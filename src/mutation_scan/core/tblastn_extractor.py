"""
tblastn-based Sequence Extractor for Dynamic Translation.

Replaces ABRicate + manual translation with translating aligner tblastn
to prevent frameshift artifacts ("The Alanine Trap").

Workflow:
1. Use tblastn to align protein reference against genomic DNA
2. Extract perfectly translated, in-frame amino acid sequence (column 13)
3. Remove alignment gaps and output clean protein FASTA
"""

import logging
import re
import subprocess
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class TblastnSequenceExtractor:
    """
    Dynamic sequence extraction using tblastn for translating alignment.

    Features:
    - Queries protein references against genomic DNA (wildtype or clinical)
    - Preserves reading frame through native translation
    - Extracts perfectly aligned amino acid sequences
    - Removes alignment gaps for clean protein output
    - Graceful handling of partial/missing genes
    """

    def __init__(
        self,
        genomes_dir: Path,
        refs_dir: Path,
        output_dir: Path,
        tblastn_binary: str = "tblastn",
        uniprot_taxid: Optional[str] = None,
    ):
        """
        Initialize TblastnSequenceExtractor.

        Args:
            genomes_dir: Directory containing clinical genome .fna files
            refs_dir: Directory containing wild-type reference .faa files
            output_dir: Directory for extracted protein sequences
            tblastn_binary: Path to tblastn executable (default: assume in PATH)
            uniprot_taxid: Optional UniProt Taxonomy ID for reference auto-fetching (e.g., 83333 for E. coli K-12)
        """
        self.genomes_dir = Path(genomes_dir)
        self.refs_dir = Path(refs_dir)
        self.output_dir = Path(output_dir)
        self.tblastn_binary = tblastn_binary
        self.uniprot_taxid = uniprot_taxid

        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Verify dependencies
        self._verify_tblastn_available()
        logger.info(f"Initialized TblastnSequenceExtractor")
        logger.info(f"  Genomes: {self.genomes_dir}")
        logger.info(f"  References: {self.refs_dir}")
        logger.info(f"  Output: {self.output_dir}")
        if self.uniprot_taxid:
            logger.info(f"  UniProt TaxID (for auto-fetch): {self.uniprot_taxid}")
        else:
            logger.info(f"  UniProt TaxID: None (using local refs only)")

    def _verify_tblastn_available(self) -> None:
        """Check if tblastn is available in PATH or as provided."""
        try:
            result = subprocess.run(
                [self.tblastn_binary, "-version"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if result.returncode == 0:
                logger.info(f"tblastn verified: {result.stdout.splitlines()[0]}")
                return
        except FileNotFoundError:
            logger.error(f"tblastn not found: {self.tblastn_binary}")
            logger.info("Install BLAST+ package and ensure tblastn is in PATH")
            raise EnvironmentError(
                f"tblastn binary not found at: {self.tblastn_binary}"
            )
        except Exception as e:
            logger.error(f"Failed to verify tblastn: {e}")
            raise

    def extract_with_tblastn(
        self,
        genome_id: str,
        target_genes: Optional[List[str]] = None,
    ) -> Tuple[int, int]:
        """
        Extract sequences for target genes from a clinical genome using tblastn.

        Workflow:
        1. For each reference protein (.faa), run tblastn against clinical genome .fna
        2. Parse output format 6 (tabular) to extract translated sequence (column 13)
        3. Remove gaps (-) from aligned sequence
        4. Save as individual protein FASTA for variant calling

        Args:
            genome_id: Genome identifier (e.g., 'GCF_000005845')
            target_genes: List of specific genes to extract (e.g., ['acrA', 'gyrA'])
                         If None, extract all available reference genes

        Returns:
            Tuple of (successfully_extracted, failed)
        """
        logger.info(f"Extracting sequences from {genome_id} using tblastn")

        # Locate genome FASTA
        genome_fna = self.genomes_dir / f"{genome_id}.fna"
        if not genome_fna.exists():
            logger.error(f"Genome file not found: {genome_fna}")
            return 0, 0

        # Get list of reference proteins
        ref_files = list(self.refs_dir.glob("*.faa"))
        if not ref_files:
            logger.warning(f"No reference protein files found in {self.refs_dir}")
            return 0, 0

        # Filter by target genes if specified
        if target_genes:
            target_genes_lower = [g.lower() for g in target_genes]
            ref_files = [
                f for f in ref_files
                if any(tg in f.stem.lower() for tg in target_genes_lower)
            ]
            if not ref_files:
                logger.warning(
                    f"No reference files match target genes: {target_genes}"
                )
                return 0, 0

        success_count = 0
        fail_count = 0

        # Snakemake may clean output directories after upstream failures.
        # Recreate defensively before writing files for this genome.
        self.output_dir.mkdir(parents=True, exist_ok=True)

        for ref_faa in ref_files:
            gene_name = ref_faa.stem.replace("_WT", "")

            try:
                # Run tblastn alignment
                translated_seq = self._run_tblastn_alignment(
                    ref_faa, genome_fna, gene_name
                )

                if translated_seq:
                    # Save as protein FASTA
                    output_faa = self.output_dir / f"{genome_id}_{gene_name}.faa"
                    record = SeqRecord(
                        seq=Seq(translated_seq),
                        id=f"{genome_id}_{gene_name}",
                        description=f"tblastn-extracted {gene_name} from {genome_id}",
                    )
                    SeqIO.write(record, str(output_faa), "fasta")
                    logger.info(f"  ✓ {gene_name}: {len(translated_seq)} aa")
                    success_count += 1
                else:
                    logger.warning(f"  ✗ {gene_name}: No alignment found")
                    fail_count += 1

            except Exception as e:
                logger.error(f"  ✗ {gene_name}: {e}")
                fail_count += 1

        logger.info(
            f"Extraction complete: {success_count} success, {fail_count} failures"
        )
        return success_count, fail_count

    def _run_tblastn_alignment(
        self,
        ref_faa: Path,
        genome_fna: Path,
        gene_name: str,
    ) -> Optional[str]:
        """
        Execute tblastn and extract translated sequence from column 13.

        Command:
            tblastn -query <ref_faa> -subject <genome_fna> \\
              -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq" \\
              -max_target_seqs 1 -max_hsps 1 > temp_output.txt

        Output Format 6 columns:
            1. qseqid: Query sequence ID
            2. sseqid: Subject sequence ID (contig)
            3. pident: Percent identity
            4. length: Alignment length
            5. mismatch: Number of mismatches
            6. gapopen: Number of gap openings
            7. qstart: Query start
            8. qend: Query end
            9. sstart: Subject start
            10. send: Subject end
            11. evalue: Expect value
            12. bitscore: Bit score
            13. sseq: **ALIGNED SUBJECT SEQUENCE (translated ORF)**

        Args:
            ref_faa: Path to reference protein FASTA
            genome_fna: Path to subject genome FASTA (nucleotides)
            gene_name: Gene name for logging

        Returns:
            Cleaned protein sequence (gaps removed) or None if no alignment
        """
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as tmp:
            tmp_output = tmp.name

        try:
            command = [
                self.tblastn_binary,
                "-query", str(ref_faa),
                "-subject", str(genome_fna),
                "-outfmt",
                "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq",
                "-max_target_seqs", "1",
                "-max_hsps", "1",
                "-out", tmp_output,
            ]

            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=True,
                timeout=60,
            )

            # Parse output
            with open(tmp_output, "r") as f:
                content = f.read().strip()

            if not content:
                # No alignment found
                return None

            # Extract column 13 (subject translated sequence)
            columns = content.split("\t")
            if len(columns) < 13:
                logger.warning(f"Unexpected output format for {gene_name}")
                return None

            aligned_seq = columns[12]  # Index 12 = column 13 (0-based)

            # Remove alignment gaps
            cleaned_seq = aligned_seq.replace("-", "")

            return cleaned_seq if cleaned_seq else None

        except subprocess.TimeoutExpired:
            logger.error(f"tblastn timeout for {gene_name}")
            return None
        except subprocess.CalledProcessError as e:
            logger.warning(
                f"tblastn failed for {gene_name} (exit {e.returncode}): "
                f"{(e.stderr or '').strip()}"
            )
            return None
        except Exception as e:
            logger.error(f"Error running tblastn for {gene_name}: {e}")
            return None
        finally:
            # Cleanup temp file
            Path(tmp_output).unlink(missing_ok=True)

    def extract_all_genomes(
        self,
        genome_ids: List[str],
        target_genes: Optional[List[str]] = None,
    ) -> pd.DataFrame:
        """
        Extract sequences from multiple clinical genomes.

        Returns DataFrame with extraction statistics.

        Args:
            genome_ids: List of genome IDs to process
            target_genes: Optional list of genes to extract

        Returns:
            DataFrame with columns: Genome, Gene, Success (bool)
        """
        # Ensure all required reference sequences exist (fetch from UniProt if needed)
        if target_genes:
            self._ensure_references_exist(target_genes)
        
        total = len(genome_ids)
        gene_count = len(target_genes) if target_genes else "all"
        logger.info(
            f"Starting extraction for {total} genomes across {gene_count} targets..."
        )

        results = []
        fully_successful = 0
        partial_or_failed = 0

        for idx, genome_id in enumerate(genome_ids, start=1):
            logger.info(f"[{idx}/{total}] Extracting {genome_id}...")
            try:
                success, fail = self.extract_with_tblastn(
                    genome_id, target_genes=target_genes
                )
            except Exception as e:
                # Graceful degradation: continue to next genome.
                logger.warning(
                    f"  -> Failed to process genome {genome_id}: {e}"
                )
                success, fail = 0, 1

            results.append({
                "Genome": genome_id,
                "Extracted": success,
                "Failed": fail,
            })

            if fail == 0 and success > 0:
                fully_successful += 1
            else:
                partial_or_failed += 1

        logger.info("=" * 50)
        logger.info(
            f"Extraction Complete: {fully_successful} fully successful, "
            f"{partial_or_failed} partial/failed."
        )
        logger.info("=" * 50)

        return pd.DataFrame(results)

    def _ensure_references_exist(self, target_genes: List[str]) -> None:
        """
        Ensure all required reference sequences exist, auto-fetching from UniProt if needed.
        
        Dynamic Reference Fetching Strategy:
        - If --uniprot-taxid is provided: Attempt to fetch missing references from UniProt
        - If --uniprot-taxid is omitted: Rely strictly on local files in refs/ directory
        
        Args:
            target_genes: List of gene names to check/fetch
        """
        self.refs_dir.mkdir(parents=True, exist_ok=True)
        
        if not self.uniprot_taxid:
            logger.info("No --uniprot-taxid provided. Relying strictly on local reference FASTA files in refs/ directory.")
            return
        
        logger.info(f"Checking references for {len(target_genes)} genes...")
        for gene in target_genes:
            ref_path = self.refs_dir / f"{gene.lower()}.fasta"
            if ref_path.exists():
                logger.debug(f"Reference exists: {ref_path.name}")
                continue
            
            logger.info(f"Reference for {gene} missing. Auto-fetching from UniProt for TaxID {self.uniprot_taxid}...")
            try:
                # Query UniProt for the reviewed (canonical) reference protein
                query = f"gene:{gene}+AND+taxonomy_id:{self.uniprot_taxid}+AND+reviewed:true"
                url = f"https://rest.uniprot.org/uniprotkb/search?query={urllib.parse.quote(query)}&format=fasta&size=1"
                
                req = urllib.request.Request(url, headers={'Accept': 'text/plain'})
                with urllib.request.urlopen(req, timeout=10) as response:
                    fasta_data = response.read().decode('utf-8').strip()
                
                if fasta_data:
                    with open(ref_path, 'w') as f:
                        f.write(fasta_data)
                    logger.info(f"Successfully saved canonical reference for {gene} to {ref_path.name}")
                else:
                    logger.warning(f"Could not find Reviewed UniProt reference for {gene} (TaxID: {self.uniprot_taxid}).")
                
                # UniProt rate limiting (be respectful)
                time.sleep(0.5)
                
            except urllib.error.HTTPError as e:
                logger.error(f"HTTP error fetching reference for {gene}: {e.code} {e.reason}")
            except urllib.error.URLError as e:
                logger.error(f"Network error fetching reference for {gene}: {e}")
            except Exception as e:
                logger.error(f"Failed to fetch reference for {gene}: {e}")

