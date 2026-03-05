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
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
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
    ):
        """
        Initialize TblastnSequenceExtractor.

        Args:
            genomes_dir: Directory containing clinical genome .fna files
            refs_dir: Directory containing wild-type reference .faa files
            output_dir: Directory for extracted protein sequences
            tblastn_binary: Path to tblastn executable (default: assume in PATH)
        """
        self.genomes_dir = Path(genomes_dir)
        self.refs_dir = Path(refs_dir)
        self.output_dir = Path(output_dir)
        self.tblastn_binary = tblastn_binary

        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Verify dependencies
        self._verify_tblastn_available()
        logger.info(f"Initialized TblastnSequenceExtractor")
        logger.info(f"  Genomes: {self.genomes_dir}")
        logger.info(f"  References: {self.refs_dir}")
        logger.info(f"  Output: {self.output_dir}")

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
                        seq=translated_seq,
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
                timeout=60,
            )

            if result.returncode != 0:
                logger.warning(f"tblastn failed for {gene_name}: {result.stderr}")
                return None

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
        logger.info(f"Processing {len(genome_ids)} genomes")

        results = []

        for genome_id in genome_ids:
            success, fail = self.extract_with_tblastn(
                genome_id, target_genes=target_genes
            )

            results.append({
                "Genome": genome_id,
                "Extracted": success,
                "Failed": fail,
            })

        return pd.DataFrame(results)
