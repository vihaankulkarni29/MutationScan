"""
Phase 1a: Protein Extraction

Extracts protein sequences from genomic DNA using tblastn against user-supplied
or auto-fetched wild-type references.

Snakemake Context Injection:
- Input:  genomes_dir, targets_file
- Output: proteins_dir, refs_dir
- Params: uniprot_taxid, reference_seed_dir (optional)
"""

import logging
import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

from mutation_scan.core.tblastn_extractor import (
    TblastnSequenceExtractor,
    load_target_genes,
    missing_reference_genes,
    resolve_reference,
    seed_references_from_dir,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# ---------------------------------------------------------
# SNAKEMAKE CONTEXT INJECTION
# ---------------------------------------------------------
genomes_dir = Path(snakemake.input.genomes_dir)
targets_file = Path(snakemake.input.targets_file)
RESULTS_DIR = Path(snakemake.params.out_dir)
os.makedirs(RESULTS_DIR, exist_ok=True)

proteins_dir = Path(snakemake.output.proteins_dir)
refs_dir = Path(snakemake.output.refs_dir)
uniprot_taxid = str(snakemake.params.uniprot_taxid or "").strip()
reference_seed_dir = str(snakemake.params.reference_seed_dir or "").strip()
skip_extraction = snakemake.params.skip_extraction

# ---------------------------------------------------------
# SANITY CHECKS & SETUP
# ---------------------------------------------------------
if not genomes_dir.exists():
    logger.error("CRITICAL: Genomes directory not found: %s", genomes_dir)
    sys.exit(1)

if not targets_file.exists():
    logger.error("CRITICAL: Targets file not found: %s", targets_file)
    sys.exit(1)

proteins_dir.mkdir(parents=True, exist_ok=True)
refs_dir.mkdir(parents=True, exist_ok=True)

target_genes = load_target_genes(targets_file)
if not target_genes:
    logger.error(
        "CRITICAL: No target genes loaded from %s. "
        "Add gene names (one per line) or see config/targets.txt.example.",
        targets_file,
    )
    sys.exit(1)

logger.info("Phase 1a Configuration:")
logger.info("  Genomes Dir: %s", genomes_dir)
logger.info("  Targets File: %s", targets_file)
logger.info("  Target Genes: %s", target_genes)
logger.info("  Proteins Output: %s", proteins_dir)
logger.info("  References Dir: %s", refs_dir)
logger.info("  Reference Seed Dir: %s", reference_seed_dir or "None")
logger.info("  UniProt TaxID: %s", uniprot_taxid or "None")
logger.info("  Skip Extraction: %s", skip_extraction)

extractor = TblastnSequenceExtractor(
    genomes_dir=genomes_dir,
    refs_dir=refs_dir,
    output_dir=proteins_dir,
    tblastn_binary="tblastn",
    uniprot_taxid=uniprot_taxid or None,
)

seeded = seed_references_from_dir(
    refs_dir=refs_dir,
    seed_dir=Path(reference_seed_dir) if reference_seed_dir else Path("."),
    target_genes=target_genes,
)
if seeded:
    logger.info("Seeded %s reference(s) from %s", seeded, reference_seed_dir)

extractor._ensure_references_exist(target_genes)

missing_refs = missing_reference_genes(refs_dir, target_genes)
if missing_refs:
    logger.error("CRITICAL: Missing reference proteins for: %s", missing_refs)
    logger.error(
        "Provide references via reference_seed_dir (files named {gene}_WT.faa or "
        "{gene}.faa), set uniprot_taxid for auto-fetch, or copy files into %s.",
        refs_dir,
    )
    sys.exit(1)

# ---------------------------------------------------------
# EARLY EXIT: skip_extraction when outputs already exist
# ---------------------------------------------------------
if skip_extraction:
    protein_files = list(proteins_dir.glob("*.faa"))
    refs_ready = all(resolve_reference(refs_dir, gene) is not None for gene in target_genes)
    if protein_files and refs_ready:
        logger.info(
            "SKIP_EXTRACTION=True and outputs exist (%s protein files). Skipping.",
            len(protein_files),
        )
        marker_file = proteins_dir / ".proteins_extracted"
        marker_file.touch()
        logger.info("Phase 1a Complete (skipped)!")
        sys.exit(0)

    logger.warning("SKIP_EXTRACTION=True but outputs are incomplete. Proceeding with extraction...")

genome_files = list(genomes_dir.glob("*.fna"))
if not genome_files:
    logger.error("CRITICAL: No .fna files found in %s", genomes_dir)
    sys.exit(1)

logger.info("Found %s genome files to extract from", len(genome_files))

genome_ids = [path.stem for path in genome_files]
extractor.extract_all_genomes(genome_ids, target_genes)

protein_files = list(proteins_dir.glob("*.faa"))
if not protein_files:
    logger.error("CRITICAL: No .faa files were extracted to %s", proteins_dir)
    sys.exit(1)

logger.info("Extraction complete: %s protein files generated", len(protein_files))

marker_file = proteins_dir / ".proteins_extracted"
marker_file.touch()
logger.info("Marker file created: %s", marker_file)
logger.info("Phase 1a Complete!")
logger.info("  Proteins extracted: %s", len(protein_files))
