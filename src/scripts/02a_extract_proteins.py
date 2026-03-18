"""
Phase 1a: Protein Extraction

This script:
1. Extracts protein sequences from genomic DNA using tblastn
2. Writes extracted protein FASTA files for downstream variant calling

Snakemake Context Injection:
- Input:  genomes_dir, targets_file
- Output: proteins_dir, refs_dir
- Params: uniprot_taxid (optional)
"""

import logging
import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

from mutation_scan.core.tblastn_extractor import TblastnSequenceExtractor

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
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
uniprot_taxid = snakemake.params.uniprot_taxid

# ---------------------------------------------------------
# SANITY CHECKS & SETUP
# ---------------------------------------------------------
if not genomes_dir.exists():
    logger.error(f"CRITICAL: Genomes directory not found: {genomes_dir}")
    sys.exit(1)

if not targets_file.exists():
    logger.error(f"CRITICAL: Targets file not found: {targets_file}")
    sys.exit(1)

proteins_dir.mkdir(parents=True, exist_ok=True)
refs_dir.mkdir(parents=True, exist_ok=True)

logger.info("Phase 1a Configuration:")
logger.info(f"  Genomes Dir: {genomes_dir}")
logger.info(f"  Targets File: {targets_file}")
logger.info(f"  Proteins Output: {proteins_dir}")
logger.info(f"  References Dir: {refs_dir}")
logger.info(f"  UniProt TaxID: {uniprot_taxid if uniprot_taxid else 'None (local refs only)'}")

# ---------------------------------------------------------
# STEP 1.4: LOAD TARGET GENES
# ---------------------------------------------------------
with open(targets_file, 'r') as f:
    target_genes = [line.strip() for line in f if line.strip()]

if not target_genes:
    logger.error(f"CRITICAL: No target genes loaded from {targets_file}")
    sys.exit(1)

logger.info(f"Target genes loaded: {target_genes}")

# ---------------------------------------------------------
# STEP 1.4: EXTRACT PROTEINS VIA TBLASTN
# ---------------------------------------------------------
logger.info("Step 1.4: Initializing Autonomous TblastnSequenceExtractor...")

extractor = TblastnSequenceExtractor(
    genomes_dir=genomes_dir,
    refs_dir=refs_dir,
    output_dir=proteins_dir,
    tblastn_binary="tblastn",
    uniprot_taxid=uniprot_taxid if uniprot_taxid else None
)

# ---------------------------------------------------------
# THE BLAST SHIELD: Pre-flight Reference Validation
# ---------------------------------------------------------
missing_refs = []
for gene in target_genes:
    ref_fasta = refs_dir / f"{gene}.fasta"
    ref_faa = refs_dir / f"{gene}_WT.faa"

    # Check if either valid file exists and has content
    fasta_valid = ref_fasta.exists() and ref_fasta.stat().st_size > 0
    faa_valid = ref_faa.exists() and ref_faa.stat().st_size > 0

    if not (fasta_valid or faa_valid):
        missing_refs.append(gene)

if missing_refs:
    logging.error(f"CRITICAL: Missing or empty reference proteins for: {missing_refs}")
    logging.error("Auto-fetch failed or local files are missing. Cannot proceed with extraction.")
    sys.exit(1)
# ---------------------------------------------------------

# Validate that genomes exist
genome_files = list(genomes_dir.glob("*.fna"))
if not genome_files:
    logger.error(f"CRITICAL: No .fna files found in {genomes_dir}")
    sys.exit(1)

logger.info(f"Found {len(genome_files)} genome files to extract from")

# Extract proteins from all genomes
genome_ids = [f.stem for f in genome_files]
logger.info(f"Extracting target genes from {len(genome_ids)} genomes...")
extractor.extract_all_genomes(genome_ids, target_genes)

# Validate extraction output
protein_files = list(proteins_dir.glob("*.faa"))
if not protein_files:
    logger.error(f"CRITICAL: No .faa files were extracted to {proteins_dir}")
    sys.exit(1)

logger.info(f"Extraction complete: {len(protein_files)} protein files generated")

# ---------------------------------------------------------
# COMPLETION
# ---------------------------------------------------------
logger.info("Phase 1a Complete!")
logger.info(f"  Proteins extracted: {len(protein_files)}")
