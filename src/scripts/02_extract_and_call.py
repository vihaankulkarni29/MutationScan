"""
Phase 1b: Protein Extraction & Variant Calling

This script:
1. Extracts protein sequences from genomic DNA using tblastn
2. Calls variants by comparing to wild-type references
3. Generates genomics report for downstream analysis

Snakemake Context Injection:
- Input:  genomes_dir (may be local or downloaded)
                targets_file (list of target genes)
- Output: proteins_dir, refs_dir, mutations_csv
- Params: uniprot_taxid (optional, for auto-fetch of references)
"""

import logging
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

from mutation_scan.core.tblastn_extractor import TblastnSequenceExtractor
from mutation_scan.analysis.variant_caller import VariantCaller

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ---------------------------------------------------------
# SNAKEMAKE CONTEXT INJECTION
# ---------------------------------------------------------
genomes_dir = Path(snakemake.input.genomes_dir)
targets_file = Path(snakemake.input.targets_file)

proteins_dir = Path(snakemake.output.proteins_dir)
mutations_csv = Path(snakemake.output.mutations_csv)

refs_dir = Path(snakemake.params.refs_dir)
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

logger.info(f"Phase 1b Configuration:")
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
# STEP 1.5: CALL VARIANTS FROM PROTEINS
# ---------------------------------------------------------
logger.info("Step 1.5: Initializing VariantCaller for mutation detection...")

caller = VariantCaller(
    refs_dir=refs_dir,
    enable_ml=False  # Disable ML for Phase 1 (can be enabled later in Phase 2-5)
)

logger.info("Calling variants from extracted proteins...")
mutations_df = caller.call_variants(
    proteins_dir=proteins_dir,
    output_csv=mutations_csv
)

if mutations_df.empty:
    logger.warning("No mutations found across the cohort.")
else:
    logger.info(f"Variant calling complete: {len(mutations_df)} mutations identified")

# ---------------------------------------------------------
# COMPLETION
# ---------------------------------------------------------
logger.info(f"Phase 1b Complete!")
logger.info(f"  Proteins extracted: {len(protein_files)}")
logger.info(f"  Mutations called: {len(mutations_df)}")
logger.info(f"  Genomics report saved: {mutations_csv}")
