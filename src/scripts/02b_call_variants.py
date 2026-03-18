"""
Phase 1b: Variant Calling

This script:
1. Reads pre-extracted protein FASTA files
2. Calls variants by comparing against wild-type references
3. Writes 1_genomics_report.csv for downstream analysis

Snakemake Context Injection:
- Input:  proteins_dir, refs_dir
- Output: report (or mutations_csv)
"""

import logging
import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

from mutation_scan.analysis.variant_caller import VariantCaller

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# ---------------------------------------------------------
# SNAKEMAKE CONTEXT INJECTION
# ---------------------------------------------------------
proteins_dir = Path(snakemake.input.proteins_dir)
refs_dir = Path(snakemake.input.refs_dir)
RESULTS_DIR = Path(snakemake.params.out_dir)
os.makedirs(RESULTS_DIR, exist_ok=True)

mutations_csv = Path(getattr(snakemake.output, "mutations_csv", snakemake.output.report))

# ---------------------------------------------------------
# SANITY CHECKS
# ---------------------------------------------------------
if not proteins_dir.exists():
    logger.error(f"CRITICAL: Proteins directory not found: {proteins_dir}")
    sys.exit(1)

if not refs_dir.exists():
    logger.error(f"CRITICAL: References directory not found: {refs_dir}")
    sys.exit(1)

logger.info("Phase 1b Configuration:")
logger.info(f"  Proteins Dir: {proteins_dir}")
logger.info(f"  References Dir: {refs_dir}")
logger.info(f"  Genomics Report: {mutations_csv}")

# ---------------------------------------------------------
# STEP 1.5: CALL VARIANTS FROM PROTEINS
# ---------------------------------------------------------
logger.info("Step 1.5: Initializing VariantCaller for mutation detection...")

caller = VariantCaller(
    refs_dir=refs_dir,
    enable_ml=False
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
logger.info("Phase 1b Complete!")
logger.info(f"  Mutations called: {len(mutations_df)}")
logger.info(f"  Genomics report saved: {mutations_csv}")
