######## snakemake preamble start (automatically inserted, do not edit) ########
import sys;sys.path.extend(['C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages', 'C:\\Users\\Vihaan\\Desktop\\MutationScan', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Scripts', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\DLLs', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow', 'C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages', 'C:\\Users\\Vihaan\\Documents\\AutoDock\\src', 'C:\\Users\\Vihaan\\Documents\\FastaAAExtractor\\src', 'C:\\Users\\Vihaan\\Desktop\\MutationScan\\src', 'C:\\Users\\Vihaan\\MutationScan\\subscan\\src', 'C:\\Users\\Vihaan\\Documents\\wildtype-aligner\\src', 'C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages\\win32', 'C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages\\win32\\lib', 'C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages\\Pythonwin', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib\\site-packages', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib\\site-packages\\win32', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib\\site-packages\\win32\\lib', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib\\site-packages\\Pythonwin', 'C:\\Users\\Vihaan\\Desktop\\MutationScan\\src\\scripts', 'C:\\Users\\Vihaan\\Desktop\\MutationScan\\src/scripts']);import pickle;from snakemake import script;script.snakemake = pickle.loads(b'\x80\x04\x95Q\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c\x14data/results/genomes\x94\x8c\x16config/acr_targets.txt\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\x0bgenomes_dir\x94K\x00N\x86\x94\x8c\x0ctargets_file\x94K\x01N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x15h\x06\x8c\x0eAttributeGuard\x94\x93\x94)\x81\x94}\x94\x8c\x04name\x94h\x15sbh\x16h\x18)\x81\x94}\x94h\x1bh\x16sbh\x0fh\nh\x11h\x0bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c\x15data/results/proteins\x94\x8c"data/results/1_genomics_report.csv\x94e}\x94(h\r}\x94(\x8c\x0cproteins_dir\x94K\x00N\x86\x94\x8c\rmutations_csv\x94K\x01N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x18)\x81\x94}\x94h\x1bh\x15sbh\x16h\x18)\x81\x94}\x94h\x1bh\x16sbh&h"h(h#ub\x8c\r_params_store\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(\x8c\x11data/results/refs\x94\x8c\x0583333\x94e}\x94(h\r}\x94(\x8c\x08refs_dir\x94K\x00N\x86\x94\x8c\runiprot_taxid\x94K\x01N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x18)\x81\x94}\x94h\x1bh\x15sbh\x16h\x18)\x81\x94}\x94h\x1bh\x16sbh7h3h9h4ub\x8c\r_params_types\x94}\x94\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x13]\x94(h\x15h\x16eh\x15h\x18)\x81\x94}\x94h\x1bh\x15sbh\x16h\x18)\x81\x94}\x94h\x1bh\x16sbub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c"C:\\Users\\Vihaan\\AppData\\Local\\Temp\x94e}\x94(h\r}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x18)\x81\x94}\x94h\x1bh\x15sbh\x16h\x18)\x81\x94}\x94h\x1bh\x16sbhUK\x01hWK\x01hYhRub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x13]\x94(h\x15h\x16eh\x15h\x18)\x81\x94}\x94h\x1bh\x15sbh\x16h\x18)\x81\x94}\x94h\x1bh\x16sbub\x8c\x06config\x94}\x94(\x8c\x05email\x94\x8c\x13default@example.com\x94\x8c\x07api_key\x94\x8c\x00\x94\x8c\tinput_csv\x94\x8c-temp_data_collection/pure_indian_metadata.csv\x94\x8c\rlocal_genomes\x94hp\x8c\x0efilter_country\x94\x8c\x05India\x94\x8c\x0ffilter_min_year\x94M\xdf\x07\x8c\x0ctargets_file\x94\x8c\x16config/acr_targets.txt\x94\x8c\runiprot_taxid\x94h4\x8c\nlocal_refs\x94hpu\x8c\x04rule\x94\x8c\x10extract_and_call\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c0C:\\Users\\Vihaan\\Desktop\\MutationScan\\src/scripts\x94ub.');del script;from snakemake.logging import logger;from snakemake.script import snakemake;__real_file__ = __file__; __file__ = 'C:\\Users\\Vihaan\\Desktop\\MutationScan\\src\\scripts\\02_extract_and_call.py';
######## snakemake preamble end #########
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
