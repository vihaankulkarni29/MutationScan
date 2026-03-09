######## snakemake preamble start (automatically inserted, do not edit) ########
import sys;sys.path.extend(['C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages', 'C:\\Users\\Vihaan\\Desktop\\MutationScan', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Scripts', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\DLLs', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow', 'C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages', 'C:\\Users\\Vihaan\\Documents\\AutoDock\\src', 'C:\\Users\\Vihaan\\Documents\\FastaAAExtractor\\src', 'C:\\Users\\Vihaan\\Desktop\\MutationScan\\src', 'C:\\Users\\Vihaan\\MutationScan\\subscan\\src', 'C:\\Users\\Vihaan\\Documents\\wildtype-aligner\\src', 'C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages\\win32', 'C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages\\win32\\lib', 'C:\\Users\\Vihaan\\AppData\\Roaming\\Python\\Python313\\site-packages\\Pythonwin', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib\\site-packages', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib\\site-packages\\win32', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib\\site-packages\\win32\\lib', 'C:\\Users\\Vihaan\\miniconda3\\envs\\amrflow\\Lib\\site-packages\\Pythonwin', 'C:\\Users\\Vihaan\\Desktop\\MutationScan\\src\\scripts', 'C:\\Users\\Vihaan\\Desktop\\MutationScan\\src/scripts']);import pickle;from snakemake import script;script.snakemake = pickle.loads(b'\x80\x04\x95\x12\x05\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94\x8c-temp_data_collection/pure_indian_metadata.csv\x94a}\x94(\x8c\x06_names\x94}\x94\x8c\x07raw_csv\x94K\x00N\x86\x94s\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x12h\x06\x8c\x0eAttributeGuard\x94\x93\x94)\x81\x94}\x94\x8c\x04name\x94h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbh\x0eh\nub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c!data/results/curated_metadata.csv\x94\x8c\x14data/results/genomes\x94e}\x94(h\x0c}\x94(\x8c\x0bcurated_csv\x94K\x00N\x86\x94\x8c\x0bgenomes_dir\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbh#h\x1fh%h ub\x8c\r_params_store\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94(\x8c\x13default@example.com\x94\x8c\x00\x94e}\x94(h\x0c}\x94(\x8c\x05email\x94K\x00N\x86\x94\x8c\x07api_key\x94K\x01N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbh4h0h6h1ub\x8c\r_params_types\x94}\x94\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c"C:\\Users\\Vihaan\\AppData\\Local\\Temp\x94e}\x94(h\x0c}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbhRK\x01hTK\x01hVhOub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94}\x94(h\x0c}\x94h\x10]\x94(h\x12h\x13eh\x12h\x15)\x81\x94}\x94h\x18h\x12sbh\x13h\x15)\x81\x94}\x94h\x18h\x13sbub\x8c\x06config\x94}\x94(\x8c\x05email\x94h0\x8c\x07api_key\x94h1\x8c\tinput_csv\x94\x8c-temp_data_collection/pure_indian_metadata.csv\x94\x8c\rlocal_genomes\x94h1\x8c\x0efilter_country\x94\x8c\x05India\x94\x8c\x0ffilter_min_year\x94M\xdf\x07\x8c\x0ctargets_file\x94\x8c\x16config/acr_targets.txt\x94\x8c\runiprot_taxid\x94\x8c\x0583333\x94\x8c\nlocal_refs\x94h1u\x8c\x04rule\x94\x8c\x0facquire_genomes\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8c0C:\\Users\\Vihaan\\Desktop\\MutationScan\\src/scripts\x94ub.');del script;from snakemake.logging import logger;from snakemake.script import snakemake;__real_file__ = __file__; __file__ = 'C:\\Users\\Vihaan\\Desktop\\MutationScan\\src\\scripts\\01_acquire_genomes.py';
######## snakemake preamble end #########
import logging
import sys
from pathlib import Path

# Adjust path so Snakemake can find the mutation_scan package
sys.path.append(str(Path(__file__).resolve().parents[2]))

from mutation_scan.core.ingestion_engine import GenomicIngestionEngine
from mutation_scan.core.metadata_interrogator import MetadataInterrogator
from mutation_scan.core.universal_downloader import UniversalGenomeDownloader

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Snakemake context injection
input_csv = snakemake.input.raw_csv
curated_csv = snakemake.output.curated_csv
genomes_dir = Path(snakemake.output.genomes_dir)

email = snakemake.params.email
api_key = snakemake.params.api_key

results_dir = Path(curated_csv).parent
results_dir.mkdir(parents=True, exist_ok=True)
genomes_dir.mkdir(parents=True, exist_ok=True)

logging.info(f"Step 1.1: Ingesting and resolving IDs from {input_csv}...")
ingestion_engine = GenomicIngestionEngine(email=email, api_key=api_key)
resolved_df = ingestion_engine.route_and_resolve_input(input_csv)

logging.info("Step 1.2: Interrogating and filtering metadata...")
interrogator = MetadataInterrogator(email=email, api_key=api_key, output_dir=results_dir)
curated_df, rejected_df = interrogator.interrogate_and_filter(resolved_df)

if not curated_df.empty:
    logging.info("Step 1.3: Downloading genomes via ThreadPool (Idempotent)...")
    downloader = UniversalGenomeDownloader(api_key=api_key, genomes_dir=genomes_dir)
    downloader.download_curated_genomes(str(curated_csv))
else:
    logging.error("No genomes passed the filters. Halting.")
    sys.exit(1)
