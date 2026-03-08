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
