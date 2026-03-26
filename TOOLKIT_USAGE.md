# MutationScan Support Toolkit Usage

The consolidated utility CLI is:

- `utility scripts/mutationscan_support_toolkit.py`

Run help:

```bash
python "utility scripts/mutationscan_support_toolkit.py" --help
```

## 1) Download genomes from BV-BRC REST

```bash
python "utility scripts/mutationscan_support_toolkit.py" download-rest \
  --csv-file "data/output/Ciproflaxcin_Run/BVBRC_genome_amr.csv" \
  --output-dir "data/output/Ciproflaxcin_Run/genomes_full" \
  --genome-id-column "Genome ID" \
  --api-limit 5000 \
  --retries 4 \
  --timeout 45
```

## 2) Fetch BV-BRC metadata and enrich CSV

```bash
python "utility scripts/mutationscan_support_toolkit.py" fetch-metadata \
  --input-csv "data/output/Ciprofloxacin_Run/BVBRC_genome_amr_Cipro.csv" \
  --genome-id-column "Genome ID" \
  --output-metadata-csv "data/output/Ciprofloxacin_Run/BVBRC_genome_amr_Cipro_genome_metadata.csv" \
  --output-enriched-csv "data/output/Ciprofloxacin_Run/BVBRC_genome_amr_Cipro_with_metadata.csv" \
  --failed-log "data/output/Ciprofloxacin_Run/BVBRC_genome_amr_Cipro_metadata_failed_ids.txt"
```

## 3) Build geospatial mutation matrix

```bash
python "utility scripts/mutationscan_support_toolkit.py" geospatial-matrix \
  --metadata-csvs \
    "data/output/Tigecycline_Run/BVBRC_genome_amr_Tigecycline_with_metadata.csv" \
    "data/output/Tetracycline_Run/BVBRC_genome_amr_tetracycline_with_metadata.csv" \
    "data/output/Ciprofloxacin_Run/BVBRC_genome_amr_Cipro_with_metadata.csv" \
    "data/output/Chloramphenicol_Run/BVBRC_genome_amr_Chloramphenicol_with_metadata.csv" \
  --genomics-reports \
    "data/output/Tigecycline_Run/1_genomics_report.csv" \
    "data/output/Tetracycline_Run/1_genomics_report.csv" \
    "data/output/Ciprofloxacin_Run/1_genomics_report.csv" \
    "data/output/Chloramphenicol_Run/1_genomics_report.csv" \
  --regulatory-genes marR acrR \
  --output-dir "data/output/geospatial/Master"
```

## 4) Generate geospatial heatmap

```bash
python "utility scripts/mutationscan_support_toolkit.py" geospatial-heatmap \
  --input-csv "data/output/geospatial/Master/Merged_Regulatory_Geospatial_AMR.csv" \
  --output-matrix-csv "data/output/geospatial/Master/Geospatial_Mutation_Matrix_Improved.csv" \
  --output-plot "data/presentation/geospatial_heatmap_improved.jpg" \
  --top-n 15
```

## 5) Assembly quality control and extraction-ready subset

```bash
python "utility scripts/mutationscan_support_toolkit.py" qc-genomes \
  --input-dir "data/output/Ciproflaxcin_Run/genomes_full" \
  --output-summary-csv "data/output/Ciproflaxcin_Run/genome_quality_summary_full.csv" \
  --output-ready-dir "data/output/Ciproflaxcin_Run/genomes_extraction_ready"
```

## 6) Presentation plots

```bash
python "utility scripts/mutationscan_support_toolkit.py" presentation-plots \
  --runs Tetracycline_Run Chloramphenicol_Run Tigecycline_Run \
  --base-output-dir "data/output" \
  --output-dir "data/presentation"
```
