# GenomeExtractor Module (NCBI Datasets API v2)

## Overview

The refactored `GenomeExtractor` module provides modern, high-performance genome downloading from NCBI using the **Datasets API v2**. It replaces the legacy Entrez.efetch implementation with batch processing, comprehensive error handling, and clinical metadata extraction.

## Key Features

✅ **Batch Processing**: Download multiple genomes in a single operation (vs. sequential requests)  
✅ **Full Assemblies**: Includes chromosomes + plasmids (captures resistance genes in plasmids)  
✅ **Clinical Metadata**: Extracts host, collection date, location, isolation source, strain info  
✅ **Dual-Mode Input**: Search mode (query string) or file mode (accession list)  
✅ **QC Checks**: Validates coverage (≥90%) and minimum length thresholds  
✅ **Comprehensive Logging**: Detailed error tracking to `genome_extractor.log`  
✅ **Anti-Hallucination**: Defaults missing fields to "N/A", never invents data  

## API Reference

### Class: `NCBIDatasetsGenomeDownloader`

Main class for genome downloading and metadata extraction.

```python
from mutation_scan.genome_extractor import NCBIDatasetsGenomeDownloader
from pathlib import Path

downloader = NCBIDatasetsGenomeDownloader(
    email="your.email@example.com",  # REQUIRED
    api_key=None,                     # Optional NCBI API key
    output_dir=Path("data/genomes"),  # Output directory
    log_file=Path("data/logs/genome_extractor.log")  # Log file
)
```

#### Constructor Parameters

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `email` | str | ✓ | Your email (required by NCBI policy) |
| `api_key` | str | Optional | NCBI API key (enables faster/priority requests) |
| `output_dir` | Path | Optional | Directory to save genomes (default: `data/genomes`) |
| `log_file` | Path | Optional | Path to error log (default: `data/logs/genome_extractor.log`) |

#### Method: `search_accessions()`

Search NCBI Assembly database and resolve to Accession IDs.

```python
accessions = downloader.search_accessions(
    query="Escherichia coli AND Antibiotic Resistant",
    max_results=100
)
# Returns: ["GCF_000005845.2", "GCF_000007045.1", ...]
```

**Parameters:**
- `query` (str): Search query (uses Entrez syntax)
- `max_results` (int): Maximum accessions to return (default: 100)

**Returns:**
- List[str]: Assembly Accessions in GCF_/GCA_ format

**Raises:**
- Exception: If search fails after 3 retry attempts

**Details:**
- Uses Bio.Entrez.esearch() → Bio.Entrez.esummary()
- Automatically retries with exponential backoff on failure
- Logs all attempts to `genome_extractor.log`
- Returns empty list if no results found

---

#### Method: `read_accession_file()`

Read accession list from a text file.

```python
accessions = downloader.read_accession_file(Path("accessions.txt"))
# accessions.txt format (one per line):
# GCF_000005845.2
# GCF_000007045.1
```

**Parameters:**
- `filepath` (Path): Path to accession list file

**Returns:**
- List[str]: Accessions from file (skips empty lines and comments)

**Raises:**
- FileNotFoundError: If file doesn't exist

**Details:**
- Skips empty lines and lines starting with `#`
- Strips whitespace automatically

---

#### Method: `download_batch()`

Batch download genomes from list of accessions.

```python
successful, failed = downloader.download_batch(accessions)
# Returns: (50, 2) - 50 successful, 2 failed
```

**Parameters:**
- `accessions` (List[str]): List of Assembly Accessions

**Returns:**
- Tuple[int, int]: (successful_downloads, failed_downloads)

**Behavior:**
- Downloads full assemblies (chromosome + plasmids)
- Automatically processes zip files
- Extracts metadata from data_report.jsonl
- Saves individual `{Accession}.fasta` files
- Generates `metadata_master.csv`
- Comprehensive logging to `genome_extractor.log`

---

### Class: `GenomeProcessor`

Validate and QC check downloaded genomes.

```python
from mutation_scan.genome_extractor import GenomeProcessor

processor = GenomeProcessor(
    min_coverage=90.0,    # 90% coverage threshold
    min_length=1000000    # 1MB minimum length
)
```

#### Method: `validate_genome()`

Validate FASTA format and length.

```python
is_valid, message = processor.validate_genome(Path("GCF_000005845.2.fasta"))
# Returns: (True, "Valid FASTA: 4641652bp")
```

**Returns:**
- Tuple[bool, str]: (is_valid, status_message)

---

#### Method: `calculate_coverage()`

Calculate sequence coverage percentage.

```python
coverage = processor.calculate_coverage(Path("genome.fasta"))
# Returns: 95.5 (95.5% coverage)
```

**Parameters:**
- `filepath` (Path): Path to FASTA file
- `reference_length` (int): Optional reference length (bp)

**Returns:**
- float: Coverage percentage

---

## Output Format

### FASTA Files

Individual files named with strict convention:

```
data/genomes/
├── GCF_000005845.2.fasta
├── GCF_000007045.1.fasta
├── GCF_000009065.1.fasta
└── metadata_master.csv
```

Each file contains:
- Full chromosome sequence
- All plasmid sequences (concatenated)
- Standard NCBI FASTA headers

### Metadata CSV

`metadata_master.csv` contains clinical metadata:

| Column | Description | Example |
|--------|-------------|---------|
| Accession | Assembly Accession | GCF_000005845.2 |
| Organism Name | Species name | Escherichia coli |
| Strain | Strain identifier | K-12 substr. MG1655 |
| Collection Date | Sample collection date | 1998 |
| Host | Host organism | Human |
| Isolation Source | Sample source | Fecal samples |
| Geo Location | Geographic location | USA |
| QC Status | QC result | PASS / QC_FAIL |
| Downloaded | Download success | True / False |

**Default Value**: "N/A" for any missing field (anti-hallucination)

### Log File

`genome_extractor.log` contains detailed error tracking:

```
2024-01-15 14:32:10 - mutation_scan.genome_extractor.entrez_handler - INFO - Initialized NCBIDatasetsGenomeDownloader for user@example.com
2024-01-15 14:32:11 - mutation_scan.genome_extractor.entrez_handler - INFO - Searching NCBI Assembly for: Escherichia coli
2024-01-15 14:32:15 - mutation_scan.genome_extractor.entrez_handler - INFO - Found 45 Assembly UIDs. Resolving to accessions...
2024-01-15 14:32:20 - mutation_scan.genome_extractor.entrez_handler - INFO - Resolved 42 valid Assembly Accessions
2024-01-15 14:32:25 - mutation_scan.genome_extractor.entrez_handler - ERROR - Download failed after 3 retries for GCF_000005845.2
```

## Usage Examples

### Example 1: Search and Download

```python
from mutation_scan.genome_extractor import NCBIDatasetsGenomeDownloader
from pathlib import Path

# Initialize
downloader = NCBIDatasetsGenomeDownloader(
    email="researcher@university.edu",
    output_dir=Path("data/genomes")
)

# Search for genomes
query = "Staphylococcus aureus AND Methicillin Resistant"
accessions = downloader.search_accessions(query, max_results=50)
print(f"Found {len(accessions)} MRSA genomes")

# Download
successful, failed = downloader.download_batch(accessions)
print(f"Downloaded {successful}, failed {failed}")
```

### Example 2: Download from File

```python
# File: accessions.txt
# GCF_000005845.2
# GCF_000007045.1
# GCF_000009065.1

downloader = NCBIDatasetsGenomeDownloader(
    email="researcher@university.edu"
)

accessions = downloader.read_accession_file(Path("accessions.txt"))
successful, failed = downloader.download_batch(accessions)
```

### Example 3: Validate Downloaded Genomes

```python
from mutation_scan.genome_extractor import GenomeProcessor

processor = GenomeProcessor(
    min_coverage=90.0,
    min_length=1000000
)

genome_file = Path("data/genomes/GCF_000005845.2.fasta")
is_valid, message = processor.validate_genome(genome_file)

if is_valid:
    print(f"✓ {message}")
    metadata = processor.extract_metadata(genome_file)
    print(f"  Sequences: {metadata['sequences']}")
else:
    print(f"✗ {message}")
```

## Configuration

Edit `config/config.yaml`:

```yaml
ncbi:
  email: "your.email@example.com"  # REQUIRED
  api_key: null                     # Optional
  batch_size: 50                    # Genomes per batch
  include_plasmids: true            # Include plasmid sequences

genome_extraction:
  input_mode: "search"              # search or file
  input_query: "Escherichia coli"   # For search mode
  input_file: "accessions.txt"      # For file mode
  output_dir: "data/genomes"
  max_genomes: 100
  min_coverage: 90
  min_length: 1000000
```

## Error Handling

The module implements robust error handling:

1. **API Failures**: Automatic retry with exponential backoff (up to 3 attempts)
2. **Zip Processing**: Graceful handling of corrupted downloads
3. **Missing Metadata**: Anti-hallucination (defaults to "N/A")
4. **QC Failures**: Logged explicitly, not silently skipped
5. **Network Issues**: Timeouts and connection errors caught and logged

All errors are logged to `genome_extractor.log` with:
- Timestamp
- Error type
- Specific accession/file affected
- Suggested remediation (when applicable)

## Performance

- **Search**: ~5-10 seconds per query (searches NCBI Assembly database)
- **Download**: ~2-5 seconds per genome (depends on size and network)
- **Batch Processing**: Up to 50 genomes simultaneously (configurable)
- **Metadata Extraction**: Inline during download (no additional overhead)

## Dependencies

- `biopython>=1.79`: Bio.Entrez for searching
- `ncbi-datasets-pylib>=16.0.0`: NCBI Datasets API v2
- `requests>=2.26.0`: HTTP requests for API calls
- `pandas>=1.3.0`: CSV output
- `python-requests-cache>=1.0.0`: Response caching (optional, improves performance)

## Troubleshooting

### Q: "403 Forbidden" error
**A**: Ensure your email is set correctly in the config. NCBI requires a valid email.

### Q: "Empty results for query"
**A**: Check your search query syntax. Use NCBI Assembly database query syntax (https://www.ncbi.nlm.nih.gov/assembly/)

### Q: "QC_FAIL for downloaded genome"
**A**: The sequence coverage is below 90%. Check the log file for details. May indicate corrupted download or partial assembly.

### Q: "No .fna file in zip"
**A**: The NCBI Datasets API response format may have changed. Check log and consider updating the module.

## Future Enhancements

- [ ] Parallel downloads (thread pool)
- [ ] Resume incomplete downloads
- [ ] Support for draft genomes (GCA_)
- [ ] Integration with SRA for raw reads
- [ ] Support for eukaryotic genomes
- [ ] Automatic taxonomy classification

## Support

For issues, questions, or contributions:
- Check `genome_extractor.log` for detailed error messages
- Review examples in `examples/genome_extractor_example.py`
- Consult NCBI Datasets API documentation: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/
