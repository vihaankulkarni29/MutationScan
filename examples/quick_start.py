"""
Quick-start guide for the refactored GenomeExtractor module.

Copy this file and modify it for your specific use case.
"""

from pathlib import Path
from mutation_scan.core import NCBIDatasetsGenomeDownloader, GenomeProcessor

# ============================================================================
# CONFIGURATION
# ============================================================================

# 1. Set your NCBI email (REQUIRED - NCBI policy)
NCBI_EMAIL = "your.email@example.com"

# 2. Optional: Set your NCBI API key (faster requests)
NCBI_API_KEY = None

# 3. Set output directory
OUTPUT_DIR = Path("data/genomes")
LOG_FILE = Path("data/logs/genome_extractor.log")

# ============================================================================
# MODE A: QUERY-BASED SEARCH
# ============================================================================

def download_by_query(query: str, max_results: int = 10):
    """
    Download genomes by search query.
    
    Example queries:
    - "Escherichia coli"
    - "Staphylococcus aureus AND Methicillin Resistant"
    - "Mycobacterium tuberculosis AND Drug Resistant"
    - "Salmonella enterica AND Antibiotic"
    """
    downloader = NCBIDatasetsGenomeDownloader(
        email=NCBI_EMAIL,
        api_key=NCBI_API_KEY,
        output_dir=OUTPUT_DIR,
        log_file=LOG_FILE,
    )
    
    print(f"Searching for: {query}")
    accessions = downloader.search_accessions(query, max_results=max_results)
    print(f"Found {len(accessions)} genomes: {accessions[:5]}...")
    
    if accessions:
        print(f"Downloading {len(accessions)} genomes...")
        successful, failed = downloader.download_batch(accessions)
        print(f"✓ Successfully downloaded: {successful}")
        print(f"✗ Failed downloads: {failed}")
        
        # Check metadata
        metadata_file = OUTPUT_DIR / "metadata_master.csv"
        if metadata_file.exists():
            print(f"✓ Metadata saved to: {metadata_file}")


# ============================================================================
# MODE B: FILE-BASED DOWNLOAD
# ============================================================================

def download_from_file(accession_file: Path):
    """
    Download genomes from a list file.
    
    File format (accessions.txt):
    GCF_000005845.2
    GCF_000007045.1
    GCF_000009065.1
    
    Lines starting with # are comments
    Empty lines are skipped
    """
    downloader = NCBIDatasetsGenomeDownloader(
        email=NCBI_EMAIL,
        api_key=NCBI_API_KEY,
        output_dir=OUTPUT_DIR,
        log_file=LOG_FILE,
    )
    
    print(f"Reading accessions from: {accession_file}")
    accessions = downloader.read_accession_file(accession_file)
    print(f"Loaded {len(accessions)} accessions")
    
    print(f"Downloading {len(accessions)} genomes...")
    successful, failed = downloader.download_batch(accessions)
    print(f"✓ Successfully downloaded: {successful}")
    print(f"✗ Failed downloads: {failed}")


# ============================================================================
# VALIDATE DOWNLOADED GENOMES
# ============================================================================

def validate_genomes():
    """
    Validate and QC-check downloaded genomes.
    """
    processor = GenomeProcessor(
        min_coverage=90.0,    # 90% minimum coverage
        min_length=1000000,   # 1MB minimum length
    )
    
    fasta_files = list(OUTPUT_DIR.glob("*.fasta"))
    print(f"Found {len(fasta_files)} FASTA files")
    
    passed = 0
    failed = 0
    
    for fasta_file in fasta_files:
        is_valid, message = processor.validate_genome(fasta_file)
        
        if is_valid:
            passed += 1
            print(f"✓ {fasta_file.name}: {message}")
            
            # Extract metadata
            metadata = processor.extract_metadata(fasta_file)
            print(f"  Sequences: {metadata['sequences']}")
        else:
            failed += 1
            print(f"✗ {fasta_file.name}: {message}")
    
    print(f"\nValidation Summary: {passed} PASS, {failed} FAIL")


# ============================================================================
# INSPECT METADATA
# ============================================================================

def show_metadata():
    """
    Display extracted metadata in CSV format.
    """
    metadata_file = OUTPUT_DIR / "metadata_master.csv"
    
    if not metadata_file.exists():
        print(f"Metadata file not found: {metadata_file}")
        print("Run download first to generate metadata.")
        return
    
    try:
        import pandas as pd
        df = pd.read_csv(metadata_file)
        print(f"\nMetadata for {len(df)} genomes:")
        print(df.to_string())
        print(f"\nMetadata file: {metadata_file}")
    except ImportError:
        print("pandas not installed. Cannot display metadata.")


# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("MutationScan GenomeExtractor - Quick Start")
    print("=" * 70)
    
    # STEP 1: Update your NCBI email at the top of this file
    if NCBI_EMAIL == "your.email@example.com":
        print("\n⚠️  WARNING: Please set your NCBI email before running!")
        print("   Edit NCBI_EMAIL = 'your.email@example.com' at the top of this file")
        exit(1)
    
    # STEP 2: Create data directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_FILE.parent.mkdir(parents=True, exist_ok=True)
    print(f"\n✓ Output directory: {OUTPUT_DIR}")
    print(f"✓ Log file: {LOG_FILE}")
    
    # STEP 3: Choose your download mode
    print("\n" + "=" * 70)
    print("SELECT MODE:")
    print("=" * 70)
    print("1. Query-based search (recommended for beginners)")
    print("2. File-based batch (for large lists)")
    print("3. Validate existing genomes")
    print("4. Show metadata")
    print("0. Exit")
    
    choice = input("\nEnter choice (0-4): ").strip()
    
    if choice == "1":
        # Example query
        query = input("Enter search query (default: 'Escherichia coli'): ").strip()
        if not query:
            query = "Escherichia coli"
        max_results = int(input("Max results (default: 10): ").strip() or "10")
        download_by_query(query, max_results)
        
    elif choice == "2":
        accession_file = Path(input("Enter accession file path: ").strip())
        if accession_file.exists():
            download_from_file(accession_file)
        else:
            print(f"File not found: {accession_file}")
            
    elif choice == "3":
        validate_genomes()
        
    elif choice == "4":
        show_metadata()
        
    else:
        print("Exiting...")
    
    print("\n" + "=" * 70)
    print("For more examples, see: examples/genome_extractor_example.py")
    print("For full API documentation, see: docs/GENOME_EXTRACTOR_API.md")
    print("=" * 70)
