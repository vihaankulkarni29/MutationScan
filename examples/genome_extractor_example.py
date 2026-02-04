"""
Example usage of the refactored GenomeExtractor with NCBI Datasets API v2.

This script demonstrates:
- Mode A: Query-based genome search and download
- Mode B: Batch download from accession file
- Metadata extraction
- Error handling
"""

import sys
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mutation_scan.core import NCBIDatasetsGenomeDownloader, GenomeProcessor


def example_search_mode():
    """
    Example 1: Search-based mode
    Query NCBI for genomes and download them.
    """
    print("=" * 60)
    print("EXAMPLE 1: Search Mode (Query-based)")
    print("=" * 60)
    
    # Initialize downloader with your NCBI email
    downloader = NCBIDatasetsGenomeDownloader(
        email="your.email@example.com",  # REQUIRED: Your email
        api_key=None,  # Optional: Your NCBI API key
        output_dir=Path("data/genomes/search_mode"),
        log_file=Path("data/logs/genome_extractor_search.log"),
    )
    
    # Example query: Antibiotic-resistant E. coli
    query = 'Escherichia coli AND "antibiotic resistant"'
    
    try:
        print(f"\nSearching for: {query}")
        accessions = downloader.search_accessions(query, max_results=10)
        print(f"Found {len(accessions)} genomes: {accessions[:3]}...")
        
        if accessions:
            print(f"\nDownloading {len(accessions)} genomes...")
            successful, failed = downloader.download_batch(accessions)
            print(f"✓ Downloaded: {successful} | ✗ Failed: {failed}")
            
    except Exception as e:
        print(f"ERROR: {e}")


def example_file_mode():
    """
    Example 2: File-based mode
    Download genomes from a list of accessions.
    """
    print("\n" + "=" * 60)
    print("EXAMPLE 2: File Mode (Accession list)")
    print("=" * 60)
    
    # Create sample accession file
    accession_file = Path("data/accessions.txt")
    accession_file.parent.mkdir(parents=True, exist_ok=True)
    
    sample_accessions = [
        "GCF_000005845.2",  # E. coli K-12
        "GCF_000007045.1",  # Helicobacter pylori
        "GCF_000009065.1",  # Vibrio cholerae
    ]
    
    with open(accession_file, 'w') as f:
        f.write("\n".join(sample_accessions))
    
    print(f"Created accession file: {accession_file}")
    print(f"Accessions: {sample_accessions}")
    
    # Initialize downloader
    downloader = NCBIDatasetsGenomeDownloader(
        email="your.email@example.com",  # REQUIRED: Your email
        api_key=None,
        output_dir=Path("data/genomes/file_mode"),
        log_file=Path("data/logs/genome_extractor_file.log"),
    )
    
    try:
        print("\nReading accessions from file...")
        accessions = downloader.read_accession_file(accession_file)
        print(f"Loaded {len(accessions)} accessions")
        
        print(f"\nDownloading {len(accessions)} genomes...")
        successful, failed = downloader.download_batch(accessions)
        print(f"✓ Downloaded: {successful} | ✗ Failed: {failed}")
        
    except Exception as e:
        print(f"ERROR: {e}")


def example_validate_genomes():
    """
    Example 3: Validate and process downloaded genomes
    """
    print("\n" + "=" * 60)
    print("EXAMPLE 3: Genome Validation & QC")
    print("=" * 60)
    
    processor = GenomeProcessor(min_coverage=90.0, min_length=1000000)
    
    genome_dir = Path("data/genomes/file_mode")
    fasta_files = list(genome_dir.glob("*.fasta"))
    
    if not fasta_files:
        print("No FASTA files found. Run Example 2 first.")
        return
    
    print(f"\nValidating {len(fasta_files)} genomes...")
    
    for fasta_file in fasta_files[:3]:  # Check first 3
        is_valid, message = processor.validate_genome(fasta_file)
        status = "✓" if is_valid else "✗"
        print(f"{status} {fasta_file.name}: {message}")
        
        if is_valid:
            metadata = processor.extract_metadata(fasta_file)
            print(f"   - Sequences: {metadata['sequences']}")
            print(f"   - File size: {metadata['file_size']} bytes")


def example_metadata_inspection():
    """
    Example 4: Inspect extracted metadata
    """
    print("\n" + "=" * 60)
    print("EXAMPLE 4: Metadata Inspection")
    print("=" * 60)
    
    metadata_file = Path("data/genomes/file_mode/metadata_master.csv")
    
    if not metadata_file.exists():
        print(f"Metadata file not found: {metadata_file}")
        print("Run Example 2 first to generate metadata.")
        return
    
    try:
        import pandas as pd
        df = pd.read_csv(metadata_file)
        print(f"\nLoaded metadata for {len(df)} genomes:")
        print(df.to_string())
    except ImportError:
        print("pandas not installed. Cannot read CSV.")


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("MutationScan GenomeExtractor Examples")
    print("NCBI Datasets API v2")
    print("=" * 60)
    
    # Uncomment the examples you want to run
    
    # Example 1: Query-based search
    # example_search_mode()
    
    # Example 2: File-based download
    # example_file_mode()
    
    # Example 3: Validate genomes
    # example_validate_genomes()
    
    # Example 4: Inspect metadata
    # example_metadata_inspection()
    
    print("\n" + "=" * 60)
    print("QUICK START:")
    print("=" * 60)
    print("""
    1. Edit your NCBI email in the examples above
    2. Uncomment the examples you want to run
    3. Run: python examples/genome_extractor_example.py
    
    IMPORTANT:
    - Set your NCBI email (required by NCBI)
    - Optionally add your NCBI API key for faster requests
    - Ensure internet connection for NCBI API access
    - Output will be saved to data/genomes/
    - Logs will be saved to data/logs/genome_extractor.log
    """)
