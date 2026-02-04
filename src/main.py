#!/usr/bin/env python3
"""
MutationScan Pipeline Orchestrator (main.py)

This script executes the end-to-end bacterial genome analysis pipeline:
1. Download genomes from NCBI (GenomeExtractor)
2. Find resistance genes (GeneFinder)
3. Extract and translate sequences (SequenceExtractor)
4. Call variants and predict resistance (VariantCaller)
5. Visualize mutations on 3D structures (PyMOLVisualizer)

Author: Senior Python Software Architect (Bioinformatics Specialization)
Platform: Windows-safe (uses pathlib.Path, shutil.which for dependency checks)
"""

import argparse
import logging
import shutil
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

# Import MutationScan modules
from mutation_scan.genome_extractor import NCBIDatasetsGenomeDownloader
from mutation_scan.gene_finder import GeneFinder
from mutation_scan.sequence_extractor import SequenceExtractor
from mutation_scan.variant_caller import VariantCaller
from mutation_scan.visualizer import PyMOLVisualizer


def print_startup_banner() -> None:
    """
    Display platform-specific startup banner with deployment instructions.
    
    Guides users to use Docker on Windows before attempting native execution.
    """
    import platform
    
    os_name = platform.system()
    is_windows = (os_name == "Windows")
    
    print("\n" + "="*70)
    print("         MUTATIONSCAN v1.0 - Clinical AMR Pipeline")
    print("="*70)
    print(f"Platform Detected: {os_name}\n")
    
    if is_windows:
        print("[WINDOWS DEPLOYMENT GUIDE]")
        print("━"*70)
        print("⚠️  Native Windows execution requires WSL or Docker.\n")
        print("RECOMMENDED: Use Docker for one-click deployment:\n")
        print("  Option 1 - Docker Compose (Easiest):")
        print("    $env:NCBI_EMAIL = 'your.email@example.com'")
        print("    docker-compose up\n")
        print("  Option 2 - Docker Run (Custom Parameters):")
        print("    docker run -v ${PWD}/data:/app/data mutationscan:v1 \\")
        print("      --email your@email.com --query 'E. coli' --limit 10\n")
        print("ALTERNATIVE: Install dependencies via WSL:")
        print("    wsl --install")
        print("    wsl -e sudo apt update")
        print("    wsl -e sudo apt install abricate ncbi-blast+")
        print("━"*70 + "\n")
    else:
        print("[UNIX/LINUX DEPLOYMENT]")
        print("━"*70)
        print("✓ Platform supported for native execution.\n")
        print("Required Dependencies:")
        print("  • ABRicate (resistance gene database)")
        print("  • NCBI BLAST+ (sequence alignment)")
        print("  • PyMOL (optional - 3D visualization)\n")
        print("Install via package manager:")
        print("  Ubuntu/Debian: sudo apt install abricate ncbi-blast+ pymol")
        print("  macOS:         brew install brewsci/bio/abricate blast pymol")
        print("━"*70 + "\n")


def setup_logging(log_file: Path = Path("data/logs/pipeline.log")) -> None:
    """
    Configure logging for the pipeline.
    
    Writes to both console (INFO) and file (DEBUG).
    
    Args:
        log_file: Path to log file
    """
    log_file.parent.mkdir(parents=True, exist_ok=True)
    
    # Root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    
    # Console handler (INFO)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )
    console_handler.setFormatter(console_formatter)
    
    # File handler (DEBUG)
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    file_handler.setFormatter(file_formatter)
    
    root_logger.addHandler(console_handler)
    root_logger.addHandler(file_handler)
    
    logging.info("="*70)
    logging.info("MutationScan Pipeline Started")
    logging.info("="*70)


def check_dependencies() -> dict:
    """
    Verify external tools are installed and accessible.
    
    Uses shutil.which() to check PATH for binaries.
    Fails gracefully with clear instructions on Windows.
    
    Returns:
        Dictionary with tool names and their status
        
    Raises:
        EnvironmentError: If critical dependencies are missing
    """
    logging.info("Checking system dependencies...")
    
    dependencies = {
        'abricate': shutil.which('abricate'),
        'blastn': shutil.which('blastn'),
        'pymol': shutil.which('pymol')
    }
    
    # Check ABRicate (critical for gene finding)
    if not dependencies['abricate']:
        logging.error("ABRicate not found in PATH")
        raise EnvironmentError(
            "\n[CRITICAL] ABRicate not found.\n"
            "On Windows, ABRicate is typically not available natively.\n"
            "SOLUTION: Run MutationScan via Docker:\n"
            "  docker build -t mutationscan .\n"
            "  docker run -v $(pwd)/data:/app/data mutationscan --email YOUR_EMAIL --query \"E. coli\"\n"
            "\nAlternatively, install ABRicate via WSL:\n"
            "  wsl --install\n"
            "  wsl -e sudo apt update\n"
            "  wsl -e sudo apt install abricate\n"
        )
    
    # Check BLAST (critical for gene finding)
    if not dependencies['blastn']:
        logging.warning(
            "BLAST+ not found in PATH. Gene finding will rely solely on ABRicate.\n"
            "Install BLAST+ for housekeeping gene detection:\n"
            "  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"
        )
    
    # Check PyMOL (optional for visualization)
    if not dependencies['pymol']:
        logging.warning(
            "PyMOL not found in PATH. Visualization will be skipped.\n"
            "Install PyMOL for 3D structure visualization:\n"
            "  conda install -c schrodinger pymol-bundle\n"
            "  OR download from: https://pymol.org/2/"
        )
    
    logging.info("Dependency check complete:")
    for tool, path in dependencies.items():
        status = f"[OK] {path}" if path else "[MISSING]"
        logging.info(f"  {tool}: {status}")
    
    return dependencies


def step1_download_genomes(
    email: str,
    api_key: Optional[str],
    query: str,
    limit: int,
    output_dir: Path
) -> int:
    """
    Step 1: Download genomes from NCBI.
    
    Args:
        email: NCBI email (required)
        api_key: Optional NCBI API key
        query: Search query (e.g., "E. coli")
        limit: Maximum genomes to download
        output_dir: Output directory for FASTA files
        
    Returns:
        Number of genomes successfully downloaded
        
    Raises:
        Exception: If download fails
    """
    logging.info("")
    logging.info("="*70)
    logging.info("STEP 1: Download Genomes from NCBI")
    logging.info("="*70)
    
    try:
        downloader = NCBIDatasetsGenomeDownloader(
            email=email,
            api_key=api_key,
            output_dir=output_dir,
            log_file=Path("data/logs/genome_extractor.log")
        )
        
        # Search for accessions
        logging.info(f"Searching NCBI for: '{query}'")
        accessions = downloader.search_accessions(query=query, max_results=limit)
        
        if not accessions:
            raise ValueError(f"No genomes found for query: {query}")
        
        logging.info(f"Found {len(accessions)} accessions")
        
        # Download batch
        logging.info("Downloading genome FASTA files...")
        success_count, fail_count = downloader.download_batch(accessions)
        
        logging.info(f"Download complete: {success_count} successful, {fail_count} failed")
        
        if success_count == 0:
            raise RuntimeError("All genome downloads failed. Check logs.")
        
        return success_count
        
    except Exception as e:
        logging.error(f"Step 1 failed: {type(e).__name__}: {e}")
        raise


def step2_find_genes(genomes_dir: Path) -> pd.DataFrame:
    """
    Step 2: Find resistance genes using ABRicate.
    
    Args:
        genomes_dir: Directory containing genome FASTA files
        
    Returns:
        DataFrame with gene coordinates
        
    Raises:
        Exception: If gene finding fails or returns empty results
    """
    logging.info("")
    logging.info("="*70)
    logging.info("STEP 2: Find Resistance Genes (ABRicate)")
    logging.info("="*70)
    
    try:
        gene_finder = GeneFinder(abricate_db="card")
        
        # Find all genome FASTA files
        genome_files = list(genomes_dir.glob("*.fasta")) + list(genomes_dir.glob("*.fna"))
        
        if not genome_files:
            raise FileNotFoundError(f"No genome FASTA files found in {genomes_dir}")
        
        logging.info(f"Scanning {len(genome_files)} genome files for resistance genes...")
        
        # Concatenate results from all genomes
        all_genes = []
        for genome_file in genome_files:
            logging.info(f"  Scanning: {genome_file.name}")
            genes_df = gene_finder.find_resistance_genes(genome_file)
            
            if not genes_df.empty:
                # Add accession column
                accession = genome_file.stem
                genes_df['Accession'] = accession
                all_genes.append(genes_df)
                logging.info(f"    Found {len(genes_df)} genes")
            else:
                logging.warning(f"    No genes found in {genome_file.name}")
        
        # Combine all results
        if not all_genes:
            raise ValueError("No resistance genes found in any genome. Pipeline cannot continue.")
        
        combined_df = pd.concat(all_genes, ignore_index=True)
        logging.info(f"Total genes found: {len(combined_df)}")
        
        # Validate DataFrame has required columns
        required_cols = ['Gene', 'Contig', 'Start', 'End', 'Strand', 'Accession']
        missing_cols = [col for col in required_cols if col not in combined_df.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        return combined_df
        
    except Exception as e:
        logging.error(f"Step 2 failed: {type(e).__name__}: {e}")
        raise


def step3_extract_sequences(
    genes_df: pd.DataFrame,
    genomes_dir: Path,
    output_dir: Path
) -> int:
    """
    Step 3: Extract and translate gene sequences.
    
    Args:
        genes_df: DataFrame with gene coordinates from Step 2
        genomes_dir: Directory containing genome FASTA files
        output_dir: Output directory for protein FASTA files
        
    Returns:
        Number of sequences successfully extracted
        
    Raises:
        Exception: If extraction fails
    """
    logging.info("")
    logging.info("="*70)
    logging.info("STEP 3: Extract and Translate Sequences")
    logging.info("="*70)
    
    try:
        extractor = SequenceExtractor(genomes_dir=genomes_dir)
        
        total_success = 0
        total_fail = 0
        
        # Group by accession
        for accession, group_df in genes_df.groupby('Accession'):
            logging.info(f"Processing {accession}: {len(group_df)} genes")
            
            success, fail = extractor.extract_sequences(
                genes_df=group_df,
                accession=accession,
                output_dir=output_dir,
                translate=True
            )
            
            total_success += success
            total_fail += fail
            logging.info(f"  Extracted: {success} succeeded, {fail} failed")
        
        logging.info(f"Total extracted: {total_success} proteins")
        
        if total_success == 0:
            raise RuntimeError("No sequences extracted. Pipeline cannot continue.")
        
        return total_success
        
    except Exception as e:
        logging.error(f"Step 3 failed: {type(e).__name__}: {e}")
        raise


def step4_call_variants(
    proteins_dir: Path,
    refs_dir: Path,
    output_csv: Path,
    enable_ml: bool = True,
    ml_models_dir: Optional[Path] = None,
    antibiotic: str = "Ciprofloxacin"
) -> pd.DataFrame:
    """
    Step 4: Call variants and predict resistance.
    
    Args:
        proteins_dir: Directory containing protein FASTA files
        refs_dir: Directory containing wild-type reference sequences
        output_csv: Path to save mutation report CSV
        enable_ml: Whether to enable ML predictions
        ml_models_dir: Path to ML models directory
        antibiotic: Antibiotic name for ML predictions
        
    Returns:
        DataFrame with mutation report
        
    Raises:
        Exception: If variant calling fails
    """
    logging.info("")
    logging.info("="*70)
    logging.info("STEP 4: Call Variants and Predict Resistance")
    logging.info("="*70)
    
    try:
        caller = VariantCaller(
            refs_dir=refs_dir,
            enable_ml=enable_ml,
            ml_models_dir=ml_models_dir,
            antibiotic=antibiotic
        )
        
        # Check if ML models exist
        if enable_ml and ml_models_dir:
            model_files = list(ml_models_dir.glob("*.pkl"))
            if model_files:
                logging.info(f"ML models found: {len(model_files)} file(s)")
            else:
                logging.warning("ML enabled but no .pkl models found in models/")
        
        # Find all protein FASTA files
        protein_files = list(proteins_dir.glob("*.faa"))
        
        if not protein_files:
            raise FileNotFoundError(f"No protein FASTA files found in {proteins_dir}")
        
        logging.info(f"Calling variants for {len(protein_files)} protein files...")
        
        # Call variants
        mutations_df = caller.call_variants(
            query_dir=proteins_dir,
            output_csv=output_csv
        )
        
        if mutations_df.empty:
            logging.warning("No mutations detected")
        else:
            logging.info(f"Detected {len(mutations_df)} mutations")
            
            # Count by status
            if 'Status' in mutations_df.columns:
                status_counts = mutations_df['Status'].value_counts()
                for status, count in status_counts.items():
                    logging.info(f"  {status}: {count}")
        
        logging.info(f"Mutation report saved to: {output_csv}")
        
        return mutations_df
        
    except Exception as e:
        logging.error(f"Step 4 failed: {type(e).__name__}: {e}")
        raise


def step5_visualize_mutations(
    mutation_csv: Path,
    output_dir: Path,
    pymol_available: bool
) -> int:
    """
    Step 5: Visualize mutations on 3D structures (optional).
    
    Args:
        mutation_csv: Path to mutation report CSV
        output_dir: Output directory for PNG images
        pymol_available: Whether PyMOL is installed
        
    Returns:
        Number of visualizations created
        
    Raises:
        Exception: If visualization fails
    """
    logging.info("")
    logging.info("="*70)
    logging.info("STEP 5: Visualize Mutations (3D Structures)")
    logging.info("="*70)
    
    if not pymol_available:
        logging.warning("PyMOL not available. Skipping visualization.")
        return 0
    
    try:
        visualizer = PyMOLVisualizer(output_dir=output_dir)
        
        # Generate visualizations
        result = visualizer.visualize_mutations(
            mutation_csv=mutation_csv,
            filter_status=["Resistant", "Predicted High Risk", "VUS"]
        )
        
        total_images = sum(len(images) for images in result.values())
        logging.info(f"Generated {total_images} visualization images")
        
        return total_images
        
    except Exception as e:
        logging.error(f"Step 5 failed: {type(e).__name__}: {e}")
        logging.warning("Visualization failed but pipeline will continue")
        return 0


def main():
    """
    Main entry point for MutationScan pipeline.
    """
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="MutationScan: End-to-end bacterial genome analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage (download 10 E. coli genomes)
  python main.py --email your.email@example.com
  
  # Custom query with visualization
  python main.py --email your.email@example.com --query "Klebsiella pneumoniae" --limit 5 --visualize
  
  # With NCBI API key for faster downloads
  python main.py --email your.email@example.com --api-key YOUR_API_KEY --limit 20
        """
    )
    
    parser.add_argument(
        '--email',
        required=True,
        help='Email address (required by NCBI policy)'
    )
    
    parser.add_argument(
        '--api-key',
        default=None,
        help='NCBI API key (optional, for faster downloads)'
    )
    
    parser.add_argument(
        '--query',
        default='Escherichia coli',
        help='NCBI search query (default: "Escherichia coli")'
    )
    
    parser.add_argument(
        '--limit',
        type=int,
        default=10,
        help='Maximum number of genomes to download (default: 10)'
    )
    
    parser.add_argument(
        '--visualize',
        action='store_true',
        help='Generate 3D structure visualizations (requires PyMOL)'
    )
    
    parser.add_argument(
        '--no-ml',
        action='store_true',
        help='Disable ML predictions for unknown mutations'
    )
    
    parser.add_argument(
        '--antibiotic',
        default='Ciprofloxacin',
        help='Antibiotic name for ML predictions (default: Ciprofloxacin)'
    )
    
    args = parser.parse_args()
    
    # Show startup banner with platform-specific instructions
    print_startup_banner()
    
    # Setup logging
    setup_logging()
    
    # Log configuration
    logging.info(f"Configuration:")
    logging.info(f"  Email: {args.email}")
    logging.info(f"  Query: {args.query}")
    logging.info(f"  Limit: {args.limit}")
    logging.info(f"  Visualize: {args.visualize}")
    logging.info(f"  ML Enabled: {not args.no_ml}")
    logging.info(f"  Antibiotic: {args.antibiotic}")
    
    # Check dependencies
    try:
        deps = check_dependencies()
    except EnvironmentError as e:
        logging.error(str(e))
        sys.exit(1)
    
    # Define paths (Windows-safe with pathlib)
    genomes_dir = Path("data") / "genomes"
    proteins_dir = Path("data") / "proteins"
    refs_dir = Path("data") / "refs"
    results_dir = Path("data") / "results"
    viz_dir = results_dir / "visualizations"
    ml_models_dir = Path("models")
    mutation_csv = results_dir / "mutation_report.csv"
    
    # Create directories
    for directory in [genomes_dir, proteins_dir, refs_dir, results_dir, viz_dir]:
        directory.mkdir(parents=True, exist_ok=True)
    
    try:
        # STEP 1: Download genomes
        genome_count = step1_download_genomes(
            email=args.email,
            api_key=args.api_key,
            query=args.query,
            limit=args.limit,
            output_dir=genomes_dir
        )
        
        # STEP 2: Find genes
        genes_df = step2_find_genes(genomes_dir=genomes_dir)
        
        # STEP 3: Extract sequences
        seq_count = step3_extract_sequences(
            genes_df=genes_df,
            genomes_dir=genomes_dir,
            output_dir=proteins_dir
        )
        
        # STEP 4: Call variants
        mutations_df = step4_call_variants(
            proteins_dir=proteins_dir,
            refs_dir=refs_dir,
            output_csv=mutation_csv,
            enable_ml=not args.no_ml,
            ml_models_dir=ml_models_dir if (not args.no_ml) else None,
            antibiotic=args.antibiotic
        )
        
        # STEP 5: Visualize (optional)
        if args.visualize:
            viz_count = step5_visualize_mutations(
                mutation_csv=mutation_csv,
                output_dir=viz_dir,
                pymol_available=deps['pymol'] is not None
            )
        
        # Pipeline summary
        logging.info("")
        logging.info("="*70)
        logging.info("PIPELINE COMPLETE")
        logging.info("="*70)
        logging.info(f"Genomes downloaded: {genome_count}")
        logging.info(f"Genes found: {len(genes_df)}")
        logging.info(f"Sequences extracted: {seq_count}")
        logging.info(f"Mutations detected: {len(mutations_df)}")
        if args.visualize:
            logging.info(f"Visualizations created: {viz_count}")
        logging.info(f"Results saved to: {results_dir}")
        logging.info("="*70)
        
    except Exception as e:
        logging.error("")
        logging.error("="*70)
        logging.error("PIPELINE FAILED")
        logging.error("="*70)
        logging.error(f"Error: {type(e).__name__}: {e}")
        logging.error("Check logs for details: data/logs/pipeline.log")
        logging.error("="*70)
        sys.exit(1)


if __name__ == "__main__":
    main()
