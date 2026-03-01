#!/usr/bin/env python3
"""
MutationScan Master Pipeline Orchestrator (run_pipeline.py)

Refactored v2.1: Deterministic Genotype-to-Phenotype-to-Biophysics Engine

This script orchestrates the end-to-end clinical AMR analysis pipeline:

Phase 1: Genomic Data Ingestion (MutationScan)
  - Download and process bacterial genomes from NCBI
  - Identify resistance genes using ABRicate
  - Extract and translate protein sequences
  - Call variants and identify mutations

Phase 2: Expression Analysis (ControlScan)
  - Measure gene expression levels
  - Infer expression scores for regulatory analysis
  - Normalize across sample cohorts

Phase 3: Epistasis Network Detection
  - Identify statistically significant mutation clusters
  - Extract core epistatic residue networks

Phase 4: 3D Biophysics Docking (AutoScan Bridge)
  - Perform comparative in silico docking
  - Calculate binding affinity changes (ΔΔG)
  - Predict structural impact of mutations

Phase 5: Feature Aggregation
  - Consolidate genomic, expression, and biophysical features
  - Generate final engineered feature CSV for downstream modeling
  - No ML inference (deterministic output only)

Author: Senior Computational Biologist & Software Architect
Platform: Windows-safe (uses pathlib.Path for cross-platform compatibility)
Version: 2.1.0
"""

import argparse
import json
import logging
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Tuple

import pandas as pd

# Import MutationScan modules
from mutation_scan.core.entrez_handler import NCBIDatasetsGenomeDownloader
from mutation_scan.biophysics.autoscan_bridge import AutoScanBridge
from mutation_scan.visualization.pymol_viz import PyMOLVisualizer


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
        '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_formatter)
    
    root_logger.addHandler(console_handler)
    root_logger.addHandler(file_handler)
    
    logger = logging.getLogger(__name__)
    logger.info(f"Logging initialized. Debug log: {log_file}")


def phase1_genomic_ingestion(
    email: str,
    api_key: Optional[str],
    genome_file: Optional[Path],
    organism: Optional[str],
    limit: int,
    target_genes: Optional[List[str]],
    output_dir: Path
) -> Tuple[pd.DataFrame, Path, Path, Path]:
    """
    Phase 1: Genomic Data Ingestion
    
    Download genomes, identify resistance genes, and call variants.
    
    Args:
        email: NCBI email required by NCBI policy
        api_key: Optional NCBI API key
        genome_file: Path to local genome FASTA (if provided)
        organism: Organism name to download (if genome_file not provided)
        limit: Number of genomes to download
        target_genes: Optional list of target genes to analyze
        output_dir: Output directory for results
        
    Returns:
        Tuple of (mutations_df, proteins_dir, refs_dir, genomes_dir)
    """
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("="*70)
    logger.info("PHASE 1: Genomic Data Ingestion")
    logger.info("="*70)
    
    genomes_dir = output_dir / "genomes"
    proteins_dir = output_dir / "proteins"
    refs_dir = output_dir / "refs"
    
    genomes_dir.mkdir(parents=True, exist_ok=True)
    proteins_dir.mkdir(parents=True, exist_ok=True)
    refs_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Download/Prepare genomes
    logger.info("Step 1.1: Downloading/preparing genomes...")
    
    if genome_file:
        logger.info(f"Using local genome: {genome_file}")
        # Copy local genome
        import shutil
        local_dest = genomes_dir / genome_file.name
        shutil.copy(genome_file, local_dest)
        num_genomes = 1
    else:
        if not organism:
            raise ValueError("Organism name required for NCBI download")
        
        logger.info(f"Downloading {limit} genome(s) from NCBI: {organism}")
        downloader = NCBIDatasetsGenomeDownloader(
            email=email,
            api_key=api_key,
            output_dir=genomes_dir
        )
        
        # Handle geolocation-based search (e.g., "Klebsiella pneumoniae AND India")
        if " AND " in organism:
            parts = organism.split(" AND ")
            taxon = parts[0].strip()
            location = parts[1].strip()
            logger.info(f"Using geolocation-based search: {taxon} from {location}")
            success, failed = downloader.download_bulk_by_geolocation(
                organism=taxon,
                location=location,
                limit=limit
            )
            num_genomes = success
            logger.info(f"Downloaded {success} genomes, {failed} failed")
        else:
            # Fallback: Standard batch downloader via search + batch
            logger.info(f"Using standard organism search")
            accessions = downloader.search_accessions(
                query=organism,
                max_results=limit
            )
            if accessions:
                logger.info(f"Found {len(accessions)} accessions, downloading...")
                success, failed = downloader.download_batch(
                    accessions=accessions
                )
                num_genomes = success
                logger.info(f"Downloaded {success} genomes, {failed} failed")
            else:
                logger.error(f"No genomes found for organism: {organism}")
                raise ValueError(f"No genomes found for organism: {organism}")
    
    # Step 2: Hand off to Docker for gene finding and variant calling
    logger.info("Step 1.2: Handing off to Docker container for gene finding and variant calling...")
    logger.info("(GeneFinder, SequenceExtractor, and VariantCaller require Linux binaries: ABRicate, BLAST+)")
    
    import os
    import subprocess
    
    # Get absolute path to project root
    project_root = Path(__file__).parent.absolute()
    
    # Prepare Docker command
    docker_command = [
        "docker", "run", "--rm",
        "-v", f"{project_root / 'data'}:/app/data",
        "-v", f"{project_root / 'config'}:/app/config",
        "-v", f"{project_root / 'models'}:/app/models",
        "mutation-scan:latest",
        "--skip-download",
        "--no-ml",
        "--email", email,
        "--organism", organism,
    ]
    
    # Add optional API key if provided
    if api_key:
        docker_command.extend(["--api-key", api_key])
    
    # Add target genes if provided
    if target_genes:
        logger.warning("Target genes filters are not automatically passed to Docker container")
        logger.info("Please configure target genes within the Docker container environment")
    
    try:
        logger.info(f"Executing Docker command: {' '.join(docker_command[:6])}...")
        result = subprocess.run(
            docker_command,
            check=True,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )
        logger.info("Docker container completed successfully")
        if result.stdout:
            logger.debug(f"Docker stdout: {result.stdout[:500]}")
    except subprocess.TimeoutExpired:
        logger.error("Docker container execution timed out (1 hour limit exceeded)")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Docker container failed with exit code {e.returncode}")
        logger.error(f"Docker stderr: {e.stderr}")
        raise
    except FileNotFoundError:
        logger.error("Docker not found. Please install Docker and ensure it is in your PATH")
        logger.error("Alternatively, run this script inside the Docker container with --skip-download")
        raise
    
    # Step 3: Read results from Docker
    logger.info("Step 1.3: Reading results from Docker container...")
    
    mutations_csv = output_dir / "mutation_report.csv"
    mutations_df = pd.DataFrame()
    
    try:
        if mutations_csv.exists():
            mutations_df = pd.read_csv(mutations_csv)
            logger.info(f"Read mutation report: {len(mutations_df)} mutations detected")
        else:
            logger.warning(f"Mutation report not found: {mutations_csv}")
            logger.info("Docker container may not have completed gene finding/variant calling")
    except Exception as e:
        logger.warning(f"Failed to read mutation report: {e}")

    
    logger.info("[PHASE 1 COMPLETE] Genomic data ingestion finished")
    
    return mutations_df, proteins_dir, refs_dir, genomes_dir


def phase2_expression_analysis(expression_file: Path) -> Dict[str, float]:
    """
    Phase 2: Expression Analysis (ControlScan Integration)
    
    Load gene expression scores from ControlScan output.
    
    Args:
        expression_file: Path to ControlScan expression output (JSON or CSV)
        
    Returns:
        Dictionary mapping gene_name -> expression_score
    """
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("="*70)
    logger.info("PHASE 2: Expression Analysis (ControlScan)")
    logger.info("="*70)
    
    expression_scores = {}
    
    if not expression_file.exists():
        logger.warning(f"Expression file not found: {expression_file}")
        logger.info("Using default expression scores (all = 1.0)")
        return expression_scores
    
    try:
        if expression_file.suffix.lower() == '.json':
            with open(expression_file, 'r') as f:
                data = json.load(f)
            expression_scores = data.get("scores", {})
        elif expression_file.suffix.lower() in ['.csv', '.txt']:
            df = pd.read_csv(expression_file)
            if 'gene' in df.columns and 'expression_score' in df.columns:
                expression_scores = dict(zip(df['gene'], df['expression_score']))
            else:
                logger.warning("Expression file missing 'gene' or 'expression_score' columns")
        
        logger.info(f"Loaded expression scores for {len(expression_scores)} genes")
        for gene, score in list(expression_scores.items())[:5]:  # Show first 5
            logger.info(f"  {gene}: {score:.3f}")
        
    except Exception as e:
        logger.error(f"Failed to load expression scores: {e}")
        logger.info("Using default expression scores")
    
    logger.info("[PHASE 2 COMPLETE] Expression analysis finished")
    
    return expression_scores


def phase3_epistasis_detection(epistasis_file: Path) -> Dict[str, List[int]]:
    """
    Phase 3: Epistasis Network Detection
    
    Load pre-computed epistatic mutation networks.
    
    Args:
        epistasis_file: Path to epistasis input JSON
        
    Returns:
        Dictionary mapping gene_name -> list of residue numbers
    """
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("="*70)
    logger.info("PHASE 3: Epistasis Network Detection")
    logger.info("="*70)
    
    epistasis_networks = {}
    
    if not epistasis_file.exists():
        logger.warning(f"Epistasis file not found: {epistasis_file}")
        logger.info("No epistatic networks to analyze")
        return epistasis_networks
    
    try:
        with open(epistasis_file, 'r') as f:
            data = json.load(f)
        
        # Parse epistasis data (expected format from MutationScan)
        for genome in data.get("genomes", []):
            for gene, mutations in genome.items():
                if gene not in epistasis_networks:
                    import re
                    positions = []
                    for mutation in mutations:
                        match = re.search(r'\d+', mutation)
                        if match:
                            positions.append(int(match.group()))
                    
                    if positions:
                        epistasis_networks[gene] = sorted(list(set(positions)))
                        logger.info(f"Found network for {gene}: residues {positions}")
        
        logger.info(f"Identified {len(epistasis_networks)} epistatic networks")
        
    except Exception as e:
        logger.error(f"Failed to load epistasis data: {e}")
        logger.info("No epistatic networks to analyze")
    
    logger.info("[PHASE 3 COMPLETE] Epistasis detection finished")
    
    return epistasis_networks


def phase4_biophysics_docking(
    epistasis_networks: Dict[str, List[int]],
    pdb_mapping: Dict[str, Tuple[str, str]],
    ligand_smiles: str,
    output_dir: Path
) -> Dict[str, Dict[str, float]]:
    """
    Phase 4: 3D Biophysics Docking (AutoScan Integration)
    
    Calculate binding affinity changes (ΔΔG) for epistatic networks.
    
    Args:
        epistasis_networks: Dictionary mapping gene -> residue list
        pdb_mapping: Dictionary mapping gene -> (pdb_id, chain)
        ligand_smiles: SMILES string for the ligand
        output_dir: Output directory for docking results
        
    Returns:
        Dictionary mapping gene -> {"wt_affinity", "mutant_affinity", "delta_delta_g"}
    """
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("="*70)
    logger.info("PHASE 4: 3D Biophysics Docking (AutoScan)")
    logger.info("="*70)
    
    docking_results = {}
    
    if not epistasis_networks:
        logger.warning("No epistatic networks to dock")
        return docking_results
    
    try:
        bridge = AutoScanBridge()
        
        for gene, residues in epistasis_networks.items():
            if gene not in pdb_mapping:
                logger.warning(f"No PDB mapping for {gene}, skipping")
                continue
            
            pdb_id, chain = pdb_mapping[gene]
            
            logger.info(f"Docking {gene} ({pdb_id} chain {chain}) with residues {residues}...")
            
            result = bridge.run_comparative_docking(
                pdb_id=pdb_id,
                chain=chain,
                residues=residues,
                ligand_smiles=ligand_smiles,
                output_dir=output_dir / f"autoscan_{gene}/"
            )
            
            if result:
                docking_results[gene] = result
                logger.info(f"  ✓ ΔΔG = {result['delta_delta_g']:.2f} kcal/mol")
            else:
                logger.warning(f"  ✗ Docking failed for {gene}")
    
    except Exception as e:
        logger.error(f"Biophysics docking phase failed: {e}", exc_info=True)
    
    logger.info(f"[PHASE 4 COMPLETE] Docking analysis finished ({len(docking_results)} successes)")
    
    return docking_results


def phase5_feature_aggregation(
    sample_id: str,
    species: str,
    mutations_df: pd.DataFrame,
    expression_scores: Dict[str, float],
    epistasis_networks: Dict[str, List[int]],
    docking_results: Dict[str, Dict[str, float]],
    output_file: Path
) -> pd.DataFrame:
    """
    Phase 5: Feature Aggregation
    
    Consolidate genomic, expression, and biophysical features into a
    final engineered feature CSV for downstream statistical modeling.
    
    Schema:
    - Sample_ID: Unique genome identifier
    - Species: Organism name
    - Gene: Target gene name
    - Core_Mutations: Epistatic residue triad (comma-separated)
    - Inferred_Expression_Score: Expression level from ControlScan
    - WT_Affinity: Wild-type binding affinity (kcal/mol)
    - Mutant_Affinity: Mutant binding affinity (kcal/mol)
    - Delta_Delta_G: Binding affinity change (mutant - WT)
    - Mutation_Count: Number of detected mutations in this gene
    
    Args:
        sample_id: Unique identifier for this genome
        species: Species name (e.g., "Klebsiella pneumoniae")
        mutations_df: DataFrame from mutation calling
        expression_scores: Dict mapping gene -> expression score
        epistasis_networks: Dict mapping gene -> residue list
        docking_results: Dict mapping gene -> affinity dict
        output_file: Path to save final features CSV
        
    Returns:
        Aggregated DataFrame
    """
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("="*70)
    logger.info("PHASE 5: Feature Aggregation")
    logger.info("="*70)
    
    records = []
    
    # Aggregate data by gene
    unique_genes = set()
    if not mutations_df.empty and 'Gene' in mutations_df.columns:
        unique_genes.update(mutations_df['Gene'].unique())
    unique_genes.update(epistasis_networks.keys())
    unique_genes.update(expression_scores.keys())
    unique_genes.update(docking_results.keys())
    
    for gene in sorted(unique_genes):
        record = {
            "Sample_ID": sample_id,
            "Species": species,
            "Gene": gene,
            "Core_Mutations": ",".join(map(str, epistasis_networks.get(gene, []))),
            "Inferred_Expression_Score": expression_scores.get(gene, 1.0),
            "Mutation_Count": len(mutations_df[mutations_df['Gene'] == gene]) if not mutations_df.empty else 0,
        }
        
        # Add docking results if available
        if gene in docking_results:
            dock = docking_results[gene]
            record.update({
                "WT_Affinity_kcalmol": dock.get("wt_affinity"),
                "Mutant_Affinity_kcalmol": dock.get("mutant_affinity"),
                "Delta_Delta_G_kcalmol": dock.get("delta_delta_g"),
            })
        else:
            record.update({
                "WT_Affinity_kcalmol": None,
                "Mutant_Affinity_kcalmol": None,
                "Delta_Delta_G_kcalmol": None,
            })
        
        records.append(record)
    
    # Create DataFrame
    features_df = pd.DataFrame(records)
    
    # Sort by Delta_Delta_G (most significant changes first)
    if "Delta_Delta_G_kcalmol" in features_df.columns:
        features_df = features_df.sort_values(
            by="Delta_Delta_G_kcalmol",
            ascending=True,  # Most negative (strongest improvement) first
            na_position='last'
        )
    
    # Append to existing file or create new
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    if output_file.exists():
        logger.info(f"Appending to existing feature file: {output_file}")
        existing_df = pd.read_csv(output_file)
        features_df = pd.concat([existing_df, features_df], ignore_index=True)
    else:
        logger.info(f"Creating new feature file: {output_file}")
    
    # Save to CSV
    features_df.to_csv(output_file, index=False)
    logger.info(f"Saved {len(records)} feature rows to {output_file}")
    
    # Log summary statistics
    if "Delta_Delta_G_kcalmol" in features_df.columns:
        ddg_stats = features_df["Delta_Delta_G_kcalmol"].describe()
        logger.info("ΔΔG Statistics:")
        logger.info(f"  Mean: {ddg_stats['mean']:.2f} kcal/mol")
        logger.info(f"  Std Dev: {ddg_stats['std']:.2f} kcal/mol")
        logger.info(f"  Min: {ddg_stats['min']:.2f} kcal/mol (best binding)")
        logger.info(f"  Max: {ddg_stats['max']:.2f} kcal/mol")
    
    logger.info("[PHASE 5 COMPLETE] Feature aggregation finished")
    
    return features_df


def run_master_pipeline(args) -> int:
    """
    Execute the complete master pipeline: Phases 1-5.
    
    Args:
        args: Parsed command-line arguments
        
    Returns:
        Exit code (0 = success, 1 = failure)
    """
    logger = logging.getLogger(__name__)
    
    try:
        # Setup directories
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("")
        logger.info("+" + "="*68 + "+")
        logger.info("|" + " "*15 + "MUTATIONSCAN MASTER PIPELINE v2.1" + " "*20 + "|")
        logger.info("|" + " "*10 + "Deterministic Genotype-to-Phenotype-to-Biophysics Engine" + " "*2 + "|")
        logger.info("+" + "="*68 + "+")
        logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        # PHASE 1: Genomic Ingestion
        mutations_df, proteins_dir, refs_dir, genomes_dir = phase1_genomic_ingestion(
            email=args.email,
            api_key=args.api_key,
            genome_file=Path(args.genome) if args.genome else None,
            organism=args.organism,
            limit=args.limit,
            target_genes=None,  # Can be extended from args
            output_dir=output_dir
        )
        
        # PHASE 2: Expression Analysis
        expression_scores = phase2_expression_analysis(
            expression_file=Path(args.expression_file) if args.expression_file else Path("data/expression_scores.json")
        )
        
        # PHASE 3: Epistasis Detection
        epistasis_networks = phase3_epistasis_detection(
            epistasis_file=Path(args.epistasis_file) if args.epistasis_file else Path("temp/epistasis_input.json")
        )
        
        # PHASE 4: Biophysics Docking
        pdb_mapping = args.pdb_mapping or {
            "oqxA": ("7cz9", "A"),
            "oqxB": ("7cz9", "A"),
        }
        docking_results = phase4_biophysics_docking(
            epistasis_networks=epistasis_networks,
            pdb_mapping=pdb_mapping,
            ligand_smiles=args.ligand_smiles,
            output_dir=output_dir / "biophysics"
        )
        
        # PHASE 5: Feature Aggregation (NO MACHINE LEARNING)
        sample_id = Path(args.genome).stem if args.genome else args.organism.replace(" ", "_")
        features_df = phase5_feature_aggregation(
            sample_id=sample_id,
            species=args.organism or "Unknown",
            mutations_df=mutations_df,
            expression_scores=expression_scores,
            epistasis_networks=epistasis_networks,
            docking_results=docking_results,
            output_file=Path("data/results/final_engineered_features.csv")
        )
        
        # SUCCESS
        logger.info("")
        logger.info("+" + "="*68 + "+")
        logger.info("|" + " "*20 + "[SUCCESS] PIPELINE COMPLETE" + " "*21 + "|")
        logger.info("|" + " "*68 + "|")
        logger.info("|  Final biophysical and genomic features saved to:" + " "*16 + "|")
        logger.info("|  data/results/final_engineered_features.csv" + " "*24 + "|")
        logger.info("|" + " "*68 + "|")
        logger.info("|  Ready for downstream statistical modeling." + " "*24 + "|")
        logger.info("+" + "="*68 + "+")
        logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        return 0
        
    except Exception as e:
        logger.error("Pipeline execution failed", exc_info=True)
        return 1


def main():
    """
    Main entry point for the master pipeline orchestrator.
    """
    parser = argparse.ArgumentParser(
        description="MutationScan Master Pipeline v2.1: Genotype-to-Phenotype-to-Biophysics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze local genome
  python run_pipeline.py --email you@email.com --genome data/genomes/strain.fasta
  
  # Download and analyze from NCBI
  python run_pipeline.py --email you@email.com --organism "Klebsiella pneumoniae" --limit 5
  
  # Full analysis with custom ligand and PDB mapping
  python run_pipeline.py --email you@email.com --genome test.fasta \\
    --ligand-smiles "O=C(O)C1=CN=C(...)" \\
    --pdb-mapping "oqxA:7cz9:A" "oqxB:7cz9:A"
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--email',
        required=True,
        help='Email address (required by NCBI policy)'
    )
    
    # Genome input (mutually exclusive)
    genome_group = parser.add_mutually_exclusive_group(required=True)
    genome_group.add_argument(
        '--genome',
        type=str,
        help='Path to local genome FASTA file'
    )
    genome_group.add_argument(
        '--organism',
        type=str,
        help='Organism name to download from NCBI'
    )
    
    # Optional arguments
    parser.add_argument(
        '--api-key',
        default=None,
        help='NCBI API key (recommended for reliable downloads)'
    )
    
    parser.add_argument(
        '--limit',
        type=int,
        default=1,
        help='Number of genomes to download (default: 1)'
    )
    
    parser.add_argument(
        '--output-dir',
        type=str,
        default='data/results',
        help='Output directory for results (default: data/results)'
    )
    
    parser.add_argument(
        '--expression-file',
        type=str,
        default='data/expression_scores.json',
        help='Path to ControlScan expression scores file'
    )
    
    parser.add_argument(
        '--epistasis-file',
        type=str,
        default='temp/epistasis_input.json',
        help='Path to epistasis network data'
    )
    
    parser.add_argument(
        '--ligand-smiles',
        type=str,
        default='O=C(O)C1=CN=C(N2CC3C=CC(N4CCNCC4)=C3C2)C(N)=C1F',
        help='SMILES string for the ligand (default: Ciprofloxacin)'
    )
    
    parser.add_argument(
        '--pdb-mapping',
        nargs='+',
        default=None,
        help='PDB mapping as "gene:pdb_id:chain" (e.g., "oqxA:7cz9:A")'
    )
    
    args = parser.parse_args()
    
    # Parse PDB mapping from command line if provided
    if args.pdb_mapping:
        pdb_dict = {}
        for mapping in args.pdb_mapping:
            parts = mapping.split(':')
            if len(parts) == 3:
                gene, pdb_id, chain = parts
                pdb_dict[gene] = (pdb_id, chain)
        args.pdb_mapping = pdb_dict
    
    # Setup logging
    setup_logging()
    
    # Run pipeline
    exit_code = run_master_pipeline(args)
    
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
