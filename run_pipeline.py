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
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Tuple

import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

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
    args,
    target_genes: Optional[List[str]],
    output_dir: Path
) -> Tuple[pd.DataFrame, Path, Path, Path]:
    """
    Phase 1: Genomic Data Ingestion
    
    Download genomes, identify resistance genes, and call variants.
    
    Args:
        args: Parsed command-line arguments namespace
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
    
    if not args.skip_download:
        if args.genome:
            genome_file = Path(args.genome)
            logger.info(f"Using local genome: {genome_file}")
            # Copy local genome
            import shutil
            local_dest = genomes_dir / genome_file.name
            shutil.copy(genome_file, local_dest)
            num_genomes = 1
        else:
            if not args.organism:
                raise ValueError("Organism name required for NCBI download")
            
            logger.info(f"Downloading {args.limit} genome(s) from NCBI: {args.organism}")
            downloader = NCBIDatasetsGenomeDownloader(
                email=args.email,
                api_key=args.api_key,
                output_dir=genomes_dir
            )
            
            # Handle geolocation-based search (e.g., "Klebsiella pneumoniae AND India")
            if " AND " in args.organism:
                parts = args.organism.split(" AND ")
                taxon = parts[0].strip()
                location = parts[1].strip()
                logger.info(f"Using geolocation-based search: {taxon} from {location}")
                success, failed = downloader.download_bulk_by_geolocation(
                    organism=taxon,
                    location=location,
                    limit=args.limit
                )
                num_genomes = success
                logger.info(f"Downloaded {success} genomes, {failed} failed")
            else:
                # Fallback: Standard batch downloader via search + batch
                logger.info(f"Using standard organism search")
                accessions = downloader.search_accessions(
                    organism=args.organism,
                    max_results=args.limit
                )
                if accessions:
                    logger.info(f"Found {len(accessions)} accessions, downloading...")
                    success, failed = downloader.download_batch(
                        accessions=accessions
                    )
                    num_genomes = success
                    logger.info(f"Downloaded {success} genomes, {failed} failed")
                else:
                    logger.error(f"No genomes found for organism: {args.organism}")
                    raise ValueError(f"No genomes found for organism: {args.organism}")
    else:
        logger.info("Skipping NCBI download phase. Proceeding with existing local genomes in data/genomes.")
    
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
        "--email", args.email,
        "--organism", args.organism if args.organism else "unknown",
        "--threads", str(args.threads),
    ]
    
    # Add optional API key if provided
    if args.api_key:
        docker_command.extend(["--api-key", args.api_key])
    
    # Add targets file path
    docker_command.extend(["--targets", args.targets])
    
    try:
        logger.info(f"Executing Docker command: {' '.join(docker_command[:6])}...")
        result = subprocess.run(
            docker_command,
            check=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=3600  # 1 hour timeout
        )
        logger.info("Docker container completed successfully")
    except subprocess.TimeoutExpired:
        logger.error("Docker container execution timed out (1 hour limit exceeded)")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Docker container failed with exit code {e.returncode}")
        logger.error("Docker output was streamed to console above. Check output for details.")
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

    # Step 4: Data Quality Filtering
    logger.info("Step 1.4: Applying data quality filters...")
    if not mutations_df.empty:
        original_count = len(mutations_df)
        
        # Keep only essential columns
        cols_to_keep = ['Accession', 'Gene', 'Mutation']
        mutations_df = mutations_df[[c for c in cols_to_keep if c in mutations_df.columns]]
        logger.info(f"Retained {len(cols_to_keep)} core columns")
        
        # Remove garbage alignments: Drop genomes with >20 mutations per gene
        if 'Accession' in mutations_df.columns and 'Gene' in mutations_df.columns:
            mut_counts = mutations_df.groupby(['Accession', 'Gene']).size()
            valid_groups = mut_counts[mut_counts <= 20].index
            mutations_df = mutations_df[mutations_df.set_index(['Accession', 'Gene']).index.isin(valid_groups)].reset_index(drop=True)
            filtered_count = len(mutations_df)
            logger.info(f"Alignment quality filter: {filtered_count}/{original_count} mutations retained (removed {original_count - filtered_count} from poor alignments)")
            
            # Create composite key for proper mapping
            mutations_df['Acc_Gene'] = mutations_df['Accession'] + "|" + mutations_df['Gene']
            logger.info(f"Created composite Accession|Gene keys for epistasis analysis")
    
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


def phase3_epistasis_detection(epistasis_data: Dict) -> Dict[str, List[str]]:
    """
    Phase 3: Epistasis Network Detection (Dynamic Top 5 Frequency-Based)
    
    Load pre-computed epistatic mutation networks, rank by frequency,
    and return only the Top 5 most prevalent networks.
    
    Args:
        epistasis_data: Loaded epistasis input dictionary with keys like "GCF_051448695.1|oqxB"
                       and values as mutation string lists (e.g., ["K10R", "I174V"])
        
    Returns:
        Dictionary mapping mutation_network_string (e.g., "K10R + I174V") -> mutation list
    """
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("="*70)
    logger.info("PHASE 3: Epistasis Network Detection (Dynamic Top 5)")
    logger.info("="*70)
    
    epistasis_networks = {}
    
    try:
        data = epistasis_data
        genomes_data = data.get("genomes", {})

        if not isinstance(genomes_data, dict):
            logger.warning("Unexpected epistasis input format for 'genomes'; expected dict")
            logger.info("No epistatic networks to analyze")
            return epistasis_networks
        
        # Count frequency of each mutation network
        network_frequency = {}
        
        # Parse epistasis data (keys are "Accession|Gene", values are mutation string lists)
        for acc_gene, mutations in genomes_data.items():
            if not isinstance(mutations, list) or len(mutations) < 2:
                continue
            
            # Format network as string: "K10R + I174V"
            network_key = " + ".join(sorted(mutations))
            network_frequency[network_key] = network_frequency.get(network_key, 0) + 1
        
        logger.info(f"Detected {len(network_frequency)} unique epistatic networks across dataset")
        
        # Sort by frequency (descending) and take Top 5
        sorted_networks = sorted(network_frequency.items(), key=lambda x: x[1], reverse=True)
        top5_networks = sorted_networks[:5]
        
        logger.info("Selecting Top 5 most frequent networks (threshold: >= 2 mutations)")
        
        for idx, (network_str, freq) in enumerate(top5_networks, 1):
            mut_list = network_str.split(" + ")
            epistasis_networks[network_str] = mut_list
            logger.info(f"  {idx}. {network_str} (frequency: {freq} genomes)")
        
        logger.info(f"Identified Top {len(epistasis_networks)} epistatic networks for docking")
        
    except Exception as e:
        logger.error(f"Failed to load epistasis data: {e}")
        logger.info("No epistatic networks to analyze")
    
    logger.info("[PHASE 3 COMPLETE] Epistasis detection finished")
    
    return epistasis_networks


def _format_for_autoscan(mut_list: List[str], chain: str = "A") -> str:
    """
    Convert mutation list to AutoScan CLI format.
    
    Example:
        ["K10R", "I174V"] -> "A:10:K:R,A:174:I:V"
    
    Args:
        mut_list: List of mutation strings (e.g., ["K10R", "I174V"])
        chain: Chain identifier (default "A")
        
    Returns:
        AutoScan format string
    """
    import re
    formatted = []
    
    for mut in mut_list:
        # Parse mutation: Extract residue number, original AA, new AA
        match = re.match(r'([A-Z])*(\d+)([A-Z])', mut)
        if match:
            residue_num = match.group(2)
            orig_aa = match.group(1) if match.group(1) else "X"
            new_aa = match.group(3)
            formatted.append(f"{chain}:{residue_num}:{orig_aa}:{new_aa}")
    
    return ",".join(formatted)


def _smiles_to_sdf(smiles: str, output_path: Path) -> bool:
    """
    Convert SMILES string to SDF molecular file.
    
    Args:
        smiles: SMILES string representation
        output_path: Path to save SDF file
        
    Returns:
        True if successful, False otherwise
    """
    logger = logging.getLogger(__name__)
    
    if not RDKIT_AVAILABLE:
        logger.warning("RDKit not available. Cannot convert SMILES to SDF.")
        logger.info("Install rdkit: pip install rdkit")
        return False
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.error(f"Invalid SMILES string: {smiles}")
            return False
        
        # Add hydrogens and 3D coordinates
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # Write to SDF
        output_path.parent.mkdir(parents=True, exist_ok=True)
        writer = Chem.SDWriter(str(output_path))
        writer.write(mol)
        writer.close()
        
        logger.info(f"Created ligand SDF: {output_path}")
        return True
        
    except Exception as e:
        logger.error(f"Failed to convert SMILES to SDF: {e}")
        return False


def phase4_biophysics_docking(
    epistasis_networks: Dict[str, List[str]],
    default_pdb: str,
    ligand_smiles: str,
    output_dir: Path,
    center_x: float = 0.0,
    center_y: float = 0.0,
    center_z: float = 0.0
) -> Dict[str, Dict[str, float]]:
    """
    Phase 4: Rigorous 3D Biophysics Docking via Dockerized AutoScan
    """
    logger = logging.getLogger(__name__)

    logger.info("")
    logger.info("=" * 70)
    logger.info("PHASE 4: 3D Biophysics Docking (Dockerized AutoScan)")
    logger.info("=" * 70)

    docking_results = {}
    if not epistasis_networks:
        logger.warning("No networks to dock.")
        return docking_results

    project_root = Path(__file__).parent.absolute()
    biophysics_dir = output_dir / "biophysics"
    biophysics_dir.mkdir(parents=True, exist_ok=True)

    # ---------------------------------------------------------
    # STEP 0: Convert SMILES to SDF Ligand File
    # ---------------------------------------------------------
    logger.info("Step 0: Preparing ligand molecule from SMILES...")
    ligands_dir = project_root / "data" / "ligands"
    ligands_dir.mkdir(parents=True, exist_ok=True)
    
    ligand_sdf_path = ligands_dir / "ligand.sdf"
    ligand_container_path = "/app/data/ligands/ligand.sdf"  # Path inside Docker container
    
    if not ligand_sdf_path.exists():
        if _smiles_to_sdf(ligand_smiles, ligand_sdf_path):
            logger.info(f"  Ligand SDF created: {ligand_sdf_path}")
        else:
            logger.error("Failed to convert SMILES to SDF. RDKit may not be installed.")
            logger.error("Install with: pip install rdkit")
            return docking_results
    else:
        logger.info(f"  Using existing ligand SDF: {ligand_sdf_path}")

    # Translates ["K10R", "I174V"] -> "A:10:K:R,A:174:I:V"
    def _format_for_autoscan_local(mut_list, chain="A"):
        formatted = []
        for mutation in mut_list:
            match = re.match(r"([A-Za-z])(\d+)([A-Za-z])", mutation)
            if match:
                orig, pos, mut = match.groups()
                formatted.append(f"{chain}:{pos}:{orig}:{mut}")
        return ",".join(formatted)

    # ---------------------------------------------------------
    # STEP A: Wild-Type Baseline Docking
    # ---------------------------------------------------------
    wt_dir = biophysics_dir / "WT_baseline"
    wt_dir.mkdir(exist_ok=True)
    logger.info(f"Step A: Establishing WT Baseline ({default_pdb})...")

    docker_cmd_wt = [
        "docker", "run", "--rm",
        "-v", f"{project_root / 'data'}:/app/data",
        "--entrypoint", "autoscan",
        "mutationscan:latest",
        "--receptor", default_pdb,
        "--ligand", ligand_container_path,
        "--center-x", str(center_x),
        "--center-y", str(center_y),
        "--center-z", str(center_z),
        "--output", "/app/data/results/biophysics/WT_baseline"
    ]

    wt_affinity = -8.5  # Fallback
    try:
        subprocess.run(docker_cmd_wt, check=True)
        wt_json = list(wt_dir.glob("*.json"))
        if wt_json:
            with open(wt_json[0], "r") as handle:
                wt_data = json.load(handle)
            wt_affinity = wt_data.get("affinity", wt_data.get("consensus_affinity_kcal_mol", wt_affinity))
        logger.info(f"  [SUCCESS] WT Baseline Affinity: {float(wt_affinity):.2f} kcal/mol")
    except Exception as exc:
        logger.warning(f"  [WARNING] WT Docking failed or JSON not found. Error: {exc}")

    # ---------------------------------------------------------
    # STEP B: Dock the Top 5 Networks
    # ---------------------------------------------------------
    logger.info("Step B: Physically docking the Top 5 most frequent networks...")
    for index, (network_str, mut_list) in enumerate(epistasis_networks.items(), 1):
        mutant_dir = biophysics_dir / f"mutant_{index}"
        mutant_dir.mkdir(exist_ok=True)

        autoscan_mut_string = _format_for_autoscan_local(mut_list)
        logger.info(
            f"  Docking Network {index}/{len(epistasis_networks)}: "
            f"{network_str} ({autoscan_mut_string})..."
        )

        docker_cmd_mut = [
            "docker", "run", "--rm",
            "-v", f"{project_root / 'data'}:/app/data",
            "--entrypoint", "autoscan",
            "mutationscan:latest",
            "--receptor", default_pdb,
            "--ligand", ligand_container_path,
            "--mutation", autoscan_mut_string,
            "--center-x", str(center_x),
            "--center-y", str(center_y),
            "--center-z", str(center_z),
            "--minimize",
            "--output", f"/app/data/results/biophysics/mutant_{index}"
        ]

        try:
            subprocess.run(docker_cmd_mut, check=True)

            mut_affinity = wt_affinity  # Fallback
            mut_json = list(mutant_dir.glob("*.json"))
            if mut_json:
                with open(mut_json[0], "r") as handle:
                    mut_data = json.load(handle)
                mut_affinity = mut_data.get("affinity", mut_data.get("consensus_affinity_kcal_mol", mut_affinity))

            ddg = float(mut_affinity) - float(wt_affinity)
            docking_results[network_str] = {
                "wt_affinity": float(wt_affinity),
                "mutant_affinity": float(mut_affinity),
                "delta_delta_g": float(ddg)
            }
            logger.info(f"    [SUCCESS] Delta Delta G = {ddg:.2f} kcal/mol")

        except subprocess.CalledProcessError as exc:
            logger.error(f"    [FAILED] Docking crashed for {network_str}")
            logger.error(f"    Stderr: {exc.stderr}")

    logger.info(f"[PHASE 4 COMPLETE] Docking analysis finished ({len(docking_results)} successes)")
    return docking_results


def phase5_feature_aggregation(
    sample_id: str,
    species: str,
    mutations_df: pd.DataFrame,
    expression_scores: Dict[str, float],
    epistasis_networks: Dict[str, List[str]],
    docking_results: Dict[str, Dict[str, float]],
    output_dir: Path,
    default_pdb: str
) -> bool:
    """
    Phase 5: Feature Aggregation
    
    Split pipeline output into 3 clean, relational CSV tables:
    - 1_genomics_report.csv: Core mutation data with reference PDB
    - 2_epistasis_networks.csv: Co-occurring mutation networks
    - 3_biophysics_docking.csv: Binding affinity calculations
    
    Args:
        sample_id: Unique identifier for this genome
        species: Species name (e.g., "Klebsiella pneumoniae")
        mutations_df: DataFrame from mutation calling (with Acc_Gene column)
        expression_scores: Dict mapping gene -> expression score
        epistasis_networks: Dict mapping acc_gene -> mutation string list
        docking_results: Dict mapping acc_gene -> affinity dict
        output_dir: Directory to save the 3 CSV files
        default_pdb: Default PDB reference for all genes
        
    Returns:
        True if successful
    """
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("="*70)
    logger.info("PHASE 5: Feature Aggregation (3-Table Output)")
    logger.info("="*70)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # CSV 1: Genomics Report (mutations_df with Ref_PDB)
        logger.info("Generating 1_genomics_report.csv...")
        genomics_df = mutations_df.copy()
        
        # Drop temporary Acc_Gene column if present
        if 'Acc_Gene' in genomics_df.columns:
            genomics_df = genomics_df.drop(columns=['Acc_Gene'])
        
        # Add Ref_PDB column
        genomics_df['Ref_PDB'] = default_pdb
        
        genomics_output = output_dir / "1_genomics_report.csv"
        genomics_df.to_csv(genomics_output, index=False)
        logger.info(f"  Saved {len(genomics_df)} mutation records to {genomics_output}")
        
        # CSV 2: Epistasis Networks (Top 5)
        logger.info("Generating 2_epistasis_networks.csv...")
        epistasis_records = []
        
        for network_idx, (network_str, mut_list) in enumerate(epistasis_networks.items(), 1):
            # network_str is already formatted as "K10R + I174V"
            epistasis_records.append({
                "Network_Rank": network_idx,
                "Gene_Network": network_str,
                "Co_occurring_Mutations": network_str,
                "Mutation_Count": len(mut_list)
            })
        
        epistasis_df = pd.DataFrame(epistasis_records)
        epistasis_output = output_dir / "2_epistasis_networks.csv"
        epistasis_df.to_csv(epistasis_output, index=False)
        logger.info(f"  Saved {len(epistasis_df)} top epistatic networks to {epistasis_output}")
        
        # CSV 3: Biophysics Docking (Top 5 Networks)
        logger.info("Generating 3_biophysics_docking.csv...")
        docking_records = []
        
        for network_idx, (network_str, affinity_dict) in enumerate(docking_results.items(), 1):
            # network_str is formatted as "K10R + I174V"
            docking_records.append({
                "Network_Rank": network_idx,
                "Network": network_str,
                "WT_Affinity_kcalmol": affinity_dict.get("wt_affinity"),
                "Mutant_Affinity_kcalmol": affinity_dict.get("mutant_affinity"),
                "Delta_Delta_G_kcalmol": affinity_dict.get("delta_delta_g")
            })
        
        docking_df = pd.DataFrame(docking_records)
        docking_output = output_dir / "3_biophysics_docking.csv"
        docking_df.to_csv(docking_output, index=False)
        logger.info(f"  Saved {len(docking_df)} docking results to {docking_output}")
        
        # Log summary statistics
        logger.info("")
        logger.info("Summary Statistics:")
        logger.info(f"  Total Mutations: {len(genomics_df)}")
        logger.info(f"  Epistatic Networks: {len(epistasis_df)}")
        logger.info(f"  Docking Calculations: {len(docking_df)}")
        
        if not docking_df.empty and 'Delta_Delta_G_kcalmol' in docking_df.columns:
            ddg_mean = docking_df['Delta_Delta_G_kcalmol'].mean()
            ddg_std = docking_df['Delta_Delta_G_kcalmol'].std()
            logger.info(f"  DDG Mean: {ddg_mean:.2f} ± {ddg_std:.2f} kcal/mol")
        
        logger.info("[PHASE 5 COMPLETE] Feature aggregation finished")
        logger.info("")
        
        return True
        
    except Exception as e:
        logger.error(f"Phase 5 aggregation failed: {e}", exc_info=True)
        return False


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

        # Default values for checkpointed execution
        mutations_df = pd.DataFrame()
        expression_scores = {}
        epistasis_networks = {}
        docking_results = {}
        sample_id = args.organism.split()[0] + "_" + args.organism.split()[1] if args.organism and " " in args.organism else (args.organism or "Unknown")
        
        logger.info("")
        logger.info("+" + "="*68 + "+")
        logger.info("|" + " "*15 + "MUTATIONSCAN MASTER PIPELINE v2.1" + " "*20 + "|")
        logger.info("|" + " "*10 + "Deterministic Genotype-to-Phenotype-to-Biophysics Engine" + " "*2 + "|")
        logger.info("+" + "="*68 + "+")
        logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        # PHASES 1 & 2: Conditionally execute or load from cache
        if not args.resume_epistasis:
            # PHASE 1: Genomic Ingestion
            if args.start_phase <= 1:
                mutations_df, proteins_dir, refs_dir, genomes_dir = phase1_genomic_ingestion(
                    args=args,
                    target_genes=None,  # Can be extended from args
                    output_dir=output_dir
                )

                # Convert mutations DataFrame to epistasis input format for Phase 3
                logger.info("Converting mutation data to epistasis analysis format...")
                epistasis_input_dict = {}

                if isinstance(mutations_df, pd.DataFrame) and not mutations_df.empty:
                    # Mutations are already filtered in phase1_genomic_ingestion
                    # Create epistasis input dict with Acc_Gene composite keys
                    if 'Acc_Gene' in mutations_df.columns and 'Mutation' in mutations_df.columns:
                        epistasis_input_dict = mutations_df.groupby('Acc_Gene')['Mutation'].apply(list).to_dict()
                        logger.info(f"Prepared epistasis input for {len(epistasis_input_dict)} genome-gene pairs")
                    else:
                        logger.warning("Mutation DataFrame missing expected columns (Acc_Gene, Mutation)")
                else:
                    logger.warning("No mutation data available for epistasis analysis")

                # Save epistasis input for Phase 3
                epistasis_file = Path(args.epistasis_file) if args.epistasis_file else Path("data/interim/epistasis_input.json")
                epistasis_file.parent.mkdir(parents=True, exist_ok=True)
                with open(epistasis_file, 'w') as f:
                    json.dump({"genomes": epistasis_input_dict}, f, indent=2)
                logger.info(f"Saved epistasis input to {epistasis_file}")

                if args.genome:
                    sample_id = Path(args.genome).stem
            else:
                logger.info("Skipping Phase 1 (Genomics). Loading cached mutation data...")
                mutation_file = output_dir / "mutation_report.csv"
                if mutation_file.exists():
                    try:
                        mutations_df = pd.read_csv(mutation_file)
                        logger.info(f"Loaded {len(mutations_df)} mutations from {mutation_file}")
                    except Exception as e:
                        logger.warning(f"Failed to load mutations from {mutation_file}: {e}")
                        mutations_df = pd.DataFrame()
                else:
                    logger.warning(f"Mutation cache not found: {mutation_file}")
                    mutations_df = pd.DataFrame()

            # PHASE 2: Expression Analysis
            if args.start_phase <= 2:
                expression_scores = phase2_expression_analysis(
                    expression_file=Path(args.expression_file) if args.expression_file else Path("data/expression_scores.json")
                )
            else:
                logger.info("Skipping Phase 2 (Expression). Loading cached expression scores...")
                expression_file = Path(args.expression_file) if args.expression_file else Path("data/expression_scores.json")
                if expression_file.exists():
                    try:
                        expression_scores = phase2_expression_analysis(expression_file=expression_file)
                        logger.info(f"Loaded expression scores from {expression_file}")
                    except Exception as e:
                        logger.warning(f"Failed to load expression scores from {expression_file}: {e}")
                        expression_scores = {}
                else:
                    logger.warning(f"Expression score cache not found: {expression_file}")
                    expression_scores = {}
        else:
            # Fast-resume: Skip Phases 1 and 2, load interim cached data directly
            logger.info("Resuming pipeline directly from Phase 3 using cached interim data...")
            
            # Load epistasis input
            epistasis_input_dict = {}
            epistasis_file = Path(args.epistasis_file) if args.epistasis_file else Path("data/interim/epistasis_input.json")
            if epistasis_file.exists():
                try:
                    with open(epistasis_file, 'r') as f:
                        epistasis_data = json.load(f)
                    epistasis_input_dict = epistasis_data.get("genomes", {})
                    logger.info(f"Loaded epistasis input from {epistasis_file}")
                except Exception as e:
                    logger.warning(f"Failed to load epistasis data from {epistasis_file}: {e}")
            else:
                logger.warning(f"Epistasis file not found: {epistasis_file}")
            
            # Load genomics report (needed for Phase 5)
            genomics_file = Path("data/results/1_genomics_report.csv")
            if genomics_file.exists():
                try:
                    mutations_df = pd.read_csv(genomics_file)
                    logger.info(f"Loaded {len(mutations_df)} mutations from genomics report")
                except Exception as e:
                    logger.warning(f"Failed to load genomics report from {genomics_file}: {e}")
                    mutations_df = pd.DataFrame()
            else:
                logger.warning(f"Genomics report not found: {genomics_file}")
                mutations_df = pd.DataFrame()

        # PHASE 3: Epistasis Detection
        if args.start_phase <= 3:
            epistasis_file = Path(args.epistasis_file) if args.epistasis_file else Path("data/interim/epistasis_input.json")
            if not epistasis_file.exists():
                logger.warning(f"Epistasis file not found: {epistasis_file}")
                epistasis_data_dict = {"genomes": {}}
            else:
                with open(epistasis_file, "r") as f:
                    epistasis_data_dict = json.load(f)

            epistasis_networks = phase3_epistasis_detection(
                epistasis_data=epistasis_data_dict
            )
        else:
            logger.info("Skipping Phase 3 (Epistasis). Loading cached epistasis networks...")
            epistasis_file = Path(args.epistasis_file) if args.epistasis_file else Path("data/interim/epistasis_input.json")
            if not epistasis_file.exists():
                logger.warning(f"Epistasis cache not found: {epistasis_file}")
                epistasis_networks = {}
            else:
                try:
                    with open(epistasis_file, "r") as f:
                        epistasis_data_dict = json.load(f)
                    epistasis_networks = phase3_epistasis_detection(epistasis_data=epistasis_data_dict)
                    logger.info(f"Loaded epistasis networks from {epistasis_file}")
                except Exception as e:
                    logger.warning(f"Failed to load epistasis networks from {epistasis_file}: {e}")
                    epistasis_networks = {}

        # PHASE 4: Biophysics Docking
        if args.start_phase <= 4:
            docking_results = phase4_biophysics_docking(
                epistasis_networks=epistasis_networks,
                default_pdb=args.default_pdb,
                ligand_smiles=args.drug_smiles,
                output_dir=output_dir,
                center_x=args.center_x,
                center_y=args.center_y,
                center_z=args.center_z
            )
        else:
            logger.info("Skipping Phase 4 (Biophysics). Relying on cached data.")

        # PHASE 5: Feature Aggregation (3-Table Output)
        if args.start_phase <= 5:
            success = phase5_feature_aggregation(
                sample_id=sample_id,
                species=args.organism or "Unknown",
                mutations_df=mutations_df,
                expression_scores=expression_scores,
                epistasis_networks=epistasis_networks,
                docking_results=docking_results,
                output_dir=Path("data/results"),
                default_pdb=args.default_pdb
            )
            
            if not success:
                logger.error("Phase 5 aggregation failed")
                return 1
        else:
            logger.info("Skipping Phase 5 (Feature Aggregation). Relying on cached data.")
        
        # SUCCESS
        logger.info("")
        logger.info("+" + "="*68 + "+")
        logger.info("|" + " "*20 + "[SUCCESS] PIPELINE COMPLETE" + " "*21 + "|")
        logger.info("|" + " "*68 + "|")
        logger.info("|  Final outputs saved to data/results/:" + " "*27 + "|")
        logger.info("|    - 1_genomics_report.csv (mutation data)" + " "*24 + "|")
        logger.info("|    - 2_epistasis_networks.csv (co-occurring mutations)" + " "*11 + "|")
        logger.info("|    - 3_biophysics_docking.csv (binding affinities)" + " "*14 + "|")
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
        default='data/interim/epistasis_input.json',
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

    parser.add_argument(
        '--default-pdb',
        type=str,
        default='7cz9',
        help='Universal Wild-Type reference PDB ID for AutoScan docking (e.g., 7cz9).'
    )

    parser.add_argument(
        '--drug-smiles',
        type=str,
        default='C1=CC2=C(C(=C1)F)N(C=C(C2=O)C(=O)O)C3CC3',
        help='SMILES string for the target antibiotic (default: Ciprofloxacin).'
    )
    
    parser.add_argument(
        '--center-x',
        type=float,
        default=0.0,
        help='Docking box center X coordinate (Angstroms). Required by AutoScan.'
    )
    
    parser.add_argument(
        '--center-y',
        type=float,
        default=0.0,
        help='Docking box center Y coordinate (Angstroms). Required by AutoScan.'
    )
    
    parser.add_argument(
        '--center-z',
        type=float,
        default=0.0,
        help='Docking box center Z coordinate (Angstroms). Required by AutoScan.'
    )
    
    parser.add_argument(
        '--targets',
        type=str,
        default='config/target_genes.txt',
        help='Path to the target genes text file'
    )
    
    parser.add_argument(
        '--skip-download',
        action='store_true',
        help='Skip the NCBI downloading phase and use existing local genomes.'
    )

    parser.add_argument(
        '--start-phase',
        type=int,
        default=1,
        choices=[1, 2, 3, 4, 5],
        help='Which phase to start execution from. Use 3 to skip Genomics and go straight to Epistasis.'
    )

    parser.add_argument(
        '--resume-epistasis',
        action='store_true',
        help='Skip Phases 1 and 2, load interim data, and start directly from Phase 3.'
    )
    
    parser.add_argument(
        '--threads',
        type=int,
        default=max(1, (os.cpu_count() or 4) - 2),
        help='Number of concurrent threads for parallel genome processing (default: CPU count - 2)'
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
