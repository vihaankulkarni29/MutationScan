#!/usr/bin/env python3
"""
run_top5_docking.py

Rigorous physical docking performance on the top 5 epistatic mutation networks
identified in the clinical AMR dataset using AutoScan (molecular docking).

This script performs baseline (WT) docking, then mutant docking for each of the 5 
networks, calculates Delta_Delta_G, and generates a detailed CSV report.

Usage:
    python run_top5_docking.py [--no-minimize] [--consensus]

Author: Senior Structural Biologist & Bioinformatics Architect
Date: March 2, 2026
"""

import argparse
import json
import logging
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# HARDCODED TOP 5 EPISTATIC NETWORKS (Clinical Dataset)
# ============================================================================

TOP5_NETWORKS = [
    {
        "network_id": 1,
        "name": "Core Triad",
        "description": "Minimal high-impact network",
        "mutation_string": "A:10:K:R,A:174:I:V,A:970:L:A",
        "mutations": ["A:10:K:R", "A:174:I:V", "A:970:L:A"]
    },
    {
        "network_id": 2,
        "name": "Triad + S540G",
        "description": "Core triad with promoter mutation",
        "mutation_string": "A:10:K:R,A:174:I:V,A:540:S:G,A:970:L:A",
        "mutations": ["A:10:K:R", "A:174:I:V", "A:540:S:G", "A:970:L:A"]
    },
    {
        "network_id": 3,
        "name": "Triad + Y782F + P1050L",
        "description": "Core triad with distant contact mutations",
        "mutation_string": "A:10:K:R,A:174:I:V,A:782:Y:F,A:970:L:A,A:1050:P:L",
        "mutations": ["A:10:K:R", "A:174:I:V", "A:782:Y:F", "A:970:L:A", "A:1050:P:L"]
    },
    {
        "network_id": 4,
        "name": "Triad + R342S + S540G + I961M",
        "description": "Hyper-resistant pentamer with stabilizing mutations",
        "mutation_string": "A:10:K:R,A:174:I:V,A:342:R:S,A:540:S:G,A:961:I:M,A:970:L:A",
        "mutations": ["A:10:K:R", "A:174:I:V", "A:342:R:S", "A:540:S:G", "A:961:I:M", "A:970:L:A"]
    },
    {
        "network_id": 5,
        "name": "Hyper-Mutator Anchor (D100N)",
        "description": "Single anchor mutation with pleiotropic effects",
        "mutation_string": "A:100:D:N",
        "mutations": ["A:100:D:N"]
    }
]


# ============================================================================
# BASELINE PARAMETERS (OqxB / PDB 7cz9)
# ============================================================================

BASELINE_CONFIG = {
    "receptor_pdb": "7cz9",              # PDB ID or local path
    "ligand_smiles": "C1=CC2=C(C(=C1)F)N(C=C(C2=O)C(=O)O)C3CC3",  # Ciprofloxacin
    "ligand_name": "ciprofloxacin",
    "chain": "A",                         # Active site chain
    "binding_pocket": "oqxB_pocket",      # Pocket name key
    
    # Grid box parameters (OqxB binding pocket coordinates)
    # These are estimated based on known QNRAV binding site
    "grid_center": {
        "x": 15.2,
        "y": 22.8,
        "z": 18.5
    }
}


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def setup_directories() -> Path:
    """Create necessary output directories."""
    output_dir = Path("data/results")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create subdirectories for docking results
    (output_dir / "top5_docking").mkdir(exist_ok=True)
    
    return output_dir


def download_pdb(pdb_id: str, output_path: Path) -> bool:
    """
    Download PDB file from RCSB using curl.
    
    Args:
        pdb_id: PDB identifier (e.g., "7cz9")
        output_path: Where to save the PDB file
        
    Returns:
        True if successful, False otherwise
    """
    pdb_id_lower = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id_lower}.pdb"
    
    try:
        logger.info(f"Downloading {pdb_id} from RCSB PDB...")
        result = subprocess.run(
            ["curl", "-s", "-o", str(output_path), url],
            timeout=30,
            capture_output=True,
            text=True
        )
        
        if result.returncode == 0 and output_path.exists():
            logger.info(f"  ✓ Downloaded {pdb_id} to {output_path}")
            return True
        else:
            logger.error(f"  ✗ Failed to download {pdb_id}")
            return False
            
    except FileNotFoundError:
        logger.error("curl not found. Install curl or provide PDB file manually.")
        return False
    except Exception as e:
        logger.error(f"Download failed: {e}")
        return False


def create_ligand_pdbqt(smiles: str, output_path: Path) -> bool:
    """
    Create PDBQT from SMILES using OpenBabel.
    
    Args:
        smiles: SMILES string for ligand
        output_path: Where to save PDBQT
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Create temporary PDB from SMILES
        temp_pdb = output_path.with_suffix(".pdb")
        
        logger.info(f"Converting SMILES to ligand structure...")
        # Use obabel to convert SMILES to PDB
        result = subprocess.run(
            [
                "obabel", "-ismi", "-O", str(temp_pdb),
                "-xh"  # Add hydrogens
            ],
            input=smiles,
            text=True,
            capture_output=True,
            timeout=30
        )
        
        if result.returncode != 0:
            logger.error(f"OpenBabel SMILES conversion failed: {result.stderr}")
            return False
        
        logger.info(f"  ✓ Created ligand PDB: {temp_pdb}")
        
        # Convert PDB to PDBQT using obabel
        result = subprocess.run(
            [
                "obabel", str(temp_pdb), "-O", str(output_path),
                "-xh", "-xn", "-xpdbqt"
            ],
            capture_output=True,
            timeout=30
        )
        
        if result.returncode == 0 and output_path.exists():
            logger.info(f"  ✓ Converted to PDBQT: {output_path}")
            return True
        else:
            logger.error(f"OpenBabel PDB->PDBQT conversion failed: {result.stderr}")
            return False
            
    except FileNotFoundError:
        logger.error("OpenBabel not found. Install via: conda install -c conda-forge openbabel")
        return False
    except Exception as e:
        logger.error(f"Ligand preparation failed: {e}")
        return False


def run_autoscan_dock(
    receptor: Path,
    ligand: Path,
    center_x: float,
    center_y: float,
    center_z: float,
    mutation: Optional[str] = None,
    output_json: Optional[Path] = None,
    minimize: bool = False,
    consensus: bool = False
) -> Optional[Dict]:
    """
    Execute AutoScan docking via CLI subprocess.
    
    Args:
        receptor: Receptor file path (PDB or PDBQT)
        ligand: Ligand file path (PDB or PDBQT)
        center_x, center_y, center_z: Grid box center coordinates
        mutation: Mutation string (e.g., "A:87:D:G") or None for WT
        output_json: Path to save JSON output
        minimize: Whether to run energy minimization
        consensus: Whether to run consensus scoring
        
    Returns:
        Parsed JSON result dict or None if failed
    """
    
    # Build command
    cmd = [
        "autoscan", "dock",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center-x", str(center_x),
        "--center-y", str(center_y),
        "--center-z", str(center_z)
    ]
    
    # Add optional flags
    if mutation:
        cmd.extend(["--mutation", mutation])
    
    if minimize:
        cmd.append("--minimize")
    
    if consensus:
        cmd.append("--use-consensus")
        cmd.extend(["--consensus-method", "weighted"])
    
    if output_json:
        cmd.extend(["--output", str(output_json)])
    
    # Run docking
    try:
        logger.info(f"Running AutoScan: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300  # 5 min timeout
        )
        
        if result.returncode != 0:
            logger.error(f"AutoScan failed: {result.stderr}")
            return None
        
        # Parse output
        if output_json and output_json.exists():
            with open(output_json, 'r') as f:
                output_data = json.load(f)
            return output_data
        else:
            # Try parsing stdout
            for line in result.stdout.split('\n'):
                if 'binding_affinity_kcal_mol' in line.lower() or 'affinity' in line.lower():
                    logger.info(f"  {line.strip()}")
            
            # If output_json not specified, try to extract from stdout
            try:
                # AutoScan should return JSON on stdout
                output_data = json.loads(result.stdout)
                return output_data
            except json.JSONDecodeError:
                logger.error(f"Could not parse AutoScan output: {result.stdout}")
                return None
        
    except FileNotFoundError:
        logger.error("AutoScan CLI not found. Install via: pip install autoscan")
        return None
    except subprocess.TimeoutExpired:
        logger.error("AutoScan docking timed out (>5 min)")
        return None
    except Exception as e:
        logger.error(f"Docking execution failed: {e}")
        return None


def extract_binding_affinity(docking_result: Dict) -> Optional[float]:
    """
    Extract binding affinity from AutoScan JSON output.
    
    Args:
        docking_result: Parsed JSON from AutoScan
        
    Returns:
        Binding affinity in kcal/mol or None if not found
    """
    if not docking_result:
        return None
    
    # Check consensus mode first
    if docking_result.get("consensus_mode") and "consensus_affinity_kcal_mol" in docking_result:
        return float(docking_result["consensus_affinity_kcal_mol"])
    
    # Standard mode
    if "binding_affinity_kcal_mol" in docking_result:
        return float(docking_result["binding_affinity_kcal_mol"])
    
    # Fallback checks
    if "affinity" in docking_result:
        return float(docking_result["affinity"])
    
    logger.warning(f"Could not extract binding affinity from: {docking_result}")
    return None


# ============================================================================
# MAIN DOCKING WORKFLOW
# ============================================================================

def run_top5_docking(args):
    """
    Main workflow: baseline docking + 5 mutant dockings + analysis.
    
    Args:
        args: Command-line arguments
    """
    
    logger.info("="*75)
    logger.info("TOP 5 EPISTATIC NETWORKS: RIGOROUS PHYSICAL DOCKING (AutoScan)")
    logger.info("="*75)
    logger.info(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("")
    
    # Setup
    output_dir = setup_directories()
    docking_subdir = output_dir / "top5_docking"
    data_dir = Path("data/proteins")  # or local PDB storage
    
    # Prepare receptor (PDB)
    receptor_pdb = data_dir / f"{BASELINE_CONFIG['receptor_pdb']}.pdb"
    if not receptor_pdb.exists():
        logger.info(f"PDB file not found: {receptor_pdb}")
        logger.info("Attempting to download from RCSB...")
        data_dir.mkdir(parents=True, exist_ok=True)
        if not download_pdb(BASELINE_CONFIG['receptor_pdb'], receptor_pdb):
            logger.error("Failed to obtain PDB file. Exiting.")
            return 1
    else:
        logger.info(f"Using local PDB: {receptor_pdb}")
    
    # Prepare ligand (PDBQT)
    ligand_pdbqt = docking_subdir / f"{BASELINE_CONFIG['ligand_name']}.pdbqt"
    if not ligand_pdbqt.exists():
        logger.info(f"Preparing ligand from SMILES...")
        if not create_ligand_pdbqt(BASELINE_CONFIG['ligand_smiles'], ligand_pdbqt):
            logger.error("Failed to prepare ligand. Exiting.")
            return 1
    else:
        logger.info(f"Using cached ligand: {ligand_pdbqt}")
    
    # Grid box coordinates
    grid_center = BASELINE_CONFIG['grid_center']
    logger.info(f"Grid box center: ({grid_center['x']}, {grid_center['y']}, {grid_center['z']})")
    
    logger.info("")
    logger.info("-"*75)
    logger.info("STEP A: BASELINE DOCKING (Wild-Type)")
    logger.info("-"*75)
    
    # Baseline docking
    wt_output_json = docking_subdir / "baseline_wt.json"
    wt_result = run_autoscan_dock(
        receptor=receptor_pdb,
        ligand=ligand_pdbqt,
        center_x=grid_center['x'],
        center_y=grid_center['y'],
        center_z=grid_center['z'],
        mutation=None,
        output_json=wt_output_json,
        minimize=args.minimize,
        consensus=args.consensus
    )
    
    wt_affinity = extract_binding_affinity(wt_result)
    if wt_affinity is None:
        logger.error("Failed to obtain WT binding affinity. Exiting.")
        return 1
    
    logger.info(f"✓ WT Binding Affinity: {wt_affinity:.2f} kcal/mol")
    logger.info("")
    
    # Storage for results
    docking_results = []
    
    # Step B: Mutant docking
    logger.info("-"*75)
    logger.info("STEP B: MUTANT DOCKING (Top 5 Networks)")
    logger.info("-"*75)
    
    for network in TOP5_NETWORKS:
        network_id = network["network_id"]
        network_name = network["name"]
        mut_string = network["mutation_string"]
        num_mutations = len(network["mutations"])
        
        logger.info(f"\nNetwork {network_id}: {network_name} ({num_mutations} mutations)")
        logger.info(f"  Mutations: {', '.join(network['mutations'])}")
        
        # Run docking
        mutant_output_json = docking_subdir / f"network_{network_id}_mutant.json"
        mutant_result = run_autoscan_dock(
            receptor=receptor_pdb,
            ligand=ligand_pdbqt,
            center_x=grid_center['x'],
            center_y=grid_center['y'],
            center_z=grid_center['z'],
            mutation=mut_string,
            output_json=mutant_output_json,
            minimize=args.minimize,
            consensus=args.consensus
        )
        
        mutant_affinity = extract_binding_affinity(mutant_result)
        if mutant_affinity is None:
            logger.warning(f"  ⚠ Failed to dock network {network_id}")
            continue
        
        # Step C: Calculate Delta_Delta_G
        delta_delta_g = mutant_affinity - wt_affinity
        
        logger.info(f"  WT Affinity: {wt_affinity:.2f} kcal/mol")
        logger.info(f"  Mutant Affinity: {mutant_affinity:.2f} kcal/mol")
        logger.info(f"  ΔΔG: {delta_delta_g:+.2f} kcal/mol", end="")
        
        # Interpret ΔΔG
        if delta_delta_g > 0.5:
            logger.info(" [DESTABILIZING] ✗")
        elif delta_delta_g < -0.5:
            logger.info(" [STABILIZING] ✓")
        else:
            logger.info(" [NEUTRAL]")
        
        # Store result
        docking_results.append({
            "Network_ID": network_id,
            "Network_Name": network_name,
            "Description": network["description"],
            "Mutation_Count": num_mutations,
            "Mutation_String": mut_string,
            "WT_Affinity_kcalmol": round(wt_affinity, 3),
            "Mutant_Affinity_kcalmol": round(mutant_affinity, 3),
            "Delta_Delta_G_kcalmol": round(delta_delta_g, 3),
            "Docking_Status": "SUCCESS"
        })
    
    # Step D: Output generation
    logger.info("")
    logger.info("-"*75)
    logger.info("STEP D: RESULTS EXPORT")
    logger.info("-"*75)
    
    # Create DataFrame
    df = pd.DataFrame(docking_results)
    
    # Sort by |ΔΔG| (most significant changes first)
    df['abs_ddg'] = df['Delta_Delta_G_kcalmol'].abs()
    df_sorted = df.sort_values('abs_ddg', ascending=False).drop(columns=['abs_ddg'])
    
    # Save to CSV
    output_csv = output_dir / "top5_physical_docking.csv"
    df_sorted.to_csv(output_csv, index=False)
    logger.info(f"✓ Saved results to: {output_csv}")
    
    # Print formatted table
    logger.info("")
    logger.info("="*75)
    logger.info("SUMMARY TABLE: TOP 5 EPISTATIC NETWORKS - BINDING AFFINITIES")
    logger.info("="*75)
    
    # Format table columns
    print("\n")
    print("┌─────────────────────────────────────────────────────────────────────────────┐")
    print("│ Network │ Name                         │ Mutations │ WT_Aff │ Mut_Aff │ ΔΔG │")
    print("├─────────────────────────────────────────────────────────────────────────────┤")
    
    for _, row in df_sorted.iterrows():
        net_id = int(row['Network_ID'])
        name = str(row['Network_Name'])[:25].ljust(25)
        n_mut = int(row['Mutation_Count'])
        wt_aff = float(row['WT_Affinity_kcalmol'])
        mut_aff = float(row['Mutant_Affinity_kcalmol'])
        ddg = float(row['Delta_Delta_G_kcalmol'])
        
        # Color coding
        if ddg > 0.5:
            ddg_symbol = "↑"  # Destabilizing
        elif ddg < -0.5:
            ddg_symbol = "↓"  # Stabilizing
        else:
            ddg_symbol = "–"  # Neutral
        
        print(f"│ {net_id:^7} │ {name} │ {n_mut:^9} │ {wt_aff:^6.2f} │ {mut_aff:^7.2f} │ {ddg_symbol}{ddg:^5.2f} │")
    
    print("└─────────────────────────────────────────────────────────────────────────────┘\n")
    
    # Statistics
    logger.info("STATISTICAL SUMMARY:")
    ddg_values = df_sorted['Delta_Delta_G_kcalmol'].values
    logger.info(f"  Mean ΔΔG: {ddg_values.mean():+.2f} kcal/mol")
    logger.info(f"  Std Dev: {ddg_values.std():.2f} kcal/mol")
    logger.info(f"  Min: {ddg_values.min():+.2f} kcal/mol")
    logger.info(f"  Max: {ddg_values.max():+.2f} kcal/mol")
    
    # Interpretation
    stabilizing = (ddg_values < -0.5).sum()
    destabilizing = (ddg_values > 0.5).sum()
    neutral = len(ddg_values) - stabilizing - destabilizing
    
    logger.info(f"\nCLINICAL INTERPRETATION:")
    logger.info(f"  Stabilizing networks (ΔΔG < -0.5): {stabilizing}/{len(ddg_values)}")
    logger.info(f"  Neutral networks: {neutral}/{len(ddg_values)}")
    logger.info(f"  Destabilizing networks (ΔΔG > 0.5): {destabilizing}/{len(ddg_values)}")
    
    logger.info("")
    logger.info("="*75)
    logger.info("[SUCCESS] Physical docking analysis complete")
    logger.info("="*75)
    logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    return 0


# ============================================================================
# CLI ENTRY POINT
# ============================================================================

def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Top 5 Epistatic Networks: Rigorous Physical Docking (AutoScan)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard docking
  python run_top5_docking.py
  
  # With energy minimization (slower but more accurate)
  python run_top5_docking.py --minimize
  
  # With consensus scoring (requires multiple docking engines)
  python run_top5_docking.py --consensus
  
  # Full rigorous mode
  python run_top5_docking.py --minimize --consensus
        """
    )
    
    parser.add_argument(
        "--minimize",
        action="store_true",
        help="Enable OpenMM energy minimization for mutant structures (slower, more accurate)"
    )
    
    parser.add_argument(
        "--consensus",
        action="store_true",
        help="Enable consensus scoring from multiple docking engines (requires Vina, AutoDock, etc.)"
    )
    
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="Skip PDB download if not found locally"
    )
    
    args = parser.parse_args()
    
    return run_top5_docking(args)


if __name__ == "__main__":
    sys.exit(main())
