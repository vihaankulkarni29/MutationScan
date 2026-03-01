"""
AutoScan Bridge Module - Population Genomics to 3D Biophysics

This module connects MutationScan's population-level mutation analysis to the
AutoScan 3D biophysics engine. It converts identified epistatic mutation networks
into in silico mutagenesis configurations and calculates binding affinity changes
(ΔΔG) for specific ligands (e.g., antibiotics).

Workflow:
1. Accept PDB ID, chain, and residue list from MutationScan
2. Generate AutoScan-compatible configuration (YAML)
3. Execute AutoScan for WT and mutant docking
4. Parse results and return ΔΔG values

Author: Senior Bioinformatics Software Architect
"""

import json
import logging
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Optional, Dict, List

import yaml

logger = logging.getLogger(__name__)


class AutoScanBridge:
    """
    Bridge between MutationScan (population genomics) and AutoScan (3D biophysics).
    
    Converts epistatic mutation networks into binding affinity predictions by
    running comparative docking simulations (WT vs mutant).
    
    Attributes:
        autoscan_path: Path or command to AutoScan executable (e.g., "python -m autoscan.main")
    """
    
    def __init__(self, autoscan_path: str = "python -m autoscan.main"):
        """
        Initialize the AutoScan bridge.
        
        Args:
            autoscan_path: Path or command to AutoScan executable.
                          Can be a full path or module invocation.
                          Default: "python -m autoscan.main"
                          
        Raises:
            ValueError: If autoscan_path is empty
        """
        if not autoscan_path or not isinstance(autoscan_path, str):
            raise ValueError("autoscan_path must be a non-empty string")
        
        self.autoscan_path = autoscan_path
        logger.info(f"Initialized AutoScanBridge with: {self.autoscan_path}")
    
    def _generate_mutation_config(
        self,
        pdb_id: str,
        chain: str,
        residues: List[int],
        output_dir: Path
    ) -> Path:
        """
        Generate YAML configuration file for AutoScan mutagenesis simulation.
        
        The configuration specifies the WT PDB, chain, and residues to be
        mutated for in silico docking analysis.
        
        Args:
            pdb_id: PDB ID (e.g., "7cz9")
            chain: Chain identifier (e.g., "A")
            residues: List of residue numbers to mutate (e.g., [10, 174, 970])
            output_dir: Directory to save the configuration file
            
        Returns:
            Path to the generated YAML configuration file
            
        Raises:
            ValueError: If inputs are invalid
        """
        # Validate inputs
        if not pdb_id or not isinstance(pdb_id, str):
            raise ValueError("pdb_id must be a non-empty string")
        if not chain or not isinstance(chain, str) or len(chain) != 1:
            raise ValueError("chain must be a single character")
        if not residues or not isinstance(residues, list):
            raise ValueError("residues must be a non-empty list")
        if not all(isinstance(r, int) and r > 0 for r in residues):
            raise ValueError("All residues must be positive integers")
        
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Build configuration dictionary
        config = {
            "metadata": {
                "source": "MutationScan Population Analysis",
                "pdb_id": pdb_id.lower(),
                "chain": chain.upper(),
                "mutation_count": len(residues),
                "generated_at": time.strftime("%Y-%m-%d %H:%M:%S")
            },
            "mutagenesis": {
                "pdb_id": pdb_id.lower(),
                "chain": chain.upper(),
                "residues": sorted(residues),  # Sort for reproducibility
                "strategy": "single_to_wild_type"  # Mutate residues back to identify WT interaction impact
            },
            "docking": {
                "engine": "AutoDock Vina",
                "minimize_md": True,
                "exhaustiveness": 8,
                "num_modes": 9
            },
            "output": {
                "save_poses": True,
                "save_trajectories": False  # Set to True if MD trajectories needed
            }
        }
        
        # Generate timestamped filename for reproducibility
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        config_file = output_dir / f"autoscan_config_{pdb_id}_{chain}_{timestamp}.yaml"
        
        try:
            with open(config_file, 'w') as f:
                yaml.dump(config, f, default_flow_style=False, sort_keys=False)
            
            logger.info(f"Generated AutoScan config: {config_file}")
            logger.debug(f"  PDB: {pdb_id}, Chain: {chain}, Residues: {len(residues)}")
            
            return config_file
            
        except Exception as e:
            logger.error(f"Failed to generate config file {config_file}: {e}")
            raise
    
    def run_comparative_docking(
        self,
        pdb_id: str,
        chain: str,
        residues: List[int],
        ligand_smiles: str,
        output_dir: Optional[Path] = None
    ) -> Optional[Dict[str, float]]:
        """
        Execute comparative docking analysis (WT vs mutant).
        
        Runs AutoScan to calculate binding affinity changes (ΔΔG) by comparing
        WT docking to in silico mutant docking.
        
        Args:
            pdb_id: PDB ID for the target protein (e.g., "7cz9")
            chain: Chain identifier (e.g., "A")
            residues: List of residue numbers to mutate (e.g., [10, 174, 970])
            ligand_smiles: SMILES string for the ligand (e.g., ciprofloxacin SMILES)
            output_dir: Optional directory for AutoScan output. 
                       If None, uses temporary directory.
            
        Returns:
            Dictionary with keys:
            - "wt_affinity": WT binding affinity (kcal/mol, negative = better binding)
            - "mutant_affinity": Mutant binding affinity (kcal/mol)
            - "delta_delta_g": ΔΔG value (mutant - WT; negative = improved binding)
            
            Returns None if docking fails (logged as error).
            
        Raises:
            ValueError: If input parameters are invalid
        """
        # Validate inputs
        if not ligand_smiles or not isinstance(ligand_smiles, str):
            raise ValueError("ligand_smiles must be a non-empty string")
        
        # Create output directory
        if output_dir is None:
            output_dir = Path(tempfile.mkdtemp(prefix="autoscan_"))
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        logger.info("")
        logger.info("="*70)
        logger.info("AUTOSCAN BRIDGE: Comparative Docking Analysis")
        logger.info("="*70)
        logger.info(f"PDB: {pdb_id}, Chain: {chain}, Residues: {residues}")
        logger.info(f"Ligand SMILES: {ligand_smiles}")
        logger.info(f"Output directory: {output_dir}")
        
        try:
            # Step A: Generate mutation configuration
            logger.info("Step A: Generating mutation configuration...")
            config_file = self._generate_mutation_config(
                pdb_id=pdb_id,
                chain=chain,
                residues=residues,
                output_dir=output_dir
            )
            
            # Step B: Execute AutoScan
            logger.info("Step B: Executing AutoScan pipeline...")
            self._execute_autoscan(
                config_file=config_file,
                ligand_smiles=ligand_smiles,
                output_dir=output_dir
            )
            
            # Step C: Parse results
            logger.info("Step C: Parsing AutoScan results...")
            results = self._parse_autoscan_results(output_dir)
            
            # Step D: Return results
            if results:
                logger.info("Step D: Results compiled successfully")
                logger.info(f"  WT Affinity: {results['wt_affinity']:.2f} kcal/mol")
                logger.info(f"  Mutant Affinity: {results['mutant_affinity']:.2f} kcal/mol")
                logger.info(f"  ΔΔG: {results['delta_delta_g']:.2f} kcal/mol")
                
                return results
            else:
                logger.error("Failed to parse AutoScan results")
                return None
                
        except subprocess.CalledProcessError as e:
            logger.error(f"AutoScan execution failed with exit code {e.returncode}")
            logger.error(f"Command: {e.cmd}")
            logger.error(f"stdout: {e.stdout}")
            logger.error(f"stderr: {e.stderr}")
            return None
            
        except Exception as e:
            logger.error(f"Comparative docking failed: {e}", exc_info=True)
            return None
    
    def _execute_autoscan(
        self,
        config_file: Path,
        ligand_smiles: str,
        output_dir: Path,
        timeout: int = 3600
    ) -> None:
        """
        Execute AutoScan subprocess with configuration and ligand.
        
        Args:
            config_file: Path to YAML configuration file
            ligand_smiles: SMILES string for the ligand
            output_dir: Output directory for results
            timeout: Timeout in seconds (default: 1 hour)
            
        Raises:
            subprocess.CalledProcessError: If AutoScan execution fails
        """
        # Build command
        cmd = self.autoscan_path.split() + [
            "--config", str(config_file),
            "--ligand", ligand_smiles,
            "--output-dir", str(output_dir),
            "--minimize-md"  # Enable MD minimization for more accurate predictions
        ]
        
        logger.debug(f"Executing: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=True  # Raise CalledProcessError if exit code != 0
            )
            
            logger.debug(f"AutoScan stdout: {result.stdout}")
            if result.stderr:
                logger.warning(f"AutoScan stderr: {result.stderr}")
                
        except subprocess.TimeoutExpired:
            logger.error(f"AutoScan timed out after {timeout}s")
            raise
        except subprocess.CalledProcessError as e:
            logger.error(f"AutoScan failed with exit code {e.returncode}")
            raise
    
    def _parse_autoscan_results(self, output_dir: Path) -> Optional[Dict[str, float]]:
        """
        Parse AutoScan output to extract binding affinity values.
        
        Looks for either:
        1. benchmark_result.json (preferred)
        2. benchmark_result.txt (fallback)
        
        Expected JSON format:
        {
            "wt_affinity": -8.5,
            "mutant_affinity": -7.2,
            "delta_delta_g": 1.3
        }
        
        Args:
            output_dir: Directory containing AutoScan output
            
        Returns:
            Dictionary with affinity values, or None if parsing fails
        """
        output_dir = Path(output_dir)
        
        # Try JSON first (preferred)
        json_file = output_dir / "benchmark_result.json"
        if json_file.exists():
            try:
                with open(json_file, 'r') as f:
                    results = json.load(f)
                
                logger.debug(f"Parsed results from {json_file}")
                
                # Validate required keys
                required_keys = {"wt_affinity", "mutant_affinity", "delta_delta_g"}
                if not required_keys.issubset(results.keys()):
                    logger.error(f"Missing required keys in {json_file}: {required_keys - set(results.keys())}")
                    return None
                
                return {
                    "wt_affinity": float(results["wt_affinity"]),
                    "mutant_affinity": float(results["mutant_affinity"]),
                    "delta_delta_g": float(results["delta_delta_g"])
                }
                
            except json.JSONDecodeError as e:
                logger.error(f"Failed to parse JSON from {json_file}: {e}")
                return None
            except (ValueError, KeyError) as e:
                logger.error(f"Invalid data in {json_file}: {e}")
                return None
        
        # Try text file (fallback)
        txt_file = output_dir / "benchmark_result.txt"
        if txt_file.exists():
            try:
                with open(txt_file, 'r') as f:
                    lines = f.readlines()
                
                results = {}
                for line in lines:
                    line = line.strip()
                    if "wt_affinity" in line.lower():
                        # Parse line like: "WT Affinity: -8.5 kcal/mol"
                        value = float(line.split(':')[-1].split()[0])
                        results["wt_affinity"] = value
                    elif "mutant_affinity" in line.lower():
                        value = float(line.split(':')[-1].split()[0])
                        results["mutant_affinity"] = value
                    elif "delta_delta_g" in line.lower() or "ddg" in line.lower():
                        value = float(line.split(':')[-1].split()[0])
                        results["delta_delta_g"] = value
                
                if len(results) == 3:
                    logger.debug(f"Parsed results from {txt_file}")
                    return results
                else:
                    logger.error(f"Could not extract all required values from {txt_file}")
                    return None
                    
            except (ValueError, IndexError) as e:
                logger.error(f"Failed to parse text from {txt_file}: {e}")
                return None
        
        logger.error(f"No result files found in {output_dir}")
        logger.debug(f"Available files: {list(output_dir.glob('*'))}")
        return None
