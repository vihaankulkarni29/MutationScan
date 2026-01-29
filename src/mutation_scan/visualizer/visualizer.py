"""
PyMOL Visualizer module for MutationScan.

This module automates PyMOL to generate 3D protein structure visualizations
with mutations highlighted in color.

CRITICAL ARCHITECTURE NOTES:
- This module writes code that writes code (Python → PyMOL .pml script → PNG)
- Uses subprocess to call PyMOL CLI (not pymol Python library)
- Headless rendering with -c flag (no GUI required)
- Groups multiple mutations on same structure

ANTI-HALLUCINATION RULES:
- Never use the pymol Python library (not standard installation)
- Always use subprocess to call: pymol -c -q script.pml
- Skip visualization if PDB_ID is "N/A"
- Validate PDB IDs are 4 characters long
- Log errors but do not crash pipeline
"""

import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd

logger = logging.getLogger(__name__)


class PyMOLVisualizer:
    """
    Automate PyMOL to visualize mutations on 3D protein structures.
    
    This class takes mutation_report.csv from VariantCaller and generates
    high-quality PNG images showing mutations highlighted on protein structures.
    
    Key Features:
    - Headless PyMOL rendering (no GUI required)
    - Groups mutations by gene (multiple mutations on same structure)
    - Color-coded by status: Resistant (red), VUS (orange)
    - High-quality ray-traced rendering
    - Error handling (invalid PDB IDs logged but don't crash)
    - .pml script generation for reproducibility
    """

    def __init__(self, output_dir: Path, pymol_path: str = "pymol"):
        """
        Initialize PyMOLVisualizer.

        Args:
            output_dir: Directory to save PNG images and .pml scripts
            pymol_path: Path to PyMOL executable (default: "pymol" assumes in PATH)

        Raises:
            FileNotFoundError: If PyMOL is not installed or not in PATH
        """
        self.output_dir = Path(output_dir)
        self.pymol_path = pymol_path
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Verify PyMOL installation
        self._verify_pymol_installation()
        
        logger.info(f"Initialized PyMOLVisualizer with output_dir: {self.output_dir}")

    def _verify_pymol_installation(self) -> None:
        """
        Verify PyMOL is installed and accessible.

        Raises:
            FileNotFoundError: If PyMOL is not found
        """
        try:
            result = subprocess.run(
                [self.pymol_path, "-c", "-e", "quit"],
                capture_output=True,
                timeout=10
            )
            logger.info("PyMOL installation verified")
        except FileNotFoundError:
            raise FileNotFoundError(
                f"PyMOL not found at '{self.pymol_path}'. "
                "Please install PyMOL or provide correct path."
            )
        except subprocess.TimeoutExpired:
            logger.warning("PyMOL verification timed out (may still work)")

    def visualize_mutations(
        self,
        mutation_csv: Path,
        filter_status: List[str] = ["Resistant", "VUS"]
    ) -> Dict[str, List[Path]]:
        """
        Generate PyMOL visualizations for all mutations in CSV.

        Args:
            mutation_csv: Path to mutation_report.csv from VariantCaller
            filter_status: List of statuses to visualize (default: ["Resistant", "VUS"])

        Returns:
            Dictionary mapping gene → list of generated PNG paths

        Example:
            {
                "gyrA": [Path("gyrA_S83L_D87N.png")],
                "parC": [Path("parC_S80I.png")]
            }
        """
        # Load mutation report
        if not mutation_csv.exists():
            raise FileNotFoundError(f"Mutation CSV not found: {mutation_csv}")
        
        mutations_df = pd.read_csv(mutation_csv)
        
        # Filter by status
        filtered = mutations_df[mutations_df['Status'].isin(filter_status)]
        
        if filtered.empty:
            logger.warning(f"No mutations with status {filter_status} found")
            return {}
        
        logger.info(f"Processing {len(filtered)} mutations for visualization")
        
        # Group by gene and PDB ID
        grouped = self._group_mutations_by_structure(filtered)
        
        # Generate visualizations
        results = {}
        
        for (gene, pdb_id), mutations_list in grouped.items():
            logger.info(f"Visualizing {gene} (PDB: {pdb_id}) with {len(mutations_list)} mutations")
            
            try:
                png_path = self._visualize_structure(
                    gene=gene,
                    pdb_id=pdb_id,
                    mutations=mutations_list
                )
                
                if png_path:
                    if gene not in results:
                        results[gene] = []
                    results[gene].append(png_path)
                    
            except Exception as e:
                logger.error(f"Failed to visualize {gene} ({pdb_id}): {e}")
                continue
        
        logger.info(f"Generated {sum(len(v) for v in results.values())} visualizations")
        return results

    def _group_mutations_by_structure(
        self,
        mutations_df: pd.DataFrame
    ) -> Dict[Tuple[str, str], List[Dict]]:
        """
        Group mutations by gene and PDB ID.

        Args:
            mutations_df: Filtered mutations DataFrame

        Returns:
            Dictionary mapping (gene, pdb_id) → list of mutation dicts
        """
        grouped = {}
        
        for _, row in mutations_df.iterrows():
            gene = row['Gene']
            pdb_id = row['Reference_PDB']
            
            # Skip if no PDB ID available
            if pd.isna(pdb_id) or pdb_id == "N/A":
                logger.debug(f"Skipping {gene} {row['Mutation']} (no PDB ID)")
                continue
            
            # Validate PDB ID format (4 characters)
            if len(str(pdb_id)) != 4:
                logger.warning(f"Invalid PDB ID '{pdb_id}' for {gene} (must be 4 characters)")
                continue
            
            key = (gene, pdb_id)
            
            if key not in grouped:
                grouped[key] = []
            
            grouped[key].append({
                'mutation': row['Mutation'],
                'status': row['Status'],
                'phenotype': row['Phenotype'],
                'accession': row['Accession']
            })
        
        return grouped

    def _visualize_structure(
        self,
        gene: str,
        pdb_id: str,
        mutations: List[Dict]
    ) -> Optional[Path]:
        """
        Generate PyMOL visualization for a single structure with mutations.

        Args:
            gene: Gene name (e.g., "gyrA")
            pdb_id: PDB ID (e.g., "3NUU")
            mutations: List of mutation dictionaries

        Returns:
            Path to generated PNG file, or None if failed
        """
        # Generate output filename
        mutation_str = "_".join([m['mutation'] for m in mutations])
        output_name = f"{gene}_{mutation_str}_{pdb_id}"
        pml_file = self.output_dir / f"{output_name}.pml"
        png_file = self.output_dir / f"{output_name}.png"
        
        # Generate PyMOL script
        pml_script = self._generate_pml_script(
            pdb_id=pdb_id,
            gene=gene,
            mutations=mutations,
            output_png=png_file
        )
        
        # Save .pml script
        with open(pml_file, 'w') as f:
            f.write(pml_script)
        
        logger.debug(f"Generated PyMOL script: {pml_file}")
        
        # Execute PyMOL
        success = self._run_pymol(pml_file)
        
        if success and png_file.exists():
            logger.info(f"Generated visualization: {png_file}")
            return png_file
        else:
            logger.error(f"PyMOL failed to generate {png_file}")
            return None

    def _generate_pml_script(
        self,
        pdb_id: str,
        gene: str,
        mutations: List[Dict],
        output_png: Path
    ) -> str:
        """
        Generate PyMOL .pml script content.

        This is the "code that writes code" method.

        Args:
            pdb_id: PDB structure ID
            gene: Gene name
            mutations: List of mutation dicts
            output_png: Output PNG path

        Returns:
            PyMOL script as string
        """
        script_lines = [
            "# PyMOL Script Generated by MutationScan Visualizer",
            f"# Gene: {gene}",
            f"# PDB: {pdb_id}",
            f"# Mutations: {', '.join([m['mutation'] for m in mutations])}",
            "",
            "# 1. FETCH STRUCTURE (SYNCHRONOUS)",
            f"fetch {pdb_id}, async=0",
            "",
            "# 2. CLEAN VIEW",
            "hide all",
            "show cartoon",
            "color white, all",
            "",
            "# 3. HIGHLIGHT MUTATIONS",
        ]
        
        # Add selection and coloring for each mutation
        for i, mutation_dict in enumerate(mutations):
            mutation = mutation_dict['mutation']
            status = mutation_dict['status']
            
            # Parse mutation (e.g., "S83L" → position 83)
            position = self._extract_position(mutation)
            
            if position is None:
                logger.warning(f"Could not parse position from {mutation}")
                continue
            
            # Color based on status
            color = "red" if status == "Resistant" else "orange"
            
            selection_name = f"mut_{i+1}"
            mutation_label = mutation
            
            script_lines.extend([
                f"",
                f"# Mutation {i+1}: {mutation} ({status})",
                f"select {selection_name}, resi {position}",
                f"show spheres, {selection_name}",
                f"color {color}, {selection_name}",
                f"label {selection_name} and name CA, \"{mutation_label}\"",
            ])
        
        # Calculate zoom target (average of all mutation positions if multiple)
        if mutations:
            positions = [self._extract_position(m['mutation']) for m in mutations]
            positions = [p for p in positions if p is not None]
            
            if positions:
                # If we have mutations, zoom to show them prominently
                script_lines.extend([
                    "",
                    "# 4. CAMERA SETUP",
                    f"zoom resi {'+'.join(map(str, positions))}, 10",
                ])
            else:
                script_lines.extend([
                    "",
                    "# 4. CAMERA SETUP",
                    "zoom",
                ])
        else:
            script_lines.extend([
                "",
                "# 4. CAMERA SETUP",
                "zoom",
            ])
        
        # Add rendering commands
        script_lines.extend([
            "",
            "# 5. RENDER HIGH-QUALITY IMAGE",
            "set ray_shadows, 0",
            "set antialias, 2",
            f"png {str(output_png)}, width=1200, height=1200, ray=1",
            "",
            "# 6. QUIT",
            "quit"
        ])
        
        return "\n".join(script_lines)

    def _extract_position(self, mutation: str) -> Optional[int]:
        """
        Extract position number from mutation string.

        Args:
            mutation: Mutation string (e.g., "S83L", "D87N")

        Returns:
            Position as integer, or None if parsing fails

        Examples:
            "S83L" → 83
            "D87N" → 87
            "K43R" → 43
        """
        try:
            # Remove first character (reference AA) and last character (mutant AA)
            # Middle part should be the position
            position_str = mutation[1:-1]
            return int(position_str)
        except (ValueError, IndexError):
            return None

    def _run_pymol(self, pml_file: Path) -> bool:
        """
        Execute PyMOL script in headless mode.

        Args:
            pml_file: Path to .pml script

        Returns:
            True if successful, False otherwise
        """
        try:
            # Run PyMOL in headless mode
            # -c: command line only (no GUI)
            # -q: quiet mode (less verbose)
            result = subprocess.run(
                [self.pymol_path, "-c", "-q", str(pml_file)],
                capture_output=True,
                text=True,
                timeout=60  # 60 second timeout
            )
            
            if result.returncode == 0:
                logger.debug(f"PyMOL executed successfully: {pml_file}")
                return True
            else:
                logger.error(f"PyMOL failed with return code {result.returncode}")
                logger.error(f"STDERR: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            logger.error(f"PyMOL execution timed out for {pml_file}")
            return False
        except Exception as e:
            logger.error(f"PyMOL execution error: {e}")
            return False

    def get_summary(self, results: Dict[str, List[Path]]) -> Dict:
        """
        Generate summary statistics from visualization results.

        Args:
            results: Output from visualize_mutations()

        Returns:
            Dictionary with summary statistics
        """
        total_genes = len(results)
        total_images = sum(len(paths) for paths in results.values())
        
        summary = {
            'total_genes_visualized': total_genes,
            'total_images_generated': total_images,
            'genes': list(results.keys()),
            'files': {gene: [str(p) for p in paths] for gene, paths in results.items()}
        }
        
        return summary


# =============================================================================
# USAGE EXAMPLE
# =============================================================================
if __name__ == "__main__":
    """
    Example usage of PyMOLVisualizer.
    
    Typical workflow:
    1. Run GenomeExtractor → Download genomes
    2. Run GeneFinder → Find resistance genes
    3. Run SequenceExtractor → Extract and translate sequences
    4. Run VariantCaller → Identify mutations
    5. Run PyMOLVisualizer → Generate 3D visualizations
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Initialize visualizer
    visualizer = PyMOLVisualizer(
        output_dir=Path("data/results/visualizations")
    )
    
    # Visualize mutations from VariantCaller output
    results = visualizer.visualize_mutations(
        mutation_csv=Path("data/results/mutation_report.csv"),
        filter_status=["Resistant", "VUS"]
    )
    
    # Print summary
    summary = visualizer.get_summary(results)
    print("\n=== VISUALIZATION SUMMARY ===")
    print(f"Genes visualized: {summary['total_genes_visualized']}")
    print(f"Images generated: {summary['total_images_generated']}")
    
    for gene, files in summary['files'].items():
        print(f"\n{gene}:")
        for f in files:
            print(f"  - {f}")
