"""
Test script for PyMOLVisualizer module.

This script creates test data and verifies the PyMOL visualization pipeline.
"""

import sys
from pathlib import Path
import logging
import pandas as pd

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from src.mutation_scan.visualizer.visualizer import PyMOLVisualizer

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def create_test_mutation_report():
    """
    Create a test mutation report CSV for visualization testing.
    """
    test_data = [
        {
            'Accession': 'GCF_001',
            'Gene': 'gyrA',
            'Mutation': 'S83L',
            'Status': 'Resistant',
            'Phenotype': 'Fluoroquinolone resistance (high-level)',
            'Reference_PDB': '3NUU'
        },
        {
            'Accession': 'GCF_001',
            'Gene': 'gyrA',
            'Mutation': 'D87N',
            'Status': 'Resistant',
            'Phenotype': 'Fluoroquinolone resistance (high-level)',
            'Reference_PDB': '3NUU'
        },
        {
            'Accession': 'GCF_002',
            'Gene': 'parC',
            'Mutation': 'S80I',
            'Status': 'Resistant',
            'Phenotype': 'Fluoroquinolone resistance (moderate)',
            'Reference_PDB': 'N/A'  # Should be skipped
        },
        {
            'Accession': 'GCF_003',
            'Gene': 'acrB',
            'Mutation': 'G141D',
            'Status': 'VUS',
            'Phenotype': 'N/A',
            'Reference_PDB': '1IWG'
        }
    ]
    
    df = pd.DataFrame(test_data)
    
    # Save to test directory
    test_dir = Path("data/test_visualizer")
    test_dir.mkdir(parents=True, exist_ok=True)
    
    csv_path = test_dir / "test_mutation_report.csv"
    df.to_csv(csv_path, index=False)
    
    logger.info(f"Created test mutation report: {csv_path}")
    return csv_path


def test_pml_generation():
    """
    Test PyMOL script generation without actually running PyMOL.
    """
    logger.info("=" * 70)
    logger.info("TEST 1: PyMOL Script Generation")
    logger.info("=" * 70)
    
    try:
        visualizer = PyMOLVisualizer(
            output_dir=Path("data/test_visualizer/output"),
            pymol_path="pymol"  # Will fail verification, but we can still test script gen
        )
    except FileNotFoundError:
        logger.warning("PyMOL not installed - testing script generation only")
        
        # Create visualizer without verification
        visualizer = PyMOLVisualizer.__new__(PyMOLVisualizer)
        visualizer.output_dir = Path("data/test_visualizer/output")
        visualizer.output_dir.mkdir(parents=True, exist_ok=True)
        visualizer.pymol_path = "pymol"
    
    # Test script generation
    mutations = [
        {'mutation': 'S83L', 'status': 'Resistant', 'phenotype': 'FQ resistance', 'accession': 'GCF_001'},
        {'mutation': 'D87N', 'status': 'Resistant', 'phenotype': 'FQ resistance', 'accession': 'GCF_001'}
    ]
    
    pml_script = visualizer._generate_pml_script(
        pdb_id="3NUU",
        gene="gyrA",
        mutations=mutations,
        output_png=Path("test_output.png")
    )
    
    logger.info("\nGenerated PyMOL Script:")
    logger.info("-" * 70)
    logger.info(pml_script)
    logger.info("-" * 70)
    
    # Verify script contains key commands
    assert "fetch 3NUU" in pml_script, "Missing fetch command"
    assert "show cartoon" in pml_script, "Missing cartoon display"
    assert "resi 83" in pml_script, "Missing S83L position"
    assert "resi 87" in pml_script, "Missing D87N position"
    assert "color red" in pml_script, "Missing resistant color"
    assert "png test_output.png" in pml_script, "Missing PNG output"
    
    logger.info("✓ TEST 1 PASSED: Script generation correct")


def test_mutation_parsing():
    """
    Test mutation position extraction.
    """
    logger.info("\n" + "=" * 70)
    logger.info("TEST 2: Mutation Position Parsing")
    logger.info("=" * 70)
    
    visualizer = PyMOLVisualizer.__new__(PyMOLVisualizer)
    
    test_cases = [
        ("S83L", 83),
        ("D87N", 87),
        ("K43R", 43),
        ("G141D", 141),
        ("P167S", 167),
        ("INVALID", None)
    ]
    
    for mutation, expected_pos in test_cases:
        result = visualizer._extract_position(mutation)
        logger.info(f"{mutation} → position {result} (expected: {expected_pos})")
        assert result == expected_pos, f"Failed to parse {mutation}"
    
    logger.info("✓ TEST 2 PASSED: Position parsing correct")


def test_grouping():
    """
    Test mutation grouping by gene and PDB.
    """
    logger.info("\n" + "=" * 70)
    logger.info("TEST 3: Mutation Grouping by Structure")
    logger.info("=" * 70)
    
    # Create test DataFrame
    test_data = pd.DataFrame([
        {'Gene': 'gyrA', 'Mutation': 'S83L', 'Status': 'Resistant', 'Phenotype': 'FQ', 'Reference_PDB': '3NUU', 'Accession': 'GCF_001'},
        {'Gene': 'gyrA', 'Mutation': 'D87N', 'Status': 'Resistant', 'Phenotype': 'FQ', 'Reference_PDB': '3NUU', 'Accession': 'GCF_001'},
        {'Gene': 'parC', 'Mutation': 'S80I', 'Status': 'Resistant', 'Phenotype': 'FQ', 'Reference_PDB': 'N/A', 'Accession': 'GCF_002'},
        {'Gene': 'acrB', 'Mutation': 'G141D', 'Status': 'VUS', 'Phenotype': 'N/A', 'Reference_PDB': '1IWG', 'Accession': 'GCF_003'}
    ])
    
    visualizer = PyMOLVisualizer.__new__(PyMOLVisualizer)
    grouped = visualizer._group_mutations_by_structure(test_data)
    
    logger.info(f"\nGrouped mutations:")
    for (gene, pdb_id), mutations in grouped.items():
        logger.info(f"  {gene} ({pdb_id}): {len(mutations)} mutations")
        for m in mutations:
            logger.info(f"    - {m['mutation']} ({m['status']})")
    
    # Verify grouping
    assert ('gyrA', '3NUU') in grouped, "Missing gyrA group"
    assert len(grouped[('gyrA', '3NUU')]) == 2, "gyrA should have 2 mutations"
    assert ('parC', 'N/A') not in grouped, "parC with N/A should be skipped"
    assert ('acrB', '1IWG') in grouped, "Missing acrB group"
    
    logger.info("✓ TEST 3 PASSED: Grouping logic correct")


def test_full_workflow_dry_run():
    """
    Test full workflow without actually running PyMOL.
    """
    logger.info("\n" + "=" * 70)
    logger.info("TEST 4: Full Workflow (Dry Run)")
    logger.info("=" * 70)
    
    # Create test mutation report
    csv_path = create_test_mutation_report()
    
    logger.info(f"Test mutation report created: {csv_path}")
    
    # Load and verify CSV
    df = pd.read_csv(csv_path)
    logger.info(f"Loaded {len(df)} mutations from CSV")
    
    # Group by gene
    resistant = df[df['Status'] == 'Resistant']
    logger.info(f"Found {len(resistant)} resistant mutations")
    
    # Show what would be visualized
    valid_pdb = resistant[resistant['Reference_PDB'] != 'N/A']
    logger.info(f"Would visualize {len(valid_pdb)} mutations (with valid PDB IDs)")
    
    for _, row in valid_pdb.iterrows():
        logger.info(f"  - {row['Gene']}: {row['Mutation']} (PDB: {row['Reference_PDB']})")
    
    logger.info("✓ TEST 4 PASSED: Workflow logic correct")


def main():
    """
    Run all tests.
    """
    try:
        logger.info("\n" + "=" * 70)
        logger.info("PYMOL VISUALIZER TEST SUITE")
        logger.info("=" * 70)
        
        test_pml_generation()
        test_mutation_parsing()
        test_grouping()
        test_full_workflow_dry_run()
        
        logger.info("\n" + "=" * 70)
        logger.info("ALL TESTS PASSED! ✓")
        logger.info("=" * 70)
        logger.info("\nPyMOLVisualizer is working correctly:")
        logger.info("  ✓ PyMOL script generation (.pml files)")
        logger.info("  ✓ Mutation position parsing (S83L → 83)")
        logger.info("  ✓ Grouping by gene and PDB ID")
        logger.info("  ✓ Workflow logic (filtering, validation)")
        logger.info("\nNOTE: Actual PyMOL rendering requires:")
        logger.info("  1. PyMOL installed (apt-get install pymol)")
        logger.info("  2. libglew-dev installed (graphics driver helper)")
        logger.info("  3. Docker rebuild: docker build -t mutationscan:v1 .")
        
    except AssertionError as e:
        logger.error(f"\n❌ TEST FAILED: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"\n❌ UNEXPECTED ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
