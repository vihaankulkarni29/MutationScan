"""
Test configuration and fixtures for SubScan Analyzer
"""

import pytest
import tempfile
import os
from pathlib import Path


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files"""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_alignment_data():
    """Sample alignment data for testing"""
    return {
        "reference_seq": ">reference\nMKRISIPSRRLRLPPPLERSHRLLPSPPPSLMSHHHHHHH",
        "mutated_seq": ">genome_1\nMKRISIPSRRLRLPPPLERSHRLLPSPPPSLMSHHAHHHH",
        "expected_mutations": [
            {
                "gene": "test_gene",
                "genome_id": "genome_1",
                "position": 34,
                "reference_aa": "H",
                "mutated_aa": "A",
                "mutation_type": "substitution",
            }
        ],
    }


@pytest.fixture
def sample_manifest():
    """Sample alignment manifest for testing"""
    return {
        "pipeline_step": "Aligner",
        "alignments": {
            "test_gene": {
                "alignment_file_path": "/path/to/test_gene_alignment.water",
                "gene_family": "test_gene",
            }
        },
    }
