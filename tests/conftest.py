"""Test fixtures and utilities."""

import pytest
from pathlib import Path

@pytest.fixture
def sample_dna_sequence():
    """Provide sample DNA sequence for testing."""
    return "ATGAAAGCGTAACGTAAA"

@pytest.fixture
def sample_protein_sequence():
    """Provide sample protein sequence for testing."""
    return "MKAYA"

@pytest.fixture
def sample_config():
    """Provide sample configuration dictionary."""
    return {
        "ncbi": {"email": "test@example.com"},
        "abricate": {"database": "card", "min_identity": 90},
        "blast": {"evalue": 1e-5}
    }
