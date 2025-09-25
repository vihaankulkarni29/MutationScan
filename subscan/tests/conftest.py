"""
Pytest configuration and fixtures for MutationScan testing.

This module provides shared fixtures and configuration for testing the
MutationScan pipeline components.
"""

import os
import tempfile
import pytest
import json
from pathlib import Path
from typing import Dict, Any, List


@pytest.fixture(scope="session")
def test_data_dir():
    """Get the test data directory path."""
    return Path(__file__).parent / "data"


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def sample_genome_manifest():
    """Create a sample genome manifest for testing."""
    return {
        "pipeline_step": "genome_harvesting",
        "timestamp": "2025-09-24T10:00:00",
        "genomes": {
            "NZ_CP107554": {
                "accession": "NZ_CP107554",
                "organism": "Escherichia coli",
                "assembly_level": "Complete Genome",
                "genome_size": 5231428,
                "gc_content": 50.8,
                "quality_score": 95.2,
                "fasta_path": "genomes/NZ_CP107554.fasta",
                "download_date": "2025-09-24"
            },
            "NZ_CP107555": {
                "accession": "NZ_CP107555",
                "organism": "Klebsiella pneumoniae",
                "assembly_level": "Complete Genome",
                "genome_size": 5644321,
                "gc_content": 57.3,
                "quality_score": 92.1,
                "fasta_path": "genomes/NZ_CP107555.fasta",
                "download_date": "2025-09-24"
            }
        },
        "summary": {
            "total_genomes": 2,
            "successful_downloads": 2,
            "failed_downloads": 0,
            "average_quality_score": 93.65
        }
    }


@pytest.fixture
def sample_annotation_manifest():
    """Create a sample annotation manifest for testing."""
    return {
        "pipeline_step": "amr_annotation",
        "timestamp": "2025-09-24T10:30:00",
        "input_manifest": "genome_manifest.json",
        "genomes": {
            "NZ_CP107554": {
                "accession": "NZ_CP107554",
                "fasta_path": "genomes/NZ_CP107554.fasta",
                "abricate_results": "annotations/NZ_CP107554_abricate.tab",
                "amr_genes_found": 5,
                "gene_families": ["blaOXA", "aac", "qnr"]
            },
            "NZ_CP107555": {
                "accession": "NZ_CP107555",
                "fasta_path": "genomes/NZ_CP107555.fasta",
                "abricate_results": "annotations/NZ_CP107555_abricate.tab",
                "amr_genes_found": 3,
                "gene_families": ["blaSHV", "aac"]
            }
        },
        "summary": {
            "total_genomes": 2,
            "genomes_with_amr": 2,
            "unique_gene_families": 4,
            "total_amr_genes": 8
        }
    }


@pytest.fixture
def sample_abricate_output():
    """Create sample ABRicate output for testing."""
    return """#FILE	SEQUENCE	START	END	STRAND	GENE	COVERAGE	COVERAGE_MAP	GAPS	%COVERAGE	%IDENTITY	DATABASE	ACCESSION	PRODUCT	RESISTANCE
NZ_CP107554.fasta	NZ_CP107554.1	1234567	1235890	+	blaOXA-1	1-324/324	===============	0/0	100.00	99.38	card	ARO:3000015	OXA beta-lactamase	BETA-LACTAM
NZ_CP107554.fasta	NZ_CP107554.1	2345678	2346234	+	aac(3)-IV	1-557/557	===============	0/0	100.00	98.92	card	ARO:3000456	aminoglycoside acetyltransferase	AMINOGLYCOSIDE
NZ_CP107555.fasta	NZ_CP107555.1	3456789	3457234	-	blaSHV-1	1-445/445	===============	0/0	100.00	99.55	card	ARO:3000789	SHV beta-lactamase	BETA-LACTAM"""


@pytest.fixture
def sample_mutation_data():
    """Create sample mutation data for testing."""
    return [
        {
            "Accession Number": "NZ_CP107554",
            "Gene Name": "blaOXA-1",
            "Gene Family": "blaOXA",
            "Position": 123,
            "Reference": "A",
            "Alternative": "G",
            "Mutation Type": "SNP",
            "Codon Change": "AAG>AGG",
            "AA Change": "K41R",
            "Frequency": 1
        },
        {
            "Accession Number": "NZ_CP107555",
            "Gene Name": "blaSHV-1",
            "Gene Family": "blaSHV",
            "Position": 456,
            "Reference": "C",
            "Alternative": "T",
            "Mutation Type": "SNP",
            "Codon Change": "CTG>TTG",
            "AA Change": "L152F",
            "Frequency": 1
        }
    ]


@pytest.fixture
def sample_fasta_content():
    """Create sample FASTA content for testing."""
    return """>NZ_CP107554.1 Escherichia coli strain test, complete genome
ATGAAAAAACTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
AAAAAACTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGTAG

>NZ_CP107555.1 Klebsiella pneumoniae strain test, complete genome
ATGAAAAAACTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
AAAAAACTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG
CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGTAG"""


@pytest.fixture
def mock_ncbi_response():
    """Create mock NCBI API response for testing."""
    return {
        "esearchresult": {
            "idlist": ["123456789", "987654321"],
            "count": "2"
        },
        "esummary": {
            "123456789": {
                "accessionversion": "NZ_CP107554.1",
                "organism": "Escherichia coli",
                "assemblyaccession": "GCF_123456789.1",
                "assemblylevel": "Complete Genome"
            },
            "987654321": {
                "accessionversion": "NZ_CP107555.1", 
                "organism": "Klebsiella pneumoniae",
                "assemblyaccession": "GCF_987654321.1",
                "assemblylevel": "Complete Genome"
            }
        }
    }


@pytest.fixture(autouse=True)
def setup_test_environment(monkeypatch):
    """Set up test environment variables."""
    # Disable external API calls during testing
    monkeypatch.setenv("MUTATIONSCAN_TESTING", "1")
    # Set test data paths
    test_dir = Path(__file__).parent
    monkeypatch.setenv("MUTATIONSCAN_TEST_DATA", str(test_dir / "data"))


def create_test_manifest_file(temp_dir: Path, manifest_data: Dict[str, Any], filename: str = "test_manifest.json") -> Path:
    """Helper function to create a test manifest file."""
    manifest_path = temp_dir / filename
    with open(manifest_path, 'w', encoding='utf-8') as f:
        json.dump(manifest_data, f, indent=2)
    return manifest_path


def create_test_fasta_file(temp_dir: Path, content: str, filename: str = "test.fasta") -> Path:
    """Helper function to create a test FASTA file."""
    fasta_path = temp_dir / filename
    with open(fasta_path, 'w', encoding='utf-8') as f:
        f.write(content)
    return fasta_path