"""
Integration tests for the refactored GenomeExtractor module.

Tests:
1. GenomeDownloader initialization
2. Accession file reading
3. Metadata parsing
4. GenomeProcessor validation
5. Error handling
"""

import json
import tempfile
import unittest
from io import BytesIO
from pathlib import Path
from unittest.mock import patch, MagicMock

# Note: Adjust imports based on your project structure
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from mutation_scan.genome_extractor import (
    NCBIDatasetsGenomeDownloader,
    GenomeProcessor,
)


class TestGenomeDownloader(unittest.TestCase):
    """Test NCBIDatasetsGenomeDownloader class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)
        self.log_file = self.output_dir / "test.log"

    def tearDown(self):
        """Clean up temporary files."""
        self.temp_dir.cleanup()

    def test_initialization(self):
        """Test downloader initialization."""
        downloader = NCBIDatasetsGenomeDownloader(
            email="test@example.com",
            api_key="test_key",
            output_dir=self.output_dir,
            log_file=self.log_file,
        )

        self.assertEqual(downloader.email, "test@example.com")
        self.assertEqual(downloader.api_key, "test_key")
        self.assertTrue(self.output_dir.exists())
        self.assertTrue(self.log_file.parent.exists())

    def test_initialization_without_email(self):
        """Test that initialization fails without email."""
        with self.assertRaises(ValueError):
            NCBIDatasetsGenomeDownloader(
                email="",
                output_dir=self.output_dir,
            )

    def test_read_accession_file(self):
        """Test reading accession file."""
        downloader = NCBIDatasetsGenomeDownloader(
            email="test@example.com",
            output_dir=self.output_dir,
        )

        # Create test file
        accession_file = self.output_dir / "accessions.txt"
        accessions_data = "GCF_000005845.2\nGCF_000007045.1\n# Comment line\nGCF_000009065.1\n\n"
        accession_file.write_text(accessions_data)

        # Read file
        accessions = downloader.read_accession_file(accession_file)

        self.assertEqual(len(accessions), 3)
        self.assertIn("GCF_000005845.2", accessions)
        self.assertIn("GCF_000007045.1", accessions)
        self.assertIn("GCF_000009065.1", accessions)

    def test_read_accession_file_not_found(self):
        """Test error when accession file not found."""
        downloader = NCBIDatasetsGenomeDownloader(
            email="test@example.com",
            output_dir=self.output_dir,
        )

        with self.assertRaises(FileNotFoundError):
            downloader.read_accession_file(Path("nonexistent.txt"))

    def test_parse_jsonl_metadata(self):
        """Test parsing JSONL metadata."""
        downloader = NCBIDatasetsGenomeDownloader(
            email="test@example.com",
            output_dir=self.output_dir,
        )

        # Create test JSONL
        jsonl_data = json.dumps({
            "organism_name": "Escherichia coli",
            "strain": "K-12",
            "collection_date": "2020-01-15",
            "host": "Human",
            "isolation_source": "Fecal",
            "country": "USA",
        })

        metadata = downloader._parse_jsonl_metadata(jsonl_data)

        self.assertEqual(metadata["Organism Name"], "Escherichia coli")
        self.assertEqual(metadata["Strain"], "K-12")
        self.assertEqual(metadata["Collection Date"], "2020-01-15")
        self.assertEqual(metadata["Host"], "Human")

    def test_parse_jsonl_missing_fields(self):
        """Test JSONL parsing with missing fields (anti-hallucination)."""
        downloader = NCBIDatasetsGenomeDownloader(
            email="test@example.com",
            output_dir=self.output_dir,
        )

        # JSONL with missing fields
        jsonl_data = json.dumps({
            "organism_name": "Escherichia coli",
            # strain, collection_date, host, isolation_source are missing
        })

        metadata = downloader._parse_jsonl_metadata(jsonl_data)

        # Defaults should be "N/A"
        self.assertEqual(metadata["Strain"], "N/A")
        self.assertEqual(metadata["Collection Date"], "N/A")
        self.assertEqual(metadata["Host"], "N/A")

    @patch("requests.get")
    def test_search_accessions(self, mock_get):
        """Test searching for accessions."""
        downloader = NCBIDatasetsGenomeDownloader(
            email="test@example.com",
            output_dir=self.output_dir,
        )

        # Mock esearch response
        search_response = MagicMock()
        search_response.json.return_value = {
            "esearchresult": {
                "idlist": ["12345", "12346"],
            }
        }

        # Mock esummary response
        summary_response = MagicMock()
        summary_response.json.return_value = {
            "result": {
                "assembly": [
                    {"assemblyaccession": "GCF_000005845.2"},
                    {"assemblyaccession": "GCF_000007045.1"},
                ]
            }
        }

        # Set up side effects
        mock_get.side_effect = [search_response, summary_response]

        # Search
        accessions = downloader.search_accessions("Escherichia coli", max_results=10)

        self.assertEqual(len(accessions), 2)
        self.assertIn("GCF_000005845.2", accessions)
        self.assertIn("GCF_000007045.1", accessions)


class TestGenomeProcessor(unittest.TestCase):
    """Test GenomeProcessor class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)

    def tearDown(self):
        """Clean up temporary files."""
        self.temp_dir.cleanup()

    def test_initialization(self):
        """Test processor initialization."""
        processor = GenomeProcessor(
            min_coverage=85.0,
            min_length=900000,
        )

        self.assertEqual(processor.min_coverage, 85.0)
        self.assertEqual(processor.min_length, 900000)

    def test_validate_valid_genome(self):
        """Test validation of valid genome."""
        processor = GenomeProcessor(min_coverage=90.0, min_length=1000)

        # Create test FASTA
        fasta_file = self.output_dir / "test.fasta"
        fasta_content = ">NC_000001\n" + "A" * 2000 + "\n"
        fasta_file.write_text(fasta_content)

        is_valid, message = processor.validate_genome(fasta_file)

        self.assertTrue(is_valid)
        self.assertIn("2000", message)

    def test_validate_invalid_genome_short(self):
        """Test validation of too-short genome."""
        processor = GenomeProcessor(min_coverage=90.0, min_length=1000000)

        # Create short FASTA
        fasta_file = self.output_dir / "short.fasta"
        fasta_content = ">NC_000001\n" + "A" * 1000 + "\n"
        fasta_file.write_text(fasta_content)

        is_valid, message = processor.validate_genome(fasta_file)

        self.assertFalse(is_valid)
        self.assertIn("too short", message.lower())

    def test_validate_invalid_file_not_found(self):
        """Test validation when file not found."""
        processor = GenomeProcessor()

        is_valid, message = processor.validate_genome(Path("nonexistent.fasta"))

        self.assertFalse(is_valid)
        self.assertIn("not found", message.lower())

    def test_calculate_coverage(self):
        """Test coverage calculation."""
        processor = GenomeProcessor()

        # Create test FASTA
        fasta_file = self.output_dir / "coverage.fasta"
        fasta_content = ">NC_000001\n" + "A" * 100 + "\n"
        fasta_file.write_text(fasta_content)

        # Without reference length
        coverage = processor.calculate_coverage(fasta_file)
        self.assertEqual(coverage, 100.0)

        # With reference length
        coverage = processor.calculate_coverage(fasta_file, reference_length=100)
        self.assertEqual(coverage, 100.0)

    def test_extract_metadata(self):
        """Test metadata extraction from FASTA."""
        processor = GenomeProcessor()

        # Create test FASTA
        fasta_file = self.output_dir / "metadata.fasta"
        fasta_content = (
            ">NC_000001 [organism=Escherichia coli] [strain=K-12]\n"
            + "A" * 1000 + "\n"
            + ">NC_000002\n"
            + "G" * 500 + "\n"
        )
        fasta_file.write_text(fasta_content)

        metadata = processor.extract_metadata(fasta_file)

        self.assertEqual(metadata["sequences"], 2)
        self.assertIn("NC_000001", metadata["first_header"])
        self.assertGreater(metadata["file_size"], 0)


class TestIntegration(unittest.TestCase):
    """Integration tests for full workflow."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.TemporaryDirectory()
        self.output_dir = Path(self.temp_dir.name)

    def tearDown(self):
        """Clean up temporary files."""
        self.temp_dir.cleanup()

    def test_full_workflow(self):
        """Test full download and validation workflow."""
        downloader = NCBIDatasetsGenomeDownloader(
            email="test@example.com",
            output_dir=self.output_dir,
        )

        processor = GenomeProcessor()

        # Create test accession file
        accession_file = self.output_dir / "accessions.txt"
        accession_file.write_text("GCF_000005845.2\nGCF_000007045.1\n")

        # Read accessions
        accessions = downloader.read_accession_file(accession_file)
        self.assertEqual(len(accessions), 2)

        # Validate (would normally download, but we test the validation part)
        # Create a mock FASTA for validation
        fasta_file = self.output_dir / "GCF_000005845.2.fasta"
        fasta_file.write_text(">NC_000001\n" + "A" * 1000000 + "\n")

        is_valid, message = processor.validate_genome(fasta_file)
        self.assertTrue(is_valid)


if __name__ == "__main__":
    unittest.main()
