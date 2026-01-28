"""Unit tests for sequence extractor module."""

import pytest
from mutation_scan.sequence_extractor import SequenceTranslator

class TestSequenceTranslator:
    """Test suite for SequenceTranslator class."""

    def setup_method(self):
        """Setup test fixtures."""
        self.translator = SequenceTranslator()

    def test_translate_standard_sequence(self, sample_dna_sequence):
        """Test translation of standard DNA sequence."""
        # Expected: ATG -> M, AAA -> K, GCG -> A, TAA -> *, CGT -> R, AAA -> K
        protein = self.translator.translate(sample_dna_sequence)
        assert isinstance(protein, str)
        assert len(protein) > 0
        assert protein[0] == "M"  # Start codon

    def test_translate_with_frame_zero(self, sample_dna_sequence):
        """Test translation with frame 0 (default)."""
        protein = self.translator.translate(sample_dna_sequence, frame=0)
        assert isinstance(protein, str)

    def test_translate_with_frame_one(self, sample_dna_sequence):
        """Test translation with frame 1."""
        protein = self.translator.translate(sample_dna_sequence, frame=1)
        assert isinstance(protein, str)

    def test_translate_with_frame_two(self, sample_dna_sequence):
        """Test translation with frame 2."""
        protein = self.translator.translate(sample_dna_sequence, frame=2)
        assert isinstance(protein, str)

    def test_translate_empty_sequence(self):
        """Test translation of empty sequence."""
        protein = self.translator.translate("")
        assert protein == ""

    def test_translate_lowercase_input(self):
        """Test that lowercase DNA is handled."""
        dna = "atgaaagcg"
        protein = self.translator.translate(dna)
        assert isinstance(protein, str)

    def test_genetic_code_completeness(self):
        """Test that genetic code is complete."""
        assert len(self.translator.STANDARD_GENETIC_CODE) == 64
