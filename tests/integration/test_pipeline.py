"""Integration tests for pipeline workflows."""

import pytest
from mutation_scan.core import SequenceTranslator
from mutation_scan.analysis import SequenceAligner, VariantDetector

class TestSequenceToVariantPipeline:
    """Test complete sequence analysis pipeline."""

    def test_translate_and_align_workflow(self):
        """Test workflow: translate DNA -> align -> detect variants."""
        # DNA sequences (simplified for testing)
        dna1 = "ATGAAAGCG"
        dna2 = "ATGAAAGCGTAA"
        
        # Translate
        translator = SequenceTranslator()
        protein1 = translator.translate(dna1)
        protein2 = translator.translate(dna2)
        
        # Align
        aligner = SequenceAligner()
        alignment = aligner.global_alignment(protein1, protein2)
        
        # Results should be valid
        assert isinstance(alignment, dict)
        assert "score" in alignment
        assert "alignment1" in alignment

    def test_variant_detection_from_alignment(self):
        """Test variant detection from sequence alignment."""
        # Aligned sequences with intentional differences
        aligned_seq1 = "MKACDEFG"
        aligned_seq2 = "MKTCDEFG"
        
        detector = VariantDetector()
        variants = detector.detect_variants(aligned_seq1, aligned_seq2, ref_seq="seq1")
        
        assert isinstance(variants, list)
