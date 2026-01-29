"""
Variant Calling module for MutationScan.

This module identifies point mutations (substitutions) by performing global alignment
between query proteins (from SequenceExtractor) and wild-type reference proteins.

CRITICAL MUTATION CALLING LOGIC:
- Uses Bio.Align.PairwiseAligner (Python-native, no external dependencies like EMBOSS Needle)
- Global alignment with BLOSUM62 matrix
- "Residue Counter Algorithm": Position counter increments ONLY when reference is NOT a gap
- Mutation format: {Ref_AA}{Position}{Query_AA} (e.g., S83L means Serine→Leucine at position 83)
- Position 83 = 83rd amino acid in reference, NOT 83rd character in alignment

ANTI-HALLUCINATION RULES:
- Never count gaps as positions
- Never crash on partial proteins (align what you have)
- Status = "Resistant" if mutation in resistance_db.json
- If not in DB, optionally route to ML predictor (Module 6)
- Default phenotype to "N/A" if not in database or prediction fails
- Always label prediction_source in output for transparency
"""

import importlib
import inspect
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import json

import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


class VariantCaller:
    """
    Identify point mutations by aligning query proteins against wild-type references.
    
    This class takes protein FASTA files from SequenceExtractor and compares them
    against wild-type references to identify substitution mutations (SNPs at protein level).
    
    Key Features:
    - Python-native global alignment (Bio.Align.PairwiseAligner)
    - BLOSUM62 scoring matrix (industry standard for protein alignment)
    - Residue counter algorithm (gap-aware position tracking)
    - Resistance interpretation via resistance_db.json
    - CSV output: Accession, Gene, Mutation, Status, Phenotype, Reference_PDB
    - Helper: _generate_dummy_references() for instant testing
    """

    def __init__(
        self,
        refs_dir: Path,
        resistance_db_path: Optional[Path] = None,
        enable_ml: bool = True,
        ml_models_dir: Optional[Path] = None,
        antibiotic: str = "Ciprofloxacin"
    ):
        """
        Initialize VariantCaller.

        Args:
            refs_dir: Path to directory containing wild-type reference FASTA files
                     Format: {GeneName}_WT.faa (e.g., gyrA_WT.faa)
            resistance_db_path: Path to resistance_db.json (optional)
                              If None, defaults to data/refs/resistance_db.json
            enable_ml: If True, use Module 6 ML predictor for unknown mutations
            ml_models_dir: Optional path to ML model directory
            antibiotic: Antibiotic name passed to ML predictor (default: Ciprofloxacin)

        Raises:
            FileNotFoundError: If refs_dir doesn't exist
        """
        self.refs_dir = Path(refs_dir)
        
        if not self.refs_dir.exists():
            raise FileNotFoundError(f"References directory not found: {self.refs_dir}")
        
        if not self.refs_dir.is_dir():
            raise NotADirectoryError(f"Not a directory: {self.refs_dir}")
        
        # Set resistance DB path
        if resistance_db_path is None:
            self.resistance_db_path = self.refs_dir / "resistance_db.json"
        else:
            self.resistance_db_path = Path(resistance_db_path)
        
        # Load resistance database
        self.resistance_db = self._load_resistance_db()

        # ML predictor settings (Module 6)
        self.enable_ml = enable_ml
        self.ml_models_dir = Path(ml_models_dir) if ml_models_dir else None
        self.antibiotic = antibiotic
        self._ml_predictor = None
        self._ml_predictor_error: Optional[Exception] = None
        
        # Initialize PairwiseAligner with BLOSUM62
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'  # Global alignment (Needleman-Wunsch)
        self.aligner.substitution_matrix = self._get_blosum62()
        
        # Set gap penalties to prefer substitutions over indels
        # This ensures S83L is detected as substitution, not indel
        self.aligner.open_gap_score = -10.0
        self.aligner.extend_gap_score = -0.5
        
        logger.info(f"Initialized VariantCaller with refs_dir: {self.refs_dir}")
        logger.info(f"Loaded {len(self.resistance_db)} genes in resistance database")

    def _get_blosum62(self):
        """
        Get BLOSUM62 substitution matrix from Biopython.
        
        Returns:
            Substitution matrix for PairwiseAligner
        """
        try:
            from Bio.Align import substitution_matrices
            blosum62 = substitution_matrices.load("BLOSUM62")
            logger.debug("Loaded BLOSUM62 matrix")
            return blosum62
        except Exception as e:
            logger.warning(f"Failed to load BLOSUM62: {e}. Using default scoring.")
            return None

    def _load_resistance_db(self) -> Dict:
        """
        Load resistance database from JSON file.
        
        Returns:
            Dictionary mapping gene -> list of mutation dicts
            Format: {
                "gyrA": [
                    {"mutation": "S83L", "phenotype": "Fluoroquinolone resistance", "pdb": "3NUU"},
                    {"mutation": "D87N", "phenotype": "Fluoroquinolone resistance", "pdb": "3NUU"}
                ]
            }
        """
        if not self.resistance_db_path.exists():
            logger.warning(f"Resistance DB not found: {self.resistance_db_path}. Returning empty DB.")
            return {}
        
        try:
            with open(self.resistance_db_path, 'r') as f:
                db = json.load(f)
            logger.info(f"Loaded resistance database: {self.resistance_db_path}")
            return db
        except Exception as e:
            logger.error(f"Failed to load resistance DB: {e}")
            return {}

    def call_variants(
        self,
        proteins_dir: Path,
        output_csv: Path
    ) -> pd.DataFrame:
        """
        Process all .faa files in proteins_dir, call mutations, and write CSV report.

        Args:
            proteins_dir: Directory containing .faa files from SequenceExtractor
                         Format: ACCESSION_GeneName.faa
            output_csv: Path to output CSV file

        Returns:
            DataFrame with columns:
            Accession, Gene, Mutation, Status, Phenotype, Reference_PDB,
            prediction_score, prediction_source

        Example:
            Input: data/proteins/GCF_001_gyrA.faa
            Output: mutation_report.csv with rows like:
                    Accession,Gene,Mutation,Status,Phenotype,Reference_PDB
                    GCF_001,gyrA,S83L,Resistant,Fluoroquinolone resistance,3NUU
        """
        proteins_dir = Path(proteins_dir)
        output_csv = Path(output_csv)
        
        if not proteins_dir.exists():
            raise FileNotFoundError(f"Proteins directory not found: {proteins_dir}")
        
        # Find all .faa files
        faa_files = list(proteins_dir.glob("*.faa"))
        
        if not faa_files:
            logger.warning(f"No .faa files found in {proteins_dir}")
            # Return empty DataFrame
            df = pd.DataFrame(columns=[
                'Accession', 'Gene', 'Mutation', 'Status', 'Phenotype', 'Reference_PDB',
                'prediction_score', 'prediction_source'
            ])
            df.to_csv(output_csv, index=False)
            return df
        
        logger.info(f"Processing {len(faa_files)} protein files")
        
        all_mutations = []
        
        for faa_file in faa_files:
            # Parse filename: ACCESSION_GeneName.faa
            filename = faa_file.stem  # Remove .faa extension
            
            try:
                # Expected format: ACCESSION_GeneName
                parts = filename.split('_', 1)  # Split on first underscore only
                
                if len(parts) != 2:
                    logger.warning(f"Skipping invalid filename format: {faa_file.name} (expected ACCESSION_GeneName.faa)")
                    continue
                
                accession, gene_name = parts
                
                # Call variants for this protein
                mutations = self._call_variants_single(faa_file, accession, gene_name)
                all_mutations.extend(mutations)
                
            except Exception as e:
                logger.error(f"Failed to process {faa_file.name}: {e}")
                continue
        
        # Create DataFrame
        df = pd.DataFrame(all_mutations)
        
        # Ensure columns exist (even if empty)
        if df.empty:
            df = pd.DataFrame(columns=[
                'Accession', 'Gene', 'Mutation', 'Status', 'Phenotype', 'Reference_PDB',
                'prediction_score', 'prediction_source'
            ])
        
        # Save to CSV
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_csv, index=False)
        
        logger.info(f"Saved mutation report: {output_csv}")
        logger.info(f"Total mutations identified: {len(df)}")
        
        # Summary statistics
        if not df.empty:
            resistant_count = (df['Status'] == 'Resistant').sum()
            vus_count = (df['Status'] == 'VUS').sum()
            logger.info(f"Resistant mutations: {resistant_count}, VUS: {vus_count}")
        
        return df

    def _call_variants_single(
        self,
        faa_file: Path,
        accession: str,
        gene_name: str
    ) -> List[Dict]:
        """
        Call variants for a single protein file.

        Args:
            faa_file: Path to .faa file
            accession: Genome accession
            gene_name: Gene name

        Returns:
            List of mutation dictionaries
        """
        # Load wild-type reference
        ref_file = self.refs_dir / f"{gene_name}_WT.faa"
        
        if not ref_file.exists():
            logger.warning(f"Wild-type reference not found for {gene_name}: {ref_file}. Skipping.")
            return []
        
        try:
            # Load reference sequence
            ref_record = SeqIO.read(ref_file, "fasta")
            
            # Load query sequence
            query_record = SeqIO.read(faa_file, "fasta")
            
            # Perform global alignment
            mutations = self._align_and_call_mutations(
                query_record,
                ref_record,
                accession,
                gene_name
            )
            
            return mutations
            
        except Exception as e:
            logger.error(f"Error calling variants for {faa_file.name}: {e}")
            return []

    def _align_and_call_mutations(
        self,
        query: SeqRecord,
        reference: SeqRecord,
        accession: str,
        gene_name: str
    ) -> List[Dict]:
        """
        Perform global alignment and call mutations using the Residue Counter Algorithm.

        CRITICAL: The "Residue Counter Algorithm"
        ----------------------------------------
        - Iterate through aligned sequences character-by-character
        - Maintain reference_position counter (starts at 1)
        - Rule: Increment counter ONLY if reference character is NOT a gap (-)
        - If Ref_AA != Query_AA AND Ref_AA != '-' AND Query_AA != '-':
            -> Record substitution: {Ref_AA}{reference_position}{Query_AA}

        Example:
            Reference:  M  K  T  -  I  A
            Query:      M  R  T  V  I  A
            Position:   1  2  3     4  5

            Mutations:
            - K2R (position 2: K→R)
            - Position 4 skipped because reference has gap
            - No mutation at position 1, 3, 4, 5 (identical)

        Args:
            query: Query protein sequence
            reference: Wild-type reference sequence
            accession: Genome accession
            gene_name: Gene name

        Returns:
            List of mutation dictionaries with keys:
            - Accession
            - Gene
            - Mutation (e.g., "S83L")
            - Status ("Resistant", or ML-predicted risk)
            - Phenotype (e.g., "Fluoroquinolone resistance" or "N/A")
            - Reference_PDB (e.g., "3NUU" or "N/A")
            - prediction_score (0.0 - 1.0)
            - prediction_source ("Clinical DB" or "AI Model")
        """
        try:
            # Perform global alignment
            # CRITICAL: Order matters! align(reference.seq, query.seq) NOT align(query.seq, reference.seq)
            # alignment[0] = reference, alignment[1] = query
            alignments = list(self.aligner.align(reference.seq, query.seq))
            
            if not alignments:
                logger.warning(f"No alignment found for {accession}_{gene_name}")
                return []
            
            # Take best alignment (first one, highest score)
            alignment = alignments[0]
            
            # Extract aligned sequences
            # alignment[0] = reference (target in Biopython terminology)
            # alignment[1] = query
            aligned_ref = str(alignment[0])    # Reference sequence with gaps
            aligned_query = str(alignment[1])  # Query sequence with gaps
            
            # Sanity check
            if len(aligned_query) != len(aligned_ref):
                logger.error(f"Alignment length mismatch for {accession}_{gene_name}")
                return []
            
            # Log alignment quality
            score = alignment.score
            identities = sum(1 for a, b in zip(aligned_query, aligned_ref) if a == b and a != '-')
            length = len(aligned_ref.replace('-', ''))  # Reference length without gaps
            identity_percent = (identities / length) * 100 if length > 0 else 0
            
            logger.info(f"Alignment {accession}_{gene_name}: Score={score:.1f}, Identity={identity_percent:.1f}%")
            
            # Call mutations using Residue Counter Algorithm
            mutations = []
            reference_position = 0  # Will increment to 1 on first non-gap residue
            
            for i in range(len(aligned_ref)):
                ref_aa = aligned_ref[i]
                query_aa = aligned_query[i]
                
                # Increment position counter ONLY if reference is not a gap
                if ref_aa != '-':
                    reference_position += 1
                
                # Check for substitution
                # Conditions:
                # 1. Reference and Query are different
                # 2. Reference is NOT a gap (we have a real position)
                # 3. Query is NOT a gap (we have a real amino acid)
                if ref_aa != query_aa and ref_aa != '-' and query_aa != '-':
                    mutation_str = f"{ref_aa}{reference_position}{query_aa}"
                    
                    # Interpret mutation
                    status, phenotype, pdb, prediction_score, prediction_source = self._interpret_mutation(
                        gene_name,
                        mutation_str
                    )
                    
                    mutations.append({
                        'Accession': accession,
                        'Gene': gene_name,
                        'Mutation': mutation_str,
                        'Status': status,
                        'Phenotype': phenotype,
                        'Reference_PDB': pdb,
                        'prediction_score': prediction_score,
                        'prediction_source': prediction_source
                    })
                    
                    logger.debug(f"Found mutation: {mutation_str} ({status})")
            
            logger.info(f"Identified {len(mutations)} mutations in {accession}_{gene_name}")
            
            return mutations
            
        except Exception as e:
            logger.error(f"Alignment failed for {accession}_{gene_name}: {e}")
            return []

    def _interpret_mutation(
        self,
        gene_name: str,
        mutation: str
    ) -> Tuple[str, str, str, Optional[float], str]:
        """
        Interpret mutation using resistance_db.json.

        Args:
            gene_name: Gene name (e.g., "gyrA")
            mutation: Mutation string (e.g., "S83L")

        Returns:
            Tuple of (status, phenotype, pdb_id, prediction_score, prediction_source)
            - status: "Resistant" if in DB, else ML prediction or "Unknown"
            - phenotype: Description from DB, else "N/A"
            - pdb_id: PDB ID from DB, else "N/A"
            - prediction_score: 0.0-1.0 for ML predictions, 1.0 for DB hits
            - prediction_source: "Clinical DB" or "AI Model"
        """
        # Check if gene exists in resistance DB
        if gene_name not in self.resistance_db:
            return self._fallback_to_ml(mutation)
        
        # Check if mutation exists for this gene
        gene_mutations = self.resistance_db[gene_name]
        
        for entry in gene_mutations:
            if entry.get('mutation') == mutation:
                return (
                    "Resistant",
                    entry.get('phenotype', 'N/A'),
                    entry.get('pdb', 'N/A'),
                    1.0,
                    "Clinical DB"
                )
        
        # Mutation not in database
        return self._fallback_to_ml(mutation)

    def _fallback_to_ml(self, mutation: str) -> Tuple[str, str, str, Optional[float], str]:
        """
        Fallback to ML predictor (Module 6) for unknown mutations.

        Args:
            mutation: Mutation string (e.g., "S83L")

        Returns:
            Tuple of (status, phenotype, pdb_id, prediction_score, prediction_source)
        """
        if not self.enable_ml:
            return ("VUS", "N/A", "N/A", None, "Clinical DB")

        prediction = self._predict_with_ml(mutation)

        if prediction.get("success"):
            risk_level = prediction.get("risk_level", "Unknown")
            resistance_prob = prediction.get("resistance_prob")

            return (
                f"Predicted {risk_level} Risk",
                f"Predicted resistance risk for {self.antibiotic}",
                "N/A",
                resistance_prob,
                "AI Model"
            )

        return (
            "Unknown (Parse Failed)",
            "N/A",
            "N/A",
            None,
            "AI Model"
        )

    def _get_ml_predictor(self):
        """
        Lazily import and initialize the ML predictor (Module 6).

        Returns:
            ResistancePredictor instance or None if unavailable
        """
        if not self.enable_ml:
            return None

        if self._ml_predictor is not None:
            return self._ml_predictor

        if self._ml_predictor_error is not None:
            return None

        try:
            module = importlib.import_module("mutation_scan.ml_predictor.inference")
            predictor_cls = getattr(module, "ResistancePredictor")
        except Exception as e:
            self._ml_predictor_error = e
            logger.warning(f"ML predictor unavailable: {e}")
            return None

        try:
            self._ml_predictor = self._init_ml_predictor(predictor_cls)
            return self._ml_predictor
        except Exception as e:
            self._ml_predictor_error = e
            logger.warning(f"Failed to initialize ML predictor: {e}")
            return None

    def _init_ml_predictor(self, predictor_cls):
        """
        Initialize ML predictor, attempting common constructor parameter names.
        """
        if self.ml_models_dir:
            try:
                sig = inspect.signature(predictor_cls)
                if "model_dir" in sig.parameters:
                    return predictor_cls(model_dir=self.ml_models_dir)
                if "models_dir" in sig.parameters:
                    return predictor_cls(models_dir=self.ml_models_dir)
                if "model_path" in sig.parameters:
                    return predictor_cls(model_path=self.ml_models_dir)
            except (TypeError, ValueError):
                pass

        return predictor_cls()

    def _predict_with_ml(self, mutation: str) -> Dict:
        """
        Run ML prediction for an unknown mutation.

        Returns:
            Dict with keys: success, resistance_prob, risk_level
        """
        predictor = self._get_ml_predictor()

        if predictor is None:
            return {"success": False, "error": "ML predictor unavailable"}

        try:
            try:
                result = predictor.predict(mutation, antibiotic=self.antibiotic)
            except TypeError:
                result = predictor.predict(mutation)

            if not isinstance(result, dict):
                return {"success": False, "error": "Unexpected predictor output"}

            success = result.get("success", True)
            resistance_prob = result.get("resistance_prob")
            risk_level = result.get("risk_level")

            if success and resistance_prob is not None and risk_level:
                return {
                    "success": True,
                    "resistance_prob": resistance_prob,
                    "risk_level": risk_level
                }

            return {"success": False, "error": "Missing prediction fields"}

        except Exception as e:
            logger.warning(f"ML prediction failed for {mutation}: {e}")
            return {"success": False, "error": str(e)}

    def _generate_dummy_references(self) -> None:
        """
        Generate placeholder E. coli K12 wild-type references for testing.
        
        This helper creates a minimal gyrA_WT.faa file so users can test
        the module immediately without needing real references.
        
        Creates:
        - data/refs/gyrA_WT.faa (E. coli K12 GyrA wild-type)
        
        Note: This is for TESTING ONLY. Production use requires real references.
        """
        # E. coli K12 GyrA protein (NCBI Reference: NP_414550.1)
        # First 100 amino acids (truncated for demo purposes)
        gyrA_wt_seq = (
            "MSDLAREITPVNIEEELKSSYLDYAMSVIVGRALPDVRDGLKPVHRRVLYAMNVLGND"
            "WNKAYKKSARVVGDVIGKYHPHGDSA"
        )
        
        gyrA_record = SeqRecord(
            Seq(gyrA_wt_seq),
            id="gyrA_WT_K12",
            description="E. coli K12 GyrA wild-type (partial sequence for testing)"
        )
        
        # Create refs directory if it doesn't exist
        self.refs_dir.mkdir(parents=True, exist_ok=True)
        
        # Write gyrA_WT.faa
        gyrA_file = self.refs_dir / "gyrA_WT.faa"
        SeqIO.write(gyrA_record, gyrA_file, "fasta")
        
        logger.info(f"Generated dummy reference: {gyrA_file}")
        
        # Also create a minimal resistance_db.json
        dummy_db = {
            "gyrA": [
                {"mutation": "S83L", "phenotype": "Fluoroquinolone resistance", "pdb": "3NUU"},
                {"mutation": "D87N", "phenotype": "Fluoroquinolone resistance", "pdb": "3NUU"},
                {"mutation": "D87G", "phenotype": "Fluoroquinolone resistance", "pdb": "3NUU"}
            ],
            "parC": [
                {"mutation": "S80I", "phenotype": "Fluoroquinolone resistance", "pdb": "N/A"},
                {"mutation": "E84K", "phenotype": "Fluoroquinolone resistance", "pdb": "N/A"}
            ]
        }
        
        db_file = self.refs_dir / "resistance_db.json"
        with open(db_file, 'w') as f:
            json.dump(dummy_db, f, indent=2)
        
        logger.info(f"Generated dummy resistance DB: {db_file}")
        logger.info("Dummy references ready for testing!")

    def get_available_references(self) -> List[str]:
        """
        Get list of available wild-type references.

        Returns:
            List of gene names (without _WT.faa suffix)
        """
        ref_files = list(self.refs_dir.glob("*_WT.faa"))
        gene_names = [f.stem.replace('_WT', '') for f in ref_files]
        
        logger.info(f"Found {len(gene_names)} wild-type references: {gene_names}")
        return gene_names

    def get_mutation_summary(self, mutations_df: pd.DataFrame) -> Dict:
        """
        Generate summary statistics from mutation report.

        Args:
            mutations_df: DataFrame from call_variants()

        Returns:
            Dictionary with summary statistics
        """
        if mutations_df.empty:
            return {
                'total_mutations': 0,
                'resistant_mutations': 0,
                'vus_mutations': 0,
                'unique_genes': 0,
                'unique_accessions': 0
            }
        
        summary = {
            'total_mutations': len(mutations_df),
            'resistant_mutations': (mutations_df['Status'] == 'Resistant').sum(),
            'vus_mutations': (mutations_df['Status'] == 'VUS').sum(),
            'unique_genes': mutations_df['Gene'].nunique(),
            'unique_accessions': mutations_df['Accession'].nunique(),
            'mutations_by_gene': mutations_df['Gene'].value_counts().to_dict()
        }
        
        return summary


# =============================================================================
# USAGE EXAMPLE
# =============================================================================
if __name__ == "__main__":
    """
    Example usage of VariantCaller.
    
    Typical workflow:
    1. Run GenomeExtractor → Download genomes to data/genomes/
    2. Run GeneFinder → Generate genes.csv with coordinates
    3. Run SequenceExtractor → Extract and translate sequences to data/proteins/
    4. Run VariantCaller → Compare against wild-type, identify mutations
    """
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Initialize VariantCaller
    caller = VariantCaller(
        refs_dir=Path("data/refs")
    )
    
    # Generate dummy references for testing (FIRST TIME ONLY)
    # Uncomment this if you need to create test references
    # caller._generate_dummy_references()
    
    # Check available references
    refs = caller.get_available_references()
    print(f"\nAvailable wild-type references: {refs}")
    
    # Call variants
    mutations_df = caller.call_variants(
        proteins_dir=Path("data/proteins"),
        output_csv=Path("data/results/mutation_report.csv")
    )
    
    # Print summary
    print("\n=== MUTATION REPORT SUMMARY ===")
    print(mutations_df.to_string(index=False))
    
    # Get statistics
    summary = caller.get_mutation_summary(mutations_df)
    print("\n=== STATISTICS ===")
    print(f"Total mutations: {summary['total_mutations']}")
    print(f"Resistant: {summary['resistant_mutations']}")
    print(f"VUS: {summary['vus_mutations']}")
    print(f"Unique genes: {summary['unique_genes']}")
    print(f"Unique accessions: {summary['unique_accessions']}")
    
    if 'mutations_by_gene' in summary:
        print("\nMutations by gene:")
        for gene, count in summary['mutations_by_gene'].items():
            print(f"  {gene}: {count}")
