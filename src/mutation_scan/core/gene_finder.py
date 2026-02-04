"""
Hybrid Gene Detection module using ABRicate and BLASTn.

This module provides functionality to detect both resistance genes (via ABRicate)
and housekeeping genes (via local BLASTn) from bacterial genome sequences.
"""

import logging
import subprocess
import tempfile
from io import StringIO
from pathlib import Path
from typing import Optional, List

import pandas as pd

logger = logging.getLogger(__name__)


class GeneFinder:
    """
    Hybrid gene detection using ABRicate (resistance) and BLASTn (housekeeping).
    
    Returns standardized DataFrames with columns:
    - Gene: Gene name (e.g., 'gyrA', 'blaNDM-1')
    - Contig: Contig ID where gene was found
    - Start: Start position (integer)
    - End: End position (integer)
    - Strand: '+' or '-'
    - Identity: Percent identity (float)
    - Source: 'ABRicate' or 'BLAST'
    """

    @staticmethod
    def load_target_genes(target_file: Path) -> List[str]:
        """
        Load target gene names from a text file.
        
        Args:
            target_file: Path to target genes file (one gene per line, # for comments)
        
        Returns:
            List of gene names (case-insensitive, lowercased)
        
        Raises:
            FileNotFoundError: If target file doesn't exist
        """
        target_file = Path(target_file)
        
        if not target_file.exists():
            logger.error(f"Target genes file not found: {target_file}")
            raise FileNotFoundError(f"Target genes file not found: {target_file}")
        
        genes = []
        with open(target_file, 'r') as f:
            for line in f:
                line = line.strip()
                # Skip empty lines and comments
                if line and not line.startswith('#'):
                    genes.append(line.lower())
        
        logger.info(f"Loaded {len(genes)} target genes from {target_file.name}")
        return genes

    def __init__(
        self, 
        abricate_db: str = "card",
        target_genes: Optional[List[str]] = None
    ):
        """
        Initialize GeneFinder.

        Args:
            abricate_db: ABRicate database to use (default: 'card')
            target_genes: Optional list of gene names to filter results (case-insensitive)
                         If None, all detected genes are returned.

        Raises:
            EnvironmentError: If ABRicate is not installed
        """
        self.abricate_db = abricate_db
        self.target_genes = [g.lower() for g in target_genes] if target_genes else None
        self.abricate_path = self._find_abricate()
        
        if not self.abricate_path:
            raise EnvironmentError(
                "ABRicate not found. Please install ABRicate and ensure it's in PATH.\n"
                "Install instructions: https://github.com/tseemann/abricate"
            )
        
        logger.info(f"Initialized GeneFinder with ABRicate: {self.abricate_path}")
        logger.info(f"Using database: {self.abricate_db}")
        if self.target_genes:
            logger.info(f"Filtering for {len(self.target_genes)} target genes: {', '.join(self.target_genes)}")

    def _find_abricate(self) -> Optional[str]:
        """
        Find ABRicate executable in system PATH.

        Returns:
            Path to abricate executable, or None if not found
        """
        possible_paths = [
            'abricate',  # In PATH
            'wsl abricate',  # WSL installation
            '/usr/local/bin/abricate',  # System install
            Path.home() / 'miniconda3' / 'bin' / 'abricate',  # Miniconda
            Path.home() / 'miniconda3' / 'envs' / 'abricate-env' / 'bin' / 'abricate',  # Conda env
            Path.home() / 'anaconda3' / 'bin' / 'abricate',  # Anaconda
        ]

        for path in possible_paths:
            try:
                path_str = str(path) if isinstance(path, Path) else path
                cmd = path_str.split() + ['--version']
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    logger.debug(f"Found ABRicate at: {path_str}")
                    return path_str
            except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
                continue

        return None

    def find_resistance_genes(self, fasta_file: Path) -> pd.DataFrame:
        """
        Find antibiotic resistance genes using ABRicate.

        Args:
            fasta_file: Path to input FASTA file

        Returns:
            DataFrame with standardized columns (Gene, Contig, Start, End, Strand, Identity, Source)
            Returns empty DataFrame on error

        Raises:
            FileNotFoundError: If FASTA file doesn't exist
        """
        fasta_file = Path(fasta_file)
        
        if not fasta_file.exists():
            logger.error(f"FASTA file not found: {fasta_file}")
            raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
        
        if fasta_file.stat().st_size == 0:
            logger.warning(f"Empty FASTA file: {fasta_file}")
            return self._empty_dataframe()
        
        logger.info(f"Running ABRicate on {fasta_file.name}")
        
        try:
            # Run ABRicate
            cmd = self.abricate_path.split() + ['--db', self.abricate_db, str(fasta_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode != 0:
                logger.error(f"ABRicate failed: {result.stderr.strip()}")
                return self._empty_dataframe()
            
            # Parse output
            return self._parse_abricate_output(result.stdout)
            
        except subprocess.TimeoutExpired:
            logger.error(f"ABRicate timed out for {fasta_file}")
            return self._empty_dataframe()
        except Exception as e:
            logger.error(f"Error running ABRicate: {e}")
            return self._empty_dataframe()

    def _parse_abricate_output(self, output: str) -> pd.DataFrame:
        """
        Parse ABRicate tab-delimited output to standardized DataFrame.

        ABRicate columns: #FILE, SEQUENCE, START, END, STRAND, GENE, COVERAGE, 
                         COVERAGE_MAP, GAPS, %COVERAGE, %IDENTITY, DATABASE, ACCESSION, PRODUCT, RESISTANCE

        Args:
            output: Raw ABRicate output text

        Returns:
            Standardized DataFrame
        """
        if not output or output.strip().startswith('#FILE') and output.count('\n') <= 1:
            logger.debug("No resistance genes found")
            return self._empty_dataframe()
        
        try:
            # Parse tab-delimited output (preserve header starting with #FILE)
            df = pd.read_csv(StringIO(output), sep='\t')
            # Strip leading '#' from header columns (e.g., #FILE -> FILE)
            df = df.rename(columns={col: col.lstrip('#') for col in df.columns})
            
            if df.empty:
                return self._empty_dataframe()
            
            # Standardize column names (ABRicate format)
            column_mapping = {
                'GENE': 'Gene',
                'SEQUENCE': 'Contig',
                'START': 'Start',
                'END': 'End',
                'STRAND': 'Strand',
                '%IDENTITY': 'Identity',
            }
            
            # Rename columns
            df = df.rename(columns=column_mapping)
            
            # Add source
            df['Source'] = 'ABRicate'
            
            # Select and reorder standardized columns
            standard_cols = ['Gene', 'Contig', 'Start', 'End', 'Strand', 'Identity', 'Source']
            df = df[standard_cols]
            
            # Ensure correct data types
            df['Start'] = df['Start'].astype(int)
            df['End'] = df['End'].astype(int)
            df['Identity'] = df['Identity'].astype(float)
            
            # Filter by target genes if specified
            if self.target_genes:
                before_count = len(df)
                df = df[df['Gene'].str.lower().isin(self.target_genes)]
                logger.info(f"Filtered {before_count} genes to {len(df)} target genes")
            else:
                logger.info(f"Found {len(df)} resistance genes")
            
            return df
            
        except Exception as e:
            logger.error(f"Error parsing ABRicate output: {e}")
            return self._empty_dataframe()

    def find_housekeeping_genes(
        self, 
        fasta_file: Path, 
        reference_db: Path,
        min_identity: float = 80.0,
        min_coverage: float = 90.0
    ) -> pd.DataFrame:
        """
        Find housekeeping genes using local BLASTn.

        Args:
            fasta_file: Path to input FASTA file (query)
            reference_db: Path to reference database (subject FASTA)
            min_identity: Minimum percent identity threshold (default: 80%)
            min_coverage: Minimum coverage threshold (default: 90%)

        Returns:
            DataFrame with standardized columns (Gene, Contig, Start, End, Strand, Identity, Source)
            Returns empty DataFrame on error

        Raises:
            FileNotFoundError: If files don't exist
            EnvironmentError: If blastn is not installed
        """
        fasta_file = Path(fasta_file)
        reference_db = Path(reference_db)
        
        if not fasta_file.exists():
            logger.error(f"Query FASTA not found: {fasta_file}")
            raise FileNotFoundError(f"Query FASTA not found: {fasta_file}")
        
        if not reference_db.exists():
            logger.error(f"Reference database not found: {reference_db}")
            raise FileNotFoundError(f"Reference database not found: {reference_db}")
        
        # Check if blastn is available
        if not self._check_blastn():
            raise EnvironmentError(
                "blastn not found. Please install BLAST+ and ensure it's in PATH.\n"
                "Install instructions: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download"
            )
        
        logger.info(f"Running BLASTn for housekeeping genes on {fasta_file.name}")
        
        try:
            # Run blastn with tabular output (outfmt 6)
            # Added 'slen' (Subject Length) for true coverage calculation
            cmd = [
                'blastn',
                '-query', str(fasta_file),
                '-subject', str(reference_db),
                '-outfmt', '6 qseqid sseqid pident length qstart qend sstart send sstrand slen',
                '-evalue', '1e-5',
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            
            if result.returncode != 0:
                logger.error(f"BLASTn failed: {result.stderr.strip()}")
                return self._empty_dataframe()
            
            # Parse output
            return self._parse_blastn_output(
                result.stdout, 
                min_identity=min_identity, 
                min_coverage=min_coverage
            )
            
        except subprocess.TimeoutExpired:
            logger.error(f"BLASTn timed out for {fasta_file}")
            return self._empty_dataframe()
        except Exception as e:
            logger.error(f"Error running BLASTn: {e}")
            return self._empty_dataframe()

    def _check_blastn(self) -> bool:
        """Check if blastn is available in PATH."""
        try:
            result = subprocess.run(
                ['blastn', '-version'], 
                capture_output=True, 
                text=True, 
                timeout=10
            )
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
            return False

    def _parse_blastn_output(
        self, 
        output: str, 
        min_identity: float = 80.0,
        min_coverage: float = 90.0
    ) -> pd.DataFrame:
        """
        Parse BLASTn tabular output to standardized DataFrame.

        BLAST outfmt 6 columns: qseqid, sseqid, pident, length, qstart, qend, sstart, send, sstrand, slen

        Args:
            output: Raw BLASTn output text
            min_identity: Minimum percent identity
            min_coverage: Minimum coverage percentage

        Returns:
            Standardized DataFrame
        """
        if not output or not output.strip():
            logger.debug("No housekeeping genes found")
            return self._empty_dataframe()
        
        try:
            # Parse tab-delimited output (includes slen for coverage calculation)
            columns = ['qseqid', 'sseqid', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'sstrand', 'slen']
            df = pd.read_csv(StringIO(output), sep='\t', names=columns, comment='#')
            
            if df.empty:
                return self._empty_dataframe()
            
            # Filter by identity
            df = df[df['pident'] >= min_identity]
            
            # Calculate true coverage percentage (alignment length / subject length)
            df['coverage_pct'] = (df['length'] / df['slen']) * 100
            df = df[df['coverage_pct'] >= min_coverage]
            
            if df.empty:
                return self._empty_dataframe()
            
            # Standardize columns
            df['Gene'] = df['sseqid']  # Subject sequence ID (gene name from reference)
            df['Contig'] = df['qseqid']  # Query sequence ID (contig name)
            
            # Use QUERY coordinates (qstart/qend) for contig positions
            # Ensure Start < End for Python slicing (regardless of strand)
            df['Start'] = df[['qstart', 'qend']].min(axis=1)
            df['End'] = df[['qstart', 'qend']].max(axis=1)
            
            # Strand standardization (BLAST returns 'plus'/'minus', we use +/-)
            df['Strand'] = df['sstrand'].apply(lambda x: '-' if 'minus' in str(x) else '+')
            
            df['Identity'] = df['pident']
            df['Source'] = 'BLAST'
            
            # Select standardized columns
            standard_cols = ['Gene', 'Contig', 'Start', 'End', 'Strand', 'Identity', 'Source']
            df = df[standard_cols]
            
            # Ensure correct data types
            df['Start'] = df['Start'].astype(int)
            df['End'] = df['End'].astype(int)
            df['Identity'] = df['Identity'].astype(float)
            
            logger.info(f"Found {len(df)} housekeeping genes after QC")
            return df
            
        except Exception as e:
            logger.error(f"Error parsing BLASTn output: {e}")
            return self._empty_dataframe()

    def find_all_genes(
        self, 
        fasta_file: Path, 
        reference_db: Optional[Path] = None
    ) -> pd.DataFrame:
        """
        Find both resistance and housekeeping genes in a single call.

        Args:
            fasta_file: Path to input FASTA file
            reference_db: Optional path to housekeeping gene reference database

        Returns:
            Combined DataFrame with all genes found
        """
        logger.info(f"Running hybrid gene detection on {fasta_file.name}")
        
        # Find resistance genes
        resistance_df = self.find_resistance_genes(fasta_file)
        
        # Find housekeeping genes if reference provided
        if reference_db and Path(reference_db).exists():
            housekeeping_df = self.find_housekeeping_genes(fasta_file, reference_db)
            
            # Combine results
            combined_df = pd.concat([resistance_df, housekeeping_df], ignore_index=True)
            logger.info(f"Total genes found: {len(combined_df)}")
            return combined_df
        else:
            logger.info(f"Total genes found: {len(resistance_df)} (resistance only)")
            return resistance_df

    def _empty_dataframe(self) -> pd.DataFrame:
        """Return empty DataFrame with standardized columns."""
        return pd.DataFrame(columns=['Gene', 'Contig', 'Start', 'End', 'Strand', 'Identity', 'Source'])
