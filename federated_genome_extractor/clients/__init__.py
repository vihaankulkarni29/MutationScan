"""
Clients package for federated genome extractor.

This package contains client modules for different bioinformatics databases:
- NCBIClient: NCBI Entrez/GenBank genome extraction
- BVBRCClient: BV-BRC (Bacterial and Viral Bioinformatics Resource Center)
- EnteroBaseClient: EnteroBase specialized database  
- PATRICClient: PATRIC pathogen database

Each client provides standardized methods for genome search, download,
and metadata extraction from their respective databases.
"""

# Make clients easily importable
try:
    from .ncbi_client import NCBIClient
    from .bvbrc_client import BVBRCClient
    
    __all__ = ['NCBIClient', 'BVBRCClient']
    
    # Optional clients that might have additional dependencies
    try:
        from .enterobase_client import EnteroBaseClient
        __all__.append('EnteroBaseClient')
    except ImportError:
        pass
        
    try:
        from .patric_client import PATRICClient  
        __all__.append('PATRICClient')
    except ImportError:
        pass
        
except ImportError as e:
    # Graceful degradation if core dependencies are missing
    __all__ = []