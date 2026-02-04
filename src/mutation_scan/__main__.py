"""
Entry point for running MutationScan as a module.

Usage:
    python -m mutation_scan.main [OPTIONS]
    
This allows running the package from anywhere in the source tree.
"""

from mutation_scan.main import main

if __name__ == "__main__":
    main()
