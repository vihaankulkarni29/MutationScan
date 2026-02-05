"""
Setup configuration for MutationScan package.

Defines package metadata, dependencies, and installation configuration.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
README_PATH = Path(__file__).parent / "README.md"
long_description = README_PATH.read_text(encoding="utf-8") if README_PATH.exists() else ""

setup(
    name="mutation_scan",
    version="1.0.0",
    author="Vihaan Kulkarni",
    author_email="contact@example.com",
    description="Democratized Bioinformatics Pipeline for AMR Detection",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vihaankulkarni29/MutationScan",
    project_urls={
        "Bug Tracker": "https://github.com/vihaankulkarni29/MutationScan/issues",
        "Documentation": "https://github.com/vihaankulkarni29/MutationScan/docs",
        "Source Code": "https://github.com/vihaankulkarni29/MutationScan",
    },
    license="MIT",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.9,<4",
    install_requires=[
        "biopython>=1.79",
        "pyyaml>=5.4",
        "requests>=2.26.0",
        "pandas>=1.3.0",
        "numpy>=1.21.0",
        "matplotlib>=3.4.0",
        "scipy>=1.7.0",
        "scikit-learn>=1.0.0",
        "joblib>=1.1.0",
        "ncbi-datasets-pylib>=16.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.2.5",
            "pytest-cov>=2.12.1",
            "black>=21.7b0",
            "flake8>=3.9.2",
            "mypy>=0.910",
            "isort>=5.9.3",
        ],
        "docs": [
            "sphinx>=4.0.0",
            "sphinx-rtd-theme>=0.5.0",
        ],
        "pymol": [
            "pymol-open-source>=2.4.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "mutation-scan=mutation_scan.main:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    zip_safe=False,
    keywords=[
        "bioinformatics",
        "antimicrobial-resistance",
        "amr",
        "genomics",
        "sequence-analysis",
        "mutation-detection",
    ],
)
