from setuptools import setup, find_packages

# Fallback setup.py for older pip versions
# Primary configuration is in pyproject.toml

setup(
    name="federated_genome_extractor",
    version="1.0.0",
    description="Federated genome extraction tool for downloading genomes from multiple bioinformatics databases",
    long_description=open("README.md", "r", encoding="utf-8").read() if __file__ else "",
    long_description_content_type="text/markdown",
    author="Vihaan Kulkarni",
    author_email="vihaankulkarni29@gmail.com",
    url="https://github.com/vihaankulkarni29/federated_genome_extractor",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.79,<2.0",
        "requests>=2.25.0,<3.0",
        "urllib3>=1.26.0,<3.0", 
        "tqdm>=4.62.0,<5.0",
        "pandas>=1.3.0,<3.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0", 
            "black>=23.0.0",
            "flake8>=6.0.0",
            "mypy>=1.0.0",
        ]
    },
    entry_points={
        "console_scripts": [
            "federated-genome-extractor=harvester:main",
        ],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research", 
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Operating System :: OS Independent",
    ],
    keywords="bioinformatics genomics ncbi genome-extraction antimicrobial-resistance federated-search",
    license="MIT",
)
