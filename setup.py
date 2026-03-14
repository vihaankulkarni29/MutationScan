from setuptools import setup, find_packages

setup(
    name="mutation_scan",
    version="2.1.0",
    description="Deterministic Genotype-to-Phenotype-to-Biophysics Engine",
    author="Vihaan Kulkarni",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "pandas",
        "biopython",
        "networkx",
        "matplotlib",
        "pyyaml",
        "requests"
    ],
    python_requires=">=3.10",
)
