from setuptools import setup, find_packages

setup(
    name="subscan-analyzer",
    version="1.0.0",
    description="Mutation analysis and detection tool for the MutationScan pipeline",
    author="Vihaan Kulkarni",
    author_email="vihaankulkarni29@gmail.com",
    url="https://github.com/vihaankulkarni29/SubScan",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "biopython>=1.79,<2.0",
        "pandas>=1.3.0,<3.0",
    ],
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    license="MIT",
)