from setuptools import setup, find_packages

setup(
    name="abricate-automator",
    version="1.0.0",
    description="Automated AMR gene annotation using ABRicate for the MutationScan pipeline",
    author="Vihaan Kulkarni",
    author_email="vihaankulkarni29@gmail.com",
    url="https://github.com/vihaankulkarni29/ABRicate-Automator",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "pandas>=1.3.0,<3.0",
        "biopython>=1.79,<2.0",
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