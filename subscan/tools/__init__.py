"""
SubScan Pipeline Tools Package

Command-line interface tools for the SubScan bioinformatics pipeline.
This package contains all seven "domino" tools that form the complete
mutation analysis workflow.

Tools:
    run_harvester: NCBI genome extraction (Domino 1)
    run_annotator: AMR gene annotation with ABRicate (Domino 2)  
    run_extractor: Protein sequence extraction (Domino 3)
    run_aligner: Wild-type sequence alignment (Domino 4)
    run_analyzer: Mutation detection and analysis (Domino 5)
    run_cooccurrence_analyzer: Co-occurrence pattern analysis (Domino 6)
    run_reporter: Interactive HTML report generation (Domino 7)
    run_pipeline: Master orchestrator for complete pipeline

Author: MutationScan Development Team
"""
