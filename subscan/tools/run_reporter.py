#!/usr/bin/env python3
"""
MutationScan Domino 7: Interactive HTML Report Generator

This module implements the final domino in the MutationScan pipeline, responsible
for generating comprehensive, interactive HTML reports from complete pipeline results.
It serves as the culmination of the AMR analysis workflow, creating professional
dashboards that enable researchers to explore and visualize mutation patterns.

The reporter aggregates data from all upstream domino tools and creates sophisticated
interactive visualizations including mutation frequency plots, co-occurrence heatmaps,
sequence alignment viewers, and comprehensive summary statistics. The output is a
self-contained HTML file optimized for scientific presentation and data sharing.

Key Features:
- Comprehensive data aggregation from all pipeline stages
- Interactive Plotly visualizations for mutation analysis
- MSA viewer integration for sequence alignment exploration  
- Professional dashboard design for scientific presentation
- Single-file HTML output for easy sharing and archival
- Responsive design optimized for different screen sizes

Usage:
    python run_reporter.py --manifest cooccurrence_manifest.json --output-dir ./reports

Integration:
    - Input: cooccurrence_manifest.json from Co-occurrence Analyzer (Domino 6)
    - Output: Interactive HTML dashboard + supporting visualization files
    - End Product: Complete AMR analysis report for researchers

Visualization Components:
    - Genome metadata summary tables
    - AMR gene distribution charts
    - Mutation frequency heatmaps
    - Co-occurrence pattern networks
    - Interactive sequence alignment viewers
    - Statistical summary panels

Dependencies:
    - pandas for data manipulation and analysis
    - plotly for interactive visualization generation
    - HTML templating for dashboard assembly
    - subscan.reporting module for specialized report functions

Author: MutationScan Development Team
Version: 1.0.0
"""

import argparse
import os
import sys
from typing import Dict, Any, List, Optional, Tuple
from pathlib import Path
import pandas as pd

# Add src directory to Python path for imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "src"))


def create_argument_parser() -> argparse.ArgumentParser:
    """
    Create and configure command-line argument parser for report generation.
    
    Sets up argument parsing for the reporter with validation for final pipeline
    manifest input, output directory specification, and report customization
    options for generating comprehensive HTML analysis dashboards.

    Returns:
        argparse.ArgumentParser: Configured parser for reporter arguments with:
            - manifest: Path to cooccurrence_manifest.json from final pipeline run
            - output_dir: Directory for HTML report and supporting files
            - report_name: Custom name for the output HTML file
            - verbose: Enable detailed progress and debugging output

    Example:
        >>> parser = create_argument_parser()
        >>> args = parser.parse_args(['--manifest', 'final_manifest.json', '--output-dir', 'reports/'])
    """
    parser = argparse.ArgumentParser(
        description="SubScan Domino 7: HTML Report Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_reporter.py --manifest cooccurrence_manifest.json --output-dir ./reports
  python run_reporter.py --manifest results/final_analysis.json --output-dir ./dashboard

Pipeline Flow:
  Complete Pipeline Results → cooccurrence_manifest.json → Domino 7 (Reporter) → Interactive HTML Dashboard

Features:
  - Comprehensive data aggregation from all pipeline stages
  - Interactive Plotly visualizations for mutation frequencies
  - MSA viewer integration for sequence alignment exploration
  - Professional dashboard design for scientific presentation
  - Single-file HTML output for easy sharing and archival

Input Requirements:
  - cooccurrence_manifest.json: Final manifest from complete pipeline run
  - All referenced result files must be accessible at their manifest paths
        """,
    )

    # Required arguments
    parser.add_argument(
        "--manifest",
        type=str,
        required=True,
        help="Path to the cooccurrence_manifest.json file from the complete pipeline run",
    )

    parser.add_argument(
        "--output-dir",
        type=str,
        required=True,
        help="Directory where the HTML report will be generated",
    )

    # Optional arguments
    parser.add_argument(
        "--report-name",
        type=str,
        default="subscan_final_report.html",
        help="Name of the output HTML file (default: subscan_final_report.html)",
    )

    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging for detailed progress tracking",
    )

    parser.add_argument(
        "--open-browser",
        action="store_true",
        help="Open the generated report in the default web browser after creation",
    )

    return parser


def validate_arguments(args: argparse.Namespace) -> None:
    """Validate command-line arguments using shared utilities"""
    # Import shared utilities
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), "src"))
    from subscan.utils import (
        validate_manifest_file,
        validate_output_directory,
        print_validation_success,
    )

    try:
        # Use shared validation functions
        validate_manifest_file(args.manifest)
        validate_output_directory(args.output_dir, create=True)

        if args.verbose:
            print_validation_success(
                args.manifest, args.output_dir, report_name=args.report_name
            )

    except (FileNotFoundError, OSError, PermissionError) as e:
        print(f"❌ Error: {e}")
        sys.exit(1)


def main() -> int:
    """
    Main entry point for the HTML Report Generator.
    
    This function orchestrates the complete report generation workflow:
    1. Parses and validates command-line arguments
    2. Initializes the ReportGenerator with the pipeline manifest
    3. Aggregates data from all pipeline stages
    4. Prepares interactive visualization components
    5. Generates the final HTML dashboard
    
    Returns:
        int: Exit code (0 for success, 1 for error)
        
    Raises:
        SystemExit: If argument validation fails or critical errors occur
    """

    print("SubScan Pipeline - Domino 7: The HTML Report Generator")
    print("=" * 60)

    # Parse and validate arguments
    parser = create_argument_parser()
    args = parser.parse_args()
    validate_arguments(args)

    print(f"📊 Input manifest: {args.manifest}")
    print(f"📁 Output directory: {args.output_dir}")
    print(f"📄 Report file: {args.report_name}")
    print()

    # TODO Phase 2: Core Data Aggregation Logic
    # Create ReportGenerator instance and aggregate data
    from subscan.reporting import ReportGenerator

    try:
        # Initialize the report generator
        report_generator = ReportGenerator(args.manifest)

        # Aggregate all pipeline data
        aggregated_data = report_generator.aggregate_data()

        # Prepare interactive visualization data
        plotly_data = report_generator.prepare_plotly_data()
        msa_data = report_generator.prepare_msa_data()

        print(f"✅ Data aggregation completed successfully")
        print(f"   - Total genomes: {aggregated_data['statistics']['total_genomes']}")
        print(
            f"   - Total mutations: {aggregated_data['statistics']['total_mutations']}"
        )
        print(
            f"   - Total co-occurrence pairs: {aggregated_data['statistics']['total_cooccurrence_pairs']}"
        )
        print(
            f"   - Pipeline completion: {aggregated_data['statistics']['data_completeness']:.1f}%"
        )

    except Exception as e:
        print(f"❌ Data aggregation failed: {e}")
        return 1

    # Phase 3: Interactive Visualization Methods (completed in Phase 2)
    print("✅ Phase 3: Interactive visualization data prepared")

    # Phase 4: Final Rendering with Jinja2 Integration
    try:
        # Determine output file path
        # Use provided report name within the specified output directory
        html_output_path = os.path.join(args.output_dir, args.report_name)

        # Generate the HTML report
        final_report_path = report_generator.generate_html_report(html_output_path)

        # Get report summary for CLI output
        report_summary = report_generator.generate_report_summary()

        print(f"\n🎉 HTML Report Generation Complete!")
        print(f"📄 Report saved to: {final_report_path}")
        print(f"� Pipeline completion: {report_summary['pipeline_completion']:.1f}%")
        print(f"🧬 Data processed:")
        print(f"   - Genomes: {report_summary['data_summary']['genomes_processed']:,}")
        print(
            f"   - Mutations: {report_summary['data_summary']['mutations_detected']:,}"
        )
        print(
            f"   - Co-occurrence pairs: {report_summary['data_summary']['cooccurrence_pairs']:,}"
        )
        print(
            f"   - Gene families: {report_summary['data_summary']['gene_families']:,}"
        )
        print(f"📈 Interactive components:")
        print(
            f"   - Charts: {report_summary['visualization_components']['interactive_charts']}"
        )
        print(
            f"   - Data tables: {report_summary['visualization_components']['data_tables']}"
        )
        print(
            f"   - MSA alignments: {report_summary['visualization_components']['msa_alignments']}"
        )

        if args.open_browser:
            try:
                import webbrowser

                webbrowser.open(f"file://{os.path.abspath(final_report_path)}")
                print(f"🌐 Opened report in default browser")
            except Exception as e:
                print(f"⚠️  Could not open browser: {e}")

        print(f"\n✅ SubScan HTML Report Generator completed successfully!")
        return 0

    except Exception as e:
        print(f"❌ HTML report generation failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
