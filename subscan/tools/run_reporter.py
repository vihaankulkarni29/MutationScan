#!/usr/bin/env python3
"""
SubScan Pipeline - Domino 7: The HTML Report Generator

This script serves as the final domino in the SubScan bioinformatics pipeline,
generating comprehensive, interactive HTML reports from the complete pipeline results.

Purpose:
- Reads cooccurrence_manifest.json from the complete pipeline run
- Aggregates all analytical results (metadata, mutations, co-occurrence data)
- Generates a professional, interactive HTML dashboard for researchers
- Provides visual exploration from high-level summaries to detailed sequence alignments

Author: SubScan Pipeline Development Team
"""

import argparse
import os
import sys
import pandas as pd

# Add src directory to Python path for imports
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(__file__)), 'src'))


def create_argument_parser():
    """Create and configure the command-line argument parser"""
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

    return parser.parse_args()


def validate_arguments(args):
    """Validate command-line arguments"""
    # Check manifest file exists
    if not os.path.isfile(args.manifest):
        print(f"❌ Error: Manifest file not found: {args.manifest}")
        sys.exit(1)

    # Validate manifest is JSON
    if not args.manifest.lower().endswith(".json"):
        print(f"⚠️  Warning: Manifest file should be JSON: {args.manifest}")

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    # Validate output directory is writable
    if not os.access(args.output_dir, os.W_OK):
        print(f"❌ Error: Output directory is not writable: {args.output_dir}")
        sys.exit(1)

    if args.verbose:
        print(f"✅ Manifest file: {args.manifest}")
        print(f"✅ Output directory: {args.output_dir}")
        print(f"✅ Report name: {args.report_name}")


def main():
    """Main entry point for the HTML Report Generator"""

    print("SubScan Pipeline - Domino 7: The HTML Report Generator")
    print("=" * 60)

    # Parse and validate arguments
    args = create_argument_parser()
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
        if args.output:
            html_output_path = args.output
        else:
            # Generate default output filename
            timestamp = pd.Timestamp.now().strftime("%Y%m%d_%H%M%S")
            html_output_path = os.path.join(
                args.output_dir, f"subscan_report_{timestamp}.html"
            )

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
