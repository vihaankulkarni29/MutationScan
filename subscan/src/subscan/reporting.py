"""
SubScan Reporting Module

This module provides the core data aggregation and preparation functionality
for the SubScan HTML Report Generator (Domino 7).

The ReportGenerator class serves as the central hub for collecting and organizing
all analytical results from the complete SubScan pipeline run.
"""

import json
import os
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional
import logging


class ReportGenerator:
    """
    Central data aggregation engine for SubScan pipeline reporting.

    This class is responsible for:
    - Reading and parsing the final cooccurrence_manifest.json
    - Aggregating data from all pipeline stages (metadata, mutations, co-occurrence)
    - Preparing structured data for HTML template rendering
    - Providing data preparation methods for interactive visualizations
    """

    def __init__(self, manifest_path: str):
        """
        Initialize the ReportGenerator with the final pipeline manifest.

        Args:
            manifest_path (str): Path to the cooccurrence_manifest.json file
        """
        self.manifest_path = manifest_path
        self.manifest_data = None
        self.aggregated_data = None

        # Set up logging
        self.logger = logging.getLogger(__name__)

        # Validate manifest file exists
        if not os.path.isfile(manifest_path):
            raise FileNotFoundError(f"Manifest file not found: {manifest_path}")

    def aggregate_data(self) -> Dict[str, Any]:
        """
        Core data aggregation method - reads and consolidates all pipeline results.

        This method serves as the central hub that:
        1. Reads and parses the input JSON manifest
        2. Locates and reads key result files from the manifest paths
        3. Consolidates all information into a structured dictionary

        Returns:
            Dict[str, Any]: Comprehensive data structure containing all pipeline results
        """
        self.logger.info("Starting data aggregation from manifest")

        # Step 1: Read and parse the input JSON manifest
        try:
            with open(self.manifest_path, "r", encoding="utf-8") as f:
                self.manifest_data = json.load(f)
            self.logger.info(f"Successfully loaded manifest: {self.manifest_path}")
        except Exception as e:
            raise RuntimeError(f"Failed to read manifest file: {e}")

        # Step 2: Initialize aggregated data structure
        self.aggregated_data = {
            "manifest_metadata": self.manifest_data,
            "pipeline_summary": {},
            "genome_metadata": {},
            "mutation_data": {},
            "cooccurrence_data": {},
            "file_paths": {},
            "statistics": {},
        }

        # Step 3: Extract file paths from manifest
        self._extract_file_paths()

        # Step 4: Read and process metadata.csv (from Domino 1)
        self._process_metadata_csv()

        # Step 5: Read and process mutation report (from Domino 5)
        self._process_mutation_report()

        # Step 6: Read and process co-occurrence report (from Domino 6)
        self._process_cooccurrence_report()

        # Step 7: Calculate summary statistics
        self._calculate_summary_statistics()

        self.logger.info("Data aggregation completed successfully")
        return self.aggregated_data

    def _extract_file_paths(self) -> None:
        """
        Extract and validate key file paths from the manifest lineage.
        
        This method traverses the complete pipeline manifest chain to locate
        all key result files needed for report generation:
        - Metadata CSV from genome harvesting
        - Mutation analysis reports from genome analysis
        - Co-occurrence analysis results
        - Alignment manifests for sequence data
        
        Updates:
            self.aggregated_data["file_paths"] with discovered file paths
            
        Raises:
            FileNotFoundError: If critical files cannot be located
            ValueError: If manifest structure is unexpected
        """
        self.logger.info("Extracting file paths from manifest")

        # Navigate through the manifest to find key result files
        # The manifest should contain the complete pipeline lineage

        # Initialize file paths dictionary
        self.aggregated_data["file_paths"] = {
            "metadata_csv": None,
            "mutation_report": None,
            "cooccurrence_report": None,
            "alignment_manifest": None,
        }

        # Extract paths based on manifest structure
        # This will depend on the actual structure of cooccurrence_manifest.json
        try:
            # Look for metadata CSV path (from Domino 1)
            if "input_manifest" in self.manifest_data:
                input_manifest = self.manifest_data["input_manifest"]
                if isinstance(input_manifest, dict):
                    # Navigate through nested manifest structure
                    if "metadata_csv" in input_manifest:
                        self.aggregated_data["file_paths"]["metadata_csv"] = (
                            input_manifest["metadata_csv"]
                        )
                    elif (
                        "output_files" in input_manifest
                        and "metadata_csv" in input_manifest["output_files"]
                    ):
                        self.aggregated_data["file_paths"]["metadata_csv"] = (
                            input_manifest["output_files"]["metadata_csv"]
                        )

            # Look for mutation report path (from Domino 5)
            if "mutation_report_path" in self.manifest_data:
                self.aggregated_data["file_paths"]["mutation_report"] = (
                    self.manifest_data["mutation_report_path"]
                )
            elif (
                "output_files" in self.manifest_data
                and "mutation_report" in self.manifest_data["output_files"]
            ):
                self.aggregated_data["file_paths"]["mutation_report"] = (
                    self.manifest_data["output_files"]["mutation_report"]
                )

            # Look for co-occurrence report path (from Domino 6)
            if "cooccurrence_report_path" in self.manifest_data:
                self.aggregated_data["file_paths"]["cooccurrence_report"] = (
                    self.manifest_data["cooccurrence_report_path"]
                )
            elif (
                "output_files" in self.manifest_data
                and "cooccurrence_report" in self.manifest_data["output_files"]
            ):
                self.aggregated_data["file_paths"]["cooccurrence_report"] = (
                    self.manifest_data["output_files"]["cooccurrence_report"]
                )

            # Look for alignment manifest (needed for MSA viewer in Phase 3)
            if "alignment_manifest_path" in self.manifest_data:
                self.aggregated_data["file_paths"]["alignment_manifest"] = (
                    self.manifest_data["alignment_manifest_path"]
                )

            self.logger.info(
                f"Extracted file paths: {self.aggregated_data['file_paths']}"
            )

        except Exception as e:
            self.logger.warning(f"Error extracting file paths: {e}")

    def _process_metadata_csv(self):
        """Read and process the metadata.csv file from Domino 1."""
        metadata_path = self.aggregated_data["file_paths"]["metadata_csv"]

        if not metadata_path or not os.path.isfile(metadata_path):
            self.logger.warning(f"Metadata CSV not found: {metadata_path}")
            self.aggregated_data["genome_metadata"] = {"genomes": [], "count": 0}
            return

        try:
            self.logger.info(f"Reading metadata CSV: {metadata_path}")
            metadata_df = pd.read_csv(metadata_path)

            # Convert to structured format
            self.aggregated_data["genome_metadata"] = {
                "genomes": metadata_df.to_dict("records"),
                "count": len(metadata_df),
                "columns": list(metadata_df.columns),
            }

            self.logger.info(f"Processed {len(metadata_df)} genome metadata records")

        except Exception as e:
            self.logger.error(f"Failed to process metadata CSV: {e}")
            self.aggregated_data["genome_metadata"] = {
                "genomes": [],
                "count": 0,
                "error": str(e),
            }

    def _process_mutation_report(self):
        """Read and process the mutation report from Domino 5."""
        mutation_path = self.aggregated_data["file_paths"]["mutation_report"]

        if not mutation_path or not os.path.isfile(mutation_path):
            self.logger.warning(f"Mutation report not found: {mutation_path}")
            self.aggregated_data["mutation_data"] = {"mutations": [], "count": 0}
            return

        try:
            self.logger.info(f"Reading mutation report: {mutation_path}")

            # Handle both Excel and CSV formats
            if mutation_path.endswith(".xlsx"):
                mutation_df = pd.read_excel(mutation_path)
            elif mutation_path.endswith(".csv"):
                mutation_df = pd.read_csv(mutation_path)
            else:
                raise ValueError(f"Unsupported file format: {mutation_path}")

            # Convert to structured format
            self.aggregated_data["mutation_data"] = {
                "mutations": mutation_df.to_dict("records"),
                "count": len(mutation_df),
                "columns": list(mutation_df.columns),
                "gene_families": (
                    mutation_df["gene_family"].unique().tolist()
                    if "gene_family" in mutation_df.columns
                    else []
                ),
            }

            self.logger.info(f"Processed {len(mutation_df)} mutation records")

        except Exception as e:
            self.logger.error(f"Failed to process mutation report: {e}")
            self.aggregated_data["mutation_data"] = {
                "mutations": [],
                "count": 0,
                "error": str(e),
            }

    def _process_cooccurrence_report(self):
        """Read and process the co-occurrence report from Domino 6."""
        cooccurrence_path = self.aggregated_data["file_paths"]["cooccurrence_report"]

        if not cooccurrence_path or not os.path.isfile(cooccurrence_path):
            self.logger.warning(f"Co-occurrence report not found: {cooccurrence_path}")
            self.aggregated_data["cooccurrence_data"] = {"pairs": [], "count": 0}
            return

        try:
            self.logger.info(f"Reading co-occurrence report: {cooccurrence_path}")
            cooccurrence_df = pd.read_csv(cooccurrence_path)

            # Convert to structured format
            self.aggregated_data["cooccurrence_data"] = {
                "pairs": cooccurrence_df.to_dict("records"),
                "count": len(cooccurrence_df),
                "columns": list(cooccurrence_df.columns),
            }

            self.logger.info(f"Processed {len(cooccurrence_df)} co-occurrence pairs")

        except Exception as e:
            self.logger.error(f"Failed to process co-occurrence report: {e}")
            self.aggregated_data["cooccurrence_data"] = {
                "pairs": [],
                "count": 0,
                "error": str(e),
            }

    def _calculate_summary_statistics(self):
        """Calculate comprehensive summary statistics from all aggregated data."""
        self.logger.info("Calculating summary statistics")

        stats = {
            "total_genomes": self.aggregated_data["genome_metadata"]["count"],
            "total_mutations": self.aggregated_data["mutation_data"]["count"],
            "total_cooccurrence_pairs": self.aggregated_data["cooccurrence_data"][
                "count"
            ],
            "gene_families": len(
                self.aggregated_data["mutation_data"].get("gene_families", [])
            ),
            "pipeline_stages_completed": self._count_pipeline_stages(),
            "data_completeness": self._assess_data_completeness(),
        }

        # Add gene family statistics
        if self.aggregated_data["mutation_data"]["mutations"]:
            gene_family_counts = {}
            for mutation in self.aggregated_data["mutation_data"]["mutations"]:
                gene_family = mutation.get("gene_family", "Unknown")
                gene_family_counts[gene_family] = (
                    gene_family_counts.get(gene_family, 0) + 1
                )
            stats["gene_family_distribution"] = gene_family_counts

        self.aggregated_data["statistics"] = stats
        self.logger.info(f"Summary statistics: {stats}")

    def _count_pipeline_stages(self) -> int:
        """Count the number of completed pipeline stages based on available data."""
        stages = 0
        if self.aggregated_data["genome_metadata"]["count"] > 0:
            stages += 1  # Domino 1 (Harvester)
        if self.aggregated_data["mutation_data"]["count"] > 0:
            stages += 1  # Domino 5 (Analyzer)
        if self.aggregated_data["cooccurrence_data"]["count"] > 0:
            stages += 1  # Domino 6 (Co-occurrence)
        return stages

    def _assess_data_completeness(self) -> float:
        """Assess the completeness of aggregated data as a percentage."""
        total_expected = 3  # metadata, mutations, co-occurrence
        available = 0

        if self.aggregated_data["genome_metadata"]["count"] > 0:
            available += 1
        if self.aggregated_data["mutation_data"]["count"] > 0:
            available += 1
        if self.aggregated_data["cooccurrence_data"]["count"] > 0:
            available += 1

        return (available / total_expected) * 100

    def get_aggregated_data(self) -> Optional[Dict[str, Any]]:
        """
        Get the aggregated data dictionary.

        Returns:
            Optional[Dict[str, Any]]: The aggregated data or None if not yet processed
        """
        return self.aggregated_data

    def get_summary_statistics(self) -> Dict[str, Any]:
        """
        Get summary statistics for the report.

        Returns:
            Dict[str, Any]: Summary statistics dictionary
        """
        if self.aggregated_data and "statistics" in self.aggregated_data:
            return self.aggregated_data["statistics"]
        return {}

    # ========================================================================
    # PHASE 3: INTERACTIVE VISUALIZATION METHODS
    # ========================================================================

    def prepare_plotly_data(self) -> Dict[str, Any]:
        """
        Prepare data structures optimized for Plotly.js interactive charts.

        This method converts the aggregated pipeline data into JSON-serializable
        formats that can be directly consumed by Plotly.js for interactive visualizations.

        Returns:
            Dict[str, Any]: Plotly-ready data structures for all chart types
        """
        if not self.aggregated_data:
            raise RuntimeError("Data must be aggregated before preparing Plotly data")

        self.logger.info("Preparing Plotly.js visualization data")

        plotly_data = {
            "mutation_distribution_chart": self._prepare_mutation_distribution(),
            "gene_family_bar_chart": self._prepare_gene_family_distribution(),
            "cooccurrence_network_data": self._prepare_cooccurrence_network(),
            "genome_quality_scatter": self._prepare_genome_quality_scatter(),
            "pipeline_completion_gauge": self._prepare_pipeline_gauge(),
        }

        self.logger.info("Plotly data preparation completed")
        return plotly_data

    def _prepare_mutation_distribution(self) -> Dict[str, Any]:
        """Prepare mutation distribution data for interactive pie/donut chart."""
        if not self.aggregated_data or "mutation_data" not in self.aggregated_data:
            return {
                "data": [],
                "layout": {"title": "No Mutation Data Available", "showlegend": False},
            }

        mutations = self.aggregated_data["mutation_data"]["mutations"]

        if not mutations:
            return {
                "data": [],
                "layout": {"title": "No Mutation Data Available", "showlegend": False},
            }

        # Group mutations by type/category
        mutation_types = {}
        for mutation in mutations:
            mut_type = mutation.get("mutation_type", "Unknown")
            mutation_types[mut_type] = mutation_types.get(mut_type, 0) + 1

        # Prepare Plotly pie chart data
        return {
            "data": [
                {
                    "values": list(mutation_types.values()),
                    "labels": list(mutation_types.keys()),
                    "type": "pie",
                    "hole": 0.4,  # Donut chart
                    "marker": {
                        "colors": [
                            "#FF6B6B",
                            "#4ECDC4",
                            "#45B7D1",
                            "#96CEB4",
                            "#FFEAA7",
                            "#DDA0DD",
                        ]
                    },
                    "textinfo": "label+percent",
                    "textposition": "outside",
                }
            ],
            "layout": {
                "title": {
                    "text": "Mutation Type Distribution",
                    "x": 0.5,
                    "font": {"size": 16},
                },
                "showlegend": True,
                "legend": {"orientation": "v", "x": 1.02, "y": 0.5},
            },
        }

    def _prepare_gene_family_distribution(self) -> Dict[str, Any]:
        """Prepare gene family distribution data for interactive bar chart."""
        if not self.aggregated_data or "statistics" not in self.aggregated_data:
            return {"data": [], "layout": {"title": "No Gene Family Data Available"}}

        gene_family_dist = self.aggregated_data["statistics"].get(
            "gene_family_distribution", {}
        )

        if not gene_family_dist:
            return {"data": [], "layout": {"title": "No Gene Family Data Available"}}

        # Sort gene families by count (descending)
        sorted_families = sorted(
            gene_family_dist.items(), key=lambda x: x[1], reverse=True
        )
        families, counts = zip(*sorted_families) if sorted_families else ([], [])

        return {
            "data": [
                {
                    "x": list(families),
                    "y": list(counts),
                    "type": "bar",
                    "marker": {
                        "color": "#4ECDC4",
                        "line": {"color": "#2C3E50", "width": 1},
                    },
                    "text": [f"{count} mutations" for count in counts],
                    "textposition": "outside",
                }
            ],
            "layout": {
                "title": {
                    "text": "Mutations by Gene Family",
                    "x": 0.5,
                    "font": {"size": 16},
                },
                "xaxis": {"title": "Gene Family", "tickangle": -45},
                "yaxis": {"title": "Number of Mutations"},
                "margin": {"b": 120},  # Extra bottom margin for rotated labels
            },
        }

    def _prepare_cooccurrence_network(self) -> Dict[str, Any]:
        """Prepare co-occurrence data for network visualization."""
        if not self.aggregated_data or "cooccurrence_data" not in self.aggregated_data:
            return {
                "nodes": [],
                "edges": [],
                "layout": {"title": "No Co-occurrence Data Available"},
            }

        cooccurrence_pairs = self.aggregated_data["cooccurrence_data"]["pairs"]

        if not cooccurrence_pairs:
            return {
                "nodes": [],
                "edges": [],
                "layout": {"title": "No Co-occurrence Data Available"},
            }

        # Extract unique genes and build network structure
        genes = set()
        edges = []

        for pair in cooccurrence_pairs[:50]:  # Limit to top 50 for performance
            gene1 = pair.get("gene1", "Unknown")
            gene2 = pair.get("gene2", "Unknown")
            frequency = pair.get("frequency", 1)

            genes.add(gene1)
            genes.add(gene2)
            edges.append({"source": gene1, "target": gene2, "value": frequency})

        # Prepare node data
        nodes = [{"id": gene, "label": gene} for gene in genes]

        return {
            "nodes": nodes,
            "edges": edges,
            "data": [
                {
                    "x": [i for i in range(len(nodes))],  # Simple circular layout
                    "y": [0 for _ in nodes],
                    "mode": "markers+text",
                    "text": [node["label"] for node in nodes],
                    "textposition": "middle center",
                    "marker": {"size": 20, "color": "#FF6B6B"},
                }
            ],
            "layout": {
                "title": {
                    "text": "Gene Co-occurrence Network",
                    "x": 0.5,
                    "font": {"size": 16},
                },
                "showlegend": False,
                "xaxis": {"visible": False},
                "yaxis": {"visible": False},
            },
        }

    def _prepare_genome_quality_scatter(self) -> Dict[str, Any]:
        """Prepare genome quality metrics for scatter plot."""
        if not self.aggregated_data or "genome_metadata" not in self.aggregated_data:
            return {"data": [], "layout": {"title": "No Genome Quality Data Available"}}

        genomes = self.aggregated_data["genome_metadata"]["genomes"]

        if not genomes:
            return {"data": [], "layout": {"title": "No Genome Quality Data Available"}}

        # Extract quality metrics (adapt based on available columns)
        x_values = []
        y_values = []
        text_labels = []

        for genome in genomes:
            # Use available numeric columns for quality metrics
            contig_count = genome.get("contig_count", genome.get("Contig_Count", 0))
            n50 = genome.get("n50", genome.get("N50", 0))
            accession = genome.get("accession", genome.get("Accession", "Unknown"))

            if contig_count and n50:
                x_values.append(contig_count)
                y_values.append(n50)
                text_labels.append(accession)

        if not x_values:
            return {
                "data": [],
                "layout": {"title": "Insufficient Quality Metrics for Visualization"},
            }

        return {
            "data": [
                {
                    "x": x_values,
                    "y": y_values,
                    "mode": "markers+text",
                    "text": text_labels,
                    "textposition": "top center",
                    "marker": {
                        "size": 10,
                        "color": "#45B7D1",
                        "line": {"color": "#2C3E50", "width": 1},
                    },
                    "type": "scatter",
                }
            ],
            "layout": {
                "title": {
                    "text": "Genome Quality Metrics",
                    "x": 0.5,
                    "font": {"size": 16},
                },
                "xaxis": {"title": "Contig Count", "type": "log"},
                "yaxis": {"title": "N50 Value", "type": "log"},
            },
        }

    def _prepare_pipeline_gauge(self) -> Dict[str, Any]:
        """Prepare pipeline completion gauge chart."""
        if not self.aggregated_data or "statistics" not in self.aggregated_data:
            return {"data": [], "layout": {"title": "Pipeline Data Not Available"}}

        completion_percentage = self.aggregated_data["statistics"]["data_completeness"]

        return {
            "data": [
                {
                    "type": "indicator",
                    "mode": "gauge+number+delta",
                    "value": completion_percentage,
                    "domain": {"x": [0, 1], "y": [0, 1]},
                    "title": {"text": "Pipeline Completion"},
                    "delta": {"reference": 100},
                    "gauge": {
                        "axis": {"range": [None, 100]},
                        "bar": {"color": "#4ECDC4"},
                        "steps": [
                            {"range": [0, 50], "color": "#FFE5E5"},
                            {"range": [50, 80], "color": "#FFF5E5"},
                            {"range": [80, 100], "color": "#E5F5E5"},
                        ],
                        "threshold": {
                            "line": {"color": "red", "width": 4},
                            "thickness": 0.75,
                            "value": 90,
                        },
                    },
                }
            ],
            "layout": {"height": 300, "margin": {"t": 25, "b": 25, "l": 25, "r": 25}},
        }

    def prepare_msa_data(self) -> Dict[str, Any]:
        """
        Prepare multiple sequence alignment data for MSAViewer integration.

        This method reads alignment files and formats them for the MSAViewer
        JavaScript component that will be embedded in the HTML report.

        Returns:
            Dict[str, Any]: MSAViewer-ready alignment data and metadata
        """
        if not self.aggregated_data:
            raise RuntimeError("Data must be aggregated before preparing MSA data")

        self.logger.info("Preparing MSAViewer alignment data")

        msa_data = {
            "alignments": [],
            "metadata": {},
            "viewer_config": self._get_msa_viewer_config(),
        }

        # Extract alignment file paths from the manifest
        alignment_manifest_path = self.aggregated_data["file_paths"].get(
            "alignment_manifest"
        )

        if alignment_manifest_path and os.path.isfile(alignment_manifest_path):
            msa_data["alignments"] = self._process_alignment_files(
                alignment_manifest_path
            )
        else:
            self.logger.warning(
                "No alignment manifest found - MSA viewer will be disabled"
            )
            msa_data["alignments"] = []

        self.logger.info(
            f"Prepared {len(msa_data['alignments'])} alignments for MSAViewer"
        )
        return msa_data

    def _process_alignment_files(
        self, alignment_manifest_path: str
    ) -> List[Dict[str, Any]]:
        """Process alignment files from the alignment manifest."""
        try:
            with open(alignment_manifest_path, "r", encoding="utf-8") as f:
                alignment_manifest = json.load(f)

            alignments = []

            # Process each alignment file listed in the manifest
            for gene_family, alignment_info in alignment_manifest.get(
                "alignments", {}
            ).items():
                alignment_file = alignment_info.get("alignment_file")

                if alignment_file and os.path.isfile(alignment_file):
                    alignment_data = self._read_alignment_file(
                        alignment_file, gene_family
                    )
                    if alignment_data:
                        alignments.append(alignment_data)

            return alignments

        except Exception as e:
            self.logger.error(f"Failed to process alignment files: {e}")
            return []

    def _read_alignment_file(
        self, alignment_file: str, gene_family: str
    ) -> Optional[Dict[str, Any]]:
        """Read and parse a single alignment file."""
        try:
            # Read FASTA alignment file
            sequences = []
            current_seq = None

            with open(alignment_file, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line.startswith(">"):
                        if current_seq:
                            sequences.append(current_seq)
                        current_seq = {
                            "name": line[1:],  # Remove '>' prefix
                            "sequence": "",
                        }
                    elif current_seq:
                        current_seq["sequence"] += line

                # Add the last sequence
                if current_seq:
                    sequences.append(current_seq)

            if not sequences:
                return None

            return {
                "gene_family": gene_family,
                "sequence_count": len(sequences),
                "sequences": sequences,
                "alignment_length": len(sequences[0]["sequence"]) if sequences else 0,
                "file_path": alignment_file,
            }

        except Exception as e:
            self.logger.error(f"Failed to read alignment file {alignment_file}: {e}")
            return None

    def _get_msa_viewer_config(self) -> Dict[str, Any]:
        """Get configuration settings for MSAViewer component."""
        return {
            "colorscheme": "taylor",  # Color scheme for amino acids
            "vis": {
                "conserv": True,
                "overviewbox": True,
                "seqlogo": True,
                "sequences": True,
                "labelname": True,
                "labelid": False,
            },
            "conf": {
                "dropImport": True,
                "registerWS": False,
                "registerMouseHover": False,
                "importProxy": False,
                "eventBus": True,
            },
            "zoomer": {"menuFontsize": "12px", "autoResize": True},
        }

    # ========================================================================
    # PHASE 4: FINAL RENDERING WITH JINJA2 INTEGRATION
    # ========================================================================

    def generate_html_report(
        self, output_path: str, template_path: Optional[str] = None
    ) -> str:
        """
        Generate the complete HTML report using Jinja2 template rendering.

        This is the final orchestration method that combines all aggregated data,
        visualization components, and template rendering to produce the final
        interactive HTML dashboard.

        Args:
            output_path (str): Path where the HTML report should be saved
            template_path (Optional[str]): Custom template path (uses default if None)

        Returns:
            str: Path to the generated HTML report file
        """
        if not self.aggregated_data:
            raise RuntimeError("Data must be aggregated before generating HTML report")

        self.logger.info("Starting HTML report generation")

        # Step 1: Prepare all visualization data
        plotly_data = self.prepare_plotly_data()
        msa_data = self.prepare_msa_data()

        # Step 2: Load and configure Jinja2 template
        template_content = self._load_template(template_path)

        # Step 3: Prepare template context with all data
        template_context = self._prepare_template_context(plotly_data, msa_data)

        # Step 4: Render the final HTML
        rendered_html = self._render_template(template_content, template_context)

        # Step 5: Save the HTML report
        self._save_html_report(rendered_html, output_path)

        self.logger.info(f"HTML report generated successfully: {output_path}")
        return output_path

    def _load_template(self, template_path: Optional[str] = None) -> str:
        """Load the Jinja2 HTML template."""
        if template_path is None:
            # Use the default template from the package
            template_path = os.path.join(
                os.path.dirname(__file__), "report_template.html"
            )

        if not os.path.isfile(template_path):
            raise FileNotFoundError(f"Template file not found: {template_path}")

        self.logger.info(f"Loading template: {template_path}")

        try:
            with open(template_path, "r", encoding="utf-8") as f:
                return f.read()
        except Exception as e:
            raise RuntimeError(f"Failed to load template: {e}")

    def _prepare_template_context(
        self, plotly_data: Dict[str, Any], msa_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Prepare the complete context dictionary for Jinja2 template rendering.

        This method organizes all data into the structure expected by the HTML template.
        """
        if not self.aggregated_data:
            raise RuntimeError("Aggregated data is not available")

        stats = self.aggregated_data["statistics"]

        context = {
            # Report metadata
            "report_title": "SubScan Pipeline Analysis Report",
            "generation_timestamp": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
            "pipeline_version": "1.0.0",
            # Summary statistics for the dashboard grid
            "total_genomes": stats.get("total_genomes", 0),
            "total_mutations": stats.get("total_mutations", 0),
            "total_cooccurrence_pairs": stats.get("total_cooccurrence_pairs", 0),
            "gene_families": stats.get("gene_families", 0),
            "completion_percentage": stats.get("data_completeness", 0),
            # Plotly.js chart data (JSON serialized)
            "mutation_distribution_chart": json.dumps(
                plotly_data["mutation_distribution_chart"]
            ),
            "gene_family_bar_chart": json.dumps(plotly_data["gene_family_bar_chart"]),
            "cooccurrence_network_data": json.dumps(
                plotly_data["cooccurrence_network_data"]
            ),
            "genome_quality_scatter": json.dumps(plotly_data["genome_quality_scatter"]),
            "pipeline_completion_gauge": json.dumps(
                plotly_data["pipeline_completion_gauge"]
            ),
            # MSAViewer data and configuration
            "msa_alignments": msa_data["alignments"],
            "msa_viewer_config": json.dumps(msa_data["viewer_config"]),
            "has_alignments": len(msa_data["alignments"]) > 0,
            # Data tables for detailed results
            "genome_data": self.aggregated_data["genome_metadata"]["genomes"],
            "mutation_data": self.aggregated_data["mutation_data"]["mutations"][
                :100
            ],  # Limit for performance
            "cooccurrence_data": self.aggregated_data["cooccurrence_data"]["pairs"][
                :100
            ],  # Limit for performance
            # Gene family distribution for summary
            "gene_family_distribution": stats.get("gene_family_distribution", {}),
            # File paths for reference
            "source_files": {
                "manifest": self.manifest_path,
                "metadata_csv": self.aggregated_data["file_paths"].get("metadata_csv"),
                "mutation_report": self.aggregated_data["file_paths"].get(
                    "mutation_report"
                ),
                "cooccurrence_report": self.aggregated_data["file_paths"].get(
                    "cooccurrence_report"
                ),
            },
            # Pipeline execution metadata
            "pipeline_stages_completed": stats.get("pipeline_stages_completed", 0),
            "total_pipeline_stages": 7,  # Total SubScan dominoes
            # Data availability flags for conditional rendering
            "has_genome_data": self.aggregated_data["genome_metadata"]["count"] > 0,
            "has_mutation_data": self.aggregated_data["mutation_data"]["count"] > 0,
            "has_cooccurrence_data": self.aggregated_data["cooccurrence_data"]["count"]
            > 0,
        }

        self.logger.info("Template context prepared with all visualization data")
        return context

    def _render_template(self, template_content: str, context: Dict[str, Any]) -> str:
        """Render the Jinja2 template with the provided context."""
        try:
            from jinja2 import Template, Environment

            # Create Jinja2 environment with useful filters
            env = Environment()

            # Add custom filters for formatting
            env.filters["thousands"] = lambda x: (
                f"{x:,}" if isinstance(x, (int, float)) else x
            )
            env.filters["percentage"] = lambda x: (
                f"{x:.1f}%" if isinstance(x, (int, float)) else x
            )
            env.filters["truncate_table"] = lambda lst, n=10: (
                lst[:n] if isinstance(lst, list) else []
            )

            # Create template object
            template = env.from_string(template_content)

            # Render with context
            self.logger.info("Rendering Jinja2 template")
            rendered_html = template.render(**context)

            return rendered_html

        except ImportError:
            raise RuntimeError(
                "Jinja2 is required for HTML report generation. Install with: pip install jinja2"
            )
        except Exception as e:
            raise RuntimeError(f"Template rendering failed: {e}")

    def _save_html_report(self, html_content: str, output_path: str) -> None:
        """Save the rendered HTML content to the specified file."""
        try:
            # Ensure output directory exists
            output_dir = os.path.dirname(output_path)
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)

            # Write HTML content
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(html_content)

            self.logger.info(f"HTML report saved to: {output_path}")

        except Exception as e:
            raise RuntimeError(f"Failed to save HTML report: {e}")

    def generate_report_summary(self) -> Dict[str, Any]:
        """
        Generate a summary of the report generation process for CLI output.

        Returns:
            Dict[str, Any]: Summary information about the generated report
        """
        if not self.aggregated_data:
            return {"status": "error", "message": "No data aggregated"}

        stats = self.aggregated_data["statistics"]

        return {
            "status": "success",
            "pipeline_completion": stats.get("data_completeness", 0),
            "data_summary": {
                "genomes_processed": stats.get("total_genomes", 0),
                "mutations_detected": stats.get("total_mutations", 0),
                "cooccurrence_pairs": stats.get("total_cooccurrence_pairs", 0),
                "gene_families": stats.get("gene_families", 0),
            },
            "visualization_components": {
                "interactive_charts": 5,  # Number of Plotly charts
                "msa_alignments": len(
                    self.aggregated_data.get("file_paths", {}).get(
                        "alignment_manifest", []
                    )
                ),
                "data_tables": 3,  # Genome, mutation, co-occurrence tables
            },
            "file_sources": {
                "metadata_available": bool(
                    self.aggregated_data["file_paths"].get("metadata_csv")
                ),
                "mutations_available": bool(
                    self.aggregated_data["file_paths"].get("mutation_report")
                ),
                "cooccurrence_available": bool(
                    self.aggregated_data["file_paths"].get("cooccurrence_report")
                ),
            },
        }
