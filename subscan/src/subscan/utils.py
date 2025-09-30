"""
SubScan Pipeline Utilities

This module provides shared utility functions used across all domino tools
to eliminate code redundancy and improve maintainability.

Author: MutationScan Development Team
"""

import json
import os
from typing import Dict, Any, List, Optional
import logging


class ManifestError(Exception):
    """Custom exception for manifest-related errors."""


def load_json_manifest(
    manifest_path: str, required_keys: Optional[List[str]] = None
) -> Dict[str, Any]:
    """
    Load and validate a JSON manifest file with comprehensive error handling.

    This function provides a centralized way to load manifest files across all domino tools,
    with consistent error handling and validation patterns.

    Args:
        manifest_path: Path to the JSON manifest file
        required_keys: List of keys that must be present in the manifest

    Returns:
        Dict[str, Any]: Parsed manifest data

    Raises:
        FileNotFoundError: If manifest file doesn't exist
        ManifestError: If manifest format is invalid or missing required keys
        PermissionError: If file cannot be read due to permissions
    """
    if not os.path.isfile(manifest_path):
        raise FileNotFoundError(f"Manifest file not found: {manifest_path}")

    try:
        with open(manifest_path, "r", encoding="utf-8") as f:
            manifest_data = json.load(f)
    except PermissionError as e:
        raise PermissionError(
            f"Cannot read manifest file (permission denied): {manifest_path}"
        ) from e
    except json.JSONDecodeError as e:
        raise ManifestError(
            f"Invalid JSON in manifest file '{manifest_path}': {e}"
        ) from e
    except Exception as e:
        raise ManifestError(
            f"Failed to load manifest file '{manifest_path}': {e}"
        ) from e

    # Validate required keys if specified
    if required_keys:
        missing_keys = [key for key in required_keys if key not in manifest_data]
        if missing_keys:
            raise ManifestError(
                f"Manifest '{manifest_path}' missing required keys: {missing_keys}"
            )

    return manifest_data


def save_json_manifest(manifest_data: Dict[str, Any], output_path: str) -> None:
    """
    Save manifest data to a JSON file with proper formatting and error handling.

    Args:
        manifest_data: Dictionary containing manifest data
        output_path: Path where the manifest should be saved

    Raises:
        PermissionError: If file cannot be written due to permissions
        OSError: If directory cannot be created or file cannot be written
    """
    # Ensure output directory exists
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        try:
            os.makedirs(output_dir, exist_ok=True)
        except OSError as e:
            raise OSError(f"Cannot create output directory '{output_dir}': {e}") from e

    try:
        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(manifest_data, f, indent=2, ensure_ascii=False)
    except PermissionError as e:
        raise PermissionError(
            f"Cannot write to manifest file (permission denied): {output_path}"
        ) from e
    except Exception as e:
        raise OSError(f"Failed to save manifest file '{output_path}': {e}") from e


def extract_genomes_from_manifest(
    manifest_data: Dict[str, Any],
) -> List[Dict[str, Any]]:
    """
    Extract genome list from manifest data, handling different manifest formats.

    This function standardizes how genome data is extracted from manifests across
    different domino tools, which may have varying manifest structures.

    Args:
        manifest_data: Parsed manifest dictionary

    Returns:
        List[Dict[str, Any]]: List of genome dictionaries

    Raises:
        ManifestError: If genomes cannot be found in expected locations
    """
    # Try different common manifest structures
    if "output_files" in manifest_data and "genomes" in manifest_data["output_files"]:
        # Domino 1 (Harvester) format
        genomes = manifest_data["output_files"]["genomes"]
    elif "genomes" in manifest_data:
        # Direct genomes format
        genomes = manifest_data["genomes"]
    elif (
        "analysis_results" in manifest_data
        and "genomes" in manifest_data["analysis_results"]
    ):
        # Analysis results format (from analyzer)
        genomes = manifest_data["analysis_results"]["genomes"]
    else:
        raise ManifestError(
            "Cannot find genomes in manifest. Expected keys: 'genomes', 'output_files.genomes', or 'analysis_results.genomes'"
        )

    if not isinstance(genomes, list):
        raise ManifestError(f"Expected genomes to be a list, got: {type(genomes)}")

    return genomes


def setup_logging(tool_name: str, log_level: str = "INFO") -> logging.Logger:
    """
    Set up standardized logging configuration for domino tools.

    Args:
        tool_name: Name of the domino tool (used as logger name)
        log_level: Logging level (DEBUG, INFO, WARNING, ERROR)

    Returns:
        logging.Logger: Configured logger instance
    """
    logger = logging.getLogger(tool_name)
    logger.setLevel(getattr(logging, log_level.upper()))

    # Create console handler if not already exists
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(getattr(logging, log_level.upper()))

        # Create formatter
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        handler.setFormatter(formatter)

        logger.addHandler(handler)

    return logger


def validate_file_exists(file_path: str, file_description: str = "File") -> None:
    """
    Validate that a file exists with descriptive error messages.

    Args:
        file_path: Path to the file to check
        file_description: Human-readable description of what the file is

    Raises:
        FileNotFoundError: If file doesn't exist
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"{file_description} not found: {file_path}")


def create_output_directory(output_dir: str) -> None:
    """
    Create output directory with proper error handling.

    Args:
        output_dir: Directory path to create

    Raises:
        OSError: If directory cannot be created
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as e:
        raise OSError(f"Cannot create output directory '{output_dir}': {e}") from e


def format_duration(seconds: float) -> str:
    """
    Format duration in seconds to human-readable format.

    Args:
        seconds: Duration in seconds

    Returns:
        str: Formatted duration string (e.g., "1h 23m 45s" or "2m 34s")
    """
    hours = int(seconds // 3600)
    minutes = int((seconds % 3600) // 60)
    secs = int(seconds % 60)

    if hours > 0:
        return f"{hours}h {minutes:02d}m {secs:02d}s"
    return f"{minutes}m {secs:02d}s"
