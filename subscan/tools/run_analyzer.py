#!/usr/bin/env python3
"""Wrapper script for Domino 5 (Analyzer) to integrate external subscan-analyzer package.

This lightweight adapter allows the orchestrator to call a stable interface
without depending on internal repo layout of the analyzer project.

It attempts to import the `subscan_analyzer` entry point from the installed
`subscan-analyzer` package (installed via git+ URL in pyproject). If unavailable,
it provides a clear error with remediation guidance.
"""
import sys
import argparse
import importlib
import json
import os
from datetime import datetime


def create_argument_parser():
    parser = argparse.ArgumentParser(
        description="Domino 5: Analyzer - Mutation detection & analysis wrapper",
    )
    parser.add_argument("--manifest", required=True, help="Path to alignment_manifest.json from Domino 4")
    parser.add_argument("--output-dir", required=True, help="Directory to write analyzer outputs")
    parser.add_argument("--threads", type=int, default=4, help="Worker threads")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--dry-run", action="store_true", help="Generate placeholder analysis outputs without heavy processing")
    return parser


def load_analyzer_main():
    candidates = [
        ("subscan_analyzer.cli", "main"),
        ("subscan_analyzer", "main"),
    ]
    for module_name, attr in candidates:
        try:
            mod = importlib.import_module(module_name)
            if hasattr(mod, attr):
                return getattr(mod, attr)
        except ModuleNotFoundError:
            continue
    return None


def run_dry_run(manifest_path: str, output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    synthetic_manifest = {
        "stage": 5,
        "domino": "Analyzer",
        "input_manifest": manifest_path,
        "generated_at": datetime.utcnow().isoformat() + "Z",
        "placeholder": True,
        "notes": "Synthetic analysis manifest generated in --dry-run mode",
        "summary": {"mutations_detected": 0, "genomes_processed": 0},
    }
    out_manifest = os.path.join(output_dir, "analysis_manifest.json")
    with open(out_manifest, "w", encoding="utf-8") as f:
        json.dump(synthetic_manifest, f, indent=2)
    # create a placeholder excel/csv to satisfy downstream expectations
    with open(os.path.join(output_dir, "subscan_mutation_report.xlsx"), "wb") as f:
        f.write(b"Placeholder XLSX (dry-run)")
    return out_manifest


def main():
    parser = create_argument_parser()
    args = parser.parse_args()

    if args.dry_run:
        manifest = run_dry_run(args.manifest, args.output_dir)
        if args.verbose:
            print(f"[dry-run] Created placeholder analysis manifest at {manifest}")
        return 0

    analyzer_entry = load_analyzer_main()
    if analyzer_entry is None:
        print(
            "ERROR: Could not locate the subscan-analyzer package entry point.\n"
            "Ensure installation succeeded (git dependency) and that your environment is active.\n"
            "Try: pip install -e .  (from the project root).",
            file=sys.stderr,
        )
        return 2

    # Delegate to analyzer
    try:
        return analyzer_entry([
            "--manifest", args.manifest,
            "--output-dir", args.output_dir,
            "--threads", str(args.threads),
        ] + (["--verbose"] if args.verbose else []))
    except Exception as e:  # noqa: BLE001
        print(f"Analyzer execution failed: {e}", file=sys.stderr)
        return 3


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
