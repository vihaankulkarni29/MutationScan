#!/usr/bin/env python3
"""MutationScan Pipeline Orchestrator (clean, ASCII-safe)

Provides two modes:
  1. Full execution of the 7 domino tools.
  2. --dry-run structural validation producing placeholder manifests + report.

Domino chain & expected manifests:
  1 Harvester     -> genome_manifest.json
  2 Annotator     -> annotation_manifest.json
  3 Extractor     -> protein_manifest.json
  4 Aligner       -> alignment_manifest.json
  5 Analyzer      -> analysis_manifest.json
  6 CoOccurrence  -> cooccurrence_manifest.json
  7 Reporter      -> final_report.html (or mutation_analysis_report.html)
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import List
import shutil

TOTAL_STAGES = 7

class DominoToolError(Exception):
    """Raised when a domino tool returns non-zero exit."""
    pass

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="MutationScan pipeline orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python run_pipeline.py --accessions acc.txt --gene-list genes.txt --email you@lab.org --output-dir run --sepi-species 'Escherichia coli'\n"
            "  python run_pipeline.py --accessions acc.txt --gene-list genes.txt --email you@lab.org --output-dir dry --sepi-species 'Escherichia coli' --dry-run"
        ),
    )
    req = p.add_argument_group("Required")
    req.add_argument("--accessions", required=True)
    req.add_argument("--gene-list", required=True)
    req.add_argument("--email", required=True)
    req.add_argument("--output-dir", required=True)
    ref = p.add_mutually_exclusive_group(required=True)
    ref.add_argument("--user-reference-dir")
    ref.add_argument("--sepi-species")
    opt = p.add_argument_group("Optional")
    opt.add_argument("--threads", type=int, default=4)
    opt.add_argument("--verbose", action="store_true")
    opt.add_argument("--open-report", action="store_true")
    opt.add_argument("--dry-run", action="store_true")
    return p

def validate_inputs(a) -> None:
    if not os.path.isfile(a.accessions):
        raise FileNotFoundError(f"Accessions file missing: {a.accessions}\nPlease provide a valid file with NCBI accession numbers, one per line.")
    if not os.path.isfile(a.gene_list):
        raise FileNotFoundError(f"Gene list file missing: {a.gene_list}\nPlease provide a valid file with gene names, one per line.")
    if a.user_reference_dir and not os.path.isdir(a.user_reference_dir):
        raise FileNotFoundError(f"Reference dir missing: {a.user_reference_dir}\nPlease provide a valid directory containing reference sequences.")
    if "@" not in a.email:
        raise ValueError("Invalid email (must contain '@'). Please enter a valid email address for NCBI queries.")
    # Validate accession file format
    with open(a.accessions, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]
        if not lines:
            raise ValueError(f"Accessions file '{a.accessions}' is empty. Please add at least one NCBI accession number.")
        for i, line in enumerate(lines, 1):
            if not (line.startswith("NZ_") or line.startswith("CP") or line.startswith("BA") or line.isalnum()):
                print(f"Warning: Line {i} in accessions file may not be a valid NCBI accession: '{line}'")
    # Validate gene list format
    with open(a.gene_list, "r", encoding="utf-8") as f:
        genes = [line.strip() for line in f if line.strip()]
        if not genes:
            raise ValueError(f"Gene list file '{a.gene_list}' is empty. Please add at least one gene name.")
        for i, gene in enumerate(genes, 1):
            if not gene.isalnum():
                print(f"Warning: Line {i} in gene list may not be a valid gene name: '{gene}'")

def make_run_dir(base: str) -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    rd = os.path.join(base, f"mutationscan_run_{ts}")
    os.makedirs(rd, exist_ok=True)
    return rd

def ensure_manifest(out_dir: str, name: str) -> str:
    path = os.path.join(out_dir, name)
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Expected manifest not found: {path}")
    return path

def run_tool(script: Path, args_list: List[str], label: str, step: int, verbose: bool) -> None:
    print(f"\n=== Domino {step}/{TOTAL_STAGES}: {label} ===")
    cmd = [sys.executable, str(script)] + args_list
    if verbose:
        print("  Command:", " ".join(cmd))
    start = time.time()
    try:
        subprocess.run(cmd, check=True, text=True)
        print(f"✅ {label} completed successfully in {time.time() - start:.1f}s.")
    except subprocess.CalledProcessError as e:
        print(f"❌ ERROR in {label} (Domino {step}): Tool failed with exit code {e.returncode}.", file=sys.stderr)
        print(f"  Command: {' '.join(cmd)}", file=sys.stderr)
        print("  Please check the output directory and logs for details.", file=sys.stderr)
        print("  Refer to the README troubleshooting section for help.", file=sys.stderr)
        raise DominoToolError(f"{label} failed (exit {e.returncode})") from e

def dry_run(run_dir: str) -> str:
    print("Dry-run: generating placeholder manifests")
    stages = [
        ("01_harvester_results", "genome_manifest.json"),
        ("02_annotator_results", "annotation_manifest.json"),
        ("03_extractor_results", "protein_manifest.json"),
        ("04_aligner_results", "alignment_manifest.json"),
        ("05_analyzer_results", "analysis_manifest.json"),
        ("06_cooccurrence_results", "cooccurrence_manifest.json"),
    ]
    prev = None
    for idx, (folder, manifest) in enumerate(stages, start=1):
        od = os.path.join(run_dir, folder)
        os.makedirs(od, exist_ok=True)
        mp = os.path.join(od, manifest)
        with open(mp, "w", encoding="utf-8") as fh:
            json.dump({
                "stage": idx,
                "manifest": manifest,
                "placeholder": True,
                "input_manifest": prev,
                "timestamp": datetime.utcnow().isoformat() + "Z"
            }, fh, indent=2)
        print(f"  [dry-run] {manifest}")
        prev = mp
    reporter_dir = os.path.join(run_dir, "07_reporter_results")
    os.makedirs(reporter_dir, exist_ok=True)
    report = os.path.join(reporter_dir, "final_report.html")
    with open(report, "w", encoding="utf-8") as fh:
        fh.write("<html><body><h1>MutationScan Dry-Run</h1><p>Placeholder.</p></body></html>")
    print("  [dry-run] final_report.html")
    print("Dry-run complete")
    return report

def execute_pipeline(a, run_dir: str) -> str:
    if a.dry_run:
        return dry_run(run_dir)
    tools = Path(__file__).parent
    current = None
    domino_labels = [
        "Harvester (Genome Download)",
        "Annotator (AMR Gene Identification)",
        "Extractor (Sequence Extraction)",
        "Aligner (Reference Alignment)",
        "Analyzer (Mutation Analysis)",
        "CoOccurrence (Pattern Analysis)",
        "Reporter (HTML Dashboard)"
    ]
    # 1 Harvester
    out1 = os.path.join(run_dir, "01_harvester_results")
    run_tool(tools / "run_harvester.py", ["--accessions", a.accessions, "--email", a.email, "--output-dir", out1], domino_labels[0], 1, a.verbose)
    current = ensure_manifest(out1, "genome_manifest.json")
    # 2 Annotator
    out2 = os.path.join(run_dir, "02_annotator_results")
    run_tool(tools / "run_annotator.py", ["--manifest", current, "--output-dir", out2, "--threads", str(a.threads)], domino_labels[1], 2, a.verbose)
    current = ensure_manifest(out2, "annotation_manifest.json")
    # 3 Extractor
    out3 = os.path.join(run_dir, "03_extractor_results")
    # Pre-flight: ensure FastaAAExtractor is available (console or module)
    if not shutil.which("fasta_aa_extractor"):
        # Try module help to see if import works
        try:
            cp = subprocess.run([sys.executable, "-m", "fasta_aa_extractor", "--help"], capture_output=True, text=True)
            if cp.returncode != 0:
                print("\nERROR: FastaAAExtractor is not available in the current environment.")
                print("Please install it so the extractor can run:")
                print("  pip install -e .  # from the extractor project")
                print("Or ensure the console script 'fasta_aa_extractor' is on PATH.")
                raise DominoToolError(3, "Extractor (Sequence Extraction)")
        except FileNotFoundError:
            print("\nERROR: Python executable not found to attempt 'python -m fasta_aa_extractor'.")
            raise DominoToolError(3, "Extractor (Sequence Extraction)")
    run_tool(tools / "run_extractor.py", ["--manifest", current, "--gene-list", a.gene_list, "--output-dir", out3, "--threads", str(a.threads)], domino_labels[2], 3, a.verbose)
    current = ensure_manifest(out3, "protein_manifest.json")
    # 4 Aligner
    out4 = os.path.join(run_dir, "04_aligner_results")
    align_args = ["--manifest", current, "--output-dir", out4, "--threads", str(a.threads)]
    if a.sepi_species:
        align_args += ["--sepi-species", a.sepi_species]
    else:
        align_args += ["--user-reference-dir", a.user_reference_dir]
    run_tool(tools / "run_aligner.py", align_args, domino_labels[3], 4, a.verbose)
    current = ensure_manifest(out4, "alignment_manifest.json")
    # 5 Analyzer
    out5 = os.path.join(run_dir, "05_analyzer_results")
    run_tool(tools / "run_analyzer.py", ["--manifest", current, "--output-dir", out5, "--threads", str(a.threads)], domino_labels[4], 5, a.verbose)
    current = ensure_manifest(out5, "analysis_manifest.json")
    # 6 CoOccurrence
    out6 = os.path.join(run_dir, "06_cooccurrence_results")
    run_tool(tools / "run_cooccurrence_analyzer.py", ["--manifest", current, "--output-dir", out6, "--threads", str(a.threads)], domino_labels[5], 6, a.verbose)
    current = ensure_manifest(out6, "cooccurrence_manifest.json")
    # 7 Reporter
    out7 = os.path.join(run_dir, "07_reporter_results")
    run_tool(tools / "run_reporter.py", ["--manifest", current, "--output-dir", out7], domino_labels[6], 7, a.verbose)
    for name in ("final_report.html", "mutation_analysis_report.html"):
        candidate = os.path.join(out7, name)
        if os.path.isfile(candidate):
            print(f"\n🎉 All domino stages completed! Your results are ready: {candidate}")
            return candidate
    print(f"\n⚠️  Final report not found. Check {out7} for output files.")
    return out7

def main():
    p = build_parser()
    a = p.parse_args()
    print("MutationScan Pipeline Orchestrator")
    print("=" * 50)
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("Mode:", "DRY-RUN" if a.dry_run else "FULL")
    start = time.time()
    try:
        print("Validating inputs...")
        validate_inputs(a)

        # Pre-flight check for all domino tool scripts
        domino_scripts = [
            "run_harvester.py",
            "run_annotator.py",
            "run_extractor.py",
            "run_aligner.py",
            "run_analyzer.py",
            "run_cooccurrence_analyzer.py",
            "run_reporter.py",
        ]
        missing = []
        tools_dir = Path(__file__).parent
        for script in domino_scripts:
            if not (tools_dir / script).is_file():
                missing.append(str(tools_dir / script))
        # Check for external tool: federated_genome_extractor
        federated_extractor_path = tools_dir.parent / "federated_genome_extractor" / "harvester.py"
        if not federated_extractor_path.is_file():
            missing.append(str(federated_extractor_path))
        if missing:
            print("\nERROR: The following required domino tool scripts are missing:")
            for m in missing:
                print(f"  - {m}")
            print("\nPlease ensure all domino tools are installed and available in the expected locations.")
            print("Refer to the README for installation instructions.")
            sys.exit(10)

        run_dir = make_run_dir(a.output_dir)
        print("Run directory:", run_dir)
        report = execute_pipeline(a, run_dir)
        elapsed = time.time() - start
        h, rem = divmod(elapsed, 3600)
        m, s = divmod(rem, 60)
        print("\nPipeline complete")
        print(f"Runtime: {int(h):02d}:{int(m):02d}:{int(s):02d}")
        if os.path.isfile(report):
            print("Report:", report)
            if a.open_report and not a.dry_run:
                try:
                    import webbrowser  # type: ignore
                    webbrowser.open(f"file://{os.path.abspath(report)}")
                except Exception:
                    pass
        else:
            print("Report directory:", report)
    except DominoToolError as e:
        print(str(e), file=sys.stderr)
        sys.exit(3)
    except (FileNotFoundError, ValueError) as e:
        print(f"Configuration error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("Interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:  # noqa: BLE001
        print(f"Unexpected error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(2)

if __name__ == "__main__":  # pragma: no cover
    main()
