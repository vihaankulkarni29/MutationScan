#!/usr/bin/env python3
"""
Domino 3 (Extractor) Smoke Test

This quick test verifies that FastaAAExtractor is available and can produce an output
FASTA on a tiny synthetic input. It works with either the real CLI or the local mock.

Exit codes:
  0 - success, extractor works
  1 - failure, extractor unavailable or run failed
"""
import os
import sys
import tempfile
import shutil
import subprocess
from pathlib import Path
from typing import Optional


def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)


def run_cmd(cmd: list[str]) -> tuple[int, str, str]:
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    tools_dir = Path(__file__).resolve().parent

    # Prefer console script, then python -m fallback, else use local mock
    extractor_exe = which("fasta_aa_extractor")
    module_fallback = [sys.executable, "-m", "fasta_aa_extractor", "--help"]
    mock_script = tools_dir / "mock_fasta_aa_extractor.py"

    # Check availability
    print("Checking FastaAAExtractor availability...")
    if extractor_exe:
        rc, out, err = run_cmd([extractor_exe, "--help"])
        if rc != 0:
            print("[FAIL] Console script exists but --help failed")
            print(err)
            return 1
        print("[OK] Console script 'fasta_aa_extractor' is available")
        runner = [extractor_exe]
    else:
        rc, out, err = run_cmd(module_fallback)
        if rc == 0:
            print("[OK] Module runner 'python -m fasta_aa_extractor' is available")
            runner = [sys.executable, "-m", "fasta_aa_extractor"]
        elif mock_script.is_file():
            print("[WARN] Real extractor not found; using local mock for smoke test")
            runner = [sys.executable, str(mock_script)]
        else:
            print("[FAIL] FastaAAExtractor not found (neither console script nor module); mock is also unavailable")
            return 1

    # Prepare tiny synthetic inputs
    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        genome = td_path / "genome.fna"
        tsv = td_path / "coords.tsv"
        genes = td_path / "genes.txt"
        outfaa = td_path / "out.faa"

        genome.write_text(">seq1\nATGAAATTTAA\n")
        tsv.write_text("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\n" \
                        "genome.fna\tseq1\t1\t9\t+\tblaTEM\n")
        genes.write_text("blaTEM\n")

        cmd = [*runner,
               "--genome", str(genome),
               "--coordinates", str(tsv),
               "--genes", str(genes),
               "--output", str(outfaa)]

        print("Running extractor smoke test...")
        rc, out, err = run_cmd(cmd)
        if rc != 0:
            print("[FAIL] Extractor returned non-zero exit code")
            if err:
                print(err)
            return 1

        if not outfaa.is_file():
            print("[FAIL] Output FASTA not created")
            return 1

        print("[OK] Output FASTA created:", outfaa)
        # Quick content check
        content = outfaa.read_text(errors="ignore")
        if ">" not in content:
            print("[FAIL] Output FASTA missing sequence headers")
            return 1

        print("[PASS] Extractor smoke test successful")
        return 0


if __name__ == "__main__":
    sys.exit(main())
