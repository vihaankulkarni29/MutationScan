"""Production entrypoint that delegates execution to Snakemake.

Usage:
    python -m mutation_scan [snakemake args]
"""

import shutil
import subprocess
import sys


def main() -> int:
    snakemake_bin = shutil.which("snakemake")
    if not snakemake_bin:
        sys.stderr.write(
            "snakemake is not available on PATH. Run the pipeline with "
            "`python -m snakemake ...` or install snakemake in the active environment.\n"
        )
        return 1

    cmd = [snakemake_bin, *sys.argv[1:]]
    return subprocess.call(cmd)


if __name__ == "__main__":
    raise SystemExit(main())
