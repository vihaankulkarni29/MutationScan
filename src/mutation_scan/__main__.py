"""``python -m mutation_scan`` -- the same program as the ``mutationscan`` command.

Kept as a one-line delegation so there is exactly one command-line implementation
to maintain, and so the tool is usable without the console script being on PATH.
"""

from .cli import main

if __name__ == "__main__":
    raise SystemExit(main())
