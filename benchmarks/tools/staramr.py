"""StarAMR wrapper.

StarAMR is a pipeline over ResFinder (acquired
genes), PointFinder (chromosomal point mutations)
and PlasmidFinder, and writes one TSV per
database into its output directory. The two that
matter for mutation benchmarking are
``resfinder_results.tsv`` (gene presence) and
``pointfinder_results.tsv`` (point mutations such
as ``gyrA S83L``) -- the latter is the closest
thing any of the external tools produce to
MutationScan's substitution calls, so it is the
one the mutation-level metrics lean on.

Both files are tab-delimited with a ``#``-prefixed
header; the parser normalises column names and
matches aliases, as the other wrappers do.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from ..schemas import MutationRecord, MutationType, ToolOutput
from .base import Tool

__all__ = ["StarAMR"]

_GENE_ALIASES = ("gene", "gene_name", "gene_symbol")
_MUTATION_ALIASES = ("mutation", "mutation_description", "change", "variant")
_RESISTANCE_ALIASES = ("resistance", "drug_class", "class")


def _normalise(column: str) -> str:
    return (
        column.strip()
        .lstrip("#")
        .strip()
        .lower()
        .replace(" ", "_")
        .replace("%", "pct_")
        .replace("-", "_")
    )


def _pick(row: dict[str, str], aliases: tuple[str, ...]) -> str:
    for alias in aliases:
        if alias in row and row[alias]:
            return row[alias]
    return ""


class StarAMR(Tool):
    """Wrap ``staramr``."""

    def __init__(self, binary: str = "staramr", **kwargs: Any) -> None:
        super().__init__(name="StarAMR", binary=binary, **kwargs)

    def command(self, inputs: dict[str, Path], output_dir: Path) -> list[str]:
        genome = inputs["genome"]
        out = output_dir / "staramr"
        cmd = [
            self.binary,
            "search",
            "--output_dir",
            str(out),
            str(genome),
        ]
        cmd += self.extra_args
        return cmd

    def parse(self, output_dir: Path, **context: Any) -> ToolOutput:
        search_dir = output_dir / "staramr"
        records: list[MutationRecord] = []
        if not search_dir.is_dir():
            return ToolOutput(tool=self.name, mutations=records)

        for filename, is_point in (
            ("resfinder_results.tsv", False),
            ("pointfinder_results.tsv", True),
        ):
            path = search_dir / filename
            if not path.is_file():
                continue
            records.extend(self._parse_tsv(path, is_point, context))
        return ToolOutput(tool=self.name, mutations=records)

    def _parse_tsv(
        self,
        path: Path,
        is_point: bool,
        context: dict[str, Any],
    ) -> list[MutationRecord]:
        records: list[MutationRecord] = []
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for raw in reader:
                row = {_normalise(key): value for key, value in raw.items() if key}
                gene = _pick(row, _GENE_ALIASES)
                if not gene:
                    continue
                mutation = _pick(row, _MUTATION_ALIASES)
                if not is_point or not mutation:
                    mutation = gene
                records.append(
                    MutationRecord(
                        accession=str(context.get("genome", "")),
                        gene=gene,
                        mutation=mutation,
                        mutation_type=(MutationType.SNV if is_point else MutationType.OTHER),
                        ref_aa="",
                        position=None,
                        alt_aa="",
                    )
                )
        return records
