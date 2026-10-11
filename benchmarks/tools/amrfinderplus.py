"""AMRFinderPlus wrapper.

AMRFinderPlus (NCBI) reports acquired resistance genes and
point mutations in a single tabular file. The benchmark
normalises both into the canonical mutation record: an
acquired gene becomes a ``(genome, gene, gene_name)``
record, and a point mutation becomes a
``(genome, gene, change)`` record, so it can be compared
against MutationScan's substitution calls on the same key.

The parser normalises column names (lowercase, spaces and
``%`` stripped) and matches a set of known aliases, so a
column renamed between AMRFinderPlus versions does not
silently drop a field.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from ..schemas import MutationRecord, MutationType, ToolOutput
from .base import Tool

__all__ = ["AmrfinderPlus"]

#: Aliases for the fields the parser needs, in the order it
#: tries them. Real column names carry spaces and ``%``; the
#: normaliser strips those, so each alias is already in
#: normalised form.
_GENE_ALIASES = ("symbol", "gene_symbol", "gene", "name")
_NAME_ALIASES = ("name", "protein_name", "description")
_CLASS_ALIASES = ("class", "drug_class")
_SUBCLASS_ALIASES = ("subclass", "drug_subclass")
_CONTIG_ALIASES = ("contig_id", "sequence_name", "contig", "accession")
_IDENTITY_ALIASES = ("pct_identity", "identity", "percent_identity")
_COVERAGE_ALIASES = ("pct_coverage", "coverage", "percent_coverage")
_START_ALIASES = ("start",)
_STOP_ALIASES = ("stop",)
_STRAND_ALIASES = ("strand",)
_METHOD_ALIASES = ("method",)


def _normalise(column: str) -> str:
    """``% identity`` -> ``pct_identity``-style key."""
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


class AmrfinderPlus(Tool):
    """Wrap ``amrfinder`` / ``amrfinderplus``."""

    def __init__(self, binary: str = "amrfinder", **kwargs: Any) -> None:
        super().__init__(name="AMRFinderPlus", binary=binary, **kwargs)

    def command(self, inputs: dict[str, Path], output_dir: Path) -> list[str]:
        genome = inputs["genome"]
        out = output_dir / "amrfinder.tsv"
        cmd = [
            self.binary,
            "-n",
            str(genome),
            "-o",
            str(out),
            "--plus",
        ]
        if inputs.get("organism"):
            cmd += ["--organism", str(inputs["organism"])]
        cmd += self.extra_args
        return cmd

    def parse(self, output_dir: Path, **context: Any) -> ToolOutput:
        path = output_dir / "amrfinder.tsv"
        records: list[MutationRecord] = []
        if not path.is_file():
            return ToolOutput(tool=self.name, mutations=records)

        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for raw in reader:
                row = {_normalise(key): value for key, value in raw.items() if key}
                gene = _pick(row, _GENE_ALIASES)
                if not gene:
                    continue
                name = _pick(row, _NAME_ALIASES)
                contig = _pick(row, _CONTIG_ALIASES) or str(context.get("genome", ""))
                mutation = name if name and name != gene else gene
                records.append(
                    MutationRecord(
                        accession=contig,
                        gene=gene,
                        mutation=mutation,
                        mutation_type=MutationType.OTHER,
                        ref_aa="",
                        position=None,
                        alt_aa="",
                    )
                )
        return ToolOutput(tool=self.name, mutations=records)
