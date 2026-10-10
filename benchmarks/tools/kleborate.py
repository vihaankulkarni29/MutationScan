"""Kleborate wrapper.

Kleborate is a *per-genome* summary: one row per
sample, with species, phylogroup, sequence type,
capsule/O-locus, virulence and resistance scores,
and a wide block of resistance-gene columns. It is
the benchmark's source of lineage labels (ST and
phylogroup) and of acquired-gene calls, but it
calls no point mutations, so its mutation records
are gene-presence records.

The resistance-gene block is version-dependent and
wide, so the parser does not hard-code its column
names: it treats every column that carries a
``yes``/allele value (and is not a known metadata
column) as a detected gene.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from ..schemas import LineageRecord, MutationRecord, MutationType, ToolOutput
from .base import Tool

__all__ = ["Kleborate"]

#: Columns that describe the sample, not a gene.
_METADATA_COLUMNS = {
    "strain",
    "species",
    "phylogroup",
    "st",
    "sequence_type",
    "kl_locus",
    "kl_locus_length",
    "kl_allele",
    "kl_match",
    "o_locus",
    "o_allele",
    "o_match",
    "k_var",
    "k_var_match",
    "virulence_score",
    "resistance_score",
    "yersiniabactin",
    "colibactin",
    "aerobactin",
    "salmochelin",
    "rmpa",
    "rmpa2",
    "peg_344",
    "iro",
    "iuc",
    "clb",
    "ybt",
}

#: Values that mean "this gene is present".
_PRESENT_VALUES = {"yes", "y", "true", "1", "present", "detected"}


def _normalise(column: str) -> str:
    return column.strip().lstrip("#").strip().lower().replace(" ", "_").replace("-", "_")


class Kleborate(Tool):
    """Wrap ``kleborate``."""

    def __init__(self, binary: str = "kleborate", **kwargs: Any) -> None:
        super().__init__(name="Kleborate", binary=binary, **kwargs)

    def command(self, inputs: dict[str, Path], output_dir: Path) -> list[str]:
        genome = inputs["genome"]
        out = output_dir / "kleborate.tsv"
        cmd = [self.binary, str(genome), "-o", str(out)]
        cmd += self.extra_args
        return cmd

    def parse(self, output_dir: Path, **context: Any) -> ToolOutput:
        path = output_dir / "kleborate.tsv"
        mutations: list[MutationRecord] = []
        lineages: list[LineageRecord] = []
        if not path.is_file():
            return ToolOutput(tool=self.name, mutations=mutations, lineages=lineages)

        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for raw in reader:
                row = {_normalise(key): value for key, value in raw.items() if key}
                strain = row.get("strain", "") or str(context.get("genome", ""))

                st = row.get("st", "") or row.get("sequence_type", "")
                phylogroup = row.get("phylogroup", "")
                if st:
                    lineages.append(LineageRecord(accession=strain, lineage=f"ST{st}"))
                if phylogroup:
                    lineages.append(LineageRecord(accession=strain, lineage=phylogroup))

                for column, value in row.items():
                    if column in _METADATA_COLUMNS or not value:
                        continue
                    if value.strip().lower() not in _PRESENT_VALUES:
                        continue
                    mutations.append(
                        MutationRecord(
                            accession=strain,
                            gene=column,
                            mutation=column,
                            mutation_type=MutationType.OTHER,
                            ref_aa="",
                            position=None,
                            alt_aa="",
                        )
                    )
        return ToolOutput(tool=self.name, mutations=mutations, lineages=lineages)
