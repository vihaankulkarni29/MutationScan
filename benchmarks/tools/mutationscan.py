"""MutationScan wrapper.

MutationScan is the reference implementation the
benchmark scores against the external tools, so its
wrapper is also the harness's self-test: running
MutationScan through the same ``Tool.run`` path as
AMRFinderPlus et al. proves the metric code is
symmetric before it is pointed at a foreign output.

It runs the full pipeline and reads the four CSVs
the pipeline writes -- ``mutations.csv``,
``epistasis_networks.csv``, ``dn_ds.csv`` and
``genome_mutation_map.csv`` -- into the canonical
records. A ``lineage_manifest`` input turns on the
lineage annotation, so the lineage-concordance
metric has something to score.
"""

from __future__ import annotations

import csv
import re
import shutil
import sys
from pathlib import Path
from typing import Any

from ..schemas import (
    DnDsRecord,
    EpistasisRecord,
    LineageRecord,
    MutationRecord,
    MutationType,
    ToolOutput,
)
from .base import Tool

__all__ = ["MutationScanTool"]

#: A substitution looks like ``S83L`` or ``A12V``: one or two
#: amino-acid letters around a position.
_SUBSTITUTION = re.compile(r"^[A-Z*]\d+[A-Z*]$")

_JOB_NAME = "benchmark"


def _classify(mutation: str) -> MutationType:
    if mutation == "LoF":
        return MutationType.LoF
    if _SUBSTITUTION.match(mutation):
        return MutationType.SNV
    return MutationType.OTHER


def _float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _int(value: Any, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


class MutationScanTool(Tool):
    """Wrap ``mutationscan run``.

    The console script is not always on ``PATH`` (an
    editable install may skip it), so the wrapper
    prefers ``python -m mutation_scan``, which works
    wherever the package is importable -- which is
    everywhere the benchmark can be run at all.
    """

    def __init__(self, binary: str = "mutationscan", **kwargs: Any) -> None:
        super().__init__(name="MutationScan", binary=binary, **kwargs)

    @property
    def available(self) -> bool:
        if shutil.which(self.binary):
            return True
        import importlib.util

        return importlib.util.find_spec("mutation_scan") is not None

    def _argv(self, args: list[str]) -> list[str]:
        if shutil.which(self.binary):
            return [self.binary, *args]
        return [sys.executable, "-m", "mutation_scan", *args]

    def command(self, inputs: dict[str, Path], output_dir: Path) -> list[str]:
        genomes = inputs["genomes"]
        references = inputs["references"]
        args = [
            "run",
            "--genomes",
            str(genomes),
            "--references",
            str(references),
            "--out",
            str(output_dir),
            "--job-name",
            _JOB_NAME,
        ]
        if inputs.get("lineage_manifest"):
            args += ["--lineage-manifest", str(inputs["lineage_manifest"])]
        if inputs.get("threads"):
            args += ["--threads", str(inputs["threads"])]
        args += self.extra_args
        return self._argv(args)

    def parse(self, output_dir: Path, **context: Any) -> ToolOutput:
        run_dir = output_dir / _JOB_NAME
        mutations = self._parse_mutations(run_dir / "mutations.csv")
        epistasis = self._parse_epistasis(run_dir / "epistasis_networks.csv")
        dn_ds = self._parse_dn_ds(run_dir / "dn_ds.csv")
        lineages = self._parse_lineages(run_dir / "genome_mutation_map.csv")
        return ToolOutput(
            tool=self.name,
            mutations=mutations,
            epistasis=epistasis,
            dn_ds=dn_ds,
            lineages=lineages,
        )

    def _parse_mutations(self, path: Path) -> list[MutationRecord]:
        records: list[MutationRecord] = []
        if not path.is_file():
            return records
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                gene = (row.get("Gene") or "").strip()
                mutation = (row.get("Mutation") or "").strip()
                if not gene or not mutation:
                    continue
                records.append(
                    MutationRecord(
                        accession=(row.get("Accession") or "").strip(),
                        gene=gene,
                        mutation=mutation,
                        mutation_type=_classify(mutation),
                        ref_aa=(row.get("Ref_AA") or "").strip(),
                        position=_int(row.get("Position")) if row.get("Position") else None,
                        alt_aa=(row.get("Alt_AA") or "").strip(),
                    )
                )
        return records

    def _parse_epistasis(self, path: Path) -> list[EpistasisRecord]:
        records: list[EpistasisRecord] = []
        if not path.is_file():
            return records
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                variant_a = (row.get("Variant_A") or "").strip()
                variant_b = (row.get("Variant_B") or "").strip()
                if not variant_a or not variant_b:
                    continue
                records.append(
                    EpistasisRecord(
                        variant_a=variant_a,
                        variant_b=variant_b,
                        gene_a=(row.get("Gene_A") or "").strip(),
                        gene_b=(row.get("Gene_B") or "").strip(),
                        n_both=_int(row.get("N_Both")),
                        q_value=_float(row.get("Q_value"), 1.0),
                        sign=(row.get("Sign") or "neutral").strip(),
                    )
                )
        return records

    def _parse_dn_ds(self, path: Path) -> list[DnDsRecord]:
        records: list[DnDsRecord] = []
        if not path.is_file():
            return records
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                gene = (row.get("Gene") or "").strip()
                if not gene:
                    continue
                dn_ds = row.get("dN_dS")
                records.append(
                    DnDsRecord(
                        gene=gene,
                        n_sites=_int(row.get("N_Sites")),
                        s_sites=_int(row.get("S_Sites")),
                        dn=_float(row.get("dN")),
                        ds=_float(row.get("dS")),
                        dn_ds=_float(dn_ds) if dn_ds not in (None, "", "inf") else None,
                        n_genomes=_int(row.get("N_Genomes_Evaluated")),
                    )
                )
        return records

    def _parse_lineages(self, path: Path) -> list[LineageRecord]:
        records: list[LineageRecord] = []
        if not path.is_file():
            return records
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                accession = (row.get("Accession") or "").strip()
                lineage = (row.get("Lineage") or "").strip()
                if accession and lineage:
                    records.append(LineageRecord(accession=accession, lineage=lineage))
        return records
