"""Mykrobe wrapper.

``mykrobe predict`` writes a JSON document whose
``samples[].genotype.calls`` list carries one entry per
gene/variant the tool genotyped. Only calls whose
genotype is ``MUT`` (or otherwise non-reference) are
mutations; ``WT`` calls are the tool asserting the
wild type, which is exactly the negative evidence a
specificity score needs.

The parser therefore keeps both: ``MUT`` calls become
mutation records, and the full call list is retained in
``metadata`` so the harness can derive true negatives
from the genes Mykrobe examined and found clean.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from ..schemas import MutationRecord, MutationType, ToolOutput
from .base import Tool, ToolError

__all__ = ["Mykrobe"]

#: Genotype strings that mean "the variant is present".
_MUT_GENOTYPES = {"MUT", "mut", "MUTANT", "mutant", "1", 1, "ALT", "alt"}


class Mykrobe(Tool):
    """Wrap ``mykrobe predict``."""

    def __init__(self, binary: str = "mykrobe", **kwargs: Any) -> None:
        super().__init__(name="Mykrobe", binary=binary, **kwargs)

    def command(self, inputs: dict[str, Path], output_dir: Path) -> list[str]:
        genome = inputs["genome"]
        out = output_dir / "mykrobe.json"
        species = str(inputs.get("species", "auto"))
        cmd = [
            self.binary,
            "predict",
            str(genome),
            species,
            "--output",
            str(out),
            "--format",
            "json",
        ]
        cmd += self.extra_args
        return cmd

    def parse(self, output_dir: Path, **context: Any) -> ToolOutput:
        path = output_dir / "mykrobe.json"
        records: list[MutationRecord] = []
        metadata: dict[str, Any] = {}
        if not path.is_file():
            return ToolOutput(tool=self.name, mutations=records, metadata=metadata)

        try:
            with path.open(encoding="utf-8") as handle:
                document = json.load(handle)
        except (json.JSONDecodeError, OSError) as exc:
            raise ToolError(f"Mykrobe: unreadable output {path}: {exc}") from exc

        samples = document.get("samples", [])
        if isinstance(samples, dict):
            samples = [samples]

        for sample in samples:
            accession = sample.get("sample_name", "") or str(context.get("genome", ""))
            genotype = sample.get("genotype", {})
            calls = genotype.get("calls", [])
            metadata.setdefault("calls", []).extend(
                {
                    "gene": call.get("gene", ""),
                    "variant": call.get("variant", ""),
                    "genotype": call.get("genotype", ""),
                }
                for call in calls
            )
            for call in calls:
                if call.get("genotype") not in _MUT_GENOTYPES:
                    continue
                gene = str(call.get("gene", ""))
                variant = str(call.get("variant", ""))
                if not gene:
                    continue
                records.append(
                    MutationRecord(
                        accession=accession,
                        gene=gene,
                        mutation=variant or gene,
                        mutation_type=MutationType.OTHER,
                        ref_aa="",
                        position=None,
                        alt_aa="",
                    )
                )
        return ToolOutput(tool=self.name, mutations=records, metadata=metadata)
