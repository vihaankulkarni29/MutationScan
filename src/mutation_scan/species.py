"""Species-aware reference routing.

Why this exists
---------------
A cohort is usually more than one species, but a single reference panel is not
valid across species. Aligning an *E. coli* K-12 reference protein against a
*Klebsiella* genome lands near 80-90% identity -- above an 80% gate -- and then
every ordinary difference between the two species is reported as an acquired
mutation. Measured on a 570-genome three-species cohort, 97.7% of the reported
mutations came from pairs below 95% identity, and 94.6% of the mutations
attributed to the E. coli genes came from *Klebsiella* genomes aligned to
E. coli references.

This is not a paralog problem. Reciprocal best hit would not catch it, and would
in fact *pass* those pairs: E. coli acrB and K. pneumoniae acrB are true
orthologs. The fix has to be to compare each genome against a reference from its
own species, and to say so rather than silently substituting another species'
protein.

Configuration
--------------
::

    species:
      escherichia_coli:
        genome_pattern: "Escherichia_coli_*"
        references: "E_coli_k12"
      klebsiella_pneumoniae:
        genome_pattern: "Klebsiella_pneumoniae_*"
        references: "K_pneumoniae_MGH 78578"

``references`` is a path relative to the run's ``references_dir``. Genomes are
matched by ``fnmatch`` on their filename stem, so no annotation file is needed --
but the mapping is explicit and checked rather than inferred, because a wrong
guess here is the exact failure this module exists to prevent.

A gene with no reference in a genome's own species is reported as unavailable
for that genome. It is never quietly substituted with another species' protein.
"""

from __future__ import annotations

import fnmatch
import logging
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

__all__ = [
    "SPECIES_UNASSIGNED",
    "SpeciesConfigError",
    "SpeciesPanel",
    "SpeciesProfile",
    "parse_species_panel",
]

#: Recorded for a genome no species pattern matched.
SPECIES_UNASSIGNED = "unassigned"

_VALID_KEYS = frozenset({"genome_pattern", "references"})


class SpeciesConfigError(ValueError):
    """Raised when the ``species`` config block is malformed."""


@dataclass(frozen=True)
class SpeciesProfile:
    """One species: how to recognise its genomes, and where its references live."""

    key: str
    genome_pattern: str
    #: Directory holding this species' references, resolved against the run's
    #: ``references_dir``.
    references_dir: Path


@dataclass(frozen=True)
class SpeciesPanel:
    """Every species in a run, plus the mapping from genomes to species.

    Empty ``profiles`` means the run is not species-aware, and every genome is
    compared against the single flat reference set. That is the historical
    behaviour and stays available, but it is only correct for a single-species
    cohort.
    """

    profiles: tuple[SpeciesProfile, ...] = ()

    @property
    def enabled(self) -> bool:
        return bool(self.profiles)

    @property
    def keys(self) -> tuple[str, ...]:
        return tuple(profile.key for profile in self.profiles)

    def profile_for(self, species_key: str) -> SpeciesProfile | None:
        for profile in self.profiles:
            if profile.key == species_key:
                return profile
        return None

    def species_for(self, genome: Path | str) -> str:
        """Which species a genome belongs to, or :data:`SPECIES_UNASSIGNED`.

        Patterns are matched against the filename stem. When two species patterns
        both match, the more specific one wins -- the pattern with more literal
        characters -- so a narrow pattern can carve an exception out of a broad
        one without the user having to order the config.
        """
        stem = Path(genome).stem
        best: tuple[int, str] | None = None
        for profile in self.profiles:
            if not fnmatch.fnmatch(stem, profile.genome_pattern):
                continue
            # More literal characters means a more specific pattern, so it wins a
            # tie. Counting wildcards would not do: "Alpha_*" and "*" both hold
            # exactly one, and they are not equally specific.
            literal = sum(1 for char in profile.genome_pattern if char not in "*?[")
            if best is None or literal > best[0]:
                best = (literal, profile.key)
        return best[1] if best is not None else SPECIES_UNASSIGNED

    def assign(self, genomes: Iterable[Path | str]) -> dict[str, str]:
        """Map every genome path to a species key."""
        return {str(genome): self.species_for(genome) for genome in genomes}

    def grouped(self, genomes: Iterable[Path | str]) -> dict[str, list[Path]]:
        """Group genomes by species, dropping the unassigned ones."""
        groups: dict[str, list[Path]] = {}
        for genome in genomes:
            path = Path(genome)
            groups.setdefault(self.species_for(path), []).append(path)
        return groups

    def unassigned(self, genomes: Iterable[Path | str]) -> list[Path]:
        return [
            Path(genome) for genome in genomes if self.species_for(genome) == SPECIES_UNASSIGNED
        ]


def parse_species_panel(
    value: Any,
    references_dir: Path | str,
) -> SpeciesPanel:
    """Parse the ``species`` config block into a :class:`SpeciesPanel`.

    Args:
        value: the raw mapping from the config file. ``None`` or ``{}`` yields a
            disabled panel, which preserves single-reference-set behaviour.
        references_dir: the run's reference root, against which each species'
            ``references`` path is resolved.

    Raises:
        SpeciesConfigError: the block is not a mapping, or an entry is missing a
            key or uses one that is not recognised.
    """
    if value is None or value == {} or value == []:
        return SpeciesPanel()
    if not isinstance(value, Mapping):
        raise SpeciesConfigError("'species' must be a mapping of species name to its settings.")

    root = Path(references_dir)
    profiles: list[SpeciesProfile] = []
    seen: set[str] = set()

    for key, settings in value.items():
        name = str(key).strip()
        if not name:
            raise SpeciesConfigError("A 'species' key must not be empty.")
        if name in seen:
            raise SpeciesConfigError(f"Duplicate species {name!r} in 'species'.")
        seen.add(name)
        if name == SPECIES_UNASSIGNED:
            raise SpeciesConfigError(
                f"{SPECIES_UNASSIGNED!r} is reserved for genomes that match no species."
            )
        if not isinstance(settings, Mapping):
            raise SpeciesConfigError(f"'species.{name}' must be a mapping.")

        unknown = sorted(set(settings) - _VALID_KEYS)
        if unknown:
            raise SpeciesConfigError(
                f"Unknown key(s) in 'species.{name}': {', '.join(unknown)}. "
                f"Valid keys: {', '.join(sorted(_VALID_KEYS))}."
            )

        pattern = str(settings.get("genome_pattern", "")).strip()
        if not pattern:
            raise SpeciesConfigError(f"'species.{name}.genome_pattern' must not be empty.")

        references = str(settings.get("references", "")).strip()
        if not references:
            raise SpeciesConfigError(f"'species.{name}.references' must not be empty.")

        reference_path = Path(references).expanduser()
        if not reference_path.is_absolute():
            reference_path = root / reference_path

        profiles.append(
            SpeciesProfile(
                key=name,
                genome_pattern=pattern,
                references_dir=reference_path,
            )
        )

    if not profiles:
        raise SpeciesConfigError("'species' was given but lists no species.")

    return SpeciesPanel(profiles=tuple(profiles))


def union_genes(
    panel: SpeciesPanel, genes_by_species: Mapping[str, Sequence[str]]
) -> tuple[str, ...]:
    """Every gene any species provides, sorted, de-duplicated case-insensitively.

    The union rather than the intersection: a gene that only some species have is
    still worth analysing for those species, and the pairs it cannot serve are
    reported as unavailable instead of being dropped from the run.
    """
    names: dict[str, str] = {}
    for genes in genes_by_species.values():
        for gene in genes:
            name = str(gene).strip()
            if name:
                names.setdefault(name.lower(), name)
    return tuple(names[key] for key in sorted(names))
