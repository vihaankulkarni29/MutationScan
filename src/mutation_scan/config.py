"""Run configuration for MutationScan.

The config is deliberately small: where the genomes are, where the reference
proteins are, and the quality gates to apply. Target genes are normally *not*
configured -- they are discovered from the reference filenames (see
:func:`mutation_scan.references.discover_targets`), which is what keeps the tool
free of any hardcoded gene, organism or drug.

All run outputs live under a single folder, ``{output_root}/{job_name}/``, and
every path in that folder is derived from :class:`Config` properties so that the
CLI and the Snakemake workflow cannot drift apart.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml

__all__ = [
    "Config",
    "CooccurrenceConfig",
    "ConfigError",
    "VALID_COOCCURRENCE_LEVELS",
    "load_config",
]


class ConfigError(ValueError):
    """Raised when a run configuration is missing, malformed or contradictory."""


#: Resolutions at which co-occurrence can be counted.
VALID_COOCCURRENCE_LEVELS = ("gene", "variant")

_DEFAULTS: dict[str, Any] = {
    "job_name": "default_run",
    "genomes_dir": "data/genomes",
    "references_dir": "data/references",
    "output_root": "data/output",
    "targets": [],
    "uniprot_taxid": "",
    "min_identity_percent": 80.0,
    "min_coverage_percent": 80.0,
    "tblastn_binary": "tblastn",
    "threads": 1,
    "cooccurrence": {"levels": ["gene"], "min_count": 1},
}


@dataclass(frozen=True)
class CooccurrenceConfig:
    """How to count co-occurrence. Counting only -- no scoring, no statistics."""

    levels: tuple[str, ...] = ("gene",)
    min_count: int = 1


@dataclass(frozen=True)
class Config:
    """A fully resolved, validated run configuration."""

    job_name: str
    genomes_dir: Path
    references_dir: Path
    output_root: Path
    targets: tuple[str, ...]
    uniprot_taxid: str
    min_identity_percent: float
    min_coverage_percent: float
    tblastn_binary: str
    threads: int
    cooccurrence: CooccurrenceConfig

    # -- derived run layout -------------------------------------------------
    # One folder per run; nothing else writes outside it.

    @property
    def run_dir(self) -> Path:
        return self.output_root / self.job_name

    @property
    def refs_dir(self) -> Path:
        """Run-local copy of the references actually used (provenance)."""
        return self.run_dir / "refs"

    @property
    def proteins_dir(self) -> Path:
        return self.run_dir / "proteins"

    @property
    def manifest_csv(self) -> Path:
        return self.run_dir / "extraction_manifest.csv"

    @property
    def mutations_csv(self) -> Path:
        return self.run_dir / "mutations.csv"

    @property
    def variant_qc_csv(self) -> Path:
        """One row per (genome, gene) pair considered by variant calling.

        Distinguishes "not mutated" from "not evaluated", and is what supplies
        honest denominators to the co-occurrence stage.
        """
        return self.run_dir / "variant_qc.csv"

    @property
    def cooccurrence_genes_csv(self) -> Path:
        return self.run_dir / "cooccurrence_genes.csv"

    @property
    def cooccurrence_matrix_csv(self) -> Path:
        return self.run_dir / "cooccurrence_genes_matrix.csv"

    @property
    def cooccurrence_variants_csv(self) -> Path:
        return self.run_dir / "cooccurrence_variants.csv"

    @property
    def run_summary_json(self) -> Path:
        return self.run_dir / "run_summary.json"

    def outputs(self) -> tuple[Path, ...]:
        """Files a complete run is expected to produce."""
        paths = [
            self.manifest_csv,
            self.mutations_csv,
            self.variant_qc_csv,
            self.cooccurrence_genes_csv,
            self.cooccurrence_matrix_csv,
            self.run_summary_json,
        ]
        if "variant" in self.cooccurrence.levels:
            paths.append(self.cooccurrence_variants_csv)
        return tuple(paths)

    # -- construction -------------------------------------------------------

    @classmethod
    def from_mapping(
        cls,
        mapping: Mapping[str, Any] | None,
        base_dir: Path | str | None = None,
    ) -> "Config":
        """Build a config from a plain mapping (YAML contents, CLI args, ...).

        Relative paths are resolved against ``base_dir`` (default: the current
        working directory, which is what Snakemake sets to the workdir).
        """
        raw = dict(_DEFAULTS)
        raw["cooccurrence"] = dict(_DEFAULTS["cooccurrence"])

        supplied = dict(mapping or {})
        unknown = sorted(set(supplied) - set(_DEFAULTS))
        if unknown:
            raise ConfigError(
                f"Unknown config key(s): {', '.join(unknown)}. "
                f"Valid keys: {', '.join(sorted(_DEFAULTS))}."
            )

        cooccurrence_raw = supplied.pop("cooccurrence", None) or {}
        if not isinstance(cooccurrence_raw, Mapping):
            raise ConfigError("'cooccurrence' must be a mapping of levels/min_count.")
        unknown_co = sorted(set(cooccurrence_raw) - set(_DEFAULTS["cooccurrence"]))
        if unknown_co:
            raise ConfigError(
                f"Unknown cooccurrence key(s): {', '.join(unknown_co)}. "
                f"Valid keys: {', '.join(sorted(_DEFAULTS['cooccurrence']))}."
            )
        raw["cooccurrence"].update(cooccurrence_raw)
        # Drop keys explicitly set to null so YAML blanks fall back to defaults.
        raw.update({k: v for k, v in supplied.items() if v is not None})

        root = Path(base_dir) if base_dir is not None else Path.cwd()

        def _resolve(key: str) -> Path:
            value = str(raw[key]).strip()
            if not value:
                raise ConfigError(f"'{key}' must not be empty.")
            path = Path(value).expanduser()
            return path if path.is_absolute() else (root / path)

        job_name = str(raw["job_name"]).strip()
        if not job_name:
            raise ConfigError("'job_name' must not be empty.")
        if Path(job_name).name != job_name or job_name in {".", ".."}:
            raise ConfigError(
                f"'job_name' must be a single folder name, got {job_name!r}."
            )

        config = cls(
            job_name=job_name,
            genomes_dir=_resolve("genomes_dir"),
            references_dir=_resolve("references_dir"),
            output_root=_resolve("output_root"),
            targets=_clean_targets(raw["targets"]),
            uniprot_taxid=str(raw["uniprot_taxid"]).strip(),
            min_identity_percent=_percent(raw["min_identity_percent"], "min_identity_percent"),
            min_coverage_percent=_percent(raw["min_coverage_percent"], "min_coverage_percent"),
            tblastn_binary=str(raw["tblastn_binary"]).strip() or "tblastn",
            threads=_positive_int(raw["threads"], "threads"),
            cooccurrence=CooccurrenceConfig(
                levels=_clean_levels(raw["cooccurrence"]["levels"]),
                min_count=_positive_int(raw["cooccurrence"]["min_count"], "cooccurrence.min_count"),
            ),
        )
        return config

    def with_overrides(self, **overrides: Any) -> "Config":
        """Return a copy with fields replaced (used by the CLI)."""
        return replace(self, **{k: v for k, v in overrides.items() if v is not None})

    def require_inputs(self) -> None:
        """Fail loudly, and early, on missing inputs.

        Kept separate from parsing so that ``--help``, dry-runs and unit tests
        can build a config without the data being present.
        """
        for label, path in (
            ("genomes_dir", self.genomes_dir),
            ("references_dir", self.references_dir),
        ):
            if not path.exists():
                raise ConfigError(f"{label} does not exist: {path}")
            if not path.is_dir():
                raise ConfigError(f"{label} is not a directory: {path}")

    def to_dict(self) -> dict[str, Any]:
        """Round-trippable plain-data view, for run_summary.json provenance."""
        return {
            "job_name": self.job_name,
            "genomes_dir": str(self.genomes_dir),
            "references_dir": str(self.references_dir),
            "output_root": str(self.output_root),
            "targets": list(self.targets),
            "uniprot_taxid": self.uniprot_taxid,
            "min_identity_percent": self.min_identity_percent,
            "min_coverage_percent": self.min_coverage_percent,
            "tblastn_binary": self.tblastn_binary,
            "threads": self.threads,
            "cooccurrence": {
                "levels": list(self.cooccurrence.levels),
                "min_count": self.cooccurrence.min_count,
            },
        }


def load_config(path: Path | str, base_dir: Path | str | None = None) -> Config:
    """Load and validate a YAML config file."""
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f"Config file not found: {path}")

    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"Invalid YAML in {path}: {exc}") from exc

    if data is None:
        data = {}
    if not isinstance(data, Mapping):
        raise ConfigError(f"Config root must be a mapping, got {type(data).__name__}: {path}")

    return Config.from_mapping(data, base_dir=base_dir)


# -- validation helpers -----------------------------------------------------


def _percent(value: Any, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"'{label}' must be a number, got {value!r}.") from exc
    if not 0.0 <= number <= 100.0:
        raise ConfigError(f"'{label}' must be between 0 and 100, got {number}.")
    return number


def _positive_int(value: Any, label: str) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"'{label}' must be an integer, got {value!r}.") from exc
    if number < 1:
        raise ConfigError(f"'{label}' must be >= 1, got {number}.")
    return number


def _clean_targets(value: Any) -> tuple[str, ...]:
    if isinstance(value, str):
        items: Iterable[Any] = [part for part in value.replace(",", " ").split()]
    elif isinstance(value, Iterable):
        items = value
    else:
        raise ConfigError(f"'targets' must be a list of gene names, got {value!r}.")

    cleaned: list[str] = []
    for item in items:
        name = str(item).strip()
        if not name:
            continue
        if name not in cleaned:
            cleaned.append(name)
    return tuple(cleaned)


def _clean_levels(value: Any) -> tuple[str, ...]:
    if isinstance(value, str):
        items: Iterable[Any] = value.replace(",", " ").split()
    elif isinstance(value, Iterable):
        items = value
    else:
        raise ConfigError(f"'cooccurrence.levels' must be a list, got {value!r}.")

    cleaned: list[str] = []
    for item in items:
        level = str(item).strip().lower()
        if not level:
            continue
        if level not in VALID_COOCCURRENCE_LEVELS:
            raise ConfigError(
                f"Unknown cooccurrence level {level!r}. "
                f"Valid levels: {', '.join(VALID_COOCCURRENCE_LEVELS)}."
            )
        if level not in cleaned:
            cleaned.append(level)

    if not cleaned:
        raise ConfigError("'cooccurrence.levels' must list at least one level.")
    # Stable, predictable ordering regardless of how the user wrote it.
    return tuple(l for l in VALID_COOCCURRENCE_LEVELS if l in cleaned)
