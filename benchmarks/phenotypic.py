"""Phenotypic concordance: does a called mutation predict drug resistance?

The strongest claim a variant caller can make is clinical: genomes
carrying the mutations it calls should be the genomes the laboratory
reports as resistant. That is the question this module answers, and it
needs two things the pipeline does not otherwise use:

* **MIC data** -- a measured minimum inhibitory concentration per
  (genome, drug). BV-BRC (formerly PATRIC) serves these; each also
  publishes which breakpoint standard the value is interpreted under.
* **Breakpoints** -- the MIC thresholds that dichotomise a
  concentration into a phenotype. CLSI M100 is the US standard and is
  paywalled; **EUCAST publishes the same breakpoints open-access**,
  which is why this module defaults to the EUCAST format. A CLSI table
  parses exactly the same way, so switching standards is a matter of
  which file is loaded.

The measurement is deliberately split in three, because they answer
different questions and a reader should not have to take one on trust
from the others:

1. **Gene-drug association** (unsupervised) -- for each (gene, drug)
   pair, is carrying a mutation in that gene associated with the
   phenotype? Two-sided Fisher exact on a 2x2 table, Benjamini-Hochberg
   across the gene x drug family. This needs no curated knowledge at
   all: it discovers from the cohort which genes matter for which
   drugs, so it is the check that MutationScan's calls carry *any*
   phenotypic signal.
2. **Burden-vs-MIC correlation** (unsupervised) -- does the number of
   resistance-associated mutations a genome carries track its MIC?
   Spearman per drug. This is the continuous form of the claim: it
   uses the whole dilution series rather than a dichotomised
   threshold, so it does not depend on the breakpoint being right.
3. **Concordance against curated rules** (supervised) -- given a
   hand-curated gene -> drug resistance rule set, the very-major and
   major error rates. VME (reported susceptible, actually resistant)
   is the dangerous direction; ME (reported resistant, actually
   susceptible) is the inconvenient one. Reported separately because
   they are not interchangeable.

A MIC value is interpreted only against the breakpoint for its
organism and drug. An organism/drug pair with no breakpoint means the
phenotype is *unknown*, which is reported as such rather than defaulted
to susceptible -- the distinction the whole module exists to preserve.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from collections.abc import Iterable, Iterator, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from scipy.stats import fisher_exact, spearmanr  # noqa: E402

logger = logging.getLogger(__name__)

__all__ = [
    "MicRecord",
    "Breakpoint",
    "GeneDrugAssociation",
    "PhenotypeRules",
    "ConcordanceResult",
    "SUSCEPTIBLE",
    "INTERMEDIATE",
    "RESISTANT",
    "UNKNOWN",
    "classify_mic",
    "read_mic_records",
    "read_breakpoints",
    "read_rules",
    "read_carriers_from_genome_map",
    "gene_drug_associations",
    "burden_correlation",
    "concordance",
    "write_report",
]

#: Phenotype categories. ``UNKNOWN`` is a real outcome, not a blank:
#: it means a MIC exists but no breakpoint was given for that
#: organism/drug, so the value cannot be dichotomised at all.
SUSCEPTIBLE = "S"
INTERMEDIATE = "I"
RESISTANT = "R"
UNKNOWN = "?"

#: Above this FDR, a gene-drug association is not reported as one.
DEFAULT_ALPHA = 0.05


# ---------------------------------------------------------------------------
# Schemas
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class MicRecord:
    """One measured MIC, in µg/mL, for one (genome, drug).

    ``organism`` is the key the breakpoint table is looked up under;
    when it is empty the breakpoint is matched on the drug alone.
    """

    accession: str
    drug: str
    mic: float
    organism: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        if not math.isfinite(self.mic) or self.mic <= 0:
            raise ValueError(
                f"MIC must be a positive finite number, got {self.mic!r} "
                f"for ({self.accession}, {self.drug})"
            )

    @property
    def log2_mic(self) -> float:
        """MIC on the log2 scale the dilution series actually lives on."""
        return math.log2(self.mic)

    @property
    def key(self) -> tuple[str, str]:
        return (self.accession, self.drug)


@dataclass(frozen=True)
class Breakpoint:
    """The MIC thresholds for one (organism, drug) pair.

    Both standards are expressed the same way: a maximum MIC that is
    still susceptible, an optional intermediate band above that, and a
    minimum MIC that is resistant. ``i_max`` being absent means there
    is no I category, so anything above ``s_max`` is resistant -- the
    EUCAST convention for drugs with no intermediate zone.

    ``standard`` records which document the numbers came from, so a
    report can say EUCAST or CLSI rather than implying one table is
    interchangeable with the other.
    """

    organism: str
    drug: str
    s_max: float | None = None
    i_max: float | None = None
    r_min: float | None = None
    standard: str = "EUCAST"


@dataclass
class GeneDrugAssociation:
    """One Fisher test: gene mutation carriage vs. resistant phenotype."""

    gene: str
    drug: str
    n_carrier: int = 0
    n_resistant: int = 0
    n_both: int = 0
    n_evaluable: int = 0
    odds_ratio: float = 0.0
    p_value: float = 1.0
    q_value: float = 1.0
    #: The FDR threshold this test was judged against, so a result can
    #: be read without consulting a module-level default that may not
    #: be the threshold the caller asked for.
    alpha: float = DEFAULT_ALPHA

    @property
    def significant(self) -> bool:
        return self.q_value < self.alpha

    def as_row(self) -> dict[str, Any]:
        return {
            "gene": self.gene,
            "drug": self.drug,
            "n_carrier": self.n_carrier,
            "n_resistant": self.n_resistant,
            "n_both": self.n_both,
            "n_evaluable": self.n_evaluable,
            "odds_ratio": (None if not math.isfinite(self.odds_ratio) else _round(self.odds_ratio)),
            "p_value": _format_p(self.p_value),
            "q_value": _format_p(self.q_value),
            "significant": self.significant,
        }


@dataclass
class PhenotypeRules:
    """A curated gene -> drug resistance mapping.

    ``resistance_genes`` maps a gene to the set of drugs a mutation in
    it is *known* to confer resistance to. This is the supervised half
    of the module: the association test discovers which genes matter,
    and this table states what the literature already says, so the two
    can be compared.
    """

    resistance_genes: Mapping[str, frozenset[str]] = field(default_factory=dict)

    def predicts_resistance(self, gene: str, drug: str) -> bool:
        """Whether a mutation in ``gene`` is known to confer resistance.

        Both sides are matched case-insensitively: a MIC table's
        ``ciprofloxacin`` and a rule table's ``Ciprofloxacin`` are the
        same drug, and a case mismatch should not read as a gene the
        rule set has simply never heard of.
        """
        return drug.strip().lower() in self.resistance_genes.get(gene, frozenset())

    def drugs_for(self, genes: Iterable[str]) -> set[str]:
        """Every drug this genome's mutations should confer resistance to."""
        predicted: set[str] = set()
        for gene in genes:
            predicted |= self.resistance_genes.get(gene, frozenset())
        return predicted


@dataclass
class ConcordanceResult:
    """How well a curated rule set tracks the measured phenotype."""

    n_evaluable: int = 0
    n_unknown_phenotype: int = 0
    agreement: float = 0.0
    #: Reported susceptible, actually resistant. The dangerous one.
    very_major_error: float = 0.0
    #: Reported resistant, actually susceptible. The inconvenient one.
    major_error: float = 0.0
    #: S vs I disagreements with no resistance implication.
    minor_error: float = 0.0
    n_very_major: int = 0
    n_major: int = 0
    n_minor: int = 0
    per_drug: dict[str, dict[str, Any]] = field(default_factory=dict)

    def as_row(self) -> dict[str, Any]:
        return {
            "n_evaluable": self.n_evaluable,
            "n_unknown_phenotype": self.n_unknown_phenotype,
            "agreement": _round(self.agreement),
            "very_major_error": _round(self.very_major_error),
            "major_error": _round(self.major_error),
            "minor_error": _round(self.minor_error),
            "n_very_major": self.n_very_major,
            "n_major": self.n_major,
            "n_minor": self.n_minor,
            "per_drug": self.per_drug,
        }


# ---------------------------------------------------------------------------
# Loaders
# ---------------------------------------------------------------------------


def _normalise_header(name: str) -> str:
    return name.strip().strip("#").strip().lower().replace(" ", "_")


def _is_comment(line: str) -> bool:
    """Whether a line is documentation rather than data.

    A benchmark table is far more useful if its provenance, its
    convention and its caveats can travel *inside* the file, and a
    ``#`` prefix is the least surprising way to do that. The loader
    has to skip them or the first line of every documented table
    reads as a header.
    """
    return line.lstrip().startswith("#")


def _is_blank(line: str) -> bool:
    """Whether a line carries no fields at all.

    ``csv.DictReader`` already drops blank rows, but *finding the
    header* has to skip them too: a documented table naturally has a
    blank line between its prose and its column names, and reading
    that as the header yields an empty one.
    """
    return not line.strip()


def _header_line(path: Path) -> str:
    """The first line of *path* that is neither a comment nor blank."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for line in handle:
            if not _is_comment(line) and not _is_blank(line):
                return line
    return ""


def _data_lines(path: Path) -> Iterator[str]:
    """Yield the lines of *path* that are neither comments nor blank."""
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for line in handle:
            if not _is_comment(line) and not _is_blank(line):
                yield line


def _dialect_for(path: Path) -> Any:
    """Pick a tab or comma dialect from the file's header line.

    ``csv.Sniffer`` raises "Could not determine delimiter" on exactly
    the tables this loader most needs to accept: a handful of rows, an
    empty optional column, or a single data line. Counting separators
    in the header is deterministic instead, and a header is guaranteed
    to hold every column name so both separators appear at least once
    when present.
    """
    header = _header_line(path)
    tabs = header.count("\t")
    commas = header.count(",")
    if tabs == 0 and commas == 0:
        # A single-column file: either dialect reads it the same way.
        return csv.excel_tab
    return csv.excel_tab if tabs >= commas else csv.excel


def _to_float(value: Any) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    # MIC tables use "<=0.25" and ">64" for off-dilution readings, and
    # "/" for a censored pair. All three collapse to a comparable
    # concentration; the inequality is recorded as a flag upstream if
    # a caller needs it, not silently reinterpreted here.
    text = text.lstrip("<=").lstrip(">=").strip()
    try:
        parsed = float(text)
    except ValueError:
        return None
    return parsed if parsed > 0 else None


def read_mic_records(path: Path | str) -> list[MicRecord]:
    """Read MIC rows from a TSV/CSV.

    Recognised columns: ``accession``/``genome``/``genome_id``,
    ``drug``/``antibiotic``, ``mic``/``mic_value``, ``organism``,
    ``source``. Header names are matched case-insensitively and ignore
    surrounding whitespace and spaces, so a BV-BRC export and a hand
    curated table parse with the same code.
    """
    path = Path(path)
    dialect = _dialect_for(path)
    reader = csv.DictReader(_data_lines(path), dialect=dialect)
    if reader.fieldnames is None:
        return []
    columns = {_normalise_header(name): name for name in reader.fieldnames}

    def column(*aliases: str) -> str | None:
        for alias in aliases:
            if alias in columns:
                return columns[alias]
        return None

    accession_col = column("accession", "genome", "genome_id", "assembly")
    drug_col = column("drug", "antibiotic", "compound")
    mic_col = column("mic", "mic_value", "mic_ug_ml")
    organism_col = column("organism", "species")
    source_col = column("source", "lab", "study")

    missing = [
        label
        for label, col in (("accession", accession_col), ("drug", drug_col), ("mic", mic_col))
        if col is None
    ]
    if missing:
        raise SystemExit(
            f"{path}: MIC table needs {', '.join(missing)} column(s); "
            f"found {sorted(reader.fieldnames)}"
        )

    records: list[MicRecord] = []
    for row in reader:
        accession = (row.get(accession_col) or "").strip()
        drug = (row.get(drug_col) or "").strip()
        mic = _to_float(row.get(mic_col))
        if not accession or not drug or mic is None:
            # An unparseable MIC is skipped, not coerced to 0: the
            # alternative is inventing a susceptible genome.
            continue
        organism = (row.get(organism_col) or "").strip() if organism_col else ""
        source = (row.get(source_col) or "").strip() if source_col else ""
        records.append(
            MicRecord(accession=accession, drug=drug, mic=mic, organism=organism, source=source)
        )
    return records


def read_breakpoints(path: Path | str) -> dict[tuple[str, str], Breakpoint]:
    """Read a breakpoint table, keyed by (organism, drug).

    Recognised columns: ``organism`` (optional), ``drug``,
    ``s_max``/``s``/``susceptible_mic``, ``i_max``/``i``,
    ``r_min``/``r``/``resistant_mic``, ``standard``.

    A breakpoint with neither an ``s_max`` nor an ``r_min`` is
    dropped rather than interpreted: with no threshold there is
    nothing to classify against, and defaulting to "susceptible"
    would manufacture agreement.
    """
    path = Path(path)
    dialect = _dialect_for(path)
    reader = csv.DictReader(_data_lines(path), dialect=dialect)
    if reader.fieldnames is None:
        return {}
    columns = {_normalise_header(name): name for name in reader.fieldnames}

    def column(*aliases: str) -> str | None:
        for alias in aliases:
            if alias in columns:
                return columns[alias]
        return None

    drug_col = column("drug", "antibiotic", "compound")
    s_col = column("s_max", "s", "susceptible_mic", "s_breakpoint")
    i_col = column("i_max", "i", "intermediate_mic", "i_breakpoint")
    r_col = column("r_min", "r", "resistant_mic", "r_breakpoint")
    organism_col = column("organism", "species")
    standard_col = column("standard", "guideline", "source")

    if drug_col is None or (s_col is None and r_col is None):
        raise SystemExit(
            f"{path}: breakpoint table needs a 'drug' column and at least one of "
            f"'s_max' / 'r_min'; found {sorted(reader.fieldnames)}"
        )

    table: dict[tuple[str, str], Breakpoint] = {}
    for row in reader:
        drug = (row.get(drug_col) or "").strip()
        if not drug:
            continue
        organism = (row.get(organism_col) or "").strip() if organism_col else ""
        standard = (row.get(standard_col) or "").strip() if standard_col else "EUCAST"
        breakpoint = Breakpoint(
            organism=organism,
            drug=drug,
            s_max=_to_float(row.get(s_col)) if s_col else None,
            i_max=_to_float(row.get(i_col)) if i_col else None,
            r_min=_to_float(row.get(r_col)) if r_col else None,
            standard=standard or "EUCAST",
        )
        if breakpoint.s_max is None and breakpoint.r_min is None:
            logger.warning("skipping %s/%s: no threshold to classify against", organism, drug)
            continue
        table[(organism.lower(), drug.lower())] = breakpoint
    return table


def read_rules(path: Path | str) -> PhenotypeRules:
    """Read a curated gene -> drug resistance mapping.

    Expected columns: ``gene`` and ``drug``, one resistance
    relationship per row. An optional ``allele`` column is read but not
    required: this module's genotype is gene-level, matching the
    collapsed ``LoF`` convention the pipeline uses everywhere else.
    """
    path = Path(path)
    resistance: dict[str, set[str]] = {}
    seen: set[tuple[str, str]] = set()
    reader = csv.DictReader(_data_lines(path), dialect=_dialect_for(path))
    if reader.fieldnames is None:
        raise SystemExit(f"{path}: empty rules table")
    columns = {_normalise_header(name): name for name in reader.fieldnames}
    gene_col = columns.get("gene")
    drug_col = columns.get("drug")
    if gene_col is None or drug_col is None:
        raise SystemExit(f"{path}: rules table needs 'gene' and 'drug' columns")
    for row in reader:
        gene = (row.get(gene_col) or "").strip()
        drug = (row.get(drug_col) or "").strip()
        if not gene or not drug:
            continue
        pair = (gene, drug.lower())
        if pair in seen:
            continue
        seen.add(pair)
        resistance.setdefault(gene, set()).add(drug)
    return PhenotypeRules(
        resistance_genes={
            gene: frozenset(drug.lower() for drug in drugs)
            for gene, drugs in sorted(resistance.items())
        }
    )


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------


def classify_mic(mic: float, breakpoint: Breakpoint | None) -> str:
    """Dichotomise a MIC into one of S / I / R / ?.

    ``?`` is returned when there is no breakpoint, when there is a
    threshold but the MIC falls inside a band the standard leaves
    undefined, and when a MIC is non-positive. A caller that treats
    ``?`` as susceptible will report perfect agreement against a
    cohort with no breakpoints, which is why it is a distinct outcome
    rather than a default.

    The two standard conventions are both supported, and the
    difference between them matters:

    * **Two-sided** (CLSI M100 form): ``s_max``, ``i_max`` and
      ``r_min`` all given. A MIC between ``i_max`` and ``r_min`` is a
      gap the standard does not define, so it is ``?``.
    * **EUCAST form**: only ``s_max``, as ``S <= x, R > x``, with no
      intermediate category at all. Here anything above ``s_max`` is
      resistant -- the absence of an ``i_max`` is the statement that
      there is no middle band, not that the standard is silent.

    ``r_min`` alone is the third shape, where anything below it is
    susceptible.
    """
    if breakpoint is None:
        return UNKNOWN
    if not math.isfinite(mic) or mic <= 0:
        return UNKNOWN

    s_max = breakpoint.s_max
    i_max = breakpoint.i_max
    r_min = breakpoint.r_min

    if s_max is not None and mic <= s_max:
        return SUSCEPTIBLE
    if r_min is not None and mic >= r_min:
        return RESISTANT
    if i_max is not None and mic <= i_max:
        # Above s_max but at or below i_max: the intermediate band.
        return INTERMEDIATE

    # Not susceptible, not resistant, not in an intermediate band. The
    # reading now depends on which convention the breakpoint came from.
    if i_max is None:
        if s_max is not None and r_min is None:
            # EUCAST's "S <= x, R > x": no middle band exists.
            return RESISTANT
        if s_max is None and r_min is not None:
            return SUSCEPTIBLE
    # Above i_max and below r_min, or a gap with no i_max to explain
    # it: the standard does not say, and guessing invents a phenotype.
    return UNKNOWN


def _lookup_breakpoint(
    table: Mapping[tuple[str, str], Breakpoint],
    record: MicRecord,
) -> Breakpoint | None:
    organism = record.organism.lower()
    drug = record.drug.lower()
    # (organism, drug) is the intended lookup. The empty-organism entry
    # is the fallback for a table that is not species-specific, and it
    # is only consulted when the organism-specific lookup missed, so a
    # per-species table can keep one shared row for drugs it applies to
    # all of them.
    if organism:
        match = table.get((organism, drug))
        if match is not None:
            return match
    return table.get(("", drug))


# ---------------------------------------------------------------------------
# Unsupervised: gene-drug association
# ---------------------------------------------------------------------------


def gene_drug_associations(
    carriers: Mapping[str, Mapping[str, bool]],
    records: Sequence[MicRecord],
    breakpoints: Mapping[tuple[str, str], Breakpoint],
    alpha: float = DEFAULT_ALPHA,
) -> list[GeneDrugAssociation]:
    """Fisher-test gene mutation carriage against resistant phenotype.

    ``carriers`` maps accession -> {gene: carries a called mutation}.
    ``records`` are the measured MICs, and ``breakpoints`` the
    thresholds that turn them into phenotypes.

    The family is every (gene, drug) pair with at least one carrier and
    one resistant genome, so each test is reported with the denominator
    it was actually computed on. Benjamini-Hochberg is applied across
    that whole family, because the number of tests is genes x drugs and
    an uncorrected Fisher at n=20 will happily call significance on
    chance alone.

    Genes with no breakpoint for their drug are simply absent: their
    phenotype is unknown, and an unknown phenotype cannot be tested
    for association.
    """
    # One phenotype per record, with a recognised organism.
    phenotyped: list[tuple[str, str, str]] = []
    for record in records:
        phenotype = classify_mic(record.mic, _lookup_breakpoint(breakpoints, record))
        if phenotype in (SUSCEPTIBLE, RESISTANT):
            phenotyped.append((record.accession, record.drug, phenotype))

    by_drug: dict[str, list[tuple[str, str]]] = {}
    for accession, drug, phenotype in phenotyped:
        by_drug.setdefault(drug.lower(), []).append((accession, phenotype))

    # A genome contributes to a drug's test only if its phenotype is
    # known; a genome with no mutation calls and no MIC contributes
    # nothing and is not counted as a negative.
    associations: list[GeneDrugAssociation] = []
    for drug, rows in sorted(by_drug.items()):
        evaluable = sorted(rows)
        genomes = [accession for accession, _ in evaluable]
        resistant = {accession for accession, phenotype in evaluable if phenotype == RESISTANT}

        # Genes are the union across all genomes that have a phenotype
        # for this drug, so a gene mutated only in genomes missing a
        # MIC for this drug does not create a test with no data.
        genes_in_drug = {
            gene
            for accession in genomes
            for gene, carrier in carriers.get(accession, {}).items()
            if carrier
        }
        if not genes_in_drug or not resistant:
            continue

        for gene in sorted(genes_in_drug):
            carrier_genomes = {
                accession for accession in genomes if carriers.get(accession, {}).get(gene, False)
            }
            if not carrier_genomes:
                continue
            n_both = len(carrier_genomes & resistant)
            # 2x2: [carrier & resistant, carrier & susceptible]
            #       [non-carrier & resistant, non-carrier & susceptible]
            n_carrier_s = len(carrier_genomes) - n_both
            n_non_carrier_r = len(resistant) - n_both
            n_non_carrier_s = len(evaluable) - n_both - n_carrier_s - n_non_carrier_r

            odds_ratio, p_value = fisher_exact(
                [[n_both, n_carrier_s], [n_non_carrier_r, n_non_carrier_s]],
                alternative="two-sided",
            )

            associations.append(
                GeneDrugAssociation(
                    gene=gene,
                    drug=drug,
                    n_carrier=len(carrier_genomes),
                    n_resistant=len(resistant),
                    n_both=n_both,
                    n_evaluable=len(evaluable),
                    odds_ratio=float(odds_ratio),
                    p_value=float(p_value),
                    q_value=1.0,
                )
            )

    _apply_bh(associations, alpha)
    return associations


def _apply_bh(associations: Sequence[GeneDrugAssociation], alpha: float) -> None:
    """Benjamini-Hochberg across the (gene, drug) family, in place.

    Sets each item's ``q_value`` and records the threshold it was
    judged against, so ``significant`` reflects what the caller asked
    for rather than a module-level default.
    """
    if not associations:
        return
    p_values = [item.p_value for item in associations]
    order = sorted(range(len(p_values)), key=lambda i: p_values[i])
    n = len(p_values)
    running = p_values[order[-1]]
    q_values: list[float] = [1.0] * n
    for rank in range(n, 0, -1):
        index = order[rank - 1]
        running = min(running, p_values[index] * n / rank)
        q_values[index] = running
    for item, q_value in zip(associations, q_values, strict=True):
        item.q_value = min(max(q_value, 0.0), 1.0)
        item.alpha = alpha


# ---------------------------------------------------------------------------
# Unsupervised: burden vs MIC
# ---------------------------------------------------------------------------


@dataclass
class BurdenCorrelation:
    """Spearman between resistance-mutation count and MIC, per drug."""

    drug: str
    n: int
    rho: float = 0.0
    p_value: float = 1.0
    n_carriers: int = 0

    @property
    def significant(self) -> bool:
        return self.p_value < DEFAULT_ALPHA and self.n >= 3

    def as_row(self) -> dict[str, Any]:
        return {
            "drug": self.drug,
            "n": self.n,
            "n_carriers": self.n_carriers,
            "rho": _round(self.rho),
            "p_value": _format_p(self.p_value),
            "significant": self.significant,
        }


def burden_correlation(
    carriers: Mapping[str, Mapping[str, bool]],
    records: Sequence[MicRecord],
    breakpoints: Mapping[tuple[str, str], Breakpoint],
    resistance_genes: Mapping[str, frozenset[str]] | None = None,
) -> list[BurdenCorrelation]:
    """Spearman between mutation burden and MIC, one test per drug.

    Uses the MIC concentration itself, on the log2 scale, so it does
    not depend on a breakpoint being correct -- a drug whose thresholds
    are disputed can still show whether genotype burden tracks the
    phenotype. ``carriers`` is accession -> {gene: bool}, and burden is
    the number of genes for which the genome carries a mutation.

    ``resistance_genes`` optionally narrows which genes count toward
    the burden. An empty mapping counts as *no filter*, not as "no
    gene qualifies": the caller derives it from a Fisher family that
    may have found nothing significant, and an empty result must not
    read as a cohort with no carriers at all.
    """
    # The associations test tells us which genes matter; feeding that
    # set in removes genes with no phenotypic relevance from the
    # burden, which would otherwise dilute it with noise.
    gene_filter: frozenset[str] | None = frozenset(resistance_genes or ())

    by_drug: dict[str, list[MicRecord]] = {}
    for record in records:
        if _lookup_breakpoint(breakpoints, record) is not None:
            by_drug.setdefault(record.drug.lower(), []).append(record)

    results: list[BurdenCorrelation] = []
    for drug, rows in sorted(by_drug.items()):
        burdens: list[float] = []
        log_mics: list[float] = []
        n_carriers = 0
        for record in rows:
            genes = [
                gene
                for gene, carrier in carriers.get(record.accession, {}).items()
                if carrier and (not gene_filter or gene in gene_filter)
            ]
            if genes:
                n_carriers += 1
            burdens.append(float(len(genes)))
            log_mics.append(record.log2_mic)

        if len(rows) < 2:
            results.append(BurdenCorrelation(drug=drug, n=len(rows), n_carriers=n_carriers))
            continue

        # Spearman needs variance in both variables; a cohort where
        # every genome carries the same genes is not a test, it is a
        # constant, and it must not be reported as one.
        if len(set(burdens)) < 2 or len(set(log_mics)) < 2:
            results.append(BurdenCorrelation(drug=drug, n=len(rows), n_carriers=n_carriers))
            continue

        rho, p_value = spearmanr(burdens, log_mics)
        results.append(
            BurdenCorrelation(
                drug=drug,
                n=len(rows),
                rho=float(rho),
                p_value=float(p_value),
                n_carriers=n_carriers,
            )
        )
    return results


# ---------------------------------------------------------------------------
# Supervised: concordance against curated rules
# ---------------------------------------------------------------------------


def concordance(
    carriers: Mapping[str, Mapping[str, bool]],
    records: Sequence[MicRecord],
    breakpoints: Mapping[tuple[str, str], Breakpoint],
    rules: PhenotypeRules,
) -> ConcordanceResult:
    """Score a curated rule set against the measured phenotype.

    For each (genome, drug) with a known phenotype, the rule set
    predicts R when any of the genome's mutated genes is a known
    resistance gene for that drug, and S otherwise. An I phenotype is
    never predicted, so an observed I counts as a minor error against
    a predicted S and a VME/ME only when the rule set called R.

    VME and ME are reported separately because they are not
    interchangeable: VME means a resistant infection would be treated
    as susceptible, which is the direction that matters clinically.
    """
    result = ConcordanceResult()
    per_drug: dict[str, dict[str, int]] = {}

    for record in records:
        phenotype = classify_mic(record.mic, _lookup_breakpoint(breakpoints, record))
        if phenotype == UNKNOWN:
            result.n_unknown_phenotype += 1
            continue

        mutated = {gene for gene, carrier in carriers.get(record.accession, {}).items() if carrier}
        predicted = (
            RESISTANT
            if any(rules.predicts_resistance(gene, record.drug.lower()) for gene in mutated)
            else SUSCEPTIBLE
        )

        drug = record.drug.lower()
        stats = per_drug.setdefault(
            drug,
            {"n": 0, "correct": 0, "vme": 0, "me": 0, "minor": 0},
        )
        stats["n"] += 1
        result.n_evaluable += 1

        if predicted == phenotype:
            stats["correct"] += 1
            continue

        if predicted == RESISTANT and phenotype == SUSCEPTIBLE:
            # Reported resistant, measured susceptible.
            result.n_major += 1
            stats["me"] += 1
        elif predicted == SUSCEPTIBLE and phenotype == RESISTANT:
            # Reported susceptible, measured resistant.
            result.n_very_major += 1
            stats["vme"] += 1
        else:
            # S vs I, or R vs I.
            result.n_minor += 1
            stats["minor"] += 1

    if result.n_evaluable:
        total = result.n_evaluable
        correct = total - result.n_major - result.n_very_major - result.n_minor
        result.agreement = correct / total
        result.very_major_error = result.n_very_major / total
        result.major_error = result.n_major / total
        result.minor_error = result.n_minor / total

    result.per_drug = {
        drug: {
            "n": stats["n"],
            "agreement": _round(stats["correct"] / stats["n"]) if stats["n"] else 0.0,
            "very_major_error": _round(stats["vme"] / stats["n"]) if stats["n"] else 0.0,
            "major_error": _round(stats["me"] / stats["n"]) if stats["n"] else 0.0,
            "minor_error": _round(stats["minor"] / stats["n"]) if stats["n"] else 0.0,
        }
        for drug, stats in sorted(per_drug.items())
    }
    return result


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------


def _round(value: float, precision: int = 4) -> float:
    if value is None or not math.isfinite(value):
        return 0.0
    return round(float(value), precision)


def _odds(value: float) -> str:
    """Render an odds ratio, keeping ``inf`` as what it is.

    A 2x2 table with a zero marginal gives an infinite odds ratio, and
    that is a real, computed result -- perfect separation, not a
    missing number. Rendering it as ``--`` would imply the metric
    could not be computed, which is a different claim.
    """
    if value is None:
        return "--"
    if math.isinf(value):
        return "inf"
    if not math.isfinite(value):
        return "--"
    return f"{value:.3f}"


def _format_p(value: float) -> str:
    if value is None or not math.isfinite(value):
        return "1.000"
    if 0 < value < 1e-4:
        return f"{value:.2e}"
    return f"{value:.4f}"


def _fmt(value: Any, precision: int = 3) -> str:
    if value is None:
        return "--"
    if isinstance(value, float) and not math.isfinite(value):
        return "--"
    if isinstance(value, float):
        return f"{value:.{precision}f}"
    return str(value)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------


def write_report(
    associations: Sequence[GeneDrugAssociation],
    burdens: Sequence[BurdenCorrelation],
    concordance_result: ConcordanceResult | None,
    output_dir: Path,
) -> dict[str, str]:
    """Write the phenotypic concordance report as Markdown and JSON.

    The three analyses are written separately, and the unsupervised
    ones come first: they are the evidence that MutationScan's calls
    carry phenotypic signal at all, before any curated table has had a
    chance to shape the answer.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("# Phenotypic concordance\n")
    lines.append(
        "Does a called mutation predict the phenotype the laboratory "
        "measures? Three measurements, in increasing order of how much "
        "they assume about the biology.\n"
    )

    lines.append("## Gene-drug association (unsupervised)\n")
    lines.append(
        "Carrying a mutation in a gene vs. being phenotypically "
        "resistant, Fisher exact two-sided, Benjamini-Hochberg across "
        "the gene x drug family. No curated knowledge used.\n"
    )
    lines.append("| gene | drug | carriers | resistant | both | n | odds ratio | p | q | sig |")
    lines.append("|------|------|----------|-----------|------|---|------------|---|---|-----|")
    if associations:
        for item in associations:
            lines.append(
                f"| {item.gene} | {item.drug} | {item.n_carrier} | {item.n_resistant} "
                f"| {item.n_both} | {item.n_evaluable} | {_odds(item.odds_ratio)} "
                f"| {_format_p(item.p_value)} | {_format_p(item.q_value)} "
                f"| {'**yes**' if item.significant else 'no'} |"
            )
    else:
        lines.append("| _none_ | | | | | | | | | |")
    n_significant = sum(1 for item in associations if item.significant)
    lines.append("")
    lines.append(
        f"**{n_significant} of {len(associations)} tests significant** at FDR < {DEFAULT_ALPHA}.\n"
    )

    lines.append("## Burden vs. MIC (unsupervised)\n")
    lines.append(
        "Spearman between the number of mutated genes a genome carries "
        "and its MIC on the log2 scale. Uses the concentration, not a "
        "dichotomised threshold, so it does not depend on the "
        "breakpoint being right.\n"
    )
    lines.append("| drug | n | carriers | rho | p | sig |")
    lines.append("|------|---|----------|-----|---|-----|")
    if burdens:
        for burden in burdens:
            lines.append(
                f"| {burden.drug} | {burden.n} | {burden.n_carriers} | {_fmt(burden.rho)} "
                f"| {_format_p(burden.p_value)} | {'**yes**' if burden.significant else 'no'} |"
            )
    else:
        lines.append("| _none_ | | | | | |")
    lines.append("")

    if concordance_result is not None:
        result = concordance_result
        lines.append("## Concordance vs. curated rules\n")
        lines.append(
            "VME = reported susceptible, measured resistant. ME = "
            "reported resistant, measured susceptible. They are not "
            "interchangeable.\n"
        )
        lines.append("| metric | value |")
        lines.append("|--------|-------|")
        lines.append(f"| evaluable | {result.n_evaluable} |")
        lines.append(f"| unknown phenotype | {result.n_unknown_phenotype} |")
        lines.append(f"| agreement | {_fmt(result.agreement)} |")
        lines.append(
            f"| very major error | {_fmt(result.very_major_error)} ({result.n_very_major}) |"
        )
        lines.append(f"| major error | {_fmt(result.major_error)} ({result.n_major}) |")
        lines.append(f"| minor error | {_fmt(result.minor_error)} ({result.n_minor}) |")
        lines.append("")
        if result.per_drug:
            lines.append("| drug | n | agreement | VME | ME | minor |")
            lines.append("|------|---|-----------|-----|----|-------|")
            for drug, stats in result.per_drug.items():
                lines.append(
                    f"| {drug} | {stats['n']} | {_fmt(stats['agreement'])} "
                    f"| {_fmt(stats['very_major_error'])} | {_fmt(stats['major_error'])} "
                    f"| {_fmt(stats['minor_error'])} |"
                )
            lines.append("")

    lines.append("## Data sources\n")
    lines.append(
        "- MIC: BV-BRC (formerly PATRIC) AMR phenotype service, or a "
        "curated table. `benchmarks/phenotypic.py` reads accession, "
        "drug, MIC and organism."
    )
    lines.append(
        "- Breakpoints: EUCAST `clinical_breakpoints` (open access) or "
        "CLSI M100 (paywalled). Both parse the same way."
    )
    lines.append("- Rules: a curated gene -> drug mapping, one row per resistance relationship.")

    markdown_path = output_dir / "phenotypic.md"
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    payload = {
        "gene_drug_associations": [item.as_row() for item in associations],
        "burden_correlations": [item.as_row() for item in burdens],
        "concordance": concordance_result.as_row() if concordance_result is not None else None,
    }
    json_path = output_dir / "phenotypic.json"
    json_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return {"markdown": str(markdown_path), "json": str(json_path)}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def read_carriers_from_genome_map(
    path: Path | str,
    mutation_column: str = "Mutated_Genes",
) -> dict[str, dict[str, bool]]:
    """Read gene-mutation carriage from a ``genome_mutation_map.csv``.

    This is the natural source: it is what MutationScan already writes,
    with one row per genome and a semicolon-joined ``Mutated_Genes``
    list. Feeding that straight in means the whole loop is
    ``mutationscan run ... && python -m benchmarks.phenotypic ...``
    with no intermediate format.

    A genome listed with no mutated genes comes back as an empty map,
    which scores as "carries nothing" rather than being dropped -- a
    genome the pipeline looked at and found clean is a result.
    """
    path = Path(path)
    dialect = _dialect_for(path)
    reader = csv.DictReader(_data_lines(path), dialect=dialect)
    if reader.fieldnames is None:
        raise SystemExit(f"{path}: empty genome map")
    columns = {_normalise_header(name): name for name in reader.fieldnames}
    accession_col = columns.get("accession")
    mutation_col = columns.get(mutation_column.lower())
    if accession_col is None or mutation_col is None:
        raise SystemExit(
            f"{path}: genome map needs 'accession' and '{mutation_column}' columns; "
            f"found {sorted(reader.fieldnames)}"
        )

    carriers: dict[str, dict[str, bool]] = {}
    for row in reader:
        accession = (row.get(accession_col) or "").strip()
        if not accession:
            continue
        genes = [gene.strip() for gene in (row.get(mutation_col) or "").split(";") if gene.strip()]
        carriers[accession] = dict.fromkeys(genes, True)
    return carriers


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mic", required=True, help="TSV/CSV of accession, drug, MIC")
    parser.add_argument(
        "--breakpoints",
        required=True,
        help="TSV/CSV of organism, drug, s_max, i_max, r_min",
    )
    parser.add_argument("--rules", help="optional curated gene -> drug resistance table")
    parser.add_argument(
        "--genome-map",
        help="genome_mutation_map.csv from a MutationScan run (preferred over --carriers)",
    )
    parser.add_argument(
        "--carriers",
        help="JSON mapping accession -> {gene: bool} (e.g. from genome_mutation_map.csv)",
    )
    parser.add_argument("--output", default="benchmarks/results/phenotypic")
    parser.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
    )

    if args.genome_map and args.carriers:
        parser.error("--genome-map and --carriers are alternative sources; pick one")

    records = read_mic_records(Path(args.mic))
    logger.info("read %d MIC records from %s", len(records), args.mic)

    breakpoints = read_breakpoints(Path(args.breakpoints))
    logger.info("read %d breakpoints", len(breakpoints))

    carriers: dict[str, dict[str, bool]] = {}
    if args.genome_map:
        carriers = read_carriers_from_genome_map(args.genome_map)
        logger.info("read carriers for %d accessions from %s", len(carriers), args.genome_map)
    elif args.carriers:
        with Path(args.carriers).open(encoding="utf-8") as handle:
            carriers = json.load(handle)
        logger.info("read carriers for %d accessions", len(carriers))

    associations = gene_drug_associations(carriers, records, breakpoints, alpha=args.alpha)
    n_significant = sum(1 for item in associations if item.significant)
    logger.info(
        "gene-drug association: %d tests, %d significant at FDR < %.2f",
        len(associations),
        n_significant,
        args.alpha,
    )

    # Feed the significant genes into the burden calculation, so the
    # burden counts phenotypically relevant genes rather than any gene.
    # An empty set means no gene was significant, which is *not* the
    # same as "filter every gene out" -- that would report zero
    # carriers for a cohort that plainly has some, and make the burden
    # column read as a null result.
    inverted: dict[str, set[str]] = {}
    for item in associations:
        if item.significant:
            inverted.setdefault(item.gene, set()).add(item.drug)
    burdens = burden_correlation(
        carriers,
        records,
        breakpoints,
        resistance_genes={gene: frozenset(drugs) for gene, drugs in inverted.items()},
    )

    concordance_result = None
    if args.rules:
        rules = read_rules(Path(args.rules))
        concordance_result = concordance(carriers, records, breakpoints, rules)
        logger.info(
            "concordance: agreement=%.3f VME=%.3f ME=%.3f",
            concordance_result.agreement,
            concordance_result.very_major_error,
            concordance_result.major_error,
        )

    paths = write_report(associations, burdens, concordance_result, Path(args.output))
    print(f"report: {paths['markdown']}")
    print(f"report: {paths['json']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
