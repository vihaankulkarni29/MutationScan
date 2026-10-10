"""Benchmark datasets: download real cohorts, or synthesise one.

Three sources, in increasing order of how well we know
the answer:

* **PATRIC** and **NCTC 3000** -- real genomes with
  curated resistance annotations. Downloading them
  needs network access and an API key (PATRIC), so the
  functions here are explicit about that and cache
  whatever they fetch.
* **Synthetic spike-in** -- genomes built from the
  project's own reference proteins with *known*
  mutations, *known* co-mutation structure and *known*
  silent changes. This is the only source where
  sensitivity, specificity and dN/dS all have exact
  ground truth, so it is the default for CI.

The synthetic generator is the interesting half. It
back-translates each reference protein to a CDS,
introduces a chosen mix of non-synonymous substitutions
(the mutations under test) and synonymous ones (which
give dS something to be non-zero), and records every
change it made. Because the genome *is* the CDS, tblastn
finds it at coordinates 1..N and the pipeline's own
dN/dS stage reads back exactly the sequence that was
written -- so a disagreement between the pipeline's dN/dS
and the recorded ground truth is a real bug, not a
coordinate artefact.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from mutation_scan.evolution import (  # noqa: E402
    GENETIC_CODE,
    back_translate,
    codon_sites,
)

__all__ = [
    "generate_synthetic",
    "download_patric",
    "download_nctc",
]

#: Amino acids that a substitution may target, weighted
#: toward chemically distinct changes so the spike-in
#: exercises the severity scorer.
_SUBSTITUTABLE = "ACDEFGHIKLMNPQRSTVWY"


def _codons_for(aa: str) -> list[str]:
    return [c for c, a in GENETIC_CODE.items() if a == aa]


def _nearby_codon(ref_codon: str, alt_aa: str, rng: random.Random) -> str:
    """A codon for ``alt_aa`` differing from ``ref_codon`` by one base
    where possible, so the spike-in looks like a point mutation."""
    candidates = _codons_for(alt_aa)
    if not candidates:
        return ref_codon
    one_off = [c for c in candidates if sum(a != b for a, b in zip(c, ref_codon, strict=True)) == 1]
    pool = one_off or candidates
    return rng.choice(pool)


def _synonymous_codon(codon: str, rng: random.Random) -> str | None:
    """A different codon encoding the same amino acid, if one exists."""
    aa = GENETIC_CODE.get(codon)
    if aa is None:
        return None
    alternatives = [c for c in _codons_for(aa) if c != codon]
    return rng.choice(alternatives) if alternatives else None


def _count_sites(reference_cds: str) -> tuple[float, float]:
    n_sites = 0.0
    s_sites = 0.0
    for i in range(0, len(reference_cds), 3):
        n, s = codon_sites(reference_cds[i : i + 3])
        n_sites += n
        s_sites += s
    return n_sites, s_sites


def _differences(ref_cds: str, query_cds: str) -> tuple[int, int]:
    """Non-synonymous and synonymous nucleotide differences."""
    n_diffs = 0
    s_diffs = 0
    for i in range(0, min(len(ref_cds), len(query_cds)), 3):
        ref_codon = ref_cds[i : i + 3]
        query_codon = query_cds[i : i + 3]
        if len(ref_codon) < 3 or len(query_codon) < 3:
            continue
        ref_aa = GENETIC_CODE.get(ref_codon)
        query_aa = GENETIC_CODE.get(query_codon)
        if ref_aa is None or query_aa is None:
            continue
        diff = sum(a != b for a, b in zip(ref_codon, query_codon, strict=True))
        if diff == 0:
            continue
        if ref_aa == query_aa:
            s_diffs += diff
        else:
            n_diffs += diff
    return n_diffs, s_diffs


def generate_synthetic(
    references_dir: Path | str,
    output_dir: Path | str,
    n_genomes: int = 50,
    seed: int = 42,
    mutation_rate: float = 0.4,
    epistatic_pairs: Sequence[tuple[str, str]] | None = None,
) -> dict[str, Any]:
    """Build a synthetic cohort with exact ground truth.

    Args:
        references_dir: directory of ``<gene>_WT.faa`` reference proteins.
        output_dir: where genomes, the lineage manifest and the ground
            truth JSON are written.
        n_genomes: cohort size.
        seed: RNG seed, so a run is reproducible.
        mutation_rate: probability a gene is mutated in a genome.
        epistatic_pairs: gene pairs to co-mutate, defining the
            ground-truth epistasis. Defaults to a few pairs drawn
            from the reference set.

    Returns:
        The ground-truth dictionary, also written to
        ``<output_dir>/ground_truth.json``.
    """
    references_dir = Path(references_dir)
    output_dir = Path(output_dir)
    rng = random.Random(seed)

    genes: dict[str, str] = {}
    for path in sorted(references_dir.glob("*_WT.faa")):
        gene = path.name[: -len("_WT.faa")]
        genes[gene] = _read_fasta(path)
    if not genes:
        raise SystemExit(f"no *_WT.faa references found in {references_dir}")

    gene_names = sorted(genes)
    if epistatic_pairs is None:
        epistatic_pairs = _default_epistatic_pairs(gene_names, rng)

    genomes_dir = output_dir / "genomes"
    # A previous run may have left genomes behind (a smaller
    # cohort reusing this directory would otherwise be scored
    # against stale assemblies that are not in the ground truth).
    if genomes_dir.is_dir():
        for stale in genomes_dir.glob("*.fna"):
            stale.unlink()
    genomes_dir.mkdir(parents=True, exist_ok=True)

    mutations: set[tuple[str, str, str]] = set()
    lineages: dict[str, str] = {}
    genome_cds: dict[str, dict[str, str]] = {}
    co_mutated: set[frozenset[str]] = set()

    for index in range(n_genomes):
        accession = f"SIM_{index:04d}"
        lineage = f"ST{rng.randint(1, 20)}"
        lineages[accession] = lineage

        mutated_genes: set[str] = set()
        cds_by_gene: dict[str, str] = {}

        for gene in gene_names:
            reference_protein = genes[gene]
            reference_cds = back_translate(reference_protein)
            codons = [reference_cds[i : i + 3] for i in range(0, len(reference_cds), 3)]

            is_epistatic = any(gene in pair for pair in epistatic_pairs)
            will_mutate = rng.random() < mutation_rate
            if is_epistatic:
                will_mutate = will_mutate or rng.random() < 0.6

            if not will_mutate:
                cds_by_gene[gene] = reference_cds
                continue

            mutated_genes.add(gene)
            n_changes = rng.randint(1, 3)
            positions = sorted(rng.sample(range(len(codons)), min(n_changes, len(codons))))
            for position in positions:
                ref_codon = codons[position]
                ref_aa = GENETIC_CODE.get(ref_codon, "")
                alt_aa = rng.choice([a for a in _SUBSTITUTABLE if a != ref_aa])
                alt_codon = _nearby_codon(ref_codon, alt_aa, rng)
                codons[position] = alt_codon
                mutations.add(
                    (
                        accession,
                        gene,
                        f"{ref_aa}{position + 1}{alt_aa}",
                    )
                )

            # A silent change or two, so dS is non-zero.
            for _ in range(rng.randint(0, 2)):
                silent = rng.randrange(len(codons))
                alt = _synonymous_codon(codons[silent], rng)
                if alt:
                    codons[silent] = alt

            cds_by_gene[gene] = "".join(codons)

        # Record co-mutation of the defined epistatic pairs.
        for gene_a, gene_b in epistatic_pairs:
            if gene_a in mutated_genes and gene_b in mutated_genes:
                co_mutated.add(frozenset((gene_a, gene_b)))

        genome_cds[accession] = cds_by_gene
        _write_genome(genomes_dir / f"{accession}.fna", cds_by_gene)

    # Ground-truth dN/dS from the known sequences.
    dn_ds: dict[str, float] = {}
    for gene in gene_names:
        reference_cds = back_translate(genes[gene])
        n_sites, s_sites = _count_sites(reference_cds)
        n_total = 0
        s_total = 0
        for cds_by_gene in genome_cds.values():
            cds = cds_by_gene.get(gene)
            if not cds:
                continue
            n_diffs, s_diffs = _differences(reference_cds, cds)
            n_total += n_diffs
            s_total += s_diffs
        dn = n_total / n_sites if n_sites else 0.0
        ds = s_total / s_sites if s_sites else 0.0
        dn_ds[gene] = dn / ds if ds else float("inf")

    # Ground-truth epistasis: pairs that co-occur more often than
    # chance. The synthetic design makes the defined pairs co-mutate,
    # so they are the ground truth.
    epistasis: set[frozenset[str]] = set()
    for gene_a, gene_b in epistatic_pairs:
        epistasis.add(frozenset((gene_a, gene_b)))

    _write_lineage_manifest(output_dir / "lineage_manifest.tsv", lineages)

    ground_truth = {
        "mutations": sorted([list(m) for m in mutations], key=lambda x: (x[0], x[1], x[2])),
        "epistasis": sorted([sorted(p) for p in epistasis], key=lambda x: (x[0], x[1])),
        "dn_ds": {
            gene: (value if value != float("inf") else "inf") for gene, value in dn_ds.items()
        },
        "lineages": lineages,
        "known_mutations": [],
        "negatives": [],
    }
    with (output_dir / "ground_truth.json").open("w", encoding="utf-8") as handle:
        json.dump(ground_truth, handle, indent=2)

    return ground_truth


def _default_epistatic_pairs(gene_names: list[str], rng: random.Random) -> list[tuple[str, str]]:
    if len(gene_names) < 2:
        return []
    pairs = list(zip(gene_names, gene_names[1:], strict=False))
    return pairs[: min(3, len(pairs))]


def _read_fasta(path: Path) -> str:
    from mutation_scan.fasta import read_protein_sequence

    record = read_protein_sequence(path, fallback_id=path.stem)
    return str(record.seq)


def _write_genome(path: Path, cds_by_gene: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for gene, cds in cds_by_gene.items():
            handle.write(f">{gene}\n")
            for i in range(0, len(cds), 60):
                handle.write(cds[i : i + 60] + "\n")


def _write_lineage_manifest(path: Path, lineages: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        handle.write("Accession\tLineage\n")
        for accession, lineage in sorted(lineages.items()):
            handle.write(f"{accession}\t{lineage}\n")


def download_patric(output_dir: Path | str, **kwargs: Any) -> Path:
    """Download a PATRIC AMR cohort.

    PATRIC's API needs a key and serves large files, so this
    is a stub that records the intent and points at the
    manual download. It raises rather than silently producing
    an empty directory that a benchmark would score as
    "found nothing".
    """
    raise NotImplementedError(
        "PATRIC download requires an API key and manual curation; "
        "see the benchmarking plan. Use generate_synthetic for CI."
    )


def download_nctc(output_dir: Path | str, **kwargs: Any) -> Path:
    """Download the NCTC 3000 cohort.

    The NCTC 3000 assemblies are served by the Wellcome Sanger
    Institute and are large; this is a stub for the same reason
    as :func:`download_patric`.
    """
    raise NotImplementedError(
        "NCTC 3000 download is a large manual fetch; "
        "see the benchmarking plan. Use generate_synthetic for CI."
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--references",
        default="data/references",
        help="directory of <gene>_WT.faa reference proteins",
    )
    parser.add_argument(
        "--output",
        default="benchmarks/data/synthetic",
        help="output directory",
    )
    parser.add_argument("--genomes", type=int, default=50)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args(argv)

    truth = generate_synthetic(args.references, args.output, n_genomes=args.genomes, seed=args.seed)
    print(
        f"wrote {args.output}: {len(truth['mutations'])} mutations, "
        f"{len(truth['epistasis'])} epistatic pairs, "
        f"{len(truth['dn_ds'])} genes with dN/dS"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
