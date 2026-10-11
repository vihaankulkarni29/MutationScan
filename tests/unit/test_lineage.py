"""Tests for lineage tagging (population-structure control)."""

from __future__ import annotations

import pandas as pd

from mutation_scan.lineage import (
    annotate_epistasis,
    annotate_genome_map,
    classify_clonality,
    lineage_distribution,
    read_lineage_manifest,
)


def _manifest(accession_lineage: list[tuple[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Accession": [a for a, _ in accession_lineage],
            "Lineage": [lin for _, lin in accession_lineage],
        }
    )


def test_read_lineage_manifest_empty() -> None:
    assert read_lineage_manifest(None) == {}
    assert read_lineage_manifest("nonexistent.tsv") == {}


def test_read_lineage_manifest_valid(tmp_path) -> None:
    p = tmp_path / "lineage.tsv"
    p.write_text("Accession\tLineage\nG1\tST131\nG2\tST11\n")
    mapping = read_lineage_manifest(p)
    assert mapping == {"G1": "ST131", "G2": "ST11"}


def test_read_lineage_manifest_case_insensitive_columns(tmp_path) -> None:
    p = tmp_path / "lineage.csv"
    p.write_text("ACCESSION,LINEAGE\nG1,ST131\nG2,ST11\n")
    mapping = read_lineage_manifest(p)
    assert mapping == {"G1": "ST131", "G2": "ST11"}


def test_read_lineage_manifest_whitespace_stripped(tmp_path) -> None:
    p = tmp_path / "lineage.tsv"
    p.write_text("Accession\tLineage\n G1 \t ST131 \n")
    mapping = read_lineage_manifest(p)
    assert mapping == {"G1": "ST131"}


def test_lineage_distribution_empty_map() -> None:
    assert lineage_distribution(["G1", "G2"], {}) == {}


def test_lineage_distribution_counts() -> None:
    m = {"G1": "ST131", "G2": "ST131", "G3": "ST11"}
    assert lineage_distribution(["G1", "G2", "G3"], m) == {"ST131": 2, "ST11": 1}


def test_lineage_distribution_unassigned() -> None:
    m = {"G1": "ST131"}
    assert lineage_distribution(["G1", "G2"], m) == {"ST131": 1, "unassigned": 1}


def test_classify_clonality_empty_map() -> None:
    assert classify_clonality(["G1"], {}) == "Unknown"
    assert classify_clonality([], {}) == "Unknown"


def test_classify_clonality_clonal() -> None:
    m = {"G1": "ST131", "G2": "ST131"}
    assert classify_clonality(["G1", "G2"], m) == "Clonal"


def test_classify_clonality_convergent() -> None:
    m = {"G1": "ST131", "G2": "ST11"}
    assert classify_clonality(["G1", "G2"], m) == "Convergent"


def test_classify_clonality_mixed_unassigned() -> None:
    m = {"G1": "ST131", "G2": "ST11"}
    assert classify_clonality(["G1", "G2", "G3"], m) == "Convergent"
    assert classify_clonality(["G3", "G4"], m) == "Unknown"


def test_annotate_genome_map() -> None:
    gm = pd.DataFrame(
        {
            "Accession": ["G1", "G2"],
            "N_Mutated_Genes": [1, 2],
            "Mutated_Genes": ["acrB", "acrB;marR"],
            "N_Mutations": [1, 2],
            "N_LoF": [0, 0],
            "Mutation_Map": ["acrB:S83L", "acrB:S83L | marR:LoF"],
        }
    )
    m = {"G1": "ST131", "G2": "ST11"}
    out = annotate_genome_map(gm, m)
    assert list(out["Lineage"]) == ["ST131", "ST11"]


def test_annotate_genome_map_empty_map() -> None:
    gm = pd.DataFrame({"Accession": ["G1"], "N_Mutated_Genes": [1]})
    out = annotate_genome_map(gm, {})
    assert "Lineage" not in out.columns or out["Lineage"].isna().all()


def test_annotate_epistasis_no_lineage() -> None:
    nets = pd.DataFrame({"Variant_A": ["acrB:S83L"], "Variant_B": ["marR:LoF"]})
    muts = pd.DataFrame(
        {
            "Accession": ["G1", "G2"],
            "Gene": ["acrB", "marR"],
            "Mutation": ["S83L", "LoF"],
        }
    )
    out = annotate_epistasis(nets, muts, {})
    assert "Lineage_Distribution" not in out.columns


def test_annotate_epistasis_clonal() -> None:
    nets = pd.DataFrame({"Variant_A": ["acrB:S83L"], "Variant_B": ["marR:LoF"]})
    muts = pd.DataFrame(
        {
            "Accession": ["G1", "G1", "G2", "G2"],
            "Gene": ["acrB", "marR", "acrB", "marR"],
            "Mutation": ["S83L", "LoF", "S83L", "LoF"],
        }
    )
    m = {"G1": "ST131", "G2": "ST131"}
    out = annotate_epistasis(nets, muts, m)
    assert out.loc[0, "Lineage_Distribution"] == "ST131:2"
    assert out.loc[0, "Clonality"] == "Clonal"


def test_annotate_epistasis_convergent() -> None:
    nets = pd.DataFrame({"Variant_A": ["acrB:S83L"], "Variant_B": ["marR:LoF"]})
    muts = pd.DataFrame(
        {
            "Accession": ["G1", "G1", "G2", "G2"],
            "Gene": ["acrB", "marR", "acrB", "marR"],
            "Mutation": ["S83L", "LoF", "S83L", "LoF"],
        }
    )
    m = {"G1": "ST131", "G2": "ST11"}
    out = annotate_epistasis(nets, muts, m)
    assert out.loc[0, "Clonality"] == "Convergent"


def test_annotate_epistasis_mixed_unassigned() -> None:
    nets = pd.DataFrame({"Variant_A": ["acrB:S83L"], "Variant_B": ["marR:LoF"]})
    muts = pd.DataFrame(
        {
            "Accession": ["G1", "G2", "G3"],
            "Gene": ["acrB", "acrB", "marR"],
            "Mutation": ["S83L", "S83L", "LoF"],
        }
    )
    m = {"G1": "ST131", "G2": "ST131"}  # G3 unassigned
    out = annotate_epistasis(nets, muts, m)
    # Both G1 and G2 carry acrB:S83L; only G3 carries marR:LoF -> intersection empty
    assert out.loc[0, "Clonality"] == "Unknown"
