"""Tests for run configuration parsing and validation."""

from __future__ import annotations

import pytest

from mutation_scan.config import Config, ConfigError, load_config


def test_defaults_are_usable(tmp_path):
    config = Config.from_mapping({}, base_dir=tmp_path)

    assert config.job_name == "default_run"
    assert config.min_identity_percent == 80.0
    assert config.min_coverage_percent == 80.0
    assert config.cooccurrence.levels == ("gene",)
    assert config.cooccurrence.min_count == 1
    assert config.targets == ()


def test_relative_paths_resolve_against_base_dir(tmp_path):
    config = Config.from_mapping(
        {"genomes_dir": "genomes", "references_dir": "refs", "output_root": "out"},
        base_dir=tmp_path,
    )

    assert config.genomes_dir == tmp_path / "genomes"
    assert config.references_dir == tmp_path / "refs"
    assert config.output_root == tmp_path / "out"


def test_absolute_paths_are_left_alone(tmp_path):
    absolute = (tmp_path / "elsewhere").resolve()
    config = Config.from_mapping({"genomes_dir": str(absolute)}, base_dir=tmp_path)

    assert config.genomes_dir == absolute


def test_run_layout_is_a_single_folder(tmp_path):
    config = Config.from_mapping({"job_name": "my_run", "output_root": "out"}, base_dir=tmp_path)
    run_dir = tmp_path / "out" / "my_run"

    assert config.run_dir == run_dir
    assert config.proteins_dir == run_dir / "proteins"
    assert config.refs_dir == run_dir / "refs"
    assert config.mutations_csv == run_dir / "mutations.csv"
    # Every declared output lives inside the run folder.
    for path in config.outputs():
        assert run_dir in path.parents


def test_variant_level_adds_its_output(tmp_path):
    gene_only = Config.from_mapping({}, base_dir=tmp_path)
    assert gene_only.cooccurrence_variants_csv not in gene_only.outputs()

    both = Config.from_mapping({"cooccurrence": {"levels": ["variant", "gene"]}}, base_dir=tmp_path)
    # Order is normalized regardless of how the user wrote it.
    assert both.cooccurrence.levels == ("gene", "variant")
    assert both.cooccurrence_variants_csv in both.outputs()


def test_targets_accept_list_or_string(tmp_path):
    from_list = Config.from_mapping({"targets": ["geneA", " geneB ", "geneA"]}, base_dir=tmp_path)
    assert from_list.targets == ("geneA", "geneB")

    from_string = Config.from_mapping({"targets": "geneA, geneB"}, base_dir=tmp_path)
    assert from_string.targets == ("geneA", "geneB")


def test_null_values_fall_back_to_defaults(tmp_path):
    # A YAML key written as `uniprot_taxid:` parses to None.
    config = Config.from_mapping({"uniprot_taxid": None, "targets": None}, base_dir=tmp_path)

    assert config.uniprot_taxid == ""
    assert config.targets == ()


@pytest.mark.parametrize(
    "mapping",
    [
        {"min_identity_percent": 101},
        {"min_identity_percent": -1},
        {"min_coverage_percent": "abc"},
        {"threads": 0},
        {"cooccurrence": {"min_count": 0}},
        {"cooccurrence": {"levels": []}},
        {"cooccurrence": {"levels": ["genome"]}},
        {"job_name": ""},
        {"job_name": "nested/run"},
        {"genomes_dir": ""},
        {"min_identity_pct": 90},  # typo in a key name must not pass silently
        {"cooccurrence": {"level": ["gene"]}},
    ],
)
def test_invalid_configs_are_rejected(mapping, tmp_path):
    with pytest.raises(ConfigError):
        Config.from_mapping(mapping, base_dir=tmp_path)


def test_require_inputs_reports_missing_directories(tmp_path):
    config = Config.from_mapping({}, base_dir=tmp_path)

    with pytest.raises(ConfigError, match="genomes_dir"):
        config.require_inputs()

    config.genomes_dir.mkdir(parents=True)
    with pytest.raises(ConfigError, match="references_dir"):
        config.require_inputs()

    config.references_dir.mkdir(parents=True)
    config.require_inputs()  # both present -> no error


def test_require_inputs_rejects_a_file_where_a_directory_is_expected(tmp_path):
    config = Config.from_mapping({}, base_dir=tmp_path)
    config.genomes_dir.parent.mkdir(parents=True, exist_ok=True)
    config.genomes_dir.write_text("not a directory", encoding="utf-8")

    with pytest.raises(ConfigError, match="not a directory"):
        config.require_inputs()


def test_load_config_reads_yaml(tmp_path):
    config_file = tmp_path / "config.yaml"
    config_file.write_text(
        "job_name: yaml_run\nmin_identity_percent: 95\ncooccurrence:\n  min_count: 3\n",
        encoding="utf-8",
    )

    config = load_config(config_file, base_dir=tmp_path)

    assert config.job_name == "yaml_run"
    assert config.min_identity_percent == 95.0
    assert config.cooccurrence.min_count == 3
    # Unspecified nested keys keep their defaults.
    assert config.cooccurrence.levels == ("gene",)


def test_load_config_errors_are_actionable(tmp_path):
    missing = tmp_path / "nope.yaml"
    with pytest.raises(ConfigError, match="not found"):
        load_config(missing)

    broken = tmp_path / "broken.yaml"
    broken.write_text("job_name: [unclosed\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="Invalid YAML"):
        load_config(broken, base_dir=tmp_path)

    not_a_mapping = tmp_path / "list.yaml"
    not_a_mapping.write_text("- a\n- b\n", encoding="utf-8")
    with pytest.raises(ConfigError, match="must be a mapping"):
        load_config(not_a_mapping, base_dir=tmp_path)


def test_empty_yaml_yields_defaults(tmp_path):
    empty = tmp_path / "empty.yaml"
    empty.write_text("# nothing here\n", encoding="utf-8")

    assert load_config(empty, base_dir=tmp_path).job_name == "default_run"


def test_shipped_config_is_valid():
    """The config we ship must parse -- the old pipeline shipped a broken one."""
    from pathlib import Path

    repo_root = Path(__file__).resolve().parents[2]
    config = load_config(repo_root / "config" / "config.yaml", base_dir=repo_root)

    assert config.job_name
    assert config.min_identity_percent > 0
