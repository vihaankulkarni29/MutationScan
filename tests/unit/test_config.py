"""Tests for run configuration parsing and validation."""

from __future__ import annotations

import pytest

from mutation_scan.config import _DEFAULTS, Config, ConfigError, load_config


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


# -- nested defaults must not leak between builds ----------------------------


def test_epistasis_settings_do_not_leak_into_later_configs(tmp_path):
    """Regression: `dict(_DEFAULTS)` is a shallow copy.

    ``raw["epistasis"]`` used to alias the module-level defaults, so building a
    config with ``epistasis.threads: 8`` rewrote the default for the rest of the
    process. Every later Config then silently inherited it -- the same class of
    cross-call state leak as the pool initializer taking no arguments.
    """
    tuned = Config.from_mapping({"epistasis": {"threads": 8, "min_count": 4}}, tmp_path)
    assert tuned.epistasis.threads == 8
    assert tuned.epistasis.min_count == 4

    plain = Config.from_mapping({}, tmp_path)

    assert plain.epistasis.threads == 1
    assert plain.epistasis.min_count == 1
    assert _DEFAULTS["epistasis"]["threads"] is None
    assert _DEFAULTS["epistasis"]["min_count"] == 1


def test_epistasis_threads_inherit_the_run_width_by_default(tmp_path):
    """One inheritance rule, in one place, so both entry points agree.

    The file default used to be a hard 1: `--threads 8` and `snakemake --cores 8`
    both widened extraction and variant calling while epistasis silently stayed
    serial, and nothing in the run's output said so.
    """
    wide = Config.from_mapping({"threads": 6}, tmp_path)
    assert wide.epistasis.threads == 6

    narrow = Config.from_mapping({"threads": 2}, tmp_path)
    assert narrow.epistasis.threads == 2


def test_explicit_epistasis_threads_still_win(tmp_path):
    config = Config.from_mapping({"threads": 6, "epistasis": {"threads": 2}}, tmp_path)

    assert config.threads == 6
    assert config.epistasis.threads == 2


def test_cooccurrence_settings_do_not_leak_into_later_configs(tmp_path):
    tuned = Config.from_mapping({"cooccurrence": {"min_count": 9}}, tmp_path)
    assert tuned.cooccurrence.min_count == 9

    plain = Config.from_mapping({}, tmp_path)

    assert plain.cooccurrence.min_count == 1
    assert _DEFAULTS["cooccurrence"]["min_count"] == 1


# -- CLI thread wiring -------------------------------------------------------


def _run_config(tmp_path, *extra):
    from mutation_scan.cli import _config_from_args, build_parser

    genomes = tmp_path / "genomes"
    references = tmp_path / "refs"
    genomes.mkdir(exist_ok=True)
    references.mkdir(exist_ok=True)
    args = build_parser().parse_args(
        [
            "run",
            "--genomes",
            str(genomes),
            "--references",
            str(references),
            "--out",
            str(tmp_path / "out"),
            *extra,
        ]
    )
    return _config_from_args(args)


def test_epistasis_threads_default_to_the_thread_count(tmp_path):
    """`--threads 8` must widen epistasis too, or the flag quietly does nothing.

    It used to leave epistasis at its config default of 1, with no command-line
    way to change it.
    """
    config = _run_config(tmp_path, "--threads", "8")

    assert config.threads == 8
    assert config.epistasis.threads == 8


def test_explicit_epistasis_threads_win(tmp_path):
    config = _run_config(tmp_path, "--threads", "8", "--epistasis-threads", "2")

    assert config.threads == 8
    assert config.epistasis.threads == 2


def test_omitting_both_leaves_epistasis_at_its_default(tmp_path):
    config = _run_config(tmp_path)

    assert config.threads == 1
    assert config.epistasis.threads == 1


def test_epistasis_threads_do_not_disturb_the_other_settings(tmp_path):
    """Widening the stage must not clobber fdr_threshold or min_count."""
    tuned = _run_config(tmp_path, "--threads", "4")
    assert tuned.epistasis.fdr_threshold == 0.05
    assert tuned.epistasis.min_count == 1


# -- the README must document a config the tool accepts ----------------------


def test_the_readme_config_block_still_loads(tmp_path):
    """Keep the documented schema honest.

    The README used to show `epistasis.threads: 1`, which parsed fine and was
    quietly wrong: it documented the exact value that left the stage serial
    while the rest of the run went parallel. A doc block that the loader rejects
    would be obvious; one that loads but means something stale is not, which is
    why this also asserts the *resolved* values rather than just parseability.
    """
    import re
    from pathlib import Path

    from mutation_scan.config import load_config

    readme = Path(__file__).resolve().parents[2] / "README.md"
    blocks = re.findall(r"```yaml\n(.*?)```", readme.read_text(encoding="utf-8"), re.S)
    assert blocks, "no yaml config block found in the README"

    block = tmp_path / "readme.yaml"
    block.write_text(blocks[0], encoding="utf-8")

    config = load_config(block, base_dir=tmp_path)

    # The documented defaults must all resolve to the same width, or the README
    # is telling a user to run a narrower stage than it implies.
    assert config.epistasis.threads == config.threads
    assert config.epistasis.fdr_threshold > 0
    assert config.epistasis.min_count >= 1
