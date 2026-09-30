"""Tests for the single source of truth about a run's parallelism.

The point of :mod:`mutation_scan.parallel_plan` is that one function answers
"how many workers will this stage actually use", and the run summary, the
`config-check` report and the workflow all read that same answer. These tests
pin the properties that make it worth trusting.
"""

from __future__ import annotations

import pytest

from mutation_scan import epistasis
from mutation_scan.config import Config
from mutation_scan.parallel import _THREAD_CEILING, available_cores, resolve_workers
from mutation_scan.parallel_plan import (
    _EPIGASTIS_SERIAL_PAIR_LIMIT,
    Stage,
    describe,
    format_plan,
)

STAGES = ("extraction", "variant calling", "co-occurrence", "epistasis")


def _config(tmp_path, **overrides):
    return Config.from_mapping(overrides, tmp_path)


class TestPlanCoversEveryStage:
    def test_every_stage_of_a_run_appears(self, tmp_path):
        names = {stage.name for stage in describe(_config(tmp_path))}

        assert names == set(STAGES)

    def test_the_plan_is_in_run_order(self, tmp_path):
        assert [stage.name for stage in describe(_config(tmp_path))] == list(STAGES)


class TestRequestedVersusEffective:
    def test_a_run_asks_for_n_and_gets_n(self, tmp_path):
        config = _config(tmp_path, threads=4, variant_threads=3)

        plan = {stage.name: stage for stage in describe(config)}

        assert plan["extraction"].requested == 4
        assert plan["extraction"].effective == 4
        assert plan["variant calling"].requested == 3
        assert plan["variant calling"].effective == 3

    def test_a_clamping_stage_reports_the_reduction(self, tmp_path):
        """Asking for more than the ceiling must be visible, not silent.

        A thread stage is hard-capped to stay inside the file-descriptor budget.
        A user who asks for 500 and gets 64 has still been given less than they
        asked for, and the report has to say so.
        """
        config = _config(tmp_path, threads=500, variant_threads=1)

        plan = {stage.name: stage for stage in describe(config)}

        assert plan["extraction"].requested == 500
        assert plan["extraction"].effective == _THREAD_CEILING
        assert plan["extraction"].clamped
        # Both can be true at once: the cap of 64 is itself above the core count
        # on most machines, and that is worth reporting alongside the clamp.
        if _THREAD_CEILING > available_cores():
            assert plan["extraction"].oversubscribed

    def test_oversubscription_is_flagged_but_not_clamped(self, tmp_path):
        """A knowingly-wide request is honoured, and reported as a mistake.

        Deliberately distinct from `clamped`: this stage runs at full width and
        is merely slower, which is a different problem from one that quietly ran
        at a fraction of the width asked for.
        """
        config = _config(tmp_path, threads=available_cores() + 4, variant_threads=1)

        plan = {stage.name: stage for stage in describe(config)}

        assert plan["extraction"].effective == available_cores() + 4
        assert plan["extraction"].oversubscribed
        assert not plan["extraction"].clamped

    def test_a_matching_request_is_neither_clamped_nor_oversubscribed(self, tmp_path):
        config = _config(tmp_path, threads=2, variant_threads=2)

        plan = {stage.name: stage for stage in describe(config)}

        assert not plan["extraction"].clamped
        assert not plan["extraction"].oversubscribed

    def test_serial_stages_are_marked_not_chosen_by_the_user(self, tmp_path):
        plan = {stage.name: stage for stage in describe(_config(tmp_path))}

        assert plan["co-occurrence"].not_parallelised


class TestEpistasisThresholdCannotDrift:
    def test_the_duplicated_threshold_matches_the_real_one(self):
        """`parallel_plan` deliberately does not import epistasis.

        Importing it would pull scipy into every `config-check`, which has to
        stay fast. The cost of that choice is a copied constant, so this test
        fails the moment the two disagree -- which would make the plan promise a
        serial fallback the stage does not actually have.
        """
        assert _EPIGASTIS_SERIAL_PAIR_LIMIT == epistasis._SERIAL_PAIR_LIMIT

    def test_the_epistasis_stage_advertises_its_serial_fallback(self, tmp_path):
        plan = {stage.name: stage for stage in describe(_config(tmp_path))}

        assert plan["epistasis"].serial_below == epistasis._SERIAL_PAIR_LIMIT


class TestRendering:
    def test_plan_renders_one_line_per_stage(self, tmp_path):
        lines = format_plan(describe(_config(tmp_path))).splitlines()

        assert len(lines) == len(STAGES)

    def test_a_clamped_stage_is_flagged_in_the_text(self, tmp_path):
        config = _config(tmp_path, threads=500, variant_threads=1)

        text = format_plan(describe(config))

        assert "clamped" in text

    def test_an_oversubscribed_stage_is_flagged_in_the_text(self, tmp_path):
        config = _config(tmp_path, threads=available_cores() + 4, variant_threads=1)

        text = format_plan(describe(config))

        assert "oversubscribed" in text

    def test_a_healthy_plan_is_not_flagged(self, tmp_path):
        config = _config(tmp_path, threads=2, variant_threads=2)

        text = format_plan(describe(config))

        assert "clamped" not in text
        assert "oversubscribed" not in text

    def test_the_serial_fallback_threshold_is_visible(self, tmp_path):
        text = format_plan(describe(_config(tmp_path)))

        assert str(_EPIGASTIS_SERIAL_PAIR_LIMIT) in text


class TestSerialisation:
    def test_a_stage_round_trips_through_to_dict(self):
        stage = Stage("x", "process", 8, 4, serial_below=99)

        payload = stage.to_dict()

        assert payload["stage"] == "x"
        assert payload["workload"] == "process"
        assert payload["requested_workers"] == 8
        assert payload["effective_workers"] == 4
        assert payload["clamped"] is True
        assert payload["serial_below_items"] == 99

    def test_optional_fields_are_omitted_when_unset(self):
        payload = Stage("x", "thread", 1, 1).to_dict()

        assert "serial_below_items" not in payload
        assert "not_parallelised" not in payload

    @pytest.mark.parametrize(
        ("workload", "expected"),
        [("process", "process"), ("thread", "thread")],
    )
    def test_workload_names_are_reported_verbatim(self, workload, expected):
        assert Stage("x", workload, 1, 1).kind == expected


class TestPlanMatchesWhatTheStagesWillActuallyDo:
    def test_effective_matches_resolve_workers(self, tmp_path):
        """The plan must not re-implement the clamp; it must call it.

        A second copy of the clamping rule is exactly the kind of drift this
        module exists to prevent.
        """
        config = _config(tmp_path, threads=3, variant_threads=5)

        plan = {stage.name: stage for stage in describe(config)}

        assert plan["extraction"].effective == resolve_workers(3, "thread")
        assert plan["variant calling"].effective == resolve_workers(5, "process")
