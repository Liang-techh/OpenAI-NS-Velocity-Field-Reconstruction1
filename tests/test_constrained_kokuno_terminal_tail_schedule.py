from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_outer_reserved_patch_schedule import (
    KokunoOuterReservedPatchSchedule,
)
from openai_ns_reconstruction.kokuno_terminal_tail_schedule import (
    KokunoTerminalTailSchedule,
)


def test_terminal_release_q_chain_replays_source_endpoint_contract() -> None:
    schedule = KokunoTerminalTailSchedule()

    assert schedule.Q_release == pytest.approx(
        (schedule.lambda_outer - schedule.h) / (1.0 - schedule.lambda_outer)
    )
    assert schedule.Q_after_steep_hold == pytest.approx(
        schedule.Q_after_first_transition
        + (1.0 - schedule.h) * schedule.steep_hold_length,
        rel=0.0,
        abs=2.0e-12,
    )
    assert schedule.Q_in > schedule.Q_p > 0.0
    assert schedule.release_hold_length == pytest.approx(
        math.log(schedule.Q_in / schedule.Q_p) / (1.0 - schedule.h),
        rel=2.0e-13,
    )

    terminal = schedule.terminal_Q(np.asarray([0.0, 1.0, 2.0, 3.0]))
    assert terminal[0] == pytest.approx(schedule.Q_p, rel=2.0e-12, abs=2.0e-14)
    assert terminal[-1] == 0.0
    assert np.all(np.diff(terminal) < 0.0)

    slope = schedule.terminal_log_slope(np.linspace(0.0, 3.0, 4097))
    assert np.min(slope) >= 0.0
    assert np.max(slope) < schedule.h / 4.0


def test_terminal_scale_chain_retains_exact_source_factors() -> None:
    schedule = KokunoTerminalTailSchedule()
    outer = schedule.outer_schedule

    assert schedule.log_X_p == pytest.approx(outer.log_X_w + outer.T_w)
    assert schedule.log_e_b == pytest.approx(
        outer.log_e_w - (0.5 + outer.lambda_outer) * outer.T_w
    )
    assert schedule.log_X_rel == pytest.approx(
        schedule.log_X_p
        + 13.0 / outer.lambda_outer
        + schedule.T_f
        + 30.0 * math.log(1.0 / outer.lambda_outer)
    )

    log_after_first_release = schedule.log_e_rel - 1.0 - 0.5 * outer.lambda_outer
    log_after_steep = log_after_first_release - 1.5 * schedule.steep_hold_length
    assert log_after_steep - log_after_first_release == pytest.approx(
        6.0 * math.log(outer.h), rel=0.0, abs=2.0e-13
    )

    f0 = float(schedule.terminal_multiplier(np.asarray(0.0)))
    assert (
        schedule.log_c_inf
        - schedule.A * schedule.log_X_tail
        + math.log(f0)
    ) == pytest.approx(schedule.log_e_tail_start, rel=0.0, abs=2.0e-13)
    assert schedule.log_X_K - schedule.log_X_tail == pytest.approx(0.2)
    assert schedule.log_X_b - schedule.log_X_tail == pytest.approx(3.0)


def test_default_existence_style_tail_is_kept_in_log_space_fail_closed() -> None:
    schedule = KokunoTerminalTailSchedule()
    report = schedule.log_scale_report()

    assert schedule.log_X_tail > math.log(np.finfo(float).max)
    assert report["X_tail"] is None
    assert report["c_inf"] is not None
    assert report["repair_ratios"]["nominal_X_K_inverse"] is None
    assert report["repair_ratios"]["log_X_star_over_X_K"] < 0.0
    assert report["repair_ratios"]["log_e_star_over_e_K"] > 0.0

    with pytest.raises(OverflowError, match="log_scale_report"):
        schedule.materialize_heat_inputs()


def test_moderate_autonomous_schedule_materializes_terminal_heat_inputs() -> None:
    outer = KokunoOuterReservedPatchSchedule(
        M_d=1.0,
        lambda_outer=0.2,
        h=0.005,
    )
    schedule = KokunoTerminalTailSchedule(outer_schedule=outer)
    inputs = schedule.materialize_heat_inputs()

    assert inputs["h"] == outer.h
    assert inputs["rho_o"] == pytest.approx(schedule.c_o * outer.h)
    assert math.isfinite(inputs["X_tail"]) and inputs["X_tail"] > 0.0
    assert math.isfinite(inputs["c_inf"]) and inputs["c_inf"] > 0.0
    assert math.log(inputs["X_tail"]) == pytest.approx(schedule.log_X_tail)
    assert math.log(inputs["c_inf"]) == pytest.approx(schedule.log_c_inf)


def test_autonomous_choices_are_guarded_by_source_slope_inequalities() -> None:
    with pytest.raises(ValueError, match="interpolation slope band"):
        KokunoTerminalTailSchedule(T_f=50.0)

    with pytest.raises(ValueError, match="terminal bound"):
        KokunoTerminalTailSchedule(c_o=0.07)


def test_terminal_schedule_payload_round_trip_and_truth_mutation_fail_closed(tmp_path) -> None:
    schedule = KokunoTerminalTailSchedule()
    target = tmp_path / "tail_schedule.json"
    schedule.save_json(target)
    loaded = KokunoTerminalTailSchedule.load_json(target)

    assert loaded.sha256 == schedule.sha256
    assert loaded.to_payload() == schedule.to_payload()

    payload = copy.deepcopy(schedule.to_payload())
    assert payload["truth_boundary"]["pde_validated"] is False
    assert payload["truth_boundary"]["paper_exact"] is False
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="hash or content mismatch"):
        KokunoTerminalTailSchedule.from_payload(payload)


def test_source_vs_autonomous_provenance_is_explicit() -> None:
    schedule = KokunoTerminalTailSchedule()
    payload = schedule.to_payload()
    truth = payload["truth_boundary"]

    assert truth["source_post_pulse_scale_chain_executable"] is True
    assert truth["source_release_Q_equation_executable"] is True
    assert truth["source_X_tail_relation_executable"] is True
    assert truth["source_c_inf_continuity_relation_executable"] is True
    assert truth["autonomous_T_f"] is True
    assert truth["autonomous_c_o"] is True
    assert truth["source_hidden_numeric_choices_recovered"] is False
    assert truth["actual_heat_discrepancy_applied_to_patch"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
