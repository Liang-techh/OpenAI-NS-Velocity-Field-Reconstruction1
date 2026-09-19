from __future__ import annotations

import copy
import inspect
from dataclasses import replace

import pytest

import openai_ns_reconstruction.kokuno_finite_correction_gain_ladder as gain_ladder
from openai_ns_reconstruction.kokuno_finite_correction_damping import (
    _analytic_before,
    _analytic_direction,
    _regression_points,
)
from openai_ns_reconstruction.kokuno_finite_correction_gain_ladder import (
    PROJECT_DERIVATIVE_STEPS,
    admit_project_derivative_ladder_gain,
    audit_project_derivative_ladder_gain,
    deterministic_receipt,
    truth_boundary,
)


def _run():
    held_in, held_out, update = _regression_points()
    return audit_project_derivative_ladder_gain(
        _analytic_before(), _analytic_direction(), held_in, held_out, update
    )


def test_selected_step_improves_actual_defect_across_project_derivative_ladder() -> None:
    report = _run()
    assert report.selected_alpha == pytest.approx(0.5, abs=0.0)
    assert report.alpha_selected_at_spatial_step == pytest.approx(0.005, abs=0.0)
    assert report.project_derivative_steps == (0.02, 0.01, 0.005)
    assert report.heldout_used_for_alpha_selection is False
    assert report.alpha_retuned_across_derivative_steps is False
    assert report.all_derivative_levels_gain_guard_passed is True
    for level in report.levels:
        assert level.held_in_rms_ratio == pytest.approx(0.2, abs=2e-9)
        assert level.held_out_rms_ratio == pytest.approx(0.2, abs=2e-9)
        assert level.held_in_max_ratio == pytest.approx(0.2, abs=2e-9)
        assert level.held_out_max_ratio == pytest.approx(0.2, abs=2e-9)
        assert level.correction_divergence_max <= 2e-12
        assert level.finite_step_gain_guard_passed is True


def test_public_admission_rejects_if_any_registered_step_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    original = gain_ladder.audit_finite_correction_transition

    def one_failed_level(*args, **kwargs):
        report = original(*args, **kwargs)
        if kwargs.get("spatial_step") == PROJECT_DERIVATIVE_STEPS[0]:
            return replace(report, finite_step_gain_guard_passed=False)
        return report

    monkeypatch.setattr(gain_ladder, "audit_finite_correction_transition", one_failed_level)
    held_in, held_out, update = _regression_points()
    with pytest.raises(ValueError, match="derivative-ladder gain audit"):
        admit_project_derivative_ladder_gain(
            _analytic_before(), _analytic_direction(), held_in, held_out, update
        )


def test_project_protocol_drift_is_fail_closed() -> None:
    payload = gain_ladder._load_project_constraints()
    bad_steps = copy.deepcopy(payload)
    bad_steps["validation"]["derivative_steps"] = [0.02, 0.01, 0.004]
    with pytest.raises(ValueError, match="derivative ladder drifted"):
        gain_ladder._validate_project_constraints(bad_steps)

    bad_gate = copy.deepcopy(payload)
    bad_gate["validation"]["thresholds"]["pde_residual_L2"] = 0.01
    with pytest.raises(ValueError, match="scientific threshold drifted"):
        gain_ladder._validate_project_constraints(bad_gate)


def test_public_api_cannot_replace_ladder_threshold_or_residual() -> None:
    parameters = set(inspect.signature(audit_project_derivative_ladder_gain).parameters)
    forbidden = {
        "spatial_step",
        "derivative_steps",
        "damping_grid",
        "residual",
        "defect",
        "stress",
        "target",
        "gain",
        "score",
        "normalized_score",
        "threshold",
        "momentum_threshold",
        "divergence_threshold",
    }
    assert not (parameters & forbidden)


def test_truth_boundary_does_not_launder_sampled_rms_into_final_pde_l2() -> None:
    boundary = truth_boundary()
    assert boundary["project_derivative_ladder_gain_audit_executable"] is True
    assert boundary["heldout_used_for_alpha_selection"] is False
    assert boundary["alpha_retuned_across_derivative_steps"] is False
    assert boundary["project_thresholds_fail_closed_against_config_drift"] is True
    assert boundary["sampled_rms_relabelled_as_volume_weighted_project_l2"] is False
    assert boundary["formal_kokuno_correction_cycle_gain_bound_claimed"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed_for_real_candidate"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_max_gate"] == 1e-3
    assert boundary["final_normalized_momentum_l2_gate"] == 1e-3
    assert boundary["final_normalized_divergence_max_gate"] == 1e-5
    assert boundary["final_normalized_divergence_l2_gate"] == 1e-5


def test_deterministic_receipt_records_one_frozen_alpha_across_ladder() -> None:
    receipt = deterministic_receipt()
    assert receipt["schema"] == "kokuno-a3-finite-correction-gain-ladder-v1"
    report = receipt["analytic_mechanics_regression"]["report"]
    assert report["selected_alpha"] == 0.5
    assert report["project_derivative_steps"] == [0.02, 0.01, 0.005]
    assert report["heldout_used_for_alpha_selection"] is False
    assert report["alpha_retuned_across_derivative_steps"] is False
    assert report["all_derivative_levels_gain_guard_passed"] is True
    assert max(level["held_out_rms_ratio"] for level in report["levels"]) == pytest.approx(
        0.2, abs=2e-9
    )
