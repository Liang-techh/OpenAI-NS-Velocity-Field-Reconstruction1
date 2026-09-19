from __future__ import annotations

import inspect

import pytest

from openai_ns_reconstruction.kokuno_finite_correction_damping import (
    DAMPING_GRID,
    _analytic_before,
    _analytic_direction,
    _regression_points,
    admit_damped_correction_search,
    audit_damped_correction_search,
    deterministic_receipt,
    truth_boundary,
)


def _run():
    held_in, held_out, update = _regression_points()
    return audit_damped_correction_search(
        _analytic_before(),
        _analytic_direction(),
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )


def test_oversized_direction_is_rescued_by_heldin_only_damping() -> None:
    report = _run()
    assert report.damping_grid == (1.0, 0.5, 0.25, 0.125)
    assert report.selected_alpha == pytest.approx(0.5, abs=0.0)
    assert report.selected_trial_index == 1
    assert report.heldout_used_for_alpha_selection is False
    assert report.heldout_trial_evaluation_count == 1
    assert report.trials[0].held_in_rms_ratio == pytest.approx(1.4, abs=2e-9)
    assert report.trials[0].held_in_eligible is False
    assert report.trials[1].held_in_rms_ratio == pytest.approx(0.2, abs=2e-9)
    assert report.trials[1].held_in_eligible is True
    assert report.trials[2].held_in_rms_ratio == pytest.approx(0.4, abs=2e-9)
    assert report.trials[3].held_in_rms_ratio == pytest.approx(0.7, abs=2e-9)
    assert report.selected_transition.held_out_rms_ratio == pytest.approx(0.2, abs=2e-9)
    assert report.selected_transition.finite_step_gain_guard_passed is True
    assert report.selected_step_admitted is True


def test_public_admission_accepts_selected_step_only_after_heldout() -> None:
    held_in, held_out, update = _regression_points()
    report = admit_damped_correction_search(
        _analytic_before(),
        _analytic_direction(),
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    assert report.selected_alpha == 0.5
    assert report.selected_step_admitted is True


def test_grid_is_frozen_and_cannot_be_threshold_tuned() -> None:
    held_in, held_out, update = _regression_points()
    with pytest.raises(ValueError, match="frozen damping grid"):
        audit_damped_correction_search(
            _analytic_before(),
            _analytic_direction(),
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=1.0e-3,
            damping_grid=(1.0, 0.4, 0.2),
        )
    assert DAMPING_GRID == (1.0, 0.5, 0.25, 0.125)


def test_heldin_heldout_overlap_is_rejected() -> None:
    held_in, held_out, update = _regression_points()
    held_out[0] = held_in[0]
    with pytest.raises(ValueError, match="disjoint"):
        audit_damped_correction_search(
            _analytic_before(),
            _analytic_direction(),
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_api_has_no_residual_gain_score_or_threshold_inputs() -> None:
    parameters = set(inspect.signature(audit_damped_correction_search).parameters)
    forbidden = {
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


def test_truth_boundary_stays_scientifically_fail_closed() -> None:
    boundary = truth_boundary()
    assert boundary["heldin_only_frozen_grid_damping_search_executable"] is True
    assert boundary["alpha_selection_uses_heldout"] is False
    assert boundary["heldout_evaluated_only_after_alpha_selection"] is True
    assert boundary["heldout_failure_causes_rejection_without_retuning"] is True
    assert boundary["formal_kokuno_correction_cycle_gain_bound_claimed"] is False
    assert boundary["real_full_candidate_defect_consumed"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["residual_reduction_claimed_for_real_candidate"] is False
    assert boundary["pde_validated"] is False
    assert boundary["final_normalized_momentum_gate"] == 1e-3
    assert boundary["final_normalized_divergence_gate"] == 1e-5


def test_deterministic_receipt_records_no_heldout_tuning() -> None:
    receipt = deterministic_receipt()
    report = receipt["analytic_mechanics_regression"]["report"]
    assert receipt["schema"] == "kokuno-a3-heldin-damped-correction-search-v1"
    assert report["selected_alpha"] == 0.5
    assert report["heldout_used_for_alpha_selection"] is False
    assert report["heldout_trial_evaluation_count"] == 1
    assert report["selected_step_admitted"] is True
