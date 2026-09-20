import inspect

import pytest

from openai_ns_reconstruction.kokuno_correction_cycle_gain_diagnostic import (
    admit_source_structured_finite_cycle,
    source_correction_stage_gain,
    truth_boundary,
)
from openai_ns_reconstruction.kokuno_finite_correction_cycle_ledger import (
    _analytic_cycle,
    _regression_points,
)


def test_stage_zero_matches_corrected_reader_arithmetic() -> None:
    stage = source_correction_stage_gain(0)
    assert stage.kappa_s == pytest.approx(0.00001, abs=0.0)
    assert stage.sigma_j == pytest.approx(0.2)
    assert stage.B_j == pytest.approx(0.7)
    assert stage.C_j == pytest.approx(1.2)
    assert stage.H_j == pytest.approx(1.19998)
    assert stage.covariance_error_gains == pytest.approx(
        (0.17998, 0.49997, 0.19997)
    )
    assert stage.covariance_min_gain == pytest.approx(0.17998)
    assert stage.wave_residual_gain == pytest.approx(0.39999)
    assert stage.tangential_mean_gain == pytest.approx(0.17)
    assert stage.compatibility_defect_gain == pytest.approx(0.89996)
    assert stage.H_minus_B_gain == pytest.approx(0.49998)
    assert stage.signed_increment_budget == pytest.approx(0.69999)
    assert stage.temporal_increment_budget == pytest.approx(1.19998)
    assert stage.next_B == pytest.approx(0.8)
    assert stage.next_C == pytest.approx(1.3)
    assert stage.required_stage_gain == pytest.approx(0.1)
    assert stage.source_bound_passed is True


def test_combined_gate_recomputes_actual_defect_and_preserves_source_gain() -> None:
    held_in, held_out, update = _regression_points()
    states, corrections = _analytic_cycle((1.0, 0.5, 0.25))
    report = admit_source_structured_finite_cycle(
        states,
        corrections,
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )
    assert report.step_count == 2
    assert report.source_gain_bounds_passed is True
    assert report.empirical_actual_defect_cycle_passed is True
    assert report.admitted is True
    assert report.minimum_source_wave_gain == pytest.approx(0.39999)
    assert report.minimum_source_mean_gain == pytest.approx(0.17)
    assert report.minimum_source_defect_gain == pytest.approx(0.89996)
    assert report.ledger.held_in_rms_contraction_factors == pytest.approx((0.5, 0.5))
    assert report.ledger.held_out_rms_contraction_factors == pytest.approx((0.5, 0.5))
    assert report.ledger.cumulative_held_in_rms_ratio == pytest.approx(0.25)
    assert report.ledger.cumulative_held_out_rms_ratio == pytest.approx(0.25)


def test_empirically_expanding_correction_is_rejected_before_source_label_can_help() -> None:
    held_in, held_out, update = _regression_points()
    states, corrections = _analytic_cycle((1.0, 0.5, 0.75))
    with pytest.raises(ValueError, match="actual-defect gain guard"):
        admit_source_structured_finite_cycle(
            states,
            corrections,
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )


def test_source_stage_relabel_and_public_gain_knobs_fail_closed() -> None:
    with pytest.raises(ValueError, match="non-negative"):
        source_correction_stage_gain(-1)

    held_in, held_out, update = _regression_points()
    later_states, later_corrections = _analytic_cycle((1.0, 0.5, 0.25))
    # Dropping stage 0 creates an otherwise consecutive transition starting at
    # j=1.  The combined gate must not let callers relabel a cycle as a later,
    # formally easier source stage.
    with pytest.raises(ValueError, match="start at stage index 0"):
        admit_source_structured_finite_cycle(
            later_states[1:],
            later_corrections[1:],
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )

    signature = inspect.signature(admit_source_structured_finite_cycle)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "debt",
        "delta_c",
        "inverse_rhs",
        "delta_y",
        "delta_a",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "required_gain",
        "kappa_s",
        "sigma_j",
        "B",
        "C",
        "damping_grid",
    }
    assert forbidden.isdisjoint(signature.parameters)

    boundary = truth_boundary()
    assert boundary["source_exponent_gain_arithmetic_executable"] is True
    assert boundary["source_gain_composed_with_actual_defect_ledger"] is True
    assert boundary["divergent_empirical_step_rejected"] is True
    assert boundary["caller_supplied_residual_allowed"] is False
    assert boundary["caller_supplied_source_gain_allowed"] is False
    assert boundary["formal_kokuno_correction_cycle_theorem_claimed"] is False
    assert boundary["real_candidate_finite_correction_cycle_run"] is False
    assert boundary["heldout_normalized_ns_residual_assessed"] is False
    assert boundary["pde_validated"] is False
