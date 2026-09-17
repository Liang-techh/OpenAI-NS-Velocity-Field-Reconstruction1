import pytest

from openai_ns_reconstruction.constrained_eq45_supported_collar_vorticity import (
    governed_supported_seed,
)
from openai_ns_reconstruction.constrained_eq45_supported_early_radial_capacity import (
    audit_supported_early_radial_capacity,
)


def test_supported_early_radial_capacity_static_span_is_selective_but_ill_conditioned():
    report = audit_supported_early_radial_capacity(governed_supported_seed())

    baseline = report["baseline_feature"]
    assert 0.05 < baseline[0] < 0.07
    assert 0.003 < baseline[1] < 0.005
    assert baseline[2] < 5e-6

    assert report["numerical_rank"] == 2
    assert 250.0 < report["condition_number"] < 350.0
    assert report["response_refinement_relative_change"] < 2e-4

    projection = report["early_only_projection"]
    assert projection["projection_fraction"] > 0.999
    assert projection["target_residual_fraction"] < 1e-3
    assert projection["projected_late_leak_l2"] < 1e-3
    # The near-perfect static signature is obtained by a large cancellation,
    # not by a well-conditioned low-amplitude direction.
    direction = projection["least_squares_direction"]
    assert abs(direction[0]) > 100.0
    assert abs(direction[1]) > 500.0

    trial = report["bounded_static_local_trial"]
    assert max(abs(value) for value in trial["coefficient_deltas"].values()) == pytest.approx(0.02)
    assert 5e-4 < trial["relative_early_reduction"] < 7e-4
    assert abs(trial["actual_delta"][1]) < 1e-8
    assert abs(trial["actual_delta"][2]) < 3e-8
    assert trial["linearization_relative_error"] < 5e-3


def test_supported_early_radial_capacity_favors_existing_phi10_affine_followup():
    report = audit_supported_early_radial_capacity(governed_supported_seed())
    screen = report["affine_temporal_screen"]

    assert screen["tau"] == pytest.approx([-1.0, 0.0, 1.0])
    assert screen["recommended_existing_mode_for_temporal_followup"] == "(1, 0)"
    assert screen["production_temporal_candidate_implemented_here"] is False

    phi10 = screen["modes"]["(1, 0)"]
    phi12 = screen["modes"]["(1, 2)"]
    assert phi10["early_to_late_selectivity"] > 1e4
    assert phi10["same_local_budget_efficiency_vs_static_cancellation"] > 20.0
    assert 1.3 < phi10["linearized_abs_slope_to_zero_early_feature"] < 1.6
    assert phi10["linearized_zero_early_within_bound"] is True
    assert phi10["max_abs_slope_under_existing_bound"] == pytest.approx(3.7)

    assert phi12["early_to_late_selectivity"] > 3e3
    assert phi12["linearized_abs_slope_to_zero_early_feature"] > 9.0
    assert phi12["linearized_zero_early_within_bound"] is False
    assert phi12["max_abs_slope_under_existing_bound"] == pytest.approx(3.9)


def test_supported_early_radial_capacity_rejects_bad_inputs():
    child = governed_supported_seed()
    with pytest.raises(ValueError):
        audit_supported_early_radial_capacity(child, coefficient_steps=(0.02, 0.04))
    with pytest.raises(ValueError):
        audit_supported_early_radial_capacity(child, times=(0.25, 0.50, 0.80))
    with pytest.raises(ValueError):
        audit_supported_early_radial_capacity(child, resolution=19)


def test_supported_early_radial_capacity_does_not_mutate_or_promote_candidate():
    child = governed_supported_seed()
    before = child.sha256
    report = audit_supported_early_radial_capacity(child, resolution=25)

    assert child.sha256 == before
    assert report["parameter_count_added"] == 0
    assert report["existing_parameters_audited"] == 2
    assert report["interpretation"]["new_spatial_basis_added"] is False
    assert report["interpretation"]["public_visual_target_fitted"] is False

    truth = report["truth_boundary"]
    assert truth["velocity_changed"] is False
    assert truth["production_coefficients_changed"] is False
    assert truth["new_basis_added"] is False
    assert truth["temporal_mode_implemented"] is False
    assert truth["forcing_or_pressure_refit"] is False
    assert truth["visualization_ready"] is False
    assert truth["visual_correspondence_verified"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert truth["blowup_proved"] is False
