import pytest

from openai_ns_reconstruction.constrained_eq45_axial_envelope_joint_renorm_screen import (
    TRUTH_BOUNDARY,
    audit_axial_envelope_joint_renorm_screen,
)


def _small_report():
    return audit_axial_envelope_joint_renorm_screen(
        alphas=(0.0, 0.25, 0.50, 1.0),
        quadrature_orders=(24, 32),
        morphology_times=(0.25, 0.50, 0.75),
        morphology_grid_size=21,
    )


def test_joint_renorm_screen_restores_energy_and_preserves_truth_boundary():
    report = _small_report()
    assert report["screen_protocol"]["alphas"] == [0.0, 0.25, 0.50, 1.0]
    assert report["screen_protocol"]["new_basis_shapes_added"] == 0
    assert report["screen_protocol"]["parameter_count_added_for_screen"] == 0
    assert report["screen_protocol"]["production_alpha_selected"] is False
    assert len(report["rows"]) == 4
    assert report["energy_quadrature_max_relative_refinement_change"] >= 0.0

    baseline = report["rows"][0]
    assert baseline["baseline_parent_replay_without_common_scale"] is True
    assert baseline["common_velocity_scale"] == 1.0
    assert baseline["reference_energy_gate_pass"]
    assert abs(baseline["jointly_renormalized_reference_energy"] - 1.0) <= 1.0e-3
    assert max(abs(value) for value in baseline["q90_gain_over_Zp_by_time"]) < 1.0e-14
    assert max(abs(value) for value in baseline["q99_gain_over_Zp_by_time"]) < 1.0e-14

    for row in report["rows"]:
        assert row["reference_energy_gate_pass"]
        assert row["validation_energy_range_all_pass"]
        assert row["profile_bound_screen"]["simple_common_scale_representation_preflight_passed"]
        assert row["representation_energy_preflight_passed"]
        assert row["positive_common_scale_morphology_invariance"]
        assert row["navier_stokes_balance_invariant_under_common_scale"] is False
        assert len(row["q90_gain_over_Zp_by_time"]) == 3
        assert len(row["q99_gain_over_Zp_by_time"]) == 3
    for row in report["rows"][1:]:
        assert row["baseline_parent_replay_without_common_scale"] is False
        assert abs(row["jointly_renormalized_reference_energy"] - 1.0) < 1.0e-10
        assert row["common_velocity_scale"] < 1.0

    checks = report["max_alpha_structure_checks_before_common_scale"]
    assert checks["alpha_zero_parent_replay_max_abs_velocity_error"] < 1.0e-12
    assert checks["outside_physical_support_max_abs_velocity"] < 1.0e-12
    assert checks["cartesian_fd_divergence_max_abs"] < 1.0e-5
    assert checks["representative_core_signs_all_pass"]

    assert report["held_out_pde_residual"]["evaluated"] is False
    assert report["smallest_mover_is_production_selection"] is False
    assert report["truth_boundary"] == TRUTH_BOUNDARY
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False


def test_joint_renorm_screen_rejects_nonmonotone_or_out_of_range_alpha():
    with pytest.raises(ValueError, match="strictly increasing"):
        audit_axial_envelope_joint_renorm_screen(
            alphas=(0.0, 0.5, 0.25),
            quadrature_orders=(16, 24),
            morphology_times=(0.5,),
            morphology_grid_size=17,
        )
    with pytest.raises(ValueError, match="alpha"):
        audit_axial_envelope_joint_renorm_screen(
            alphas=(0.0, 1.01),
            quadrature_orders=(16, 24),
            morphology_times=(0.5,),
            morphology_grid_size=17,
        )
