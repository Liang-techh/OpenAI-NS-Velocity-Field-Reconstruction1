import pytest

from openai_ns_reconstruction.constrained_eq45_axial_envelope_small_alpha_screen import (
    TRUTH_BOUNDARY,
    audit_axial_envelope_small_alpha_screen,
)


def _small_report():
    return audit_axial_envelope_small_alpha_screen(
        alphas=(0.0, 0.10, 0.25),
        quadrature_orders=(24, 32),
        morphology_times=(0.25, 0.50, 0.75),
        morphology_grid_size=21,
    )


def test_small_alpha_screen_is_bounded_and_truthful():
    report = _small_report()
    assert report["screen_protocol"]["alphas"] == [0.0, 0.10, 0.25]
    assert report["screen_protocol"]["new_basis_shapes_added"] == 0
    assert report["screen_protocol"]["parameter_count_added_for_screen"] == 0
    assert len(report["energy_rows"]) == 3
    assert len(report["morphology_rows"]) == 9
    assert len(report["diagnostic_time_rows"]) == 3
    assert report["energy_quadrature_max_relative_refinement_change"] >= 0.0

    baseline = report["diagnostic_time_rows"][0]
    assert baseline["alpha"] == 0.0
    assert abs(baseline["q90_gain_over_Zp"]) < 1.0e-14
    assert abs(baseline["q99_gain_over_Zp"]) < 1.0e-14
    assert abs(baseline["axial_rms_relative_change"]) < 1.0e-14
    assert abs(baseline["radial_rms_relative_change"]) < 1.0e-14
    assert abs(baseline["collar_enstrophy_fraction_absolute_change"]) < 1.0e-14

    checks = report["max_alpha_structure_checks"]
    assert checks["alpha_zero_parent_replay_max_abs_velocity_error"] < 1.0e-12
    assert checks["outside_physical_support_max_abs_velocity"] < 1.0e-12
    assert checks["cartesian_fd_divergence_max_abs"] < 1.0e-5
    assert checks["representative_core_signs_all_pass"]

    assert report["held_out_pde_residual"]["evaluated"] is False
    assert report["smallest_mover_is_production_selection"] is False
    assert report["truth_boundary"] == TRUTH_BOUNDARY
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False


def test_screen_rejects_alpha_outside_bounded_line():
    with pytest.raises(ValueError, match="alpha<=0.25"):
        audit_axial_envelope_small_alpha_screen(
            alphas=(0.0, 0.30),
            quadrature_orders=(16, 24),
            morphology_times=(0.5,),
            morphology_grid_size=17,
        )
