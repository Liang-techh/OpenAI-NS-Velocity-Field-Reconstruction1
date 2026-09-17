import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_bipolar_compact_poloidal_capacity import (
    COMPACT_POLOIDAL,
    PHI03,
    _basis_structure_checks,
    _compact_poloidal_basis_velocity,
    _source_field,
    _with_phi_delta,
    audit_compact_poloidal_axial_capacity,
)


def test_compact_poloidal_basis_is_parity_compatible_divergence_free_and_collar_zero():
    field = _source_field()
    checks = _basis_structure_checks(field)
    assert checks["radial_even_max_abs_error"] < 1e-12
    assert checks["axial_odd_max_abs_error"] < 1e-12
    assert checks["swirl_response_rms"] < 1e-12
    assert checks["support_flank_max_abs_velocity"] < 1e-14
    assert checks["cartesian_fd_divergence_max_abs"] < 1e-7

    points = np.array(
        [
            [0.45, 0.0, 0.0],
            [0.25, 0.0, 0.40],
            [0.25, 0.0, -0.40],
        ],
        dtype=float,
    )
    response = _compact_poloidal_basis_velocity(field, points)
    assert response[0, 0] < 0.0  # positive coefficient: inward midplane radial response
    assert response[1, 2] > 0.0
    assert response[2, 2] < 0.0
    assert response[1, 2] == pytest.approx(-response[2, 2], abs=1e-14)


def test_compact_poloidal_screen_adds_one_stable_target_free_direction():
    report = audit_compact_poloidal_axial_capacity(
        coefficient_steps=(0.02, 0.01),
        grid_size=17,
        diagnostic_trial_coefficient=0.20,
    )
    velocity = report["public_velocity_response"]
    assert velocity["normalized_rank"] == 3
    assert np.isfinite(velocity["normalized_condition_number"])
    assert velocity["normalized_condition_number"] < 100.0
    assert velocity["finite_difference_refinement_relative_change"] < 1e-6
    assert velocity["compact_novelty_outside_Phi01_Phi03_span"] > 0.05

    morphology = report["vorticity_morphology"]
    assert morphology["response_rank"]["normalized_rank"] >= 2
    assert morphology["finite_difference_refinement_relative_change"] < 0.10
    assert morphology["diagnostic_axial_rms_span_over_Zp"] > 0.0
    assert morphology["diagnostic_outer_axial_enstrophy_fraction_span"] > 0.0
    assert report["parameter_growth_if_materialized"]["fixed_spatial_basis_shapes_added"] == 1
    assert report["parameter_growth_if_materialized"]["scalar_coefficients_added"] == 1
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False


def test_phi03_bound_and_malformed_inputs_fail_closed():
    field = _source_field()
    with pytest.raises(ValueError, match="coefficient_limit"):
        _with_phi_delta(field, PHI03, -0.01)
    with pytest.raises(ValueError):
        _compact_poloidal_basis_velocity(field, np.zeros((2, 2)))
    with pytest.raises(ValueError):
        audit_compact_poloidal_axial_capacity(coefficient_steps=(0.01, 0.02))
    with pytest.raises(ValueError):
        audit_compact_poloidal_axial_capacity(grid_size=16)
    assert COMPACT_POLOIDAL == "COMPACT_C4_ODD_Z_POLOIDAL"
