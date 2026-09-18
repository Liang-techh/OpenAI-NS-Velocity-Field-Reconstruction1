import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_bipolar_axial_cap_poloidal_capacity import (
    AXIAL_CAP_POLOIDAL,
    CAP_Q_INNER,
    CAP_Q_OUTER,
    _axial_cap_basis_velocity,
    _basis_axial_localization,
    _basis_structure_checks,
    _cap_velocity,
    audit_axial_cap_poloidal_capacity,
)
from openai_ns_reconstruction.constrained_eq45_bipolar_compact_poloidal_capacity import _source_field


def test_axial_cap_basis_is_banded_divergence_free_and_parity_compatible():
    field = _source_field()
    checks = _basis_structure_checks(field)
    assert checks["radial_even_max_abs_error"] < 1e-12
    assert checks["axial_odd_max_abs_error"] < 1e-12
    assert checks["swirl_response_rms"] < 1e-12
    assert checks["support_flank_max_abs_velocity"] < 1e-14
    assert checks["midplane_max_abs_velocity"] < 1e-14
    assert checks["inner_axial_band_max_abs_velocity"] < 1e-14
    assert checks["cartesian_fd_divergence_max_abs"] < 1e-7

    points = np.array([[0.45, 0.0, 1.12], [0.45, 0.0, -1.12]], dtype=float)
    response = _axial_cap_basis_velocity(field, points)
    assert np.linalg.norm(response) > 0.0
    assert response[0, 0] == pytest.approx(response[1, 0], abs=1e-14)
    assert response[0, 2] == pytest.approx(-response[1, 2], abs=1e-14)
    assert CAP_Q_INNER == pytest.approx(0.55)
    assert CAP_Q_OUTER == pytest.approx(0.95)


def test_axial_cap_support_moves_response_outward_of_polynomial_shoulder():
    field = _source_field()
    localization = _basis_axial_localization(field, grid_size=21)
    assert localization["cap_abs_z_energy_centroid_over_Zp"] > localization[
        "shoulder_abs_z_energy_centroid_over_Zp"
    ]
    assert localization["cap_to_shoulder_centroid_ratio"] > 1.05
    assert localization["cap_outer_065_response_energy_fraction"] > localization[
        "shoulder_outer_065_response_energy_fraction"
    ]
    assert localization["cap_outer_075_response_energy_fraction"] > 0.0


def test_axial_cap_screen_adds_one_stable_target_free_direction():
    report = audit_axial_cap_poloidal_capacity(
        coefficient_steps=(0.02, 0.01),
        grid_size=21,
        diagnostic_trial_coefficient=0.20,
    )
    velocity = report["public_velocity_response"]
    assert velocity["normalized_rank"] == 5
    assert np.isfinite(velocity["normalized_condition_number"])
    assert velocity["normalized_condition_number"] < 1.0e4
    assert velocity["finite_difference_refinement_relative_change"] < 1e-6
    assert velocity["cap_novelty_outside_Phi01_Phi03_center_compact_shoulder_span"] > 0.05

    morphology = report["vorticity_morphology"]
    assert morphology["response_rank"]["normalized_rank"] >= 3
    assert np.isfinite(morphology["response_rank"]["normalized_condition_number"])
    assert morphology["finite_difference_refinement_relative_change"] < 0.10
    assert morphology["diagnostic_axial_rms_span_over_Zp"] > 0.0
    assert morphology["diagnostic_axial_q90_span_over_Zp"] >= 0.0
    assert morphology["diagnostic_axial_q99_span_over_Zp"] >= 0.0
    assert morphology["diagnostic_outer_065_span"] >= 0.0
    assert morphology["diagnostic_outer_075_span"] >= 0.0
    assert report["parameter_growth_if_materialized"]["fixed_spatial_basis_shapes_added"] == 1
    assert report["parameter_growth_if_materialized"]["scalar_coefficients_added"] == 1
    assert report["truth_boundary"]["velocity_changed"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False


def test_axial_cap_inputs_fail_closed():
    field = _source_field()
    with pytest.raises(ValueError):
        _axial_cap_basis_velocity(field, np.zeros((2, 2)))
    with pytest.raises(ValueError, match="implementation guard"):
        _cap_velocity(field, np.zeros((1, 3)), 0.5, 5.0)
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_capacity(coefficient_steps=(0.01, 0.02))
    with pytest.raises(ValueError):
        audit_axial_cap_poloidal_capacity(grid_size=16)
    assert AXIAL_CAP_POLOIDAL == "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
