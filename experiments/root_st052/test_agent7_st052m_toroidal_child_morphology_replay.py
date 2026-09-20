from __future__ import annotations

import numpy as np

import agent7_st052m_toroidal_child_morphology_replay as audit


def test_projection_rule_strengthens_midpoint_swirl_without_radial_axial_leakage() -> None:
    parent_fn, _, _ = audit._build_reservoir_parent()
    beta, calibration = audit.derive_beta(parent_fn)
    assert calibration["coefficient_scan_performed"] is False
    assert calibration["optimizer_used"] is False
    assert calibration["source_numeric_target_used"] is False
    assert calibration["fresh_714_data_used"] is False
    assert calibration["nontrivial"] is True

    fingerprint = audit.morphology_fingerprint(parent_fn, beta)
    assert fingerprint["finite"] is True
    assert fingerprint["midpoint_swirl_rms_strictly_increased"] is True
    assert fingerprint["radial_axial_invariant_all_times"] is True
    midpoint = fingerprint["by_time"]["0.500"]
    assert midpoint["u_theta_rms_delta"] > audit.SWIRL_INCREASE_ABS_MIN
    assert midpoint["radial_component_change_max_abs"] <= audit.ZERO_GATE
    assert midpoint["axial_component_change_max_abs"] <= audit.ZERO_GATE


def test_child_is_exact_parent_at_start_and_toroidal_increment_is_pure() -> None:
    parent_fn, _, _ = audit._build_reservoir_parent()
    beta, _ = audit.derive_beta(parent_fn)
    points = audit._probe_points()

    parent_start = np.asarray(parent_fn(points, 0.25), dtype=float)
    child_start = np.asarray(audit.child_velocity(parent_fn, points, 0.25, beta), dtype=float)
    assert np.max(np.abs(child_start - parent_start)) <= audit.ZERO_GATE

    time = 0.50
    parent_cyl = audit.toroidal.cylindrical_components(points, parent_fn(points, time))
    child_cyl = audit.toroidal.cylindrical_components(
        points, audit.child_velocity(parent_fn, points, time, beta)
    )
    np.testing.assert_allclose(child_cyl[:, 0], parent_cyl[:, 0], atol=audit.ZERO_GATE, rtol=0.0)
    np.testing.assert_allclose(child_cyl[:, 2], parent_cyl[:, 2], atol=audit.ZERO_GATE, rtol=0.0)
    assert np.linalg.norm(child_cyl[:, 1] - parent_cyl[:, 1]) > 0.0


def test_structure_capacity_and_truth_boundary() -> None:
    parent_fn, _, _ = audit._build_reservoir_parent()
    beta, _ = audit.derive_beta(parent_fn)
    structure = audit.structure_check(parent_fn, beta)
    assert structure["passes"] is True
    assert structure["symbolic_toroidal_divergence_exact"] is True
    assert structure["poloidal_toroidal_capacity"]["passes"] is True
    assert structure["start_identity_max_abs"] <= audit.ZERO_GATE
    assert structure["outside_unit_toroidal_max_abs"] <= audit.ZERO_GATE
    assert structure["toroidal_correction_cartesian_fd_divergence_max"] <= audit.DIVERGENCE_MAX

    assert audit.TRUTH["basis_dimension_increment"] == 1
    assert audit.TRUTH["second_poloidal_basis_added"] is False
    assert audit.TRUTH["held_out_pde_residual_evaluated"] is False
    assert audit.TRUTH["pressure_or_force_changed"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
    assert audit.TRUTH["paper_exact"] is False
    assert audit.TRUTH["openai_field_identified"] is False
