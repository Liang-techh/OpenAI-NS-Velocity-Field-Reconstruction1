from __future__ import annotations

import numpy as np

import agent7_st052m_total_child_axial_stretch_fingerprint as audit


def test_preregistered_probe_contract_is_frozen() -> None:
    assert audit.TASK_ID == "CR003-ST052M-TOTAL-CHILD-AXIAL-STRETCH-FINGERPRINT-113"
    assert audit.PREREG_ISSUE == 793
    assert audit.SOURCE_CHILD_PR == 775
    assert audit.SOURCE_CHILD_HEAD == "93887a59729d22113badf2ae4dac2f7d868e2703"
    assert audit.SOURCE_SCOPE_PR == 786
    assert audit.SOURCE_SCOPE_HEAD == "f3c1e5a2a86bcc09c9d1f16e9f865c5d57a94ae3"
    assert audit.RELATED_FRESH_HOLDOUT_PR == 714
    assert audit.ABS_Z_BANDS == (1.55, 1.75)
    assert audit.PROBE_TIMES == (0.375, 0.50, 0.625, 0.75)
    assert len(audit.RADII) == 55
    assert abs(float(audit.RADII[0]) - 0.20) < 1.0e-15
    assert abs(float(audit.RADII[-1]) - 1.55) < 1.0e-15
    assert len(audit.ANGLES) == 8


def test_ring_cloud_balances_both_z_signs_and_angles() -> None:
    points, signs = audit._ring_cloud(1.55)
    assert points.shape == (55 * 2 * 8, 3)
    assert signs.shape == (55 * 2 * 8,)
    assert np.sum(signs < 0.0) == 55 * 8
    assert np.sum(signs > 0.0) == 55 * 8
    assert np.allclose(np.abs(points[:, 2]), 1.55)
    radius = np.sqrt(points[:, 0] ** 2 + points[:, 1] ** 2)
    assert np.isclose(radius.min(), 0.20)
    assert np.isclose(radius.max(), 1.55)


def test_crossing_helper_recovers_manufactured_central_positive_profile() -> None:
    radii = audit.RADII
    profile = 1.10 - radii
    crossing = audit.first_positive_to_nonpositive_crossing(radii, profile)
    assert crossing is not None
    assert abs(crossing - 1.10) < 1.0e-12
    assert audit.CROSSING_RADIUS_MIN < crossing < audit.CROSSING_RADIUS_MAX
    assert audit._contiguous_central_positive_fraction(profile) > 0.0
    assert audit.first_positive_to_nonpositive_crossing(radii, -np.ones_like(radii)) is None


def test_truth_boundary_forbids_retuning_and_scientific_promotion() -> None:
    assert audit.TRUTH["frozen_775_total_child_reconstructed"] is True
    assert audit.TRUTH["historical_development_paths_used_only_to_rederive_frozen_775_alpha"] is True
    for key in (
        "velocity_formula_changed_this_increment",
        "candidate_coefficient_changed_this_increment",
        "basis_dimension_changed",
        "new_spatial_basis_added",
        "new_temporal_basis_added",
        "fresh_714_path_data_used",
        "fresh_714_path_data_used_for_retuning",
        "parameter_scan_performed",
        "optimization_performed",
        "post_result_retuning_performed",
        "pressure_or_force_changed",
        "held_out_pde_residual_evaluated",
        "public_image_numeric_target_used",
        "pixel_similarity_objective_used",
        "canonical_velocity_changed",
        "saved_velocity_changed",
        "production_candidate_selected",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert audit.TRUTH[key] is False
