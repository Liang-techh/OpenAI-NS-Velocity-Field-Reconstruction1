from __future__ import annotations

import numpy as np

import agent7_st052m_compact_toroidal_swirl_preflight as audit


def test_symbolic_toroidal_divergence_and_truth_boundary() -> None:
    symbolic = audit.symbolic_toroidal_divergence()
    assert symbolic["exact_identity"] is True
    assert symbolic["divergence"] == "0"

    report = audit.build_report()
    assert report["task_id"] == "CR003-ST052M-COMPACT-TOROIDAL-SWIRL-PREFLIGHT-119"
    assert report["prereg_issue"] == 836
    assert report["source_parent"] == {
        "pr": 828,
        "head": "034a5e04df8900d8def05099fdbd1a287b43b4a0",
    }
    assert report["candidate_velocity_changed"] is False
    assert report["canonical_velocity_changed"] is False
    assert report["basis_dimension_changed"] is False
    assert report["toroidal_coefficient_selected"] is False
    assert report["new_temporal_basis_added"] is False
    assert report["parameter_scan_performed"] is False
    assert report["optimization_performed"] is False
    assert report["fresh_714_path_data_used"] is False
    assert report["pressure_or_force_changed"] is False
    assert report["held_out_pde_residual_evaluated"] is False
    assert report["public_image_numeric_target_used"] is False
    assert report["source_numeric_swirl_target_used"] is False
    assert report["closer_visualization_delivery_established"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False


def test_toroidal_channel_is_pure_swirl_and_compact() -> None:
    selectivity = audit.channel_selectivity()
    assert selectivity["toroidal_u_r_max_abs"] <= audit.ZERO_GATE
    assert selectivity["toroidal_u_z_max_abs"] <= audit.ZERO_GATE
    assert selectivity["toroidal_u_theta_rms"] > audit.NONTRIVIAL_RMS_MIN
    assert selectivity["poloidal_u_theta_max_abs"] <= audit.ZERO_GATE
    assert selectivity["poloidal_rz_rms"] > audit.NONTRIVIAL_RMS_MIN
    assert selectivity["passes"] is True

    support = audit.support_and_axis_check()
    assert support["toroidal_outside_max_abs"] <= audit.ZERO_GATE
    assert support["poloidal_outside_max_abs"] <= audit.ZERO_GATE
    assert support["toroidal_axis_transverse_max_abs"] <= audit.ZERO_GATE
    assert support["passes"] is True


def test_toroidal_and_poloidal_response_columns_are_independent() -> None:
    ident = audit.response_identifiability()
    assert ident["normalized_column_rank"] == 2
    assert abs(ident["normalized_column_cosine"]) <= audit.COSINE_MAX_ABS
    assert ident["normalized_column_condition"] <= audit.CONDITION_MAX
    assert ident["poloidal_response_norm"] > 0.0
    assert ident["toroidal_response_norm"] > 0.0
    assert ident["passes"] is True


def test_axis_and_cylindrical_sign_convention_with_simple_points() -> None:
    points = np.asarray(
        [
            [0.0, 0.0, 1.55],
            [0.55, 0.0, 1.55],
            [0.0, 0.55, 1.55],
            [-0.55, 0.0, -1.55],
        ],
        dtype=float,
    )
    toroidal = audit.toroidal_unit_velocity(points)
    cyl = audit.cylindrical_components(points, toroidal)
    assert np.max(np.abs(cyl[:, 0])) <= audit.ZERO_GATE
    assert np.max(np.abs(cyl[:, 2])) <= audit.ZERO_GATE
    assert np.max(np.abs(toroidal[0])) <= audit.ZERO_GATE
    assert cyl[1, 1] > 0.0
    assert cyl[2, 1] > 0.0
    assert cyl[3, 1] > 0.0
