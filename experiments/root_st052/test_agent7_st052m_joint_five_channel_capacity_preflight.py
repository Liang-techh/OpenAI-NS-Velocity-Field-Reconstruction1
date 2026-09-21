from __future__ import annotations

import math

import numpy as np

import agent7_st052m_joint_five_channel_capacity_preflight as audit


def test_source_lock_and_protocol_are_frozen() -> None:
    audit._assert_source_lock()
    assert audit.TASK_ID == "CR003-ST052M-JOINT-FIVE-CHANNEL-CAPACITY-123"
    assert audit.PREREG_ISSUE == 1031
    assert audit.SOURCE_PARENT_PR == 862
    assert audit.SOURCE_PARENT_HEAD == "f5dd27a47b03fc6239a2e6186cef78c4f3d3e024"
    assert audit.CHANNELS == (
        "amplitude",
        "radial_shape",
        "axial_turnover",
        "temporal_curvature",
        "toroidal_swirl",
    )
    assert audit.RANK_TARGET == 5
    assert audit.COSINE_MAX_ABS == 0.995
    assert audit.CONDITION_MAX == 25.0


def test_joint_identifiability_receipt_is_finite_and_complete() -> None:
    receipt = audit.joint_identifiability()
    assert receipt["channel_order"] == list(audit.CHANNELS)
    assert receipt["response_times"] == [0.375, 0.5, 0.625, 0.75]
    assert 1 <= receipt["normalized_column_rank"] <= 5
    assert len(receipt["singular_values"]) == 5
    assert len(receipt["pairwise_cosines"]) == 10
    assert set(receipt["column_norms"]) == set(audit.CHANNELS)
    assert all(math.isfinite(float(v)) and float(v) > 0.0 for v in receipt["column_norms"].values())
    assert all(math.isfinite(float(v)) for v in receipt["singular_values"])
    assert math.isfinite(float(receipt["normalized_condition"]))
    assert math.isfinite(float(receipt["max_abs_pairwise_cosine"]))
    assert isinstance(receipt["passes"], bool)


def test_structural_guards_preserve_temporal_endpoints_and_channel_selectivity() -> None:
    guards = audit.structural_guards()
    assert guards["temporal_endpoint_max_abs"] == 0.0
    assert guards["temporal_endpoints_exactly_zero"] is True
    assert guards["toroidal_rz_leakage_max_abs"] <= audit.TOROIDAL_LEAKAGE_MAX
    assert guards["toroidal_u_theta_rms"] > 1.0e-8
    assert guards["poloidal_u_theta_leakage_max_abs"] <= audit.TOROIDAL_LEAKAGE_MAX
    assert guards["toroidal_channel_selective"] is True
    assert guards["poloidal_channel_selective"] is True
    assert guards["all_pass"] is True


def test_visualization_fingerprints_are_target_free_and_toroidal_is_azimuthal() -> None:
    fingerprints = audit.visualization_fingerprint_sensitivities()
    assert fingerprints["autonomous_target_free_fingerprints"] is True
    assert fingerprints["direct_candidate_change"] is False
    assert fingerprints["direct_fingerprint_improvement"] == 0.0
    channels = fingerprints["channels"]
    assert set(channels) == set(audit.CHANNELS)

    tor = channels["toroidal_swirl"]["normalized_component_energy_fractions"]
    assert tor["u_theta"] > 1.0 - 1.0e-12
    assert tor["u_r"] <= 1.0e-12
    assert tor["u_z"] <= 1.0e-12

    for name in ("amplitude", "radial_shape", "axial_turnover", "temporal_curvature"):
        frac = channels[name]["normalized_component_energy_fractions"]
        assert frac["u_theta"] <= 1.0e-12
        assert np.isclose(frac["u_r"] + frac["u_theta"] + frac["u_z"], 1.0, rtol=1e-12, atol=1e-12)


def test_report_keeps_truth_boundary_fail_closed() -> None:
    report = audit.build_report()
    assert report["frozen_protocol"]["new_basis_dimension"] == 0
    assert report["frozen_protocol"]["coefficient_selected"] is None
    assert report["decision"]["additional_basis_dimension_justified_by_capacity_only"] is False
    assert report["decision"]["actual_velocity_changed"] is False
    assert report["decision"]["direct_visualization_fingerprint_improvement"] == 0.0
    assert report["held_out_pde_residual"]["evaluated"] is False
    assert report["held_out_pde_residual"]["st006_comparison_performed"] is False
    assert report["candidate_velocity_changed"] is False
    assert report["basis_dimension_changed"] is False
    assert report["pressure_or_force_changed"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
