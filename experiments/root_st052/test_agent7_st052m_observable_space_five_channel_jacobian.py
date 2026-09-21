from __future__ import annotations

import math

import numpy as np

import agent7_st052m_observable_space_five_channel_jacobian as audit


def test_source_lock_and_preregistered_protocol_are_frozen() -> None:
    audit._assert_source_lock()
    assert audit.TASK_ID == "CR003-ST052M-OBSERVABLE-SPACE-FIVE-CHANNEL-JACOBIAN-127"
    assert audit.PREREG_ISSUE == 1067
    assert audit.SOURCE_PARENT_PR == 1060
    assert audit.SOURCE_PARENT_HEAD == "637a86d84928c543526d0a7068b699ca90fe38a8"
    assert audit.TIMES == (0.25, 0.50, 0.75)
    assert audit.RING_RADII == (0.40, 0.70, 1.00, 1.30)
    assert audit.RING_ABS_Z == (0.30, 0.80)
    assert audit.AZIMUTH_COUNT == 16
    assert audit.VORTICITY_GRID == 25
    assert audit.EPSILONS == (1.0e-3, 5.0e-4)
    assert audit.PRIMARY_EPSILON == 5.0e-4
    assert audit.RANK_TARGET == 5
    assert audit.CONDITION_MAX == 25.0
    assert audit.COSINE_MAX_ABS == 0.995
    assert audit.DERIVATIVE_DRIFT_MAX == 0.05
    assert audit.DERIVATIVE_COSINE_MIN == 0.999


def test_observable_vector_covers_streamline_core_and_vorticity_families() -> None:
    base_fn, _, _ = audit._build_field_and_tangents()
    obs = audit.observable_vector(base_fn)
    assert len(obs) == 163
    assert all(math.isfinite(float(v)) for v in obs.values())
    groups = {audit._group_for_name(name) for name in obs}
    assert groups == set(audit.OBSERVABLE_GROUPS)
    assert sum(name.startswith("radial_ring/") for name in obs) == 48
    assert sum(name.startswith("swirl_ring/") for name in obs) == 48
    assert sum(name.startswith("axial_ring/") for name in obs) == 48
    assert sum(name.startswith("core_speed_width/") for name in obs) == 12
    assert sum(name.startswith("vorticity_shape/") for name in obs) == 7


def test_existing_tangents_preserve_frozen_endpoint_semantics() -> None:
    guards = audit.structural_endpoint_guards()
    assert guards["all_current_channels_exact_start_identity"] is True
    assert guards["temporal_curvature_exact_end_identity"] is True
    assert max(guards["start_tangent_max_abs_by_channel"].values()) <= 1.0e-12
    assert guards["temporal_curvature_end_max_abs"] <= 1.0e-12


def test_observable_space_jacobian_receipts_scientific_failure_instead_of_hiding_it() -> None:
    receipt = audit.observable_space_jacobian()
    assert receipt["observable_count"] == 163
    assert receipt["observable_count_by_group"] == {
        "radial_ring": 48,
        "swirl_ring": 48,
        "axial_ring": 48,
        "core_speed_width": 12,
        "vorticity_shape": 7,
    }
    assert receipt["rank_target"] == 5
    assert 0 <= receipt["normalized_column_rank"] <= 5
    assert isinstance(receipt["capacity_gate_passes"], bool)
    assert isinstance(receipt["all_derivative_stability_gates_pass"], bool)
    if receipt["normalized_condition"] is not None:
        assert math.isfinite(float(receipt["normalized_condition"]))
    if receipt["max_abs_pairwise_cosine"] is not None:
        assert math.isfinite(float(receipt["max_abs_pairwise_cosine"]))

    for channel in audit._joint().CHANNELS:
        row = receipt["channels"][channel]
        assert math.isfinite(float(row["scaled_derivative_l2"]))
        fractions = row["observable_group_energy_fractions"]
        assert set(fractions) == set(audit.OBSERVABLE_GROUPS)
        assert all(math.isfinite(float(v)) and float(v) >= 0.0 for v in fractions.values())
        if row["scaled_derivative_l2"] > np.finfo(float).tiny:
            assert abs(sum(fractions.values()) - 1.0) <= 2.0e-12
        assert len(row["dominant_observables"]) == 8

        stability = receipt["derivative_stability"][channel]
        assert stability["coarse_epsilon"] == 1.0e-3
        assert stability["fine_epsilon"] == 5.0e-4
        assert isinstance(stability["passes"], bool)
        if stability["relative_column_drift"] is not None:
            assert math.isfinite(float(stability["relative_column_drift"]))
        if stability["coarse_fine_cosine"] is not None:
            assert math.isfinite(float(stability["coarse_fine_cosine"]))


def test_report_keeps_basis_growth_delivery_and_pde_truth_fail_closed() -> None:
    report = audit.build_report()
    protocol = report["frozen_protocol"]
    assert protocol["new_basis_dimension"] == 0
    assert protocol["coefficient_selected"] is None
    assert protocol["public_openai_numeric_target"] is None
    assert report["decision"]["additional_basis_dimension_justified_by_this_audit"] is False
    assert report["decision"]["actual_velocity_changed"] is False
    assert report["decision"]["direct_visualization_fingerprint_improvement"] == 0.0
    assert report["decision"]["closer_visualization_delivery_established"] is False
    assert report["held_out_pde_residual"]["evaluated"] is False
    assert report["held_out_pde_residual"]["st006_comparison_performed"] is False
    assert report["candidate_velocity_changed"] is False
    assert report["basis_dimension_changed"] is False
    assert report["pressure_or_force_changed"] is False
    assert report["public_image_numeric_target_used"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False
