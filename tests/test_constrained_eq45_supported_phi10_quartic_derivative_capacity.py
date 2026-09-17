import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_derivative_capacity import (
    audit_supported_phi10_quartic_derivative_capacity,
    balanced_nullspace_coefficient,
    quartic_phi10_delta,
    quartic_phi10_delta_dtau,
    quartic_nullspace,
)


def test_quartic_nullspace_preserves_screened_snapshots_and_adds_one_well_conditioned_degree():
    report = audit_supported_phi10_quartic_derivative_capacity()

    assert report["task_id"] == "CR003-EQ45-SUPPORTED-PHI10-QUARTIC-DERIVATIVE-CAPACITY-025"
    assert report["time_design"]["rank"] == 4
    assert report["time_design"]["basis_dimension_increment"] == 1
    assert report["time_design"]["condition_number"] < 9.0

    coefficient = report["nullspace_coefficient"]
    assert coefficient == pytest.approx(-0.2831460674157303, rel=0.0, abs=1.0e-14)
    np.testing.assert_allclose(quartic_nullspace(np.array([-1.0, 0.0, 0.5, 1.0])), 0.0, atol=0.0)
    np.testing.assert_allclose(
        quartic_phi10_delta(np.array([-1.0, 0.0, 0.5, 1.0])),
        np.array([-1.4, 0.0, 0.0, 0.0]),
        rtol=0.0,
        atol=2.0e-15,
    )

    for row in report["anchor_public_velocity_checks"]:
        assert row["public_velocity_rms_difference"] < 2.0e-13
        assert row["public_velocity_max_abs_difference"] < 1.0e-12

    extrema = report["schedule_extrema"]
    assert extrema["minimum_coefficient"] >= -extrema["coefficient_limit"]
    assert extrema["maximum_coefficient"] <= extrema["coefficient_limit"]
    assert report["truth_boundary"]["new_spatial_basis_added"] is False
    assert report["truth_boundary"]["pde_validated"] is False
    assert report["truth_boundary"]["visual_correspondence_verified"] is False


def test_balanced_quartic_reduces_return_slope_and_off_keyframe_velocity_collateral():
    report = audit_supported_phi10_quartic_derivative_capacity()

    capacity = report["return_frame_derivative_capacity"]
    assert capacity["quartic_to_cubic_l2_ratio"] == pytest.approx(
        0.29981267559834457, rel=0.0, abs=2.0e-14
    )
    assert capacity["fractional_l2_reduction"] > 0.70

    late = report["late_half_localization"]
    assert late["quartic_to_cubic_ratio"] < 0.29
    assert late["quartic_max_abs_delta"] < 0.0064

    rows = {row["time"]: row for row in report["off_keyframe_public_velocity_checks"]}
    assert rows[0.375]["quartic_to_cubic_public_velocity_delta_rms_ratio"] < 0.71
    assert rows[0.5625]["quartic_to_cubic_public_velocity_delta_rms_ratio"] < 0.25
    assert rows[0.6875]["quartic_to_cubic_public_velocity_delta_rms_ratio"] < 0.07

    early_cost = report["early_endpoint_temporal_cost"]
    assert 1.27 < early_cost["quartic_to_cubic_abs_ratio"] < 1.29
    assert report["truth_boundary"]["pde_objective_used_to_choose_nullspace_coefficient"] is False
    assert report["truth_boundary"]["production_temporal_shape_promoted"] is False


def test_balanced_coefficient_is_stationary_for_return_node_slope_objective_and_inputs_fail_closed():
    coefficient = balanced_nullspace_coefficient()
    return_tau = np.array([0.0, 0.5, 1.0])
    slopes = np.asarray(
        quartic_phi10_delta_dtau(return_tau, nullspace_coefficient=coefficient), dtype=float
    )

    step = 1.0e-6
    minus = np.asarray(
        quartic_phi10_delta_dtau(return_tau, nullspace_coefficient=coefficient - step), dtype=float
    )
    plus = np.asarray(
        quartic_phi10_delta_dtau(return_tau, nullspace_coefficient=coefficient + step), dtype=float
    )
    objective = float(np.dot(slopes, slopes))
    assert objective < float(np.dot(minus, minus))
    assert objective < float(np.dot(plus, plus))

    with pytest.raises(ValueError):
        balanced_nullspace_coefficient(early_delta=0.0)
    with pytest.raises(ValueError):
        quartic_nullspace(1.01)
    with pytest.raises(ValueError):
        quartic_phi10_delta(0.0, nullspace_coefficient=np.inf)
