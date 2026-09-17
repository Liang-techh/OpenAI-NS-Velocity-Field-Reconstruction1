import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_compact_quartic_blend_temporal_derivative import (
    AUDIT_TIMES,
    BLEND_WEIGHTS,
    PROBE_SEED,
    REFERENCE_STEP,
    TIME_LEVELS,
    audit_blend_temporal_derivative,
)


@pytest.fixture(scope="module")
def report():
    return audit_blend_temporal_derivative()


def test_three_level_public_time_derivative_converges_for_endpoints_and_midpoint(report):
    assert report["probe_seed"] == PROBE_SEED
    assert report["probe_count"] == 24
    assert report["blend_weights"] == list(BLEND_WEIGHTS)
    assert report["audit_times"] == list(AUDIT_TIMES)
    assert report["time_levels"] == list(TIME_LEVELS)
    assert report["reference_step"] == REFERENCE_STEP
    assert report["probe_region_counts"] == {
        "axial_collar": 6,
        "corner_collar": 4,
        "plateau": 8,
        "radial_collar": 6,
    }

    for key in ("0.000000", "0.500000", "1.000000"):
        convergence = report["convergence"][key]
        assert convergence["reference_step"] == REFERENCE_STEP
        assert convergence["reference_operator"] == "fourth_order_public_velocity_time_difference"
        errors = [level["sampled_rms"] for level in convergence["levels"]]
        normalized = [level["normalized_rms"] for level in convergence["levels"]]
        assert errors[2] < errors[1] < errors[0]
        assert normalized[2] < normalized[1] < normalized[0]
        assert normalized[2] < 0.05
        assert all(1.5 < order < 2.5 for order in convergence["observed_orders"])
        assert all(np.isfinite(level["sampled_max_vector"]) for level in convergence["levels"])


def test_midpoint_public_velocity_and_time_derivative_are_affine_between_endpoints(report):
    rows = report["midpoint_endpoint_affinity"]
    assert [row["time"] for row in rows] == list(AUDIT_TIMES)
    for row in rows:
        assert row["instantaneous_velocity_affine_rms_error"] < 1e-12
        assert row["instantaneous_velocity_affine_max_vector_error"] < 1e-11
        assert row["time_derivative_affine_rms_error"] < 1e-9
        assert row["time_derivative_affine_max_vector_error"] < 1e-8
        assert row["time_derivative_affine_normalized_rms_error"] < 1e-8
        assert np.isfinite(row["endpoint_time_derivative_span_rms"])

    active = [
        row for row in rows if row["time"] in (0.25, 0.3125, 0.375, 0.4375)
    ]
    assert all(row["endpoint_time_derivative_span_rms"] > 1e-8 for row in active)


def test_sign_mutation_calibration_remains_sensitive(report):
    calibration = report["mutation_calibration"]
    assert calibration["good_sampled_rms_error"] < 1e-10
    assert calibration["mutated_sampled_rms_error"] > 3.0
    assert calibration["mutated_sampled_max_vector_error"] > 5.0


def test_truth_boundary_keeps_blend_derivative_audit_out_of_pde_selection(report):
    truth = report["truth_boundary"]
    assert truth["serialized_candidates_reloaded"] is True
    assert truth["public_velocity_only"] is True
    assert truth["blend_weight_selected"] is False
    assert truth["training_loss_reused"] is False
    assert truth["pressure_fitted"] is False
    assert truth["forcing_fitted"] is False
    assert truth["vorticity_or_momentum_residual_evaluated"] is False
    assert truth["formal_pde_gate_assessed"] is False
    assert truth["visualization_candidate_only"] is True
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False


@pytest.mark.parametrize(
    "kwargs",
    [
        {"weights": (0.0, 1.0)},
        {"weights": (0.0, 0.5, 0.5, 1.0)},
        {"weights": (0.0, 0.25, 1.0)},
        {"weights": (0.0, 0.5, 1.1)},
        {"probe_seed": -1},
    ],
)
def test_invalid_audit_contract_fails_closed(kwargs):
    with pytest.raises(ValueError):
        audit_blend_temporal_derivative(**kwargs)
