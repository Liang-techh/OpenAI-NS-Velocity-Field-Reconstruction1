import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_temporal_derivative import (
    AUDIT_TIMES,
    PROBE_COUNT,
    PROBE_SEED,
    REFERENCE_STEP,
    STATIC_RETURN_TIMES,
    TIME_LEVELS,
    _manufactured_sign_mutation_calibration,
    audit_reloaded_candidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)


EXPECTED_QUARTIC_SHA256 = "03fdae73170469ae1160297b489b1b83a27adddfc941119e116a21b98bd15049"


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    candidate = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(base=base)
    path = tmp_path_factory.mktemp("quartic-ut") / "candidate.json"
    candidate.save_json(path)
    loaded = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate.load_json(path)
    assert loaded.sha256 == candidate.sha256 == EXPECTED_QUARTIC_SHA256
    return audit_reloaded_candidate(loaded)


def test_three_level_public_quartic_time_derivative_converges(report):
    assert report["probe_seed"] == PROBE_SEED
    assert report["probe_count"] == PROBE_COUNT
    assert report["audit_times"] == list(AUDIT_TIMES)
    assert report["time_levels"] == list(TIME_LEVELS)
    assert report["reference_step"] == REFERENCE_STEP
    assert report["probe_region_counts"] == {
        "axial_collar": 6,
        "corner_collar": 4,
        "plateau": 8,
        "radial_collar": 6,
    }

    for key in ("quartic_convergence", "cubic_convergence", "base_convergence"):
        convergence = report[key]
        assert convergence["reference_step"] == REFERENCE_STEP
        assert convergence["reference_operator"] == "fourth_order_public_velocity_time_difference"
        errors = [level["sampled_rms"] for level in convergence["levels"]]
        normalized = [level["normalized_rms"] for level in convergence["levels"]]
        assert errors[2] < errors[1] < errors[0]
        assert normalized[2] < normalized[1] < normalized[0]
        assert normalized[2] < 0.05
        assert all(1.5 < order < 2.5 for order in convergence["observed_orders"])
        assert all(np.isfinite(level["sampled_max_vector"]) for level in convergence["levels"])


def test_quartic_balancing_reduces_public_dynamic_collateral_at_static_returns(report):
    comparisons = report["anchor_dynamic_comparisons"]
    rows = comparisons["static_returns"]
    assert [row["time"] for row in rows] == list(STATIC_RETURN_TIMES)
    for row in rows:
        assert row["quartic_static_velocity_rms_difference"] == 0.0
        assert row["quartic_static_velocity_max_vector_difference"] == 0.0
        assert row["cubic_static_velocity_rms_difference"] == 0.0
        assert row["cubic_static_velocity_max_vector_difference"] == 0.0
        assert row["quartic_static_time_derivative_rms_difference"] > 1e-7
        assert row["cubic_static_time_derivative_rms_difference"] > 1e-7
        assert 0.0 < row["quartic_to_cubic_dynamic_collateral_ratio"] < 0.6

    aggregate = comparisons["aggregate_static_return_dynamic_collateral"]
    assert aggregate["quartic_rms"] > 0.0
    assert aggregate["cubic_rms"] > 0.0
    assert 0.0 < aggregate["quartic_to_cubic_ratio"] < 0.6


def test_same_early_snapshot_still_has_different_quartic_and_cubic_time_derivatives(report):
    early = report["anchor_dynamic_comparisons"]["early_same_snapshot"]
    assert early["time"] == 0.25
    assert early["quartic_cubic_velocity_rms_difference"] == 0.0
    assert early["quartic_cubic_velocity_max_vector_difference"] == 0.0
    assert early["quartic_cubic_time_derivative_rms_difference"] > 1e-5
    assert early["quartic_cubic_time_derivative_difference_over_cubic_rms"] > 0.0

    off_keyframes = report["off_keyframe_derivative_comparisons"]
    assert [row["time"] for row in off_keyframes] == [0.375, 0.5625, 0.6875]
    for row in off_keyframes:
        assert row["quartic_cubic_velocity_rms_difference"] > 0.0
        assert row["quartic_cubic_time_derivative_rms_difference"] > 0.0
        assert np.isfinite(row["quartic_static_time_derivative_rms_difference"])
        assert np.isfinite(row["cubic_static_time_derivative_rms_difference"])


def test_manufactured_sign_mutation_is_detected():
    calibration = _manufactured_sign_mutation_calibration()
    assert calibration["good_sampled_rms_error"] < 1e-10
    assert calibration["mutated_sampled_rms_error"] > 3.0
    assert calibration["mutated_sampled_max_vector_error"] > 5.0


def test_truth_boundary_keeps_quartic_derivative_audit_out_of_formal_pde_gate(report):
    truth = report["truth_boundary"]
    assert truth["serialized_candidate_reloaded"] is True
    assert truth["public_velocity_only"] is True
    assert truth["training_loss_reused"] is False
    assert truth["pressure_fitted"] is False
    assert truth["forcing_fitted"] is False
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
