import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_temporal_derivative import (
    AUDIT_TIMES,
    PROBE_COUNT,
    PROBE_SEED,
    REFERENCE_STEP,
    TIME_LEVELS,
    _manufactured_sign_mutation_calibration,
    audit_reloaded_candidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)


@pytest.fixture(scope="module")
def report(tmp_path_factory):
    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    candidate = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(base=base)
    path = tmp_path_factory.mktemp("cubic-ut") / "candidate.json"
    candidate.save_json(path)
    loaded = Eq45SupportedPhi10CubicLocalizedTemporalCandidate.load_json(path)
    assert loaded.sha256 == candidate.sha256
    return audit_reloaded_candidate(loaded)


def test_three_level_public_time_derivative_converges_against_independent_reference(report):
    assert report["probe_seed"] == PROBE_SEED
    assert report["probe_count"] == PROBE_COUNT
    assert report["audit_times"] == list(AUDIT_TIMES)
    assert report["time_levels"] == list(TIME_LEVELS)
    assert report["probe_region_counts"] == {
        "axial_collar": 6,
        "corner_collar": 4,
        "plateau": 8,
        "radial_collar": 6,
    }

    for key in ("cubic_convergence", "quadratic_convergence", "base_convergence"):
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


def test_snapshot_identity_does_not_hide_temporal_derivative_difference(report):
    comparisons = report["same_snapshot_dynamic_comparisons"]
    expected_identity_pairs = {
        (0.25, "cubic", "quadratic"),
        (0.50, "cubic", "base"),
        (0.50, "quadratic", "base"),
        (0.625, "cubic", "base"),
        (0.75, "cubic", "base"),
        (0.75, "quadratic", "base"),
    }
    seen = set()
    for row in comparisons:
        key = (row["time"], row["left"], row["right"])
        seen.add(key)
        assert row["instantaneous_velocity_rms_difference"] == 0.0
        assert row["instantaneous_velocity_max_vector_difference"] == 0.0
        assert row["time_derivative_rms_difference"] > 1e-6
        assert np.isfinite(row["time_derivative_difference_over_right_rms"])
    assert seen == expected_identity_pairs

    cubic_at_0625 = next(
        row
        for row in comparisons
        if row["time"] == 0.625 and row["left"] == "cubic" and row["right"] == "base"
    )
    assert cubic_at_0625["time_derivative_rms_difference"] > 1e-5


def test_manufactured_sign_mutation_is_detected():
    calibration = _manufactured_sign_mutation_calibration()
    assert calibration["good_sampled_rms_error"] < 1e-10
    assert calibration["mutated_sampled_rms_error"] > 3.0
    assert calibration["mutated_sampled_max_vector_error"] > 5.0


def test_truth_boundary_keeps_temporal_audit_out_of_formal_pde_gate(report):
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
