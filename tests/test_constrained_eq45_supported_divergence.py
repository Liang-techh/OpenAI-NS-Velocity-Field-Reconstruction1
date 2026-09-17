import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_eq45_supported_candidate import (
    Eq45SupportedVelocityCandidate,
)
from openai_ns_reconstruction.constrained_eq45_supported_divergence import (
    CLAIM_SCOPE,
    PublicDivergenceMutation,
    audit_public_velocity_divergence,
    audit_supported_candidate_artifact_divergence,
)


ROOT = Path(__file__).resolve().parents[1]
PARENT_SEED = ROOT / "artifacts" / "constrained" / "eq45_velocity_candidate_seed.json"
CONSTRAINTS = ROOT / "configs" / "constraints.json"


def _probes():
    # Every point stays away from the r=1.6/2.0 and |z|=1.6/2.0 interfaces by
    # more than the coarsest registered derivative step (0.02).
    points = np.array(
        [
            # untouched plateau
            [0.43, 0.27, 0.31],
            [-0.62, 0.18, -0.47],
            [0.71, -0.39, 0.55],
            [-0.38, -0.64, -0.72],
            # radial taper collar, axial plateau
            [1.71, 0.19, 0.35],
            [-1.77, 0.16, -0.42],
            [0.24, 1.83, 0.58],
            [-0.31, -1.74, -0.66],
            # axial taper collar, radial plateau
            [0.48, 0.21, 1.72],
            [-0.55, 0.32, -1.81],
            [0.73, -0.24, 1.86],
            [-0.69, -0.34, -1.69],
            # simultaneous radial + axial taper collar
            [1.70, 0.18, 1.71],
            [-1.76, 0.12, -1.74],
            [0.17, 1.82, 1.68],
            [-0.23, -1.73, -1.84],
        ],
        dtype=float,
    )
    times = np.array(
        [
            0.25, 0.50, 0.75, 0.50,
            0.25, 0.50, 0.75, 0.50,
            0.25, 0.50, 0.75, 0.50,
            0.25, 0.50, 0.75, 0.50,
        ],
        dtype=float,
    )
    labels = np.array(
        ["plateau"] * 4
        + ["radial_collar"] * 4
        + ["axial_collar"] * 4
        + ["corner_collar"] * 4,
        dtype=object,
    )
    return points, times, labels


def _roundtrip_child(tmp_path):
    parent = Eq45VelocityCandidate.load_json(PARENT_SEED)
    child = Eq45SupportedVelocityCandidate(parent=parent)
    path = tmp_path / "eq45_supported_candidate.json"
    child.save_json(path)
    loaded = Eq45SupportedVelocityCandidate.load_json(path)
    assert loaded.sha256 == child.sha256
    return path, loaded


def test_supported_child_three_level_divergence_revalidation(tmp_path, capsys):
    artifact, child = _roundtrip_child(tmp_path)
    points, times, labels = _probes()

    report = audit_supported_candidate_artifact_divergence(
        artifact,
        CONSTRAINTS,
        points,
        times,
        labels,
    )

    assert report["claim_scope"] == CLAIM_SCOPE
    assert report["artifact_reloaded"] is True
    assert report["candidate_sha256"] == child.sha256
    assert report["parent_sha256"] == child.parent_sha256
    assert report["velocity_access"] == "at_points_only"
    assert [level["step"] for level in report["levels"]] == [0.02, 0.01, 0.005]
    assert report["registered_thresholds"] == {
        "divergence_max": 1e-5,
        "divergence_L2": 1e-5,
    }
    assert report["regions"] == [
        "plateau",
        "radial_collar",
        "axial_collar",
        "corner_collar",
    ]

    for level in report["levels"]:
        assert np.isfinite(level["max_abs"])
        assert np.isfinite(level["rms"])
        assert np.isfinite(level["normalized_rms"])
        assert set(level["by_region"]) == set(report["regions"])
        for metrics in level["by_region"].values():
            assert np.isfinite(metrics["max_abs"])
            assert np.isfinite(metrics["rms"])
            assert np.isfinite(metrics["normalized_rms"])

    # Preserve the preregistered numeric thresholds without pretending that a
    # 16-probe RMS is the registered volume-weighted spatial L2 or global max.
    finest = report["levels"][-1]
    assert report["sampled_threshold_indicators"] == {
        "sampled_max_below_registered_max": finest["max_abs"] <= 1e-5,
        "sampled_rms_below_registered_L2_number": finest["rms"] <= 1e-5,
    }
    assert report["sampled_rms_only"] is True
    assert report["global_domain_max_assessed"] is False
    assert report["volume_weighted_L2_assessed"] is False
    assert report["cr001_divergence_gate_assessed"] is False
    assert report["full_momentum_residual_assessed"] is False
    assert report["physical_support_validated"] is False
    assert report["pde_validated"] is False
    assert report["visualization_ready"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False

    # Print deterministic production numerics into exact-head Actions logs.
    with capsys.disabled():
        print("CR006_SUPPORTED_DIVERGENCE_REPORT=" + json.dumps(report, sort_keys=True))


def test_supported_child_divergence_mutation_is_detected(tmp_path, capsys):
    _, child = _roundtrip_child(tmp_path)
    points, times, labels = _probes()
    pristine = audit_public_velocity_divergence(
        child,
        points,
        times,
        labels,
        steps=(0.02, 0.01, 0.005),
    )
    mutated = audit_public_velocity_divergence(
        PublicDivergenceMutation(child, epsilon=0.03),
        points,
        times,
        labels,
        steps=(0.02, 0.01, 0.005),
    )

    pristine_finest = pristine["levels"][-1]["rms"]
    mutated_finest = mutated["levels"][-1]["rms"]
    assert mutated_finest > 0.02
    assert mutated_finest > 10.0 * max(pristine_finest, 1e-12)
    assert abs(mutated_finest - 0.03) < 0.01

    with capsys.disabled():
        print(
            "CR006_SUPPORTED_DIVERGENCE_MUTATION="
            + json.dumps(
                {
                    "epsilon": 0.03,
                    "pristine_finest_rms": pristine_finest,
                    "mutated_finest_rms": mutated_finest,
                    "detected": True,
                },
                sort_keys=True,
            )
        )


def test_supported_divergence_audit_fails_closed(tmp_path):
    _, child = _roundtrip_child(tmp_path)
    points, times, labels = _probes()

    with pytest.raises(ValueError, match="at least three"):
        audit_public_velocity_divergence(
            child, points, times, labels, steps=(0.02, 0.01)
        )
    with pytest.raises(ValueError, match="strictly decreasing"):
        audit_public_velocity_divergence(
            child, points, times, labels, steps=(0.02, 0.01, 0.01)
        )
    with pytest.raises(ValueError, match="one region label"):
        audit_public_velocity_divergence(
            child, points, times, labels[:-1], steps=(0.02, 0.01, 0.005)
        )
    with pytest.raises(ValueError, match="epsilon"):
        PublicDivergenceMutation(child, epsilon=0.0)
