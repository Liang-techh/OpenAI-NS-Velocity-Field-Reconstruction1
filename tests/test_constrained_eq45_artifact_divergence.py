from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_artifact_divergence import (
    audit_eq45_candidate_artifact_divergence,
    audit_public_velocity_divergence,
)
from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate


ROOT = Path(__file__).resolve().parents[1]
SEED_ARTIFACT = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"


def _heldout_points_times(seed=914103, count=128):
    rng = np.random.default_rng(seed)
    points = np.column_stack(
        (
            rng.uniform(-0.40, 0.40, count),
            rng.uniform(-0.40, 0.40, count),
            rng.uniform(-0.16, 0.16, count),
        )
    )
    times = rng.uniform(0.31, 0.69, count)
    return points, times


def test_reloaded_eq45_artifact_divergence_contracts_across_three_steps():
    points, times = _heldout_points_times()
    report = audit_eq45_candidate_artifact_divergence(
        SEED_ARTIFACT,
        points,
        times,
        steps=(0.06, 0.03, 0.015),
    )

    levels = report["levels"]
    rms = [level["rms"] for level in levels]
    maximum = [level["max_abs"] for level in levels]

    assert report["artifact_reloaded"] is True
    assert report["velocity_access"] == "at_points_only"
    assert report["point_count"] == 128
    assert report["candidate_sha256"] == Eq45VelocityCandidate.load_json(SEED_ARTIFACT).sha256
    assert rms[0] > rms[1] > rms[2]
    assert maximum[0] > maximum[1] > maximum[2]
    assert all(order is not None and order > 1.5 for order in report["observed_rms_orders"])
    assert levels[-1]["normalized_rms"] < levels[0]["normalized_rms"]
    assert report["full_momentum_residual_assessed"] is False
    assert report["physical_support_validated"] is False
    assert report["pde_validated"] is False


class _DivergenceMutation:
    def __init__(self, base, epsilon):
        self.base = base
        self.epsilon = float(epsilon)

    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        values = np.array(self.base.at_points(points, time), dtype=float, copy=True)
        values[:, 0] += self.epsilon * points[:, 0]
        return values


def test_public_velocity_mutation_produces_noncontracting_divergence_floor():
    points, times = _heldout_points_times(seed=914107, count=96)
    base = Eq45VelocityCandidate.load_json(SEED_ARTIFACT)
    epsilon = 0.125
    mutated = _DivergenceMutation(base, epsilon)

    report = audit_public_velocity_divergence(
        mutated,
        points,
        times,
        steps=(0.06, 0.03, 0.015),
    )
    rms = [level["rms"] for level in report["levels"]]

    assert rms[-1] > 0.10
    assert abs(rms[-1] - epsilon) < 0.02
    assert rms[-1] / rms[0] > 0.8
    assert report["pde_validated"] is False


def test_audit_uses_only_public_at_points_interface():
    class PublicOnlyRotation:
        @property
        def candidate(self):
            raise AssertionError("candidate internals must not be touched")

        def at_points(self, points, time):
            points = np.asarray(points, dtype=float)
            return np.column_stack((-points[:, 1], points[:, 0], np.zeros(points.shape[0])))

    points, times = _heldout_points_times(seed=914109, count=32)
    report = audit_public_velocity_divergence(
        PublicOnlyRotation(), points, times, steps=(0.08, 0.04, 0.02)
    )
    assert report["levels"][-1]["max_abs"] < 1e-12
    assert report["levels"][-1]["rms"] < 1e-12


def test_divergence_audit_fails_closed_on_invalid_inputs():
    points, times = _heldout_points_times(seed=914111, count=8)
    field = Eq45VelocityCandidate.load_json(SEED_ARTIFACT)

    with pytest.raises(ValueError, match="at least three"):
        audit_public_velocity_divergence(field, points, times, steps=(0.04, 0.02))
    with pytest.raises(ValueError, match="strictly decreasing"):
        audit_public_velocity_divergence(field, points, times, steps=(0.02, 0.04, 0.01))
    with pytest.raises(ValueError, match="shape"):
        audit_public_velocity_divergence(field, points[:, :2], times)

    class BadField:
        def at_points(self, points, time):
            return np.full((len(points), 3), np.nan)

    with pytest.raises(ValueError, match="finite"):
        audit_public_velocity_divergence(BadField(), points, times)
