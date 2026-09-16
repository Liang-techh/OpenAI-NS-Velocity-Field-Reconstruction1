from dataclasses import dataclass

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_pressure_holdout import (
    audit_projected_pressure,
)


@dataclass(frozen=True)
class LinearPressureToy:
    pressure_coefficients: tuple[float, float] = (0.0, 0.0)

    def velocity(self, points, time):
        return np.zeros_like(np.asarray(points, dtype=float))

    def pressure_basis(self, points, time):
        x = np.asarray(points, dtype=float)
        return np.column_stack((x[:, 0], x[:, 1]))

    def pressure(self, points, time):
        return (self.pressure_basis(points, time)
                @ np.asarray(self.pressure_coefficients))


def constant_force(vector):
    value = np.asarray(vector, dtype=float)

    def force(points, time):
        return np.broadcast_to(value, np.asarray(points).shape).copy()

    return force


def test_holdout_audit_recovers_pressure_and_preserves_velocity():
    train_rng = np.random.default_rng(41)
    holdout_rng = np.random.default_rng(42)
    train_x = train_rng.uniform(-0.8, 0.8, (64, 3))
    train_t = train_rng.uniform(0.3, 0.7, 64)
    holdout_x = holdout_rng.uniform(-0.8, 0.8, (64, 3))
    holdout_t = np.repeat([0.35, 0.65], 32)

    audit = audit_projected_pressure(
        LinearPressureToy(),
        constant_force([0.3, -0.2, 0.0]),
        train_x,
        train_t,
        holdout_x,
        holdout_t,
        0.01,
        training_seed=41,
        holdout_seed=42,
    )

    assert np.allclose(audit.projected_coefficients, [0.3, -0.2], atol=2e-10)
    assert audit.training_rms_after < 1e-10
    assert audit.heldout_rms_after < 1e-9
    assert audit.heldout_rms_ratio < 1e-8
    assert audit.velocity_max_change == 0.0


def test_holdout_audit_does_not_hide_unrepresentable_force_channel():
    train_rng = np.random.default_rng(51)
    holdout_rng = np.random.default_rng(52)
    train_x = train_rng.uniform(-0.7, 0.7, (40, 3))
    train_t = train_rng.uniform(0.3, 0.7, 40)
    holdout_x = holdout_rng.uniform(-0.7, 0.7, (40, 3))
    holdout_t = np.repeat([0.4, 0.6], 20)

    audit = audit_projected_pressure(
        LinearPressureToy(),
        constant_force([0.0, 0.0, 0.4]),
        train_x,
        train_t,
        holdout_x,
        holdout_t,
        0.01,
        training_seed=51,
        holdout_seed=52,
    )

    assert audit.heldout_max_after == pytest.approx(0.4, rel=0, abs=1e-10)
    assert audit.heldout_rms_after == pytest.approx(0.4, rel=0, abs=1e-10)
    assert audit.velocity_max_change == 0.0


def test_holdout_audit_rejects_seed_reuse_and_bad_holdout_shape():
    x = np.zeros((4, 3))
    t = np.full(4, 0.5)
    force = constant_force([0.0, 0.0, 0.0])

    with pytest.raises(ValueError, match='seeds'):
        audit_projected_pressure(
            LinearPressureToy(), force, x, t, x, t, 0.01,
            training_seed=7, holdout_seed=7,
        )
    with pytest.raises(ValueError, match='times'):
        audit_projected_pressure(
            LinearPressureToy(), force, x, t, x, t[:3], 0.01,
            training_seed=7, holdout_seed=8,
        )
