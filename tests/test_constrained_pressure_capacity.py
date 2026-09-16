from dataclasses import dataclass

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_pressure_capacity import (
    diagnose_pressure_capacity,
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
        return (
            self.pressure_basis(points, time)
            @ np.asarray(self.pressure_coefficients)
        )


def constant_force(vector):
    value = np.asarray(vector, dtype=float)

    def force(points, time):
        return np.broadcast_to(value, np.asarray(points).shape).copy()

    return force


def time_shift_force(points, time):
    x = np.asarray(points, dtype=float)
    t = np.asarray(time, dtype=float)
    coefficient = 0.2 + t
    if coefficient.ndim == 0:
        first = np.full(len(x), float(coefficient))
    else:
        first = np.broadcast_to(coefficient, (len(x),))
    return np.column_stack((first, np.zeros(len(x)), np.zeros(len(x))))


def samples(seed, n, time):
    rng = np.random.default_rng(seed)
    return rng.uniform(-0.8, 0.8, (n, 3)), np.full(n, time)


def test_representable_pressure_has_no_capacity_gap():
    train_x, train_t = samples(101, 48, 0.4)
    holdout_x, holdout_t = samples(102, 48, 0.6)
    result = diagnose_pressure_capacity(
        LinearPressureToy(),
        constant_force([0.3, -0.2, 0.0]),
        train_x,
        train_t,
        holdout_x,
        holdout_t,
        0.01,
        training_seed=101,
        holdout_seed=102,
    )
    assert np.allclose(result.training_coefficients, [0.3, -0.2], atol=2e-10)
    assert np.allclose(result.capacity_coefficients, [0.3, -0.2], atol=2e-10)
    assert result.frozen_holdout_rms < 1e-9
    assert result.capacity_holdout_rms < 1e-9
    assert result.training_rank == 2
    assert result.holdout_rank == 2
    assert result.velocity_max_change == 0.0
    assert result.capacity_fraction_of_frozen is None
    assert result.recoverable_fraction_of_frozen is None
    assert "not independent validation" in result.truth_boundary


def test_unrepresentable_channel_persists_at_capacity_ceiling():
    train_x, train_t = samples(111, 40, 0.4)
    holdout_x, holdout_t = samples(112, 40, 0.6)
    result = diagnose_pressure_capacity(
        LinearPressureToy(),
        constant_force([0.0, 0.0, 0.4]),
        train_x,
        train_t,
        holdout_x,
        holdout_t,
        0.01,
        training_seed=111,
        holdout_seed=112,
    )
    assert result.frozen_holdout_rms == pytest.approx(0.4, abs=1e-10)
    assert result.capacity_holdout_rms == pytest.approx(0.4, abs=1e-10)
    assert result.capacity_fraction_of_frozen == pytest.approx(1.0, abs=1e-10)
    assert result.recoverable_fraction_of_frozen == pytest.approx(0.0, abs=1e-10)


def test_capacity_probe_exposes_coefficient_generalization_gap():
    train_x, train_t = samples(121, 64, 0.35)
    holdout_x, holdout_t = samples(122, 64, 0.65)
    result = diagnose_pressure_capacity(
        LinearPressureToy(),
        time_shift_force,
        train_x,
        train_t,
        holdout_x,
        holdout_t,
        0.01,
        training_seed=121,
        holdout_seed=122,
    )
    # Training sees coefficient 0.55; held-out sees 0.85.
    assert result.training_coefficients[0] == pytest.approx(0.55, abs=2e-10)
    assert result.capacity_coefficients[0] == pytest.approx(0.85, abs=2e-10)
    assert result.frozen_holdout_rms == pytest.approx(0.30, abs=2e-9)
    assert result.capacity_holdout_rms < 1e-9
    assert result.recoverable_fraction_of_frozen > 0.999999
    assert result.velocity_max_change == 0.0


def test_fail_closed_seed_reuse_and_bad_bounds():
    x, t = samples(131, 8, 0.5)
    with pytest.raises(ValueError, match="seeds"):
        diagnose_pressure_capacity(
            LinearPressureToy(), constant_force([0, 0, 0]),
            x, t, x, t, 0.01, training_seed=3, holdout_seed=3
        )
    with pytest.raises(ValueError, match="invalid"):
        diagnose_pressure_capacity(
            LinearPressureToy(), constant_force([0, 0, 0]),
            x, t, x, t, 0.01, training_seed=3, holdout_seed=4,
            lower=1.0, upper=-1.0
        )
