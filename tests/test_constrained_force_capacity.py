import numpy as np
import pytest

from openai_ns_reconstruction.constrained_force import RestrictedForce
from openai_ns_reconstruction.constrained_force_capacity import (
    diagnose_fixed_candidate_force_capacity,
    fit_restricted_force_capacity,
)


def _points(seed=1702, n=96):
    return np.random.default_rng(seed).uniform(-1.2, 1.2, size=(n, 3))


def test_exact_recovery_inside_preregistered_force_family():
    points = _points()
    time = 0.5
    target = RestrictedForce(2.3, 1.7)(points, time)
    result = fit_restricted_force_capacity(target, points, time)
    assert result["design_rank"] == 2
    assert result["coefficients"]["a"] == pytest.approx(2.3, abs=2e-10)
    assert result["coefficients"]["c"] == pytest.approx(1.7, abs=2e-10)
    assert result["capacity_rms"] < 1e-10
    assert result["pde_validated"] is False
    assert "not independent validation" in result["scope"]


def test_orthogonal_residual_floor_is_not_erased_by_restricted_force():
    points = _points(seed=17, n=128)
    time = 0.5
    fa = RestrictedForce(1.0, 0.0)(points, time)
    fc = RestrictedForce(0.0, 1.0)(points, time)
    matrix = np.column_stack((fa.ravel(), fc.ravel()))

    raw = np.column_stack((
        np.ones(len(points)),
        np.linspace(-1.0, 1.0, len(points)),
        np.cos(np.arange(len(points))),
    )).ravel()
    raw = raw - matrix @ np.linalg.lstsq(matrix, raw, rcond=None)[0]
    raw = raw / np.sqrt(np.mean(raw * raw)) * 0.4
    floor = raw.reshape((-1, 3))

    base = RestrictedForce(0.8, 1.1)(points, time) + floor
    result = fit_restricted_force_capacity(base, points, time)
    assert result["coefficients"]["a"] == pytest.approx(0.8, abs=2e-10)
    assert result["coefficients"]["c"] == pytest.approx(1.1, abs=2e-10)
    assert result["capacity_rms"] == pytest.approx(np.sqrt(3.0) * 0.4, rel=2e-10)
    assert result["recoverable_fraction"] < 1.0


def test_preregistered_upper_bound_caps_force_capacity():
    points = _points(seed=27)
    time = 0.5
    base = 12.0 * RestrictedForce(1.0, 0.0)(points, time)
    result = fit_restricted_force_capacity(base, points, time)
    assert result["coefficients"]["a"] == pytest.approx(10.0, abs=1e-8)
    assert "a:upper" in result["active_bounds"]
    assert result["capacity_rms"] > 1e-4


def test_fixed_candidate_wrapper_uses_independent_operator_but_stays_capacity_only():
    points = _points(seed=99, n=24)

    def velocity(x, t):
        x = np.asarray(x, dtype=float)
        return np.zeros_like(x)

    def pressure(x, t):
        x = np.asarray(x, dtype=float)
        return np.zeros(x.shape[:-1], dtype=float)

    result = diagnose_fixed_candidate_force_capacity(
        velocity,
        pressure,
        points,
        0.5,
        nu=0.01,
        step=0.005,
    )
    assert result["velocity_and_pressure_frozen"] is True
    assert result["force_fit_uses_diagnostic_samples"] is True
    assert result["pde_validated"] is False
    assert result["before_rms"] == pytest.approx(0.0, abs=1e-15)
    assert result["capacity_rms"] == pytest.approx(0.0, abs=1e-15)


def test_fail_closed_inputs_and_bounds():
    points = _points(n=8)
    base = np.zeros_like(points)
    with pytest.raises(ValueError):
        fit_restricted_force_capacity(base, points, 0.5, lower=-1.0)
    with pytest.raises(ValueError):
        fit_restricted_force_capacity(base, points, 0.5, upper=11.0)
    bad = base.copy()
    bad[0, 0] = np.nan
    with pytest.raises(ValueError):
        fit_restricted_force_capacity(bad, points, 0.5)
