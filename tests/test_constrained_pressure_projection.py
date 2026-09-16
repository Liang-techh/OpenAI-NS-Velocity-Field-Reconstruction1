from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pytest

from openai_ns_reconstruction.constrained_pressure_projection import (
    apply_projected_pressure,
    pressure_gradient_matrix,
    project_pressure_coefficients,
)


@dataclass(frozen=True)
class LinearPressureToy:
    pressure_coefficients: tuple[float, float] = (0.0, 0.0)

    def velocity(self, points, time):
        x = np.asarray(points, dtype=float)
        return np.zeros_like(x)

    def pressure_basis(self, points, time):
        x = np.asarray(points, dtype=float)
        return np.column_stack((x[:, 0], x[:, 1]))

    def pressure(self, points, time):
        return self.pressure_basis(points, time) @ np.asarray(self.pressure_coefficients)


def constant_force(vector):
    value = np.asarray(vector, dtype=float)
    def force(points, time):
        return np.broadcast_to(value, np.asarray(points).shape).copy()
    return force


def test_exact_linear_pressure_projection_recovers_manufactured_coefficients():
    rng = np.random.default_rng(23017)
    x = rng.uniform(-0.8, 0.8, (64, 3))
    t = rng.uniform(0.3, 0.7, len(x))
    target = np.array([0.3, -0.2])
    result = project_pressure_coefficients(
        LinearPressureToy(), constant_force([target[0], target[1], 0.0]), x, t, 0.01,
        lower=-1.0, upper=1.0,
    )
    assert np.allclose(result.coefficients, target, atol=2e-10)
    assert result.rank == 2
    assert result.residual_rms_after < 1e-10
    assert result.residual_rms_after < result.residual_rms_before * 1e-8
    assert result.active_lower == result.active_upper == 0


def test_projection_respects_bounds_and_fail_closed_inputs():
    x = np.array([[0.2, -0.1, 0.3], [-0.4, 0.5, -0.2], [0.7, 0.6, 0.1]])
    t = np.array([0.35, 0.5, 0.65])
    result = project_pressure_coefficients(
        LinearPressureToy(), constant_force([1.5, 0.0, 0.0]), x, t, 0.01,
        lower=-1.0, upper=1.0,
    )
    assert result.coefficients[0] == pytest.approx(1.0, abs=1e-8)
    assert result.active_upper == 1
    assert result.residual_rms_after < result.residual_rms_before
    with pytest.raises(ValueError):
        pressure_gradient_matrix(LinearPressureToy(), np.zeros((3, 2)), t)
    with pytest.raises(ValueError):
        project_pressure_coefficients(LinearPressureToy((2.0, 0.0)), constant_force([0, 0, 0]), x, t, 0.01)
    with pytest.raises(ValueError):
        project_pressure_coefficients(LinearPressureToy(), constant_force([0, 0, 0]), x, t, -0.01)


def test_current_poloidal_anchor_fixed_velocity_projection_never_worsens_l2_training_residual():
    artifact = Path('artifacts/constrained/poloidal_anchor/candidate.json')
    if not artifact.exists():
        pytest.skip('active constrained artifact is not present in this checkout')
    from openai_ns_reconstruction.constrained_poloidal import PoloidalCandidate

    candidate = PoloidalCandidate.load(artifact)
    rng = np.random.default_rng(62051)
    x = rng.uniform(-1.5, 1.5, (48, 3))
    t = rng.uniform(0.28, 0.72, len(x))
    result = project_pressure_coefficients(candidate, candidate.force, x, t, 0.01)
    projected = apply_projected_pressure(candidate, result)

    assert len(result.coefficients) == 27
    assert result.rank > 0
    assert result.residual_rms_after <= result.residual_rms_before * (1 + 1e-10)
    assert np.max(np.abs(result.coefficients)) <= 1.0 + 1e-12
    assert np.allclose(
        candidate.velocity(x[:8], t[:8]), projected.velocity(x[:8], t[:8]), atol=0.0, rtol=0.0
    )
    assert candidate.force is projected.force
