import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_vorticity_residual import (
    audit_frozen_eq45,
    audit_vorticity_equation,
    vorticity_equation_residual,
)


class RigidRotation:
    def at_points(self, points, time):
        p = np.asarray(points, dtype=float)
        x, y, _ = p.T
        return np.column_stack((-y, x, np.zeros_like(x)))


class SinusoidalShear:
    def at_points(self, points, time):
        p = np.asarray(points, dtype=float)
        y = p[:, 1]
        return np.column_stack((np.sin(y), np.zeros_like(y), np.zeros_like(y)))


class ExactShearForce:
    def __init__(self, nu, sign=1.0):
        self.nu = float(nu)
        self.sign = float(sign)

    def __call__(self, points, time):
        p = np.asarray(points, dtype=float)
        y = p[:, 1]
        return np.column_stack(
            (self.sign * self.nu * np.sin(y), np.zeros_like(y), np.zeros_like(y))
        )


def test_rigid_rotation_has_zero_pressure_free_vorticity_residual():
    points = np.array(
        [[0.12, -0.21, 0.07], [-0.18, 0.09, -0.11], [0.04, 0.16, 0.19]],
        dtype=float,
    )
    residual, terms = vorticity_equation_residual(
        RigidRotation(),
        points,
        0.5,
        spatial_step=0.01,
        time_step=0.005,
        nu=0.01,
    )
    np.testing.assert_allclose(terms["omega"], np.array([[0.0, 0.0, 2.0]] * 3), atol=2e-12)
    assert np.max(np.linalg.norm(residual, axis=1)) < 1e-9


def test_sinusoidal_shear_checks_diffusion_and_force_curl_sign():
    nu = 0.01
    points = np.array(
        [[0.0, -0.7, 0.0], [0.1, -0.2, -0.1], [-0.1, 0.35, 0.2], [0.0, 0.8, -0.2]],
        dtype=float,
    )
    zero, _ = vorticity_equation_residual(
        SinusoidalShear(),
        points,
        0.5,
        spatial_step=0.004,
        time_step=0.003,
        nu=nu,
    )
    expected = np.column_stack(
        (np.zeros(len(points)), np.zeros(len(points)), -nu * np.cos(points[:, 1]))
    )
    np.testing.assert_allclose(zero, expected, atol=2e-6, rtol=2e-4)

    closed, _ = vorticity_equation_residual(
        SinusoidalShear(),
        points,
        0.5,
        spatial_step=0.004,
        time_step=0.003,
        nu=nu,
        force=ExactShearForce(nu, sign=1.0),
    )
    reversed_force, _ = vorticity_equation_residual(
        SinusoidalShear(),
        points,
        0.5,
        spatial_step=0.004,
        time_step=0.003,
        nu=nu,
        force=ExactShearForce(nu, sign=-1.0),
    )
    assert np.max(np.linalg.norm(closed, axis=1)) < 3e-6
    assert np.sqrt(np.mean(np.sum(reversed_force**2, axis=1))) > 1.8 * np.sqrt(
        np.mean(np.sum(zero**2, axis=1))
    )


def test_frozen_eq45_audit_uses_registered_levels_and_keeps_truth_boundary():
    report = audit_frozen_eq45(sample_count=4, seed=914131)
    assert report["schema"] == "eq45_vorticity_equation_crosscheck_v1"
    assert report["spatial_steps"] == [0.02, 0.01, 0.005]
    assert len(report["rows"]) == 3
    assert report["restricted_force_comparison"]["a"] == 0.0
    assert report["restricted_force_comparison"]["c"] == pytest.approx(2.19114231)
    assert report["finest_zero_force_rms"] > 1e-3
    assert report["finest_restricted_force_rms"] > 1e-3
    for row in report["rows"]:
        assert np.isfinite(row["zero_force_rms"])
        assert np.isfinite(row["restricted_force_rms"])
        assert row["zero_force_rms"] > 0.0
        assert row["restricted_force_rms"] > 0.0
    assert report["forcing_fitted_in_this_audit"] is False
    assert report["physical_support_validated"] is False
    assert report["visualization_ready"] is False
    assert report["pde_validated"] is False
    assert report["paper_exact"] is False
    assert report["openai_field_identified"] is False


def test_audit_fails_closed_on_insufficient_levels_and_bad_velocity_shape():
    points = np.array([[0.1, 0.1, 0.1], [-0.1, -0.1, -0.1]])
    with pytest.raises(ValueError, match="at least three"):
        audit_vorticity_equation(
            RigidRotation(), points, 0.5, spatial_steps=(0.02, 0.01), time_step=0.005
        )

    class BadField:
        def at_points(self, points, time):
            return np.zeros((len(points), 2))

    with pytest.raises(ValueError, match="shape"):
        vorticity_equation_residual(
            BadField(),
            points,
            0.5,
            spatial_step=0.01,
            time_step=0.005,
            nu=0.01,
        )
