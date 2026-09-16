import numpy as np
import pytest

from openai_ns_reconstruction.constrained_velocity_grid_consistency import (
    _trilinear,
    audit_velocity_grid_consistency,
    refinement_summary,
)


class SmoothField:
    def at_points(self, points, time):
        p = np.asarray(points, dtype=float)
        x, y, z = p.T
        return np.column_stack((
            x + 0.2*y + 0.1*time*np.sin(z),
            y - 0.3*z + 0.1*np.sin(x*y),
            z + 0.25*x + 0.05*time*np.cos(y),
        ))

    def grid(self, x, y, z, times):
        tt, xx, yy, zz = np.meshgrid(np.asarray(times), x, y, z, indexing='ij')
        pts = np.stack((xx, yy, zz), axis=-1)
        vals = []
        for i, time in enumerate(np.asarray(times)):
            vals.append(self.at_points(pts[i].reshape(-1, 3), float(time)).reshape(pts.shape[1:]))
        return np.stack(vals, axis=0)


def test_uniform_trilinear_interpolation_is_exact_for_affine_vector_field():
    axes = np.linspace(-1, 1, 7)
    xx, yy, zz = np.meshgrid(axes, axes, axes, indexing='ij')
    values = np.stack((xx + 2*yy, yy - 3*zz, 0.5*xx + zz), axis=-1)
    points = np.array([[0.123, -0.456, 0.789], [-0.77, 0.18, -0.31]])
    got = _trilinear(values, points, ((-1, 1), (-1, 1), (-1, 1)))
    expected = np.column_stack((points[:,0] + 2*points[:,1], points[:,1] - 3*points[:,2], 0.5*points[:,0] + points[:,2]))
    np.testing.assert_allclose(got, expected, atol=2e-15, rtol=0)


def test_three_level_offgrid_audit_detects_resolution_improvement():
    rows = audit_velocity_grid_consistency(
        SmoothField(), grid_sizes=(7, 13, 25), times=(0.35, 0.65),
        sample_count=256, seed=12345,
    )
    summary = refinement_summary(rows)
    assert len(rows) == 6
    for data in summary.values():
        assert data['coarse_to_fine_ratio'] < 0.30
        assert data['finest_relative_rms'] < 5e-4
        assert data['interpretation'].endswith('not PDE evidence')


def test_audit_fail_closes_on_invalid_grid_contract():
    with pytest.raises(ValueError):
        audit_velocity_grid_consistency(SmoothField(), grid_sizes=(9, 9, 17))


def test_packaged_velocity_field_is_audited_through_public_api():
    from openai_ns_reconstruction.velocity_components import VelocityField
    field = VelocityField()
    rows = audit_velocity_grid_consistency(
        field, grid_sizes=(5, 9, 13), times=(0.5,), sample_count=32, seed=20260916,
    )
    assert [row.grid_size for row in rows] == [5, 9, 13]
    assert all(np.isfinite(row.max_vector_error) for row in rows)
    assert all(np.isfinite(row.relative_rms_error) for row in rows)
