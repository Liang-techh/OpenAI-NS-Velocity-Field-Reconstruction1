from __future__ import annotations

import numpy as np
import pytest

import agent7_st048s_piola_swirl_gain_screen as screen


class _SyntheticAxisymmetric:
    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        # u_r=-.2 r, u_theta=.3 r, u_z=.1 z
        return np.column_stack((
            -0.2 * x - 0.3 * y,
            -0.2 * y + 0.3 * x,
            0.1 * z,
        ))


def _cylindrical(points, velocity):
    points = np.asarray(points, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    x, y = points[:, 0], points[:, 1]
    r = np.hypot(x, y)
    rx, ry = x / r, y / r
    ur = rx * velocity[:, 0] + ry * velocity[:, 1]
    ut = -ry * velocity[:, 0] + rx * velocity[:, 1]
    return ur, ut, velocity[:, 2]


def test_swirl_gain_changes_only_azimuthal_component():
    points = np.array([
        [0.6, 0.0, 0.3],
        [0.0, 0.9, -0.3],
        [0.8 / np.sqrt(2.0), 0.8 / np.sqrt(2.0), 0.1],
    ])
    parent = _SyntheticAxisymmetric()
    base = parent.at_points(points, 0.5)
    child = screen.SwirlGainField(parent, 0.1).at_points(points, 0.5)
    ur0, ut0, uz0 = _cylindrical(points, base)
    ur1, ut1, uz1 = _cylindrical(points, child)
    assert np.allclose(ur1, ur0, rtol=0.0, atol=1e-14)
    assert np.allclose(uz1, uz0, rtol=0.0, atol=1e-14)
    assert np.allclose(ut1, 1.1 * ut0, rtol=0.0, atol=1e-14)


def test_kappa_grid_and_truth_boundary_fail_closed():
    assert screen._validate_kappas(screen.DEFAULT_KAPPAS) == screen.DEFAULT_KAPPAS
    with pytest.raises(ValueError):
        screen._validate_kappa(-1e-3)
    with pytest.raises(ValueError):
        screen._validate_kappa(0.151)
    with pytest.raises(ValueError):
        screen._validate_kappas((0.0, 0.05, 0.025))
    assert screen.TRUTH_BOUNDARY["new_spatial_basis_added"] is False
    assert screen.TRUTH_BOUNDARY["held_out_pde_residual_evaluated"] is False
    assert screen.TRUTH_BOUNDARY["pde_validated"] is False
    assert screen.TRUTH_BOUNDARY["visual_correspondence_verified"] is False


def test_real_screen_smoke_and_sensitivity_rank():
    report = screen.audit_swirl_gain_screen(
        kappas=(0.0, 0.05),
        energy_orders=(16, 24),
        grid_size=17,
    )
    assert report["parent_id"] == screen.PARENT_ID
    assert report["fixed_temporal_geometry"]["gamma"] == screen.BASE_GAMMA
    assert len(report["rows"]) == 2
    assert report["rows"][0]["preflight_pass"]
    assert report["rows"][1]["preflight_pass"]
    assert min(report["rows"][1]["mean_angular_rate_relative_changes_vs_kappa0"]) > 0.0
    sensitivity = report["three_column_velocity_sensitivity"]
    assert sensitivity["parameter_count"] == 3
    assert sensitivity["normalized_rank"] == 3
    assert np.isfinite(sensitivity["normalized_condition_number"])
    assert sensitivity["normalized_condition_number"] < 50.0
