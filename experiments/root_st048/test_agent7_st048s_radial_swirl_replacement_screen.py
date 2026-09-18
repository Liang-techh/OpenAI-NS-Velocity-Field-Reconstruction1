from __future__ import annotations

import numpy as np
import pytest

import agent7_st048s_radial_swirl_replacement_screen as screen


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


def test_ring_profile_is_compact_smooth_shape():
    radii = np.array([0.0, 0.6, 0.9, 1.2, np.sqrt(2.0), 2.0, 2.2])
    values = screen.radial_ring_profile(radii)
    assert values[0] == 0.0
    assert values[-2] == 0.0
    assert values[-1] == 0.0
    assert np.all(values[1:4] > 0.0)
    assert values[1] < values[2] < values[3] < 1.0
    assert values[4] == pytest.approx(1.0, abs=2e-15)


def test_radial_swirl_gain_changes_only_azimuthal_component():
    points = np.array([
        [0.6, 0.0, 0.3],
        [0.0, 0.9, -0.3],
        [1.2 / np.sqrt(2.0), 1.2 / np.sqrt(2.0), 0.1],
    ])
    parent = _SyntheticAxisymmetric()
    base = parent.at_points(points, 0.5)
    child = screen.RadialSwirlGainField(parent, 0.1).at_points(points, 0.5)
    ur0, ut0, uz0 = _cylindrical(points, base)
    ur1, ut1, uz1 = _cylindrical(points, child)
    expected = 1.0 + 0.1 * screen.radial_ring_profile(np.hypot(points[:, 0], points[:, 1]))
    assert np.allclose(ur1, ur0, rtol=0.0, atol=1e-14)
    assert np.allclose(uz1, uz0, rtol=0.0, atol=1e-14)
    assert np.allclose(ut1, expected * ut0, rtol=0.0, atol=1e-14)


def test_grid_truth_real_smoke_and_sensitivity_rank():
    assert screen._validate_ring_gains(screen.DEFAULT_RING_GAINS) == screen.DEFAULT_RING_GAINS
    with pytest.raises(ValueError):
        screen._validate_ring_gain(-1e-3)
    with pytest.raises(ValueError):
        screen._validate_ring_gain(0.201)
    with pytest.raises(ValueError):
        screen._validate_ring_gains((0.0, 0.1, 0.075))
    assert screen.TRUTH_BOUNDARY["canonical_velocity_changed"] is False
    assert screen.TRUTH_BOUNDARY["held_out_pde_residual_evaluated"] is False
    assert screen.TRUTH_BOUNDARY["material_path_integration_performed"] is False
    assert screen.TRUTH_BOUNDARY["visual_correspondence_verified"] is False
    assert screen.TRUTH_BOUNDARY["pde_validated"] is False

    report = screen.audit_radial_swirl_replacement_screen(
        ring_gains=(0.0, 0.10),
        energy_orders=(16, 24),
        grid_size=17,
    )
    assert report["parent_id"] == screen.PARENT_ID
    assert len(report["rows"]) == 2
    assert report["rows"][0]["preflight_pass"]
    assert report["rows"][1]["preflight_pass"]
    assert min(report["rows"][1]["mean_angular_rate_relative_changes_vs_kappa0"]) > 0.0
    sensitivity = report["four_column_velocity_sensitivity"]
    assert sensitivity["parameter_count_in_diagnostic"] == 4
    assert sensitivity["production_recommendation_parameter_count"] == 3
    assert sensitivity["normalized_rank"] == 4
    assert np.isfinite(sensitivity["normalized_condition_number"])
    assert sensitivity["normalized_condition_number"] < 100.0
    assert sensitivity["ring_column_novelty_fraction_outside_beta_gamma_global_span"] > 1e-3
