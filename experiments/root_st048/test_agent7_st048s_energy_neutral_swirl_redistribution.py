from __future__ import annotations

import numpy as np
import pytest

import agent7_st048s_energy_neutral_swirl_redistribution as screen


class _SyntheticAxisymmetric:
    def at_points(self, points, time):
        points = np.asarray(points, dtype=float)
        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        # u_r=-.2 r, u_theta=.3 r*(1-.15 r^2), u_z=.1 z
        swirl = 0.3 * (1.0 - 0.15 * (x * x + y * y))
        return np.column_stack((
            -0.2 * x - swirl * y,
            -0.2 * y + swirl * x,
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


def test_profiles_route_inner_mid_positive_and_outer_negative_after_balance():
    parent = _SyntheticAxisymmetric()
    balance = screen._balance_coefficient(parent, order=32)
    alpha = balance["outer_balance_coefficient_alpha"]
    assert alpha > 0.0
    assert balance["balanced_moment_relative_defect"] < 1e-14
    radii = np.array([0.0, 0.6, 0.9, 1.2, 2.0])
    profile = screen.redistribution_profile(radii, alpha)
    assert profile[0] == 0.0
    assert profile[-1] == 0.0
    assert profile[1] > 0.0
    assert profile[2] > 0.0
    assert profile[3] < 0.0


def test_redistribution_changes_only_swirl_before_common_normalization():
    parent = _SyntheticAxisymmetric()
    alpha = screen._balance_coefficient(parent, order=32)["outer_balance_coefficient_alpha"]
    points = np.array([
        [0.6, 0.0, 0.3],
        [0.0, 0.9, -0.3],
        [1.2 / np.sqrt(2.0), 1.2 / np.sqrt(2.0), 0.1],
    ])
    base = parent.at_points(points, 0.5)
    child = screen.EnergyNeutralSwirlRedistributionField(parent, 0.025, alpha).at_points(points, 0.5)
    ur0, ut0, uz0 = _cylindrical(points, base)
    ur1, ut1, uz1 = _cylindrical(points, child)
    r = np.hypot(points[:, 0], points[:, 1])
    expected = 1.0 + 0.025 * screen.redistribution_profile(r, alpha)
    assert np.allclose(ur1, ur0, rtol=0.0, atol=1e-14)
    assert np.allclose(uz1, uz0, rtol=0.0, atol=1e-14)
    assert np.allclose(ut1, expected * ut0, rtol=0.0, atol=1e-14)
    assert expected[0] > 1.0 and expected[1] > 1.0 and expected[2] < 1.0


def test_grid_truth_real_smoke_and_sensitivity_rank():
    assert screen._validate_gains(screen.DEFAULT_GAINS) == screen.DEFAULT_GAINS
    with pytest.raises(ValueError):
        screen._validate_gain(-1e-3)
    with pytest.raises(ValueError):
        screen._validate_gain(0.051)
    with pytest.raises(ValueError):
        screen._validate_gains((0.0, 0.025, 0.015))
    assert screen.TRUTH_BOUNDARY["canonical_velocity_changed"] is False
    assert screen.TRUTH_BOUNDARY["global_kappa_stacked"] is False
    assert screen.TRUTH_BOUNDARY["held_out_pde_residual_evaluated"] is False
    assert screen.TRUTH_BOUNDARY["material_path_integration_performed"] is False
    assert screen.TRUTH_BOUNDARY["visual_correspondence_verified"] is False
    assert screen.TRUTH_BOUNDARY["pde_validated"] is False

    report = screen.audit_energy_neutral_swirl_redistribution(
        gains=(0.0, 0.025),
        energy_orders=(16, 24),
        grid_size=17,
    )
    assert report["parent_id"] == screen.PARENT_ID
    assert len(report["rows"]) == 2
    assert report["energy_balance_fine"]["balanced_moment_relative_defect"] < 1e-12
    assert report["profile_design"]["seed_radius_profile_values"]["0.6"] > 0.0
    assert report["profile_design"]["seed_radius_profile_values"]["0.9"] > 0.0
    assert report["profile_design"]["seed_radius_profile_values"]["1.2"] < 0.0
    assert report["rows"][0]["preflight_pass"]
    assert report["rows"][1]["preflight_pass"]
    assert min(report["rows"][1]["seed_radius_angular_rate_relative_changes"]["0.6"]) > 0.0
    assert max(report["rows"][1]["seed_radius_angular_rate_relative_changes"]["1.2"]) < 0.0
    sensitivity = report["four_column_velocity_sensitivity"]
    assert sensitivity["parameter_count_in_diagnostic"] == 4
    assert sensitivity["production_recommendation_parameter_count"] == 3
    assert sensitivity["normalized_rank"] == 4
    assert np.isfinite(sensitivity["normalized_condition_number"])
    assert sensitivity["normalized_condition_number"] < 100.0
    assert sensitivity["redistribution_column_novelty_fraction_outside_beta_gamma_global_span"] > 1e-3
