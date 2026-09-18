from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import agent7_st051b_frozen_redistribution_transfer as screen


def test_preregistered_frozen_transfer_contract():
    assert screen.PREREG_ISSUE == 468
    assert screen.PARENT_ID == "ST051-B"
    assert screen.PARENT_HEAD == "4b784f1b8457af2ead49295631d834d4e882000b"
    assert screen.SOURCE_TRANSFER_PR == 458
    assert screen.SOURCE_PARENT_ID == "ST050R-C"
    assert screen.SOURCE_ALPHA == 2.520520814687742
    assert screen.GAIN == 0.025
    assert screen.INNER_WINDOW == (0.30, 1.05)
    assert screen.OUTER_WINDOW == (0.95, 1.85)
    assert screen.GRID_SIZE == 33
    assert screen.CRITERIA == {
        "inner_gain_floor": 0.02,
        "mid_gain_floor": 0.01,
        "outer_gain_ceiling": 0.005,
        "normalization_deviation_max": 0.005,
        "inward_speed_loss_max": 0.005,
        "axial_rms_abs_change_max": 0.02,
        "radial_rms_growth_max": 0.03,
        "collar_ratio_max": 1.25,
        "support_max_abs": 1e-12,
        "divergence_fd_max": 1e-5,
    }


def test_frozen_profile_has_intended_inner_mid_outer_signs_and_support():
    vals = screen.h_profile(np.array([0.0, 0.6, 0.9, 1.2, 2.0]))
    assert vals[0] == 0.0
    assert vals[1] > 0.0
    assert vals[2] > 0.0
    assert vals[3] < 0.0
    assert vals[4] == 0.0


def test_compact_bump_derivative_matches_finite_difference():
    r = np.array([0.41, 0.63, 0.82, 1.08, 1.36, 1.68])
    eps = 1e-6
    for window in (screen.INNER_WINDOW, screen.OUTER_WINDOW):
        analytic = screen.compact_bump_derivative(r, *window)
        numeric = (screen.compact_bump(r + eps, *window) - screen.compact_bump(r - eps, *window)) / (2.0 * eps)
        assert np.allclose(analytic, numeric, rtol=2e-7, atol=2e-8)


def test_clean_rule_is_fail_closed():
    row = dict(
        angular_gain_r06=0.021,
        angular_gain_r09=0.011,
        angular_gain_r12=-0.01,
        normalization=0.999,
        inward_speed_loss=0.001,
        axial_rms_relative_changes=[0.0, 0.001, -0.001],
        radial_rms_relative_changes=[-0.01, 0.0, 0.01],
        collar_ratios=[1.0, 1.0, 1.01],
        core_signs_pass=True,
        support_max_abs=0.0,
        divergence_fd_max=1e-8,
    )
    assert screen.clean_rule(row)
    row["angular_gain_r06"] = 0.019
    assert not screen.clean_rule(row)
