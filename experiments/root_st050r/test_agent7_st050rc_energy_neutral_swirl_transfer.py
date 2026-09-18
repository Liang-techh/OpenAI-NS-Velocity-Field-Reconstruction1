from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import agent7_st050rc_energy_neutral_swirl_transfer as screen


def test_preregistered_contract_is_frozen():
    assert screen.PREREG_ISSUE == 457
    assert screen.PARENT_ID == "ST050R-C"
    assert screen.PARENT_HEAD == "dc1d5476ea9979566b774293add76196083288c5"
    assert screen.GAINS == (0.0, 0.015, 0.025, 0.035, 0.05)
    assert screen.INNER_WINDOW == (0.30, 1.05)
    assert screen.OUTER_WINDOW == (0.95, 1.85)
    assert screen.CRITERIA == {
        "inner_gain_floor": 0.02,
        "mid_gain_floor": 0.01,
        "outer_gain_ceiling": 0.005,
        "normalization_deviation_max": 0.005,
        "inward_speed_loss_max": 0.005,
        "axial_rms_abs_change_max": 0.02,
        "radial_rms_growth_max": 0.03,
        "collar_ratio_max": 1.25,
    }


def test_compact_bump_derivative_matches_finite_difference():
    r = np.array([0.41, 0.63, 0.82, 1.08, 1.36, 1.68])
    eps = 1e-6
    for window in (screen.INNER_WINDOW, screen.OUTER_WINDOW):
        analytic = screen.compact_bump_derivative(r, *window)
        numeric = (screen.compact_bump(r + eps, *window) - screen.compact_bump(r - eps, *window)) / (2 * eps)
        assert np.allclose(analytic, numeric, rtol=2e-7, atol=2e-8)


def test_balance_profile_has_expected_localization_signs():
    # The sign structure is independent of the eventual candidate-derived alpha as long as alpha>0.
    alpha = 2.0
    assert screen.h_profile(np.array([0.6]), alpha)[0] > 0
    assert screen.h_profile(np.array([1.2]), alpha)[0] < 0
    assert screen.h_profile(np.array([0.0, 2.0]), alpha).tolist() == [0.0, 0.0]
