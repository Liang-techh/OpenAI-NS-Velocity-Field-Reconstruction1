from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import agent7_st052m_frozen_redistribution as screen


def test_preregistered_contract_is_frozen():
    assert screen.PREREG_ISSUE == 527
    assert screen.PARENT_ID == "ST052-M"
    assert screen.PARENT_HEAD == "b3b8bfdbe1077f9ec967d158602951997d81e17d"
    assert screen.UPSTREAM_PARENT_ID == "ST051-B"
    assert screen.SOURCE_ALPHA == 2.520520814687742
    assert screen.GAIN == 0.05
    assert screen.INNER_WINDOW == (0.30, 1.05)
    assert screen.OUTER_WINDOW == (0.95, 1.85)
    assert screen.TIMES == (0.25, 0.50, 0.75)
    assert screen.RADII == (0.6, 0.9, 1.2)
    assert screen.CRITERIA["inner_gain_floor"] == 0.04
    assert screen.CRITERIA["mid_gain_floor"] == 0.025
    assert screen.CRITERIA["outer_gain_ceiling"] == 0.0
    assert screen.CRITERIA["response_condition_max"] == 5.0
    assert screen.CRITERIA["response_abs_cosine_max"] == 0.90


def test_profile_has_frozen_inner_mid_outer_signs_and_support():
    vals = screen.h_profile(np.array([0.0, 0.6, 0.9, 1.2, 2.0]))
    assert vals[0] == 0.0
    assert vals[1] > 0.0
    assert vals[2] > 0.0
    assert vals[3] < 0.0
    assert vals[4] == 0.0


def test_clean_rule_is_fail_closed_in_capacity_and_rank():
    row = {
        "angular_gain_r06": 0.041,
        "angular_gain_r09": 0.026,
        "angular_gain_r12": -0.01,
        "normalization": 1.001,
        "axial_rms_relative_changes": [0.0, 0.001, -0.001],
        "radial_rms_relative_changes": [-0.01, 0.0, 0.01],
        "collar_ratios": [1.0, 1.0, 1.01],
        "core_signs_pass": True,
        "support_max_abs": 0.0,
        "divergence_fd_max": 1e-8,
        "response_rank": 2,
        "response_condition_number": 2.0,
        "response_cosine": 0.2,
    }
    assert screen.clean_rule(row)
    bad = dict(row)
    bad["angular_gain_r06"] = 0.039
    assert not screen.clean_rule(bad)
    bad = dict(row)
    bad["response_rank"] = 1
    assert not screen.clean_rule(bad)
    bad = dict(row)
    bad["response_cosine"] = 0.91
    assert not screen.clean_rule(bad)
