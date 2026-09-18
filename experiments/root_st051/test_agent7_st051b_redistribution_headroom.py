from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import agent7_st051b_redistribution_headroom as screen


def test_preregistered_headroom_contract_is_frozen():
    assert screen.PREREG_ISSUE == 479
    assert screen.PARENT_ID == "ST051-B"
    assert screen.PARENT_HEAD == "4b784f1b8457af2ead49295631d834d4e882000b"
    assert screen.SOURCE_PROFILE_PR == 458
    assert screen.TRANSFER_PR == 469
    assert screen.SOURCE_ALPHA == 2.520520814687742
    assert screen.INNER_WINDOW == (0.30, 1.05)
    assert screen.OUTER_WINDOW == (0.95, 1.85)
    assert screen.GAINS == (0.0, 0.015, 0.025, 0.035, 0.05)
    assert screen.GRID_SIZES == (25, 33)
    assert screen.CRITERIA == {
        "inner_gain_floor": 0.04,
        "mid_gain_floor": 0.025,
        "outer_gain_ceiling": 0.0,
        "normalization_deviation_max": 0.005,
        "inward_speed_loss_max": 0.005,
        "axial_rms_abs_change_max": 0.02,
        "radial_rms_growth_max": 0.03,
        "collar_ratio_max": 1.25,
        "support_max_abs": 1e-12,
        "divergence_fd_max": 1e-5,
    }


def test_clean_headroom_rule_is_fail_closed_and_requires_gain_above_existing_crossing():
    row = dict(
        gain=0.05,
        angular_gain_r06=0.041,
        angular_gain_r09=0.026,
        angular_gain_r12=-0.01,
        normalization=1.002,
        inward_speed_loss=0.0,
        axial_rms_relative_changes=[0.0, -0.001, 0.001],
        radial_rms_relative_changes=[-0.01, 0.0, 0.01],
        collar_ratios=[1.0, 1.0, 1.01],
        core_signs_pass=True,
        support_max_abs=0.0,
        divergence_fd_max=1e-8,
    )
    assert screen.clean_headroom_rule(row)
    row["gain"] = 0.025
    assert not screen.clean_headroom_rule(row)
    row["gain"] = 0.05
    row["angular_gain_r09"] = 0.024
    assert not screen.clean_headroom_rule(row)


def test_no_new_basis_dimension_is_declared_by_design():
    assert len(screen.GAINS) == 5
    assert screen.GAINS[-1] == 0.05
    assert screen.SOURCE_ALPHA == screen.base.SOURCE_ALPHA
    assert screen.INNER_WINDOW == screen.base.INNER_WINDOW
    assert screen.OUTER_WINDOW == screen.base.OUTER_WINDOW
