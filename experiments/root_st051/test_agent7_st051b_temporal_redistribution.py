from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import agent7_st051b_temporal_redistribution as screen


def test_temporal_redistribution_contract_is_frozen():
    assert screen.PREREG_ISSUE == 488
    assert screen.PARENT_ID == "ST051-B"
    assert screen.PARENT_HEAD == "4b784f1b8457af2ead49295631d834d4e882000b"
    assert screen.STACK_BASE_PR == 480
    assert screen.STACK_BASE_HEAD == "187157d377b14bae56415648ac100342f2b8bbb0"
    assert screen.SOURCE_PROFILE_PR == 458
    assert screen.TRANSFER_PR == 469
    assert screen.SOURCE_ALPHA == 2.520520814687742
    assert screen.INNER_WINDOW == (0.30, 1.05)
    assert screen.OUTER_WINDOW == (0.95, 1.85)
    assert screen.BASE_GAIN == 0.025
    assert screen.GAMMAS == (0.0, 0.010, 0.015, 0.020, 0.025)
    assert screen.GRID_SIZE == 33
    assert screen.SENSITIVITY_EPS == 1e-4


def test_temporal_shape_keeps_reference_frame_frozen():
    for gamma in screen.GAMMAS:
        assert screen.time_shape(0.25) == 0.0
        assert screen.gain_at_time(gamma, 0.25) == screen.BASE_GAIN
    assert screen.time_shape(0.50) == 0.5
    assert screen.time_shape(0.75) == 1.0
    assert screen.gain_at_time(0.025, 0.75) == 0.05


def test_clean_rule_requires_independent_well_conditioned_temporal_capacity():
    angular = {
        str(0.25): {
            str(0.6): {"additional_gain_from_static": 0.0},
            str(0.9): {"additional_gain_from_static": 0.0},
            str(1.2): {"additional_gain_from_static": 0.0},
        },
        str(0.5): {
            str(0.6): {"additional_gain_from_static": 0.01},
            str(0.9): {"additional_gain_from_static": 0.006},
            str(1.2): {"additional_gain_from_static": -0.02},
        },
        str(0.75): {
            str(0.6): {"additional_gain_from_static": 0.021},
            str(0.9): {"additional_gain_from_static": 0.011},
            str(1.2): {"additional_gain_from_static": -0.04},
        },
    }
    row = dict(
        gamma=0.02,
        angular=angular,
        early_field_max_abs_delta=0.0,
        axial_rms_relative_changes_vs_static=[0.0, 0.001, 0.002],
        radial_rms_relative_changes_vs_static=[0.0, 0.001, 0.002],
        collar_ratios_vs_static=[1.0, 1.01, 1.02],
        energies=[1.0, 0.8, 0.7],
        reference_energy_relative_error=0.0,
        core_signs_pass=True,
        support_max_abs=0.0,
        divergence_fd_max=1e-8,
    )
    sensitivity = dict(
        normalized_rank=2,
        normalized_condition_number=3.0,
        normalized_column_cosine=0.4,
    )
    assert screen.clean_temporal_rule(row, sensitivity)
    sensitivity["normalized_rank"] = 1
    assert not screen.clean_temporal_rule(row, sensitivity)
    sensitivity["normalized_rank"] = 2
    row["gamma"] = 0.0
    assert not screen.clean_temporal_rule(row, sensitivity)


def test_no_new_spatial_basis_is_introduced_by_design():
    assert screen.SOURCE_ALPHA == screen.base.SOURCE_ALPHA
    assert screen.INNER_WINDOW == screen.base.INNER_WINDOW
    assert screen.OUTER_WINDOW == screen.base.OUTER_WINDOW
    assert max(screen.GAMMAS) == 0.025
