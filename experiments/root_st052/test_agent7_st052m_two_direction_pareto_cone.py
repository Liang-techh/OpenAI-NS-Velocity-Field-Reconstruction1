from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import agent7_st052m_two_direction_pareto_cone as audit


def test_swirl_column_is_centered_in_normalized_sign_coordinate():
    report = {
        "oriented_incremental_desirabilities": {
            "plus": {
                "aspect_gain_increment": 0.4,
                "tip_thinning_increment": -0.2,
                "turns_fidelity_increment": 0.6,
                "axial_pair_fidelity_increment": 0.0,
            },
            "minus": {
                "aspect_gain_increment": -0.2,
                "tip_thinning_increment": 0.4,
                "turns_fidelity_increment": -0.4,
                "axial_pair_fidelity_increment": 0.0,
            },
        }
    }
    got = audit._swirl_column(report)
    np.testing.assert_allclose(got, [0.3, -0.3, 0.5, 0.0])


def test_shoulder_column_converts_per_lambda_to_normalized_coordinate():
    report = {
        "oriented_desirability_sensitivities": {
            "aspect_desirability": 2.0,
            "tip_thinning_desirability": -3.0,
            "path_turns_desirability": 4.0,
            "path_pair_axial_desirability": 5.0,
        }
    }
    got = audit._shoulder_column(report)
    np.testing.assert_allclose(
        got,
        audit.shoulder.EPSILON * np.asarray([2.0, -3.0, 4.0, 5.0]),
    )


def test_pareto_solver_finds_common_improving_direction_when_present():
    j = np.asarray(
        [
            [0.20, 0.00],
            [0.00, 0.20],
            [0.10, 0.10],
            [0.00, 0.00],
        ],
        dtype=float,
    )
    result = audit.solve_pareto_cone(j)
    assert result["target_free_local_pareto_direction"] is not None


def test_pareto_solver_rejects_conflicting_cone():
    # Any gain in metric 0 from x0 worsens metric 1 by the same amount, while
    # x1 is unable to change those rows.  With improvement required above the
    # same tolerance used for non-worsening, no Pareto witness exists.
    j = np.asarray(
        [
            [0.20, 0.00],
            [-0.20, 0.00],
            [0.00, 0.10],
            [0.00, -0.10],
        ],
        dtype=float,
    )
    result = audit.solve_pareto_cone(j)
    assert result["target_free_local_pareto_direction"] is None


def test_truth_boundary_is_fail_closed():
    assert audit.TRUTH["canonical_velocity_changed"] is False
    assert audit.TRUTH["combined_velocity_child_constructed"] is False
    assert audit.TRUTH["held_out_pde_residual_evaluated"] is False
    assert audit.TRUTH["visual_correspondence_verified"] is False
    assert audit.TRUTH["pde_validated"] is False
