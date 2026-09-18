import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_spacetime_damped_quadratic_gain import (
    solve_shared_spacetime_damping,
)


def _cell(time, z, target_value, *, l1=0.1):
    target = np.tile(np.array([[target_value, 0.0]], dtype=float), (4, 1))
    linear = np.tile(np.array([[1.0, 0.0]], dtype=float), (4, 1))
    quadratic = np.tile(np.array([[1.0, 0.0]], dtype=float), (4, 1))
    return {
        "time": time,
        "z": z,
        "target_stress": target,
        "active_target_mask": np.ones(4, dtype=bool),
        "measurement": {
            "linear_covariance_change_theta_axial": linear,
            "self_covariance_theta_axial": quadratic,
            "exact_covariance_change_theta_axial": linear + quadratic,
            "quadratic_identity_passed": True,
            "aggregate_l1_update": l1,
        },
    }


def test_shared_exact_quartic_recovers_common_half_step():
    rows = [_cell(0.375, 0.06, 0.75), _cell(0.625, 0.10, 0.75)]
    report = solve_shared_spacetime_damping(rows)

    assert report["cell_count"] == 2
    assert report["shared_selected_damping"] == pytest.approx(0.5, abs=1e-12)
    assert report["shared_best_relative_stress_residual_rms"] < 1e-12
    assert report["shared_damped_aggregate_l1_update"] == pytest.approx(0.05)
    assert report["shared_damping_improves_every_cell"] is True
    for cell in report["cells"]:
        assert cell["individually_optimal_damping"] == pytest.approx(0.5, abs=1e-12)
        assert cell["shared_damping_relative_stress_residual_rms"] < 1e-12


def test_shared_damping_is_one_compromise_not_per_cell_retuning():
    rows = [_cell(0.375, 0.06, 0.75), _cell(0.625, 0.10, 0.24)]
    report = solve_shared_spacetime_damping(rows)

    local = [cell["individually_optimal_damping"] for cell in report["cells"]]
    shared = report["shared_selected_damping"]
    assert abs(local[0] - local[1]) > 0.1
    assert 0.0 < shared < 1.0
    assert abs(shared - local[0]) > 1e-4
    assert abs(shared - local[1]) > 1e-4
    assert report["objective_weighting"].startswith("equal weight per active radial node")


def test_shared_damping_rejects_cell_specific_step_budget():
    rows = [_cell(0.375, 0.06, 0.75, l1=0.1), _cell(0.625, 0.10, 0.75, l1=0.09)]
    with pytest.raises(ValueError, match="same proposed bounded update"):
        solve_shared_spacetime_damping(rows)


def test_shared_damping_rejects_broken_quadratic_identity():
    rows = [_cell(0.375, 0.06, 0.75), _cell(0.625, 0.10, 0.75)]
    rows[1]["measurement"]["exact_covariance_change_theta_axial"] = np.zeros((4, 2))
    with pytest.raises(ValueError, match="inconsistent with linear\+self"):
        solve_shared_spacetime_damping(rows)
