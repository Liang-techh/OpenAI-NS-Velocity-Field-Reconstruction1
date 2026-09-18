import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_spacetime_damping_holdout import (
    evaluate_fixed_damping_on_disjoint_cells,
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


def test_fixed_damping_transfers_without_validation_retuning():
    rows = [_cell(0.3125, 0.07, 0.75), _cell(0.6875, 0.09, 0.75)]
    report = evaluate_fixed_damping_on_disjoint_cells(
        rows,
        damping=0.5,
        calibration_pairs={(0.375, 0.06), (0.625, 0.10)},
    )

    assert report["fixed_damping"] == pytest.approx(0.5)
    assert report["aggregate_fixed_damping_relative_stress_residual_rms"] < 1e-12
    assert report["fixed_damped_aggregate_l1_update"] == pytest.approx(0.05)
    assert report["validation_transfer_preflight_passed"] is True
    assert report["heldout_damping_reoptimized"] is False
    assert report["validation_objective_used_for_selection"] is False
    for cell in report["cells"]:
        assert cell["fixed_damping_relative_stress_residual_rms"] < 1e-12
        assert "individually_optimal_damping" not in cell


def test_fixed_damping_records_heldout_failure_without_retuning():
    rows = [_cell(0.3125, 0.07, 0.24), _cell(0.6875, 0.09, 0.24)]
    report = evaluate_fixed_damping_on_disjoint_cells(
        rows,
        damping=0.5,
        calibration_pairs={(0.375, 0.06), (0.625, 0.10)},
    )

    assert report["fixed_damping"] == pytest.approx(0.5)
    assert report["aggregate_fixed_damping_relative_stress_residual_rms"] > 1.0
    assert report["validation_transfer_preflight_passed"] is False
    assert report["heldout_damping_reoptimized"] is False


def test_fixed_damping_rejects_overlap_with_calibration_cells():
    rows = [_cell(0.375, 0.06, 0.75), _cell(0.6875, 0.09, 0.75)]
    with pytest.raises(ValueError, match="disjoint from damping calibration"):
        evaluate_fixed_damping_on_disjoint_cells(
            rows,
            damping=0.5,
            calibration_pairs={(0.375, 0.06), (0.625, 0.10)},
        )


def test_fixed_damping_rejects_cell_specific_update_budget():
    rows = [_cell(0.3125, 0.07, 0.75, l1=0.1), _cell(0.6875, 0.09, 0.75, l1=0.09)]
    with pytest.raises(ValueError, match="same frozen update"):
        evaluate_fixed_damping_on_disjoint_cells(
            rows,
            damping=0.5,
            calibration_pairs={(0.375, 0.06), (0.625, 0.10)},
        )


def test_fixed_damping_rejects_broken_quadratic_identity():
    rows = [_cell(0.3125, 0.07, 0.75), _cell(0.6875, 0.09, 0.75)]
    rows[1]["measurement"]["exact_covariance_change_theta_axial"] = np.zeros((4, 2))
    with pytest.raises(ValueError, match=r"inconsistent with linear\+self"):
        evaluate_fixed_damping_on_disjoint_cells(
            rows,
            damping=0.5,
            calibration_pairs={(0.375, 0.06), (0.625, 0.10)},
        )
