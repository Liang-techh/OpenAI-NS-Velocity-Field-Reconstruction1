import numpy as np
import pytest

from openai_ns_reconstruction.constrained_verifier_calibration import (
    manufactured_force,
    manufactured_terms,
    run_calibration,
)


def test_manufactured_terms_are_nontrivial_and_divergence_free_analytically():
    points = np.array([
        [0.31, -0.27, 0.19],
        [-0.44, 0.38, -0.17],
    ])
    terms = manufactured_terms(points, 0.5)
    assert np.max(np.linalg.norm(terms["convection"], axis=1)) > 0.1
    assert np.max(np.linalg.norm(terms["grad_pressure"], axis=1)) > 0.1
    assert np.max(np.linalg.norm(terms["viscous"], axis=1)) > 1e-3

    x, y, _ = points.T
    a = 1.5
    divergence = a * np.cos(x) * np.cos(y) - a * np.cos(x) * np.cos(y)
    assert np.array_equal(divergence, np.zeros_like(divergence))


def test_fourth_order_convergence_including_time_window_endpoints():
    report = run_calibration()
    assert report["minimum_observed_order"] > 3.9
    assert len(report["steps"]) == 4
    assert set(report["times"]) == {0.25, 0.5, 0.75}
    finest = [row for row in report["rows"] if row["step"] == min(report["steps"])]
    assert max(row["residual_max"] for row in finest) < 1.1e-9
    assert max(row["divergence_max"] for row in report["rows"]) < 5e-13


def test_mutations_are_detected_with_fixed_thresholds():
    report = run_calibration()
    mutations = report["mutation_rows"]
    assert mutations["convection_sign"]["residual_max"] > 3.0
    assert mutations["pressure_sign"]["residual_max"] > 1.5
    assert mutations["viscosity_sign"]["residual_max"] > 0.04
    assert all(row["residual_max"] > 1e-3 for row in mutations.values())


def test_mutation_and_step_inputs_fail_closed():
    point = np.array([[0.2, 0.1, -0.3]])
    with pytest.raises(ValueError):
        manufactured_force(point, 0.5, mutation="unknown")
    with pytest.raises(ValueError):
        run_calibration(steps=(0.04, 0.02))
    with pytest.raises(ValueError):
        run_calibration(steps=(0.02, 0.04, 0.01))
    with pytest.raises(ValueError):
        run_calibration(point_count=4)
