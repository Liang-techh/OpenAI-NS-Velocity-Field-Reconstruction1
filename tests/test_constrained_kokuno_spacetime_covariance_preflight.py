import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_missing_covariance_column_target import (
    PROFILE_Z,
    build_missing_covariance_column_target,
)
from openai_ns_reconstruction.kokuno_spacetime_covariance_preflight import (
    REFERENCE_BUDGET_Z,
    SPATIAL_SCREEN_Z,
    evaluate_spacetime_velocity_column_preflight,
)


def _basis(points):
    points = np.asarray(points, dtype=float)
    x, y, _ = points.T
    r = np.hypot(x, y)
    er = np.column_stack((x / r, y / r, np.zeros_like(r)))
    et = np.column_stack((-y / r, x / r, np.zeros_like(r)))
    ez = np.column_stack((np.zeros_like(r), np.zeros_like(r), np.ones_like(r)))
    return r, er, et, ez


def _theta_unit_column(points, time, phase):
    r, er, et, _ = _basis(points)
    c = np.cos(phase)
    return (r * c)[:, None] * er + (2.0 * r * c)[:, None] * et


def _axial_unit_column(points, time, phase):
    r, er, _, ez = _basis(points)
    c = np.cos(phase)
    return (r * c)[:, None] * er + (2.0 * r * c)[:, None] * ez


def _synthetic_cells(amplitude=0.25):
    radii = np.linspace(0.1, 0.4, 9)
    current = np.column_stack((2.0 * amplitude * radii**2, np.zeros_like(radii)))
    second = np.column_stack((np.zeros_like(radii), 2.0 * amplitude * radii**2))
    target = 0.01 * current + 0.02 * second
    cells = []
    for z_scale, z in ((1.0, 0.06), (1.2, 0.10)):
        for time_scale, time in ((1.0, 0.4), (1.1, 0.6)):
            receipt = build_missing_covariance_column_target(
                radii,
                z_scale * time_scale * target,
                current,
            )
            cells.append((time, z, receipt))
    return cells


def test_spacetime_independent_velocity_column_passes_every_cell():
    audit = evaluate_spacetime_velocity_column_preflight(
        _synthetic_cells(),
        _axial_unit_column,
        amplitude=0.25,
        angular_count=8,
        phase_count=16,
        coefficient_budget=0.05,
    )
    assert audit["cell_count"] == 4
    assert audit["all_cells_bounded_inverse_preflight_passed"] is True
    assert audit["finite_cycle_rerun_allowed"] is True
    assert audit["rank2_required_nodes_total"] == audit["required_nodes_total"]
    assert all(
        item["finite_cycle_rerun_allowed_for_cell"] is True
        for item in audit["cells"]
    )


def test_spacetime_duplicate_velocity_column_is_rejected_everywhere():
    audit = evaluate_spacetime_velocity_column_preflight(
        _synthetic_cells(),
        _theta_unit_column,
        amplitude=0.25,
        angular_count=8,
        phase_count=16,
        coefficient_budget=0.05,
    )
    assert audit["finite_cycle_rerun_allowed"] is False
    assert audit["required_nodes_total"] > 0
    assert audit["rank2_required_nodes_total"] == 0
    assert all(
        item["finite_cycle_rerun_allowed_for_cell"] is False
        for item in audit["cells"]
    )


def test_spacetime_preflight_rejects_duplicate_pairs_and_mismatched_grids():
    cells = _synthetic_cells()
    with pytest.raises(ValueError, match="unique"):
        evaluate_spacetime_velocity_column_preflight(
            [cells[0], cells[0]],
            _axial_unit_column,
            amplitude=0.25,
            coefficient_budget=0.05,
        )

    bad_radii = np.linspace(0.11, 0.41, 9)
    current = np.column_stack((bad_radii**2, np.zeros_like(bad_radii)))
    target = current + np.column_stack((np.zeros_like(bad_radii), bad_radii**2))
    bad = build_missing_covariance_column_target(bad_radii, target, current)
    with pytest.raises(ValueError, match="same radial grid"):
        evaluate_spacetime_velocity_column_preflight(
            [cells[0], (0.7, 0.12, bad)],
            _axial_unit_column,
            amplitude=0.25,
            coefficient_budget=0.05,
        )


def test_default_axial_screen_is_fixed_core_safe_and_contains_reference():
    assert SPATIAL_SCREEN_Z == (0.06, float(PROFILE_Z), 0.10)
    assert REFERENCE_BUDGET_Z == float(PROFILE_Z)
    assert list(SPATIAL_SCREEN_Z) == sorted(SPATIAL_SCREEN_Z)
