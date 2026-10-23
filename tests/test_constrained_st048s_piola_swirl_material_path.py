import copy

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_st048s_piola_swirl_material_path import (
    SWIRL_KAPPA,
    TRUTH_BOUNDARY,
    TemporalPiolaSwirlVelocity,
    _apply_swirl_gain_xy,
    _fixed_seed_table,
    audit_truth_boundary,
    beta_at_time,
    measure_material_paths,
)


def _cylindrical_components(points, velocity):
    points = np.asarray(points, dtype=float)
    velocity = np.asarray(velocity, dtype=float)
    radius = np.hypot(points[:, 0], points[:, 1])
    rx = points[:, 0] / radius
    ry = points[:, 1] / radius
    radial = rx * velocity[:, 0] + ry * velocity[:, 1]
    azimuthal = -ry * velocity[:, 0] + rx * velocity[:, 1]
    return radial, azimuthal, velocity[:, 2]


def test_temporal_schedule_and_swirl_gain_are_exactly_bounded():
    assert beta_at_time(0.25) == pytest.approx(0.1)
    assert beta_at_time(0.50) == pytest.approx(0.075)
    assert beta_at_time(0.75) == pytest.approx(0.1)
    with pytest.raises(ValueError):
        beta_at_time(0.2)
    with pytest.raises(ValueError):
        beta_at_time(0.5, gamma=0.2)

    points = np.array([[0.6, 0.0, 0.1], [0.0, 0.9, -0.2], [-0.8, 0.8, 0.3]])
    velocity = np.array([[0.2, 0.3, 0.4], [-0.5, 0.25, -0.1], [0.7, -0.2, 0.6]])
    before = _cylindrical_components(points, velocity)
    after_velocity = _apply_swirl_gain_xy(points, velocity, SWIRL_KAPPA)
    after = _cylindrical_components(points, after_velocity)
    assert np.allclose(after[0], before[0], rtol=0.0, atol=2e-16)
    assert np.allclose(after[1], (1.0 + SWIRL_KAPPA) * before[1], rtol=2e-15, atol=2e-16)
    assert np.array_equal(after[2], before[2])
    with pytest.raises(ValueError):
        _apply_swirl_gain_xy(points, velocity, 0.16)


def test_temporal_piola_swirl_wrapper_is_nontrivial_and_finite():
    def parent(points, time):
        points = np.asarray(points)
        return np.column_stack((
            -0.1 * points[:, 0] - 0.2 * points[:, 1],
            0.2 * points[:, 0] - 0.1 * points[:, 1],
            0.05 * points[:, 2] + 0.01 * float(time),
        ))

    field = TemporalPiolaSwirlVelocity(parent, kappa=SWIRL_KAPPA)
    points = np.array([[0.6, 0.0, -0.3], [0.0, 0.9, 0.3], [0.8, -0.4, 0.2]])
    values = field(points, 0.5)
    assert values.shape == points.shape
    assert np.isfinite(values).all()
    assert np.linalg.norm(values) > 0.0
    with pytest.raises(ValueError):
        TemporalPiolaSwirlVelocity(parent, kappa=-0.01)


def test_frozen_seed_contract_and_material_path_oracle():
    seeds, metadata = _fixed_seed_table()
    assert seeds.shape == (48, 3)
    assert len(metadata) == 48
    assert np.allclose(np.unique(np.round(np.hypot(seeds[:, 0], seeds[:, 1]), 12)), [0.6, 0.9, 1.2])
    assert np.allclose(np.unique(seeds[:, 2]), [-0.3, 0.3])

    decay = 0.2
    omega = 0.4
    stretch = 0.1

    def linear_field(points, _time):
        points = np.asarray(points, dtype=float)
        x, y, z = points.T
        return np.column_stack((
            -decay * x - omega * y,
            omega * x - decay * y,
            stretch * z,
        ))

    report = measure_material_paths(linear_field)
    dt = 0.5
    expected_turns = abs(omega * dt) / (2.0 * np.pi)
    expected_sep_ratio = np.exp(stretch * dt)
    expected_mean_radius_change = np.mean(
        np.hypot(seeds[:, 0], seeds[:, 1]) * (np.exp(-decay * dt) - 1.0)
    )
    assert report["path_count"] == 48
    assert report["paired_material_line_count"] == 24
    assert report["inward_path_count"] == 48
    assert report["outward_path_count"] == 0
    assert report["pair_axial_separation_growth_count"] == 24
    assert report["pair_axial_separation_shrink_count"] == 0
    assert report["mean_absolute_turns"] == pytest.approx(expected_turns, rel=2e-10, abs=2e-12)
    assert report["mean_pair_axial_separation_ratio"] == pytest.approx(expected_sep_ratio, rel=2e-10)
    assert report["mean_radius_change"] == pytest.approx(expected_mean_radius_change, rel=2e-10)


def test_truth_boundary_rejects_any_promotion():
    report = {
        "task_id": "CR-A9-048",
        "candidate": {"pde_validated": False},
        "agent7_formula_crosscheck": {"crosscheck_passed": True},
        "truth_boundary": copy.deepcopy(TRUTH_BOUNDARY),
    }
    audit_truth_boundary(report)
    bad = copy.deepcopy(report)
    bad["truth_boundary"]["visualization_ready"] = True
    with pytest.raises(ValueError):
        audit_truth_boundary(bad)
    bad = copy.deepcopy(report)
    bad["truth_boundary"]["held_out_pde_residual_evaluated"] = True
    with pytest.raises(ValueError):
        audit_truth_boundary(bad)
