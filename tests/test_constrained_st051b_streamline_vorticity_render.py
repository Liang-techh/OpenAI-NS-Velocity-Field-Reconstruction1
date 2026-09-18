import numpy as np
import pytest

from openai_ns_reconstruction.constrained_st051b_streamline_vorticity_render import (
    REFERENCE_TIMES,
    analyze_pair,
    sample_velocity_grid,
    seed_points,
    vorticity,
)


def _solid_rotation(points, time):
    points = np.asarray(points, dtype=float)
    return np.column_stack((-points[:, 1], points[:, 0], np.zeros(len(points))))


def test_seed_contract_is_frozen_48_paths():
    seeds = seed_points()
    assert seeds.shape == (48, 3)
    radii = np.hypot(seeds[:, 0], seeds[:, 1])
    assert sorted(set(np.round(radii, 12))) == [0.6, 0.9, 1.2]
    assert sorted(set(np.round(seeds[:, 2], 12))) == [-0.3, 0.3]


def test_vorticity_recovers_solid_rotation():
    axis, fields = sample_velocity_grid(_solid_rotation, resolution=9)
    spacing = float(axis[1] - axis[0])
    omega, magnitude = vorticity(fields[0], spacing)
    assert np.max(np.abs(omega[..., 0])) < 1e-12
    assert np.max(np.abs(omega[..., 1])) < 1e-12
    assert np.max(np.abs(omega[..., 2] - 2.0)) < 1e-12
    assert np.max(np.abs(magnitude - 2.0)) < 1e-12


def test_target_free_pair_analysis_preserves_truth_boundary():
    axis, parent = sample_velocity_grid(_solid_rotation, resolution=9)
    child = 1.05 * parent
    report = analyze_pair(axis, parent, axis, child)
    assert report["contract"]["streamline_count"] == 48
    assert report["contract"]["image_or_openai_numeric_target_used"] is False
    assert [row["time"] for row in report["per_time"]] == list(REFERENCE_TIMES)
    for row in report["per_time"]:
        parent_turns = row["parent"]["streamlines"]["mean_absolute_turns"]
        child_turns = row["child"]["streamlines"]["mean_absolute_turns"]
        assert child_turns == pytest.approx(parent_turns, rel=2e-5)
    assert all(value is False for value in report["truth"].values())
