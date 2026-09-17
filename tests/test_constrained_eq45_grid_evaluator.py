from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate


ROOT = Path(__file__).resolve().parents[1]
SEED_ARTIFACT = ROOT / "artifacts/constrained/eq45_velocity_candidate_seed.json"


def test_reloaded_eq45_candidate_grid_matches_direct_public_velocity(tmp_path):
    candidate = Eq45VelocityCandidate.load_json(SEED_ARTIFACT)

    x = np.array([-0.75, 0.0, 0.55])
    y = np.array([-0.45, 0.25])
    z = np.array([-0.8, 0.0, 0.65])
    times = np.array([0.25, 0.50, 0.75])

    grid = candidate.grid(x, y, z, times)
    assert grid.shape == (times.size, x.size, y.size, z.size, 3)
    assert np.all(np.isfinite(grid))
    assert np.max(np.linalg.norm(grid, axis=-1)) > 1e-6

    tt, xx, yy, zz = np.meshgrid(times, x, y, z, indexing="ij")
    points = np.stack((xx, yy, zz), axis=-1)
    direct = candidate.at_points(points, tt)
    assert np.array_equal(grid, direct)

    saved = tmp_path / "eq45_candidate.json"
    candidate.save_json(saved)
    restored = Eq45VelocityCandidate.load_json(saved)
    assert restored.sha256 == candidate.sha256
    assert np.array_equal(restored.grid(x, y, z, times), grid)

    truth = restored.to_dict()["truth_boundary"]
    assert truth["velocity_export_ready"] is True
    assert truth["visualization_ready"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False


def test_eq45_grid_fail_closed_on_malformed_axes_and_time():
    candidate = Eq45VelocityCandidate.seed()
    axis = np.array([-0.5, 0.0, 0.5])
    times = np.array([0.25, 0.50, 0.75])

    with pytest.raises(ValueError, match="nonempty 1D"):
        candidate.grid(np.array([]), axis, axis, times)
    with pytest.raises(ValueError, match="nonempty 1D"):
        candidate.grid(axis.reshape(1, -1), axis, axis, times)
    with pytest.raises(ValueError, match="finite"):
        candidate.grid(axis, np.array([-0.5, np.nan, 0.5]), axis, times)
    with pytest.raises(ValueError, match="declared interval"):
        candidate.grid(axis, axis, axis, np.array([0.25, 0.751]))
