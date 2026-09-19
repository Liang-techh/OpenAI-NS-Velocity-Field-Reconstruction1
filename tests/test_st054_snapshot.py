from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat

from openai_ns_reconstruction.st054_snapshot import (
    ST054_REFERENCE_ATOL,
    ST054_SOURCE_COMMIT,
    ST054_TESTED_MAT_SHA256,
    available_st054_models,
    load_st054,
)


ROOT = Path(__file__).resolve().parents[1]
MAT = ROOT / "visualization" / "matlab" / "data" / "st054_models.mat"


def _models():
    payload = loadmat(MAT, simplify_cells=True)
    models = payload["models"]
    if isinstance(models, dict):
        return [models]
    return np.asarray(models, dtype=object).reshape(-1).tolist()


def _model(model_id: str):
    return next(model for model in _models() if str(model["id"]) == model_id)


@pytest.mark.parametrize("model_id", available_st054_models())
def test_velocity_replays_frozen_source_reference(model_id: str):
    field = load_st054(model_id)
    source = _model(model_id)
    points = np.asarray(source["ref_points"], dtype=float)
    times = np.asarray(source["ref_times"], dtype=float).reshape(-1)
    expected = np.asarray(source["ref_velocity"], dtype=float)
    assert expected.shape == (times.size, points.shape[0], 3)

    errors = []
    for index, time in enumerate(times):
        actual = field.velocity(points[:, 0], points[:, 1], points[:, 2], float(time))
        errors.append(float(np.max(np.abs(actual - expected[index]))))
    assert max(errors) <= ST054_REFERENCE_ATOL


@pytest.mark.parametrize("model_id", available_st054_models())
def test_scalar_and_broadcast_api(model_id: str):
    field = load_st054(model_id)
    scalar = field.velocity(0.1, 0.0, 0.1, 0.5)
    assert scalar.shape == (3,)
    assert np.isfinite(scalar).all()

    x = np.array([[0.0], [0.1]])
    y = np.array([0.0, 0.2, -0.2])
    z = 0.1
    value = field.velocity(x, y, z, 0.5)
    assert value.shape == (2, 3, 3)
    assert np.isfinite(value).all()


@pytest.mark.parametrize("model_id", available_st054_models())
def test_loader_binds_native_matlab_tested_payload_and_truth_boundary(model_id: str):
    field = load_st054(model_id)
    assert field.source_commit == ST054_SOURCE_COMMIT
    assert field.mat_sha256 == ST054_TESTED_MAT_SHA256
    assert field.nu == 0.01
    assert (field.tmin, field.tmax) == (0.25, 0.75)
    assert (field.support_radius, field.support_z) == (2.0, 2.0)


@pytest.mark.parametrize("bad_time", [0.249999, 0.750001, np.nan, np.inf])
def test_velocity_rejects_invalid_time(bad_time: float):
    with pytest.raises(ValueError):
        load_st054("ST054-Q2").velocity(0.0, 0.0, 0.0, bad_time)


def test_loader_rejects_unknown_model():
    with pytest.raises(ValueError, match="Unknown ST054 model"):
        load_st054("ST054-HIDDEN")


def test_loader_rejects_mutated_mat_bytes(tmp_path: Path):
    mutated = tmp_path / "st054_models.mat"
    data = bytearray(MAT.read_bytes())
    data[-1] ^= 1
    mutated.write_bytes(data)
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_st054("ST054-Q2", mutated)
