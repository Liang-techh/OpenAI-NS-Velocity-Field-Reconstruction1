from pathlib import Path

import numpy as np
import pytest
from scipy.io import netcdf_file

from openai_ns_reconstruction import constrained_st048s_visualization_grid as mod


def _affine_velocity(points: np.ndarray, time: float) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    return np.column_stack(
        (
            points[:, 0] + time,
            points[:, 1] - 2.0 * time,
            points[:, 2] + 3.0 * time,
        )
    )


def test_cross_language_grid_roundtrip_and_component_order(tmp_path: Path):
    path = tmp_path / "grid.nc"
    result = mod._export_velocity_netcdf(_affine_velocity, path, resolution=5)
    loaded = mod.load_st048s_netcdf(path)

    assert result["grid_sha256"] == loaded["grid_sha256"]
    assert loaded["u"].shape == (3, 5, 5, 5)
    center = 2
    assert loaded["x"][center] == 0.0
    assert loaded["y"][center] == 0.0
    assert loaded["z"][center] == 0.0
    assert loaded["u"][0, center, center, center] == pytest.approx(0.25)
    assert loaded["v"][0, center, center, center] == pytest.approx(-0.50)
    assert loaded["w"][0, center, center, center] == pytest.approx(0.75)
    assert loaded["speed"][0, center, center, center] == pytest.approx(
        np.sqrt(0.25**2 + 0.50**2 + 0.75**2)
    )


def test_truth_state_laundering_is_rejected(tmp_path: Path):
    path = tmp_path / "grid.nc"
    mod._export_velocity_netcdf(_affine_velocity, path, resolution=5)
    with netcdf_file(str(path), mode="a") as nc:
        nc.pde_validated = 1
    with pytest.raises(ValueError, match="truth-state promotion"):
        mod.load_st048s_netcdf(path)


def test_velocity_tamper_is_rejected(tmp_path: Path):
    path = tmp_path / "grid.nc"
    mod._export_velocity_netcdf(_affine_velocity, path, resolution=5)
    with netcdf_file(str(path), mode="a") as nc:
        nc.variables["u"].data[0, 2, 2, 2] += 0.125
        u = np.asarray(nc.variables["u"].data)
        v = np.asarray(nc.variables["v"].data)
        w = np.asarray(nc.variables["w"].data)
        nc.variables["speed"].data[:] = np.sqrt(u * u + v * v + w * w)
    with pytest.raises(ValueError, match="checksum mismatch"):
        mod.load_st048s_netcdf(path)


def test_export_rejects_zero_field_bad_resolution_and_overwrite(tmp_path: Path):
    path = tmp_path / "grid.nc"

    def zero(points: np.ndarray, _time: float) -> np.ndarray:
        return np.zeros_like(points)

    with pytest.raises(ValueError, match="exact-zero"):
        mod._export_velocity_netcdf(zero, path, resolution=5)
    with pytest.raises(ValueError, match="odd"):
        mod._export_velocity_netcdf(_affine_velocity, path, resolution=6)

    mod._export_velocity_netcdf(_affine_velocity, path, resolution=5)
    with pytest.raises(FileExistsError):
        mod._export_velocity_netcdf(_affine_velocity, path, resolution=5)


def test_registered_identity_constants_are_fail_closed():
    assert mod.TASK_ID == "CR-A9-045"
    assert mod.BASE_MAIN_SHA == "f0193d66c9d92948b4820ebcb70263673995b324"
    assert mod.SOURCE_PR == 390
    assert mod.SOURCE_HEAD_SHA == "97695a86f85ce68fb4ae70c41fc81c904d655183"
    assert mod.CANDIDATE_ID == "ST048-S"
    assert mod.CANDIDATE_RAW_SHA256 == "6d9ce8407e29aca20d32599332ae3ec654f45678e17783179267665428ec8e09"
    assert np.array_equal(mod.REFERENCE_TIMES, np.asarray([0.25, 0.50, 0.75]))
