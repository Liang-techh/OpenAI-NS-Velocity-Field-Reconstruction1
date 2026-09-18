from pathlib import Path

import numpy as np
import pytest
from scipy.io import netcdf_file

from research_baseline import CANDIDATE_SHA256, load_best
from research_baseline.netcdf import export_st006_netcdf, load_st006_netcdf


def test_st006_netcdf_round_trip_matches_callable(tmp_path: Path):
    path = tmp_path / "st006.nc"
    receipt = export_st006_netcdf(path, resolution=5)
    grid = load_st006_netcdf(path)

    assert receipt["candidate_sha256"] == CANDIDATE_SHA256
    assert receipt["grid_sha256"] == grid.grid_sha256
    assert receipt["dimension_order"] == "time,x,y,z"
    assert receipt["component_order"] == "u,v,w"
    assert receipt["velocity_export_ready"] is True
    assert grid.u.shape == (3, 5, 5, 5)
    assert np.array_equal(grid.time, np.array([0.25, 0.50, 0.75]))
    assert grid.x[0] == grid.y[0] == grid.z[0] == -2.0
    assert grid.x[-1] == grid.y[-1] == grid.z[-1] == 2.0

    field = load_best()
    index = 3
    point = np.array([grid.x[index], grid.y[2], grid.z[1]])
    expected = np.asarray(field.at_points(point, 0.50), dtype=float)
    actual = np.array([grid.u[1, index, 2, 1], grid.v[1, index, 2, 1], grid.w[1, index, 2, 1]])
    assert np.allclose(actual, expected, rtol=0.0, atol=2.0e-13)
    assert np.allclose(grid.speed, np.sqrt(grid.u**2 + grid.v**2 + grid.w**2), rtol=0.0, atol=2.0e-14)

    for array in (grid.time, grid.x, grid.y, grid.z, grid.u, grid.v, grid.w, grid.speed):
        assert array.flags.writeable is False
    assert grid.visualization_ready is False
    assert grid.visual_correspondence_verified is False
    assert grid.pde_validated is False
    assert grid.openai_field_identified is False


def test_st006_netcdf_has_named_cross_language_variables(tmp_path: Path):
    path = tmp_path / "st006.nc"
    export_st006_netcdf(path, resolution=5)

    with netcdf_file(str(path), mode="r", mmap=False) as nc:
        assert set(nc.variables) == {"time", "x", "y", "z", "u", "v", "w", "speed"}
        assert nc.variables["u"].dimensions == ("time", "x", "y", "z")
        assert nc.variables["v"].dimensions == ("time", "x", "y", "z")
        assert nc.variables["w"].dimensions == ("time", "x", "y", "z")
        assert bytes(nc.component_order).decode("ascii") == "u,v,w"
        assert bytes(nc.delivery_semantics).decode("ascii") == "visualization_only_not_pde_evidence"


def test_st006_netcdf_rejects_checksum_tamper_even_with_consistent_speed(tmp_path: Path):
    path = tmp_path / "st006.nc"
    export_st006_netcdf(path, resolution=5)

    with netcdf_file(str(path), mode="a", mmap=False) as nc:
        i = (1, 2, 2, 2)
        u = float(nc.variables["u"][i]) + 0.125
        v = float(nc.variables["v"][i])
        w = float(nc.variables["w"][i])
        nc.variables["u"][i] = u
        nc.variables["speed"][i] = np.sqrt(u * u + v * v + w * w)

    with pytest.raises(ValueError, match="checksum"):
        load_st006_netcdf(path)


def test_st006_netcdf_rejects_truth_state_laundering(tmp_path: Path):
    path = tmp_path / "st006.nc"
    export_st006_netcdf(path, resolution=5)

    with netcdf_file(str(path), mode="a", mmap=False) as nc:
        nc.pde_validated = 1

    with pytest.raises(ValueError, match="pde_validated"):
        load_st006_netcdf(path)


def test_st006_netcdf_fail_closed_inputs_and_overwrite(tmp_path: Path):
    with pytest.raises(ValueError, match="odd"):
        export_st006_netcdf(tmp_path / "even.nc", resolution=6)
    with pytest.raises(ValueError, match=r"\.nc"):
        export_st006_netcdf(tmp_path / "wrong.bin", resolution=5)

    path = tmp_path / "st006.nc"
    export_st006_netcdf(path, resolution=5)
    with pytest.raises(FileExistsError):
        export_st006_netcdf(path, resolution=5)

    second = export_st006_netcdf(path, resolution=5, overwrite=True)
    assert second["candidate_sha256"] == CANDIDATE_SHA256
    assert load_st006_netcdf(path).grid_sha256 == second["grid_sha256"]
