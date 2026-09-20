from __future__ import annotations

import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction import load_st054
from openai_ns_reconstruction.st054_vtk import (
    _vtk_scalar_rows,
    _vtk_vector_rows,
    export_st054_vtk,
    sample_st054_grid,
    verify_vtk_export,
    verify_vtk_manifest,
)


def _read_vectors_and_speed(path, point_count):
    lines = path.read_text(encoding="ascii").splitlines()
    vector_index = lines.index("VECTORS velocity double")
    vectors = np.array(
        [
            [float(x) for x in row.split()]
            for row in lines[vector_index + 1 : vector_index + 1 + point_count]
        ],
        dtype=float,
    )
    scalar_index = lines.index("SCALARS speed double 1")
    assert lines[scalar_index + 1] == "LOOKUP_TABLE default"
    speed = np.array(
        [float(x) for x in lines[scalar_index + 2 : scalar_index + 2 + point_count]],
        dtype=float,
    )
    return vectors, speed


def test_vtk_point_order_is_x_fastest():
    velocity = np.zeros((2, 3, 2, 3), dtype=float)
    scalar = np.zeros((2, 3, 2), dtype=float)
    for i in range(2):
        for j in range(3):
            for k in range(2):
                marker = 100.0 * k + 10.0 * j + i
                velocity[i, j, k] = (marker, marker + 0.1, marker + 0.2)
                scalar[i, j, k] = marker
    expected = [0, 1, 10, 11, 20, 21, 100, 101, 110, 111, 120, 121]
    assert _vtk_vector_rows(velocity)[:, 0].tolist() == expected
    assert _vtk_scalar_rows(scalar).tolist() == expected


@pytest.mark.parametrize("model_id", ["ST054-Q2", "ST054-M3"])
def test_small_export_round_trips_exact_text_values(tmp_path, model_id):
    vtk_path = tmp_path / f"{model_id}.vtk"
    manifest = export_st054_vtk(vtk_path, model_id, 0.503, grid_size=5)
    manifest_path = vtk_path.with_suffix(".vtk.json")
    loaded = verify_vtk_export(vtk_path, manifest_path)
    assert loaded == manifest
    assert manifest["sampling_only"] is True
    assert manifest["grid_interpolation_used"] is False
    assert manifest["truth_boundary"]["vtk_export_ready"] is True
    assert manifest["truth_boundary"]["visualization_ready"] is False
    assert manifest["truth_boundary"]["pde_validated"] is False
    assert manifest["grid"]["dimensions"] == [5, 5, 5]

    x, y, z, velocity, speed = sample_st054_grid(model_id, 0.503, 5)
    assert np.array_equal(x, y)
    assert np.array_equal(y, z)
    vectors_text, speed_text = _read_vectors_and_speed(vtk_path, 5**3)
    assert np.array_equal(vectors_text, _vtk_vector_rows(velocity))
    assert np.array_equal(speed_text, _vtk_scalar_rows(speed))

    field = load_st054(model_id)
    ix, iy, iz = 3, 1, 4
    point_velocity = field.velocity(x[ix], y[iy], z[iz], 0.503)
    vtk_offset = iz * 25 + iy * 5 + ix
    assert np.max(np.abs(vectors_text[vtk_offset] - point_velocity)) <= 1.0e-12


def test_export_guards_and_truth_boundary_fail_closed(tmp_path):
    with pytest.raises(ValueError):
        export_st054_vtk(tmp_path / "bad.vtk", "not-a-model", 0.5, grid_size=5)
    with pytest.raises(ValueError):
        export_st054_vtk(tmp_path / "bad.vtk", "ST054-Q2", 0.2, grid_size=5)
    with pytest.raises(ValueError):
        export_st054_vtk(tmp_path / "bad.vtk", "ST054-Q2", 0.5, grid_size=2)
    with pytest.raises(ValueError):
        export_st054_vtk(tmp_path / "bad.txt", "ST054-Q2", 0.5, grid_size=5)

    vtk_path = tmp_path / "once.vtk"
    manifest = export_st054_vtk(vtk_path, "ST054-Q2", 0.5, grid_size=3)
    with pytest.raises(FileExistsError):
        export_st054_vtk(vtk_path, "ST054-Q2", 0.5, grid_size=3)

    promoted = copy.deepcopy(manifest)
    promoted["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="pde_validated"):
        verify_vtk_manifest(promoted)

    manifest_path = vtk_path.with_suffix(".vtk.json")
    on_disk = json.loads(manifest_path.read_text(encoding="utf-8"))
    on_disk["truth_boundary"]["visualization_ready"] = True
    manifest_path.write_text(json.dumps(on_disk), encoding="utf-8")
    with pytest.raises(ValueError, match="visualization_ready"):
        verify_vtk_export(vtk_path, manifest_path)


def test_vtk_byte_tamper_rejected(tmp_path):
    vtk_path = tmp_path / "tamper.vtk"
    export_st054_vtk(vtk_path, "ST054-M3", 0.5, grid_size=3)
    manifest_path = vtk_path.with_suffix(".vtk.json")
    vtk_path.write_text(vtk_path.read_text(encoding="ascii") + "\n", encoding="ascii")
    with pytest.raises(ValueError, match="byte identity"):
        verify_vtk_export(vtk_path, manifest_path)
