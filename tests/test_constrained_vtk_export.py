import hashlib
import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_vtk_export import export_velocity_vtk_series


def _field(points: np.ndarray, time: float) -> np.ndarray:
    x, y, z = points.T
    return np.column_stack((x + time, 10.0 * y, -z + 2.0 * time))


def _vector_rows(text: str, count: int) -> np.ndarray:
    lines = text.splitlines()
    start = lines.index("VECTORS velocity double") + 1
    return np.array(
        [[float(value) for value in line.split()] for line in lines[start : start + count]]
    )


def test_vtk_series_preserves_xyz_point_order_and_truth_boundary(tmp_path):
    result = export_velocity_vtk_series(
        _field,
        output_dir=tmp_path / "vtk",
        x=[-1.0, 2.0],
        y=[-2.0, 0.5],
        z=[-3.0, 1.0],
        times=[0.25, 0.75],
        candidate_id="manufactured_vtk_smoke",
    )

    manifest = json.loads(result.manifest_path.read_text())
    assert manifest["claim_scope"] == "visualization_export_only"
    assert manifest["axis_order"] == "x,y,z with x-fastest VTK point ordering"
    assert manifest["claims"]["visualization_export_ready"] is True
    assert manifest["claims"]["pde_validated"] is False
    assert manifest["claims"]["paper_exact"] is False
    assert manifest["claims"]["openai_field_identified"] is False
    assert manifest["claims"]["blowup_proved"] is False

    vtk_path = result.vtk_paths[0]
    text = vtk_path.read_text(encoding="ascii")
    assert "DATASET RECTILINEAR_GRID" in text
    assert "DIMENSIONS 2 2 2" in text
    assert "POINT_DATA 8" in text
    assert "SCALARS speed double 1" in text

    X, Y, Z = np.meshgrid(
        np.array([-1.0, 2.0]),
        np.array([-2.0, 0.5]),
        np.array([-3.0, 1.0]),
        indexing="ij",
    )
    points = np.column_stack(
        (X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F"))
    )
    expected = _field(points, 0.25)
    actual = _vector_rows(text, 8)
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=0.0)

    expected_sha = hashlib.sha256(vtk_path.read_bytes()).hexdigest()
    assert manifest["files"][0]["sha256"] == expected_sha


def test_vtk_series_is_deterministic(tmp_path):
    kwargs = dict(
        x=[-1.0, 0.0, 1.0],
        y=[-1.0, 1.0],
        z=[-2.0, 2.0],
        times=[0.5],
        candidate_id="deterministic_field",
    )
    first = export_velocity_vtk_series(_field, output_dir=tmp_path / "a", **kwargs)
    second = export_velocity_vtk_series(_field, output_dir=tmp_path / "b", **kwargs)
    assert first.vtk_paths[0].read_bytes() == second.vtk_paths[0].read_bytes()
    assert first.manifest_path.read_bytes() == second.manifest_path.read_bytes()


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"x": [0.0, 0.0]}, "x must be strictly increasing"),
        ({"times": [0.5, 0.25]}, "times must be strictly increasing"),
        ({"candidate_id": "   "}, "candidate_id must be nonempty"),
    ],
)
def test_vtk_series_fails_closed_on_bad_schema(tmp_path, kwargs, message):
    call = dict(
        output_dir=tmp_path / "bad",
        x=[-1.0, 1.0],
        y=[-1.0, 1.0],
        z=[-1.0, 1.0],
        times=[0.25],
        candidate_id="candidate",
    )
    call.update(kwargs)
    with pytest.raises(ValueError, match=message):
        export_velocity_vtk_series(_field, **call)


def test_vtk_series_rejects_bad_velocity_shape_and_nonfinite(tmp_path):
    def bad_shape(points: np.ndarray, time: float) -> np.ndarray:
        return np.zeros((points.shape[0], 2))

    with pytest.raises(ValueError, match="must return shape"):
        export_velocity_vtk_series(
            bad_shape,
            output_dir=tmp_path / "shape",
            x=[-1.0, 1.0],
            y=[-1.0, 1.0],
            z=[-1.0, 1.0],
            times=[0.25],
            candidate_id="bad-shape",
        )

    def bad_value(points: np.ndarray, time: float) -> np.ndarray:
        out = np.zeros((points.shape[0], 3))
        out[0, 0] = np.nan
        return out

    with pytest.raises(ValueError, match="non-finite"):
        export_velocity_vtk_series(
            bad_value,
            output_dir=tmp_path / "nan",
            x=[-1.0, 1.0],
            y=[-1.0, 1.0],
            z=[-1.0, 1.0],
            times=[0.25],
            candidate_id="bad-value",
        )
