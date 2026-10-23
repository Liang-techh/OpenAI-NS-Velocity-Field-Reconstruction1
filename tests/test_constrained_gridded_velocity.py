import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_gridded_velocity import StructuredVelocityField


def affine_velocity(points, time):
    p = np.asarray(points, dtype=float)
    x, y, z = p.T
    t = float(time)
    return np.column_stack((
        1.0 + 2.0*t + 3.0*x - y + 0.5*z,
        -t + x + 2.0*y,
        4.0 + z - 2.0*x,
    ))


def make_field():
    return StructuredVelocityField.sample(
        affine_velocity,
        x=[-1.0, -0.2, 0.7, 1.3],
        y=[-0.8, 0.1, 1.1],
        z=[-1.2, -0.3, 0.4, 1.0],
        times=[0.25, 0.5, 0.75],
        candidate_id="manufactured-affine-v1",
    )


def test_linear_space_time_interpolation_is_exact_for_affine_field():
    field = make_field()
    rng = np.random.default_rng(90617)
    points = np.column_stack((
        rng.uniform(-1.0, 1.3, 257),
        rng.uniform(-0.8, 1.1, 257),
        rng.uniform(-1.2, 1.0, 257),
    ))
    for time in (0.31, 0.43, 0.69):
        got = field.at_points(points, time)
        expected = affine_velocity(points, time)
        np.testing.assert_allclose(got, expected, rtol=0.0, atol=4e-15)


def test_npz_round_trip_preserves_callable_grid_and_truth_boundary(tmp_path):
    field = make_field()
    path = tmp_path / "velocity_grid.npz"
    digest = field.save_npz(path)
    assert len(digest) == 64

    loaded = StructuredVelocityField.load_npz(path)
    assert loaded.candidate_id == "manufactured-affine-v1"
    assert loaded.metadata["claim_scope"] == "visualization_interpolation_only"
    assert loaded.metadata["velocity_export_ready"] is True
    assert loaded.metadata["visualization_ready"] is False
    assert loaded.metadata["pde_validated"] is False
    assert loaded.metadata["paper_exact"] is False
    assert loaded.metadata["openai_field_identified"] is False

    xq = [-0.9, 0.0, 1.0]
    yq = [-0.7, 0.4]
    zq = [-1.0, 0.2, 0.8]
    tq = [0.3, 0.7]
    grid = loaded.grid(xq, yq, zq, tq)
    assert grid.shape == (2, 3, 2, 3, 3)
    X, Y, Z = np.meshgrid(xq, yq, zq, indexing="ij")
    pts = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    for i, t in enumerate(tq):
        np.testing.assert_allclose(grid[i].reshape(-1, 3), affine_velocity(pts, t), rtol=0.0, atol=4e-15)


def test_fail_closed_on_bad_axes_values_bounds_and_metadata(tmp_path):
    with pytest.raises(ValueError, match="strictly increasing"):
        StructuredVelocityField.sample(
            affine_velocity,
            x=[-1.0, 0.0, 0.0], y=[0.0, 1.0], z=[0.0, 1.0], times=[0.25, 0.75],
            candidate_id="bad",
        )

    field = make_field()
    with pytest.raises(ValueError):
        field.at_points(np.array([[2.0, 0.0, 0.0]]), 0.5)
    with pytest.raises(ValueError):
        field.at_points(np.array([[0.0, 0.0, 0.0]]), 0.9)

    path = tmp_path / "tampered.npz"
    meta = dict(field.metadata)
    meta["pde_validated"] = True
    np.savez_compressed(
        path,
        x=field.x,
        y=field.y,
        z=field.z,
        times=field.times,
        values=field.values,
        metadata_json=np.asarray(json.dumps(meta)),
    )
    with pytest.raises(ValueError, match="pde_validated"):
        StructuredVelocityField.load_npz(path)

    def bad_velocity(points, time):
        out = affine_velocity(points, time)
        out[0, 0] = np.nan
        return out

    with pytest.raises(FloatingPointError, match="non-finite"):
        StructuredVelocityField.sample(
            bad_velocity,
            x=[-1.0, 0.0], y=[-1.0, 0.0], z=[-1.0, 0.0], times=[0.25, 0.75],
            candidate_id="nan-field",
        )
