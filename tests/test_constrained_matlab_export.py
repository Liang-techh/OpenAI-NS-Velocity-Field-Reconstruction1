import numpy as np
import pytest
from scipy.io import loadmat, savemat

from openai_ns_reconstruction.constrained_matlab_export import (
    export_velocity_mat,
    load_velocity_mat,
    sample_velocity_grid,
)


def manufactured_velocity(points, time):
    x, y, z = points.T
    return np.column_stack((x + time, y - 2.0 * time, z + x * y))


def test_matlab_round_trip_preserves_grid_and_velocity(tmp_path):
    x = np.linspace(-1.0, 1.0, 5)
    y = np.linspace(-0.8, 0.8, 4)
    z = np.linspace(-0.6, 0.6, 3)
    times = np.array([0.25, 0.50, 0.75])

    path = export_velocity_mat(
        tmp_path / "velocity_bundle",
        manufactured_velocity,
        x=x,
        y=y,
        z=z,
        times=times,
        candidate_id="manufactured-smoke",
        metadata={"task_id": "CR-A9-003"},
    )
    loaded = load_velocity_mat(path)
    expected = sample_velocity_grid(manufactured_velocity, x, y, z, times)

    assert path.suffix == ".mat"
    assert loaded["candidate_id"] == "manufactured-smoke"
    assert loaded["claim_scope"] == "visualization_export_only"
    for name in ("x", "y", "z", "times", "U", "V", "W"):
        np.testing.assert_allclose(loaded[name], expected[name], rtol=0.0, atol=0.0)

    raw = loadmat(path, squeeze_me=True)
    assert int(raw["pde_validated"]) == 0
    assert int(raw["paper_exact"]) == 0
    assert int(raw["openai_field_identified"]) == 0
    assert str(raw["meta_task_id"]) == "CR-A9-003"


def test_export_matches_matlab_meshgrid_axis_order(tmp_path):
    x = np.array([-1.0, 2.0])
    y = np.array([-3.0, 4.0, 5.0])
    z = np.array([0.1, 0.7])
    times = np.array([0.5])

    path = export_velocity_mat(
        tmp_path / "axes.mat",
        manufactured_velocity,
        x=x,
        y=y,
        z=z,
        times=times,
        candidate_id="axis-smoke",
    )
    loaded = load_velocity_mat(path)

    X, Y, Z = np.meshgrid(x, y, z, indexing="xy")
    np.testing.assert_allclose(loaded["U"][..., 0], X + 0.5)
    np.testing.assert_allclose(loaded["V"][..., 0], Y - 1.0)
    np.testing.assert_allclose(loaded["W"][..., 0], Z + X * Y)


def test_fail_closed_on_bad_inputs_and_tampered_bundle(tmp_path):
    with pytest.raises(ValueError, match="strictly increasing"):
        sample_velocity_grid(
            manufactured_velocity,
            np.array([0.0, 0.0]),
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([0.5]),
        )

    def bad_velocity(points, time):
        return np.zeros((points.shape[0], 2))

    with pytest.raises(ValueError, match="must return shape"):
        sample_velocity_grid(
            bad_velocity,
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([0.0, 1.0]),
            np.array([0.5]),
        )

    good = export_velocity_mat(
        tmp_path / "good.mat",
        manufactured_velocity,
        x=np.array([0.0, 1.0]),
        y=np.array([0.0, 1.0]),
        z=np.array([0.0, 1.0]),
        times=np.array([0.5]),
        candidate_id="tamper-smoke",
    )
    raw = loadmat(good, squeeze_me=False)
    raw["axis_order"] = "wrong"
    tampered = tmp_path / "tampered.mat"
    savemat(tampered, {k: v for k, v in raw.items() if not k.startswith("__")})
    with pytest.raises(ValueError, match="axis_order"):
        load_velocity_mat(tampered)
