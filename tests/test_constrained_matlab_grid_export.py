from pathlib import Path

import numpy as np
import pytest
from scipy.io import loadmat, savemat

from openai_ns_reconstruction.constrained_matlab_grid_export import (
    export_velocity_grid_mat,
    load_velocity_grid_mat,
)


SHA = "1" * 64


def sample_grid():
    x = np.array([-1.0, 0.0, 1.0])
    y = np.array([-0.5, 0.5])
    z = np.array([-2.0, -0.25, 0.75, 2.0])
    t = np.array([0.25, 0.5, 0.75])
    tt, xx, yy, zz = np.meshgrid(t, x, y, z, indexing="ij")
    velocity = np.stack(
        (
            xx + 2.0 * tt,
            yy - tt + 0.1 * zz,
            zz + 0.25 * xx + tt,
        ),
        axis=-1,
    )
    return x, y, z, t, velocity


def export(tmp_path: Path, **kwargs):
    x, y, z, t, velocity = sample_grid()
    args = dict(
        x=x,
        y=y,
        z=z,
        t=t,
        velocity=velocity,
        candidate_sha256=SHA,
        candidate_provenance="candidate:test-fixture",
        grid_provenance="grid:fixed-cartesian-fixture",
    )
    args.update(kwargs)
    path = tmp_path / "field.mat"
    receipt = export_velocity_grid_mat(path, **args)
    return path, receipt, args


def test_roundtrip_preserves_axes_velocity_and_identity(tmp_path):
    path, receipt, args = export(tmp_path)
    loaded = load_velocity_grid_mat(path)
    assert receipt.roundtrip_verified is True
    assert receipt.velocity_shape == args["velocity"].shape
    assert np.array_equal(loaded.x, args["x"])
    assert np.array_equal(loaded.y, args["y"])
    assert np.array_equal(loaded.z, args["z"])
    assert np.array_equal(loaded.t, args["t"])
    assert np.array_equal(loaded.velocity, args["velocity"])
    assert loaded.candidate_sha256 == SHA


def test_mat_variables_are_matlab_facing_txyz_components(tmp_path):
    path, _, args = export(tmp_path)
    raw = loadmat(path, appendmat=False, squeeze_me=False, chars_as_strings=True)
    expected = args["velocity"].shape[:-1]
    assert raw["u_txyz"].shape == expected
    assert raw["v_txyz"].shape == expected
    assert raw["w_txyz"].shape == expected
    assert np.array_equal(raw["u_txyz"], args["velocity"][..., 0])
    assert np.array_equal(raw["v_txyz"], args["velocity"][..., 1])
    assert np.array_equal(raw["w_txyz"], args["velocity"][..., 2])
    assert raw["x"].shape == (3, 1)
    assert raw["t"].shape == (3, 1)


def test_truth_flags_remain_false_and_arrays_readonly(tmp_path):
    path, receipt, _ = export(tmp_path)
    loaded = load_velocity_grid_mat(path)
    assert all(value is False for value in receipt.truth_flags.values())
    assert all(value is False for value in loaded.truth_flags.values())
    assert loaded.velocity.flags.writeable is False
    assert loaded.x.flags.writeable is False
    with pytest.raises(ValueError):
        loaded.velocity[0, 0, 0, 0, 0] = 9.0


def test_diagnostics_match_direct_speed(tmp_path):
    path, receipt, args = export(tmp_path)
    loaded = load_velocity_grid_mat(path)
    speed = np.linalg.norm(args["velocity"], axis=-1)
    expected_max = np.max(speed, axis=(1, 2, 3))
    expected_rms = np.sqrt(np.mean(speed * speed, axis=(1, 2, 3)))
    assert np.allclose(loaded.per_time_max_speed, expected_max)
    assert np.allclose(loaded.per_time_rms_speed, expected_rms)
    assert receipt.global_max_speed == pytest.approx(float(np.max(speed)))


def test_refuses_overwrite_without_explicit_permission(tmp_path):
    path, _, args = export(tmp_path)
    with pytest.raises(FileExistsError):
        export_velocity_grid_mat(path, **args)
    receipt = export_velocity_grid_mat(path, **args, overwrite=True)
    assert receipt.roundtrip_verified is True


@pytest.mark.parametrize("axis", ["x", "y", "z", "t"])
def test_rejects_nonmonotone_axes(tmp_path, axis):
    _, _, args = export(tmp_path)
    values = np.array(args[axis], copy=True)
    if values.size == 1:
        pytest.skip("fixture axis has one entry")
    values[1] = values[0]
    args[axis] = values
    with pytest.raises(ValueError, match="strictly increasing"):
        export_velocity_grid_mat(tmp_path / f"bad-{axis}.mat", **args)


def test_rejects_nonfinite_axis(tmp_path):
    x, y, z, t, velocity = sample_grid()
    x[1] = np.nan
    with pytest.raises(ValueError, match="finite"):
        export_velocity_grid_mat(
            tmp_path / "bad.mat",
            x=x, y=y, z=z, t=t, velocity=velocity,
            candidate_sha256=SHA,
            candidate_provenance="candidate:x",
            grid_provenance="grid:x",
        )


def test_rejects_wrong_velocity_shape(tmp_path):
    x, y, z, t, velocity = sample_grid()
    with pytest.raises(ValueError, match="velocity has shape"):
        export_velocity_grid_mat(
            tmp_path / "bad.mat",
            x=x, y=y, z=z, t=t, velocity=velocity[..., :2],
            candidate_sha256=SHA,
            candidate_provenance="candidate:x",
            grid_provenance="grid:x",
        )


def test_rejects_nonfinite_velocity(tmp_path):
    x, y, z, t, velocity = sample_grid()
    velocity[0, 0, 0, 0, 0] = np.inf
    with pytest.raises(ValueError, match="finite"):
        export_velocity_grid_mat(
            tmp_path / "bad.mat",
            x=x, y=y, z=z, t=t, velocity=velocity,
            candidate_sha256=SHA,
            candidate_provenance="candidate:x",
            grid_provenance="grid:x",
        )


def test_rejects_zero_time_slice(tmp_path):
    x, y, z, t, velocity = sample_grid()
    velocity[1, ...] = 0.0
    with pytest.raises(ValueError, match="time slice"):
        export_velocity_grid_mat(
            tmp_path / "bad.mat",
            x=x, y=y, z=z, t=t, velocity=velocity,
            candidate_sha256=SHA,
            candidate_provenance="candidate:x",
            grid_provenance="grid:x",
        )


@pytest.mark.parametrize("sha", ["", "A" * 64, "1" * 63, "z" * 64])
def test_rejects_bad_candidate_sha(tmp_path, sha):
    x, y, z, t, velocity = sample_grid()
    with pytest.raises(ValueError):
        export_velocity_grid_mat(
            tmp_path / "bad.mat",
            x=x, y=y, z=z, t=t, velocity=velocity,
            candidate_sha256=sha,
            candidate_provenance="candidate:x",
            grid_provenance="grid:x",
        )


@pytest.mark.parametrize("key", ["candidate_provenance", "grid_provenance"])
def test_rejects_missing_provenance(tmp_path, key):
    _, _, args = export(tmp_path)
    args[key] = "   "
    with pytest.raises(ValueError, match="non-empty"):
        export_velocity_grid_mat(tmp_path / f"bad-{key}.mat", **args)


def test_rejects_non_mat_suffix_and_missing_parent(tmp_path):
    x, y, z, t, velocity = sample_grid()
    common = dict(
        x=x, y=y, z=z, t=t, velocity=velocity,
        candidate_sha256=SHA,
        candidate_provenance="candidate:x",
        grid_provenance="grid:x",
    )
    with pytest.raises(ValueError, match="end in .mat"):
        export_velocity_grid_mat(tmp_path / "field.npz", **common)
    with pytest.raises(FileNotFoundError):
        export_velocity_grid_mat(tmp_path / "missing" / "field.mat", **common)


def test_loader_rejects_truth_promotion(tmp_path):
    path, _, _ = export(tmp_path)
    raw = loadmat(path, appendmat=False, squeeze_me=False, chars_as_strings=True)
    clean = {key: value for key, value in raw.items() if not key.startswith("__")}
    clean["truth_visualization_ready"] = np.array([[1]], dtype=np.uint8)
    savemat(path, clean, appendmat=False, format="5")
    with pytest.raises(ValueError, match="forbidden"):
        load_velocity_grid_mat(path)


def test_loader_rejects_layout_mutation(tmp_path):
    path, _, _ = export(tmp_path)
    raw = loadmat(path, appendmat=False, squeeze_me=False, chars_as_strings=True)
    clean = {key: value for key, value in raw.items() if not key.startswith("__")}
    clean["layout"] = "x,y,z,time,component"
    savemat(path, clean, appendmat=False, format="5")
    with pytest.raises(ValueError, match="layout"):
        load_velocity_grid_mat(path)
