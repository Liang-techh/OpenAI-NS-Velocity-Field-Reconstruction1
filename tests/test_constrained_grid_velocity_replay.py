import numpy as np
import pytest

from openai_ns_reconstruction.constrained_grid_velocity_replay import (
    GridVelocityReplay,
    metadata_json,
)


SHA = "1" * 64
PROV = "unit-test frozen velocity grid"


def analytic(t, x, y, z):
    return np.stack(
        (
            1.0 + t + 2.0 * x - 3.0 * y + 0.5 * z,
            -0.5 + 2.0 * t - x + 0.25 * y + z,
            0.75 - t + 0.5 * x + 2.0 * y - 1.5 * z,
        ),
        axis=-1,
    )


def field():
    t = np.array([0.25, 0.5, 0.75])
    x = np.array([-1.0, 0.0, 1.0])
    y = np.array([-1.0, 0.5, 1.0])
    z = np.array([-1.0, 0.25, 1.0])
    T, X, Y, Z = np.meshgrid(t, x, y, z, indexing="ij")
    values = analytic(T, X, Y, Z)
    return GridVelocityReplay(
        x=x,
        y=y,
        z=z,
        t=t,
        velocity_txyzc=values,
        candidate_sha256=SHA,
        provenance=PROV,
    )


def test_linear_replay_is_exact_inside_cells():
    f = field()
    x = np.array([-0.7, 0.2, 0.8])
    y = np.array([-0.4, 0.0, 0.7])
    z = np.array([-0.2, 0.4, 0.9])
    t = np.array([0.30, 0.60, 0.70])
    np.testing.assert_allclose(f.velocity(x, y, z, t), analytic(t, x, y, z), rtol=0, atol=2e-14)


def test_grid_nodes_component_order_and_scalar_shape():
    f = field()
    got = f.velocity(0.0, 0.5, 0.25, 0.5)
    assert got.shape == (3,)
    np.testing.assert_allclose(got, analytic(0.5, 0.0, 0.5, 0.25))
    assert f.metadata()["component_order"] == ["u", "v", "w"]
    assert f.metadata()["layout"] == "time,x,y,z,component"


def test_at_points_broadcasts_time_and_matches_public_signature():
    f = field()
    points = np.array([[-0.7, -0.4, -0.2], [0.2, 0.0, 0.4], [0.8, 0.7, 0.9]])
    times = np.array([0.30, 0.60, 0.70])
    np.testing.assert_allclose(f.at_points(points, times), analytic(times, points[:, 0], points[:, 1], points[:, 2]))


def test_save_load_round_trip_and_readonly_arrays(tmp_path):
    f = field()
    out = f.save(tmp_path / "velocity_grid.npz")
    g = GridVelocityReplay.load(out)
    assert g.grid_sha256 == f.grid_sha256
    assert g.candidate_sha256 == SHA
    assert g.provenance == PROV
    np.testing.assert_array_equal(g.velocity_txyzc, f.velocity_txyzc)
    assert not g.x.flags.writeable
    assert not g.velocity_txyzc.flags.writeable
    np.testing.assert_allclose(g.velocity(0.2, 0.0, 0.4, 0.6), f.velocity(0.2, 0.0, 0.4, 0.6))
    with pytest.raises(FileExistsError):
        f.save(out)


def test_tampered_grid_sha_fails_closed(tmp_path):
    f = field()
    path = f.save(tmp_path / "grid.npz")
    with np.load(path, allow_pickle=False) as data:
        payload = {k: np.array(data[k], copy=True) for k in data.files}
    payload["grid_sha256"] = np.asarray("0" * 64)
    with path.open("wb") as handle:
        np.savez_compressed(handle, **payload)
    with pytest.raises(ValueError, match="grid_sha256 mismatch"):
        GridVelocityReplay.load(path)


def test_outside_grid_and_nonfinite_query_fail_closed():
    f = field()
    with pytest.raises(ValueError, match="extrapolation is forbidden"):
        f.velocity(1.01, 0.0, 0.0, 0.5)
    with pytest.raises(ValueError, match="finite"):
        f.velocity(np.nan, 0.0, 0.0, 0.5)


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda kw: kw.update(x=[0.0, -1.0]), "strictly increasing"),
        (lambda kw: kw.update(candidate_sha256="BAD"), "64 lowercase"),
        (lambda kw: kw.update(provenance="  "), "non-empty"),
        (lambda kw: kw.update(velocity_txyzc=np.zeros((3, 3, 3, 3, 3))), "nontrivial"),
        (lambda kw: kw.update(velocity_txyzc=np.zeros((3, 3, 3, 3, 2))), "shape"),
    ],
)
def test_invalid_inputs_fail_closed(change, message):
    f = field()
    kw = dict(
        x=f.x,
        y=f.y,
        z=f.z,
        t=f.t,
        velocity_txyzc=f.velocity_txyzc,
        candidate_sha256=SHA,
        provenance=PROV,
    )
    change(kw)
    with pytest.raises(ValueError, match=message):
        GridVelocityReplay(**kw)


def test_nonfinite_grid_fails_closed():
    f = field()
    bad = np.array(f.velocity_txyzc, copy=True)
    bad[0, 0, 0, 0, 0] = np.inf
    with pytest.raises(ValueError, match="finite"):
        GridVelocityReplay(
            x=f.x, y=f.y, z=f.z, t=f.t, velocity_txyzc=bad,
            candidate_sha256=SHA, provenance=PROV,
        )


def test_truth_boundary_is_explicit_and_json_stable():
    f = field()
    m = f.metadata()
    assert m["velocity_export_ready"] is True
    assert m["visualization_interpolation_only"] is True
    assert m["visualization_ready"] is False
    assert m["visual_correspondence_verified"] is False
    assert m["pde_validated"] is False
    assert m["interpolated_derivatives_valid_for_pde_acceptance"] is False
    assert m["openai_field_identified"] is False
    assert m["extrapolation_allowed"] is False
    text = metadata_json(f)
    assert text == metadata_json(f)
    assert '"pde_validated":false' in text
