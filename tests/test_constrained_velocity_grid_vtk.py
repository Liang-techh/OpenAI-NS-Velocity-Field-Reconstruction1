import json

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_velocity_grid_vtk import (
    export_velocity_grid_vtk,
    load_velocity_grid_vtk_manifest,
)


SHA = "a" * 64
PROVENANCE = {"task_id": "CR-A9-027", "candidate": "synthetic-affine-regression"}


def sample_grid():
    x = np.array([-1.0, 0.5, 2.0])
    y = np.array([-0.5, 1.0])
    z = np.array([0.25, 1.25])
    t = np.array([0.25, 0.75])
    tt, xx, yy, zz = np.meshgrid(t, x, y, z, indexing="ij")
    velocity = np.stack(
        [1.0 + xx + 2.0 * tt, -2.0 + 3.0 * yy - tt, 0.5 + 4.0 * zz + xx],
        axis=-1,
    )
    return x, y, z, t, velocity


def test_export_preserves_axes_components_and_vtk_point_order(tmp_path):
    x, y, z, t, velocity = sample_grid()
    out = tmp_path / "vtk"
    manifest = export_velocity_grid_vtk(
        out,
        x=x,
        y=y,
        z=z,
        t=t,
        velocity_txyz=velocity,
        candidate_sha256=SHA,
        provenance=PROVENANCE,
    )
    assert manifest["dataset"] == "RECTILINEAR_GRID"
    assert manifest["component_order"] == ["u", "v", "w"]
    assert manifest["vtk_point_order"] == "x-fastest, then y, then z"
    lines = (out / "velocity_t000.vtk").read_text(encoding="ascii").splitlines()
    assert "DIMENSIONS 3 2 2" in lines
    assert "POINT_DATA 12" in lines
    start = lines.index("VECTORS velocity double") + 1
    rows = np.array([[float(v) for v in line.split()] for line in lines[start : start + 12]])
    expected = np.transpose(velocity[0], (2, 1, 0, 3)).reshape(-1, 3)
    np.testing.assert_array_equal(rows, expected)


def test_manifest_roundtrip_keeps_false_scientific_claims(tmp_path):
    x, y, z, t, velocity = sample_grid()
    out = tmp_path / "vtk"
    written = export_velocity_grid_vtk(
        out,
        x=x,
        y=y,
        z=z,
        t=t,
        velocity_txyz=velocity,
        candidate_sha256=SHA,
        provenance=PROVENANCE,
    )
    loaded = load_velocity_grid_vtk_manifest(out)
    assert loaded == written
    assert len(loaded["grid_sha256"]) == 64
    assert loaded["truth_boundary"]["velocity_export_ready"] is True
    assert loaded["truth_boundary"]["registered_nontriviality_gate_assessed"] is False
    assert loaded["truth_boundary"]["vtk_derivatives_valid_for_pde_acceptance"] is False
    for key in (
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert loaded["truth_boundary"][key] is False


def test_tampered_vtk_file_fails_closed(tmp_path):
    x, y, z, t, velocity = sample_grid()
    out = tmp_path / "vtk"
    export_velocity_grid_vtk(out, x=x, y=y, z=z, t=t, velocity_txyz=velocity, candidate_sha256=SHA, provenance=PROVENANCE)
    target = out / "velocity_t001.vtk"
    target.write_text(target.read_text(encoding="ascii") + "# tamper\n", encoding="ascii")
    with pytest.raises(ValueError, match="integrity check failed"):
        load_velocity_grid_vtk_manifest(out)


def test_promoted_pde_claim_fails_closed(tmp_path):
    x, y, z, t, velocity = sample_grid()
    out = tmp_path / "vtk"
    export_velocity_grid_vtk(out, x=x, y=y, z=z, t=t, velocity_txyz=velocity, candidate_sha256=SHA, provenance=PROVENANCE)
    path = out / "manifest.json"
    manifest = json.loads(path.read_text())
    manifest["truth_boundary"]["pde_validated"] = True
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="pde_validated"):
        load_velocity_grid_vtk_manifest(out)


def test_exact_zero_time_slice_fails_closed(tmp_path):
    x, y, z, t, velocity = sample_grid()
    velocity[0] = 0.0
    with pytest.raises(ValueError, match="exact-zero"):
        export_velocity_grid_vtk(tmp_path / "vtk", x=x, y=y, z=z, t=t, velocity_txyz=velocity, candidate_sha256=SHA, provenance=PROVENANCE)


@pytest.mark.parametrize(
    "change, message",
    [
        (lambda x, y, z, t, v: (x[::-1], y, z, t, v), "strictly increasing"),
        (lambda x, y, z, t, v: (x, np.array([-0.5, np.nan]), z, t, v), "finite"),
        (lambda x, y, z, t, v: (x, y, z, t, v[..., :2]), "shape"),
    ],
)
def test_bad_axes_or_shape_fail_closed(tmp_path, change, message):
    x, y, z, t, velocity = change(*sample_grid())
    with pytest.raises(ValueError, match=message):
        export_velocity_grid_vtk(tmp_path / "vtk", x=x, y=y, z=z, t=t, velocity_txyz=velocity, candidate_sha256=SHA, provenance=PROVENANCE)


def test_bad_identity_nonfinite_velocity_and_existing_destination_fail_closed(tmp_path):
    x, y, z, t, velocity = sample_grid()
    bad = velocity.copy()
    bad[0, 0, 0, 0, 0] = np.inf
    with pytest.raises(ValueError, match="finite"):
        export_velocity_grid_vtk(tmp_path / "a", x=x, y=y, z=z, t=t, velocity_txyz=bad, candidate_sha256=SHA, provenance=PROVENANCE)
    with pytest.raises(ValueError, match="64 lowercase"):
        export_velocity_grid_vtk(tmp_path / "b", x=x, y=y, z=z, t=t, velocity_txyz=velocity, candidate_sha256="bad", provenance=PROVENANCE)
    with pytest.raises(ValueError, match="non-empty"):
        export_velocity_grid_vtk(tmp_path / "c", x=x, y=y, z=z, t=t, velocity_txyz=velocity, candidate_sha256=SHA, provenance={})
    existing = tmp_path / "d"
    existing.mkdir()
    with pytest.raises(FileExistsError, match="refusing"):
        export_velocity_grid_vtk(existing, x=x, y=y, z=z, t=t, velocity_txyz=velocity, candidate_sha256=SHA, provenance=PROVENANCE)
