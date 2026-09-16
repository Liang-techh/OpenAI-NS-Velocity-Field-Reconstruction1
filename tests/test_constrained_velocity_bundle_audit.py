import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.velocity_bundle_audit import audit_velocity_bundle


class FakeField:
    def __init__(self, candidate_path):
        self.path = Path(candidate_path)
        self.sha256 = hashlib.sha256(self.path.read_bytes()).hexdigest()

    def grid(self, x, y, z, times):
        tt, xx, yy, zz = np.meshgrid(times, x, y, z, indexing="ij")
        return np.stack((xx + tt, yy - 2 * tt, zz + xx * yy), axis=-1)

    def metadata(self):
        return {
            "family": "test_velocity_v1",
            "candidate_sha256": self.sha256,
            "components": ["u", "v", "w"],
            "coordinates": "right-handed Cartesian x,y,z",
            "units": "dimensionless; no physical unit or OpenAI scene calibration claimed",
            "time_interval": [0.25, 0.75],
            "compact_support": "test-only",
            "visual_correspondence": "not verified",
            "pde_acceptance": "not evaluated",
            "grid_layout": ["time", "x", "y", "z", "component"],
        }


def make_bundle(tmp_path):
    candidate = tmp_path / "candidate.json"
    candidate.write_text('{"family":"test_velocity_v1"}\n', encoding="utf-8")
    field = FakeField(candidate)
    out = tmp_path / "bundle"
    out.mkdir()
    axes = np.linspace(-1, 1, 5)
    times = np.array([0.25, 0.5, 0.75])
    values = field.grid(axes, axes, axes, times)
    grid_path = out / "grid.npz"
    np.savez_compressed(
        grid_path,
        x=axes,
        y=axes,
        z=axes,
        times=times,
        u=values[..., 0],
        v=values[..., 1],
        w=values[..., 2],
    )
    meta = field.metadata()
    meta["grid_shape"] = list(values.shape)
    meta["grid_npz_sha256"] = hashlib.sha256(grid_path.read_bytes()).hexdigest()
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return field, out


def test_reproducible_bundle_round_trip(tmp_path):
    field, out = make_bundle(tmp_path)
    report = audit_velocity_bundle(field, out)
    assert report["status"] == "passed"
    assert report["max_abs_regeneration_error"] == 0.0
    assert report["grid_shape"] == [3, 5, 5, 5, 3]
    assert report["truth_boundary"] == {
        "velocity_export_ready": True,
        "visualization_ready": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_correspondence_verified": False,
        "note": "artifact reproducibility only; stronger claims require independent evidence",
    }


def test_tampered_grid_fails_before_regeneration(tmp_path):
    field, out = make_bundle(tmp_path)
    with (out / "grid.npz").open("ab") as handle:
        handle.write(b"tamper")
    with pytest.raises(ValueError, match="grid.npz SHA256 mismatch"):
        audit_velocity_bundle(field, out)


def test_candidate_identity_and_public_metadata_fail_closed(tmp_path):
    field, out = make_bundle(tmp_path)
    meta_path = out / "metadata.json"
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    meta["candidate_sha256"] = "0" * 64
    meta_path.write_text(json.dumps(meta), encoding="utf-8")
    with pytest.raises(ValueError, match="candidate SHA256 mismatch"):
        audit_velocity_bundle(field, out)


def test_checked_in_velocity_bundle_regenerates_through_public_api():
    from openai_ns_reconstruction.velocity_components import VelocityField

    repo_root = Path(__file__).resolve().parents[1]
    artifact_dir = repo_root / "artifacts" / "visual" / "velocity_api"
    report = audit_velocity_bundle(VelocityField(), artifact_dir)
    assert report["status"] == "passed"
    assert report["candidate_sha256"] == "a8c9e4b4536409e91c3cd06937e625f2085c475ad7f8c717ccd23c83d804e381"
    assert report["grid_npz_sha256"] == "3ac3a589ab2a94a78621ec3c93f8e6881681c6b8c80c09beb4b58420ee9a7016"
    assert report["grid_shape"] == [3, 17, 17, 17, 3]
    assert report["reference_times"] == [0.25, 0.5, 0.75]
