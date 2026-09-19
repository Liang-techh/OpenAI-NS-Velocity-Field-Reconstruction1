from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.st052_grid_export import (
    MANIFEST_FILENAME,
    MAT_FILENAME,
    NPZ_FILENAME,
    export_velocity_grid,
    sample_velocity_grid,
    verify_velocity_grid_export,
)


class FakeWholeCandidate:
    candidate_id = "ST052-M-linear-temporal-child-v1"
    identity_sha256 = "a" * 64
    manifest = {
        "truth_boundary": {
            "whole_child_bundle_materialized": True,
            "standalone_package_parent_runtime_ready": False,
            "held_out_temporal_child_pde_residual_evaluated": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "production_candidate_selected": False,
            "pde_validated": False,
            "source_correspondence_verified": False,
        }
    }

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.stack((x + t, y - 2.0 * t, z + x * y + 0.5 * t), axis=-1)


def test_sample_layout_is_time_xyz_components():
    arrays = sample_velocity_grid(
        FakeWholeCandidate(), grid_points=5, times=(0.25, 0.5, 0.75), domain=(-1.0, 1.0)
    )
    assert arrays["u"].shape == (3, 5, 5, 5)
    assert arrays["v"].shape == arrays["u"].shape
    assert arrays["w"].shape == arrays["u"].shape
    assert np.allclose(arrays["u"][0], arrays["x"][:, None, None] + 0.25)


def test_npz_and_mat_roundtrip_exactly(tmp_path: Path):
    manifest = export_velocity_grid(
        FakeWholeCandidate(), tmp_path, grid_points=5, times=(0.25, 0.5, 0.75)
    )
    checked = verify_velocity_grid_export(tmp_path)
    assert checked["npz_mat_exact_equal"] is True
    assert checked["npz_payload_sha256"] == manifest["grid_payload_sha256"]
    assert checked["mat_payload_sha256"] == manifest["grid_payload_sha256"]
    assert manifest["grid"]["velocity_component_shape"] == [3, 5, 5, 5]
    assert manifest["truth_boundary"]["grid_export_materialized"] is True
    assert manifest["truth_boundary"]["visualization_ready"] is False
    assert manifest["truth_boundary"]["pde_validated"] is False
    assert manifest["execution_environment"]["dependency_runtime_identity_closed"] is False
    assert manifest["matlab_compatibility"]["actual_matlab_runtime_executed"] is False


def test_raw_file_tamper_fails_closed(tmp_path: Path):
    export_velocity_grid(FakeWholeCandidate(), tmp_path, grid_points=4, times=(0.25, 0.75))
    target = tmp_path / NPZ_FILENAME
    target.write_bytes(target.read_bytes() + b"tamper")
    with pytest.raises(ValueError, match="checksum mismatch"):
        verify_velocity_grid_export(tmp_path)


def test_manifest_payload_tamper_fails_closed(tmp_path: Path):
    export_velocity_grid(FakeWholeCandidate(), tmp_path, grid_points=4, times=(0.25, 0.75))
    path = tmp_path / MANIFEST_FILENAME
    data = json.loads(path.read_text())
    data["grid_payload_sha256"] = "0" * 64
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="payload checksum mismatch"):
        verify_velocity_grid_export(tmp_path)


def test_truth_boundary_is_required(tmp_path: Path):
    class Bad(FakeWholeCandidate):
        manifest = {"truth_boundary": {**FakeWholeCandidate.manifest["truth_boundary"], "pde_validated": True}}

    with pytest.raises(ValueError, match="pde_validated=false"):
        export_velocity_grid(Bad(), tmp_path, grid_points=3, times=(0.25, 0.75))


def test_refuses_overwrite(tmp_path: Path):
    export_velocity_grid(FakeWholeCandidate(), tmp_path, grid_points=3, times=(0.25, 0.75))
    assert (tmp_path / NPZ_FILENAME).is_file()
    assert (tmp_path / MAT_FILENAME).is_file()
    with pytest.raises(ValueError, match="refusing to overwrite"):
        export_velocity_grid(FakeWholeCandidate(), tmp_path, grid_points=3, times=(0.25, 0.75))
