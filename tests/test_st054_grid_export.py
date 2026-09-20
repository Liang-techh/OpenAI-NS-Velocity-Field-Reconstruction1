from __future__ import annotations

import json

import numpy as np
import pytest
from scipy.io import loadmat

from openai_ns_reconstruction.st054_grid_export import (
    _ARRAY_ORDER,
    _canonical_payload_sha256,
    export_st054_grid,
    sample_st054_grid,
    verify_st054_grid_export,
)
from openai_ns_reconstruction.st054_snapshot import ST054_REFERENCE_ATOL


@pytest.mark.parametrize("model_id", ["ST054-Q2", "ST054-M3"])
def test_small_grid_npz_mat_exact_and_callable_bound(tmp_path, model_id):
    receipt = export_st054_grid(
        tmp_path,
        model_id,
        grid_size=5,
        times=(0.25, 0.5, 0.75),
    )
    assert receipt["passed"] is True
    assert receipt["pde_validated"] is False
    assert receipt["visualization_ready"] is False
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["callable_node_max_component_error"] <= ST054_REFERENCE_ATOL

    stem = model_id.lower().replace("-", "_") + "_velocity_grid"
    npz_path = tmp_path / f"{stem}.npz"
    mat_path = tmp_path / f"{stem}.mat"
    manifest_path = tmp_path / f"{stem}.json"
    manifest = json.loads(manifest_path.read_text())
    assert manifest["grid_size"] == 5
    assert manifest["times"] == [0.25, 0.5, 0.75]
    assert manifest["sampling_only"] is True
    assert manifest["grid_interpolation_used"] is False
    assert manifest["truth_boundary"]["pde_validated"] is False

    with np.load(npz_path, allow_pickle=False) as npz:
        arrays = {name: np.asarray(npz[name]) for name in _ARRAY_ORDER}
    mat = loadmat(mat_path, simplify_cells=True)
    for name in _ARRAY_ORDER:
        right = np.asarray(mat[name]).squeeze() if name in ("x", "y", "z", "t") else np.asarray(mat[name])
        assert np.array_equal(arrays[name], right)
    assert arrays["u"].shape == (3, 5, 5, 5)
    assert arrays["v"].shape == arrays["u"].shape
    assert arrays["w"].shape == arrays["u"].shape
    assert _canonical_payload_sha256(arrays) == manifest["canonical_payload_sha256"]


def test_sample_rejects_noncanonical_time_requests():
    with pytest.raises(ValueError, match=r"\[0.25, 0.75\]"):
        sample_st054_grid(grid_size=3, times=(0.2, 0.5))
    with pytest.raises(ValueError, match="strictly increasing"):
        sample_st054_grid(grid_size=3, times=(0.5, 0.5))
    with pytest.raises(ValueError, match="integer >= 3"):
        sample_st054_grid(grid_size=2, times=(0.5,))


def test_export_refuses_overwrite_and_truth_boundary_tamper(tmp_path):
    receipt = export_st054_grid(
        tmp_path,
        "ST054-Q2",
        grid_size=3,
        times=(0.25, 0.75),
    )
    assert receipt["passed"] is True
    with pytest.raises(FileExistsError):
        export_st054_grid(
            tmp_path,
            "ST054-Q2",
            grid_size=3,
            times=(0.25, 0.75),
        )

    stem = "st054_q2_velocity_grid"
    manifest_path = tmp_path / f"{stem}.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["truth_boundary"]["pde_validated"] = True
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="truth boundary"):
        verify_st054_grid_export(
            tmp_path / f"{stem}.npz",
            tmp_path / f"{stem}.mat",
            manifest_path,
        )
