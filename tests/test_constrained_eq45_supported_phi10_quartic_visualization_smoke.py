import json
from pathlib import Path

import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_supported_phi10_quartic_visualization_smoke import (
    governed_quartic_candidate,
    meridional_velocity_slice,
    write_quartic_visualization_smoke,
)


EXPECTED_QUARTIC_SHA = "03fdae73170469ae1160297b489b1b83a27adddfc941119e116a21b98bd15049"


def test_quartic_visualization_smoke_writes_public_slice_bundle(tmp_path: Path):
    manifest = write_quartic_visualization_smoke(
        tmp_path,
        times=(0.25,),
        grid_size=17,
        extent=2.0,
        write_png=True,
    )

    assert manifest["candidate_sha256"] == EXPECTED_QUARTIC_SHA
    assert manifest["times"] == [0.25]
    assert manifest["truth_boundary"]["visualization_smoke_generated"] is True
    assert manifest["truth_boundary"]["visualization_ready"] is False
    assert manifest["truth_boundary"]["visual_correspondence_verified"] is False
    assert manifest["truth_boundary"]["pde_validated"] is False

    saved_manifest = json.loads(
        (tmp_path / "quartic_eq45_visualization_smoke_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    assert saved_manifest == manifest

    entry = manifest["files"][0]
    npz_path = tmp_path / entry["npz"]
    png_path = tmp_path / entry["png"]
    assert npz_path.is_file()
    assert png_path.is_file()
    assert png_path.stat().st_size > 0

    with np.load(npz_path) as data:
        assert str(data["candidate_sha256"]) == EXPECTED_QUARTIC_SHA
        assert float(data["time"]) == 0.25
        for name in ("u", "v", "w", "speed", "poloidal_speed"):
            assert data[name].shape == (17, 17)
            assert np.all(np.isfinite(data[name]))
        expected_speed = np.sqrt(data["u"] ** 2 + data["v"] ** 2 + data["w"] ** 2)
        assert np.allclose(data["speed"], expected_speed, rtol=0.0, atol=1.0e-14)


def test_quartic_midpoint_slice_is_exact_static_supported_snapshot():
    candidate = governed_quartic_candidate()
    data = meridional_velocity_slice(candidate, 0.50, grid_size=9, extent=1.25)

    x = np.asarray(data["x"], dtype=float)
    z = np.asarray(data["z"], dtype=float)
    xx, zz = np.meshgrid(x, z, indexing="xy")
    yy = np.zeros_like(xx)
    expected = candidate.base.velocity_xyz(xx, yy, zz, 0.50)
    actual = np.stack((data["u"], data["v"], data["w"]), axis=-1)

    assert candidate.sha256 == EXPECTED_QUARTIC_SHA
    assert np.array_equal(actual, expected)
    assert np.max(np.asarray(data["speed"], dtype=float)) > 0.0


def test_quartic_visualization_smoke_exports_reference_times_without_png(tmp_path: Path):
    manifest = write_quartic_visualization_smoke(
        tmp_path,
        times=(0.25, 0.50, 0.75),
        grid_size=9,
        write_png=False,
    )
    assert manifest["times"] == [0.25, 0.50, 0.75]
    assert len(manifest["files"]) == 3
    assert all(entry["png"] is None for entry in manifest["files"])
    assert all((tmp_path / entry["npz"]).is_file() for entry in manifest["files"])
    assert manifest["files"][0]["max_speed"] > 0.0


@pytest.mark.parametrize(
    "kwargs",
    [
        {"grid_size": 8},
        {"grid_size": 7},
        {"extent": 0.0},
        {"extent": float("nan")},
        {"times": ()},
        {"times": (0.25, 0.25)},
        {"times": (0.24,)},
        {"times": (float("nan"),)},
    ],
)
def test_quartic_visualization_smoke_rejects_invalid_contract(tmp_path: Path, kwargs):
    with pytest.raises(ValueError):
        write_quartic_visualization_smoke(tmp_path, write_png=False, **kwargs)
