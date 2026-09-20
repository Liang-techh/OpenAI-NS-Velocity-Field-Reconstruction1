"""Deterministic sampled-grid export for the published continuous ST054 velocity field.

This is a delivery adapter only.  It samples the checksum-bound continuous
``ST054Snapshot.velocity`` callable onto a Cartesian grid for convenient Python /
MATLAB visualization.  It does not define a new continuous field, does not add
scientific freedom, and does not upgrade the PDE or visual-correspondence truth
boundary.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.io import loadmat, savemat

from .st054_snapshot import (
    ST054_REFERENCE_ATOL,
    ST054_SOURCE_COMMIT,
    ST054_TESTED_MAT_SHA256,
    available_st054_models,
    load_st054,
)

DEFAULT_GRID_SIZE = 33
DEFAULT_TIMES = (0.25, 0.375, 0.5, 0.625, 0.75)
_ARRAY_ORDER = ("x", "y", "z", "t", "u", "v", "w")
_SCHEMA = "st054_sampled_velocity_grid_v1"


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _canonical_payload_sha256(arrays: dict[str, np.ndarray]) -> str:
    """Hash numerical identity independent of NPZ/MAT container metadata."""
    h = hashlib.sha256()
    h.update((_SCHEMA + "\n").encode())
    for name in _ARRAY_ORDER:
        array = np.ascontiguousarray(np.asarray(arrays[name]))
        h.update(name.encode() + b"\0")
        h.update(array.dtype.str.encode() + b"\0")
        h.update(",".join(str(v) for v in array.shape).encode() + b"\0")
        h.update(array.tobytes(order="C"))
    return h.hexdigest()


def _validate_grid_inputs(grid_size: int, times: Iterable[float]) -> np.ndarray:
    if int(grid_size) != grid_size or int(grid_size) < 3:
        raise ValueError("grid_size must be an integer >= 3")
    t = np.asarray(tuple(times), dtype=float)
    if t.ndim != 1 or t.size == 0 or not np.isfinite(t).all():
        raise ValueError("times must be a nonempty finite 1-D sequence")
    if np.any((t < 0.25) | (t > 0.75)):
        raise ValueError("all times must lie in [0.25, 0.75]")
    if np.any(np.diff(t) <= 0.0):
        raise ValueError("times must be strictly increasing")
    return t


def sample_st054_grid(
    model_id: str = "ST054-Q2",
    *,
    grid_size: int = DEFAULT_GRID_SIZE,
    times: Iterable[float] = DEFAULT_TIMES,
) -> dict[str, np.ndarray]:
    """Sample the published continuous field on ``[-2,2]^3``.

    Velocity arrays use the fixed layout ``[time, x, y, z]``.
    """
    t = _validate_grid_inputs(grid_size, times)
    field = load_st054(model_id)
    axis = np.linspace(-2.0, 2.0, int(grid_size), dtype=np.float64)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    velocity = np.empty((t.size, axis.size, axis.size, axis.size, 3), dtype=np.float64)
    for index, time in enumerate(t):
        velocity[index] = field.velocity(xx, yy, zz, float(time))
    return {
        "x": axis.copy(),
        "y": axis.copy(),
        "z": axis.copy(),
        "t": t.astype(np.float64, copy=True),
        "u": velocity[..., 0],
        "v": velocity[..., 1],
        "w": velocity[..., 2],
    }


def _mat_payload(arrays: dict[str, np.ndarray], manifest: dict[str, Any]) -> dict[str, Any]:
    payload: dict[str, Any] = {name: arrays[name] for name in _ARRAY_ORDER}
    payload.update(
        {
            "schema": _SCHEMA,
            "model_id": manifest["model_id"],
            "source_commit": ST054_SOURCE_COMMIT,
            "source_mat_sha256": ST054_TESTED_MAT_SHA256,
            "nu": 0.01,
            "support_radius": 2.0,
            "support_z": 2.0,
            "pde_validated": np.uint8(0),
            "visual_correspondence_verified": np.uint8(0),
        }
    )
    return payload


def export_st054_grid(
    output_dir: str | Path,
    model_id: str = "ST054-Q2",
    *,
    grid_size: int = DEFAULT_GRID_SIZE,
    times: Iterable[float] = DEFAULT_TIMES,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write exact-value NPZ and MATLAB-v5 copies plus a fail-closed manifest."""
    if model_id not in available_st054_models():
        raise ValueError(f"unknown model_id {model_id!r}")
    arrays = sample_st054_grid(model_id, grid_size=grid_size, times=times)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = model_id.lower().replace("-", "_") + "_velocity_grid"
    npz_path = out / f"{stem}.npz"
    mat_path = out / f"{stem}.mat"
    manifest_path = out / f"{stem}.json"
    for path in (npz_path, mat_path, manifest_path):
        if path.exists() and not overwrite:
            raise FileExistsError(f"refusing to overwrite existing export: {path}")

    payload_sha = _canonical_payload_sha256(arrays)
    manifest: dict[str, Any] = {
        "schema": _SCHEMA,
        "model_id": model_id,
        "source_repo": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
        "source_commit": ST054_SOURCE_COMMIT,
        "source_mat_sha256": ST054_TESTED_MAT_SHA256,
        "continuous_source": "openai_ns_reconstruction.load_st054(...).velocity(x,y,z,t)",
        "sampling_only": True,
        "grid_interpolation_used": False,
        "domain": {"physical": "R^3", "evaluation_cube": [-2.0, 2.0]},
        "nu": 0.01,
        "support": {"radius": 2.0, "z_abs": 2.0},
        "layout": {"u": ["time", "x", "y", "z"], "v": ["time", "x", "y", "z"], "w": ["time", "x", "y", "z"]},
        "grid_size": int(np.asarray(arrays["x"]).size),
        "times": np.asarray(arrays["t"]).tolist(),
        "canonical_payload_sha256": payload_sha,
        "engineering_callable_atol": ST054_REFERENCE_ATOL,
        "truth_boundary": {
            "pde_validated": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "source_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "external_migration": {
            "classification": "direct migration / public API only",
            "source_repo": "scipy/scipy",
            "source_commit": "2c2998c8d9d5fc63d3e6e3d330db46ee476379b2",
            "license": "BSD-3-Clause",
            "scope": "scipy.io.savemat/loadmat container I/O only; no SciPy implementation code copied",
        },
    }

    np.savez_compressed(npz_path, **arrays)
    savemat(mat_path, _mat_payload(arrays, manifest), do_compression=True)
    manifest["npz_sha256"] = _sha256_file(npz_path)
    manifest["mat_sha256"] = _sha256_file(mat_path)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    receipt = verify_st054_grid_export(npz_path, mat_path, manifest_path)
    if not receipt["passed"]:
        raise RuntimeError("ST054 grid export verification failed")
    return receipt


def verify_st054_grid_export(
    npz_path: str | Path,
    mat_path: str | Path,
    manifest_path: str | Path,
) -> dict[str, Any]:
    """Verify container hashes, exact NPZ/MAT arrays, and callable node parity."""
    npz_path = Path(npz_path)
    mat_path = Path(mat_path)
    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != _SCHEMA:
        raise ValueError("unexpected ST054 grid export schema")
    if manifest.get("source_commit") != ST054_SOURCE_COMMIT:
        raise ValueError("ST054 source commit drift")
    if manifest.get("source_mat_sha256") != ST054_TESTED_MAT_SHA256:
        raise ValueError("ST054 source MAT identity drift")
    truth = manifest.get("truth_boundary", {})
    required_false = (
        "pde_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    )
    if any(truth.get(key) is not False for key in required_false):
        raise ValueError("ST054 sampled export truth boundary was unexpectedly promoted")
    if _sha256_file(npz_path) != manifest.get("npz_sha256"):
        raise ValueError("NPZ byte checksum mismatch")
    if _sha256_file(mat_path) != manifest.get("mat_sha256"):
        raise ValueError("MAT byte checksum mismatch")

    with np.load(npz_path, allow_pickle=False) as npz:
        npz_arrays = {name: np.asarray(npz[name]) for name in _ARRAY_ORDER}
    mat = loadmat(mat_path, simplify_cells=True)
    mat_arrays = {name: np.asarray(mat[name]).squeeze() if name in ("x", "y", "z", "t") else np.asarray(mat[name]) for name in _ARRAY_ORDER}
    for name in _ARRAY_ORDER:
        if npz_arrays[name].shape != mat_arrays[name].shape or not np.array_equal(npz_arrays[name], mat_arrays[name]):
            raise ValueError(f"NPZ/MAT numerical mismatch for {name}")
    payload_sha = _canonical_payload_sha256(npz_arrays)
    if payload_sha != manifest.get("canonical_payload_sha256"):
        raise ValueError("canonical numerical payload checksum mismatch")

    x, y, z, t = (npz_arrays[name] for name in ("x", "y", "z", "t"))
    u, v, w = (npz_arrays[name] for name in ("u", "v", "w"))
    expected_shape = (t.size, x.size, y.size, z.size)
    if any(array.shape != expected_shape for array in (u, v, w)):
        raise ValueError("velocity layout mismatch")
    if not all(np.isfinite(array).all() for array in npz_arrays.values()):
        raise ValueError("non-finite value in sampled export")

    field = load_st054(str(manifest["model_id"]))
    indices = sorted(set((0, x.size // 2, x.size - 1)))
    max_error = 0.0
    for it, time in enumerate(t):
        ii, jj, kk = np.meshgrid(indices, indices, indices, indexing="ij")
        xi, yi, zi = x[ii], y[jj], z[kk]
        expected = field.velocity(xi, yi, zi, float(time))
        saved = np.stack((u[it, ii, jj, kk], v[it, ii, jj, kk], w[it, ii, jj, kk]), axis=-1)
        max_error = max(max_error, float(np.max(np.abs(expected - saved))))
    passed = bool(max_error <= ST054_REFERENCE_ATOL)
    return {
        "schema": _SCHEMA + "_receipt",
        "passed": passed,
        "model_id": manifest["model_id"],
        "canonical_payload_sha256": payload_sha,
        "npz_sha256": manifest["npz_sha256"],
        "mat_sha256": manifest["mat_sha256"],
        "callable_node_max_component_error": max_error,
        "engineering_callable_atol": ST054_REFERENCE_ATOL,
        "pde_validated": False,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
    }


def _main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--model", choices=available_st054_models(), default="ST054-Q2")
    parser.add_argument("--grid-size", type=int, default=DEFAULT_GRID_SIZE)
    parser.add_argument("--times", type=float, nargs="+", default=list(DEFAULT_TIMES))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    receipt = export_st054_grid(
        args.output_dir,
        args.model,
        grid_size=args.grid_size,
        times=args.times,
        overwrite=args.overwrite,
    )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    _main()
