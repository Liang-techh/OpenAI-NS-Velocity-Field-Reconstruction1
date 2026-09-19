"""Checksum-bound NPZ/MAT grid export for the frozen ST052 temporal child.

This is a delivery adapter, not a scientific validator.  It samples an already
loaded whole candidate on a frozen Cartesian visualization grid and emits both
Python ``.npz`` and MATLAB v5 ``.mat`` files plus a manifest.  The manifest
binds the exported numerical payload to the whole-candidate identity while
keeping PDE, visual-correspondence, and dependency-runtime truth boundaries
explicitly fail-closed.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path
import sys
from typing import Any, Iterable

import numpy as np
import scipy
from scipy.io import loadmat, savemat

SCHEMA = "st052-linear-temporal-grid-export/v1"
TASK_ID = "CR-A9-071"
DEFAULT_TIMES = (0.25, 0.375, 0.5, 0.625, 0.75)
DEFAULT_GRID_POINTS = 33
DEFAULT_DOMAIN = (-2.0, 2.0)
NPZ_FILENAME = "st052_velocity_grid.npz"
MAT_FILENAME = "st052_velocity_grid.mat"
MANIFEST_FILENAME = "st052_velocity_grid_manifest.json"
_ARRAY_NAMES = ("x", "y", "z", "t", "u", "v", "w")

# Engineering replay tolerance only.  This reuses the already-frozen 5e-12
# exact-source whole-candidate parity tolerance from the ST052 capsule lane.
# It is not a PDE, visualization, or source-correspondence acceptance gate.
# NPZ <-> MAT serialized numerical payload equality remains exact below.
CALLABLE_GRID_PARITY_ATOL = 5e-12


def callable_grid_parity_passes(max_component_error: float) -> bool:
    """Return whether scalar-vs-batched callable replay is within the fixed software tolerance."""
    error = float(max_component_error)
    return bool(np.isfinite(error) and 0.0 <= error <= CALLABLE_GRID_PARITY_ATOL)


def _canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False
    ).encode("utf-8")


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _f64(array: Any) -> np.ndarray:
    out = np.asarray(array, dtype="<f8")
    if not np.all(np.isfinite(out)):
        raise ValueError("grid export contains non-finite values")
    return np.ascontiguousarray(out)


def _payload_sha256(arrays: dict[str, np.ndarray], *, whole_candidate_identity: str) -> str:
    ordered: dict[str, np.ndarray] = {name: _f64(arrays[name]) for name in _ARRAY_NAMES}
    descriptor = {
        "schema": SCHEMA,
        "whole_candidate_identity_sha256": whole_candidate_identity,
        "array_order": list(_ARRAY_NAMES),
        "arrays": {
            name: {"dtype": "<f8", "shape": list(ordered[name].shape)}
            for name in _ARRAY_NAMES
        },
    }
    h = hashlib.sha256()
    h.update(_canonical_json_bytes(descriptor))
    for name in _ARRAY_NAMES:
        encoded = name.encode("ascii")
        h.update(len(encoded).to_bytes(2, "big"))
        h.update(encoded)
        raw = ordered[name].tobytes(order="C")
        h.update(len(raw).to_bytes(8, "big"))
        h.update(raw)
    return h.hexdigest()


def _runtime_receipt() -> dict[str, Any]:
    return {
        "python_implementation": platform.python_implementation(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "dependency_runtime_identity_closed": False,
        "binary_or_blas_identity_bound": False,
    }


def _candidate_truth(candidate: Any) -> dict[str, Any]:
    manifest = getattr(candidate, "manifest", None)
    if not isinstance(manifest, dict):
        raise ValueError("whole candidate manifest is required for export")
    truth = manifest.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("whole candidate truth boundary is missing")
    if truth.get("whole_child_bundle_materialized") is not True:
        raise ValueError("export requires a materialized whole-child bundle")
    for key in (
        "pde_validated",
        "visualization_ready",
        "visual_correspondence_verified",
        "production_candidate_selected",
        "source_correspondence_verified",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"whole candidate must keep {key}=false at this export stage")
    return truth


def sample_velocity_grid(
    candidate: Any,
    *,
    grid_points: int = DEFAULT_GRID_POINTS,
    times: Iterable[float] = DEFAULT_TIMES,
    domain: tuple[float, float] = DEFAULT_DOMAIN,
) -> dict[str, np.ndarray]:
    """Sample ``candidate.velocity`` on a Cartesian ``ij`` grid.

    Returned velocity components use shape ``(nt, nx, ny, nz)``.  Coordinates
    and times are float64 one-dimensional arrays.
    """
    _candidate_truth(candidate)
    n = int(grid_points)
    if n < 3:
        raise ValueError("grid_points must be >= 3")
    lo, hi = map(float, domain)
    if not np.isfinite(lo) or not np.isfinite(hi) or not lo < hi:
        raise ValueError("domain must be finite with lower < upper")
    tt = _f64(tuple(float(v) for v in times))
    if tt.ndim != 1 or len(tt) == 0 or np.any(np.diff(tt) <= 0):
        raise ValueError("times must be a nonempty strictly increasing sequence")
    if tt[0] < 0.25 or tt[-1] > 0.75:
        raise ValueError("ST052 export times must stay inside [0.25, 0.75]")

    x = _f64(np.linspace(lo, hi, n, dtype=float))
    y = x.copy()
    z = x.copy()
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    u = np.empty((len(tt), n, n, n), dtype="<f8")
    v = np.empty_like(u)
    w = np.empty_like(u)
    for i, time in enumerate(tt):
        velocity = _f64(candidate.velocity(xx, yy, zz, float(time)))
        if velocity.shape != xx.shape + (3,):
            raise ValueError(
                f"candidate velocity returned shape {velocity.shape}, expected {xx.shape + (3,)}"
            )
        u[i] = velocity[..., 0]
        v[i] = velocity[..., 1]
        w[i] = velocity[..., 2]
    return {"x": x, "y": y, "z": z, "t": tt, "u": u, "v": v, "w": w}


def export_velocity_grid(
    candidate: Any,
    out_dir: str | Path,
    *,
    grid_points: int = DEFAULT_GRID_POINTS,
    times: Iterable[float] = DEFAULT_TIMES,
    domain: tuple[float, float] = DEFAULT_DOMAIN,
) -> dict[str, Any]:
    """Write NPZ + MATLAB v5 exports and a checksum/truth-boundary manifest."""
    truth = _candidate_truth(candidate)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name in (NPZ_FILENAME, MAT_FILENAME, MANIFEST_FILENAME):
        if (out_dir / name).exists():
            raise ValueError(f"refusing to overwrite existing export file: {name}")

    whole_id = str(getattr(candidate, "identity_sha256", ""))
    candidate_id = str(getattr(candidate, "candidate_id", ""))
    if len(whole_id) != 64 or not candidate_id:
        raise ValueError("candidate ID and 64-hex whole-candidate identity are required")

    arrays = sample_velocity_grid(
        candidate, grid_points=grid_points, times=times, domain=domain
    )
    payload_sha = _payload_sha256(arrays, whole_candidate_identity=whole_id)

    npz_path = out_dir / NPZ_FILENAME
    mat_path = out_dir / MAT_FILENAME
    np.savez_compressed(npz_path, **arrays)
    savemat(
        mat_path,
        {
            **arrays,
            "whole_candidate_identity_sha256": whole_id,
            "grid_payload_sha256": payload_sha,
            "array_layout": "u,v,w: [time,x,y,z]; meshgrid indexing=ij",
        },
        do_compression=True,
        oned_as="row",
    )

    runtime = _runtime_receipt()
    manifest = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "candidate_id": candidate_id,
        "whole_candidate_identity_sha256": whole_id,
        "grid_payload_sha256": payload_sha,
        "grid": {
            "domain": [float(arrays["x"][0]), float(arrays["x"][-1])],
            "shape_xyz": [len(arrays["x"]), len(arrays["y"]), len(arrays["z"])],
            "times": arrays["t"].tolist(),
            "velocity_component_shape": list(arrays["u"].shape),
            "array_layout": "u,v,w: [time,x,y,z]",
            "meshgrid_indexing": "ij",
            "dtype": "float64 little-endian canonical payload",
        },
        "files": {
            NPZ_FILENAME: {"sha256": _sha256_file(npz_path), "format": "NumPy compressed NPZ"},
            MAT_FILENAME: {"sha256": _sha256_file(mat_path), "format": "MATLAB v5 MAT via scipy.io.savemat"},
        },
        "external_method": {
            "source_repo": "scipy/scipy",
            "api": "scipy.io.savemat/loadmat",
            "license": "BSD-3-Clause",
            "classification": "direct migration / public API only",
            "source_commit_bound": False,
            "resolved_package_version": runtime["scipy_version"],
            "scope": "MAT-file serialization/deserialization only; no upstream implementation copied",
        },
        "execution_environment": runtime,
        "matlab_compatibility": {
            "mat_v5_written": True,
            "scipy_loadmat_roundtrip_required": True,
            "actual_matlab_runtime_executed": False,
        },
        "truth_boundary": {
            "whole_child_bundle_materialized": truth["whole_child_bundle_materialized"],
            "grid_export_materialized": True,
            "grid_export_payload_checksum_bound_to_whole_candidate": True,
            "dependency_runtime_identity_closed": False,
            "standalone_package_parent_runtime_ready": bool(
                truth.get("standalone_package_parent_runtime_ready", False)
            ),
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "production_candidate_selected": False,
            "held_out_temporal_child_pde_residual_evaluated": bool(
                truth.get("held_out_temporal_child_pde_residual_evaluated", False)
            ),
            "pde_validated": False,
            "source_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    (out_dir / MANIFEST_FILENAME).write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return manifest


def verify_velocity_grid_export(out_dir: str | Path) -> dict[str, Any]:
    """Verify file bytes and exact numerical NPZ/MAT payload round-trip."""
    out_dir = Path(out_dir)
    manifest = json.loads((out_dir / MANIFEST_FILENAME).read_text(encoding="utf-8"))
    if manifest.get("schema") != SCHEMA or manifest.get("task_id") != TASK_ID:
        raise ValueError("unsupported ST052 grid-export manifest")
    whole_id = str(manifest.get("whole_candidate_identity_sha256", ""))
    files = manifest.get("files")
    if not isinstance(files, dict):
        raise ValueError("grid-export file receipt missing")
    for name in (NPZ_FILENAME, MAT_FILENAME):
        path = out_dir / name
        record = files.get(name)
        if not path.is_file() or not isinstance(record, dict):
            raise ValueError(f"missing export file: {name}")
        if record.get("sha256") != _sha256_file(path):
            raise ValueError(f"export file checksum mismatch: {name}")

    with np.load(out_dir / NPZ_FILENAME, allow_pickle=False) as data:
        npz_arrays = {name: _f64(data[name]) for name in _ARRAY_NAMES}
    raw_mat = loadmat(out_dir / MAT_FILENAME, squeeze_me=False, struct_as_record=False)
    mat_arrays = {
        "x": _f64(raw_mat["x"].reshape(-1)),
        "y": _f64(raw_mat["y"].reshape(-1)),
        "z": _f64(raw_mat["z"].reshape(-1)),
        "t": _f64(raw_mat["t"].reshape(-1)),
        "u": _f64(raw_mat["u"]),
        "v": _f64(raw_mat["v"]),
        "w": _f64(raw_mat["w"]),
    }
    npz_payload = _payload_sha256(npz_arrays, whole_candidate_identity=whole_id)
    mat_payload = _payload_sha256(mat_arrays, whole_candidate_identity=whole_id)
    expected = manifest.get("grid_payload_sha256")
    if npz_payload != expected or mat_payload != expected:
        raise ValueError("NPZ/MAT canonical grid payload checksum mismatch")
    for name in _ARRAY_NAMES:
        if not np.array_equal(npz_arrays[name], mat_arrays[name]):
            raise ValueError(f"NPZ/MAT numerical payload differs: {name}")

    truth = manifest.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("grid-export truth boundary missing")
    for key in (
        "dependency_runtime_identity_closed",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "openai_field_identified",
    ):
        if truth.get(key) is not False:
            raise ValueError(f"grid-export truth boundary must keep {key}=false")
    return {
        "manifest": manifest,
        "npz_payload_sha256": npz_payload,
        "mat_payload_sha256": mat_payload,
        "npz_mat_exact_equal": True,
    }