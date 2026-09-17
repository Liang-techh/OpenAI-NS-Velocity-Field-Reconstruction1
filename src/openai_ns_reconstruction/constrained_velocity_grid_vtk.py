from __future__ import annotations

from collections.abc import Mapping
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

import numpy as np


_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_FALSE_CLAIMS = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def _axis(name: str, values: Any, minimum: int) -> np.ndarray:
    arr = np.asarray(values, dtype=np.float64)
    if arr.ndim != 1 or arr.size < minimum:
        raise ValueError(f"{name} must be a 1-D axis with at least {minimum} entries")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    if arr.size > 1 and not np.all(np.diff(arr) > 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    return np.ascontiguousarray(arr)


def _identity(candidate_sha256: str, provenance: Mapping[str, Any]) -> tuple[str, dict[str, Any]]:
    candidate_sha = str(candidate_sha256).strip().lower()
    if _SHA256_RE.fullmatch(candidate_sha) is None:
        raise ValueError("candidate_sha256 must be 64 lowercase hexadecimal characters")
    if not isinstance(provenance, Mapping) or not provenance:
        raise ValueError("provenance must be a non-empty mapping")
    try:
        encoded = json.dumps(provenance, sort_keys=True, separators=(",", ":"), allow_nan=False)
        provenance_json = json.loads(encoded)
    except (TypeError, ValueError) as exc:
        raise ValueError("provenance must be finite JSON data") from exc
    return candidate_sha, provenance_json


def _file_hash(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _grid_hash(arrays: list[np.ndarray], candidate_sha: str, provenance: dict[str, Any]) -> str:
    digest = sha256(b"constrained-vtk-grid-v1")
    for arr in arrays:
        item = np.ascontiguousarray(arr, dtype="<f8")
        digest.update(np.asarray(item.shape, dtype="<i8").tobytes())
        digest.update(item.tobytes())
    digest.update(candidate_sha.encode("ascii"))
    digest.update(json.dumps(provenance, sort_keys=True, separators=(",", ":")).encode())
    return digest.hexdigest()


def _numbers(values: np.ndarray, width: int = 6) -> list[str]:
    flat = np.asarray(values, dtype=np.float64).reshape(-1)
    return [" ".join(f"{float(v):.17g}" for v in flat[i : i + width]) for i in range(0, flat.size, width)]


def _vtk_text(x: np.ndarray, y: np.ndarray, z: np.ndarray, velocity: np.ndarray, time: float) -> str:
    nx, ny, nz = x.size, y.size, z.size
    vectors = np.transpose(velocity, (2, 1, 0, 3)).reshape(-1, 3)
    speed = np.linalg.norm(vectors, axis=1)
    lines = [
        "# vtk DataFile Version 3.0",
        f"truth-bounded velocity t={time:.17g}",
        "ASCII",
        "DATASET RECTILINEAR_GRID",
        f"DIMENSIONS {nx} {ny} {nz}",
        f"X_COORDINATES {nx} double",
        *_numbers(x),
        f"Y_COORDINATES {ny} double",
        *_numbers(y),
        f"Z_COORDINATES {nz} double",
        *_numbers(z),
        f"POINT_DATA {nx * ny * nz}",
        "VECTORS velocity double",
    ]
    lines.extend(" ".join(f"{float(v):.17g}" for v in row) for row in vectors)
    lines.extend(["SCALARS speed double 1", "LOOKUP_TABLE default", *_numbers(speed)])
    return "\n".join(lines) + "\n"


def export_velocity_grid_vtk(
    output_directory: str | Path,
    *,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    velocity_txyz: Any,
    candidate_sha256: str,
    provenance: Mapping[str, Any],
) -> dict[str, Any]:
    """Export frozen ``(t,x,y,z,[u,v,w])`` samples to legacy VTK rectilinear grids."""

    x_arr = _axis("x", x, 2)
    y_arr = _axis("y", y, 2)
    z_arr = _axis("z", z, 2)
    t_arr = _axis("t", t, 1)
    candidate_sha, provenance_json = _identity(candidate_sha256, provenance)
    velocity = np.asarray(velocity_txyz, dtype=np.float64)
    expected = (t_arr.size, x_arr.size, y_arr.size, z_arr.size, 3)
    if velocity.shape != expected:
        raise ValueError(f"velocity_txyz must have shape {expected}")
    if not np.all(np.isfinite(velocity)):
        raise ValueError("velocity_txyz must be finite")
    velocity = np.ascontiguousarray(velocity)
    speeds = np.linalg.norm(velocity.reshape(t_arr.size, -1, 3), axis=2)
    if np.any(np.max(speeds, axis=1) == 0.0):
        raise ValueError("exact-zero time slices are rejected; the registered energy gate remains separate")

    output = Path(output_directory)
    if output.exists():
        raise FileExistsError("refusing to overwrite an existing VTK export directory")
    output.mkdir(parents=True)

    files = []
    for index, time_value in enumerate(t_arr):
        filename = f"velocity_t{index:03d}.vtk"
        path = output / filename
        path.write_text(_vtk_text(x_arr, y_arr, z_arr, velocity[index], float(time_value)), encoding="ascii")
        files.append(
            {
                "index": index,
                "time": float(time_value),
                "filename": filename,
                "sha256": _file_hash(path),
                "max_speed": float(np.max(speeds[index])),
                "rms_speed": float(np.sqrt(np.mean(speeds[index] ** 2))),
            }
        )

    truth = {key: False for key in _FALSE_CLAIMS}
    truth.update(
        {
            "velocity_export_ready": True,
            "registered_nontriviality_gate_assessed": False,
            "vtk_derivatives_valid_for_pde_acceptance": False,
        }
    )
    manifest = {
        "schema": "constrained_velocity_grid_vtk_manifest_v1",
        "format": "VTK legacy ASCII",
        "dataset": "RECTILINEAR_GRID",
        "component_order": ["u", "v", "w"],
        "input_layout": "time,x,y,z,component",
        "vtk_point_order": "x-fastest, then y, then z",
        "axes": {"x": x_arr.tolist(), "y": y_arr.tolist(), "z": z_arr.tolist(), "t": t_arr.tolist()},
        "candidate_sha256": candidate_sha,
        "grid_sha256": _grid_hash([x_arr, y_arr, z_arr, t_arr, velocity], candidate_sha, provenance_json),
        "provenance": provenance_json,
        "files": files,
        "truth_boundary": truth,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def load_velocity_grid_vtk_manifest(output_directory: str | Path) -> dict[str, Any]:
    """Load the export manifest and fail closed on file or truth-state tampering."""

    output = Path(output_directory)
    try:
        manifest = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("missing or invalid VTK manifest") from exc
    if manifest.get("schema") != "constrained_velocity_grid_vtk_manifest_v1":
        raise ValueError("unexpected VTK manifest schema")
    if manifest.get("component_order") != ["u", "v", "w"] or manifest.get("vtk_point_order") != "x-fastest, then y, then z":
        raise ValueError("VTK ordering metadata was altered")
    _identity(manifest.get("candidate_sha256", ""), manifest.get("provenance", {}))
    truth = manifest.get("truth_boundary", {})
    for key in _FALSE_CLAIMS:
        if truth.get(key) is not False:
            raise ValueError(f"truth flag {key} was altered")
    if truth.get("registered_nontriviality_gate_assessed") is not False or truth.get("vtk_derivatives_valid_for_pde_acceptance") is not False:
        raise ValueError("VTK export cannot promote energy or PDE evidence")
    files = manifest.get("files")
    times = manifest.get("axes", {}).get("t", [])
    if not isinstance(files, list) or len(files) != len(times):
        raise ValueError("VTK file list does not match time axis")
    for index, item in enumerate(files):
        filename = f"velocity_t{index:03d}.vtk"
        path = output / filename
        if item.get("filename") != filename or not path.is_file() or _file_hash(path) != item.get("sha256"):
            raise ValueError(f"VTK file integrity check failed: {filename}")
    return manifest
