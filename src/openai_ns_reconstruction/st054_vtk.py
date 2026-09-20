"""Portable legacy-VTK handoff for the published continuous ST054 velocity snapshots.

This is a delivery adapter.  It samples the immutable continuous
``load_st054(...).velocity(x,y,z,t)`` callable onto a Cartesian grid and writes
an ASCII VTK legacy ``STRUCTURED_POINTS`` file that ParaView/VTK-style tools can
consume.  It does not fit or interpolate a new continuous field and it does not
perform Navier--Stokes acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from .st054_snapshot import ST054_REFERENCE_ATOL, available_st054_models, load_st054


SCHEMA = "st054_legacy_vtk_v1"
TASK_ID = "CR-A9-077"
DEFAULT_GRID_SIZE = 33
VTK_SPEC = "https://docs.vtk.org/en/latest/vtk_file_formats/vtk_legacy_file_format.html"
VTK_SOURCE = {
    "source_repo": "Kitware/VTK",
    "source_commit": "aba012833b82d219168fcb1e7066dfcba85a4118",
    "license": "BSD-3-Clause",
    "classification": "suitable_reimplementation",
    "migration_scope": "legacy VTK STRUCTURED_POINTS file-format contract only",
    "copied_implementation_code": False,
    "difference": "repository-local sampler/writer; no VTK implementation source is copied",
}
MESHIO_SOURCE = {
    "source_repo": "nschloe/meshio",
    "source_commit": "b2ee99842e119901349fdeee06b5bf61e01f450a",
    "version": "5.3.5",
    "license": "MIT",
    "classification": "direct_verification_ci_only",
    "migration_scope": "independent legacy-VTK read-back in CI only",
    "runtime_dependency": False,
    "copied_implementation_code": False,
}
_HARD_FALSE = (
    "visualization_ready",
    "visual_correspondence_verified",
    "source_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_payload_sha256(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: float,
    velocity: np.ndarray,
    speed: np.ndarray,
) -> str:
    digest = hashlib.sha256()
    for name, value in (
        ("x", np.asarray(x, dtype="<f8")),
        ("y", np.asarray(y, dtype="<f8")),
        ("z", np.asarray(z, dtype="<f8")),
        ("t", np.asarray([t], dtype="<f8")),
        ("velocity", np.asarray(velocity, dtype="<f8")),
        ("speed", np.asarray(speed, dtype="<f8")),
    ):
        array = np.ascontiguousarray(value)
        digest.update(name.encode("ascii") + b"\0")
        digest.update(str(array.shape).encode("ascii") + b"\0")
        digest.update(array.dtype.str.encode("ascii") + b"\0")
        digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def _grid_axis(grid_size: int) -> np.ndarray:
    size = int(grid_size)
    if size != grid_size or not (3 <= size <= 129):
        raise ValueError("grid_size must be an integer in [3, 129]")
    return np.linspace(-2.0, 2.0, size, dtype=float)


def sample_st054_grid(
    model_id: str,
    time: float,
    grid_size: int = DEFAULT_GRID_SIZE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Sample one published continuous ST054 field on the declared Cartesian box."""
    field = load_st054(model_id)
    t = float(time)
    if not np.isfinite(t) or not (field.tmin <= t <= field.tmax):
        raise ValueError(f"time must lie in [{field.tmin}, {field.tmax}]")
    axis = _grid_axis(grid_size)
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    velocity = np.asarray(field.velocity(xx, yy, zz, t), dtype=float)
    expected_shape = (axis.size, axis.size, axis.size, 3)
    if velocity.shape != expected_shape:
        raise RuntimeError(f"unexpected velocity shape {velocity.shape}; expected {expected_shape}")
    if not np.isfinite(velocity).all():
        raise RuntimeError("continuous ST054 evaluator produced non-finite velocity")
    speed = np.linalg.norm(velocity, axis=-1)
    return axis.copy(), axis.copy(), axis.copy(), velocity, speed


def _vtk_vector_rows(velocity: np.ndarray) -> np.ndarray:
    """Return legacy-VTK point ordering: x varies fastest, then y, then z."""
    values = np.asarray(velocity, dtype=float)
    if values.ndim != 4 or values.shape[-1] != 3:
        raise ValueError("velocity must have shape [nx, ny, nz, 3]")
    return np.ascontiguousarray(values.transpose(2, 1, 0, 3)).reshape(-1, 3)


def _vtk_scalar_rows(values: np.ndarray) -> np.ndarray:
    data = np.asarray(values, dtype=float)
    if data.ndim != 3:
        raise ValueError("scalar field must have shape [nx, ny, nz]")
    return np.ascontiguousarray(data.transpose(2, 1, 0)).reshape(-1)


def _write_legacy_structured_points(
    path: Path,
    *,
    model_id: str,
    time: float,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    velocity: np.ndarray,
    speed: np.ndarray,
) -> None:
    nx, ny, nz = len(x), len(y), len(z)
    if velocity.shape != (nx, ny, nz, 3) or speed.shape != (nx, ny, nz):
        raise ValueError("grid/value shape mismatch")
    spacing = (
        float(x[1] - x[0]),
        float(y[1] - y[0]),
        float(z[1] - z[0]),
    )
    if not (
        np.allclose(np.diff(x), spacing[0], rtol=0.0, atol=1e-15)
        and np.allclose(np.diff(y), spacing[1], rtol=0.0, atol=1e-15)
        and np.allclose(np.diff(z), spacing[2], rtol=0.0, atol=1e-15)
    ):
        raise ValueError("legacy STRUCTURED_POINTS export requires uniform axes")
    vectors = _vtk_vector_rows(velocity)
    scalars = _vtk_scalar_rows(speed)
    with path.open("x", encoding="ascii", newline="\n") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write(f"{model_id} continuous ST054 sample at t={time:.17g}; delivery-only\n")
        handle.write("ASCII\n")
        handle.write("DATASET STRUCTURED_POINTS\n")
        handle.write(f"DIMENSIONS {nx} {ny} {nz}\n")
        handle.write(f"ORIGIN {x[0]:.17g} {y[0]:.17g} {z[0]:.17g}\n")
        handle.write(f"SPACING {spacing[0]:.17g} {spacing[1]:.17g} {spacing[2]:.17g}\n")
        handle.write(f"POINT_DATA {nx * ny * nz}\n")
        handle.write("VECTORS velocity double\n")
        for ux, uy, uz in vectors:
            handle.write(f"{ux:.17g} {uy:.17g} {uz:.17g}\n")
        handle.write("SCALARS speed double 1\n")
        handle.write("LOOKUP_TABLE default\n")
        for value in scalars:
            handle.write(f"{value:.17g}\n")


def verify_vtk_manifest(manifest: dict[str, Any]) -> None:
    if manifest.get("schema") != SCHEMA or manifest.get("task_id") != TASK_ID:
        raise ValueError("unexpected VTK export schema/task")
    if manifest.get("sampling_only") is not True:
        raise ValueError("VTK handoff must remain sampling-only")
    if manifest.get("grid_interpolation_used") is not False:
        raise ValueError("grid interpolation cannot be promoted as the continuous field")
    if manifest.get("continuous_authority") != "load_st054(...).velocity(x,y,z,t)":
        raise ValueError("continuous field authority drift")
    truth = manifest.get("truth_boundary")
    if not isinstance(truth, dict):
        raise ValueError("missing truth boundary")
    for key in _HARD_FALSE:
        if truth.get(key) is not False:
            raise ValueError(f"truth-boundary promotion rejected: {key}")


def verify_vtk_export(vtk_path: str | Path, manifest_path: str | Path) -> dict[str, Any]:
    vtk_file = Path(vtk_path)
    manifest_file = Path(manifest_path)
    manifest = json.loads(manifest_file.read_text(encoding="utf-8"))
    verify_vtk_manifest(manifest)
    if _sha256(vtk_file) != manifest.get("vtk_sha256"):
        raise ValueError("VTK byte identity mismatch")

    with vtk_file.open("r", encoding="ascii") as handle:
        header = [handle.readline().rstrip("\n") for _ in range(10)]
    if header[0] != "# vtk DataFile Version 3.0":
        raise ValueError("unexpected legacy VTK version header")
    if header[2] != "ASCII" or header[3] != "DATASET STRUCTURED_POINTS":
        raise ValueError("unexpected legacy VTK dataset encoding/type")
    dims = tuple(int(item) for item in header[4].split()[1:])
    expected_dims = tuple(int(item) for item in manifest["grid"]["dimensions"])
    if dims != expected_dims:
        raise ValueError("VTK dimensions do not match manifest")
    point_count = int(header[7].split()[1])
    if point_count != int(np.prod(expected_dims)):
        raise ValueError("VTK point count does not match grid dimensions")
    if header[8] != "VECTORS velocity double":
        raise ValueError("velocity vector declaration missing")
    return manifest


def export_st054_vtk(
    output_path: str | Path,
    model_id: str,
    time: float,
    *,
    grid_size: int = DEFAULT_GRID_SIZE,
    overwrite: bool = False,
) -> dict[str, Any]:
    """Write one sampled ST054 velocity to legacy VTK plus a JSON provenance manifest."""
    if model_id not in available_st054_models():
        raise ValueError(f"unknown ST054 model: {model_id!r}")
    output = Path(output_path)
    if output.suffix.lower() != ".vtk":
        raise ValueError("output_path must end in .vtk")
    manifest_path = output.with_suffix(output.suffix + ".json")
    if (output.exists() or manifest_path.exists()) and not overwrite:
        raise FileExistsError("refusing to overwrite an existing VTK export or manifest")
    if overwrite:
        output.unlink(missing_ok=True)
        manifest_path.unlink(missing_ok=True)
    output.parent.mkdir(parents=True, exist_ok=True)

    x, y, z, velocity, speed = sample_st054_grid(model_id, time, grid_size)
    field = load_st054(model_id)
    _write_legacy_structured_points(
        output,
        model_id=model_id,
        time=float(time),
        x=x,
        y=y,
        z=z,
        velocity=velocity,
        speed=speed,
    )
    manifest: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "model_id": model_id,
        "physical_time": float(time),
        "format": {
            "family": "VTK legacy",
            "version": "3.0",
            "dataset": "STRUCTURED_POINTS",
            "encoding": "ASCII",
            "spec": VTK_SPEC,
            "point_order": "x-fastest, then y, then z",
            "point_vectors": ["velocity"],
            "point_scalars": ["speed"],
        },
        "grid": {
            "dimensions": [int(x.size), int(y.size), int(z.size)],
            "origin": [float(x[0]), float(y[0]), float(z[0])],
            "spacing": [
                float(x[1] - x[0]),
                float(y[1] - y[0]),
                float(z[1] - z[0]),
            ],
            "evaluation_box": [[-2.0, 2.0], [-2.0, 2.0], [-2.0, 2.0]],
        },
        "continuous_authority": "load_st054(...).velocity(x,y,z,t)",
        "sampling_only": True,
        "grid_interpolation_used": False,
        "engineering_node_replay_atol": ST054_REFERENCE_ATOL,
        "source": {
            "source_commit": field.source_commit,
            "source_mat_sha256": field.mat_sha256,
            "nu": float(field.nu),
            "support_radius": float(field.support_radius),
            "support_z": float(field.support_z),
            "time_interval": [float(field.tmin), float(field.tmax)],
        },
        "payload_sha256": _canonical_payload_sha256(
            x, y, z, float(time), velocity, speed
        ),
        "vtk_sha256": _sha256(output),
        "external_result_screening": [VTK_SOURCE, MESHIO_SOURCE],
        "truth_boundary": {
            "vtk_export_ready": True,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "source_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    verify_vtk_manifest(manifest)
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    verify_vtk_export(output, manifest_path)
    return manifest


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Sample a published continuous ST054 velocity into legacy VTK."
    )
    parser.add_argument("--model", choices=available_st054_models(), required=True)
    parser.add_argument("--time", type=float, required=True)
    parser.add_argument("--grid-size", type=int, default=DEFAULT_GRID_SIZE)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main() -> int:
    args = _parser().parse_args()
    manifest = export_st054_vtk(
        args.out,
        args.model,
        args.time,
        grid_size=args.grid_size,
        overwrite=args.overwrite,
    )
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
