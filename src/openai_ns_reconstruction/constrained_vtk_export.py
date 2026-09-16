"""Truth-bounded VTK export for frozen 3D velocity callables.

This module is visualization/serialization infrastructure only. Export success
does not establish Navier--Stokes validity, divergence-free structure,
correspondence to OpenAI's hidden numerical field, or blow-up.

The writer deliberately implements only the small legacy ASCII VTK
RECTILINEAR_GRID subset needed for a vector velocity field plus speed. It does
not copy meshio/VTK implementation code.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

VelocityCallable = Callable[[np.ndarray, float], np.ndarray]

_SCHEMA_VERSION = 1
_CLAIM_SCOPE = "visualization_export_only"
_AXIS_ORDER = "x,y,z with x-fastest VTK point ordering"
_VTK_DATASET = "RECTILINEAR_GRID"


@dataclass(frozen=True)
class VTKSeriesExport:
    output_dir: Path
    manifest_path: Path
    vtk_paths: tuple[Path, ...]


def _strict_axis(name: str, values: Iterable[float]) -> np.ndarray:
    axis = np.asarray(tuple(values), dtype=float)
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError(f"{name} must be a 1D axis with at least two points")
    if not np.all(np.isfinite(axis)):
        raise ValueError(f"{name} must contain only finite values")
    if not np.all(np.diff(axis) > 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    return axis


def _strict_times(values: Iterable[float]) -> np.ndarray:
    times = np.asarray(tuple(values), dtype=float)
    if times.ndim != 1 or times.size < 1:
        raise ValueError("times must be a nonempty 1D sequence")
    if not np.all(np.isfinite(times)):
        raise ValueError("times must contain only finite values")
    if times.size > 1 and not np.all(np.diff(times) > 0.0):
        raise ValueError("times must be strictly increasing")
    return times


def _sample_velocity(
    velocity: VelocityCallable,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    time: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    points = np.column_stack(
        (X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F"))
    )
    values = np.asarray(velocity(points, float(time)), dtype=float)
    if values.shape != (points.shape[0], 3):
        raise ValueError(
            "velocity(points, time) must return shape "
            f"({points.shape[0]}, 3), got {values.shape}"
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("velocity returned non-finite values")

    shape = (x.size, y.size, z.size)
    U = values[:, 0].reshape(shape, order="F")
    V = values[:, 1].reshape(shape, order="F")
    W = values[:, 2].reshape(shape, order="F")
    return U, V, W


def _write_numbers(handle, values: np.ndarray, per_line: int = 6) -> None:
    flat = np.asarray(values, dtype=float).ravel()
    for start in range(0, flat.size, per_line):
        chunk = flat[start : start + per_line]
        handle.write(" ".join(f"{float(v):.17g}" for v in chunk))
        handle.write("\n")


def write_rectilinear_velocity_vtk(
    path: str | Path,
    *,
    x: Iterable[float],
    y: Iterable[float],
    z: Iterable[float],
    U: np.ndarray,
    V: np.ndarray,
    W: np.ndarray,
    title: str = "constrained velocity visualization export",
) -> Path:
    """Write one legacy ASCII VTK RECTILINEAR_GRID velocity snapshot.

    Array convention: U[i,j,k] == u(x[i], y[j], z[k]) (same for V and W).
    VTK point data are emitted with x varying fastest.
    """

    x_axis = _strict_axis("x", x)
    y_axis = _strict_axis("y", y)
    z_axis = _strict_axis("z", z)
    shape = (x_axis.size, y_axis.size, z_axis.size)

    components = []
    for name, value in (("U", U), ("V", V), ("W", W)):
        array = np.asarray(value, dtype=float)
        if array.shape != shape:
            raise ValueError(f"{name} must have shape {shape}, got {array.shape}")
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{name} must contain only finite values")
        components.append(array)

    U_arr, V_arr, W_arr = components
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    safe_title = " ".join(str(title).split())[:255]

    u_flat = U_arr.ravel(order="F")
    v_flat = V_arr.ravel(order="F")
    w_flat = W_arr.ravel(order="F")
    speed = np.sqrt(u_flat * u_flat + v_flat * v_flat + w_flat * w_flat)

    with target.open("w", encoding="ascii", newline="\n") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write(f"{safe_title}\n")
        handle.write("ASCII\n")
        handle.write(f"DATASET {_VTK_DATASET}\n")
        handle.write(f"DIMENSIONS {x_axis.size} {y_axis.size} {z_axis.size}\n")

        handle.write(f"X_COORDINATES {x_axis.size} double\n")
        _write_numbers(handle, x_axis)
        handle.write(f"Y_COORDINATES {y_axis.size} double\n")
        _write_numbers(handle, y_axis)
        handle.write(f"Z_COORDINATES {z_axis.size} double\n")
        _write_numbers(handle, z_axis)

        n_points = int(np.prod(shape))
        handle.write(f"POINT_DATA {n_points}\n")
        handle.write("VECTORS velocity double\n")
        for ux, uy, uz in zip(u_flat, v_flat, w_flat, strict=True):
            handle.write(f"{ux:.17g} {uy:.17g} {uz:.17g}\n")

        handle.write("SCALARS speed double 1\n")
        handle.write("LOOKUP_TABLE default\n")
        _write_numbers(handle, speed)

    return target


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def export_velocity_vtk_series(
    velocity: VelocityCallable,
    *,
    output_dir: str | Path,
    x: Iterable[float],
    y: Iterable[float],
    z: Iterable[float],
    times: Iterable[float],
    candidate_id: str,
) -> VTKSeriesExport:
    """Sample a frozen callable and export a ParaView-readable VTK series.

    The manifest is explicitly truth-bounded. ``visualization_export_ready``
    means only that deterministic files were produced in the declared schema.
    """

    if not str(candidate_id).strip():
        raise ValueError("candidate_id must be nonempty")

    x_axis = _strict_axis("x", x)
    y_axis = _strict_axis("y", y)
    z_axis = _strict_axis("z", z)
    time_axis = _strict_times(times)

    directory = Path(output_dir)
    directory.mkdir(parents=True, exist_ok=True)

    records: list[dict[str, object]] = []
    paths: list[Path] = []

    for index, time in enumerate(time_axis):
        U, V, W = _sample_velocity(velocity, x_axis, y_axis, z_axis, float(time))
        filename = f"velocity_t{index:03d}.vtk"
        vtk_path = write_rectilinear_velocity_vtk(
            directory / filename,
            x=x_axis,
            y=y_axis,
            z=z_axis,
            U=U,
            V=V,
            W=W,
            title=f"{candidate_id} t={float(time):.17g} {_CLAIM_SCOPE}",
        )
        paths.append(vtk_path)
        records.append(
            {
                "index": index,
                "time": float(time),
                "file": filename,
                "sha256": _sha256(vtk_path),
            }
        )

    manifest = {
        "schema_version": _SCHEMA_VERSION,
        "format": "legacy_vtk_ascii_rectilinear_grid",
        "dataset": _VTK_DATASET,
        "candidate_id": str(candidate_id),
        "claim_scope": _CLAIM_SCOPE,
        "axis_order": _AXIS_ORDER,
        "grid_shape": [int(x_axis.size), int(y_axis.size), int(z_axis.size)],
        "x": x_axis.tolist(),
        "y": y_axis.tolist(),
        "z": z_axis.tolist(),
        "times": time_axis.tolist(),
        "files": records,
        "claims": {
            "visualization_export_ready": True,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "truth_boundary": (
            "VTK export success verifies serialization/schema only. It does not "
            "validate Navier-Stokes, divergence, support, visual correspondence "
            "to OpenAI public media, hidden-field identity, or blow-up."
        ),
    }

    manifest_path = directory / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return VTKSeriesExport(
        output_dir=directory,
        manifest_path=manifest_path,
        vtk_paths=tuple(paths),
    )
