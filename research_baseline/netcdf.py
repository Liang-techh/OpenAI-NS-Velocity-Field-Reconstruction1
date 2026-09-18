"""Truth-bounded NetCDF interchange for the retained ST006 velocity field.

This module is a delivery/visualization adapter only.  It does not change the
published velocity, fit OpenAI imagery, or provide PDE-acceptance evidence.
In particular, derivatives taken from an exported grid are not a substitute
for the repository's independent full-domain validation operator.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import uuid

import numpy as np
from scipy.io import netcdf_file

from . import CANDIDATE_SHA256, load_best

SCHEMA = "st006_velocity_netcdf_v1"
REFERENCE_TIMES = np.array([0.25, 0.50, 0.75], dtype=np.float64)
EVALUATION_BOUNDS = (-2.0, 2.0)
_FALSE_TRUTH_FLAGS = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


@dataclass(frozen=True)
class ST006NetCDFGrid:
    """Validated, immutable arrays loaded from one governed NetCDF file."""

    time: np.ndarray
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    speed: np.ndarray
    grid_sha256: str
    candidate_sha256: str = CANDIDATE_SHA256
    visualization_ready: bool = False
    visual_correspondence_verified: bool = False
    pde_validated: bool = False
    paper_exact: bool = False
    openai_field_identified: bool = False
    blowup_proved: bool = False


def _validate_resolution(resolution: int) -> int:
    if isinstance(resolution, bool) or not isinstance(resolution, (int, np.integer)):
        raise ValueError("resolution must be an odd integer")
    resolution = int(resolution)
    if resolution < 5 or resolution > 129 or resolution % 2 == 0:
        raise ValueError("resolution must be odd and lie in [5,129]")
    return resolution


def _hash_array(hasher: "hashlib._Hash", name: str, value: np.ndarray) -> None:
    array = np.ascontiguousarray(np.asarray(value, dtype="<f8"))
    hasher.update(name.encode("ascii"))
    hasher.update(b"\0")
    hasher.update(np.asarray(array.shape, dtype="<i8").tobytes())
    hasher.update(array.tobytes(order="C"))


def _new_grid_hasher(time: np.ndarray, x: np.ndarray, y: np.ndarray, z: np.ndarray):
    hasher = hashlib.sha256()
    hasher.update(SCHEMA.encode("ascii"))
    hasher.update(b"\0")
    hasher.update(CANDIDATE_SHA256.encode("ascii"))
    hasher.update(b"\0time,x,y,z;u,v,w\0")
    _hash_array(hasher, "time", time)
    _hash_array(hasher, "x", x)
    _hash_array(hasher, "y", y)
    _hash_array(hasher, "z", z)
    return hasher


def _text(value) -> str:
    if isinstance(value, bytes):
        return value.decode("ascii")
    if isinstance(value, np.bytes_):
        return bytes(value).decode("ascii")
    return str(value)


def _readonly(value: np.ndarray) -> np.ndarray:
    value = np.asarray(value)
    value.setflags(write=False)
    return value


def export_st006_netcdf(
    path: str | os.PathLike[str],
    *,
    resolution: int = 33,
    overwrite: bool = False,
) -> dict:
    """Export ST006 on the registered box at fixed reference times.

    The file is NetCDF-3 64-bit-offset and stores named dimensions plus
    separate ``u``, ``v`` and ``w`` variables, so MATLAB ``ncread`` and
    Python/SciPy can consume the same artifact without component-order
    inference.  The result remains a visualization/delivery grid only.
    """

    resolution = _validate_resolution(resolution)
    target = Path(path)
    if target.suffix.lower() != ".nc":
        raise ValueError("NetCDF output path must end in .nc")
    if target.exists() and not overwrite:
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)

    axis = np.linspace(EVALUATION_BOUNDS[0], EVALUATION_BOUNDS[1], resolution, dtype=np.float64)
    time = REFERENCE_TIMES.copy()
    field = load_best()
    if field.sha256 != CANDIDATE_SHA256 or field.pde_validated is not False:
        raise ValueError("Unexpected retained candidate identity or truth state")

    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    hasher = _new_grid_hasher(time, axis, axis, axis)
    maximum_abs_velocity = 0.0
    temporary = target.with_name(f".{target.name}.{uuid.uuid4().hex}.tmp")

    try:
        with netcdf_file(str(temporary), mode="w", version=2) as nc:
            nc.createDimension("time", len(time))
            nc.createDimension("x", resolution)
            nc.createDimension("y", resolution)
            nc.createDimension("z", resolution)

            time_var = nc.createVariable("time", "d", ("time",))
            x_var = nc.createVariable("x", "d", ("x",))
            y_var = nc.createVariable("y", "d", ("y",))
            z_var = nc.createVariable("z", "d", ("z",))
            u_var = nc.createVariable("u", "d", ("time", "x", "y", "z"))
            v_var = nc.createVariable("v", "d", ("time", "x", "y", "z"))
            w_var = nc.createVariable("w", "d", ("time", "x", "y", "z"))
            speed_var = nc.createVariable("speed", "d", ("time", "x", "y", "z"))

            time_var[:] = time
            x_var[:] = axis
            y_var[:] = axis
            z_var[:] = axis
            time_var.long_name = "registered dimensionless time"
            x_var.long_name = "x coordinate"
            y_var.long_name = "y coordinate"
            z_var.long_name = "z coordinate"
            for variable, label in ((u_var, "u"), (v_var, "v"), (w_var, "w")):
                variable.long_name = f"Cartesian velocity component {label}"
                variable.units = "dimensionless"
            speed_var.long_name = "Cartesian velocity magnitude"
            speed_var.units = "dimensionless"

            for index, current_time in enumerate(time):
                velocity = np.asarray(field.velocity(xx, yy, zz, float(current_time)), dtype=np.float64)
                if velocity.shape != (resolution, resolution, resolution, 3):
                    raise ValueError("Retained velocity returned an unexpected shape")
                if not np.isfinite(velocity).all():
                    raise ValueError("Retained velocity returned nonfinite values")
                speed = np.linalg.norm(velocity, axis=-1)
                u_var[index] = velocity[..., 0]
                v_var[index] = velocity[..., 1]
                w_var[index] = velocity[..., 2]
                speed_var[index] = speed
                _hash_array(hasher, f"u[{index}]", velocity[..., 0])
                _hash_array(hasher, f"v[{index}]", velocity[..., 1])
                _hash_array(hasher, f"w[{index}]", velocity[..., 2])
                maximum_abs_velocity = max(maximum_abs_velocity, float(np.max(np.abs(velocity))))

            if maximum_abs_velocity == 0.0:
                raise ValueError("Exact-zero sampled velocity cannot be exported")

            digest = hasher.hexdigest()
            nc.schema = SCHEMA
            nc.candidate = "ST006"
            nc.candidate_sha256 = CANDIDATE_SHA256
            nc.grid_sha256 = digest
            nc.dimension_order = "time,x,y,z"
            nc.component_order = "u,v,w"
            nc.registered_box = "[-2,2]^3"
            nc.registered_time_window = "[0.25,0.75]"
            nc.reference_times = "0.25,0.50,0.75"
            nc.delivery_semantics = "visualization_only_not_pde_evidence"
            nc.velocity_export_ready = 1
            for flag in _FALSE_TRUTH_FLAGS:
                setattr(nc, flag, 0)

        os.replace(temporary, target)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise

    return {
        "path": str(target),
        "schema": SCHEMA,
        "candidate": "ST006",
        "candidate_sha256": CANDIDATE_SHA256,
        "grid_sha256": digest,
        "resolution": resolution,
        "times": tuple(float(value) for value in time),
        "dimension_order": "time,x,y,z",
        "component_order": "u,v,w",
        "velocity_export_ready": True,
        **{flag: False for flag in _FALSE_TRUTH_FLAGS},
    }


def load_st006_netcdf(path: str | os.PathLike[str]) -> ST006NetCDFGrid:
    """Load and fail-closed validate one ST006 NetCDF interchange artifact."""

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)

    with netcdf_file(str(source), mode="r", mmap=False) as nc:
        if _text(getattr(nc, "schema", "")) != SCHEMA:
            raise ValueError("Unexpected NetCDF schema")
        if _text(getattr(nc, "candidate", "")) != "ST006":
            raise ValueError("Unexpected candidate label")
        if _text(getattr(nc, "candidate_sha256", "")) != CANDIDATE_SHA256:
            raise ValueError("Candidate identity mismatch")
        if _text(getattr(nc, "dimension_order", "")) != "time,x,y,z":
            raise ValueError("Unexpected dimension order")
        if _text(getattr(nc, "component_order", "")) != "u,v,w":
            raise ValueError("Unexpected component order")
        if _text(getattr(nc, "delivery_semantics", "")) != "visualization_only_not_pde_evidence":
            raise ValueError("Unexpected delivery semantics")
        if int(getattr(nc, "velocity_export_ready", 0)) != 1:
            raise ValueError("Velocity export readiness was not recorded")
        for flag in _FALSE_TRUTH_FLAGS:
            if int(getattr(nc, flag, 1)) != 0:
                raise ValueError(f"Unsupported truth-state promotion: {flag}")

        required = {"time", "x", "y", "z", "u", "v", "w", "speed"}
        if set(nc.variables) != required:
            raise ValueError("Unexpected NetCDF variable set")
        if nc.variables["u"].dimensions != ("time", "x", "y", "z"):
            raise ValueError("Unexpected u dimensions")
        for name in ("v", "w", "speed"):
            if nc.variables[name].dimensions != nc.variables["u"].dimensions:
                raise ValueError(f"Unexpected {name} dimensions")

        time = np.array(nc.variables["time"].data, dtype=np.float64, copy=True)
        x = np.array(nc.variables["x"].data, dtype=np.float64, copy=True)
        y = np.array(nc.variables["y"].data, dtype=np.float64, copy=True)
        z = np.array(nc.variables["z"].data, dtype=np.float64, copy=True)
        u = np.array(nc.variables["u"].data, dtype=np.float64, copy=True)
        v = np.array(nc.variables["v"].data, dtype=np.float64, copy=True)
        w = np.array(nc.variables["w"].data, dtype=np.float64, copy=True)
        speed = np.array(nc.variables["speed"].data, dtype=np.float64, copy=True)
        stored_digest = _text(getattr(nc, "grid_sha256", ""))

    if not np.array_equal(time, REFERENCE_TIMES):
        raise ValueError("Reference times were changed")
    if any(array.ndim != 1 or not np.isfinite(array).all() for array in (x, y, z)):
        raise ValueError("Coordinate arrays must be finite one-dimensional axes")
    if not (len(x) == len(y) == len(z)):
        raise ValueError("Spatial axes must use one common resolution")
    resolution = _validate_resolution(len(x))
    for name, axis in (("x", x), ("y", y), ("z", z)):
        if not np.all(np.diff(axis) > 0.0):
            raise ValueError(f"{name} axis is not strictly increasing")
        if axis[0] != EVALUATION_BOUNDS[0] or axis[-1] != EVALUATION_BOUNDS[1]:
            raise ValueError(f"{name} axis leaves the registered evaluation box")

    expected_shape = (len(time), resolution, resolution, resolution)
    for name, array in (("u", u), ("v", v), ("w", w), ("speed", speed)):
        if array.shape != expected_shape or not np.isfinite(array).all():
            raise ValueError(f"{name} has invalid shape or values")
    recomputed_speed = np.sqrt(u * u + v * v + w * w)
    tolerance = 64.0 * np.finfo(np.float64).eps * max(1.0, float(np.max(recomputed_speed)))
    if float(np.max(np.abs(speed - recomputed_speed))) > tolerance:
        raise ValueError("Stored speed is inconsistent with [u,v,w]")
    if max(float(np.max(np.abs(u))), float(np.max(np.abs(v))), float(np.max(np.abs(w)))) == 0.0:
        raise ValueError("Exact-zero sampled velocity is not a valid delivery artifact")

    hasher = _new_grid_hasher(time, x, y, z)
    for index in range(len(time)):
        _hash_array(hasher, f"u[{index}]", u[index])
        _hash_array(hasher, f"v[{index}]", v[index])
        _hash_array(hasher, f"w[{index}]", w[index])
    if len(stored_digest) != 64 or hasher.hexdigest() != stored_digest:
        raise ValueError("NetCDF grid checksum mismatch")

    return ST006NetCDFGrid(
        time=_readonly(time),
        x=_readonly(x),
        y=_readonly(y),
        z=_readonly(z),
        u=_readonly(u),
        v=_readonly(v),
        w=_readonly(w),
        speed=_readonly(speed),
        grid_sha256=stored_digest,
    )
