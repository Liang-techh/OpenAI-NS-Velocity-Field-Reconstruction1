"""Governed Python visualization-grid bundle for the retained ST006 field.

This module samples only the published ST006 callable on the registered box and
three declared reference times.  The resulting grid is a visualization/export
artifact, not Navier--Stokes acceptance evidence.
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import tempfile

import numpy as np

from . import CANDIDATE_SHA256, load_best

SCHEMA = "st006_visualization_grid_v1"
LAYOUT = "time,x,y,z,component"
COMPONENT_ORDER = "u,v,w"
CLAIM_SCOPE = "retained_st006_visualization_sampling_only"
REFERENCE_TIMES = np.array([0.25, 0.50, 0.75], dtype=np.float64)
BOX = (-2.0, 2.0)
TRUTH_FLAGS = (
    "visualization_ready",
    "visual_correspondence_verified",
    "pde_validated",
    "paper_exact",
    "openai_field_identified",
    "blowup_proved",
)


def _readonly(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a)
    a.setflags(write=False)
    return a


def _validate_resolution(resolution: int) -> int:
    if isinstance(resolution, bool) or not isinstance(resolution, (int, np.integer)):
        raise ValueError("resolution must be an odd integer")
    resolution = int(resolution)
    if resolution < 5 or resolution > 129 or resolution % 2 == 0:
        raise ValueError("resolution must be odd and lie in [5,129]")
    return resolution


def _grid_sha256(x: np.ndarray, y: np.ndarray, z: np.ndarray, times: np.ndarray,
                 velocity: np.ndarray) -> str:
    h = hashlib.sha256()
    for text in (SCHEMA, CANDIDATE_SHA256, LAYOUT, COMPONENT_ORDER, CLAIM_SCOPE):
        h.update(text.encode("utf-8")); h.update(b"\0")
    for array in (x, y, z, times, velocity):
        a = np.ascontiguousarray(array, dtype="<f8")
        h.update(str(a.shape).encode("ascii")); h.update(b"\0"); h.update(a.tobytes(order="C"))
    return h.hexdigest()


def sample_reference_grid(resolution: int = 17) -> dict:
    """Sample ST006 on the fixed ``[-2,2]^3`` box at ``.25/.50/.75``.

    No time, camera, registration, threshold, forcing, or candidate parameter is
    fitted here.  Odd grids are required so the physical origin is represented.
    """
    resolution = _validate_resolution(resolution)
    x = np.linspace(BOX[0], BOX[1], resolution, dtype=np.float64)
    y = x.copy(); z = x.copy(); times = REFERENCE_TIMES.copy()
    xx, yy, zz = np.meshgrid(x, y, z, indexing="ij")
    points = np.stack((xx, yy, zz), axis=-1)
    field = load_best()
    velocity = np.stack([field.at_points(points, float(t)) for t in times], axis=0)
    if velocity.shape != (3, resolution, resolution, resolution, 3):
        raise ValueError("published velocity returned an unexpected grid shape")
    if not np.isfinite(velocity).all():
        raise ValueError("published velocity grid contains non-finite values")
    if not np.any(velocity != 0.0):
        raise ValueError("refusing an exact-zero visualization grid")
    speed = np.linalg.norm(velocity, axis=-1)
    digest = _grid_sha256(x, y, z, times, velocity)
    return {
        "schema": SCHEMA,
        "candidate_sha256": CANDIDATE_SHA256,
        "grid_sha256": digest,
        "layout": LAYOUT,
        "component_order": COMPONENT_ORDER,
        "claim_scope": CLAIM_SCOPE,
        "x": _readonly(x), "y": _readonly(y), "z": _readonly(z),
        "times": _readonly(times),
        "velocity": _readonly(velocity),
        "speed": _readonly(speed),
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def write_visualization_bundle(path, resolution: int = 17, *, overwrite: bool = False) -> dict:
    """Write one pickle-free compressed NPZ visualization bundle atomically."""
    path = Path(path)
    if path.suffix.lower() != ".npz":
        raise ValueError("visualization bundle path must end in .npz")
    if path.exists() and not overwrite:
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    grid = sample_reference_grid(resolution)
    payload = {k: v for k, v in grid.items() if isinstance(v, np.ndarray)}
    payload.update({
        "schema": np.array(grid["schema"]),
        "candidate_sha256": np.array(grid["candidate_sha256"]),
        "grid_sha256": np.array(grid["grid_sha256"]),
        "layout": np.array(grid["layout"]),
        "component_order": np.array(grid["component_order"]),
        "claim_scope": np.array(grid["claim_scope"]),
    })
    for flag in TRUTH_FLAGS:
        payload[flag] = np.array(False, dtype=np.bool_)
    fd, tmp_name = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp.npz", dir=path.parent)
    os.close(fd)
    tmp = Path(tmp_name)
    try:
        np.savez_compressed(tmp, **payload)
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()
    return {k: grid[k] for k in (
        "schema", "candidate_sha256", "grid_sha256", "layout", "component_order",
        "claim_scope", *TRUTH_FLAGS,
    )}


def _scalar(data, key: str):
    value = data[key]
    if value.shape != ():
        raise ValueError(f"{key} must be a scalar")
    return value.item()


def load_visualization_bundle(path) -> dict:
    """Load and fail-close on identity, layout, hash, or truth-state drift."""
    required = {
        "schema", "candidate_sha256", "grid_sha256", "layout", "component_order", "claim_scope",
        "x", "y", "z", "times", "velocity", "speed", *TRUTH_FLAGS,
    }
    with np.load(Path(path), allow_pickle=False) as data:
        if set(data.files) != required:
            raise ValueError("unexpected visualization bundle fields")
        scalars = {k: _scalar(data, k) for k in (
            "schema", "candidate_sha256", "grid_sha256", "layout", "component_order", "claim_scope",
            *TRUTH_FLAGS,
        )}
        arrays = {k: np.array(data[k], dtype=np.float64, copy=True) for k in (
            "x", "y", "z", "times", "velocity", "speed",
        )}
    if scalars["schema"] != SCHEMA or scalars["candidate_sha256"] != CANDIDATE_SHA256:
        raise ValueError("visualization bundle candidate/schema identity mismatch")
    if scalars["layout"] != LAYOUT or scalars["component_order"] != COMPONENT_ORDER:
        raise ValueError("visualization bundle layout/component order mismatch")
    if scalars["claim_scope"] != CLAIM_SCOPE:
        raise ValueError("visualization bundle claim scope mismatch")
    if any(bool(scalars[flag]) for flag in TRUTH_FLAGS):
        raise ValueError("unsupported visualization/scientific claim promotion")
    x, y, z, times = arrays["x"], arrays["y"], arrays["z"], arrays["times"]
    for name, axis in (("x", x), ("y", y), ("z", z)):
        if axis.ndim != 1 or len(axis) < 5 or len(axis) % 2 == 0 or not np.isfinite(axis).all():
            raise ValueError(f"invalid {name} axis")
        if not np.all(np.diff(axis) > 0) or axis[0] != BOX[0] or axis[-1] != BOX[1]:
            raise ValueError(f"{name} axis left the registered visualization box")
    if not (np.array_equal(x, y) and np.array_equal(x, z)):
        raise ValueError("reference visualization axes must use one shared resolution")
    if not np.array_equal(times, REFERENCE_TIMES):
        raise ValueError("reference visualization times changed")
    velocity, speed = arrays["velocity"], arrays["speed"]
    n = len(x)
    if velocity.shape != (3, n, n, n, 3) or speed.shape != (3, n, n, n):
        raise ValueError("visualization array shape/layout mismatch")
    if not np.isfinite(velocity).all() or not np.isfinite(speed).all() or not np.any(velocity != 0.0):
        raise ValueError("invalid visualization velocity/speed values")
    if not np.array_equal(speed, np.linalg.norm(velocity, axis=-1)):
        raise ValueError("visualization speed does not match [u,v,w]")
    actual = _grid_sha256(x, y, z, times, velocity)
    if scalars["grid_sha256"] != actual:
        raise ValueError("visualization grid checksum mismatch")
    result = {**scalars, **arrays}
    for key in arrays:
        result[key] = _readonly(result[key])
    return result
