"""Truth-preserving MATLAB export for frozen velocity fields.

This module is a visualization/export utility only.  A successful round trip
means that sampled velocity values were written and loaded consistently; it is
not Navier--Stokes validation, visual correspondence evidence, or identification
of OpenAI's underlying numerical field.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable, Mapping

import numpy as np
from scipy.io import loadmat, savemat


VelocityCallable = Callable[[np.ndarray, float], np.ndarray]
_REQUIRED_ARRAYS = ("x", "y", "z", "times", "U", "V", "W")
_AXIS_ORDER = "MATLAB meshgrid [y,x,z,time]"
_CLAIM_SCOPE = "visualization_export_only"


def _strict_vector(name: str, values: np.ndarray, *, minimum_size: int = 2) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size < minimum_size:
        raise ValueError(f"{name} must contain at least {minimum_size} values")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must be finite")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    return arr


def _matlab_string(value: object) -> str:
    arr = np.asarray(value)
    if arr.size == 0:
        return ""
    item = arr.squeeze()
    if isinstance(item, np.ndarray):
        item = item.item()
    return str(item)


def sample_velocity_grid(
    velocity: VelocityCallable,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    times: np.ndarray,
) -> dict[str, np.ndarray]:
    """Sample ``velocity(points, time)`` on a MATLAB-``meshgrid`` compatible grid.

    Returned ``U,V,W`` have shape ``(len(y), len(x), len(z), len(times))`` so
    that in MATLAB ``[X,Y,Z] = meshgrid(x,y,z)`` matches ``U(:,:,:,k)``.
    """

    x = _strict_vector("x", x)
    y = _strict_vector("y", y)
    z = _strict_vector("z", z)
    times = _strict_vector("times", times, minimum_size=1)

    X, Y, Z = np.meshgrid(x, y, z, indexing="xy")
    points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    shape = X.shape + (times.size,)

    U = np.empty(shape, dtype=float)
    V = np.empty(shape, dtype=float)
    W = np.empty(shape, dtype=float)

    for k, time in enumerate(times):
        values = np.asarray(velocity(points, float(time)), dtype=float)
        if values.shape != (points.shape[0], 3):
            raise ValueError(
                "velocity(points, time) must return shape "
                f"({points.shape[0]}, 3), got {values.shape}"
            )
        if not np.all(np.isfinite(values)):
            raise ValueError("velocity output must be finite")

        U[..., k] = values[:, 0].reshape(X.shape)
        V[..., k] = values[:, 1].reshape(X.shape)
        W[..., k] = values[:, 2].reshape(X.shape)

    return {"x": x, "y": y, "z": z, "times": times, "U": U, "V": V, "W": W}


def export_velocity_mat(
    path: str | Path,
    velocity: VelocityCallable,
    *,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    times: np.ndarray,
    candidate_id: str,
    metadata: Mapping[str, str | int | float] | None = None,
) -> Path:
    """Write a frozen velocity sample bundle as a MATLAB v5 ``.mat`` file.

    The file is intended for visualization and cross-language interchange.
    Export success must not be reported as PDE validation.
    """

    candidate_id = str(candidate_id).strip()
    if not candidate_id:
        raise ValueError("candidate_id must be non-empty")

    grid = sample_velocity_grid(velocity, x, y, z, times)
    payload: dict[str, object] = {
        **grid,
        "schema_version": np.array([[1]], dtype=np.int64),
        "candidate_id": candidate_id,
        "axis_order": _AXIS_ORDER,
        "claim_scope": _CLAIM_SCOPE,
        "pde_validated": np.array([[0]], dtype=np.uint8),
        "paper_exact": np.array([[0]], dtype=np.uint8),
        "openai_field_identified": np.array([[0]], dtype=np.uint8),
    }

    if metadata:
        for raw_key, raw_value in metadata.items():
            key = str(raw_key)
            if not key or not key.replace("_", "a").isalnum() or key[0].isdigit():
                raise ValueError(f"invalid MATLAB metadata key: {key!r}")
            full_key = f"meta_{key}"
            if full_key in payload:
                raise ValueError(f"metadata key collides with reserved field: {full_key}")
            if not isinstance(raw_value, (str, int, float, np.integer, np.floating)):
                raise TypeError(f"unsupported metadata value for {key!r}")
            payload[full_key] = raw_value

    target = Path(path)
    if target.suffix.lower() != ".mat":
        target = target.with_suffix(".mat")
    target.parent.mkdir(parents=True, exist_ok=True)

    savemat(target, payload, format="5", do_compression=True, oned_as="row")
    return target


def load_velocity_mat(path: str | Path) -> dict[str, object]:
    """Load and validate a bundle written by :func:`export_velocity_mat`."""

    raw = loadmat(Path(path), squeeze_me=True, struct_as_record=False)
    missing = [name for name in _REQUIRED_ARRAYS if name not in raw]
    if missing:
        raise ValueError(f"missing required MATLAB fields: {missing}")

    x = np.asarray(raw["x"], dtype=float).reshape(-1)
    y = np.asarray(raw["y"], dtype=float).reshape(-1)
    z = np.asarray(raw["z"], dtype=float).reshape(-1)
    times = np.asarray(raw["times"], dtype=float).reshape(-1)
    _strict_vector("x", x)
    _strict_vector("y", y)
    _strict_vector("z", z)
    _strict_vector("times", times, minimum_size=1)

    expected = (y.size, x.size, z.size, times.size)
    arrays: dict[str, np.ndarray] = {}
    for name in ("U", "V", "W"):
        arr = np.asarray(raw[name], dtype=float)
        # scipy.squeeze_me removes a singleton time dimension, so restore it.
        if times.size == 1 and arr.shape == expected[:-1]:
            arr = arr[..., None]
        if arr.shape != expected:
            raise ValueError(f"{name} has shape {arr.shape}, expected {expected}")
        if not np.all(np.isfinite(arr)):
            raise ValueError(f"{name} must be finite")
        arrays[name] = arr

    if _matlab_string(raw.get("axis_order", "")) != _AXIS_ORDER:
        raise ValueError("unexpected or missing axis_order")
    if _matlab_string(raw.get("claim_scope", "")) != _CLAIM_SCOPE:
        raise ValueError("unexpected or missing claim_scope")
    candidate_id = _matlab_string(raw.get("candidate_id", "")).strip()
    if not candidate_id:
        raise ValueError("missing candidate_id")

    return {
        "x": x,
        "y": y,
        "z": z,
        "times": times,
        **arrays,
        "candidate_id": candidate_id,
        "axis_order": _AXIS_ORDER,
        "claim_scope": _CLAIM_SCOPE,
    }
