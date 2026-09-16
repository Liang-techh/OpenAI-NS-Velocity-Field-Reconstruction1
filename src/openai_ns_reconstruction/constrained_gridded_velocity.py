"""Truth-bounded structured-grid velocity sampling and interpolation.

This module turns a frozen callable velocity field into a portable rectilinear
space-time grid and reconstructs a callable ``velocity(points, time)`` interface
from that saved grid. Interpolation is a visualization/delivery convenience only:
it is never evidence that the underlying candidate satisfies Navier--Stokes.
"""
from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import json

import numpy as np
from scipy.interpolate import RegularGridInterpolator


_SCHEMA_VERSION = 1
_AXIS_ORDER = "t,x,y,z,component"
_CLAIM_SCOPE = "visualization_interpolation_only"


def _strict_axis(values, name: str, *, minimum_size: int = 2) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 1 or arr.size < minimum_size:
        raise ValueError(f"{name} must be a one-dimensional axis with at least {minimum_size} points")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(np.diff(arr) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    out = np.array(arr, copy=True)
    out.setflags(write=False)
    return out


def _metadata(candidate_id: str) -> dict:
    if not isinstance(candidate_id, str) or not candidate_id.strip():
        raise ValueError("candidate_id must be a nonempty string")
    return {
        "schema_version": _SCHEMA_VERSION,
        "candidate_id": candidate_id,
        "axis_order": _AXIS_ORDER,
        "interpolation": "scipy.interpolate.RegularGridInterpolator(method=linear,bounds_error=True)",
        "claim_scope": _CLAIM_SCOPE,
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def _validate_metadata(meta: dict) -> dict:
    required = {
        "schema_version": _SCHEMA_VERSION,
        "axis_order": _AXIS_ORDER,
        "claim_scope": _CLAIM_SCOPE,
        "velocity_export_ready": True,
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    if not isinstance(meta, dict):
        raise ValueError("metadata must be a JSON object")
    for key, expected in required.items():
        if meta.get(key) != expected:
            raise ValueError(f"metadata field {key!r} must equal {expected!r}")
    candidate_id = meta.get("candidate_id")
    if not isinstance(candidate_id, str) or not candidate_id.strip():
        raise ValueError("metadata candidate_id must be nonempty")
    interpolation = meta.get("interpolation")
    if interpolation != "scipy.interpolate.RegularGridInterpolator(method=linear,bounds_error=True)":
        raise ValueError("unexpected interpolation contract")
    return dict(meta)


@dataclass(frozen=True)
class StructuredVelocityField:
    """Frozen rectilinear space-time velocity grid with strict linear interpolation."""

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    times: np.ndarray
    values: np.ndarray
    metadata: dict

    def __post_init__(self):
        x = _strict_axis(self.x, "x")
        y = _strict_axis(self.y, "y")
        z = _strict_axis(self.z, "z")
        times = _strict_axis(self.times, "times")
        values = np.asarray(self.values, dtype=float)
        expected = (times.size, x.size, y.size, z.size, 3)
        if values.shape != expected:
            raise ValueError(f"values must have shape {expected}, got {values.shape}")
        if not np.all(np.isfinite(values)):
            raise ValueError("values must contain only finite velocity components")
        meta = _validate_metadata(self.metadata)

        frozen_values = np.array(values, copy=True)
        frozen_values.setflags(write=False)
        object.__setattr__(self, "x", x)
        object.__setattr__(self, "y", y)
        object.__setattr__(self, "z", z)
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "values", frozen_values)
        object.__setattr__(self, "metadata", meta)
        object.__setattr__(self, "_interpolator", RegularGridInterpolator(
            (times, x, y, z),
            frozen_values,
            method="linear",
            bounds_error=True,
        ))

    @property
    def candidate_id(self) -> str:
        return self.metadata["candidate_id"]

    @classmethod
    def sample(
        cls,
        velocity,
        *,
        x,
        y,
        z,
        times,
        candidate_id: str,
    ) -> "StructuredVelocityField":
        """Sample a frozen ``velocity(points,time)->[u,v,w]`` callable on a grid."""
        x = _strict_axis(x, "x")
        y = _strict_axis(y, "y")
        z = _strict_axis(z, "z")
        times = _strict_axis(times, "times")
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
        rows = []
        for time in times:
            out = np.asarray(velocity(points, float(time)), dtype=float)
            if out.shape != (points.shape[0], 3):
                raise ValueError("velocity callable must return shape (n_points,3)")
            if not np.all(np.isfinite(out)):
                raise FloatingPointError("velocity callable returned non-finite values")
            rows.append(out.reshape(x.size, y.size, z.size, 3))
        values = np.stack(rows, axis=0)
        return cls(x=x, y=y, z=z, times=times, values=values, metadata=_metadata(candidate_id))

    def at_points(self, points, time: float) -> np.ndarray:
        """Evaluate interpolated velocity at finite points for one in-range time."""
        pts = np.asarray(points, dtype=float)
        if pts.ndim != 2 or pts.shape[1] != 3 or not np.all(np.isfinite(pts)):
            raise ValueError("points must be a finite array with shape (n,3)")
        if not np.isfinite(time):
            raise ValueError("time must be finite")
        q = np.column_stack((np.full(pts.shape[0], float(time)), pts))
        return np.asarray(self._interpolator(q), dtype=float)

    def __call__(self, points, time: float) -> np.ndarray:
        return self.at_points(points, time)

    def grid(self, x, y, z, times) -> np.ndarray:
        """Evaluate the interpolated field on a requested rectilinear grid."""
        xq = _strict_axis(x, "x")
        yq = _strict_axis(y, "y")
        zq = _strict_axis(z, "z")
        tq = _strict_axis(times, "times")
        X, Y, Z = np.meshgrid(xq, yq, zq, indexing="ij")
        xyz = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
        out = []
        for time in tq:
            out.append(self.at_points(xyz, float(time)).reshape(xq.size, yq.size, zq.size, 3))
        return np.stack(out, axis=0)

    def save_npz(self, path) -> str:
        """Save the frozen grid without pickle and return the file SHA256."""
        path = Path(path)
        metadata_json = json.dumps(self.metadata, sort_keys=True, separators=(",", ":"))
        np.savez_compressed(
            path,
            x=self.x,
            y=self.y,
            z=self.z,
            times=self.times,
            values=self.values,
            metadata_json=np.asarray(metadata_json),
        )
        return sha256(path.read_bytes()).hexdigest()

    @classmethod
    def load_npz(cls, path) -> "StructuredVelocityField":
        """Load and fail-close validate a saved grid bundle."""
        path = Path(path)
        with np.load(path, allow_pickle=False) as data:
            required = {"x", "y", "z", "times", "values", "metadata_json"}
            if set(data.files) != required:
                raise ValueError(f"grid bundle must contain exactly {sorted(required)}")
            raw_meta = data["metadata_json"]
            if raw_meta.shape != ():
                raise ValueError("metadata_json must be a scalar string")
            try:
                meta = json.loads(str(raw_meta.item()))
            except (TypeError, ValueError, json.JSONDecodeError) as exc:
                raise ValueError("metadata_json is not valid JSON") from exc
            return cls(
                x=np.array(data["x"], copy=True),
                y=np.array(data["y"], copy=True),
                z=np.array(data["z"], copy=True),
                times=np.array(data["times"], copy=True),
                values=np.array(data["values"], copy=True),
                metadata=meta,
            )
