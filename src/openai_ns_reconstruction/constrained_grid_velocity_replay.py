"""Truth-bounded callable replay of a sampled time-varying Cartesian velocity grid.

This adapter is for visualization/delivery only.  Linear interpolation of a
sampled grid is not an independent PDE, divergence, or source-correspondence
validator and must not be used to promote those states.
"""

from __future__ import annotations

from hashlib import sha256
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Any

import numpy as np
from scipy.interpolate import RegularGridInterpolator

_SCHEMA_VERSION = 1
_COMPONENT_ORDER = ("u", "v", "w")
_GRID_LAYOUT = "time,x,y,z,component"
_SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
_REQUIRED_KEYS = {
    "schema_version",
    "x",
    "y",
    "z",
    "t",
    "velocity_txyzc",
    "candidate_sha256",
    "provenance",
    "grid_sha256",
}


def _validated_axis(name: str, values: Any) -> np.ndarray:
    axis = np.asarray(values, dtype=np.float64)
    if axis.ndim != 1 or axis.size < 2:
        raise ValueError(f"{name} must be a 1-D axis with at least two points")
    if not np.all(np.isfinite(axis)):
        raise ValueError(f"{name} must be finite")
    if not np.all(np.diff(axis) > 0.0):
        raise ValueError(f"{name} must be strictly increasing")
    axis = np.array(axis, dtype=np.float64, copy=True)
    axis.setflags(write=False)
    return axis


def _validated_text(name: str, value: Any) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value.strip()


def _grid_digest(
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    velocity: np.ndarray,
    candidate_sha256: str,
    provenance: str,
) -> str:
    h = sha256()
    h.update(f"grid_velocity_replay_v{_SCHEMA_VERSION}".encode("ascii"))
    for name, arr in (("x", x), ("y", y), ("z", z), ("t", t), ("velocity", velocity)):
        canonical = np.ascontiguousarray(arr, dtype="<f8")
        h.update(name.encode("ascii"))
        h.update(np.asarray(canonical.shape, dtype="<i8").tobytes())
        h.update(canonical.tobytes())
    h.update(candidate_sha256.encode("ascii"))
    h.update(provenance.encode("utf-8"))
    return h.hexdigest()


class GridVelocityReplay:
    """Multilinear callable view of a frozen ``(t,x,y,z,3)`` velocity grid."""

    def __init__(
        self,
        *,
        x: Any,
        y: Any,
        z: Any,
        t: Any,
        velocity_txyzc: Any,
        candidate_sha256: str,
        provenance: str,
    ) -> None:
        self.x = _validated_axis("x", x)
        self.y = _validated_axis("y", y)
        self.z = _validated_axis("z", z)
        self.t = _validated_axis("t", t)
        candidate_sha256 = _validated_text("candidate_sha256", candidate_sha256)
        if _SHA256_RE.fullmatch(candidate_sha256) is None:
            raise ValueError("candidate_sha256 must be 64 lowercase hexadecimal characters")
        self.candidate_sha256 = candidate_sha256
        self.provenance = _validated_text("provenance", provenance)

        velocity = np.asarray(velocity_txyzc, dtype=np.float64)
        expected = (self.t.size, self.x.size, self.y.size, self.z.size, 3)
        if velocity.shape != expected:
            raise ValueError(f"velocity_txyzc must have shape {expected}, got {velocity.shape}")
        if not np.all(np.isfinite(velocity)):
            raise ValueError("velocity_txyzc must be finite")
        if not np.any(np.linalg.norm(velocity, axis=-1) > 0.0):
            raise ValueError("velocity_txyzc must be nontrivial; the exact zero field is rejected")
        self.velocity_txyzc = np.array(velocity, dtype=np.float64, copy=True)
        self.velocity_txyzc.setflags(write=False)

        self.grid_sha256 = _grid_digest(
            self.x,
            self.y,
            self.z,
            self.t,
            self.velocity_txyzc,
            self.candidate_sha256,
            self.provenance,
        )
        self._interpolator = RegularGridInterpolator(
            (self.t, self.x, self.y, self.z),
            self.velocity_txyzc,
            method="linear",
            bounds_error=True,
            fill_value=np.nan,
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Return interpolated Cartesian ``[...,3]`` values in ``[u,v,w]`` order."""
        xb, yb, zb, tb = np.broadcast_arrays(
            np.asarray(x, dtype=np.float64),
            np.asarray(y, dtype=np.float64),
            np.asarray(z, dtype=np.float64),
            np.asarray(t, dtype=np.float64),
        )
        for name, arr in (("x", xb), ("y", yb), ("z", zb), ("t", tb)):
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"query {name} must be finite")
        q = np.column_stack((tb.ravel(), xb.ravel(), yb.ravel(), zb.ravel()))
        try:
            out = np.asarray(self._interpolator(q), dtype=np.float64)
        except ValueError as exc:
            raise ValueError("query lies outside the frozen grid; extrapolation is forbidden") from exc
        out = out.reshape(xb.shape + (3,))
        if not np.all(np.isfinite(out)):
            raise RuntimeError("interpolator returned non-finite velocity")
        return out

    __call__ = velocity

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        """Evaluate ``velocity`` at an array whose final dimension is ``(x,y,z)``."""
        points = np.asarray(points_xyz, dtype=np.float64)
        if points.ndim < 1 or points.shape[-1] != 3:
            raise ValueError("points_xyz must have final dimension 3")
        if not np.all(np.isfinite(points)):
            raise ValueError("points_xyz must be finite")
        lead = points.shape[:-1]
        tb = np.broadcast_to(np.asarray(t, dtype=np.float64), lead)
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], tb)

    def metadata(self) -> dict[str, Any]:
        return {
            "schema_version": _SCHEMA_VERSION,
            "layout": _GRID_LAYOUT,
            "component_order": list(_COMPONENT_ORDER),
            "interpolation": "scipy.interpolate.RegularGridInterpolator(method='linear', bounds_error=True)",
            "candidate_sha256": self.candidate_sha256,
            "grid_sha256": self.grid_sha256,
            "provenance": self.provenance,
            "velocity_export_ready": True,
            "visualization_interpolation_only": True,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
            "interpolated_derivatives_valid_for_pde_acceptance": False,
            "extrapolation_allowed": False,
        }

    def save(self, path: str | os.PathLike[str], *, overwrite: bool = False) -> Path:
        target = Path(path)
        if target.exists() and not overwrite:
            raise FileExistsError(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema_version": np.asarray(_SCHEMA_VERSION, dtype=np.int64),
            "x": self.x,
            "y": self.y,
            "z": self.z,
            "t": self.t,
            "velocity_txyzc": self.velocity_txyzc,
            "candidate_sha256": np.asarray(self.candidate_sha256),
            "provenance": np.asarray(self.provenance),
            "grid_sha256": np.asarray(self.grid_sha256),
        }
        fd, tmp_name = tempfile.mkstemp(prefix=target.name + ".", suffix=".tmp", dir=target.parent)
        os.close(fd)
        tmp = Path(tmp_name)
        try:
            with tmp.open("wb") as handle:
                np.savez_compressed(handle, **payload)
            loaded = type(self).load(tmp)
            if loaded.grid_sha256 != self.grid_sha256:
                raise RuntimeError("saved grid failed identity round-trip")
            os.replace(tmp, target)
        finally:
            if tmp.exists():
                tmp.unlink()
        return target

    @classmethod
    def load(cls, path: str | os.PathLike[str]) -> "GridVelocityReplay":
        source = Path(path)
        with np.load(source, allow_pickle=False) as data:
            if set(data.files) != _REQUIRED_KEYS:
                raise ValueError("grid artifact keys do not match the governed schema")
            schema = int(np.asarray(data["schema_version"]).item())
            if schema != _SCHEMA_VERSION:
                raise ValueError(f"unsupported schema_version {schema}")
            stored_grid_sha = str(np.asarray(data["grid_sha256"]).item())
            obj = cls(
                x=np.array(data["x"], copy=True),
                y=np.array(data["y"], copy=True),
                z=np.array(data["z"], copy=True),
                t=np.array(data["t"], copy=True),
                velocity_txyzc=np.array(data["velocity_txyzc"], copy=True),
                candidate_sha256=str(np.asarray(data["candidate_sha256"]).item()),
                provenance=str(np.asarray(data["provenance"]).item()),
            )
        if _SHA256_RE.fullmatch(stored_grid_sha) is None or stored_grid_sha != obj.grid_sha256:
            raise ValueError("grid_sha256 mismatch; artifact bytes or metadata changed")
        return obj


def metadata_json(field: GridVelocityReplay) -> str:
    """Stable compact JSON for manifests without serializing NumPy objects."""
    return json.dumps(field.metadata(), sort_keys=True, separators=(",", ":"))
