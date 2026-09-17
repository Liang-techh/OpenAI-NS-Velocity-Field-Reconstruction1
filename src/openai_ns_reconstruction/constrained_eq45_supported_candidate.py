"""Serializable Eq. (4.5) candidate with a physical-space exterior connection.

This module composes three already-governed pieces without changing any Eq45
profile coefficient:

* :class:`Eq45VelocityCandidate` supplies the callable Eq. (4.5) field;
* :func:`eq45_candidate_cylindrical_state` exposes its matching streamfunction
  and cylindrical components; and
* :class:`AxisymmetricPhysicalTaper` applies the C4 streamfunction-level
  exterior connection in physical ``(r,z)`` coordinates.

The resulting object is a *child candidate*: the exterior collar changes
``[u,v,w]``, so it receives a new identity and must be revalidated on its own.
The construction gives a callable/saveable velocity with exact zero extension
at/outside the declared taper cylinder and exact equality with the parent on
the taper plateau.  Those representation properties do not by themselves
establish visualization correspondence, Navier--Stokes validity, paper
exactness, or identification of an OpenAI hidden field.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_axisymmetric_physical_taper import AxisymmetricPhysicalTaper
from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_cylindrical_state import eq45_candidate_cylindrical_state


SCHEMA = "eq45_supported_velocity_candidate_v1"

_CLASSIFICATION = {
    "velocity_interface": "user_requirement",
    "eq45_parent": "parent_provenance",
    "physical_support_transform": "autonomous_design",
    "taper_parameters": "autonomous_design",
    "openai_hidden_profiles": "pending_unknown",
}

_TRUTH_BOUNDARY = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
    "requires_physical_support_connection": False,
}


@dataclass(frozen=True)
class Eq45SupportedVelocityCandidate:
    """Eq45 child candidate after the governed physical support transform."""

    parent: Eq45VelocityCandidate
    taper: AxisymmetricPhysicalTaper = AxisymmetricPhysicalTaper()

    def __post_init__(self) -> None:
        if not isinstance(self.parent, Eq45VelocityCandidate):
            raise TypeError("parent must be Eq45VelocityCandidate")
        if not isinstance(self.taper, AxisymmetricPhysicalTaper):
            raise TypeError("taper must be AxisymmetricPhysicalTaper")

    @property
    def time_start(self) -> float:
        return float(self.parent.time_start)

    @property
    def time_end(self) -> float:
        return float(self.parent.time_end)

    @property
    def parent_sha256(self) -> str:
        return self.parent.sha256

    def velocity(self, points, time) -> np.ndarray:
        """Evaluate the support-connected candidate as Cartesian ``[...,3]``."""
        points_arr = np.asarray(points, dtype=float)
        state = eq45_candidate_cylindrical_state(self.parent, points_arr, time)
        out = self.taper.apply_cylindrical(
            points_arr,
            **state.taper_inputs(),
        )
        if not np.all(np.isfinite(out)):
            raise RuntimeError("support-connected Eq45 velocity became nonfinite")
        return out

    def at_points(self, points, time) -> np.ndarray:
        """Public validator/visualizer alias for :meth:`velocity`."""
        return self.velocity(points, time)

    __call__ = velocity

    def velocity_xyz(self, x, y, z, time) -> np.ndarray:
        """Broadcast x/y/z/time and evaluate the support-connected field."""
        arrays = [np.asarray(value, dtype=float) for value in (x, y, z, time)]
        if not all(np.all(np.isfinite(value)) for value in arrays):
            raise ValueError("x, y, z and time must be finite")
        try:
            x_arr, y_arr, z_arr, time_arr = np.broadcast_arrays(*arrays)
        except ValueError as exc:
            raise ValueError("x, y, z and time must be broadcast-compatible") from exc
        points = np.stack((x_arr, y_arr, z_arr), axis=-1)
        return self.velocity(points, time_arr)

    def grid(self, x, y, z, times) -> np.ndarray:
        """Return ``(time,x,y,z,component)`` samples for visualization/export."""
        axes = [np.asarray(value, dtype=float) for value in (times, x, y, z)]
        if any(axis.ndim != 1 or axis.size == 0 for axis in axes):
            raise ValueError("grid axes and times must be nonempty 1D arrays")
        if not all(np.all(np.isfinite(axis)) for axis in axes):
            raise ValueError("grid axes and times must be finite")
        tt, xx, yy, zz = np.meshgrid(*axes, indexing="ij")
        return self.velocity_xyz(xx, yy, zz, tt)

    @property
    def sha256(self) -> str:
        payload = json.dumps(
            self.to_dict(), sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode("utf-8")
        return hashlib.sha256(payload).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        """Return a fail-closed child artifact including explicit parent identity."""
        return {
            "schema": SCHEMA,
            "classification": dict(_CLASSIFICATION),
            "parent_sha256": self.parent.sha256,
            "parent_candidate": self.parent.to_dict(),
            "physical_taper": self.taper.to_dict(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "Eq45SupportedVelocityCandidate":
        if not isinstance(data, Mapping):
            raise ValueError("supported candidate payload must be an object")
        if data.get("schema") != SCHEMA:
            raise ValueError(f"schema must be {SCHEMA!r}")
        if data.get("classification") != _CLASSIFICATION:
            raise ValueError("supported candidate source classification was modified")
        if data.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("supported candidate truth-boundary metadata was modified")
        if "parent_candidate" not in data or "physical_taper" not in data:
            raise ValueError("supported candidate payload is missing a required component")

        parent = Eq45VelocityCandidate.from_dict(data["parent_candidate"])
        claimed_parent_sha = data.get("parent_sha256")
        if claimed_parent_sha != parent.sha256:
            raise ValueError("parent_sha256 does not match the embedded parent candidate")
        taper = AxisymmetricPhysicalTaper.from_dict(data["physical_taper"])
        return cls(parent=parent, taper=taper)

    def save_json(self, path: str | Path) -> None:
        target = Path(path)
        target.write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(cls, path: str | Path) -> "Eq45SupportedVelocityCandidate":
        target = Path(path)
        return cls.from_dict(json.loads(target.read_text(encoding="utf-8")))
