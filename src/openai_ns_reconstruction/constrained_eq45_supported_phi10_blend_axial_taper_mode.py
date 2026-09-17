"""Bounded axial-taper representation for the supported Eq45 Phi10 blend.

Agent-7 PR #193 showed that the existing compact--quartic ``Phi(1,0)`` blend
and ``Phi(1,2)`` profile coefficient are both strongly suppressed at the
extreme top/bottom physical-support collar.  This module materializes the
smallest representation response to that diagnosed bottleneck: vary only the
axial *plateau* of the already-governed divergence-preserving physical taper.

The Eq. (4.5) parent, compact--quartic temporal schedule, blend weight, radial
taper, and outer support faces are unchanged.  The caller supplies one bounded
``axial_plateau_q`` in ``[0.64, 0.81]``.  For the governed half-height ``H=2``
this moves the exact-identity axial plateau from ``|z|<=1.6`` to at most
``|z|<=1.8`` while keeping exact zero at/outside ``|z|=2``.  No value is
selected here and no public image, PDE residual, force, pressure, or hidden
OpenAI parameter enters the construction.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_axisymmetric_physical_taper import AxisymmetricPhysicalTaper
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)


SCHEMA = "eq45_supported_phi10_blend_axial_taper_candidate_v1"
TASK_ID = "CR003-EQ45-SUPPORTED-BLEND-AXIAL-TAPER-MODE-034"
EXPECTED_SUPPORTED_BASE_SHA256 = (
    "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
)
BASE_AXIAL_PLATEAU_Q = 0.64
MAX_AXIAL_PLATEAU_Q = 0.81
AXIAL_PLATEAU_Q_BOUNDS = (BASE_AXIAL_PLATEAU_Q, MAX_AXIAL_PLATEAU_Q)

_CLASSIFICATION = {
    "velocity_interface": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "base_blend_candidate": "parent_provenance",
    "physical_support_transform": "autonomous_design",
    "axial_taper_plateau_control": "autonomous_design",
    "axial_taper_bound": "autonomous_design",
    "openai_hidden_profiles": "pending_unknown",
}

_TRUTH_BOUNDARY = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "axial_taper_control_materialized": True,
    "axial_taper_value_selected_by_this_module": False,
    "canonical_velocity_changed": False,
    "spatial_profile_basis_grown": False,
    "temporal_family_grown": False,
    "radial_taper_changed": False,
    "physical_support_outer_faces_changed": False,
    "force_or_pressure_fitted": False,
    "pde_objective_used_to_choose_axial_taper": False,
    "public_image_fitted": False,
    "physical_support_validated": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validate_axial_plateau_q(value: float) -> float:
    checked = float(value)
    if not np.isfinite(checked) or not (
        BASE_AXIAL_PLATEAU_Q <= checked <= MAX_AXIAL_PLATEAU_Q
    ):
        raise ValueError(
            "axial_plateau_q must be finite and lie in the autonomous trial bound "
            f"[{BASE_AXIAL_PLATEAU_Q}, {MAX_AXIAL_PLATEAU_Q}]"
        )
    return checked


@dataclass(frozen=True)
class Eq45SupportedPhi10BlendAxialTaperCandidate:
    """Compact--quartic Eq45 blend with one caller-declared axial taper control."""

    base: Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate
    axial_plateau_q: float = BASE_AXIAL_PLATEAU_Q

    def __post_init__(self) -> None:
        if not isinstance(
            self.base, Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate
        ):
            raise TypeError(
                "base must be Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate"
            )
        object.__setattr__(
            self, "axial_plateau_q", _validate_axial_plateau_q(self.axial_plateau_q)
        )

        supported = self.base.base
        if supported.sha256 != EXPECTED_SUPPORTED_BASE_SHA256:
            raise ValueError("axial taper family requires the frozen governed supported base")
        taper = supported.taper
        if (
            taper.radial_support != 2.0
            or taper.axial_half_height != 2.0
            or taper.radial_plateau_q != 0.64
            or taper.axial_plateau_q != BASE_AXIAL_PLATEAU_Q
        ):
            raise ValueError("base physical taper no longer matches the governed support contract")

    @property
    def time_start(self) -> float:
        return float(self.base.time_start)

    @property
    def time_end(self) -> float:
        return float(self.base.time_end)

    @property
    def base_sha256(self) -> str:
        return self.base.sha256

    @property
    def supported_base_sha256(self) -> str:
        return self.base.base.sha256

    @property
    def axial_identity_half_height(self) -> float:
        return float(self.base.base.taper.axial_half_height) * float(
            np.sqrt(self.axial_plateau_q)
        )

    @property
    def taper(self) -> AxisymmetricPhysicalTaper:
        return replace(
            self.base.base.taper,
            axial_plateau_q=float(self.axial_plateau_q),
        )

    def snapshot(self, time: float) -> Eq45SupportedVelocityCandidate:
        """Return one scalar-time blend snapshot with the axial taper replaced."""
        values = np.asarray(time, dtype=float)
        if values.ndim != 0:
            raise ValueError("snapshot time must be a scalar")
        supported_snapshot = self.base.snapshot(float(values))
        return Eq45SupportedVelocityCandidate(
            parent=supported_snapshot.parent,
            taper=self.taper,
        )

    def _validated_inputs(self, points, time) -> tuple[np.ndarray, np.ndarray]:
        points_arr = np.asarray(points, dtype=float)
        if points_arr.ndim == 0 or points_arr.shape[-1] != 3:
            raise ValueError("points must have shape (..., 3)")
        if not np.all(np.isfinite(points_arr)):
            raise ValueError("points must be finite")
        time_arr = np.asarray(time, dtype=float)
        if not np.all(np.isfinite(time_arr)):
            raise ValueError("time must be finite")
        try:
            time_arr = np.broadcast_to(time_arr, points_arr.shape[:-1])
        except ValueError as exc:
            raise ValueError("time must broadcast to points.shape[:-1]") from exc
        self.base.tau(time_arr)
        return points_arr, time_arr

    def velocity(self, points, time) -> np.ndarray:
        """Evaluate Cartesian ``[...,3] = [u,v,w]`` through the retapered snapshots."""
        points_arr, time_arr = self._validated_inputs(points, time)
        flat_points = points_arr.reshape((-1, 3))
        flat_times = np.asarray(time_arr, dtype=float).reshape(-1)
        flat_out = np.empty_like(flat_points)
        for scalar_time in np.unique(flat_times):
            mask = flat_times == scalar_time
            flat_out[mask] = self.snapshot(float(scalar_time)).at_points(
                flat_points[mask], float(scalar_time)
            )
        output = flat_out.reshape(points_arr.shape)
        if not np.all(np.isfinite(output)):
            raise RuntimeError("axial-taper Eq45 velocity became nonfinite")
        return output

    def at_points(self, points, time) -> np.ndarray:
        return self.velocity(points, time)

    __call__ = velocity

    def velocity_xyz(self, x, y, z, time) -> np.ndarray:
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
        """Return regular samples in ``(time,x,y,z,component)`` order."""
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
        return {
            "schema": SCHEMA,
            "task_id": TASK_ID,
            "classification": dict(_CLASSIFICATION),
            "base_sha256": self.base.sha256,
            "base_blend_candidate": self.base.to_dict(),
            "axial_taper_control": {
                "parameter": "AxisymmetricPhysicalTaper.axial_plateau_q",
                "value": float(self.axial_plateau_q),
                "bounds": list(AXIAL_PLATEAU_Q_BOUNDS),
                "base_value": BASE_AXIAL_PLATEAU_Q,
                "value_selected": False,
                "axial_identity_half_height": self.axial_identity_half_height,
                "outer_axial_half_height": float(self.taper.axial_half_height),
                "resolved_taper": self.taper.to_dict(),
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @classmethod
    def from_dict(
        cls, data: Mapping[str, Any]
    ) -> "Eq45SupportedPhi10BlendAxialTaperCandidate":
        if not isinstance(data, Mapping):
            raise ValueError("axial taper candidate payload must be an object")
        if data.get("schema") != SCHEMA:
            raise ValueError(f"schema must be {SCHEMA!r}")
        if data.get("task_id") != TASK_ID:
            raise ValueError("axial taper task identity was modified")
        if data.get("classification") != _CLASSIFICATION:
            raise ValueError("candidate source classification was modified")
        if data.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("candidate truth-boundary metadata was modified")
        if "base_blend_candidate" not in data or "axial_taper_control" not in data:
            raise ValueError("axial taper candidate payload is missing a required component")

        base = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate.from_dict(
            data["base_blend_candidate"]
        )
        if data.get("base_sha256") != base.sha256:
            raise ValueError("base_sha256 does not match the embedded blend candidate")

        control = data["axial_taper_control"]
        if not isinstance(control, Mapping):
            raise ValueError("axial_taper_control must be an object")
        if control.get("parameter") != "AxisymmetricPhysicalTaper.axial_plateau_q":
            raise ValueError("axial taper parameter identity was modified")
        if control.get("bounds") != list(AXIAL_PLATEAU_Q_BOUNDS):
            raise ValueError("axial taper bounds were modified")
        if control.get("base_value") != BASE_AXIAL_PLATEAU_Q:
            raise ValueError("axial taper base value was modified")
        if control.get("value_selected") is not False:
            raise ValueError("axial taper selection provenance was modified")
        if "value" not in control:
            raise ValueError("axial taper value is missing")

        candidate = cls(base=base, axial_plateau_q=float(control["value"]))
        expected_control = candidate.to_dict()["axial_taper_control"]
        if dict(control) != expected_control:
            raise ValueError("serialized axial taper control does not match resolved candidate")
        return candidate

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(
        cls, path: str | Path
    ) -> "Eq45SupportedPhi10BlendAxialTaperCandidate":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
