"""Callable supported Eq45 candidate on the screened compact--quartic Phi10 blend.

Agent-7 PR #184 screened one scalar interpolation degree between the existing
balanced-quartic and slope-capped compact-C2 temporal schedules on the same
``Phi(1,0)`` spatial mode.  It deliberately selected no blend weight.  This
module only materializes that already-screened representation so downstream
visualization and independent validation can evaluate an explicitly chosen,
bounded weight through one stable ``velocity(x,y,z,t)->[u,v,w]`` object.

No image, PDE residual, force, pressure, or hidden OpenAI parameter chooses the
weight here.  The class is a visualization/research candidate only.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_c2_compact_temporal_capacity import (
    DEFAULT_EARLY_DELTA,
    compact_phi10_snapshot,
    compact_window_parameters,
)
from .constrained_eq45_supported_phi10_compact_quartic_blend_capacity import (
    blended_phi10_delta,
    blended_phi10_snapshot,
)
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)


SCHEMA = "eq45_supported_phi10_compact_quartic_blend_temporal_candidate_v1"
_MODE = (1, 0)
_BLEND_FORMULA = "(1-lambda)*delta_quartic(tau)+lambda*delta_compact(tau)"
_TAU_DEFINITION = "delivery_interval_endpoints_map_to_minus1_plus1"
_COMPACT_RULE = "earliest_c2_return_under_balanced_quartic_peak_slope_cap"
_ENDPOINT_FAMILIES = ("derivative_balanced_quartic", "slope_capped_c2_compact")

_CLASSIFICATION = {
    "velocity_interface": "user_requirement",
    "eq45_backbone": "public_source_fact",
    "base_supported_candidate": "parent_provenance",
    "phi10_mode": "autonomous_design",
    "compact_quartic_blend_family": "autonomous_design",
    "blend_weight": "autonomous_design",
    "early_delta": "autonomous_design",
    "openai_hidden_profiles": "pending_unknown",
}

_TRUTH_BOUNDARY = {
    "callable_serializable": True,
    "velocity_export_ready": True,
    "visualization_candidate_only": True,
    "canonical_velocity_changed": False,
    "screened_temporal_blend_materialized": True,
    "blend_weight_selected_by_this_module": False,
    "physical_support_connection_implemented": True,
    "physical_support_validated": False,
    "spatial_basis_grown": False,
    "force_or_pressure_fitted": False,
    "pde_objective_used_to_choose_blend": False,
    "public_image_fitted": False,
    "production_temporal_shape_promoted": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validate_weight(value: float) -> float:
    weight = float(value)
    if not np.isfinite(weight) or not (0.0 <= weight <= 1.0):
        raise ValueError("blend_weight must be finite and lie in [0,1]")
    return weight


@dataclass(frozen=True)
class Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate:
    """Support-connected Eq45 field with one explicit compact--quartic blend weight."""

    base: Eq45SupportedVelocityCandidate
    blend_weight: float
    early_delta: float = DEFAULT_EARLY_DELTA

    def __post_init__(self) -> None:
        if not isinstance(self.base, Eq45SupportedVelocityCandidate):
            raise TypeError("base must be Eq45SupportedVelocityCandidate")
        if _MODE not in self.profile_basis.mode_indices:
            raise ValueError("base profile basis must contain Phi(1,0)")
        object.__setattr__(self, "blend_weight", _validate_weight(self.blend_weight))
        if not np.isscalar(self.early_delta) or not np.isfinite(self.early_delta):
            raise ValueError("early_delta must be a finite scalar")
        if float(self.early_delta) == 0.0:
            raise ValueError("early_delta must be nonzero; keep the static base candidate otherwise")

        # Reuse the already-audited endpoint families instead of inventing a
        # second bound proof.  The compact trajectory is a monotone smootherstep
        # between early_delta and zero; the quartic object validates its complete
        # schedule.  Their convex coefficient blend therefore remains inside the
        # same coefficient interval for every delivery time.
        Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
            base=self.base, early_delta=float(self.early_delta)
        )
        compact_phi10_snapshot(
            self.base, self.time_start, early_delta=float(self.early_delta)
        )
        compact_phi10_snapshot(
            self.base, self.time_midpoint, early_delta=float(self.early_delta)
        )

    @property
    def profile_basis(self):
        return self.base.parent.profile_basis

    @property
    def mode_index(self) -> int:
        return self.profile_basis.mode_indices.index(_MODE)

    @property
    def midpoint_coefficient(self) -> float:
        return float(self.profile_basis.phi_coefficients[self.mode_index])

    @property
    def time_start(self) -> float:
        return float(self.base.time_start)

    @property
    def time_end(self) -> float:
        return float(self.base.time_end)

    @property
    def time_midpoint(self) -> float:
        return 0.5 * (self.time_start + self.time_end)

    @property
    def base_sha256(self) -> str:
        return self.base.sha256

    @property
    def compact_parameters(self) -> dict[str, float]:
        return compact_window_parameters(early_delta=float(self.early_delta))

    def tau(self, time):
        values = np.asarray(time, dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("time must be finite")
        if np.any((values < self.time_start) | (values > self.time_end)):
            raise ValueError(
                f"time must lie in declared interval [{self.time_start}, {self.time_end}]"
            )
        half_width = 0.5 * (self.time_end - self.time_start)
        return (values - self.time_midpoint) / half_width

    def coefficient_at_tau(self, tau):
        return self.midpoint_coefficient + blended_phi10_delta(
            tau,
            blend_weight=float(self.blend_weight),
            early_delta=float(self.early_delta),
        )

    def coefficient_at(self, time):
        return self.coefficient_at_tau(self.tau(time))

    def snapshot(self, time: float) -> Eq45SupportedVelocityCandidate:
        values = np.asarray(time, dtype=float)
        if values.ndim != 0:
            raise ValueError("snapshot time must be a scalar")
        return blended_phi10_snapshot(
            self.base,
            float(values),
            blend_weight=float(self.blend_weight),
            early_delta=float(self.early_delta),
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
        self.tau(time_arr)
        return points_arr, time_arr

    def velocity(self, points, time) -> np.ndarray:
        """Evaluate Cartesian ``[...,3] = [u,v,w]`` from the blended temporal field."""
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
            raise RuntimeError("compact-quartic blend Eq45 velocity became nonfinite")
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
        compact = self.compact_parameters
        return {
            "schema": SCHEMA,
            "classification": dict(_CLASSIFICATION),
            "base_sha256": self.base.sha256,
            "base_supported_candidate": self.base.to_dict(),
            "temporal_mode": {
                "family": "phi",
                "mode": list(_MODE),
                "early_delta": float(self.early_delta),
                "blend_weight": float(self.blend_weight),
                "blend_weight_selected": False,
                "endpoint_families": list(_ENDPOINT_FAMILIES),
                "blend_formula": _BLEND_FORMULA,
                "tau_definition": _TAU_DEFINITION,
                "compact_window_rule": _COMPACT_RULE,
                "compact_return_tau": float(compact["return_tau"]),
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @classmethod
    def from_dict(
        cls, data: Mapping[str, Any]
    ) -> "Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate":
        if not isinstance(data, Mapping):
            raise ValueError("compact-quartic blend candidate payload must be an object")
        if data.get("schema") != SCHEMA:
            raise ValueError(f"schema must be {SCHEMA!r}")
        if data.get("classification") != _CLASSIFICATION:
            raise ValueError("candidate source classification was modified")
        if data.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("candidate truth-boundary metadata was modified")
        if "base_supported_candidate" not in data or "temporal_mode" not in data:
            raise ValueError("candidate payload is missing a required component")

        base = Eq45SupportedVelocityCandidate.from_dict(data["base_supported_candidate"])
        if data.get("base_sha256") != base.sha256:
            raise ValueError("base_sha256 does not match the embedded supported candidate")
        temporal = data["temporal_mode"]
        if not isinstance(temporal, Mapping):
            raise ValueError("temporal_mode must be an object")
        if temporal.get("family") != "phi" or temporal.get("mode") != list(_MODE):
            raise ValueError("temporal mode must remain exactly Phi(1,0)")
        if temporal.get("endpoint_families") != list(_ENDPOINT_FAMILIES):
            raise ValueError("temporal endpoint families were modified")
        if temporal.get("blend_formula") != _BLEND_FORMULA:
            raise ValueError("temporal blend formula was modified")
        if temporal.get("tau_definition") != _TAU_DEFINITION:
            raise ValueError("temporal tau definition was modified")
        if temporal.get("compact_window_rule") != _COMPACT_RULE:
            raise ValueError("compact window rule was modified")
        if temporal.get("blend_weight_selected") is not False:
            raise ValueError("blend weight selection provenance was modified")
        if "early_delta" not in temporal or "blend_weight" not in temporal:
            raise ValueError("temporal_mode is missing a required coefficient")

        early_delta = float(temporal["early_delta"])
        expected_return = float(compact_window_parameters(early_delta=early_delta)["return_tau"])
        supplied_return = float(temporal.get("compact_return_tau", np.nan))
        if not np.isfinite(supplied_return) or supplied_return != expected_return:
            raise ValueError("serialized compact return does not match the deterministic rule")
        return cls(
            base=base,
            blend_weight=float(temporal["blend_weight"]),
            early_delta=early_delta,
        )

    def save_json(self, path: str | Path) -> None:
        Path(path).write_text(
            json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    @classmethod
    def load_json(
        cls, path: str | Path
    ) -> "Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate":
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))
