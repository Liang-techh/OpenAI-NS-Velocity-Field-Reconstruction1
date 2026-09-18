"""Fail-closed segmented assembly for executable Kokuno leading-profile stages.

This module is the first repository object that places the independently reconstructed
Kokuno leading-profile pieces behind one candidate-facing velocity interface without
inventing the still-missing radial stages.  It composes

* the source-aligned nonlinear core/reference continuation through its published
  short logarithmic continuation, and
* the source I2 reserved-patch power profile with the corrected heat discrepancy
  applied by ``KokunoHighPrecisionHierarchicalHeatRepair``.

The I2 continuous-moment backend used here is the Agent-1 high-precision tanh--sinh
construction independently audited on repository PR #339 (exact audit head
``abf4152b981fba12a493f1fc9db0495661995505``) with a different arbitrary-
precision quadrature.  That audit is repository evidence, not a statement that the
field is OpenAI paper-exact.

Crucially, this object does *not* bridge the enormous unreconstructed radial gap
between the reference continuation and I2, and it does not pretend that the terminal
heat factor is already a source-matched global velocity.  Requests in unsupported
regions fail closed.  This makes missing source stages machine-visible while giving
the existing candidate evaluator a single ``velocity(x,y,z,t)->[...,3]`` entry point
for every reconstructed segment.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_high_precision_heat_repair import KokunoHighPrecisionHierarchicalHeatRepair
from .kokuno_reference_continuation import KokunoReferenceContinuationCandidate
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-segmented-leading-assembly-v1"
AUDIT_PR = 339
AUDIT_HEAD = "abf4152b981fba12a493f1fc9db0495661995505"

_REFERENCE_STAGE = "reference_continuation"
_I2_STAGE = "audited_repaired_I2"
_UNSUPPORTED_STAGE = "unreconstructed"

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "single_velocity_router_executable": True,
    "reference_continuation_stage_executable": True,
    "audited_high_precision_I2_repair_integrated": True,
    "unsupported_radial_gaps_fail_closed": True,
    "autonomous_gap_bridge_inserted": False,
    "heat_exterior_formula_available_elsewhere": True,
    "source_matched_heat_exterior_bound_here": False,
    "unified_profile_derivative_contract_completed": False,
    "global_leading_profile_reconstructed": False,
    "complete_kokuno_composite_velocity": False,
    "formal_full_domain_pde_gate_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


class UnreconstructedKokunoStageError(ValueError):
    """Raised when a requested point lies in a not-yet-reconstructed source stage."""


@dataclass(frozen=True)
class KokunoSegmentedLeadingAssembly:
    """Route reconstructed Kokuno leading-profile segments without filling gaps.

    The default component parameters are the autonomous, source-admissible choices
    already pinned by their component classes.  They are checked for a common
    ``h`` and a common Kokuno-native coordinate convention before any evaluation.
    """

    reference: KokunoReferenceContinuationCandidate = field(
        default_factory=KokunoReferenceContinuationCandidate
    )
    repaired_i2: KokunoHighPrecisionHierarchicalHeatRepair = field(
        default_factory=KokunoHighPrecisionHierarchicalHeatRepair
    )

    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.reference, KokunoReferenceContinuationCandidate):
            raise TypeError("reference must be a KokunoReferenceContinuationCandidate")
        if not isinstance(self.repaired_i2, KokunoHighPrecisionHierarchicalHeatRepair):
            raise TypeError(
                "repaired_i2 must be a KokunoHighPrecisionHierarchicalHeatRepair"
            )
        h_reference = float(self.reference.h)
        h_i2 = float(self.repaired_i2.outer_schedule.h)
        if not math.isclose(h_reference, h_i2, rel_tol=0.0, abs_tol=1.0e-15):
            raise ValueError("segmented assembly requires the same h in every stage")
        object.__setattr__(
            self, "_coordinates", KokunoNativeSimilarityCoordinates(h=h_reference)
        )

        low, high = self.i2_log_interval
        if not math.log(self.reference.X2) < low < high:
            raise ValueError("stage ordering is inconsistent: reference must precede I2")

    @property
    def h(self) -> float:
        return float(self.reference.h)

    @property
    def i2_log_interval(self) -> tuple[float, float]:
        low, high = self.repaired_i2.outer_schedule.reserved_log_intervals()["I2"]
        return float(low), float(high)

    @property
    def heat_log_X_K(self) -> float:
        return float(self.repaired_i2.target.terminal_schedule.log_X_K)

    def stage_for_similarity(self, X: Any) -> np.ndarray:
        """Classify finite nonnegative similarity radii into executable stages."""

        values = np.asarray(X, dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("X must contain only finite nonnegative values")
        stages = np.full(values.shape, _UNSUPPORTED_STAGE, dtype=object)
        stages[values <= self.reference.X2] = _REFERENCE_STAGE
        positive = values > 0.0
        log_X = np.full(values.shape, -math.inf, dtype=float)
        log_X[positive] = np.log(values[positive])
        low, high = self.i2_log_interval
        i2 = (log_X > low) & (log_X < high)
        stages[i2] = _I2_STAGE
        return stages

    def coverage_report(self) -> dict[str, Any]:
        """Return source-stage coverage in log-X coordinates without overflow."""

        low, high = self.i2_log_interval
        heat_log = self.heat_log_X_K
        return {
            "reference": {
                "stage": _REFERENCE_STAGE,
                "X_min": 0.0,
                "X_max": float(self.reference.X2),
                "log_X_max": math.log(float(self.reference.X2)),
                "executable_velocity": True,
            },
            "repaired_I2": {
                "stage": _I2_STAGE,
                "log_X_open_interval": [low, high],
                "log_X_star": float(self.repaired_i2.outer_schedule.log_X_star),
                "continuous_moment_backend": "arbitrary-precision tanh-sinh",
                "independent_audit_pr": AUDIT_PR,
                "independent_audit_head": AUDIT_HEAD,
                "executable_velocity": True,
            },
            "terminal_heat": {
                "log_X_K": heat_log,
                "log_X_b": float(
                    self.repaired_i2.target.terminal_schedule.log_X_b
                ),
                "matched_velocity_routed_here": False,
                "reason": (
                    "terminal heat formula exists as a component, but the intervening "
                    "outer stages and a source-matched global pressure/velocity splice "
                    "are not reconstructed in this assembly"
                ),
            },
            "unreconstructed_log_X_gaps": [
                [math.log(float(self.reference.X2)), low],
                [high, heat_log],
            ],
            "global_coverage_complete": False,
        }

    def _broadcast_inputs(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        arrays = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        if not all(np.all(np.isfinite(array)) for array in arrays):
            raise ValueError("physical inputs must contain only finite values")
        return arrays  # type: ignore[return-value]

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Evaluate reconstructed stages through one fail-closed velocity interface."""

        x_array, y_array, z_array, t_array = self._broadcast_inputs(x, y, z, t)
        coordinates = self._coordinates.evaluate(x_array, y_array, z_array, t_array)
        stages = self.stage_for_similarity(coordinates["X"])
        unsupported = stages == _UNSUPPORTED_STAGE
        if np.any(unsupported):
            log_values = np.full(stages.shape, -math.inf, dtype=float)
            X = np.asarray(coordinates["X"], dtype=float)
            positive = X > 0.0
            log_values[positive] = np.log(X[positive])
            bad = log_values[unsupported]
            raise UnreconstructedKokunoStageError(
                "requested point lies in an unreconstructed Kokuno radial stage; "
                f"unsupported_count={bad.size}, log_X_range="
                f"[{float(np.min(bad)):.17g},{float(np.max(bad)):.17g}]"
            )

        shape = x_array.shape
        x_flat = x_array.reshape(-1)
        y_flat = y_array.reshape(-1)
        z_flat = z_array.reshape(-1)
        t_flat = t_array.reshape(-1)
        stage_flat = stages.reshape(-1)
        output = np.empty((x_flat.size, 3), dtype=float)

        reference_indices = np.flatnonzero(stage_flat == _REFERENCE_STAGE)
        if reference_indices.size:
            values = self.reference.velocity(
                x_flat[reference_indices],
                y_flat[reference_indices],
                z_flat[reference_indices],
                t_flat[reference_indices],
            )
            output[reference_indices] = np.asarray(values, dtype=float).reshape(-1, 3)

        i2_indices = np.flatnonzero(stage_flat == _I2_STAGE)
        if i2_indices.size:
            values = self.repaired_i2.velocity(
                x_flat[i2_indices],
                y_flat[i2_indices],
                z_flat[i2_indices],
                t_flat[i2_indices],
            )
            output[i2_indices] = np.asarray(values, dtype=float).reshape(-1, 3)

        return output.reshape(shape + (3,))

    def pressure(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Evaluate pressure only where the reconstructed reference pressure exists.

        The I2 velocity repair does not yet carry the globally matched pressure
        integration constant.  Refusing it here prevents a local ``Pi_X=F^2``
        primitive from being mislabeled as the global pressure.
        """

        x_array, y_array, z_array, t_array = self._broadcast_inputs(x, y, z, t)
        coordinates = self._coordinates.evaluate(x_array, y_array, z_array, t_array)
        stages = self.stage_for_similarity(coordinates["X"])
        if np.any(stages != _REFERENCE_STAGE):
            raise UnreconstructedKokunoStageError(
                "global matched pressure is not reconstructed outside the reference stage"
            )
        return np.asarray(
            self.reference.pressure(x_array, y_array, z_array, t_array), dtype=float
        )

    def at_points(self, points_xyz: Any, t: Any) -> np.ndarray:
        points = np.asarray(points_xyz, dtype=float)
        if points.ndim == 0 or points.shape[-1] != 3:
            raise ValueError("points_xyz must have final dimension 3")
        if not np.all(np.isfinite(points)):
            raise ValueError("points_xyz must contain only finite values")
        return self.velocity(points[..., 0], points[..., 1], points[..., 2], t)

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "independent_repository_audit": {
                "pr": AUDIT_PR,
                "head": AUDIT_HEAD,
                "scope": "continuous high-precision three-moment I2 repair only",
            },
            "reference": self.reference.to_payload(),
            "repaired_i2": self.repaired_i2.to_payload(),
            "coverage": self.coverage_report(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(
            _canonical_json(self._unsigned_payload()).encode("utf-8")
        ).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSegmentedLeadingAssembly":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected segmented leading-assembly schema")
        reference_payload = payload.get("reference")
        repaired_payload = payload.get("repaired_i2")
        if not isinstance(reference_payload, dict) or not isinstance(repaired_payload, dict):
            raise ValueError("serialized leading-profile components are missing")
        obj = cls(
            reference=KokunoReferenceContinuationCandidate.from_payload(reference_payload),
            repaired_i2=KokunoHighPrecisionHierarchicalHeatRepair.from_payload(
                repaired_payload
            ),
        )
        if obj.to_payload() != payload:
            raise ValueError("segmented leading-assembly payload hash or content mismatch")
        return obj

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoSegmentedLeadingAssembly":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
