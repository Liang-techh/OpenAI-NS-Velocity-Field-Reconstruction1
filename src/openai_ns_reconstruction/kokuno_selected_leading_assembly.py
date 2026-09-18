"""Fail-closed selected Kokuno leading assembly through the I2 heat patch.

This module composes already executable leading-profile pieces without claiming
that the unpublished Kokuno shear loop was recovered.  It routes

* the reconstructed nonlinear reference continuation;
* the RF40b/RF40c outer/base schedule;
* the explicitly autonomous source-form finite-frequency insertion from #391;
* the eta-smooth PA.17 I1 repair family from #399; and
* the independently audited high-precision I2 heat repair from #337/#339.

A subtle point is essential for an executable velocity: the compact modulation
changes the prefix moment M=int U dX.  Even where the local modulation has
returned to E,U equal to the RF40 base, incompressibility retains a radial
memory through M/X until the I1 repair cancels it.  Merely adding the local I1
velocity correction would therefore be wrong.  This assembly carries that
prefix-M contribution through the post-modulation gap, through I1, and (at the
small interpolation residual left by the autonomous smooth family) through I2.

The still-unreconstructed inner-to-outer join and the later I3/I4 overlays are
left fail-closed.  Pressure is exposed only on the already matched reference
stage.  This is an executable provenance-labelled candidate segment, not a
paper-exact field and not independent Navier--Stokes validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_eta_smooth_shear_repair_family import KokunoEtaSmoothShearRepairFamily
from .kokuno_high_precision_heat_repair import KokunoHighPrecisionHierarchicalHeatRepair
from .kokuno_outer_base_schedule import KokunoOuterBaseSchedule
from .kokuno_reference_continuation import KokunoReferenceContinuationCandidate
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-selected-leading-assembly-v1"
I2_AUDIT_PR = 339
I2_AUDIT_HEAD = "abf4152b981fba12a493f1fc9db0495661995505"

_REFERENCE = "reference_continuation"
_OUTER_BASE = "rf40_outer_base"
_MODULATION = "selected_autonomous_modulation"
_POST_MOD_MEMORY = "post_modulation_prefix_memory"
_I1 = "selected_eta_smooth_I1_repair"
_POST_I1_MEMORY = "post_I1_residual_memory"
_I2 = "audited_repaired_I2"
_POST_I2_BASE = "post_I2_base_before_I3"
_UNSUPPORTED = "unreconstructed"

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "single_velocity_router_executable": True,
    "reference_continuation_stage_executable": True,
    "rf40_outer_base_routed": True,
    "selected_autonomous_modulation_routed": True,
    "modulation_prefix_M_memory_carried": True,
    "eta_smooth_I1_repair_routed": True,
    "audited_high_precision_I2_repair_routed": True,
    "source_hidden_loop_parameters_recovered": False,
    "source_admissible_loop_reconstructed": False,
    "actual_source_admissible_loop_discrepancy_supplied": False,
    "I1_cone_repair_applied_to_actual_source_modulation": False,
    "I3_positive_order_moment_correction_applied": False,
    "I4_mean_correction_applied": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
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


def _finite_array(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    return array


class UnreconstructedSelectedKokunoStageError(ValueError):
    """Raised when the selected assembly reaches a still-missing source stage."""


@dataclass(frozen=True)
class KokunoSelectedLeadingAssembly:
    """Compose the selected autonomous leading route without filling source gaps."""

    reference: KokunoReferenceContinuationCandidate = field(
        default_factory=KokunoReferenceContinuationCandidate
    )
    i1_family: KokunoEtaSmoothShearRepairFamily = field(
        default_factory=KokunoEtaSmoothShearRepairFamily
    )
    repaired_i2: KokunoHighPrecisionHierarchicalHeatRepair = field(
        default_factory=KokunoHighPrecisionHierarchicalHeatRepair
    )

    _outer_base: KokunoOuterBaseSchedule = field(init=False, repr=False, compare=False)
    _coordinates: KokunoNativeSimilarityCoordinates = field(init=False, repr=False, compare=False)
    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)
    _modulation_M_unit: float = field(init=False, repr=False, compare=False)
    _repair_M_row: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.reference, KokunoReferenceContinuationCandidate):
            raise TypeError("reference must be a KokunoReferenceContinuationCandidate")
        if not isinstance(self.i1_family, KokunoEtaSmoothShearRepairFamily):
            raise TypeError("i1_family must be a KokunoEtaSmoothShearRepairFamily")
        if not isinstance(self.repaired_i2, KokunoHighPrecisionHierarchicalHeatRepair):
            raise TypeError("repaired_i2 must be a KokunoHighPrecisionHierarchicalHeatRepair")

        schedule = self.i1_family.modulation.outer_schedule
        if schedule.to_payload() != self.repaired_i2.outer_schedule.to_payload():
            raise ValueError("I1 and I2 components must use the same outer schedule")
        if not math.isclose(float(self.reference.h), float(schedule.h), rel_tol=0.0, abs_tol=1.0e-15):
            raise ValueError("reference and outer stages must use the same h")

        outer_base = KokunoOuterBaseSchedule(outer_schedule=schedule)
        if not math.log(float(self.reference.X2)) < outer_base.log_X_R:
            raise ValueError("reference continuation must precede the RF40 outer schedule")

        modulation = self.i1_family.modulation
        mod_left, mod_right = modulation.support_log_interval
        i1_low, i1_high = schedule.reserved_log_intervals()["I1"]
        i2_low, i2_high = schedule.reserved_log_intervals()["I2"]
        i3_low, _ = schedule.reserved_log_intervals()["I3"]
        if not (
            outer_base.log_X_R < mod_left < mod_right < i1_low < i1_high
            <= i2_low < i2_high < i3_low <= outer_base.log_X_end
        ):
            raise ValueError("selected outer stage ordering is inconsistent")

        nodes, weights = leggauss(int(modulation.quadrature_order))
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)

        # For the selected insertion, normalized Delta M is exactly f(eta)
        # times its eta=0 value.  Cache that unit value once so candidate
        # evaluation does not replay the full five-moment quadrature.
        modulation_M_unit = float(
            np.asarray(modulation.normalized_discrepancy(0.0), dtype=float)[0]
        )
        if not math.isfinite(modulation_M_unit) or modulation_M_unit == 0.0:
            raise RuntimeError("selected modulation must have a finite nonzero M discrepancy")

        # PA.17's M row is exactly linear and eta-independent.  Cache its row
        # so the residual prefix memory and its analytic eta derivative can be
        # propagated without re-solving or finite differences.
        repair = self.i1_family.repair
        repair_M_row = np.asarray(
            repair.coefficient_jacobian(np.zeros(5, dtype=float), f_eta=1.0)[0],
            dtype=float,
        )
        repair_M_row.setflags(write=False)

        object.__setattr__(self, "_outer_base", outer_base)
        object.__setattr__(self, "_coordinates", KokunoNativeSimilarityCoordinates(h=float(schedule.h)))
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)
        object.__setattr__(self, "_modulation_M_unit", modulation_M_unit)
        object.__setattr__(self, "_repair_M_row", repair_M_row)

    @property
    def h(self) -> float:
        return float(self.reference.h)

    @property
    def outer_base(self) -> KokunoOuterBaseSchedule:
        return self._outer_base

    @property
    def outer_schedule(self):
        return self.i1_family.modulation.outer_schedule

    def _intervals(self) -> dict[str, tuple[float, float]]:
        reserved = self.outer_schedule.reserved_log_intervals()
        return {
            "modulation": tuple(float(v) for v in self.i1_family.modulation.support_log_interval),
            "I1": tuple(float(v) for v in reserved["I1"]),
            "I2": tuple(float(v) for v in reserved["I2"]),
            "I3": tuple(float(v) for v in reserved["I3"]),
            "I4": tuple(float(v) for v in reserved["I4"]),
        }

    def stage_for_similarity(self, X: Any) -> np.ndarray:
        values = np.asarray(X, dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values < 0.0):
            raise ValueError("X must contain only finite nonnegative values")
        stages = np.full(values.shape, _UNSUPPORTED, dtype=object)
        stages[values <= float(self.reference.X2)] = _REFERENCE

        positive = values > 0.0
        log_X = np.full(values.shape, -math.inf, dtype=float)
        log_X[positive] = np.log(values[positive])
        mod_left, mod_right = self._intervals()["modulation"]
        i1_low, i1_high = self._intervals()["I1"]
        i2_low, i2_high = self._intervals()["I2"]
        i3_low, _ = self._intervals()["I3"]

        mask = (log_X >= self.outer_base.log_X_R) & (log_X < mod_left)
        stages[mask] = _OUTER_BASE
        mask = (log_X >= mod_left) & (log_X <= mod_right)
        stages[mask] = _MODULATION
        mask = (log_X > mod_right) & (log_X <= i1_low)
        stages[mask] = _POST_MOD_MEMORY
        mask = (log_X > i1_low) & (log_X < i1_high)
        stages[mask] = _I1
        mask = (log_X >= i1_high) & (log_X <= i2_low)
        stages[mask] = _POST_I1_MEMORY
        mask = (log_X > i2_low) & (log_X < i2_high)
        stages[mask] = _I2
        mask = (log_X >= i2_high) & (log_X <= i3_low)
        stages[mask] = _POST_I2_BASE
        return stages

    def coverage_report(self) -> dict[str, Any]:
        intervals = self._intervals()
        return {
            "reference": {
                "X_max": float(self.reference.X2),
                "log_X_max": math.log(float(self.reference.X2)),
                "velocity": True,
                "matched_pressure": True,
            },
            "unreconstructed_inner_outer_join_log_X": [
                math.log(float(self.reference.X2)),
                self.outer_base.log_X_R,
            ],
            "selected_outer": {
                "log_X_start": self.outer_base.log_X_R,
                "log_X_stop_before_I3": intervals["I3"][0],
                "RF40_base": True,
                "selected_autonomous_modulation_log_X": list(intervals["modulation"]),
                "I1_eta_smooth_repair_log_X": list(intervals["I1"]),
                "I2_audited_heat_repair_log_X": list(intervals["I2"]),
                "modulation_prefix_M_memory_carried": True,
                "I2_independent_audit_pr": I2_AUDIT_PR,
                "I2_independent_audit_head": I2_AUDIT_HEAD,
            },
            "unreconstructed_later_overlays": {
                "I3": list(intervals["I3"]),
                "I4": list(intervals["I4"]),
                "terminal_heat_log_X_K": float(
                    self.repaired_i2.target.terminal_schedule.log_X_K
                ),
            },
            "global_coverage_complete": False,
            "global_pressure_matched": False,
        }

    def _modulation_M_normalized(self, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        f = np.asarray(self.outer_schedule.source_f(eta_array), dtype=float)
        f_eta = np.asarray(self.outer_schedule.source_f_eta(eta_array), dtype=float)
        return self._modulation_M_unit * f, self._modulation_M_unit * f_eta

    def _post_i1_M_normalized(self, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        eta_array = _finite_array(eta, "eta")
        mod_M, mod_M_eta = self._modulation_M_normalized(eta_array)
        coefficients = np.asarray(self.i1_family.coefficients(eta_array), dtype=float)
        coefficient_eta = np.asarray(self.i1_family.coefficient_eta(eta_array), dtype=float)
        repair_M = np.tensordot(coefficients, self._repair_M_row, axes=([-1], [0]))
        repair_M_eta = np.tensordot(coefficient_eta, self._repair_M_row, axes=([-1], [0]))
        return mod_M + repair_M, mod_M_eta + repair_M_eta

    def moment_memory_report(self, eta: Any) -> dict[str, Any]:
        eta_array = _finite_array(eta, "eta")
        mod_M, mod_M_eta = self._modulation_M_normalized(eta_array)
        residual_M, residual_M_eta = self._post_i1_M_normalized(eta_array)
        return {
            "eta": eta_array.tolist(),
            "normalized_modulation_M": np.asarray(mod_M).tolist(),
            "normalized_modulation_M_eta": np.asarray(mod_M_eta).tolist(),
            "normalized_post_I1_M_residual": np.asarray(residual_M).tolist(),
            "normalized_post_I1_M_residual_eta": np.asarray(residual_M_eta).tolist(),
            "I1_family_closure_tolerance": float(self.i1_family.closure_tolerance),
            "post_I1_residual_is_interpolation_error_not_source_target": True,
        }

    def _moment_scale_over_x(self, log_X: float) -> float:
        exponent = (
            float(self.i1_family.repair.log_X_1)
            + float(self.i1_family.repair.log_e_1)
            - float(log_X)
        )
        if exponent < math.log(np.nextafter(0.0, 1.0)):
            return 0.0
        if exponent > math.log(np.finfo(float).max):
            raise OverflowError("prefix-moment scale/X is outside float64 range")
        return math.exp(exponent)

    def _full_modulation_memory_over_x(self, log_X: float, eta: float) -> tuple[float, float]:
        scale = self._moment_scale_over_x(log_X)
        mod_M, mod_M_eta = self._modulation_M_normalized(float(eta))
        return scale * float(np.asarray(mod_M)), scale * float(np.asarray(mod_M_eta))

    def _post_i1_memory_over_x(self, log_X: float, eta: float) -> tuple[float, float]:
        scale = self._moment_scale_over_x(log_X)
        residual_M, residual_M_eta = self._post_i1_M_normalized(float(eta))
        return scale * float(np.asarray(residual_M)), scale * float(np.asarray(residual_M_eta))

    def _partial_modulation_memory_over_x(self, log_X: float, eta: float) -> tuple[float, float]:
        left, right = self._intervals()["modulation"]
        upper = min(float(log_X), right)
        if upper <= left:
            return 0.0, 0.0
        modulation = self.i1_family.modulation
        panels = max(
            1,
            int(math.ceil((upper - left) * float(modulation.frequency) * float(modulation.panels_per_period))),
        )
        edges = np.linspace(left, upper, panels + 1)
        points: list[np.ndarray] = []
        weights: list[np.ndarray] = []
        for a, b in zip(edges[:-1], edges[1:]):
            half = 0.5 * (b - a)
            centre = 0.5 * (b + a)
            points.append(centre + half * self._nodes)
            weights.append(half * self._weights)
        y = np.concatenate(points)
        w = np.concatenate(weights)
        profile = modulation.profile_values_logX(y, np.full_like(y, float(eta)))
        delta_U = np.asarray(profile["U"], dtype=float)
        m_over_x = float(np.dot(w, delta_U * np.exp(y - float(log_X))))
        f = float(np.asarray(self.outer_schedule.source_f(float(eta))))
        f_eta = float(np.asarray(self.outer_schedule.source_f_eta(float(eta))))
        m_eta_over_x = m_over_x * f_eta / f
        return m_over_x, m_eta_over_x

    def _memory_v0(self, eta: float, m_over_x: float, m_eta_over_x: float, delta_U: float = 0.0) -> float:
        D = 0.5 - self.h
        d = 1.0 - eta * eta
        L = 1.0 - 2.0 * self.h * eta * eta
        return (
            2.0 * eta * delta_U
            - 2.0 * D * eta * m_over_x
            - d * m_eta_over_x
        ) / L

    def _velocity_from_profiles(
        self, x: float, y: float, q: float, F: float, U: float, v0: float
    ) -> np.ndarray:
        swirl = q ** (-1.0 - self.h) * F
        radial = v0 / (2.0 * q)
        axial = q ** (-0.5 - self.h) * U
        return np.asarray(
            [radial * x - swirl * y, radial * y + swirl * x, axial],
            dtype=float,
        )

    def _outer_velocity_scalar(self, x: float, y: float, z: float, t: float, stage: str) -> np.ndarray:
        coordinates = self._coordinates.evaluate(x, y, z, t)
        X = float(np.asarray(coordinates["X"]))
        eta = float(np.asarray(coordinates["eta"]))
        q = float(np.asarray(coordinates["q"]))
        log_X = math.log(X)

        if stage == _OUTER_BASE:
            return np.asarray(self.outer_base.velocity(x, y, z, t), dtype=float).reshape(3)

        if stage == _MODULATION:
            base = self.outer_base.profile_values_logX(log_X, eta)
            mod = self.i1_family.modulation.profile_values_logX(log_X, eta)
            m_over_x, m_eta_over_x = self._partial_modulation_memory_over_x(log_X, eta)
            delta_U = float(np.asarray(mod["U"])) - float(np.asarray(base["U"]))
            v0 = float(np.asarray(base["v0"])) + self._memory_v0(
                eta, m_over_x, m_eta_over_x, delta_U=delta_U
            )
            return self._velocity_from_profiles(
                x,
                y,
                q,
                float(np.asarray(mod["F"])),
                float(np.asarray(mod["U"])),
                v0,
            )

        if stage == _POST_MOD_MEMORY:
            base = self.outer_base.profile_values_logX(log_X, eta)
            m_over_x, m_eta_over_x = self._full_modulation_memory_over_x(log_X, eta)
            v0 = float(np.asarray(base["v0"])) + self._memory_v0(
                eta, m_over_x, m_eta_over_x
            )
            return self._velocity_from_profiles(
                x,
                y,
                q,
                float(np.asarray(base["F"])),
                float(np.asarray(base["U"])),
                v0,
            )

        if stage == _I1:
            base = self.outer_base.profile_values_logX(log_X, eta)
            repair = self.i1_family.profile_correction_logX(log_X, eta)
            m_over_x, m_eta_over_x = self._full_modulation_memory_over_x(log_X, eta)
            v0 = (
                float(np.asarray(base["v0"]))
                + float(np.asarray(repair["delta_v0"]))
                + self._memory_v0(eta, m_over_x, m_eta_over_x)
            )
            return self._velocity_from_profiles(
                x,
                y,
                q,
                float(np.asarray(base["F"])) + float(np.asarray(repair["delta_F"])),
                float(np.asarray(base["U"])) + float(np.asarray(repair["delta_U"])),
                v0,
            )

        if stage in (_POST_I1_MEMORY, _POST_I2_BASE):
            base = self.outer_base.profile_values_logX(log_X, eta)
            m_over_x, m_eta_over_x = self._post_i1_memory_over_x(log_X, eta)
            v0 = float(np.asarray(base["v0"])) + self._memory_v0(
                eta, m_over_x, m_eta_over_x
            )
            return self._velocity_from_profiles(
                x,
                y,
                q,
                float(np.asarray(base["F"])),
                float(np.asarray(base["U"])),
                v0,
            )

        if stage == _I2:
            # Isolate the audited Decimal heat correction from its local
            # pure-swirl patch base, then place that correction on the full
            # RF40 base so the latter's incompressibility-derived v0 is not
            # silently discarded by the local patch helper.
            local = self.repaired_i2.velocity_decimal(x, y, z, t)
            patch = np.asarray(
                self.outer_schedule.patch_velocity(x, y, z, t), dtype=float
            ).reshape(3)
            base = np.asarray(self.outer_base.velocity(x, y, z, t), dtype=float).reshape(3)
            with_delta = np.asarray(
                [
                    float(Decimal(repr(float(base[index]))) + (local[index] - Decimal(repr(float(patch[index])))))
                    for index in range(3)
                ],
                dtype=float,
            )
            m_over_x, m_eta_over_x = self._post_i1_memory_over_x(log_X, eta)
            delta_v0 = self._memory_v0(eta, m_over_x, m_eta_over_x)
            radial = delta_v0 / (2.0 * q)
            with_delta[0] += radial * x
            with_delta[1] += radial * y
            return with_delta

        raise UnreconstructedSelectedKokunoStageError(
            f"selected outer stage {stage!r} is not executable"
        )

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        arrays = np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            _finite_array(z, "z"),
            _finite_array(t, "t"),
        )
        coordinates = self._coordinates.evaluate(*arrays)
        stages = self.stage_for_similarity(coordinates["X"])
        unsupported = stages == _UNSUPPORTED
        if np.any(unsupported):
            X = np.asarray(coordinates["X"], dtype=float)
            log_values = np.full(X.shape, -math.inf, dtype=float)
            positive = X > 0.0
            log_values[positive] = np.log(X[positive])
            bad = log_values[unsupported]
            raise UnreconstructedSelectedKokunoStageError(
                "requested point lies in an unreconstructed selected Kokuno stage; "
                f"unsupported_count={bad.size}, log_X_range="
                f"[{float(np.min(bad)):.17g},{float(np.max(bad)):.17g}]"
            )

        shape = arrays[0].shape
        flat = [array.reshape(-1) for array in arrays]
        stage_flat = stages.reshape(-1)
        output = np.empty((flat[0].size, 3), dtype=float)
        for index, stage in enumerate(stage_flat):
            if stage == _REFERENCE:
                value = self.reference.velocity(
                    float(flat[0][index]),
                    float(flat[1][index]),
                    float(flat[2][index]),
                    float(flat[3][index]),
                )
                output[index] = np.asarray(value, dtype=float).reshape(3)
            else:
                output[index] = self._outer_velocity_scalar(
                    float(flat[0][index]),
                    float(flat[1][index]),
                    float(flat[2][index]),
                    float(flat[3][index]),
                    str(stage),
                )
        return output.reshape(shape + (3,))

    def pressure(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        arrays = np.broadcast_arrays(
            _finite_array(x, "x"),
            _finite_array(y, "y"),
            _finite_array(z, "z"),
            _finite_array(t, "t"),
        )
        coordinates = self._coordinates.evaluate(*arrays)
        stages = self.stage_for_similarity(coordinates["X"])
        if np.any(stages != _REFERENCE):
            raise UnreconstructedSelectedKokunoStageError(
                "global matched pressure is not reconstructed on selected outer stages"
            )
        return np.asarray(self.reference.pressure(*arrays), dtype=float)

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
            "I2_independent_repository_audit": {
                "pr": I2_AUDIT_PR,
                "head": I2_AUDIT_HEAD,
                "scope": "continuous high-precision three-moment I2 repair only",
            },
            "reference": self.reference.to_payload(),
            "i1_family": self.i1_family.to_payload(),
            "repaired_i2": self.repaired_i2.to_payload(),
            "assembly_numerics": {
                "partial_modulation_prefix_M": (
                    "composite Gauss-Legendre using the selected modulation's declared order/panels-per-period"
                ),
                "full_modulation_M": "exact f(eta) factor times cached eta=0 normalized discrepancy",
                "I1_M_row": "cached exact-linear PA.17 M row with analytic coefficient_eta",
                "post_I1_memory": "carry residual instead of silently resetting it to zero",
            },
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoSelectedLeadingAssembly":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno selected-leading assembly schema")
        reference_payload = payload.get("reference")
        i1_payload = payload.get("i1_family")
        i2_payload = payload.get("repaired_i2")
        if not all(isinstance(value, dict) for value in (reference_payload, i1_payload, i2_payload)):
            raise ValueError("serialized selected-leading components are missing")
        obj = cls(
            reference=KokunoReferenceContinuationCandidate.from_payload(reference_payload),
            i1_family=KokunoEtaSmoothShearRepairFamily.from_payload(i1_payload),
            repaired_i2=KokunoHighPrecisionHierarchicalHeatRepair.from_payload(i2_payload),
        )
        if obj.to_payload() != payload:
            raise ValueError("selected-leading assembly payload hash or content mismatch")
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
    def load_json(cls, path: str | Path) -> "KokunoSelectedLeadingAssembly":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
