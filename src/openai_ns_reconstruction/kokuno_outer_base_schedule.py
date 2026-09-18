"""Executable Kokuno outer/base radial schedule before cone and repair overlays.

This module implements the explicit source schedule recorded in the corrected
2026-09-09 Kokuno reconstruction (RF40b/RF40c).  It is deliberately the
*base* outer profile only.  Later cone modulation and the four reserved-patch
corrections are overlays on this backbone and are not silently inserted here.

With ``x=X/X_R`` and ``f=(1+eta^2)^(-1)``, the source temporary reference
branch is

    U = 4 eta,
    E = P_* f x^(1/10),                 x <= 1.

At ``X_R`` the logarithmic coordinate is reset at each stage.  Writing
``l = d_y log(sqrt(2X) E)``, RF40b prescribes consecutively

    length 1:     U=4 eta,       l=(3/5)(1-sigma(y));
    length T_d:   U=k(y) eta,    l=0,
                  k(y)=4[1-sigma(log(1+y)/M_d)];
    length 1:     U=0,           l=-lambda sigma(y).

The next RF40c stage has length ``T_w`` and keeps ``U=0,l=-lambda``.
The endpoint amplitudes are exactly the relations already implemented by
``KokunoOuterReservedPatchSchedule``:

    P_1=P_* exp(-1/5),
    e_w=P_1 exp(-T_d/2) exp(-1/2-lambda/2),
    X_w=X_R exp(T_d+2).

The source flat step is reused from the reference-continuation module.  Only
its primitive on (0,1) needs numerical quadrature; this implementation uses a
fixed Gauss--Legendre rule and pins the exact endpoint integral 1/2 using the
source symmetry sigma(1-s)=1-sigma(s).

The object exposes the low-dimensional profile and analytic first derivatives,
the incompressibility-derived ``v0``, and a native-coordinate velocity API.
It is useful as the backbone on which the still-missing cone / I1 / I2 / I3 /
I4 overlays can act.  It is NOT routed as a final global Kokuno candidate and
is not independent Navier--Stokes validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_reference_continuation import source_smooth_step
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-outer-base-schedule-v1"

_REFERENCE = "temporary_reference"
_FIRST_TURN = "first_l_turn"
_AXIAL_SHUTDOWN = "axial_shutdown"
_LAMBDA_TURN = "lambda_turn"
_POWER = "power_law"
_OUTSIDE = "outside_source_base_schedule"

_SOURCE_FORMULAS = {
    "temporary_reference": (
        "x=X/X_R; y=log x; U=4 eta; E=P_* (1+eta^2)^(-1) exp(y/10) for y<=0"
    ),
    "stage_reset": "the logarithmic coordinate y is reset to zero at each stage left endpoint",
    "first_turn": "length 1; U=4 eta; l=(3/5)(1-sigma(y))",
    "axial_shutdown": (
        "length T_d; U=k(y) eta; l=0; "
        "k(y)=4[1-sigma(log(1+y)/M_d)]"
    ),
    "lambda_turn": "length 1; U=0; l=-lambda*sigma(y)",
    "power_law": "length T_w; U=0; l=-lambda",
    "log_E_evolution": "d_y log E=l-1/2",
    "f": "f(eta)=(1+eta^2)^(-1)",
    "endpoint_amplitudes": (
        "P1=P_* exp(-1/5); e_w=P1 exp(-T_d/2) exp(-1/2-lambda/2)"
    ),
    "incompressibility": (
        "V0=(2 eta X U-2D eta M-d M_eta)/L; v0=V0/X"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "outer_base_schedule_executable": True,
    "rf40b_stage_resets_executable": True,
    "rf40c_power_stage_executable": True,
    "analytic_first_profile_derivatives": True,
    "native_velocity_api_executable": True,
    "autonomous_quadrature_only": True,
    "autonomous_outer_parameter_values_inherited": True,
    "cone_modulation_completed": False,
    "I1_cone_repair_applied": False,
    "I2_heat_compensation_applied_here": False,
    "I3_positive_order_moment_correction_applied": False,
    "I4_mean_correction_applied": False,
    "global_pressure_matched": False,
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


def _source_step_derivative(s: Any) -> np.ndarray:
    """Analytic derivative of the source flat step on its open unit interval."""

    values = _finite_array(s, "s")
    out = np.zeros_like(values, dtype=float)
    mask = (values > 0.0) & (values < 1.0)
    if np.any(mask):
        local = values[mask]
        sigma = source_smooth_step(local)
        logit_prime = 2.0 / np.power(local, 3) + 2.0 / np.power(1.0 - local, 3)
        out[mask] = sigma * (1.0 - sigma) * logit_prime
    return out


@dataclass(frozen=True)
class KokunoOuterBaseSchedule:
    """Vectorized RF40b/RF40c base-profile evaluator.

    ``outer_schedule`` carries the source scale relations and the explicit
    autonomous admissible numerical choices already pinned by Agent 1.  The
    only new numerical choice is the fixed quadrature order used to integrate
    the source flat step and the stage-2 prefix moment.  It changes numerical
    evaluation accuracy, not the source profile definition.
    """

    outer_schedule: KokunoOuterReservedPatchSchedule = field(
        default_factory=KokunoOuterReservedPatchSchedule
    )
    quadrature_points: int = 64

    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)
    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be a KokunoOuterReservedPatchSchedule")
        if isinstance(self.quadrature_points, bool) or not isinstance(
            self.quadrature_points, (int, np.integer)
        ):
            raise TypeError("quadrature_points must be an integer")
        order = int(self.quadrature_points)
        if not 24 <= order <= 192:
            raise ValueError("quadrature_points must lie in [24,192]")
        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)
        object.__setattr__(self, "quadrature_points", order)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)
        object.__setattr__(
            self,
            "_coordinates",
            KokunoNativeSimilarityCoordinates(h=float(self.outer_schedule.h)),
        )

    @property
    def h(self) -> float:
        return float(self.outer_schedule.h)

    @property
    def log_X_R(self) -> float:
        return float(self.outer_schedule.log_X_R)

    @property
    def log_X_w(self) -> float:
        return float(self.outer_schedule.log_X_w)

    @property
    def log_X_end(self) -> float:
        return self.log_X_w + float(self.outer_schedule.T_w)

    @property
    def log_stage1_end(self) -> float:
        return self.log_X_R + 1.0

    @property
    def log_stage2_end(self) -> float:
        return self.log_X_R + 1.0 + float(self.outer_schedule.T_d)

    def source_interval_report(self) -> dict[str, Any]:
        return {
            "temporary_reference_log_X": [None, self.log_X_R],
            "first_l_turn_log_X": [self.log_X_R, self.log_stage1_end],
            "axial_shutdown_log_X": [self.log_stage1_end, self.log_stage2_end],
            "lambda_turn_log_X": [self.log_stage2_end, self.log_X_w],
            "power_law_log_X": [self.log_X_w, self.log_X_end],
            "reserved_log_intervals": {
                key: list(value)
                for key, value in self.outer_schedule.reserved_log_intervals().items()
            },
            "final_overlay_corrections_applied": False,
        }

    def stage_for_log_X(self, log_X: Any) -> np.ndarray:
        values = _finite_array(log_X, "log_X")
        out = np.full(values.shape, _OUTSIDE, dtype=object)
        out[values <= self.log_X_R] = _REFERENCE
        mask = (values > self.log_X_R) & (values <= self.log_stage1_end)
        out[mask] = _FIRST_TURN
        mask = (values > self.log_stage1_end) & (values <= self.log_stage2_end)
        out[mask] = _AXIAL_SHUTDOWN
        mask = (values > self.log_stage2_end) & (values <= self.log_X_w)
        out[mask] = _LAMBDA_TURN
        mask = (values > self.log_X_w) & (values <= self.log_X_end)
        out[mask] = _POWER
        return out

    def _step_primitive(self, y: np.ndarray) -> np.ndarray:
        """Return integral_0^y sigma(s) ds, using exact endpoint values."""

        values = np.asarray(y, dtype=float)
        out = np.empty_like(values)
        out[values <= 0.0] = 0.0
        out[values >= 1.0] = 0.5
        mask = (values > 0.0) & (values < 1.0)
        if np.any(mask):
            flat = values[mask].reshape(-1)
            local = np.empty_like(flat)
            for i, upper in enumerate(flat):
                points = 0.5 * upper * (self._nodes + 1.0)
                local[i] = 0.5 * upper * float(
                    np.dot(self._weights, source_smooth_step(points))
                )
            out[mask] = local.reshape(values[mask].shape)
        return out

    def _stage2_k_and_Dk(self, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        local = np.maximum(np.asarray(y, dtype=float), 0.0)
        r = np.log1p(local) / float(self.outer_schedule.M_d)
        sigma = source_smooth_step(r)
        dsigma = _source_step_derivative(r)
        k = 4.0 * (1.0 - sigma)
        Dk = -4.0 * dsigma / (
            float(self.outer_schedule.M_d) * (1.0 + local)
        )
        return k, Dk

    def _stage2_m_ratio(self, y: np.ndarray) -> np.ndarray:
        """Return B(X)/X where M=eta*B during the axial-shutdown stage."""

        values = np.asarray(y, dtype=float)
        flat = values.reshape(-1)
        out = np.empty_like(flat)
        for i, upper in enumerate(flat):
            if upper <= 0.0:
                out[i] = 4.0
                continue
            points = 0.5 * upper * (self._nodes + 1.0)
            k, _ = self._stage2_k_and_Dk(points)
            weighted = k * np.exp(points - upper)
            integral = 0.5 * upper * float(np.dot(self._weights, weighted))
            out[i] = 4.0 * math.exp(-upper) + integral
        return out.reshape(values.shape)

    def _base_state_from_log_X(
        self, log_X: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return log(E/f), k=U/eta, D_X k, l, and B/X."""

        log_values = np.asarray(log_X, dtype=float)
        stages = self.stage_for_log_X(log_values)
        if np.any(stages == _OUTSIDE):
            bad = log_values[stages == _OUTSIDE]
            raise ValueError(
                "outer base schedule ends after the RF40c T_w stage; "
                f"unsupported log_X range [{float(np.min(bad)):.17g},"
                f"{float(np.max(bad)):.17g}]"
            )

        log_E_over_f = np.empty_like(log_values)
        k = np.empty_like(log_values)
        Dk = np.empty_like(log_values)
        ell = np.empty_like(log_values)
        m_ratio = np.empty_like(log_values)

        mask = stages == _REFERENCE
        if np.any(mask):
            y = log_values[mask] - self.log_X_R
            log_E_over_f[mask] = float(self.outer_schedule.log_P_star) + 0.1 * y
            k[mask] = 4.0
            Dk[mask] = 0.0
            ell[mask] = 0.6
            m_ratio[mask] = 4.0

        mask = stages == _FIRST_TURN
        if np.any(mask):
            y = log_values[mask] - self.log_X_R
            sigma = source_smooth_step(y)
            primitive = self._step_primitive(y)
            log_E_over_f[mask] = (
                float(self.outer_schedule.log_P_star)
                + 0.1 * y
                - 0.6 * primitive
            )
            k[mask] = 4.0
            Dk[mask] = 0.0
            ell[mask] = 0.6 * (1.0 - sigma)
            m_ratio[mask] = 4.0

        mask = stages == _AXIAL_SHUTDOWN
        if np.any(mask):
            y = log_values[mask] - self.log_stage1_end
            local_k, local_Dk = self._stage2_k_and_Dk(y)
            log_E_over_f[mask] = float(self.outer_schedule.log_P1) - 0.5 * y
            k[mask] = local_k
            Dk[mask] = local_Dk
            ell[mask] = 0.0
            m_ratio[mask] = self._stage2_m_ratio(y)

        end_ratio = float(
            self._stage2_m_ratio(np.asarray(float(self.outer_schedule.T_d))).item()
        )

        mask = stages == _LAMBDA_TURN
        if np.any(mask):
            y = log_values[mask] - self.log_stage2_end
            sigma = source_smooth_step(y)
            primitive = self._step_primitive(y)
            log_E_over_f[mask] = (
                float(self.outer_schedule.log_P1)
                - 0.5 * float(self.outer_schedule.T_d)
                - 0.5 * y
                - float(self.outer_schedule.lambda_outer) * primitive
            )
            k[mask] = 0.0
            Dk[mask] = 0.0
            ell[mask] = -float(self.outer_schedule.lambda_outer) * sigma
            m_ratio[mask] = end_ratio * np.exp(-y)

        mask = stages == _POWER
        if np.any(mask):
            y = log_values[mask] - self.log_X_w
            exponent = -0.5 - float(self.outer_schedule.lambda_outer)
            log_E_over_f[mask] = float(self.outer_schedule.log_e_w) + exponent * y
            k[mask] = 0.0
            Dk[mask] = 0.0
            ell[mask] = -float(self.outer_schedule.lambda_outer)
            m_ratio[mask] = end_ratio * np.exp(-1.0 - y)

        return log_E_over_f, k, Dk, ell, m_ratio

    @staticmethod
    def _checked_exp(log_values: np.ndarray, name: str) -> np.ndarray:
        max_log = math.log(np.finfo(float).max)
        min_log = math.log(np.nextafter(0.0, 1.0))
        if np.any(log_values > max_log) or np.any(log_values < min_log):
            raise OverflowError(
                f"{name} is outside float range; use source_interval_report/log-scale data"
            )
        return np.exp(log_values)

    def profile_values_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate the source base profile using log-X as the radial input.

        Logarithmic radial derivatives ``D_X_*`` are returned directly so the
        schedule remains usable even when an autonomous source-admissible scale
        would place absolute X outside ordinary floating-point range.
        """

        log_array, eta_array = np.broadcast_arrays(
            _finite_array(log_X, "log_X"), _finite_array(eta, "eta")
        )
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")

        log_E_over_f, k, Dk, ell, m_ratio = self._base_state_from_log_X(log_array)
        f = self.outer_schedule.source_f(eta_array)
        f_eta = self.outer_schedule.source_f_eta(eta_array)
        log_f = np.log(f)
        log_E = log_E_over_f + log_f
        E = self._checked_exp(log_E, "E")
        log_F = log_E - 0.5 * (math.log(2.0) + log_array)
        F = self._checked_exp(log_F, "F")

        f_ratio = f_eta / f
        U = eta_array * k
        U_eta = k
        D_X_U = eta_array * Dk
        D_X_E = (ell - 0.5) * E
        E_eta = f_ratio * E
        D_X_F = (ell - 1.0) * F
        F_eta = f_ratio * F

        L = 1.0 - 2.0 * self.h * eta_array * eta_array
        v0 = 2.0 * eta_array * eta_array * k / L - m_ratio
        D_X_m_ratio = k - m_ratio
        D_X_v0 = 2.0 * eta_array * eta_array * Dk / L - D_X_m_ratio
        v0_eta = 4.0 * eta_array * k / (L * L)

        return {
            "stage": self.stage_for_log_X(log_array),
            "log_X": log_array,
            "ell": ell,
            "k": k,
            "m_ratio": m_ratio,
            "E": E,
            "D_X_E": D_X_E,
            "E_eta": E_eta,
            "F": F,
            "D_X_F": D_X_F,
            "F_eta": F_eta,
            "U": U,
            "D_X_U": D_X_U,
            "U_eta": U_eta,
            "v0": v0,
            "D_X_v0": D_X_v0,
            "v0_eta": v0_eta,
            "Pi_X": F * F,
        }

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate values plus ordinary X derivatives on materialized X>0."""

        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array <= 0.0):
            raise ValueError("X must be strictly positive")
        values = self.profile_values_logX(np.log(X_array), eta_array)
        values["X"] = X_array
        values["E_X"] = values["D_X_E"] / X_array
        values["F_X"] = values["D_X_F"] / X_array
        values["U_X"] = values["D_X_U"] / X_array
        values["v0_X"] = values["D_X_v0"] / X_array
        return values

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Evaluate the RF40 base velocity in Kokuno native coordinates.

        This is a source-stage velocity, not the final cone/correction-overlaid
        field.  Calls outside the explicit RF40c base schedule fail closed.
        """

        coordinates = self._coordinates.evaluate(x, y, z, t)
        profiles = self.profile_values(coordinates["X"], coordinates["eta"])
        x_array, y_array = np.broadcast_arrays(_finite_array(x, "x"), _finite_array(y, "y"))
        x_array = np.broadcast_to(x_array, np.shape(profiles["F"]))
        y_array = np.broadcast_to(y_array, np.shape(profiles["F"]))
        q = np.asarray(coordinates["q"], dtype=float)
        swirl_factor = np.power(q, -1.0 - self.h) * profiles["F"]
        radial_factor = profiles["v0"] / (2.0 * q)
        axial_factor = np.power(q, -0.5 - self.h) * profiles["U"]
        return np.stack(
            (
                radial_factor * x_array - swirl_factor * y_array,
                radial_factor * y_array + swirl_factor * x_array,
                axial_factor,
            ),
            axis=-1,
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
            "source_formulas": dict(_SOURCE_FORMULAS),
            "outer_schedule": self.outer_schedule.to_payload(),
            "quadrature_points": self.quadrature_points,
            "autonomous_numerics": {
                "flat_step_primitive": "fixed Gauss-Legendre; exact 0 and 1 endpoint primitives pinned",
                "stage2_prefix_moment": "fixed Gauss-Legendre in reset log-X coordinate",
            },
            "intervals": self.source_interval_report(),
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoOuterBaseSchedule":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno outer base-schedule schema")
        schedule_payload = payload.get("outer_schedule")
        if not isinstance(schedule_payload, dict):
            raise ValueError("outer_schedule is missing")
        obj = cls(
            outer_schedule=KokunoOuterReservedPatchSchedule.from_payload(schedule_payload),
            quadrature_points=int(payload.get("quadrature_points")),
        )
        if obj.to_payload() != payload:
            raise ValueError("outer base-schedule payload hash or content mismatch")
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
    def load_json(cls, path: str | Path) -> "KokunoOuterBaseSchedule":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
