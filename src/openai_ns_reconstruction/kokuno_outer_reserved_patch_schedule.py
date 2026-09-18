"""Executable Kokuno outer reserved-patch scale binding.

The corrected 2026-09-09 reconstruction does not publish one hidden numerical
outer schedule, but it *does* publish the scale relations that generate the
four reserved power-law patches.  In the notation of the pinned reader,

    T_d = exp(M_d) + 10,
    T_w = 60 log(1/lambda),
    X_R = 110 (C P_*)^10,
    X_w = X_R exp(T_d + 2),

and the second reserved interval is

    I_2 = X_w exp((T_w-20, T_w-15)).

The power-law background on every reserved interval is exactly

    U = 0,
    E = c_patch f(eta) X^(-1/2-lambda),
    f(eta) = (1+eta^2)^(-1),

with

    P_1 = P_* exp(-1/5),
    e_w = P_1 exp(-T_d/2) exp(-1/2-lambda/2),
    c_patch = e_w X_w^(1/2+lambda).

The two one-half integrals above follow from the source flat step symmetry
sigma(1-s)=1-sigma(s).  This module keeps the source free choices explicit.
For a repository-ready local repair coordinate it makes one *autonomous*
choice: X_* is the geometric centre of I_2.  Then

    e_* = c_patch X_*^(-1/2-lambda)

is no longer guessed independently.  These are precisely the physical scales
needed by ``KokunoHeatDiscrepancyRepair`` and by the V2 heat-discrepancy target.

This is a local leading-profile component, not the full outer construction.
The source terminal-tail scales (X_tail, c_inf, rho_o), the actual three-bump
compensation, the global velocity/pressure assembly and the held-out PDE gate
remain separate dependencies.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-outer-reserved-patch-schedule-v1"

# The existing repository heat-repair primitive uses autonomous bumps supported
# in 0.8 <= x_* <= 1.7.  The source only requires ordered compact bumps inside
# I_2; keeping this entire implementation support inside I_2 is a fail-closed
# compatibility guard, not a claimed source constant.
_REPAIR_XSTAR_MIN = 0.80
_REPAIR_XSTAR_MAX = 1.70
_DEFAULT_I2_LOG_OFFSET = -17.5

_SOURCE_FORMULAS = {
    "outer_order": (
        "T_d=exp(M_d)+10, P_*>exp(T_d), then lambda,h and a sufficiently large X_R"
    ),
    "continuation_scale": "X_R=110*(C*P_*)^10",
    "long_power_stage": "T_w=60*log(1/lambda), X_w=X_R*exp(T_d+2)",
    "reserved_intervals": (
        "I1=X_w*exp((T_w-25,T_w-20)); I2=X_w*exp((T_w-20,T_w-15)); "
        "I3=X_w*exp((T_w-14,T_w-9)); I4=X_w*exp((T_w-8,T_w-3))"
    ),
    "power_background": (
        "U=0; E=c_patch*(1+eta^2)^(-1)*X^(-1/2-lambda) on the reserved intervals"
    ),
    "amplitudes": (
        "P1=P_* exp(-1/5); e_w=P1*exp(-T_d/2)*exp(-1/2-lambda/2); "
        "c_patch=e_w*X_w^(1/2+lambda)"
    ),
    "repair_scaling": (
        "x_*=X/X_*; E=e_* f x_*^(-1/2-lambda); e_*=c_patch*X_*^(-1/2-lambda)"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_reserved_interval_formulas_executable": True,
    "source_reserved_patch_amplitude_formula_executable": True,
    "source_power_background_profile_executable": True,
    "source_X_star_e_star_relation_executable": True,
    "autonomous_outer_parameter_values": True,
    "autonomous_I2_repair_center": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_terminal_tail_schedule_bound": False,
    "actual_heat_discrepancy_applied_to_patch": False,
    "heat_compensation_completed": False,
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


def _exp_if_float(log_value: float) -> float | None:
    max_log = math.log(np.finfo(float).max)
    min_log = math.log(np.nextafter(0.0, 1.0))
    if log_value > max_log or log_value < min_log:
        return None
    value = math.exp(log_value)
    return value if math.isfinite(value) and value > 0.0 else None


@dataclass(frozen=True)
class KokunoOuterReservedPatchSchedule:
    """Source-aligned reserved-patch scales with explicit autonomous choices.

    ``M_d``, ``log_p_star_margin``, ``C`` and ``lambda_outer`` instantiate
    source free choices.  They are not recovered OpenAI/Kokuno hidden numbers.
    ``repair_log_offset`` selects an interior point of I_2; the default -17.5
    is the logarithmic midpoint of the source interval (-20,-15).
    """

    M_d: float = 3.0
    log_p_star_margin: float = 1.0
    C: float = 2.0
    lambda_outer: float = 0.05
    h: float = 0.005
    repair_log_offset: float = _DEFAULT_I2_LOG_OFFSET

    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        M_d = float(self.M_d)
        margin = float(self.log_p_star_margin)
        C = float(self.C)
        lam = float(self.lambda_outer)
        h = float(self.h)
        offset = float(self.repair_log_offset)
        if not math.isfinite(M_d) or not (0.0 < M_d <= 6.0):
            raise ValueError("M_d must be finite and lie in (0,6]")
        if not math.isfinite(margin) or not (0.0 < margin <= 20.0):
            raise ValueError("log_p_star_margin must lie in (0,20]")
        if not math.isfinite(C) or not (1.0 < C <= 1.0e3):
            raise ValueError("C must be finite and lie in (1,1e3]")
        if not math.isfinite(lam) or not (0.0 < lam < math.exp(-25.0 / 60.0)):
            raise ValueError("lambda_outer must satisfy 0<lambda<exp(-25/60)")
        if not math.isfinite(h) or not (0.0 < h < 1.0e-2):
            raise ValueError("h must satisfy the source construction bound 0<h<1e-2")
        if not math.isfinite(offset) or not (-20.0 < offset < -15.0):
            raise ValueError("repair_log_offset must lie strictly inside I2 offsets (-20,-15)")

        # Keep the complete autonomous bump support used by #271 inside I_2.
        if offset + math.log(_REPAIR_XSTAR_MIN) <= -20.0:
            raise ValueError("repair support would cross the lower I2 edge")
        if offset + math.log(_REPAIR_XSTAR_MAX) >= -15.0:
            raise ValueError("repair support would cross the upper I2 edge")

        object.__setattr__(self, "M_d", M_d)
        object.__setattr__(self, "log_p_star_margin", margin)
        object.__setattr__(self, "C", C)
        object.__setattr__(self, "lambda_outer", lam)
        object.__setattr__(self, "h", h)
        object.__setattr__(self, "repair_log_offset", offset)
        object.__setattr__(self, "_coordinates", KokunoNativeSimilarityCoordinates(h=h))

    @property
    def T_d(self) -> float:
        return math.exp(self.M_d) + 10.0

    @property
    def log_P_star(self) -> float:
        # Autonomous finite instantiation of the source strict inequality
        # log(P_*) > T_d.
        return self.T_d + self.log_p_star_margin

    @property
    def T_w(self) -> float:
        return 60.0 * math.log(1.0 / self.lambda_outer)

    @property
    def log_X_R(self) -> float:
        return math.log(110.0) + 10.0 * (math.log(self.C) + self.log_P_star)

    @property
    def log_X_w(self) -> float:
        return self.log_X_R + self.T_d + 2.0

    @property
    def log_P1(self) -> float:
        # sigma(1-s)=1-sigma(s), hence int_0^1 sigma = 1/2 exactly.
        return self.log_P_star - 0.2

    @property
    def log_e_w(self) -> float:
        return self.log_P1 - 0.5 * self.T_d - 0.5 - 0.5 * self.lambda_outer

    @property
    def log_c_patch(self) -> float:
        return self.log_e_w + (0.5 + self.lambda_outer) * self.log_X_w

    def reserved_log_intervals(self) -> dict[str, tuple[float, float]]:
        base = self.log_X_w + self.T_w
        return {
            "I1": (base - 25.0, base - 20.0),
            "I2": (base - 20.0, base - 15.0),
            "I3": (base - 14.0, base - 9.0),
            "I4": (base - 8.0, base - 3.0),
        }

    @property
    def log_X_star(self) -> float:
        return self.log_X_w + self.T_w + self.repair_log_offset

    @property
    def log_e_star(self) -> float:
        return self.log_c_patch - (0.5 + self.lambda_outer) * self.log_X_star

    def repair_scales(self) -> dict[str, float]:
        """Materialize the source-derived ``X_*`` and ``e_*`` when representable."""

        X_star = _exp_if_float(self.log_X_star)
        e_star = _exp_if_float(self.log_e_star)
        if X_star is None or e_star is None:
            raise OverflowError(
                "reserved-patch scales are outside float range; use log_scale_report()"
            )
        return {"X_star": X_star, "e_star": e_star}

    def log_scale_report(self) -> dict[str, Any]:
        intervals = self.reserved_log_intervals()
        return {
            "T_d": self.T_d,
            "log_P_star": self.log_P_star,
            "T_w": self.T_w,
            "log_X_R": self.log_X_R,
            "log_X_w": self.log_X_w,
            "log_P1": self.log_P1,
            "log_e_w": self.log_e_w,
            "log_c_patch": self.log_c_patch,
            "reserved_log_intervals": {key: list(value) for key, value in intervals.items()},
            "log_X_star": self.log_X_star,
            "log_e_star": self.log_e_star,
            "X_star": _exp_if_float(self.log_X_star),
            "e_star": _exp_if_float(self.log_e_star),
        }

    @staticmethod
    def source_f(eta: Any) -> np.ndarray:
        values = _finite_array(eta, "eta")
        if np.any(np.abs(values) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        return 1.0 / (1.0 + values * values)

    @staticmethod
    def source_f_eta(eta: Any) -> np.ndarray:
        values = _finite_array(eta, "eta")
        if np.any(np.abs(values) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        return -2.0 * values / np.power(1.0 + values * values, 2)

    def patch_profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate the exact reserved-patch base ``E,F,U,v0`` and first partials."""

        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any(X_array <= 0.0):
            raise ValueError("reserved-patch profile requires X>0")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("reserved-patch profile requires |eta|<=1")

        log_X = np.log(X_array)
        low, high = self.reserved_log_intervals()["I2"]
        if np.any(log_X <= low) or np.any(log_X >= high):
            raise ValueError("reserved-patch profile evaluation is restricted to the open I2 interval")

        x_star = np.exp(log_X - self.log_X_star)
        f = self.source_f(eta_array)
        f_eta = self.source_f_eta(eta_array)
        exponent = -0.5 - self.lambda_outer
        e_star = math.exp(self.log_e_star)
        E = e_star * f * np.power(x_star, exponent)
        E_X = exponent * E / X_array
        E_eta = e_star * f_eta * np.power(x_star, exponent)
        sqrt_2X = np.sqrt(2.0 * X_array)
        F = E / sqrt_2X
        F_X = -(1.0 + self.lambda_outer) * F / X_array
        F_eta = E_eta / sqrt_2X
        zeros = np.zeros_like(F)
        return {
            "x_star": x_star,
            "E": E,
            "E_X": E_X,
            "E_eta": E_eta,
            "F": F,
            "F_X": F_X,
            "F_eta": F_eta,
            "U": zeros,
            "U_X": zeros,
            "U_eta": zeros,
            "v0": zeros,
            "Pi_X": F * F,
        }

    def patch_velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        """Evaluate the source pure-swirl leading velocity on the reserved I2 patch.

        This is deliberately domain-restricted.  It is not a global candidate and
        refuses points whose similarity radius leaves I2.
        """

        coordinates = self._coordinates.evaluate(x, y, z, t)
        profiles = self.patch_profile_values(coordinates["X"], coordinates["eta"])
        factor = np.power(coordinates["q"], -1.0 - self.h) * profiles["F"]
        x_array, y_array = np.broadcast_arrays(_finite_array(x, "x"), _finite_array(y, "y"))
        x_array = np.broadcast_to(x_array, np.shape(factor))
        y_array = np.broadcast_to(y_array, np.shape(factor))
        zeros = np.zeros_like(factor)
        return np.stack((-y_array * factor, x_array * factor, zeros), axis=-1)

    def _unsigned_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "provenance": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": dict(_SOURCE_FORMULAS),
            "parameters": {
                "M_d": self.M_d,
                "log_p_star_margin": self.log_p_star_margin,
                "C": self.C,
                "lambda_outer": self.lambda_outer,
                "h": self.h,
                "repair_log_offset": self.repair_log_offset,
            },
            "autonomous_numerics": {
                "P_star_instantiation": "log(P_*)=T_d+log_p_star_margin",
                "repair_scale_choice": "X_* is the geometric centre of I2 by default",
                "repair_x_star_support": [_REPAIR_XSTAR_MIN, _REPAIR_XSTAR_MAX],
            },
            "derived_log_scales": self.log_scale_report(),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self._unsigned_payload()).encode("utf-8")).hexdigest()

    def to_payload(self) -> dict[str, Any]:
        payload = self._unsigned_payload()
        payload["sha256"] = self.sha256
        return payload

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoOuterReservedPatchSchedule":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno outer reserved-patch schedule schema")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("parameters are missing")
        obj = cls(**params)
        if obj.to_payload() != payload:
            raise ValueError("outer reserved-patch schedule payload hash or content mismatch")
        return obj

    def save(self, path: str | Path) -> None:
        Path(path).write_text(_canonical_json(self.to_payload()) + "\n", encoding="utf-8")

    @classmethod
    def load(cls, path: str | Path) -> "KokunoOuterReservedPatchSchedule":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
