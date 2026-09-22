"""Source-ideal Kokuno exterior Q_s transport and matching-length target.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
blob 205a99807302e21a51c5eaf223390c0dfc42bcd0
corrected 2026-09-09 reconstruction.

The corrected reconstruction states that, after the relative-swirl correction,
its *ideal* leading state has

    I = X H/(1-lambda),   U=M=J=0,

and therefore at the start of the exterior transition

    Q_s = -1 + (1-h) I/(XH) = (lambda-h)/(1-lambda).

It then gives

    Q_s' + (1+l) Q_s = -l-h,

through the first unit transition, the l=-1 hold of length 4 log(1/h),
and the second unit transition.  Before the terminal multiplier interval the
source holds l=-h until Q_s reaches the terminal target Q_p.

This module executes that SOURCE-IDEAL scalar transport and the corresponding
ideal matching length.  It intentionally does NOT relabel the result as the
current repository candidate's Q_s: the current A1 lineage preserves its
actual post-pulse M/M_eta history instead of imposing the source-ideal
U=M=J=0 identity exactly.  Hence this module does not extend Cartesian
velocity and does not close the current-Q_s blocker.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_postswirl_release1 import (
    _sigma,
    _sigma_integral,
)
from .kokuno_public_terminal_multiplier_target import (
    KokunoPublicTerminalMultiplierTarget,
)

SCHEMA = "kokuno-public-ideal-exterior-qs-schedule-v1"
PARENT_EXACT_HEAD = "fb81b8e90272caa9e8bbd075f05b94dbb1195310"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
QS_QUADRATURE_ORDER = 192
RELEASE_LENGTH = 1.0

_SOURCE_FORMULAS = {
    "ideal_exterior_entry": (
        "after the correction I=XH/(1-lambda), U=M=J=0, hence "
        "Q_s=(lambda-h)/(1-lambda)"
    ),
    "transport": "Q_s'+(1+l)Q_s=-l-h",
    "release1": "l=-lambda-(1-lambda)*sigma(s), 0<=s<=1",
    "hold": "l=-1 for T_hold=4*log(1/h), so Q_s'=1-h",
    "release2": "l=-1+(1-h)*sigma(s), 0<=s<=1",
    "terminal_target_order": "hold l=-h before terminal interval until Q_s reaches Q_p",
    "matching_bridge": (
        "Q_s(y)=Q_in*exp(-(1-h)y), "
        "T_match=(1-h)^(-1)*log(Q_in/Q_p)"
    ),
}

_NUMERICAL_REALIZATION = {
    "quadrature": "fixed 192-point Gauss-Legendre for the two unit-transition Duhamel integrals",
    "sigma_integral": "reuse A1 fixed public-step quadrature implementation",
    "arithmetic": "binary64 NumPy/Python arithmetic",
    "terminal_target": "consume exact #1213 Q_p without changing c_o or its quadrature",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "public_source_ideal_exterior_qs_schedule_materialized": True,
    "public_source_ideal_qs_release2_endpoint_materialized": True,
    "public_source_ideal_l_minus_h_matching_length_materialized": True,
    "current_q_s_release2_endpoint_materialized": False,
    "current_l_minus_h_matching_bridge_materialized": False,
    "current_cartesian_terminal_multiplier_composed": False,
    "source_terminal_multiplier_materialized": False,
    "source_exterior_heat_replacement_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "unified_global_cartesian_velocity_export_ready": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "same_protocol_comparable_to_st006": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}

_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(QS_QUADRATURE_ORDER)


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _integral_0_upper(func: Callable[[np.ndarray], np.ndarray], upper: Any) -> np.ndarray:
    """Fixed Gauss-Legendre integral of a vectorized function from 0 to upper."""
    upper_arr = _finite(upper, "upper integration limit")
    tol = 64.0 * np.finfo(float).eps
    if np.any((upper_arr < -tol) | (upper_arr > 1.0 + tol)):
        raise ValueError("unit-transition integration limit must lie in [0,1]")
    clipped = np.clip(upper_arr, 0.0, 1.0)
    flat = clipped.reshape(-1)
    out = np.empty_like(flat)
    for i, u in enumerate(flat):
        if u == 0.0:
            out[i] = 0.0
            continue
        x = 0.5 * u * (_GL_NODES + 1.0)
        vals = np.asarray(func(x), dtype=float)
        if vals.shape != x.shape or np.any(~np.isfinite(vals)):
            raise RuntimeError("Q_s transport quadrature integrand is invalid")
        out[i] = 0.5 * u * float(np.dot(_GL_WEIGHTS, vals))
    return out.reshape(clipped.shape)


@dataclass(frozen=True)
class KokunoPublicIdealExteriorQsSchedule:
    """Executable source-ideal Q_s schedule; not current-candidate Q_s evidence."""

    parent: KokunoPublicTerminalMultiplierTarget = field(
        default_factory=KokunoPublicTerminalMultiplierTarget,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPublicTerminalMultiplierTarget):
            raise TypeError("parent must be KokunoPublicTerminalMultiplierTarget")
        truth = self.parent.truth_boundary
        if not truth["public_terminal_multiplier_target_api_materialized"]:
            raise ValueError("exact A1 #1213 terminal Q_p target is required")
        if truth["current_q_s_release2_endpoint_materialized"]:
            raise ValueError("parent unexpectedly claims current Q_s materialization")
        if not (0.0 < self.h_value < 0.01):
            raise ValueError("h must satisfy the public construction range 0<h<1/100")
        if not (0.0 < self.lambda_value < 1.0):
            raise ValueError("lambda must lie in (0,1)")
        if not self.q_in_source_ideal > self.q_p:
            raise ValueError("source-ideal matching bridge requires Q_in>Q_p")

    @property
    def release2(self):
        return self.parent.parent

    @property
    def hold(self):
        return self.release2.parent

    @property
    def release1(self):
        return self.hold.parent

    @property
    def h_value(self) -> float:
        return float(self.release2.h_value)

    @property
    def lambda_value(self) -> float:
        return float(self.release2.lambda_value)

    @property
    def q_p(self) -> float:
        return float(self.parent.q_p)

    @property
    def source_formulas(self) -> dict[str, str]:
        return dict(_SOURCE_FORMULAS)

    @property
    def numerical_realization(self) -> dict[str, str]:
        return dict(_NUMERICAL_REALIZATION)

    @property
    def truth_boundary(self) -> dict[str, bool]:
        truth = dict(self.parent.truth_boundary)
        truth.update(_TRUTH_UPDATES)
        return truth

    @property
    def q_exterior_entry_source_ideal(self) -> float:
        lam = self.lambda_value
        return (lam - self.h_value) / (1.0 - lam)

    @property
    def hold_length(self) -> float:
        return 4.0 * math.log(1.0 / self.h_value)

    def _release1_A(self, s: Any) -> np.ndarray:
        ss = _finite(s, "release1 s")
        return (1.0 - self.lambda_value) * (ss - _sigma_integral(ss))

    def _release2_A(self, s: Any) -> np.ndarray:
        ss = _finite(s, "release2 s")
        return (1.0 - self.h_value) * _sigma_integral(ss)

    def q_release1(self, s: Any) -> np.ndarray:
        """Source-ideal Q_s on the first one-unit exterior transition."""
        ss = _finite(s, "release1 s")
        A = self._release1_A(ss)
        lam, h = self.lambda_value, self.h_value

        def integrand(v: np.ndarray) -> np.ndarray:
            Av = self._release1_A(v)
            ell = -lam - (1.0 - lam) * _sigma(v)
            b = -ell - h
            return np.exp(Av) * b

        source = _integral_0_upper(integrand, ss)
        out = np.exp(-A) * (self.q_exterior_entry_source_ideal + source)
        if np.any(~np.isfinite(out)) or np.any(out <= 0.0):
            raise RuntimeError("source-ideal release1 Q_s became invalid")
        return out

    @property
    def q_release1_end_source_ideal(self) -> float:
        return float(np.asarray(self.q_release1(1.0)))

    def q_hold(self, s: Any) -> np.ndarray:
        """Source-ideal Q_s on the constant l=-1 hold."""
        ss = _finite(s, "hold s")
        tol = 128.0 * np.finfo(float).eps * max(1.0, self.hold_length)
        if np.any((ss < -tol) | (ss > self.hold_length + tol)):
            raise ValueError("hold coordinate must lie in [0,4 log(1/h)]")
        sc = np.clip(ss, 0.0, self.hold_length)
        out = self.q_release1_end_source_ideal + (1.0 - self.h_value) * sc
        if np.any(~np.isfinite(out)) or np.any(out <= 0.0):
            raise RuntimeError("source-ideal hold Q_s became invalid")
        return out

    @property
    def q_hold_end_source_ideal(self) -> float:
        return float(np.asarray(self.q_hold(self.hold_length)))

    def q_release2(self, s: Any) -> np.ndarray:
        """Source-ideal Q_s on the second one-unit -1 -> -h transition."""
        ss = _finite(s, "release2 s")
        A = self._release2_A(ss)
        h = self.h_value

        def integrand(v: np.ndarray) -> np.ndarray:
            Av = self._release2_A(v)
            sig = _sigma(v)
            b = (1.0 - h) * (1.0 - sig)
            return np.exp(Av) * b

        source = _integral_0_upper(integrand, ss)
        out = np.exp(-A) * (self.q_hold_end_source_ideal + source)
        if np.any(~np.isfinite(out)) or np.any(out <= 0.0):
            raise RuntimeError("source-ideal release2 Q_s became invalid")
        return out

    @property
    def q_in_source_ideal(self) -> float:
        return float(np.asarray(self.q_release2(1.0)))

    @property
    def ideal_matching_length(self) -> float:
        value = math.log(self.q_in_source_ideal / self.q_p) / (1.0 - self.h_value)
        if not (math.isfinite(value) and value > 0.0):
            raise RuntimeError("source-ideal matching length must be finite and positive")
        return value

    def q_matching_bridge_source_ideal(self, y: Any) -> np.ndarray:
        """Ideal l=-h exponential bridge target, not current-candidate Q_s."""
        yy = _finite(y, "matching y")
        T = self.ideal_matching_length
        tol = 128.0 * np.finfo(float).eps * max(1.0, T)
        if np.any((yy < -tol) | (yy > T + tol)):
            raise ValueError("matching coordinate must lie in [0,T_match]")
        yc = np.clip(yy, 0.0, T)
        out = self.q_in_source_ideal * np.exp(-(1.0 - self.h_value) * yc)
        if np.any(~np.isfinite(out)) or np.any(out <= 0.0):
            raise RuntimeError("source-ideal matching Q_s became invalid")
        return out

    def release1_ode_rhs(self, s: Any) -> np.ndarray:
        ss = _finite(s, "release1 s")
        q = self.q_release1(ss)
        ell = -self.lambda_value - (1.0 - self.lambda_value) * _sigma(ss)
        return -ell - self.h_value - (1.0 + ell) * q

    def release2_ode_rhs(self, s: Any) -> np.ndarray:
        ss = _finite(s, "release2 s")
        q = self.q_release2(ss)
        ell = -1.0 + (1.0 - self.h_value) * _sigma(ss)
        return -ell - self.h_value - (1.0 + ell) * q

    def schedule_report(self) -> dict[str, Any]:
        q0 = self.q_exterior_entry_source_ideal
        q1 = self.q_release1_end_source_ideal
        qh = self.q_hold_end_source_ideal
        qin = self.q_in_source_ideal
        T = self.ideal_matching_length
        qend = float(np.asarray(self.q_matching_bridge_source_ideal(T)))
        return {
            "schema": "kokuno-agent1-public-ideal-exterior-qs-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "lambda_current": self.lambda_value,
            "h_current": self.h_value,
            "source_ideal_assumptions": ["I=XH/(1-lambda)", "U=0", "M=0", "J=0"],
            "q_exterior_entry_source_ideal": q0,
            "q_release1_end_source_ideal": q1,
            "hold_length": self.hold_length,
            "q_hold_end_source_ideal": qh,
            "q_in_source_ideal": qin,
            "q_p": self.q_p,
            "ideal_matching_length": T,
            "q_matching_end_source_ideal": qend,
            "matching_end_relative_error_to_q_p": abs(qend - self.q_p) / self.q_p,
            "current_candidate_q_s_claimed": False,
            "current_cartesian_bridge_composed": False,
            "truth_boundary": self.truth_boundary,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source_repository": SOURCE_REPOSITORY,
            "source_commit": SOURCE_COMMIT,
            "source_blob": SOURCE_BLOB,
            "source_path": SOURCE_PATH,
            "source_release": SOURCE_RELEASE,
            "source_release_date": SOURCE_RELEASE_DATE,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "lambda_current": self.lambda_value,
            "h_current": self.h_value,
            "q_s_quadrature_order": QS_QUADRATURE_ORDER,
            "source_formulas": self.source_formulas,
            "numerical_realization": self.numerical_realization,
            "truth_boundary": self.truth_boundary,
        }

    @property
    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]):
        obj = cls()
        if dict(payload) != obj.configuration():
            raise ValueError("configuration/provenance drift")
        return obj

    @classmethod
    def load_configuration(cls, path: str | Path):
        return cls.from_configuration(json.loads(Path(path).read_text()))
