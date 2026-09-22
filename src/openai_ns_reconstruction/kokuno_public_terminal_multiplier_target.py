"""Public Kokuno terminal-multiplier target algebra after the -1 -> -h release.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
blob 205a99807302e21a51c5eaf223390c0dfc42bcd0
corrected 2026-09-09 reconstruction.

The corrected source does NOT apply the terminal multiplier immediately after
the one-unit l=-1 -> -h release.  It first defines

    f_o(y) = 1 - rho_o psi_o(y),      0 <= y <= 3,
    psi_o(y) = 1 - sigma((y-1)/2),
    rho_o = c_o h > 0,

with c_o fixed sufficiently small that 0 <= f_o'/f_o < h/4, and the terminal
target

    Q_p = int_0^3 exp(int_0^v (1+l(s)) ds) (f_o'(v)/f_o(v)) dv,
    l = -h + f_o'/f_o.

The current Q_s must first decay on an intervening l=-h bridge until it reaches
Q_p.  This module materializes ONLY the terminal multiplier / Q_p target
algebra.  It does not manufacture current Q_s, does not compose the bridge, and
does not extend the Cartesian candidate past A1 #1204.

The public source does not expose a unique numerical c_o.  This repository
freezes c_o=1/64 as an autonomous, residual-independent realization and verifies
the source inequality against the inherited public flat step.  It is not
claimed source-exact.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_postswirl_release1 import _sigma, _sigma_prime
from .kokuno_pa16_current_cartesian_postswirl_release2 import (
    KokunoPA16CurrentCartesianPostSwirlRelease2,
)

SCHEMA = "kokuno-public-terminal-multiplier-target-v1"
PARENT_EXACT_HEAD = "18e75e6e2df43206147db34788c68ee6a7fa0f37"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"

TERMINAL_LENGTH = 3.0
C_O_AUTONOMOUS = 1.0 / 64.0
QP_QUADRATURE_ORDER = 192
# A deliberately loose repository certificate for the inherited flat step.
# The actual maximum is 8 at z=1/2; the regression scans the implemented
# derivative independently.  This is numerical-realization metadata, not a
# source constant.
SIGMA_PRIME_CERTIFICATE_BOUND = 9.0

_SOURCE_FORMULAS = {
    "terminal_multiplier": "f_o(y)=1-rho_o*psi_o(y), 0<=y<=3",
    "terminal_cutoff": "psi_o(y)=1-sigma((y-1)/2)",
    "terminal_amplitude": "rho_o=c_o*h>0 with c_o fixed so 0<=f_o'/f_o<h/4",
    "terminal_slope": "l=-h+f_o'/f_o",
    "terminal_target": (
        "Q_p=int_0^3 exp(int_0^v(1+l(s))ds)*(f_o'(v)/f_o(v))dv"
    ),
    "bridge_order": (
        "after -1->-h release, hold l=-h until current Q_s reaches Q_p; "
        "only then apply terminal interval"
    ),
    "bridge_solution": (
        "on intervening l=-h bridge: Q_s(y)=Q_in*exp(-(1-h)y), "
        "T_match=(1-h)^(-1)*log(Q_in/Q_p)"
    ),
}

_NUMERICAL_REALIZATION = {
    "c_o": "1/64 fixed repository-autonomous choice; not source-exact",
    "c_o_selection": (
        "chosen before any residual evaluation solely to satisfy the public "
        "0<=f_o'/f_o<h/4 inequality with a wide margin"
    ),
    "q_p_quadrature": "fixed 192-point Gauss-Legendre on [0,3]",
    "profile_arithmetic": "binary64 NumPy/Python arithmetic",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "public_terminal_multiplier_target_api_materialized": True,
    "repository_autonomous_terminal_c_o_materialized": True,
    "source_exact_terminal_c_o_recovered": False,
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

_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(QP_QUADRATURE_ORDER)


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


def _sigma_second(value: Any) -> np.ndarray:
    """Analytic second derivative of the inherited public flat step."""
    z = _finite(value, "sigma argument")
    out = np.zeros_like(z, dtype=float)
    inside = (z > 0.0) & (z < 1.0)
    if np.any(inside):
        x = z[inside]
        sig = _sigma(x)
        lp = 2.0 / (x**3) + 2.0 / ((1.0 - x) ** 3)
        lpp = -6.0 / (x**4) + 6.0 / ((1.0 - x) ** 4)
        out[inside] = sig * (1.0 - sig) * (
            (1.0 - 2.0 * sig) * lp * lp + lpp
        )
    return out


@dataclass(frozen=True)
class KokunoPublicTerminalMultiplierTarget:
    """Executable public terminal-multiplier / Q_p target algebra."""

    parent: KokunoPA16CurrentCartesianPostSwirlRelease2 = field(
        default_factory=KokunoPA16CurrentCartesianPostSwirlRelease2,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianPostSwirlRelease2):
            raise TypeError("parent must be KokunoPA16CurrentCartesianPostSwirlRelease2")
        truth = self.parent.truth_boundary
        if not truth["source_l_minus1_to_minus_h_transition_materialized"]:
            raise ValueError("exact A1 #1204 -1 -> -h source transition is required")
        if not truth["current_cartesian_l_minus1_to_minus_h_composed"]:
            raise ValueError("exact A1 #1204 Cartesian release2 composition is required")
        h = self.h_value
        if not (0.0 < h < 0.01):
            raise ValueError("current h must satisfy the public construction range 0<h<1/100")
        if not (0.0 < self.rho_o < 1.0):
            raise ValueError("autonomous rho_o must lie in (0,1)")
        # Sufficient bound under the repository-fixed derivative certificate:
        # f' <= rho*(1/2)*B and f >= 1-rho.
        certified = (
            0.5
            * self.rho_o
            * SIGMA_PRIME_CERTIFICATE_BOUND
            / (1.0 - self.rho_o)
        )
        if not certified < h / 4.0:
            raise ValueError("autonomous c_o does not satisfy the public slope inequality")

    @property
    def h_value(self) -> float:
        return float(self.parent.h_value)

    @property
    def c_o(self) -> float:
        return C_O_AUTONOMOUS

    @property
    def rho_o(self) -> float:
        return self.c_o * self.h_value

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

    def _broadcast_terminal(self, y: Any) -> np.ndarray:
        arr = _finite(y, "terminal y")
        tol = 128.0 * np.finfo(float).eps
        if np.any((arr < -tol) | (arr > TERMINAL_LENGTH + tol)):
            raise ValueError("terminal multiplier requires 0<=y<=3")
        return np.clip(arr, 0.0, TERMINAL_LENGTH)

    def terminal_multiplier(self, y: Any) -> dict[str, np.ndarray]:
        """Vectorized f_o, slope and analytic y derivatives on 0<=y<=3."""
        yy = self._broadcast_terminal(y)
        z = 0.5 * (yy - 1.0)
        sig = _sigma(z)
        sig_p = _sigma_prime(z)
        sig_pp = _sigma_second(z)

        psi = 1.0 - sig
        rho = self.rho_o
        f = 1.0 - rho * psi
        f_prime = 0.5 * rho * sig_p
        f_second = 0.25 * rho * sig_pp
        dlogf = f_prime / f
        d2logf = f_second / f - dlogf * dlogf
        ell = -self.h_value + dlogf
        ell_prime = d2logf

        arrays = (psi, f, f_prime, f_second, dlogf, d2logf, ell, ell_prime)
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("terminal multiplier became non-finite")
        if np.any(f <= 0.0):
            raise RuntimeError("terminal multiplier lost positivity")
        if np.any(f_prime < -64.0 * np.finfo(float).eps):
            raise RuntimeError("terminal multiplier lost monotonicity")

        return {
            "y": yy,
            "psi_o": psi,
            "f_o": f,
            "f_o_prime": f_prime,
            "f_o_second": f_second,
            "D_y_log_f_o": dlogf,
            "D_y2_log_f_o": d2logf,
            "ell": ell,
            "D_y_ell": ell_prime,
        }

    def terminal_logH_ratio(self, y: Any) -> np.ndarray:
        """log(H/H_entry) for the public terminal multiplier stage."""
        p = self.terminal_multiplier(y)
        yy = np.asarray(p["y"], dtype=float)
        f = np.asarray(p["f_o"], dtype=float)
        f0 = 1.0 - self.rho_o
        return -self.h_value * yy + np.log(f / f0)

    def terminal_logE_ratio(self, y: Any) -> np.ndarray:
        """log(E/E_entry) because H=sqrt(2X)E."""
        yy = self._broadcast_terminal(y)
        return self.terminal_logH_ratio(yy) - 0.5 * yy

    def q_p_direct(self) -> float:
        """Evaluate the public Q_p definition with its exact integrating factor."""
        v = 0.5 * TERMINAL_LENGTH * (_GL_NODES + 1.0)
        p = self.terminal_multiplier(v)
        f = np.asarray(p["f_o"], dtype=float)
        dlogf = np.asarray(p["D_y_log_f_o"], dtype=float)
        f0 = 1.0 - self.rho_o
        log_if = (1.0 - self.h_value) * v + np.log(f / f0)
        integrand = np.exp(log_if) * dlogf
        return 0.5 * TERMINAL_LENGTH * float(np.dot(_GL_WEIGHTS, integrand))

    def q_p_reduced(self) -> float:
        """Equivalent stable identity Q_p=(1/f(0))*int exp((1-h)v) f'(v) dv."""
        v = 0.5 * TERMINAL_LENGTH * (_GL_NODES + 1.0)
        p = self.terminal_multiplier(v)
        fp = np.asarray(p["f_o_prime"], dtype=float)
        f0 = 1.0 - self.rho_o
        integrand = np.exp((1.0 - self.h_value) * v) * fp / f0
        return 0.5 * TERMINAL_LENGTH * float(np.dot(_GL_WEIGHTS, integrand))

    @property
    def q_p(self) -> float:
        value = self.q_p_reduced()
        if not (math.isfinite(value) and value > 0.0):
            raise RuntimeError("terminal target Q_p must be finite and positive")
        return value

    def source_inequality_report(self) -> dict[str, float | bool]:
        # The public flat step has its active derivative on z in (0,1).
        z = np.linspace(0.0, 1.0, 200001, dtype=float)
        sigma_prime_max = float(np.max(_sigma_prime(z)))
        # Since y maps to z=(y-1)/2, max f' is rho/2 max sigma'.
        # The exact maximum ratio is also measured on the terminal grid.
        y = np.linspace(0.0, TERMINAL_LENGTH, 200001, dtype=float)
        dlogf = np.asarray(self.terminal_multiplier(y)["D_y_log_f_o"], dtype=float)
        actual_max = float(np.max(dlogf))
        certified_max = (
            0.5
            * self.rho_o
            * SIGMA_PRIME_CERTIFICATE_BOUND
            / (1.0 - self.rho_o)
        )
        limit = self.h_value / 4.0
        return {
            "sigma_prime_sampled_max": sigma_prime_max,
            "sigma_prime_certificate_bound": SIGMA_PRIME_CERTIFICATE_BOUND,
            "actual_sampled_max_fprime_over_f": actual_max,
            "certified_max_fprime_over_f": certified_max,
            "public_upper_limit_h_over_4": limit,
            "sampled_inequality_satisfied": bool(actual_max < limit),
            "certificate_inequality_satisfied": bool(certified_max < limit),
        }

    def target_report(self) -> dict[str, Any]:
        q_direct = self.q_p_direct()
        q_reduced = self.q_p_reduced()
        inequality = self.source_inequality_report()
        p0 = self.terminal_multiplier(0.0)
        p3 = self.terminal_multiplier(3.0)
        return {
            "schema": "kokuno-agent1-terminal-multiplier-target-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "h_current": self.h_value,
            "c_o_autonomous": self.c_o,
            "rho_o": self.rho_o,
            "terminal_length": TERMINAL_LENGTH,
            "q_p_direct": q_direct,
            "q_p_reduced": q_reduced,
            "q_p_relative_replay_error": abs(q_direct - q_reduced) / q_reduced,
            "q_p_over_h": q_reduced / self.h_value,
            "f_o_start": float(np.asarray(p0["f_o"])),
            "f_o_end": float(np.asarray(p3["f_o"])),
            "ell_start": float(np.asarray(p0["ell"])),
            "ell_end": float(np.asarray(p3["ell"])),
            "source_inequality": inequality,
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
            "h_current": self.h_value,
            "c_o_autonomous": self.c_o,
            "rho_o": self.rho_o,
            "terminal_length": TERMINAL_LENGTH,
            "q_p_quadrature_order": QP_QUADRATURE_ORDER,
            "sigma_prime_certificate_bound": SIGMA_PRIME_CERTIFICATE_BOUND,
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
