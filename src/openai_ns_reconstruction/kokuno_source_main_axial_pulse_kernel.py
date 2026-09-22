"""Executable low-dimensional kernel for Kokuno's public main axial pulse.

Pinned provenance: KokunoYumeto corrected 2026-09-09 reconstruction,
commit 143f6773feb424ad9ed3a8d116653200f20346b7.

At the pulse-entry seam X_p the public reconstruction uses

    y = log(X/X_p),  xi = lambda*y,
    E = E_p exp[-(1/2+lambda)y],
    U = E R_b,
    R_b = Amp(eta) R_0(xi) + c1(eta) beta1(y) + c2(eta) beta2(y),
    R_0(xi) = phi(xi) [1-sigma(xi-10)],
    phi(xi) = int_0^xi sigma(v/0.02) dv.

This bounded Agent-1 increment materializes only the public main kernel R_0
and an explicitly repository-autonomous principal-amplitude proxy

    A_principal = sqrt((1-exp(-26))/(4 K_b)),
    K_b = int_0^13 exp(-2 xi) R_0(xi)^2 dxi,

obtained by retaining the public principal spin-moment term.  The source exact
Amp(eta) root and c1/c2 end compensators are NOT materialized.  The executable
range therefore stops at xi=11, where R_0 vanishes, before those compensators.

This module is a reusable F/E/U profile contract for composition with the
current #1088 Cartesian candidate.  It is not itself a global velocity field,
not paper/OpenAI exact, and not independent PDE validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any, Mapping

import numpy as np


SCHEMA = "kokuno-source-main-axial-pulse-kernel-v1"
PARENT_EXACT_HEAD = "bd85ed48323feb3b33b9c404a078630245dfe783"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected 208-page reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
QUADRATURE_ORDER = 96
MAIN_XI_MAX = 11.0

_SOURCE_FORMULAS = {
    "pulse_coordinate": "y=log(X/X_p), xi=lambda*y",
    "pulse_E": "E=E_p exp[-(1/2+lambda)y]",
    "pulse_U": "U=E R_b",
    "pulse_ratio": "R_b=Amp(eta) R_0(lambda y)+c1(eta) beta1(y)+c2(eta) beta2(y)",
    "main_kernel": "R_0(xi)=phi(xi)[1-sigma(xi-10)]",
    "phi": "phi(xi)=int_0^xi sigma(v/0.02) dv",
    "sigma": "0 below 0, 1 above 1, standard C-infinity quotient on (0,1)",
    "principal_spin_term": "4 K_b Amp^2 approximately equals 1-exp(-26)",
    "K_b": "K_b=int_0^13 exp(-2 xi) R_0(xi)^2 dxi",
    "F": "F=E/sqrt(2X)",
}

_NUMERICAL_REALIZATION = {
    "principal_amplitude": (
        "repository-autonomous A_principal=sqrt((1-exp(-26))/(4 K_b)); "
        "not the source exact Amp(eta) root"
    ),
    "quadrature": (
        "fixed 96-point Gauss-Legendre panels; phi uses quadrature only for "
        "0<xi<0.02; K_b uses panels [0,0.02],[0.02,10],[10,11]"
    ),
    "domain_stop": "fail closed for xi<0 or xi>11; c1/c2 bumps are absent",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_main_pulse_kernel_formula_materialized": True,
    "source_smooth_cutoff_materialized": True,
    "source_phi_materialized": True,
    "repository_autonomous_principal_amplitude_materialized": True,
    "source_exact_amplitude_root_materialized": False,
    "source_pulse_end_MJ_corrections_materialized": False,
    "source_hidden_parameters_recovered": False,
    "current_cartesian_pulse_velocity_composed": False,
    "matched_global_pressure_materialized": False,
    "restricted_forcing_materialized": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def smooth_sigma(z: Any) -> np.ndarray:
    """Public C-infinity step sigma, evaluated stably and vectorized."""
    x = _finite(z, "z")
    out = np.empty_like(x, dtype=float)
    lo, hi = x <= 0.0, x >= 1.0
    mid = ~(lo | hi)
    out[lo], out[hi] = 0.0, 1.0
    if np.any(mid):
        q = x[mid]
        log_ratio = -1.0 / (q * q) + 1.0 / ((1.0 - q) ** 2)
        vals = np.empty_like(log_ratio)
        nonneg = log_ratio >= 0.0
        vals[nonneg] = 1.0 / (1.0 + np.exp(-log_ratio[nonneg]))
        e = np.exp(log_ratio[~nonneg])
        vals[~nonneg] = e / (1.0 + e)
        out[mid] = vals
    return out


def smooth_sigma_prime(z: Any) -> np.ndarray:
    """Analytic derivative of smooth_sigma."""
    x = _finite(z, "z")
    out = np.zeros_like(x, dtype=float)
    mid = (x > 0.0) & (x < 1.0)
    if np.any(mid):
        q = x[mid]
        s = smooth_sigma(q)
        out[mid] = s * (1.0 - s) * (2.0 / q**3 + 2.0 / (1.0 - q) ** 3)
    return out


def _gauss_legendre_integral(func, a: float, b: float) -> float:
    if not (math.isfinite(a) and math.isfinite(b) and b >= a):
        raise ValueError("quadrature bounds must be finite and ordered")
    if b == a:
        return 0.0
    nodes, weights = np.polynomial.legendre.leggauss(QUADRATURE_ORDER)
    half, mid = 0.5 * (b - a), 0.5 * (b + a)
    x = half * nodes + mid
    return float(half * np.sum(weights * np.asarray(func(x), dtype=float)))


def phi(xi: Any) -> np.ndarray:
    """Public phi=int_0^xi sigma(v/0.02)dv, materialized for xi>=0."""
    x = _finite(xi, "xi")
    if np.any(x < 0.0):
        raise ValueError("phi is materialized only for xi>=0")
    out = np.zeros_like(x, dtype=float)
    startup = (x > 0.0) & (x < 0.02)
    if np.any(startup):
        flat, out_flat = x.reshape(-1), out.reshape(-1)
        for i in np.flatnonzero(startup.reshape(-1)):
            upper = float(flat[i])
            out_flat[i] = _gauss_legendre_integral(
                lambda v: smooth_sigma(v / 0.02), 0.0, upper
            )
    late = x >= 0.02
    # sigma(s)+sigma(1-s)=1, hence int_0^0.02 sigma(v/0.02)dv=0.01.
    out[late] = x[late] - 0.01
    return out


def main_kernel_R0(xi: Any) -> np.ndarray:
    """Public source main-pulse kernel R_0(xi)."""
    x = _finite(xi, "xi")
    if np.any((x < 0.0) | (x > MAIN_XI_MAX)):
        raise ValueError(f"main source pulse kernel requires 0<=xi<={MAIN_XI_MAX}")
    return phi(x) * (1.0 - smooth_sigma(x - 10.0))


def main_kernel_R0_prime(xi: Any) -> np.ndarray:
    """Analytic xi derivative of R_0."""
    x = _finite(xi, "xi")
    if np.any((x < 0.0) | (x > MAIN_XI_MAX)):
        raise ValueError(f"main source pulse kernel requires 0<=xi<={MAIN_XI_MAX}")
    return (
        smooth_sigma(x / 0.02) * (1.0 - smooth_sigma(x - 10.0))
        - phi(x) * smooth_sigma_prime(x - 10.0)
    )


def principal_Kb() -> float:
    """Deterministic repository realization of the public principal K_b integral."""
    def integrand(x: np.ndarray) -> np.ndarray:
        r = main_kernel_R0(x)
        return np.exp(-2.0 * x) * r * r

    return sum(
        _gauss_legendre_integral(integrand, a, b)
        for a, b in ((0.0, 0.02), (0.02, 10.0), (10.0, 11.0))
    )


def principal_amplitude() -> float:
    """Autonomous principal-term amplitude; not source-exact Amp(eta)."""
    kb = principal_Kb()
    if not (kb > 0.0 and math.isfinite(kb)):
        raise RuntimeError("principal K_b must be positive and finite")
    return math.sqrt((1.0 - math.exp(-26.0)) / (4.0 * kb))


class KokunoSourceMainAxialPulseKernel:
    """Vectorized F/E/U contract for the main pulse with autonomous A0."""

    @property
    def K_b(self) -> float:
        return principal_Kb()

    @property
    def A_principal(self) -> float:
        return principal_amplitude()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return dict(_TRUTH_BOUNDARY)

    def profile_values(
        self, X: Any, eta: Any, *, X_p: float, E_entry: Any, lambda_value: float
    ) -> dict[str, np.ndarray]:
        if not (math.isfinite(X_p) and X_p > 0.0):
            raise ValueError("X_p must be finite and positive")
        if not (math.isfinite(lambda_value) and lambda_value > 0.0):
            raise ValueError("lambda_value must be finite and positive")
        X_arr, eta_arr, E0_arr = np.broadcast_arrays(
            _finite(X, "X"), _finite(eta, "eta"), _finite(E_entry, "E_entry")
        )
        if np.any(X_arr < X_p):
            raise ValueError("main pulse requires X>=X_p")
        if np.any(E0_arr <= 0.0):
            raise ValueError("E_entry must be strictly positive")
        y = np.log(X_arr / X_p)
        xi = lambda_value * y
        if np.any(xi > MAIN_XI_MAX + 2.0e-13):
            raise ValueError("main pulse stops at xi=11 before source end compensators")
        xi = np.minimum(xi, MAIN_XI_MAX)
        decay = np.exp(-(0.5 + lambda_value) * y)
        E = E0_arr * decay
        F = E / np.sqrt(2.0 * X_arr)
        R0 = main_kernel_R0(xi)
        U = E * self.A_principal * R0
        return {
            "X": X_arr,
            "eta": eta_arr,
            "y": y,
            "xi": xi,
            "R0": R0,
            "amplitude_principal": np.full_like(X_arr, self.A_principal),
            "E_source_main_pulse_principal": E,
            "F_source_main_pulse_principal": F,
            "U_source_main_pulse_principal": U,
        }

    def radial_derivatives(
        self, X: Any, eta: Any, *, X_p: float, E_entry: Any, lambda_value: float
    ) -> dict[str, np.ndarray]:
        values = self.profile_values(
            X, eta, X_p=X_p, E_entry=E_entry, lambda_value=lambda_value
        )
        X_arr = np.asarray(values["X"], dtype=float)
        E = np.asarray(values["E_source_main_pulse_principal"], dtype=float)
        F = np.asarray(values["F_source_main_pulse_principal"], dtype=float)
        U = np.asarray(values["U_source_main_pulse_principal"], dtype=float)
        xi = np.asarray(values["xi"], dtype=float)
        a = 0.5 + lambda_value
        return {
            **values,
            "E_X_source_main_pulse_principal": -a * E / X_arr,
            "F_X_source_main_pulse_principal": -(1.0 + lambda_value) * F / X_arr,
            "U_X_source_main_pulse_principal": (
                -a * U
                + E * self.A_principal * lambda_value * main_kernel_R0_prime(xi)
            ) / X_arr,
        }

    def eta_derivatives(
        self,
        X: Any,
        eta: Any,
        *,
        X_p: float,
        E_entry: Any,
        E_entry_eta: Any,
        lambda_value: float,
    ) -> dict[str, np.ndarray]:
        values = self.profile_values(
            X, eta, X_p=X_p, E_entry=E_entry, lambda_value=lambda_value
        )
        X_arr, _, E0_eta = np.broadcast_arrays(
            np.asarray(values["X"], dtype=float),
            _finite(E_entry, "E_entry"),
            _finite(E_entry_eta, "E_entry_eta"),
        )
        decay = np.exp(-(0.5 + lambda_value) * np.asarray(values["y"], dtype=float))
        E_eta = E0_eta * decay
        F_eta = E_eta / np.sqrt(2.0 * X_arr)
        U_eta = E_eta * self.A_principal * np.asarray(values["R0"], dtype=float)
        return {
            **values,
            "E_eta_source_main_pulse_principal": E_eta,
            "F_eta_source_main_pulse_principal": F_eta,
            "U_eta_source_main_pulse_principal": U_eta,
        }

    def to_configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": SOURCE_RELEASE,
                "release_date": SOURCE_RELEASE_DATE,
            },
            "numerical_realization": {
                "quadrature_order": QUADRATURE_ORDER,
                "main_xi_max": MAIN_XI_MAX,
                "amplitude_model": "repository_autonomous_principal_spin_term",
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "KokunoSourceMainAxialPulseKernel":
        if dict(payload) != cls().to_configuration():
            raise ValueError("main-pulse configuration/source/truth boundary is immutable")
        return cls()

    @property
    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_configuration()).encode()).hexdigest()

    def receipt(self) -> dict[str, Any]:
        kb, amp = self.K_b, self.A_principal
        return {
            "schema": SCHEMA,
            "semantic_sha256": self.semantic_sha256,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "source": self.to_configuration()["source"],
            "source_formulas": dict(_SOURCE_FORMULAS),
            "numerical_realization": dict(_NUMERICAL_REALIZATION),
            "principal_amplitude_diagnostic": {
                "K_b": kb,
                "A_principal": amp,
                "principal_identity_error": abs(
                    4.0 * kb * amp * amp - (1.0 - math.exp(-26.0))
                ),
            },
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }
