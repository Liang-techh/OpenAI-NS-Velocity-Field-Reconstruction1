"""Executable Appendix-B boundary profile at ``X_i=110``.

Pinned provenance is the corrected KokunoYumeto 2026-09-09 reconstruction,
commit 143f6773feb424ad9ed3a8d116653200f20346b7.  The public continuation
prescribes a reference profile, stress activation with

    kappa = 1 - (1-kappa_0) sigma(y/t1),
    D_X U = -kappa X n_{s,r}/2,
    D_X log(phi) = -kappa p_{1,r}/2,

then keeps kappa=kappa_0 through X=100.  Between 100 and 110 it first shuts
off the axial prescription on a short logarithmic interval, then interpolates
``a`` to 0.8 while keeping ``D_X U=0``; the final interval has
``a=0.8, l=0.6, D_X U=0``.  The source leaves kappa_0 and the two final widths
as existence choices.  This module makes explicit autonomous finite choices
for them and does not call those choices paper-exact.

This increment exposes the actual profile boundary values required by PA.16:

    ell_i(eta) = log(C E(X_i,eta)),   G_i(eta) = U(X_i,eta),
    E = sqrt(2X) F.

The incoming five prefix moments and the source T_sh lower-bound certificate
remain separate dependencies; therefore this class is not yet a complete
inner-to-outer join or a global leading velocity field.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from .kokuno_reference_continuation import (
    KokunoReferenceContinuationCandidate,
    source_smooth_step,
)

SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-appendix-b-boundary-v1"
X_I = 110.0
X_FINAL_START = 100.0
A_FINAL = 0.8

_SOURCE_FORMULAS = {
    "activation": (
        "e_a=(1-kappa_0)sigma(y/t1), kappa=1-e_a; "
        "a=kappa*p1_r; D_X U=-kappa*X*n_s_r/2; "
        "D_X log(phi)=-kappa*p1_r/2"
    ),
    "reference_primitives": (
        "D_X p1_r=X*S_q_r/L-l_r*p1_r; D_X n_s_r+n_s_r=S_n_r/L"
    ),
    "final_join": (
        "from X=100 shut off the axial prescription on a short log interval; "
        "then D_X U=0 and interpolate a to 0.8; final a=0.8,l=0.6"
    ),
    "boundary": "X_i=110; ell_i=log(C*E(X_i,eta)); G_i=U(X_i,eta)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "appendix_B_Xi_profile_executable": True,
    "source_reference_primitive_equations_executable": True,
    "source_activation_equations_executable": True,
    "autonomous_kappa0_and_final_widths": True,
    "source_hidden_numeric_choices_recovered": False,
    "actual_upstream_appendix_B_boundary_values_executable_for_selected_realization": True,
    "incoming_five_moment_discrepancy_bound": False,
    "source_T_sh_lower_bound_verified": False,
    "inner_to_outer_join_completed": False,
    "global_pressure_matched": False,
    "global_leading_profile_reconstructed": False,
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


@dataclass(frozen=True)
class _ReferenceEtaSlice:
    eta: float
    X2: float
    F2: float
    U2: float
    F_eta2: float
    U_eta2: float
    Pi2: float
    Pi_eta2: float
    M2: float
    M_eta2: float


@dataclass(frozen=True)
class AppendixBBoundaryValue:
    eta: float
    ell_i: float
    G_i: float
    F_i: float
    E_i: float
    p1_reference_i: float
    n_s_reference_i: float
    radial_log_slope_F_i: float
    radial_log_slope_U_i: float


@dataclass(frozen=True)
class KokunoAppendixBBoundary:
    """Selected executable realization of the public Appendix-B exit profile."""

    reference: KokunoReferenceContinuationCandidate = field(
        default_factory=KokunoReferenceContinuationCandidate
    )
    kappa0: float = 0.01
    axial_shutdown_log_width: float = 0.02
    angular_settle_log_width: float = 0.02
    rtol: float = 2.0e-9
    atol: float = 2.0e-11
    max_step: float = 0.04
    eta_derivative_step: float = 2.0e-4

    def __post_init__(self) -> None:
        if not isinstance(self.reference, KokunoReferenceContinuationCandidate):
            raise TypeError("reference must be a KokunoReferenceContinuationCandidate")
        k0 = float(self.kappa0)
        w1 = float(self.axial_shutdown_log_width)
        w2 = float(self.angular_settle_log_width)
        rtol = float(self.rtol)
        atol = float(self.atol)
        max_step = float(self.max_step)
        eta_step = float(self.eta_derivative_step)
        if not math.isfinite(k0) or not 0.0 < k0 < 0.5:
            raise ValueError("kappa0 must satisfy 0<kappa0<1/2")
        available = math.log(X_I / X_FINAL_START)
        if not all(math.isfinite(v) and v > 0.0 for v in (w1, w2)):
            raise ValueError("final logarithmic widths must be positive and finite")
        if w1 + w2 >= available:
            raise ValueError("final logarithmic widths must fit strictly between X=100 and X_i=110")
        if not (0.0 < rtol <= 1.0e-6 and 0.0 < atol <= 1.0e-8):
            raise ValueError("ODE tolerances are outside the guarded finite-realization range")
        if not math.isfinite(max_step) or not 0.0 < max_step <= 0.1:
            raise ValueError("max_step must lie in (0,0.1]")
        if not math.isfinite(eta_step) or not 1.0e-6 <= eta_step <= 1.0e-2:
            raise ValueError("eta_derivative_step must lie in [1e-6,1e-2]")
        object.__setattr__(self, "kappa0", k0)
        object.__setattr__(self, "axial_shutdown_log_width", w1)
        object.__setattr__(self, "angular_settle_log_width", w2)
        object.__setattr__(self, "rtol", rtol)
        object.__setattr__(self, "atol", atol)
        object.__setattr__(self, "max_step", max_step)
        object.__setattr__(self, "eta_derivative_step", eta_step)

    @property
    def X0(self) -> float:
        return self.reference.X0

    @property
    def y100(self) -> float:
        return math.log(X_FINAL_START / self.X0)

    @property
    def y_i(self) -> float:
        return math.log(X_I / self.X0)

    @property
    def y_axial_end(self) -> float:
        return self.y100 + self.axial_shutdown_log_width

    @property
    def y_angular_end(self) -> float:
        return self.y_axial_end + self.angular_settle_log_width

    def geometry_report(self) -> dict[str, Any]:
        return {
            "X0": self.X0,
            "X_reference_freeze": self.reference.X2,
            "X_final_start": X_FINAL_START,
            "X_i": X_I,
            "kappa0": self.kappa0,
            "axial_shutdown_log_width": self.axial_shutdown_log_width,
            "angular_settle_log_width": self.angular_settle_log_width,
            "final_constant_a_log_width": self.y_i - self.y_angular_end,
            "source_free_choices_are_autonomous": True,
        }

    def _pi_eta(self, X: float, eta: float) -> float:
        h = self.eta_derivative_step
        if eta <= -1.0 + 2.0 * h:
            f0 = float(self.reference.Pi(X, eta))
            f1 = float(self.reference.Pi(X, eta + h))
            f2 = float(self.reference.Pi(X, eta + 2.0 * h))
            return (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h)
        if eta >= 1.0 - 2.0 * h:
            f0 = float(self.reference.Pi(X, eta))
            f1 = float(self.reference.Pi(X, eta - h))
            f2 = float(self.reference.Pi(X, eta - 2.0 * h))
            return (3.0 * f0 - 4.0 * f1 + f2) / (2.0 * h)
        return (
            float(self.reference.Pi(X, eta + h))
            - float(self.reference.Pi(X, eta - h))
        ) / (2.0 * h)

    def _reference_slice(self, eta: float) -> _ReferenceEtaSlice:
        eta = float(eta)
        X2 = float(self.reference.X2)
        F2 = float(self.reference.F(X2, eta))
        U2 = float(self.reference.U(X2, eta))
        F_eta2 = float(self.reference.F_eta(X2, eta))
        U_eta2 = float(self.reference.U_eta(X2, eta))
        Pi2 = float(self.reference.Pi(X2, eta))
        Pi_eta2 = self._pi_eta(X2, eta)
        avg2, avg_eta2 = self.reference._average_U_and_eta(X2, eta)
        return _ReferenceEtaSlice(
            eta=eta,
            X2=X2,
            F2=F2,
            U2=U2,
            F_eta2=F_eta2,
            U_eta2=U_eta2,
            Pi2=Pi2,
            Pi_eta2=Pi_eta2,
            M2=X2 * float(avg2),
            M_eta2=X2 * float(avg_eta2),
        )

    def _reference_state(self, X: float, eta: float, sl: _ReferenceEtaSlice) -> dict[str, float]:
        h = float(self.reference.h)
        A = 0.5 + h
        D = 0.5 - h
        d = 1.0 - eta * eta
        L = 1.0 - 2.0 * h * eta * eta
        if X >= sl.X2:
            F = sl.F2
            U = sl.U2
            F_X = 0.0
            U_X = 0.0
            F_eta = sl.F_eta2
            U_eta = sl.U_eta2
            Pi = sl.Pi2 + (X - sl.X2) * F * F
            Pi_eta = sl.Pi_eta2 + (X - sl.X2) * 2.0 * F * F_eta
            M = sl.M2 + (X - sl.X2) * U
            M_eta = sl.M_eta2 + (X - sl.X2) * U_eta
        else:
            F = float(self.reference.F(X, eta))
            U = float(self.reference.U(X, eta))
            F_X = float(self.reference.F_X(X, eta))
            U_X = float(self.reference.U_X(X, eta))
            F_eta = float(self.reference.F_eta(X, eta))
            U_eta = float(self.reference.U_eta(X, eta))
            Pi = float(self.reference.Pi(X, eta))
            Pi_eta = self._pi_eta(X, eta)
            avg, avg_eta = self.reference._average_U_and_eta(X, eta)
            M = X * float(avg)
            M_eta = X * float(avg_eta)
        if not F > 0.0:
            raise RuntimeError("Appendix-B reference requires positive F")
        W = 1.0 - (2.0 * D * eta * M + d * M_eta) / X
        Hc = D * eta + d * U
        l = 1.0 + X * F_X / F
        logE_eta = F_eta / F
        S_q = -W * l - h * (1.0 - 2.0 * eta * U) - Hc * logE_eta
        S_n = (
            -W * X * U_X
            - A * (1.0 - 2.0 * eta * U) * U
            - Hc * U_eta
            - d * Pi_eta
            + 4.0 * A * eta * Pi
            + 2.0 * eta * X * F * F
        )
        return {
            "F": F,
            "U": U,
            "F_X": F_X,
            "U_X": U_X,
            "l": l,
            "S_q": S_q,
            "S_n": S_n,
            "L": L,
        }

    def _rhs(self, y: float, state: np.ndarray, eta: float, sl: _ReferenceEtaSlice) -> np.ndarray:
        X = self.X0 * math.exp(float(y))
        ref = self._reference_state(X, eta, sl)
        p1, n_s = float(state[0]), float(state[1])
        dp1 = X * ref["S_q"] / ref["L"] - ref["l"] * p1
        dn_s = ref["S_n"] / ref["L"] - n_s

        if y < self.y100:
            if y < self.reference.log_transition_width:
                s = y / self.reference.log_transition_width
                kappa = 1.0 - (1.0 - self.kappa0) * float(source_smooth_step(np.asarray(s)))
            else:
                kappa = self.kappa0
            a = kappa * p1
            dlogF = -0.5 * a
            dU = -0.5 * kappa * X * n_s
        elif y < self.y_axial_end:
            s = (y - self.y100) / self.axial_shutdown_log_width
            beta = 1.0 - float(source_smooth_step(np.asarray(s)))
            a = self.kappa0 * p1
            dlogF = -0.5 * a
            dU = -0.5 * self.kappa0 * beta * X * n_s
        elif y < self.y_angular_end:
            s = (y - self.y_axial_end) / self.angular_settle_log_width
            sigma = float(source_smooth_step(np.asarray(s)))
            a = (1.0 - sigma) * self.kappa0 * p1 + sigma * A_FINAL
            dlogF = -0.5 * a
            dU = 0.0
        else:
            dlogF = -0.5 * A_FINAL
            dU = 0.0
        return np.asarray([dp1, dn_s, dlogF, dU], dtype=float)

    def _solve_eta(self, eta: float) -> AppendixBBoundaryValue:
        eta = float(eta)
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        sl = self._reference_slice(eta)
        F0 = float(self.reference.F(self.X0, eta))
        U0 = float(self.reference.U(self.X0, eta))
        F_X0 = float(self.reference.F_X(self.X0, eta))
        U_X0 = float(self.reference.U_X(self.X0, eta))
        p10 = -2.0 * self.X0 * F_X0 / F0
        n0 = -2.0 * U_X0
        initial = np.asarray([p10, n0, math.log(F0), U0], dtype=float)
        solved = solve_ivp(
            lambda yy, zz: self._rhs(yy, zz, eta, sl),
            (0.0, self.y_i),
            initial,
            method="DOP853",
            rtol=self.rtol,
            atol=self.atol,
            max_step=self.max_step,
        )
        if not solved.success or solved.y.size == 0:
            raise RuntimeError(f"Appendix-B boundary integration failed: {solved.message}")
        p1_i, n_i, logF_i, G_i = (float(v) for v in solved.y[:, -1])
        F_i = math.exp(logF_i)
        E_i = math.sqrt(2.0 * X_I) * F_i
        ell_i = math.log(float(self.reference.C) * E_i)
        values = np.asarray([ell_i, G_i, F_i, E_i, p1_i, n_i])
        if not np.all(np.isfinite(values)) or F_i <= 0.0 or E_i <= 0.0:
            raise RuntimeError("Appendix-B boundary integration produced non-finite data")
        return AppendixBBoundaryValue(
            eta=eta,
            ell_i=ell_i,
            G_i=G_i,
            F_i=F_i,
            E_i=E_i,
            p1_reference_i=p1_i,
            n_s_reference_i=n_i,
            radial_log_slope_F_i=-0.5 * A_FINAL,
            radial_log_slope_U_i=0.0,
        )

    def boundary_values(self, eta: Any) -> dict[str, np.ndarray]:
        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        out = {name: np.empty_like(eta_array, dtype=float) for name in (
            "ell_i", "G_i", "F_i", "E_i", "p1_reference_i", "n_s_reference_i"
        )}
        for index, value in enumerate(eta_array.reshape(-1)):
            item = self._solve_eta(float(value))
            for name in out:
                out[name].reshape(-1)[index] = float(getattr(item, name))
        return out

    def boundary_eta_derivatives(self, eta: Any) -> dict[str, np.ndarray]:
        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        step = self.eta_derivative_step
        d_ell = np.empty_like(eta_array, dtype=float)
        d_G = np.empty_like(eta_array, dtype=float)
        for index, value in enumerate(eta_array.reshape(-1)):
            e = float(value)
            if e <= -1.0 + step:
                a, b = self._solve_eta(e), self._solve_eta(e + step)
                de = step
            elif e >= 1.0 - step:
                a, b = self._solve_eta(e - step), self._solve_eta(e)
                de = step
            else:
                a, b = self._solve_eta(e - step), self._solve_eta(e + step)
                de = 2.0 * step
            d_ell.reshape(-1)[index] = (b.ell_i - a.ell_i) / de
            d_G.reshape(-1)[index] = (b.G_i - a.G_i) / de
        return {"ell_i_eta": d_ell, "G_i_eta": d_G}

    def to_payload(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "parameters": {
                "kappa0": self.kappa0,
                "axial_shutdown_log_width": self.axial_shutdown_log_width,
                "angular_settle_log_width": self.angular_settle_log_width,
                "rtol": self.rtol,
                "atol": self.atol,
                "max_step": self.max_step,
                "eta_derivative_step": self.eta_derivative_step,
                "origin": "autonomous_source_compatible_existence_choices_not_hidden_parameters",
            },
            "formula_map": dict(_SOURCE_FORMULAS),
            "geometry": self.geometry_report(),
            "dependencies": {"reference_sha256": self.reference.sha256},
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

    @property
    def sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_payload()).encode("utf-8")).hexdigest()

    def save_json(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        payload = self.to_payload()
        payload["sha256"] = self.sha256
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return target

    @classmethod
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoAppendixBBoundary":
        if not isinstance(payload, dict):
            raise ValueError("Appendix-B boundary payload must be a JSON object")
        unsigned = {k: v for k, v in payload.items() if k != "sha256"}
        if payload.get("schema") != SCHEMA:
            raise ValueError("unsupported Appendix-B boundary schema")
        if payload.get("formula_map") != _SOURCE_FORMULAS:
            raise ValueError("Appendix-B source formula map changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("Appendix-B truth-boundary metadata changed")
        if "sha256" in payload:
            expected = hashlib.sha256(_canonical_json(unsigned).encode("utf-8")).hexdigest()
            if payload["sha256"] != expected:
                raise ValueError("Appendix-B boundary payload hash mismatch")
        params = dict(payload["parameters"])
        params.pop("origin", None)
        obj = cls(**params)
        expected_payload = obj.to_payload()
        if unsigned != expected_payload:
            raise ValueError("Appendix-B boundary payload dependencies or geometry changed")
        return obj

    @classmethod
    def load_json(cls, path: str | Path) -> "KokunoAppendixBBoundary":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
