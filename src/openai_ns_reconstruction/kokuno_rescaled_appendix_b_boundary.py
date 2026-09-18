"""Log-stable selected-pressure Appendix-B trajectory through ``X_i=110``.

This module carries :class:`KokunoSourceRescaledReferenceContinuation` through
the public Appendix-B stress-activation/final-join equations in KokunoYumeto's
corrected 2026-09-09 reconstruction.  It is deliberately a *selected numerical
realization*: the public source leaves ``kappa_0`` and two short final
logarithmic widths as existence choices, while the upstream source fixed point
and complex-domain contraction certificate are still unresolved.

The implementation keeps ``log(F)`` primary.  At the selected Appendix-A
pressure scale, ordinary binary64 ``F`` legitimately underflows away from the
real-axis phase maximum; ratios such as ``F_eta/F`` are therefore evaluated as
``partial_eta log(F)`` rather than by dividing two underflowed numbers.  The
selected path also propagates

    M(X,eta)  = integral_0^X U(s,eta) ds,
    Pi_X      = F^2,

so a callable native-coordinate velocity can retain the incompressibility
prefix memory instead of dropping it at the activation stage.

This is not the source fixed point, not hidden-parameter recovery, not a
completed PA.16 join, and not an independent Navier--Stokes validation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.integrate import solve_ivp

from .kokuno_reference_continuation import source_smooth_step
from .kokuno_rescaled_reference_continuation import (
    KokunoSourceRescaledReferenceContinuation,
)
from .kokuno_similarity_coordinates import KokunoNativeSimilarityCoordinates


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-source-rescaled-appendix-b-boundary-v1"

X_I = 110.0
X_FINAL_START = 100.0
A_FINAL = 0.8

_SOURCE_FORMULAS = {
    "reference_primitives": (
        "D_X p1_r=X*S_q,r/L-l_r*p1_r; "
        "D_X n_s,r+n_s,r=S_n,r/L"
    ),
    "activation": (
        "e_a=(1-kappa_0)sigma(y/t1), kappa=1-e_a; "
        "a=kappa*p1_r; D_X U=-kappa*X*n_s,r/2; "
        "D_X log(phi)=-kappa*p1_r/2"
    ),
    "prefix_pressure": "D_X M=X*U; D_X Pi=X*F^2 in y=log(X/X0)",
    "final_join": (
        "from X=100 shut off the axial prescription on a short log interval; "
        "then keep D_X U=0 and interpolate a to 0.8; "
        "final a=0.8, l=1-a/2=0.6"
    ),
    "boundary": (
        "X_i=110; ell_i=log(C*E_i)=log(C)+log(E_i); G_i=U_i"
    ),
    "incompressibility": (
        "v0=(2 eta U-2D eta M/X-(1-eta^2)M_eta/X)/(1-2h eta^2)"
    ),
    "cartesian_velocity": (
        "u1=x*v0/(2q)-y*q^(-1-h)*F; "
        "u2=y*v0/(2q)+x*q^(-1-h)*F; u3=q^(-1/2-h)*U"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "selected_pressure_appendix_B_activation_executable": True,
    "selected_pressure_appendix_B_Xi_profile_executable": True,
    "selected_pressure_appendix_B_velocity_executable": True,
    "log_amplitude_trajectory_executable": True,
    "prefix_M_and_pressure_identity_carried": True,
    "autonomous_kappa0_and_final_widths": True,
    "continued_eta_derivatives_are_numerical": True,
    "selected_real_axis_C_normalization_is_autonomous": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_fixed_point_solved": False,
    "source_contraction_threshold_verified": False,
    "source_complex_C_bound_verified": False,
    "actual_source_appendix_B_trajectory_reconstructed": False,
    "actual_source_incoming_five_moments_bound": False,
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
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _exp_log(log_value: float, name: str) -> float:
    value = float(log_value)
    if not math.isfinite(value):
        raise OverflowError(f"{name} logarithm became non-finite")
    log_max = math.log(np.finfo(float).max)
    log_tiny = math.log(np.finfo(float).tiny)
    if value > log_max:
        raise OverflowError(f"{name} exceeds binary64 range")
    if value < log_tiny:
        return 0.0
    return math.exp(value)


def _square_from_log(log_value: float) -> float:
    return _exp_log(2.0 * float(log_value), "F^2")


@dataclass(frozen=True)
class _ReferenceSlice:
    eta: float
    X2: float
    log_F2: float
    log_F_eta2: float
    F2_sq: float
    U2: float
    U_eta2: float
    Pi2: float
    Pi_eta2: float
    M2: float
    M_eta2: float


@dataclass(frozen=True)
class AppendixBSelectedBoundaryValue:
    eta: float
    ell_i: float
    G_i: float
    log_F_i: float
    F_i: float
    E_i: float
    M_i: float
    Pi_i: float
    p1_reference_i: float
    n_s_reference_i: float
    radial_log_slope_F_i: float
    radial_log_slope_U_i: float


@dataclass(frozen=True)
class KokunoSourceRescaledAppendixBBoundary:
    """Selected-pressure executable Appendix-B activation through ``X_i``."""

    reference: KokunoSourceRescaledReferenceContinuation = field(
        default_factory=KokunoSourceRescaledReferenceContinuation
    )
    kappa0: float = 0.01
    axial_shutdown_log_width: float = 0.02
    angular_settle_log_width: float = 0.02
    rtol: float = 2.0e-8
    atol: float = 2.0e-10
    max_step: float = 0.25
    eta_derivative_step: float = 2.0e-4

    _coordinates: KokunoNativeSimilarityCoordinates = field(
        init=False, repr=False, compare=False
    )
    _slice_cache: dict[float, _ReferenceSlice] = field(
        init=False, repr=False, compare=False
    )
    _state_cache: dict[tuple[float, float], np.ndarray] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.reference, KokunoSourceRescaledReferenceContinuation):
            raise TypeError(
                "reference must be a KokunoSourceRescaledReferenceContinuation"
            )
        k0 = float(self.kappa0)
        w1 = float(self.axial_shutdown_log_width)
        w2 = float(self.angular_settle_log_width)
        rtol = float(self.rtol)
        atol = float(self.atol)
        max_step = float(self.max_step)
        eta_step = float(self.eta_derivative_step)
        if not math.isfinite(k0) or not 0.0 < k0 < 0.5:
            raise ValueError("kappa0 must satisfy 0<kappa0<1/2")
        if not all(math.isfinite(v) and v > 0.0 for v in (w1, w2)):
            raise ValueError("final logarithmic widths must be positive and finite")
        if w1 + w2 >= math.log(X_I / X_FINAL_START):
            raise ValueError(
                "final logarithmic widths must fit strictly between X=100 and X_i=110"
            )
        if not (0.0 < rtol <= 1.0e-6 and 0.0 < atol <= 1.0e-8):
            raise ValueError(
                "ODE tolerances are outside the guarded finite-realization range"
            )
        if not math.isfinite(max_step) or not 0.02 <= max_step <= 0.5:
            raise ValueError("max_step must lie in [0.02,0.5]")
        if not math.isfinite(eta_step) or not 1.0e-6 <= eta_step <= 1.0e-2:
            raise ValueError("eta_derivative_step must lie in [1e-6,1e-2]")
        object.__setattr__(self, "kappa0", k0)
        object.__setattr__(self, "axial_shutdown_log_width", w1)
        object.__setattr__(self, "angular_settle_log_width", w2)
        object.__setattr__(self, "rtol", rtol)
        object.__setattr__(self, "atol", atol)
        object.__setattr__(self, "max_step", max_step)
        object.__setattr__(self, "eta_derivative_step", eta_step)
        object.__setattr__(
            self, "_coordinates", KokunoNativeSimilarityCoordinates(h=self.h)
        )
        object.__setattr__(self, "_slice_cache", {})
        object.__setattr__(self, "_state_cache", {})

    @property
    def h(self) -> float:
        return float(self.reference.h)

    @property
    def D(self) -> float:
        return 0.5 - self.h

    @property
    def X0(self) -> float:
        return float(self.reference.X0)

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

    @property
    def log_C(self) -> float:
        return float(self.reference.seed.log_C_real_axis)

    def geometry_report(self) -> dict[str, Any]:
        return {
            "X0": self.X0,
            "X_reference_freeze": self.reference.X2,
            "X_final_start": X_FINAL_START,
            "X_i": X_I,
            "log_X0": math.log(self.X0),
            "log_X_i": math.log(X_I),
            "source_y_i": self.y_i,
            "kappa0": self.kappa0,
            "axial_shutdown_log_width": self.axial_shutdown_log_width,
            "angular_settle_log_width": self.angular_settle_log_width,
            "final_constant_a_log_width": self.y_i - self.y_angular_end,
            "source_free_choices_are_autonomous": True,
        }

    def _eta_fd(self, fn: Callable[[float], float], eta: float) -> float:
        h = self.eta_derivative_step
        if eta <= -1.0 + 2.0 * h:
            f0, f1, f2 = fn(eta), fn(eta + h), fn(eta + 2.0 * h)
            return (-3.0 * f0 + 4.0 * f1 - f2) / (2.0 * h)
        if eta >= 1.0 - 2.0 * h:
            f0, f1, f2 = fn(eta), fn(eta - h), fn(eta - 2.0 * h)
            return (3.0 * f0 - 4.0 * f1 + f2) / (2.0 * h)
        return (fn(eta + h) - fn(eta - h)) / (2.0 * h)

    def _make_slice(self, eta: float) -> _ReferenceSlice:
        eta = float(eta)
        cached = self._slice_cache.get(eta)
        if cached is not None:
            return cached
        X2 = float(self.reference.X2)
        values = self.reference.profile_values(X2, eta)
        scalar = {key: float(np.asarray(value)) for key, value in values.items()}
        log_F_eta2 = self._eta_fd(
            lambda ee: float(np.asarray(
                self.reference.profile_values(X2, ee)["log_F"]
            )),
            eta,
        )
        Pi_eta2 = self._eta_fd(
            lambda ee: float(np.asarray(self.reference.Pi(X2, ee))), eta
        )
        M_eta2 = self._eta_fd(
            lambda ee: float(np.asarray(
                self.reference.profile_values(X2, ee)["M"]
            )),
            eta,
        )
        item = _ReferenceSlice(
            eta=eta,
            X2=X2,
            log_F2=scalar["log_F"],
            log_F_eta2=log_F_eta2,
            F2_sq=_square_from_log(scalar["log_F"]),
            U2=scalar["U"],
            U_eta2=scalar["U_eta"],
            Pi2=scalar["Pi"],
            Pi_eta2=Pi_eta2,
            M2=scalar["M"],
            M_eta2=M_eta2,
        )
        self._slice_cache[eta] = item
        return item

    def _reference_state(
        self, X: float, eta: float, sl: _ReferenceSlice
    ) -> dict[str, float]:
        X = float(X)
        eta = float(eta)
        d = 1.0 - eta * eta
        L = 1.0 - 2.0 * self.h * eta * eta
        A = 0.5 + self.h
        H_coeff = self.D * eta

        if X >= sl.X2:
            log_F = sl.log_F2
            F_sq = sl.F2_sq
            U = sl.U2
            U_eta = sl.U_eta2
            D_X_log_F = 0.0
            D_X_U = 0.0
            log_F_eta = sl.log_F_eta2
            Pi = sl.Pi2 + (X - sl.X2) * F_sq
            Pi_eta = sl.Pi_eta2 + (X - sl.X2) * 2.0 * F_sq * log_F_eta
            M = sl.M2 + (X - sl.X2) * U
            M_eta = sl.M_eta2 + (X - sl.X2) * U_eta
        else:
            values = self.reference.profile_values(X, eta)
            log_F = float(np.asarray(values["log_F"]))
            F_sq = _square_from_log(log_F)
            U = float(np.asarray(values["U"]))
            U_eta = float(np.asarray(values["U_eta"]))
            D_X_log_F = float(np.asarray(values["D_X_log_F"]))
            D_X_U = float(np.asarray(values["D_X_U"]))
            Pi = float(np.asarray(values["Pi"]))
            M = float(np.asarray(values["M"]))
            log_F_eta = self._eta_fd(
                lambda ee: float(np.asarray(
                    self.reference.profile_values(X, ee)["log_F"]
                )),
                eta,
            )
            Pi_eta = self._eta_fd(
                lambda ee: float(np.asarray(self.reference.Pi(X, ee))), eta
            )
            M_eta = self._eta_fd(
                lambda ee: float(np.asarray(
                    self.reference.profile_values(X, ee)["M"]
                )),
                eta,
            )

        W = 1.0 - (2.0 * self.D * eta * M + d * M_eta) / X
        Hc = H_coeff + d * U
        l = 1.0 + D_X_log_F
        S_q = -W * l - self.h * (1.0 - 2.0 * eta * U) - Hc * log_F_eta
        S_n = (
            -W * D_X_U
            - A * (1.0 - 2.0 * eta * U) * U
            - Hc * U_eta
            - d * Pi_eta
            + 4.0 * A * eta * Pi
            + 2.0 * eta * X * F_sq
        )
        out = {
            "log_F": log_F,
            "F_sq": F_sq,
            "U": U,
            "U_eta": U_eta,
            "D_X_log_F": D_X_log_F,
            "D_X_U": D_X_U,
            "log_F_eta": log_F_eta,
            "Pi": Pi,
            "Pi_eta": Pi_eta,
            "M": M,
            "M_eta": M_eta,
            "l": l,
            "S_q": S_q,
            "S_n": S_n,
            "L": L,
        }
        if not all(math.isfinite(v) for v in out.values()):
            raise OverflowError("reference primitive state became non-finite")
        return out

    def _initial_state(self, eta: float) -> np.ndarray:
        values = self.reference.profile_values(self.X0, eta)
        log_F = float(np.asarray(values["log_F"]))
        U = float(np.asarray(values["U"]))
        M = float(np.asarray(values["M"]))
        Pi = float(np.asarray(values["Pi"]))
        D_X_log_F = float(np.asarray(values["D_X_log_F"]))
        D_X_U = float(np.asarray(values["D_X_U"]))
        p1 = -2.0 * D_X_log_F
        n_s = -2.0 * D_X_U / self.X0
        state = np.asarray([p1, n_s, log_F, U, M, Pi], dtype=float)
        if np.any(~np.isfinite(state)):
            raise OverflowError("selected Appendix-B initial state became non-finite")
        return state

    def _radial_controls(
        self, y: float, p1: float, n_s: float
    ) -> tuple[float, float]:
        if y < self.y100:
            if y < self.reference.log_transition_width:
                s = y / self.reference.log_transition_width
                kappa = 1.0 - (1.0 - self.kappa0) * float(
                    source_smooth_step(np.asarray(s))
                )
            else:
                kappa = self.kappa0
            return -0.5 * kappa * p1, -0.5 * kappa * self._X_from_y(y) * n_s

        if y < self.y_axial_end:
            s = (y - self.y100) / self.axial_shutdown_log_width
            beta = 1.0 - float(source_smooth_step(np.asarray(s)))
            return -0.5 * self.kappa0 * p1, -0.5 * self.kappa0 * beta * self._X_from_y(y) * n_s

        if y < self.y_angular_end:
            s = (y - self.y_axial_end) / self.angular_settle_log_width
            sigma = float(source_smooth_step(np.asarray(s)))
            a = (1.0 - sigma) * self.kappa0 * p1 + sigma * A_FINAL
            return -0.5 * a, 0.0

        return -0.5 * A_FINAL, 0.0

    def _X_from_y(self, y: float) -> float:
        value = self.X0 * math.exp(float(y))
        if not math.isfinite(value) or value <= 0.0:
            raise OverflowError("selected Appendix-B radial coordinate became invalid")
        return value

    def _rhs(
        self, y: float, state: np.ndarray, eta: float, sl: _ReferenceSlice
    ) -> np.ndarray:
        X = self._X_from_y(y)
        ref = self._reference_state(X, eta, sl)
        p1, n_s, log_F, U, _, _ = (float(v) for v in state)
        dp1 = X * ref["S_q"] / ref["L"] - ref["l"] * p1
        dn_s = ref["S_n"] / ref["L"] - n_s
        dlogF, dU = self._radial_controls(y, p1, n_s)
        F_sq = _square_from_log(log_F)
        dM = X * U
        dPi = X * F_sq
        out = np.asarray([dp1, dn_s, dlogF, dU, dM, dPi], dtype=float)
        if np.any(~np.isfinite(out)):
            raise OverflowError("selected Appendix-B ODE derivative became non-finite")
        return out

    def _solve_state(self, y: float, eta: float) -> np.ndarray:
        y = float(y)
        eta = float(eta)
        if not math.isfinite(y) or not 0.0 <= y <= self.y_i:
            raise ValueError("source y must lie in [0,y_i]")
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        key = (y, eta)
        cached = self._state_cache.get(key)
        if cached is not None:
            return cached.copy()

        initial = self._initial_state(eta)
        if y == 0.0:
            self._state_cache[key] = initial.copy()
            return initial

        sl = self._make_slice(eta)
        solved = solve_ivp(
            lambda yy, zz: self._rhs(yy, zz, eta, sl),
            (0.0, y),
            initial,
            method="DOP853",
            rtol=self.rtol,
            atol=self.atol,
            max_step=self.max_step,
        )
        if not solved.success or solved.y.size == 0:
            raise RuntimeError(
                f"selected-pressure Appendix-B integration failed: {solved.message}"
            )
        final = np.asarray(solved.y[:, -1], dtype=float)
        if np.any(~np.isfinite(final)):
            raise OverflowError(
                "selected-pressure Appendix-B integration produced non-finite data"
            )
        self._state_cache[key] = final.copy()
        return final

    def _selected_scalar(self, X: float, eta: float) -> dict[str, float]:
        X = float(X)
        eta = float(eta)
        if not math.isfinite(X) or not 0.0 <= X <= X_I:
            raise ValueError("selected Appendix-B profile requires 0<=X<=110")
        if not math.isfinite(eta) or not -1.0 <= eta <= 1.0:
            raise ValueError("eta must lie in [-1,1]")
        if X <= self.X0:
            values = self.reference.profile_values(X, eta)
            out = {key: float(np.asarray(value)) for key, value in values.items()}
            log_F_eta = self._eta_fd(
                lambda ee: float(np.asarray(
                    self.reference.profile_values(X, ee)["log_F"]
                )),
                eta,
            ) if X > 0.0 else float(
                self.reference.seed.rescaling_lambda
                * np.asarray(self.reference.seed.axis_state(np.asarray(eta))["zeta_star"])
            )
            M_eta = self._eta_fd(
                lambda ee: float(np.asarray(
                    self.reference.profile_values(X, ee)["M"]
                )),
                eta,
            ) if X > 0.0 else 0.0
            out["log_F_eta"] = log_F_eta
            out["M_eta"] = M_eta
            out["p1_reference"] = -2.0 * out["D_X_log_F"]
            if X > 0.0:
                out["n_s_reference"] = -2.0 * out["D_X_U"] / X
            else:
                axis = self.reference.seed.axis_state(np.asarray(eta))
                B = -float(np.asarray(axis["Z_star"])) / (
                    2.0 * float(np.asarray(axis["L"]))
                )
                out["n_s_reference"] = -2.0 * B
            out["source_a"] = -2.0 * out["D_X_log_F"]
            out["source_l"] = 1.0 + out["D_X_log_F"]
            return out

        y = math.log(X / self.X0)
        state = self._solve_state(y, eta)
        p1, n_s, log_F, U, M, Pi = (float(v) for v in state)
        sl = self._make_slice(eta)
        rhs = self._rhs(y, state, eta, sl)
        D_X_log_F = float(rhs[2])
        D_X_U = float(rhs[3])
        log_F_eta = self._eta_fd(
            lambda ee: float(self._solve_state(y, ee)[2]), eta
        )
        U_eta = self._eta_fd(
            lambda ee: float(self._solve_state(y, ee)[3]), eta
        )
        M_eta = self._eta_fd(
            lambda ee: float(self._solve_state(y, ee)[4]), eta
        )
        F = _exp_log(log_F, "F")
        F_sq = _square_from_log(log_F)
        E = _exp_log(0.5 * math.log(2.0 * X) + log_F, "E")
        D_X_F = F * D_X_log_F
        D_X_E = 0.5 * E + math.sqrt(2.0 * X) * D_X_F
        d = 1.0 - eta * eta
        L = 1.0 - 2.0 * self.h * eta * eta
        v0 = (
            2.0 * eta * U
            - 2.0 * self.D * eta * M / X
            - d * M_eta / X
        ) / L
        out = {
            "log_F": log_F,
            "F": F,
            "E": E,
            "U": U,
            "v0": v0,
            "Pi": Pi,
            "Pi_X": F_sq,
            "M": M,
            "M_eta": M_eta,
            "D_X_log_F": D_X_log_F,
            "D_X_F": D_X_F,
            "D_X_E": D_X_E,
            "D_X_U": D_X_U,
            "log_F_eta": log_F_eta,
            "F_eta": F * log_F_eta,
            "U_eta": U_eta,
            "p1_reference": p1,
            "n_s_reference": n_s,
            "source_a": -2.0 * D_X_log_F,
            "source_l": 1.0 + D_X_log_F,
        }
        if not all(math.isfinite(value) for value in out.values()):
            raise OverflowError("selected Appendix-B profile produced non-finite values")
        return out

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        X_array, eta_array = np.broadcast_arrays(
            _finite_array(X, "X"), _finite_array(eta, "eta")
        )
        if np.any((X_array < 0.0) | (X_array > X_I)):
            raise ValueError("selected Appendix-B profile requires 0<=X<=110")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        keys = (
            "log_F", "F", "E", "U", "v0", "Pi", "Pi_X", "M", "M_eta",
            "D_X_log_F", "D_X_F", "D_X_E", "D_X_U", "log_F_eta",
            "F_eta", "U_eta", "p1_reference", "n_s_reference",
            "source_a", "source_l",
        )
        out = {key: np.empty_like(X_array, dtype=float) for key in keys}
        for index, (xx, ee) in enumerate(
            zip(X_array.reshape(-1), eta_array.reshape(-1), strict=True)
        ):
            values = self._selected_scalar(float(xx), float(ee))
            for key in keys:
                out[key].reshape(-1)[index] = values[key]
        if any(np.any(~np.isfinite(value)) for value in out.values()):
            raise OverflowError("selected Appendix-B batch profile became non-finite")
        return out

    def profile_values_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        log_array, eta_array = np.broadcast_arrays(
            _finite_array(log_X, "log_X"), _finite_array(eta, "eta")
        )
        if np.any(log_array > math.log(X_I)):
            raise ValueError("selected Appendix-B log_X exceeds log(110)")
        with np.errstate(under="ignore"):
            X = np.exp(log_array)
        return self.profile_values(X, eta_array)

    def F(self, X: Any, eta: Any) -> np.ndarray:
        return self.profile_values(X, eta)["F"]

    def U(self, X: Any, eta: Any) -> np.ndarray:
        return self.profile_values(X, eta)["U"]

    def v0(self, X: Any, eta: Any) -> np.ndarray:
        return self.profile_values(X, eta)["v0"]

    def Pi(self, X: Any, eta: Any) -> np.ndarray:
        return self.profile_values(X, eta)["Pi"]

    def boundary_values(self, eta: Any) -> dict[str, np.ndarray]:
        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        names = (
            "ell_i", "G_i", "log_F_i", "F_i", "E_i", "M_i", "Pi_i",
            "p1_reference_i", "n_s_reference_i",
            "radial_log_slope_F_i", "radial_log_slope_U_i",
        )
        out = {name: np.empty_like(eta_array, dtype=float) for name in names}
        for index, value in enumerate(eta_array.reshape(-1)):
            e = float(value)
            state = self._solve_state(self.y_i, e)
            p1, n_s, log_F, U, M, Pi = (float(v) for v in state)
            rhs = self._rhs(self.y_i, state, e, self._make_slice(e))
            log_E = 0.5 * math.log(2.0 * X_I) + log_F
            item = AppendixBSelectedBoundaryValue(
                eta=e,
                ell_i=self.log_C + log_E,
                G_i=U,
                log_F_i=log_F,
                F_i=_exp_log(log_F, "F"),
                E_i=_exp_log(log_E, "E"),
                M_i=M,
                Pi_i=Pi,
                p1_reference_i=p1,
                n_s_reference_i=n_s,
                radial_log_slope_F_i=float(rhs[2]),
                radial_log_slope_U_i=float(rhs[3]),
            )
            for name in names:
                out[name].reshape(-1)[index] = float(getattr(item, name))
        return out

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        x_array, y_array, z_array, t_array = np.broadcast_arrays(
            _finite_array(x, "x"), _finite_array(y, "y"),
            _finite_array(z, "z"), _finite_array(t, "t")
        )
        coordinates = self._coordinates.evaluate(x_array, y_array, z_array, t_array)
        profiles = self.profile_values(coordinates["X"], coordinates["eta"])
        q = np.asarray(coordinates["q"], dtype=float)
        swirl_factor = np.power(q, -1.0 - self.h) * profiles["F"]
        radial_factor = profiles["v0"] / (2.0 * q)
        axial_factor = np.power(q, -0.5 - self.h) * profiles["U"]
        result = np.stack(
            (
                radial_factor * x_array - swirl_factor * y_array,
                radial_factor * y_array + swirl_factor * x_array,
                axial_factor,
            ),
            axis=-1,
        )
        if np.any(~np.isfinite(result)):
            raise OverflowError("selected Appendix-B velocity became non-finite")
        return result

    __call__ = velocity

    def report(self) -> dict[str, Any]:
        eta = float(self.reference.seed.phase_stationary_eta)
        boundary = self.boundary_values(np.asarray([eta, 0.0]))
        return {
            "pressure_scale": self.reference.pressure_scale,
            "rescaling_lambda": self.reference.rescaling_lambda,
            "log_C_real_axis": self.log_C,
            "geometry": self.geometry_report(),
            "sample_eta": [eta, 0.0],
            "sample_ell_i": boundary["ell_i"].tolist(),
            "sample_G_i": boundary["G_i"].tolist(),
            "sample_log_F_i": boundary["log_F_i"].tolist(),
            "all_boundary_outputs_finite": bool(
                all(np.all(np.isfinite(value)) for value in boundary.values())
            ),
            "truth_boundary": dict(_TRUTH_BOUNDARY),
        }

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
            "reference": self.reference.to_payload(),
            "parameters": {
                "kappa0": self.kappa0,
                "axial_shutdown_log_width": self.axial_shutdown_log_width,
                "angular_settle_log_width": self.angular_settle_log_width,
                "rtol": self.rtol,
                "atol": self.atol,
                "max_step": self.max_step,
                "eta_derivative_step": self.eta_derivative_step,
            },
            "autonomous_numerics": {
                "kappa0_and_final_widths": (
                    "explicit source-compatible existence choices, not hidden parameters"
                ),
                "ODE": "scipy DOP853 in source y=log(X/X0)",
                "eta_derivatives": (
                    "second-order centered/one-sided finite differences on the selected path"
                ),
                "amplitude": (
                    "log(F) primary; binary64 F/E may underflow while log values remain finite"
                ),
                "C_normalization": (
                    "selected real-axis phase maximum inherited from the rescaled core seed"
                ),
            },
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
    def from_payload(
        cls, payload: dict[str, Any]
    ) -> "KokunoSourceRescaledAppendixBBoundary":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected source-rescaled Appendix-B boundary schema")
        if payload.get("source_formulas") != _SOURCE_FORMULAS:
            raise ValueError("source-rescaled Appendix-B formulas changed")
        if payload.get("truth_boundary") != _TRUTH_BOUNDARY:
            raise ValueError("source-rescaled Appendix-B truth-boundary metadata changed")
        claimed = payload.get("sha256")
        unsigned = {key: value for key, value in payload.items() if key != "sha256"}
        expected = hashlib.sha256(
            _canonical_json(unsigned).encode("utf-8")
        ).hexdigest()
        if claimed != expected:
            raise ValueError("source-rescaled Appendix-B payload SHA-256 mismatch")
        params = payload.get("parameters")
        if not isinstance(params, dict):
            raise ValueError("missing source-rescaled Appendix-B parameters")
        obj = cls(
            reference=KokunoSourceRescaledReferenceContinuation.from_payload(
                payload.get("reference")
            ),
            kappa0=float(params["kappa0"]),
            axial_shutdown_log_width=float(params["axial_shutdown_log_width"]),
            angular_settle_log_width=float(params["angular_settle_log_width"]),
            rtol=float(params["rtol"]),
            atol=float(params["atol"]),
            max_step=float(params["max_step"]),
            eta_derivative_step=float(params["eta_derivative_step"]),
        )
        if obj.sha256 != claimed:
            raise ValueError("source-rescaled Appendix-B replay changed SHA-256")
        return obj

    def save_json(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return output

    @classmethod
    def load_json(
        cls, path: str | Path
    ) -> "KokunoSourceRescaledAppendixBBoundary":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
