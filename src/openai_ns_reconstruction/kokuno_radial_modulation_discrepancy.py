"""Executable finite-frequency radial insertion and five-moment discrepancy.

The corrected Kokuno reconstruction (2026-09-09, Appendix C) inserts an
admissible shear loop into the unchanged logarithmic radial coordinate through
periodic zero-mean antiderivatives ``A`` and ``B``:

    E_N = E exp(A(X,eta,N log X)/N),
    U_N = U + B(X,eta,N log X)/N.

With phase held fixed in slow derivatives, the source identities are

    a_N = a_L - 2 D_X A/N,
    b_N = exp(-A/N) [ b_L + 2 D_X B/(N E) ],

where ``a=1-2 D_X log E`` and ``b=2 D_X U/E``.  The insertion changes the
five prefix moments

    M   = int U dx,
    I   = int sqrt(2x) E dx,
    J   = int U sqrt(2x) E dx,
    S   = int (U^2-E^2/2) dx,
    C_p = int E^2/(2x) dx.

Those discrepancies are the input to the later PA.17 I1 repair already
implemented by :mod:`kokuno_shear_moment_repair`.

This module makes that source insertion contract executable on the existing
RF40 power-law backbone.  The source proves existence of a particular
admissible loop but does not publish one hidden numerical loop realization.
For an executable diagnostic this module therefore uses an explicitly
AUTONOMOUS smooth periodic antiderivative pair on a compact subinterval before
I1:

    A = alpha chi(log X) cos(2 pi phi),
    B = beta E chi(log X) sin(2 pi phi),       phi in R/Z.

It has the source-required zero phase means and compact radial support, so it
exercises the exact insertion and moment formulas.  It is *not* claimed to be
the recovered source admissible cone loop.  Consequently the generated PA.17
target is an executable source-form diagnostic target, not the paper-exact
target and not independent PDE validation.

The discrepancy is evaluated directly in log-X and normalized by the I1 scales
``X_1,e_1`` used by ``KokunoShearMomentRepair``.  This avoids forming enormous
absolute powers of X on the existence-style outer schedule.  Row order is the
PA.17 repair order ``(M,J,I,S,C_p)``.
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

from .kokuno_outer_base_schedule import KokunoOuterBaseSchedule
from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_shear_moment_repair import KokunoShearMomentRepair, ShearMomentSolveResult


SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
CORRECTED_RELEASE = "zenodo:22678406"
CORRECTED_RELEASE_DATE = "2026-09-09"
SCHEMA = "kokuno-radial-modulation-discrepancy-v1"
MOMENT_ORDER = ("M", "J", "I", "S", "C_p")

# The source chooses a compact modulation interval before the first reserved
# repair patch but does not expose one hidden numerical interval.  These
# offsets are autonomous and lie strictly inside the already executable RF40
# power stage, with a gap before I1.
_DEFAULT_SUPPORT_LEFT_OFFSET = -35.3
_DEFAULT_SUPPORT_RIGHT_OFFSET = -29.1

_SOURCE_FORMULAS = {
    "insertion": "E_N=E exp(A(X,eta,N log X)/N); U_N=U+B(X,eta,N log X)/N",
    "phase_chain_rule": "X d_X = D_X + N d_phi after phi=N log X",
    "shear_a": "a_N=a_L-2 D_X A/N",
    "shear_b": "b_N=exp(-A/N)*(b_L+2 D_X B/(N E))",
    "moments": (
        "M=int U dx; I=int sqrt(2x)E dx; J=int U sqrt(2x)E dx; "
        "S=int(U^2-E^2/2)dx; C_p=int E^2/(2x)dx"
    ),
    "exact_increment_integrands": (
        "u; sqrt(2X)e; H0*u+U0*sqrt(2X)e+sqrt(2X)u*e; "
        "2U0*u+u^2-E0*e-e^2/2; (E0*e+e^2/2)/X"
    ),
    "source_bound": "all five fixed-eta parameter-derivative discrepancies are O_m(N^-1)",
    "repair_dependency": "negative pre-I1 discrepancy is restored on reserved I1 by PA.17",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "source_finite_frequency_insertion_formula_executable": True,
    "source_shear_chain_rule_executable": True,
    "source_five_prefix_moment_discrepancy_executable": True,
    "stable_I1_normalized_discrepancy_executable": True,
    "PA17_target_from_autonomous_insertion_executable": True,
    "autonomous_antiderivative_realization": True,
    "autonomous_modulation_interval": True,
    "source_hidden_loop_parameters_recovered": False,
    "source_admissible_loop_reconstructed": False,
    "actual_source_admissible_loop_discrepancy_supplied": False,
    "eta_smooth_repair_coefficient_family_reconstructed": False,
    "cone_modulation_completed": False,
    "I1_cone_repair_applied_to_actual_modulation": False,
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


def _cinfty_window(log_X: np.ndarray, left: float, right: float) -> np.ndarray:
    values = np.asarray(log_X, dtype=float)
    centre = 0.5 * (left + right)
    half = 0.5 * (right - left)
    s = (values - centre) / half
    out = np.zeros_like(values)
    active = np.abs(s) < 1.0
    if np.any(active):
        sa = s[active]
        out[active] = np.exp(1.0 - 1.0 / (1.0 - sa * sa))
    return out


def _cinfty_window_dlogX(log_X: np.ndarray, left: float, right: float) -> np.ndarray:
    values = np.asarray(log_X, dtype=float)
    centre = 0.5 * (left + right)
    half = 0.5 * (right - left)
    s = (values - centre) / half
    out = np.zeros_like(values)
    active = np.abs(s) < 1.0
    if np.any(active):
        sa = s[active]
        bump = np.exp(1.0 - 1.0 / (1.0 - sa * sa))
        out[active] = bump * (-2.0 * sa / (1.0 - sa * sa) ** 2) / half
    return out


@dataclass(frozen=True)
class KokunoRadialModulationDiscrepancy:
    """Source-form radial modulation with an autonomous periodic loop probe."""

    outer_schedule: KokunoOuterReservedPatchSchedule = field(
        default_factory=KokunoOuterReservedPatchSchedule
    )
    frequency: int = 16
    alpha: float = 5.0e-4
    beta: float = 5.0e-4
    support_left_offset: float = _DEFAULT_SUPPORT_LEFT_OFFSET
    support_right_offset: float = _DEFAULT_SUPPORT_RIGHT_OFFSET
    quadrature_order: int = 8
    panels_per_period: int = 3

    _outer_base: KokunoOuterBaseSchedule = field(init=False, repr=False, compare=False)
    _repair: KokunoShearMomentRepair = field(init=False, repr=False, compare=False)
    _nodes: np.ndarray = field(init=False, repr=False, compare=False)
    _weights: np.ndarray = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be a KokunoOuterReservedPatchSchedule")
        if isinstance(self.frequency, bool) or not isinstance(self.frequency, (int, np.integer)):
            raise TypeError("frequency must be an integer")
        frequency = int(self.frequency)
        if not 4 <= frequency <= 128:
            raise ValueError("frequency must lie in [4,128]")
        alpha = float(self.alpha)
        beta = float(self.beta)
        if not math.isfinite(alpha) or not (0.0 < abs(alpha) <= 1.0e-2):
            raise ValueError("alpha must be finite, nonzero, and |alpha|<=1e-2")
        if not math.isfinite(beta) or not (0.0 < abs(beta) <= 1.0e-2):
            raise ValueError("beta must be finite, nonzero, and |beta|<=1e-2")
        left_offset = float(self.support_left_offset)
        right_offset = float(self.support_right_offset)
        if not math.isfinite(left_offset) or not math.isfinite(right_offset):
            raise ValueError("support offsets must be finite")
        if not left_offset < right_offset:
            raise ValueError("support_left_offset must be smaller than support_right_offset")
        if isinstance(self.quadrature_order, bool) or not isinstance(
            self.quadrature_order, (int, np.integer)
        ):
            raise TypeError("quadrature_order must be an integer")
        order = int(self.quadrature_order)
        if not 4 <= order <= 16:
            raise ValueError("quadrature_order must lie in [4,16]")
        if isinstance(self.panels_per_period, bool) or not isinstance(
            self.panels_per_period, (int, np.integer)
        ):
            raise TypeError("panels_per_period must be an integer")
        panels_per_period = int(self.panels_per_period)
        if not 2 <= panels_per_period <= 8:
            raise ValueError("panels_per_period must lie in [2,8]")

        base = KokunoOuterBaseSchedule(outer_schedule=self.outer_schedule)
        repair = KokunoShearMomentRepair(outer_schedule=self.outer_schedule)
        support_left = base.log_X_end + left_offset
        support_right = base.log_X_end + right_offset
        i1_left, _ = self.outer_schedule.reserved_log_intervals()["I1"]
        if support_left <= base.log_X_w:
            raise ValueError("autonomous modulation support must lie inside the RF40 power stage")
        if support_right >= i1_left:
            raise ValueError("autonomous modulation support must end strictly before source I1")
        stage_probe = base.stage_for_log_X(
            np.asarray([support_left + 1.0e-9, support_right - 1.0e-9])
        )
        if np.any(stage_probe != "power_law"):
            raise ValueError("autonomous modulation support must remain in the power-law stage")

        nodes, weights = leggauss(order)
        nodes = np.asarray(nodes, dtype=float)
        weights = np.asarray(weights, dtype=float)
        nodes.setflags(write=False)
        weights.setflags(write=False)

        object.__setattr__(self, "frequency", frequency)
        object.__setattr__(self, "alpha", alpha)
        object.__setattr__(self, "beta", beta)
        object.__setattr__(self, "support_left_offset", left_offset)
        object.__setattr__(self, "support_right_offset", right_offset)
        object.__setattr__(self, "quadrature_order", order)
        object.__setattr__(self, "panels_per_period", panels_per_period)
        object.__setattr__(self, "_outer_base", base)
        object.__setattr__(self, "_repair", repair)
        object.__setattr__(self, "_nodes", nodes)
        object.__setattr__(self, "_weights", weights)

    @property
    def support_log_interval(self) -> tuple[float, float]:
        return (
            self._outer_base.log_X_end + self.support_left_offset,
            self._outer_base.log_X_end + self.support_right_offset,
        )

    @property
    def lambda_outer(self) -> float:
        return float(self.outer_schedule.lambda_outer)

    @property
    def base_a(self) -> float:
        # On RF40 power law E ~ X^(-1/2-lambda).
        return 2.0 + 2.0 * self.lambda_outer

    def source_interval_report(self) -> dict[str, Any]:
        left, right = self.support_log_interval
        i1 = self.outer_schedule.reserved_log_intervals()["I1"]
        return {
            "RF40_power_log_interval": [self._outer_base.log_X_w, self._outer_base.log_X_end],
            "autonomous_modulation_support_log_X": [left, right],
            "source_I1_log_interval": list(i1),
            "gap_to_I1_log_units": float(i1[0] - right),
            "support_is_autonomous": True,
            "source_admissible_loop_reconstructed": False,
        }

    @staticmethod
    def _phase(log_X: np.ndarray, frequency: int) -> np.ndarray:
        # The source antiderivatives are 1-periodic in phi=N log X.
        return np.remainder(float(frequency) * np.asarray(log_X, dtype=float), 1.0)

    def _window(self, log_X: np.ndarray) -> np.ndarray:
        left, right = self.support_log_interval
        return _cinfty_window(log_X, left, right)

    def _window_dlogX(self, log_X: np.ndarray) -> np.ndarray:
        left, right = self.support_log_interval
        return _cinfty_window_dlogX(log_X, left, right)

    def antiderivatives(self, log_X: Any, eta: Any, *, phase: Any | None = None) -> dict[str, np.ndarray]:
        """Return the autonomous A,B realization and phase derivatives.

        ``B_over_E`` is returned instead of materializing B separately because
        the source identity uses ``D_X B/E`` and the existence-style E can be
        extremely small.  The actual B is ``E * B_over_E``.
        """

        log_array, eta_array = np.broadcast_arrays(
            _finite_array(log_X, "log_X"), _finite_array(eta, "eta")
        )
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        if phase is None:
            phase_array = self._phase(log_array, self.frequency)
        else:
            phase_array = np.broadcast_to(_finite_array(phase, "phase"), log_array.shape)
        window = self._window(log_array)
        window_y = self._window_dlogX(log_array)
        angle = 2.0 * math.pi * phase_array
        sine = np.sin(angle)
        cosine = np.cos(angle)
        A = self.alpha * window * cosine
        B_over_E = self.beta * window * sine
        A_phi = -2.0 * math.pi * self.alpha * window * sine
        B_over_E_phi = 2.0 * math.pi * self.beta * window * cosine
        A_slow = self.alpha * window_y * cosine
        # D_X(B)/E with phase held fixed.  On the power law D_X E/E=-p.
        p = 0.5 + self.lambda_outer
        D_X_B_over_E = self.beta * (window_y - p * window) * sine
        return {
            "phase": phase_array,
            "window": window,
            "window_dlogX": window_y,
            "A": A,
            "A_phi": A_phi,
            "D_X_A_slow": A_slow,
            "B_over_E": B_over_E,
            "B_over_E_phi": B_over_E_phi,
            "D_X_B_over_E": D_X_B_over_E,
        }

    def loop_shear(self, log_X: Any, eta: Any, *, phase: Any | None = None) -> dict[str, np.ndarray]:
        """Return base and autonomous loop shear before finite-N insertion."""

        anti = self.antiderivatives(log_X, eta, phase=phase)
        a = np.full(np.shape(anti["A"]), self.base_a, dtype=float)
        b = np.zeros_like(a)
        a_L = a - 2.0 * anti["A_phi"]
        b_L = b + 2.0 * anti["B_over_E_phi"]
        return {"a": a, "b": b, "a_L": a_L, "b_L": b_L, **anti}

    def profile_values_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate exact source-form E_N,U_N and analytic first derivatives.

        This method is restricted to the RF40 power stage.  It exposes the
        profile contract needed for later global assembly, but deliberately
        does not fabricate the incompressibility-derived radial component: the
        eta-smooth PA.17 repaired family has not yet been assembled.
        """

        log_array, eta_array = np.broadcast_arrays(
            _finite_array(log_X, "log_X"), _finite_array(eta, "eta")
        )
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        stages = self._outer_base.stage_for_log_X(log_array)
        if np.any(stages != "power_law"):
            raise ValueError("radial modulation profile is restricted to the RF40 power-law stage")

        base = self._outer_base.profile_values_logX(log_array, eta_array)
        anti = self.antiderivatives(log_array, eta_array)
        inv_N = 1.0 / float(self.frequency)
        exp_factor = np.exp(anti["A"] * inv_N)
        E0 = np.asarray(base["E"], dtype=float)
        F0 = np.asarray(base["F"], dtype=float)
        D_X_E0 = np.asarray(base["D_X_E"], dtype=float)
        D_X_F0 = np.asarray(base["D_X_F"], dtype=float)
        E0_eta = np.asarray(base["E_eta"], dtype=float)
        F0_eta = np.asarray(base["F_eta"], dtype=float)

        E = E0 * exp_factor
        F = F0 * exp_factor
        U = E0 * anti["B_over_E"] * inv_N
        # Total log-radial derivative after phi=N log X.  For A/N the fast
        # derivative contributes A_phi exactly, while the slow derivative
        # retains 1/N.
        d_A_over_N = anti["A_phi"] + anti["D_X_A_slow"] * inv_N
        D_X_E = exp_factor * (D_X_E0 + E0 * d_A_over_N)
        D_X_F = exp_factor * (D_X_F0 + F0 * d_A_over_N)
        # B=E0*(B/E0).  Fast phase differentiation contributes N*(B/E)_phi.
        D_X_U = (
            D_X_E0 * anti["B_over_E"] * inv_N
            + E0 * anti["D_X_B_over_E"] * inv_N
            + E0 * anti["B_over_E_phi"]
            - D_X_E0 * anti["B_over_E"] * inv_N
        )
        # The previous line intentionally writes the slow B derivative through
        # the stored D_X B / E quantity; simplify to the exact source identity.
        D_X_U = E0 * (
            anti["B_over_E_phi"] + anti["D_X_B_over_E"] * inv_N
        )
        ratio_eta = np.divide(E0_eta, E0, out=np.zeros_like(E0_eta), where=E0 != 0.0)
        E_eta = E * ratio_eta
        F_eta = F * np.divide(F0_eta, F0, out=np.zeros_like(F0_eta), where=F0 != 0.0)
        U_eta = U * ratio_eta

        a_L = self.base_a - 2.0 * anti["A_phi"]
        b_L = 2.0 * anti["B_over_E_phi"]
        a_N = a_L - 2.0 * anti["D_X_A_slow"] * inv_N
        b_N = np.exp(-anti["A"] * inv_N) * (
            b_L + 2.0 * anti["D_X_B_over_E"] * inv_N
        )
        a_direct = 1.0 - 2.0 * np.divide(D_X_E, E)
        b_direct = 2.0 * np.divide(D_X_U, E)

        return {
            "stage": stages,
            "log_X": log_array,
            "phase": anti["phase"],
            "window": anti["window"],
            "A": anti["A"],
            "B_over_E0": anti["B_over_E"],
            "E0": E0,
            "F0": F0,
            "U0": np.asarray(base["U"], dtype=float),
            "E": E,
            "D_X_E": D_X_E,
            "E_eta": E_eta,
            "F": F,
            "D_X_F": D_X_F,
            "F_eta": F_eta,
            "U": U,
            "D_X_U": D_X_U,
            "U_eta": U_eta,
            "Pi_X": F * F,
            "a_L": a_L,
            "b_L": b_L,
            "a_N": a_N,
            "b_N": b_N,
            "a_direct": a_direct,
            "b_direct": b_direct,
        }

    def _quadrature(self) -> tuple[np.ndarray, np.ndarray]:
        left, right = self.support_log_interval
        width = right - left
        panels = max(
            8,
            int(math.ceil(width * float(self.frequency) * float(self.panels_per_period))),
        )
        edges = np.linspace(left, right, panels + 1)
        points: list[np.ndarray] = []
        weights: list[np.ndarray] = []
        for a, b in zip(edges[:-1], edges[1:]):
            half = 0.5 * (b - a)
            centre = 0.5 * (b + a)
            points.append(centre + half * self._nodes)
            weights.append(half * self._weights)
        return np.concatenate(points), np.concatenate(weights)

    def normalized_discrepancy(self, eta: Any) -> np.ndarray:
        """Return exact five-moment changes normalized to PA.17 I1 scales.

        Output final axis is ``(M,J,I,S,C_p)``.  Integration is performed in
        ``y=log X`` with ``dX=X dy``.  The power-law relation to the I1 scale is
        used analytically, so no enormous absolute X powers are materialized.
        """

        eta_array = _finite_array(eta, "eta")
        if np.any(np.abs(eta_array) > 1.0):
            raise ValueError("eta must lie in [-1,1]")
        original_shape = eta_array.shape
        eta_flat = eta_array.reshape(-1)
        y, weights = self._quadrature()
        anti = self.antiderivatives(y, np.zeros_like(y))
        phase = np.asarray(anti["phase"], dtype=float)
        window = np.asarray(anti["window"], dtype=float)
        angle = 2.0 * math.pi * phase
        sine = np.sin(angle)
        cosine = np.cos(angle)
        A = self.alpha * window * cosine
        expm1_factor = np.expm1(A / float(self.frequency))

        log_X_1 = float(self._repair.log_X_1)
        r = np.exp(y - log_X_1)
        p = 0.5 + self.lambda_outer
        f = np.asarray(self.outer_schedule.source_f(eta_flat), dtype=float)
        E_hat = f[:, None] * np.power(r[None, :], -p)
        # B=beta*E0*window*sin, so u=(U_N-U0)/e1 is below.
        u_hat = (
            self.beta
            * E_hat
            * window[None, :]
            * sine[None, :]
            / float(self.frequency)
        )
        e_hat = E_hat * expm1_factor[None, :]
        root2 = math.sqrt(2.0)
        r1 = r[None, :]
        r32 = np.power(r, 1.5)[None, :]

        # On this autonomous support the RF40 power stage has U0=0 exactly.
        integrand_M = r1 * u_hat
        integrand_J = root2 * r32 * (E_hat * u_hat + u_hat * e_hat)
        integrand_I = root2 * r32 * e_hat
        integrand_S = r1 * (u_hat * u_hat - E_hat * e_hat - 0.5 * e_hat * e_hat)
        integrand_Cp = E_hat * e_hat + 0.5 * e_hat * e_hat
        rows = np.stack(
            (integrand_M, integrand_J, integrand_I, integrand_S, integrand_Cp),
            axis=1,
        )
        values = np.tensordot(rows, weights, axes=([-1], [0]))
        return np.asarray(values, dtype=float).reshape(original_shape + (5,))

    def repair_target_normalized(self, eta: Any) -> np.ndarray:
        """Return the negative pre-I1 discrepancy in PA.17 row order."""

        return -self.normalized_discrepancy(eta)

    def solve_repair_at_eta(self, eta: float) -> ShearMomentSolveResult:
        """Route one scalar autonomous discrepancy into the existing PA.17 inverse.

        This is intentionally pointwise in eta.  It does not construct or claim
        the source-required smooth coefficient family ``c(eta)``.
        """

        eta_value = float(eta)
        if not math.isfinite(eta_value) or abs(eta_value) > 1.0:
            raise ValueError("eta must be finite and lie in [-1,1]")
        target = self.repair_target_normalized(eta_value)
        f_eta = float(np.asarray(self.outer_schedule.source_f(eta_value)).item())
        return self._repair.solve(target, f_eta=f_eta)

    def discrepancy_report(self, eta: float) -> dict[str, Any]:
        eta_value = float(eta)
        discrepancy = np.asarray(self.normalized_discrepancy(eta_value), dtype=float)
        target = -discrepancy
        solution = self.solve_repair_at_eta(eta_value)
        return {
            "eta": eta_value,
            "moment_order": list(MOMENT_ORDER),
            "normalized_discrepancy": [float(v) for v in discrepancy],
            "normalized_repair_target": [float(v) for v in target],
            "pointwise_repair_coefficients": [float(v) for v in solution.coefficients],
            "pointwise_repair_max_abs_residual": float(solution.max_abs_residual),
            "pointwise_repair_success": bool(solution.success),
            "eta_smooth_family_claimed": False,
            "source_admissible_loop_claimed": False,
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
            "moment_order": list(MOMENT_ORDER),
            "outer_schedule": self.outer_schedule.to_payload(),
            "frequency": self.frequency,
            "alpha": self.alpha,
            "beta": self.beta,
            "support_left_offset": self.support_left_offset,
            "support_right_offset": self.support_right_offset,
            "quadrature_order": self.quadrature_order,
            "panels_per_period": self.panels_per_period,
            "autonomous_realization": {
                "A": "alpha*chi(logX)*cos(2*pi*phi)",
                "B": "beta*E0*chi(logX)*sin(2*pi*phi)",
                "chi": "C-infinity compact bump on autonomous pre-I1 power-stage interval",
                "phase": "phi=N log X mod 1",
                "claim": "diagnostic source-form loop realization; not recovered source admissible loop",
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
    def from_payload(cls, payload: dict[str, Any]) -> "KokunoRadialModulationDiscrepancy":
        if not isinstance(payload, dict) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected Kokuno radial-modulation schema")
        schedule_payload = payload.get("outer_schedule")
        if not isinstance(schedule_payload, dict):
            raise ValueError("outer_schedule is missing")
        obj = cls(
            outer_schedule=KokunoOuterReservedPatchSchedule.from_payload(schedule_payload),
            frequency=int(payload.get("frequency")),
            alpha=float(payload.get("alpha")),
            beta=float(payload.get("beta")),
            support_left_offset=float(payload.get("support_left_offset")),
            support_right_offset=float(payload.get("support_right_offset")),
            quadrature_order=int(payload.get("quadrature_order")),
            panels_per_period=int(payload.get("panels_per_period")),
        )
        if obj.to_payload() != payload:
            raise ValueError("radial-modulation payload hash or content mismatch")
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
    def load_json(cls, path: str | Path) -> "KokunoRadialModulationDiscrepancy":
        return cls.from_payload(json.loads(Path(path).read_text(encoding="utf-8")))
