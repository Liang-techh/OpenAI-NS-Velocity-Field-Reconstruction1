"""Bind the exact current A1 I1 candidate through the corrected-source I2 heat repair.

This module is one narrow Kokuno Agent-1 composition increment.  Its parent is
``KokunoPA16CurrentCartesianI1FrozenGate`` (A1 PR #1051 exact head
``ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5``).  The parent already carries
the repository-autonomous finite-N insertion and PA.17 I1 repair on the current
PA.10 -> PA.16 -> RF40 lineage, with the I1 construction gate frozen at 5e-7.

The corrected public reconstruction later repairs the three heat-replacement
moment discrepancies on reserved interval I2 using

    delta E = e_* sum_i c_i(eta) beta_i(X/X_*),    delta U = 0,

with the three normalized rows ``(C_p, S, I_sub)``.  The repository already has
an arbitrary-precision realization of that local map in
``KokunoHighPrecisionHierarchicalHeatRepair``; Agent-4 PR #339 independently
audited its continuous moments with an implementation-distinct quadrature.
This adapter reuses that executable component rather than constructing another
heat solver.

The new seam is *current-lineage* composition.  The I1 exit generally leaves a
small candidate-side incompressibility-memory remainder.  I2 is pure swirl and
therefore does not change physical M=int U dX or M_eta.  We measure the actual
I1-exit difference from the current RF40 base and carry it as

    Delta(M/X)(X)     = (X_I1_end/X) Delta(M/X)_I1_end,
    Delta(M_eta/X)(X) = (X_I1_end/X) Delta(M_eta/X)_I1_end.

The heat correction is then added only to E (and F=E/sqrt(2X)); v0 is rebuilt
from the full current M/M_eta history.  The candidate fails closed at the I2
upper endpoint.  I3/I4, terminal/global assembly, matched pressure/forcing and
held-out Navier--Stokes validation remain separate tasks.

KokunoYumeto's corrected 2026-09-09 reconstruction is provenance only.  The
three concrete compact bumps and high-precision numerical backend remain
repository-autonomous choices; this module does not recover hidden OpenAI
parameters and does not claim a paper-exact field.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass, field
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_high_precision_heat_repair import (
    KokunoHighPrecisionHierarchicalHeatRepair,
    _BUMP_CENTERS as _I2_BUMP_CENTERS,
    _BUMP_HALF_WIDTH as _I2_BUMP_HALF_WIDTH,
)
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa16_current_cartesian_i1_frozen_gate import (
    I1_CLOSURE_TOLERANCE,
    KokunoPA16CurrentCartesianI1FrozenGate,
)

SCHEMA = "kokuno-pa16-current-cartesian-i2-heat-repair-v1"
PARENT_EXACT_HEAD = "ea59dc305b4265fcb0bb0f948c2a45f4477d1ae5"
I2_AUDIT_PR = 339
I2_AUDIT_HEAD = "abf4152b981fba12a493f1fc9db0495661995505"

_SOURCE_FORMULAS = {
    "I2_background": "U=0; E=e_* f(eta) (X/X_*)^(-1/2-lambda)",
    "I2_repair": "delta E=e_* sum_i c_i(eta) beta_i(X/X_*); delta U=0",
    "I2_rows": "normalized row order (C_p,S,I_sub)",
    "I2_target": (
        "target=-(Delta C_p/e_*^2, Delta S/(X_*e_*^2), "
        "Delta I_sub/(X_*^(3/2)e_*))"
    ),
    "F": "F=E/sqrt(2X)",
    "incompressibility_primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "post_I1_memory": (
        "I2 has delta U=0, so the physical I1-exit Delta M and Delta M_eta "
        "are constants and their normalized ratios decay as X_I1_end/X"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact #1051 frozen-gate current I1 candidate",
    "I2_solver": (
        "reuse KokunoHighPrecisionHierarchicalHeatRepair: arbitrary-precision "
        "continuous moment tensor with repository-autonomous compact bumps"
    ),
    "I2_independent_audit": (
        "Agent-4 #339 / abf4152b... independently re-integrated the continuous "
        "three-row moments; that local audit is not a PDE validation"
    ),
    "float_profile_view": (
        "the selected Decimal I2 correction and analytic bump derivative are "
        "materialized to float64 only at the candidate API boundary"
    ),
    "memory_handoff": (
        "measure actual total-minus-RF40 M/X and M_eta/X at exact current I1 exit; "
        "carry the corresponding physical memory through I2 without resetting it"
    ),
    "new_residual_tuning": "none",
    "domain_stop": "fail closed at reserved I2 upper endpoint before I3",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "exact_current_I1_frozen_gate_parent_consumed": True,
    "current_I1_exit_M_memory_carried_through_I2": True,
    "existing_high_precision_I2_heat_repair_consumed": True,
    "independent_local_I2_moment_audit_exists": True,
    "current_I2_overlay_applied": True,
    "cartesian_velocity_executable_through_I2": True,
    "velocity_interface_vectorized": True,
    "configuration_serializable": True,
    "source_hidden_loop_parameters_recovered": False,
    "source_admissible_loop_reconstructed": False,
    "paper_exact_I2_bump_parameters_recovered": False,
    "current_I3_overlay_applied": False,
    "current_I4_overlay_applied": False,
    "source_terminal_tail_schedule_bound_into_current_velocity": False,
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


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


def _D(value: float | str | Decimal) -> Decimal:
    if isinstance(value, Decimal):
        return value
    if isinstance(value, str):
        return Decimal(value)
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("cannot convert a non-finite value to Decimal")
    return Decimal(repr(number))


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianI2HeatRepair:
    """Current A1 Cartesian leading candidate through the reserved I2 patch."""

    parent: KokunoPA16CurrentCartesianI1FrozenGate = field(
        default_factory=KokunoPA16CurrentCartesianI1FrozenGate,
        repr=False,
        compare=False,
    )

    _i2: KokunoHighPrecisionHierarchicalHeatRepair = field(
        init=False, repr=False, compare=False
    )
    _solution_cache: dict[float, Any] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianI1FrozenGate):
            raise TypeError("parent must be KokunoPA16CurrentCartesianI1FrozenGate")
        if self.parent.closure_tolerance != I1_CLOSURE_TOLERANCE:
            raise ValueError("current I2 lineage requires the exact frozen 5e-7 I1 gate")

        i2 = KokunoHighPrecisionHierarchicalHeatRepair()
        schedule = self.parent.pre_i1.modulation.outer_schedule
        if i2.outer_schedule.to_payload() != schedule.to_payload():
            raise ValueError(
                "audited high-precision I2 repair and current lineage must use the same outer schedule"
            )

        i1_low, i1_high = schedule.reserved_log_intervals()["I1"]
        i2_low, i2_high = schedule.reserved_log_intervals()["I2"]
        i3_low, _ = schedule.reserved_log_intervals()["I3"]
        if not math.isclose(
            float(i1_high), self.parent.log_X_I1_end, rel_tol=0.0, abs_tol=2.0e-12
        ):
            raise ValueError("parent I1 endpoint and current outer schedule disagree")
        if not float(i1_high) <= float(i2_low) < float(i2_high) <= float(i3_low):
            raise ValueError("reserved I1/I2/I3 ordering is inconsistent")
        if float(i2_high) >= math.log(np.finfo(float).max):
            raise ValueError("current I2 endpoint is not materializable in float64 X")

        object.__setattr__(self, "_i2", i2)
        object.__setattr__(self, "_solution_cache", {})

    @property
    def i2(self) -> KokunoHighPrecisionHierarchicalHeatRepair:
        return self._i2

    @property
    def current(self):
        return self.parent.current

    @property
    def geometry(self):
        return self.parent.geometry

    @property
    def A(self) -> float:
        return float(self.parent.A)

    @property
    def D(self) -> float:
        return float(self.parent.D)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.parent.eta_interval)

    @property
    def log_X_I1_end(self) -> float:
        return float(self.parent.log_X_I1_end)

    @property
    def X_I1_end(self) -> float:
        return float(self.parent.X_I1_end)

    @property
    def log_X_I2_start(self) -> float:
        return float(self.i2.outer_schedule.reserved_log_intervals()["I2"][0])

    @property
    def log_X_I2_end(self) -> float:
        return float(self.i2.outer_schedule.reserved_log_intervals()["I2"][1])

    @property
    def X_I2_start(self) -> float:
        return float(math.exp(self.log_X_I2_start))

    @property
    def X_I2_end(self) -> float:
        return float(math.exp(self.log_X_I2_end))

    @property
    def log_X_I3_start(self) -> float:
        return float(self.i2.outer_schedule.reserved_log_intervals()["I3"][0])

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.parent.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_I2_end)):
            raise ValueError(
                "X must lie in the current I2-repair domain "
                f"[0,{self.X_I2_end:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _solution(self, eta: float):
        key = float(eta)
        cached = self._solution_cache.get(key)
        if cached is None:
            cached = self.i2.solve(key)
            self._solution_cache[key] = cached
        return cached

    @staticmethod
    def _decimal_bump_and_derivative(
        x_star: Decimal, center_text: str
    ) -> tuple[Decimal, Decimal]:
        center = Decimal(center_text)
        half_width = Decimal(_I2_BUMP_HALF_WIDTH)
        s = (x_star - center) / half_width
        if abs(s) >= 1:
            return Decimal(0), Decimal(0)
        one_minus = Decimal(1) - s * s
        beta = (Decimal(1) - Decimal(1) / one_minus).exp()
        beta_x = beta * (-Decimal(2) * s / (one_minus * one_minus)) / half_width
        return beta, beta_x

    def _i2_delta_decimal(
        self, X: float, eta: float
    ) -> tuple[Decimal, Decimal, Decimal, Decimal]:
        """Return (delta_E, delta_E_X, delta_F, delta_F_X) in Decimal.

        The source bump derivative is analytic.  Decimal is retained until the
        candidate API boundary because the audited I2 hierarchy spans hundreds
        of decades and the physical correction can be sub-ulp in float64.
        """

        X_value = float(X)
        eta_value = float(eta)
        if not (self.X_I2_start < X_value < self.X_I2_end):
            zero = Decimal(0)
            return zero, zero, zero, zero
        log_x_star = math.log(X_value) - float(self.i2.outer_schedule.log_X_star)
        x_star_float = math.exp(log_x_star)
        centers = tuple(float(value) for value in _I2_BUMP_CENTERS)
        half_width = float(_I2_BUMP_HALF_WIDTH)
        if not (
            min(centers) - half_width < x_star_float < max(centers) + half_width
        ):
            zero = Decimal(0)
            return zero, zero, zero, zero

        solution = self._solution(eta_value)
        with localcontext() as ctx:
            ctx.prec = int(solution.precision_digits)
            x_star = _D(x_star_float)
            correction = Decimal(0)
            correction_x = Decimal(0)
            for coefficient, center in zip(
                solution.coefficient_decimals(), _I2_BUMP_CENTERS
            ):
                beta, beta_x = self._decimal_bump_and_derivative(x_star, center)
                correction += coefficient * beta
                correction_x += coefficient * beta_x

            e_star = _D(float(self.i2.outer_schedule.log_e_star)).exp()
            inv_X_star = (-_D(float(self.i2.outer_schedule.log_X_star))).exp()
            delta_E = e_star * correction
            delta_E_X = e_star * inv_X_star * correction_x
            Xd = _D(X_value)
            root = (Decimal(2) * Xd).sqrt()
            delta_F = delta_E / root
            delta_F_X = delta_E_X / root - delta_E / (Decimal(2) * Xd * root)
            return delta_E, delta_E_X, delta_F, delta_F_X

    def _i1_exit_overlay_memory(
        self, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        X = np.full(unique.shape, self.X_I1_end, dtype=float)
        total = self.parent.similarity_profile_values(X, unique)
        base = self.current.similarity_profile_values(X, unique)
        delta = np.asarray(total["M_over_X_current_i1_repair"], dtype=float) - np.asarray(
            base["M_over_X_current_rf40_power_law"], dtype=float
        )
        delta_eta = np.asarray(
            total["M_eta_over_X_current_i1_repair"], dtype=float
        ) - np.asarray(base["M_eta_over_X_current_rf40_power_law"], dtype=float)
        return (
            delta[inverse].reshape(values.shape),
            delta_eta[inverse].reshape(values.shape),
        )

    def _carried_i1_memory(
        self, X: np.ndarray, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        dm, dm_eta = self._i1_exit_overlay_memory(eta)
        factor = self.X_I1_end / np.asarray(X, dtype=float)
        return factor * dm, factor * dm_eta

    def similarity_profile_values(
        self, X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)

        F = np.empty_like(xf)
        U = np.empty_like(xf)
        E = np.empty_like(xf)
        m_ratio = np.empty_like(xf)
        m_eta_ratio = np.empty_like(xf)
        i1_memory = np.zeros_like(xf)
        i1_memory_eta = np.zeros_like(xf)
        delta_E = np.zeros_like(xf)
        delta_F = np.zeros_like(xf)
        v0 = np.empty_like(xf)
        region = np.empty(xf.shape, dtype=object)

        inherited = xf <= self.X_I1_end
        if np.any(inherited):
            p = self.parent.similarity_profile_values(xf[inherited], ef[inherited])
            F[inherited] = np.asarray(p["F_current_i1_repair"], dtype=float)
            U[inherited] = np.asarray(p["U_current_i1_repair"], dtype=float)
            E[inherited] = np.asarray(p["E_current_i1_repair"], dtype=float)
            m_ratio[inherited] = np.asarray(
                p["M_over_X_current_i1_repair"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_i1_repair"], dtype=float
            )
            v0[inherited] = np.asarray(p["v0_current_i1_repair"], dtype=float)
            region[inherited] = np.asarray(p["region"], dtype=object)

        post_i1 = ~inherited
        if np.any(post_i1):
            Xo = xf[post_i1]
            etao = ef[post_i1]
            base = self.current.similarity_profile_values(Xo, etao)
            dm, dm_eta = self._carried_i1_memory(Xo, etao)

            base_F = np.asarray(base["F_current_rf40_power_law"], dtype=float)
            base_U = np.asarray(base["U_current_rf40_power_law"], dtype=float)
            base_E = np.asarray(base["E_current_rf40_power_law"], dtype=float)
            dE = np.zeros_like(Xo)
            dF = np.zeros_like(Xo)
            for index, (Xv, etav) in enumerate(zip(Xo, etao)):
                dE_dec, _, dF_dec, _ = self._i2_delta_decimal(float(Xv), float(etav))
                dE[index] = float(dE_dec)
                dF[index] = float(dF_dec)

            mo = np.asarray(base["M_over_X_current_rf40_power_law"], dtype=float) + dm
            meo = np.asarray(
                base["M_eta_over_X_current_rf40_power_law"], dtype=float
            ) + dm_eta
            Uo = base_U
            axis = self.geometry.physical_profiles.axis_profiles.values(
                np.zeros_like(etao), etao
            )
            L = np.asarray(axis["L"], dtype=float)
            d = np.asarray(axis["d"], dtype=float)
            v0o = (
                2.0 * etao * Uo
                - 2.0 * self.D * etao * mo
                - d * meo
            ) / L

            F[post_i1] = base_F + dF
            U[post_i1] = Uo
            E[post_i1] = base_E + dE
            m_ratio[post_i1] = mo
            m_eta_ratio[post_i1] = meo
            i1_memory[post_i1] = dm
            i1_memory_eta[post_i1] = dm_eta
            delta_E[post_i1] = dE
            delta_F[post_i1] = dF
            v0[post_i1] = v0o
            logX = np.log(Xo)
            region[post_i1] = np.where(
                (logX > self.log_X_I2_start) & (logX < self.log_X_I2_end),
                "current_high_precision_I2_heat_repair",
                "current_post_I1_memory_before_I2",
            )

        arrays = (F, U, E, m_ratio, m_eta_ratio, i1_memory, i1_memory_eta, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current I2 profile produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current I2 profile lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "F_current_i2_repair": F.reshape(shape),
            "U_current_i2_repair": U.reshape(shape),
            "E_current_i2_repair": E.reshape(shape),
            "M_over_X_current_i2_repair": m_ratio.reshape(shape),
            "M_eta_over_X_current_i2_repair": m_eta_ratio.reshape(shape),
            "delta_M_over_X_from_I1_exit": i1_memory.reshape(shape),
            "delta_M_eta_over_X_from_I1_exit": i1_memory_eta.reshape(shape),
            "delta_E_high_precision_I2_float_view": delta_E.reshape(shape),
            "delta_F_high_precision_I2_float_view": delta_F.reshape(shape),
            "v0_current_i2_repair": v0.reshape(shape),
        }

    def similarity_radial_derivatives(
        self, X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast_similarity(X, eta)
        if np.any(X_arr <= 0.0):
            raise ValueError("similarity_radial_derivatives requires X>0")
        shape = X_arr.shape
        xf = X_arr.reshape(-1)
        ef = eta_arr.reshape(-1)
        F_X = np.empty_like(xf)
        U_X = np.empty_like(xf)
        E_X = np.empty_like(xf)
        delta_F_X = np.zeros_like(xf)
        delta_E_X = np.zeros_like(xf)

        inherited = xf <= self.X_I1_end
        if np.any(inherited):
            d0 = self.parent.similarity_radial_derivatives(
                xf[inherited], ef[inherited]
            )
            F_X[inherited] = np.asarray(d0["F_current_i1_repair_X"], dtype=float)
            U_X[inherited] = np.asarray(d0["U_current_i1_repair_X"], dtype=float)
            E_X[inherited] = np.asarray(d0["E_current_i1_repair_X"], dtype=float)

        post_i1 = ~inherited
        if np.any(post_i1):
            Xo = xf[post_i1]
            etao = ef[post_i1]
            base = self.current.similarity_radial_derivatives(Xo, etao)
            dEX = np.zeros_like(Xo)
            dFX = np.zeros_like(Xo)
            for index, (Xv, etav) in enumerate(zip(Xo, etao)):
                _, dEX_dec, _, dFX_dec = self._i2_delta_decimal(float(Xv), float(etav))
                dEX[index] = float(dEX_dec)
                dFX[index] = float(dFX_dec)
            F_X[post_i1] = np.asarray(
                base["F_current_rf40_power_law_X"], dtype=float
            ) + dFX
            U_X[post_i1] = np.asarray(
                base["U_current_rf40_power_law_X"], dtype=float
            )
            E_X[post_i1] = np.asarray(
                base["E_current_rf40_power_law_X"], dtype=float
            ) + dEX
            delta_F_X[post_i1] = dFX
            delta_E_X[post_i1] = dEX

        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_current_i2_repair_X": F_X.reshape(shape),
            "U_current_i2_repair_X": U_X.reshape(shape),
            "E_current_i2_repair_X": E_X.reshape(shape),
            "delta_F_high_precision_I2_X_float_view": delta_F_X.reshape(shape),
            "delta_E_high_precision_I2_X_float_view": delta_E_X.reshape(shape),
        }

    def values(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        tol = 128.0 * np.finfo(float).eps * self.X_I2_end
        if np.any(X > self.X_I2_end + tol):
            raise ValueError("Cartesian point lies beyond the current I2 domain")
        X = np.minimum(X, self.X_I2_end)
        eta = np.asarray(coords["eta"], dtype=float)
        p = self.similarity_profile_values(X, eta)

        xb, yb, qb = np.broadcast_arrays(
            _finite(x, "x"), _finite(y, "y"), np.asarray(coords["q"], dtype=float)
        )
        F = np.asarray(p["F_current_i2_repair"], dtype=float)
        U = np.asarray(p["U_current_i2_repair"], dtype=float)
        v0 = np.asarray(p["v0_current_i2_repair"], dtype=float)
        radial = v0 / (2.0 * qb)
        swirl = qb ** (-self.A - 0.5) * F
        axial = qb ** (-self.A) * U
        u = radial * xb - swirl * yb
        v = radial * yb + swirl * xb
        w = axial
        return {
            **coords,
            **p,
            "u": np.asarray(u, dtype=float),
            "v": np.asarray(v, dtype=float),
            "w": np.asarray(w, dtype=float),
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        vals = self.values(x, y, z, t)
        return np.stack((vals["u"], vals["v"], vals["w"]), axis=-1)

    def i2_decimal_correction_report(self, X: float, eta: float) -> dict[str, Any]:
        X_value = float(X)
        eta_value = float(eta)
        if not self.X_I2_start < X_value < self.X_I2_end:
            raise ValueError("X must lie strictly inside reserved I2")
        dE, dEX, dF, dFX = self._i2_delta_decimal(X_value, eta_value)
        solution = self._solution(eta_value)
        return {
            "X": X_value,
            "eta": eta_value,
            "delta_E_decimal": str(dE),
            "delta_E_X_decimal": str(dEX),
            "delta_F_decimal": str(dF),
            "delta_F_X_decimal": str(dFX),
            "solution_sha256": solution.sha256,
            "precision_digits": solution.precision_digits,
            "max_relative_residual": solution.max_relative_residual,
            "coefficients": list(solution.coefficients),
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_i1_frozen_gate": self.parent.configuration(),
            "bound_high_precision_i2": {
                "sha256": self.i2.sha256,
                "payload": self.i2.to_payload(),
                "independent_audit_pr": I2_AUDIT_PR,
                "independent_audit_head": I2_AUDIT_HEAD,
                "mutable": False,
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianI2HeatRepair":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current I2 heat-repair schema")
        parent_payload = payload.get("parent_i1_frozen_gate")
        bound = payload.get("bound_high_precision_i2")
        if not isinstance(parent_payload, Mapping) or not isinstance(bound, Mapping):
            raise ValueError("missing current I2 parent/dependency configuration")
        if bound.get("mutable") is not False:
            raise ValueError("bound high-precision I2 dependency must be immutable")
        if bound.get("independent_audit_pr") != I2_AUDIT_PR:
            raise ValueError("unexpected I2 independent-audit PR")
        if bound.get("independent_audit_head") != I2_AUDIT_HEAD:
            raise ValueError("unexpected I2 independent-audit head")

        parent = KokunoPA16CurrentCartesianI1FrozenGate.from_configuration(parent_payload)
        obj = cls(parent=parent)
        if bound.get("sha256") != obj.i2.sha256:
            raise ValueError("serialized I2 semantic SHA differs from the bound audited component")
        if bound.get("payload") != obj.i2.to_payload():
            raise ValueError("serialized I2 payload differs from the bound audited component")
        return obj

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianI2HeatRepair":
        return cls.from_configuration(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "i2_sha256": self.i2.sha256,
            "i2_independent_audit": {
                "pr": I2_AUDIT_PR,
                "head": I2_AUDIT_HEAD,
            },
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "truth_boundary": _TRUTH_BOUNDARY,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta0 = np.asarray(0.0)
        exit_total = self.parent.similarity_profile_values(
            np.asarray(self.X_I1_end), eta0
        )
        exit_base = self.current.similarity_profile_values(
            np.asarray(self.X_I1_end), eta0
        )
        delta_m = float(
            np.asarray(exit_total["M_over_X_current_i1_repair"])
            - np.asarray(exit_base["M_over_X_current_rf40_power_law"])
        )
        delta_m_eta = float(
            np.asarray(exit_total["M_eta_over_X_current_i1_repair"])
            - np.asarray(exit_base["M_eta_over_X_current_rf40_power_law"])
        )
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "semantic_sha256": self.semantic_sha256,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "corrected_release": CORRECTED_RELEASE,
                "corrected_release_date": CORRECTED_RELEASE_DATE,
            },
            "intervals": {
                "log_X_I1_end": self.log_X_I1_end,
                "log_X_I2_start": self.log_X_I2_start,
                "log_X_I2_end": self.log_X_I2_end,
                "log_X_I3_start": self.log_X_I3_start,
            },
            "current_I1_exit_memory_eta0": {
                "delta_M_over_X": delta_m,
                "delta_M_eta_over_X": delta_m_eta,
                "physical_memory_carried_constant_through_I2": True,
            },
            "high_precision_I2": {
                "sha256": self.i2.sha256,
                "independent_audit_pr": I2_AUDIT_PR,
                "independent_audit_head": I2_AUDIT_HEAD,
                "changes_U_or_M": False,
            },
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "the finite-N source loop and concrete I2 bumps remain repository-autonomous candidate choices",
                "the Agent-4 I2 audit certifies only the local continuous three-row moment seam, not Navier-Stokes",
                "float64 velocity can round the Decimal I2 swirl increment below one ulp; use i2_decimal_correction_report for the local correction receipt",
                "I3/I4, terminal/global leading assembly, matched pressure/forcing and held-out full NS residual remain unassessed",
            ],
        }


def _main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--config-output", required=True)
    args = parser.parse_args()
    candidate = KokunoPA16CurrentCartesianI2HeatRepair()
    candidate.save_configuration(args.config_output)
    target = Path(args.output)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps(candidate.report(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    _main()
