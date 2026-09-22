"""Compose the current A1 Cartesian leading field through Kokuno's pulse-end bumps.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
navier-stokes/navier_stokes_workbench.tex
corrected 2026-09-09 reconstruction.

This bounded Agent-1 increment stacks exactly on A1 #1124.  The parent already
materializes the actual current-lineage pulse-entry quantities

    m_p = M(X_p)/(X_p E_p),
    j_p = J(X_p)/(X_p H_p E_p),

without assuming J=0, and #1116 materializes the public two-row end solve

    R_b = Amp R_0(lambda y) + c1 beta1(y) + c2 beta2(y),
    y = log(X/X_p),

with public row slopes s1=1/2-lambda and s2=1/2-2lambda.  Here those pieces are
finally composed into the current Cartesian leading candidate over the two end
bumps, extending the callable log-X field from the principal endpoint xi=11 to
the full public pulse endpoint xi=13.

Scientific boundary
-------------------
The executable candidate keeps the repository-autonomous A_principal from
#1094/#1107 and the repository-autonomous pointwise C-infinity bump realization
from #1116.  Neither is promoted to hidden/source-exact data.  This module does
not add pressure, forcing, a terminal/exterior tail, a complete held-out NS
residual, or PDE validation.

Current primitive continuity
----------------------------
The exact #1107 current M/M_eta state at xi=11 is preserved.  For xi>11 the
main R_0 kernel is already zero and only the two end bumps contribute.  Rather
than reconstructing M from a detached idealized row, the current primitive is
continued as

    M/X = exp(y_11-y) (M/X)_11
          + E_p exp(s1*y1-y) sum_j c_j C_1j(y),

where C_ij(y)=int exp(s_i(s-y1)) beta_j(s) ds is the source-prescribed row-scaled
cumulative bump moment.  The eta derivative is propagated analytically from
E_p, c1, c2 and their eta jets.  This gives the implementation-distinct identity

    d_logX(M/X) = U - M/X

inside each bump while maintaining an exact seam with #1107 at xi=11.

The public M/J two-row state is also exposed in the same scaled coordinates.
Its endpoint closure is algebraic evidence for this autonomous realization; it
is not a complete Navier-Stokes validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_current_pulse_entry_moments import KokunoCurrentPulseEntryMoments
from .kokuno_public_pulse_end_compensator import (
    MAIN_XI_END,
    PULSE_XI_END,
    QUADRATURE_ORDER,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_RELEASE,
    SOURCE_RELEASE_DATE,
    SOURCE_REPOSITORY,
)

SCHEMA = "kokuno-pa16-current-cartesian-pulse-end-compensated-v1"
PARENT_EXACT_HEAD = "4be97cbeb3b78b4d6a162447ae8fea73d0c2c477"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "pulse_coordinate": "y=log(X/X_p), xi=lambda*y",
    "pulse_E": "E=E_p exp[-(1/2+lambda)y]",
    "pulse_F": "F=E/sqrt(2X)",
    "pulse_U": "U=E R_b",
    "pulse_ratio": "R_b=Amp(eta)R_0(lambda*y)+c1(eta)beta1(y)+c2(eta)beta2(y)",
    "row_slopes": "s1=1/2-lambda, s2=1/2-2lambda",
    "bump_centers": "y1=13/lambda-3, y2=13/lambda-1",
    "bump_width": "beta1,beta2 have identical width 0.3",
    "M": "M(X)=int_0^X U dx",
    "J": "J(X)=int_0^X U H dx, H=sqrt(2X)E",
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1124 current pulse-entry M/J identity",
    "amplitude": "repository-autonomous A_principal inherited unchanged from #1094/#1107",
    "bump_shape": "repository-autonomous pointwise C-infinity bump inherited unchanged from #1116",
    "coordinate": "overflow-safe log-X/log-F continuation inherited from #1107",
    "cumulative_bump_quadrature": (
        f"fixed {QUADRATURE_ORDER}-point Gauss-Legendre on the exact compact bump support; "
        "full-support values reuse #1116's precomputed row-scaled matrix"
    ),
    "current_M_policy": "preserve exact #1107 M/M_eta at xi=11, then add only end-bump increments",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "current_cartesian_end_compensation_composed": True,
    "current_pulse_endpoint_xi13_materialized": True,
    "public_MJ_end_compensator_algebra_consumed": True,
    "public_end_bump_formula_candidate_materialized": True,
    "current_endpoint_MJ_closure_numerically_assessed": True,
    "source_exact_amplitude_root_materialized": False,
    "source_exact_bump_shape_recovered": False,
    "source_hidden_parameters_recovered": False,
    "source_terminal_tail_schedule_bound_into_current_velocity": False,
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


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    return out


_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(QUADRATURE_ORDER)
_GL_NODES = np.asarray(_GL_NODES, dtype=float)
_GL_WEIGHTS = np.asarray(_GL_WEIGHTS, dtype=float)
_GL_NODES.setflags(write=False)
_GL_WEIGHTS.setflags(write=False)


@dataclass(frozen=True)
class KokunoPA16CurrentCartesianPulseEndCompensated:
    """Current leading Cartesian candidate through the full public pulse endpoint xi=13."""

    parent: KokunoCurrentPulseEntryMoments = field(
        default_factory=KokunoCurrentPulseEntryMoments,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoCurrentPulseEntryMoments):
            raise TypeError("parent must be KokunoCurrentPulseEntryMoments")
        if not self.parent.truth_boundary["current_lineage_J_entry_materialized"]:
            raise ValueError("pulse-end composition requires current-lineage J entry")
        if self.parent.truth_boundary["current_cartesian_end_compensation_composed"]:
            raise ValueError("parent identity unexpectedly already composes pulse-end correction")
        if self.lambda_value != self.compensator.lambda_value:
            raise ValueError("pulse-entry and end-compensator lambda mismatch")
        if not (MAIN_XI_END < PULSE_XI_END):
            raise ValueError("invalid public pulse endpoint ordering")
        if not self.log_radius_q1_pulse_end < math.log(np.finfo(float).max):
            raise ValueError("public pulse endpoint has non-representable physical radius")

    @property
    def leading(self):
        return self.parent.leading

    @property
    def compensator(self):
        return self.parent.compensator

    @property
    def geometry(self):
        return self.leading.geometry

    @property
    def A(self) -> float:
        return float(self.leading.A)

    @property
    def D(self) -> float:
        return float(self.leading.D)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.leading.eta_interval)

    @property
    def lambda_value(self) -> float:
        return float(self.leading.lambda_value)

    @property
    def log_X_p(self) -> float:
        return float(self.leading.log_X_p)

    @property
    def log_X_main_end(self) -> float:
        return float(self.leading.log_X_source_end)

    @property
    def log_X_pulse_end(self) -> float:
        return self.log_X_p + PULSE_XI_END / self.lambda_value

    @property
    def log_radius_q1_pulse_end(self) -> float:
        return 0.5 * (_LOG2 + self.log_X_pulse_end)

    def log_X_from_xi(self, xi: float) -> float:
        value = float(xi)
        if not math.isfinite(value) or not (0.0 <= value <= PULSE_XI_END):
            raise ValueError(f"xi must lie in [0,{PULSE_XI_END:g}]")
        return self.log_X_p + value / self.lambda_value

    def similarity_coordinates_logX(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.leading.similarity_coordinates_logX(x, y, z, t)

    def _broadcast_log_similarity(
        self, log_X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_pulse_end))
        if np.any(log_arr > self.log_X_pulse_end + tol):
            raise ValueError("log_X lies beyond the public xi=13 pulse endpoint")
        log_arr = np.minimum(log_arr, self.log_X_pulse_end)
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return log_arr, eta_arr

    def _partial_scaled_bump_integral_scalar(self, y: float, row: int, column: int) -> float:
        if row not in (0, 1) or column not in (0, 1):
            raise ValueError("row and column must be 0 or 1")
        slope = (self.compensator.s1, self.compensator.s2)[row]
        center = (self.compensator.y1, self.compensator.y2)[column]
        beta = (self.compensator.beta1, self.compensator.beta2)[column]
        half = 0.5 * self.compensator.width
        lo, hi = center - half, center + half
        if y <= lo:
            return 0.0
        if y >= hi:
            return float(self.compensator.matrix_scaled[row, column])
        half_interval = 0.5 * (y - lo)
        mid = 0.5 * (y + lo)
        points = half_interval * _GL_NODES + mid
        integrand = np.exp(slope * (points - self.compensator.y1)) * beta(points)
        return float(half_interval * np.dot(_GL_WEIGHTS, integrand))

    def _cumulative_scaled_bump_integrals(self, y: np.ndarray) -> np.ndarray:
        values = np.asarray(y, dtype=float)
        out = np.empty(values.shape + (2, 2), dtype=float)
        flat = values.reshape(-1)
        flat_out = out.reshape((-1, 2, 2))
        for i, value in enumerate(flat):
            for row in range(2):
                for column in range(2):
                    flat_out[i, row, column] = self._partial_scaled_bump_integral_scalar(
                        float(value), row, column
                    )
        return out

    def _coefficient_state(self, eta: np.ndarray) -> dict[str, np.ndarray]:
        inputs = self.parent.compensator_inputs(eta)
        sol, c1_eta, c2_eta = self.compensator.solve_eta_jet(
            inputs["amplitude"],
            inputs["m_entry_ratio"],
            inputs["j_entry_ratio"],
            inputs["amplitude_eta"],
            inputs["m_entry_ratio_eta"],
            inputs["j_entry_ratio_eta"],
        )
        return {
            "amplitude": np.asarray(inputs["amplitude"], dtype=float),
            "m_entry": np.asarray(inputs["m_entry_ratio"], dtype=float),
            "j_entry": np.asarray(inputs["j_entry_ratio"], dtype=float),
            "m_entry_eta": np.asarray(inputs["m_entry_ratio_eta"], dtype=float),
            "j_entry_eta": np.asarray(inputs["j_entry_ratio_eta"], dtype=float),
            "c1": np.asarray(sol.c1, dtype=float),
            "c2": np.asarray(sol.c2, dtype=float),
            "c1_eta": np.asarray(c1_eta, dtype=float),
            "c2_eta": np.asarray(c2_eta, dtype=float),
            "row1_residual": np.asarray(sol.residual_row1, dtype=float),
            "row2_residual": np.asarray(sol.residual_row2, dtype=float),
        }

    def _end_profile_logX(
        self, log_X: np.ndarray, eta: np.ndarray
    ) -> dict[str, np.ndarray]:
        y = np.asarray(log_X, dtype=float) - self.log_X_p
        xi = self.lambda_value * y
        if np.any((xi <= MAIN_XI_END - 2.0e-13) | (xi > PULSE_XI_END + 2.0e-13)):
            raise ValueError("pulse-end profile requires 11<xi<=13")
        xi = np.minimum(xi, PULSE_XI_END)

        E_entry, E_entry_eta, _, _ = self.leading._pulse_entry_state(eta)
        E_entry = np.asarray(E_entry, dtype=float)
        E_entry_eta = np.asarray(E_entry_eta, dtype=float)
        decay = np.exp(-(0.5 + self.lambda_value) * y)
        E = E_entry * decay
        E_eta = E_entry_eta * decay
        log_E = np.log(E_entry) - (0.5 + self.lambda_value) * y
        log_F = log_E - 0.5 * (_LOG2 + log_X)
        F = np.exp(log_F)

        coeff = self._coefficient_state(eta)
        c1, c2 = coeff["c1"], coeff["c2"]
        c1_eta, c2_eta = coeff["c1_eta"], coeff["c2_eta"]
        beta1 = self.compensator.beta1(y)
        beta2 = self.compensator.beta2(y)
        ratio = c1 * beta1 + c2 * beta2
        ratio_eta = c1_eta * beta1 + c2_eta * beta2
        U = E * ratio
        U_eta = E_eta * ratio + E * ratio_eta

        integrals = self._cumulative_scaled_bump_integrals(y)
        row1_corr = c1 * integrals[..., 0, 0] + c2 * integrals[..., 0, 1]
        row1_corr_eta = (
            c1_eta * integrals[..., 0, 0] + c2_eta * integrals[..., 0, 1]
        )
        row2_corr = c1 * integrals[..., 1, 0] + c2 * integrals[..., 1, 1]

        main_end_log = np.full_like(eta, self.log_X_main_end, dtype=float)
        main = self.leading.similarity_profile_values_logX(main_end_log, eta)
        M_main = np.asarray(
            main["M_over_X_current_leading_with_main_pulse_logX"], dtype=float
        )
        M_eta_main = np.asarray(
            main["M_eta_over_X_current_leading_with_main_pulse_logX"], dtype=float
        )
        y_main = MAIN_XI_END / self.lambda_value
        carry = np.exp(y_main - y)
        row1_factor = np.exp(self.compensator.s1 * self.compensator.y1 - y)
        M_ratio = carry * M_main + E_entry * row1_factor * row1_corr
        M_eta_ratio = (
            carry * M_eta_main
            + row1_factor
            * (E_entry_eta * row1_corr + E_entry * row1_corr_eta)
        )

        m_scaled = (
            coeff["m_entry"] * math.exp(-self.compensator.s1 * self.compensator.y1)
            + coeff["amplitude"] * self.compensator.main_moments_scaled[0]
            + row1_corr
        )
        j_scaled = (
            coeff["j_entry"] * math.exp(-self.compensator.s2 * self.compensator.y1)
            + coeff["amplitude"] * self.compensator.main_moments_scaled[1]
            + row2_corr
        )

        axis = self.geometry.physical_profiles.axis_profiles.values(
            np.zeros_like(eta), eta
        )
        L = np.asarray(axis["L"], dtype=float)
        d = np.asarray(axis["d"], dtype=float)
        v0 = (
            2.0 * eta * U
            - 2.0 * self.D * eta * M_ratio
            - d * M_eta_ratio
        ) / L

        arrays = (
            E,
            F,
            U,
            U_eta,
            M_ratio,
            M_eta_ratio,
            v0,
            m_scaled,
            j_scaled,
        )
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("pulse-end compensated profile produced non-finite values")
        if np.any(E <= 0.0) or np.any(F < 0.0):
            raise RuntimeError("pulse-end compensated profile lost positive E/nonnegative F")

        return {
            "log_X": np.asarray(log_X, dtype=float),
            "eta": np.asarray(eta, dtype=float),
            "y": y,
            "xi": xi,
            "E": E,
            "E_eta": E_eta,
            "F": F,
            "log_F": log_F,
            "U": U,
            "U_eta": U_eta,
            "M_over_X": M_ratio,
            "M_eta_over_X": M_eta_ratio,
            "v0": v0,
            "c1": c1,
            "c2": c2,
            "c1_eta": c1_eta,
            "c2_eta": c2_eta,
            "public_M_row_scaled": m_scaled,
            "public_J_row_scaled": j_scaled,
            "public_row1_solver_residual": coeff["row1_residual"],
            "public_row2_solver_residual": coeff["row2_residual"],
        }

    def similarity_profile_values_logX(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_log_similarity(log_X, eta)
        shape = log_arr.shape
        lf, ef = log_arr.reshape(-1), eta_arr.reshape(-1)

        F = np.empty_like(lf)
        log_F = np.empty_like(lf)
        U = np.empty_like(lf)
        E = np.empty_like(lf)
        M_ratio = np.empty_like(lf)
        M_eta_ratio = np.empty_like(lf)
        v0 = np.empty_like(lf)
        xi = np.zeros_like(lf)
        public_m = np.zeros_like(lf)
        public_j = np.zeros_like(lf)
        region = np.empty(lf.shape, dtype=object)

        inherited = lf <= self.log_X_main_end
        if np.any(inherited):
            p = self.leading.similarity_profile_values_logX(lf[inherited], ef[inherited])
            F[inherited] = np.asarray(
                p["F_current_leading_with_main_pulse_logX"], dtype=float
            )
            log_F[inherited] = np.asarray(
                p["log_F_current_leading_with_main_pulse_logX"], dtype=float
            )
            U[inherited] = np.asarray(
                p["U_current_leading_with_main_pulse_logX"], dtype=float
            )
            E[inherited] = np.asarray(
                p["E_current_leading_with_main_pulse_logX"], dtype=float
            )
            M_ratio[inherited] = np.asarray(
                p["M_over_X_current_leading_with_main_pulse_logX"], dtype=float
            )
            M_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_leading_with_main_pulse_logX"], dtype=float
            )
            v0[inherited] = np.asarray(
                p["v0_current_leading_with_main_pulse_logX"], dtype=float
            )
            xi[inherited] = np.asarray(
                p["xi_current_main_pulse_logX"], dtype=float
            )
            region[inherited] = np.asarray(p["region"], dtype=object)

        active = ~inherited
        if np.any(active):
            p = self._end_profile_logX(lf[active], ef[active])
            F[active] = p["F"]
            log_F[active] = p["log_F"]
            U[active] = p["U"]
            E[active] = p["E"]
            M_ratio[active] = p["M_over_X"]
            M_eta_ratio[active] = p["M_eta_over_X"]
            v0[active] = p["v0"]
            xi[active] = p["xi"]
            public_m[active] = p["public_M_row_scaled"]
            public_j[active] = p["public_J_row_scaled"]
            region[active] = "current_cartesian_pulse_end_compensated_logX"

        arrays = (F, log_F, U, E, M_ratio, M_eta_ratio, v0)
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise RuntimeError("current pulse-end Cartesian candidate produced non-finite values")
        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "xi_current_pulse_end": xi.reshape(shape),
            "F_current_leading_with_pulse_end": F.reshape(shape),
            "log_F_current_leading_with_pulse_end": log_F.reshape(shape),
            "U_current_leading_with_pulse_end": U.reshape(shape),
            "E_current_leading_with_pulse_end": E.reshape(shape),
            "M_over_X_current_leading_with_pulse_end": M_ratio.reshape(shape),
            "M_eta_over_X_current_leading_with_pulse_end": M_eta_ratio.reshape(shape),
            "v0_current_leading_with_pulse_end": v0.reshape(shape),
            "public_M_row_scaled_current_pulse_end": public_m.reshape(shape),
            "public_J_row_scaled_current_pulse_end": public_j.reshape(shape),
        }

    def end_eta_jet_logX(self, log_X: Any, eta: Any) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_log_similarity(log_X, eta)
        if np.any(log_arr <= self.log_X_main_end):
            raise ValueError("end eta-jet is exposed only on 11<xi<=13")
        p = self._end_profile_logX(log_arr, eta_arr)
        return {
            "log_X": log_arr,
            "eta": eta_arr,
            "U_eta_current_leading_with_pulse_end": p["U_eta"],
            "M_eta_over_X_current_leading_with_pulse_end": p["M_eta_over_X"],
            "c1_eta": p["c1_eta"],
            "c2_eta": p["c2_eta"],
        }

    def similarity_log_radial_derivatives(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_log_similarity(log_X, eta)
        shape = log_arr.shape
        lf, ef = log_arr.reshape(-1), eta_arr.reshape(-1)
        F_D = np.empty_like(lf)
        U_D = np.empty_like(lf)
        E_D = np.empty_like(lf)
        log_F_D = np.empty_like(lf)

        inherited = lf <= self.log_X_main_end
        if np.any(inherited):
            d = self.leading.similarity_log_radial_derivatives(lf[inherited], ef[inherited])
            F_D[inherited] = np.asarray(
                d["F_DlogX_current_leading_with_main_pulse"], dtype=float
            )
            U_D[inherited] = np.asarray(
                d["U_DlogX_current_leading_with_main_pulse"], dtype=float
            )
            E_D[inherited] = np.asarray(
                d["E_DlogX_current_leading_with_main_pulse"], dtype=float
            )
            log_F_D[inherited] = np.asarray(
                d["log_F_DlogX_current_leading_with_main_pulse"], dtype=float
            )

        active = ~inherited
        if np.any(active):
            p = self._end_profile_logX(lf[active], ef[active])
            y = np.asarray(p["y"], dtype=float)
            E = np.asarray(p["E"], dtype=float)
            F = np.asarray(p["F"], dtype=float)
            U = np.asarray(p["U"], dtype=float)
            c1 = np.asarray(p["c1"], dtype=float)
            c2 = np.asarray(p["c2"], dtype=float)
            ratio_y = (
                c1 * self.compensator.beta1_prime(y)
                + c2 * self.compensator.beta2_prime(y)
            )
            a = 0.5 + self.lambda_value
            E_D[active] = -a * E
            F_D[active] = -(1.0 + self.lambda_value) * F
            log_F_D[active] = -(1.0 + self.lambda_value)
            U_D[active] = -a * U + E * ratio_y

        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_DlogX_current_leading_with_pulse_end": F_D.reshape(shape),
            "log_F_DlogX_current_leading_with_pulse_end": log_F_D.reshape(shape),
            "U_DlogX_current_leading_with_pulse_end": U_D.reshape(shape),
            "E_DlogX_current_leading_with_pulse_end": E_D.reshape(shape),
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        coords = self.similarity_coordinates_logX(x, y, z, t)
        xb = np.asarray(coords["x"], dtype=float)
        yb = np.asarray(coords["y"], dtype=float)
        zb = np.asarray(coords["z"], dtype=float)
        tb = np.asarray(coords["t"], dtype=float)
        q = np.asarray(coords["q"], dtype=float)
        eta = np.asarray(coords["eta"], dtype=float)
        radius = np.asarray(coords["radius"], dtype=float)
        log_X = np.asarray(coords["log_X"], dtype=float)

        out = np.empty(radius.shape + (3,), dtype=float)
        axis = radius == 0.0
        if np.any(axis):
            out[axis] = self.leading.velocity(
                xb[axis], yb[axis], zb[axis], tb[axis]
            )

        nonaxis = ~axis
        if np.any(nonaxis):
            p = self.similarity_profile_values_logX(log_X[nonaxis], eta[nonaxis])
            log_F = np.asarray(p["log_F_current_leading_with_pulse_end"], dtype=float)
            U = np.asarray(p["U_current_leading_with_pulse_end"], dtype=float)
            v0 = np.asarray(p["v0_current_leading_with_pulse_end"], dtype=float)
            unit_x = xb[nonaxis] / radius[nonaxis]
            unit_y = yb[nonaxis] / radius[nonaxis]
            radial_radius = v0 * radius[nonaxis] / (2.0 * q[nonaxis])
            log_swirl_radius = (
                (-self.A - 0.5) * np.log(q[nonaxis])
                + log_F
                + np.log(radius[nonaxis])
            )
            swirl_radius = np.exp(log_swirl_radius)
            u1 = radial_radius * unit_x - swirl_radius * unit_y
            u2 = radial_radius * unit_y + swirl_radius * unit_x
            axial = q[nonaxis] ** (-self.A) * U
            out[nonaxis] = np.stack((u1, u2, axial), axis=-1)

        if np.any(~np.isfinite(out)):
            raise RuntimeError("pulse-end compensated Cartesian velocity became non-finite")
        return out

    @property
    def truth_boundary(self) -> dict[str, bool]:
        truth = dict(self.parent.truth_boundary)
        truth.update(_TRUTH_UPDATES)
        return truth

    @property
    def source_formulas(self) -> dict[str, str]:
        return dict(_SOURCE_FORMULAS)

    @property
    def numerical_realization(self) -> dict[str, str]:
        return dict(_NUMERICAL_REALIZATION)

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_current_pulse_entry_moments": self.parent.configuration(),
            "bound_scope": {
                "parent_exact_head": PARENT_EXACT_HEAD,
                "source": {
                    "repository": SOURCE_REPOSITORY,
                    "commit": SOURCE_COMMIT,
                    "blob": SOURCE_BLOB,
                    "path": SOURCE_PATH,
                    "release": SOURCE_RELEASE,
                    "release_date": SOURCE_RELEASE_DATE,
                },
                "pulse_xi_interval": [MAIN_XI_END, PULSE_XI_END],
                "amplitude_role": "repository_autonomous_principal_not_source_exact",
                "bump_shape_role": "repository_autonomous_not_source_exact",
                "current_M_policy": "preserve_xi11_parent_then_integrate_end_bumps",
                "current_J_policy": "consume_actual_current_entry_J_in_public_two_row_algebra",
                "mutable": False,
            },
            "truth_boundary": self.truth_boundary,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianPulseEndCompensated":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current Cartesian pulse-end schema")
        parent_payload = payload.get("parent_current_pulse_entry_moments")
        if not isinstance(parent_payload, Mapping):
            raise ValueError("missing current pulse-entry parent configuration")
        candidate = cls(parent=KokunoCurrentPulseEntryMoments.from_configuration(parent_payload))
        if _canonical_json(dict(payload)) != _canonical_json(candidate.configuration()):
            raise ValueError("serialized pulse-end source/truth binding changed")
        return candidate

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianPulseEndCompensated":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "source_formulas": self.source_formulas,
            "numerical_realization": self.numerical_realization,
            "domain": {
                "log_X_p": self.log_X_p,
                "lambda": self.lambda_value,
                "main_xi_end": MAIN_XI_END,
                "pulse_xi_end": PULSE_XI_END,
                "log_X_pulse_end": self.log_X_pulse_end,
            },
            "truth_boundary": self.truth_boundary,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    def end_closure_report(
        self, eta: Any = (-0.83, -0.61, -0.37, 0.19, 0.73)
    ) -> dict[str, Any]:
        eta_arr = _finite(eta, "eta")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        end_log = np.full(eta_arr.shape, self.log_X_pulse_end, dtype=float)
        end = self._end_profile_logX(end_log, eta_arr)
        main_log = np.full(eta_arr.shape, self.log_X_main_end, dtype=float)
        main = self.leading.similarity_profile_values_logX(main_log, eta_arr)
        main_M = np.asarray(
            main["M_over_X_current_leading_with_main_pulse_logX"], dtype=float
        )
        end_M = np.asarray(end["M_over_X"], dtype=float)
        reduction = np.abs(end_M) / np.maximum(np.abs(main_M), 1.0e-300)
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "semantic_sha256": self.semantic_sha256,
            "eta": eta_arr.tolist(),
            "lambda_value": self.lambda_value,
            "main_xi_end": MAIN_XI_END,
            "pulse_xi_end": PULSE_XI_END,
            "log_X_pulse_end": self.log_X_pulse_end,
            "max_abs_public_M_row_scaled_endpoint": float(
                np.max(np.abs(end["public_M_row_scaled"]))
            ),
            "max_abs_public_J_row_scaled_endpoint": float(
                np.max(np.abs(end["public_J_row_scaled"]))
            ),
            "max_abs_solver_row1_residual": float(
                np.max(np.abs(end["public_row1_solver_residual"]))
            ),
            "max_abs_solver_row2_residual": float(
                np.max(np.abs(end["public_row2_solver_residual"]))
            ),
            "max_abs_actual_current_M_over_X_endpoint": float(np.max(np.abs(end_M))),
            "max_actual_M_endpoint_to_xi11_ratio": float(np.max(reduction)),
            "max_abs_actual_current_M_eta_over_X_endpoint": float(
                np.max(np.abs(end["M_eta_over_X"]))
            ),
            "truth_boundary": self.truth_boundary,
            "limitations": [
                "A_principal remains repository-autonomous rather than source-exact Amp(eta)",
                "the pointwise beta1/beta2 shape remains the frozen repository-autonomous realization from #1116",
                "public two-row endpoint closure is algebraic evidence for this realization, not complete NS validation",
                "the current physical M path is seam-preserving from #1107; J is tracked in the public scaled end-row algebra rather than a new global physical-J state object",
                "terminal/exterior heat, global leading velocity, matched pressure, restricted forcing, and held-out complete NS residual remain absent",
            ],
        }


__all__ = ["KokunoPA16CurrentCartesianPulseEndCompensated"]
