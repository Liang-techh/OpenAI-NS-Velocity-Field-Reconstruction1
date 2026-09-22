"""Overflow-safe log-X Cartesian composition of the current Kokuno main pulse.

This bounded Agent-1 increment stacks exactly on PR #1100 head
``2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156``.  #1100 already composes the
current pulse-entry state with Kokuno's public main-pulse kernel and preserves
the real current M/M_eta history, but its Cartesian path forms the similarity
radius X explicitly.  For the frozen autonomous lambda=0.05 lineage, the
public endpoint xi=11 has finite physical radius while X itself exceeds the
float64 range.

Here the same candidate is evaluated in

    log X = 2 log(hypot(x,y)) - log(2 q),

so no x^2+y^2 or explicit X is formed in the active large-X pulse.  The source
profile remains

    y=log X-log X_p,  xi=lambda y,
    E=E_p exp[-(1/2+lambda)y],
    F=E exp[-(log 2+log X)/2],
    U=E A_principal R_0(xi),

and the normalized primitive carry is

    M/X = exp(log X_p-log X)(M_p/X_p) + Delta(M/X),

with the same public/source-coordinate integral used by #1100.  No new
physical parameter is introduced.  In particular A_principal remains the
repository-autonomous principal-term approximation from #1094; this module
does not recover source-exact Amp(eta), c1/c2 end compensators, pressure,
forcing, or any independent PDE validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_main_pulse import (
    KokunoPA16CurrentCartesianMainPulse,
)
from .kokuno_source_main_axial_pulse_kernel import (
    MAIN_XI_MAX,
    QUADRATURE_ORDER,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_RELEASE,
    SOURCE_RELEASE_DATE,
    SOURCE_REPOSITORY,
    main_kernel_R0,
    main_kernel_R0_prime,
)

SCHEMA = "kokuno-pa16-current-cartesian-main-pulse-logx-v1"
PARENT_EXACT_HEAD = "2a6c58ff6bf6d1ba6b875bca0e6abf9d1118e156"
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "overflow_safe_coordinate": "log X=2 log hypot(x,y)-log(2q)",
    "pulse_coordinate": "y=log X-log X_p, xi=lambda*y",
    "pulse_E": "E=E_p exp[-(1/2+lambda)y]",
    "pulse_F": "F=E exp[-(log 2+log X)/2]",
    "pulse_U": "U=E A_principal R_0(xi)",
    "primitive_carry": "M/X=exp(log X_p-log X)(M_p/X_p)+Delta(M/X)",
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact #1100 current Cartesian main-pulse prefix identity",
    "coordinate_change": "evaluate active pulse in log X; no new source/physical parameter",
    "primitive_quadrature": f"fixed {QUADRATURE_ORDER}-point Gauss-Legendre rule in xi, identical order to #1100",
    "source_endpoint": "materialize the public principal main-kernel through xi=11 inclusively",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "overflow_safe_logX_similarity_materialized": True,
    "full_source_xi_11_current_cartesian_materialized": True,
    "source_exact_amplitude_root_materialized": False,
    "source_pulse_end_MJ_corrections_materialized": False,
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
class KokunoPA16CurrentCartesianMainPulseLogX:
    """Current Cartesian leading candidate through the complete xi<=11 main kernel."""

    parent: KokunoPA16CurrentCartesianMainPulse = field(
        default_factory=KokunoPA16CurrentCartesianMainPulse,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianMainPulse):
            raise TypeError("parent must be KokunoPA16CurrentCartesianMainPulse")
        if not (0.0 < self.lambda_value < 1.0):
            raise ValueError("current outer lambda must lie in (0,1)")
        if not self.log_X_source_end > self.log_X_p:
            raise ValueError("invalid public main-pulse log-X interval")
        if not self.log_X_source_end > math.log(np.finfo(float).max):
            raise ValueError(
                "this schema is specifically the overflow-safe continuation beyond finite float64 X"
            )
        # The physical cylindrical radius at q=1 must itself remain representable;
        # otherwise log-X would not cure the Cartesian representation problem.
        if not self.log_radius_q1_source_end < math.log(np.finfo(float).max):
            raise ValueError("public main-pulse endpoint has non-representable physical radius")

    @property
    def current(self):
        return self.parent.current

    @property
    def outer_schedule(self):
        return self.parent.outer_schedule

    @property
    def geometry(self):
        return self.parent.geometry

    @property
    def kernel(self):
        return self.parent.kernel

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
    def X_p(self) -> float:
        return float(self.parent.X_p)

    @property
    def log_X_p(self) -> float:
        return float(self.parent.log_X_p)

    @property
    def lambda_value(self) -> float:
        return float(self.parent.lambda_value)

    @property
    def log_X_source_end(self) -> float:
        return self.log_X_p + MAIN_XI_MAX / self.lambda_value

    @property
    def log_radius_q1_source_end(self) -> float:
        return 0.5 * (_LOG2 + self.log_X_source_end)

    def log_X_from_xi(self, xi: float) -> float:
        value = float(xi)
        if not math.isfinite(value) or not (0.0 <= value <= MAIN_XI_MAX):
            raise ValueError(f"xi must lie in [0,{MAIN_XI_MAX:g}]")
        return self.log_X_p + value / self.lambda_value

    def similarity_coordinates_logX(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        xb, yb, zb, tb = np.broadcast_arrays(
            _finite(x, "x"), _finite(y, "y"), _finite(z, "z"), _finite(t, "t")
        )
        # Reuse the exact current q/eta solver while setting the dummy radial
        # coordinate to zero so the parent's finite-X arithmetic cannot overflow.
        base = self.parent.similarity_coordinates(
            np.zeros_like(xb), np.zeros_like(yb), zb, tb
        )
        q = np.asarray(base["q"], dtype=float)
        eta = np.asarray(base["eta"], dtype=float)
        if np.any(~np.isfinite(q)) or np.any(q <= 0.0):
            raise RuntimeError("current similarity geometry produced invalid q")
        radius = np.hypot(xb, yb)
        log_X = np.full(radius.shape, -np.inf, dtype=float)
        nonaxis = radius > 0.0
        if np.any(nonaxis):
            log_X[nonaxis] = (
                2.0 * np.log(radius[nonaxis])
                - _LOG2
                - np.log(q[nonaxis])
            )
        return {
            "q": q,
            "eta": eta,
            "radius": radius,
            "log_X": log_X,
            "x": xb,
            "y": yb,
            "z": zb,
            "t": tb,
        }

    def _broadcast_log_similarity(
        self, log_X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        log_arr, eta_arr = np.broadcast_arrays(
            _finite(log_X, "log_X"), _finite(eta, "eta")
        )
        tol = 256.0 * np.finfo(float).eps * max(1.0, abs(self.log_X_source_end))
        if np.any(log_arr > self.log_X_source_end + tol):
            raise ValueError("log_X lies beyond the public xi=11 main-kernel endpoint")
        log_arr = np.minimum(log_arr, self.log_X_source_end)
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return log_arr, eta_arr

    def _pulse_entry_state(
        self, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        return self.parent._pulse_entry_state(eta)

    def _normalized_pulse_integral(self, xi: np.ndarray) -> np.ndarray:
        values = np.asarray(xi, dtype=float)
        if np.any((values < 0.0) | (values > MAIN_XI_MAX + 2.0e-13)):
            raise ValueError("xi lies outside the public main-kernel interval")
        values = np.minimum(values, MAIN_XI_MAX)
        out = np.zeros_like(values)
        flat, out_flat = values.reshape(-1), out.reshape(-1)
        lam = self.lambda_value
        for i, upper in enumerate(flat):
            if upper <= 0.0:
                continue
            z = 0.5 * upper * (_GL_NODES + 1.0)
            exponent = ((0.5 - lam) * z - upper) / lam
            integrand = np.exp(exponent) * main_kernel_R0(z)
            out_flat[i] = 0.5 * upper * float(np.dot(_GL_WEIGHTS, integrand))
        return out

    def _active_profile_logX(
        self, log_X: np.ndarray, eta: np.ndarray
    ) -> dict[str, np.ndarray]:
        E_entry, E_entry_eta, m_entry, m_eta_entry = self._pulse_entry_state(eta)
        y = np.asarray(log_X, dtype=float) - self.log_X_p
        xi = self.lambda_value * y
        if np.any((xi < -2.0e-13) | (xi > MAIN_XI_MAX + 2.0e-13)):
            raise ValueError("active log-X profile lies outside 0<=xi<=11")
        xi = np.clip(xi, 0.0, MAIN_XI_MAX)
        decay = np.exp(-(0.5 + self.lambda_value) * y)
        E = E_entry * decay
        F = E * np.exp(-0.5 * (_LOG2 + log_X))
        R0 = main_kernel_R0(xi)
        U = E * self.kernel.A_principal * R0
        norm_int = self._normalized_pulse_integral(xi)
        amp_over_lambda = self.kernel.A_principal / self.lambda_value
        delta_m = E_entry * amp_over_lambda * norm_int
        delta_m_eta = E_entry_eta * amp_over_lambda * norm_int
        carry = np.exp(self.log_X_p - log_X)
        m_ratio = carry * m_entry + delta_m
        m_eta_ratio = carry * m_eta_entry + delta_m_eta
        axis = self.geometry.physical_profiles.axis_profiles.values(
            np.zeros_like(eta), eta
        )
        L = np.asarray(axis["L"], dtype=float)
        d = np.asarray(axis["d"], dtype=float)
        v0 = (
            2.0 * eta * U
            - 2.0 * self.D * eta * m_ratio
            - d * m_eta_ratio
        ) / L
        return {
            "log_X": np.asarray(log_X, dtype=float),
            "eta": np.asarray(eta, dtype=float),
            "y": y,
            "xi": xi,
            "R0": R0,
            "E_entry": E_entry,
            "E_entry_eta": E_entry_eta,
            "F": F,
            "U": U,
            "E": E,
            "M_over_X": m_ratio,
            "M_eta_over_X": m_eta_ratio,
            "v0": v0,
            "delta_M_over_X": delta_m,
            "delta_M_eta_over_X": delta_m_eta,
        }

    def similarity_profile_values_logX(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        log_arr, eta_arr = self._broadcast_log_similarity(log_X, eta)
        shape = log_arr.shape
        lf, ef = log_arr.reshape(-1), eta_arr.reshape(-1)
        F = np.empty_like(lf)
        U = np.empty_like(lf)
        E = np.empty_like(lf)
        m_ratio = np.empty_like(lf)
        m_eta_ratio = np.empty_like(lf)
        v0 = np.empty_like(lf)
        delta_m = np.zeros_like(lf)
        delta_m_eta = np.zeros_like(lf)
        xi = np.zeros_like(lf)
        region = np.empty(lf.shape, dtype=object)

        inherited = lf <= self.log_X_p
        if np.any(inherited):
            X = np.exp(lf[inherited])
            p = self.parent.similarity_profile_values(X, ef[inherited])
            F[inherited] = np.asarray(p["F_current_leading_with_main_pulse"], dtype=float)
            U[inherited] = np.asarray(p["U_current_leading_with_main_pulse"], dtype=float)
            E[inherited] = np.asarray(p["E_current_leading_with_main_pulse"], dtype=float)
            m_ratio[inherited] = np.asarray(
                p["M_over_X_current_leading_with_main_pulse"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_leading_with_main_pulse"], dtype=float
            )
            v0[inherited] = np.asarray(
                p["v0_current_leading_with_main_pulse"], dtype=float
            )
            region[inherited] = np.asarray(p["region"], dtype=object)

        active = ~inherited
        if np.any(active):
            p = self._active_profile_logX(lf[active], ef[active])
            F[active] = p["F"]
            U[active] = p["U"]
            E[active] = p["E"]
            m_ratio[active] = p["M_over_X"]
            m_eta_ratio[active] = p["M_eta_over_X"]
            v0[active] = p["v0"]
            delta_m[active] = p["delta_M_over_X"]
            delta_m_eta[active] = p["delta_M_eta_over_X"]
            xi[active] = p["xi"]
            region[active] = "current_cartesian_main_pulse_principal_logX_full"

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("log-X current main-pulse candidate produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("log-X current main-pulse candidate lost positive F/E")
        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "xi_current_main_pulse_logX": xi.reshape(shape),
            "F_current_leading_with_main_pulse_logX": F.reshape(shape),
            "U_current_leading_with_main_pulse_logX": U.reshape(shape),
            "E_current_leading_with_main_pulse_logX": E.reshape(shape),
            "M_over_X_current_leading_with_main_pulse_logX": m_ratio.reshape(shape),
            "M_eta_over_X_current_leading_with_main_pulse_logX": m_eta_ratio.reshape(shape),
            "v0_current_leading_with_main_pulse_logX": v0.reshape(shape),
            "delta_M_over_X_main_pulse_logX": delta_m.reshape(shape),
            "delta_M_eta_over_X_main_pulse_logX": delta_m_eta.reshape(shape),
        }

    def similarity_log_radial_derivatives(
        self, log_X: Any, eta: Any
    ) -> dict[str, np.ndarray]:
        """Return D_X=X*d/dX=d/d(log X) derivatives without forming X."""
        log_arr, eta_arr = self._broadcast_log_similarity(log_X, eta)
        shape = log_arr.shape
        lf, ef = log_arr.reshape(-1), eta_arr.reshape(-1)
        F_D = np.empty_like(lf)
        U_D = np.empty_like(lf)
        E_D = np.empty_like(lf)
        inherited = lf <= self.log_X_p
        if np.any(inherited):
            X = np.exp(lf[inherited])
            d = self.parent.similarity_radial_derivatives(X, ef[inherited])
            F_D[inherited] = X * np.asarray(
                d["F_current_leading_with_main_pulse_X"], dtype=float
            )
            U_D[inherited] = X * np.asarray(
                d["U_current_leading_with_main_pulse_X"], dtype=float
            )
            E_D[inherited] = X * np.asarray(
                d["E_current_leading_with_main_pulse_X"], dtype=float
            )
        active = ~inherited
        if np.any(active):
            p = self._active_profile_logX(lf[active], ef[active])
            xi = np.asarray(p["xi"], dtype=float)
            E = np.asarray(p["E"], dtype=float)
            F = np.asarray(p["F"], dtype=float)
            U = np.asarray(p["U"], dtype=float)
            a = 0.5 + self.lambda_value
            E_D[active] = -a * E
            F_D[active] = -(1.0 + self.lambda_value) * F
            U_D[active] = (
                -a * U
                + E
                * self.kernel.A_principal
                * self.lambda_value
                * main_kernel_R0_prime(xi)
            )
        return {
            "log_X": lf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_DlogX_current_leading_with_main_pulse": F_D.reshape(shape),
            "U_DlogX_current_leading_with_main_pulse": U_D.reshape(shape),
            "E_DlogX_current_leading_with_main_pulse": E_D.reshape(shape),
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        coords = self.similarity_coordinates_logX(x, y, z, t)
        xb = np.asarray(coords["x"], dtype=float)
        yb = np.asarray(coords["y"], dtype=float)
        q = np.asarray(coords["q"], dtype=float)
        eta = np.asarray(coords["eta"], dtype=float)
        radius = np.asarray(coords["radius"], dtype=float)
        log_X = np.asarray(coords["log_X"], dtype=float)

        F = np.empty_like(radius)
        U = np.empty_like(radius)
        v0 = np.empty_like(radius)
        axis = radius == 0.0
        if np.any(axis):
            p = self.parent.similarity_profile_values(
                np.zeros(np.count_nonzero(axis), dtype=float), eta[axis]
            )
            F[axis] = np.asarray(p["F_current_leading_with_main_pulse"], dtype=float)
            U[axis] = np.asarray(p["U_current_leading_with_main_pulse"], dtype=float)
            v0[axis] = np.asarray(p["v0_current_leading_with_main_pulse"], dtype=float)
        nonaxis = ~axis
        if np.any(nonaxis):
            p = self.similarity_profile_values_logX(log_X[nonaxis], eta[nonaxis])
            F[nonaxis] = np.asarray(
                p["F_current_leading_with_main_pulse_logX"], dtype=float
            )
            U[nonaxis] = np.asarray(
                p["U_current_leading_with_main_pulse_logX"], dtype=float
            )
            v0[nonaxis] = np.asarray(
                p["v0_current_leading_with_main_pulse_logX"], dtype=float
            )

        radial = v0 / (2.0 * q)
        swirl = q ** (-self.A - 0.5) * F
        axial = q ** (-self.A) * U
        out = np.stack(
            (radial * xb - swirl * yb, radial * yb + swirl * xb, axial), axis=-1
        )
        if np.any(~np.isfinite(out)):
            raise RuntimeError("overflow-safe current Cartesian velocity became non-finite")
        return out

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_finite_X_main_pulse": self.parent.configuration(),
            "bound_scope": {
                "source_commit": SOURCE_COMMIT,
                "source_role": "current_cartesian_main_pulse_principal_logX_through_xi_11",
                "amplitude_role": "repository_autonomous_principal_not_source_exact",
                "full_source_xi_11_current_cartesian_materialized": True,
                "source_pulse_end_MJ_corrections_materialized": False,
                "primitive_quadrature_order": QUADRATURE_ORDER,
                "mutable": False,
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianMainPulseLogX":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected log-X current Cartesian main-pulse schema")
        parent_payload = payload.get("parent_finite_X_main_pulse")
        bound = payload.get("bound_scope")
        if not isinstance(parent_payload, Mapping) or not isinstance(bound, Mapping):
            raise ValueError("missing log-X parent or source/truth binding")
        expected = {
            "source_commit": SOURCE_COMMIT,
            "source_role": "current_cartesian_main_pulse_principal_logX_through_xi_11",
            "amplitude_role": "repository_autonomous_principal_not_source_exact",
            "full_source_xi_11_current_cartesian_materialized": True,
            "source_pulse_end_MJ_corrections_materialized": False,
            "primitive_quadrature_order": QUADRATURE_ORDER,
            "mutable": False,
        }
        if dict(bound) != expected:
            raise ValueError("serialized log-X current main-pulse source/truth binding changed")
        parent = KokunoPA16CurrentCartesianMainPulse.from_configuration(parent_payload)
        return cls(parent=parent)

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA16CurrentCartesianMainPulseLogX":
        return cls.from_configuration(json.loads(Path(path).read_text(encoding="utf-8")))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": SOURCE_RELEASE,
                "release_date": SOURCE_RELEASE_DATE,
            },
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "domain": {
                "log_X_p": self.log_X_p,
                "lambda": self.lambda_value,
                "source_main_xi_max": MAIN_XI_MAX,
                "log_X_source_end": self.log_X_source_end,
            },
            "truth_boundary": self.truth_boundary,
        }
        return hashlib.sha256(_canonical_json(payload).encode()).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        truth = dict(self.parent.truth_boundary)
        truth.update(_TRUTH_UPDATES)
        return truth

    def receipt(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "semantic_sha256": self.semantic_sha256,
            "source_commit": SOURCE_COMMIT,
            "log_X_p": self.log_X_p,
            "lambda": self.lambda_value,
            "source_main_xi_max": MAIN_XI_MAX,
            "log_X_source_end": self.log_X_source_end,
            "float64_log_X_max": math.log(np.finfo(float).max),
            "log_radius_q1_source_end": self.log_radius_q1_source_end,
            "A_principal": self.kernel.A_principal,
            "truth_boundary": self.truth_boundary,
            "warning": (
                "Full public principal main-kernel Cartesian materialization through xi=11 "
                "via log-X numerics only; source-exact Amp(eta), c1/c2 M/J end closure, "
                "global pressure/forcing and independent PDE validation remain absent."
            ),
        }


__all__ = ["KokunoPA16CurrentCartesianMainPulseLogX"]
