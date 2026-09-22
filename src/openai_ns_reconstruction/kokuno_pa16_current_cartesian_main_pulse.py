"""Compose the current Kokuno Cartesian leading field into the public main pulse.

This is one bounded Kokuno Agent-1 increment stacked exactly on A1 PR #1094
head ``dc33914d3dd60c1dc4c8dfaab0760ef071385d54``.  It consumes two already
separated contracts:

* #1088: the exact current Cartesian leading candidate through pulse entry X_p,
  including the real current M/M_eta incompressibility history;
* #1094: Kokuno's public main axial-pulse F/E/U kernel with the explicitly
  repository-autonomous principal-amplitude approximation.

For X>X_p, with y=log(X/X_p) and xi=lambda*y, #1094 supplies

    E = E_p exp[-(1/2+lambda)y],
    F = E/sqrt(2X),
    U = E A_principal R_0(xi).

The pulse-entry leading shape is f(eta)=(1+eta^2)^(-1), hence

    E_{p,eta} = -2 eta/(1+eta^2) E_p.

The current primitive history is not reset at X_p.  Writing M_p for the exact
#1088 primitive at the seam, this module adds the pulse contribution and keeps
it in normalized form to avoid large-X overflow,

    Delta(M/X)
      = E_p A_principal/lambda
        int_0^xi exp(((1/2-lambda)z-xi)/lambda) R_0(z) dz,

and the identical formula with E_{p,eta} for Delta(M_eta/X).  Then

    M/X = (X_p/X)(M_p/X_p) + Delta(M/X),

with the analogous M_eta/X, and v0 is rebuilt from the public
incompressibility formula.  The integral is evaluated by a fixed Gauss-
Legendre rule; this is deterministic numerical bookkeeping, not a residual
fit or an additional physical parameter.

A numerical boundary matters for this exact current lineage.  Its autonomous
X_p is already so large that float64 cannot represent X_p*exp(11/lambda), even
though the standalone #1094 low-dimensional source kernel is defined through
xi=11.  This Cartesian composition therefore exposes the largest safely
materializable finite-X prefix and fails closed beyond it.  It does NOT claim
to have materialized the missing final xi-prefix, the source-exact Amp(eta),
the c1/c2 M/J end compensators, a terminal/global velocity, pressure, forcing,
or independent PDE validation.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_post_i4_pulse_entry_bridge import (
    KokunoPA16CurrentCartesianPostI4PulseEntryBridge,
)
from .kokuno_source_main_axial_pulse_kernel import (
    MAIN_XI_MAX,
    QUADRATURE_ORDER,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_RELEASE,
    SOURCE_RELEASE_DATE,
    SOURCE_REPOSITORY,
    KokunoSourceMainAxialPulseKernel,
    main_kernel_R0,
)


SCHEMA = "kokuno-pa16-current-cartesian-main-pulse-v1"
PARENT_EXACT_HEAD = "dc33914d3dd60c1dc4c8dfaab0760ef071385d54"
_FLOAT_X_SAFETY_DIVISOR = 4.0

_SOURCE_FORMULAS = {
    "pulse_coordinate": "y=log(X/X_p), xi=lambda*y",
    "pulse_E": "E=E_p exp[-(1/2+lambda)y]",
    "pulse_F": "F=E/sqrt(2X)",
    "pulse_U": "U=E R_b; current bounded realization uses A_principal R_0(xi)",
    "entry_eta_jet": "E_p proportional to (1+eta^2)^(-1), so E_p,eta=-2 eta E_p/(1+eta^2)",
    "primitive": "M(X,eta)=int_0^X U(s,eta) ds",
    "normalized_pulse_primitive": (
        "Delta(M/X)=E_p A/lambda int_0^xi "
        "exp(((1/2-lambda)z-xi)/lambda) R_0(z) dz"
    ),
    "normalized_pulse_primitive_eta": (
        "replace E_p by E_p,eta in the same integral because the autonomous "
        "principal amplitude is eta-independent"
    ),
    "radial_profile": "v0=(2 eta U-2 D eta M/X-d M_eta/X)/L",
    "cartesian_velocity": (
        "u1=(v0/(2q))x-q^(-A-1/2)Fy; "
        "u2=(v0/(2q))y+q^(-A-1/2)Fx; u3=q^(-A)U"
    ),
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact #1088 current Cartesian pulse-entry candidate through #1094 ancestry",
    "pulse_kernel": "consume exact #1094 public R_0/F/E/U kernel contract",
    "amplitude": (
        "repository-autonomous #1094 A_principal; explicitly not the source-exact Amp(eta) root"
    ),
    "primitive_quadrature": (
        f"fixed {QUADRATURE_ORDER}-point Gauss-Legendre integration in xi of the "
        "stable normalized pulse primitive"
    ),
    "float_X_guard": (
        "stop before float64 X overflow and before 2X overflow in F=E/sqrt(2X); "
        f"safety divisor={_FLOAT_X_SAFETY_DIVISOR:g}"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "exact_current_pulse_entry_parent_consumed": True,
    "source_main_pulse_kernel_formula_consumed": True,
    "repository_autonomous_principal_amplitude_consumed": True,
    "current_cartesian_main_pulse_prefix_materialized": True,
    "current_pulse_M_and_M_eta_history_integrated": True,
    "current_v0_rebuilt_from_full_primitive_history": True,
    "velocity_interface_vectorized": True,
    "configuration_serializable": True,
    "source_exact_amplitude_root_materialized": False,
    "source_hidden_parameters_recovered": False,
    "source_pulse_end_MJ_corrections_materialized": False,
    "full_source_xi_11_current_cartesian_materialized": False,
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
class KokunoPA16CurrentCartesianMainPulse:
    """Current Cartesian leading candidate through the finite-X main-pulse prefix."""

    parent: KokunoPA16CurrentCartesianPostI4PulseEntryBridge = field(
        default_factory=KokunoPA16CurrentCartesianPostI4PulseEntryBridge,
        repr=False,
        compare=False,
    )
    kernel: KokunoSourceMainAxialPulseKernel = field(
        default_factory=KokunoSourceMainAxialPulseKernel,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA16CurrentCartesianPostI4PulseEntryBridge):
            raise TypeError(
                "parent must be KokunoPA16CurrentCartesianPostI4PulseEntryBridge"
            )
        if not isinstance(self.kernel, KokunoSourceMainAxialPulseKernel):
            raise TypeError("kernel must be KokunoSourceMainAxialPulseKernel")
        if not (0.0 < self.lambda_value < 1.0):
            raise ValueError("current outer lambda must lie in (0,1)")
        if not (self.log_X_p < self.log_X_materializable_end):
            raise ValueError("no finite-X main-pulse prefix is materializable")
        if not (0.0 < self.xi_materializable_max <= MAIN_XI_MAX):
            raise ValueError("invalid materializable main-pulse xi bound")
        # This exact lineage is expected to hit the float64-X ceiling before
        # xi=11.  Keep that limitation explicit rather than silently changing
        # the public low-dimensional kernel's source domain.
        if not (self.xi_materializable_max < MAIN_XI_MAX):
            raise ValueError(
                "this frozen current-lineage schema is specifically the finite-X prefix; "
                "a full-xi lineage requires a new schema"
            )

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
        # #1088 exposes the exact current RF40 power-law object as ``current``;
        # bind the pulse coordinate to that already-frozen autonomous source choice.
        return float(self.current.lambda_outer)

    @property
    def log_X_materializable_end(self) -> float:
        source_end = self.log_X_p + MAIN_XI_MAX / self.lambda_value
        float_end = math.log(np.finfo(float).max / _FLOAT_X_SAFETY_DIVISOR)
        return min(source_end, float_end)

    @property
    def X_materializable_end(self) -> float:
        return math.exp(self.log_X_materializable_end)

    @property
    def xi_materializable_max(self) -> float:
        return self.lambda_value * (self.log_X_materializable_end - self.log_X_p)

    def X_from_xi(self, xi: float) -> float:
        value = float(xi)
        if not math.isfinite(value) or not (0.0 <= value <= self.xi_materializable_max):
            raise ValueError(
                f"xi must lie in [0,{self.xi_materializable_max:.17g}] for the current finite-X lineage"
            )
        return math.exp(self.log_X_p + value / self.lambda_value)

    def similarity_coordinates(
        self, x: Any, y: Any, z: Any, t: Any
    ) -> dict[str, np.ndarray]:
        return self.parent.similarity_coordinates(x, y, z, t)

    def _broadcast_similarity(
        self, X: Any, eta: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_materializable_end)):
            raise ValueError(
                "X must lie in the current Cartesian main-pulse materialization domain "
                f"[0,{self.X_materializable_end:.17e}]"
            )
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return X_arr, eta_arr

    def _pulse_entry_state(
        self, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        values = np.asarray(eta, dtype=float)
        flat = values.reshape(-1)
        unique, inverse = np.unique(flat, return_inverse=True)
        X = np.full(unique.shape, self.X_p, dtype=float)
        p = self.parent.similarity_profile_values(X, unique)
        E = np.asarray(p["E_current_leading_to_pulse_entry"], dtype=float)
        m_ratio = np.asarray(
            p["M_over_X_current_leading_to_pulse_entry"], dtype=float
        )
        m_eta_ratio = np.asarray(
            p["M_eta_over_X_current_leading_to_pulse_entry"], dtype=float
        )
        E_eta = -2.0 * unique * E / (1.0 + unique * unique)
        return tuple(
            arr[inverse].reshape(values.shape)
            for arr in (E, E_eta, m_ratio, m_eta_ratio)
        )

    def _normalized_pulse_integral(self, xi: np.ndarray) -> np.ndarray:
        """Return stable normalized integral multiplying E_entry*A/lambda."""
        values = np.asarray(xi, dtype=float)
        if np.any((values < 0.0) | (values > self.xi_materializable_max + 2.0e-13)):
            raise ValueError("xi lies outside the current finite-X pulse prefix")
        out = np.zeros_like(values)
        flat = values.reshape(-1)
        out_flat = out.reshape(-1)
        lam = self.lambda_value
        for i, upper in enumerate(flat):
            if upper <= 0.0:
                continue
            z = 0.5 * upper * (_GL_NODES + 1.0)
            exponent = ((0.5 - lam) * z - upper) / lam
            integrand = np.exp(exponent) * main_kernel_R0(z)
            out_flat[i] = 0.5 * upper * float(np.dot(_GL_WEIGHTS, integrand))
        return out

    def _pulse_profile(
        self, X: np.ndarray, eta: np.ndarray
    ) -> dict[str, np.ndarray]:
        E_entry, E_entry_eta, m_entry, m_eta_entry = self._pulse_entry_state(eta)
        k = self.kernel.profile_values(
            X,
            eta,
            X_p=self.X_p,
            E_entry=E_entry,
            lambda_value=self.lambda_value,
        )
        xi = np.asarray(k["xi"], dtype=float)
        norm_int = self._normalized_pulse_integral(xi)
        amp_over_lambda = self.kernel.A_principal / self.lambda_value
        delta_m = E_entry * amp_over_lambda * norm_int
        delta_m_eta = E_entry_eta * amp_over_lambda * norm_int
        carry = self.X_p / np.asarray(X, dtype=float)
        m_ratio = carry * m_entry + delta_m
        m_eta_ratio = carry * m_eta_entry + delta_m_eta
        U = np.asarray(k["U_source_main_pulse_principal"], dtype=float)

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
            **k,
            "E_entry": E_entry,
            "E_entry_eta": E_entry_eta,
            "M_over_X_entry": m_entry,
            "M_eta_over_X_entry": m_eta_entry,
            "delta_M_over_X_main_pulse": delta_m,
            "delta_M_eta_over_X_main_pulse": delta_m_eta,
            "M_over_X_current_main_pulse": m_ratio,
            "M_eta_over_X_current_main_pulse": m_eta_ratio,
            "v0_current_main_pulse": v0,
        }

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
        v0 = np.empty_like(xf)
        delta_m = np.zeros_like(xf)
        delta_m_eta = np.zeros_like(xf)
        xi = np.zeros_like(xf)
        region = np.empty(xf.shape, dtype=object)

        inherited = xf <= self.X_p
        if np.any(inherited):
            p = self.parent.similarity_profile_values(xf[inherited], ef[inherited])
            F[inherited] = np.asarray(p["F_current_leading_to_pulse_entry"], dtype=float)
            U[inherited] = np.asarray(p["U_current_leading_to_pulse_entry"], dtype=float)
            E[inherited] = np.asarray(p["E_current_leading_to_pulse_entry"], dtype=float)
            m_ratio[inherited] = np.asarray(
                p["M_over_X_current_leading_to_pulse_entry"], dtype=float
            )
            m_eta_ratio[inherited] = np.asarray(
                p["M_eta_over_X_current_leading_to_pulse_entry"], dtype=float
            )
            v0[inherited] = np.asarray(
                p["v0_current_leading_to_pulse_entry"], dtype=float
            )
            region[inherited] = np.asarray(p["region"], dtype=object)

        pulse = ~inherited
        if np.any(pulse):
            p = self._pulse_profile(xf[pulse], ef[pulse])
            F[pulse] = np.asarray(p["F_source_main_pulse_principal"], dtype=float)
            U[pulse] = np.asarray(p["U_source_main_pulse_principal"], dtype=float)
            E[pulse] = np.asarray(p["E_source_main_pulse_principal"], dtype=float)
            m_ratio[pulse] = np.asarray(p["M_over_X_current_main_pulse"], dtype=float)
            m_eta_ratio[pulse] = np.asarray(
                p["M_eta_over_X_current_main_pulse"], dtype=float
            )
            v0[pulse] = np.asarray(p["v0_current_main_pulse"], dtype=float)
            delta_m[pulse] = np.asarray(p["delta_M_over_X_main_pulse"], dtype=float)
            delta_m_eta[pulse] = np.asarray(
                p["delta_M_eta_over_X_main_pulse"], dtype=float
            )
            xi[pulse] = np.asarray(p["xi"], dtype=float)
            region[pulse] = "current_cartesian_main_pulse_principal_prefix"

        arrays = (F, U, E, m_ratio, m_eta_ratio, v0)
        if any(np.any(~np.isfinite(arr)) for arr in arrays):
            raise RuntimeError("current Cartesian main-pulse candidate produced non-finite values")
        if np.any(F <= 0.0) or np.any(E < 0.0):
            raise RuntimeError("current Cartesian main-pulse candidate lost positive F/E")
        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "region": region.reshape(shape),
            "xi_current_main_pulse": xi.reshape(shape),
            "F_current_leading_with_main_pulse": F.reshape(shape),
            "U_current_leading_with_main_pulse": U.reshape(shape),
            "E_current_leading_with_main_pulse": E.reshape(shape),
            "M_over_X_current_leading_with_main_pulse": m_ratio.reshape(shape),
            "M_eta_over_X_current_leading_with_main_pulse": m_eta_ratio.reshape(shape),
            "v0_current_leading_with_main_pulse": v0.reshape(shape),
            "delta_M_over_X_main_pulse": delta_m.reshape(shape),
            "delta_M_eta_over_X_main_pulse": delta_m_eta.reshape(shape),
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

        inherited = xf <= self.X_p
        if np.any(inherited):
            p = self.parent.similarity_radial_derivatives(xf[inherited], ef[inherited])
            F_X[inherited] = np.asarray(
                p["F_current_leading_to_pulse_entry_X"], dtype=float
            )
            U_X[inherited] = np.asarray(
                p["U_current_leading_to_pulse_entry_X"], dtype=float
            )
            E_X[inherited] = np.asarray(
                p["E_current_leading_to_pulse_entry_X"], dtype=float
            )

        pulse = ~inherited
        if np.any(pulse):
            E_entry, _, _, _ = self._pulse_entry_state(ef[pulse])
            p = self.kernel.radial_derivatives(
                xf[pulse],
                ef[pulse],
                X_p=self.X_p,
                E_entry=E_entry,
                lambda_value=self.lambda_value,
            )
            F_X[pulse] = np.asarray(
                p["F_X_source_main_pulse_principal"], dtype=float
            )
            U_X[pulse] = np.asarray(
                p["U_X_source_main_pulse_principal"], dtype=float
            )
            E_X[pulse] = np.asarray(
                p["E_X_source_main_pulse_principal"], dtype=float
            )

        return {
            "X": xf.reshape(shape),
            "eta": ef.reshape(shape),
            "F_current_leading_with_main_pulse_X": F_X.reshape(shape),
            "U_current_leading_with_main_pulse_X": U_X.reshape(shape),
            "E_current_leading_with_main_pulse_X": E_X.reshape(shape),
        }

    def velocity(self, x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
        coords = self.similarity_coordinates(x, y, z, t)
        X = np.asarray(coords["X"], dtype=float)
        tol = 128.0 * np.finfo(float).eps * self.X_materializable_end
        if np.any(X > self.X_materializable_end + tol):
            raise ValueError("Cartesian point lies beyond the finite-X current main-pulse domain")
        X = np.minimum(X, self.X_materializable_end)
        eta = np.asarray(coords["eta"], dtype=float)
        p = self.similarity_profile_values(X, eta)

        xb, yb, qb = np.broadcast_arrays(
            _finite(x, "x"), _finite(y, "y"), np.asarray(coords["q"], dtype=float)
        )
        F = np.asarray(p["F_current_leading_with_main_pulse"], dtype=float)
        U = np.asarray(p["U_current_leading_with_main_pulse"], dtype=float)
        v0 = np.asarray(p["v0_current_leading_with_main_pulse"], dtype=float)
        radial = v0 / (2.0 * qb)
        swirl = qb ** (-self.A - 0.5) * F
        axial = qb ** (-self.A) * U
        return np.stack(
            (radial * xb - swirl * yb, radial * yb + swirl * xb, axial), axis=-1
        )

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_pulse_entry": self.parent.configuration(),
            "main_pulse_kernel": self.kernel.to_configuration(),
            "bound_scope": {
                "source_commit": SOURCE_COMMIT,
                "source_role": "current_cartesian_main_pulse_principal_finite_X_prefix",
                "amplitude_role": "repository_autonomous_principal_not_source_exact",
                "full_source_xi_11_current_cartesian_materialized": False,
                "source_pulse_end_MJ_corrections_materialized": False,
                "float_X_safety_divisor": _FLOAT_X_SAFETY_DIVISOR,
                "primitive_quadrature_order": QUADRATURE_ORDER,
                "mutable": False,
            },
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA16CurrentCartesianMainPulse":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected current Cartesian main-pulse schema")
        parent_payload = payload.get("parent_pulse_entry")
        kernel_payload = payload.get("main_pulse_kernel")
        bound = payload.get("bound_scope")
        if not all(isinstance(v, Mapping) for v in (parent_payload, kernel_payload, bound)):
            raise ValueError("missing current pulse parent, kernel, or scope binding")
        expected_bound = {
            "source_commit": SOURCE_COMMIT,
            "source_role": "current_cartesian_main_pulse_principal_finite_X_prefix",
            "amplitude_role": "repository_autonomous_principal_not_source_exact",
            "full_source_xi_11_current_cartesian_materialized": False,
            "source_pulse_end_MJ_corrections_materialized": False,
            "float_X_safety_divisor": _FLOAT_X_SAFETY_DIVISOR,
            "primitive_quadrature_order": QUADRATURE_ORDER,
            "mutable": False,
        }
        if dict(bound) != expected_bound:
            raise ValueError("serialized current main-pulse source/truth binding changed")
        parent = KokunoPA16CurrentCartesianPostI4PulseEntryBridge.from_configuration(
            parent_payload
        )
        kernel = KokunoSourceMainAxialPulseKernel.from_configuration(kernel_payload)
        return cls(parent=parent, kernel=kernel)

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
    ) -> "KokunoPA16CurrentCartesianMainPulse":
        return cls.from_configuration(
            json.loads(Path(path).read_text(encoding="utf-8"))
        )

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "kernel_semantic_sha256": self.kernel.semantic_sha256,
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
                "xi_materializable_max": self.xi_materializable_max,
                "log_X_materializable_end": self.log_X_materializable_end,
            },
            "truth_boundary": _TRUTH_BOUNDARY,
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return dict(_TRUTH_BOUNDARY)

    def receipt(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_exact_head": PARENT_EXACT_HEAD,
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "kernel_semantic_sha256": self.kernel.semantic_sha256,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": SOURCE_RELEASE,
                "release_date": SOURCE_RELEASE_DATE,
            },
            "domain": {
                "X_p": self.X_p,
                "log_X_p": self.log_X_p,
                "lambda": self.lambda_value,
                "source_main_xi_max": MAIN_XI_MAX,
                "xi_materializable_max": self.xi_materializable_max,
                "X_materializable_end": self.X_materializable_end,
                "log_X_materializable_end": self.log_X_materializable_end,
            },
            "principal_amplitude": {
                "K_b": self.kernel.K_b,
                "A_principal": self.kernel.A_principal,
                "origin": "repository_autonomous_principal_spin_term_not_source_exact",
            },
            "truth_boundary": self.truth_boundary,
        }
