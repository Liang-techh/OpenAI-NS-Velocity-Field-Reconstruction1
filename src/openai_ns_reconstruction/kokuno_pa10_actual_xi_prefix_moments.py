"""Current-lineage Kokuno PA.10 prefix moments through ``X_i=110``.

This module binds the executable Agent-1 actual profile chain

    natural PA.10 -> stress activation -> fixed-kappa continuation
      -> final bridge -> X_i=110

behind one low-dimensional ``F(X,eta), U(X,eta), E(X,eta)`` evaluator on
``0 <= X <= X_i`` and propagates the five public prefix moments

    M   = int U dX,
    I   = int H dX,                  H = sqrt(2X) E = 2 X F,
    J   = int U H dX,
    S   = int (U^2-E^2/2) dX,
    C_p = int E^2/(2X) dX = int F^2 dX.

The corrected KokunoYumeto 2026-09-09 reconstruction then uses the PA.15
normalization

    (M,I,J,S,C_p)
      = (X_R Mhat, X_R^(3/2) Ihat, X_R^(3/2) Jhat,
         X_R Shat, Cphat),

with ``X_R=110(CP_*)^10`` and ideal pair

    U_0 = 4 eta,
    E_0 = P_* (1+eta^2)^(-1) (X/X_R)^(1/10).

The source publishes the scale hierarchy but not one hidden numerical outer
choice.  This module therefore uses ``KokunoOuterReservedPatchSchedule`` with
the *same* source-normalization ``C`` as the current A1 chain and its explicit
repository-autonomous source-hierarchy values for the remaining outer
parameters.  No hidden OpenAI/Kokuno parameter is recovered.

The resulting discrepancy is an executable *candidate-side* upstream PA.15
input at ``X_i``.  It is not the source-prepared Appendix-B fixed point:
current kappa_0 and reference pressure/stress are still autonomous.  T_sh,
PA.16 repair, the inner-to-outer join, global Cartesian velocity, matched
pressure/forcing and held-out Navier-Stokes validation remain separate and
fail closed.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_outer_reserved_patch_schedule import KokunoOuterReservedPatchSchedule
from .kokuno_pa10_actual_final_bridge_xi110 import KokunoPA10ActualFinalBridgeToXi110
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa10_source_c_normalized_physical_center import SELECTED_SOURCE_C


SCHEMA = "kokuno-pa10-actual-xi-prefix-moments-v1"
PREFIX_QUADRATURE_ORDER = 48

_SOURCE_FORMULAS = {
    "five_prefix_moments": (
        "M=int U dX; I=int H dX; J=int U H dX; "
        "S=int(U^2-E^2/2)dX; C_p=int E^2/(2X)dX; H=sqrt(2X)E=2XF"
    ),
    "PA15_scaling": (
        "(M,I,J,S,C_p)=(X_R*Mhat,X_R^(3/2)*Ihat,"
        "X_R^(3/2)*Jhat,X_R*Shat,Cphat)"
    ),
    "outer_scale": "X_R=110(C P_*)^10; x=X/X_R",
    "ideal_pair": "U_0=4 eta; E_0=P_* f x^(1/10); f=(1+eta^2)^(-1)",
    "incoming_discrepancy": (
        "current-lineage actual PA15 prefix at X_i minus ideal-pair PA15 prefix at X_i"
    ),
}

_NUMERICAL_REALIZATION = {
    "outer_schedule": (
        "repository-autonomous KokunoOuterReservedPatchSchedule with the same C "
        "as the current A1 source-normalized profile; no hidden-number recovery"
    ),
    "prefix_quadrature": (
        f"fixed Gauss-Legendre order {PREFIX_QUADRATURE_ORDER} independently on "
        "[0,X_0], [X_0,X_1], [X_1,100], [100,110]"
    ),
    "profile_chain": (
        "public executable A1 profile values from the current #940 ancestry; "
        "no residual or forcing enters the moment calculation"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "current_lineage_profile_0_to_Xi_executable": True,
    "current_lineage_radial_derivatives_executable": True,
    "public_five_prefix_moment_map_executable": True,
    "public_PA15_scaling_executable": True,
    "candidate_side_upstream_five_moment_discrepancy_at_Xi_materialized": True,
    "candidate_side_PA16_input_tuple_materialized": True,
    "outer_scale_is_repository_autonomous_source_hierarchy_realization": True,
    "source_hidden_numeric_choices_recovered": False,
    "source_prepared_appendixA_pressure_stress_materialized": False,
    "source_admitted_kappa0_materialized": False,
    "source_prepared_upstream_five_moment_discrepancy_materialized": False,
    "source_T_sh_lower_bound_verified": False,
    "five_moment_repair_applied": False,
    "inner_to_outer_join_completed": False,
    "outer_global_leading_velocity_materialized": False,
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


@lru_cache(maxsize=8)
def _legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(order, int) or order < 16:
        raise ValueError("quadrature order must be an integer >=16")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


def _default_outer_schedule() -> KokunoOuterReservedPatchSchedule:
    return KokunoOuterReservedPatchSchedule(C=SELECTED_SOURCE_C)


@dataclass(frozen=True)
class KokunoPA10ActualXiPrefixMoments:
    """Unify the current A1 profile and compute its PA.15 prefix discrepancy."""

    bridge: KokunoPA10ActualFinalBridgeToXi110 = field(
        default_factory=KokunoPA10ActualFinalBridgeToXi110,
        repr=False,
        compare=False,
    )
    outer_schedule: KokunoOuterReservedPatchSchedule = field(
        default_factory=_default_outer_schedule,
        repr=False,
        compare=False,
    )
    quadrature_order: int = PREFIX_QUADRATURE_ORDER

    def __post_init__(self) -> None:
        if not isinstance(self.bridge, KokunoPA10ActualFinalBridgeToXi110):
            raise TypeError("bridge must be KokunoPA10ActualFinalBridgeToXi110")
        if not isinstance(self.outer_schedule, KokunoOuterReservedPatchSchedule):
            raise TypeError("outer_schedule must be KokunoOuterReservedPatchSchedule")
        order = int(self.quadrature_order)
        if isinstance(self.quadrature_order, bool) or order != self.quadrature_order:
            raise TypeError("quadrature_order must be an integer")
        if not 16 <= order <= 128:
            raise ValueError("quadrature_order must lie in [16,128]")
        object.__setattr__(self, "quadrature_order", order)
        if not math.isclose(
            float(self.outer_schedule.C), self.C, rel_tol=0.0, abs_tol=1e-14
        ):
            raise ValueError("outer schedule C must equal the current A1 source normalization C")
        min_log = math.log(np.finfo(float).tiny)
        if -1.5 * float(self.outer_schedule.log_X_R) <= min_log:
            raise ValueError(
                "outer schedule makes PA.15 X_R^(-3/2) underflow in binary64; "
                "a future log-scaled adapter is required"
            )
        if abs(self.C - SELECTED_SOURCE_C) > 0.0:
            raise ValueError("current A1 source-normalized C identity drifted")

    @property
    def C(self) -> float:
        return float(self.bridge.C)

    @property
    def X_i(self) -> float:
        return float(self.bridge.X_i)

    @property
    def X_b_ref(self) -> float:
        return float(self.bridge.X_b_ref)

    @property
    def X_0(self) -> float:
        return float(self.bridge.actual_x100.activation.X_0)

    @property
    def X_1(self) -> float:
        return float(self.bridge.actual_x100.activation.X_1)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.bridge.eta_interval)

    @property
    def log_X_R(self) -> float:
        return float(self.outer_schedule.log_X_R)

    @property
    def log_x_i(self) -> float:
        return math.log(self.X_i) - self.log_X_R

    @property
    def log_P_star(self) -> float:
        return float(self.outer_schedule.log_P_star)

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any((X_arr < 0.0) | (X_arr > self.X_i)):
            raise ValueError(f"X must lie in the current prefix domain [0,{self.X_i}]")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in the current executable interval [{lo},{hi}]")
        return X_arr, eta_arr

    def profile_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return vectorized current-lineage ``F,U,E`` on ``0<=X<=X_i``."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        F = np.empty_like(x)
        U = np.empty_like(x)

        natural = self.bridge.actual_x100.activation.reference.source_normalized
        activation = self.bridge.actual_x100.activation
        fixed = self.bridge.actual_x100

        m0 = x <= self.X_0
        m1 = (x > self.X_0) & (x <= self.X_1)
        m2 = (x > self.X_1) & (x < self.X_b_ref)
        m3 = x >= self.X_b_ref

        if np.any(m0):
            vals = natural.values(x[m0], e[m0])
            F[m0] = vals["F_0"]
            U[m0] = vals["U_0"]
        if np.any(m1):
            vals = activation.values(x[m1], e[m1])
            F[m1] = vals["F_activation"]
            U[m1] = vals["U_activation"]
        if np.any(m2):
            vals = fixed.values(x[m2], e[m2])
            F[m2] = vals["F_fixed_kappa"]
            U[m2] = vals["U_fixed_kappa"]
        if np.any(m3):
            vals = self.bridge.values(x[m3], e[m3])
            F[m3] = vals["F_final_bridge"]
            U[m3] = vals["U_final_bridge"]

        E = np.sqrt(2.0 * x) * F
        if np.any(~np.isfinite(F)) or np.any(~np.isfinite(U)) or np.any(~np.isfinite(E)):
            raise RuntimeError("current prefix evaluator produced non-finite values")
        if np.any(F <= 0.0):
            raise RuntimeError("current prefix evaluator lost positive F")
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "F_actual_prefix": F.reshape(shape),
            "U_actual_prefix": U.reshape(shape),
            "E_actual_prefix": E.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return piecewise analytic/source-ODE first radial derivatives for ``X>0``."""
        X_arr, eta_arr = self._broadcast(X, eta)
        if np.any(X_arr <= 0.0):
            raise ValueError("radial_derivatives requires X>0")
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        F_X = np.empty_like(x)
        U_X = np.empty_like(x)

        natural = self.bridge.actual_x100.activation.reference.source_normalized
        activation = self.bridge.actual_x100.activation
        fixed = self.bridge.actual_x100

        m0 = x <= self.X_0
        m1 = (x > self.X_0) & (x <= self.X_1)
        m2 = (x > self.X_1) & (x < self.X_b_ref)
        m3 = x >= self.X_b_ref

        if np.any(m0):
            vals = natural.derivatives(x[m0], e[m0])
            F_X[m0] = vals["F_0_X"]
            U_X[m0] = vals["U_0_X"]
        if np.any(m1):
            vals = activation.radial_derivatives(x[m1], e[m1])
            F_X[m1] = vals["F_activation_X"]
            U_X[m1] = vals["U_activation_X"]
        if np.any(m2):
            vals = fixed.radial_derivatives(x[m2], e[m2])
            F_X[m2] = vals["F_fixed_kappa_X"]
            U_X[m2] = vals["U_fixed_kappa_X"]
        if np.any(m3):
            vals = self.bridge.radial_derivatives(x[m3], e[m3])
            F_X[m3] = vals["F_final_bridge_X"]
            U_X[m3] = vals["U_final_bridge_X"]

        profile = self.profile_values(x, e)
        F = profile["F_actual_prefix"]
        E_X = F / np.sqrt(2.0 * x) + np.sqrt(2.0 * x) * F_X
        return {
            "F_actual_prefix_X": F_X.reshape(shape),
            "U_actual_prefix_X": U_X.reshape(shape),
            "E_actual_prefix_X": E_X.reshape(shape),
        }

    def moment_densities(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return public physical five-moment densities in ``(M,I,J,S,C_p)`` order."""
        vals = self.profile_values(X, eta)
        X_arr = vals["X"]
        F = vals["F_actual_prefix"]
        U = vals["U_actual_prefix"]
        H = 2.0 * X_arr * F
        return {
            "M_density": U,
            "I_density": H,
            "J_density": U * H,
            "S_density": U * U - X_arr * F * F,
            "C_p_density": F * F,
        }

    def _integrate_segment(
        self, a: float, b: float, eta: float, *, order: int
    ) -> np.ndarray:
        if not b > a:
            raise ValueError("moment segment must have b>a")
        nodes, weights = _legendre_rule(order)
        half = 0.5 * (b - a)
        mid = 0.5 * (a + b)
        X = mid + half * nodes
        eta_nodes = np.full_like(X, eta, dtype=float)
        d = self.moment_densities(X, eta_nodes)
        matrix = np.stack(
            [
                d["M_density"],
                d["I_density"],
                d["J_density"],
                d["S_density"],
                d["C_p_density"],
            ],
            axis=0,
        )
        return half * (matrix @ weights)

    def physical_prefix_moments_at_Xi(
        self, eta: float, *, order: int | None = None
    ) -> np.ndarray:
        """Integrate the current actual five physical moments from ``0`` to ``X_i``."""
        eta_f = float(eta)
        lo, hi = self.eta_interval
        if not math.isfinite(eta_f) or not lo <= eta_f <= hi:
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        q = self.quadrature_order if order is None else int(order)
        if q < 16:
            raise ValueError("quadrature order must be >=16")
        segments = (
            (0.0, self.X_0),
            (self.X_0, self.X_1),
            (self.X_1, self.X_b_ref),
            (self.X_b_ref, self.X_i),
        )
        out = np.zeros(5, dtype=float)
        for a, b in segments:
            out += self._integrate_segment(a, b, eta_f, order=q)
        if np.any(~np.isfinite(out)):
            raise RuntimeError("current actual prefix moments are non-finite")
        return out

    def actual_PA15_scaled_moments(
        self, eta: float, *, order: int | None = None
    ) -> np.ndarray:
        physical = self.physical_prefix_moments_at_Xi(eta, order=order)
        log_X_R = self.log_X_R
        scale_1 = math.exp(-log_X_R)
        scale_3_2 = math.exp(-1.5 * log_X_R)
        scaled = np.asarray(
            [
                physical[0] * scale_1,
                physical[1] * scale_3_2,
                physical[2] * scale_3_2,
                physical[3] * scale_1,
                physical[4],
            ],
            dtype=float,
        )
        if np.any(~np.isfinite(scaled)):
            raise RuntimeError("PA.15 scaling produced non-finite current moments")
        return scaled

    def ideal_PA15_scaled_moments(self, eta: float) -> np.ndarray:
        """Return the analytic ideal-pair PA.15 prefix through ``x_i=X_i/X_R``."""
        eta_f = float(eta)
        lo, hi = self.eta_interval
        if not math.isfinite(eta_f) or not lo <= eta_f <= hi:
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        f = 1.0 / (1.0 + eta_f * eta_f)
        lx = self.log_x_i
        lp = self.log_P_star
        x = math.exp(lx)
        I_hat = (5.0 * math.sqrt(2.0) / 8.0) * f * math.exp(lp + 1.6 * lx)
        E_term = (5.0 / 12.0) * f * f * math.exp(2.0 * lp + 1.2 * lx)
        C_p_hat = 2.5 * f * f / (self.C * self.C)
        out = np.asarray(
            [
                4.0 * eta_f * x,
                I_hat,
                4.0 * eta_f * I_hat,
                16.0 * eta_f * eta_f * x - E_term,
                C_p_hat,
            ],
            dtype=float,
        )
        if np.any(~np.isfinite(out)):
            raise RuntimeError("analytic ideal PA.15 prefix is non-finite")
        return out

    def incoming_PA15_discrepancy(
        self, eta: Any, *, order: int | None = None
    ) -> np.ndarray:
        """Vectorized current candidate-side ``actual-ideal`` PA.15 discrepancy."""
        eta_arr = _finite(eta, "eta")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        out = np.empty(eta_arr.shape + (5,), dtype=float)
        flat = out.reshape(-1, 5)
        for index, value in enumerate(eta_arr.reshape(-1)):
            actual = self.actual_PA15_scaled_moments(float(value), order=order)
            ideal = self.ideal_PA15_scaled_moments(float(value))
            flat[index] = actual - ideal
        if np.any(~np.isfinite(out)):
            raise RuntimeError("incoming PA.15 discrepancy is non-finite")
        return out

    def pa16_input_at_eta(
        self, eta: float, *, order: int | None = None
    ) -> dict[str, Any]:
        """Expose the exact current data shape consumed by the existing PA.16 join.

        This method does *not* choose ``T_sh`` and does not execute the repair.
        """
        eta_f = float(eta)
        handoff = self.bridge.handoff_at_Xi(np.asarray(eta_f))
        discrepancy = self.incoming_PA15_discrepancy(np.asarray(eta_f), order=order)
        return {
            "eta": eta_f,
            "ell_i": float(np.asarray(handoff["ell_i"])),
            "G_i": float(np.asarray(handoff["G_i"])),
            "incoming_scaled_discrepancy": tuple(
                float(v) for v in np.asarray(discrepancy).reshape(5)
            ),
            "outer_schedule_sha256": self.outer_schedule.sha256,
            "source_T_sh_lower_bound_verified": False,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "bridge_semantic_sha256": self.bridge.semantic_sha256,
            "outer_schedule": self.outer_schedule.to_payload(),
            "quadrature_order": self.quadrature_order,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10ActualXiPrefixMoments":
        if not isinstance(payload, Mapping) or payload.get("schema") != SCHEMA:
            raise ValueError("unexpected actual-Xi-prefix-moments schema")
        schedule_payload = payload.get("outer_schedule")
        if not isinstance(schedule_payload, dict):
            raise ValueError("outer_schedule payload is missing")
        obj = cls(
            outer_schedule=KokunoOuterReservedPatchSchedule.from_payload(schedule_payload),
            quadrature_order=int(payload.get("quadrature_order")),
        )
        if obj.configuration() != dict(payload):
            raise ValueError("actual-Xi-prefix-moments configuration does not replay exactly")
        return obj

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA10ActualXiPrefixMoments":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_formulas": _SOURCE_FORMULAS,
            "numerical_realization": _NUMERICAL_REALIZATION,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta_probe = np.asarray([-0.5, 0.0, 0.5])
        discrepancy = self.incoming_PA15_discrepancy(eta_probe)
        pa16 = [self.pa16_input_at_eta(float(v)) for v in eta_probe]
        return {
            "schema": SCHEMA,
            "source": {
                "repository": SOURCE_REPOSITORY,
                "commit": SOURCE_COMMIT,
                "path": SOURCE_PATH,
                "release": CORRECTED_RELEASE,
                "release_date": CORRECTED_RELEASE_DATE,
            },
            "source_formulas": copy.deepcopy(_SOURCE_FORMULAS),
            "numerical_realization": copy.deepcopy(_NUMERICAL_REALIZATION),
            "configuration": self.configuration(),
            "geometry": {
                "X_0": self.X_0,
                "X_1": self.X_1,
                "X_b_ref": self.X_b_ref,
                "X_i": self.X_i,
                "log_X_R": self.log_X_R,
                "log_x_i": self.log_x_i,
            },
            "eta_probe": eta_probe.tolist(),
            "incoming_PA15_discrepancy": discrepancy.tolist(),
            "pa16_inputs": pa16,
            "semantic_sha256": self.semantic_sha256,
            "truth_boundary": self.truth_boundary,
        }
