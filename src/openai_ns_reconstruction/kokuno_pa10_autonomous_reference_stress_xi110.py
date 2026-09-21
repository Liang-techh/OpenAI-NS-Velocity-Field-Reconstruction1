"""Extend the current Kokuno PA.10 reference stress from X=100 to X_i=110.

Pinned public provenance: ``KokunoYumeto/yang-mills-interacting-workbench`` at
``143f6773feb424ad9ed3a8d116653200f20346b7``, corrected 2026-09-09 reader.

The public Appendix-B continuation makes the reference ``phi_r`` and ``U_r``
constant in ``y=log(X/X_0)`` after the first reference collar, while pressure
and the prefix primitives continue from the axis by their original integral
definitions.  In particular, on the already-flat reference interval,

    F_r(X,eta) = F_r(100,eta),
    U_r(X,eta) = U_r(100,eta),
    E_r = sqrt(2X) F_r,
    Pi_X = F_r^2,
    M(X,eta) = M(100,eta) + (X-100) U_r(100,eta),
    ell_r = X partial_X log(2 X F_r) = 1.

The corrected reconstruction keeps using the reference stress primitives in the
actual final bridge.  They obey

    D_X p_{1,r} = X S_{q,r}/L - ell_r p_{1,r},
    D_X N_{s,r} + N_{s,r} = S_{n,r},
    n_{s,r} = N_{s,r}/L,
    p_{2,r} = X n_{s,r}/E_r,

with ``D_X=X partial_X`` and the public ``S_q``/``S_n`` formulas.

This module is deliberately a prerequisite seam only: it extends the *current*
repository-autonomous-pressure reference stress from Agent-1 #918 through the
public ``X_i=110``.  It does not execute the actual ``100<X<110`` shutdown /
``a -> 0.8`` bridge.  The absolute axis pressure inherited from #913 remains a
repository-autonomous admissible realization, not the Appendix-A source datum.
No residual, forcing, image fit, amplitude collapse, or PDE threshold is used.
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

from .kokuno_pa10_autonomous_reference_ns_x100 import (
    KokunoPA10AutonomousReferenceNsToX100,
)
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa10_reference_p1_x100 import ETA_FD_STEP


SCHEMA = "kokuno-pa10-autonomous-reference-stress-xi110-v1"
SOURCE_X_B_REF = 100.0
SOURCE_X_I = 110.0
STRESS_QUADRATURE_ORDER = 64

_SOURCE_FORMULAS = {
    "reference_plateau": (
        "after the first reference collar phi_r,U_r are constant in y=log(X/X_0); "
        "therefore F_r,U_r are X-constant and ell_r=1"
    ),
    "reference_pressure": "Pi_X=F_r^2",
    "reference_prefix": "M(X)=M(100)+(X-100)U_r(100)",
    "angular_source": "S_q=-W ell-h(1-2eta U)-H_c (log E)_eta",
    "axial_source": (
        "S_n=-W X U_X-A(1-2eta U)U-H_c U_eta-d Pi_eta"
        "+4A eta Pi+2eta X Pi_X"
    ),
    "p1_ode": "D_X p_{1,r}=X S_{q,r}/L-ell_r p_{1,r}",
    "Ns_ode": "D_X N_{s,r}+N_{s,r}=S_{n,r}",
    "stress_definitions": "n_{s,r}=N_{s,r}/L; p_{2,r}=X n_{s,r}/E_r",
    "final_reference_domain": "X_b^ref=100 through X_i=110",
}

_NUMERICAL_REALIZATION = {
    "pressure_dependency": (
        "inherits Agent-1 #913 autonomous admissible Pi0 seed through #918; "
        "not source-prepared Appendix-A Pi0"
    ),
    "eta_derivative": f"fixed centered FD4 with h_eta={ETA_FD_STEP:.17g}",
    "stress_radial_integral": (
        f"fixed Gauss-Legendre order {STRESS_QUADRATURE_ORDER} on [100,X]"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_reference_plateau_formula_executable_to_Xi": True,
    "public_reference_stress_odes_executable_to_Xi": True,
    "autonomous_pressure_reference_stress_pair_to_Xi_executable": True,
    "candidate_side_reference_Xi_handoff_executable": True,
    "source_prepared_appendixA_Pi0_materialized": False,
    "source_prepared_reference_nsr_materialized": False,
    "source_prepared_full_reference_stress_pair_materialized": False,
    "source_admitted_global_kappa0_smallness": False,
    "actual_final_interpolation_100_to_Xi_materialized": False,
    "actual_G_i_at_Xi_materialized": False,
    "actual_ell_i_at_Xi_materialized": False,
    "actual_upstream_five_moment_discrepancy_materialized": False,
    "actual_inner_to_outer_bridge_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
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


@lru_cache(maxsize=4)
def _legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(order, int) or order < 16:
        raise ValueError("quadrature order must be an integer >=16")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


@dataclass(frozen=True)
class KokunoPA10AutonomousReferenceStressToXi110:
    """Continue the current autonomous-pressure reference stress to ``X_i=110``."""

    parent: KokunoPA10AutonomousReferenceNsToX100 = field(
        default_factory=KokunoPA10AutonomousReferenceNsToX100,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPA10AutonomousReferenceNsToX100):
            raise TypeError("parent must be KokunoPA10AutonomousReferenceNsToX100")
        if abs(self.parent.X_b_ref - SOURCE_X_B_REF) > 0.0:
            raise ValueError("parent reference cutoff must be exactly X_b^ref=100")
        if not SOURCE_X_I > SOURCE_X_B_REF:
            raise ValueError("source Xi must lie beyond the reference cutoff")

    @property
    def X_b_ref(self) -> float:
        return SOURCE_X_B_REF

    @property
    def X_i(self) -> float:
        return SOURCE_X_I

    @property
    def h(self) -> float:
        return float(self.parent.h)

    @property
    def A(self) -> float:
        return float(self.parent.A)

    @property
    def D(self) -> float:
        return float(self.parent.D)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.parent.eta_interval

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        x_hi = math.nextafter(self.X_i, math.inf)
        if np.any(X_arr < self.X_b_ref) or np.any(X_arr > x_hi):
            raise ValueError(
                f"X must lie in the reference-stress Xi extension [{self.X_b_ref},{self.X_i}]"
            )
        margin = 2.0 * ETA_FD_STEP
        lo, hi = self.eta_interval
        if np.any(eta_arr < lo + margin) or np.any(eta_arr > hi - margin):
            raise ValueError(
                "eta lies too close to the executable boundary for the frozen centered eta derivative"
            )
        return np.minimum(X_arr, self.X_i), eta_arr

    def _eta_derivative(self, evaluator, eta: np.ndarray) -> np.ndarray:
        h = ETA_FD_STEP
        return (
            -evaluator(eta + 2.0 * h)
            + 8.0 * evaluator(eta + h)
            - 8.0 * evaluator(eta - h)
            + evaluator(eta - 2.0 * h)
        ) / (12.0 * h)

    def _x100_data(self, eta: np.ndarray) -> dict[str, np.ndarray]:
        """Return all frozen reference data needed to continue beyond X=100."""
        X100 = np.full_like(eta, self.X_b_ref, dtype=float)
        reference = self.parent.reference.values(X100, eta)
        kinematics = self.parent.p1_reference.source_terms(X100, eta)
        ns = self.parent.values(X100, eta)
        p1 = self.parent.p1_reference.values(X100, eta)["p1_reference"]
        Pi = self.parent.pressure.pressure(X100, eta)
        Pi_eta = self.parent.pressure.eta_derivative(X100, eta)

        def F100(e: np.ndarray) -> np.ndarray:
            x = np.full_like(e, self.X_b_ref, dtype=float)
            return self.parent.reference.values(x, e)["F_reference"]

        def U100(e: np.ndarray) -> np.ndarray:
            x = np.full_like(e, self.X_b_ref, dtype=float)
            return self.parent.reference.values(x, e)["U_reference"]

        F_eta = self._eta_derivative(F100, eta)
        U_eta = self._eta_derivative(U100, eta)
        return {
            "F100": reference["F_reference"],
            "U100": reference["U_reference"],
            "F100_eta": F_eta,
            "U100_eta": U_eta,
            "M100": kinematics["M_reference"],
            "M100_eta": kinematics["M_reference_eta"],
            "Pi100": Pi,
            "Pi100_eta": Pi_eta,
            "p1_100": p1,
            "Ns100": ns["N_s_reference_autonomous"],
            "ns100": ns["n_s_reference_autonomous"],
            "p2_100": ns["p2_reference_autonomous"],
        }

    def reference_state(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate the public frozen-reference kinematics/sources on [100,110]."""
        X_arr, eta_arr = self._broadcast(X, eta)
        data = self._x100_data(eta_arr)
        F = data["F100"]
        U = data["U100"]
        F_eta = data["F100_eta"]
        U_eta = data["U100_eta"]
        delta = X_arr - self.X_b_ref
        M = data["M100"] + delta * U
        M_eta = data["M100_eta"] + delta * U_eta
        Pi = data["Pi100"] + delta * F * F
        Pi_X = F * F
        Pi_eta = data["Pi100_eta"] + 2.0 * delta * F * F_eta
        E = np.sqrt(2.0 * X_arr) * F
        if np.any(F <= 0.0) or np.any(E <= 0.0):
            raise RuntimeError("reference F/E must remain positive through Xi")

        d = 1.0 - eta_arr * eta_arr
        L = 1.0 - 2.0 * self.h * eta_arr * eta_arr
        if np.any(L <= 0.0):
            raise RuntimeError("source L lost positivity on the executable eta interval")
        W = 1.0 - (2.0 * self.D * eta_arr * M + d * M_eta) / X_arr
        H_c = self.D * eta_arr + d * U
        ell = np.ones_like(X_arr, dtype=float)
        log_E_eta = F_eta / F
        S_q = -W * ell - self.h * (1.0 - 2.0 * eta_arr * U) - H_c * log_E_eta
        S_n = (
            -self.A * (1.0 - 2.0 * eta_arr * U) * U
            - H_c * U_eta
            - d * Pi_eta
            + 4.0 * self.A * eta_arr * Pi
            + 2.0 * eta_arr * X_arr * Pi_X
        )
        H = 2.0 * X_arr * F
        if np.any(~np.isfinite(S_q)) or np.any(~np.isfinite(S_n)):
            raise RuntimeError("reference Xi source propagation produced non-finite terms")
        return {
            "X": X_arr,
            "eta": eta_arr,
            "F_reference": F,
            "U_reference": U,
            "E_reference": E,
            "F_reference_eta": F_eta,
            "U_reference_eta": U_eta,
            "M_reference": M,
            "M_reference_eta": M_eta,
            "Pi_reference_autonomous": Pi,
            "Pi_reference_autonomous_X": Pi_X,
            "Pi_reference_autonomous_eta": Pi_eta,
            "d": d,
            "L": L,
            "W_reference": W,
            "H_c_reference": H_c,
            "ell_reference": ell,
            "log_E_reference_eta": log_E_eta,
            "H_reference": H,
            "S_q_reference": S_q,
            "S_n_reference_autonomous": S_n,
        }

    def values(
        self, X: Any, eta: Any, *, order: int = STRESS_QUADRATURE_ORDER
    ) -> dict[str, np.ndarray]:
        """Return reference ``p1,N_s,n_s,p2`` and frozen profiles through Xi."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        start = self._x100_data(e)
        nodes, weights = _legendre_rule(order)
        half = 0.5 * (x - self.X_b_ref)
        mid = 0.5 * (x + self.X_b_ref)
        X_nodes = mid[:, None] + half[:, None] * nodes[None, :]
        eta_nodes = np.broadcast_to(e[:, None], X_nodes.shape)
        terms = self.reference_state(X_nodes, eta_nodes)

        p1_integrand = terms["H_reference"] * terms["S_q_reference"] / terms["L"]
        p1_integral = half * np.sum(weights[None, :] * p1_integrand, axis=1)
        H100 = 2.0 * self.X_b_ref * start["F100"]
        H = 2.0 * x * start["F100"]
        p1 = (H100 * start["p1_100"] + p1_integral) / H

        Ns_integral = half * np.sum(
            weights[None, :] * terms["S_n_reference_autonomous"], axis=1
        )
        Ns = (self.X_b_ref * start["Ns100"] + Ns_integral) / x
        L = 1.0 - 2.0 * self.h * e * e
        ns = Ns / L
        E = np.sqrt(2.0 * x) * start["F100"]
        p2 = x * ns / E
        Pi = start["Pi100"] + (x - self.X_b_ref) * start["F100"] ** 2
        if any(np.any(~np.isfinite(v)) for v in (p1, Ns, ns, p2, Pi)):
            raise RuntimeError("reference stress Xi continuation became non-finite")
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "F_reference": start["F100"].reshape(shape),
            "U_reference": start["U100"].reshape(shape),
            "E_reference": E.reshape(shape),
            "Pi_reference_autonomous": Pi.reshape(shape),
            "p1_reference": p1.reshape(shape),
            "N_s_reference_autonomous": Ns.reshape(shape),
            "n_s_reference_autonomous": ns.reshape(shape),
            "p2_reference_autonomous": p2.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return first radial derivatives directly from the public reference ODEs."""
        X_arr, eta_arr = self._broadcast(X, eta)
        values = self.values(X_arr, eta_arr)
        terms = self.reference_state(X_arr, eta_arr)
        p1 = values["p1_reference"]
        Ns = values["N_s_reference_autonomous"]
        p1_X = terms["S_q_reference"] / terms["L"] - p1 / X_arr
        Ns_X = (terms["S_n_reference_autonomous"] - Ns) / X_arr
        ns_X = Ns_X / terms["L"]
        F = values["F_reference"]
        return {
            "F_reference_X": np.zeros_like(X_arr, dtype=float),
            "U_reference_X": np.zeros_like(X_arr, dtype=float),
            "E_reference_X": F / np.sqrt(2.0 * X_arr),
            "Pi_reference_autonomous_X": F * F,
            "p1_reference_X": p1_X,
            "N_s_reference_autonomous_X": Ns_X,
            "n_s_reference_autonomous_X": ns_X,
        }

    def handoff_at_Xi(self, eta: Any) -> dict[str, np.ndarray]:
        eta_arr = _finite(eta, "eta")
        X = np.full_like(eta_arr, self.X_i, dtype=float)
        return {
            **self.values(X, eta_arr),
            **self.radial_derivatives(X, eta_arr),
            **self.reference_state(X, eta_arr),
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_semantic_sha256": self.parent.semantic_sha256,
            "source_X_b_ref": self.X_b_ref,
            "source_X_i": self.X_i,
            "eta_fd_step": ETA_FD_STEP,
            "stress_quadrature_order": STRESS_QUADRATURE_ORDER,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10AutonomousReferenceStressToXi110":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError("reference-stress Xi configuration is frozen; serialized data differ")
        return cls()

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoPA10AutonomousReferenceStressToXi110":
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
        eta_probe = np.asarray([-0.70, -0.25, 0.0, 0.25, 0.70])
        start = self.values(np.full_like(eta_probe, self.X_b_ref), eta_probe)
        end = self.handoff_at_Xi(eta_probe)
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
            "geometry": {"X_b_ref": self.X_b_ref, "X_i": self.X_i},
            "X100_probe": {
                "eta": eta_probe.tolist(),
                "p1": start["p1_reference"].tolist(),
                "n_s": start["n_s_reference_autonomous"].tolist(),
            },
            "Xi_probe": {
                "eta": eta_probe.tolist(),
                "F_reference": end["F_reference"].tolist(),
                "U_reference": end["U_reference"].tolist(),
                "p1": end["p1_reference"].tolist(),
                "n_s": end["n_s_reference_autonomous"].tolist(),
                "p2": end["p2_reference_autonomous"].tolist(),
            },
            "semantic_sha256": self.semantic_sha256,
            "truth_boundary": self.truth_boundary,
        }
