"""Autonomous-pressure continuation of the public PA.10 reference ``n_s`` to X=100.

Pinned public provenance: ``KokunoYumeto/yang-mills-interacting-workbench`` at
``143f6773feb424ad9ed3a8d116653200f20346b7``, corrected 2026-09-09 reader.

For a leading profile the corrected reconstruction defines

    d = 1-eta^2,                 L = 1-2 h eta^2,
    A = 1/2+h,                  D = 1/2-h,
    W = 1-(2 D eta M + d M_eta)/X,
    H_c = D eta + d U,

and the axial source

    S_n = -W X U_X - A(1-2 eta U)U - H_c U_eta
          - d Pi_eta + 4 A eta Pi + 2 eta X Pi_X.

The regular stress primitive obeys

    D_X N_s + N_s = S_n,        D_X = X partial_X,
    n_s = N_s/L,

and on the natural stress-free segment the source gives
``n_s=-2 U_X``.  This module propagates that source ODE from the already
materialized natural/reference handoff ``X_0=4/Lambda`` through ``X=100``.

The pressure used here is deliberately the Agent-1 *repository-autonomous*
admissible axis-pressure seed from #913 plus the source-determined pressure
increment.  Therefore the resulting ``n_s`` is an executable candidate-side
continuation, not the Appendix-A source-prepared ``n_s`` and not matched/global
pressure evidence.  No residual, held-out sample, forcing, image fit, or PDE
threshold is used to choose any parameter in this module.

Production ``U_eta`` uses one frozen centered fourth-order finite difference,
and the radial N_s propagation uses fixed Gauss--Legendre quadrature.  Those
are repository numerical realizations; the displayed ``S_n`` formula, stress
ODE, ``n_s=N_s/L``, and natural stress-free initial identity are public-source
formulas.
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

from .kokuno_pa10_autonomous_axis_pressure_seed_x100 import (
    KokunoPA10AutonomousAxisPressureSeedToX100,
)
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa10_reference_p1_x100 import (
    ETA_FD_STEP,
    KokunoPA10ReferenceP1ToX100,
)


SCHEMA = "kokuno-pa10-autonomous-reference-ns-x100-v1"
NS_QUADRATURE_ORDER = 64

_SOURCE_FORMULAS = {
    "kinematics": (
        "d=1-eta^2; L=1-2h eta^2; A=1/2+h; D=1/2-h; "
        "W=1-(2D eta M+d M_eta)/X; H_c=D eta+dU"
    ),
    "axial_source": (
        "S_n=-W X U_X-A(1-2eta U)U-H_c U_eta-d Pi_eta"
        "+4A eta Pi+2eta X Pi_X"
    ),
    "Ns_ode": "D_X N_s+N_s=S_n; D_X=X partial_X",
    "ns_definition": "n_s=N_s/L",
    "stress_free_initial": "n_s(X_0,eta)=-2 U_X(X_0,eta)",
    "p2_definition": "p_2=X n_s/E",
}

_NUMERICAL_REALIZATION = {
    "pressure_dependency": (
        "Agent-1 autonomous admissible Pi0 seed plus source-determined C_p; "
        "not source-prepared Appendix-A Pi0"
    ),
    "U_eta": f"fixed centered FD4 with h_eta={ETA_FD_STEP:.17g}",
    "Ns_radial_integral": f"fixed Gauss-Legendre order {NS_QUADRATURE_ORDER}",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_Sn_formula_executable": True,
    "public_Ns_ode_executable_to_X100": True,
    "autonomous_pressure_reference_ns_to_X100_executable": True,
    "autonomous_pressure_reference_p2_to_X100_executable": True,
    "autonomous_pressure_reference_stress_pair_p1r_nsr_executable": True,
    "source_prepared_appendixA_Pi0_materialized": False,
    "source_prepared_reference_nsr_materialized": False,
    "source_prepared_full_reference_stress_pair_materialized": False,
    "source_large_P_star_regime_admitted": False,
    "outer_tail_pressure_compatibility_verified": False,
    "selected_kappa0_global_source_smallness_admitted": False,
    "kappa0_continuation_to_X100_materialized": False,
    "final_interpolation_to_Xi_materialized": False,
    "actual_inner_to_outer_bridge_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
    "cartesian_matched_pressure_gradient_materialized": False,
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


@lru_cache(maxsize=4)
def _legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(order, int) or order < 16:
        raise ValueError("quadrature order must be an integer >=16")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


@dataclass(frozen=True)
class KokunoPA10AutonomousReferenceNsToX100:
    """Propagate the public reference ``n_s`` ODE using the autonomous pressure seam."""

    p1_reference: KokunoPA10ReferenceP1ToX100 = field(
        default_factory=KokunoPA10ReferenceP1ToX100,
        repr=False,
        compare=False,
    )
    pressure: KokunoPA10AutonomousAxisPressureSeedToX100 = field(
        default_factory=KokunoPA10AutonomousAxisPressureSeedToX100,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.p1_reference, KokunoPA10ReferenceP1ToX100):
            raise TypeError("p1_reference must be KokunoPA10ReferenceP1ToX100")
        if not isinstance(self.pressure, KokunoPA10AutonomousAxisPressureSeedToX100):
            raise TypeError(
                "pressure must be KokunoPA10AutonomousAxisPressureSeedToX100"
            )
        p1_profile_sha = self.p1_reference.reference.semantic_sha256
        pressure_profile_sha = self.pressure.pressure_increment.reference.semantic_sha256
        if p1_profile_sha != pressure_profile_sha:
            raise ValueError(
                "p1/reference and pressure/reference profile identities differ; "
                "the axial stress continuation must use one reference profile"
            )
        if abs(self.pressure.P_star - 1.0) > 0.0:
            raise ValueError("autonomous pressure P_star identity drifted")

    @property
    def reference(self):
        return self.p1_reference.reference

    @property
    def X_0(self) -> float:
        return float(self.p1_reference.X_0)

    @property
    def X_b_ref(self) -> float:
        return float(self.p1_reference.X_b_ref)

    @property
    def h(self) -> float:
        return float(self.p1_reference.h)

    @property
    def D(self) -> float:
        return float(self.p1_reference.D)

    @property
    def A(self) -> float:
        return 0.5 + self.h

    @property
    def eta_interval(self) -> tuple[float, float]:
        # The autonomous axis seed is intentionally defined only on [-1,1].
        return (-1.0, 1.0)

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        x_hi = math.nextafter(self.X_b_ref, math.inf)
        if np.any(X_arr < self.X_0) or np.any(X_arr > x_hi):
            raise ValueError(
                f"X must lie in the autonomous n_s propagation domain "
                f"[{self.X_0},{self.X_b_ref}]"
            )
        margin = 2.0 * ETA_FD_STEP
        if np.any(eta_arr < -1.0 + margin) or np.any(eta_arr > 1.0 - margin):
            raise ValueError(
                "eta lies too close to [-1,1] endpoints for the frozen centered eta derivative"
            )
        return np.minimum(X_arr, self.X_b_ref), eta_arr

    def _eta_derivative(self, evaluator, X: np.ndarray, eta: np.ndarray) -> np.ndarray:
        h = ETA_FD_STEP
        return (
            -evaluator(X, eta + 2.0 * h)
            + 8.0 * evaluator(X, eta + h)
            - 8.0 * evaluator(X, eta - h)
            + evaluator(X, eta - 2.0 * h)
        ) / (12.0 * h)

    def source_terms(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate the public axial-source formula with the autonomous pressure."""
        X_arr, eta_arr = self._broadcast(X, eta)
        reference_values = self.reference.values(X_arr, eta_arr)
        reference_derivatives = self.reference.radial_derivatives(X_arr, eta_arr)
        kinematics = self.p1_reference.source_terms(X_arr, eta_arr)

        def U_eval(x: np.ndarray, e: np.ndarray) -> np.ndarray:
            return self.reference.values(x, e)["U_reference"]

        U = reference_values["U_reference"]
        U_X = reference_derivatives["U_reference_X"]
        U_eta = self._eta_derivative(U_eval, X_arr, eta_arr)

        Pi = self.pressure.pressure(X_arr, eta_arr)
        Pi_X = self.pressure.radial_derivative(X_arr, eta_arr)
        Pi_eta = self.pressure.eta_derivative(X_arr, eta_arr)

        d = 1.0 - eta_arr * eta_arr
        W = kinematics["W_reference"]
        H_c = kinematics["H_c_reference"]
        L = kinematics["L"]
        S_n = (
            -W * X_arr * U_X
            - self.A * (1.0 - 2.0 * eta_arr * U) * U
            - H_c * U_eta
            - d * Pi_eta
            + 4.0 * self.A * eta_arr * Pi
            + 2.0 * eta_arr * X_arr * Pi_X
        )
        if np.any(~np.isfinite(S_n)):
            raise RuntimeError("autonomous-pressure S_n produced a non-finite value")
        return {
            "X": X_arr,
            "eta": eta_arr,
            "d": d,
            "L": L,
            "W_reference": W,
            "H_c_reference": H_c,
            "U_reference": U,
            "U_reference_X": U_X,
            "U_reference_eta": U_eta,
            "Pi_reference_autonomous": Pi,
            "Pi_reference_autonomous_X": Pi_X,
            "Pi_reference_autonomous_eta": Pi_eta,
            "S_n_reference_autonomous": S_n,
        }

    def initial_values(self, eta: Any) -> dict[str, np.ndarray]:
        """Return the source stress-free handoff at ``X_0=4/Lambda``.

        This is a value handoff only.  Because the subsequent pressure datum is
        repository-autonomous rather than the source-prepared Appendix-A datum,
        this module does not claim the autonomous ``S_n`` derivative at X_0 is
        the paper's stress-free derivative.
        """
        eta_arr = _finite(eta, "eta")
        margin = 2.0 * ETA_FD_STEP
        if np.any(eta_arr < -1.0 + margin) or np.any(eta_arr > 1.0 - margin):
            raise ValueError("eta outside autonomous n_s numerical-derivative interior")
        X = np.full_like(eta_arr, self.X_0, dtype=float)
        source = self.reference.collar.source_normalized
        U_X = source.derivatives(X, eta_arr)["U_0_X"]
        L = 1.0 - 2.0 * self.h * eta_arr * eta_arr
        n_s = -2.0 * U_X
        N_s = L * n_s
        return {
            "X": X,
            "eta": eta_arr,
            "L": L,
            "U_0_X": U_X,
            "n_s_initial": n_s,
            "N_s_initial": N_s,
        }

    def values(
        self, X: Any, eta: Any, *, order: int = NS_QUADRATURE_ORDER
    ) -> dict[str, np.ndarray]:
        """Return ``N_s,n_s,p_2`` on ``X_0<=X<=100``.

        From ``D_X N_s+N_s=S_n`` one has

            X N_s(X) = X_0 N_s(X_0) + integral_X0^X S_n(s) ds.
        """
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        initial = self.initial_values(e)

        nodes, weights = _legendre_rule(order)
        half = 0.5 * (x - self.X_0)
        mid = 0.5 * (x + self.X_0)
        X_nodes = mid[:, None] + half[:, None] * nodes[None, :]
        eta_nodes = np.broadcast_to(e[:, None], X_nodes.shape)
        S_nodes = self.source_terms(X_nodes, eta_nodes)["S_n_reference_autonomous"]
        integral = half * np.sum(weights[None, :] * S_nodes, axis=1)

        N_s = (self.X_0 * initial["N_s_initial"] + integral) / x
        L = 1.0 - 2.0 * self.h * e * e
        n_s = N_s / L
        E = self.reference.values(x, e)["E_reference"]
        if np.any(E <= 0.0):
            raise RuntimeError("reference E must remain positive for p2=X n_s/E")
        p2 = x * n_s / E
        if np.any(~np.isfinite(N_s)) or np.any(~np.isfinite(n_s)) or np.any(~np.isfinite(p2)):
            raise RuntimeError("autonomous reference stress propagation became non-finite")
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "N_s_reference_autonomous": N_s.reshape(shape),
            "n_s_reference_autonomous": n_s.reshape(shape),
            "p2_reference_autonomous": p2.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return source-ODE radial derivatives of ``N_s`` and ``n_s``."""
        X_arr, eta_arr = self._broadcast(X, eta)
        values = self.values(X_arr, eta_arr)
        terms = self.source_terms(X_arr, eta_arr)
        N_s = values["N_s_reference_autonomous"]
        L = terms["L"]
        N_s_X = (terms["S_n_reference_autonomous"] - N_s) / X_arr
        n_s_X = N_s_X / L
        return {
            "N_s_reference_autonomous_X": N_s_X,
            "n_s_reference_autonomous_X": n_s_X,
        }

    def handoff_at_X100(self, eta: Any) -> dict[str, np.ndarray]:
        eta_arr = _finite(eta, "eta")
        X = np.full_like(eta_arr, self.X_b_ref, dtype=float)
        return {
            **self.values(X, eta_arr),
            **self.radial_derivatives(X, eta_arr),
            **self.source_terms(X, eta_arr),
            "p1_reference": self.p1_reference.values(X, eta_arr)["p1_reference"],
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "p1_reference_semantic_sha256": self.p1_reference.semantic_sha256,
            "autonomous_pressure_semantic_sha256": self.pressure.semantic_sha256,
            "eta_fd_step": ETA_FD_STEP,
            "Ns_quadrature_order": NS_QUADRATURE_ORDER,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10AutonomousReferenceNsToX100":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError(
                "autonomous-reference-ns configuration is frozen; serialized data differ"
            )
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
    ) -> "KokunoPA10AutonomousReferenceNsToX100":
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
        initial = self.initial_values(eta_probe)
        handoff = self.handoff_at_X100(eta_probe)
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
            "semantic_sha256": self.semantic_sha256,
            "initial_probe": {
                "eta": eta_probe.tolist(),
                "n_s_initial": initial["n_s_initial"].tolist(),
                "N_s_initial": initial["N_s_initial"].tolist(),
            },
            "X100_handoff_probe": {
                "eta": eta_probe.tolist(),
                "p1_reference": handoff["p1_reference"].tolist(),
                "n_s_reference_autonomous": handoff["n_s_reference_autonomous"].tolist(),
                "p2_reference_autonomous": handoff["p2_reference_autonomous"].tolist(),
                "S_n_reference_autonomous": handoff["S_n_reference_autonomous"].tolist(),
                "n_s_reference_autonomous_X": handoff[
                    "n_s_reference_autonomous_X"
                ].tolist(),
            },
            "truth_boundary": self.truth_boundary,
        }

    def save_report(self, path: str | Path) -> dict[str, Any]:
        payload = self.report()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload


def _main() -> None:
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--config-output", type=Path)
    args = parser.parse_args()
    stress = KokunoPA10AutonomousReferenceNsToX100()
    payload = stress.save_report(args.output)
    if args.config_output is not None:
        stress.save_configuration(args.config_output)
    print("semantic_sha256=", payload["semantic_sha256"])
    print("initial_probe=", payload["initial_probe"])
    print("X100_handoff_probe=", payload["X100_handoff_probe"])


if __name__ == "__main__":
    _main()
