"""Executable public reference angular-stress primitive ``p_{1,r}`` to X=100.

Pinned provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.

For the reference leading profile the corrected reconstruction defines

    H = 2 X F,
    ell = X d_X log H,
    W = 1 - (2 D eta M + d M_eta)/X,
    H_c = D eta + d U,
    S_q = -W ell - h(1-2 eta U) - H_c (log E)_eta,

and the regular angular-stress primitive by

    D_X p_{1,r} = X S_{q,r}/L - ell_r p_{1,r},
    D_X = X d_X.

Since ``D_X log H=ell``, the source ODE has the integrating-factor form

    p1(X) = [H(X0)p1(X0) + int_{X0}^X H(s) S_q(s)/L ds] / H(X),

with the stress-free natural initial datum

    p1(X0) = -2 D_X log F_r(X0).

This module materializes exactly that *p1* seam on top of the public reference
profile already continued through ``X_b^ref=100``.  The source needs eta
parameter derivatives inside ``S_q``.  The current reference-collar API does
not yet expose analytic eta jets, so this increment uses one frozen fourth-
order centered eta-difference realization (repository numerical realization,
not hidden source data) while keeping the radial p1 ODE/integrating factor
source-exact.  ``n_{s,r}``, the fixed-kappa actual continuation, pressure,
forcing, global velocity and PDE validation remain deliberately out of scope.
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

from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa10_reference_plateau_x100 import KokunoPA10ReferencePlateauToX100


SCHEMA = "kokuno-pa10-reference-p1-x100-v1"
ETA_FD_STEP = 2.0e-5
M_QUADRATURE_ORDER = 40
P1_QUADRATURE_ORDER = 48

_SOURCE_FORMULAS = {
    "reference_kinematics": (
        "H=2XF; ell=X partial_X log H; W=1-(2D eta M+d M_eta)/X; "
        "H_c=D eta+dU"
    ),
    "angular_source": "S_q=-W ell-h(1-2eta U)-H_c (log E)_eta",
    "p1_ode": "D_X p_{1,r}=X S_{q,r}/L-ell_r p_{1,r}; D_X=X partial_X",
    "p1_integrating_factor": (
        "p1(X)=[H(X0)p1(X0)+int_X0^X H S_q/L ds]/H(X)"
    ),
    "stress_free_initial": "p1(X0)=-2 D_X log F_r(X0)",
}

_NUMERICAL_REALIZATION = {
    "eta_derivative": "fixed fourth-order centered difference",
    "eta_step": ETA_FD_STEP,
    "M_quadrature": "Gauss-Legendre",
    "M_quadrature_order": M_QUADRATURE_ORDER,
    "p1_quadrature": "Gauss-Legendre",
    "p1_quadrature_order": P1_QUADRATURE_ORDER,
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_reference_p1_ode_executable_to_X100": True,
    "reference_p1_vectorized_to_X100": True,
    "reference_p1_radial_derivative_executable": True,
    "reference_Sq_executable": True,
    "reference_eta_jets_analytic": False,
    "reference_eta_jets_fixed_fd_realization": True,
    "reference_nsr_materialized": False,
    "reference_full_stress_pair_p1r_nsr_materialized": False,
    "selected_kappa0_global_source_smallness_admitted": False,
    "kappa0_continuation_to_X100_materialized": False,
    "final_interpolation_to_Xi_materialized": False,
    "actual_inner_to_outer_bridge_materialized": False,
    "outer_global_leading_velocity_materialized": False,
    "matched_global_pressure_materialized": False,
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


@lru_cache(maxsize=8)
def _legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(order, int) or order < 16:
        raise ValueError("quadrature order must be an integer >=16")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


@dataclass(frozen=True)
class KokunoPA10ReferenceP1ToX100:
    """Materialize the public reference ``p_{1,r}`` primitive to ``X=100``."""

    reference: KokunoPA10ReferencePlateauToX100 = field(
        default_factory=KokunoPA10ReferencePlateauToX100,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.reference, KokunoPA10ReferencePlateauToX100):
            raise TypeError("reference must be KokunoPA10ReferencePlateauToX100")
        if not self.X_0 > 0.0:
            raise ValueError("source X_0 must be positive")

    @property
    def X_0(self) -> float:
        return float(self.reference.X_0)

    @property
    def X_b_ref(self) -> float:
        return float(self.reference.X_b_ref)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.reference.eta_interval

    @property
    def axis_domain(self):
        return self.reference.collar.source_normalized.axis_domain

    @property
    def h(self) -> float:
        return float(self.axis_domain.h)

    @property
    def D(self) -> float:
        return float(self.axis_domain.D)

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        x_hi = math.nextafter(self.X_b_ref, math.inf)
        if np.any(X_arr < self.X_0) or np.any(X_arr > x_hi):
            raise ValueError(f"X must lie in the p1 propagation domain [{self.X_0},{self.X_b_ref}]")
        lo, hi = self.eta_interval
        margin = 2.0 * ETA_FD_STEP
        if np.any(eta_arr < lo + margin) or np.any(eta_arr > hi - margin):
            raise ValueError(
                "eta lies too close to the executable reference boundary for the frozen centered eta derivative"
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

    def _M_reference_flat(
        self, X: np.ndarray, eta: np.ndarray, *, order: int = M_QUADRATURE_ORDER
    ) -> np.ndarray:
        """Return M=int_0^X U_r dx, reusing the analytic natural M at X0."""
        if X.ndim != 1 or eta.ndim != 1 or X.shape != eta.shape:
            raise ValueError("M integration expects matching flat vectors")
        if X.size == 0:
            return np.empty(0, dtype=float)
        source = self.reference.collar.source_normalized
        x0 = np.full_like(X, self.X_0)
        initial = source.values(x0, eta)["M_0"]
        nodes, weights = _legendre_rule(order)
        half = 0.5 * (X - self.X_0)
        mid = 0.5 * (X + self.X_0)
        X_nodes = mid[:, None] + half[:, None] * nodes[None, :]
        eta_nodes = np.broadcast_to(eta[:, None], X_nodes.shape)
        U_nodes = self.reference.values(X_nodes, eta_nodes)["U_reference"]
        return initial + half * np.sum(weights[None, :] * U_nodes, axis=1)

    def M_reference(self, X: Any, eta: Any) -> np.ndarray:
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        result = self._M_reference_flat(X_arr.reshape(-1), eta_arr.reshape(-1))
        return result.reshape(shape)

    def source_terms(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate the public reference kinematics and angular source ``S_q``."""
        X_arr, eta_arr = self._broadcast(X, eta)
        values = self.reference.values(X_arr, eta_arr)
        slopes = self.reference.radial_log_slopes(X_arr, eta_arr)
        M = self.M_reference(X_arr, eta_arr)

        def M_eval(x: np.ndarray, e: np.ndarray) -> np.ndarray:
            shape = np.broadcast(x, e).shape
            xb, eb = np.broadcast_arrays(x, e)
            return self._M_reference_flat(xb.reshape(-1), eb.reshape(-1)).reshape(shape)

        M_eta = self._eta_derivative(M_eval, X_arr, eta_arr)

        def log_E_eval(x: np.ndarray, e: np.ndarray) -> np.ndarray:
            E = self.reference.values(x, e)["E_reference"]
            if np.any(E <= 0.0):
                raise RuntimeError("reference E must be positive for log derivative")
            return np.log(E)

        log_E_eta = self._eta_derivative(log_E_eval, X_arr, eta_arr)
        U = values["U_reference"]
        ell = slopes["l_reference"]
        d = 1.0 - eta_arr * eta_arr
        L = 1.0 - 2.0 * self.h * eta_arr * eta_arr
        if np.any(L <= 0.0):
            raise RuntimeError("source L lost positivity on the executable eta interval")
        W = 1.0 - (2.0 * self.D * eta_arr * M + d * M_eta) / X_arr
        H_c = self.D * eta_arr + d * U
        S_q = -W * ell - self.h * (1.0 - 2.0 * eta_arr * U) - H_c * log_E_eta
        H = 2.0 * X_arr * values["F_reference"]
        return {
            "M_reference": M,
            "M_reference_eta": M_eta,
            "log_E_reference_eta": log_E_eta,
            "W_reference": W,
            "H_c_reference": H_c,
            "H_reference": H,
            "ell_reference": ell,
            "L": L,
            "S_q_reference": S_q,
        }

    def initial_p1(self, eta: Any) -> np.ndarray:
        eta_arr = _finite(eta, "eta")
        lo, hi = self.eta_interval
        margin = 2.0 * ETA_FD_STEP
        if np.any(eta_arr < lo + margin) or np.any(eta_arr > hi - margin):
            raise ValueError("eta outside p1 numerical-derivative interior")
        X = np.full_like(eta_arr, self.X_0, dtype=float)
        source = self.reference.collar.source_normalized
        values = source.values(X, eta_arr)
        derivatives = source.derivatives(X, eta_arr)
        return -2.0 * X * derivatives["F_0_X"] / values["F_0"]

    def values(
        self, X: Any, eta: Any, *, order: int = P1_QUADRATURE_ORDER
    ) -> dict[str, np.ndarray]:
        """Return vectorized ``p_{1,r}`` using the public integrating factor."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        p0 = self.initial_p1(e)
        H0 = 2.0 * self.X_0 * self.reference.values(
            np.full_like(e, self.X_0), e
        )["F_reference"]

        nodes, weights = _legendre_rule(order)
        half = 0.5 * (x - self.X_0)
        mid = 0.5 * (x + self.X_0)
        X_nodes = mid[:, None] + half[:, None] * nodes[None, :]
        eta_nodes = np.broadcast_to(e[:, None], X_nodes.shape)
        terms = self.source_terms(X_nodes, eta_nodes)
        integrand = terms["H_reference"] * terms["S_q_reference"] / terms["L"]
        integral = half * np.sum(weights[None, :] * integrand, axis=1)
        final = self.reference.values(x, e)
        H = 2.0 * x * final["F_reference"]
        p1 = (H0 * p0 + integral) / H
        if np.any(~np.isfinite(p1)):
            raise RuntimeError("reference p1 propagation produced a non-finite value")
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "p1_reference": p1.reshape(shape),
        }

    def radial_derivative(self, X: Any, eta: Any) -> np.ndarray:
        """Return ``partial_X p1`` directly from the public first-order ODE."""
        X_arr, eta_arr = self._broadcast(X, eta)
        p1 = self.values(X_arr, eta_arr)["p1_reference"]
        terms = self.source_terms(X_arr, eta_arr)
        return terms["S_q_reference"] / terms["L"] - terms["ell_reference"] * p1 / X_arr

    def handoff_at_X100(self, eta: Any) -> dict[str, np.ndarray]:
        eta_arr = _finite(eta, "eta")
        X = np.full_like(eta_arr, self.X_b_ref, dtype=float)
        values = self.values(X, eta_arr)
        return {
            **values,
            "p1_reference_X": self.radial_derivative(X, eta_arr),
            **self.source_terms(X, eta_arr),
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_reference_semantic_sha256": self.reference.semantic_sha256,
            "eta_fd_step": ETA_FD_STEP,
            "M_quadrature_order": M_QUADRATURE_ORDER,
            "p1_quadrature_order": P1_QUADRATURE_ORDER,
        }

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "KokunoPA10ReferenceP1ToX100":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError("reference-p1 configuration is frozen; serialized data differ")
        return cls()

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA10ReferenceP1ToX100":
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
        eta = np.asarray([-0.75, -0.25, 0.0, 0.25, 0.75])
        x0 = self.values(np.full_like(eta, self.X_0), eta)["p1_reference"]
        x100 = self.handoff_at_X100(eta)
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
            "geometry": {"X_0": self.X_0, "X_b_ref": self.X_b_ref},
            "probe": {
                "eta": eta.tolist(),
                "p1_X0": x0.tolist(),
                "p1_X100": x100["p1_reference"].tolist(),
                "p1_X100_X": x100["p1_reference_X"].tolist(),
                "S_q_X100": x100["S_q_reference"].tolist(),
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
    primitive = KokunoPA10ReferenceP1ToX100()
    payload = primitive.save_report(args.output)
    if args.config_output is not None:
        primitive.save_configuration(args.config_output)
    print("semantic_sha256=", payload["semantic_sha256"])
    print("geometry=", payload["geometry"])
    print("probe=", payload["probe"])


if __name__ == "__main__":
    _main()
