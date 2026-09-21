"""Executable fixed-kappa continuation of the Kokuno PA.10 leading profile to X=100.

Pinned public provenance: ``KokunoYumeto/yang-mills-interacting-workbench`` at
``143f6773feb424ad9ed3a8d116653200f20346b7``, corrected 2026-09-09 reader.

After the first stress-activation collar the corrected reconstruction keeps
``kappa=kappa_0`` and continues the actual profile with the reference stress
primitives

    D_X U = -kappa_0 X n_{s,r}/2,
    partial_y log(phi) = -kappa_0 p_{1,r}/2,
    y = log(X/X_0),                    D_X = X partial_X,

so, because the normalization C is constant,

    partial_X log(F) = -kappa_0 p_{1,r}/(2X),
    partial_X U      = -kappa_0 n_{s,r}/2.

This module starts exactly from the already materialized Agent-1 activation
endpoint ``X_1=X_0 exp(t_1)`` and propagates those two equations through the
public reference cutoff ``X_b^ref=100``.  It consumes the current executable
reference ``p_{1,r}`` and the explicitly repository-autonomous-pressure
``n_{s,r}`` from #918.

Consequently this is a candidate-side executable continuation, not the
source-prepared Appendix-A continuation: the numerical ``Pi_0`` underlying the
current ``n_{s,r}`` is autonomous, and the frozen ``kappa_0=0.1`` has not been
admitted by the source's later global smallness/cone proof.  No residual,
forcing, held-out sample, image fit, or PDE threshold is used to choose any
parameter here.
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
from .kokuno_pa10_stress_activation_collar import (
    SELECTED_KAPPA0,
    KokunoPA10StressActivationCollar,
)


SCHEMA = "kokuno-pa10-fixed-kappa-x100-v1"
CONTINUATION_QUADRATURE_ORDER = 48

_SOURCE_FORMULAS = {
    "continuation_range": "X_1<=X<=X_b^ref=100 after the activation collar",
    "fixed_kappa": "kappa=kappa_0",
    "axial_profile": "D_X U=-kappa_0 X n_{s,r}/2; D_X=X partial_X",
    "swirl_profile": "partial_y log(phi)=-kappa_0 p_{1,r}/2; y=log(X/X_0)",
    "normalized_swirl": "C constant => partial_X log(F)=-kappa_0 p_{1,r}/(2X)",
    "azimuthal_profile": "E=sqrt(2X) F",
    "actual_shear": "s_shear=kappa_0 (p_{1,r}, p_{2,r} E_r/E)",
}

_NUMERICAL_REALIZATION = {
    "kappa_0": (
        "frozen repository-autonomous kappa_0=0.1 inherited from the activation collar; "
        "not residual-tuned and not source-smallness-admitted"
    ),
    "reference_n_s": (
        "Agent-1 autonomous-pressure reference n_s from #918; not source-prepared Appendix-A n_s"
    ),
    "radial_integration": (
        f"fixed Gauss-Legendre order {CONTINUATION_QUADRATURE_ORDER} on [X_1,X]"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_fixed_kappa_continuation_formula_executable": True,
    "activation_endpoint_consumed": True,
    "autonomous_pressure_reference_stress_pair_consumed": True,
    "candidate_side_fixed_kappa_continuation_to_X100_executable": True,
    "candidate_side_X100_profile_handoff_executable": True,
    "selected_kappa0_is_repository_autonomous": True,
    "selected_kappa0_chosen_from_ns_residual": False,
    "selected_kappa0_global_source_smallness_admitted": False,
    "source_prepared_appendixA_Pi0_materialized": False,
    "source_prepared_reference_nsr_materialized": False,
    "source_prepared_full_reference_stress_pair_materialized": False,
    "source_exact_fixed_kappa_continuation_to_X100_materialized": False,
    "final_interpolation_to_Xi_materialized": False,
    "actual_G_i_at_Xi_materialized": False,
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
class KokunoPA10FixedKappaContinuationToX100:
    """Continue the actual PA.10 profile at frozen ``kappa_0`` through X=100."""

    activation: KokunoPA10StressActivationCollar = field(
        default_factory=KokunoPA10StressActivationCollar,
        repr=False,
        compare=False,
    )
    reference_stress: KokunoPA10AutonomousReferenceNsToX100 = field(
        default_factory=KokunoPA10AutonomousReferenceNsToX100,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.activation, KokunoPA10StressActivationCollar):
            raise TypeError("activation must be KokunoPA10StressActivationCollar")
        if not isinstance(self.reference_stress, KokunoPA10AutonomousReferenceNsToX100):
            raise TypeError(
                "reference_stress must be KokunoPA10AutonomousReferenceNsToX100"
            )
        # #918's reference plateau must descend from exactly the same #880 collar
        # that supplies the actual activation endpoint.
        stress_collar_sha = self.reference_stress.reference.collar.semantic_sha256
        if stress_collar_sha != self.activation.reference.semantic_sha256:
            raise ValueError(
                "activation and reference-stress backbones differ; fixed-kappa continuation "
                "requires one shared reference collar identity"
            )
        if abs(SELECTED_KAPPA0 - 0.1) > 0.0:
            raise ValueError("registered autonomous kappa_0 identity drifted")
        if not self.X_1 < self.X_b_ref:
            raise ValueError("fixed-kappa continuation requires X_1 < X_b^ref")

    @property
    def X_1(self) -> float:
        return float(self.activation.X_1)

    @property
    def X_b_ref(self) -> float:
        return float(self.reference_stress.X_b_ref)

    @property
    def kappa_0(self) -> float:
        return float(SELECTED_KAPPA0)

    @property
    def eta_interval(self) -> tuple[float, float]:
        # The activation collar is the narrower executable source interval and
        # therefore controls the admissible continuation interval.
        return self.activation.eta_interval

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        x_lo = math.nextafter(self.X_1, -math.inf)
        x_hi = math.nextafter(self.X_b_ref, math.inf)
        if np.any((X_arr < x_lo) | (X_arr > x_hi)):
            raise ValueError(
                f"X must lie in the fixed-kappa continuation domain "
                f"[{self.X_1},{self.X_b_ref}]"
            )
        eta_lo, eta_hi = self.eta_interval
        if np.any((eta_arr < eta_lo) | (eta_arr > eta_hi)):
            raise ValueError(
                f"eta must lie in the executable activation interval [{eta_lo},{eta_hi}]"
            )
        return np.clip(X_arr, self.X_1, self.X_b_ref), eta_arr

    def endpoint_values(self, eta: Any) -> dict[str, np.ndarray]:
        """Return the exact actual activation handoff at ``X_1``."""
        eta_arr = _finite(eta, "eta")
        eta_lo, eta_hi = self.eta_interval
        if np.any((eta_arr < eta_lo) | (eta_arr > eta_hi)):
            raise ValueError("eta outside activation endpoint interval")
        X = np.full_like(eta_arr, self.X_1, dtype=float)
        values = self.activation.values(X, eta_arr)
        derivatives = self.activation.radial_derivatives(X, eta_arr)
        return {**values, **derivatives}

    def reference_stress_data(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return candidate-side reference p1, n_s, p2 and reference E on [X1,100]."""
        X_arr, eta_arr = self._broadcast(X, eta)
        ns = self.reference_stress.values(X_arr, eta_arr)
        p1 = self.reference_stress.p1_reference.values(X_arr, eta_arr)["p1_reference"]
        E_r = self.reference_stress.reference.values(X_arr, eta_arr)["E_reference"]
        return {
            "p1_reference": p1,
            "n_s_reference_autonomous": ns["n_s_reference_autonomous"],
            "p2_reference_autonomous": ns["p2_reference_autonomous"],
            "E_reference": E_r,
        }

    def _continuation_integrals(
        self, X: np.ndarray, eta: np.ndarray, *, order: int = CONTINUATION_QUADRATURE_ORDER
    ) -> tuple[np.ndarray, np.ndarray]:
        """Integrate the two displayed fixed-kappa profile ODEs from X1 to X."""
        if X.ndim != 1 or eta.ndim != 1 or X.shape != eta.shape:
            raise ValueError("continuation integration expects matching flat vectors")
        if X.size == 0:
            return np.empty(0), np.empty(0)
        nodes, weights = _legendre_rule(order)
        half = 0.5 * (X - self.X_1)
        mid = 0.5 * (X + self.X_1)
        X_nodes = mid[:, None] + half[:, None] * nodes[None, :]
        eta_nodes = np.broadcast_to(eta[:, None], X_nodes.shape)
        stress = self.reference_stress_data(X_nodes, eta_nodes)
        d_log_F_dX = -0.5 * self.kappa_0 * stress["p1_reference"] / X_nodes
        d_U_dX = -0.5 * self.kappa_0 * stress["n_s_reference_autonomous"]
        delta_log_F = half * np.sum(weights[None, :] * d_log_F_dX, axis=1)
        delta_U = half * np.sum(weights[None, :] * d_U_dX, axis=1)
        return delta_log_F, delta_U

    def values(
        self, X: Any, eta: Any, *, order: int = CONTINUATION_QUADRATURE_ORDER
    ) -> dict[str, np.ndarray]:
        """Return vectorized actual ``F,U,E`` on ``X_1<=X<=100``."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        endpoint = self.endpoint_values(e)
        delta_log_F, delta_U = self._continuation_integrals(x, e, order=order)
        F = endpoint["F_activation"] * np.exp(delta_log_F)
        U = endpoint["U_activation"] + delta_U
        E = np.sqrt(2.0 * x) * F
        if np.any(F <= 0.0) or np.any(E <= 0.0):
            raise RuntimeError("fixed-kappa continuation lost positive F/E")
        if np.any(~np.isfinite(F)) or np.any(~np.isfinite(U)) or np.any(~np.isfinite(E)):
            raise RuntimeError("fixed-kappa continuation produced non-finite profile values")
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "F_fixed_kappa": F.reshape(shape),
            "U_fixed_kappa": U.reshape(shape),
            "E_fixed_kappa": E.reshape(shape),
            "delta_log_F_from_activation_endpoint": delta_log_F.reshape(shape),
            "delta_U_from_activation_endpoint": delta_U.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return analytic first radial derivatives from the public continuation ODE."""
        X_arr, eta_arr = self._broadcast(X, eta)
        values = self.values(X_arr, eta_arr)
        stress = self.reference_stress_data(X_arr, eta_arr)
        F = values["F_fixed_kappa"]
        F_X = -0.5 * self.kappa_0 * stress["p1_reference"] * F / X_arr
        U_X = -0.5 * self.kappa_0 * stress["n_s_reference_autonomous"]
        E_X = F / np.sqrt(2.0 * X_arr) + np.sqrt(2.0 * X_arr) * F_X
        return {
            "F_fixed_kappa_X": F_X,
            "U_fixed_kappa_X": U_X,
            "E_fixed_kappa_X": E_X,
        }

    def shear_schedule(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return the displayed fixed-kappa actual shear schedule."""
        X_arr, eta_arr = self._broadcast(X, eta)
        actual = self.values(X_arr, eta_arr)
        stress = self.reference_stress_data(X_arr, eta_arr)
        a = self.kappa_0 * stress["p1_reference"]
        second = (
            self.kappa_0
            * stress["p2_reference_autonomous"]
            * stress["E_reference"]
            / actual["E_fixed_kappa"]
        )
        return {
            "kappa": np.full_like(X_arr, self.kappa_0, dtype=float),
            "a": a,
            "shear_first": a,
            "shear_second": second,
        }

    def handoff_at_X100(self, eta: Any) -> dict[str, np.ndarray]:
        eta_arr = _finite(eta, "eta")
        X = np.full_like(eta_arr, self.X_b_ref, dtype=float)
        return {
            **self.values(X, eta_arr),
            **self.radial_derivatives(X, eta_arr),
            **self.shear_schedule(X, eta_arr),
            **self.reference_stress_data(X, eta_arr),
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "activation_semantic_sha256": self.activation.semantic_sha256,
            "reference_stress_semantic_sha256": self.reference_stress.semantic_sha256,
            "selected_kappa0": self.kappa_0,
            "continuation_quadrature_order": CONTINUATION_QUADRATURE_ORDER,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10FixedKappaContinuationToX100":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError(
                "fixed-kappa continuation configuration is frozen; serialized data differ"
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
    ) -> "KokunoPA10FixedKappaContinuationToX100":
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
        endpoint = self.endpoint_values(eta_probe)
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
            "geometry": {
                "X_1": self.X_1,
                "X_b_ref": self.X_b_ref,
                "selected_kappa0": self.kappa_0,
            },
            "activation_endpoint_probe": {
                "eta": eta_probe.tolist(),
                "F_activation": endpoint["F_activation"].tolist(),
                "U_activation": endpoint["U_activation"].tolist(),
            },
            "X100_handoff_probe": {
                "eta": eta_probe.tolist(),
                "F_fixed_kappa": handoff["F_fixed_kappa"].tolist(),
                "U_fixed_kappa": handoff["U_fixed_kappa"].tolist(),
                "E_fixed_kappa": handoff["E_fixed_kappa"].tolist(),
                "p1_reference": handoff["p1_reference"].tolist(),
                "n_s_reference_autonomous": handoff[
                    "n_s_reference_autonomous"
                ].tolist(),
                "shear_first": handoff["shear_first"].tolist(),
                "shear_second": handoff["shear_second"].tolist(),
            },
            "semantic_sha256": self.semantic_sha256,
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
    continuation = KokunoPA10FixedKappaContinuationToX100()
    payload = continuation.save_report(args.output)
    if args.config_output is not None:
        continuation.save_configuration(args.config_output)
    print("semantic_sha256=", payload["semantic_sha256"])
    print("geometry=", payload["geometry"])
    print("X100_handoff_probe=", payload["X100_handoff_probe"])


if __name__ == "__main__":
    _main()
