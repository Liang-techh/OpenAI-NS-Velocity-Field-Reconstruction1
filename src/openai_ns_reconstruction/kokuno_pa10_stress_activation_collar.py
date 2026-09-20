"""Executable first stress-activation collar for the Kokuno leading profile.

Pinned public provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.

On the source activation interval y=log(X/X_0) in [0,t_1], the reference
profile is still the natural PA.10 profile.  The corrected reconstruction gives

    e_a = (1-kappa_0) sigma(y/t_1),   kappa = 1-e_a,
    a = kappa p_{1,r},
    D_X U = -kappa X n_{s,r}/2,
    partial_y log(phi) = -kappa p_{1,r}/2,

with the stress-free reference identities

    p_{1,r} = -2 D_X log(phi_r),
    n_{s,r} = -2 partial_X U_r,
    p_{2,r} = X n_{s,r}/E_r.

Because C is constant, the same logarithmic equation holds for F=phi/C.  This
module therefore turns the public activation prescription into executable
F/U/E values, analytic radial derivatives, and the source shear schedule on
this first collar.  It consumes the exact #880 reference-collar implementation
and does not extrapolate beyond X_0 exp(t_1).

The source chooses kappa_0 only after proof-dependent bounds such as
kappa_0 V_max < 1/4 have been fixed.  Those global reference bounds are not yet
materialized in this repository.  We consequently freeze kappa_0=0.1 as an
explicit repository-autonomous realization, selected before residual work and
satisfying only the source's displayed 0<kappa_0<1/2 requirement.  It is not
claimed to be the paper/OpenAI hidden value or to satisfy the later global
smallness/cone proof.
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
from .kokuno_pa10_reference_continuation_collar import (
    QUADRATURE_ORDER,
    SELECTED_T1,
    KokunoPA10ReferenceContinuationCollar,
    smooth_step,
)


SCHEMA = "kokuno-pa10-stress-activation-collar-v1"
SELECTED_KAPPA0 = 0.1

_SOURCE_FORMULAS = {
    "activation_coordinate": "y=log(X/X_0), 0<y<t_1",
    "activation_flat_factor": "e_a=(1-kappa_0) sigma(y/t_1)",
    "activation_kappa": "kappa=1-e_a",
    "angular_shear": "a=kappa p_{1,r}",
    "axial_profile_slope": "D_X U=-kappa X n_{s,r}/2",
    "swirl_profile_slope": "partial_y log phi=-kappa p_{1,r}/2",
    "reference_p1": "p_{1,r}=-2 D_X log(phi_r)",
    "reference_ns": "n_{s,r}=-2 partial_X U_r",
    "reference_p2": "p_{2,r}=X n_{s,r}/E_r",
    "actual_shear": "s_shear=kappa (p_{1,r},p_{2,r} E_r/E)",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_stress_activation_formula_executable": True,
    "stress_activation_materialized": True,
    "stress_activation_profile_values_executable": True,
    "stress_activation_analytic_radial_derivatives_executable": True,
    "stress_activation_shear_schedule_executable": True,
    "selected_kappa0_is_repository_autonomous": True,
    "selected_kappa0_chosen_from_ns_residual": False,
    "selected_kappa0_displayed_range_satisfied": True,
    "selected_kappa0_global_source_smallness_admitted": False,
    "stress_activation_cone_admissibility_independently_certified": False,
    "kappa0_continuation_to_X100_materialized": False,
    "final_interpolation_to_Xi_materialized": False,
    "actual_G_i_at_Xi_materialized": False,
    "actual_upstream_five_moment_discrepancy_materialized": False,
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


@lru_cache(maxsize=4)
def _legendre_rule(order: int) -> tuple[np.ndarray, np.ndarray]:
    if not isinstance(order, int) or order < 16:
        raise ValueError("quadrature order must be an integer >=16")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    nodes.setflags(write=False)
    weights.setflags(write=False)
    return nodes, weights


@dataclass(frozen=True)
class KokunoPA10StressActivationCollar:
    """Source-governed first activation collar stacked on the reference profile."""

    reference: KokunoPA10ReferenceContinuationCollar = field(
        default_factory=KokunoPA10ReferenceContinuationCollar,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.reference, KokunoPA10ReferenceContinuationCollar):
            raise TypeError("reference must be KokunoPA10ReferenceContinuationCollar")
        if not (0.0 < SELECTED_KAPPA0 < 0.5):
            raise ValueError("selected kappa_0 violates the displayed source range")

    @property
    def X_0(self) -> float:
        return float(self.reference.X_0)

    @property
    def X_1(self) -> float:
        return float(self.reference.X_1)

    @property
    def Lambda(self) -> float:
        return float(self.reference.Lambda)

    @property
    def C(self) -> float:
        return float(self.reference.C)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return self.reference.eta_interval

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        # Allow tiny binary64 endpoint roundoff only through nextafter bounds.
        x_lo = math.nextafter(self.X_0, -math.inf)
        x_hi = math.nextafter(self.X_1, math.inf)
        if np.any((X_arr < x_lo) | (X_arr > x_hi)):
            raise ValueError(
                f"X must lie in the first stress-activation collar [{self.X_0},{self.X_1}]"
            )
        eta_lo, eta_hi = self.eta_interval
        if np.any((eta_arr < eta_lo) | (eta_arr > eta_hi)):
            raise ValueError(
                f"eta must lie in the executable natural interval [{eta_lo},{eta_hi}]"
            )
        return np.clip(X_arr, self.X_0, self.X_1), eta_arr

    def reference_stress_data(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return the exact stress-free reference p1,n_s,p2,v_r on 0<=y<=t1."""
        X_arr, eta_arr = self._broadcast(X, eta)
        natural = self.reference.source_normalized.values(X_arr, eta_arr)
        natural_d = self.reference.source_normalized.derivatives(X_arr, eta_arr)
        Phi = natural["Phi_0"]
        E = natural["E_0"]
        if np.any(Phi <= 0.0) or np.any(E <= 0.0):
            raise RuntimeError("reference Phi and E must remain positive on activation collar")
        p1 = -2.0 * X_arr * natural_d["Phi_0_X"] / Phi
        ns = -2.0 * natural_d["U_0_X"]
        if np.any(p1 <= 0.0):
            raise RuntimeError("reference p_{1,r} lost its source-positive sign")
        p2 = X_arr * ns / E
        vr = p1 + p2 * p2 / p1
        return {
            "p1_reference": p1,
            "ns_reference": ns,
            "p2_reference": p2,
            "v_reference": vr,
            "E_reference": E,
        }

    def activation_schedule(self, X: Any, eta: Any = 0.0) -> dict[str, np.ndarray]:
        X_arr, eta_arr = self._broadcast(X, eta)
        y = np.log(X_arr / self.X_0)
        s = np.clip(y / SELECTED_T1, 0.0, 1.0)
        sigma = smooth_step(s)
        e_a = (1.0 - SELECTED_KAPPA0) * sigma
        kappa = 1.0 - e_a
        # Enforce exact source endpoint values despite log roundoff.
        e_a = np.where(X_arr <= self.X_0, 0.0, e_a)
        kappa = np.where(X_arr <= self.X_0, 1.0, kappa)
        e_a = np.where(X_arr >= self.X_1, 1.0 - SELECTED_KAPPA0, e_a)
        kappa = np.where(X_arr >= self.X_1, SELECTED_KAPPA0, kappa)
        return {
            "X": X_arr,
            "eta": eta_arr,
            "y": y,
            "sigma": sigma,
            "e_a": e_a,
            "kappa": kappa,
        }

    def _activation_integrals(
        self, X: np.ndarray, eta: np.ndarray, *, order: int = QUADRATURE_ORDER
    ) -> tuple[np.ndarray, np.ndarray]:
        """Integrate the public exact-difference equations from y=0 to y(X)."""
        if X.ndim != 1 or eta.ndim != 1 or X.shape != eta.shape:
            raise ValueError("activation integration expects matching flat vectors")
        if X.size == 0:
            return np.empty(0), np.empty(0)
        nodes, weights = _legendre_rule(order)
        y = np.log(X / self.X_0)
        half = 0.5 * y
        mid = 0.5 * y
        y_nodes = mid[:, None] + half[:, None] * nodes[None, :]
        y_weights = half[:, None] * weights[None, :]
        X_nodes = self.X_0 * np.exp(y_nodes)
        eta_nodes = np.broadcast_to(eta[:, None], X_nodes.shape)
        stress = self.reference_stress_data(X_nodes, eta_nodes)
        e_a = (1.0 - SELECTED_KAPPA0) * smooth_step(y_nodes / SELECTED_T1)
        d_log_F_difference_dy = 0.5 * e_a * stress["p1_reference"]
        d_U_difference_dy = 0.5 * e_a * X_nodes * stress["ns_reference"]
        return (
            np.sum(y_weights * d_log_F_difference_dy, axis=1),
            np.sum(y_weights * d_U_difference_dy, axis=1),
        )

    def values(
        self, X: Any, eta: Any, *, order: int = QUADRATURE_ORDER
    ) -> dict[str, np.ndarray]:
        """Return vectorized actual F,U,E on the first stress-activation collar."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        reference = self.reference.values(x, e)
        delta_log_F, delta_U = self._activation_integrals(x, e, order=order)
        F = reference["F_reference"] * np.exp(delta_log_F)
        U = reference["U_reference"] + delta_U
        E = np.sqrt(2.0 * x) * F
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "F_activation": F.reshape(shape),
            "U_activation": U.reshape(shape),
            "E_activation": E.reshape(shape),
            "delta_log_F_from_reference": delta_log_F.reshape(shape),
            "delta_U_from_reference": delta_U.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return analytic X derivatives from the public activation ODE."""
        X_arr, eta_arr = self._broadcast(X, eta)
        actual = self.values(X_arr, eta_arr)
        stress = self.reference_stress_data(X_arr, eta_arr)
        schedule = self.activation_schedule(X_arr, eta_arr)
        kappa = schedule["kappa"]
        F_X = -0.5 * kappa * stress["p1_reference"] * actual["F_activation"] / X_arr
        U_X = -0.5 * kappa * stress["ns_reference"]
        return {
            "F_activation_X": F_X,
            "U_activation_X": U_X,
        }

    def shear_schedule(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return the source actual-shear schedule on the materialized collar."""
        X_arr, eta_arr = self._broadcast(X, eta)
        stress = self.reference_stress_data(X_arr, eta_arr)
        actual = self.values(X_arr, eta_arr)
        schedule = self.activation_schedule(X_arr, eta_arr)
        kappa = schedule["kappa"]
        a = kappa * stress["p1_reference"]
        second = kappa * stress["p2_reference"] * stress["E_reference"] / actual["E_activation"]
        t_s = (
            stress["p2_reference"]
            / stress["p1_reference"]
            * stress["E_reference"]
            / actual["E_activation"]
        )
        return {
            "a": a,
            "shear_first": a,
            "shear_second": second,
            "t_s": t_s,
            "kappa": kappa,
            "e_a": schedule["e_a"],
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "selected_kappa0": SELECTED_KAPPA0,
            "selected_t1": SELECTED_T1,
            "quadrature_order": QUADRATURE_ORDER,
            "parent_reference_semantic_sha256": self.reference.semantic_sha256,
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10StressActivationCollar":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError(
                "stress-activation configuration is frozen; serialized data differ "
                "from the registered source/autonomous realization"
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
    ) -> "KokunoPA10StressActivationCollar":
        return cls.from_configuration(json.loads(Path(path).read_text()))

    @property
    def semantic_sha256(self) -> str:
        payload = {
            "schema": SCHEMA,
            "source_commit": SOURCE_COMMIT,
            "source_formulas": _SOURCE_FORMULAS,
            "configuration": self.configuration(),
        }
        return hashlib.sha256(_canonical_json(payload).encode("utf-8")).hexdigest()

    @property
    def truth_boundary(self) -> dict[str, bool]:
        return copy.deepcopy(_TRUTH_BOUNDARY)

    def report(self) -> dict[str, Any]:
        eta_probe = np.asarray([-0.75, -0.25, 0.0, 0.25, 0.75])
        X0 = np.full_like(eta_probe, self.X_0)
        X1 = np.full_like(eta_probe, self.X_1)
        start = self.values(X0, eta_probe)
        end = self.values(X1, eta_probe)
        reference_end = self.reference.values(X1, eta_probe)
        schedule0 = self.activation_schedule(X0, eta_probe)
        schedule1 = self.activation_schedule(X1, eta_probe)
        stress_end = self.reference_stress_data(X1, eta_probe)
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
            "configuration": self.configuration(),
            "semantic_sha256": self.semantic_sha256,
            "geometry": {
                "X_0": self.X_0,
                "X_1": self.X_1,
                "selected_t1": SELECTED_T1,
                "selected_kappa0": SELECTED_KAPPA0,
                "displayed_kappa0_range_satisfied": 0.0 < SELECTED_KAPPA0 < 0.5,
            },
            "endpoint_probe": {
                "eta": eta_probe.tolist(),
                "start_delta_log_F_max": float(np.max(np.abs(start["delta_log_F_from_reference"]))),
                "start_delta_U_max": float(np.max(np.abs(start["delta_U_from_reference"]))),
                "end_delta_log_F_max": float(np.max(np.abs(end["delta_log_F_from_reference"]))),
                "end_delta_U_max": float(np.max(np.abs(end["delta_U_from_reference"]))),
                "end_F_reference_ratio": (end["F_activation"] / reference_end["F_reference"]).tolist(),
                "start_kappa": schedule0["kappa"].tolist(),
                "end_kappa": schedule1["kappa"].tolist(),
                "end_reference_v_r": stress_end["v_reference"].tolist(),
            },
            "truth_boundary": self.truth_boundary,
        }
