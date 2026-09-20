"""Executable first reference-continuation collar for the Kokuno leading profiles.

Pinned public provenance: KokunoYumeto/yang-mills-interacting-workbench at
143f6773feb424ad9ed3a8d116653200f20346b7, corrected 2026-09-09 reader.

The public continuation derivation fixes

    X_0 = 4/Lambda,   X_b^ref = 100,   X_i = 110,
    y = log(X/X_0),

and the exact flat smooth step

    sigma(s) = exp(-1/s^2) /
               [exp(-1/s^2) + exp(-1/(1-s)^2)], 0<s<1,

extended by 0/1 outside [0,1].  For 0<t_1<=bar t with
4 exp(2 bar t)<4.1 it retains the natural analytic profile through y=t_1,
then on t_1<y<2t_1 prescribes

    d_y log(phi_r) =
        [1-sigma((y-t_1)/t_1)] D_X log(phi_nat),
    d_y U_r =
        [1-sigma((y-t_1)/t_1)] D_X U_nat,

and makes the reference profiles constant in y after 2t_1.

This module materializes only that first collar.  Because C is constant,
the same logarithmic equation holds for F=phi/C.  The collar consumes the
already executable source-C-normalized PA.10 center as the natural profile,
uses a fixed repository choice t_1=0.005 satisfying the displayed geometric
source bound, and evaluates the two source ODE integrals by deterministic
Gauss-Legendre quadrature.  t_1 is not fitted to a residual and is not claimed
to satisfy every later smallness inequality in the proof.

The actual stress activation, kappa_0 continuation through X=100, final
interpolation to X_i=110, shaping to X_sep, -8/-7 axial restoration, five-
moment repair, matched pressure/forcing, and global Cartesian field are not
materialized here.  This is a source-governed executable profile increment,
not PDE validation and not an OpenAI/paper-exact field claim.
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
from .kokuno_pa10_source_c_normalized_physical_center import (
    KokunoPA10SourceCNormalizedPhysicalCenter,
)


SCHEMA = "kokuno-pa10-reference-continuation-collar-v1"
SOURCE_X0_TIMES_LAMBDA = 4.0
SOURCE_ANALYTIC_Y_UPPER_TIMES_LAMBDA = 4.1
SOURCE_X_B_REF = 100.0
SOURCE_X_I = 110.0
SELECTED_T1 = 0.005
QUADRATURE_ORDER = 64

_SOURCE_FORMULAS = {
    "X0": "X_0=4/Lambda",
    "reference_cutoff": "X_b^ref=100",
    "Xi": "X_i=110",
    "coordinate": "y=log(X/X_0)",
    "smooth_step": (
        "sigma(s)=exp(-1/s^2)/(exp(-1/s^2)+exp(-1/(1-s)^2)), "
        "0<s<1; sigma=0 for s<=0 and 1 for s>=1"
    ),
    "geometric_bound": "0<t_1<=bar t and 4 exp(2 bar t)<4.1",
    "F_collar": (
        "d_y log F_r=(1-sigma((y-t_1)/t_1)) D_X log F_nat "
        "(equivalent to the public phi equation because C is constant)"
    ),
    "U_collar": (
        "d_y U_r=(1-sigma((y-t_1)/t_1)) D_X U_nat"
    ),
    "after_collar": "reference F_r,U_r are constant in y after 2 t_1",
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_first_reference_collar_formula_executable": True,
    "exact_public_smooth_step_executable": True,
    "source_X0_4_over_Lambda_executable": True,
    "source_reference_cutoff_100_recorded": True,
    "source_Xi_110_recorded": True,
    "selected_t1_is_repository_autonomous": True,
    "selected_t1_chosen_from_ns_residual": False,
    "selected_t1_satisfies_displayed_geometric_source_bound": True,
    "selected_t1_full_later_smallness_admitted": False,
    "natural_to_reference_value_and_radial_jet_handoff_executable": True,
    "stress_activation_materialized": False,
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


def smooth_step(s: Any) -> np.ndarray:
    """Return the exact source flat step with stable endpoint evaluation."""
    z = _finite(s, "s")
    out = np.empty_like(z, dtype=float)
    low = z <= 0.0
    high = z >= 1.0
    interior = ~(low | high)
    out[low] = 0.0
    out[high] = 1.0
    if np.any(interior):
        q = z[interior]
        log_a = -1.0 / (q * q)
        one_minus = 1.0 - q
        log_b = -1.0 / (one_minus * one_minus)
        d = log_b - log_a
        vals = np.empty_like(d)
        nonnegative = d >= 0.0
        if np.any(nonnegative):
            e = np.exp(-d[nonnegative])
            vals[nonnegative] = e / (1.0 + e)
        if np.any(~nonnegative):
            e = np.exp(d[~nonnegative])
            vals[~nonnegative] = 1.0 / (1.0 + e)
        out[interior] = vals
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
class KokunoPA10ReferenceContinuationCollar:
    """First public radial reference collar applied to the executable PA.10 center."""

    source_normalized: KokunoPA10SourceCNormalizedPhysicalCenter = field(
        default_factory=KokunoPA10SourceCNormalizedPhysicalCenter,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(
            self.source_normalized, KokunoPA10SourceCNormalizedPhysicalCenter
        ):
            raise TypeError(
                "source_normalized must be KokunoPA10SourceCNormalizedPhysicalCenter"
            )
        if not self.geometric_source_bound_satisfied:
            raise ValueError("selected t_1 violates the displayed source collar bound")
        inner_max = float(self.source_normalized.physical_profiles.source_X_interval[1])
        if self.X_2 > inner_max:
            raise ValueError(
                "selected t_1 reaches beyond the executable natural PA.10 interval"
            )

    @property
    def Lambda(self) -> float:
        return float(self.source_normalized.Lambda)

    @property
    def C(self) -> float:
        return float(self.source_normalized.C)

    @property
    def X_0(self) -> float:
        return SOURCE_X0_TIMES_LAMBDA / self.Lambda

    @property
    def X_1(self) -> float:
        return self.X_0 * math.exp(SELECTED_T1)

    @property
    def X_2(self) -> float:
        return self.X_0 * math.exp(2.0 * SELECTED_T1)

    @property
    def geometric_t1_upper_bound(self) -> float:
        return 0.5 * math.log(
            SOURCE_ANALYTIC_Y_UPPER_TIMES_LAMBDA / SOURCE_X0_TIMES_LAMBDA
        )

    @property
    def geometric_source_bound_satisfied(self) -> bool:
        return bool(
            0.0 < SELECTED_T1 < self.geometric_t1_upper_bound
            and SOURCE_X0_TIMES_LAMBDA * math.exp(2.0 * SELECTED_T1)
            < SOURCE_ANALYTIC_Y_UPPER_TIMES_LAMBDA
        )

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(
            float(v)
            for v in self.source_normalized.physical_profiles.eta_interval
        )

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        if np.any(X_arr < 0.0) or np.any(X_arr > self.X_2):
            raise ValueError(
                f"X must lie in the first reference-collar domain [0,{self.X_2}]"
            )
        eta_lo, eta_hi = self.eta_interval
        if np.any(eta_arr < eta_lo) or np.any(eta_arr > eta_hi):
            raise ValueError(
                f"eta must lie in the executable natural interval [{eta_lo},{eta_hi}]"
            )
        return X_arr, eta_arr

    def _transition_integrals(
        self, X: np.ndarray, eta: np.ndarray, *, order: int = QUADRATURE_ORDER
    ) -> tuple[np.ndarray, np.ndarray]:
        """Integrate the two displayed source collar ODEs from X_1 to X."""
        if X.ndim != 1 or eta.ndim != 1 or X.shape != eta.shape:
            raise ValueError("transition integration expects matching flat vectors")
        if np.any(X < self.X_1) or np.any(X > self.X_2):
            raise ValueError("transition integration is restricted to [X_1,X_2]")
        if X.size == 0:
            return np.empty(0), np.empty(0)

        nodes, weights = _legendre_rule(order)
        y = np.log(X / self.X_0)
        half = 0.5 * (y - SELECTED_T1)
        mid = 0.5 * (y + SELECTED_T1)
        y_nodes = mid[:, None] + half[:, None] * nodes[None, :]
        y_weights = half[:, None] * weights[None, :]
        X_nodes = self.X_0 * np.exp(y_nodes)
        eta_nodes = np.broadcast_to(eta[:, None], X_nodes.shape)

        natural = self.source_normalized.values(X_nodes, eta_nodes)
        natural_d = self.source_normalized.derivatives(X_nodes, eta_nodes)
        F_nat = natural["F_0"]
        if np.any(F_nat <= 0.0):
            raise RuntimeError(
                "natural F must stay strictly positive on the source collar"
            )
        gate = 1.0 - smooth_step(
            (y_nodes - SELECTED_T1) / SELECTED_T1
        )
        d_log_F_dy = gate * X_nodes * natural_d["F_0_X"] / F_nat
        d_U_dy = gate * X_nodes * natural_d["U_0_X"]
        return (
            np.sum(y_weights * d_log_F_dy, axis=1),
            np.sum(y_weights * d_U_dy, axis=1),
        )

    def values(
        self, X: Any, eta: Any, *, order: int = QUADRATURE_ORDER
    ) -> dict[str, np.ndarray]:
        """Return vectorized F_r,U_r,E_r on the first public reference collar."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        F = np.empty_like(x)
        U = np.empty_like(x)

        natural_mask = x <= self.X_1
        if np.any(natural_mask):
            natural = self.source_normalized.values(x[natural_mask], e[natural_mask])
            F[natural_mask] = natural["F_0"]
            U[natural_mask] = natural["U_0"]

        transition_mask = ~natural_mask
        if np.any(transition_mask):
            xt = x[transition_mask]
            et = e[transition_mask]
            start = self.source_normalized.values(
                np.full_like(xt, self.X_1), et
            )
            if np.any(start["F_0"] <= 0.0):
                raise RuntimeError("natural F at X_1 must be strictly positive")
            int_log_F, int_U = self._transition_integrals(
                xt, et, order=order
            )
            F[transition_mask] = start["F_0"] * np.exp(int_log_F)
            U[transition_mask] = start["U_0"] + int_U

        E = np.sqrt(2.0 * x) * F
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "F_reference": F.reshape(shape),
            "E_reference": E.reshape(shape),
            "U_reference": U.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return analytic radial derivatives implied by the source collar ODE."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        F_X = np.empty_like(x)
        U_X = np.empty_like(x)

        natural_mask = x <= self.X_1
        if np.any(natural_mask):
            natural_d = self.source_normalized.derivatives(
                x[natural_mask], e[natural_mask]
            )
            F_X[natural_mask] = natural_d["F_0_X"]
            U_X[natural_mask] = natural_d["U_0_X"]

        transition_mask = ~natural_mask
        if np.any(transition_mask):
            xt = x[transition_mask]
            et = e[transition_mask]
            ref = self.values(xt, et)
            natural = self.source_normalized.values(xt, et)
            natural_d = self.source_normalized.derivatives(xt, et)
            if np.any(natural["F_0"] <= 0.0):
                raise RuntimeError("natural F must stay positive in transition")
            y = np.log(xt / self.X_0)
            gate = 1.0 - smooth_step((y - SELECTED_T1) / SELECTED_T1)
            F_X[transition_mask] = (
                ref["F_reference"]
                * gate
                * natural_d["F_0_X"]
                / natural["F_0"]
            )
            U_X[transition_mask] = gate * natural_d["U_0_X"]

        return {
            "F_reference_X": F_X.reshape(shape),
            "U_reference_X": U_X.reshape(shape),
        }

    def handoff_at_collar_end(self, eta: Any) -> dict[str, np.ndarray]:
        """Return the value/zero-radial-jet handoff at y=2t_1."""
        eta_arr = _finite(eta, "eta")
        X = np.full_like(eta_arr, self.X_2, dtype=float)
        values = self.values(X, eta_arr)
        derivatives = self.radial_derivatives(X, eta_arr)
        return {**values, **derivatives}

    def configuration(self) -> dict[str, Any]:
        upstream = self.source_normalized.physical_profiles.configuration()
        return {
            "schema": SCHEMA,
            "selected_t1": SELECTED_T1,
            "quadrature_order": QUADRATURE_ORDER,
            "source_X0_times_Lambda": SOURCE_X0_TIMES_LAMBDA,
            "source_X_b_ref": SOURCE_X_B_REF,
            "source_X_i": SOURCE_X_I,
            "upstream_physical_configuration_sha256": hashlib.sha256(
                _canonical_json(upstream).encode("utf-8")
            ).hexdigest(),
        }

    @classmethod
    def from_configuration(
        cls, payload: Mapping[str, Any]
    ) -> "KokunoPA10ReferenceContinuationCollar":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError(
                "reference-collar configuration is frozen; serialized data differ "
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
    ) -> "KokunoPA10ReferenceContinuationCollar":
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
        handoff = self.handoff_at_collar_end(eta_probe)
        radial_max = max(
            float(np.max(np.abs(handoff["F_reference_X"]))),
            float(np.max(np.abs(handoff["U_reference_X"]))),
        )
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
                "Lambda": self.Lambda,
                "X_0": self.X_0,
                "X_1": self.X_1,
                "X_2": self.X_2,
                "geometric_t1_upper_bound": self.geometric_t1_upper_bound,
                "selected_t1": SELECTED_T1,
                "geometric_source_bound_satisfied": self.geometric_source_bound_satisfied,
            },
            "handoff_probe": {
                "eta": eta_probe.tolist(),
                "F_reference": handoff["F_reference"].tolist(),
                "U_reference": handoff["U_reference"].tolist(),
                "max_abs_radial_derivative_at_X2": radial_max,
            },
            "truth_boundary": self.truth_boundary,
        }
