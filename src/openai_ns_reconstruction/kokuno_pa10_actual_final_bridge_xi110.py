"""Executable candidate-side Kokuno final bridge from ``X=100`` to ``X_i=110``.

Pinned public provenance: ``KokunoYumeto/yang-mills-interacting-workbench`` at
``143f6773feb424ad9ed3a8d116653200f20346b7``, corrected 2026-09-09 reader.

The corrected reconstruction prescribes the following final Appendix-B bridge.
Starting at ``X=100`` it first multiplies the axial prescription by a smooth
``beta`` decreasing from one to zero on a short logarithmic interval.  It then
holds ``D_X U=0`` and interpolates the angular prescription ``a`` from
``kappa_0 p_{1,r}`` to ``0.8`` on a second short logarithmic interval.  On the
remaining interval through ``X_i=110`` it fixes

    a=0.8,     l=0.6,     D_X U=0.

With ``y=log(X/100)`` and the same profile convention used by the preceding
Agent-1 continuation, the executable profile equations are

    partial_y log(F) = -a/2,
    partial_y U = -kappa_0 beta X n_{s,r}/2       (first interval),
    partial_y U = 0                                (later intervals),
    E = sqrt(2X) F,
    l = 1 + partial_y log(F).

The two short logarithmic widths are source existence choices rather than
published numeric constants.  This repository freezes both to ``0.02`` before
any residual evaluation.  They are explicit autonomous realization choices,
not recovered hidden Kokuno/OpenAI parameters.  The bridge consumes the current
repository-autonomous-pressure reference stress from Agent-1 #933 and the
current autonomous ``kappa_0=0.1`` actual profile from #925.  Consequently this
module materializes a *candidate-side* ``G_i``/``ell_i`` handoff; it does not
claim the source's global smallness/cone inequalities, Appendix-A pressure
identity, a global Cartesian velocity, or PDE validation.
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

from .kokuno_pa10_autonomous_reference_stress_xi110 import (
    SOURCE_X_B_REF,
    SOURCE_X_I,
    KokunoPA10AutonomousReferenceStressToXi110,
)
from .kokuno_pa10_fixed_kappa_x100 import KokunoPA10FixedKappaContinuationToX100
from .kokuno_pa10_physical_center_profile_contract import (
    CORRECTED_RELEASE,
    CORRECTED_RELEASE_DATE,
    SOURCE_COMMIT,
    SOURCE_PATH,
    SOURCE_REPOSITORY,
)
from .kokuno_pa10_reference_continuation_collar import smooth_step
from .kokuno_pa10_reference_p1_x100 import ETA_FD_STEP


SCHEMA = "kokuno-pa10-actual-final-bridge-xi110-v1"
AXIAL_SHUTDOWN_LOG_WIDTH = 0.02
ANGULAR_SETTLE_LOG_WIDTH = 0.02
A_FINAL = 0.8
L_FINAL = 0.6
BRIDGE_QUADRATURE_ORDER = 48

_SOURCE_FORMULAS = {
    "bridge_coordinate": "y=log(X/100), 100<=X<=X_i=110",
    "first_interval": (
        "multiply the axial prescription by a smooth beta from one to zero; "
        "a=kappa_0 p_{1,r}; partial_y log(F)=-a/2; "
        "partial_y U=-kappa_0 beta X n_{s,r}/2"
    ),
    "second_interval": (
        "hold D_X U=0 and interpolate a convexly from kappa_0 p_{1,r} to 0.8; "
        "partial_y log(F)=-a/2"
    ),
    "final_interval": "a=0.8; l=0.6; D_X U=0 through X_i=110",
    "azimuthal_profile": "E=sqrt(2X) F",
    "Xi_boundary": "G_i(eta)=U(X_i,eta); ell_i(eta)=log(C E(X_i,eta))",
}

_NUMERICAL_REALIZATION = {
    "axial_shutdown_log_width": (
        "repository-autonomous 0.02 existence choice, fixed before residual evaluation"
    ),
    "angular_settle_log_width": (
        "repository-autonomous 0.02 existence choice, fixed before residual evaluation"
    ),
    "kappa_0": (
        "inherits repository-autonomous kappa_0=0.1 from the current activation/fixed-kappa chain"
    ),
    "reference_stress": (
        "inherits #933 autonomous-pressure p1/n_s through Xi; not source-prepared Appendix-A stress"
    ),
    "profile_integration": (
        f"piecewise fixed Gauss-Legendre order {BRIDGE_QUADRATURE_ORDER} in y"
    ),
}

_TRUTH_BOUNDARY = {
    "public_reconstruction_source": True,
    "public_final_bridge_formula_executable": True,
    "candidate_side_actual_final_interpolation_100_to_Xi_materialized": True,
    "candidate_side_actual_G_i_at_Xi_materialized": True,
    "candidate_side_actual_ell_i_at_Xi_materialized": True,
    "autonomous_final_widths_fixed_before_residual": True,
    "selected_kappa0_is_repository_autonomous": True,
    "selected_kappa0_chosen_from_ns_residual": False,
    "source_prepared_appendixA_Pi0_materialized": False,
    "source_prepared_full_reference_stress_pair_materialized": False,
    "source_admitted_global_kappa0_smallness": False,
    "final_bridge_cone_admissibility_independently_certified": False,
    "source_hidden_final_widths_recovered": False,
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
class KokunoPA10ActualFinalBridgeToXi110:
    """Continue the current actual PA.10 candidate from ``X=100`` to ``X_i=110``."""

    actual_x100: KokunoPA10FixedKappaContinuationToX100 = field(
        default_factory=KokunoPA10FixedKappaContinuationToX100,
        repr=False,
        compare=False,
    )
    reference_stress: KokunoPA10AutonomousReferenceStressToXi110 = field(
        default_factory=KokunoPA10AutonomousReferenceStressToXi110,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.actual_x100, KokunoPA10FixedKappaContinuationToX100):
            raise TypeError("actual_x100 must be KokunoPA10FixedKappaContinuationToX100")
        if not isinstance(
            self.reference_stress, KokunoPA10AutonomousReferenceStressToXi110
        ):
            raise TypeError(
                "reference_stress must be KokunoPA10AutonomousReferenceStressToXi110"
            )
        if abs(self.actual_x100.X_b_ref - SOURCE_X_B_REF) > 0.0:
            raise ValueError("actual parent must end exactly at X=100")
        if (
            self.reference_stress.parent.semantic_sha256
            != self.actual_x100.reference_stress.semantic_sha256
        ):
            raise ValueError(
                "actual X100 and Xi reference-stress branches do not share one #918 identity"
            )
        if not (
            AXIAL_SHUTDOWN_LOG_WIDTH > 0.0
            and ANGULAR_SETTLE_LOG_WIDTH > 0.0
            and self.y_angular_end < self.y_i
        ):
            raise ValueError("autonomous final bridge widths must fit strictly between 100 and 110")
        if abs((1.0 - 0.5 * A_FINAL) - L_FINAL) > 1e-15:
            raise ValueError("registered final a=0.8 and l=0.6 identity drifted")

    @property
    def X_b_ref(self) -> float:
        return float(SOURCE_X_B_REF)

    @property
    def X_i(self) -> float:
        return float(SOURCE_X_I)

    @property
    def kappa_0(self) -> float:
        return float(self.actual_x100.kappa_0)

    @property
    def C(self) -> float:
        return float(self.actual_x100.activation.C)

    @property
    def y_i(self) -> float:
        return math.log(self.X_i / self.X_b_ref)

    @property
    def y_axial_end(self) -> float:
        return float(AXIAL_SHUTDOWN_LOG_WIDTH)

    @property
    def y_angular_end(self) -> float:
        return float(AXIAL_SHUTDOWN_LOG_WIDTH + ANGULAR_SETTLE_LOG_WIDTH)

    @property
    def X_axial_end(self) -> float:
        return self.X_b_ref * math.exp(self.y_axial_end)

    @property
    def X_angular_end(self) -> float:
        return self.X_b_ref * math.exp(self.y_angular_end)

    @property
    def eta_interval(self) -> tuple[float, float]:
        lo_a, hi_a = self.actual_x100.eta_interval
        lo_r, hi_r = self.reference_stress.eta_interval
        margin = 2.0 * ETA_FD_STEP
        return max(lo_a, lo_r + margin), min(hi_a, hi_r - margin)

    def _broadcast(self, X: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        X_arr, eta_arr = np.broadcast_arrays(_finite(X, "X"), _finite(eta, "eta"))
        x_lo = math.nextafter(self.X_b_ref, -math.inf)
        x_hi = math.nextafter(self.X_i, math.inf)
        if np.any((X_arr < x_lo) | (X_arr > x_hi)):
            raise ValueError(f"X must lie in the actual final bridge [{self.X_b_ref},{self.X_i}]")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError(
                f"eta must lie in the final-bridge executable interval [{lo},{hi}]"
            )
        return np.clip(X_arr, self.X_b_ref, self.X_i), eta_arr

    def _reference_values(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        return self.reference_stress.values(X, eta)

    def schedule(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return the public piecewise ``beta/a/l/D_XU`` bridge prescription."""
        X_arr, eta_arr = self._broadcast(X, eta)
        y = np.log(X_arr / self.X_b_ref)
        stress = self._reference_values(X_arr, eta_arr)
        p1 = stress["p1_reference"]
        ns = stress["n_s_reference_autonomous"]

        beta = np.zeros_like(y, dtype=float)
        sigma_a = np.ones_like(y, dtype=float)
        a = np.full_like(y, A_FINAL, dtype=float)
        dU_dy = np.zeros_like(y, dtype=float)

        first = y <= self.y_axial_end
        if np.any(first):
            s1 = np.clip(y[first] / AXIAL_SHUTDOWN_LOG_WIDTH, 0.0, 1.0)
            beta[first] = 1.0 - smooth_step(s1)
            sigma_a[first] = 0.0
            a[first] = self.kappa_0 * p1[first]
            dU_dy[first] = -0.5 * self.kappa_0 * beta[first] * X_arr[first] * ns[first]

        second = (y > self.y_axial_end) & (y <= self.y_angular_end)
        if np.any(second):
            s2 = np.clip(
                (y[second] - self.y_axial_end) / ANGULAR_SETTLE_LOG_WIDTH,
                0.0,
                1.0,
            )
            sigma = smooth_step(s2)
            sigma_a[second] = sigma
            a[second] = (1.0 - sigma) * self.kappa_0 * p1[second] + sigma * A_FINAL

        final = y > self.y_angular_end
        if np.any(final):
            beta[final] = 0.0
            sigma_a[final] = 1.0
            a[final] = A_FINAL

        # Enforce exact public endpoint states despite logarithm roundoff.
        at_start = X_arr <= self.X_b_ref
        beta = np.where(at_start, 1.0, beta)
        sigma_a = np.where(at_start, 0.0, sigma_a)
        a = np.where(at_start, self.kappa_0 * p1, a)
        dU_dy = np.where(at_start, -0.5 * self.kappa_0 * X_arr * ns, dU_dy)

        at_or_after_axial = X_arr >= self.X_axial_end
        beta = np.where(at_or_after_axial, 0.0, beta)
        dU_dy = np.where(at_or_after_axial, 0.0, dU_dy)
        at_or_after_angular = X_arr >= self.X_angular_end
        sigma_a = np.where(at_or_after_angular, 1.0, sigma_a)
        a = np.where(at_or_after_angular, A_FINAL, a)

        dlogF_dy = -0.5 * a
        l_profile = 1.0 + dlogF_dy
        D_X_U = dU_dy
        return {
            "X": X_arr,
            "eta": eta_arr,
            "y": y,
            "beta": beta,
            "angular_sigma": sigma_a,
            "a": a,
            "l_profile": l_profile,
            "dlogF_dy": dlogF_dy,
            "D_X_U": D_X_U,
        }

    def _integrate_piece(
        self,
        y_lo: float,
        y_hi: float,
        eta: float,
        *,
        order: int,
    ) -> tuple[float, float]:
        if y_hi <= y_lo:
            return 0.0, 0.0
        nodes, weights = _legendre_rule(order)
        half = 0.5 * (y_hi - y_lo)
        mid = 0.5 * (y_hi + y_lo)
        y_nodes = mid + half * nodes
        X_nodes = self.X_b_ref * np.exp(y_nodes)
        eta_nodes = np.full_like(X_nodes, eta, dtype=float)
        sched = self.schedule(X_nodes, eta_nodes)
        delta_log_F = half * float(np.dot(weights, sched["dlogF_dy"]))
        delta_U = half * float(np.dot(weights, sched["D_X_U"]))
        return delta_log_F, delta_U

    def _bridge_integrals(
        self,
        X: np.ndarray,
        eta: np.ndarray,
        *,
        order: int = BRIDGE_QUADRATURE_ORDER,
    ) -> tuple[np.ndarray, np.ndarray]:
        if X.ndim != 1 or eta.ndim != 1 or X.shape != eta.shape:
            raise ValueError("final-bridge integration expects matching flat vectors")
        dlog = np.zeros_like(X, dtype=float)
        dU = np.zeros_like(X, dtype=float)
        for index, (xx, ee) in enumerate(zip(X, eta, strict=True)):
            y = math.log(float(xx) / self.X_b_ref)
            total_log = 0.0
            total_U = 0.0
            stop1 = min(y, self.y_axial_end)
            qlog, qU = self._integrate_piece(0.0, stop1, float(ee), order=order)
            total_log += qlog
            total_U += qU
            if y > self.y_axial_end:
                stop2 = min(y, self.y_angular_end)
                qlog, qU = self._integrate_piece(
                    self.y_axial_end, stop2, float(ee), order=order
                )
                total_log += qlog
                total_U += qU
            if y > self.y_angular_end:
                # The final source interval is exactly constant: a=0.8, D_XU=0.
                total_log += -0.5 * A_FINAL * (y - self.y_angular_end)
            dlog[index] = total_log
            dU[index] = total_U
        return dlog, dU

    def x100_handoff(self, eta: Any) -> dict[str, np.ndarray]:
        eta_arr = _finite(eta, "eta")
        lo, hi = self.eta_interval
        if np.any((eta_arr < lo) | (eta_arr > hi)):
            raise ValueError("eta outside final-bridge executable interval")
        return self.actual_x100.handoff_at_X100(eta_arr)

    def values(
        self, X: Any, eta: Any, *, order: int = BRIDGE_QUADRATURE_ORDER
    ) -> dict[str, np.ndarray]:
        """Return vectorized actual ``F,U,E`` on ``100<=X<=110``."""
        X_arr, eta_arr = self._broadcast(X, eta)
        shape = X_arr.shape
        x = X_arr.reshape(-1)
        e = eta_arr.reshape(-1)
        start = self.x100_handoff(e)
        delta_log_F, delta_U = self._bridge_integrals(x, e, order=order)
        F = start["F_fixed_kappa"] * np.exp(delta_log_F)
        U = start["U_fixed_kappa"] + delta_U
        E = np.sqrt(2.0 * x) * F
        if np.any(F <= 0.0) or np.any(E <= 0.0):
            raise RuntimeError("actual final bridge lost positive F/E")
        if any(np.any(~np.isfinite(v)) for v in (F, U, E)):
            raise RuntimeError("actual final bridge produced non-finite profile values")
        return {
            "X": x.reshape(shape),
            "eta": e.reshape(shape),
            "F_final_bridge": F.reshape(shape),
            "U_final_bridge": U.reshape(shape),
            "E_final_bridge": E.reshape(shape),
            "delta_log_F_from_X100": delta_log_F.reshape(shape),
            "delta_U_from_X100": delta_U.reshape(shape),
        }

    def radial_derivatives(self, X: Any, eta: Any) -> dict[str, np.ndarray]:
        """Return first radial derivatives directly from the public bridge schedule."""
        X_arr, eta_arr = self._broadcast(X, eta)
        values = self.values(X_arr, eta_arr)
        sched = self.schedule(X_arr, eta_arr)
        F = values["F_final_bridge"]
        F_X = F * sched["dlogF_dy"] / X_arr
        U_X = sched["D_X_U"] / X_arr
        E_X = F / np.sqrt(2.0 * X_arr) + np.sqrt(2.0 * X_arr) * F_X
        return {
            "F_final_bridge_X": F_X,
            "U_final_bridge_X": U_X,
            "E_final_bridge_X": E_X,
        }

    def handoff_at_Xi(self, eta: Any) -> dict[str, np.ndarray]:
        eta_arr = _finite(eta, "eta")
        X = np.full_like(eta_arr, self.X_i, dtype=float)
        values = self.values(X, eta_arr)
        derivatives = self.radial_derivatives(X, eta_arr)
        sched = self.schedule(X, eta_arr)
        ell_i = np.log(self.C * values["E_final_bridge"])
        G_i = values["U_final_bridge"]
        return {
            **values,
            **derivatives,
            **sched,
            "G_i": G_i,
            "ell_i": ell_i,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "actual_x100_semantic_sha256": self.actual_x100.semantic_sha256,
            "reference_stress_xi_semantic_sha256": self.reference_stress.semantic_sha256,
            "source_X_b_ref": self.X_b_ref,
            "source_X_i": self.X_i,
            "selected_kappa0": self.kappa_0,
            "axial_shutdown_log_width": AXIAL_SHUTDOWN_LOG_WIDTH,
            "angular_settle_log_width": ANGULAR_SETTLE_LOG_WIDTH,
            "bridge_quadrature_order": BRIDGE_QUADRATURE_ORDER,
            "a_final": A_FINAL,
            "l_final": L_FINAL,
        }

    @classmethod
    def from_configuration(cls, payload: Mapping[str, Any]) -> "KokunoPA10ActualFinalBridgeToXi110":
        if not isinstance(payload, Mapping):
            raise TypeError("configuration must be a mapping")
        expected = cls().configuration()
        if dict(payload) != expected:
            raise ValueError("actual final-bridge configuration is frozen; serialized data differ")
        return cls()

    def save_configuration(self, path: str | Path) -> dict[str, Any]:
        payload = self.configuration()
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        return payload

    @classmethod
    def load_configuration(cls, path: str | Path) -> "KokunoPA10ActualFinalBridgeToXi110":
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
        eta_probe = np.asarray([-0.65, -0.25, 0.0, 0.25, 0.65])
        start = self.x100_handoff(eta_probe)
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
            "geometry": {
                "X_b_ref": self.X_b_ref,
                "X_axial_end": self.X_axial_end,
                "X_angular_end": self.X_angular_end,
                "X_i": self.X_i,
                "remaining_final_log_width": self.y_i - self.y_angular_end,
            },
            "X100_probe": {
                "eta": eta_probe.tolist(),
                "F": start["F_fixed_kappa"].tolist(),
                "U": start["U_fixed_kappa"].tolist(),
            },
            "Xi_probe": {
                "eta": eta_probe.tolist(),
                "F": end["F_final_bridge"].tolist(),
                "U": end["U_final_bridge"].tolist(),
                "E": end["E_final_bridge"].tolist(),
                "G_i": end["G_i"].tolist(),
                "ell_i": end["ell_i"].tolist(),
                "a": end["a"].tolist(),
                "l_profile": end["l_profile"].tolist(),
            },
            "semantic_sha256": self.semantic_sha256,
            "truth_boundary": self.truth_boundary,
        }
