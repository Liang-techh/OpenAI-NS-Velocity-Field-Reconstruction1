"""Deterministic witness certificate for the current Kokuno scalar bridge question.

This module stacks on A1 #1233's executable current ``l=-h`` transport.  The
corrected 2026-09-09 reconstruction uses a scalar log-X matching interval in
its source-ideal specialization, but the current repository candidate carries
eta-dependent ``M/M_eta`` and therefore has an eta-resolved ``Q_s`` state.

A single scalar terminal time ``T`` could be valid only if every admissible eta
hits the same target ``Q_p`` at that same ``T``.  To *refute* such a scalar
bridge it is enough to exhibit two eta values whose target roots are each
unique and whose root brackets are disjoint.  This module implements exactly
that one-way certificate:

* a fixed Chebyshev-Lobatto witness grid on the inherited eta interval;
* a sufficient analytic strict-monotonicity test for each fixed eta;
* fixed-mechanics bracketing and bisection of the unique ``Q_s=Q_p`` root;
* a deterministic search for two certified disjoint root brackets.

The monotonicity test follows directly from A1 #1233's closed form.  With
``a=1-h``, ``q0=Q_s(0,eta)`` and ``p=P(eta)``,

    Q_s(y)=exp(-a y) [q0-p+p exp(-h y)]

and hence

    exp(a y) d_y Q_s = -a(q0-p) - p exp(-h y).

If ``p>=0`` it is sufficient that ``q0-p>0``.  If ``p<0`` it is sufficient
that ``a*q0+h*p>0``; in that case the right-hand side is already negative at
y=0 and decreases thereafter.  Together with ``q0>Q_p>0`` these conditions
make the positive target crossing unique.

A positive witness result disproves one common scalar bridge for the current
*executable repository candidate*.  A negative witness result is inconclusive:
it does not establish a common bridge.  In either case this finite witness does
not establish full-eta target totality, smoothness of eta->T(eta), or an
eta-dependent Cartesian terminal surface.  Those remain governed by CR002
#1229/#1241.  This is representation mechanics, not PDE validation and not a
paper-exact/OpenAI-field claim.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_current_exterior_lminus_h_bridge import (
    KokunoCurrentExteriorLMinusHBridge,
)

SCHEMA = "kokuno-current-bridge-common-scalar-refutation-v1"
PARENT_EXACT_HEAD = "6ae44fdac17483b1fd042dca6f5bdf06701ae068"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"

# Fixed representation-mechanics choices.  They are deliberately not public
# scientific tuning knobs and do not define a terminal Cartesian geometry.
_WITNESS_COUNT = 65
_ROOT_EXPANSIONS = 16
_ROOT_BISECTIONS = 96
_MARGIN_REL = 2.0 ** -40

_SOURCE_FORMULAS = {
    "current_bridge": (
        "Q_s(y,eta)=exp(-(1-h)y)[Q_0(eta)-P(eta)(1-exp(-hy))]"
    ),
    "scaled_derivative": (
        "exp((1-h)y)D_y Q_s=-(1-h)(Q_0-P)-P*exp(-hy)"
    ),
    "target": "hold l=-h until Q_s reaches the pre-existing terminal Q_p",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1233 current eta-resolved l=-h transport",
    "eta_witness_grid": "fixed 65-point Chebyshev-Lobatto grid on inherited eta interval",
    "root_mechanics": "fixed 16 bracket expansions plus 96 binary64 bisections",
    "sign_margin": "fixed 2^-40 relative engineering margin for one-way monotonicity witness",
    "promotion_policy": (
        "a disjoint unique-root pair may refute one common scalar bridge; "
        "absence of such a pair does not establish a common bridge"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "current_common_scalar_bridge_witness_check_materialized": True,
    "full_eta_interval_entry_above_qp_established": False,
    "full_eta_interval_target_totality_established": False,
    "full_eta_interval_root_uniqueness_established": False,
    "full_eta_interval_transversality_established": False,
    "smooth_eta_target_time_map_established": False,
    "full_eta_interval_common_scalar_bridge_length_established": False,
    "current_l_minus_h_matching_bridge_materialized": False,
    "eta_dependent_cartesian_matching_boundary_materialized": False,
    "current_cartesian_terminal_multiplier_composed": False,
    "outer_global_leading_velocity_materialized": False,
    "unified_global_cartesian_velocity_export_ready": False,
    "heldout_ns_residual_assessed": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_json(payload: Mapping[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def _strict_decrease_margin(q0: np.ndarray, p0: np.ndarray, h: float) -> np.ndarray:
    """Return a sufficient positive margin for strict decrease on y>=0.

    For p>=0 the scaled derivative is largest as y->infinity, giving the
    condition ``(1-h)(q0-p)>0``.  For p<0 it is largest at y=0, giving
    ``(1-h)q0+h p>0``.  A positive returned value therefore certifies the
    corresponding exact closed-form curve is strictly decreasing.
    """

    q = np.asarray(q0, dtype=float)
    p = np.asarray(p0, dtype=float)
    a = 1.0 - float(h)
    return np.where(p >= 0.0, a * (q - p), a * q + float(h) * p)


@dataclass(frozen=True)
class KokunoCurrentBridgeCommonScalarRefutation:
    """One-way common-scalar refutation witness for exact A1 #1233."""

    parent: KokunoCurrentExteriorLMinusHBridge = field(
        default_factory=KokunoCurrentExteriorLMinusHBridge,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoCurrentExteriorLMinusHBridge):
            raise TypeError("parent must be KokunoCurrentExteriorLMinusHBridge")
        truth = self.parent.truth_boundary
        if not truth["current_l_minus_h_q_s_transport_materialized"]:
            raise ValueError("exact A1 #1233 current l=-h transport is required")
        if truth["current_l_minus_h_matching_bridge_materialized"]:
            raise ValueError("parent unexpectedly already claims a matching bridge")
        if truth["full_eta_interval_common_scalar_bridge_length_established"]:
            raise ValueError("parent unexpectedly already claims one common scalar length")

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.parent.eta_interval)

    @property
    def h_value(self) -> float:
        return float(self.parent.h_value)

    @property
    def q_p(self) -> float:
        return float(self.parent.q_p)

    @property
    def source_formulas(self) -> dict[str, str]:
        return dict(_SOURCE_FORMULAS)

    @property
    def numerical_realization(self) -> dict[str, str]:
        return dict(_NUMERICAL_REALIZATION)

    @property
    def truth_boundary(self) -> dict[str, bool]:
        truth = dict(self.parent.truth_boundary)
        truth.update(_TRUTH_UPDATES)
        return truth

    def witness_eta(self) -> np.ndarray:
        """Return the fixed non-tunable Chebyshev-Lobatto witness grid."""

        lo, hi = self.eta_interval
        k = np.arange(_WITNESS_COUNT, dtype=float)
        nodes = np.cos(math.pi * k / float(_WITNESS_COUNT - 1))
        eta = 0.5 * (lo + hi) + 0.5 * (hi - lo) * nodes
        return np.sort(np.asarray(eta, dtype=float))

    def unique_root_certificate(self, eta: Any) -> dict[str, np.ndarray]:
        """Sufficient fixed-eta uniqueness/transversality certificate.

        This is intentionally a one-way certificate.  ``certified=False`` does
        not mean the root is absent or non-unique; it means only that this
        bounded sufficient condition did not prove it.
        """

        ee = np.asarray(eta, dtype=float)
        lo, hi = self.eta_interval
        if np.any(~np.isfinite(ee)) or np.any((ee < lo) | (ee > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        entry = self.parent.entry_state(ee)
        q0 = np.asarray(entry["Q_s_entry"], dtype=float)
        p0 = np.asarray(entry["P_entry"], dtype=float)
        decrease_margin = _strict_decrease_margin(q0, p0, self.h_value)
        entry_margin = q0 - self.q_p
        asymptotic_coeff = q0 - p0
        scale = np.maximum.reduce(
            [
                np.ones_like(q0),
                np.abs(q0),
                np.abs(p0),
                np.full_like(q0, abs(self.q_p)),
            ]
        )
        numerical_floor = _MARGIN_REL * scale
        certified = (
            (entry_margin > numerical_floor)
            & (decrease_margin > numerical_floor)
            & (asymptotic_coeff > numerical_floor)
        )
        return {
            "eta": ee,
            "Q_s_entry": q0,
            "P_entry": p0,
            "entry_minus_Q_p": entry_margin,
            "strict_decrease_margin": decrease_margin,
            "asymptotic_positive_coefficient": asymptotic_coeff,
            "numerical_margin_floor": numerical_floor,
            "unique_positive_target_root_certified": certified,
        }

    def unique_root_bracket(self, eta: Any) -> dict[str, np.ndarray]:
        """Bracket unique target roots at eta values passing the certificate."""

        cert = self.unique_root_certificate(eta)
        ee = np.asarray(cert["eta"], dtype=float)
        ok = np.asarray(cert["unique_positive_target_root_certified"], dtype=bool)
        if np.any(~ok):
            raise RuntimeError("unique-root sufficient condition failed at requested eta")

        low = np.zeros_like(ee, dtype=float)
        high = np.ones_like(ee, dtype=float)

        def residual(yv: np.ndarray) -> np.ndarray:
            return np.asarray(self.parent.state(yv, ee)["Q_s"], dtype=float) - self.q_p

        f_low = residual(low)
        if np.any(f_low <= 0.0):
            raise RuntimeError("certified bridge entry no longer lies above Q_p")
        f_high = residual(high)
        for _ in range(_ROOT_EXPANSIONS):
            mask = f_high > 0.0
            if not np.any(mask):
                break
            high = np.where(mask, 2.0 * high, high)
            f_high = residual(high)
        if np.any(f_high > 0.0):
            raise RuntimeError("failed to bracket certified unique Q_p root")

        for _ in range(_ROOT_BISECTIONS):
            mid = 0.5 * (low + high)
            f_mid = residual(mid)
            above = f_mid > 0.0
            low = np.where(above, mid, low)
            high = np.where(above, high, mid)

        # Floating-point bisection can stagnate once the endpoints become
        # adjacent representable numbers.  Preserve the sign-oriented bracket
        # semantics rather than claiming arbitrary extra digits.
        f_low = residual(low)
        f_high = residual(high)
        if np.any(f_low < -8.0 * np.finfo(float).eps * max(1.0, abs(self.q_p))):
            raise RuntimeError("lower target bracket lost its above-target sign")
        if np.any(f_high > 8.0 * np.finfo(float).eps * max(1.0, abs(self.q_p))):
            raise RuntimeError("upper target bracket lost its below-target sign")
        return {
            **cert,
            "root_lower": low,
            "root_upper": high,
            "root_midpoint": 0.5 * (low + high),
            "root_width": high - low,
            "lower_residual": f_low,
            "upper_residual": f_high,
        }

    def refutation_report(self) -> dict[str, Any]:
        """Return a deterministic one-way common-scalar bridge certificate."""

        eta = self.witness_eta()
        cert = self.unique_root_certificate(eta)
        ok = np.asarray(cert["unique_positive_target_root_certified"], dtype=bool)
        certified_eta = eta[ok]

        pair: dict[str, Any] | None = None
        gap_lower_bound = 0.0
        certified_count = int(np.count_nonzero(ok))
        if certified_count >= 2:
            brackets = self.unique_root_bracket(certified_eta)
            lower = np.asarray(brackets["root_lower"], dtype=float)
            upper = np.asarray(brackets["root_upper"], dtype=float)
            midpoint = np.asarray(brackets["root_midpoint"], dtype=float)
            i_lo = int(np.argmin(midpoint))
            i_hi = int(np.argmax(midpoint))
            if i_lo != i_hi and upper[i_lo] < lower[i_hi]:
                gap_lower_bound = float(lower[i_hi] - upper[i_lo])
                pair = {
                    "eta_earlier": float(certified_eta[i_lo]),
                    "eta_later": float(certified_eta[i_hi]),
                    "earlier_root_bracket": [float(lower[i_lo]), float(upper[i_lo])],
                    "later_root_bracket": [float(lower[i_hi]), float(upper[i_hi])],
                    "disjoint_gap_lower_bound": gap_lower_bound,
                    "earlier_strict_decrease_margin": float(
                        np.asarray(brackets["strict_decrease_margin"])[i_lo]
                    ),
                    "later_strict_decrease_margin": float(
                        np.asarray(brackets["strict_decrease_margin"])[i_hi]
                    ),
                }

        refuted = pair is not None
        return {
            "schema": "kokuno-agent1-current-common-scalar-refutation-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "eta_interval": list(self.eta_interval),
            "witness_count": int(_WITNESS_COUNT),
            "unique_root_certified_count": certified_count,
            "q_p": self.q_p,
            "common_scalar_bridge_refuted_by_disjoint_unique_root_brackets": refuted,
            "disjoint_root_pair": pair,
            "disjoint_gap_lower_bound": gap_lower_bound,
            "finite_witness_is_sufficient_only_for_refutation": True,
            "absence_of_refutation_does_not_establish_common_scalar_bridge": True,
            "full_eta_interval_entry_above_qp_established": False,
            "full_eta_interval_target_totality_established": False,
            "full_eta_interval_root_uniqueness_established": False,
            "full_eta_interval_transversality_established": False,
            "smooth_eta_target_time_map_established": False,
            "full_eta_interval_common_scalar_bridge_length_established": False,
            "current_l_minus_h_matching_bridge_materialized": False,
            "eta_dependent_cartesian_matching_boundary_materialized": False,
            "truth_boundary": self.truth_boundary,
        }

    def configuration(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "source_repository": SOURCE_REPOSITORY,
            "source_commit": SOURCE_COMMIT,
            "source_blob": SOURCE_BLOB,
            "source_path": SOURCE_PATH,
            "source_release": SOURCE_RELEASE,
            "source_release_date": SOURCE_RELEASE_DATE,
            "witness_count": _WITNESS_COUNT,
            "root_expansions": _ROOT_EXPANSIONS,
            "root_bisections": _ROOT_BISECTIONS,
            "relative_sign_margin": _MARGIN_REL,
            "truth_boundary": self.truth_boundary,
        }

    @property
    def semantic_sha256(self) -> str:
        return hashlib.sha256(_canonical_json(self.configuration()).encode("utf-8")).hexdigest()

    def save_configuration(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(
            json.dumps(self.configuration(), indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        return target

    @classmethod
    def load_configuration(
        cls, path: str | Path
    ) -> "KokunoCurrentBridgeCommonScalarRefutation":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        candidate = cls()
        if payload != candidate.configuration():
            raise ValueError("current bridge scalar-refutation configuration/provenance mismatch")
        return candidate
