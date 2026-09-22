"""Current-candidate Kokuno ``l=-h`` exterior ``Q_s`` bridge transport.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
``navier-stokes/navier_stokes_workbench.tex``
blob 205a99807302e21a51c5eaf223390c0dfc42bcd0
corrected 2026-09-09 reconstruction.

After the public ``l:-1 -> -h`` release, the corrected reconstruction holds
``l=-h`` until the exterior radial quantity reaches the terminal target ``Q_p``.
For the source-ideal specialization ``M=0`` this reduces to a homogeneous
exponential.  The current A1 lineage deliberately carries the actual post-pulse
``M/M_eta`` state, so its general transport is instead

    Q_s' + (1-h) Q_s = h (W-1),

with ``U=0`` and ``(log E)_eta=0`` on the bridge.  Physical ``M,M_eta`` remain
constant, hence if

    P(eta) = 2 D eta (M/X)_0 + (1-eta^2) (M_eta/X)_0 = 1-W_0,

then

    W(y)-1 = -P(eta) exp(-y)

and the exact current bridge transport is

    Q_s(y,eta) = exp(-(1-h)y)
                   [Q_0(eta) - P(eta) (1-exp(-h y))].

This module materializes that eta-resolved current transport and deterministic
pointwise target-time diagnostics.  It intentionally does *not* replace the
eta-resolved state by a mean/min/max/RMS scalar bridge length, does not define an
eta-dependent Cartesian terminal surface, and does not extend the public
Cartesian velocity.  Those are separate representation choices governed by
CR002 #1229.

This is executable candidate-state transport, not a paper-exact field and not
Navier--Stokes validation.  No pressure, forcing, residual, optimizer, hidden
parameter, or held-out tuning enters here.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_current_exterior_qs_release2_state import (
    KokunoCurrentExteriorQsRelease2State,
)

SCHEMA = "kokuno-current-exterior-lminus-h-bridge-v1"
PARENT_EXACT_HEAD = "8503f2ede2dd604392198b31ffdeb79321822a30"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"

# Fixed internal numerical mechanics for pointwise diagnostics only.  These are
# not caller-tunable scientific parameters and do not define a global bridge.
_ROOT_EXPANSIONS = 12
_ROOT_BISECTIONS = 96

_SOURCE_FORMULAS = {
    "general_transport": (
        "D_logX Q_s+(1+l)Q_s=-W*l-h(1-2etaU)-H_c*(log E)_eta"
    ),
    "bridge_specialization": (
        "l=-h, U=0, (log E)_eta=0 => Q_s'+(1-h)Q_s=h(W-1)"
    ),
    "carried_primitives": (
        "physical M,M_eta constant for U=0; M/X and M_eta/X decay as exp(-y)"
    ),
    "current_closed_form": (
        "Q_s(y,eta)=exp(-(1-h)y)[Q_0(eta)-P(eta)(1-exp(-h y))], "
        "P=2D eta(M/X)_0+(1-eta^2)(M_eta/X)_0=1-W_0"
    ),
    "terminal_order": "hold l=-h until Q_s reaches the already-defined Q_p",
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1225 current release2-end Q_s state",
    "bridge_arithmetic": "analytic closed form in binary64 NumPy/Python arithmetic",
    "pointwise_target_time": (
        "deterministic bracket expansion plus 96 fixed bisections; diagnostic only"
    ),
    "scalar_bridge_policy": (
        "no mean/min/max/RMS or representative-eta bridge length is promoted"
    ),
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "current_q_s_release2_endpoint_materialized": True,
    "current_l_minus_h_q_s_transport_materialized": True,
    "current_eta_resolved_pointwise_target_time_diagnostic_materialized": True,
    "full_eta_interval_common_scalar_bridge_length_established": False,
    "current_l_minus_h_matching_bridge_materialized": False,
    "eta_dependent_cartesian_matching_boundary_materialized": False,
    "current_cartesian_terminal_multiplier_composed": False,
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


def _semantic_value(obj: Any) -> str:
    value = getattr(obj, "semantic_sha256")
    return str(value() if callable(value) else value)


def _finite_nonnegative(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any(out < 0.0):
        raise ValueError(f"{name} must be nonnegative")
    return out


@dataclass(frozen=True)
class KokunoCurrentExteriorLMinusHBridge:
    """Exact eta-resolved current ``l=-h`` radial transport after A1 #1225."""

    parent: KokunoCurrentExteriorQsRelease2State = field(
        default_factory=KokunoCurrentExteriorQsRelease2State,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoCurrentExteriorQsRelease2State):
            raise TypeError("parent must be KokunoCurrentExteriorQsRelease2State")
        truth = self.parent.truth_boundary
        if not truth["current_q_s_release2_endpoint_materialized"]:
            raise ValueError("exact A1 #1225 current Q_s endpoint is required")
        if truth["current_l_minus_h_matching_bridge_materialized"]:
            raise ValueError("parent unexpectedly already claims the matching bridge")
        if not (0.0 < self.h_value < 0.01):
            raise ValueError("current h must satisfy the public range 0<h<1/100")
        if not (math.isfinite(self.q_p) and self.q_p > 0.0):
            raise ValueError("terminal Q_p must be finite and positive")

    @property
    def h_value(self) -> float:
        return float(self.parent.h_value)

    @property
    def D(self) -> float:
        return float(self.parent.D)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.parent.eta_interval)

    @property
    def q_p(self) -> float:
        return float(self.parent.parent.q_p)

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

    def _broadcast(self, y: Any, eta: Any) -> tuple[np.ndarray, np.ndarray]:
        yy = _finite_nonnegative(y, "bridge y")
        ee = np.asarray(eta, dtype=float)
        if np.any(~np.isfinite(ee)):
            raise ValueError("eta must contain only finite values")
        lo, hi = self.eta_interval
        if np.any((ee < lo) | (ee > hi)):
            raise ValueError(f"eta must lie in [{lo},{hi}]")
        return np.broadcast_arrays(yy, ee)

    def entry_state(self, eta: Any) -> dict[str, np.ndarray]:
        """Return the exact A1 #1225 bridge-entry state and ``P=1-W``."""
        ee = np.asarray(eta, dtype=float)
        state = self.parent.state(ee)
        q0 = np.asarray(state["Q_s_current_release2_end"], dtype=float)
        m0 = np.asarray(state["M_over_X"], dtype=float)
        me0 = np.asarray(state["M_eta_over_X"], dtype=float)
        w0 = np.asarray(state["W"], dtype=float)
        d = 1.0 - ee * ee
        p_direct = 2.0 * self.D * ee * m0 + d * me0
        p_from_w = 1.0 - w0
        scale = np.maximum.reduce(
            [np.abs(p_direct), np.abs(p_from_w), np.full_like(p_direct, 1.0e-300)]
        )
        if np.any(np.abs(p_direct - p_from_w) > 8.0e-14 * scale + 1.0e-300):
            raise RuntimeError("current bridge P no longer replays 1-W from #1225")
        arrays = (q0, m0, me0, w0, p_direct)
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("current bridge entry state became non-finite")
        return {
            "eta": ee,
            "Q_s_entry": q0,
            "M_over_X_entry": m0,
            "M_eta_over_X_entry": me0,
            "W_entry": w0,
            "P_entry": p_direct,
        }

    def state(self, y: Any, eta: Any) -> dict[str, np.ndarray]:
        """Evaluate the exact current ``l=-h`` transport at ``(y,eta)``."""
        yy, ee = self._broadcast(y, eta)
        entry = self.entry_state(ee)
        q0 = np.asarray(entry["Q_s_entry"], dtype=float)
        m0 = np.asarray(entry["M_over_X_entry"], dtype=float)
        me0 = np.asarray(entry["M_eta_over_X_entry"], dtype=float)
        p0 = np.asarray(entry["P_entry"], dtype=float)

        e1 = np.exp(-yy)
        eh = np.exp(-self.h_value * yy)
        ea = np.exp(-(1.0 - self.h_value) * yy)
        m = m0 * e1
        me = me0 * e1
        p = p0 * e1
        w = 1.0 - p
        q = ea * (q0 - p0 * (1.0 - eh))
        dq = -(1.0 - self.h_value) * q - self.h_value * p

        arrays = (m, me, p, w, q, dq)
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("current l=-h bridge state became non-finite")
        return {
            "y": yy,
            "eta": ee,
            "M_over_X": m,
            "M_eta_over_X": me,
            "P": p,
            "W": w,
            "Q_s": q,
            "D_y_Q_s": dq,
            "Q_p": np.full_like(q, self.q_p, dtype=float),
        }

    def pointwise_target_time(self, eta: Any) -> np.ndarray:
        """Diagnostic eta-wise first target time; does not define a global bridge."""
        ee = np.asarray(eta, dtype=float)
        lo_eta, hi_eta = self.eta_interval
        if np.any(~np.isfinite(ee)) or np.any((ee < lo_eta) | (ee > hi_eta)):
            raise ValueError(f"eta must lie in [{lo_eta},{hi_eta}]")
        entry = self.entry_state(ee)
        q0 = np.asarray(entry["Q_s_entry"], dtype=float)
        if np.any(q0 <= self.q_p):
            raise RuntimeError("current bridge entry must lie strictly above terminal Q_p")

        low = np.zeros_like(q0, dtype=float)
        # Start from the public source-ideal matching scale, but only as a
        # deterministic search bracket.  It is never promoted as current time.
        q_ideal = float(self.parent.parent.q_in_source_ideal)
        source_scale = math.log(q_ideal / self.q_p) / (1.0 - self.h_value)
        high = np.full_like(q0, max(1.0, source_scale), dtype=float)

        def residual(yv: np.ndarray) -> np.ndarray:
            return np.asarray(self.state(yv, ee)["Q_s"], dtype=float) - self.q_p

        f_high = residual(high)
        for _ in range(_ROOT_EXPANSIONS):
            mask = f_high > 0.0
            if not np.any(mask):
                break
            high = np.where(mask, 2.0 * high, high)
            f_high = residual(high)
        if np.any(f_high > 0.0):
            raise RuntimeError("failed to bracket current pointwise Q_p target")

        for _ in range(_ROOT_BISECTIONS):
            mid = 0.5 * (low + high)
            f_mid = residual(mid)
            above = f_mid > 0.0
            low = np.where(above, mid, low)
            high = np.where(above, high, mid)
        return 0.5 * (low + high)

    def pointwise_target_report(
        self,
        eta: Any = (-0.83, -0.61, -0.37, 0.0, 0.19, 0.73),
    ) -> dict[str, Any]:
        ee = np.asarray(eta, dtype=float)
        times = np.asarray(self.pointwise_target_time(ee), dtype=float)
        hit = np.asarray(self.state(times, ee)["Q_s"], dtype=float)
        residual = hit - self.q_p
        q0 = np.asarray(self.entry_state(ee)["Q_s_entry"], dtype=float)
        source_ideal_time = math.log(
            float(self.parent.parent.q_in_source_ideal) / self.q_p
        ) / (1.0 - self.h_value)
        return {
            "schema": "kokuno-agent1-current-lminus-h-pointwise-target-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "eta": ee.tolist(),
            "q_s_entry": q0.tolist(),
            "q_p": self.q_p,
            "pointwise_target_time": times.tolist(),
            "pointwise_target_time_min": float(np.min(times)),
            "pointwise_target_time_max": float(np.max(times)),
            "pointwise_target_time_spread": float(np.max(times) - np.min(times)),
            "pointwise_target_residual_max_abs": float(np.max(np.abs(residual))),
            "source_ideal_reference_time": float(source_ideal_time),
            "finite_probe_only": True,
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
    def load_configuration(cls, path: str | Path) -> "KokunoCurrentExteriorLMinusHBridge":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        candidate = cls()
        if payload != candidate.configuration():
            raise ValueError("current l=-h bridge configuration/provenance mismatch")
        return candidate
