"""Current-candidate Kokuno exterior ``Q_s`` state at the release2 endpoint.

Pinned public provenance
------------------------
KokunoYumeto/yang-mills-interacting-workbench
commit 143f6773feb424ad9ed3a8d116653200f20346b7
``navier-stokes/navier_stokes_workbench.tex``
blob 205a99807302e21a51c5eaf223390c0dfc42bcd0
corrected 2026-09-09 reconstruction.

The corrected reconstruction gives the general regular radial identity

    W = 1 - 2 D eta M/X - d M_eta/X,

    Q_s = -W + ((1-h) I - D eta I_eta - d J_eta
                + 2(h-D) eta J)/(X H).

The source later simplifies this only after imposing its ideal state
``I=XH/(1-lambda), U=M=J=0``.  A1 #1217 deliberately kept that source-ideal
schedule separate because the current repository candidate preserves its actual
post-pulse ``M/M_eta`` history.  This module closes exactly that bookkeeping
seam at the end of the public ``l:-1 -> -h`` release.  It does *not* yet build
the subsequent ``l=-h`` matching bridge.

The current angular state comes from A1 #1168/#1171.  The relative-swirl solve
sets the endpoint angular condition up to its actual algebraic row residual;
we transport that residual stably rather than subtracting two O(1) values.
The current ``J,J_eta`` state is reconstructed from the already-executable
A1 #1133 pulse-end row, including the real frozen-I1 ``J`` entry and analytic
eta jets.  Since ``U=0`` after the pulse, physical ``J`` remains constant
through eta flattening, relative swirl, release1, the ``l=-1`` hold, and
release2.  ``M/M_eta`` are read directly from the exact #1204 release2 endpoint.

This is executable current-candidate state, not a paper-exact field and not
Navier--Stokes validation.  No pressure, forcing, residual, optimizer or
hidden/source parameter is introduced here.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .kokuno_pa16_current_cartesian_pulse_end_compensated import (
    KokunoPA16CurrentCartesianPulseEndCompensated,
)
from .kokuno_public_ideal_exterior_qs_schedule import (
    KokunoPublicIdealExteriorQsSchedule,
)

SCHEMA = "kokuno-current-exterior-qs-release2-state-v1"
PARENT_EXACT_HEAD = "3978560078104bff7c9d5556bd9e658e4874c963"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_BLOB = "205a99807302e21a51c5eaf223390c0dfc42bcd0"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_RELEASE = "corrected-208-page-reconstruction"
SOURCE_RELEASE_DATE = "2026-09-09"
_LOG2 = math.log(2.0)

_SOURCE_FORMULAS = {
    "radial_identity": (
        "W=1-2D*eta*M/X-d*M_eta/X; "
        "Q_s=-W+((1-h)I-D*eta*I_eta-d*J_eta+2(h-D)*eta*J)/(XH)"
    ),
    "angular_cumulative": "H=sqrt(2X)E; I=int_0^X H dx; r_I=I/(XH)",
    "angular_transport": "D_logX r_I+(1+l)r_I=1",
    "pulse_J": "J=int_0^X U H dx",
    "postpulse": "U=0 after the pulse; physical M,M_eta,J,J_eta remain constant",
    "general_q_transport": (
        "D_logX Q_s+(1+l)Q_s=-W*l-h(1-2etaU)-H_c*(log E)_eta"
    ),
    "source_ideal_specialization": (
        "I=XH/(1-lambda), U=M=J=0 => Q_s=(lambda-h)/(1-lambda) at exterior entry"
    ),
}

_NUMERICAL_REALIZATION = {
    "parent": "consume exact A1 #1217 source-ideal schedule and exact #1204 current release2 identity",
    "angular_state": (
        "consume exact #1168 current angular target and #1154 row solve; retain its actual row residual"
    ),
    "pulse_J_state": (
        "reconstruct exact #1133 pulse-end J/J_eta normalization from its public two-row algebra and jets"
    ),
    "current_primitives": "read exact #1204 M/X and M_eta/X at the release2 endpoint",
    "arithmetic": "binary64 state algebra; no residual-driven or held-out tuning",
    "new_residual_tuning_parameters": "none",
}

_TRUTH_UPDATES = {
    "current_q_s_release2_endpoint_materialized": True,
    "current_q_s_uses_general_radial_identity": True,
    "current_absolute_I_state_consumed": True,
    "current_actual_M_M_eta_state_consumed": True,
    "current_actual_J_J_eta_state_consumed": True,
    "source_ideal_q_s_relabelled_as_current": False,
    "current_l_minus_h_matching_bridge_materialized": False,
    "current_cartesian_terminal_multiplier_composed": False,
    "source_terminal_multiplier_materialized": False,
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


def _finite_eta(value: Any, interval: tuple[float, float]) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if np.any(~np.isfinite(out)):
        raise ValueError("eta must contain only finite values")
    lo, hi = interval
    if np.any((out < lo) | (out > hi)):
        raise ValueError(f"eta must lie in [{lo},{hi}]")
    return out


@dataclass(frozen=True)
class KokunoCurrentExteriorQsRelease2State:
    """Exact-current bookkeeping for ``Q_s`` at the #1204 release2 endpoint."""

    parent: KokunoPublicIdealExteriorQsSchedule = field(
        default_factory=KokunoPublicIdealExteriorQsSchedule,
        repr=False,
        compare=False,
    )
    _pulse_end: KokunoPA16CurrentCartesianPulseEndCompensated = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.parent, KokunoPublicIdealExteriorQsSchedule):
            raise TypeError("parent must be KokunoPublicIdealExteriorQsSchedule")
        truth = self.parent.truth_boundary
        if not truth["public_source_ideal_exterior_qs_schedule_materialized"]:
            raise ValueError("exact A1 #1217 source-ideal schedule is required")
        if truth["current_q_s_release2_endpoint_materialized"]:
            raise ValueError("parent unexpectedly already claims current Q_s")
        rtruth = self.release2.truth_boundary
        if not rtruth["current_cartesian_l_minus1_to_minus_h_composed"]:
            raise ValueError("exact A1 #1204 current release2 candidate is required")
        if not self.split.target.truth_boundary[
            "absolute_current_r_I_at_flattening_endpoint_materialized"
        ]:
            raise ValueError("exact A1 #1168 absolute current angular state is required")

        node: Any = self.split.target.hold
        pulse = None
        for _ in range(12):
            if isinstance(node, KokunoPA16CurrentCartesianPulseEndCompensated):
                pulse = node
                break
            node = getattr(node, "parent", None)
            if node is None:
                break
        if pulse is None:
            raise ValueError("current A1 ancestry no longer exposes exact #1133 pulse-end state")
        if not pulse.truth_boundary["current_pulse_endpoint_xi13_materialized"]:
            raise ValueError("pulse ancestor does not materialize the current xi=13 endpoint")
        object.__setattr__(self, "_pulse_end", pulse)

        # At the end of the relative-swirl hold both compact bumps and their
        # eta jets must be exactly inactive.  This makes XH eta-independent on
        # the subsequent public exterior stages, so I_eta/(XH)=d_eta r_I.
        end = self.split.log_X_hold_end - self.split.log_X_flatten_end
        if float(np.asarray(self.split.compensator.beta1(end))) != 0.0:
            raise ValueError("first relative-swirl bump did not close before exterior entry")
        if float(np.asarray(self.split.compensator.beta2(end))) != 0.0:
            raise ValueError("second relative-swirl bump did not close before exterior entry")

    @property
    def release2(self):
        return self.parent.release2

    @property
    def split(self):
        return self.release2.split

    @property
    def pulse_end(self) -> KokunoPA16CurrentCartesianPulseEndCompensated:
        return self._pulse_end

    @property
    def lambda_value(self) -> float:
        return float(self.release2.lambda_value)

    @property
    def h_value(self) -> float:
        return float(self.release2.h_value)

    @property
    def D(self) -> float:
        return float(self.release2.D)

    @property
    def eta_interval(self) -> tuple[float, float]:
        return tuple(float(v) for v in self.release2.eta_interval)

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

    def _angular_relative_endpoint(self, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Return actual r_I and I_eta/(XH) after the relative-swirl solve.

        #1168 chooses a scaled target that would put r_I exactly at
        1/(1-lambda).  We retain the *actual* #1154 angular row residual and
        transport only that small quantity back to the physical endpoint.
        """
        target, target_eta = self.split.target.target_with_eta(eta)
        sol, c1_eta, c2_eta = self.split.compensator.solve_target_jet(
            target, target_eta
        )
        row = self.split.compensator.angular_row_scaled
        residual = row[0] * np.asarray(sol.c1) + row[1] * np.asarray(sol.c2) - target
        residual_eta = row[0] * np.asarray(c1_eta) + row[1] * np.asarray(c2_eta) - target_eta
        # The row is scaled at y1; undo that scaling only over the bounded
        # distance from y1 to the hold endpoint.  This avoids O(1)-O(1)
        # cancellation in r_I itself.
        k = 1.0 - self.lambda_value
        attenuation = math.exp(-k * (self.split.hold_length - self.split.compensator.y1))
        r_star = 1.0 / k
        r_entry = r_star + attenuation * residual
        r_entry_eta = attenuation * residual_eta
        return np.asarray(r_entry, dtype=float), np.asarray(r_entry_eta, dtype=float)

    def _angular_release2_endpoint(self, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Transport the actual angular-row residual through the exterior schedule."""
        r_entry, r_entry_eta = self._angular_relative_endpoint(eta)
        entry_profile = self.split.split_profile_logX(
            np.full_like(eta, self.split.log_X_hold_end, dtype=float), eta
        )
        end_profile = self.release2.profile_logX(
            np.full_like(eta, self.release2.log_X_release2_end, dtype=float), eta
        )
        log_xh_entry = (
            _LOG2
            + 2.0 * self.split.log_X_hold_end
            + np.asarray(entry_profile["log_F_base"], dtype=float)
        )
        log_xh_end = (
            _LOG2
            + 2.0 * self.release2.log_X_release2_end
            + np.asarray(end_profile["log_F"], dtype=float)
        )
        homogeneous = np.exp(log_xh_entry - log_xh_end)
        if np.any(~np.isfinite(homogeneous)) or np.any(homogeneous <= 0.0):
            raise RuntimeError("current angular residual transport became invalid")

        # For the source-ideal state, Q=-1+(1-h)r at release2 end.  #1217
        # already executes Q with the same public l schedule, so this is an
        # implementation-distinct way to recover the ideal r endpoint.
        r_ideal_end = (self.parent.q_in_source_ideal + 1.0) / (1.0 - self.h_value)
        r_star = 1.0 / (1.0 - self.lambda_value)
        r_end = r_ideal_end + (r_entry - r_star) * homogeneous
        r_end_eta = r_entry_eta * homogeneous
        return np.asarray(r_end, dtype=float), np.asarray(r_end_eta, dtype=float)

    def _pulse_j_normalized_with_eta(
        self, eta: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return J/(X_p H_p E_p), its eta jet, and the scaled-row replay."""
        inputs = self.pulse_end.parent.compensator_inputs(eta)
        comp = self.pulse_end.compensator
        sol, c1_eta, c2_eta = comp.solve_eta_jet(
            inputs["amplitude"],
            inputs["m_entry_ratio"],
            inputs["j_entry_ratio"],
            inputs["amplitude_eta"],
            inputs["m_entry_ratio_eta"],
            inputs["j_entry_ratio_eta"],
        )
        matrix = np.asarray(comp.matrix_scaled, dtype=float)
        main = np.asarray(comp.main_moments_scaled, dtype=float)
        exp_scale = math.exp(-comp.s2 * comp.y1)
        j_scaled = (
            np.asarray(inputs["j_entry_ratio"], dtype=float) * exp_scale
            + np.asarray(inputs["amplitude"], dtype=float) * main[1]
            + np.asarray(sol.c1, dtype=float) * matrix[1, 0]
            + np.asarray(sol.c2, dtype=float) * matrix[1, 1]
        )
        j_scaled_eta = (
            np.asarray(inputs["j_entry_ratio_eta"], dtype=float) * exp_scale
            + np.asarray(inputs["amplitude_eta"], dtype=float) * main[1]
            + np.asarray(c1_eta, dtype=float) * matrix[1, 0]
            + np.asarray(c2_eta, dtype=float) * matrix[1, 1]
        )
        unscale = math.exp(comp.s2 * comp.y1)
        j_norm = j_scaled * unscale
        j_norm_eta = j_scaled_eta * unscale

        # Cross-check against the exact #1133 public pulse-end row value.
        pulse_values = self.pulse_end.similarity_profile_values_logX(
            np.full_like(eta, self.pulse_end.log_X_pulse_end, dtype=float), eta
        )
        row_replay = np.asarray(pulse_values["public_J_row_scaled"], dtype=float)
        scale = np.maximum.reduce(
            [np.abs(j_scaled), np.abs(row_replay), np.full_like(j_scaled, 1.0e-300)]
        )
        if np.any(np.abs(j_scaled - row_replay) > 2.0e-11 * scale + 1.0e-300):
            raise RuntimeError("current pulse J row no longer replays exact #1133")
        return np.asarray(j_norm), np.asarray(j_norm_eta), np.asarray(j_scaled)

    def _j_release2_ratios(self, eta: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return current J/(XH), J_eta/(XH), and pulse scaled-row replay."""
        j_norm, j_norm_eta, j_scaled = self._pulse_j_normalized_with_eta(eta)
        E_p, E_p_eta, _, _ = self.pulse_end.leading._pulse_entry_state(eta)
        E_p = np.asarray(E_p, dtype=float)
        E_p_eta = np.asarray(E_p_eta, dtype=float)
        if np.any(E_p <= 0.0) or np.any(~np.isfinite(E_p)):
            raise RuntimeError("current pulse-entry E is invalid")

        end = self.release2.profile_logX(
            np.full_like(eta, self.release2.log_X_release2_end, dtype=float), eta
        )
        log_E_end = np.asarray(end["log_F"], dtype=float) + 0.5 * (
            _LOG2 + self.release2.log_X_release2_end
        )
        log_factor = (
            1.5 * (self.pulse_end.log_X_p - self.release2.log_X_release2_end)
            + 2.0 * np.log(E_p)
            - log_E_end
        )
        factor = np.exp(log_factor)
        if np.any(~np.isfinite(factor)):
            raise RuntimeError("current pulse-to-exterior J scale became non-finite")
        J_ratio = factor * j_norm
        J_eta_ratio = factor * (
            j_norm_eta + 2.0 * j_norm * (E_p_eta / E_p)
        )
        return np.asarray(J_ratio), np.asarray(J_eta_ratio), np.asarray(j_scaled)

    def state(self, eta: Any) -> dict[str, np.ndarray]:
        """Materialize the exact-current release2-end radial ``Q_s`` identity."""
        values = _finite_eta(eta, self.eta_interval)
        profile = self.release2.profile_logX(
            np.full_like(values, self.release2.log_X_release2_end, dtype=float), values
        )
        M_ratio = np.asarray(profile["M_over_X"], dtype=float)
        M_eta_ratio = np.asarray(profile["M_eta_over_X"], dtype=float)
        r_I, I_eta_ratio = self._angular_release2_endpoint(values)
        J_ratio, J_eta_ratio, j_scaled = self._j_release2_ratios(values)

        d = 1.0 - values * values
        W = 1.0 - 2.0 * self.D * values * M_ratio - d * M_eta_ratio
        term_W = -W
        term_I = (1.0 - self.h_value) * r_I
        term_I_eta = -self.D * values * I_eta_ratio
        term_J_eta = -d * J_eta_ratio
        term_J = 2.0 * (self.h_value - self.D) * values * J_ratio
        q_current = term_W + term_I + term_I_eta + term_J_eta + term_J
        q_ideal = np.full_like(q_current, self.parent.q_in_source_ideal, dtype=float)
        q_delta = q_current - q_ideal

        arrays = (
            M_ratio, M_eta_ratio, r_I, I_eta_ratio, J_ratio, J_eta_ratio,
            W, term_W, term_I, term_I_eta, term_J_eta, term_J, q_current, q_delta,
        )
        if any(np.any(~np.isfinite(a)) for a in arrays):
            raise RuntimeError("current release2-end Q_s state became non-finite")

        return {
            "eta": values,
            "M_over_X": M_ratio,
            "M_eta_over_X": M_eta_ratio,
            "I_over_XH": r_I,
            "I_eta_over_XH": I_eta_ratio,
            "J_over_XH": J_ratio,
            "J_eta_over_XH": J_eta_ratio,
            "W": W,
            "term_minus_W": term_W,
            "term_I": term_I,
            "term_I_eta": term_I_eta,
            "term_J_eta": term_J_eta,
            "term_J": term_J,
            "Q_s_current_release2_end": q_current,
            "Q_s_source_ideal_release2_end": q_ideal,
            "Q_s_current_minus_source_ideal": q_delta,
            "pulse_public_J_row_scaled_replay": j_scaled,
        }

    def report(self, eta: Any = (-0.83, -0.61, -0.37, 0.0, 0.19, 0.73)) -> dict[str, Any]:
        values = _finite_eta(eta, self.eta_interval)
        state = self.state(values)
        q = np.asarray(state["Q_s_current_release2_end"], dtype=float)
        delta = np.asarray(state["Q_s_current_minus_source_ideal"], dtype=float)
        return {
            "schema": "kokuno-agent1-current-exterior-qs-release2-report-v1",
            "semantic_sha256": self.semantic_sha256,
            "parent_semantic_sha256": _semantic_value(self.parent),
            "eta": values.tolist(),
            "q_s_current_release2_end": q.tolist(),
            "q_s_source_ideal_release2_end": float(self.parent.q_in_source_ideal),
            "q_s_current_minus_source_ideal": delta.tolist(),
            "max_abs_current_minus_source_ideal": float(np.max(np.abs(delta))),
            "min_current_q_s": float(np.min(q)),
            "max_current_q_s": float(np.max(q)),
            "q_p": float(self.parent.q_p),
            "current_l_minus_h_matching_bridge_materialized": False,
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
    def load_configuration(cls, path: str | Path) -> "KokunoCurrentExteriorQsRelease2State":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        candidate = cls()
        if payload != candidate.configuration():
            raise ValueError("current exterior Q_s configuration/provenance mismatch")
        return candidate
