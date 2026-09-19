"""Multi-step finite correction-cycle ledger for Kokuno Agent 3.

This module extends the single-transition actual-defect audit in
``kokuno_finite_correction_cycle_contract`` to a finite sequence

    u^(k+1) = u^k + delta_u_k.

Every adjacent transition is re-audited from raw velocity, velocity-time-
derivative, pressure, and restricted-forcing providers.  The ledger never
accepts a caller-supplied residual, defect, stress, target, gain, normalized
score, or scientific threshold.  It preserves one disjoint held-in/held-out
split across the finite cycle and records empirical contraction factors from
actual recomputed defects.

The empirical factors are an engineering divergence/rejection diagnostic, not
Kokuno's formal correction-cycle gain theorem.  The deterministic receipt uses
an analytic nonzero mechanics regression only; it is not evidence that a real
Kokuno composite candidate has run a correction cycle or reached the final
1e-3 Navier--Stokes gate.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from .kokuno_finite_correction_cycle_contract import (
    CorrectionFieldProvider,
    FiniteCorrectionCycleReport,
    audit_finite_correction_transition,
)
from .kokuno_same_cycle_defect_contract import (
    CycleIdentity,
    RestrictedForcingProvider,
    ScalarFieldProvider,
    VectorFieldProvider,
)


@dataclass(frozen=True)
class CycleStateProviders:
    """Raw providers for one assembled finite-correction cycle state."""

    velocity: VectorFieldProvider
    velocity_dt: VectorFieldProvider
    pressure: ScalarFieldProvider
    restricted_forcing: RestrictedForcingProvider

    @property
    def identity(self) -> CycleIdentity:
        identity = self.velocity.identity
        for name, other in (
            ("velocity_dt", self.velocity_dt.identity),
            ("pressure", self.pressure.identity),
            ("restricted_forcing", self.restricted_forcing.identity),
        ):
            if other != identity:
                raise ValueError(f"state {name} does not share one cycle identity")
        return identity


@dataclass(frozen=True)
class FiniteCorrectionCycleLedgerReport:
    """Audited finite sequence and empirical actual-defect contraction ledger."""

    cycle_id: str
    state_count: int
    step_count: int
    step_reports: tuple[FiniteCorrectionCycleReport, ...]
    held_in_rms_contraction_factors: tuple[float, ...]
    held_out_rms_contraction_factors: tuple[float, ...]
    held_in_max_contraction_factors: tuple[float, ...]
    held_out_max_contraction_factors: tuple[float, ...]
    cumulative_held_in_rms_ratio: float
    cumulative_held_out_rms_ratio: float
    cumulative_held_in_max_ratio: float
    cumulative_held_out_max_ratio: float
    worst_step_held_out_rms_ratio: float
    cumulative_correction_rms_budget: float
    maximum_correction_rms: float
    maximum_correction_divergence_max: float
    every_step_nontrivial: bool
    every_step_actual_defect_gain_guard_passed: bool
    empirical_cycle_contraction_passed: bool

    def to_receipt(self) -> dict[str, object]:
        payload = asdict(self)
        payload["step_reports"] = [step.to_receipt() for step in self.step_reports]
        return payload


def _as_points(
    points: Sequence[Sequence[float]] | np.ndarray, *, label: str
) -> np.ndarray:
    arr = np.asarray(points, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 4 or arr.shape[0] == 0:
        raise ValueError(f"{label} must have shape (N, 4) with N >= 1")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{label} must contain only finite values")
    rows = [tuple(float(v) for v in row) for row in arr]
    if len(set(rows)) != len(rows):
        raise ValueError(f"{label} contains duplicate points")
    return arr


def _ratio(final: float, initial: float) -> float:
    if initial == 0.0:
        return 1.0 if final == 0.0 else float("inf")
    return float(final / initial)


def _validate_cycle_shape(
    states: Sequence[CycleStateProviders],
    corrections: Sequence[CorrectionFieldProvider],
) -> tuple[CycleStateProviders, ...]:
    state_tuple = tuple(states)
    correction_tuple = tuple(corrections)
    if len(state_tuple) < 2:
        raise ValueError("finite cycle requires at least two states")
    if len(correction_tuple) != len(state_tuple) - 1:
        raise ValueError("finite cycle requires exactly one correction per transition")

    identities = [state.identity for state in state_tuple]
    cycle_id = identities[0].cycle_id
    for index, identity in enumerate(identities):
        if identity.cycle_id != cycle_id:
            raise ValueError("all finite-cycle states must share one cycle_id")
        if identity.cycle_index != identities[0].cycle_index + index:
            raise ValueError("finite-cycle state indices must be consecutive")
        if index > 0 and identity.state_token == identities[index - 1].state_token:
            raise ValueError("adjacent finite-cycle states must change state_token")

    for index, correction in enumerate(correction_tuple):
        if correction.from_identity != identities[index]:
            raise ValueError("correction does not bind the exact from-state identity")
        if correction.to_identity != identities[index + 1]:
            raise ValueError("correction does not bind the exact to-state identity")
    return state_tuple


def audit_finite_correction_cycle(
    states: Sequence[CycleStateProviders],
    corrections: Sequence[CorrectionFieldProvider],
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    update_check_points: Sequence[Sequence[float]] | np.ndarray,
    *,
    viscosity: float,
    spatial_step: float,
) -> FiniteCorrectionCycleLedgerReport:
    """Recompute every adjacent actual-defect transition in a finite cycle."""

    state_tuple = _validate_cycle_shape(states, corrections)
    correction_tuple = tuple(corrections)
    held_in = _as_points(held_in_points, label="held_in_points")
    held_out = _as_points(held_out_points, label="held_out_points")
    update_points = _as_points(update_check_points, label="update_check_points")
    held_in_rows = {tuple(float(v) for v in row) for row in held_in}
    held_out_rows = {tuple(float(v) for v in row) for row in held_out}
    if held_in_rows & held_out_rows:
        raise ValueError("held-in and held-out point sets must be disjoint")

    step_reports: list[FiniteCorrectionCycleReport] = []
    for index, correction in enumerate(correction_tuple):
        before = state_tuple[index]
        after = state_tuple[index + 1]
        step_reports.append(
            audit_finite_correction_transition(
                before.velocity,
                before.velocity_dt,
                before.pressure,
                before.restricted_forcing,
                after.velocity,
                after.velocity_dt,
                after.pressure,
                after.restricted_forcing,
                correction,
                held_in,
                held_out,
                update_points,
                viscosity=viscosity,
                spatial_step=spatial_step,
            )
        )

    first = step_reports[0]
    last = step_reports[-1]
    held_in_rms_factors = tuple(step.held_in_rms_ratio for step in step_reports)
    held_out_rms_factors = tuple(step.held_out_rms_ratio for step in step_reports)
    held_in_max_factors = tuple(step.held_in_max_ratio for step in step_reports)
    held_out_max_factors = tuple(step.held_out_max_ratio for step in step_reports)
    every_step_gain = all(step.finite_step_gain_guard_passed for step in step_reports)
    every_step_nontrivial = all(step.correction_nontrivial for step in step_reports)

    return FiniteCorrectionCycleLedgerReport(
        cycle_id=state_tuple[0].identity.cycle_id,
        state_count=len(state_tuple),
        step_count=len(step_reports),
        step_reports=tuple(step_reports),
        held_in_rms_contraction_factors=held_in_rms_factors,
        held_out_rms_contraction_factors=held_out_rms_factors,
        held_in_max_contraction_factors=held_in_max_factors,
        held_out_max_contraction_factors=held_out_max_factors,
        cumulative_held_in_rms_ratio=_ratio(
            last.held_in_after_rms, first.held_in_before_rms
        ),
        cumulative_held_out_rms_ratio=_ratio(
            last.held_out_after_rms, first.held_out_before_rms
        ),
        cumulative_held_in_max_ratio=_ratio(
            last.held_in_after_max, first.held_in_before_max
        ),
        cumulative_held_out_max_ratio=_ratio(
            last.held_out_after_max, first.held_out_before_max
        ),
        worst_step_held_out_rms_ratio=max(held_out_rms_factors),
        cumulative_correction_rms_budget=float(
            sum(step.correction_vector_rms for step in step_reports)
        ),
        maximum_correction_rms=max(step.correction_vector_rms for step in step_reports),
        maximum_correction_divergence_max=max(
            step.correction_divergence_max for step in step_reports
        ),
        every_step_nontrivial=every_step_nontrivial,
        every_step_actual_defect_gain_guard_passed=every_step_gain,
        empirical_cycle_contraction_passed=bool(
            every_step_gain
            and every_step_nontrivial
            and last.held_in_after_rms < first.held_in_before_rms
            and last.held_out_after_rms < first.held_out_before_rms
        ),
    )


def admit_finite_correction_cycle(*args: object, **kwargs: object) -> FiniteCorrectionCycleLedgerReport:
    """Fail closed unless every recomputed actual-defect transition contracts."""

    report = audit_finite_correction_cycle(*args, **kwargs)
    if not report.every_step_actual_defect_gain_guard_passed:
        raise ValueError("at least one finite correction step failed the actual-defect gain guard")
    if not report.empirical_cycle_contraction_passed:
        raise ValueError("finite correction cycle failed the empirical contraction guard")
    return report


def truth_boundary() -> dict[str, object]:
    return {
        "finite_correction_cycle_ledger_executable": True,
        "each_step_actual_defect_recomputed_from_raw_providers": True,
        "one_disjoint_heldin_heldout_split_preserved_across_cycle": True,
        "empirical_actual_defect_contraction_diagnostic_executable": True,
        "divergent_intermediate_step_rejected": True,
        "caller_supplied_residual_allowed": False,
        "caller_supplied_defect_allowed": False,
        "caller_supplied_stress_allowed": False,
        "caller_supplied_target_allowed": False,
        "caller_supplied_gain_allowed": False,
        "caller_supplied_normalized_score_allowed": False,
        "caller_supplied_scientific_threshold_allowed": False,
        "formal_kokuno_correction_cycle_gain_bound_claimed": False,
        "restricted_forcing_semantics_independently_validated_here": False,
        "real_full_candidate_defect_consumed": False,
        "full_same_cycle_composite_requested_stress_materialized": False,
        "candidate_finite_head_mean_debt_materialized": False,
        "signed_mean_inverse_input_ready": False,
        "public_velocity_correction_materialized": False,
        "real_candidate_finite_correction_cycle_run": False,
        "heldout_normalized_ns_residual_assessed": False,
        "residual_reduction_claimed_for_real_candidate": False,
        "pde_validated": False,
        "final_normalized_momentum_gate": 1.0e-3,
        "final_normalized_divergence_gate": 1.0e-5,
    }


def _analytic_state(scale: float, index: int) -> CycleStateProviders:
    identity = CycleIdentity("analytic-finite-cycle-ledger-v1", index, f"scale-{scale:g}")
    velocity = VectorFieldProvider(
        identity,
        f"analytic:u={scale:g}*t*y e_x",
        lambda x, y, z, t, s=scale: (s * t * y, 0.0, 0.0),
    )
    velocity_dt = VectorFieldProvider(
        identity,
        f"analytic:u_t={scale:g}*y e_x",
        lambda x, y, z, t, s=scale: (s * y, 0.0, 0.0),
    )
    pressure = ScalarFieldProvider(
        identity,
        "analytic:p=0",
        lambda x, y, z, t: 0.0,
    )
    forcing = RestrictedForcingProvider(
        identity,
        "analytic:fixed-zero-forcing",
        "regression-zero-forcing; fixed before residual evaluation",
        lambda x, y, z, t: (0.0, 0.0, 0.0),
    )
    return CycleStateProviders(velocity, velocity_dt, pressure, forcing)


def _analytic_cycle(scales: Sequence[float]) -> tuple[tuple[CycleStateProviders, ...], tuple[CorrectionFieldProvider, ...]]:
    if len(scales) < 2:
        raise ValueError("analytic cycle needs at least two scales")
    states = tuple(_analytic_state(float(scale), index) for index, scale in enumerate(scales))
    corrections: list[CorrectionFieldProvider] = []
    for index, (before_scale, after_scale) in enumerate(zip(scales[:-1], scales[1:])):
        delta = float(after_scale - before_scale)
        corrections.append(
            CorrectionFieldProvider(
                states[index].identity,
                states[index + 1].identity,
                f"analytic:delta_scale={delta:g}",
                lambda x, y, z, t, d=delta: (d * t * y, 0.0, 0.0),
                lambda x, y, z, t, d=delta: (d * y, 0.0, 0.0),
            )
        )
    return states, tuple(corrections)


def _regression_points() -> tuple[list[tuple[float, float, float, float]], list[tuple[float, float, float, float]], list[tuple[float, float, float, float]]]:
    held_in = [
        (-0.31, 0.27, 0.14, 0.37),
        (0.22, -0.19, -0.33, 0.63),
        (0.17, 0.41, -0.12, 0.46),
    ]
    held_out = [
        (0.11, 0.36, -0.21, 0.41),
        (-0.28, -0.17, 0.29, 0.59),
        (0.33, 0.23, 0.18, 0.54),
    ]
    update = [
        (-0.13, 0.29, 0.08, 0.35),
        (0.26, -0.24, -0.15, 0.51),
        (0.19, 0.31, 0.27, 0.67),
    ]
    return held_in, held_out, update


def deterministic_receipt() -> dict[str, object]:
    held_in, held_out, update = _regression_points()
    states, corrections = _analytic_cycle((1.0, 0.5, 0.25))
    report = admit_finite_correction_cycle(
        states,
        corrections,
        held_in,
        held_out,
        update,
        viscosity=0.01,
        spatial_step=1.0e-3,
    )

    bad_states, bad_corrections = _analytic_cycle((1.0, 0.5, 0.75))
    expansion_rejected = False
    try:
        admit_finite_correction_cycle(
            bad_states,
            bad_corrections,
            held_in,
            held_out,
            update,
            viscosity=0.01,
            spatial_step=1.0e-3,
        )
    except ValueError:
        expansion_rejected = True

    return {
        "schema": "kokuno-a3-finite-correction-cycle-ledger-v1",
        "task": "KOKUNO-A3-FINITE-CYCLE-LEDGER-048",
        "provenance": {
            "parent_agent3_pr": 637,
            "parent_agent3_head": "f7830a34784127a9ad7894047039b0ea1e8795eb",
            "purpose": "multi-step actual-defect finite-cycle ledger while full Agent-1 composite remains unavailable",
        },
        "analytic_mechanics_regression": {
            "scales": [1.0, 0.5, 0.25],
            "expected_step_rms_ratios": [0.5, 0.5],
            "expected_cumulative_rms_ratio": 0.25,
            "report": report.to_receipt(),
            "expansion_mutation_scales": [1.0, 0.5, 0.75],
            "expansion_mutation_rejected": expansion_rejected,
        },
        "truth_boundary": truth_boundary(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = deterministic_receipt()
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
