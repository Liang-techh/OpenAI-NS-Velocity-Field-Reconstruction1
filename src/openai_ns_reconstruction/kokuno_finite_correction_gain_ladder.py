"""Derivative-ladder robustness audit for one Kokuno Agent-3 correction step.

The correction direction and damping factor are selected exactly once by the existing
held-in-only finite search, at the project's finest registered Cartesian derivative
step.  After that selection is frozen, the same transition is replayed on the full
project derivative ladder.  Held-out data never selects or retunes ``alpha``.

This is a repository engineering robustness diagnostic for the finite correction
cycle.  It is not Kokuno's formal gain theorem, not a volume-weighted full-domain
PDE validator, and not evidence that the real leading+oscillatory composite has
already been corrected.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from .kokuno_finite_correction_cycle_contract import (
    FiniteCorrectionCycleReport,
    audit_finite_correction_transition,
)
from .kokuno_finite_correction_cycle_ledger import CycleStateProviders
from .kokuno_finite_correction_damping import (
    CorrectionDirectionProvider,
    _analytic_before,
    _analytic_direction,
    _build_trial,
    _regression_points,
    audit_damped_correction_search,
)

PROJECT_VISCOSITY = 0.01
PROJECT_DERIVATIVE_STEPS: tuple[float, ...] = (0.02, 0.01, 0.005)
FINAL_MOMENTUM_MAX_GATE = 1.0e-3
FINAL_MOMENTUM_L2_GATE = 1.0e-3
FINAL_DIVERGENCE_MAX_GATE = 1.0e-5
FINAL_DIVERGENCE_L2_GATE = 1.0e-5


def _constraints_path() -> Path:
    return Path(__file__).resolve().parents[2] / "configs" / "constraints.json"


def _load_project_constraints() -> dict[str, object]:
    payload = json.loads(_constraints_path().read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("project constraints must decode to an object")
    return payload


def _validate_project_constraints(payload: Mapping[str, object]) -> None:
    """Fail closed if the registered derivative/gate protocol has drifted."""

    if float(payload.get("nu", float("nan"))) != PROJECT_VISCOSITY:
        raise ValueError("project viscosity drifted from the frozen correction audit")
    if payload.get("units") != "dimensionless; L_ref=U_ref=1, residual scale=U_ref^2/L_ref=1":
        raise ValueError("project residual normalization contract drifted")
    validation = payload.get("validation")
    if not isinstance(validation, Mapping):
        raise ValueError("project validation contract is missing")
    steps = tuple(float(x) for x in validation.get("derivative_steps", ()))
    if steps != PROJECT_DERIVATIVE_STEPS:
        raise ValueError("project derivative ladder drifted")
    norms = validation.get("norms")
    if norms != [
        "max Euclidean vector norm",
        "volume-weighted L2 spatial norm at each time",
    ]:
        raise ValueError("project validation norm contract drifted")
    thresholds = validation.get("thresholds")
    if not isinstance(thresholds, Mapping):
        raise ValueError("project validation thresholds are missing")
    expected = {
        "pde_residual_max": FINAL_MOMENTUM_MAX_GATE,
        "pde_residual_L2": FINAL_MOMENTUM_L2_GATE,
        "divergence_max": FINAL_DIVERGENCE_MAX_GATE,
        "divergence_L2": FINAL_DIVERGENCE_L2_GATE,
    }
    for key, value in expected.items():
        if float(thresholds.get(key, float("nan"))) != value:
            raise ValueError(f"project scientific threshold drifted: {key}")
    failure_policy = str(validation.get("failure_policy", ""))
    if "changing thresholds requires a new experiment version" not in failure_policy:
        raise ValueError("project threshold-change failure policy drifted")


@dataclass(frozen=True)
class GainLadderLevel:
    spatial_step: float
    held_in_rms_ratio: float
    held_in_max_ratio: float
    held_out_rms_ratio: float
    held_out_max_ratio: float
    correction_divergence_rms: float
    correction_divergence_max: float
    finite_step_gain_guard_passed: bool


@dataclass(frozen=True)
class CorrectionGainLadderReport:
    cycle_id: str
    from_cycle_index: int
    selected_alpha: float
    alpha_selected_at_spatial_step: float
    project_derivative_steps: tuple[float, ...]
    levels: tuple[GainLadderLevel, ...]
    heldout_used_for_alpha_selection: bool
    alpha_retuned_across_derivative_steps: bool
    all_derivative_levels_gain_guard_passed: bool
    maximum_held_out_rms_ratio: float
    maximum_held_out_max_ratio: float

    def to_receipt(self) -> dict[str, object]:
        payload = asdict(self)
        payload["levels"] = [asdict(level) for level in self.levels]
        return payload


def _level_from_transition(
    spatial_step: float, transition: FiniteCorrectionCycleReport
) -> GainLadderLevel:
    return GainLadderLevel(
        spatial_step=float(spatial_step),
        held_in_rms_ratio=transition.held_in_rms_ratio,
        held_in_max_ratio=transition.held_in_max_ratio,
        held_out_rms_ratio=transition.held_out_rms_ratio,
        held_out_max_ratio=transition.held_out_max_ratio,
        correction_divergence_rms=transition.correction_divergence_rms,
        correction_divergence_max=transition.correction_divergence_max,
        finite_step_gain_guard_passed=bool(transition.finite_step_gain_guard_passed),
    )


def audit_project_derivative_ladder_gain(
    before: CycleStateProviders,
    direction: CorrectionDirectionProvider,
    held_in_points: Sequence[Sequence[float]] | np.ndarray,
    held_out_points: Sequence[Sequence[float]] | np.ndarray,
    update_check_points: Sequence[Sequence[float]] | np.ndarray,
) -> CorrectionGainLadderReport:
    """Freeze alpha on held-in finest-step data, then replay one step on all h."""

    _validate_project_constraints(_load_project_constraints())
    finest_step = PROJECT_DERIVATIVE_STEPS[-1]
    search = audit_damped_correction_search(
        before,
        direction,
        held_in_points,
        held_out_points,
        update_check_points,
        viscosity=PROJECT_VISCOSITY,
        spatial_step=finest_step,
    )
    selected_alpha = float(search.selected_alpha)
    after, correction = _build_trial(before, direction, selected_alpha)

    levels: list[GainLadderLevel] = []
    for spatial_step in PROJECT_DERIVATIVE_STEPS:
        transition = audit_finite_correction_transition(
            before.velocity,
            before.velocity_dt,
            before.pressure,
            before.restricted_forcing,
            after.velocity,
            after.velocity_dt,
            after.pressure,
            after.restricted_forcing,
            correction,
            held_in_points,
            held_out_points,
            update_check_points,
            viscosity=PROJECT_VISCOSITY,
            spatial_step=spatial_step,
        )
        levels.append(_level_from_transition(spatial_step, transition))

    all_passed = all(level.finite_step_gain_guard_passed for level in levels)
    return CorrectionGainLadderReport(
        cycle_id=before.identity.cycle_id,
        from_cycle_index=before.identity.cycle_index,
        selected_alpha=selected_alpha,
        alpha_selected_at_spatial_step=finest_step,
        project_derivative_steps=PROJECT_DERIVATIVE_STEPS,
        levels=tuple(levels),
        heldout_used_for_alpha_selection=bool(search.heldout_used_for_alpha_selection),
        alpha_retuned_across_derivative_steps=False,
        all_derivative_levels_gain_guard_passed=bool(all_passed),
        maximum_held_out_rms_ratio=max(level.held_out_rms_ratio for level in levels),
        maximum_held_out_max_ratio=max(level.held_out_max_ratio for level in levels),
    )


def admit_project_derivative_ladder_gain(*args: object, **kwargs: object) -> CorrectionGainLadderReport:
    """Reject a correction if its actual-defect gain fails at any registered h."""

    report = audit_project_derivative_ladder_gain(*args, **kwargs)
    if not report.all_derivative_levels_gain_guard_passed:
        raise ValueError("selected correction fails the frozen derivative-ladder gain audit")
    return report


def truth_boundary() -> dict[str, object]:
    return {
        "project_derivative_ladder_gain_audit_executable": True,
        "project_viscosity": PROJECT_VISCOSITY,
        "project_derivative_steps": list(PROJECT_DERIVATIVE_STEPS),
        "alpha_selected_at_finest_registered_step": True,
        "heldout_used_for_alpha_selection": False,
        "alpha_retuned_across_derivative_steps": False,
        "actual_defect_recomputed_from_raw_providers_at_each_step": True,
        "project_thresholds_fail_closed_against_config_drift": True,
        "sampled_rms_relabelled_as_volume_weighted_project_l2": False,
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
        "final_normalized_momentum_max_gate": FINAL_MOMENTUM_MAX_GATE,
        "final_normalized_momentum_l2_gate": FINAL_MOMENTUM_L2_GATE,
        "final_normalized_divergence_max_gate": FINAL_DIVERGENCE_MAX_GATE,
        "final_normalized_divergence_l2_gate": FINAL_DIVERGENCE_L2_GATE,
    }


def deterministic_receipt() -> dict[str, object]:
    held_in, held_out, update = _regression_points()
    report = audit_project_derivative_ladder_gain(
        _analytic_before(), _analytic_direction(), held_in, held_out, update
    )
    return {
        "schema": "kokuno-a3-finite-correction-gain-ladder-v1",
        "task": "KOKUNO-A3-FINITE-CORRECTION-GAIN-LADDER-050",
        "provenance": {
            "parent_agent3_pr": 650,
            "parent_agent3_head": "2a743997d4883967c6b7930e150df773c86d922e",
            "source_classification": (
                "Kokuno corrected finite-correction/gain structure motivates rejection; "
                "the derivative ladder and held-in-only damping policy are repository engineering"
            ),
        },
        "analytic_mechanics_regression": {
            "description": (
                "u=(t*y,0,0), fixed f=0, oversized delta scale -2.4; "
                "selected alpha=0.5 gives exact defect scale 0.2"
            ),
            "report": report.to_receipt(),
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
