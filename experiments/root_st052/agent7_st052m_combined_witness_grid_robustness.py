"""Audit grid robustness of the razor-thin morphology miss from Agent-7 PR #652.

Issue #658 freezes this diagnostic before evaluation.  The exact nonlinear child
from #652 is not retuned or promoted.  We only ask whether the `41^3` smooth-tip
Pareto miss is numerically resolved relative to the grid dependence of the same
threshold-free vorticity observable.

This is expression-capacity / numerical-observable robustness evidence only.
It is not PDE validation, visual-correspondence verification, or identification
of an exact OpenAI field.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_combined_pareto_witness as witness

TASK_ID = "CR003-ST052M-COMBINED-WITNESS-GRID-ROBUSTNESS-093"
PREREG_ISSUE = 658
SOURCE_WITNESS_PR = 652
SOURCE_WITNESS_HEAD = "b0d98d0d5fdb3ca6a438fdf6ac0dbaaef7781faf"
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = witness.SOURCE_TEMPORAL_HEAD
GRID_RESOLUTIONS = (33, 41, 49)
MID_TIME = 0.50
PARETO_TOL = 5.0e-5
MARGIN_RESOLUTION_FRACTION = 0.25

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "witness_retuned": False,
    "parameter_grid_scan_performed": False,
    "optimization_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "parent_pde_receipt_transferred": False,
    "public_image_numeric_target_used": False,
    "visual_acceptance_threshold_defined": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "source_correspondence_verified": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def morphology_metrics(velocity_fn, time: float, resolution: int) -> dict:
    axis, velocity = witness.base.morph.prior.sample_velocity_grid(
        velocity_fn, float(time), int(resolution)
    )
    spacing = float(axis[1] - axis[0])
    _, omega = witness.base.morph.prior.vorticity(velocity, spacing)
    return witness.base.morph.enstrophy_moment_metrics(omega, axis)


def compare_metrics(control: dict, candidate: dict) -> dict:
    return witness.base.morph.compare_metrics(control, candidate)


def resolution_record(control_fn, ramp_fn, child_fn, resolution: int) -> dict:
    control_abs = morphology_metrics(control_fn, MID_TIME, resolution)
    linear_abs = morphology_metrics(ramp_fn, MID_TIME, resolution)
    child_abs = morphology_metrics(child_fn, MID_TIME, resolution)
    linear_rel = compare_metrics(control_abs, linear_abs)
    child_rel = compare_metrics(control_abs, child_abs)
    return {
        "resolution": int(resolution),
        "absolute": {
            "control": control_abs,
            "linear": linear_abs,
            "child": child_abs,
        },
        "relative_to_control": {
            "linear": linear_rel,
            "child": child_rel,
        },
        "increments_child_vs_linear": {
            "aspect_gain_increment": float(
                child_rel["full_aspect_ratio"] - linear_rel["full_aspect_ratio"]
            ),
            "tip_thinning_increment": float(
                linear_rel["smooth_tip_radial_rms"]
                - child_rel["smooth_tip_radial_rms"]
            ),
            "central_radial_rms_increment": float(
                child_rel["central_radial_rms"] - linear_rel["central_radial_rms"]
            ),
        },
    }


def _all_finite(obj) -> bool:
    if isinstance(obj, dict):
        return all(_all_finite(v) for v in obj.values())
    if isinstance(obj, (list, tuple)):
        return all(_all_finite(v) for v in obj)
    if isinstance(obj, (float, int, np.floating, np.integer)):
        return bool(np.isfinite(float(obj)))
    return True


def robustness_decision(records: dict[int, dict]) -> dict:
    tips = {
        n: float(records[n]["increments_child_vs_linear"]["tip_thinning_increment"])
        for n in GRID_RESOLUTIONS
    }
    aspects = {
        n: float(records[n]["increments_child_vs_linear"]["aspect_gain_increment"])
        for n in GRID_RESOLUTIONS
    }
    margin_41 = abs(tips[41] + PARETO_TOL)
    drift_41_to_49 = abs(tips[49] - tips[41])
    same_aspect_sign = bool(
        all(v > 0.0 for v in aspects.values())
        or all(v < 0.0 for v in aspects.values())
    )
    finite = _all_finite(records)
    both_fine_below_boundary = bool(tips[41] < -PARETO_TOL and tips[49] < -PARETO_TOL)
    drift_resolves_margin = bool(
        drift_41_to_49 <= MARGIN_RESOLUTION_FRACTION * margin_41
    )
    resolved = bool(
        finite
        and both_fine_below_boundary
        and drift_resolves_margin
        and same_aspect_sign
    )
    return {
        "frozen_pareto_tolerance": PARETO_TOL,
        "tip_increment_by_resolution": tips,
        "aspect_increment_by_resolution": aspects,
        "distance_of_41_tip_from_boundary": margin_41,
        "tip_drift_41_to_49": drift_41_to_49,
        "allowed_drift_fraction_of_41_margin": MARGIN_RESOLUTION_FRACTION,
        "allowed_drift_for_resolution_claim": MARGIN_RESOLUTION_FRACTION * margin_41,
        "both_41_and_49_below_boundary": both_fine_below_boundary,
        "drift_small_enough_to_resolve_41_margin": drift_resolves_margin,
        "aspect_increment_same_sign_all_grids": same_aspect_sign,
        "all_morphology_values_finite": finite,
        "pareto_rejection_grid_resolved": resolved,
    }


def run(out: Path) -> dict:
    # Re-run the exact #652 evaluator once as a reference identity check.  Its
    # negative Pareto verdict remains binding regardless of this grid audit.
    reference_path = Path(out).with_name("source_witness_reference.json")
    source_reference = witness.run(reference_path)
    if source_reference["combined_child_target_free_pareto_verified"] is not False:
        raise RuntimeError("source #652 negative Pareto verdict drifted")
    if abs(witness.SWIRL_A - 0.065899695471146) > 0.0:
        raise RuntimeError("#652 frozen swirl coordinate drifted")
    if abs(witness.SHOULDER_LAMBDA - (-0.02)) > 1.0e-15:
        raise RuntimeError("#652 frozen shoulder coordinate drifted")

    field, raw = witness.base.morph.prior.replay_st052.reconstruct()
    redistribution_scale, _, _ = witness.base.morph.prior.base.child_scale(field, raw)
    energy_solve = witness.base.morph.prior.comp.solve_energy_beta(
        field, raw, redistribution_scale
    )
    if not energy_solve["root_exists"]:
        raise RuntimeError("#559 frozen energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - witness.base.EXPECTED_BETA) > 2.0e-10:
        raise RuntimeError("#559 beta identity drift")

    control_fn = lambda p, t: witness.base.control_velocity(
        field, raw, p, t, redistribution_scale
    )
    ramp_fn = lambda p, t: witness.base.ramp_velocity(
        field, raw, p, t, redistribution_scale, beta
    )
    epsilon, normalization = witness.swirl.derive_epsilon(ramp_fn)
    child_fn = lambda p, t: witness.combined_velocity(
        field, raw, p, t, redistribution_scale, beta, epsilon
    )

    records = {
        n: resolution_record(control_fn, ramp_fn, child_fn, n)
        for n in GRID_RESOLUTIONS
    }
    decision = robustness_decision(records)

    # Bind the 41^3 record back to the exact source workflow result.  This is
    # numerical identity evidence only, not a new acceptance gate.
    ref_inc = source_reference["nonlinear_oriented_desirability_increments_vs_linear"]
    record_41 = records[41]["increments_child_vs_linear"]
    reference_replay = {
        "source_pareto_verified": False,
        "source_tip_increment": float(ref_inc["tip_thinning_increment"]),
        "audit_tip_increment_41": float(record_41["tip_thinning_increment"]),
        "tip_increment_abs_difference": float(
            abs(record_41["tip_thinning_increment"] - ref_inc["tip_thinning_increment"])
        ),
        "source_aspect_increment": float(ref_inc["aspect_gain_increment"]),
        "audit_aspect_increment_41": float(record_41["aspect_gain_increment"]),
        "aspect_increment_abs_difference": float(
            abs(record_41["aspect_gain_increment"] - ref_inc["aspect_gain_increment"])
        ),
    }

    if decision["pareto_rejection_grid_resolved"]:
        routing = (
            "retain #652 rejection as grid-resolved morphology/trajectory tradeoff; "
            "keep #587 minimal family and do not grow basis without a new target-specific visual defect"
        )
    else:
        routing = (
            "retain #652 formal rejection, but its tolerance-edge tip miss is not grid-resolved; "
            "do not add/retune basis from this margin; keep #587 minimal family and route future "
            "selection to a fixed target-specific visual objective"
        )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_witness": {"pr": SOURCE_WITNESS_PR, "head": SOURCE_WITNESS_HEAD},
        "source_temporal_parent": {"pr": SOURCE_TEMPORAL_PR, "head": SOURCE_TEMPORAL_HEAD},
        "frozen_child": {
            "swirl_a": witness.SWIRL_A,
            "shoulder_b": witness.SHOULDER_B,
            "shoulder_lambda": witness.SHOULDER_LAMBDA,
            "compact_swirl_epsilon": float(epsilon),
            "effective_compact_swirl_amplitude": float(witness.SWIRL_A * epsilon),
            "time_envelope": "g(t)=2*(t-.25)",
            "parameter_grid_scan_performed": False,
            "optimization_performed": False,
        },
        "grid_resolutions": list(GRID_RESOLUTIONS),
        "time": MID_TIME,
        "normalization": normalization,
        "morphology_by_resolution": {str(k): v for k, v in records.items()},
        "source_41_replay_identity": reference_replay,
        "grid_robustness_decision": decision,
        "source_652_formal_rejection_remains_binding": True,
        "combined_child_promoted": False,
        "basis_growth_justified_by_this_audit": False,
        "routing": routing,
        **TRUTH,
    }

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    r = run(args.out)
    print("source_652_formal_rejection_remains_binding=", r["source_652_formal_rejection_remains_binding"])
    print("grid_robustness_decision=", r["grid_robustness_decision"])
    print("source_41_replay_identity=", r["source_41_replay_identity"])
    print("routing=", r["routing"])


if __name__ == "__main__":
    main()
