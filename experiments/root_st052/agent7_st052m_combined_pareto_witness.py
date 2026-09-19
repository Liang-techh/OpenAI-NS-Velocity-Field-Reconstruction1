"""Verify one preregistered nonlinear combination of two screened ST052-M directions.

Issue #651 freezes one point before this child is evaluated:

    compact-swirl coordinate a = 0.0658997
    shoulder-timing coordinate b = -1.0  -> lambda = -0.02

The parent is the exact Agent-7 #587 linear-ramp child.  The shoulder-timing
coordinate reuses #634 and the compact pure-swirl coordinate reuses #621.  No
new basis, scan, fit, image target, pressure, forcing, or common rescaling is
introduced.  This module asks only whether the *actual nonlinear child* obeys
the same target-free Pareto rule that motivated the first-order #642 cone
routing.

Expression-capacity / morphology-trajectory evidence only.  This is not PDE
validation, visual-correspondence verification, or identification of an exact
OpenAI field.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_two_direction_pareto_cone as cone

swirl = cone.swirl
shoulder = cone.shoulder
base = shoulder.base

TASK_ID = "CR003-ST052M-COMBINED-PARETO-WITNESS-092"
PREREG_ISSUE = 651
SOURCE_CONE_PR = 642
SOURCE_CONE_HEAD = "f523c221847e1a587ec9af4ec053f6f864ce90c4"
SOURCE_SWIRL_PR = 621
SOURCE_SWIRL_HEAD = "d01ea33e4e18c0fe8f3a40baaf22b233a66fe96a"
SOURCE_SHOULDER_PR = 634
SOURCE_SHOULDER_HEAD = "9c628326cad1b8637b13016ffefb945901c4dd53"
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
SOURCE_RENDER_PR = 601
SOURCE_RENDER_HEAD = "02fb726e03268d17f3ecd1a7d503282a225d50eb"

SWIRL_A = 0.0658997
SHOULDER_B = -1.0
SHOULDER_LAMBDA = SHOULDER_B * shoulder.EPSILON
PARETO_TOL = cone.PARETO_TOL
MID_TIME = 0.50
LATE_TIME = 0.75
GRID_RESOLUTION = 41
ENERGY_DRIFT_MAX = 0.005
EARLY_IDENTITY_TOL = 1.0e-11
DIVERGENCE_MAX = 1.0e-5
EARLY_SEED = 9176511
EARLY_PROBES = 96

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "screened_combined_velocity_child_constructed": True,
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


def add_compact_swirl(phased_velocity: np.ndarray, points: np.ndarray, time: float, epsilon: float) -> np.ndarray:
    """Apply the one frozen fractional #621 compact-swirl coordinate."""
    pts = np.asarray(points, dtype=float)
    phased = np.asarray(phased_velocity, dtype=float)
    correction = (
        SWIRL_A
        * float(epsilon)
        * base.activation(float(time))
        * swirl.swirl_correction(pts)
    )
    if phased.shape != correction.shape:
        raise ValueError("phased velocity and compact-swirl correction shape mismatch")
    return phased + correction


def combined_velocity(field, raw, points, time, redistribution_scale, beta, epsilon):
    """Exact one-point nonlinear child frozen in issue #651."""
    pts = np.asarray(points, dtype=float)
    phased = shoulder.phased_velocity(
        field,
        raw,
        pts,
        float(time),
        redistribution_scale,
        beta,
        SHOULDER_LAMBDA,
    )
    return add_compact_swirl(phased, pts, float(time), epsilon)


def morphology_metrics(velocity_fn, time: float) -> dict:
    axis, velocity = base.morph.prior.sample_velocity_grid(
        velocity_fn, float(time), GRID_RESOLUTION
    )
    spacing = float(axis[1] - axis[0])
    _, omega = base.morph.prior.vorticity(velocity, spacing)
    return base.morph.enstrophy_moment_metrics(omega, axis)


def morphology_relative(control_fn, velocity_fn, time: float) -> dict:
    control = morphology_metrics(control_fn, time)
    candidate = morphology_metrics(velocity_fn, time)
    return base.morph.compare_metrics(control, candidate)


def oriented_increments(linear_morph: dict, child_morph: dict, linear_paths: dict, child_paths: dict) -> dict:
    """Four frozen desirability increments; positive means better."""
    return {
        "aspect_gain_increment": float(
            child_morph["full_aspect_ratio"] - linear_morph["full_aspect_ratio"]
        ),
        "tip_thinning_increment": float(
            linear_morph["smooth_tip_radial_rms"] - child_morph["smooth_tip_radial_rms"]
        ),
        "turns_fidelity_increment": float(
            child_paths["mean_absolute_turns_relative"]
            - linear_paths["mean_absolute_turns_relative"]
        ),
        "axial_pair_fidelity_increment": float(
            child_paths["mean_pair_axial_separation_change_relative"]
            - linear_paths["mean_pair_axial_separation_change_relative"]
        ),
    }


def nonlinear_pareto_rule(increments: dict, structural_pass: bool) -> dict:
    values = np.asarray(
        [
            increments["aspect_gain_increment"],
            increments["tip_thinning_increment"],
            increments["turns_fidelity_increment"],
            increments["axial_pair_fidelity_increment"],
        ],
        dtype=float,
    )
    nonworsening = bool(np.all(values >= -PARETO_TOL - 1.0e-12))
    improving = bool(np.any(values > PARETO_TOL))
    return {
        "tolerance": PARETO_TOL,
        "all_desirabilities_nonworsening": nonworsening,
        "strict_improvement_any": improving,
        "passes": bool(structural_pass and nonworsening and improving),
    }


def early_identity(ramp_fn, child_fn) -> float:
    rng = np.random.default_rng(EARLY_SEED)
    radius = rng.uniform(0.0, 1.9, size=EARLY_PROBES)
    angle = rng.uniform(0.0, 2.0 * np.pi, size=EARLY_PROBES)
    z = rng.uniform(-1.9, 1.9, size=EARLY_PROBES)
    pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
    return float(
        np.max(
            np.abs(
                np.asarray(child_fn(pts, 0.25), float)
                - np.asarray(ramp_fn(pts, 0.25), float)
            )
        )
    )


def run(out: Path) -> dict:
    # Reconstruct exactly the same ST052-M / redistribution / energy-only beta
    # lineage used by #587, #621 and #634.
    field, raw = base.morph.prior.replay_st052.reconstruct()
    redistribution_scale, _, _ = base.morph.prior.base.child_scale(field, raw)
    energy_solve = base.morph.prior.comp.solve_energy_beta(field, raw, redistribution_scale)
    if not energy_solve["root_exists"]:
        raise RuntimeError("#559 frozen energy root disappeared")
    beta = float(energy_solve["beta"])
    if abs(beta - base.EXPECTED_BETA) > 2.0e-10:
        raise RuntimeError("#559 beta identity drift")

    control_fn = lambda p, t: base.control_velocity(field, raw, p, t, redistribution_scale)
    ramp_fn = lambda p, t: base.ramp_velocity(field, raw, p, t, redistribution_scale, beta)

    epsilon, normalization = swirl.derive_epsilon(ramp_fn)
    child_fn = lambda p, t: combined_velocity(
        field, raw, p, t, redistribution_scale, beta, epsilon
    )

    linear_mid = morphology_relative(control_fn, ramp_fn, MID_TIME)
    child_mid = morphology_relative(control_fn, child_fn, MID_TIME)
    linear_late = morphology_relative(control_fn, ramp_fn, LATE_TIME)
    child_late = morphology_relative(control_fn, child_fn, LATE_TIME)

    paths_control = base.measure_paths(control_fn)
    paths_linear = base.measure_paths(ramp_fn)
    paths_child = base.measure_paths(child_fn)
    path_relative = {
        "linear": base.compare_paths(paths_linear, paths_control),
        "child": base.compare_paths(paths_child, paths_control),
    }

    increments = oriented_increments(
        linear_mid,
        child_mid,
        path_relative["linear"],
        path_relative["child"],
    )

    identity = early_identity(ramp_fn, child_fn)
    structure = base.structure_preflight(child_fn, seed=9176512)
    e_linear = base.axisymmetric_energy_at_time(ramp_fn, MID_TIME)
    e_child = base.axisymmetric_energy_at_time(child_fn, MID_TIME)
    energy_drift = float(e_child / e_linear - 1.0)

    structural_pass = bool(
        identity <= EARLY_IDENTITY_TOL
        and structure["support_max_abs"] <= base.CRITERIA["support_max_abs"]
        and structure["divergence_fd_max"] <= DIVERGENCE_MAX
        and abs(energy_drift) <= ENERGY_DRIFT_MAX
    )
    pareto = nonlinear_pareto_rule(increments, structural_pass)
    verified = bool(pareto["passes"])

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "sources": {
            "first_order_cone": {"pr": SOURCE_CONE_PR, "head": SOURCE_CONE_HEAD},
            "compact_swirl": {"pr": SOURCE_SWIRL_PR, "head": SOURCE_SWIRL_HEAD},
            "shoulder_temporal": {"pr": SOURCE_SHOULDER_PR, "head": SOURCE_SHOULDER_HEAD},
            "linear_temporal": {"pr": SOURCE_TEMPORAL_PR, "head": SOURCE_TEMPORAL_HEAD},
            "fixed_render": {"pr": SOURCE_RENDER_PR, "head": SOURCE_RENDER_HEAD},
        },
        "frozen_witness": {
            "swirl_a": SWIRL_A,
            "shoulder_b": SHOULDER_B,
            "shoulder_lambda": SHOULDER_LAMBDA,
            "compact_swirl_epsilon": epsilon,
            "effective_compact_swirl_amplitude": float(SWIRL_A * epsilon),
            "time_envelope": "g(t)=2*(t-.25)",
            "parameter_grid_scan_performed": False,
            "optimization_performed": False,
        },
        "normalization": normalization,
        "midtime_morphology_relative_to_control": {
            "linear": linear_mid,
            "child": child_mid,
        },
        "late_morphology_relative_to_control_descriptive_only": {
            "linear": linear_late,
            "child": child_late,
        },
        "path_relative_to_control": path_relative,
        "nonlinear_oriented_desirability_increments_vs_linear": increments,
        "early_identity_max_abs": identity,
        "structure_preflight": structure,
        "midpoint_energy": {
            "linear": e_linear,
            "child": e_child,
            "child_relative_to_linear": energy_drift,
        },
        "structural_gates_pass": structural_pass,
        "nonlinear_pareto_rule": pareto,
        "combined_child_target_free_pareto_verified": verified,
        "routing": (
            "route this exact one-point child to independent fixed-render/delivery comparison; do not add a third basis"
            if verified
            else "keep #587 as the minimal family; do not retune the failed one-point witness"
        ),
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
    print("combined_child_target_free_pareto_verified=", r["combined_child_target_free_pareto_verified"])
    print("witness=", r["frozen_witness"])
    print("increments=", r["nonlinear_oriented_desirability_increments_vs_linear"])
    print("midtime_morphology=", r["midtime_morphology_relative_to_control"])
    print("paths=", r["path_relative_to_control"])
    print("late_morphology=", r["late_morphology_relative_to_control_descriptive_only"])
    print("structure=", r["structure_preflight"])
    print("energy=", r["midpoint_energy"])
    print("routing=", r["routing"])


if __name__ == "__main__":
    main()
