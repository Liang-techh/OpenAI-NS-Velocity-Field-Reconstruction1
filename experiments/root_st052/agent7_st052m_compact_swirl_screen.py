"""Screen one axis-regular compact pure-swirl basis around the ST052-M linear ramp.

Preregistered in issue #620 before evaluation. The frozen parent is Agent-7 PR
#587's single linear time law. This module adds exactly one compact axisymmetric
pure-swirl correction

    C_swirl = (-y*S(r^2,z^2), x*S(r^2,z^2), 0),

with normalized fourth-power polynomial bumps, radial support r<1.35 and axial
support 0.70<|z|<1.45. The already-accepted linear activation g(t)=2(t-.25)
is reused; no second time degree is introduced. Amplitude is not fitted: at
t=.50 on the frozen 41^3 grid it is derived once so the actual perturbation RMS
is exactly 0.5% of the baseline linear-ramp velocity RMS. Exactly the two signs
are evaluated; there is no amplitude/window/time-law scan.

This is expression-capacity evidence only, not PDE validation, visual
correspondence verification, or identification of an OpenAI field.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_vector_potential_poloidal_screen as prior

base = prior.base

TASK_ID = "CR003-ST052M-COMPACT-SWIRL-089"
PREREG_ISSUE = 620
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
SOURCE_VECTOR_POTENTIAL_PR = 613
SOURCE_VECTOR_POTENTIAL_HEAD = "f211e706caafc58687345ffa922a723380d4c472"
SOURCE_RENDER_PR = 601

GRID_RESOLUTION = 41
MID_TIME = 0.50
R_SUPPORT = 1.35
Z_SUPPORT_INNER = 0.70
Z_SUPPORT_OUTER = 1.45
PERTURBATION_RMS_FRACTION = 0.005
PARETO_TOL = 5.0e-5
DIVERGENCE_STEP = 1.0e-5
ENERGY_DRIFT_MAX = 0.01
CONDITION_MAX = 12.0
PROBE_SEED = 9176201
PROBES_PER_TIME = 16
ACTIVE_TIMES = (0.375, 0.50, 0.625)

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "candidate_basis_changed": False,
    "screened_extra_compact_swirl_basis": True,
    "parameter_scan_performed": False,
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


def swirl_correction(points: np.ndarray) -> np.ndarray:
    """Axis-regular compact pure swirl (-y*S,x*S,0), analytically divergence-free."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    r2 = x * x + y * y
    z2 = z * z
    radial, _ = prior._poly_bump(r2, 0.0, R_SUPPORT**2)
    axial, _ = prior._poly_bump(z2, Z_SUPPORT_INNER**2, Z_SUPPORT_OUTER**2)
    scalar = radial * axial
    return np.column_stack((-y * scalar, x * scalar, np.zeros_like(scalar)))


def swirl_velocity(ramp_fn, points, time: float, epsilon: float, sign: int) -> np.ndarray:
    if int(sign) not in (-1, 1):
        raise ValueError("only the two preregistered signs are allowed")
    pts = np.asarray(points, dtype=float)
    return np.asarray(ramp_fn(pts, float(time)), dtype=float) + (
        float(sign)
        * float(epsilon)
        * base.activation(time)
        * swirl_correction(pts)
    )


def derive_epsilon(ramp_fn) -> tuple[float, dict]:
    axis, velocity = base.morph.prior.sample_velocity_grid(ramp_fn, MID_TIME, GRID_RESOLUTION)
    pts = prior._grid_points(axis)
    envelope = base.activation(MID_TIME)
    correction = envelope * swirl_correction(pts).reshape(velocity.shape)
    baseline_rms = float(np.sqrt(np.mean(np.sum(np.square(velocity), axis=-1))))
    correction_rms = float(np.sqrt(np.mean(np.sum(np.square(correction), axis=-1))))
    if correction_rms <= np.finfo(float).tiny or baseline_rms <= np.finfo(float).tiny:
        raise RuntimeError("degenerate RMS normalization")
    epsilon = float(PERTURBATION_RMS_FRACTION * baseline_rms / correction_rms)
    return epsilon, {
        "baseline_velocity_rms": baseline_rms,
        "unscaled_midtime_enveloped_correction_rms": correction_rms,
        "midtime_activation": envelope,
        "derived_epsilon": epsilon,
        "scaled_correction_to_velocity_rms_ratio": float(epsilon * correction_rms / baseline_rms),
    }


def _morphology_metrics(velocity_fn, time: float = MID_TIME) -> dict:
    axis, velocity = base.morph.prior.sample_velocity_grid(velocity_fn, time, GRID_RESOLUTION)
    spacing = float(axis[1] - axis[0])
    _, omega = base.morph.prior.vorticity(velocity, spacing)
    return base.morph.enstrophy_moment_metrics(omega, axis)


def response_diagnostic(control_fn, ramp_fn, plus_fn, minus_fn, epsilon: float) -> dict:
    rng = np.random.default_rng(PROBE_SEED)
    baseline_blocks = []
    swirl_blocks = []
    for time in ACTIVE_TIMES:
        radius = rng.uniform(0.2, 1.30, size=PROBES_PER_TIME)
        angle = rng.uniform(0.0, 2.0 * np.pi, size=PROBES_PER_TIME)
        z = rng.uniform(-1.40, 1.40, size=PROBES_PER_TIME)
        pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
        u0 = np.asarray(control_fn(pts, time), float)
        ur = np.asarray(ramp_fn(pts, time), float)
        up = np.asarray(plus_fn(pts, time), float)
        um = np.asarray(minus_fn(pts, time), float)
        baseline_blocks.append((ur - u0).ravel())
        swirl_blocks.append(((up - um) / (2.0 * epsilon)).ravel())
    a = np.concatenate(baseline_blocks)
    b = np.concatenate(swirl_blocks)
    na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    if na <= np.finfo(float).tiny or nb <= np.finfo(float).tiny:
        raise RuntimeError("degenerate response column")
    matrix = np.column_stack((a / na, b / nb))
    singular = np.linalg.svd(matrix, compute_uv=False)
    return {
        "rank": int(np.linalg.matrix_rank(matrix, tol=1.0e-10)),
        "singular_values": singular.tolist(),
        "condition_number": float(singular[0] / singular[-1]),
        "column_cosine": float(np.dot(matrix[:, 0], matrix[:, 1])),
        "baseline_response_norm": na,
        "swirl_response_norm_per_epsilon": nb,
    }


def structure_preflight(ramp_fn, plus_fn, minus_fn) -> dict:
    rng = np.random.default_rng(9176202)
    central = rng.uniform(-1.20, 1.20, size=(64, 3))
    central[:, 2] = rng.uniform(-0.65, 0.65, size=64)
    central_errors = {}
    for label, fn in (("plus", plus_fn), ("minus", minus_fn)):
        central_errors[label] = float(
            np.max(np.abs(np.asarray(fn(central, MID_TIME)) - np.asarray(ramp_fn(central, MID_TIME))))
        )

    outside = np.array(
        [
            [1.36, 0.0, 1.0], [-1.36, 0.0, -1.0], [0.0, 1.36, 1.0],
            [0.0, -1.36, -1.0], [0.8, 0.0, 0.69], [0.8, 0.0, 1.46],
            [0.8, 0.0, -1.46], [0.0, 0.0, 1.0],
        ], dtype=float,
    )
    correction_outside = float(np.max(np.abs(swirl_correction(outside))))

    pts = rng.uniform(-1.25, 1.25, size=(36, 3))
    div = {}
    for label, fn in (("plus", plus_fn), ("minus", minus_fn)):
        div_max = 0.0
        for time in (0.50, 0.625):
            current = np.zeros(len(pts))
            for axis in range(3):
                shift = np.zeros(3)
                shift[axis] = DIVERGENCE_STEP
                up = np.asarray(fn(pts + shift, time), float)
                um = np.asarray(fn(pts - shift, time), float)
                current += (up[:, axis] - um[:, axis]) / (2.0 * DIVERGENCE_STEP)
            div_max = max(div_max, float(np.max(np.abs(current))))
        div[label] = div_max
    return {
        "central_velocity_perturbation_max_abs": central_errors,
        "outside_swirl_correction_max_abs": correction_outside,
        "divergence_fd_max": div,
    }


def early_identity(ramp_fn, plus_fn, minus_fn) -> dict:
    rng = np.random.default_rng(9176203)
    pts = rng.uniform(-1.8, 1.8, size=(96, 3))
    return {
        label: float(np.max(np.abs(fn(pts, 0.25) - ramp_fn(pts, 0.25))))
        for label, fn in (("plus", plus_fn), ("minus", minus_fn))
    }


def _pareto(oriented: dict, structural_pass: bool) -> dict:
    result = {}
    for label in ("plus", "minus"):
        values = np.asarray(list(oriented[label].values()), dtype=float)
        nonworsening = bool(np.all(values >= -PARETO_TOL))
        improvement = bool(np.any(values > PARETO_TOL))
        result[label] = {
            "nonworsening_all": nonworsening,
            "strict_improvement_any": improvement,
            "passes": bool(structural_pass and nonworsening and improvement),
        }
    direction = "+epsilon" if result["plus"]["passes"] else (
        "-epsilon" if result["minus"]["passes"] else None
    )
    return {"tolerance": PARETO_TOL, "by_sign": result, "target_free_pareto_direction": direction}


def run(out: Path) -> dict:
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
    epsilon, normalization = derive_epsilon(ramp_fn)
    plus_fn = lambda p, t: swirl_velocity(ramp_fn, p, t, epsilon, +1)
    minus_fn = lambda p, t: swirl_velocity(ramp_fn, p, t, epsilon, -1)

    metrics_control = _morphology_metrics(control_fn)
    metrics_ramp = _morphology_metrics(ramp_fn)
    metrics_plus = _morphology_metrics(plus_fn)
    metrics_minus = _morphology_metrics(minus_fn)
    relative = {
        "linear": base.morph.compare_metrics(metrics_control, metrics_ramp),
        "plus": base.morph.compare_metrics(metrics_control, metrics_plus),
        "minus": base.morph.compare_metrics(metrics_control, metrics_minus),
    }

    path_control = base.measure_paths(control_fn)
    path_ramp = base.measure_paths(ramp_fn)
    path_plus = base.measure_paths(plus_fn)
    path_minus = base.measure_paths(minus_fn)
    path_relative = {
        "linear": base.compare_paths(path_ramp, path_control),
        "plus": base.compare_paths(path_plus, path_control),
        "minus": base.compare_paths(path_minus, path_control),
    }

    oriented = {}
    for label in ("plus", "minus"):
        oriented[label] = {
            "aspect_gain_increment": float(relative[label]["full_aspect_ratio"] - relative["linear"]["full_aspect_ratio"]),
            "tip_thinning_increment": float(relative["linear"]["smooth_tip_radial_rms"] - relative[label]["smooth_tip_radial_rms"]),
            "turns_fidelity_increment": float(path_relative[label]["mean_absolute_turns_relative"] - path_relative["linear"]["mean_absolute_turns_relative"]),
            "axial_pair_fidelity_increment": float(path_relative[label]["mean_pair_axial_separation_change_relative"] - path_relative["linear"]["mean_pair_axial_separation_change_relative"]),
        }

    identity = early_identity(ramp_fn, plus_fn, minus_fn)
    structure = structure_preflight(ramp_fn, plus_fn, minus_fn)
    response = response_diagnostic(control_fn, ramp_fn, plus_fn, minus_fn, epsilon)
    energy_linear = base.axisymmetric_energy_at_time(ramp_fn, MID_TIME)
    energy = {}
    for label, fn in (("plus", plus_fn), ("minus", minus_fn)):
        value = base.axisymmetric_energy_at_time(fn, MID_TIME)
        energy[label] = {
            "value": value,
            "relative_to_linear": float(value / energy_linear - 1.0),
        }

    structural_pass = bool(
        all(v <= 1.0e-11 for v in identity.values())
        and all(v <= 1.0e-11 for v in structure["central_velocity_perturbation_max_abs"].values())
        and structure["outside_swirl_correction_max_abs"] <= 1.0e-12
        and all(v <= 1.0e-5 for v in structure["divergence_fd_max"].values())
        and all(abs(v["relative_to_linear"]) <= ENERGY_DRIFT_MAX for v in energy.values())
        and response["rank"] >= 2
        and response["condition_number"] <= CONDITION_MAX
    )
    pareto = _pareto(oriented, structural_pass)
    justified = bool(pareto["target_free_pareto_direction"] is not None)

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_temporal_pr": SOURCE_TEMPORAL_PR,
        "source_temporal_head": SOURCE_TEMPORAL_HEAD,
        "source_vector_potential_pr": SOURCE_VECTOR_POTENTIAL_PR,
        "source_vector_potential_head": SOURCE_VECTOR_POTENTIAL_HEAD,
        "source_render_pr": SOURCE_RENDER_PR,
        "frozen_screen": {
            "grid_resolution": GRID_RESOLUTION,
            "mid_time": MID_TIME,
            "radial_support_max": R_SUPPORT,
            "axial_support_abs": [Z_SUPPORT_INNER, Z_SUPPORT_OUTER],
            "perturbation_rms_fraction": PERTURBATION_RMS_FRACTION,
            "time_envelope": "g(t)=2*(t-.25)",
            "signs": [-1, 1],
            "parameter_scan_performed": False,
            "held_out_pde_residual_evaluated": False,
        },
        "normalization": normalization,
        "morphology_relative_to_control": relative,
        "path_relative_to_control": path_relative,
        "oriented_incremental_desirabilities": oriented,
        "early_identity": identity,
        "structure_preflight": structure,
        "midtime_energy": {"linear": energy_linear, **energy},
        "velocity_response_diagnostic": response,
        "structural_gates_pass": structural_pass,
        "local_pareto_rule": pareto,
        "compact_swirl_basis_growth_justified": justified,
        **TRUTH,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    print("epsilon=", report["normalization"]["derived_epsilon"])
    print("compact_swirl_basis_growth_justified=", report["compact_swirl_basis_growth_justified"])
    print("pareto=", report["local_pareto_rule"])
    print("oriented=", report["oriented_incremental_desirabilities"])
    print("morphology=", report["morphology_relative_to_control"])
    print("paths=", report["path_relative_to_control"])
    print("response=", report["velocity_response_diagnostic"])
    print("structure=", report["structure_preflight"])


if __name__ == "__main__":
    main()
