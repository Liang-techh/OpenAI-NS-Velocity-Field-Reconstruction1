"""Screen one endpoint-zero temporal curvature direction around the ST052-M ramp.

Preregistered in issue #603 before evaluation.  The frozen base field is the
Agent-7 PR #587 linear-ramp child.  This module does not optimize or select a
new coefficient.  It evaluates exactly lambda=+/-0.02 for

    g_lambda(t) = 2*(t-.25) + lambda*16*(t-.25)*(.75-t),

so the t=.25 control endpoint and t=.75 static-#559 endpoint are preserved
exactly.  The scientific question is whether this second temporal direction
has a target-free first-order Pareto sign across morphology and material-path
observables.  Independent velocity-response rank is recorded separately: an
independent direction is not, by itself, a reason to grow the basis.

This is expression-capacity evidence only.  It is not PDE validation, visual
correspondence verification, or identification of the OpenAI field.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_linear_temporal_activation as base

TASK_ID = "CR003-ST052M-TEMPORAL-CURVATURE-SENSITIVITY-087"
PREREG_ISSUE = 603
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
SOURCE_STATIC_PR = 559
SOURCE_STATIC_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
SOURCE_CONVERGENCE_PR = 594
SOURCE_RENDER_PR = 601

EPSILON = 0.02
ZERO_TOL = 1.0e-6
MID_TIME = 0.50
ACTIVE_TIMES = (0.375, 0.50, 0.625)
GRID_RESOLUTION = 41
PROBE_SEED = 9176031
PROBES_PER_TIME = 16
ENDPOINT_SEED = 9176032
ENDPOINT_PROBES = 64

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "candidate_basis_changed": False,
    "screened_extra_temporal_basis": True,
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


def curvature_basis(time: float) -> float:
    time = float(time)
    if not (base.TIME_INTERVAL[0] - 1.0e-12 <= time <= base.TIME_INTERVAL[1] + 1.0e-12):
        raise ValueError("time outside frozen interval")
    return float(16.0 * (time - 0.25) * (0.75 - time))


def curved_activation(time: float, lam: float) -> float:
    if abs(float(lam)) > EPSILON + 1.0e-15:
        raise ValueError("only the preregistered +/-epsilon screen is allowed")
    return float(base.activation(time) + float(lam) * curvature_basis(time))


def curved_velocity(field, raw, points, time, redistribution_scale, beta, lam):
    g = curved_activation(time, lam)
    u = base.morph.prior.taper.taper_velocity(
        field,
        raw,
        np.asarray(points, dtype=float),
        float(time),
        redistribution_scale=redistribution_scale,
        tau=base.TAPER_TAU * g,
        taper_scale=1.0,
    )
    return base.morph.prior.comp.apply_swirl_compensation(points, u, beta * g)


def morphology_relative(velocity_fn, control_fn, time: float = MID_TIME) -> dict:
    axis, uc = base.morph.prior.sample_velocity_grid(control_fn, time, GRID_RESOLUTION)
    _, uu = base.morph.prior.sample_velocity_grid(velocity_fn, time, GRID_RESOLUTION)
    spacing = float(axis[1] - axis[0])
    _, oc = base.morph.prior.vorticity(uc, spacing)
    _, ou = base.morph.prior.vorticity(uu, spacing)
    mc = base.morph.enstrophy_moment_metrics(oc, axis)
    mu = base.morph.enstrophy_moment_metrics(ou, axis)
    return base.morph.compare_metrics(mc, mu)


def response_diagnostic(control_fn, ramp_fn, plus_fn, minus_fn) -> dict:
    rng = np.random.default_rng(PROBE_SEED)
    linear_blocks = []
    curvature_blocks = []
    for time in ACTIVE_TIMES:
        radius = rng.uniform(0.2, 1.45, size=PROBES_PER_TIME)
        angle = rng.uniform(0.0, 2.0 * np.pi, size=PROBES_PER_TIME)
        z = rng.uniform(-1.45, 1.45, size=PROBES_PER_TIME)
        pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
        u0 = np.asarray(control_fn(pts, time), float)
        linear_blocks.append((np.asarray(ramp_fn(pts, time), float) - u0).ravel())
        curvature_blocks.append(
            ((np.asarray(plus_fn(pts, time), float) - np.asarray(minus_fn(pts, time), float))
             / (2.0 * EPSILON)).ravel()
        )
    linear = np.concatenate(linear_blocks)
    curvature = np.concatenate(curvature_blocks)
    n_linear = float(np.linalg.norm(linear))
    n_curvature = float(np.linalg.norm(curvature))
    if n_linear <= np.finfo(float).tiny or n_curvature <= np.finfo(float).tiny:
        raise RuntimeError("degenerate response column")
    matrix = np.column_stack((linear / n_linear, curvature / n_curvature))
    singular = np.linalg.svd(matrix, compute_uv=False)
    return {
        "rank": int(np.linalg.matrix_rank(matrix, tol=1.0e-10)),
        "singular_values": singular.tolist(),
        "condition_number": float(singular[0] / singular[-1]),
        "column_cosine": float(np.dot(matrix[:, 0], matrix[:, 1])),
        "linear_response_norm": n_linear,
        "curvature_response_norm_per_lambda": n_curvature,
    }


def endpoint_identity(control_fn, static_fn, plus_fn, minus_fn) -> dict:
    rng = np.random.default_rng(ENDPOINT_SEED)
    radius = rng.uniform(0.0, 1.9, size=ENDPOINT_PROBES)
    angle = rng.uniform(0.0, 2.0 * np.pi, size=ENDPOINT_PROBES)
    z = rng.uniform(-1.9, 1.9, size=ENDPOINT_PROBES)
    pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
    early = np.asarray(control_fn(pts, 0.25), float)
    late = np.asarray(static_fn(pts, 0.75), float)
    out = {}
    for label, fn in (("plus", plus_fn), ("minus", minus_fn)):
        out[label] = {
            "early_control_max_abs": float(np.max(np.abs(np.asarray(fn(pts, 0.25), float) - early))),
            "late_static_max_abs": float(np.max(np.abs(np.asarray(fn(pts, 0.75), float) - late))),
        }
    return out


def activation_range() -> dict:
    times = np.linspace(0.25, 0.75, 401)
    out = {}
    for label, lam in (("plus", EPSILON), ("minus", -EPSILON)):
        values = np.asarray([curved_activation(t, lam) for t in times], float)
        out[label] = {"min": float(values.min()), "max": float(values.max())}
    return out


def _central_derivative(plus: float, minus: float) -> float:
    return float((float(plus) - float(minus)) / (2.0 * EPSILON))


def local_pareto(oriented: dict) -> dict:
    values = np.asarray(list(oriented.values()), dtype=float)
    plus_ok = bool(np.all(values >= -ZERO_TOL) and np.any(values > ZERO_TOL))
    minus_values = -values
    minus_ok = bool(np.all(minus_values >= -ZERO_TOL) and np.any(minus_values > ZERO_TOL))
    direction = "+lambda" if plus_ok else ("-lambda" if minus_ok else None)
    return {
        "zero_tolerance": ZERO_TOL,
        "plus_lambda_nonworsening_all_and_improving_one": plus_ok,
        "minus_lambda_nonworsening_all_and_improving_one": minus_ok,
        "target_free_local_pareto_direction": direction,
    }


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
    static_fn = lambda p, t: base.static_velocity(field, raw, p, t, redistribution_scale, beta)
    ramp_fn = lambda p, t: base.ramp_velocity(field, raw, p, t, redistribution_scale, beta)
    plus_fn = lambda p, t: curved_velocity(field, raw, p, t, redistribution_scale, beta, EPSILON)
    minus_fn = lambda p, t: curved_velocity(field, raw, p, t, redistribution_scale, beta, -EPSILON)

    morphology = {
        "linear": morphology_relative(ramp_fn, control_fn),
        "plus": morphology_relative(plus_fn, control_fn),
        "minus": morphology_relative(minus_fn, control_fn),
    }

    paths_control = base.measure_paths(control_fn)
    paths_linear = base.measure_paths(ramp_fn)
    paths_plus = base.measure_paths(plus_fn)
    paths_minus = base.measure_paths(minus_fn)
    path_relative = {
        "linear": base.compare_paths(paths_linear, paths_control),
        "plus": base.compare_paths(paths_plus, paths_control),
        "minus": base.compare_paths(paths_minus, paths_control),
    }

    raw_sensitivities = {
        "aspect_change_per_lambda": _central_derivative(
            morphology["plus"]["full_aspect_ratio"], morphology["minus"]["full_aspect_ratio"]
        ),
        "tip_radial_change_per_lambda": _central_derivative(
            morphology["plus"]["smooth_tip_radial_rms"], morphology["minus"]["smooth_tip_radial_rms"]
        ),
        "path_turns_relative_per_lambda": _central_derivative(
            path_relative["plus"]["mean_absolute_turns_relative"],
            path_relative["minus"]["mean_absolute_turns_relative"],
        ),
        "path_pair_axial_relative_per_lambda": _central_derivative(
            path_relative["plus"]["mean_pair_axial_separation_change_relative"],
            path_relative["minus"]["mean_pair_axial_separation_change_relative"],
        ),
    }
    oriented = {
        "aspect_desirability": raw_sensitivities["aspect_change_per_lambda"],
        "tip_thinning_desirability": -raw_sensitivities["tip_radial_change_per_lambda"],
        "path_turns_desirability": raw_sensitivities["path_turns_relative_per_lambda"],
        "path_pair_axial_desirability": raw_sensitivities["path_pair_axial_relative_per_lambda"],
    }
    pareto = local_pareto(oriented)

    response = response_diagnostic(control_fn, ramp_fn, plus_fn, minus_fn)
    endpoints = endpoint_identity(control_fn, static_fn, plus_fn, minus_fn)
    activation_bounds = activation_range()
    structure = {
        "plus": base.structure_preflight(plus_fn),
        "minus": base.structure_preflight(minus_fn),
    }

    endpoint_max = max(
        endpoints["plus"]["early_control_max_abs"],
        endpoints["plus"]["late_static_max_abs"],
        endpoints["minus"]["early_control_max_abs"],
        endpoints["minus"]["late_static_max_abs"],
    )
    support_max = max(structure["plus"]["support_max_abs"], structure["minus"]["support_max_abs"])
    divergence_max = max(
        structure["plus"]["divergence_fd_max"], structure["minus"]["divergence_fd_max"]
    )
    structural_preflight_pass = bool(
        endpoint_max <= 1.0e-11
        and support_max <= base.CRITERIA["support_max_abs"]
        and divergence_max <= base.CRITERIA["divergence_fd_max"]
        and activation_bounds["plus"]["min"] >= -1.0e-12
        and activation_bounds["plus"]["max"] <= 1.0 + 1.0e-12
        and activation_bounds["minus"]["min"] >= -1.0e-12
        and activation_bounds["minus"]["max"] <= 1.0 + 1.0e-12
    )
    response_independent = bool(response["rank"] >= 2 and response["condition_number"] <= 8.0)
    growth_justified = bool(
        structural_preflight_pass
        and response_independent
        and pareto["target_free_local_pareto_direction"] is not None
    )

    report = {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_temporal_pr": SOURCE_TEMPORAL_PR,
        "source_temporal_head": SOURCE_TEMPORAL_HEAD,
        "source_static_pr": SOURCE_STATIC_PR,
        "source_static_head": SOURCE_STATIC_HEAD,
        "source_convergence_pr": SOURCE_CONVERGENCE_PR,
        "source_render_pr": SOURCE_RENDER_PR,
        "frozen_screen": {
            "epsilon": EPSILON,
            "lambda_values": [-EPSILON, EPSILON],
            "base_activation": "2*(t-0.25)",
            "curvature_basis": "16*(t-0.25)*(0.75-t)",
            "redistribution_gain": base.morph.prior.base.GAIN,
            "taper_tau_static_endpoint": base.TAPER_TAU,
            "shoulder_beta_static_endpoint": beta,
            "post_transform_common_scale": 1.0,
            "grid_resolution": GRID_RESOLUTION,
            "morphology_time": MID_TIME,
            "parameter_scan_performed": False,
        },
        "activation_bounds": activation_bounds,
        "endpoint_identity": endpoints,
        "morphology_relative_to_control": morphology,
        "path_relative_to_control": path_relative,
        "raw_central_sensitivities": raw_sensitivities,
        "oriented_desirability_sensitivities": oriented,
        "local_pareto_rule": pareto,
        "velocity_response_diagnostic": response,
        "structure_preflight": structure,
        "structural_preflight_pass": structural_preflight_pass,
        "response_independent": response_independent,
        "temporal_basis_growth_justified": growth_justified,
        **TRUTH,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    r = run(args.out)
    print("temporal_basis_growth_justified=", r["temporal_basis_growth_justified"])
    print("pareto=", r["local_pareto_rule"])
    print("oriented_sensitivities=", r["oriented_desirability_sensitivities"])
    print("response=", r["velocity_response_diagnostic"])
    print("morphology=", r["morphology_relative_to_control"])
    print("paths=", r["path_relative_to_control"])
    print("structure=", r["structure_preflight"])


if __name__ == "__main__":
    main()
