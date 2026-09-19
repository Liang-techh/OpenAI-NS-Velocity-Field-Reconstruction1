"""Screen one compact vector-potential poloidal basis around the ST052-M linear ramp.

Preregistered in issue #612 before evaluation. The frozen parent is Agent-7 PR
#587's single linear time law. This module adds exactly one divergence-free,
off-midplane poloidal correction generated as curl(A), with

    A=(-y*f, x*f, 0),
    f=R(r^2) Z(z^2),

and an endpoint-zero temporal envelope h(t)=16(t-.25)(.75-t). Its amplitude is
not fitted: at t=.50 on the frozen 41^3 grid it is derived once so the RMS
perturbation equals 0.5% of the baseline linear-ramp velocity RMS. Exactly the
two signs are evaluated. No spatial window, amplitude, time profile, image
numeric target, pressure, or forcing is optimized.

This is expression-capacity evidence only, not PDE validation, visual
correspondence verification, or identification of an OpenAI field.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

import agent7_st052m_linear_temporal_activation as base

TASK_ID = "CR003-ST052M-VECTOR-POTENTIAL-POLOIDAL-088"
PREREG_ISSUE = 612
SOURCE_TEMPORAL_PR = 587
SOURCE_TEMPORAL_HEAD = "0b93819095f6c8576a7571bdc2d2fbef4154944d"
SOURCE_STATIC_PR = 559
SOURCE_STATIC_HEAD = "39b106ad8cb8df2064cabead3a12682089575e74"
SOURCE_CURVATURE_PR = 604

GRID_RESOLUTION = 41
MID_TIME = 0.50
R_SUPPORT = 1.55
Z_SUPPORT_INNER = 0.70
Z_SUPPORT_OUTER = 1.45
PERTURBATION_RMS_FRACTION = 0.005
PARETO_TOL = 5.0e-5
DIVERGENCE_STEP = 1.0e-5
ENERGY_DRIFT_MAX = 0.01
CONDITION_MAX = 12.0
PROBE_SEED = 9176121
PROBES_PER_TIME = 16
ACTIVE_TIMES = (0.375, 0.50, 0.625)

TRUTH = {
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "production_candidate_selected": False,
    "candidate_basis_changed": False,
    "screened_extra_vector_potential_basis": True,
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


def endpoint_zero_envelope(time: float) -> float:
    time = float(time)
    if not (base.TIME_INTERVAL[0] - 1.0e-12 <= time <= base.TIME_INTERVAL[1] + 1.0e-12):
        raise ValueError("time outside frozen interval")
    return float(16.0 * (time - 0.25) * (0.75 - time))


def _poly_bump(q: np.ndarray, a: float, b: float) -> tuple[np.ndarray, np.ndarray]:
    """Normalized ((q-a)(b-q))^4 bump and analytic d/dq, zero outside."""
    q = np.asarray(q, dtype=float)
    value = np.zeros_like(q)
    deriv = np.zeros_like(q)
    mask = (q > a) & (q < b)
    if not np.any(mask):
        return value, deriv
    p = (q[mask] - a) * (b - q[mask])
    pmax = ((b - a) * 0.5) ** 2
    value[mask] = (p / pmax) ** 4
    deriv[mask] = 4.0 * p**3 * (a + b - 2.0 * q[mask]) / pmax**4
    return value, deriv


def vector_potential_correction(points: np.ndarray) -> np.ndarray:
    """Analytic curl of A=(-y*f,x*f,0); axisymmetric and divergence-free."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    r2 = x * x + y * y
    z2 = z * z
    radial, d_radial = _poly_bump(r2, 0.0, R_SUPPORT**2)
    axial, d_axial = _poly_bump(z2, Z_SUPPORT_INNER**2, Z_SUPPORT_OUTER**2)
    f = radial * axial
    f_r2 = d_radial * axial
    f_z = radial * d_axial * (2.0 * z)
    return np.column_stack(
        (-x * f_z, -y * f_z, 2.0 * f + 2.0 * r2 * f_r2)
    )


def vector_velocity(ramp_fn, points, time: float, epsilon: float, sign: int) -> np.ndarray:
    if int(sign) not in (-1, 1):
        raise ValueError("only the two preregistered signs are allowed")
    pts = np.asarray(points, dtype=float)
    return np.asarray(ramp_fn(pts, float(time)), dtype=float) + (
        float(sign) * float(epsilon) * endpoint_zero_envelope(time)
        * vector_potential_correction(pts)
    )


def _grid_points(axis: np.ndarray) -> np.ndarray:
    xx, yy, zz = np.meshgrid(axis, axis, axis, indexing="ij")
    return np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))


def derive_epsilon(ramp_fn) -> tuple[float, dict]:
    axis, velocity = base.morph.prior.sample_velocity_grid(ramp_fn, MID_TIME, GRID_RESOLUTION)
    pts = _grid_points(axis)
    correction = vector_potential_correction(pts).reshape(velocity.shape)
    baseline_rms = float(np.sqrt(np.mean(np.sum(np.square(velocity), axis=-1))))
    correction_rms = float(np.sqrt(np.mean(np.sum(np.square(correction), axis=-1))))
    if correction_rms <= np.finfo(float).tiny or baseline_rms <= np.finfo(float).tiny:
        raise RuntimeError("degenerate RMS normalization")
    epsilon = float(PERTURBATION_RMS_FRACTION * baseline_rms / correction_rms)
    return epsilon, {
        "baseline_velocity_rms": baseline_rms,
        "unscaled_correction_rms": correction_rms,
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
    vector_blocks = []
    for time in ACTIVE_TIMES:
        radius = rng.uniform(0.2, 1.45, size=PROBES_PER_TIME)
        angle = rng.uniform(0.0, 2.0 * np.pi, size=PROBES_PER_TIME)
        z = rng.uniform(-1.40, 1.40, size=PROBES_PER_TIME)
        pts = np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))
        u0 = np.asarray(control_fn(pts, time), float)
        ur = np.asarray(ramp_fn(pts, time), float)
        up = np.asarray(plus_fn(pts, time), float)
        um = np.asarray(minus_fn(pts, time), float)
        baseline_blocks.append((ur - u0).ravel())
        vector_blocks.append(((up - um) / (2.0 * epsilon)).ravel())
    a = np.concatenate(baseline_blocks)
    b = np.concatenate(vector_blocks)
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
        "vector_response_norm_per_epsilon": nb,
    }


def structure_preflight(ramp_fn, plus_fn, minus_fn) -> dict:
    rng = np.random.default_rng(9176122)
    central = rng.uniform(-1.25, 1.25, size=(64, 3))
    central[:, 2] = rng.uniform(-0.65, 0.65, size=64)
    central_errors = {}
    for label, fn in (("plus", plus_fn), ("minus", minus_fn)):
        central_errors[label] = float(
            np.max(np.abs(np.asarray(fn(central, MID_TIME)) - np.asarray(ramp_fn(central, MID_TIME))))
        )

    outside = np.array(
        [
            [1.56, 0.0, 1.0], [-1.56, 0.0, -1.0], [0.0, 1.56, 1.0],
            [0.0, -1.56, -1.0], [0.8, 0.0, 0.69], [0.8, 0.0, 1.46],
            [0.8, 0.0, -1.46],
        ], dtype=float,
    )
    correction_outside = float(np.max(np.abs(vector_potential_correction(outside))))

    pts = rng.uniform(-1.35, 1.35, size=(36, 3))
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
        "outside_vector_correction_max_abs": correction_outside,
        "divergence_fd_max": div,
    }


def endpoint_identity(ramp_fn, plus_fn, minus_fn) -> dict:
    rng = np.random.default_rng(9176123)
    pts = rng.uniform(-1.8, 1.8, size=(96, 3))
    out = {}
    for label, fn in (("plus", plus_fn), ("minus", minus_fn)):
        out[label] = {
            "early_ramp_max_abs": float(np.max(np.abs(fn(pts, 0.25) - ramp_fn(pts, 0.25)))),
            "late_ramp_max_abs": float(np.max(np.abs(fn(pts, 0.75) - ramp_fn(pts, 0.75)))),
        }
    return out


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
    plus_fn = lambda p, t: vector_velocity(ramp_fn, p, t, epsilon, +1)
    minus_fn = lambda p, t: vector_velocity(ramp_fn, p, t, epsilon, -1)

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

    endpoint = endpoint_identity(ramp_fn, plus_fn, minus_fn)
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
        all(v["early_ramp_max_abs"] <= 1.0e-11 and v["late_ramp_max_abs"] <= 1.0e-11 for v in endpoint.values())
        and all(v <= 1.0e-11 for v in structure["central_velocity_perturbation_max_abs"].values())
        and structure["outside_vector_correction_max_abs"] <= 1.0e-12
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
        "source_static_pr": SOURCE_STATIC_PR,
        "source_static_head": SOURCE_STATIC_HEAD,
        "source_curvature_pr": SOURCE_CURVATURE_PR,
        "frozen_screen": {
            "grid_resolution": GRID_RESOLUTION,
            "mid_time": MID_TIME,
            "radial_support_max": R_SUPPORT,
            "axial_support_abs": [Z_SUPPORT_INNER, Z_SUPPORT_OUTER],
            "perturbation_rms_fraction": PERTURBATION_RMS_FRACTION,
            "time_envelope": "16*(t-.25)*(.75-t)",
            "signs": [-1, 1],
            "parameter_scan_performed": False,
            "held_out_pde_residual_evaluated": False,
        },
        "normalization": normalization,
        "morphology_relative_to_control": relative,
        "path_relative_to_control": path_relative,
        "oriented_incremental_desirabilities": oriented,
        "endpoint_identity": endpoint,
        "structure_preflight": structure,
        "midtime_energy": {"linear": energy_linear, **energy},
        "velocity_response_diagnostic": response,
        "structural_gates_pass": structural_pass,
        "local_pareto_rule": pareto,
        "vector_potential_basis_growth_justified": justified,
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
    print("vector_potential_basis_growth_justified=", report["vector_potential_basis_growth_justified"])
    print("pareto=", report["local_pareto_rule"])
    print("oriented=", report["oriented_incremental_desirabilities"])
    print("morphology=", report["morphology_relative_to_control"])
    print("paths=", report["path_relative_to_control"])
    print("response=", report["velocity_response_diagnostic"])
    print("structure=", report["structure_preflight"])


if __name__ == "__main__":
    main()
