"""Observable-space Jacobian for the five existing ST052-M morphology controls.

Preregistered in issue #1067 and stacked exactly on Agent-7 PR #1060 head.
This increment changes no selected candidate and adds no basis.  It maps the
already-existing five local velocity tangents into deterministic,
renderer-independent morphology observables so future discrepancies can be
routed to the smallest existing control before basis growth.

The observable protocol is autonomous repository diagnostics only.  No public
OpenAI numerical target, camera registration, pixel loss, pressure/forcing
refit, held-out PDE residual, or candidate promotion is introduced.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_space_time_local_five_channel_capacity as parent
import agent7_st052m_threshold_free_morphology as morph

TASK_ID = "CR003-ST052M-OBSERVABLE-SPACE-FIVE-CHANNEL-JACOBIAN-127"
PREREG_ISSUE = 1067
SOURCE_PARENT_PR = 1060
SOURCE_PARENT_HEAD = "637a86d84928c543526d0a7068b699ca90fe38a8"

TIMES = (0.25, 0.50, 0.75)
RING_RADII = (0.40, 0.70, 1.00, 1.30)
RING_ABS_Z = (0.30, 0.80)
AZIMUTH_COUNT = 16
CORE_RADII = tuple(float(v) for v in np.linspace(0.20, 1.60, 15))
CORE_SPEED_RADII = (0.40, 0.80, 1.20)
VORTICITY_TIME = 0.50
VORTICITY_GRID = 25
EPSILONS = (1.0e-3, 5.0e-4)
PRIMARY_EPSILON = 5.0e-4

RANK_TARGET = 5
CONDITION_MAX = 25.0
COSINE_MAX_ABS = 0.995
DERIVATIVE_DRIFT_MAX = 0.05
DERIVATIVE_COSINE_MIN = 0.999

OBSERVABLE_GROUPS = (
    "radial_ring",
    "swirl_ring",
    "axial_ring",
    "core_speed_width",
    "vorticity_shape",
)

VORTICITY_KEYS = (
    "full_axial_rms",
    "full_radial_rms",
    "full_aspect_ratio",
    "smooth_tip_radial_rms",
    "smooth_tip_axial_rms",
    "smooth_tip_enstrophy_fraction",
    "central_radial_rms",
)

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "coefficient_selected": False,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "st006_comparison_performed": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "renderer_or_camera_fit_used": False,
    "observable_space_sensitivities_recorded": True,
    "direct_visualization_fingerprint_improvement": 0.0,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _joint():
    return parent.time_local.regional.joint


def _assert_source_lock() -> None:
    parent._assert_source_lock()
    if parent.TASK_ID != "CR003-ST052M-SPACE-TIME-LOCAL-FIVE-CHANNEL-CAPACITY-126":
        raise RuntimeError("#1060 parent task identity drifted")
    if parent.PREREG_ISSUE != 1057:
        raise RuntimeError("#1060 preregistration identity drifted")
    if parent.SOURCE_PARENT_PR != 1050:
        raise RuntimeError("#1060 parent lineage drifted")
    if tuple(_joint().CHANNELS) != (
        "amplitude",
        "radial_shape",
        "axial_turnover",
        "temporal_curvature",
        "toroidal_swirl",
    ):
        raise RuntimeError("five-channel order drifted")
    if tuple(EPSILONS) != (1.0e-3, 5.0e-4) or PRIMARY_EPSILON != 5.0e-4:
        raise RuntimeError("preregistered finite-difference steps drifted")


def _build_field_and_tangents() -> tuple[Callable, dict[str, Callable], dict[str, Any]]:
    """Build the unchanged live field and five existing local tangent callables."""
    _assert_source_lock()
    joint = _joint()
    parent_fn, s_alphas, s_calibrations, s_children, _k_alphas, k_children = (
        joint.axial.build_local_families()
    )
    live_alpha = float(s_alphas[joint.axial.S_LIVE])
    live_fn = s_children[joint.axial.S_LIVE]

    def amplitude(points: np.ndarray, time: float) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        return np.asarray(live_fn(pts, time), dtype=float) - np.asarray(
            parent_fn(pts, time), dtype=float
        )

    def radial_shape(points: np.ndarray, time: float) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        lo = np.asarray(k_children[joint.axial.radial_shape.K_MINUS](pts, time), dtype=float)
        hi = np.asarray(k_children[joint.axial.radial_shape.K_PLUS](pts, time), dtype=float)
        return (hi - lo) / (
            joint.axial.radial_shape.K_PLUS - joint.axial.radial_shape.K_MINUS
        )

    def axial_turnover(points: np.ndarray, time: float) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        lo = np.asarray(s_children[joint.axial.S_MINUS](pts, time), dtype=float)
        hi = np.asarray(s_children[joint.axial.S_PLUS](pts, time), dtype=float)
        return (hi - lo) / (joint.axial.S_PLUS - joint.axial.S_MINUS)

    def temporal_curvature(points: np.ndarray, time: float) -> np.ndarray:
        return np.asarray(
            joint.temporal.temporal_curvature_tangent(
                np.asarray(points, dtype=float), float(time), live_alpha
            ),
            dtype=float,
        )

    def toroidal_swirl(points: np.ndarray, time: float) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        return float(joint.toroidal.activation_g1(float(time))) * np.asarray(
            joint.toroidal.toroidal_unit_velocity(pts), dtype=float
        )

    tangents = {
        "amplitude": amplitude,
        "radial_shape": radial_shape,
        "axial_turnover": axial_turnover,
        "temporal_curvature": temporal_curvature,
        "toroidal_swirl": toroidal_swirl,
    }
    metadata = {
        "live_alpha": live_alpha,
        "live_alpha_relative_mismatch_vs_exact_775": float(
            s_calibrations[joint.axial.S_LIVE]["relative_mismatch_vs_exact_775"]
        ),
        "channel_order": list(_joint().CHANNELS),
    }
    return live_fn, tangents, metadata


def _azimuths() -> np.ndarray:
    return np.arange(AZIMUTH_COUNT, dtype=float) * (2.0 * np.pi / AZIMUTH_COUNT)


def _ring_points(radius: float, z: float) -> np.ndarray:
    angles = _azimuths()
    return np.column_stack(
        (
            float(radius) * np.cos(angles),
            float(radius) * np.sin(angles),
            np.full(AZIMUTH_COUNT, float(z), dtype=float),
        )
    )


def _cylindrical(points: np.ndarray, velocity: np.ndarray) -> np.ndarray:
    return _joint().toroidal.cylindrical_components(
        np.asarray(points, dtype=float), np.asarray(velocity, dtype=float)
    )


def _core_speed_profile(field_fn: Callable, time: float) -> tuple[np.ndarray, float]:
    means: list[float] = []
    for radius in CORE_RADII:
        pts = _ring_points(float(radius), 0.0)
        vel = np.asarray(field_fn(pts, float(time)), dtype=float)
        means.append(float(np.mean(np.linalg.norm(vel, axis=1))))
    profile = np.asarray(means, dtype=float)
    if not np.all(np.isfinite(profile)) or np.any(profile < 0.0):
        raise RuntimeError("nonfinite core-speed profile")
    total = float(np.sum(profile))
    if total <= np.finfo(float).tiny:
        raise RuntimeError("degenerate core-speed profile")
    radii = np.asarray(CORE_RADII, dtype=float)
    weighted_rms = float(np.sqrt(np.sum(radii * radii * profile) / total))
    return profile, weighted_rms


def _vorticity_observables(field_fn: Callable) -> dict[str, float]:
    axis, velocity = morph.prior.sample_velocity_grid(
        field_fn, VORTICITY_TIME, VORTICITY_GRID
    )
    spacing = float(axis[1] - axis[0])
    _omega, magnitude = morph.prior.vorticity(velocity, spacing)
    metrics = morph.enstrophy_moment_metrics(magnitude, axis)
    out = {key: float(metrics[key]) for key in VORTICITY_KEYS}
    if not all(np.isfinite(v) for v in out.values()):
        raise RuntimeError("nonfinite vorticity morphology observable")
    return out


def observable_vector(field_fn: Callable) -> dict[str, float]:
    """Return the frozen renderer-independent morphology observable vector."""
    out: dict[str, float] = {}
    for time in TIMES:
        for abs_z in RING_ABS_Z:
            for sign in (-1.0, 1.0):
                z = float(sign * abs_z)
                for radius in RING_RADII:
                    pts = _ring_points(float(radius), z)
                    vel = np.asarray(field_fn(pts, float(time)), dtype=float)
                    cyl = _cylindrical(pts, vel)
                    suffix = f"t={time:.3f}/z={z:+.2f}/r={radius:.2f}"
                    out[f"radial_ring/{suffix}"] = float(np.mean(cyl[:, 0]))
                    out[f"swirl_ring/{suffix}"] = float(np.mean(cyl[:, 1]))
                    out[f"axial_ring/{suffix}"] = float(sign * np.mean(cyl[:, 2]))

        profile, weighted_rms = _core_speed_profile(field_fn, float(time))
        out[f"core_speed_width/t={time:.3f}/speed_weighted_rms_radius"] = weighted_rms
        for radius in CORE_SPEED_RADII:
            idx = int(np.argmin(np.abs(np.asarray(CORE_RADII) - float(radius))))
            if abs(CORE_RADII[idx] - float(radius)) > 1.0e-12:
                raise RuntimeError("frozen core-speed probe radius not represented")
            out[f"core_speed_width/t={time:.3f}/mean_speed_r={radius:.2f}"] = float(
                profile[idx]
            )

    for key, value in _vorticity_observables(field_fn).items():
        out[f"vorticity_shape/t={VORTICITY_TIME:.3f}/{key}"] = float(value)

    if not out or not all(np.isfinite(v) for v in out.values()):
        raise RuntimeError("invalid morphology observable vector")
    return out


def _baseline_speed_scale(field_fn: Callable) -> float:
    samples: list[np.ndarray] = []
    for time in TIMES:
        for abs_z in RING_ABS_Z:
            for sign in (-1.0, 1.0):
                for radius in RING_RADII:
                    pts = _ring_points(float(radius), float(sign * abs_z))
                    samples.append(np.asarray(field_fn(pts, float(time)), dtype=float))
        for radius in CORE_RADII:
            pts = _ring_points(float(radius), 0.0)
            samples.append(np.asarray(field_fn(pts, float(time)), dtype=float))
    values = np.vstack(samples)
    scale = float(np.sqrt(np.mean(np.sum(values * values, axis=1))))
    if not np.isfinite(scale) or scale <= np.finfo(float).tiny:
        raise RuntimeError("invalid baseline RMS-speed scale")
    return scale


def _observable_scale(name: str, baseline_speed_scale: float) -> float:
    if name.startswith(("radial_ring/", "swirl_ring/", "axial_ring/")):
        return float(baseline_speed_scale)
    if name.startswith("core_speed_width/") and "/mean_speed_" in name:
        return float(baseline_speed_scale)
    return 1.0


def _perturbed_field(
    base_fn: Callable, tangent_fn: Callable, signed_epsilon: float
) -> Callable:
    def field(points: np.ndarray, time: float) -> np.ndarray:
        pts = np.asarray(points, dtype=float)
        return np.asarray(base_fn(pts, float(time)), dtype=float) + float(
            signed_epsilon
        ) * np.asarray(tangent_fn(pts, float(time)), dtype=float)

    return field


def _scaled_derivative(
    base_fn: Callable,
    tangent_fn: Callable,
    epsilon: float,
    names: list[str],
    scales: np.ndarray,
) -> np.ndarray:
    plus = observable_vector(_perturbed_field(base_fn, tangent_fn, float(epsilon)))
    minus = observable_vector(_perturbed_field(base_fn, tangent_fn, -float(epsilon)))
    if list(plus) != names or list(minus) != names:
        raise RuntimeError("observable ordering drifted under perturbation")
    raw = np.asarray(
        [(plus[name] - minus[name]) / (2.0 * float(epsilon)) for name in names],
        dtype=float,
    )
    scaled = raw / scales
    if not np.all(np.isfinite(scaled)):
        raise RuntimeError("nonfinite observable derivative")
    return scaled


def _group_for_name(name: str) -> str:
    prefix = name.split("/", 1)[0]
    if prefix not in OBSERVABLE_GROUPS:
        raise RuntimeError(f"unregistered observable group: {prefix}")
    return prefix


def _column_summary(
    names: list[str], derivative: np.ndarray
) -> dict[str, Any]:
    vector = np.asarray(derivative, dtype=float)
    power = vector * vector
    total = float(np.sum(power))
    if total <= np.finfo(float).tiny:
        fractions = {group: 0.0 for group in OBSERVABLE_GROUPS}
    else:
        fractions = {
            group: float(
                np.sum(
                    [
                        power[i]
                        for i, name in enumerate(names)
                        if _group_for_name(name) == group
                    ]
                )
                / total
            )
            for group in OBSERVABLE_GROUPS
        }
    order = np.argsort(-np.abs(vector))
    dominant = [
        {"observable": names[int(i)], "scaled_derivative": float(vector[int(i)])}
        for i in order[: min(8, len(order))]
    ]
    return {
        "scaled_derivative_l2": float(np.linalg.norm(vector)),
        "observable_group_energy_fractions": fractions,
        "dominant_observables": dominant,
    }


def observable_space_jacobian() -> dict[str, Any]:
    base_fn, tangents, metadata = _build_field_and_tangents()
    baseline = observable_vector(base_fn)
    names = list(baseline)
    baseline_speed = _baseline_speed_scale(base_fn)
    scales = np.asarray(
        [_observable_scale(name, baseline_speed) for name in names], dtype=float
    )
    if np.any(scales <= 0.0) or not np.all(np.isfinite(scales)):
        raise RuntimeError("invalid observable scales")

    derivatives: dict[float, dict[str, np.ndarray]] = {}
    for epsilon in EPSILONS:
        derivatives[float(epsilon)] = {
            channel: _scaled_derivative(
                base_fn, tangents[channel], float(epsilon), names, scales
            )
            for channel in _joint().CHANNELS
        }

    fine = derivatives[PRIMARY_EPSILON]
    normalized_columns: list[np.ndarray] = []
    fine_norms: dict[str, float] = {}
    stability: dict[str, Any] = {}
    summaries: dict[str, Any] = {}
    for channel in _joint().CHANNELS:
        vf = fine[channel]
        vc = derivatives[EPSILONS[0]][channel]
        nf = float(np.linalg.norm(vf))
        nc = float(np.linalg.norm(vc))
        fine_norms[channel] = nf
        if nf <= np.finfo(float).tiny:
            normalized_columns.append(np.zeros_like(vf))
            rel_drift = None
            cosine = None
            stable = False
        else:
            normalized_columns.append(vf / nf)
            rel_drift = float(np.linalg.norm(vf - vc) / nf)
            cosine = (
                None
                if nc <= np.finfo(float).tiny
                else float(np.dot(vf, vc) / (nf * nc))
            )
            stable = bool(
                rel_drift <= DERIVATIVE_DRIFT_MAX
                and cosine is not None
                and cosine >= DERIVATIVE_COSINE_MIN
            )
        stability[channel] = {
            "coarse_epsilon": EPSILONS[0],
            "fine_epsilon": PRIMARY_EPSILON,
            "coarse_scaled_derivative_l2": nc,
            "fine_scaled_derivative_l2": nf,
            "relative_column_drift": rel_drift,
            "coarse_fine_cosine": cosine,
            "drift_gate": DERIVATIVE_DRIFT_MAX,
            "cosine_gate": DERIVATIVE_COSINE_MIN,
            "passes": stable,
        }
        summaries[channel] = _column_summary(names, vf)

    matrix = np.column_stack(normalized_columns)
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    condition = (
        None
        if singular[-1] <= np.finfo(float).tiny
        else float(singular[0] / singular[-1])
    )
    gram = matrix.T @ matrix
    cosines: dict[str, float | None] = {}
    for i, j in itertools.combinations(range(len(_joint().CHANNELS)), 2):
        ni = fine_norms[_joint().CHANNELS[i]]
        nj = fine_norms[_joint().CHANNELS[j]]
        key = f"{_joint().CHANNELS[i]}__{_joint().CHANNELS[j]}"
        cosines[key] = None if min(ni, nj) <= np.finfo(float).tiny else float(gram[i, j])
    finite_cosines = [abs(v) for v in cosines.values() if v is not None]
    max_abs_cosine = None if not finite_cosines else float(max(finite_cosines))

    capacity_pass = bool(
        rank == RANK_TARGET
        and condition is not None
        and condition <= CONDITION_MAX
        and max_abs_cosine is not None
        and max_abs_cosine < COSINE_MAX_ABS
    )
    stability_pass = bool(all(row["passes"] for row in stability.values()))

    baseline_by_group = {
        group: int(sum(_group_for_name(name) == group for name in names))
        for group in OBSERVABLE_GROUPS
    }
    return {
        **metadata,
        "observable_count": int(len(names)),
        "observable_count_by_group": baseline_by_group,
        "observable_names": names,
        "baseline_rms_speed_scale": baseline_speed,
        "baseline_observables": {name: float(baseline[name]) for name in names},
        "primary_epsilon": PRIMARY_EPSILON,
        "column_norms_in_scaled_observable_space": fine_norms,
        "normalized_column_rank": rank,
        "rank_target": RANK_TARGET,
        "normalized_condition": condition,
        "condition_gate": CONDITION_MAX,
        "pairwise_cosines": cosines,
        "max_abs_pairwise_cosine": max_abs_cosine,
        "cosine_abs_gate": COSINE_MAX_ABS,
        "singular_values": [float(v) for v in singular],
        "capacity_gate_passes": capacity_pass,
        "derivative_stability": stability,
        "all_derivative_stability_gates_pass": stability_pass,
        "channels": summaries,
    }


def structural_endpoint_guards() -> dict[str, Any]:
    base_fn, tangents, _ = _build_field_and_tangents()
    del base_fn
    probe = np.vstack(
        [
            _ring_points(0.70, 0.30),
            _ring_points(1.00, 0.80),
            _ring_points(1.30, -0.80),
        ]
    )
    start_max = {
        channel: float(np.max(np.abs(tangent(probe, 0.25))))
        for channel, tangent in tangents.items()
    }
    temporal_end = float(
        np.max(np.abs(tangents["temporal_curvature"](probe, 0.75)))
    )
    return {
        "start_tangent_max_abs_by_channel": start_max,
        "all_current_channels_exact_start_identity": bool(
            max(start_max.values()) <= 1.0e-12
        ),
        "temporal_curvature_end_max_abs": temporal_end,
        "temporal_curvature_exact_end_identity": bool(temporal_end <= 1.0e-12),
        "tolerance": 1.0e-12,
    }


def build_report() -> dict[str, Any]:
    jacobian = observable_space_jacobian()
    guards = structural_endpoint_guards()
    scientific_pass = bool(
        jacobian["capacity_gate_passes"]
        and jacobian["all_derivative_stability_gates_pass"]
        and guards["all_current_channels_exact_start_identity"]
        and guards["temporal_curvature_exact_end_identity"]
    )
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "frozen_protocol": {
            "channels": list(_joint().CHANNELS),
            "times": list(TIMES),
            "ring_radii": list(RING_RADII),
            "ring_abs_z": list(RING_ABS_Z),
            "azimuth_count": AZIMUTH_COUNT,
            "core_radii": list(CORE_RADII),
            "core_speed_radii": list(CORE_SPEED_RADII),
            "vorticity_time": VORTICITY_TIME,
            "vorticity_grid": VORTICITY_GRID,
            "epsilons": list(EPSILONS),
            "primary_epsilon": PRIMARY_EPSILON,
            "rank_target": RANK_TARGET,
            "condition_gate": CONDITION_MAX,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "derivative_drift_gate": DERIVATIVE_DRIFT_MAX,
            "derivative_cosine_gate": DERIVATIVE_COSINE_MIN,
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "observable_space_jacobian": jacobian,
        "structural_endpoint_guards": guards,
        "decision": {
            "observable_space_capacity_and_linearization_pass": scientific_pass,
            "additional_basis_dimension_justified_by_this_audit": False,
            "if_pass_next_action": (
                "route a later candidate-bound morphology discrepancy to the smallest "
                "existing channel with leverage in its observable family"
            ),
            "if_fail_next_action": (
                "inspect/reparameterize the weak observable-space direction or one "
                "specifically missing observable family before any basis growth"
            ),
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "st006_comparison_performed": False,
            "reason": "diagnostic-only observable Jacobian; no matched pressure/forcing child",
        },
        **TRUTH,
    }


def run(out: Path) -> dict[str, Any]:
    report = build_report()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.out)
    jac = report["observable_space_jacobian"]
    compact = {
        "observable_count": jac["observable_count"],
        "rank": jac["normalized_column_rank"],
        "condition": jac["normalized_condition"],
        "max_abs_pairwise_cosine": jac["max_abs_pairwise_cosine"],
        "capacity_gate_passes": jac["capacity_gate_passes"],
        "all_derivative_stability_gates_pass": jac[
            "all_derivative_stability_gates_pass"
        ],
        "channel_group_energy_fractions": {
            channel: row["observable_group_energy_fractions"]
            for channel, row in jac["channels"].items()
        },
        "decision": report["decision"],
    }
    print(json.dumps(compact, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
