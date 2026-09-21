"""Map the five existing ST052-M morphology controls into signed defect coordinates.

Preregistered in issue #1077 and stacked exactly on Agent-7 PR #1069 head.
This is a diagnostic/routing increment only: the selected candidate, basis
size, coefficients, pressure, forcing and scientific validation state are not
changed.  The coordinates below are autonomous renderer-independent summaries
of observables already frozen in #1069; they are not numerical targets inferred
from the public OpenAI visualization.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_observable_space_five_channel_jacobian as parent

TASK_ID = "CR003-ST052M-DEFECT-COORDINATE-CONTROLLABILITY-128"
PREREG_ISSUE = 1077
SOURCE_PARENT_PR = 1069
SOURCE_PARENT_HEAD = "a564fbc8bfa6a5c630b5839fa208e41cb7bc72d7"

SUPPORT_RADIUS = 2.0
EPSILONS = (1.0e-3, 5.0e-4)
PRIMARY_EPSILON = 5.0e-4
RANK_TARGET = 5
CONDITION_MAX = 25.0
COSINE_MAX_ABS = 0.995
DERIVATIVE_DRIFT_MAX = 0.05
DERIVATIVE_COSINE_MIN = 0.999
ROW_RESPONSE_FLOOR = 1.0e-12

DEFECT_COORDINATES = (
    "inward_radial_index_t050",
    "swirl_to_poloidal_ratio_t050",
    "axial_vorticity_aspect_t050",
    "core_width_t050",
    "tip_radial_thickness_t050",
    "angular_rotation_radial_variation_t050",
    "core_width_change_025_to_075",
    "core_speed_change_025_to_075",
)

COORDINATE_SEMANTICS = {
    "inward_radial_index_t050": {
        "positive_direction": "stronger inward radial motion",
        "public_observable_relation": "inward-spiraling trajectories",
    },
    "swirl_to_poloidal_ratio_t050": {
        "positive_direction": "more swirl relative to radial/axial motion",
        "public_observable_relation": "vortex/swirl presence; autonomous ratio diagnostic",
    },
    "axial_vorticity_aspect_t050": {
        "positive_direction": "more axially elongated enstrophy distribution",
        "public_observable_relation": "axial stretching/elongation",
    },
    "core_width_t050": {
        "positive_direction": "thicker speed-weighted midplane core",
        "public_observable_relation": "shrinking central region while speed increases",
    },
    "tip_radial_thickness_t050": {
        "positive_direction": "wider/blunter smooth-tip enstrophy region",
        "public_observable_relation": "autonomous tip-morphology diagnostic",
    },
    "angular_rotation_radial_variation_t050": {
        "positive_direction": "more spatial/radial variation in angular rotation",
        "public_observable_relation": "spatial variation in angular rotation; radius-dependent circulation",
    },
    "core_width_change_025_to_075": {
        "positive_direction": "core widens from start to end; negative means contraction",
        "public_observable_relation": "shrinking central region while speed increases",
    },
    "core_speed_change_025_to_075": {
        "positive_direction": "core speed increases from start to end",
        "public_observable_relation": "shrinking central region while speed increases",
    },
}

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
    "defect_coordinate_sensitivities_recorded": True,
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


def _assert_source_lock() -> None:
    parent._assert_source_lock()
    if parent.TASK_ID != "CR003-ST052M-OBSERVABLE-SPACE-FIVE-CHANNEL-JACOBIAN-127":
        raise RuntimeError("#1069 parent task identity drifted")
    if parent.PREREG_ISSUE != 1067:
        raise RuntimeError("#1069 parent preregistration identity drifted")
    if parent.SOURCE_PARENT_PR != 1060:
        raise RuntimeError("#1069 lineage drifted")
    if tuple(parent.EPSILONS) != EPSILONS or parent.PRIMARY_EPSILON != PRIMARY_EPSILON:
        raise RuntimeError("#1069 finite-difference protocol drifted")
    if tuple(parent._joint().CHANNELS) != (
        "amplitude",
        "radial_shape",
        "axial_turnover",
        "temporal_curvature",
        "toroidal_swirl",
    ):
        raise RuntimeError("five-channel order drifted")


def _ring_values(observables: dict[str, float], group: str, time: float) -> np.ndarray:
    prefix = f"{group}/t={time:.3f}/"
    values = np.asarray(
        [float(value) for name, value in observables.items() if name.startswith(prefix)],
        dtype=float,
    )
    if values.size != 16:
        raise RuntimeError(f"expected 16 {group} ring summaries at t={time:.3f}, got {values.size}")
    if not np.all(np.isfinite(values)):
        raise RuntimeError(f"nonfinite {group} ring summaries")
    return values


def _ring_radius_from_name(name: str) -> float:
    marker = "/r="
    if marker not in name:
        raise RuntimeError(f"missing radius in observable name: {name}")
    return float(name.rsplit(marker, 1)[1])


def _angular_rotation_variation(
    observables: dict[str, float], baseline_speed_scale: float
) -> float:
    prefix = "swirl_ring/t=0.500/"
    rows = [
        (name, float(value))
        for name, value in observables.items()
        if name.startswith(prefix)
    ]
    if len(rows) != 16:
        raise RuntimeError("unexpected t=.50 swirl-ring count")
    omega_star = np.asarray(
        [
            SUPPORT_RADIUS * value / (_ring_radius_from_name(name) * baseline_speed_scale)
            for name, value in rows
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(omega_star)):
        raise RuntimeError("nonfinite dimensionless angular-speed samples")
    return float(np.std(omega_star, ddof=0))


def _core_width(observables: dict[str, float], time: float) -> float:
    return float(
        observables[f"core_speed_width/t={time:.3f}/speed_weighted_rms_radius"]
    )


def _mean_core_probe_speed(observables: dict[str, float], time: float) -> float:
    return float(
        np.mean(
            [
                observables[f"core_speed_width/t={time:.3f}/mean_speed_r={radius:.2f}"]
                for radius in parent.CORE_SPEED_RADII
            ]
        )
    )


def defect_coordinate_vector(
    field_fn: Callable, baseline_speed_scale: float
) -> dict[str, float]:
    """Return the frozen target-free signed morphology defect coordinates."""
    if not np.isfinite(baseline_speed_scale) or baseline_speed_scale <= np.finfo(float).tiny:
        raise RuntimeError("invalid frozen baseline RMS-speed scale")
    obs = parent.observable_vector(field_fn)

    radial = _ring_values(obs, "radial_ring", 0.50)
    swirl = _ring_values(obs, "swirl_ring", 0.50)
    axial = _ring_values(obs, "axial_ring", 0.50)
    poloidal_rms = float(np.sqrt(np.mean(radial * radial + axial * axial)))
    swirl_rms = float(np.sqrt(np.mean(swirl * swirl)))
    if poloidal_rms <= np.finfo(float).tiny:
        raise RuntimeError("degenerate poloidal-ring RMS")

    width_start = _core_width(obs, 0.25)
    width_mid = _core_width(obs, 0.50)
    width_end = _core_width(obs, 0.75)
    speed_start = _mean_core_probe_speed(obs, 0.25)
    speed_end = _mean_core_probe_speed(obs, 0.75)

    out = {
        "inward_radial_index_t050": float(-np.mean(radial) / baseline_speed_scale),
        "swirl_to_poloidal_ratio_t050": float(swirl_rms / poloidal_rms),
        "axial_vorticity_aspect_t050": float(
            obs["vorticity_shape/t=0.500/full_aspect_ratio"]
        ),
        "core_width_t050": float(width_mid / SUPPORT_RADIUS),
        "tip_radial_thickness_t050": float(
            obs["vorticity_shape/t=0.500/smooth_tip_radial_rms"] / SUPPORT_RADIUS
        ),
        "angular_rotation_radial_variation_t050": _angular_rotation_variation(
            obs, baseline_speed_scale
        ),
        "core_width_change_025_to_075": float(
            (width_end - width_start) / SUPPORT_RADIUS
        ),
        "core_speed_change_025_to_075": float(
            (speed_end - speed_start) / baseline_speed_scale
        ),
    }
    if tuple(out) != DEFECT_COORDINATES:
        raise RuntimeError("defect-coordinate ordering drifted")
    if not all(np.isfinite(value) for value in out.values()):
        raise RuntimeError("nonfinite defect coordinate")
    return out


def _coordinate_derivative(
    base_fn: Callable,
    tangent_fn: Callable,
    epsilon: float,
    baseline_speed_scale: float,
) -> np.ndarray:
    plus = defect_coordinate_vector(
        parent._perturbed_field(base_fn, tangent_fn, float(epsilon)),
        baseline_speed_scale,
    )
    minus = defect_coordinate_vector(
        parent._perturbed_field(base_fn, tangent_fn, -float(epsilon)),
        baseline_speed_scale,
    )
    if tuple(plus) != DEFECT_COORDINATES or tuple(minus) != DEFECT_COORDINATES:
        raise RuntimeError("defect-coordinate ordering drifted under perturbation")
    derivative = np.asarray(
        [
            (plus[name] - minus[name]) / (2.0 * float(epsilon))
            for name in DEFECT_COORDINATES
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(derivative)):
        raise RuntimeError("nonfinite defect-coordinate derivative")
    return derivative


def controllability_audit() -> dict[str, Any]:
    _assert_source_lock()
    base_fn, tangents, metadata = parent._build_field_and_tangents()
    baseline_speed = parent._baseline_speed_scale(base_fn)
    baseline = defect_coordinate_vector(base_fn, baseline_speed)

    derivatives: dict[float, dict[str, np.ndarray]] = {}
    for epsilon in EPSILONS:
        derivatives[float(epsilon)] = {
            channel: _coordinate_derivative(
                base_fn,
                tangents[channel],
                float(epsilon),
                baseline_speed,
            )
            for channel in parent._joint().CHANNELS
        }

    fine = derivatives[PRIMARY_EPSILON]
    coarse = derivatives[EPSILONS[0]]
    normalized_columns: list[np.ndarray] = []
    column_norms: dict[str, float] = {}
    stability: dict[str, Any] = {}
    for channel in parent._joint().CHANNELS:
        vf = fine[channel]
        vc = coarse[channel]
        nf = float(np.linalg.norm(vf))
        nc = float(np.linalg.norm(vc))
        column_norms[channel] = nf
        if nf <= np.finfo(float).tiny:
            normalized_columns.append(np.zeros_like(vf))
            drift = None
            cosine = None
            stable = False
        else:
            normalized_columns.append(vf / nf)
            drift = float(np.linalg.norm(vf - vc) / nf)
            cosine = (
                None
                if nc <= np.finfo(float).tiny
                else float(np.dot(vf, vc) / (nf * nc))
            )
            stable = bool(
                drift <= DERIVATIVE_DRIFT_MAX
                and cosine is not None
                and cosine >= DERIVATIVE_COSINE_MIN
            )
        stability[channel] = {
            "coarse_epsilon": EPSILONS[0],
            "fine_epsilon": PRIMARY_EPSILON,
            "coarse_derivative_l2": nc,
            "fine_derivative_l2": nf,
            "relative_column_drift": drift,
            "coarse_fine_cosine": cosine,
            "drift_gate": DERIVATIVE_DRIFT_MAX,
            "cosine_gate": DERIVATIVE_COSINE_MIN,
            "passes": stable,
        }

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
    channels = tuple(parent._joint().CHANNELS)
    for i, j in itertools.combinations(range(len(channels)), 2):
        key = f"{channels[i]}__{channels[j]}"
        cosines[key] = (
            None
            if min(column_norms[channels[i]], column_norms[channels[j]])
            <= np.finfo(float).tiny
            else float(gram[i, j])
        )
    finite_cosines = [abs(value) for value in cosines.values() if value is not None]
    max_abs_cosine = None if not finite_cosines else float(max(finite_cosines))

    coordinate_rows: dict[str, Any] = {}
    for row_index, name in enumerate(DEFECT_COORDINATES):
        raw = np.asarray([fine[channel][row_index] for channel in channels], dtype=float)
        raw_norm = float(np.linalg.norm(raw))
        aligned = np.asarray(matrix[row_index, :], dtype=float)
        aligned_power = aligned * aligned
        aligned_total = float(np.sum(aligned_power))
        if aligned_total <= np.finfo(float).tiny:
            fractions = {channel: 0.0 for channel in channels}
            dominant = None
        else:
            fractions = {
                channel: float(aligned_power[index] / aligned_total)
                for index, channel in enumerate(channels)
            }
            dominant = channels[int(np.argmax(np.abs(aligned)))]
        coordinate_rows[name] = {
            **COORDINATE_SEMANTICS[name],
            "baseline_value": float(baseline[name]),
            "fine_raw_derivative_by_channel": {
                channel: float(raw[index]) for index, channel in enumerate(channels)
            },
            "fine_raw_row_l2": raw_norm,
            "structurally_unresponsive_at_numeric_floor": bool(
                raw_norm <= ROW_RESPONSE_FLOOR
            ),
            "normalized_column_alignment_by_channel": {
                channel: float(aligned[index]) for index, channel in enumerate(channels)
            },
            "normalized_alignment_energy_fraction_by_channel": fractions,
            "dominant_aligned_existing_channel": dominant,
        }

    all_rows_respond = bool(
        all(
            not row["structurally_unresponsive_at_numeric_floor"]
            for row in coordinate_rows.values()
        )
    )
    geometry_pass = bool(
        rank == RANK_TARGET
        and condition is not None
        and condition <= CONDITION_MAX
        and max_abs_cosine is not None
        and max_abs_cosine < COSINE_MAX_ABS
    )
    stability_pass = bool(all(row["passes"] for row in stability.values()))

    return {
        **metadata,
        "defect_coordinate_count": len(DEFECT_COORDINATES),
        "defect_coordinate_order": list(DEFECT_COORDINATES),
        "baseline_rms_speed_scale": float(baseline_speed),
        "support_radius_scale": SUPPORT_RADIUS,
        "baseline_defect_coordinates": {
            name: float(baseline[name]) for name in DEFECT_COORDINATES
        },
        "column_norms_in_defect_coordinate_space": column_norms,
        "normalized_column_rank": rank,
        "rank_target": RANK_TARGET,
        "normalized_condition": condition,
        "condition_gate": CONDITION_MAX,
        "pairwise_cosines": cosines,
        "max_abs_pairwise_cosine": max_abs_cosine,
        "cosine_abs_gate": COSINE_MAX_ABS,
        "singular_values": [float(value) for value in singular],
        "compressed_routing_geometry_gate_passes": geometry_pass,
        "derivative_stability": stability,
        "all_derivative_stability_gates_pass": stability_pass,
        "all_coordinates_have_numeric_response": all_rows_respond,
        "numeric_row_response_floor": ROW_RESPONSE_FLOOR,
        "coordinates": coordinate_rows,
    }


def build_report() -> dict[str, Any]:
    audit = controllability_audit()
    scientific_routing_pass = bool(
        audit["compressed_routing_geometry_gate_passes"]
        and audit["all_derivative_stability_gates_pass"]
        and audit["all_coordinates_have_numeric_response"]
    )
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "frozen_protocol": {
            "channels": list(parent._joint().CHANNELS),
            "defect_coordinates": list(DEFECT_COORDINATES),
            "support_radius_scale": SUPPORT_RADIUS,
            "epsilons": list(EPSILONS),
            "primary_epsilon": PRIMARY_EPSILON,
            "rank_target": RANK_TARGET,
            "condition_gate": CONDITION_MAX,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "derivative_drift_gate": DERIVATIVE_DRIFT_MAX,
            "derivative_cosine_gate": DERIVATIVE_COSINE_MIN,
            "numeric_row_response_floor": ROW_RESPONSE_FLOOR,
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "defect_coordinate_controllability": audit,
        "decision": {
            "compressed_defect_routing_and_linearization_pass": scientific_routing_pass,
            "additional_basis_dimension_justified_by_this_audit": False,
            "candidate_mutation_authorized_by_this_audit": False,
            "if_coordinate_has_response": (
                "wait for an actual candidate-bound public-observable discrepancy, then route "
                "that named defect to the smallest existing aligned channel"
            ),
            "if_coordinate_is_unresponsive_or_geometry_is_ill_conditioned": (
                "record the structural routing obstruction; only after an actual discrepancy "
                "selects that defect may one minimal reparameterization or targeted basis be tested"
            ),
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "st006_comparison_performed": False,
            "reason": "diagnostic-only defect-coordinate Jacobian; no candidate change",
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
    audit = report["defect_coordinate_controllability"]
    compact = {
        "defect_coordinate_count": audit["defect_coordinate_count"],
        "rank": audit["normalized_column_rank"],
        "condition": audit["normalized_condition"],
        "max_abs_pairwise_cosine": audit["max_abs_pairwise_cosine"],
        "routing_geometry_gate_passes": audit["compressed_routing_geometry_gate_passes"],
        "derivative_stability_passes": audit["all_derivative_stability_gates_pass"],
        "all_coordinates_have_numeric_response": audit[
            "all_coordinates_have_numeric_response"
        ],
        "dominant_existing_channel_by_coordinate": {
            name: row["dominant_aligned_existing_channel"]
            for name, row in audit["coordinates"].items()
        },
        "decision": report["decision"],
    }
    print(json.dumps(compact, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
