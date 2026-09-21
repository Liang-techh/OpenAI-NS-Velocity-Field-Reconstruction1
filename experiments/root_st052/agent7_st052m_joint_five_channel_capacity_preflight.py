"""Joint capacity audit for the existing minimal ST052-M morphology controls.

Preregistered in issue #1031 and stacked exactly on Agent-7 PR #862 head.
This increment deliberately adds no basis and changes no selected candidate.
It asks one bounded question before any further basis growth: are the already
available low-dimensional morphology controls jointly identifiable on one
common frozen spacetime response cloud?

The five frozen tangent directions are:

1. live scalar amplitude of the existing outer-reservoir child;
2. recalibrated radial-envelope shape k from #810;
3. same-dimension axial-turnover coordinate s from #862;
4. endpoint-preserving temporal curvature g2 from #818;
5. compact toroidal/swirl channel from the #837/#855 lineage.

A PASS only means these controls are jointly nonredundant enough that blind
basis growth is not justified by capacity alone.  A FAIL would instead route
the next Agent-7 step toward reparameterizing the redundant directions before
adding dimensions.  Neither outcome changes velocity(x,y,z,t), establishes
visual correspondence, or evaluates the complete Navier-Stokes residual.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_axial_turnover_shape_preflight as axial
import agent7_st052m_compact_toroidal_swirl_preflight as toroidal
import agent7_st052m_endpoint_temporal_curvature_preflight as temporal

TASK_ID = "CR003-ST052M-JOINT-FIVE-CHANNEL-CAPACITY-123"
PREREG_ISSUE = 1031
SOURCE_PARENT_PR = 862
SOURCE_PARENT_HEAD = "f5dd27a47b03fc6239a2e6186cef78c4f3d3e024"

SOURCE_RADIAL_PR = 810
SOURCE_TEMPORAL_PR = 818
SOURCE_TOROIDAL_PR = 837
SOURCE_TOROIDAL_VORTICITY_PR = 855
SOURCE_AXIAL_PR = 862

CHANNELS = (
    "amplitude",
    "radial_shape",
    "axial_turnover",
    "temporal_curvature",
    "toroidal_swirl",
)
RANK_TARGET = 5
COSINE_MAX_ABS = 0.995
CONDITION_MAX = 25.0
TOROIDAL_LEAKAGE_MAX = 1.0e-12

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
    "visualization_fingerprint_sensitivities_recorded": True,
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
    axial._assert_source_lock()
    temporal._assert_source_lock()
    toroidal._assert_source_lock()
    if axial.TASK_ID != "CR003-ST052M-AXIAL-TURNOVER-SHAPE-PREFLIGHT-122":
        raise RuntimeError("#862 task identity drifted")
    if axial.PREREG_ISSUE != 861:
        raise RuntimeError("#862 preregistration drifted")
    if temporal.TASK_ID != "CR003-ST052M-ENDPOINT-TEMPORAL-CURVATURE-PREFLIGHT-116":
        raise RuntimeError("#818 task identity drifted")
    if temporal.PREREG_ISSUE != 817:
        raise RuntimeError("#818 preregistration drifted")
    if toroidal.TASK_ID != "CR003-ST052M-COMPACT-TOROIDAL-SWIRL-PREFLIGHT-119":
        raise RuntimeError("#837 task identity drifted")
    if toroidal.PREREG_ISSUE != 836:
        raise RuntimeError("#837 preregistration drifted")
    if tuple(axial.RESPONSE_TIMES) != tuple(temporal.shape_audit.RESPONSE_TIMES):
        raise RuntimeError("axial/temporal response-time semantics drifted")
    if tuple(axial.RESPONSE_TIMES) != tuple(toroidal.RESPONSE_TIMES):
        raise RuntimeError("axial/toroidal response-time semantics drifted")


def _normalize_column(values: np.ndarray, name: str) -> tuple[np.ndarray, float]:
    arr = np.asarray(values, dtype=float).ravel()
    norm = float(np.linalg.norm(arr))
    if not np.isfinite(norm) or norm <= np.finfo(float).tiny:
        raise RuntimeError(f"zero/nonfinite response column: {name}")
    return arr / norm, norm


def _build_tangent_blocks() -> tuple[np.ndarray, dict[str, list[np.ndarray]], dict[str, Any]]:
    """Evaluate all five frozen tangent directions on one exact response cloud."""
    _assert_source_lock()
    parent_fn, s_alphas, s_calibrations, s_children, k_alphas, k_children = (
        axial.build_local_families()
    )
    points = axial.radial_shape._response_points()
    live_alpha = float(s_alphas[axial.S_LIVE])
    toroidal_spatial = toroidal.toroidal_unit_velocity(points)

    blocks: dict[str, list[np.ndarray]] = {name: [] for name in CHANNELS}
    for time in axial.RESPONSE_TIMES:
        t = float(time)
        parent = np.asarray(parent_fn(points, t), dtype=float)
        live = np.asarray(s_children[axial.S_LIVE](points, t), dtype=float)
        k_minus = np.asarray(k_children[axial.radial_shape.K_MINUS](points, t), dtype=float)
        k_plus = np.asarray(k_children[axial.radial_shape.K_PLUS](points, t), dtype=float)
        s_minus = np.asarray(s_children[axial.S_MINUS](points, t), dtype=float)
        s_plus = np.asarray(s_children[axial.S_PLUS](points, t), dtype=float)

        blocks["amplitude"].append(live - parent)
        blocks["radial_shape"].append(
            (k_plus - k_minus) / (axial.radial_shape.K_PLUS - axial.radial_shape.K_MINUS)
        )
        blocks["axial_turnover"].append(
            (s_plus - s_minus) / (axial.S_PLUS - axial.S_MINUS)
        )
        blocks["temporal_curvature"].append(
            temporal.temporal_curvature_tangent(points, t, live_alpha)
        )
        blocks["toroidal_swirl"].append(toroidal.activation_g1(t) * toroidal_spatial)

    for name, channel_blocks in blocks.items():
        if len(channel_blocks) != len(axial.RESPONSE_TIMES):
            raise RuntimeError(f"incomplete response blocks for {name}")
        if any(np.asarray(block).shape != points.shape for block in channel_blocks):
            raise RuntimeError(f"response-shape drift for {name}")
        if not all(np.all(np.isfinite(block)) for block in channel_blocks):
            raise RuntimeError(f"nonfinite response for {name}")

    metadata = {
        "point_count_per_time": int(points.shape[0]),
        "response_times": [float(v) for v in axial.RESPONSE_TIMES],
        "live_alpha": live_alpha,
        "live_alpha_relative_mismatch_vs_exact_775": float(
            s_calibrations[axial.S_LIVE]["relative_mismatch_vs_exact_775"]
        ),
        "radial_alphas": {f"k={k:.2f}": float(v) for k, v in sorted(k_alphas.items())},
        "axial_alphas": {f"s={s:.2f}": float(v) for s, v in sorted(s_alphas.items())},
    }
    return points, blocks, metadata


def joint_identifiability() -> dict[str, Any]:
    """Rank/conditioning audit of all five already-existing morphology controls."""
    _, blocks, metadata = _build_tangent_blocks()
    normalized: dict[str, np.ndarray] = {}
    norms: dict[str, float] = {}
    for name in CHANNELS:
        vector = np.concatenate([np.asarray(v, dtype=float).ravel() for v in blocks[name]])
        normalized[name], norms[name] = _normalize_column(vector, name)

    matrix = np.column_stack([normalized[name] for name in CHANNELS])
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0.0 else float("inf")
    gram = matrix.T @ matrix

    cosines: dict[str, float] = {}
    for i, j in itertools.combinations(range(len(CHANNELS)), 2):
        cosines[f"{CHANNELS[i]}__{CHANNELS[j]}"] = float(gram[i, j])
    max_abs_cosine = float(max(abs(v) for v in cosines.values()))
    passed = bool(
        rank == RANK_TARGET
        and condition <= CONDITION_MAX
        and max_abs_cosine < COSINE_MAX_ABS
    )
    return {
        **metadata,
        "channel_order": list(CHANNELS),
        "flattened_component_count": int(matrix.shape[0]),
        "column_norms": norms,
        "normalized_column_rank": rank,
        "rank_target": RANK_TARGET,
        "normalized_condition": condition,
        "condition_gate": CONDITION_MAX,
        "pairwise_cosines": cosines,
        "max_abs_pairwise_cosine": max_abs_cosine,
        "cosine_abs_gate": COSINE_MAX_ABS,
        "singular_values": [float(v) for v in singular],
        "passes": passed,
    }


def visualization_fingerprint_sensitivities() -> dict[str, Any]:
    """Record target-free cylindrical/component and temporal sensitivity fingerprints."""
    points, blocks, _ = _build_tangent_blocks()
    repeated_points = np.vstack([points for _ in axial.RESPONSE_TIMES])
    result: dict[str, Any] = {}

    for name in CHANNELS:
        vectors = np.vstack(blocks[name])
        cylindrical = toroidal.cylindrical_components(repeated_points, vectors)
        component_rms = np.sqrt(np.mean(cylindrical * cylindrical, axis=0))
        total_rms = float(np.sqrt(np.mean(np.sum(cylindrical * cylindrical, axis=1))))
        per_time_vector_rms = {
            f"{float(time):.3f}": float(
                np.sqrt(np.mean(np.sum(np.asarray(block, dtype=float) ** 2, axis=1)))
            )
            for time, block in zip(axial.RESPONSE_TIMES, blocks[name])
        }
        result[name] = {
            "cylindrical_component_rms": {
                "u_r": float(component_rms[0]),
                "u_theta": float(component_rms[1]),
                "u_z": float(component_rms[2]),
            },
            "total_vector_rms": total_rms,
            "per_time_vector_rms": per_time_vector_rms,
            "normalized_component_energy_fractions": {
                "u_r": float(component_rms[0] ** 2 / max(total_rms**2, np.finfo(float).tiny)),
                "u_theta": float(
                    component_rms[1] ** 2 / max(total_rms**2, np.finfo(float).tiny)
                ),
                "u_z": float(component_rms[2] ** 2 / max(total_rms**2, np.finfo(float).tiny)),
            },
        }

    return {
        "autonomous_target_free_fingerprints": True,
        "direct_candidate_change": False,
        "direct_fingerprint_improvement": 0.0,
        "channels": result,
    }


def structural_guards() -> dict[str, Any]:
    """Recheck the two representation-specific invariants on the common cloud."""
    points = axial.radial_shape._response_points()
    live_alpha = float(axial.build_local_families()[1][axial.S_LIVE])

    temporal_start = temporal.temporal_curvature_tangent(points, 0.25, live_alpha)
    temporal_end = temporal.temporal_curvature_tangent(points, 0.75, live_alpha)
    temporal_endpoint_max = float(
        max(np.max(np.abs(temporal_start)), np.max(np.abs(temporal_end)))
    )

    toroidal_unit = toroidal.toroidal_unit_velocity(points)
    toroidal_cyl = toroidal.cylindrical_components(points, toroidal_unit)
    toroidal_rz_leakage = float(
        max(np.max(np.abs(toroidal_cyl[:, 0])), np.max(np.abs(toroidal_cyl[:, 2])))
    )
    toroidal_theta_rms = float(np.sqrt(np.mean(toroidal_cyl[:, 1] ** 2)))

    amplitude_mid = axial.unit_time_correction_s(points, 0.50, axial.S_LIVE)
    amplitude_cyl = toroidal.cylindrical_components(points, amplitude_mid)
    poloidal_theta_leakage = float(np.max(np.abs(amplitude_cyl[:, 1])))

    temporal_exact = bool(temporal_endpoint_max == 0.0)
    toroidal_selective = bool(
        toroidal_rz_leakage <= TOROIDAL_LEAKAGE_MAX and toroidal_theta_rms > 1.0e-8
    )
    poloidal_selective = bool(poloidal_theta_leakage <= TOROIDAL_LEAKAGE_MAX)
    return {
        "temporal_endpoint_max_abs": temporal_endpoint_max,
        "temporal_endpoints_exactly_zero": temporal_exact,
        "toroidal_rz_leakage_max_abs": toroidal_rz_leakage,
        "toroidal_u_theta_rms": toroidal_theta_rms,
        "poloidal_u_theta_leakage_max_abs": poloidal_theta_leakage,
        "leakage_gate": TOROIDAL_LEAKAGE_MAX,
        "toroidal_channel_selective": toroidal_selective,
        "poloidal_channel_selective": poloidal_selective,
        "all_pass": bool(temporal_exact and toroidal_selective and poloidal_selective),
    }


def build_report() -> dict[str, Any]:
    identifiability = joint_identifiability()
    fingerprints = visualization_fingerprint_sensitivities()
    guards = structural_guards()
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_assets": {
            "radial_shape_pr": SOURCE_RADIAL_PR,
            "temporal_curvature_pr": SOURCE_TEMPORAL_PR,
            "toroidal_swirl_pr": SOURCE_TOROIDAL_PR,
            "toroidal_vorticity_pr": SOURCE_TOROIDAL_VORTICITY_PR,
            "axial_turnover_pr": SOURCE_AXIAL_PR,
        },
        "frozen_protocol": {
            "channels": list(CHANNELS),
            "response_times": [float(v) for v in axial.RESPONSE_TIMES],
            "rank_target": RANK_TARGET,
            "condition_gate": CONDITION_MAX,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "toroidal_leakage_gate": TOROIDAL_LEAKAGE_MAX,
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "joint_identifiability": identifiability,
        "visualization_fingerprint_sensitivities": fingerprints,
        "structural_guards": guards,
        "decision": {
            "existing_low_dimensional_controls_jointly_identifiable": bool(
                identifiability["passes"] and guards["all_pass"]
            ),
            "additional_basis_dimension_justified_by_capacity_only": False,
            "if_pass_next_action": "wait_for_candidate_bound_frozen_ST052_morphology_discrepancy",
            "if_fail_next_action": "reparameterize_redundant_existing_channels_before_basis_growth",
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "reason": "capacity-only tangent audit; no selected child and no matched pressure/restricted forcing rebuild",
            "st006_comparison_performed": False,
        },
        **TRUTH,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
