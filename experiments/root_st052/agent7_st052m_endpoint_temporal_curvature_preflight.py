"""Preflight one endpoint-preserving temporal basis for the frozen ST052-M reservoir channel.

Frozen in issue #817 and stacked directly on exact Agent-7 PR #810 head.
This is one expression-capacity increment only: it changes no candidate velocity
and selects no temporal coefficient.

The live outer-reservoir child uses

    g1(t) = 2 * (t - 0.25)

on the registered interval [0.25, 0.75].  This audit adds only a hypothetical
response direction

    g2(t) = 16 * (t - 0.25) * (0.75 - t),

which vanishes at both registered endpoints and equals one at t=0.50.  The
question is whether that time-curvature direction is locally independent of the
already-audited scalar-amplitude and recalibrated radial-shape directions.

All morphology quantities are autonomous repository fingerprints.  No OpenAI
numeric time law, core radius, image target, or PDE claim is introduced here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_total_shape_coordinate_sensitivity as shape_audit

TASK_ID = "CR003-ST052M-ENDPOINT-TEMPORAL-CURVATURE-PREFLIGHT-116"
PREREG_ISSUE = 817
SOURCE_PARENT_PR = 810
SOURCE_PARENT_HEAD = "b077eaffabfe66ce0b15ac29601a82bfddeb4cf9"
SOURCE_FROZEN_CHILD_PR = 775
SOURCE_FROZEN_CHILD_HEAD = "93887a59729d22113badf2ae4dac2f7d868e2703"
RELATED_FRESH_HOLDOUT_PR = 714

TIME_START = 0.25
TIME_END = 0.75
TIME_MID = 0.50
IDENTIFIABILITY_COSINE_MAX_ABS = 0.995
IDENTIFIABILITY_CONDITION_MAX = 25.0

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "second_spatial_poloidal_basis_added": False,
    "temporal_basis_dimension_changed": False,
    "temporal_coefficient_selected": False,
    "temporal_basis_capacity_only": True,
    "parameter_scan_performed": False,
    "optimization_performed": False,
    "post_result_retuning_performed": False,
    "fresh_714_path_data_used": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "total_field_visualization_fingerprints_recorded": False,
    "temporal_tangent_visualization_fingerprints_recorded": True,
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
    shape_audit._assert_source_lock()
    if shape_audit.TASK_ID != "CR003-ST052M-TOTAL-SHAPE-COORDINATE-SENSITIVITY-115":
        raise RuntimeError("#810 task identity drifted")
    if shape_audit.PREREG_ISSUE != 809:
        raise RuntimeError("#810 preregistration identity drifted")
    if shape_audit.SOURCE_CHILD_PR != SOURCE_FROZEN_CHILD_PR:
        raise RuntimeError("frozen #775 child identity drifted")
    if shape_audit.SOURCE_CHILD_HEAD != SOURCE_FROZEN_CHILD_HEAD:
        raise RuntimeError("frozen #775 child head drifted")
    if tuple(shape_audit.K_VALUES) != (0.90, 1.00, 1.10):
        raise RuntimeError("#810 local radial-shape family drifted")
    if tuple(shape_audit.RESPONSE_TIMES) != (0.375, 0.50, 0.625, 0.75):
        raise RuntimeError("#810 response times drifted")
    if RELATED_FRESH_HOLDOUT_PR != 714:
        raise RuntimeError("fresh holdout identity drifted")


def temporal_curvature(time: float) -> float:
    """Frozen endpoint-preserving autonomous time shape g2(t)."""
    t = float(time)
    if t < TIME_START - 1.0e-15 or t > TIME_END + 1.0e-15:
        raise ValueError("time outside registered interval [0.25,0.75]")
    return float(16.0 * (t - TIME_START) * (TIME_END - t))


def temporal_curvature_tangent(points: np.ndarray, time: float, alpha_live: float) -> np.ndarray:
    """Unit coefficient tangent alpha_live*g2(t)*C_live(x)."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    spatial = shape_audit.reservoir_correction_k(pts, shape_audit.K_LIVE)
    out = float(alpha_live) * temporal_curvature(float(time)) * spatial
    if not np.all(np.isfinite(out)):
        raise RuntimeError("nonfinite temporal-curvature tangent")
    return out


def _normalize_column(values: np.ndarray) -> tuple[np.ndarray, float]:
    arr = np.asarray(values, dtype=float).ravel()
    norm = float(np.linalg.norm(arr))
    if not np.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("response column is zero or nonfinite")
    return arr / norm, norm


def three_column_identifiability() -> dict[str, Any]:
    """Audit amplitude, radial-shape and endpoint-curvature tangent independence."""
    _assert_source_lock()
    parent_fn, alphas, calibrations, children = shape_audit.build_local_family()
    live_alpha = float(alphas[shape_audit.K_LIVE])
    points = shape_audit._response_points()

    amplitude_blocks: list[np.ndarray] = []
    shape_blocks: list[np.ndarray] = []
    temporal_blocks: list[np.ndarray] = []
    for time in shape_audit.RESPONSE_TIMES:
        parent = np.asarray(parent_fn(points, float(time)), dtype=float)
        live = np.asarray(children[shape_audit.K_LIVE](points, float(time)), dtype=float)
        minus = np.asarray(children[shape_audit.K_MINUS](points, float(time)), dtype=float)
        plus = np.asarray(children[shape_audit.K_PLUS](points, float(time)), dtype=float)
        amplitude_blocks.append((live - parent).ravel())
        shape_blocks.append(
            ((plus - minus) / (shape_audit.K_PLUS - shape_audit.K_MINUS)).ravel()
        )
        temporal_blocks.append(
            temporal_curvature_tangent(points, float(time), live_alpha).ravel()
        )

    amplitude, amplitude_norm = _normalize_column(np.concatenate(amplitude_blocks))
    shape, shape_norm = _normalize_column(np.concatenate(shape_blocks))
    temporal, temporal_norm = _normalize_column(np.concatenate(temporal_blocks))
    matrix = np.column_stack((amplitude, shape, temporal))
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix))
    condition = float(singular[0] / singular[-1]) if singular[-1] > 0.0 else float("inf")
    gram = matrix.T @ matrix
    cosines = {
        "amplitude_shape": float(gram[0, 1]),
        "amplitude_temporal": float(gram[0, 2]),
        "shape_temporal": float(gram[1, 2]),
    }
    max_abs_cosine = float(max(abs(v) for v in cosines.values()))
    passes = bool(
        rank == 3
        and condition <= IDENTIFIABILITY_CONDITION_MAX
        and max_abs_cosine < IDENTIFIABILITY_COSINE_MAX_ABS
    )
    return {
        "live_alpha": live_alpha,
        "live_alpha_replay_relative_mismatch": float(
            calibrations[shape_audit.K_LIVE]["relative_mismatch_vs_exact_775"]
        ),
        "response_point_count": int(points.shape[0]),
        "response_times": [float(v) for v in shape_audit.RESPONSE_TIMES],
        "column_norms": {
            "amplitude": amplitude_norm,
            "radial_shape": shape_norm,
            "temporal_curvature": temporal_norm,
        },
        "normalized_column_rank": rank,
        "normalized_three_column_condition": condition,
        "normalized_pairwise_cosines": cosines,
        "max_abs_pairwise_cosine": max_abs_cosine,
        "singular_values": [float(v) for v in singular],
        "passes": passes,
    }


def _temporal_tangent_field(alpha_live: float) -> Callable[[np.ndarray, float], np.ndarray]:
    return lambda pts, time: temporal_curvature_tangent(pts, float(time), alpha_live)


def temporal_morphology_fingerprints(alpha_live: float) -> dict[str, Any]:
    """Record autonomous visible-tip / axial-stretch response of the new time tangent."""
    _assert_source_lock()
    tangent_fn = _temporal_tangent_field(float(alpha_live))
    by_time: dict[str, Any] = {}

    axial_radii = np.asarray(shape_audit.total_axial.RADII, dtype=float)
    core_mask = axial_radii <= 0.90 + 1.0e-15
    annulus_mask = (axial_radii >= 1.20 - 1.0e-15) & (axial_radii <= 1.50 + 1.0e-15)
    if not np.any(core_mask) or not np.any(annulus_mask):
        raise RuntimeError("registered axial-stretch masks are empty")

    for time in (0.25, 0.375, 0.50, 0.625, 0.75):
        radial_bands: dict[str, Any] = {}
        axial_bands: dict[str, Any] = {}
        for abs_z in (1.55, 1.75):
            radial_points = shape_audit.replay._probe_points(float(abs_z))
            radial = shape_audit.replay.screen.loc._radial_velocity(
                tangent_fn, radial_points, float(time)
            )
            radial_bands[f"|z|={abs_z:.2f}"] = {
                "mean_radial_tangent": float(np.mean(radial)),
                "rms_radial_tangent": float(np.sqrt(np.mean(np.asarray(radial, dtype=float) ** 2))),
                "max_abs_radial_tangent": float(np.max(np.abs(radial))),
            }

            profile = np.asarray(
                shape_audit.total_axial.axial_stretch_profile(
                    tangent_fn, float(abs_z), float(time)
                ),
                dtype=float,
            )
            axial_bands[f"|z|={abs_z:.2f}"] = {
                "mean_core_away_from_midplane_axial_tangent": float(np.mean(profile[core_mask])),
                "mean_annular_away_from_midplane_axial_tangent": float(np.mean(profile[annulus_mask])),
                "max_abs_profile_tangent": float(np.max(np.abs(profile))),
            }

        by_time[f"{time:.3f}"] = {
            "g2": temporal_curvature(float(time)),
            "visible_radial": radial_bands,
            "axial_stretch": axial_bands,
        }

    endpoint_points = shape_audit._response_points()
    endpoint_zero = {
        "t=0.25_max_abs": float(
            np.max(np.abs(temporal_curvature_tangent(endpoint_points, 0.25, float(alpha_live))))
        ),
        "t=0.75_max_abs": float(
            np.max(np.abs(temporal_curvature_tangent(endpoint_points, 0.75, float(alpha_live))))
        ),
    }
    return {
        "time_profile": {
            f"{time:.3f}": temporal_curvature(float(time))
            for time in (0.25, 0.375, 0.50, 0.625, 0.75)
        },
        "endpoint_zero": endpoint_zero,
        "by_time": by_time,
        "fingerprints_are_autonomous_not_openai_numeric_targets": True,
    }


def build_report() -> dict[str, Any]:
    identifiability = three_column_identifiability()
    morphology = temporal_morphology_fingerprints(float(identifiability["live_alpha"]))
    endpoint_exact = bool(
        morphology["endpoint_zero"]["t=0.25_max_abs"] == 0.0
        and morphology["endpoint_zero"]["t=0.75_max_abs"] == 0.0
    )
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "source_frozen_child": {
            "pr": SOURCE_FROZEN_CHILD_PR,
            "head": SOURCE_FROZEN_CHILD_HEAD,
        },
        "time_basis": {
            "live": "g1(t)=2*(t-0.25)",
            "new_capacity_direction": "g2(t)=16*(t-0.25)*(0.75-t)",
            "registered_interval": [TIME_START, TIME_END],
            "g2_start": temporal_curvature(TIME_START),
            "g2_mid": temporal_curvature(TIME_MID),
            "g2_end": temporal_curvature(TIME_END),
            "coefficient_selected": None,
            "public_openai_numeric_time_target": None,
        },
        "three_column_identifiability": identifiability,
        "temporal_morphology_fingerprints": morphology,
        "decision": {
            "endpoint_preserving_temporal_curvature_direction_identifiable": bool(
                identifiability["passes"]
            ),
            "endpoint_field_exactly_preserved_by_unit_time_tangent": endpoint_exact,
            "temporal_coefficient_selected": False,
            "second_spatial_poloidal_basis_justified_by_this_increment": False,
            "actual_velocity_changed": False,
            "closer_visualization_delivery_established": False,
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
