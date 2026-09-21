"""Regional localization of the existing five-channel ST052-M morphology capacity.

Preregistered in issue #1041 and stacked exactly on Agent-7 PR #1033 head
c692a445e4875bcc8965d649804a3b1a33a0ecec.

This increment deliberately adds no basis and changes no candidate.  It reuses
exactly the five frozen tangent directions and response cloud from #1033, then
asks a narrower representation question: where in space does each existing
control have leverage, and does any broad morphology-relevant region expose a
local rank/conditioning obstruction?

The regions are target-free sensitivity localizers derived only from medians of
the already-frozen response cloud.  They are not numerical targets inferred
from the public OpenAI visualization.  A regional failure routes future work
toward reparameterizing locally redundant existing directions; it does not
authorize blind basis growth.  A pass likewise does not establish visual
correspondence or PDE validity.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_joint_five_channel_capacity_preflight as joint

TASK_ID = "CR003-ST052M-REGIONAL-FIVE-CHANNEL-CAPACITY-124"
PREREG_ISSUE = 1041
SOURCE_PARENT_PR = 1033
SOURCE_PARENT_HEAD = "c692a445e4875bcc8965d649804a3b1a33a0ecec"

REGIONS = ("inner_core", "outer_radial", "midplane", "axial_tip")
MIN_REGION_FRACTION = 0.25
RANK_TARGET = joint.RANK_TARGET
CONDITION_MAX = joint.CONDITION_MAX
COSINE_MAX_ABS = joint.COSINE_MAX_ABS
LEAKAGE_MAX = joint.TOROIDAL_LEAKAGE_MAX

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
    "regional_visualization_sensitivity_localized": True,
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
    joint._assert_source_lock()
    if joint.TASK_ID != "CR003-ST052M-JOINT-FIVE-CHANNEL-CAPACITY-123":
        raise RuntimeError("#1033 task identity drifted")
    if joint.PREREG_ISSUE != 1031:
        raise RuntimeError("#1033 preregistration drifted")
    if tuple(joint.CHANNELS) != (
        "amplitude",
        "radial_shape",
        "axial_turnover",
        "temporal_curvature",
        "toroidal_swirl",
    ):
        raise RuntimeError("#1033 channel set drifted")
    if joint.RANK_TARGET != 5 or joint.CONDITION_MAX != 25.0:
        raise RuntimeError("#1033 rank/conditioning gate drifted")
    if joint.COSINE_MAX_ABS != 0.995:
        raise RuntimeError("#1033 cosine gate drifted")


def _region_masks(points: np.ndarray) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    """Build the four preregistered median-derived spatial masks."""
    xyz = np.asarray(points, dtype=float)
    if xyz.ndim != 2 or xyz.shape[1] != 3 or xyz.shape[0] < 4:
        raise RuntimeError("invalid frozen response cloud")
    if not np.all(np.isfinite(xyz)):
        raise RuntimeError("nonfinite frozen response cloud")

    radius = np.hypot(xyz[:, 0], xyz[:, 1])
    abs_z = np.abs(xyz[:, 2])
    r50 = float(np.median(radius))
    z50 = float(np.median(abs_z))

    masks = {
        "inner_core": radius <= r50,
        "outer_radial": radius >= r50,
        "midplane": abs_z <= z50,
        "axial_tip": abs_z >= z50,
    }
    n = int(xyz.shape[0])
    for name in REGIONS:
        mask = np.asarray(masks[name], dtype=bool)
        if mask.shape != (n,):
            raise RuntimeError(f"region mask shape drift: {name}")
        if float(np.mean(mask)) < MIN_REGION_FRACTION:
            raise RuntimeError(f"region too small under frozen protocol: {name}")

    return masks, {"r50": r50, "abs_z50": z50}


def _flatten_channel(blocks: list[np.ndarray], mask: np.ndarray | None = None) -> np.ndarray:
    selected: list[np.ndarray] = []
    for block in blocks:
        arr = np.asarray(block, dtype=float)
        selected.append(arr if mask is None else arr[mask])
    return np.concatenate([arr.ravel() for arr in selected])


def _safe_normalized_columns(
    region_vectors: dict[str, np.ndarray],
) -> tuple[np.ndarray, dict[str, float], dict[str, bool]]:
    columns: list[np.ndarray] = []
    norms: dict[str, float] = {}
    nonzero: dict[str, bool] = {}
    for name in joint.CHANNELS:
        vector = np.asarray(region_vectors[name], dtype=float).ravel()
        norm = float(np.linalg.norm(vector))
        if not np.isfinite(norm):
            raise RuntimeError(f"nonfinite regional norm: {name}")
        norms[name] = norm
        is_nonzero = bool(norm > np.finfo(float).tiny)
        nonzero[name] = is_nonzero
        columns.append(vector / norm if is_nonzero else np.zeros_like(vector))
    return np.column_stack(columns), norms, nonzero


def _regional_component_fingerprint(
    points: np.ndarray,
    blocks: dict[str, list[np.ndarray]],
    mask: np.ndarray,
) -> dict[str, Any]:
    region_points = np.asarray(points, dtype=float)[mask]
    repeated_points = np.vstack([region_points for _ in joint.axial.RESPONSE_TIMES])
    out: dict[str, Any] = {}
    labels = ("u_r", "u_theta", "u_z")

    for name in joint.CHANNELS:
        vectors = np.vstack([np.asarray(block, dtype=float)[mask] for block in blocks[name]])
        cylindrical = joint.toroidal.cylindrical_components(repeated_points, vectors)
        rms = np.sqrt(np.mean(cylindrical * cylindrical, axis=0))
        total = float(np.sqrt(np.mean(np.sum(cylindrical * cylindrical, axis=1))))
        denom = max(total * total, np.finfo(float).tiny)
        fractions = {
            label: float(value * value / denom) for label, value in zip(labels, rms)
        }
        dominant = labels[int(np.argmax(rms))]
        out[name] = {
            "cylindrical_component_rms": {
                label: float(value) for label, value in zip(labels, rms)
            },
            "normalized_component_energy_fractions": fractions,
            "dominant_component": dominant,
            "total_vector_rms": total,
        }
    return out


def _audit_one_region(
    points: np.ndarray,
    blocks: dict[str, list[np.ndarray]],
    mask: np.ndarray,
    full_norms: dict[str, float],
) -> dict[str, Any]:
    region_vectors = {name: _flatten_channel(blocks[name], mask) for name in joint.CHANNELS}
    matrix, raw_norms, nonzero = _safe_normalized_columns(region_vectors)
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    condition: float | None
    if singular.size == 0 or singular[-1] <= np.finfo(float).tiny:
        condition = None
    else:
        condition = float(singular[0] / singular[-1])

    gram = matrix.T @ matrix
    cosines: dict[str, float | None] = {}
    finite_cosines: list[float] = []
    for i, j in itertools.combinations(range(len(joint.CHANNELS)), 2):
        left = joint.CHANNELS[i]
        right = joint.CHANNELS[j]
        key = f"{left}__{right}"
        if nonzero[left] and nonzero[right]:
            value = float(gram[i, j])
            cosines[key] = value
            finite_cosines.append(abs(value))
        else:
            cosines[key] = None
    max_abs_cosine = float(max(finite_cosines)) if finite_cosines else None

    energy_fraction = {
        name: float((raw_norms[name] / full_norms[name]) ** 2)
        if full_norms[name] > np.finfo(float).tiny
        else 0.0
        for name in joint.CHANNELS
    }
    top_two = sorted(joint.CHANNELS, key=lambda name: raw_norms[name], reverse=True)[:2]
    all_nonzero = bool(all(nonzero.values()))
    passed = bool(
        all_nonzero
        and rank == RANK_TARGET
        and condition is not None
        and condition <= CONDITION_MAX
        and max_abs_cosine is not None
        and max_abs_cosine < COSINE_MAX_ABS
    )

    fingerprints = _regional_component_fingerprint(points, blocks, mask)
    toroidal = fingerprints["toroidal_swirl"]["cylindrical_component_rms"]
    toroidal_rz_leak = float(max(abs(toroidal["u_r"]), abs(toroidal["u_z"])))
    poloidal_theta_leak = float(
        max(
            abs(fingerprints[name]["cylindrical_component_rms"]["u_theta"])
            for name in joint.CHANNELS
            if name != "toroidal_swirl"
        )
    )

    return {
        "point_count": int(np.count_nonzero(mask)),
        "point_fraction": float(np.mean(mask)),
        "flattened_component_count": int(matrix.shape[0]),
        "raw_column_norms": raw_norms,
        "channel_nonzero": nonzero,
        "regional_energy_fraction_of_full_channel": energy_fraction,
        "normalized_column_rank": rank,
        "rank_target": RANK_TARGET,
        "normalized_condition": condition,
        "condition_gate": CONDITION_MAX,
        "pairwise_cosines": cosines,
        "max_abs_pairwise_cosine": max_abs_cosine,
        "cosine_abs_gate": COSINE_MAX_ABS,
        "singular_values": [float(value) for value in singular],
        "top_two_channels_by_raw_norm": list(top_two),
        "component_fingerprints": fingerprints,
        "toroidal_rz_leakage_rms_max": toroidal_rz_leak,
        "poloidal_theta_leakage_rms_max": poloidal_theta_leak,
        "selectivity_guard_pass": bool(
            toroidal_rz_leak <= LEAKAGE_MAX and poloidal_theta_leak <= LEAKAGE_MAX
        ),
        "passes": passed,
    }


def regional_capacity_localization() -> dict[str, Any]:
    """Localize joint rank, conditioning and component leverage in four broad regions."""
    _assert_source_lock()
    points, blocks, metadata = joint._build_tangent_blocks()
    masks, thresholds = _region_masks(points)
    full_norms = {
        name: float(np.linalg.norm(_flatten_channel(blocks[name]))) for name in joint.CHANNELS
    }
    if any((not np.isfinite(value)) or value <= np.finfo(float).tiny for value in full_norms.values()):
        raise RuntimeError("zero/nonfinite full-cloud channel norm")

    regions = {
        name: _audit_one_region(points, blocks, masks[name], full_norms) for name in REGIONS
    }
    failing = [name for name in REGIONS if not regions[name]["passes"]]
    selectivity_failures = [
        name for name in REGIONS if not regions[name]["selectivity_guard_pass"]
    ]
    return {
        **metadata,
        "region_order": list(REGIONS),
        "median_thresholds": thresholds,
        "minimum_region_fraction": MIN_REGION_FRACTION,
        "full_cloud_channel_norms": full_norms,
        "regions": regions,
        "all_regions_pass": bool(not failing),
        "failing_regions": failing,
        "all_region_selectivity_guards_pass": bool(not selectivity_failures),
        "selectivity_guard_failures": selectivity_failures,
    }


def build_report() -> dict[str, Any]:
    localization = regional_capacity_localization()
    parent_guards = joint.structural_guards()
    all_mechanical_guards = bool(
        parent_guards["all_pass"] and localization["all_region_selectivity_guards_pass"]
    )
    if localization["all_regions_pass"]:
        next_action = (
            "wait_for_candidate_bound_frozen_ST052_morphology_discrepancy_and_choose_"
            "the_smallest_existing_control"
        )
    else:
        next_action = (
            "inspect_or_reparameterize_locally_redundant_existing_channels_before_"
            "any_basis_growth"
        )

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "frozen_protocol": {
            "channels": list(joint.CHANNELS),
            "regions": list(REGIONS),
            "region_rule": "median-derived masks on exact #1033 response cloud",
            "minimum_region_fraction": MIN_REGION_FRACTION,
            "rank_target": RANK_TARGET,
            "condition_gate": CONDITION_MAX,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "regional_capacity_localization": localization,
        "parent_structural_guards": parent_guards,
        "mechanical_guards_pass": all_mechanical_guards,
        "decision": {
            "local_representation_obstruction_regions": localization["failing_regions"],
            "additional_basis_dimension_justified_by_this_audit": False,
            "next_action": next_action,
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "st006_comparison_performed": False,
            "st006_momentum_max_reference": 0.1082289305112118,
            "st006_momentum_volume_l2_reference": 0.10758432876230622,
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
    loc = report["regional_capacity_localization"]
    print("all_regions_pass=", loc["all_regions_pass"])
    print("failing_regions=", loc["failing_regions"])
    for name in REGIONS:
        region = loc["regions"][name]
        print(
            name,
            "rank=", region["normalized_column_rank"],
            "condition=", region["normalized_condition"],
            "max_abs_cosine=", region["max_abs_pairwise_cosine"],
            "top_two=", region["top_two_channels_by_raw_norm"],
        )


if __name__ == "__main__":
    main()
