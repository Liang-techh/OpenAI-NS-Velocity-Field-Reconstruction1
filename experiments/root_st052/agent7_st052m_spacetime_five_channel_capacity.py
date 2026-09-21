"""Space-time localization of existing ST052-M morphology capacity.

Preregistered in issue #1057 and stacked exactly on Agent-7 PR #1050 head
3b94159939ddc48811069a34bc84e5fe18e2fbb2.

This increment adds no basis and changes no candidate. It intersects #1042's
frozen spatial masks with #1050's frozen time slices so a future candidate-bound
morphology discrepancy can be routed to the smallest existing control with
actual leverage in the affected region and time.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_time_local_five_channel_capacity as timeaudit

TASK_ID = "CR003-ST052M-SPACE-TIME-LOCAL-FIVE-CHANNEL-CAPACITY-126"
PREREG_ISSUE = 1057
SOURCE_PARENT_PR = 1050
SOURCE_PARENT_HEAD = "3b94159939ddc48811069a34bc84e5fe18e2fbb2"

REGIONS = timeaudit.regional.REGIONS
TIME_SLICES = timeaudit.TIME_SLICES
ACTIVE_RELATIVE_NORM_MIN = timeaudit.ACTIVE_RELATIVE_NORM_MIN
CONDITION_MAX = timeaudit.CONDITION_MAX
COSINE_MAX_ABS = timeaudit.COSINE_MAX_ABS
LEAKAGE_MAX = timeaudit.LEAKAGE_MAX

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
    "space_time_visualization_sensitivity_localized": True,
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
    timeaudit._assert_source_lock()
    if timeaudit.TASK_ID != "CR003-ST052M-TIME-LOCAL-FIVE-CHANNEL-CAPACITY-125":
        raise RuntimeError("#1050 task identity drifted")
    if timeaudit.PREREG_ISSUE != 1048:
        raise RuntimeError("#1050 preregistration drifted")
    if timeaudit.SOURCE_PARENT_PR != 1042:
        raise RuntimeError("#1050 parent identity drifted")
    if tuple(REGIONS) != ("inner_core", "outer_radial", "midplane", "axial_tip"):
        raise RuntimeError("regional mask set drifted")
    if tuple(TIME_SLICES) != (0.25, 0.375, 0.50, 0.625, 0.75):
        raise RuntimeError("time-slice set drifted")
    if ACTIVE_RELATIVE_NORM_MIN != 1.0e-12:
        raise RuntimeError("activity gate drifted")
    if CONDITION_MAX != 25.0 or COSINE_MAX_ABS != 0.995 or LEAKAGE_MAX != 1.0e-12:
        raise RuntimeError("capacity/selectivity gates drifted")


def _norm(block: np.ndarray) -> float:
    value = float(np.linalg.norm(np.asarray(block, dtype=float).ravel()))
    if not np.isfinite(value):
        raise RuntimeError("nonfinite response norm")
    return value


def _component_fingerprint(points: np.ndarray, block: np.ndarray) -> dict[str, Any]:
    return timeaudit._component_fingerprint(np.asarray(points, dtype=float), np.asarray(block, dtype=float))


def _audit_cell(
    time: float,
    region: str,
    points: np.ndarray,
    mask: np.ndarray,
    blocks: dict[str, np.ndarray],
    full_norms: dict[str, float],
) -> dict[str, Any]:
    channels = tuple(timeaudit.regional.joint.CHANNELS)
    global_norms = {name: _norm(blocks[name]) for name in channels}
    global_relative = {name: float(global_norms[name] / full_norms[name]) for name in channels}
    global_active = tuple(name for name in channels if global_relative[name] > ACTIVE_RELATIVE_NORM_MIN)
    expected_inactive = tuple(timeaudit.EXPECTED_INACTIVE[timeaudit._time_key(time)])
    expected_global_active = tuple(name for name in channels if name not in expected_inactive)
    global_semantics_match = bool(global_active == expected_global_active)

    region_points = np.asarray(points, dtype=float)[mask]
    local_blocks = {name: np.asarray(blocks[name], dtype=float)[mask] for name in channels}
    local_norms = {name: _norm(local_blocks[name]) for name in channels}
    local_relative = {name: float(local_norms[name] / full_norms[name]) for name in channels}
    local_active = tuple(
        name for name in global_active if local_relative[name] > ACTIVE_RELATIVE_NORM_MIN
    )
    locally_negligible = tuple(name for name in global_active if name not in local_active)

    if local_active:
        columns = []
        for name in local_active:
            norm = local_norms[name]
            columns.append(local_blocks[name].ravel() / norm)
        matrix = np.column_stack(columns)
        singular = np.linalg.svd(matrix, compute_uv=False)
        rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
        condition = float(singular[0] / singular[-1]) if singular[-1] > 0.0 else None
        gram = matrix.T @ matrix
        cosines: dict[str, float] = {}
        for i, j in itertools.combinations(range(len(local_active)), 2):
            cosines[f"{local_active[i]}__{local_active[j]}"] = float(gram[i, j])
        max_abs_cosine = float(max(abs(v) for v in cosines.values())) if cosines else 0.0
        local_subspace_pass = bool(
            rank == len(local_active)
            and condition is not None
            and condition <= CONDITION_MAX
            and max_abs_cosine < COSINE_MAX_ABS
        )
        singular_values = [float(v) for v in singular]
    else:
        rank = 0
        condition = None
        cosines = {}
        max_abs_cosine = None
        singular_values = []
        local_subspace_pass = bool(timeaudit._time_key(time) == "0.250")

    fingerprints = {
        name: _component_fingerprint(region_points, local_blocks[name]) for name in channels
    }
    tor = fingerprints["toroidal_swirl"]["cylindrical_component_rms"]
    toroidal_rz_leakage = float(max(abs(tor["u_r"]), abs(tor["u_z"])))
    poloidal_theta_leakage = float(
        max(
            abs(fingerprints[name]["cylindrical_component_rms"]["u_theta"])
            for name in channels
            if name != "toroidal_swirl"
        )
    )
    selectivity_pass = bool(
        toroidal_rz_leakage <= LEAKAGE_MAX and poloidal_theta_leakage <= LEAKAGE_MAX
    )
    top_two = sorted(global_active, key=lambda name: local_norms[name], reverse=True)[:2]
    endpoint_empty_expected = bool(timeaudit._time_key(time) == "0.250" and not global_active)
    local_obstruction = bool(locally_negligible or (global_active and not local_subspace_pass))
    passes = bool(
        global_semantics_match
        and selectivity_pass
        and (
            endpoint_empty_expected
            or (not locally_negligible and local_subspace_pass)
        )
    )

    return {
        "time": float(time),
        "region": region,
        "point_count": int(np.count_nonzero(mask)),
        "point_fraction": float(np.mean(mask)),
        "global_slice_raw_norms": global_norms,
        "global_slice_relative_to_full_spacetime_norm": global_relative,
        "expected_global_active_channels": list(expected_global_active),
        "global_active_channels": list(global_active),
        "global_activity_matches_preregistered_time_semantics": global_semantics_match,
        "regional_raw_norms": local_norms,
        "regional_relative_to_full_spacetime_norm": local_relative,
        "regional_active_channels": list(local_active),
        "locally_negligible_globally_active_channels": list(locally_negligible),
        "active_relative_norm_gate": ACTIVE_RELATIVE_NORM_MIN,
        "active_subspace_rank": rank,
        "active_subspace_target_rank": int(len(local_active)),
        "active_subspace_condition": condition,
        "condition_gate": CONDITION_MAX,
        "pairwise_cosines_active_subspace": cosines,
        "max_abs_pairwise_cosine_active_subspace": max_abs_cosine,
        "cosine_abs_gate": COSINE_MAX_ABS,
        "singular_values_active_subspace": singular_values,
        "active_subspace_pass": local_subspace_pass,
        "component_fingerprints": fingerprints,
        "top_two_existing_channels_by_regional_raw_norm": list(top_two),
        "toroidal_rz_leakage_rms_max": toroidal_rz_leakage,
        "poloidal_theta_leakage_rms_max": poloidal_theta_leakage,
        "selectivity_guard_pass": selectivity_pass,
        "expected_empty_start_endpoint": endpoint_empty_expected,
        "local_representation_obstruction": local_obstruction,
        "passes_preregistered_cell_gate": passes,
    }


def spacetime_capacity_localization() -> dict[str, Any]:
    """Audit existing-control leverage in each frozen region×time cell."""
    _assert_source_lock()
    joint = timeaudit.regional.joint
    points, full_blocks, metadata = joint._build_tangent_blocks()
    masks, thresholds = timeaudit.regional._region_masks(points)
    full_norms = {
        name: float(
            np.linalg.norm(
                np.concatenate([np.asarray(v, dtype=float).ravel() for v in full_blocks[name]])
            )
        )
        for name in joint.CHANNELS
    }
    if any((not np.isfinite(v)) or v <= np.finfo(float).tiny for v in full_norms.values()):
        raise RuntimeError("zero/nonfinite full-spacetime channel norm")

    parent_fn, s_alphas, _, s_children, _, k_children = joint.axial.build_local_families()
    toroidal_spatial = joint.toroidal.toroidal_unit_velocity(points)

    cells: dict[str, Any] = {}
    for time in TIME_SLICES:
        blocks = timeaudit._build_blocks_at_time(
            float(time), points, parent_fn, s_alphas, s_children, k_children, toroidal_spatial
        )
        for region in REGIONS:
            key = f"t={timeaudit._time_key(time)}::{region}"
            cells[key] = _audit_cell(
                float(time), region, points, masks[region], blocks, full_norms
            )

    failing = [key for key, value in cells.items() if not value["passes_preregistered_cell_gate"]]
    obstructions = [key for key, value in cells.items() if value["local_representation_obstruction"]]
    unexpected_activity = [
        key
        for key, value in cells.items()
        if not value["global_activity_matches_preregistered_time_semantics"]
    ]
    selectivity_failures = [
        key for key, value in cells.items() if not value["selectivity_guard_pass"]
    ]
    return {
        **metadata,
        "region_order": list(REGIONS),
        "time_order": [float(v) for v in TIME_SLICES],
        "cell_count": int(len(cells)),
        "median_thresholds": thresholds,
        "full_spacetime_channel_norms": full_norms,
        "cells": cells,
        "all_cells_pass": bool(not failing),
        "failing_cells": failing,
        "local_representation_obstruction_cells": obstructions,
        "unexpected_global_activity_cells": unexpected_activity,
        "selectivity_guard_failures": selectivity_failures,
    }


def build_report() -> dict[str, Any]:
    loc = spacetime_capacity_localization()
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "frozen_protocol": {
            "channels": list(timeaudit.regional.joint.CHANNELS),
            "regions": list(REGIONS),
            "times": [float(v) for v in TIME_SLICES],
            "region_rule": "reuse #1042 median-derived masks on exact frozen response cloud",
            "active_relative_norm_gate": ACTIVE_RELATIVE_NORM_MIN,
            "condition_gate": CONDITION_MAX,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "selectivity_leakage_gate": LEAKAGE_MAX,
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "space_time_capacity_localization": loc,
        "decision": {
            "local_representation_obstruction_cells": loc["local_representation_obstruction_cells"],
            "additional_basis_dimension_justified_by_this_audit": False,
            "routing_rule": (
                "use_candidate_bound_discrepancy_to_choose_smallest_existing_control_with_"
                "local_leverage; if the relevant cell is obstructed, inspect targeted_"
                "reparameterization_or_one_specific_channel_before_any_basis_growth"
            ),
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
    loc = report["space_time_capacity_localization"]
    print("all_cells_pass=", loc["all_cells_pass"])
    print("failing_cells=", loc["failing_cells"])
    print("local_representation_obstruction_cells=", loc["local_representation_obstruction_cells"])
    for key, cell in loc["cells"].items():
        print(
            key,
            "active=", cell["regional_active_channels"],
            "rank=", cell["active_subspace_rank"],
            "condition=", cell["active_subspace_condition"],
            "top_two=", cell["top_two_existing_channels_by_regional_raw_norm"],
        )


if __name__ == "__main__":
    main()
