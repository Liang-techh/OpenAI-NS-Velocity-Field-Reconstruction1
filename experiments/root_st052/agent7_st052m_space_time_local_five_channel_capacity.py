"""Space-time-local capacity audit for the existing five ST052-M morphology controls.

Preregistered in issue #1057 and stacked exactly on Agent-7 PR #1050 head
3b94159939ddc48811069a34bc84e5fe18e2fbb2.

#1033 audited the same five frozen tangent directions globally, #1042 localized
that capacity in broad morphology-relevant spatial regions, and #1050 localized
it by registered time.  This increment intersects those two localizations: at a
registered time, does each already-existing control retain leverage in the
specific core/radial/midplane/tip region where a future candidate-bound visual
discrepancy could occur?

No candidate, coefficient or basis is changed.  Scientific rank/conditioning
failure is preserved in the receipt rather than hidden as a CI failure.  The
regions and times are target-free routing coordinates, not numerical targets
inferred from the public OpenAI visualization.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

import agent7_st052m_time_local_five_channel_capacity as time_local

TASK_ID = "CR003-ST052M-SPACE-TIME-LOCAL-FIVE-CHANNEL-CAPACITY-126"
PREREG_ISSUE = 1057
SOURCE_PARENT_PR = 1050
SOURCE_PARENT_HEAD = "3b94159939ddc48811069a34bc84e5fe18e2fbb2"

REGIONS = tuple(time_local.regional.REGIONS)
TIME_SLICES = tuple(time_local.TIME_SLICES)
ACTIVE_RELATIVE_NORM_MIN = time_local.ACTIVE_RELATIVE_NORM_MIN
CONDITION_MAX = time_local.CONDITION_MAX
COSINE_MAX_ABS = time_local.COSINE_MAX_ABS
LEAKAGE_MAX = time_local.LEAKAGE_MAX

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
    "space_time_local_visualization_sensitivity_audited": True,
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


def _cell_key(time: float, region: str) -> str:
    return f"{time_local._time_key(time)}/{region}"


def _assert_source_lock() -> None:
    time_local._assert_source_lock()
    if time_local.TASK_ID != "CR003-ST052M-TIME-LOCAL-FIVE-CHANNEL-CAPACITY-125":
        raise RuntimeError("#1050 task identity drifted")
    if time_local.PREREG_ISSUE != 1048:
        raise RuntimeError("#1050 preregistration drifted")
    if time_local.SOURCE_PARENT_PR != 1042:
        raise RuntimeError("#1050 parent identity drifted")
    if tuple(time_local.regional.REGIONS) != (
        "inner_core",
        "outer_radial",
        "midplane",
        "axial_tip",
    ):
        raise RuntimeError("#1042 region set drifted")
    if TIME_SLICES != (0.25, 0.375, 0.50, 0.625, 0.75):
        raise RuntimeError("#1050 time grid drifted")
    if ACTIVE_RELATIVE_NORM_MIN != 1.0e-12:
        raise RuntimeError("#1050 activity gate drifted")
    if CONDITION_MAX != 25.0 or COSINE_MAX_ABS != 0.995:
        raise RuntimeError("joint conditioning gates drifted")


def _normalized_subspace(
    active: tuple[str, ...],
    region_vectors: dict[str, np.ndarray],
) -> tuple[int, float | None, list[float], dict[str, float | None], float | None]:
    """Return rank/condition/SVD/cosines without dropping tiny active columns."""
    if not active:
        return 0, None, [], {}, None

    columns: list[np.ndarray] = []
    nonzero: dict[str, bool] = {}
    for name in active:
        vector = np.asarray(region_vectors[name], dtype=float).ravel()
        norm = float(np.linalg.norm(vector))
        is_nonzero = bool(np.isfinite(norm) and norm > np.finfo(float).tiny)
        nonzero[name] = is_nonzero
        columns.append(vector / norm if is_nonzero else np.zeros_like(vector))

    matrix = np.column_stack(columns)
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.linalg.matrix_rank(matrix, tol=1.0e-10))
    condition = (
        float(singular[0] / singular[-1])
        if singular.size and singular[-1] > np.finfo(float).tiny
        else None
    )
    gram = matrix.T @ matrix
    cosines: dict[str, float | None] = {}
    finite: list[float] = []
    for i, j in itertools.combinations(range(len(active)), 2):
        key = f"{active[i]}__{active[j]}"
        if nonzero[active[i]] and nonzero[active[j]]:
            value = float(gram[i, j])
            cosines[key] = value
            finite.append(abs(value))
        else:
            cosines[key] = None
    max_abs = float(max(finite)) if finite else 0.0
    return rank, condition, [float(v) for v in singular], cosines, max_abs


def _component_fingerprints(
    points: np.ndarray,
    blocks: dict[str, np.ndarray],
    mask: np.ndarray,
) -> dict[str, Any]:
    region_points = np.asarray(points, dtype=float)[mask]
    out: dict[str, Any] = {}
    labels = ("u_r", "u_theta", "u_z")
    for name in time_local.regional.joint.CHANNELS:
        fp = time_local._component_fingerprint(region_points, np.asarray(blocks[name])[mask])
        rms = fp["cylindrical_component_rms"]
        fp["dominant_component"] = labels[int(np.argmax([rms[label] for label in labels]))]
        out[name] = fp
    return out


def _audit_cell(
    time: float,
    region: str,
    points: np.ndarray,
    mask: np.ndarray,
    blocks: dict[str, np.ndarray],
    full_norms: dict[str, float],
    global_active: tuple[str, ...],
) -> dict[str, Any]:
    channels = tuple(time_local.regional.joint.CHANNELS)
    region_vectors = {name: np.asarray(blocks[name], dtype=float)[mask].ravel() for name in channels}
    raw_norms = {name: float(np.linalg.norm(region_vectors[name])) for name in channels}
    relative_norms = {name: float(raw_norms[name] / full_norms[name]) for name in channels}
    locally_negligible = tuple(
        name for name in global_active if relative_norms[name] <= ACTIVE_RELATIVE_NORM_MIN
    )

    rank, condition, singular, cosines, max_abs_cosine = _normalized_subspace(
        global_active, region_vectors
    )
    active_subspace_pass = bool(
        not locally_negligible
        and rank == len(global_active)
        and condition is not None
        and condition <= CONDITION_MAX
        and max_abs_cosine is not None
        and max_abs_cosine < COSINE_MAX_ABS
    ) if global_active else bool(time_local._time_key(time) == "0.250")

    fingerprints = _component_fingerprints(points, blocks, mask)
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
    expected_endpoint_obstruction = bool(
        time_local._time_key(time) == "0.250" and not global_active
    )
    representation_obstruction = bool(expected_endpoint_obstruction or not active_subspace_pass)
    top_two = sorted(channels, key=lambda name: raw_norms[name], reverse=True)[:2]

    return {
        "time": float(time),
        "region": region,
        "point_count": int(np.count_nonzero(mask)),
        "point_fraction": float(np.mean(mask)),
        "global_active_channels": list(global_active),
        "global_active_channel_count": int(len(global_active)),
        "raw_regional_column_norms": raw_norms,
        "regional_norm_over_full_spacetime_channel_norm": relative_norms,
        "active_relative_norm_gate": ACTIVE_RELATIVE_NORM_MIN,
        "locally_negligible_globally_active_channels": list(locally_negligible),
        "active_subspace_rank": rank,
        "active_subspace_rank_target": int(len(global_active)),
        "active_subspace_condition": condition,
        "condition_gate": CONDITION_MAX,
        "active_subspace_singular_values": singular,
        "pairwise_cosines_active_subspace": cosines,
        "max_abs_pairwise_cosine_active_subspace": max_abs_cosine,
        "cosine_abs_gate": COSINE_MAX_ABS,
        "active_subspace_pass": active_subspace_pass,
        "component_fingerprints": fingerprints,
        "top_two_channels_by_raw_regional_response": list(top_two),
        "toroidal_rz_leakage_rms_max": toroidal_rz_leakage,
        "poloidal_theta_leakage_rms_max": poloidal_theta_leakage,
        "selectivity_guard_pass": selectivity_pass,
        "expected_start_endpoint_obstruction": expected_endpoint_obstruction,
        "local_representation_obstruction": representation_obstruction,
        "passes_preregistered_cell_gate": bool(active_subspace_pass and selectivity_pass),
    }


def space_time_local_capacity() -> dict[str, Any]:
    """Audit the existing five controls in all 20 frozen region×time cells."""
    _assert_source_lock()
    joint = time_local.regional.joint
    points, full_blocks, metadata = joint._build_tangent_blocks()
    masks, thresholds = time_local.regional._region_masks(points)
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

    by_time: dict[str, Any] = {}
    obstruction_cells: list[str] = []
    unexpected_failure_cells: list[str] = []
    expected_endpoint_cells: list[str] = []
    cell_count = 0

    for time in TIME_SLICES:
        key = time_local._time_key(time)
        blocks = time_local._build_blocks_at_time(
            float(time), points, parent_fn, s_alphas, s_children, k_children, toroidal_spatial
        )
        slice_norms = {
            name: float(np.linalg.norm(np.asarray(blocks[name], dtype=float).ravel()))
            for name in joint.CHANNELS
        }
        slice_relative = {name: float(slice_norms[name] / full_norms[name]) for name in joint.CHANNELS}
        global_active = tuple(
            name for name in joint.CHANNELS if slice_relative[name] > ACTIVE_RELATIVE_NORM_MIN
        )
        global_inactive = tuple(name for name in joint.CHANNELS if name not in global_active)
        expected_inactive = tuple(time_local.EXPECTED_INACTIVE[key])
        inactive_matches = bool(global_inactive == expected_inactive)

        regions: dict[str, Any] = {}
        for region in REGIONS:
            row = _audit_cell(
                float(time), region, points, masks[region], blocks, full_norms, global_active
            )
            regions[region] = row
            cell_count += 1
            cell_key = _cell_key(float(time), region)
            if row["local_representation_obstruction"]:
                obstruction_cells.append(cell_key)
            if row["expected_start_endpoint_obstruction"]:
                expected_endpoint_cells.append(cell_key)
            elif not row["passes_preregistered_cell_gate"]:
                unexpected_failure_cells.append(cell_key)

        by_time[key] = {
            "time": float(time),
            "slice_column_norms": slice_norms,
            "slice_norm_over_full_spacetime_channel_norm": slice_relative,
            "global_active_channels": list(global_active),
            "global_inactive_channels": list(global_inactive),
            "expected_global_inactive_channels": list(expected_inactive),
            "global_inactive_set_matches_preregistered_expectation": inactive_matches,
            "regions": regions,
        }

    unexpected_inactive_times = [
        key
        for key, row in by_time.items()
        if not row["global_inactive_set_matches_preregistered_expectation"]
    ]
    locally_negligible_by_cell = {
        _cell_key(float(key), region): row["locally_negligible_globally_active_channels"]
        for key, time_row in by_time.items()
        for region, row in time_row["regions"].items()
        if row["locally_negligible_globally_active_channels"]
    }

    return {
        **metadata,
        "region_order": list(REGIONS),
        "time_order": [float(v) for v in TIME_SLICES],
        "cell_count": cell_count,
        "median_region_thresholds": thresholds,
        "full_spacetime_channel_norms": full_norms,
        "active_relative_norm_gate": ACTIVE_RELATIVE_NORM_MIN,
        "by_time": by_time,
        "expected_start_endpoint_obstruction_cells": expected_endpoint_cells,
        "local_representation_obstruction_cells": obstruction_cells,
        "unexpected_capacity_failure_cells": unexpected_failure_cells,
        "locally_negligible_channels_by_cell": locally_negligible_by_cell,
        "unexpected_global_inactive_time_slices": unexpected_inactive_times,
        "all_global_time_activity_semantics_match": bool(not unexpected_inactive_times),
        "all_nonendpoint_cells_pass_preregistered_gates": bool(not unexpected_failure_cells),
    }


def build_report() -> dict[str, Any]:
    capacity = space_time_local_capacity()
    if capacity["unexpected_capacity_failure_cells"]:
        next_action = (
            "inspect_or_reparameterize_the_recorded_region_time_existing_channels_before_any_basis_growth"
        )
    else:
        next_action = (
            "wait_for_candidate_bound_morphology_discrepancy_and_route_to_the_smallest_existing_"
            "active_control_in_the_matching_region_time_cell"
        )
    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "source_parent": {"pr": SOURCE_PARENT_PR, "head": SOURCE_PARENT_HEAD},
        "frozen_protocol": {
            "channels": list(time_local.regional.joint.CHANNELS),
            "regions": list(REGIONS),
            "times": [float(v) for v in TIME_SLICES],
            "cell_count": len(REGIONS) * len(TIME_SLICES),
            "region_rule": "reuse #1042 median-derived masks on the exact frozen response cloud",
            "time_activity_rule": "reuse #1050 slice_norm/full_spacetime_norm > 1e-12",
            "condition_gate": CONDITION_MAX,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "selectivity_leakage_gate": LEAKAGE_MAX,
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "space_time_local_capacity": capacity,
        "decision": {
            "expected_start_endpoint_obstruction_cells": capacity[
                "expected_start_endpoint_obstruction_cells"
            ],
            "unexpected_local_representation_obstruction_cells": capacity[
                "unexpected_capacity_failure_cells"
            ],
            "additional_basis_dimension_justified_by_this_audit": False,
            "if_future_discrepancy_at_t025": (
                "inspect_parent_family_or_endpoint_nonpreserving_time_coordinate_before_basis_growth"
            ),
            "if_future_discrepancy_in_nonendpoint_cell": (
                "use_the_smallest_existing_well_conditioned_active_control_before_new_basis"
            ),
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
    cap = report["space_time_local_capacity"]
    print("all_nonendpoint_cells_pass=", cap["all_nonendpoint_cells_pass_preregistered_gates"])
    print("unexpected_capacity_failure_cells=", cap["unexpected_capacity_failure_cells"])
    print("locally_negligible_channels_by_cell=", cap["locally_negligible_channels_by_cell"])
    for key in (time_local._time_key(v) for v in TIME_SLICES):
        for region in REGIONS:
            row = cap["by_time"][key]["regions"][region]
            print(
                f"{key}/{region}",
                "active=", row["global_active_channels"],
                "local_negligible=", row["locally_negligible_globally_active_channels"],
                "rank=", row["active_subspace_rank"],
                "condition=", row["active_subspace_condition"],
                "max_abs_cosine=", row["max_abs_pairwise_cosine_active_subspace"],
            )


if __name__ == "__main__":
    main()
