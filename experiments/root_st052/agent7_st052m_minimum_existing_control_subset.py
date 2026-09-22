"""Find the smallest existing ST052-M control subset that preserves defect-space capacity.

Preregistered in issue #1130.  This is a target-free representation-capacity audit
stacked on Agent-7 #1123.  It consumes the already-frozen #1078 8x5
morphology-defect Jacobian and does not add a basis, select coefficients, or
change the candidate velocity.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

import agent7_st052m_defect_coordinate_controllability as parent

TASK_ID = "CR003-ST052M-MINIMUM-EXISTING-CONTROL-SUBSET-134"
PREREG_ISSUE = 1130
STACK_BASE_PR = 1123
STACK_BASE_HEAD = "216481c9c86b0982e3ab8c5b3d54a9b3641e8c01"
SCIENTIFIC_PARENT_PR = 1078
SCIENTIFIC_PARENT_TASK = "CR003-ST052M-DEFECT-COORDINATE-CONTROLLABILITY-128"

CHANNELS = (
    "amplitude",
    "radial_shape",
    "axial_turnover",
    "temporal_curvature",
    "toroidal_swirl",
)
CONDITION_MAX = 25.0
COSINE_MAX_ABS = 0.995
ROW_RESPONSE_FLOOR = 1.0e-12
RANK_TOL = 1.0e-10

TRUTH = {
    "candidate_velocity_changed": False,
    "canonical_velocity_changed": False,
    "saved_velocity_changed": False,
    "basis_dimension_changed": False,
    "new_spatial_basis_added": False,
    "new_temporal_basis_added": False,
    "coefficient_selected": False,
    "optimization_performed": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "st006_comparison_performed": False,
    "public_image_numeric_target_used": False,
    "pixel_similarity_objective_used": False,
    "renderer_or_camera_fit_used": False,
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
    if parent.TASK_ID != SCIENTIFIC_PARENT_TASK:
        raise RuntimeError("#1078 task identity drifted")
    if parent.PREREG_ISSUE != 1077:
        raise RuntimeError("#1078 preregistration identity drifted")
    if tuple(parent.DEFECT_COORDINATES) != (
        "inward_radial_index_t050",
        "swirl_to_poloidal_ratio_t050",
        "axial_vorticity_aspect_t050",
        "core_width_t050",
        "tip_radial_thickness_t050",
        "angular_rotation_radial_variation_t050",
        "core_width_change_025_to_075",
        "core_speed_change_025_to_075",
    ):
        raise RuntimeError("#1078 defect-coordinate order drifted")
    if tuple(parent.parent._joint().CHANNELS) != CHANNELS:
        raise RuntimeError("#1078 five-channel order drifted")
    if parent.CONDITION_MAX != CONDITION_MAX:
        raise RuntimeError("#1078 condition gate drifted")
    if parent.COSINE_MAX_ABS != COSINE_MAX_ABS:
        raise RuntimeError("#1078 cosine gate drifted")
    if parent.ROW_RESPONSE_FLOOR != ROW_RESPONSE_FLOOR:
        raise RuntimeError("#1078 row-response floor drifted")


def _matrix(audit: dict[str, Any], subset: tuple[str, ...], key: str) -> np.ndarray:
    rows = []
    for coordinate in parent.DEFECT_COORDINATES:
        mapping = audit["coordinates"][coordinate][key]
        rows.append([float(mapping[channel]) for channel in subset])
    matrix = np.asarray(rows, dtype=float)
    if matrix.shape != (len(parent.DEFECT_COORDINATES), len(subset)):
        raise RuntimeError("unexpected subset matrix shape")
    if not np.all(np.isfinite(matrix)):
        raise RuntimeError("nonfinite subset matrix")
    return matrix


def _pairwise_cosines(normalized: np.ndarray, subset: tuple[str, ...]) -> dict[str, float | None]:
    norms = np.linalg.norm(normalized, axis=0)
    out: dict[str, float | None] = {}
    for i, j in itertools.combinations(range(len(subset)), 2):
        key = f"{subset[i]}__{subset[j]}"
        if min(float(norms[i]), float(norms[j])) <= np.finfo(float).tiny:
            out[key] = None
        else:
            out[key] = float(np.dot(normalized[:, i], normalized[:, j]) / (norms[i] * norms[j]))
    return out


def subset_metrics(audit: dict[str, Any], subset: Iterable[str]) -> dict[str, Any]:
    subset = tuple(subset)
    if not subset:
        raise ValueError("subset must be non-empty")
    if len(set(subset)) != len(subset) or any(channel not in CHANNELS for channel in subset):
        raise ValueError("invalid channel subset")

    normalized = _matrix(audit, subset, "normalized_column_alignment_by_channel")
    raw = _matrix(audit, subset, "fine_raw_derivative_by_channel")
    singular = np.linalg.svd(normalized, compute_uv=False)
    rank = int(np.linalg.matrix_rank(normalized, tol=RANK_TOL))
    condition = None if singular[-1] <= np.finfo(float).tiny else float(singular[0] / singular[-1])
    cosines = _pairwise_cosines(normalized, subset)
    finite_cosines = [abs(value) for value in cosines.values() if value is not None]
    max_abs_cosine = 0.0 if len(subset) == 1 else (None if not finite_cosines else float(max(finite_cosines)))

    row_l2 = np.linalg.norm(raw, axis=1)
    responsive = {
        coordinate: bool(float(row_l2[index]) > ROW_RESPONSE_FLOOR)
        for index, coordinate in enumerate(parent.DEFECT_COORDINATES)
    }
    stability = {
        channel: bool(audit["derivative_stability"][channel]["passes"])
        for channel in subset
    }

    rank_pass = rank == len(subset)
    condition_pass = condition is not None and condition <= CONDITION_MAX
    cosine_pass = max_abs_cosine is not None and max_abs_cosine < COSINE_MAX_ABS
    response_pass = all(responsive.values())
    stability_pass = all(stability.values())

    failed = []
    if not stability_pass:
        failed.append("derivative_stability")
    if not rank_pass:
        failed.append("full_subset_rank")
    if not condition_pass:
        failed.append("condition")
    if not cosine_pass:
        failed.append("pairwise_cosine")
    if not response_pass:
        failed.append("coordinate_response")

    return {
        "channels": list(subset),
        "size": len(subset),
        "rank": rank,
        "rank_target": len(subset),
        "singular_values": [float(value) for value in singular],
        "condition": condition,
        "condition_gate": CONDITION_MAX,
        "pairwise_cosines": cosines,
        "max_abs_pairwise_cosine": max_abs_cosine,
        "cosine_abs_gate": COSINE_MAX_ABS,
        "channel_derivative_stability": stability,
        "all_included_channels_stable": stability_pass,
        "raw_response_l2_by_coordinate": {
            coordinate: float(row_l2[index])
            for index, coordinate in enumerate(parent.DEFECT_COORDINATES)
        },
        "coordinate_has_response": responsive,
        "all_coordinates_have_response": response_pass,
        "row_response_floor": ROW_RESPONSE_FLOOR,
        "passes_all_frozen_gates": bool(
            stability_pass and rank_pass and condition_pass and cosine_pass and response_pass
        ),
        "failed_gates": failed,
    }


def analyze_subsets(audit: dict[str, Any]) -> dict[str, Any]:
    by_size: dict[str, list[dict[str, Any]]] = {}
    minimum_size: int | None = None
    minimum_subsets: list[dict[str, Any]] = []

    for size in range(1, len(CHANNELS) + 1):
        rows = [subset_metrics(audit, subset) for subset in itertools.combinations(CHANNELS, size)]
        by_size[str(size)] = rows
        passing = [row for row in rows if row["passes_all_frozen_gates"]]
        if passing and minimum_size is None:
            minimum_size = size
            minimum_subsets = passing

    leave_one_out = {
        omitted: subset_metrics(audit, tuple(channel for channel in CHANNELS if channel != omitted))
        for omitted in CHANNELS
    }

    return {
        "channel_order": list(CHANNELS),
        "defect_coordinate_order": list(parent.DEFECT_COORDINATES),
        "all_nonempty_subsets_evaluated": 2 ** len(CHANNELS) - 1,
        "subsets_by_size": by_size,
        "minimum_admissible_subset_size": minimum_size,
        "minimum_admissible_subsets": [row["channels"] for row in minimum_subsets],
        "minimum_admissible_subset_metrics": minimum_subsets,
        "leave_one_out_four_channel_audit": leave_one_out,
        "existing_five_channel_family_locally_overcomplete_for_these_coordinates": bool(
            minimum_size is not None and minimum_size < len(CHANNELS)
        ),
        "all_five_channels_required_under_frozen_gates": bool(minimum_size == len(CHANNELS)),
        "no_admissible_existing_subset": minimum_size is None,
    }


def build_report() -> dict[str, Any]:
    _assert_source_lock()
    parent_audit = parent.controllability_audit()
    subset_audit = analyze_subsets(parent_audit)
    minimum = subset_audit["minimum_admissible_subset_size"]
    if minimum is None:
        interpretation = (
            "existing controls have a stability/conditioning/coverage obstruction under the frozen gates; "
            "do not infer that a new basis is required before a concrete discrepancy and obstruction diagnosis"
        )
    elif minimum < len(CHANNELS):
        interpretation = (
            "the five-control family is locally overcomplete for these autonomous defect coordinates; "
            "this is evidence against blind basis growth, not authorization to delete a production control"
        )
    else:
        interpretation = (
            "all five existing controls are locally required under the frozen gates; a sixth basis still needs "
            "a concrete candidate-bound discrepancy that the existing controls cannot repair"
        )

    return {
        "task_id": TASK_ID,
        "prereg_issue": PREREG_ISSUE,
        "stack_base": {"pr": STACK_BASE_PR, "head": STACK_BASE_HEAD},
        "scientific_parent": {"pr": SCIENTIFIC_PARENT_PR, "task": SCIENTIFIC_PARENT_TASK},
        "frozen_protocol": {
            "channels": list(CHANNELS),
            "defect_coordinates": list(parent.DEFECT_COORDINATES),
            "subset_count": 31,
            "condition_gate": CONDITION_MAX,
            "cosine_abs_gate": COSINE_MAX_ABS,
            "row_response_floor": ROW_RESPONSE_FLOOR,
            "source_jacobian": "#1078 fine defect-coordinate derivatives and normalized columns",
            "new_basis_dimension": 0,
            "coefficient_selected": None,
            "public_openai_numeric_target": None,
        },
        "minimum_existing_control_subset_audit": subset_audit,
        "decision": {
            "interpretation": interpretation,
            "new_basis_dimension_justified_by_this_audit": False,
            "candidate_mutation_authorized_by_this_audit": False,
            "actual_velocity_changed": False,
            "direct_visualization_fingerprint_improvement": 0.0,
            "closer_visualization_delivery_established": False,
        },
        "held_out_pde_residual": {
            "evaluated": False,
            "st006_comparison_performed": False,
            "reason": "subset audit only; no candidate velocity was changed",
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
    audit = report["minimum_existing_control_subset_audit"]
    print(json.dumps({
        "minimum_admissible_subset_size": audit["minimum_admissible_subset_size"],
        "minimum_admissible_subsets": audit["minimum_admissible_subsets"],
        "leave_one_out_failed_gates": {
            omitted: row["failed_gates"]
            for omitted, row in audit["leave_one_out_four_channel_audit"].items()
        },
        "decision": report["decision"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
