"""Joint target-free capacity audit for the two existing interior swirl redistributions.

This increment does not add another basis mode.  It asks whether the already
screened fixed radial and axial redistributions span two genuinely independent
visualization-facing directions before any mixed r-z swirl mode is considered.

With the inherited compact interior swirl response

    B = (r/Rp) (1-(r/Rp)^2)^5_+ (1-(z/Zp)^2)^5_+,

the two fixed directions are

    B_r = B * (1-(13/2)(r/Rp)^2),
    B_z = B * (1-23(z/Zp)^2).

The constants are analytic rather than fitted.  Under the separable cylindrical
response-energy measure B^2 r dr dz, each direction is orthogonal to B and the
two redistributions are mutually orthogonal.  The numerical audit below uses the
public [u,v,w] evaluator and target-free swirl-energy morphology fingerprints to
check whether this structural independence survives in the delivered candidate.

No coefficient is selected, no public image or residual map is fitted, and no
perturbed-candidate PDE residual is evaluated.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_f20_capacity import F10, F20, _embedded_field
from .constrained_eq45_bipolar_interior_compact_swirl_capacity import (
    COMPACT_MODE,
    _compact_basis_velocity,
    _plateau_half_widths,
    _response,
    _rings,
)
from .constrained_eq45_bipolar_interior_axial_swirl_redistribution_capacity import (
    AXIAL_REDISTRIBUTION_MODE,
    ORTHOGONALIZER as AXIAL_ORTHOGONALIZER,
    _axial_redistribution_basis_velocity,
)

TASK_ID = "CR003-BIPOLAR-INTERIOR-SWIRL-2D-JACOBIAN-040"
RADIAL_REDISTRIBUTION_MODE = "INTERIOR_C4_SWIRL_REDISTRIBUTION"
RADIAL_ORTHOGONALIZER = 13.0 / 2.0
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_mode_added": False,
    "radial_redistribution_coefficient_selected": False,
    "axial_redistribution_coefficient_selected": False,
    "residual_map_used_to_fit_basis_shape": False,
    "public_image_used_to_fit_basis_shape": False,
    "perturbed_candidate_pde_residual_evaluated": False,
    "force_or_pressure_changed": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _radial_redistribution_basis_velocity(field, points) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    base = _compact_basis_velocity(field, points)
    radial_plateau, _ = _plateau_half_widths(field)
    radius = np.hypot(points[:, 0], points[:, 1])
    factor = 1.0 - RADIAL_ORTHOGONALIZER * (radius / radial_plateau) ** 2
    return base * factor[:, None]


def _joint_velocity(
    field,
    points: np.ndarray,
    time: float,
    *,
    radial_coefficient: float = 0.0,
    axial_coefficient: float = 0.0,
) -> np.ndarray:
    radial_coefficient = float(radial_coefficient)
    axial_coefficient = float(axial_coefficient)
    if not np.isfinite(radial_coefficient) or not np.isfinite(axial_coefficient):
        raise ValueError("redistribution coefficients must be finite")
    limit = float(field.parent.profile_basis.coefficient_limit)
    if abs(radial_coefficient) > limit or abs(axial_coefficient) > limit:
        raise ValueError("diagnostic redistribution coefficient exceeds inherited implementation guard")
    points = np.asarray(points, dtype=float)
    base = field.at_points(points, float(time))
    return (
        base
        + radial_coefficient * _radial_redistribution_basis_velocity(field, points)
        + axial_coefficient * _axial_redistribution_basis_velocity(field, points)
    )


def _probe_groups() -> dict[str, np.ndarray]:
    return {
        "center": _rings((0.25, 0.55, 0.85), (-0.15, 0.0, 0.15)),
        "radial_shoulder": _rings((0.80, 1.00, 1.20), (-0.35, 0.0, 0.35)),
        "axial_shoulder": _rings((0.25, 0.55, 0.85), (-0.75, -0.55, 0.55, 0.75)),
        "interior_corner": _rings((0.80, 1.05), (-0.80, -0.60, 0.60, 0.80)),
        "radial_flank": _rings((1.65, 1.75, 1.85), (-0.35, 0.0, 0.35)),
        "axial_flank": _rings((0.35, 0.75, 1.15), (-1.85, -1.70, 1.70, 1.85)),
    }


def _joint_response(field, points: np.ndarray, time: float, step: float, axis: str) -> np.ndarray:
    if axis == "radial":
        plus = _joint_velocity(field, points, time, radial_coefficient=step)
        minus = _joint_velocity(field, points, time, radial_coefficient=-step)
    elif axis == "axial":
        plus = _joint_velocity(field, points, time, axial_coefficient=step)
        minus = _joint_velocity(field, points, time, axial_coefficient=-step)
    else:
        raise ValueError("axis must be radial or axial")
    return (plus - minus) / (2.0 * step)


def _response_matrix(field, points: np.ndarray, *, time: float, step: float) -> np.ndarray:
    columns = [
        _response(field, F10, step, points, time).reshape(-1),
        _response(field, F20, step, points, time).reshape(-1),
        _response(field, COMPACT_MODE, step, points, time).reshape(-1),
        _joint_response(field, points, time, step, "radial").reshape(-1),
        _joint_response(field, points, time, step, "axial").reshape(-1),
    ]
    return np.column_stack(columns)


def _rank_condition(matrix: np.ndarray) -> dict[str, Any]:
    matrix = np.asarray(matrix, dtype=float)
    norms = np.linalg.norm(matrix, axis=0)
    if np.any(norms <= 0.0):
        raise RuntimeError("response matrix contains a zero column")
    normalized = matrix / norms
    singular = np.linalg.svd(normalized, compute_uv=False)
    tol = max(normalized.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.sum(singular > tol))
    condition = float(singular[0] / singular[-1]) if rank == normalized.shape[1] else float("inf")
    return {
        "column_norms": [float(value) for value in norms],
        "normalized_singular_values": [float(value) for value in singular],
        "normalized_rank": rank,
        "normalized_condition_number": condition,
    }


def _novelty(existing: np.ndarray, proposed: np.ndarray) -> float:
    coefficients, *_ = np.linalg.lstsq(existing, proposed, rcond=None)
    projection = existing @ coefficients
    return float(np.linalg.norm(proposed - projection) / max(np.linalg.norm(proposed), 1e-300))


def _cosine(left: np.ndarray, right: np.ndarray) -> float:
    denom = np.linalg.norm(left) * np.linalg.norm(right)
    if denom <= 0.0:
        raise RuntimeError("cannot compute cosine for zero response")
    return float(np.dot(left, right) / denom)


def _fingerprint_grid(field) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    radial_plateau, axial_plateau = _plateau_half_widths(field)
    radius_axis = np.linspace(0.02, 1.98, 96)
    axial_axis = np.linspace(-1.98, 1.98, 129)
    rr, zz = np.meshgrid(radius_axis, axial_axis, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    return points, rr.ravel(), zz.ravel(), float(radial_plateau), float(axial_plateau)


def _morphology_metrics(
    velocity: np.ndarray,
    radius: np.ndarray,
    axial: np.ndarray,
    radial_plateau: float,
    axial_plateau: float,
) -> dict[str, float]:
    velocity = np.asarray(velocity, dtype=float)
    swirl = velocity[:, 1]
    energy = radius * swirl * swirl
    total = float(np.sum(energy))
    if not total > 0.0:
        raise RuntimeError("zero swirl energy in morphology audit")
    abs_z = np.abs(axial)
    radial_inner = radius < radial_plateau / np.sqrt(RADIAL_ORTHOGONALIZER)
    axial_center = abs_z < axial_plateau / np.sqrt(AXIAL_ORTHOGONALIZER)
    collar = (radius >= radial_plateau) | (abs_z >= axial_plateau)
    return {
        "radial_centroid_over_Rp": float(np.sum(radius * energy) / total / radial_plateau),
        "abs_z_centroid_over_Zp": float(np.sum(abs_z * energy) / total / axial_plateau),
        "inner_radial_swirl_energy_fraction": float(np.sum(energy[radial_inner]) / total),
        "central_axial_swirl_energy_fraction": float(np.sum(energy[axial_center]) / total),
        "physical_plateau_collar_swirl_energy": float(np.sum(energy[collar])),
        "total_swirl_energy_weight": total,
    }


_MORPHOLOGY_KEYS = (
    "radial_centroid_over_Rp",
    "abs_z_centroid_over_Zp",
    "inner_radial_swirl_energy_fraction",
    "central_axial_swirl_energy_fraction",
)


def _metric_vector(metrics: dict[str, float]) -> np.ndarray:
    return np.asarray([metrics[key] for key in _MORPHOLOGY_KEYS], dtype=float)


def _morphology_jacobian(field, *, time: float, step: float) -> tuple[np.ndarray, dict[str, float]]:
    points, radius, axial, radial_plateau, axial_plateau = _fingerprint_grid(field)
    base = _morphology_metrics(
        field.at_points(points, float(time)), radius, axial, radial_plateau, axial_plateau
    )
    columns = []
    for axis in ("radial", "axial"):
        if axis == "radial":
            plus_velocity = _joint_velocity(field, points, time, radial_coefficient=step)
            minus_velocity = _joint_velocity(field, points, time, radial_coefficient=-step)
        else:
            plus_velocity = _joint_velocity(field, points, time, axial_coefficient=step)
            minus_velocity = _joint_velocity(field, points, time, axial_coefficient=-step)
        plus = _morphology_metrics(plus_velocity, radius, axial, radial_plateau, axial_plateau)
        minus = _morphology_metrics(minus_velocity, radius, axial, radial_plateau, axial_plateau)
        columns.append((_metric_vector(plus) - _metric_vector(minus)) / (2.0 * step))
    return np.column_stack(columns), base


def _joint_corner_fingerprint(field, *, time: float, coefficient: float) -> dict[str, Any]:
    points, radius, axial, radial_plateau, axial_plateau = _fingerprint_grid(field)
    base = _morphology_metrics(
        field.at_points(points, float(time)), radius, axial, radial_plateau, axial_plateau
    )
    rows: dict[str, dict[str, float]] = {}
    for radial_sign, axial_sign in itertools.product((-1.0, 1.0), repeat=2):
        label = f"r{int(radial_sign):+d}_z{int(axial_sign):+d}"
        velocity = _joint_velocity(
            field,
            points,
            time,
            radial_coefficient=radial_sign * coefficient,
            axial_coefficient=axial_sign * coefficient,
        )
        rows[label] = _morphology_metrics(velocity, radius, axial, radial_plateau, axial_plateau)
    radial_values = [row["radial_centroid_over_Rp"] for row in rows.values()]
    axial_values = [row["abs_z_centroid_over_Zp"] for row in rows.values()]
    collar_denom = max(abs(base["physical_plateau_collar_swirl_energy"]), 1e-300)
    max_collar_change = max(
        abs(row["physical_plateau_collar_swirl_energy"] - base["physical_plateau_collar_swirl_energy"])
        / collar_denom
        for row in rows.values()
    )
    return {
        "symmetric_diagnostic_coefficient": float(coefficient),
        "coefficient_is_candidate_bound": False,
        "base": base,
        "corners": rows,
        "radial_centroid_span": float(max(radial_values) - min(radial_values)),
        "abs_z_centroid_span": float(max(axial_values) - min(axial_values)),
        "max_relative_physical_plateau_collar_energy_change": float(max_collar_change),
    }


def audit_interior_swirl_2d_jacobian(
    *,
    time: float = 0.5,
    response_steps: Iterable[float] = (0.02, 0.01),
    morphology_steps: Iterable[float] = (0.02, 0.01),
    joint_trial_coefficient: float = 0.10,
) -> dict[str, Any]:
    response_steps = tuple(float(value) for value in response_steps)
    morphology_steps = tuple(float(value) for value in morphology_steps)
    if len(response_steps) != 2 or not response_steps[0] > response_steps[1] > 0.0:
        raise ValueError("response_steps must be two decreasing positive values")
    if len(morphology_steps) != 2 or not morphology_steps[0] > morphology_steps[1] > 0.0:
        raise ValueError("morphology_steps must be two decreasing positive values")
    if not np.isfinite(joint_trial_coefficient) or joint_trial_coefficient <= 0.0:
        raise ValueError("joint_trial_coefficient must be finite and positive")

    field = _embedded_field()
    groups = _probe_groups()
    points = np.concatenate(tuple(groups.values()))

    response_coarse = _response_matrix(field, points, time=float(time), step=response_steps[0])
    response_fine = _response_matrix(field, points, time=float(time), step=response_steps[1])
    response_refinement = float(
        np.linalg.norm(response_coarse - response_fine) / max(np.linalg.norm(response_fine), 1e-300)
    )
    response_summary = _rank_condition(response_fine)

    radial = response_fine[:, 3]
    axial = response_fine[:, 4]
    radial_novelty = _novelty(response_fine[:, :3], radial)
    axial_novelty_after_radial = _novelty(response_fine[:, :4], axial)
    pair_summary = _rank_condition(response_fine[:, 3:5])
    pair_cosine = _cosine(radial, axial)

    morphology_coarse, base_metrics = _morphology_jacobian(
        field, time=float(time), step=morphology_steps[0]
    )
    morphology_fine, _ = _morphology_jacobian(field, time=float(time), step=morphology_steps[1])
    morphology_refinement = float(
        np.linalg.norm(morphology_coarse - morphology_fine)
        / max(np.linalg.norm(morphology_fine), 1e-300)
    )
    morphology_summary = _rank_condition(morphology_fine)
    centroid_jacobian = morphology_fine[:2, :]
    centroid_singular = np.linalg.svd(centroid_jacobian, compute_uv=False)
    centroid_rank = int(np.linalg.matrix_rank(centroid_jacobian))
    centroid_condition = (
        float(centroid_singular[0] / centroid_singular[-1])
        if centroid_rank == 2 and centroid_singular[-1] > 0.0
        else float("inf")
    )
    radial_primary_fraction = float(
        abs(morphology_fine[0, 0])
        / max(abs(morphology_fine[0, 0]) + abs(morphology_fine[1, 0]), 1e-300)
    )
    axial_primary_fraction = float(
        abs(morphology_fine[1, 1])
        / max(abs(morphology_fine[0, 1]) + abs(morphology_fine[1, 1]), 1e-300)
    )

    stop_on_capacity_ground = bool(
        pair_summary["normalized_rank"] == 2
        and morphology_summary["normalized_rank"] == 2
        and np.isfinite(pair_summary["normalized_condition_number"])
        and np.isfinite(morphology_summary["normalized_condition_number"])
    )

    radial_plateau, axial_plateau = _plateau_half_widths(field)
    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "time": float(time),
        "basis_contract": {
            "parent": COMPACT_MODE,
            "radial_mode": RADIAL_REDISTRIBUTION_MODE,
            "axial_mode": AXIAL_REDISTRIBUTION_MODE,
            "radial_shape": "B*(1-(13/2)*(r/Rp)^2)",
            "axial_shape": "B*(1-23*(z/Zp)^2)",
            "radial_plateau_half_width": float(radial_plateau),
            "axial_plateau_half_width": float(axial_plateau),
            "radial_orthogonalizer": RADIAL_ORTHOGONALIZER,
            "axial_orthogonalizer": AXIAL_ORTHOGONALIZER,
            "analytic_parent_radial_inner_product": 0.0,
            "analytic_parent_axial_inner_product": 0.0,
            "analytic_radial_axial_inner_product": 0.0,
            "analytic_measure": "separable cylindrical response-energy measure B^2 r dr dz",
            "residual_fitted_shape_parameters": 0,
            "public_image_fitted_shape_parameters": 0,
        },
        "public_velocity_response": {
            **response_summary,
            "response_steps": list(response_steps),
            "finite_difference_refinement_relative_change": response_refinement,
            "radial_novelty_outside_F10_F20_compact_span": radial_novelty,
            "axial_novelty_outside_F10_F20_compact_radial_span": axial_novelty_after_radial,
            "radial_axial_pair": {
                **pair_summary,
                "cosine": pair_cosine,
            },
        },
        "morphology_jacobian": {
            "metric_order": list(_MORPHOLOGY_KEYS),
            "control_order": [RADIAL_REDISTRIBUTION_MODE, AXIAL_REDISTRIBUTION_MODE],
            "base_metrics": base_metrics,
            "coarse_step": morphology_steps[0],
            "fine_step": morphology_steps[1],
            "fine_jacobian": [[float(value) for value in row] for row in morphology_fine],
            "finite_difference_refinement_relative_change": morphology_refinement,
            "normalized_column_singular_values": morphology_summary["normalized_singular_values"],
            "normalized_rank": morphology_summary["normalized_rank"],
            "normalized_condition_number": morphology_summary["normalized_condition_number"],
            "column_norms": morphology_summary["column_norms"],
            "centroid_jacobian": [[float(value) for value in row] for row in centroid_jacobian],
            "centroid_singular_values": [float(value) for value in centroid_singular],
            "centroid_rank": centroid_rank,
            "centroid_condition_number": centroid_condition,
            "radial_control_primary_centroid_fraction": radial_primary_fraction,
            "axial_control_primary_centroid_fraction": axial_primary_fraction,
        },
        "joint_visualization_fingerprint": _joint_corner_fingerprint(
            field, time=float(time), coefficient=float(joint_trial_coefficient)
        ),
        "parameter_growth_this_increment": {
            "new_spatial_basis_shapes_added": 0,
            "new_scalar_coefficients_selected": 0,
            "shape_parameters_fitted": 0,
            "temporal_degrees_added": 0,
        },
        "routing": {
            "rank_obstruction_to_two_mode_interior_swirl_control_observed": not stop_on_capacity_ground,
            "mixed_rz_swirl_mode_justified_by_capacity_rank_alone": False if stop_on_capacity_ground else None,
            "stop_interior_swirl_basis_growth_on_capacity_ground": stop_on_capacity_ground,
            "interpretation": (
                "If both existing controls remain rank-2 in public velocity and target-free morphology, "
                "there is no capacity-rank reason to add a mixed r-z swirl mode. Use the existing radial "
                "control only for an independently established core-vs-shoulder mismatch and the axial "
                "control only for an independently established center-vs-cap mismatch. Gross axial-tip "
                "or vorticity-core reach should continue to the poloidal/support lane rather than more swirl basis growth."
            ),
        },
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/bipolar_interior_swirl_2d_jacobian/report.json",
    )
    args = parser.parse_args()
    report = audit_interior_swirl_2d_jacobian()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
