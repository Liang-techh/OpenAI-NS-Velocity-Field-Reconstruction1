"""Target-free capacity screen for one interior swirl-redistribution basis.

The already-screened compact swirl mode changes the amount/concentration of
interior swirl with one sign.  This audit asks whether one additional *fixed*
sign-changing radial shape can instead move swirl response between the core and
interior shoulder without touching the physical support collars.

With s=r/Rp and q=z/Zp, define

    B = s (1-s^2)^5_+ (1-q^2)^5_+
    B_red = B (1 - (13/2) s^2).

The factor 13/2 is not fitted to a residual or image.  It is the analytic value
that makes B_red cylindrically L2-orthogonal to B in radius under the response
energy measure r dr:

    int B^2 (1-c s^2) r dr = 0  =>  c = Beta(2,11)/Beta(3,11)=13/2.

This is representation/visualization-capacity evidence only.  It does not
select a coefficient, evaluate a perturbed-candidate PDE residual, or promote a
new canonical velocity field.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_f20_capacity import F10, F20, _embedded_field
from .constrained_eq45_bipolar_interior_compact_swirl_capacity import (
    COMPACT_MODE,
    _compact_basis_velocity,
    _compact_velocity,
    _plateau_half_widths,
    _probe_groups,
    _response,
    _response_summary,
    _rings,
    _theta_unit,
)

TASK_ID = "CR003-BIPOLAR-INTERIOR-SWIRL-REDISTRIBUTION-CAPACITY-038"
REDISTRIBUTION_MODE = "INTERIOR_C4_SWIRL_REDISTRIBUTION"
ORTHOGONALIZER = 13.0 / 2.0
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "redistribution_coefficient_selected": False,
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


def _redistribution_basis_velocity(field, points) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    base = _compact_basis_velocity(field, points)
    radial_plateau, _ = _plateau_half_widths(field)
    r = np.hypot(points[:, 0], points[:, 1])
    factor = 1.0 - ORTHOGONALIZER * (r / radial_plateau) ** 2
    return base * factor[:, None]


def _redistribution_velocity(field, points, time: float, coefficient: float) -> np.ndarray:
    coefficient = float(coefficient)
    limit = float(field.parent.profile_basis.coefficient_limit)
    if not np.isfinite(coefficient) or abs(coefficient) > limit:
        raise ValueError("redistribution coefficient exceeds inherited coefficient_limit")
    return field.at_points(points, float(time)) + coefficient * _redistribution_basis_velocity(field, points)


def _mode_response(field, mode, step: float, points: np.ndarray, time: float) -> np.ndarray:
    if mode == REDISTRIBUTION_MODE:
        plus = _redistribution_velocity(field, points, time, step)
        minus = _redistribution_velocity(field, points, time, -step)
        return (plus - minus) / (2.0 * step)
    return _response(field, mode, step, points, time)


def _response_matrix(field, groups, *, time: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    points = np.concatenate(tuple(groups.values()))
    modes = (F10, F20, COMPACT_MODE, REDISTRIBUTION_MODE)
    columns = [_mode_response(field, mode, step, points, time).reshape(-1) for mode in modes]
    return points, np.column_stack(columns)


def _signed_swirl_mean(points: np.ndarray, velocity: np.ndarray) -> float:
    theta = _theta_unit(points)
    return float(np.mean(np.sum(np.asarray(velocity) * theta, axis=1)))


def _sign_localization(field) -> dict[str, Any]:
    radial_plateau, _ = _plateau_half_widths(field)
    inner = _rings((0.25, 0.45), (-0.35, 0.0, 0.35))
    shoulder = _rings((0.80, 1.00, 1.20), (-0.35, 0.0, 0.35))
    inner_response = _redistribution_basis_velocity(field, inner)
    shoulder_response = _redistribution_basis_velocity(field, shoulder)
    sign_change_fraction = float(np.sqrt(1.0 / ORTHOGONALIZER))
    return {
        "analytic_sign_change_r_over_Rp": sign_change_fraction,
        "analytic_sign_change_radius": float(radial_plateau * sign_change_fraction),
        "inner_signed_swirl_mean": _signed_swirl_mean(inner, inner_response),
        "shoulder_signed_swirl_mean": _signed_swirl_mean(shoulder, shoulder_response),
        "inner_total_rms": float(np.sqrt(np.mean(np.sum(inner_response * inner_response, axis=1)))),
        "shoulder_total_rms": float(np.sqrt(np.mean(np.sum(shoulder_response * shoulder_response, axis=1)))),
    }


def _visual_fingerprint(field, coefficient: float, *, time: float) -> dict[str, Any]:
    radial_plateau, axial_plateau = _plateau_half_widths(field)
    r = np.linspace(0.02, 1.98, 96)
    z = np.linspace(-1.98, 1.98, 129)
    rr, zz = np.meshgrid(r, z, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    radius = rr.ravel()
    weights = radius
    sign_radius = radial_plateau / np.sqrt(ORTHOGONALIZER)

    def metrics(velocity: np.ndarray) -> dict[str, float]:
        swirl = velocity[:, 1]
        energy = weights * swirl * swirl
        total = float(np.sum(energy))
        if not total > 0.0:
            raise RuntimeError("visual fingerprint encountered zero swirl energy")
        collar = (radius >= radial_plateau) | (np.abs(zz.ravel()) >= axial_plateau)
        inner = radius < sign_radius
        return {
            "swirl_energy_weight": total,
            "weighted_radial_centroid": float(np.sum(radius * energy) / total),
            "weighted_radial_rms": float(np.sqrt(np.sum(radius * radius * energy) / total)),
            "inner_swirl_energy_fraction": float(np.sum(energy[inner]) / total),
            "physical_plateau_collar_swirl_energy": float(np.sum(energy[collar])),
        }

    base_velocity = field.at_points(points, float(time))
    plus_velocity = _redistribution_velocity(field, points, time, coefficient)
    minus_velocity = _redistribution_velocity(field, points, time, -coefficient)
    base = metrics(base_velocity)
    plus = metrics(plus_velocity)
    minus = metrics(minus_velocity)
    collar_denom = max(abs(base["physical_plateau_collar_swirl_energy"]), 1e-300)
    max_collar_change = max(
        abs(plus["physical_plateau_collar_swirl_energy"] - base["physical_plateau_collar_swirl_energy"]),
        abs(minus["physical_plateau_collar_swirl_energy"] - base["physical_plateau_collar_swirl_energy"]),
    ) / collar_denom
    return {
        "time": float(time),
        "symmetric_diagnostic_coefficient": float(coefficient),
        "coefficient_is_candidate_bound": False,
        "grid_shape_r_z": [int(r.size), int(z.size)],
        "base": base,
        "plus": plus,
        "minus": minus,
        "max_relative_physical_plateau_collar_energy_change": float(max_collar_change),
        "plus_minus_radial_centroid_span": float(abs(plus["weighted_radial_centroid"] - minus["weighted_radial_centroid"])),
        "plus_minus_inner_energy_fraction_span": float(abs(plus["inner_swirl_energy_fraction"] - minus["inner_swirl_energy_fraction"])),
        "interpretation": (
            "Target-free morphology leverage only: a nonzero span means the fixed mode can redistribute "
            "interior swirl radially without changing velocity in the physical taper collars."
        ),
    }


def audit_interior_swirl_redistribution_capacity(
    *,
    time: float = 0.5,
    coefficient_steps: Iterable[float] = (0.02, 0.01),
    visual_trial_coefficient: float = 0.10,
) -> dict[str, Any]:
    steps = tuple(float(value) for value in coefficient_steps)
    if len(steps) != 2 or not steps[0] > steps[1] > 0.0:
        raise ValueError("coefficient_steps must be two decreasing positive values")
    field = _embedded_field()
    groups = _probe_groups()
    radial_plateau, axial_plateau = _plateau_half_widths(field)

    points, coarse = _response_matrix(field, groups, time=float(time), step=steps[0])
    _, fine = _response_matrix(field, groups, time=float(time), step=steps[1])
    refinement = float(np.linalg.norm(coarse - fine) / max(np.linalg.norm(fine), 1e-300))

    column_norms = np.linalg.norm(fine, axis=0)
    if np.any(column_norms <= 0.0):
        raise RuntimeError("a screened response column is numerically zero")
    normalized = fine / column_norms
    singular_values = np.linalg.svd(normalized, compute_uv=False)
    tol = max(normalized.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tol))
    condition = float(singular_values[0] / singular_values[-1]) if rank == 4 else float("inf")

    existing = fine[:, :3]
    redistribution = fine[:, 3]
    projection_coefficients, *_ = np.linalg.lstsq(existing, redistribution, rcond=None)
    projection = existing @ projection_coefficients
    novelty = float(
        np.linalg.norm(redistribution - projection) / max(np.linalg.norm(redistribution), 1e-300)
    )

    modes = ("F10", "F20", COMPACT_MODE, REDISTRIBUTION_MODE)
    overall = {
        name: _response_summary(points, fine[:, j].reshape(len(points), 3))
        for j, name in enumerate(modes)
    }
    regional: dict[str, Any] = {}
    cursor = 0
    for region_name, region_points in groups.items():
        n = len(region_points)
        block = fine[3 * cursor : 3 * (cursor + n), :]
        regional[region_name] = {
            name: _response_summary(region_points, block[:, j].reshape(n, 3))
            for j, name in enumerate(modes)
        }
        cursor += n

    red_energy_by_region = {
        name: float(regional[name][REDISTRIBUTION_MODE]["total_rms"] ** 2 * len(groups[name]))
        for name in groups
    }
    red_energy_total = sum(red_energy_by_region.values())
    collar_energy = red_energy_by_region["radial_flank"] + red_energy_by_region["axial_flank"]

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "time": float(time),
        "basis_definition": {
            "mode": REDISTRIBUTION_MODE,
            "classification": "autonomous_design",
            "radial_plateau_half_width": radial_plateau,
            "axial_plateau_half_width": axial_plateau,
            "shape": "B*(1-(13/2)*(r/Rp)^2), B=(r/Rp)*(1-(r/Rp)^2)^5_+*(1-(z/Zp)^2)^5_+ in u_theta",
            "orthogonalizer": ORTHOGONALIZER,
            "orthogonality_measure": "cylindrical response L2 radial measure r dr",
            "analytic_parent_inner_product": 0.0,
            "residual_fitted_shape_parameters": 0,
            "public_image_fitted_shape_parameters": 0,
        },
        "coefficient_steps": list(steps),
        "response_column_norms": [float(v) for v in column_norms],
        "normalized_response_singular_values": [float(v) for v in singular_values],
        "normalized_response_rank": rank,
        "normalized_response_condition_number": condition,
        "finite_difference_refinement_relative_change": refinement,
        "redistribution_novelty_outside_F10_F20_compact_span": novelty,
        "projection_coefficients_onto_F10_F20_compact": [float(v) for v in projection_coefficients],
        "overall_response": overall,
        "regional_response": regional,
        "redistribution_collar_response_energy_fraction": float(collar_energy / max(red_energy_total, 1e-300)),
        "sign_localization": _sign_localization(field),
        "visualization_fingerprint_leverage": _visual_fingerprint(
            field, float(visual_trial_coefficient), time=float(time)
        ),
        "parameter_growth": {
            "fixed_spatial_basis_shapes_added_if_materialized": 1,
            "scalar_coefficients_added_if_materialized": 1,
            "shape_parameters_fitted": 0,
            "temporal_degrees_added": 0,
        },
        "routing_contract": (
            "Keep this fixed redistribution direction only if a later materialization/validation lane needs "
            "core-versus-shoulder swirl control. Do not infer PDE improvement or select its coefficient here."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/bipolar_interior_swirl_redistribution_capacity/report.json",
    )
    args = parser.parse_args()
    report = audit_interior_swirl_redistribution_capacity()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
