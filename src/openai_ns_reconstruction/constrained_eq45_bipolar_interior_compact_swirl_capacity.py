"""Target-free capacity screen for one interior-localized compact swirl basis.

The basis shape is fixed *before* any residual or public-image comparison from the
existing physical support/plateau geometry.  It is an additive axisymmetric
swirl-only response

    d u_theta = a * (r/Rp) * (1-(r/Rp)^2)^5_+ * (1-(z/Zp)^2)^5_+,

where ``Rp`` and ``Zp`` are the already-declared radial/axial identity-plateau
half widths of :class:`AxisymmetricPhysicalTaper`.  The factor ``r`` makes the
Cartesian field regular on the axis, and the C4 envelopes make the correction
identically zero before the support collars begin.

This module is representation/visualization-capacity evidence only.  It does
not fit ``a``, does not evaluate a perturbed-candidate PDE residual, and does not
promote a new canonical velocity field.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_f20_capacity import (
    F10,
    F20,
    _embedded_field,
    _with_swirl_delta,
)

TASK_ID = "CR003-BIPOLAR-INTERIOR-COMPACT-SWIRL-CAPACITY-037"
COMPACT_MODE = "INTERIOR_C4_SWIRL"
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "compact_coefficient_selected": False,
    "residual_map_used_to_fit_basis_shape": False,
    "perturbed_candidate_pde_residual_evaluated": False,
    "force_or_pressure_changed": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
}


def _rings(radii, z_values, *, n_theta: int = 12) -> np.ndarray:
    theta = 2.0 * np.pi * (np.arange(n_theta) + 0.5) / n_theta
    rows = []
    for radius in radii:
        for z in z_values:
            rows.append(
                np.column_stack(
                    (
                        float(radius) * np.cos(theta),
                        float(radius) * np.sin(theta),
                        np.full(n_theta, float(z), dtype=float),
                    )
                )
            )
    return np.concatenate(rows)


def _probe_groups() -> dict[str, np.ndarray]:
    return {
        "core": _rings((0.25, 0.45, 0.65), (-0.35, 0.0, 0.35)),
        "interior": _rings((0.80, 1.00, 1.20), (-0.35, 0.0, 0.35)),
        "radial_flank": _rings((1.65, 1.75, 1.85), (-0.35, 0.0, 0.35)),
        "axial_flank": _rings((0.35, 0.75, 1.15), (-1.85, -1.70, 1.70, 1.85)),
    }


def _theta_unit(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    r = np.hypot(points[:, 0], points[:, 1])
    if np.any(r <= 1e-14):
        raise ValueError("theta direction is undefined on the axis")
    return np.column_stack((-points[:, 1] / r, points[:, 0] / r, np.zeros(len(points))))


def _plateau_half_widths(field) -> tuple[float, float]:
    taper = field.taper
    radial = float(taper.radial_support * np.sqrt(taper.radial_plateau_q))
    axial = float(taper.axial_half_height * np.sqrt(taper.axial_plateau_q))
    if not radial > 0.0 or not axial > 0.0:
        raise ValueError("invalid physical taper plateau geometry")
    return radial, axial


def _compact_basis_velocity(field, points) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must be finite")

    radial_plateau, axial_plateau = _plateau_half_widths(field)
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    r = np.hypot(x, y)
    qr = (r / radial_plateau) ** 2
    qz = (z / axial_plateau) ** 2
    br = np.where(qr < 1.0, (1.0 - qr) ** 5, 0.0)
    bz = np.where(qz < 1.0, (1.0 - qz) ** 5, 0.0)
    u_theta = (r / radial_plateau) * br * bz

    cos_theta = np.divide(x, r, out=np.ones_like(r), where=r > 0.0)
    sin_theta = np.divide(y, r, out=np.zeros_like(r), where=r > 0.0)
    out = np.column_stack((-u_theta * sin_theta, u_theta * cos_theta, np.zeros_like(r)))
    out[r == 0.0] = 0.0
    if not np.all(np.isfinite(out)):
        raise RuntimeError("compact swirl basis became nonfinite")
    return out


def _compact_velocity(field, points, time: float, coefficient: float) -> np.ndarray:
    coefficient = float(coefficient)
    limit = float(field.parent.profile_basis.coefficient_limit)
    if not np.isfinite(coefficient) or abs(coefficient) > limit:
        raise ValueError("compact coefficient exceeds inherited coefficient_limit")
    return field.at_points(points, float(time)) + coefficient * _compact_basis_velocity(field, points)


def _response(field, mode, step: float, points: np.ndarray, time: float) -> np.ndarray:
    if mode == COMPACT_MODE:
        plus = _compact_velocity(field, points, time, step)
        minus = _compact_velocity(field, points, time, -step)
    else:
        plus_field = _with_swirl_delta(field, mode, step)
        minus_field = _with_swirl_delta(field, mode, -step)
        plus = plus_field.at_points(points, time)
        minus = minus_field.at_points(points, time)
    return (plus - minus) / (2.0 * step)


def _rms(values) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _response_summary(points: np.ndarray, response: np.ndarray) -> dict[str, float]:
    theta = _theta_unit(points)
    swirl = np.sum(response * theta, axis=1)
    total = np.sqrt(np.sum(response * response, axis=1))
    return {
        "swirl_rms": _rms(swirl),
        "total_rms": _rms(total),
        "swirl_fraction": float(np.linalg.norm(swirl) / max(np.linalg.norm(response), 1e-300)),
    }


def _spatial_response_matrix(field, groups, *, time: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    points = np.concatenate(tuple(groups.values()))
    columns = []
    for mode in (F10, F20, COMPACT_MODE):
        columns.append(_response(field, mode, step, points, time).reshape(-1))
    return points, np.column_stack(columns)


def _visual_swirl_fingerprint(field, coefficient: float, *, time: float = 0.5) -> dict[str, Any]:
    radial_plateau, axial_plateau = _plateau_half_widths(field)
    r = np.linspace(0.02, 1.98, 96)
    z = np.linspace(-1.98, 1.98, 129)
    rr, zz = np.meshgrid(r, z, indexing="ij")
    points = np.column_stack((rr.ravel(), np.zeros(rr.size), zz.ravel()))
    weights = rr.ravel()

    def metrics(velocity: np.ndarray) -> dict[str, float]:
        # At y=0, x>0 the Cartesian v component is the cylindrical swirl u_theta.
        swirl = velocity[:, 1]
        energy = weights * swirl * swirl
        total = float(np.sum(energy))
        if not total > 0.0:
            raise RuntimeError("visual fingerprint encountered zero swirl energy")
        radial_rms = float(np.sqrt(np.sum((rr.ravel() ** 2) * energy) / total))
        axial_rms = float(np.sqrt(np.sum((zz.ravel() ** 2) * energy) / total))
        collar = (rr.ravel() >= radial_plateau) | (np.abs(zz.ravel()) >= axial_plateau)
        collar_energy = float(np.sum(energy[collar]))
        return {
            "swirl_energy_weight": total,
            "weighted_radial_rms": radial_rms,
            "weighted_axial_rms": axial_rms,
            "physical_plateau_collar_swirl_energy": collar_energy,
            "physical_plateau_collar_swirl_energy_fraction": collar_energy / total,
        }

    base_velocity = field.at_points(points, time)
    plus_velocity = _compact_velocity(field, points, time, coefficient)
    minus_velocity = _compact_velocity(field, points, time, -coefficient)
    base = metrics(base_velocity)
    plus = metrics(plus_velocity)
    minus = metrics(minus_velocity)
    denom = max(abs(base["physical_plateau_collar_swirl_energy"]), 1e-300)
    collar_change = max(
        abs(plus["physical_plateau_collar_swirl_energy"] - base["physical_plateau_collar_swirl_energy"]),
        abs(minus["physical_plateau_collar_swirl_energy"] - base["physical_plateau_collar_swirl_energy"]),
    ) / denom
    radial_leverage = max(
        abs(plus["weighted_radial_rms"] - base["weighted_radial_rms"]),
        abs(minus["weighted_radial_rms"] - base["weighted_radial_rms"]),
    ) / max(abs(base["weighted_radial_rms"]), 1e-300)
    return {
        "time": float(time),
        "symmetric_trial_coefficient": float(coefficient),
        "grid_shape_r_z": [int(r.size), int(z.size)],
        "base": base,
        "plus": plus,
        "minus": minus,
        "max_relative_physical_plateau_collar_energy_change": float(collar_change),
        "max_fractional_weighted_radial_rms_leverage": float(radial_leverage),
        "interpretation": (
            "This is target-free morphology leverage, not an OpenAI-image score: the compact mode can "
            "change interior swirl concentration while leaving velocities in the physical taper collars unchanged."
        ),
    }


def audit_interior_compact_swirl_capacity(
    *,
    time: float = 0.5,
    coefficient_steps: Iterable[float] = (0.02, 0.01),
    visual_trial_coefficient: float = 0.25,
) -> dict[str, Any]:
    steps = tuple(float(value) for value in coefficient_steps)
    if len(steps) != 2 or not steps[0] > steps[1] > 0.0:
        raise ValueError("coefficient_steps must be two decreasing positive values")
    field = _embedded_field()
    groups = _probe_groups()
    radial_plateau, axial_plateau = _plateau_half_widths(field)

    points, coarse = _spatial_response_matrix(field, groups, time=float(time), step=steps[0])
    _, fine = _spatial_response_matrix(field, groups, time=float(time), step=steps[1])
    refinement = float(np.linalg.norm(coarse - fine) / max(np.linalg.norm(fine), 1e-300))

    column_norms = np.linalg.norm(fine, axis=0)
    if np.any(column_norms <= 0.0):
        raise RuntimeError("a screened response column is numerically zero")
    normalized = fine / column_norms
    singular_values = np.linalg.svd(normalized, compute_uv=False)
    tol = max(normalized.shape) * np.finfo(float).eps * singular_values[0]
    rank = int(np.sum(singular_values > tol))
    condition = float(singular_values[0] / singular_values[-1]) if rank == 3 else float("inf")

    existing = fine[:, :2]
    compact = fine[:, 2]
    projection_coefficients, *_ = np.linalg.lstsq(existing, compact, rcond=None)
    projection = existing @ projection_coefficients
    novelty = float(np.linalg.norm(compact - projection) / max(np.linalg.norm(compact), 1e-300))

    region_rows: dict[str, Any] = {}
    cursor = 0
    mode_names = ("F10", "F20", COMPACT_MODE)
    for region_name, region_points in groups.items():
        n = len(region_points)
        block = fine[3 * cursor : 3 * (cursor + n), :]
        region_rows[region_name] = {}
        for j, name in enumerate(mode_names):
            response = block[:, j].reshape(n, 3)
            region_rows[region_name][name] = _response_summary(region_points, response)
        cursor += n

    overall = {}
    for j, name in enumerate(mode_names):
        overall[name] = _response_summary(points, fine[:, j].reshape(len(points), 3))

    compact_energy_by_region = {
        name: float(region_rows[name][COMPACT_MODE]["total_rms"] ** 2 * len(groups[name]))
        for name in groups
    }
    compact_energy_total = sum(compact_energy_by_region.values())
    collar_energy = compact_energy_by_region["radial_flank"] + compact_energy_by_region["axial_flank"]

    visual = _visual_swirl_fingerprint(
        field,
        float(visual_trial_coefficient),
        time=float(time),
    )
    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "time": float(time),
        "basis_definition": {
            "mode": COMPACT_MODE,
            "classification": "autonomous_design",
            "radial_plateau_half_width": radial_plateau,
            "axial_plateau_half_width": axial_plateau,
            "shape": "(r/Rp)*(1-(r/Rp)^2)^5_+*(1-(z/Zp)^2)^5_+ in u_theta",
            "residual_fitted_shape_parameters": 0,
        },
        "coefficient_steps": list(steps),
        "response_column_norms": [float(v) for v in column_norms],
        "normalized_response_singular_values": [float(v) for v in singular_values],
        "normalized_response_rank": rank,
        "normalized_response_condition_number": condition,
        "finite_difference_refinement_relative_change": refinement,
        "compact_novelty_outside_F10_F20_span": novelty,
        "projection_coefficients_onto_F10_F20": [float(v) for v in projection_coefficients],
        "overall_response": overall,
        "regional_response": region_rows,
        "compact_collar_response_energy_fraction": float(collar_energy / max(compact_energy_total, 1e-300)),
        "visualization_fingerprint_leverage": visual,
        "parameter_growth": {
            "fixed_spatial_basis_shapes_added_if_materialized": 1,
            "scalar_coefficients_added_if_materialized": 1,
            "shape_parameters_fitted": 0,
            "temporal_degrees_added": 0,
        },
        "routing_contract": (
            "Keep this one-mode block only if its independent interior localization is useful to the next "
            "materialization/optimization lane. Do not infer PDE improvement or choose its coefficient from this screen."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/bipolar_interior_compact_swirl_capacity/report.json",
    )
    args = parser.parse_args()
    report = audit_interior_compact_swirl_capacity()
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
