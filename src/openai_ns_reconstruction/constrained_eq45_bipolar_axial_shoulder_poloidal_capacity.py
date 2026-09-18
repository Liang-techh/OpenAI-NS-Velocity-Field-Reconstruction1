"""Target-free axial-shoulder poloidal capacity screen for the capped bipolar field.

The preceding compact-poloidal screen found a useful bulk axial-thickness knob
but did not move q90/q99 axial vorticity reach.  This increment therefore tests
one fixed follow-up direction with no fitted shape parameters:

    A_theta = .5*r*q^3*(1-(r/Rp)^2)^6_+*(1-q^2)^6_+,
    q = z/Zp.

Relative to the integrated q-weighted compact mode, the extra q^2 suppresses
the midplane and shifts response toward axial shoulders while preserving the
same axis regularity, bipolar parity and C4 velocity support inside the existing
physical plateau.  No residual map or public image is used to choose the shape.

This module is representation/visualization-capacity evidence only.  It does not
select a coefficient, change the saved candidate, fit pressure/forcing, evaluate
a perturbed-candidate PDE residual, or establish OpenAI visual correspondence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_compact_poloidal_capacity import (
    COMPACT_POLOIDAL,
    PHI01,
    PHI03,
    _compact_poloidal_basis_velocity,
    _compact_velocity,
    _cylindrical_components,
    _novelty,
    _probe_groups,
    _rank_condition,
    _response,
    _source_field,
    _weighted_quantile,
    _with_phi_delta,
)
from .constrained_eq45_bipolar_interior_compact_swirl_capacity import _plateau_half_widths, _rings

TASK_ID = "CR003-BIPOLAR-AXIAL-SHOULDER-POLOIDAL-CAPACITY-042"
AXIAL_SHOULDER_POLOIDAL = "AXIAL_SHOULDER_C4_ODD_Z_POLOIDAL"
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_materialized": False,
    "axial_shoulder_poloidal_coefficient_selected": False,
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


def _axial_shoulder_basis_velocity(field, points) -> np.ndarray:
    """Evaluate curl(A_theta e_theta) for the fixed q^3 shoulder-weighted mode."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")
    if not np.all(np.isfinite(points)):
        raise ValueError("points must be finite")

    radial_plateau, axial_plateau = _plateau_half_widths(field)
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    radius = np.hypot(x, y)
    s = (radius / radial_plateau) ** 2
    q = z / axial_plateau

    radial_inside = s < 1.0
    axial_inside = np.abs(q) < 1.0
    active = radial_inside & axial_inside
    one_minus_s = np.where(radial_inside, 1.0 - s, 0.0)
    one_minus_q2 = np.where(axial_inside, 1.0 - q * q, 0.0)
    ar5 = one_minus_s**5
    ar6 = one_minus_s**6
    az5 = one_minus_q2**5
    az6 = one_minus_q2**6

    # b(q)=q^3(1-q^2)^6, b'(q)=3q^2(1-q^2)^5(1-5q^2).
    b = q**3 * az6
    bprime = 3.0 * q * q * az5 * (1.0 - 5.0 * q * q)
    u_r = -0.5 * radius * ar6 * bprime / axial_plateau
    u_z = ar5 * (1.0 - 7.0 * s) * b
    u_r = np.where(active, u_r, 0.0)
    u_z = np.where(active, u_z, 0.0)

    cos_theta = np.divide(x, radius, out=np.ones_like(radius), where=radius > 0.0)
    sin_theta = np.divide(y, radius, out=np.zeros_like(radius), where=radius > 0.0)
    out = np.column_stack((u_r * cos_theta, u_r * sin_theta, u_z))
    if not np.all(np.isfinite(out)):
        raise RuntimeError("axial-shoulder poloidal basis became nonfinite")
    return out


def _shoulder_velocity(field, points, time: float, coefficient: float) -> np.ndarray:
    coefficient = float(coefficient)
    if not np.isfinite(coefficient):
        raise ValueError("axial-shoulder coefficient must be finite")
    limit = float(field.parent.profile_basis.coefficient_limit)
    if abs(coefficient) > limit:
        raise ValueError("axial-shoulder coefficient exceeds inherited implementation guard")
    return field.at_points(np.asarray(points, dtype=float), float(time)) + coefficient * _axial_shoulder_basis_velocity(field, points)


def _shoulder_response(field, step: float, points: np.ndarray, time: float) -> np.ndarray:
    step = float(step)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    plus = _shoulder_velocity(field, points, time, step)
    minus = _shoulder_velocity(field, points, time, -step)
    return (plus - minus) / (2.0 * step)


def _response_matrix(field, groups, *, time: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    points = np.concatenate(tuple(groups.values()))
    columns = [
        _response(field, PHI01, step, points, time).reshape(-1),
        _response(field, PHI03, step, points, time).reshape(-1),
        _response(field, COMPACT_POLOIDAL, step, points, time).reshape(-1),
        _shoulder_response(field, step, points, time).reshape(-1),
    ]
    return points, np.column_stack(columns)


def _basis_structure_checks(field) -> dict[str, float]:
    positive = _rings((0.20, 0.45, 0.80, 1.20), (0.20, 0.55, 1.10), n_theta=8)
    mirrored = positive.copy()
    mirrored[:, 2] *= -1.0
    vp = _axial_shoulder_basis_velocity(field, positive)
    vm = _axial_shoulder_basis_velocity(field, mirrored)
    urp, utp, uzp = _cylindrical_components(positive, vp)
    urm, utm, uzm = _cylindrical_components(mirrored, vm)

    flank = np.concatenate((_probe_groups()["radial_flank"], _probe_groups()["axial_flank"]))
    midplane = _rings((0.25, 0.55, 0.85, 1.15), (0.0,), n_theta=8)

    check = np.array(
        [[0.31, 0.22, 0.17], [0.62, -0.27, -0.41], [0.88, 0.31, 0.67], [1.02, -0.36, -0.73]],
        dtype=float,
    )
    h = 1.0e-5
    divergence = np.zeros(len(check), dtype=float)
    for axis in range(3):
        offset = np.zeros_like(check)
        offset[:, axis] = h
        plus = _axial_shoulder_basis_velocity(field, check + offset)
        minus = _axial_shoulder_basis_velocity(field, check - offset)
        divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * h)

    return {
        "radial_even_max_abs_error": float(np.max(np.abs(urp - urm))),
        "axial_odd_max_abs_error": float(np.max(np.abs(uzp + uzm))),
        "swirl_response_rms": float(np.sqrt(np.mean(np.concatenate((utp, utm)) ** 2))),
        "support_flank_max_abs_velocity": float(np.max(np.abs(_axial_shoulder_basis_velocity(field, flank)))),
        "midplane_max_abs_velocity": float(np.max(np.abs(_axial_shoulder_basis_velocity(field, midplane)))),
        "cartesian_fd_divergence_max_abs": float(np.max(np.abs(divergence))),
    }


def _basis_axial_localization(field, *, grid_size: int = 31) -> dict[str, float]:
    if grid_size < 17 or grid_size % 2 == 0:
        raise ValueError("grid_size must be odd and >=17")
    radial_plateau, axial_plateau = _plateau_half_widths(field)
    axis = np.linspace(-radial_plateau, radial_plateau, grid_size)
    zaxis = np.linspace(-axial_plateau, axial_plateau, grid_size)
    x, y, z = np.meshgrid(axis, axis, zaxis, indexing="ij")
    inside = (x * x + y * y) < radial_plateau**2
    points = np.column_stack((x[inside], y[inside], z[inside]))

    def metrics(response: np.ndarray) -> tuple[float, float]:
        energy = np.sum(response * response, axis=1)
        total = float(np.sum(energy))
        if total <= np.finfo(float).tiny:
            raise RuntimeError("basis response is numerically inactive")
        abs_z = np.abs(points[:, 2])
        centroid = float(np.sum(abs_z * energy) / total / axial_plateau)
        outer = float(np.sum(energy[abs_z >= 0.50 * axial_plateau]) / total)
        return centroid, outer

    center = metrics(_compact_poloidal_basis_velocity(field, points))
    shoulder = metrics(_axial_shoulder_basis_velocity(field, points))
    return {
        "center_compact_abs_z_energy_centroid_over_Zp": center[0],
        "shoulder_abs_z_energy_centroid_over_Zp": shoulder[0],
        "centroid_ratio_shoulder_over_center": shoulder[0] / max(center[0], 1e-300),
        "center_compact_outer_half_response_energy_fraction": center[1],
        "shoulder_outer_half_response_energy_fraction": shoulder[1],
        "outer_half_fraction_ratio_shoulder_over_center": shoulder[1] / max(center[1], 1e-300),
    }


def _velocity_for_mode(field, points: np.ndarray, time: float, mode, coefficient: float) -> np.ndarray:
    if mode is None:
        return field.at_points(points, float(time))
    if mode == AXIAL_SHOULDER_POLOIDAL:
        return _shoulder_velocity(field, points, time, coefficient)
    if mode == COMPACT_POLOIDAL:
        return _compact_velocity(field, points, time, coefficient)
    if mode in (PHI01, PHI03):
        return _with_phi_delta(field, mode, coefficient).at_points(points, float(time))
    raise ValueError(f"unknown morphology mode {mode}")


def _vorticity_morphology(
    field,
    *,
    time: float,
    mode=None,
    coefficient: float = 0.0,
    grid_size: int = 33,
    box_half_width: float = 1.95,
) -> dict[str, float]:
    if grid_size < 17 or grid_size % 2 == 0:
        raise ValueError("grid_size must be odd and >=17")
    radial_plateau, axial_plateau = _plateau_half_widths(field)
    axis = np.linspace(-float(box_half_width), float(box_half_width), int(grid_size))
    spacing = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))
    velocity = _velocity_for_mode(field, points, time, mode, coefficient)
    velocity = np.asarray(velocity, dtype=float).reshape(grid_size, grid_size, grid_size, 3)

    u, v, w = velocity[..., 0], velocity[..., 1], velocity[..., 2]
    two_h = 2.0 * spacing
    dw_dy = (w[1:-1, 2:, 1:-1] - w[1:-1, :-2, 1:-1]) / two_h
    dv_dz = (v[1:-1, 1:-1, 2:] - v[1:-1, 1:-1, :-2]) / two_h
    du_dz = (u[1:-1, 1:-1, 2:] - u[1:-1, 1:-1, :-2]) / two_h
    dw_dx = (w[2:, 1:-1, 1:-1] - w[:-2, 1:-1, 1:-1]) / two_h
    dv_dx = (v[2:, 1:-1, 1:-1] - v[:-2, 1:-1, 1:-1]) / two_h
    du_dy = (u[1:-1, 2:, 1:-1] - u[1:-1, :-2, 1:-1]) / two_h
    omega_sq = (dw_dy - dv_dz) ** 2 + (du_dz - dw_dx) ** 2 + (dv_dx - du_dy) ** 2
    if not np.all(np.isfinite(omega_sq)):
        raise RuntimeError("vorticity reconstruction became nonfinite")
    total = float(np.sum(omega_sq))
    if total <= np.finfo(float).tiny:
        raise RuntimeError("vorticity is numerically inactive")

    xi = x[1:-1, 1:-1, 1:-1]
    yi = y[1:-1, 1:-1, 1:-1]
    zi = z[1:-1, 1:-1, 1:-1]
    radius_sq = xi * xi + yi * yi
    radial_rms = float(np.sqrt(np.sum(radius_sq * omega_sq) / total))
    axial_rms = float(np.sqrt(np.sum(zi * zi * omega_sq) / total))
    abs_z = np.abs(zi)
    axial_q90 = _weighted_quantile(abs_z, omega_sq, 0.90)
    axial_q99 = _weighted_quantile(abs_z, omega_sq, 0.99)
    outer_half = float(np.sum(omega_sq[abs_z >= 0.50 * axial_plateau]) / total)
    outer_065 = float(np.sum(omega_sq[abs_z >= 0.65 * axial_plateau]) / total)
    outer_075 = float(np.sum(omega_sq[abs_z >= 0.75 * axial_plateau]) / total)
    collar = (np.sqrt(radius_sq) >= radial_plateau) | (abs_z >= axial_plateau)
    collar_fraction = float(np.sum(omega_sq[collar]) / total)
    return {
        "radial_rms_over_Rp": radial_rms / radial_plateau,
        "axial_rms_over_Zp": axial_rms / axial_plateau,
        "aspect_z_over_r": axial_rms / max(radial_rms, 1e-300),
        "axial_q90_over_Zp": axial_q90 / axial_plateau,
        "axial_q99_over_Zp": axial_q99 / axial_plateau,
        "outer_half_enstrophy_fraction": outer_half,
        "outer_065_enstrophy_fraction": outer_065,
        "outer_075_enstrophy_fraction": outer_075,
        "physical_plateau_collar_enstrophy_fraction": collar_fraction,
    }


def _morphology_vector(metrics: dict[str, float]) -> np.ndarray:
    return np.array(
        [
            metrics["radial_rms_over_Rp"],
            metrics["axial_rms_over_Zp"],
            metrics["aspect_z_over_r"],
            metrics["outer_half_enstrophy_fraction"],
            metrics["outer_065_enstrophy_fraction"],
            metrics["outer_075_enstrophy_fraction"],
        ],
        dtype=float,
    )


def _morphology_response(field, mode, step: float, *, time: float, grid_size: int) -> np.ndarray:
    if mode == PHI03:
        base = _morphology_vector(_vorticity_morphology(field, time=time, grid_size=grid_size))
        plus = _morphology_vector(
            _vorticity_morphology(field, time=time, mode=mode, coefficient=step, grid_size=grid_size)
        )
        return (plus - base) / step
    plus = _morphology_vector(
        _vorticity_morphology(field, time=time, mode=mode, coefficient=step, grid_size=grid_size)
    )
    minus = _morphology_vector(
        _vorticity_morphology(field, time=time, mode=mode, coefficient=-step, grid_size=grid_size)
    )
    return (plus - minus) / (2.0 * step)


def audit_axial_shoulder_poloidal_capacity(
    *,
    time: float = 0.5,
    coefficient_steps: Iterable[float] = (0.02, 0.01),
    grid_size: int = 33,
    diagnostic_trial_coefficient: float = 0.25,
) -> dict[str, Any]:
    steps = tuple(float(value) for value in coefficient_steps)
    if len(steps) != 2 or not steps[0] > steps[1] > 0.0:
        raise ValueError("coefficient_steps must contain two decreasing positive values")
    trial = float(diagnostic_trial_coefficient)
    if not np.isfinite(trial) or not 0.0 < trial <= 1.0:
        raise ValueError("diagnostic_trial_coefficient must lie in (0,1]")

    field = _source_field()
    groups = _probe_groups()
    points, coarse = _response_matrix(field, groups, time=float(time), step=steps[0])
    _, fine = _response_matrix(field, groups, time=float(time), step=steps[1])
    response_refinement = float(np.linalg.norm(coarse - fine) / max(np.linalg.norm(fine), 1e-300))
    response_rank = _rank_condition(fine)
    shoulder_novelty = _novelty(fine[:, :3], fine[:, 3])
    per_vector_rms = np.linalg.norm(fine, axis=0) / np.sqrt(len(points))

    morphology_by_step = []
    modes = (PHI01, PHI03, COMPACT_POLOIDAL, AXIAL_SHOULDER_POLOIDAL)
    for step in steps:
        morphology_by_step.append(
            np.column_stack(
                [_morphology_response(field, mode, step, time=float(time), grid_size=grid_size) for mode in modes]
            )
        )
    morph_coarse, morph_fine = morphology_by_step
    morphology_refinement = float(
        np.linalg.norm(morph_coarse - morph_fine) / max(np.linalg.norm(morph_fine), 1e-300)
    )
    morphology_rank = _rank_condition(morph_fine)

    base_morph = _vorticity_morphology(field, time=float(time), grid_size=grid_size)
    plus_morph = _vorticity_morphology(
        field, time=float(time), mode=AXIAL_SHOULDER_POLOIDAL, coefficient=trial, grid_size=grid_size
    )
    minus_morph = _vorticity_morphology(
        field, time=float(time), mode=AXIAL_SHOULDER_POLOIDAL, coefficient=-trial, grid_size=grid_size
    )

    axial_efficiency = np.abs(morph_fine[1, :]) / np.maximum(per_vector_rms, 1e-300)
    outer065_efficiency = np.abs(morph_fine[4, :]) / np.maximum(per_vector_rms, 1e-300)
    outer075_efficiency = np.abs(morph_fine[5, :]) / np.maximum(per_vector_rms, 1e-300)
    radial_efficiency = np.abs(morph_fine[0, :]) / np.maximum(per_vector_rms, 1e-300)
    radial_plateau, axial_plateau = _plateau_half_widths(field)

    names = ("Phi01", "Phi03", COMPACT_POLOIDAL, AXIAL_SHOULDER_POLOIDAL)
    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "time": float(time),
        "basis_definition": {
            "mode": AXIAL_SHOULDER_POLOIDAL,
            "representation": "axisymmetric_vector_potential_A_theta_equivalently_streamfunction",
            "classification": "autonomous_design",
            "shape": "A_theta=.5*r*q^3*(1-(r/Rp)^2)^6_+*(1-(z/Zp)^2)^6_+",
            "radial_plateau_half_width": radial_plateau,
            "axial_plateau_half_width": axial_plateau,
            "shape_parameters_fitted": 0,
            "parity": "u_r even in z; u_z odd in z; u_theta=0",
            "design_reason": "extra q^2 suppresses the midplane and shifts capacity toward axial shoulders",
        },
        "coefficient_steps": list(steps),
        "basis_structure_checks": _basis_structure_checks(field),
        "basis_axial_localization": _basis_axial_localization(field),
        "public_velocity_response": {
            **response_rank,
            "finite_difference_refinement_relative_change": response_refinement,
            "shoulder_novelty_outside_Phi01_Phi03_center_compact_span": shoulder_novelty,
            "per_vector_response_rms": {name: float(value) for name, value in zip(names, per_vector_rms)},
        },
        "vorticity_morphology": {
            "grid_size": int(grid_size),
            "fingerprint_rows": [
                "radial_rms_over_Rp",
                "axial_rms_over_Zp",
                "aspect_z_over_r",
                "outer_half_enstrophy_fraction",
                "outer_065_enstrophy_fraction",
                "outer_075_enstrophy_fraction",
            ],
            "finest_response_columns": {
                name: [float(v) for v in morph_fine[:, index]] for index, name in enumerate(names)
            },
            "response_rank": morphology_rank,
            "finite_difference_refinement_relative_change": morphology_refinement,
            "axial_rms_leverage_per_public_velocity_rms": {
                name: float(value) for name, value in zip(names, axial_efficiency)
            },
            "outer_065_leverage_per_public_velocity_rms": {
                name: float(value) for name, value in zip(names, outer065_efficiency)
            },
            "outer_075_leverage_per_public_velocity_rms": {
                name: float(value) for name, value in zip(names, outer075_efficiency)
            },
            "radial_rms_leverage_per_public_velocity_rms": {
                name: float(value) for name, value in zip(names, radial_efficiency)
            },
            "baseline": base_morph,
            "diagnostic_plus": plus_morph,
            "diagnostic_minus": minus_morph,
            "diagnostic_trial_coefficient": trial,
            "diagnostic_axial_rms_span_over_Zp": float(abs(plus_morph["axial_rms_over_Zp"] - minus_morph["axial_rms_over_Zp"])),
            "diagnostic_axial_q90_span_over_Zp": float(abs(plus_morph["axial_q90_over_Zp"] - minus_morph["axial_q90_over_Zp"])),
            "diagnostic_axial_q99_span_over_Zp": float(abs(plus_morph["axial_q99_over_Zp"] - minus_morph["axial_q99_over_Zp"])),
            "diagnostic_outer_065_span": float(abs(plus_morph["outer_065_enstrophy_fraction"] - minus_morph["outer_065_enstrophy_fraction"])),
            "diagnostic_outer_075_span": float(abs(plus_morph["outer_075_enstrophy_fraction"] - minus_morph["outer_075_enstrophy_fraction"])),
            "diagnostic_radial_rms_span_over_Rp": float(abs(plus_morph["radial_rms_over_Rp"] - minus_morph["radial_rms_over_Rp"])),
        },
        "parameter_growth_if_materialized": {
            "fixed_spatial_basis_shapes_added": 1,
            "scalar_coefficients_added": 1,
            "shape_parameters_fitted": 0,
            "temporal_degrees_added": 0,
        },
        "routing_contract": (
            "Promote this one shoulder-weighted compact poloidal direction only if its target-free axial-reach "
            "leverage is materially stronger than the integrated center-weighted compact mode and later governed "
            "PDE screening accepts an explicit bounded coefficient. If q90/q99 and outer-tip fingerprints remain "
            "effectively immobile, stop polynomial q-weight growth and route the remaining tip defect to a support-"
            "localized poloidal/support-geometry experiment rather than adding more swirl or high-degree monomials."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/bipolar_axial_shoulder_poloidal_capacity/report.json")
    parser.add_argument("--grid-size", type=int, default=33)
    args = parser.parse_args()
    report = audit_axial_shoulder_poloidal_capacity(grid_size=args.grid_size)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
