"""Target-free compact poloidal capacity screen for the capped bipolar field.

The preceding Agent-7 joint swirl screen found no local velocity-rank reason to
keep growing the swirl family, while its axial swirl redistribution had weak
gross axial morphology leverage.  This increment therefore tests exactly one
minimal *poloidal* alternative, without fitting a residual map or public image.

The proposed correction is generated as ``curl(A_theta e_theta)`` with

    A_theta = 0.5 * r * q * (1-s)^6_+ * (1-q^2)^6_+,
    s = (r/Rp)^2, q = z/Zp,

where ``Rp`` and ``Zp`` are the already-declared physical-taper plateau half
widths.  It is axis-regular, divergence-free by construction, odd in z in the
axial component, even in z in the radial component, and identically zero before
the physical support collars begin.  Positive coefficient gives inward radial
motion on the midplane and opposite-sign axial motion above/below it, matching
the parity class of the current bipolar seed.

This module is expression/visualization-capacity evidence only.  It does not
select a coefficient, change the saved candidate, fit pressure/forcing, evaluate
a perturbed-candidate PDE residual, or establish OpenAI visual correspondence.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .constrained_eq45_bipolar_f20_capacity import _source_field
from .constrained_eq45_bipolar_interior_compact_swirl_capacity import (
    _plateau_half_widths,
    _rings,
)

TASK_ID = "CR003-BIPOLAR-COMPACT-POLOIDAL-AXIAL-CAPACITY-041"
PHI01 = (0, 1)
PHI03 = (0, 3)
COMPACT_POLOIDAL = "COMPACT_C4_ODD_Z_POLOIDAL"
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_materialized": False,
    "compact_poloidal_coefficient_selected": False,
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


def _with_phi_delta(field, mode: tuple[int, int], delta: float):
    basis = field.parent.profile_basis
    if mode not in basis.mode_indices:
        raise ValueError(f"missing Phi mode {mode}")
    index = basis.mode_indices.index(mode)
    coefficients = list(basis.phi_coefficients)
    coefficients[index] += float(delta)
    if not np.isfinite(coefficients[index]) or abs(coefficients[index]) > basis.coefficient_limit:
        raise ValueError("Phi perturbation exceeds inherited coefficient_limit")
    return replace(
        field,
        parent=replace(
            field.parent,
            profile_basis=replace(basis, phi_coefficients=tuple(coefficients)),
        ),
    )


def _compact_poloidal_basis_velocity(field, points) -> np.ndarray:
    """Evaluate the fixed curl(A_theta e_theta) response in Cartesian form."""
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

    # For A_theta = .5*r*Ar(s)*B(q), B=q*(1-q^2)^6:
    # u_r = -d_z A_theta and u_z = (1/r)d_r(r A_theta).
    b = q * az6
    bprime = az5 * (1.0 - 13.0 * q * q)
    u_r = -0.5 * radius * ar6 * bprime / axial_plateau
    u_z = ar5 * (1.0 - 7.0 * s) * b
    u_r = np.where(active, u_r, 0.0)
    u_z = np.where(active, u_z, 0.0)

    cos_theta = np.divide(x, radius, out=np.ones_like(radius), where=radius > 0.0)
    sin_theta = np.divide(y, radius, out=np.zeros_like(radius), where=radius > 0.0)
    out = np.column_stack((u_r * cos_theta, u_r * sin_theta, u_z))
    if not np.all(np.isfinite(out)):
        raise RuntimeError("compact poloidal basis became nonfinite")
    return out


def _compact_velocity(field, points, time: float, coefficient: float) -> np.ndarray:
    coefficient = float(coefficient)
    if not np.isfinite(coefficient):
        raise ValueError("compact poloidal coefficient must be finite")
    # This is only an implementation guard for the diagnostic, not a selected bound.
    limit = float(field.parent.profile_basis.coefficient_limit)
    if abs(coefficient) > limit:
        raise ValueError("compact poloidal coefficient exceeds inherited implementation guard")
    return field.at_points(points, float(time)) + coefficient * _compact_poloidal_basis_velocity(field, points)


def _probe_groups() -> dict[str, np.ndarray]:
    return {
        "midplane_core": _rings((0.25, 0.45, 0.65), (-0.15, 0.0, 0.15)),
        "axial_interior": _rings((0.25, 0.45, 0.65), (-1.30, -0.95, 0.95, 1.30)),
        "radial_interior": _rings((0.85, 1.05, 1.25), (-0.55, 0.0, 0.55)),
        "interior_corner": _rings((0.85, 1.10), (-1.15, -0.80, 0.80, 1.15)),
        "radial_flank": _rings((1.68, 1.78, 1.88), (-0.45, 0.0, 0.45)),
        "axial_flank": _rings((0.35, 0.75, 1.15), (-1.88, -1.72, 1.72, 1.88)),
    }


def _response(field, mode, step: float, points: np.ndarray, time: float) -> np.ndarray:
    step = float(step)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    if mode == COMPACT_POLOIDAL:
        plus = _compact_velocity(field, points, time, step)
        minus = _compact_velocity(field, points, time, -step)
        return (plus - minus) / (2.0 * step)
    if mode == PHI01:
        plus = _with_phi_delta(field, PHI01, step).at_points(points, time)
        minus = _with_phi_delta(field, PHI01, -step).at_points(points, time)
        return (plus - minus) / (2.0 * step)
    if mode == PHI03:
        # The capped candidate has Phi(0,3)=-4 at the inherited bound.  Use the
        # feasible inward one-sided derivative rather than silently widening it.
        plus = _with_phi_delta(field, PHI03, step).at_points(points, time)
        base = field.at_points(points, time)
        return (plus - base) / step
    raise ValueError(f"unknown response mode {mode}")


def _response_matrix(field, groups, *, time: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    points = np.concatenate(tuple(groups.values()))
    columns = [
        _response(field, PHI01, step, points, time).reshape(-1),
        _response(field, PHI03, step, points, time).reshape(-1),
        _response(field, COMPACT_POLOIDAL, step, points, time).reshape(-1),
    ]
    return points, np.column_stack(columns)


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
        "column_norms": [float(v) for v in norms],
        "normalized_singular_values": [float(v) for v in singular],
        "normalized_rank": rank,
        "normalized_condition_number": condition,
    }


def _novelty(existing: np.ndarray, proposed: np.ndarray) -> float:
    coefficients, *_ = np.linalg.lstsq(existing, proposed, rcond=None)
    orthogonal = proposed - existing @ coefficients
    return float(np.linalg.norm(orthogonal) / max(np.linalg.norm(proposed), 1e-300))


def _cylindrical_components(points: np.ndarray, velocity: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    radius = np.hypot(points[:, 0], points[:, 1])
    cos_theta = np.divide(points[:, 0], radius, out=np.ones_like(radius), where=radius > 0.0)
    sin_theta = np.divide(points[:, 1], radius, out=np.zeros_like(radius), where=radius > 0.0)
    radial = velocity[:, 0] * cos_theta + velocity[:, 1] * sin_theta
    swirl = -velocity[:, 0] * sin_theta + velocity[:, 1] * cos_theta
    return radial, swirl, velocity[:, 2]


def _basis_structure_checks(field) -> dict[str, float]:
    positive = _rings((0.20, 0.45, 0.80, 1.20), (0.20, 0.55, 1.10), n_theta=8)
    mirrored = positive.copy()
    mirrored[:, 2] *= -1.0
    vp = _compact_poloidal_basis_velocity(field, positive)
    vm = _compact_poloidal_basis_velocity(field, mirrored)
    urp, utp, uzp = _cylindrical_components(positive, vp)
    urm, utm, uzm = _cylindrical_components(mirrored, vm)

    radial_even_error = float(np.max(np.abs(urp - urm)))
    axial_odd_error = float(np.max(np.abs(uzp + uzm)))
    swirl_rms = float(np.sqrt(np.mean(np.concatenate((utp, utm)) ** 2)))

    flank = np.concatenate((_probe_groups()["radial_flank"], _probe_groups()["axial_flank"]))
    flank_max = float(np.max(np.abs(_compact_poloidal_basis_velocity(field, flank))))

    # Independent Cartesian finite-difference divergence check away from the axis
    # and away from the compact support transition.
    check = np.array(
        [[0.31, 0.22, 0.17], [0.62, -0.27, -0.41], [0.88, 0.31, 0.67], [1.02, -0.36, -0.73]],
        dtype=float,
    )
    h = 1.0e-5
    divergence = np.zeros(len(check), dtype=float)
    for axis in range(3):
        offset = np.zeros_like(check)
        offset[:, axis] = h
        plus = _compact_poloidal_basis_velocity(field, check + offset)
        minus = _compact_poloidal_basis_velocity(field, check - offset)
        derivative = (plus[:, axis] - minus[:, axis]) / (2.0 * h)
        divergence += derivative

    return {
        "radial_even_max_abs_error": radial_even_error,
        "axial_odd_max_abs_error": axial_odd_error,
        "swirl_response_rms": swirl_rms,
        "support_flank_max_abs_velocity": flank_max,
        "cartesian_fd_divergence_max_abs": float(np.max(np.abs(divergence))),
    }


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, fraction: float) -> float:
    values = np.asarray(values, dtype=float).ravel()
    weights = np.asarray(weights, dtype=float).ravel()
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    total = float(np.sum(weights))
    if not total > 0.0:
        raise RuntimeError("weighted quantile encountered zero weight")
    cumulative = np.cumsum(weights)
    index = int(np.searchsorted(cumulative, float(fraction) * total, side="left"))
    return float(values[min(index, len(values) - 1)])


def _vorticity_morphology(
    field,
    *,
    time: float,
    mode=None,
    coefficient: float = 0.0,
    grid_size: int = 25,
    box_half_width: float = 1.95,
) -> dict[str, float]:
    if grid_size < 17 or grid_size % 2 == 0:
        raise ValueError("grid_size must be odd and >=17")
    radial_plateau, axial_plateau = _plateau_half_widths(field)
    axis = np.linspace(-float(box_half_width), float(box_half_width), int(grid_size))
    spacing = float(axis[1] - axis[0])
    x, y, z = np.meshgrid(axis, axis, axis, indexing="ij")
    points = np.column_stack((x.ravel(), y.ravel(), z.ravel()))

    if mode is None:
        velocity = field.at_points(points, float(time))
    elif mode == COMPACT_POLOIDAL:
        velocity = _compact_velocity(field, points, float(time), float(coefficient))
    elif mode in (PHI01, PHI03):
        velocity = _with_phi_delta(field, mode, float(coefficient)).at_points(points, float(time))
    else:
        raise ValueError(f"unknown morphology mode {mode}")
    velocity = np.asarray(velocity, dtype=float).reshape(grid_size, grid_size, grid_size, 3)

    u, v, w = velocity[..., 0], velocity[..., 1], velocity[..., 2]
    two_h = 2.0 * spacing
    dw_dy = (w[1:-1, 2:, 1:-1] - w[1:-1, :-2, 1:-1]) / two_h
    dv_dz = (v[1:-1, 1:-1, 2:] - v[1:-1, 1:-1, :-2]) / two_h
    du_dz = (u[1:-1, 1:-1, 2:] - u[1:-1, 1:-1, :-2]) / two_h
    dw_dx = (w[2:, 1:-1, 1:-1] - w[:-2, 1:-1, 1:-1]) / two_h
    dv_dx = (v[2:, 1:-1, 1:-1] - v[:-2, 1:-1, 1:-1]) / two_h
    du_dy = (u[1:-1, 2:, 1:-1] - u[1:-1, :-2, 1:-1]) / two_h
    omega_sq = (
        (dw_dy - dv_dz) ** 2
        + (du_dz - dw_dx) ** 2
        + (dv_dx - du_dy) ** 2
    )
    if not np.all(np.isfinite(omega_sq)):
        raise RuntimeError("vorticity reconstruction became nonfinite")
    total = float(np.sum(omega_sq))
    if not total > np.finfo(float).tiny:
        raise RuntimeError("vorticity is numerically inactive")

    xi = x[1:-1, 1:-1, 1:-1]
    yi = y[1:-1, 1:-1, 1:-1]
    zi = z[1:-1, 1:-1, 1:-1]
    radius_sq = xi * xi + yi * yi
    radial_rms = float(np.sqrt(np.sum(radius_sq * omega_sq) / total))
    axial_rms = float(np.sqrt(np.sum(zi * zi * omega_sq) / total))
    axial_q90 = _weighted_quantile(np.abs(zi), omega_sq, 0.90)
    axial_q99 = _weighted_quantile(np.abs(zi), omega_sq, 0.99)
    outer_axial = np.abs(zi) >= 0.50 * axial_plateau
    outer_fraction = float(np.sum(omega_sq[outer_axial]) / total)
    collar = (np.sqrt(radius_sq) >= radial_plateau) | (np.abs(zi) >= axial_plateau)
    collar_fraction = float(np.sum(omega_sq[collar]) / total)
    return {
        "radial_rms": radial_rms,
        "axial_rms": axial_rms,
        "aspect_z_over_r": axial_rms / max(radial_rms, 1e-300),
        "axial_q90": axial_q90,
        "axial_q99": axial_q99,
        "outer_axial_enstrophy_fraction": outer_fraction,
        "physical_plateau_collar_enstrophy_fraction": collar_fraction,
        "radial_rms_over_Rp": radial_rms / radial_plateau,
        "axial_rms_over_Zp": axial_rms / axial_plateau,
        "axial_q90_over_Zp": axial_q90 / axial_plateau,
        "axial_q99_over_Zp": axial_q99 / axial_plateau,
    }


def _morphology_vector(metrics: dict[str, float]) -> np.ndarray:
    return np.array(
        [
            metrics["radial_rms_over_Rp"],
            metrics["axial_rms_over_Zp"],
            metrics["aspect_z_over_r"],
            metrics["outer_axial_enstrophy_fraction"],
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


def audit_compact_poloidal_axial_capacity(
    *,
    time: float = 0.5,
    coefficient_steps: Iterable[float] = (0.02, 0.01),
    grid_size: int = 25,
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
    compact_novelty = _novelty(fine[:, :2], fine[:, 2])

    per_vector_rms = np.linalg.norm(fine, axis=0) / np.sqrt(len(points))
    morphology_by_step = []
    for step in steps:
        matrix = np.column_stack(
            [
                _morphology_response(field, PHI01, step, time=float(time), grid_size=grid_size),
                _morphology_response(field, PHI03, step, time=float(time), grid_size=grid_size),
                _morphology_response(field, COMPACT_POLOIDAL, step, time=float(time), grid_size=grid_size),
            ]
        )
        morphology_by_step.append(matrix)
    morph_coarse, morph_fine = morphology_by_step
    morphology_refinement = float(
        np.linalg.norm(morph_coarse - morph_fine) / max(np.linalg.norm(morph_fine), 1e-300)
    )
    morphology_rank = _rank_condition(morph_fine)

    base_morph = _vorticity_morphology(field, time=float(time), grid_size=grid_size)
    plus_morph = _vorticity_morphology(
        field,
        time=float(time),
        mode=COMPACT_POLOIDAL,
        coefficient=trial,
        grid_size=grid_size,
    )
    minus_morph = _vorticity_morphology(
        field,
        time=float(time),
        mode=COMPACT_POLOIDAL,
        coefficient=-trial,
        grid_size=grid_size,
    )

    axial_efficiency = np.abs(morph_fine[1, :]) / np.maximum(per_vector_rms, 1e-300)
    tip_efficiency = np.abs(morph_fine[3, :]) / np.maximum(per_vector_rms, 1e-300)
    radial_efficiency = np.abs(morph_fine[0, :]) / np.maximum(per_vector_rms, 1e-300)

    structure = _basis_structure_checks(field)
    radial_plateau, axial_plateau = _plateau_half_widths(field)
    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "time": float(time),
        "basis_definition": {
            "mode": COMPACT_POLOIDAL,
            "representation": "axisymmetric_vector_potential_A_theta_equivalently_streamfunction",
            "classification": "autonomous_design",
            "shape": "A_theta=.5*r*q*(1-(r/Rp)^2)^6_+*(1-(z/Zp)^2)^6_+",
            "radial_plateau_half_width": radial_plateau,
            "axial_plateau_half_width": axial_plateau,
            "shape_parameters_fitted": 0,
            "parity": "u_r even in z; u_z odd in z; u_theta=0",
        },
        "existing_comparison_modes": {
            "Phi01": "central bounded derivative",
            "Phi03": "one-sided feasible derivative because baseline coefficient is at -4 bound",
        },
        "coefficient_steps": list(steps),
        "basis_structure_checks": structure,
        "public_velocity_response": {
            **response_rank,
            "finite_difference_refinement_relative_change": response_refinement,
            "compact_novelty_outside_Phi01_Phi03_span": compact_novelty,
            "per_vector_response_rms": {
                "Phi01": float(per_vector_rms[0]),
                "Phi03": float(per_vector_rms[1]),
                COMPACT_POLOIDAL: float(per_vector_rms[2]),
            },
        },
        "vorticity_morphology": {
            "grid_size": int(grid_size),
            "fingerprint_rows": [
                "radial_rms_over_Rp",
                "axial_rms_over_Zp",
                "aspect_z_over_r",
                "outer_axial_enstrophy_fraction_abs_z_ge_half_Zp",
            ],
            "finest_response_columns": {
                "Phi01": [float(v) for v in morph_fine[:, 0]],
                "Phi03": [float(v) for v in morph_fine[:, 1]],
                COMPACT_POLOIDAL: [float(v) for v in morph_fine[:, 2]],
            },
            "response_rank": morphology_rank,
            "finite_difference_refinement_relative_change": morphology_refinement,
            "axial_rms_leverage_per_public_velocity_rms": {
                "Phi01": float(axial_efficiency[0]),
                "Phi03": float(axial_efficiency[1]),
                COMPACT_POLOIDAL: float(axial_efficiency[2]),
            },
            "outer_axial_fraction_leverage_per_public_velocity_rms": {
                "Phi01": float(tip_efficiency[0]),
                "Phi03": float(tip_efficiency[1]),
                COMPACT_POLOIDAL: float(tip_efficiency[2]),
            },
            "radial_rms_leverage_per_public_velocity_rms": {
                "Phi01": float(radial_efficiency[0]),
                "Phi03": float(radial_efficiency[1]),
                COMPACT_POLOIDAL: float(radial_efficiency[2]),
            },
            "baseline": base_morph,
            "diagnostic_plus": plus_morph,
            "diagnostic_minus": minus_morph,
            "diagnostic_trial_coefficient": trial,
            "diagnostic_axial_rms_span_over_Zp": float(
                abs(plus_morph["axial_rms_over_Zp"] - minus_morph["axial_rms_over_Zp"])
            ),
            "diagnostic_axial_q90_span_over_Zp": float(
                abs(plus_morph["axial_q90_over_Zp"] - minus_morph["axial_q90_over_Zp"])
            ),
            "diagnostic_outer_axial_enstrophy_fraction_span": float(
                abs(
                    plus_morph["outer_axial_enstrophy_fraction"]
                    - minus_morph["outer_axial_enstrophy_fraction"]
                )
            ),
            "diagnostic_radial_rms_span_over_Rp": float(
                abs(plus_morph["radial_rms_over_Rp"] - minus_morph["radial_rms_over_Rp"])
            ),
        },
        "parameter_growth_if_materialized": {
            "fixed_spatial_basis_shapes_added": 1,
            "scalar_coefficients_added": 1,
            "shape_parameters_fitted": 0,
            "temporal_degrees_added": 0,
        },
        "routing_contract": (
            "This screen only asks whether one compact parity-compatible poloidal direction supplies useful "
            "axial/vorticity morphology capacity beyond existing odd-Phi controls. Materialize it only if the "
            "measured morphology leverage is useful and a later governed optimization gives an explicit bound; "
            "otherwise keep it zero and do not continue indiscriminate basis growth."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/bipolar_compact_poloidal_capacity/report.json",
    )
    parser.add_argument("--grid-size", type=int, default=25)
    args = parser.parse_args()
    report = audit_compact_poloidal_axial_capacity(grid_size=args.grid_size)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
