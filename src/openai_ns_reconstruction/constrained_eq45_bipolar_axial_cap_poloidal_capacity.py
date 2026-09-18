"""Target-free support-adjacent axial-cap poloidal capacity screen.

The preceding axial-shoulder screen showed that multiplying the same centered
compact envelope by higher odd powers of q moves response energy outward but
still leaves the gross q90/q99 vorticity reach unchanged.  This increment tests
one different minimal direction: move the *support of the response itself* into
fixed symmetric axial-cap bands while keeping a curl-generated, axis-regular
poloidal correction.

With s=(r/Rp)^2 and q=z/Zp, define

    y = (q^2-c2)/d2,
    g(q) = q*(1-y^2)^6_+,
    A_theta = .5*r*(1-s)^6_+*g(q),

where the fixed autonomous band is 0.55 < |q| < 0.95.  Thus the correction is
zero around the midplane, is concentrated immediately inside the existing
physical plateau edge, and is also zero before the physical taper collar starts.
No residual map, public image, or fitted shape parameter chooses this band.

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

from .constrained_eq45_bipolar_axial_shoulder_poloidal_capacity import (
    AXIAL_SHOULDER_POLOIDAL,
    _axial_shoulder_basis_velocity,
    _morphology_response,
    _morphology_vector,
    _shoulder_response,
    _vorticity_morphology,
)
from .constrained_eq45_bipolar_compact_poloidal_capacity import (
    COMPACT_POLOIDAL,
    PHI01,
    PHI03,
    _compact_poloidal_basis_velocity,
    _cylindrical_components,
    _novelty,
    _probe_groups,
    _rank_condition,
    _response,
    _source_field,
)
from .constrained_eq45_bipolar_interior_compact_swirl_capacity import _plateau_half_widths, _rings

TASK_ID = "CR003-BIPOLAR-AXIAL-CAP-POLOIDAL-CAPACITY-043"
AXIAL_CAP_POLOIDAL = "AXIAL_CAP_BANDED_C4_ODD_Z_POLOIDAL"
CAP_Q_INNER = 0.55
CAP_Q_OUTER = 0.95
CAP_Q2_CENTER = 0.5 * (CAP_Q_INNER**2 + CAP_Q_OUTER**2)
CAP_Q2_HALF_WIDTH = 0.5 * (CAP_Q_OUTER**2 - CAP_Q_INNER**2)
TRUTH_BOUNDARY = {
    "velocity_changed": False,
    "candidate_artifact_changed": False,
    "new_basis_materialized": False,
    "axial_cap_poloidal_coefficient_selected": False,
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


def _cap_envelope(q: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return g(q), g'(q), and the fixed symmetric axial-cap activity mask."""
    q = np.asarray(q, dtype=float)
    y = (q * q - CAP_Q2_CENTER) / CAP_Q2_HALF_WIDTH
    active = np.abs(y) < 1.0
    one_minus_y2 = np.where(active, 1.0 - y * y, 0.0)
    h5 = one_minus_y2**5
    h6 = one_minus_y2**6
    g = q * h6
    # y' = 2q/d2 and d[q(1-y^2)^6]/dq
    # = (1-y^2)^6 - 24 q^2 y/d2 (1-y^2)^5.
    gprime = np.where(
        active,
        h6 - 24.0 * q * q * y * h5 / CAP_Q2_HALF_WIDTH,
        0.0,
    )
    return g, gprime, active


def _axial_cap_basis_velocity(field, points) -> np.ndarray:
    """Evaluate curl(A_theta e_theta) for the fixed support-adjacent cap mode."""
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
    one_minus_s = np.where(radial_inside, 1.0 - s, 0.0)
    ar5 = one_minus_s**5
    ar6 = one_minus_s**6
    g, gprime, axial_active = _cap_envelope(q)
    active = radial_inside & axial_active

    # For A_theta=.5*r*Ar(s)*g(q):
    # u_r=-d_z A_theta, u_z=(1/r)d_r(r A_theta).
    u_r = -0.5 * radius * ar6 * gprime / axial_plateau
    u_z = ar5 * (1.0 - 7.0 * s) * g
    u_r = np.where(active, u_r, 0.0)
    u_z = np.where(active, u_z, 0.0)

    cos_theta = np.divide(x, radius, out=np.ones_like(radius), where=radius > 0.0)
    sin_theta = np.divide(y, radius, out=np.zeros_like(radius), where=radius > 0.0)
    out = np.column_stack((u_r * cos_theta, u_r * sin_theta, u_z))
    if not np.all(np.isfinite(out)):
        raise RuntimeError("axial-cap poloidal basis became nonfinite")
    return out


def _cap_velocity(field, points, time: float, coefficient: float) -> np.ndarray:
    coefficient = float(coefficient)
    if not np.isfinite(coefficient):
        raise ValueError("axial-cap coefficient must be finite")
    # Diagnostic implementation guard only; this does not select a materialized bound.
    limit = float(field.parent.profile_basis.coefficient_limit)
    if abs(coefficient) > limit:
        raise ValueError("axial-cap coefficient exceeds inherited implementation guard")
    points = np.asarray(points, dtype=float)
    return field.at_points(points, float(time)) + coefficient * _axial_cap_basis_velocity(field, points)


def _cap_response(field, step: float, points: np.ndarray, time: float) -> np.ndarray:
    step = float(step)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    plus = _cap_velocity(field, points, time, step)
    minus = _cap_velocity(field, points, time, -step)
    return (plus - minus) / (2.0 * step)


def _response_matrix(field, groups, *, time: float, step: float) -> tuple[np.ndarray, np.ndarray]:
    points = np.concatenate(tuple(groups.values()))
    columns = [
        _response(field, PHI01, step, points, time).reshape(-1),
        _response(field, PHI03, step, points, time).reshape(-1),
        _response(field, COMPACT_POLOIDAL, step, points, time).reshape(-1),
        _shoulder_response(field, step, points, time).reshape(-1),
        _cap_response(field, step, points, time).reshape(-1),
    ]
    return points, np.column_stack(columns)


def _basis_structure_checks(field) -> dict[str, float]:
    # These z values lie well inside the fixed cap band for Zp=1.6.
    positive = _rings((0.20, 0.45, 0.80, 1.20), (0.95, 1.15, 1.35), n_theta=8)
    mirrored = positive.copy()
    mirrored[:, 2] *= -1.0
    vp = _axial_cap_basis_velocity(field, positive)
    vm = _axial_cap_basis_velocity(field, mirrored)
    urp, utp, uzp = _cylindrical_components(positive, vp)
    urm, utm, uzm = _cylindrical_components(mirrored, vm)

    flank = np.concatenate((_probe_groups()["radial_flank"], _probe_groups()["axial_flank"]))
    midplane = _rings((0.25, 0.55, 0.85, 1.15), (0.0,), n_theta=8)
    inner_axial = _rings((0.25, 0.55, 0.85, 1.15), (-0.70, -0.35, 0.35, 0.70), n_theta=8)

    check = np.array(
        [[0.31, 0.22, 1.02], [0.62, -0.27, -1.08], [0.88, 0.31, 1.24], [1.02, -0.36, -1.30]],
        dtype=float,
    )
    h = 1.0e-5
    divergence = np.zeros(len(check), dtype=float)
    for axis in range(3):
        offset = np.zeros_like(check)
        offset[:, axis] = h
        plus = _axial_cap_basis_velocity(field, check + offset)
        minus = _axial_cap_basis_velocity(field, check - offset)
        divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * h)

    return {
        "radial_even_max_abs_error": float(np.max(np.abs(urp - urm))),
        "axial_odd_max_abs_error": float(np.max(np.abs(uzp + uzm))),
        "swirl_response_rms": float(np.sqrt(np.mean(np.concatenate((utp, utm)) ** 2))),
        "support_flank_max_abs_velocity": float(np.max(np.abs(_axial_cap_basis_velocity(field, flank)))),
        "midplane_max_abs_velocity": float(np.max(np.abs(_axial_cap_basis_velocity(field, midplane)))),
        "inner_axial_band_max_abs_velocity": float(np.max(np.abs(_axial_cap_basis_velocity(field, inner_axial)))),
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

    def metrics(response: np.ndarray) -> tuple[float, float, float]:
        energy = np.sum(response * response, axis=1)
        total = float(np.sum(energy))
        if total <= np.finfo(float).tiny:
            raise RuntimeError("basis response is numerically inactive")
        abs_z = np.abs(points[:, 2])
        centroid = float(np.sum(abs_z * energy) / total / axial_plateau)
        outer065 = float(np.sum(energy[abs_z >= 0.65 * axial_plateau]) / total)
        outer075 = float(np.sum(energy[abs_z >= 0.75 * axial_plateau]) / total)
        return centroid, outer065, outer075

    center = metrics(_compact_poloidal_basis_velocity(field, points))
    shoulder = metrics(_axial_shoulder_basis_velocity(field, points))
    cap = metrics(_axial_cap_basis_velocity(field, points))
    return {
        "center_compact_abs_z_energy_centroid_over_Zp": center[0],
        "shoulder_abs_z_energy_centroid_over_Zp": shoulder[0],
        "cap_abs_z_energy_centroid_over_Zp": cap[0],
        "cap_to_shoulder_centroid_ratio": cap[0] / max(shoulder[0], 1e-300),
        "center_outer_065_response_energy_fraction": center[1],
        "shoulder_outer_065_response_energy_fraction": shoulder[1],
        "cap_outer_065_response_energy_fraction": cap[1],
        "cap_to_shoulder_outer_065_ratio": cap[1] / max(shoulder[1], 1e-300),
        "center_outer_075_response_energy_fraction": center[2],
        "shoulder_outer_075_response_energy_fraction": shoulder[2],
        "cap_outer_075_response_energy_fraction": cap[2],
    }


class _CapFieldProxy:
    """Delegate field metadata while adding one cap coefficient to at_points."""

    def __init__(self, field, coefficient: float):
        self._field = field
        self._coefficient = float(coefficient)

    def __getattr__(self, name):
        return getattr(self._field, name)

    def at_points(self, points, time: float):
        return _cap_velocity(self._field, points, time, self._coefficient)


def _cap_vorticity_morphology(field, *, time: float, coefficient: float, grid_size: int) -> dict[str, float]:
    proxy = _CapFieldProxy(field, coefficient)
    return _vorticity_morphology(proxy, time=time, grid_size=grid_size)


def _cap_morphology_response(field, step: float, *, time: float, grid_size: int) -> np.ndarray:
    plus = _morphology_vector(_cap_vorticity_morphology(field, time=time, coefficient=step, grid_size=grid_size))
    minus = _morphology_vector(_cap_vorticity_morphology(field, time=time, coefficient=-step, grid_size=grid_size))
    return (plus - minus) / (2.0 * step)


def audit_axial_cap_poloidal_capacity(
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
    if grid_size < 17 or grid_size % 2 == 0:
        raise ValueError("grid_size must be odd and >=17")

    field = _source_field()
    groups = _probe_groups()
    points, coarse = _response_matrix(field, groups, time=float(time), step=steps[0])
    _, fine = _response_matrix(field, groups, time=float(time), step=steps[1])
    response_refinement = float(np.linalg.norm(coarse - fine) / max(np.linalg.norm(fine), 1e-300))
    response_rank = _rank_condition(fine)
    cap_novelty = _novelty(fine[:, :4], fine[:, 4])
    per_vector_rms = np.linalg.norm(fine, axis=0) / np.sqrt(len(points))

    existing_modes = (PHI01, PHI03, COMPACT_POLOIDAL, AXIAL_SHOULDER_POLOIDAL)
    morphology_by_step = []
    for step in steps:
        existing = [
            _morphology_response(field, mode, step, time=float(time), grid_size=grid_size)
            for mode in existing_modes
        ]
        existing.append(_cap_morphology_response(field, step, time=float(time), grid_size=grid_size))
        morphology_by_step.append(np.column_stack(existing))
    morph_coarse, morph_fine = morphology_by_step
    morphology_refinement = float(
        np.linalg.norm(morph_coarse - morph_fine) / max(np.linalg.norm(morph_fine), 1e-300)
    )
    morphology_rank = _rank_condition(morph_fine)

    base_morph = _vorticity_morphology(field, time=float(time), grid_size=grid_size)
    plus_morph = _cap_vorticity_morphology(field, time=float(time), coefficient=trial, grid_size=grid_size)
    minus_morph = _cap_vorticity_morphology(field, time=float(time), coefficient=-trial, grid_size=grid_size)

    axial_efficiency = np.abs(morph_fine[1, :]) / np.maximum(per_vector_rms, 1e-300)
    outer065_efficiency = np.abs(morph_fine[4, :]) / np.maximum(per_vector_rms, 1e-300)
    outer075_efficiency = np.abs(morph_fine[5, :]) / np.maximum(per_vector_rms, 1e-300)
    radial_efficiency = np.abs(morph_fine[0, :]) / np.maximum(per_vector_rms, 1e-300)
    radial_plateau, axial_plateau = _plateau_half_widths(field)
    names = ("Phi01", "Phi03", COMPACT_POLOIDAL, AXIAL_SHOULDER_POLOIDAL, AXIAL_CAP_POLOIDAL)

    return {
        "task_id": TASK_ID,
        "source_candidate_sha256": field.sha256,
        "time": float(time),
        "basis_definition": {
            "mode": AXIAL_CAP_POLOIDAL,
            "representation": "axisymmetric_vector_potential_A_theta_equivalently_streamfunction",
            "classification": "autonomous_design",
            "shape": "A_theta=.5*r*(1-(r/Rp)^2)^6_+*q*(1-y(q)^2)^6_+",
            "y_of_q": "(q^2-c2)/d2",
            "fixed_abs_q_band": [CAP_Q_INNER, CAP_Q_OUTER],
            "q2_center": CAP_Q2_CENTER,
            "q2_half_width": CAP_Q2_HALF_WIDTH,
            "radial_plateau_half_width": radial_plateau,
            "axial_plateau_half_width": axial_plateau,
            "shape_parameters_fitted": 0,
            "parity": "u_r even in z; u_z odd in z; u_theta=0",
            "design_reason": (
                "move response support into symmetric axial-cap bands instead of multiplying the same centered "
                "envelope by another odd polynomial"
            ),
        },
        "coefficient_steps": list(steps),
        "basis_structure_checks": _basis_structure_checks(field),
        "basis_axial_localization": _basis_axial_localization(field),
        "public_velocity_response": {
            **response_rank,
            "finite_difference_refinement_relative_change": response_refinement,
            "cap_novelty_outside_Phi01_Phi03_center_compact_shoulder_span": cap_novelty,
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
            "diagnostic_axial_rms_span_over_Zp": float(
                abs(plus_morph["axial_rms_over_Zp"] - minus_morph["axial_rms_over_Zp"])
            ),
            "diagnostic_axial_q90_span_over_Zp": float(
                abs(plus_morph["axial_q90_over_Zp"] - minus_morph["axial_q90_over_Zp"])
            ),
            "diagnostic_axial_q99_span_over_Zp": float(
                abs(plus_morph["axial_q99_over_Zp"] - minus_morph["axial_q99_over_Zp"])
            ),
            "diagnostic_outer_065_span": float(
                abs(plus_morph["outer_065_enstrophy_fraction"] - minus_morph["outer_065_enstrophy_fraction"])
            ),
            "diagnostic_outer_075_span": float(
                abs(plus_morph["outer_075_enstrophy_fraction"] - minus_morph["outer_075_enstrophy_fraction"])
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
            "If this fixed cap-localized direction materially moves q90/q99 or outer-tip morphology while remaining "
            "well conditioned, hand Agent 1/2 only this one bounded coefficient for candidate-specific PDE/energy "
            "screening. If q90/q99 still remain effectively immobile despite moving response support into the caps, "
            "stop additive interior/cap basis growth for the gross tip problem and investigate the base support/taper "
            "geometry instead. Do not respond by adding higher odd-q degree or more swirl modes."
        ),
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="artifacts/bipolar_axial_cap_poloidal_capacity/report.json")
    parser.add_argument("--grid-size", type=int, default=33)
    args = parser.parse_args()
    report = audit_axial_cap_poloidal_capacity(grid_size=args.grid_size)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
