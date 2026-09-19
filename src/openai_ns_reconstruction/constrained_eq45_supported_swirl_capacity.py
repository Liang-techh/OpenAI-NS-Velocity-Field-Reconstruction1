"""Post-support expression-capacity audit for existing Eq45 swirl channels.

This module is deliberately target-free.  It asks whether two *already present*
bounded profile coefficients, ``F(1,0)`` and ``F(1,2)``, remain useful and
independent controls of swirl/poloidal morphology after the production
physical-support transform.  It does not fit an image, add a basis function,
or change a saved velocity candidate.

Every diagnostic below samples Cartesian ``[u,v,w]`` only through the public
``at_points`` evaluators of the parent and support-connected child.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Iterable

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate


MODES = ((1, 0), (1, 2))
REGIONS = ("plateau", "radial_collar", "axial_collar")
DEFAULT_TIMES = (0.25, 0.50, 0.75)
DEFAULT_STEPS = (0.04, 0.02)


def _cylindrical_points(radii: Iterable[float], zs: Iterable[float]) -> np.ndarray:
    angles = (0.0, 0.5 * np.pi, np.pi, 1.5 * np.pi)
    rows = []
    for r in radii:
        for z in zs:
            for angle in angles:
                rows.append((r * np.cos(angle), r * np.sin(angle), z))
    return np.asarray(rows, dtype=float)


def _probe_sets() -> dict[str, np.ndarray]:
    # Default physical taper has an exact identity plateau at r,|z| <= 1.6
    # and support faces at r,|z| = 2.  All probes stay away from interfaces.
    return {
        "plateau": _cylindrical_points((0.60, 1.00), (-0.80, -0.35, 0.35, 0.80)),
        "radial_collar": _cylindrical_points((1.72, 1.88), (-0.90, -0.40, 0.40, 0.90)),
        "axial_collar": _cylindrical_points((0.60, 1.00), (-1.88, -1.72, 1.72, 1.88)),
    }


def _replace_swirl_mode(
    candidate: Eq45VelocityCandidate,
    mode: tuple[int, int],
    delta: float,
) -> Eq45VelocityCandidate:
    basis = candidate.profile_basis
    try:
        index = basis.mode_indices.index(mode)
    except ValueError as exc:
        raise ValueError(f"candidate does not contain swirl mode {mode}") from exc
    coefficients = list(basis.swirl_coefficients)
    coefficients[index] += float(delta)
    new_basis = replace(basis, swirl_coefficients=tuple(coefficients))
    return replace(candidate, profile_basis=new_basis)


def _perturbed_child(
    child: Eq45SupportedVelocityCandidate,
    mode: tuple[int, int],
    delta: float,
) -> Eq45SupportedVelocityCandidate:
    return Eq45SupportedVelocityCandidate(
        parent=_replace_swirl_mode(child.parent, mode, delta),
        taper=child.taper,
    )


def _cylindrical_components(points: np.ndarray, velocity: np.ndarray) -> tuple[np.ndarray, ...]:
    x = points[:, 0]
    y = points[:, 1]
    r = np.hypot(x, y)
    if np.any(r <= 0.0):
        raise ValueError("swirl-capacity probes must stay off the symmetry axis")
    u = velocity[:, 0]
    v = velocity[:, 1]
    radial = (x * u + y * v) / r
    swirl = (-y * u + x * v) / r
    axial = velocity[:, 2]
    return radial, swirl, axial


def _log_swirl_poloidal_ratio(points: np.ndarray, velocity: np.ndarray) -> float:
    radial, swirl, axial = _cylindrical_components(points, velocity)
    swirl_rms = float(np.sqrt(np.mean(swirl * swirl)))
    poloidal_rms = float(np.sqrt(np.mean(radial * radial + axial * axial)))
    if not np.isfinite(swirl_rms) or not np.isfinite(poloidal_rms):
        raise RuntimeError("nonfinite morphology norm")
    if swirl_rms <= 1e-14 or poloidal_rms <= 1e-14:
        raise RuntimeError("probe set has a degenerate swirl/poloidal channel")
    return float(np.log(swirl_rms / poloidal_rms))


def _morphology_vector(
    candidate: Eq45SupportedVelocityCandidate,
    probes: dict[str, np.ndarray],
    times: tuple[float, ...],
) -> np.ndarray:
    rows = []
    for region in REGIONS:
        points = probes[region]
        for time in times:
            velocity = candidate.at_points(points, time)
            rows.append(_log_swirl_poloidal_ratio(points, velocity))
    return np.asarray(rows, dtype=float)


def _velocity_response(
    candidate,
    plus,
    minus,
    probes: dict[str, np.ndarray],
    times: tuple[float, ...],
    step: float,
) -> dict[str, np.ndarray]:
    del candidate  # the public plus/minus evaluators are the only values consumed here
    out: dict[str, np.ndarray] = {}
    for region in REGIONS:
        pieces = []
        points = probes[region]
        for time in times:
            vp = plus.at_points(points, time)
            vm = minus.at_points(points, time)
            pieces.append((vp - vm) / (2.0 * step))
        out[region] = np.concatenate(pieces, axis=0)
    return out


def _response_matrix(
    child: Eq45SupportedVelocityCandidate,
    probes: dict[str, np.ndarray],
    times: tuple[float, ...],
    step: float,
) -> tuple[np.ndarray, dict[tuple[int, int], dict[str, np.ndarray]]]:
    columns = []
    velocity_responses: dict[tuple[int, int], dict[str, np.ndarray]] = {}
    for mode in MODES:
        plus = _perturbed_child(child, mode, step)
        minus = _perturbed_child(child, mode, -step)
        fp = _morphology_vector(plus, probes, times)
        fm = _morphology_vector(minus, probes, times)
        columns.append((fp - fm) / (2.0 * step))
        velocity_responses[mode] = _velocity_response(
            child, plus, minus, probes, times, step
        )
    return np.column_stack(columns), velocity_responses


def _parent_velocity_response(
    parent: Eq45VelocityCandidate,
    mode: tuple[int, int],
    probes: dict[str, np.ndarray],
    times: tuple[float, ...],
    step: float,
) -> dict[str, np.ndarray]:
    plus = _replace_swirl_mode(parent, mode, step)
    minus = _replace_swirl_mode(parent, mode, -step)
    return _velocity_response(parent, plus, minus, probes, times, step)


def _relative_difference(a: np.ndarray, b: np.ndarray) -> float:
    denominator = float(np.linalg.norm(b.ravel()))
    if denominator <= 1e-14:
        raise RuntimeError("reference response is numerically zero")
    return float(np.linalg.norm((a - b).ravel()) / denominator)


def audit_supported_swirl_capacity(
    child: Eq45SupportedVelocityCandidate,
    *,
    times: tuple[float, ...] = DEFAULT_TIMES,
    coefficient_steps: tuple[float, float] = DEFAULT_STEPS,
) -> dict:
    """Measure existing post-support swirl-control capacity without fitting.

    The response matrix uses the derivative of a dimensionless morphology
    feature, ``log(RMS(u_theta) / RMS(u_poloidal))``, in the identity plateau,
    radial collar, and axial collar at each requested time.
    """
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("child must be Eq45SupportedVelocityCandidate")
    times = tuple(float(t) for t in times)
    if not times or not all(np.isfinite(t) for t in times):
        raise ValueError("times must be a nonempty finite tuple")
    if any(t < child.time_start or t > child.time_end for t in times):
        raise ValueError("times must lie inside the child time interval")
    if len(coefficient_steps) != 2:
        raise ValueError("coefficient_steps must contain coarse and fine steps")
    coarse, fine = (float(v) for v in coefficient_steps)
    if not np.isfinite(coarse) or not np.isfinite(fine) or not (coarse > fine > 0.0):
        raise ValueError("coefficient steps must satisfy finite coarse > fine > 0")

    # The screen is intentionally limited to existing modes.  Reject instead of
    # silently growing the basis.
    for mode in MODES:
        if mode not in child.parent.profile_basis.mode_indices:
            raise ValueError(f"parent is missing required existing mode {mode}")
        index = child.parent.profile_basis.mode_indices.index(mode)
        value = child.parent.profile_basis.swirl_coefficients[index]
        limit = child.parent.profile_basis.coefficient_limit
        if abs(value) + coarse > limit:
            raise ValueError(f"coefficient step leaves the registered bound for mode {mode}")

    probes = _probe_sets()
    coarse_matrix, _ = _response_matrix(child, probes, times, coarse)
    fine_matrix, child_velocity = _response_matrix(child, probes, times, fine)

    singular_values = np.linalg.svd(fine_matrix, compute_uv=False)
    tolerance = np.finfo(float).eps * max(fine_matrix.shape) * singular_values[0]
    rank = int(np.count_nonzero(singular_values > tolerance))
    condition = float(singular_values[0] / singular_values[-1])
    refinement = float(
        np.linalg.norm(fine_matrix - coarse_matrix)
        / max(np.linalg.norm(fine_matrix), np.finfo(float).tiny)
    )

    f10 = fine_matrix[:, 0]
    f12 = fine_matrix[:, 1]
    projection = f10 * (float(np.dot(f10, f12)) / float(np.dot(f10, f10)))
    novelty = float(np.linalg.norm(f12 - projection) / np.linalg.norm(f12))

    region_slices = {
        region: slice(i * len(times), (i + 1) * len(times))
        for i, region in enumerate(REGIONS)
    }
    region_response_share = {}
    for column_index, mode in enumerate(MODES):
        column = fine_matrix[:, column_index]
        total = float(np.linalg.norm(column))
        region_response_share[str(mode)] = {
            region: float(np.linalg.norm(column[sl]) / total)
            for region, sl in region_slices.items()
        }

    swirl_fraction = {}
    plateau_mismatch = {}
    collar_change = {}
    for mode in MODES:
        parent_velocity = _parent_velocity_response(child.parent, mode, probes, times, fine)
        mode_key = str(mode)
        plateau_mismatch[mode_key] = _relative_difference(
            child_velocity[mode]["plateau"], parent_velocity["plateau"]
        )
        collar_change[mode_key] = {
            region: _relative_difference(child_velocity[mode][region], parent_velocity[region])
            for region in ("radial_collar", "axial_collar")
        }

        all_response = np.concatenate(
            [child_velocity[mode][region] for region in REGIONS], axis=0
        )
        all_points = np.concatenate(
            [
                np.tile(probes[region], (len(times), 1))
                for region in REGIONS
            ],
            axis=0,
        )
        radial, swirl, axial = _cylindrical_components(all_points, all_response)
        total_norm = float(np.linalg.norm(all_response.ravel()))
        theta_norm = float(np.linalg.norm(swirl))
        poloidal_norm = float(np.sqrt(np.linalg.norm(radial) ** 2 + np.linalg.norm(axial) ** 2))
        swirl_fraction[mode_key] = {
            "swirl_response_fraction": theta_norm / total_norm,
            "poloidal_cross_talk_fraction": poloidal_norm / total_norm,
        }

    return {
        "schema": "eq45_supported_swirl_capacity_v1",
        "parent_sha256": child.parent.sha256,
        "child_sha256": child.sha256,
        "modes": [list(mode) for mode in MODES],
        "parameter_count_added": 0,
        "existing_parameters_audited": 2,
        "times": list(times),
        "coefficient_steps": [coarse, fine],
        "morphology_feature": "log_rms_swirl_to_poloidal_ratio",
        "regions": list(REGIONS),
        "singular_values": [float(v) for v in singular_values],
        "numerical_rank": rank,
        "condition_number": condition,
        "response_refinement_relative_change": refinement,
        "f12_novelty_outside_f10_span": novelty,
        "response_norms": {
            str(mode): float(np.linalg.norm(fine_matrix[:, i]))
            for i, mode in enumerate(MODES)
        },
        "region_response_share": region_response_share,
        "swirl_poloidal_cross_talk": swirl_fraction,
        "plateau_supported_vs_parent_response_relative_difference": plateau_mismatch,
        "collar_supported_vs_parent_response_relative_difference": collar_change,
        "interpretation": {
            "new_basis_required_by_this_audit": False,
            "recommended_first_swirl_control": "F(1,0)",
            "recommended_axially_shaped_swirl_control_if_needed": "F(1,2)",
            "public_visual_target_fitted": False,
        },
        "truth_boundary": {
            "velocity_changed": False,
            "production_coefficients_changed": False,
            "new_basis_added": False,
            "forcing_or_pressure_refit": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
