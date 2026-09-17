"""Screen one quartic temporal nullspace degree for supported Eq45 Phi(1,0).

The cubic-localized Phi(1,0) schedule already preserves the selected early field
and returns exactly to the static supported field at t=.50/.625/.75.  Independent
public-output validation nevertheless shows that snapshot identity does not make
the time derivative identical at those return frames.

This module asks one representation-only question: can exactly one additional
temporal degree reduce that return-frame derivative collateral without adding a
spatial basis mode or changing any of the four already-screened snapshots?

The extra degree is the quartic nullspace

    g(tau) = (tau + 1) tau (tau - 1/2) (tau - 1),

which vanishes at tau=-1,0,1/2,1.  We add ``a*g`` to the existing cubic and
choose ``a`` by a target-free least-squares criterion that minimizes the three
coefficient slopes at the static-return nodes tau=0,1/2,1.  This is a temporal
capacity screen, not a PDE fit or public-image fit.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_cubic_temporal_capacity import (
    DEFAULT_PROBES,
    cubic_phi10_delta,
    cubic_phi10_delta_dtau,
    cubic_phi10_snapshot,
)
from .constrained_eq45_supported_phi10_quadratic_temporal_capacity import DEFAULT_EARLY_DELTA


RETURN_TAU = np.array([0.0, 0.5, 1.0], dtype=float)
TIME_DESIGN_TAU = np.linspace(-1.0, 1.0, 9)
OFF_KEYFRAME_TIMES = (0.375, 0.5625, 0.6875)

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "screen_trial_velocity_changed": True,
    "new_spatial_basis_added": False,
    "cubic_temporal_basis_dimension": 3,
    "screen_temporal_basis_dimension": 4,
    "temporal_basis_dimension_increment": 1,
    "force_or_pressure_fitted": False,
    "pde_objective_used_to_choose_nullspace_coefficient": False,
    "public_image_fitted": False,
    "production_temporal_shape_promoted": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _tau(base: Eq45SupportedVelocityCandidate, time: float) -> float:
    value = float(time)
    if not np.isfinite(value):
        raise ValueError("time must be finite")
    if value < base.time_start or value > base.time_end:
        raise ValueError("time must lie inside the candidate delivery interval")
    midpoint = 0.5 * (base.time_start + base.time_end)
    half_width = 0.5 * (base.time_end - base.time_start)
    return (value - midpoint) / half_width


def quartic_nullspace(tau: float | np.ndarray):
    tau_arr = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(tau_arr)) or np.any((tau_arr < -1.0) | (tau_arr > 1.0)):
        raise ValueError("tau must be finite and lie in [-1,1]")
    return (tau_arr + 1.0) * tau_arr * (tau_arr - 0.5) * (tau_arr - 1.0)


def quartic_nullspace_dtau(tau: float | np.ndarray):
    tau_arr = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(tau_arr)) or np.any((tau_arr < -1.0) | (tau_arr > 1.0)):
        raise ValueError("tau must be finite and lie in [-1,1]")
    return 4.0 * tau_arr**3 - 1.5 * tau_arr**2 - 2.0 * tau_arr + 0.5


def balanced_nullspace_coefficient(*, early_delta: float = DEFAULT_EARLY_DELTA) -> float:
    """Least-squares coefficient minimizing slopes at tau=0,1/2,1.

    No PDE residual or image target enters this calculation.  It only minimizes
    the Euclidean norm of the temporal profile coefficient derivative at the
    three frames where the cubic displacement itself is exactly zero.
    """
    early = float(early_delta)
    if not np.isfinite(early) or early == 0.0:
        raise ValueError("early_delta must be finite and nonzero")
    cubic_slopes = np.asarray(
        cubic_phi10_delta_dtau(RETURN_TAU, early_delta=early), dtype=float
    )
    null_slopes = np.asarray(quartic_nullspace_dtau(RETURN_TAU), dtype=float)
    denominator = float(np.dot(null_slopes, null_slopes))
    if denominator <= 0.0:
        raise RuntimeError("quartic nullspace derivative unexpectedly degenerate")
    return float(-np.dot(null_slopes, cubic_slopes) / denominator)


def quartic_phi10_delta(
    tau: float | np.ndarray,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
    nullspace_coefficient: float | None = None,
):
    early = float(early_delta)
    if not np.isfinite(early) or early == 0.0:
        raise ValueError("early_delta must be finite and nonzero")
    coefficient = (
        balanced_nullspace_coefficient(early_delta=early)
        if nullspace_coefficient is None
        else float(nullspace_coefficient)
    )
    if not np.isfinite(coefficient):
        raise ValueError("nullspace_coefficient must be finite")
    return np.asarray(cubic_phi10_delta(tau, early_delta=early), dtype=float) + coefficient * np.asarray(
        quartic_nullspace(tau), dtype=float
    )


def quartic_phi10_delta_dtau(
    tau: float | np.ndarray,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
    nullspace_coefficient: float | None = None,
):
    early = float(early_delta)
    coefficient = (
        balanced_nullspace_coefficient(early_delta=early)
        if nullspace_coefficient is None
        else float(nullspace_coefficient)
    )
    if not np.isfinite(coefficient):
        raise ValueError("nullspace_coefficient must be finite")
    return np.asarray(cubic_phi10_delta_dtau(tau, early_delta=early), dtype=float) + coefficient * np.asarray(
        quartic_nullspace_dtau(tau), dtype=float
    )


def _phi10_midpoint(base: Eq45SupportedVelocityCandidate) -> tuple[int, float]:
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    basis = base.parent.profile_basis
    mode = (1, 0)
    if mode not in basis.mode_indices:
        raise RuntimeError("base profile basis does not contain Phi(1,0)")
    index = basis.mode_indices.index(mode)
    return index, float(basis.phi_coefficients[index])


def _delta_polynomial(early_delta: float, nullspace_coefficient: float) -> np.poly1d:
    early = float(early_delta)
    a = float(nullspace_coefficient)
    cubic = (early / -3.0) * np.poly1d([1.0, -1.5, 0.5, 0.0])
    nullspace = a * np.poly1d([1.0, -0.5, -1.0, 0.5, 0.0])
    return cubic + nullspace


def _interval_extrema(
    early_delta: float,
    nullspace_coefficient: float,
    lower: float,
    upper: float,
) -> tuple[float, float, float]:
    polynomial = _delta_polynomial(early_delta, nullspace_coefficient)
    roots = np.roots(np.polyder(polynomial))
    locations = [float(lower), float(upper)]
    for root in roots:
        if abs(float(np.imag(root))) <= 1.0e-12:
            value = float(np.real(root))
            if lower <= value <= upper:
                locations.append(value)
    values = np.asarray([float(np.polyval(polynomial, value)) for value in locations])
    index = int(np.argmax(np.abs(values)))
    return float(values[index]), float(locations[index]), float(np.max(np.abs(values)))


def _schedule_extrema(
    base: Eq45SupportedVelocityCandidate,
    early_delta: float,
    nullspace_coefficient: float,
) -> dict[str, float]:
    _, midpoint = _phi10_midpoint(base)
    polynomial = _delta_polynomial(early_delta, nullspace_coefficient)
    roots = np.roots(np.polyder(polynomial))
    locations = [-1.0, 1.0]
    for root in roots:
        if abs(float(np.imag(root))) <= 1.0e-12:
            value = float(np.real(root))
            if -1.0 <= value <= 1.0:
                locations.append(value)
    deltas = np.asarray([float(np.polyval(polynomial, value)) for value in locations])
    coefficients = midpoint + deltas
    limit = float(base.parent.profile_basis.coefficient_limit)
    minimum = float(np.min(coefficients))
    maximum = float(np.max(coefficients))
    if minimum < -limit or maximum > limit:
        raise ValueError("quartic Phi(1,0) schedule exceeds the existing profile coefficient bound")
    return {
        "minimum_coefficient": minimum,
        "maximum_coefficient": maximum,
        "coefficient_limit": limit,
        "max_abs_delta": float(np.max(np.abs(deltas))),
    }


def quartic_phi10_snapshot(
    base: Eq45SupportedVelocityCandidate,
    time: float,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
    nullspace_coefficient: float | None = None,
) -> Eq45SupportedVelocityCandidate:
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    early = float(early_delta)
    a = (
        balanced_nullspace_coefficient(early_delta=early)
        if nullspace_coefficient is None
        else float(nullspace_coefficient)
    )
    _schedule_extrema(base, early, a)
    index, midpoint = _phi10_midpoint(base)
    coefficient = midpoint + float(
        quartic_phi10_delta(_tau(base, time), early_delta=early, nullspace_coefficient=a)
    )
    phi = list(base.parent.profile_basis.phi_coefficients)
    phi[index] = coefficient
    basis = replace(base.parent.profile_basis, phi_coefficients=tuple(phi))
    parent = replace(base.parent, profile_basis=basis)
    return Eq45SupportedVelocityCandidate(parent=parent, taper=base.taper)


def _rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def _time_design() -> dict[str, object]:
    tau = np.asarray(TIME_DESIGN_TAU, dtype=float)
    design = np.column_stack((tau, tau**2, tau**3, tau**4))
    singular_values = np.linalg.svd(design, compute_uv=False)
    return {
        "sample_tau": tau.tolist(),
        "basis": ["tau", "tau_squared", "tau_cubed", "tau_fourth"],
        "rank": int(np.linalg.matrix_rank(design)),
        "singular_values": singular_values.tolist(),
        "condition_number": float(singular_values[0] / singular_values[-1]),
        "cubic_temporal_basis_dimension": 3,
        "screen_temporal_basis_dimension": 4,
        "basis_dimension_increment": 1,
    }


def audit_supported_phi10_quartic_derivative_capacity(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
    probes: np.ndarray = DEFAULT_PROBES,
) -> dict[str, object]:
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    early = float(early_delta)
    a = balanced_nullspace_coefficient(early_delta=early)
    _schedule_extrema(child, early, a)

    point_array = np.asarray(probes, dtype=float)
    if point_array.ndim != 2 or point_array.shape[1] != 3 or point_array.shape[0] < 3:
        raise ValueError("probes must have shape (n,3) with n>=3")
    if not np.all(np.isfinite(point_array)):
        raise ValueError("probes must be finite")

    cubic_return_slopes = np.asarray(
        cubic_phi10_delta_dtau(RETURN_TAU, early_delta=early), dtype=float
    )
    quartic_return_slopes = np.asarray(
        quartic_phi10_delta_dtau(
            RETURN_TAU, early_delta=early, nullspace_coefficient=a
        ),
        dtype=float,
    )
    return_slope_ratio = float(
        np.linalg.norm(quartic_return_slopes) / np.linalg.norm(cubic_return_slopes)
    )
    cubic_start = float(cubic_phi10_delta_dtau(-1.0, early_delta=early))
    quartic_start = float(
        quartic_phi10_delta_dtau(-1.0, early_delta=early, nullspace_coefficient=a)
    )

    _, cubic_late_location, cubic_late_max = _interval_extrema(early, 0.0, 0.0, 1.0)
    _, quartic_late_location, quartic_late_max = _interval_extrema(early, a, 0.0, 1.0)

    anchor_rows: list[dict[str, float]] = []
    for time in (0.25, 0.50, 0.625, 0.75):
        quartic = quartic_phi10_snapshot(child, time, early_delta=early, nullspace_coefficient=a)
        if time == 0.25:
            reference = cubic_phi10_snapshot(child, time, early_delta=early)
            reference_name = "cubic"
        else:
            reference = child
            reference_name = "static"
        difference = quartic.at_points(point_array, time) - reference.at_points(point_array, time)
        anchor_rows.append(
            {
                "time": float(time),
                "tau": float(_tau(child, time)),
                "reference": reference_name,
                "public_velocity_rms_difference": _rms(difference),
                "public_velocity_max_abs_difference": float(np.max(np.abs(difference))),
            }
        )

    off_keyframe_rows: list[dict[str, float]] = []
    for time in OFF_KEYFRAME_TIMES:
        tau = _tau(child, time)
        cubic = cubic_phi10_snapshot(child, time, early_delta=early)
        quartic = quartic_phi10_snapshot(
            child, time, early_delta=early, nullspace_coefficient=a
        )
        cubic_difference = cubic.at_points(point_array, time) - child.at_points(point_array, time)
        quartic_difference = quartic.at_points(point_array, time) - child.at_points(point_array, time)
        cubic_rms = _rms(cubic_difference)
        quartic_rms = _rms(quartic_difference)
        off_keyframe_rows.append(
            {
                "time": float(time),
                "tau": float(tau),
                "cubic_coefficient_delta": float(cubic_phi10_delta(tau, early_delta=early)),
                "quartic_coefficient_delta": float(
                    quartic_phi10_delta(tau, early_delta=early, nullspace_coefficient=a)
                ),
                "cubic_public_velocity_delta_rms": cubic_rms,
                "quartic_public_velocity_delta_rms": quartic_rms,
                "quartic_to_cubic_public_velocity_delta_rms_ratio": float(
                    quartic_rms / max(cubic_rms, np.finfo(float).tiny)
                ),
            }
        )

    return {
        "schema": "eq45_supported_phi10_quartic_derivative_capacity_v1",
        "task_id": "CR003-EQ45-SUPPORTED-PHI10-QUARTIC-DERIVATIVE-CAPACITY-025",
        "base_candidate_sha": child.sha256,
        "mode": {"family": "Phi", "index": [1, 0]},
        "early_delta": early,
        "nullspace_coefficient": a,
        "selection_objective": (
            "minimize the Euclidean norm of d(delta Phi10)/d tau at the three "
            "static-return nodes tau=0,0.5,1; no PDE residual or image target"
        ),
        "time_design": _time_design(),
        "schedule_extrema": _schedule_extrema(child, early, a),
        "return_frame_derivative_capacity": {
            "tau": RETURN_TAU.tolist(),
            "cubic_d_delta_dtau": cubic_return_slopes.tolist(),
            "quartic_d_delta_dtau": quartic_return_slopes.tolist(),
            "quartic_to_cubic_l2_ratio": return_slope_ratio,
            "fractional_l2_reduction": 1.0 - return_slope_ratio,
        },
        "early_endpoint_temporal_cost": {
            "cubic_d_delta_dtau": cubic_start,
            "quartic_d_delta_dtau": quartic_start,
            "quartic_to_cubic_abs_ratio": abs(quartic_start) / abs(cubic_start),
        },
        "late_half_localization": {
            "cubic_max_abs_delta": cubic_late_max,
            "cubic_max_abs_delta_tau": cubic_late_location,
            "quartic_max_abs_delta": quartic_late_max,
            "quartic_max_abs_delta_tau": quartic_late_location,
            "quartic_to_cubic_ratio": quartic_late_max / cubic_late_max,
        },
        "anchor_public_velocity_checks": anchor_rows,
        "off_keyframe_public_velocity_checks": off_keyframe_rows,
        "routing_recommendation": (
            "One more temporal degree has substantial localization headroom without spatial basis "
            "growth, but it trades that gain for a sharper early endpoint slope. If public visual "
            "fingerprints need still tighter early-only localization, this quartic direction is the "
            "next minimal representation to test; otherwise stop temporal growth and keep the cubic."
        ),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
