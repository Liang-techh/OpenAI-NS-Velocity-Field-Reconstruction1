"""Screen a compact C2 temporal window for supported Eq45 Phi(1,0).

The derivative-balanced quartic Phi(1,0) candidate lowers late-time morphology
collateral, but it still carries nonzero coefficient slope at exact static-return
snapshots and its largest temporal slope occurs at the selected early endpoint.
Repeated polynomial-degree growth is not automatically the smallest response.

This module screens a different temporal family without adding a spatial mode.
The selected early displacement is multiplied by the quintic smootherstep
remainder ``q(s)=1-10s^3+15s^4-6s^5`` and then set identically to zero after a
deterministic return time.  The return is chosen target-free: it is the earliest
return whose peak ``|d(delta Phi10)/d tau|`` does not exceed the already-screened
quartic peak.  The compact trial therefore preserves the selected early snapshot,
returns exactly to the static supported field with zero first and second
coefficient derivatives, and adds no fitted temporal coefficient.

The audit compares public ``at_points(...)->[u,v,w]`` excursions against the
materialized quartic field at fixed off-keyframes.  No PDE residual, force,
pressure, public image, or hidden OpenAI parameter chooses the window.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_cubic_temporal_capacity import DEFAULT_PROBES
from .constrained_eq45_supported_phi10_quadratic_temporal_capacity import DEFAULT_EARLY_DELTA
from .constrained_eq45_supported_phi10_quartic_derivative_capacity import (
    balanced_nullspace_coefficient,
)
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)


_MODE = (1, 0)
_SMOOTHERSTEP_PEAK_SLOPE = 15.0 / 8.0
_SMOOTHERSTEP_SLOPE_ENERGY = 10.0 / 7.0
_SMOOTHERSTEP_PEAK_CURVATURE = 10.0 / np.sqrt(3.0)
DEFAULT_OFF_KEYFRAME_TIMES = (
    0.28125,
    0.3125,
    0.34375,
    0.375,
    0.4,
    0.4375,
    0.5625,
    0.6875,
)

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "screen_trial_velocity_changed": True,
    "new_spatial_basis_added": False,
    "new_fitted_temporal_parameter_added": False,
    "compact_temporal_family_screened": True,
    "force_or_pressure_fitted": False,
    "pde_objective_used_to_choose_window": False,
    "public_image_fitted": False,
    "production_temporal_shape_promoted": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _validate_early_delta(early_delta: float) -> float:
    early = float(early_delta)
    if not np.isfinite(early) or early == 0.0:
        raise ValueError("early_delta must be finite and nonzero")
    return early


def _tau(base: Eq45SupportedVelocityCandidate, time: float) -> float:
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    value = float(time)
    if not np.isfinite(value):
        raise ValueError("time must be finite")
    if value < base.time_start or value > base.time_end:
        raise ValueError("time must lie inside the candidate delivery interval")
    midpoint = 0.5 * (float(base.time_start) + float(base.time_end))
    half_width = 0.5 * (float(base.time_end) - float(base.time_start))
    return (value - midpoint) / half_width


def _time_from_tau(base: Eq45SupportedVelocityCandidate, tau: float) -> float:
    value = float(tau)
    if not np.isfinite(value) or not (-1.0 <= value <= 1.0):
        raise ValueError("tau must be finite and lie in [-1,1]")
    midpoint = 0.5 * (float(base.time_start) + float(base.time_end))
    half_width = 0.5 * (float(base.time_end) - float(base.time_start))
    return midpoint + half_width * value


def smootherstep_remainder(s: float | np.ndarray):
    """Return the C2 quintic remainder for normalized progress ``s in [0,1]``."""
    values = np.asarray(s, dtype=float)
    if not np.all(np.isfinite(values)) or np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("s must be finite and lie in [0,1]")
    return 1.0 - 10.0 * values**3 + 15.0 * values**4 - 6.0 * values**5


def smootherstep_remainder_ds(s: float | np.ndarray):
    values = np.asarray(s, dtype=float)
    if not np.all(np.isfinite(values)) or np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("s must be finite and lie in [0,1]")
    return -30.0 * values**2 + 60.0 * values**3 - 30.0 * values**4


def smootherstep_remainder_d2s(s: float | np.ndarray):
    values = np.asarray(s, dtype=float)
    if not np.all(np.isfinite(values)) or np.any((values < 0.0) | (values > 1.0)):
        raise ValueError("s must be finite and lie in [0,1]")
    return -60.0 * values + 180.0 * values**2 - 120.0 * values**3


def _quartic_delta_derivative_polynomial(
    early_delta: float, nullspace_coefficient: float
) -> np.poly1d:
    early = _validate_early_delta(early_delta)
    coefficient = float(nullspace_coefficient)
    if not np.isfinite(coefficient):
        raise ValueError("nullspace_coefficient must be finite")
    cubic_scale = early / -3.0
    return np.poly1d(
        [
            4.0 * coefficient,
            3.0 * cubic_scale - 1.5 * coefficient,
            -3.0 * cubic_scale - 2.0 * coefficient,
            0.5 * cubic_scale + 0.5 * coefficient,
        ]
    )


def _polynomial_peak_abs(poly: np.poly1d, lower: float, upper: float) -> tuple[float, float]:
    roots = np.roots(np.polyder(poly))
    locations = [float(lower), float(upper)]
    for root in roots:
        if abs(float(np.imag(root))) <= 1.0e-12:
            value = float(np.real(root))
            if lower <= value <= upper:
                locations.append(value)
    magnitudes = [abs(float(np.polyval(poly, value))) for value in locations]
    index = int(np.argmax(magnitudes))
    return float(magnitudes[index]), float(locations[index])


def _polynomial_square_integral(poly: np.poly1d, lower: float, upper: float) -> float:
    antiderivative = np.polyint(np.polymul(poly, poly))
    return float(np.polyval(antiderivative, upper) - np.polyval(antiderivative, lower))


def compact_window_parameters(*, early_delta: float = DEFAULT_EARLY_DELTA) -> dict[str, float]:
    """Return the earliest C2 return allowed by the quartic peak-slope cap."""
    early = _validate_early_delta(early_delta)
    nullspace = float(balanced_nullspace_coefficient(early_delta=early))
    quartic_derivative = _quartic_delta_derivative_polynomial(early, nullspace)
    quartic_peak_slope, quartic_peak_slope_tau = _polynomial_peak_abs(
        quartic_derivative, -1.0, 1.0
    )
    compact_numerator = abs(early) * _SMOOTHERSTEP_PEAK_SLOPE
    window_tau_length = compact_numerator / quartic_peak_slope
    if not (0.0 < window_tau_length < 1.0):
        raise ValueError(
            "slope-capped smootherstep does not return before tau=0 for this early_delta"
        )
    return_tau = -1.0 + window_tau_length

    quartic_slope_energy = _polynomial_square_integral(quartic_derivative, -1.0, 1.0)
    compact_slope_energy = early * early * _SMOOTHERSTEP_SLOPE_ENERGY / window_tau_length

    quartic_second = np.polyder(quartic_derivative)
    quartic_peak_curvature, quartic_peak_curvature_tau = _polynomial_peak_abs(
        quartic_second, -1.0, 1.0
    )
    compact_peak_curvature = (
        abs(early) * _SMOOTHERSTEP_PEAK_CURVATURE / (window_tau_length**2)
    )

    return {
        "early_delta": early,
        "quartic_nullspace_coefficient": nullspace,
        "quartic_peak_abs_d_delta_dtau": quartic_peak_slope,
        "quartic_peak_abs_d_delta_dtau_tau": quartic_peak_slope_tau,
        "window_tau_length": window_tau_length,
        "return_tau": return_tau,
        "compact_peak_abs_d_delta_dtau": compact_numerator / window_tau_length,
        "compact_to_quartic_peak_slope_ratio": compact_numerator
        / (window_tau_length * quartic_peak_slope),
        "quartic_integral_d_delta_dtau_squared": quartic_slope_energy,
        "compact_integral_d_delta_dtau_squared": compact_slope_energy,
        "compact_to_quartic_slope_energy_ratio": compact_slope_energy / quartic_slope_energy,
        "quartic_peak_abs_d2_delta_dtau2": quartic_peak_curvature,
        "quartic_peak_abs_d2_delta_dtau2_tau": quartic_peak_curvature_tau,
        "compact_peak_abs_d2_delta_dtau2": compact_peak_curvature,
        "compact_to_quartic_peak_curvature_ratio": compact_peak_curvature
        / quartic_peak_curvature,
    }


def compact_phi10_delta(
    tau: float | np.ndarray, *, early_delta: float = DEFAULT_EARLY_DELTA
):
    """Return the slope-capped C2 compact Phi(1,0) displacement."""
    early = _validate_early_delta(early_delta)
    values = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(values)) or np.any((values < -1.0) | (values > 1.0)):
        raise ValueError("tau must be finite and lie in [-1,1]")
    parameters = compact_window_parameters(early_delta=early)
    length = float(parameters["window_tau_length"])
    return_tau = float(parameters["return_tau"])
    progress = np.clip((values + 1.0) / length, 0.0, 1.0)
    active = values <= return_tau
    return np.where(active, early * smootherstep_remainder(progress), 0.0)


def compact_phi10_delta_dtau(
    tau: float | np.ndarray, *, early_delta: float = DEFAULT_EARLY_DELTA
):
    early = _validate_early_delta(early_delta)
    values = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(values)) or np.any((values < -1.0) | (values > 1.0)):
        raise ValueError("tau must be finite and lie in [-1,1]")
    parameters = compact_window_parameters(early_delta=early)
    length = float(parameters["window_tau_length"])
    return_tau = float(parameters["return_tau"])
    progress = np.clip((values + 1.0) / length, 0.0, 1.0)
    active = (values > -1.0) & (values < return_tau)
    return np.where(
        active,
        early * smootherstep_remainder_ds(progress) / length,
        0.0,
    )


def compact_phi10_delta_d2tau(
    tau: float | np.ndarray, *, early_delta: float = DEFAULT_EARLY_DELTA
):
    early = _validate_early_delta(early_delta)
    values = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(values)) or np.any((values < -1.0) | (values > 1.0)):
        raise ValueError("tau must be finite and lie in [-1,1]")
    parameters = compact_window_parameters(early_delta=early)
    length = float(parameters["window_tau_length"])
    return_tau = float(parameters["return_tau"])
    progress = np.clip((values + 1.0) / length, 0.0, 1.0)
    active = (values > -1.0) & (values < return_tau)
    return np.where(
        active,
        early * smootherstep_remainder_d2s(progress) / (length * length),
        0.0,
    )


def _phi10_midpoint(base: Eq45SupportedVelocityCandidate) -> tuple[int, float]:
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    basis = base.parent.profile_basis
    if _MODE not in basis.mode_indices:
        raise RuntimeError("base profile basis does not contain Phi(1,0)")
    index = basis.mode_indices.index(_MODE)
    return index, float(basis.phi_coefficients[index])


def compact_phi10_snapshot(
    base: Eq45SupportedVelocityCandidate,
    time: float,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
) -> Eq45SupportedVelocityCandidate:
    """Return one supported snapshot of the compact temporal trial."""
    tau = _tau(base, time)
    index, midpoint = _phi10_midpoint(base)
    delta = float(compact_phi10_delta(tau, early_delta=early_delta))
    coefficient = midpoint + delta
    limit = float(base.parent.profile_basis.coefficient_limit)
    if not np.isfinite(coefficient) or coefficient < -limit or coefficient > limit:
        raise ValueError("compact Phi(1,0) schedule exceeds the existing coefficient bound")
    phi = list(base.parent.profile_basis.phi_coefficients)
    phi[index] = coefficient
    basis = replace(base.parent.profile_basis, phi_coefficients=tuple(phi))
    parent = replace(base.parent, profile_basis=basis)
    return Eq45SupportedVelocityCandidate(parent=parent, taper=base.taper)


def _rms(values: np.ndarray) -> float:
    array = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(array * array)))


def audit_supported_phi10_c2_compact_temporal_capacity(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
    probes: np.ndarray = DEFAULT_PROBES,
    off_keyframe_times=DEFAULT_OFF_KEYFRAME_TIMES,
) -> dict[str, object]:
    """Compare the slope-capped compact window with the materialized quartic field."""
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    early = _validate_early_delta(early_delta)
    point_array = np.asarray(probes, dtype=float)
    if point_array.ndim != 2 or point_array.shape[1] != 3 or point_array.shape[0] < 3:
        raise ValueError("probes must have shape (n,3) with n>=3")
    if not np.all(np.isfinite(point_array)):
        raise ValueError("probes must be finite")

    parameters = compact_window_parameters(early_delta=early)
    return_time = _time_from_tau(child, float(parameters["return_tau"]))
    quartic = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
        base=child, early_delta=early
    )

    early_compact = compact_phi10_snapshot(child, child.time_start, early_delta=early)
    early_difference = early_compact.at_points(point_array, child.time_start) - quartic.at_points(
        point_array, child.time_start
    )

    times = tuple(float(value) for value in off_keyframe_times)
    if not times or not all(np.isfinite(value) for value in times):
        raise ValueError("off_keyframe_times must contain finite values")
    if any(value <= child.time_start or value >= child.time_end for value in times):
        raise ValueError("off_keyframe_times must lie strictly inside the delivery interval")

    rows: list[dict[str, float | bool]] = []
    for time in times:
        tau = _tau(child, time)
        static_velocity = child.at_points(point_array, time)
        quartic_velocity = quartic.at_points(point_array, time)
        compact_snapshot = compact_phi10_snapshot(child, time, early_delta=early)
        compact_velocity = compact_snapshot.at_points(point_array, time)
        quartic_rms = _rms(quartic_velocity - static_velocity)
        compact_rms = _rms(compact_velocity - static_velocity)
        rows.append(
            {
                "time": time,
                "tau": tau,
                "after_compact_return": bool(time >= return_time),
                "quartic_coefficient_delta": float(
                    quartic.coefficient_at(time) - quartic.midpoint_coefficient
                ),
                "compact_coefficient_delta": float(compact_phi10_delta(tau, early_delta=early)),
                "quartic_public_velocity_delta_rms": quartic_rms,
                "compact_public_velocity_delta_rms": compact_rms,
                "compact_to_quartic_public_velocity_delta_rms_ratio": float(
                    compact_rms / max(quartic_rms, np.finfo(float).tiny)
                ),
                "compact_public_velocity_max_abs_from_static": float(
                    np.max(np.abs(compact_velocity - static_velocity))
                ),
            }
        )

    return {
        "schema": "eq45_supported_phi10_c2_compact_temporal_capacity_v1",
        "task_id": "CR003-EQ45-SUPPORTED-PHI10-C2-COMPACT-WINDOW-CAPACITY-029",
        "claim_scope": "target_free_temporal_representation_capacity_screen",
        "base_supported_sha256": child.sha256,
        "quartic_candidate_sha256": quartic.sha256,
        "mode": {"family": "Phi", "index": [1, 0]},
        "selection_rule": (
            "earliest C2 smootherstep return whose peak abs d(delta Phi10)/d tau does not "
            "exceed the already-screened derivative-balanced quartic peak"
        ),
        "parameters": {**parameters, "return_time": return_time},
        "endpoint_checks": {
            "early_public_velocity_rms_difference_from_quartic": _rms(early_difference),
            "early_public_velocity_max_abs_difference_from_quartic": float(
                np.max(np.abs(early_difference))
            ),
            "compact_d_delta_dtau_at_start": float(
                compact_phi10_delta_dtau(-1.0, early_delta=early)
            ),
            "compact_d2_delta_dtau2_at_start": float(
                compact_phi10_delta_d2tau(-1.0, early_delta=early)
            ),
            "compact_d_delta_dtau_at_return": float(
                compact_phi10_delta_dtau(parameters["return_tau"], early_delta=early)
            ),
            "compact_d2_delta_dtau2_at_return": float(
                compact_phi10_delta_d2tau(parameters["return_tau"], early_delta=early)
            ),
        },
        "off_keyframe_public_velocity_checks": rows,
        "routing_recommendation": (
            "This compact C2 family removes all post-return coefficient and public-velocity "
            "collateral without exceeding the quartic peak slope, but it spends more integrated "
            "slope energy and peak curvature and reshapes the early transition. Treat it as the "
            "next non-polynomial temporal representation to compare visually/PDE-wise, not as a "
            "promotion and not as evidence for more spatial basis growth."
        ),
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
