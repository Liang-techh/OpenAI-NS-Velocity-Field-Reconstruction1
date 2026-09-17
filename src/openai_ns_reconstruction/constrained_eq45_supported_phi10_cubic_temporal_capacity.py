"""Screen one extra temporal degree for a more localized supported Eq45 Phi(1,0) trial.

The existing early-localized quadratic schedule preserves the successful early
Phi(1,0) endpoint and returns to the static supported field at the midpoint and
late endpoint.  Independent morphology audit #154 nevertheless measures a
small, resolved collateral at t=0.625, where the quadratic coefficient still
has a +0.175 excursion.

This module asks one minimal representation question: can exactly one more
*temporal* degree remove that known intermediate collateral without adding any
spatial basis term?  The screened cubic displacement is

    delta c(tau) = d_early * tau * (tau - 0.5) * (tau - 1) / (-3),

so delta c(-1)=d_early while delta c(0)=delta c(0.5)=delta c(1)=0.  For the
already-selected d_early=-1.4 it therefore keeps the successful t=0.25 field
but returns exactly to the static supported field at t=0.50, 0.625, and 0.75.

This is target-free capacity evidence only.  It does not fit a public image,
change the canonical candidate, validate the PDE, or identify an OpenAI field.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_phi10_quadratic_temporal_capacity import (
    DEFAULT_EARLY_DELTA,
    quadratic_phi10_delta,
    quadratic_phi10_snapshot,
)


DEFAULT_PROBE_TIMES = (0.25, 0.375, 0.50, 0.625, 0.6875, 0.75)
DEFAULT_PROBES = np.array(
    [
        [0.37, 0.11, -0.22],
        [1.10, -0.25, 0.42],
        [1.72, 0.00, 0.30],
        [1.68, 0.18, 1.72],
        [0.70, -0.60, -1.72],
        [-1.72, 0.10, 0.60],
    ],
    dtype=float,
)

_TRUTH_BOUNDARY = {
    "canonical_velocity_changed": False,
    "screen_trial_velocity_changed": True,
    "new_spatial_basis_added": False,
    "quadratic_temporal_basis_dimension": 2,
    "screen_temporal_basis_dimension": 3,
    "temporal_basis_dimension_increment": 1,
    "production_temporal_shape_promoted": False,
    "openai_image_fitted": False,
    "physical_support_validated": False,
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


def cubic_phi10_delta(
    tau: float | np.ndarray,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
):
    """Return the cubic Phi(1,0) displacement with roots at tau=0,0.5,1."""
    early = float(early_delta)
    if not np.isfinite(early) or early == 0.0:
        raise ValueError("early_delta must be finite and nonzero")
    tau_arr = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(tau_arr)) or np.any((tau_arr < -1.0) | (tau_arr > 1.0)):
        raise ValueError("tau must be finite and lie in [-1,1]")
    product = tau_arr * (tau_arr - 0.5) * (tau_arr - 1.0)
    return early * (product / -3.0)


def cubic_phi10_delta_dtau(
    tau: float | np.ndarray,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
):
    """Return d(delta Phi10)/d tau for the cubic screen."""
    early = float(early_delta)
    if not np.isfinite(early) or early == 0.0:
        raise ValueError("early_delta must be finite and nonzero")
    tau_arr = np.asarray(tau, dtype=float)
    if not np.all(np.isfinite(tau_arr)) or np.any((tau_arr < -1.0) | (tau_arr > 1.0)):
        raise ValueError("tau must be finite and lie in [-1,1]")
    return (early / -3.0) * (3.0 * tau_arr * tau_arr - 3.0 * tau_arr + 0.5)


def _quadratic_delta_dtau(tau: float, early_delta: float) -> float:
    return float(0.5 * early_delta * (2.0 * float(tau) - 1.0))


def _phi10_midpoint_coefficient(base: Eq45SupportedVelocityCandidate) -> tuple[int, float]:
    basis = base.parent.profile_basis
    mode = (1, 0)
    if mode not in basis.mode_indices:
        raise RuntimeError("base profile basis does not contain Phi(1,0)")
    index = basis.mode_indices.index(mode)
    return index, float(basis.phi_coefficients[index])


def cubic_phi10_coefficient(
    base: Eq45SupportedVelocityCandidate,
    time: float,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
) -> float:
    """Return the bounded Phi(1,0) coefficient for the cubic screen."""
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    _, midpoint = _phi10_midpoint_coefficient(base)
    return midpoint + float(cubic_phi10_delta(_tau(base, time), early_delta=early_delta))


def _schedule_extrema(
    base: Eq45SupportedVelocityCandidate,
    early_delta: float,
) -> dict[str, object]:
    _, midpoint = _phi10_midpoint_coefficient(base)
    root_lo = (3.0 - np.sqrt(3.0)) / 6.0
    root_hi = (3.0 + np.sqrt(3.0)) / 6.0
    tau_values = np.array([-1.0, 0.0, root_lo, 0.5, root_hi, 1.0], dtype=float)
    deltas = np.asarray(cubic_phi10_delta(tau_values, early_delta=early_delta), dtype=float)
    coefficients = midpoint + deltas
    limit = float(base.parent.profile_basis.coefficient_limit)
    lo = float(np.min(coefficients))
    hi = float(np.max(coefficients))
    if lo < -limit or hi > limit:
        raise ValueError("cubic Phi(1,0) schedule exceeds the existing profile coefficient bound")
    return {
        "sample_tau": tau_values.tolist(),
        "sample_delta": deltas.tolist(),
        "minimum_coefficient": lo,
        "maximum_coefficient": hi,
        "coefficient_limit": limit,
        "late_half_max_abs_delta": float(np.max(np.abs(deltas[2:]))),
    }


def cubic_phi10_snapshot(
    base: Eq45SupportedVelocityCandidate,
    time: float,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
) -> Eq45SupportedVelocityCandidate:
    """Materialize one static supported snapshot of the cubic screen."""
    if not isinstance(base, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    early = float(early_delta)
    _schedule_extrema(base, early)
    index, _ = _phi10_midpoint_coefficient(base)
    coefficient = cubic_phi10_coefficient(base, float(time), early_delta=early)
    phi = list(base.parent.profile_basis.phi_coefficients)
    phi[index] = coefficient
    basis = replace(base.parent.profile_basis, phi_coefficients=tuple(phi))
    parent = replace(base.parent, profile_basis=basis)
    return Eq45SupportedVelocityCandidate(parent=parent, taper=base.taper)


def _time_design() -> dict[str, object]:
    tau = np.array([-1.0, -0.5, 0.0, 0.5, 1.0], dtype=float)
    design = np.column_stack((tau, tau * tau, tau * tau * tau))
    singular_values = np.linalg.svd(design, compute_uv=False)
    return {
        "sample_tau": tau.tolist(),
        "basis": ["tau", "tau_squared", "tau_cubed"],
        "rank": int(np.linalg.matrix_rank(design)),
        "singular_values": singular_values.tolist(),
        "condition_number": float(singular_values[0] / singular_values[-1]),
        "quadratic_temporal_basis_dimension": 2,
        "screen_temporal_basis_dimension": 3,
        "basis_dimension_increment": 1,
    }


def _delta_metrics(
    candidate: Eq45SupportedVelocityCandidate,
    baseline: Eq45SupportedVelocityCandidate,
    probes: np.ndarray,
    time: float,
) -> dict[str, float]:
    delta = candidate.at_points(probes, time) - baseline.at_points(probes, time)
    return {
        "rms": float(np.sqrt(np.mean(delta * delta))),
        "max_abs": float(np.max(np.abs(delta))),
    }


def audit_supported_phi10_cubic_temporal_capacity(
    base: Eq45SupportedVelocityCandidate | None = None,
    *,
    early_delta: float = DEFAULT_EARLY_DELTA,
    probe_times: tuple[float, ...] = DEFAULT_PROBE_TIMES,
    probes: np.ndarray = DEFAULT_PROBES,
) -> dict[str, object]:
    """Compare the existing quadratic schedule with one cubic temporal degree."""
    child = governed_supported_seed() if base is None else base
    if not isinstance(child, Eq45SupportedVelocityCandidate):
        raise TypeError("base must be Eq45SupportedVelocityCandidate")
    early = float(early_delta)
    if not np.isfinite(early) or early == 0.0:
        raise ValueError("early_delta must be finite and nonzero")

    times = tuple(float(value) for value in probe_times)
    if len(times) < 6 or not all(np.isfinite(value) for value in times):
        raise ValueError("probe_times must contain at least six finite times")
    if times != tuple(sorted(times)):
        raise ValueError("probe_times must be increasing")
    if times[0] != child.time_start or times[-1] != child.time_end:
        raise ValueError("probe_times must span the full delivery interval")
    if 0.50 not in times or 0.625 not in times:
        raise ValueError("probe_times must include t=0.50 and t=0.625")

    point_array = np.asarray(probes, dtype=float)
    if point_array.ndim != 2 or point_array.shape[1] != 3 or point_array.shape[0] < 3:
        raise ValueError("probes must have shape (n,3) with n>=3")
    if not np.all(np.isfinite(point_array)):
        raise ValueError("probes must be finite")

    extrema = _schedule_extrema(child, early)
    rows: list[dict[str, object]] = []
    for time in times:
        tau = _tau(child, time)
        quadratic = quadratic_phi10_snapshot(child, time, early_delta=early)
        cubic = cubic_phi10_snapshot(child, time, early_delta=early)
        quadratic_delta = _delta_metrics(quadratic, child, point_array, time)
        cubic_delta = _delta_metrics(cubic, child, point_array, time)
        cross = cubic.at_points(point_array, time) - quadratic.at_points(point_array, time)
        quadratic_rms = quadratic_delta["rms"]
        rows.append(
            {
                "time": time,
                "tau": tau,
                "quadratic_coefficient_delta": float(quadratic_phi10_delta(tau, early_delta=early)),
                "cubic_coefficient_delta": float(cubic_phi10_delta(tau, early_delta=early)),
                "quadratic_public_velocity_delta": quadratic_delta,
                "cubic_public_velocity_delta": cubic_delta,
                "cubic_to_quadratic_public_velocity_delta_rms_ratio": (
                    None if quadratic_rms <= 1.0e-15 else float(cubic_delta["rms"] / quadratic_rms)
                ),
                "cubic_quadratic_public_velocity_max_abs_difference": float(
                    np.max(np.abs(cross))
                ),
            }
        )

    row_by_time = {float(row["time"]): row for row in rows}
    early_match = row_by_time[child.time_start]
    static_return = {
        str(time): {
            "cubic_public_velocity_delta_rms": row_by_time[time][
                "cubic_public_velocity_delta"
            ]["rms"],
            "cubic_public_velocity_delta_max_abs": row_by_time[time][
                "cubic_public_velocity_delta"
            ]["max_abs"],
        }
        for time in (0.50, 0.625, child.time_end)
    }

    quadratic_late_max = abs(float(quadratic_phi10_delta(0.5, early_delta=early)))
    cubic_late_max = float(extrema["late_half_max_abs_delta"])
    quadratic_start_derivative = abs(_quadratic_delta_dtau(-1.0, early))
    cubic_start_derivative = abs(float(cubic_phi10_delta_dtau(-1.0, early_delta=early)))

    return {
        "schema": "eq45_supported_phi10_cubic_temporal_capacity_v1",
        "task_id": "CR003-EQ45-SUPPORTED-PHI10-CUBIC-TEMPORAL-CAPACITY-022",
        "base_candidate_sha": child.candidate_sha,
        "mode": {"family": "Phi", "index": [1, 0]},
        "early_delta": early,
        "time_design": _time_design(),
        "schedule_extrema": extrema,
        "probe_rows": rows,
        "early_endpoint": {
            "cubic_quadratic_public_velocity_max_abs_difference": early_match[
                "cubic_quadratic_public_velocity_max_abs_difference"
            ]
        },
        "static_return": static_return,
        "late_half_localization": {
            "quadratic_max_abs_coefficient_delta": quadratic_late_max,
            "cubic_max_abs_coefficient_delta": cubic_late_max,
            "cubic_to_quadratic_ratio": cubic_late_max / quadratic_late_max,
        },
        "temporal_sharpness": {
            "quadratic_abs_d_delta_dtau_at_start": quadratic_start_derivative,
            "cubic_abs_d_delta_dtau_at_start": cubic_start_derivative,
            "cubic_to_quadratic_start_derivative_ratio": (
                cubic_start_derivative / quadratic_start_derivative
            ),
            "pde_implication": (
                "The cubic localizes morphology more aggressively but increases the start-endpoint "
                "coefficient time derivative; PDE residual must be independently rerun if this shape "
                "is materialized."
            ),
        },
        "visualization_fingerprint_implication": {
            "source": "PR #154 t=0.625 whole-domain morphology audit",
            "field_identity_statement": (
                "Because the cubic Phi(1,0) displacement is exactly zero at t=0.625, its full "
                "supported [u,v,w] field is exactly the static baseline there; therefore every "
                "deterministic morphology fingerprint of that field returns to the static value."
            ),
            "quadratic_known_intermediate_collateral_removed_by_identity": True,
        },
        "pde_validation_rerun": False,
        "truth_boundary": dict(_TRUTH_BOUNDARY),
    }
