"""Curl/time-commutation diagnostic for the frozen Kokuno Agent-2 velocity.

This module does **not** create or retune an oscillatory candidate.  It consumes
only the already-frozen public Agent-2 providers

``velocity_osc(x,y,z,t)`` and ``velocity_osc_dt(x,y,z,t)``

and checks the numerical seam

    d_t(curl u_osc) = curl(d_t u_osc).

The two public fields inherit the corrected Kokuno localized-wave / complete-curl
structure, while the finite-difference operators below are repository diagnostics.
The frozen time modulation is repository-autonomous candidate data, not a recovered
Kokuno hidden time law.  None of the errors reported here is a Navier--Stokes
momentum residual or an independent Agent-4 validation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_oscillatory_time_derivative import velocity_osc_dt
from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A2-PUBLIC-VORTICITY-TIME-COMMUTATION-046"
SCHEMA = "kokuno-a2-public-vorticity-time-commutation-v1"
PARENT_AGENT2_PR = 669
PARENT_AGENT2_HEAD = "f40436f2eddb77db94b07921fefa1e09ab94b134"
VELOCITY_PR = 561
VELOCITY_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
VELOCITY_DT_PR = 579
VELOCITY_DT_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"

# Deliberately different from #669's centered FD6 spatial diagnostic.
_FD4_OFFSETS = (-2, -1, 1, 2)
_FD4_COEFFICIENTS = (1.0, -8.0, 8.0, -1.0)
_FD4_DENOMINATOR = 12.0

# Independent time differentiation of the *numerically reconstructed curl*.
_FD6_TIME_OFFSETS = (-3, -2, -1, 1, 2, 3)
_FD6_TIME_COEFFICIENTS = (-1.0, 9.0, -45.0, 45.0, -9.0, 1.0)
_FD6_TIME_DENOMINATOR = 60.0

Provider = Callable[[Any, Any, Any, Any], np.ndarray]


def _broadcast_xyzt(x: Any, y: Any, z: Any, t: Any) -> tuple[np.ndarray, ...]:
    arrays = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    if not all(np.all(np.isfinite(a)) for a in arrays):
        raise ValueError("x, y, z and t must be finite and broadcastable")
    return tuple(arrays)


def _validate_spatial_step(step: float) -> float:
    h = float(step)
    if not math.isfinite(h) or not (1.0e-5 <= h <= 5.0e-2):
        raise ValueError("spatial_step must be finite and lie in [1e-5, 5e-2]")
    return h


def _validate_time_step(step: float) -> float:
    dt = float(step)
    if not math.isfinite(dt) or not (1.0e-5 <= dt <= 5.0e-2):
        raise ValueError("time_step must be finite and lie in [1e-5, 5e-2]")
    return dt


def _provider_value(
    provider: Provider,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
) -> np.ndarray:
    value = np.asarray(provider(x, y, z, t), dtype=float)
    expected = x.shape + (3,)
    if value.shape != expected:
        raise RuntimeError(
            f"vector provider returned shape {value.shape}, expected {expected}"
        )
    if not np.all(np.isfinite(value)):
        raise RuntimeError("vector provider returned non-finite values")
    return value


def _fd4_velocity_gradient(
    provider: Provider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    """Return ``[..., component, derivative_axis]`` with centered FD4."""
    h = _validate_spatial_step(spatial_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)
    base = (xb, yb, zb)
    derivatives: list[np.ndarray] = []
    for axis in range(3):
        accum = np.zeros(xb.shape + (3,), dtype=float)
        for offset, coefficient in zip(
            _FD4_OFFSETS, _FD4_COEFFICIENTS, strict=True
        ):
            shifted = [a for a in base]
            shifted[axis] = shifted[axis] + offset * h
            accum += coefficient * _provider_value(
                provider, shifted[0], shifted[1], shifted[2], tb
            )
        derivatives.append(accum / (_FD4_DENOMINATOR * h))
    return np.stack(derivatives, axis=-1)


def _curl_fd4(
    provider: Provider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
) -> np.ndarray:
    grad = _fd4_velocity_gradient(
        provider, x, y, z, t, spatial_step=spatial_step
    )
    return np.stack(
        (
            grad[..., 2, 1] - grad[..., 1, 2],
            grad[..., 0, 2] - grad[..., 2, 0],
            grad[..., 1, 0] - grad[..., 0, 1],
        ),
        axis=-1,
    )


def _time_fd6_of_curl_fd4(
    provider: Provider,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float,
    time_step: float,
) -> np.ndarray:
    """Apply centered time-FD6 to an FD4 Cartesian curl of ``provider``."""
    h = _validate_spatial_step(spatial_step)
    dt = _validate_time_step(time_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    field = default_field()
    if np.any(tb - 3.0 * dt < field.time_min) or np.any(
        tb + 3.0 * dt > field.time_max
    ):
        raise ValueError("time_step stencil leaves the registered candidate interval")

    accum = np.zeros(xb.shape + (3,), dtype=float)
    for offset, coefficient in zip(
        _FD6_TIME_OFFSETS, _FD6_TIME_COEFFICIENTS, strict=True
    ):
        accum += coefficient * _curl_fd4(
            provider,
            xb,
            yb,
            zb,
            tb + offset * dt,
            spatial_step=h,
        )
    return accum / (_FD6_TIME_DENOMINATOR * dt)


def evaluate_vorticity_time_commutation(
    x: Any,
    y: Any,
    z: Any,
    t: Any,
    *,
    spatial_step: float = 0.0045,
    time_step: float = 0.004,
) -> dict[str, np.ndarray | float | str]:
    """Evaluate ``d_t curl(u)`` versus ``curl(u_t)`` for the frozen candidate.

    Route A uses centered time-FD6 applied to a centered Cartesian FD4 curl of
    public ``velocity_osc``.  Route B applies the same FD4 curl to the exact
    frozen-candidate public ``velocity_osc_dt``.  The returned commutator is a
    derivative-consistency diagnostic, not a PDE residual.
    """
    h = _validate_spatial_step(spatial_step)
    dt = _validate_time_step(time_step)
    xb, yb, zb, tb = _broadcast_xyzt(x, y, z, t)

    time_of_curl = _time_fd6_of_curl_fd4(
        velocity_osc,
        xb,
        yb,
        zb,
        tb,
        spatial_step=h,
        time_step=dt,
    )
    curl_of_time = _curl_fd4(
        velocity_osc_dt,
        xb,
        yb,
        zb,
        tb,
        spatial_step=h,
    )
    commutator = time_of_curl - curl_of_time
    return {
        "vorticity_dt_time_fd6_of_curl_fd4": time_of_curl,
        "vorticity_dt_curl_fd4_of_velocity_dt": curl_of_time,
        "commutator": commutator,
        "spatial_step": h,
        "time_step": dt,
        "spatial_operator": "centered_cartesian_fd4",
        "time_operator": "centered_time_fd6",
    }


def _vector_rms(value: np.ndarray) -> float:
    a = np.asarray(value, dtype=float)
    return float(np.sqrt(np.mean(np.sum(a * a, axis=-1))))


def _new_offgrid_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Fresh support-interior cloud, distinct from #579 and #669 receipts."""
    radius = np.asarray(
        [
            0.38,
            0.47,
            0.55,
            0.64,
            0.72,
            0.80,
            0.88,
            0.97,
            1.05,
            1.13,
            0.42,
            0.51,
            0.60,
            0.75,
            0.84,
            0.93,
            1.02,
            1.15,
        ],
        dtype=float,
    )
    angle = np.asarray(
        [
            0.29,
            1.03,
            1.71,
            2.43,
            3.02,
            -2.29,
            -1.61,
            -0.88,
            -0.27,
            0.48,
            1.26,
            2.04,
            2.71,
            -2.76,
            -1.98,
            -1.24,
            -0.57,
            0.76,
        ],
        dtype=float,
    )
    z = np.asarray(
        [
            -1.11,
            -0.77,
            -0.42,
            -0.07,
            0.29,
            0.66,
            1.01,
            1.24,
            -1.00,
            -0.64,
            -0.26,
            0.13,
            0.48,
            0.81,
            1.07,
            -1.27,
            0.59,
            0.03,
        ],
        dtype=float,
    )
    t = np.asarray(([0.35] * 6) + ([0.50] * 6) + ([0.65] * 6), dtype=float)
    return radius * np.cos(angle), radius * np.sin(angle), z, t


def build_receipt() -> dict[str, Any]:
    """Build the preregistered curl/time-commutation and spatial-stability receipt."""
    x, y, z, t = _new_offgrid_cloud()
    spatial_step = 0.0045
    time_steps = (0.016, 0.008, 0.004)

    time_levels = [
        evaluate_vorticity_time_commutation(
            x,
            y,
            z,
            t,
            spatial_step=spatial_step,
            time_step=dt,
        )
        for dt in time_steps
    ]
    reference = np.asarray(
        time_levels[-1]["vorticity_dt_curl_fd4_of_velocity_dt"], dtype=float
    )
    reference_rms = _vector_rms(reference)
    reference_max = float(np.max(np.linalg.norm(reference, axis=-1)))

    time_rows: list[dict[str, float]] = []
    absolute_rms_errors: list[float] = []
    for dt, level in zip(time_steps, time_levels, strict=True):
        error = np.asarray(level["commutator"], dtype=float)
        error_rms = _vector_rms(error)
        error_max = float(np.max(np.linalg.norm(error, axis=-1)))
        absolute_rms_errors.append(error_rms)
        time_rows.append(
            {
                "time_step": dt,
                "absolute_rms_error": error_rms,
                "absolute_max_error": error_max,
                "relative_rms_error": error_rms / max(reference_rms, 1.0e-300),
                "relative_max_error": error_max / max(reference_max, 1.0e-300),
            }
        )
    time_refinement = [
        absolute_rms_errors[0] / max(absolute_rms_errors[1], np.finfo(float).tiny),
        absolute_rms_errors[1] / max(absolute_rms_errors[2], np.finfo(float).tiny),
    ]

    spatial_steps = (0.018, 0.009, 0.0045)
    spatial_values = [
        _curl_fd4(velocity_osc_dt, x, y, z, t, spatial_step=h)
        for h in spatial_steps
    ]
    spatial_diff_rms = [
        _vector_rms(spatial_values[0] - spatial_values[1]),
        _vector_rms(spatial_values[1] - spatial_values[2]),
    ]
    spatial_refinement = spatial_diff_rms[0] / max(
        spatial_diff_rms[1], np.finfo(float).tiny
    )

    exterior = evaluate_vorticity_time_commutation(
        np.asarray([0.0, 1.60, 0.60]),
        np.asarray([0.0, 0.0, 0.0]),
        np.asarray([0.0, 0.0, 2.16]),
        np.asarray([0.50, 0.50, 0.50]),
        spatial_step=spatial_step,
        time_step=time_steps[-1],
    )
    exterior_abs_max = max(
        float(np.max(np.abs(np.asarray(exterior[key], dtype=float))))
        for key in (
            "vorticity_dt_time_fd6_of_curl_fd4",
            "vorticity_dt_curl_fd4_of_velocity_dt",
            "commutator",
        )
    )

    guards = {
        "vorticity_dt_nontrivial_rms_min": 1.0e-8,
        "finest_commutator_relative_rms_max": 5.0e-7,
        "finest_commutator_relative_max_max": 1.0e-6,
        "time_fd6_refinement_ratio_min": 20.0,
        "spatial_fd4_refinement_ratio_min": 8.0,
        "support_exterior_absolute_max": 1.0e-12,
    }
    failed: list[str] = []
    if reference_rms < guards["vorticity_dt_nontrivial_rms_min"]:
        failed.append("vorticity_dt_nontrivial_rms")
    if time_rows[-1]["relative_rms_error"] > guards[
        "finest_commutator_relative_rms_max"
    ]:
        failed.append("finest_commutator_relative_rms")
    if time_rows[-1]["relative_max_error"] > guards[
        "finest_commutator_relative_max_max"
    ]:
        failed.append("finest_commutator_relative_max")
    if any(ratio < guards["time_fd6_refinement_ratio_min"] for ratio in time_refinement):
        failed.append("time_fd6_refinement_ratio")
    if spatial_refinement < guards["spatial_fd4_refinement_ratio_min"]:
        failed.append("spatial_fd4_refinement_ratio")
    if exterior_abs_max > guards["support_exterior_absolute_max"]:
        failed.append("support_exterior_absolute_max")

    identity_payload = {
        "schema": SCHEMA,
        "task": TASK,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "velocity_pr": VELOCITY_PR,
        "velocity_head": VELOCITY_HEAD,
        "velocity_dt_pr": VELOCITY_DT_PR,
        "velocity_dt_head": VELOCITY_DT_HEAD,
        "spatial_operator": "centered_cartesian_fd4",
        "time_operator": "centered_time_fd6",
        "spatial_step_for_commutator": spatial_step,
        "time_steps": list(time_steps),
        "spatial_stability_steps": list(spatial_steps),
    }
    identity_sha256 = hashlib.sha256(
        json.dumps(identity_payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    return {
        **identity_payload,
        "identity_sha256": identity_sha256,
        "sample_count": int(t.size),
        "verification_times": sorted({float(v) for v in t}),
        "vorticity_dt_reference_rms": reference_rms,
        "vorticity_dt_reference_sampled_max": reference_max,
        "time_commutation_ladder": {
            "rows": time_rows,
            "successive_rms_refinement_ratios": time_refinement,
        },
        "spatial_stability_ladder": {
            "steps": list(spatial_steps),
            "successive_difference_rms": spatial_diff_rms,
            "successive_difference_refinement_ratio": spatial_refinement,
        },
        "support_exterior_absolute_max": exterior_abs_max,
        "guards": guards,
        "failed_guards": failed,
        "self_diagnostic_passed": failed == [],
        "provenance": {
            "source_structure": (
                "inherits corrected-Kokuno localized oscillatory complete-curl structure "
                "from the frozen upstream public velocity"
            ),
            "time_modulation": "repository_autonomous_candidate_choice",
            "fd4_spatial_operator": "repository_numerical_diagnostic",
            "fd6_time_operator": "repository_numerical_diagnostic",
        },
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "velocity_time_derivative_changed": False,
            "oscillatory_coefficients_retuned": False,
            "source_formula_changed": False,
            "numerical_curl_time_commutation_diagnostic_only": True,
            "independent_agent4_validation_replaced": False,
            "independent_agent4_vector_potential_audit_required": True,
            "agent3_momentum_defect_or_correction_duplicated": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = build_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
