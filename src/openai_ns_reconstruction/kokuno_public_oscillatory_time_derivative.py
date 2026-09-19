"""Exact time derivative for the frozen admitted Kokuno A2 oscillatory candidate.

The admitted public field remains the Agent-2 #561 provider
``kokuno_public_z_pullback_velocity.velocity_osc``.  This module does not retune
or rewrite that field.  It differentiates only the repository-autonomous time
modulation already frozen in the candidate.

For that candidate the source covariance inverse gives

    y_+ = (n-c)/(2 h_+),   y_- = (n+c)/(2 h_-),   a_sigma=sqrt(y_sigma),

because ``T_N=-n A_c`` and ``T_K=c u_star``.  The source phase, support,
complete-curl basis, h_sigma, background, and public-z pullback are all time
independent in this frozen realization.  Each sign-resolved complete-curl field
is therefore linear in ``a_sigma(t)`` and

    d_t u_sigma = (d_t a_sigma / a_sigma) u_sigma,
    d_t a_+/a_+ = .5 (n'-c')/(n-c),
    d_t a_-/a_- = .5 (n'+c')/(n+c).

The n(t), c(t) modulation itself is a repository-autonomous numerical choice,
not a recovered Kokuno hidden time law.  The resulting derivative is exact for
this candidate only; it is not a paper-exact derivative or an NS validation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_z_pullback_velocity import default_field, velocity_osc

TASK = "KOKUNO-A2-PUBLIC-TIME-DERIVATIVE-041"
SCHEMA = "kokuno-a2-public-oscillatory-time-derivative-v1"
PARENT_AGENT2_PR = 570
PARENT_AGENT2_HEAD = "134095115fe88e4e7c285a02d8592fd81451b44c"
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
FD6_STEPS = (0.01, 0.005, 0.0025)

VelocityEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def target_log_amplitude_derivative(t: Any) -> np.ndarray:
    """Return exact ``(a_+'/a_+, a_-'/a_-)`` for the frozen candidate time law."""
    field = default_field()
    tt = _finite(t, "t")
    if np.any(tt < field.time_min) or np.any(tt > field.time_max):
        raise ValueError("t is outside the registered candidate interval")

    duration = field.time_max - field.time_min
    tau = (tt - field.time_min) / duration
    angle = 2.0 * math.pi * tau
    omega = 2.0 * math.pi / duration

    normal = field.normal_target_ratio * (
        1.0 + field.time_modulation * np.sin(angle)
    )
    cross = field.cross_target_ratio * (
        1.0 + field.time_modulation * np.cos(angle)
    )
    d_normal = (
        field.normal_target_ratio
        * field.time_modulation
        * omega
        * np.cos(angle)
    )
    d_cross = (
        -field.cross_target_ratio
        * field.time_modulation
        * omega
        * np.sin(angle)
    )

    plus_gap = normal - cross
    minus_gap = normal + cross
    if np.any(plus_gap <= 0.0) or np.any(minus_gap <= 0.0):
        raise RuntimeError("frozen covariance target lost its positive two-sign gap")

    return np.stack(
        (
            0.5 * (d_normal - d_cross) / plus_gap,
            0.5 * (d_normal + d_cross) / minus_gap,
        ),
        axis=-1,
    )


def evaluate_velocity_osc_dt(x: Any, y: Any, z: Any, t: Any) -> dict[str, Any]:
    """Evaluate the exact candidate time derivative by beta/sign and in total.

    No finite difference is used in this production derivative.  The public
    sign-resolved field returned by the admitted provider is multiplied by the
    exact logarithmic derivative of its signed covariance amplitude.
    """
    field = default_field()
    evaluated = field.evaluate(x, y, z, t)
    by_sign = np.asarray(evaluated["velocity_cartesian_by_beta_sign"], dtype=float)
    sample_shape = by_sign.shape[:-3]
    tt = np.broadcast_to(_finite(t, "t"), sample_shape)
    log_derivative = target_log_amplitude_derivative(tt)

    dt_by_sign = by_sign * log_derivative[..., None, :, None]
    dt_by_beta = np.sum(dt_by_sign, axis=-2)
    dt_total = np.sum(dt_by_beta, axis=-2)

    return {
        "velocity_dt_cartesian_total": dt_total,
        "velocity_dt_cartesian_by_beta": dt_by_beta,
        "velocity_dt_cartesian_by_beta_sign": dt_by_sign,
        "amplitude_log_derivative_by_sign": log_derivative,
        "support_mask": np.asarray(evaluated["support_mask"], dtype=bool),
        "beta_labels": evaluated["beta_labels"],
        "sign_labels": evaluated["sign_labels"],
        "coordinate_contract": evaluated["coordinate_contract"],
        "time_interval": evaluated["time_interval"],
        "derivative_contract": (
            "exact derivative of the frozen repository-autonomous covariance-target "
            "time modulation; source phase/support/complete-curl basis are time independent"
        ),
        "velocity_candidate_changed": False,
        "source_formula_changed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def velocity_osc_dt(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return ``partial_t velocity_osc`` as a Cartesian ``[...,3]`` array."""
    return evaluate_velocity_osc_dt(x, y, z, t)["velocity_dt_cartesian_total"]


def _fd6_time_derivative(
    evaluator: VelocityEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    """Independent centered sixth-order time derivative for verification only."""
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")

    f_m3 = np.asarray(evaluator(x, y, z, t - 3.0 * step), dtype=float)
    f_m2 = np.asarray(evaluator(x, y, z, t - 2.0 * step), dtype=float)
    f_m1 = np.asarray(evaluator(x, y, z, t - step), dtype=float)
    f_p1 = np.asarray(evaluator(x, y, z, t + step), dtype=float)
    f_p2 = np.asarray(evaluator(x, y, z, t + 2.0 * step), dtype=float)
    f_p3 = np.asarray(evaluator(x, y, z, t + 3.0 * step), dtype=float)
    return (
        -f_m3 + 9.0 * f_m2 - 45.0 * f_m1
        + 45.0 * f_p1 - 9.0 * f_p2 + f_p3
    ) / (60.0 * step)


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Deterministic support-interior cloud, separate from Agent-4 audit points."""
    radii = (0.43, 0.68, 0.93, 1.17)
    z_values = (-1.08, -0.31, 0.74)
    times = (0.32, 0.50, 0.68)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    index = 0
    for time in times:
        for z in z_values:
            for radius in radii:
                theta = 0.219 + (index + 1) * golden
                rows.append(
                    (
                        radius * math.cos(theta),
                        radius * math.sin(theta),
                        z,
                        time,
                    )
                )
                index += 1
    array = np.asarray(rows, dtype=float)
    return array[:, 0], array[:, 1], array[:, 2], array[:, 3]


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def verification_receipt() -> dict[str, Any]:
    """Build a target-free FD6 check of the analytic candidate derivative."""
    x, y, z, t = _verification_cloud()
    exact = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float)
    exact_rms = _vector_rms(exact)
    exact_max = float(np.max(np.linalg.norm(exact, axis=-1)))
    rows = []
    for step in FD6_STEPS:
        numerical = _fd6_time_derivative(velocity_osc, x, y, z, t, step)
        error = numerical - exact
        error_rms = _vector_rms(error)
        error_max = float(np.max(np.linalg.norm(error, axis=-1)))
        rows.append(
            {
                "step": step,
                "absolute_rms_error": error_rms,
                "absolute_max_error": error_max,
                "relative_rms_error": error_rms / max(exact_rms, 1.0e-300),
                "relative_max_error": error_max / max(exact_max, 1.0e-300),
            }
        )

    field = default_field()
    exterior = np.asarray(
        velocity_osc_dt(
            np.asarray((0.0, 1.8, 0.8)),
            np.asarray((0.0, 0.0, 0.0)),
            np.asarray((0.0, 0.0, 2.1)),
            np.asarray((0.5, 0.5, 0.5)),
        ),
        dtype=float,
    )
    exterior_max = float(np.max(np.abs(exterior)))

    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "admitted_agent2_pr": ADMITTED_AGENT2_PR,
        "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
        "public_provider": "openai_ns_reconstruction.kokuno_public_z_pullback_velocity:velocity_osc",
        "derivative_provider": "openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative:velocity_osc_dt",
        "sample_count": int(t.size),
        "verification_times": sorted({float(v) for v in t}),
        "fd6_steps": list(FD6_STEPS),
        "analytic_derivative_rms": exact_rms,
        "analytic_derivative_sampled_max": exact_max,
        "fd6_comparison": rows,
        "support_exterior_absolute_max": exterior_max,
        "time_law": {
            "normal_target_ratio": field.normal_target_ratio,
            "cross_target_ratio": field.cross_target_ratio,
            "time_modulation": field.time_modulation,
            "time_interval": [field.time_min, field.time_max],
            "provenance": "repository_autonomous_candidate_time_modulation",
        },
        "derivation": (
            "T_N=-n*A_c, T_K=c*u_star imply y_plus=(n-c)/(2*h_plus), "
            "y_minus=(n+c)/(2*h_minus); frozen basis gives dt(u_sigma)=(dt a_sigma/a_sigma)u_sigma"
        ),
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "source_formula_changed": False,
            "autonomous_time_law_derivative_executable": True,
            "independent_fd6_time_difference_compared": True,
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
    receipt = verification_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
