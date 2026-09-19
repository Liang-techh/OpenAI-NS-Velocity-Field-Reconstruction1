"""Independent Agent-4 audit of the public Kokuno oscillatory time derivative.

This module treats both public APIs as black boxes:

    velocity_osc(x,y,z,t)
    velocity_osc_dt(x,y,z,t)

It does not read the signed covariance amplitudes, phase/curl basis, target-law
helper, or Agent-2 verification FD6 implementation.  A fresh held-out cloud is
checked with an independent centered FD4 time derivative, and a separate
Gauss--Legendre fundamental-theorem check integrates the public derivative over
an interval and compares it with the public velocity endpoint difference.

This is derivative-interface evidence needed by a future independent momentum
evaluator.  It is not a full Navier--Stokes residual and does not change the
project gates (normalized momentum max/L2 <= 1e-3; divergence max/L2 <= 1e-5).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
from numpy.polynomial.legendre import leggauss

from .kokuno_public_oscillatory_time_derivative import velocity_osc_dt
from .kokuno_public_z_pullback_velocity import velocity_osc

TASK = "KOKUNO-A4-PUBLIC-TIME-DERIVATIVE-AUDIT-042"
SCHEMA = "kokuno-agent4-public-time-derivative-audit-v1"
AUDITED_AGENT2_PR = 579
AUDITED_AGENT2_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
ADMITTED_VELOCITY_AGENT2_PR = 561
ADMITTED_VELOCITY_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
PARENT_AGENT4_PR = 572
PARENT_AGENT4_HEAD = "9385af81ea4f4e3b2b2dfdf3862534a15a8b72e8"
SEED = 9173261
FD4_STEPS = (0.012, 0.006, 0.003)
QUADRATURE_ORDERS = (16, 32, 64)
TIME_INTEGRAL_INTERVAL = (0.29, 0.71)

# Frozen before the exact-head Actions evaluation.  These are interface guards,
# not substitutes for the final full-candidate PDE acceptance thresholds.
GUARDS = {
    "finest_fd4_relative_rms": 5.0e-6,
    "finest_fd4_relative_max": 1.0e-5,
    "minimum_fd4_refinement_ratio": 8.0,
    "finest_integral_relative_rms": 1.0e-9,
    "finest_integral_relative_max": 1.0e-9,
    "support_exterior_absolute_max": 1.0e-12,
    "minimum_derivative_rms": 1.0e-8,
    "minimum_sign_flip_integral_mismatch": 5.0e-1,
}

VelocityEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _fd4_time_derivative(
    evaluator: VelocityEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    """Centered fourth-order derivative using only public field values."""
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    f_m2 = np.asarray(evaluator(x, y, z, t - 2.0 * step), dtype=float)
    f_m1 = np.asarray(evaluator(x, y, z, t - step), dtype=float)
    f_p1 = np.asarray(evaluator(x, y, z, t + step), dtype=float)
    f_p2 = np.asarray(evaluator(x, y, z, t + 2.0 * step), dtype=float)
    return (f_m2 - 8.0 * f_m1 + 8.0 * f_p1 - f_p2) / (12.0 * step)


def _heldout_cloud(count: int = 36) -> tuple[np.ndarray, np.ndarray]:
    """Fresh off-grid support-interior spacetime points for local FD4 checks."""
    rng = np.random.default_rng(SEED)
    r = rng.uniform(0.36, 1.14, size=count)
    theta = rng.uniform(-math.pi, math.pi, size=count)
    z = rng.uniform(-1.22, 1.22, size=count)
    t = rng.uniform(0.31, 0.69, size=count)
    points = np.stack((r * np.cos(theta), r * np.sin(theta), z), axis=-1)
    return points, t


def _integral_points(count: int = 12) -> np.ndarray:
    """Deterministic spatial points, disjoint from the local FD4 cloud."""
    rng = np.random.default_rng(SEED + 97)
    r = rng.uniform(0.40, 1.08, size=count)
    theta = rng.uniform(-math.pi, math.pi, size=count)
    z = rng.uniform(-1.05, 1.05, size=count)
    return np.stack((r * np.cos(theta), r * np.sin(theta), z), axis=-1)


def _integrate_public_dt(points: np.ndarray, order: int) -> np.ndarray:
    """Integrate the public derivative in time with independent Gauss--Legendre."""
    if order < 2:
        raise ValueError("order must be at least two")
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n,3)")

    lo, hi = TIME_INTEGRAL_INTERVAL
    nodes, weights = leggauss(order)
    mid = 0.5 * (lo + hi)
    half = 0.5 * (hi - lo)
    times = mid + half * nodes

    n = points.shape[0]
    x = np.repeat(points[:, 0], order)
    y = np.repeat(points[:, 1], order)
    z = np.repeat(points[:, 2], order)
    t = np.tile(times, n)
    values = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float).reshape(n, order, 3)
    return half * np.einsum("j,ijc->ic", weights, values)


def _support_exterior_absolute_max() -> float:
    x = np.asarray((0.0, 0.08, 1.55, 0.72, 0.91, -0.63), dtype=float)
    y = np.asarray((0.0, 0.00, 0.00, 0.00, 0.00, 0.00), dtype=float)
    z = np.asarray((0.0, 0.4, -0.3, 2.08, -2.21, 2.35), dtype=float)
    t = np.asarray((0.33, 0.41, 0.49, 0.55, 0.63, 0.67), dtype=float)
    values = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float)
    return float(np.max(np.abs(values), initial=0.0))


def run_audit() -> dict[str, Any]:
    points, times = _heldout_cloud()
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    exact = np.asarray(velocity_osc_dt(x, y, z, times), dtype=float)
    exact_rms = _vector_rms(exact)
    exact_max = float(np.max(np.linalg.norm(exact, axis=-1), initial=0.0))

    fd_rows: list[dict[str, float]] = []
    for step in FD4_STEPS:
        numerical = _fd4_time_derivative(velocity_osc, x, y, z, times, step)
        error = numerical - exact
        abs_rms = _vector_rms(error)
        abs_max = float(np.max(np.linalg.norm(error, axis=-1), initial=0.0))
        fd_rows.append(
            {
                "step": float(step),
                "absolute_rms_error": abs_rms,
                "absolute_max_error": abs_max,
                "relative_rms_error": abs_rms / max(exact_rms, 1.0e-300),
                "relative_max_error": abs_max / max(exact_max, 1.0e-300),
            }
        )
    fd_refinement = [
        fd_rows[i]["relative_rms_error"]
        / max(fd_rows[i + 1]["relative_rms_error"], 1.0e-300)
        for i in range(len(fd_rows) - 1)
    ]

    integral_points = _integral_points()
    lo, hi = TIME_INTEGRAL_INTERVAL
    n = integral_points.shape[0]
    start = np.asarray(
        velocity_osc(
            integral_points[:, 0], integral_points[:, 1], integral_points[:, 2], np.full(n, lo)
        ),
        dtype=float,
    )
    end = np.asarray(
        velocity_osc(
            integral_points[:, 0], integral_points[:, 1], integral_points[:, 2], np.full(n, hi)
        ),
        dtype=float,
    )
    delta = end - start
    delta_rms = _vector_rms(delta)
    delta_max = float(np.max(np.linalg.norm(delta, axis=-1), initial=0.0))

    integral_rows: list[dict[str, float]] = []
    finest_integral: np.ndarray | None = None
    for order in QUADRATURE_ORDERS:
        integral = _integrate_public_dt(integral_points, order)
        mismatch = integral - delta
        mismatch_rms = _vector_rms(mismatch)
        mismatch_max = float(np.max(np.linalg.norm(mismatch, axis=-1), initial=0.0))
        integral_rows.append(
            {
                "order": int(order),
                "absolute_rms_mismatch": mismatch_rms,
                "absolute_max_mismatch": mismatch_max,
                "relative_rms_mismatch": mismatch_rms / max(delta_rms, 1.0e-300),
                "relative_max_mismatch": mismatch_max / max(delta_max, 1.0e-300),
            }
        )
        finest_integral = integral
    assert finest_integral is not None

    # External sign-flip mutation of the public derivative integral.  If the
    # checker is sensitive, -integral cannot satisfy the endpoint identity.
    sign_flip_mismatch = _vector_rms(-finest_integral - delta) / max(delta_rms, 1.0e-300)
    support_max = _support_exterior_absolute_max()
    finest_fd = fd_rows[-1]
    finest_integral_row = integral_rows[-1]

    passed = bool(
        finest_fd["relative_rms_error"] <= GUARDS["finest_fd4_relative_rms"]
        and finest_fd["relative_max_error"] <= GUARDS["finest_fd4_relative_max"]
        and min(fd_refinement) >= GUARDS["minimum_fd4_refinement_ratio"]
        and finest_integral_row["relative_rms_mismatch"]
        <= GUARDS["finest_integral_relative_rms"]
        and finest_integral_row["relative_max_mismatch"]
        <= GUARDS["finest_integral_relative_max"]
        and support_max <= GUARDS["support_exterior_absolute_max"]
        and exact_rms >= GUARDS["minimum_derivative_rms"]
        and sign_flip_mismatch >= GUARDS["minimum_sign_flip_integral_mismatch"]
    )

    failed: list[str] = []
    if finest_fd["relative_rms_error"] > GUARDS["finest_fd4_relative_rms"]:
        failed.append("finest_fd4_relative_rms")
    if finest_fd["relative_max_error"] > GUARDS["finest_fd4_relative_max"]:
        failed.append("finest_fd4_relative_max")
    if min(fd_refinement) < GUARDS["minimum_fd4_refinement_ratio"]:
        failed.append("minimum_fd4_refinement_ratio")
    if finest_integral_row["relative_rms_mismatch"] > GUARDS["finest_integral_relative_rms"]:
        failed.append("finest_integral_relative_rms")
    if finest_integral_row["relative_max_mismatch"] > GUARDS["finest_integral_relative_max"]:
        failed.append("finest_integral_relative_max")
    if support_max > GUARDS["support_exterior_absolute_max"]:
        failed.append("support_exterior_absolute_max")
    if exact_rms < GUARDS["minimum_derivative_rms"]:
        failed.append("minimum_derivative_rms")
    if sign_flip_mismatch < GUARDS["minimum_sign_flip_integral_mismatch"]:
        failed.append("minimum_sign_flip_integral_mismatch")

    return {
        "task": TASK,
        "schema": SCHEMA,
        "audited_agent2_pr": AUDITED_AGENT2_PR,
        "audited_agent2_head": AUDITED_AGENT2_HEAD,
        "admitted_velocity_agent2_pr": ADMITTED_VELOCITY_AGENT2_PR,
        "admitted_velocity_agent2_head": ADMITTED_VELOCITY_AGENT2_HEAD,
        "parent_agent4_pr": PARENT_AGENT4_PR,
        "parent_agent4_head": PARENT_AGENT4_HEAD,
        "seed": SEED,
        "local_point_count": int(points.shape[0]),
        "integral_point_count": int(integral_points.shape[0]),
        "fd4_steps": list(FD4_STEPS),
        "quadrature_orders": list(QUADRATURE_ORDERS),
        "time_integral_interval": list(TIME_INTEGRAL_INTERVAL),
        "guards_frozen_before_actions": GUARDS,
        "public_derivative_rms": exact_rms,
        "public_derivative_sampled_max": exact_max,
        "fd4_comparison": fd_rows,
        "fd4_refinement_ratios": fd_refinement,
        "fundamental_theorem": {
            "endpoint_delta_rms": delta_rms,
            "endpoint_delta_sampled_max": delta_max,
            "quadrature_comparison": integral_rows,
            "sign_flip_mutation_relative_rms": sign_flip_mismatch,
        },
        "support_exterior_absolute_max": support_max,
        "failed_guards": failed,
        "public_time_derivative_preflight_passed": passed,
        "independence": {
            "public_values_only": True,
            "agent2_signed_amplitudes_read": False,
            "agent2_target_derivative_helper_read": False,
            "agent2_fd6_verifier_reused": False,
            "independent_local_operator": "centered_fd4",
            "independent_global_operator": "gauss_legendre_fundamental_theorem",
        },
        "truth_boundary": {
            "time_derivative_interface_only": True,
            "velocity_candidate_changed": False,
            "pressure_assessed": False,
            "forcing_assessed": False,
            "heldout_ns_momentum_residual_assessed": False,
            "formal_full_domain_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "final_momentum_gate": 1.0e-3,
            "final_divergence_gate": 1.0e-5,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = run_audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
