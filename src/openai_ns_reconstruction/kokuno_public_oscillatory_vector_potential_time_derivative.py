"""Exact time derivative of the frozen public Kokuno A2 vector potential.

This module is a non-retuning extension of Agent-2 #616.  The admitted
oscillatory velocity remains the Agent-2 #561 public-z pullback field and its
public time derivative remains Agent-2 #579.  The public gauge potential from
#616 is linear in the same two signed covariance amplitudes and, for the frozen
repository candidate, its phase, support, complete-curl coefficient and band
scales are time independent.  Consequently

    d_t A_sigma = (d_t a_sigma / a_sigma) A_sigma.

The signed amplitude law is the repository-autonomous candidate modulation
already documented by #579.  This derivative is exact for that frozen
candidate; it is not a recovered Kokuno time law and is not paper-exact.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_oscillatory_time_derivative import (
    target_log_amplitude_derivative,
    velocity_osc_dt,
)
from .kokuno_public_oscillatory_vector_potential import (
    evaluate_vector_potential_osc,
    vector_potential_osc,
)

TASK = "KOKUNO-A2-PUBLIC-VECTOR-POTENTIAL-DT-043"
SCHEMA = "kokuno-a2-public-oscillatory-vector-potential-dt-v1"
PARENT_AGENT2_PR = 616
PARENT_AGENT2_HEAD = "aa79cb5eb601ad179c79fb8a4f3e21a093d94c31"
ADMITTED_AGENT2_PR = 561
ADMITTED_AGENT2_HEAD = "732800ce4990464b49c8aa32d0dff4580f6684d4"
PUBLIC_TIME_DERIVATIVE_PR = 579
PUBLIC_TIME_DERIVATIVE_HEAD = "6c8e71c800a17c0e9df1802042f9feebf9dbce04"
TIME_FD6_STEPS = (0.01, 0.005, 0.0025)
SPACE_FD6_STEPS = (0.018, 0.009, 0.0045)

VectorEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def _finite(value: Any, name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} must contain only finite real values")
    return out


def evaluate_vector_potential_osc_dt(x: Any, y: Any, z: Any, t: Any) -> dict[str, Any]:
    """Evaluate the exact candidate ``partial_t A_osc`` by beta/sign and total."""
    evaluated = evaluate_vector_potential_osc(x, y, z, t)
    by_sign = np.asarray(
        evaluated["vector_potential_cartesian_by_beta_sign"], dtype=float
    )
    sample_shape = by_sign.shape[:-3]
    tt = np.broadcast_to(_finite(t, "t"), sample_shape)
    log_derivative = target_log_amplitude_derivative(tt)

    dt_by_sign = by_sign * log_derivative[..., None, :, None]
    dt_by_beta = np.sum(dt_by_sign, axis=-2)
    dt_total = np.sum(dt_by_beta, axis=-2)

    return {
        "vector_potential_dt_cartesian_total": dt_total,
        "vector_potential_dt_cartesian_by_beta": dt_by_beta,
        "vector_potential_dt_cartesian_by_beta_sign": dt_by_sign,
        "amplitude_log_derivative_by_sign": log_derivative,
        "support_mask": np.asarray(evaluated["support_mask"], dtype=bool),
        "beta_labels": evaluated["beta_labels"],
        "sign_labels": evaluated["sign_labels"],
        "coordinate_contract": evaluated["coordinate_contract"],
        "derivative_contract": (
            "exact derivative of the frozen repository-autonomous signed-amplitude "
            "modulation; phase/support/gauge complete-curl basis are time independent"
        ),
        "velocity_candidate_changed": False,
        "vector_potential_candidate_changed": False,
        "source_formula_changed": False,
        "paper_exact": False,
        "pde_validated": False,
    }


def vector_potential_osc_dt(x: Any, y: Any, z: Any, t: Any) -> np.ndarray:
    """Return ``partial_t A_osc`` as a Cartesian ``[...,3]`` array."""
    return evaluate_vector_potential_osc_dt(x, y, z, t)[
        "vector_potential_dt_cartesian_total"
    ]


def _fd6_time_derivative(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    return (
        -np.asarray(evaluator(x, y, z, t - 3.0 * step), dtype=float)
        + 9.0 * np.asarray(evaluator(x, y, z, t - 2.0 * step), dtype=float)
        - 45.0 * np.asarray(evaluator(x, y, z, t - step), dtype=float)
        + 45.0 * np.asarray(evaluator(x, y, z, t + step), dtype=float)
        - 9.0 * np.asarray(evaluator(x, y, z, t + 2.0 * step), dtype=float)
        + np.asarray(evaluator(x, y, z, t + 3.0 * step), dtype=float)
    ) / (60.0 * step)


def _fd6_axis_derivative(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
    axis: int,
) -> np.ndarray:
    if not math.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]

    def shifted(multiplier: float) -> np.ndarray:
        args = [value.copy() for value in coords]
        args[axis] = args[axis] + multiplier * step
        return np.asarray(evaluator(args[0], args[1], args[2], t), dtype=float)

    return (
        -shifted(-3.0)
        + 9.0 * shifted(-2.0)
        - 45.0 * shifted(-1.0)
        + 45.0 * shifted(1.0)
        - 9.0 * shifted(2.0)
        + shifted(3.0)
    ) / (60.0 * step)


def _fd6_curl(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    d_dx = _fd6_axis_derivative(evaluator, x, y, z, t, step, 0)
    d_dy = _fd6_axis_derivative(evaluator, x, y, z, t, step, 1)
    d_dz = _fd6_axis_derivative(evaluator, x, y, z, t, step, 2)
    return np.stack(
        (
            d_dy[..., 2] - d_dz[..., 1],
            d_dz[..., 0] - d_dx[..., 2],
            d_dx[..., 1] - d_dy[..., 0],
        ),
        axis=-1,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    radii = (0.50, 0.77, 1.04)
    z_values = (-0.55, 0.61)
    times = (0.34, 0.50, 0.66)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    index = 0
    for time in times:
        for z in z_values:
            for radius in radii:
                theta = 0.173 + (index + 1) * golden
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


def _error_row(approx: np.ndarray, exact: np.ndarray, step: float) -> dict[str, float]:
    error = np.asarray(approx, dtype=float) - np.asarray(exact, dtype=float)
    exact_rms = _vector_rms(exact)
    exact_max = float(np.max(np.linalg.norm(exact, axis=-1)))
    abs_rms = _vector_rms(error)
    abs_max = float(np.max(np.linalg.norm(error, axis=-1)))
    return {
        "step": float(step),
        "absolute_rms_error": abs_rms,
        "absolute_max_error": abs_max,
        "relative_rms_error": abs_rms / max(exact_rms, 1.0e-300),
        "relative_max_error": abs_max / max(exact_max, 1.0e-300),
    }


def _ratios(rows: list[dict[str, float]]) -> list[float]:
    return [
        rows[i]["relative_rms_error"] / max(rows[i + 1]["relative_rms_error"], 1.0e-300)
        for i in range(len(rows) - 1)
    ]


def verification_receipt() -> dict[str, Any]:
    """Cross-check exact ``A_t`` in time and through ``curl(A_t)=u_t``."""
    x, y, z, t = _verification_cloud()
    exact_at = np.asarray(vector_potential_osc_dt(x, y, z, t), dtype=float)
    exact_ut = np.asarray(velocity_osc_dt(x, y, z, t), dtype=float)

    time_rows = [
        _error_row(
            _fd6_time_derivative(vector_potential_osc, x, y, z, t, step),
            exact_at,
            step,
        )
        for step in TIME_FD6_STEPS
    ]
    curl_rows = [
        _error_row(_fd6_curl(vector_potential_osc_dt, x, y, z, t, step), exact_ut, step)
        for step in SPACE_FD6_STEPS
    ]
    time_ratios = _ratios(time_rows)
    curl_ratios = _ratios(curl_rows)

    exterior = np.asarray(
        vector_potential_osc_dt(
            np.asarray((0.0, 1.8, 0.8)),
            np.asarray((0.0, 0.0, 0.0)),
            np.asarray((0.0, 0.0, 2.1)),
            np.asarray((0.5, 0.5, 0.5)),
        ),
        dtype=float,
    )
    exterior_max = float(np.max(np.abs(exterior)))

    guards = {
        "time_finest_relative_rms": time_rows[-1]["relative_rms_error"] <= 5.0e-8,
        "time_finest_relative_max": time_rows[-1]["relative_max_error"] <= 1.0e-7,
        "time_refinement": min(time_ratios) >= 20.0,
        "curl_finest_relative_rms": curl_rows[-1]["relative_rms_error"] <= 5.0e-6,
        "curl_finest_relative_max": curl_rows[-1]["relative_max_error"] <= 1.0e-5,
        "curl_refinement": min(curl_ratios) >= 20.0,
        "nontrivial_vector_potential_dt": _vector_rms(exact_at) >= 1.0e-8,
        "nontrivial_velocity_dt": _vector_rms(exact_ut) >= 1.0e-8,
        "support_exterior": exterior_max <= 1.0e-12,
    }

    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "admitted_agent2_pr": ADMITTED_AGENT2_PR,
        "admitted_agent2_head": ADMITTED_AGENT2_HEAD,
        "public_time_derivative_pr": PUBLIC_TIME_DERIVATIVE_PR,
        "public_time_derivative_head": PUBLIC_TIME_DERIVATIVE_HEAD,
        "public_vector_potential_provider": (
            "openai_ns_reconstruction.kokuno_public_oscillatory_vector_potential:vector_potential_osc"
        ),
        "public_vector_potential_dt_provider": (
            "openai_ns_reconstruction.kokuno_public_oscillatory_vector_potential_time_derivative:vector_potential_osc_dt"
        ),
        "public_velocity_dt_provider": (
            "openai_ns_reconstruction.kokuno_public_oscillatory_time_derivative:velocity_osc_dt"
        ),
        "sample_count": int(t.size),
        "verification_times": sorted({float(v) for v in t}),
        "time_fd6_steps": list(TIME_FD6_STEPS),
        "space_fd6_steps": list(SPACE_FD6_STEPS),
        "vector_potential_dt_rms": _vector_rms(exact_at),
        "velocity_dt_rms": _vector_rms(exact_ut),
        "time_fd6_comparison": time_rows,
        "time_rms_refinement_ratios": time_ratios,
        "curl_dt_fd6_comparison": curl_rows,
        "curl_dt_rms_refinement_ratios": curl_ratios,
        "support_exterior_absolute_max": exterior_max,
        "guards": guards,
        "failed_guards": [name for name, passed in guards.items() if not passed],
        "derivation": (
            "the frozen public gauge potential and admitted complete-curl velocity are linear "
            "in the same signed amplitudes; therefore dt A_sigma=(dt a_sigma/a_sigma) A_sigma "
            "and continuum differentiation commutes as curl(dt A)=dt curl(A)=dt u"
        ),
        "provenance": {
            "source_localized_vector_potential_formula_reused": True,
            "signed_amplitude_inverse_formula_reused": True,
            "time_modulation": "repository_autonomous_candidate_choice",
            "public_gauge_Q_scaling": "repository_autonomous_public_coordinate_realization",
        },
        "truth_boundary": {
            "velocity_candidate_changed": False,
            "vector_potential_candidate_changed": False,
            "source_formula_changed": False,
            "public_vector_potential_time_derivative_materialized": True,
            "agent2_self_check_only": True,
            "independent_agent4_vector_potential_audit_required": True,
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
