"""Public potential-time seam for the frozen A2 signed-amplitude correction.

This Agent-2 module exposes the time derivative of the already-materialized
correction vector potential and checks the two commutation identities needed by
an independent complete-curl audit:

    d_t(delta A) = delta A_t,
    curl(delta A_t) = delta u_t.

The underlying correction remains vector-potential first. No Agent-3 defect,
mean/radial inverse, pressure, forcing, target, gain, or residual enters this
API. The corrected Kokuno reconstruction supplies structural provenance only;
the public-z pullback, compact spline lift, autonomous time modulation, and
finite-difference audit are repository realizations rather than paper-exact
hidden data.
"""
from __future__ import annotations

import argparse
import inspect
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_public_signed_amplitude_complete_curl import (
    SignedRadialAmplitudeCorrection,
    correction_vector_potential,
    correction_velocity_dt,
    evaluate_signed_amplitude_correction,
    profile_from_delta_a,
)

TASK = "KOKUNO-A2-SIGNED-AMPLITUDE-POTENTIAL-TIME-SEAM-051"
SCHEMA = "kokuno-a2-signed-amplitude-potential-time-seam-v1"
PARENT_AGENT2_PR = 734
PARENT_AGENT2_HEAD = "c535eeec3267869f630c1327506b3683d0417969"
SOURCE_CORRECTED_READER_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_CORRECTED_READER_DATE = "2026-09-09"
TIME_FD6_STEPS = (0.012, 0.006, 0.003)
CURL_FD4_STEPS = (0.012, 0.006, 0.003)

VectorEvaluator = Callable[[Any, Any, Any, Any], np.ndarray]


def correction_vector_potential_dt(
    correction: SignedRadialAmplitudeCorrection,
    x: Any,
    y: Any,
    z: Any,
    t: Any,
) -> np.ndarray:
    """Return the frozen-candidate time derivative of ``delta A``.

    Time dependence is exactly the repository-autonomous signed-amplitude lift
    already used by the parent correction kernel. This function exposes that
    existing quantity; it does not introduce a new time law.
    """
    return np.asarray(
        evaluate_signed_amplitude_correction(correction, x, y, z, t)[
            "vector_potential_dt_cartesian_total"
        ],
        dtype=float,
    )


def _finite(value: Any, label: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{label} must contain only finite real values")
    return out


def _fd4_axis_derivative(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
    axis: int,
) -> np.ndarray:
    h = float(step)
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError("step must be positive and finite")
    coords = [np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float)]

    def shifted(multiplier: float) -> np.ndarray:
        args = [coord.copy() for coord in coords]
        args[axis] = args[axis] + multiplier * h
        return np.asarray(evaluator(args[0], args[1], args[2], t), dtype=float)

    return (
        shifted(-2.0)
        - 8.0 * shifted(-1.0)
        + 8.0 * shifted(1.0)
        - shifted(2.0)
    ) / (12.0 * h)


def _fd4_curl(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    d_dx = _fd4_axis_derivative(evaluator, x, y, z, t, step, 0)
    d_dy = _fd4_axis_derivative(evaluator, x, y, z, t, step, 1)
    d_dz = _fd4_axis_derivative(evaluator, x, y, z, t, step, 2)
    return np.stack(
        (
            d_dy[..., 2] - d_dz[..., 1],
            d_dz[..., 0] - d_dx[..., 2],
            d_dx[..., 1] - d_dy[..., 0],
        ),
        axis=-1,
    )


def _fd6_time_derivative(
    evaluator: VectorEvaluator,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    t: np.ndarray,
    step: float,
) -> np.ndarray:
    h = float(step)
    if not math.isfinite(h) or h <= 0.0:
        raise ValueError("step must be positive and finite")
    return (
        -np.asarray(evaluator(x, y, z, t - 3.0 * h), dtype=float)
        + 9.0 * np.asarray(evaluator(x, y, z, t - 2.0 * h), dtype=float)
        - 45.0 * np.asarray(evaluator(x, y, z, t - h), dtype=float)
        + 45.0 * np.asarray(evaluator(x, y, z, t + h), dtype=float)
        - 9.0 * np.asarray(evaluator(x, y, z, t + 2.0 * h), dtype=float)
        + np.asarray(evaluator(x, y, z, t + 3.0 * h), dtype=float)
    ) / (60.0 * h)


def _vector_rms(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(arr * arr, axis=-1))))


def _verification_profile() -> SignedRadialAmplitudeCorrection:
    radii = np.linspace(0.29, 1.21, 19)
    s = (radii - radii[0]) / (radii[-1] - radii[0])
    bump = np.sin(math.pi * s) ** 6
    delta_a = np.stack(
        (
            0.011 * bump * (1.0 + 0.11 * np.cos(2.0 * math.pi * s)),
            -0.0075 * bump * (1.0 - 0.08 * np.sin(2.0 * math.pi * s)),
        ),
        axis=-1,
    )
    delta_a[[0, -1], :] = 0.0
    return profile_from_delta_a(
        radii,
        delta_a,
        reference_time=0.50,
        producer_kind="analytic-regression-not-agent3-real-candidate",
        provenance=(
            "fresh signed radial mechanics profile for potential-time/curl-time "
            "commutation only; no Agent-3 defect or residual input"
        ),
        source_mean_amplitude_differential_certified=False,
    )


def _verification_cloud() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    radii = (0.43, 0.69, 0.95)
    z_values = (-0.57, 0.61)
    times = (0.38, 0.50, 0.62)
    golden = math.pi * (3.0 - math.sqrt(5.0))
    rows: list[tuple[float, float, float, float]] = []
    idx = 0
    for time in times:
        for z_value in z_values:
            for radius in radii:
                theta = 0.271 + (idx + 1) * golden
                rows.append(
                    (
                        radius * math.cos(theta),
                        radius * math.sin(theta),
                        z_value,
                        time,
                    )
                )
                idx += 1
    arr = np.asarray(rows, dtype=float)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]


def verification_receipt() -> dict[str, Any]:
    correction = _verification_profile()
    x, y, z, t = _verification_cloud()

    potential = lambda xx, yy, zz, tt: correction_vector_potential(
        correction, xx, yy, zz, tt
    )
    potential_dt = lambda xx, yy, zz, tt: correction_vector_potential_dt(
        correction, xx, yy, zz, tt
    )
    velocity_dt = lambda xx, yy, zz, tt: correction_velocity_dt(
        correction, xx, yy, zz, tt
    )

    exact_A_dt = potential_dt(x, y, z, t)
    exact_u_dt = velocity_dt(x, y, z, t)
    A_dt_rms = _vector_rms(exact_A_dt)
    A_dt_max = float(np.max(np.linalg.norm(exact_A_dt, axis=-1)))
    u_dt_rms = _vector_rms(exact_u_dt)
    u_dt_max = float(np.max(np.linalg.norm(exact_u_dt, axis=-1)))

    time_rows: list[dict[str, float]] = []
    for step in TIME_FD6_STEPS:
        numerical = _fd6_time_derivative(potential, x, y, z, t, step)
        error = numerical - exact_A_dt
        abs_rms = _vector_rms(error)
        abs_max = float(np.max(np.linalg.norm(error, axis=-1)))
        time_rows.append(
            {
                "step": step,
                "relative_rms_error": abs_rms / max(A_dt_rms, 1.0e-300),
                "relative_max_error": abs_max / max(A_dt_max, 1.0e-300),
            }
        )
    time_refinement = [
        time_rows[i]["relative_rms_error"]
        / max(time_rows[i + 1]["relative_rms_error"], 1.0e-300)
        for i in range(len(time_rows) - 1)
    ]

    curl_rows: list[dict[str, float]] = []
    for step in CURL_FD4_STEPS:
        numerical = _fd4_curl(potential_dt, x, y, z, t, step)
        error = numerical - exact_u_dt
        abs_rms = _vector_rms(error)
        abs_max = float(np.max(np.linalg.norm(error, axis=-1)))
        curl_rows.append(
            {
                "step": step,
                "relative_rms_error": abs_rms / max(u_dt_rms, 1.0e-300),
                "relative_max_error": abs_max / max(u_dt_max, 1.0e-300),
            }
        )
    curl_refinement = [
        curl_rows[i]["relative_rms_error"]
        / max(curl_rows[i + 1]["relative_rms_error"], 1.0e-300)
        for i in range(len(curl_rows) - 1)
    ]

    exterior_x = np.asarray((0.0, 0.18, 1.33, 0.74), dtype=float)
    exterior_y = np.zeros_like(exterior_x)
    exterior_z = np.asarray((0.0, 0.0, 0.0, 2.05), dtype=float)
    exterior_t = np.full_like(exterior_x, 0.50)
    exterior_A_dt = potential_dt(exterior_x, exterior_y, exterior_z, exterior_t)
    exterior_u_dt = velocity_dt(exterior_x, exterior_y, exterior_z, exterior_t)
    exterior_max = float(
        max(
            np.max(np.abs(exterior_A_dt), initial=0.0),
            np.max(np.abs(exterior_u_dt), initial=0.0),
        )
    )

    guards = {
        "nontrivial_vector_potential_dt": A_dt_rms >= 1.0e-8,
        "nontrivial_velocity_dt": u_dt_rms >= 1.0e-8,
        "time_finest_relative_rms": time_rows[-1]["relative_rms_error"] <= 5.0e-8,
        "time_finest_relative_max": time_rows[-1]["relative_max_error"] <= 1.0e-7,
        "time_refinement": min(time_refinement) >= 20.0,
        "curl_finest_relative_rms": curl_rows[-1]["relative_rms_error"] <= 2.0e-5,
        "curl_finest_relative_max": curl_rows[-1]["relative_max_error"] <= 5.0e-5,
        "curl_refinement": min(curl_refinement) >= 8.0,
        "support_exterior": exterior_max <= 1.0e-12,
    }

    return {
        "task": TASK,
        "schema": SCHEMA,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
        "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        "sample_count": int(t.size),
        "vector_potential_dt_rms": A_dt_rms,
        "velocity_dt_rms": u_dt_rms,
        "time_fd6_steps": list(TIME_FD6_STEPS),
        "potential_time_derivative_comparison": time_rows,
        "potential_time_refinement_ratios": time_refinement,
        "curl_fd4_steps": list(CURL_FD4_STEPS),
        "curl_A_dt_vs_velocity_dt": curl_rows,
        "curl_rms_refinement_ratios": curl_refinement,
        "support_exterior_absolute_max": exterior_max,
        "guards": guards,
        "failed_guards": [name for name, passed in guards.items() if not passed],
        "provenance": {
            "source_structure": (
                "signed complete-curl/vector-potential structure and coefficient-level "
                "localization from the pinned corrected Kokuno reconstruction"
            ),
            "repository_realization": (
                "public-z pullback, compact quintic radial lift, autonomous signed time "
                "modulation, and FD6/FD4 audit operators"
            ),
            "paper_exact_hidden_time_law_claimed": False,
            "agent3_lane_reimplemented": False,
        },
        "truth_boundary": truth_boundary(),
    }


def truth_boundary() -> dict[str, Any]:
    signature = inspect.signature(correction_vector_potential_dt)
    forbidden = {
        "residual",
        "defect",
        "stress",
        "pressure",
        "forcing",
        "target",
        "gain",
        "normalized_score",
        "scientific_threshold",
        "delta_y",
    }
    return {
        "source_corrected_reader_commit": SOURCE_CORRECTED_READER_COMMIT,
        "source_corrected_reader_date": SOURCE_CORRECTED_READER_DATE,
        "parent_agent2_pr": PARENT_AGENT2_PR,
        "parent_agent2_head": PARENT_AGENT2_HEAD,
        "public_vector_potential_time_derivative_executable": True,
        "potential_time_commutation_auditable": True,
        "curl_time_commutation_auditable": True,
        "vector_potential_first_complete_curl_reused": True,
        "forbidden_public_parameters_absent": forbidden.isdisjoint(signature.parameters),
        "candidate_amplitude_phase_support_retuned": False,
        "agent3_mean_radial_chain_reimplemented": False,
        "source_agent2_complete_curl_certified_for_full_candidate": False,
        "independent_agent4_correction_vector_potential_audit_required": True,
        "real_agent3_delta_a_bound": False,
        "real_full_candidate_correction_cycle_run": False,
        "heldout_ns_momentum_residual_assessed": False,
        "residual_reduction_claimed": False,
        "paper_exact": False,
        "pde_validated": False,
        "openai_field_identified": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = verification_receipt()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if receipt["failed_guards"]:
        raise SystemExit(
            "preregistered signed-amplitude potential-time guards failed: "
            + ", ".join(receipt["failed_guards"])
        )


if __name__ == "__main__":
    main()
