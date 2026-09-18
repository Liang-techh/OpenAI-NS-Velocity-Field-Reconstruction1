"""Map public oscillatory velocity tangents into Kokuno mean-covariance columns.

Agent 3 owns the zero-harmonic mean/stress correction layer. Agent 2 owns the
projected-pulse and complete-curl construction. This module joins those layers
without reimplementing Agent 2 mathematics: a caller supplies phase-dependent
public velocity values ``w`` and a phase-dependent public tangent
``s = partial_xi w``.

For either cylindrical covariance component the product rule gives

    partial_xi <w_r w_theta> = <s_r w_theta + w_r s_theta>,
    partial_xi <w_r w_z>     = <s_r w_z     + w_r s_z>.

The resulting two-component response can be sent directly to Agent 3's bounded
two-column inverse preflight. It is not itself a velocity correction and does
not authorize another finite correction cycle.

The real calibration in this module uses only the already-public Agent-2
complete-curl amplitude tangent. Because that tangent stays in the existing
column, it must reproduce the current rank-one covariance derivative and be
rejected as a fake second direction. A future actual D_r/D_z public tangent can
use the same interface once Agent 2 materializes one.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from math import pi
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    PROFILE_Z,
    MissingCovarianceColumnTarget,
    build_missing_covariance_column_target,
    generate_actual_core_report as generate_missing_target_report,
)
from .kokuno_second_column_bounded_inverse import (
    current_signed_coefficient_budget,
    evaluate_bounded_second_column,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)

TASK = "KOKUNO-A3-COVARIANCE-TANGENT-PREFLIGHT-011"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "Signed covariance and exact signed differential inverse"

PhaseVelocityEvaluator = Callable[[np.ndarray, float, float], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("values must have shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _wrap_phase(value: float) -> float:
    wrapped = (float(value) + pi) % (2.0 * pi) - pi
    if wrapped == -pi and value > 0.0:
        return pi
    return wrapped


def _ring_points(radii: np.ndarray, *, angular_count: int, z: float) -> np.ndarray:
    radii = np.asarray(radii, dtype=float)
    if radii.ndim != 1 or len(radii) < 3 or not np.isfinite(radii).all():
        raise ValueError("radii must be a finite one-dimensional array")
    if np.any(radii <= 0.0) or np.any(np.diff(radii) <= 0.0):
        raise ValueError("radii must be positive and strictly increasing")
    if angular_count < 8 or angular_count % 2:
        raise ValueError("angular_count must be an even integer >= 8")
    if not np.isfinite(z):
        raise ValueError("z must be finite")
    angles = 2.0 * pi * np.arange(angular_count, dtype=float) / angular_count
    rr, aa = np.meshgrid(radii, angles, indexing="ij")
    return np.stack(
        (rr * np.cos(aa), rr * np.sin(aa), np.full_like(rr, float(z))),
        axis=-1,
    )


def _cylindrical_velocity(
    values: np.ndarray, points: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.asarray(values, dtype=float)
    points = np.asarray(points, dtype=float)
    if values.shape != points.shape or values.shape[-1] != 3:
        raise ValueError("velocity and points must have matching (...,3) shape")
    if not np.isfinite(values).all() or not np.isfinite(points).all():
        raise ValueError("velocity and points must be finite")
    radius = np.hypot(points[..., 0], points[..., 1])
    if np.any(radius <= 1.0e-12):
        raise ValueError("cylindrical covariance is undefined on the axis")
    ur = (points[..., 0] * values[..., 0] + points[..., 1] * values[..., 1]) / radius
    utheta = (
        -points[..., 1] * values[..., 0] + points[..., 0] * values[..., 1]
    ) / radius
    return ur, utheta, values[..., 2]


def complete_curl_phase_evaluator(
    correction: KokunoCompleteCurlCorrection,
) -> PhaseVelocityEvaluator:
    """Adapt Agent 2's public complete-curl object to a phase evaluator.

    ``phase_offset`` is added to the object's frozen phase. No amplitude,
    polarization, support, or wave-vector parameter is changed here.
    """
    if not isinstance(correction, KokunoCompleteCurlCorrection):
        raise TypeError("correction must be KokunoCompleteCurlCorrection")

    def evaluate(points: np.ndarray, time: float, phase_offset: float) -> np.ndarray:
        if not np.isfinite([time, phase_offset]).all():
            raise ValueError("time and phase_offset must be finite")
        probe = replace(
            correction,
            phase=_wrap_phase(float(correction.phase) + float(phase_offset)),
        )
        values = np.asarray(probe.at_points(points, float(time)), dtype=float)
        expected = np.asarray(points).shape
        if values.shape != expected or expected[-1:] != (3,):
            raise ValueError("public complete-curl evaluator returned the wrong shape")
        if not np.isfinite(values).all():
            raise ValueError("public complete-curl evaluator returned nonfinite values")
        return values

    return evaluate


def measure_phase_mean_covariance_tangent(
    base_velocity: PhaseVelocityEvaluator,
    tangent_velocity: PhaseVelocityEvaluator,
    radii: np.ndarray,
    *,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> np.ndarray:
    """Measure the derivative of ``(<w_r w_theta>, <w_r w_z>)``.

    Both the base field and its tangent are evaluated through caller-supplied
    public velocity interfaces on the same physical-angle and phase quadrature.
    The tangent is not inferred by finite differences inside this function.
    """
    if not callable(base_velocity) or not callable(tangent_velocity):
        raise TypeError("base_velocity and tangent_velocity must be callable")
    if phase_count < 4 or phase_count % 2:
        raise ValueError("phase_count must be an even integer >= 4")
    if not np.isfinite([time, z]).all():
        raise ValueError("time and z must be finite")

    points = _ring_points(radii, angular_count=angular_count, z=z)
    flat = points.reshape(-1, 3)
    accumulator = np.zeros((len(radii), 2), dtype=float)
    shifts = 2.0 * pi * np.arange(phase_count, dtype=float) / phase_count

    for shift in shifts:
        base = np.asarray(base_velocity(flat, float(time), float(shift)), dtype=float)
        tangent = np.asarray(
            tangent_velocity(flat, float(time), float(shift)), dtype=float
        )
        if base.shape != flat.shape or tangent.shape != flat.shape:
            raise ValueError("phase velocity evaluators must return shape (N,3)")
        if not np.isfinite(base).all() or not np.isfinite(tangent).all():
            raise ValueError("phase velocity evaluators must return finite values")

        base = base.reshape(points.shape)
        tangent = tangent.reshape(points.shape)
        wr, wtheta, wz = _cylindrical_velocity(base, points)
        sr, stheta, sz = _cylindrical_velocity(tangent, points)
        accumulator[:, 0] += np.mean(sr * wtheta + wr * stheta, axis=1)
        accumulator[:, 1] += np.mean(sr * wz + wr * sz, axis=1)

    return accumulator / float(phase_count)


def evaluate_covariance_tangent_preflight(
    receipt: MissingCovarianceColumnTarget,
    second_response: np.ndarray,
    *,
    coefficient_budget: float | None = None,
) -> dict[str, Any]:
    """Route one measured public covariance tangent through the #302 guard."""
    if not isinstance(receipt, MissingCovarianceColumnTarget):
        raise TypeError("receipt must be MissingCovarianceColumnTarget")
    response = np.asarray(second_response, dtype=float)
    if response.shape != receipt.target_stress.shape or not np.isfinite(response).all():
        raise ValueError("second_response must be finite and match target shape")
    budget = (
        current_signed_coefficient_budget(receipt)
        if coefficient_budget is None
        else float(coefficient_budget)
    )
    bounded = evaluate_bounded_second_column(
        receipt,
        response,
        coefficient_budget=budget,
    )
    return {
        "coefficient_budget": budget,
        "response_vector_rms": _vector_rms(response[receipt.active_target_mask]),
        "bounded_inverse": bounded,
        "finite_cycle_rerun_allowed": bool(bounded["bounded_inverse_preflight_passed"]),
        "interpretation": (
            "A measured public velocity tangent may enter a two-column finite cycle only "
            "after this covariance/stress inverse preflight passes. Passing is necessary, "
            "not sufficient, for Navier-Stokes residual reduction."
        ),
    }


def _receipt_from_target_report(report: dict[str, Any]) -> MissingCovarianceColumnTarget:
    target = report["missing_second_column_target"]
    inputs = report["inputs"]
    return build_missing_covariance_column_target(
        np.asarray(target["radii"], dtype=float),
        np.asarray(target["target_stress_theta_axial"], dtype=float),
        np.asarray(target["current_response_theta_axial"], dtype=float),
        target_floor_fraction=float(inputs["target_floor_fraction"]),
        response_floor_fraction=float(inputs["response_floor_fraction"]),
    )


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/covariance_tangent_preflight_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/covariance_tangent_preflight_target.json",
) -> dict[str, Any]:
    """Calibrate the tangent bridge on the existing physical rank-one column."""
    target_report = generate_missing_target_report(output=target_output)
    receipt = _receipt_from_target_report(target_report)

    base = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    amplitude_tangent = replace(base, amplitude=1.0)
    measured = measure_phase_mean_covariance_tangent(
        complete_curl_phase_evaluator(base),
        complete_curl_phase_evaluator(amplitude_tangent),
        receipt.radii,
        time=PROFILE_TIME,
        z=PROFILE_Z,
        angular_count=int(target_report["inputs"]["angular_count"]),
        phase_count=int(target_report["inputs"]["phase_count"]),
    )
    active = receipt.active_target_mask
    calibration_relative_error = _vector_rms(
        measured[active] - receipt.current_response[active]
    ) / max(_vector_rms(receipt.current_response[active]), np.finfo(float).tiny)

    preflight = evaluate_covariance_tangent_preflight(receipt, measured)
    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "section": SOURCE_SECTION,
            "source_structure": (
                "Kokuno uses a two-column signed covariance inverse. The product-rule "
                "covariance tangent implemented here is repository algebra for consuming "
                "a public velocity tangent; it is not quoted as a new source formula."
            ),
        },
        "inputs": {
            "leading_candidate_sha256": target_report["inputs"]["leading_candidate_sha256"],
            "oscillatory_amplitude": float(base.amplitude),
            "oscillatory_phase": float(base.phase),
            "profile_time": PROFILE_TIME,
            "profile_z": PROFILE_Z,
            "profile_annulus": list(PROFILE_ANNULUS),
            "angular_count": int(target_report["inputs"]["angular_count"]),
            "phase_count": int(target_report["inputs"]["phase_count"]),
            "surrogate_defect_used": False,
        },
        "amplitude_tangent_calibration": {
            "response_vector_rms": _vector_rms(measured[active]),
            "current_response_vector_rms": _vector_rms(
                receipt.current_response[active]
            ),
            "relative_error_to_existing_covariance_derivative": calibration_relative_error,
            "bounded_preflight": preflight,
            "expected_routing": (
                "reject: an amplitude tangent is the existing covariance direction, "
                "not an independent second column"
            ),
        },
        "agent2_handoff": {
            "projected_pulse_parameter_sensitivity_api_exists": True,
            "actual_source_D_r_public_velocity_tangent_available": False,
            "actual_source_D_z_public_velocity_tangent_available": False,
            "second_public_covariance_column_available": False,
            "next_admissible_input": (
                "Agent-2 public phase-dependent [u,v,w] tangent from an actual source "
                "D_r or D_z path after the source forcing/background seam is instantiated"
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_target_consumed": True,
            "surrogate_defect_used": False,
            "covariance_tangent_bridge_executable": True,
            "amplitude_tangent_used_as_second_physical_column": False,
            "new_oscillatory_column_constructed": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "finite_cycle_rerun_allowed": bool(
                preflight["finite_cycle_rerun_allowed"]
            ),
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Preflight a public Kokuno velocity tangent as a second covariance column"
    )
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/covariance_tangent_preflight_report.json",
    )
    parser.add_argument(
        "--target-output",
        default="artifacts/kokuno_agent3/covariance_tangent_preflight_target.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(
        output=args.output,
        target_output=args.target_output,
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
