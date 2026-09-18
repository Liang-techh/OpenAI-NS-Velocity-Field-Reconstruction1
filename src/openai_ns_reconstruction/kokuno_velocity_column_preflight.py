"""Preflight an independently materialized oscillatory velocity column for mean correction.

Agent 3 previously required a caller-supplied public velocity tangent in order to
form a two-component mean-covariance response.  Agent 2 now has a source-real
m=+/-1 complete-curl pair materializer (PR #321) whose natural handoff is a
physical velocity column rather than an already differentiated amplitude
interface.

For any phase-dependent unit column ``w`` and a positive multiplicative
amplitude ``a`` the repository algebra is

    d/da <(a w)_r (a w)_theta> = 2 a <w_r w_theta>,
    d/da <(a w)_r (a w)_z>     = 2 a <w_r w_z>.

This module turns that identity into an adapter for Agent 3's existing bounded
two-column inverse preflight.  It does not reconstruct Agent 2's pulse/curl
machinery and does not claim the displayed derivative as an additional Kokuno
source formula.  A source-real pair still needs actual source mode data and a
phase/physical-coordinate evaluator before it can be scored here.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_covariance_tangent_preflight import (
    complete_curl_phase_evaluator,
    evaluate_covariance_tangent_preflight,
    measure_phase_mean_covariance_tangent,
)
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    PROFILE_Z,
    MissingCovarianceColumnTarget,
    build_missing_covariance_column_target,
    generate_actual_core_report as generate_missing_target_report,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)

TASK = "KOKUNO-A3-SOURCE-REAL-PAIR-COLUMN-PREFLIGHT-012"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_SECTION = "Signed covariance / real conjugate oscillatory pair"
AGENT2_REAL_PAIR_PR = 321
AGENT2_REAL_PAIR_HEAD = "1c5cf93b7468441a0f735336837a4ecdcb8e5afd"
AGENT2_REAL_PAIR_DEDICATED_RUN = 35307264120

PhaseVelocityEvaluator = Callable[[np.ndarray, float, float], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("values must have shape (N,2)")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def scaled_phase_velocity(
    unit_velocity: PhaseVelocityEvaluator,
    amplitude: float,
) -> PhaseVelocityEvaluator:
    """Return ``amplitude * unit_velocity`` with fail-closed output checks."""
    if not callable(unit_velocity):
        raise TypeError("unit_velocity must be callable")
    amplitude = float(amplitude)
    if not np.isfinite(amplitude) or amplitude <= 0.0:
        raise ValueError("amplitude must be finite and positive")

    def evaluate(points: np.ndarray, time: float, phase_offset: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must be a finite array with shape (N,3)")
        if not np.isfinite([time, phase_offset]).all():
            raise ValueError("time and phase_offset must be finite")
        values = np.asarray(unit_velocity(points, float(time), float(phase_offset)), dtype=float)
        if values.shape != points.shape or not np.isfinite(values).all():
            raise ValueError("unit_velocity must return finite shape-(N,3) values")
        return amplitude * values

    return evaluate


def measure_multiplicative_covariance_response(
    unit_velocity: PhaseVelocityEvaluator,
    radii: np.ndarray,
    *,
    amplitude: float = 1.0,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> np.ndarray:
    """Measure the covariance derivative of one multiplicative velocity column.

    The unit column is treated as an independently materialized velocity basis.
    The amplitude derivative is evaluated by the exact product-rule bridge
    already used by Agent 3; no finite-difference amplitude perturbation is
    introduced here.
    """
    amplitude = float(amplitude)
    base = scaled_phase_velocity(unit_velocity, amplitude)
    return measure_phase_mean_covariance_tangent(
        base,
        unit_velocity,
        np.asarray(radii, dtype=float),
        time=float(time),
        z=float(z),
        angular_count=int(angular_count),
        phase_count=int(phase_count),
    )


def evaluate_velocity_column_preflight(
    receipt: MissingCovarianceColumnTarget,
    unit_velocity: PhaseVelocityEvaluator,
    *,
    amplitude: float = 1.0,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
    coefficient_budget: float | None = None,
) -> dict[str, Any]:
    """Measure one public velocity column and route it through the bounded inverse."""
    if not isinstance(receipt, MissingCovarianceColumnTarget):
        raise TypeError("receipt must be MissingCovarianceColumnTarget")
    response = measure_multiplicative_covariance_response(
        unit_velocity,
        receipt.radii,
        amplitude=amplitude,
        time=time,
        z=z,
        angular_count=angular_count,
        phase_count=phase_count,
    )
    routed = evaluate_covariance_tangent_preflight(
        receipt,
        response,
        coefficient_budget=coefficient_budget,
    )
    return {
        "amplitude": float(amplitude),
        "response_vector_rms": _vector_rms(response[receipt.active_target_mask]),
        "response_theta_axial": response.tolist(),
        "bounded_preflight": routed,
        "finite_cycle_rerun_allowed": bool(routed["finite_cycle_rerun_allowed"]),
        "interpretation": (
            "An independently materialized velocity column can be screened without an "
            "explicit amplitude-tangent API. Passing remains only a prerequisite for "
            "a later two-column velocity correction and held-out NS residual cycle."
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
    output: str | Path = "artifacts/kokuno_agent3/velocity_column_preflight_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/velocity_column_preflight_target.json",
) -> dict[str, Any]:
    """Calibrate the new column adapter on the existing real rank-one response."""
    target_report = generate_missing_target_report(output=target_output)
    receipt = _receipt_from_target_report(target_report)

    unit = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=1.0,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    measured = measure_multiplicative_covariance_response(
        complete_curl_phase_evaluator(unit),
        receipt.radii,
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        time=PROFILE_TIME,
        z=PROFILE_Z,
        angular_count=int(target_report["inputs"]["angular_count"]),
        phase_count=int(target_report["inputs"]["phase_count"]),
    )
    active = receipt.active_target_mask
    reference_scale = max(
        _vector_rms(receipt.current_response[active]),
        np.finfo(float).tiny,
    )
    calibration_relative_error = (
        _vector_rms(measured[active] - receipt.current_response[active])
        / reference_scale
    )
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
                "Kokuno uses real conjugate oscillatory pairs and a two-column signed "
                "covariance inverse. The multiplicative covariance derivative used by "
                "this adapter is repository algebra, not an additional quoted source formula."
            ),
        },
        "inputs": {
            "leading_candidate_sha256": target_report["inputs"]["leading_candidate_sha256"],
            "oscillatory_amplitude": ROUTED_OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": ROUTED_OSCILLATORY_PHASE,
            "profile_time": PROFILE_TIME,
            "profile_z": PROFILE_Z,
            "profile_annulus": list(PROFILE_ANNULUS),
            "angular_count": int(target_report["inputs"]["angular_count"]),
            "phase_count": int(target_report["inputs"]["phase_count"]),
            "surrogate_defect_used": False,
        },
        "multiplicative_column_calibration": {
            "response_vector_rms": _vector_rms(measured[active]),
            "current_response_vector_rms": _vector_rms(receipt.current_response[active]),
            "relative_error_to_existing_covariance_derivative": calibration_relative_error,
            "bounded_preflight": preflight,
            "expected_routing": (
                "reject as a second direction: this calibration intentionally reuses the "
                "existing physical column and only validates the new velocity-column adapter"
            ),
        },
        "agent2_real_pair_handoff": {
            "pr": AGENT2_REAL_PAIR_PR,
            "exact_head": AGENT2_REAL_PAIR_HEAD,
            "dedicated_ci_run": AGENT2_REAL_PAIR_DEDICATED_RUN,
            "source_real_conjugate_pair_materializer_available": True,
            "physical_cartesian_pointwise_velocity_available_from_supplied_mode_data": True,
            "explicit_amplitude_tangent_required_by_agent3": False,
            "source_actual_pulse_forcing_instantiated": False,
            "source_actual_background_path_instantiated": False,
            "public_xyz_t_velocity_correction_materialized": False,
            "genuinely_independent_second_source_pair_available": False,
            "next_admissible_input": (
                "a source-motivated second real-pair phase evaluator built from actual mode "
                "data; Agent 3 can now measure its covariance column directly before a finite cycle"
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_target_consumed": True,
            "surrogate_defect_used": False,
            "velocity_column_covariance_adapter_executable": True,
            "agent2_real_pair_code_reimplemented": False,
            "new_oscillatory_column_constructed": False,
            "second_public_covariance_column_available": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "finite_cycle_rerun_allowed": bool(preflight["finite_cycle_rerun_allowed"]),
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
        description="Preflight an independently materialized Kokuno velocity covariance column"
    )
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/velocity_column_preflight_report.json",
    )
    parser.add_argument(
        "--target-output",
        default="artifacts/kokuno_agent3/velocity_column_preflight_target.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output, target_output=args.target_output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
