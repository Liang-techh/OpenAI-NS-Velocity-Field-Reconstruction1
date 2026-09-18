"""Bridge Agent-2 bounded family coordinates into Agent-3 mean covariance machinery.

Agent-2 PR #400 exposes two repository-autonomous coordinates on an already
physical, Q-scaled slow-label family,

    m_beta = 1 + delta_common + delta_band g_beta,

with dimensionless fractional-multiplier units.  Agent-3 PR #393, however,
freezes the existing correction budget in the older physical-amplitude
coordinate.  Reusing the same numeric budget in the new coordinate would be a
unit error.  This module makes the conversion explicit and measures the full
cross-term-aware covariance tangents before calling the existing bounded inverse.

For assembled W=sum_beta w_beta and a public velocity tangent S_j,

    d_j C_theta = <S_j,r W_theta + W_r S_j,theta>,
    d_j C_z     = <S_j,r W_z     + W_r S_j,z>.

For the common coordinate S_common=W, hence d_common C=2 C exactly.  If a
one-column calibration is materialized as W=a*w_unit, the old physical-amplitude
budget B_a and the fractional budget B_delta represent the same velocity update
only when

    B_delta = B_a / |a|.

That conversion is repository algebra, not a Kokuno source constant or a new
scientific threshold.  A future differently normalized physical family must
supply its own explicit normalization before this bridge may declare budget-unit
compatibility.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from math import pi
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

from .kokuno_bounded_family_coefficients import KokunoBoundedFamilyCoefficientCoordinates
from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_covariance_tangent_preflight import complete_curl_phase_evaluator, _ring_points
from .kokuno_family_bounded_inverse import (
    _receipt_from_target_report,
    evaluate_family_bounded_inverse,
)
from .kokuno_family_covariance_cross_terms import single_column_family_evaluator
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    PROFILE_Z,
    generate_actual_core_report as generate_missing_target_report,
)
from .kokuno_second_column_bounded_inverse import current_signed_coefficient_budget
from .kokuno_signed_covariance_inverse import (
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)

TASK = "KOKUNO-A3-BOUNDED-COORDINATE-UNIT-BRIDGE-020"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
AGENT2_COORDINATE_PR = 400
AGENT2_COORDINATE_HEAD = "d4464acc3e7072f473746df9b02df823cdbcb8c4"
AGENT3_INVERSE_PR = 393
AGENT3_INVERSE_HEAD = "bef9206d3aac61f2de8f0c479dfa74cc40938746"
EXPECTED_PARAMETER_NAMES = ("delta_common", "delta_band")
EXPECTED_PARAMETER_UNIT = "dimensionless_fractional_multiplier_of_Q_scaled_physical_velocity"

FamilyPhaseEvaluator = Callable[[np.ndarray, float, float], Mapping[str, Any]]


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.shape[-1:] != (2,) or not np.isfinite(values).all():
        raise ValueError("values must be finite with final dimension 2")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))


def convert_amplitude_budget_to_fractional(
    frozen_amplitude_budget: float,
    reference_physical_amplitude: float,
    *,
    maximum_fractional_budget: float = 0.5,
) -> float:
    """Convert a frozen physical-amplitude update budget into fractional units.

    The conversion is valid only when the supplied physical family is known to
    equal ``reference_physical_amplitude * unit_family``.  It does not infer a
    normalization for a future source family.
    """
    budget = float(frozen_amplitude_budget)
    amplitude = float(reference_physical_amplitude)
    maximum = float(maximum_fractional_budget)
    if not np.isfinite([budget, amplitude, maximum]).all():
        raise ValueError("budget/amplitude/maximum must be finite")
    if budget <= 0.0 or amplitude == 0.0:
        raise ValueError("frozen budget must be positive and reference amplitude nonzero")
    if not (0.0 < maximum <= 0.5):
        raise ValueError("maximum_fractional_budget must lie in (0,0.5]")
    converted = budget / abs(amplitude)
    if converted > maximum + 8.0 * np.finfo(float).eps:
        raise ValueError(
            "converted fractional budget exceeds the bounded-coordinate positivity contract"
        )
    return float(converted)


def covariance_tangent_from_velocity_and_tangents(
    base_velocity: np.ndarray,
    tangent_velocity: np.ndarray,
    *,
    sample_axis: int = 0,
) -> np.ndarray:
    """Return d(<W_r W_theta>,<W_r W_z>) for one or more velocity tangents."""
    base = np.asarray(base_velocity, dtype=float)
    tangent = np.asarray(tangent_velocity, dtype=float)
    if base.ndim < 2 or base.shape[-1] != 3 or not np.isfinite(base).all():
        raise ValueError("base_velocity must be finite with final dimension 3")
    if tangent.shape[:-2] != base.shape[:-1] or tangent.shape[-1] != 3:
        raise ValueError("tangent_velocity must have shape base_sample_shape+(parameters,3)")
    if not np.isfinite(tangent).all():
        raise ValueError("tangent_velocity must be finite")
    axis = int(sample_axis)
    if axis < 0:
        axis += base.ndim - 1
    if not 0 <= axis < base.ndim - 1:
        raise ValueError("sample_axis must index a base sample dimension")

    wr = base[..., 0]
    wtheta = base[..., 1]
    wz = base[..., 2]
    sr = tangent[..., 0]
    stheta = tangent[..., 1]
    sz = tangent[..., 2]
    theta_part = np.mean(sr * wtheta[..., None] + wr[..., None] * stheta, axis=axis)
    axial_part = np.mean(sr * wz[..., None] + wr[..., None] * sz, axis=axis)
    return np.stack((theta_part, axial_part), axis=-1)


def measure_bounded_coordinate_covariance(
    family_velocity: FamilyPhaseEvaluator,
    radii: np.ndarray,
    *,
    coordinate_contract: KokunoBoundedFamilyCoefficientCoordinates,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Measure #400 coordinate covariance tangents including label cross terms."""
    if not callable(family_velocity):
        raise TypeError("family_velocity must be callable")
    if not isinstance(coordinate_contract, KokunoBoundedFamilyCoefficientCoordinates):
        raise TypeError("coordinate_contract must be KokunoBoundedFamilyCoefficientCoordinates")
    radii = np.asarray(radii, dtype=float)
    if radii.ndim != 1 or len(radii) < 3 or not np.isfinite(radii).all():
        raise ValueError("radii must be a finite one-dimensional array")
    if angular_count < 8 or angular_count % 2:
        raise ValueError("angular_count must be an even integer >= 8")
    if phase_count < 4 or phase_count % 2:
        raise ValueError("phase_count must be an even integer >= 4")
    if not np.isfinite([time, z]).all():
        raise ValueError("time and z must be finite")

    points = _ring_points(radii, angular_count=angular_count, z=float(z))
    flat = points.reshape(-1, 3)
    shifts = 2.0 * pi * np.arange(phase_count, dtype=float) / phase_count
    covariance = np.zeros((len(radii), 2), dtype=float)
    response = np.zeros((len(radii), 2, 2), dtype=float)
    labels: tuple[str, ...] | None = None
    band_active: bool | None = None

    for shift in shifts:
        sample = dict(family_velocity(flat, float(time), float(shift)))
        out = coordinate_contract.consume_physical_family(sample)
        names = tuple(str(x) for x in out["parameter_names"])
        units = tuple(str(x) for x in out["parameter_units"])
        if names != EXPECTED_PARAMETER_NAMES:
            raise ValueError("unexpected bounded-coordinate parameter names")
        if units != (EXPECTED_PARAMETER_UNIT, EXPECTED_PARAMETER_UNIT):
            raise ValueError("unexpected bounded-coordinate parameter units")
        current_labels = tuple(str(x) for x in out["beta_labels"])
        if labels is None:
            labels = current_labels
            band_active = bool(out["band_contrast_active"])
        elif current_labels != labels or bool(out["band_contrast_active"]) != band_active:
            raise ValueError("bounded family labels/active-band state changed across phase samples")

        base = np.asarray(out["velocity_physical_cylindrical_total_base"], dtype=float)
        tangents = np.asarray(out["velocity_parameter_tangent_cylindrical"], dtype=float)
        base = base.reshape(len(radii), angular_count, 3)
        tangents = tangents.reshape(len(radii), angular_count, 2, 3)
        covariance[:, 0] += np.mean(base[..., 0] * base[..., 1], axis=1)
        covariance[:, 1] += np.mean(base[..., 0] * base[..., 2], axis=1)
        response += covariance_tangent_from_velocity_and_tangents(
            base, tangents, sample_axis=1
        )

    covariance /= float(phase_count)
    response /= float(phase_count)
    common_error = response[:, 0, :] - 2.0 * covariance
    common_relative_error = _vector_rms(common_error) / max(
        2.0 * _vector_rms(covariance), np.finfo(float).tiny
    )
    if common_relative_error > 5.0e-12:
        raise RuntimeError("delta_common covariance tangent failed the exact quadratic identity")

    assert labels is not None and band_active is not None
    return {
        "radii": radii.tolist(),
        "beta_labels": list(labels),
        "parameter_names": list(EXPECTED_PARAMETER_NAMES),
        "parameter_units": [EXPECTED_PARAMETER_UNIT, EXPECTED_PARAMETER_UNIT],
        "repository_coefficient_unit_mapping_available": True,
        "source_coefficient_unit_mapping_available": False,
        "band_contrast_active": band_active,
        "base_covariance_theta_axial": covariance.tolist(),
        "coordinate_covariance_jacobian_theta_axial": response.tolist(),
        "common_quadratic_identity_relative_error": float(common_relative_error),
    }


def _one_label_real_family() -> FamilyPhaseEvaluator:
    unit = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=1.0,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    wrapped = single_column_family_evaluator(
        complete_curl_phase_evaluator(unit),
        physical_amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        label="existing_complete_curl",
    )

    def evaluate(points: np.ndarray, time: float, phase_offset: float) -> Mapping[str, Any]:
        out = dict(wrapped(points, time, phase_offset))
        # The numeric ell label is an autonomous calibration label only.  A
        # one-band family has no delta_band direction regardless of its value.
        out["beta_labels"] = ((5, ("existing_complete_curl", 0, 0)),)
        return out

    return evaluate


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/bounded_coordinate_covariance_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/bounded_coordinate_covariance_target.json",
) -> dict[str, Any]:
    """Run the real-defect one-band calibration through #400 -> #393."""
    target_report = generate_missing_target_report(output=target_output)
    receipt = _receipt_from_target_report(target_report)
    frozen_amplitude_budget = current_signed_coefficient_budget(receipt)
    fractional_budget = convert_amplitude_budget_to_fractional(
        frozen_amplitude_budget,
        ROUTED_OSCILLATORY_AMPLITUDE,
    )
    coordinates = KokunoBoundedFamilyCoefficientCoordinates(max_l1_update=fractional_budget)
    measurement = measure_bounded_coordinate_covariance(
        _one_label_real_family(),
        receipt.radii,
        coordinate_contract=coordinates,
        time=float(target_report["inputs"]["profile_time"]),
        z=float(target_report["inputs"]["profile_z"]),
    )
    jacobian = np.asarray(measurement["coordinate_covariance_jacobian_theta_axial"], dtype=float)
    expected_common = ROUTED_OSCILLATORY_AMPLITUDE * receipt.current_response
    common_calibration_error = _vector_rms(jacobian[:, 0, :] - expected_common) / max(
        _vector_rms(expected_common), np.finfo(float).tiny
    )
    band_rms = _vector_rms(jacobian[:, 1, :])

    inverse = evaluate_family_bounded_inverse(
        receipt,
        jacobian,
        coefficient_labels=EXPECTED_PARAMETER_NAMES,
        current_coefficient_budget=frozen_amplitude_budget,
        additional_family_l1_budget=fractional_budget,
        budget_unit_matches_jacobian=True,
        coefficient_unit=EXPECTED_PARAMETER_UNIT,
    )

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "source_scope": (
                "Kokuno motivates the signed covariance/stress correction structure. "
                "The two bounded coefficient coordinates and this unit conversion are repository engineering."
            ),
        },
        "handoff": {
            "agent2_coordinate_pr": AGENT2_COORDINATE_PR,
            "agent2_coordinate_head": AGENT2_COORDINATE_HEAD,
            "agent3_inverse_pr": AGENT3_INVERSE_PR,
            "agent3_inverse_head": AGENT3_INVERSE_HEAD,
        },
        "inputs": {
            "leading_candidate_sha256": target_report["inputs"]["leading_candidate_sha256"],
            "oscillatory_amplitude": ROUTED_OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": ROUTED_OSCILLATORY_PHASE,
            "profile_time": target_report["inputs"]["profile_time"],
            "profile_z": target_report["inputs"]["profile_z"],
            "profile_annulus": target_report["inputs"]["profile_annulus"],
            "surrogate_defect_used": False,
            "frozen_physical_amplitude_budget": frozen_amplitude_budget,
            "converted_fractional_budget": fractional_budget,
            "budget_conversion": "B_fractional=B_physical_amplitude/abs(reference_physical_amplitude)",
            "budget_changed": False,
        },
        "measurement": measurement,
        "calibration": {
            "expected_common_response_scale_vs_old_amplitude_response": ROUTED_OSCILLATORY_AMPLITUDE,
            "common_coordinate_relative_error_to_scaled_existing_response": float(common_calibration_error),
            "band_coordinate_response_vector_rms": float(band_rms),
            "one_band_band_coordinate_inactive": not bool(measurement["band_contrast_active"]),
        },
        "bounded_inverse": inverse,
        "truth_boundary": {
            "real_candidate_defect_target_consumed": True,
            "surrogate_defect_used": False,
            "repository_coefficient_unit_mapping_consumed": True,
            "source_coefficient_unit_mapping_available": False,
            "actual_source_multiband_family_consumed": False,
            "genuinely_independent_second_covariance_column_ready": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
            "finite_correction_cycle_rerun_allowed": False,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
        "interpretation": (
            "The new repository coordinate units can now be consumed without silently reusing a physical-amplitude "
            "budget as a dimensionless number.  The actual one-band control remains rank deficient, so no finite "
            "correction cycle is authorized."
        ),
    }
    if inverse["family_bounded_inverse_preflight_passed"]:
        raise RuntimeError("one-band bounded-coordinate control unexpectedly passed the rank-two preflight")
    if common_calibration_error > 5.0e-11 or band_rms > 5.0e-14:
        raise RuntimeError("real one-band bounded-coordinate calibration failed")

    output_path = Path(output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/bounded_coordinate_covariance_report.json",
    )
    parser.add_argument(
        "--target-output",
        default="artifacts/kokuno_agent3/bounded_coordinate_covariance_target.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output, target_output=args.target_output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
