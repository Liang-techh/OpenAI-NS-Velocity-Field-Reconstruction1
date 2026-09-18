"""Retain overlap cross terms when a Kokuno real-pair family enters mean correction.

Agent 2 PR #380 exposes a supplied-data localized real-pair *family* with
physical cylindrical velocities by slow label ``beta`` and their assembled
sum.  Agent 3 owns the zero-harmonic covariance/stress layer.  For an assembled
field

    W = sum_beta w_beta,

the mean covariance is quadratic, so overlapping labels cannot be screened by
summing only diagonal/self covariances.  The repository algebra is

    C_theta = <W_r W_theta>,       C_z = <W_r W_z>,

and, for a label coefficient ``a_beta`` with public tangent ``s_beta``,

    d_beta C_theta = <s_beta,r W_theta + W_r s_beta,theta>,
    d_beta C_z     = <s_beta,r W_z     + W_r s_beta,z>.

If ``a_beta`` is the dimensionless multiplier of the already-materialized
label contribution, then ``s_beta=w_beta``.  In that special case the exact
quadratic identity

    sum_beta d_beta C = 2 C

must hold.  A different source coefficient normalization can be represented by
caller-supplied tangent scales, but its units must be kept explicit before any
existing signed-coefficient budget is reused.

These product-rule and self/cross decompositions are repository algebra.  They
are not claimed as new Kokuno source formulas.  This module consumes an
Agent-2-style by-label velocity interface and does not reimplement pulse,
partition, localization, or complete-curl construction.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from math import pi
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_covariance_tangent_preflight import (
    _cylindrical_velocity,
    _ring_points,
    complete_curl_phase_evaluator,
)
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    PROFILE_Z,
    generate_actual_core_report as generate_missing_target_report,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)

TASK = "KOKUNO-A3-FAMILY-COVARIANCE-CROSS-TERMS-018"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
AGENT2_FAMILY_PR = 380
AGENT2_FAMILY_HEAD = "43f021e904dbb09c9bc2aba1d9ea2954083937e1"

FamilyPhaseEvaluator = Callable[[np.ndarray, float, float], Mapping[str, Any]]
PhaseVelocityEvaluator = Callable[[np.ndarray, float, float], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    if values.shape[-1:] != (2,) or not np.isfinite(values).all():
        raise ValueError("values must be finite with final dimension 2")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))


def _validate_scales(scales: Sequence[float] | np.ndarray | None, beta_count: int) -> np.ndarray:
    if scales is None:
        return np.ones(beta_count, dtype=float)
    out = np.asarray(scales, dtype=float)
    if out.shape != (beta_count,) or not np.isfinite(out).all():
        raise ValueError("coefficient_tangent_scales must be finite with shape (beta_count,)")
    if np.any(out == 0.0):
        raise ValueError("coefficient_tangent_scales must be nonzero")
    return out


def _validate_family_sample(
    sample: Mapping[str, Any],
    sample_count: int,
    *,
    expected_beta_count: int | None = None,
    total_rtol: float = 2.0e-12,
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...]]:
    if not isinstance(sample, Mapping):
        raise TypeError("family evaluator must return a mapping")
    required = {
        "velocity_physical_cylindrical_by_beta",
        "velocity_physical_cylindrical_total",
    }
    if not required.issubset(sample):
        raise ValueError("family evaluator is missing physical cylindrical by-beta/total velocity")
    by_beta = np.asarray(sample["velocity_physical_cylindrical_by_beta"], dtype=float)
    total = np.asarray(sample["velocity_physical_cylindrical_total"], dtype=float)
    if by_beta.ndim != 3 or by_beta.shape[0] != sample_count or by_beta.shape[-1] != 3:
        raise ValueError("by-beta velocity must have shape (N,beta_count,3)")
    if total.shape != (sample_count, 3):
        raise ValueError("total velocity must have shape (N,3)")
    if not np.isfinite(by_beta).all() or not np.isfinite(total).all():
        raise ValueError("family velocities must be finite")
    beta_count = by_beta.shape[1]
    if beta_count < 1:
        raise ValueError("family must contain at least one beta label")
    if expected_beta_count is not None and beta_count != expected_beta_count:
        raise ValueError("family beta count changed across phase samples")
    summed = np.sum(by_beta, axis=1)
    scale = max(1.0, float(np.max(np.abs(total))), float(np.max(np.abs(summed))))
    if not np.allclose(total, summed, rtol=0.0, atol=total_rtol * scale):
        raise ValueError("declared total velocity does not equal the sum of by-beta velocities")
    labels_raw = sample.get("beta_labels")
    if labels_raw is None:
        labels = tuple(f"beta_{index}" for index in range(beta_count))
    else:
        if not isinstance(labels_raw, Sequence) or isinstance(labels_raw, (str, bytes)):
            raise ValueError("beta_labels must be a sequence when supplied")
        labels = tuple(str(label) for label in labels_raw)
        if len(labels) != beta_count or len(set(labels)) != beta_count:
            raise ValueError("beta_labels must be unique and match beta_count")
    return by_beta, total, labels


def single_column_family_evaluator(
    velocity: PhaseVelocityEvaluator,
    *,
    physical_amplitude: float,
    label: str = "existing_complete_curl",
) -> FamilyPhaseEvaluator:
    """Wrap one public Cartesian velocity column as a one-label physical family.

    The returned label contribution includes ``physical_amplitude``.  To recover
    the derivative with respect to that physical amplitude, call
    :func:`measure_family_covariance_jacobian` with tangent scale ``1/amplitude``.
    """
    if not callable(velocity):
        raise TypeError("velocity must be callable")
    physical_amplitude = float(physical_amplitude)
    if not np.isfinite(physical_amplitude) or physical_amplitude == 0.0:
        raise ValueError("physical_amplitude must be finite and nonzero")
    label = str(label)
    if not label:
        raise ValueError("label must be nonempty")

    def evaluate(points: np.ndarray, time: float, phase_offset: float) -> Mapping[str, Any]:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must have finite shape (N,3)")
        cart = np.asarray(velocity(points, float(time), float(phase_offset)), dtype=float)
        if cart.shape != points.shape or not np.isfinite(cart).all():
            raise ValueError("velocity must return finite shape-(N,3) Cartesian values")
        cart = physical_amplitude * cart
        ur, utheta, uz = _cylindrical_velocity(cart, points)
        cyl = np.stack((ur, utheta, uz), axis=-1)
        return {
            "beta_labels": (label,),
            "velocity_physical_cylindrical_by_beta": cyl[:, None, :],
            "velocity_physical_cylindrical_total": cyl,
        }

    return evaluate


def measure_family_covariance_jacobian(
    family_velocity: FamilyPhaseEvaluator,
    radii: np.ndarray,
    *,
    coefficient_tangent_scales: Sequence[float] | np.ndarray | None = None,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
    total_rtol: float = 2.0e-12,
) -> dict[str, Any]:
    """Measure assembled covariance, self/cross split, and per-label Jacobian.

    ``coefficient_tangent_scales=None`` means the coefficient for every label is
    its dimensionless multiplier, so ``s_beta=w_beta`` and the Euler quadratic
    identity is enforced.  Non-unit scales change coefficient units and disable
    that particular identity; the scale vector is returned verbatim so later
    bounded-inverse code cannot silently forget the normalization change.
    """
    if not callable(family_velocity):
        raise TypeError("family_velocity must be callable")
    radii = np.asarray(radii, dtype=float)
    if radii.ndim != 1 or len(radii) < 3 or not np.isfinite(radii).all():
        raise ValueError("radii must be a finite one-dimensional array")
    if phase_count < 4 or phase_count % 2:
        raise ValueError("phase_count must be an even integer >= 4")
    if angular_count < 8 or angular_count % 2:
        raise ValueError("angular_count must be an even integer >= 8")
    if not np.isfinite([time, z, total_rtol]).all() or total_rtol <= 0.0:
        raise ValueError("time/z must be finite and total_rtol positive")

    points = _ring_points(radii, angular_count=angular_count, z=float(z))
    flat = points.reshape(-1, 3)
    shifts = 2.0 * pi * np.arange(phase_count, dtype=float) / phase_count

    total_covariance: np.ndarray | None = None
    diagonal_covariance: np.ndarray | None = None
    jacobian: np.ndarray | None = None
    labels: tuple[str, ...] | None = None
    scales: np.ndarray | None = None

    for shift in shifts:
        sample = family_velocity(flat, float(time), float(shift))
        by_beta, total, current_labels = _validate_family_sample(
            sample,
            len(flat),
            expected_beta_count=None if labels is None else len(labels),
            total_rtol=total_rtol,
        )
        if labels is None:
            labels = current_labels
            scales = _validate_scales(coefficient_tangent_scales, len(labels))
            total_covariance = np.zeros((len(radii), 2), dtype=float)
            diagonal_covariance = np.zeros((len(radii), 2), dtype=float)
            jacobian = np.zeros((len(radii), len(labels), 2), dtype=float)
        elif current_labels != labels:
            raise ValueError("beta labels changed across phase samples")

        assert scales is not None
        assert total_covariance is not None
        assert diagonal_covariance is not None
        assert jacobian is not None

        beta_count = len(labels)
        by_beta = by_beta.reshape(len(radii), angular_count, beta_count, 3)
        total = total.reshape(len(radii), angular_count, 3)
        wr = by_beta[..., 0]
        wtheta = by_beta[..., 1]
        wz = by_beta[..., 2]
        Wr = total[..., 0]
        Wtheta = total[..., 1]
        Wz = total[..., 2]

        total_covariance[:, 0] += np.mean(Wr * Wtheta, axis=1)
        total_covariance[:, 1] += np.mean(Wr * Wz, axis=1)
        diagonal_covariance[:, 0] += np.sum(np.mean(wr * wtheta, axis=1), axis=1)
        diagonal_covariance[:, 1] += np.sum(np.mean(wr * wz, axis=1), axis=1)

        sr = wr * scales[None, None, :]
        stheta = wtheta * scales[None, None, :]
        sz = wz * scales[None, None, :]
        jacobian[:, :, 0] += np.mean(
            sr * Wtheta[:, :, None] + Wr[:, :, None] * stheta,
            axis=1,
        )
        jacobian[:, :, 1] += np.mean(
            sr * Wz[:, :, None] + Wr[:, :, None] * sz,
            axis=1,
        )

    assert labels is not None
    assert scales is not None
    assert total_covariance is not None
    assert diagonal_covariance is not None
    assert jacobian is not None
    divisor = float(phase_count)
    total_covariance /= divisor
    diagonal_covariance /= divisor
    jacobian /= divisor
    cross_covariance = total_covariance - diagonal_covariance

    total_norm = max(_vector_rms(total_covariance), np.finfo(float).tiny)
    cross_fraction = _vector_rms(cross_covariance) / total_norm
    unit_multiplier_semantics = bool(np.array_equal(scales, np.ones_like(scales)))
    if unit_multiplier_semantics:
        euler_error = np.sum(jacobian, axis=1) - 2.0 * total_covariance
        euler_relative_error = _vector_rms(euler_error) / max(
            2.0 * total_norm, np.finfo(float).tiny
        )
        if euler_relative_error > 5.0e-12:
            raise RuntimeError("family covariance Jacobian failed the quadratic Euler identity")
    else:
        euler_relative_error = None

    return {
        "radii": radii.tolist(),
        "beta_labels": list(labels),
        "coefficient_tangent_scales": scales.tolist(),
        "coefficient_semantics": (
            "dimensionless label multipliers" if unit_multiplier_semantics
            else "caller-scaled label tangents; coefficient units must be preserved downstream"
        ),
        "total_covariance_theta_axial": total_covariance.tolist(),
        "diagonal_self_covariance_theta_axial": diagonal_covariance.tolist(),
        "cross_covariance_theta_axial": cross_covariance.tolist(),
        "cross_covariance_relative_vector_rms": float(cross_fraction),
        "per_label_covariance_jacobian_theta_axial": jacobian.tolist(),
        "quadratic_euler_identity_checked": unit_multiplier_semantics,
        "quadratic_euler_relative_error": euler_relative_error,
        "bounded_inverse_ready": False,
        "bounded_inverse_blocker": (
            "A source coefficient-unit mapping and an actually instantiated independent "
            "family are still required before reusing Agent-3 signed-coefficient budgets."
        ),
    }


def _target_current_response(target_report: Mapping[str, Any]) -> np.ndarray:
    target = target_report["missing_second_column_target"]
    return np.asarray(target["current_response_theta_axial"], dtype=float)


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/family_covariance_cross_terms_report.json",
    target_output: str | Path = "artifacts/kokuno_agent3/family_covariance_cross_terms_target.json",
) -> dict[str, Any]:
    """Calibrate the family contract on the existing one-label physical column.

    This is a real-defect calibration only.  It deliberately does not fabricate
    the missing multi-label source family: one label has zero cross covariance.
    The tangent scale ``1/a`` maps the materialized contribution ``a*w`` back to
    the existing physical-amplitude derivative used by Agent 3 #322.
    """
    target_report = generate_missing_target_report(output=target_output)
    radii = np.asarray(target_report["missing_second_column_target"]["radii"], dtype=float)
    unit = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=1.0,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    family = single_column_family_evaluator(
        complete_curl_phase_evaluator(unit),
        physical_amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
    )
    measurement = measure_family_covariance_jacobian(
        family,
        radii,
        coefficient_tangent_scales=[1.0 / ROUTED_OSCILLATORY_AMPLITUDE],
        time=PROFILE_TIME,
        z=PROFILE_Z,
        angular_count=int(target_report["inputs"]["angular_count"]),
        phase_count=int(target_report["inputs"]["phase_count"]),
    )
    jacobian = np.asarray(
        measurement["per_label_covariance_jacobian_theta_axial"], dtype=float
    )[:, 0, :]
    reference = _target_current_response(target_report)
    active = np.asarray(
        target_report["missing_second_column_target"]["active_target_mask"], dtype=bool
    )
    scale = max(_vector_rms(reference[active]), np.finfo(float).tiny)
    relative_error = _vector_rms((jacobian - reference)[active]) / scale
    cross = np.asarray(measurement["cross_covariance_theta_axial"], dtype=float)

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "source_scope": "slow-label real-pair family composition and signed covariance correction",
            "repository_algebra": (
                "assembled covariance is quadratic in the sum over beta; the per-label "
                "product-rule Jacobian and self/cross decomposition are Agent-3 algebra"
            ),
        },
        "agent2_family_handoff": {
            "pr": AGENT2_FAMILY_PR,
            "exact_head": AGENT2_FAMILY_HEAD,
            "by_beta_interface": "velocity_physical_cylindrical_by_beta",
            "total_interface": "velocity_physical_cylindrical_total",
            "supplied_data_family_composition_available": True,
            "concrete_source_partition_bumps_reconstructed": False,
            "source_actual_partition_labels_instantiated": False,
            "source_actual_background_path_instantiated": False,
            "genuinely_independent_second_covariance_column_ready": False,
        },
        "inputs": {
            "leading_candidate_sha256": target_report["inputs"]["leading_candidate_sha256"],
            "profile_time": PROFILE_TIME,
            "profile_z": PROFILE_Z,
            "profile_annulus": list(PROFILE_ANNULUS),
            "oscillatory_amplitude": ROUTED_OSCILLATORY_AMPLITUDE,
            "oscillatory_phase": ROUTED_OSCILLATORY_PHASE,
            "surrogate_defect_used": False,
        },
        "one_label_real_calibration": {
            "coefficient_tangent_scale": 1.0 / ROUTED_OSCILLATORY_AMPLITUDE,
            "relative_error_to_existing_current_response": float(relative_error),
            "cross_covariance_vector_rms": _vector_rms(cross),
            "expected_cross_covariance": 0.0,
            "interpretation": (
                "The new family Jacobian reproduces the already-measured one-column "
                "physical-amplitude response. This validates the handoff only; it is not "
                "a new second direction or a residual improvement."
            ),
        },
        "family_measurement": measurement,
        "truth_boundary": {
            "real_candidate_defect_target_consumed": True,
            "surrogate_defect_used": False,
            "agent2_family_code_reimplemented": False,
            "family_cross_term_contract_executable": True,
            "actual_source_multilabel_family_consumed": False,
            "cross_label_source_covariance_measured": False,
            "source_coefficient_unit_mapping_available": False,
            "second_public_covariance_column_available": False,
            "bounded_inverse_rerun_allowed": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
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
        description="Audit overlap cross terms in Agent-2-style Kokuno real-pair families"
    )
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/family_covariance_cross_terms_report.json",
    )
    parser.add_argument(
        "--target-output",
        default="artifacts/kokuno_agent3/family_covariance_cross_terms_target.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output, target_output=args.target_output)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
