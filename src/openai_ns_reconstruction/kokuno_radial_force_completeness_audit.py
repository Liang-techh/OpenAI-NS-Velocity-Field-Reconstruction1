"""Audit the source-required radial force from Agent-3 compact axial stress.

The existing Agent-3 radial reconstruction uses the corrected-reader moment
complement for the theta/e=2 and z/e=1 channels, but its second-column routing
has so far screened only those two tangential stress channels.  The corrected
reader's full symmetric stress tensor has an additional physical radial force

    (div T)_r = d_z sigma_1,

where ``sigma_1`` is the e=1 axial radial stress.  That term is explicitly
nonzero in the source and therefore must not be silently dropped when deciding
whether a future two-column correction is safe to materialize.

This module consumes the same *real* phase-mean Navier--Stokes defect as the
existing Agent-3 machinery, reconstructs ``sigma_1`` at neighboring physical z
slices, and differentiates it on a frozen centered-step ladder.  The derivative
ladder is an autonomous numerical audit, not a Navier--Stokes acceptance gate.
No second oscillatory column is invented, no forcing or pressure is fitted, and
no finite correction cycle is rerun here.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_missing_covariance_column_target import (
    DEFAULT_ANGULAR_COUNT,
    DEFAULT_PHASE_COUNT,
    DEFAULT_RADIAL_COUNT,
    PROFILE_Z,
)
from .kokuno_radial_stress import (
    CompactRadialStressInverse,
    sample_real_phase_mean_radial_profiles,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)

TASK = "KOKUNO-A3-RADIAL-FORCE-COMPLETENESS-AUDIT-016"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_FORMULAS = "R33-R34 and R41"

# Frozen before measuring the real field.  These are ordinary numerical
# differentiation controls, not scientific acceptance thresholds.
Z_DERIVATIVE_STEP_LADDER = (0.02, 0.01, 0.005)
FINE_PAIR_RELATIVE_STABILITY_TOLERANCE = 5.0e-2


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _vector_rms(*components: np.ndarray) -> float:
    arrays = [np.asarray(component, dtype=float) for component in components]
    if not arrays or any(array.shape != arrays[0].shape for array in arrays):
        raise ValueError("components must be nonempty arrays with matching shapes")
    total = np.zeros_like(arrays[0], dtype=float)
    for array in arrays:
        total += array * array
    return float(np.sqrt(np.mean(total)))


def _validate_z_steps(steps: Iterable[float]) -> tuple[float, ...]:
    values = tuple(float(step) for step in steps)
    if len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("z derivative ladder must contain at least two finite steps")
    if any(step <= 0.0 for step in values):
        raise ValueError("z derivative steps must be positive")
    if any(right >= left for left, right in zip(values, values[1:])):
        raise ValueError("z derivative steps must be strictly decreasing")
    return values


def centered_radial_force(
    sigma_plus: np.ndarray,
    sigma_minus: np.ndarray,
    step: float,
) -> np.ndarray:
    """Return the centered approximation to the source term ``d_z sigma_1``."""
    plus = np.asarray(sigma_plus, dtype=float)
    minus = np.asarray(sigma_minus, dtype=float)
    step = float(step)
    if plus.shape != minus.shape or plus.ndim != 1 or len(plus) < 3:
        raise ValueError("sigma_plus and sigma_minus must be matching 1-D arrays")
    if not np.isfinite(plus).all() or not np.isfinite(minus).all():
        raise ValueError("stress arrays must be finite")
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    return (plus - minus) / (2.0 * step)


def audit_radial_force_derivative_ladder(
    steps: Iterable[float],
    derivatives: Iterable[np.ndarray],
    *,
    active_mask: np.ndarray | None = None,
    fine_pair_relative_stability_tolerance: float = FINE_PAIR_RELATIVE_STABILITY_TOLERANCE,
) -> dict[str, Any]:
    """Check that the reconstructed radial-force derivative is numerically stable."""
    checked_steps = _validate_z_steps(steps)
    arrays = tuple(np.asarray(value, dtype=float) for value in derivatives)
    if len(arrays) != len(checked_steps):
        raise ValueError("one derivative array is required for every z step")
    if not arrays or arrays[0].ndim != 1:
        raise ValueError("derivative arrays must be one-dimensional")
    shape = arrays[0].shape
    if any(array.shape != shape or not np.isfinite(array).all() for array in arrays):
        raise ValueError("derivative arrays must be finite and have one common shape")

    if active_mask is None:
        mask = np.ones(shape, dtype=bool)
    else:
        mask = np.asarray(active_mask, dtype=bool)
        if mask.shape != shape or not np.any(mask):
            raise ValueError("active_mask must select at least one derivative entry")

    tolerance = float(fine_pair_relative_stability_tolerance)
    if not np.isfinite(tolerance) or not 0.0 < tolerance < 0.5:
        raise ValueError("fine-pair stability tolerance must lie in (0,0.5)")

    levels = []
    for step, derivative in zip(checked_steps, arrays):
        selected = derivative[mask]
        levels.append(
            {
                "z_step": step,
                "radial_force_rms": _rms(selected),
                "radial_force_max_abs": float(np.max(np.abs(selected))),
            }
        )

    comparisons = []
    for coarse_step, coarse, fine_step, fine in zip(
        checked_steps,
        arrays,
        checked_steps[1:],
        arrays[1:],
    ):
        fine_scale = max(_rms(fine[mask]), np.finfo(float).tiny)
        comparisons.append(
            {
                "coarse_z_step": coarse_step,
                "fine_z_step": fine_step,
                "relative_rms_difference": _rms((coarse - fine)[mask]) / fine_scale,
            }
        )

    finest_difference = float(comparisons[-1]["relative_rms_difference"])
    return {
        "z_steps": list(checked_steps),
        "levels": levels,
        "comparisons": comparisons,
        "fine_pair_relative_stability_tolerance": tolerance,
        "finest_pair_relative_rms_difference": finest_difference,
        "radial_force_derivative_stability_preflight_passed": bool(
            finest_difference <= tolerance
        ),
        "interpretation": (
            "This checks only numerical stability of the source-required d_z sigma_1 "
            "side effect. It is not a PDE residual threshold and does not authorize a "
            "correction cycle by itself."
        ),
    }


def _annular_window(radii: np.ndarray, r_inner: float, r_outer: float) -> np.ndarray:
    radii = np.asarray(radii, dtype=float)
    s = 2.0 * (radii - r_inner) / (r_outer - r_inner) - 1.0
    inside = np.abs(s) < 1.0
    base = np.where(inside, 1.0 - s * s, 0.0)
    return base**5


def sample_real_radial_mean_defect(
    base_velocity: Callable[[np.ndarray, float], np.ndarray],
    correction: KokunoCompleteCurlCorrection,
    contract: PhaseMeanDefectContract,
    radii: np.ndarray,
    *,
    time: float,
    z: float,
    angular_count: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Ring-average the radial component of the same real phase-mean defect."""
    radii = np.asarray(radii, dtype=float)
    if radii.ndim != 1 or len(radii) < 9 or np.any(radii <= 0.0):
        raise ValueError("radii must be a positive 1-D grid with at least 9 nodes")
    if not np.isfinite(radii).all() or np.any(np.diff(radii) <= 0.0):
        raise ValueError("radii must be finite and strictly increasing")
    if angular_count < 8 or angular_count % 2:
        raise ValueError("angular_count must be an even integer >= 8")

    angles = 2.0 * np.pi * np.arange(angular_count, dtype=float) / angular_count
    rr, aa = np.meshgrid(radii, angles, indexing="ij")
    points = np.stack(
        (rr * np.cos(aa), rr * np.sin(aa), np.full_like(rr, float(z))), axis=-1
    )
    flat_points = points.reshape(-1, 3)
    sample = contract.evaluate(base_velocity, correction, flat_points, float(time))
    increment = np.asarray(sample.mean_defect_increment, dtype=float)
    radius = np.hypot(flat_points[:, 0], flat_points[:, 1])
    e_r = np.column_stack(
        (
            flat_points[:, 0] / radius,
            flat_points[:, 1] / radius,
            np.zeros(len(flat_points)),
        )
    )
    radial = np.einsum("ni,ni->n", increment, e_r).reshape(
        len(radii), angular_count
    )
    raw = np.mean(radial, axis=1)
    window = _annular_window(radii, float(radii[0]), float(radii[-1]))
    return raw, raw * window


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/radial_force_completeness_audit_report.json",
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
    z_steps: Iterable[float] = Z_DERIVATIVE_STEP_LADDER,
) -> dict[str, Any]:
    """Measure the omitted radial-force channel on the current real correction target."""
    checked_steps = _validate_z_steps(z_steps)
    if radial_count < 17:
        raise ValueError("radial_count must be at least 17")
    if angular_count < 8 or angular_count % 2:
        raise ValueError("angular_count must be an even integer >= 8")
    if phase_count < 4 or phase_count % 2:
        raise ValueError("phase_count must be an even integer >= 4")

    leading = KokunoLeadingCoreSeriesCandidate()
    correction = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=int(phase_count))
    r_inner, r_outer = (float(PROFILE_ANNULUS[0]), float(PROFILE_ANNULUS[1]))
    radii = np.linspace(r_inner, r_outer, int(radial_count))
    active_mask = (radii > r_inner) & (radii < r_outer)

    theta_profile, axial_profile, center_projection = sample_real_phase_mean_radial_profiles(
        leading.at_points,
        correction,
        contract,
        time=float(time),
        z=float(z),
        r_inner=r_inner,
        r_outer=r_outer,
        radial_count=int(radial_count),
        angular_count=int(angular_count),
    )
    theta_inverse = CompactRadialStressInverse(theta_profile)
    axial_inverse = CompactRadialStressInverse(axial_profile)

    tangential_theta_force = (
        -theta_inverse.source(radii)
        + theta_inverse.normalized_bump(radii) * theta_inverse.moment
    )
    tangential_axial_force = (
        -axial_inverse.source(radii)
        + axial_inverse.normalized_bump(radii) * axial_inverse.moment
    )

    derivative_levels: list[np.ndarray] = []
    slice_receipts: list[dict[str, Any]] = []
    for step in checked_steps:
        if not (correction.center[2] - correction.half_widths[2] < z - step < z + step < correction.center[2] + correction.half_widths[2]):
            raise ValueError("z derivative ladder must stay inside correction support")
        _, axial_plus, projection_plus = sample_real_phase_mean_radial_profiles(
            leading.at_points,
            correction,
            contract,
            time=float(time),
            z=float(z + step),
            r_inner=r_inner,
            r_outer=r_outer,
            radial_count=int(radial_count),
            angular_count=int(angular_count),
        )
        _, axial_minus, projection_minus = sample_real_phase_mean_radial_profiles(
            leading.at_points,
            correction,
            contract,
            time=float(time),
            z=float(z - step),
            r_inner=r_inner,
            r_outer=r_outer,
            radial_count=int(radial_count),
            angular_count=int(angular_count),
        )
        sigma_plus = CompactRadialStressInverse(axial_plus).stress(radii)
        sigma_minus = CompactRadialStressInverse(axial_minus).stress(radii)
        radial_force = centered_radial_force(sigma_plus, sigma_minus, step)
        derivative_levels.append(radial_force)
        slice_receipts.append(
            {
                "z_step": step,
                "z_minus": float(z - step),
                "z_plus": float(z + step),
                "plus_full_mean_defect_increment_rms": projection_plus[
                    "full_mean_defect_increment_rms"
                ],
                "minus_full_mean_defect_increment_rms": projection_minus[
                    "full_mean_defect_increment_rms"
                ],
            }
        )

    derivative_audit = audit_radial_force_derivative_ladder(
        checked_steps,
        derivative_levels,
        active_mask=active_mask,
    )
    finest_radial_force = derivative_levels[-1]
    raw_radial_defect, gated_radial_defect = sample_real_radial_mean_defect(
        leading.at_points,
        correction,
        contract,
        radii,
        time=float(time),
        z=float(z),
        angular_count=int(angular_count),
    )

    radial_force_rms = _rms(finest_radial_force[active_mask])
    radial_force_max = float(np.max(np.abs(finest_radial_force[active_mask])))
    raw_radial_rms = _rms(raw_radial_defect[active_mask])
    gated_radial_rms = _rms(gated_radial_defect[active_mask])
    tangential_force_rms = _vector_rms(
        tangential_theta_force[active_mask], tangential_axial_force[active_mask]
    )
    full_stress_force_rms = _vector_rms(
        finest_radial_force[active_mask],
        tangential_theta_force[active_mask],
        tangential_axial_force[active_mask],
    )
    scale = max(tangential_force_rms, np.finfo(float).tiny)
    nonzero_floor = 1024.0 * np.finfo(float).eps * max(full_stress_force_rms, 1.0)
    radial_force_nonzero = bool(radial_force_max > nonzero_floor)

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "formulas": SOURCE_FORMULAS,
            "source_structure": (
                "R33-R34 give div(T)_r=d_z sigma_1 in addition to the theta and z "
                "moment-complement rows. R41 shows the same radial term survives the "
                "physical-to-chart map with the source epsilon gain; it is not zero."
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading.sha256,
            "oscillatory_amplitude": float(correction.amplitude),
            "oscillatory_phase": float(correction.phase),
            "time": float(time),
            "z": float(z),
            "profile_annulus": [r_inner, r_outer],
            "radial_count": int(radial_count),
            "angular_count": int(angular_count),
            "phase_count": int(phase_count),
            "z_derivative_steps": list(checked_steps),
            "surrogate_defect_used": False,
        },
        "center_real_defect_projection": center_projection,
        "center_radial_stress_metrics": {
            "theta_e2": theta_inverse.metrics(),
            "axial_e1": axial_inverse.metrics(),
        },
        "z_derivative_slice_receipts": slice_receipts,
        "radial_force_derivative_audit": derivative_audit,
        "full_tensor_force_completeness": {
            "finest_radial_force_rms": radial_force_rms,
            "finest_radial_force_max_abs": radial_force_max,
            "raw_ring_radial_mean_defect_rms": raw_radial_rms,
            "gated_ring_radial_mean_defect_rms": gated_radial_rms,
            "tangential_reconstructed_force_vector_rms": tangential_force_rms,
            "full_reconstructed_force_vector_rms_with_radial_term": full_stress_force_rms,
            "radial_force_over_tangential_force_rms": radial_force_rms / scale,
            "radial_force_over_raw_ring_radial_defect_rms": radial_force_rms
            / max(raw_radial_rms, np.finfo(float).tiny),
            "radial_force_over_gated_ring_radial_defect_rms": radial_force_rms
            / max(gated_radial_rms, np.finfo(float).tiny),
            "radial_force_machine_nonzero": radial_force_nonzero,
        },
        "routing": {
            "two_tangential_channel_stress_is_full_tensor_force": False,
            "radial_force_side_effect_accounted": True,
            "radial_force_side_effect_numerically_stable": derivative_audit[
                "radial_force_derivative_stability_preflight_passed"
            ],
            "second_public_covariance_column_available": False,
            "finite_correction_cycle_rerun_allowed": False,
            "next_required": (
                "When Agent 2 supplies a genuinely independent public second velocity "
                "column, keep the existing theta/z bounded-inverse screen but also carry "
                "the source-required radial d_z sigma_1 channel into the materialized "
                "correction and judge the complete held-in/held-out NS residual. A "
                "two-channel covariance fit alone must not authorize the cycle."
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_consumed": True,
            "surrogate_defect_used": False,
            "source_radial_force_channel_measured": True,
            "new_oscillatory_column_constructed": False,
            "agent2_complete_curl_reimplemented": False,
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
    path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/radial_force_completeness_audit_report.json",
    )
    args = parser.parse_args()
    report = generate_actual_core_report(output=args.output)
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
