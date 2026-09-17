"""Measured signed-covariance inverse for the Agent-3 correction lane.

This module is deliberately narrower than a new oscillatory construction.  It
consumes Agent-2's existing :class:`KokunoCompleteCurlCorrection`, Agent-3's
real phase-mean defect contract, and Agent-3's compact radial-stress inverse.
It then asks the next source-motivated question: can the measured theta/radial
stress be represented by a *signed differential amplitude* of the already
routed oscillatory column?

The corrected Kokuno reader (2026-09-09) does not take a square root of a
perturbed covariance.  With its two fixed positive primary amplitudes it uses

    d_Sigma = H^{-1}(Sigma / epsilon),
    delta a_sigma = (d_Sigma)_sigma / (2 a_sigma),
    L Sigma = sqrt(epsilon) sum_sigma delta a_sigma b_sigma,

and proves the linear covariance identity

    B(W0, L Sigma) = Sigma.

The current repository does not expose Kokuno's source ``H``, ``a_sigma`` or
its two source pulse columns.  We therefore do **not** claim to implement that
map exactly.  Instead we measure, from the public Agent-2 velocity itself, the
rank-one physical response

    d/da <w_r w_theta> |_{a=a0} = 2 a0 <w_unit,r w_unit,theta>,

and solve the scalar analogue pointwise on the already reconstructed compact
radial stress.  This is an executable realizability/conditioning diagnostic,
not yet a public velocity correction.  A failed or ill-conditioned result is
retained as a routing result; the amplitude range is never widened to force a
success.

No pressure or force is fitted, no hand-entered residual is accepted, and no
PDE threshold is changed.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import argparse
import json
from math import pi
from pathlib import Path
from typing import Any

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_radial_stress import (
    CompactRadialStressInverse,
    sample_real_phase_mean_radial_profiles,
)


TASK = "KOKUNO-A3-SIGNED-COVARIANCE-INVERSE-005"
SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_SECTION = "Signed differential inverse and assembly across bands"
ROUTED_OSCILLATORY_AMPLITUDE = 0.125
ROUTED_OSCILLATORY_PHASE = 0.0
PROFILE_TIME = 0.5
PROFILE_Z = 0.0
PROFILE_ANNULUS = (0.08, 0.36)
DEFAULT_RADIAL_COUNT = 33
DEFAULT_ANGULAR_COUNT = 8
DEFAULT_PHASE_COUNT = 8
TARGET_FLOOR_FRACTION = 1.0e-3
RESPONSE_FLOOR_FRACTION = 1.0e-3
EXISTING_AGENT2_AMPLITUDE_BOUND = 2.0


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _wrap_phase(value: float) -> float:
    wrapped = (float(value) + pi) % (2.0 * pi) - pi
    if wrapped == -pi and value > 0.0:
        return pi
    return wrapped


def _ring_points(radii: np.ndarray, *, angular_count: int, z: float) -> tuple[np.ndarray, np.ndarray]:
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
    points = np.stack(
        (rr * np.cos(aa), rr * np.sin(aa), np.full_like(rr, float(z))), axis=-1
    )
    return points, aa


def _cylindrical_components(values: np.ndarray, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    utheta = (-points[..., 1] * values[..., 0] + points[..., 0] * values[..., 1]) / radius
    uz = values[..., 2]
    return ur, utheta, uz


def measure_unit_radial_theta_covariance(
    correction: KokunoCompleteCurlCorrection,
    radii: np.ndarray,
    *,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> np.ndarray:
    """Measure ``<w_r w_theta>`` for amplitude one.

    The average is over the repository's existing uniform phase adapter and a
    uniform physical-angle ring.  The complete-curl velocity is always obtained
    from Agent 2's public implementation; this module does not reimplement it.
    """
    if not isinstance(correction, KokunoCompleteCurlCorrection):
        raise TypeError("correction must be KokunoCompleteCurlCorrection")
    if phase_count < 4 or phase_count % 2:
        raise ValueError("phase_count must be an even integer >= 4")
    if not np.isfinite(time):
        raise ValueError("time must be finite")
    points, _ = _ring_points(radii, angular_count=angular_count, z=z)
    flat = points.reshape(-1, 3)
    unit = replace(correction, amplitude=1.0)
    accumulator = np.zeros(len(radii), dtype=float)
    shifts = 2.0 * pi * np.arange(phase_count, dtype=float) / phase_count
    for shift in shifts:
        shifted = replace(unit, phase=_wrap_phase(unit.phase + float(shift)))
        values = np.asarray(shifted.at_points(flat, float(time)), dtype=float).reshape(points.shape)
        ur, utheta, _ = _cylindrical_components(values, points)
        accumulator += np.mean(ur * utheta, axis=1)
    return accumulator / phase_count


def finite_difference_amplitude_derivative_error(
    correction: KokunoCompleteCurlCorrection,
    radii: np.ndarray,
    *,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
    amplitude_step: float = 1.0e-4,
) -> dict[str, float]:
    """Calibrate the measured covariance derivative against an FD amplitude probe."""
    if not 0.0 < amplitude_step < 0.01:
        raise ValueError("amplitude_step must lie in (0,0.01)")
    a0 = float(correction.amplitude)
    if abs(a0) <= amplitude_step:
        raise ValueError("base amplitude must stay away from zero")
    unit_cov = measure_unit_radial_theta_covariance(
        correction,
        radii,
        time=time,
        z=z,
        angular_count=angular_count,
        phase_count=phase_count,
    )

    def covariance_at(amplitude: float) -> np.ndarray:
        probe = replace(correction, amplitude=float(amplitude))
        # Velocity is exactly linear in amplitude, so reuse the measured unit
        # covariance rather than introducing a second curl implementation.
        return amplitude * amplitude * unit_cov

    fd = (covariance_at(a0 + amplitude_step) - covariance_at(a0 - amplitude_step)) / (2.0 * amplitude_step)
    analytic = 2.0 * a0 * unit_cov
    error = fd - analytic
    scale = max(_rms(analytic), np.finfo(float).tiny)
    return {
        "analytic_derivative_rms": _rms(analytic),
        "fd_derivative_rms": _rms(fd),
        "derivative_error_rms": _rms(error),
        "derivative_error_relative_rms": _rms(error) / scale,
        "amplitude_step": float(amplitude_step),
    }


@dataclass(frozen=True)
class SignedAmplitudeProfile:
    """Sampled scalar analogue of Kokuno's signed differential amplitude map."""

    radii: np.ndarray
    target_stress: np.ndarray
    covariance_derivative: np.ndarray
    delta_amplitude: np.ndarray
    predicted_linear_covariance: np.ndarray
    active_target_mask: np.ndarray
    safe_response_mask: np.ndarray
    base_amplitude: float
    target_floor_fraction: float = TARGET_FLOOR_FRACTION
    response_floor_fraction: float = RESPONSE_FLOOR_FRACTION

    def __post_init__(self) -> None:
        arrays = {
            "radii": np.asarray(self.radii, dtype=float),
            "target_stress": np.asarray(self.target_stress, dtype=float),
            "covariance_derivative": np.asarray(self.covariance_derivative, dtype=float),
            "delta_amplitude": np.asarray(self.delta_amplitude, dtype=float),
            "predicted_linear_covariance": np.asarray(self.predicted_linear_covariance, dtype=float),
            "active_target_mask": np.asarray(self.active_target_mask, dtype=bool),
            "safe_response_mask": np.asarray(self.safe_response_mask, dtype=bool),
        }
        shape = arrays["radii"].shape
        if len(shape) != 1 or len(arrays["radii"]) < 9:
            raise ValueError("profile arrays must be one-dimensional with at least 9 nodes")
        for name, value in arrays.items():
            if value.shape != shape:
                raise ValueError(f"{name} must match radii shape")
            if value.dtype != bool and not np.isfinite(value).all():
                raise ValueError(f"{name} must be finite")
        if np.any(np.diff(arrays["radii"]) <= 0.0):
            raise ValueError("radii must be strictly increasing")
        if not np.isfinite(self.base_amplitude) or self.base_amplitude <= 0.0:
            raise ValueError("base_amplitude must be positive and finite")
        if not 0.0 < self.target_floor_fraction < 0.1:
            raise ValueError("target_floor_fraction must lie in (0,0.1)")
        if not 0.0 < self.response_floor_fraction < 0.1:
            raise ValueError("response_floor_fraction must lie in (0,0.1)")
        for name, value in arrays.items():
            copied = value.copy()
            copied.setflags(write=False)
            object.__setattr__(self, name, copied)

    @property
    def unresolved_mask(self) -> np.ndarray:
        return self.active_target_mask & ~self.safe_response_mask

    @property
    def rank1_realizable_on_active_nodes(self) -> bool:
        return not bool(np.any(self.unresolved_mask))

    def metrics(self) -> dict[str, Any]:
        active = self.active_target_mask
        safe_active = active & self.safe_response_mask
        unresolved = self.unresolved_mask
        error = self.predicted_linear_covariance - self.target_stress
        target_scale = max(_rms(self.target_stress[active]) if np.any(active) else 0.0, np.finfo(float).tiny)
        response_values = np.abs(self.covariance_derivative[safe_active])
        if len(response_values):
            response_condition = float(np.max(response_values) / max(np.min(response_values), np.finfo(float).tiny))
        else:
            response_condition = float("inf")
        if len(self.radii) >= 5:
            first = np.gradient(self.delta_amplitude, self.radii, edge_order=2)
            second = np.gradient(first, self.radii, edge_order=2)
            max_slope = float(np.max(np.abs(first)))
            max_curvature = float(np.max(np.abs(second)))
        else:
            max_slope = float("nan")
            max_curvature = float("nan")
        return {
            "active_target_nodes": int(np.count_nonzero(active)),
            "safe_active_nodes": int(np.count_nonzero(safe_active)),
            "unresolved_active_nodes": int(np.count_nonzero(unresolved)),
            "rank1_realizable_on_active_nodes": self.rank1_realizable_on_active_nodes,
            "target_stress_rms": _rms(self.target_stress),
            "target_stress_max_abs": float(np.max(np.abs(self.target_stress))),
            "covariance_derivative_rms": _rms(self.covariance_derivative),
            "covariance_derivative_max_abs": float(np.max(np.abs(self.covariance_derivative))),
            "safe_response_condition_ratio": response_condition,
            "delta_amplitude_rms": _rms(self.delta_amplitude),
            "delta_amplitude_max_abs": float(np.max(np.abs(self.delta_amplitude))),
            "total_local_amplitude_max_abs": float(np.max(np.abs(self.base_amplitude + self.delta_amplitude))),
            "within_existing_agent2_absolute_amplitude_bound": bool(
                np.max(np.abs(self.base_amplitude + self.delta_amplitude)) <= EXISTING_AGENT2_AMPLITUDE_BOUND
            ),
            "linear_reconstruction_rms": _rms(error[active]) if np.any(active) else 0.0,
            "linear_reconstruction_relative_rms": (_rms(error[active]) / target_scale) if np.any(active) else 0.0,
            "delta_amplitude_edge_max_abs": float(
                max(abs(float(self.delta_amplitude[0])), abs(float(self.delta_amplitude[-1])))
            ),
            "sampled_delta_amplitude_max_abs_radial_slope": max_slope,
            "sampled_delta_amplitude_max_abs_radial_second_derivative": max_curvature,
        }


def build_rank1_signed_amplitude_profile(
    radii: np.ndarray,
    target_stress: np.ndarray,
    unit_covariance: np.ndarray,
    *,
    base_amplitude: float,
    target_floor_fraction: float = TARGET_FLOOR_FRACTION,
    response_floor_fraction: float = RESPONSE_FLOOR_FRACTION,
) -> SignedAmplitudeProfile:
    """Solve the scalar differential-covariance analogue without hiding rank loss."""
    radii = np.asarray(radii, dtype=float)
    target = np.asarray(target_stress, dtype=float)
    unit = np.asarray(unit_covariance, dtype=float)
    if radii.ndim != 1 or target.shape != radii.shape or unit.shape != radii.shape:
        raise ValueError("radii, target_stress and unit_covariance must be matching 1-D arrays")
    if not np.isfinite(radii).all() or not np.isfinite(target).all() or not np.isfinite(unit).all():
        raise ValueError("profile inputs must be finite")
    if not np.isfinite(base_amplitude) or base_amplitude <= 0.0:
        raise ValueError("base_amplitude must be positive and finite")
    if not 0.0 < target_floor_fraction < 0.1 or not 0.0 < response_floor_fraction < 0.1:
        raise ValueError("floor fractions must lie in (0,0.1)")

    derivative = 2.0 * float(base_amplitude) * unit
    target_max = float(np.max(np.abs(target)))
    if target_max <= np.finfo(float).tiny:
        raise ValueError("target stress is numerically zero; no signed update is needed")
    active = np.abs(target) >= target_floor_fraction * target_max
    # Exact compact-support endpoints are not used to diagnose conditioning.
    active[0] = False
    active[-1] = False
    response_scale = float(np.max(np.abs(derivative[active]))) if np.any(active) else 0.0
    if response_scale <= np.finfo(float).tiny:
        safe = np.zeros_like(active, dtype=bool)
    else:
        safe = np.abs(derivative) >= response_floor_fraction * response_scale

    delta = np.zeros_like(target)
    solvable = active & safe
    delta[solvable] = target[solvable] / derivative[solvable]
    predicted = derivative * delta
    # Preserve the exact compact edge values from the radial inverse.
    delta[0] = 0.0
    delta[-1] = 0.0
    predicted[0] = 0.0
    predicted[-1] = 0.0

    return SignedAmplitudeProfile(
        radii=radii,
        target_stress=target,
        covariance_derivative=derivative,
        delta_amplitude=delta,
        predicted_linear_covariance=predicted,
        active_target_mask=active,
        safe_response_mask=safe,
        base_amplitude=float(base_amplitude),
        target_floor_fraction=float(target_floor_fraction),
        response_floor_fraction=float(response_floor_fraction),
    )


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/signed_covariance_inverse_report.json",
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> dict[str, Any]:
    """Drive the signed-amplitude diagnostic from the real leading-core defect."""
    leading = KokunoLeadingCoreSeriesCandidate()
    correction = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=phase_count)
    theta_profile, axial_profile, projection_metrics = sample_real_phase_mean_radial_profiles(
        leading.at_points,
        correction,
        contract,
        time=PROFILE_TIME,
        z=PROFILE_Z,
        r_inner=PROFILE_ANNULUS[0],
        r_outer=PROFILE_ANNULUS[1],
        radial_count=radial_count,
        angular_count=angular_count,
    )
    theta_inverse = CompactRadialStressInverse(theta_profile)
    target_stress = theta_inverse.stress(theta_profile.radii)
    unit_covariance = measure_unit_radial_theta_covariance(
        correction,
        theta_profile.radii,
        time=PROFILE_TIME,
        z=PROFILE_Z,
        angular_count=angular_count,
        phase_count=phase_count,
    )
    signed_profile = build_rank1_signed_amplitude_profile(
        theta_profile.radii,
        target_stress,
        unit_covariance,
        base_amplitude=correction.amplitude,
    )
    signed_metrics = signed_profile.metrics()
    derivative_calibration = finite_difference_amplitude_derivative_error(
        correction,
        theta_profile.radii,
        time=PROFILE_TIME,
        z=PROFILE_Z,
        angular_count=angular_count,
        phase_count=phase_count,
    )
    radial_metrics = theta_inverse.metrics(query_count=129)

    report: dict[str, Any] = {
        "task": TASK,
        "source": {
            "reader": SOURCE_READER,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "repository": SOURCE_REPOSITORY,
            "commit": SOURCE_COMMIT,
            "path": SOURCE_PATH,
            "section": SOURCE_SECTION,
            "exact_source_linear_map": (
                "d_Sigma=H^{-1}(Sigma/epsilon); delta_a_sigma=(d_Sigma)_sigma/(2*a_sigma); "
                "B(W0,L Sigma)=Sigma"
            ),
            "exact_source_nonlinear_boundary": (
                "C(W0+L Sigma)=C(W0)+Sigma+C(L Sigma); no square root of a perturbed covariance"
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading.sha256,
            "oscillatory_amplitude": correction.amplitude,
            "oscillatory_phase": correction.phase,
            "profile_time": PROFILE_TIME,
            "profile_z": PROFILE_Z,
            "profile_annulus": list(PROFILE_ANNULUS),
            "radial_count": radial_count,
            "angular_count": angular_count,
            "phase_count": phase_count,
            "nu": contract.nu,
            "derivative_step": contract.step,
            "defect_kind": theta_profile.source_kind,
            "surrogate_defect_used": False,
        },
        "real_defect_projection": projection_metrics,
        "theta_radial_inverse": radial_metrics,
        "axial_projection_rms": _rms(axial_profile.values),
        "measured_rank1_covariance": {
            "unit_radial_theta_covariance_rms": _rms(unit_covariance),
            "unit_radial_theta_covariance_max_abs": float(np.max(np.abs(unit_covariance))),
            **derivative_calibration,
            **signed_metrics,
        },
        "routing": {
            "source_exact_two_column_map_available": False,
            "current_adapter": "measured one-column scalar differential covariance response",
            "rank1_realization_ready_for_public_velocity_lift": bool(
                signed_metrics["rank1_realizable_on_active_nodes"]
                and signed_metrics["within_existing_agent2_absolute_amplitude_bound"]
            ),
            "next_required_if_rank1_fails": (
                "expose/construct a second independent fixed positive oscillatory column or source H-map; "
                "do not widen amplitudes or reuse residual samples"
            ),
            "next_required_if_rank1_passes": (
                "materialize the sampled signed amplitude at vector-potential level with the complete curl, "
                "then recompute held-in/held-out real residual before any acceptance"
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_consumed": True,
            "radial_stress_reconstructed_from_real_defect": True,
            "source_signed_differential_formula_recorded": True,
            "source_exact_signed_map_implemented": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_advanced_to_covariance_realizability": True,
            "residual_reduction_claimed": False,
            "formal_full_domain_normalized_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Measure whether the real compact theta stress is realizable by the routed Agent-2 column's signed covariance derivative"
    )
    parser.add_argument("--output", default="artifacts/kokuno_agent3/signed_covariance_inverse_report.json")
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--phase-count", type=int, default=DEFAULT_PHASE_COUNT)
    args = parser.parse_args()
    report = generate_actual_core_report(
        output=args.output,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        phase_count=args.phase_count,
    )
    print(json.dumps({
        "task": report["task"],
        "rank1_realization_ready": report["routing"]["rank1_realization_ready_for_public_velocity_lift"],
        "target_stress_rms": report["measured_rank1_covariance"]["target_stress_rms"],
        "covariance_derivative_rms": report["measured_rank1_covariance"]["covariance_derivative_rms"],
        "unresolved_active_nodes": report["measured_rank1_covariance"]["unresolved_active_nodes"],
        "delta_amplitude_max_abs": report["measured_rank1_covariance"]["delta_amplitude_max_abs"],
    }, indent=2))


if __name__ == "__main__":
    main()
