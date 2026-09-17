"""Audit whether Agent-2 phase copies can supply Kokuno's missing second covariance column.

The corrected Kokuno reader uses two *distinct fixed positive* pulse columns in the
signed differential inverse.  In its local (N,K) frame,

    H_sigma = h_sigma (-A_c N - sigma u_* K + e_sigma),

so, before the small column errors, the reference response matrix is

    H_ref = [[-A_c h_+, -A_c h_-],
             [-u_* h_+,  u_* h_-]],

with determinant ``-2 A_c u_* h_+ h_-``.  Thus the source inverse depends on a
genuinely rank-two covariance map.  PR #253 measured only one physical response
column from Agent 2's already-routed complete curl, and PR #263 showed that the
resulting one-column signed correction does not reduce the real held-out mean
defect after public complete-curl materialization.

This module tests one tempting but invalid shortcut without constructing a new
oscillatory field: can two phase offsets of the *existing* Agent-2 column act as
the missing source columns after the same phase averaging used by the mean-defect
contract?  It consumes the real Agent-1 + Agent-2 phase-mean defect, reconstructs
both compact theta/e=2 and axial/e=1 radial stresses, and measures the two-component
covariance response ``(<w_r w_theta>, <w_r w_z>)`` only through Agent 2's public
``at_points`` velocity.

Uniform phase averaging makes a phase translate of one periodic column the same
quadratic covariance column.  The report therefore fails closed if the measured
phase-orbit matrix is not rank two.  No polarization, wave vector, amplitude,
pressure, force, residual threshold, or velocity correction is selected here.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import argparse
import json
from math import pi
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_defect import PhaseMeanDefectContract
from .kokuno_radial_stress import CompactRadialStressInverse, sample_real_phase_mean_radial_profiles
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    PROFILE_TIME,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
)


TASK = "KOKUNO-A3-PHASE-ORBIT-COVARIANCE-RANK-007"
SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_SECTION = "Signed covariance and exact signed differential inverse"
PROFILE_Z = 0.08
PHASE_OFFSETS = (0.0, 0.5 * pi)
DEFAULT_RADIAL_COUNT = 33
DEFAULT_ANGULAR_COUNT = 8
DEFAULT_PHASE_COUNT = 8
RANK_RELATIVE_TOLERANCE = 1.0e-10


def _rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _wrap_phase(value: float) -> float:
    wrapped = (float(value) + pi) % (2.0 * pi) - pi
    if wrapped == -pi and value > 0.0:
        return pi
    return wrapped


def source_reference_response_matrix(
    *,
    A_c: float,
    u_star: float,
    h_plus: float,
    h_minus: float,
) -> np.ndarray:
    """Return the source's error-free two-column reference matrix in (N,K)."""
    values = np.asarray([A_c, u_star, h_plus, h_minus], dtype=float)
    if not np.isfinite(values).all() or np.any(values <= 0.0):
        raise ValueError("A_c, u_star, h_plus and h_minus must be positive and finite")
    return np.asarray(
        [
            [-A_c * h_plus, -A_c * h_minus],
            [-u_star * h_plus, u_star * h_minus],
        ],
        dtype=float,
    )


def source_reference_determinant(
    *,
    A_c: float,
    u_star: float,
    h_plus: float,
    h_minus: float,
) -> float:
    """Exact determinant ``-2 A_c u_* h_+ h_-`` from the corrected reader."""
    source_reference_response_matrix(
        A_c=A_c, u_star=u_star, h_plus=h_plus, h_minus=h_minus
    )
    return float(-2.0 * A_c * u_star * h_plus * h_minus)


def _ring_points(
    radii: np.ndarray,
    *,
    angular_count: int,
    z: float,
) -> np.ndarray:
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


def _cylindrical_velocity(values: np.ndarray, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    return ur, utheta, values[..., 2]


def measure_phase_mean_covariance_vector(
    correction: KokunoCompleteCurlCorrection,
    radii: np.ndarray,
    *,
    phase_offset: float,
    time: float = PROFILE_TIME,
    z: float = PROFILE_Z,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
) -> np.ndarray:
    """Return ``(<w_r w_theta>, <w_r w_z>)`` for an amplitude-one phase copy.

    Both physical-angle and phase averages use public Agent-2 velocity values.
    The caller-supplied ``phase_offset`` only translates the same existing
    column; this function never changes its polarization, wave vector or support.
    """
    if not isinstance(correction, KokunoCompleteCurlCorrection):
        raise TypeError("correction must be KokunoCompleteCurlCorrection")
    if phase_count < 4 or phase_count % 2:
        raise ValueError("phase_count must be an even integer >= 4")
    if not np.isfinite([phase_offset, time, z]).all():
        raise ValueError("phase_offset, time and z must be finite")

    points = _ring_points(radii, angular_count=angular_count, z=z)
    flat = points.reshape(-1, 3)
    unit = replace(correction, amplitude=1.0)
    accumulator = np.zeros((len(radii), 2), dtype=float)
    shifts = 2.0 * pi * np.arange(phase_count, dtype=float) / phase_count
    for shift in shifts:
        probe = replace(
            unit,
            phase=_wrap_phase(unit.phase + float(phase_offset) + float(shift)),
        )
        velocity = np.asarray(probe.at_points(flat, float(time)), dtype=float).reshape(points.shape)
        ur, utheta, uz = _cylindrical_velocity(velocity, points)
        accumulator[:, 0] += np.mean(ur * utheta, axis=1)
        accumulator[:, 1] += np.mean(ur * uz, axis=1)
    return accumulator / float(phase_count)


@dataclass(frozen=True)
class ResponseRankAudit:
    singular_values: np.ndarray
    rank: int
    condition: float
    relative_target_projection_rms: float
    coefficients: np.ndarray

    def __post_init__(self) -> None:
        singular = np.asarray(self.singular_values, dtype=float)
        coefficients = np.asarray(self.coefficients, dtype=float)
        if singular.ndim != 1 or coefficients.ndim != 1:
            raise ValueError("singular_values and coefficients must be one-dimensional")
        if not np.isfinite(singular).all() or not np.isfinite(coefficients).all():
            raise ValueError("rank-audit arrays must be finite")
        if self.rank < 0 or self.rank > len(singular):
            raise ValueError("invalid response rank")
        if not np.isfinite(self.relative_target_projection_rms) or self.relative_target_projection_rms < 0.0:
            raise ValueError("projection error must be finite and nonnegative")
        if not (np.isfinite(self.condition) or np.isinf(self.condition)):
            raise ValueError("condition must be finite or infinity")
        singular = singular.copy()
        coefficients = coefficients.copy()
        singular.setflags(write=False)
        coefficients.setflags(write=False)
        object.__setattr__(self, "singular_values", singular)
        object.__setattr__(self, "coefficients", coefficients)


def analyze_response_matrix(
    matrix: np.ndarray,
    target: np.ndarray,
    *,
    relative_rank_tolerance: float = RANK_RELATIVE_TOLERANCE,
) -> ResponseRankAudit:
    """Measure numerical rank and least-squares target coverage without selection."""
    matrix = np.asarray(matrix, dtype=float)
    target = np.asarray(target, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] < 1 or target.shape != (matrix.shape[0],):
        raise ValueError("matrix must be (M,K) and target must be (M,)")
    if not np.isfinite(matrix).all() or not np.isfinite(target).all():
        raise ValueError("matrix and target must be finite")
    if not 0.0 < relative_rank_tolerance < 1.0e-2:
        raise ValueError("relative_rank_tolerance must lie in (0,1e-2)")
    singular = np.linalg.svd(matrix, compute_uv=False)
    if len(singular) == 0 or singular[0] <= np.finfo(float).tiny:
        raise ValueError("response matrix is numerically zero")
    rank = int(np.count_nonzero(singular > relative_rank_tolerance * singular[0]))
    condition = float(singular[0] / singular[-1]) if rank == matrix.shape[1] else float("inf")
    coefficients, _, _, _ = np.linalg.lstsq(matrix, target, rcond=None)
    residual = matrix @ coefficients - target
    target_scale = max(_rms(target), np.finfo(float).tiny)
    return ResponseRankAudit(
        singular_values=singular,
        rank=rank,
        condition=condition,
        relative_target_projection_rms=_rms(residual) / target_scale,
        coefficients=coefficients,
    )


def _stack_response(response: np.ndarray) -> np.ndarray:
    response = np.asarray(response, dtype=float)
    if response.ndim != 2 or response.shape[1] != 2 or len(response) < 3:
        raise ValueError("response must have shape (radial_count,2)")
    # The compact stress is exactly zero on the two annulus edges, and those
    # endpoint rows carry no rank information.
    return response[1:-1].reshape(-1)


def generate_actual_core_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/phase_orbit_covariance_rank_report.json",
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
    phase_offsets: Iterable[float] = PHASE_OFFSETS,
) -> dict[str, Any]:
    """Audit the real two-component stress against phase copies of one curl column."""
    offsets = tuple(float(v) for v in phase_offsets)
    if len(offsets) != 2 or not np.isfinite(offsets).all():
        raise ValueError("exactly two finite phase offsets are required")

    leading = KokunoLeadingCoreSeriesCandidate()
    correction = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=phase_count)
    theta_profile, axial_profile, projection = sample_real_phase_mean_radial_profiles(
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
    axial_inverse = CompactRadialStressInverse(axial_profile)
    radii = theta_profile.radii
    target = np.stack(
        (theta_inverse.stress(radii), axial_inverse.stress(radii)), axis=-1
    )

    covariances = [
        measure_phase_mean_covariance_vector(
            correction,
            radii,
            phase_offset=offset,
            time=PROFILE_TIME,
            z=PROFILE_Z,
            angular_count=angular_count,
            phase_count=phase_count,
        )
        for offset in offsets
    ]
    response_columns = [2.0 * correction.amplitude * covariance for covariance in covariances]
    matrix = np.stack([_stack_response(column) for column in response_columns], axis=1)
    target_vector = _stack_response(target)
    audit = analyze_response_matrix(matrix, target_vector)
    column_scale = max(_rms(response_columns[0]), np.finfo(float).tiny)
    phase_orbit_difference = _rms(response_columns[1] - response_columns[0]) / column_scale

    source_demo = source_reference_response_matrix(
        A_c=1.0, u_star=1.0, h_plus=1.0, h_minus=1.0
    )
    source_demo_singular = np.linalg.svd(source_demo, compute_uv=False)

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
            "signed_inverse": (
                "d_Sigma=H^{-1}(Sigma/epsilon); delta_a_sigma=(d_Sigma)_sigma/(2*a_sigma); "
                "B(W0,L Sigma)=Sigma"
            ),
            "reference_columns": "H_sigma=h_sigma*(-A_c*N-sigma*u_*K+e_sigma)",
            "reference_determinant": "det(H_ref)=-2*A_c*u_*h_+*h_-",
            "nonlinear_remainder_retained": "C(W0+L Sigma)=epsilon*T0+Sigma+C(L Sigma)",
        },
        "inputs": {
            "leading_candidate_sha256": leading.sha256,
            "oscillatory_amplitude": float(correction.amplitude),
            "oscillatory_phase": float(correction.phase),
            "profile_time": PROFILE_TIME,
            "profile_z": PROFILE_Z,
            "profile_annulus": list(PROFILE_ANNULUS),
            "radial_count": radial_count,
            "angular_count": angular_count,
            "phase_count": phase_count,
            "phase_offsets": list(offsets),
            "nu": contract.nu,
            "derivative_step": contract.step,
            "surrogate_defect_used": False,
            "new_oscillatory_column_constructed": False,
        },
        "real_defect_projection": projection,
        "real_compact_stress": {
            "theta_stress_rms": _rms(target[:, 0]),
            "theta_stress_max_abs": float(np.max(np.abs(target[:, 0]))),
            "axial_stress_rms": _rms(target[:, 1]),
            "axial_stress_max_abs": float(np.max(np.abs(target[:, 1]))),
            "two_component_stress_rms": _rms(target),
            "theta_radial_identity_fd_rms": theta_inverse.metrics(query_count=129)["radial_identity_fd_rms"],
            "axial_radial_identity_fd_rms": axial_inverse.metrics(query_count=129)["radial_identity_fd_rms"],
        },
        "source_rank2_reference_calibration": {
            "unit_parameter_matrix": source_demo.tolist(),
            "determinant": float(np.linalg.det(source_demo)),
            "singular_values": [float(v) for v in source_demo_singular],
            "rank": int(np.linalg.matrix_rank(source_demo)),
            "used_as_candidate_data": False,
        },
        "existing_agent2_phase_orbit": {
            "response_column_0_rms": _rms(response_columns[0]),
            "response_column_1_rms": _rms(response_columns[1]),
            "phase_orbit_column_relative_difference_rms": phase_orbit_difference,
            "singular_values": [float(v) for v in audit.singular_values],
            "rank": int(audit.rank),
            "condition": None if np.isinf(audit.condition) else float(audit.condition),
            "condition_is_infinite": bool(np.isinf(audit.condition)),
            "real_target_projection_relative_rms": float(audit.relative_target_projection_rms),
            "least_squares_coefficients": [float(v) for v in audit.coefficients],
            "rank2_required_by_source_signed_inverse": True,
            "phase_orbit_supplies_second_independent_column": bool(audit.rank >= 2),
        },
        "routing": {
            "phase_copy_accepted_as_source_second_column": False,
            "second_independent_covariance_column_available": bool(audit.rank >= 2),
            "finite_correction_cycle_rerun_allowed_from_this_increment": False,
            "next_required": (
                "consume a genuinely independent source-motivated Agent-2 pulse/covariance column "
                "after it exists; do not use phase translation of the current column, widen the "
                "rejected signed amplitude, or rerun a finite correction cycle from this rank audit"
            ),
        },
        "truth_boundary": {
            "real_candidate_defect_consumed": True,
            "theta_and_axial_radial_stress_reconstructed": True,
            "source_rank2_signed_map_requirement_recorded": True,
            "agent2_complete_curl_reimplemented": False,
            "new_oscillatory_basis_or_polarization_selected": False,
            "public_velocity_correction_materialized": False,
            "finite_correction_cycle_run": False,
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
        description="Audit whether phase copies of the routed Agent-2 curl provide Kokuno's required second covariance column"
    )
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/phase_orbit_covariance_rank_report.json",
    )
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
    result = report["existing_agent2_phase_orbit"]
    print(
        json.dumps(
            {
                "task": report["task"],
                "rank": result["rank"],
                "singular_values": result["singular_values"],
                "phase_orbit_column_relative_difference_rms": result[
                    "phase_orbit_column_relative_difference_rms"
                ],
                "real_target_projection_relative_rms": result[
                    "real_target_projection_relative_rms"
                ],
                "second_independent_covariance_column_available": report["routing"][
                    "second_independent_covariance_column_available"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
