"""Materialize Agent-3's signed covariance profile through Agent-2's complete curl.

This is one finite correction-cycle increment downstream of
``kokuno_signed_covariance_inverse``.  The corrected Kokuno reader uses signed
*differential amplitudes* to realize a compact covariance/stress target.  PR
#253 already measured the rank-one physical analogue from the real leading-core
mean defect and found a small, well-conditioned signed amplitude profile.

The missing operation is a public velocity.  We deliberately do not reimplement
Agent 2's curl kernel.  Let ``A_unit`` and ``w_unit=curl(A_unit)`` be Agent 2's
existing public vector potential and complete curl at amplitude one.  For a
compact radial signed amplitude ``alpha(r)`` we materialize

    A_delta = alpha(r) A_unit,
    delta_u = curl(A_delta)
            = alpha w_unit + grad(alpha) x A_unit.

The second term is mandatory: multiplying an already-computed velocity by a
spatial amplitude would drop the cutoff/amplitude-gradient curl remainder and
would not be divergence free.  The radial amplitude is a fixed-degree compact
C4 Chebyshev adapter to the sampled signed profile from PR #253.  Its degree is
fixed before the residual screen and is not tuned on held-in or held-out data.

The report then recomputes the real phase-mean nonlinear NS defect with the
independent CR006 finite-difference operator on fresh held-in and disjoint
held-out samples.  No pressure or force is fitted, no Agent-2 phase/amplitude is
optimized, and no project threshold is changed.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import argparse
import json
from math import pi
from pathlib import Path
from typing import Any, Callable

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval, chebvander

from .constrained_validation import residual
from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from .kokuno_mean_defect import (
    PhaseMeanDefectContract,
    cylindrical_theta_component,
)
from .kokuno_radial_stress import (
    CompactRadialStressInverse,
    sample_real_phase_mean_radial_profiles,
)
from .kokuno_signed_covariance_inverse import (
    PROFILE_ANNULUS,
    PROFILE_TIME,
    PROFILE_Z,
    ROUTED_OSCILLATORY_AMPLITUDE,
    ROUTED_OSCILLATORY_PHASE,
    SignedAmplitudeProfile,
    build_rank1_signed_amplitude_profile,
    measure_unit_radial_theta_covariance,
)


TASK = "KOKUNO-A3-SIGNED-AMPLITUDE-CURL-CYCLE-006"
SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_REPOSITORY = "KokunoYumeto/yang-mills-interacting-workbench"
SOURCE_COMMIT = "143f6773feb424ad9ed3a8d116653200f20346b7"
SOURCE_PATH = "navier-stokes/navier_stokes_workbench.tex"
SOURCE_SECTION = "Signed differential inverse and finite correction cycle"
CHEBYSHEV_DEGREE = 10
PROFILE_FIT_RELATIVE_RMS_LIMIT = 0.10
TRAIN_SEED = 9_173_021
HELD_OUT_SEED = 9_173_027
DIVERGENCE_SEED = 9_173_033
TRAIN_COUNT = 16
HELD_OUT_COUNT = 32
DIVERGENCE_COUNT = 96
TRAIN_TIME = 0.5
HELD_OUT_TIMES = (0.375, 0.5, 0.625)
CORE_SAMPLE_RADIAL_BOUNDS = (0.09, 0.22)
CORE_SAMPLE_Z_HALF_WIDTH = 0.12
DEFAULT_RADIAL_COUNT = 33
DEFAULT_ANGULAR_COUNT = 8
DEFAULT_PHASE_COUNT = 8
PointVelocity = Callable[[np.ndarray, float], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))


def _vector_max(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.max(np.linalg.norm(values, axis=-1)))


def _scalar_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(values * values)))


def _wrap_phase(value: float) -> float:
    wrapped = (float(value) + pi) % (2.0 * pi) - pi
    if wrapped == -pi and value > 0.0:
        return pi
    return wrapped


def _zero_pressure(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros(len(points), dtype=float)


def _zero_force(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros_like(points, dtype=float)


def _zero_velocity(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros_like(points, dtype=float)


@dataclass(frozen=True)
class CompactChebyshevSignedAmplitude:
    """C4 compact radial adapter for the sampled signed differential amplitude."""

    r_inner: float
    r_outer: float
    coefficients: tuple[float, ...]
    degree: int
    sampled_fit_rms: float
    sampled_fit_relative_rms: float
    design_condition: float

    def __post_init__(self) -> None:
        values = np.asarray(
            [self.r_inner, self.r_outer, self.sampled_fit_rms,
             self.sampled_fit_relative_rms, self.design_condition],
            dtype=float,
        )
        if not np.isfinite(values).all():
            raise ValueError("compact signed-amplitude metadata must be finite")
        if not 0.0 < self.r_inner < self.r_outer:
            raise ValueError("radial support must satisfy 0 < r_inner < r_outer")
        if self.degree < 0 or len(self.coefficients) != self.degree + 1:
            raise ValueError("coefficient count must equal degree + 1")
        coeff = np.asarray(self.coefficients, dtype=float)
        if not np.isfinite(coeff).all():
            raise ValueError("Chebyshev coefficients must be finite")
        if self.sampled_fit_rms < 0.0 or self.sampled_fit_relative_rms < 0.0:
            raise ValueError("fit errors must be nonnegative")
        if self.design_condition < 1.0:
            raise ValueError("design condition must be at least one")

    def _scaled(self, radius: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        radius = np.asarray(radius, dtype=float)
        if not np.isfinite(radius).all():
            raise ValueError("radius must be finite")
        scale = 2.0 / (self.r_outer - self.r_inner)
        s = scale * (radius - self.r_inner) - 1.0
        inside = (radius > self.r_inner) & (radius < self.r_outer)
        return s, inside, np.asarray(scale)

    def value_and_radial_derivative(self, radius: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        radius = np.asarray(radius, dtype=float)
        s, inside, scale = self._scaled(radius)
        base = np.where(inside, 1.0 - s * s, 0.0)
        gate = base**5
        gate_s = np.where(inside, -10.0 * s * base**4, 0.0)
        coeff = np.asarray(self.coefficients, dtype=float)
        poly = chebval(s, coeff)
        dcoeff = chebder(coeff)
        poly_s = chebval(s, dcoeff) if len(dcoeff) else np.zeros_like(s)
        value = np.where(inside, gate * poly, 0.0)
        derivative = np.where(
            inside,
            (gate_s * poly + gate * poly_s) * float(scale),
            0.0,
        )
        return value, derivative

    def __call__(self, radius: np.ndarray) -> np.ndarray:
        return self.value_and_radial_derivative(radius)[0]

    def metrics(self, radii: np.ndarray, sampled_target: np.ndarray) -> dict[str, float]:
        radii = np.asarray(radii, dtype=float)
        sampled_target = np.asarray(sampled_target, dtype=float)
        predicted, derivative = self.value_and_radial_derivative(radii)
        error = predicted - sampled_target
        scale = max(_scalar_rms(sampled_target), np.finfo(float).tiny)
        dense = np.linspace(self.r_inner, self.r_outer, 513)
        dense_value, dense_derivative = self.value_and_radial_derivative(dense)
        return {
            "degree": int(self.degree),
            "sampled_fit_rms": _scalar_rms(error),
            "sampled_fit_relative_rms": _scalar_rms(error) / scale,
            "design_condition": float(self.design_condition),
            "sampled_max_abs": float(np.max(np.abs(predicted))),
            "dense_max_abs": float(np.max(np.abs(dense_value))),
            "dense_max_abs_radial_derivative": float(np.max(np.abs(dense_derivative))),
            "inner_value": float(self(np.asarray([self.r_inner]))[0]),
            "outer_value": float(self(np.asarray([self.r_outer]))[0]),
            "inner_radial_derivative": float(
                self.value_and_radial_derivative(np.asarray([self.r_inner]))[1][0]
            ),
            "outer_radial_derivative": float(
                self.value_and_radial_derivative(np.asarray([self.r_outer]))[1][0]
            ),
        }


def fit_compact_chebyshev_signed_amplitude(
    sampled: SignedAmplitudeProfile,
    *,
    degree: int = CHEBYSHEV_DEGREE,
) -> CompactChebyshevSignedAmplitude:
    """Fit one fixed-degree C4 radial amplitude without looking at residual samples."""
    if not isinstance(sampled, SignedAmplitudeProfile):
        raise TypeError("sampled must be SignedAmplitudeProfile")
    if degree < 1 or degree >= len(sampled.radii) - 2:
        raise ValueError("degree must lie in [1, radial_count-3]")
    if not sampled.rank1_realizable_on_active_nodes:
        raise ValueError("cannot materialize a rank-deficient signed amplitude profile")

    radii = np.asarray(sampled.radii, dtype=float)
    target = np.asarray(sampled.delta_amplitude, dtype=float)
    r_inner = float(radii[0])
    r_outer = float(radii[-1])
    s = 2.0 * (radii - r_inner) / (r_outer - r_inner) - 1.0
    gate = np.maximum(1.0 - s * s, 0.0) ** 5
    design = gate[:, None] * chebvander(s, degree)
    coefficients, _, rank, singular = np.linalg.lstsq(design, target, rcond=None)
    if rank != degree + 1:
        raise RuntimeError("compact Chebyshev amplitude fit lost rank")
    positive = singular[singular > np.finfo(float).tiny]
    if len(positive) != degree + 1:
        raise RuntimeError("compact Chebyshev amplitude fit has singular design")
    condition = float(np.max(positive) / np.min(positive))
    predicted = design @ coefficients
    fit_rms = _scalar_rms(predicted - target)
    target_rms = max(_scalar_rms(target), np.finfo(float).tiny)
    return CompactChebyshevSignedAmplitude(
        r_inner=r_inner,
        r_outer=r_outer,
        coefficients=tuple(float(v) for v in coefficients),
        degree=int(degree),
        sampled_fit_rms=fit_rms,
        sampled_fit_relative_rms=fit_rms / target_rms,
        design_condition=condition,
    )


@dataclass(frozen=True)
class KokunoSignedAmplitudeCurlCorrection:
    """Signed radial amplitude materialized at Agent-2 vector-potential level."""

    primary: KokunoCompleteCurlCorrection
    radial_amplitude: CompactChebyshevSignedAmplitude

    def __post_init__(self) -> None:
        if not isinstance(self.primary, KokunoCompleteCurlCorrection):
            raise TypeError("primary must be KokunoCompleteCurlCorrection")
        if not isinstance(self.radial_amplitude, CompactChebyshevSignedAmplitude):
            raise TypeError("radial_amplitude must be CompactChebyshevSignedAmplitude")

    @property
    def unit(self) -> KokunoCompleteCurlCorrection:
        return replace(self.primary, amplitude=1.0)

    def with_phase(self, phase: float) -> "KokunoSignedAmplitudeCurlCorrection":
        return replace(self, primary=replace(self.primary, phase=_wrap_phase(phase)))

    def _amplitude_gradient(self, x, y, z, time) -> tuple[np.ndarray, np.ndarray]:
        x, y, z, time = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(time, dtype=float),
        )
        del time
        radius = np.hypot(x, y)
        amplitude, radial = self.radial_amplitude.value_and_radial_derivative(radius)
        gradient = np.zeros(amplitude.shape + (3,), dtype=float)
        mask = radius > 0.0
        gradient[..., 0] = np.where(mask, radial * x / np.where(mask, radius, 1.0), 0.0)
        gradient[..., 1] = np.where(mask, radial * y / np.where(mask, radius, 1.0), 0.0)
        return amplitude, gradient

    def vector_potential(self, x, y, z, time) -> np.ndarray:
        amplitude, _ = self._amplitude_gradient(x, y, z, time)
        return amplitude[..., None] * self.unit.vector_potential(x, y, z, time)

    def velocity(self, x, y, z, time) -> np.ndarray:
        amplitude, gradient = self._amplitude_gradient(x, y, z, time)
        unit = self.unit
        return (
            amplitude[..., None] * unit.velocity(x, y, z, time)
            + np.cross(gradient, unit.vector_potential(x, y, z, time))
        )

    __call__ = velocity

    def at_points(self, points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must be a finite (N,3) array")
        values = self.velocity(points[:, 0], points[:, 1], points[:, 2], time)
        if values.shape != points.shape or not np.isfinite(values).all():
            raise RuntimeError("signed complete-curl correction returned invalid velocity")
        return values

    def metadata(self) -> dict[str, Any]:
        return {
            "family": "kokuno_signed_radial_amplitude_complete_curl_adapter_v1",
            "source_structure": (
                "signed differential amplitude is materialized at vector-potential level; "
                "curl(alpha A)=alpha curl(A)+grad(alpha)xA"
            ),
            "primary_amplitude": float(self.primary.amplitude),
            "primary_phase": float(self.primary.phase),
            "radial_support": [
                float(self.radial_amplitude.r_inner),
                float(self.radial_amplitude.r_outer),
            ],
            "chebyshev_degree": int(self.radial_amplitude.degree),
            "source_exact_two_column_map": False,
            "agent2_curl_reimplemented": False,
        }


@dataclass(frozen=True)
class PhaseMeanCurlCycleSample:
    points: np.ndarray
    time: float
    base_momentum: np.ndarray
    baseline_phase_momentum: np.ndarray
    corrected_phase_momentum: np.ndarray
    corrected_phase_divergence: np.ndarray

    @property
    def before_increment(self) -> np.ndarray:
        return np.mean(self.baseline_phase_momentum, axis=0) - self.base_momentum

    @property
    def after_increment(self) -> np.ndarray:
        return np.mean(self.corrected_phase_momentum, axis=0) - self.base_momentum

    @property
    def corrected_mean_momentum(self) -> np.ndarray:
        return np.mean(self.corrected_phase_momentum, axis=0)


def _evaluate_phase_mean_curl_cycle(
    leading: KokunoLeadingCoreSeriesCandidate,
    primary: KokunoCompleteCurlCorrection,
    signed: KokunoSignedAmplitudeCurlCorrection,
    contract: PhaseMeanDefectContract,
    points: np.ndarray,
    time: float,
) -> PhaseMeanCurlCycleSample:
    baseline = contract.evaluate(leading.at_points, primary, points, time)
    phases = 2.0 * pi * np.arange(contract.phase_count, dtype=float) / contract.phase_count
    corrected_momentum: list[np.ndarray] = []
    corrected_divergence: list[np.ndarray] = []
    for shift in phases:
        phase = _wrap_phase(primary.phase + float(shift))
        shifted_primary = replace(primary, phase=phase)
        shifted_signed = signed.with_phase(phase)

        def composite_velocity(x: np.ndarray, t: float) -> np.ndarray:
            x = np.asarray(x, dtype=float)
            values = (
                np.asarray(leading.at_points(x, t), dtype=float)
                + np.asarray(shifted_primary.at_points(x, t), dtype=float)
                + np.asarray(shifted_signed.at_points(x, t), dtype=float)
            )
            if values.shape != x.shape or not np.isfinite(values).all():
                raise RuntimeError("corrected phase composite returned invalid velocity")
            return values

        result = residual(
            composite_velocity,
            _zero_pressure,
            _zero_force,
            points,
            float(time),
            nu=contract.nu,
            step=contract.step,
            time_bounds=contract.time_bounds,
        )
        corrected_momentum.append(np.asarray(result["momentum"], dtype=float))
        corrected_divergence.append(np.asarray(result["divergence"], dtype=float))

    return PhaseMeanCurlCycleSample(
        points=np.asarray(points, dtype=float).copy(),
        time=float(time),
        base_momentum=np.asarray(baseline.base_momentum, dtype=float),
        baseline_phase_momentum=np.asarray(baseline.phase_momentum, dtype=float),
        corrected_phase_momentum=np.stack(corrected_momentum, axis=0),
        corrected_phase_divergence=np.stack(corrected_divergence, axis=0),
    )


def _sample_core_points(*, seed: int, count: int) -> np.ndarray:
    if count <= 0:
        raise ValueError("count must be positive")
    r0, r1 = CORE_SAMPLE_RADIAL_BOUNDS
    rng = np.random.default_rng(seed)
    radius = np.sqrt(rng.uniform(r0 * r0, r1 * r1, count))
    angle = rng.uniform(0.0, 2.0 * pi, count)
    z = rng.uniform(-CORE_SAMPLE_Z_HALF_WIDTH, CORE_SAMPLE_Z_HALF_WIDTH, count)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def _assert_stencil_safe(
    leading: KokunoLeadingCoreSeriesCandidate,
    points: np.ndarray,
    times: tuple[float, ...],
    step: float,
) -> None:
    for time in times:
        for dt in (-step, 0.0, step):
            leading.coordinates(points[:, 0], points[:, 1], points[:, 2], time + dt)
        for axis in range(3):
            for multiple in (-2.0, -1.0, 1.0, 2.0):
                shifted = np.asarray(points, dtype=float).copy()
                shifted[:, axis] += multiple * step
                leading.coordinates(shifted[:, 0], shifted[:, 1], shifted[:, 2], time)


def _pressure_at_points(leading: KokunoLeadingCoreSeriesCandidate):
    def pressure(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        return np.asarray(
            leading.pressure(points[:, 0], points[:, 1], points[:, 2], time),
            dtype=float,
        )
    return pressure


def _fixed_pressure_gradient(
    leading: KokunoLeadingCoreSeriesCandidate,
    points: np.ndarray,
    time: float,
    contract: PhaseMeanDefectContract,
) -> np.ndarray:
    result = residual(
        _zero_velocity,
        _pressure_at_points(leading),
        _zero_force,
        points,
        float(time),
        nu=contract.nu,
        step=contract.step,
        time_bounds=contract.time_bounds,
    )
    return np.asarray(result["momentum"], dtype=float)


def _increment_metrics(values: np.ndarray, points: np.ndarray) -> dict[str, float]:
    theta = cylindrical_theta_component(values, points)
    return {
        "rms": _vector_rms(values),
        "max": _vector_max(values),
        "theta_rms": _scalar_rms(theta),
        "theta_max": float(np.max(np.abs(theta))),
    }


def _independent_divergence(
    correction: KokunoSignedAmplitudeCurlCorrection,
    points: np.ndarray,
    *,
    time: float = TRAIN_TIME,
    step: float = 1.0e-5,
) -> np.ndarray:
    divergence = np.zeros(len(points), dtype=float)
    for axis in range(3):
        offset = np.eye(3)[axis] * step
        plus = correction.at_points(points + offset, time)
        minus = correction.at_points(points - offset, time)
        divergence += ((plus - minus) / (2.0 * step))[:, axis]
    return divergence


def _build_actual_signed_correction(
    *,
    radial_count: int,
    angular_count: int,
    phase_count: int,
) -> tuple[
    KokunoLeadingCoreSeriesCandidate,
    KokunoCompleteCurlCorrection,
    PhaseMeanDefectContract,
    SignedAmplitudeProfile,
    CompactChebyshevSignedAmplitude,
    KokunoSignedAmplitudeCurlCorrection,
    dict[str, Any],
]:
    leading = KokunoLeadingCoreSeriesCandidate()
    primary = replace(
        KokunoCompleteCurlCorrection(),
        amplitude=ROUTED_OSCILLATORY_AMPLITUDE,
        phase=ROUTED_OSCILLATORY_PHASE,
    )
    contract = PhaseMeanDefectContract(phase_count=phase_count)
    theta_profile, axial_profile, projection = sample_real_phase_mean_radial_profiles(
        leading.at_points,
        primary,
        contract,
        time=PROFILE_TIME,
        z=PROFILE_Z,
        r_inner=PROFILE_ANNULUS[0],
        r_outer=PROFILE_ANNULUS[1],
        radial_count=radial_count,
        angular_count=angular_count,
    )
    inverse = CompactRadialStressInverse(theta_profile)
    target_stress = inverse.stress(theta_profile.radii)
    unit_covariance = measure_unit_radial_theta_covariance(
        primary,
        theta_profile.radii,
        time=PROFILE_TIME,
        z=PROFILE_Z,
        angular_count=angular_count,
        phase_count=phase_count,
    )
    sampled = build_rank1_signed_amplitude_profile(
        theta_profile.radii,
        target_stress,
        unit_covariance,
        base_amplitude=primary.amplitude,
    )
    compact = fit_compact_chebyshev_signed_amplitude(sampled)
    signed = KokunoSignedAmplitudeCurlCorrection(primary=primary, radial_amplitude=compact)
    adapter_metrics = compact.metrics(sampled.radii, sampled.delta_amplitude)
    predicted_covariance = sampled.covariance_derivative * compact(sampled.radii)
    active = sampled.active_target_mask
    covariance_error = predicted_covariance - sampled.target_stress
    target_scale = max(
        _scalar_rms(sampled.target_stress[active]) if np.any(active) else 0.0,
        np.finfo(float).tiny,
    )
    adapter_metrics.update({
        "fitted_linear_covariance_reconstruction_rms": (
            _scalar_rms(covariance_error[active]) if np.any(active) else 0.0
        ),
        "fitted_linear_covariance_reconstruction_relative_rms": (
            _scalar_rms(covariance_error[active]) / target_scale if np.any(active) else 0.0
        ),
        "sampled_signed_delta_amplitude_max_abs": float(
            np.max(np.abs(sampled.delta_amplitude))
        ),
        "sampled_signed_delta_amplitude_rms": _scalar_rms(sampled.delta_amplitude),
        "axial_ring_projection_rms": _scalar_rms(axial_profile.values),
        "real_defect_projection": projection,
        "radial_inverse": inverse.metrics(query_count=129),
    })
    return leading, primary, contract, sampled, compact, signed, adapter_metrics


def generate_signed_amplitude_curl_cycle_report(
    *,
    output: str | Path = "artifacts/kokuno_agent3/signed_amplitude_curl_cycle_report.json",
    radial_count: int = DEFAULT_RADIAL_COUNT,
    angular_count: int = DEFAULT_ANGULAR_COUNT,
    phase_count: int = DEFAULT_PHASE_COUNT,
    train_count: int = TRAIN_COUNT,
    held_out_count: int = HELD_OUT_COUNT,
) -> dict[str, Any]:
    """Materialize the signed amplitude and run one nonselecting correction cycle."""
    if train_count <= 0 or held_out_count <= 0:
        raise ValueError("sample counts must be positive")
    (
        leading,
        primary,
        contract,
        sampled,
        compact,
        signed,
        adapter_metrics,
    ) = _build_actual_signed_correction(
        radial_count=radial_count,
        angular_count=angular_count,
        phase_count=phase_count,
    )

    train_points = _sample_core_points(seed=TRAIN_SEED, count=train_count)
    held_out_points = _sample_core_points(seed=HELD_OUT_SEED, count=held_out_count)
    _assert_stencil_safe(leading, train_points, (TRAIN_TIME,), contract.step)
    _assert_stencil_safe(leading, held_out_points, HELD_OUT_TIMES, contract.step)

    train_sample = _evaluate_phase_mean_curl_cycle(
        leading, primary, signed, contract, train_points, TRAIN_TIME
    )
    train_before = _increment_metrics(train_sample.before_increment, train_points)
    train_after = _increment_metrics(train_sample.after_increment, train_points)
    training_ratio = train_after["rms"] / max(train_before["rms"], np.finfo(float).tiny)

    rows: list[dict[str, Any]] = []
    all_before: list[np.ndarray] = []
    all_after: list[np.ndarray] = []
    all_pressure_before: list[np.ndarray] = []
    all_pressure_after: list[np.ndarray] = []
    all_correction: list[np.ndarray] = []
    for time in HELD_OUT_TIMES:
        sample = _evaluate_phase_mean_curl_cycle(
            leading, primary, signed, contract, held_out_points, time
        )
        before = _increment_metrics(sample.before_increment, held_out_points)
        after = _increment_metrics(sample.after_increment, held_out_points)
        grad_p = _fixed_pressure_gradient(leading, held_out_points, time, contract)
        pressure_before = np.mean(sample.baseline_phase_momentum, axis=0) + grad_p
        pressure_after = sample.corrected_mean_momentum + grad_p
        correction_values = signed.at_points(held_out_points, time)
        rows.append({
            "time": float(time),
            "mean_defect_before": before,
            "mean_defect_after": after,
            "mean_defect_rms_ratio": after["rms"] / max(before["rms"], np.finfo(float).tiny),
            "theta_mean_defect_rms_ratio": after["theta_rms"] / max(before["theta_rms"], np.finfo(float).tiny),
            "pressure_inclusive_raw_phase_mean_operator_rms_before": _vector_rms(pressure_before),
            "pressure_inclusive_raw_phase_mean_operator_rms_after": _vector_rms(pressure_after),
            "pressure_inclusive_raw_phase_mean_operator_rms_ratio": (
                _vector_rms(pressure_after) / max(_vector_rms(pressure_before), np.finfo(float).tiny)
            ),
            "signed_correction_rms": _vector_rms(correction_values),
            "signed_correction_max": _vector_max(correction_values),
            "corrected_phase_divergence_max": float(np.max(np.abs(sample.corrected_phase_divergence))),
        })
        all_before.append(sample.before_increment)
        all_after.append(sample.after_increment)
        all_pressure_before.append(pressure_before)
        all_pressure_after.append(pressure_after)
        all_correction.append(correction_values)

    before_all = np.concatenate(all_before, axis=0)
    after_all = np.concatenate(all_after, axis=0)
    pressure_before_all = np.concatenate(all_pressure_before, axis=0)
    pressure_after_all = np.concatenate(all_pressure_after, axis=0)
    correction_all = np.concatenate(all_correction, axis=0)
    theta_before = cylindrical_theta_component(before_all, np.tile(held_out_points, (len(HELD_OUT_TIMES), 1)))
    theta_after = cylindrical_theta_component(after_all, np.tile(held_out_points, (len(HELD_OUT_TIMES), 1)))

    held_out_ratio = _vector_rms(after_all) / max(_vector_rms(before_all), np.finfo(float).tiny)
    theta_ratio = _scalar_rms(theta_after) / max(_scalar_rms(theta_before), np.finfo(float).tiny)
    pressure_ratio = _vector_rms(pressure_after_all) / max(
        _vector_rms(pressure_before_all), np.finfo(float).tiny
    )

    divergence_points = _sample_core_points(seed=DIVERGENCE_SEED, count=DIVERGENCE_COUNT)
    divergence = _independent_divergence(signed, divergence_points)
    divergence_rms = _scalar_rms(divergence)
    divergence_max = float(np.max(np.abs(divergence)))

    checks = {
        "compact_profile_fit_relative_rms_le_10pct": bool(
            compact.sampled_fit_relative_rms <= PROFILE_FIT_RELATIVE_RMS_LIMIT
        ),
        "held_in_total_mean_defect_improves": bool(training_ratio < 1.0),
        "held_out_total_mean_defect_improves": bool(held_out_ratio < 1.0),
        "held_out_theta_mean_defect_does_not_worsen": bool(theta_ratio <= 1.0),
        "held_out_pressure_inclusive_raw_operator_nonworsening_0p1pct": bool(
            pressure_ratio <= 1.001
        ),
        "independent_signed_correction_divergence_max_le_1e-5": bool(
            divergence_max <= 1.0e-5
        ),
        "signed_correction_nontrivial": bool(_vector_rms(correction_all) > 1.0e-12),
    }
    accepted = bool(all(checks.values()))

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
            "signed_differential_formula": (
                "d_Sigma=H^{-1}(Sigma/epsilon); delta_a_sigma=(d_Sigma)_sigma/(2*a_sigma)"
            ),
            "velocity_materialization_identity": (
                "A_delta=alpha(r)A_unit; curl(A_delta)=alpha curl(A_unit)+grad(alpha)xA_unit"
            ),
        },
        "inputs": {
            "leading_candidate_sha256": leading.sha256,
            "primary_oscillatory_amplitude": primary.amplitude,
            "primary_oscillatory_phase": primary.phase,
            "profile_annulus": list(PROFILE_ANNULUS),
            "profile_time": PROFILE_TIME,
            "profile_z": PROFILE_Z,
            "radial_count": radial_count,
            "angular_count": angular_count,
            "phase_count": phase_count,
            "chebyshev_degree_fixed_before_residual_screen": CHEBYSHEV_DEGREE,
            "profile_fit_relative_rms_limit": PROFILE_FIT_RELATIVE_RMS_LIMIT,
            "nu": contract.nu,
            "derivative_step": contract.step,
            "train_seed": TRAIN_SEED,
            "held_out_seed": HELD_OUT_SEED,
            "divergence_seed": DIVERGENCE_SEED,
            "train_count": train_count,
            "held_out_count": held_out_count,
            "held_out_times": list(HELD_OUT_TIMES),
            "surrogate_defect_used": False,
            "residual_points_used_to_fit_signed_profile": False,
        },
        "sampled_signed_covariance_profile": sampled.metrics(),
        "compact_amplitude_adapter": adapter_metrics,
        "public_signed_curl_correction": signed.metadata(),
        "held_in": {
            "time": TRAIN_TIME,
            "before": train_before,
            "after": train_after,
            "rms_ratio": training_ratio,
            "used_for_parameter_selection": False,
        },
        "held_out": {
            "rows": rows,
            "aggregate_mean_defect_rms_before": _vector_rms(before_all),
            "aggregate_mean_defect_rms_after": _vector_rms(after_all),
            "aggregate_mean_defect_rms_ratio": held_out_ratio,
            "aggregate_theta_rms_before": _scalar_rms(theta_before),
            "aggregate_theta_rms_after": _scalar_rms(theta_after),
            "aggregate_theta_rms_ratio": theta_ratio,
            "pressure_inclusive_raw_phase_mean_operator_rms_before": _vector_rms(pressure_before_all),
            "pressure_inclusive_raw_phase_mean_operator_rms_after": _vector_rms(pressure_after_all),
            "pressure_inclusive_raw_phase_mean_operator_rms_ratio": pressure_ratio,
            "signed_correction_rms": _vector_rms(correction_all),
            "signed_correction_max": _vector_max(correction_all),
        },
        "independent_signed_correction_divergence": {
            "seed": DIVERGENCE_SEED,
            "count": DIVERGENCE_COUNT,
            "step": 1.0e-5,
            "rms": divergence_rms,
            "max": divergence_max,
        },
        "finite_cycle_decision": {
            "checks": checks,
            "accepted_for_next_cycle": accepted,
            "no_post_result_amplitude_or_degree_widening": True,
        },
        "truth_boundary": {
            "real_candidate_defect_consumed": True,
            "radial_stress_reconstructed_from_real_defect": True,
            "public_velocity_correction_materialized": True,
            "agent2_complete_curl_reimplemented": False,
            "source_exact_two_column_signed_map_implemented": False,
            "source_exact_annular_hierarchy_implemented": False,
            "finite_correction_cycle_run": True,
            "residual_reduction_claimed": bool(held_out_ratio < 1.0),
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
        description="Materialize Agent-3 signed radial amplitude through Agent-2 complete curl and run one finite residual cycle"
    )
    parser.add_argument(
        "--output",
        default="artifacts/kokuno_agent3/signed_amplitude_curl_cycle_report.json",
    )
    parser.add_argument("--radial-count", type=int, default=DEFAULT_RADIAL_COUNT)
    parser.add_argument("--angular-count", type=int, default=DEFAULT_ANGULAR_COUNT)
    parser.add_argument("--phase-count", type=int, default=DEFAULT_PHASE_COUNT)
    parser.add_argument("--train-count", type=int, default=TRAIN_COUNT)
    parser.add_argument("--held-out-count", type=int, default=HELD_OUT_COUNT)
    args = parser.parse_args()
    report = generate_signed_amplitude_curl_cycle_report(
        output=args.output,
        radial_count=args.radial_count,
        angular_count=args.angular_count,
        phase_count=args.phase_count,
        train_count=args.train_count,
        held_out_count=args.held_out_count,
    )
    print(json.dumps({
        "task": report["task"],
        "compact_fit_relative_rms": report["compact_amplitude_adapter"]["sampled_fit_relative_rms"],
        "held_in_ratio": report["held_in"]["rms_ratio"],
        "held_out_ratio": report["held_out"]["aggregate_mean_defect_rms_ratio"],
        "held_out_theta_ratio": report["held_out"]["aggregate_theta_rms_ratio"],
        "pressure_operator_ratio": report["held_out"]["pressure_inclusive_raw_phase_mean_operator_rms_ratio"],
        "signed_correction_divergence_max": report["independent_signed_correction_divergence"]["max"],
        "accepted_for_next_cycle": report["finite_cycle_decision"]["accepted_for_next_cycle"],
    }, indent=2))


if __name__ == "__main__":
    main()
