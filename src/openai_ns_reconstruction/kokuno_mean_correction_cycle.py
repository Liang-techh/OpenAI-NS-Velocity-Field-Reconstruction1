"""One bounded mean-velocity correction screen driven by real phase-mean defect.

KokunoYumeto's corrected 2026-09-09 reader sends compact radial stress through
a signed covariance-amplitude map; it does not identify stress with mean
velocity. This repository does not yet expose that full amplitude hierarchy.
We therefore use the reconstructed theta/e=2 stress only as a fixed *shape* for
one compact, axisymmetric mean-swirl velocity correction. The lift is an
autonomous engineering adapter, not a paper-exact identity.

A scalar amplitude is selected on held-in samples only. The frozen nonzero
correction is then evaluated on disjoint held-out samples/times through the real
public candidate and the independent CR006 finite-difference operator. No force
is fitted and no acceptance threshold is changed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import argparse
import json
from math import pi
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.interpolate import CubicSpline

from .constrained_validation import residual
from .kokuno_complete_curl import KokunoCompleteCurlCorrection
from .kokuno_mean_defect import PhaseMeanDefectContract, cylindrical_theta_component
from .kokuno_radial_stress import CompactRadialStressInverse, sample_real_phase_mean_radial_profiles

SOURCE_READER = "KokunoYumeto corrected 208-page reconstruction"
SOURCE_EDITION = "2026.09.09-consolidated"
SOURCE_DOI = "10.5281/zenodo.22678406"
SOURCE_RADIAL_FORMULA = "C12 / R6-R9 weighted moment-complement radial stress"
SOURCE_CYCLE_CONTEXT = "finite correction cycle: source stress -> signed covariance-amplitude map"
ENGINEERING_LIFT = "normalized theta/e=2 radial stress -> compact axisymmetric mean-swirl velocity"
AMPLITUDE_GRID = (-0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04)
TRAIN_SEED = 9172813
HELD_OUT_SEED = 9172819
TRAIN_COUNT = 64
HELD_OUT_COUNT = 128
TRAIN_TIME = 0.5
HELD_OUT_TIMES = (0.375, 0.5, 0.625)
AXIAL_HALF_WIDTH = 0.58
RADIAL_SAMPLE_BOUNDS = (0.16, 0.68)
Z_SAMPLE_HALF_WIDTH = 0.44
PointVelocity = Callable[[np.ndarray, float], np.ndarray]


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))


def _vector_max(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.max(np.linalg.norm(values, axis=-1)))


def _zero_pressure(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros(len(points), dtype=float)


def _zero_force(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros_like(points, dtype=float)


def _c4_centered_window(values: np.ndarray, half_width: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if not np.isfinite(half_width) or half_width <= 0.0:
        raise ValueError("half_width must be finite and positive")
    s = values / half_width
    inside = np.abs(s) < 1.0
    base = np.where(inside, 1.0 - s * s, 0.0)
    return base**5


@dataclass(frozen=True)
class StressShapedMeanSwirlCorrection:
    """Steady compact axisymmetric pure-swirl engineering correction."""

    inverse: CompactRadialStressInverse
    amplitude: float
    axial_half_width: float = AXIAL_HALF_WIDTH
    normalization_count: int = 513
    _stress_scale: float = field(init=False, repr=False, compare=False)
    _radial_spline: CubicSpline = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.inverse.profile.component != "theta" or self.inverse.profile.exponent != 2:
            raise ValueError("mean-swirl lift requires the theta/e=2 radial stress")
        if not np.isfinite(self.amplitude) or abs(self.amplitude) > 0.04:
            raise ValueError("amplitude must lie in the preregistered [-0.04, 0.04] range")
        if not np.isfinite(self.axial_half_width) or not 0.1 <= self.axial_half_width <= 1.0:
            raise ValueError("axial_half_width must lie in [0.1, 1.0]")
        if self.normalization_count < 65:
            raise ValueError("normalization_count must be >= 65")
        radii = np.linspace(self.inverse.r_inner, self.inverse.r_outer, self.normalization_count)
        stress = self.inverse.stress(radii)
        scale = float(np.max(np.abs(stress)))
        if not np.isfinite(scale) or scale <= 1.0e-14:
            raise ValueError("radial stress is too small to define a normalized correction shape")
        object.__setattr__(self, "_stress_scale", scale)
        object.__setattr__(self, "_radial_spline", CubicSpline(radii, stress / scale, extrapolate=False))

    @property
    def stress_scale(self) -> float:
        return self._stress_scale

    def at_points(self, points: np.ndarray, time: float) -> np.ndarray:
        del time
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError("points must be a finite (N,3) array")
        radius = np.hypot(points[:, 0], points[:, 1])
        radial_shape = np.zeros(len(points), dtype=float)
        inside = (radius > self.inverse.r_inner) & (radius < self.inverse.r_outer)
        if np.any(inside):
            radial_shape[inside] = self._radial_spline(radius[inside])
        axial_shape = _c4_centered_window(points[:, 2], self.axial_half_width)
        scalar = self.amplitude * radial_shape * axial_shape
        out = np.zeros_like(points)
        nonaxis = radius > 1.0e-14
        out[nonaxis, 0] = -scalar[nonaxis] * points[nonaxis, 1] / radius[nonaxis]
        out[nonaxis, 1] = scalar[nonaxis] * points[nonaxis, 0] / radius[nonaxis]
        return out

    __call__ = at_points

    def metadata(self) -> dict:
        return {
            "family": "kokuno_agent3_stress_shaped_mean_swirl_v1",
            "amplitude": float(self.amplitude),
            "amplitude_bound": 0.04,
            "axial_half_width": float(self.axial_half_width),
            "radial_support": [self.inverse.r_inner, self.inverse.r_outer],
            "stress_scale": self.stress_scale,
            "source_stress_formula": SOURCE_RADIAL_FORMULA,
            "engineering_lift": ENGINEERING_LIFT,
            "source_exact_velocity_lift": False,
        }


def _sample_cylindrical_points(*, seed: int, count: int,
                               radial_bounds: tuple[float, float] = RADIAL_SAMPLE_BOUNDS,
                               z_half_width: float = Z_SAMPLE_HALF_WIDTH) -> np.ndarray:
    if count <= 0:
        raise ValueError("count must be positive")
    r0, r1 = map(float, radial_bounds)
    if not 0.0 < r0 < r1 or not 0.0 < z_half_width:
        raise ValueError("invalid cylindrical sampling domain")
    rng = np.random.default_rng(seed)
    radius = np.sqrt(rng.uniform(r0 * r0, r1 * r1, count))
    angle = rng.uniform(0.0, 2.0 * pi, count)
    z = rng.uniform(-z_half_width, z_half_width, count)
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def _compose(base_velocity: PointVelocity,
             correction: StressShapedMeanSwirlCorrection) -> PointVelocity:
    def velocity(points: np.ndarray, time: float) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        base = np.asarray(base_velocity(points, time), dtype=float)
        delta = correction.at_points(points, time)
        if base.shape != points.shape:
            raise ValueError("base velocity must return (N,3) values")
        return base + delta
    return velocity


def _phase_mean_after(base_velocity: PointVelocity, oscillatory: KokunoCompleteCurlCorrection,
                      mean_correction: StressShapedMeanSwirlCorrection, points: np.ndarray,
                      time: float, contract: PhaseMeanDefectContract,
                      original_base_momentum: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    corrected_base = _compose(base_velocity, mean_correction)
    sample = contract.evaluate(corrected_base, oscillatory, points, float(time))
    increment_relative_original_base = sample.mean_momentum - original_base_momentum
    return increment_relative_original_base, sample.mean_momentum, float(np.max(np.abs(sample.phase_divergence)))


def _increment_metrics(values: np.ndarray, points: np.ndarray) -> dict[str, float]:
    theta = cylindrical_theta_component(values, points)
    return {
        "rms": _vector_rms(values),
        "max": _vector_max(values),
        "theta_rms": float(np.sqrt(np.mean(theta * theta))),
        "theta_max": float(np.max(np.abs(theta))),
    }


def _independent_divergence(correction: StressShapedMeanSwirlCorrection,
                            points: np.ndarray, *, step: float = 1.0e-5) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    divergence = np.zeros(len(points), dtype=float)
    for axis in range(3):
        delta = np.eye(3)[axis] * step
        plus = correction.at_points(points + delta, TRAIN_TIME)
        minus = correction.at_points(points - delta, TRAIN_TIME)
        divergence += ((plus - minus) / (2.0 * step))[:, axis]
    return divergence


def _train_amplitude(base_velocity: PointVelocity, oscillatory: KokunoCompleteCurlCorrection,
                     inverse: CompactRadialStressInverse, contract: PhaseMeanDefectContract,
                     points: np.ndarray) -> tuple[float, list[dict], dict]:
    baseline = contract.evaluate(base_velocity, oscillatory, points, TRAIN_TIME)
    before_metrics = _increment_metrics(baseline.mean_defect_increment, points)
    rows: list[dict] = []
    for amplitude in AMPLITUDE_GRID:
        mean_correction = StressShapedMeanSwirlCorrection(inverse=inverse, amplitude=float(amplitude))
        after, _, phase_divergence_max = _phase_mean_after(
            base_velocity, oscillatory, mean_correction, points, TRAIN_TIME, contract,
            baseline.base_momentum)
        metrics = _increment_metrics(after, points)
        rows.append({
            "amplitude": float(amplitude),
            **metrics,
            "ratio_to_zero_mean_defect_rms": metrics["rms"] / max(before_metrics["rms"], np.finfo(float).tiny),
            "composite_phase_divergence_max": phase_divergence_max,
        })
    selected = min(rows, key=lambda row: row["rms"])
    return float(selected["amplitude"]), rows, before_metrics


def generate_current_candidate_mean_cycle_report(
    *,
    candidate_path: str | Path = "artifacts/bipolar_joint_capped/candidate.json",
    output: str | Path = "artifacts/kokuno_agent3/mean_correction_cycle_report.json",
    radial_count: int = 65,
    angular_count: int = 16,
) -> dict:
    """Run one bounded stress-shaped mean-correction screen on a real artifact."""
    from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate

    candidate_path = Path(candidate_path)
    candidate = Eq45SupportedVelocityCandidate.load_json(candidate_path)
    oscillatory = KokunoCompleteCurlCorrection()
    contract = PhaseMeanDefectContract(phase_count=8)
    theta_profile, _, projection = sample_real_phase_mean_radial_profiles(
        candidate.at_points, oscillatory, contract, time=TRAIN_TIME, z=0.0,
        r_inner=0.12, r_outer=0.72, radial_count=radial_count, angular_count=angular_count)
    inverse = CompactRadialStressInverse(theta_profile)

    train_points = _sample_cylindrical_points(seed=TRAIN_SEED, count=TRAIN_COUNT)
    held_out_points = _sample_cylindrical_points(seed=HELD_OUT_SEED, count=HELD_OUT_COUNT)
    selected_amplitude, training_rows, training_before = _train_amplitude(
        candidate.at_points, oscillatory, inverse, contract, train_points)
    selected = StressShapedMeanSwirlCorrection(inverse=inverse, amplitude=selected_amplitude)

    held_out_rows = []
    aggregate_before_increment, aggregate_after_increment = [], []
    aggregate_before_operator, aggregate_after_operator = [], []
    for time in HELD_OUT_TIMES:
        before_sample = contract.evaluate(candidate.at_points, oscillatory, held_out_points, time)
        before_increment = before_sample.mean_defect_increment
        after_increment, after_mean_operator, phase_divergence_max = _phase_mean_after(
            candidate.at_points, oscillatory, selected, held_out_points, time, contract,
            before_sample.base_momentum)
        before_metrics = _increment_metrics(before_increment, held_out_points)
        after_metrics = _increment_metrics(after_increment, held_out_points)
        before_operator_rms = _vector_rms(before_sample.mean_momentum)
        after_operator_rms = _vector_rms(after_mean_operator)
        correction_values = selected.at_points(held_out_points, time)
        held_out_rows.append({
            "time": float(time),
            "mean_defect_before": before_metrics,
            "mean_defect_after": after_metrics,
            "mean_defect_rms_ratio": after_metrics["rms"] / max(before_metrics["rms"], np.finfo(float).tiny),
            "raw_phase_mean_operator_rms_before": before_operator_rms,
            "raw_phase_mean_operator_rms_after": after_operator_rms,
            "raw_phase_mean_operator_rms_ratio": after_operator_rms / max(before_operator_rms, np.finfo(float).tiny),
            "selected_correction_rms": _vector_rms(correction_values),
            "selected_correction_max": _vector_max(correction_values),
            "composite_phase_divergence_max": phase_divergence_max,
        })
        aggregate_before_increment.append(before_increment)
        aggregate_after_increment.append(after_increment)
        aggregate_before_operator.append(before_sample.mean_momentum)
        aggregate_after_operator.append(after_mean_operator)

    before_all = np.concatenate(aggregate_before_increment, axis=0)
    after_all = np.concatenate(aggregate_after_increment, axis=0)
    before_operator_all = np.concatenate(aggregate_before_operator, axis=0)
    after_operator_all = np.concatenate(aggregate_after_operator, axis=0)
    heldout_before_rms, heldout_after_rms = _vector_rms(before_all), _vector_rms(after_all)
    heldout_ratio = heldout_after_rms / max(heldout_before_rms, np.finfo(float).tiny)
    raw_before_rms, raw_after_rms = _vector_rms(before_operator_all), _vector_rms(after_operator_all)
    raw_ratio = raw_after_rms / max(raw_before_rms, np.finfo(float).tiny)

    div_points = _sample_cylindrical_points(seed=9172829, count=192,
                                            radial_bounds=(0.15, 0.69), z_half_width=0.46)
    correction_divergence = _independent_divergence(selected, div_points)
    correction_divergence_max = float(np.max(np.abs(correction_divergence)))
    correction_divergence_rms = float(np.sqrt(np.mean(correction_divergence * correction_divergence)))

    train_selected_row = next(row for row in training_rows if row["amplitude"] == selected_amplitude)
    train_improved = bool(train_selected_row["rms"] < training_before["rms"])
    heldout_improved = bool(heldout_after_rms < heldout_before_rms)
    raw_operator_nonworsening = bool(raw_ratio <= 1.001)
    divergence_guard = bool(correction_divergence_max <= 1.0e-5)
    nontrivial = bool(abs(selected_amplitude) > 0.0 and selected.stress_scale > 1.0e-14)
    accepted_for_next_cycle = bool(
        train_improved and heldout_improved and raw_operator_nonworsening and divergence_guard and nontrivial)

    report = {
        "task": "KOKUNO-A3-MEAN-CORRECTION-CYCLE-003",
        "source": {
            "reader": SOURCE_READER,
            "edition": SOURCE_EDITION,
            "doi": SOURCE_DOI,
            "radial_formula": SOURCE_RADIAL_FORMULA,
            "cycle_context": SOURCE_CYCLE_CONTEXT,
            "important_source_boundary": "source radial stress is realized through a signed covariance-amplitude map; it is not identified there with mean velocity",
            "source_cycle_gain_diagnostic_only": {
                "wave_residual_gain_exponent": 0.39999,
                "tangential_mean_gain_exponent": 0.17,
                "defect_gain_exponent": 0.89996,
                "used_as_numeric_contraction_gate": False,
            },
        },
        "candidate_path": str(candidate_path),
        "candidate_sha256": candidate.sha256,
        "oscillatory_correction": oscillatory.metadata(),
        "real_defect_contract": {
            "nu": contract.nu,
            "derivative_step": contract.step,
            "phase_count": contract.phase_count,
            "time_bounds": list(contract.time_bounds),
            "pressure_for_increment": "unchanged terms cancel by subtracting original base operator; no pressure fitted",
            "force_for_increment": "unchanged terms cancel by subtracting original base operator; no force fitted",
        },
        "stress_shape_input": {
            "projection": projection,
            "theta_e2_metrics": inverse.metrics(),
            "radial_count": radial_count,
            "angular_count": angular_count,
            "profile_time": TRAIN_TIME,
            "profile_z": 0.0,
        },
        "engineering_mean_velocity_lift": selected.metadata(),
        "selection": {
            "objective": "held-in phase-mean nonlinear defect increment RMS relative to original base operator",
            "amplitude_grid": list(AMPLITUDE_GRID),
            "zero_amplitude_is_baseline_not_a_selectable_success": True,
            "train_seed": TRAIN_SEED,
            "train_count": TRAIN_COUNT,
            "train_time": TRAIN_TIME,
            "training_before": training_before,
            "training_rows": training_rows,
            "selected_amplitude": selected_amplitude,
            "selected_amplitude_at_grid_boundary": bool(abs(selected_amplitude) == max(abs(x) for x in AMPLITUDE_GRID)),
            "training_improved": train_improved,
        },
        "held_out": {
            "seed": HELD_OUT_SEED,
            "count": HELD_OUT_COUNT,
            "times": list(HELD_OUT_TIMES),
            "points_reused_from_training": False,
            "rows": held_out_rows,
            "aggregate_mean_defect_rms_before": heldout_before_rms,
            "aggregate_mean_defect_rms_after": heldout_after_rms,
            "aggregate_mean_defect_rms_ratio": heldout_ratio,
            "aggregate_raw_phase_mean_operator_rms_before": raw_before_rms,
            "aggregate_raw_phase_mean_operator_rms_after": raw_after_rms,
            "aggregate_raw_phase_mean_operator_rms_ratio": raw_ratio,
            "heldout_mean_defect_improved": heldout_improved,
            "raw_operator_nonworsening_0p1pct_guard": raw_operator_nonworsening,
        },
        "correction_checks": {
            "independent_divergence_seed": 9172829,
            "independent_divergence_points": 192,
            "independent_divergence_fd_step": 1.0e-5,
            "correction_divergence_rms": correction_divergence_rms,
            "correction_divergence_max": correction_divergence_max,
            "divergence_guard_1e-5": divergence_guard,
            "nontrivial_correction": nontrivial,
            "radial_support": [inverse.r_inner, inverse.r_outer],
            "axial_support": [-selected.axial_half_width, selected.axial_half_width],
        },
        "cycle_decision": {
            "accepted_for_next_cycle": accepted_for_next_cycle,
            "criteria": [
                "nonzero fixed-grid amplitude improves held-in mean-defect increment RMS",
                "same frozen amplitude improves aggregate disjoint held-out mean-defect increment RMS",
                "held-out raw phase-mean zero-pressure/zero-force operator does not worsen by >0.1%",
                "independent correction divergence max <= 1e-5",
            ],
            "if_rejected": "retain failed bounded correction attempt; do not widen the amplitude grid in this run",
        },
        "truth_boundary": {
            "consumes_real_candidate_artifact": True,
            "consumes_real_phase_mean_defect": True,
            "surrogate_defect_used": False,
            "stress_to_velocity_lift_is_autonomous_engineering_adapter": True,
            "source_exact_velocity_lift": False,
            "mean_velocity_correction_evaluated": True,
            "finite_correction_cycle_one_step_screen_run": True,
            "full_pressure_force_model_fitted": False,
            "formal_normalized_pde_gate_assessed": False,
            "normalized_ns_residual_le_1e-3_claimed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Run one bounded stress-shaped Kokuno Agent-3 mean correction cycle")
    parser.add_argument("--candidate", default="artifacts/bipolar_joint_capped/candidate.json")
    parser.add_argument("--output", default="artifacts/kokuno_agent3/mean_correction_cycle_report.json")
    parser.add_argument("--radial-count", type=int, default=65)
    parser.add_argument("--angular-count", type=int, default=16)
    args = parser.parse_args()
    report = generate_current_candidate_mean_cycle_report(
        candidate_path=args.candidate, output=args.output,
        radial_count=args.radial_count, angular_count=args.angular_count)
    print(json.dumps({
        "candidate_sha256": report["candidate_sha256"],
        "selected_amplitude": report["selection"]["selected_amplitude"],
        "training_improved": report["selection"]["training_improved"],
        "held_out_aggregate_ratio": report["held_out"]["aggregate_mean_defect_rms_ratio"],
        "raw_operator_ratio": report["held_out"]["aggregate_raw_phase_mean_operator_rms_ratio"],
        "correction_checks": report["correction_checks"],
        "cycle_decision": report["cycle_decision"],
    }, indent=2))


if __name__ == "__main__":
    main()
