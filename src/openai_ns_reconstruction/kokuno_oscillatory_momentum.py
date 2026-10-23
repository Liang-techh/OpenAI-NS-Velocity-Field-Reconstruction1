"""Truth-bounded preliminary momentum screen for the Agent-2 complete curl.

This module answers one narrow question: with every non-amplitude choice in the
existing :class:`KokunoCompleteCurlCorrection` frozen, how does adding that
correction change the *raw* zero-pressure/zero-force momentum operator on a real
repository velocity artifact?

The correction remains the autonomous Cartesian surrogate documented in
``kokuno_complete_curl``.  KokunoYumeto's corrected 2026-09-09 reader motivates
complete curls and potential-level localization, but it does not supply the
Cartesian affine phase, C4 box support, parameter values, or the amplitude grid
used here.  Nothing in this module is a paper-exact field identification.

For a base velocity ``b`` and correction ``w`` we evaluate the exact algebraic
increment

    N(b+w)-N(b)
      = w_t + (b . grad)w + (w . grad)b + (w . grad)w - nu Delta w,

with Agent-2's analytic ``w_t``, ``grad w`` and ``Delta w``.  Only ``grad b`` is
finite-differenced.  The resulting hybrid momentum is independently compared
with ``constrained_validation.residual`` applied to the full composite velocity.
That comparison is a derivative-consistency check, not the registered PDE gate.

The current report deliberately uses the capped-bipolar candidate only as a
``temporary_engineering_bridge_not_kokuno_leading`` until Agent 1 exposes a full
3-D Kokuno leading field.  Pressure and forcing are fixed to zero in both sides
of the comparison.  No amplitude is selected or promoted from this held-out
screen, and any phase-mean/zero-harmonic burden remains Agent 3's lane.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import argparse
import json
from pathlib import Path
from typing import Callable, Sequence

import numpy as np

from .constrained_validation import residual
from .kokuno_complete_curl import KokunoCompleteCurlCorrection


PointVelocity = Callable[[np.ndarray, float], np.ndarray]
DEFAULT_AMPLITUDES = (0.0, 0.0625, 0.125, 0.25)
DEFAULT_TIMES = (0.375, 0.5, 0.625)
DEFAULT_STEPS = (0.02, 0.01, 0.005)


def _zero_pressure(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros(len(points), dtype=float)


def _zero_force(points: np.ndarray, time: float) -> np.ndarray:
    del time
    return np.zeros_like(points, dtype=float)


def _vector_rms(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt(np.mean(np.sum(values * values, axis=-1))))


def _vector_max(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.max(np.linalg.norm(values, axis=-1)))


def _validate_points(points: np.ndarray) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0:
        raise ValueError("points must be a nonempty (N,3) array")
    if not np.isfinite(points).all():
        raise ValueError("points must be finite")
    return points


def _base_jacobian(
    base_velocity: PointVelocity,
    points: np.ndarray,
    time: float,
    step: float,
) -> np.ndarray:
    """Fourth-order Cartesian Jacobian of only the base velocity."""
    points = _validate_points(points)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("step must be positive and finite")
    jacobian = np.empty((len(points), 3, 3), dtype=float)
    for axis in range(3):
        delta = np.eye(3, dtype=float)[axis] * step
        um2, um1, up1, up2 = [
            np.asarray(base_velocity(points + offset * delta, float(time)), dtype=float)
            for offset in (-2.0, -1.0, 1.0, 2.0)
        ]
        if any(value.shape != points.shape for value in (um2, um1, up1, up2)):
            raise ValueError("base_velocity must return an (N,3) array")
        jacobian[:, :, axis] = (um2 - 8.0 * um1 + 8.0 * up1 - up2) / (12.0 * step)
    if not np.isfinite(jacobian).all():
        raise RuntimeError("base Jacobian became nonfinite")
    return jacobian


@dataclass(frozen=True)
class OscillatoryMomentumSample:
    """One direct-vs-analytic comparison at a fixed amplitude/time/step."""

    amplitude: float
    time: float
    step: float
    base_momentum: np.ndarray
    analytic_increment: np.ndarray
    hybrid_composite_momentum: np.ndarray
    direct_composite_momentum: np.ndarray
    direct_composite_divergence: np.ndarray
    analytic_correction_divergence: np.ndarray

    @property
    def agreement_error(self) -> np.ndarray:
        return self.direct_composite_momentum - self.hybrid_composite_momentum

    def metrics(self) -> dict[str, float]:
        base_rms = _vector_rms(self.base_momentum)
        composite_rms = _vector_rms(self.direct_composite_momentum)
        return {
            "amplitude": float(self.amplitude),
            "time": float(self.time),
            "step": float(self.step),
            "base_momentum_rms": base_rms,
            "base_momentum_max": _vector_max(self.base_momentum),
            "composite_momentum_rms": composite_rms,
            "composite_momentum_max": _vector_max(self.direct_composite_momentum),
            "composite_minus_base_rms": composite_rms - base_rms,
            "composite_over_base_rms": composite_rms / max(base_rms, np.finfo(float).tiny),
            "analytic_increment_rms": _vector_rms(self.analytic_increment),
            "analytic_increment_max": _vector_max(self.analytic_increment),
            "direct_vs_hybrid_rms": _vector_rms(self.agreement_error),
            "direct_vs_hybrid_max": _vector_max(self.agreement_error),
            "direct_composite_divergence_max": float(np.max(np.abs(self.direct_composite_divergence))),
            "analytic_correction_divergence_max": float(np.max(np.abs(self.analytic_correction_divergence))),
        }


@dataclass(frozen=True)
class OscillatoryMomentumIncrementContract:
    """Analytic correction increment with an independent full-composite replay."""

    nu: float = 0.01
    time_bounds: tuple[float, float] = (0.25, 0.75)

    def __post_init__(self) -> None:
        if not np.isfinite([self.nu, *self.time_bounds]).all() or self.nu <= 0.0:
            raise ValueError("nu and time bounds must be finite, with nu > 0")
        if self.time_bounds[0] >= self.time_bounds[1]:
            raise ValueError("time_bounds must be increasing")

    def evaluate(
        self,
        base_velocity: PointVelocity,
        correction: KokunoCompleteCurlCorrection,
        points: np.ndarray,
        time: float,
        *,
        step: float,
    ) -> OscillatoryMomentumSample:
        if not callable(base_velocity):
            raise TypeError("base_velocity must be callable")
        if not isinstance(correction, KokunoCompleteCurlCorrection):
            raise TypeError("correction must be KokunoCompleteCurlCorrection")
        points = _validate_points(points)
        if not np.isfinite(time) or not self.time_bounds[0] <= time <= self.time_bounds[1]:
            raise ValueError("time must lie inside time_bounds")
        if not np.isfinite(step) or step <= 0.0 or 4.0 * step > self.time_bounds[1] - self.time_bounds[0]:
            raise ValueError("step is invalid for the declared time interval")

        base_values = np.asarray(base_velocity(points, float(time)), dtype=float)
        if base_values.shape != points.shape or not np.isfinite(base_values).all():
            raise ValueError("base_velocity must return finite (N,3) values")

        base_result = residual(
            base_velocity,
            _zero_pressure,
            _zero_force,
            points,
            float(time),
            nu=self.nu,
            step=float(step),
            time_bounds=self.time_bounds,
        )
        base_jacobian = _base_jacobian(base_velocity, points, float(time), float(step))

        x, y, z = points.T
        wave = correction.at_points(points, float(time))
        wave_t = correction.time_derivative(x, y, z, float(time))
        wave_jacobian = correction.spatial_jacobian(x, y, z, float(time))
        wave_laplacian = correction.laplacian(x, y, z, float(time))
        wave_divergence = correction.divergence(x, y, z, float(time))

        analytic_increment = (
            wave_t
            + np.einsum("nij,nj->ni", wave_jacobian, base_values)
            + np.einsum("nij,nj->ni", base_jacobian, wave)
            + np.einsum("nij,nj->ni", wave_jacobian, wave)
            - self.nu * wave_laplacian
        )
        base_momentum = np.asarray(base_result["momentum"], dtype=float)
        hybrid = base_momentum + analytic_increment

        def composite_velocity(sample_points: np.ndarray, sample_time: float) -> np.ndarray:
            base = np.asarray(base_velocity(sample_points, sample_time), dtype=float)
            delta = correction.at_points(sample_points, sample_time)
            if base.shape != sample_points.shape or delta.shape != sample_points.shape:
                raise ValueError("base and correction must return (N,3) values")
            return base + delta

        direct = residual(
            composite_velocity,
            _zero_pressure,
            _zero_force,
            points,
            float(time),
            nu=self.nu,
            step=float(step),
            time_bounds=self.time_bounds,
        )
        direct_momentum = np.asarray(direct["momentum"], dtype=float)
        if not all(np.isfinite(value).all() for value in (analytic_increment, hybrid, direct_momentum, wave_divergence)):
            raise RuntimeError("momentum screen produced nonfinite values")

        return OscillatoryMomentumSample(
            amplitude=float(correction.amplitude),
            time=float(time),
            step=float(step),
            base_momentum=base_momentum,
            analytic_increment=analytic_increment,
            hybrid_composite_momentum=hybrid,
            direct_composite_momentum=direct_momentum,
            direct_composite_divergence=np.asarray(direct["divergence"], dtype=float),
            analytic_correction_divergence=np.asarray(wave_divergence, dtype=float),
        )


def _held_out_points(
    correction: KokunoCompleteCurlCorrection,
    *,
    seed: int,
    count: int,
) -> np.ndarray:
    """Fresh interior support points, independent of prior Agent-2/3/4 seeds."""
    if count <= 0:
        raise ValueError("count must be positive")
    center = np.asarray(correction.center, dtype=float)
    widths = np.asarray(correction.half_widths, dtype=float)
    rng = np.random.default_rng(seed)
    return rng.uniform(center - 0.72 * widths, center + 0.72 * widths, size=(count, 3))


def screen_frozen_amplitudes(
    base_velocity: PointVelocity,
    points: np.ndarray,
    *,
    amplitudes: Sequence[float] = DEFAULT_AMPLITUDES,
    times: Sequence[float] = DEFAULT_TIMES,
    steps: Sequence[float] = DEFAULT_STEPS,
    correction_template: KokunoCompleteCurlCorrection | None = None,
    contract: OscillatoryMomentumIncrementContract | None = None,
) -> list[dict[str, float]]:
    """Evaluate a predeclared amplitude grid without selecting a winner."""
    points = _validate_points(points)
    template = correction_template or KokunoCompleteCurlCorrection()
    contract = contract or OscillatoryMomentumIncrementContract()
    amplitudes = tuple(float(value) for value in amplitudes)
    times = tuple(float(value) for value in times)
    steps = tuple(float(value) for value in steps)
    if not amplitudes or not times or not steps:
        raise ValueError("amplitudes, times and steps must be nonempty")
    if any(not np.isfinite(value) or abs(value) > 2.0 for value in amplitudes):
        raise ValueError("screen amplitudes must lie in the existing [-2,2] bound")

    rows: list[dict[str, float]] = []
    for amplitude in amplitudes:
        correction = replace(template, amplitude=amplitude)
        for time in times:
            for step in steps:
                rows.append(
                    contract.evaluate(
                        base_velocity,
                        correction,
                        points,
                        time,
                        step=step,
                    ).metrics()
                )
    return rows


def generate_current_candidate_report(
    *,
    candidate_path: str | Path = "artifacts/bipolar_joint_capped/candidate.json",
    output: str | Path = "artifacts/kokuno_agent2/oscillatory_momentum_screen.json",
    seed: int = 9172791,
    point_count: int = 96,
) -> dict:
    """Run the frozen screen on the temporary real engineering bridge."""
    from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate

    candidate_path = Path(candidate_path)
    candidate = Eq45SupportedVelocityCandidate.load_json(candidate_path)
    template = KokunoCompleteCurlCorrection()
    points = _held_out_points(template, seed=seed, count=point_count)
    rows = screen_frozen_amplitudes(candidate.at_points, points, correction_template=template)

    finest = min(DEFAULT_STEPS)
    finest_rows = [row for row in rows if row["step"] == finest]
    by_amplitude: dict[str, dict[str, float]] = {}
    for amplitude in DEFAULT_AMPLITUDES:
        selected = [row for row in finest_rows if row["amplitude"] == amplitude]
        by_amplitude[str(amplitude)] = {
            "mean_base_momentum_rms": float(np.mean([row["base_momentum_rms"] for row in selected])),
            "mean_composite_momentum_rms": float(np.mean([row["composite_momentum_rms"] for row in selected])),
            "mean_composite_over_base_rms": float(np.mean([row["composite_over_base_rms"] for row in selected])),
            "max_direct_vs_hybrid_rms": float(np.max([row["direct_vs_hybrid_rms"] for row in selected])),
            "max_analytic_correction_divergence": float(np.max([row["analytic_correction_divergence_max"] for row in selected])),
        }

    report = {
        "task": "K2-OSC-003",
        "source": {
            "reader": "KokunoYumeto corrected 208-page reconstruction",
            "edition": "2026.09.09-consolidated",
            "doi": "10.5281/zenodo.22678406",
            "source_structure_used": [
                "oscillatory velocity is represented by a complete curl",
                "potential-level localization retains the cutoff-gradient curl remainder",
            ],
            "autonomous_screen_choices": [
                "Agent-2 Cartesian affine phase/frame",
                "Agent-2 separable compact C4 box support",
                "frozen amplitude grid [0, 0.0625, 0.125, 0.25]",
                "zero pressure and zero forcing for this preliminary difference screen",
            ],
        },
        "bridge": {
            "candidate_path": str(candidate_path),
            "candidate_sha256": candidate.sha256,
            "classification": "temporary_engineering_bridge_not_kokuno_leading",
        },
        "sampling": {
            "seed": int(seed),
            "point_count": int(point_count),
            "times": list(DEFAULT_TIMES),
            "derivative_steps": list(DEFAULT_STEPS),
            "amplitudes": list(DEFAULT_AMPLITUDES),
            "nu": 0.01,
            "points_are_optimization_samples": False,
            "amplitude_grid_frozen_before_evaluation": True,
        },
        "operator": {
            "analytic_correction_terms": ["u_t", "spatial_jacobian", "laplacian"],
            "base_jacobian": "fourth-order centered finite difference",
            "direct_cross_check": "constrained_validation.residual on the full composite velocity",
            "pressure": "zero/fixed on both base and composite",
            "forcing": "zero/fixed on both base and composite",
        },
        "rows": rows,
        "finest_step_summary": by_amplitude,
        "truth_boundary": {
            "amplitude_selected": False,
            "candidate_promoted": False,
            "temporary_bridge_is_kokuno_leading": False,
            "mean_correction_applied": False,
            "pressure_or_forcing_fitted": False,
            "formal_pde_gate_assessed": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Preliminary frozen-amplitude momentum screen for the Agent-2 complete curl"
    )
    parser.add_argument("--candidate", default="artifacts/bipolar_joint_capped/candidate.json")
    parser.add_argument("--output", default="artifacts/kokuno_agent2/oscillatory_momentum_screen.json")
    parser.add_argument("--seed", type=int, default=9172791)
    parser.add_argument("--points", type=int, default=96)
    args = parser.parse_args()
    report = generate_current_candidate_report(
        candidate_path=args.candidate,
        output=args.output,
        seed=args.seed,
        point_count=args.points,
    )
    print(json.dumps({
        "candidate_sha256": report["bridge"]["candidate_sha256"],
        "sampling": report["sampling"],
        "finest_step_summary": report["finest_step_summary"],
        "truth_boundary": report["truth_boundary"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
