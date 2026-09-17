"""Restricted-force cross-check for the frozen C2 compact Phi(1,0) screen.

Agent 7 PR #175 chooses the compact temporal window without using a PDE residual,
force, pressure or public image.  This module freezes that window first and asks one
narrow CR005 question: how does the already-preregistered two-parameter
``RestrictedForce(a,c)`` behave on the resulting public support-connected velocity
field relative to the static supported field and the derivative-balanced quartic
trial?

Only ``0 <= a,c <= 10`` are fitted, by the already-checked bounded linear variable
projection on fixed stratified training probes.  Coefficients are frozen on the
disjoint three-level holdout ladder.  No pressure basis, new force direction,
residual-defined force, temporal coefficient, spatial coefficient, support taper,
validation threshold or visual target is fitted here.

This is a pressure-free vorticity-equation routing diagnostic, not the preregistered
full-momentum acceptance gate.  It cannot establish PDE validity, visual
correspondence, paper exactness, OpenAI-field identity or blow-up.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_force_projection import (
    FIT_MAX_ITER,
    FIT_TOL,
    FORCE_BOUNDS,
    project_supported_restricted_force,
)
from .constrained_eq45_supported_phi10_c2_compact_temporal_capacity import (
    compact_phi10_snapshot,
    compact_window_parameters,
)
from .constrained_eq45_supported_phi10_force_crosscheck import _project_temporal_trial
from .constrained_eq45_supported_phi10_quadratic_temporal_capacity import (
    DEFAULT_EARLY_DELTA,
)
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)

TASK_ID = "CR005-EQ45-SUPPORTED-PHI10-C2-COMPACT-FORCE-CROSSCHECK-030"
FROZEN_EARLY_DELTA = float(DEFAULT_EARLY_DELTA)


@dataclass(frozen=True)
class _CompactTemporalField:
    """Thin public-velocity wrapper around #175's already-screened snapshot map."""

    base: Eq45SupportedVelocityCandidate
    early_delta: float = FROZEN_EARLY_DELTA

    def __post_init__(self) -> None:
        if not isinstance(self.base, Eq45SupportedVelocityCandidate):
            raise TypeError("base must be Eq45SupportedVelocityCandidate")
        value = float(self.early_delta)
        if not np.isfinite(value) or value == 0.0:
            raise ValueError("early_delta must be finite and nonzero")
        # Fail closed on the frozen schedule before the PDE operator is touched.
        compact_window_parameters(early_delta=value)
        compact_phi10_snapshot(self.base, self.time_start, early_delta=value)

    @property
    def time_start(self) -> float:
        return float(self.base.time_start)

    @property
    def time_end(self) -> float:
        return float(self.base.time_end)

    @property
    def time_midpoint(self) -> float:
        return 0.5 * (self.time_start + self.time_end)

    def snapshot(self, time: float) -> Eq45SupportedVelocityCandidate:
        time_array = np.asarray(time, dtype=float)
        if time_array.ndim != 0:
            raise ValueError("snapshot time must be scalar")
        return compact_phi10_snapshot(
            self.base, float(time_array), early_delta=float(self.early_delta)
        )

    def at_points(self, points, time) -> np.ndarray:
        points_array = np.asarray(points, dtype=float)
        if points_array.ndim == 0 or points_array.shape[-1] != 3:
            raise ValueError("points must have shape (...,3)")
        if not np.all(np.isfinite(points_array)):
            raise ValueError("points must be finite")
        time_array = np.asarray(time, dtype=float)
        if not np.all(np.isfinite(time_array)):
            raise ValueError("time must be finite")
        try:
            time_array = np.broadcast_to(time_array, points_array.shape[:-1])
        except ValueError as exc:
            raise ValueError("time must broadcast to points.shape[:-1]") from exc
        if np.any((time_array < self.time_start) | (time_array > self.time_end)):
            raise ValueError("time lies outside the candidate delivery interval")

        flat_points = points_array.reshape((-1, 3))
        flat_times = np.asarray(time_array, dtype=float).reshape(-1)
        flat_output = np.empty_like(flat_points)
        for scalar_time in np.unique(flat_times):
            mask = flat_times == scalar_time
            flat_output[mask] = self.snapshot(float(scalar_time)).at_points(
                flat_points[mask], float(scalar_time)
            )
        output = flat_output.reshape(points_array.shape)
        if not np.all(np.isfinite(output)):
            raise FloatingPointError("compact temporal field returned nonfinite velocity")
        return output

    velocity = at_points
    __call__ = at_points

    @property
    def sha256(self) -> str:
        parameters = compact_window_parameters(early_delta=float(self.early_delta))
        payload = {
            "schema": "eq45_supported_phi10_c2_compact_screen_field_v1",
            "base_sha256": self.base.sha256,
            "early_delta": float(self.early_delta),
            "window_tau_length": float(parameters["window_tau_length"]),
            "return_tau": float(parameters["return_tau"]),
        }
        encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()


def _fractional_change(new: float, old: float) -> float:
    return float((new - old) / max(abs(old), np.finfo(float).tiny))


def compare_c2_compact_phi10_force(
    *, constraints_path: str | Path | None = None
) -> dict:
    """Compare static, quartic and frozen compact fields under one force contract."""
    root = Path(__file__).resolve().parents[2]
    constraints_path = Path(constraints_path or root / "configs/constraints.json")
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise AssertionError("active forcing convention drifted from preregistration")

    static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    quartic = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
        base=static, early_delta=FROZEN_EARLY_DELTA
    )
    compact = _CompactTemporalField(base=static, early_delta=FROZEN_EARLY_DELTA)
    parameters = compact_window_parameters(early_delta=FROZEN_EARLY_DELTA)

    # Freeze velocity identities before any PDE/force calculation.
    static_sha = static.sha256
    quartic_sha = quartic.sha256
    compact_sha = compact.sha256

    nu = float(constraints["nu"])
    steps = tuple(float(value) for value in constraints["validation"]["derivative_steps"])
    static_report = project_supported_restricted_force(
        static, nu=nu, derivative_steps=steps
    )
    quartic_report = _project_temporal_trial(
        quartic, nu=nu, derivative_steps=steps
    )
    compact_report = _project_temporal_trial(
        compact, nu=nu, derivative_steps=steps
    )

    if (static.sha256, quartic.sha256, compact.sha256) != (
        static_sha,
        quartic_sha,
        compact_sha,
    ):
        raise AssertionError("velocity identity changed during restricted-force projection")

    probes = np.asarray(
        [
            [0.37, 0.11, -0.22],
            [0.95, -0.35, 0.48],
            [1.45, 0.15, -0.62],
            [1.72, 0.00, 0.30],
        ],
        dtype=float,
    )
    return_time = 0.5 * (static.time_start + static.time_end) + 0.5 * (
        static.time_end - static.time_start
    ) * float(parameters["return_tau"])
    public_identity_checks = {
        "early_compact_equals_quartic": bool(
            np.array_equal(
                compact.at_points(probes, static.time_start),
                quartic.at_points(probes, static.time_start),
            )
        ),
        "post_return_t04375_compact_equals_static": bool(
            np.array_equal(
                compact.at_points(probes, 0.4375),
                static.at_points(probes, 0.4375),
            )
        ),
        "post_return_t05625_compact_equals_static": bool(
            np.array_equal(
                compact.at_points(probes, 0.5625),
                static.at_points(probes, 0.5625),
            )
        ),
        "post_return_t06875_compact_equals_static": bool(
            np.array_equal(
                compact.at_points(probes, 0.6875),
                static.at_points(probes, 0.6875),
            )
        ),
        "active_t0375_compact_differs_from_static": bool(
            not np.array_equal(
                compact.at_points(probes, 0.375),
                static.at_points(probes, 0.375),
            )
        ),
        "active_t0375_compact_differs_from_quartic": bool(
            not np.array_equal(
                compact.at_points(probes, 0.375),
                quartic.at_points(probes, 0.375),
            )
        ),
    }
    if not all(public_identity_checks.values()):
        raise AssertionError("compact public-velocity identity contract drifted")

    static_finest = static_report["holdout_rows"][-1]
    quartic_finest = quartic_report["holdout_rows"][-1]
    compact_finest = compact_report["holdout_rows"][-1]

    static_zero = float(static_finest["before"]["rms"])
    quartic_zero = float(quartic_finest["before"]["rms"])
    compact_zero = float(compact_finest["before"]["rms"])
    static_forced = float(static_finest["after"]["rms"])
    quartic_forced = float(quartic_finest["after"]["rms"])
    compact_forced = float(compact_finest["after"]["rms"])

    return {
        "schema": "eq45_supported_phi10_c2_compact_force_crosscheck_v1",
        "task_id": TASK_ID,
        "claim_scope": "bounded_restricted_force_tradeoff_for_frozen_c2_compact_temporal_screen",
        "velocity_schedule_frozen_before_pde_fit": True,
        "early_delta": FROZEN_EARLY_DELTA,
        "window_tau_length": float(parameters["window_tau_length"]),
        "return_tau": float(parameters["return_tau"]),
        "return_time": float(return_time),
        "window_selection_rule": "earliest_C2_return_under_quartic_peak_coefficient_slope_cap",
        "pde_objective_used_to_choose_temporal_schedule": False,
        "mode": {"family": "phi", "index": [1, 0]},
        "static_supported_sha256": static_sha,
        "quartic_trial_sha256": quartic_sha,
        "compact_screen_sha256": compact_sha,
        "public_identity_checks": public_identity_checks,
        "nu": nu,
        "derivative_steps": list(steps),
        "training_random_seed": None,
        "training_randomness_used": False,
        "training_probe_contract": "fixed_stratified_4x4_plateau_radial_axial_corner",
        "holdout_probe_contract": "agent3_supported_vorticity_fixed_16_region_probes",
        "fit_and_holdout_separate": True,
        "holdout_force_refit": False,
        "force_family": "preregistered_restricted_two_parameter_family",
        "force_bounds": list(FORCE_BOUNDS),
        "force_solver": "scipy.optimize.lsq_linear_bounded_linear_variable_projection",
        "force_solver_max_iterations": FIT_MAX_ITER,
        "force_solver_tolerance": FIT_TOL,
        "pressure_fitted": False,
        "new_force_direction_added": False,
        "residual_defined_force_allowed": False,
        "temporal_schedule_refit_on_pde": False,
        "static_projection": {
            "fit": static_report["fit"],
            "holdout_rows": static_report["holdout_rows"],
        },
        "quartic_projection": quartic_report,
        "compact_projection": compact_report,
        "finest_holdout_comparison": {
            "spatial_step": float(compact_finest["spatial_step"]),
            "static_zero_force_rms": static_zero,
            "static_projected_force_rms": static_forced,
            "quartic_zero_force_rms": quartic_zero,
            "quartic_projected_force_rms": quartic_forced,
            "compact_zero_force_rms": compact_zero,
            "compact_projected_force_rms": compact_forced,
            "compact_vs_static_zero_force_fractional_change": _fractional_change(
                compact_zero, static_zero
            ),
            "compact_vs_static_projected_force_fractional_change": _fractional_change(
                compact_forced, static_forced
            ),
            "compact_vs_quartic_zero_force_fractional_change": _fractional_change(
                compact_zero, quartic_zero
            ),
            "compact_vs_quartic_projected_force_fractional_change": _fractional_change(
                compact_forced, quartic_forced
            ),
            "compact_force_reduction_fraction": float(
                compact_finest["rms_reduction_fraction"]
            ),
        },
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "fixed-probe pressure-free vorticity RMS is not the preregistered "
            "volume-weighted full-momentum L2/global maximum and pressure is not bound"
        ),
        "visualization_candidate_only": True,
        "canonical_velocity_changed": False,
        "production_temporal_shape_promoted": False,
        "visual_correspondence_verified": False,
        "visualization_ready": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = compare_c2_compact_phi10_force()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
