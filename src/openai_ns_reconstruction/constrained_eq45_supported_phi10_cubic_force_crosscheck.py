"""Restricted-force cross-check for the frozen cubic-localized supported Phi(1,0) trial.

The cubic temporal shape is fixed before this module evaluates any PDE diagnostic.
Only the preregistered ``RestrictedForce(a,c)`` family with ``0 <= a,c <= 10`` is
projected, using the already-checked supported-child vorticity operator, fixed
stratified training probes, and disjoint held-out derivative ladder.  No pressure
basis, force direction, temporal coefficient, spatial basis coefficient, taper,
or acceptance threshold is fitted here.

This is a routing diagnostic for the final callable ``[u,v,w]``.  It asks whether
the sharper cubic temporal localization materialized in PR #159 has a different
pressure-free vorticity obstruction than the previously checked quadratic trial.
It is not a full momentum validation and cannot promote PDE validity, visual
correspondence, paper exactness, or identification of an OpenAI field.
"""
from __future__ import annotations

import argparse
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
from .constrained_eq45_supported_phi10_cubic_temporal_capacity import (
    DEFAULT_EARLY_DELTA,
    cubic_phi10_delta_dtau,
)
from .constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_early_localized_temporal_mode import (
    Eq45SupportedPhi10EarlyLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_force_crosscheck import _project_temporal_trial

FROZEN_EARLY_DELTA = float(DEFAULT_EARLY_DELTA)


def _fractional_change(new: float, old: float) -> float:
    return float((new - old) / max(abs(old), np.finfo(float).tiny))


def _quadratic_delta_dtau_at_start(early_delta: float) -> float:
    # delta_q(tau) = 0.5 * early_delta * tau * (tau - 1)
    return float(0.5 * early_delta * (2.0 * -1.0 - 1.0))


def compare_cubic_localized_phi10_force(
    *, constraints_path: str | Path | None = None
) -> dict:
    """Compare static, quadratic and cubic fields under one frozen force contract."""
    root = Path(__file__).resolve().parents[2]
    constraints_path = Path(constraints_path or root / "configs/constraints.json")
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise AssertionError("active forcing convention drifted from preregistration")

    static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    quadratic = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(
        base=static, early_delta=FROZEN_EARLY_DELTA
    )
    cubic = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(
        base=static, early_delta=FROZEN_EARLY_DELTA
    )

    # Freeze all velocity-side choices before looking at the PDE diagnostic.
    static_sha = static.sha256
    quadratic_sha = quadratic.sha256
    cubic_sha = cubic.sha256

    nu = float(constraints["nu"])
    steps = tuple(float(value) for value in constraints["validation"]["derivative_steps"])
    static_report = project_supported_restricted_force(
        static, nu=nu, derivative_steps=steps
    )
    quadratic_report = _project_temporal_trial(
        quadratic, nu=nu, derivative_steps=steps
    )
    cubic_report = _project_temporal_trial(
        cubic, nu=nu, derivative_steps=steps
    )

    # Fail closed if a helper mutated or replaced a frozen candidate.
    if (static.sha256, quadratic.sha256, cubic.sha256) != (
        static_sha,
        quadratic_sha,
        cubic_sha,
    ):
        raise AssertionError("velocity candidate identity changed during force projection")

    probes = np.asarray(
        [
            [0.37, 0.11, -0.22],
            [0.95, -0.35, 0.48],
            [1.45, 0.15, -0.62],
            [1.72, 0.00, 0.30],
        ],
        dtype=float,
    )
    public_identity_checks = {
        "early_cubic_equals_quadratic": bool(
            np.array_equal(
                cubic.at_points(probes, cubic.time_start),
                quadratic.at_points(probes, quadratic.time_start),
            )
        ),
        "midpoint_cubic_equals_static": bool(
            np.array_equal(
                cubic.at_points(probes, cubic.time_midpoint),
                static.at_points(probes, cubic.time_midpoint),
            )
        ),
        "t0625_cubic_equals_static": bool(
            np.array_equal(
                cubic.at_points(probes, 0.625),
                static.at_points(probes, 0.625),
            )
        ),
        "late_cubic_equals_static": bool(
            np.array_equal(
                cubic.at_points(probes, cubic.time_end),
                static.at_points(probes, cubic.time_end),
            )
        ),
    }
    if not all(public_identity_checks.values()):
        raise AssertionError("cubic public-velocity identity contract drifted")

    static_finest = static_report["holdout_rows"][-1]
    quadratic_finest = quadratic_report["holdout_rows"][-1]
    cubic_finest = cubic_report["holdout_rows"][-1]

    static_zero = float(static_finest["before"]["rms"])
    quadratic_zero = float(quadratic_finest["before"]["rms"])
    cubic_zero = float(cubic_finest["before"]["rms"])
    static_forced = float(static_finest["after"]["rms"])
    quadratic_forced = float(quadratic_finest["after"]["rms"])
    cubic_forced = float(cubic_finest["after"]["rms"])

    cubic_start_derivative = abs(
        float(cubic_phi10_delta_dtau(-1.0, early_delta=FROZEN_EARLY_DELTA))
    )
    quadratic_start_derivative = abs(
        _quadratic_delta_dtau_at_start(FROZEN_EARLY_DELTA)
    )

    return {
        "schema": "eq45_supported_phi10_cubic_force_crosscheck_v1",
        "task_id": "CR005-EQ45-SUPPORTED-PHI10-CUBIC-FORCE-CROSSCHECK-024",
        "claim_scope": "bounded_restricted_force_tradeoff_for_frozen_cubic_visualization_trial",
        "velocity_schedule_frozen_before_pde_fit": True,
        "early_delta": FROZEN_EARLY_DELTA,
        "mode": {"family": "phi", "index": [1, 0]},
        "static_supported_sha256": static_sha,
        "quadratic_trial_sha256": quadratic_sha,
        "cubic_trial_sha256": cubic_sha,
        "public_identity_checks": public_identity_checks,
        "nu": nu,
        "derivative_steps": list(steps),
        "training_random_seed": None,
        "training_randomness_used": False,
        "training_probe_contract": "fixed_stratified_4x4_plateau_radial_axial_corner",
        "holdout_probe_contract": "agent3_supported_vorticity_fixed_16_region_probes",
        "fit_and_holdout_separate": True,
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
        "quadratic_projection": quadratic_report,
        "cubic_projection": cubic_report,
        "temporal_sharpness_context": {
            "quadratic_abs_d_delta_dtau_at_start": quadratic_start_derivative,
            "cubic_abs_d_delta_dtau_at_start": cubic_start_derivative,
            "cubic_to_quadratic_start_derivative_ratio": float(
                cubic_start_derivative / quadratic_start_derivative
            ),
        },
        "finest_holdout_comparison": {
            "spatial_step": float(cubic_finest["spatial_step"]),
            "static_zero_force_rms": static_zero,
            "static_projected_force_rms": static_forced,
            "quadratic_zero_force_rms": quadratic_zero,
            "quadratic_projected_force_rms": quadratic_forced,
            "cubic_zero_force_rms": cubic_zero,
            "cubic_projected_force_rms": cubic_forced,
            "quadratic_vs_static_zero_force_fractional_change": _fractional_change(
                quadratic_zero, static_zero
            ),
            "quadratic_vs_static_projected_force_fractional_change": _fractional_change(
                quadratic_forced, static_forced
            ),
            "cubic_vs_static_zero_force_fractional_change": _fractional_change(
                cubic_zero, static_zero
            ),
            "cubic_vs_static_projected_force_fractional_change": _fractional_change(
                cubic_forced, static_forced
            ),
            "cubic_vs_quadratic_zero_force_fractional_change": _fractional_change(
                cubic_zero, quadratic_zero
            ),
            "cubic_vs_quadratic_projected_force_fractional_change": _fractional_change(
                cubic_forced, quadratic_forced
            ),
            "cubic_force_reduction_fraction": float(
                cubic_finest["rms_reduction_fraction"]
            ),
        },
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "fixed-probe pressure-free vorticity RMS is not the preregistered "
            "volume-weighted full-momentum L2/global maximum and pressure is not bound"
        ),
        "visualization_candidate_only": True,
        "canonical_velocity_changed": False,
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
    report = compare_cubic_localized_phi10_force()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
