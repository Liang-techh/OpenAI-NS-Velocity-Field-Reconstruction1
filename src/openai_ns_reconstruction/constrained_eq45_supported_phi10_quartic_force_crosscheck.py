"""Restricted-force cross-check for the frozen derivative-balanced quartic Phi(1,0) trial.

The quartic temporal schedule is fixed by the target-free derivative-balancing rule
from the representation lane before this module evaluates any PDE diagnostic.  This
module only reuses the preregistered ``RestrictedForce(a,c)`` family with
``0 <= a,c <= 10`` and the already checked independent supported-child vorticity
operator.  Force coefficients are fitted on fixed stratified training probes and
then frozen on a disjoint three-level held-out derivative ladder.

No pressure basis, new force direction, temporal coefficient, spatial coefficient,
taper parameter, acceptance threshold, or visual target is fitted here.  The report
is a routing diagnostic for the final callable ``[u,v,w]`` candidate; it is not a
full momentum validation and cannot establish PDE validity, visual correspondence,
paper exactness, identification of an OpenAI field, or blow-up.
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
from .constrained_eq45_supported_phi10_cubic_temporal_mode import (
    Eq45SupportedPhi10CubicLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_force_crosscheck import _project_temporal_trial
from .constrained_eq45_supported_phi10_quartic_derivative_capacity import (
    DEFAULT_EARLY_DELTA,
    balanced_nullspace_coefficient,
)
from .constrained_eq45_supported_phi10_quartic_temporal_mode import (
    Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate,
)

FROZEN_EARLY_DELTA = float(DEFAULT_EARLY_DELTA)
FROZEN_NULLSPACE_COEFFICIENT = float(
    balanced_nullspace_coefficient(early_delta=FROZEN_EARLY_DELTA)
)


def _fractional_change(new: float, old: float) -> float:
    return float((new - old) / max(abs(old), np.finfo(float).tiny))


def compare_quartic_balanced_phi10_force(
    *, constraints_path: str | Path | None = None
) -> dict:
    """Compare static, cubic and quartic fields under one frozen force contract."""
    root = Path(__file__).resolve().parents[2]
    constraints_path = Path(constraints_path or root / "configs/constraints.json")
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise AssertionError("active forcing convention drifted from preregistration")

    static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    cubic = Eq45SupportedPhi10CubicLocalizedTemporalCandidate(
        base=static, early_delta=FROZEN_EARLY_DELTA
    )
    quartic = Eq45SupportedPhi10QuarticDerivativeBalancedTemporalCandidate(
        base=static, early_delta=FROZEN_EARLY_DELTA
    )

    if quartic.nullspace_coefficient != FROZEN_NULLSPACE_COEFFICIENT:
        raise AssertionError("quartic derivative-balancing rule drifted")

    # Freeze every velocity-side choice before the PDE diagnostic is evaluated.
    static_sha = static.sha256
    cubic_sha = cubic.sha256
    quartic_sha = quartic.sha256

    nu = float(constraints["nu"])
    steps = tuple(float(value) for value in constraints["validation"]["derivative_steps"])
    static_report = project_supported_restricted_force(
        static, nu=nu, derivative_steps=steps
    )
    cubic_report = _project_temporal_trial(
        cubic, nu=nu, derivative_steps=steps
    )
    quartic_report = _project_temporal_trial(
        quartic, nu=nu, derivative_steps=steps
    )

    if (static.sha256, cubic.sha256, quartic.sha256) != (
        static_sha,
        cubic_sha,
        quartic_sha,
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
        "early_quartic_equals_cubic": bool(
            np.array_equal(
                quartic.at_points(probes, quartic.time_start),
                cubic.at_points(probes, cubic.time_start),
            )
        ),
        "midpoint_quartic_equals_static": bool(
            np.array_equal(
                quartic.at_points(probes, quartic.time_midpoint),
                static.at_points(probes, quartic.time_midpoint),
            )
        ),
        "t0625_quartic_equals_static": bool(
            np.array_equal(
                quartic.at_points(probes, 0.625),
                static.at_points(probes, 0.625),
            )
        ),
        "late_quartic_equals_static": bool(
            np.array_equal(
                quartic.at_points(probes, quartic.time_end),
                static.at_points(probes, quartic.time_end),
            )
        ),
        "off_keyframe_quartic_differs_from_cubic": bool(
            any(
                not np.array_equal(
                    quartic.at_points(probes, time),
                    cubic.at_points(probes, time),
                )
                for time in (0.375, 0.5625, 0.6875)
            )
        ),
    }
    if not all(public_identity_checks.values()):
        raise AssertionError("quartic public-velocity identity contract drifted")

    static_finest = static_report["holdout_rows"][-1]
    cubic_finest = cubic_report["holdout_rows"][-1]
    quartic_finest = quartic_report["holdout_rows"][-1]

    static_zero = float(static_finest["before"]["rms"])
    cubic_zero = float(cubic_finest["before"]["rms"])
    quartic_zero = float(quartic_finest["before"]["rms"])
    static_forced = float(static_finest["after"]["rms"])
    cubic_forced = float(cubic_finest["after"]["rms"])
    quartic_forced = float(quartic_finest["after"]["rms"])

    return {
        "schema": "eq45_supported_phi10_quartic_force_crosscheck_v1",
        "task_id": "CR005-EQ45-SUPPORTED-PHI10-QUARTIC-FORCE-CROSSCHECK-027",
        "claim_scope": "bounded_restricted_force_tradeoff_for_frozen_derivative_balanced_quartic_trial",
        "velocity_schedule_frozen_before_pde_fit": True,
        "early_delta": FROZEN_EARLY_DELTA,
        "nullspace_coefficient": FROZEN_NULLSPACE_COEFFICIENT,
        "nullspace_rule": "least_squares_minimize_return_node_coefficient_slope",
        "pde_objective_used_to_choose_temporal_schedule": False,
        "mode": {"family": "phi", "index": [1, 0]},
        "static_supported_sha256": static_sha,
        "cubic_trial_sha256": cubic_sha,
        "quartic_trial_sha256": quartic_sha,
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
        "cubic_projection": cubic_report,
        "quartic_projection": quartic_report,
        "finest_holdout_comparison": {
            "spatial_step": float(quartic_finest["spatial_step"]),
            "static_zero_force_rms": static_zero,
            "static_projected_force_rms": static_forced,
            "cubic_zero_force_rms": cubic_zero,
            "cubic_projected_force_rms": cubic_forced,
            "quartic_zero_force_rms": quartic_zero,
            "quartic_projected_force_rms": quartic_forced,
            "cubic_vs_static_zero_force_fractional_change": _fractional_change(
                cubic_zero, static_zero
            ),
            "cubic_vs_static_projected_force_fractional_change": _fractional_change(
                cubic_forced, static_forced
            ),
            "quartic_vs_static_zero_force_fractional_change": _fractional_change(
                quartic_zero, static_zero
            ),
            "quartic_vs_static_projected_force_fractional_change": _fractional_change(
                quartic_forced, static_forced
            ),
            "quartic_vs_cubic_zero_force_fractional_change": _fractional_change(
                quartic_zero, cubic_zero
            ),
            "quartic_vs_cubic_projected_force_fractional_change": _fractional_change(
                quartic_forced, cubic_forced
            ),
            "quartic_force_reduction_fraction": float(
                quartic_finest["rms_reduction_fraction"]
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
    report = compare_quartic_balanced_phi10_force()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
