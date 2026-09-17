"""Restricted-force cross-check for the early-localized supported Phi(1,0) trial.

This module combines two already-checked dependencies: the serializable quadratic
Phi(1,0) temporal candidate and Agent-2's independent supported-child vorticity /
restricted-force projection.  The quadratic schedule is frozen before any PDE
calculation.  Only the preregistered RestrictedForce(a,c), 0 <= a,c <= 10, is
projected on the existing fixed training probes and then frozen on the disjoint
holdout probes.

The result is a routing diagnostic for the final callable [u,v,w].  It does not
fit pressure, add a force direction, re-optimize the temporal schedule, relax a
threshold, fit an OpenAI image, or promote solver convergence to PDE validity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_phi10_early_localized_temporal_mode import (
    Eq45SupportedPhi10EarlyLocalizedTemporalCandidate,
)
from .constrained_eq45_supported_phi10_force_crosscheck import (
    FROZEN_SLOPE,
    FORCE_BOUNDS,
    _project_temporal_trial,
    compare_supported_phi10_temporal_force,
)
from .constrained_eq45_supported_temporal_mode import (
    Eq45SupportedAffineTemporalModeCandidate,
)

FROZEN_EARLY_DELTA = -1.4


def _fractional_change(new: float, old: float) -> float:
    return float((new - old) / max(abs(old), np.finfo(float).tiny))


def compare_early_localized_phi10_force(
    *, constraints_path: str | Path | None = None
) -> dict:
    """Compare static, affine and early-localized fields under one force contract."""
    root = Path(__file__).resolve().parents[2]
    constraints_path = Path(constraints_path or root / "configs/constraints.json")
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise AssertionError("active forcing convention drifted from preregistration")

    # Reuse the complete checked affine/static projection from PR #144 rather than
    # reimplementing its operator, probe sets or force projection.
    affine_report = compare_supported_phi10_temporal_force(
        slope=FROZEN_SLOPE, constraints_path=constraints_path
    )

    static = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    if static.sha256 != affine_report["static_supported_sha256"]:
        raise AssertionError("static supported identity drifted across dependencies")
    quadratic = Eq45SupportedPhi10EarlyLocalizedTemporalCandidate(
        base=static, early_delta=FROZEN_EARLY_DELTA
    )

    nu = float(constraints["nu"])
    steps = tuple(float(value) for value in constraints["validation"]["derivative_steps"])
    quadratic_report = _project_temporal_trial(
        quadratic,
        nu=nu,
        derivative_steps=steps,
    )

    # Fail closed on the intended endpoint relationship between the two visual
    # trials and the static supported field.  These are public evaluator checks,
    # not coefficient-only checks.
    affine = Eq45SupportedAffineTemporalModeCandidate(
        base=static, family="phi", mode_i=1, mode_j=0, slope=FROZEN_SLOPE
    )
    probes = np.asarray(
        [
            [0.37, 0.11, -0.22],
            [0.95, -0.35, 0.48],
            [1.45, 0.15, -0.62],
            [1.72, 0.00, 0.30],
        ],
        dtype=float,
    )
    endpoint_checks = {
        "early_equals_affine": bool(
            np.array_equal(
                quadratic.at_points(probes, quadratic.time_start),
                affine.at_points(probes, affine.time_start),
            )
        ),
        "midpoint_equals_static": bool(
            np.array_equal(
                quadratic.at_points(probes, quadratic.time_midpoint),
                static.at_points(probes, quadratic.time_midpoint),
            )
        ),
        "late_equals_static": bool(
            np.array_equal(
                quadratic.at_points(probes, quadratic.time_end),
                static.at_points(probes, quadratic.time_end),
            )
        ),
    }
    if not all(endpoint_checks.values()):
        raise AssertionError("early-localized public velocity endpoint identity drifted")

    static_finest = affine_report["static_projection"]["holdout_rows"][-1]
    affine_finest = affine_report["temporal_projection"]["holdout_rows"][-1]
    quadratic_finest = quadratic_report["holdout_rows"][-1]

    static_zero = float(static_finest["before"]["rms"])
    static_forced = float(static_finest["after"]["rms"])
    affine_zero = float(affine_finest["before"]["rms"])
    affine_forced = float(affine_finest["after"]["rms"])
    quadratic_zero = float(quadratic_finest["before"]["rms"])
    quadratic_forced = float(quadratic_finest["after"]["rms"])

    return {
        "schema": "eq45_supported_phi10_early_localized_force_crosscheck_v1",
        "claim_scope": "bounded_restricted_force_tradeoff_for_frozen_early_localized_visualization_trial",
        "schedule_frozen_before_pde_fit": True,
        "early_delta": FROZEN_EARLY_DELTA,
        "affine_reference_slope": FROZEN_SLOPE,
        "mode": {"family": "phi", "index": [1, 0]},
        "static_supported_sha256": static.sha256,
        "affine_trial_sha256": affine_report["temporal_trial_sha256"],
        "early_localized_trial_sha256": quadratic.sha256,
        "endpoint_checks": endpoint_checks,
        "nu": nu,
        "derivative_steps": list(steps),
        "force_family": "preregistered_restricted_two_parameter_family",
        "force_bounds": list(FORCE_BOUNDS),
        "pressure_fitted": False,
        "new_force_direction_added": False,
        "residual_defined_force_allowed": False,
        "temporal_schedule_refit_on_pde": False,
        "static_projection": affine_report["static_projection"],
        "affine_projection": affine_report["temporal_projection"],
        "early_localized_projection": quadratic_report,
        "finest_holdout_comparison": {
            "spatial_step": float(quadratic_finest["spatial_step"]),
            "static_zero_force_rms": static_zero,
            "static_projected_force_rms": static_forced,
            "affine_zero_force_rms": affine_zero,
            "affine_projected_force_rms": affine_forced,
            "early_localized_zero_force_rms": quadratic_zero,
            "early_localized_projected_force_rms": quadratic_forced,
            "affine_vs_static_zero_force_fractional_change": _fractional_change(
                affine_zero, static_zero
            ),
            "affine_vs_static_projected_force_fractional_change": _fractional_change(
                affine_forced, static_forced
            ),
            "early_localized_vs_static_zero_force_fractional_change": _fractional_change(
                quadratic_zero, static_zero
            ),
            "early_localized_vs_static_projected_force_fractional_change": _fractional_change(
                quadratic_forced, static_forced
            ),
            "early_localized_vs_affine_zero_force_fractional_change": _fractional_change(
                quadratic_zero, affine_zero
            ),
            "early_localized_vs_affine_projected_force_fractional_change": _fractional_change(
                quadratic_forced, affine_forced
            ),
            "early_localized_force_reduction_fraction": float(
                quadratic_finest["rms_reduction_fraction"]
            ),
        },
        "formal_pde_gate_assessed": False,
        "reason_formal_pde_gate_unassessed": (
            "fixed-probe pressure-free vorticity RMS is not the preregistered "
            "volume-weighted full-momentum L2/global maximum and pressure is not bound"
        ),
        "visualization_candidate_only": True,
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
    report = compare_early_localized_phi10_force()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
