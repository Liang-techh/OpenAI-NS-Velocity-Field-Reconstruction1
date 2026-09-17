"""Bounded restricted-force screen across the frozen compact--quartic Phi10 blend.

Agent 7 PR #184 preregistered the target-free blend grid and Agent 1 PR #186
materialized the callable blend candidate.  This CR005 increment does not choose a
blend weight.  It freezes the five already-screened weights before any PDE
calculation, fits only the preregistered ``RestrictedForce(a,c)`` family on the
existing fixed stratified training probes, and freezes those coefficients on the
independent three-level supported-child holdout.

The report is a pressure-free vorticity-equation routing diagnostic.  It is intended
to tell the final ``[u,v,w]`` selector whether visual interpolation between the
quartic and compact endpoints has a materially different PDE-side cost.  It does not
fit pressure, add a force direction, fit a public image, select a blend weight, or
promote solver convergence to PDE validity.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from .constrained_eq45_candidate import Eq45VelocityCandidate
from .constrained_eq45_supported_candidate import Eq45SupportedVelocityCandidate
from .constrained_eq45_supported_force_projection import (
    FIT_MAX_ITER,
    FIT_TOL,
    FORCE_BOUNDS,
)
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)
from .constrained_eq45_supported_phi10_force_crosscheck import _project_temporal_trial

TASK_ID = "CR005-EQ45-SUPPORTED-PHI10-COMPACT-QUARTIC-BLEND-FORCE-SCREEN-033"
FROZEN_BLEND_WEIGHTS = (0.0, 0.25, 0.5, 0.75, 1.0)
FROZEN_EARLY_DELTA = -1.4


def _fractional_change(new: float, old: float) -> float:
    return float((new - old) / max(abs(old), np.finfo(float).tiny))


def _candidate_roundtrip(candidate):
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "candidate.json"
        candidate.save_json(path)
        return type(candidate).load_json(path)


def screen_compact_quartic_blend_force(
    *, constraints_path: str | Path | None = None
) -> dict:
    """Evaluate the frozen blend grid under one bounded force-projection contract."""
    root = Path(__file__).resolve().parents[2]
    constraints_path = Path(constraints_path or root / "configs/constraints.json")
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise AssertionError("active forcing convention drifted from preregistration")

    base = Eq45SupportedVelocityCandidate(parent=Eq45VelocityCandidate.seed())
    nu = float(constraints["nu"])
    steps = tuple(float(value) for value in constraints["validation"]["derivative_steps"])

    # Freeze every public velocity identity before touching a PDE residual.
    candidates = tuple(
        Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
            base=base,
            blend_weight=weight,
            early_delta=FROZEN_EARLY_DELTA,
        )
        for weight in FROZEN_BLEND_WEIGHTS
    )
    frozen_shas = tuple(candidate.sha256 for candidate in candidates)

    rows = []
    for weight, candidate, frozen_sha in zip(FROZEN_BLEND_WEIGHTS, candidates, frozen_shas):
        projection = _project_temporal_trial(
            candidate,
            nu=nu,
            derivative_steps=steps,
        )
        if candidate.sha256 != frozen_sha:
            raise AssertionError("blend candidate identity changed during force projection")
        finest = projection["holdout_rows"][-1]
        rows.append(
            {
                "blend_weight": float(weight),
                "candidate_sha256": frozen_sha,
                "fit": projection["fit"],
                "holdout_rows": projection["holdout_rows"],
                "finest_zero_force_rms": float(finest["before"]["rms"]),
                "finest_projected_force_rms": float(finest["after"]["rms"]),
                "finest_zero_force_max": float(finest["before"]["max"]),
                "finest_projected_force_max": float(finest["after"]["max"]),
                "finest_force_reduction_fraction": float(finest["rms_reduction_fraction"]),
            }
        )

    # Endpoint references are descriptive only; they do not select an interior weight.
    quartic_zero = float(rows[0]["finest_zero_force_rms"])
    quartic_forced = float(rows[0]["finest_projected_force_rms"])
    compact_zero = float(rows[-1]["finest_zero_force_rms"])
    compact_forced = float(rows[-1]["finest_projected_force_rms"])
    for row in rows:
        row["vs_quartic_zero_force_fractional_change"] = _fractional_change(
            float(row["finest_zero_force_rms"]), quartic_zero
        )
        row["vs_quartic_projected_force_fractional_change"] = _fractional_change(
            float(row["finest_projected_force_rms"]), quartic_forced
        )
        row["vs_compact_zero_force_fractional_change"] = _fractional_change(
            float(row["finest_zero_force_rms"]), compact_zero
        )
        row["vs_compact_projected_force_fractional_change"] = _fractional_change(
            float(row["finest_projected_force_rms"]), compact_forced
        )

    # Direct final-deliverable checks use an interior weight but do not promote it.
    midpoint_candidate = candidates[FROZEN_BLEND_WEIGHTS.index(0.5)]
    roundtripped = _candidate_roundtrip(midpoint_candidate)
    if roundtripped.sha256 != midpoint_candidate.sha256:
        raise AssertionError("blend candidate serialization changed identity")
    axes = np.asarray([-1.0, 0.0, 1.0], dtype=float)
    times = np.asarray([0.25, 0.375, 0.5], dtype=float)
    grid = roundtripped.grid(axes, axes, axes, times)
    expected_shape = (len(times), len(axes), len(axes), len(axes), 3)
    if grid.shape != expected_shape or not np.all(np.isfinite(grid)):
        raise AssertionError("blend candidate failed direct finite [u,v,w] grid export")
    if not bool(np.any(np.abs(grid) > 0.0)):
        raise AssertionError("blend candidate exported an identically zero velocity grid")

    return {
        "schema": "eq45_supported_phi10_compact_quartic_blend_force_screen_v1",
        "task_id": TASK_ID,
        "claim_scope": "bounded_restricted_force_tradeoff_curve_on_preregistered_blend_grid",
        "velocity_schedule_frozen_before_pde_fit": True,
        "blend_weights": [float(value) for value in FROZEN_BLEND_WEIGHTS],
        "blend_weight_selected": False,
        "blend_weight_selection_rule": None,
        "early_delta": FROZEN_EARLY_DELTA,
        "mode": {"family": "phi", "index": [1, 0]},
        "base_supported_sha256": base.sha256,
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
        "public_image_fitted": False,
        "rows": rows,
        "endpoint_finest_holdout": {
            "quartic_zero_force_rms": quartic_zero,
            "quartic_projected_force_rms": quartic_forced,
            "compact_zero_force_rms": compact_zero,
            "compact_projected_force_rms": compact_forced,
            "compact_vs_quartic_zero_force_fractional_change": _fractional_change(
                compact_zero, quartic_zero
            ),
            "compact_vs_quartic_projected_force_fractional_change": _fractional_change(
                compact_forced, quartic_forced
            ),
        },
        "direct_velocity_export_check": {
            "checked_blend_weight": 0.5,
            "candidate_sha256": midpoint_candidate.sha256,
            "serialization_roundtrip_exact": True,
            "grid_shape": list(grid.shape),
            "grid_finite": True,
            "grid_nonzero": True,
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
    report = screen_compact_quartic_blend_force()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
