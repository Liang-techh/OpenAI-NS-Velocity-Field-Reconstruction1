"""Restricted-force PDE-side screen for the bounded axial-taper control.

Agent 7 PR #193 isolated the extreme top/bottom support collar as a representation
bottleneck, and Agent 1 PR #194 materialized the smallest available control:
``AxisymmetricPhysicalTaper.axial_plateau_q`` in ``[0.64, 0.81]``.  This CR005
increment freezes one representative compact--quartic blend weight and a small
predeclared axial-plateau grid *before* any PDE calculation.  It then reuses the
already-checked support-connected vorticity-equation operator and projects only the
preregistered bounded ``RestrictedForce(a,c)`` family.

The purpose is routing: quantify the PDE-side cost/benefit of moving the axial
identity plateau while keeping the final public ``[u,v,w]`` directly evaluable.
No taper value is selected here.  Pressure is not fitted, no force direction is
added, and a sampled vorticity residual is not a formal PDE acceptance test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from .constrained_eq45_supported_collar_vorticity import governed_supported_seed
from .constrained_eq45_supported_force_projection import (
    FIT_MAX_ITER,
    FIT_TOL,
    FORCE_BOUNDS,
)
from .constrained_eq45_supported_phi10_blend_axial_taper_mode import (
    Eq45SupportedPhi10BlendAxialTaperCandidate,
)
from .constrained_eq45_supported_phi10_compact_quartic_blend_temporal_mode import (
    Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate,
)
from .constrained_eq45_supported_phi10_force_crosscheck import _project_temporal_trial

TASK_ID = "CR005-EQ45-SUPPORTED-BLEND-AXIAL-TAPER-FORCE-SCREEN-035"
FROZEN_BLEND_WEIGHT = 0.5
FROZEN_EARLY_DELTA = -1.4
FROZEN_IDENTITY_HALF_HEIGHTS = (1.60, 1.65, 1.70, 1.75, 1.80)
FROZEN_AXIAL_PLATEAU_Q = tuple((height / 2.0) ** 2 for height in FROZEN_IDENTITY_HALF_HEIGHTS)


def _fractional_change(new: float, old: float) -> float:
    return float((new - old) / max(abs(old), np.finfo(float).tiny))


def _roundtrip(candidate):
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "candidate.json"
        candidate.save_json(path)
        return type(candidate).load_json(path)


def screen_blend_axial_taper_force(
    *, constraints_path: str | Path | None = None
) -> dict:
    """Screen the preregistered axial-taper grid under one force contract."""
    root = Path(__file__).resolve().parents[2]
    constraints_path = Path(constraints_path or root / "configs/constraints.json")
    constraints = json.loads(constraints_path.read_text(encoding="utf-8"))
    if constraints["forcing"]["mode"] != "restricted_two_parameter_family":
        raise AssertionError("active forcing convention drifted from preregistration")

    supported = governed_supported_seed()
    blend = Eq45SupportedPhi10CompactQuarticBlendTemporalCandidate(
        base=supported,
        blend_weight=FROZEN_BLEND_WEIGHT,
        early_delta=FROZEN_EARLY_DELTA,
    )
    nu = float(constraints["nu"])
    steps = tuple(float(value) for value in constraints["validation"]["derivative_steps"])

    # Candidate identities and the trial grid are fixed before any residual is evaluated.
    candidates = tuple(
        Eq45SupportedPhi10BlendAxialTaperCandidate(
            base=blend,
            axial_plateau_q=q,
        )
        for q in FROZEN_AXIAL_PLATEAU_Q
    )
    frozen_shas = tuple(candidate.sha256 for candidate in candidates)

    axial_response_points = np.asarray(
        [
            [0.35, 0.0, 1.68],
            [0.80, 0.0, 1.75],
            [1.25, 0.0, 1.82],
            [0.35, 0.0, -1.68],
            [0.80, 0.0, -1.75],
            [1.25, 0.0, -1.82],
        ],
        dtype=float,
    )
    response_time = 0.3125
    baseline_velocity = candidates[0].at_points(axial_response_points, response_time)
    baseline_velocity_rms = float(np.sqrt(np.mean(np.sum(baseline_velocity**2, axis=1))))

    rows = []
    for height, q, candidate, frozen_sha in zip(
        FROZEN_IDENTITY_HALF_HEIGHTS,
        FROZEN_AXIAL_PLATEAU_Q,
        candidates,
        frozen_shas,
    ):
        projection = _project_temporal_trial(
            candidate,
            nu=nu,
            derivative_steps=steps,
        )
        if candidate.sha256 != frozen_sha:
            raise AssertionError("axial-taper candidate identity changed during force projection")
        finest = projection["holdout_rows"][-1]
        velocity = candidate.at_points(axial_response_points, response_time)
        velocity_delta_rms = float(
            np.sqrt(np.mean(np.sum((velocity - baseline_velocity) ** 2, axis=1)))
        )
        rows.append(
            {
                "axial_identity_half_height": float(height),
                "axial_plateau_q": float(q),
                "candidate_sha256": frozen_sha,
                "fit": projection["fit"],
                "holdout_rows": projection["holdout_rows"],
                "finest_zero_force_rms": float(finest["before"]["rms"]),
                "finest_projected_force_rms": float(finest["after"]["rms"]),
                "finest_zero_force_max": float(finest["before"]["max"]),
                "finest_projected_force_max": float(finest["after"]["max"]),
                "finest_force_reduction_fraction": float(finest["rms_reduction_fraction"]),
                "finest_by_region": finest["by_region"],
                "public_axial_probe_velocity_delta_rms_vs_q064": velocity_delta_rms,
                "public_axial_probe_velocity_relative_delta_vs_q064": velocity_delta_rms
                / max(baseline_velocity_rms, np.finfo(float).tiny),
            }
        )

    base_zero = float(rows[0]["finest_zero_force_rms"])
    base_forced = float(rows[0]["finest_projected_force_rms"])
    base_regions = rows[0]["finest_by_region"]
    for row in rows:
        row["vs_q064_zero_force_fractional_change"] = _fractional_change(
            float(row["finest_zero_force_rms"]), base_zero
        )
        row["vs_q064_projected_force_fractional_change"] = _fractional_change(
            float(row["finest_projected_force_rms"]), base_forced
        )
        region_changes = {}
        for region, metrics in row["finest_by_region"].items():
            region_changes[region] = {
                "zero_force_rms_fractional_change": _fractional_change(
                    float(metrics["before"]["rms"]),
                    float(base_regions[region]["before"]["rms"]),
                ),
                "projected_force_rms_fractional_change": _fractional_change(
                    float(metrics["after"]["rms"]),
                    float(base_regions[region]["after"]["rms"]),
                ),
            }
        row["vs_q064_by_region"] = region_changes

    # One interior candidate must remain a replayable direct [u,v,w] deliverable.
    export_index = 2
    export_candidate = candidates[export_index]
    reloaded = _roundtrip(export_candidate)
    if reloaded.sha256 != export_candidate.sha256:
        raise AssertionError("axial-taper candidate serialization changed identity")
    axes = np.asarray([-1.0, 0.0, 1.0], dtype=float)
    z_axis = np.asarray([-1.75, 0.0, 1.75], dtype=float)
    times = np.asarray([0.25, 0.375, 0.5], dtype=float)
    grid = reloaded.grid(axes, axes, z_axis, times)
    expected_shape = (len(times), len(axes), len(axes), len(z_axis), 3)
    if grid.shape != expected_shape or not np.all(np.isfinite(grid)):
        raise AssertionError("axial-taper candidate failed finite [u,v,w] grid export")
    if not bool(np.any(np.abs(grid) > 0.0)):
        raise AssertionError("axial-taper candidate exported an identically zero velocity grid")

    return {
        "schema": "eq45_supported_blend_axial_taper_force_screen_v1",
        "task_id": TASK_ID,
        "claim_scope": "bounded_restricted_force_tradeoff_curve_on_preregistered_axial_taper_grid",
        "velocity_family_frozen_before_pde_fit": True,
        "blend_weight": FROZEN_BLEND_WEIGHT,
        "blend_weight_selected_by_this_module": False,
        "early_delta": FROZEN_EARLY_DELTA,
        "axial_identity_half_heights": [float(value) for value in FROZEN_IDENTITY_HALF_HEIGHTS],
        "axial_plateau_q_values": [float(value) for value in FROZEN_AXIAL_PLATEAU_Q],
        "axial_taper_value_selected": False,
        "axial_taper_selection_rule": None,
        "base_supported_sha256": supported.sha256,
        "base_blend_sha256": blend.sha256,
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
        "taper_refit_on_pde": False,
        "public_image_fitted": False,
        "rows": rows,
        "baseline_q064_finest_holdout": {
            "zero_force_rms": base_zero,
            "projected_force_rms": base_forced,
        },
        "direct_velocity_export_check": {
            "checked_axial_identity_half_height": float(FROZEN_IDENTITY_HALF_HEIGHTS[export_index]),
            "checked_axial_plateau_q": float(FROZEN_AXIAL_PLATEAU_Q[export_index]),
            "candidate_sha256": export_candidate.sha256,
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
        "production_axial_taper_promoted": False,
        "visual_correspondence_verified": False,
        "visualization_ready": False,
        "physical_support_validated": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = screen_blend_axial_taper_force()
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
