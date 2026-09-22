"""Audit whether the frozen morphology times alias the time-skewed basis probe.

Preregistered in issue #1165 and stacked on Agent-7 PR #1156. This module is
analysis-only: it changes no candidate, coefficient, basis, pressure, forcing,
or scientific readiness state.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

import agent7_st052m_representation_archetype_left_null_preflight as rep

TASK_ID = "CR003-ST052M-TEMPORAL-SAMPLING-ALIAS-AUDIT-138"
PREREG_ISSUE = 1165
STACK_BASE_PR = 1156
STACK_BASE_HEAD = "eb65d27b0b7bfb12fbb8010858f1a116902b84be"

OBSERVATION_TIMES = (0.25, 0.50, 0.75)
OFFGRID_WITNESS_TIMES = (0.375, 0.625)
ZERO_TOL = 1.0e-15
OFFGRID_NONZERO_FLOOR = 1.0e-6

TRUTH = {
    "candidate_velocity_changed": False,
    "basis_dimension_changed": False,
    "coefficient_selected": False,
    "pressure_or_force_changed": False,
    "held_out_pde_residual_evaluated": False,
    "public_image_numeric_target_used": False,
    "renderer_or_camera_fit_used": False,
    "direct_visualization_fingerprint_improvement": 0.0,
    "closer_visualization_delivery_established": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _scientific_parent():
    return rep._scientific_parent()


def _assert_source_lock() -> None:
    rep._assert_source_lock()
    if rep.TASK_ID != "CR003-ST052M-REPRESENTATION-ARCHETYPE-LEFT-NULL-PREFLIGHT-136":
        raise RuntimeError("#1147 representation parent identity drifted")
    if rep.PREREG_ISSUE != 1146 or rep.STACK_BASE_PR != 1139:
        raise RuntimeError("#1147 preregistration/lineage drifted")
    scientific = _scientific_parent()
    if scientific.TASK_ID != "CR003-ST052M-DEFECT-COORDINATE-CONTROLLABILITY-128":
        raise RuntimeError("#1078 scientific parent identity drifted")
    expected = (
        "inward_radial_index_t050",
        "swirl_to_poloidal_ratio_t050",
        "axial_vorticity_aspect_t050",
        "core_width_t050",
        "tip_radial_thickness_t050",
        "angular_rotation_radial_variation_t050",
        "core_width_change_025_to_075",
        "core_speed_change_025_to_075",
    )
    if tuple(scientific.DEFECT_COORDINATES) != expected:
        raise RuntimeError("#1078 defect-coordinate time support drifted")


def _sample_points() -> np.ndarray:
    # Interior points where the frozen compact poloidal unit is nonzero.
    return np.asarray(
        [
            [0.25, 0.10, 0.40],
            [-0.40, 0.20, -0.55],
            [0.60, -0.25, 0.75],
        ],
        dtype=float,
    )


def audit_activation(activation: Callable[[float], float]) -> dict[str, Any]:
    observed = {f"{t:.3f}": float(activation(t)) for t in OBSERVATION_TIMES}
    witness = {f"{t:.3f}": float(activation(t)) for t in OFFGRID_WITNESS_TIMES}
    obs_zero = max(abs(v) for v in observed.values()) <= ZERO_TOL
    witness_nonzero = min(abs(v) for v in witness.values()) >= OFFGRID_NONZERO_FLOOR
    opposite_sign = witness["0.375"] * witness["0.625"] < 0.0
    return {
        "observation_activation": observed,
        "offgrid_witness_activation": witness,
        "zero_at_all_observation_times": bool(obs_zero),
        "offgrid_nonzero": bool(witness_nonzero),
        "offgrid_opposite_sign": bool(opposite_sign),
        "temporal_sampling_alias_detected": bool(obs_zero and witness_nonzero and opposite_sign),
    }


def build_report() -> dict[str, Any]:
    _assert_source_lock()
    scalar = audit_activation(rep._activation_skew)
    tangent = rep._probe_tangents()["time_skewed_poloidal_extension"]
    points = _sample_points()

    field_obs = {}
    field_witness = {}
    for t in OBSERVATION_TIMES:
        field_obs[f"{t:.3f}"] = float(np.max(np.abs(tangent(points, t))))
    for t in OFFGRID_WITNESS_TIMES:
        field_witness[f"{t:.3f}"] = float(np.max(np.abs(tangent(points, t))))

    field_zero = max(field_obs.values()) <= ZERO_TOL
    field_nonzero = min(field_witness.values()) >= OFFGRID_NONZERO_FLOOR
    alias = bool(scalar["temporal_sampling_alias_detected"] and field_zero and field_nonzero)

    consequence = (
        "exclude_time_skewed_probe_from_zero_derivative_capacity_ranking_until_offgrid_temporal_morphology_coordinate_exists"
        if alias
        else "no_alias_established"
    )
    return {
        "task_id": TASK_ID,
        "preregistered_issue": PREREG_ISSUE,
        "stack_base_pr": STACK_BASE_PR,
        "stack_base_head": STACK_BASE_HEAD,
        "frozen_observation_times": list(OBSERVATION_TIMES),
        "frozen_offgrid_witness_times": list(OFFGRID_WITNESS_TIMES),
        "activation_audit": scalar,
        "field_tangent_max_abs_at_observation_times": field_obs,
        "field_tangent_max_abs_at_offgrid_witness_times": field_witness,
        "field_zero_at_observation_times": bool(field_zero),
        "field_nonzero_offgrid": bool(field_nonzero),
        "temporal_sampling_alias_detected": alias,
        "routing_consequence": consequence,
        "minimum_followup": (
            "add_one_preregistered_offgrid temporal morphology coordinate before judging time-dependent coefficient capacity"
            if alias
            else "none"
        ),
        "truth": dict(TRUTH),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_report()
    text = json.dumps(report, indent=2, sort_keys=True)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
