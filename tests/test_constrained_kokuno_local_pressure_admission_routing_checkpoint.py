from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_local_pressure_admission_routing_checkpoint import (
    FORMAL_GATES,
    build_checkpoint,
    validate_checkpoint,
)


def test_v48_admits_only_independent_local_pressure_preflight() -> None:
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)

    states = checkpoint["states"]
    assert states["selected_pressure_local_coordinate_independent_preflight_passed"] is True
    assert states["selected_pressure_local_coordinate_independently_admitted"] is True
    assert states["leading_ready"] is False
    assert states["source_pressure_radius_one_ball_bounds_ready"] is False
    assert states["source_R1_R2_MK_ready"] is False
    assert states["global_matched_pressure_ready"] is False
    assert states["global_leading_velocity_ready"] is False
    assert states["correction_ready"] is False
    assert states["candidate_artifact_instantiated"] is False
    assert states["pde_validated"] is False

    metrics = checkpoint["independent_local_pressure_metrics"]
    assert metrics["fd8_finest_relative_rms"] <= 5.0e-6
    assert metrics["fd8_finest_relative_max"] <= 1.0e-5
    assert min(metrics["fd8_refinement_ratios"]) >= 8.0
    assert metrics["fundamental_theorem_finest_relative_rms"] <= 1.0e-7
    assert metrics["fundamental_theorem_finest_relative_max"] <= 5.0e-7
    assert min(metrics["fundamental_theorem_refinement_ratios"]) >= 4.0
    assert metrics["physical_chain_rule_relative_rms"] <= 5.0e-6
    assert metrics["physical_chain_rule_relative_max"] <= 1.0e-5
    assert metrics["mutation_sign_flip_relative_rms"] >= 0.5
    assert checkpoint["formal_gates_unchanged"] == FORMAL_GATES


def test_v48_keeps_unresolved_newer_siblings_fail_closed() -> None:
    checkpoint = build_checkpoint()
    newer = checkpoint["upstream"]["newer_sibling_work"]
    assert newer["agent1_pressure_ball"]["admitted_by_this_checkpoint"] is False
    assert newer["agent3_finite_cycle"]["admitted_by_this_checkpoint"] is False
    assert newer["agent4_pressure_ball_audit"]["resolved_at_checkpoint_freeze"] is False
    assert newer["agent4_pressure_ball_audit"]["admitted_by_this_checkpoint"] is False


def test_v48_rejects_premature_leading_or_pde_promotion() -> None:
    for key in ("leading_ready", "global_matched_pressure_ready", "pde_validated"):
        checkpoint = build_checkpoint()
        checkpoint["states"][key] = True
        with pytest.raises(ValueError, match="premature downstream promotion"):
            validate_checkpoint(checkpoint)


def test_v48_rejects_local_pressure_pass_demotion() -> None:
    checkpoint = build_checkpoint()
    checkpoint["states"]["selected_pressure_local_coordinate_independently_admitted"] = False
    with pytest.raises(ValueError, match="required admitted state"):
        validate_checkpoint(checkpoint)


def test_v48_rejects_gate_or_hash_mutation() -> None:
    checkpoint = build_checkpoint()
    bad_gate = copy.deepcopy(checkpoint)
    bad_gate["formal_gates_unchanged"]["held_out_normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError, match="checkpoint SHA256 mismatch|formal PDE gates changed"):
        validate_checkpoint(bad_gate)

    bad_hash = copy.deepcopy(checkpoint)
    bad_hash["checkpoint_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="checkpoint SHA256 mismatch"):
        validate_checkpoint(bad_hash)


def test_v48_rejects_unresolved_pressure_ball_promotion() -> None:
    checkpoint = build_checkpoint()
    checkpoint["upstream"]["newer_sibling_work"]["agent4_pressure_ball_audit"][
        "resolved_at_checkpoint_freeze"
    ] = True
    with pytest.raises(ValueError, match="checkpoint SHA256 mismatch|unresolved pressure-ball audit"):
        validate_checkpoint(checkpoint)
