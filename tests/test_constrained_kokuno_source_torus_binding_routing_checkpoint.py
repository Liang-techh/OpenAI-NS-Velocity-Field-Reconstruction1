from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.kokuno_source_torus_binding_independent import generate_report
from openai_ns_reconstruction.kokuno_source_torus_binding_routing_checkpoint import (
    AGENT4_RECEIPT,
    SCHEMA,
    ST006_BASELINE,
    TASK_ID,
    build_checkpoint,
    checkpoint_sha256,
    validate_checkpoint,
)


def test_checkpoint_identity_and_hash_are_deterministic() -> None:
    first = build_checkpoint()
    second = build_checkpoint()
    assert first == second
    assert first["schema"] == SCHEMA
    assert first["task_id"] == TASK_ID
    assert first["checkpoint_sha256"] == checkpoint_sha256(first)
    validate_checkpoint(first)


def test_displayed_source_torus_binding_is_bound_without_pde_promotion() -> None:
    checkpoint = build_checkpoint()
    states = checkpoint["states"]
    a2 = checkpoint["upstream"]["agent2"]
    a4 = checkpoint["upstream"]["agent4"]
    assert states["displayed_source_torus_constants_bound"] is True
    assert states["displayed_source_torus_binding_independently_audited"] is True
    assert a2["displayed_source_torus_constants_bound"] is True
    assert a4["base_agent2_head"] == a2["head"]
    assert a4["local_structural_preflight_passed"] is True
    assert a4["formal_full_domain_pde_gate_assessed"] is False
    assert checkpoint["states"]["pde_validated"] is False


def test_agent4_public_report_still_passes_frozen_source_binding_guards() -> None:
    report = generate_report()
    summary = report["summary"]
    guards = AGENT4_RECEIPT["local_guards"]
    assert report["local_structural_preflight_passed"] is True
    assert summary["worst_finest_primary_relative_rms"] <= guards["finest_relative_rms"]
    assert summary["worst_phase_embedding_max_abs_error"] <= guards["phase_embedding_max_abs"]
    assert summary["worst_source_binding_relative_error"] <= guards["source_binding_relative"]
    assert summary["worst_refinement_ratio"] >= guards["minimum_refinement_ratio"]
    assert summary["weakest_wrong_sign_mutation_relative_rms"] >= guards[
        "wrong_sign_mutation_relative_rms_floor"
    ]
    assert summary["weakest_exponent_mutation_relative_rms"] >= guards[
        "exponent_mutation_relative_rms_floor"
    ]
    assert report["formal_gates"]["formal_full_domain_pde_gate_assessed"] is False
    assert report["formal_gates"]["pde_validated"] is False


def test_missing_scientific_objects_remain_fail_closed() -> None:
    states = build_checkpoint()["states"]
    assert states["oscillatory_machinery_ready"] is True
    assert states["displayed_source_torus_constants_bound"] is True
    assert states["displayed_source_torus_binding_independently_audited"] is True
    assert states["shared_spacetime_quadratic_damping_ready"] is True
    for key in (
        "leading_ready",
        "actual_positive_order_background_bound",
        "actual_auxiliary_torus_mode_family_bound",
        "public_q_scaled_by_beta_total_velocity_ready",
        "oscillatory_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
    ):
        assert states[key] is False


def test_st006_and_formal_gates_are_unchanged_and_noncomparable() -> None:
    checkpoint = build_checkpoint()
    baseline = checkpoint["baseline_vs_kokuno"]["st006"]
    assert baseline == ST006_BASELINE
    assert baseline["momentum_sampled_max"] == pytest.approx(0.1082289305112118, rel=0, abs=0)
    assert baseline["momentum_volume_l2"] == pytest.approx(0.10758432876230622, rel=0, abs=0)
    assert checkpoint["formal_gates"]["held_out_normalized_momentum_max"] == 1.0e-3
    assert checkpoint["formal_gates"]["held_out_divergence_max"] == 1.0e-5
    assert checkpoint["formal_gates"]["changed_this_round"] is False
    assert baseline["comparison_to_agent4_local_source_binding_audit"] == "not directly comparable"


def test_truth_boundary_tampering_is_rejected() -> None:
    checkpoint = build_checkpoint()

    promoted = deepcopy(checkpoint)
    promoted["states"]["pde_validated"] = True
    promoted["checkpoint_sha256"] = checkpoint_sha256(promoted)
    with pytest.raises(ValueError, match="pde_validated"):
        validate_checkpoint(promoted)

    erased_audit = deepcopy(checkpoint)
    erased_audit["states"]["displayed_source_torus_binding_independently_audited"] = False
    erased_audit["checkpoint_sha256"] = checkpoint_sha256(erased_audit)
    with pytest.raises(ValueError, match="independently audited"):
        validate_checkpoint(erased_audit)

    relaxed = deepcopy(checkpoint)
    relaxed["formal_gates"]["held_out_normalized_momentum_max"] = 2.0e-3
    relaxed["checkpoint_sha256"] = checkpoint_sha256(relaxed)
    with pytest.raises(ValueError, match="momentum gate"):
        validate_checkpoint(relaxed)

    instantiated = deepcopy(checkpoint)
    instantiated["candidate_artifact_contract"]["instantiated"] = True
    instantiated["checkpoint_sha256"] = checkpoint_sha256(instantiated)
    with pytest.raises(ValueError, match="candidate artifact"):
        validate_checkpoint(instantiated)
