from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_oscillatory_admission_routing_checkpoint import (
    AGENT2_HEAD,
    AGENT4_ARTIFACT_ID,
    AGENT4_HEAD,
    FORMAL_GATES,
    build_checkpoint,
    load_checkpoint,
    save_checkpoint,
    validate_checkpoint,
)


def test_v40_records_only_independent_oscillatory_component_promotion() -> None:
    payload = build_checkpoint()
    states = payload["states"]

    assert payload["schema"] == "kokuno-agent5-oscillatory-admission-routing-checkpoint-v40"
    assert states["leading_ready"] is False
    assert states["public_provenance_labelled_xyz_t_oscillatory_velocity_ready"] is True
    assert states["public_oscillatory_project_support_passed"] is True
    assert states["public_oscillatory_covariance_rank_passed"] is True
    assert states["public_oscillatory_local_divergence_passed"] is True
    assert states["public_oscillatory_preflight_passed"] is True
    assert states["oscillatory_ready"] is True

    assert states["correction_receipt_identity_pinned_to_current_pass"] is False
    assert states["correction_ingest_allowed"] is False
    assert states["same_cycle_requested_stress_materialized"] is False
    assert states["candidate_numeric_finite_head_mean_debt_materialized"] is False
    assert states["correction_ready"] is False
    assert states["candidate_artifact_instantiated"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False

    osc = payload["typed_handoffs"]["oscillatory"]
    assert osc["candidate_head"] == AGENT2_HEAD
    assert osc["validator_head"] == AGENT4_HEAD
    assert osc["validator_artifact_id"] == AGENT4_ARTIFACT_ID
    assert osc["admitted"] is True
    assert osc["scope"] == "oscillatory_component_preflight_only"
    assert payload["formal_gates"] == FORMAL_GATES


def test_v40_frozen_metrics_satisfy_unchanged_local_guards() -> None:
    payload = build_checkpoint()
    audit = payload["upstream"]["agent4_public_z_independent_audit"]
    passes = audit["derived_guard_passes"]
    assert all(passes.values())
    assert audit["scientific_guards_changed_from_agent4_543"] is False
    assert audit["scientific_verdict"] == "PASS_local_oscillatory_preflight_only"

    metrics = audit["metrics"]
    assert metrics["relative_divergence_rms_by_step"]["0.005"] < 2.0e-5
    assert metrics["finest_relative_divergence_max"] < 1.0e-4
    assert min(metrics["relative_rms_refinement_ratios"]) > 3.0
    assert min(metrics["covariance_rank_ratio_by_phase_resolution"].values()) > 2.0e-2
    assert metrics["project_axial_exterior_absolute_max"] == 0.0


def test_v40_roundtrip(tmp_path) -> None:
    path = tmp_path / "checkpoint.json"
    written = save_checkpoint(path)
    loaded = load_checkpoint(path)
    assert loaded == written
    assert loaded["checkpoint_sha256"] == written["checkpoint_sha256"]


def test_v40_rejects_pde_promotion() -> None:
    payload = build_checkpoint()
    mutated = copy.deepcopy(payload)
    mutated["states"]["pde_validated"] = True
    mutated.pop("checkpoint_sha256")
    from openai_ns_reconstruction.kokuno_oscillatory_admission_routing_checkpoint import _canonical_sha256

    mutated["checkpoint_sha256"] = _canonical_sha256(mutated)
    with pytest.raises(ValueError, match="downstream state must remain closed: pde_validated"):
        validate_checkpoint(mutated)


def test_v40_rejects_correction_ingest_before_current_pass_identity_pin() -> None:
    payload = build_checkpoint()
    mutated = copy.deepcopy(payload)
    mutated["states"]["correction_ingest_allowed"] = True
    mutated.pop("checkpoint_sha256")
    from openai_ns_reconstruction.kokuno_oscillatory_admission_routing_checkpoint import _canonical_sha256

    mutated["checkpoint_sha256"] = _canonical_sha256(mutated)
    with pytest.raises(ValueError, match="correction ingest cannot open"):
        validate_checkpoint(mutated)


def test_v40_rejects_validator_identity_laundering() -> None:
    payload = build_checkpoint()
    mutated = copy.deepcopy(payload)
    mutated["upstream"]["agent4_public_z_independent_audit"]["audited_agent2_head"] = "0" * 40
    mutated.pop("checkpoint_sha256")
    from openai_ns_reconstruction.kokuno_oscillatory_admission_routing_checkpoint import _canonical_sha256

    mutated["checkpoint_sha256"] = _canonical_sha256(mutated)
    with pytest.raises(ValueError, match="oscillatory audit identity mismatch"):
        validate_checkpoint(mutated)


def test_v40_truth_boundary_keeps_component_pass_separate_from_full_ns() -> None:
    payload = build_checkpoint()
    truth = payload["truth_boundary"]
    assert truth["oscillatory_component_preflight_only"] is True
    assert truth["kokuno_hidden_background_recovered"] is False
    assert truth["paper_exact"] is False
    assert truth["free_forcing_used_to_cancel_residual"] is False
    assert truth["thresholds_relaxed_after_result"] is False
    assert truth["agent4_component_preflight_counted_as_full_pde_validation"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
