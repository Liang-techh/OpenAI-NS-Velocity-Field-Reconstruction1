from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_source_compatible_partition_routing_checkpoint import (
    AGENT4_PARTITION_AUDIT_RECEIPT,
    AUTONOMOUS_ROUTE_POLICY,
    FIXED_GATES,
    PIPELINE_FRONTIER,
    STATES,
    UPSTREAM,
    _sha,
    build_checkpoint,
    validate_checkpoint,
    write_bundle,
)


def _resign(payload: dict) -> None:
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    payload["checkpoint_sha256"] = _sha(unsigned)


def test_checkpoint_binds_independent_partition_pass() -> None:
    checkpoint = build_checkpoint()
    assert validate_checkpoint(checkpoint) == checkpoint
    receipt = checkpoint["agent4_partition_audit_receipt"]
    assert receipt == AGENT4_PARTITION_AUDIT_RECEIPT
    assert receipt["source_pr"] == 394
    assert receipt["all_local_checks_passed"] is True
    assert receipt["worst_partition_max_abs"] == pytest.approx(4.440892098500626e-16)
    assert receipt["worst_differentiated_closure_max_abs"] == pytest.approx(1.3877787807814457e-16)
    assert receipt["maximum_active_ell_span"] == 1
    assert max(receipt["finest_radial_relative_rms_by_origin"]) < 1.0e-6
    assert max(receipt["finest_axial_relative_rms_by_origin"]) < 1.0e-6
    assert min(min(row) for row in receipt["radial_refinement_ratios_by_origin"]) >= 8.0
    assert min(min(row) for row in receipt["axial_refinement_ratios_by_origin"]) >= 8.0
    assert receipt["formal_full_domain_pde_gate_assessed"] is False


def test_autonomous_scaffold_truth_boundary_is_explicit() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["autonomous_route_policy"] == AUTONOMOUS_ROUTE_POLICY
    policy = checkpoint["autonomous_route_policy"]
    assert policy["kokuno_derived_autonomous_choices_allowed"] is True
    assert policy["source_compatible_partition_may_be_used_as_scaffold"] is True
    assert policy["source_compatible_partition_may_be_called_recovered_source"] is False
    assert policy["autonomous_agent1_loop_may_be_used_as_kokuno_derived_scaffold"] is True
    assert policy["autonomous_agent1_loop_may_be_called_source_hidden_loop"] is False
    assert policy["paper_exact"] is False

    a2 = checkpoint["upstream"]["agent2"]
    assert a2["latest_pr"] == 392
    assert a2["source_compatible_partition_realization_ready"] is True
    assert a2["partition_independently_audited"] is True
    assert a2["partition_is_repository_autonomous"] is True
    assert a2["source_actual_partition_recovered"] is False
    assert a2["public_xyz_t_by_beta_velocity_ready"] is False


def test_latest_agent1_and_agent3_siblings_are_bound_without_laundering() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["upstream"] == UPSTREAM
    a1 = checkpoint["upstream"]["agent1"]
    assert a1["latest_pr"] == 391
    assert a1["autonomous_diagnostic_loop_available"] is True
    assert a1["source_admissible_loop_reconstructed"] is False
    assert a1["global_leading_profile_reconstructed"] is False
    assert a1["consumed_in_executable_ancestry"] is False

    a3 = checkpoint["upstream"]["agent3"]
    assert a3["latest_pr"] == 393
    assert a3["family_bounded_inverse_machinery_ready"] is True
    assert a3["nodes_requiring_second_direction"] == 25
    assert a3["duplicate_control_rank_two_coverage"] == "0/25"
    assert a3["missing_relative_stress_response"] == pytest.approx(0.210194502136004)
    assert a3["finite_correction_cycle_rerun_allowed"] is False
    assert a3["consumed_in_executable_ancestry"] is False


def test_pipeline_frontier_clears_only_partition_seam() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["states"] == STATES
    assert checkpoint["pipeline_frontier"] == PIPELINE_FRONTIER
    assert checkpoint["states"]["source_compatible_partition_realization_ready"] is True
    assert checkpoint["states"]["source_compatible_partition_independently_audited"] is True
    assert checkpoint["states"]["source_actual_partition_recovered"] is False
    assert checkpoint["states"]["oscillatory_machinery_ready"] is True
    assert checkpoint["states"]["oscillatory_ready"] is False
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["candidate_artifact_instantiated"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert all(stage["ready"] is False for stage in checkpoint["pipeline_frontier"])
    ingest = next(stage for stage in checkpoint["pipeline_frontier"] if stage["stage"] == "source_profile_ingest")
    assert ingest["status"] == "source_compatible_autonomous_partition_independently_cleared"


def test_checkpoint_rejects_source_recovery_or_global_promotion() -> None:
    mutated = copy.deepcopy(build_checkpoint())
    mutated["upstream"]["agent2"]["source_actual_partition_recovered"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="upstream routing receipt changed|autonomous partition was promoted"):
        validate_checkpoint(mutated)

    for key in (
        "source_actual_partition_recovered",
        "leading_ready",
        "oscillatory_ready",
        "public_kokuno_derived_oscillatory_xyz_t_velocity_ready",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
    ):
        mutated = copy.deepcopy(build_checkpoint())
        mutated["states"][key] = True
        _resign(mutated)
        with pytest.raises(ValueError, match="routing states changed|global readiness was falsely promoted"):
            validate_checkpoint(mutated)


def test_checkpoint_rejects_threshold_or_finite_cycle_change() -> None:
    mutated = copy.deepcopy(build_checkpoint())
    mutated["fixed_gates"]["held_out_normalized_full_momentum_max"] = 2.0e-3
    _resign(mutated)
    with pytest.raises(ValueError, match="fixed gates changed"):
        validate_checkpoint(mutated)

    mutated = copy.deepcopy(build_checkpoint())
    mutated["upstream"]["agent3"]["finite_correction_cycle_rerun_allowed"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="upstream routing receipt changed|finite correction cycle"):
        validate_checkpoint(mutated)


def test_bundle_roundtrip_retains_receipts_and_truth_boundary(tmp_path) -> None:
    paths = write_bundle(tmp_path / "bundle")
    checkpoint = json.loads(paths["checkpoint"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    validate_checkpoint(checkpoint)
    assert manifest["checkpoint_payload_sha256"] == checkpoint["checkpoint_sha256"]
    assert manifest["bound_agent4_artifact_id"] == 10545247621
    assert manifest["bound_agent4_artifact_digest"] == "sha256:6b2cc78a2f0ba6977b938fdafb6aa62a4059964727e14ad9fa681c7246c09765"
    assert manifest["bound_agent3_sibling_artifact_id"] == 10545117446
    assert manifest["source_compatible_partition_independently_audited"] is True
    assert manifest["source_actual_partition_recovered"] is False
    assert manifest["candidate_artifact_instantiated"] is False
    assert manifest["pde_validated"] is False
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False
