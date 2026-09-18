from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_source_wave_shell_pass_routing_checkpoint import (
    AGENT4_PASS_RECEIPT,
    FIXED_GATES,
    HISTORICAL_AGENT4_REJECTION,
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


def test_checkpoint_binds_new_agent4_pass_without_rewriting_old_rejection() -> None:
    checkpoint = build_checkpoint()
    assert validate_checkpoint(checkpoint) == checkpoint
    assert checkpoint["agent4_pass_receipt"] == AGENT4_PASS_RECEIPT
    assert checkpoint["historical_agent4_rejection"] == HISTORICAL_AGENT4_REJECTION
    assert checkpoint["agent4_pass_receipt"]["source_pr"] == 383
    assert checkpoint["agent4_pass_receipt"]["source_wave_shell_roundoff_ladder_preflight_passed"] is True
    assert checkpoint["agent4_pass_receipt"]["finest"]["step"] == pytest.approx(1.0e-3)
    assert checkpoint["agent4_pass_receipt"]["finest"]["curl_relative_rms"] == pytest.approx(2.755788271978995e-10)
    assert checkpoint["agent4_pass_receipt"]["finest"]["divergence_max_abs"] == pytest.approx(3.580782588711177e-09)
    assert checkpoint["agent4_pass_receipt"]["divergence_refinement_ratios"] == pytest.approx(
        [15.995372546512321, 16.024253341545737]
    )
    assert checkpoint["historical_agent4_rejection"]["source_wave_shell_independent_cartesian_preflight_passed"] is False
    assert checkpoint["historical_agent4_rejection"]["failed_checks"] == ["divergence_refinement_ratio"]
    assert checkpoint["historical_agent4_rejection"]["overwritten"] is False


def test_latest_sibling_handoffs_are_bound_but_not_laundered() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["upstream"] == UPSTREAM
    a2 = checkpoint["upstream"]["agent2"]
    assert a2["latest_pr"] == 380
    assert a2["dedicated_run"] == 35335492487
    assert a2["standard_run"] == 35335492324
    assert a2["localized_real_pair_family_ready"] is True
    assert a2["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert a2["genuinely_independent_second_covariance_column_ready"] is False
    assert a2["consumed_in_executable_ancestry"] is False
    a3 = checkpoint["upstream"]["agent3"]
    assert a3["latest_pr"] == 382
    assert a3["dedicated_run"] == 35336161949
    assert a3["artifact_id"] == 10542888440
    assert a3["family_cross_term_contract_executable"] is True
    assert a3["bounded_inverse_ready"] is False
    assert a3["finite_correction_cycle_rerun_allowed"] is False
    assert a3["consumed_in_executable_ancestry"] is False


def test_pipeline_frontier_clears_only_local_preflight() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["states"] == STATES
    assert checkpoint["pipeline_frontier"] == PIPELINE_FRONTIER
    assert checkpoint["states"]["source_wave_shell_independent_cartesian_preflight_passed"] is True
    assert checkpoint["states"]["historical_source_wave_shell_rejection_preserved"] is True
    assert checkpoint["states"]["oscillatory_machinery_ready"] is True
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["states"]["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert checkpoint["states"]["genuinely_independent_second_covariance_column_ready"] is False
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["candidate_artifact_instantiated"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert all(stage["ready"] is False for stage in checkpoint["pipeline_frontier"])
    osc = next(stage for stage in checkpoint["pipeline_frontier"] if stage["stage"] == "oscillatory_augmentation")
    assert osc["status"] == "source_shell_preflight_cleared_waiting_for_actual_source_inputs"


def test_checkpoint_rejects_rewriting_historical_failure() -> None:
    mutated = copy.deepcopy(build_checkpoint())
    mutated["historical_agent4_rejection"]["source_wave_shell_independent_cartesian_preflight_passed"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="historical Agent-4 rejection changed|historical source-shell rejection was overwritten"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_global_or_pde_promotion() -> None:
    for key in (
        "public_source_oscillatory_xyz_t_velocity_ready",
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
    with pytest.raises(ValueError, match="upstream routing receipt changed|finite correction cycle was falsely authorized"):
        validate_checkpoint(mutated)


def test_bundle_roundtrip_retains_pass_and_sibling_provenance(tmp_path) -> None:
    paths = write_bundle(tmp_path / "bundle")
    checkpoint = json.loads(paths["checkpoint"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    validate_checkpoint(checkpoint)
    assert manifest["checkpoint_payload_sha256"] == checkpoint["checkpoint_sha256"]
    assert manifest["bound_agent4_artifact_id"] == 10542893497
    assert manifest["bound_agent4_artifact_digest"] == "sha256:4350fbeba26ca8a54ee08386f48c33709af034fd96df443eb4a9ec904e1e91a1"
    assert manifest["bound_agent3_sibling_artifact_id"] == 10542888440
    assert manifest["historical_agent4_rejection_preserved"] is True
    assert manifest["pde_validated"] is False
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False
