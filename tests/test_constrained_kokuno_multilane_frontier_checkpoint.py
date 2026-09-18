from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_multilane_frontier_checkpoint import (
    CANDIDATE_ARTIFACT_CONTRACT,
    FIXED_GATES,
    LATEST,
    PIPELINE_STAGES,
    STATES,
    _sha,
    build_checkpoint,
    validate_checkpoint,
    write_bundle,
)


def _resign(payload: dict) -> None:
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    payload["checkpoint_sha256"] = _sha(unsigned)


def test_frontier_binds_latest_verified_a1_a2_without_promoting_candidate() -> None:
    checkpoint = build_checkpoint()
    assert validate_checkpoint(checkpoint) == checkpoint
    assert checkpoint["latest"] == LATEST
    assert checkpoint["latest"]["agent1"]["pr"] == 368
    assert checkpoint["latest"]["agent1"]["outer_base_schedule_ready"] is True
    assert checkpoint["latest"]["agent1"]["global_leading_profile_reconstructed"] is False
    assert checkpoint["latest"]["agent2"]["pr"] == 369
    assert checkpoint["latest"]["agent2"]["slow_squared_partition_family_contract_ready"] is True
    assert checkpoint["latest"]["agent2"]["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert checkpoint["states"]["pde_validated"] is False


def test_pipeline_and_candidate_contract_fail_closed() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["pipeline_stages"] == PIPELINE_STAGES
    assert checkpoint["candidate_artifact_contract"] == CANDIDATE_ARTIFACT_CONTRACT
    assert checkpoint["candidate_artifact_contract"]["status"] == "reserved_not_instantiated"
    assert all(stage["ready"] is False for stage in checkpoint["pipeline_stages"])
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False


def test_source_wave_shell_a4_dependency_cannot_be_laundered() -> None:
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    mutated["latest"]["agent2"]["source_wave_shell_independent_cartesian_preflight_passed"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="source-wave-shell independent preflight was falsely promoted"):
        validate_checkpoint(mutated)


def test_finite_cycle_remains_blocked_without_second_column() -> None:
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    mutated["latest"]["agent3"]["finite_correction_cycle_rerun_allowed"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="finite correction cycle was falsely authorized"):
        validate_checkpoint(mutated)


def test_incomplete_pipeline_stage_cannot_be_promoted() -> None:
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    mutated["pipeline_stages"][1]["ready"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="incomplete pipeline stage was promoted"):
        validate_checkpoint(mutated)


def test_fixed_gate_change_is_rejected() -> None:
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    mutated["fixed_gates"]["held_out_normalized_full_momentum_max"] = 2.0e-3
    _resign(mutated)
    with pytest.raises(ValueError, match="fixed gates changed"):
        validate_checkpoint(mutated)


def test_bundle_roundtrip_is_deterministic_and_not_a_candidate(tmp_path) -> None:
    paths = write_bundle(tmp_path / "bundle")
    checkpoint = json.loads(paths["checkpoint"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    validate_checkpoint(checkpoint)
    assert manifest["checkpoint_payload_sha256"] == checkpoint["checkpoint_sha256"]
    assert manifest["pde_validated"] is False
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False
    assert checkpoint["truth_boundary"]["candidate_artifact_fabricated"] is False
