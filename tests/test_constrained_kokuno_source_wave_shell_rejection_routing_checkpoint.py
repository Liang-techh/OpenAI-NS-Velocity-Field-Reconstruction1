from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_source_wave_shell_rejection_routing_checkpoint import (
    AGENT4_VERIFIED_RECEIPT,
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


def test_checkpoint_retains_exact_agent4_rejection_receipt() -> None:
    checkpoint = build_checkpoint()
    assert validate_checkpoint(checkpoint) == checkpoint
    receipt = checkpoint["agent4_verified_receipt"]
    assert receipt == AGENT4_VERIFIED_RECEIPT
    assert receipt["source_pr"] == 371
    assert receipt["source_wave_shell_independent_cartesian_preflight_assessed"] is True
    assert receipt["source_wave_shell_independent_cartesian_preflight_passed"] is False
    assert receipt["failed_checks"] == ["divergence_refinement_ratio"]
    assert receipt["finest"]["step"] == pytest.approx(2.5e-4)
    assert receipt["finest"]["curl_relative_rms"] == pytest.approx(1.3245337639271412e-12)
    assert receipt["finest"]["divergence_max_abs"] == pytest.approx(2.1924837789035134e-10)
    assert receipt["curl_refinement_ratios"] == pytest.approx([15.994123110028152, 14.5964801840349])
    assert receipt["divergence_refinement_ratios"] == pytest.approx([11.789475096192975, 1.193958216011321])


def test_latest_a3_spacetime_envelope_is_bound_as_verified_sibling() -> None:
    checkpoint = build_checkpoint()
    agent3 = checkpoint["upstream"]["agent3"]
    assert agent3["pr"] == 370
    assert agent3["dedicated_run"] == 35331611617
    assert agent3["standard_run"] == 35331611700
    assert agent3["artifact_id"] == 10542155114
    assert agent3["radial_force_spacetime_envelope_ready"] is True
    assert agent3["all_cells_derivative_stable"] is True
    assert agent3["cell_count"] == 9
    assert agent3["radial_force_rms_min"] == pytest.approx(6.688247650141098e-05)
    assert agent3["radial_force_rms_max"] == pytest.approx(1.0448444107009646e-04)
    assert agent3["radial_to_tangential_force_rms_ratio_max"] == pytest.approx(0.688005978982281)
    assert agent3["fine_pair_relative_rms_difference_max"] == pytest.approx(0.001227751423269186)
    assert agent3["finite_correction_cycle_rerun_allowed"] is False


def test_pipeline_frontier_and_truth_states_remain_fail_closed() -> None:
    checkpoint = build_checkpoint()
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["states"] == STATES
    assert checkpoint["upstream"] == UPSTREAM
    assert checkpoint["pipeline_frontier"] == PIPELINE_FRONTIER
    assert checkpoint["states"]["outer_base_backbone_ready"] is True
    assert checkpoint["states"]["slow_squared_partition_family_contract_ready"] is True
    assert checkpoint["states"]["radial_force_spacetime_envelope_ready"] is True
    assert checkpoint["states"]["source_wave_shell_independent_cartesian_preflight_passed"] is False
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["states"]["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert checkpoint["states"]["genuinely_independent_second_covariance_column_ready"] is False
    assert checkpoint["states"]["correction_ready"] is False
    assert checkpoint["states"]["candidate_artifact_instantiated"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["formal_full_domain_pde_gate_assessed"] is False
    assert checkpoint["states"]["pde_validated"] is False
    assert all(stage["ready"] is False for stage in checkpoint["pipeline_frontier"])


def test_checkpoint_rejects_false_source_shell_promotion() -> None:
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    mutated["agent4_verified_receipt"]["source_wave_shell_independent_cartesian_preflight_passed"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="Agent-4 verified receipt changed|source-wave-shell rejection was falsely promoted"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_changed_failure_reason() -> None:
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    mutated["agent4_verified_receipt"]["failed_checks"] = []
    _resign(mutated)
    with pytest.raises(ValueError, match="Agent-4 verified receipt changed|rejection reason changed"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_threshold_change() -> None:
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    mutated["fixed_gates"]["held_out_normalized_full_momentum_max"] = 2e-3
    _resign(mutated)
    with pytest.raises(ValueError, match="fixed gates changed"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_a3_finite_cycle_promotion() -> None:
    checkpoint = build_checkpoint()
    mutated = copy.deepcopy(checkpoint)
    mutated["upstream"]["agent3"]["finite_correction_cycle_rerun_allowed"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="upstream verification receipt changed|finite correction cycle was falsely authorized"):
        validate_checkpoint(mutated)


def test_bundle_roundtrip_retains_artifact_provenance(tmp_path) -> None:
    paths = write_bundle(tmp_path / "bundle")
    checkpoint = json.loads(paths["checkpoint"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["manifest"].read_text(encoding="utf-8"))
    validate_checkpoint(checkpoint)
    assert manifest["checkpoint_payload_sha256"] == checkpoint["checkpoint_sha256"]
    assert manifest["bound_agent4_artifact_id"] == 10541526034
    assert manifest["bound_agent4_artifact_digest"] == "sha256:5f24d5b153b8d501b4c8f2cb312604392d1a922dd42723c15cf1644fcf1c8a92"
    assert manifest["bound_agent3_artifact_id"] == 10542155114
    assert manifest["bound_agent3_artifact_digest"] == "sha256:8cc8de6e34510183eaa18ba67d4b4e37ede022b69aa7016a7bbcb615f105b3a0"
    assert manifest["pde_validated"] is False
    assert checkpoint["truth_boundary"]["local_absolute_error_used_to_override_failed_refinement_guard"] is False
    assert checkpoint["truth_boundary"]["roundoff_sensitive_upstream_verdict_recomputed_by_agent5"] is False
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False
