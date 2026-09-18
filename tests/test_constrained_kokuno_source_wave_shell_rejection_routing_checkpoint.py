from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_source_wave_shell_cartesian_independent import (
    build_report as build_agent4_report,
)
from openai_ns_reconstruction.kokuno_source_wave_shell_rejection_routing_checkpoint import (
    FIXED_GATES,
    STATES,
    _sha,
    _validate_agent4_report,
    build_checkpoint,
    validate_checkpoint,
    write_bundle,
)


@pytest.fixture(scope="module")
def agent4_report() -> dict:
    return build_agent4_report()


def _resign(payload: dict) -> None:
    unsigned = dict(payload)
    unsigned.pop("checkpoint_sha256", None)
    payload["checkpoint_sha256"] = _sha(unsigned)


def test_checkpoint_retains_agent4_source_shell_rejection(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    assert validate_checkpoint(checkpoint) == checkpoint
    audit = checkpoint["agent4_source_wave_shell_audit"]
    assert audit["source_wave_shell_independent_cartesian_preflight_assessed"] is True
    assert audit["source_wave_shell_independent_cartesian_preflight_passed"] is False
    assert audit["failed_checks"] == ["divergence_refinement_ratio"]
    assert audit["finest"]["step"] == pytest.approx(2.5e-4)
    assert audit["finest"]["curl_relative_rms"] == pytest.approx(1.3245337639271412e-12, rel=1e-3)
    assert audit["finest"]["divergence_max_abs"] == pytest.approx(2.1924837789035134e-10, rel=1e-3)
    assert audit["curl_refinement_ratios"] == pytest.approx([15.994123110028152, 14.5964801840349], rel=1e-3)
    assert audit["divergence_refinement_ratios"] == pytest.approx([11.789475096192975, 1.193958216011321], rel=1e-3)


def test_rejection_is_only_refinement_guard(agent4_report: dict) -> None:
    _validate_agent4_report(agent4_report)
    checks = agent4_report["local_guards"]["checks"]
    assert checks["divergence_refinement_ratio"] is False
    assert all(value is True for key, value in checks.items() if key != "divergence_refinement_ratio")
    assert agent4_report["local_guards"]["thresholds"]["minimum_refinement_ratio"] == pytest.approx(2.5)
    assert agent4_report["local_guards"]["thresholds"]["finest_divergence_max_abs_max"] == pytest.approx(2e-4)


def test_pipeline_frontier_and_sibling_truths_remain_blocked(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["states"] == STATES
    assert checkpoint["upstream"]["agent1"]["pr"] == 368
    assert checkpoint["upstream"]["agent1"]["outer_base_schedule_ready"] is True
    assert checkpoint["upstream"]["agent1"]["leading_ready"] is False
    assert checkpoint["upstream"]["agent2"]["pr"] == 369
    assert checkpoint["upstream"]["agent2"]["consumed_in_executable_ancestry"] is True
    assert checkpoint["upstream"]["agent2"]["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert checkpoint["upstream"]["agent3"]["finite_correction_cycle_rerun_allowed"] is False
    assert all(stage["ready"] is False for stage in checkpoint["pipeline_frontier"])
    assert checkpoint["states"]["velocity_export_ready"] is False
    assert checkpoint["states"]["pde_validated"] is False


def test_checkpoint_rejects_false_local_promotion(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    mutated = copy.deepcopy(checkpoint)
    mutated["agent4_source_wave_shell_audit"]["source_wave_shell_independent_cartesian_preflight_passed"] = True
    _resign(mutated)
    with pytest.raises(ValueError, match="source-wave-shell rejection was falsely promoted"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_changed_failure_reason(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    mutated = copy.deepcopy(checkpoint)
    mutated["agent4_source_wave_shell_audit"]["failed_checks"] = []
    _resign(mutated)
    with pytest.raises(ValueError, match="rejection reason changed"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_threshold_change(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    mutated = copy.deepcopy(checkpoint)
    mutated["fixed_gates"]["held_out_normalized_full_momentum_max"] = 2e-3
    _resign(mutated)
    with pytest.raises(ValueError, match="fixed gates changed"):
        validate_checkpoint(mutated)


def test_agent4_rejection_must_not_be_erased(agent4_report: dict) -> None:
    mutated = copy.deepcopy(agent4_report)
    mutated["local_guards"]["source_wave_shell_independent_cartesian_preflight_passed"] = True
    with pytest.raises(ValueError, match="expected source-wave-shell local rejection is absent"):
        _validate_agent4_report(mutated)


def test_bundle_roundtrip_retains_truth_boundary(tmp_path, agent4_report: dict) -> None:
    paths = write_bundle(tmp_path / "bundle")
    checkpoint = json.loads(paths["checkpoint"].read_text(encoding="utf-8"))
    report = json.loads(paths["report"].read_text(encoding="utf-8"))
    validate_checkpoint(checkpoint)
    _validate_agent4_report(report)
    assert checkpoint["truth_boundary"]["local_absolute_error_used_to_override_failed_refinement_guard"] is False
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False
    assert checkpoint["states"]["candidate_artifact_instantiated"] is False
