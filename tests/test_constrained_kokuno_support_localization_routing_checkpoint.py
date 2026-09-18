from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_support_localized_cartesian_independent import (
    build_report as build_agent4_report,
)
from openai_ns_reconstruction.kokuno_support_localization_routing_checkpoint import (
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


def test_checkpoint_preserves_agent4_local_rejection(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    assert validate_checkpoint(checkpoint) == checkpoint
    audit = checkpoint["agent4_independent_support_localization_audit"]
    finest = audit["finest"]
    assert audit["structural_preflight_passed"] is False
    assert finest["step"] == pytest.approx(5.0e-4)
    assert finest["curl_relative_rms"] == pytest.approx(1.9695954432183738e-05)
    assert finest["divergence_max_abs"] == pytest.approx(3.777363724997359e-04)
    assert finest["regions"]["off_grid"]["divergence_max_abs"] == pytest.approx(
        2.0182680496128174e-05
    )
    assert finest["regions"]["small_radius"]["divergence_max_abs"] == pytest.approx(
        3.777363724997359e-04
    )
    checks = audit["local_guards"]["checks"]
    assert checks["finest_cartesian_curl_relative_rms"] is True
    assert checks["finest_cartesian_divergence_max_abs"] is False
    assert checks["curl_refinement_ratio"] is True
    assert checks["divergence_refinement_ratio"] is True
    assert checkpoint["states"]["support_localized_cartesian_preflight_passed"] is False
    assert checkpoint["states"]["pde_validated"] is False


def test_checkpoint_binds_lane_state_without_false_promotion(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    assert checkpoint["fixed_gates"] == FIXED_GATES
    assert checkpoint["upstream"]["agent2"]["consumed_in_executable_ancestry"] is True
    assert checkpoint["upstream"]["agent2"]["public_source_oscillatory_xyz_t_velocity_ready"] is False
    assert checkpoint["upstream"]["agent3"]["finite_correction_cycle_rerun_allowed"] is False
    assert checkpoint["upstream"]["agent3"]["duplicate_control_rank2_required_nodes"] == 0
    assert checkpoint["upstream"]["agent3"]["finest_nodes_requiring_second_direction"] == 99
    assert checkpoint["baseline_vs_kokuno"][1]["directly_comparable_to_ST006"] is False
    assert checkpoint["states"] == STATES
    assert checkpoint["states"]["leading_ready"] is False
    assert checkpoint["states"]["velocity_export_ready"] is False


def test_checkpoint_rejects_rehashed_truth_promotion(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    mutated = copy.deepcopy(checkpoint)
    mutated["truth_boundary"]["pde_validated"] = True
    unsigned = dict(mutated)
    unsigned.pop("checkpoint_sha256")
    mutated["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError, match="truth boundary promoted"):
        validate_checkpoint(mutated)


def test_checkpoint_rejects_rehashed_threshold_change(agent4_report: dict) -> None:
    checkpoint = build_checkpoint(agent4_report)
    mutated = copy.deepcopy(checkpoint)
    mutated["fixed_gates"]["held_out_normalized_full_momentum_max"] = 2.0e-3
    unsigned = dict(mutated)
    unsigned.pop("checkpoint_sha256")
    mutated["checkpoint_sha256"] = _sha(unsigned)
    with pytest.raises(ValueError, match="fixed gates changed"):
        validate_checkpoint(mutated)


def test_agent4_failed_guard_cannot_be_laundered_as_pass(agent4_report: dict) -> None:
    promoted = copy.deepcopy(agent4_report)
    promoted["local_guards"]["structural_preflight_passed"] = True
    with pytest.raises(ValueError, match="rejection disappeared"):
        _validate_agent4_report(promoted)


def test_bundle_round_trip_retains_scientific_failure(tmp_path) -> None:
    paths = write_bundle(tmp_path / "bundle")
    checkpoint = json.loads(paths["checkpoint"].read_text(encoding="utf-8"))
    validate_checkpoint(checkpoint)
    assert checkpoint["agent4_independent_support_localization_audit"][
        "structural_preflight_passed"
    ] is False
    assert checkpoint["truth_boundary"]["pde_validated"] is False
