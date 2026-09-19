from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_time_derivative_routing_checkpoint import (
    AGENT3_MEAN_STRESS_HEAD,
    AGENT4_TIME_DERIVATIVE_ARTIFACT_ID,
    FORMAL_GATES,
    SCHEMA,
    TIME_DERIVATIVE_METRICS,
    _canonical_sha256,
    _time_derivative_passes,
    build_checkpoint,
    validate_checkpoint,
)


def test_v42_admits_only_independently_validated_time_derivative_handoff() -> None:
    payload = build_checkpoint()
    assert payload["schema"] == SCHEMA
    assert all(
        payload["upstream"]["agent4_independent_time_derivative_audit"][
            "derived_guard_passes"
        ].values()
    )
    states = payload["states"]
    assert states["oscillatory_ready"] is True
    assert states["public_oscillatory_time_derivative_handoff_ready"] is True
    assert states["public_oscillatory_time_derivative_independently_validated"] is True
    assert states["correction_ingest_allowed"] is True
    assert states["oscillatory_requested_stress_component_materialized"] is False
    assert states["same_cycle_requested_stress_materialized"] is False
    assert states["candidate_numeric_finite_head_mean_debt_materialized"] is False
    assert states["correction_ready"] is False
    assert states["velocity_export_ready"] is False
    assert states["pde_validated"] is False


def test_v42_records_a3_exact_head_generation_failure_without_promotion() -> None:
    payload = build_checkpoint()
    a3 = payload["upstream"]["agent3_latest_mean_stress_attempt"]
    assert a3["head"] == AGENT3_MEAN_STRESS_HEAD
    assert a3["dedicated_workflow_conclusion"] == "failure"
    assert a3["focused_tests_passed"] is True
    assert a3["receipt_generated"] is False
    assert a3["artifact_uploaded"] is False
    assert a3["scientific_promotion_allowed"] is False
    assert "receipt_generation_failed" in a3["blocker"]


def test_v42_fails_closed_if_independent_derivative_guard_is_mutated() -> None:
    metrics = copy.deepcopy(TIME_DERIVATIVE_METRICS)
    metrics["fd4_relative_rms"][-1] = 1.0e-3
    verdict = _time_derivative_passes(metrics)
    assert verdict["fd4_relative_rms"] is False
    with pytest.raises(ValueError, match="time-derivative handoff rejected"):
        build_checkpoint(metrics)


def test_v42_roundtrip_sha_and_truth_boundary_are_fail_closed() -> None:
    payload = build_checkpoint()
    encoded = json.dumps(payload, sort_keys=True)
    roundtrip = json.loads(encoded)
    validate_checkpoint(roundtrip)
    assert roundtrip["checkpoint_sha256"] == _canonical_sha256(roundtrip)
    assert roundtrip["formal_gates_unchanged"] == FORMAL_GATES
    assert (
        roundtrip["upstream"]["agent4_independent_time_derivative_audit"]["artifact_id"]
        == AGENT4_TIME_DERIVATIVE_ARTIFACT_ID
    )

    bad = copy.deepcopy(roundtrip)
    bad["states"]["pde_validated"] = True
    bad["checkpoint_sha256"] = _canonical_sha256(bad)
    with pytest.raises(ValueError, match="pde_validated"):
        validate_checkpoint(bad)


def test_v42_rejects_laundered_a3_artifact() -> None:
    payload = build_checkpoint()
    bad = copy.deepcopy(payload)
    bad["upstream"]["agent3_latest_mean_stress_attempt"]["artifact_uploaded"] = True
    bad["checkpoint_sha256"] = _canonical_sha256(bad)
    with pytest.raises(ValueError, match="laundered"):
        validate_checkpoint(bad)
