from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_oscillatory_audit_identity_admission import (
    ACTUAL_AUDITED_AGENT2_HEAD,
    ACTUAL_AUDITOR_AGENT4_HEAD,
    AUDIT_ARTIFACT_ID,
    AUDIT_ARTIFACT_ZIP_DIGEST,
    AUDIT_RECEIPT_CANONICAL_SHA256,
    AUDIT_RECEIPT_RAW_SHA256,
    AUDIT_WORKFLOW_RUN_ID,
    PROTOCOL_ORIGIN_AGENT2_HEAD,
    _canonical_payload_bytes,
    _sha256_bytes,
    evaluate_pinned_oscillatory_audit_admission,
)

RECEIPT = Path("artifacts/constrained/kokuno_agent3/agent4_543_public_oscillatory_audit.json")


def test_exact_agent4_543_receipt_is_identity_bound_but_scientifically_rejected() -> None:
    raw = RECEIPT.read_bytes()
    report = evaluate_pinned_oscillatory_audit_admission(raw)

    assert _sha256_bytes(raw) == AUDIT_RECEIPT_RAW_SHA256
    payload = json.loads(raw)
    assert _sha256_bytes(_canonical_payload_bytes(payload)) == AUDIT_RECEIPT_CANONICAL_SHA256

    assert report["execution_identity_bound"] is True
    assert report["protocol_origin_agent2_head"] == PROTOCOL_ORIGIN_AGENT2_HEAD
    assert report["protocol_origin_semantics"] == "validator_protocol_origin_metadata_only"
    assert report["actual_audited_agent2_head"] == ACTUAL_AUDITED_AGENT2_HEAD
    assert report["actual_audited_agent2_head"] != report["protocol_origin_agent2_head"]
    assert report["actual_auditor_agent4_head"] == ACTUAL_AUDITOR_AGENT4_HEAD
    assert report["audit_workflow_run_id"] == AUDIT_WORKFLOW_RUN_ID
    assert report["audit_artifact_id"] == AUDIT_ARTIFACT_ID
    assert report["audit_artifact_zip_digest"] == AUDIT_ARTIFACT_ZIP_DIGEST

    assert report["scientific_failed_guards"] == (
        "finest_relative_divergence_rms",
        "finest_relative_divergence_max",
        "minimum_divergence_refinement_ratio",
        "minimum_covariance_rank_ratio",
    )
    assert report["independent_public_oscillatory_preflight_passed"] is False
    assert report["correction_ingest_allowed"] is False
    assert report["same_cycle_requested_stress_materialization_allowed"] is False
    assert report["candidate_finite_head_mean_debt_materialization_allowed"] is False
    assert report["finite_correction_cycle_rerun_allowed"] is False
    assert report["future_pass_requires_new_agent3_identity_pin"] is True


def test_modified_scientific_receipt_cannot_inherit_agent4_543_execution_identity() -> None:
    payload = json.loads(RECEIPT.read_text(encoding="utf-8"))
    payload["divergence"]["rows"][-1]["relative_rms"] = 1.0e-6
    payload["divergence"]["rows"][-1]["relative_max"] = 1.0e-6
    payload["divergence"]["relative_rms_refinement_ratios"] = [4.0, 4.0]
    payload["divergence"]["local_divergence_passed"] = True
    for row in payload["phase_mean_covariance"]["rows"]:
        row["rank_ratio"] = 0.03
    payload["phase_mean_covariance"]["covariance_rank_preflight_passed"] = True
    payload["local_curl_covariance_preflight_passed"] = True
    payload["public_oscillatory_preflight_passed"] = True
    modified = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")

    with pytest.raises(ValueError, match="receipt bytes"):
        evaluate_pinned_oscillatory_audit_admission(modified)


def test_semantically_equivalent_reserialization_does_not_count_as_exact_artifact_member() -> None:
    payload = json.loads(RECEIPT.read_text(encoding="utf-8"))
    compact = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    assert json.loads(compact) == payload

    with pytest.raises(ValueError, match="receipt bytes"):
        evaluate_pinned_oscillatory_audit_admission(compact)


def test_pinned_execution_keeps_all_downstream_scientific_claims_closed() -> None:
    report = evaluate_pinned_oscillatory_audit_admission(RECEIPT.read_bytes())
    for key in (
        "requested_stress_actual_state_values_materialized",
        "finite_head_mean_debt_materialized",
        "real_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "finite_correction_cycle_rerun_allowed",
        "finite_correction_cycle_run",
        "heldout_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
        "surrogate_defect_used",
    ):
        assert report[key] is False


def test_non_bytes_receipt_is_rejected() -> None:
    with pytest.raises(TypeError, match="must be bytes"):
        evaluate_pinned_oscillatory_audit_admission({})  # type: ignore[arg-type]
