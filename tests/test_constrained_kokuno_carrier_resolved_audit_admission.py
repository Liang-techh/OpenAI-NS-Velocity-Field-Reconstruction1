from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_carrier_resolved_audit_admission import (
    ACTUAL_AUDITED_AGENT2_HEAD,
    ACTUAL_AUDITOR_AGENT4_HEAD,
    AUDIT_ARTIFACT_ID,
    AUDIT_ARTIFACT_ZIP_DIGEST,
    AUDIT_RECEIPT_CANONICAL_SHA256,
    AUDIT_RECEIPT_RAW_SHA256,
    AUDIT_WORKFLOW_RUN_ID,
    _canonical_payload_bytes,
    _recompute_scientific_verdict,
    _sha256_bytes,
    evaluate_carrier_resolved_audit_admission,
)

RECEIPT = Path(
    "artifacts/constrained/kokuno_agent3/agent4_553_carrier_resolved_independent_audit.json"
)


def test_exact_agent4_553_receipt_binds_execution_and_partial_scientific_passes() -> None:
    raw = RECEIPT.read_bytes()
    report = evaluate_carrier_resolved_audit_admission(raw)

    assert _sha256_bytes(raw) == AUDIT_RECEIPT_RAW_SHA256
    payload = json.loads(raw)
    assert _sha256_bytes(_canonical_payload_bytes(payload)) == AUDIT_RECEIPT_CANONICAL_SHA256

    assert report["execution_identity_bound"] is True
    assert report["actual_audited_agent2_head"] == ACTUAL_AUDITED_AGENT2_HEAD
    assert report["actual_auditor_agent4_head"] == ACTUAL_AUDITOR_AGENT4_HEAD
    assert report["audit_workflow_run_id"] == AUDIT_WORKFLOW_RUN_ID
    assert report["audit_artifact_id"] == AUDIT_ARTIFACT_ID
    assert report["audit_artifact_zip_digest"] == AUDIT_ARTIFACT_ZIP_DIGEST

    assert report["project_support_preflight_passed"] is True
    assert report["independent_candidate_covariance_rank_preflight_passed"] is True
    assert report["candidate_second_covariance_direction_independently_observed"] is True
    assert report["covariance_input_structurally_available"] is True
    assert report["minimum_covariance_rank_ratio"] == pytest.approx(0.09026569603651097)
    assert report["maximum_covariance_resolution_drift"] == pytest.approx(
        3.288617254673393e-13
    )


def test_exact_agent4_553_receipt_still_rejects_overall_oscillatory_handoff() -> None:
    report = evaluate_carrier_resolved_audit_admission(RECEIPT.read_bytes())

    assert report["finest_relative_divergence_rms"] == pytest.approx(5.451453500506105e-05)
    assert report["finest_relative_divergence_max"] == pytest.approx(6.930909589100082e-05)
    assert report["divergence_refinement_ratios"] == pytest.approx(
        (1.5199992545044203, 1.0233425481785872)
    )
    assert report["divergence_max_passed"] is True
    assert report["divergence_rms_passed"] is False
    assert report["divergence_refinement_passed"] is False
    assert report["failed_guards"] == (
        "finest_relative_divergence_rms",
        "minimum_divergence_refinement_ratio",
    )
    assert report["public_oscillatory_preflight_passed"] is False


def test_covariance_pass_does_not_open_any_correction_or_pde_claim() -> None:
    report = evaluate_carrier_resolved_audit_admission(RECEIPT.read_bytes())
    for key in (
        "correction_ingest_allowed",
        "same_cycle_requested_stress_materialization_allowed",
        "candidate_finite_head_mean_debt_materialization_allowed",
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
    assert report["future_overall_pass_requires_new_agent3_identity_pin"] is True


def test_recomputed_science_does_not_trust_declared_pass_bits() -> None:
    payload = json.loads(RECEIPT.read_text(encoding="utf-8"))
    payload["divergence"]["local_divergence_passed"] = True
    payload["public_oscillatory_preflight_passed"] = True
    scientific = _recompute_scientific_verdict(payload)
    assert scientific["local_divergence_passed"] is False
    assert scientific["public_oscillatory_preflight_passed"] is False
    assert scientific["covariance_rank_preflight_passed"] is True


def test_modified_or_reserialized_receipt_cannot_inherit_execution_identity() -> None:
    payload = json.loads(RECEIPT.read_text(encoding="utf-8"))
    payload["phase_mean_covariance"]["rows"][0]["rank_ratio"] = 0.5
    modified = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    with pytest.raises(ValueError, match="receipt bytes"):
        evaluate_carrier_resolved_audit_admission(modified)

    original = json.loads(RECEIPT.read_text(encoding="utf-8"))
    compact = json.dumps(original, sort_keys=True, separators=(",", ":")).encode("utf-8")
    assert json.loads(compact) == original
    with pytest.raises(ValueError, match="receipt bytes"):
        evaluate_carrier_resolved_audit_admission(compact)


def test_non_bytes_receipt_is_rejected() -> None:
    with pytest.raises(TypeError, match="must be bytes"):
        evaluate_carrier_resolved_audit_admission({})  # type: ignore[arg-type]
