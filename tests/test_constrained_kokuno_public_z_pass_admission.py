from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_public_z_pass_admission import (
    ACTUAL_AUDITED_AGENT2_HEAD,
    ACTUAL_AUDITOR_AGENT4_HEAD,
    AUDIT_ARTIFACT_ID,
    AUDIT_ARTIFACT_ZIP_DIGEST,
    AUDIT_RECEIPT_CANONICAL_SHA256,
    AUDIT_RECEIPT_RAW_SHA256,
    AUDIT_WORKFLOW_RUN_ID,
    evaluate_public_z_pass_admission,
)

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "artifacts/constrained/kokuno_agent3/agent4_563_public_z_pullback_independent_audit.json"


def _raw() -> bytes:
    return AUDIT.read_bytes()


def test_exact_public_z_pass_receipt_is_identity_bound() -> None:
    result = evaluate_public_z_pass_admission(_raw())
    assert result["execution_identity_bound"] is True
    assert result["actual_audited_agent2_head"] == ACTUAL_AUDITED_AGENT2_HEAD
    assert result["actual_auditor_agent4_head"] == ACTUAL_AUDITOR_AGENT4_HEAD
    assert result["audit_workflow_run_id"] == AUDIT_WORKFLOW_RUN_ID
    assert result["audit_artifact_id"] == AUDIT_ARTIFACT_ID
    assert result["audit_artifact_zip_digest"] == AUDIT_ARTIFACT_ZIP_DIGEST
    assert result["correction_receipt_identity_pinned_to_current_pass"] is True


def test_exact_receipt_recomputes_all_frozen_local_guards_as_pass() -> None:
    result = evaluate_public_z_pass_admission(_raw())
    assert result["failed_guards"] == ()
    assert result["local_divergence_passed"] is True
    assert result["covariance_rank_preflight_passed"] is True
    assert result["project_support_preflight_passed"] is True
    assert result["public_oscillatory_preflight_passed"] is True
    assert result["finest_relative_divergence_rms"] == pytest.approx(1.767246633658725e-07)
    assert result["finest_relative_divergence_max"] == pytest.approx(2.0538436904463344e-07)
    assert result["divergence_refinement_ratios"] == pytest.approx(
        (15.824371829440508, 15.955526532124916)
    )
    assert result["minimum_covariance_rank_ratio"] == pytest.approx(0.09276278453864097)


def test_pass_opens_only_same_cycle_defect_materialization_stage() -> None:
    result = evaluate_public_z_pass_admission(_raw())
    assert result["correction_ingest_allowed"] is True
    assert result["same_cycle_requested_stress_materialization_allowed"] is True
    assert result["candidate_finite_head_mean_debt_materialization_allowed"] is True

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
        assert result[key] is False


def test_exact_raw_and_canonical_hashes_are_pinned() -> None:
    raw = _raw()
    assert hashlib.sha256(raw).hexdigest() == AUDIT_RECEIPT_RAW_SHA256
    payload = json.loads(raw.decode("utf-8"))
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    assert hashlib.sha256(canonical).hexdigest() == AUDIT_RECEIPT_CANONICAL_SHA256


def test_receipt_byte_or_value_substitution_fails_closed() -> None:
    raw = _raw()
    with pytest.raises(ValueError, match="receipt bytes"):
        evaluate_public_z_pass_admission(raw + b" ")

    payload = json.loads(raw.decode("utf-8"))
    payload["divergence"]["rows"][-1]["relative_rms"] = 0.0
    modified = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    with pytest.raises(ValueError, match="receipt bytes"):
        evaluate_public_z_pass_admission(modified)
