from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.kokuno_independent_radial_stress_admission import (
    ACTUAL_AUDITOR_AGENT4_HEAD,
    AUDIT_ARTIFACT_ID,
    AUDIT_ARTIFACT_ZIP_DIGEST,
    AUDIT_RECEIPT_CANONICAL_SHA256,
    AUDIT_RECEIPT_RAW_SHA256,
    AUDIT_WORKFLOW_RUN_ID,
    evaluate_independent_radial_stress_admission,
)

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "artifacts/constrained/kokuno_agent3/agent4_595_independent_radial_stress_audit.json"


def _raw() -> bytes:
    return AUDIT.read_bytes()


def test_exact_independent_radial_stress_receipt_is_identity_bound() -> None:
    result = evaluate_independent_radial_stress_admission(_raw())
    assert result["execution_identity_bound"] is True
    assert result["actual_auditor_agent4_head"] == ACTUAL_AUDITOR_AGENT4_HEAD
    assert result["audit_workflow_run_id"] == AUDIT_WORKFLOW_RUN_ID
    assert result["audit_artifact_id"] == AUDIT_ARTIFACT_ID
    assert result["audit_artifact_zip_digest"] == AUDIT_ARTIFACT_ZIP_DIGEST


def test_independent_fd6_cross_audit_recomputes_as_pass() -> None:
    result = evaluate_independent_radial_stress_admission(_raw())
    assert result["failed_guards"] == ()
    assert result["independent_compact_radial_stress_cross_audit_passed"] is True
    assert result["compact_radial_stress_operator_independently_admitted"] is True
    assert result["theta_e2"]["finest_relative_rms"] == pytest.approx(0.00010215080513304861)
    assert result["theta_e2"]["finest_relative_max"] == pytest.approx(0.00013808983925326713)
    assert result["theta_e2"]["refinement_ratios"] == pytest.approx(
        (12.029826646162858, 14.923789334921668)
    )
    assert result["axial_e1"]["finest_relative_rms"] == pytest.approx(0.0001864344008886948)
    assert result["axial_e1"]["finest_relative_max"] == pytest.approx(0.00023610373350348252)
    assert result["axial_e1"]["refinement_ratios"] == pytest.approx(
        (14.554891779445125, 14.786694125947507)
    )
    assert result["theta_e2"]["finest_sign_flip_mutation_relative_rms"] > 1.99
    assert result["axial_e1"]["finest_sign_flip_mutation_relative_rms"] > 1.99


def test_pass_promotes_only_operator_reuse_not_correction_cycle() -> None:
    result = evaluate_independent_radial_stress_admission(_raw())
    assert result["full_composite_radial_stress_execution_allowed_when_actual_defect_available"] is True
    assert result["real_oscillatory_component_defect_independently_consumed"] is True
    for key in (
        "full_same_cycle_composite_requested_stress_materialized",
        "agent1_leading_cross_terms_included",
        "matched_pressure_included",
        "restricted_forcing_included",
        "candidate_finite_head_mean_debt_materialized",
        "real_full_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "finite_correction_cycle_rerun_allowed",
        "finite_correction_cycle_run",
        "heldout_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
        "paper_exact",
        "blowup_proved",
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


def test_receipt_byte_or_threshold_substitution_fails_closed() -> None:
    raw = _raw()
    with pytest.raises(ValueError, match="receipt bytes"):
        evaluate_independent_radial_stress_admission(raw + b" ")

    payload = json.loads(raw.decode("utf-8"))
    payload["frozen_guards"]["finest_relative_rms_max"] = 1.0
    modified = (json.dumps(payload, indent=2, sort_keys=True) + "\n").encode("utf-8")
    with pytest.raises(ValueError, match="receipt bytes"):
        evaluate_independent_radial_stress_admission(modified)
