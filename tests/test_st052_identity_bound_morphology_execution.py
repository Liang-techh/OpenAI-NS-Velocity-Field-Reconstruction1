from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction import st052_identity_bound_morphology_execution as mod


def _manifest() -> dict:
    return {
        "candidate_id": mod.CANDIDATE_ID,
        "whole_candidate_identity_sha256": "a" * 64,
        "runtime": {
            "exact_source_runtime_identity_sha256": mod.SOURCE_RUNTIME_IDENTITY_SHA256,
        },
        "truth_boundary": {
            "visualization_ready": False,
            "pde_validated": False,
            "source_correspondence_verified": False,
        },
    }


def _contract() -> dict:
    return {
        "candidate_id": mod.CANDIDATE_ID,
        "decision": "accept_authenticated_exact_source_checkout_as_external_visualization_runtime_dependency",
        "required_source_runtime": {
            "source_head": mod.SOURCE_HEAD,
            "source_tree": mod.SOURCE_TREE,
            "source_runtime_identity_sha256": mod.SOURCE_RUNTIME_IDENTITY_SHA256,
        },
        "runtime_policy": {
            "exact_source_checkout_dependency_explicitly_accepted": True,
            "resolved_dependency_versions_are_part_of_frozen_candidate_identity": False,
        },
        "truth_boundary": {
            "velocity_export_ready": False,
            "visualization_ready": False,
            "pde_validated": False,
        },
    }


def _receipt() -> dict:
    payload = mod.build_velocity_identity_payload(
        _manifest(), _contract(), dependency_contract_sha256="b" * 64
    )
    velocity_sha = mod._canonical_sha256(payload)
    candidate_sha = "a" * 64
    morphology = {
        "candidate_identity_bound": True,
        "cylindrical_morphology_diagnostic_ready": True,
        "candidate_identity": {
            "candidate_id": mod.CANDIDATE_ID,
            "candidate_sha256": candidate_sha,
            "velocity_identity_sha256": velocity_sha,
        },
        "protocol": {
            "source_numeric_targets_used": False,
            "renderer_or_camera_used": False,
            "pixel_loss_used": False,
        },
        "visualization_ready": False,
        "visual_correspondence_verified": False,
        "source_correspondence_verified": False,
        "pde_validated": False,
        "paper_exact": False,
        "openai_field_identified": False,
        "blowup_proved": False,
    }
    receipt = {
        "schema": mod.SCHEMA,
        "task_id": mod.TASK_ID,
        "provenance": {
            "main_morphology_base": mod.MAIN_MORPHOLOGY_BASE,
            "constrained_runtime_acceptance_merge": mod.CONSTRAINED_RUNTIME_ACCEPTANCE_MERGE,
            "source_head": mod.SOURCE_HEAD,
            "source_tree": mod.SOURCE_TREE,
        },
        "candidate_identity": {
            "candidate_id": mod.CANDIDATE_ID,
            "candidate_sha256": candidate_sha,
            "velocity_identity_payload": payload,
            "velocity_identity_sha256": velocity_sha,
        },
        "execution_runtime": {"resolved_versions_identity_bound": False},
        "morphology_receipt": morphology,
        "truth_boundary": dict(mod.TRUTH_BOUNDARY),
    }
    receipt["receipt_sha256"] = mod._canonical_sha256(receipt)
    return receipt


def test_velocity_identity_binds_runtime_acceptance_contract_bytes() -> None:
    first = mod.build_velocity_identity_payload(
        _manifest(), _contract(), dependency_contract_sha256="b" * 64
    )
    second = mod.build_velocity_identity_payload(
        _manifest(), _contract(), dependency_contract_sha256="c" * 64
    )
    assert mod._canonical_sha256(first) != mod._canonical_sha256(second)
    assert first["whole_candidate_identity_sha256"] == "a" * 64
    assert first["resolved_runtime_versions_identity_bound"] is False


def test_runtime_source_or_truth_drift_fails_closed() -> None:
    contract = _contract()
    contract["required_source_runtime"]["source_head"] = "0" * 40
    with pytest.raises(ValueError, match="source identity drift"):
        mod.build_velocity_identity_payload(
            _manifest(), contract, dependency_contract_sha256="b" * 64
        )

    contract = _contract()
    contract["truth_boundary"]["visualization_ready"] = True
    with pytest.raises(ValueError, match="truth boundary promoted"):
        mod.build_velocity_identity_payload(
            _manifest(), contract, dependency_contract_sha256="b" * 64
        )


def test_execution_receipt_binds_exact_candidate_and_rejects_promotions() -> None:
    receipt = _receipt()
    mod.validate_execution_receipt(receipt)

    bad = copy.deepcopy(receipt)
    bad.pop("receipt_sha256")
    bad["morphology_receipt"]["candidate_identity"]["candidate_sha256"] = "f" * 64
    bad["receipt_sha256"] = mod._canonical_sha256(bad)
    with pytest.raises(ValueError, match="different candidate identity"):
        mod.validate_execution_receipt(bad)

    bad = copy.deepcopy(receipt)
    bad.pop("receipt_sha256")
    bad["truth_boundary"]["pde_validated"] = True
    bad["receipt_sha256"] = mod._canonical_sha256(bad)
    with pytest.raises(ValueError, match="truth boundary drift"):
        mod.validate_execution_receipt(bad)

    bad = copy.deepcopy(receipt)
    bad.pop("receipt_sha256")
    bad["morphology_receipt"]["protocol"]["source_numeric_targets_used"] = True
    bad["receipt_sha256"] = mod._canonical_sha256(bad)
    with pytest.raises(ValueError, match="forbidden source/render target"):
        mod.validate_execution_receipt(bad)
