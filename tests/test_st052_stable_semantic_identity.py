from __future__ import annotations

from copy import deepcopy
import hashlib

import pytest

from openai_ns_reconstruction import st052_stable_semantic_identity as stable


def _sha(label: str) -> str:
    return hashlib.sha256(label.encode("utf-8")).hexdigest()


def _candidate_manifest() -> dict:
    identity = {
        "schema": stable.LEGACY_WHOLE_SCHEMA,
        "parent_candidate_id": stable.PARENT_CANDIDATE_ID,
        "parent_source_head": stable.bridge.SOURCE_HEAD,
        "parent_replay_identity_sha256": _sha("parent-replay-semantics"),
        "parent_candidate_file_sha256": _sha("candidate-materialization-A"),
        "parent_validation_file_sha256": _sha("validation-materialization-A"),
        "parent_manifest_file_sha256": _sha("manifest-materialization-A"),
        "temporal_transform_spec_sha256": _sha("temporal-transform-semantics"),
        "static_transform_spec_sha256": _sha("static-temporal-semantics"),
        "exact_source_runtime_identity_sha256": stable.bridge.SOURCE_RUNTIME_IDENTITY_SHA256,
    }
    return {
        "schema": stable.LEGACY_WHOLE_SCHEMA,
        "task_id": "CR-A9-070",
        "candidate_id": stable.CANDIDATE_ID,
        "whole_candidate_identity_sha256": stable.canonical_sha256(identity),
        "identity_payload": identity,
        "files": {
            "parent_candidate.json": identity["parent_candidate_file_sha256"],
            "parent_validation.json": identity["parent_validation_file_sha256"],
            "parent_replay_manifest.json": identity["parent_manifest_file_sha256"],
            "temporal_transform_spec.json": _sha("temporal-spec-file-A"),
        },
        "runtime": {
            "exact_source_runtime_identity_sha256": stable.bridge.SOURCE_RUNTIME_IDENTITY_SHA256,
            "exact_source_runtime_identity_closed": True,
            "historical_module_cache_isolated": True,
        },
        "truth_boundary": {
            "velocity_export_ready": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "source_correspondence_verified": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def _dependency_contract() -> dict:
    return {
        "schema": stable.DEPENDENCY_SCHEMA,
        "task_id": stable.DEPENDENCY_TASK_ID,
        "candidate_id": stable.CANDIDATE_ID,
        "decision": stable.DEPENDENCY_DECISION,
        "rationale": "prose is evidence, not callable semantics",
        "required_source_runtime": {
            "source_head": stable.bridge.SOURCE_HEAD,
            "source_tree": stable.bridge.SOURCE_TREE,
            "recipe_path": stable.SOURCE_RECIPE_PATH,
            "recipe_git_blob_sha1": stable.SOURCE_RECIPE_GIT_BLOB_SHA1,
            "source_runtime_identity_sha256": stable.bridge.SOURCE_RUNTIME_IDENTITY_SHA256,
            "worktree_must_authenticate": True,
            "historical_module_cache_must_be_isolated": True,
        },
        "runtime_policy": {
            "vendoring_parent_into_installable_package_required_for_visualization_candidate": False,
            "exact_source_checkout_dependency_explicitly_accepted": True,
            "normal_package_plus_exact_source_checkout_load_path_accepted": True,
            "integration_smoke_must_record_resolved_python_numpy_scipy_sympy_versions": True,
            "resolved_dependency_versions_are_part_of_frozen_candidate_identity": False,
            "portable_bitwise_identity_across_dependency_builds_claimed": False,
            "standalone_package_parent_runtime_ready": False,
        },
        "truth_boundary": {
            "external_runtime_dependency_accepted_for_visualization_delivery": True,
            "standalone_package_parent_runtime_ready": False,
            "velocity_export_ready": False,
            "visualization_ready": False,
            "visual_correspondence_verified": False,
            "pde_validated": False,
            "paper_exact": False,
            "openai_field_identified": False,
            "blowup_proved": False,
        },
    }


def _receipt(manifest=None, contract=None, *, contract_file="contract-A") -> dict:
    return stable.build_receipt(
        _candidate_manifest() if manifest is None else manifest,
        _dependency_contract() if contract is None else contract,
        dependency_contract_file_sha256=_sha(contract_file),
    )


def test_materialization_hash_drift_does_not_change_semantic_candidate_or_velocity_identity():
    first_manifest = _candidate_manifest()
    first = _receipt(first_manifest)

    second_manifest = deepcopy(first_manifest)
    identity = second_manifest["identity_payload"]
    replacements = {
        "parent_candidate_file_sha256": _sha("candidate-materialization-B"),
        "parent_validation_file_sha256": _sha("validation-materialization-B"),
        "parent_manifest_file_sha256": _sha("manifest-materialization-B"),
    }
    identity.update(replacements)
    second_manifest["files"]["parent_candidate.json"] = replacements[
        "parent_candidate_file_sha256"
    ]
    second_manifest["files"]["parent_validation.json"] = replacements[
        "parent_validation_file_sha256"
    ]
    second_manifest["files"]["parent_replay_manifest.json"] = replacements[
        "parent_manifest_file_sha256"
    ]
    second_manifest["files"]["temporal_transform_spec.json"] = _sha(
        "temporal-spec-file-B"
    )
    second_manifest["whole_candidate_identity_sha256"] = stable.canonical_sha256(identity)
    second = _receipt(second_manifest, contract_file="contract-B")

    assert (
        first["stable_identity"]["candidate_semantic_identity_sha256"]
        == second["stable_identity"]["candidate_semantic_identity_sha256"]
    )
    assert (
        first["stable_identity"]["velocity_semantic_identity_sha256"]
        == second["stable_identity"]["velocity_semantic_identity_sha256"]
    )
    assert (
        first["materialization_evidence"]["legacy_whole_candidate_identity_sha256"]
        != second["materialization_evidence"]["legacy_whole_candidate_identity_sha256"]
    )
    assert first["materialization_evidence_sha256"] != second["materialization_evidence_sha256"]


def test_contract_prose_or_raw_file_sha_drift_does_not_change_normalized_runtime_semantics():
    contract_a = _dependency_contract()
    contract_b = deepcopy(contract_a)
    contract_b["rationale"] = "different prose, identical accepted runtime semantics"
    first = _receipt(contract=contract_a, contract_file="contract-file-A")
    second = _receipt(contract=contract_b, contract_file="contract-file-B")
    assert (
        first["stable_identity"]["runtime_acceptance_semantics_sha256"]
        == second["stable_identity"]["runtime_acceptance_semantics_sha256"]
    )
    assert (
        first["stable_identity"]["velocity_semantic_identity_sha256"]
        == second["stable_identity"]["velocity_semantic_identity_sha256"]
    )
    assert first["materialization_evidence_sha256"] != second["materialization_evidence_sha256"]


def test_true_parent_replay_semantic_drift_creates_new_candidate_and_velocity_identity():
    changed = _candidate_manifest()
    changed["identity_payload"]["parent_replay_identity_sha256"] = _sha(
        "different-parent-replay-semantics"
    )
    changed["whole_candidate_identity_sha256"] = stable.canonical_sha256(
        changed["identity_payload"]
    )
    first = _receipt()
    second = _receipt(changed)
    assert (
        first["stable_identity"]["candidate_semantic_identity_sha256"]
        != second["stable_identity"]["candidate_semantic_identity_sha256"]
    )
    assert (
        first["stable_identity"]["velocity_semantic_identity_sha256"]
        != second["stable_identity"]["velocity_semantic_identity_sha256"]
    )


def test_temporal_semantic_drift_creates_new_identity_instead_of_reusing_old_hash():
    changed = _candidate_manifest()
    changed["identity_payload"]["temporal_transform_spec_sha256"] = _sha(
        "different-temporal-transform"
    )
    changed["whole_candidate_identity_sha256"] = stable.canonical_sha256(
        changed["identity_payload"]
    )
    assert (
        _receipt()["stable_identity"]["candidate_semantic_identity_sha256"]
        != _receipt(changed)["stable_identity"]["candidate_semantic_identity_sha256"]
    )


def test_exact_source_runtime_drift_fails_closed():
    changed = _candidate_manifest()
    changed["identity_payload"]["exact_source_runtime_identity_sha256"] = _sha(
        "foreign-runtime"
    )
    changed["runtime"]["exact_source_runtime_identity_sha256"] = changed[
        "identity_payload"
    ]["exact_source_runtime_identity_sha256"]
    changed["whole_candidate_identity_sha256"] = stable.canonical_sha256(
        changed["identity_payload"]
    )
    with pytest.raises(ValueError, match="runtime semantic identity drift"):
        _receipt(changed)


def test_runtime_acceptance_semantic_drift_fails_closed():
    changed = _dependency_contract()
    changed["runtime_policy"][
        "resolved_dependency_versions_are_part_of_frozen_candidate_identity"
    ] = True
    with pytest.raises(ValueError, match="runtime-dependency policy drift"):
        _receipt(contract=changed)


def test_truth_promotion_fails_closed():
    changed = _candidate_manifest()
    changed["truth_boundary"]["visualization_ready"] = True
    with pytest.raises(ValueError, match="visualization_ready"):
        _receipt(changed)


def test_receipt_rejects_laundering_materialization_evidence_into_stable_identity():
    receipt = _receipt()
    receipt["materialization_evidence"]["included_in_stable_candidate_identity"] = True
    receipt["materialization_evidence_sha256"] = stable.canonical_sha256(
        receipt["materialization_evidence"]
    )
    receipt["receipt_sha256"] = None
    with pytest.raises(ValueError, match="leaked into stable candidate identity"):
        stable.validate_receipt(receipt)
