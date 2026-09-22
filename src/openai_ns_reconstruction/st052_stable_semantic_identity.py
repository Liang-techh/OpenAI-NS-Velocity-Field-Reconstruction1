"""Stable semantic identity for the frozen ST052-M temporal child.

The legacy whole-child capsule intentionally authenticates every materialized
file, but its composite identity also includes regenerated candidate/validation/
manifest file SHA-256 values.  The exact ST052 replay contract explicitly says
those raw bytes can change while the frozen mathematical replay identity stays
fixed.  This module adds a *new* identity schema that separates those roles:

* semantic candidate identity: parent replay semantics + temporal-transform
  semantics + exact source-runtime semantics;
* semantic callable identity: semantic candidate identity + normalized accepted
  runtime-dependency semantics + the exact constrained callable implementation;
* materialization evidence: legacy whole-child/file/contract checksums, still
  authenticated and retained, but not folded into either semantic digest.

This is an identity/governance adapter only.  It does not modify the ST052
velocity, pressure, forcing, morphology, residual, or scientific thresholds.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import re
import subprocess
from typing import Any

from . import st052_identity_bound_morphology_execution as bridge

SCHEMA = "st052-stable-semantic-identity-split/v1"
CANDIDATE_SCHEMA = "st052-linear-temporal-semantic-candidate/v1"
VELOCITY_SCHEMA = "st052-stable-callable-velocity-identity/v1"
TASK_ID = "CR-A9-105"
CANDIDATE_ID = bridge.CANDIDATE_ID
PARENT_CANDIDATE_ID = "ST052-M"
LEGACY_WHOLE_SCHEMA = "st052-linear-temporal-whole-candidate-capsule/v3"
DEPENDENCY_SCHEMA = "st052m-runtime-dependency-acceptance/v1"
DEPENDENCY_TASK_ID = "A8-DELIVERY-01"
DEPENDENCY_DECISION = (
    "accept_authenticated_exact_source_checkout_as_external_visualization_runtime_dependency"
)
SOURCE_RECIPE_PATH = "experiments/root_st052/recipe.json"
SOURCE_RECIPE_GIT_BLOB_SHA1 = "e30c769052379f72afeee46ca264482884cc5ac7"
CONSTRAINED_RUNTIME_IMPLEMENTATION_REPOSITORY = (
    "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1"
)
CONSTRAINED_RUNTIME_IMPLEMENTATION_COMMIT = bridge.CONSTRAINED_RUNTIME_ACCEPTANCE_MERGE
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_SHA1 = re.compile(r"^[0-9a-f]{40}$")

TRUTH_BOUNDARY = {
    "candidate_changed": False,
    "velocity_coefficients_changed": False,
    "pressure_changed": False,
    "forcing_changed": False,
    "scientific_threshold_changed": False,
    "legacy_materialization_identity_relabelled_stable": False,
    "materialization_file_hashes_identity_bound": False,
    "constrained_callable_implementation_identity_bound": True,
    "stable_semantic_identity_ready": True,
    "velocity_export_ready": False,
    "visualization_ready": False,
    "visual_correspondence_verified": False,
    "source_correspondence_verified": False,
    "pde_validated": False,
    "paper_exact": False,
    "openai_field_identified": False,
    "blowup_proved": False,
}


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    return hashlib.sha256(_canonical_bytes(value)).hexdigest()


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_dict(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase 64-hex SHA-256")
    return value


def _require_sha1(value: Any, label: str) -> str:
    if not isinstance(value, str) or _SHA1.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase 40-hex SHA-1")
    return value


def _require_false(mapping: dict[str, Any], key: str, label: str) -> None:
    if mapping.get(key) is not False:
        raise ValueError(f"{label}.{key} must be exactly false")


def _authenticate_constrained_runtime_checkout(constrained_root: Path) -> str:
    """Bind callable execution to the exact clean constrained implementation commit."""
    try:
        head = subprocess.run(
            ["git", "-C", str(constrained_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "-C", str(constrained_root), "status", "--porcelain"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError) as exc:
        raise ValueError("unable to authenticate constrained callable runtime checkout") from exc
    _require_sha1(head, "constrained callable runtime HEAD")
    if head != CONSTRAINED_RUNTIME_IMPLEMENTATION_COMMIT:
        raise ValueError("constrained callable runtime implementation commit drift")
    if status:
        raise ValueError("constrained callable runtime worktree must be clean")
    return head


def candidate_semantic_payload(candidate_manifest: dict[str, Any]) -> dict[str, Any]:
    """Extract only output-defining semantics from a verified legacy bundle manifest."""
    if candidate_manifest.get("schema") != LEGACY_WHOLE_SCHEMA:
        raise ValueError("unexpected legacy ST052 whole-child schema")
    if candidate_manifest.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("unexpected ST052 whole-child candidate id")

    identity = _require_dict(candidate_manifest.get("identity_payload"), "identity_payload")
    runtime = _require_dict(candidate_manifest.get("runtime"), "runtime")
    truth = _require_dict(candidate_manifest.get("truth_boundary"), "truth_boundary")

    # The old digest must still be internally authentic: it is retained as
    # materialization evidence, not silently discarded.
    legacy_id = _require_sha256(
        candidate_manifest.get("whole_candidate_identity_sha256"),
        "whole_candidate_identity_sha256",
    )
    if canonical_sha256(identity) != legacy_id:
        raise ValueError("legacy whole-candidate identity checksum mismatch")

    if identity.get("parent_candidate_id") != PARENT_CANDIDATE_ID:
        raise ValueError("parent candidate id drift")
    if identity.get("parent_source_head") != bridge.SOURCE_HEAD:
        raise ValueError("parent source head drift")
    replay_id = _require_sha256(
        identity.get("parent_replay_identity_sha256"),
        "parent_replay_identity_sha256",
    )
    temporal_id = _require_sha256(
        identity.get("temporal_transform_spec_sha256"),
        "temporal_transform_spec_sha256",
    )
    static_temporal_id = _require_sha256(
        identity.get("static_transform_spec_sha256"),
        "static_transform_spec_sha256",
    )
    runtime_id = _require_sha256(
        identity.get("exact_source_runtime_identity_sha256"),
        "exact_source_runtime_identity_sha256",
    )
    if runtime_id != bridge.SOURCE_RUNTIME_IDENTITY_SHA256:
        raise ValueError("exact source runtime semantic identity drift")
    if runtime.get("exact_source_runtime_identity_sha256") != runtime_id:
        raise ValueError("manifest runtime identity disagrees with identity payload")
    if runtime.get("exact_source_runtime_identity_closed") is not True:
        raise ValueError("exact source runtime identity is not closed")
    if runtime.get("historical_module_cache_isolated") is not True:
        raise ValueError("historical module-cache isolation is not asserted")

    for key in (
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "source_correspondence_verified",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_false(truth, key, "legacy truth_boundary")

    return {
        "schema": CANDIDATE_SCHEMA,
        "candidate_id": CANDIDATE_ID,
        "parent": {
            "candidate_id": PARENT_CANDIDATE_ID,
            "source_head": bridge.SOURCE_HEAD,
            "replay_identity_sha256": replay_id,
        },
        "temporal_transform": {
            "spec_sha256": temporal_id,
            "static_spec_sha256": static_temporal_id,
        },
        "exact_source_runtime": {
            "source_head": bridge.SOURCE_HEAD,
            "source_tree": bridge.SOURCE_TREE,
            "identity_sha256": runtime_id,
        },
    }


def normalized_runtime_acceptance(dependency_contract: dict[str, Any]) -> dict[str, Any]:
    """Normalize only the accepted runtime semantics, excluding prose/file bytes."""
    if dependency_contract.get("schema") != DEPENDENCY_SCHEMA:
        raise ValueError("runtime-dependency contract schema drift")
    if dependency_contract.get("task_id") != DEPENDENCY_TASK_ID:
        raise ValueError("runtime-dependency contract task drift")
    if dependency_contract.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("runtime-dependency candidate drift")
    if dependency_contract.get("decision") != DEPENDENCY_DECISION:
        raise ValueError("runtime-dependency acceptance decision drift")

    required = _require_dict(
        dependency_contract.get("required_source_runtime"), "required_source_runtime"
    )
    policy = _require_dict(dependency_contract.get("runtime_policy"), "runtime_policy")
    truth = _require_dict(dependency_contract.get("truth_boundary"), "dependency truth_boundary")

    expected_required = {
        "source_head": bridge.SOURCE_HEAD,
        "source_tree": bridge.SOURCE_TREE,
        "recipe_path": SOURCE_RECIPE_PATH,
        "recipe_git_blob_sha1": SOURCE_RECIPE_GIT_BLOB_SHA1,
        "source_runtime_identity_sha256": bridge.SOURCE_RUNTIME_IDENTITY_SHA256,
        "worktree_must_authenticate": True,
        "historical_module_cache_must_be_isolated": True,
    }
    for key, expected in expected_required.items():
        if required.get(key) != expected:
            raise ValueError(f"runtime-dependency required-source drift: {key}")
    _require_sha1(required.get("recipe_git_blob_sha1"), "recipe_git_blob_sha1")
    _require_sha256(
        required.get("source_runtime_identity_sha256"),
        "source_runtime_identity_sha256",
    )

    expected_policy = {
        "vendoring_parent_into_installable_package_required_for_visualization_candidate": False,
        "exact_source_checkout_dependency_explicitly_accepted": True,
        "normal_package_plus_exact_source_checkout_load_path_accepted": True,
        "integration_smoke_must_record_resolved_python_numpy_scipy_sympy_versions": True,
        "resolved_dependency_versions_are_part_of_frozen_candidate_identity": False,
        "portable_bitwise_identity_across_dependency_builds_claimed": False,
        "standalone_package_parent_runtime_ready": False,
    }
    for key, expected in expected_policy.items():
        if policy.get(key) is not expected:
            raise ValueError(f"runtime-dependency policy drift: {key}")

    if truth.get("external_runtime_dependency_accepted_for_visualization_delivery") is not True:
        raise ValueError("external runtime dependency acceptance truth is absent")
    for key in (
        "standalone_package_parent_runtime_ready",
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        _require_false(truth, key, "dependency truth_boundary")

    return {
        "schema": "st052-normalized-runtime-acceptance-semantics/v1",
        "candidate_id": CANDIDATE_ID,
        "decision": DEPENDENCY_DECISION,
        "required_source_runtime": expected_required,
        "runtime_policy": expected_policy,
    }


def velocity_semantic_payload(
    candidate_semantics: dict[str, Any], dependency_contract: dict[str, Any]
) -> dict[str, Any]:
    candidate_id = canonical_sha256(candidate_semantics)
    acceptance = normalized_runtime_acceptance(dependency_contract)
    return {
        "schema": VELOCITY_SCHEMA,
        "candidate_id": CANDIDATE_ID,
        "candidate_semantic_identity_sha256": candidate_id,
        "runtime_acceptance_semantics_sha256": canonical_sha256(acceptance),
        "exact_source_runtime": {
            "source_head": bridge.SOURCE_HEAD,
            "source_tree": bridge.SOURCE_TREE,
            "source_runtime_identity_sha256": bridge.SOURCE_RUNTIME_IDENTITY_SHA256,
        },
        "constrained_callable_runtime": {
            "repository": CONSTRAINED_RUNTIME_IMPLEMENTATION_REPOSITORY,
            "commit": CONSTRAINED_RUNTIME_IMPLEMENTATION_COMMIT,
        },
        "callable": "St052LinearTemporalWholeCandidate.velocity(x,y,z,t)->[...,3]",
        "resolved_runtime_versions_identity_bound": False,
    }


def materialization_evidence(
    candidate_manifest: dict[str, Any], *, dependency_contract_file_sha256: str
) -> dict[str, Any]:
    identity = _require_dict(candidate_manifest.get("identity_payload"), "identity_payload")
    files = _require_dict(candidate_manifest.get("files"), "files")
    legacy_id = _require_sha256(
        candidate_manifest.get("whole_candidate_identity_sha256"),
        "whole_candidate_identity_sha256",
    )
    if canonical_sha256(identity) != legacy_id:
        raise ValueError("legacy materialization identity is not self-consistent")

    volatile = {}
    for key in (
        "parent_candidate_file_sha256",
        "parent_validation_file_sha256",
        "parent_manifest_file_sha256",
    ):
        volatile[key] = _require_sha256(identity.get(key), key)
    file_hashes: dict[str, str] = {}
    for name, digest in sorted(files.items()):
        if not isinstance(name, str):
            raise ValueError("materialization filename must be a string")
        file_hashes[name] = _require_sha256(digest, f"files[{name!r}]")

    return {
        "schema": "st052-materialization-evidence/v1",
        "legacy_whole_candidate_identity_sha256": legacy_id,
        "regenerated_parent_file_sha256": volatile,
        "bundle_file_sha256": file_hashes,
        "runtime_dependency_contract_file_sha256": _require_sha256(
            dependency_contract_file_sha256,
            "dependency_contract_file_sha256",
        ),
        "included_in_stable_candidate_identity": False,
        "included_in_stable_velocity_identity": False,
        "purpose": "byte-level authentication and reproducibility evidence only",
    }


def build_receipt(
    candidate_manifest: dict[str, Any],
    dependency_contract: dict[str, Any],
    *,
    dependency_contract_file_sha256: str,
) -> dict[str, Any]:
    candidate_payload = candidate_semantic_payload(candidate_manifest)
    candidate_id = canonical_sha256(candidate_payload)
    acceptance = normalized_runtime_acceptance(dependency_contract)
    velocity_payload = velocity_semantic_payload(candidate_payload, dependency_contract)
    velocity_id = canonical_sha256(velocity_payload)
    evidence = materialization_evidence(
        candidate_manifest,
        dependency_contract_file_sha256=dependency_contract_file_sha256,
    )

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "candidate_id": CANDIDATE_ID,
        "stable_identity": {
            "candidate_semantic_payload": candidate_payload,
            "candidate_semantic_identity_sha256": candidate_id,
            "velocity_semantic_payload": velocity_payload,
            "velocity_semantic_identity_sha256": velocity_id,
            "runtime_acceptance_semantics": acceptance,
            "runtime_acceptance_semantics_sha256": canonical_sha256(acceptance),
        },
        "materialization_evidence": evidence,
        "materialization_evidence_sha256": canonical_sha256(evidence),
        "migration_boundary": {
            "classification": "internal_semantic_identity_schema_upgrade",
            "source_repo": "Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1",
            "source_commit": bridge.CONSTRAINED_RUNTIME_ACCEPTANCE_MERGE,
            "license": "repository-internal; no third-party implementation migrated",
            "migration_scope": (
                "identity semantics only: split stable frozen-field/callable identity from regenerated byte checksums"
            ),
            "difference_from_legacy": (
                "legacy whole_candidate_identity_sha256 binds regenerated parent candidate/validation/manifest file SHA-256 values; "
                "the new candidate semantic identity excludes those raw file hashes while the callable identity additionally binds the exact constrained implementation commit"
            ),
            "legacy_identity_reused_as_stable_identity": False,
        },
        "direct_contribution": (
            "makes repeated exact-source rematerializations of the same frozen ST052 [u,v,w] addressable by one semantic candidate identity while binding callable identity to the exact implementation commit and retaining byte-level evidence separately"
        ),
        "remaining_limits": [
            "this schema does not itself rerun or admit the failed #1095 axial-aspect receipt",
            "ST052 still lacks the separately governed repository-unified velocity API promotion",
            "no public-source visual correspondence verdict is made",
            "no compatible pressure/restricted forcing or fresh complete-NS validation is performed",
        ],
        "truth_boundary": dict(TRUTH_BOUNDARY),
    }
    receipt["receipt_sha256"] = canonical_sha256(receipt)
    validate_receipt(receipt)
    return receipt


def validate_receipt(receipt: dict[str, Any]) -> None:
    if receipt.get("schema") != SCHEMA or receipt.get("task_id") != TASK_ID:
        raise ValueError("unexpected stable-identity receipt schema/task")
    if receipt.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("stable-identity candidate id drift")
    stable = _require_dict(receipt.get("stable_identity"), "stable_identity")
    evidence = _require_dict(receipt.get("materialization_evidence"), "materialization_evidence")
    migration = _require_dict(receipt.get("migration_boundary"), "migration_boundary")
    truth = _require_dict(receipt.get("truth_boundary"), "truth_boundary")

    candidate_payload = _require_dict(
        stable.get("candidate_semantic_payload"), "candidate_semantic_payload"
    )
    candidate_id = _require_sha256(
        stable.get("candidate_semantic_identity_sha256"),
        "candidate_semantic_identity_sha256",
    )
    if canonical_sha256(candidate_payload) != candidate_id:
        raise ValueError("stable candidate semantic digest mismatch")

    velocity_payload = _require_dict(
        stable.get("velocity_semantic_payload"), "velocity_semantic_payload"
    )
    velocity_id = _require_sha256(
        stable.get("velocity_semantic_identity_sha256"),
        "velocity_semantic_identity_sha256",
    )
    if canonical_sha256(velocity_payload) != velocity_id:
        raise ValueError("stable velocity semantic digest mismatch")
    if velocity_payload.get("candidate_semantic_identity_sha256") != candidate_id:
        raise ValueError("stable velocity/candidate identity linkage mismatch")

    acceptance = _require_dict(
        stable.get("runtime_acceptance_semantics"), "runtime_acceptance_semantics"
    )
    acceptance_id = _require_sha256(
        stable.get("runtime_acceptance_semantics_sha256"),
        "runtime_acceptance_semantics_sha256",
    )
    if canonical_sha256(acceptance) != acceptance_id:
        raise ValueError("runtime acceptance semantic digest mismatch")
    if velocity_payload.get("runtime_acceptance_semantics_sha256") != acceptance_id:
        raise ValueError("velocity/runtime-acceptance semantic linkage mismatch")

    constrained_runtime = _require_dict(
        velocity_payload.get("constrained_callable_runtime"),
        "constrained_callable_runtime",
    )
    if constrained_runtime.get("repository") != CONSTRAINED_RUNTIME_IMPLEMENTATION_REPOSITORY:
        raise ValueError("constrained callable runtime implementation repository drift")
    implementation_commit = _require_sha1(
        constrained_runtime.get("commit"),
        "constrained callable runtime implementation commit",
    )
    if implementation_commit != CONSTRAINED_RUNTIME_IMPLEMENTATION_COMMIT:
        raise ValueError("constrained callable runtime implementation commit drift")

    if evidence.get("included_in_stable_candidate_identity") is not False:
        raise ValueError("materialization evidence leaked into stable candidate identity")
    if evidence.get("included_in_stable_velocity_identity") is not False:
        raise ValueError("materialization evidence leaked into stable velocity identity")
    if receipt.get("materialization_evidence_sha256") != canonical_sha256(evidence):
        raise ValueError("materialization evidence checksum mismatch")
    if migration.get("legacy_identity_reused_as_stable_identity") is not False:
        raise ValueError("legacy materialization identity was relabelled stable")
    if truth != TRUTH_BOUNDARY:
        raise ValueError("stable-identity truth boundary drift")

    recorded = receipt.get("receipt_sha256")
    if recorded is not None:
        unsigned = dict(receipt)
        unsigned.pop("receipt_sha256")
        if recorded != canonical_sha256(unsigned):
            raise ValueError("stable-identity receipt checksum mismatch")


def execute(*, constrained_root: str | Path, bundle_dir: str | Path) -> dict[str, Any]:
    """Verify the real bundle/acceptance contract, then emit the stable identity split."""
    constrained_root = Path(constrained_root).resolve()
    bundle_dir = Path(bundle_dir).resolve()
    _authenticate_constrained_runtime_checkout(constrained_root)
    bridge._install_constrained_package_path(constrained_root)
    acceptance_mod = importlib.import_module(
        "openai_ns_reconstruction.constrained_st052_runtime_dependency_acceptance"
    )
    whole_mod = importlib.import_module("openai_ns_reconstruction.st052_linear_temporal_capsule")

    # These calls retain the old byte-level and runtime-contract authentication.
    acceptance_mod.audit(constrained_root)
    candidate_manifest = whole_mod.verify_bundle(bundle_dir)
    contract_path = constrained_root / bridge.DEPENDENCY_CONTRACT_REL
    dependency_contract = json.loads(contract_path.read_text(encoding="utf-8"))
    return build_receipt(
        candidate_manifest,
        dependency_contract,
        dependency_contract_file_sha256=_sha256_file(contract_path),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--constrained-root", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    receipt = execute(constrained_root=args.constrained_root, bundle_dir=args.bundle_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "candidate_semantic_identity_sha256": receipt["stable_identity"][
                    "candidate_semantic_identity_sha256"
                ],
                "velocity_semantic_identity_sha256": receipt["stable_identity"][
                    "velocity_semantic_identity_sha256"
                ],
                "legacy_whole_candidate_identity_sha256": receipt[
                    "materialization_evidence"
                ]["legacy_whole_candidate_identity_sha256"],
                "materialization_file_hashes_identity_bound": receipt["truth_boundary"][
                    "materialization_file_hashes_identity_bound"
                ],
                "visualization_ready": receipt["truth_boundary"]["visualization_ready"],
                "pde_validated": receipt["truth_boundary"]["pde_validated"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
