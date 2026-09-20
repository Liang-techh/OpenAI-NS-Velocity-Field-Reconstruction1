"""Audit the explicit external-runtime acceptance for frozen ST052-M delivery.

A8-DELIVERY-01 deliberately does not vendor or reimplement the historical
ST052-M parent.  Instead it accepts the already authenticated exact #508
checkout as an explicit runtime dependency for visualization-candidate delivery
and the next end-to-end integration smoke.

This module is governance/delivery plumbing only.  It must not promote the
experimental ST052 candidate to standalone-package export readiness, visual
correspondence, PDE validity, paper exactness, or OpenAI-field identity.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from . import st052_linear_temporal_capsule as whole_capsule
from . import st052_source_runtime_identity as source_runtime

CONTRACT_REL = Path("configs/st052m_runtime_dependency_acceptance.json")
SCHEMA = "st052m-runtime-dependency-acceptance/v1"
TASK_ID = "A8-DELIVERY-01"
CANDIDATE_ID = "ST052-M-linear-temporal-child-v1"


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def _load_contract(repo_root: str | Path) -> dict[str, Any]:
    path = Path(repo_root) / CONTRACT_REL
    obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("ST052 runtime dependency contract must be a JSON object")
    return obj


def audit(repo_root: str | Path = ".") -> dict[str, Any]:
    """Fail closed unless the accepted dependency exactly matches live runtime semantics."""
    contract = _load_contract(repo_root)
    if contract.get("schema") != SCHEMA:
        raise ValueError("runtime dependency contract schema mismatch")
    if contract.get("task_id") != TASK_ID:
        raise ValueError("runtime dependency contract task mismatch")
    if contract.get("candidate_id") != CANDIDATE_ID:
        raise ValueError("runtime dependency contract candidate mismatch")
    if contract.get("decision") != "accept_authenticated_exact_source_checkout_as_external_visualization_runtime_dependency":
        raise ValueError("runtime dependency acceptance decision changed")

    required = contract.get("required_source_runtime")
    policy = contract.get("runtime_policy")
    truth = contract.get("truth_boundary")
    if not all(isinstance(v, dict) for v in (required, policy, truth)):
        raise ValueError("runtime dependency contract sections malformed")

    expected_runtime_sha = source_runtime.source_runtime_identity_sha256()
    expected = {
        "source_head": source_runtime.SOURCE_HEAD,
        "source_tree": source_runtime.SOURCE_TREE,
        "recipe_path": source_runtime.SOURCE_RECIPE_PATH,
        "recipe_git_blob_sha1": source_runtime.SOURCE_RECIPE_GIT_BLOB_SHA1,
        "source_runtime_identity_sha256": expected_runtime_sha,
    }
    for key, value in expected.items():
        if required.get(key) != value:
            raise ValueError(f"accepted exact-source runtime drift: {key}")
    if required.get("worktree_must_authenticate") is not True:
        raise ValueError("exact-source worktree authentication cannot be relaxed")
    if required.get("historical_module_cache_must_be_isolated") is not True:
        raise ValueError("historical module-cache isolation cannot be relaxed")

    if policy.get("exact_source_checkout_dependency_explicitly_accepted") is not True:
        raise ValueError("external exact-source dependency is not explicitly accepted")
    if policy.get("normal_package_plus_exact_source_checkout_load_path_accepted") is not True:
        raise ValueError("normal package + exact source checkout load path must remain explicit")
    if policy.get("vendoring_parent_into_installable_package_required_for_visualization_candidate") is not False:
        raise ValueError("this contract intentionally rejects vendoring as a visualization blocker")
    if policy.get("integration_smoke_must_record_resolved_python_numpy_scipy_sympy_versions") is not True:
        raise ValueError("next integration smoke must record the resolved numerical runtime")
    if policy.get("resolved_dependency_versions_are_part_of_frozen_candidate_identity") is not False:
        raise ValueError("runtime receipt must not silently mutate the frozen candidate identity")
    if policy.get("portable_bitwise_identity_across_dependency_builds_claimed") is not False:
        raise ValueError("portable bitwise runtime identity remains unproved")
    if policy.get("standalone_package_parent_runtime_ready") is not False:
        raise ValueError("standalone package parent runtime is still not ready")

    live = whole_capsule.TRUTH_BOUNDARY
    required_live_true = {
        "whole_child_bundle_materialized": True,
        "whole_child_save_load_ready_with_exact_source_runtime": True,
        "exact_source_runtime_identity_closed": True,
        "historical_module_cache_isolated": True,
    }
    for key, value in required_live_true.items():
        if live.get(key) is not value:
            raise ValueError(f"live whole-candidate positive fact changed: {key}")
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
        if live.get(key) is not False:
            raise ValueError(f"live whole-candidate truth state must remain false: {key}")
        if truth.get(key) is not False:
            raise ValueError(f"acceptance contract truth state must remain false: {key}")

    if truth.get("external_runtime_dependency_accepted_for_visualization_delivery") is not True:
        raise ValueError("accepted dependency must be machine-visible")
    if contract.get("unblocks") != ["A8-DELIVERY-02"]:
        raise ValueError("runtime acceptance may unblock only the next delivery smoke")

    contract_sha = hashlib.sha256(_canonical_bytes(contract)).hexdigest()
    return {
        "schema": SCHEMA,
        "task_id": TASK_ID,
        "candidate_id": CANDIDATE_ID,
        "contract_sha256": contract_sha,
        "exact_source_runtime_identity_sha256": expected_runtime_sha,
        "external_runtime_dependency_accepted_for_visualization_delivery": True,
        "standalone_package_parent_runtime_ready": False,
        "velocity_export_ready": False,
        "visualization_ready": False,
        "pde_validated": False,
        "next_task": "A8-DELIVERY-02",
    }


if __name__ == "__main__":
    print(json.dumps(audit(Path(__file__).resolve().parents[2]), indent=2, sort_keys=True))
