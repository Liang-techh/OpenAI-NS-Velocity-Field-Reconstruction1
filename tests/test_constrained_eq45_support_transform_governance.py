from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_support_transform_governance import (
    audit_documents,
    audit_eq45_support_transform_evidence_scope,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(relative: str) -> dict:
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _documents() -> tuple[dict, dict, dict, dict]:
    return (
        _load("configs/eq45_support_transform_evidence_scope.json"),
        _load("configs/constraints.json"),
        _load("artifacts/constrained/eq45_velocity_candidate_seed.json"),
        _load("configs/delivery_state_contract.json"),
    )


def test_current_support_transform_scope_is_truth_bounded() -> None:
    report = audit_eq45_support_transform_evidence_scope(ROOT)
    assert report["parent_velocity_export_ready"] is True
    assert report["parent_requires_physical_support_connection"] is True
    assert report["production_support_transform_materialized"] is False
    assert report["child_candidate_identity"] == "pending_unknown"
    assert report["child_global_evidence_must_be_reestablished"] is True
    assert report["pde_failure_blocks_child_velocity_export"] is False
    assert report["physical_support_validated"] is False
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False
    assert report["paper_exact"] is False
    assert report["velocity_changed"] is False


def test_rejects_parent_validation_as_child_or_provenance() -> None:
    contract, constraints, candidate, delivery = _documents()

    inherited = deepcopy(contract)
    inherited["evidence_inheritance"]["provenance_only"].append("pde_validated")
    with pytest.raises(ValueError, match="cannot be inherited as provenance"):
        audit_documents(inherited, constraints, candidate, delivery)

    missing_recheck = deepcopy(contract)
    missing_recheck["evidence_inheritance"][
        "must_be_reestablished_for_transformed_candidate"
    ].remove("visual_resolution_stability")
    with pytest.raises(ValueError, match="revalidation set drifted"):
        audit_documents(missing_recheck, constraints, candidate, delivery)

    child_inherits_parent = deepcopy(contract)
    child_inherits_parent["evidence_inheritance"][
        "must_not_be_inherited_as_child_validation"
    ].remove("parent_pde_residual_metrics")
    with pytest.raises(ValueError, match="forbidden child-evidence inheritance set drifted"):
        audit_documents(child_inherits_parent, constraints, candidate, delivery)


def test_rejects_premature_support_promotion_or_source_reclassification() -> None:
    contract, constraints, candidate, delivery = _documents()

    materialized = deepcopy(contract)
    materialized["current_transform_state"]["production_support_transform_materialized"] = True
    materialized["current_transform_state"]["child_candidate_identity"] = "child-sha"
    with pytest.raises(ValueError, match="now requires child audit"):
        audit_documents(materialized, constraints, candidate, delivery)

    promoted = deepcopy(contract)
    promoted["truth_boundary"]["physical_support_validated"] = True
    with pytest.raises(ValueError, match="cannot promote readiness"):
        audit_documents(promoted, constraints, candidate, delivery)

    reclassified = deepcopy(contract)
    reclassified["classification"]["support_transform_design"] = "public_source_fact"
    with pytest.raises(ValueError, match="source classification drifted"):
        audit_documents(reclassified, constraints, candidate, delivery)


def test_rejects_free_force_or_pde_export_gate_drift() -> None:
    contract, constraints, candidate, delivery = _documents()

    free_force = deepcopy(constraints)
    free_force["forcing"]["mode"] = "residual_defined_pointwise_force"
    with pytest.raises(ValueError, match="forcing family drifted"):
        audit_documents(contract, free_force, candidate, delivery)

    blocks_export = deepcopy(contract)
    blocks_export["policy"]["pde_failure_blocks_child_velocity_export"] = True
    with pytest.raises(ValueError, match="pde_failure_blocks_child_velocity_export"):
        audit_documents(blocks_export, constraints, candidate, delivery)

    parent_claims_support = deepcopy(candidate)
    parent_claims_support["truth_boundary"]["physical_support_validated"] = True
    with pytest.raises(ValueError, match="parent candidate support/delivery state drifted"):
        audit_documents(contract, constraints, parent_claims_support, delivery)
