from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_st052m_parent_capsule_runtime_scope_governance import (
    CAPSULE_REL,
    CONSTRAINTS_REL,
    CONTRACT_REL,
    STATUS_REL,
    WORKFLOW_REL,
    audit_repository,
    validate_contract_data,
)

ROOT = Path(__file__).resolve().parents[1]


def _inputs():
    contract = json.loads((ROOT / CONTRACT_REL).read_text())
    constraints = json.loads((ROOT / CONSTRAINTS_REL).read_text())
    status = json.loads((ROOT / STATUS_REL).read_text())
    capsule_source = (ROOT / CAPSULE_REL).read_text()
    workflow_source = (ROOT / WORKFLOW_REL).read_text()
    return contract, constraints, status, capsule_source, workflow_source


def test_live_repository_audit_passes_and_keeps_delivery_stages_separate():
    receipt = audit_repository(ROOT)
    assert receipt == {
        "schema": "st052m-parent-capsule-runtime-scope-governance/v1",
        "candidate_id": "ST052-M",
        "builder_source_head_verified": False,
        "workflow_exact_source_head_bound": True,
        "parent_runtime_ready": False,
        "whole_child_save_load_ready": False,
        "experimental_st052_velocity_export_ready": False,
        "canonical_velocity_export_ready": True,
        "pde_validated": False,
    }


def test_rejects_builder_source_head_laundering():
    contract, constraints, status, capsule, workflow = _inputs()
    bad = copy.deepcopy(contract)
    bad["observed_live_scope"]["source_head_verified_inside_builder"] = True
    with pytest.raises(ValueError, match="source_head_verified_inside_builder"):
        validate_contract_data(bad, constraints, status, capsule, workflow)


def test_rejects_workflow_artifact_as_parent_runtime_promotion():
    contract, constraints, status, capsule, workflow = _inputs()
    bad = copy.deepcopy(contract)
    bad["truth_states"]["parent_runtime_ready"] = True
    with pytest.raises(ValueError, match="parent_runtime_ready"):
        validate_contract_data(bad, constraints, status, capsule, workflow)


def test_rejects_workflow_artifact_as_whole_child_save_load():
    contract, constraints, status, capsule, workflow = _inputs()
    bad = copy.deepcopy(contract)
    bad["truth_states"]["whole_child_save_load_ready"] = True
    with pytest.raises(ValueError, match="whole_child_save_load_ready"):
        validate_contract_data(bad, constraints, status, capsule, workflow)


def test_rejects_experimental_velocity_export_promotion_from_parent_capsule():
    contract, constraints, status, capsule, workflow = _inputs()
    bad = copy.deepcopy(contract)
    bad["truth_states"]["experimental_st052_velocity_export_ready"] = True
    with pytest.raises(ValueError, match="experimental_st052_velocity_export_ready"):
        validate_contract_data(bad, constraints, status, capsule, workflow)


def test_rejects_exact_checkout_ref_drift():
    contract, constraints, status, capsule, workflow = _inputs()
    changed_workflow = workflow.replace(
        "ref: b3b8bfdbe1077f9ec967d158602951997d81e17d",
        "ref: 0000000000000000000000000000000000000000",
    )
    with pytest.raises(ValueError, match="exact-source checkout ref drift"):
        validate_contract_data(contract, constraints, status, capsule, changed_workflow)


def test_rejects_cr001_threshold_drift():
    contract, constraints, status, capsule, workflow = _inputs()
    bad_constraints = copy.deepcopy(constraints)
    bad_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="pde_residual_max drift"):
        validate_contract_data(contract, bad_constraints, status, capsule, workflow)


def test_rejects_visual_or_pde_truth_promotion():
    contract, constraints, status, capsule, workflow = _inputs()
    for key in ("visual_correspondence_verified", "pde_validated", "paper_exact", "openai_field_identified"):
        bad = copy.deepcopy(contract)
        bad["truth_states"][key] = True
        with pytest.raises(ValueError, match=key):
            validate_contract_data(bad, constraints, status, capsule, workflow)
