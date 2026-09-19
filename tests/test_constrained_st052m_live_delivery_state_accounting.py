from __future__ import annotations

from copy import deepcopy
import json

import pytest

from openai_ns_reconstruction.constrained_st052m_live_delivery_state_accounting import (
    CONTRACT_PATH,
    ROOT,
    _audit_cr001,
    audit_repository,
    load_contract,
    validate_contract,
)


def test_live_repository_accounting_audit_passes() -> None:
    receipt = audit_repository()
    assert receipt["st052_source_runtime_backed_callable_save_load_ready"] is True
    assert receipt["project_status_st052_materialization_fields_stale"] is True
    assert receipt["canonical_eq45_velocity_export_ready"] is True
    assert receipt["st052_velocity_export_ready"] is False
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["pde_validated"] is False
    assert receipt["paper_exact"] is False
    assert receipt["openai_field_identified"] is False
    assert receipt["cr001_unchanged"] is True


@pytest.mark.parametrize(
    "key",
    [
        "st052_velocity_export_ready",
        "st052_production_candidate_selected",
        "st052_visualization_ready",
        "st052_visual_correspondence_verified",
        "st052_held_out_temporal_child_pde_residual_evaluated",
        "st052_pde_validated",
        "st052_source_correspondence_verified",
        "st052_paper_exact",
        "st052_openai_field_identified",
        "st052_blowup_proved",
    ],
)
def test_premature_st052_claim_promotion_fails_closed(key: str) -> None:
    contract = deepcopy(load_contract())
    contract["claim_states"][key] = True
    with pytest.raises(ValueError, match="premature ST052 claim promotion"):
        validate_contract(contract)


def test_source_runtime_backed_callable_fact_cannot_be_erased() -> None:
    contract = deepcopy(load_contract())
    contract["claim_states"]["st052_source_runtime_backed_callable_save_load_ready"] = False
    with pytest.raises(ValueError, match="callable/save-load fact missing"):
        validate_contract(contract)


def test_whole_child_materialization_fact_cannot_be_reverted_to_old_status() -> None:
    contract = deepcopy(load_contract())
    contract["integrated_st052_runtime"]["whole_child_bundle_materialized"] = False
    with pytest.raises(ValueError, match="integrated positive fact missing"):
        validate_contract(contract)


def test_source_classes_must_remain_distinct() -> None:
    contract = deepcopy(load_contract())
    contract["source_classification"][3]["classification"] = "autonomous_design"
    with pytest.raises(ValueError, match="source classification classes drifted"):
        validate_contract(contract)


def test_cr001_momentum_threshold_mutation_is_rejected() -> None:
    contract = load_contract()
    constraints = json.loads((ROOT / "configs" / "constraints.json").read_text(encoding="utf-8"))
    mutated = deepcopy(constraints)
    mutated["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="CR001 threshold drifted: pde_residual_max"):
        _audit_cr001(contract, mutated)


def test_cr001_free_force_mutation_is_rejected() -> None:
    contract = load_contract()
    constraints = json.loads((ROOT / "configs" / "constraints.json").read_text(encoding="utf-8"))
    mutated = deepcopy(constraints)
    mutated["forcing"]["restriction"] = "Allow residual-defined pointwise force."
    with pytest.raises(ValueError, match="free residual-defined forcing prohibition missing"):
        _audit_cr001(contract, mutated)


def test_transitional_status_observation_is_explicitly_versioned() -> None:
    contract = load_contract(CONTRACT_PATH)
    observation = contract["status_reporting_observation"]
    assert observation["project_status_st052_materialization_fields_stale"] is True
    assert "deliberately revised" in observation["upgrade_rule"]
