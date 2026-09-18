from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_render_readiness_scope_governance import (
    audit_render_readiness_scope,
    load_contract,
    load_delivery_state_contract,
)


def test_render_readiness_scope_contract_passes():
    result = audit_render_readiness_scope()
    assert result["status"] == "governance_pass"
    assert result["integrated_render_pr"] == 490
    assert result["render_artifact_ready"] is True
    assert result["visualization_ready"] is False
    assert result["visual_correspondence_verified"] is False
    assert result["instantaneous_streamlines_are_material_paths"] is False
    assert result["sampled_vorticity_is_pde_acceptance_evidence"] is False
    assert result["canonical_velocity_export_ready"] is True
    assert result["canonical_thresholds_unchanged"] is True
    assert result["pde_validated"] is False


def test_rejects_render_artifact_promotion_to_visualization_ready():
    contract = deepcopy(load_contract())
    contract["integrated_pr490"]["unqualified_visualization_ready"] = True
    with pytest.raises(ValueError, match="unqualified_visualization_ready"):
        audit_render_readiness_scope(contract=contract)


def test_rejects_instantaneous_streamline_material_path_laundering():
    contract = deepcopy(load_contract())
    contract["integrated_pr490"]["streamline_may_be_claimed_as_material_path"] = True
    with pytest.raises(ValueError, match="streamline_may_be_claimed_as_material_path"):
        audit_render_readiness_scope(contract=contract)


def test_rejects_sampled_vorticity_pde_laundering():
    contract = deepcopy(load_contract())
    contract["integrated_pr490"]["sampled_vorticity_may_be_claimed_as_pde_acceptance_evidence"] = True
    with pytest.raises(ValueError, match="sampled_vorticity_may_be_claimed_as_pde_acceptance_evidence"):
        audit_render_readiness_scope(contract=contract)


def test_rejects_missing_visual_resolution_sanity_evidence():
    delivery_state = deepcopy(load_delivery_state_contract())
    delivery_state["states"]["visualization_ready"]["positive_evidence"] = [
        "callable_velocity",
        "supported_visualization_entry_point",
    ]
    with pytest.raises(ValueError, match="visualization-ready evidence"):
        audit_render_readiness_scope(delivery_state_contract=delivery_state)


def test_rejects_target_free_render_as_visual_correspondence():
    contract = deepcopy(load_contract())
    contract["readiness_semantics"]["visual_correspondence_verified"]["target_free_render_is_sufficient"] = True
    with pytest.raises(ValueError, match="target_free_render_is_sufficient"):
        audit_render_readiness_scope(contract=contract)


def test_rejects_cr001_threshold_drift():
    contract = deepcopy(load_contract())
    contract["canonical_cr001"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="PDE max"):
        audit_render_readiness_scope(contract=contract)
