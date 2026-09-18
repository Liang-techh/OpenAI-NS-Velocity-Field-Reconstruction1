from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st048s_material_path_evidence_governance import (
    audit_material_path_evidence_governance,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_material_path_governance_passes_on_registered_contract():
    result = audit_material_path_evidence_governance()
    assert result["status"] == "pass"
    assert result["material_path_protocol_registered"] is True
    assert result["canonical_velocity_unchanged"] is True
    assert result["live_route_unchanged"] is True
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False


def test_material_path_protocol_drift_fails_closed():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["material_path_contract"]["rtol"] = 1e-6
    with pytest.raises(ValueError, match="pathline rtol"):
        audit_material_path_evidence_governance(scope=mutated)


def test_proxy_cannot_be_promoted_to_material_path_evidence():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["observable_semantics"]["fixed_probe_eulerian_proxy_is_material_path_evidence"] = True
    with pytest.raises(ValueError, match="fixed_probe_eulerian_proxy_is_material_path_evidence"):
        audit_material_path_evidence_governance(scope=mutated)


def test_ring_proxy_failure_cannot_select_global_kappa():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["pr437_proxy_receipt"]["proxy_failure_promotes_global_kappa_to_production"] = True
    with pytest.raises(ValueError, match="proxy_failure_promotes_global_kappa_to_production"):
        audit_material_path_evidence_governance(scope=mutated)


def test_openai_material_trajectory_truth_stays_pending():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["source_classification"]["openai_hidden_particle_seeds_or_material_trajectories"] = "public_source_fact"
    with pytest.raises(ValueError, match="openai_hidden_particle_seeds_or_material_trajectories"):
        audit_material_path_evidence_governance(scope=mutated)


def test_material_path_tradeoff_cannot_be_relabelled_visual_correspondence():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["pr435_diagnostic_receipt"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError, match="PR435 visual_correspondence_verified"):
        audit_material_path_evidence_governance(scope=mutated)


def test_material_path_evidence_cannot_replace_canonical_route():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["promotion_contract"]["material_path_or_proxy_improvement_may_replace_live_scientific_route"] = True
    with pytest.raises(ValueError, match="material_path_or_proxy_improvement_may_replace_live_scientific_route"):
        audit_material_path_evidence_governance(scope=mutated)


def test_cr001_threshold_drift_is_rejected():
    constraints = load_constraints()
    mutated = deepcopy(constraints)
    mutated["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="pde_residual_max"):
        audit_material_path_evidence_governance(constraints=mutated)


def test_project_truth_state_promotion_is_rejected():
    status = load_project_status()
    mutated = deepcopy(status)
    mutated["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="project state pde_validated"):
        audit_material_path_evidence_governance(project_status=mutated)
