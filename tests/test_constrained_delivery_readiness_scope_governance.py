from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_delivery_readiness_scope_governance import (
    audit_delivery_readiness_scope,
    load_constraints,
    load_contract,
    load_delivery_state_contract,
    load_grid_identity_contract,
    load_project_status,
)


def _inputs():
    return (
        load_contract(),
        load_delivery_state_contract(),
        load_grid_identity_contract(),
        load_constraints(),
        load_project_status(),
    )


def test_delivery_readiness_scope_contract_passes_current_live_state():
    result = audit_delivery_readiness_scope()
    assert result["status"] == "governance_pass"
    assert result["observed_pr"] == 476
    assert result["unqualified_velocity_export_ready_requires_unified_candidate_delivery"] is True
    assert result["sampled_grid_has_separate_readiness_scope"] is True
    assert result["pr476_unqualified_readiness_claim_contract_compatible"] is False
    assert result["pr476_candidate_velocity_export_ready"] is False
    assert result["pr476_sampled_grid_export_ready"] == "pending_exact_source_ci"
    assert result["canonical_velocity_export_ready"] is True
    assert result["canonical_thresholds_unchanged"] is True
    assert result["pde_validated"] is False
    assert result["visual_correspondence_verified"] is False


def test_sampled_grid_cannot_borrow_unqualified_velocity_export_ready():
    contract, delivery_state, grid_identity, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["observed_pr476"]["unqualified_velocity_export_ready_claim_is_contract_compatible"] = True
    with pytest.raises(ValueError, match="PR476 readiness compatibility"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)

    bad = deepcopy(contract)
    bad["readiness_semantics"]["unqualified_velocity_export_ready"]["sampled_grid_alone_is_sufficient"] = True
    with pytest.raises(ValueError, match="sampled_grid_alone_is_sufficient"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)


def test_candidate_local_points_time_callable_is_not_the_unified_velocity_api():
    contract, delivery_state, grid_identity, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["readiness_semantics"]["unqualified_velocity_export_ready"][
        "auxiliary_points_time_callable_alone_is_sufficient"
    ] = True
    with pytest.raises(ValueError, match="auxiliary_points_time_callable_alone_is_sufficient"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)

    bad = deepcopy(contract)
    bad["observed_pr476"]["registers_repository_unified_velocity_api"] = True
    with pytest.raises(ValueError, match="PR476 unified API registration"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)


def test_existing_grid_identity_contract_remains_fail_closed():
    contract, delivery_state, grid_identity, constraints, project_status = _inputs()
    bad_grid = deepcopy(grid_identity)
    bad_grid["identity_layers"]["sampled_grid"]["sampled_grid_is_continuum_velocity_api"] = True
    with pytest.raises(ValueError, match="existing grid/continuum separation"):
        audit_delivery_readiness_scope(contract, delivery_state, bad_grid, constraints, project_status)


def test_open_grid_export_cannot_replace_canonical_candidate_or_api():
    contract, delivery_state, grid_identity, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["canonical_callable_delivery"]["open_pr476_changes_canonical_velocity_api"] = True
    with pytest.raises(ValueError, match="PR476 canonical API separation"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)

    bad = deepcopy(contract)
    bad["canonical_callable_delivery"]["velocity_api"] = "st051b_grid.interpolate"
    with pytest.raises(ValueError, match="canonical velocity API"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)


def test_cr001_thresholds_and_scientific_truth_states_cannot_be_laundered():
    contract, delivery_state, grid_identity, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["canonical_cr001"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="PDE max"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)

    bad = deepcopy(contract)
    bad["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="contract truth pde_validated"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)

    bad = deepcopy(contract)
    bad["source_classification"]["openai_hidden_numerical_velocity"] = "public_source_fact"
    with pytest.raises(ValueError, match="source class openai_hidden_numerical_velocity"):
        audit_delivery_readiness_scope(bad, delivery_state, grid_identity, constraints, project_status)
