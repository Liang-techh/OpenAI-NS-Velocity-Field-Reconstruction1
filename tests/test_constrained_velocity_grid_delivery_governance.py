from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_velocity_grid_delivery_governance import (
    audit_velocity_grid_delivery_identity,
    load_constraints,
    load_contract,
    load_project_status,
)


def _inputs():
    return load_contract(), load_constraints(), load_project_status()


def test_velocity_grid_delivery_identity_contract_passes_current_live_state():
    result = audit_velocity_grid_delivery_identity()
    assert result["status"] == "governance_pass"
    assert result["analytic_candidate"] == "ST048-S"
    assert result["sampled_grid_is_separate_identity"] is True
    assert result["sampled_grid_is_not_continuum_velocity_api"] is True
    assert result["st048_canonical_promotion"] is False
    assert result["grid_derivatives_are_pde_evidence"] is False
    assert result["canonical_thresholds_unchanged"] is True
    assert result["pde_validated"] is False
    assert result["visual_correspondence_verified"] is False


def test_grid_identity_cannot_be_laundered_into_candidate_identity():
    contract, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["identity_layers"]["sampled_grid"]["grid_sha_is_candidate_sha"] = True
    with pytest.raises(ValueError, match="grid/candidate SHA separation"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)


def test_sampled_grid_cannot_be_promoted_to_continuum_api_or_exact_openai_field():
    contract, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["identity_layers"]["sampled_grid"]["sampled_grid_is_continuum_velocity_api"] = True
    with pytest.raises(ValueError, match="grid/continuum separation"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)

    bad = deepcopy(contract)
    bad["grid_evidence_scope"]["sampled_grid_may_be_called_exact_openai_field"] = True
    with pytest.raises(ValueError, match="sampled_grid_may_be_called_exact_openai_field"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)


def test_open_grid_export_cannot_silently_replace_canonical_velocity_delivery():
    contract, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["identity_layers"]["canonical_callable_delivery"][
        "open_unintegrated_grid_export_changes_canonical_velocity_api"
    ] = True
    with pytest.raises(ValueError, match="grid/API promotion separation"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)

    bad = deepcopy(contract)
    bad["identity_layers"]["canonical_callable_delivery"]["velocity_api"] = "st048_grid.interpolate"
    with pytest.raises(ValueError, match="canonical velocity API"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)


def test_grid_sampling_cannot_replace_cr001_validation_protocol():
    contract, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["grid_evidence_scope"]["grid_nodes_may_replace_cr001_held_out_sample"] = True
    with pytest.raises(ValueError, match="grid_nodes_may_replace_cr001_held_out_sample"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)

    bad = deepcopy(contract)
    bad["canonical_cr001"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="PDE max"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)


def test_source_classes_and_scientific_truth_states_fail_closed():
    contract, constraints, project_status = _inputs()
    bad = deepcopy(contract)
    bad["source_classification"]["openai_hidden_numerical_velocity"] = "public_source_fact"
    with pytest.raises(ValueError, match="source class openai_hidden_numerical_velocity"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)

    bad = deepcopy(contract)
    bad["truth_boundary"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError, match="contract truth visual_correspondence_verified"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)

    bad = deepcopy(contract)
    bad["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="contract truth pde_validated"):
        audit_velocity_grid_delivery_identity(bad, constraints, project_status)
