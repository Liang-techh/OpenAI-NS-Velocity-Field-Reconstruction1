from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st051b_sampled_grid_piola_governance import (
    audit_sampled_grid_piola_governance,
    load_contract,
    load_constraints,
    load_project_status,
)


def test_sampled_grid_piola_governance_passes_frozen_contract():
    receipt = audit_sampled_grid_piola_governance()
    assert receipt["audit"] == "passed"
    assert receipt["sampled_piola_is_distinct_numerical_representation"] is True
    assert receipt["continuum_divergence_transfer_allowed"] is False
    assert receipt["sampled_energy_is_cr001_energy_acceptance"] is False
    assert receipt["nested_17_is_continuum_convergence"] is False
    assert receipt["canonical_velocity_export_ready"] is True
    assert receipt["pde_validated"] is False


def test_rejects_continuum_divergence_laundering_into_sampled_interpolant():
    contract = load_contract()
    contract["representation_layers"]["sampled_piola_child"][
        "continuum_divergence_identity_is_discrete_divergence_certificate"
    ] = True
    with pytest.raises(ValueError, match="continuum-to-discrete divergence transfer"):
        audit_sampled_grid_piola_governance(contract=contract)


def test_rejects_nested_subsample_being_promoted_to_continuum_convergence():
    contract = load_contract()
    contract["representation_layers"]["nested_coarse_view"]["continuum_convergence_study"] = True
    with pytest.raises(ValueError, match="coarse continuum convergence"):
        audit_sampled_grid_piola_governance(contract=contract)


def test_rejects_sampled_energy_restoration_as_cr001_energy_acceptance():
    contract = load_contract()
    contract["observed_pr499_evidence"]["sampled_energy_restoration_is_cr001_reference_energy_validation"] = True
    with pytest.raises(ValueError, match="sampled energy/CR001 energy transfer"):
        audit_sampled_grid_piola_governance(contract=contract)


def test_rejects_beta_only_identity_and_production_promotion():
    contract = load_contract()
    contract["representation_layers"]["sampled_piola_child"]["beta_alone_identifies_child"] = True
    with pytest.raises(ValueError, match="beta-only child identity"):
        audit_sampled_grid_piola_governance(contract=contract)

    contract = load_contract()
    contract["observed_pr499_evidence"]["capacity_pass_selects_production_beta"] = True
    with pytest.raises(ValueError, match="capacity production-beta selection"):
        audit_sampled_grid_piola_governance(contract=contract)


def test_rejects_cr001_threshold_drift_and_free_force_laundering():
    constraints = load_constraints()
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="momentum max"):
        audit_sampled_grid_piola_governance(constraints=constraints)

    constraints = load_constraints()
    constraints["forcing"]["restriction"] = "fit arbitrary pointwise force"
    with pytest.raises(ValueError, match="free-force guard"):
        audit_sampled_grid_piola_governance(constraints=constraints)


def test_rejects_canonical_delivery_state_promotion():
    status = deepcopy(load_project_status())
    status["states"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError, match="project-status visual_correspondence_verified"):
        audit_sampled_grid_piola_governance(project_status=status)
