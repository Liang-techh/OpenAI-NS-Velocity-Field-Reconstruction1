from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_redistribution_gain_headroom_governance import (
    audit_redistribution_gain_headroom,
    load_constraints,
    load_contract,
    load_cross_backbone_contract,
    load_delivery_state_contract,
    load_project_status,
)


def _inputs():
    return (
        load_contract(),
        load_cross_backbone_contract(),
        load_delivery_state_contract(),
        load_constraints(),
        load_project_status(),
    )


def test_redistribution_gain_headroom_contract_passes_current_live_state():
    result = audit_redistribution_gain_headroom()
    assert result["status"] == "governance_pass"
    assert result["headroom_pr"] == 480
    assert result["lower_gain_with_material_paths"] == 0.025
    assert result["higher_gain_capacity_crossing"] == 0.05
    assert result["higher_gain_production_selected"] is False
    assert result["higher_gain_material_path_evidence_exists"] is False
    assert result["lower_gain_evidence_transfer_allowed"] is False
    assert result["higher_gain_requires_new_identity_and_sha"] is True
    assert result["canonical_velocity_export_ready"] is True
    assert result["canonical_thresholds_unchanged"] is True
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False


def test_clean_headroom_crossing_cannot_select_production_gain_or_bound():
    contract, cross, delivery, constraints, project = _inputs()

    bad = deepcopy(contract)
    bad["gain_evidence_semantics"]["clean_headroom_crossing_selects_production_gain"] = True
    with pytest.raises(ValueError, match="clean_headroom_crossing_selects_production_gain"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["gain_evidence_semantics"]["largest_preregistered_clean_gain_is_production_upper_bound"] = True
    with pytest.raises(ValueError, match="largest_preregistered_clean_gain_is_production_upper_bound"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)


def test_lower_gain_material_path_evidence_cannot_be_inherited_or_extrapolated():
    contract, cross, delivery, constraints, project = _inputs()

    bad = deepcopy(contract)
    bad["gain_evidence_semantics"]["lower_gain_material_paths_may_be_inherited_by_higher_gain"] = True
    with pytest.raises(ValueError, match="lower_gain_material_paths_may_be_inherited_by_higher_gain"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["gain_evidence_semantics"]["lower_gain_material_paths_may_be_linearly_extrapolated_to_higher_gain"] = True
    with pytest.raises(ValueError, match="lower_gain_material_paths_may_be_linearly_extrapolated_to_higher_gain"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["evidence_by_gain"]["gain_0_05"]["material_path_evidence_exists"] = True
    with pytest.raises(ValueError, match="higher-gain truth material_path_evidence_exists"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)


def test_lower_gain_grid_or_parent_pde_receipt_cannot_be_relabelled_at_higher_gain():
    contract, cross, delivery, constraints, project = _inputs()

    bad = deepcopy(contract)
    bad["gain_evidence_semantics"]["lower_gain_sampled_grid_may_be_relabelled_as_higher_gain"] = True
    with pytest.raises(ValueError, match="lower_gain_sampled_grid_may_be_relabelled_as_higher_gain"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["gain_evidence_semantics"]["parent_or_lower_gain_pde_receipt_may_be_inherited"] = True
    with pytest.raises(ValueError, match="parent_or_lower_gain_pde_receipt_may_be_inherited"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["evidence_by_gain"]["gain_0_05"]["fresh_full_momentum_receipt_exists"] = True
    with pytest.raises(ValueError, match="higher-gain truth fresh_full_momentum_receipt_exists"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)


def test_gain_change_requires_new_materialized_identity_and_fresh_replay():
    contract, cross, delivery, constraints, project = _inputs()

    bad = deepcopy(contract)
    bad["fixed_transform_identity"]["changing_gain_requires_new_candidate_sha256"] = False
    with pytest.raises(ValueError, match="new candidate SHA rule"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["materialization_requirements_for_gain_0_05"]["fresh_material_paths_required_before_material_path_claim"] = False
    with pytest.raises(ValueError, match="fresh_material_paths_required_before_material_path_claim"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["materialization_requirements_for_gain_0_05"][
        "fresh_compatible_pressure_and_restricted_force_required_before_pde_claim"
    ] = False
    with pytest.raises(ValueError, match="fresh_compatible_pressure_and_restricted_force_required_before_pde_claim"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)


def test_proxy_or_grid_stability_cannot_launder_visual_or_openai_identity():
    contract, cross, delivery, constraints, project = _inputs()

    bad = deepcopy(contract)
    bad["gain_evidence_semantics"]["sampled_q90_q99_stability_is_continuum_morphology_certificate"] = True
    with pytest.raises(ValueError, match="sampled_q90_q99_stability_is_continuum_morphology_certificate"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["gain_evidence_semantics"]["eulerian_proxy_improvement_is_visual_correspondence_certificate"] = True
    with pytest.raises(ValueError, match="eulerian_proxy_improvement_is_visual_correspondence_certificate"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["source_classification"]["openai_hidden_swirl_redistribution_gain"] = "public_source_fact"
    with pytest.raises(ValueError, match="source class openai_hidden_swirl_redistribution_gain"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)


def test_cr001_thresholds_and_canonical_delivery_cannot_drift():
    contract, cross, delivery, constraints, project = _inputs()

    bad = deepcopy(contract)
    bad["canonical_cr001"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="PDE max"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["canonical_delivery_and_truth"]["velocity_api"] = "st051b_gain_005.velocity"
    with pytest.raises(ValueError, match="canonical velocity API"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)

    bad = deepcopy(contract)
    bad["canonical_delivery_and_truth"]["pde_validated"] = True
    with pytest.raises(ValueError, match="contract truth pde_validated"):
        audit_redistribution_gain_headroom(bad, cross, delivery, constraints, project)
