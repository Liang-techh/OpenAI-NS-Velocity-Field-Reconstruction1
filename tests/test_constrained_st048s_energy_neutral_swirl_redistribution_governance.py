from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st048s_energy_neutral_swirl_redistribution_governance import (
    audit_energy_neutral_swirl_redistribution_governance,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_redistribution_governance_passes_registered_contract():
    result = audit_energy_neutral_swirl_redistribution_governance()
    assert result["status"] == "pass"
    assert result["first_order_energy_neutrality_only"] is True
    assert result["diagnostic_crossing_only"] is True
    assert result["material_path_gain_verified"] is False
    assert result["canonical_velocity_unchanged"] is True
    assert result["live_route_unchanged"] is True
    assert result["visual_correspondence_verified"] is False
    assert result["pde_validated"] is False


def test_first_order_energy_neutrality_cannot_be_upgraded_to_finite_energy_invariance():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["energy_semantics"]["first_order_energy_neutral_means_exact_energy_invariant_for_finite_gain"] = True
    with pytest.raises(ValueError, match="first_order_energy_neutral_means_exact_energy_invariant_for_finite_gain"):
        audit_energy_neutral_swirl_redistribution_governance(scope=mutated)


def test_balance_coefficient_cannot_be_relabelled_source_recovered():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["source_classification"]["reference_time_balance_coefficient_alpha"] = "public_source_fact"
    with pytest.raises(ValueError, match="reference_time_balance_coefficient_alpha"):
        audit_energy_neutral_swirl_redistribution_governance(scope=mutated)


def test_diagnostic_crossing_cannot_select_production_gain():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["redistribution_definition"]["production_coefficient_selected"] = True
    with pytest.raises(ValueError, match="production coefficient selection"):
        audit_energy_neutral_swirl_redistribution_governance(scope=mutated)


def test_fixed_probe_gain_cannot_be_promoted_to_material_path_gain():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["promotion_contract"]["fixed_probe_angular_rate_gain_may_be_called_material_path_gain"] = True
    with pytest.raises(ValueError, match="fixed_probe_angular_rate_gain_may_be_called_material_path_gain"):
        audit_energy_neutral_swirl_redistribution_governance(scope=mutated)


def test_st006_radius_attribution_cannot_be_promoted_to_openai_truth():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["radius_attribution_semantics"]["st006_is_openai_visual_truth"] = True
    with pytest.raises(ValueError, match="ST006/OpenAI truth separation"):
        audit_energy_neutral_swirl_redistribution_governance(scope=mutated)


def test_rank_novelty_cannot_be_promoted_to_pde_leverage():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["sensitivity_semantics"]["rank_four_proves_pde_leverage"] = True
    with pytest.raises(ValueError, match="rank_four_proves_pde_leverage"):
        audit_energy_neutral_swirl_redistribution_governance(scope=mutated)


def test_materialized_child_cannot_drop_multiplier_positivity_guard():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["representation_and_axis_contract"]["materialized_child_requires_multiplier_finite_and_positive_on_support"] = False
    with pytest.raises(ValueError, match="materialized_child_requires_multiplier_finite_and_positive_on_support"):
        audit_energy_neutral_swirl_redistribution_governance(scope=mutated)


def test_cr001_threshold_drift_is_rejected():
    constraints = load_constraints()
    mutated = deepcopy(constraints)
    mutated["validation"]["thresholds"]["pde_residual_L2"] = 0.01
    with pytest.raises(ValueError, match="pde_residual_L2"):
        audit_energy_neutral_swirl_redistribution_governance(constraints=mutated)


def test_project_truth_state_promotion_is_rejected():
    status = load_project_status()
    mutated = deepcopy(status)
    mutated["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="project state pde_validated"):
        audit_energy_neutral_swirl_redistribution_governance(project_status=mutated)
