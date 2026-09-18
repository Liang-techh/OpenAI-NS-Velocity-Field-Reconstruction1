from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_cross_backbone_swirl_redistribution_governance import (
    audit_cross_backbone_swirl_redistribution_governance,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_cross_backbone_governance_passes_on_registered_contract():
    result = audit_cross_backbone_swirl_redistribution_governance()
    assert result["audit"] == "pass"
    assert result["same_numeric_gain_is_portable_across_parents"] is False
    assert result["st050rc_material_path_receipt_exists_at_audited_head"] is False
    assert result["canonical_velocity_unchanged"] is True
    assert result["cr001_unchanged"] is True


def test_parent_conditioned_alpha_cannot_be_erased():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["evidence_snapshots"]["st050rc"]["balance_alpha"] = mutated["evidence_snapshots"]["st048s"]["balance_alpha"]
    with pytest.raises(ValueError, match="ST050 alpha|alpha difference"):
        audit_cross_backbone_swirl_redistribution_governance(scope=mutated)


def test_same_numeric_gain_cannot_be_declared_portable():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["representation_identity"]["same_numeric_kappa_is_portable_production_parameter_across_parent_families"] = True
    with pytest.raises(ValueError, match="same_numeric_kappa_is_portable"):
        audit_cross_backbone_swirl_redistribution_governance(scope=mutated)


def test_rebalanced_cross_parent_comparison_cannot_be_single_factor():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["cross_backbone_comparability"]["cross_parent_child_difference_attributable_only_to_backbone_when_alpha_rebalanced"] = True
    with pytest.raises(ValueError, match="cross_parent_child_difference_attributable_only_to_backbone"):
        audit_cross_backbone_swirl_redistribution_governance(scope=mutated)


def test_two_clean_crossings_cannot_select_production_gain():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["cross_backbone_comparability"]["same_clean_crossing_on_two_parents_selects_production_gain"] = True
    with pytest.raises(ValueError, match="same_clean_crossing_on_two_parents_selects_production_gain"):
        audit_cross_backbone_swirl_redistribution_governance(scope=mutated)


def test_green_standard_ci_cannot_replace_failed_exact_source_replay():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["ci_evidence_semantics"]["green_standard_ci_may_replace_missing_dedicated_scientific_receipt"] = True
    with pytest.raises(ValueError, match="green_standard_ci_may_replace"):
        audit_cross_backbone_swirl_redistribution_governance(scope=mutated)


def test_failed_replay_cannot_be_relabelled_scientific_negative_result():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["ci_evidence_semantics"]["failed_dedicated_replay_before_numerical_receipt_is_scientific_negative_result"] = True
    with pytest.raises(ValueError, match="scientific_negative_result"):
        audit_cross_backbone_swirl_redistribution_governance(scope=mutated)


def test_openai_hidden_swirl_profile_stays_pending():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["source_classification"]["openai_hidden_radial_swirl_redistribution"] = "public_source_fact"
    with pytest.raises(ValueError, match="openai_hidden_radial_swirl_redistribution"):
        audit_cross_backbone_swirl_redistribution_governance(scope=mutated)


def test_cr001_threshold_drift_is_rejected():
    constraints = load_constraints()
    mutated = deepcopy(constraints)
    mutated["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="pde_residual_max"):
        audit_cross_backbone_swirl_redistribution_governance(constraints=mutated)


def test_project_pde_promotion_is_rejected():
    status = load_project_status()
    mutated = deepcopy(status)
    mutated["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="project state pde_validated"):
        audit_cross_backbone_swirl_redistribution_governance(project_status=mutated)
