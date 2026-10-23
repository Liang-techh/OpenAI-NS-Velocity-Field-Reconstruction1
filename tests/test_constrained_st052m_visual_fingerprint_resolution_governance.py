from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_visual_fingerprint_resolution_governance import (
    audit_contract,
    audit_repository,
    load_constraints,
    load_contract,
    load_delivery,
    load_project_status,
)


def _inputs():
    return load_contract(), load_constraints(), load_delivery(), load_project_status()


def _audit_mutated(mutator):
    contract, constraints, delivery, status = _inputs()
    mutated = deepcopy(contract)
    mutator(mutated)
    with pytest.raises(ValueError):
        audit_contract(mutated, constraints, delivery, status)


def test_repository_visual_fingerprint_resolution_governance_passes():
    receipt = audit_repository()
    assert receipt["status"] == "PASS"
    assert receipt["scientific_verdict"] == "FAIL"
    assert receipt["fine_grid_magnitude_guards_pass"] is True
    assert receipt["cross_resolution_sign_guard_pass"] is False
    assert receipt["stable_tip_band_sign"] is True
    assert receipt["top_cloud_radial_resolution_stable"] is False
    assert receipt["continuum_convergence_certified"] is False
    assert receipt["child_velocity_export_ready"] is False
    assert receipt["canonical_velocity_export_ready"] is True
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["pde_validated"] is False


def test_fine_grid_pass_cannot_override_preregistered_aggregate_fail():
    _audit_mutated(lambda c: c["preregistered_rule"].__setitem__("aggregate_clean_verdict", True))
    _audit_mutated(lambda c: c["resolution_semantics"].__setitem__("fine_grid_pass_can_override_preregistered_aggregate_fail", True))
    _audit_mutated(lambda c: c["scientific_scope"].__setitem__("scientific_verdict", "PASS"))


def test_two_fixed_grids_cannot_be_promoted_to_continuum_or_resolution_stability():
    _audit_mutated(lambda c: c["resolution_semantics"].__setitem__("two_grid_comparison_is_continuum_convergence_certificate", True))
    _audit_mutated(lambda c: c["resolution_semantics"].__setitem__("two_grid_comparison_is_resolution_stability_certificate_for_top_cloud_radial_metric", True))
    _audit_mutated(lambda c: c["resolution_semantics"].__setitem__("top_cloud_radial_tightening_has_same_desired_sign_on_both_grids_all_times", True))


def test_observed_coarse_fine_sign_disagreement_is_fail_closed():
    _audit_mutated(lambda c: c["coarse_grid_evidence"]["top_radial_relative_changes"].__setitem__(2, -0.00475))
    _audit_mutated(lambda c: c["coarse_grid_evidence"].__setitem__("top_cloud_radial_sign_flips_at_t_0p75", False))
    _audit_mutated(lambda c: c["preregistered_rule"].__setitem__("desired_sign_agreement_required_on_both_grids", False))


def test_stable_tip_band_sign_remains_descriptive_not_correspondence_or_pde_evidence():
    _audit_mutated(lambda c: c["resolution_semantics"].__setitem__("stable_tip_band_sign_is_visual_correspondence_evidence", True))
    _audit_mutated(lambda c: c["resolution_semantics"].__setitem__("stable_tip_band_sign_is_pde_evidence", True))
    _audit_mutated(lambda c: c["scientific_scope"].__setitem__("visual_correspondence_verified", True))
    _audit_mutated(lambda c: c["scientific_scope"].__setitem__("pde_validated", True))


def test_target_free_protocol_cannot_be_relabelled_openai_acceptance():
    _audit_mutated(lambda c: c["protocol"].__setitem__("openai_image_or_numeric_target_used", True))
    _audit_mutated(lambda c: c["protocol"].__setitem__("visual_acceptance_threshold_against_openai_used", True))
    _audit_mutated(lambda c: c["scientific_scope"].__setitem__("production_visual_fingerprint_selected", True))


def test_visual_receipt_does_not_create_delivery_identity():
    _audit_mutated(lambda c: c["child_delivery_scope"].__setitem__("versioned_child_materialized", True))
    _audit_mutated(lambda c: c["child_delivery_scope"].__setitem__("save_load_contract_exists", True))
    _audit_mutated(lambda c: c["child_delivery_scope"].__setitem__("velocity_export_ready", True))
    _audit_mutated(lambda c: c["canonical_delivery"].__setitem__("visual_fingerprint_child_replaces_canonical_candidate", True))


def test_cr001_thresholds_force_and_validation_contract_remain_unchanged():
    _audit_mutated(lambda c: c["canonical_cr001"].__setitem__("pde_residual_max", 0.01))
    _audit_mutated(lambda c: c["canonical_cr001"].__setitem__("divergence_max", 1e-4))
    _audit_mutated(lambda c: c["canonical_cr001"].__setitem__("validation_seed", 9170000))
    _audit_mutated(lambda c: c["canonical_cr001"].__setitem__("free_or_residual_defined_force_allowed", True))
    _audit_mutated(lambda c: c["canonical_cr001"].__setitem__("threshold_relaxation_under_same_experiment_allowed", True))


def test_independent_truth_states_remain_separate():
    _audit_mutated(lambda c: c["independent_truth_states"].__setitem__("visualization_ready", True))
    _audit_mutated(lambda c: c["independent_truth_states"].__setitem__("visual_correspondence_verified", True))
    _audit_mutated(lambda c: c["independent_truth_states"].__setitem__("pde_validated", True))
    _audit_mutated(lambda c: c["independent_truth_states"].__setitem__("paper_exact", True))
