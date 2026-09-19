from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_localized_piola_energy_root_governance import (
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


def test_repository_energy_root_governance_passes():
    receipt = audit_repository()
    assert receipt["status"] == "PASS"
    assert receipt["energy_root_found"] is False
    assert receipt["scientific_result"] == "negative_capacity_evidence_for_frozen_family_only"
    assert receipt["canonical_velocity_export_ready"] is True
    assert receipt["pde_validated"] is False
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["paper_exact"] is False


def test_green_execution_cannot_become_scientific_root_pass():
    _audit_mutated(lambda c: c["scientific_failure_scope"].__setitem__("scientific_energy_neutral_crossing_found", True))
    _audit_mutated(lambda c: c["scientific_failure_scope"].__setitem__("green_ci_is_scientific_pass", True))
    _audit_mutated(lambda c: c["energy_root_experiment"].__setitem__("root_exists_in_preregistered_bracket", True))


def test_no_root_cannot_be_promoted_to_global_impossibility():
    _audit_mutated(lambda c: c["scientific_failure_scope"].__setitem__("no_root_proves_no_energy_neutral_localized_piola_transform_exists", True))
    _audit_mutated(lambda c: c["scientific_failure_scope"].__setitem__("no_root_proves_no_other_divergence_preserving_compensation_can_work", True))
    _audit_mutated(lambda c: c["scientific_failure_scope"].__setitem__("no_root_identifies_any_openai_hidden_parameter", True))


def test_scaled_taper_evidence_cannot_transfer_to_nonexistent_root_child():
    _audit_mutated(lambda c: c["scaled_taper_evidence"].__setitem__("material_path_evidence_transfers_to_hypothetical_unscaled_energy_root_child", True))
    _audit_mutated(lambda c: c["scaled_taper_evidence"].__setitem__("eulerian_morphology_evidence_transfers_to_hypothetical_unscaled_energy_root_child", True))


def test_order64_root_screen_cannot_replace_registered_cr001_energy_ladder():
    _audit_mutated(lambda c: c["formal_cr001_energy_distinction"].__setitem__("order_64_root_screen_is_registered_reference_energy_acceptance", True))
    _audit_mutated(lambda c: c["formal_cr001_energy_distinction"].__setitem__("registered_quadrature_orders_per_axis", [64]))


def test_failed_root_experiment_cannot_materialize_or_replace_canonical_delivery():
    _audit_mutated(lambda c: c["energy_root_experiment"].__setitem__("versioned_child_materialized", True))
    _audit_mutated(lambda c: c["energy_root_experiment"].__setitem__("velocity_export_ready", True))
    _audit_mutated(lambda c: c["canonical_delivery"].__setitem__("energy_root_experiment_replaces_canonical_candidate", True))


def test_piola_or_render_evidence_cannot_promote_pde_or_visual_truth():
    _audit_mutated(lambda c: c["representation_and_pde_scope"].__setitem__("piola_divergence_structure_is_navier_stokes_invariance", True))
    _audit_mutated(lambda c: c["representation_and_pde_scope"].__setitem__("target_free_morphology_or_render_can_promote_pde_validated", True))
    _audit_mutated(lambda c: c["representation_and_pde_scope"].__setitem__("target_free_morphology_or_render_can_promote_visual_correspondence_verified", True))


def test_post_result_family_change_requires_new_experiment_identity():
    _audit_mutated(lambda c: c["scientific_failure_scope"].__setitem__("changing_sign_window_bracket_or_compensation_family_after_result_requires_new_experiment_identity", False))
    _audit_mutated(lambda c: c["scientific_failure_scope"].__setitem__("same_experiment_may_relax_or_expand_bracket_after_observing_failure", True))
