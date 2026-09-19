from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_local_swirl_energy_root_governance import (
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


def test_repository_local_swirl_energy_root_governance_passes():
    receipt = audit_repository()
    assert receipt["status"] == "PASS"
    assert receipt["energy_root_found"] is True
    assert receipt["beta"] == pytest.approx(0.08837490297155456)
    assert receipt["frozen_material_path_aggregate_change"] == "reported_zero_at_receipt_precision"
    assert receipt["child_velocity_export_ready"] is False
    assert receipt["canonical_velocity_export_ready"] is True
    assert receipt["pde_validated"] is False
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["paper_exact"] is False


def test_energy_root_cannot_be_promoted_to_portable_or_openai_parameter():
    _audit_mutated(lambda c: c["root_semantics"].__setitem__("beta_is_portable_across_changed_parent_or_transform_identity", True))
    _audit_mutated(lambda c: c["root_semantics"].__setitem__("beta_is_a_production_coefficient", True))
    _audit_mutated(lambda c: c["root_semantics"].__setitem__("beta_is_an_openai_hidden_parameter", True))
    _audit_mutated(lambda c: c["root_semantics"].__setitem__("reusing_beta_without_rechecking_energy_on_changed_representation_allowed", True))


def test_order64_energy_root_cannot_replace_registered_cr001_energy_ladder():
    _audit_mutated(lambda c: c["formal_cr001_energy_distinction"].__setitem__("order_64_root_screen_is_registered_reference_energy_acceptance", True))
    _audit_mutated(lambda c: c["formal_cr001_energy_distinction"].__setitem__("registered_quadrature_orders_per_axis", [64]))


def test_theta_only_divergence_structure_cannot_transfer_parent_pde_receipts():
    _audit_mutated(lambda c: c["representation_scope"].__setitem__("theta_only_divergence_structure_is_navier_stokes_invariance", True))
    _audit_mutated(lambda c: c["representation_scope"].__setitem__("parent_pressure_receipt_may_transfer", True))
    _audit_mutated(lambda c: c["representation_scope"].__setitem__("parent_restricted_force_receipt_may_transfer", True))
    _audit_mutated(lambda c: c["representation_scope"].__setitem__("parent_momentum_receipt_may_transfer", True))


def test_central_identity_cannot_be_promoted_to_whole_field_identity():
    _audit_mutated(lambda c: c["representation_scope"].__setitem__("central_identity_implies_entire_velocity_field_identity", True))
    _audit_mutated(lambda c: c["representation_scope"].__setitem__("post_transform_common_scale_is_unity", False))


def test_zero_change_frozen_paths_cannot_be_generalized():
    _audit_mutated(lambda c: c["material_path_evidence"].__setitem__("zero_change_on_frozen_paths_proves_global_trajectory_identity", True))
    _audit_mutated(lambda c: c["material_path_evidence"].__setitem__("zero_change_on_frozen_paths_proves_velocity_field_identity", True))
    _audit_mutated(lambda c: c["material_path_evidence"].__setitem__("zero_change_on_frozen_paths_negates_eulerian_morphology_change", True))
    _audit_mutated(lambda c: c["material_path_evidence"].__setitem__("material_path_evidence_is_pde_evidence", True))


def test_target_free_morphology_cannot_promote_visual_or_pde_truth():
    _audit_mutated(lambda c: c["eulerian_capacity_evidence"].__setitem__("can_promote_visual_correspondence_verified", True))
    _audit_mutated(lambda c: c["eulerian_capacity_evidence"].__setitem__("can_promote_pde_validated", True))
    _audit_mutated(lambda c: c["eulerian_capacity_evidence"].__setitem__("public_openai_numeric_target_used", True))


def test_successful_root_and_paths_do_not_create_delivery_identity():
    _audit_mutated(lambda c: c["child_delivery_scope"].__setitem__("versioned_child_materialized", True))
    _audit_mutated(lambda c: c["child_delivery_scope"].__setitem__("save_load_contract_exists", True))
    _audit_mutated(lambda c: c["child_delivery_scope"].__setitem__("velocity_export_ready", True))
    _audit_mutated(lambda c: c["canonical_delivery"].__setitem__("local_swirl_energy_child_replaces_canonical_candidate", True))


def test_cr001_and_independent_truth_states_remain_fail_closed():
    _audit_mutated(lambda c: c["canonical_cr001"].__setitem__("pde_residual_max", 0.01))
    _audit_mutated(lambda c: c["canonical_cr001"].__setitem__("free_or_residual_defined_force_allowed", True))
    _audit_mutated(lambda c: c["independent_truth_states"].__setitem__("pde_validated", True))
    _audit_mutated(lambda c: c["independent_truth_states"].__setitem__("visual_correspondence_verified", True))
