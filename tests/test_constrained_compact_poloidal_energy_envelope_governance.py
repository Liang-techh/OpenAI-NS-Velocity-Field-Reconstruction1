from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_compact_poloidal_energy_envelope_governance import (
    audit_compact_poloidal_energy_envelope_scope,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_current_compact_poloidal_energy_envelope_scope_passes():
    report = audit_compact_poloidal_energy_envelope_scope()
    assert report["status"] == "governance_pass"
    assert report["mode_name"] == "COMPACT_C4_ODD_Z_POLOIDAL"
    assert report["upstream_pr"] == 309
    assert report["canonical_thresholds_unchanged"] is True
    assert report["diagnostic_energy_envelope_selected_as_bound"] is False
    assert report["normalization_root_selected_as_coefficient"] is False
    assert report["joint_renormalization_requires_new_candidate_validation"] is True
    assert report["small_morphology_span_is_capacity_evidence_only"] is True
    assert report["live_next_integration_task_preserved"] is True
    assert report["velocity_export_ready"] is True
    assert report["visual_correspondence_verified"] is False
    assert report["pde_validated"] is False


def test_rejects_source_laundering_and_evidence_class_promotion():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["source_classification"]["compact_poloidal_mode"] = "public_source_fact"
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["evidence_relation_labels_are_source_classes"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(scope=mutated)


def test_rejects_diagnostic_envelope_or_root_promotion():
    scope = load_scope()
    for key in (
        "diagnostic_envelope_is_selected_materialization_bound",
        "diagnostic_envelope_selects_coefficient",
        "diagnostic_envelope_is_complete_feasible_coefficient_interval",
        "diagnostic_exact_normalization_root_selects_coefficient",
        "fixed_parent_energy_admissibility_transfers_through_joint_renormalization",
        "validation_time_energy_range_is_full_candidate_acceptance",
        "quadrature_refinement_is_pde_validation",
        "previous_capacity_trial_invalidated_by_energy_envelope",
    ):
        mutated = deepcopy(scope)
        mutated["energy_envelope_evidence_only"][key] = True
        with pytest.raises(ValueError):
            audit_compact_poloidal_energy_envelope_scope(scope=mutated)


def test_rejects_morphology_to_pde_or_visual_promotion():
    scope = load_scope()
    for key in (
        "small_span_proves_no_pde_leverage",
        "zero_q90_q99_span_proves_mode_useless",
        "morphology_span_is_momentum_residual_jacobian",
        "morphology_span_is_visual_correspondence",
        "morphology_span_selects_coefficient_or_sign",
    ):
        mutated = deepcopy(scope)
        mutated["morphology_semantics"][key] = True
        with pytest.raises(ValueError):
            audit_compact_poloidal_energy_envelope_scope(scope=mutated)


def test_rejects_silent_routing_or_energy_gate_relaxation():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["routing_scope"]["energy_envelope_may_automatically_define_materialization_bound"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["routing_scope"]["exact_energy_root_may_automatically_define_materialized_child"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["routing_scope"]["reference_energy_gate_may_be_relaxed_for_larger_move"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["routing_scope"]["additional_basis_growth_allowed_before_materialization_screen"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(scope=mutated)


def test_rejects_missing_identity_fresh_validation_or_joint_renormalization_recheck():
    scope = load_scope()
    for key in (
        "explicit_autonomous_coefficient_bound_required",
        "explicit_coefficient_value_required",
        "new_representation_family_identity_required",
        "new_candidate_sha_required",
        "fresh_full_per_component_momentum_divergence_required",
        "joint_renormalization_requires_fresh_nonlinear_pde_validation",
        "model_selection_data_may_not_be_reused_as_acceptance_data",
    ):
        mutated = deepcopy(scope)
        mutated["future_nonzero_materialization"][key] = False
        with pytest.raises(ValueError):
            audit_compact_poloidal_energy_envelope_scope(scope=mutated)


def test_rejects_cr001_and_project_status_drift():
    constraints = load_constraints()
    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(constraints=mutated_constraints)

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["forcing"]["parameters"]["a"] = [0.0, 20.0]
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(constraints=mutated_constraints)

    status = load_project_status()
    mutated_status = deepcopy(status)
    mutated_status["next_integration_task"] = "add another compact poloidal basis"
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["integration_policy"]["avoid"].remove("coefficient_selection_from_capacity_or_morphology_alone")
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(project_status=mutated_status)

    mutated_status = deepcopy(status)
    mutated_status["states"]["pde_validated"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_energy_envelope_scope(project_status=mutated_status)
