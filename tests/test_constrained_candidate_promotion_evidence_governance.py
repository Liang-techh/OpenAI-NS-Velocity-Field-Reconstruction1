from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_candidate_promotion_evidence_governance import (
    audit_candidate_promotion_evidence,
    load_constraints,
    load_contract,
    load_project_status,
)


def test_current_candidate_promotion_lifecycle_contract_passes():
    report = audit_candidate_promotion_evidence()
    assert report["status"] == "governance_pass"
    assert report["retained_baseline"] == "ST006"
    assert report["retained_baseline_pde_validated"] is False
    assert report["current_challengers_are_unintegrated"] is True
    assert report["direct_st006_comparison_requires_replay"] is True
    assert report["parent_paired_claim_is_parent_local_only"] is True
    assert report["canonical_delivery_is_separate_from_pde_benchmark"] is True
    assert report["canonical_thresholds_unchanged"] is True
    assert report["velocity_export_ready"] is True
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False


def test_rejects_fifth_source_class_or_measurement_laundering():
    contract = load_contract()

    mutated = deepcopy(contract)
    mutated["allowed_source_classes"].append("project_derived_measurement")
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["source_classification"]["challenger_branch_measurements"] = "public_source_fact"
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["relation_labels_are_source_classes"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)


def test_rejects_parent_local_improvement_as_repository_or_scientific_promotion():
    contract = load_contract()
    for key in (
        "allows_repository_st006_superiority_claim",
        "allows_retained_baseline_replacement",
        "allows_canonical_velocity_replacement",
        "allows_pde_validation_claim",
        "allows_visual_correspondence_claim",
        "allows_paper_exact_claim",
        "allows_openai_field_identification_claim",
    ):
        mutated = deepcopy(contract)
        mutated["paired_parent_claim_scope"][key] = True
        with pytest.raises(ValueError):
            audit_candidate_promotion_evidence(contract=mutated)


def test_rejects_lifecycle_laundering_from_open_challenger():
    contract = load_contract()

    for key in (
        "open_or_draft_pr_changes_retained_baseline",
        "protocol_scoped_superiority_changes_retained_baseline",
        "protocol_scoped_superiority_changes_canonical_velocity",
        "retained_baseline_replacement_implies_pde_validated",
    ):
        mutated = deepcopy(contract)
        mutated["lifecycle_separation"][key] = True
        with pytest.raises(ValueError):
            audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["current_challenger_snapshot"]["st048"]["retained_baseline_replaced"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["current_challenger_snapshot"]["st047_e"]["pde_validated"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)


def test_rejects_promotion_without_explicit_identity_receipt_or_status_update():
    contract = load_contract()
    for key in (
        "retained_baseline_replacement_requires_explicit_integration_or_promotion",
        "retained_baseline_replacement_requires_immutable_candidate_identity",
        "retained_baseline_replacement_requires_governed_validation_receipt",
        "retained_baseline_replacement_requires_status_document_update",
        "canonical_velocity_replacement_requires_explicit_api_or_candidate_status_update",
        "scientific_acceptance_requires_all_unchanged_formal_gates",
    ):
        mutated = deepcopy(contract)
        mutated["lifecycle_separation"][key] = False
        with pytest.raises(ValueError):
            audit_candidate_promotion_evidence(contract=mutated)


def test_rejects_protocol_identity_weakening():
    contract = load_contract()

    mutated = deepcopy(contract)
    mutated["repository_baseline"]["direct_superiority_claim_requires_identical_protocol"] = False
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["repository_baseline"]["required_protocol_identity_fields"].remove("sample_points_and_seed")
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["repository_baseline"]["unknown_protocol_field_fails_closed"] = False
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["canonical_cr001_comparability"]["directly_comparable_to_st006_without_replay"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)


def test_rejects_st006_record_or_cr001_protocol_drift():
    contract = load_contract()

    mutated = deepcopy(contract)
    mutated["repository_baseline"]["retained_momentum_sampled_max"] = 0.09
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    constraints = load_constraints()
    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["seed"] = 9172801
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(constraints=mutated_constraints)

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["held_out_points"] = 2048
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(constraints=mutated_constraints)


def test_rejects_threshold_force_delivery_or_truth_state_laundering():
    constraints = load_constraints()
    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(constraints=mutated_constraints)

    mutated_constraints = deepcopy(constraints)
    mutated_constraints["forcing"]["parameters"]["c"] = [0.0, 20.0]
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(constraints=mutated_constraints)

    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["delivery_separation"]["benchmark_challenge_changes_velocity_api_automatically"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["scientific_truth_boundary"]["retained_baseline_replacement_implies_pde_validated"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    status = load_project_status()
    mutated_status = deepcopy(status)
    mutated_status["states"]["paper_exact"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(project_status=mutated_status)


def test_rejects_claim_that_governance_changed_the_science_object_or_route():
    contract = load_contract()
    for key in (
        "velocity_changed",
        "candidate_changed",
        "pressure_changed",
        "forcing_changed",
        "optimizer_changed",
        "sampling_changed",
        "norm_definition_changed",
        "threshold_changed",
        "active_route_changed",
    ):
        mutated = deepcopy(contract)
        mutated["mutation_scope"][key] = True
        with pytest.raises(ValueError):
            audit_candidate_promotion_evidence(contract=mutated)
