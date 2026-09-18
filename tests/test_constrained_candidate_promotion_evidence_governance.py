from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_candidate_promotion_evidence_governance import (
    audit_candidate_promotion_evidence,
    load_constraints,
    load_contract,
    load_project_status,
)


def test_current_candidate_promotion_contract_passes():
    report = audit_candidate_promotion_evidence()
    assert report["status"] == "governance_pass"
    assert report["st006_seed"] == 9172801
    assert report["cr001_seed"] == 914027
    assert report["direct_st006_comparison_requires_replay"] is True
    assert report["parent_paired_claim_is_parent_local_only"] is True
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
    mutated["source_classification"]["st006_measured_residual_record"] = "public_source_fact"
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
        "allows_pde_validation_claim",
        "allows_visual_correspondence_claim",
        "allows_paper_exact_claim",
        "allows_openai_field_identification_claim",
    ):
        mutated = deepcopy(contract)
        mutated["paired_parent_claim_scope"][key] = True
        with pytest.raises(ValueError):
            audit_candidate_promotion_evidence(contract=mutated)


def test_rejects_protocol_identity_weakening():
    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["repository_baseline_promotion"]["direct_superiority_claim_requires_identical_protocol"] = False
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["repository_baseline_promotion"]["required_protocol_identity_fields"].remove("sample_points_and_seed")
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["repository_baseline_promotion"]["unknown_protocol_field_fails_closed"] = False
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    mutated = deepcopy(contract)
    mutated["canonical_cr001_comparability"]["directly_comparable_to_st006_without_replay"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)


def test_rejects_st006_record_or_cr001_protocol_drift():
    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["repository_baseline_promotion"]["retained_momentum_sampled_max"] = 0.09
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


def test_rejects_threshold_force_or_truth_state_laundering():
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
    mutated["scientific_truth_boundary"]["repository_baseline_improvement_implies_pde_validated"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(contract=mutated)

    status = load_project_status()
    mutated_status = deepcopy(status)
    mutated_status["states"]["paper_exact"] = True
    with pytest.raises(ValueError):
        audit_candidate_promotion_evidence(project_status=mutated_status)


def test_rejects_claim_that_governance_changed_the_science_object():
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
    ):
        mutated = deepcopy(contract)
        mutated["mutation_scope"][key] = True
        with pytest.raises(ValueError):
            audit_candidate_promotion_evidence(contract=mutated)
