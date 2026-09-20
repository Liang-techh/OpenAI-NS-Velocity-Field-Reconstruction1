from __future__ import annotations

from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_kokuno_cartesian_center_definition_evaluator_identity_scope import (
    GovernanceError,
    audit,
    evaluator_identity_witness,
    load_constraints,
    load_contract,
)


def test_exact_upstream_definition_evaluator_identity_audit_passes():
    result = audit()
    assert result["upstream_pr"] == 803
    assert result["upstream_head"] == "0666fff748d7c2659774b0aacba70042d045aab5"
    assert result["CR001_unchanged"] is True

    structure = result["source_structure"]
    assert structure["default_bisection_iterations"] == 80
    assert structure["q_inverse_uses_configured_iteration_count"] is True
    assert structure["field_configuration_binds_q_inverse_iterations"] is True
    assert structure["field_sha256_hashes_complete_field_configuration"] is True
    assert structure["save_load_preserves_evaluator_realization"] is True

    witness = result["evaluator_identity_witness"]
    assert witness["same_physical_profile_configuration"] is True
    assert witness["same_registered_time_interval"] is True
    assert witness["same_source_q_inverse_equation"] is True
    assert witness["same_q_inverse_method"] is True
    assert witness["different_evaluator_iteration_count"] is True
    assert witness["different_current_field_sha256"] is True
    assert witness["velocity_bitwise_equality_claimed"] is False
    assert witness["public_source_definition_difference_claimed"] is False

    state = result["identity_state"]
    assert state["inner_velocity_callable"] is True
    assert state["current_field_sha_is_exact_evaluator_replay_identity"] is True
    assert state["source_mathematical_definition_digest_exposed"] is False
    assert state["solver_iterations_are_public_source_data"] is False
    assert state["unified_kokuno_velocity_export_ready_promoted"] is False
    assert state["visual_correspondence_verified"] is False
    assert state["pde_validated"] is False
    assert state["paper_exact"] is False
    assert state["openai_field_identified"] is False


def test_behavioral_witness_changes_digest_without_claiming_velocity_equivalence():
    witness = evaluator_identity_witness()
    assert witness["left_iterations"] == 80
    assert witness["right_iterations"] == 96
    assert witness["different_current_field_sha256"] is True
    assert witness["velocity_bitwise_equality_claimed"] is False


def test_public_source_laundering_of_iteration_count_fails_closed():
    contract = load_contract()
    contract["source_classification"]["public_source_fact"].append(
        "DEFAULT_BISECTION_ITERATIONS=80 is a public-source parameter"
    )
    with pytest.raises(GovernanceError, match="laundered into public-source facts"):
        audit(contract=contract)


def test_field_sha_cannot_be_promoted_to_pure_source_definition_identity():
    contract = load_contract()
    contract["identity_roles"][
        "current_field_sha_is_source_mathematical_definition_identity"
    ] = True
    with pytest.raises(GovernanceError, match="identity laundering"):
        audit(contract=contract)


def test_cross_implementation_validation_may_not_require_same_evaluator_digest():
    contract = load_contract()
    contract["cross_implementation_validation_rules"][
        "independent_validator_may_use_different_q_solver_realization"
    ] = False
    with pytest.raises(GovernanceError, match="cross-implementation rule lost"):
        audit(contract=contract)


def test_exact_executable_artifact_may_not_ignore_evaluator_identity():
    contract = load_contract()
    contract["cross_implementation_validation_rules"][
        "different_evaluator_sha_may_be_called_same_exact_executable_artifact"
    ] = True
    with pytest.raises(GovernanceError, match="cross-implementation laundering"):
        audit(contract=contract)


def test_unified_export_and_scientific_states_cannot_be_promoted_by_inner_identity_audit():
    contract = load_contract()
    contract["delivery_claim_boundary"]["unified_kokuno_velocity_export_ready"] = True
    with pytest.raises(GovernanceError, match="premature delivery/scientific promotion"):
        audit(contract=contract)

    contract = load_contract()
    contract["delivery_claim_boundary"]["pde_validated"] = True
    with pytest.raises(GovernanceError, match="premature delivery/scientific promotion"):
        audit(contract=contract)


def test_cr001_threshold_relaxation_and_free_force_routes_fail_closed():
    constraints = load_constraints()
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.002
    with pytest.raises(GovernanceError, match="threshold drift"):
        audit(constraints=constraints)

    constraints = load_constraints()
    constraints["forcing"]["restriction"] = "Only a,c may be fitted."
    with pytest.raises(GovernanceError, match="residual-defined force allowed"):
        audit(constraints=constraints)


def test_source_mutation_dropping_iteration_binding_fails_closed():
    source_path = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "openai_ns_reconstruction"
        / "kokuno_pa10_cartesian_center_velocity.py"
    )
    source = source_path.read_text()
    mutated = source.replace(
        '"iterations": int(self.bisection_iterations)',
        '"iterations": 80',
        1,
    )
    assert mutated != source
    with pytest.raises(
        GovernanceError, match="field_configuration identity binding missing"
    ):
        audit(source_text=mutated)


def test_claim_boundary_mutation_fails_closed():
    contract = load_contract()
    contract["claim_boundaries"][
        "field_sha_difference_proves_different_public_source_formula"
    ] = True
    with pytest.raises(GovernanceError, match="claim boundary promoted"):
        audit(contract=contract)
