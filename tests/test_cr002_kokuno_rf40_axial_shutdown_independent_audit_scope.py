from __future__ import annotations

import copy

from openai_ns_reconstruction.audit_kokuno_rf40_axial_shutdown_independent_audit_scope import (
    CANONICAL_CONSTRAINTS_GIT_BLOB_SHA1,
    EXPECTED_CR001,
    audit_default,
    audit_scope,
    load_canonical_constraints,
    load_contract,
)


def _mutated_contract():
    return copy.deepcopy(load_contract())


def test_default_scope_contract_is_fail_closed_and_replays_canonical_constraints():
    assert audit_default() == []
    canonical, blob_sha1 = load_canonical_constraints()
    assert blob_sha1 == CANONICAL_CONSTRAINTS_GIT_BLOB_SHA1
    assert canonical["nu"] == EXPECTED_CR001["nu"]
    assert canonical["validation"]["thresholds"]["pde_residual_max"] == EXPECTED_CR001["momentum_max"]
    assert canonical["validation"]["thresholds"]["divergence_max"] == EXPECTED_CR001["divergence_max"]


def test_four_way_provenance_is_explicit_and_nonempty():
    contract = load_contract()
    provenance = contract["four_way_provenance"]
    assert set(provenance) == {
        "user_requirements",
        "public_source_facts",
        "autonomous_repository_facts",
        "pending_or_unknown",
    }
    assert all(provenance[key] for key in provenance)
    assert contract["scope_logic_witness"]["classification"] == "autonomous_mechanics_only"
    assert contract["scope_logic_witness"]["not_a_public_source_fact"] is True
    assert contract["scope_logic_witness"]["not_candidate_numerical_evidence"] is True


def test_first_turn_audit_cannot_be_laundered_into_axial_shutdown_coverage():
    contract = _mutated_contract()
    contract["machine_locked_distinctions"]["agent4_995_covers_axial_shutdown_X1_to_X2"] = True
    errors = audit_scope(contract)
    assert any("agent4_995_covers_axial_shutdown_X1_to_X2" in error for error in errors)

    contract = _mutated_contract()
    contract["machine_locked_distinctions"]["first_turn_a4_995_audit_may_be_relabelled_as_axial_shutdown_audit"] = True
    errors = audit_scope(contract)
    assert any("first_turn_a4_995_audit_may_be_relabelled_as_axial_shutdown_audit" in error for error in errors)

    contract = _mutated_contract()
    contract["upstream_scope"]["agent4_995"]["axial_shutdown_probes_X1_to_X2"] = True
    errors = audit_scope(contract)
    assert any("axial_shutdown_probes_X1_to_X2" in error for error in errors)


def test_a1_production_checks_cannot_be_laundered_into_independent_cartesian_validation():
    contract = _mutated_contract()
    contract["machine_locked_distinctions"][
        "axial_shutdown_leading_independent_cartesian_public_velocity_audit_available"
    ] = True
    errors = audit_scope(contract)
    assert any("axial_shutdown_leading_independent_cartesian_public_velocity_audit_available" in error for error in errors)

    contract = _mutated_contract()
    contract["machine_locked_distinctions"][
        "source_coordinate_handoff_or_derivative_replay_may_be_relabelled_as_independent_axial_shutdown_cartesian_validation"
    ] = True
    errors = audit_scope(contract)
    assert any("source_coordinate_handoff_or_derivative_replay" in error for error in errors)

    contract = _mutated_contract()
    contract["upstream_scope"]["agent1_993"]["implementation_distinct_axial_shutdown_cartesian_public_velocity_audit"] = True
    errors = audit_scope(contract)
    assert any("implementation_distinct_axial_shutdown_cartesian_public_velocity_audit" in error for error in errors)


def test_non_consumed_sibling_identity_cannot_be_silently_promoted_to_composite():
    contract = _mutated_contract()
    contract["upstream_scope"]["agent5_996"]["agent1_993_recorded_as_non_consumed_sibling"] = False
    errors = audit_scope(contract)
    assert any("agent1_993_recorded_as_non_consumed_sibling" in error for error in errors)

    contract = _mutated_contract()
    contract["upstream_scope"]["agent5_996"][
        "leading_plus_oscillatory_velocity_through_axial_shutdown_materialized"
    ] = True
    errors = audit_scope(contract)
    assert any("leading_plus_oscillatory_velocity_through_axial_shutdown_materialized" in error for error in errors)

    contract = _mutated_contract()
    contract["machine_locked_distinctions"][
        "axial_shutdown_leading_materialization_implies_axial_shutdown_composite_materialization"
    ] = True
    errors = audit_scope(contract)
    assert any("axial_shutdown_leading_materialization_implies_axial_shutdown_composite_materialization" in error for error in errors)


def test_scientific_and_delivery_promotions_stay_fail_closed():
    contract = _mutated_contract()
    contract["kokuno_route_state"]["pde_validated"] = True
    errors = audit_scope(contract)
    assert any("kokuno_route_state.pde_validated" in error for error in errors)

    contract = _mutated_contract()
    contract["kokuno_route_state"]["velocity_export_ready"] = True
    errors = audit_scope(contract)
    assert any("kokuno_route_state.velocity_export_ready" in error for error in errors)

    contract = _mutated_contract()
    contract["canonical_delivery_independence"]["velocity_export_ready"] = False
    errors = audit_scope(contract)
    assert any("canonical_delivery_independence.velocity_export_ready" in error for error in errors)

    contract = _mutated_contract()
    contract["scientific_changes"]["scientific_thresholds_changed"] = True
    errors = audit_scope(contract)
    assert any("scientific_changes.scientific_thresholds_changed" in error for error in errors)


def test_cr001_snapshot_cannot_be_weakened():
    for key, bad in (
        ("momentum_max", 0.01),
        ("momentum_l2", 0.01),
        ("divergence_max", 1e-4),
        ("divergence_l2", 1e-4),
        ("held_out_points", 512),
        ("residual_defined_free_forcing_forbidden", False),
        ("candidate_collapse_forbidden", False),
        ("post_hoc_threshold_relaxation_forbidden", False),
    ):
        contract = _mutated_contract()
        contract["cr001_snapshot"][key] = bad
        errors = audit_scope(contract)
        assert "mismatch:cr001_snapshot" in errors


def test_canonical_constraints_replay_detects_contract_or_source_drift():
    contract = load_contract()
    canonical, blob_sha1 = load_canonical_constraints()
    assert audit_scope(contract, canonical, blob_sha1) == []

    mutated_canonical = copy.deepcopy(canonical)
    mutated_canonical["validation"]["thresholds"]["pde_residual_max"] = 0.01
    errors = audit_scope(contract, mutated_canonical, blob_sha1)
    assert "mismatch:canonical_constraints_replay" in errors

    errors = audit_scope(contract, canonical, "0" * 40)
    assert "mismatch:canonical_constraints_replay" in errors


def test_required_forbidden_promotions_cannot_be_deleted():
    contract = _mutated_contract()
    contract["forbidden_promotions"].remove("A4 #995 first-turn audit -> axial-shutdown audit")
    errors = audit_scope(contract)
    assert any(error.startswith("missing:forbidden_promotions:") for error in errors)


def test_future_promotion_contract_requires_actual_axial_shutdown_public_velocity_scope():
    contract = load_contract()
    requirements = contract["future_promotion_requirements"][
        "axial_shutdown_leading_independent_cartesian_public_velocity_audit_available"
    ]
    joined = " ".join(requirements)
    assert "serialized/reloaded A1 #993 public Cartesian velocity" in joined
    assert "X_1 < X <= X_2" in joined
    assert "both sides of X_1" in joined
    assert "implementation-distinct operator" in joined
    assert "whole-domain [24,48,96]" in joined


def test_scope_logic_is_set_coverage_not_source_or_candidate_numerics():
    contract = load_contract()
    witness = contract["scope_logic_witness"]
    assert "S_first={X <= X_1}" in witness["statement"]
    assert "S_shutdown={X_1 < X <= X_2}" in witness["statement"]
    assert witness["not_a_public_source_fact"] is True
    assert witness["not_candidate_numerical_evidence"] is True
