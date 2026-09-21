"""Regressions for CR002 RF40 power-law composite/audit scope governance."""
from __future__ import annotations

import copy
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_rf40_power_law_composite_audit_scope import (
    audit_scope,
    load_canonical_constraints,
    load_contract,
)

ROOT = Path(__file__).resolve().parents[1]
CANONICAL = ROOT / "configs" / "constraints.json"


def _baseline():
    return load_contract()


def _mutate(path, value):
    payload = copy.deepcopy(_baseline())
    cur = payload
    for key in path[:-1]:
        cur = cur[key]
    cur[path[-1]] = value
    return payload


def _has(errors, token):
    return any(token in err for err in errors)


def test_contract_baseline_passes_without_repository_replay():
    assert audit_scope(_baseline()) == []


@pytest.mark.skipif(not CANONICAL.exists(), reason="repository canonical constraints not present")
def test_contract_matches_live_canonical_constraints():
    canonical, blob = load_canonical_constraints(CANONICAL)
    assert audit_scope(_baseline(), canonical, blob) == []


def test_base_head_drift_fails_closed():
    errors = audit_scope(_mutate(("exact_base", "head"), "0" * 40))
    assert _has(errors, "exact_base")


def test_active_integration_head_drift_fails_closed():
    errors = audit_scope(_mutate(("active_integration_reference", "head"), "1" * 40))
    assert _has(errors, "active_integration_reference")


def test_four_way_provenance_bucket_removal_fails_closed():
    payload = _baseline()
    del payload["four_way_provenance"]["pending_or_unknown"]
    errors = audit_scope(payload)
    assert _has(errors, "four_way_provenance")


def test_agent1_identity_drift_fails_closed():
    errors = audit_scope(_mutate(("upstream_scope", "agent1_1005", "source_blob"), "2" * 40))
    assert _has(errors, "agent1_1005.source_blob")


def test_agent1_leading_cannot_be_promoted_to_power_law_composite():
    errors = audit_scope(_mutate(("upstream_scope", "agent1_1005", "power_law_composite_materialized"), True))
    assert _has(errors, "power_law_composite_materialized")


def test_lambda_turn_a2_cannot_be_relabelled_as_power_law_composite():
    errors = audit_scope(_mutate(("upstream_scope", "agent2_1004", "power_law_composite_materialized"), True))
    assert _has(errors, "agent2_1004.power_law_composite_materialized")


def test_lambda_turn_a4_cannot_claim_strict_power_law_probes():
    errors = audit_scope(_mutate(("upstream_scope", "agent4_1007", "strict_power_law_probes_X3_to_X4"), True))
    assert _has(errors, "strict_power_law_probes_X3_to_X4")


def test_a5_cannot_register_unmaterialized_power_law_composite():
    errors = audit_scope(_mutate(("upstream_scope", "agent5_1008", "leading_plus_oscillatory_velocity_through_power_law_materialized"), True))
    assert _has(errors, "leading_plus_oscillatory_velocity_through_power_law_materialized")


def test_machine_scope_cannot_promote_power_law_composite():
    errors = audit_scope(_mutate(("machine_locked_distinctions", "power_law_composite_materialized"), True))
    assert _has(errors, "power_law_composite_materialized")


def test_machine_scope_cannot_promote_independent_power_law_audit():
    errors = audit_scope(_mutate(("machine_locked_distinctions", "power_law_composite_independent_cartesian_public_velocity_audit_available"), True))
    assert _has(errors, "power_law_composite_independent_cartesian_public_velocity_audit_available")


def test_a1_self_checks_cannot_be_relabelled_as_independent_composite_validation():
    errors = audit_scope(_mutate(("machine_locked_distinctions", "agent1_1005_production_or_source_coordinate_checks_may_be_relabelled_as_independent_power_law_composite_validation"), True))
    assert _has(errors, "agent1_1005_production_or_source_coordinate_checks_may_be_relabelled")


def test_lambda_turn_evidence_cannot_be_relabelled_as_power_law_evidence():
    errors = audit_scope(_mutate(("machine_locked_distinctions", "lambda_turn_composite_or_audit_may_be_relabelled_as_power_law_composite_or_audit"), True))
    assert _has(errors, "lambda_turn_composite_or_audit_may_be_relabelled")


def test_cr001_threshold_mutation_fails_closed():
    errors = audit_scope(_mutate(("cr001_snapshot", "momentum_max"), 0.002))
    assert _has(errors, "cr001_snapshot")


def test_free_residual_forcing_enablement_fails_closed():
    errors = audit_scope(_mutate(("cr001_snapshot", "residual_defined_free_forcing_forbidden"), False))
    assert _has(errors, "cr001_snapshot")


def test_canonical_eq45_delivery_cannot_be_downgraded_by_kokuno_gap():
    errors = audit_scope(_mutate(("canonical_delivery_independence", "velocity_export_ready"), False))
    assert _has(errors, "canonical_delivery_independence.velocity_export_ready")


def test_kokuno_export_cannot_be_promoted_by_power_law_leading_only():
    errors = audit_scope(_mutate(("kokuno_route_state", "velocity_export_ready"), True))
    assert _has(errors, "kokuno_route_state.velocity_export_ready")


def test_scientific_byte_change_claim_fails_closed():
    errors = audit_scope(_mutate(("scientific_changes", "velocity_or_profile_bytes_changed"), True))
    assert _has(errors, "scientific_changes.velocity_or_profile_bytes_changed")


def test_required_forbidden_promotion_cannot_be_removed():
    payload = _baseline()
    payload["forbidden_promotions"].remove(
        "A2 #1004 lambda-turn composite -> power-law composite"
    )
    errors = audit_scope(payload)
    assert _has(errors, "forbidden_promotions")


def test_future_promotion_requirements_must_stay_nonempty():
    payload = _baseline()
    payload["future_promotion_requirements"]["power_law_composite_ready"] = []
    errors = audit_scope(payload)
    assert _has(errors, "future_promotion_requirements.power_law_composite_ready")


def test_scope_witness_is_mechanics_only_not_candidate_evidence():
    payload = _baseline()
    witness = payload["scope_logic_witness"]
    assert witness["classification"] == "autonomous_mechanics_only"
    assert witness["not_a_public_source_fact"] is True
    assert witness["not_candidate_numerical_evidence"] is True
    inner = {0, 1, 2, 3}
    strict_outer = {4, 5}
    assert inner.isdisjoint(strict_outer)


def test_live_canonical_mutation_is_detected_when_supplied():
    payload = _baseline()
    canonical = {
        "nu": 0.02,
        "domain": {
            "physical": "R^3",
            "evaluation_box": [[-2, 2], [-2, 2], [-2, 2]],
            "support": "r < 2 and abs(z) < 2",
            "time_interval": [0.25, 0.75],
        },
        "forcing": {
            "mode": "restricted_two_parameter_family",
            "restriction": "Only a,c may be fitted. No residual-dependent basis or pointwise free force.",
        },
        "nontriviality": {
            "reference_energy": 1.0,
            "reference_energy_abs_tolerance": 0.001,
            "enforcement": "normalize initial velocity parameters, reject collapsed candidates",
        },
        "validation": {
            "seed": 914027,
            "held_out_points": 4096,
            "derivative_steps": [0.02, 0.01, 0.005],
            "quadrature_orders_per_axis": [24, 48, 96],
            "thresholds": {
                "pde_residual_max": 0.001,
                "pde_residual_L2": 0.001,
                "divergence_max": 0.00001,
                "divergence_L2": 0.00001,
            },
            "failure_policy": "retain failed results; changing thresholds requires a new experiment version",
        },
    }
    errors = audit_scope(payload, canonical, "6c559e42895a606e2ef025ade4cb448966d75814")
    assert _has(errors, "live_canonical_constraints")
