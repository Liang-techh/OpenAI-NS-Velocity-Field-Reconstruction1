from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_current_partial_radial_stress_ingest_contract import (
    AGENT2_LATEST_DIFFERENTIAL_SIBLING,
    AGENT3_RADIAL_STRESS,
    AGENT4_RADIAL_STRESS_AUDIT,
    FINAL_GATE,
    PARENT_A5,
    READINESS,
    ST006_BASELINE,
    build_contract,
    enforce_contract,
    validate_contract,
)

EXACT_HEAD = "a" * 40


def _payload():
    return build_contract(EXACT_HEAD)


def _rehash(payload):
    # Deliberately rebuild through the public constructor and then transplant the
    # mutation. Tests that exercise semantic drift do not need a valid digest;
    # digest failure is itself fail-closed behavior.
    return payload


def test_current_partial_radial_stress_contract_accepts_frozen_state():
    payload = _payload()
    assert validate_contract(payload) == []
    enforce_contract(payload)
    assert payload["parent_a5"] == PARENT_A5
    assert payload["agent3_radial_stress"] == AGENT3_RADIAL_STRESS
    assert payload["agent4_radial_stress_audit"] == AGENT4_RADIAL_STRESS_AUDIT
    assert payload["readiness"] == READINESS
    assert payload["final_gate"] == FINAL_GATE
    assert payload["st006_baseline"] == ST006_BASELINE


def test_exact_lineage_and_independent_reference_are_pinned():
    payload = _payload()
    a3 = payload["agent3_radial_stress"]
    a4 = payload["agent4_radial_stress_audit"]
    assert a3["head"] == "0a777668050ff9e73d6c8db1f72e5fabfa9b04a6"
    assert a3["source_blob"] == "9afa50c573186b86cd86fa7ffe44a5986a808b6d"
    assert a4["head"] == "2b015feccb5d55735ced8763fcb05a563a4cf12e"
    assert a4["audited_agent3_head"] == a3["head"]
    assert a4["uses_production_radial_inverse_as_reference"] is False
    assert a4["reference_path"] == (
        "piecewise cubic interpolation plus order-8 Gauss-Legendre cell integration"
    )


def test_latest_a2_differential_is_sibling_not_retroactive_dependency():
    payload = _payload()
    a2 = payload["agent2_latest_differential_sibling"]
    assert a2 == AGENT2_LATEST_DIFFERENTIAL_SIBLING
    assert a2["pr"] == 975
    assert a2["consumed_by_agent3_976"] is False
    assert a2["replaces_independent_agent4_validation"] is False


def test_truth_boundary_advances_only_scoped_radial_inverse():
    truth = _payload()["truth_boundary"]
    assert truth["current_partial_nonlinear_radial_inverse_performed"] is True
    assert truth["current_partial_quadratic_radial_stress_materialized"] is True
    assert truth["independent_a4_radial_stress_audit_registered"] is True
    assert truth["independent_a4_radial_stress_audit_admitted"] is False
    for key in (
        "velocity_beyond_xh_materialized",
        "outer_global_leading_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "scoped_current_partial_radial_stress_authorized_as_correction_target",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "complete_candidate_api_ready",
        "heldout_normalized_ns_residual_assessed",
        "same_protocol_st006_comparison_available_now",
        "residual_reduction_claimed",
        "scientific_admission",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False


def test_readiness_and_final_gates_are_unchanged():
    payload = _payload()
    assert payload["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert payload["final_gate"]["normalized_momentum_sampled_max"] == 1e-3
    assert payload["final_gate"]["normalized_momentum_volume_l2"] == 1e-3
    assert payload["final_gate"]["divergence_sampled_max"] == 1e-5
    assert payload["final_gate"]["divergence_volume_l2"] == 1e-5
    assert payload["final_gate"]["canonical_volume_quadrature_ladder"] == [24, 48, 96]
    assert payload["frozen_science"]["residual_defined_free_forcing_forbidden"] is True


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("truth_boundary", "complete_ns_defect_materialized"), True),
        (("truth_boundary", "scoped_current_partial_radial_stress_authorized_as_correction_target"), True),
        (("truth_boundary", "real_agent3_ns_correction_velocity_materialized"), True),
        (("truth_boundary", "scientific_admission"), True),
        (("readiness", "correction_ready"), True),
        (("readiness", "pde_validated"), True),
        (("final_gate", "normalized_momentum_sampled_max"), 2e-3),
        (("frozen_science", "residual_defined_free_forcing_forbidden"), False),
        (("agent4_radial_stress_audit", "uses_production_radial_inverse_as_reference"), True),
        (("agent4_radial_stress_audit", "exact_cartesian_axis_evidence"), True),
        (("agent2_latest_differential_sibling", "consumed_by_agent3_976"), True),
    ],
)
def test_contract_fails_closed_on_scientific_scope_mutations(path, value):
    payload = copy.deepcopy(_payload())
    payload[path[0]][path[1]] = value
    assert validate_contract(_rehash(payload))
    with pytest.raises(AssertionError):
        enforce_contract(payload)


def test_contract_rejects_invalid_or_drifted_exact_head():
    with pytest.raises(ValueError):
        build_contract("not-a-commit")
    payload = _payload()
    payload["agent5_exact_head"] = "NOT-HEX"
    assert "contract_sha256_mismatch" in validate_contract(payload)
    assert "agent5_exact_head_invalid" in validate_contract(payload)
