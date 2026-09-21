from __future__ import annotations

import copy

import openai_ns_reconstruction.kokuno_a5_current_partial_composite_mean_ingest_contract as mod


HEAD = "1" * 40


def _rehash(payload):
    out = copy.deepcopy(payload)
    out.pop("contract_sha256", None)
    out["contract_sha256"] = mod._sha256(out)
    return out


def test_contract_accepts_exact_registered_seam() -> None:
    payload = mod.build_contract(HEAD)
    assert mod.validate_contract(payload) == []
    assert payload["parent_a5"]["pr"] == 968
    assert payload["parent_a5"]["head"] == "5558d0859d142914ae8921a4286f56077fe2232c"
    assert payload["agent2_partial_composite"]["pr"] == 970
    assert payload["agent2_partial_composite"]["head"] == "3a6405bbd3d10b8c3b38078f0c989c45e5d407b4"
    assert payload["agent3_partial_nonlinear_mean"]["pr"] == 971
    assert payload["agent3_partial_nonlinear_mean"]["head"] == "bdb02ad487a175ab62748687c3172127d09aeb9c"
    assert payload["agent4_leading_only_audit"]["pr"] == 966


def test_narrow_truth_advance_is_partial_composite_and_mean_only() -> None:
    payload = mod.build_contract(HEAD)
    truth = payload["truth_boundary"]
    for key in (
        "current_cartesian_leading_velocity_materialized_through_xh",
        "current_partial_leading_plus_oscillatory_velocity_materialized_through_xh",
        "current_partial_mixed_nonlinear_mean_materialized",
        "current_partial_quadratic_nonlinear_mean_materialized",
        "current_partial_aggregate_nonlinear_mean_materialized",
    ):
        assert truth[key] is True

    for key in (
        "velocity_beyond_xh_materialized",
        "outer_global_leading_velocity_materialized",
        "global_compact_support_completed",
        "independent_a4_composite_audit_available",
        "independent_a4_nonlinear_mean_audit_available",
        "radial_inverse_performed_on_current_partial_mean",
        "scoped_nonlinear_mean_authorized_as_correction_target",
        "real_agent3_ns_correction_velocity_materialized",
        "matched_cartesian_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "complete_candidate_api_ready",
        "same_protocol_full_ns_residual_available",
        "same_protocol_st006_comparison_available_now",
        "current_candidate_eligible_for_full_ns_validation",
        "scientific_admission",
        "pde_validated" if "pde_validated" in truth else "paper_exact",
    ):
        if key in truth:
            assert truth[key] is False, key

    assert payload["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_a4_leading_audit_cannot_be_laundered_into_composite_audit() -> None:
    payload = mod.build_contract(HEAD)
    payload["agent4_leading_only_audit"]["audits_agent2_partial_composite"] = True
    payload = _rehash(payload)
    errors = mod.validate_contract(payload)
    assert "agent4_leading_only_audit_drift" in errors
    assert "a4_composite_scope_laundered" in errors


def test_a4_absence_cannot_be_promoted_in_truth_boundary() -> None:
    payload = mod.build_contract(HEAD)
    payload["truth_boundary"]["independent_a4_composite_audit_available"] = True
    payload = _rehash(payload)
    errors = mod.validate_contract(payload)
    assert "truth_boundary_drift" in errors
    assert "premature_promotion_independent_a4_composite_audit_available" in errors


def test_partial_mean_cannot_be_authorized_as_correction_target() -> None:
    payload = mod.build_contract(HEAD)
    payload["agent3_partial_nonlinear_mean"]["authorized_as_correction_target"] = True
    payload = _rehash(payload)
    errors = mod.validate_contract(payload)
    assert "agent3_partial_nonlinear_mean_drift" in errors
    assert "a3_correction_authorization_laundered" in errors


def test_partial_mean_cannot_be_relabelled_complete_ns_defect() -> None:
    payload = mod.build_contract(HEAD)
    payload["truth_boundary"]["complete_ns_defect_materialized"] = True
    payload = _rehash(payload)
    errors = mod.validate_contract(payload)
    assert "truth_boundary_drift" in errors
    assert "premature_promotion_complete_ns_defect_materialized" in errors


def test_pde_readiness_cannot_be_promoted() -> None:
    payload = mod.build_contract(HEAD)
    payload["readiness"]["pde_validated"] = True
    payload = _rehash(payload)
    errors = mod.validate_contract(payload)
    assert "readiness_drift" in errors
    assert "readiness_promoted" in errors


def test_fixed_final_gates_cannot_be_relaxed() -> None:
    payload = mod.build_contract(HEAD)
    payload["final_gate"]["normalized_momentum_volume_l2"] = 2.0e-3
    payload = _rehash(payload)
    errors = mod.validate_contract(payload)
    assert "final_gate_drift" in errors


def test_residual_defined_free_forcing_remains_forbidden() -> None:
    payload = mod.build_contract(HEAD)
    payload["frozen_science"]["residual_defined_free_forcing_forbidden"] = False
    payload = _rehash(payload)
    errors = mod.validate_contract(payload)
    assert "frozen_science_drift" in errors
    assert "free_forcing_firewall_removed" in errors


def test_a2_to_a3_exact_lineage_is_required() -> None:
    payload = mod.build_contract(HEAD)
    payload["agent3_partial_nonlinear_mean"]["consumed_agent2_composite_head"] = "2" * 40
    payload = _rehash(payload)
    errors = mod.validate_contract(payload)
    assert "agent3_partial_nonlinear_mean_drift" in errors
    assert "a2_a3_lineage_mismatch" in errors


def test_contract_sha_detects_unrehash_mutation() -> None:
    payload = mod.build_contract(HEAD)
    payload["st006_baseline"]["momentum_sampled_max"] = 0.0
    assert "contract_sha256_mismatch" in mod.validate_contract(payload)
