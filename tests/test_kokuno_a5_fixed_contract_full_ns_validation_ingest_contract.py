from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.kokuno_a5_fixed_contract_full_ns_validation_ingest_contract import (
    AGENT3,
    AGENT4,
    FINAL_GATE,
    PARENT_A5,
    QUADRATURE_LADDER,
    ST006_BASELINE,
    _contract_payload_without_sha,
    _sha256,
    build_contract,
    validate_contract,
)


def _fresh() -> dict:
    return build_contract(exact_head="a" * 40)


def _resign(contract: dict) -> dict:
    contract["contract_sha256"] = _sha256(_contract_payload_without_sha(contract))
    return contract


def _must_reject(mutator) -> None:
    contract = copy.deepcopy(_fresh())
    mutator(contract)
    _resign(contract)
    with pytest.raises(ValueError):
        validate_contract(contract)


def test_contract_builds_and_hash_round_trips() -> None:
    contract = _fresh()
    validate_contract(contract)
    assert contract["contract_sha256"] == _sha256(_contract_payload_without_sha(contract))


def test_exact_provenance_and_unresolved_ci_are_snapshot_bound() -> None:
    contract = _fresh()
    assert contract["parent_a5"] == PARENT_A5
    a3 = contract["agent3_fixed_physical_contract_handoff"]
    a4 = contract["agent4_blackbox_full_ns_validator_handoff"]
    assert a3["pr"] == 895
    assert a3["head"] == AGENT3["head"]
    assert a3["source_blob"] == "0cd41dfcc998064c37075b6086a246e93831ebfd"
    assert a3["observed_ci"]["dedicated"]["status"] == "queued"
    assert a3["observed_ci"]["dedicated"]["conclusion"] is None
    assert a4["pr"] == 896
    assert a4["head"] == AGENT4["head"]
    assert a4["source_blob"] == "a9e46daf60d909a63261dd10821e1a08101dd872"
    assert a4["observed_ci"]["dedicated"]["status"] == "queued"
    assert a4["observed_ci"]["dedicated"]["conclusion"] is None


def test_fixed_physical_contract_and_velocity_only_correction_firewall() -> None:
    a3 = _fresh()["agent3_fixed_physical_contract_handoff"]
    assert a3["api"] == "run_fixed_physical_contract_finite_correction_stage"
    assert a3["physical_contract_sha256_required"] is True
    assert a3["before_after_same_physical_contract_sha256"] is True
    assert a3["velocity_only_correction"] is True
    assert a3["pressure_frozen_through_correction"] is True
    assert a3["forcing_frozen_through_correction"] is True
    assert a3["held_in_held_out_disjoint_required"] is True
    assert a3["residual_defined_forcing_forbidden"] is True
    assert a3["joint_pressure_or_forcing_updates_require_separate_preregistered_stage"] is True
    assert a3["scientific_admission"] is False


def test_blackbox_validator_and_canonical_l2_firewall() -> None:
    a4 = _fresh()["agent4_blackbox_full_ns_validator_handoff"]
    assert a4["public_candidate_surfaces"] == [
        "velocity",
        "pressure",
        "forcing",
        "validation_metadata",
    ]
    assert a4["construction_derivatives_allowed"] is False
    assert a4["construction_defect_tensors_allowed"] is False
    assert a4["training_state_allowed"] is False
    assert a4["operator"]["viscosity"] == 0.01
    assert a4["operator"]["derivative_step_ladder"] == [0.02, 0.01, 0.005]
    assert a4["heldout"]["seed"] == 914027
    assert a4["heldout"]["sample_count"] == 4096
    assert a4["canonical_quadrature_ladder"] == QUADRATURE_LADDER == [24, 48, 96]
    assert a4["quadrature_ladder_assessed_required_for_pde_validated"] is True
    assert a4["scientific_admission"] is False


def test_current_candidate_remains_ineligible_and_readiness_unchanged() -> None:
    contract = _fresh()
    seam = contract["integration_seam"]
    for key in (
        "global_join_present",
        "matched_pressure_present",
        "restricted_forcing_present",
        "restricted_forcing_preregistered",
        "real_correction_velocity_present",
        "complete_candidate_api_ready",
        "current_candidate_eligible_for_full_ns_validation",
        "same_protocol_st006_comparison_available_now",
    ):
        assert seam[key] is False
    assert contract["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }


def test_final_gates_and_st006_baseline_are_frozen() -> None:
    contract = _fresh()
    assert contract["final_gate"] == FINAL_GATE == {
        "momentum_sampled_max": 1.0e-3,
        "momentum_volume_l2": 1.0e-3,
        "divergence_sampled_max": 1.0e-5,
        "divergence_volume_l2": 1.0e-5,
    }
    assert contract["baseline"]["st006"] == ST006_BASELINE
    assert ST006_BASELINE["momentum_sampled_max"] == 0.1082289305112118
    assert ST006_BASELINE["momentum_volume_l2"] == 0.10758432876230622


@pytest.mark.parametrize(
    "mutator",
    [
        lambda c: c["parent_a5"].__setitem__("head", "0" * 40),
        lambda c: c["agent3_fixed_physical_contract_handoff"].__setitem__("head", "0" * 40),
        lambda c: c["agent4_blackbox_full_ns_validator_handoff"].__setitem__("head", "0" * 40),
        lambda c: c["agent3_fixed_physical_contract_handoff"].__setitem__("scientific_admission", True),
        lambda c: c["agent4_blackbox_full_ns_validator_handoff"].__setitem__("scientific_admission", True),
        lambda c: c["agent3_fixed_physical_contract_handoff"].__setitem__("velocity_only_correction", False),
        lambda c: c["agent3_fixed_physical_contract_handoff"].__setitem__("pressure_frozen_through_correction", False),
        lambda c: c["agent3_fixed_physical_contract_handoff"].__setitem__("forcing_frozen_through_correction", False),
        lambda c: c["agent3_fixed_physical_contract_handoff"].__setitem__("residual_defined_forcing_forbidden", False),
        lambda c: c["agent4_blackbox_full_ns_validator_handoff"].__setitem__("construction_derivatives_allowed", True),
        lambda c: c["agent4_blackbox_full_ns_validator_handoff"].__setitem__("construction_defect_tensors_allowed", True),
        lambda c: c["agent4_blackbox_full_ns_validator_handoff"].__setitem__("canonical_quadrature_ladder", [24, 48]),
        lambda c: c["agent4_blackbox_full_ns_validator_handoff"].__setitem__("quadrature_ladder_assessed_required_for_pde_validated", False),
        lambda c: c["integration_seam"].__setitem__("correction_and_validator_physical_contract_sha256_must_match", False),
        lambda c: c["integration_seam"].__setitem__("global_join_present", True),
        lambda c: c["integration_seam"].__setitem__("matched_pressure_present", True),
        lambda c: c["integration_seam"].__setitem__("restricted_forcing_present", True),
        lambda c: c["integration_seam"].__setitem__("real_correction_velocity_present", True),
        lambda c: c["integration_seam"].__setitem__("complete_candidate_api_ready", True),
        lambda c: c["integration_seam"].__setitem__("current_candidate_eligible_for_full_ns_validation", True),
        lambda c: c["integration_seam"].__setitem__("same_protocol_st006_comparison_available_now", True),
        lambda c: c["readiness"].__setitem__("pde_validated", True),
        lambda c: c["readiness"].__setitem__("correction_ready", True),
        lambda c: c["final_gate"].__setitem__("momentum_volume_l2", 1.001e-3),
        lambda c: c["final_gate"].__setitem__("divergence_volume_l2", 1.001e-5),
        lambda c: c["baseline"]["st006"].__setitem__("momentum_volume_l2", 0.01),
        lambda c: c["truth_boundary"].__setitem__("sampled_rms_is_not_canonical_volume_l2", False),
        lambda c: c["truth_boundary"].__setitem__("monte_carlo_volume_l2_does_not_replace_canonical_quadrature_ladder", False),
        lambda c: c["truth_boundary"].__setitem__("residual_defined_forcing_forbidden", False),
        lambda c: c["truth_boundary"].__setitem__("kokuno_reconstruction_is_paper_exact", True),
    ],
)
def test_fail_closed_semantic_mutations(mutator) -> None:
    _must_reject(mutator)


def test_hash_tamper_is_rejected_without_resigning() -> None:
    contract = _fresh()
    contract["contract_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="contract SHA mismatch"):
        validate_contract(contract)
