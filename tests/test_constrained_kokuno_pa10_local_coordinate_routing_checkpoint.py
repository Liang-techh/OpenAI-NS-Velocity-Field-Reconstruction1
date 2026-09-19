import copy

import pytest

from openai_ns_reconstruction.kokuno_pa10_local_coordinate_routing_checkpoint import (
    AGENT4_ARTIFACT_DIGEST,
    AGENT4_ARTIFACT_ID,
    AGENT4_REAUDIT_HEAD,
    EXPECTED_FAILED_GUARDS,
    FORMAL_GATES,
    _canonical_sha256,
    build_checkpoint,
    validate_checkpoint,
)


@pytest.fixture(scope="module")
def checkpoint():
    value = build_checkpoint()
    validate_checkpoint(value)
    return value


def test_v46_admits_only_local_pressure_repair_and_keeps_eta_global_identity_closed(checkpoint):
    states = checkpoint["states"]
    assert states["selected_pressure_primitive_interface_exposed"] is True
    assert states["selected_pressure_stationary_local_value_path_repair_confirmed"] is True
    assert states["selected_pressure_local_p_eta_fd_guard_passed"] is True
    assert states["selected_pressure_old_global_phase_subtraction_blocker_closed"] is True
    assert states["selected_pressure_local_coordinate_contract_required"] is True

    assert states["selected_pressure_eta_fundamental_theorem_guard_passed"] is False
    assert states["selected_pressure_current_eta_api_independent_preflight_passed"] is False
    assert states["leading_ready"] is False
    assert states["same_cycle_requested_stress_materialized"] is False
    assert states["public_velocity_correction_materialized"] is False
    assert states["correction_ready"] is False
    assert states["candidate_artifact_instantiated"] is False
    assert states["velocity_export_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False

    assert states["oscillatory_ready"] is True
    assert states["correction_ingest_allowed"] is True
    assert states["compact_radial_stress_operator_independently_admitted"] is True
    assert states["independent_spacetime_radial_stress_generalization_admitted"] is True
    assert checkpoint["formal_gates_unchanged"] == FORMAL_GATES


def test_v46_pins_agent4_reaudit_identity_and_exact_remaining_failed_guards(checkpoint):
    a4 = checkpoint["upstream"]["agent4_repaired_pressure_reaudit"]
    assert a4["head"] == AGENT4_REAUDIT_HEAD
    assert a4["dedicated_workflow_run"] == 35438452570
    assert a4["artifact_id"] == AGENT4_ARTIFACT_ID
    assert a4["artifact_zip_digest"] == AGENT4_ARTIFACT_DIGEST
    assert a4["repaired_public_p_eta_independent_preflight_passed"] is False
    assert set(a4["failed_guards"]) == EXPECTED_FAILED_GUARDS


def test_v46_freezes_local_fd_repair_and_new_eta_coordinate_floor(checkpoint):
    a4 = checkpoint["upstream"]["agent4_repaired_pressure_reaudit"]
    levels = a4["p_eta_fd6_levels"]
    ratios = a4["p_eta_fd6_refinement_ratios"]
    fto = a4["eta_fundamental_theorem_levels"]

    rms = [row["relative_rms"] for row in levels]
    assert rms == pytest.approx(
        [1.7863971549e-1, 8.9919454005e-3, 3.8369244790e-4],
        rel=5e-7,
    )
    assert levels[-1]["relative_max"] == pytest.approx(4.5837622294e-4, rel=5e-7)
    assert ratios[0] > 19.0
    assert ratios[1] > 23.0
    assert a4["p_eta_public_reference"]["rms"] > 3.0e14
    assert a4["p_eta_sign_flip_relative_rms"] > 1.9

    fto_rms = [row["relative_rms"] for row in fto]
    assert fto_rms == pytest.approx(
        [1.04116586e-4, 7.22338148e-5, 5.48304239e-5],
        rel=5e-7,
    )
    assert fto[-1]["relative_max"] == pytest.approx(7.00226922e-5, rel=5e-7)

    diagnosis = checkpoint["pressure_coordinate_diagnosis"]
    assert diagnosis["coordinate_resolution_hypothesis_only"] is True
    assert diagnosis["formula_error_claimed"] is False
    assert diagnosis["observed_fto_floor_is_same_order_as_ulp_over_width"] is True
    assert diagnosis["ulp_over_layer_width"] == pytest.approx(4.3396e-5, rel=2e-4)
    assert diagnosis["finest_eta_fundamental_theorem_relative_rms"] > diagnosis[
        "eta_fundamental_theorem_rms_gate"
    ]


def test_v46_routes_typed_local_coordinate_chain_rule_without_claiming_pass(checkpoint):
    leading = checkpoint["typed_handoffs"]["leading"]
    assert leading["ready"] is False
    assert leading["current_selected_pressure_api"] == "evaluate(Y, eta)"
    assert leading["next_required_api"] == "evaluate_local(Y, s) or equivalent typed compensated coordinate"
    assert leading["local_coordinate_definition"] == "s=sqrt(Lambda)*(eta-eta_*)"
    assert leading["required_chain_rule"] == "p_eta=sqrt(Lambda)*p_s"


@pytest.mark.parametrize(
    "key",
    [
        "leading_ready",
        "selected_pressure_eta_fundamental_theorem_guard_passed",
        "selected_pressure_current_eta_api_independent_preflight_passed",
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "real_full_candidate_defect_consumed",
        "signed_mean_inverse_input_ready",
        "public_velocity_correction_materialized",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ],
)
def test_v46_rejects_overpromotion_even_if_resigned(checkpoint, key):
    tampered = copy.deepcopy(checkpoint)
    tampered["states"][key] = True
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError):
        validate_checkpoint(tampered)


def test_v46_rejects_dropping_local_coordinate_requirement_even_if_resigned(checkpoint):
    tampered = copy.deepcopy(checkpoint)
    tampered["states"]["selected_pressure_local_coordinate_contract_required"] = False
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="required admitted state changed"):
        validate_checkpoint(tampered)


def test_v46_rejects_agent4_artifact_substitution_even_if_resigned(checkpoint):
    tampered = copy.deepcopy(checkpoint)
    tampered["upstream"]["agent4_repaired_pressure_reaudit"]["artifact_id"] += 1
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="wrong Agent-4 repaired-pressure artifact"):
        validate_checkpoint(tampered)


def test_v46_rejects_failed_guard_laundering_even_if_resigned(checkpoint):
    tampered = copy.deepcopy(checkpoint)
    tampered["upstream"]["agent4_repaired_pressure_reaudit"]["failed_guards"] = []
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="failed-guard set changed"):
        validate_checkpoint(tampered)


def test_v46_rejects_threshold_change_even_if_resigned(checkpoint):
    tampered = copy.deepcopy(checkpoint)
    tampered["formal_gates_unchanged"]["held_out_normalized_momentum_max"] = 2.0e-3
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="formal final gates changed"):
        validate_checkpoint(tampered)
