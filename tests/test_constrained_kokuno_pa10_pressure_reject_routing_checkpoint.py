import copy

import pytest

from openai_ns_reconstruction.kokuno_pa10_pressure_reject_routing_checkpoint import (
    AGENT3_SPACETIME_ADMISSION_HEAD,
    AGENT4_PRESSURE_ARTIFACT_ID,
    AGENT4_PRESSURE_AUDIT_HEAD,
    EXPECTED_FAILED_GUARDS,
    FORMAL_GATES,
    _canonical_sha256,
    build_checkpoint,
    validate_checkpoint,
)


@pytest.fixture(scope="module")
def checkpoint():
    """Run the expensive independent Agent-4 audit once for this regression module."""
    value = build_checkpoint()
    validate_checkpoint(value)
    return value


def test_v45_routes_scientific_pressure_reject_without_weakening_prior_handoffs(checkpoint):
    states = checkpoint["states"]

    assert states["selected_pressure_primitive_interface_exposed"] is True
    assert states["selected_pressure_primitive_independent_preflight_passed"] is False
    assert states["selected_pressure_public_value_path_repair_required"] is True
    assert states["oscillatory_ready"] is True
    assert states["correction_ingest_allowed"] is True
    assert states["compact_radial_stress_operator_independently_admitted"] is True
    assert states["independent_spacetime_radial_stress_generalization_admitted"] is True
    assert states["finite_correction_cycle_radial_operator_reuse_allowed_without_retuning"] is True

    assert states["leading_ready"] is False
    assert states["same_cycle_requested_stress_materialized"] is False
    assert states["signed_mean_inverse_input_ready"] is False
    assert states["public_velocity_correction_materialized"] is False
    assert states["correction_ready"] is False
    assert states["candidate_artifact_instantiated"] is False
    assert states["velocity_export_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False
    assert checkpoint["formal_gates_unchanged"] == FORMAL_GATES


def test_v45_pins_latest_agent3_and_agent4_execution_identity(checkpoint):
    a3 = checkpoint["upstream"]["agent3_spacetime_radial_stress_admission"]
    a4 = checkpoint["upstream"]["agent4_selected_pressure_audit"]

    assert a3["head"] == AGENT3_SPACETIME_ADMISSION_HEAD
    assert a3["dedicated_workflow_run"] == 35435350975
    assert a3["standard_workflow_run"] == 35435350972
    assert a3["artifact_id"] == 10581334240
    assert a3["artifact_zip_digest"] == (
        "sha256:68dbb652b633092b766d7ec7ea089d937ed26c4d05d2b45ba4d29fc3077b5d3a"
    )

    assert a4["head"] == AGENT4_PRESSURE_AUDIT_HEAD
    assert a4["dedicated_workflow_run"] == 35435721364
    assert a4["artifact_id"] == AGENT4_PRESSURE_ARTIFACT_ID
    assert a4["artifact_zip_digest"] == (
        "sha256:408892549b47495142aab4c1055a7a4a40cca035a75e189f87757b68398b614e"
    )
    assert a4["selected_pressure_primitive_independent_preflight_passed"] is False
    assert set(a4["failed_guards"]) == EXPECTED_FAILED_GUARDS


def test_v45_freezes_substantive_p_eta_failure_and_separates_roundoff_failures(checkpoint):
    a4 = checkpoint["upstream"]["agent4_selected_pressure_audit"]
    p_eta = a4["p_eta"]
    p_y = a4["p_Y"]
    phi_eta = a4["Phi_eta"]

    errors = [row["relative_rms"] for row in p_eta["levels"]]
    assert errors == pytest.approx(
        [0.9991143298529059, 0.9039058221188769, 1.226775129471845]
    )
    assert p_eta["levels"][-1]["relative_max"] == pytest.approx(1.5890750527047752)

    # Freeze the scientifically relevant fact rather than a rounded PR-body
    # representation of derived ratios: the first refinement barely improves,
    # while the second gets worse.  Also verify the receipt ratios are exactly
    # the ratios implied by the live audit's own RMS ladder.
    expected_ratios = [errors[0] / errors[1], errors[1] / errors[2]]
    assert p_eta["refinement_ratios"] == pytest.approx(expected_ratios, rel=1e-12, abs=1e-12)
    assert 1.0 < p_eta["refinement_ratios"][0] < 1.2
    assert p_eta["refinement_ratios"][1] < 1.0
    assert a4["nontrivial_reference_rms"]["p_eta"] > 1.0e15

    assert p_y["levels"][-1]["relative_rms"] < 1.0e-12
    assert phi_eta["levels"][-1]["relative_rms"] < 1.0e-5
    interpretation = checkpoint["pressure_audit_interpretation"]
    assert interpretation["p_Y_refinement_guard_failed_at_machine_floor"] is True
    assert interpretation["Y_Phi_Y_stationary_center_cloud_not_binary64_identifiable"] is True
    assert interpretation["p_eta_is_substantive_blocker"] is True
    assert interpretation["p_eta_failure_attributed_to_formula_error"] is False


@pytest.mark.parametrize(
    "key",
    [
        "leading_ready",
        "selected_pressure_primitive_independent_preflight_passed",
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
def test_v45_rejects_overpromotion_even_if_resigned(checkpoint, key):
    tampered = copy.deepcopy(checkpoint)
    tampered["states"][key] = True
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError):
        validate_checkpoint(tampered)


def test_v45_rejects_pressure_audit_identity_substitution_even_if_resigned(checkpoint):
    tampered = copy.deepcopy(checkpoint)
    tampered["upstream"]["agent4_selected_pressure_audit"]["artifact_id"] += 1
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="wrong Agent-4 pressure audit artifact"):
        validate_checkpoint(tampered)


def test_v45_rejects_failed_guard_laundering_even_if_resigned(checkpoint):
    tampered = copy.deepcopy(checkpoint)
    tampered["upstream"]["agent4_selected_pressure_audit"]["failed_guards"] = []
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="pressure failed-guard set changed"):
        validate_checkpoint(tampered)


def test_v45_rejects_threshold_change_even_if_resigned(checkpoint):
    tampered = copy.deepcopy(checkpoint)
    tampered["formal_gates_unchanged"]["held_out_normalized_momentum_max"] = 2.0e-3
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="formal final gates changed"):
        validate_checkpoint(tampered)
