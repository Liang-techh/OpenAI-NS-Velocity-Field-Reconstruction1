import copy

import pytest

from openai_ns_reconstruction.kokuno_independent_radial_stress_routing_checkpoint import (
    FORMAL_GATES,
    _canonical_sha256,
    build_checkpoint,
    validate_checkpoint,
)


def test_v44_admits_only_independently_audited_radial_stress_operator():
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    states = checkpoint["states"]
    assert states["oscillatory_ready"] is True
    assert states["correction_ingest_allowed"] is True
    assert states["oscillatory_requested_stress_component_materialized"] is True
    assert states["oscillatory_component_finite_head_mean_debt_materialized"] is True
    assert states["independent_compact_radial_stress_cross_audit_passed"] is True
    assert states["compact_radial_stress_operator_independently_admitted"] is True
    assert states["full_composite_radial_stress_execution_allowed_when_actual_defect_available"] is True
    assert states["same_cycle_requested_stress_materialized"] is False
    assert states["candidate_numeric_finite_head_mean_debt_materialized"] is False
    assert states["signed_mean_inverse_input_ready"] is False
    assert states["public_velocity_correction_materialized"] is False
    assert states["correction_ready"] is False
    assert states["candidate_artifact_instantiated"] is False
    assert states["velocity_export_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False
    assert checkpoint["formal_gates_unchanged"] == FORMAL_GATES


def test_v44_pins_exact_independent_radial_stress_execution_and_metrics():
    checkpoint = build_checkpoint()
    radial = checkpoint["upstream"]["agent3_independent_radial_stress_admission"]
    assert radial["head"] == "a97f0882ab0b1fcc0d4b9575b5366fb44b35c233"
    assert radial["dedicated_workflow_run"] == 35432615191
    assert radial["standard_workflow_run"] == 35432615180
    assert radial["artifact_id"] == 10581435351
    assert radial["artifact_zip_digest"] == (
        "sha256:be1c19d3257086622f90ef5a8b7be9181dcbe9628d568eb8a4f9bfeac213ec6c"
    )
    assert radial["agent4_audit_head"] == "19296acad4f84ba05d3b095bb8c01bf5a2c95891"
    assert radial["agent4_dedicated_workflow_run"] == 35431599990
    assert radial["agent4_standard_workflow_run"] == 35431599977
    assert radial["failed_guards"] == []
    assert radial["theta_e2"]["finest_relative_rms"] == pytest.approx(
        0.00010215080513304861
    )
    assert radial["axial_e1"]["finest_relative_rms"] == pytest.approx(
        0.0001864344008886948
    )
    assert radial["theta_e2"]["refinement_ratios"] == pytest.approx(
        [12.029826646162858, 14.923789334921668]
    )
    assert radial["axial_e1"]["refinement_ratios"] == pytest.approx(
        [14.554891779445125, 14.786694125947507]
    )


def test_v44_keeps_agent1_diagnostic_algebra_fail_closed_as_source_leading():
    checkpoint = build_checkpoint()
    leading = checkpoint["upstream"]["agent1_latest_leading"]
    assert leading["exact_displayed_R1_R2_monomial_bookkeeping_executable"] is True
    assert leading["one_factor_at_a_time_lipschitz_bookkeeping_executable"] is True
    assert leading["source_R1_radius_one_ball_norm_machine_bound"] is False
    assert leading["source_R2_radius_one_ball_norm_machine_bound"] is False
    assert leading["source_R1_radius_one_ball_lipschitz_machine_bound"] is False
    assert leading["source_R2_radius_one_ball_lipschitz_machine_bound"] is False
    assert leading["source_operator_M_K_machine_bound"] is False
    assert leading["source_B0_T_sh_machine_bound"] is False
    assert leading["global_leading_velocity_pressure_ready"] is False


@pytest.mark.parametrize(
    "key",
    [
        "leading_ready",
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
def test_v44_rejects_downstream_overpromotion_even_if_resigned(key):
    checkpoint = build_checkpoint()
    tampered = copy.deepcopy(checkpoint)
    tampered["states"][key] = True
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError):
        validate_checkpoint(tampered)


def test_v44_rejects_threshold_change_even_if_resigned():
    checkpoint = build_checkpoint()
    tampered = copy.deepcopy(checkpoint)
    tampered["formal_gates_unchanged"]["held_out_normalized_momentum_max"] = 2.0e-3
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="formal final gates changed"):
        validate_checkpoint(tampered)


def test_v44_rejects_radial_stress_identity_substitution_even_if_resigned():
    checkpoint = build_checkpoint()
    tampered = copy.deepcopy(checkpoint)
    radial = tampered["upstream"]["agent3_independent_radial_stress_admission"]
    radial["artifact_id"] += 1
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="wrong Agent-3 radial-stress admission artifact"):
        validate_checkpoint(tampered)


def test_v44_rejects_component_to_full_correction_laundering_even_if_resigned():
    checkpoint = build_checkpoint()
    tampered = copy.deepcopy(checkpoint)
    tampered["truth_boundary"][
        "component_radial_stress_laundered_as_full_same_cycle_correction"
    ] = True
    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="cannot be laundered as full correction"):
        validate_checkpoint(tampered)
