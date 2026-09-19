import copy

import pytest

from openai_ns_reconstruction.kokuno_oscillatory_mean_debt_routing_checkpoint import (
    FORMAL_GATES,
    build_checkpoint,
    validate_checkpoint,
)


def test_v43_promotes_only_component_mean_debt():
    checkpoint = build_checkpoint()
    validate_checkpoint(checkpoint)
    states = checkpoint["states"]
    assert states["oscillatory_ready"] is True
    assert states["public_oscillatory_time_derivative_independently_validated"] is True
    assert states["correction_ingest_allowed"] is True
    assert states["oscillatory_requested_stress_component_materialized"] is True
    assert states["oscillatory_component_finite_head_mean_debt_materialized"] is True
    assert states["same_cycle_requested_stress_materialized"] is False
    assert states["candidate_numeric_finite_head_mean_debt_materialized"] is False
    assert states["signed_mean_inverse_input_ready"] is False
    assert states["correction_ready"] is False
    assert states["candidate_artifact_instantiated"] is False
    assert states["velocity_export_ready"] is False
    assert states["formal_full_domain_pde_gate_assessed"] is False
    assert states["pde_validated"] is False
    assert checkpoint["formal_gates_unchanged"] == FORMAL_GATES


def test_v43_pins_real_a3_receipt_and_not_formal_missing_weight():
    checkpoint = build_checkpoint()
    a3 = checkpoint["upstream"]["agent3_oscillatory_component_mean_debt"]
    assert a3["dedicated_workflow_conclusion"] == "success"
    assert a3["standard_workflow_conclusion"] == "success"
    assert a3["artifact_id"] == 10580380247
    assert a3["artifact_zip_digest"].startswith("sha256:")
    assert a3["requested_stress_rms"] > 0.0
    assert a3["oscillatory_component_debt_rms"] > 0.0
    assert a3["identity_closure_max_abs"] <= 1e-12
    assert a3["real_oscillatory_self_defect_component_consumed"] is True
    assert a3["autonomous_factor_is_formal_theorem_missing_weight"] is False
    assert a3["full_same_cycle_composite_requested_stress_materialized"] is False
    assert a3["candidate_full_finite_head_mean_debt_materialized"] is False


def test_v43_keeps_latest_in_progress_radial_audit_fail_closed():
    checkpoint = build_checkpoint()
    later = checkpoint["upstream"]["agent3_later_radial_closure_audit"]
    assert later["status_at_checkpoint_freeze"] == "in_progress"
    assert later["scientific_promotion_allowed"] is False


@pytest.mark.parametrize(
    "key",
    [
        "same_cycle_requested_stress_materialized",
        "candidate_numeric_finite_head_mean_debt_materialized",
        "signed_mean_inverse_input_ready",
        "correction_ready",
        "finite_correction_cycle_run",
        "candidate_artifact_instantiated",
        "velocity_export_ready",
        "formal_full_domain_pde_gate_assessed",
        "heldout_ns_residual_assessed",
        "pde_validated",
    ],
)
def test_v43_rejects_downstream_overpromotion(key):
    checkpoint = build_checkpoint()
    tampered = copy.deepcopy(checkpoint)
    tampered["states"][key] = True
    # Re-sign to prove the semantic guard, not only hash checking, rejects it.
    from openai_ns_reconstruction.kokuno_oscillatory_mean_debt_routing_checkpoint import _canonical_sha256

    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError):
        validate_checkpoint(tampered)


def test_v43_rejects_threshold_change_even_if_resigned():
    checkpoint = build_checkpoint()
    tampered = copy.deepcopy(checkpoint)
    tampered["formal_gates_unchanged"]["held_out_normalized_momentum_max"] = 2e-3
    from openai_ns_reconstruction.kokuno_oscillatory_mean_debt_routing_checkpoint import _canonical_sha256

    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="formal final gates changed"):
        validate_checkpoint(tampered)


def test_v43_rejects_component_to_full_laundering_even_if_resigned():
    checkpoint = build_checkpoint()
    tampered = copy.deepcopy(checkpoint)
    a3 = tampered["upstream"]["agent3_oscillatory_component_mean_debt"]
    a3["full_same_cycle_composite_requested_stress_materialized"] = True
    from openai_ns_reconstruction.kokuno_oscillatory_mean_debt_routing_checkpoint import _canonical_sha256

    tampered["checkpoint_sha256"] = _canonical_sha256(tampered)
    with pytest.raises(ValueError, match="laundered"):
        validate_checkpoint(tampered)
