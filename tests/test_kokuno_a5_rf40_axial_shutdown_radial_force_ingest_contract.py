from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_rf40_axial_shutdown_radial_stress_ingest_contract as parent
from openai_ns_reconstruction import kokuno_a5_rf40_axial_shutdown_radial_force_ingest_contract as subject


def test_receipt_is_deterministic_and_enforced() -> None:
    first = subject.materialize_registration_receipt()
    second = subject.materialize_registration_receipt()
    assert first == second
    assert first["schema"] == subject.SCHEMA_NAME
    assert first["task_id"] == subject.TASK_ID
    subject.enforce_registration_receipt(first)


def test_exact_parent_and_agent3_radial_force_provenance_is_frozen() -> None:
    r = subject.materialize_registration_receipt()["registration"]
    assert r["parent_a5"]["pr"] == 1019
    assert r["parent_a5"]["head"] == "b5b9b30f6a209e13f22eb8aa3b7f0a251637cc42"
    assert r["parent_a5"]["source_blob"] == "b5d3c845dd8daffc7cd3dc2634d7532acd1435b6"
    a3 = r["agent3_axial_shutdown_radial_force"]
    assert a3["pr"] == 1022
    assert a3["head"] == "93200f28deae4372c2338e0620b785543af892d4"
    assert a3["source_blob"] == "222e75ab5ecb1f0625ed3b25d4e2defe6f780b8f"
    assert a3["test_blob"] == "5c4aeb46ad229e94c43a3a5ebc94746363eea01a"
    assert a3["workflow_blob"] == "a70cf9b5d04d368ad886cce3fe4a35f0521a4c30"
    assert a3["parent_stress_pr"] == 1015
    assert a3["parent_stress_head"] == parent.AGENT3_AXIAL_SHUTDOWN_RADIAL_STRESS["head"]


def test_radial_force_advance_and_no_a4_force_audit_are_both_locked() -> None:
    r = subject.materialize_registration_receipt()["registration"]
    a3 = r["agent3_axial_shutdown_radial_force"]
    a4 = r["agent4_stress_audit_only"]
    truth = r["truth_boundary"]
    assert a3["source_formula"] == "(div T)_r = partial_z sigma_1"
    assert a3["z_derivative_step_ladder"] == [0.02, 0.01, 0.005]
    assert a3["radial_force_through_X_2_materialized"] is True
    assert a3["independent_third_radial_moment_inverse_introduced"] is False
    assert a3["recorded_radial_mean_equated_to_radial_force"] is False
    assert a3["authorized_as_ns_correction_target"] is False
    assert a3["consumes_parent_x4_composite"] is False
    assert a3["candidate_residual_evidence"] is False
    assert a4["pr"] == 1017
    assert a4["audited_agent3_pr"] == 1015
    assert a4["audits_agent3_1022_radial_force"] is False
    assert a4["radial_force_independent_audit_authority"] is False
    assert truth["current_nonlinear_radial_stress_through_axial_shutdown_materialized"] is True
    assert truth["current_nonlinear_radial_force_through_axial_shutdown_materialized"] is True
    assert truth["agent3_1022_axial_shutdown_radial_force_registered"] is True
    assert truth["agent4_1017_stress_audit_does_not_audit_radial_force"] is True
    assert truth["agent4_dedicated_axial_shutdown_radial_force_audit_present"] is False
    assert truth["agent4_independent_axial_shutdown_radial_force_audit_registered"] is False
    assert truth["agent4_independent_axial_shutdown_radial_force_audit_admitted"] is False
    assert truth["scoped_axial_shutdown_radial_force_authorized_as_ns_correction_target"] is False


def test_final_gates_baseline_and_readiness_are_unchanged() -> None:
    r = subject.materialize_registration_receipt()["registration"]
    assert r["final_gate"] == parent.FINAL_GATE
    assert r["final_gate"]["normalized_momentum_sampled_max"] == 1.0e-3
    assert r["final_gate"]["normalized_momentum_volume_l2"] == 1.0e-3
    assert r["final_gate"]["divergence_sampled_max"] == 1.0e-5
    assert r["final_gate"]["divergence_volume_l2"] == 1.0e-5
    assert r["final_gate"]["canonical_volume_quadrature_ladder"] == [24, 48, 96]
    assert r["st006_baseline"] == parent.ST006_BASELINE
    assert r["readiness"] == {
        "leading_ready": False,
        "oscillatory_ready": True,
        "correction_ready": False,
        "velocity_export_ready": False,
        "pde_validated": False,
    }
    assert r["frozen_science"]["residual_defined_free_forcing_forbidden"] is True
    truth = r["truth_boundary"]
    assert truth["complete_ns_defect_materialized"] is False
    assert truth["real_agent3_ns_correction_velocity_materialized"] is False
    assert truth["real_candidate_finite_correction_cycle_run"] is False
    assert truth["heldout_normalized_ns_residual_assessed"] is False
    assert truth["same_protocol_st006_comparison_available_now"] is False
    assert truth["residual_reduction_claimed"] is False
    assert truth["pde_validated"] is False


@pytest.mark.parametrize(
    "mutator",
    [
        lambda x: x["registration"]["agent3_axial_shutdown_radial_force"].__setitem__("head", "0" * 40),
        lambda x: x["registration"]["agent3_axial_shutdown_radial_force"].__setitem__("authorized_as_ns_correction_target", True),
        lambda x: x["registration"]["agent3_axial_shutdown_radial_force"].__setitem__("consumes_parent_x4_composite", True),
        lambda x: x["registration"]["agent4_stress_audit_only"].__setitem__("audits_agent3_1022_radial_force", True),
        lambda x: x["registration"]["truth_boundary"].__setitem__("agent4_dedicated_axial_shutdown_radial_force_audit_present", True),
        lambda x: x["registration"]["truth_boundary"].__setitem__("scoped_axial_shutdown_radial_force_authorized_as_ns_correction_target", True),
        lambda x: x["registration"]["readiness"].__setitem__("correction_ready", True),
        lambda x: x["registration"]["readiness"].__setitem__("pde_validated", True),
        lambda x: x["registration"]["final_gate"].__setitem__("normalized_momentum_sampled_max", 2.0e-3),
    ],
)
def test_mutations_fail_closed(mutator) -> None:
    receipt = copy.deepcopy(subject.materialize_registration_receipt())
    mutator(receipt)
    registration = receipt["registration"]
    receipt["registration_sha256"] = subject._sha256(registration)
    with pytest.raises(ValueError):
        subject.enforce_registration_receipt(receipt)


def test_write_roundtrip(tmp_path) -> None:
    path = subject.write_registration_receipt(tmp_path / "receipt.json")
    loaded = json.loads(path.read_text(encoding="utf-8"))
    subject.enforce_registration_receipt(loaded)
