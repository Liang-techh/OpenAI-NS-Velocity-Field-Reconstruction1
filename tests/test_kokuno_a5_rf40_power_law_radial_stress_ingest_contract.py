from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_rf40_power_law_radial_stress_ingest_contract as reg


def _rehash(receipt: dict) -> dict:
    out = copy.deepcopy(receipt)
    out["registration_sha256"] = reg._sha256(out["registration"])
    return out


def test_receipt_registers_exact_x4_stress_and_matching_independent_audit_without_science_promotion() -> None:
    receipt = reg.materialize_registration_receipt()
    reg.enforce_registration_receipt(receipt)
    reg.enforce_registration_receipt(json.loads(json.dumps(receipt, sort_keys=True)))

    r = receipt["registration"]
    stress = r["agent3_power_law_stress"]
    mean_audit = r["agent4_power_law_mean_audit"]
    stress_audit = r["agent4_power_law_stress_audit"]
    truth = r["truth_boundary"]

    assert r["parent_a5"]["pr"] == 1030
    assert stress["pr"] == 1037
    assert stress["parent_agent3_pr"] == 1028
    assert mean_audit["pr"] == 1029
    assert mean_audit["audited_agent3_pr"] == stress["parent_agent3_pr"]
    assert mean_audit["audits_agent3_1037_radial_stress"] is False
    assert stress_audit["pr"] == 1038
    assert stress_audit["audited_agent3_pr"] == stress["pr"]
    assert stress_audit["audited_agent3_head"] == stress["head"]
    assert stress_audit["audited_agent3_source_blob"] == stress["source_blob"]
    assert stress_audit["implementation_distinct"] is True
    assert stress_audit["independent_operator"] != stress_audit["production_operator"]
    assert stress["radial_force_requires_later_partial_z_sigma_1"] is True

    assert truth["current_nonlinear_mean_through_power_law_materialized"] is True
    assert truth["agent3_1037_power_law_radial_stress_registered"] is True
    assert truth["current_nonlinear_radial_stress_through_power_law_materialized"] is True
    assert truth["agent4_dedicated_power_law_radial_stress_audit_present"] is True
    assert truth["agent4_independent_power_law_radial_stress_audit_registered"] is True
    assert truth["agent4_independent_power_law_radial_stress_audit_admitted"] is False
    assert truth["current_nonlinear_radial_force_through_power_law_materialized"] is False
    assert truth["scoped_power_law_radial_stress_authorized_as_ns_correction_target"] is False
    assert r["readiness"]["correction_ready"] is False
    assert r["readiness"]["velocity_export_ready"] is False
    assert r["readiness"]["pde_validated"] is False


def test_rejects_wrong_stress_identity_even_if_rehashed() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent3_power_law_stress"]["head"] = "0" * 40
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_parent_mean_lineage_or_mean_audit_scope_drift() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent3_power_law_stress"]["parent_agent3_head"] = "1" * 40
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))

    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent4_power_law_mean_audit"]["audits_agent3_1037_radial_stress"] = True
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_matching_stress_audit_target_or_operator_drift() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent4_power_law_stress_audit"]["audited_agent3_head"] = "2" * 40
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))

    receipt = reg.materialize_registration_receipt()
    audit = receipt["registration"]["agent4_power_law_stress_audit"]
    audit["independent_operator"] = audit["production_operator"]
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_scoped_stress_or_audit_as_complete_ns_evidence() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent3_power_law_stress"]["authorized_as_ns_correction_target"] = True
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))

    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent4_power_law_stress_audit"]["complete_ns_residual_evidence"] = True
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_premature_stress_audit_admission_force_or_correction_promotion() -> None:
    for key in (
        "agent4_independent_power_law_radial_stress_audit_admitted",
        "current_nonlinear_radial_force_through_power_law_materialized",
        "scoped_power_law_radial_stress_authorized_as_ns_correction_target",
    ):
        receipt = reg.materialize_registration_receipt()
        receipt["registration"]["truth_boundary"][key] = True
        with pytest.raises(ValueError):
            reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_complete_ns_or_pde_promotion() -> None:
    for key in (
        "matched_cartesian_pressure_materialized",
        "restricted_forcing_materialized",
        "complete_ns_defect_materialized",
        "real_agent3_ns_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "pde_validated",
    ):
        receipt = reg.materialize_registration_receipt()
        receipt["registration"]["truth_boundary"][key] = True
        with pytest.raises(ValueError):
            reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_final_gate_baseline_or_free_forcing_drift() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["final_gate"]["normalized_momentum_sampled_max"] = 2e-3
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))

    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["st006_baseline"]["momentum_sampled_max"] = 0.01
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))

    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["frozen_science"]["residual_defined_free_forcing_forbidden"] = False
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))
