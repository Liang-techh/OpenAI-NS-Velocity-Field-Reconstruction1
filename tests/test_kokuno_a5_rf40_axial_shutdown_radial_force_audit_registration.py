from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import (
    kokuno_a5_rf40_axial_shutdown_radial_force_audit_registration as reg,
)


def _rehash(receipt: dict) -> dict:
    out = copy.deepcopy(receipt)
    out["registration_sha256"] = reg._sha256(out["registration"])
    return out


def test_materialized_receipt_registers_exact_a4_audit_without_science_promotion() -> None:
    receipt = reg.materialize_registration_receipt()
    reg.enforce_registration_receipt(receipt)
    roundtripped = json.loads(json.dumps(receipt, sort_keys=True))
    reg.enforce_registration_receipt(roundtripped)

    r = receipt["registration"]
    a3 = r["agent3_radial_force"]
    a4 = r["agent4_radial_force_audit"]
    truth = r["truth_boundary"]

    assert r["parent_a5"]["pr"] == 1024
    assert a3["pr"] == 1022
    assert a4["pr"] == 1023
    assert a4["audited_agent3_pr"] == a3["pr"]
    assert a4["audited_agent3_head"] == a3["head"]
    assert a4["audited_agent3_source_blob"] == a3["source_blob"]
    assert a4["implementation_distinct"] is True
    assert a4["independent_operator"] != a4["production_operator"]

    assert truth["agent4_dedicated_axial_shutdown_radial_force_audit_present"] is True
    assert truth["agent4_independent_axial_shutdown_radial_force_audit_registered"] is True
    assert truth["agent4_independent_axial_shutdown_radial_force_audit_admitted"] is False
    assert truth["scoped_axial_shutdown_radial_force_authorized_as_ns_correction_target"] is False
    assert r["readiness"]["correction_ready"] is False
    assert r["readiness"]["velocity_export_ready"] is False
    assert r["readiness"]["pde_validated"] is False
    assert r["frozen_science"]["residual_defined_free_forcing_forbidden"] is True


def test_rejects_wrong_a4_identity_even_if_receipt_is_rehashed() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent4_radial_force_audit"]["head"] = "0" * 40
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_a4_target_drift_even_if_receipt_is_rehashed() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent4_radial_force_audit"]["audited_agent3_head"] = "1" * 40
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_collapse_to_production_operator() -> None:
    receipt = reg.materialize_registration_receipt()
    a4 = receipt["registration"]["agent4_radial_force_audit"]
    a4["independent_operator"] = a4["production_operator"]
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_scoped_audit_as_complete_ns_evidence() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["agent4_radial_force_audit"]["complete_ns_residual_evidence"] = True
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_scoped_force_correction_authorization() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["truth_boundary"][
        "scoped_axial_shutdown_radial_force_authorized_as_ns_correction_target"
    ] = True
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_premature_a4_scientific_admission() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["truth_boundary"][
        "agent4_independent_axial_shutdown_radial_force_audit_admitted"
    ] = True
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_final_gate_or_baseline_drift() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["final_gate"]["normalized_momentum_max"] = 2.0e-3
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))

    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["st006_baseline"]["momentum_sampled_max"] = 0.01
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))


def test_rejects_free_forcing_or_pde_promotion() -> None:
    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["frozen_science"]["residual_defined_free_forcing_forbidden"] = False
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))

    receipt = reg.materialize_registration_receipt()
    receipt["registration"]["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        reg.enforce_registration_receipt(_rehash(receipt))
