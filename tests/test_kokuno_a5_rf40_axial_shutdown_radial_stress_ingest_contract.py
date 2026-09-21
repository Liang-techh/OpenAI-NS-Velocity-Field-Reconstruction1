from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction import kokuno_a5_current_rf40_power_law_composite_ingest_contract as parent
from openai_ns_reconstruction import kokuno_a5_rf40_axial_shutdown_radial_stress_ingest_contract as subject


def test_receipt_is_deterministic_and_enforced() -> None:
    first = subject.materialize_registration_receipt()
    second = subject.materialize_registration_receipt()
    assert first == second
    assert first["schema"] == subject.SCHEMA_NAME
    assert first["task_id"] == subject.TASK_ID
    subject.enforce_registration_receipt(first)


def test_exact_parent_and_sibling_provenance_is_frozen() -> None:
    r = subject.materialize_registration_receipt()["registration"]
    assert r["parent_a5"]["pr"] == 1013
    assert r["parent_a5"]["head"] == "df75e6994ca4ae1fa405d56ef162b5a89d5e6d58"
    assert r["parent_a5"]["source_blob"] == "3d7c80f5f8e7f46f3da39402ea163f13b9120a8e"
    assert r["agent3_axial_shutdown_radial_stress"]["pr"] == 1015
    assert r["agent3_axial_shutdown_radial_stress"]["head"] == "2e7683f06bea810881afbe39d8f0cf64815ddd56"
    assert r["agent3_axial_shutdown_radial_stress"]["source_blob"] == "9c6082e6dc5edacf7b6a64e87549cbec142b3d4a"
    assert r["agent4_axial_shutdown_radial_stress_audit"]["pr"] == 1017
    assert r["agent4_axial_shutdown_radial_stress_audit"]["head"] == "cc3109bd17d518f3414d972f2303a20f2dabc1b1"
    assert r["agent4_axial_shutdown_radial_stress_audit"]["source_blob"] == "c886a7e9d222313125f5e242ca9cc9828babdb23"


def test_lineage_and_scope_remain_fail_closed() -> None:
    r = subject.materialize_registration_receipt()["registration"]
    a3 = r["agent3_axial_shutdown_radial_stress"]
    a4 = r["agent4_axial_shutdown_radial_stress_audit"]
    truth = r["truth_boundary"]
    assert a3["parent_mean_pr"] == 1011
    assert a3["agent2_composite_pr"] == 999
    assert a3["agent1_leading_pr"] == 993
    assert a3["consumes_parent_x4_composite"] is False
    assert a3["radial_force_through_X_2_materialized"] is False
    assert a4["audited_agent3_head"] == a3["head"]
    assert a4["uses_a3_production_radial_inverse_for_reference"] is False
    assert a4["momentum_or_full_ns_evidence"] is False
    assert truth["current_nonlinear_radial_stress_through_axial_shutdown_materialized"] is True
    assert truth["current_nonlinear_radial_force_through_axial_shutdown_materialized"] is False
    assert truth["current_nonlinear_radial_stress_through_power_law_materialized"] is False
    assert truth["agent4_1017_independent_axial_shutdown_radial_stress_audit_registered"] is True
    assert truth["agent4_1017_independent_axial_shutdown_radial_stress_audit_admitted"] is False


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


@pytest.mark.parametrize(
    "mutator",
    [
        lambda x: x["registration"]["agent3_axial_shutdown_radial_stress"].__setitem__("head", "0" * 40),
        lambda x: x["registration"]["agent4_axial_shutdown_radial_stress_audit"].__setitem__("head", "1" * 40),
        lambda x: x["registration"]["truth_boundary"].__setitem__("current_nonlinear_radial_force_through_axial_shutdown_materialized", True),
        lambda x: x["registration"]["truth_boundary"].__setitem__("agent4_1017_independent_axial_shutdown_radial_stress_audit_admitted", True),
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
