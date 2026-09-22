from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_main_pulse_amplitude_scope import (
    CONTRACT_PATH,
    CONSTRAINTS_PATH,
    PROJECT_STATUS_PATH,
    SOURCE_PATH,
    amplitude_dependence_witness,
    audit_contract,
)

ROOT = Path(__file__).resolve().parents[1]


def _live_inputs():
    contract = json.loads((ROOT / CONTRACT_PATH).read_text())
    constraints = json.loads((ROOT / CONSTRAINTS_PATH).read_text())
    status = json.loads((ROOT / PROJECT_STATUS_PATH).read_text())
    source_text = (ROOT / SOURCE_PATH).read_text()
    return contract, source_text, constraints, status


def test_live_repository_scope_is_fail_closed_and_clean():
    contract, source_text, constraints, status = _live_inputs()
    assert audit_contract(contract, source_text, constraints, status) == []


def test_mechanics_witness_is_nonvacuous_and_not_scientific_evidence():
    w = amplitude_dependence_witness()
    assert w["same_public_shape_inputs"] is True
    assert w["different_executable_U"] is True
    assert w["U_a"] == pytest.approx(0.6)
    assert w["U_b"] == pytest.approx(0.675)
    assert w["not_source_numeric_evidence"] is True
    assert w["not_openai_numeric_evidence"] is True
    assert w["not_pde_evidence"] is True


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("current_truth", "source_exact_amplitude_root_materialized"), True),
        (("current_truth", "source_pulse_end_MJ_corrections_materialized"), True),
        (("current_truth", "source_hidden_parameters_recovered"), True),
        (("current_truth", "current_cartesian_pulse_velocity_composed"), True),
        (("current_truth", "terminal_global_velocity_materialized"), True),
        (("current_truth", "velocity_export_ready_for_kokuno_route"), True),
        (("current_truth", "visual_correspondence_verified"), True),
        (("current_truth", "pde_validated"), True),
        (("current_truth", "paper_exact"), True),
        (("current_truth", "openai_field_identified"), True),
        (("current_truth", "blowup_proved"), True),
        (("scope", "candidate_bytes_changed"), True),
        (("scope", "forcing_changed"), True),
        (("scope", "scientific_threshold_changed"), True),
    ],
)
def test_truth_and_scope_promotions_fail_closed(path, value):
    contract, source_text, constraints, status = _live_inputs()
    mutated = copy.deepcopy(contract)
    target = mutated
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    assert audit_contract(mutated, source_text, constraints, status)


def test_autonomous_amplitude_cannot_be_laundered_into_public_source_bucket():
    contract, source_text, constraints, status = _live_inputs()
    mutated = copy.deepcopy(contract)
    autonomous = mutated["provenance_classes"]["autonomous_design"]
    moved = next(x for x in autonomous if "A_principal=sqrt" in x)
    autonomous.remove(moved)
    mutated["provenance_classes"]["public_source_fact"].append(moved)
    violations = audit_contract(mutated, source_text, constraints, status)
    assert any("laundered" in x or "classification" in x for x in violations)


def test_source_exact_amplitude_must_remain_pending_until_materialized():
    contract, source_text, constraints, status = _live_inputs()
    mutated = copy.deepcopy(contract)
    mutated["provenance_classes"]["pending_unknown"] = [
        x
        for x in mutated["provenance_classes"]["pending_unknown"]
        if "full source Amp(eta) root" not in x
    ]
    assert audit_contract(mutated, source_text, constraints, status)


def test_amplitude_replacement_must_create_new_identity_and_revalidate_evidence():
    contract, source_text, constraints, status = _live_inputs()
    mutated = copy.deepcopy(contract)
    mutated["promotion_requirements"]["evidence_transfer_after_amplitude_replacement"] = [
        "shared R_0 permits automatic receipt inheritance"
    ]
    violations = audit_contract(mutated, source_text, constraints, status)
    assert any("new identity" in x or "evidence-transfer" in x for x in violations)


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("thresholds", "pde_residual_max", 0.01),
        ("thresholds", "pde_residual_L2", 0.01),
        ("thresholds", "divergence_max", 1e-4),
        ("thresholds", "divergence_L2", 1e-4),
    ],
)
def test_cr001_threshold_relaxation_is_rejected(section, key, value):
    contract, source_text, constraints, status = _live_inputs()
    mutated_constraints = copy.deepcopy(constraints)
    mutated_constraints["validation"][section][key] = value
    assert audit_contract(contract, source_text, mutated_constraints, status)


def test_residual_defined_free_force_drift_is_rejected():
    contract, source_text, constraints, status = _live_inputs()
    mutated_constraints = copy.deepcopy(constraints)
    mutated_constraints["forcing"]["mode"] = "residual_defined_pointwise"
    mutated_constraints["forcing"]["restriction"] = "fit f=R pointwise"
    violations = audit_contract(contract, source_text, mutated_constraints, status)
    assert any("forcing" in x.lower() for x in violations)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("velocity_export_ready", False),
        ("visualization_ready", True),
        ("visual_correspondence_verified", True),
        ("pde_validated", True),
        ("paper_exact", True),
        ("openai_field_identified", True),
    ],
)
def test_canonical_eq45_independent_state_cannot_be_rewritten_by_kokuno_scope(key, value):
    contract, source_text, constraints, status = _live_inputs()
    mutated_status = copy.deepcopy(status)
    mutated_status["states"][key] = value
    assert audit_contract(contract, source_text, constraints, mutated_status)


def test_source_token_drift_is_rejected():
    contract, source_text, constraints, status = _live_inputs()
    mutated_source = source_text.replace(
        '"source_exact_amplitude_root_materialized": False',
        '"source_exact_amplitude_root_materialized": True',
        1,
    )
    violations = audit_contract(contract, mutated_source, constraints, status)
    assert any("source token" in x for x in violations)
