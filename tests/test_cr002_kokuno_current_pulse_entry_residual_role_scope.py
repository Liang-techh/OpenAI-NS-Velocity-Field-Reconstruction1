from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_current_pulse_entry_residual_role_scope import (
    assert_contract,
    audit_contract,
    bookkeeping_projection_witness,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "configs" / "kokuno_current_pulse_entry_residual_role_scope.json"
CONSTRAINTS = ROOT / "configs" / "constraints.json"
CANDIDATE_MODULE = (
    ROOT / "src" / "openai_ns_reconstruction" / "kokuno_current_pulse_entry_moments.py"
)
PROJECT_STATUS = ROOT / "project_status.json"


def _payload() -> dict:
    return json.loads(CONTRACT.read_text())


def _audit_temp(tmp_path: Path, payload: dict) -> list[str]:
    path = tmp_path / "contract.json"
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return audit_contract(path, CONSTRAINTS, CANDIDATE_MODULE, PROJECT_STATUS)


def test_exact_parent_contract_audits_clean() -> None:
    assert_contract(CONTRACT, CONSTRAINTS, CANDIDATE_MODULE, PROJECT_STATUS)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("role_separation", "construction_closure_residual_is_public_source_data"), True),
        (("role_separation", "construction_closure_residual_is_forcing"), True),
        (("role_separation", "construction_closure_residual_is_heldout_ns_residual"), True),
        (("role_separation", "construction_closure_residual_may_authorize_force_fitting"), True),
        (("role_separation", "residual_defined_free_forcing_allowed"), True),
        (("role_separation", "current_cartesian_end_compensation_composed"), True),
        (("role_separation", "unified_global_kokuno_velocity_export_ready"), True),
        (("independent_delivery_states", "kokuno_visual_correspondence_verified"), True),
        (("independent_delivery_states", "kokuno_pde_validated"), True),
        (("independent_delivery_states", "kokuno_paper_exact"), True),
        (("independent_delivery_states", "kokuno_openai_field_identified"), True),
        (("cr001_invariants", "threshold_relaxation_this_increment"), True),
        (("cr001_invariants", "forcing_family_change_this_increment"), True),
    ],
)
def test_forbidden_truth_and_shortcut_mutations_fail_closed(
    tmp_path: Path, path: tuple[str, str], value: object
) -> None:
    payload = _payload()
    payload[path[0]][path[1]] = value
    errors = _audit_temp(tmp_path, payload)
    assert errors, f"mutation {path}={value!r} did not fail closed"


def test_reclassifying_construction_residual_as_source_fact_fails_closed(tmp_path: Path) -> None:
    payload = _payload()
    payload["role_separation"][
        "construction_closure_residual_classification"
    ] = "public_source_fact"
    errors = _audit_temp(tmp_path, payload)
    assert any("classification" in error for error in errors)


def test_relaxing_registered_momentum_gate_fails_closed(tmp_path: Path) -> None:
    payload = _payload()
    payload["cr001_invariants"]["thresholds"]["pde_residual_max"] = 0.1
    errors = _audit_temp(tmp_path, payload)
    assert any("CR001 thresholds" in error for error in errors)


def test_mechanics_witness_is_nonvacuous_and_information_losing() -> None:
    payload = _payload()
    witness = payload["mechanics_witness"]
    receipt = bookkeeping_projection_witness(witness["residual_a"], witness["residual_b"])
    assert receipt["same_mj_projection"] is True
    assert receipt["same_full_residual"] is False
    assert receipt["norm_a"] != receipt["norm_b"]


def test_mechanics_witness_is_not_source_or_pde_evidence() -> None:
    payload = _payload()
    witness = payload["mechanics_witness"]
    assert witness["classification"] == "autonomous_mechanics_only"
    assert payload["role_separation"]["construction_closure_residual_is_forcing"] is False
    assert (
        payload["role_separation"]["construction_closure_residual_is_heldout_ns_residual"]
        is False
    )
    assert payload["independent_delivery_states"]["kokuno_pde_validated"] is False


def test_canonical_eq45_delivery_stays_independent() -> None:
    payload = _payload()
    states = payload["independent_delivery_states"]
    assert states["canonical_eq45_velocity_export_ready"] is True
    assert states["kokuno_pulse_entry_route_velocity_export_ready"] is False
