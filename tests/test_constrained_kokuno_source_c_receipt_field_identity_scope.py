from __future__ import annotations

import copy
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_kokuno_source_c_receipt_field_identity_scope import (
    GovernanceError,
    audit,
    behavioral_identity_witness,
    load_constraints,
    load_contract,
    validate_contract,
    validate_cr001,
    validate_source_structure,
)


def _source_text() -> str:
    root = Path(__file__).resolve().parents[1]
    return (
        root
        / "src"
        / "openai_ns_reconstruction"
        / "kokuno_pa10_source_c_normalized_physical_center.py"
    ).read_text()


def test_live_source_c_receipt_field_identity_audit() -> None:
    receipt = audit()
    assert receipt["schema"] == "cr002-kokuno-source-c-receipt-field-identity-scope-v1"
    assert receipt["upstream_pr"] == 796
    assert receipt["source_structure"]["selected_source_C"] == 1000.0
    assert receipt["source_structure"]["sha256_is_report_receipt_sha256"] is True
    assert receipt["identity_state"]["field_semantics_callable"] is True
    assert receipt["identity_state"]["current_sha256_is_evidence_receipt_identity"] is True
    assert receipt["identity_state"]["stable_field_semantic_digest_exposed_by_796"] is False
    assert receipt["identity_state"]["receipt_sha256_is_candidate_sha256"] is False
    assert receipt["identity_state"]["pde_validated"] is False
    assert receipt["CR001_unchanged"] is True


def test_certificate_chunking_is_not_field_semantics() -> None:
    witness = behavioral_identity_witness()
    assert witness == {
        "same_physical_profile_configuration": True,
        "different_full_evidence_configuration": True,
        "same_values_on_frozen_probe": True,
        "same_derivatives_on_frozen_probe": True,
        "left_chunk_size": 1024,
        "right_chunk_size": 2048,
        "certificate_executed": False,
    }


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("identity_roles", "current_sha256_is_stable_field_semantic_identity"), True),
        (("identity_roles", "current_sha256_is_candidate_sha256"), True),
        (("identity_roles", "certificate_execution_parameters_are_field_parameters"), True),
        (("claim_boundaries", "receipt_sha256_is_field_identity"), True),
        (("claim_boundaries", "source_C_machine_self_certificate_implies_pde_validated"), True),
        (("truth_boundary", "velocity_export_ready_promoted"), True),
        (("truth_boundary", "openai_field_identified"), True),
    ],
)
def test_contract_rejects_identity_or_scientific_laundering(
    path: tuple[str, str], value: bool
) -> None:
    payload = copy.deepcopy(load_contract())
    payload[path[0]][path[1]] = value
    with pytest.raises(GovernanceError):
        validate_contract(payload)


def test_contract_rejects_moving_source_c_out_of_field_identity() -> None:
    payload = copy.deepcopy(load_contract())
    payload["identity_roles"]["future_field_identity_must_bind_selected_C"] = False
    with pytest.raises(GovernanceError):
        validate_contract(payload)


def test_source_audit_rejects_certificate_knob_in_field_evaluation() -> None:
    source = _source_text()
    mutated = source.replace(
        "return self.physical_profiles.values(X, eta)",
        "return self.physical_profiles.values(X, eta) if self.chunk_size else {}",
        1,
    )
    assert mutated != source
    with pytest.raises(GovernanceError):
        validate_source_structure(mutated)


def test_source_audit_rejects_sha256_decoupled_from_receipt_without_new_governance() -> None:
    source = _source_text()
    mutated = source.replace(
        'return str(self.report()["receipt_sha256"])',
        'return str(self.configuration()["physical_profile_configuration"])',
        1,
    )
    assert mutated != source
    with pytest.raises(GovernanceError):
        validate_source_structure(mutated)


def test_cr001_threshold_relaxation_fails_closed() -> None:
    constraints = copy.deepcopy(load_constraints())
    constraints["validation"]["thresholds"]["pde_residual_max"] = 0.002
    with pytest.raises(GovernanceError):
        validate_cr001(constraints, load_contract())


def test_cr001_free_force_escape_fails_closed() -> None:
    constraints = copy.deepcopy(load_constraints())
    constraints["forcing"]["restriction"] = "Any pointwise residual-defined force may be fitted."
    with pytest.raises(GovernanceError):
        validate_cr001(constraints, load_contract())


def test_cr001_amplitude_collapse_escape_fails_closed() -> None:
    constraints = copy.deepcopy(load_constraints())
    constraints["nontriviality"]["enforcement"] = "allow collapsed candidates"
    with pytest.raises(GovernanceError):
        validate_cr001(constraints, load_contract())
