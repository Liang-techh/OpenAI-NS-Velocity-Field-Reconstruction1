from __future__ import annotations

import json
from copy import deepcopy
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_eq45_supported_energy_failure_governance import (
    audit_documents,
    audit_repository,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def _documents() -> tuple[dict, dict, dict, dict]:
    return (
        _load("configs/eq45_supported_energy_failure_scope.json"),
        _load("configs/constraints.json"),
        _load("artifacts/constrained/eq45_supported_energy_audit.json"),
        _load("configs/delivery_state_contract.json"),
    )


def test_supported_energy_failure_contract_matches_checked_child() -> None:
    report = audit_repository(ROOT)

    assert report["candidate_sha256"] == (
        "2fdff812c22131d56eac1b7d6e455187ff3207500c6385aa9abf282a8e3d7b1d"
    )
    assert report["reference_energy"] == 0.9054304573980382
    assert report["reference_energy_abs_error"] == 0.09456954260196182
    assert report["reference_energy_gate"] == "fail"
    assert report["quadrature_gate"] == "pass"
    assert report["per_time_energy_range_gate"] == "pass"
    assert report["velocity_export_blocked_by_energy_failure"] is False
    assert report["implicit_renormalization_allowed"] is False
    assert report["renormalization_requires_new_candidate_identity"] is True
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("reference_energy_failure_blocks_velocity_export", True),
        ("implicit_export_time_amplitude_renormalization_allowed", True),
        ("amplitude_rescale_preserves_candidate_identity", True),
        ("amplitude_rescale_is_metadata_only", True),
        ("renormalized_child_may_inherit_parent_validation_without_recheck", True),
        ("renormalization_factor_must_be_serialized", False),
    ],
)
def test_policy_mutations_fail_closed(key: str, value: object) -> None:
    contract, constraints, energy_audit, delivery = _documents()
    mutated = deepcopy(contract)
    mutated["policy"][key] = value

    with pytest.raises(ValueError):
        audit_documents(mutated, constraints, energy_audit, delivery)


def test_reference_energy_failure_cannot_be_reclassified_as_pass() -> None:
    contract, constraints, energy_audit, delivery = _documents()
    mutated = deepcopy(contract)
    mutated["measured_supported_child_state"]["reference_energy_gate_pass"] = True
    mutated["measured_supported_child_state"]["cr001_energy_normalization_status"] = "pass"

    with pytest.raises(ValueError):
        audit_documents(mutated, constraints, energy_audit, delivery)


def test_reference_tolerance_cannot_be_relaxed_in_place() -> None:
    contract, constraints, energy_audit, delivery = _documents()
    mutated = deepcopy(contract)
    mutated["registered_energy_contract"]["reference_energy_abs_tolerance"] = 0.1

    with pytest.raises(ValueError):
        audit_documents(mutated, constraints, energy_audit, delivery)


def test_renormalized_child_revalidation_cannot_drop_divergence() -> None:
    contract, constraints, energy_audit, delivery = _documents()
    mutated = deepcopy(contract)
    mutated["renormalized_child_requires_revalidation"].remove("divergence_validation")

    with pytest.raises(ValueError):
        audit_documents(mutated, constraints, energy_audit, delivery)


def test_scientific_state_cannot_be_promoted_by_energy_governance() -> None:
    contract, constraints, energy_audit, delivery = _documents()
    mutated = deepcopy(contract)
    mutated["truth_boundary"]["pde_validated"] = True

    with pytest.raises(ValueError):
        audit_documents(mutated, constraints, energy_audit, delivery)
