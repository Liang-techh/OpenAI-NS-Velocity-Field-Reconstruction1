from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_pre_i1_autonomous_modulation_scope import (
    CONTRACT_ID,
    audit_contract,
    audit_repository,
)

ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "configs" / "kokuno_pre_i1_autonomous_modulation_scope.json"


def _contract() -> dict:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def _identity_sha(payload: dict) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    return hashlib.sha256(raw).hexdigest()


def test_contract_id_and_baseline_are_clean() -> None:
    contract = _contract()
    assert contract["contract_id"] == CONTRACT_ID
    assert audit_contract(contract) == []


def test_repository_audit_passes_on_exact_ancestry() -> None:
    assert audit_repository(ROOT) == []


@pytest.mark.parametrize(
    "key",
    [
        "source_hidden_loop_parameters_recovered",
        "source_admissible_loop_reconstructed",
        "autonomous_modulation_equals_hidden_source_loop",
        "a2_1010_consumes_a1_1035_modulated_identity",
        "a4_1012_validates_a1_1035_modulated_velocity",
        "a3_1037_stress_applies_to_a1_1035_modulated_identity",
        "a4_1038_validates_a1_1035_modulated_identity",
        "leading_plus_oscillatory_modulated_pre_i1_materialized",
        "modulated_pre_i1_independent_cartesian_audit_available",
        "modulated_pre_i1_correction_mean_stress_force_materialized",
        "outer_global_leading_velocity_materialized",
        "kokuno_velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
    ],
)
def test_forbidden_promotions_fail_closed(key: str) -> None:
    mutated = copy.deepcopy(_contract())
    mutated["required_false"][key] = True
    assert audit_contract(mutated)


@pytest.mark.parametrize(
    "key",
    [
        "a1_1035_pre_i1_leading_cartesian_velocity_callable",
        "a1_1035_configuration_serializable",
        "a1_1035_semantic_identity_available",
        "a1_1035_repository_autonomous_modulation_bound",
        "a1_1035_prefix_incompressibility_memory_carried",
        "a1_1035_fails_closed_at_reserved_i1_boundary",
        "corrected_public_source_provenance_pinned",
        "existing_unmodulated_a2_a3_a4_evidence_remains_scoped_to_its_own_identity",
    ],
)
def test_required_scope_facts_cannot_be_erased(key: str) -> None:
    mutated = copy.deepcopy(_contract())
    mutated["allowed_true"][key] = False
    assert audit_contract(mutated)


def test_autonomous_parameters_cannot_be_relabelled_as_source_recovery() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["source_classification"]["repository_autonomous_realization"][
        "classification"
    ] = "public_source_parameter_recovery"
    assert audit_contract(mutated)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("a2_1010_consumed_leading_head", "ddf71bab2a2038bab2fba5b854ad713c6b6c7e95"),
        ("a3_1037_head", "0" * 40),
        ("a4_1038_head", "f" * 40),
    ],
)
def test_lineage_evidence_drift_is_rejected(key: str, value: str) -> None:
    mutated = copy.deepcopy(_contract())
    mutated["lineage_evidence"][key] = value
    assert audit_contract(mutated)


def test_domain_overlap_does_not_equal_candidate_identity_mechanics_witness() -> None:
    # Autonomous mechanics-only witness: identical declared coordinate coverage is
    # intentionally insufficient for identity equality when the leading lineage changes.
    common = {
        "domain": "same pre-I1 coordinate interval",
        "velocity_api": "velocity(x,y,z,t)",
        "time_interval": [0.25, 0.75],
    }
    unmodulated = {**common, "leading_head": "2c76ebdc41d6c566f43a1305034ba2ff9dce410b"}
    modulated = {**common, "leading_head": "ddf71bab2a2038bab2fba5b854ad713c6b6c7e95"}
    assert _identity_sha(unmodulated) != _identity_sha(modulated)


def test_source_exactness_requires_independent_evidence() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["promotion_requirements"]["future_source_exactness"] = [
        "reuse the autonomous loop and call it source exact"
    ]
    assert audit_contract(mutated)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("nu", 0.02),
        ("forcing_mode", "free_residual_force"),
        ("reference_energy", 0.0),
        ("momentum_max_l2_threshold", 0.01),
        ("divergence_max_l2_threshold", 0.001),
        ("residual_defined_free_force_forbidden", False),
        ("candidate_collapse_forbidden", False),
        ("post_hoc_threshold_relaxation_forbidden", False),
    ],
)
def test_cr001_snapshot_mutations_are_rejected(key: str, value: object) -> None:
    mutated = copy.deepcopy(_contract())
    mutated["cr001_snapshot"][key] = value
    assert audit_contract(mutated)


def test_canonical_eq45_delivery_cannot_be_downgraded_by_kokuno_incompleteness() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["canonical_delivery_independence"]["velocity_export_ready"] = False
    assert audit_contract(mutated)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("visualization_ready", True),
        ("visual_correspondence_verified", True),
        ("pde_validated", True),
        ("paper_exact", True),
        ("openai_field_identified", True),
    ],
)
def test_canonical_truth_promotions_fail_closed(key: str, value: bool) -> None:
    mutated = copy.deepcopy(_contract())
    mutated["canonical_delivery_independence"][key] = value
    assert audit_contract(mutated)


def test_queued_ci_cannot_be_rewritten_as_scientific_pass() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["ci_observation_at_creation"]["a1_1035_observed_status"] = "success"
    mutated["ci_observation_at_creation"]["scientific_pass_claimed_from_these_observations"] = True
    assert audit_contract(mutated)


def test_live_integration_ref_is_locked_separately_from_stale_doc_snapshots() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["active_integration_snapshot"]["live_head"] = "ca37d25d19b97d893030194ebd6364160ae4355e"
    assert audit_contract(mutated)


def test_unreviewed_truth_keys_fail_closed() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["required_false"]["visual_similarity_means_pde_valid"] = False
    assert audit_contract(mutated)
