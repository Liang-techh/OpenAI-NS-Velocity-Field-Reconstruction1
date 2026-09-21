from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_x4_mean_audit_registration_scope import (
    CONTRACT_REL,
    audit_contract,
    audit_repository,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def _contract() -> dict:
    return json.loads((REPO_ROOT / CONTRACT_REL).read_text(encoding="utf-8"))


def test_contract_audits_clean() -> None:
    assert audit_contract(_contract()) == []


def test_repository_binding_audits_clean() -> None:
    assert audit_repository(REPO_ROOT) == []


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("current_repository_evidence", "agent5_1030_consumes_agent3_1028", True),
        ("current_repository_evidence", "agent5_1030_consumes_agent4_1029", True),
        ("current_repository_evidence", "agent5_1030_registers_agent4_1029", True),
        ("required_false", "x4_compact_radial_stress_materialized", True),
        ("required_false", "x4_radial_force_partial_z_sigma1_materialized", True),
        ("required_false", "x4_mean_audit_means_x4_stress_or_force_audit", True),
        ("required_false", "x4_mean_audit_means_complete_ns_defect", True),
        ("required_false", "x4_mean_audit_means_authorized_correction_target", True),
        ("required_false", "cartesian_correction_velocity_materialized", True),
        ("required_false", "complete_candidate_api_materialized_on_kokuno_route", True),
        ("required_false", "pde_validated", True),
        ("required_false", "visual_correspondence_verified", True),
        ("required_false", "paper_exact", True),
        ("required_false", "openai_field_identified", True),
    ],
)
def test_forbidden_promotions_fail_closed(section: str, key: str, value: object) -> None:
    mutated = copy.deepcopy(_contract())
    mutated[section][key] = value
    errors = audit_contract(mutated)
    assert errors, f"mutation {section}.{key} was not rejected"


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("nu", 0.02),
        ("forcing_mode", "residual_defined_free_force"),
        ("reference_energy", 0.0),
        ("reference_energy_abs_tolerance", 0.01),
        ("validation_seed", 1),
        ("held_out_points", 1024),
        ("derivative_steps", [0.02, 0.01]),
        ("quadrature_orders_per_axis", [24, 48]),
        ("momentum_max_l2_threshold", 0.002),
        ("divergence_max_l2_threshold", 0.00002),
        ("residual_defined_free_force_forbidden", False),
        ("candidate_collapse_forbidden", False),
        ("post_hoc_threshold_relaxation_forbidden", False),
    ],
)
def test_cr001_mutations_fail_closed(key: str, value: object) -> None:
    mutated = copy.deepcopy(_contract())
    mutated["cr001_snapshot"][key] = value
    assert audit_contract(mutated), f"CR001 mutation {key} was not rejected"


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
def test_canonical_delivery_state_laundering_fails_closed(key: str, value: object) -> None:
    mutated = copy.deepcopy(_contract())
    mutated["canonical_delivery_independence"][key] = value
    assert audit_contract(mutated), f"canonical state mutation {key} was not rejected"


def test_x4_mean_audit_identity_mutation_fails_closed() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["current_repository_evidence"]["agent4_x4_mean_audit_head"] = "0" * 40
    assert audit_contract(mutated)


def test_a5_checkpoint_identity_mutation_fails_closed() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["base_a5"]["head"] = "f" * 40
    assert audit_contract(mutated)


def test_public_source_laundering_fails_closed() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["provenance_classes"]["public_source_fact"] = [
        "X4 compact radial stress is a public OpenAI fact"
    ]
    assert audit_contract(mutated)


def test_scoped_audit_cannot_be_rewritten_as_stress_force_validation() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["evidence_semantics"]["representation_rule"] = (
        "Independent X4 nonlinear-mean validation also validates compact stress and radial force."
    )
    assert audit_contract(mutated)


def test_future_force_promotion_requires_compact_stress_and_complete_defect_firewall() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["promotion_requirements"]["future_x4_radial_force_materialization"] = [
        "consume exact A3 #1028"
    ]
    errors = audit_contract(mutated)
    assert any("compact-stress" in error for error in errors)
    assert any("complete defect" in error for error in errors)


def test_queued_ci_cannot_be_rewritten_as_scientific_pass() -> None:
    mutated = copy.deepcopy(_contract())
    mutated["ci_observation_at_creation"]["scientific_pass_claimed_from_these_observations"] = True
    assert audit_contract(mutated)
