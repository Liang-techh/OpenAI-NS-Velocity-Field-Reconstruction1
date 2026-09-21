from __future__ import annotations

from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_i2_float64_observability_scope import (
    EXPECTED_CONSTRAINTS_BLOB,
    EXPECTED_PARENT_HEAD,
    EXPECTED_UPSTREAM_HEAD,
    audit_contract,
    git_blob_sha,
    load_contract,
    repository_live_violations,
    subulp_decimal_witness,
)


def _set_path(payload, dotted_path, value):
    node = payload
    parts = dotted_path.split(".")
    for part in parts[:-1]:
        node = node[part]
    node[parts[-1]] = value


def test_contract_baseline_is_fail_closed_clean():
    contract = load_contract()
    assert audit_contract(contract) == []


def test_subulp_witness_is_nonvacuous_and_mechanics_only():
    witness = subulp_decimal_witness()
    assert witness["classification"] == (
        "autonomous_mechanics_only_not_kokuno_or_openai_candidate_data"
    )
    assert witness["high_precision_delta_nonzero"] is True
    assert witness["float64_values_equal"] is True
    assert witness["semantic_payload_sha_changed"] is True
    assert float(witness["synthetic_delta"]) < witness["float64_ulp_at_base"]


@pytest.mark.parametrize(
    ("path", "mutated"),
    [
        ("schema_version", 2),
        ("contract_id", "laundered"),
        ("scope", "scientific_admission"),
        ("target.upstream_pr", 1065),
        ("target.upstream_head", "0" * 40),
        ("target.parent_i1_pr", 1043),
        ("target.parent_i1_head", "1" * 40),
        (
            "target.module_path",
            "src/openai_ns_reconstruction/not_the_i2_candidate.py",
        ),
        ("representation_truth.high_precision_i2_correction_receipt_nonzero", False),
        ("representation_truth.decimal_retained_until_float64_api_boundary", False),
        ("representation_truth.public_float64_velocity_callable_through_i2", False),
        ("representation_truth.public_float64_delta_from_i1_parent_verified", True),
        ("representation_truth.public_float64_delta_probe_protocol_registered", True),
        ("representation_truth.public_float64_visual_effect_verified", True),
        ("representation_truth.high_precision_nonzero_implies_public_float64_delta", True),
        ("representation_truth.semantic_identity_progress_implies_public_float64_delta", True),
        ("representation_truth.callable_through_i2_implies_visual_correspondence", True),
        ("cr001_freeze.nu", 0.02),
        ("cr001_freeze.momentum_max", 0.01),
        ("cr001_freeze.divergence_max", 0.001),
        ("cr001_freeze.validation_seed", 7),
        ("cr001_freeze.held_out_points", 256),
        ("cr001_freeze.residual_defined_free_force_forbidden", False),
        ("cr001_freeze.candidate_collapse_forbidden", False),
        ("cr001_freeze.posthoc_threshold_relaxation_forbidden", False),
        ("independent_project_states.canonical_eq45_velocity_export_ready", False),
        ("independent_project_states.canonical_eq45_visualization_ready", True),
        ("independent_project_states.canonical_eq45_visual_correspondence_verified", True),
        ("independent_project_states.canonical_eq45_pde_validated", True),
        ("independent_project_states.kokuno_current_route_velocity_export_ready", True),
        ("independent_project_states.kokuno_current_route_visual_correspondence_verified", True),
        ("independent_project_states.kokuno_current_route_pde_validated", True),
        ("independent_project_states.kokuno_current_route_paper_exact", True),
        ("independent_project_states.kokuno_current_route_openai_field_identified", True),
    ],
)
def test_protected_scope_mutations_fail_closed(path, mutated):
    contract = load_contract()
    _set_path(contract, path, mutated)
    assert audit_contract(contract), path


def test_missing_provenance_class_fails_closed():
    contract = load_contract()
    del contract["provenance_classes"]["pending_or_unknown"]
    violations = audit_contract(contract)
    assert any("four governed classes" in item for item in violations)


def test_empty_provenance_class_fails_closed():
    contract = load_contract()
    contract["provenance_classes"]["public_source_facts"] = []
    violations = audit_contract(contract)
    assert any("public_source_facts" in item for item in violations)


def test_missing_forbidden_implication_fails_closed():
    contract = load_contract()
    contract["forbidden_implications"].remove(
        "high_precision_nonzero_correction -> nonzero_public_float64_velocity_delta"
    )
    violations = audit_contract(contract)
    assert any("forbidden_implications missing" in item for item in violations)


def test_public_float64_promotion_protocol_is_explicit_and_negative_evidence_safe():
    contract = load_contract()
    requirements = " ".join(
        contract["promotion_requirements"][
            "public_float64_delta_from_i1_parent_verified"
        ]
    ).lower()
    for token in ("save", "reload", "#1061", "#1051", "float64", "ulp", "nonzero", "negative"):
        assert token in requirements
    assert "do not change" in requirements


def test_exact_identity_pins_are_current():
    contract = load_contract()
    assert contract["target"]["upstream_head"] == EXPECTED_UPSTREAM_HEAD
    assert contract["target"]["parent_i1_head"] == EXPECTED_PARENT_HEAD
    assert contract["cr001_freeze"]["constraints_blob_sha"] == EXPECTED_CONSTRAINTS_BLOB


def test_repository_live_replay():
    assert repository_live_violations() == []


def test_canonical_constraints_git_blob_replay():
    root = Path(__file__).resolve().parents[1]
    data = (root / "configs/constraints.json").read_bytes()
    assert git_blob_sha(data) == EXPECTED_CONSTRAINTS_BLOB
