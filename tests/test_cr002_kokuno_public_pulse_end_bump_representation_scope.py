from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_public_pulse_end_bump_representation_scope import (
    CONTRACT_PATH,
    audit_scope,
    autonomous_bump_identity_witness,
)


ROOT = Path(__file__).resolve().parents[1]


def _contract() -> dict:
    return json.loads((ROOT / CONTRACT_PATH).read_text(encoding="utf-8"))


def test_repository_live_scope_audit_passes() -> None:
    assert audit_scope(ROOT) == []


def test_autonomous_bump_identity_witness_is_nonvacuous_and_closes_each_system() -> None:
    witness = autonomous_bump_identity_witness()
    assert witness["matrix_linf_delta"] > 1.0e-4
    assert witness["coefficient_linf_delta"] > 1.0e-3
    assert witness["shape1_closure_max"] <= 1.0e-12
    assert witness["shape2_closure_max"] <= 1.0e-12


@pytest.mark.parametrize(
    ("path", "value"),
    [
        (("source_vs_realization", "repository_autonomous_pointwise_bump_materialized"), False),
        (("source_vs_realization", "source_exact_pointwise_bump_recovered"), True),
        (("source_vs_realization", "current_cartesian_end_compensation_composed"), True),
        (("evidence_transfer", "same_public_centers_width_imply_same_pointwise_bump_identity"), True),
        (("evidence_transfer", "same_public_algebra_implies_same_numerical_c1_c2_across_bump_realizations"), True),
        (("evidence_transfer", "coefficient_receipts_transfer_across_pointwise_bump_replacement"), True),
        (("evidence_transfer", "fresh_semantic_identity_required_if_pointwise_bump_changes"), False),
        (("evidence_transfer", "fresh_save_load_and_numerical_replay_required_if_pointwise_bump_changes"), False),
        (("evidence_transfer", "a2_pr1117_receipts_cover_pr1116_end_compensation"), True),
        (("truth_states", "source_exact_bump_shape_recovered"), True),
        (("truth_states", "current_lineage_J_entry_materialized"), True),
        (("truth_states", "current_cartesian_end_compensation_composed"), True),
        (("truth_states", "kokuno_velocity_export_ready"), True),
        (("truth_states", "visual_correspondence_verified"), True),
        (("truth_states", "pde_validated"), True),
        (("truth_states", "paper_exact"), True),
        (("truth_states", "openai_field_identified"), True),
        (("canonical_cr001", "momentum_max"), 0.01),
        (("canonical_cr001", "divergence_max"), 0.001),
        (("canonical_cr001", "residual_defined_free_forcing_forbidden"), False),
        (("canonical_delivery", "velocity_export_ready"), False),
        (("canonical_delivery", "pde_validated"), True),
    ],
)
def test_forbidden_scope_mutations_fail_closed(path: tuple[str, str], value: object) -> None:
    payload = copy.deepcopy(_contract())
    section, key = path
    payload[section][key] = value
    violations = audit_scope(ROOT, payload)
    assert violations, f"mutation {section}.{key} unexpectedly passed"


def test_wrong_exact_source_identity_fails_closed() -> None:
    payload = copy.deepcopy(_contract())
    payload["exact_stack"]["source_module_git_blob_sha1"] = "0" * 40
    assert audit_scope(ROOT, payload)


def test_a2_consumption_cannot_be_promoted_without_new_lineage() -> None:
    payload = copy.deepcopy(_contract())
    payload["exact_stack"]["current_a2_consumes_base_pr1116"] = True
    assert audit_scope(ROOT, payload)


@pytest.mark.parametrize(
    "bucket",
    ["user_requirements", "public_source_facts", "autonomous_design", "pending_unknown"],
)
def test_four_way_provenance_buckets_are_required(bucket: str) -> None:
    payload = copy.deepcopy(_contract())
    payload["four_way_provenance"][bucket] = []
    assert audit_scope(ROOT, payload)


def test_required_forbidden_implications_cannot_be_removed() -> None:
    payload = copy.deepcopy(_contract())
    payload["promotion_rule"]["forbidden_implications"] = [
        item
        for item in payload["promotion_rule"]["forbidden_implications"]
        if item != "same_centers_width -> transferable_c1_c2_receipt"
    ]
    assert audit_scope(ROOT, payload)
