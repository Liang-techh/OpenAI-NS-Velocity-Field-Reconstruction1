from __future__ import annotations

import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.audit_kokuno_rf40_lambda_turn_source_provenance_semantic_scope import (
    SOURCE_PROVENANCE,
    assert_scope,
    audit_scope,
    semantic_scope_witness,
)


CONTRACT_PATH = (
    Path(__file__).resolve().parents[1]
    / "configs"
    / "kokuno_rf40_lambda_turn_source_provenance_semantic_scope.json"
)


def _contract():
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def _assert_rejected(mutated):
    errors = audit_scope(mutated)
    assert errors
    with pytest.raises(ValueError):
        assert_scope(mutated)


def test_clean_contract_is_admitted():
    contract = _contract()
    assert audit_scope(contract) == []
    assert_scope(contract)


def test_four_way_provenance_stays_explicit():
    provenance = _contract()["four_way_provenance"]
    assert set(provenance) == {
        "user_requirements",
        "public_source_facts",
        "autonomous_repository_facts",
        "pending_or_unknown",
    }
    assert all(provenance[key] for key in provenance)


def test_exact_a1_lambda_turn_identity_is_pinned():
    a1 = _contract()["upstream_scope"]["agent1_998"]
    assert a1["head"] == "43b295444b1e9558222d757cb551e04385eddcb5"
    assert a1["source_blob"] == "07743c30978360e305a8863e05e7ed322d818b32"
    assert a1["semantic_payload_fields_observed"] == [
        "schema",
        "source_formulas",
        "numerical_realization",
        "truth_boundary",
        "configuration",
    ]


def test_report_and_ci_source_provenance_do_not_become_semantic_binding():
    locked = _contract()["machine_locked_distinctions"]
    assert locked["a1_998_report_carries_external_source_provenance"] is True
    assert locked["a1_998_dedicated_ci_pins_external_source_blob"] is True
    assert locked["a1_998_semantic_sha_binds_external_source_provenance"] is False
    assert locked["ci_source_blob_pin_implies_serialized_semantic_source_binding"] is False
    assert locked["report_source_block_implies_serialized_semantic_source_binding"] is False
    assert locked["source_provenance_semantic_identity_closed"] is False


def test_source_provenance_mutation_is_rejected():
    mutated = _contract()
    mutated["corrected_public_source_provenance"]["commit"] = "f" * 40
    _assert_rejected(mutated)


def test_false_semantic_source_binding_promotion_is_rejected():
    mutated = _contract()
    mutated["machine_locked_distinctions"][
        "a1_998_semantic_sha_binds_external_source_provenance"
    ] = True
    _assert_rejected(mutated)


def test_false_source_identity_closure_is_rejected():
    mutated = _contract()
    mutated["machine_locked_distinctions"]["source_provenance_semantic_identity_closed"] = True
    _assert_rejected(mutated)


@pytest.mark.parametrize(
    "key",
    [
        "kokuno_lambda_turn_velocity_export_ready",
        "kokuno_pde_validated",
        "kokuno_paper_exact",
        "kokuno_openai_field_identified",
        "kokuno_blowup_proved",
        "residual_defined_free_force_allowed",
        "candidate_collapse_allowed",
        "post_hoc_threshold_relaxation_allowed",
    ],
)
def test_forbidden_promotions_fail_closed(key):
    mutated = _contract()
    mutated["forbidden_promotions"][key] = True
    _assert_rejected(mutated)


def test_lambda_turn_composite_cannot_be_inferred_from_a1_callable():
    mutated = _contract()
    mutated["upstream_scope"]["agent5_1002"]["lambda_turn_composite_materialized"] = True
    _assert_rejected(mutated)


def test_cr001_threshold_drift_is_rejected():
    mutated = _contract()
    mutated["cr001_snapshot"]["momentum_max_threshold"] = 0.01
    _assert_rejected(mutated)


def test_free_force_drift_is_rejected_even_if_thresholds_unchanged():
    mutated = _contract()
    mutated["cr001_snapshot"]["residual_defined_free_force_allowed"] = True
    _assert_rejected(mutated)


def test_eq45_delivery_cannot_be_downgraded_by_kokuno_incompleteness():
    mutated = _contract()
    mutated["canonical_delivery_independent_state"]["velocity_export_ready"] = False
    _assert_rejected(mutated)


def test_semantic_scope_witness_detects_only_bound_source_mutation():
    witness = semantic_scope_witness()
    assert witness["classification"] == "autonomous_mechanics_only"
    assert witness["not_a_public_source_fact"] is True
    assert witness["not_candidate_numerical_evidence"] is True
    assert witness["unbound_digest_a"] == witness["unbound_digest_b"]
    assert witness["unbound_digest_detects_source_mutation"] is False
    assert witness["bound_digest_a"] != witness["bound_digest_b"]
    assert witness["bound_digest_detects_source_mutation"] is True


def test_future_source_identity_closure_requires_full_tuple_and_replay():
    requirements = _contract()["future_promotion_requirements"][
        "source_provenance_semantic_identity_closed"
    ]
    joined = " ".join(requirements)
    for token in ("repository", "commit", "source path", "source blob"):
        assert token in joined
    assert "recompute" in joined
    assert "fail closed" in joined
    assert "paper exactness" in joined
    assert "OpenAI-field identity" in joined


def test_original_openai_paper_exactness_remains_unverified_here():
    assert (
        _contract()["corrected_public_source_provenance"][
            "original_openai_paper_exact_correspondence_independently_verified_here"
        ]
        is False
    )


def test_source_provenance_block_matches_frozen_expected_tuple():
    assert _contract()["corrected_public_source_provenance"] == SOURCE_PROVENANCE
