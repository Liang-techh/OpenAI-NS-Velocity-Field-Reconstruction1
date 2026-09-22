from __future__ import annotations

import copy
import json

import pytest

from openai_ns_reconstruction.kokuno_release1_decimal_precision_governance import (
    CONSTRAINTS_PATH,
    IMPLEMENTATION_PATH,
    GovernanceError,
    audit,
    audit_payload,
    load_contract,
    mechanics_only_precision_witness,
)


def _inputs():
    contract = load_contract()
    impl_bytes = IMPLEMENTATION_PATH.read_bytes()
    constraints_bytes = CONSTRAINTS_PATH.read_bytes()
    return (
        contract,
        impl_bytes.decode("utf-8"),
        impl_bytes,
        json.loads(constraints_bytes.decode("utf-8")),
        constraints_bytes,
    )


def test_exact_a1_release1_precision_scope_passes() -> None:
    report = audit()
    assert report["status"] == "pass"
    assert report["decimal_output_encoding_materialized"] is True
    assert report["end_to_end_96digit_arithmetic_materialized"] is False
    assert report["global_velocity_export_promoted"] is False
    assert report["pde_promoted"] is False


def test_decimal_embedding_does_not_recover_sub_ulp_distinction() -> None:
    witness = mechanics_only_precision_witness()
    assert witness["scope"] == "autonomous_representation_mechanics_only"
    assert witness["exact_delta"] == "1E-30"
    assert witness["binary64_values_equal"] is True
    assert witness["decimal_embeddings_equal"] is True
    assert witness["embedding_error_nonzero"] is True


@pytest.mark.parametrize(
    "key",
    [
        "release1_end_to_end_96digit_arithmetic_materialized",
        "parent_1179_sub_epsilon_survival_evidence_transfers_to_new_release1_arithmetic",
        "binary64_to_decimal_embedding_recovers_discarded_sub_ulp_information",
        "source_exact_numerical_precision_identified",
        "unified_global_cartesian_velocity_export_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ],
)
def test_forbidden_precision_or_scientific_promotions_fail_closed(key: str) -> None:
    contract, text, impl_bytes, constraints, constraints_bytes = _inputs()
    mutated = copy.deepcopy(contract)
    mutated["machine_state"][key] = True
    with pytest.raises(GovernanceError):
        audit_payload(mutated, text, impl_bytes, constraints, constraints_bytes)


def test_decimal_output_must_not_be_relabelled_as_end_to_end_decimal_arithmetic() -> None:
    contract, text, impl_bytes, constraints, constraints_bytes = _inputs()
    mutated = copy.deepcopy(contract)
    mutated["representation_identity"]["new_release1_cartesian_assembly_arithmetic"] = "decimal96"
    with pytest.raises(GovernanceError):
        audit_payload(mutated, text, impl_bytes, constraints, constraints_bytes)


def test_binary64_before_decimal_implementation_mechanics_are_pinned() -> None:
    contract, text, impl_bytes, constraints, constraints_bytes = _inputs()
    mutated_text = text.replace(
        "out[idx] = _embed_decimal(out_float[idx])",
        "out[idx] = out_float[idx]",
        1,
    )
    with pytest.raises(GovernanceError):
        audit_payload(contract, mutated_text, impl_bytes, constraints, constraints_bytes)


def test_source_classification_cannot_launder_autonomous_precision_as_public_fact() -> None:
    contract, text, impl_bytes, constraints, constraints_bytes = _inputs()
    mutated = copy.deepcopy(contract)
    autonomous_items = list(mutated["source_classification"]["autonomous_design"])
    mutated["source_classification"]["public_source_fact"].extend(autonomous_items)
    mutated["source_classification"]["autonomous_design"] = []
    with pytest.raises(GovernanceError):
        audit_payload(mutated, text, impl_bytes, constraints, constraints_bytes)


def test_cr001_threshold_relaxation_fails_closed() -> None:
    contract, text, impl_bytes, constraints, constraints_bytes = _inputs()
    mutated = copy.deepcopy(contract)
    mutated["cr001_unchanged"]["pde_residual_max"] = 0.002
    with pytest.raises(GovernanceError):
        audit_payload(mutated, text, impl_bytes, constraints, constraints_bytes)
