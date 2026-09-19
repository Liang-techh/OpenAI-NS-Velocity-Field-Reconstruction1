from __future__ import annotations

import copy

import pytest

from openai_ns_reconstruction.constrained_st052m_whole_candidate_runtime_binding_governance import (
    GovernanceError,
    audit_contract,
    audit_default_contract,
    load_contract,
)


def _contract():
    return load_contract()


def test_default_contract_passes():
    audited = audit_default_contract()
    assert audited["snapshot"]["whole_candidate_identity_sha256"] == (
        "333615d89423a22995d9356b29b604d1c241c1b40a631bd66892c2701c536064"
    )
    assert audited["delivery_states"]["experimental_st052_velocity_export_ready"] is False


def test_rejects_source_head_authentication_laundering():
    obj = copy.deepcopy(_contract())
    obj["external_runtime_binding"]["source_git_head_verified_on_load"] = True
    with pytest.raises(GovernanceError, match="source_git_head_verified_on_load"):
        audit_contract(obj)


def test_rejects_runtime_tree_binding_laundering():
    obj = copy.deepcopy(_contract())
    obj["external_runtime_binding"]["source_git_tree_sha_bound_by_bundle"] = True
    with pytest.raises(GovernanceError, match="source_git_tree_sha_bound_by_bundle"):
        audit_contract(obj)


def test_rejects_probe_binding_as_whole_domain_identity():
    obj = copy.deepcopy(_contract())
    obj["external_runtime_binding"]["parent_probe_binding_is_whole_domain_function_identity"] = True
    with pytest.raises(GovernanceError, match="whole_domain"):
        audit_contract(obj)


def test_rejects_bundle_hash_as_complete_function_identity():
    obj = copy.deepcopy(_contract())
    obj["identity_scope"]["whole_candidate_identity_sha256_is_complete_executable_function_identity"] = True
    with pytest.raises(GovernanceError, match="complete_executable"):
        audit_contract(obj)


def test_rejects_exact_ci_parity_generalization():
    obj = copy.deepcopy(_contract())
    obj["identity_scope"]["dedicated_ci_parity_generalizes_to_arbitrary_exact_source_root_argument"] = True
    with pytest.raises(GovernanceError, match="arbitrary_exact_source_root"):
        audit_contract(obj)


def test_rejects_premature_st052_export_promotion():
    obj = copy.deepcopy(_contract())
    obj["delivery_states"]["experimental_st052_velocity_export_ready"] = True
    with pytest.raises(GovernanceError, match="experimental_st052_velocity_export_ready"):
        audit_contract(obj)


def test_rejects_pde_or_visual_truth_promotion():
    for key in ("pde_validated", "visual_correspondence_verified", "paper_exact", "openai_field_identified"):
        obj = copy.deepcopy(_contract())
        obj["delivery_states"][key] = True
        with pytest.raises(GovernanceError, match=key):
            audit_contract(obj)


def test_rejects_runtime_closure_prerequisite_removal():
    obj = copy.deepcopy(_contract())
    obj["promotion_requirements"]["loader_must_fail_on_runtime_identity_mismatch"] = False
    with pytest.raises(GovernanceError, match="loader_must_fail"):
        audit_contract(obj)


def test_rejects_cr001_threshold_or_seed_drift():
    for key, value in (("momentum_max", 0.01), ("validation_seed", 9175291)):
        obj = copy.deepcopy(_contract())
        obj["cr001_frozen"][key] = value
        with pytest.raises(GovernanceError, match="CR001 drift"):
            audit_contract(obj)


def test_rejects_source_classification_laundering():
    obj = copy.deepcopy(_contract())
    obj["source_classification"]["openai_hidden_velocity_coefficients"] = "public_source_fact"
    with pytest.raises(GovernanceError, match="hidden-field"):
        audit_contract(obj)
