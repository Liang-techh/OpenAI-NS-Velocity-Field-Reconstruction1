from __future__ import annotations

from openai_ns_reconstruction.constrained_st052m_transitive_runtime_graph_governance import (
    audit,
    audit_contract,
    load_contract,
    mutated,
    observe_upstream_loader,
)


def test_registered_contract_matches_upstream_runtime_state() -> None:
    report = audit()
    assert report["passed"], report["errors"]
    truth = report["truth_boundary"]
    assert truth["exact_source_checkout_identity_closed"] is True
    assert truth["direct_replay_import_cache_isolation_verified"] is True
    assert truth["transitive_runtime_import_cache_isolation_verified"] is False
    assert truth["whole_executable_runtime_identity_closed"] is False
    assert truth["experimental_st052_velocity_export_ready"] is False


def test_live_upstream_has_direct_guard_but_not_transitive_guard() -> None:
    observed = observe_upstream_loader()
    assert observed["authenticates_source_before_import"] is True
    assert observed["identity_guard_checks_git_head"] is True
    assert observed["identity_guard_checks_git_tree"] is True
    assert observed["evicts_direct_replay_module"] is True
    assert observed["uses_bare_replay_import"] is True
    assert observed["transitive_cache_eviction_present"] is False
    assert observed["transitive_origin_check_present"] is False


def test_transitive_cache_promotion_fails_closed() -> None:
    contract = load_contract()
    promoted = mutated(
        contract,
        "allowed_states",
        "transitive_runtime_import_cache_isolation_verified",
        value=True,
    )
    report = audit_contract(promoted)
    assert report["passed"] is False
    assert any("transitive_runtime_import_cache_isolation_verified" in e for e in report["errors"])


def test_whole_runtime_identity_promotion_fails_closed() -> None:
    contract = load_contract()
    promoted = mutated(
        contract,
        "allowed_states",
        "whole_executable_runtime_identity_closed",
        value=True,
    )
    report = audit_contract(promoted)
    assert report["passed"] is False
    assert any("whole_executable_runtime_identity_closed" in e for e in report["errors"])


def test_velocity_export_promotion_fails_closed() -> None:
    contract = load_contract()
    promoted = mutated(
        contract,
        "allowed_states",
        "experimental_st052_velocity_export_ready",
        value=True,
    )
    report = audit_contract(promoted)
    assert report["passed"] is False
    assert any("experimental_st052_velocity_export_ready" in e for e in report["errors"])


def test_transitive_graph_contract_cannot_claim_origin_checks() -> None:
    contract = load_contract()
    promoted = mutated(
        contract,
        "runtime_graph_scope",
        "transitive_module_origin_checks_verified",
        value=True,
    )
    report = audit_contract(promoted)
    assert report["passed"] is False
    assert any("transitive origin verification" in e for e in report["errors"])


def test_cr001_threshold_drift_is_rejected() -> None:
    contract = load_contract()
    relaxed = mutated(
        contract,
        "cr001_nonmutation",
        "pde_residual_max",
        value=0.01,
    )
    report = audit_contract(relaxed)
    assert report["passed"] is False
    assert "CR001 scientific contract drift" in report["errors"]


def test_source_classification_drift_is_rejected() -> None:
    contract = load_contract()
    relabeled = mutated(
        contract,
        "source_classification",
        "pending_unknown",
        value=["hidden OpenAI coefficients only"],
    )
    report = audit_contract(relabeled)
    assert report["passed"] is False
    assert any("transitive runtime identity must remain pending" in e for e in report["errors"])
