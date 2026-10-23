from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_dependency_runtime_identity_governance import (
    audit_contract,
    load_contract,
)


def test_live_dependency_runtime_scope_audit_passes_without_promotion():
    report = audit_contract()
    assert report["passed"] is True
    assert report["source_runtime_code_identity_closed"] is True
    assert report["historical_module_cache_isolated"] is True
    assert report["whole_child_save_load_ready_with_authenticated_exact_source_runtime"] is True
    assert report["external_dependency_runtime_identity_bound"] is False
    assert report["complete_executable_environment_identity_closed"] is False
    assert report["standalone_reproducible_velocity_identity_ready"] is False
    assert report["velocity_export_ready"] is False
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False


@pytest.mark.parametrize(
    "key",
    [
        "python_micro_version_bound_in_whole_candidate_identity",
        "python_implementation_build_bound_in_whole_candidate_identity",
        "numpy_resolved_version_bound_in_whole_candidate_identity",
        "scipy_resolved_version_bound_in_whole_candidate_identity",
        "sympy_resolved_version_bound_in_whole_candidate_identity",
        "binary_wheel_or_blas_build_identity_bound",
        "external_dependency_runtime_identity_bound",
        "complete_executable_environment_identity_closed",
        "exact_ci_receipt_portable_across_dependency_versions",
    ],
)
def test_dependency_environment_promotion_laundering_is_rejected(key: str):
    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["execution_environment_scope"][key] = True
    with pytest.raises(ValueError):
        audit_contract(mutated)


@pytest.mark.parametrize(
    "key",
    [
        "standalone_package_parent_runtime_ready",
        "standalone_reproducible_velocity_identity_ready",
        "velocity_export_ready",
        "visualization_ready",
        "visual_correspondence_verified",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ],
)
def test_scientific_or_delivery_state_laundering_is_rejected(key: str):
    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["states"][key] = True
    with pytest.raises(ValueError):
        audit_contract(mutated)


def test_source_classification_drift_is_rejected():
    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["source_classification"]["paper_fact"] = mutated["source_classification"].pop(
        "autonomous_design"
    )
    with pytest.raises(ValueError, match="classification"):
        audit_contract(mutated)


def test_cr001_threshold_relaxation_is_rejected():
    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["cr001_snapshot"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="CR001"):
        audit_contract(mutated)


def test_cr001_forcing_expansion_is_rejected():
    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["cr001_snapshot"]["forcing_c_bounds"] = [0.0, 100.0]
    with pytest.raises(ValueError, match="CR001"):
        audit_contract(mutated)


def test_positive_integrated_runtime_facts_cannot_be_erased():
    contract = load_contract()
    mutated = deepcopy(contract)
    mutated["integrated_positive_facts"]["historical_source_module_graph_cache_isolated"] = False
    with pytest.raises(ValueError):
        audit_contract(mutated)
