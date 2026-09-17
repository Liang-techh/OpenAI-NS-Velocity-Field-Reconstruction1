from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_bipolar_error_map_routing_governance import (
    audit_bipolar_error_map_routing_scope,
    audit_bipolar_error_map_routing_scope_file,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(path: str) -> dict:
    with (ROOT / path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _inputs() -> tuple[dict, dict, dict, dict]:
    scope = _load("configs/bipolar_error_map_routing_scope.json")
    constraints = _load("configs/constraints.json")
    delivery = _load("configs/delivery_state_contract.json")
    candidate = _load(scope["source_candidate"]["artifact"])
    return scope, constraints, delivery, candidate


def test_current_error_map_routing_scope_passes() -> None:
    report = audit_bipolar_error_map_routing_scope_file(repo_root=ROOT)
    assert report["scope_pass"] is True
    assert report["error_map_role"] == "routing_diagnostic_only"
    assert report["f20_selected"] is False
    assert report["interior_localized_mode_materialized"] is False
    assert report["model_selection_samples_reusable_as_independent_validation"] is False
    assert report["fresh_or_pristine_validation_required_after_freeze"] is True
    assert report["pde_validated"] is False


def test_residual_map_cannot_fit_spatial_envelope() -> None:
    scope, constraints, delivery, candidate = _inputs()
    mutated = copy.deepcopy(scope)
    mutated["model_selection_boundary"]["residual_map_may_fit_spatial_envelope"] = True
    with pytest.raises(ValueError, match="may not fit correction envelope"):
        audit_bipolar_error_map_routing_scope(mutated, constraints, delivery, candidate)


def test_model_selection_samples_cannot_be_reused_as_independent_validation() -> None:
    scope, constraints, delivery, candidate = _inputs()
    mutated = copy.deepcopy(scope)
    mutated["model_selection_boundary"]["consumed_samples_may_count_as_independent_validation"] = True
    with pytest.raises(ValueError, match="laundered into validation"):
        audit_bipolar_error_map_routing_scope(mutated, constraints, delivery, candidate)


def test_frozen_candidate_requires_fresh_or_pristine_validation() -> None:
    scope, constraints, delivery, candidate = _inputs()
    mutated = copy.deepcopy(scope)
    mutated["model_selection_boundary"]["required_after_candidate_freeze"] = "reuse_overlap_map"
    with pytest.raises(ValueError, match="fresh validation requirement drift"):
        audit_bipolar_error_map_routing_scope(mutated, constraints, delivery, candidate)


def test_routing_evidence_cannot_select_f20() -> None:
    scope, constraints, delivery, candidate = _inputs()
    mutated = copy.deepcopy(scope)
    mutated["states"]["f20_selected"] = True
    with pytest.raises(ValueError, match="premature state promotion: f20_selected"):
        audit_bipolar_error_map_routing_scope(mutated, constraints, delivery, candidate)


def test_routing_evidence_ordering_is_pinned() -> None:
    scope, constraints, delivery, candidate = _inputs()
    mutated = copy.deepcopy(scope)
    mutated["upstream_routing_evidence"]["f20_vs_f10"]["response_error_cosine_f20"] = 0.8
    with pytest.raises(ValueError, match="routing evidence drift: F20 error cosine"):
        audit_bipolar_error_map_routing_scope(mutated, constraints, delivery, candidate)


def test_registered_pde_threshold_cannot_be_relaxed() -> None:
    scope, constraints, delivery, candidate = _inputs()
    mutated = copy.deepcopy(scope)
    mutated["registered_contract"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="registered threshold drift: pde_residual_max"):
        audit_bipolar_error_map_routing_scope(mutated, constraints, delivery, candidate)


def test_source_classification_vocabulary_cannot_expand() -> None:
    scope, constraints, delivery, candidate = _inputs()
    mutated = copy.deepcopy(scope)
    mutated["classification_vocabulary"].append("parent_provenance")
    with pytest.raises(ValueError, match="scope source vocabulary drift"):
        audit_bipolar_error_map_routing_scope(mutated, constraints, delivery, candidate)
