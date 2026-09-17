from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_bipolar_parent_provenance_governance import (
    audit_bipolar_parent_provenance_scope,
    audit_bipolar_parent_provenance_scope_file,
)


ROOT = Path(__file__).resolve().parents[1]
SCOPE_PATH = ROOT / "configs" / "bipolar_parent_provenance_scope.json"


def _load(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _inputs() -> tuple[dict, dict, dict, dict, dict]:
    scope = _load(SCOPE_PATH)
    candidate = _load(ROOT / scope["inputs"]["candidate"])
    delivery = _load(ROOT / scope["inputs"]["delivery_state_contract"])
    constraints = _load(ROOT / scope["inputs"]["constraints"])
    status = _load(ROOT / scope["inputs"]["project_status"])
    return scope, candidate, delivery, constraints, status


def test_current_bipolar_parent_provenance_scope_passes() -> None:
    result = audit_bipolar_parent_provenance_scope_file(SCOPE_PATH, repo_root=ROOT)
    assert result["scope_pass"] is True
    assert result["legacy_parent_relation"] == "parent_provenance"
    assert result["legacy_parent_relation_is_source_classification"] is False
    assert result["canonical_source_classes"] == [
        "autonomous_design",
        "pending_unknown",
        "public_source_fact",
        "user_requirement",
    ]
    assert result["velocity_export_ready"] is True
    assert result["pde_validated"] is False
    assert result["paper_exact"] is False


def test_parent_relation_cannot_be_relabelled_as_public_source() -> None:
    scope, candidate, delivery, constraints, status = _inputs()
    mutant = copy.deepcopy(candidate)
    mutant["classification"]["eq45_parent"] = "public_source_fact"
    with pytest.raises(ValueError, match="parent relation changed"):
        audit_bipolar_parent_provenance_scope(scope, mutant, delivery, constraints, status)


def test_parent_provenance_cannot_become_a_fifth_source_class() -> None:
    scope, candidate, delivery, constraints, status = _inputs()
    mutant = copy.deepcopy(delivery)
    mutant["classification_vocabulary"]["parent_provenance"] = "incorrect fifth class"
    with pytest.raises(ValueError, match="canonical source vocabulary drift"):
        audit_bipolar_parent_provenance_scope(scope, candidate, mutant, constraints, status)


def test_real_information_class_must_stay_canonical() -> None:
    scope, candidate, delivery, constraints, status = _inputs()
    mutant = copy.deepcopy(candidate)
    mutant["classification"]["velocity_interface"] = "user_required"
    with pytest.raises(ValueError, match="artifact source classification drift"):
        audit_bipolar_parent_provenance_scope(scope, mutant, delivery, constraints, status)


def test_bipolar_profile_choice_cannot_be_laundered_as_public_fact() -> None:
    scope, candidate, delivery, constraints, status = _inputs()
    mutant = copy.deepcopy(scope)
    mutant["information_classification"]["bipolar_phi01_profile_choice"] = "public_source_fact"
    with pytest.raises(ValueError, match="information classification drift"):
        audit_bipolar_parent_provenance_scope(mutant, candidate, delivery, constraints, status)


def test_pde_claim_and_registered_threshold_drift_fail_closed() -> None:
    scope, candidate, delivery, constraints, status = _inputs()

    candidate_mutant = copy.deepcopy(candidate)
    candidate_mutant["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="candidate truth-state promotion: pde_validated"):
        audit_bipolar_parent_provenance_scope(
            scope, candidate_mutant, delivery, constraints, status
        )

    constraint_mutant = copy.deepcopy(constraints)
    constraint_mutant["validation"]["thresholds"]["pde_residual_L2"] = 0.01
    with pytest.raises(ValueError, match="registered threshold drift: pde_residual_L2"):
        audit_bipolar_parent_provenance_scope(
            scope, candidate, delivery, constraint_mutant, status
        )
