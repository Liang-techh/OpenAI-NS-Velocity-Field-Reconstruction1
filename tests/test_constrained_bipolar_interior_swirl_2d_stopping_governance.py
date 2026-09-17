from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_bipolar_interior_swirl_2d_stopping_governance import (
    audit_bipolar_interior_swirl_2d_stopping_scope,
    audit_repository,
)
from openai_ns_reconstruction.constrained_eq45_bipolar_interior_swirl_2d_jacobian import (
    TRUTH_BOUNDARY as UPSTREAM_TRUTH_BOUNDARY,
)

ROOT = Path(__file__).resolve().parents[1]


def _payloads():
    scope = json.loads((ROOT / "configs/bipolar_interior_swirl_2d_stopping_scope.json").read_text())
    constraints = json.loads((ROOT / "configs/constraints.json").read_text())
    delivery = json.loads((ROOT / "configs/delivery_state_contract.json").read_text())
    return scope, constraints, delivery


def test_current_joint_stopping_scope_passes():
    result = audit_repository()
    assert result["task_id"] == "CR002-BIPOLAR-INTERIOR-SWIRL-2D-STOPPING-SCOPE-035"
    assert result["registered_validation_points"] == 4096
    assert result["pde_validated"] is False
    assert UPSTREAM_TRUTH_BOUNDARY["new_basis_mode_added"] is False
    assert UPSTREAM_TRUTH_BOUNDARY["radial_redistribution_coefficient_selected"] is False
    assert UPSTREAM_TRUTH_BOUNDARY["axial_redistribution_coefficient_selected"] is False
    assert UPSTREAM_TRUTH_BOUNDARY["perturbed_candidate_pde_residual_evaluated"] is False
    assert UPSTREAM_TRUTH_BOUNDARY["pde_validated"] is False


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("joint_capacity_semantics", "public_velocity_rank_is_residual_jacobian_rank", True),
        ("joint_capacity_semantics", "morphology_rank_implies_visual_correspondence", True),
        ("joint_capacity_semantics", "raw_response_orthogonality_implies_finite_child_energy_neutrality", True),
        ("stopping_rule_semantics", "global_no_mixed_mode_theorem", True),
        ("stopping_rule_semantics", "candidate_selection_resolved", True),
        ("materialization_semantics", "diagnostic_plus_minus_0p10_selects_sign_or_value", True),
        ("materialization_semantics", "inherited_abs_4_implementation_guard_is_preregistered_bound_for_autonomous_modes", True),
        ("truth_boundary", "pde_validated", True),
        ("truth_boundary", "visual_correspondence_verified", True),
    ],
)
def test_claim_laundering_mutations_fail_closed(section, key, value):
    scope, constraints, delivery = _payloads()
    mutated = copy.deepcopy(scope)
    mutated[section][key] = value
    with pytest.raises(ValueError):
        audit_bipolar_interior_swirl_2d_stopping_scope(mutated, constraints, delivery)


def test_threshold_relaxation_fails_closed():
    scope, constraints, delivery = _payloads()
    mutated = copy.deepcopy(scope)
    mutated["registered_contract"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="threshold drift"):
        audit_bipolar_interior_swirl_2d_stopping_scope(mutated, constraints, delivery)


def test_source_vocabulary_collapse_fails_closed():
    scope, constraints, delivery = _payloads()
    mutated = copy.deepcopy(scope)
    mutated["classification_vocabulary"].remove("pending_unknown")
    with pytest.raises(ValueError, match="source vocabulary"):
        audit_bipolar_interior_swirl_2d_stopping_scope(mutated, constraints, delivery)
