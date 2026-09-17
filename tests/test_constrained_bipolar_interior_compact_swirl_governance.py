from copy import deepcopy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_bipolar_interior_compact_swirl_governance import (
    audit_bipolar_interior_compact_swirl_scope,
    audit_bipolar_interior_compact_swirl_scope_file,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(relative):
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _inputs():
    scope = _load("configs/bipolar_interior_compact_swirl_capacity_scope.json")
    constraints = _load("configs/constraints.json")
    delivery = _load("configs/delivery_state_contract.json")
    source = _load(scope["source_candidate"]["artifact"])
    return scope, constraints, delivery, source


def test_interior_compact_swirl_scope_passes_current_state():
    receipt = audit_bipolar_interior_compact_swirl_scope_file(repo_root=ROOT)
    assert receipt["scope_pass"] is True
    assert receipt["mode_origin"] == "autonomous_design"
    assert receipt["is_existing_eq45_profile_basis_coefficient"] is False
    assert receipt["diagnostic_trial_selected"] is False
    assert receipt["nonzero_materialization_requires_new_identity"] is True
    assert receipt["pde_validated"] is False


def test_rejects_laundering_physical_space_mode_into_eq45_source_basis():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["representation_semantics"]["is_eq45_public_source_profile_mode"] = True
    with pytest.raises(ValueError, match="public Eq45 mode"):
        audit_bipolar_interior_compact_swirl_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["representation_semantics"]["is_existing_eq45_profile_basis_coefficient"] = True
    with pytest.raises(ValueError, match="profile coefficient"):
        audit_bipolar_interior_compact_swirl_scope(mutated, constraints, delivery, source)


def test_rejects_promoting_symmetric_trial_to_selection_or_bound():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["diagnostic_trial_semantics"]["trial_coefficients_are_selected"] = True
    with pytest.raises(ValueError, match="trial promotion"):
        audit_bipolar_interior_compact_swirl_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["diagnostic_trial_semantics"]["trial_coefficients_define_materialization_bound"] = True
    with pytest.raises(ValueError, match="trial promotion"):
        audit_bipolar_interior_compact_swirl_scope(mutated, constraints, delivery, source)


def test_rejects_parent_identity_reuse_or_scientific_promotion():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["representation_semantics"]["parent_candidate_identity_may_be_reused_for_nonzero_child"] = True
    with pytest.raises(ValueError, match="identity reuse"):
        audit_bipolar_interior_compact_swirl_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["states"]["interior_compact_nonzero_candidate_materialized"] = True
    with pytest.raises(ValueError, match="premature state promotion"):
        audit_bipolar_interior_compact_swirl_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="premature state promotion"):
        audit_bipolar_interior_compact_swirl_scope(mutated, constraints, delivery, source)


def test_rejects_registered_gate_or_force_relaxation():
    scope, constraints, delivery, source = _inputs()
    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="registered threshold drift"):
        audit_bipolar_interior_compact_swirl_scope(scope, mutated_constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["registered_contract"]["residual_defined_free_force_allowed"] = True
    with pytest.raises(ValueError, match="free-force route enabled"):
        audit_bipolar_interior_compact_swirl_scope(mutated, constraints, delivery, source)
