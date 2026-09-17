from copy import deepcopy
import json
import math
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_bipolar_interior_swirl_redistribution_governance import (
    audit_bipolar_interior_swirl_redistribution_scope,
    audit_bipolar_interior_swirl_redistribution_scope_file,
)
from openai_ns_reconstruction.constrained_eq45_bipolar_interior_swirl_redistribution_capacity import (
    ORTHOGONALIZER,
    REDISTRIBUTION_MODE,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(relative):
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _inputs():
    scope = _load("configs/bipolar_interior_swirl_redistribution_capacity_scope.json")
    constraints = _load("configs/constraints.json")
    delivery = _load("configs/delivery_state_contract.json")
    source = _load(scope["source_candidate"]["artifact"])
    return scope, constraints, delivery, source


def test_redistribution_scope_passes_current_stacked_state():
    receipt = audit_bipolar_interior_swirl_redistribution_scope_file(repo_root=ROOT)
    assert receipt["scope_pass"] is True
    assert receipt["mode"] == REDISTRIBUTION_MODE == "INTERIOR_C4_SWIRL_REDISTRIBUTION"
    assert receipt["mode_origin"] == "autonomous_design"
    assert receipt["raw_basis_pair_orthogonality_only"] is True
    assert receipt["full_candidate_energy_neutrality_claimed"] is False
    assert receipt["technical_guard_is_materialization_bound"] is False
    assert receipt["nonzero_materialization_requires_new_identity"] is True
    assert receipt["pde_validated"] is False


def test_scope_binds_agent7_analytic_orthogonalizer_without_promoting_energy_neutrality():
    scope, constraints, delivery, source = _inputs()
    assert math.isclose(ORTHOGONALIZER, 13.0 / 2.0, rel_tol=0.0, abs_tol=0.0)
    assert math.isclose(
        scope["upstream_capacity_evidence"]["orthogonalizer"],
        ORTHOGONALIZER,
        rel_tol=0.0,
        abs_tol=0.0,
    )
    assert math.isclose(
        scope["upstream_capacity_evidence"]["analytic_sign_change_r_over_Rp"],
        math.sqrt(2.0 / 13.0),
        rel_tol=1e-12,
    )

    mutated = deepcopy(scope)
    mutated["orthogonality_semantics"]["orthogonality_implies_full_candidate_energy_neutrality"] = True
    with pytest.raises(ValueError, match="orthogonality/collar inference promoted"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["orthogonality_semantics"]["orthogonality_implies_finite_coefficient_energy_neutrality"] = True
    with pytest.raises(ValueError, match="orthogonality/collar inference promoted"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)


def test_rejects_theta_only_nonlinear_momentum_inference_or_hidden_quadratic_energy_term():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["orthogonality_semantics"]["pure_swirl_velocity_response_implies_theta_only_nonlinear_momentum_impact"] = True
    with pytest.raises(ValueError, match="orthogonality/collar inference promoted"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["orthogonality_semantics"]["finite_nonzero_coefficient_has_positive_quadratic_energy_term"] = False
    with pytest.raises(ValueError, match="quadratic energy term hidden"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)


def test_rejects_borrowed_capacity_guard_or_trial_as_materialization_bound():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["coefficient_semantics"]["technical_guard_is_preregistered_materialization_bound"] = True
    with pytest.raises(ValueError, match="technical guard laundered"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["coefficient_semantics"]["trial_coefficients_define_materialization_bound"] = True
    with pytest.raises(ValueError, match="diagnostic coefficient promotion"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["coefficient_semantics"]["trial_coefficients_are_selected"] = True
    with pytest.raises(ValueError, match="diagnostic coefficient promotion"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)


def test_rejects_parent_identity_reuse_or_scientific_promotion():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["representation_semantics"]["parent_candidate_identity_may_be_reused_for_nonzero_child"] = True
    with pytest.raises(ValueError, match="identity reuse"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["states"]["interior_redistribution_nonzero_candidate_materialized"] = True
    with pytest.raises(ValueError, match="premature state promotion"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="premature state promotion"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)


def test_rejects_registered_gate_or_force_relaxation():
    scope, constraints, delivery, source = _inputs()
    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_L2"] = 0.01
    with pytest.raises(ValueError, match="registered threshold drift"):
        audit_bipolar_interior_swirl_redistribution_scope(scope, mutated_constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["registered_contract"]["residual_defined_free_force_allowed"] = True
    with pytest.raises(ValueError, match="free-force route enabled"):
        audit_bipolar_interior_swirl_redistribution_scope(mutated, constraints, delivery, source)
