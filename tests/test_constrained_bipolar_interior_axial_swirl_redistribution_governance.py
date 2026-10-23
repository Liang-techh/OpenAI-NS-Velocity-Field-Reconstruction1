from copy import deepcopy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_bipolar_interior_axial_swirl_redistribution_governance import (
    audit_bipolar_interior_axial_swirl_redistribution_scope,
    audit_bipolar_interior_axial_swirl_redistribution_scope_file,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(relative):
    return json.loads((ROOT / relative).read_text(encoding="utf-8"))


def _inputs():
    scope = _load("configs/bipolar_interior_axial_swirl_redistribution_scope.json")
    constraints = _load("configs/constraints.json")
    delivery = _load("configs/delivery_state_contract.json")
    source = _load(scope["source_candidate"]["artifact"])
    return scope, constraints, delivery, source


def test_axial_swirl_redistribution_scope_passes_current_state():
    receipt = audit_bipolar_interior_axial_swirl_redistribution_scope_file(repo_root=ROOT)
    assert receipt["scope_pass"] is True
    assert receipt["mode"] == "INTERIOR_C4_SWIRL_AXIAL_REDISTRIBUTION"
    assert receipt["mode_origin"] == "autonomous_design"
    assert receipt["raw_axial_orthogonality_only"] is True
    assert receipt["global_swirl_positivity_preregistered"] is False
    assert receipt["diagnostic_coefficient_selected"] is False
    assert receipt["nonzero_materialization_requires_new_identity"] is True
    assert receipt["pde_validated"] is False


def test_rejects_orthogonality_laundering_into_full_candidate_or_energy_neutrality():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["orthogonality_semantics"]["orthogonality_is_to_complete_candidate_velocity"] = True
    with pytest.raises(ValueError, match="raw orthogonality laundered"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["orthogonality_semantics"]["orthogonality_implies_finite_child_energy_unchanged"] = True
    with pytest.raises(ValueError, match="energy neutrality"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)


def test_rejects_invented_global_swirl_positivity_or_required_counter_rotation():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["local_rotation_sign_semantics"]["global_u_theta_positivity_is_preregistered"] = True
    with pytest.raises(ValueError, match="invented global swirl-positivity gate"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["local_rotation_sign_semantics"]["capacity_sign_change_is_not_a_required_global_counter_rotation"] = False
    with pytest.raises(ValueError, match="promoted to requirement"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)


def test_rejects_inherited_guard_or_diagnostic_probe_as_materialization_bound():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["coefficient_semantics"]["inherited_guard_is_preregistered_bound_for_this_autonomous_mode"] = True
    with pytest.raises(ValueError, match="inherited guard laundered"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["coefficient_semantics"]["diagnostic_coefficients_define_materialization_bound"] = True
    with pytest.raises(ValueError, match="diagnostic coefficients promoted"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)


def test_rejects_parent_identity_reuse_or_scientific_promotion():
    scope, constraints, delivery, source = _inputs()
    mutated = deepcopy(scope)
    mutated["representation_semantics"]["parent_candidate_identity_may_be_reused_for_nonzero_child"] = True
    with pytest.raises(ValueError, match="identity reuse"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="premature state promotion"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)


def test_rejects_registered_gate_or_force_relaxation():
    scope, constraints, delivery, source = _inputs()
    mutated_constraints = deepcopy(constraints)
    mutated_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="registered threshold drift"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(scope, mutated_constraints, delivery, source)

    mutated = deepcopy(scope)
    mutated["registered_contract"]["residual_defined_free_force_allowed"] = True
    with pytest.raises(ValueError, match="free-force route enabled"):
        audit_bipolar_interior_axial_swirl_redistribution_scope(mutated, constraints, delivery, source)
