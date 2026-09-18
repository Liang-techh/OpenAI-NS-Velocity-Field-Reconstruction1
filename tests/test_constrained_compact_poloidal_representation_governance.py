from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_compact_poloidal_representation_governance import (
    audit_compact_poloidal_scope,
    load_constraints,
    load_project_status,
    load_scope,
)


def test_current_compact_poloidal_scope_passes():
    report = audit_compact_poloidal_scope()
    assert report["status"] == "governance_pass"
    assert report["mode_name"] == "COMPACT_C4_ODD_Z_POLOIDAL"
    assert report["canonical_thresholds_unchanged"] is True
    assert report["nonzero_child_materialized"] is False
    assert report["coefficient_selected"] is False
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False


def test_rejects_source_laundering_and_bound_promotion():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["source_classification"]["compact_poloidal_shape"] = "public_source_fact"
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["capacity_evidence_only"]["diagnostic_coefficients_define_materialization_bound"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["capacity_evidence_only"]["inherited_profile_coefficient_limit_is_materialization_bound"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)


def test_rejects_capacity_to_scientific_claim_promotion():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["capacity_evidence_only"]["analytic_or_fd_divergence_response_is_registered_divergence_acceptance"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["capacity_evidence_only"]["morphology_leverage_is_visual_correspondence"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["truth_boundary"]["visual_correspondence_verified"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)


def test_rejects_missing_new_identity_or_fresh_revalidation():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["future_nonzero_materialization"]["new_candidate_sha_required"] = False
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["future_nonzero_materialization"]["new_representation_family_identity_required"] = False
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["future_nonzero_materialization"]["fresh_full_per_component_momentum_divergence_required"] = False
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["future_nonzero_materialization"]["model_selection_data_may_not_be_reused_as_acceptance_data"] = False
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)


def test_rejects_cr001_threshold_or_force_drift():
    constraints = load_constraints()
    mutated = deepcopy(constraints)
    mutated["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(constraints=mutated)

    mutated = deepcopy(constraints)
    mutated["forcing"]["parameters"]["c"] = [0.0, 20.0]
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(constraints=mutated)

    mutated = deepcopy(constraints)
    mutated["validation"]["held_out_points"] = 2048
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(constraints=mutated)


def test_rejects_tip_geometry_or_project_state_laundering():
    scope = load_scope()
    mutated = deepcopy(scope)
    mutated["morphology_scope"]["diagnostic_q90_q99_axial_reach_moved"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    mutated = deepcopy(scope)
    mutated["morphology_scope"]["tip_geometry_solution_claimed"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(scope=mutated)

    status = load_project_status()
    mutated_status = deepcopy(status)
    mutated_status["states"]["paper_exact"] = True
    with pytest.raises(ValueError):
        audit_compact_poloidal_scope(project_status=mutated_status)
