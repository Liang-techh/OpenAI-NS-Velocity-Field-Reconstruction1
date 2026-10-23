from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_derived_challenger_lineage_governance import (
    audit_derived_challenger_lineage,
    load_constraints,
    load_contract,
    load_project_status,
    validate_parent_lineage_record,
)


def _inputs():
    return load_contract(), load_constraints(), load_project_status()


def test_current_derived_challenger_lineage_contract_passes():
    report = audit_derived_challenger_lineage()
    assert report["status"] == "governance_pass"
    assert report["observed_child_pr"] == 407
    assert report["observed_parent_pr"] == 398
    assert report["parent_identity_is_immutable"] is True
    assert report["branch_tip_alone_is_insufficient"] is True
    assert report["stacked_child_is_protocol_scoped_only"] is True
    assert report["parent_evidence_transfers_automatically"] is False
    assert report["canonical_velocity_unchanged"] is True
    assert report["canonical_thresholds_unchanged"] is True
    assert report["velocity_export_ready"] is True
    assert report["pde_validated"] is False
    assert report["visual_correspondence_verified"] is False


def test_branch_name_alone_is_not_a_parent_identity():
    with pytest.raises(ValueError, match="immutable parent git commit required"):
        validate_parent_lineage_record({
            "parent_branch": "agent7/cr003-st048-morphology-transfer-060",
            "parent_representation_identity": "ST048-S",
            "parent_candidate_sha256": "6" * 64,
        })

    with pytest.raises(ValueError, match="parent representation identity required"):
        validate_parent_lineage_record({
            "parent_branch": "some-branch",
            "parent_git_commit": "1" * 40,
            "parent_candidate_sha256": "6" * 64,
        })

    with pytest.raises(ValueError, match="parent candidate sha256 required"):
        validate_parent_lineage_record({
            "parent_branch": "some-branch",
            "parent_git_commit": "1" * 40,
            "parent_representation_identity": "ST048-S",
        })


def test_parent_or_child_status_cannot_transfer_through_stacked_lineage():
    contract, constraints, status = _inputs()
    for key in (
        "parent_validation_receipt_transfers_automatically",
        "parent_pressure_or_forcing_transfers_automatically",
        "parent_retained_baseline_status_transfers_automatically",
        "parent_canonical_delivery_status_transfers_automatically",
        "parent_pde_status_transfers_automatically",
        "parent_visual_status_transfers_automatically",
        "parent_paper_exact_status_transfers_automatically",
        "parent_openai_identity_status_transfers_automatically",
    ):
        bad = deepcopy(contract)
        bad["lineage_contract"][key] = True
        with pytest.raises(ValueError):
            audit_derived_challenger_lineage(bad, constraints, status)

    bad = deepcopy(contract)
    bad["current_child_snapshot"]["parent_full_momentum_may_be_inherited"] = True
    with pytest.raises(ValueError):
        audit_derived_challenger_lineage(bad, constraints, status)


def test_unintegrated_parent_dependency_must_be_resolved_before_promotion():
    contract, constraints, status = _inputs()
    bad = deepcopy(contract)
    bad["unintegrated_parent_dependency"]["promotion_requires_explicit_dependency_resolution"] = False
    with pytest.raises(ValueError, match="dependency resolution requirement"):
        audit_derived_challenger_lineage(bad, constraints, status)

    bad = deepcopy(contract)
    bad["unintegrated_parent_dependency"]["downstream_green_ci_promotes_parent"] = True
    with pytest.raises(ValueError, match="downstream_green_ci_promotes_parent"):
        audit_derived_challenger_lineage(bad, constraints, status)

    bad = deepcopy(contract)
    bad["unintegrated_parent_dependency"]["branch_tip_only_dependency_is_acceptable_for_promotion"] = True
    with pytest.raises(ValueError, match="branch_tip_only_dependency_is_acceptable_for_promotion"):
        audit_derived_challenger_lineage(bad, constraints, status)


def test_capacity_or_divergence_preflight_cannot_promote_child_or_pde_state():
    contract, constraints, status = _inputs()
    for key in (
        "open_or_stacked_child_replaces_retained_baseline",
        "open_or_stacked_child_replaces_canonical_velocity",
        "capacity_crossing_selects_production_parameter",
        "divergence_preservation_implies_pde_validated",
    ):
        bad = deepcopy(contract)
        bad["child_promotion_boundary"][key] = True
        with pytest.raises(ValueError):
            audit_derived_challenger_lineage(bad, constraints, status)

    bad = deepcopy(contract)
    bad["current_child_snapshot"]["production_beta_selected"] = True
    with pytest.raises(ValueError, match="production_beta_selected"):
        audit_derived_challenger_lineage(bad, constraints, status)

    bad = deepcopy(contract)
    bad["current_child_snapshot"]["pde_validated"] = True
    with pytest.raises(ValueError, match="pde_validated"):
        audit_derived_challenger_lineage(bad, constraints, status)


def test_source_classes_delivery_route_and_cr001_gates_fail_closed():
    contract, constraints, status = _inputs()
    bad = deepcopy(contract)
    bad["allowed_source_classes"].append("project_derived_measurement")
    with pytest.raises(ValueError, match="source-class vocabulary"):
        audit_derived_challenger_lineage(bad, constraints, status)

    bad = deepcopy(contract)
    bad["source_classification"]["openai_hidden_coordinate_map"] = "public_source_fact"
    with pytest.raises(ValueError, match="source class openai_hidden_coordinate_map"):
        audit_derived_challenger_lineage(bad, constraints, status)

    bad_status = deepcopy(status)
    bad_status["velocity_api"] = "experimental.st048_warp:velocity"
    with pytest.raises(ValueError, match="canonical velocity API"):
        audit_derived_challenger_lineage(contract, constraints, bad_status)

    bad_status = deepcopy(status)
    bad_status["active_scientific_route"] = "promote_st048_warp"
    with pytest.raises(ValueError, match="active scientific route"):
        audit_derived_challenger_lineage(contract, constraints, bad_status)

    bad_constraints = deepcopy(constraints)
    bad_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError, match="PDE max"):
        audit_derived_challenger_lineage(contract, bad_constraints, status)

    bad_constraints = deepcopy(constraints)
    bad_constraints["forcing"]["parameters"]["c"] = [0.0, 20.0]
    with pytest.raises(ValueError, match="force bounds"):
        audit_derived_challenger_lineage(contract, bad_constraints, status)
