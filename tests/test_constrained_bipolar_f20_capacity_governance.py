import copy
import json
from pathlib import Path

import pytest

from openai_ns_reconstruction.constrained_bipolar_f20_capacity_governance import (
    audit_bipolar_f20_capacity_scope,
)
from openai_ns_reconstruction.constrained_eq45_bipolar_f20_capacity import (
    audit_bipolar_f20_capacity,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(path):
    with (ROOT / path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _inputs():
    scope = _load("configs/bipolar_f20_capacity_scope.json")
    constraints = _load(scope["inputs"]["constraints"])
    delivery = _load(scope["inputs"]["delivery_state_contract"])
    candidate = _load(scope["inputs"]["source_candidate"])
    report = audit_bipolar_f20_capacity()
    return scope, constraints, delivery, candidate, report


def test_f20_capacity_scope_preserves_selection_pde_and_energy_boundaries():
    scope, constraints, delivery, candidate, report = _inputs()
    result = audit_bipolar_f20_capacity_scope(
        scope, constraints, delivery, candidate, report
    )
    assert result["scope_pass"] is True
    assert result["f20_mode"] == [2, 0]
    assert result["f20_source_classification"] == "autonomous_design"
    assert result["f20_nonzero_candidate_materialized"] is False
    assert result["candidate_selection_resolved"] is False
    assert result["held_out_pde_residual_evaluated"] is False
    assert result["pde_validated"] is False
    assert result["visual_correspondence_verified"] is False

    mutations = []

    mutated = copy.deepcopy(scope)
    mutated["information_source"]["f20_basis_mode"] = "public_source_fact"
    mutations.append(mutated)

    mutated = copy.deepcopy(scope)
    mutated["forbidden_inferences"][
        "swirl_pure_velocity_response_implies_theta_only_momentum_change"
    ] = True
    mutations.append(mutated)

    mutated = copy.deepcopy(scope)
    mutated["forbidden_inferences"][
        "swirl_pure_velocity_response_implies_energy_unchanged"
    ] = True
    mutations.append(mutated)

    mutated = copy.deepcopy(scope)
    mutated["capacity_evidence_contract"]["coefficient_selected"] = True
    mutations.append(mutated)

    mutated = copy.deepcopy(scope)
    mutated["materialization_state"]["f20_nonzero_candidate_materialized"] = True
    mutated["materialization_state"]["f20_nonzero_candidate_sha256"] = "0" * 64
    mutated["materialization_state"]["f20_nonzero_velocity_export_ready"] = True
    mutations.append(mutated)

    mutated = copy.deepcopy(scope)
    mutated["requirements_before_nonzero_f20_promotion"][
        "reference_energy_revalidated"
    ] = False
    mutations.append(mutated)

    mutated = copy.deepcopy(scope)
    mutated["requirements_before_nonzero_f20_promotion"][
        "parent_energy_or_pde_results_not_inherited"
    ] = False
    mutations.append(mutated)

    for bad_scope in mutations:
        with pytest.raises(ValueError):
            audit_bipolar_f20_capacity_scope(
                bad_scope, constraints, delivery, candidate, report
            )


def test_f20_capacity_scope_rejects_registered_threshold_or_report_promotion():
    scope, constraints, delivery, candidate, report = _inputs()

    bad_constraints = copy.deepcopy(constraints)
    bad_constraints["validation"]["thresholds"]["pde_residual_max"] = 0.01
    with pytest.raises(ValueError):
        audit_bipolar_f20_capacity_scope(
            scope, bad_constraints, delivery, candidate, report
        )

    bad_report = copy.deepcopy(report)
    bad_report["truth_boundary"]["held_out_pde_residual_evaluated"] = True
    with pytest.raises(ValueError):
        audit_bipolar_f20_capacity_scope(
            scope, constraints, delivery, candidate, bad_report
        )

    bad_report = copy.deepcopy(report)
    bad_report["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError):
        audit_bipolar_f20_capacity_scope(
            scope, constraints, delivery, candidate, bad_report
        )

    bad_report = copy.deepcopy(report)
    bad_report["source_candidate_sha256"] = "f" * 64
    with pytest.raises(ValueError):
        audit_bipolar_f20_capacity_scope(
            scope, constraints, delivery, candidate, bad_report
        )
