from __future__ import annotations

from copy import deepcopy
import json

import pytest

from openai_ns_reconstruction.constrained_st052m_axis_representation_audit import (
    ROOT,
    _audit_cr001,
    audit_repository,
    load_contract,
    mutation_rejects_irregular_parent,
    validate_contract,
)


def test_live_st052_axis_representation_audit_passes() -> None:
    receipt = audit_repository()
    assert receipt["exact_axis_division_safe"] is True
    assert receipt["conditional_transverse_O_r_preserved"] is True
    assert receipt["hard_radius_guard_present"] is True
    assert receipt["hard_radius_guard"] == 1.0e-14
    assert receipt["analytic_C1_or_smoother_axis_regularity_proved"] is False
    assert receipt["whole_parent_axis_regularity_independently_proved_here"] is False
    assert receipt["whole_candidate_axis_regularity_theorem_proved"] is False
    assert receipt["pde_validated"] is False
    assert receipt["visual_correspondence_verified"] is False
    assert receipt["cr001_unchanged"] is True

    synthetic = receipt["synthetic_regression"]
    assert synthetic["parent_precondition"]["passed"] is True
    assert synthetic["max_axis_transverse"]["static"] <= 1.0e-14
    assert synthetic["max_axis_transverse"]["temporal"] <= 1.0e-14
    assert synthetic["max_transverse_over_r"]["static"] < 20.0
    assert synthetic["max_transverse_over_r"]["temporal"] < 20.0
    assert synthetic["static_guard_effect_in_tangential_over_r"] > 0.001
    assert synthetic["temporal_guard_effect_in_tangential_over_r"] > 0.001


def test_parent_axis_prerequisite_is_real_not_assumed_away() -> None:
    assert mutation_rejects_irregular_parent() is True


@pytest.mark.parametrize(
    "key",
    [
        "analytic_C1_or_smoother_axis_regularity_proved",
        "whole_parent_axis_regularity_independently_proved_here",
        "whole_candidate_axis_regularity_theorem_proved",
    ],
)
def test_analytic_axis_claims_fail_closed(key: str) -> None:
    contract = deepcopy(load_contract())
    contract["axis_representation"][key] = True
    with pytest.raises(ValueError, match="premature analytic axis claim"):
        validate_contract(contract)


def test_radius_guard_cannot_be_relabelled_as_source_fact() -> None:
    contract = deepcopy(load_contract())
    contract["axis_representation"]["hard_radius_guard_is_autonomous_numerical_design"] = False
    with pytest.raises(ValueError, match="radius guard must remain autonomous design"):
        validate_contract(contract)


def test_pde_claim_cannot_be_promoted_by_axis_safety() -> None:
    contract = deepcopy(load_contract())
    contract["claim_states"]["pde_validated"] = True
    with pytest.raises(ValueError, match="premature claim promotion"):
        validate_contract(contract)


def test_source_classes_remain_separate() -> None:
    contract = deepcopy(load_contract())
    contract["source_classification"][3]["classification"] = "public_source_fact"
    with pytest.raises(ValueError, match="source classification classes drifted"):
        validate_contract(contract)


def test_cr001_threshold_mutation_is_rejected() -> None:
    contract = load_contract()
    constraints = json.loads((ROOT / "configs" / "constraints.json").read_text(encoding="utf-8"))
    mutated = deepcopy(constraints)
    mutated["validation"]["thresholds"]["pde_residual_L2"] = 0.01
    with pytest.raises(ValueError, match="CR001 threshold drifted: pde_residual_L2"):
        _audit_cr001(contract, mutated)


def test_cr001_free_force_mutation_is_rejected() -> None:
    contract = load_contract()
    constraints = json.loads((ROOT / "configs" / "constraints.json").read_text(encoding="utf-8"))
    mutated = deepcopy(constraints)
    mutated["forcing"]["restriction"] = "Allow pointwise residual-defined free force."
    with pytest.raises(ValueError, match="free residual-defined forcing prohibition missing"):
        _audit_cr001(contract, mutated)
