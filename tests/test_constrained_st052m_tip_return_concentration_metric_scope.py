from __future__ import annotations

from copy import deepcopy

import pytest

from openai_ns_reconstruction.constrained_st052m_tip_return_concentration_metric_scope import (
    GovernanceError,
    audit,
    load_contract,
    repository_root,
)


def test_live_metric_scope_audits() -> None:
    result = audit()
    assert result["audited_pr"] == 750
    assert result["audited_head"] == "376223a3bb2dae0d53756ea48da26937a95b6e58"
    assert result["diagnostic_scope"] == "one_dimensional_axial_response_shape"
    assert result["absolute_or_3d_concentration_verified"] is False
    assert result["total_child_return_flow_concentration_verified"] is False
    assert result["pde_validated"] is False


@pytest.mark.parametrize(
    ("section", "key", "value"),
    [
        ("metric_semantics", "recalibrated_alpha_94_used", True),
        ("metric_semantics", "cylindrical_volume_jacobian_used", True),
        ("metric_semantics", "cartesian_or_cylindrical_3d_volume_integral_used", True),
        ("metric_semantics", "nonlinear_child_velocity_evaluated", True),
        ("metric_semantics", "absolute_correction_magnitude_compared", True),
        ("metric_semantics", "three_dimensional_correction_energy_compared", True),
        ("claim_states", "absolute_correction_magnitude_verified_by_750", True),
        ("claim_states", "three_dimensional_correction_energy_concentration_verified", True),
        ("claim_states", "nonlinear_child_total_field_return_flow_concentration_verified", True),
        ("claim_states", "nonlinear_child_740_passed_at_contract_freeze", True),
        ("claim_states", "visual_correspondence_verified", True),
        ("claim_states", "pde_validated", True),
        ("claim_states", "paper_exact", True),
        ("claim_states", "openai_field_identified", True),
    ],
)
def test_metric_or_claim_promotion_fails_closed(section: str, key: str, value: bool) -> None:
    contract = deepcopy(load_contract())
    contract[section][key] = value
    with pytest.raises(GovernanceError):
        audit(contract=contract)


@pytest.mark.parametrize(
    "key",
    [
        "axial_shape_concentration_implies_larger_absolute_physical_return_velocity",
        "axial_shape_concentration_implies_three_dimensional_correction_energy_concentration",
        "axial_shape_concentration_implies_total_child_return_flow_concentration",
        "axial_shape_concentration_implies_nonlinear_child_passes_740",
        "axial_shape_concentration_implies_visual_correspondence",
        "axial_shape_concentration_implies_pde_validation",
    ],
)
def test_forbidden_inference_cannot_be_promoted(key: str) -> None:
    contract = deepcopy(load_contract())
    contract["interpretation_rule"][key] = True
    with pytest.raises(GovernanceError):
        audit(contract=contract)


def test_source_mutation_coefficient_injection_fails_closed() -> None:
    root = repository_root()
    source = (root / "experiments/root_st052/agent7_st052m_tip_return_concentration_audit.py").read_text(
        encoding="utf-8"
    )
    needle = "    z = np.linspace(relocation.AXIAL_Z_INNER, relocation.AXIAL_Z_OUTER, int(n))"
    mutated = source.replace(needle, "    alpha = 1.0\n" + needle, 1)
    assert mutated != source
    with pytest.raises(GovernanceError):
        audit(audit_source=mutated)


def test_source_mutation_pde_promotion_fails_closed() -> None:
    root = repository_root()
    source = (root / "experiments/root_st052/agent7_st052m_tip_return_concentration_audit.py").read_text(
        encoding="utf-8"
    )
    mutated = source.replace('"pde_validated": False', '"pde_validated": True', 1)
    assert mutated != source
    with pytest.raises(GovernanceError):
        audit(audit_source=mutated)


def test_cr001_free_force_collapse_and_threshold_drift_fail_closed() -> None:
    contract = deepcopy(load_contract())
    contract["cr001_lock"]["residual_defined_or_pointwise_free_force_allowed"] = True
    with pytest.raises(GovernanceError):
        audit(contract=contract)

    contract = deepcopy(load_contract())
    contract["cr001_lock"]["amplitude_collapse_success_allowed"] = True
    with pytest.raises(GovernanceError):
        audit(contract=contract)

    contract = deepcopy(load_contract())
    contract["cr001_lock"]["pde_residual_max"] = 0.01
    with pytest.raises(GovernanceError):
        audit(contract=contract)


def test_later_physical_claim_requirements_cannot_be_weakened() -> None:
    contract = deepcopy(load_contract())
    contract["interpretation_rule"]["later_absolute_or_3d_claim_requires_alpha_bound_and_physical_measure"] = False
    with pytest.raises(GovernanceError):
        audit(contract=contract)

    contract = deepcopy(load_contract())
    contract["interpretation_rule"]["later_total_field_claim_requires_parent_plus_correction_evaluation"] = False
    with pytest.raises(GovernanceError):
        audit(contract=contract)
