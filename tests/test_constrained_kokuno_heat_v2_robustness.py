import json

import pytest

from openai_ns_reconstruction.kokuno_heat_discrepancy_v2_robustness import (
    PARAMETER_CASES,
    QUADRATURE_LEVELS,
    run_audit,
)


@pytest.fixture(scope="module")
def report():
    return run_audit()


def test_corrected_heat_discrepancy_is_robust_under_fixed_parameter_perturbations(report):
    assert [case["name"] for case in PARAMETER_CASES] == [
        "lower_scale",
        "baseline",
        "upper_scale",
    ]
    assert [level["label"] for level in QUADRATURE_LEVELS] == [
        "coarse",
        "medium",
        "fine",
    ]
    assert report["structural_preflight_passed"] is True
    assert all(report["checks"].values())


def test_independent_audit_detects_the_removed_factor_mutation(report):
    summary = report["summary"]
    assert summary["missing_copy_mutation_min_cp_relative_error"] >= 1.0e-4
    assert summary["tail_fix_vs_independent_half_first_order_max_relative"] <= 5.0e-5
    assert summary["unchanged_delta_s_i_max_abs"] == 0.0
    assert summary["endpoint_eta_pm1_max_abs"] == 0.0


def test_resolution_and_truth_boundaries_remain_fail_closed(report):
    summary = report["summary"]
    assert summary["medium_to_fine_max_relative"] <= 5.0e-5
    assert report["st006_directly_comparable"] is False
    assert report["fixed_project_references"] == {
        "normalized_full_momentum": 1.0e-3,
        "divergence_max": 1.0e-5,
        "changed": False,
    }
    assert report["formal_full_domain_pde_gate_assessed"] is False
    assert report["normalized_ns_residual_le_1e-3_claimed"] is False
    assert report["pde_validated"] is False
    json.dumps(report, allow_nan=False)
