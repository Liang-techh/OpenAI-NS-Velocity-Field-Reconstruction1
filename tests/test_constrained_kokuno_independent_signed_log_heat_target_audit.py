from __future__ import annotations

import math

import pytest

from openai_ns_reconstruction.kokuno_independent_signed_log_heat_target_audit import (
    FORMAL_DIVERGENCE_GATE,
    FORMAL_MOMENTUM_GATE,
    LOCAL_GUARDS,
    run_audit,
)


@pytest.fixture(scope="module")
def report() -> dict:
    return run_audit()


def test_independent_signed_log_heat_target_preflight_passes(report: dict) -> None:
    assert report["structural_preflight_passed"] is True
    assert all(report["guards"].values())

    metrics = report["metrics"]
    assert metrics["max_public_log_abs_error"] <= LOCAL_GUARDS["max_public_log_abs_error"]
    assert (
        metrics["max_medium_to_fine_log_change"]
        <= LOCAL_GUARDS["max_medium_to_fine_log_change"]
    )
    assert (
        metrics["max_materialized_relative_error"]
        <= LOCAL_GUARDS["max_materialized_relative_error"]
    )
    assert (
        metrics["minimum_old_tail_bug_log_shift"]
        >= LOCAL_GUARDS["minimum_old_tail_bug_log_shift"]
    )


def test_default_case_preserves_underflow_hierarchy_and_nonzero_signs(report: dict) -> None:
    default = next(case for case in report["cases"] if case["name"] == "default_extreme")
    target = default["finest_eta_0p2"]

    assert target["sign"] == [1, -1, 1]
    assert target["underflow_channels"] == ["C_p", "S"]
    assert target["span_decades"] >= LOCAL_GUARDS["minimum_default_span_decades"]

    log10_abs = target["log10_abs"]
    assert log10_abs[0] == pytest.approx(-664.7, abs=0.2)
    assert log10_abs[1] == pytest.approx(-464.7, abs=0.2)
    assert log10_abs[2] == pytest.approx(-243.3, abs=0.2)


def test_small_z_oracle_is_far_inside_preregistered_regime(report: dict) -> None:
    metrics = report["metrics"]
    assert metrics["maximum_tested_log_z"] <= LOCAL_GUARDS["maximum_tested_log_z"]

    for case in report["cases"]:
        for level in case["levels"]:
            for row in level["rows"]:
                omitted = row["first_omitted_relative_log_scale"]
                assert math.isfinite(omitted)
                assert omitted < -80.0


def test_old_terminal_cp_factor_bug_is_detected(report: dict) -> None:
    for case in report["cases"]:
        assert (
            case["old_tail_bug_cp_log_shift"]
            >= LOCAL_GUARDS["minimum_old_tail_bug_log_shift"]
        )


def test_formal_pde_truth_boundary_stays_false(report: dict) -> None:
    assert report["formal_full_domain_pde_gate_assessed"] is False
    assert report["normalized_ns_residual_le_1e-3_claimed"] is False
    assert report["pde_validated"] is False
    assert report["formal_momentum_gate"] == FORMAL_MOMENTUM_GATE == 1.0e-3
    assert report["formal_divergence_gate"] == FORMAL_DIVERGENCE_GATE == 1.0e-5
    assert report["st006_context"]["directly_comparable_to_this_local_audit"] is False
