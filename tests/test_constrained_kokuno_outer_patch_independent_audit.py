import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_outer_patch_independent_audit import (
    LOCAL_NORMALIZED_DIVERGENCE_TOLERANCE,
    MUTATION_MIN_NORMALIZED_DIVERGENCE,
    Q_RELATION_RELATIVE_TOLERANCE,
    SCALE_LAW_RELATIVE_TOLERANCE,
    SCALE_LOG_ABSOLUTE_TOLERANCE,
    VELOCITY_RELATIVE_TOLERANCE,
    _q_bisection,
    _reference_scales,
    audit_schedule,
    run_audit,
)
from openai_ns_reconstruction.kokuno_outer_reserved_patch_schedule import (
    KokunoOuterReservedPatchSchedule,
)


def test_high_precision_scale_reconstruction_is_independent_and_consistent():
    schedule = KokunoOuterReservedPatchSchedule()
    reference = _reference_scales(schedule)
    public = schedule.log_scale_report()

    keys = (
        "T_d",
        "log_P_star",
        "T_w",
        "log_X_R",
        "log_X_w",
        "log_P1",
        "log_e_w",
        "log_c_patch",
        "log_X_star",
        "log_e_star",
    )
    assert max(abs(reference[key] - float(public[key])) for key in keys) <= SCALE_LOG_ABSOLUTE_TOLERANCE
    assert reference["I2_low"] < reference["log_X_star"] < reference["I2_high"]


def test_bisection_inverts_source_relation_without_production_fixed_point():
    h = 0.006
    for t in (0.375, 0.5, 0.625):
        for eta in (-0.53, -0.17, 0.0, 0.31, 0.54):
            q_expected = (1.0 - t) / (1.0 - eta * eta)
            z = q_expected ** (0.5 - h) * eta
            q = _q_bisection(z, t, h)
            relation = q - z * z * q ** (2.0 * h) - (1.0 - t)
            assert q == pytest.approx(q_expected, rel=3.0e-14, abs=3.0e-14)
            assert abs(relation) / max(q, 1.0 - t) <= Q_RELATION_RELATIVE_TOLERANCE


def test_black_box_patch_audit_passes_and_detects_radial_mutation():
    schedule = KokunoOuterReservedPatchSchedule()
    result = audit_schedule(schedule, seed=9173061, samples=9)

    assert result["passed"] is True
    assert result["max_public_vs_independent_log_scale_abs_error"] <= SCALE_LOG_ABSOLUTE_TOLERANCE
    assert result["max_public_vs_independent_velocity_relative_error"] <= VELOCITY_RELATIVE_TOLERANCE
    assert result["max_independent_q_relation_relative_residual"] <= Q_RELATION_RELATIVE_TOLERANCE
    assert result["radial_scale_law_relative_error"] <= SCALE_LAW_RELATIVE_TOLERANCE
    finest = result["normalized_divergence_by_relative_step"]["2.0e-04"]["max"]
    assert finest <= LOCAL_NORMALIZED_DIVERGENCE_TOLERANCE
    assert result["mutation_finest_normalized_divergence_min"] >= MUTATION_MIN_NORMALIZED_DIVERGENCE
    assert result["public_speed_min"] > 0.0
    assert np.isfinite(result["public_speed_max"])


def test_full_independent_matrix_is_truth_bounded_and_not_a_pde_gate():
    report = run_audit(seed=9173061, samples_per_case=6)

    assert report["structural_preflight_passed"] is True
    assert report["parameter_case_count"] == 3
    assert all(case["passed"] for case in report["cases"])
    assert report["formal_project_gate"]["normalized_full_momentum_threshold"] == pytest.approx(1.0e-3)
    assert report["formal_project_gate"]["divergence_max_threshold"] == pytest.approx(1.0e-5)
    assert report["formal_project_gate"]["assessed"] is False
    assert report["st006_cross_route_boundary"]["directly_comparable"] is False
    truth = report["truth_boundary"]
    assert truth["training_loss_used"] is False
    assert truth["candidate_internal_derivatives_used"] is False
    assert truth["free_residual_defined_forcing_used"] is False
    assert truth["complete_kokuno_composite_velocity"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_bad_sample_count_fails_closed():
    with pytest.raises(ValueError):
        audit_schedule(KokunoOuterReservedPatchSchedule(), seed=1, samples=5)
