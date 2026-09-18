import numpy as np

from openai_ns_reconstruction.kokuno_correction_cycle_gain_guard import (
    ComparableFullDomainResidual,
    ST006_MOMENTUM_MAX,
    ST006_VOLUME_L2,
    empirical_exponent_gain,
    evaluate_correction_cycle_gain_guard,
    source_exponent_gain_margins,
)


def _cycle_report(
    *,
    held_in_ratio=0.95,
    held_out_ratio=0.94,
    theta_ratio=0.97,
    pressure_ratio=1.0005,
    correction_rms=1.0e-3,
    divergence_max=1.0e-7,
):
    return {
        "inputs": {"surrogate_defect_used": False},
        "held_in": {"rms_ratio": held_in_ratio},
        "held_out": {
            "aggregate_mean_defect_rms_ratio": held_out_ratio,
            "aggregate_theta_rms_ratio": theta_ratio,
            "pressure_inclusive_raw_phase_mean_operator_rms_ratio": pressure_ratio,
            "signed_correction_rms": correction_rms,
        },
        "independent_signed_correction_divergence": {"max": divergence_max},
        "truth_boundary": {
            "real_candidate_defect_consumed": True,
            "finite_correction_cycle_run": True,
        },
    }


def test_source_gain_formulas_are_exponent_margins_not_ratios():
    receipt = source_exponent_gain_margins()
    np.testing.assert_allclose(receipt["wave_residual_exponent_gain"], 0.39999)
    np.testing.assert_allclose(receipt["full_tangential_mean_exponent_gain"], 0.17)
    np.testing.assert_allclose(receipt["defect_exponent_gain"], 0.89996)
    assert receipt["all_source_displayed_gains_exceed_required"]
    assert "not measured residual_before/residual_after" in receipt["interpretation"]


def test_empirical_exponent_conversion_requires_a_real_scale():
    q_scale = 100.0
    ratio = q_scale ** (-0.2)
    np.testing.assert_allclose(empirical_exponent_gain(1.0, ratio, q_scale), 0.2)
    with np.testing.assert_raises(ValueError):
        empirical_exponent_gain(1.0, ratio, 1.0)


def test_measured_worsening_rejects_cycle_despite_positive_source_margins():
    audit = evaluate_correction_cycle_gain_guard(
        _cycle_report(held_in_ratio=1.05, held_out_ratio=1.15, theta_ratio=1.08)
    )
    assert audit["source_exponent_gain_audit"]["all_source_displayed_gains_exceed_required"]
    assert not audit["accepted_for_next_cycle"]
    assert not audit["finite_cycle_checks"]["held_in_total_mean_defect_improves"]
    assert not audit["finite_cycle_checks"]["held_out_total_mean_defect_improves"]
    assert not audit["finite_cycle_checks"]["held_out_theta_mean_defect_nonworsening"]
    assert not audit["source_class_comparable_to_this_cycle"]
    assert audit["empirical_source_exponent_comparison"] is None
    assert not audit["st006_comparison"]["directly_comparable"]


def test_source_class_comparison_fails_closed_without_q_scale():
    with np.testing.assert_raises(ValueError):
        evaluate_correction_cycle_gain_guard(
            _cycle_report(), source_class_comparable=True
        )
    with np.testing.assert_raises(ValueError):
        evaluate_correction_cycle_gain_guard(
            _cycle_report(), source_class_comparable=False, source_q_scale=100.0
        )


def test_explicit_comparable_full_domain_receipt_reports_st006_ratios_only():
    full = ComparableFullDomainResidual(
        momentum_max=0.9 * ST006_MOMENTUM_MAX,
        volume_l2=0.8 * ST006_VOLUME_L2,
        seed=12345,
        point_count=4096,
        finest_spatial_step=0.005,
    )
    audit = evaluate_correction_cycle_gain_guard(
        _cycle_report(), full_domain_receipt=full
    )
    comparison = audit["st006_comparison"]
    assert comparison["directly_comparable"]
    np.testing.assert_allclose(comparison["momentum_max_ratio_to_st006"], 0.9)
    np.testing.assert_allclose(comparison["volume_l2_ratio_to_st006"], 0.8)
    assert comparison["beats_st006_on_both_reported_momentum_metrics"]
    assert not comparison["repository_level_pde_improvement_claimed"]


def test_guard_refuses_surrogate_or_unexecuted_cycle():
    surrogate = _cycle_report()
    surrogate["inputs"]["surrogate_defect_used"] = True
    with np.testing.assert_raises(ValueError):
        evaluate_correction_cycle_gain_guard(surrogate)

    unexecuted = _cycle_report()
    unexecuted["truth_boundary"]["finite_correction_cycle_run"] = False
    with np.testing.assert_raises(ValueError):
        evaluate_correction_cycle_gain_guard(unexecuted)
