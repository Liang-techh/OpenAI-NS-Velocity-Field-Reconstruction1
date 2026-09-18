import numpy as np

from openai_ns_reconstruction.kokuno_agent4_signed_complete_curl_covariance_audit import (
    MAX_DUPLICATED_SIGN_RANK_RATIO,
    MAX_RESOLUTION_RELATIVE_DRIFT,
    MIN_PHYSICAL_RANK_RATIO,
    MIN_VELOCITY_RMS,
    RADII,
    RESOLUTIONS,
    SEED,
    TARGET_STEP,
    _covariance_from_public_cartesian,
    build_report,
)


def test_cartesian_to_cylindrical_covariance_oracle_is_independent_of_theta_rotation():
    theta = np.linspace(0.11, 5.93, 96).reshape(12, 8)
    ur = 0.7 + 0.13 * np.cos(3.0 * theta)
    ut = -0.4 + 0.09 * np.sin(2.0 * theta)
    uz = 0.2 + 0.07 * np.cos(theta)
    c = np.cos(theta)
    s = np.sin(theta)
    cart = np.stack((ur * c - ut * s, ur * s + ut * c, uz), axis=-1)

    covariance, rms = _covariance_from_public_cartesian(cart, theta)
    expected = np.array([np.mean(ur * ut), np.mean(ur * uz)])
    np.testing.assert_allclose(covariance, expected, rtol=2e-15, atol=2e-15)
    assert rms > 0.0


def test_signed_complete_curl_covariance_audit_emits_frozen_fail_closed_receipt():
    report = build_report()
    assert report["seed"] == SEED == 9173191
    assert tuple(report["resolutions"]) == RESOLUTIONS == (12, 24, 48)
    assert tuple(report["radii"]) == RADII == (0.08, 0.35, 0.90)
    assert report["target_step"] == TARGET_STEP == 2e-4
    assert report["guards_frozen_before_execution"] == {
        "minimum_physical_rank_ratio": MIN_PHYSICAL_RANK_RATIO,
        "maximum_resolution_relative_drift": MAX_RESOLUTION_RELATIVE_DRIFT,
        "minimum_velocity_rms": MIN_VELOCITY_RMS,
        "maximum_duplicated_sign_mutation_rank_ratio": MAX_DUPLICATED_SIGN_RANK_RATIO,
    }

    oracle = report["oracle"]
    assert oracle["candidate_consumption"] == "public velocity_physical_cartesian_total values only"
    assert oracle["complete_curl_helper_reused_by_oracle"] is False
    assert oracle["returned_candidate_amplitudes_or_tangents_read_by_oracle"] is False

    metrics = report["metrics"]
    for value in metrics.values():
        assert np.isfinite(value)
        assert value >= 0.0
    assert len(report["cases"]) == 3
    assert all(len(case["resolutions"]) == 3 for case in report["cases"])
    assert all(len(case["resolution_stability"]) == 2 for case in report["cases"])

    # The scientific verdict is deliberately data-driven.  A REJECT is still a
    # valid audit artifact and must not be hidden by changing these guards.
    assert isinstance(report["local_physical_covariance_preflight_passed"], bool)

    truth = report["truth_boundary"]
    assert truth["supplied_manufactured_physical_signed_covariance_rank_two_assessed"] is True
    for key in (
        "actual_source_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_source_signed_rectangles_or_modes_bound",
        "actual_source_partition_labels_instantiated",
        "public_source_bound_xyz_t_oscillatory_velocity_ready",
        "actual_source_physical_covariance_rank_two_assessed",
        "genuinely_independent_second_covariance_column_ready",
        "correction_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert truth[key] is False

    gates = report["formal_project_gates_unchanged"]
    assert gates == {
        "normalized_momentum_max_and_L2": 1e-3,
        "divergence_max_and_L2": 1e-5,
        "assessed_here": False,
    }
