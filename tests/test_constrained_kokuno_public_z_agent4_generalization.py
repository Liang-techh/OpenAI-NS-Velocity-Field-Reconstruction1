from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_public_z_agent4_generalization import (
    AUDITED_AGENT2_HEAD,
    FD6_STEPS,
    GUARDS,
    PARENT_AGENT4_HEAD,
    SEED,
    _fd6_jacobian,
    _heldout_stratified_points,
    run_audit,
)


def test_fd6_jacobian_is_independently_calibrated_on_degree_six_polynomial() -> None:
    points = np.asarray(
        [
            [0.21, -0.17, 0.33],
            [-0.38, 0.29, -0.24],
            [0.47, 0.13, -0.19],
        ],
        dtype=float,
    )
    times = np.asarray([0.31, 0.49, 0.67], dtype=float)

    def field(x, y, z, t):
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        z = np.asarray(z, dtype=float)
        t = np.asarray(t, dtype=float)
        return np.stack(
            (
                x**6 + 2.0 * y + 0.5 * t,
                y**6 - 3.0 * z + 0.25 * t,
                z**6 + 4.0 * x - 0.75 * t,
            ),
            axis=-1,
        )

    jac = _fd6_jacobian(points, times, field, 0.017)
    expected = np.zeros_like(jac)
    expected[:, 0, 0] = 6.0 * points[:, 0] ** 5
    expected[:, 0, 1] = 2.0
    expected[:, 1, 1] = 6.0 * points[:, 1] ** 5
    expected[:, 1, 2] = -3.0
    expected[:, 2, 0] = 4.0
    expected[:, 2, 2] = 6.0 * points[:, 2] ** 5
    np.testing.assert_allclose(jac, expected, rtol=2.0e-10, atol=2.0e-11)


def test_stratified_cloud_is_fresh_and_stencils_remain_inside_support() -> None:
    points, times, labels = _heldout_stratified_points()
    assert SEED == 9173251
    assert points.shape == (48, 3)
    assert times.shape == (48,)
    assert set(labels.tolist()) == {"core", "radial_collar", "axial_collar"}
    assert all(int(np.sum(labels == name)) == 16 for name in set(labels.tolist()))

    max_probe = 3.0 * max(FD6_STEPS)
    radius = np.hypot(points[:, 0], points[:, 1])
    radial = labels == "radial_collar"
    axial = labels == "axial_collar"
    assert np.min(radius[radial]) - max_probe > 0.15
    assert np.max(radius[radial]) + max_probe < 1.35
    assert np.max(np.abs(points[axial, 2])) + max_probe < 2.0
    assert np.all((times >= 0.25) & (times <= 0.75))


def test_report_is_fail_closed_and_scientific_verdict_is_recomputed() -> None:
    report = run_audit()
    assert report["audited_agent2_head"] == AUDITED_AGENT2_HEAD
    assert report["parent_agent4_head"] == PARENT_AGENT4_HEAD
    assert report["fd6_steps"] == list(FD6_STEPS)
    assert report["guards_frozen_before_actions"] == GUARDS
    assert report["point_count"] == 48

    finest = report["divergence"]["rows"][-1]
    refinements = report["divergence"]["refinement_ratios"]
    strata = report["divergence"]["finest_stratum_relative_rms"]
    support = report["support_exterior"]
    perturb = report["parameter_perturbations"]
    expected_pass = bool(
        finest["relative_rms"] <= GUARDS["finest_relative_divergence_rms"]
        and finest["relative_max"] <= GUARDS["finest_relative_divergence_max"]
        and min(refinements) >= GUARDS["minimum_divergence_refinement_ratio"]
        and max(strata.values()) <= GUARDS["maximum_stratum_finest_relative_divergence_rms"]
        and support["absolute_max"] <= GUARDS["support_exterior_absolute_max"]
        and report["heldout_velocity_rms"] >= GUARDS["minimum_nontrivial_velocity_rms"]
        and min(perturb.values()) >= GUARDS["minimum_parameter_perturbation_relative_change"]
        and report["divergence"]["mutation_relative_rms"]
        >= GUARDS["minimum_divergence_mutation_relative_rms"]
    )
    assert report["public_z_generalization_passed"] is expected_pass

    truth = report["truth_boundary"]
    assert truth["different_derivative_order_from_parent_agent4_fd4"] is True
    assert truth["candidate_internal_derivative_helper_used"] is False
    assert truth["pressure_or_forcing_fitted"] is False
    assert truth["formal_full_domain_pde_gate_assessed"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
    assert truth["formal_momentum_normalized_max_l2_gate"] == 1.0e-3
    assert truth["formal_divergence_max_l2_gate"] == 1.0e-5
