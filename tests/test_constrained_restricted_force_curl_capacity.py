import numpy as np
import pytest

from openai_ns_reconstruction.constrained_eq45_candidate import Eq45VelocityCandidate
from openai_ns_reconstruction.constrained_restricted_force_curl_capacity import (
    DEFAULT_DERIVATIVE_STEPS,
    DEFAULT_FIT_SEED,
    DEFAULT_HOLDOUT_SEED,
    audit_default_eq45_restricted_force_capacity,
    audit_restricted_force_curl_capacity,
    deterministic_probe_cloud,
)


def test_eq45_restricted_force_fit_is_bounded_separate_and_low_capacity():
    candidate = Eq45VelocityCandidate.seed()
    report = audit_default_eq45_restricted_force_capacity(candidate)

    assert report.fit_seed == DEFAULT_FIT_SEED
    assert report.holdout_seed == DEFAULT_HOLDOUT_SEED
    assert report.fit_point_count == report.holdout_point_count == 16
    assert report.force_bounds == (0.0, 10.0)
    assert report.fit_and_holdout_separate is True
    assert report.velocity_changed is False
    assert report.pressure_fitted is False
    assert report.forcing_fitted is True
    assert report.pde_validated is False

    fit = report.fit
    assert fit.design_rank == 2
    assert fit.design_condition == pytest.approx(22.7168663284, rel=3e-6)
    assert fit.singular_values == pytest.approx((7.54993736, 0.33234942), rel=3e-6)
    assert fit.a < 1e-8
    assert fit.c == pytest.approx(2.19114231, rel=3e-6)
    assert fit.active_mask == (-1, 0)
    assert fit.nonlinear_function_evaluations == 0
    assert fit.curl_rms_before == pytest.approx(17.1276738407, rel=3e-6)
    assert fit.curl_rms_after == pytest.approx(16.6208546589, rel=3e-6)
    assert fit.recoverable_rms_fraction == pytest.approx(0.0295906605, rel=3e-6)

    assert tuple(level.spatial_step for level in report.holdout_levels) == DEFAULT_DERIVATIVE_STEPS
    expected = (
        (23.7703689102, 23.3515394560, 0.0176198130),
        (24.0659630536, 23.6514378128, 0.0172245441),
        (24.1409212848, 23.7274736268, 0.0171264242),
    )
    for level, (before, after, fraction) in zip(report.holdout_levels, expected):
        assert level.curl_rms_before == pytest.approx(before, rel=3e-6)
        assert level.curl_rms_after == pytest.approx(after, rel=3e-6)
        assert level.recoverable_rms_fraction == pytest.approx(fraction, rel=3e-6)
        assert level.curl_rms_after > 20.0

    # The gain survives refinement but remains tiny: the registered two-parameter
    # force cannot explain more than 2% of held-out curl RMS on this probe contract.
    assert all(0.0 < level.recoverable_rms_fraction < 0.02 for level in report.holdout_levels)


def test_report_truth_boundary_does_not_promote_capacity_to_pde_validation():
    report = audit_default_eq45_restricted_force_capacity(Eq45VelocityCandidate.seed())
    payload = report.to_dict()

    assert payload["force_family"] == "preregistered_restricted_two_parameter_family"
    assert payload["pressure_role"] == "scalar_gradient_cannot_change_curl_obstruction"
    assert payload["velocity_changed"] is False
    assert payload["pressure_fitted"] is False
    assert payload["forcing_fitted"] is True
    assert payload["pde_validated"] is False
    assert "capacity evidence only" in payload["interpretation"]


def test_fit_and_holdout_cannot_be_collapsed_or_constraint_drifted():
    candidate = Eq45VelocityCandidate.seed()
    points, times = deterministic_probe_cloud(DEFAULT_FIT_SEED, count=4)
    other_points, other_times = deterministic_probe_cloud(DEFAULT_HOLDOUT_SEED, count=4)

    with pytest.raises(ValueError, match="disjoint"):
        audit_restricted_force_curl_capacity(
            candidate,
            points,
            times,
            points.copy(),
            times.copy(),
            fit_seed=DEFAULT_FIT_SEED,
            holdout_seed=DEFAULT_HOLDOUT_SEED,
        )

    with pytest.raises(ValueError, match="preregistered"):
        audit_restricted_force_curl_capacity(
            candidate,
            points,
            times,
            other_points,
            other_times,
            nu=0.02,
        )

    with pytest.raises(ValueError, match="at least three"):
        audit_restricted_force_curl_capacity(
            candidate,
            points,
            times,
            other_points,
            other_times,
            derivative_steps=(0.02, 0.01),
        )

    with pytest.raises(ValueError, match="distinct"):
        audit_restricted_force_curl_capacity(
            candidate,
            points,
            times,
            other_points,
            other_times,
            fit_seed=17,
            holdout_seed=17,
        )


def test_probe_generation_is_reproducible_and_independent_by_seed():
    first_points, first_times = deterministic_probe_cloud(DEFAULT_FIT_SEED)
    replay_points, replay_times = deterministic_probe_cloud(DEFAULT_FIT_SEED)
    holdout_points, holdout_times = deterministic_probe_cloud(DEFAULT_HOLDOUT_SEED)

    assert np.array_equal(first_points, replay_points)
    assert np.array_equal(first_times, replay_times)
    assert not np.array_equal(first_points, holdout_points)
    assert not np.array_equal(first_times, holdout_times)
    assert np.max(np.abs(first_points)) <= 0.25
    assert np.all((first_times >= 0.38) & (first_times <= 0.62))
