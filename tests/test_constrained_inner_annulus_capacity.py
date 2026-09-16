from pathlib import Path

import numpy as np

from openai_ns_reconstruction.constrained_inner_annulus_capacity import (
    FROZEN_COLLAR_POINT,
    INNER_FULL_RADIUS,
    INNER_ZERO_RADIUS,
    OLD_CORRECTION_ZERO_RADIUS,
    default_summary,
    guarded_swirl_basis,
    moment_cancelled_basis,
    moment_ratios,
    reach_rank,
)
from openai_ns_reconstruction.constrained_momentum_budget import angular_moment
from openai_ns_reconstruction.constrained_quintic_swirl import QuinticSwirlCandidate


CANDIDATE = Path("artifacts/constrained/quintic_swirl/candidate.json")


def test_new_modes_reach_exact_old_blind_region_but_leave_core_open():
    candidate = QuinticSwirlCandidate.load(CANDIDATE)
    h = 0.005
    p = FROZEN_COLLAR_POINT[0]
    stencil = [p]
    for axis in range(3):
        direction = np.eye(3)[axis]
        stencil.extend(p + k * h * direction for k in (-2, -1, 1, 2))
    stencil = np.asarray(stencil)

    # The current quintic corrections and all their fourth-order FD stencil
    # values are exactly zero here, so their coefficient Jacobian cannot repair
    # the recorded frozen-collar residual.
    old = candidate.correction_basis(stencil, 0.75)
    assert np.max(np.abs(old)) == 0.0
    assert np.max(np.hypot(stencil[:, 0], stencil[:, 1])) < OLD_CORRECTION_ZERO_RADIUS

    ratios = moment_ratios(candidate.parent, 48)
    new = moment_cancelled_basis(candidate.parent, FROZEN_COLLAR_POINT, 0.75, ratios=ratios)
    raw = guarded_swirl_basis(FROZEN_COLLAR_POINT)
    np.testing.assert_array_equal(new, raw)
    assert abs(new[0, 1, 0]) > 0.3
    assert abs(new[0, 1, 1]) > 0.3

    times = np.linspace(0.25, 0.75, 21)
    tau = 1.0 - times
    core = np.column_stack((0.1 * np.sqrt(tau), np.zeros(len(times)), 0.1 * tau**0.495))
    assert np.max(np.hypot(core[:, 0], core[:, 1])) < INNER_ZERO_RADIUS
    assert np.max(np.abs(moment_cancelled_basis(candidate.parent, core, times, ratios=ratios))) == 0.0

    exterior = np.array([[2.0, 0.0, 0.0], [0.4, 0.0, 2.0], [2.2, 0.0, 1.0]])
    assert np.max(np.abs(guarded_swirl_basis(exterior))) == 0.0
    assert INNER_ZERO_RADIUS < INNER_FULL_RADIUS < OLD_CORRECTION_ZERO_RADIUS


def test_two_mode_growth_increases_inner_annulus_reach_rank_with_mild_conditioning():
    one = reach_rank(mode_count=1)
    two = reach_rank(mode_count=2)
    assert one["rank"] == 1
    assert one["condition_number"] == 1.0
    assert two["rank"] == 2
    assert two["condition_number"] < 5.0

    summary = default_summary()
    assert summary["rank_by_basis_dimension"]["2"]["rank"] == 2
    assert min(abs(v) for v in summary["frozen_point_y_components"]) > 0.3


def test_modes_are_pure_swirl_divergence_free_and_moment_cancel_numerically():
    candidate = QuinticSwirlCandidate.load(CANDIDATE)
    ratios = moment_ratios(candidate.parent, 48)

    rng = np.random.default_rng(7702)
    radius = rng.uniform(0.11, 0.34, 32)
    theta = rng.uniform(0.0, 2.0 * np.pi, 32)
    z = rng.uniform(0.3, 1.2, 32)
    points = np.column_stack((radius * np.cos(theta), radius * np.sin(theta), z))
    step = 2.0e-6
    for mode in range(2):
        divergence = np.zeros(len(points))
        for axis in range(3):
            direction = np.eye(3)[axis] * step
            plus = guarded_swirl_basis(points + direction)[..., :, mode]
            minus = guarded_swirl_basis(points - direction)[..., :, mode]
            divergence += (plus[:, axis] - minus[:, axis]) / (2.0 * step)
        assert np.max(np.abs(divergence)) < 2.0e-7

        def raw_field(p, t, selected=mode):
            del t
            return guarded_swirl_basis(p)[..., :, selected]

        def cancelled_field(p, t, selected=mode):
            return moment_cancelled_basis(candidate.parent, p, t, ratios=ratios)[..., :, selected]

        raw_moment = angular_moment(raw_field, 0.5, 48)
        cancelled_moment = angular_moment(cancelled_field, 0.5, 48)
        assert abs(cancelled_moment) <= 1.0e-9 * max(1.0, abs(raw_moment))


def test_invalid_shapes_and_basis_dimensions_fail_closed():
    with np.testing.assert_raises(ValueError):
        guarded_swirl_basis(np.zeros((4, 2)))
    with np.testing.assert_raises(ValueError):
        guarded_swirl_basis(np.array([[np.nan, 0.0, 0.0]]))
    with np.testing.assert_raises(ValueError):
        reach_rank(mode_count=0)
    with np.testing.assert_raises(ValueError):
        reach_rank(mode_count=3)
    with np.testing.assert_raises(TypeError):
        reach_rank(mode_count=True)
