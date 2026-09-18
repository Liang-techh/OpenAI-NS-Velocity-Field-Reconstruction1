import numpy as np
import pytest

from openai_ns_reconstruction.constrained_streamline_seed_plan import (
    plan_activity_aware_streamline_seeds,
)


def helical_velocity(points, time):
    p = np.asarray(points, dtype=float)
    return np.column_stack((-1.4 * p[:, 1], 1.4 * p[:, 0], 0.25 + 0.2 * p[:, 2]))


def test_exact_200_deterministic_well_spaced_seeds_inside_declared_support():
    kwargs = dict(
        time=0.5,
        count=200,
        support_radius=2.0,
        support_half_height=2.0,
        visualization_seed=914113,
        relative_activity_floor=0.02,
    )
    a = plan_activity_aware_streamline_seeds(helical_velocity, **kwargs)
    b = plan_activity_aware_streamline_seeds(helical_velocity, **kwargs)

    assert a.seeds.shape == (200, 3)
    assert a.speeds.shape == (200,)
    np.testing.assert_array_equal(a.seeds, b.seeds)
    np.testing.assert_array_equal(a.speeds, b.speeds)
    assert not a.seeds.flags.writeable
    assert not a.speeds.flags.writeable

    r = np.hypot(a.seeds[:, 0], a.seeds[:, 1])
    assert np.max(r) <= 2.0 * 0.95 + 1e-12
    assert np.max(np.abs(a.seeds[:, 2])) <= 2.0 * 0.95 + 1e-12
    assert a.min_normalized_separation > 0.05
    assert np.min(a.speeds) >= a.activity_floor
    assert a.claim_scope == "visualization_seed_selection_only"
    assert not a.visualization_ready
    assert not a.visual_correspondence_verified
    assert not a.pde_validated
    assert not a.paper_exact
    assert not a.openai_field_identified
    assert not a.blowup_proved


def test_activity_floor_rejects_slow_center_without_changing_truth_boundary():
    def pure_swirl(points, time):
        p = np.asarray(points, dtype=float)
        return np.column_stack((-p[:, 1], p[:, 0], np.zeros(len(p))))

    plan = plan_activity_aware_streamline_seeds(
        pure_swirl,
        0.5,
        count=64,
        support_radius=1.0,
        support_half_height=1.0,
        visualization_seed=914117,
        relative_activity_floor=0.65,
        pool_factor=8,
    )
    r = np.hypot(plan.seeds[:, 0], plan.seeds[:, 1])
    assert np.min(plan.speeds) >= plan.activity_floor
    assert np.min(r) >= plan.activity_floor - 1e-12
    assert plan.pde_validated is False


def test_fail_closed_zero_malformed_nonfinite_and_bad_settings():
    def zero(points, time):
        return np.zeros_like(points, dtype=float)

    with pytest.raises(ValueError, match="inactive"):
        plan_activity_aware_streamline_seeds(
            zero, 0.5, count=16, support_radius=1.0, support_half_height=1.0
        )

    def bad_shape(points, time):
        return np.zeros((len(points), 2))

    with pytest.raises(ValueError, match="shape"):
        plan_activity_aware_streamline_seeds(
            bad_shape, 0.5, count=16, support_radius=1.0, support_half_height=1.0
        )

    def nan_velocity(points, time):
        out = np.zeros_like(points, dtype=float)
        out[0, 0] = np.nan
        return out

    with pytest.raises(ValueError, match="non-finite"):
        plan_activity_aware_streamline_seeds(
            nan_velocity, 0.5, count=16, support_radius=1.0, support_half_height=1.0
        )

    with pytest.raises(ValueError):
        plan_activity_aware_streamline_seeds(
            helical_velocity, np.nan, count=16, support_radius=1.0, support_half_height=1.0
        )
    with pytest.raises(ValueError):
        plan_activity_aware_streamline_seeds(
            helical_velocity, 0.5, count=0, support_radius=1.0, support_half_height=1.0
        )
    with pytest.raises(ValueError):
        plan_activity_aware_streamline_seeds(
            helical_velocity, 0.5, count=16, support_radius=-1.0, support_half_height=1.0
        )
