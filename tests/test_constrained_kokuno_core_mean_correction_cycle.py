import numpy as np

from openai_ns_reconstruction.kokuno_core_mean_correction_cycle import (
    AMPLITUDE_GRID,
    CORE_SAMPLE_RADIAL_BOUNDS,
    CORE_SAMPLE_Z_HALF_WIDTH,
    HELD_OUT_SEED,
    TRAIN_SEED,
    _assert_stencil_safe,
    _decision,
    _sample_core_points,
)
from openai_ns_reconstruction.kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate


def test_core_sampling_is_deterministic_disjoint_and_stencil_safe():
    train_a = _sample_core_points(seed=TRAIN_SEED, count=12)
    train_b = _sample_core_points(seed=TRAIN_SEED, count=12)
    held_out = _sample_core_points(seed=HELD_OUT_SEED, count=18)
    assert np.array_equal(train_a, train_b)
    assert not np.array_equal(train_a, held_out[:12])

    radius = np.hypot(train_a[:, 0], train_a[:, 1])
    assert np.all(radius >= CORE_SAMPLE_RADIAL_BOUNDS[0])
    assert np.all(radius <= CORE_SAMPLE_RADIAL_BOUNDS[1])
    assert np.all(np.abs(train_a[:, 2]) <= CORE_SAMPLE_Z_HALF_WIDTH)

    leading = KokunoLeadingCoreSeriesCandidate()
    _assert_stencil_safe(leading, train_a, (0.375, 0.5, 0.625), 0.005)
    _assert_stencil_safe(leading, held_out, (0.375, 0.5, 0.625), 0.005)


def test_core_cycle_keeps_pr232_amplitude_grid_without_widening():
    assert AMPLITUDE_GRID == (-0.04, -0.02, -0.01, -0.005, 0.005, 0.01, 0.02, 0.04)
    assert 0.0 not in AMPLITUDE_GRID


def test_cycle_decision_rejects_aggregate_gain_when_theta_worsens():
    accepted, checks = _decision(
        training_ratio=0.98,
        held_out_ratio=0.96,
        theta_ratio=1.02,
        pressure_operator_ratio=0.999,
        correction_divergence_max=1.0e-8,
        amplitude=0.005,
        stress_scale=1.0e-5,
    )
    assert not accepted
    assert checks["held_in_total_mean_defect_improves"]
    assert checks["held_out_total_mean_defect_improves"]
    assert not checks["held_out_theta_mean_defect_does_not_worsen"]


def test_cycle_decision_requires_nontrivial_bounded_divergence_correction():
    accepted, checks = _decision(
        training_ratio=0.98,
        held_out_ratio=0.97,
        theta_ratio=0.99,
        pressure_operator_ratio=1.0,
        correction_divergence_max=2.0e-5,
        amplitude=0.005,
        stress_scale=1.0e-5,
    )
    assert not accepted
    assert not checks["independent_correction_divergence_max_le_1e-5"]

    accepted_zero, zero_checks = _decision(
        training_ratio=0.98,
        held_out_ratio=0.97,
        theta_ratio=0.99,
        pressure_operator_ratio=1.0,
        correction_divergence_max=1.0e-8,
        amplitude=0.0,
        stress_scale=1.0e-5,
    )
    assert not accepted_zero
    assert not zero_checks["correction_nontrivial"]
