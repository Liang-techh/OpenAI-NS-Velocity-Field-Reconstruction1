from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_strict_inner_differentiable_candidate_independent_audit import (
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_NORMALIZED_GATE,
    RICHARDSON_TIME_STEPS,
    component_error_metrics,
    frozen_gate_failures,
    public_contract,
    refinement_ratio,
    richardson_time_level,
)


def _poly_velocity(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack(
        (
            x + 2.0 * t + 0.5 * t**2 - 0.25 * t**3,
            y - 0.7 * t + 0.3 * t**2 + 0.2 * t**4,
            z + 0.4 * t - 0.1 * t**2 + 0.05 * t**3,
        ),
        axis=-1,
    )


def _poly_velocity_dt(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    del x, y, z
    return np.stack(
        (
            2.0 + t - 0.75 * t**2,
            -0.7 + 0.6 * t + 0.8 * t**3,
            0.4 - 0.2 * t + 0.15 * t**2,
        ),
        axis=-1,
    )


def _fixture_cloud():
    x = np.asarray((0.12, -0.21, 0.08, -0.17, 0.03))
    y = np.asarray((-0.09, 0.14, -0.18, 0.06, 0.22))
    z = np.asarray((0.04, -0.16, 0.19, -0.07, 0.11))
    t = np.asarray((0.41, 0.46, 0.50, 0.55, 0.59))
    return x, y, z, t


def test_richardson_time_derivative_locks_fixed_cartesian_orientation():
    x, y, z, t = _fixture_cloud()
    expected = _poly_velocity_dt(x, y, z, t)
    levels = tuple(
        richardson_time_level(_poly_velocity, x, y, z, t, step=h)
        for h in RICHARDSON_TIME_STEPS
    )
    for level in levels:
        np.testing.assert_allclose(level.velocity_dt, expected, atol=3.0e-12, rtol=3.0e-12)
    assert frozen_gate_failures(expected, levels) == []


def test_frozen_mutations_are_detected_without_moving_gates():
    x, y, z, t = _fixture_cloud()
    production = _poly_velocity_dt(x, y, z, t)
    levels = tuple(
        richardson_time_level(_poly_velocity, x, y, z, t, step=h)
        for h in RICHARDSON_TIME_STEPS
    )

    failures_scale = frozen_gate_failures(0.99 * production, levels)
    assert "fine_relative_rms" in failures_scale
    assert "fine_relative_sampled_max" in failures_scale

    flipped = production.copy()
    flipped[..., 2] *= -1.0
    failures_flip = frozen_gate_failures(flipped, levels)
    assert "fine_relative_rms" in failures_flip
    assert "fine_relative_sampled_max" in failures_flip


def test_component_metrics_report_all_cartesian_time_components():
    x, y, z, t = _fixture_cloud()
    expected = _poly_velocity_dt(x, y, z, t)
    observed = expected.copy()
    observed[:, 1] += 1.0e-6
    metrics = component_error_metrics(expected, observed)
    assert [m["component"] for m in metrics] == ["u_t", "v_t", "w_t"]
    assert metrics[0]["absolute_rms"] == pytest.approx(0.0)
    assert metrics[1]["absolute_rms"] == pytest.approx(1.0e-6)
    assert metrics[2]["absolute_rms"] == pytest.approx(0.0)


def test_refinement_requires_ordered_three_level_ladder():
    x, y, z, t = _fixture_cloud()
    levels = tuple(
        richardson_time_level(_poly_velocity, x, y, z, t, step=h)
        for h in RICHARDSON_TIME_STEPS
    )
    coarse_to_medium, medium_to_fine, ratio = refinement_ratio(*levels)
    assert coarse_to_medium >= 0.0
    assert medium_to_fine >= 0.0
    assert ratio >= 0.0
    with pytest.raises(ValueError, match="ordered"):
        refinement_ratio(levels[2], levels[1], levels[0])


def test_public_contract_keeps_final_pde_boundary_fail_closed():
    contract = public_contract()
    assert contract["heldout_seed"] == 9173531
    assert contract["richardson_time_steps"] == list(RICHARDSON_TIME_STEPS)
    assert contract["frozen_gates"]["final_momentum_normalized_max_l2"] == FINAL_MOMENTUM_NORMALIZED_GATE == 1.0e-3
    assert contract["frozen_gates"]["final_divergence_max_l2"] == FINAL_DIVERGENCE_GATE == 1.0e-5
    assert contract["independent_reference"]["numerical_reference_uses"] == "reloaded public velocity(x,y,z,t) only"
    assert contract["independent_reference"]["uses_agent1_analytic_time_derivative_as_reference"] is False
    assert contract["independent_reference"]["uses_agent2_fd4_verifier_as_reference"] is False
    truth = contract["truth_boundary"]
    assert truth["strict_inner_velocity_dt_independently_assessed"] is True
    assert truth["matched_pressure_included"] is False
    assert truth["restricted_forcing_included"] is False
    assert truth["correction_velocity_included"] is False
    assert truth["heldout_ns_momentum_residual_assessed"] is False
    assert truth["pde_validated"] is False
