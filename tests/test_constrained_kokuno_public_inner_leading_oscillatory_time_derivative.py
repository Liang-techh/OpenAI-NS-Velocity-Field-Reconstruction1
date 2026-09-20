from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_time_derivative import (
    AGENT1_HEAD,
    AGENT1_PR,
    PARENT_AGENT2_HEAD,
    PARENT_AGENT2_PR,
    evaluate_inner_leading_oscillatory_time_derivative,
    public_contract,
)


class ManufacturedInnerBackend:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(z, dtype=float), np.asarray(t, dtype=float),
        )
        return np.stack(
            (x * t + 0.25 * y, y * t * t - 0.5 * z, z * np.sin(t) + 0.2 * x),
            axis=-1,
        )

    def velocity_dt(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(z, dtype=float), np.asarray(t, dtype=float),
        )
        return np.stack((x, 2.0 * y * t, z * np.cos(t)), axis=-1)


def test_typed_composite_and_derivative_close_exactly():
    backend = ManufacturedInnerBackend()
    x = np.asarray((0.31, -0.47, 0.66, -0.24))
    y = np.asarray((0.42, 0.28, -0.36, -0.51))
    z = np.asarray((-0.72, 0.41, 0.63, -0.18))
    t = np.asarray((0.44, 0.48, 0.52, 0.56))
    result = evaluate_inner_leading_oscillatory_time_derivative(backend, x, y, z, t)

    np.testing.assert_array_equal(
        result.inner_plus_oscillatory_velocity,
        result.inner_leading_velocity + result.oscillatory_velocity,
    )
    np.testing.assert_array_equal(
        result.inner_plus_oscillatory_velocity_dt,
        result.inner_leading_velocity_dt + result.oscillatory_velocity_dt,
    )
    assert result.inner_plus_oscillatory_velocity.shape == (4, 3)
    assert result.inner_plus_oscillatory_velocity_dt.shape == (4, 3)
    assert np.all(np.isfinite(result.inner_plus_oscillatory_velocity_dt))


def test_production_derivative_matches_independent_fd4_of_composite():
    backend = ManufacturedInnerBackend()
    x = np.asarray((0.35, -0.52, 0.71))
    y = np.asarray((0.44, 0.33, -0.29))
    z = np.asarray((-0.61, 0.37, 0.58))
    t = np.asarray((0.46, 0.50, 0.54))
    result = evaluate_inner_leading_oscillatory_time_derivative(backend, x, y, z, t)

    def total(tt):
        r = evaluate_inner_leading_oscillatory_time_derivative(backend, x, y, z, tt)
        return r.inner_plus_oscillatory_velocity

    h = 1.0e-3
    fd4 = (total(t - 2*h) - 8*total(t - h) + 8*total(t + h) - total(t + 2*h)) / (12*h)
    np.testing.assert_allclose(
        result.inner_plus_oscillatory_velocity_dt,
        fd4,
        rtol=2.0e-6,
        atol=2.0e-8,
    )


def test_backend_contract_fails_closed():
    class MissingDt:
        def velocity(self, x, y, z, t):
            shape = np.broadcast(np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(t)).shape
            return np.zeros(shape + (3,))

    with pytest.raises(TypeError, match="velocity_dt"):
        evaluate_inner_leading_oscillatory_time_derivative(MissingDt(), 0.2, 0.3, 0.1, 0.5)

    class WrongShape(ManufacturedInnerBackend):
        def velocity_dt(self, x, y, z, t):
            return np.zeros((2,))

    with pytest.raises(ValueError, match="shape"):
        evaluate_inner_leading_oscillatory_time_derivative(WrongShape(), 0.2, 0.3, 0.1, 0.5)


def test_public_truth_boundary_remains_inner_only():
    contract = public_contract()
    assert contract["parent_agent2_pr"] == PARENT_AGENT2_PR == 812
    assert contract["parent_agent2_head"] == PARENT_AGENT2_HEAD
    assert contract["agent1_inner_velocity_dt_pr"] == AGENT1_PR == 811
    assert contract["agent1_inner_velocity_dt_head"] == AGENT1_HEAD
    assert contract["forbidden_public_inputs_present"] == []

    signature = inspect.signature(evaluate_inner_leading_oscillatory_time_derivative)
    assert list(signature.parameters) == ["inner_leading_backend", "x", "y", "z", "t"]

    truth = contract["truth_boundary"]
    assert truth["inner_plus_oscillatory_velocity_dt_executable"] is True
    assert truth["production_composite_dt_uses_finite_difference"] is False
    for key in (
        "agent1_exact_head_ci_assumed_passed",
        "agent4_inner_field_independent_admission_assumed",
        "inner_leading_is_final_corrected_fixed_point",
        "global_leading_velocity_materialized",
        "outer_join_materialized",
        "matched_pressure_included",
        "restricted_forcing_included",
        "agent3_mean_radial_chain_reimplemented",
        "complete_kokuno_candidate_assembled",
        "heldout_ns_residual_assessed",
        "residual_reduction_claimed",
        "paper_exact",
        "pde_validated",
    ):
        assert truth[key] is False
