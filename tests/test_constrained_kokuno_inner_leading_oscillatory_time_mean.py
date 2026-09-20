from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_inner_leading_oscillatory_time_mean import (
    AGENT2_SOURCE_BLOB_SHA,
    ANGULAR_ORDERS,
    ExactAgent2InnerLeadingOscillatoryTimeBackend,
    _materialize_with_evaluator,
    materialize_agent2_inner_leading_oscillatory_time_mean,
    public_contract,
    truth_boundary,
)


INNER_MEAN = np.asarray([0.09, -0.03, 0.04], dtype=float)
OSC_MEAN = np.asarray([-0.02, 0.05, 0.07], dtype=float)
TOTAL_MEAN = INNER_MEAN + OSC_MEAN


def _cylindrical_vector(theta: np.ndarray, mean: np.ndarray, phase: float) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)
    radial = mean[0] + 0.03 * np.cos(2.0 * theta + phase)
    tangential = mean[1] + 0.02 * np.sin(2.0 * theta - phase)
    axial = mean[2] + 0.015 * np.cos(2.0 * theta - 0.5 * phase)
    return np.stack(
        [radial * c - tangential * s, radial * s + tangential * c, axial],
        axis=-1,
    )


class DummyInnerBackend:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(z, dtype=float), np.asarray(t, dtype=float),
        )
        return np.stack([0.1 + 0.0 * x, 0.2 + 0.0 * y, 0.3 + 0.0 * z], axis=-1)

    def velocity_dt(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(z, dtype=float), np.asarray(t, dtype=float),
        )
        theta = np.arctan2(y, x)
        return _cylindrical_vector(theta, INNER_MEAN, 0.3)


def _good_evaluator(inner, x, y, z, t):
    dt_inner = np.asarray(inner.velocity_dt(x, y, z, t), dtype=float)
    theta = np.arctan2(np.asarray(y, dtype=float), np.asarray(x, dtype=float))
    dt_osc = _cylindrical_vector(theta, OSC_MEAN, -0.4)
    return SimpleNamespace(
        inner_leading_velocity_dt=dt_inner,
        oscillatory_velocity_dt=dt_osc,
        inner_plus_oscillatory_velocity_dt=dt_inner + dt_osc,
    )


def _witness(evaluator=_good_evaluator):
    return _materialize_with_evaluator(
        evaluator,
        DummyInnerBackend(),
        np.asarray([0.55, 0.9, 1.25], dtype=float),
        np.asarray([-0.25, 0.0, 0.3], dtype=float),
        np.asarray([0.46, 0.50, 0.54], dtype=float),
        backend=None,
        backend_kind="test_fixture",
    )


def test_rotating_cylindrical_projection_recovers_registered_time_means():
    witness = _witness()
    expected_inner = np.broadcast_to(INNER_MEAN, (3, 3))
    expected_osc = np.broadcast_to(OSC_MEAN, (3, 3))
    expected_total = np.broadcast_to(TOTAL_MEAN, (3, 3))

    np.testing.assert_allclose(
        witness.mean_inner_leading_velocity_dt_cylindrical,
        expected_inner,
        rtol=0.0,
        atol=2.0e-14,
    )
    np.testing.assert_allclose(
        witness.mean_oscillatory_velocity_dt_cylindrical,
        expected_osc,
        rtol=0.0,
        atol=2.0e-14,
    )
    np.testing.assert_allclose(
        witness.mean_inner_plus_oscillatory_velocity_dt_cylindrical,
        expected_total,
        rtol=0.0,
        atol=2.0e-14,
    )
    assert witness.additive_mean_closure_absolute_max <= 2.0e-14
    assert max(witness.successive_total_dt_mean_relative_differences) <= 2.0e-14
    assert witness.total_dt_mean_rms > 0.0
    assert witness.total_dt_mean_to_full_rms_ratio <= 1.0 + 1.0e-12


def test_public_contract_freezes_ladder_and_keeps_truth_boundary_fail_closed():
    contract = public_contract()
    assert contract["angular_orders"] == list(ANGULAR_ORDERS) == [32, 64, 128]
    assert contract["forbidden_public_inputs_present"] == []
    assert contract["final_gates_unchanged"] == {
        "normalized_momentum_max_l2": 1.0e-3,
        "normalized_divergence_max_l2": 1.0e-5,
    }
    assert set(inspect.signature(materialize_agent2_inner_leading_oscillatory_time_mean).parameters) == {
        "backend", "inner_leading_backend", "radius", "z", "t"
    }

    boundary = truth_boundary()
    for key in (
        "inner_leading_is_final_corrected_fixed_point",
        "outer_global_leading_join_materialized",
        "global_leading_velocity_materialized",
        "inner_leading_oscillatory_velocity_dt_mean_is_complete_ns_defect",
        "inner_leading_self_advection_included_here",
        "inner_leading_oscillatory_cross_advection_included_here",
        "pressure_gradient_included",
        "restricted_forcing_included",
        "real_agent3_delta_a_bound",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
        "paper_exact",
        "blowup_proved",
    ):
        assert boundary[key] is False


def test_public_materializer_requires_exact_backend_and_inner_velocity_dt():
    exact = ExactAgent2InnerLeadingOscillatoryTimeBackend(
        evaluator=_good_evaluator,
        source_file="fixture",
        source_blob_sha=AGENT2_SOURCE_BLOB_SHA,
    )

    class VelocityOnly:
        def velocity(self, x, y, z, t):
            shape = np.broadcast(
                np.asarray(x), np.asarray(y), np.asarray(z), np.asarray(t)
            ).shape
            return np.zeros(shape + (3,), dtype=float)

    with pytest.raises(TypeError, match="velocity_dt"):
        materialize_agent2_inner_leading_oscillatory_time_mean(
            exact, VelocityOnly(), 0.7, 0.0, 0.5
        )

    drifted = ExactAgent2InnerLeadingOscillatoryTimeBackend(
        evaluator=_good_evaluator,
        source_file="fixture",
        source_blob_sha="0" * 40,
    )
    with pytest.raises(ValueError, match="provenance drifted"):
        materialize_agent2_inner_leading_oscillatory_time_mean(
            drifted, DummyInnerBackend(), 0.7, 0.0, 0.5
        )


def test_payload_shape_nonfinite_and_additive_closure_fail_closed():
    def bad_shape(inner, x, y, z, t):
        good = _good_evaluator(inner, x, y, z, t)
        return SimpleNamespace(
            inner_leading_velocity_dt=good.inner_leading_velocity_dt[..., :2],
            oscillatory_velocity_dt=good.oscillatory_velocity_dt,
            inner_plus_oscillatory_velocity_dt=good.inner_plus_oscillatory_velocity_dt,
        )

    with pytest.raises(ValueError, match="unexpected ring shape"):
        _witness(bad_shape)

    def nonfinite(inner, x, y, z, t):
        good = _good_evaluator(inner, x, y, z, t)
        total = np.array(good.inner_plus_oscillatory_velocity_dt, copy=True)
        total[..., 0] = np.nan
        return SimpleNamespace(
            inner_leading_velocity_dt=good.inner_leading_velocity_dt,
            oscillatory_velocity_dt=good.oscillatory_velocity_dt,
            inner_plus_oscillatory_velocity_dt=total,
        )

    with pytest.raises(ValueError, match="non-finite"):
        _witness(nonfinite)

    def closure_drift(inner, x, y, z, t):
        good = _good_evaluator(inner, x, y, z, t)
        total = np.array(good.inner_plus_oscillatory_velocity_dt, copy=True)
        total[..., 2] += 1.0e-4
        return SimpleNamespace(
            inner_leading_velocity_dt=good.inner_leading_velocity_dt,
            oscillatory_velocity_dt=good.oscillatory_velocity_dt,
            inner_plus_oscillatory_velocity_dt=total,
        )

    with pytest.raises(ValueError, match="additive closure"):
        _witness(closure_drift)


def test_ring_domain_contract_is_inherited_fail_closed():
    with pytest.raises(ValueError, match="0 < r < 2"):
        _materialize_with_evaluator(
            _good_evaluator,
            DummyInnerBackend(),
            np.asarray([0.0]),
            np.asarray([0.0]),
            np.asarray([0.5]),
            backend=None,
            backend_kind="test_fixture",
        )
