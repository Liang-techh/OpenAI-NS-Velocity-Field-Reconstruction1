from __future__ import annotations

import inspect
from types import SimpleNamespace

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_inner_leading_oscillatory_full_advection_mean import (
    AGENT2_INNER_LEADING_FD4_STEP,
    AGENT2_OSCILLATORY_FD6_STEP,
    AGENT2_SOURCE_BLOB_SHA,
    ANGULAR_ORDERS,
    ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend,
    _materialize_with_evaluator,
    materialize_agent2_inner_leading_oscillatory_full_advection_mean,
    public_contract,
    truth_boundary,
)


INNER_SELF_MEAN = np.asarray([0.16, -0.03, 0.12], dtype=float)
INNER_OSC_MEAN = np.asarray([0.21, -0.06, 0.14], dtype=float)
OSC_INNER_MEAN = np.asarray([-0.08, 0.02, 0.11], dtype=float)
OSC_SELF_MEAN = np.asarray([0.07, 0.04, -0.02], dtype=float)
MIXED_MEAN = INNER_OSC_MEAN + OSC_INNER_MEAN
TOTAL_MEAN = INNER_SELF_MEAN + MIXED_MEAN + OSC_SELF_MEAN


def _cylindrical_vector(theta: np.ndarray, mean: np.ndarray, phase: float) -> np.ndarray:
    theta = np.asarray(theta, dtype=float)
    c = np.cos(theta)
    s = np.sin(theta)
    radial = mean[0] + 0.025 * np.cos(2.0 * theta + phase)
    tangential = mean[1] + 0.018 * np.sin(2.0 * theta - phase)
    axial = mean[2] + 0.012 * np.cos(2.0 * theta + 0.5 * phase)
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
        return np.stack([0.2 + 0.0 * x, -0.1 + 0.0 * y, 0.3 + 0.0 * z], axis=-1)

    def self_advection(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(z, dtype=float), np.asarray(t, dtype=float),
        )
        theta = np.arctan2(y, x)
        return _cylindrical_vector(theta, INNER_SELF_MEAN, 0.2)


def _good_evaluator(
    inner,
    x,
    y,
    z,
    t,
    *,
    inner_leading_spatial_step,
    oscillatory_spatial_step,
):
    assert inner_leading_spatial_step == AGENT2_INNER_LEADING_FD4_STEP
    assert oscillatory_spatial_step == AGENT2_OSCILLATORY_FD6_STEP
    theta = np.arctan2(np.asarray(y, dtype=float), np.asarray(x, dtype=float))
    inner_self = np.asarray(inner.self_advection(x, y, z, t), dtype=float)
    inner_osc = _cylindrical_vector(theta, INNER_OSC_MEAN, 0.4)
    osc_inner = _cylindrical_vector(theta, OSC_INNER_MEAN, -0.5)
    osc_self = _cylindrical_vector(theta, OSC_SELF_MEAN, 0.8)
    mixed = inner_osc + osc_inner
    return SimpleNamespace(
        inner_leading_self_advection=inner_self,
        inner_advects_oscillation=inner_osc,
        oscillation_advects_inner=osc_inner,
        mixed_cross_advection=mixed,
        oscillatory_self_advection=osc_self,
        inner_plus_oscillatory_self_advection=inner_self + mixed + osc_self,
        inner_leading_spatial_step=inner_leading_spatial_step,
        oscillatory_spatial_step=oscillatory_spatial_step,
    )


def _witness(evaluator=_good_evaluator):
    return _materialize_with_evaluator(
        evaluator,
        DummyInnerBackend(),
        np.asarray([0.55, 0.90, 1.25], dtype=float),
        np.asarray([-0.25, 0.0, 0.30], dtype=float),
        np.asarray([0.46, 0.50, 0.54], dtype=float),
        backend=None,
        backend_kind="test_fixture",
    )


def test_rotating_projection_recovers_all_registered_full_advection_means():
    witness = _witness()
    expected = {
        "inner_self": np.broadcast_to(INNER_SELF_MEAN, (3, 3)),
        "inner_osc": np.broadcast_to(INNER_OSC_MEAN, (3, 3)),
        "osc_inner": np.broadcast_to(OSC_INNER_MEAN, (3, 3)),
        "mixed": np.broadcast_to(MIXED_MEAN, (3, 3)),
        "osc_self": np.broadcast_to(OSC_SELF_MEAN, (3, 3)),
        "total": np.broadcast_to(TOTAL_MEAN, (3, 3)),
    }
    observed = {
        "inner_self": witness.mean_inner_leading_self_advection_cylindrical,
        "inner_osc": witness.mean_inner_advects_oscillation_cylindrical,
        "osc_inner": witness.mean_oscillation_advects_inner_cylindrical,
        "mixed": witness.mean_mixed_cross_advection_cylindrical,
        "osc_self": witness.mean_oscillatory_self_advection_cylindrical,
        "total": witness.mean_inner_plus_oscillatory_full_advection_cylindrical,
    }
    for key in expected:
        np.testing.assert_allclose(observed[key], expected[key], rtol=0.0, atol=2.0e-14)

    assert witness.mixed_term_mean_closure_absolute_max <= 2.0e-14
    assert witness.full_decomposition_mean_closure_absolute_max <= 2.0e-14
    assert max(witness.successive_total_advection_mean_relative_differences) <= 2.0e-14
    assert witness.total_advection_mean_rms > 0.0
    assert witness.total_advection_mean_to_full_rms_ratio <= 1.0 + 1.0e-12


def test_public_contract_freezes_upstream_and_keeps_truth_boundary_fail_closed():
    contract = public_contract()
    assert contract["angular_orders"] == list(ANGULAR_ORDERS) == [32, 64, 128]
    assert contract["forbidden_public_inputs_present"] == []
    assert contract["agent2_source_blob_sha"] == AGENT2_SOURCE_BLOB_SHA
    assert contract["frozen_agent2_derivative_steps"] == {
        "inner_leading_fd4": 1.0e-3,
        "oscillatory_fd6": 1.0e-3,
    }
    assert contract["final_gates_unchanged"] == {
        "normalized_momentum_max_l2": 1.0e-3,
        "normalized_divergence_max_l2": 1.0e-5,
    }
    assert set(
        inspect.signature(
            materialize_agent2_inner_leading_oscillatory_full_advection_mean
        ).parameters
    ) == {"backend", "inner_leading_backend", "radius", "z", "t"}

    boundary = truth_boundary()
    for key in (
        "inner_leading_is_final_corrected_fixed_point",
        "outer_global_leading_join_materialized",
        "global_leading_velocity_materialized",
        "strict_inner_full_advection_mean_is_complete_ns_defect",
        "inner_plus_oscillatory_velocity_dt_mean_included_here",
        "base_viscous_term_included_here",
        "correction_transport_mean_included_here",
        "pressure_gradient_included",
        "restricted_forcing_included",
        "real_agent3_delta_a_bound",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert boundary[key] is False
    assert boundary["inner_leading_self_advection_included_here"] is True
    assert boundary["inner_leading_oscillatory_cross_advection_included_here"] is True
    assert boundary["oscillatory_self_advection_mean_included_here"] is True


def test_public_materializer_requires_exact_backend_and_self_advection():
    exact = ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend(
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

    with pytest.raises(TypeError, match="self_advection"):
        materialize_agent2_inner_leading_oscillatory_full_advection_mean(
            exact, VelocityOnly(), 0.7, 0.0, 0.5
        )

    drifted = ExactAgent2InnerLeadingOscillatoryFullAdvectionBackend(
        evaluator=_good_evaluator,
        source_file="fixture",
        source_blob_sha="0" * 40,
    )
    with pytest.raises(ValueError, match="provenance drifted"):
        materialize_agent2_inner_leading_oscillatory_full_advection_mean(
            drifted, DummyInnerBackend(), 0.7, 0.0, 0.5
        )


def test_shape_nonfinite_mixed_and_full_closure_fail_closed():
    def bad_shape(inner, x, y, z, t, **kwargs):
        good = _good_evaluator(inner, x, y, z, t, **kwargs)
        good.inner_leading_self_advection = good.inner_leading_self_advection[..., :2]
        return good

    with pytest.raises(ValueError, match="unexpected ring shape"):
        _witness(bad_shape)

    def nonfinite(inner, x, y, z, t, **kwargs):
        good = _good_evaluator(inner, x, y, z, t, **kwargs)
        total = np.array(good.inner_plus_oscillatory_self_advection, copy=True)
        total[..., 0] = np.nan
        good.inner_plus_oscillatory_self_advection = total
        return good

    with pytest.raises(ValueError, match="non-finite"):
        _witness(nonfinite)

    def mixed_drift(inner, x, y, z, t, **kwargs):
        good = _good_evaluator(inner, x, y, z, t, **kwargs)
        mixed = np.array(good.mixed_cross_advection, copy=True)
        mixed[..., 1] += 1.0e-4
        good.mixed_cross_advection = mixed
        good.inner_plus_oscillatory_self_advection = (
            good.inner_leading_self_advection + mixed + good.oscillatory_self_advection
        )
        return good

    with pytest.raises(ValueError, match="mixed-advection closure"):
        _witness(mixed_drift)

    def full_drift(inner, x, y, z, t, **kwargs):
        good = _good_evaluator(inner, x, y, z, t, **kwargs)
        total = np.array(good.inner_plus_oscillatory_self_advection, copy=True)
        total[..., 2] += 1.0e-4
        good.inner_plus_oscillatory_self_advection = total
        return good

    with pytest.raises(ValueError, match="full-advection decomposition closure"):
        _witness(full_drift)


def test_derivative_step_and_ring_domain_drift_fail_closed():
    def step_drift(inner, x, y, z, t, **kwargs):
        good = _good_evaluator(inner, x, y, z, t, **kwargs)
        good.inner_leading_spatial_step = 2.0e-3
        return good

    with pytest.raises(ValueError, match="FD4 step drifted"):
        _witness(step_drift)

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
