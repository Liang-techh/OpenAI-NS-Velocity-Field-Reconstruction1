from __future__ import annotations

from dataclasses import dataclass
import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_current_partial_nonlinear_mean_attribution import (
    COMPOSITION_CLOSURE_GATE,
    FIXED_LEADING_SPATIAL_STEP,
    ExactCurrentPartialNonlinearBackend,
    _components_from_interfaces,
    materialize_current_partial_nonlinear_mean_attribution,
    truth_boundary,
)


@dataclass(frozen=True)
class _OscEval:
    velocity: np.ndarray
    velocity_jacobian: np.ndarray


def _leading(x, y, z, t):
    del t
    x, y, z = np.broadcast_arrays(x, y, z)
    return np.stack((x + 2.0 * y, -y + 0.5 * z, 3.0 * z - x), axis=-1)


def _osc_diff(points, time):
    del time
    p = np.asarray(points, dtype=float)
    x, y, z = p[..., 0], p[..., 1], p[..., 2]
    velocity = np.stack((y, z, x), axis=-1)
    jac = np.zeros(p.shape[:-1] + (3, 3), dtype=float)
    jac[..., 0, 1] = 1.0
    jac[..., 1, 2] = 1.0
    jac[..., 2, 0] = 1.0
    return _OscEval(velocity=velocity, velocity_jacobian=jac)


def _composite(x, y, z, t):
    lead = _leading(x, y, z, t)
    p = np.stack(np.broadcast_arrays(x, y, z), axis=-1)
    return lead + _osc_diff(p, t).velocity


def test_internal_formula_adapter_reconstructs_three_nonlinear_pieces():
    x = np.asarray([0.31, -0.42, 0.53])
    y = np.asarray([-0.17, 0.28, 0.19])
    z = np.asarray([0.11, -0.23, 0.07])
    t = np.asarray([0.41, 0.47, 0.53])

    a, b, q, mixed, aggregate = _components_from_interfaces(
        _leading, _composite, _osc_diff, x, y, z, t
    )

    leading = _leading(x, y, z, t)
    osc = np.stack((y, z, x), axis=-1)
    j_osc = np.asarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    j_lead = np.asarray([[1.0, 2.0, 0.0], [0.0, -1.0, 0.5], [-1.0, 0.0, 3.0]])
    expected_a = np.einsum("ij,...j->...i", j_osc, leading)
    expected_b = np.einsum("ij,...j->...i", j_lead, osc)
    expected_q = np.einsum("ij,...j->...i", j_osc, osc)

    assert np.max(np.abs(a - expected_a)) < 2.0e-12
    assert np.max(np.abs(b - expected_b)) < 2.0e-12
    assert np.max(np.abs(q - expected_q)) < 2.0e-12
    assert np.max(np.abs(mixed - (expected_a + expected_b))) < 2.0e-12
    assert np.max(np.abs(aggregate - (expected_a + expected_b + expected_q))) < 2.0e-12


def test_internal_formula_adapter_fails_on_composition_drift():
    def bad_composite(x, y, z, t):
        out = _composite(x, y, z, t).copy()
        out[..., 0] += 10.0 * COMPOSITION_CLOSURE_GATE
        return out

    with pytest.raises(RuntimeError, match="no longer closes"):
        _components_from_interfaces(
            _leading, bad_composite, _osc_diff, 0.4, -0.2, 0.1, 0.5
        )


def test_exact_binder_rejects_unpinned_backend():
    class FakeField:
        leading_backend = object()

        def velocity(self, x, y, z, t):
            return _leading(x, y, z, t)

        def semantic_payload(self):
            return {}

    with pytest.raises(ValueError, match="module identity"):
        ExactCurrentPartialNonlinearBackend.bind(FakeField(), _osc_diff)


def test_public_contract_has_no_scientific_tuning_knobs():
    params = inspect.signature(materialize_current_partial_nonlinear_mean_attribution).parameters
    assert list(params) == ["backend", "radius", "z", "t"]
    for forbidden in (
        "residual", "defect", "mean", "stress", "pressure", "forcing", "gain",
        "damping", "spatial_step", "derivative_step", "angular_order", "viscosity",
        "scientific_threshold", "correction",
    ):
        assert forbidden not in params
    assert FIXED_LEADING_SPATIAL_STEP == 1.0e-3


def test_truth_boundary_stays_scoped_and_fail_closed():
    b = truth_boundary()
    assert b["current_partial_leading_plus_oscillation_consumed"] is True
    assert b["current_partial_mixed_nonlinear_mean_materialized"] is True
    assert b["current_partial_quadratic_nonlinear_mean_materialized"] is True
    assert b["agent2_oscillatory_jacobian_reimplemented_by_agent3"] is False
    assert b["caller_tunable_leading_spatial_step"] is False
    assert b["forbidden_public_parameters_absent"] is True
    for key in (
        "radial_inverse_performed_in_this_increment",
        "pressure_gradient_included",
        "restricted_forcing_included",
        "velocity_beyond_Xh_materialized",
        "outer_global_leading_velocity_materialized",
        "complete_ns_defect",
        "scoped_nonlinear_mean_authorized_as_correction_target",
        "mean_correction_velocity_materialized",
        "real_candidate_finite_correction_cycle_run",
        "heldout_normalized_ns_residual_assessed",
        "residual_reduction_claimed",
        "same_protocol_comparable_to_st006",
        "paper_exact",
        "pde_validated",
        "blowup_proved",
    ):
        assert b[key] is False, key
    assert b["final_normalized_momentum_gate"] == 1.0e-3
    assert b["final_normalized_divergence_gate"] == 1.0e-5
