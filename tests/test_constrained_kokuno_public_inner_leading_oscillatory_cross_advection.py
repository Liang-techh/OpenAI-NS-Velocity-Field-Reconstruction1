from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_cross_advection import (
    directional_fd4_cross_advection,
    evaluate_inner_leading_oscillatory_cross_advection,
    public_contract,
)


class PolynomialInnerLeading:
    """Manufactured backend; not Kokuno scientific evidence."""

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack(
            (
                x + 2.0 * y + 0.5 * z + 0.2 * t,
                0.3 * x * x - y + z - 0.1 * t,
                z + x * y + 0.15 * t,
            ),
            axis=-1,
        )

    @staticmethod
    def jacobian(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        out = np.empty(x.shape + (3, 3), dtype=float)
        out[..., 0, :] = np.stack(
            (np.ones_like(x), 2.0 * np.ones_like(x), 0.5 * np.ones_like(x)), axis=-1
        )
        out[..., 1, :] = np.stack(
            (0.6 * x, -np.ones_like(x), np.ones_like(x)), axis=-1
        )
        out[..., 2, :] = np.stack((y, x, np.ones_like(x)), axis=-1)
        return out


class BadShapeBackend:
    def velocity(self, x, y, z, t):
        x = np.asarray(x, dtype=float)
        return np.zeros(x.shape + (2,), dtype=float)


def _points():
    radius = np.asarray([0.48, 0.57, 0.66, 0.75, 0.84, 0.93], dtype=float)
    angle = np.asarray([0.31, 1.08, 1.83, 2.54, -2.37, -1.12], dtype=float)
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)
    z = np.asarray([-0.82, -0.43, -0.08, 0.29, 0.61, 0.94], dtype=float)
    t = np.asarray([0.46, 0.48, 0.50, 0.52, 0.54, 0.56], dtype=float)
    return x, y, z, t


def _vector_rms(a):
    a = np.asarray(a, dtype=float)
    return float(np.sqrt(np.mean(np.sum(a * a, axis=-1))))


def test_manufactured_inner_jacobian_and_cross_assembly():
    backend = PolynomialInnerLeading()
    x, y, z, t = _points()
    result = evaluate_inner_leading_oscillatory_cross_advection(
        backend,
        x,
        y,
        z,
        t,
        inner_leading_spatial_step=1.0e-3,
        oscillatory_spatial_step=1.0e-3,
    )

    expected_j = backend.jacobian(x, y, z, t)
    assert np.max(np.abs(result.inner_leading_jacobian - expected_j)) < 2.0e-11
    np.testing.assert_array_equal(
        result.inner_plus_oscillatory_velocity,
        result.inner_leading_velocity + result.oscillatory_velocity,
    )
    np.testing.assert_array_equal(
        result.cross_advection,
        result.inner_advects_oscillation + result.oscillation_advects_inner,
    )

    expected_inner_advects_osc = np.einsum(
        "...ij,...j->...i", result.oscillatory_jacobian, result.inner_leading_velocity
    )
    expected_osc_advects_inner = np.einsum(
        "...ij,...j->...i", expected_j, result.oscillatory_velocity
    )
    assert np.max(np.abs(result.inner_advects_oscillation - expected_inner_advects_osc)) == 0.0
    assert np.max(np.abs(result.oscillation_advects_inner - expected_osc_advects_inner)) < 1.0e-11
    assert _vector_rms(result.cross_advection) > 1.0e-10


def test_directional_fd4_is_independent_consistency_check():
    backend = PolynomialInnerLeading()
    x, y, z, t = _points()
    production = evaluate_inner_leading_oscillatory_cross_advection(
        backend, x, y, z, t, inner_leading_spatial_step=1.0e-3, oscillatory_spatial_step=1.0e-3
    ).cross_advection
    direct = directional_fd4_cross_advection(
        backend, x, y, z, t, directional_step=1.0e-3
    )
    relative = _vector_rms(direct - production) / max(_vector_rms(production), np.finfo(float).tiny)
    assert relative < 5.0e-3


def test_contract_truth_boundary_and_forbidden_inputs():
    contract = public_contract()
    assert contract["forbidden_public_inputs_present"] == []
    truth = contract["truth_boundary"]
    assert truth["inner_plus_oscillatory_velocity_executable"] is True
    assert truth["inner_oscillatory_cross_advection_executable"] is True
    assert truth["agent1_exact_head_ci_assumed_passed"] is False
    assert truth["agent4_inner_field_independent_admission_assumed"] is False
    assert truth["inner_leading_is_final_corrected_fixed_point"] is False
    assert truth["global_leading_velocity_materialized"] is False
    assert truth["matched_pressure_included"] is False
    assert truth["restricted_forcing_included"] is False
    assert truth["agent3_mean_radial_chain_reimplemented"] is False
    assert truth["complete_kokuno_candidate_assembled"] is False
    assert truth["heldout_ns_residual_assessed"] is False
    assert truth["paper_exact"] is False
    assert truth["pde_validated"] is False


def test_fail_closed_backend_shape_and_steps():
    x, y, z, t = _points()
    with pytest.raises(ValueError, match="shape"):
        evaluate_inner_leading_oscillatory_cross_advection(BadShapeBackend(), x, y, z, t)
    with pytest.raises(ValueError, match="inner_leading_spatial_step"):
        evaluate_inner_leading_oscillatory_cross_advection(
            PolynomialInnerLeading(), x, y, z, t, inner_leading_spatial_step=0.0
        )
    with pytest.raises(ValueError, match="oscillatory_spatial_step"):
        evaluate_inner_leading_oscillatory_cross_advection(
            PolynomialInnerLeading(), x, y, z, t, oscillatory_spatial_step=np.inf
        )
