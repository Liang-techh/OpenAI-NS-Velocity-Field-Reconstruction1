from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_transport_delta import (
    INDEPENDENT_FD4_STEPS,
    VISCOSITY,
    evaluate_inner_leading_oscillatory_transport_delta,
    independent_fd4_inner_transport,
    independent_fd4_oscillatory_transport_delta,
    public_contract,
)


class PolynomialInnerBackend:
    """Manufactured strict-inner backend; mechanics fixture only."""

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack(
            (
                x + 0.17 * t + 0.04 * y * z,
                -0.8 * y + 0.13 * t + 0.03 * x * z,
                0.21 * z - 0.24 * t + 0.025 * x * y,
            ),
            axis=-1,
        )

    def velocity_dt(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        del y, z, t
        return np.broadcast_to(
            np.asarray((0.17, 0.13, -0.24)), x.shape + (3,)
        ).copy()

    def self_advection(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        u = self.velocity(x, y, z, t)
        jac = np.empty(x.shape + (3, 3), dtype=float)
        jac[..., 0, 0] = 1.0
        jac[..., 0, 1] = 0.04 * z
        jac[..., 0, 2] = 0.04 * y
        jac[..., 1, 0] = 0.03 * z
        jac[..., 1, 1] = -0.8
        jac[..., 1, 2] = 0.03 * x
        jac[..., 2, 0] = 0.025 * y
        jac[..., 2, 1] = 0.025 * x
        jac[..., 2, 2] = 0.21
        return np.einsum("...ij,...j->...i", jac, u)

    def velocity_laplacian(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        del y, z, t
        return np.zeros(x.shape + (3,), dtype=float)


def test_production_transport_delta_closes_before_after_decomposition() -> None:
    backend = PolynomialInnerBackend()
    x = np.asarray((0.54, 0.76, 0.91), dtype=float)
    y = np.asarray((0.19, -0.23, 0.29), dtype=float)
    z = np.asarray((-0.38, 0.14, 0.49), dtype=float)
    t = np.asarray((0.46, 0.50, 0.54), dtype=float)

    out = evaluate_inner_leading_oscillatory_transport_delta(
        backend, x, y, z, t
    )
    assert out.viscosity == VISCOSITY == 0.01
    assert out.oscillatory_transport_increment.shape == (3, 3)
    assert np.all(np.isfinite(out.oscillatory_transport_increment))
    np.testing.assert_allclose(
        out.inner_plus_oscillatory_transport,
        out.inner_leading_transport + out.oscillatory_transport_increment,
        rtol=0.0,
        atol=5.0e-13,
    )
    np.testing.assert_allclose(
        out.oscillatory_transport_increment,
        out.oscillatory_time_increment
        + out.oscillatory_nonlinear_increment
        + out.oscillatory_viscous_increment,
        rtol=0.0,
        atol=5.0e-13,
    )
    np.testing.assert_array_equal(
        out.inner_plus_oscillatory_velocity,
        out.inner_leading_velocity + out.oscillatory_velocity,
    )


def test_independent_inner_and_delta_fd4_are_finite_and_stabilize() -> None:
    backend = PolynomialInnerBackend()
    x = np.asarray((0.57, 0.81), dtype=float)
    y = np.asarray((0.20, -0.25), dtype=float)
    z = np.asarray((-0.33, 0.45), dtype=float)
    t = np.asarray((0.47, 0.53), dtype=float)

    inner_values = [
        independent_fd4_inner_transport(backend, x, y, z, t, step=h)
        for h in INDEPENDENT_FD4_STEPS
    ]
    delta_values = [
        independent_fd4_oscillatory_transport_delta(
            backend, x, y, z, t, step=h
        )
        for h in INDEPENDENT_FD4_STEPS
    ]
    assert all(value.shape == (2, 3) for value in inner_values + delta_values)
    assert all(np.all(np.isfinite(value)) for value in inner_values + delta_values)

    inner_coarse = np.linalg.norm(inner_values[0] - inner_values[1])
    inner_fine = np.linalg.norm(inner_values[1] - inner_values[2])
    delta_coarse = np.linalg.norm(delta_values[0] - delta_values[1])
    delta_fine = np.linalg.norm(delta_values[1] - delta_values[2])
    assert inner_fine < inner_coarse
    assert delta_fine < delta_coarse


def test_missing_required_backend_method_fails_closed() -> None:
    class IncompleteBackend:
        def velocity(self, x, y, z, t):
            return PolynomialInnerBackend().velocity(x, y, z, t)

        def velocity_dt(self, x, y, z, t):
            return PolynomialInnerBackend().velocity_dt(x, y, z, t)

        def velocity_laplacian(self, x, y, z, t):
            return PolynomialInnerBackend().velocity_laplacian(x, y, z, t)

    with pytest.raises(TypeError, match="self_advection"):
        evaluate_inner_leading_oscillatory_transport_delta(
            IncompleteBackend(), 0.72, 0.18, 0.12, 0.5
        )


def test_public_contract_keeps_transport_delta_below_residual_boundary() -> None:
    contract = public_contract()
    assert contract["schema"] == "kokuno-a2-inner-leading-oscillatory-transport-delta-v1"
    assert contract["parent_agent2_pr"] == 840
    assert contract["viscosity"] == 0.01
    assert contract["forbidden_public_inputs_present"] == []

    signature = inspect.signature(evaluate_inner_leading_oscillatory_transport_delta)
    for forbidden in (
        "residual", "defect", "pressure", "forcing", "target", "gain",
        "delta_a", "viscosity", "nu", "scientific_threshold", "spatial_step",
    ):
        assert forbidden not in signature.parameters

    truth = contract["truth_boundary"]
    assert truth["strict_inner_before_after_transport_executable"] is True
    assert truth["oscillatory_transport_increment_executable"] is True
    assert truth["viscosity_fixed_not_tunable"] is True
    for key in (
        "diagnostic_is_complete_ns_residual",
        "matched_pressure_included",
        "restricted_forcing_included",
        "correction_velocity_included",
        "inner_leading_is_final_corrected_fixed_point",
        "global_leading_velocity_materialized",
        "outer_join_materialized",
        "agent3_mean_radial_chain_reimplemented",
        "complete_kokuno_candidate_assembled",
        "heldout_ns_residual_assessed",
        "st006_same_protocol_comparison_valid",
        "residual_reduction_claimed",
        "paper_exact",
        "pde_validated",
    ):
        assert truth[key] is False
