from __future__ import annotations

import inspect

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_public_inner_leading_oscillatory_transport import (
    INDEPENDENT_FD4_STEPS,
    VISCOSITY,
    evaluate_inner_leading_oscillatory_transport,
    independent_fd4_inner_plus_oscillatory_transport,
    public_contract,
)


class PolynomialInnerBackend:
    """Manufactured strict-inner backend; mechanics fixture, not Kokuno evidence."""

    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack(
            (
                x + 0.2 * t + 0.05 * y * z,
                -y + 0.1 * t + 0.04 * x * z,
                0.15 * z - 0.3 * t + 0.03 * x * y,
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
        return np.broadcast_to(np.asarray((0.2, 0.1, -0.3)), x.shape + (3,)).copy()

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
        jac[..., 0, 1] = 0.05 * z
        jac[..., 0, 2] = 0.05 * y
        jac[..., 1, 0] = 0.04 * z
        jac[..., 1, 1] = -1.0
        jac[..., 1, 2] = 0.04 * x
        jac[..., 2, 0] = 0.03 * y
        jac[..., 2, 1] = 0.03 * x
        jac[..., 2, 2] = 0.15
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


def test_transport_assembly_and_fixed_viscosity() -> None:
    backend = PolynomialInnerBackend()
    x = np.asarray((0.51, 0.72, 0.93), dtype=float)
    y = np.asarray((0.18, -0.27, 0.31), dtype=float)
    z = np.asarray((-0.42, 0.16, 0.53), dtype=float)
    t = np.asarray((0.43, 0.50, 0.57), dtype=float)

    out = evaluate_inner_leading_oscillatory_transport(backend, x, y, z, t)
    assert out.transport.shape == (3, 3)
    assert np.all(np.isfinite(out.transport))
    assert out.viscosity == VISCOSITY == 0.01
    np.testing.assert_array_equal(
        out.inner_plus_oscillatory_velocity,
        out.inner_leading_velocity + out.oscillatory_velocity,
    )
    np.testing.assert_array_equal(
        out.inner_plus_oscillatory_velocity_dt,
        out.inner_leading_velocity_dt + out.oscillatory_velocity_dt,
    )
    np.testing.assert_array_equal(
        out.inner_plus_oscillatory_laplacian,
        out.inner_leading_laplacian + out.oscillatory_laplacian,
    )
    np.testing.assert_allclose(
        out.transport,
        out.inner_plus_oscillatory_velocity_dt
        + out.inner_plus_oscillatory_self_advection
        - VISCOSITY * out.inner_plus_oscillatory_laplacian,
        rtol=0.0,
        atol=1.0e-14,
    )


def test_independent_fd4_transport_is_finite_and_stabilizes() -> None:
    backend = PolynomialInnerBackend()
    x = np.asarray((0.58, 0.79), dtype=float)
    y = np.asarray((0.21, -0.24), dtype=float)
    z = np.asarray((-0.34, 0.47), dtype=float)
    t = np.asarray((0.47, 0.53), dtype=float)

    values = [
        independent_fd4_inner_plus_oscillatory_transport(
            backend, x, y, z, t, step=h
        )
        for h in INDEPENDENT_FD4_STEPS
    ]
    assert all(value.shape == (2, 3) for value in values)
    assert all(np.all(np.isfinite(value)) for value in values)
    coarse = np.linalg.norm(values[0] - values[1])
    fine = np.linalg.norm(values[1] - values[2])
    assert fine < coarse


def test_missing_laplacian_backend_fails_closed() -> None:
    class IncompleteBackend:
        def velocity(self, x, y, z, t):
            return PolynomialInnerBackend().velocity(x, y, z, t)

        def velocity_dt(self, x, y, z, t):
            return PolynomialInnerBackend().velocity_dt(x, y, z, t)

        def self_advection(self, x, y, z, t):
            return PolynomialInnerBackend().self_advection(x, y, z, t)

    with pytest.raises(TypeError, match="velocity_laplacian"):
        evaluate_inner_leading_oscillatory_transport(
            IncompleteBackend(), 0.7, 0.2, 0.1, 0.5
        )


def test_public_contract_keeps_scope_fail_closed() -> None:
    contract = public_contract()
    assert contract["schema"] == "kokuno-a2-inner-leading-oscillatory-transport-v1"
    assert contract["viscosity"] == 0.01
    assert contract["forbidden_public_inputs_present"] == []
    signature = inspect.signature(evaluate_inner_leading_oscillatory_transport)
    for forbidden in (
        "residual", "defect", "pressure", "forcing", "target", "gain",
        "delta_a", "viscosity", "nu", "scientific_threshold", "spatial_step",
    ):
        assert forbidden not in signature.parameters

    truth = contract["truth_boundary"]
    assert truth["strict_inner_transport_executable"] is True
    assert truth["viscosity_fixed_not_tunable"] is True
    for key in (
        "agent1_exact_head_ci_assumed_passed",
        "agent4_inner_field_independent_admission_assumed",
        "agent4_a2_oscillatory_independent_admission_assumed",
        "inner_leading_is_final_corrected_fixed_point",
        "global_leading_velocity_materialized",
        "outer_join_materialized",
        "matched_pressure_included",
        "restricted_forcing_included",
        "agent3_mean_radial_chain_reimplemented",
        "correction_velocity_included",
        "complete_kokuno_candidate_assembled",
        "complete_ns_residual",
        "heldout_ns_residual_assessed",
        "residual_reduction_claimed",
        "paper_exact",
        "pde_validated",
    ):
        assert truth[key] is False
