from __future__ import annotations

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_a4_oscillatory_transport_delta_independent_audit import (
    FD8_STEPS,
    before_after_fd8_transport,
    frozen_gate_failures,
    public_contract,
    relative_rms,
)

A_INNER = np.asarray(((1.0, 0.2, 0.0), (0.1, -1.0, 0.0), (0.0, 0.0, 0.0)))
A_OSC = np.asarray(((0.0, 0.15, 0.0), (-0.15, 0.0, 0.05), (0.0, -0.05, 0.0)))
B_INNER = np.asarray((0.17, 0.13, -0.24))
B_OSC = np.asarray((0.03, -0.02, 0.04))


def _affine_velocity(matrix, temporal):
    def provider(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float), np.asarray(y, dtype=float),
            np.asarray(z, dtype=float), np.asarray(t, dtype=float),
        )
        xyz = np.stack((x, y, z), axis=-1)
        return np.einsum("ij,...j->...i", matrix, xyz) + t[..., None] * temporal
    return provider


INNER = _affine_velocity(A_INNER, B_INNER)
TOTAL = _affine_velocity(A_INNER + A_OSC, B_INNER + B_OSC)


def _analytic_transport(matrix, temporal, velocity):
    return temporal + np.einsum("ij,...j->...i", matrix, velocity)


def _fixture():
    x = np.asarray((0.41, 0.62, 0.83))
    y = np.asarray((0.12, -0.21, 0.31))
    z = np.asarray((-0.27, 0.18, 0.09))
    t = np.asarray((0.45, 0.50, 0.55))
    levels = tuple(
        before_after_fd8_transport(INNER, TOTAL, x, y, z, t, step=h)
        for h in FD8_STEPS
    )
    inner_velocity = INNER(x, y, z, t)
    total_velocity = TOTAL(x, y, z, t)
    exact_inner = _analytic_transport(A_INNER, B_INNER, inner_velocity)
    exact_total = _analytic_transport(
        A_INNER + A_OSC, B_INNER + B_OSC, total_velocity
    )
    return levels, inner_velocity, total_velocity, exact_inner, exact_total - exact_inner


def test_fd8_before_after_recovers_affine_transport_and_divergence() -> None:
    levels, inner_velocity, total_velocity, exact_inner, exact_delta = _fixture()
    fine = levels[-1]
    assert relative_rms(exact_inner, fine.inner.transport) < 1.0e-8
    assert relative_rms(exact_delta, fine.delta_transport) < 1.0e-7
    np.testing.assert_allclose(fine.inner.divergence, 0.0, atol=1.0e-9, rtol=0.0)
    np.testing.assert_allclose(fine.total.divergence, 0.0, atol=1.0e-9, rtol=0.0)
    assert frozen_gate_failures(
        exact_inner, exact_delta, levels,
        inner_velocity=inner_velocity, total_velocity=total_velocity,
    ) == []


def test_frozen_gate_detects_delta_mutations() -> None:
    levels, inner_velocity, total_velocity, exact_inner, exact_delta = _fixture()
    scaled = frozen_gate_failures(
        exact_inner, 0.90 * exact_delta, levels,
        inner_velocity=inner_velocity, total_velocity=total_velocity,
    )
    reversed_delta = frozen_gate_failures(
        exact_inner, -exact_delta, levels,
        inner_velocity=inner_velocity, total_velocity=total_velocity,
    )
    assert "delta_fine_relative_rms" in scaled
    assert "delta_fine_relative_rms" in reversed_delta


def test_invalid_velocity_shape_fails_closed() -> None:
    def bad_velocity(x, y, z, t):
        x, y, z, t = np.broadcast_arrays(x, y, z, t)
        return np.zeros(x.shape + (2,))
    with pytest.raises(ValueError, match="velocity provider must return shape"):
        before_after_fd8_transport(
            bad_velocity, TOTAL, 0.4, 0.2, -0.1, 0.5, step=FD8_STEPS[-1]
        )


def test_public_contract_keeps_delta_below_full_ns_boundary() -> None:
    contract = public_contract()
    assert contract["schema"] == "kokuno-a4-oscillatory-transport-delta-independent-audit-v1"
    assert contract["agent2_delta_pr"] == 849
    assert contract["viscosity"] == 0.01
    assert contract["independent_reference"]["uses_public_inner_velocity_only"] is True
    assert contract["independent_reference"]["uses_public_total_velocity_only"] is True
    assert contract["independent_reference"]["uses_agent2_fd4_verifier"] is False
    assert contract["norm_scope"]["oscillation_rms_ratio_is_observation_only"] is True
    assert contract["norm_scope"]["same_protocol_st006_comparison_valid"] is False
    truth = contract["truth_boundary"]
    assert truth["strict_inner_before_after_transport_independently_assessed"] is True
    assert truth["oscillatory_transport_delta_independently_assessed"] is True
    for key in (
        "global_leading_velocity_materialized", "outer_join_materialized",
        "matched_pressure_included", "restricted_forcing_included",
        "correction_velocity_included", "complete_kokuno_candidate_assembled",
        "leading_only_ns_residual_assessed",
        "leading_plus_oscillatory_ns_residual_assessed",
        "after_correction_ns_residual_assessed",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed", "boundary_support_gate_assessed",
        "pde_validated", "paper_exact",
    ):
        assert truth[key] is False
