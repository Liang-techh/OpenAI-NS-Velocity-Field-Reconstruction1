from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_a4_inner_leading_oscillatory_full_advection_independent_audit import (
    FD6_STEPS,
    cartesian_fd6_self_advection,
    frozen_gate_failures,
    public_contract,
)


def _polynomial_velocity(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack(
        (
            1.0 + x + 2.0 * y + 3.0 * z + 0.25 * t,
            2.0 - x * x + y * z - 0.5 * t,
            -1.0 + 0.5 * x * y + z**3 + 0.125 * t,
        ),
        axis=-1,
    )


def _analytic_self_advection(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    u = _polynomial_velocity(x, y, z, t)
    jac = np.empty(x.shape + (3, 3), dtype=float)
    jac[..., 0, 0] = 1.0
    jac[..., 0, 1] = 2.0
    jac[..., 0, 2] = 3.0
    jac[..., 1, 0] = -2.0 * x
    jac[..., 1, 1] = z
    jac[..., 1, 2] = y
    jac[..., 2, 0] = 0.5 * y
    jac[..., 2, 1] = 0.5 * x
    jac[..., 2, 2] = 3.0 * z * z
    return np.einsum("...j,...ij->...i", u, jac)


def test_fd6_cartesian_advection_matches_manufactured_polynomial():
    x = np.asarray([-0.31, -0.07, 0.19, 0.42])
    y = np.asarray([0.22, -0.26, 0.11, -0.18])
    z = np.asarray([-0.24, 0.09, 0.27, -0.05])
    t = np.asarray([0.41, 0.48, 0.53, 0.59])
    exact = _analytic_self_advection(x, y, z, t)
    level = cartesian_fd6_self_advection(
        _polynomial_velocity, x, y, z, t, step=4.0e-4
    )
    assert np.max(np.abs(level.self_advection - exact)) < 2.0e-10


def test_frozen_gate_detects_scaled_and_component_sign_mutations():
    x = np.asarray([-0.31, -0.07, 0.19, 0.42])
    y = np.asarray([0.22, -0.26, 0.11, -0.18])
    z = np.asarray([-0.24, 0.09, 0.27, -0.05])
    t = np.asarray([0.41, 0.48, 0.53, 0.59])
    exact = _analytic_self_advection(x, y, z, t)
    levels = tuple(
        cartesian_fd6_self_advection(_polynomial_velocity, x, y, z, t, step=h)
        for h in FD6_STEPS
    )
    velocity = _polynomial_velocity(x, y, z, t)
    assert frozen_gate_failures(exact, levels, total_velocity=velocity) == []

    scaled = 0.99 * exact
    assert "fine_relative_rms" in frozen_gate_failures(
        scaled, levels, total_velocity=velocity
    )

    flipped = exact.copy()
    flipped[..., 2] *= -1.0
    failures = frozen_gate_failures(flipped, levels, total_velocity=velocity)
    assert "fine_relative_rms" in failures
    assert "fine_relative_sampled_max" in failures


def test_truth_boundary_stays_fail_closed():
    contract = public_contract()
    assert contract["fd6_steps"] == list(FD6_STEPS)
    assert contract["independent_reference"] == {
        "uses_public_total_velocity_only": True,
        "cartesian_fd6": True,
        "uses_agent1_self_advection": False,
        "uses_a2_cross_advection": False,
        "uses_a2_oscillatory_self_advection": False,
        "uses_a2_directional_fd4_verifier": False,
    }
    truth = contract["truth_boundary"]
    assert truth["strict_inner_leading_plus_oscillatory_advection_assessed"] is True
    for key in (
        "global_leading_velocity_materialized",
        "outer_join_materialized",
        "matched_pressure_included",
        "restricted_forcing_included",
        "viscous_laplacian_included",
        "correction_velocity_included",
        "complete_kokuno_candidate_assembled",
        "leading_only_ns_residual_assessed",
        "leading_plus_oscillatory_ns_residual_assessed",
        "after_correction_ns_residual_assessed",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
    ):
        assert truth[key] is False
