from __future__ import annotations

from dataclasses import replace

import numpy as np

from openai_ns_reconstruction.kokuno_a4_inner_leading_oscillatory_transport_independent_audit import (
    FD8_STEPS,
    VISCOSITY,
    cartesian_fd8_transport,
    frozen_gate_failures,
    public_contract,
    relative_rms,
    transport_at_viscosity,
)


def manufactured_velocity(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack(
        (
            y**3 + 0.1 * z**2 + 0.2 * t * x,
            x**3 + 0.05 * z**3 - 0.2 * t * y,
            0.4 * x * y + 0.3 * t**2,
        ),
        axis=-1,
    )


def manufactured_exact(x, y, z, t):
    x, y, z, t = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    u = manufactured_velocity(x, y, z, t)
    dt = np.stack((0.2 * x, -0.2 * y, 0.6 * t), axis=-1)
    jac = np.empty(x.shape + (3, 3), dtype=float)
    jac[..., 0, 0] = 0.2 * t
    jac[..., 0, 1] = 3.0 * y**2
    jac[..., 0, 2] = 0.2 * z
    jac[..., 1, 0] = 3.0 * x**2
    jac[..., 1, 1] = -0.2 * t
    jac[..., 1, 2] = 0.15 * z**2
    jac[..., 2, 0] = 0.4 * y
    jac[..., 2, 1] = 0.4 * x
    jac[..., 2, 2] = 0.0
    lap = np.stack((6.0 * y + 0.2, 6.0 * x + 0.3 * z, np.zeros_like(x)), axis=-1)
    adv = np.einsum("...j,...ij->...i", u, jac)
    transport = dt + adv - VISCOSITY * lap
    return dt, jac, lap, adv, transport


def test_fd8_operator_recovers_manufactured_polynomial() -> None:
    x = np.asarray((0.17, -0.24, 0.31, -0.11), dtype=float)
    y = np.asarray((-0.22, 0.19, 0.13, 0.28), dtype=float)
    z = np.asarray((0.14, -0.27, 0.21, -0.18), dtype=float)
    t = np.asarray((0.43, 0.47, 0.52, 0.56), dtype=float)
    expected_dt, expected_jac, expected_lap, expected_adv, expected_transport = (
        manufactured_exact(x, y, z, t)
    )
    level = cartesian_fd8_transport(
        manufactured_velocity, x, y, z, t, step=2.0e-3
    )
    np.testing.assert_allclose(level.velocity_dt, expected_dt, rtol=0.0, atol=2.0e-11)
    np.testing.assert_allclose(level.jacobian, expected_jac, rtol=0.0, atol=2.0e-11)
    np.testing.assert_allclose(level.laplacian, expected_lap, rtol=0.0, atol=2.0e-8)
    np.testing.assert_allclose(level.self_advection, expected_adv, rtol=0.0, atol=3.0e-11)
    np.testing.assert_allclose(level.transport, expected_transport, rtol=0.0, atol=3.0e-10)
    assert float(np.max(np.abs(level.divergence))) <= 2.0e-11


def test_frozen_gates_accept_exact_polynomial_and_detect_mutations() -> None:
    x = np.asarray((0.16, 0.21, 0.27, 0.32, 0.24, 0.19), dtype=float)
    y = np.asarray((0.18, -0.23, 0.11, -0.15, 0.29, -0.12), dtype=float)
    z = np.asarray((-0.17, 0.22, -0.26, 0.13, 0.09, 0.28), dtype=float)
    t = np.asarray((0.45, 0.47, 0.49, 0.51, 0.53, 0.55), dtype=float)
    levels = tuple(
        cartesian_fd8_transport(manufactured_velocity, x, y, z, t, step=h)
        for h in FD8_STEPS
    )
    production = manufactured_exact(x, y, z, t)[-1]
    total_velocity = manufactured_velocity(x, y, z, t)
    assert frozen_gate_failures(production, levels, total_velocity=total_velocity) == []

    scale_failures = frozen_gate_failures(
        0.90 * production, levels, total_velocity=total_velocity
    )
    assert "fine_relative_rms" in scale_failures

    flipped = production.copy()
    flipped[..., 2] *= -1.0
    assert frozen_gate_failures(flipped, levels, total_velocity=total_velocity)

    mutated_fine = replace(levels[-1], divergence=levels[-1].divergence + 1.0e-3)
    divergence_failures = frozen_gate_failures(
        production,
        (levels[0], levels[1], mutated_fine),
        total_velocity=total_velocity,
    )
    assert "divergence_sampled_max" in divergence_failures
    assert "divergence_sampled_rms" in divergence_failures


def test_viscosity_perturbation_is_diagnostic_not_retuning() -> None:
    level = cartesian_fd8_transport(
        manufactured_velocity, 0.23, -0.17, 0.19, 0.51, step=FD8_STEPS[-1]
    )
    nominal = transport_at_viscosity(level, VISCOSITY)
    plus = transport_at_viscosity(level, VISCOSITY * 1.001)
    minus = transport_at_viscosity(level, VISCOSITY * 0.999)
    assert relative_rms(nominal, plus) > 0.0
    assert relative_rms(nominal, minus) > 0.0
    np.testing.assert_allclose(level.transport, nominal, rtol=0.0, atol=0.0)


def test_public_contract_keeps_full_pde_gate_fail_closed() -> None:
    contract = public_contract()
    assert contract["schema"] == "kokuno-a4-inner-leading-oscillatory-transport-independent-audit-v1"
    assert contract["agent2_transport_pr"] == 840
    assert contract["viscosity"] == 0.01
    assert contract["fd8_steps"] == list(FD8_STEPS)
    independent = contract["independent_reference"]
    assert independent["uses_public_total_velocity_only"] is True
    assert independent["centered_fd8_time"] is True
    assert independent["cartesian_fd8_jacobian"] is True
    assert independent["cartesian_fd8_laplacian"] is True
    for key in (
        "uses_agent1_velocity_dt",
        "uses_agent1_self_advection",
        "uses_agent1_velocity_laplacian",
        "uses_agent2_production_subterms",
        "uses_agent2_fd4_directional_verifier",
    ):
        assert independent[key] is False

    scope = contract["norm_scope"]
    assert scope["divergence_l2_is_heldout_sampled_rms"] is True
    assert scope["whole_domain_volume_l2_assessed"] is False
    assert scope["complete_ns_residual_assessed"] is False

    truth = contract["truth_boundary"]
    assert truth["strict_inner_leading_plus_oscillatory_transport_assessed"] is True
    assert truth["strict_inner_composite_divergence_assessed"] is True
    for key in (
        "global_leading_velocity_materialized",
        "outer_join_materialized",
        "matched_pressure_included",
        "restricted_forcing_included",
        "correction_velocity_included",
        "complete_kokuno_candidate_assembled",
        "leading_only_ns_residual_assessed",
        "leading_plus_oscillatory_ns_residual_assessed",
        "after_correction_ns_residual_assessed",
        "heldout_ns_momentum_residual_assessed",
        "formal_full_domain_pde_gate_assessed",
        "boundary_support_gate_assessed",
        "pde_validated",
        "paper_exact",
    ):
        assert truth[key] is False
