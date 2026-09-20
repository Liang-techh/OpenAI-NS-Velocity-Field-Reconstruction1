from __future__ import annotations

import numpy as np

from openai_ns_reconstruction.kokuno_a4_strict_inner_transport_artifact_independent_audit import (
    FINAL_DIVERGENCE_GATE,
    FINAL_MOMENTUM_NORMALIZED_GATE,
    HELDOUT_SEED,
    MEDIUM_TO_FINE_RELATIVE_CHANGE_GATE,
    NU,
    SPACE_STEPS,
    TIME_STEPS,
    gate_failures,
    public_contract,
    three_resolution_levels,
    transport_level,
    vector_rms,
)


def _quadratic_divergence_free_velocity(x, y, z, t):
    x, y, z, _ = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack((y * y, z * z, x * x), axis=-1)


def _quadratic_exact_transport(x, y, z, t, nu=NU):
    x, y, z, _ = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )
    return np.stack(
        (
            2.0 * y * z * z - 2.0 * nu,
            2.0 * z * x * x - 2.0 * nu,
            2.0 * x * y * y - 2.0 * nu,
        ),
        axis=-1,
    )


def _divergence_mutation_velocity(x, y, z, t):
    base = _quadratic_divergence_free_velocity(x, y, z, t)
    x = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
        np.asarray(t, dtype=float),
    )[0]
    out = np.array(base, copy=True)
    out[..., 0] += 1.0e-3 * x
    return out


def test_transport_operator_recovers_quadratic_manufactured_field():
    x = np.asarray([-0.31, -0.09, 0.17, 0.42])
    y = np.asarray([0.23, -0.27, 0.35, -0.12])
    z = np.asarray([-0.25, 0.16, 0.32, -0.39])
    t = np.asarray([0.46, 0.49, 0.52, 0.54])
    level = transport_level(
        _quadratic_divergence_free_velocity,
        x,
        y,
        z,
        t,
        space_step=1.5e-3,
        time_step=1.2e-3,
    )
    expected = _quadratic_exact_transport(x, y, z, t)
    np.testing.assert_allclose(level.transport, expected, rtol=0.0, atol=3.0e-10)
    np.testing.assert_allclose(level.divergence, 0.0, rtol=0.0, atol=2.0e-12)
    np.testing.assert_allclose(level.laplacian, 2.0, rtol=0.0, atol=3.0e-9)


def test_three_resolution_gate_is_stable_and_detects_divergence_mutation():
    rng = np.random.default_rng(HELDOUT_SEED)
    x = rng.uniform(-0.4, 0.4, 64)
    y = rng.uniform(-0.4, 0.4, 64)
    z = rng.uniform(-0.4, 0.4, 64)
    t = rng.uniform(0.46, 0.54, 64)
    levels = three_resolution_levels(_quadratic_divergence_free_velocity, x, y, z, t)
    assert gate_failures(levels) == []
    assert vector_rms(levels[-1].velocity) > 1.0e-8
    assert levels[-1].normalized_rms >= 0.0
    assert levels[-1].normalized_sampled_max >= 0.0

    mutated = three_resolution_levels(_divergence_mutation_velocity, x, y, z, t)
    failures = gate_failures(mutated)
    assert "fine_divergence_sampled_max" in failures
    assert "fine_divergence_sampled_rms" in failures


def test_public_contract_freezes_operator_and_truth_boundary():
    contract = public_contract()
    assert contract["heldout_seed"] == HELDOUT_SEED
    assert contract["space_steps"] == list(SPACE_STEPS)
    assert contract["time_steps"] == list(TIME_STEPS)
    assert contract["nu"] == NU == 0.01
    independent = contract["independent_reference"]
    assert independent["artifact_must_be_saved_and_reloaded"] is True
    assert independent["numerical_source"] == "reloaded public velocity(x,y,z,t) only"
    assert independent["agent1_analytic_derivatives_used"] is False
    assert independent["agent2_production_derivatives_used"] is False
    assert independent["pressure_fit_used"] is False
    assert independent["forcing_fit_used"] is False
    gates = contract["frozen_gates"]
    assert gates["medium_to_fine_relative_change"] == MEDIUM_TO_FINE_RELATIVE_CHANGE_GATE
    assert gates["final_momentum_normalized"] == FINAL_MOMENTUM_NORMALIZED_GATE == 1.0e-3
    assert gates["final_divergence"] == FINAL_DIVERGENCE_GATE == 1.0e-5
    truth = contract["truth_boundary"]
    assert truth["strict_inner_transport_precursor_assessed"] is True
    assert truth["transport_precursor_is_complete_ns_residual"] is False
    assert truth["matched_pressure_included"] is False
    assert truth["restricted_forcing_included"] is False
    assert truth["agent3_correction_velocity_included"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False
