from __future__ import annotations

import inspect

import numpy as np

import openai_ns_reconstruction.kokuno_a4_pa10_cartesian_center_velocity_dt_independent_audit as audit


class _ManufacturedVelocityProvider:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        return np.stack(
            (
                x + t**4,
                y - 2.0 * t**3,
                z + 3.0 * t**2,
            ),
            axis=-1,
        )


def test_public_velocity_only_fd4_time_derivative_orientation() -> None:
    provider = _ManufacturedVelocityProvider()
    x = np.asarray([0.13, -0.27, 0.41])
    y = np.asarray([-0.31, 0.19, 0.07])
    z = np.asarray([0.23, -0.11, 0.37])
    t = np.asarray([0.37, 0.51, 0.63])

    actual = audit._fd4_velocity_dt(provider, x, y, z, t, 1.0e-3)
    expected = np.stack((4.0 * t**3, -6.0 * t**2, 6.0 * t), axis=-1)
    assert np.max(np.abs(actual - expected)) < 2.0e-11


def test_independent_audit_avoids_agent1_derivative_oracles() -> None:
    source = inspect.getsource(audit)
    forbidden_calls = (
        ".coordinate_time_derivatives(",
        ".v0_eta(",
        ".physical_profiles.derivatives(",
        ".report(",
    )
    for forbidden in forbidden_calls:
        assert forbidden not in source

    assert "provider.velocity(x, y, z, t - 2.0 * step)" in source
    assert "provider.velocity(x, y, z, t + 2.0 * step)" in source
    assert audit.TIME_STEPS == (4.0e-4, 2.0e-4, 1.0e-4)


def test_pa10_cartesian_center_velocity_dt_independent_audit() -> None:
    result = audit.run_independent_audit()
    payload = result.payload

    assert payload["schema"] == audit.SCHEMA
    assert payload["seed"] == 9173471
    assert payload["sample_count_random"] == 1024
    assert payload["sample_count_total"] == 1029
    assert payload["time_steps"] == [4.0e-4, 2.0e-4, 1.0e-4]

    assert payload["frozen_gates"]["final_normalized_momentum_max_l2"] == 1.0e-3
    assert payload["frozen_gates"]["final_divergence_max_l2"] == 1.0e-5

    truth = payload["truth_boundary"]
    assert truth["agent4_public_velocity_only_fd4_time_derivative_used"] is True
    assert truth["agent4_three_time_resolutions_used"] is True
    assert truth["source_complex_C_normalization_independently_admitted"] is False
    assert truth["global_cartesian_spacetime_leading_velocity_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["complete_restricted_forcing_materialized"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False

    assert set(payload["measurements"]["by_step"]) == {
        "0.0004000",
        "0.0002000",
        "0.0001000",
    }
    for step_payload in payload["measurements"]["by_step"].values():
        assert set(step_payload["component_errors"]) == {"u_t", "v_t", "w_t"}

    assert payload["gates"]["mutations_detected"] is True
    assert payload["passed"] is True
    assert result.passed is True
    assert len(payload["receipt_sha256"]) == 64
