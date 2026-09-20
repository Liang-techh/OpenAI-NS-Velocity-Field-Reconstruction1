from __future__ import annotations

import inspect

import numpy as np

import openai_ns_reconstruction.kokuno_a4_pa10_cartesian_center_spatial_derivatives_independent_audit as audit


class _ManufacturedVelocityProvider:
    def velocity(self, x, y, z, t):
        x, y, z, t = np.broadcast_arrays(
            np.asarray(x, dtype=float),
            np.asarray(y, dtype=float),
            np.asarray(z, dtype=float),
            np.asarray(t, dtype=float),
        )
        del t
        return np.stack(
            (
                x**5 + 2.0 * y**3 - z**2,
                -3.0 * x**2 + y**4 + 0.5 * z**3,
                x * y + 2.0 * y * z + z**5,
            ),
            axis=-1,
        )


def _manufactured_jacobian(x, y, z):
    x, y, z = np.broadcast_arrays(
        np.asarray(x, dtype=float),
        np.asarray(y, dtype=float),
        np.asarray(z, dtype=float),
    )
    return np.stack(
        (
            np.stack((5.0 * x**4, 6.0 * y**2, -2.0 * z), axis=-1),
            np.stack((-6.0 * x, 4.0 * y**3, 1.5 * z**2), axis=-1),
            np.stack((y, x + 2.0 * z, 2.0 * y + 5.0 * z**4), axis=-1),
        ),
        axis=-2,
    )


def test_public_velocity_only_fd6_spatial_orientation_and_curl() -> None:
    provider = _ManufacturedVelocityProvider()
    x = np.asarray([0.13, -0.27, 0.41])
    y = np.asarray([-0.31, 0.19, 0.07])
    z = np.asarray([0.23, -0.11, 0.37])
    t = np.asarray([0.37, 0.51, 0.63])

    actual = audit._fd6_jacobian(provider, x, y, z, t, 1.0e-3)
    expected = _manufactured_jacobian(x, y, z)
    assert np.max(np.abs(actual - expected)) < 2.0e-10

    actual_curl = audit._curl_from_jacobian(actual)
    expected_curl = audit._curl_from_jacobian(expected)
    assert np.max(np.abs(actual_curl - expected_curl)) < 3.0e-10


def test_independent_reference_avoids_agent1_spatial_derivative_oracles() -> None:
    source = inspect.getsource(audit._fd6_axis)
    assert "provider.velocity(" in source
    forbidden_calls = (
        ".velocity_jacobian(",
        ".divergence(",
        ".vorticity(",
        ".coordinate_spatial_derivatives(",
        ".physical_profiles.derivatives(",
    )
    for forbidden in forbidden_calls:
        assert forbidden not in source

    assert audit.SPACE_STEPS == (1.2e-3, 6.0e-4, 3.0e-4)
    assert audit.SEED == 9173481


def test_pa10_cartesian_center_spatial_derivative_independent_audit_contract() -> None:
    result = audit.run_independent_audit()
    payload = result.payload

    assert payload["schema"] == audit.SCHEMA
    assert payload["seed"] == 9173481
    assert payload["sample_count_random"] == 384
    assert payload["sample_count_total"] == 389
    assert payload["space_steps"] == [1.2e-3, 6.0e-4, 3.0e-4]

    frozen = payload["frozen_gates"]
    assert frozen["independent_divergence_sampled_max"] == 1.0e-5
    assert frozen["independent_divergence_sampled_l2_rms"] == 1.0e-5
    assert frozen["final_normalized_momentum_max_l2"] == 1.0e-3
    assert frozen["final_divergence_max_l2"] == 1.0e-5

    assert set(payload["measurements"]["by_step"]) == {
        "0.0012000",
        "0.0006000",
        "0.0003000",
    }
    for step_payload in payload["measurements"]["by_step"].values():
        assert len(step_payload["component_errors"]) == 9

    truth = payload["truth_boundary"]
    assert truth["agent4_public_velocity_only_fd6_spatial_derivative_used"] is True
    assert truth["agent4_three_spatial_resolutions_used"] is True
    assert truth["source_complex_C_normalization_independently_admitted"] is False
    assert truth["global_cartesian_spacetime_leading_velocity_materialized"] is False
    assert truth["matched_global_pressure_materialized"] is False
    assert truth["complete_restricted_forcing_materialized"] is False
    assert truth["leading_only_ns_residual_assessed"] is False
    assert truth["leading_plus_oscillatory_ns_residual_assessed"] is False
    assert truth["after_correction_ns_residual_assessed"] is False
    assert truth["pde_validated"] is False

    assert payload["gates"]["mutations_detected"] is True
    assert truth["inner_cartesian_center_spatial_derivatives_independently_admitted"] is payload["passed"]
    assert result.passed is payload["passed"]
    assert len(payload["receipt_sha256"]) == 64
