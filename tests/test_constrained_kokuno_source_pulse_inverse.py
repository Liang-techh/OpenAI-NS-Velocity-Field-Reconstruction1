from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_complete_curl import (
    KokunoSourceCompleteCurlContract,
)
from openai_ns_reconstruction.kokuno_source_pulse_inverse import (
    KokunoProjectedPulseInverseContract,
)


def _manufactured_path(contract: KokunoProjectedPulseInverseContract, count: int):
    v = np.linspace(0.0, 1.0, count)
    alpha = 0.7
    n = np.stack(
        (np.cos(alpha * v), np.sin(alpha * v), np.zeros_like(v)), axis=-1
    )
    n_prime = alpha * np.stack(
        (-np.sin(alpha * v), np.cos(alpha * v), np.zeros_like(v)), axis=-1
    )
    tangent = np.stack(
        (-np.sin(alpha * v), np.cos(alpha * v), np.zeros_like(v)), axis=-1
    )
    tangent_prime = -alpha * n
    g = v * (1.0 + 0.2 * v)
    g_prime = 1.0 + 0.4 * v
    ez = np.array([0.0, 0.0, 1.0])
    exact = g[:, None] * tangent + 0.15j * g[:, None] * ez
    exact_prime = (
        g_prime[:, None] * tangent
        + g[:, None] * tangent_prime
        + 0.15j * g_prime[:, None] * ez
    )

    K = np.zeros((count, 3, 3), dtype=np.complex128)
    K[:, 0, 0] = 0.10 + 0.03 * v
    K[:, 0, 1] = -0.07
    K[:, 1, 0] = 0.04
    K[:, 1, 1] = -0.08 + 0.02 * v
    K[:, 2, 0] = 0.03
    K[:, 2, 2] = 0.05

    d = contract.diffusion(n)
    forcing = -(
        exact_prime
        + np.einsum("...ij,...j->...i", K, exact)
        + (contract.m**2 * d)[:, None] * exact
    )
    return v, n, n_prime, K, forcing, exact, exact_prime


def test_projected_pulse_inverse_refines_on_rotating_manufactured_solution():
    contract = KokunoProjectedPulseInverseContract(epsilon=0.2, m=2)
    errors = []
    projection_corrections = []
    for count in (33, 65, 129):
        v, n, n_prime, K, forcing, exact, _ = _manufactured_path(contract, count)
        solved = contract.solve_path(v, n, n_prime, K, forcing)
        t = solved["t_m"]
        errors.append(float(np.max(np.linalg.norm(t - exact, axis=-1))))
        projection_corrections.append(float(solved["max_projection_correction_norm"]))
        assert solved["max_constraint_defect_abs"] < 2.0e-15

    # Tabulated coefficient interpolation is the limiting second-order step.
    assert errors[0] / errors[1] > 3.5
    assert errors[1] / errors[2] > 3.5
    assert errors[2] < 4.0e-6
    assert projection_corrections[2] < projection_corrections[1] < projection_corrections[0]
    assert projection_corrections[2] < 3.0e-10


def test_pressure_recovers_full_mode_equation_and_active_pulse_feeds_source_curl():
    contract = KokunoProjectedPulseInverseContract(epsilon=0.2, m=2)
    v, n, n_prime, K, forcing, exact, exact_prime = _manufactured_path(contract, 129)

    pressure = contract.pressure(n, n_prime, K, forcing, exact)
    np.testing.assert_allclose(pressure, 0.0, atol=3.0e-15)
    residual = contract.full_equation_residual(
        exact_prime, n, n_prime, K, forcing, exact
    )
    assert float(np.max(np.linalg.norm(residual, axis=-1))) < 3.0e-14

    solved = contract.solve_path(v, n, n_prime, K, forcing)
    # The source pulse starts from zero entrance data.  #281's coefficient
    # contract is applied on the active interior path where t_m is nonzero.
    active_n = n[1:]
    active_t = solved["t_m"][1:]
    C = contract.complete_curl_coefficients(active_n, active_t)
    curl = KokunoSourceCompleteCurlContract(epsilon=contract.epsilon, m=contract.m)
    leading = curl.phase_derivative_leading(active_n, C)
    np.testing.assert_allclose(leading, active_t, rtol=2.0e-10, atol=2.0e-12)
    assert np.all(np.isfinite(C.real)) and np.all(np.isfinite(C.imag))


def test_zero_source_path_and_input_guards_fail_closed():
    contract = KokunoProjectedPulseInverseContract(epsilon=0.25, m=1)
    path = np.linspace(0.0, 1.0, 17)
    n = np.tile(np.array([1.0, 0.2, -0.1]), (path.size, 1))
    n_prime = np.zeros_like(n)
    K = np.zeros((path.size, 3, 3))
    forcing = np.zeros((path.size, 3), dtype=np.complex128)
    solved = contract.solve_path(path, n, n_prime, K, forcing)
    np.testing.assert_array_equal(solved["t_m"], 0.0)
    np.testing.assert_array_equal(solved["pi_m"], 0.0)
    assert solved["max_constraint_defect_abs"] == 0.0
    assert solved["max_projection_correction_norm"] == 0.0

    with pytest.raises(ValueError, match="strictly increasing"):
        contract.solve_path(path[::-1], n, n_prime, K, forcing)
    with pytest.raises(ValueError, match="shape"):
        contract.solve_path(path, n[:-1], n_prime, K, forcing)
    with pytest.raises(ValueError, match="t0 must satisfy"):
        contract.solve_path(path, n, n_prime, K, forcing, t0=n[0])


def test_projected_pulse_payload_round_trip_and_truth_boundary(tmp_path):
    contract = KokunoProjectedPulseInverseContract(epsilon=0.125, m=-3)
    path = contract.save_json(tmp_path / "projected_pulse_inverse.json")
    loaded = KokunoProjectedPulseInverseContract.load_json(path)
    assert loaded.sha256 == contract.sha256
    assert loaded.k_m == loaded.k * loaded.m

    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["source_projected_pulse_equation_executable"] is True
    assert truth["source_complete_curl_coefficient_bridge_executable"] is True
    assert truth["source_actual_pulse_forcing_instantiated"] is False
    assert truth["source_pulse_inverse_reconstructed"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False
    assert "not directly comparable" in payload["routing"]["st006_comparison"]

    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoProjectedPulseInverseContract.from_payload(payload)
