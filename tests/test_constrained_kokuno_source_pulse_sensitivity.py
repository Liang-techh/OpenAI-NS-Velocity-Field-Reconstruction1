from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_pulse_sensitivity import (
    KokunoProjectedPulseSensitivityContract,
)


def _manufactured_parameter_path(
    contract: KokunoProjectedPulseSensitivityContract,
    count: int,
    *,
    xi: float = 0.37,
):
    path = np.linspace(0.0, 1.0, count)
    beta = 0.4
    cb = np.cos(beta)
    sb = np.sin(beta)
    e = np.tile(np.array([np.cos(xi), np.sin(xi), 0.0]), (count, 1))
    p = np.tile(np.array([-np.sin(xi), np.cos(xi), 0.0]), (count, 1))
    ez = np.array([0.0, 0.0, 1.0])
    n = e
    n_xi = p
    n_prime = np.zeros_like(n)
    n_prime_xi = np.zeros_like(n)
    h = np.sin(0.7 * path)
    h_prime = 0.7 * np.cos(0.7 * path)
    direction = cb * p + sb * ez
    direction_xi = -cb * e
    exact = h[:, None] * direction
    exact_xi = h[:, None] * direction_xi
    exact_prime = h_prime[:, None] * direction
    exact_prime_xi = h_prime[:, None] * direction_xi

    K = np.zeros((count, 3, 3), dtype=np.complex128)
    K_xi = np.zeros_like(K)
    d = contract.pulse.diffusion(n)
    forcing = -(exact_prime + (contract.m**2 * d)[:, None] * exact)
    forcing_xi = -(exact_prime_xi + (contract.m**2 * d)[:, None] * exact_xi)
    exact_C_xi = (1j / contract.k_m) * h[:, None] * sb * e
    return (
        path,
        n,
        n_prime,
        K,
        forcing,
        n_xi,
        n_prime_xi,
        K_xi,
        forcing_xi,
        exact,
        exact_xi,
        exact_C_xi,
    )


def test_projected_pulse_sensitivity_refines_and_preserves_differentiated_constraint():
    contract = KokunoProjectedPulseSensitivityContract(epsilon=0.2, m=2)
    t_errors = []
    sensitivity_errors = []
    coefficient_errors = []
    for count in (33, 65, 129):
        data = _manufactured_parameter_path(contract, count)
        solved = contract.solve_path(*data[:9])
        exact = data[9]
        exact_xi = data[10]
        exact_C_xi = data[11]
        t_errors.append(float(np.max(np.linalg.norm(solved["t_m"] - exact, axis=-1))))
        sensitivity_errors.append(
            float(np.max(np.linalg.norm(solved["t_m_xi"] - exact_xi, axis=-1)))
        )
        coefficient_errors.append(
            float(np.max(np.linalg.norm(solved["C_m_xi"] - exact_C_xi, axis=-1)))
        )
        assert solved["max_pulse_constraint_defect_abs"] < 2.0e-15
        assert solved["max_sensitivity_constraint_defect_abs"] < 2.0e-15

    assert t_errors[0] / t_errors[1] > 3.8
    assert t_errors[1] / t_errors[2] > 3.8
    assert sensitivity_errors[0] / sensitivity_errors[1] > 3.8
    assert sensitivity_errors[1] / sensitivity_errors[2] > 3.8
    assert coefficient_errors[0] / coefficient_errors[1] > 3.8
    assert coefficient_errors[1] / coefficient_errors[2] > 3.8
    assert t_errors[2] < 1.7e-6
    assert sensitivity_errors[2] < 1.6e-6
    assert coefficient_errors[2] < 1.1e-7


def test_sensitivity_rhs_matches_independent_parameter_difference_of_parent_equation():
    contract = KokunoProjectedPulseSensitivityContract(epsilon=0.2, m=2)
    rng = np.random.default_rng(9173011)
    n = np.array([1.1, -0.4, 0.7])
    n_xi = np.array([0.2, 0.1, -0.3])
    n_prime = np.array([0.1, 0.2, -0.1])
    n_prime_xi = np.array([-0.03, 0.04, 0.02])
    K = rng.normal(size=(3, 3)) + 0.1j * rng.normal(size=(3, 3))
    K_xi = 0.2 * rng.normal(size=(3, 3))
    forcing = rng.normal(size=3) + 0.2j * rng.normal(size=3)
    forcing_xi = 0.2 * rng.normal(size=3) + 0.1j * rng.normal(size=3)
    t_m = rng.normal(size=3) + 0.1j * rng.normal(size=3)
    t_m_xi = rng.normal(size=3) + 0.1j * rng.normal(size=3)

    delta = 1.0e-6
    plus = contract.pulse.projected_rhs(
        n + delta * n_xi,
        n_prime + delta * n_prime_xi,
        K + delta * K_xi,
        forcing + delta * forcing_xi,
        t_m + delta * t_m_xi,
    )
    minus = contract.pulse.projected_rhs(
        n - delta * n_xi,
        n_prime - delta * n_prime_xi,
        K - delta * K_xi,
        forcing - delta * forcing_xi,
        t_m - delta * t_m_xi,
    )
    finite_difference = (plus - minus) / (2.0 * delta)
    analytic = contract.sensitivity_rhs(
        n,
        n_prime,
        K,
        forcing,
        n_xi,
        n_prime_xi,
        K_xi,
        forcing_xi,
        t_m,
        t_m_xi,
    )
    np.testing.assert_allclose(analytic, finite_difference, rtol=2.0e-9, atol=2.0e-9)


def test_complete_curl_coefficient_sensitivity_matches_centered_parameter_difference():
    contract = KokunoProjectedPulseSensitivityContract(epsilon=0.125, m=-3)
    n = np.array([1.2, -0.3, 0.5])
    n_xi = np.array([0.2, 0.4, -0.1])
    raw_t = np.array([0.4 + 0.1j, -0.2j, 0.3 - 0.05j])
    raw_s = np.array([-0.1j, 0.2 + 0.05j, -0.07])
    t = contract.pulse.project_transverse(n, raw_t)
    _, _, _, _, P, P_xi = contract._projector_data(n, n_xi)
    s = P @ raw_s + P_xi @ raw_t
    analytic = contract.coefficient_sensitivity(n, n_xi, t, s)

    delta = 1.0e-6
    n_plus = n + delta * n_xi
    n_minus = n - delta * n_xi
    t_plus = contract.pulse.project_transverse(n_plus, raw_t + delta * raw_s)
    t_minus = contract.pulse.project_transverse(n_minus, raw_t - delta * raw_s)
    C_plus = contract.pulse.complete_curl_coefficients(n_plus, t_plus)
    C_minus = contract.pulse.complete_curl_coefficients(n_minus, t_minus)
    finite_difference = (C_plus - C_minus) / (2.0 * delta)
    np.testing.assert_allclose(analytic, finite_difference, rtol=3.0e-9, atol=3.0e-10)
    assert abs(np.dot(n, s) + np.dot(n_xi, t)) < 3.0e-16


def test_sensitivity_payload_round_trip_truth_boundary_and_guards(tmp_path):
    contract = KokunoProjectedPulseSensitivityContract(epsilon=0.125, m=3)
    path = contract.save_json(tmp_path / "pulse_sensitivity.json")
    loaded = KokunoProjectedPulseSensitivityContract.load_json(path)
    assert loaded.sha256 == contract.sha256

    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["generic_parameter_sensitivity_operator_executable"] is True
    assert truth["complete_curl_coefficient_sensitivity_executable"] is True
    assert truth["source_actual_D_r_path_instantiated"] is False
    assert truth["source_actual_D_z_path_instantiated"] is False
    assert truth["source_pulse_parameter_derivatives_reconstructed"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False
    assert "not directly comparable" in payload["routing"]["st006_comparison"]

    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoProjectedPulseSensitivityContract.from_payload(payload)

    data = list(_manufactured_parameter_path(contract, 17))
    data[5] = data[5][:-1]
    with pytest.raises(ValueError, match="n_phi_xi must have shape"):
        contract.solve_path(*data[:9])
