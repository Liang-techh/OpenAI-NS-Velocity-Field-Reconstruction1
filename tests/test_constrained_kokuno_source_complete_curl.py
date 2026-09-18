from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_complete_curl import (
    KokunoSourceCompleteCurlContract,
)


def _phase_and_covector(contract, R, theta, Z):
    R, theta, Z = np.broadcast_arrays(
        np.asarray(R, dtype=float), np.asarray(theta, dtype=float), np.asarray(Z, dtype=float)
    )
    p = 2.0
    p_z = 0.7
    x0 = 1.2
    pulse = 0.35
    a, b, c = 0.11, -0.07, 0.09
    H = a * R * R + b * R * Z + c * Z * Z
    H_R = 2.0 * a * R + b * Z
    H_Z = b * R + 2.0 * c * Z
    phase = p * theta + p_z * Z / contract.epsilon + x0 * R - pulse * H
    n = np.stack(
        (
            x0 - pulse * H_R,
            p / R,
            p_z - contract.epsilon * pulse * H_Z,
        ),
        axis=-1,
    )
    return phase, n


def _transverse_t(contract, R, theta, Z):
    _, n = _phase_and_covector(contract, R, theta, Z)
    R, Z = np.broadcast_arrays(np.asarray(R, dtype=float), np.asarray(Z, dtype=float))
    raw = np.stack(
        (
            0.8 + 0.10 * R,
            -0.35 + 0.06 * Z,
            0.25 + 0.04 * R * Z,
        ),
        axis=-1,
    ).astype(np.complex128)
    n2 = np.sum(n * n, axis=-1)
    return raw - n * (np.sum(n * raw, axis=-1) / n2)[..., None]


def _coefficient_fn(contract, R, theta, Z):
    _, n = _phase_and_covector(contract, R, theta, Z)
    return contract.coefficient(n, _transverse_t(contract, R, theta, Z))


def _fd4(fn, x, h):
    return (fn(x - 2 * h) - 8 * fn(x - h) + 8 * fn(x + h) - fn(x + 2 * h)) / (12 * h)


def _coefficient_derivatives(contract, R, theta, Z, h=2.0e-5):
    D_r = _fd4(lambda RR: _coefficient_fn(contract, RR, theta, Z), R, h)
    partial_Z = _fd4(lambda ZZ: _coefficient_fn(contract, R, theta, ZZ), Z, h)
    return D_r, contract.epsilon * partial_Z


def _potential(contract, R, theta, Z):
    phase, _ = _phase_and_covector(contract, R, theta, Z)
    C = _coefficient_fn(contract, R, theta, Z)
    return contract.mode_vector_potential(phase, C)


def _public_velocity(contract, R, theta, Z):
    phase, n = _phase_and_covector(contract, R, theta, Z)
    t = _transverse_t(contract, R, theta, Z)
    D_r_C, D_z_C = _coefficient_derivatives(contract, R, theta, Z)
    return contract.mode_velocity(R, phase, n, t, D_r_C, D_z_C)["velocity"]


def _independent_cylindrical_curl(contract, R, theta, Z, h):
    A = _potential(contract, R, theta, Z)
    d_theta = _fd4(lambda th: _potential(contract, R, th, Z), theta, h)
    D_r = _fd4(lambda RR: _potential(contract, RR, theta, Z), R, h)
    D_z = contract.epsilon * _fd4(lambda ZZ: _potential(contract, R, theta, ZZ), Z, h)
    return np.array(
        [
            d_theta[2] / R - D_z[1],
            D_z[0] - D_r[2],
            D_r[1] + A[1] / R - d_theta[0] / R,
        ],
        dtype=np.complex128,
    )


def _independent_divergence(contract, R, theta, Z, h):
    u = _public_velocity(contract, R, theta, Z)
    D_r_ur = _fd4(lambda RR: _public_velocity(contract, RR, theta, Z)[0], R, h)
    d_theta_utheta = _fd4(
        lambda th: _public_velocity(contract, R, th, Z)[1], theta, h
    )
    D_z_uz = contract.epsilon * _fd4(
        lambda ZZ: _public_velocity(contract, R, theta, ZZ)[2], Z, h
    )
    return D_r_ur + u[0] / R + d_theta_utheta / R + D_z_uz


def test_source_coefficient_recovers_transverse_leading_term_and_rejects_bad_input():
    contract = KokunoSourceCompleteCurlContract(epsilon=0.2, m=2)
    R = np.array([0.7, 0.9, 1.1])
    theta = np.array([0.2, -0.4, 0.8])
    Z = np.array([-0.15, 0.05, 0.2])
    _, n = _phase_and_covector(contract, R, theta, Z)
    t = _transverse_t(contract, R, theta, Z)
    C = contract.coefficient(n, t)
    np.testing.assert_allclose(contract.phase_derivative_leading(n, C), t, rtol=3e-14, atol=3e-14)
    with pytest.raises(ValueError, match="transverse"):
        contract.coefficient(n, t + 0.1 * n)
    with pytest.raises(ValueError, match="R>0"):
        contract.remainder(0.0, C[0], C[0], C[0])


def test_source_complete_mode_matches_independent_cylindrical_curl_and_divergence_refines():
    contract = KokunoSourceCompleteCurlContract(epsilon=0.2, m=1)
    probes = [(0.76, 0.31, -0.12), (1.08, -0.44, 0.17)]
    steps = [0.02, 0.01, 0.005]
    curl_errors = []
    div_errors = []
    for h in steps:
        ce = []
        de = []
        for R, theta, Z in probes:
            expected = _public_velocity(contract, R, theta, Z)
            independent = _independent_cylindrical_curl(contract, R, theta, Z, h)
            ce.append(np.linalg.norm(independent - expected))
            de.append(abs(_independent_divergence(contract, R, theta, Z, h)))
        curl_errors.append(max(ce))
        div_errors.append(max(de))
    assert curl_errors[2] < 2.0e-7
    assert div_errors[2] < 2.0e-7
    assert curl_errors[1] / curl_errors[2] > 8.0
    assert curl_errors[0] / curl_errors[1] > 8.0
    assert div_errors[1] / div_errors[2] > 7.0
    assert div_errors[0] / div_errors[1] > 7.0


def test_batch_complete_amplitude_keeps_longitudinal_remainder_and_shapes():
    contract = KokunoSourceCompleteCurlContract(epsilon=0.3, m=-1)
    R = np.array([0.72, 0.95, 1.18])
    theta = np.array([0.1, 0.3, -0.2])
    Z = np.array([-0.12, 0.02, 0.16])
    phase, n = _phase_and_covector(contract, R, theta, Z)
    t = _transverse_t(contract, R, theta, Z)
    D_r_C = np.stack([_coefficient_derivatives(contract, r, th, z)[0] for r, th, z in zip(R, theta, Z)])
    D_z_C = np.stack([_coefficient_derivatives(contract, r, th, z)[1] for r, th, z in zip(R, theta, Z)])
    out = contract.mode_velocity(R, phase, n, t, D_r_C, D_z_C)
    assert out["velocity"].shape == (3, 3)
    assert out["vector_potential"].shape == (3, 3)
    assert np.all(np.isfinite(out["velocity"]))
    assert np.linalg.norm(out["r_m"]) > 1.0e-6
    longitudinal = np.sum(n * out["a_m"], axis=-1)
    # Source explicitly keeps the complete amplitude's longitudinal remainder.
    assert np.max(np.abs(longitudinal)) > 1.0e-8
    np.testing.assert_allclose(np.sum(n * t, axis=-1), 0.0, atol=2e-14)


def test_phase_binding_and_truth_metadata_fail_closed(tmp_path):
    class Phase:
        epsilon = 0.125
        m = 3

    contract = KokunoSourceCompleteCurlContract.from_phase_contract(Phase())
    assert contract.epsilon == 0.125
    assert contract.m == 3
    assert contract.k_m == contract.k * 3
    path = contract.save_json(tmp_path / "source_complete_curl.json")
    loaded = KokunoSourceCompleteCurlContract.load_json(path)
    assert loaded.sha256 == contract.sha256
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["truth_boundary"]["source_complete_curl_remainder_executable"] is True
    assert payload["truth_boundary"]["source_pulse_inverse_reconstructed"] is False
    assert payload["truth_boundary"]["actual_background_V_G_mapping_completed"] is False
    assert payload["truth_boundary"]["public_velocity_correction_materialized"] is False
    assert payload["truth_boundary"]["pde_validated"] is False
    payload["truth_boundary"]["pde_validated"] = True
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoSourceCompleteCurlContract.from_payload(payload)
