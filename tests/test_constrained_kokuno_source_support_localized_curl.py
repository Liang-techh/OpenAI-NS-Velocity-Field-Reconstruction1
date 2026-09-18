from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_complete_curl import (
    KokunoSourceCompleteCurlContract,
)
from openai_ns_reconstruction.kokuno_source_support_localized_curl import (
    KokunoSourceSupportLocalizedCurl,
)


def _phase_and_covector(contract, R, theta, Z):
    R, theta, Z = np.broadcast_arrays(
        np.asarray(R, dtype=float), np.asarray(theta, dtype=float), np.asarray(Z, dtype=float)
    )
    p = 1.5
    p_z = 0.65
    x0 = 1.1
    pulse = 0.28
    a, b, c = 0.09, -0.06, 0.07
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
    R, theta, Z = np.broadcast_arrays(
        np.asarray(R, dtype=float), np.asarray(theta, dtype=float), np.asarray(Z, dtype=float)
    )
    # Source slow amplitudes are independent of the fast angular phase.  Keep
    # this manufactured check in that same class so an independent theta
    # derivative tests only exp(i*k_m*Phi), exactly as the complete-curl
    # contract assumes.
    raw = np.stack(
        (
            0.75 + 0.08 * R,
            -0.31 + 0.05 * Z,
            0.22 + 0.03 * R * Z,
        ),
        axis=-1,
    ).astype(np.complex128)
    n2 = np.sum(n * n, axis=-1)
    return raw - n * (np.sum(n * raw, axis=-1) / n2)[..., None]


def _coefficient(contract, R, theta, Z):
    _, n = _phase_and_covector(contract, R, theta, Z)
    return contract.coefficient(n, _transverse_t(contract, R, theta, Z))


def _fd4(fn, x, h):
    return (fn(x - 2 * h) - 8 * fn(x - h) + 8 * fn(x + h) - fn(x + 2 * h)) / (12 * h)


def _coefficient_derivatives(contract, R, theta, Z, h=1.5e-5):
    D_r = _fd4(lambda RR: _coefficient(contract, RR, theta, Z), R, h)
    D_z = contract.epsilon * _fd4(lambda ZZ: _coefficient(contract, R, theta, ZZ), Z, h)
    return D_r, D_z


def _eta(contract, R, Z):
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    return 0.55 + 0.08 * R - 0.05 * Z + 0.03 * R * Z


def _eta_derivatives(contract, R, Z):
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    D_r = 0.08 + 0.03 * Z
    D_z = contract.epsilon * (-0.05 + 0.03 * R)
    return D_r, D_z


def _localized_public(localizer, R, theta, Z):
    contract = localizer.complete_curl
    phase, n = _phase_and_covector(contract, R, theta, Z)
    t = _transverse_t(contract, R, theta, Z)
    D_r_C, D_z_C = _coefficient_derivatives(contract, R, theta, Z)
    D_r_eta, D_z_eta = _eta_derivatives(contract, R, Z)
    return localizer.localized_mode(
        R,
        phase,
        n,
        t,
        D_r_C,
        D_z_C,
        _eta(contract, R, Z),
        D_r_eta,
        D_z_eta,
    )["velocity"]


def _localized_potential(localizer, R, theta, Z):
    contract = localizer.complete_curl
    phase, _ = _phase_and_covector(contract, R, theta, Z)
    C = _coefficient(contract, R, theta, Z)
    return _eta(contract, R, Z) * C * np.exp(1j * contract.k_m * phase)


def _independent_curl(localizer, R, theta, Z, h):
    contract = localizer.complete_curl
    A = _localized_potential(localizer, R, theta, Z)
    d_theta = _fd4(lambda th: _localized_potential(localizer, R, th, Z), theta, h)
    D_r = _fd4(lambda RR: _localized_potential(localizer, RR, theta, Z), R, h)
    D_z = contract.epsilon * _fd4(
        lambda ZZ: _localized_potential(localizer, R, theta, ZZ), Z, h
    )
    return np.array(
        [
            d_theta[2] / R - D_z[1],
            D_z[0] - D_r[2],
            D_r[1] + A[1] / R - d_theta[0] / R,
        ],
        dtype=np.complex128,
    )


def _independent_divergence(localizer, R, theta, Z, h):
    contract = localizer.complete_curl
    u = _localized_public(localizer, R, theta, Z)
    D_r_ur = _fd4(lambda RR: _localized_public(localizer, RR, theta, Z)[0], R, h)
    d_theta_utheta = _fd4(
        lambda th: _localized_public(localizer, R, th, Z)[1], theta, h
    )
    D_z_uz = contract.epsilon * _fd4(
        lambda ZZ: _localized_public(localizer, R, theta, ZZ)[2], Z, h
    )
    return D_r_ur + u[0] / R + d_theta_utheta / R + D_z_uz


def test_localized_mode_matches_independent_curl_and_divergence_refines():
    localizer = KokunoSourceSupportLocalizedCurl(epsilon=0.2, m=1)
    probes = [(0.74, 0.27, -0.11), (1.06, -0.39, 0.16)]
    steps = [0.02, 0.01, 0.005]
    curl_errors = []
    divergence_errors = []
    for h in steps:
        ce = []
        de = []
        for R, theta, Z in probes:
            expected = _localized_public(localizer, R, theta, Z)
            independent = _independent_curl(localizer, R, theta, Z, h)
            ce.append(np.linalg.norm(independent - expected))
            de.append(abs(_independent_divergence(localizer, R, theta, Z, h)))
        curl_errors.append(max(ce))
        divergence_errors.append(max(de))
    assert curl_errors[2] < 3.0e-7
    assert divergence_errors[2] < 3.0e-7
    assert curl_errors[0] / curl_errors[1] > 7.0
    assert curl_errors[1] / curl_errors[2] > 7.0
    assert divergence_errors[0] / divergence_errors[1] > 6.0
    assert divergence_errors[1] / divergence_errors[2] > 6.0


def test_support_gradient_terms_are_retained_and_not_post_velocity_cutoff():
    localizer = KokunoSourceSupportLocalizedCurl(epsilon=0.2, m=1)
    contract = localizer.complete_curl
    R, theta, Z = 0.91, 0.37, 0.12
    phase, n = _phase_and_covector(contract, R, theta, Z)
    t = _transverse_t(contract, R, theta, Z)
    D_r_C, D_z_C = _coefficient_derivatives(contract, R, theta, Z)
    eta = _eta(contract, R, Z)
    D_r_eta, D_z_eta = _eta_derivatives(contract, R, Z)
    out = localizer.localized_mode(
        R, phase, n, t, D_r_C, D_z_C, eta, D_r_eta, D_z_eta
    )
    base = contract.mode_velocity(R, phase, n, t, D_r_C, D_z_C)["velocity"]
    naive = eta * base
    assert np.linalg.norm(out["velocity"] - naive) > 1.0e-4
    np.testing.assert_allclose(out["leading_local"], eta * t, rtol=2e-13, atol=2e-13)


def test_exact_support_zero_and_batch_product_rule():
    localizer = KokunoSourceSupportLocalizedCurl(epsilon=0.25, m=-1)
    contract = localizer.complete_curl
    R = np.array([0.76, 0.94, 1.12])
    theta = np.array([0.1, -0.2, 0.4])
    Z = np.array([-0.1, 0.03, 0.14])
    phase, n = _phase_and_covector(contract, R, theta, Z)
    t = _transverse_t(contract, R, theta, Z)
    D_r_C = np.stack([
        _coefficient_derivatives(contract, r, th, z)[0]
        for r, th, z in zip(R, theta, Z)
    ])
    D_z_C = np.stack([
        _coefficient_derivatives(contract, r, th, z)[1]
        for r, th, z in zip(R, theta, Z)
    ])
    zero = localizer.localized_mode(
        R, phase, n, t, D_r_C, D_z_C, 0.0, 0.0, 0.0
    )
    np.testing.assert_array_equal(zero["vector_potential"], np.zeros((3, 3), dtype=np.complex128))
    np.testing.assert_array_equal(zero["velocity"], np.zeros((3, 3), dtype=np.complex128))

    eta = _eta(contract, R, Z)
    D_r_eta, D_z_eta = _eta_derivatives(contract, R, Z)
    data = localizer.localized_coefficient_data(
        n, t, D_r_C, D_z_C, eta, D_r_eta, D_z_eta
    )
    C = contract.coefficient(n, t)
    np.testing.assert_allclose(data["C_local"], eta[:, None] * C)
    np.testing.assert_allclose(
        data["D_r_C_local"], eta[:, None] * D_r_C + D_r_eta[:, None] * C
    )
    np.testing.assert_allclose(
        data["D_z_C_local"], eta[:, None] * D_z_C + D_z_eta[:, None] * C
    )
    with pytest.raises(ValueError, match="squared partition"):
        localizer.localized_coefficient_data(n, t, D_r_C, D_z_C, 1.01, 0.0, 0.0)


def test_truth_boundary_and_serialization_fail_closed(tmp_path):
    localizer = KokunoSourceSupportLocalizedCurl(epsilon=0.125, m=2)
    path = localizer.save_json(tmp_path / "source_support_localized_curl.json")
    loaded = KokunoSourceSupportLocalizedCurl.load_json(path)
    assert loaded.sha256 == localizer.sha256
    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["source_squared_partition_structure_executable"] is True
    assert truth["vector_potential_level_localization_executable"] is True
    assert truth["support_gradient_complete_curl_terms_retained"] is True
    assert truth["concrete_source_partition_bumps_reconstructed"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False
    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        KokunoSourceSupportLocalizedCurl.from_payload(payload)
