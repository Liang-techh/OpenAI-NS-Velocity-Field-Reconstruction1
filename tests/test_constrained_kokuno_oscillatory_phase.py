import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract


def base_values(R, Z):
    # Manufactured smooth V=R*F with nontrivial R/Z dependence.
    F = 1.0 + 0.2 * R + 0.1 * Z + 0.03 * R * Z
    F_R = 0.2 + 0.03 * Z
    F_Z = 0.1 + 0.03 * R
    V = R * F
    V_R = F + R * F_R
    V_Z = R * F_Z
    G = -0.2 + 0.3 * R - 0.15 * Z + 0.04 * R * Z
    G_R = 0.3 + 0.04 * Z
    G_Z = -0.15 + 0.04 * R
    return V, V_R, V_Z, G, G_R, G_Z


def test_source_phase_covector_and_quotient_identity():
    c = KokunoOscillatoryPhaseContract(p=2, p_z=-0.75, x0=1.25, epsilon=0.2, m=-1)
    R = np.array([0.7, 1.1, 1.6])
    Z = np.array([-0.2, 0.05, 0.3])
    theta = np.array([0.1, -0.4, 1.2])
    pulse_v = np.array([0.3, 0.5, 0.8])
    V, V_R, V_Z, G, G_R, G_Z = base_values(R, Z)
    out = c.evaluate_from_V(R, theta, Z, pulse_v, V, V_R, V_Z, G, G_R, G_Z)

    F = 1.0 + 0.2 * R + 0.1 * Z + 0.03 * R * Z
    F_R = 0.2 + 0.03 * Z
    F_Z = 0.1 + 0.03 * R
    H = c.p * F + c.p_z * G
    H_R = c.p * F_R + c.p_z * G_R
    H_Z = c.p * F_Z + c.p_z * G_Z
    expected_phi = c.p * theta + c.p_z * Z / c.epsilon + c.x0 * R - pulse_v * H
    expected_n = np.stack(
        (c.x0 - pulse_v * H_R, c.p / R, c.p_z - c.epsilon * pulse_v * H_Z), axis=-1
    )

    np.testing.assert_allclose(out["F"], F, rtol=0, atol=2e-15)
    np.testing.assert_allclose(out["F_R"], F_R, rtol=0, atol=3e-15)
    np.testing.assert_allclose(out["F_Z"], F_Z, rtol=0, atol=2e-15)
    np.testing.assert_allclose(out["Phi"], expected_phi, rtol=0, atol=2e-15)
    np.testing.assert_allclose(out["n_Phi"], expected_n, rtol=0, atol=3e-15)
    assert out["k"] == 3
    assert out["k_m"] == -3


def test_covector_matches_independent_phase_finite_differences():
    c = KokunoOscillatoryPhaseContract(p=3, p_z=0.8, x0=-0.7, epsilon=0.17, m=1)
    R = np.array([0.8, 1.2, 1.7])
    Z = np.array([-0.25, 0.1, 0.35])
    theta = np.array([0.2, -0.1, 0.7])
    pulse_v = np.array([0.4, 0.55, 0.7])
    V, V_R, V_Z, G, G_R, G_Z = base_values(R, Z)
    out = c.evaluate_from_V(R, theta, Z, pulse_v, V, V_R, V_Z, G, G_R, G_Z)

    def phase(Rq, thetaq, Zq):
        Vq, _, _, Gq, _, _ = base_values(Rq, Zq)
        return c.phase_from_V(Rq, thetaq, Zq, pulse_v, Vq, Gq)

    errors = []
    for h in (2.0e-3, 1.0e-3, 5.0e-4):
        dR = (phase(R + h, theta, Z) - phase(R - h, theta, Z)) / (2.0 * h)
        dtheta = (phase(R, theta + h, Z) - phase(R, theta - h, Z)) / (2.0 * h)
        dZ = (phase(R, theta, Z + h) - phase(R, theta, Z - h)) / (2.0 * h)
        fd_n = np.stack((dR, dtheta / R, c.epsilon * dZ), axis=-1)
        errors.append(float(np.max(np.abs(fd_n - out["n_Phi"]))))
    assert errors[-1] < 2.0e-9
    assert errors[-1] <= errors[0] + 1e-12


def test_batch_axis_guard_bounds_and_nonfinite_rejection():
    c = KokunoOscillatoryPhaseContract()
    R = np.array([[0.6], [1.1]])
    Z = np.array([[-0.2, 0.2]])
    V, V_R, V_Z, G, G_R, G_Z = base_values(R, Z)
    out = c.evaluate_from_V(R, 0.3, Z, 0.4, V, V_R, V_Z, G, G_R, G_Z)
    assert out["Phi"].shape == (2, 2)
    assert out["n_Phi"].shape == (2, 2, 3)
    assert np.all(out["n_Phi_norm_sq"] > 0)

    with pytest.raises(ValueError, match="R>0"):
        c.phase_from_V(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    with pytest.raises(ValueError, match="finite"):
        c.phase_from_V(1.0, np.nan, 0.0, 0.0, 1.0, 0.0)
    with pytest.raises(ValueError):
        KokunoOscillatoryPhaseContract(p=0)
    with pytest.raises(ValueError):
        KokunoOscillatoryPhaseContract(epsilon=0.0)
    with pytest.raises(ValueError):
        KokunoOscillatoryPhaseContract(m=0)


def test_truth_boundary_and_fail_closed_roundtrip(tmp_path):
    c = KokunoOscillatoryPhaseContract(p=2, p_z=1.25, x0=-0.5, epsilon=0.125, m=2)
    path = c.save_json(tmp_path / "phase.json")
    loaded = KokunoOscillatoryPhaseContract.load_json(path)
    assert loaded.sha256 == c.sha256
    payload = json.loads(path.read_text())
    truth = payload["truth_boundary"]
    assert truth["source_phase_formula_executable"] is True
    assert truth["source_covector_formula_executable"] is True
    assert truth["agent1_profile_to_source_V_G_mapping_completed"] is False
    assert truth["complete_curl_velocity_changed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    payload["truth_boundary"]["paper_exact"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoOscillatoryPhaseContract.load_json(path)
