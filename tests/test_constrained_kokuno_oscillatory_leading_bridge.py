import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from openai_ns_reconstruction.kokuno_oscillatory_leading_bridge import (
    KokunoAgent1LeadingChartPhaseBridge,
)


def candidate():
    return KokunoLeadingCoreSeriesCandidate(
        h=0.005,
        j0=0.02,
        sigma=0.5,
        Lambda=10.0,
        C=2.0,
        maxdegree=14,
        eta_nodes=257,
        quadrature_points=16,
    )


def test_chart_scaling_matches_public_leading_velocity():
    c = candidate()
    bridge = KokunoAgent1LeadingChartPhaseBridge(ell=4)
    R = np.array([0.42, 0.50, 0.58])
    Z = np.array([-0.05, 0.0, 0.05])
    T = np.array([0.78, 0.82, 0.86])
    fields = bridge.leading_chart_fields(c, R, Z, T)

    Q = 2.0 ** (-bridge.ell)
    A = 0.5 + c.h
    D = 0.5 - c.h
    theta = np.array([0.2, -0.5, 0.9])
    r = np.sqrt(Q) * R
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    z = (Q**D) * Z
    t = 1.0 - Q * T
    u = c.velocity(x, y, z, t)
    u_theta = -np.sin(theta) * u[..., 0] + np.cos(theta) * u[..., 1]

    np.testing.assert_allclose((Q**A) * u_theta, fields["V"], rtol=2e-11, atol=2e-12)
    np.testing.assert_allclose((Q**A) * u[..., 2], fields["G"], rtol=2e-11, atol=2e-12)
    np.testing.assert_allclose(fields["F_chart"], fields["V"] / R, rtol=0, atol=2e-15)


def test_chart_R_Z_derivatives_match_independent_finite_differences():
    c = candidate()
    bridge = KokunoAgent1LeadingChartPhaseBridge(ell=32)
    R = np.array([0.40, 0.48, 0.56])
    Z = np.array([-0.08, 0.02, 0.09])
    T = np.array([0.80, 0.84, 0.88])
    analytic = bridge.leading_chart_fields(c, R, Z, T)

    errors = []
    for h in (2.0e-3, 1.0e-3, 5.0e-4):
        plus_R = bridge.leading_chart_fields(c, R + h, Z, T)
        minus_R = bridge.leading_chart_fields(c, R - h, Z, T)
        plus_Z = bridge.leading_chart_fields(c, R, Z + h, T)
        minus_Z = bridge.leading_chart_fields(c, R, Z - h, T)
        fd = {
            "V_R": (plus_R["V"] - minus_R["V"]) / (2.0 * h),
            "G_R": (plus_R["G"] - minus_R["G"]) / (2.0 * h),
            "V_Z": (plus_Z["V"] - minus_Z["V"]) / (2.0 * h),
            "G_Z": (plus_Z["G"] - minus_Z["G"]) / (2.0 * h),
        }
        errors.append(
            max(float(np.max(np.abs(fd[name] - analytic[name]))) for name in fd)
        )

    assert errors[-1] < 2.0e-4
    assert errors[-1] < errors[0]


def test_source_phase_covector_composes_on_leading_chart_and_fd_checks():
    c = candidate()
    bridge = KokunoAgent1LeadingChartPhaseBridge(ell=64)
    contract = bridge.phase_contract(c, p=2, p_z=-0.8, x0=0.7, m=1)
    R = np.array([0.40, 0.48, 0.56])
    Z = np.array([-0.06, 0.01, 0.07])
    T = np.array([0.82, 0.84, 0.86])
    theta = np.array([0.2, -0.3, 0.6])
    pulse_v = np.array([0.3, 0.5, 0.7])
    out = bridge.evaluate_phase_on_leading(c, contract, R, theta, Z, T, pulse_v)

    def phase(Rq, thetaq, Zq):
        fields = bridge.leading_chart_fields(c, Rq, Zq, T)
        return contract.phase_from_V(
            Rq, thetaq, Zq, pulse_v, fields["V"], fields["G"]
        )

    errors = []
    for h in (2.0e-3, 1.0e-3, 5.0e-4):
        dR = (phase(R + h, theta, Z) - phase(R - h, theta, Z)) / (2.0 * h)
        dtheta = (phase(R, theta + h, Z) - phase(R, theta - h, Z)) / (2.0 * h)
        dZ = (phase(R, theta, Z + h) - phase(R, theta, Z - h)) / (2.0 * h)
        fd_n = np.stack((dR, dtheta / R, contract.epsilon * dZ), axis=-1)
        errors.append(float(np.max(np.abs(fd_n - out["n_Phi"]))))

    assert errors[-1] < 2.0e-4
    assert errors[-1] < errors[0]


def test_bounds_truth_boundary_and_fail_closed_phase_epsilon(tmp_path):
    c = candidate()
    bridge = KokunoAgent1LeadingChartPhaseBridge(ell=32)
    contract = bridge.phase_contract(c)
    assert contract.epsilon == pytest.approx(bridge.epsilon(c), rel=2e-15)

    wrong = type(contract)(epsilon=0.25)
    with pytest.raises(ValueError, match="epsilon"):
        bridge.evaluate_phase_on_leading(c, wrong, 0.5, 0.0, 0.0, 0.8, 0.3)
    with pytest.raises(ValueError, match="R>0"):
        bridge.leading_chart_fields(c, 0.0, 0.0, 0.8)
    with pytest.raises(ValueError, match="T>0"):
        bridge.leading_chart_fields(c, 0.5, 0.0, 0.0)
    with pytest.raises(ValueError):
        KokunoAgent1LeadingChartPhaseBridge(ell=0)

    path = bridge.save_json(tmp_path / "bridge.json", c)
    payload = json.loads(path.read_text())
    truth = payload["truth_boundary"]
    assert truth["agent1_leading_profile_to_source_leading_V_G_mapping_completed"] is True
    assert truth["source_phase_composed_with_agent1_leading_profile"] is True
    assert truth["actual_background_V_G_mapping_completed"] is False
    assert truth["positive_order_background_corrections_included"] is False
    assert truth["complete_curl_velocity_changed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert payload["sha256"] == bridge.sha256(c)
