import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from openai_ns_reconstruction.kokuno_source_dz_pulse_bridge import (
    KokunoLeadingSourceDzPulseBridge,
)
from openai_ns_reconstruction.kokuno_source_pulse_inverse import (
    KokunoProjectedPulseInverseContract,
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


def setup_contract():
    c = candidate()
    bridge = KokunoLeadingSourceDzPulseBridge(ell=64, dz_step=2.0e-4)
    phase = bridge.leading_bridge.phase_contract(c, p=2, p_z=-0.8, x0=0.7, m=1)
    return c, bridge, phase


def test_source_n_prime_and_K_are_bound_on_fixed_slow_chart():
    c, bridge, phase = setup_contract()
    v = np.linspace(0.15, 0.85, 17)
    geometry = bridge.source_geometry(c, phase, 0.48, 0.025, 0.84, v)

    # The corrected reader gives n_Phi' as the v derivative at fixed slow chart.
    independent_v_derivative = np.gradient(geometry["n_Phi"], v, axis=0, edge_order=2)
    np.testing.assert_allclose(
        independent_v_derivative,
        geometry["n_Phi_prime"],
        rtol=2.0e-11,
        atol=2.0e-11,
    )

    # K is a background object and therefore constant along this fixed-slow-point v path.
    np.testing.assert_allclose(
        geometry["K"],
        np.broadcast_to(geometry["K"][0], geometry["K"].shape),
        rtol=0.0,
        atol=2.0e-14,
    )
    F = geometry["F"]
    F_R = geometry["F_R"]
    G_R = geometry["G_R"]
    np.testing.assert_allclose(geometry["K"][:, 0, 1], -2.0 * F, rtol=0.0, atol=1.0e-14)
    np.testing.assert_allclose(
        geometry["K"][:, 1, 0], 2.0 * F + 0.48 * F_R, rtol=0.0, atol=1.0e-14
    )
    np.testing.assert_allclose(geometry["K"][:, 2, 0], G_R, rtol=0.0, atol=1.0e-14)


def _forcing(v):
    base_direction = np.array([1.0 + 0.2j, -0.7 + 0.1j, 0.5 - 0.3j])
    dz_direction = np.array([-0.3 + 0.4j, 0.6 - 0.2j, 0.2 + 0.5j])
    f = (0.004 + 0.002 * v)[:, None] * base_direction[None, :]
    D_z_f = (0.0004 + 0.0002 * v)[:, None] * dz_direction[None, :]
    return f, D_z_f


def _parent_C(parent, geometry, forcing):
    solved = parent.solve_path(
        geometry["pulse_v"],
        geometry["n_Phi"],
        geometry["n_Phi_prime"],
        geometry["K"],
        forcing,
    )
    t = solved["t_m"]
    C = np.zeros_like(t)
    active = np.linalg.norm(t, axis=-1) > 0.0
    C[active] = parent.complete_curl_coefficients(geometry["n_Phi"][active], t[active])
    return C


def test_Dz_C_m_matches_independent_shifted_parent_solves_and_refines():
    c, bridge, phase = setup_contract()
    v = np.linspace(0.1, 0.9, 65)
    R, Z0, T = 0.48, 0.025, 0.84
    f0, D_z_f = _forcing(v)
    solved = bridge.solve_Dz_coefficient(c, phase, R, Z0, T, v, f0, D_z_f)
    reference = solved["D_z_C_m"]
    parent = KokunoProjectedPulseInverseContract(epsilon=phase.epsilon, m=phase.m)

    errors = []
    for delta in (2.0e-3, 1.0e-3, 5.0e-4):
        plus_geometry = bridge.source_geometry(c, phase, R, Z0 + delta, T, v)
        minus_geometry = bridge.source_geometry(c, phase, R, Z0 - delta, T, v)
        # By construction D_z f = epsilon partial_Z f for this independent check.
        f_plus = f0 + (delta / phase.epsilon) * D_z_f
        f_minus = f0 - (delta / phase.epsilon) * D_z_f
        C_plus = _parent_C(parent, plus_geometry, f_plus)
        C_minus = _parent_C(parent, minus_geometry, f_minus)
        fd = phase.epsilon * (C_plus - C_minus) / (2.0 * delta)
        scale = max(1.0e-12, float(np.linalg.norm(reference[1:])))
        errors.append(float(np.linalg.norm(fd[1:] - reference[1:]) / scale))

    assert errors[-1] < 2.0e-5
    assert errors[1] < 0.4 * errors[0]
    assert errors[2] < 0.4 * errors[1]
    assert np.max(solved["pulse_constraint_defect_abs"]) < 1.0e-10
    assert np.max(solved["sensitivity_constraint_defect_abs"]) < 1.0e-10
    assert np.linalg.norm(reference[1:]) > 0.0


def test_source_Dz_geometry_matches_separate_centered_difference():
    c, bridge, phase = setup_contract()
    v = np.linspace(0.2, 0.8, 13)
    R, Z0, T = 0.5, -0.03, 0.86
    geometry = bridge.source_geometry(c, phase, R, Z0, T, v)

    # Independent lower-order check: deliberately use a different step and stencil.
    delta = 7.5e-5
    plus = bridge._raw_geometry(c, phase, R, Z0 + delta, T, v)
    minus = bridge._raw_geometry(c, phase, R, Z0 - delta, T, v)
    for name, dz_name in (
        ("n_Phi", "D_z_n_Phi"),
        ("n_Phi_prime", "D_z_n_Phi_prime"),
        ("K", "D_z_K"),
    ):
        independent = phase.epsilon * (plus[name] - minus[name]) / (2.0 * delta)
        np.testing.assert_allclose(independent, geometry[dz_name], rtol=2.0e-6, atol=2.0e-8)


def test_fail_closed_provenance_and_bounds(tmp_path):
    c, bridge, phase = setup_contract()
    v = np.linspace(0.2, 0.8, 9)
    with pytest.raises(ValueError, match="strictly increasing"):
        bridge.source_geometry(c, phase, 0.5, 0.0, 0.84, v[::-1])
    with pytest.raises(ValueError, match="R>0"):
        bridge.source_geometry(c, phase, 0.0, 0.0, 0.84, v)
    with pytest.raises(ValueError):
        KokunoLeadingSourceDzPulseBridge(ell=64, dz_step=1.0e-8)

    wrong_phase = type(phase)(p=phase.p, p_z=phase.p_z, x0=phase.x0, epsilon=0.25, m=1)
    with pytest.raises(ValueError, match="epsilon"):
        bridge.source_geometry(c, wrong_phase, 0.5, 0.0, 0.84, v)

    path = bridge.save_json(tmp_path / "dz_bridge.json", c, phase)
    payload = json.loads(path.read_text())
    truth = payload["truth_boundary"]
    assert truth["source_D_z_operator_bound_to_leading_background"] is True
    assert truth["leading_only_D_z_C_m_executable_with_caller_supplied_forcing_derivative"] is True
    assert truth["source_actual_pulse_forcing_instantiated"] is False
    assert truth["source_actual_background_path_instantiated"] is False
    assert truth["source_actual_D_z_pulse_path_instantiated"] is False
    assert truth["actual_background_V_G_mapping_completed"] is False
    assert truth["positive_order_background_corrections_included"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert payload["sha256"] == bridge.sha256(c, phase)
    loaded = type(bridge).load_json(path, c, phase)
    assert loaded == bridge

    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        type(bridge).from_payload(payload, c, phase)
