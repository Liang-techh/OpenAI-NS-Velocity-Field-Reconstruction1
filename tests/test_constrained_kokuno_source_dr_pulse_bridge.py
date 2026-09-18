import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_leading_core_series import KokunoLeadingCoreSeriesCandidate
from openai_ns_reconstruction.kokuno_source_dr_pulse_bridge import (
    KokunoLeadingSourceDrPulseBridge,
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
    bridge = KokunoLeadingSourceDrPulseBridge(ell=64, dr_step=1.0e-4)
    phase = bridge.leading_bridge.phase_contract(c, p=2, p_z=-0.8, x0=0.7, m=1)
    return c, bridge, phase


def _forcing(v):
    base_direction = np.array([1.0 + 0.2j, -0.7 + 0.1j, 0.5 - 0.3j])
    dr_direction = np.array([0.25 - 0.35j, -0.45 + 0.15j, 0.3 + 0.4j])
    f = (0.004 + 0.002 * v)[:, None] * base_direction[None, :]
    D_r_f = (0.0005 + 0.00015 * v)[:, None] * dr_direction[None, :]
    return f, D_r_f


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


def _shifted_parent_C(bridge, candidate, phase, parent, R0, Z, T, v, f0, D_r_f, shift):
    geometry = bridge.source_geometry(candidate, phase, R0 + shift, Z, T, v)
    # This manufactured forcing is auxiliary-independent, so source D_r f is
    # exactly partial_R f and a direct R shift gives an independent comparator.
    forcing = f0 + shift * D_r_f
    return _parent_C(parent, geometry, forcing)


def test_source_Dr_geometry_reduces_to_partial_R_only_on_declared_sector():
    c, bridge, phase = setup_contract()
    v = np.linspace(0.2, 0.8, 13)
    R0, Z, T = 0.5, -0.03, 0.86
    geometry = bridge.source_geometry(c, phase, R0, Z, T, v)

    assert geometry["leading_geometry_auxiliary_independent"] is True
    assert "L_i=0" in geometry["D_r_geometry_reduction"]
    assert "M_i" in geometry["D_r_definition"]

    # Independent comparator uses a different step and stencil from the stored
    # FD4 derivative, so this is not self-comparison of the implementation.
    delta = 6.0e-5
    plus = bridge._raw_geometry(c, phase, R0 + delta, Z, T, v)
    minus = bridge._raw_geometry(c, phase, R0 - delta, Z, T, v)
    for name, dr_name in (
        ("n_Phi", "D_r_n_Phi"),
        ("n_Phi_prime", "D_r_n_Phi_prime"),
        ("K", "D_r_K"),
    ):
        independent = (plus[name] - minus[name]) / (2.0 * delta)
        np.testing.assert_allclose(independent, geometry[dr_name], rtol=2.0e-5, atol=1.0e-7)


def test_Dr_C_m_matches_independent_shifted_parent_solves_and_refines():
    c, bridge, phase = setup_contract()
    v = np.linspace(0.1, 0.9, 65)
    R0, Z, T = 0.48, 0.025, 0.84
    f0, D_r_f = _forcing(v)
    solved = bridge.solve_Dr_coefficient(c, phase, R0, Z, T, v, f0, D_r_f)
    reference = solved["D_r_C_m"]
    parent = KokunoProjectedPulseInverseContract(epsilon=phase.epsilon, m=phase.m)

    errors = []
    for delta in (8.0e-4, 4.0e-4, 2.0e-4):
        C_m2 = _shifted_parent_C(bridge, c, phase, parent, R0, Z, T, v, f0, D_r_f, -2.0 * delta)
        C_m1 = _shifted_parent_C(bridge, c, phase, parent, R0, Z, T, v, f0, D_r_f, -delta)
        C_p1 = _shifted_parent_C(bridge, c, phase, parent, R0, Z, T, v, f0, D_r_f, delta)
        C_p2 = _shifted_parent_C(bridge, c, phase, parent, R0, Z, T, v, f0, D_r_f, 2.0 * delta)
        fd4 = (C_m2 - 8.0 * C_m1 + 8.0 * C_p1 - C_p2) / (12.0 * delta)
        scale = max(1.0e-12, float(np.linalg.norm(reference[1:])))
        errors.append(float(np.linalg.norm(fd4[1:] - reference[1:]) / scale))

    assert errors[-1] < 2.0e-4, f"fourth-order D_r C_m errors={errors}"
    assert errors[1] < 0.3 * errors[0], f"fourth-order D_r C_m errors={errors}"
    assert errors[2] < 0.3 * errors[1], f"fourth-order D_r C_m errors={errors}"
    assert np.max(solved["pulse_constraint_defect_abs"]) < 1.0e-10
    assert np.max(solved["sensitivity_constraint_defect_abs"]) < 1.0e-10
    assert np.linalg.norm(reference[1:]) > 0.0


def test_fail_closed_axis_provenance_and_bounds(tmp_path):
    c, bridge, phase = setup_contract()
    v = np.linspace(0.2, 0.8, 9)
    with pytest.raises(ValueError, match="strictly increasing"):
        bridge.source_geometry(c, phase, 0.5, 0.0, 0.84, v[::-1])
    with pytest.raises(ValueError, match="2\\*dr_step"):
        bridge.source_geometry(c, phase, 1.5e-4, 0.0, 0.84, v)
    with pytest.raises(ValueError):
        KokunoLeadingSourceDrPulseBridge(ell=64, dr_step=1.0e-8)

    wrong_phase = type(phase)(p=phase.p, p_z=phase.p_z, x0=phase.x0, epsilon=0.25, m=1)
    with pytest.raises(ValueError, match="epsilon"):
        bridge.source_geometry(c, wrong_phase, 0.5, 0.0, 0.84, v)

    path = bridge.save_json(tmp_path / "dr_bridge.json", c, phase)
    payload = json.loads(path.read_text())
    truth = payload["truth_boundary"]
    assert truth["source_D_r_operator_decomposed"] is True
    assert truth["leading_geometry_auxiliary_independent"] is True
    assert truth["source_D_r_geometry_reduces_to_partial_R_on_leading_only_bridge"] is True
    assert truth["leading_only_D_r_C_m_executable_with_caller_supplied_forcing_derivative"] is True
    assert truth["source_actual_pulse_forcing_instantiated"] is False
    assert truth["source_actual_background_path_instantiated"] is False
    assert truth["source_actual_D_r_pulse_path_instantiated"] is False
    assert truth["source_actual_auxiliary_dependent_D_r_path_instantiated"] is False
    assert truth["actual_background_V_G_mapping_completed"] is False
    assert truth["positive_order_background_corrections_included"] is False
    assert truth["public_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False
    assert "D_r=partial_R+" in payload["source"]["source_formulas"]["source_D_r"]
    assert payload["sha256"] == bridge.sha256(c, phase)
    loaded = type(bridge).load_json(path, c, phase)
    assert loaded == bridge

    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        type(bridge).from_payload(payload, c, phase)
