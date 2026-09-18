import json
from math import ceil, sqrt

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_oscillatory_phase import KokunoOscillatoryPhaseContract
from openai_ns_reconstruction.kokuno_source_phase_label import KokunoSourcePhaseLabelSelection


def test_source_phase_label_formulas_and_phase_contract_bridge():
    s = KokunoSourcePhaseLabelSelection(
        epsilon=0.25,
        lambda0=1.0,
        u_star=2.0,
        L_s=1.5,
        sigma=-1,
        R0=1.25,
        g0=(1.2, -0.4),
        K=(0.3, 0.9),
    )
    k = ceil(s.epsilon ** -0.5)
    B_s = sqrt(s.lambda0 / (s.epsilon * k * k * (1.0 + s.u_star**2) ** 1.5))
    g0 = np.array(s.g0)
    K = np.array(s.K)
    q = B_s * (K - s.sigma * s.u_star * g0 / (s.L_s * np.dot(g0, g0)))
    tilde_p = s.R0 * q[0]
    j = s.nearest_nonzero_integer(k * tilde_p)

    assert s.k == k
    np.testing.assert_allclose(s.B_s, B_s, rtol=0, atol=2e-15)
    np.testing.assert_allclose(s.x0, s.sigma * B_s * s.u_star / 2.0, rtol=0, atol=2e-15)
    np.testing.assert_allclose(s.q, q, rtol=0, atol=2e-15)
    np.testing.assert_allclose(s.tilde_p, tilde_p, rtol=0, atol=2e-15)
    assert s.j == j != 0
    assert s.p == pytest.approx(j / k, abs=0.0)
    assert s.p_z == pytest.approx(q[1], abs=2e-15)
    assert s.k * s.p == pytest.approx(float(s.j), abs=2e-15)

    phase = s.phase_contract(m=-2)
    assert phase.p == s.p
    assert phase.j == s.j
    assert phase.p_z == s.p_z
    assert phase.x0 == s.x0
    assert phase.epsilon == s.epsilon
    assert phase.m == -2


def test_phase_contract_accepts_source_rational_grid_and_rejects_off_grid():
    c = KokunoOscillatoryPhaseContract(p=0.25, epsilon=0.0625)
    assert c.k == 4
    assert c.j == 1
    assert c.p == 0.25

    # Legacy integer labels remain a subset of the source-compatible grid.
    integer = KokunoOscillatoryPhaseContract(p=2, epsilon=0.2)
    assert integer.k == 3
    assert integer.j == 6

    with pytest.raises(ValueError, match="p=j/k"):
        KokunoOscillatoryPhaseContract(p=0.3, epsilon=0.0625)
    with pytest.raises(ValueError, match="p=j/k"):
        KokunoOscillatoryPhaseContract(p=0.0, epsilon=0.25)


def test_nearest_nonzero_integer_autonomous_tie_convention_is_explicit():
    nearest = KokunoSourcePhaseLabelSelection.nearest_nonzero_integer
    assert nearest(1.49) == 1
    assert nearest(1.50) == 2
    assert nearest(-1.49) == -1
    assert nearest(-1.50) == -2
    assert nearest(0.49) == 1
    assert nearest(-0.49) == -1
    assert nearest(0.0) == 1
    with pytest.raises(ValueError, match="finite"):
        nearest(np.nan)


def test_phase_label_truth_boundary_roundtrip_and_tamper_rejection(tmp_path):
    s = KokunoSourcePhaseLabelSelection()
    path = s.save_json(tmp_path / "phase-label.json")
    loaded = KokunoSourcePhaseLabelSelection.load_json(path)
    assert loaded.sha256 == s.sha256

    payload = json.loads(path.read_text())
    truth = payload["truth_boundary"]
    assert truth["source_phase_label_selection_formulas_executable"] is True
    assert truth["source_requires_p_equals_j_over_k"] is True
    assert truth["autonomous_deterministic_tie_rule_used"] is True
    assert truth["source_actual_discrete_label_recovered"] is False
    assert truth["source_actual_background_path_instantiated"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    payload["truth_boundary"]["paper_exact"] = True
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="truth_boundary metadata changed"):
        KokunoSourcePhaseLabelSelection.load_json(path)


def test_phase_label_input_guards():
    with pytest.raises(ValueError, match="sigma"):
        KokunoSourcePhaseLabelSelection(sigma=0)
    with pytest.raises(ValueError, match="g0"):
        KokunoSourcePhaseLabelSelection(g0=(0.0, 0.0))
    with pytest.raises(ValueError, match="2-vector"):
        KokunoSourcePhaseLabelSelection(K=(1.0, 2.0, 3.0))
    with pytest.raises(ValueError, match="epsilon"):
        KokunoSourcePhaseLabelSelection(epsilon=0.0)
