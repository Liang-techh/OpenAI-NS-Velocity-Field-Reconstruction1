import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_real_pair import (
    KokunoSourceRealConjugatePair,
)


def _mode_data():
    R = np.array([1.2, 1.5, 1.8])
    theta = np.array([0.1, -0.4, 0.7])
    phase = np.array([0.2, 0.8, -0.3])
    n = np.broadcast_to(np.array([1.0, 2.0, 3.0]), (3, 3)).copy()
    # n dot [2,-1,0] = 0, including for a complex scalar amplitude.
    scalar = np.array([0.3 + 0.2j, -0.1 + 0.4j, 0.25 - 0.15j])
    t_plus = scalar[:, None] * np.array([2.0, -1.0, 0.0])[None, :]
    D_r_C_plus = np.array(
        [
            [0.01 + 0.02j, -0.03 + 0.01j, 0.02 - 0.04j],
            [0.02 - 0.01j, 0.01 + 0.03j, -0.04 + 0.02j],
            [-0.01 + 0.04j, 0.02 - 0.02j, 0.03 + 0.01j],
        ]
    )
    D_z_C_plus = np.array(
        [
            [-0.02 + 0.01j, 0.04 + 0.02j, 0.01 - 0.03j],
            [0.03 + 0.01j, -0.02 + 0.04j, 0.02 + 0.01j],
            [0.01 - 0.02j, 0.03 + 0.02j, -0.01 + 0.04j],
        ]
    )
    return R, theta, phase, n, t_plus, D_r_C_plus, D_z_C_plus


def test_real_pair_equals_direct_plus_minus_and_two_real_plus():
    pair = KokunoSourceRealConjugatePair(epsilon=0.25, h=0.005)
    R, _, phase, n, t_plus, D_r_C_plus, D_z_C_plus = _mode_data()

    out = pair.normalized_pair(R, phase, n, t_plus, D_r_C_plus, D_z_C_plus)

    np.testing.assert_allclose(
        out["minus_velocity"], np.conjugate(out["plus_velocity"]), rtol=0.0, atol=2e-13
    )
    np.testing.assert_allclose(
        out["minus_vector_potential"],
        np.conjugate(out["plus_vector_potential"]),
        rtol=0.0,
        atol=2e-13,
    )
    direct = out["plus_velocity"] + out["minus_velocity"]
    np.testing.assert_allclose(out["velocity_cylindrical"], direct.real, rtol=0.0, atol=2e-13)
    np.testing.assert_allclose(
        out["velocity_cylindrical"], 2.0 * out["plus_velocity"].real, rtol=0.0, atol=2e-13
    )
    assert np.max(np.abs(direct.imag)) < 2e-13
    assert np.linalg.norm(out["velocity_cylindrical"]) > 0.0


def test_physical_scaling_and_cartesian_rotation_are_exact():
    pair = KokunoSourceRealConjugatePair(epsilon=0.25, h=0.005)
    R, theta, phase, n, t_plus, D_r_C_plus, D_z_C_plus = _mode_data()
    Q = np.array([0.25, 0.5, 0.75])

    normalized = pair.normalized_pair(R, phase, n, t_plus, D_r_C_plus, D_z_C_plus)
    out = pair.physical_pair(
        Q, R, theta, phase, n, t_plus, D_r_C_plus, D_z_C_plus
    )

    expected_cyl = Q[:, None] ** (-pair.A) * normalized["velocity_cylindrical"]
    np.testing.assert_allclose(out["velocity_physical_cylindrical"], expected_cyl, rtol=2e-15, atol=2e-15)

    c = np.cos(theta)
    s = np.sin(theta)
    expected_cart = np.stack(
        (
            expected_cyl[:, 0] * c - expected_cyl[:, 1] * s,
            expected_cyl[:, 0] * s + expected_cyl[:, 1] * c,
            expected_cyl[:, 2],
        ),
        axis=-1,
    )
    np.testing.assert_allclose(out["velocity_physical_cartesian"], expected_cart, rtol=2e-15, atol=2e-15)
    np.testing.assert_allclose(out["physical_scale"], Q ** (-pair.A), rtol=2e-15, atol=0.0)


def test_batch_broadcast_guards_and_source_h_range():
    pair = KokunoSourceRealConjugatePair(epsilon=0.25, h=0.005)
    R, _, phase, n, t_plus, D_r_C_plus, D_z_C_plus = _mode_data()

    out = pair.physical_pair(
        0.5,
        R,
        0.25,
        phase,
        n,
        t_plus,
        D_r_C_plus,
        D_z_C_plus,
    )
    assert out["velocity_physical_cartesian"].shape == (3, 3)
    assert np.all(np.isfinite(out["velocity_physical_cartesian"]))

    with pytest.raises(ValueError, match="Q must lie"):
        pair.physical_pair(0.0, R, 0.0, phase, n, t_plus, D_r_C_plus, D_z_C_plus)
    with pytest.raises(ValueError, match="Q must lie"):
        pair.physical_pair(1.1, R, 0.0, phase, n, t_plus, D_r_C_plus, D_z_C_plus)
    with pytest.raises(ValueError, match="source range"):
        KokunoSourceRealConjugatePair(epsilon=0.25, h=0.01)
    with pytest.raises(ValueError, match="R>0"):
        pair.normalized_pair(0.0, 0.0, np.array([1.0, 2.0, 3.0]), np.array([2.0, -1.0, 0.0]), np.zeros(3), np.zeros(3))


def test_payload_roundtrip_and_truth_boundary_fail_closed(tmp_path):
    pair = KokunoSourceRealConjugatePair(epsilon=0.125, h=0.004)
    path = pair.save_json(tmp_path / "real_pair.json")
    payload = json.loads(path.read_text())

    restored = KokunoSourceRealConjugatePair.load_json(path)
    assert restored == pair
    assert restored.sha256 == pair.sha256
    assert payload["source"]["commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    assert payload["truth_boundary"]["real_conjugate_pair_materialized_from_caller_mode_data"] is True
    assert payload["truth_boundary"]["source_actual_pulse_forcing_instantiated"] is False
    assert payload["truth_boundary"]["public_xyz_t_velocity_correction_materialized"] is False
    assert payload["truth_boundary"]["paper_exact"] is False
    assert payload["truth_boundary"]["openai_field_identified"] is False

    tampered = json.loads(path.read_text())
    tampered["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary metadata changed"):
        KokunoSourceRealConjugatePair.from_payload(tampered)
