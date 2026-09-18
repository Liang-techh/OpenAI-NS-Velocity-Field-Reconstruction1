import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_homogeneous_pulse import (
    KokunoSourceHomogeneousPrimaryPulse,
)


def _pulse():
    return KokunoSourceHomogeneousPrimaryPulse(
        epsilon=0.25,
        u_star=2.0,
        lambda0=0.7,
        c0=-0.4,
        L_s=1.2,
        sigma=1,
    )


def _geometry(v):
    v = np.asarray(v, dtype=float)
    n = np.stack(
        (
            0.35 + 0.12 * v,
            1.10 + 0.05 * v,
            -0.75 + 0.08 * v,
        ),
        axis=-1,
    )
    n_prime = np.broadcast_to(np.array([0.12, 0.05, 0.08]), n.shape).copy()
    R = 0.8
    F = 0.22 + 0.03 * v
    F_R = -0.04 + 0.01 * v
    G_R = 0.06 - 0.015 * v
    K = np.zeros((v.size, 3, 3), dtype=float)
    K[:, 0, 1] = -2.0 * F
    K[:, 1, 0] = 2.0 * F + R * F_R
    K[:, 2, 0] = G_R
    return n, n_prime, K


def test_source_envelope_uses_exact_antiderivative_and_peaks_at_midpoint():
    pulse = _pulse()
    v = np.array([0.0, 0.21, pulse.L_s / 2.0, 0.87, pulse.L_s])
    P = pulse.envelope(v)

    assert np.all(P > 0.0)
    assert np.max(P) <= 1.0 + 2.0e-15
    np.testing.assert_allclose(P[2], 1.0, rtol=0.0, atol=2.0e-15)

    point = 0.41
    step = 2.0e-6
    derivative = (
        pulse.log_envelope(point + step) - pulse.log_envelope(point - step)
    ) / (2.0 * step)
    np.testing.assert_allclose(derivative, pulse.a_net(point), rtol=2.0e-9, atol=2.0e-10)


def test_moving_plane_is_transverse_and_coordinates_invert_B():
    pulse = _pulse()
    v = np.linspace(0.0, pulse.L_s, 17)
    n, _, _ = _geometry(v)
    plane = pulse.moving_plane(v, n)

    defect = np.einsum("ni,nij->nj", n, plane["B"])
    assert np.max(np.abs(defect)) < 2.0e-15

    z = np.stack(
        (
            pulse.envelope(v) * (1.0 + 0.05 * v),
            0.08 * pulse.envelope(v) * np.sin(1.3 * v),
        ),
        axis=-1,
    )
    t = np.einsum("nij,nj->ni", plane["B"], z)
    recovered = pulse.plane_coordinates(v, n, t)
    np.testing.assert_allclose(recovered["z"].real, z, rtol=2.0e-14, atol=2.0e-14)
    assert np.max(np.abs(recovered["z"].imag)) == 0.0


def _solve(nodes):
    pulse = _pulse()
    v = np.linspace(0.0, pulse.L_s, nodes)
    n, n_prime, K = _geometry(v)
    return pulse, pulse.solve_path(v, n, n_prime, K)


def test_homogeneous_source_path_refines_and_feeds_positive_mode_curl():
    pulse33, out33 = _solve(33)
    _, out65 = _solve(65)
    _, out129 = _solve(129)
    _, out257 = _solve(257)

    t33 = out33["t_real_cosine"][-1]
    t65 = out65["t_real_cosine"][-1]
    t129 = out129["t_real_cosine"][-1]
    t257 = out257["t_real_cosine"][-1]
    differences = [
        float(np.linalg.norm(t33 - t65)),
        float(np.linalg.norm(t65 - t129)),
        float(np.linalg.norm(t129 - t257)),
    ]
    assert differences[1] < 0.45 * differences[0], differences
    assert differences[2] < 0.45 * differences[1], differences

    for out in (out33, out65, out129, out257):
        assert out["max_constraint_defect_abs"] < 2.0e-10
        assert out["initial_coordinate_error_abs"] < 2.0e-12
        np.testing.assert_allclose(
            out["t_plus"], 0.5 * out["t_real_cosine"], rtol=0.0, atol=0.0
        )
        assert np.max(np.abs(out["t_plus"].imag)) == 0.0
        assert np.linalg.norm(out["C_plus"][1:]) > 0.0
        assert out["forcing_origin"] == "source homogeneous primary pulse: f=0"

    np.testing.assert_allclose(
        out33["source_initial_z"],
        np.array([pulse33.envelope(0.0), 0.0]),
        rtol=0.0,
        atol=2.0e-15,
    )


def test_guards_serialization_and_truth_boundary_fail_closed(tmp_path):
    pulse = _pulse()
    v = np.linspace(0.0, pulse.L_s, 9)
    n, n_prime, K = _geometry(v)

    with pytest.raises(ValueError, match="v=0 to v=L_s"):
        pulse.solve_path(v + 0.02, n, n_prime, K)

    bad_n = n.copy()
    bad_n[3, 1:] = 0.0
    with pytest.raises(ValueError, match="tangential norm"):
        pulse.solve_path(v, bad_n, n_prime, K)

    with pytest.raises(ValueError, match="strictly negative"):
        KokunoSourceHomogeneousPrimaryPulse(c0=0.1)
    with pytest.raises(ValueError, match="sigma"):
        KokunoSourceHomogeneousPrimaryPulse(sigma=0)

    path = pulse.save_json(tmp_path / "homogeneous_pulse.json")
    payload = json.loads(path.read_text())
    restored = KokunoSourceHomogeneousPrimaryPulse.load_json(path)
    assert restored == pulse
    assert restored.sha256 == pulse.sha256
    assert payload["source"]["commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    truth = payload["truth_boundary"]
    assert truth["source_real_homogeneous_primary_pulse_equations_executable"] is True
    assert truth["source_complete_curl_C_plus_executable"] is True
    assert truth["source_actual_background_path_instantiated"] is False
    assert truth["source_homogeneous_D_r_C_plus_instantiated"] is False
    assert truth["source_homogeneous_D_z_C_plus_instantiated"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    assert truth["openai_field_identified"] is False

    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary metadata changed"):
        KokunoSourceHomogeneousPrimaryPulse.from_payload(payload)
