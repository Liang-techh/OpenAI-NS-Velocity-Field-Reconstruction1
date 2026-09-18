import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_homogeneous_spatial_sensitivity import (
    KokunoSourceHomogeneousSpatialSensitivity,
)


def _contract(label="D_z"):
    return KokunoSourceHomogeneousSpatialSensitivity(
        derivative_label=label,
        epsilon=0.25,
        u_star=2.0,
        lambda0=0.7,
        c0=-0.4,
        L_s=1.2,
        sigma=1,
    )


def _geometry(v, xi=0.0):
    v = np.asarray(v, dtype=float)
    n = np.stack(
        (
            0.35 + 0.12 * v + xi * (0.03 + 0.01 * v),
            1.10 + 0.05 * v + xi * (-0.04 + 0.02 * v),
            -0.75 + 0.08 * v + xi * (0.025 - 0.015 * v),
        ),
        axis=-1,
    )
    n_prime = np.stack(
        (
            np.full_like(v, 0.12 + 0.01 * xi),
            np.full_like(v, 0.05 + 0.02 * xi),
            np.full_like(v, 0.08 - 0.015 * xi),
        ),
        axis=-1,
    )
    n_xi = np.stack(
        (0.03 + 0.01 * v, -0.04 + 0.02 * v, 0.025 - 0.015 * v),
        axis=-1,
    )
    n_prime_xi = np.broadcast_to(np.array([0.01, 0.02, -0.015]), n.shape).copy()

    R = 0.8
    F = 0.22 + 0.03 * v + xi * (0.015 - 0.005 * v)
    F_R = -0.04 + 0.01 * v + xi * (0.012 + 0.004 * v)
    G_R = 0.06 - 0.015 * v + xi * (-0.02 + 0.006 * v)
    F_xi = 0.015 - 0.005 * v
    F_R_xi = 0.012 + 0.004 * v
    G_R_xi = -0.02 + 0.006 * v
    K = np.zeros((v.size, 3, 3), dtype=float)
    K_xi = np.zeros_like(K)
    K[:, 0, 1] = -2.0 * F
    K[:, 1, 0] = 2.0 * F + R * F_R
    K[:, 2, 0] = G_R
    K_xi[:, 0, 1] = -2.0 * F_xi
    K_xi[:, 1, 0] = 2.0 * F_xi + R * F_R_xi
    K_xi[:, 2, 0] = G_R_xi
    return n, n_prime, K, n_xi, n_prime_xi, K_xi


def _fourth_difference(solve, step):
    return (
        -solve(2.0 * step)
        + 8.0 * solve(step)
        - 8.0 * solve(-step)
        + solve(-2.0 * step)
    ) / (12.0 * step)


def test_moving_plane_spatial_derivative_matches_independent_finite_difference():
    contract = _contract("D_r")
    v = np.linspace(0.0, contract.L_s, 41)
    n, _, _, n_xi, _, _ = _geometry(v)
    out = contract.moving_plane_sensitivity(v, n, n_xi)
    step = 2.0e-5

    def B_at(xi):
        n_shift, *_ = _geometry(v, xi)
        return contract.pulse.moving_plane(v, n_shift)["B"]

    fd = _fourth_difference(B_at, step)
    np.testing.assert_allclose(out["B_xi"], fd, rtol=2.0e-7, atol=2.0e-9)
    assert float(out["differentiated_transverse_defect_abs"]) < 2.0e-14


def _solve(nodes, label="D_z"):
    contract = _contract(label)
    v = np.linspace(0.0, contract.L_s, nodes)
    geometry = _geometry(v)
    return contract, v, geometry, contract.solve_path(v, *geometry)


def test_nonzero_entrance_t_and_C_sensitivities_match_full_resolve_finite_difference():
    contract, v, geometry, out = _solve(129)
    n, n_prime, K, *_ = geometry
    pulse = contract.pulse
    base = pulse.solve_path(v, n, n_prime, K)
    np.testing.assert_allclose(out["t_plus"], base["t_plus"], rtol=2.0e-13, atol=2.0e-13)
    np.testing.assert_allclose(out["C_plus"], base["C_plus"], rtol=2.0e-13, atol=2.0e-13)

    step = 2.0e-5

    def shifted(key):
        def solve(xi):
            n_s, n_p_s, K_s, *_ = _geometry(v, xi)
            return pulse.solve_path(v, n_s, n_p_s, K_s)[key]
        return solve

    fd_t = _fourth_difference(shifted("t_plus"), step)
    fd_C = _fourth_difference(shifted("C_plus"), step)
    t_error = float(np.max(np.abs(out["t_plus_xi"] - fd_t)))
    C_error = float(np.max(np.abs(out["C_plus_xi"] - fd_C)))
    print({"max_t_plus_xi_fd_error": t_error, "max_C_plus_xi_fd_error": C_error})
    np.testing.assert_allclose(out["t_plus_xi"], fd_t, rtol=3.0e-6, atol=3.0e-8)
    np.testing.assert_allclose(out["C_plus_xi"], fd_C, rtol=4.0e-6, atol=4.0e-8)
    assert out["max_pulse_constraint_defect_abs"] < 3.0e-10
    assert out["max_sensitivity_constraint_defect_abs"] < 3.0e-10
    assert out["entrance_sensitivity_constraint_defect_abs"] < 3.0e-12
    assert np.linalg.norm(out["C_plus_xi"]) > 0.0


def test_spatial_sensitivity_refines_and_keeps_source_coordinate_contract():
    _, _, _, out33 = _solve(33)
    _, _, _, out65 = _solve(65)
    _, _, _, out129 = _solve(129)
    _, _, _, out257 = _solve(257)
    values = [out["C_plus_xi"][-1] for out in (out33, out65, out129, out257)]
    differences = [float(np.linalg.norm(a - b)) for a, b in zip(values[:-1], values[1:])]
    print({"endpoint_C_plus_xi_refinement_differences": differences})
    assert differences[1] < 0.46 * differences[0], differences
    assert differences[2] < 0.46 * differences[1], differences
    for out in (out33, out65, out129, out257):
        assert out["spatial_coordinate_contract"] == "source band rectangle has D_r v=D_z v=0"
        assert out["forcing_origin"].endswith("spatial derivative f_xi=0")
        assert out["max_sensitivity_constraint_defect_abs"] < 3.0e-10


def test_bounds_serialization_and_truth_boundary_fail_closed(tmp_path):
    with pytest.raises(ValueError, match="D_r.*D_z"):
        _contract("partial_R")

    for label in ("D_r", "D_z"):
        contract = _contract(label)
        path = contract.save_json(tmp_path / f"homogeneous_{label}_sensitivity.json")
        payload = json.loads(path.read_text())
        restored = KokunoSourceHomogeneousSpatialSensitivity.load_json(path)
        assert restored == contract
        assert restored.sha256 == contract.sha256
        assert payload["source"]["commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
        truth = payload["truth_boundary"]
        assert truth["source_spatial_pulse_coordinate_invariance_used"] is True
        assert truth["source_homogeneous_C_plus_sensitivity_executable"] is True
        assert truth["source_actual_background_path_instantiated"] is False
        assert truth[f"source_actual_{label}_path_instantiated"] is False
        assert truth["public_xyz_t_velocity_correction_materialized"] is False
        assert truth["pde_validated"] is False
        assert truth["paper_exact"] is False
        assert truth["openai_field_identified"] is False

    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary metadata changed"):
        KokunoSourceHomogeneousSpatialSensitivity.from_payload(payload)
