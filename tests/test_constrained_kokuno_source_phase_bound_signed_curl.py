import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_band_covering import KokunoSourceBandCovering
from openai_ns_reconstruction.kokuno_source_phase_bound_signed_curl import (
    KokunoSourcePhaseBoundSignedCurlFamily,
)


def _phase_inputs():
    return {
        "R": np.array([0.72, 0.91]),
        "Z": np.array([0.11, -0.07]),
        "theta": 0.43,
        "v": np.array([[0.16, 0.22], [0.19, 0.14]]),
        "beta_labels": ((40, "a"), (41, "b")),
        "R0": np.array([0.75, 0.88]),
        "F0": np.array([1.2, 1.1]),
        "a": np.array([3.1, 3.3]),
        "b_s": np.array([0.28, -0.24]),
        "u_star": np.array([2.0, 2.2]),
        "L_s": np.array([0.62, 0.58]),
        "F": np.array([1.15, 1.04]),
        "G": np.array([0.31, -0.26]),
        "F_R": np.array([0.12, -0.09]),
        "G_R": np.array([-0.08, 0.11]),
        "F_Z": np.array([0.05, -0.04]),
        "G_Z": np.array([0.07, 0.06]),
    }


def _complete_inputs(family, phase_inputs):
    frame = family.phase_frame(**phase_inputs)
    n = frame["n_phi"]
    base = np.stack((n[..., 1], -n[..., 0], np.zeros_like(n[..., 0])), axis=-1)
    t_plus = (1.0 + 0.17j) * base.astype(np.complex128)
    D_r_C = np.zeros_like(t_plus)
    D_z_C = np.zeros_like(t_plus)

    c0 = frame["c0_by_beta"]
    u_star = np.asarray(phase_inputs["u_star"], dtype=float)
    A_c = -c0 * np.sqrt(1.0 + u_star * u_star)
    h_plus = np.array([1.1, 0.95])
    h_minus = np.array([0.9, 1.05])
    p = np.array([0.40, 0.36])
    q = np.array([0.08, -0.06])
    return {
        "t_plus_prototype": t_plus,
        "D_r_C_plus_prototype": D_r_C,
        "D_z_C_plus_prototype": D_z_C,
        "eta": np.array([0.6, 0.8]),
        "D_r_eta": np.array([0.08, -0.06]),
        "D_z_eta": np.array([0.04, -0.03]),
        "A_c": A_c,
        "h_plus": h_plus,
        "h_minus": h_minus,
        "T_N": -A_c * p,
        "T_K": u_star * q,
        "direction_gap_eta": 0.1,
    }


def test_displayed_phase_gradient_matches_independent_fd4():
    family = KokunoSourcePhaseBoundSignedCurlFamily(h=0.005)
    args = _phase_inputs()
    frame = family.phase_frame(**args)

    assert frame["sigma_labels"] == ("sigma_plus", "sigma_minus")
    assert np.all(frame["kp_by_beta_sign"] != 0)
    for j, ell in enumerate((40, 41)):
        schedule = KokunoSourceBandCovering(ell, family.h)
        k = frame["k_by_beta"][j]
        assert k == int(np.ceil(schedule.epsilon ** -0.5))
        assert np.all(frame["rounding_error_by_beta_sign"][j] <= 1.0 / k + 1e-15)
        assert 1.0 <= schedule.epsilon * k * k <= 4.0

        for sigma_index in range(2):
            p = frame["p_by_beta_sign"][j, sigma_index]
            pz = frame["p_z_by_beta_sign"][j, sigma_index]
            x0 = frame["x0_by_beta_sign"][j, sigma_index]
            vv = args["v"][j, sigma_index]
            R0 = args["R"][j]
            Z0 = args["Z"][j]
            th0 = args["theta"]
            F = args["F"][j]
            G = args["G"][j]
            FR = args["F_R"][j]
            GR = args["G_R"][j]
            FZ = args["F_Z"][j]
            GZ = args["G_Z"][j]

            def phi(R, th, Z):
                FF = F + FR * (R - R0) + FZ * (Z - Z0)
                GG = G + GR * (R - R0) + GZ * (Z - Z0)
                return p * th + pz * Z / schedule.epsilon + x0 * R - vv * (p * FF + pz * GG)

            step = 2.0e-5

            def fd4(which):
                def eval_at(offset):
                    R, th, Z = R0, th0, Z0
                    if which == "R":
                        R += offset
                    elif which == "theta":
                        th += offset
                    else:
                        Z += offset
                    return phi(R, th, Z)
                return (
                    -eval_at(2 * step)
                    + 8 * eval_at(step)
                    - 8 * eval_at(-step)
                    + eval_at(-2 * step)
                ) / (12 * step)

            expected = np.array(
                [
                    fd4("R"),
                    fd4("theta") / R0,
                    schedule.epsilon * fd4("Z"),
                ]
            )
            np.testing.assert_allclose(
                frame["n_phi"][j, sigma_index], expected, rtol=3e-10, atol=3e-10
            )


def test_phase_bound_bridge_equals_direct_signed_complete_curl_call():
    family = KokunoSourcePhaseBoundSignedCurlFamily(h=0.005)
    phase_args = _phase_inputs()
    extra = _complete_inputs(family, phase_args)

    bridged = family.physical_family(**phase_args, **extra)
    frame = family.phase_frame(**phase_args)
    direct = family.complete_curl_family.physical_family(
        R=phase_args["R"],
        theta=phase_args["theta"],
        phase=frame["phase"],
        n_phi=frame["n_phi"],
        t_plus_prototype=extra["t_plus_prototype"],
        D_r_C_plus_prototype=extra["D_r_C_plus_prototype"],
        D_z_C_plus_prototype=extra["D_z_C_plus_prototype"],
        eta=extra["eta"],
        D_r_eta=extra["D_r_eta"],
        D_z_eta=extra["D_z_eta"],
        beta_labels=phase_args["beta_labels"],
        A_c=extra["A_c"],
        u_star=phase_args["u_star"],
        h_plus=extra["h_plus"],
        h_minus=extra["h_minus"],
        T_N=extra["T_N"],
        T_K=extra["T_K"],
        direction_gap_eta=extra["direction_gap_eta"],
    )
    np.testing.assert_array_equal(
        bridged["velocity_physical_cartesian_total"],
        direct["velocity_physical_cartesian_total"],
    )
    np.testing.assert_array_equal(
        bridged["velocity_physical_cartesian_by_beta_sign"],
        direct["velocity_physical_cartesian_by_beta_sign"],
    )
    assert bridged["source_phase_formula_bound"]
    assert bridged["source_phase_gradient_formula_bound"]
    assert not bridged["actual_positive_order_background_bound"]
    assert not bridged["public_xyz_t_velocity_correction_materialized"]


def test_phase_formula_retains_source_sign_and_nonfree_carrier():
    family = KokunoSourcePhaseBoundSignedCurlFamily()
    args = _phase_inputs()
    frame = family.phase_frame(**args)
    np.testing.assert_allclose(frame["x0_by_beta_sign"][:, 0], -frame["x0_by_beta_sign"][:, 1])
    reconstructed = frame["kp_by_beta_sign"] / frame["k_by_beta"][:, None]
    np.testing.assert_array_equal(frame["p_by_beta_sign"], reconstructed)
    assert np.all(frame["B_s_by_beta"] > 0.0)
    assert np.all(frame["c0_by_beta"] < 0.0)
    assert np.all(frame["v_s_by_beta"] > 2.0)


def test_phase_frame_fails_closed_on_axis_cone_and_rectangle_domain():
    family = KokunoSourcePhaseBoundSignedCurlFamily()
    args = _phase_inputs()

    bad_axis = dict(args)
    bad_axis["R"] = np.array([0.0, 0.91])
    with pytest.raises(ValueError, match="R>0"):
        family.phase_frame(**bad_axis)

    bad_cone = dict(args)
    bad_cone["a"] = np.array([1.5, 3.3])
    bad_cone["b_s"] = np.array([0.0, -0.24])
    with pytest.raises(ValueError, match="v_s"):
        family.phase_frame(**bad_cone)

    bad_v = dict(args)
    vv = np.array(args["v"], copy=True)
    vv[0, 0] = args["L_s"][0] + 0.01
    bad_v["v"] = vv
    with pytest.raises(ValueError, match="0<=v<=L_s"):
        family.phase_frame(**bad_v)


def test_source_formula_vs_caller_data_truth_boundary_and_serialization(tmp_path):
    family = KokunoSourcePhaseBoundSignedCurlFamily()
    payload = family.to_payload()
    assert payload["truth_boundary"]["source_phase_formula_bound"]
    assert payload["truth_boundary"]["source_phase_gradient_formula_bound"]
    assert "caller supplied" in payload["caller_contract"]["background"]
    for key in (
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_bound",
        "actual_auxiliary_torus_mode_family_bound",
        "source_actual_partition_labels_instantiated",
        "public_xyz_t_velocity_correction_materialized",
        "actual_source_physical_covariance_rank_two_assessed",
        "genuinely_independent_second_covariance_column_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated",
        "paper_exact",
        "openai_field_identified",
        "blowup_proved",
    ):
        assert payload["truth_boundary"][key] is False

    path = family.save_json(tmp_path / "phase_bound.json")
    loaded = KokunoSourcePhaseBoundSignedCurlFamily.load_json(path)
    assert loaded.sha256 == family.sha256
    mutated = json.loads(path.read_text())
    mutated["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary metadata changed"):
        KokunoSourcePhaseBoundSignedCurlFamily.from_payload(mutated)
