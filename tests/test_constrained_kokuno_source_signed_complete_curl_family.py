import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_band_covering import KokunoSourceBandCovering
from openai_ns_reconstruction.kokuno_source_signed_complete_curl_family import (
    KokunoSourceSignedCompleteCurlFamily,
)
from openai_ns_reconstruction.kokuno_source_support_localized_curl import (
    KokunoSourceSupportLocalizedCurl,
)


def _fixture():
    beta_labels = ((5, "a"), (6, "b"))
    eta = np.array([0.6, 0.8])
    D_r_eta = np.array([0.08, -0.06])
    D_z_eta = np.array([0.04, -0.03])
    phase = np.array([[0.17, -0.23], [0.31, 0.07]])
    n_phi = np.zeros((2, 2, 3))
    n_phi[...] = np.array([1.0, 0.0, 1.0])
    t_plus = np.empty((2, 2, 3), dtype=np.complex128)
    t_plus[:, 0, :] = np.array([1.0 + 0.2j, 0.0 + 0.3j, -1.0 - 0.2j])
    t_plus[:, 1, :] = np.array([0.0, 1.0 + 0.1j, 0.0])
    t_plus[1] *= 0.83
    D_r_C = np.zeros_like(t_plus)
    D_z_C = np.zeros_like(t_plus)
    A_c = np.array([2.0, 2.2])
    u_star = np.array([3.0, 2.8])
    h_plus = np.array([1.1, 0.95])
    h_minus = np.array([0.9, 1.05])
    p = np.array([0.40, 0.36])
    q = np.array([0.08, -0.06])
    T_N = -A_c * p
    T_K = u_star * q
    return {
        "R": np.array([0.72, 0.91]),
        "theta": 0.43,
        "phase": phase,
        "n_phi": n_phi,
        "t_plus_prototype": t_plus,
        "D_r_C_plus_prototype": D_r_C,
        "D_z_C_plus_prototype": D_z_C,
        "eta": eta,
        "D_r_eta": D_r_eta,
        "D_z_eta": D_z_eta,
        "beta_labels": beta_labels,
        "A_c": A_c,
        "u_star": u_star,
        "h_plus": h_plus,
        "h_minus": h_minus,
        "T_N": T_N,
        "T_K": T_K,
        "D_r_A_c": np.array([0.03, -0.02]),
        "D_r_u_star": np.array([0.04, 0.01]),
        "D_r_h_plus": np.array([0.05, -0.03]),
        "D_r_h_minus": np.array([-0.02, 0.04]),
        "D_r_T_N": np.array([0.12, -0.09]),
        "D_r_T_K": np.array([0.08, 0.06]),
        "D_z_A_c": np.array([-0.02, 0.01]),
        "D_z_u_star": np.array([0.03, -0.04]),
        "D_z_h_plus": np.array([-0.04, 0.02]),
        "D_z_h_minus": np.array([0.03, -0.01]),
        "D_z_T_N": np.array([-0.10, 0.07]),
        "D_z_T_K": np.array([0.06, -0.05]),
        "direction_gap_eta": 0.1,
    }


def test_signed_amplitude_gradients_enter_complete_curl_exactly():
    family = KokunoSourceSignedCompleteCurlFamily(h=0.005)
    args = _fixture()
    result = family.physical_family(**args)
    j = 0
    sigma = 0
    schedule = KokunoSourceBandCovering(args["beta_labels"][j][0], family.h)
    localizer = KokunoSourceSupportLocalizedCurl(
        epsilon=schedule.epsilon, m=1, partition_atol=family.partition_atol
    )
    n = args["n_phi"][j, sigma]
    t_proto = args["t_plus_prototype"][j, sigma]
    C_proto = localizer.complete_curl.coefficient(n, t_proto)
    a = result["reference_amplitudes"][j, sigma]
    D_r_a = result["D_r_reference_amplitudes"][j, sigma]
    D_z_a = result["D_z_reference_amplitudes"][j, sigma]
    g = np.sqrt(schedule.epsilon) * a
    D_r_g = np.sqrt(schedule.epsilon) * D_r_a
    D_z_g = np.sqrt(schedule.epsilon) * D_z_a
    t_scaled = g * t_proto
    D_r_C_scaled = g * args["D_r_C_plus_prototype"][j, sigma] + D_r_g * C_proto
    D_z_C_scaled = g * args["D_z_C_plus_prototype"][j, sigma] + D_z_g * C_proto
    plus = localizer.localized_mode(
        args["R"][j], args["phase"][j, sigma], n, t_scaled, D_r_C_scaled,
        D_z_C_scaled, args["eta"][j], args["D_r_eta"][j], args["D_z_eta"][j]
    )
    minus = KokunoSourceSupportLocalizedCurl(
        epsilon=schedule.epsilon, m=-1, partition_atol=family.partition_atol
    ).localized_mode(
        args["R"][j], args["phase"][j, sigma], n, np.conjugate(t_scaled),
        np.conjugate(D_r_C_scaled), np.conjugate(D_z_C_scaled), args["eta"][j],
        args["D_r_eta"][j], args["D_z_eta"][j]
    )
    expected = schedule.Q ** (-family.A) * (plus["velocity"] + minus["velocity"]).real
    np.testing.assert_allclose(
        result["velocity_physical_cylindrical_by_beta_sign"][j, sigma],
        expected, rtol=2e-13, atol=2e-13
    )

    dropped_plus = localizer.localized_mode(
        args["R"][j], args["phase"][j, sigma], n, t_scaled,
        g * args["D_r_C_plus_prototype"][j, sigma],
        g * args["D_z_C_plus_prototype"][j, sigma], args["eta"][j],
        args["D_r_eta"][j], args["D_z_eta"][j]
    )
    dropped_minus = KokunoSourceSupportLocalizedCurl(
        epsilon=schedule.epsilon, m=-1, partition_atol=family.partition_atol
    ).localized_mode(
        args["R"][j], args["phase"][j, sigma], n, np.conjugate(t_scaled),
        np.conjugate(g * args["D_r_C_plus_prototype"][j, sigma]),
        np.conjugate(g * args["D_z_C_plus_prototype"][j, sigma]), args["eta"][j],
        args["D_r_eta"][j], args["D_z_eta"][j]
    )
    dropped = schedule.Q ** (-family.A) * (
        dropped_plus["velocity"] + dropped_minus["velocity"]
    ).real
    relative = np.linalg.norm(expected - dropped) / max(np.linalg.norm(expected), 1e-15)
    assert relative > 1e-4


def test_sign_beta_total_contract_and_source_scales():
    family = KokunoSourceSignedCompleteCurlFamily(h=0.005)
    result = family.physical_family(**_fixture())
    assert result["velocity_physical_cartesian_by_beta_sign"].shape == (2, 2, 3)
    np.testing.assert_allclose(
        result["velocity_physical_cartesian_by_beta"],
        np.sum(result["velocity_physical_cartesian_by_beta_sign"], axis=-2),
        rtol=0.0, atol=2e-12
    )
    expected_total = np.sum(np.take(
        result["velocity_physical_cartesian_by_beta"], result["canonical_beta_indices"], axis=-2
    ), axis=-2)
    np.testing.assert_array_equal(result["velocity_physical_cartesian_total"], expected_total)
    for j, ell in enumerate((5, 6)):
        schedule = KokunoSourceBandCovering(ell, family.h)
        assert result["Q_by_beta"][j] == schedule.Q
        assert result["epsilon_by_beta"][j] == schedule.epsilon
    assert result["reference_covariance_rank_two"]
    assert not result["physical_complete_curl_covariance_rank_two_assessed"]
    assert not result["genuinely_independent_second_covariance_column_ready"]


def test_semantic_sigma_swap_with_tk_sign_flip_preserves_total():
    family = KokunoSourceSignedCompleteCurlFamily(h=0.005)
    args = _fixture()
    baseline = family.physical_family(**args)["velocity_physical_cartesian_total"]
    swapped = dict(args)
    for key in ("phase", "n_phi", "t_plus_prototype", "D_r_C_plus_prototype",
                "D_z_C_plus_prototype"):
        swapped[key] = np.asarray(args[key])[:, ::-1, ...]
    swapped["h_plus"] = args["h_minus"]
    swapped["h_minus"] = args["h_plus"]
    swapped["D_r_h_plus"] = args["D_r_h_minus"]
    swapped["D_r_h_minus"] = args["D_r_h_plus"]
    swapped["D_z_h_plus"] = args["D_z_h_minus"]
    swapped["D_z_h_minus"] = args["D_z_h_plus"]
    swapped["T_K"] = -args["T_K"]
    swapped["D_r_T_K"] = -args["D_r_T_K"]
    swapped["D_z_T_K"] = -args["D_z_T_K"]
    candidate = family.physical_family(**swapped)["velocity_physical_cartesian_total"]
    np.testing.assert_allclose(candidate, baseline, rtol=3e-13, atol=3e-13)


def test_axis_partition_and_provenance_fail_closed(tmp_path):
    family = KokunoSourceSignedCompleteCurlFamily()
    args = _fixture()
    bad_axis = dict(args)
    bad_axis["R"] = np.array([0.0, 0.91])
    with pytest.raises(ValueError, match="R>0"):
        family.physical_family(**bad_axis)
    bad_partition = dict(args)
    bad_partition["eta"] = np.array([0.61, 0.8])
    with pytest.raises(ValueError, match="squared partition"):
        family.physical_family(**bad_partition)
    path = family.save_json(tmp_path / "signed_complete_curl.json")
    loaded = KokunoSourceSignedCompleteCurlFamily.load_json(path)
    assert loaded.sha256 == family.sha256
    payload = json.loads(path.read_text())
    payload["truth_boundary"]["paper_exact"] = True
    with pytest.raises(ValueError, match="truth_boundary metadata changed"):
        KokunoSourceSignedCompleteCurlFamily.from_payload(payload)


def test_truth_boundary_does_not_promote_missing_source_data():
    truth = KokunoSourceSignedCompleteCurlFamily().to_payload()["truth_boundary"]
    assert truth["signed_amplitude_gradient_complete_curl_terms_retained"]
    assert truth["q_scaled_sign_by_beta_total_physical_family_executable"]
    for key in (
        "actual_positive_order_background_bound",
        "actual_source_h_sigma_pulse_integrals_bound",
        "actual_signed_auxiliary_rectangles_bound",
        "actual_auxiliary_torus_mode_family_bound",
        "source_actual_partition_labels_instantiated",
        "public_xyz_t_velocity_correction_materialized",
        "genuinely_independent_second_covariance_column_ready",
        "formal_full_domain_pde_gate_assessed",
        "pde_validated", "paper_exact", "openai_field_identified", "blowup_proved",
    ):
        assert truth[key] is False
