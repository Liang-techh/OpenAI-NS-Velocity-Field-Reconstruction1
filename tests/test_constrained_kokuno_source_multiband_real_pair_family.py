from __future__ import annotations

import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_band_covering import KokunoSourceBandCovering
from openai_ns_reconstruction.kokuno_source_localized_real_pair_family import (
    KokunoSourceLocalizedRealPairFamily,
)
from openai_ns_reconstruction.kokuno_source_multiband_real_pair_family import (
    KokunoSourceMultiBandRealPairFamily,
)
from openai_ns_reconstruction.kokuno_source_support_localized_curl import (
    KokunoSourceSupportLocalizedCurl,
)


def _inputs(labels=((5, (0, 0, 0)), (6, (1, 0, 0)))):
    R = np.array([0.82, 1.01, 1.19], dtype=float)
    theta = np.array([0.17, -0.31, 0.42], dtype=float)
    alpha = 0.23 + 0.15 * R
    eta = np.stack((np.cos(alpha), np.sin(alpha)), axis=-1)
    tangent = np.stack((-np.sin(alpha), np.cos(alpha)), axis=-1)
    D_r_eta = 0.15 * tangent
    D_z_eta = np.zeros_like(eta)

    phase = np.stack((0.91 * R + 0.13, 1.07 * R - 0.27), axis=-1)
    n_phi = np.zeros((3, 2, 3), dtype=float)
    n_phi[..., 0] = 1.0
    t_plus = np.empty((3, 2, 3), dtype=np.complex128)
    t_plus[:, 0, :] = np.array([0.0, 0.24 + 0.07j, -0.11 + 0.04j])
    t_plus[:, 1, :] = np.array([0.0, -0.18 + 0.05j, 0.16 - 0.03j])
    D_r_C_plus = np.zeros_like(t_plus)
    D_z_C_plus = np.zeros_like(t_plus)
    return {
        "R": R,
        "theta": theta,
        "phase": phase,
        "n_phi": n_phi,
        "t_plus": t_plus,
        "D_r_C_plus": D_r_C_plus,
        "D_z_C_plus": D_z_C_plus,
        "eta": eta,
        "D_r_eta": D_r_eta,
        "D_z_eta": D_z_eta,
        "beta_labels": labels,
    }


def _high_band_inputs():
    """Three high bands make summation-order roundoff macroscopically visible."""
    R = np.array([0.73, 0.91, 1.08, 1.26], dtype=float)
    theta = np.array([0.13, -0.29, 0.47, -0.18], dtype=float)
    alpha = 0.37 + 0.11 * R
    beta = -0.21 + 0.07 * R
    ca, sa = np.cos(alpha), np.sin(alpha)
    cb, sb = np.cos(beta), np.sin(beta)
    eta = np.stack((ca, sa * cb, sa * sb), axis=-1)
    D_r_eta = np.stack(
        (
            -0.11 * sa,
            0.11 * ca * cb - 0.07 * sa * sb,
            0.11 * ca * sb + 0.07 * sa * cb,
        ),
        axis=-1,
    )
    D_z_eta = np.zeros_like(eta)

    phase = np.stack(
        (0.83 * R + 0.19, -1.11 * R + 0.07, 1.37 * R - 0.31), axis=-1
    )
    n_phi = np.zeros((R.size, 3, 3), dtype=float)
    n_phi[..., 0] = 1.0
    t_plus = np.empty((R.size, 3, 3), dtype=np.complex128)
    t_plus[:, 0, :] = np.array([0.0, 0.31 + 0.09j, -0.22 + 0.04j])
    t_plus[:, 1, :] = np.array([0.0, -0.27 + 0.06j, 0.18 - 0.08j])
    t_plus[:, 2, :] = np.array([0.0, 0.14 - 0.11j, 0.29 + 0.03j])
    D_r_C_plus = np.zeros_like(t_plus)
    D_z_C_plus = np.zeros_like(t_plus)
    return {
        "R": R,
        "theta": theta,
        "phase": phase,
        "n_phi": n_phi,
        "t_plus": t_plus,
        "D_r_C_plus": D_r_C_plus,
        "D_z_C_plus": D_z_C_plus,
        "eta": eta,
        "D_r_eta": D_r_eta,
        "D_z_eta": D_z_eta,
        "beta_labels": ((221, (2, 0, 0)), (223, (0, 1, 0)), (225, (0, 0, 2))),
    }


def _permute_beta(data, permutation):
    out = dict(data)
    for key in ("phase", "eta", "D_r_eta", "D_z_eta"):
        out[key] = data[key][..., permutation]
    for key in ("n_phi", "t_plus", "D_r_C_plus", "D_z_C_plus"):
        out[key] = data[key][..., permutation, :]
    out["beta_labels"] = tuple(data["beta_labels"][j] for j in permutation)
    return out


def _evaluate(labels=((5, (0, 0, 0)), (6, (1, 0, 0))), h=0.004):
    contract = KokunoSourceMultiBandRealPairFamily(h=h)
    data = _inputs(labels)
    return contract, data, contract.physical_family(**data)


def test_two_distinct_bands_use_their_own_source_Q_and_epsilon():
    h = 0.004
    _, _, out = _evaluate(h=h)
    expected = [KokunoSourceBandCovering(5, h), KokunoSourceBandCovering(6, h)]
    np.testing.assert_allclose(out["Q_source_by_label"], [b.Q for b in expected], rtol=0, atol=2e-16)
    np.testing.assert_allclose(
        out["epsilon_source_by_label"], [b.epsilon for b in expected], rtol=0, atol=2e-16
    )
    assert out["epsilon_source_by_label"][0] != out["epsilon_source_by_label"][1]
    assert tuple(out["active_ell_bands"]) == (5, 6)
    assert out["distinct_dyadic_band_columns_materialized"] is True
    assert out["genuinely_independent_second_covariance_column_ready"] is False
    assert out["agent3_rank_screen_still_required"] is True


def test_multiband_result_matches_explicit_per_label_complete_curls():
    h = 0.004
    contract, data, out = _evaluate(h=h)
    manual = []
    for j, ell in enumerate((5, 6)):
        band = KokunoSourceBandCovering(ell, h)
        plus = KokunoSourceSupportLocalizedCurl(epsilon=band.epsilon, m=1).localized_mode(
            data["R"], data["phase"][:, j], data["n_phi"][:, j, :],
            data["t_plus"][:, j, :], data["D_r_C_plus"][:, j, :],
            data["D_z_C_plus"][:, j, :], data["eta"][:, j],
            data["D_r_eta"][:, j], data["D_z_eta"][:, j],
        )
        minus = KokunoSourceSupportLocalizedCurl(epsilon=band.epsilon, m=-1).localized_mode(
            data["R"], data["phase"][:, j], data["n_phi"][:, j, :],
            np.conjugate(data["t_plus"][:, j, :]),
            np.conjugate(data["D_r_C_plus"][:, j, :]),
            np.conjugate(data["D_z_C_plus"][:, j, :]), data["eta"][:, j],
            data["D_r_eta"][:, j], data["D_z_eta"][:, j],
        )
        pair = (plus["velocity"] + minus["velocity"]).real
        manual.append((band.Q ** (-contract.A)) * pair)
    manual = np.stack(manual, axis=-2)
    np.testing.assert_allclose(
        out["velocity_physical_cylindrical_by_beta"], manual, rtol=2e-13, atol=2e-13
    )
    np.testing.assert_allclose(
        out["velocity_physical_cylindrical_total"], np.sum(manual, axis=-2), rtol=0, atol=2e-13
    )
    np.testing.assert_allclose(
        np.sum(out["velocity_physical_cartesian_by_band"], axis=-2),
        out["velocity_physical_cartesian_total"], rtol=0, atol=2e-12,
    )


def test_simultaneous_beta_permutation_is_bitwise_invariant_at_high_Q_scaling():
    contract = KokunoSourceMultiBandRealPairFamily(h=0.009)
    data = _high_band_inputs()
    baseline = contract.physical_family(**data)
    permuted = contract.physical_family(**_permute_beta(data, (2, 0, 1)))

    assert tuple(baseline["active_ell_bands"]) == (221, 223, 225)
    assert baseline["deterministic_permutation_invariant_aggregation"] is True
    assert baseline["deterministic_aggregation_beta_order"] == permuted[
        "deterministic_aggregation_beta_order"
    ]
    np.testing.assert_array_equal(
        baseline["velocity_physical_cartesian_by_band"],
        permuted["velocity_physical_cartesian_by_band"],
    )
    np.testing.assert_array_equal(
        baseline["velocity_physical_cartesian_total"],
        permuted["velocity_physical_cartesian_total"],
    )
    np.testing.assert_array_equal(
        baseline["velocity_physical_cylindrical_total"],
        permuted["velocity_physical_cylindrical_total"],
    )


def test_shared_epsilon_legacy_contract_cannot_encode_two_source_band_scales():
    h = 0.004
    eps5 = KokunoSourceBandCovering(5, h).epsilon
    eps6 = KokunoSourceBandCovering(6, h).epsilon
    legacy = KokunoSourceLocalizedRealPairFamily(epsilon=eps5, h=h)
    _, _, out = _evaluate(h=h)
    assert legacy.epsilon == eps5
    assert eps5 != eps6
    np.testing.assert_allclose(out["epsilon_source_by_label"], [eps5, eps6], rtol=0, atol=2e-16)


def test_agent3_handoff_exposes_distinct_band_columns_but_not_rank_claim():
    contract, data, out = _evaluate()
    handoff = contract.agent3_handoff(out, theta=data["theta"], delta_common=0.02, delta_band=-0.01)
    assert tuple(handoff["active_ell_bands"]) == (5, 6)
    assert handoff["rank_screen_columns_cartesian"].shape == (3, 2, 3)
    assert handoff["parameter_count"] == 2
    assert len(handoff["parameter_units"]) == 2
    assert handoff["bounded_coordinates"]["band_contrast_active"] is True
    assert handoff["genuinely_independent_second_covariance_column_ready"] is False
    assert handoff["agent3_rank_screen_still_required"] is True


def test_single_band_phase_copy_and_too_wide_band_window_fail_closed():
    contract = KokunoSourceMultiBandRealPairFamily(h=0.004)
    same_band = _inputs(labels=((5, (0, 0, 0)), (5, (1, 0, 0))))
    with pytest.raises(ValueError, match="same-band phase copy"):
        contract.physical_family(**same_band)

    wide = _inputs(labels=((5, (0, 0, 0)), (10, (1, 0, 0))))
    with pytest.raises(ValueError, match="active interacting band window"):
        contract.physical_family(**wide)


def test_axis_and_zero_band_collapse_fail_closed():
    contract = KokunoSourceMultiBandRealPairFamily(h=0.004)
    axis = _inputs()
    axis["R"] = np.array([0.0, 1.01, 1.19])
    with pytest.raises(ValueError, match="strictly away"):
        contract.physical_family(**axis)

    zero_second = _inputs()
    zero_second["t_plus"] = zero_second["t_plus"].copy()
    zero_second["t_plus"][:, 1, :] = 0.0
    with pytest.raises((ValueError, RuntimeError)):
        contract.physical_family(**zero_second)


def test_truth_boundary_and_serialization_fail_closed(tmp_path):
    contract = KokunoSourceMultiBandRealPairFamily(h=0.004, coefficient_max_l1_update=0.125)
    path = contract.save_json(tmp_path / "multiband.json")
    loaded = KokunoSourceMultiBandRealPairFamily.load_json(path)
    assert loaded.sha256 == contract.sha256

    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["source_per_band_Q_and_epsilon_schedule_executable"] is True
    assert truth["distinct_dyadic_band_columns_materialized_from_supplied_modes"] is True
    assert truth["deterministic_permutation_invariant_aggregation"] is True
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_auxiliary_torus_mode_family_bound"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["genuinely_independent_second_covariance_column_ready"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    payload["truth_boundary"]["paper_exact"] = True
    payload.pop("sha256")
    with pytest.raises(ValueError, match="truth_boundary"):
        KokunoSourceMultiBandRealPairFamily.from_payload(payload)
