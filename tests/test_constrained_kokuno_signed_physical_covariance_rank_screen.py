import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_signed_physical_covariance_rank_screen import (
    KokunoSignedPhysicalCovarianceRankScreen,
)
from openai_ns_reconstruction.kokuno_source_signed_complete_curl_family import (
    KokunoSourceSignedCompleteCurlFamily,
)


def _synthetic_result(*, duplicate_signs: bool = False):
    sample_count = 8
    plus = np.array([1.0, 1.0, 0.0])
    minus = plus.copy() if duplicate_signs else np.array([1.0, 0.0, 1.0])
    by = np.zeros((sample_count, 2, 2, 3), dtype=float)
    by[:, 0, 0, :] = 0.6 * plus
    by[:, 1, 0, :] = 0.4 * plus
    by[:, 0, 1, :] = 0.35 * minus
    by[:, 1, 1, :] = 0.65 * minus
    total = np.sum(by, axis=(-3, -2))
    return {
        "beta_labels": ((5, "a"), (6, "b")),
        "sigma_labels": ("sigma_plus", "sigma_minus"),
        "velocity_physical_cylindrical_by_beta_sign": by,
        "velocity_physical_cylindrical_total": total,
        "reference_covariance_rank_two": True,
        "physical_complete_curl_covariance_rank_two_assessed": False,
        "genuinely_independent_second_covariance_column_ready": False,
    }


def test_exact_product_rule_rank_two_and_cross_interaction_identity():
    result = KokunoSignedPhysicalCovarianceRankScreen().evaluate(
        _synthetic_result(), averaging_axes=(0,)
    )
    np.testing.assert_allclose(
        result["covariance_response_jacobian"],
        np.array([[3.0, 1.0], [1.0, 3.0]]),
        rtol=0.0,
        atol=2e-15,
    )
    np.testing.assert_allclose(result["singular_values"], np.array([4.0, 2.0]))
    assert result["rank_two_cells"] == 1
    assert result["total_cells"] == 1
    assert result["local_supplied_signed_physical_covariance_rank_two"]
    assert result["minimum_singular_value_ratio"] == pytest.approx(0.5)
    assert result["maximum_quadratic_homogeneity_relative_error"] < 1e-15
    assert result["cross_sign_and_cross_beta_terms_retained"] is True
    assert result["genuinely_independent_second_covariance_column_ready"] is False
    assert result["finite_correction_cycle_rerun_allowed"] is False


def test_duplicate_physical_sign_column_cannot_manufacture_rank_two():
    result = KokunoSignedPhysicalCovarianceRankScreen().evaluate(
        _synthetic_result(duplicate_signs=True), averaging_axes=(0,)
    )
    assert result["rank_two_cells"] == 0
    assert result["rank_two_fraction"] == 0.0
    assert not result["local_supplied_signed_physical_covariance_rank_two"]


def test_total_or_sigma_semantic_drift_fails_closed():
    screen = KokunoSignedPhysicalCovarianceRankScreen()
    bad_total = _synthetic_result()
    bad_total["velocity_physical_cylindrical_total"] = (
        bad_total["velocity_physical_cylindrical_total"].copy()
    )
    bad_total["velocity_physical_cylindrical_total"][0, 0] += 1e-4
    with pytest.raises(RuntimeError, match="do not reproduce"):
        screen.evaluate(bad_total, averaging_axes=(0,))

    bad_sigma = _synthetic_result()
    bad_sigma["sigma_labels"] = ("+", "-")
    with pytest.raises(ValueError, match="source sigma ordering"):
        screen.evaluate(bad_sigma, averaging_axes=(0,))


def _supplied_signed_complete_curl_family():
    radial = np.array([0.78, 0.90, 1.02], dtype=float)
    quadrature_phase = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False, dtype=float)
    R, theta = np.meshgrid(radial, quadrature_phase, indexing="ij")

    alpha = 0.55 + 0.04 * (R - 0.90)
    eta = np.stack((np.cos(alpha), np.sin(alpha)), axis=-1)
    tangent = np.stack((-np.sin(alpha), np.cos(alpha)), axis=-1)
    D_r_eta = 0.04 * tangent
    D_z_eta = np.zeros_like(eta)

    phase = np.empty(R.shape + (2, 2), dtype=float)
    phase[..., 0, 0] = quadrature_phase[None, :] + 0.17 + 0.20 * R
    phase[..., 0, 1] = quadrature_phase[None, :] - 0.23 + 0.20 * R
    phase[..., 1, 0] = quadrature_phase[None, :] + 0.31 - 0.10 * R
    phase[..., 1, 1] = quadrature_phase[None, :] + 0.07 - 0.10 * R

    n_phi = np.zeros(R.shape + (2, 2, 3), dtype=float)
    n_phi[...] = np.array([1.0, 0.0, 1.0])
    t_plus = np.empty(R.shape + (2, 2, 3), dtype=np.complex128)
    prototypes = np.array(
        [
            [[1.0 + 0.2j, 0.0 + 0.3j, -1.0 - 0.2j],
             [0.0 + 0.0j, 1.0 + 0.1j, 0.0 + 0.0j]],
            [[0.83 + 0.166j, 0.0 + 0.249j, -0.83 - 0.166j],
             [0.0 + 0.0j, 0.83 + 0.083j, 0.0 + 0.0j]],
        ],
        dtype=np.complex128,
    )
    t_plus[...] = prototypes
    D_r_C = np.zeros_like(t_plus)
    D_z_C = np.zeros_like(t_plus)

    A_c = np.array([2.0, 2.2])
    u_star = np.array([3.0, 2.8])
    h_plus = np.array([1.1, 0.95])
    h_minus = np.array([0.9, 1.05])
    p = np.array([0.40, 0.36])
    q = np.array([0.08, -0.06])

    family = KokunoSourceSignedCompleteCurlFamily(h=0.005)
    return family.physical_family(
        R=R,
        theta=theta,
        phase=phase,
        n_phi=n_phi,
        t_plus_prototype=t_plus,
        D_r_C_plus_prototype=D_r_C,
        D_z_C_plus_prototype=D_z_C,
        eta=eta,
        D_r_eta=D_r_eta,
        D_z_eta=D_z_eta,
        beta_labels=((5, "a"), (6, "b")),
        A_c=A_c,
        u_star=u_star,
        h_plus=h_plus,
        h_minus=h_minus,
        T_N=-A_c * p,
        T_K=u_star * q,
        D_r_A_c=np.array([0.03, -0.02]),
        D_r_u_star=np.array([0.04, 0.01]),
        D_r_h_plus=np.array([0.05, -0.03]),
        D_r_h_minus=np.array([-0.02, 0.04]),
        D_r_T_N=np.array([0.12, -0.09]),
        D_r_T_K=np.array([0.08, 0.06]),
        D_z_A_c=np.array([-0.02, 0.01]),
        D_z_u_star=np.array([0.03, -0.04]),
        D_z_h_plus=np.array([-0.04, 0.02]),
        D_z_h_minus=np.array([0.03, -0.01]),
        D_z_T_N=np.array([-0.10, 0.07]),
        D_z_T_K=np.array([0.06, -0.05]),
        direction_gap_eta=0.1,
    )


def test_agent2_signed_complete_curl_handoff_is_screenable_without_promotion():
    physical = _supplied_signed_complete_curl_family()
    result = KokunoSignedPhysicalCovarianceRankScreen().evaluate(
        physical, averaging_axes=(1,)
    )
    assert result["total_cells"] == 3
    assert np.all(np.isfinite(result["singular_values"]))
    assert result["maximum_quadratic_homogeneity_relative_error"] <= 2e-12
    assert result["physical_complete_curl_covariance_rank_two_assessed"] is True
    assert result["actual_source_h_sigma_pulse_integrals_bound"] is False
    assert result["actual_signed_auxiliary_rectangles_bound"] is False
    assert result["real_candidate_defect_consumed"] is False
    assert result["genuinely_independent_second_covariance_column_ready"] is False
    assert result["heldout_ns_residual_assessed"] is False
    assert result["pde_validated"] is False
