import copy

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_source_signed_covariance_pair import (
    KokunoSourceSignedCovariancePair,
)


def _args():
    return dict(A_c=2.0, u_star=3.0, h_plus=1.0, h_minus=1.0, T_N=-2.0, T_K=0.6)


def test_reference_two_sign_inverse_has_nonzero_covariance_rank():
    pair = KokunoSourceSignedCovariancePair()
    out = pair.solve_reference_amplitudes(**_args(), direction_gap_eta=0.2)

    np.testing.assert_allclose(
        out["reference_column_matrix"],
        np.array([[-2.0, -2.0], [-3.0, 3.0]]),
        rtol=0.0,
        atol=0.0,
    )
    assert out["reference_column_determinant"] == pytest.approx(-12.0)
    np.testing.assert_allclose(out["squared_amplitudes"], [0.4, 0.6], rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(out["reconstructed_target"], [-2.0, 0.6], rtol=0.0, atol=2e-15)
    assert out["reference_covariance_rank_two"] is True
    assert out["positive_reference_coefficients"] is True
    assert out["sigma_axis_is_not_fourier_harmonic_axis"] is True
    assert out["genuinely_independent_second_covariance_column_ready"] is False


def test_vectorized_reference_inverse_and_exact_determinant_identity():
    pair = KokunoSourceSignedCovariancePair()
    A = np.array([1.5, 2.0, 2.5])
    u = np.array([2.0, 3.0, 4.0])
    hp = np.array([0.7, 1.1, 1.4])
    hm = np.array([0.9, 1.3, 1.8])
    TN = -A * np.array([1.0, 1.2, 1.4])
    TK = u * np.array([0.1, -0.2, 0.3])
    out = pair.solve_reference_amplitudes(A, u, hp, hm, TN, TK)

    expected_det = -2.0 * A * u * hp * hm
    np.testing.assert_allclose(out["reference_column_determinant"], expected_det, rtol=3e-15, atol=0.0)
    np.testing.assert_allclose(out["reconstructed_target"], np.stack((TN, TK), axis=-1), rtol=0.0, atol=3e-15)
    assert np.all(out["squared_amplitudes"] > 0.0)


def test_source_direction_gap_and_positive_cone_fail_closed():
    pair = KokunoSourceSignedCovariancePair()

    with pytest.raises(ValueError, match="positive two-sign source covariance cone"):
        pair.solve_reference_amplitudes(A_c=2.0, u_star=3.0, h_plus=1.0, h_minus=1.0, T_N=-2.0, T_K=6.0)

    # p=1, q=.8 remains in the positive cone, but violates the declared eta=.3 gap |q|<=.7 p.
    with pytest.raises(ValueError, match="strict source direction gap"):
        pair.solve_reference_amplitudes(A_c=2.0, u_star=3.0, h_plus=1.0, h_minus=1.0, T_N=-2.0, T_K=2.4, direction_gap_eta=0.3)

    with pytest.raises(ValueError, match="A_c must be strictly positive"):
        pair.solve_reference_amplitudes(A_c=0.0, u_star=3.0, h_plus=1.0, h_minus=1.0, T_N=-2.0, T_K=0.0)


def test_directional_derivative_matches_independent_centered_difference():
    pair = KokunoSourceSignedCovariancePair()
    args = dict(A_c=2.2, u_star=3.1, h_plus=0.8, h_minus=1.3, T_N=-2.8, T_K=0.5)
    direction = dict(dA_c=0.11, du_star=-0.07, dh_plus=0.05, dh_minus=-0.03, dT_N=0.13, dT_K=-0.09)
    analytic = pair.directional_derivative(**args, **direction)

    step = 2.0e-6
    plus = {name: args[name] + step * direction["d" + name] for name in args}
    minus = {name: args[name] - step * direction["d" + name] for name in args}
    plus_out = pair.solve_reference_amplitudes(**plus)
    minus_out = pair.solve_reference_amplitudes(**minus)
    fd_y = (plus_out["squared_amplitudes"] - minus_out["squared_amplitudes"]) / (2.0 * step)
    fd_a = (plus_out["amplitudes"] - minus_out["amplitudes"]) / (2.0 * step)

    np.testing.assert_allclose(analytic["d_squared_amplitudes"], fd_y, rtol=2e-9, atol=2e-10)
    np.testing.assert_allclose(analytic["d_amplitudes"], fd_a, rtol=2e-9, atol=2e-10)
    np.testing.assert_allclose(analytic["differentiated_identity_error"], 0.0, rtol=0.0, atol=2e-15)


def test_primary_tangent_assembly_keeps_sigma_axis_distinct():
    pair = KokunoSourceSignedCovariancePair()
    out = pair.solve_reference_amplitudes(**_args())
    prototypes = np.array(
        [
            [[1.0, 2.0, 0.0], [0.0, -1.0, 3.0]],
            [[-2.0, 1.0, 0.5], [4.0, 0.0, -1.0]],
        ]
    )
    epsilon = np.array([0.25, 0.36])
    amplitudes = np.broadcast_to(out["amplitudes"], (2, 2))
    actual = pair.assemble_primary_tangent_wave(prototypes, amplitudes, epsilon)
    expected = np.sqrt(epsilon)[:, None] * np.sum(amplitudes[..., :, None] * prototypes, axis=-2)
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2e-15)


def test_payload_records_sigma_not_m_and_keeps_actual_family_false(tmp_path):
    pair = KokunoSourceSignedCovariancePair()
    payload = pair.to_payload()
    truth = payload["truth_boundary"]
    assert truth["source_rectangle_sign_axis_identified"] is True
    assert truth["rectangle_sign_sigma_distinct_from_fourier_harmonic_m"] is True
    assert truth["source_reference_two_column_covariance_map_executable"] is True
    assert truth["actual_positive_order_background_bound"] is False
    assert truth["actual_source_h_sigma_pulse_integrals_bound"] is False
    assert truth["public_xyz_t_velocity_correction_materialized"] is False
    assert truth["genuinely_independent_second_covariance_column_ready"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False

    path = pair.save_json(tmp_path / "signed_pair.json")
    loaded = KokunoSourceSignedCovariancePair.load_json(path)
    assert loaded.sha256 == pair.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["public_xyz_t_velocity_correction_materialized"] = True
    with pytest.raises(ValueError, match="truth_boundary"):
        KokunoSourceSignedCovariancePair.from_payload(tampered)


def test_source_formula_provenance_is_pinned_and_not_autonomous_band_contrast():
    payload = KokunoSourceSignedCovariancePair().to_payload()
    assert payload["source"]["commit"] == "143f6773feb424ad9ed3a8d116653200f20346b7"
    formulas = payload["source"]["formulas"]
    assert "sigma=+/-" in formulas["rectangle_sign_fields"]
    assert "det(H_ref)=-2 A_c u_* h_+ h_-" == formulas["determinant"]
    assert "delta_band" not in " ".join(formulas.values())
