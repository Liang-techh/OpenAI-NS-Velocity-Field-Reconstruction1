import copy
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_rescaled_core_seed import (
    KokunoSourceRescaledCoreSeed,
)


def test_source_scale_shifted_center_stays_finite_and_nontrivial():
    seed = KokunoSourceRescaledCoreSeed()
    assert seed.pressure_scale > 1.0e13
    assert seed.rescaling_lambda > 1.0e26

    eta = np.array([-1.0, seed.phase_stationary_eta, 0.0, 1.0])
    values = seed.profile_values_scaled(np.full(eta.shape, 2.0), eta)
    for key in (
        "Phi_0",
        "u_0",
        "U",
        "U_eta",
        "v0",
        "Pi",
        "log_g",
        "log_F",
    ):
        assert np.all(np.isfinite(values[key])), key

    # The log representation remains informative even when ordinary binary64
    # angular amplitude underflows away from the real-axis phase maximum.
    assert np.max(values["log_g"]) == pytest.approx(0.0, abs=1.0e-12)
    assert np.max(values["F"]) > 0.1
    assert np.max(np.abs(values["U"])) > 1.0e-3
    assert np.all(values["Pi"] < 0.0)


def test_center_reconstructs_exact_source_shift_and_pressure_identity():
    seed = KokunoSourceRescaledCoreSeed()
    eta = np.array([seed.phase_stationary_eta, 0.25])
    Y = np.array([1.25, 3.0])
    values = seed.profile_values_scaled(Y, eta)

    expected_U = values["U_star"] + values["u_0"] / seed.rescaling_lambda
    assert np.allclose(values["U"], expected_U, rtol=2.0e-15, atol=0.0)
    assert np.allclose(values["Pi_X"], values["F"] ** 2, rtol=0.0, atol=0.0)
    assert np.all(values["physical_pressure_correction"] >= 0.0)
    assert np.all(values["physical_pressure_correction"] <= 4.1 / seed.rescaling_lambda)

    axis = seed.profile_values_scaled(np.zeros(eta.shape), eta)
    assert np.allclose(axis["Phi_0"], 1.0, rtol=0.0, atol=0.0)
    assert np.allclose(axis["U"], axis["U_star"], rtol=0.0, atol=0.0)
    assert np.allclose(axis["Pi"], axis["Pi_0"], rtol=0.0, atol=0.0)


def test_center_radial_derivatives_are_analytic_not_finite_difference_state():
    seed = KokunoSourceRescaledCoreSeed()
    eta = seed.phase_stationary_eta
    Y = 1.7
    values = seed.profile_values_scaled(np.asarray(Y), np.asarray(eta))

    eps = 2.0e-6
    plus = seed.profile_values_scaled(np.asarray(Y + eps), np.asarray(eta))
    minus = seed.profile_values_scaled(np.asarray(Y - eps), np.asarray(eta))
    finite_DF = Y * (float(plus["F"]) - float(minus["F"])) / (2.0 * eps)
    finite_DU = Y * (float(plus["U"]) - float(minus["U"])) / (2.0 * eps)
    assert float(values["D_X_F"]) == pytest.approx(finite_DF, rel=2.0e-7, abs=1.0e-12)
    assert float(values["D_X_U"]) == pytest.approx(finite_DU, rel=3.0e-5, abs=1.0e-15)


def test_vectorized_native_velocity_is_axis_regular_and_fail_closed():
    seed = KokunoSourceRescaledCoreSeed()
    axis_velocity = seed.velocity(0.0, 0.0, 0.0, 0.0)
    assert axis_velocity.shape == (3,)
    assert np.all(np.isfinite(axis_velocity))
    assert axis_velocity[0] == pytest.approx(0.0, abs=0.0)
    assert axis_velocity[1] == pytest.approx(0.0, abs=0.0)
    assert axis_velocity[2] == pytest.approx(seed.j0, rel=2.0e-15)

    eta = seed.phase_stationary_eta
    q = 1.0 / (1.0 - eta * eta)
    z = q ** seed.D * eta
    X = 2.0 / seed.rescaling_lambda
    r = math.sqrt(2.0 * q * X)
    velocity = seed.velocity(np.array([r, 0.5 * r]), 0.0, z, 0.0)
    assert velocity.shape == (2, 3)
    assert np.all(np.isfinite(velocity))
    assert np.linalg.norm(velocity[0]) > 1.0e-3

    outside_r = math.sqrt(2.0 * q * (4.2 / seed.rescaling_lambda))
    with pytest.raises(ValueError, match="0<=Y<=4.1"):
        seed.velocity(outside_r, 0.0, z, 0.0)


def test_serialization_replays_and_truth_boundary_fails_closed(tmp_path):
    seed = KokunoSourceRescaledCoreSeed()
    payload = seed.to_payload()
    replay = KokunoSourceRescaledCoreSeed.from_payload(payload)
    assert replay.sha256 == seed.sha256
    assert replay.report()["all_sample_shifted_fields_finite"] is True
    assert replay.to_payload()["truth_boundary"]["source_fixed_point_solved"] is False
    assert replay.to_payload()["truth_boundary"]["global_leading_profile_reconstructed"] is False
    assert replay.to_payload()["truth_boundary"]["paper_exact"] is False

    path = tmp_path / "rescaled-core-seed.json"
    seed.save_json(path)
    assert KokunoSourceRescaledCoreSeed.load_json(path).sha256 == seed.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_fixed_point_solved"] = True
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoSourceRescaledCoreSeed.from_payload(tampered)


def test_report_records_autonomous_scale_without_upgrading_global_claims():
    seed = KokunoSourceRescaledCoreSeed()
    report = seed.report()
    assert report["lambda_to_pressure_square_ratio"] == pytest.approx(1.0)
    assert report["rescaling_lambda"] > 1.0e26
    assert report["all_sample_shifted_fields_finite"] is True
    truth = report["truth_boundary"]
    assert truth["selected_rescaling_lambda_is_autonomous_conditioning_choice"] is True
    assert truth["source_contraction_threshold_verified"] is False
    assert truth["source_complex_C_bound_verified"] is False
    assert truth["source_pressure_datum_applied_to_selected_global_path"] is False
    assert truth["pde_validated"] is False
