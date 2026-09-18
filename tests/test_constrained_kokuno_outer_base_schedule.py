import json
import math

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_outer_base_schedule import KokunoOuterBaseSchedule


def test_rf40_endpoint_amplitudes_and_stage_values_are_source_consistent():
    profile = KokunoOuterBaseSchedule()
    schedule = profile.outer_schedule
    eta = 0.37
    f = float(schedule.source_f(eta))

    at_XR = profile.profile_values_logX(profile.log_X_R, eta)
    np.testing.assert_allclose(
        at_XR["E"], math.exp(schedule.log_P_star) * f, rtol=2e-14, atol=0.0
    )
    np.testing.assert_allclose(at_XR["U"], 4.0 * eta, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(at_XR["ell"], 0.6, rtol=0.0, atol=0.0)

    at_stage1_end = profile.profile_values_logX(profile.log_stage1_end, eta)
    np.testing.assert_allclose(
        at_stage1_end["E"], math.exp(schedule.log_P1) * f, rtol=2e-14, atol=0.0
    )
    np.testing.assert_allclose(at_stage1_end["U"], 4.0 * eta, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(at_stage1_end["ell"], 0.0, rtol=0.0, atol=0.0)

    at_stage2_end = profile.profile_values_logX(profile.log_stage2_end, eta)
    expected_stage2_E = math.exp(schedule.log_P1 - 0.5 * schedule.T_d) * f
    np.testing.assert_allclose(at_stage2_end["E"], expected_stage2_E, rtol=2e-14, atol=0.0)
    np.testing.assert_allclose(at_stage2_end["U"], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(at_stage2_end["ell"], 0.0, rtol=0.0, atol=0.0)

    at_Xw = profile.profile_values_logX(profile.log_X_w, eta)
    np.testing.assert_allclose(
        at_Xw["E"], math.exp(schedule.log_e_w) * f, rtol=2e-14, atol=0.0
    )
    np.testing.assert_allclose(at_Xw["U"], 0.0, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(
        at_Xw["ell"], -schedule.lambda_outer, rtol=0.0, atol=0.0
    )


def test_power_stage_matches_existing_reserved_patch_background():
    profile = KokunoOuterBaseSchedule()
    schedule = profile.outer_schedule
    low, high = schedule.reserved_log_intervals()["I2"]
    log_X = 0.5 * (low + high)
    eta = -0.23
    values = profile.profile_values_logX(log_X, eta)

    expected_log_E = (
        schedule.log_c_patch
        + math.log(float(schedule.source_f(eta)))
        + (-0.5 - schedule.lambda_outer) * log_X
    )
    np.testing.assert_allclose(values["E"], math.exp(expected_log_E), rtol=3e-13, atol=0.0)
    np.testing.assert_allclose(values["U"], 0.0, rtol=0.0, atol=0.0)
    assert values["stage"].item() == "power_law"


def test_logarithmic_radial_derivatives_agree_with_independent_centered_differences():
    profile = KokunoOuterBaseSchedule(quadrature_points=96)
    schedule = profile.outer_schedule
    samples = [
        profile.log_X_R + 0.37,
        profile.log_stage1_end + 2.5,
        profile.log_stage2_end + 0.41,
        profile.log_X_w + 3.0,
    ]
    eta = 0.31
    eps = 2.0e-6
    for log_X in samples:
        center = profile.profile_values_logX(log_X, eta)
        plus = profile.profile_values_logX(log_X + eps, eta)
        minus = profile.profile_values_logX(log_X - eps, eta)
        for name, derivative_name in [
            ("E", "D_X_E"),
            ("F", "D_X_F"),
            ("U", "D_X_U"),
            ("v0", "D_X_v0"),
        ]:
            fd = (float(plus[name]) - float(minus[name])) / (2.0 * eps)
            expected = float(center[derivative_name])
            scale = max(1.0e-14, abs(fd), abs(expected))
            assert abs(fd - expected) / scale < 3.0e-7

    # The source endpoint relation remains exact and is not a quadrature artifact.
    assert math.isclose(schedule.log_P1, schedule.log_P_star - 0.2, rel_tol=0.0, abs_tol=0.0)


def test_eta_derivatives_agree_with_independent_centered_differences():
    profile = KokunoOuterBaseSchedule()
    log_X = profile.log_stage1_end + 4.0
    eta = 0.28
    eps = 2.0e-6
    center = profile.profile_values_logX(log_X, eta)
    plus = profile.profile_values_logX(log_X, eta + eps)
    minus = profile.profile_values_logX(log_X, eta - eps)
    for name, derivative_name in [
        ("E", "E_eta"),
        ("F", "F_eta"),
        ("U", "U_eta"),
        ("v0", "v0_eta"),
    ]:
        fd = (float(plus[name]) - float(minus[name])) / (2.0 * eps)
        expected = float(center[derivative_name])
        scale = max(1.0e-14, abs(fd), abs(expected))
        assert abs(fd - expected) / scale < 2.0e-7


def test_prefix_moment_ratio_obeys_incompressibility_ode_on_shutdown_stage():
    profile = KokunoOuterBaseSchedule(quadrature_points=128)
    log_X = profile.log_stage1_end + 6.0
    eps = 2.0e-6
    center = profile.profile_values_logX(log_X, 0.2)
    plus = profile.profile_values_logX(log_X + eps, 0.2)
    minus = profile.profile_values_logX(log_X - eps, 0.2)
    fd = (float(plus["m_ratio"]) - float(minus["m_ratio"])) / (2.0 * eps)
    rhs = float(center["k"] - center["m_ratio"])
    assert abs(fd - rhs) / max(1.0e-14, abs(rhs)) < 2.0e-7


def test_native_velocity_replays_cartesian_profile_map_on_materialized_outer_stage():
    profile = KokunoOuterBaseSchedule()
    log_X = profile.log_stage1_end + 3.0
    X = math.exp(log_X)
    q = 0.5
    eta = 0.3
    D = 0.5 - profile.h
    t = 1.0 - q * (1.0 - eta * eta)
    z = (q**D) * eta
    x = math.sqrt(2.0 * q * X)
    y = 0.0

    values = profile.profile_values(X, eta)
    expected = np.asarray(
        [
            float(values["v0"]) * x / (2.0 * q),
            (q ** (-1.0 - profile.h)) * float(values["F"]) * x,
            (q ** (-0.5 - profile.h)) * float(values["U"]),
        ]
    )
    actual = np.asarray(profile.velocity(x, y, z, t), dtype=float)
    np.testing.assert_allclose(actual, expected, rtol=2e-11, atol=0.0)
    assert np.all(np.isfinite(actual))
    assert np.linalg.norm(actual) > 0.0


def test_schedule_fails_closed_past_explicit_Tw_base_stage():
    profile = KokunoOuterBaseSchedule()
    with pytest.raises(ValueError, match="outer base schedule ends"):
        profile.profile_values_logX(profile.log_X_end + 1.0e-3, 0.0)


def test_payload_roundtrip_and_truth_boundary_fail_closed(tmp_path):
    profile = KokunoOuterBaseSchedule()
    path = tmp_path / "outer-base.json"
    profile.save_json(path)
    loaded = KokunoOuterBaseSchedule.load_json(path)
    assert loaded.to_payload() == profile.to_payload()
    assert loaded.sha256 == profile.sha256

    payload = json.loads(path.read_text(encoding="utf-8"))
    truth = payload["truth_boundary"]
    assert truth["outer_base_schedule_executable"] is True
    assert truth["cone_modulation_completed"] is False
    assert truth["global_leading_profile_reconstructed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False
    payload["truth_boundary"]["global_leading_profile_reconstructed"] = True
    with pytest.raises(ValueError, match="payload hash or content mismatch"):
        KokunoOuterBaseSchedule.from_payload(payload)
