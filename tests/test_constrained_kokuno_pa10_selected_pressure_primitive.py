import copy
import json

import numpy as np
import pytest

from openai_ns_reconstruction.kokuno_pa10_selected_pressure_primitive import (
    KokunoPA10SelectedPressurePrimitive,
)


def test_selected_pressure_matches_existing_core_primitive_and_radial_identity():
    primitive = KokunoPA10SelectedPressurePrimitive(quadrature_points=64)
    Y = np.asarray([0.0, 0.35, 1.7, 4.1])
    eta = np.asarray([-0.75, -0.2, 0.35, 0.82])
    got = primitive.evaluate(Y, eta)
    core = primitive.core.profile_values_scaled(Y, eta)

    np.testing.assert_allclose(
        got["p"], core["scaled_pressure_correction"], rtol=3.0e-12, atol=2.0e-14
    )
    np.testing.assert_allclose(got["p_Y"], got["g"] ** 2 * got["Phi"] ** 2, rtol=0.0, atol=0.0)
    np.testing.assert_allclose(got["Y_p_Y"], Y * got["p_Y"], rtol=0.0, atol=0.0)


def test_selected_profile_derivatives_match_independent_centered_differences():
    primitive = KokunoPA10SelectedPressurePrimitive()
    Y = np.asarray([0.4, 1.2, 2.7])
    eta = np.asarray([-0.65, 0.12, 0.61])
    h = 2.0e-6
    got = primitive.evaluate(Y, eta)

    plus_y = primitive.evaluate(Y * np.exp(h), eta)["Phi"]
    minus_y = primitive.evaluate(Y * np.exp(-h), eta)["Phi"]
    fd_log_y = (plus_y - minus_y) / (2.0 * h)
    np.testing.assert_allclose(got["Y_Phi_Y"], fd_log_y, rtol=3.0e-7, atol=2.0e-10)

    plus_eta = primitive.evaluate(Y, eta + h)["Phi"]
    minus_eta = primitive.evaluate(Y, eta - h)["Phi"]
    fd_eta = (plus_eta - minus_eta) / (2.0 * h)
    np.testing.assert_allclose(got["Phi_eta"], fd_eta, rtol=3.0e-7, atol=2.0e-10)


def test_pressure_axis_initial_value_and_eta_derivative_are_exactly_zero():
    primitive = KokunoPA10SelectedPressurePrimitive(quadrature_points=64)
    eta = np.asarray([-0.9, -0.25, 0.0, 0.4, 0.95])
    got = primitive.evaluate(np.zeros_like(eta), eta)
    np.testing.assert_array_equal(got["p"], np.zeros_like(eta))
    np.testing.assert_array_equal(got["p_eta"], np.zeros_like(eta))
    np.testing.assert_array_equal(got["Y_p_Y"], np.zeros_like(eta))


def test_engineering_envelope_is_nontrivial_nested_and_fail_closed():
    primitive = KokunoPA10SelectedPressurePrimitive(quadrature_points=32)
    report = primitive.engineering_envelope(y_count=7, eta_count=9)
    for name in ("Phi", "Phi_eta", "Y_Phi_Y", "p", "p_eta", "Y_p_Y"):
        assert report["engineering_abs_envelope"][name] >= report["fine_sampled_abs_max"][name]
        assert report["engineering_abs_envelope"][name] >= report["coarse_sampled_abs_max"][name]
    assert report["engineering_abs_envelope"]["Phi"] > 0.0
    # The known selected phase maximizer is explicitly sampled so materialized
    # g underflow on generic eta nodes cannot masquerade as p == 0 everywhere.
    assert report["phase_stationary_eta_injected"] is True
    assert report["engineering_abs_envelope"]["p"] > 1.0
    assert report["engineering_abs_envelope"]["Y_p_Y"] > 1.0
    assert report["coarse_grid"]["eta_count_actual"] >= report["coarse_grid"]["eta_count_requested"]
    assert report["fine_grid"]["eta_count_actual"] >= report["fine_grid"]["eta_count_requested"]
    assert report["continuum_supremum_certified"] is False
    assert report["source_radius_one_ball_bound"] is False
    truth = primitive.truth_boundary
    assert truth["selected_center_pressure_primitive_executable"] is True
    assert truth["source_pressure_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["selected_pa16_handoff_allowed"] is False
    assert truth["pde_validated"] is False


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path):
    primitive = KokunoPA10SelectedPressurePrimitive(quadrature_points=32)
    path = tmp_path / "selected-pressure.json"
    primitive.save(path)
    restored = KokunoPA10SelectedPressurePrimitive.load(path)
    assert restored.sha256 == primitive.sha256
    np.testing.assert_allclose(
        restored.evaluate(1.3, 0.1)["p"], primitive.evaluate(1.3, 0.1)["p"]
    )

    payload = primitive.to_payload()
    bad = copy.deepcopy(payload)
    bad["truth_boundary"]["source_pressure_radius_one_ball_norm_machine_bound"] = True
    unsigned = {key: value for key, value in bad.items() if key != "sha256"}
    import hashlib

    bad["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    ).hexdigest()
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoPA10SelectedPressurePrimitive.from_payload(bad)


def test_domain_guards():
    primitive = KokunoPA10SelectedPressurePrimitive()
    with pytest.raises(ValueError):
        primitive.evaluate(-1.0e-6, 0.0)
    with pytest.raises(ValueError):
        primitive.evaluate(4.100001, 0.0)
    with pytest.raises(ValueError):
        primitive.evaluate(1.0, 1.000001)
