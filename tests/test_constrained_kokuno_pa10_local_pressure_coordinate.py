import copy
import hashlib
import json

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from openai_ns_reconstruction.kokuno_pa10_local_pressure_coordinate import (
    KokunoPA10LocalPressureCoordinate,
)


def _fd6(values, h):
    fm3, fm2, fm1, fp1, fp2, fp3 = values
    return (-fm3 + 9.0 * fm2 - 45.0 * fm1 + 45.0 * fp1 - 9.0 * fp2 + fp3) / (
        60.0 * h
    )


def test_local_coordinate_center_matches_public_pressure_center():
    local = KokunoPA10LocalPressureCoordinate(radial_quadrature_points=64)
    Y = np.asarray([0.0, 0.7, 2.0, 4.1])
    got = local.evaluate_local(Y, np.zeros_like(Y))
    public = local.primitive.evaluate(Y, np.full_like(Y, local.eta_star))

    np.testing.assert_allclose(got["p"], public["p"], rtol=2.0e-13, atol=2.0e-14)
    np.testing.assert_allclose(got["Phi"], public["Phi"], rtol=2.0e-13, atol=2.0e-14)
    np.testing.assert_array_equal(got["p_s"], np.zeros_like(Y))
    np.testing.assert_array_equal(got["p_eta"], np.zeros_like(Y))
    np.testing.assert_array_equal(got["log_g"], np.zeros_like(Y))


def test_local_layer_is_nontrivial_and_not_flattened_by_eta_materialization():
    local = KokunoPA10LocalPressureCoordinate()
    Y = np.full(5, 2.0)
    s = np.asarray([-0.8, -0.4, 0.0, 0.4, 0.8])
    got = local.evaluate_local(Y, s)

    assert got["p"][2] == pytest.approx(2.0, rel=0.0, abs=2.0e-13)
    assert np.all(got["p"][[0, 1, 3, 4]] < got["p"][2])
    assert np.all(got["log_g"][[0, 1, 3, 4]] < 0.0)
    assert float(np.max(np.abs(got["p_s"]))) > 1.0e-6

    hi, lo = local.eta_two_term(s)
    recovered = (hi - local.eta_star) + lo
    np.testing.assert_allclose(recovered, local.eta_offset(s), rtol=0.0, atol=0.0)
    assert np.any(lo != 0.0)


def test_local_pressure_derivative_has_fd6_refinement():
    local = KokunoPA10LocalPressureCoordinate(radial_quadrature_points=64)
    Y = np.asarray([0.7, 1.2, 1.9, 2.6, 3.2, 3.8])
    s = np.asarray([-0.83, -0.51, -0.19, 0.14, 0.43, 0.79])
    reference = local.evaluate_local(Y, s)["p_s"]

    errors = []
    for h in (0.04, 0.02, 0.01):
        values = [
            local.evaluate_local(Y, s + offset * h)["p"]
            for offset in (-3.0, -2.0, -1.0, 1.0, 2.0, 3.0)
        ]
        fd = _fd6(values, h)
        rms = float(np.sqrt(np.mean((fd - reference) ** 2)))
        scale = float(np.sqrt(np.mean(reference**2)))
        errors.append(rms / scale)

    assert errors[-1] <= 2.0e-5
    assert errors[0] / errors[1] >= 8.0
    assert errors[1] / errors[2] >= 8.0


def test_local_s_fundamental_theorem_closes_without_eta_roundtrip():
    local = KokunoPA10LocalPressureCoordinate(radial_quadrature_points=64)
    intervals = ((-0.91, 0.31), (-0.62, 0.77), (-0.23, 0.94))
    y_values = (0.9, 2.1, 3.7)

    rms_by_order = []
    max_by_order = []
    for order in (16, 32, 64):
        nodes, weights = leggauss(order)
        errors = []
        refs = []
        for y_value, (left, right) in zip(y_values, intervals, strict=True):
            half = 0.5 * (right - left)
            midpoint = 0.5 * (right + left)
            points = midpoint + half * nodes
            integral = half * float(
                np.dot(weights, local.evaluate_local(y_value, points)["p_s"])
            )
            endpoint = float(
                local.evaluate_local(y_value, right)["p"]
                - local.evaluate_local(y_value, left)["p"]
            )
            errors.append(integral - endpoint)
            refs.append(endpoint)
        errors = np.asarray(errors)
        refs = np.asarray(refs)
        rms_by_order.append(
            float(np.sqrt(np.mean(errors**2)) / np.sqrt(np.mean(refs**2)))
        )
        max_by_order.append(float(np.max(np.abs(errors)) / np.max(np.abs(refs))))

    assert rms_by_order[-1] <= 1.0e-8
    assert max_by_order[-1] <= 5.0e-8


def test_physical_chain_rule_is_exact_in_public_local_contract():
    local = KokunoPA10LocalPressureCoordinate()
    got = local.evaluate_local(
        np.asarray([0.5, 1.3, 2.4, 3.9]),
        np.asarray([-0.7, -0.2, 0.35, 0.88]),
    )
    np.testing.assert_array_equal(got["p_eta"], local.sqrt_lambda * got["p_s"])
    np.testing.assert_array_equal(got["Phi_eta"], local.sqrt_lambda * got["Phi_s"])


def test_receipt_and_truth_boundary_remain_fail_closed():
    local = KokunoPA10LocalPressureCoordinate(radial_quadrature_points=64)
    receipt = local.diagnostic_receipt()
    assert receipt["selfcheck_passed"] is True
    assert receipt["independent_agent4_audit"] == "required"
    assert receipt["heldout_ns_residual_assessed"] is False
    assert receipt["pde_validated"] is False
    truth = local.truth_boundary
    assert truth["selected_pressure_local_coordinate_executable"] is True
    assert truth["independent_agent4_local_coordinate_audit_required"] is True
    assert truth["source_pressure_radius_one_ball_norm_machine_bound"] is False
    assert truth["selected_pa16_handoff_allowed"] is False
    assert truth["global_pressure_matched"] is False
    assert truth["pde_validated"] is False


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path):
    local = KokunoPA10LocalPressureCoordinate(radial_quadrature_points=32)
    path = tmp_path / "local-pressure-coordinate.json"
    local.save(path)
    restored = KokunoPA10LocalPressureCoordinate.load(path)
    assert restored.sha256 == local.sha256
    np.testing.assert_allclose(
        restored.evaluate_local(1.7, 0.4)["p"], local.evaluate_local(1.7, 0.4)["p"]
    )

    bad = copy.deepcopy(local.to_payload())
    bad["truth_boundary"]["source_pressure_radius_one_ball_norm_machine_bound"] = True
    unsigned = {key: value for key, value in bad.items() if key != "sha256"}
    bad["sha256"] = hashlib.sha256(
        json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False).encode(
            "utf-8"
        )
    ).hexdigest()
    with pytest.raises(ValueError, match="truth-boundary"):
        KokunoPA10LocalPressureCoordinate.from_payload(bad)


def test_domain_guards():
    local = KokunoPA10LocalPressureCoordinate()
    with pytest.raises(ValueError):
        local.evaluate_local(-1.0e-6, 0.0)
    with pytest.raises(ValueError):
        local.evaluate_local(4.100001, 0.0)
    outside = (1.0 - abs(local.eta_star)) * local.sqrt_lambda * 1.001
    with pytest.raises(ValueError):
        local.evaluate_local(1.0, outside)
