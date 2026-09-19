from __future__ import annotations

import copy
import hashlib
import json
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_operator_primitives import (
    KokunoPA10OperatorPrimitiveBounds,
)


def test_source_convolution_and_product_constants_are_executable() -> None:
    bounds = KokunoPA10OperatorPrimitiveBounds()
    expected_c_sq = 4.0 * math.pi * math.pi / 3.0
    assert bounds.C_sq == pytest.approx(expected_c_sq, rel=2e-15)
    assert bounds.product_constant == pytest.approx(expected_c_sq**2, rel=2e-15)


def test_J_nu_vanishing_bounds_follow_displayed_coefficients() -> None:
    bounds = KokunoPA10OperatorPrimitiveBounds()
    assert bounds.j_nu_bound(1, 0) == 80.0
    assert bounds.j_nu_bound(2, 0) == 40.0
    assert bounds.j_nu_bound(2, 1) == pytest.approx(80.0 / 6.0)
    assert bounds.j_nu_bound(2, 4) < bounds.j_nu_bound(2, 1)


def test_mixed_derivative_template_keeps_rho_explicit() -> None:
    bounds = KokunoPA10OperatorPrimitiveBounds()
    rho = 0.25
    expected = (80.0 / rho) * bounds.product_constant
    assert bounds.mixed_derivative_template_constant(rho) == pytest.approx(expected)


def test_T_power_and_inverse_absolute_series_are_reproducible() -> None:
    bounds = KokunoPA10OperatorPrimitiveBounds()
    assert bounds.t_power_bound(0, 1.0) == 1.0
    assert bounds.t_power_bound(1, 1.0) == 20.0
    assert bounds.t_power_bound(2, 1.0) == pytest.approx(400.0 / 3.0)
    assert bounds.inverse_one_plus_t_bound(0.0) == 1.0

    # Independent direct recurrence for the same positive majorant at M_chi=1.
    term = 1.0
    expected = 1.0
    for k in range(1, 200):
        term *= 40.0 / (k * (k + 1.0))
        expected += term
        if term <= 1.0e-16 * expected:
            break
    actual = bounds.inverse_one_plus_t_bound(1.0)
    assert actual == pytest.approx(expected, rel=5e-15)
    assert actual > 1.0


def test_conditional_pair_envelope_propagates_without_source_promotion() -> None:
    bounds = KokunoPA10OperatorPrimitiveBounds()
    result = bounds.conditional_pair_envelopes(
        multiplier_norm_chi=0.0,
        R1_ball_norm=2.0,
        R2_ball_norm=3.0,
        R1_lipschitz_norm=5.0,
        R2_lipschitz_norm=7.0,
    )
    assert result["inverse_one_plus_T_absolute_bound"] == 1.0
    assert result["conditional_M_phi_component"] == pytest.approx(40.0)
    assert result["conditional_M_u_component"] == pytest.approx(120.0)
    assert result["conditional_M_pair_max_envelope"] == pytest.approx(120.0)
    assert result["conditional_K_phi_component"] == pytest.approx(100.0)
    assert result["conditional_K_u_component"] == pytest.approx(280.0)
    assert result["conditional_K_pair_max_envelope"] == pytest.approx(280.0)

    gate = bounds.diagnostic_fixed_point_gate_from_remainder_bounds(
        multiplier_norm_chi=0.0,
        R1_ball_norm=2.0,
        R2_ball_norm=3.0,
        R1_lipschitz_norm=5.0,
        R2_lipschitz_norm=7.0,
    )
    assert gate.operator_constant_M == pytest.approx(120.0)
    assert gate.operator_constant_K == pytest.approx(280.0)
    assert gate.conditional_threshold == pytest.approx(560.0)
    truth = gate.report()["truth_boundary"]
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["source_fixed_point_distance_machine_bound"] is False


def test_invalid_source_dependent_inputs_fail_closed() -> None:
    bounds = KokunoPA10OperatorPrimitiveBounds()
    with pytest.raises(ValueError, match="nu"):
        bounds.j_nu_bound(3)
    with pytest.raises(TypeError, match="vanishing_order"):
        bounds.j_nu_bound(2, 1.5)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="vanishing_order"):
        bounds.j_nu_bound(2, -1)
    with pytest.raises(ValueError, match="rho"):
        bounds.mixed_derivative_template_constant(0.0)
    with pytest.raises(ValueError, match="multiplier_norm_chi"):
        bounds.t_power_bound(1, -1.0)
    with pytest.raises(TypeError, match="k"):
        bounds.t_power_bound(True, 1.0)  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="R1_ball_norm"):
        bounds.conditional_pair_envelopes(
            multiplier_norm_chi=1.0,
            R1_ball_norm=-1.0,
            R2_ball_norm=1.0,
            R1_lipschitz_norm=1.0,
            R2_lipschitz_norm=1.0,
        )


def test_truth_boundary_exposes_only_universal_machine_bound_progress() -> None:
    truth = KokunoPA10OperatorPrimitiveBounds().report()["truth_boundary"]
    assert truth["source_universal_operator_primitives_executable"] is True
    assert truth["source_product_constant_executable"] is True
    assert truth["source_J_nu_bound_executable"] is True
    assert truth["source_T_power_bound_executable"] is True
    assert truth["conditional_MK_envelope_executable"] is True
    assert truth["source_rho_machine_bound"] is False
    assert truth["source_chi_multiplier_norm_machine_bound"] is False
    assert truth["source_R1_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_R2_radius_one_ball_norm_machine_bound"] is False
    assert truth["source_operator_constant_M_machine_bound"] is False
    assert truth["source_operator_constant_K_machine_bound"] is False
    assert truth["source_B0_dependencies_machine_bound"] is False
    assert truth["source_T_sh_lower_bound_verified"] is False
    assert truth["selected_pa16_handoff_allowed"] is False
    assert truth["pde_validated"] is False
    assert truth["paper_exact"] is False


def test_payload_roundtrip_and_truth_tamper_fail_closed(tmp_path) -> None:
    bounds = KokunoPA10OperatorPrimitiveBounds()
    payload = bounds.to_payload()
    replay = KokunoPA10OperatorPrimitiveBounds.from_payload(payload)
    assert replay.sha256 == bounds.sha256
    assert replay.to_payload() == payload

    path = tmp_path / "operator_primitives.json"
    bounds.save_json(path)
    loaded = KokunoPA10OperatorPrimitiveBounds.load_json(path)
    assert loaded.sha256 == bounds.sha256
    assert json.loads(path.read_text())["sha256"] == bounds.sha256

    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_operator_constant_M_machine_bound"] = True
    body = copy.deepcopy(tampered)
    body.pop("sha256")
    canonical = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False)
    tampered["sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    with pytest.raises(ValueError, match="truth boundary"):
        KokunoPA10OperatorPrimitiveBounds.from_payload(tampered)
