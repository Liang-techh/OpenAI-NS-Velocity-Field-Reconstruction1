from __future__ import annotations

import copy
import json
import math

import pytest

from openai_ns_reconstruction.kokuno_pa10_chi_multiplier import (
    KokunoPA10SelectedChiMultiplier,
)
from openai_ns_reconstruction.kokuno_pa10_remainder_ball_bounds import (
    BallFactorBound,
    KokunoPA10RemainderBallBounds,
    KokunoPA10RemainderPrimitiveBounds,
)
from openai_ns_reconstruction.kokuno_rescaled_core_seed import (
    KokunoSourceRescaledCoreSeed,
)


def _calculator() -> KokunoPA10RemainderBallBounds:
    core = KokunoSourceRescaledCoreSeed()
    return KokunoPA10RemainderBallBounds(
        h=core.h,
        A=core.A,
        D=core.D,
        rescaling_lambda=core.rescaling_lambda,
    )


def _replace(
    primitives: KokunoPA10RemainderPrimitiveBounds,
    name: str,
    value: BallFactorBound,
) -> KokunoPA10RemainderPrimitiveBounds:
    payload = primitives.to_dict()
    payload[name] = {"norm": value.norm, "lipschitz": value.lipschitz}
    return KokunoPA10RemainderPrimitiveBounds.from_dict(payload)


def test_diagnostic_fixture_exercises_exact_remainder_monomials() -> None:
    calculator = _calculator()
    primitives = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
    result = calculator.remainder_envelopes(primitives)

    assert set(result["R1_terms"]) == {
        "coefficient_times_Phi",
        "W_times_Y_Phi_Y",
        "H_c_times_Phi_eta",
    }
    assert set(result["R2_terms"]) == {
        "linear_coefficient_times_u",
        "2A_eta_Lambda_inverse_u_squared",
        "W_times_Y_u_Y",
        "H_star_times_u_eta",
        "d_Lambda_inverse_u_u_eta",
        "4A_eta_p",
        "d_times_p_eta",
        "2_eta_Y_p_Y",
    }
    assert result["R1"]["norm"] > 0.0
    assert result["R1"]["lipschitz"] > 0.0
    assert result["R2"]["norm"] > 0.0
    assert result["R2"]["lipschitz"] > 0.0
    assert all(
        math.isfinite(value)
        for name in ("R1", "R2")
        for value in result[name].values()
    )


def test_product_uses_source_algebra_and_one_factor_lipschitz_rule() -> None:
    calculator = _calculator()
    first = BallFactorBound(2.0, 3.0)
    second = BallFactorBound(5.0, 7.0)
    third = BallFactorBound(11.0, 13.0)
    bound = calculator.product(first, second, third)
    P = calculator.product_constant

    assert bound.norm == pytest.approx(P**2 * 2.0 * 5.0 * 11.0)
    assert bound.lipschitz == pytest.approx(
        P**2 * (3.0 * 5.0 * 11.0 + 7.0 * 2.0 * 11.0 + 13.0 * 2.0 * 5.0)
    )


def test_zero_moving_lipschitz_produces_zero_remainder_lipschitz() -> None:
    calculator = _calculator()
    primitives = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
    payload = primitives.to_dict()
    for value in payload.values():
        value["lipschitz"] = 0.0
    frozen = KokunoPA10RemainderPrimitiveBounds.from_dict(payload)
    result = calculator.remainder_envelopes(frozen)

    assert result["R1"]["lipschitz"] == 0.0
    assert result["R2"]["lipschitz"] == 0.0


def test_explicit_log_radial_derivative_input_is_not_inferred() -> None:
    calculator = _calculator()
    primitives = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
    base = calculator.remainder_envelopes(primitives)
    changed = _replace(primitives, "Y_Phi_Y", BallFactorBound(4.0, 6.0))
    result = calculator.remainder_envelopes(changed)

    assert result["R1_terms"]["W_times_Y_Phi_Y"]["norm"] > base["R1_terms"][
        "W_times_Y_Phi_Y"
    ]["norm"]
    assert result["R1_terms"]["W_times_Y_Phi_Y"]["lipschitz"] > base[
        "R1_terms"
    ]["W_times_Y_Phi_Y"]["lipschitz"]
    assert result["R2"] == base["R2"]


def test_pressure_seam_is_explicit_and_monotone() -> None:
    calculator = _calculator()
    primitives = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
    base = calculator.remainder_envelopes(primitives)
    changed = _replace(primitives, "p", BallFactorBound(8.0, 9.0))
    result = calculator.remainder_envelopes(changed)

    assert result["R2_terms"]["4A_eta_p"]["norm"] > base["R2_terms"][
        "4A_eta_p"
    ]["norm"]
    assert result["R2_terms"]["4A_eta_p"]["lipschitz"] > base["R2_terms"][
        "4A_eta_p"
    ]["lipschitz"]
    assert result["R1"] == base["R1"]


def test_selected_chi_can_feed_only_diagnostic_conditional_MK() -> None:
    calculator = _calculator()
    primitives = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
    selected_chi = KokunoPA10SelectedChiMultiplier()
    selected_mchi = float(
        selected_chi.to_payload()["certificate"]["selected_multiplier_norm_chi_upper"]
    )
    pair = calculator.diagnostic_conditional_MK(
        primitives,
        multiplier_norm_chi=selected_mchi,
    )

    assert pair["conditional_M_pair_max_envelope"] > 0.0
    assert pair["conditional_K_pair_max_envelope"] > 0.0
    assert math.isfinite(pair["conditional_M_pair_max_envelope"])
    assert math.isfinite(pair["conditional_K_pair_max_envelope"])
    assert calculator.truth_boundary["source_operator_constant_M_machine_bound"] is False
    assert calculator.truth_boundary["source_operator_constant_K_machine_bound"] is False


def test_payload_roundtrip_and_truth_boundary_fail_closed(tmp_path) -> None:
    calculator = _calculator()
    primitives = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture()
    selected_chi = KokunoPA10SelectedChiMultiplier()
    selected_mchi = float(
        selected_chi.to_payload()["certificate"]["selected_multiplier_norm_chi_upper"]
    )
    path = tmp_path / "remainder.json"
    calculator.save_json(path, primitives, multiplier_norm_chi=selected_mchi)
    loaded, loaded_primitives = KokunoPA10RemainderBallBounds.load_json(path)

    assert loaded == calculator
    assert loaded_primitives == primitives
    assert loaded.to_payload(
        loaded_primitives, multiplier_norm_chi=selected_mchi
    ) == calculator.to_payload(primitives, multiplier_norm_chi=selected_mchi)

    payload = json.loads(path.read_text())
    tampered = copy.deepcopy(payload)
    tampered["truth_boundary"]["source_R1_radius_one_ball_norm_machine_bound"] = True
    unsigned = copy.deepcopy(tampered)
    unsigned.pop("sha256")
    import hashlib

    canonical = json.dumps(unsigned, sort_keys=True, separators=(",", ":"), allow_nan=False)
    tampered["sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    with pytest.raises(ValueError, match="truth boundary mismatch"):
        KokunoPA10RemainderBallBounds.from_payload(tampered)


def test_invalid_primitive_or_parameter_inputs_fail_closed() -> None:
    with pytest.raises(ValueError):
        BallFactorBound(-1.0)
    with pytest.raises(ValueError):
        KokunoPA10RemainderBallBounds(h=0.6, A=1.1, D=0.1, rescaling_lambda=1.0)

    payload = KokunoPA10RemainderPrimitiveBounds.diagnostic_unit_fixture().to_dict()
    payload.pop("p_eta")
    with pytest.raises(ValueError, match="primitive-bound fields mismatch"):
        KokunoPA10RemainderPrimitiveBounds.from_dict(payload)
