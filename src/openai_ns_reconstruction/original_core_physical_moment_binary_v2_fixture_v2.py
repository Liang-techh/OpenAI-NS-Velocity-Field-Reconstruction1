"""Synthetic checks for the signed-dyadic v2 scale ledger."""

from __future__ import annotations

from typing import Any, Callable


SCHEMA = "original_core_physical_moment_binary_v2_fixture_v2"


def _endpoint(sign: int, mantissa: int, exponent: int) -> list[int]:
    return [sign, mantissa, exponent, mantissa.bit_length()]


def _positive(mantissa: int = 1, exponent: int = 0) -> list[int]:
    return _endpoint(0, mantissa, exponent)


def _negative(mantissa: int = 1, exponent: int = 0) -> list[int]:
    return _endpoint(1, mantissa, exponent)


def _row(value: list[int], derivative: list[int], family: str = "fixture-family") -> dict[str, Any]:
    return {
        "z0": {"finite": [value, value], "error_upper": _positive(0)},
        "z1": {"finite": [derivative, derivative], "error_upper": _positive(0)},
        "source_family": family,
    }


def _moments(family: str = "fixture-family") -> dict[str, Any]:
    values = {
        "H": (_positive(1, 4), _positive(1, 3)),
        "M": (_positive(1, 5), _positive(1, 2)),
        "K": (_positive(1, 6), _positive(1, 4)),
        "A": (_positive(1, 7), _negative(1, 3)),
        "B": (_positive(1, 8), _positive(1, 5)),
        "C": (_positive(1, 9), _positive(1, 6)),
    }
    return {name: _row(value, derivative, family) for name, (value, derivative) in values.items()}


def _pair_value(value: list[int]) -> tuple[int, int]:
    signed, exponent = value
    return int(signed), int(exponent)


def _cmp_pair(left: list[int], right: list[int]) -> int:
    left_signed, left_exponent = _pair_value(left)
    right_signed, right_exponent = _pair_value(right)
    if left_signed == 0:
        return 0 if right_signed == 0 else (-1 if right_signed > 0 else 1)
    if right_signed == 0:
        return 1 if left_signed > 0 else -1
    left_sign = -1 if left_signed < 0 else 1
    right_sign = -1 if right_signed < 0 else 1
    if left_sign != right_sign:
        return -1 if left_sign < right_sign else 1
    left_top = abs(left_signed).bit_length() + left_exponent
    right_top = abs(right_signed).bit_length() + right_exponent
    if left_top != right_top:
        result = 1 if left_top > right_top else -1
    else:
        shift = left_exponent - right_exponent
        if shift >= 0:
            left_m, right_m = abs(left_signed) << shift, abs(right_signed)
        else:
            left_m, right_m = abs(left_signed), abs(right_signed) << (-shift)
        result = (left_m > right_m) - (left_m < right_m)
    return result if left_sign > 0 else -result


def _contains(interval: dict[str, Any], target: list[int]) -> bool:
    return _cmp_pair(interval["lower"], target) <= 0 <= _cmp_pair(interval["upper"], target)


def _exact(interval: dict[str, Any], target: list[int]) -> bool:
    return interval["lower"] == target and interval["upper"] == target


def _raises(call: Callable[[], Any]) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def run_fixture(scale: Callable[..., Any]) -> dict[str, Any]:
    """Run hand, extreme-exponent, tail, and rejection checks."""

    moments = _moments()
    log_derivative = {
        "finite": [_negative(3, 0), _negative(1, 0)],
        "source_family": "fixture-family",
    }
    result = scale(moments, log_derivative, source_family="fixture-family")
    checks: dict[str, bool] = {}
    checks["schema"] = result.get("schema") == "original_core_physical_moment_binary_v2"
    expected = {
        "I_V": [1, 7],
        "I_V2": [1, 9],
        "I_RF": [1, 7],
        "I_RFV": [1, 9],
        "I_RF2": [1, 12],
        "I_F2": [1, 11],
    }
    checks["hand_finite_values"] = all(
        _contains(result["finite"][name], value) for name, value in expected.items()
    )
    expected_z1 = {
        "I_V": [1, 4],
        "I_V2": [-1, 5],
        "I_RF": [-5, 6],
        "I_RFV": [-11, 7],
        "I_RF2": [-47, 9],
        "I_F2": [-47, 8],
    }
    checks["hand_Z1_intervals"] = all(
        _contains(result["finite_Z1"][name], value) for name, value in expected_z1.items()
    )
    checks["scale_powers"] = all(
        result["scale"][name] == value
        for name, value in {
            "I_V": [1, 0], "I_V2": [1, 0], "I_RF": [2, 1],
            "I_RFV": [2, 1], "I_RF2": [2, 2], "I_F2": [1, 2],
        }.items()
    )
    checks["zero_errors"] = all(result["error_upper"][name] == [0, 0] for name in expected)

    noisy = _moments()
    for row in noisy.values():
        row["z0"]["error_upper"] = _positive(1, -2)
        row["z1"]["error_upper"] = _positive(1, -3)
    noisy_result = scale(noisy, log_derivative, source_family="fixture-family")
    # I_RF2 has k=2, c=16: 1/8 + 2*3/4 = 13/8, hence error 26.
    checks["nonzero_tail_k2"] = (
        _cmp_pair(noisy_result["error_upper_Z1"]["I_RF2"], [13, 1]) == 0
        and _cmp_pair(
            noisy_result["enclosure_Z1"]["I_RF2"]["lower"],
            noisy_result["finite_Z1"]["I_RF2"]["lower"],
        ) < 0
    )

    interval_moments = _moments()
    interval_moments["H"]["z0"]["finite"] = [_positive(1, 3), _positive(1, 4)]
    interval_result = scale(interval_moments, log_derivative, source_family="fixture-family")
    checks["input_interval_preserved"] = _contains(interval_result["finite"]["I_RF"], [1, 6])

    add = getattr(scale, "interval_add", None)
    mul = getattr(scale, "interval_mul", None)
    checks["extreme_helpers_exposed"] = callable(add) and callable(mul)
    if checks["extreme_helpers_exposed"]:
        huge = _positive(1, 10**18)
        one = _positive(1, 0)
        tiny = _positive(1, -10**18)
        neg_huge = _negative(1, 10**18)
        point_huge = {"lower": huge, "upper": huge}
        point_tiny = {"lower": tiny, "upper": tiny}
        point_neg_huge = {"lower": neg_huge, "upper": neg_huge}
        residual = add(add(point_huge, {"lower": one, "upper": one}), point_neg_huge)
        zero_tiny = add(
            {"lower": _positive(0), "upper": _positive(0)},
            point_tiny,
        )
        near_cancel = add(
            {"lower": _endpoint(0, 3, 600), "upper": _endpoint(0, 3, 600)},
            {"lower": _endpoint(1, 2, 600), "upper": _endpoint(1, 2, 600)},
        )
        wide_near_cancel = add(
            {
                "lower": _endpoint(0, (1 << 799) + 1, 0),
                "upper": _endpoint(0, (1 << 799) + 1, 0),
            },
            {
                "lower": _endpoint(1, 1 << 199, 600),
                "upper": _endpoint(1, 1 << 199, 600),
            },
        )
        tiny_product = mul(point_huge, point_tiny)
        checks["extreme_residual_contains_zero"] = _contains(residual, [0, 0])
        checks["extreme_residual_contains_one"] = _contains(residual, [1, 0])
        checks["zero_plus_extreme_tiny_preserved"] = _exact(zero_tiny, [1, -10**18])
        checks["near_equal_opposite_contains_positive"] = _contains(near_cancel, [1, 600])
        checks["wide_mantissa_near_cancel_contains_one"] = _contains(wide_near_cancel, [1, 0])
        checks["extreme_product_contains_one"] = _contains(tiny_product, [1, 0])

    missing = _moments()
    del missing["H"]["z1"]["error_upper"]
    checks["undeclared_error_rejected"] = _raises(lambda: scale(missing, log_derivative))
    checks["float_rejected"] = _raises(lambda: scale(moments, 0.5))
    checks["reversed_interval_rejected"] = _raises(
        lambda: scale(moments, {"finite": [_positive(1), _negative(1)]})
    )
    cross = _moments("family-A")
    cross["C"]["source_family"] = "family-B"
    checks["source_mismatch_rejected"] = _raises(
        lambda: scale(cross, {"finite": [_negative(3), _negative(1)]})
    )
    checks["scope_flags"] = (
        result["assumptions"]["epsilon_and_F0_values_materialized"] is False
        and result["assumptions"]["source_identity_certified"] is False
        and result["assumptions"]["old_errors_cancelled"] is False
    )
    return {"schema": SCHEMA, "passed": all(checks.values()), "checks": checks}


__all__ = ["SCHEMA", "run_fixture"]
