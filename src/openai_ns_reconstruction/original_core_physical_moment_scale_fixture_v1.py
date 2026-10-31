"""Synthetic checks for the coefficient-only physical scale ledger."""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_core_physical_moment_scale_fixture_v1"


def _endpoints(value: Any) -> tuple[Any, Any]:
    if hasattr(value, "lower") and hasattr(value, "upper"):
        return value.lower, value.upper
    if isinstance(value, (tuple, list)) and len(value) == 2:
        return value[0], value[1]
    if isinstance(value, dict):
        return value["lower"], value["upper"]
    raise TypeError("fixture cannot read interval endpoints")


def _exact(value: Any, target: Fraction) -> bool:
    lower, upper = _endpoints(value)
    return lower == target == upper


def _contains(value: Any, target: Fraction) -> bool:
    lower, upper = _endpoints(value)
    return lower <= target <= upper


def _raises(call: Callable[[], Any]) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _moments(source_family: str = "fixture-family") -> dict[str, Any]:
    values = {
        "H": (1, 2),
        "M": (2, 3),
        "K": (3, 4),
        "A": (5, 6),
        "B": (7, 8),
        "C": (11, 12),
    }
    return {
        name: {
            "z0": {"finite": [value, value], "error_upper": 0},
            "z1": {"finite": [derivative, derivative], "error_upper": 0},
            "source_family": source_family,
        }
        for name, (value, derivative) in values.items()
    }


def run_fixture(scale: Callable[..., Any]) -> dict[str, Any]:
    """Run hand values, uncertainty, interval, and rejection checks."""

    moments = _moments()
    log_derivative = {"finite": [-3, -1], "source_family": "fixture-family"}
    result = scale(moments, log_derivative, source_family="fixture-family")
    expected_finite = {
        "I_V": Fraction(8),
        "I_V2": Fraction(20),
        "I_RF": Fraction(8),
        "I_RFV": Fraction(24),
        "I_RF2": Fraction(112),
        "I_F2": Fraction(44),
    }
    expected_z1 = {
        "I_V": Fraction(12),
        "I_V2": Fraction(24),
        "I_RF": (-8, 8),
        "I_RFV": (-40, 8),
        "I_RF2": (-544, -96),
        "I_F2": (-216, -40),
    }
    checks: dict[str, bool] = {}
    checks["schema"] = result.get("schema") == SCHEMA.replace("_fixture", "")
    checks["finite_hand_values"] = all(
        _exact(result["finite"][name], value) for name, value in expected_finite.items()
    )
    checks["ordinary_Z1_hand_values"] = all(
        (
            _exact(result["finite_Z1"][name], value)
            if isinstance(value, Fraction)
            else _endpoints(result["finite_Z1"][name]) == tuple(Fraction(item) for item in value)
        )
        for name, value in expected_z1.items()
    )
    checks["scale_powers_separate"] = all(
        result["scale"][name] == powers
        for name, powers in {
            "I_V": (1, 0),
            "I_V2": (1, 0),
            "I_RF": (2, 1),
            "I_RFV": (2, 1),
            "I_RF2": (2, 2),
            "I_F2": (1, 2),
        }.items()
    )
    checks["negative_nontrivial_L_retained"] = (
        result["L"].lower == Fraction(-3)
        and result["L"].upper == Fraction(-1)
        and result["finite_Z1"]["I_RF"].lower < 0 < result["finite_Z1"]["I_RF"].upper
    )
    checks["zero_errors_by_default"] = all(
        result["error_upper"][name] == 0 for name in expected_finite
    )
    checks["total_is_widened_coefficient_row"] = all(
        result["total"][name] == result["enclosure"][name] for name in expected_finite
    )

    noisy_moments = _moments()
    for row in noisy_moments.values():
        row["z0"]["error_upper"] = "1/10"
        row["z1"]["error_upper"] = "1/20"
    noisy = scale(noisy_moments, log_derivative, source_family="fixture-family")
    checks["nonzero_errors_add"] = (
        noisy["error_upper"]["I_V"] == Fraction(2, 5)
        and noisy["error_upper_Z1"]["I_V"] == Fraction(1, 5)
        and noisy["error_upper_Z1"]["I_RF"] == Fraction(14, 5)
        and noisy["enclosure"]["I_RF"].lower < noisy["finite"]["I_RF"].lower
    )

    interval_moments = _moments()
    interval_moments["H"]["z0"]["finite"] = ["9/10", "11/10"]
    interval_moments["H"]["z1"]["finite"] = ["19/10", "21/10"]
    interval_result = scale(interval_moments, log_derivative, source_family="fixture-family")
    checks["input_intervals_preserved"] = (
        _contains(interval_result["finite"]["I_RF"], Fraction(36, 5))
        and _contains(interval_result["finite_Z1"]["I_RF"], Fraction(-8))
    )

    checks["float_rejected"] = _raises(
        lambda: scale(moments, 0.5)
    )
    checks["negative_error_rejected"] = _raises(
        lambda: scale(
            {**moments, "H": {**moments["H"], "z0": {"finite": [1, 1], "error_upper": -1}}},
            log_derivative,
        )
    )
    missing_error = _moments()
    del missing_error["H"]["z0"]["error_upper"]
    checks["missing_error_rejected"] = _raises(
        lambda: scale(missing_error, log_derivative)
    )
    checks["reversed_interval_rejected"] = _raises(
        lambda: scale(moments, [1, -1])
    )
    cross = _moments("family-A")
    cross["C"]["source_family"] = "family-B"
    checks["cross_source_rejected"] = _raises(
        lambda: scale(cross, log_derivative)
    )
    checks["caller_scope_flags"] = (
        result["assumptions"]["L_is_caller_supplied_interval"] is True
        and result["assumptions"]["source_identity_certified"] is False
        and result["assumptions"]["closure_proved"] is False
        and result["assumptions"]["old_errors_cancelled"] is False
        and result["assumptions"]["physical_scale_materialized"] is False
    )
    return {"schema": SCHEMA, "passed": all(checks.values()), "checks": checks}


__all__ = ["SCHEMA", "run_fixture"]
