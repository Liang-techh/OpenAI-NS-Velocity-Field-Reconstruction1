"""Independent small-rational checks for the variable-Z reserved leaf.

The fixture deliberately supplies nonzero values for all four unknown source
derivatives.  It checks the leaf's baseline chain rule and each coefficient
from the direct product rule; it does not import an ancestor producer.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields, is_dataclass
from typing import Any


def _require(condition: object, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def _field(value: Any, name: str) -> Any:
    if isinstance(value, Mapping):
        _require(name in value, "leaf result is missing " + name)
        return value[name]
    if is_dataclass(value):
        names = {item.name for item in fields(value)}
        _require(name in names, "leaf dataclass is missing " + name)
        return getattr(value, name)
    _require(False, "leaf result is not a dataclass or mapping")
    raise AssertionError(name)


def _singleton(c: Any, endpoint: Any) -> Any:
    return c.make_mpf((tuple(endpoint), tuple(endpoint)))


def _intersects(c: Any, left: Any, right: Any) -> bool:
    """Use numerical singleton comparisons, never lexicographic mpf tuples."""

    _require(hasattr(left, "_mpi_") and hasattr(right, "_mpi_"), "expected interval result")
    left_lower, left_upper = left._mpi_
    right_lower, right_upper = right._mpi_
    left_lower_value = _singleton(c, left_lower)
    left_upper_value = _singleton(c, left_upper)
    right_lower_value = _singleton(c, right_lower)
    right_upper_value = _singleton(c, right_upper)
    return bool(not (left_upper_value < right_lower_value or right_upper_value < left_lower_value))


def _excludes_zero(c: Any, value: Any) -> bool:
    _require(hasattr(value, "_mpi_"), "expected interval result")
    lower, upper = value._mpi_
    lower_value = _singleton(c, lower)
    upper_value = _singleton(c, upper)
    zero = c.mpf(0)
    return bool(upper_value < zero or zero < lower_value)


def _sensitivity_mapping(value: Any) -> dict[str, Any]:
    _require(isinstance(value, (list, tuple)), "sensitivity_coefficients is not an ordered pair tuple")
    pairs: dict[str, Any] = {}
    for index, pair in enumerate(value):
        _require(isinstance(pair, (list, tuple)) and len(pair) == 2, f"sensitivity pair {index} is invalid")
        name, coefficient = pair
        _require(isinstance(name, str), f"sensitivity pair {index} name is not text")
        _require(name not in pairs, "duplicate sensitivity name: " + name)
        pairs[name] = coefficient
    _require(set(pairs) == {"r0_Z", "r1_Z", "r2_Z", "mu_Z"}, "sensitivity names are incomplete")
    return pairs


def run_fixture(leaf: Any, c: Any) -> dict[str, Any]:
    """Run one independent rational probe through ``leaf``.

    All supplied inputs are small integers (hence rational values).  ``mu``
    and ``Z`` are zero so the probe isolates the ordinary Z1 product rule while
    the four derivative inputs remain deliberately nonzero.
    """

    evaluate = getattr(leaf, "evaluate_partial_reserved_minus_pulse", None)
    _require(callable(evaluate), "leaf evaluate_partial_reserved_minus_pulse is missing")

    def iv(number: int) -> Any:
        return c.mpf(number)

    r_i = (iv(2), iv(3), iv(5))
    r_i_Z = (iv(7), iv(11), iv(13))
    mu = iv(0)
    mu_Z = iv(17)
    theta_i = (iv(19), iv(23))
    theta_i_Z = (iv(29), iv(31))
    J = iv(37)
    J_Z = iv(41)
    P = iv(43)
    P_Z = iv(47)
    Z = iv(0)

    # This is an independent re-expression of the declared direct formula.
    exponent = -c.mpf(17) / 10 - c.mpf(5) * mu / 2
    amplitude = c.exp(exponent) / (1 + Z**2)
    amplitude_Z = -2 * Z * c.exp(exponent) / (1 + Z**2) ** 2
    own = -r_i[0] * amplitude * theta_i[0] - r_i[1] * amplitude * theta_i[1]
    c0_expected = own + r_i[2] * J + (1 - r_i[2]) * P
    baseline_expected = (
        -r_i[0] * amplitude_Z * theta_i[0]
        - r_i[1] * amplitude_Z * theta_i[1]
        - r_i[0] * amplitude * theta_i_Z[0]
        - r_i[1] * amplitude * theta_i_Z[1]
        + r_i[2] * J_Z
        + (1 - r_i[2]) * P_Z
    )
    sensitivity_expected = {
        "r0_Z": -amplitude * theta_i[0],
        "r1_Z": -amplitude * theta_i[1],
        "r2_Z": J - P,
        # The coefficient is -(5/2)*own, not -(5/2)*amplitude.
        "mu_Z": -(c.mpf(5) / 2) * own,
    }
    z1_expected = baseline_expected + sum(
        (sensitivity_expected[name] * value)
        for name, value in (
            ("r0_Z", r_i_Z[0]),
            ("r1_Z", r_i_Z[1]),
            ("r2_Z", r_i_Z[2]),
            ("mu_Z", mu_Z),
        )
    )

    result = evaluate(
        c,
        r_i,
        r_i_Z,
        mu,
        mu_Z,
        theta_i,
        theta_i_Z,
        J,
        J_Z,
        P,
        P_Z,
        Z,
    )
    c0 = _field(result, "C0")
    z1_baseline = _field(result, "Z1_baseline")
    z1 = _field(result, "Z1")
    coefficients = _sensitivity_mapping(_field(result, "sensitivity_coefficients"))
    c0_terms = _field(result, "C0_terms")
    _require(c0_terms is not None, "C0_terms is missing")

    checks = {
        "C0_intersects_manual_rational_formula": _intersects(c, c0, c0_expected),
        "Z1_baseline_intersects_manual_chain_rule": _intersects(c, z1_baseline, baseline_expected),
        "Z1_intersects_nonzero_derivative_chain_rule": _intersects(c, z1, z1_expected),
        "r0_Z_coefficient_intersects": _intersects(c, coefficients["r0_Z"], sensitivity_expected["r0_Z"]),
        "r1_Z_coefficient_intersects": _intersects(c, coefficients["r1_Z"], sensitivity_expected["r1_Z"]),
        "r2_Z_coefficient_intersects": _intersects(c, coefficients["r2_Z"], sensitivity_expected["r2_Z"]),
        "mu_Z_coefficient_uses_own": _intersects(c, coefficients["mu_Z"], sensitivity_expected["mu_Z"]),
        "r0_Z_coefficient_nonzero": _excludes_zero(c, coefficients["r0_Z"]),
        "r1_Z_coefficient_nonzero": _excludes_zero(c, coefficients["r1_Z"]),
        "r2_Z_coefficient_nonzero": _excludes_zero(c, coefficients["r2_Z"]),
        "mu_Z_coefficient_nonzero": _excludes_zero(c, coefficients["mu_Z"]),
        "all_derivative_inputs_nonzero": all(
            _excludes_zero(c, value) for value in (*r_i_Z, mu_Z)
        ),
        "C0_terms_present": True,
    }
    _require(all(checks.values()), "independent variable-Z fixture check failed")
    return {
        "schema": "reserved_m15_variable_z_fixture_v1",
        "passed": True,
        "check_count": len(checks),
        "checks": checks,
        "probe": {
            "r_i": r_i,
            "r_i_Z": r_i_Z,
            "mu": mu,
            "mu_Z": mu_Z,
            "theta_i": theta_i,
            "theta_i_Z": theta_i_Z,
            "J": J,
            "J_Z": J_Z,
            "P": P,
            "P_Z": P_Z,
            "Z": Z,
        },
        "expected": {
            "C0": c0_expected,
            "Z1_baseline": baseline_expected,
            "Z1": z1_expected,
            "sensitivity_coefficients": sensitivity_expected,
        },
        "derivative_inputs_are_nonzero": True,
        "mu_Z_coefficient_formula": "-(5/2)*own",
        "comparison_mode": "interval_intersection_only",
        "strict_error_proof": False,
        "P1_consumed": False,
        "physical_Rp_used": False,
    }
