"""Synthetic-only checks for :mod:`n1_sigma_interval_rational_v1`.

The fixture receives the candidate module from its caller.  It constructs no
graph, reads no receipt, imports no prior sigma family, and uses no floating
point reference.  The checks exercise exact endpoints, reflection, cells
crossing the support boundaries, the analytic midpoint value, proved
rational exponential bounds, and fail-closed input validation.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Callable


SCHEMA = "n1_sigma_interval_rational_v1_fixture"


def _rejected(call: Callable[[], Any]) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _all_false(checks: dict[str, bool], names: tuple[str, ...]) -> None:
    for name in names:
        checks[name] = False


def run_fixture(module: Any) -> dict[str, Any]:
    """Run independent synthetic interval checks against ``module``."""

    names = (
        "api",
        "schema",
        "rational_endpoint_types",
        "negative_endpoint",
        "zero_endpoint",
        "one_endpoint",
        "positive_endpoint",
        "half_exact",
        "cross_zero",
        "cross_half",
        "cross_one",
        "whole_support",
        "reflection_sigma",
        "reflection_derivative",
        "quarter_analytic_order",
        "independent_exp_one_bracket",
        "independent_quarter_bracket",
        "endpoint_derivative_cap",
        "derivative_nonnegative",
        "exp_zero_exact",
        "exp_large_rational_bound",
        "ordinary_x_metadata",
        "source_identity_unverified",
        "no_owner_binding",
        "no_legacy_module",
        "wrong_type_rejected",
        "reversed_interval_rejected",
        "bad_control_rejected",
        "bad_exp_input_rejected",
        "bad_interval_rejected",
        "interval_resource_rejected",
        "resource_input_bits_rejected",
        "resource_terms_rejected",
        "resource_reduction_rejected",
        "resource_power_rejected",
    )
    checks: dict[str, bool] = {name: False for name in names}
    evaluator = getattr(module, "evaluate_sigma_interval", None)
    exp_evaluator = getattr(module, "rational_exp_negative", None)
    interval_type = getattr(module, "RationalInterval", None)
    checks["api"] = (
        callable(evaluator) and callable(exp_evaluator) and interval_type is not None
    )
    selected: dict[str, Any] = {}

    if checks["api"]:
        try:
            F = Fraction
            negative = evaluator(F(-2), F(-1))
            zero = evaluator(F(0), F(0))
            one = evaluator(F(1), F(1))
            positive = evaluator(F(2), F(3))
            half = evaluator(F(1, 2), F(1, 2))
            cross_zero = evaluator(F(-1, 2), F(1, 4))
            cross_half = evaluator(F(1, 4), F(3, 4))
            cross_one = evaluator(F(3, 4), F(3, 2))
            whole = evaluator(F(-1), F(2))
            selected = {
                "schema": getattr(module, "SCHEMA", None),
                "component": getattr(module, "COMPONENT", None),
                "whole_sigma": (
                    whole.sigma.lower,
                    whole.sigma.upper,
                ),
                "whole_derivative_upper": whole.d_sigma_dx.upper,
            }

            checks["schema"] = (
                selected["schema"] == "n1_sigma_interval_rational_v1"
                and selected["component"] == "original_flat_sigma_interval_rational"
            )
            checks["rational_endpoint_types"] = all(
                isinstance(value, Fraction)
                for result in (negative, zero, one, positive, half, cross_half)
                for interval in (result.sigma, result.d_sigma_dx)
                for value in (interval.lower, interval.upper)
            )
            checks["negative_endpoint"] = (
                negative.sigma.lower == F(0)
                and negative.sigma.upper == F(0)
                and negative.d_sigma_dx.lower == F(0)
                and negative.d_sigma_dx.upper == F(0)
            )
            checks["zero_endpoint"] = (
                zero.sigma.lower == F(0)
                and zero.sigma.upper == F(0)
                and zero.d_sigma_dx.lower == F(0)
                and zero.d_sigma_dx.upper == F(0)
            )
            checks["one_endpoint"] = (
                one.sigma.lower == F(1)
                and one.sigma.upper == F(1)
                and one.d_sigma_dx.lower == F(0)
                and one.d_sigma_dx.upper == F(0)
            )
            checks["positive_endpoint"] = (
                positive.sigma.lower == F(1)
                and positive.sigma.upper == F(1)
                and positive.d_sigma_dx.lower == F(0)
                and positive.d_sigma_dx.upper == F(0)
            )
            checks["half_exact"] = (
                half.sigma.lower == F(1, 2)
                and half.sigma.upper == F(1, 2)
                and half.d_sigma_dx.lower == F(8)
                and half.d_sigma_dx.upper == F(8)
            )
            checks["cross_zero"] = (
                cross_zero.sigma.lower == F(0)
                and F(0) < cross_zero.sigma.upper < F(1)
                and cross_zero.d_sigma_dx.lower == F(0)
                and cross_zero.d_sigma_dx.upper > F(0)
            )
            checks["cross_half"] = (
                F(0) < cross_half.sigma.lower < cross_half.sigma.upper < F(1)
                and cross_half.d_sigma_dx.lower == F(0)
                and cross_half.d_sigma_dx.upper >= F(8)
            )
            checks["cross_one"] = (
                F(0) < cross_one.sigma.lower < F(1)
                and cross_one.sigma.upper == F(1)
                and cross_one.d_sigma_dx.lower == F(0)
                and cross_one.d_sigma_dx.upper > F(0)
            )
            checks["whole_support"] = (
                whole.sigma == interval_type(F(0), F(1))
                and whole.d_sigma_dx.lower == F(0)
                and whole.d_sigma_dx.upper > F(0)
            )

            left = evaluator(F(1, 7), F(2, 7))
            right = evaluator(F(5, 7), F(6, 7))
            checks["reflection_sigma"] = right.sigma == interval_type(
                F(1) - left.sigma.upper,
                F(1) - left.sigma.lower,
            )
            checks["reflection_derivative"] = (
                right.d_sigma_dx == left.d_sigma_dx
            )

            quarter = evaluator(F(1, 4), F(1, 4))
            checks["quarter_analytic_order"] = (
                F(0) < quarter.sigma.lower
                and quarter.sigma.upper < F(1, 2)
                and quarter.d_sigma_dx.lower >= F(0)
                and quarter.d_sigma_dx.upper > F(0)
            )

            # Independent Bernoulli bounds, not the leaf's Taylor/range
            # reduction, prove the comparison brackets.  For t > 0 and n >
            # t, (1-t/n)^n <= exp(-t) <= (1+t/n)^(-n).
            exp_one_lower = (F(255, 256)) ** 256
            exp_one_upper = (F(256, 257)) ** 256
            exp_one = exp_evaluator(F(1))
            checks["independent_exp_one_bracket"] = (
                F(367, 1000) < exp_one_lower
                and exp_one_upper < F(369, 1000)
                and exp_one_lower <= exp_one.lower <= exp_one.upper <= exp_one_upper
            )
            t_quarter = F(128, 9)
            q_quarter_lower = (F(1) - t_quarter / F(256)) ** 256
            q_quarter_upper = (F(1) + t_quarter / F(256)) ** (-256)
            independent_quarter_lower = q_quarter_lower / (F(1) + q_quarter_lower)
            independent_quarter_upper = q_quarter_upper / (F(1) + q_quarter_upper)
            checks["independent_quarter_bracket"] = (
                independent_quarter_lower <= quarter.sigma.lower
                <= quarter.sigma.upper <= independent_quarter_upper
            )
            endpoint_strip = evaluator(F(0), F(1, 4))
            endpoint_cap = getattr(
                module,
                "ENDPOINT_DERIVATIVE_CAP",
                getattr(module, "_ENDPOINT_DERIVATIVE_CAP", None),
            )
            checks["endpoint_derivative_cap"] = (
                isinstance(endpoint_cap, Fraction)
                and endpoint_cap > F(0)
                and endpoint_strip.d_sigma_dx.upper == endpoint_cap
            )
            checks["derivative_nonnegative"] = all(
                result.d_sigma_dx.lower >= F(0)
                for result in (negative, zero, one, positive, half, cross_zero, cross_half, cross_one, whole)
            )

            exp_zero = exp_evaluator(F(0))
            exp_large = exp_evaluator(F(1000000))
            checks["exp_zero_exact"] = exp_zero == interval_type(F(1), F(1))
            checks["exp_large_rational_bound"] = (
                F(0) <= exp_large.lower <= exp_large.upper <= F(1)
                and exp_large.upper < F(1, 1000000)
            )

            checks["ordinary_x_metadata"] = (
                cross_half.certified is True
                and cross_half.ordinary_x_derivative is True
                and cross_half.z_derivative is False
                and cross_half.derivative == cross_half.d_sigma_dx
                and cross_half.algorithm_family == "n1_sigma_interval_rational_v1"
            )
            checks["source_identity_unverified"] = (
                cross_half.source_identity_verified is False
            )
            checks["no_owner_binding"] = cross_half.owner_binding is False
            checks["no_legacy_module"] = cross_half.used_legacy_modules is False

            checks["wrong_type_rejected"] = _rejected(
                lambda: evaluator(0, F(1))
            ) and _rejected(lambda: evaluator(F(0), 1))
            checks["reversed_interval_rejected"] = _rejected(
                lambda: evaluator(F(1), F(0))
            )
            checks["bad_control_rejected"] = (
                _rejected(lambda: evaluator(F(0), F(1), series_terms=True))
                and _rejected(lambda: evaluator(F(0), F(1), max_reduction_steps=0))
                and _rejected(lambda: evaluator(F(0), F(1), fallback_power=0))
            )
            checks["bad_exp_input_rejected"] = (
                _rejected(lambda: exp_evaluator(F(-1)))
                and _rejected(lambda: exp_evaluator(1))
            )
            checks["bad_interval_rejected"] = _rejected(
                lambda: interval_type(F(1), F(0))
            )
            intermediate_limit = getattr(module, "MAX_INTERMEDIATE_NUMERATOR_BITS")
            checks["interval_resource_rejected"] = _rejected(
                lambda: interval_type(
                    F(1 << intermediate_limit), F(1 << intermediate_limit)
                )
            )

            input_limit = getattr(module, "MAX_INPUT_NUMERATOR_BITS")
            denominator_limit = getattr(module, "MAX_INPUT_DENOMINATOR_BITS")
            checks["resource_input_bits_rejected"] = _rejected(
                lambda: evaluator(F(1 << input_limit), F(1 << input_limit))
            ) and _rejected(
                lambda: evaluator(F(1, 1 << denominator_limit), F(1))
            )
            terms_limit = getattr(module, "MAX_SERIES_TERMS")
            reduction_limit = getattr(module, "MAX_REDUCTION_STEPS")
            power_limit = getattr(module, "MAX_FALLBACK_POWER")
            checks["resource_terms_rejected"] = _rejected(
                lambda: evaluator(F(0), F(1), series_terms=terms_limit + 1)
            )
            checks["resource_reduction_rejected"] = _rejected(
                lambda: evaluator(
                    F(0), F(1), max_reduction_steps=reduction_limit + 1
                )
            )
            checks["resource_power_rejected"] = _rejected(
                lambda: evaluator(F(0), F(1), fallback_power=power_limit + 1)
            )
        except Exception:
            _all_false(checks, names[1:])

    return {
        "schema": SCHEMA,
        "component": "original_flat_sigma_interval_rational",
        "passed": all(checks.values()),
        "checks": checks,
        "synthetic_only": True,
        "used_old_results": False,
        "used_old_module": False,
        "owner_binding": False,
        "selected": selected,
    }


run = run_fixture


__all__ = ["SCHEMA", "run", "run_fixture"]
