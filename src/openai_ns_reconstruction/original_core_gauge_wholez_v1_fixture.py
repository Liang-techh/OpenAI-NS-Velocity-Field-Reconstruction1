"""Independent synthetic contract for the whole-Z gauge leaf.

This fixture intentionally does not import the production leaf or any older
project module.  It uses a different small rational source and proves a
nonzero positive integral on a simple cell by exact interval arithmetic.
The admitted synthetic Z domain is also only ``[-1,1]``; no real-axis
extension is tested.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any


ZERO = Fraction(0)
ONE = Fraction(1)


def _iv(value: Any) -> tuple[Fraction, Fraction]:
    if isinstance(value, tuple):
        if len(value) != 2:
            raise ValueError("synthetic interval tuple must have two endpoints")
        lo, hi = value
    else:
        lo = hi = value
    if not isinstance(lo, Fraction):
        lo = Fraction(lo)
    if not isinstance(hi, Fraction):
        hi = Fraction(hi)
    if lo > hi:
        raise ValueError("synthetic interval is reversed")
    return lo, hi


def _add(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]):
    return a[0] + b[0], a[1] + b[1]


def _mul(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]):
    values = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(values), max(values)


def _scale(a: tuple[Fraction, Fraction], value: Fraction):
    return _mul(a, (value, value))


def _square(a: tuple[Fraction, Fraction]):
    if a[0] <= 0 <= a[1]:
        return ZERO, max(a[0] * a[0], a[1] * a[1])
    return min(a[0] * a[0], a[1] * a[1]), max(a[0] * a[0], a[1] * a[1])


def _sub(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]):
    return a[0] - b[1], a[1] - b[0]


def _div_positive(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]):
    if b[0] <= 0:
        raise ArithmeticError("synthetic denominator positivity failed")
    return _mul(a, (ONE / b[1], ONE / b[0]))


def _h_point(v: Fraction, delta: Fraction, j: Fraction) -> Fraction:
    """Independent point formula; this does not call the production leaf."""

    return ((ONE - delta) * v) / 2 + (ONE - v * v) * (4 * v + j)


def _g_point(
    v: Fraction,
    delta: Fraction,
    j: Fraction,
    sigma: Fraction,
) -> Fraction:
    """Independent exact point value used only as an enclosure witness."""

    h = _h_point(v, delta, j)
    denominator = h * h + sigma * sigma
    if denominator <= 0:
        raise ArithmeticError("synthetic point denominator is not positive")
    return (ONE - delta * v * v) * h / denominator


def _scalar_record(value: Fraction) -> dict[str, int]:
    return {"numerator": value.numerator, "denominator": value.denominator}


def _record_scalar(record: dict[str, Any]) -> Fraction:
    return Fraction(record["numerator"], record["denominator"])


def _record_interval(record: dict[str, Any]) -> tuple[Fraction, Fraction]:
    return _record_scalar(record["lower"]), _record_scalar(record["upper"])


def _synthetic_parts(cell: tuple[Fraction, Fraction]):
    """Independent rectangle enclosure for derivative and denominator."""

    j = Fraction(1, 1000)
    delta = (ZERO, Fraction(1, 10**10))
    sigma = j / 500
    z = _iv(cell)
    one_minus_delta = (ONE - delta[1], ONE - delta[0])
    h = _add(
        _scale(_mul(one_minus_delta, z), Fraction(1, 2)),
        _mul(_sub((ONE, ONE), _square(z)),
             _add(_scale(z, Fraction(4)), (j, j))),
    )
    l = _sub((ONE, ONE), _mul(delta, _square(z)))
    denominator = _add(_square(h), (sigma * sigma, sigma * sigma))
    return _div_positive(_mul(l, h), denominator), denominator


def synthetic_contract() -> dict[str, Any]:
    """Return exact independent checks for one nonzero synthetic integral."""

    left = Fraction(1, 4)
    right = Fraction(1, 2)
    j = Fraction(1, 1000)
    sigma = j / 500
    derivative, denominator = _synthetic_parts((left, right))
    width = right - left
    integral = derivative[0] * width, derivative[1] * width
    point_z = Fraction(1, 4)
    point_delta_zero = _g_point(point_z, ZERO, j, sigma)
    positive_left = Fraction(1, 8)
    positive_right = Fraction(1, 4)
    positive_derivative, positive_denominator = _synthetic_parts(
        (positive_left, positive_right)
    )
    positive_integral_lower = positive_derivative[0] * (
        positive_right - positive_left
    )
    return {
        "schema": "original_core_gauge_wholez_v1_synthetic_fixture",
        "synthetic_inputs": {
            "j": {"numerator": 1, "denominator": 1000},
            "delta_interval": {
                "lower": {"numerator": 0, "denominator": 1},
                "upper": {"numerator": 1, "denominator": 10**10},
            },
            "sigma": {"numerator": 1, "denominator": 500000},
            "Z_cell": {
                "left": {"numerator": left.numerator, "denominator": left.denominator},
                "right": {"numerator": right.numerator, "denominator": right.denominator},
            },
        },
        "derivative_interval": {
            "lower": {"numerator": derivative[0].numerator,
                      "denominator": derivative[0].denominator},
            "upper": {"numerator": derivative[1].numerator,
                      "denominator": derivative[1].denominator},
        },
        "integral_interval": {
            "lower": {"numerator": integral[0].numerator,
                      "denominator": integral[0].denominator},
            "upper": {"numerator": integral[1].numerator,
                      "denominator": integral[1].denominator},
        },
        "checks": {
            "denominator_strictly_positive": denominator[0] > 0,
            "derivative_strictly_positive": derivative[0] > 0,
            "integral_nonzero_lower_bound": integral[0] > 0,
            "all_bounds_exact_fraction": True,
            "independent_of_production_leaf": True,
        },
        "provenance_declarations": {
            "source_owner_proved": False,
            "physical_m15": False,
            "closure": False,
        },
        "independent_reference": {
            "point_z": _scalar_record(point_z),
            "point_delta": _scalar_record(ZERO),
            "point_g_delta_zero": _scalar_record(point_delta_zero),
            "positive_integral_cell": {
                "left": _scalar_record(positive_left),
                "right": _scalar_record(positive_right),
            },
            "positive_integral_lower": _scalar_record(positive_integral_lower),
            "positive_integral_lower_is_strictly_positive": positive_integral_lower > 0,
            "positive_cell_denominator_lower": _scalar_record(positive_denominator[0]),
            "G_at_zero_nonnegative_used": True,
            "delta_upper_is_only_an_interval_bound": True,
        },
    }


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run independent references against the newly supplied production leaf.

    The caller passes the imported new leaf explicitly.  This fixture never
    imports that module (or any archived producer), so its point and positive
    integral witnesses remain independent arithmetic checks.
    """

    independent = synthetic_contract()
    j = Fraction(1, 1000)
    delta_upper = Fraction(1, 10**10)
    sigma = j / 500
    z_left = Fraction(1, 4)
    z_right = Fraction(1, 2)
    root_tolerance = Fraction(1, 1 << 80)
    geometry_levels = 24
    max_panels = 128

    source = leaf.GaugeSource(
        j=j,
        delta=leaf.RationalInterval(ZERO, delta_upper),
        sigma=sigma,
    )
    owner = leaf.ActualCoreGaugeWholeZ(
        source=source,
        root_tolerance=root_tolerance,
        geometry_levels=geometry_levels,
        max_panels=max_panels,
    )
    evaluation = owner.evaluate(z_left, z_right)

    # Independent delta = 0 point witness at the left edge of the cell.
    point_value = _g_point(z_left, ZERO, j, sigma)
    g_z_lower, g_z_upper = _record_interval(evaluation["G_Z_interval"])
    g_lower, g_upper = _record_interval(evaluation["G_interval"])

    positive_left = Fraction(1, 8)
    positive_right = Fraction(1, 4)
    positive_derivative, positive_denominator = _synthetic_parts(
        (positive_left, positive_right)
    )
    positive_integral_lower = positive_derivative[0] * (
        positive_right - positive_left
    )

    synthetic_anchor_derivative_lower = (
        (ONE - delta_upper) / 2 + 4 - 12 * j * j
    )
    synthetic_anchor_left_h = _h_point(-j, delta_upper, j)
    synthetic_anchor_zero_h = _h_point(ZERO, ZERO, j)
    g_at_zero_nonnegative = (
        synthetic_anchor_derivative_lower > ZERO
        and synthetic_anchor_left_h < ZERO
        and synthetic_anchor_zero_h > ZERO
    )

    source_record = evaluation["source"]
    source_delta_lower, source_delta_upper = _record_interval(
        source_record["delta_interval"]
    )
    anchor_record = evaluation["anchor"]
    anchor_delta_values = {
        _record_scalar(anchor_record["delta_lower_root"]["delta"]),
        _record_scalar(anchor_record["delta_upper_root"]["delta"]),
    }
    admitted_domain_record = evaluation["admitted_Z_domain"]
    z_record = evaluation["Z_cell"]
    error_budget = evaluation["error_budget"]
    error_intervals_are_ordered = all(
        _record_interval(error_budget[key])[0] == ZERO
        and _record_interval(error_budget[key])[1] >= ZERO
        for key in (
            "integration_truncation",
            "anchor_uncertainty",
            "rounding",
            "G_Z_rounding",
        )
    )

    try:
        owner.evaluate(Fraction(-2), ZERO)
    except ValueError:
        illegal_cell_rejected = True
    else:
        illegal_cell_rejected = False

    checks = {
        "leaf_called": True,
        "Z_cell_exactly_requested": (
            _record_scalar(z_record["lower"]) == z_left
            and _record_scalar(z_record["upper"]) == z_right
        ),
        "admitted_domain_is_exactly_minus_one_to_one": (
            _record_scalar(admitted_domain_record["lower"]) == -ONE
            and _record_scalar(admitted_domain_record["upper"]) == ONE
        ),
        "G_Z_contains_independent_delta_zero_point": (
            g_z_lower <= point_value <= g_z_upper
        ),
        "G_cell_bounds_nonnegative_and_ordered": (
            ZERO <= g_lower <= g_upper
        ),
        "G_cell_has_strictly_positive_lower_bound": g_lower > ZERO,
        "independent_positive_integral_lower_is_strict": (
            positive_denominator[0] > ZERO and positive_integral_lower > ZERO
        ),
        "G_upper_covers_positive_integral_from_[1/8,1/4]": (
            g_upper >= positive_integral_lower
        ),
        "G_at_zero_nonnegative_used": g_at_zero_nonnegative,
        "source_j_is_synthetic_exact": (
            _record_scalar(source_record["j"]) == j
        ),
        "source_sigma_is_exact_j_over_500": (
            _record_scalar(source_record["sigma"]) == sigma
        ),
        "delta_parameter_enclosure_is_exact": (
            source_delta_lower == ZERO and source_delta_upper == delta_upper
        ),
        "delta_lower_not_selected": (
            source_record["delta_lower_is_not_selected_value"] is True
        ),
        "delta_upper_not_selected": (
            source_record["delta_upper_is_not_selected_value"] is True
        ),
        "anchor_parameter_endpoints_are_only_bounds": (
            anchor_delta_values == {ZERO, delta_upper}
            and anchor_record["root_center_not_selected"] is True
        ),
        "no_delta_or_anchor_center_selected": (
            evaluation["delta_center_selected"] is False
            and evaluation["anchor_center_selected"] is False
            and evaluation["midpoint_sample_used"] is False
        ),
        "error_fields_are_separate_and_fractional": (
            error_budget["all_endpoints_exact_fraction"] is True
            and "integration_truncation" in error_budget
            and "anchor_uncertainty" in error_budget
            and "rounding" in error_budget
            and "G_Z_rounding" in error_budget
            and error_intervals_are_ordered
        ),
        "illegal_cell_rejected": illegal_cell_rejected,
    }

    return {
        "schema": "original_core_gauge_wholez_v1_fixture_leaf_call",
        "synthetic_only": True,
        "passed": all(checks.values()),
        "independent_reference": independent,
        "leaf_call": {
            "root_tolerance": _scalar_record(root_tolerance),
            "geometry_levels": geometry_levels,
            "max_panels": max_panels,
            "Z_cell": {
                "left": _scalar_record(z_left),
                "right": _scalar_record(z_right),
            },
            "admitted_Z_domain": {
                "lower": _scalar_record(-ONE),
                "upper": _scalar_record(ONE),
            },
            "source_delta_upper_is_bound_not_value": True,
        },
        "evaluation": evaluation,
        "independent_witnesses": {
            "delta_zero_g_at_Z_left": _scalar_record(point_value),
            "positive_integral_cell": {
                "left": _scalar_record(positive_left),
                "right": _scalar_record(positive_right),
            },
            "positive_integral_lower": _scalar_record(positive_integral_lower),
            "G_at_zero_nonnegative": g_at_zero_nonnegative,
            "anchor_sign_evidence": {
                "H_at_minus_j": _scalar_record(synthetic_anchor_left_h),
                "H_at_zero": _scalar_record(synthetic_anchor_zero_h),
                "H_prime_lower": _scalar_record(synthetic_anchor_derivative_lower),
            },
        },
        "checks": checks,
    }


__all__ = ["synthetic_contract", "run_fixture"]
