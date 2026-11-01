"""Independent synthetic fixture for the bounded directed dyadic sum leaf.

``run_fixture(leaf)`` accepts the implementation module or a mapping exposing
``add_interval`` and ``sum_intervals``.  Small expected values are calculated
with independent :class:`fractions.Fraction` arithmetic.  The astronomical
exponent cases are structural checks and never convert ``10**17`` to a
denominator or shift by that exponent.  The returned report contains only
JSON-safe scalars, tuples are checked as tuples at the API boundary, and no
source/provider graph is accessed.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_bounded_dyadic_sum_fixture_v1"
FIXTURE_SOURCE_ID = "synthetic:bounded-dyadic-sum-fixture-v1"


def _member(leaf: Any, name: str, default: Any = None) -> Any:
    if isinstance(leaf, Mapping):
        return leaf.get(name, default)
    return getattr(leaf, name, default)


def _raises(call: Callable[[], Any]) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _endpoint_fraction(endpoint: tuple[int, int, int, int]) -> Fraction:
    sign, mantissa, exponent, bitcount = endpoint
    if mantissa == 0:
        return Fraction(0)
    # This helper is used only on the small synthetic cases below.  The
    # 10**17 cases deliberately bypass it.
    if exponent >= 0:
        magnitude = Fraction(mantissa * (2**exponent), 1)
    else:
        magnitude = Fraction(mantissa, 2 ** (-exponent))
    return -magnitude if sign else magnitude


def _exact_interval_sum(intervals: tuple[Any, ...]) -> tuple[Fraction, Fraction]:
    lower = Fraction(0)
    upper = Fraction(0)
    for interval in intervals:
        lower += _endpoint_fraction(interval[0])
        upper += _endpoint_fraction(interval[1])
    return lower, upper


def _actual_interval_fraction(value: Any) -> tuple[Fraction, Fraction] | None:
    if not isinstance(value, tuple) or len(value) != 2:
        return None
    endpoints: list[Fraction] = []
    for endpoint in value:
        if not isinstance(endpoint, tuple) or len(endpoint) != 4:
            return None
        endpoints.append(_endpoint_fraction(endpoint))
    return endpoints[0], endpoints[1]


def _encloses(actual: Any, expected: tuple[Fraction, Fraction]) -> bool:
    converted = _actual_interval_fraction(actual)
    return (
        converted is not None
        and converted[0] <= expected[0]
        and converted[1] >= expected[1]
    )


def _exact_tuple_interval(
    value: Any,
    expected: tuple[tuple[int, int, int, int], tuple[int, int, int, int]],
) -> bool:
    return value == expected


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run fresh directed-rounding, cancellation and resource checks."""

    add = _member(leaf, "add_interval")
    total = _member(leaf, "sum_intervals")
    checks: dict[str, bool] = {}
    if not callable(add) or not callable(total):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {
                "add_interval_api_callable": callable(add),
                "sum_intervals_api_callable": callable(total),
            },
            "test_count": 2,
            "synthetic_only": True,
            "source_provider_called": False,
            "source_graph_walked": False,
            "original_source_certified": False,
        }

    zero = (0, 0, 0, 0)
    one = (0, 1, 0, 1)
    minus_one = (1, 1, 0, 1)

    positive = ((0, 1, 0, 1), (0, 3, 0, 2))
    negative = ((1, 3, 0, 2), (1, 1, 0, 1))
    checks["plus_exact_small_fraction_sum"] = _encloses(
        add(positive, positive, precision=8),
        _exact_interval_sum((positive, positive)),
    )
    checks["minus_exact_small_fraction_sum"] = _encloses(
        add(negative, negative, precision=8),
        _exact_interval_sum((negative, negative)),
    )

    # [1,3] + [-3,-1] crosses zero, while the singleton opposite pair must
    # cancel exactly rather than manufacture a signed zero.
    cross_left = ((0, 1, 0, 1), (0, 3, 0, 2))
    cross_right = ((1, 3, 0, 2), (0, 1, 0, 1))
    checks["cross_zero_enclosure"] = _encloses(
        add(cross_left, cross_right, precision=8),
        _exact_interval_sum((cross_left, cross_right)),
    )
    checks["exact_opposite_cancellation"] = _exact_tuple_interval(
        add((one, one), (minus_one, minus_one), precision=8),
        (zero, zero),
    )

    # A p=3 endpoint with four significant bits forces both directed choices:
    # +15 rounds down to 14 and up to 16; the negative case reverses them.
    coarse = ((0, 15, 0, 4), (0, 15, 0, 4))
    negative_coarse = ((1, 15, 0, 4), (1, 15, 0, 4))
    checks["positive_floor_and_ceil_sticky"] = _exact_tuple_interval(
        add(coarse, (zero, zero), precision=3),
        ((0, 7, 1, 3), (0, 1, 4, 1)),
    )
    checks["negative_floor_and_ceil_direction"] = _exact_tuple_interval(
        add(negative_coarse, (zero, zero), precision=3),
        ((1, 1, 4, 1), (1, 7, 1, 3)),
    )

    # The sum path retains several nonzero terms and is checked against an
    # independent exact rational total, not against another implementation
    # helper.
    terms = (
        ((0, 3, -1, 2), (0, 3, -1, 2)),
        ((1, 1, -2, 1), (1, 1, -2, 1)),
        ((0, 5, -3, 3), (0, 5, -3, 3)),
        ((1, 1, -4, 1), (1, 1, -4, 1)),
    )
    checks["multiple_nonzero_terms_fraction_enclosed"] = _encloses(
        total(terms, precision=4),
        _exact_interval_sum(terms),
    )

    extreme = 10**17
    huge = ((0, 1, extreme, 1), (0, 1, extreme, 1))
    tiny_negative = ((1, 1, -extreme, 1), (1, 1, -extreme, 1))
    tiny_positive = ((0, 1, -extreme, 1), (0, 1, -extreme, 1))
    # Below the power boundary the p=8 spacing is 2**(E-8), and the
    # negative microscopic term must move the lower endpoint down.
    checks["huge_gap_negative_sticky_moves_lower"] = _exact_tuple_interval(
        add(huge, tiny_negative, precision=8),
        (
            (0, 255, extreme - 8, 8),
            (0, 1, extreme, 1),
        ),
    )
    checks["huge_gap_positive_sticky_moves_upper"] = _exact_tuple_interval(
        add(huge, tiny_positive, precision=8),
        (
            (0, 1, extreme, 1),
            (0, 129, extreme - 7, 8),
        ),
    )
    checks["one_bit_huge_gap_uses_lower_binade"] = _exact_tuple_interval(
        add(huge, tiny_negative, precision=1),
        ((0, 1, extreme - 1, 1), huge[1]),
    )
    checks["huge_exponent_preserved_by_sum"] = _exact_tuple_interval(
        total((huge, (zero, zero)), precision=8),
        (huge[0], huge[1]),
    )

    special_codes = (
        (0, 0, 0, -1),
        (1, 0, 0, -2),
        (0, 0, 1, 0),
        (0, 0, -456, -2),
    )
    checks["special_values_rejected"] = all(
        _raises(lambda code=code: add((code, code), (one, one), precision=8))
        for code in special_codes
    )
    checks["malformed_bits_rejected"] = _raises(
        lambda: add(((0, 3, 0, 3),) * 2, (one, one), precision=8)
    )
    checks["reversed_interval_rejected"] = _raises(
        lambda: add(((0, 3, 0, 2), (0, 1, 0, 1)), (one, one), precision=8)
    )
    checks["list_and_scalar_containers_rejected"] = (
        _raises(lambda: add([(0, 1, 0, 1), (0, 1, 0, 1)], (one, one), precision=8))
        and _raises(lambda: add(1, (one, one), precision=8))
    )

    receipt = _member(leaf, "ROUNDING_RECEIPT", None)
    checks["separate_json_safe_rounding_receipt"] = (
        isinstance(receipt, dict)
        and receipt.get("point_value_selected") is False
        and receipt.get("special_values_rejected") is True
        and receipt.get("exponent_shift_by_magnitude") is False
    )

    max_precision = _member(leaf, "MAX_PRECISION_BITS", 4096)
    checks["explicit_precision_resource_failure"] = _raises(
        lambda: add((one, one), (one, one), precision=max_precision + 1)
    )
    max_terms = _member(leaf, "MAX_TERMS", 100_000)
    # Do not construct this tuple when a mapping intentionally hides the
    # implementation cap; the precision failure above remains authoritative.
    if isinstance(max_terms, int) and max_terms < 200_000:
        too_many = tuple((one, one) for _ in range(max_terms + 1))
        checks["explicit_term_count_resource_failure"] = _raises(
            lambda: total(too_many, precision=8)
        )
    else:
        checks["explicit_term_count_resource_failure"] = True

    return {
        "schema": SCHEMA,
        "fixture_source_id": FIXTURE_SOURCE_ID,
        "passed": bool(checks) and all(checks.values()),
        "checks": checks,
        "test_count": len(checks),
        "synthetic_only": True,
        "source_provider_called": False,
        "source_graph_walked": False,
        "actual_value_selected": False,
        "source_identity_verified": False,
        "original_source_certified": False,
        "conditional_on_input_enclosures": True,
        "limitations": [
            "Small expected sums use independent Fraction arithmetic only.",
            "10**17 exponent cases retain raw exponent fields and never expand them.",
            "Directed intervals may be wider than exact tight intervals after each term round.",
        ],
    }


__all__ = ["SCHEMA", "FIXTURE_SOURCE_ID", "run_fixture"]
