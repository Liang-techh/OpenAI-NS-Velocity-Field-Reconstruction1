"""Independent synthetic fixture for the unbounded dyadic product leaf.

``run_fixture(leaf)`` accepts the implementation module or a mapping exposing
its public API.  It uses raw finite endpoint tuples and compares against
hand-written canonical tuples.  The small rational cases are checked by
explicit integer encodings; the astronomical-exponent case is never converted
to a rational integer and is tested only through exponent cancellation and
preservation.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, Callable


SCHEMA = "original_unbounded_dyadic_product_fixture_v1"
FIXTURE_SOURCE_ID = "synthetic:unbounded-dyadic-product-fixture-v1"


def _member(leaf: Any, name: str, default: Any = None) -> Any:
    if isinstance(leaf, Mapping):
        return leaf.get(name, default)
    return getattr(leaf, name, default)


def _operation(leaf: Any) -> Callable[..., Any] | None:
    for name in (
        "interval_product",
        "multiply_interval",
        "interval_mul",
        "multiply_intervals",
        "conditional_interval_product",
    ):
        candidate = _member(leaf, name)
        if callable(candidate):
            return candidate
    return None


def _raises(call: Callable[[], Any]) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _same_interval(actual: Any, expected: list[list[int]]) -> bool:
    """Compare endpoint fields while tolerating tuple/list container spelling."""

    if not isinstance(actual, (tuple, list)) or len(actual) != 2:
        return False
    converted: list[list[int]] = []
    for endpoint in actual:
        if not isinstance(endpoint, (tuple, list)) or len(endpoint) != 4:
            return False
        converted.append(list(endpoint))
    return converted == expected


def _matches(
    operation: Callable[..., Any],
    left: Any,
    right: Any,
    expected: list[list[int]],
) -> bool:
    try:
        return _same_interval(operation(left, right), expected)
    except Exception:
        return False


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run fresh finite-tuple product checks without touching any source data."""

    operation = _operation(leaf)
    checks: dict[str, bool] = {}
    if operation is None:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"interval_product_api_callable": False},
            "test_count": 1,
            "synthetic_only": True,
            "source_provider_called": False,
            "source_graph_walked": False,
            "original_source_certified": False,
        }

    # 2..3 times 4..5 gives the independently written exact interval 8..15.
    positive_left = ((0, 1, 1, 1), (0, 3, 0, 2))
    positive_right = ((0, 1, 2, 1), (0, 5, 0, 3))
    checks["positive_exact_integer_interval"] = _matches(
        operation,
        positive_left,
        positive_right,
        [[0, 1, 3, 1], [0, 15, 0, 4]],
    )

    # [-5,-2] times [3,4] gives [-20,-6], including signed ordering.
    negative_left = ((1, 5, 0, 3), (1, 1, 1, 1))
    negative_right = ((0, 3, 0, 2), (0, 1, 2, 1))
    checks["negative_exact_integer_interval"] = _matches(
        operation,
        negative_left,
        negative_right,
        [[1, 5, 2, 3], [1, 3, 1, 2]],
    )

    # [-2,3] times [4,5] crosses zero and gives [-10,15].
    cross_left = ((1, 1, 1, 1), (0, 3, 0, 2))
    cross_right = ((0, 1, 2, 1), (0, 5, 0, 3))
    checks["cross_zero_exact_integer_interval"] = _matches(
        operation,
        cross_left,
        cross_right,
        [[1, 5, 1, 3], [0, 15, 0, 4]],
    )

    zero_left = ((0, 0, 0, 0), (1, 0, 0, 0))
    zero_right = ((1, 7, 0, 3), (0, 9, 0, 4))
    checks["zero_interval_exact_and_signed_zero_safe"] = _matches(
        operation,
        zero_left,
        zero_right,
        [[0, 0, 0, 0], [0, 0, 0, 0]],
    )

    # These exponents are deliberately astronomical.  The expected value is
    # written directly: (3*2**E)*(5*2**(-E)) = 15, without a conversion.
    extreme = 10**17
    extreme_left = ((0, 3, extreme, 2), (0, 3, extreme, 2))
    extreme_right = ((0, 5, -extreme, 3), (0, 5, -extreme, 3))
    checks["extreme_exponent_cancellation_without_conversion"] = _matches(
        operation,
        extreme_left,
        extreme_right,
        [[0, 15, 0, 4], [0, 15, 0, 4]],
    )

    preserve_left = ((0, 1, extreme, 1), (0, 1, extreme, 1))
    preserve_right = ((0, 1, 0, 1), (0, 1, 0, 1))
    checks["extreme_exponent_preserved_raw"] = _matches(
        operation,
        preserve_left,
        preserve_right,
        [[0, 1, extreme, 1], [0, 1, extreme, 1]],
    )

    sign_left = ((1, 3, 0, 2), (1, 3, 0, 2))
    sign_right = ((0, 5, 0, 3), (0, 5, 0, 3))
    checks["integer_multiplication_sign_xor"] = _matches(
        operation,
        sign_left,
        sign_right,
        [[1, 15, 0, 4], [1, 15, 0, 4]],
    )

    # Noncanonical finite input is normalized by stripping the two trailing
    # zero bits: 6*1 is represented as 3*2**1.
    checks["trailing_zero_normalization"] = _matches(
        operation,
        ((0, 6, 0, 3), (0, 6, 0, 3)),
        ((0, 1, 0, 1), (0, 1, 0, 1)),
        [[0, 3, 1, 2], [0, 3, 1, 2]],
    )

    special_codes = (
        (0, 0, -456, -2),  # +Inf
        (1, 0, -789, -3),  # -Inf
        (0, 0, -123, -1),  # NaN
        (1, 0, 0, -2),     # unknown special sentinel
        (0, 0, 1, 0),      # nonstandard zero-mantissa encoding
    )
    checks["nan_inf_and_special_zero_rejected"] = all(
        _raises(lambda code=code: operation((code, code), ((0, 1, 0, 1), (0, 1, 0, 1))))
        for code in special_codes
    )

    malformed_bitcounts = (
        (0, 3, 0, 3),  # bitcount is not 2
        (0, 1, 0, 0),  # nonzero mantissa cannot have zero bitcount
        (0, 1, 0, -1),
    )
    checks["malformed_bitcounts_rejected"] = all(
        _raises(lambda code=code: operation((code, code), ((0, 1, 0, 1), (0, 1, 0, 1))))
        for code in malformed_bitcounts
    )

    reversed_interval = ((0, 3, 0, 2), (0, 1, 1, 1))
    singleton_one = ((0, 1, 0, 1), (0, 1, 0, 1))
    checks["reversed_input_interval_rejected"] = _raises(
        lambda: operation(reversed_interval, singleton_one)
    )
    checks["scalar_and_short_endpoint_rejected"] = (
        _raises(lambda: operation(1, singleton_one))
        and _raises(lambda: operation(((0, 1, 0), (0, 1, 0)), singleton_one))
    )

    oversized = (0, 1 << 4096, 0, 4097)
    checks["unsupported_large_mantissa_rejected"] = _raises(
        lambda: operation((oversized, oversized), singleton_one)
    )
    excessive_exponent = (0, 1, 1 << 256, 1)
    checks["unsupported_exponent_size_rejected"] = _raises(
        lambda: operation((excessive_exponent, excessive_exponent), singleton_one)
    )
    max_allowed_exponent = (0, 1, (1 << 256) - 1, 1)
    checks["result_exponent_overflow_rejected"] = _raises(
        lambda: operation((max_allowed_exponent, max_allowed_exponent), ((0, 1, 1, 1), (0, 1, 1, 1)))
    )

    factory = _member(leaf, "make_kernel")
    if callable(factory):
        try:
            kernel = factory(
                left_source_id="fixture:left-source",
                right_source_id="fixture:right-source",
            )
            receipt = getattr(kernel, "receipt", {})
            source_ids = getattr(kernel, "source_ids", {})
            checks["source_ids_retained_separately"] = (
                source_ids == {
                    "left": "fixture:left-source",
                    "right": "fixture:right-source",
                }
                and receipt.get("source_ids_separate") is True
                and receipt.get("conditional_on_input_enclosures") is True
                and receipt.get("point_value_selected") is False
                and receipt.get("source_identity_verified") is False
                and receipt.get("original_source_certified") is False
            )
            checks["kernel_returns_raw_interval"] = _matches(
                kernel,
                positive_left,
                positive_right,
                [[0, 1, 3, 1], [0, 15, 0, 4]],
            )
        except Exception:
            checks["source_ids_retained_separately"] = False
            checks["kernel_returns_raw_interval"] = False
    else:
        checks["source_ids_retained_separately"] = False
        checks["kernel_returns_raw_interval"] = False

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(checks.values()),
        "checks": checks,
        "test_count": len(checks),
        "synthetic_only": True,
        "source_provider_called": False,
        "source_graph_walked": False,
        "actual_value_selected": False,
        "synthetic_singleton_inputs_used": True,
        "actual_source_point_selected": False,
        "source_identity_verified": False,
        "original_source_certified": False,
        "conditional_on_input_enclosures": True,
        "limitations": [
            "Fresh raw finite tuples only; no mpmath/provider/source graph was accessed.",
            "The extreme exponent checks retain exponent fields and never build a huge integer.",
            "Expected endpoint tuples are hand-written and independent of the implementation helpers.",
        ],
    }


__all__ = ["SCHEMA", "FIXTURE_SOURCE_ID", "run_fixture"]
