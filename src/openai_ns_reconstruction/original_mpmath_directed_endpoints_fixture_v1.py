"""Synthetic checks for :mod:`original_mpmath_directed_endpoints_v1`.

The fixture supplies a fresh ``FakeInterval`` carrying raw endpoint tuples.
It never imports mpmath, invokes a provider, reads a source graph, or calls an
older conversion helper. A caller passes the implementation module (or a
mapping exposing its public API) to :func:`run_fixture`.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_mpmath_directed_endpoints_fixture_v1"


class FakeInterval:
    """Small interval-shaped object with the private member under test."""

    __slots__ = ("_mpi_",)

    def __init__(self, lower: Any, upper: Any):
        self._mpi_ = [lower, upper]


def _member(api: Any, name: str, default: Any = None) -> Any:
    if isinstance(api, Mapping):
        return api.get(name, default)
    return getattr(api, name, default)


def _raises(call: Callable[[], Any]) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _convertible(converter: Callable[[Any], Any], value: Any) -> Any:
    """Use the injected callable only; this helper never imports a backend."""

    return converter(value)


def run_fixture(api: Any) -> dict[str, Any]:
    """Run exact representation and resource-boundary checks."""

    make = _member(api, "make_converter")
    direct = _member(api, "convert")
    checks: dict[str, bool] = {}
    if not callable(make) or not callable(direct):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"converter_api_callable": False},
            "synthetic_only": True,
            "source_provider_called": False,
        }

    converter = make(max_bits=64)
    checks["converter_callable"] = callable(converter)

    positive = FakeInterval((0, 3, -1, 2), (0, 5, -1, 3))
    negative = FakeInterval((1, 3, -1, 2), (1, 5, -2, 3))
    zero = FakeInterval((0, 0, 0, 0), (1, 0, 0, 0))

    positive_result = _convertible(converter, positive)
    negative_result = _convertible(converter, negative)
    zero_result = _convertible(converter, zero)
    checks["positive_dyadic_exact"] = positive_result == (
        Fraction(3, 2),
        Fraction(5, 2),
    )
    checks["negative_dyadic_exact"] = negative_result == (
        Fraction(-5, 4),
        Fraction(-3, 2),
    )
    checks["zero_exact_and_signed_zero_safe"] = zero_result == (
        Fraction(0),
        Fraction(0),
    )
    checks["returns_fraction_pair"] = all(
        isinstance(item, Fraction)
        for item in positive_result + negative_result + zero_result
    )

    # The lower endpoint is deliberately negative and the upper endpoint
    # positive; the converter must retain direction and never swap values.
    checks["endpoint_direction_retained"] = positive_result[0] <= positive_result[1]
    checks["reversed_direction_rejected"] = _raises(
        lambda: converter(FakeInterval((0, 1, 0, 1), (1, 1, 0, 1)))
    )

    special_codes = (
        (0, 0, -456, -2),  # mpmath 1.3.0 +inf
        (1, 0, -789, -3),  # mpmath 1.3.0 -inf
        (0, 0, -123, -1),  # mpmath 1.3.0 nan
        (1, 0, 0, -2),  # unknown signed special code
        (0, 0, 1, 0),   # unknown zero-mantissa code
    )
    checks["all_special_codes_rejected"] = all(
        _raises(lambda code=code: converter(FakeInterval(code, code)))
        for code in special_codes
    )

    illegal_fields = (
        (2, 1, 0, 1),    # invalid sign
        (-1, 1, 0, 1),   # invalid sign
        (0, -1, 0, 1),   # negative mantissa
        (0, 3, -1, 3),   # incorrect bitcount
        (0, 1, 0, -1),   # negative bitcount / special sentinel
    )
    checks["illegal_sign_mantissa_bitcount_rejected"] = all(
        _raises(lambda code=code: converter(FakeInterval(code, code)))
        for code in illegal_fields
    )

    checks["malformed_endpoint_rejected"] = (
        _raises(lambda: converter(FakeInterval((0, 1, 0), (0, 1, 0, 1))))
        and _raises(lambda: converter(FakeInterval((0, 1, 0, 1, 0), (0, 1, 0, 1))))
    )
    checks["malformed_outer_member_rejected"] = (
        _raises(lambda: converter(type("NoPair", (), {"_mpi_": [(0, 1, 0, 1)]})()))
        and _raises(
            lambda: converter(
                type(
                    "ThreePair",
                    (),
                    {"_mpi_": [(0, 1, 0, 1), (0, 1, 0, 1), (0, 1, 0, 1)]},
                )()
            )
        )
    )

    # The shifts below must be rejected before an exact integer is built.
    wide_mantissa = (0, 1 << 70, 0, 71)
    wide_exponent = (0, 1, 1 << 70, 1)
    checks["oversized_mantissa_budget_rejected"] = _raises(
        lambda: converter(FakeInterval(wide_mantissa, wide_mantissa))
    )
    checks["oversized_exponent_budget_rejected"] = _raises(
        lambda: converter(FakeInterval(wide_exponent, wide_exponent))
    )
    checks["invalid_max_bits_rejected"] = (
        _raises(lambda: make(max_bits=0))
        and _raises(lambda: make(max_bits=65.0))
    )

    checks["unknown_version_rejected"] = (
        _raises(lambda: make(version="1.2.0"))
        and _raises(lambda: direct(zero, version="2.0.0"))
    )
    checks["direct_scalars_rejected"] = (
        _raises(lambda: converter(1))
        and _raises(lambda: converter(1.0))
        and _raises(lambda: converter("1"))
    )
    float_field = (0.0, 1, 0, 1)
    checks["float_fields_rejected"] = _raises(
        lambda: converter(FakeInterval(float_field, float_field))
    )

    receipt = getattr(converter, "receipt", {})
    checks["receipt_scope_is_unverified"] = (
        isinstance(receipt, Mapping)
        and receipt.get("version_family") == "1.3.x"
        and receipt.get("conversion_contract_recognized") is True
        and receipt.get("version_supplied_by_caller") is True
        and receipt.get("runtime_version_verified") is False
        and receipt.get("source_provenance_unknown") is True
        and receipt.get("provider_verified") is False
        and receipt.get("source_verified") is False
        and receipt.get("rounding_verified") is False
        and receipt.get("actual_value_selected") is False
    )

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(checks.values()),
        "checks": checks,
        "test_count": len(checks),
        "synthetic_only": True,
        "source_provider_called": False,
        "source_graph_walked": False,
        "actual_value_selected": False,
        "limitations": [
            "Fresh FakeInterval tuples only; no mpmath object was imported.",
            "The receipt covers representation semantics, not provider/source/rounding certification.",
        ],
    }


__all__ = ["SCHEMA", "FakeInterval", "run_fixture"]
