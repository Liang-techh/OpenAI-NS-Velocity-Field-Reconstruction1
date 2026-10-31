"""Fresh synthetic checks for :mod:`original_q4_fraction_interval_v1`.

The fixture receives an implementation module (or an object exposing the
same public names) instead of importing a provider or an archived graph.  It
is intentionally inert at import time; a caller may invoke ``run_fixture``
after freezing the implementation.  Its polynomial and all expected values
are exact rationals, and its receipt checks retain the conditional nature of
the interval result.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any


SCHEMA = "original_q4_fraction_interval_fixture_v1"


def _member(api: Any, name: str, default: Any = None) -> Any:
    if isinstance(api, Mapping):
        return api.get(name, default)
    return getattr(api, name, default)


def _raises(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _interval_bounds(value: Any) -> tuple[Any, Any] | None:
    interval = getattr(value, "interval", None)
    if interval is None:
        return None
    lower = getattr(interval, "lo", getattr(interval, "lower", None))
    upper = getattr(interval, "hi", getattr(interval, "upper", None))
    if lower is None or upper is None:
        return None
    return lower, upper


def _new_backend(factory: Any) -> Any:
    if not callable(factory):
        return None
    try:
        return factory()
    except TypeError:
        # A directly supplied backend instance is also accepted by the
        # caller-facing fixture API.
        return factory


def run_fixture(api: Any) -> dict[str, Any]:
    """Run exact arithmetic and conditional quadrature checks.

    ``api`` is caller supplied.  No provider, source graph, binder, or old
    helper is imported or called, and passing this fixture does not certify
    any same-named original source.
    """

    interval_type = _member(api, "Interval")
    integrate = _member(api, "range_integrate")
    integrate_interval = _member(api, "range_integrate_interval")
    backend_factory = _member(api, "Backend")
    extract = _member(api, "extract_component_interval")
    checks: dict[str, bool] = {}

    if not callable(interval_type) or not callable(integrate):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"interval_and_integrator_callable": False},
            "limitations": ["No injected exact interval implementation was supplied."],
        }

    def interval(lo: Any, hi: Any) -> Any:
        return interval_type(lo, hi)

    # Synthetic source function: x^3 + 2 x + 1 on [-1, 2].  Its exact
    # integral is F(2)-F(-1) = 39/4.
    expected = Fraction(39, 4)

    def polynomial(cell: Any) -> Any:
        return cell * cell * cell + 2 * cell + 1

    coarse = integrate(polynomial, -1, 2, partitions=7)
    fine = integrate(polynomial, interval(-1, -1), interval(2, 2), partitions=28)
    coarse_bounds = _interval_bounds(coarse)
    fine_bounds = _interval_bounds(fine)
    checks["synthetic_true_value_is_enclosed"] = bool(
        coarse_bounds is not None
        and coarse_bounds[0] <= expected <= coarse_bounds[1]
        and fine_bounds is not None
        and fine_bounds[0] <= expected <= fine_bounds[1]
    )
    checks["grid_refinement_narrows_enclosure"] = bool(
        coarse_bounds is not None
        and fine_bounds is not None
        and (fine_bounds[1] - fine_bounds[0]) < (coarse_bounds[1] - coarse_bounds[0])
    )
    checks["coverage_cells_are_recorded"] = (
        getattr(coarse, "coverage_count", None) == 7
        and getattr(fine, "coverage_count", None) == 28
        and getattr(coarse, "kind", None) == "conditional_enclosure"
    )
    checks["conditional_flags_do_not_claim_source_certification"] = (
        getattr(coarse, "provider_certified", None) is False
        and getattr(coarse, "actual_value_selected", None) is False
    )

    positive = interval(1, 3)
    negative = interval(-3, -1)
    checks["positive_and_negative_signs"] = (
        positive.is_positive() is True and negative.is_positive() is False
    )
    quotient = positive / interval(2, 4)
    checks["strict_positive_interval_division"] = (
        getattr(quotient, "lo", None) == Fraction(1, 4)
        and getattr(quotient, "hi", None) == Fraction(3, 2)
    )
    negative_quotient = negative / interval(2, 4)
    checks["negative_numerator_interval_division"] = (
        negative_quotient.lo == Fraction(-3, 2)
        and negative_quotient.hi == Fraction(-1, 4)
    )
    checks["division_by_zero_is_rejected"] = _raises(lambda: positive / 0)
    checks["nonpositive_interval_denominator_is_rejected"] = _raises(
        lambda: positive / interval(-1, 1)
    )

    checks["float_endpoint_is_rejected"] = _raises(lambda: interval(0.0, 1))
    checks["float_scalar_is_rejected"] = _raises(lambda: positive + 0.5)
    checks["point_only_scalar_integrand_is_rejected"] = _raises(
        lambda: integrate(
            lambda cell: (cell.lo + cell.hi) / 2,
            0,
            1,
            partitions=2,
        )
    )
    checks["nonpositive_partition_count_is_rejected"] = _raises(
        lambda: integrate(polynomial, 0, 1, partitions=0)
    )

    if callable(integrate_interval):
        receipts = []
        compatibility = integrate_interval(polynomial, 0, 1, partitions=2, receipt_sink=receipts)
        checks["compatibility_integrator_keeps_interval_tag"] = (
            getattr(compatibility, "kind", None) == "enclosure"
            and len(receipts) == 1 and receipts[0].interval == compatibility
        )
    else:
        checks["compatibility_integrator_keeps_interval_tag"] = False

    backend = _new_backend(backend_factory)
    if backend is None:
        checks["backend_exact_operations"] = False
        checks["backend_rejects_transcendental"] = False
    else:
        rational = backend.rational(3, 2)
        summed = backend.add(rational, backend.rational(1, 2))
        multiplied = backend.multiply(summed, backend.rational(2))
        checks["backend_exact_operations"] = (
            getattr(rational, "kind", None) == "enclosure"
            and getattr(summed, "kind", None) == "enclosure"
            and getattr(multiplied, "kind", None) == "enclosure"
            and getattr(multiplied, "lo", None) == Fraction(4)
            and getattr(multiplied, "hi", None) == Fraction(4)
            and backend.is_positive(rational) is True
        )
        checks["backend_rejects_transcendental"] = _raises(
            lambda: backend.unary("exp", rational)
        )
        checks["backend_rejects_float"] = _raises(lambda: backend.rational(1.0, 2))

    class _SyntheticComponent:
        kind = "enclosure"
        actual_value_selected = False
        interval = object()

    if callable(extract):
        extracted = extract(
            _SyntheticComponent(),
            lambda _raw: (Fraction(-5, 2), Fraction(7, 3)),
        )
        checks["explicit_component_endpoint_extraction"] = (
            getattr(extracted, "kind", None) == "conditional_enclosure"
            and extracted.interval.lo == Fraction(-5, 2)
            and extracted.interval.hi == Fraction(7, 3)
            and extracted.endpoint_conversion_verified is False
        )
        checks["component_float_conversion_is_rejected"] = _raises(
            lambda: extract(_SyntheticComponent(), lambda _raw: (0.0, 1.0))
        )
    else:
        checks["explicit_component_endpoint_extraction"] = False
        checks["component_float_conversion_is_rejected"] = False

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "limitations": [
            "Synthetic rational callbacks only; no original graph or provider was inspected.",
            "Conditional enclosures retain caller-supplied validity and integrability assumptions.",
            "A passing fixture is not a source, rounding, quadrature, or physical certification.",
        ],
    }


__all__ = ["SCHEMA", "run_fixture"]
