"""Independent exact coverage and accumulation guard for Q4 integrals.

The older Q4 interval runtime keeps the cells *claimed* by an injected
integrator.  This module is a separate boundary that checks the cells which
were actually passed to the callback and independently recomputes the
interval accumulation.  It deliberately does not import an interval class,
an integrator, an evaluator, or any source/graph helper.  Both are injected
by the caller.

Only exact ``int`` and :class:`fractions.Fraction` endpoints are accepted.
No endpoint is converted to a midpoint or another point sample.  A successful
result consequently certifies structural coverage and exact endpoint
accumulation only; callback source identity, whole-cell mathematical range
semantics, and certified rounding remain outside this guard's claim.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_q4_integral_coverage_guard_v1"

# This is a resource bound, not a numerical accuracy bound.  It also bounds
# the number of callback records that the guard will retain.
MAX_PARTITIONS = 100_000
_MAX_PARTITIONS = MAX_PARTITIONS
_MISSING = object()


class CoverageGuardError(RuntimeError):
    """Fail-closed error raised for an invalid or unverifiable integral."""

    def __init__(self, code: str, detail: Any = None):
        self.code = code
        self.detail = detail
        message = code if detail is None else f"{code}: {detail}"
        super().__init__(message)


# Descriptive aliases make the failure type easy to discover without making
# callers depend on a private implementation exception.
IntegralCoverageGuardError = CoverageGuardError
GuardRejected = CoverageGuardError


def _reject(code: str, detail: Any = None) -> None:
    raise CoverageGuardError(code, detail)


def _member(value: Any, name: str) -> Any:
    """Read one explicitly named member without invoking arbitrary helpers."""

    try:
        if isinstance(value, Mapping):
            if name in value:
                return value[name]
            return _MISSING
        return getattr(value, name)
    except Exception:
        return _MISSING


def _exact_fraction(value: Any, label: str) -> Fraction:
    """Normalize only exact integers and Fractions.

    ``bool`` is intentionally excluded even though it subclasses ``int``.
    No float, decimal, string, or generic numeric conversion is attempted.
    """

    if type(value) is int:
        return Fraction(value, 1)
    if isinstance(value, Fraction):
        try:
            numerator = value.numerator
            denominator = value.denominator
            if type(numerator) is not int or type(denominator) is not int:
                _reject("non_exact_interval_endpoint", label)
            if denominator <= 0:
                _reject("invalid_fraction_endpoint", label)
            return Fraction(numerator, denominator)
        except CoverageGuardError:
            raise
        except Exception as exc:
            raise CoverageGuardError(
                "invalid_fraction_endpoint", label
            ) from exc
    _reject("non_exact_interval_endpoint", label)
    raise AssertionError("unreachable")


def _checked_interval(
    value: Any,
    interval_type: type,
    code: str,
    *,
    strict_positive: bool = False,
) -> tuple[Fraction, Fraction]:
    """Validate an exact interval object and return normalized endpoints."""

    if type(value) is not interval_type:
        _reject(code, "interval type")
    kind = _member(value, "kind")
    if kind != "enclosure":
        _reject(code, "interval kind")
    raw_lo = _member(value, "lo")
    raw_hi = _member(value, "hi")
    if raw_lo is _MISSING or raw_hi is _MISSING:
        _reject(code, "missing .lo/.hi")
    lower = _exact_fraction(raw_lo, f"{code}.lo")
    upper = _exact_fraction(raw_hi, f"{code}.hi")
    if lower > upper:
        _reject(code, "reversed interval")
    if strict_positive and not lower < upper:
        _reject(code, "interval must have strictly positive width")
    return lower, upper


def _checked_bound(value: Any, interval_type: type, label: str) -> Fraction:
    """Validate one degenerate exact interval endpoint."""

    lower, upper = _checked_interval(value, interval_type, label)
    if lower != upper:
        _reject("bound_must_be_degenerate", label)
    return lower


def _checked_partitions(value: Any) -> int:
    if type(value) is not int or value <= 0:
        _reject("partitions_must_be_positive_integer")
    if value > MAX_PARTITIONS:
        _reject("partitions_exceed_limit", MAX_PARTITIONS)
    return value


def _sequence(value: Any, label: str) -> tuple[Any, ...]:
    """Accept only finite, already-materialized coverage containers."""

    if not isinstance(value, (tuple, list)):
        _reject("coverage_cells_must_be_finite_sequence", label)
    try:
        size = len(value)
    except Exception as exc:
        raise CoverageGuardError("coverage_cells_unreadable", label) from exc
    if size > MAX_PARTITIONS:
        _reject("coverage_cells_exceed_limit", label)
    return tuple(value)


def _validate_ordered_cells(
    cells: tuple[Any, ...],
    *,
    interval_type: type,
    left: Fraction,
    right: Fraction,
    partitions: int,
    label: str,
) -> tuple[tuple[Any, Fraction, Fraction], ...]:
    """Check an ordered complete partition without sorting it."""

    if len(cells) != partitions:
        _reject("coverage_count_mismatch", label)
    cursor = left
    checked: list[tuple[Any, Fraction, Fraction]] = []
    for index, cell in enumerate(cells):
        cell_left, cell_right = _checked_interval(
            cell,
            interval_type,
            "cell_not_exact_enclosure",
            strict_positive=True,
        )
        if cell_left < left or cell_right > right:
            _reject("cell_outside_requested_range", f"{label}[{index}]")
        # Requiring the next left endpoint to equal the prior right endpoint
        # rejects gaps, overlaps, and reversed traversal without sorting away
        # the evidence's original order.
        if cell_left != cursor:
            _reject("coverage_not_strictly_contiguous", f"{label}[{index}]")
        checked.append((cell, cell_left, cell_right))
        cursor = cell_right
    if cursor != right:
        _reject("coverage_does_not_reach_upper_bound", label)
    return tuple(checked)


def _same_endpoint_pairs(
    actual: tuple[tuple[Any, Fraction, Fraction], ...],
    declared: tuple[tuple[Any, Fraction, Fraction], ...],
) -> None:
    if len(actual) != len(declared):
        _reject("actual_and_declared_coverage_count_mismatch")
    for index, (actual_item, declared_item) in enumerate(zip(actual, declared)):
        _, actual_left, actual_right = actual_item
        _, declared_left, declared_right = declared_item
        if actual_left != declared_left or actual_right != declared_right:
            _reject("actual_and_declared_coverage_mismatch", index)


@dataclass(frozen=True, slots=True)
class GuardedIntegralEnclosure:
    """Conditional integral result produced after guard verification.

    The interval is the underlying integrator's interval object, unchanged.
    In particular, a deliberately wider valid interval is never narrowed to
    the independently recomputed sum.
    """

    interval: Any
    coverage_cells: tuple[Any, ...] = field(default_factory=tuple)
    kind: str = field(default="conditional_enclosure", init=False)
    provider_certified: bool = field(default=False, init=False)
    actual_value_selected: bool = field(default=False, init=False)
    structural_coverage_verified: bool = field(default=True, init=False)
    exact_accumulation_verified: bool = field(default=True, init=False)

    @property
    def coverage_count(self) -> int:
        return len(self.coverage_cells)

    @property
    def count(self) -> int:
        return self.coverage_count

    @property
    def cells(self) -> tuple[Any, ...]:
        return self.coverage_cells

    @property
    def conditional(self) -> bool:
        return True

    def receipt(self) -> dict[str, Any]:
        """Return the narrow non-certifying claim made by this guard."""

        return {
            "schema": SCHEMA,
            "kind": self.kind,
            "coverage_count": self.coverage_count,
            "structural_coverage_verified": True,
            "exact_accumulation_verified": True,
            "coverage_verified": True,
            "accumulation_verified": True,
            "actual_value_selected": False,
            "provider_certified": False,
            # These deliberately remain false.  Structural validation of the
            # cell objects does not prove what the callback mathematically
            # encloses, where it came from, or how it rounded endpoints.
            "callback_cell_structure_verified": True,
            "callback_whole_cell_range_verified": False,
            "callback_source_verified": False,
            "callback_certified_rounding_verified": False,
            "whole_cell_range_verified": False,
            "source_verified": False,
            "certified_rounding_verified": False,
            "source_identity_verified": False,
            "condition": (
                "structural coverage and exact Fraction accumulation verified "
                "only; callback whole-cell range/source/certified rounding "
                "remain unverified"
            ),
        }

    def __float__(self) -> float:
        raise TypeError("an integral enclosure cannot be converted to float")

    def __bool__(self) -> bool:
        raise TypeError("an integral enclosure has no scalar truth value")


class GuardedIntegrator:
    """Wrap an injected range integrator with exact coverage checks.

    Exact API::

        GuardedIntegrator(underlying, interval_type).integrate(
            callback, lo, hi, *, partitions
        )

    ``lo`` and ``hi`` are exact, degenerate instances of ``interval_type``;
    their endpoint values define the requested range.  The underlying object
    must expose ``integrate(callback, lo, hi, partitions=...)``.  A
    ``range_integrate`` method is accepted as a compatibility spelling when
    ``integrate`` is absent.
    """

    def __init__(self, underlying: Any, interval_type: type) -> None:
        if not isinstance(interval_type, type):
            _reject("interval_type_must_be_a_type")
        if underlying is None:
            _reject("underlying_integrator_required")
        self.underlying = underlying
        self.interval_type = interval_type

    def _method(self) -> Callable[..., Any]:
        method = _member(self.underlying, "integrate")
        if not callable(method):
            method = _member(self.underlying, "range_integrate")
        if not callable(method):
            _reject("underlying_integrator_method_required")
        return method

    def integrate(
        self,
        callback: Callable[[Any], Any],
        lo: Any,
        hi: Any,
        *,
        partitions: int,
    ) -> GuardedIntegralEnclosure:
        """Run and verify one exact full-cell integral, failing closed."""

        if not callable(callback):
            _reject("callback_must_be_callable")
        count = _checked_partitions(partitions)
        left = _checked_bound(lo, self.interval_type, "lo")
        right = _checked_bound(hi, self.interval_type, "hi")
        if not left < right:
            _reject("integration_bounds_must_be_strictly_increasing")

        events: list[tuple[Any, Any, Fraction, Fraction, Fraction, Fraction]] = []
        callback_failure: CoverageGuardError | None = None

        def tracked_callback(cell: Any) -> Any:
            nonlocal callback_failure
            if callback_failure is not None:
                raise callback_failure
            try:
                if len(events) >= count:
                    _reject("callback_count_exceeds_partitions")
                cell_left, cell_right = _checked_interval(
                    cell,
                    self.interval_type,
                    "callback_cell_not_exact_enclosure",
                    strict_positive=True,
                )
                if cell_left < left or cell_right > right:
                    _reject("callback_cell_outside_requested_range")
                returned = callback(cell)
                if _checked_interval(cell, self.interval_type, "callback_cell_changed", strict_positive=True) != (cell_left, cell_right):
                    _reject("callback_cell_changed")
                returned_left, returned_right = _checked_interval(
                    returned,
                    self.interval_type,
                    "callback_result_not_exact_enclosure",
                )
            except CoverageGuardError as exc:
                callback_failure = exc
                raise
            except Exception as exc:
                callback_failure = CoverageGuardError(
                    "callback_or_cell_validation_error", type(exc).__name__
                )
                raise callback_failure from exc
            events.append(
                (
                    cell,
                    returned,
                    cell_left,
                    cell_right,
                    returned_left,
                    returned_right,
                )
            )
            return returned

        try:
            underlying_result = self._method()(
                tracked_callback,
                lo,
                hi,
                partitions=count,
            )
        except CoverageGuardError as exc:
            if callback_failure is not None:
                raise callback_failure
            raise exc
        except Exception as exc:
            if callback_failure is not None:
                raise callback_failure
            raise CoverageGuardError(
                "underlying_integrator_error", type(exc).__name__
            ) from exc

        # An underlying integrator that catches a callback exception must not
        # be able to turn that failed call into a successful receipt.
        if callback_failure is not None:
            raise callback_failure
        if len(events) != count:
            _reject("callback_count_mismatch", len(events))

        kind = _member(underlying_result, "kind")
        if kind != "conditional_enclosure":
            _reject("underlying_result_not_conditional_enclosure")
        actual_value_selected = _member(underlying_result, "actual_value_selected")
        if actual_value_selected is not False:
            _reject("underlying_result_selected_actual_value")
        raw_interval = _member(underlying_result, "interval")
        if raw_interval is _MISSING:
            _reject("underlying_result_interval_missing")
        result_left, result_right = _checked_interval(
            raw_interval,
            self.interval_type,
            "underlying_result_interval_invalid",
        )
        raw_coverage = _member(underlying_result, "coverage_cells")
        if raw_coverage is _MISSING or raw_coverage is None:
            _reject("underlying_result_coverage_cells_missing")
        declared_cells = _sequence(raw_coverage, "underlying_result.coverage_cells")

        actual_cells = tuple(item[0] for item in events)
        actual_checked = _validate_ordered_cells(
            actual_cells,
            interval_type=self.interval_type,
            left=left,
            right=right,
            partitions=count,
            label="actual callback cells",
        )
        for event, checked in zip(events, actual_checked):
            if (event[2], event[3]) != (checked[1], checked[2]):
                _reject("callback_cell_changed_after_call")
        declared_checked = _validate_ordered_cells(
            declared_cells,
            interval_type=self.interval_type,
            left=left,
            right=right,
            partitions=count,
            label="declared coverage cells",
        )
        _same_endpoint_pairs(actual_checked, declared_checked)

        # Independently recompute the enclosure using only Fraction endpoint
        # arithmetic.  No interval operation supplied by the underlying
        # implementation participates in this check.
        accumulated_left = Fraction(0, 1)
        accumulated_right = Fraction(0, 1)
        for _, _, cell_left, cell_right, value_left, value_right in events:
            width = cell_right - cell_left
            accumulated_left += width * value_left
            accumulated_right += width * value_right
        if result_left > accumulated_left or result_right < accumulated_right:
            _reject("underlying_interval_does_not_contain_exact_accumulation")

        return GuardedIntegralEnclosure(
            interval=raw_interval,
            coverage_cells=declared_cells,
        )


# A concise alias is useful to code that names all integral result classes
# ``IntegralEnclosure`` while retaining the explicit guard name above.
GuardedIntegralResult = GuardedIntegralEnclosure


__all__ = [
    "SCHEMA",
    "MAX_PARTITIONS",
    "CoverageGuardError",
    "IntegralCoverageGuardError",
    "GuardRejected",
    "GuardedIntegralEnclosure",
    "GuardedIntegralResult",
    "GuardedIntegrator",
]
