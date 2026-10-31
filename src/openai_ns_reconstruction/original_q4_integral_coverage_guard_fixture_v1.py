"""Fresh synthetic fixture for ``original_q4_integral_coverage_guard_v1``.

This fixture owns every interval and integrator used below.  It does not
import an earlier Q4 module, graph, provider, or arithmetic helper.  The
entry point injects the new guard module (or a mapping exposing its public
names) through ``run_fixture(leaf)``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any


SCHEMA = "original_q4_integral_coverage_guard_fixture_v1"
MAX_PARTITIONS = 100_000
_MISSING = object()


def _fraction(value: Any, label: str = "endpoint") -> Fraction:
    if type(value) is int:
        return Fraction(value, 1)
    if isinstance(value, Fraction):
        return Fraction(value.numerator, value.denominator)
    raise TypeError(f"{label} must be an exact int or Fraction")


@dataclass(frozen=True, slots=True)
class FakeInterval:
    """Exact enclosure type used only by this fixture."""

    lo: Fraction | int
    hi: Fraction | int
    kind: str = field(default="enclosure", init=False)

    def __post_init__(self) -> None:
        lower = _fraction(self.lo, "FakeInterval.lo")
        upper = _fraction(self.hi, "FakeInterval.hi")
        if lower > upper:
            raise ValueError("FakeInterval lower endpoint exceeds upper endpoint")
        object.__setattr__(self, "lo", lower)
        object.__setattr__(self, "hi", upper)

    @property
    def lower(self) -> Fraction:
        return self.lo

    @property
    def upper(self) -> Fraction:
        return self.hi

    @property
    def width(self) -> Fraction:
        return self.hi - self.lo

    @staticmethod
    def _operand(value: Any) -> "FakeInterval":
        if type(value) is FakeInterval:
            return value
        scalar = _fraction(value, "FakeInterval scalar")
        return FakeInterval(scalar, scalar)

    def __add__(self, other: Any) -> "FakeInterval":
        rhs = self._operand(other)
        return FakeInterval(self.lo + rhs.lo, self.hi + rhs.hi)

    def __radd__(self, other: Any) -> "FakeInterval":
        return self + other

    def __sub__(self, other: Any) -> "FakeInterval":
        rhs = self._operand(other)
        return FakeInterval(self.lo - rhs.hi, self.hi - rhs.lo)

    def __rsub__(self, other: Any) -> "FakeInterval":
        return self._operand(other) - self

    def __mul__(self, other: Any) -> "FakeInterval":
        rhs = self._operand(other)
        products = (
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        )
        return FakeInterval(min(products), max(products))

    def __rmul__(self, other: Any) -> "FakeInterval":
        return self * other

    def __truediv__(self, other: Any) -> "FakeInterval":
        rhs = self._operand(other)
        if rhs.lo <= 0 <= rhs.hi:
            raise ZeroDivisionError("FakeInterval denominator contains zero")
        reciprocal = FakeInterval(
            min(Fraction(1, rhs.lo), Fraction(1, rhs.hi)),
            max(Fraction(1, rhs.lo), Fraction(1, rhs.hi)),
        )
        return self * reciprocal

    def __neg__(self) -> "FakeInterval":
        return FakeInterval(-self.hi, -self.lo)

    def __float__(self) -> float:
        raise TypeError("FakeInterval cannot become a point sample")

    def __bool__(self) -> bool:
        raise TypeError("FakeInterval has no scalar truth value")


@dataclass(frozen=True, slots=True)
class FakeIntegralResult:
    """Underlying result shape consumed by the independent guard."""

    interval: FakeInterval
    coverage_cells: tuple[FakeInterval, ...]
    kind: str = field(default="conditional_enclosure", init=False)
    actual_value_selected: bool = field(default=False, init=False)
    provider_certified: bool = field(default=False, init=False)

    def receipt(self) -> dict[str, Any]:
        return {
            "schema": SCHEMA,
            "kind": self.kind,
            "coverage_count": len(self.coverage_cells),
            "actual_value_selected": False,
            "provider_certified": False,
        }


def _bound(value: Any, label: str) -> Fraction:
    if type(value) is not FakeInterval:
        raise TypeError(f"{label} must be a FakeInterval")
    if value.lo != value.hi:
        raise ValueError(f"{label} must be degenerate")
    return value.lo


def _cells(lower: FakeInterval, upper: FakeInterval, partitions: int) -> tuple[FakeInterval, ...]:
    left = _bound(lower, "lower")
    right = _bound(upper, "upper")
    if not type(partitions) is int or partitions <= 0 or partitions > MAX_PARTITIONS:
        raise ValueError("invalid fixture partition count")
    if not left < right:
        raise ValueError("fixture bounds must be increasing")
    width = (right - left) / partitions
    return tuple(
        FakeInterval(left + width * index, left + width * (index + 1))
        for index in range(partitions)
    )


def _accumulate(cells: tuple[FakeInterval, ...], values: list[FakeInterval]) -> FakeInterval:
    if len(cells) != len(values):
        raise ValueError("fixture accumulation needs one value per cell")
    lower = Fraction(0, 1)
    upper = Fraction(0, 1)
    for cell, value in zip(cells, values):
        if type(value) is not FakeInterval:
            raise TypeError("fixture value must be a FakeInterval")
        width = cell.hi - cell.lo
        lower += width * value.lo
        upper += width * value.hi
    return FakeInterval(lower, upper)


class FakeRangeIntegrator:
    """Small scriptable integrator used to exercise each guard rejection."""

    def __init__(self, mode: str = "honest") -> None:
        self.mode = mode
        self.calls: list[Any] = []

    def integrate(
        self,
        callback: Any,
        lo: FakeInterval,
        hi: FakeInterval,
        partitions: int,
    ) -> FakeIntegralResult:
        base = _cells(lo, hi, partitions)
        if self.mode == "honest" or self.mode == "wide":
            call_cells = base
            declared = base
        elif self.mode == "lying_cells":
            call_cells = base
            declared = tuple(reversed(base))
        elif self.mode == "lying_partition":
            call_cells = base
            middle = (base[0].lo + base[1].hi) / 2 + base[0].width / 4
            declared = (FakeInterval(base[0].lo, middle), FakeInterval(middle, base[1].hi)) + base[2:]
        elif self.mode == "missing_callback":
            call_cells = base[:-1]
            declared = base
        elif self.mode == "overlap":
            width = base[0].width
            left = base[0].lo
            call_cells = (FakeInterval(left, left + 2 * width),) + tuple(
                FakeInterval(left + width * index, left + width * (index + 1))
                for index in range(1, partitions)
            )
            declared = call_cells
        elif self.mode == "gap":
            # Six positive cells with a deliberate gap after the first one.
            # The fixture always calls this mode with six partitions.
            half = base[0].width / 2
            left = base[0].lo
            call_cells = tuple(
                [FakeInterval(left, left + half)]
                + [
                    FakeInterval(left + 2 * half * index, left + 2 * half * (index + 1))
                    for index in range(1, partitions)
                ]
            )
            # The final item is widened to end at the requested upper bound;
            # all cells remain positive and in range, but the first gap stays.
            final = FakeInterval(call_cells[-1].lo, base[-1].hi)
            call_cells = call_cells[:-1] + (final,)
            declared = call_cells
        elif self.mode == "reverse":
            call_cells = tuple(reversed(base))
            declared = call_cells
        elif self.mode == "narrow":
            call_cells = base
            declared = base
        elif self.mode == "scalar_point":
            call_cells = base
            declared = base
        elif self.mode == "extra_callback":
            call_cells = base + (base[0],)
            declared = base
        else:
            raise ValueError(f"unknown fixture mode: {self.mode}")

        values: list[FakeInterval] = []
        for cell in call_cells:
            self.calls.append(cell)
            if self.mode == "scalar_point":
                argument: Any = (cell.lo + cell.hi) / 2
            else:
                argument = cell
            values.append(callback(argument))

        if self.mode == "narrow":
            result_interval = FakeInterval(6, 6)
        else:
            result_interval = _accumulate(call_cells, values)
            if self.mode == "wide":
                result_interval = FakeInterval(
                    result_interval.lo - 3, result_interval.hi + 3
                )
        return FakeIntegralResult(result_interval, tuple(declared))


class HonestIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("honest")


class WiderIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("wide")


class LyingCellsIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("lying_cells")


class MissingCallbackIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("missing_callback")


class OverlapIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("overlap")


class GapIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("gap")


class ReverseIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("reverse")


class NarrowResultIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("narrow")


class ScalarPointIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("scalar_point")


class ExtraCallbackIntegrator(FakeRangeIntegrator):
    def __init__(self) -> None:
        super().__init__("extra_callback")


# Alternate explicit names make the individual adversarial fixtures useful to
# callers without requiring knowledge of the mode-string implementation.
FakeIntegrator = FakeRangeIntegrator
MisreportedCoverageIntegrator = LyingCellsIntegrator
OmittedCallbackIntegrator = MissingCallbackIntegrator
OverlappingIntegrator = OverlapIntegrator
GappedIntegrator = GapIntegrator
ReversedOrderIntegrator = ReverseIntegrator
TooNarrowIntegrator = NarrowResultIntegrator
PointCallbackIntegrator = ScalarPointIntegrator


def g(cell: FakeInterval) -> FakeInterval:
    """Fresh input function ``g(x) = 2 x + 3`` on a complete cell."""

    if type(cell) is not FakeInterval:
        raise TypeError("g requires a complete FakeInterval cell")
    return 2 * cell + 3


def _member(api: Any, name: str) -> Any:
    if isinstance(api, dict):
        return api.get(name, _MISSING)
    return getattr(api, name, _MISSING)


def _raises(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _run_guard(guard_type: Any, integrator: Any, callback: Any = g) -> Any:
    guard = guard_type(integrator, FakeInterval)
    return guard.integrate(
        callback,
        FakeInterval(-2, -2),
        FakeInterval(1, 1),
        partitions=6,
    )


def run_fixture(leaf: Any) -> dict[str, Any]:
    """Run the exact synthetic guard checks through an injected ``leaf``."""

    guard_type = _member(leaf, "GuardedIntegrator")
    checks: dict[str, bool] = {}
    if not callable(guard_type):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"guard_callable": False},
            "limitations": ["The injected GuardedIntegrator API is missing."],
        }

    exact_value = Fraction(6, 1)
    try:
        honest = _run_guard(guard_type, HonestIntegrator())
    except Exception as exc:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"honest_integral_runs": False},
            "error": type(exc).__name__,
            "limitations": ["The synthetic guard could not produce its baseline result."],
        }

    checks["honest_integral_runs"] = True
    checks["exact_integral_is_enclosed"] = (
        type(honest.interval) is FakeInterval
        and honest.interval.lo <= exact_value <= honest.interval.hi
    )
    checks["coverage_is_retained"] = (
        honest.kind == "conditional_enclosure"
        and honest.actual_value_selected is False
        and len(honest.coverage_cells) == 6
        and all(type(cell) is FakeInterval for cell in honest.coverage_cells)
    )
    receipt = honest.receipt()
    checks["receipt_is_structural_and_exact_only"] = (
        receipt.get("structural_coverage_verified") is True
        and receipt.get("exact_accumulation_verified") is True
        and receipt.get("actual_value_selected") is False
        and receipt.get("callback_whole_cell_range_verified") is False
        and receipt.get("callback_source_verified") is False
        and receipt.get("callback_certified_rounding_verified") is False
        and "structural coverage" in receipt.get("condition", "")
        and "exact Fraction accumulation" in receipt.get("condition", "")
    )

    wider = _run_guard(guard_type, WiderIntegrator())
    checks["wider_underlying_interval_is_preserved"] = (
        wider.interval.lo < honest.interval.lo
        and wider.interval.hi > honest.interval.hi
        and wider.interval.lo <= exact_value <= wider.interval.hi
    )

    checks["underlying_lie_about_cells_rejected"] = _raises(
        lambda: _run_guard(guard_type, LyingCellsIntegrator())
    )
    checks["different_valid_declared_partition_rejected"] = _raises(
        lambda: _run_guard(guard_type, FakeRangeIntegrator("lying_partition"))
    )
    checks["missing_callback_rejected"] = _raises(
        lambda: _run_guard(guard_type, MissingCallbackIntegrator())
    )
    checks["overlapping_cells_rejected"] = _raises(
        lambda: _run_guard(guard_type, OverlapIntegrator())
    )
    checks["gapped_cells_rejected"] = _raises(
        lambda: _run_guard(guard_type, GapIntegrator())
    )
    checks["reverse_order_rejected"] = _raises(
        lambda: _run_guard(guard_type, ReverseIntegrator())
    )
    checks["too_narrow_integral_rejected"] = _raises(
        lambda: _run_guard(guard_type, NarrowResultIntegrator())
    )
    checks["scalar_point_callback_rejected"] = _raises(
        lambda: _run_guard(guard_type, ScalarPointIntegrator())
    )
    checks["extra_callback_rejected"] = _raises(
        lambda: _run_guard(guard_type, ExtraCallbackIntegrator())
    )
    checks["point_return_value_rejected"] = _raises(
        lambda: _run_guard(
            guard_type,
            HonestIntegrator(),
            lambda cell: (cell.lo + cell.hi) / 2,
        )
    )

    # Bounds and partition checks are intentionally explicit: the source
    # range endpoints are degenerate interval objects, never point samples.
    checks["nondegenerate_bound_rejected"] = _raises(
        lambda: guard_type(HonestIntegrator(), FakeInterval).integrate(
            g,
            FakeInterval(-2, -1),
            FakeInterval(1, 1),
            partitions=6,
        )
    )
    checks["raw_scalar_bound_rejected"] = _raises(
        lambda: guard_type(HonestIntegrator(), FakeInterval).integrate(
            g, -2, FakeInterval(1, 1), partitions=6
        )
    )
    checks["nonpositive_partitions_rejected"] = _raises(
        lambda: guard_type(HonestIntegrator(), FakeInterval).integrate(
            g,
            FakeInterval(-2, -2),
            FakeInterval(1, 1),
            partitions=0,
        )
    )

    checks["all"] = all(checks.values())
    return {
        "schema": SCHEMA,
        "passed": checks["all"],
        "checks": checks,
        "successful_result": honest,
        "wider_result": wider,
        "independent_exact_value": exact_value,
        "limitations": [
            "Only fresh rational cells and the synthetic g(x)=2x+3 callback were used.",
            "The guard does not certify callback source identity or directed rounding.",
            "A passing fixture is not an original graph or physical-source certification.",
        ],
    }


run = run_fixture


__all__ = [
    "SCHEMA",
    "FakeInterval",
    "FakeIntegralResult",
    "FakeRangeIntegrator",
    "FakeIntegrator",
    "HonestIntegrator",
    "WiderIntegrator",
    "LyingCellsIntegrator",
    "MisreportedCoverageIntegrator",
    "MissingCallbackIntegrator",
    "OmittedCallbackIntegrator",
    "OverlapIntegrator",
    "OverlappingIntegrator",
    "GapIntegrator",
    "GappedIntegrator",
    "ReverseIntegrator",
    "ReversedOrderIntegrator",
    "NarrowResultIntegrator",
    "TooNarrowIntegrator",
    "ScalarPointIntegrator",
    "PointCallbackIntegrator",
    "ExtraCallbackIntegrator",
    "g",
    "run_fixture",
    "run",
]
