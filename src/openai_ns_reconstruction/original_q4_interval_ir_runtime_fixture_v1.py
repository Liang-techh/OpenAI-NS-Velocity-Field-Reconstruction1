"""Fresh synthetic fixture for ``original_q4_interval_ir_runtime_v1``.

The fixture is self-contained and imports no prior kernel, evaluator,
adapter, lowering helper, or archived graph.  The runner injects the new
runtime callable into :func:`run_fixture`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Any


_SCHEMA = "original_q4_interval_ir_runtime_fixture_v1"
_FAMILY = "synthetic_fresh_runtime_family_20261031"
_GRAPH_HASH = "b" * 64


def _fraction(value: Any) -> Fraction:
    if type(value) is int:
        return Fraction(value, 1)
    if isinstance(value, Fraction):
        return Fraction(value.numerator, value.denominator)
    raise TypeError("fixture endpoints must be exact")


@dataclass(frozen=True, slots=True)
class FreshInterval:
    """Exact rational enclosure used only by this fixture."""

    lo: Fraction | int
    hi: Fraction | int
    kind: str = field(default="enclosure", init=False)

    def __post_init__(self) -> None:
        lower = _fraction(self.lo)
        upper = _fraction(self.hi)
        if lower > upper:
            raise ValueError("fixture interval is reversed")
        object.__setattr__(self, "lo", lower)
        object.__setattr__(self, "hi", upper)

    @property
    def width(self) -> Fraction:
        return self.hi - self.lo

    def __add__(self, other: Any) -> "FreshInterval":
        rhs = other if isinstance(other, FreshInterval) else FreshInterval(other, other)
        return FreshInterval(self.lo + rhs.lo, self.hi + rhs.hi)

    def __radd__(self, other: Any) -> "FreshInterval":
        return self + other

    def __sub__(self, other: Any) -> "FreshInterval":
        rhs = other if isinstance(other, FreshInterval) else FreshInterval(other, other)
        return FreshInterval(self.lo - rhs.hi, self.hi - rhs.lo)

    def __rsub__(self, other: Any) -> "FreshInterval":
        return FreshInterval(other, other) - self

    def __mul__(self, other: Any) -> "FreshInterval":
        rhs = other if isinstance(other, FreshInterval) else FreshInterval(other, other)
        values = (
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        )
        return FreshInterval(min(values), max(values))

    def __rmul__(self, other: Any) -> "FreshInterval":
        return self * other

    def __truediv__(self, other: Any) -> "FreshInterval":
        rhs = other if isinstance(other, FreshInterval) else FreshInterval(other, other)
        if rhs.lo <= 0:
            raise ZeroDivisionError("fixture denominator is not positive")
        return self * FreshInterval(1 / rhs.hi, 1 / rhs.lo)

    def __neg__(self) -> "FreshInterval":
        return FreshInterval(-self.hi, -self.lo)

    def __abs__(self) -> "FreshInterval":
        if self.lo >= 0:
            return self
        if self.hi <= 0:
            return FreshInterval(-self.hi, -self.lo)
        return FreshInterval(0, max(-self.lo, self.hi))

    def __float__(self) -> float:
        raise TypeError("fixture interval cannot become a point")

    def __bool__(self) -> bool:
        raise TypeError("fixture interval has no scalar truth value")


@dataclass(frozen=True, slots=True)
class FreshComponentEnclosure:
    interval: FreshInterval
    kind: str = field(default="enclosure", init=False)
    actual_value_selected: bool = field(default=False, init=False)
    certified_by_adapter: bool = field(default=False, init=False)


@dataclass(frozen=True, slots=True)
class FreshComponentExtraction:
    interval: FreshInterval
    kind: str = field(default="conditional_enclosure", init=False)
    endpoint_conversion_verified: bool = field(default=False, init=False)
    source_identity_verified: bool = field(default=False, init=False)
    actual_value_selected: bool = field(default=False, init=False)

    def receipt(self) -> dict[str, Any]:
        return {
            "schema": _SCHEMA,
            "kind": self.kind,
            "interval_type": "FreshInterval",
            "endpoint_conversion_verified": False,
            "source_identity_verified": False,
            "actual_value_selected": False,
        }


@dataclass(frozen=True, slots=True)
class FreshIntegralEnclosure:
    interval: FreshInterval
    coverage_cells: tuple[FreshInterval, ...]
    kind: str = field(default="conditional_enclosure", init=False)
    actual_value_selected: bool = field(default=False, init=False)
    provider_certified: bool = field(default=False, init=False)

    def receipt(self) -> dict[str, Any]:
        return {
            "schema": _SCHEMA,
            "kind": self.kind,
            "coverage_count": len(self.coverage_cells),
            "provider_certified": False,
            "actual_value_selected": False,
        }


class FreshBackend:
    """Exact interval arithmetic with no scalar conversion path."""

    def rational(self, numerator: int, denominator: int = 1) -> FreshInterval:
        if type(numerator) is not int or type(denominator) is not int:
            raise TypeError("fixture rational must be integer based")
        if denominator <= 0:
            raise ValueError("fixture denominator must be positive")
        value = Fraction(numerator, denominator)
        return FreshInterval(value, value)

    def add(self, left: FreshInterval, right: FreshInterval) -> FreshInterval:
        return left + right

    def multiply(self, left: FreshInterval, right: FreshInterval) -> FreshInterval:
        return left * right

    def negative(self, value: FreshInterval) -> FreshInterval:
        return -value

    def divide(self, numerator: FreshInterval, denominator: FreshInterval) -> FreshInterval:
        return numerator / denominator

    def unary(self, name: str, value: FreshInterval) -> FreshInterval:
        if name == "abs":
            return abs(value)
        if name == "exp":
            raise ValueError("fixture exp is intentionally not approximated")
        raise ValueError("fixture unary is unsupported")

    def is_positive(self, value: FreshInterval) -> bool:
        return value.lo > 0


class FreshRangeIntegrator:
    """Range integrator that passes a complete cell to every callback."""

    def integrate(
        self,
        callback: Any,
        lower: FreshInterval,
        upper: FreshInterval,
        partitions: int,
    ) -> FreshIntegralEnclosure:
        if type(partitions) is not int or partitions <= 0:
            raise ValueError("fixture partitions must be positive")
        if not isinstance(lower, FreshInterval) or not isinstance(upper, FreshInterval):
            raise TypeError("fixture bounds must be intervals")
        left = lower.lo
        right = upper.hi
        if left > right:
            raise ValueError("fixture bounds are reversed")
        width = (right - left) / partitions
        total = FreshInterval(0, 0)
        cells: list[FreshInterval] = []
        for index in range(partitions):
            cell = FreshInterval(left + width * index, left + width * (index + 1))
            value = callback(cell)
            if type(value) is not FreshInterval or value.kind != "enclosure":
                raise TypeError("fixture callback must return an interval")
            cells.append(cell)
            total = total + value * width
        return FreshIntegralEnclosure(total, tuple(cells))


class FreshOracle:
    mode = "synthetic_reference"

    def __init__(self, integrator: FreshRangeIntegrator | None = None) -> None:
        self.integrator = integrator or FreshRangeIntegrator()
        self.parameter_calls: list[dict[str, Any]] = []
        self.component_calls: list[dict[str, Any]] = []

    def parameter(self, row: dict[str, Any]) -> FreshInterval:
        self.parameter_calls.append(dict(row))
        return FreshInterval(1, 1)

    def component(
        self, row: dict[str, Any], rho: FreshInterval, z_value: FreshInterval
    ) -> FreshComponentEnclosure:
        if not isinstance(rho, FreshInterval) or not isinstance(z_value, FreshInterval):
            raise TypeError("fixture component arguments must be intervals")
        self.component_calls.append({"row": dict(row), "rho": rho, "Z": z_value})
        return FreshComponentEnclosure(rho + z_value)

    def integrate(
        self,
        callback: Any,
        lower: FreshInterval,
        upper: FreshInterval,
        partitions: int,
    ) -> FreshIntegralEnclosure:
        return self.integrator.integrate(callback, lower, upper, partitions)


def extract_component(
    component: FreshComponentEnclosure, endpoint_converter: Any
) -> FreshComponentExtraction:
    if component.kind != "enclosure" or component.actual_value_selected is not False:
        raise TypeError("fixture component tag is invalid")
    raw = endpoint_converter(component.interval)
    if type(raw) is not FreshInterval:
        raise TypeError("fixture endpoint conversion must retain exact interval")
    return FreshComponentExtraction(raw)


def endpoint_converter(raw: FreshInterval) -> FreshInterval:
    if type(raw) is not FreshInterval:
        raise TypeError("fixture endpoint is not exact")
    return raw


def _rational(numerator: int, denominator: int = 1, *, ref: str | None = None) -> dict[str, Any]:
    row: dict[str, Any] = {
        "op": "exact_rational",
        "numerator": numerator,
        "denominator": denominator,
        "source_row": {"operation": "exact_rational", "numerator": numerator, "denominator": denominator},
    }
    if ref is not None:
        row["source_ref"] = ref
    return row


def _fresh(token: str) -> dict[str, Any]:
    return {"kind": "fresh_binder", "token": token, "source_name": "same_name"}


def build_synthetic_ir() -> dict[str, Any]:
    """Build a new expression with same-name, distinct fresh tokens."""

    outer_token = "fresh_runtime_outer_20261031"
    inner_token = "fresh_runtime_inner_20261031"
    outer_ref = _fresh(outer_token)
    inner_ref = _fresh(inner_token)
    parameter = {
        "op": "current_original_source_parameter",
        "source_ref": "parameter_runtime_20261031",
        "source_row": {
            "operation": "current_original_source_parameter",
            "name": "fresh_parameter_20261031",
            "source_family": _FAMILY,
        },
    }
    component = {
        "op": "current_original_unique_analytic_core_component",
        "source_ref": "component_runtime_20261031",
        "source_row": {
            "operation": "current_original_unique_analytic_core_component",
            "source_family": _FAMILY,
            "quantity": "synthetic_component",
        },
        "arguments": {
            "rho": {"op": "sum", "terms": [outer_ref, inner_ref]},
            "Z_variable": parameter,
        },
    }
    inner_lower = {
        "op": "product",
        "factors": [_rational(0, ref="zero_inner_bound_runtime_20261031"), _fresh(outer_token)],
    }
    inner_integrand = {
        "op": "product",
        "factors": [
            {
                "op": "positive_quotient",
                "numerator": component,
                "denominator": _rational(1, ref="one_quotient_runtime_20261031"),
            },
            {
                "op": "analytic_unary",
                "function": "abs",
                "argument": {"op": "negative", "argument": _rational(1, ref="minus_one_runtime_20261031")},
            },
        ],
    }
    inner_integral = {
        "op": "definite_integral",
        "source_ref": "integral_inner_runtime_20261031",
        "binder": {"kind": "fresh_binder", "token": inner_token, "source_name": "same_name"},
        "lower": inner_lower,
        "upper": _rational(1, ref="one_inner_bound_runtime_20261031"),
        "integrand": inner_integrand,
    }
    return {
        "op": "definite_integral",
        "source_ref": "integral_outer_runtime_20261031",
        "binder": {"kind": "fresh_binder", "token": outer_token, "source_name": "same_name"},
        "lower": _rational(0, ref="zero_outer_bound_runtime_20261031"),
        "upper": _rational(1, ref="one_outer_bound_runtime_20261031"),
        "integrand": {"op": "sum", "terms": [parameter, inner_integral]},
    }


class _RawScalarOracle(FreshOracle):
    def parameter(self, row: dict[str, Any]) -> Any:
        return Fraction(1, 1)


class _RawMidpointIntegrator(FreshRangeIntegrator):
    def integrate(
        self,
        callback: Any,
        lower: FreshInterval,
        upper: FreshInterval,
        partitions: int,
    ) -> Any:
        return callback(Fraction(1, 2))


class _OriginalModeOracle(FreshOracle):
    mode = "original_source"
    source_family = _FAMILY
    graph_hash = _GRAPH_HASH


def make_fixture_inputs() -> dict[str, Any]:
    return {
        "oracle": FreshOracle(),
        "backend": FreshBackend(),
        "interval_type": FreshInterval,
        "extract_component": extract_component,
        "endpoint_converter": endpoint_converter,
        "partitions": 8,
        "free_variables": {},
        "admission": {"mode": "synthetic_reference"},
    }


def run_fixture(evaluate_interval: Any) -> dict[str, Any]:
    """Run all checks through an injected new-runtime callable."""

    if not callable(evaluate_interval):
        raise TypeError("the new interval runtime callable is required")
    expression = build_synthetic_ir()
    inputs = make_fixture_inputs()
    result = evaluate_interval(expression, **inputs)
    # Integral_0^1 [1 + Integral_0^1 (outer + inner + 1) d(inner)] d(outer) = 3.
    exact_value = Fraction(3, 1)
    interval = result.get("interval") if isinstance(result, dict) else None
    checks: dict[str, bool] = {}
    checks["conditional_status"] = (
        isinstance(result, dict)
        and result.get("status") == "conditional_enclosure"
        and result.get("actual_value_selected") is False
        and result.get("value") is None
        and result.get("certified") is False
    )
    checks["exact_value_contained"] = (
        type(interval) is FreshInterval
        and interval.lo <= exact_value <= interval.hi
    )
    components = result.get("component_extractions", []) if isinstance(result, dict) else []
    checks["complete_component_receipt"] = bool(components) and all(
        item.get("component_kind") == "enclosure"
        and item.get("extraction_kind") == "conditional_enclosure"
        and type(item.get("extracted_interval")) is FreshInterval
        and item.get("actual_value_selected") is False
        for item in components
    )
    integrals = result.get("integral_receipts", []) if isinstance(result, dict) else []
    checks["nested_same_name_tokens"] = (
        len(integrals) == 9
        and len({item.get("token") for item in integrals}) == 2
        and {item.get("source_name") for item in integrals} == {"same_name"}
    )
    checks["coverage_retained"] = bool(integrals) and all(
        item.get("coverage_count") == len(item.get("coverage_cells", ()))
        and all(type(cell) is FreshInterval for cell in item.get("coverage_cells", ()))
        for item in integrals
    )
    finer_inputs = dict(inputs)
    finer_inputs["partitions"] = 16
    finer_result = evaluate_interval(expression, **finer_inputs)
    finer_interval = finer_result.get("interval") if isinstance(finer_result, dict) else None
    checks["convergence_interval"] = (
        finer_result.get("status") == "conditional_enclosure"
        and type(finer_interval) is FreshInterval
        and finer_interval.lo <= exact_value <= finer_interval.hi
        and type(interval) is FreshInterval
        and finer_interval.width <= interval.width
    )

    bad_scalar = dict(inputs)
    bad_scalar["oracle"] = _RawScalarOracle()
    raw_scalar_result = evaluate_interval(expression, **bad_scalar)
    checks["raw_scalar_rejected"] = raw_scalar_result.get("status") == "rejected"

    bad_midpoint = dict(inputs)
    bad_midpoint["oracle"] = FreshOracle(_RawMidpointIntegrator())
    raw_midpoint_result = evaluate_interval(expression, **bad_midpoint)
    checks["raw_midpoint_rejected"] = raw_midpoint_result.get("status") == "rejected"

    original_inputs = dict(inputs)
    original_inputs["oracle"] = _OriginalModeOracle()
    original_inputs["admission"] = {
        "mode": "original_source",
        "source_family": _FAMILY,
        "graph_hash": _GRAPH_HASH,
        "free_variable_namespace": [],
    }
    original_expression = {
        "status": "lowered",
        "ir": expression,
        "source_family": _FAMILY,
        "graph_hash": _GRAPH_HASH,
    }
    original_result = evaluate_interval(original_expression, **original_inputs)
    checks["original_missing_binder_evidence_rejected"] = (
        original_result.get("status") == "rejected"
        and original_result.get("value") is None
    )

    checks["all"] = all(checks.values())
    return {
        "schema": _SCHEMA,
        "status": "passed" if checks["all"] else "failed",
        "checks": checks,
        "successful_result": result,
        "finer_result": finer_result,
        "raw_scalar_result": raw_scalar_result,
        "raw_midpoint_result": raw_midpoint_result,
        "original_result": original_result,
        "independent_exact_value": exact_value,
    }


run = run_fixture


__all__ = [
    "FreshInterval",
    "FreshComponentEnclosure",
    "FreshComponentExtraction",
    "FreshIntegralEnclosure",
    "FreshBackend",
    "FreshRangeIntegrator",
    "FreshOracle",
    "build_synthetic_ir",
    "make_fixture_inputs",
    "extract_component",
    "endpoint_converter",
    "run_fixture",
    "run",
]
