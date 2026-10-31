"""Synthetic checks for :mod:`original_q4_core_enclosure_adapter_v2`.

The fixture accepts an injected adapter class or factory.  Its interval class
and provider are local synthetic objects; no project provider, source graph,
archive, old implementation, or numerical source value is imported or
called.  ``run_fixture`` is intentionally not executed at import time.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Callable


SCHEMA = "original_q4_core_enclosure_adapter_fixture_v2"


class SyntheticInterval:
    """Fresh interval object used to ensure the adapter does not return a scalar."""

    __slots__ = ("label", "lower", "upper")

    def __init__(self, label: str, lower: Any, upper: Any):
        self.label = label
        self.lower = lower
        self.upper = upper


class _FakeProvider:
    def __init__(self, family: Any, *, mode: str = "complete"):
        self.family = family
        self.mode = mode
        self.calls: list[dict[str, Any]] = []

    def __call__(
        self,
        quantity: str,
        rho: Any,
        Z: Any,
        *,
        radial_order: int,
        axial_order: int,
        degree: int,
        depth: int,
    ) -> dict[str, Any]:
        self.calls.append(
            {
                "quantity": quantity,
                "rho": rho,
                "Z": Z,
                "radial_order": radial_order,
                "axial_order": axial_order,
                "degree": degree,
                "depth": depth,
            }
        )
        result: dict[str, Any] = {
            "quantity": quantity,
            "radial_order": radial_order,
            "axial_order": axial_order,
            "finite_polynomial_enclosure": SyntheticInterval(
                "finite", Fraction(1), Fraction(2)
            ),
            "infinite_radial_tail_upper": SyntheticInterval(
                "tail", Fraction(0), Fraction(1, 10)
            ),
            "enclosure": SyntheticInterval(
                "complete", Fraction(9, 10), Fraction(21, 10)
            ),
            "all_rows_are_ordinary_derivatives": True,
            "actual_twenty_term_nonlinear_solution_enclosed": True,
            "numerical_center_not_selected": True,
            "source_family": self.family,
            "radial_degree": degree,
            "axial_depth": depth,
        }
        if self.mode == "missing_ordinary":
            result.pop("all_rows_are_ordinary_derivatives")
        elif self.mode == "false_ordinary":
            result["all_rows_are_ordinary_derivatives"] = False
        elif self.mode == "missing_nonlinear":
            result.pop("actual_twenty_term_nonlinear_solution_enclosed")
        elif self.mode == "missing_center":
            result.pop("numerical_center_not_selected")
        elif self.mode in {
            "missing_finite",
            "missing_tail",
            "missing_enclosure",
        }:
            result.pop(
                {
                    "missing_finite": "finite_polynomial_enclosure",
                    "missing_tail": "infinite_radial_tail_upper",
                    "missing_enclosure": "enclosure",
                }[self.mode]
            )
        elif self.mode == "scalar_enclosure":
            result["enclosure"] = Fraction(3, 2)
        return result


def _factory(adapter: Any) -> Callable[..., Any] | None:
    """Accept a module, class, or direct factory without importing it."""

    for name in (
        "CoreEnclosureAdapter",
        "Q4CoreEnclosureAdapter",
        "ComponentEnclosureAdapter",
        "make_adapter",
        "build_adapter",
        "make_component_adapter",
    ):
        candidate = getattr(adapter, name, None)
        if callable(candidate):
            return candidate
    return adapter if callable(adapter) else None


def _new(
    factory: Callable[..., Any],
    provider: _FakeProvider,
    family: Any,
    descriptor: Any,
    evidence: Any,
    *,
    degree: int = 6,
    depth: int = 3,
    validator: Callable[[Any], Any] | None = None,
) -> Any:
    kwargs = {
        "provider": provider,
        "source_family": family,
        "provider_descriptor": descriptor,
        "degree": degree,
        "depth": depth,
        "interval_validator": validator,
        "source_evidence": evidence,
    }
    try:
        return factory(**kwargs)
    except TypeError:
        # Keep compatibility with a positional constructor while retaining
        # exactly the same injected values.
        return factory(
            provider,
            family,
            descriptor,
            degree,
            depth,
            validator,
            evidence,
        )


def _component(instance: Any) -> Callable[..., Any]:
    method = getattr(instance, "component", None)
    if callable(method):
        return method
    if callable(instance):
        return instance
    raise TypeError("injected adapter has no component callable")


def _record(family: Any, descriptor: Any, evidence: Any, **updates: Any) -> dict[str, Any]:
    value = {
        "quantity": "Phi",
        "radial_order": 1,
        "ordinary_slow_Z_derivative_order": 2,
        "source_family": family,
        "source_provider": descriptor,
        "source_evidence": evidence,
    }
    value.update(updates)
    return value


def _raises(call: Callable[[], Any]) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def run_fixture(adapter: Any) -> dict[str, Any]:
    """Run independent synthetic admission and forwarding checks."""

    factory = _factory(adapter)
    checks: dict[str, bool] = {}
    if factory is None:
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"adapter_factory_callable": False},
            "limitations": ["No injected adapter factory was supplied."],
        }

    family = {"name": "synthetic-core-family", "revision": 1}
    descriptor = {"name": "synthetic-component-enclosure", "revision": 1}
    evidence = {"source": "synthetic fixture only", "caller_asserted": True}

    provider = _FakeProvider(family)
    validator = lambda value: isinstance(value, SyntheticInterval) and value.label == "complete"
    instance = _new(factory, provider, family, descriptor, evidence, validator=validator)
    component = _component(instance)
    checks["provider_lazy_before_component"] = provider.calls == []
    result = component(_record(family, descriptor, evidence), Fraction(1), Fraction(2))
    checks["returns_complete_interval_object"] = (
        getattr(result, "kind", None) == "enclosure"
        and isinstance(getattr(result, "interval", None), SyntheticInterval)
        and result.interval.label == "complete"
        and getattr(result, "actual_value_selected", None) is False
    )
    checks["does_not_return_finite_or_tail"] = result.interval.label not in {"finite", "tail"}
    checks["rejects_scalar_conversion"] = _raises(lambda: float(result))
    checks["forwards_canonical_arguments"] = provider.calls == [
        {
            "quantity": "Phi",
            "rho": Fraction(1),
            "Z": Fraction(2),
            "radial_order": 1,
            "axial_order": 2,
            "degree": 6,
            "depth": 3,
        }
    ]
    receipt = getattr(instance, "last_receipt", {})
    checks["receipt_retains_uncertified_enclosure_semantics"] = (
        receipt.get("kind") == "enclosure"
        and receipt.get("global_status") == "uncertified"
        and receipt.get("actual_value_selected") is False
        and receipt.get("source_binding_independently_verified") is False
        and receipt.get("errors_independently_verified") is False
        and receipt.get("caller_claims", {}).get("interval_validator_is_proof") is False
    )

    def rejects(record: dict[str, Any], **options: Any) -> bool:
        local_provider = _FakeProvider(family, mode=options.pop("mode", "complete"))
        local = _new(
            factory,
            local_provider,
            family,
            descriptor,
            evidence,
            validator=options.pop("validator", validator),
            degree=options.pop("degree", 6),
            depth=options.pop("depth", 3),
        )
        return _raises(lambda: _component(local)(record, Fraction(1), Fraction(2)))

    checks["rejects_unsupported_quantity"] = rejects(
        _record(family, descriptor, evidence, quantity="not-a-core-component")
    )
    checks["rejects_derivative_sum_above_five"] = rejects(
        _record(family, descriptor, evidence, radial_order=4, ordinary_slow_Z_derivative_order=2)
    )
    checks["rejects_axial_order_above_depth"] = rejects(
        _record(family, descriptor, evidence, ordinary_slow_Z_derivative_order=4)
    )
    checks["rejects_degree_below_six"] = rejects(
        _record(family, descriptor, evidence), degree=5
    )
    checks["rejects_family_mismatch"] = rejects(
        _record({"name": "other-family"}, descriptor, evidence)
    )
    checks["rejects_provider_descriptor_mismatch"] = rejects(
        _record(family, {"name": "other-provider"}, evidence)
    )
    checks["rejects_source_evidence_mismatch"] = rejects(
        _record(family, descriptor, {"source": "wrong"})
    )
    checks["rejects_missing_ordinary_derivative_claim"] = rejects(
        _record(family, descriptor, evidence), mode="missing_ordinary"
    )
    checks["rejects_missing_twenty_term_claim"] = rejects(
        _record(family, descriptor, evidence), mode="missing_nonlinear"
    )
    checks["rejects_missing_center_claim"] = rejects(
        _record(family, descriptor, evidence), mode="missing_center"
    )
    checks["rejects_missing_finite_polynomial"] = rejects(
        _record(family, descriptor, evidence), mode="missing_finite"
    )
    checks["rejects_missing_tail"] = rejects(
        _record(family, descriptor, evidence), mode="missing_tail"
    )
    checks["rejects_missing_complete_enclosure"] = rejects(
        _record(family, descriptor, evidence), mode="missing_enclosure"
    )
    checks["rejects_unverified_enclosure"] = rejects(
        _record(family, descriptor, evidence), validator=lambda _value: False
    )
    checks["rejects_scalar_point_output"] = rejects(
        _record(family, descriptor, evidence), mode="scalar_enclosure"
    )

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "limitations": [
            "All source families, provider descriptors, evidence, intervals, and values are synthetic.",
            "The interval callback is caller validation, not an independently verified proof.",
            "No archived graph, old provider, finite polynomial, tail cap, or midpoint was used.",
        ],
    }


__all__ = ["SCHEMA", "SyntheticInterval", "run_fixture"]

