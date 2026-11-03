"""Synthetic, JSON-safe checks for :mod:`n1_scoped_graph_evaluator_v1`.

The fixture uses only tiny in-memory graphs.  It intentionally does not load
the real archive or any prior runtime.  Its integral callback is a small exact
mock for the named fixture rules; it checks that the evaluator passes scalar
endpoints and a body closure, not endpoint node ids.
"""

from __future__ import annotations

from collections.abc import Mapping
from fractions import Fraction
from typing import Any


def _json_safe(value: Any) -> Any:
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


class _ExactBackend:
    """Exact Fraction backend plus deliberately small fixture integrals."""

    def __init__(self) -> None:
        self.integral_calls: list[dict[str, Any]] = []
        self.integral_metadata_clean = True

    def validate_scalar(self, value: Any) -> bool:
        return isinstance(value, Fraction)

    def rational(self, numerator: int, denominator: int) -> Fraction:
        return Fraction(numerator, denominator)

    def add(self, left: Fraction, right: Fraction) -> Fraction:
        return left + right

    def mul(self, left: Fraction, right: Fraction) -> Fraction:
        return left * right

    def is_positive(self, value: Fraction) -> bool:
        return value > 0

    def divide(self, numerator: Fraction, denominator: Fraction) -> Fraction:
        return numerator / denominator

    def unary(self, function: str, value: Fraction) -> Fraction:
        if function == "abs":
            return abs(value)
        raise ValueError("fixture unary function is not installed: " + function)

    def integral(
        self,
        lower: Fraction,
        upper: Fraction,
        body: Any,
        metadata: Mapping[str, Any],
    ) -> Fraction:
        if not isinstance(lower, Fraction) or not isinstance(upper, Fraction):
            raise TypeError("integral endpoints must be scalar Fraction values")
        if not callable(body):
            raise TypeError("integral body must be a callable closure")
        forbidden = {"integrand", "lower", "upper", "variable"}
        if any(key in metadata for key in forbidden):
            self.integral_metadata_clean = False
            raise ValueError("endpoint ids leaked into integral metadata")
        self.integral_calls.append(
            {"lower": lower, "upper": upper, "body": callable(body), "metadata": dict(metadata)}
        )
        rule = metadata.get("fixture_rule")
        if rule == "identity":
            at_zero = body(Fraction(0))
            at_one = body(Fraction(1))
            if at_zero != Fraction(0) or at_one != Fraction(1):
                raise ValueError("identity body scope was not installed")
            return (upper * upper - lower * lower) / 2
        if rule == "constant":
            at_zero = body(Fraction(0))
            if body(Fraction(1)) != at_zero:
                raise ValueError("constant body changed with the binder")
            return at_zero * (upper - lower)
        if rule == "outer_plus_identity":
            at_zero = body(Fraction(0))
            at_one = body(Fraction(1))
            if at_one - at_zero != Fraction(1):
                raise ValueError("outer-plus-identity body scope was not installed")
            return at_zero * (upper - lower) + (upper * upper - lower * lower) / 2
        raise ValueError("unknown fixture integral rule")


def _rejected(call: Any) -> bool:
    try:
        call()
    except Exception:
        return True
    return False


def _substitution_value_is_outer(module: Any) -> tuple[bool, Any]:
    # value = x + y is evaluated before x is rebound.  With outer x=7,y=2,
    # the body x must receive 9 rather than recursively seeing its own value.
    graph = [
        {"operation": "exact_rational", "numerator": 7},  # 0
        {"operation": "exact_rational", "numerator": 2},  # 1
        {"operation": "bound_variable", "name": "x"},  # 2
        {"operation": "bound_variable", "name": "y"},  # 3
        {"operation": "sum", "arguments": [2, 3]},  # 4: value x+y
        {"operation": "function_substitution", "expression": 2, "variable": 2, "value": 4},  # 5
    ]
    backend = _ExactBackend()
    result = module.evaluate(
        graph,
        5,
        backend,
        free_variables={2: Fraction(7), 3: Fraction(2)},
    )
    return result["value"] == Fraction(9), result["value"]


def _binder_shadow_is_by_id(module: Any) -> tuple[bool, Any]:
    # Two distinct bound-variable ids intentionally share the name t.  The
    # inner substitution binds id 3, while its body references outer id 2.
    graph = [
        {"operation": "exact_rational", "numerator": 3},  # 0
        {"operation": "exact_rational", "numerator": 7},  # 1
        {"operation": "bound_variable", "name": "t"},  # 2 outer
        {"operation": "bound_variable", "name": "t"},  # 3 inner
        {"operation": "function_substitution", "expression": 2, "variable": 3, "value": 1},  # 4
        {"operation": "function_substitution", "expression": 4, "variable": 2, "value": 0},  # 5
    ]
    result = module.evaluate(graph, 5, _ExactBackend())
    return result["value"] == Fraction(3), result["value"]


def _integral_scope_is_installed(module: Any) -> tuple[bool, Any, _ExactBackend]:
    graph = [
        {"operation": "exact_rational", "numerator": 0},  # 0
        {"operation": "exact_rational", "numerator": 2},  # 1
        {"operation": "bound_variable", "name": "t"},  # 2
        {
            "operation": "definite_integral",
            "integrand": 2,
            "variable": "t",
            "lower": 0,
            "upper": 1,
            "fixture_rule": "identity",
        },  # 3
    ]
    backend = _ExactBackend()
    result = module.evaluate(graph, 3, backend, integral_binders={3: 2})
    passed = (
        result["value"] == Fraction(2)
        and len(backend.integral_calls) == 1
        and backend.integral_calls[0]["body"] is True
        and backend.integral_metadata_clean
    )
    return passed, result["value"], backend


def _integral_outer_endpoints_inner_frame(module: Any) -> tuple[bool, Any]:
    # The endpoint nodes are outer bound variables.  The integrand uses the
    # same bound-variable id as the sidecar, but must see a fresh inner frame;
    # an outer free value for that id is deliberately present as a trap.
    graph = [
        {"operation": "exact_rational", "numerator": 1},  # 0
        {"operation": "exact_rational", "numerator": 3},  # 1
        {"operation": "bound_variable", "name": "lower"},  # 2
        {"operation": "bound_variable", "name": "upper"},  # 3
        {"operation": "bound_variable", "name": "t"},  # 4
        {
            "operation": "definite_integral",
            "integrand": 4,
            "variable": "t",
            "lower": 2,
            "upper": 3,
            "fixture_rule": "identity",
        },  # 5
    ]
    backend = _ExactBackend()
    result = module.evaluate(
        graph,
        5,
        backend,
        free_variables={2: Fraction(1), 3: Fraction(3), 4: Fraction(99)},
        integral_binders={5: 4},
    )
    passed = (
        result["value"] == Fraction(4)
        and backend.integral_calls[0]["lower"] == Fraction(1)
        and backend.integral_calls[0]["upper"] == Fraction(3)
    )
    return passed, result["value"]


def _same_node_can_use_two_envs(module: Any) -> tuple[bool, Any]:
    graph = [
        {"operation": "exact_rational", "numerator": 1},  # 0
        {"operation": "exact_rational", "numerator": 2},  # 1
        {"operation": "bound_variable", "name": "x"},  # 2 shared body
        {"operation": "function_substitution", "expression": 2, "variable": 2, "value": 0},  # 3
        {"operation": "function_substitution", "expression": 2, "variable": 2, "value": 1},  # 4
        {"operation": "sum", "arguments": [3, 4]},  # 5
    ]
    result = module.evaluate(graph, 5, _ExactBackend())
    return result["value"] == Fraction(3), result["value"]


def _integral_sidecar_is_explicit(module: Any) -> dict[str, bool]:
    graph = [
        {"operation": "exact_rational", "numerator": 0},  # 0
        {"operation": "exact_rational", "numerator": 1},  # 1
        {"operation": "bound_variable", "name": "t"},  # 2
        {"operation": "bound_variable", "name": "t"},  # 3, duplicate name
        {"operation": "bound_variable", "name": "s"},  # 4, mismatched name
        {
            "operation": "definite_integral",
            "integrand": 2,
            "variable": "t",
            "lower": 0,
            "upper": 1,
            # With the wrong sidecar, an outer free value for id 2 would make
            # this body look constant and could silently pass old code.
            "fixture_rule": "constant",
        },  # 5
    ]
    backend = _ExactBackend()
    missing = _rejected(lambda: module.evaluate(graph, 5, backend))
    ambiguous_explicit = _rejected(
        lambda: module.evaluate(graph, 5, _ExactBackend(), integral_binders={5: 2})
    )
    wrong_same_name_id = _rejected(
        lambda: module.evaluate(
            graph,
            5,
            _ExactBackend(),
            integral_binders={5: 3},
            free_variables={2: Fraction(9)},
        )
    )
    wrong_name = _rejected(
        lambda: module.evaluate(graph, 5, _ExactBackend(), integral_binders={5: 4})
    )
    return {
        "missing_sidecar_rejected": missing,
        "same_name_ambiguity_rejected": ambiguous_explicit,
        "same_name_wrong_id_rejected": wrong_same_name_id,
        "same_name_wrong_id_with_free_value_rejected": wrong_same_name_id,
        "mismatched_name_rejected": wrong_name,
    }


def _opaque_and_positivity_fail_closed(module: Any) -> dict[str, bool]:
    bad_division_graph = [
        {"operation": "exact_rational", "numerator": 1},  # 0
        {"operation": "exact_rational", "numerator": -1},  # 1
        {"operation": "positive_quotient", "numerator": 0, "denominator": 1},  # 2
    ]
    positivity_rejected = _rejected(
        lambda: module.evaluate(bad_division_graph, 2, _ExactBackend())
    )

    unknown_graph = [
        {
            "operation": "future_opaque_operation",
            "payload": {"argument": 999, "binder": {"children": [998]}},
        }
    ]
    calls: list[tuple[int, str, Any]] = []

    def unknown_oracle(node_id: int, record: Any, category: str, env: Any) -> Fraction:
        calls.append((node_id, category, record))
        return Fraction(4)

    oracle_result = module.evaluate(
        unknown_graph,
        0,
        _ExactBackend(),
        oracle={"unknown": unknown_oracle},
    )
    oracle_success = (
        oracle_result["value"] == Fraction(4)
        and calls == [(0, "unknown", unknown_graph[0])]
        and oracle_result["opaque_nodes_used"] == [0]
    )
    oracle_missing = _rejected(lambda: module.evaluate(unknown_graph, 0, _ExactBackend()))
    variable_missing = _rejected(
        lambda: module.evaluate(
            [{"operation": "bound_variable", "name": "free"}],
            0,
            _ExactBackend(),
        )
    )
    validator_domain_rejected = _rejected(
        lambda: module.evaluate(
            [{"operation": "bound_variable", "name": "fraction_only"}],
            0,
            _ExactBackend(),
            free_variables={0: 0.5},
        )
    )
    return {
        "positive_denominator_required": positivity_rejected,
        "category_oracle_success": oracle_success,
        "missing_unknown_oracle_rejected": oracle_missing,
        "missing_free_variable_rejected": variable_missing,
        "backend_validator_rejects_real_outside_domain": validator_domain_rejected,
    }


def _cycle_rejected(module: Any) -> bool:
    cycle_graph = [{"operation": "sum", "arguments": [0]}]
    return _rejected(lambda: module.evaluate(cycle_graph, 0, _ExactBackend()))


def run_fixture(module: Any) -> dict[str, Any]:
    """Run synthetic scope checks and return only JSON-safe data."""

    checks: dict[str, bool] = {}
    selected: dict[str, Any] = {}
    evaluate = getattr(module, "evaluate", None)
    execute = getattr(module, "execute", None)
    checks["evaluator_api"] = callable(evaluate)
    checks["execute_api"] = callable(execute)
    if not checks["evaluator_api"]:
        return {
            "schema": "n1_scoped_graph_evaluator_fixture_v1",
            "passed": False,
            "checks": checks,
            "selected": selected,
            "json_safe": True,
            "synthetic_only": True,
            "used_old_results": False,
        }

    try:
        substitution_ok, substitution_value = _substitution_value_is_outer(module)
        checks["substitution_value_uses_outer_env"] = substitution_ok
        selected["substitution_value"] = substitution_value

        shadow_ok, shadow_value = _binder_shadow_is_by_id(module)
        checks["same_name_binders_are_id_scoped"] = shadow_ok
        selected["shadow_value"] = shadow_value

        integral_ok, integral_value, integral_backend = _integral_scope_is_installed(module)
        checks["integral_body_uses_fresh_scope"] = integral_ok
        checks["integral_backend_received_closure"] = bool(integral_backend.integral_calls)
        selected["integral_value"] = integral_value

        endpoint_ok, endpoint_value = _integral_outer_endpoints_inner_frame(module)
        checks["integral_outer_endpoints_inner_frame"] = endpoint_ok
        selected["integral_endpoint_scope_value"] = endpoint_value

        same_node_ok, same_node_value = _same_node_can_use_two_envs(module)
        checks["same_node_different_envs"] = same_node_ok
        selected["same_node_sum"] = same_node_value

        checks.update(_integral_sidecar_is_explicit(module))
        checks.update(_opaque_and_positivity_fail_closed(module))
        checks["cycles_rejected"] = _cycle_rejected(module)

        result_probe = module.evaluate(
            [{"operation": "exact_rational", "numerator": 1}],
            0,
            _ExactBackend(),
        )
        claims = result_probe.get("claims", {})
        provenance = result_probe.get("provenance", {})
        checks["conditional_not_certified"] = (
            result_probe.get("status") == "conditional"
            and result_probe.get("conditional") is True
            and result_probe.get("certified") is False
            and result_probe.get("claims", {}).get("certified") is False
        )
        checks["no_physical_claims"] = all(
            claims.get(name) is False
            for name in ("function_equality", "physical_m15", "Kchi", "Navier_Stokes")
        )
        checks["scope_capabilities_explicit"] = (
            provenance.get("lexical_binder_supported") is True
            and provenance.get("integral_runtime_supported") is True
            and provenance.get("substitution_runtime_supported") is True
        )
    except Exception:
        for name in (
            "substitution_value_uses_outer_env",
            "same_name_binders_are_id_scoped",
            "integral_body_uses_fresh_scope",
            "integral_backend_received_closure",
            "integral_outer_endpoints_inner_frame",
            "same_node_different_envs",
            "missing_sidecar_rejected",
            "same_name_ambiguity_rejected",
            "same_name_wrong_id_rejected",
            "same_name_wrong_id_with_free_value_rejected",
            "mismatched_name_rejected",
            "positive_denominator_required",
            "category_oracle_success",
            "missing_unknown_oracle_rejected",
            "missing_free_variable_rejected",
            "backend_validator_rejects_real_outside_domain",
            "cycles_rejected",
            "conditional_not_certified",
            "no_physical_claims",
            "scope_capabilities_explicit",
        ):
            checks.setdefault(name, False)

    checks["json_safe"] = True
    checks["synthetic_only"] = True
    return {
        "schema": "n1_scoped_graph_evaluator_fixture_v1",
        "passed": all(checks.values()),
        "checks": checks,
        "selected": _json_safe(selected),
        "json_safe": True,
        "synthetic_only": True,
        "used_old_results": False,
    }


__all__ = ["run_fixture"]
