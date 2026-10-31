"""Independent synthetic records for ``original_q4_saved_moment_join_v1``.

The fixture constructs fresh JSON-shaped graph snapshots and fresh endpoint
tuples.  It never imports, executes, evaluates, or decodes a producer/archive.
The public return is intentionally just ``{passed, checks}`` so a runner can
record it without treating the fixture as source evidence.
"""

from __future__ import annotations

import copy
from typing import Any, Callable


SCHEMA = "original_q4_saved_moment_join_fixture_v1"
FAMILY = {
    "actual_five_defect_family_sha256": "a" * 64,
    "implicit_source_sha256": "b" * 64,
    "datum_enclosure_sha256": "c" * 64,
}
MOMENT_ORDERS = {
    "M": (0, 1, 2),
    "A": (0, 1, 2),
    "B": (0, 1, 2, 3),
    "C": (0, 1, 2),
}
RHO_NAME = "original_n1_true_core_rho"
MOMENT_VARIABLE_NAME = "original_core_inlet_moment_t"


def _endpoint(sign: int, mantissa: int, exponent: int) -> list[int]:
    return [sign, mantissa, exponent, mantissa.bit_length()]


def _finite_row(name: str, order: int) -> dict[str, Any]:
    # Distinct exact tuples make accidental cross-order copying visible.
    lower = _endpoint(0, 1 + order, -4 - len(name))
    upper = _endpoint(0, 2 + order, -4 - len(name))
    return {
        "lower_exact_mpf_tuple": lower,
        "upper_exact_mpf_tuple": upper,
        "source_row_tag": "fresh-finite-" + name + "-Z" + str(order),
    }


def _tail_row(name: str, order: int) -> dict[str, Any]:
    upper = _endpoint(0, 1 + order, -9 - len(name))
    return {
        "lower_exact_mpf_tuple": _endpoint(0, 0, 0),
        "upper_exact_mpf_tuple": upper,
        "source_row_tag": "fresh-tail-" + name + "-Z" + str(order),
    }


def _ordinary_row(name: str, order: int) -> dict[str, Any]:
    return {
        "lower_exact_mpf_tuple": _endpoint(0, 1, -8 - order),
        "upper_exact_mpf_tuple": _endpoint(0, 2, -8 - order),
        "source_row_tag": "fresh-ordinary-" + name + "-Z" + str(order),
        "is_enclosure": True,
    }


def _graph() -> tuple[list[dict[str, Any]], dict[str, list[int]], dict[str, list[int]], list[dict[str, Any]]]:
    nodes: list[dict[str, Any]] = [
        {"operation": "bound_variable", "name": RHO_NAME, "fresh_handle_namespace": "rho-A"},
        {"operation": "exact_rational", "numerator": 0, "denominator": 1},
        {"operation": "exact_rational", "numerator": 1, "denominator": 1},
        {"operation": "exact_rational", "numerator": 4, "denominator": 1},
    ]
    moments: dict[str, list[int]] = {}
    inlet_moments: dict[str, list[int]] = {}
    inspections: list[dict[str, Any]] = []
    for name, orders in MOMENT_ORDERS.items():
        moments[name] = []
        inlet_moments[name] = []
        for order in orders:
            density = len(nodes)
            nodes.append({
                "operation": "fresh_moment_density",
                "moment_name": name,
                "ordinary_Z_order": order,
                "fresh_handle_namespace": "density-" + name + "-" + str(order),
            })
            integral = len(nodes)
            integral_row = {
                "operation": "definite_integral",
                "integrand": density,
                "variable": MOMENT_VARIABLE_NAME,
                "lower": 1,
                "upper": 2,
                "actual_core_normalized_moment": name,
                "ordinary_slow_Z_derivative_order": order,
                "exact_function_integral": True,
                "numerical_value_not_installed": True,
                "fresh_integral_namespace": "integral-" + name + "-" + str(order),
            }
            nodes.append(integral_row)
            substitution = len(nodes)
            nodes.append({
                "operation": "function_substitution",
                "expression": integral,
                "variable": 0,
                "value": 3,
                "Z_independent_substitution": True,
                "fresh_substitution_namespace": "substitution-" + name + "-" + str(order),
            })
            moments[name].append(integral)
            inlet_moments[name].append(substitution)
            inspections.append({
                "frontier": {
                    "operation": "definite_integral",
                    "node_id": integral,
                    "integrand_node_id": density,
                    "binder_name": MOMENT_VARIABLE_NAME,
                },
                "inspection": {
                    "status": "observed",
                    "evidence": {
                        "integral_node": integral,
                        "integrand_node": density,
                        "binder_name": MOMENT_VARIABLE_NAME,
                        "unit_interval_bounds_verified": True,
                        "integral_source_row": copy.deepcopy(integral_row),
                    },
                    "scope": {
                        "capture_safe_substitution_proved": False,
                        "numeric_function_value_computed": False,
                        "numeric_field_value_computed": False,
                        "old_errors_removed": False,
                        "terminal_closure_proved": False,
                        "conditional_on_source_operator_semantics": True,
                        "integrand_scan_complete": False,
                    },
                },
            })
    return nodes, moments, inlet_moments, inspections


def _interior() -> dict[str, Any]:
    finite: dict[str, dict[str, Any]] = {}
    tails: dict[str, dict[str, Any]] = {}
    ordinary: dict[str, dict[str, Any]] = {}
    for name, orders in MOMENT_ORDERS.items():
        finite[name] = {"rho0_Z" + str(order): _finite_row(name, order) for order in orders}
        tails[name] = {"rho0_Z" + str(order): _tail_row(name, order) for order in orders}
        ordinary[name] = {"rho0_Z" + str(order): _ordinary_row(name, order) for order in orders}
    view = {
        "rho": {
            "lower_exact_mpf_tuple": _endpoint(0, 1, 2),
            "upper_exact_mpf_tuple": _endpoint(0, 1, 2),
        },
        "Z": {
            "lower_exact_mpf_tuple": _endpoint(1, 1, 0),
            "upper_exact_mpf_tuple": _endpoint(0, 1, 0),
        },
        "radial_degree": 24,
        "axial_depth": 6,
        "finite_polynomial_is_not_complete_solution": True,
        "source_enclosures_not_resolved_point_values": True,
        "finite_integrated_axial_Taylor_coefficient_grids": finite,
        "directed_integrated_axial_Taylor_tail_bounds": tails,
        "ordinary_rho4_axial6_moment_grids": ordinary,
    }
    return {
        **copy.deepcopy(FAMILY),
        "source_family": copy.deepcopy(FAMILY),
        "current_core_interior_views": {"inlet": view},
    }


def _physical(interior: dict[str, Any]) -> dict[str, Any]:
    view = interior["current_core_interior_views"]["inlet"]
    input_moments: dict[str, dict[str, Any]] = {}
    for name, orders in MOMENT_ORDERS.items():
        zones: dict[str, Any] = {}
        for order in (0, 1):
            finite = view["finite_integrated_axial_Taylor_coefficient_grids"][name]["rho0_Z" + str(order)]
            tail = view["directed_integrated_axial_Taylor_tail_bounds"][name]["rho0_Z" + str(order)]
            zones["z" + str(order)] = {
                "finite": [
                    copy.deepcopy(finite["lower_exact_mpf_tuple"]),
                    copy.deepcopy(finite["upper_exact_mpf_tuple"]),
                ],
                "error_upper": copy.deepcopy(tail["upper_exact_mpf_tuple"]),
            }
        input_moments[name] = zones
    return {
        "status": "returned",
        "schema": "original_core_physical_moment_binary_v2",
        "family_evidence": {
            "all_three_Q4_core_family_records_match": True,
            "fixed_point_family_matches": True,
            "saved_core_family": copy.deepcopy(FAMILY),
        },
        "views": {"inlet": {"input_moments": input_moments}},
    }


def _case() -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]], dict[str, Any], dict[str, Any]]:
    q4_nodes, moments, inlet_moments, inspections = _graph()
    core_nodes = copy.deepcopy(q4_nodes)
    core = {
        **copy.deepcopy(FAMILY),
        "source_family": copy.deepcopy(FAMILY),
        "exact_graph_nodes": core_nodes,
        "core_inlet_axial6_bindings": {
            "moments": copy.deepcopy(moments),
            "inlet_moments": copy.deepcopy(inlet_moments),
        },
    }
    q4 = {
        **copy.deepcopy(FAMILY),
        "source_family": copy.deepcopy(FAMILY),
        "exact_graph_nodes": q4_nodes,
    }
    interior = _interior()
    physical = _physical(interior)
    return q4, core, inspections, interior, physical


def _result(bridge: Callable[..., Any], case: tuple[Any, ...]) -> Any:
    try:
        # The fixture stores Q4 first for its mutation helpers; the public
        # bridge takes the source-core archive first.
        return bridge(case[1], case[0], *case[2:])
    except Exception:
        return None


def _unresolved(result: Any) -> bool:
    return isinstance(result, dict) and result.get("passed") is False and result.get("status") == "unresolved"


def run_fixture(bridge: Callable[..., Any]) -> dict[str, Any]:
    checks: dict[str, bool] = {}
    if not callable(bridge):
        return {"passed": False, "checks": {"bridge_callable": False}}
    checks["bridge_callable"] = True

    positive = _case()
    good = _result(bridge, positive)
    checks["positive_matched"] = isinstance(good, dict) and good.get("status") == "matched" and good.get("passed") is True
    checks["positive_13_rows"] = checks["positive_matched"] and len(good.get("evidence", {}).get("moment_rows", ())) == 13
    checks["positive_only_MABC"] = checks["positive_matched"] and {
        row.get("name") for row in good.get("evidence", {}).get("moment_rows", ())
    } == set(MOMENT_ORDERS)
    checks["positive_factorials_2_and_6"] = checks["positive_matched"] and (
        good.get("evidence", {}).get("factorial_conversion", {}).get("M", {}).get("2") == 2
        and good.get("evidence", {}).get("factorial_conversion", {}).get("B", {}).get("3") == 6
    )
    checks["positive_raw_rows_retained"] = checks["positive_matched"] and all(
        "saved_inlet_rows" in row
        and "finite_integrated_axial_Taylor_coefficient_row" in row["saved_inlet_rows"]
        and "directed_integrated_axial_Taylor_tail_bound_row" in row["saved_inlet_rows"]
        and "ordinary_rho4_axial6_moment_enclosure_row" in row["saved_inlet_rows"]
        for row in good.get("evidence", {}).get("moment_rows", ())
    )
    checks["positive_scope_does_not_claim_capture"] = checks["positive_matched"] and (
        good.get("scope", {}).get("capture_safe_substitution_proved") is False
        and good.get("scope", {}).get("dynamic_binder_identity_proved") is False
        and good.get("scope", {}).get("source_evaluator_invoked") is False
    )

    def reject(label: str, mutate: Callable[[tuple[Any, ...]], None]) -> None:
        case = list(_case())
        mutate(case)
        checks["reject_" + label] = _unresolved(_result(bridge, tuple(case)))

    def wrong_order(case: tuple[Any, ...]) -> None:
        q4, core, _inspections, _interior, _physical_result = case
        node_id = core["core_inlet_axial6_bindings"]["moments"]["B"][3]
        q4["exact_graph_nodes"][node_id]["ordinary_slow_Z_derivative_order"] = 2
        core["exact_graph_nodes"][node_id]["ordinary_slow_Z_derivative_order"] = 2

    def wrong_row(case: tuple[Any, ...]) -> None:
        _q4, core, _inspections, _interior, _physical_result = case
        node_id = core["core_inlet_axial6_bindings"]["moments"]["A"][1]
        core["exact_graph_nodes"][node_id]["actual_core_normalized_moment"] = "M"

    def wrong_endpoint(case: tuple[Any, ...]) -> None:
        q4, _core, _inspections, _interior, _physical_result = case
        node_id = q4["core_inlet_axial6_bindings"]["moments"]["M"][0] if "core_inlet_axial6_bindings" in q4 else 5
        lower = q4["exact_graph_nodes"][node_id]["lower"]
        q4["exact_graph_nodes"][lower]["numerator"] = 7

    def family(case: tuple[Any, ...]) -> None:
        q4, _core, _inspections, _interior, _physical_result = case
        q4["source_family"]["datum_enclosure_sha256"] = "d" * 64

    def missing(case: tuple[Any, ...]) -> None:
        _q4, _core, inspections, _interior, _physical_result = case
        del inspections[-1]

    def opaque(case: tuple[Any, ...]) -> None:
        q4, core, _inspections, _interior, _physical_result = case
        inlet_id = 6
        q4["exact_graph_nodes"][inlet_id]["operation"] = "opaque_q4_scope"
        core["exact_graph_nodes"][inlet_id]["operation"] = "opaque_q4_scope"

    def missing_saved_row(case: tuple[Any, ...]) -> None:
        _q4, _core, _inspections, interior, _physical_result = case
        del interior["current_core_interior_views"]["inlet"]["ordinary_rho4_axial6_moment_grids"]["C"]["rho0_Z2"]

    def physical_family(case: tuple[Any, ...]) -> None:
        _q4, _core, _inspections, _interior, physical = case
        physical["family_evidence"]["saved_core_family"]["actual_five_defect_family_sha256"] = "e" * 64

    def physical_finite(case: tuple[Any, ...]) -> None:
        _q4, _core, _inspections, _interior, physical = case
        physical["views"]["inlet"]["input_moments"]["A"]["z0"]["finite"][0][1] += 1

    def physical_tail(case: tuple[Any, ...]) -> None:
        _q4, _core, _inspections, _interior, physical = case
        physical["views"]["inlet"]["input_moments"]["C"]["z1"]["error_upper"][1] += 1

    reject("wrong_order", wrong_order)
    reject("wrong_row", wrong_row)
    reject("wrong_endpoint", wrong_endpoint)
    reject("family", family)
    reject("missing", missing)
    reject("opaque_q4_scope", opaque)
    reject("missing_saved_row", missing_saved_row)
    reject("physical_family", physical_family)
    reject("physical_finite", physical_finite)
    reject("physical_tail", physical_tail)
    checks["reject_missing_physical_result"] = _unresolved(_result(bridge, positive[:4] + (None,)))

    return {"passed": bool(checks) and all(checks.values()), "checks": checks}


__all__ = ["SCHEMA", "run_fixture"]
