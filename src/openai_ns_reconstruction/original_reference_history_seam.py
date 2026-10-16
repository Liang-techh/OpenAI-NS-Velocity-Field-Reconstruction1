"""Archive-only reference/restoration history seam fingerprints.

This module compares the accepted graph rows at reference native ``1`` and
restoration native ``0``.  It records exact source fingerprints and never
turns a finite-N error bound into a function identity.  In particular, the
native Jacobians and the physical ``U_y``/``Q_y`` rows remain outside this
certificate.
"""

from __future__ import annotations

from typing import Any, Mapping

from .original_active_viscous_m15 import ActiveEndpointFingerprint


_CHARTS = ("reference", "restoration")
_HISTORY_COMPONENTS = ("m", "h", "k", "e", "p")
_FINITE_FIELDS = ("E_N", "V_N", "delta_E", "delta_V", "A_over_N")
_OUTER_FIELDS = ("V", "m")


def _map(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    value = _map(value)
    nested = value.get(key)
    return _map(nested) if isinstance(nested, Mapping) else value


def _handle(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    return value if isinstance(value, int) else None


def _row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _handle(value)
    if index is None:
        return None
    try:
        result = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return result if isinstance(result, Mapping) else None


def _children(row: Mapping[str, Any] | None) -> tuple[int, ...]:
    if not isinstance(row, Mapping):
        return ()
    operation = row.get("operation")
    if operation in {"sum", "product", "logical_or", "logical_and"}:
        keys = ("arguments",)
    elif operation in {
        "negative", "analytic_unary", "exact_original_source_flat_sigma",
        "fractional_part", "regular_even_radial_axis_continuation",
    }:
        keys = ("argument",)
    elif operation == "positive_quotient":
        keys = ("numerator", "denominator")
    elif operation == "original_lazy_flat_branch":
        keys = ("active_body", "Delta", "eta", "flat_value", "flat_predicate",
                "original_collar_coordinate", "original_collar_upper")
    elif operation == "original_monotone_phase_inverse":
        keys = ("phase_function", "target_phase", "angle_lower", "angle_upper", "Delta", "eta")
    elif operation in {
        "substitute_original_inverse_angle", "substitute_original_native_coordinate",
        "substitute_free_query_at_native_coordinate",
    }:
        keys = ("body", "inverse_angle", "coordinate", "value", "lower", "upper")
    elif operation == "function_substitution":
        # The bound variable is lexical scope, not a value dependency.
        keys = ("expression", "value")
    elif operation == "definite_integral":
        keys = ("integrand", "lower", "upper")
    elif operation == "exact_real_comparison":
        keys = ("left", "right")
    else:
        return ()
    result: list[int] = []
    for key in keys:
        value = row.get(key)
        if isinstance(value, (list, tuple)):
            result.extend(index for item in value if (index := _handle(item)) is not None)
        elif (index := _handle(value)) is not None:
            result.append(index)
    return tuple(result)


def _reachable(nodes: Any, root: Any) -> set[int]:
    start = _handle(root)
    if start is None:
        return set()
    result: set[int] = set()
    pending = [start]
    while pending:
        index = pending.pop()
        if index in result:
            continue
        result.add(index)
        pending.extend(child for child in _children(_row(nodes, index)) if child not in result)
    return result


def _exact_zero(nodes: Any, value: Any) -> bool:
    row = _row(nodes, value)
    return bool(row and row.get("operation") == "exact_rational" and
                row.get("numerator") == 0 and row.get("denominator") == 1)


def _exact_one(nodes: Any, value: Any) -> bool:
    row = _row(nodes, value)
    return bool(row and row.get("operation") == "exact_rational" and
                row.get("numerator") == 1 and row.get("denominator") == 1)


def _endpoint_scope(nodes: Any, charts: Mapping[str, Any], failures: list[dict[str, Any]]):
    reference = _map(charts.get("reference"))
    restoration = _map(charts.get("restoration"))
    reference_coordinate = _handle(reference.get("coordinate"))
    restoration_coordinate = _handle(restoration.get("coordinate"))
    reference_domain = reference.get("domain")
    restoration_domain = restoration.get("domain")
    if (reference_coordinate is None or restoration_coordinate is None or
            not isinstance(reference_domain, (list, tuple)) or len(reference_domain) < 2 or
            not isinstance(restoration_domain, (list, tuple)) or len(restoration_domain) < 2):
        failures.append({"kind": "endpoint_scope", "reason": "missing chart coordinate/domain"})
        return None, None
    fingerprint = ActiveEndpointFingerprint(
        nodes, {"coordinate": reference_coordinate, "loop": charts}
    )
    try:
        scope = {
            reference_coordinate: fingerprint.at(reference_domain[1], {}),
            restoration_coordinate: fingerprint.at(restoration_domain[0], {}),
        }
    except Exception as error:  # archive rows must remain an explicit unproved item
        failures.append({"kind": "endpoint_scope", "reason": type(error).__name__ + ": " + str(error)})
        return fingerprint, None
    return fingerprint, scope


def _fingerprint(fingerprint: ActiveEndpointFingerprint, node: Any, scope: Mapping[int, Any]):
    handle = _handle(node)
    if handle is None:
        return None, {"node": node, "reason": "noninteger source handle"}
    try:
        value = fingerprint.at(handle, scope)
        return value[0], None
    except Exception as error:
        return None, {"node": handle, "reason": type(error).__name__ + ": " + str(error)}


def _compare(
    fingerprint: ActiveEndpointFingerprint,
    scope: Mapping[int, Any],
    left: Any,
    right: Any,
    group: str,
    quantity: str,
    order: int,
    identity_rows: list[dict[str, Any]],
    unproved: list[dict[str, Any]],
) -> bool:
    left_fp, left_error = _fingerprint(fingerprint, left, scope)
    right_fp, right_error = _fingerprint(fingerprint, right, scope)
    if left_error or right_error:
        unproved.append({
            "group": group, "quantity": quantity, "ordinary_Z_order": order,
            "left_node": _handle(left), "right_node": _handle(right),
            "left_error": left_error, "right_error": right_error,
        })
        return False
    if left_fp != right_fp:
        unproved.append({
            "group": group, "quantity": quantity, "ordinary_Z_order": order,
            "left_node": _handle(left), "right_node": _handle(right),
            "left_endpoint_fingerprint": left_fp, "right_endpoint_fingerprint": right_fp,
            "reason": "actual source fingerprints differ; no alias installed",
        })
        return False
    identity_rows.append({
        "group": group, "quantity": quantity, "ordinary_Z_order": order,
        "left_node": _handle(left), "right_node": _handle(right),
        "canonical_endpoint_fingerprint": left_fp,
    })
    return True


def _compare_lists(
    fingerprint: ActiveEndpointFingerprint,
    scope: Mapping[int, Any],
    left: Any,
    right: Any,
    group: str,
    quantity: str,
    identity_rows: list[dict[str, Any]],
    unproved: list[dict[str, Any]],
    *,
    minimum_orders: int = 3,
) -> int:
    if not isinstance(left, (list, tuple)) or not isinstance(right, (list, tuple)):
        unproved.append({"group": group, "quantity": quantity, "reason": "missing ordered source list"})
        return 0
    if len(left) != len(right) or len(left) < minimum_orders:
        unproved.append({
            "group": group, "quantity": quantity,
            "left_length": len(left), "right_length": len(right),
            "reason": "source derivative order list is incomplete or mismatched",
        })
    count = 0
    for order, (left_node, right_node) in enumerate(zip(left, right)):
        if _compare(fingerprint, scope, left_node, right_node, group, quantity, order,
                    identity_rows, unproved):
            count += 1
    return count


def _z_free(nodes: Any, root: Any, coordinate: int) -> dict[str, Any]:
    """Check source Z-independence of geometry handles structurally."""
    visited: set[int] = set()
    unknown: list[dict[str, Any]] = []

    def visit(value: Any) -> bool:
        index = _handle(value)
        if index is None:
            unknown.append({"node": value, "reason": "noninteger geometry handle"})
            return False
        if index in visited:
            return True
        visited.add(index)
        row = _row(nodes, index)
        if row is None:
            unknown.append({"node": index, "reason": "missing geometry row"})
            return False
        operation = row.get("operation")
        if operation == "bound_variable":
            return index == coordinate
        if operation in {"exact_rational", "mathematical_pi"}:
            return True
        if operation in {"exact_positive_integer_expression", "current_original_source_parameter",
                         "original_source_parameter"}:
            return row.get("Z_independent") is True
        if operation == "function_substitution":
            return visit(row.get("value")) and visit(row.get("expression"))
        children = _children(row)
        if not children:
            unknown.append({"node": index, "operation": operation, "reason": "opaque geometry row"})
            return False
        return all(visit(child) for child in children)

    passed = visit(root)
    return {"passed": bool(passed), "root": _handle(root),
            "visited_handles": sorted(visited), "unknown": unknown[:8],
            "unknown_count": len(unknown)}


def _incoming_scope_rows(
    nodes: Any,
    fingerprint: ActiveEndpointFingerprint,
    scope: Mapping[int, Any],
    chart_name: str,
    chart: Mapping[str, Any],
    unproved: list[dict[str, Any]],
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    incoming = _map(chart.get("incoming"))
    for component in _HISTORY_COMPONENTS:
        rows = incoming.get(component)
        checks = []
        if not isinstance(rows, (list, tuple)) or len(rows) < 3:
            unproved.append({"group": "incoming_scope", "chart": chart_name,
                             "quantity": component, "reason": "missing incoming history orders"})
            result[component] = {"passed": False, "rows": checks}
            continue
        for order, handle in enumerate(rows):
            source = _row(nodes, handle)
            local = bool(
                source and source.get("operation") == "function_substitution" and
                source.get("Z_independent_substitution") is True and
                _exact_one(nodes, source.get("value")) and
                _handle(source.get("variable")) is not None and
                _handle(source.get("expression")) is not None
            )
            if not local:
                unproved.append({"group": "incoming_scope", "chart": chart_name,
                                 "quantity": component, "ordinary_Z_order": order,
                                 "node": _handle(handle),
                                 "reason": "incoming history is not an exact fixed endpoint substitution"})
            # Evaluate the actual source row so empty integral preservation is
            # measured from the accepted graph, not from metadata.
            _, error = _fingerprint(fingerprint, handle, scope)
            if error:
                local = False
                unproved.append({"group": "incoming_scope", "chart": chart_name,
                                 "quantity": component, "ordinary_Z_order": order,
                                 "node": _handle(handle), "reason": error["reason"]})
            checks.append({
                "order": order, "node": _handle(handle),
                "expression": _handle(source.get("expression")) if source else None,
                "variable": _handle(source.get("variable")) if source else None,
                "passed": local,
            })
        result[component] = {"passed": all(row["passed"] for row in checks), "rows": checks}
    return result


def _empty_integrals(nodes: Any, roots: Mapping[str, Any], empty: set[int]) -> dict[str, list[int]]:
    result: dict[str, list[int]] = {}
    for name, root in roots.items():
        values = root if isinstance(root, (list, tuple)) else (root,)
        reachable = set().union(*(_reachable(nodes, value) for value in values))
        result[name] = sorted(reachable & empty)
    return result


def _chart_value(bindings: Mapping[str, Any], chart_name: str) -> Mapping[str, Any]:
    return _map(_map(bindings.get("charts")).get(chart_name))


def build_reference_history_seam(nodes: Any, l: Mapping[str, Any], modified: Mapping[str, Any], outer: Mapping[str, Any]) -> dict[str, Any]:
    """Build the accepted reference native-1/restoration native-0 seam packet."""
    loop = _unwrap(l, "prescribed_loop_C5_bindings")
    modified_bindings = _unwrap(modified, "modified_Q4_pressure_C3_bindings")
    outer_bindings = _unwrap(outer, "prescribed_R100_Rm_C6_bindings")
    charts = _map(loop.get("charts"))
    failures: list[dict[str, Any]] = []
    missing = [name for name in _CHARTS if not isinstance(charts.get(name), Mapping)]
    if missing:
        return {
            "passed": False,
            "method": "archived_active_endpoint_fingerprint_reference_history_seam",
            "unproved": [{"kind": "charts", "reason": "missing chart(s): " + ", ".join(missing)}],
            "claims_not_made": ["native_Jacobian_equality", "U_y_equality", "Q_y_equality", "D_native_B_zero"],
        }
    fingerprint, scope = _endpoint_scope(nodes, charts, failures)
    if fingerprint is None or scope is None:
        return {
            "passed": False,
            "method": "archived_active_endpoint_fingerprint_reference_history_seam",
            "unproved": failures,
            "claims_not_made": ["native_Jacobian_equality", "U_y_equality", "Q_y_equality", "D_native_B_zero"],
        }

    reference = _chart_value(loop, "reference")
    restoration = _chart_value(loop, "restoration")
    identity_rows: list[dict[str, Any]] = []
    unproved: list[dict[str, Any]] = list(failures)

    # Shared radius and global phase are compared at the actual endpoints;
    # no phase value is substituted or reset to 0/1.
    for quantity in ("radius", "global_offset", "actual_phase"):
        _compare(fingerprint, scope, reference.get(quantity), restoration.get(quantity),
                 "loop_geometry", quantity, 0, identity_rows, unproved)

    modified_charts = _map(modified_bindings.get("charts"))
    outer_charts = _map(outer_bindings.get("charts"))
    modified_identity_counts: dict[str, int] = {}
    modified_order_availability: dict[str, dict[str, int]] = {}
    if not modified_charts or not outer_charts:
        unproved.append({"kind": "source_archives", "reason": "modified or outer chart bindings missing"})

    # Physical U, m and m_y are directly comparable.  U_y/Q_y are deliberately
    # left open because their native derivative theorem is owned by the parent.
    for quantity in ("U", "m", "m_y"):
        left = _map(_map(modified_charts.get("reference")).get("physical")).get(quantity)
        right = _map(_map(modified_charts.get("restoration")).get("physical")).get(quantity)
        modified_identity_counts[quantity] = _compare_lists(
            fingerprint, scope, left, right, "modified_physical", quantity,
            identity_rows, unproved,
        )
    for quantity in ("U", "m", "m_y", "Q_y"):
        left = _map(_map(modified_charts.get("reference")).get("physical")).get(quantity)
        right = _map(_map(modified_charts.get("restoration")).get("physical")).get(quantity)
        left_length = len(left) if isinstance(left, (list, tuple)) else 0
        right_length = len(right) if isinstance(right, (list, tuple)) else 0
        modified_order_availability[quantity] = {
            "reference_length": left_length, "restoration_length": right_length,
        }
        if min(left_length, right_length) < 3:
            unproved.append({
                "group": "modified_physical_order_availability", "quantity": quantity,
                "reason": "Qy_Z input derivative order through ordinary Z2 is unavailable",
                "reference_length": left_length, "restoration_length": right_length,
            })

    # The loop archive is the source of corrected finite-N V/m and all five
    # retained histories.  Keep every stored order while explicitly checking
    # the first three ordinary-Z orders needed by Q_y_Z.
    finite_identity_counts: dict[str, int] = {}
    for quantity in _FINITE_FIELDS:
        finite_identity_counts[quantity] = _compare_lists(
            fingerprint, scope,
            _map(reference.get("finite_N")).get(quantity),
            _map(restoration.get("finite_N")).get(quantity),
            "finite_N", quantity, identity_rows, unproved,
        )
    for component in _HISTORY_COMPONENTS:
        finite_identity_counts["finite_N." + component] = _compare_lists(
            fingerprint, scope,
            _map(_map(reference.get("finite_N")).get("five_signed_density_increments")).get(component),
            _map(_map(restoration.get("finite_N")).get("five_signed_density_increments")).get(component),
            "finite_N_signed_density", component, identity_rows, unproved,
        )
        finite_identity_counts["history." + component] = _compare_lists(
            fingerprint, scope, _map(reference.get("history")).get(component),
            _map(restoration.get("history")).get(component), "retained_history", component,
            identity_rows, unproved,
        )

    incoming_scope = {
        name: _incoming_scope_rows(nodes, fingerprint, scope, name, _chart_value(loop, name), unproved)
        for name in _CHARTS
    }

    # Verify that each outer chart uses the same radius, physical measure and
    # native phase coordinate as its loop chart.  The distinct outer leading
    # V/m sources are compared across the two charts below.
    outer_identity_counts: dict[str, int] = {}
    for chart_name in _CHARTS:
        loop_chart = _chart_value(loop, chart_name)
        outer_chart = _chart_value(outer_bindings, chart_name)
        if not outer_chart:
            unproved.append({"group": "outer_source", "chart": chart_name, "reason": "missing outer chart"})
            continue
        # Compare outer geometry handles to the corresponding loop handles.
        _compare(fingerprint, scope, outer_chart.get("radius"), loop_chart.get("radius"),
                 "outer_geometry", chart_name + ".radius", 0, identity_rows, unproved)
        _compare(fingerprint, scope, outer_chart.get("physical_measure"), loop_chart.get("physical_measure"),
                 "outer_geometry", chart_name + ".physical_measure", 0, identity_rows, unproved)
        _compare(fingerprint, scope, outer_chart.get("phase"), loop_chart.get("coordinate"),
                 "outer_geometry", chart_name + ".phase", 0, identity_rows, unproved)

    # Outer reference/restoration leading fields are a useful source-family
    # check; compare V/m and V_y through ordinary Z2 only.
    outer_ref = _chart_value(outer_bindings, "reference")
    outer_rest = _chart_value(outer_bindings, "restoration")
    for quantity in ("V", "m"):
        outer_identity_counts["cross." + quantity] = _compare_lists(
            fingerprint, scope, _map(outer_ref.get("fields")).get(quantity),
            _map(outer_rest.get("fields")).get(quantity), "outer_cross_chart", quantity,
            identity_rows, unproved,
        )
    outer_identity_counts["cross.V_y"] = _compare_lists(
        fingerprint, scope, _map(outer_ref.get("radial")).get("V_y"),
        _map(outer_rest.get("radial")).get("V_y"), "outer_cross_chart", "V_y",
        identity_rows, unproved,
    )

    geometry_Z_free = {}
    for name in _CHARTS:
        chart = _chart_value(loop, name)
        coordinate = _handle(chart.get("coordinate"))
        geometry_Z_free[name] = {
            quantity: _z_free(nodes, chart.get(quantity), coordinate)
            for quantity in ("radius", "global_offset", "actual_phase")
        }
        for quantity, check in geometry_Z_free[name].items():
            if not check["passed"]:
                unproved.append({"group": "geometry_Z_free", "chart": name,
                                 "quantity": quantity, "reason": "geometry graph is not structurally Z-free",
                                 "details": check})

    incoming_roots = {
        name + ".incoming." + component: _chart_value(loop, name).get("incoming", {}).get(component)
        for name in _CHARTS for component in _HISTORY_COMPONENTS
    }
    empty_integrals = _empty_integrals(nodes, incoming_roots, fingerprint.empty)

    required_orders = {
        "ordinary_Z_orders_checked": [0, 1, 2],
        "finite_N_V_N": all(finite_identity_counts.get("V_N", 0) >= 3 for _ in [0]),
        "finite_N_m_increment": all(finite_identity_counts.get("finite_N.m", 0) >= 3 for _ in [0]),
        "retained_history": all(finite_identity_counts.get("history." + name, 0) >= 3 for name in _HISTORY_COMPONENTS),
        "modified_m_y": modified_identity_counts.get("m_y", 0) >= 3,
        "modified_Q_y_Z2_inputs_available": all(
            modified_order_availability.get(name, {}).get("reference_length", 0) >= 3 and
            modified_order_availability.get(name, {}).get("restoration_length", 0) >= 3
            for name in ("U", "m", "m_y", "Q_y")
        ),
    }
    if not all(required_orders.values()):
        unproved.append({"group": "derivative_orders", "reason": "ordinary Z2 source rows incomplete",
                         "checks": required_orders})

    return {
        "passed": not unproved,
        "method": "archived_active_endpoint_fingerprint_reference_history_seam",
        "source_only": True,
        "endpoint_scope": {
            "reference_chart": "reference",
            "restoration_chart": "restoration",
            "reference_native_endpoint": reference.get("domain", [None, None])[1],
            "restoration_native_endpoint": restoration.get("domain", [None, None])[0],
            "coordinate_handles": [reference.get("coordinate"), restoration.get("coordinate")],
            "global_phase_not_reset": True,
        },
        "identity_rows": identity_rows,
        "identity_row_count": len(identity_rows),
        "modified_physical_identity_counts": modified_identity_counts,
        "modified_physical_order_availability": modified_order_availability,
        "finite_N_and_history_identity_counts": finite_identity_counts,
        "outer_source_identity_counts": outer_identity_counts,
        "incoming_history_scope": incoming_scope,
        "empty_integrals_preserving_incoming_histories": empty_integrals,
        "geometry_Z_free": geometry_Z_free,
        "required_orders": required_orders,
        "unproved": unproved,
        "claims_not_made": [
            "native_Jacobian_equality",
            "modified_U_y_equality",
            "modified_Q_y_equality",
            "D_native_B_zero",
            "numeric_source_field_values",
            "finite_N_error_cap_as_function",
        ],
    }


__all__ = ["build_reference_history_seam"]
