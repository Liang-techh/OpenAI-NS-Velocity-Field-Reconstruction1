"""Source-bound symbolic main/gap signed seam.

The runner supplies two already authenticated, decoded archive mappings.  This
leaf intentionally does not open an archive, execute a source function,
perform a broad graph traversal, or turn a graph handle into a number. It reads only the selected
ordinary ``Z0``/``Z1`` rows, chart query/radius handles, the three geometry
handles, the two selected ``own`` domain-edge rows, and the two selected query
rows in each graph.

The returned expression objects are JSON-ready symbolic records.  In
particular, ``source_at(expression_ref, query_ref, exact11)`` is a closure,
not a call.  The two endpoint sides remain independent all the way through
``K`` and ``K_Z``; no main/gap cancellation or source value is manufactured.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from copy import deepcopy
from typing import Any


SCHEMA = "original_f_main_gap_signed_seam_v1"
COMPONENT = SCHEMA
ORDERS = (0, 1)
CHARTS = ("main", "gap")
FIELDS = ("U", "Q", "Q_y")
F_SOURCE_KEY = "Rp_native_C5_pulse_source_bindings"
OWN_SOURCE_KEY = "own_main_gap_end_fields"
F_SOURCE_POINTER = "$.Rp_native_C5_pulse_source_bindings"
OWN_SOURCE_POINTER = "$.own_main_gap_end_fields"
FAMILY_KEYS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)

# This is the one explicitly protected historical node required by the
# source contract.  It is compared as an opaque handle only: it is never
# indexed, even when a synthetic graph happens to be long enough to contain
# it.
PROTECTED_NODE_IDS = frozenset({33982})


class FMainGapSignedSeamInputError(ValueError):
    """Raised only when the two public build arguments are not mappings."""


def _is_sequence(value: Any) -> bool:
    return isinstance(value, Sequence) and not isinstance(
        value, (str, bytes, bytearray)
    )


def _copy_json(value: Any) -> Any:
    """Copy small metadata values without copying a graph or source object."""

    if isinstance(value, Mapping):
        return {str(key): _copy_json(item) for key, item in value.items()}
    if _is_sequence(value):
        return [_copy_json(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return repr(value)


def _expr(operator: str, *arguments: Any, **metadata: Any) -> dict[str, Any]:
    """Create one unsimplified symbolic expression node."""

    result: dict[str, Any] = {
        "operator": operator,
        "arguments": list(arguments),
        "source_values_selected": False,
    }
    result.update(metadata)
    return result


def _literal(text: str) -> dict[str, Any]:
    # Literals in this module are syntax, not recovered graph values.
    return {
        "operator": "exact_literal",
        "text": text,
        "source_values_selected": False,
        "node_id_is_not_a_numeric_value": True,
    }


def _formal(name: str, **metadata: Any) -> dict[str, Any]:
    result = _expr("formal_source_symbol", _literal(name), name=name)
    result.update(metadata)
    return result


def _handle(value: Any) -> int | None:
    """Extract an opaque integer handle, never treating it as a scalar."""

    if type(value) is int:
        return value
    if isinstance(value, Mapping):
        for key in ("source_node_id", "node", "handle", "node_id", "root"):
            candidate = value.get(key)
            if type(candidate) is int:
                return candidate
    return None


def _opaque_ref(value: Any, role: str, owner: str, namespace: str) -> dict[str, Any]:
    node = _handle(value)
    if node is None:
        return {
            "status": "unresolved_source_ref",
            "kind": "opaque_source_ref",
            "owner": owner,
            "namespace": namespace,
            "source_role": role,
            "source_node_id": None,
            "node_id_is_not_a_numeric_value": True,
            "source_values_selected": False,
        }
    return {
        "status": "opaque_source_ref",
        "kind": "opaque_source_ref",
        "owner": owner,
        "namespace": namespace,
        "source_role": role,
        "source_node_id": node,
        "node_id_is_not_a_numeric_value": True,
        "source_values_selected": False,
    }


def _unresolved_ref(role: str, reason: str) -> dict[str, Any]:
    return {
        "status": "unresolved_source_ref",
        "kind": "opaque_source_ref",
        "source_role": role,
        "reason": reason,
        "source_node_id": None,
        "node_id_is_not_a_numeric_value": True,
        "source_values_selected": False,
    }


def _require_mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise FMainGapSignedSeamInputError(f"{label} must be a mapping")
    return value


def _family(archive: Mapping[str, Any], source: Mapping[str, Any]) -> dict[str, Any]:
    raw = archive.get("source_family")
    if not isinstance(raw, Mapping):
        raw = source.get("source_family")
    if not isinstance(raw, Mapping):
        return {
            "present": False,
            "complete": False,
            "fields": {key: None for key in FAMILY_KEYS},
        }
    fields = {key: _copy_json(raw.get(key)) for key in FAMILY_KEYS}
    return {
        "present": True,
        "complete": all(key in raw and raw.get(key) is not None for key in FAMILY_KEYS),
        "fields": fields,
    }


def _family_check(
    f_family: Mapping[str, Any], own_family: Mapping[str, Any]
) -> dict[str, Any]:
    left = f_family.get("fields", {})
    right = own_family.get("fields", {})
    complete = f_family.get("complete") is True and own_family.get("complete") is True
    equal = complete and all(left.get(key) == right.get(key) for key in FAMILY_KEYS)
    return {
        "passed": bool(equal),
        "fields": list(FAMILY_KEYS),
        "F": _copy_json(left),
        "own": _copy_json(right),
        "all_three_fields_present": bool(complete),
        "source_family_function_identity_proved": False,
    }


def _chart(source: Mapping[str, Any], chart: str) -> Mapping[str, Any]:
    value = source.get(chart)
    return value if isinstance(value, Mapping) else {}


def _leading(source: Mapping[str, Any], chart: str) -> Mapping[str, Any]:
    leading = source.get("leading_fields")
    if not isinstance(leading, Mapping):
        return {}
    value = leading.get(chart)
    return value if isinstance(value, Mapping) else {}


def _coordinate(source: Mapping[str, Any], chart: str) -> Mapping[str, Any]:
    coordinates = source.get("pulse_coordinates")
    if not isinstance(coordinates, Mapping):
        return {}
    value = coordinates.get(chart)
    return value if isinstance(value, Mapping) else {}


def _row_pair(value: Any) -> tuple[Any, Any, str]:
    """Read exactly the ordinary Z0/Z1 pair from a saved field record."""

    if _is_sequence(value) and len(value) >= 2:
        return value[0], value[1], "sequence[0:2]"
    if not isinstance(value, Mapping):
        return None, None, "missing"
    for key in (
        "ordinary_Z0_to_Z1",
        "ordinary_Z0_to_Z2",
        "ordinary_Z_rows",
        "rows",
        "all_saved_rows",
    ):
        nested = value.get(key)
        if _is_sequence(nested) and len(nested) >= 2:
            return nested[0], nested[1], key + "[0:2]"
        if isinstance(nested, Mapping):
            first = next(
                (nested.get(key0) for key0 in (0, "0", "Z0", "ordinary_Z0") if key0 in nested),
                None,
            )
            second = next(
                (nested.get(key1) for key1 in (1, "1", "Z1", "ordinary_Z1") if key1 in nested),
                None,
            )
            if first is not None and second is not None:
                return first, second, key + "{0,1}"
    first = next(
        (
            value.get(key)
            for key in ("ordinary_Z0", "Z0", "z0", "0", "value")
            if key in value
        ),
        None,
    )
    second = next(
        (
            value.get(key)
            for key in ("ordinary_Z1", "Z1", "z1", "1", "derivative")
            if key in value
        ),
        None,
    )
    if first is not None and second is not None:
        return first, second, "named_rows"
    return None, None, "missing"


def _selected_rows(source: Mapping[str, Any], chart: str) -> dict[str, Any]:
    leading = _leading(source, chart)
    result: dict[str, Any] = {}
    for field in FIELDS:
        first, second, selected_from = _row_pair(leading.get(field))
        result[field] = {
            "Z0": _opaque_ref(
                first,
                f"{F_SOURCE_POINTER}.leading_fields.{chart}.{field}[0]",
                "F",
                "F_exact_graph",
            ),
            "Z1": _opaque_ref(
                second,
                f"{F_SOURCE_POINTER}.leading_fields.{chart}.{field}[1]",
                "F",
                "F_exact_graph",
            ),
            "selected_from": selected_from,
            "ordinary_Z1_from_F_row_1": True,
            "row_reorder_or_owner_identity_authenticated": False,
            "owner_identity_proved": False,
            "source_values_selected": False,
        }
    return result


def _coord_refs(source: Mapping[str, Any], chart: str) -> dict[str, Any]:
    coordinate = _coordinate(source, chart)
    return {
        "query": _opaque_ref(
            coordinate.get("query"),
            f"{F_SOURCE_POINTER}.pulse_coordinates.{chart}.query",
            "F",
            "F_exact_graph",
        ),
        "radius": _opaque_ref(
            coordinate.get("radius"),
            f"{F_SOURCE_POINTER}.pulse_coordinates.{chart}.radius",
            "F",
            "F_exact_graph",
        ),
    }


def _geometry_refs(source: Mapping[str, Any]) -> dict[str, Any]:
    geometry = source.get("geometry")
    geometry = geometry if isinstance(geometry, Mapping) else {}
    return {
        name: _opaque_ref(
            geometry.get(name),
            f"{F_SOURCE_POINTER}.geometry.{name}",
            "F",
            "F_exact_graph",
        )
        for name in ("delta", "mu", "Rp")
    }


def _graph_nodes(archive: Mapping[str, Any]) -> Any:
    value = archive.get("exact_graph_nodes")
    return value


def _node_at(nodes: Any, handle: Any, label: str) -> tuple[Mapping[str, Any] | None, str]:
    """Index only one requested row; never iterate a graph."""

    node = _handle(handle)
    if node is None:
        return None, "invalid_handle"
    if node in PROTECTED_NODE_IDS:
        return None, "protected_node_not_dereferenced"
    if not isinstance(nodes, (Mapping, Sequence)) or isinstance(
        nodes, (str, bytes, bytearray)
    ):
        return None, "graph_missing"
    try:
        if isinstance(nodes, Mapping):
            if node in nodes:
                row = nodes[node]
            elif str(node) in nodes:
                row = nodes[str(node)]
            else:
                return None, "row_missing"
        else:
            if node < 0 or node >= len(nodes):
                return None, "row_missing"
            row = nodes[node]
    except (KeyError, IndexError, TypeError):
        return None, "row_missing"
    if not isinstance(row, Mapping):
        return None, "row_not_mapping"
    return row, "row_read"


def _operation(row: Mapping[str, Any] | None) -> Any:
    if not isinstance(row, Mapping):
        return None
    operation = row.get("operation", row.get("op"))
    return operation if isinstance(operation, str) else None


def _bound_signature(row: Mapping[str, Any] | None) -> tuple[Any, Any] | None:
    operation = _operation(row)
    if operation != "bound_variable" or not isinstance(row, Mapping):
        return None
    name = row.get("name", row.get("variable"))
    if not isinstance(name, str) or not name:
        return None
    return (operation, name)


def _edge_domain(chart_record: Mapping[str, Any], side: str) -> Any:
    domain = chart_record.get("domain")
    if isinstance(domain, Mapping):
        return domain.get(side)
    if _is_sequence(domain) and len(domain) >= 2:
        return domain[0 if side == "lower" else 1]
    key = "domain_" + side
    if key in chart_record:
        return chart_record.get(key)
    return chart_record.get(side)


def _own_chart(source: Mapping[str, Any], chart: str) -> Mapping[str, Any]:
    value = source.get(chart)
    return value if isinstance(value, Mapping) else {}


def _exact11(
    own_nodes: Any, edge_handle: Any, role: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Read one own edge row and admit only typed exact rational 11."""

    ref = _opaque_ref(edge_handle, role, "own", "own_exact_graph")
    node = _handle(edge_handle)
    if node is None:
        return (
            {
                "status": "unresolved_exact11",
                "kind": "exact_rational_endpoint",
                "source_ref": ref,
                "source_values_selected": False,
                "typed_exact_rational_11": False,
                "row_read": False,
            },
            {"passed": False, "reason": "edge handle is absent", "row_read": False},
        )
    row, status = _node_at(own_nodes, node, role)
    if status != "row_read":
        return (
            {
                "status": "unresolved_exact11",
                "kind": "exact_rational_endpoint",
                "source_ref": ref,
                "source_values_selected": False,
                "typed_exact_rational_11": False,
                "row_read": False,
                "protected_or_unread_reason": status,
            },
            {
                "passed": False,
                "reason": status,
                "row_read": False,
                "protected_node_not_dereferenced": status == "protected_node_not_dereferenced",
            },
        )
    operation = _operation(row)
    numerator = row.get("numerator", row.get("value"))
    denominator = row.get("denominator", 1)
    passed = (
        operation == "exact_rational"
        and type(numerator) is int
        and type(denominator) is int
        and numerator == 11
        and denominator == 1
    )
    endpoint = {
        "status": "typed_exact11" if passed else "unresolved_exact11",
        "kind": "exact_rational_endpoint",
        "source_ref": ref,
        "numerator": "11" if passed else None,
        "denominator": "1" if passed else None,
        "source_values_selected": False,
        "typed_exact_rational_11": bool(passed),
        "row_read": True,
        "row_operation": operation,
        "graph_handle_is_opaque": True,
    }
    return endpoint, {
        "passed": bool(passed),
        "reason": "typed exact_rational 11" if passed else "edge row is not exact_rational 11/1",
        "row_read": True,
        "protected_node_not_dereferenced": False,
    }


def _domain_edge_in_f(source: Mapping[str, Any], chart: str, side: str) -> Any:
    """Find an optional F edge pointer without searching the graph."""

    coordinates = _coordinate(source, chart)
    for container in (coordinates, _chart(source, chart)):
        if not isinstance(container, Mapping):
            continue
        domain = container.get("domain")
        if isinstance(domain, Mapping) and side in domain:
            return domain.get(side)
        if _is_sequence(domain) and len(domain) >= 2:
            return domain[0 if side == "lower" else 1]
        key = "domain_" + side
        if key in container:
            return container.get(key)
    return None


def _query_identity(
    f_archive: Mapping[str, Any], own_archive: Mapping[str, Any],
    f_source: Mapping[str, Any], own_source: Mapping[str, Any], chart: str,
) -> dict[str, Any]:
    f_chart = _coordinate(f_source, chart)
    own_chart = _own_chart(own_source, chart)
    f_query = _handle(f_chart.get("query"))
    own_query = _handle(own_chart.get("query"))
    same_handle = f_query is not None and own_query is not None and f_query == own_query
    f_row, f_status = _node_at(_graph_nodes(f_archive), f_query, f"F.{chart}.query")
    own_row, own_status = _node_at(_graph_nodes(own_archive), own_query, f"own.{chart}.query")
    f_signature = _bound_signature(f_row)
    own_signature = _bound_signature(own_row)
    rows_match = f_signature is not None and f_signature == own_signature
    return {
        "F_query": _opaque_ref(
            f_chart.get("query"),
            f"{F_SOURCE_POINTER}.pulse_coordinates.{chart}.query",
            "F",
            "F_exact_graph",
        ),
        "own_query": _opaque_ref(
            own_chart.get("query"),
            f"{OWN_SOURCE_POINTER}.{chart}.query",
            "own",
            "own_exact_graph",
        ),
        "same_integer_handle": bool(same_handle),
        "F_query_row_status": f_status,
        "own_query_row_status": own_status,
        "F_selected_row": {
            "operation": f_signature[0] if f_signature else None,
            "name": f_signature[1] if f_signature else None,
            "row_read": f_status == "row_read",
        },
        "own_selected_row": {
            "operation": own_signature[0] if own_signature else None,
            "name": own_signature[1] if own_signature else None,
            "row_read": own_status == "row_read",
        },
        "selected_bound_variable_rows_match": bool(rows_match),
        "query_reference_rows_match": bool(same_handle and rows_match),
        "semantic_query_identity_proved": False,
        "namespace_identity_conditional": True,
    }


def _source_at(
    expression_ref: Any, query_ref: Any, exact11: Any, role: str
) -> dict[str, Any]:
    return _expr(
        "source_at",
        expression_ref,
        query_ref,
        exact11,
        source_role=role,
        expression_ref=expression_ref,
        query_ref=query_ref,
        endpoint_ref=exact11,
        source_executed=False,
        capture_safe_endpoint_proved=False,
    )


def _partial_z(value: Any, role: str) -> dict[str, Any]:
    return _expr(
        "partial_Z_at_fixed_query",
        value,
        source_role=role,
        unresolved=True,
        fixed_query=True,
        source_expression=True,
        source_values_selected=False,
    )


def _source_atom(name: str, ref: Mapping[str, Any], role: str) -> dict[str, Any]:
    return _expr(
        "formal_source_atom",
        _literal(name),
        name=name,
        source_ref=ref,
        source_role=role,
        source_values_selected=False,
    )


def _pow(value: Any, exponent: str) -> dict[str, Any]:
    return _expr("power", value, _literal(exponent), formula_text=f"({exponent})-power")


def _prod(values: Sequence[Any]) -> dict[str, Any]:
    return _expr("product", *values)


def _sum(values: Sequence[Any]) -> dict[str, Any]:
    return _expr("sum", *values)


def _negative(value: Any) -> dict[str, Any]:
    return _expr("negative", value)


def _product_z(values: Sequence[Any], derivatives: Sequence[Any], formula_text: str) -> dict[str, Any]:
    terms: list[Any] = []
    for index in range(len(values)):
        factors = [derivatives[index] if index == j else values[j] for j in range(len(values))]
        terms.append(_prod(factors))
    return _expr("sum", *terms, formula_text=formula_text, product_rule_expanded=True)


def _radius_shape(
    f_archive: Mapping[str, Any], chart: str, coords: Mapping[str, Any], geometry: Mapping[str, Any]
) -> dict[str, Any]:
    """Verify only the direct current-radius rows needed for a chain rule."""

    radius_handle = _handle(coords.get("radius"))
    rp_handle = _handle(geometry.get("Rp", {}))
    mu_handle = _handle(geometry.get("mu", {}))
    query_handle = _handle(coords.get("query"))
    result: dict[str, Any] = {
        "chart": chart,
        "verified": False,
        "rows_read": [],
        "direct_shape": "R=Rp*exp(query/mu)",
        "source_values_selected": False,
    }
    if None in (radius_handle, rp_handle, mu_handle, query_handle):
        result["reason"] = "radius, Rp, mu, or query handle is missing"
        return result
    radius_row, radius_status = _node_at(
        _graph_nodes(f_archive), radius_handle, f"F.{chart}.radius"
    )
    result["rows_read"].append({"role": f"F.{chart}.radius", "handle": radius_handle, "status": radius_status})
    if radius_status != "row_read":
        result["reason"] = radius_status
        return result
    operation = _operation(radius_row)
    args = radius_row.get("arguments", radius_row.get("operands"))
    if operation not in {"product", "mul"} or not _is_sequence(args) or len(args) != 2:
        result["reason"] = "radius row is not a two-factor product"
        return result
    if rp_handle not in {_handle(args[0]), _handle(args[1])}:
        result["reason"] = "radius product has no geometry Rp factor"
        return result
    exp_handle = args[1] if _handle(args[0]) == rp_handle else args[0]
    exp_id = _handle(exp_handle)
    if exp_id is None:
        result["reason"] = "radius exponential handle is missing"
        return result
    exp_row, exp_status = _node_at(_graph_nodes(f_archive), exp_id, f"F.{chart}.radius.exp")
    result["rows_read"].append({"role": f"F.{chart}.radius.exp", "handle": exp_id, "status": exp_status})
    if exp_status != "row_read":
        result["reason"] = exp_status
        return result
    if _operation(exp_row) != "analytic_unary" or exp_row.get("name") != "exp":
        result["reason"] = "radius factor is not exp"
        return result
    y_id = _handle(exp_row.get("argument"))
    if y_id is None:
        result["reason"] = "exp argument is missing"
        return result
    y_row, y_status = _node_at(_graph_nodes(f_archive), y_id, f"F.{chart}.radius.exp.argument")
    result["rows_read"].append({"role": f"F.{chart}.radius.exp.argument", "handle": y_id, "status": y_status})
    if y_status != "row_read":
        result["reason"] = y_status
        return result
    y_operation = _operation(y_row)
    y_args = y_row.get("arguments", y_row.get("operands"))
    if y_operation == "positive_quotient":
        y_numerator = y_row.get("numerator")
        y_denominator = y_row.get("denominator")
        ratio_ok = (
            _handle(y_numerator) == query_handle
            and _handle(y_denominator) == mu_handle
        )
        query_match = _handle(y_numerator) == query_handle
        mu_match = _handle(y_denominator) == mu_handle
    else:
        ratio_ok = (
            y_operation in {"quotient", "div"}
            and _is_sequence(y_args)
            and len(y_args) == 2
            and _handle(y_args[0]) == query_handle
            and _handle(y_args[1]) == mu_handle
        )
        query_match = _handle(y_args[0]) == query_handle if _is_sequence(y_args) and len(y_args) >= 1 else False
        mu_match = _handle(y_args[1]) == mu_handle if _is_sequence(y_args) and len(y_args) >= 2 else False
    result["ratio_row"] = {
        "operation": y_operation,
        "query_handle_match": query_match,
        "mu_handle_match": mu_match,
    }
    result["verified"] = bool(ratio_ok)
    result["reason"] = "direct current radius shape verified" if ratio_ok else "query/mu quotient shape not verified"
    return result


def _radius_z(
    radius: Any, geometry: Mapping[str, Any], shape: Mapping[str, Any], role: str
) -> dict[str, Any]:
    unresolved = _partial_z(radius, role + ".partial_Z_at_fixed_query")
    if shape.get("verified") is not True:
        unresolved["shape_candidate"] = None
        unresolved["formula_text"] = "partial_Z_at_fixed_query(R)"
        return unresolved
    rp = _source_atom("Rp", geometry["Rp"], role + ".Rp")
    mu = _source_atom("mu", geometry["mu"], role + ".mu")
    rp_z = _partial_z(rp, role + ".Rp_Z")
    mu_z = _partial_z(mu, role + ".mu_Z")
    xi = _literal("11")
    expanded = _prod(
        [
            _expr("exp", _expr("quotient", xi, mu), formula_text="exp(11/mu)"),
            _expr(
                "subtract",
                rp_z,
                _expr(
                    "quotient",
                    _prod([xi, rp, mu_z]),
                    _pow(mu, "2"),
                    formula_text="11*Rp*mu_Z/mu^2",
                ),
                formula_text="Rp_Z-11*Rp*mu_Z/mu^2",
            ),
        ]
    )
    unresolved["shape_candidate"] = expanded
    unresolved["formula_text"] = "exp(11/mu)*(Rp_Z-11*Rp*mu_Z/mu^2)"
    unresolved["direct_radius_shape_used"] = True
    return unresolved


def _endpoint(
    chart: str,
    endpoint: str,
    f_archive: Mapping[str, Any],
    own_archive: Mapping[str, Any],
    f_source: Mapping[str, Any],
    own_source: Mapping[str, Any],
    geometry: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    own_chart = _own_chart(own_source, chart)
    side = "upper" if endpoint == "upper" else "lower"
    edge_value = _edge_domain(own_chart, side)
    exact11, edge_check = _exact11(
        _graph_nodes(own_archive),
        edge_value,
        f"{OWN_SOURCE_POINTER}.{chart}.domain.{side}",
    )
    query_check = _query_identity(f_archive, own_archive, f_source, own_source, chart)
    coords = _coord_refs(f_source, chart)
    rows = _selected_rows(f_source, chart)
    radius_shape = _radius_shape(f_archive, chart, _coordinate(f_source, chart), geometry)
    query_ref = coords["query"]
    fields: dict[str, Any] = {}
    for field in FIELDS:
        fields[field] = {
            "Z0": _source_at(
                rows[field]["Z0"], query_ref, exact11, f"F.{field}.{chart}.{endpoint}.Z0"
            ),
            "Z1": _source_at(
                rows[field]["Z1"], query_ref, exact11, f"F.{field}.{chart}.{endpoint}.ordinary_Z1"
            ),
            "ordinary_Z1_owner_identity_proved": False,
            "row_reorder_or_owner_identity_authenticated": False,
            "source_values_selected": False,
        }
    radius = _source_at(
        coords["radius"], query_ref, exact11, f"F.R.{chart}.{endpoint}.Z0"
    )
    radius_z = _radius_z(radius, geometry, radius_shape, f"F.R.{chart}.{endpoint}")
    Z = _formal("Z", source_role="ordinary_Z_coordinate")
    delta = _source_atom("delta", geometry["delta"], f"F.geometry.delta.{chart}.{endpoint}")
    delta_z = _partial_z(delta, f"F.geometry.delta_Z.{chart}.{endpoint}")
    L = _expr(
        "subtract",
        _literal("1"),
        _prod([delta, _pow(Z, "2")]),
        formula_text="L=1-delta*Z^2",
    )
    L_z = _negative(
        _sum(
            [
                _prod([delta_z, _pow(Z, "2")]),
                _prod([_literal("2"), delta, Z]),
            ]
        )
    )
    L_z["formula_text"] = "L_Z=-delta_Z*Z^2-2*delta*Z"
    R2 = _pow(radius, "2")
    Q2 = _pow(fields["Q"]["Z0"], "2")
    U2 = _pow(fields["U"]["Z0"], "2")
    R2_z = _product_z([radius, radius], [radius_z, radius_z], "(R^2)_Z=R_Z*R+R*R_Z")
    Q2_z = _product_z(
        [fields["Q"]["Z0"], fields["Q"]["Z0"]],
        [fields["Q"]["Z1"], fields["Q"]["Z1"]],
        "(Q^2)_Z=Q_Z*Q+Q*Q_Z",
    )
    U2_z = _product_z(
        [fields["U"]["Z0"], fields["U"]["Z0"]],
        [fields["U"]["Z1"], fields["U"]["Z1"]],
        "(U^2)_Z=U_Z*U+U*U_Z",
    )
    A = _prod([R2, fields["Q"]["Z0"]])
    A_z = _product_z(
        [R2, fields["Q"]["Z0"]], [R2_z, fields["Q"]["Z1"]], "A_Z=(R^2)_Z*Q+R^2*Q_Z"
    )
    B = _prod([R2, Q2])
    B_z = _product_z([R2, Q2], [R2_z, Q2_z], "B_Z=(R^2)_Z*Q^2+R^2*(Q^2)_Z")
    C = _prod([R2, U2])
    C_z = _product_z([R2, U2], [R2_z, U2_z], "C_Z=(R^2)_Z*U^2+R^2*(U^2)_Z")
    D = _prod([radius, fields["Q_y"]["Z0"]])
    D_z = _product_z(
        [radius, fields["Q_y"]["Z0"]], [radius_z, fields["Q_y"]["Z1"]], "D_Z=R_Z*Q_y+R*Q_y_Z"
    )
    K = _sum(
        [
            _expr("quotient", A, _prod([_literal("2"), L]), formula_text="R^2*Q/(2*L)"),
            _expr("quotient", B, _literal("4"), formula_text="R^2*Q^2/4"),
            _negative(
                _expr(
                    "quotient",
                    _prod([_pow(Z, "2"), C]),
                    _pow(L, "2"),
                    formula_text="Z^2*R^2*U^2/L^2",
                )
            ),
            _negative(D),
        ]
    )
    K["formula_text"] = "K=R^2*Q/(2*L)+R^2*Q^2/4-Z^2*R^2*U^2/L^2-R*Q_y"
    K_z = _sum(
        [
            _expr("quotient", A_z, _prod([_literal("2"), L]), formula_text="A_Z/(2*L)"),
            _negative(
                _expr(
                    "quotient",
                    _prod([A, L_z]),
                    _prod([_literal("2"), _pow(L, "2")]),
                    formula_text="A*L_Z/(2*L^2)",
                )
            ),
            _expr("quotient", B_z, _literal("4"), formula_text="B_Z/4"),
            _negative(
                _expr(
                    "quotient",
                    _sum([_prod([_literal("2"), Z, C]), _prod([_pow(Z, "2"), C_z])]),
                    _pow(L, "2"),
                    formula_text="(2*Z*C+Z^2*C_Z)/L^2",
                )
            ),
            _expr(
                "quotient",
                _prod([_literal("2"), _pow(Z, "2"), C, L_z]),
                _pow(L, "3"),
                formula_text="2*Z^2*C*L_Z/L^3",
            ),
            _negative(D_z),
        ]
    )
    K_z["formula_text"] = (
        "K_Z=A_Z/(2*L)-A*L_Z/(2*L^2)+B_Z/4-"
        "(2*Z*C+Z^2*C_Z)/L^2+2*Z^2*C*L_Z/L^3-D_Z"
    )
    f_edge_ref = _domain_edge_in_f(f_source, chart, side)
    endpoint = {
        "chart": chart,
        "endpoint": endpoint,
        "source_domain_edge": {
            "own": exact11,
            "own_check": edge_check,
            "F_edge_present": f_edge_ref is not None,
            "F_edge_opaque_ref": (
                _opaque_ref(f_edge_ref, f"{F_SOURCE_POINTER}.{chart}.domain.{side}", "F", "F_exact_graph")
                if f_edge_ref is not None
                else None
            ),
            "F_edge_row_read": False,
            "namespace_identity_conditional": f_edge_ref is None,
        },
        "query_identity": query_check,
        "moving_edge_term": {
            "expression": _literal("0"),
            "value": "0",
            "reason": "xi=11 is fixed for this endpoint",
            "limited_to_fixed_xi_11": True,
            "source_values_selected": False,
        },
        "radius_shape": radius_shape,
        "R": radius,
        "R_Z": radius_z,
        "fields": fields,
        "U": fields["U"]["Z0"],
        "U_Z1": fields["U"]["Z1"],
        "Q": fields["Q"]["Z0"],
        "Q_Z1": fields["Q"]["Z1"],
        "Q_y": fields["Q_y"]["Z0"],
        "Q_y_Z1": fields["Q_y"]["Z1"],
        "ordinary_Z1": {field: fields[field]["Z1"] for field in FIELDS},
        "ordinary_Z1_owner_identity_proved": False,
        "geometry": {
            "delta": delta,
            "delta_Z": delta_z,
            "mu": _source_atom("mu", geometry["mu"], f"F.geometry.mu.{chart}.{endpoint}"),
            "Rp": _source_atom("Rp", geometry["Rp"], f"F.geometry.Rp.{chart}.{endpoint}"),
            "L": L,
            "L_Z": L_z,
        },
        "intermediate": {
            "R_squared": R2,
            "R_squared_Z": R2_z,
            "Q_squared": Q2,
            "Q_squared_Z": Q2_z,
            "U_squared": U2,
            "U_squared_Z": U2_z,
            "A": A,
            "A_Z": A_z,
            "B": B,
            "B_Z": B_z,
            "C": C,
            "C_Z": C_z,
            "D": D,
            "D_Z": D_z,
        },
        "K": K,
        "K_Z": K_z,
        "source_values_selected": False,
        "source_executed": False,
    }
    structural = bool(
        edge_check.get("passed")
        and query_check.get("query_reference_rows_match")
        and all(
            rows[field][order]["status"] == "opaque_source_ref"
            for field in FIELDS
            for order in ("Z0", "Z1")
        )
        and coords["query"].get("status") == "opaque_source_ref"
        and coords["radius"].get("status") == "opaque_source_ref"
    )
    return endpoint, {
        "chart": chart,
        "endpoint": endpoint,
        "passed": structural,
        "edge_check": edge_check,
        "query_check": query_check,
        "radius_shape_verified": radius_shape.get("verified") is True,
    }


def _candidate_old_c3_reference(source: Mapping[str, Any]) -> dict[str, Any]:
    """Retain an old-C3 pointer without reading or using its endpoint rows."""

    present = False
    pointer = F_SOURCE_POINTER + ".inputs.<candidate_C3>.endpoints"
    refs = {}
    if isinstance(source.get("candidate_old_C3_endpoints"), Mapping):
        present = True
        pointer = F_SOURCE_POINTER + ".candidate_old_C3_endpoints"
    elif isinstance(source.get("inputs"), Mapping):
        inputs = source.get("inputs")
        for key in (
            "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz",
            "C3",
            "n1_pulse_pressure_C3_functions",
            "n1_pulse_pressure_C3_functions.json.gz",
        ):
            item = inputs.get(key)
            if isinstance(item, Mapping) and isinstance(item.get("endpoints"), Mapping):
                present = True
                pointer = F_SOURCE_POINTER + f".inputs[{key!r}].endpoints"
                for chart, index in (("main", 1), ("gap", 0)):
                    chart_ends = item['endpoints'].get(chart)
                    if not _is_sequence(chart_ends) or len(chart_ends) <= index or not isinstance(chart_ends[index], Mapping):
                        continue
                    fields = chart_ends[index]
                    selected = {}
                    for field in ('R_squared_q', 'R_squared_q_squared', 'R_squared_U_squared', 'R_q_y'):
                        rows = fields.get(field)
                        if _is_sequence(rows) and len(rows) >= 2:
                            selected[field] = [_opaque_ref(rows[k], f"old_C3.{chart}.{field}.Z{k}", "old_C3", "old_C3_graph") for k in (0, 1)]
                    refs[chart] = selected
                break
    return {
        "status": "candidate_only" if present else "candidate_reference_not_supplied",
        "archive_pointer": pointer,
        "present": bool(present),
        "selected_endpoint_refs": refs,
        "used_as_endpoint_source": False,
        "equal_to_new_F_endpoint_source_proved": False,
        "replacement_source": "F.leading_fields plus F.pulse_coordinates",
        "source_values_selected": False,
    }


def _obligations(
    family: Mapping[str, Any], endpoint_checks: Mapping[str, Any]
) -> list[dict[str, Any]]:
    entries = [
        ("source_family_three_fields_equal", family.get("passed") is True, "F and own source_family fields must match"),
        ("source_field_radius_owner_identity", False, "F field/radius owner identity is not inferred from equal handles"),
        ("capture_safe_endpoint", False, "source_at closures are not executed and capture safety is unproved"),
        ("ordinary_Z_commutation", False, "ordinary Z transport/commutation is outside this leaf"),
        ("numerical_interval_result", False, "no source value, numerical value, or interval is selected"),
        ("seam_cancellation", False, "main and gap terms are retained as independent expressions"),
        ("m15_closure", False, "m15 closure and error transport are outside this leaf"),
        ("ordinary_Z1_row_owner_identity", False, "F row[1] is retained but reorder/owner identity is unverified"),
        ("R_Z_source_resolution", False, "R_Z remains an unresolved formal source derivative"),
        ("delta_Z_source_resolution", False, "delta_Z remains an independent unresolved formal source derivative"),
        ("mu_nonzero", False, "the symbolic radius formula requires the same source mu to be nonzero"),
        ("L_nonzero", False, "the endpoint quotient requires L=1-delta*Z^2 to be nonzero"),
        ("radius_positive", False, "physical positivity of the two symbolic radii is not certified here"),
        (
            "selected_endpoint_checks",
            all(bool(item.get("passed")) for item in endpoint_checks.values()),
            "selected own edge and query rows must be authenticated",
        ),
    ]
    return [
        {"obligation": name, "proved": bool(passed), "status": "closed" if passed else "open", "reason": reason}
        for name, passed, reason in entries
    ]


def build(F_archive: Mapping[str, Any], own_archive: Mapping[str, Any]) -> dict[str, Any]:
    """Construct the source-bound symbolic seam from runner-authenticated maps.

    ``F_archive`` and ``own_archive`` are already decoded JSON objects.  No
    raw bytes are opened here and no source or graph evaluator is called.
    Missing or conflicting evidence yields ``status='unresolved'`` while the
    unsimplified expression records remain available for inspection.
    """

    F = _require_mapping(F_archive, "F_archive")
    own = _require_mapping(own_archive, "own_archive")
    f_source_raw = F.get(F_SOURCE_KEY)
    own_source_raw = own.get(OWN_SOURCE_KEY)
    f_source = f_source_raw if isinstance(f_source_raw, Mapping) else {}
    own_source = own_source_raw if isinstance(own_source_raw, Mapping) else {}
    f_family = _family(F, f_source)
    own_family = _family(own, own_source)
    family = _family_check(f_family, own_family)
    geometry_refs = _geometry_refs(f_source)
    endpoints: dict[str, Any] = {}
    endpoint_checks: dict[str, Any] = {}
    for chart, endpoint in (("main", "upper"), ("gap", "lower")):
        endpoint_record, check = _endpoint(
            chart,
            endpoint,
            F,
            own,
            f_source,
            own_source,
            geometry_refs,
        )
        key = f"{chart}.{endpoint}"
        endpoints[key] = endpoint_record
        endpoint_checks[key] = check
    signed_K = _expr(
        "subtract",
        endpoints["main.upper"]["K"],
        endpoints["gap.lower"]["K"],
        formula_text="signed_seam_K=K(main.upper)-K(gap.lower)",
        source_identity_proved=False,
    )
    signed_K_z = _expr(
        "subtract",
        endpoints["main.upper"]["K_Z"],
        endpoints["gap.lower"]["K_Z"],
        formula_text="signed_seam_K_Z=K_Z(main.upper)-K_Z(gap.lower)",
        source_identity_proved=False,
    )
    structural = bool(
        family.get("passed") is True
        and isinstance(f_source_raw, Mapping)
        and isinstance(own_source_raw, Mapping)
        and all(item.get("passed") is True for item in endpoint_checks.values())
        and all(
            geometry_refs[name].get("status") == "opaque_source_ref"
            for name in ("delta", "mu", "Rp")
        )
    )
    status = "constructed" if structural else "unresolved"
    obligations = _obligations(family, endpoint_checks)
    source_family_output = {
        "F": f_family,
        "own": own_family,
        "comparison": family,
    }
    result: dict[str, Any] = {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": status,
        "passed": status == "constructed",
        "constructed_means_symbolic_records_only": True,
        "source_only": True,
        "source_executed": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "interval_result": False,
        "graph_traversed": True,
        "graph_walk_scope": "selected query, own edge and current radius shape rows only",
        "graph_handles_used_as_values": False,
        "raw_sha_authenticated_by_leaf": False,
        "raw_authentication_required_of_caller": True,
        "source_family": source_family_output,
        "geometry": {
            name: {
                "source_ref": geometry_refs[name],
                "formal": _source_atom(name, geometry_refs[name], f"F.geometry.{name}"),
                "ordinary_Z_derivative": _partial_z(
                    _source_atom(name, geometry_refs[name], f"F.geometry.{name}"),
                    f"F.geometry.{name}_Z",
                ),
                "source_values_selected": False,
            }
            for name in ("delta", "mu", "Rp")
        },
        "source_rows": {
            chart: {
                "coordinates": _coord_refs(f_source, chart),
                "leading_fields": _selected_rows(f_source, chart),
                "F_query_and_own_query_same_handle": endpoints[f"{chart}.upper" if chart == "main" else f"{chart}.lower"]["query_identity"].get("same_integer_handle") is True,
            }
            for chart in CHARTS
        },
        "endpoints": endpoints,
        "signed_seam": {
            "K": signed_K,
            "K_Z": signed_K_z,
            "orientation": "main.upper-minus-gap.lower",
            "main": endpoints["main.upper"],
            "gap": endpoints["gap.lower"],
            "source_identity_proved": False,
            "seam_cancellation_proved": False,
            "source_values_selected": False,
        },
        "signed_seam_K": signed_K,
        "signed_seam_K_Z": signed_K_z,
        "candidate_old_C3_endpoints": _candidate_old_c3_reference(f_source),
        "obligations": obligations,
        "scope": {
            "source_field_radius_owner_identity": False,
            "capture_safe_endpoint": False,
            "numerical_interval_result": False,
            "seam_cancellation": False,
            "ordinary_Z_commutation": False,
            "m15_closure": False,
            "source_values_computed": False,
            "source_function_executed": False,
            "error_transport": False,
            "old_errors_removed": False,
            "full_goal_complete": False,
        },
        "audit": {
            "old_C3_endpoints_are_reference_only": True,
            "old_C3_endpoints_used_as_F": False,
            "old_C3_endpoint_equality_claimed": False,
            "protected_node_ids_not_dereferenced": sorted(PROTECTED_NODE_IDS),
            "selected_graph_rows_only": True,
            "broad_graph_traversal": False,
            "product_and_quotient_sources_retained": True,
            "main_gap_terms_cancelled": False,
            "R_Z_defaulted_zero": False,
            "delta_Z_defaulted_zero": False,
            "moving_edge_term": {
                "value": "0",
                "reason": "xi=11 is a fixed endpoint in this leaf",
                "limited_to_fixed_xi_11": True,
            },
            "endpoint_checks": endpoint_checks,
        },
    }
    return result


__all__ = [
    "CHARTS",
    "COMPONENT",
    "FIELDS",
    "FAMILY_KEYS",
    "FMainGapSignedSeamInputError",
    "ORDERS",
    "PROTECTED_NODE_IDS",
    "SCHEMA",
    "build",
]
