"""Fail-closed consumer for the native Q4 first-micro lower endpoint.

The consumer accepts an already-live Q4 object.  It does not construct a Q4
owner, call a producer/helper, evaluate a graph node, or compare a corrected
object with the native object.  The only source claim made here is the local
definition of the two ordinary initial rows::

    q4.initial["R_Q_y"].rows[k]
        = function_substitution(
            q4.charts["first_micro"]["finite_boundary"]["R_Q_y"].rows[k],
            q4.charts["first_micro"]["source_chart"]["coordinate"],
            q4.charts["first_micro"]["source_chart"]["domain"][0],
        )

for ``k = 0, 1`` in the same live function graph.  Handles are inspected as
``(graph, node)`` references and the graph is indexed only at those exact
node ids; the graph is never walked and no node expression is evaluated.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from typing import Any


SCHEMA = "original_q4_lower_substitution_v1"
UNRESOLVED_STATUS = "unresolved"
ARCHIVED_SCHEMA = "original_q4_lower_substitution_archived_evidence_v1"
ARCHIVED_KEYS = (
    "actual_common_prefix_and_native_initial_flux_ports",
    "actual_native_RQy_product_and_Z_independent_pullback",
    "original_loop_flat_branch_is_radial_not_a_Z_sample",
)

_MISSING = object()
_PRIMITIVES = (str, int, float, bool, type(None))


def _read(value: Any, key: Any) -> tuple[Any, str | None]:
    """Read one already-materialized field without calling a producer."""

    if value is _MISSING or value is None:
        return _MISSING, "parent is missing"
    if isinstance(value, Mapping):
        try:
            return (value[key], None) if key in value else (_MISSING, None)
        except Exception as exc:  # pragma: no cover - defensive mapping guard
            return _MISSING, type(exc).__name__ + ": " + str(exc)
    try:
        return value[key], None
    except (KeyError, IndexError, TypeError, AttributeError):
        pass
    except Exception as exc:  # pragma: no cover - defensive sequence guard
        return _MISSING, type(exc).__name__ + ": " + str(exc)
    try:
        return getattr(value, key), None
    except AttributeError:
        return _MISSING, None
    except Exception as exc:  # pragma: no cover - defensive property guard
        return _MISSING, type(exc).__name__ + ": " + str(exc)


def _read_index(value: Any, index: int) -> tuple[Any, str | None]:
    if value is _MISSING or value is None:
        return _MISSING, "container is missing"
    try:
        return value[index], None
    except (KeyError, IndexError, TypeError, AttributeError):
        return _MISSING, None
    except Exception as exc:  # pragma: no cover - defensive sequence guard
        return _MISSING, type(exc).__name__ + ": " + str(exc)


def _rows(jet: Any) -> tuple[tuple[Any, Any] | None, str | None]:
    """Expose only ordinary rows 0 and 1 of a FunctionJet-like object."""

    rows, detail = _read(jet, "rows")
    if rows is _MISSING:
        return None, detail or "rows field missing"
    if isinstance(rows, (str, bytes, bytearray, Mapping)):
        return None, "rows is not a sequence"
    try:
        if len(rows) < 2:
            return None, "ordinary Z0/Z1 rows are short"
        zero, zero_detail = _read_index(rows, 0)
        one, one_detail = _read_index(rows, 1)
    except Exception as exc:  # pragma: no cover - defensive sequence guard
        return None, type(exc).__name__ + ": " + str(exc)
    if zero is _MISSING or one is _MISSING:
        return None, zero_detail or one_detail or "ordinary Z0/Z1 rows are missing"
    return (zero, one), None


def _ref_parts(ref: Any, expected_graph: Any) -> tuple[dict[str, Any] | None, str | None]:
    """Read a FunctionRef's graph and node, requiring graph identity."""

    graph, graph_detail = _read(ref, "graph")
    node, node_detail = _read(ref, "node")
    if graph is _MISSING or node is _MISSING:
        return None, graph_detail or node_detail or "FunctionRef graph/node missing"
    # FunctionTransportGraph node ids are ordinary Python integer indices.
    if type(node) is not int:
        return None, "FunctionRef.node is not an integer id"
    if graph is not expected_graph:
        return None, "FunctionRef.graph is not the live q4.g object"
    return {"graph": graph, "node": node}, None


def _node_at(graph: Any, node: int) -> tuple[dict[str, Any] | None, str | None]:
    """Read exactly one graph row; this is not graph traversal."""

    nodes, detail = _read(graph, "nodes")
    if nodes is _MISSING:
        return None, detail or "graph.nodes missing"
    try:
        count = len(nodes)
    except Exception as exc:  # pragma: no cover - defensive graph guard
        return None, type(exc).__name__ + ": " + str(exc)
    if node < 0 or node >= count:
        return None, "node id is outside live graph.nodes"
    try:
        row = nodes[node]
    except Exception as exc:  # pragma: no cover - defensive graph guard
        return None, type(exc).__name__ + ": " + str(exc)
    if not isinstance(row, Mapping):
        return None, "live graph row is not a mapping"
    return dict(row), None


def _jsonable(value: Any) -> Any:
    """Make a small node row safe for metadata hashing only."""

    if isinstance(value, _PRIMITIVES):
        return value
    if isinstance(value, Mapping):
        return {str(key): _jsonable(child) for key, child in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(child) for child in value]
    # Substitution rows in the native schema contain only JSON primitives.
    # Keep an unusual extension hashable without returning the object itself.
    return {"__type__": type(value).__name__}


def _digest(row: Mapping[str, Any] | None) -> str | None:
    if row is None:
        return None
    try:
        encoded = json.dumps(
            _jsonable(row), sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False,
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()
    except Exception:  # pragma: no cover - metadata-only defensive guard
        return None


def _scalar(value: Any) -> Any:
    """Return only JSON scalar metadata; never return a live handle/object."""

    return value if isinstance(value, _PRIMITIVES) else None


def _record(
    checks: dict[str, bool],
    diagnostics: list[dict[str, Any]],
    name: str,
    passed: bool,
    *,
    reason: str | None = None,
    unresolved: bool = False,
    **extra: Any,
) -> None:
    checks[name] = bool(passed)
    item: dict[str, Any] = {
        "check": name,
        "status": "passed" if passed else (UNRESOLVED_STATUS if unresolved else "failed"),
    }
    if reason:
        item["reason"] = reason
    for key, value in extra.items():
        item[key] = _scalar(value)
    diagnostics.append(item)


def _marker(path: str, reason: str, detail: str | None = None) -> dict[str, Any]:
    result: dict[str, Any] = {"status": UNRESOLVED_STATUS, "path": path, "reason": reason}
    if detail:
        result["detail"] = detail
    return result


def _expected_substitution(
    initial_row: Mapping[str, Any] | None,
    boundary_node: int | None,
    coordinate_node: int | None,
    lower_node: int | None,
) -> dict[str, Any]:
    return {
        "operation": "function_substitution",
        "expression": boundary_node,
        "variable": coordinate_node,
        "value": lower_node,
        "Z_independent_substitution": True,
    }


def verify(q4: Any) -> dict[str, Any]:
    """Verify native Q4 first-micro initial R_Q_y substitution provenance.

    The return value is JSON-safe metadata.  A malformed or incomplete live
    object is reported as ``status='unresolved'`` and ``passed=False`` rather
    than being repaired from labels, node numbers, or attached declarations.
    """

    checks: dict[str, bool] = {}
    diagnostics: list[dict[str, Any]] = []
    unresolved: list[dict[str, Any]] = []
    matched: dict[str, dict[str, Any]] = {}
    source_ids: dict[str, Any] = {
        "coordinate": None, "lower": None, "upper": None, "sc": None, "half": None,
    }
    source_hashes: dict[str, Any] = {
        "coordinate": None, "lower": None, "upper": None, "sc": None, "half": None,
    }
    hash_checks: dict[str, bool] = {"Z0": False, "Z1": False}

    if q4 is None:
        _record(checks, diagnostics, "q4_present", False,
                reason="q4 object is None", unresolved=True)
        unresolved.append(_marker("q4", "missing_q4"))
        return {
            "schema": SCHEMA,
            "status": UNRESOLVED_STATUS,
            "passed": False,
            "checks": checks,
            "diagnostics": diagnostics,
            "unresolved": unresolved,
            "matched": matched,
            "matched_source_ids": {},
            "matched_source_hashes": {},
            "source_ids": source_ids,
            "source_hashes": source_hashes,
            "hash_checks": hash_checks,
            "source_only": True,
            "metadata_only": True,
            "numeric_values_computed": False,
            "graph_traversed": False,
            "old_owner_constructed": False,
            "old_helper_called": False,
            "cross_owner_identity_claimed": False,
            "corrected_source_identity_proved": False,
        }

    _record(checks, diagnostics, "q4_present", True)

    loop, loop_detail = _read(q4, "loop")
    q4_graph, q4_graph_detail = _read(q4, "g")
    q4_family, q4_family_detail = _read(q4, "family")
    q4_charts, q4_charts_detail = _read(q4, "charts")
    initial, initial_detail = _read(q4, "initial")
    loop_graph, loop_graph_detail = _read(loop, "g")
    loop_charts, loop_charts_detail = _read(loop, "charts")
    first_micro, first_detail = _read(q4_charts, "first_micro")
    loop_first_micro, loop_first_detail = _read(loop_charts, "first_micro")
    source_chart, source_chart_detail = _read(first_micro, "source_chart")
    finite_boundary, boundary_detail = _read(first_micro, "finite_boundary")
    boundary_rq, boundary_rq_detail = _read(finite_boundary, "R_Q_y")
    initial_rq, initial_rq_detail = _read(initial, "R_Q_y")
    coordinate, coordinate_detail = _read(source_chart, "coordinate")
    domain, domain_detail = _read(source_chart, "domain")
    lower, lower_detail = _read_index(domain, 0)
    upper, upper_detail = _read_index(domain, 1)
    sc, sc_detail = _read(loop, "sc")

    graph_identity = (
        q4_graph is not _MISSING and loop_graph is not _MISSING and q4_graph is loop_graph
    )
    _record(
        checks, diagnostics, "q4_g_is_loop_g", graph_identity,
        reason=None if graph_identity else (q4_graph_detail or loop_graph_detail or "live graph identity mismatch"),
        unresolved=q4_graph is _MISSING or loop_graph is _MISSING,
    )
    if not graph_identity:
        unresolved.append(_marker("q4.g", "q4_g_is_loop_g_failed"))

    chart_identity = (
        source_chart is not _MISSING and loop_first_micro is not _MISSING
        and source_chart is loop_first_micro
    )
    _record(
        checks, diagnostics, "q4_first_micro_source_chart_is_loop_chart", chart_identity,
        reason=None if chart_identity else (source_chart_detail or loop_first_detail or "source chart identity mismatch"),
        unresolved=source_chart is _MISSING or loop_first_micro is _MISSING,
    )
    if not chart_identity:
        unresolved.append(_marker("q4.charts.first_micro.source_chart", "source_chart_identity_failed"))

    _record(
        checks, diagnostics, "q4_first_micro_fields_present",
        first_micro is not _MISSING and finite_boundary is not _MISSING
        and boundary_rq is not _MISSING and initial_rq is not _MISSING,
        reason=first_detail or boundary_detail or boundary_rq_detail or initial_rq_detail,
        unresolved=(first_micro is _MISSING or finite_boundary is _MISSING
                    or boundary_rq is _MISSING or initial_rq is _MISSING),
    )
    if first_micro is _MISSING:
        unresolved.append(_marker("q4.charts.first_micro", "missing_first_micro_chart"))
    if boundary_rq is _MISSING:
        unresolved.append(_marker("q4.charts.first_micro.finite_boundary.R_Q_y", "missing_boundary_R_Q_y"))
    if initial_rq is _MISSING:
        unresolved.append(_marker("q4.initial.R_Q_y", "missing_initial_R_Q_y"))

    coordinate_info, coordinate_info_detail = _ref_parts(coordinate, q4_graph)
    lower_info, lower_info_detail = _ref_parts(lower, q4_graph)
    upper_info, upper_info_detail = _ref_parts(upper, q4_graph)
    for name, info, detail in (
        ("coordinate", coordinate_info, coordinate_info_detail or coordinate_detail),
        ("lower", lower_info, lower_info_detail or lower_detail),
        ("upper", upper_info, upper_info_detail or upper_detail),
    ):
        ok = info is not None
        _record(
            checks, diagnostics, "source_chart_" + name + "_is_live_function_ref", ok,
            reason=None if ok else (detail or "source chart endpoint is not a live FunctionRef"),
            unresolved=info is None,
        )
        if info is None:
            unresolved.append(_marker("source_chart." + name, "not_live_function_ref", detail))
        else:
            source_ids[name] = info["node"]
            row, row_detail = _node_at(q4_graph, info["node"])
            source_hashes[name] = _digest(row)
            _record(
                checks, diagnostics, "source_chart_" + name + "_node_in_range", row is not None,
                reason=None if row is not None else (row_detail or "source node unavailable"),
                unresolved=row is None,
                node=info["node"],
            )
            if row is None:
                unresolved.append(_marker("source_chart." + name, "source_node_unavailable", row_detail))

    # The native loop constructor defines the first-micro lower endpoint as
    # ``g.mul(self.sc, g.constant('1/2'))``.  This is checked independently of
    # the substitution row: the substitution's Z-independent flag cannot
    # certify that the endpoint really is the source-owned sc/2 node.
    sc_info, sc_info_detail = _ref_parts(sc, q4_graph)
    _record(
        checks, diagnostics, "loop_sc_is_live_function_ref", sc_info is not None,
        reason=None if sc_info is not None else (sc_info_detail or sc_detail or "loop.sc is not a live FunctionRef"),
        unresolved=sc_info is None,
    )
    sc_row = None
    sc_row_detail = None
    if sc_info is not None:
        source_ids["sc"] = sc_info["node"]
        sc_row, sc_row_detail = _node_at(q4_graph, sc_info["node"])
        source_hashes["sc"] = _digest(sc_row)
    _record(
        checks, diagnostics, "loop_sc_node_in_range", sc_row is not None,
        reason=None if sc_row is not None else (sc_row_detail or "loop.sc node unavailable"),
        unresolved=sc_row is None,
        node=None if sc_info is None else sc_info["node"],
    )
    sc_native = (
        sc_row is not None
        and sc_row.get("operation") == "current_original_source_parameter"
        and sc_row.get("name") == "sc"
    )
    _record(
        checks, diagnostics, "sc_node_is_native_source_parameter", sc_native,
        reason=None if sc_native else "live sc node metadata is not the native source parameter",
        unresolved=sc_row is None,
    )
    sc_family_matches = (
        sc_row is not None and q4_family is not _MISSING
        and sc_row.get("source_family", _MISSING) == q4_family
    )
    _record(
        checks, diagnostics, "sc_node_source_family_matches_q4_family", sc_family_matches,
        reason=None if sc_family_matches else (
            "sc source_family metadata is missing or differs from q4.family"
            if q4_family is not _MISSING else (q4_family_detail or "q4.family is missing")
        ),
        unresolved=sc_row is None or q4_family is _MISSING
                 or "source_family" not in sc_row,
    )
    sc_z_independent = sc_row is not None and sc_row.get("Z_independent") is True
    _record(
        checks, diagnostics, "sc_node_is_Z_independent", sc_z_independent,
        reason=None if sc_z_independent else "sc node Z_independent metadata is absent or false",
        unresolved=sc_row is None or "Z_independent" not in sc_row,
    )

    lower_row = None
    lower_row_detail = None
    lower_info_for_sc, lower_info_for_sc_detail = _ref_parts(lower, q4_graph)
    if lower_info_for_sc is not None:
        lower_row, lower_row_detail = _node_at(q4_graph, lower_info_for_sc["node"])
    lower_product = lower_row is not None and lower_row.get("operation") == "product"
    _record(
        checks, diagnostics, "first_micro_lower_product_operation", lower_product,
        reason=None if lower_product else (lower_row_detail or lower_info_for_sc_detail or "first-micro lower node is not product"),
        unresolved=lower_row is None,
    )
    product_arguments = lower_row.get("arguments", _MISSING) if lower_row is not None else _MISSING
    half_node = None
    if isinstance(product_arguments, (list, tuple)) and len(product_arguments) == 2:
        if sc_info is not None and product_arguments[0] == sc_info["node"]:
            half_node = product_arguments[1]
    arguments_match = (
        isinstance(product_arguments, (list, tuple))
        and sc_info is not None and len(product_arguments) == 2
        and product_arguments[0] == sc_info["node"]
    )
    _record(
        checks, diagnostics, "first_micro_lower_product_arguments", arguments_match,
        reason=None if arguments_match else "first-micro lower product arguments do not start with loop.sc",
        unresolved=lower_row is None or sc_info is None,
    )
    half_row = None
    half_row_detail = None
    if type(half_node) is int:
        source_ids["half"] = half_node
        half_row, half_row_detail = _node_at(q4_graph, half_node)
        source_hashes["half"] = _digest(half_row)
    half_exact = (
        half_row is not None
        and half_row.get("operation") == "exact_rational"
        and half_row.get("numerator") == 1
        and half_row.get("denominator") == 2
    )
    _record(
        checks, diagnostics, "first_micro_lower_half_is_exact_rational_one_half", half_exact,
        reason=None if half_exact else (half_row_detail or "sc/2 half node is not exact rational 1/2"),
        unresolved=half_node is None or half_row is None,
        node=half_node,
    )
    lower_sc_over_two = (
        lower_product and arguments_match and half_exact and sc_z_independent
        and sc_family_matches
    )
    _record(
        checks, diagnostics, "first_micro_lower_is_sc_over_two_exact", lower_sc_over_two,
        reason=None if lower_sc_over_two else "first-micro source lower is not the live native sc/2 definition",
        unresolved=lower_row is None or sc_row is None or half_row is None,
    )
    if not lower_sc_over_two:
        unresolved.append(_marker("source_chart.domain[0]", "native_sc_over_two_unresolved"))

    boundary_rows, boundary_rows_detail = _rows(boundary_rq)
    initial_rows, initial_rows_detail = _rows(initial_rq)
    boundary_rows_ok = boundary_rows is not None
    initial_rows_ok = initial_rows is not None
    _record(
        checks, diagnostics, "finite_boundary_R_Q_y_rows_Z0_Z1_present", boundary_rows_ok,
        reason=None if boundary_rows_ok else (boundary_rows_detail or "boundary ordinary rows unavailable"),
        unresolved=not boundary_rows_ok,
    )
    _record(
        checks, diagnostics, "initial_R_Q_y_rows_Z0_Z1_present", initial_rows_ok,
        reason=None if initial_rows_ok else (initial_rows_detail or "initial ordinary rows unavailable"),
        unresolved=not initial_rows_ok,
    )
    if not boundary_rows_ok:
        unresolved.append(_marker("finite_boundary.R_Q_y.rows", "ordinary_Z0_Z1_missing", boundary_rows_detail))
    if not initial_rows_ok:
        unresolved.append(_marker("initial.R_Q_y.rows", "ordinary_Z0_Z1_missing", initial_rows_detail))

    coordinate_node = coordinate_info["node"] if coordinate_info is not None else None
    lower_node = lower_info["node"] if lower_info is not None else None
    for order, label in ((0, "Z0"), (1, "Z1")):
        boundary_info = None
        initial_info = None
        boundary_row = None
        initial_row = None
        boundary_node = None
        initial_node = None
        if boundary_rows is not None:
            boundary_info, boundary_info_detail = _ref_parts(boundary_rows[order], q4_graph)
            _record(
                checks, diagnostics, label + "_boundary_ref_is_live_q4_graph",
                boundary_info is not None,
                reason=None if boundary_info is not None else (boundary_info_detail or "boundary ref graph mismatch"),
                unresolved=boundary_info is None,
            )
            if boundary_info is not None:
                boundary_node = boundary_info["node"]
                boundary_row, row_detail = _node_at(q4_graph, boundary_node)
                _record(
                    checks, diagnostics, label + "_boundary_node_in_range", boundary_row is not None,
                    reason=None if boundary_row is not None else (row_detail or "boundary node unavailable"),
                    unresolved=boundary_row is None,
                    node=boundary_node,
                )
            else:
                _record(checks, diagnostics, label + "_boundary_node_in_range", False,
                        reason="boundary FunctionRef unavailable", unresolved=True)
        else:
            _record(checks, diagnostics, label + "_boundary_ref_is_live_q4_graph", False,
                    reason="boundary ordinary rows unavailable", unresolved=True)
            _record(checks, diagnostics, label + "_boundary_node_in_range", False,
                    reason="boundary ordinary rows unavailable", unresolved=True)

        if initial_rows is not None:
            initial_info, initial_info_detail = _ref_parts(initial_rows[order], q4_graph)
            _record(
                checks, diagnostics, label + "_initial_ref_is_live_q4_graph",
                initial_info is not None,
                reason=None if initial_info is not None else (initial_info_detail or "initial ref graph mismatch"),
                unresolved=initial_info is None,
            )
            if initial_info is not None:
                initial_node = initial_info["node"]
                initial_row, row_detail = _node_at(q4_graph, initial_node)
                _record(
                    checks, diagnostics, label + "_initial_node_in_range", initial_row is not None,
                    reason=None if initial_row is not None else (row_detail or "initial node unavailable"),
                    unresolved=initial_row is None,
                    node=initial_node,
                )
            else:
                _record(checks, diagnostics, label + "_initial_node_in_range", False,
                        reason="initial FunctionRef unavailable", unresolved=True)
        else:
            _record(checks, diagnostics, label + "_initial_ref_is_live_q4_graph", False,
                    reason="initial ordinary rows unavailable", unresolved=True)
            _record(checks, diagnostics, label + "_initial_node_in_range", False,
                    reason="initial ordinary rows unavailable", unresolved=True)

        expected = _expected_substitution(initial_row, boundary_node, coordinate_node, lower_node)
        exact = True
        for field, expected_value in expected.items():
            actual_value = _MISSING if initial_row is None else initial_row.get(field, _MISSING)
            if field == "Z_independent_substitution":
                field_ok = actual_value is True and expected_value is True
            else:
                field_ok = actual_value == expected_value and actual_value is not _MISSING
            exact = exact and field_ok
            check_name = label + "_substitution_" + field
            _record(
                checks, diagnostics, check_name, field_ok,
                reason=None if field_ok else "initial node literal does not match native substitution",
                unresolved=initial_row is None or boundary_info is None
                         or coordinate_info is None or lower_info is None,
                actual=_scalar(actual_value), expected=_scalar(expected_value),
            )
        _record(
            checks, diagnostics, label + "_native_substitution_exact", exact,
            reason=None if exact else "native initial row is not the exact live Q4 substitution",
            unresolved=initial_row is None or boundary_info is None
                       or coordinate_info is None or lower_info is None,
        )

        if exact and initial_row is not None and boundary_row is not None:
            initial_hash = _digest(initial_row)
            boundary_hash = _digest(boundary_row)
            hash_checks[label] = initial_hash is not None and boundary_hash is not None
            matched[label] = {
                "ordinary_Z_order": order,
                "initial_node": initial_node,
                "boundary_node": boundary_node,
                "substitution_node": initial_node,
                "expression": _scalar(initial_row.get("expression", _MISSING)),
                "variable": _scalar(initial_row.get("variable", _MISSING)),
                "value": _scalar(initial_row.get("value", _MISSING)),
                "Z_independent_substitution": initial_row.get(
                    "Z_independent_substitution", None
                ) is True,
                "initial_node_sha256": initial_hash,
                "boundary_node_sha256": boundary_hash,
            }
        else:
            hash_checks[label] = False
        _record(
            checks, diagnostics, "runtime_hash_" + label + "_available", hash_checks[label],
            reason=None if hash_checks[label] else "exact source row hash unavailable",
            unresolved=not hash_checks[label],
        )
        if not exact:
            unresolved.append(_marker("initial.R_Q_y.rows[" + str(order) + "]",
                                     "native_substitution_mismatch"))

    both_exact = (
        checks.get("Z0_native_substitution_exact") is True
        and checks.get("Z1_native_substitution_exact") is True
    )
    _record(
        checks, diagnostics, "both_ordinary_Z_rows_native_substitution_exact", both_exact,
        reason=None if both_exact else "one or both ordinary initial rows failed literal matching",
        unresolved=not both_exact,
    )
    hashes_available = hash_checks["Z0"] and hash_checks["Z1"]
    _record(
        checks, diagnostics, "runtime_hashes_available", hashes_available,
        reason=None if hashes_available else "both matched rows require runtime SHA-256 metadata",
        unresolved=not hashes_available,
    )

    required = (
        "q4_present", "q4_g_is_loop_g", "q4_first_micro_source_chart_is_loop_chart",
        "q4_first_micro_fields_present", "source_chart_coordinate_is_live_function_ref",
        "source_chart_lower_is_live_function_ref", "source_chart_coordinate_node_in_range",
        "source_chart_lower_node_in_range", "finite_boundary_R_Q_y_rows_Z0_Z1_present",
        "initial_R_Q_y_rows_Z0_Z1_present", "Z0_native_substitution_exact",
        "Z1_native_substitution_exact", "both_ordinary_Z_rows_native_substitution_exact",
        "runtime_hashes_available", "loop_sc_is_live_function_ref", "loop_sc_node_in_range",
        "sc_node_is_native_source_parameter", "sc_node_is_Z_independent",
        "sc_node_source_family_matches_q4_family",
        "first_micro_lower_product_operation", "first_micro_lower_product_arguments",
        "first_micro_lower_half_is_exact_rational_one_half", "first_micro_lower_is_sc_over_two_exact",
    )
    passed = all(checks.get(name) is True for name in required)
    if passed:
        status = "verified"
    else:
        status = UNRESOLVED_STATUS
    return {
        "schema": SCHEMA,
        "status": status,
        "passed": passed,
        "checks": checks,
        "diagnostics": diagnostics,
        "unresolved": unresolved,
        "matched": matched,
        "matched_source_ids": {
            label: {
                "initial_node": rows.get("initial_node"),
                "boundary_node": rows.get("boundary_node"),
            }
            for label, rows in matched.items()
        },
        "matched_source_hashes": {
            label: {
                "initial_node_sha256": rows.get("initial_node_sha256"),
                "boundary_node_sha256": rows.get("boundary_node_sha256"),
            }
            for label, rows in matched.items()
        },
        "source_ids": source_ids,
        "source_hashes": source_hashes,
        "hash_checks": hash_checks,
        "hash_algorithm": "sha256",
        "ordinary_Z_orders": [0, 1],
        "source_chart": "q4.charts.first_micro.source_chart",
        "source_definition": "native Q4 live first_micro finite_boundary.R_Q_y at source_chart.domain[0]",
        "first_micro_lower_definition": "loop.sc * exact_rational(1, 2)",
        "source_only": True,
        "metadata_only": True,
        "numeric_values_computed": False,
        "graph_traversed": False,
        "old_owner_constructed": False,
        "old_helper_called": False,
        "cross_owner_identity_claimed": False,
        "corrected_source_identity_proved": False,
    }


def consume_archived_evidence(decoded_report: Any) -> dict[str, Any]:
    """Consume three already-decoded archive metadata booleans only.

    This is deliberately separate from :func:`verify`: archived evidence can
    describe a prior source scope, but it cannot establish identity with the
    current live Q4 graph or with a corrected owner.  The caller is responsible
    for authenticating the raw report hash before passing decoded metadata.
    """

    component, component_detail = _read(decoded_report, "component_result")
    audit, audit_detail = _read(component, "source_prerequisite_audit")
    checks: dict[str, bool] = {}
    diagnostics: list[dict[str, Any]] = []
    scope: dict[str, Any] = {}
    if isinstance(component, Mapping):
        # Keep the small declared domain/range fields when present; do not
        # inspect definition graphs, source nodes, or archived raw payloads.
        for key in ("ordinary_Z_orders", "source_owned_Z_cell", "domain", "range"):
            value, _ = _read(component, key)
            scope[key] = _jsonable(value) if value is not _MISSING else None

    for key in ARCHIVED_KEYS:
        value, detail = _read(audit, key)
        passed = value is True
        checks[key] = passed
        diagnostics.append({
            "check": key,
            "status": "passed" if passed else UNRESOLVED_STATUS,
            "reason": None if passed else (detail or "archived metadata field is not true"),
        })

    passed = all(checks.values()) if checks else False
    return {
        "schema": ARCHIVED_SCHEMA,
        "status": "archived_verified" if passed else UNRESOLVED_STATUS,
        "passed": passed,
        "checks": checks,
        "diagnostics": diagnostics,
        "scope": scope,
        "source_only": True,
        "metadata_only": True,
        "raw_hash_authenticated_by_consumer": False,
        "archive_graph_rebuilt": False,
        "live_q4_identity_proved": False,
        "corrected_source_identity_proved": False,
        "cross_owner_identity_claimed": False,
        "component_result_present": component is not _MISSING,
        "source_prerequisite_audit_present": audit is not _MISSING,
        "component_detail": component_detail,
        "audit_detail": audit_detail,
    }


verify_live_q4 = verify

__all__ = ["SCHEMA", "ARCHIVED_KEYS", "verify", "verify_live_q4", "consume_archived_evidence"]
