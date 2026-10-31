"""Source-only extraction of the live original Q4 lower chart.

``extract`` consumes an already-live Q4 object.  It only reads the source
fields documented by the original Q4/C5 constructors, preserves the opaque
references, and never constructs an ancestor, evaluates a graph, or compares
source and corrected owners.  Missing fields remain explicit unresolved
records.
"""

from __future__ import annotations

from typing import Any


SCHEMA = "original_q4_lower_live_bundle_v1"
UNRESOLVED_STATUS = "unresolved"
_MISSING = object()


def _unresolved(path: str, reason: str, detail: str | None = None) -> dict[str, Any]:
    result: dict[str, Any] = {"status": UNRESOLVED_STATUS, "path": path, "reason": reason}
    if detail:
        result["detail"] = detail
    return result


def _is_unresolved(value: Any) -> bool:
    return isinstance(value, dict) and value.get("status") == UNRESOLVED_STATUS


def _read(value: Any, key: Any) -> tuple[Any, str | None]:
    """Read one field without invoking a producer or graph operation."""

    if value is _MISSING or _is_unresolved(value):
        return _MISSING, "parent is unresolved"
    if value is None:
        return _MISSING, "parent is None"
    try:
        if isinstance(value, dict):
            return (value[key], None) if key in value else (_MISSING, None)
    except Exception as exc:  # pragma: no cover - defensive unusual mapping
        return _MISSING, type(exc).__name__ + ": " + str(exc)
    try:
        return value[key], None
    except (KeyError, IndexError, TypeError, AttributeError):
        pass
    except Exception as exc:  # pragma: no cover - defensive wrapper
        return _MISSING, type(exc).__name__ + ": " + str(exc)
    try:
        return getattr(value, key), None
    except AttributeError:
        return _MISSING, None
    except Exception as exc:  # pragma: no cover - defensive property
        return _MISSING, type(exc).__name__ + ": " + str(exc)


def _read_index(value: Any, index: int) -> tuple[Any, str | None]:
    if value is _MISSING or _is_unresolved(value) or value is None:
        return _MISSING, "row container is unresolved"
    try:
        return value[index], None
    except (KeyError, IndexError, TypeError, AttributeError):
        return _MISSING, None
    except Exception as exc:  # pragma: no cover - defensive sequence wrapper
        return _MISSING, type(exc).__name__ + ": " + str(exc)


def _capture(value: Any, key: Any, path: str, unresolved: list[dict[str, Any]]) -> Any:
    result, detail = _read(value, key)
    if result is not _MISSING:
        return result
    marker = _unresolved(path, "missing_field", detail)
    unresolved.append(marker)
    return marker


def _capture_optional(value: Any, key: Any, path: str) -> Any:
    result, detail = _read(value, key)
    return result if result is not _MISSING else _unresolved(
        path, "optional_field_not_exposed", detail
    )


def _ordinary_z01(value: Any, path: str, unresolved: list[dict[str, Any]]) -> dict[str, Any]:
    """Return raw ordinary Z0/Z1 rows; no handle is converted to a value."""

    labels = ["Z0", "Z1"]
    if value is _MISSING or _is_unresolved(value):
        marker = _unresolved(path, "missing_ordinary_rows")
        unresolved.append(marker)
        return {"raw": value if value is not _MISSING else marker, "rows": marker,
                "Z0": marker, "Z1": marker, "ordinary_Z_orders": [0, 1],
                "ordinary_Z_labels": labels}

    rows_value, rows_detail = _read(value, "rows")
    if rows_value is _MISSING:
        z0, z0_detail = _read(value, "Z0")
        z1, z1_detail = _read(value, "Z1")
        if z0 is _MISSING:
            z0, z0_detail = _read(value, 0)
        if z1 is _MISSING:
            z1, z1_detail = _read(value, 1)
        if z0 is not _MISSING and z1 is not _MISSING:
            return {"raw": value, "rows": (z0, z1), "Z0": z0, "Z1": z1,
                    "ordinary_Z_orders": [0, 1], "ordinary_Z_labels": labels}
        marker = _unresolved(path, "ordinary_Z0_Z1_rows_missing",
                             z0_detail or z1_detail or rows_detail)
        unresolved.append(marker)
        return {"raw": value, "rows": marker, "Z0": marker, "Z1": marker,
                "ordinary_Z_orders": [0, 1], "ordinary_Z_labels": labels}

    z0, z0_detail = _read_index(rows_value, 0)
    z1, z1_detail = _read_index(rows_value, 1)
    if z0 is _MISSING or z1 is _MISSING:
        marker = _unresolved(path, "ordinary_Z0_Z1_rows_short", z0_detail or z1_detail)
        unresolved.append(marker)
        return {"raw": value, "rows": rows_value,
                "Z0": marker if z0 is _MISSING else z0,
                "Z1": marker if z1 is _MISSING else z1,
                "ordinary_Z_orders": [0, 1], "ordinary_Z_labels": labels}
    return {"raw": value, "rows": rows_value, "Z0": z0, "Z1": z1,
            "ordinary_Z_orders": [0, 1], "ordinary_Z_labels": labels}


def _direct_graph(value: Any) -> tuple[Any, str] | None:
    """Read only a reference's graph owner, never ``graph.nodes``."""

    if value is _MISSING or _is_unresolved(value) or value is None:
        return None
    try:
        graph = getattr(value, "graph")
    except Exception:
        graph = _MISSING
    if graph is not _MISSING and graph is not None:
        return graph, "graph"
    # Owner wrappers may spell the graph field ``g``; require ``node`` so a
    # chart metadata object's ordinary ``g`` field is not mistaken for a ref.
    try:
        node = getattr(value, "node")
    except Exception:
        node = _MISSING
    if node is not _MISSING:
        try:
            graph = getattr(value, "g")
        except Exception:
            graph = _MISSING
        if graph is not _MISSING and graph is not None:
            return graph, "g"
    return None


def _iter_graphs(value: Any, path: str, seen: set[int] | None = None):
    """Walk selected containers only; this is not graph structure traversal."""

    if seen is None:
        seen = set()
    if value is _MISSING or _is_unresolved(value) or value is None:
        return
    direct = _direct_graph(value)
    if direct is not None:
        yield path, direct[0], direct[1]
        return
    identity = id(value)
    if identity in seen:
        return
    seen.add(identity)
    if isinstance(value, dict):
        for key, child in value.items():
            yield from _iter_graphs(child, path + "." + str(key), seen)
        return
    if isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            yield from _iter_graphs(child, path + "[" + str(index) + "]", seen)
        return
    rows_value, _ = _read(value, "rows")
    if rows_value is not _MISSING and rows_value is not value:
        yield from _iter_graphs(rows_value, path + ".rows", seen)


def _graph_identity_check(name: str, value: Any, expected_graph: Any,
                          unresolved: list[dict[str, Any]],
                          diagnostics: list[dict[str, Any]]) -> bool:
    if expected_graph is _MISSING or _is_unresolved(expected_graph):
        marker = _unresolved(name, "expected_graph_missing")
        unresolved.append(marker)
        diagnostics.append({"check": name, "status": UNRESOLVED_STATUS,
                            "reason": marker["reason"]})
        return False
    refs = list(_iter_graphs(value, name))
    if not refs:
        marker = _unresolved(name, "no_graph_reference_exposed")
        unresolved.append(marker)
        diagnostics.append({"check": name, "status": UNRESOLVED_STATUS,
                            "reason": marker["reason"]})
        return False
    passed = all(graph is expected_graph for _, graph, _ in refs)
    diagnostics.append({"check": name, "status": "passed" if passed else "failed",
                        "reference_count": len(refs), "identity_only": True})
    return passed


def _same_object_check(name: str, left: Any, right: Any,
                       unresolved: list[dict[str, Any]],
                       diagnostics: list[dict[str, Any]]) -> bool:
    if left is _MISSING or right is _MISSING or _is_unresolved(left) or _is_unresolved(right):
        marker = _unresolved(name, "identity_operand_missing")
        unresolved.append(marker)
        diagnostics.append({"check": name, "status": UNRESOLVED_STATUS,
                            "reason": marker["reason"]})
        return False
    passed = left is right
    diagnostics.append({"check": name, "status": "passed" if passed else "failed",
                        "identity_only": True})
    return passed


def _initial_lower_evidence(q4: Any, initial: Any, first_micro: Any,
                            source_domain: Any,
                            unresolved: list[dict[str, Any]],
                            diagnostics: list[dict[str, Any]]) -> tuple[Any, bool]:
    """Record a lower-source ref without claiming the substitution identity.

    The real Q4 constructor evaluates ``initial`` at ``chart['domain'][0]``
    but does not retain that argument in ``initial``.  Consequently a live
    object remains unresolved here, even when it carries a matching witness;
    labels, numerical equality, and graph node ids are not evidence.
    """

    lower, lower_detail = _read_index(source_domain, 0)
    if lower is _MISSING:
        marker = _unresolved("first_micro.initial_lower",
                             "source_chart_domain_lower_missing", lower_detail)
        unresolved.append(marker)
        diagnostics.append({"check": "first_micro_initial_lower_is_source_chart_domain_lower",
                            "status": UNRESOLVED_STATUS, "reason": marker["reason"]})
        return marker, False

    candidates: list[Any] = []
    for owner, key, path in (
        (q4, "initial_lower", "q4.initial_lower"),
        (q4, "initial_lower_source", "q4.initial_lower_source"),
        (q4, "initial_endpoint_lower", "q4.initial_endpoint_lower"),
        (first_micro, "initial_lower", "q4.charts.first_micro.initial_lower"),
        (initial, "lower", "q4.initial.lower"),
        (initial, "lower_source", "q4.initial.lower_source"),
        (initial, "evaluation_coordinate", "q4.initial.evaluation_coordinate"),
    ):
        candidate, _ = _read(owner, key)
        if candidate is not _MISSING:
            candidates.append(candidate)
            diagnostics.append({"check": path, "status": "observed_explicit_reference"})

    if not candidates:
        marker = _unresolved(
            "first_micro.initial_lower",
            "initial_lower_identity_not_exposed",
            "constructor endpoint argument is not retained by the live Q4 object",
        )
        unresolved.append(marker)
        diagnostics.append({"check": "first_micro_initial_lower_is_source_chart_domain_lower",
                            "status": UNRESOLVED_STATUS, "reason": marker["reason"]})
        return marker, False
    same_reference = all(candidate is lower for candidate in candidates)
    # Even an attached reference with the right identity does not prove that
    # the initial flux was actually formed by substitution at that reference.
    # The constructor's defining operation must be authenticated separately.
    marker = _unresolved("first_micro.initial_lower",
                         "initial_substitution_binding_not_inspected")
    unresolved.append(marker)
    diagnostics.append({"check": "first_micro_initial_lower_is_source_chart_domain_lower",
                        "status": UNRESOLVED_STATUS,
                        "candidate_is_domain_lower": same_reference,
                        "candidate_count": len(candidates),
                        "reason": marker["reason"]})
    return candidates[0], False


def _internal_optional_identity(value: Any, label: str) -> dict[str, Any]:
    """Inspect an optional corrected side internally, never against source."""

    if value is None:
        return {"present": False, "graph_identity": None, "status": "not_supplied",
                "cross_source_identity_checked": False}
    refs = list(_iter_graphs(value, label))
    if not refs:
        return {"present": True, "graph_identity": False, "status": UNRESOLVED_STATUS,
                "reason": "no_graph_reference_exposed", "cross_source_identity_checked": False}
    anchor = refs[0][1]
    passed = all(graph is anchor for _, graph, _ in refs)
    return {"present": True, "graph_identity": passed,
            "status": "passed" if passed else "failed", "reference_count": len(refs),
            "identity_only": True, "cross_source_identity_checked": False}


def extract(q4: Any, corrected_packet: Any = None,
            corrected_view: Any = None) -> dict[str, Any]:
    """Return a source-only bundle from an already-live Q4 instance."""

    unresolved: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []

    loop = _capture(q4, "loop", "q4.loop", unresolved)
    q4_graph = _capture(q4, "g", "q4.g", unresolved)
    q4_background = _capture(q4, "background", "q4.background", unresolved)
    q4_charts = _capture(q4, "charts", "q4.charts", unresolved)
    q4_initial = _capture(q4, "initial", "q4.initial", unresolved)

    loop_graph = _capture(loop, "g", "q4.loop.g", unresolved)
    loop_background = _capture(loop, "background", "q4.loop.background", unresolved)
    loop_charts = _capture(loop, "charts", "q4.loop.charts", unresolved)
    loop_N = _capture(loop, "N", "q4.loop.N", unresolved)

    first_micro = _capture(q4_charts, "first_micro", "q4.charts.first_micro", unresolved)
    source_chart = _capture(first_micro, "source_chart",
                            "q4.charts.first_micro.source_chart", unresolved)
    loop_first_micro = _capture(loop_charts, "first_micro",
                                "q4.loop.charts.first_micro", unresolved)
    physical = _capture(first_micro, "physical",
                        "q4.charts.first_micro.physical", unresolved)
    finite_boundary = _capture(first_micro, "finite_boundary",
                               "q4.charts.first_micro.finite_boundary", unresolved)
    finite_boundary_R_Q_y = _capture(
        finite_boundary, "R_Q_y", "q4.charts.first_micro.finite_boundary.R_Q_y", unresolved
    )
    initial_R_Q_y = _capture(q4_initial, "R_Q_y", "q4.initial.R_Q_y", unresolved)

    physical_fields = {
        key: _capture(physical, key, "q4.charts.first_micro.physical." + key, unresolved)
        for key in ("U", "m", "U_y", "Q_y")
    }
    source_chart_fields = {
        key: _capture(source_chart, key,
                      "q4.charts.first_micro.source_chart." + key, unresolved)
        for key in ("coordinate", "domain", "radius", "modified_fields")
    }
    loop_first_fields = {
        key: _capture(loop_first_micro, key,
                      "q4.loop.charts.first_micro." + key, unresolved)
        for key in ("incoming", "history")
    }
    delta = _capture(q4_background, "delta", "q4.background.delta", unresolved)
    Pstar = _capture(q4_background, "Pstar", "q4.background.Pstar", unresolved)
    initial_rows = _ordinary_z01(initial_R_Q_y, "q4.initial.R_Q_y", unresolved)
    initial_lower_source, initial_lower_identity = _initial_lower_evidence(
        q4, q4_initial, first_micro, source_chart_fields["domain"], unresolved, diagnostics
    )

    checks: dict[str, bool] = {}
    checks["q4_g_is_loop_g"] = _same_object_check(
        "q4_g_is_loop_g", q4_graph, loop_graph, unresolved, diagnostics
    )
    checks["q4_background_is_loop_background"] = _same_object_check(
        "q4_background_is_loop_background", q4_background, loop_background,
        unresolved, diagnostics
    )
    checks["q4_first_micro_source_chart_is_loop_first_micro"] = _same_object_check(
        "q4_first_micro_source_chart_is_loop_first_micro", source_chart,
        loop_first_micro, unresolved, diagnostics
    )
    checks["first_micro_initial_lower_is_source_chart_domain_lower"] = initial_lower_identity

    graph_fields = {
        "q4_loop_N_graph_identity": loop_N,
        "background_delta_graph_identity": delta,
        "background_Pstar_graph_identity": Pstar,
        "physical_U_graph_identity": physical_fields["U"],
        "physical_m_graph_identity": physical_fields["m"],
        "physical_U_y_graph_identity": physical_fields["U_y"],
        "physical_Q_y_graph_identity": physical_fields["Q_y"],
        "finite_boundary_R_Q_y_graph_identity": finite_boundary_R_Q_y,
        "initial_R_Q_y_Z0_graph_identity": initial_rows["Z0"],
        "initial_R_Q_y_Z1_graph_identity": initial_rows["Z1"],
        "source_chart_coordinate_graph_identity": source_chart_fields["coordinate"],
        "source_chart_domain_graph_identity": source_chart_fields["domain"],
        "source_chart_radius_graph_identity": source_chart_fields["radius"],
        "source_chart_modified_fields_graph_identity": source_chart_fields["modified_fields"],
        "loop_first_micro_incoming_graph_identity": loop_first_fields["incoming"],
        "loop_first_micro_history_graph_identity": loop_first_fields["history"],
    }
    for name, value in graph_fields.items():
        checks[name] = _graph_identity_check(name, value, q4_graph, unresolved, diagnostics)

    checks["first_micro_physical_fields_present"] = all(
        not _is_unresolved(value) for value in physical_fields.values()
    )
    checks["first_micro_finite_boundary_R_Q_y_present"] = not _is_unresolved(
        finite_boundary_R_Q_y
    )
    checks["initial_R_Q_y_ordinary_Z0_Z1_present"] = not any(
        _is_unresolved(initial_rows[key]) for key in ("Z0", "Z1")
    )
    checks["source_chart_fields_present"] = all(
        not _is_unresolved(value) for value in source_chart_fields.values()
    )
    checks["loop_first_micro_incoming_history_present"] = all(
        not _is_unresolved(value) for value in loop_first_fields.values()
    )

    # Keep the original error objects as-is.  Error fields are optional in the
    # live classes, so their absence is visible in the bundle but does not
    # manufacture a source failure.
    raw_errors = {
        "q4": _capture_optional(q4, "errors", "q4.errors"),
        "loop": _capture_optional(loop, "errors", "q4.loop.errors"),
        "background": _capture_optional(q4_background, "errors", "q4.background.errors"),
        "q4_first_micro": _capture_optional(first_micro, "errors",
                                              "q4.charts.first_micro.errors"),
        "loop_first_micro": _capture_optional(loop_first_micro, "errors",
                                                "q4.loop.charts.first_micro.errors"),
    }
    optional_internal = {
        "corrected_packet": _internal_optional_identity(corrected_packet, "corrected_packet"),
        "corrected_view": _internal_optional_identity(corrected_view, "corrected_view"),
    }

    first_micro_bundle = {
        "chart": first_micro, "physical": physical, "physical_fields": physical_fields,
        "finite_boundary": finite_boundary, "finite_boundary_R_Q_y": finite_boundary_R_Q_y,
        "initial_R_Q_y": initial_R_Q_y, "initial_R_Q_y_rows": initial_rows,
        "source_chart": source_chart, "source_chart_fields": source_chart_fields,
        "loop_chart": loop_first_micro, "loop_incoming": loop_first_fields["incoming"],
        "loop_history": loop_first_fields["history"],
        "initial_lower_source": initial_lower_source,
        "initial_lower_identity_proved": initial_lower_identity,
    }
    return {
        "schema": SCHEMA, "source_only": True, "metadata_only": True,
        "constructed": False, "ancestor_calls": 0, "construction_calls": 0,
        "recompute_calls": 0, "numeric_values_computed": False,
        "graph_handles_treated_as_values": False,
        "q4": q4, "loop": loop, "background": q4_background, "graph": q4_graph,
        "first_micro": first_micro_bundle, "physical_fields": physical_fields,
        "finite_boundary": {"raw": finite_boundary, "R_Q_y": finite_boundary_R_Q_y},
        "initial": {"raw": q4_initial, "R_Q_y": initial_R_Q_y,
                     "R_Q_y_rows": initial_rows},
        "source_chart": {"raw": source_chart, **source_chart_fields},
        "loop_first_micro": {"raw": loop_first_micro,
                              "incoming": loop_first_fields["incoming"],
                              "history": loop_first_fields["history"]},
        "background_fields": {"raw": q4_background, "delta": delta, "Pstar": Pstar},
        "loop_fields": {"raw": loop, "N": loop_N},
        "errors": raw_errors, "identity_checks": checks, "checks": checks,
        "unresolved": unresolved, "diagnostics": diagnostics,
        "optional_internal_views": optional_internal,
        "cross_endpoint_identity_claimed": False,
        "native_corrected_identity_proved": False, "full_goal_complete": False,
    }


extract_source_bundle = extract

__all__ = ["SCHEMA", "extract", "extract_source_bundle"]

