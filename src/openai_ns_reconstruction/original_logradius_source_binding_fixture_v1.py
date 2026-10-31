"""Independent synthetic checks for the log-radius source binding leaf.

``run_fixture(build)`` supplies tiny JSON-shaped archives directly to a
caller-owned builder.  The integer handles are labels only.  A guarded
``exact_graph_nodes`` mapping is included to make accidental graph access
observable; this fixture never imports a builder or evaluates a source DAG.
"""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from typing import Any, Callable


SCHEMA = "original_logradius_source_binding_fixture_v1"
SIDES = ("main.upper", "gap.lower")


class _ForbiddenGraph(Mapping[int, Any]):
    """A graph-shaped input that must never be indexed or traversed."""

    def __init__(self) -> None:
        self.accessed = False

    def __getitem__(self, key: int) -> Any:
        self.accessed = True
        raise AssertionError("fixture graph must not be accessed")

    def __iter__(self):
        self.accessed = True
        raise AssertionError("fixture graph must not be traversed")

    def __len__(self) -> int:
        self.accessed = True
        raise AssertionError("fixture graph must not be sized")


def _candidate(namespace: str, side: str, field: str, handle: int) -> dict[str, Any]:
    return {
        "operator": "source_at",
        "arguments": [
            {"operator": "opaque_source_ref", "namespace": namespace, "source_node_id": handle},
            {"operator": "exact_rational", "text": "11"},
        ],
        "namespace": namespace,
        "source_role": f"F.{side}.{field}",
        "source_node_id": handle,
        "source_values_selected": False,
    }


def _archives() -> tuple[dict[str, Any], dict[str, Any], _ForbiddenGraph, _ForbiddenGraph]:
    own_graph = _ForbiddenGraph()
    prior_graph = _ForbiddenGraph()
    own = {
        "own_main_gap_end_fields": {
            "main": {"log_radius": 101, "query": 102},
            "gap": {"log_radius": 201, "query": 202},
        },
        "exact_graph_nodes": own_graph,
    }
    prior = {
        "source_family": {"fixture": "logradius-v1", "source": "synthetic"},
        "endpoints": {
            "main.upper": {
                "R": _candidate("F.main.upper", "main.upper", "R", 301),
                "R_Z": _candidate("F.main.upper", "main.upper", "R_Z", 302),
            },
            "gap.lower": {
                "R": _candidate("F.gap.lower", "gap.lower", "R", 401),
                "R_Z": _candidate("F.gap.lower", "gap.lower", "R_Z", 402),
            },
        },
        "audit": {
            "old_log_radius_is_dimensionless_actual_source": True,
            "F_R_positive_same_reference_convention": True,
            "edge_substitution_capture_safe": False,
            "ordinary_Z1_source_row_semantics": True,
            "F_R_Z_is_derivative": False,
            "old_errors_removed": False,
        },
        "exact_graph_nodes": prior_graph,
    }
    return own, prior, own_graph, prior_graph


def _contains(value: Any, *needles: str) -> bool:
    text = repr(value)
    return all(needle in text for needle in needles)


def _endpoint(result: Mapping[str, Any], side: str) -> Mapping[str, Any]:
    endpoints = result.get("endpoints")
    item = endpoints.get(side) if isinstance(endpoints, Mapping) else None
    return item if isinstance(item, Mapping) else {}


def _call(build: Callable[..., Any], own: Mapping[str, Any], prior: Mapping[str, Any]) -> Mapping[str, Any] | None:
    try:
        result = build(own, prior)
    except Exception:
        return None
    return result if isinstance(result, Mapping) else None


def run_fixture(build: Callable[..., Any]) -> dict[str, Any]:
    """Run positive and fail-closed synthetic checks for one builder."""

    if not callable(build):
        raise TypeError("build must be callable")
    own, prior, own_graph, prior_graph = _archives()
    result = _call(build, own, prior)
    checks: dict[str, bool] = {}
    checks["source_reference_status"] = bool(
        result is not None and result.get("status") == "source-reference candidate"
    )
    checks["graphs_not_accessed"] = not own_graph.accessed and not prior_graph.accessed

    main = _endpoint(result or {}, "main.upper")
    gap = _endpoint(result or {}, "gap.lower")
    checks["both_independent_endpoints"] = bool(
        main.get("source_at_old_log_radius", {}).get("operator") == "source_at"
        and gap.get("source_at_old_log_radius", {}).get("operator") == "source_at"
    )

    def source_at_handles(endpoint: Mapping[str, Any], log_handle: int, query_handle: int) -> bool:
        node = endpoint.get("source_at_old_log_radius")
        args = node.get("arguments") if isinstance(node, Mapping) else None
        if not isinstance(args, list) or len(args) != 3:
            return False
        return bool(
            args[0].get("source_node_id") == log_handle
            and args[1].get("source_node_id") == query_handle
            and args[2].get("typed_exact_rational_11") is True
        )

    checks["exact_old_main_binding"] = source_at_handles(main, 101, 102)
    checks["exact_old_gap_binding"] = source_at_handles(gap, 201, 202)
    checks["D_references_both_actual_sources"] = bool(
        _contains(main.get("D_native"), "101", "102", "301")
        and _contains(gap.get("D_native"), "201", "202", "401")
    )
    checks["D_Z_references_both_R_Z_sources"] = bool(
        _contains(main.get("D_Z"), "101", "302", "301", "rho_ref_Z")
        and _contains(gap.get("D_Z"), "201", "402", "401", "rho_ref_Z")
    )
    checks["explicit_rho_terms"] = bool(
        main.get("rho_ref", {}).get("operator") == "rho_ref"
        and main.get("rho_ref_Z", {}).get("operator") == "rho_ref_Z"
        and main.get("D", {}).get("rho_ref", {}).get("operator") == "rho_ref"
        and main.get("D_Z", {}).get("rho_ref_Z", {}).get("operator") == "rho_ref_Z"
    )
    checks["sides_not_interchangeable"] = bool(
        not _contains(main.get("D_native"), "201", "202", "401")
        and not _contains(gap.get("D_native"), "101", "102", "301")
    )
    checks["old_errors_not_deleted"] = bool(
        result is not None
        and result.get("old_errors_removed") is False
        and result.get("historical_error_terms_deleted") is False
    )
    checks["not_evaluated"] = bool(
        result is not None
        and result.get("source_values_computed") is False
        and result.get("numeric_values_computed") is False
    )

    missing = deepcopy(own)
    del missing["own_main_gap_end_fields"]["gap"]["query"]
    missing_graphs = (missing.get("exact_graph_nodes"), prior.get("exact_graph_nodes"))
    checks["missing_named_path_fails_closed"] = _call(build, missing, prior) is None
    checks["missing_path_did_not_touch_graph"] = not any(graph.accessed for graph in missing_graphs if graph is not None)

    selected = deepcopy(prior)
    selected["endpoints"]["main.upper"]["R"]["source_values_selected"] = True
    checks["selected_candidate_fails_closed"] = _call(build, own, selected) is None

    checks = {name: bool(value) for name, value in checks.items()}
    return {
        "schema": SCHEMA,
        "synthetic_only": True,
        "passed": bool(all(checks.values())),
        "checks": checks,
        "source_values_computed": False,
        "graph_traversed": False,
        "old_errors_removed": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
