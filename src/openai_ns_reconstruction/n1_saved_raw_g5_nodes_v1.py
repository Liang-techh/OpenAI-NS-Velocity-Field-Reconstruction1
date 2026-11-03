"""Read-only references into the saved raw n=1 Volterra g5 graph.

The caller supplies one already authenticated and decoded graph mapping. This
leaf performs no source import, producer replay, numerical evaluation, interval
comparison, AST walk, or graph comparison. It indexes only the saved C3
g['5'] rows for the core and first systems and then selects the ordinary
Z0/Z1 positions from those rows.

The ordinary labels are a structural interpretation of the saved C3 positions
0 and 1. They are not a live-owner check and are not numerical values.
In particular, the solution-function handles, axis owner, m15, cancellation,
and a full Kchi construction are deliberately outside this consumer's claims.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


__all__ = [
    "SCHEMA",
    "ORDINARY_Z_LABELS",
    "SavedRawG5Error",
    "consume_saved_raw_g5",
    "read_saved_raw_g5",
    "extract_saved_raw_g5",
]


SCHEMA = "n1_saved_raw_g5_nodes_v1_result"
ORDINARY_Z_LABELS = ("Z0", "Z1")
_CHARTS = ("core", "first")
_SOURCE_PATH_PREFIX = ("inner_Volterra_bindings", "systems")
_MISSING = object()


class SavedRawG5Error(ValueError):
    """A saved raw graph does not satisfy this narrow consumer contract."""


def _json_copy(value: Any, label: str) -> Any:
    """Copy JSON-shaped values without interpreting graph records."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Mapping):
        copied: dict[Any, Any] = {}
        for key, item in value.items():
            if not isinstance(key, (str, int, float, bool)):
                raise SavedRawG5Error(f"{label}: non-JSON mapping key {key!r}")
            copied[key] = _json_copy(item, f"{label}[{key!r}]")
        return copied
    if isinstance(value, (list, tuple)):
        return [_json_copy(item, f"{label}[{index}]") for index, item in enumerate(value)]
    raise SavedRawG5Error(f"{label}: value is not JSON-shaped")


def _required_metadata(graph: Mapping[str, Any], names: tuple[str, ...], label: str) -> Any:
    """Read one authenticated runner metadata field without deriving it."""

    for name in names:
        value = graph.get(name, _MISSING)
        if value is not _MISSING:
            if value is None:
                raise SavedRawG5Error(f"{label}: {name} is null")
            return _json_copy(value, label)
    raise SavedRawG5Error(f"{label}: authenticated {names[0]} is required")


def _node_id(value: Any, *, label: str, node_count: int) -> int:
    """Validate one saved node index; bool is intentionally not an integer id."""

    if isinstance(value, bool) or not isinstance(value, int):
        raise SavedRawG5Error(f"{label}: node id must be an int, got {value!r}")
    if value < 0 or value >= node_count:
        raise SavedRawG5Error(
            f"{label}: node id {value!r} is outside exact_graph_nodes[0:{node_count}]"
        )
    return value


def _source_path(chart: str) -> list[str]:
    return [*_SOURCE_PATH_PREFIX, chart, "g", "5"]


def _g5_nodes(bindings: Mapping[str, Any], chart: str) -> list[Any]:
    """Directly read one encoded C3 g5 row, with no graph traversal."""

    try:
        systems = bindings["systems"]
        system = systems[chart]
        source = system["g"]
        row = source["5"]
    except (KeyError, TypeError) as exc:
        raise SavedRawG5Error(
            f"{chart}: missing direct inner_Volterra_bindings.systems.{chart}.g['5']"
        ) from exc
    if not isinstance(row, list) or len(row) != 4:
        raise SavedRawG5Error(
            f"{chart}: saved C3 g['5'] must be a node-id list of length 4"
        )
    return row


def consume_saved_raw_g5(graph: Mapping[str, Any]) -> dict[str, Any]:
    """Return exact raw-record references for saved core/first g5 Z0 and Z1.

    graph must already carry runner authentication metadata. The function
    never computes a graph hash or infers a file path, and it never substitutes
    the separate solution functions handles for the two raw g['5'] rows.
    """

    if not isinstance(graph, Mapping):
        raise SavedRawG5Error("authenticated graph mapping is required")

    source_family = _required_metadata(graph, ("source_family",), "source_family")
    graph_hash = _required_metadata(graph, ("graph_hash", "graph_sha256"), "graph_hash")
    file_path = _required_metadata(
        graph, ("file_path", "graph_file_path", "graph_path"), "file_path"
    )

    nodes = graph.get("exact_graph_nodes", _MISSING)
    if not isinstance(nodes, list) or not nodes:
        raise SavedRawG5Error("exact_graph_nodes must be a non-empty decoded list")

    bindings = graph.get("inner_Volterra_bindings", _MISSING)
    if not isinstance(bindings, Mapping):
        raise SavedRawG5Error("inner_Volterra_bindings mapping is required")

    selected: dict[str, Any] = {}
    for chart in _CHARTS:
        c3_ids = _g5_nodes(bindings, chart)
        saved_ids = [
            _node_id(c3_ids[position], label=f"{chart}.g5[{position}]", node_count=len(nodes))
            for position in range(4)
        ]

        ordinary: dict[str, Any] = {}
        for position, z_label in enumerate(ORDINARY_Z_LABELS):
            node_id = saved_ids[position]
            # This is the sole exact_graph_nodes selection performed by the
            # consumer for this chart and ordinary Z position.
            raw_record = _json_copy(nodes[node_id], f"{chart}.g5.{z_label}.raw_record")
            ordinary[z_label] = {
                "ordinary_Z_label": z_label,
                "saved_C3_position": position,
                "node_id": node_id,
                "node_path": ["exact_graph_nodes", node_id],
                "source_path": _source_path(chart),
                "raw_record": raw_record,
            }

        selected[chart] = {
            "source_path": _source_path(chart),
            "saved_C3_node_ids": saved_ids,
            "ordinary_Z": ordinary,
        }

    return {
        "schema": SCHEMA,
        "passed": True,
        "source_family": source_family,
        "graph_hash": graph_hash,
        "file_path": file_path,
        "graph_node_count": len(nodes),
        "g5": selected,
        "interpretation": {
            "ordinary_Z_labels": list(ORDINARY_Z_LABELS),
            "ordinary_Z_source": "saved C3 positions 0 and 1",
            "references_are_symbolic_graph_nodes": True,
            "point_values_computed": False,
            "live_owner_verified": False,
            "function_solution_handles_used": False,
            "axis_owner_claimed": False,
            "m15_claimed": False,
            "cancellation_claimed": False,
            "full_Kchi_claimed": False,
            "intervals_compared": False,
        },
    }


def read_saved_raw_g5(graph: Mapping[str, Any]) -> dict[str, Any]:
    """Named alias for callers that describe the operation as a read."""

    return consume_saved_raw_g5(graph)


def extract_saved_raw_g5(graph: Mapping[str, Any]) -> dict[str, Any]:
    """Named alias retained for source-consumer discovery code."""

    return consume_saved_raw_g5(graph)
