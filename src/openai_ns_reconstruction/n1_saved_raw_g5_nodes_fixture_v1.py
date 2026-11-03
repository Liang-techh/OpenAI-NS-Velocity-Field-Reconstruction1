"""Synthetic-only fixture for n1_saved_raw_g5_nodes_v1.

The fixture constructs a new fictional decoded graph in memory. It never reads
the saved archive, imports a producer, evaluates a mathematical function, or
uses an existing graph record. The supplied module is the new consumer under
review and is the only implementation called by this fixture.
"""

from __future__ import annotations

import copy
import json
from typing import Any


__all__ = ["SCHEMA", "run_fixture", "run"]


SCHEMA = "n1_saved_raw_g5_nodes_fixture_v1"


def _synthetic_graph() -> dict[str, Any]:
    """Return a fresh fictional graph with distinct core/first source rows."""

    nodes = [
        {"record": "fictional-node-0", "source": {"sector": "unused", "order": 0}},
        {"record": "fictional-node-1", "source": {"sector": "unused", "order": 1}},
        {
            "record": "fictional-core-g5-Z0",
            "source": {"chart": "core", "kind": "raw-g5", "ordinary_Z": "Z0"},
            "payload": ["core", 0, "symbolic"],
        },
        {
            "record": "fictional-first-g5-Z0",
            "source": {"chart": "first", "kind": "raw-g5", "ordinary_Z": "Z0"},
            "payload": ["first", 0, "symbolic"],
        },
        {
            "record": "fictional-core-g5-Z1",
            "source": {"chart": "core", "kind": "raw-g5", "ordinary_Z": "Z1"},
            "payload": ["core", 1, "symbolic"],
        },
        {
            "record": "fictional-core-g5-Z2",
            "source": {"chart": "core", "kind": "raw-g5", "ordinary_Z": "Z2"},
            "payload": ["core", 2, "symbolic"],
        },
        {
            "record": "fictional-core-g5-Z3",
            "source": {"chart": "core", "kind": "raw-g5", "ordinary_Z": "Z3"},
            "payload": ["core", 3, "symbolic"],
        },
        {
            "record": "fictional-first-g5-Z1",
            "source": {"chart": "first", "kind": "raw-g5", "ordinary_Z": "Z1"},
            "payload": ["first", 1, "symbolic"],
        },
        {
            "record": "fictional-first-g5-Z2",
            "source": {"chart": "first", "kind": "raw-g5", "ordinary_Z": "Z2"},
            "payload": ["first", 2, "symbolic"],
        },
        {
            "record": "fictional-first-g5-Z3",
            "source": {"chart": "first", "kind": "raw-g5", "ordinary_Z": "Z3"},
            "payload": ["first", 3, "symbolic"],
        },
    ]
    return {
        "source_family": {"name": "fictional-saved-family", "revision": 1},
        "graph_hash": "fictional-graph-sha256",
        "file_path": "fixture/fictional_saved_graph.json.gz",
        "exact_graph_nodes": nodes,
        "inner_Volterra_bindings": {
            "systems": {
                "core": {"g": {"5": [2, 4, 5, 6]}},
                "first": {"g": {"5": [3, 7, 8, 9]}},
            },
            # Deliberately unrelated handles: the consumer must not read these
            # solution-function entries as raw g5 nodes.
            "functions": {
                "core": {"F1": [0, 1, 2, 3]},
                "first": {"F1": [4, 5, 6, 7]},
            },
        },
    }


def _rejected(consumer: Any, graph: dict[str, Any]) -> bool:
    try:
        consumer(graph)
    except Exception:
        return True
    return False


def run_fixture(module: Any) -> dict[str, Any]:
    """Exercise acceptance and fail-closed cases against a fictional graph."""

    checks: dict[str, bool] = {}
    consumer = getattr(module, "consume_saved_raw_g5", None)
    checks["consumer_api"] = callable(consumer)
    selected: dict[str, Any] = {}

    if checks["consumer_api"]:
        graph = _synthetic_graph()
        try:
            selected = consumer(graph)
            checks["accepted"] = selected.get("passed") is True
            checks["schema"] = selected.get("schema") == getattr(
                module, "SCHEMA", "n1_saved_raw_g5_nodes_v1_result"
            )
            checks["core_and_first"] = set(selected.get("g5", {})) == {"core", "first"}
            checks["core_Z0_node"] = (
                selected["g5"]["core"]["ordinary_Z"]["Z0"]["node_id"] == 2
            )
            checks["core_Z1_node"] = (
                selected["g5"]["core"]["ordinary_Z"]["Z1"]["node_id"] == 4
            )
            checks["first_Z0_node"] = (
                selected["g5"]["first"]["ordinary_Z"]["Z0"]["node_id"] == 3
            )
            checks["first_Z1_node"] = (
                selected["g5"]["first"]["ordinary_Z"]["Z1"]["node_id"] == 7
            )
            checks["complete_C3_ids_retained"] = (
                selected["g5"]["core"]["saved_C3_node_ids"] == [2, 4, 5, 6]
                and selected["g5"]["first"]["saved_C3_node_ids"] == [3, 7, 8, 9]
            )
            checks["source_family_preserved"] = (
                selected.get("source_family") == graph["source_family"]
            )
            checks["graph_hash_preserved"] = (
                selected.get("graph_hash") == graph["graph_hash"]
            )
            checks["file_path_preserved"] = (
                selected.get("file_path") == graph["file_path"]
            )
            checks["core_source_record_preserved"] = (
                selected["g5"]["core"]["ordinary_Z"]["Z0"]["raw_record"]
                == graph["exact_graph_nodes"][2]
                and selected["g5"]["core"]["ordinary_Z"]["Z1"]["raw_record"]
                == graph["exact_graph_nodes"][4]
            )
            checks["first_source_record_preserved"] = (
                selected["g5"]["first"]["ordinary_Z"]["Z0"]["raw_record"]
                == graph["exact_graph_nodes"][3]
                and selected["g5"]["first"]["ordinary_Z"]["Z1"]["raw_record"]
                == graph["exact_graph_nodes"][7]
            )
            interpretation = selected.get("interpretation", {})
            checks["symbolic_not_numeric"] = (
                interpretation.get("references_are_symbolic_graph_nodes") is True
                and interpretation.get("point_values_computed") is False
                and interpretation.get("live_owner_verified") is False
            )
            checks["solution_handles_not_used"] = (
                interpretation.get("function_solution_handles_used") is False
            )
            checks["unsupported_claims_false"] = all(
                interpretation.get(name) is False
                for name in (
                    "axis_owner_claimed",
                    "m15_claimed",
                    "cancellation_claimed",
                    "full_Kchi_claimed",
                    "intervals_compared",
                )
            )
        except Exception:
            checks["accepted"] = False
            checks["schema"] = False
            checks["core_and_first"] = False

        missing_Z1 = copy.deepcopy(graph)
        missing_Z1["inner_Volterra_bindings"]["systems"]["core"]["g"]["5"] = [2, 4, 5]
        checks["rejects_missing_Z1"] = _rejected(consumer, missing_Z1)

        out_of_range = copy.deepcopy(graph)
        out_of_range["inner_Volterra_bindings"]["systems"]["first"]["g"]["5"] = [
            3,
            7,
            8,
            99,
        ]
        checks["rejects_out_of_range"] = _rejected(consumer, out_of_range)

        mixed_non_int = copy.deepcopy(graph)
        mixed_non_int["inner_Volterra_bindings"]["systems"]["first"]["g"]["5"] = [
            3,
            7,
            "8",
            9,
        ]
        checks["rejects_mixed_non_int"] = _rejected(consumer, mixed_non_int)
    else:
        for name in (
            "accepted",
            "schema",
            "core_and_first",
            "rejects_missing_Z1",
            "rejects_out_of_range",
            "rejects_mixed_non_int",
        ):
            checks[name] = False

    output: dict[str, Any] = {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "synthetic_only": True,
        "json_safe": True,
        "selected": selected,
    }
    try:
        json.dumps(output, sort_keys=True)
    except (TypeError, ValueError):
        output["json_safe"] = False
        output["passed"] = False
    return output


run = run_fixture
