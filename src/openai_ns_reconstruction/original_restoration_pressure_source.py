"""Source-only restoration/postrestore pressure continuation.

The accepted modified long/reference graph is restored from its serialized
Q4 rows and saved suffixes.  This module appends only the two downstream
restoration charts.  It binds the six physical signed I1 volume functions
and the reduced D/B/S1 source formulae from the pinned Q4 ports, then joins
the two chart increments while retaining every new seam as an uncontracted
source expression.  No accepted owner, builder, checker, producer, norm,
recurrence, quadrature, numerical value, cap, or midpoint is executed.
"""

from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from .original_signed_graph_pairs import GraphPairs


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
MODULE = "src/openai_ns_reconstruction/original_restoration_pressure_source.py"

Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
Q4_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json"
C5_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
C5_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
C5_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json"
C6_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
C6_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
C6_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json"
BRIDGE_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
GRAPH_HELPER = "src/openai_ns_reconstruction/original_signed_graph_pairs.py"

Q4_NODE_COUNT = 183700
BASE_NODE_COUNT = 185704
BASE_DIGEST = "2d6efeb13efe91dd95d4f24fa03eb4626c7685d04e5f29f2bdc6532e101994e3"
PARENT_JOINED = (185700, 185703)
PSTAR_HANDLE = 33572
PSTAR_LOG_ARGUMENT_HANDLE = 5
N_HANDLE = 18

SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}


_CHARTS: dict[str, dict[str, Any]] = {
    "restoration": {
        "native_chart": "restoration",
        "q4_chart": "restoration",
        "c6_chart": "restoration",
        "phase": 116798,
        "global_offset": 116796,
        "coordinate": 83010,
        "radius": 83014,
        "J": 1,
        "U": (177665, 177666, 177667, 177668, 177669, 177670),
        "U_y": (177671, 177672, 177673, 177674, 177675, 177676),
        "m": (177677, 177678, 177679, 177680, 177681, 177682),
        "m_y": (177689, 177690, 177691, 177692, 177693, 177694),
        "Q": (177760, 177764, 177770, 177776, 177782),
        "Q_y": (177848, 177852, 177858, 177864, 177870),
        "finite_RQy": (178007, 178008, 178009, 178010),
        "local_S1": (178037, 178038, 178039, 178040),
        "cumulative_S1": (178041, 178042, 178043, 178044),
        "incoming_S1": (171853, 171854, 171855, 171856),
        "upper_RQy": (178093, 178094, 178095, 178096),
        "incoming_m": (116764, 116765, 116766, 116767, 116768, 116769),
        "history_m": (119228, 119229, 119230, 119231, 119232, 119233),
        "numeric_unit": "Wrestoration",
        "graph_measure": "R^2*1*dq",
    },
    "postrestore": {
        "native_chart": "postrestore",
        "q4_chart": "postrestore",
        "c6_chart": "post_restore",
        "phase": 119416,
        "global_offset": 119414,
        "coordinate": 84528,
        "radius": 84533,
        "J": 1,
        "U": (183268, 183269, 183270, 183271, 183272, 183273),
        "U_y": (183274, 183275, 183276, 183277, 183278, 183279),
        "m": (183280, 183281, 183282, 183283, 183284, 183285),
        "m_y": (183292, 183293, 183294, 183295, 183296, 183297),
        "Q": (183363, 183367, 183373, 183379, 183385),
        "Q_y": (183451, 183455, 183461, 183467, 183473),
        "finite_RQy": (183610, 183611, 183612, 183613),
        "local_S1": (183640, 183641, 183642, 183643),
        "cumulative_S1": (183644, 183645, 183646, 183647),
        "incoming_S1": (178077, 178078, 178079, 178080),
        "upper_RQy": (183696, 183697, 183698, 183699),
        "incoming_m": (119382, 119383, 119384, 119385, 119386, 119387),
        "history_m": (121602, 121603, 121604, 121605, 121606, 121607),
        "numeric_unit": "Wpostrestore",
        "graph_measure": "R^2*1*dq",
    },
}

_I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")


def _require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def _gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"gzip JSON object required: {path}")
    return value


def _row(nodes: Sequence[Mapping[str, Any]], handle: int) -> Mapping[str, Any]:
    _require(type(handle) is int and 0 <= handle < len(nodes), f"bad graph handle {handle}")
    return nodes[handle]


def _pair(value: Sequence[int], label: str) -> tuple[int, int]:
    _require(
        isinstance(value, (list, tuple))
        and len(value) == 2
        and all(type(item) is int for item in value),
        f"{label} must be an ordinary Z 0 / 1 pair",
    )
    return int(value[0]), int(value[1])


def _append_saved(graph: GraphPairs, row: Mapping[str, Any]) -> int:
    _require(isinstance(row, Mapping) and isinstance(row.get("operation"), str), "serialized graph row malformed")
    return graph.node(
        row["operation"],
        **{key: value for key, value in row.items() if key != "operation"},
    )


def _function_source(path: Path, name: str) -> str:
    with path.open("r", encoding="utf-8", newline="") as stream:
        raw = stream.read()
    tree = ast.parse(raw, filename=str(path))
    candidates = [
        ast.get_source_segment(raw, node) or ""
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name
    ]
    _require(candidates, f"source function missing: {path.name}:{name}")
    return max(candidates, key=len)


def _source_guard(root: Path, path_text: str, functions: Mapping[str, Sequence[str]]) -> dict[str, Any]:
    path = root / path_text
    entries: dict[str, Any] = {}
    for name, snippets in functions.items():
        source = _function_source(path, name)
        compact = "".join(source.split())
        tokens = ["".join(token.split()) for token in snippets]
        matched = all(token in compact for token in tokens)
        _require(matched, f"source guard failed: {path.name}:{name}")
        entries[name] = {
            "present": True,
            "tokens": list(snippets),
            "matched": True,
            "source_sha256": hashlib.sha256(source.encode("utf-8")).hexdigest(),
        }
    return {
        "path": path_text,
        "functions": entries,
        "source_guard_tuple": tuple((name, tuple(snippets)) for name, snippets in functions.items()),
        "passed": True,
        "source_executed": False,
    }


def _format_guard(root: Path) -> dict[str, Any]:
    raw = (root / MODULE).read_bytes()
    _require(b"\r" not in raw, "new source must use LF line endings")
    _require(raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "new source must have one final LF")
    tree = ast.parse(raw.decode("utf-8"), filename=MODULE)
    return {
        "passed": True,
        "LF": True,
        "exactly_one_EOF_LF": True,
        "AST_dump_sha256": _digest(ast.dump(tree, include_attributes=False)),
        "source_executed": False,
    }


def _source_pins(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    packet = parent.get("actual_long_reference_pressure_source")
    _require(isinstance(packet, Mapping), "accepted long/reference source packet missing")
    inherited = packet.get("source_pins")
    _require(isinstance(inherited, Mapping), "accepted source pins missing")
    paths = tuple(sorted(inherited))
    result: dict[str, Any] = {}
    for path in paths:
        item = inherited.get(path)
        _require(isinstance(item, Mapping), f"missing accepted source pin: {path}")
        expected = item.get("expected_sha256")
        _require(isinstance(expected, str) and item.get("passed") is True, f"malformed source pin: {path}")
        actual = _sha(root / path)
        _require(actual == expected, f"changed source pin: {path}")
        result[path] = {
            "path": path,
            "expected_sha256": expected,
            "actual_sha256": actual,
            "passed": True,
            "source_executed": False,
        }
        if item.get("self_pin") is True:
            result[path]["self_pin"] = True
    parent_actual = _sha(root / PARENT)
    _require(parent_actual == PARENT_SHA256, "changed accepted parent source pin")
    result[PARENT] = {
        "path": PARENT,
        "expected_sha256": PARENT_SHA256,
        "actual_sha256": parent_actual,
        "passed": True,
        "source_executed": False,
    }
    result[MODULE] = {
        "path": MODULE,
        "expected_sha256": _sha(root / MODULE),
        "actual_sha256": _sha(root / MODULE),
        "passed": True,
        "self_pin": True,
        "source_executed": False,
    }
    return result


def _validate_parent(parent: Mapping[str, Any]) -> dict[str, Any]:
    checks = parent.get("checks")
    _require(isinstance(checks, Mapping) and checks.get("passed") is True,
             "accepted parent checks are not passed")
    packet = checks.get("packet_invariants")
    _require(isinstance(packet, Mapping)
             and packet.get("accepted_math_payload_preserved") is True
             and packet.get("all_standalone_caps_retained") is True
             and packet.get("external_reference_RQy_retained") is True,
             "accepted payload/cap/reference receipt changed")
    _require(parent.get("source_family") == SOURCE_FAMILY, "accepted source family changed")
    reserved = parent.get("actual_reserved_point_m15")
    _require(isinstance(reserved, Mapping) and reserved.get("passed") is True,
             "accepted reserved m15 packet missing")
    _require(reserved.get("old_ledger_or_error_caps_changed") is False,
             "accepted old error/cap ledger changed")
    _require(reserved.get("old_six_signed_contributions_retained") is True,
             "accepted six signed terms were not retained")
    _require(reserved.get("reserved_response_full_support_only") is True,
             "accepted reserved response support changed")
    _require(reserved.get("independent_target_consumed_once") is True,
             "accepted independent target consumption changed")
    _require(reserved.get("independent_pulse_target_enclosure_admitted") is True
             and reserved.get("reserved_response_enclosure_admitted") is True,
             "accepted target/reserved enclosure receipt changed")
    _require(reserved.get("lower_preRp_saved_errors_preserved") is True,
             "accepted pre-Rp saved errors changed")
    pressure = parent.get("actual_long_reference_pressure_source")
    correlated = parent.get("actual_long_reference_correlated_pressure")
    _require(isinstance(pressure, Mapping) and pressure.get("passed") is True, "accepted pressure source missing")
    _require(isinstance(correlated, Mapping) and correlated.get("passed") is True, "accepted correlated source missing")
    _require(pressure.get("source_family") == SOURCE_FAMILY, "pressure source family changed")
    _require(correlated.get("source_family") == SOURCE_FAMILY, "correlated source family changed")
    _require(pressure.get("source_execution") is False and correlated.get("source_executed") is False,
             "accepted source execution receipt changed")
    _require(pressure.get("accepted_ancestors_replayed") is False
             and correlated.get("accepted_ancestors_replayed") is False,
             "accepted ancestor replay receipt changed")
    _require(pressure.get("joined_source_entry") == [185300, 185303], "accepted joined entry changed")
    _require(correlated.get("joined_contracted") == [185700, 185703], "accepted contracted entry changed")
    _require(correlated.get("joined_uncontracted") == [185688, 185689], "accepted raw joined source changed")
    patch = correlated.get("graph_patch")
    _require(isinstance(patch, Mapping), "accepted correlation graph patch missing")
    _require(
        patch.get("base_node_count") == 185690
        and patch.get("base_canonical_nodes_sha256") == "8422242e6c2bce6750323cf2aafd2a5a4bce2bf406f033ddd6a18ad61ae33d28"
        and patch.get("appended_node_count") == 14
        and patch.get("full_node_count") == BASE_NODE_COUNT
        and patch.get("full_canonical_nodes_sha256") == BASE_DIGEST
        and patch.get("base_prefix_unchanged") is True,
        "accepted 185704 graph prefix changed",
    )
    limits = correlated.get("limits", {})
    _require(limits.get("standalone_caps_retained") is True and limits.get("all_old_data_retained") is True,
             "accepted cap/data retention receipt missing")
    _require(pressure.get("limits", {}).get("old_212_terms_and_errors_unchanged") is True,
             "accepted old 212/error receipt missing")
    return {
        "passed": True,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": copy.deepcopy(SOURCE_FAMILY),
        "accepted_prefix_node_count": BASE_NODE_COUNT,
        "accepted_prefix_digest": BASE_DIGEST,
        "parent_joined_contracted": list(PARENT_JOINED),
        "old_212_terms_and_errors_unchanged": True,
        "standalone_caps_retained": True,
        "all_old_data_retained": True,
        "accepted_ancestors_replayed": False,
        "older_germ_unproved_status_retained": True,
        "older_germ_used_as_proof": False,
        "accepted_numeric_computations_replayed": False,
    }


def _validate_q4(root: Path, q4: Mapping[str, Any]) -> dict[str, Any]:
    _require(q4.get("source_family") == SOURCE_FAMILY, "Q4 source family changed")
    nodes = q4.get("exact_graph_nodes")
    _require(isinstance(nodes, list) and len(nodes) == Q4_NODE_COUNT, "Q4 graph node count changed")
    bindings = q4.get("modified_Q4_pressure_C3_bindings")
    _require(isinstance(bindings, Mapping), "Q4 bindings missing")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping), "Q4 chart/endpoints missing")
    receipts: dict[str, Any] = {}
    for name, cfg in _CHARTS.items():
        chart = charts.get(name)
        endpoint = endpoints.get(name)
        _require(isinstance(chart, Mapping) and isinstance(endpoint, Mapping), f"Q4 {name} chart missing")
        physical = chart.get("physical")
        _require(isinstance(physical, Mapping), f"Q4 {name} physical ports missing")
        for field in ("U", "U_y", "m", "m_y", "Q", "Q_y"):
            _require(physical.get(field) == list(cfg[field]), f"Q4 {name} {field} changed")
        _require(chart.get("finite_boundary", {}).get("R_Q_y") == list(cfg["finite_RQy"]),
                 f"Q4 {name} finite R_Q_y changed")
        _require(chart.get("local_integrals", {}).get("S1") == list(cfg["local_S1"]),
                 f"Q4 {name} local S1 changed")
        _require(chart.get("cumulative_integrals", {}).get("S1") == list(cfg["cumulative_S1"]),
                 f"Q4 {name} cumulative S1 changed")
        _require(chart.get("incoming_integrals", {}).get("S1") == list(cfg["incoming_S1"]),
                 f"Q4 {name} incoming S1 changed")
        _require(endpoint.get("boundary", {}).get("R_Q_y") == list(cfg["upper_RQy"]),
                 f"Q4 {name} external R_Q_y changed")
        for order, handle in enumerate(cfg["finite_RQy"][:2]):
            row = _row(nodes, handle)
            _require(row.get("operation") == "product"
                     and row.get("arguments") == [cfg["Q_y"][order], cfg["radius"]],
                     f"Q4 {name} finite R_Q_y source changed")
        for order, handle in enumerate(cfg["local_S1"][:2]):
            row = _row(nodes, handle)
            _require(row.get("operation") == "definite_integral"
                     and row.get("original_chart") == name
                     and row.get("physical_radius_power") == 2
                     and row.get("physical_Jacobian_applied_once") is True
                     and row.get("ordinary_slow_Z_derivative_order") == order,
                     f"Q4 {name} local S1 source changed")
        receipts[name] = {
            "physical_ports": {field: list(cfg[field]) for field in ("U", "U_y", "m", "m_y", "Q", "Q_y")},
            "finite_RQy": list(cfg["finite_RQy"][:2]),
            "local_S1": list(cfg["local_S1"][:2]),
            "cumulative_S1": list(cfg["cumulative_S1"][:2]),
            "incoming_S1": list(cfg["incoming_S1"][:2]),
            "upper_RQy": list(cfg["upper_RQy"][:2]),
            "physical_radius_power": 2,
            "physical_Jacobian_applied_once": True,
        }
    _require(bindings.get("definition", {}).get("physical_U_mean_Uy_C5_and_Q_Qy_C4_Pstar_once") is True,
             "Q4 physical Pstar source premise missing")
    _require(bindings.get("definition", {}).get("reduced_original_S_and_signed_m13_m15_C3_through_Rm") is True,
             "Q4 reduced pressure source premise missing")
    return {
        "passed": True,
        "node_count": len(nodes),
        "common_N": N_HANDLE,
        "Pstar": PSTAR_HANDLE,
        "Pstar_log_argument": PSTAR_LOG_ARGUMENT_HANDLE,
        "charts": receipts,
        "source_executed": False,
    }


def _validate_c5(root: Path, c5: Mapping[str, Any]) -> dict[str, Any]:
    _require(c5.get("source_family") == SOURCE_FAMILY, "C5 source family changed")
    bindings = c5.get("prescribed_loop_C5_bindings")
    _require(isinstance(bindings, Mapping), "C5 bindings missing")
    _require(bindings.get("N") == N_HANDLE, "C5 common N changed")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping), "C5 charts/endpoints missing")
    result: dict[str, Any] = {}
    for name, cfg in _CHARTS.items():
        chart = charts.get(name)
        _require(isinstance(chart, Mapping), f"C5 {name} chart missing")
        _require(chart.get("actual_phase") == cfg["phase"]
                 and chart.get("global_offset") == cfg["global_offset"]
                 and chart.get("coordinate") == cfg["coordinate"]
                 and chart.get("radius") == cfg["radius"]
                 and chart.get("physical_measure") == cfg["J"]
                 and chart.get("domain") == [0, 1], f"C5 {name} geometry changed")
        incoming = chart.get("incoming")
        history = chart.get("history")
        memory = chart.get("memory")
        _require(isinstance(incoming, Mapping) and isinstance(history, Mapping) and isinstance(memory, Mapping),
                 f"C5 {name} memory ports missing")
        _require(incoming.get("m") == list(cfg["incoming_m"]) and history.get("m") == list(cfg["history_m"]),
                 f"C5 {name} m history changed")
        for memory_name in ("m", "h", "k", "e", "p"):
            _require(isinstance(incoming.get(memory_name), list)
                     and len(incoming[memory_name]) == 6
                     and all(type(item) is int for item in incoming[memory_name])
                     and isinstance(history.get(memory_name), list)
                     and len(history[memory_name]) == 6
                     and all(type(item) is int for item in history[memory_name])
                     and type(memory.get(memory_name)) is int,
                     f"C5 {name} {memory_name} memory changed")
        _require(memory.get("p") == 1, f"C5 {name} pressure memory reset")
        result[name] = {
            "actual_phase": cfg["phase"],
            "global_offset": cfg["global_offset"],
            "coordinate": cfg["coordinate"],
            "radius": cfg["radius"],
            "physical_measure": cfg["J"],
            "incoming": copy.deepcopy(dict(incoming)),
            "history": copy.deepcopy(dict(history)),
            "memory": copy.deepcopy(dict(memory)),
            "nonreset_memories": True,
            "phase_Z_independent": True,
        }
    _require(endpoints.get("reference", {}).get("m") == charts["restoration"].get("incoming", {}).get("m"),
             "reference to restoration m history reset")
    _require(endpoints.get("restoration", {}).get("m") == charts["postrestore"].get("incoming", {}).get("m"),
             "restoration to postrestore m history reset")
    definition = bindings.get("definition", {})
    for key in (
        "complete_original_rminus_to_Rm_Duhamel_correction_histories_C5",
        "correction_zero_inlet_only_at_rminus_background_nonzero_memories_unchanged",
        "original_rate_zero_pressure_memory_retained",
    ):
        _require(definition.get(key) is True, f"C5 source premise missing: {key}")
    return {
        "passed": True,
        "common_N": N_HANDLE,
        "charts": result,
        "postrestore_alias": "postrestore",
        "nonreset_incoming_histories": True,
        "source_executed": False,
    }


def _validate_c6(root: Path, c6: Mapping[str, Any]) -> dict[str, Any]:
    _require(c6.get("source_family") == SOURCE_FAMILY, "C6 source family changed")
    bindings = c6.get("prescribed_R100_Rm_C6_bindings")
    _require(isinstance(bindings, Mapping), "C6 bindings missing")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping), "C6 charts/endpoints missing")
    result: dict[str, Any] = {}
    for name, cfg in _CHARTS.items():
        c6_name = cfg["c6_chart"]
        chart = charts.get(c6_name)
        _require(isinstance(chart, Mapping), f"C6 {c6_name} chart missing")
        _require(chart.get("phase") == cfg["coordinate"]
                 and chart.get("radius") == cfg["radius"]
                 and chart.get("physical_measure") == cfg["J"]
                 and chart.get("domain") == [0, 1], f"C6 {c6_name} geometry changed")
        actual = chart.get("actual")
        moments = chart.get("moments")
        _require(isinstance(actual, Mapping) and isinstance(moments, Mapping), f"C6 {c6_name} ports missing")
        for field in ("phi", "U", "phi_y", "U_y"):
            values = actual.get(field)
            _require(isinstance(values, list) and len(values) == 7 and all(type(item) is int for item in values),
                     f"C6 {c6_name} {field} changed")
        for moment_name in ("H", "M", "K", "A", "B", "C"):
            values = moments.get(moment_name)
            _require(isinstance(values, list) and len(values) == 7 and all(type(item) is int for item in values),
                     f"C6 {c6_name} {moment_name} changed")
        result[name] = {
            "chart": c6_name,
            "phase": cfg["coordinate"],
            "radius": cfg["radius"],
            "physical_measure": cfg["J"],
            "actual_rows": {field: list(actual[field]) for field in ("phi", "U", "phi_y", "U_y")},
            "moment_row_lengths": {field: len(moments[field]) for field in ("H", "M", "K", "A", "B", "C")},
        }
    _require(isinstance(endpoints.get("Rz"), Mapping)
             and isinstance(endpoints.get("Rrestore"), Mapping)
             and isinstance(endpoints.get("Rm"), Mapping), "C6 Rz/Rrestore/Rm endpoint chain missing")
    _require(endpoints["Rrestore"].get("radius") == 84352
             and endpoints["Rm"].get("radius") == 85743, "C6 restoration endpoint radii changed")
    _require(bindings.get("definition", {}).get("all_six_nonzero_incoming_and_every_chart_memory_retained") is True,
             "C6 nonreset memory premise missing")
    return {
        "passed": True,
        "charts": result,
        "endpoint_chain": ["Rz", "Rrestore", "Rm"],
        "postrestore_alias": "post_restore",
        "source_executed": False,
    }


def _rehydrate_prefix(parent: Mapping[str, Any], q4: Mapping[str, Any]) -> tuple[GraphPairs, dict[str, Any]]:
    graph = GraphPairs(q4["exact_graph_nodes"])
    pieces = (
        ("actual_micro_macro_correlated_m15", "merged_namespace", Q4_NODE_COUNT, 987),
        ("actual_macro_switch_RQy_seam", "graph_patch", 184687, 16),
        ("actual_macro_first_switch_pressure_source", "graph_patch", 184703, 195),
        ("actual_switch_post_pressure_source", "graph_patch", 184898, 410),
        ("actual_long_reference_pressure_source", "graph_patch", 185308, 382),
        ("actual_long_reference_correlated_pressure", "graph_patch", 185690, 14),
    )
    for key, field, start, length in pieces:
        packet = parent.get(key, {}).get(field)
        _require(isinstance(packet, Mapping), f"accepted graph segment missing: {key}")
        suffix = packet.get("appended_nodes")
        _require(packet.get("base_node_count") == start and isinstance(suffix, list) and len(suffix) == length,
                 f"accepted graph segment changed: {key}")
        _require(len(graph.nodes) == start, f"accepted graph segment start changed: {key}")
        for row in suffix:
            _append_saved(graph, row)
    _require(len(graph.nodes) == BASE_NODE_COUNT and _digest(graph.nodes) == BASE_DIGEST,
             "accepted 185704 graph prefix changed")
    return graph, {
        "base_node_count": BASE_NODE_COUNT,
        "base_canonical_nodes_sha256": BASE_DIGEST,
        "base_prefix_unchanged": True,
        "full_prefix_snapshot": copy.deepcopy(graph.nodes),
        "serialized_rehydration": True,
        "segment_starts": {"q4": 0, "accepted_merged": 183700, "accepted_seam": 184687,
                            "accepted_first": 184703, "accepted_post": 184898,
                            "accepted_long_reference": 185308, "accepted_correlation": 185690,
                            "new_suffix": BASE_NODE_COUNT},
        "suffix_counts": {"accepted_merged": 987, "accepted_seam": 16, "accepted_first": 195,
                           "accepted_post": 410, "accepted_long_reference": 382,
                           "accepted_correlation": 14},
    }


def _build_six_source(graph: GraphPairs, cfg: Mapping[str, Any]) -> tuple[dict[str, list[int]], dict[str, Any]]:
    q = _pair(cfg["Q"][:2], "Q")
    q_z = _pair(cfg["Q"][1:3], "Q_Z")
    u = _pair(cfg["U"][:2], "U")
    u_z = _pair(cfg["U"][1:3], "U_Z")
    density = {
        "q": q,
        "q_Z": q_z,
        "q_squared": graph.pmul(q, q),
        "U_q_Z": graph.pmul(u, q_z),
        "U_squared": graph.pmul(u, u),
        "U_U_Z": graph.pmul(u, u_z),
    }
    weight = graph.mul(cfg["radius"], cfg["radius"], cfg["J"])
    integrals: dict[str, list[int]] = {}
    for name in _I1_NAMES:
        integrals[name] = list(
            graph.integral(
                graph.pscale(density[name], weight),
                cfg["coordinate"],
                0,
                1,
                f"actual_{cfg['native_chart']}_current_N_{name}",
                metadata={
                    "original_chart": cfg["q4_chart"],
                    "source_layer": "physical_six_I1_density",
                    "ordinary_Z_rows": [0, 1],
                    "physical_radius_power": 2,
                    "physical_Jacobian_applied_once": True,
                    "Pstar_applied_in_physical_fields_once": True,
                    "endpoint_independent_of_Z": True,
                    "physical_y_derivatives_claimed": False,
                    "raw_physical_unit": "R^2*J*dcoordinate",
                    "numerical_value_not_installed": True,
                },
            )
        )
    return integrals, {
        "passed": True,
        "density_pairs": {key: list(value) for key, value in density.items()},
        "weight_handle": weight,
        "weight_formula": f"R_{cfg['native_chart']}^2 * J_{cfg['native_chart']}",
        "weight_units": cfg["graph_measure"],
        "graph_source_raw_physical": True,
        "ordinary_Z_orders": [0, 1],
        "Pstar_applied_once": True,
        "R_not_a_mean_scale": True,
        "integral_handles": integrals,
        "numerical_values_selected": False,
    }


def _build_reduced_D(graph: GraphPairs, integrals: Mapping[str, Sequence[int]]) -> tuple[list[int], dict[str, Any]]:
    z = _pair((20, 1), "Z")
    d = _pair((184430, 184429), "d")
    h = _pair((184432, 184431), "h")
    L = _pair((184439, 184438), "L")
    L2 = _pair((184440, 184443), "L^2")
    one, _ = graph.fixed(1)
    half = 28
    two = 48
    inv_L = graph.pdiv((one, graph.constant(0)), L, "positive physical L=1-delta*Z^2")
    hZ_over_2L = graph.pscale(
        graph.pdiv(graph.pmul(h, z), L, "positive physical L=1-delta*Z^2"), half
    )
    minus_half = graph.pneg(graph.fixed(half))
    d_over_L = graph.pdiv(d, L, "positive physical L=1-delta*Z^2")
    two_hZ2_over_L2 = graph.pdiv(
        graph.pscale(graph.pmul(h, graph.pmul(z, z)), two),
        L2,
        "positive physical L=1-delta*Z^2",
    )
    two_Zd_over_L2 = graph.pdiv(
        graph.pscale(graph.pmul(z, d), two),
        L2,
        "positive physical L=1-delta*Z^2",
    )
    terms = {
        "minus_inv_L": graph.pmul(graph.pneg(inv_L), _pair(integrals["q"], "q")),
        "hZ_over_2L": graph.pmul(hZ_over_2L, _pair(integrals["q_Z"], "q_Z")),
        "minus_half": graph.pmul(minus_half, _pair(integrals["q_squared"], "q_squared")),
        "d_over_L": graph.pmul(d_over_L, _pair(integrals["U_q_Z"], "U_q_Z")),
        "two_hZ2_over_L2": graph.pmul(two_hZ2_over_L2, _pair(integrals["U_squared"], "U_squared")),
        "two_Zd_over_L2": graph.pmul(two_Zd_over_L2, _pair(integrals["U_U_Z"], "U_U_Z")),
    }
    rows = list(graph.padd(*terms.values()))
    return rows, {
        "passed": True,
        "rows": rows,
        "terms": {key: list(value) for key, value in terms.items()},
        "formula": "D=-q/L+h*Z*q_Z/(2L)-q^2/2+d*U*q_Z/L+2*h*Z^2*U^2/L^2+2*Z*d*U*U_Z/L^2",
        "graph_source_unit": "raw physical R^2*J*dcoordinate",
        "normalization_not_in_graph": True,
        "D_is_not_full_S1": True,
        "correct_two_h_Z2_term": True,
        "numerical_values_selected": False,
    }


def _build_reduced_B(graph: GraphPairs, cfg: Mapping[str, Any]) -> tuple[list[int], dict[str, Any]]:
    z = _pair((20, 1), "Z")
    u = _pair(cfg["U"][:2], "U")
    q = _pair(cfg["Q"][:2], "Q")
    L = _pair((184439, 184438), "L")
    L2 = _pair((184440, 184443), "L^2")
    radius_squared = graph.mul(cfg["radius"], cfg["radius"])
    first = graph.pscale(graph.pdiv(q, L, "positive physical L=1-delta*Z^2"), radius_squared)
    second = graph.pscale(graph.pscale(graph.pmul(q, q), radius_squared), 28)
    third = graph.pscale(
        graph.pmul(graph.pmul(z, z), graph.fixed(radius_squared), graph.pmul(u, u)), 48
    )
    third = graph.pneg(graph.pdiv(third, L2, "positive physical L=1-delta*Z^2"))
    boundary = graph.padd(first, second, third)
    upper = graph.at(boundary, cfg["coordinate"], 1)
    lower = graph.at(boundary, cfg["coordinate"], 0)
    rows = list(graph.padd(upper, graph.pneg(lower)))
    return rows, {
        "passed": True,
        "rows": rows,
        "boundary_integrand": list(boundary),
        "upper": list(upper),
        "lower": list(lower),
        "formula": "B=[R^2*Q/L+R^2*Q^2/2-2*Z^2*R^2*U^2/L^2]_{t=0}^{t=1}",
        "physical_radius_power": 2,
        "physical_Jacobian_applied": False,
        "endpoint_measure_has_no_J": True,
        "graph_source_unit": "raw physical R^2*Q endpoint",
        "numeric_endpoint_not_evaluated": True,
        "ordinary_Z_orders": [0, 1],
        "correct_minus_two_Z2_R2_U2_term": True,
        "numerical_values_selected": False,
    }


def _linear_coefficients(nodes: Sequence[Mapping[str, Any]], root: int) -> dict[int, int]:
    out: dict[int, int] = {}
    pending = [(root, 1)]
    while pending:
        handle, sign = pending.pop()
        row = _row(nodes, handle)
        operation = row.get("operation")
        if operation == "sum":
            pending.extend((child, sign) for child in row["arguments"])
        elif operation == "negative":
            pending.append((row["argument"], -sign))
        elif operation == "exact_rational" and int(row.get("numerator", 1)) == 0:
            continue
        else:
            out[handle] = out.get(handle, 0) + sign
            if out[handle] == 0:
                del out[handle]
    return out


def _build_chart(graph: GraphPairs, cfg: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    integrals, six = _build_six_source(graph, cfg)
    b_rows, b_source = _build_reduced_B(graph, cfg)
    d_rows, d_source = _build_reduced_D(graph, integrals)
    reduced = list(graph.padd(_pair(b_rows, "B"), _pair(d_rows, "D")))
    lower = graph.at(_pair(cfg["finite_RQy"][:2], f"{cfg['native_chart']} finite R_Q_y"),
                     cfg["coordinate"], graph.constant(0))
    output = {
        "passed": True,
        "native_chart": cfg["native_chart"],
        "q4_chart": cfg["q4_chart"],
        "c6_chart": cfg["c6_chart"],
        "c5_postrestore_alias": "postrestore" if cfg["native_chart"] == "postrestore" else None,
        "common_N": N_HANDLE,
        "N_materialized": False,
        "phase": cfg["phase"],
        "global_offset": cfg["global_offset"],
        "phase_Z_independent": True,
        "coordinate": cfg["coordinate"],
        "domain": [0, 1],
        "radius": cfg["radius"],
        "native_J": cfg["J"],
        "physical_ports": {key: list(cfg[key]) for key in ("U", "U_y", "m", "m_y", "Q", "Q_y")},
        "incoming_m": list(cfg["incoming_m"]),
        "history_m": list(cfg["history_m"]),
        "native_memories": ["m", "h", "k", "e", "p"],
        "native_history_moments": ["H", "M", "K", "A", "B", "C"],
        "local_S1": list(cfg["local_S1"][:2]),
        "cumulative_S1": list(cfg["cumulative_S1"][:2]),
        "incoming_S1": list(cfg["incoming_S1"][:2]),
        "finite_boundary_RQy": list(cfg["finite_RQy"][:2]),
        "lower_wrapper_t0_or_q0": list(lower),
        "lower_RQy_wrappers": list(lower),
        "lower_RQy": list(lower),
        "upper_RQy_standalone": list(cfg["upper_RQy"][:2]),
        "upper_RQy_traces": list(cfg["upper_RQy"][:2]),
        "upper_RQy": list(cfg["upper_RQy"][:2]),
        "six_I1_source": six,
        "integrals": integrals,
        "B_source": b_source,
        "D_source": d_source,
        "S1_source": {"rows": reduced, "formula": "S1=B+D", "D_is_not_full_S1": True},
        "raw_graph_unit": cfg["graph_measure"],
        "physical_R_squared_applied_once": True,
        "physical_Jacobian_applied_once": True,
        "Pstar_applied_once": True,
        "physical_R_not_mean_scale": True,
        "source_values_selected": False,
        "interval_or_cap_values_used": False,
    }
    return output, {"local_S1": list(cfg["local_S1"][:2]), "lower": list(lower), "upper": list(cfg["upper_RQy"][:2])}


def _build_join(graph: GraphPairs, charts: Mapping[str, Mapping[str, Any]], ports: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    half = {name: [graph.mul(int(value), 28) for value in ports[name]["local_S1"]] for name in _CHARTS}
    increments = {
        name: [
            graph.add(half[name][order], ports[name]["lower"][order], graph.neg(ports[name]["upper"][order]))
            for order in range(2)
        ]
        for name in _CHARTS
    }
    after_restoration = [graph.add(PARENT_JOINED[order], increments["restoration"][order]) for order in range(2)]
    after_postrestore = [graph.add(after_restoration[order], increments["postrestore"][order]) for order in range(2)]
    parent_coefficients = [_linear_coefficients(graph.nodes, root) for root in PARENT_JOINED]
    final_coefficients = [_linear_coefficients(graph.nodes, root) for root in after_postrestore]
    external_reference = [171869, 171870]
    external_rm = list(_CHARTS["postrestore"]["upper_RQy"][:2])
    _require(all(parent_coefficients[order].get(external_reference[order]) == -1 for order in range(2)),
             "accepted external reference RQy coefficient changed")
    _require(all(final_coefficients[order].get(external_reference[order]) == -1 for order in range(2)),
             "external reference RQy coefficient was not retained")
    _require(all(final_coefficients[order].get(external_rm[order]) == -1 for order in range(2)),
             "external Rm RQy coefficient was not retained")
    return {
        "passed": True,
        "ordinary_Z_orders": [0, 1],
        "parent_joined_contracted": list(PARENT_JOINED),
        "restoration_local_S1": list(charts["restoration"]["local_S1"][:2]),
        "postrestore_local_S1": list(charts["postrestore"]["local_S1"][:2]),
        "restoration_half_S1": half["restoration"],
        "postrestore_half_S1": half["postrestore"],
        "restoration_lower_RQy_t0": list(ports["restoration"]["lower"]),
        "postrestore_lower_RQy_t0": list(ports["postrestore"]["lower"]),
        "restoration_upper_RQy": list(ports["restoration"]["upper"]),
        "external_Rm_RQy": external_rm,
        "external_reference_RQy_retained": external_reference,
        "increments_uncontracted": increments,
        "joined_after_restoration_uncontracted": after_restoration,
        "joined_after_postrestore_uncontracted": after_postrestore,
        "joined_uncontracted": after_postrestore,
        "joined_contracted": None,
        "contracted_before_linearization": None,
        "seam_differences_uncontracted": True,
        "contraction_allowed": False,
        "contraction_status": "BLOCKED_UNPROVEN_RESTORATION_POSTRESTORE_SOURCE_CONGRUENCE",
        "linearization_receipt": {
            "passed": False,
            "contraction_allowed": False,
            "uncontracted_source_preserved": True,
            "standalone_caps_preserved": True,
            "reason": "actual restoration/postrestore seam source proof was not supplied",
        },
        "formula_uncontracted": (
            "parent_joined + restoration_local_S1/2 + restoration_R_Q_y(0) - restoration_R_Q_y(1) "
            "+ postrestore_local_S1/2 + postrestore_R_Q_y(0) - Rm_R_Q_y"
        ),
        "formula_contracted": None,
        "local_S1_added_once": True,
        "cumulative_S1_not_used": True,
        "external_Rm_RQy_retained": True,
        "external_reference_RQy_coefficient": [-1, -1],
        "external_Rm_RQy_coefficient": [-1, -1],
        "signed_relation_not_global_closure": True,
        "source_values_selected": False,
        "interval_or_cap_values_used": False,
    }


def restored_nodes(root: str | Path, parent: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Rehydrate only the accepted 185704-node serialized graph prefix.

    This helper reads the pinned Q4 archive and saved graph suffix rows.  It
    never calls an accepted owner or source builder and returns a detached
    serialization copy for a later append-only source extension.
    """

    root_path = Path(root)
    _require(_sha(root_path / PARENT) == PARENT_SHA256, "accepted modified-seams parent bytes changed")
    saved_parent = _json(root_path / PARENT)
    _require(_digest(saved_parent) == _digest(parent), "exact accepted modified-seams parent required")
    _validate_parent(parent)
    q4 = _gzip_json(root_path / Q4_ARCHIVE)
    _validate_q4(root_path, q4)
    graph, _ = _rehydrate_prefix(parent, q4)
    return copy.deepcopy(graph.nodes)


def build(root: str | Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    """Append the source-defined restoration and postrestore pressure terms."""

    root_path = Path(root)
    _require(isinstance(parent, Mapping), "parent report mapping required")
    _require(_sha(root_path / PARENT) == PARENT_SHA256, "accepted modified-seams parent bytes changed")
    saved_parent = _json(root_path / PARENT)
    _require(_digest(saved_parent) == _digest(parent), "exact accepted modified-seams parent required")
    parent_scope = _validate_parent(parent)
    pins = _source_pins(root_path, parent)
    q4 = _gzip_json(root_path / Q4_ARCHIVE)
    c5 = _gzip_json(root_path / C5_ARCHIVE)
    c6 = _gzip_json(root_path / C6_ARCHIVE)
    q4_check = _validate_q4(root_path, q4)
    c5_check = _validate_c5(root_path, c5)
    c6_check = _validate_c6(root_path, c6)
    guards = {
        "self_format": _format_guard(root_path),
        "Q4": _source_guard(root_path, Q4_SOURCE, {
            "physical_pressure_fields": (
                "U=a.scale(fields['V'],Pstar)",
                "Uy=a.scale(radial['V_y'],Pstar)",
                "m=a.scale(fields['m'],Pstar)",
                "my=a.sub(U,m)",
                "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
            ),
            "reduced_pressure_S": (
                "linear=div(a.add(Q,a.scale(a.mul(a.sub(one,de),a.mul(z,Qz)),'1/2'),Qy))",
                "quadratic=a.mul(Q,a.add(a.scale(Q,'1/2'),Qy))",
            ),
        }),
        "C5": _source_guard(root_path, C5_SOURCE, {
            "transport": (
                "width=g.mul(jacobian,g.sub(coordinate,lower))",
                "contribution=a.integral(body,t,lower,coordinate,actual_finite_N_signed_correction=True",
                "physical_Jacobian_applied_once=True",
                "predecessor_memory_never_reset=True",
            ),
            "loop_functions": (
                "inverse=g.node('original_monotone_phase_inverse'",
                "Aphi=lazy(",
                "Bphi=lazy(",
            ),
        }),
        "bridge": _source_guard(root_path, BRIDGE_SOURCE, {
            "generic_fields": (
                "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)",
                "radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y'])),V_y=a.scale(actual['U_y'],inverse))",
            ),
            "densities": (
                "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2)",
                "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi))",
            ),
            "histories": (
                "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase",
                "original_incoming_memory_retained=True",
                "source_measure_applied_once=True",
            ),
        }),
        "C6": _source_guard(root_path, C6_SOURCE, {
            "__init__": (
                "t=g.symbol('true_prescribed_restoration_t')",
                "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])",
                "p=g.symbol('true_prescribed_post_restore_fraction')",
                "self.add_chart('post_restore',p,R,g.one,actual,self.endpoints['Rrestore']['moments'])",
                "self.endpoints['Rm']=self.end('post_restore',g.one)",
            ),
            "add_chart": (
                "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
                "fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)",
            ),
        }),
    }
    guards["source_execution"] = False
    guards["source_values_selected"] = False
    graph, prefix = _rehydrate_prefix(parent, q4)
    chart_output: dict[str, Any] = {}
    ports: dict[str, dict[str, Any]] = {}
    for name, cfg in _CHARTS.items():
        chart_output[name], ports[name] = _build_chart(graph, cfg)
    join = _build_join(graph, chart_output, ports)
    for name in _CHARTS:
        chart_output[name]["local_half_S1"] = list(join[f"{name}_half_S1"])
    _require(graph.nodes[:BASE_NODE_COUNT] == prefix["full_prefix_snapshot"], "accepted 185704-node prefix changed")
    suffix = copy.deepcopy(graph.nodes[BASE_NODE_COUNT:])
    _require(suffix, "restoration source did not append graph rows")
    graph_patch = {
        "base_archive": Q4_ARCHIVE,
        "base_node_count": BASE_NODE_COUNT,
        "base_prefix_unchanged": True,
        "base_canonical_nodes_sha256": BASE_DIGEST,
        "base_prefix_canonical_nodes_sha256": BASE_DIGEST,
        "new_node_start": BASE_NODE_COUNT,
        "appended_nodes": suffix,
        "appended_node_count": len(suffix),
        "full_node_count": len(graph.nodes),
        "full_prefix_digest": BASE_DIGEST,
        "parent_full_canonical_nodes_sha256": BASE_DIGEST,
        "full_canonical_nodes_sha256": _digest(graph.nodes),
        "segment_starts": prefix["segment_starts"],
        "prefix_suffix_counts": prefix["suffix_counts"],
        "endpoint_handles": {
            "parent_joined_contracted": list(PARENT_JOINED),
            "parent_joined": list(PARENT_JOINED),
            "parent_external_reference_RQy": [171869, 171870],
            "restoration_lower_RQy": list(ports["restoration"]["lower"]),
            "restoration_upper_RQy": list(ports["restoration"]["upper"]),
            "postrestore_lower_RQy": list(ports["postrestore"]["lower"]),
            "external_Rm_RQy": list(ports["postrestore"]["upper"]),
            "restoration_local_S1": list(chart_output["restoration"]["local_S1"][:2]),
            "postrestore_local_S1": list(chart_output["postrestore"]["local_S1"][:2]),
            "joined_after_restoration_uncontracted": join["joined_after_restoration_uncontracted"],
            "joined_after_postrestore_uncontracted": join["joined_after_postrestore_uncontracted"],
            "joined_uncontracted": join["joined_uncontracted"],
            "new_raw_joined": join["joined_uncontracted"],
            "joined_contracted": None,
        },
    }
    chart_handles = {
        "parent_joined": list(PARENT_JOINED),
        "parent_joined_contracted": list(PARENT_JOINED),
        "restoration": {
            "lower_RQy_wrappers": list(ports["restoration"]["lower"]),
            "upper_RQy_traces": list(ports["restoration"]["upper"]),
            "local_half_S1": list(join["restoration_half_S1"]),
            "local_S1": list(chart_output["restoration"]["local_S1"][:2]),
        },
        "postrestore": {
            "lower_RQy_wrappers": list(ports["postrestore"]["lower"]),
            "upper_RQy_traces": list(ports["postrestore"]["upper"]),
            "local_half_S1": list(join["postrestore_half_S1"]),
            "local_S1": list(chart_output["postrestore"]["local_S1"][:2]),
        },
        "new_raw_joined": list(join["joined_uncontracted"]),
        "new_raw_joined_uncontracted": list(join["joined_uncontracted"]),
        "new_raw_joined_after_restoration": list(join["joined_after_restoration_uncontracted"]),
        "new_raw_joined_after_postrestore": list(join["joined_after_postrestore_uncontracted"]),
    }
    graph_patch["stable_chart_handles"] = copy.deepcopy(chart_handles)
    preservation = {
        "accepted_185704_prefix_unchanged": True,
        "old_212_terms_and_errors_unchanged": True,
        "old_six_signed_terms_retained": parent.get("actual_reserved_point_m15", {}).get("old_six_signed_contributions_retained") is True,
        "reserved_n1_own_swirl_retained": parent.get("actual_reserved_point_m15", {}).get("reserved_response_full_support_only") is True,
        "negative_full_pulse_m15_target_retained": parent.get("actual_reserved_point_m15", {}).get("independent_target_consumed_once") is True,
        "lower_preRp_saved_errors_preserved": parent.get("actual_reserved_point_m15", {}).get("lower_preRp_saved_errors_preserved") is True,
        "old_ledger_or_error_caps_unchanged": parent.get("actual_reserved_point_m15", {}).get("old_ledger_or_error_caps_changed") is False,
        "all_standalone_caps_retained": True,
        "all_error_caps_retained": True,
        "target_handles_retained": True,
        "new_seam_atoms_uncontracted": True,
    }
    source_contract = {
        "passed": True,
        "restoration_postrestore_source_only": True,
        "ordinary_Z_orders": [0, 1],
        "common_N_handle": N_HANDLE,
        "N_not_materialized": True,
        "phase_Z_independent": True,
        "physical_Pstar_applied_once": True,
        "physical_R_squared_applied_once": True,
        "physical_Jacobian_applied_once": True,
        "native_J": {name: cfg["J"] for name, cfg in _CHARTS.items()},
        "nonreset_memories": True,
        "all_six_signed_I1_sources_defined": True,
        "D_B_S1_sources_defined": True,
        "seam_differences_contracted": False,
        "external_Rm_RQy_retained": True,
        "global_scale_exp_materialized": False,
        "accepted_ancestors_replayed": False,
        "source_execution": False,
        "source_executed": False,
        "complete_pressure_numeric": False,
        "global_five_moment_closure": False,
        "full_goal_complete": False,
    }
    limits = {
        "restoration_source_defined": True,
        "postrestore_source_defined": True,
        "six_physical_signed_I1_defined": True,
        "joined_pressure_source_defined": True,
        "joined_pressure_numeric_enclosure": False,
        "source_formula_values_selected": False,
        "interval_or_cap_values_used": False,
        "new_seams_proved": False,
        "new_seams_contracted": False,
        "external_Rm_RQy_retained": True,
        "all_caps_retained": True,
        "all_old_data_retained": True,
        "old_212_terms_and_errors_unchanged": True,
        "global_closure": False,
        "complete_pressure_global_closure": False,
        "accepted_ancestors_replayed": False,
        "N_materialized": False,
        "mathematical_completion": False,
    }
    result = {
        "passed": True,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": copy.deepcopy(SOURCE_FAMILY),
        "module_sha256": _sha(root_path / MODULE),
        "source_pins": pins,
        "input_pins": copy.deepcopy(pins),
        "source_guards": guards,
        "source_guard_tuples": {
            key: value.get("source_guard_tuple")
            for key, value in guards.items()
            if isinstance(value, Mapping) and "source_guard_tuple" in value
        },
        "parent_scope": parent_scope,
        "c6_cache": c6_check,
        "c5_cache": c5_check,
        "q4_cache": q4_check,
        "charts": chart_output,
        "joined_source_entry": list(PARENT_JOINED),
        "joined_pressure_source": join,
        "stable_chart_handles": chart_handles,
        "graph_patch": graph_patch,
        "preservation": preservation,
        "source_contract": source_contract,
        "contract": copy.deepcopy(source_contract),
        "limits": limits,
        "signed_relation_not_global_closure": True,
        "accepted_ancestors_replayed": False,
        "source_execution": False,
        "source_executed": False,
        "complete_pressure_numeric": False,
        "global_five_moment_closure": False,
        "full_goal_complete": False,
        "mathematical_completion": {
            "complete_pressure_numeric": False,
            "global_five_moment_closure": False,
            "exact_heat_energy_stress_flat": False,
            "true_temporal_n_recursion": False,
            "mean_oscillatory_correction": False,
            "independent_corrected_Cartesian_NS": False,
            "full_goal_complete": False,
        },
    }
    json.dumps(result, sort_keys=True, ensure_ascii=True, allow_nan=False)
    return result


__all__ = ["build", "restored_nodes"]
