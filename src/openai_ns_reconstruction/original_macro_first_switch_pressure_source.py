"""Source-only finite-N first-switch pressure extension.

The fresh first-switch finite-N numeric packet supplies signed six-I1 and
reduced-D enclosures.  This adapter binds those rows to the full Q4
first-switch source and appends an ordinary-Z source graph for the local
first-switch volume.  It then extends the already joined micro/macro pressure
source by one first-switch segment.  No archived owner, producer, checker,
recurrence, norm, quadrature, or numerical builder is imported or executed
here.
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
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_macro_finite_N_I1.json"
PARENT_SHA256 = "066c6eaf48b49f836379b1604725f1b7f33ee9372bb2e4521652ddeaff1d7c3c"
MODULE = "src/openai_ns_reconstruction/original_macro_first_switch_pressure_source.py"

Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
Q4_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json"
LOOP_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
LOOP_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
LOOP_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json"
DENSITY_SOURCE = "src/openai_ns_reconstruction/original_finite_N_density_kernel.py"
GRAPH_HELPER = "src/openai_ns_reconstruction/original_signed_graph_pairs.py"
SEAM_SOURCE = "src/openai_ns_reconstruction/original_macro_switch_RQy_seam.py"

ORDER = 1
I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")

N_HANDLE = 18
MACRO_PHASE_HANDLE = 101909
FIRST_PHASE_HANDLE = 104537
MACRO_J_HANDLE = 66987
FIRST_J_HANDLE = 9
MACRO_RADIUS_HANDLE = 71448
FIRST_RADIUS_HANDLE = 85905
MACRO_COORDINATE_HANDLE = 70265
FIRST_COORDINATE_HANDLE = 73190
PSTAR_HANDLE = 33572
PSTAR_LOG_ARGUMENT_HANDLE = 5

FIRST_U = (148085, 148086, 148087, 148088, 148089, 148090)
FIRST_UY = (148091, 148092, 148093, 148094, 148095, 148096)
FIRST_M = (148097, 148098, 148099, 148100, 148101, 148102)
FIRST_MY = (148109, 148110, 148111, 148112, 148113, 148114)
FIRST_Q = (148180, 148184, 148190, 148196, 148202)
FIRST_QY = (148268, 148272, 148278, 148284, 148290)
FIRST_FINITE_RQY = (148427, 148428, 148429, 148430)
FIRST_LOCAL_S1 = (148458, 148459, 148460, 148461)
FIRST_CUMULATIVE_S1 = (148462, 148463, 148464, 148465)
FIRST_UPPER_RQY = (148514, 148515, 148516, 148517)
FIRST_INCOMING_M = (104503, 104504, 104505, 104506, 104507, 104508)

MACRO_UPPER_RQY = (141114, 141115, 141116, 141117)
MACRO_LOWER_RQY = (184328, 184329, 184330, 184331)
SEAM_FIRST_LOWER = (184687, 184688)
SEAM_CANCELLATION = (184700, 184702)
PARENT_JOINED = (184685, 184686)
PARENT_UNCONTRACTED = (184679, 184680)

Q4_NODE_COUNT = 183700
MERGED_NODE_COUNT = 184687
SEAM_NODE_COUNT = 184703


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
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"JSON object required: {path.name}")
    return value


def _json_gzip(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"gzip JSON object required: {path.name}")
    return value


def _source_guard(path: Path, functions: Mapping[str, tuple[str, ...]]) -> dict[str, Any]:
    """Read-only AST guard; snippets are intentionally required as tuples."""

    with path.open("r", encoding="utf-8", newline="") as stream:
        raw = stream.read()
    tree = ast.parse(raw, filename=str(path))
    found: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            found[node.name] = ast.get_source_segment(raw, node) or ""
    result: dict[str, Any] = {"path": path.as_posix(), "functions": {}}
    for name, snippets in functions.items():
        _require(isinstance(snippets, tuple), f"source guard snippets must be tuples: {name}")
        body = found.get(name, "")
        compact = "".join(body.split())
        matched = all("".join(token.split()) in compact for token in snippets)
        result["functions"][name] = {
            "present": bool(body),
            "tokens": list(snippets),
            "matched": matched,
        }
        _require(body and matched, f"source guard failed: {path.name}:{name}")
    result["passed"] = True
    result["source_executed"] = False
    return result


def _self_format_guard(root: Path) -> dict[str, Any]:
    raw = (root / MODULE).read_bytes()
    _require(b"\r" not in raw, "new source must use LF line endings")
    _require(raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
             "new source must have exactly one final LF")
    tree = ast.parse(raw.decode("utf-8"), filename=MODULE)
    return {
        "passed": True,
        "LF": True,
        "exactly_one_EOF_LF": True,
        "AST_dump_sha256": _digest(ast.dump(tree, include_attributes=False)),
        "source_executed": False,
    }


def _parent_expected(parent: Mapping[str, Any], path: str) -> str | None:
    if path.startswith(BASE):
        table = parent.get("input_hashes", {})
        if not isinstance(table, Mapping):
            return None
        return table.get(path) or table.get(Path(path).name)
    table = parent.get("implementation_hashes", {})
    if not isinstance(table, Mapping):
        return None
    return table.get(path)


def _pin_parent(root: Path, parent: Mapping[str, Any], path: str) -> dict[str, Any]:
    expected = PARENT_SHA256 if path == PARENT else _parent_expected(parent, path)
    _require(isinstance(expected, str), f"parent has no source pin for {path}")
    actual = _sha(root / path)
    _require(actual == expected, f"changed inherited source: {path}")
    return {"path": path, "expected_sha256": expected, "passed": True}


def _packet_pin(root: Path, parent: Mapping[str, Any], path: str, value: Any) -> dict[str, Any]:
    if isinstance(value, Mapping):
        expected = value.get("expected_sha256")
        self_pin = value.get("self_pin") is True
    else:
        expected = value
        self_pin = False
    _require(isinstance(expected, str), f"numeric packet pin missing: {path}")
    inherited = _parent_expected(parent, path)
    if path == PARENT:
        _require(expected == PARENT_SHA256, "numeric packet parent pin changed")
    elif inherited is not None:
        _require(expected == inherited, f"numeric packet conflicts with parent pin: {path}")
    _require(_sha(root / path) == expected, f"changed numeric packet source: {path}")
    result = {"path": path, "expected_sha256": expected, "passed": True}
    if self_pin:
        result["self_pin"] = True
    return result


def _append_saved(graph: GraphPairs, row: Mapping[str, Any]) -> int:
    _require(isinstance(row, Mapping) and isinstance(row.get("operation"), str),
             "serialized graph row malformed")
    attrs = {key: value for key, value in row.items() if key != "operation"}
    return graph.node(row["operation"], **attrs)


def _row(nodes: Sequence[Mapping[str, Any]], handle: int) -> Mapping[str, Any]:
    _require(type(handle) is int and 0 <= handle < len(nodes),
             f"graph handle outside prefix: {handle}")
    return nodes[handle]


def _has_dependency(nodes: Sequence[Mapping[str, Any]], root: int, target: int) -> bool:
    seen: set[int] = set()
    pending = [root]
    while pending:
        handle = pending.pop()
        if handle in seen:
            continue
        seen.add(handle)
        if handle == target:
            return True
        item = _row(nodes, handle)
        for key in ("argument", "numerator", "denominator", "expression", "variable",
                    "value", "integrand", "lower", "upper"):
            child = item.get(key)
            if type(child) is int and child >= 0:
                pending.append(child)
        arguments = item.get("arguments")
        if isinstance(arguments, list):
            pending.extend(child for child in arguments
                          if type(child) is int and child >= 0)
    return False


def _validate_q4(root: Path) -> dict[str, Any]:
    archive = _json_gzip(root / Q4_ARCHIVE)
    nodes = archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list) and len(nodes) == Q4_NODE_COUNT,
             "Q4 exact graph node count changed")
    bindings = archive.get("modified_Q4_pressure_C3_bindings")
    _require(isinstance(bindings, Mapping), "Q4 bindings missing")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping),
             "Q4 charts/endpoints missing")
    first = charts.get("first_switch")
    macro = charts.get("frozen_macro")
    _require(isinstance(first, Mapping) and isinstance(macro, Mapping),
             "Q4 first-switch/frozen-macro charts missing")
    _require(first.get("physical", {}).get("U") == list(FIRST_U),
             "Q4 first-switch U ports changed")
    _require(first.get("physical", {}).get("U_y") == list(FIRST_UY),
             "Q4 first-switch U_y ports changed")
    _require(first.get("physical", {}).get("m") == list(FIRST_M),
             "Q4 first-switch m ports changed")
    _require(first.get("physical", {}).get("m_y") == list(FIRST_MY),
             "Q4 first-switch m_y ports changed")
    _require(first.get("physical", {}).get("Q") == list(FIRST_Q),
             "Q4 first-switch Q ports changed")
    _require(first.get("physical", {}).get("Q_y") == list(FIRST_QY),
             "Q4 first-switch Q_y ports changed")
    _require(first.get("finite_boundary", {}).get("R_Q_y") == list(FIRST_FINITE_RQY),
             "Q4 first-switch finite R_Q_y ports changed")
    _require(first.get("local_integrals", {}).get("S1") == list(FIRST_LOCAL_S1),
             "Q4 first-switch local S1 ports changed")
    _require(first.get("cumulative_integrals", {}).get("S1") == list(FIRST_CUMULATIVE_S1),
             "Q4 first-switch cumulative S1 ports changed")
    _require(first.get("incoming_integrals", {}).get("S1")[:2] == [141098, 141099],
             "Q4 first-switch incoming S1 changed")
    _require(endpoints.get("first_switch", {}).get("boundary", {}).get("R_Q_y") == list(FIRST_UPPER_RQY),
             "Q4 first-switch upper R_Q_y ports changed")

    for order, handle in enumerate(FIRST_LOCAL_S1[:2]):
        item = _row(nodes, handle)
        _require(
            item.get("operation") == "definite_integral"
            and item.get("original_chart") == "first_switch"
            and item.get("physical_radius_power") == 2
            and item.get("physical_Jacobian_applied_once") is True
            and item.get("original_signed_reduced_pressure") is True
            and item.get("ordinary_slow_Z_derivative_order") == order,
            f"Q4 first-switch S1 source changed at Z {order}",
        )

    for order, handle in enumerate(FIRST_FINITE_RQY[:2]):
        item = _row(nodes, handle)
        _require(item.get("operation") == "product"
                 and item.get("arguments") == [FIRST_QY[order], FIRST_RADIUS_HANDLE],
                 f"Q4 first-switch R_Q_y product changed at Z {order}")
    for order, handle in enumerate(FIRST_UPPER_RQY[:2]):
        item = _row(nodes, handle)
        _require(item.get("operation") == "function_substitution"
                 and item.get("expression") == FIRST_FINITE_RQY[order]
                 and item.get("variable") == FIRST_COORDINATE_HANDLE
                 and item.get("value") == 1
                 and item.get("Z_independent_substitution") is True,
                 f"Q4 first-switch upper endpoint changed at Z {order}")

    nrow = _row(nodes, N_HANDLE)
    _require(nrow.get("operation") == "exact_positive_integer_expression"
             and nrow.get("one_common_N_for_all_charts_and_control_denominators") is True,
             "Q4 common N source changed")
    _require(nrow.get("Z_independent") is True,
             "Q4 common N gained a slow-Z dependence")
    phase = _row(nodes, FIRST_PHASE_HANDLE)
    _require(phase.get("operation") == "analytic_unary"
             and phase.get("name") == "fractional_part",
             "first-switch phase source changed")
    _require(not _has_dependency(nodes, FIRST_PHASE_HANDLE, 20),
             "first-switch phase gained a slow-Z dependency")
    pstar = _row(nodes, PSTAR_HANDLE)
    _require(pstar.get("operation") == "analytic_unary"
             and pstar.get("name") == "exp"
             and pstar.get("argument") == PSTAR_LOG_ARGUMENT_HANDLE,
             "physical Pstar source changed")
    return {
        "passed": True,
        "node_count": len(nodes),
        "physical_ports": {
            "U": list(FIRST_U),
            "U_y": list(FIRST_UY),
            "m": list(FIRST_M),
            "m_y": list(FIRST_MY),
            "Q": list(FIRST_Q),
            "Q_y": list(FIRST_QY),
        },
        "finite_boundary_R_Q_y": list(FIRST_FINITE_RQY),
        "local_S1": list(FIRST_LOCAL_S1[:2]),
        "cumulative_S1": list(FIRST_CUMULATIVE_S1[:2]),
        "upper_R_Q_y": list(FIRST_UPPER_RQY[:2]),
        "ordinary_Z_orders": [0, 1],
        "physical_measure": {
            "radius": FIRST_RADIUS_HANDLE,
            "radius_power": 2,
            "native_J": FIRST_J_HANDLE,
            "R_squared_J_applied_once": True,
        },
        "phase": {
            "handle": FIRST_PHASE_HANDLE,
            "Z_independent": True,
            "total_Z_phase_derivative_exact_zero": True,
        },
        "Pstar": {
            "handle": PSTAR_HANDLE,
            "log_argument_handle": PSTAR_LOG_ARGUMENT_HANDLE,
            "exp_source_not_literal_exp5": True,
        },
    }


def _validate_loop(root: Path) -> dict[str, Any]:
    archive = _json_gzip(root / LOOP_ARCHIVE)
    bindings = archive.get("prescribed_loop_C5_bindings")
    _require(isinstance(bindings, Mapping), "C5 bindings missing")
    _require(bindings.get("N") == N_HANDLE, "C5 common N changed")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping),
             "C5 charts/endpoints missing")
    first = charts.get("first_switch")
    macro = charts.get("frozen_macro")
    _require(isinstance(first, Mapping) and isinstance(macro, Mapping),
             "C5 first-switch/frozen-macro charts missing")
    _require(first.get("actual_phase") == FIRST_PHASE_HANDLE
             and first.get("coordinate") == FIRST_COORDINATE_HANDLE
             and first.get("domain") == [0, 1]
             and first.get("physical_measure") == FIRST_J_HANDLE
             and first.get("radius") == FIRST_RADIUS_HANDLE,
             "C5 first-switch geometry changed")
    _require(macro.get("actual_phase") == MACRO_PHASE_HANDLE
             and macro.get("coordinate") == MACRO_COORDINATE_HANDLE
             and macro.get("physical_measure") == MACRO_J_HANDLE
             and macro.get("radius") == MACRO_RADIUS_HANDLE,
             "C5 macro geometry changed")
    _require(first.get("incoming", {}).get("m")
             == endpoints.get("frozen_macro", {}).get("m"),
             "first-switch incoming m is not the frozen-macro endpoint")
    _require(first.get("incoming", {}).get("m") == list(FIRST_INCOMING_M),
             "first-switch incoming m handles changed")
    definition = bindings.get("definition", {})
    _require(definition.get("correction_zero_inlet_only_at_rminus_background_nonzero_memories_unchanged") is True,
             "C5 nonreset incoming history premise changed")
    return {
        "passed": True,
        "common_N_handle": N_HANDLE,
        "first_switch": {
            "phase": FIRST_PHASE_HANDLE,
            "coordinate": FIRST_COORDINATE_HANDLE,
            "domain": [0, 1],
            "native_J": FIRST_J_HANDLE,
            "radius": FIRST_RADIUS_HANDLE,
            "incoming_m": list(first["incoming"]["m"]),
            "history_m": list(first.get("history", {}).get("m", [])),
        },
        "frozen_macro": {
            "phase": MACRO_PHASE_HANDLE,
            "coordinate": MACRO_COORDINATE_HANDLE,
            "native_J": MACRO_J_HANDLE,
            "radius": MACRO_RADIUS_HANDLE,
            "endpoint_m": list(endpoints["frozen_macro"]["m"]),
        },
        "same_common_N": True,
        "phase_Z_independent": True,
        "nonreset_incoming_history": True,
    }


def _validate_seam(parent: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    seam = parent.get("actual_macro_switch_RQy_seam")
    _require(isinstance(seam, Mapping) and seam.get("passed") is True,
             "accepted macro/first-switch seam packet missing")
    guards = seam.get("source_guards", {})
    _require(guards.get("common_N") == N_HANDLE,
             "seam common N changed")
    _require(guards.get("native_J_distinct") == {
        "frozen_macro": MACRO_J_HANDLE,
        "first_switch": FIRST_J_HANDLE,
    }, "seam native J distinction changed")
    _require(guards.get("phase_handles") == {
        "frozen_macro": MACRO_PHASE_HANDLE,
        "first_switch": FIRST_PHASE_HANDLE,
    }, "seam phase handles changed")
    _require(guards.get("nonreset_M") == {
        "frozen_macro": True,
        "first_switch": True,
    }, "seam nonreset M premise changed")
    patch = seam.get("graph_patch")
    _require(isinstance(patch, Mapping)
             and patch.get("base_node_count") == MERGED_NODE_COUNT
             and patch.get("appended_node_count") == 16,
             "accepted seam suffix changed")
    endpoints = patch.get("endpoint_handles", {})
    _require(endpoints.get("first_switch_lower_new_s0_Z0_Z1") == list(SEAM_FIRST_LOWER),
             "accepted first-switch lower seam handles changed")
    _require(endpoints.get("signed_cancellation") == list(SEAM_CANCELLATION),
             "accepted seam cancellation handles changed")
    for k, handle in enumerate(SEAM_CANCELLATION):
        item = _row(nodes, handle)
        _require(item.get("operation") == "sum"
                 and item.get("arguments") == [
                     _negative_child(nodes, MACRO_UPPER_RQY[k]), SEAM_FIRST_LOWER[k]
                 ], "accepted macro-upper/first-lower cancellation changed")
    return {
        "passed": True,
        "common_N_handle": N_HANDLE,
        "macro_native_J": MACRO_J_HANDLE,
        "first_switch_native_J": FIRST_J_HANDLE,
        "native_J_distinct": True,
        "phase_handles": {
            "frozen_macro": MACRO_PHASE_HANDLE,
            "first_switch": FIRST_PHASE_HANDLE,
        },
        "nonreset_M": True,
        "conditional_identity_only_in_joined_expression": True,
        "seam_patch_base_count": MERGED_NODE_COUNT,
        "seam_patch_count": 16,
    }


def _negative_child(nodes: Sequence[Mapping[str, Any]], target: int) -> int:
    matches = [index for index, item in enumerate(nodes)
               if item.get("operation") == "negative" and item.get("argument") == target]
    _require(matches, f"negative endpoint trace missing for {target}")
    return matches[-1]


def _rehydrate_prefix(root: Path, parent: Mapping[str, Any], q4: Mapping[str, Any]) -> tuple[GraphPairs, dict[str, Any]]:
    q4_nodes = q4["exact_graph_nodes"]
    correlated = parent["actual_micro_macro_correlated_m15"]["merged_namespace"]["appended_nodes"]
    seam_rows = parent["actual_macro_switch_RQy_seam"]["graph_patch"]["appended_nodes"]
    _require(isinstance(correlated, list) and len(correlated) == 987,
             "complete micro/macro merged suffix is not 987 rows")
    _require(isinstance(seam_rows, list) and len(seam_rows) == 16,
             "complete macro/first-switch seam suffix is not 16 rows")
    graph = GraphPairs(q4_nodes)
    for row in correlated:
        _append_saved(graph, row)
    _require(len(graph.nodes) == MERGED_NODE_COUNT,
             "merged micro/macro prefix is not 184687 nodes")
    merged_snapshot = copy.deepcopy(graph.nodes)
    for row in seam_rows:
        _append_saved(graph, row)
    _require(len(graph.nodes) == SEAM_NODE_COUNT,
             "full current prefix is not 184703 nodes")
    _require(graph.nodes[:Q4_NODE_COUNT] == q4_nodes,
             "Q4 graph prefix changed")
    _require(graph.nodes[Q4_NODE_COUNT:MERGED_NODE_COUNT] == correlated,
             "micro/macro merged suffix changed")
    _require(graph.nodes[MERGED_NODE_COUNT:] == seam_rows,
             "macro/first-switch seam suffix changed")
    namespace = parent["actual_micro_macro_correlated_m15"]["graph_namespace"]
    safe = parent["actual_micro_macro_correlated_m15"]["safe_graphs"]
    old_digest_payload = {
        "base_sha256": _sha(root / Q4_ARCHIVE),
        "base_node_count": Q4_NODE_COUNT,
        "correlated_suffix": safe["accepted_correlated"]["appended_nodes"],
        "relocated_rqy_suffix": safe["relocated_rqy"]["appended_nodes"],
        "macro_suffix": safe["new_macro"]["appended_nodes"],
        "joined_suffix": safe["joined_expression"]["appended_nodes"],
    }
    _require(_digest(old_digest_payload) == namespace.get("namespace_digest"),
             "old 184687 namespace digest changed")
    _require(namespace.get("merged_node_count") == MERGED_NODE_COUNT,
             "old merged node count changed")
    return graph, {
        "Q4_node_count": Q4_NODE_COUNT,
        "merged_node_count": MERGED_NODE_COUNT,
        "seam_node_count": SEAM_NODE_COUNT,
        "old_namespace_digest": namespace["namespace_digest"],
        "merged_snapshot": merged_snapshot,
        "full_prefix_snapshot": copy.deepcopy(graph.nodes),
        "full_prefix_canonical_nodes_sha256": _digest(graph.nodes),
        "segment_starts": {
            "q4": 0,
            "micro_macro_merged": Q4_NODE_COUNT,
            "macro_first_switch_seam": MERGED_NODE_COUNT,
            "new_suffix": SEAM_NODE_COUNT,
        },
        "suffix_counts": {
            "micro_macro_merged": len(correlated),
            "macro_first_switch_seam": len(seam_rows),
        },
    }


def _pair(value: Sequence[int], label: str) -> tuple[int, int]:
    _require(isinstance(value, (list, tuple)) and len(value) == 2
             and all(type(item) is int for item in value),
             f"{label} must be an ordinary Z 0 / 1 pair")
    return int(value[0]), int(value[1])


def _build_six_source(graph: GraphPairs) -> tuple[dict[str, list[int]], dict[str, Any]]:
    q = _pair(FIRST_Q[:2], "first Q")
    q_z = _pair(FIRST_Q[1:3], "first Q_Z")
    u = _pair(FIRST_U[:2], "first U")
    u_z = _pair(FIRST_U[1:3], "first U_Z")
    source_pairs = {
        "q": q,
        "q_Z": q_z,
        "q_squared": graph.pmul(q, q),
        "U_q_Z": graph.pmul(u, q_z),
        "U_squared": graph.pmul(u, u),
        "U_U_Z": graph.pmul(u, u_z),
    }
    weight = graph.mul(FIRST_RADIUS_HANDLE, FIRST_RADIUS_HANDLE, FIRST_J_HANDLE)
    integrals: dict[str, list[int]] = {}
    for name in I1_NAMES:
        weighted = graph.pscale(source_pairs[name], weight)
        integral = graph.integral(
            weighted,
            FIRST_COORDINATE_HANDLE,
            0,
            1,
            f"actual_first_switch_finite_N_I1_{name}",
            metadata={
                "original_chart": "first_switch",
                "source_layer": "physical_six_I1_density",
                "ordinary_Z_rows": [0, 1],
                "physical_radius_power": 2,
                "physical_Jacobian_applied_once": True,
                "Pstar_applied_in_physical_fields_once": True,
                "endpoint_independent_of_Z": True,
                "physical_y_derivatives_claimed": False,
            },
        )
        integrals[name] = list(integral)
    return integrals, {
        "passed": True,
        "density_pairs": {name: list(pair) for name, pair in source_pairs.items()},
        "weight_handle": weight,
        "weight_formula": "R_first^2 * J_first",
        "weight_units": "R^2*J*dcoordinate",
        "graph_source_raw_physical": True,
        "numeric_packet_normalization": "Wswitch=10000*hbS",
        "integral_handles": integrals,
        "ordinary_Z_orders": [0, 1],
        "Pstar_applied_once": True,
        "R_not_a_mean_scale": True,
    }


def _build_reduced_D(graph: GraphPairs, integrals: Mapping[str, Sequence[int]]) -> tuple[list[int], dict[str, Any]]:
    z = _pair((20, 1), "Z")
    d = _pair((184430, 184429), "d")
    h = _pair((184432, 184431), "h")
    L = _pair((184439, 184438), "L")
    L2 = _pair((184440, 184443), "L^2")
    one = graph.fixed(1)
    half = 28
    two = 48

    inv_L = graph.pdiv(one, L, "same positive physical L=1-delta*Z^2")
    minus_inv_L = graph.pneg(inv_L)
    hZ_over_2L = graph.pscale(
        graph.pdiv(graph.pmul(h, z), L, "same positive physical L=1-delta*Z^2"),
        half,
    )
    minus_half = graph.pneg(graph.fixed(half))
    d_over_L = graph.pdiv(d, L, "same positive physical L=1-delta*Z^2")
    two_hZ2_over_L2 = graph.pdiv(
        graph.pscale(graph.pmul(h, graph.pmul(z, z)), two),
        L2,
        "same positive physical L=1-delta*Z^2",
    )
    two_Zd_over_L2 = graph.pdiv(
        graph.pscale(graph.pmul(z, d), two),
        L2,
        "same positive physical L=1-delta*Z^2",
    )
    terms = {
        "minus_inv_L": graph.pmul(minus_inv_L, _pair(integrals["q"], "q integral")),
        "hZ_over_2L": graph.pmul(hZ_over_2L, _pair(integrals["q_Z"], "q_Z integral")),
        "minus_half": graph.pmul(minus_half, _pair(integrals["q_squared"], "q_squared integral")),
        "d_over_L": graph.pmul(d_over_L, _pair(integrals["U_q_Z"], "U_q_Z integral")),
        "two_hZ2_over_L2": graph.pmul(two_hZ2_over_L2, _pair(integrals["U_squared"], "U_squared integral")),
        "two_Zd_over_L2": graph.pmul(two_Zd_over_L2, _pair(integrals["U_U_Z"], "U_U_Z integral")),
    }
    D = graph.padd(*terms.values())
    return list(D), {
        "passed": True,
        "rows": list(D),
        "terms": {name: list(pair) for name, pair in terms.items()},
        "coefficients": {
            "minus_inv_L": list(minus_inv_L),
            "hZ_over_2L": list(hZ_over_2L),
            "minus_half": list(minus_half),
            "d_over_L": list(d_over_L),
            "two_hZ2_over_L2": list(two_hZ2_over_L2),
            "two_Zd_over_L2": list(two_Zd_over_L2),
        },
        "formula": (
            "D=-q/L+h*Z*q_Z/(2L)-q^2/2+d*U*q_Z/L+"
            "2*h*Z^2*U^2/L^2+2*Z*d*U*U_Z/L^2"
        ),
        "graph_source_unit": "raw physical R^2*J*dcoordinate",
        "numeric_packet_unit": "Wswitch",
        "normalization_not_in_graph": True,
        "D_is_not_full_S1": True,
    }


def _build_reduced_B(graph: GraphPairs) -> tuple[list[int], dict[str, Any]]:
    """Build the endpoint B source with its physical R^2 factor only."""

    z = _pair((20, 1), "Z")
    u = _pair(FIRST_U[:2], "first U")
    q = _pair(FIRST_Q[:2], "first Q")
    L = _pair((184439, 184438), "L")
    L2 = _pair((184440, 184443), "L^2")
    radius_squared = graph.mul(FIRST_RADIUS_HANDLE, FIRST_RADIUS_HANDLE)
    half = 28
    two = 48
    first_term = graph.pscale(graph.pdiv(q, L, "same positive physical L=1-delta*Z^2"), radius_squared)
    second_term = graph.pscale(graph.pmul(q, q), radius_squared)
    second_term = graph.pscale(second_term, half)
    third_numerator = graph.pscale(
        graph.pmul(graph.pmul(z, z), graph.fixed(radius_squared), graph.pmul(u, u)),
        two,
    )
    third_term = graph.pneg(
        graph.pdiv(third_numerator, L2, "same positive physical L=1-delta*Z^2")
    )
    boundary_integrand = graph.padd(first_term, second_term, third_term)
    upper = graph.at(boundary_integrand, FIRST_COORDINATE_HANDLE, 1)
    lower = graph.at(boundary_integrand, FIRST_COORDINATE_HANDLE, 0)
    B = graph.padd(upper, graph.pneg(lower))
    return list(B), {
        "passed": True,
        "rows": list(B),
        "boundary_integrand": list(boundary_integrand),
        "upper": list(upper),
        "lower": list(lower),
        "formula": (
            "B=[R^2*Q/L+R^2*Q^2/2-2*Z^2*R^2*U^2/L^2]"
            "_{s=0}^{s=1}"
        ),
        "physical_radius_power": 2,
        "physical_Jacobian_applied": False,
        "endpoint_measure_has_no_J": True,
        "graph_source_unit": "raw physical R^2*Q endpoint",
        "numeric_endpoint_not_evaluated": True,
        "ordinary_Z_orders": [0, 1],
    }


def _numeric_packet_contract(packet: Mapping[str, Any]) -> dict[str, Any]:
    _require(packet.get("passed") is True, "first-switch numeric packet is not passed")
    _require(packet.get("parent_report_sha256") == PARENT_SHA256,
             "first-switch numeric packet has the wrong parent SHA")
    _require(isinstance(packet.get("source_pins"), Mapping),
             "first-switch numeric packet source_pins missing")
    for key in ("I1_full", "D_full", "joined_macro_first_switch_volume"):
        _require(isinstance(packet.get(key), Mapping) and packet[key].get("passed") is True,
                 f"first-switch numeric {key} missing")
    _require(packet["I1_full"].get("normalized_by_Wswitch") is True,
             "first-switch I1 packet is not in Wswitch units")
    _require(packet["D_full"].get("normalized_by_Wswitch") is True,
             "first-switch D packet is not in Wswitch units")
    source_contract = packet.get("source_contract", {})
    _require(source_contract.get("chart") == "first_switch"
             and source_contract.get("N_handle") == N_HANDLE
             and source_contract.get("physical_J_handle") == FIRST_J_HANDLE
             and source_contract.get("phase_Z_independent") is True
             and source_contract.get("physical_y_derivatives_claimed") is False
             and source_contract.get("nonreset_mean_positive_kernel") is True,
             "first-switch numeric source contract changed")
    limits = packet.get("limits", {})
    _require(limits.get("outer_RQy_retained") is True
             and limits.get("outer_RQy_evaluated") is False
             and limits.get("complete_pressure_global_closure") is False,
             "first-switch numeric scope changed")
    joined = packet["joined_macro_first_switch_volume"]
    _require(joined.get("complete_pressure_enclosed") is False
             and joined.get("outer_RQy_retained") is True,
             "joined numeric volume scope changed")
    joined_I1 = joined.get("I1", {})
    joined_D = joined.get("D", {})
    _require(joined_I1.get("normalized_by_Wmacro") is True
             and joined_D.get("normalized_by_Wmacro") is True,
             "joined numeric volume units changed")
    return {
        "passed": True,
        "parent_report_sha256": PARENT_SHA256,
        "source_pins_bound": True,
        "source_pins": copy.deepcopy(dict(packet["source_pins"])),
        "source_contract": copy.deepcopy(dict(packet.get("source_contract", {}))),
        "I1_full": copy.deepcopy(packet["I1_full"]),
        "D_full": copy.deepcopy(packet["D_full"]),
        "joined_macro_first_switch_volume": copy.deepcopy(joined),
        "physical_error_log_bounds": copy.deepcopy(
            packet.get("physical_error_log_bounds", {})
        ),
        "six_I1_difference_error_log_bounds": copy.deepcopy(
            packet.get("six_I1_difference_error_log_bounds", {})
        ),
        "signed_I1_radius": copy.deepcopy(packet.get("signed_I1_radius")),
        "I1_error_threshold": copy.deepcopy(packet.get("I1_error_threshold")),
        "limits": copy.deepcopy(packet.get("limits", {})),
        "D_is_not_full_S1": True,
        "full_outer_RQy_numeric": False,
        "complete_pressure_numeric": False,
    }


def _build_join(graph: GraphPairs, parent: Mapping[str, Any], local_s1: Sequence[int]) -> dict[str, Any]:
    joined_packet = parent["actual_micro_macro_correlated_m15"]["actual_joined_pressure_function"]
    _require(joined_packet.get("joined") == list(PARENT_JOINED),
             "parent whole micro/macro joined expression changed")
    _require(joined_packet.get("uncontracted_joined") == list(PARENT_UNCONTRACTED),
             "parent uncontracted joined expression changed")
    _require(joined_packet.get("macro_outer_RQy_retained") is True,
             "parent macro outer R_Q_y was not retained")
    for k, root in enumerate(PARENT_JOINED):
        _require(_has_dependency(graph.nodes, root, MACRO_UPPER_RQY[k]),
                 "parent joined row lost macro-upper trace")
        _require(_has_dependency(graph.nodes, root, 183935 + k),
                 "parent joined row lost first-micro lower trace")

    half_s1 = [graph.mul(local_s1[k], 28) for k in range(2)]
    first_upper_negative = [graph.neg(FIRST_UPPER_RQY[k]) for k in range(2)]
    parent_without_macro_upper = [
        graph.add(PARENT_JOINED[k], MACRO_UPPER_RQY[k]) for k in range(2)
    ]
    contracted_trace = [
        graph.add(parent_without_macro_upper[k], SEAM_CANCELLATION[k])
        for k in range(2)
    ]
    joined_contracted = [
        graph.add(contracted_trace[k], half_s1[k], first_upper_negative[k])
        for k in range(2)
    ]
    first_increment_uncontracted = [
        graph.add(half_s1[k], SEAM_FIRST_LOWER[k], first_upper_negative[k])
        for k in range(2)
    ]
    joined_uncontracted = [
        graph.add(PARENT_JOINED[k], first_increment_uncontracted[k])
        for k in range(2)
    ]
    return {
        "passed": True,
        "ordinary_Z_orders": [0, 1],
        "parent_joined": list(PARENT_JOINED),
        "parent_uncontracted": list(PARENT_UNCONTRACTED),
        "first_switch_local_S1": list(local_s1),
        "first_switch_half_S1": half_s1,
        "first_switch_upper_RQy": list(FIRST_UPPER_RQY[:2]),
        "first_switch_lower_RQy": list(SEAM_FIRST_LOWER),
        "macro_upper_RQy": list(MACRO_UPPER_RQY[:2]),
        "macro_lower_RQy_standalone": list(MACRO_LOWER_RQY),
        "seam_cancellation": list(SEAM_CANCELLATION),
        "parent_without_macro_upper": parent_without_macro_upper,
        "contracted_internal_trace": contracted_trace,
        "first_increment_uncontracted": first_increment_uncontracted,
        "joined_uncontracted": joined_uncontracted,
        "joined_contracted": joined_contracted,
        "formula_uncontracted": (
            "parent_joined + first_local_S1/2 + R_Q_y(first_switch,s=0) "
            "- R_Q_y(first_switch,s=1)"
        ),
        "formula_contracted": (
            "parent_joined with internal -R_Q_y(macro,m=1) + "
            "R_Q_y(first_switch,s=0) contracted, plus first_local_S1/2 "
            "- R_Q_y(first_switch,s=1)"
        ),
        "conditional_identity_only_in_joined_expression": True,
        "first_switch_upper_standalone_retained": True,
        "macro_lower_standalone_retained": True,
        "outer_RQy_and_caps_retained": True,
        "signed_relation_not_global_closure": True,
    }


def build(root: str | Path, parent: Mapping[str, Any], first_switch_numeric_packet: Mapping[str, Any]) -> dict[str, Any]:
    """Bind first-switch finite-N source and append a joined pressure suffix."""

    root_path = Path(root)
    _require(_sha(root_path / PARENT) == PARENT_SHA256,
             "current macro finite-N parent file changed")
    saved_parent = _json(root_path / PARENT)
    _require(_digest(saved_parent) == _digest(parent),
             "exact current macro finite-N parent required")
    _require(isinstance(first_switch_numeric_packet, Mapping),
             "first-switch numeric packet mapping required")

    pin_paths = (Q4_SOURCE, Q4_ARCHIVE, Q4_CHECK, LOOP_SOURCE, LOOP_ARCHIVE,
                 LOOP_CHECK, DENSITY_SOURCE, GRAPH_HELPER, SEAM_SOURCE)
    source_pins: dict[str, Any] = {PARENT: _pin_parent(root_path, parent, PARENT)}
    for path in pin_paths:
        source_pins[path] = _pin_parent(root_path, parent, path)
    for path, value in first_switch_numeric_packet["source_pins"].items():
        pin = _packet_pin(root_path, parent, path, value)
        previous = source_pins.get(path)
        if previous is not None:
            _require(previous["expected_sha256"] == pin["expected_sha256"],
                     f"conflicting source pin: {path}")
        source_pins[path] = pin
    source_pins[MODULE] = {
        "path": MODULE,
        "expected_sha256": _sha(root_path / MODULE),
        "passed": True,
        "self_pin": True,
    }

    q4_guards = _validate_q4(root_path)
    loop_guards = _validate_loop(root_path)
    q4_archive = _json_gzip(root_path / Q4_ARCHIVE)
    graph, prefix = _rehydrate_prefix(root_path, parent, q4_archive)
    seam_guards = _validate_seam(parent, graph.nodes)
    numeric = _numeric_packet_contract(first_switch_numeric_packet)
    ast_guards = {
        "q4_physical_fields": _source_guard(root_path / Q4_SOURCE, {
            "physical_pressure_fields": (
                "U=a.scale(fields['V'],Pstar)",
                "Uy=a.scale(radial['V_y'],Pstar)",
                "m=a.scale(fields['m'],Pstar)",
                "my=a.sub(U,m)",
                "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
            ),
        }),
        "q4_pressure_measure": _source_guard(root_path / Q4_SOURCE, {
            "__init__": (
                "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,chart['radius']))",
                "physical_measures='S0:R*Jacobian*dcoordinate;S1:R^2*Jacobian*dcoordinate'",
                "m15_minus_endpoint_sign_applied_once=True",
            ),
        }),
        "c5_finite_N": _source_guard(root_path / LOOP_SOURCE, {
            "finite_N_functions": (
                "s=a.scale(loop['B_over_Pstar'],invN)",
                "V_N=a.add(V,s)",
            ),
        }),
        "c5_transport": _source_guard(root_path / LOOP_SOURCE, {
            "transport": (
                "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
                "physical_Jacobian_applied_once=True",
                "predecessor_memory_never_reset=True",
            ),
        }),
        "c5_phase_and_history": _source_guard(root_path / LOOP_SOURCE, {
            "__init__": (
                "self.N=self.unique('exact_positive_integer_expression',None)",
                "phase=g.unary('fractional_part',g.mul(self.N,offset))",
                "incoming=self.endpoints[name]",
            ),
        }),
    }
    guards = {
        "passed": True,
        "AST": ast_guards,
        "self_format": _self_format_guard(root_path),
        "Q4": q4_guards,
        "C5": loop_guards,
        "seam": seam_guards,
        "source_execution": False,
    }

    integrals, six_source = _build_six_source(graph)
    B_rows, B_source = _build_reduced_B(graph)
    D_rows, D_source = _build_reduced_D(graph, integrals)
    reduced_S1 = list(graph.padd(_pair(B_rows, "B rows"), _pair(D_rows, "D rows")))
    join = _build_join(graph, parent, FIRST_LOCAL_S1[:2])
    new_suffix = copy.deepcopy(graph.nodes[SEAM_NODE_COUNT:])
    _require(graph.nodes[:SEAM_NODE_COUNT] == prefix["full_prefix_snapshot"],
             "full 184703-node prefix changed")
    _require(new_suffix == graph.appended[len(graph.appended) - len(new_suffix):],
             "new graph suffix serialization changed")

    source_contract = {
        "passed": True,
        "chart": "first_switch",
        "finite_N_source": True,
        "finite_N_fields_from_C5_loop": True,
        "ordinary_Z_orders": [0, 1],
        "common_N_handle": N_HANDLE,
        "N_not_materialized": True,
        "phase_handles": {
            "frozen_macro": MACRO_PHASE_HANDLE,
            "first_switch": FIRST_PHASE_HANDLE,
        },
        "phase_Z_independent": True,
        "native_J": {
            "frozen_macro": MACRO_J_HANDLE,
            "first_switch": FIRST_J_HANDLE,
            "distinct": True,
        },
        "physical_y_derivatives_claimed": False,
        "physical_Pstar_applied_once": True,
        "Pstar_handle": PSTAR_HANDLE,
        "Pstar_log_argument_handle": PSTAR_LOG_ARGUMENT_HANDLE,
        "physical_R_not_a_mean_scale": True,
        "physical_measure": "R_first^2*J_first*ds",
        "six_I1_D_unit": "raw physical graph; numeric packet Wswitch",
        "joined_volume_unit": "Wmacro",
        "full_outer_RQy_numeric": False,
        "B_endpoint_numeric": False,
        "complete_pressure_numeric": False,
        "endpoint_measure": "R_first*Q_y (no J)",
        "nonreset_incoming_m": True,
        "incoming_m_path": "C5.first_switch.incoming.m = C5.endpoints.frozen_macro.m",
        "incoming_m_handles": list(FIRST_INCOMING_M),
        "local_S1_added_once": True,
        "numeric_D_is_not_full_S1": True,
        "conditional_RQy_identity_only_in_joined_expression": True,
    }
    first_switch_source = {
        "passed": True,
        "finite_N": True,
        "finite_N_source": "Q4 first_switch physical fields with C5 B_over_Pstar/N correction",
        "physical_fields": {
            "U": list(FIRST_U),
            "U_y": list(FIRST_UY),
            "m": list(FIRST_M),
            "m_y": list(FIRST_MY),
            "Q": list(FIRST_Q),
            "Q_y": list(FIRST_QY),
        },
        "common_N": N_HANDLE,
        "phase": FIRST_PHASE_HANDLE,
        "phase_Z_independent": True,
        "incoming_m": list(FIRST_INCOMING_M),
        "nonreset_incoming_history": True,
        "six_I1": {
            "names": list(I1_NAMES),
            "source": six_source,
            "numeric_packet_I1_full": copy.deepcopy(numeric["I1_full"]),
            "graph_source_unit": "raw physical R^2*J*dcoordinate",
            "numeric_unit": "Wswitch (packet only)",
            "graph_source_not_normalized": True,
            "physical_error_radius_is_not_a_selected_value": True,
        },
        "ordinary_Z_orders": [0, 1],
        "six_I1_names": list(I1_NAMES),
        "coordinate": FIRST_COORDINATE_HANDLE,
        "domain": [0, 1],
        "physical_measure": {
            "radius": FIRST_RADIUS_HANDLE,
            "radius_power": 2,
            "length": FIRST_J_HANDLE,
            "weight": six_source["weight_handle"],
            "R_squared_length_applied_once": True,
            "Z_independent": True,
        },
        "source": copy.deepcopy(six_source["density_pairs"]),
        "integrals": copy.deepcopy(six_source["integral_handles"]),
        "endpoint_fields": {
            "lower": {
                "R_Q_y": list(SEAM_FIRST_LOWER),
                "finite_source": list(FIRST_FINITE_RQY[:2]),
                "coordinate": FIRST_COORDINATE_HANDLE,
                "value": 0,
                "J_applied": False,
            },
            "upper": {
                "R_Q_y": list(FIRST_UPPER_RQY[:2]),
                "coordinate": FIRST_COORDINATE_HANDLE,
                "value": 1,
                "J_applied": False,
            },
        },
        "formula_context": {
            "z": [20, 1],
            "delta": [25220, 0],
            "L": [184439, 184438],
            "L2": [184440, 184443],
            "two_id": 48,
        },
        "D": {
            "source": D_source,
            "numeric_packet_D_full": copy.deepcopy(numeric["D_full"]),
            "graph_source_unit": "raw physical R^2*J*dcoordinate",
            "numeric_unit": "Wswitch (packet only)",
            "graph_source_not_normalized": True,
            "numeric_D_is_not_full_S1": True,
        },
        "reduced_boundary_B": B_rows,
        "B_source": B_source,
        "reduced_volume_D": D_rows,
        "reduced_S1": reduced_S1,
        "B": {
            "source_formula": B_source["formula"],
            "boundary_F_formula": "F=Q/(2L)+Q^2/4-Z^2*U^2/L^2",
            "source_ports": {
                "Q": list(FIRST_Q[:2]),
                "U": list(FIRST_U[:2]),
                "radius": FIRST_RADIUS_HANDLE,
                "coordinate": FIRST_COORDINATE_HANDLE,
                "lower": 0,
                "upper": 1,
            },
            "numericized": False,
            "physical_R_squared_applied_once": True,
        },
        "S1": {
            "source_handles": list(FIRST_LOCAL_S1[:2]),
            "reduced_B_plus_D_source": reduced_S1,
            "measure": "R^2*J*ds",
            "graph_source_unit": "raw physical R^2*J*ds",
            "numeric_joined_unit": "Wmacro (packet only)",
            "direct_Q4_local_integral": True,
            "conditional_identity": "S1 = B + D under the Q4 FTC source identity",
            "D_is_not_full_S1": True,
            "numeric_packet_joined_volume": copy.deepcopy(
                numeric["joined_macro_first_switch_volume"]
            ),
        },
    }

    graph_patch = {
        "base_archive": Q4_ARCHIVE,
        "base_node_count": SEAM_NODE_COUNT,
        "base_prefix_unchanged": graph.nodes[:SEAM_NODE_COUNT] == prefix["full_prefix_snapshot"],
        "base_parent_merged_namespace_digest": prefix["old_namespace_digest"],
        "base_prefix_canonical_nodes_sha256": prefix["full_prefix_canonical_nodes_sha256"],
        "new_node_start": SEAM_NODE_COUNT,
        "appended_nodes": new_suffix,
        "appended_node_count": len(new_suffix),
        "full_node_count": len(graph.nodes),
        "full_canonical_nodes_sha256": _digest(graph.nodes),
        "segment_starts": prefix["segment_starts"],
        "prefix_suffix_counts": prefix["suffix_counts"],
        "endpoint_handles": {
            "macro_lower_standalone": list(MACRO_LOWER_RQY),
            "macro_upper_existing": list(MACRO_UPPER_RQY),
            "first_switch_lower_existing": list(SEAM_FIRST_LOWER),
            "first_switch_upper_standalone": list(FIRST_UPPER_RQY),
            "accepted_seam_cancellation": list(SEAM_CANCELLATION),
            "first_switch_local_S1": list(FIRST_LOCAL_S1[:2]),
            "first_switch_half_S1": join["first_switch_half_S1"],
            "joined_uncontracted": join["joined_uncontracted"],
            "joined_contracted": join["joined_contracted"],
        },
    }

    limits = {
        "first_switch_six_I1_D_numeric_packet_passed": True,
        "first_switch_full_S1_source_defined": True,
        "first_switch_D_is_not_full_S1": True,
        "joined_pressure_source_defined": True,
        "joined_pressure_numeric_enclosure": False,
        "outer_RQy_retained": True,
        "outer_RQy_evaluated": False,
        "macro_lower_standalone_retained": True,
        "first_switch_upper_standalone_retained": True,
        "all_caps_retained": True,
        "conditional_RQy_identity_only_in_joined_expression": True,
        "global_closure": False,
        "complete_pressure_global_closure": False,
        "accepted_ancestors_replayed": False,
        "old_212_terms_and_errors_unchanged": True,
        "physical_y_derivatives_claimed": False,
        "N_materialized": False,
    }
    result = {
        "passed": True,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": parent.get("source_family"),
        "source_pins": source_pins,
        "input_pins": copy.deepcopy(source_pins),
        "source_guards": guards,
        "source_contract": source_contract,
        "contract": copy.deepcopy(source_contract),
        "first_switch_numeric_packet": numeric,
        "I1_full": copy.deepcopy(numeric["I1_full"]),
        "D_full": copy.deepcopy(numeric["D_full"]),
        "joined_macro_first_switch_volume": copy.deepcopy(
            numeric["joined_macro_first_switch_volume"]
        ),
        "errors": {
            "physical_error_log_bounds": copy.deepcopy(
                numeric["physical_error_log_bounds"]
            ),
            "six_I1_difference_error_log_bounds": copy.deepcopy(
                numeric["six_I1_difference_error_log_bounds"]
            ),
        },
        "first_switch_source": first_switch_source,
        "joined_pressure_source": join,
        "graph_patch": graph_patch,
        "limits": limits,
        "signed_relation_not_global_closure": True,
        "accepted_ancestors_replayed": False,
    }
    json.dumps(result, sort_keys=True, ensure_ascii=True, allow_nan=False)
    return result


__all__ = ["build"]
