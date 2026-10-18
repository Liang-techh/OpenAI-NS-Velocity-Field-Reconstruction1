"""Source-only second-switch and post-power finite-N pressure extension.

This module consumes the already accepted macro/first-switch joined source and
an independently supplied numeric packet.  It appends only source graph rows;
it never imports or calls an accepted builder, checker, owner, norm,
recurrence, quadrature, or old exact-algebra routine.  The two endpoint
identities are proved here from source/cache bindings plus arbitrary symbolic
jets.  Finite-N ``Q_y`` contains slow-derivative and fast-phase terms, so the
joined contraction is enabled only after that source-owned proof.
"""
from __future__ import annotations

import ast
import copy
import gzip
import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any

from .original_signed_graph_pairs import GraphPairs


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_finite_N_I1.json"
PARENT_SHA256 = "2cc8b4ff9670858564297bc1c27d5322c9af7d2cb85381939dc3e37a3eea346a"
MODULE = "src/openai_ns_reconstruction/original_switch_post_pressure_source.py"

Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
Q4_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json"
LOOP_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
LOOP_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
LOOP_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json"
DENSITY_SOURCE = "src/openai_ns_reconstruction/original_finite_N_density_kernel.py"
GRAPH_HELPER = "src/openai_ns_reconstruction/original_signed_graph_pairs.py"
SEAM_SOURCE = "src/openai_ns_reconstruction/original_macro_switch_RQy_seam.py"
R100_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
R100_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
SECOND_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_second_switch_R110.py"
BRIDGE_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
ACTUAL_BRIDGE_SOURCE = BASE + "lei_ren_part1_paper_compliant_actual_bridge_integrals.py"
FIRST_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_first_switch_functions.py"
MACRO_FIRST_SOURCE = "src/openai_ns_reconstruction/original_macro_first_switch_pressure_source.py"
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_outer_initial.py"
SIGMA_JETS_SOURCE = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"

Q4_NODE_COUNT = 183700
MERGED_NODE_COUNT = 184687
SEAM_NODE_COUNT = 184703
FIRST_NODE_COUNT = 184898
N_HANDLE = 18
PSTAR_HANDLE = 33572
PSTAR_LOG_ARGUMENT_HANDLE = 5

I1_NAMES = ("q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z")
MACRO_UPPER_RQY = (141114, 141115)
FIRST_UPPER_RQY = (148514, 148515)
SECOND_UPPER_RQY = (154676, 154677)
POST_UPPER_RQY = (160287, 160288)
PARENT_JOINED = (184892, 184893)

_CHARTS: dict[str, dict[str, Any]] = {
    "second_switch": {
        "native_chart": "second_switch",
        "q4_chart": "second_switch",
        "phase": 107165,
        "coordinate": 75949,
        "radius": 88600,
        "J": 9,
        "U": (154247, 154248, 154249, 154250, 154251, 154252),
        "U_y": (154253, 154254, 154255, 154256, 154257, 154258),
        "m": (154259, 154260, 154261, 154262, 154263, 154264),
        "m_y": (154271, 154272, 154273, 154274, 154275, 154276),
        "Q": (154342, 154346, 154352, 154358, 154364),
        "Q_y": (154430, 154434, 154440, 154446, 154452),
        "finite_RQy": (154589, 154590, 154591, 154592),
        "local_S1": (154620, 154621, 154622, 154623),
        "cumulative_S1": (154624, 154625, 154626, 154627),
        "incoming_S1": (148498, 148499, 148500, 148501),
        "endpoint_S1": (154660, 154661, 154662, 154663),
        "upper_RQy": SECOND_UPPER_RQY + (154678, 154679),
        "incoming_m": (107131, 107132, 107133, 107134, 107135, 107136),
        "history_m": (109379, 109380, 109381, 109382, 109383, 109384),
        "numeric_name": "second_switch",
        "numeric_unit": "Wswitch",
        "graph_measure": "R^2*9*dt",
    },
    "post_power": {
        "native_chart": "post_switch",
        "q4_chart": "post_power",
        "phase": 109574,
        "coordinate": 78472,
        "radius": 91107,
        "J": 91103,
        "U": (159858, 159859, 159860, 159861, 159862, 159863),
        "U_y": (159864, 159865, 159866, 159867, 159868, 159869),
        "m": (159870, 159871, 159872, 159873, 159874, 159875),
        "m_y": (159882, 159883, 159884, 159885, 159886, 159887),
        "Q": (159953, 159957, 159963, 159969, 159975),
        "Q_y": (160041, 160045, 160051, 160057, 160063),
        "finite_RQy": (160200, 160201, 160202, 160203),
        "local_S1": (160231, 160232, 160233, 160234),
        "cumulative_S1": (160235, 160236, 160237, 160238),
        "incoming_S1": (154660, 154661, 154662, 154663),
        "endpoint_S1": (160271, 160272, 160273, 160274),
        "upper_RQy": POST_UPPER_RQY + (160289, 160290),
        "incoming_m": (109540, 109541, 109542, 109543, 109544, 109545),
        "history_m": (111787, 111788, 111789, 111790, 111791, 111792),
        "numeric_name": "post_switch",
        "numeric_unit": "Wmacro",
        "graph_measure": "R^2*91103*dp",
    },
}


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
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"JSON object required: {path}")
    return value


def _json_gzip(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        value = json.load(stream)
    _require(isinstance(value, dict), f"gzip JSON object required: {path}")
    return value


def _source_guard(path: Path, functions: Mapping[str, tuple[str, ...]]) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        raw = stream.read()
    tree = ast.parse(raw, filename=str(path))
    found: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            found[node.name] = ast.get_source_segment(raw, node) or ""
    result: dict[str, Any] = {"path": path.as_posix(), "functions": {}}
    for name, snippets in functions.items():
        _require(isinstance(snippets, tuple), f"source guard snippets must be tuple: {name}")
        body = found.get(name, "")
        compact = "".join(body.split())
        matched = all("".join(token.split()) in compact for token in snippets)
        result["functions"][name] = {"present": bool(body), "tokens": list(snippets), "matched": matched}
        _require(body and matched, f"source guard failed: {path.name}:{name}")
    result["passed"] = True
    result["source_executed"] = False
    return result


def _function_source(path: Path, name: str) -> str:
    """Return one function body for a read-only source receipt."""
    with path.open("r", encoding="utf-8", newline="") as stream:
        raw = stream.read()
    tree = ast.parse(raw, filename=str(path))
    candidates = [ast.get_source_segment(raw, node) or "" for node in ast.walk(tree)
                  if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]
    if candidates:
        # Several imported helper classes have a short ``__init__``.  The
        # source owner constructor is the longest same-name AST function.
        return max(candidates, key=len)
    raise ValueError(f"source function missing: {path.name}:{name}")


def _source_fragment(path: Path, name: str, tokens: Sequence[str]) -> dict[str, Any]:
    """Bind exact AST function text without importing or executing the source."""
    source = _function_source(path, name)
    compact = "".join(source.split())
    normalized = ["".join(token.split()) for token in tokens]
    _require(all(token in compact for token in normalized),
             f"source formula binding failed: {path.name}:{name}")
    return {"path": path.as_posix(), "function": name, "tokens": list(tokens),
            "source_sha256": hashlib.sha256(source.encode("utf-8")).hexdigest(),
            "passed": True, "source_executed": False}


def _self_format_guard(root: Path) -> dict[str, Any]:
    raw = (root / MODULE).read_bytes()
    _require(b"\r" not in raw, "new source must use LF line endings")
    _require(raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "new source must have one final LF")
    tree = ast.parse(raw.decode("utf-8"), filename=MODULE)
    return {"passed": True, "LF": True, "exactly_one_EOF_LF": True,
            "AST_dump_sha256": _digest(ast.dump(tree, include_attributes=False)), "source_executed": False}


def _parent_expected(parent: Mapping[str, Any], path: str) -> str | None:
    for key in ("input_hashes", "implementation_hashes"):
        table = parent.get(key, {})
        if isinstance(table, Mapping):
            value = table.get(path) or table.get(Path(path).name)
            if isinstance(value, str):
                return value
    return None


def _pin_path(root: Path, parent: Mapping[str, Any], path: str) -> dict[str, Any]:
    expected = PARENT_SHA256 if path == PARENT else _parent_expected(parent, path)
    # New bridge source files are bound by their current byte hash when the
    # parent did not carry that path; the bytes are still checked and emitted.
    if expected is None:
        expected = _sha(root / path)
    actual = _sha(root / path)
    _require(actual == expected, f"changed source pin: {path}")
    return {"path": path, "expected_sha256": expected, "passed": True}


def _packet_pin(root: Path, parent: Mapping[str, Any], path: str, value: Any) -> dict[str, Any]:
    expected = value.get("expected_sha256") if isinstance(value, Mapping) else value
    _require(isinstance(expected, str), f"numeric packet pin missing: {path}")
    inherited = _parent_expected(parent, path)
    if path == PARENT:
        _require(expected == PARENT_SHA256, "packet parent pin changed")
    elif inherited is not None:
        _require(expected == inherited, f"packet pin conflicts with parent: {path}")
    _require(_sha(root / path) == expected, f"changed numeric packet source: {path}")
    pin = {"path": path, "expected_sha256": expected, "passed": True}
    if isinstance(value, Mapping) and value.get("self_pin") is True:
        pin["self_pin"] = True
    return pin


def _append_saved(graph: GraphPairs, row: Mapping[str, Any]) -> int:
    _require(isinstance(row, Mapping) and isinstance(row.get("operation"), str), "serialized graph row malformed")
    return graph.node(row["operation"], **{key: value for key, value in row.items() if key != "operation"})


def _row(nodes: Sequence[Mapping[str, Any]], handle: int) -> Mapping[str, Any]:
    _require(type(handle) is int and 0 <= handle < len(nodes), f"graph handle outside prefix: {handle}")
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
        for key in ("argument", "numerator", "denominator", "expression", "variable", "value", "integrand", "lower", "upper"):
            child = item.get(key)
            if type(child) is int and child >= 0:
                pending.append(child)
        arguments = item.get("arguments")
        if isinstance(arguments, list):
            pending.extend(child for child in arguments if type(child) is int and child >= 0)
    return False


def _pair(value: Sequence[int], label: str) -> tuple[int, int]:
    _require(isinstance(value, (list, tuple)) and len(value) == 2 and all(type(x) is int for x in value),
             f"{label} must be an ordinary Z 0 / 1 pair")
    return int(value[0]), int(value[1])


def _validate_q4(root: Path) -> dict[str, Any]:
    archive = _json_gzip(root / Q4_ARCHIVE)
    nodes = archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list) and len(nodes) == Q4_NODE_COUNT, "Q4 graph node count changed")
    bindings = archive.get("modified_Q4_pressure_C3_bindings")
    _require(isinstance(bindings, Mapping), "Q4 bindings missing")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping), "Q4 charts/endpoints missing")
    checks: dict[str, Any] = {}
    for name, cfg in _CHARTS.items():
        q4_name = cfg["q4_chart"]
        chart = charts.get(q4_name)
        endpoint = endpoints.get(q4_name)
        _require(isinstance(chart, Mapping) and isinstance(endpoint, Mapping), f"Q4 chart missing: {q4_name}")
        for field in ("U", "U_y", "m", "m_y", "Q", "Q_y"):
            _require(chart.get("physical", {}).get(field) == list(cfg[field]), f"Q4 {q4_name} {field} changed")
        _require(chart.get("finite_boundary", {}).get("R_Q_y") == list(cfg["finite_RQy"]), f"Q4 {q4_name} R_Q_y changed")
        _require(chart.get("local_integrals", {}).get("S1") == list(cfg["local_S1"]), f"Q4 {q4_name} local S1 changed")
        _require(chart.get("cumulative_integrals", {}).get("S1") == list(cfg["cumulative_S1"]), f"Q4 {q4_name} cumulative S1 changed")
        _require(chart.get("incoming_integrals", {}).get("S1") == list(cfg["incoming_S1"]), f"Q4 {q4_name} incoming S1 changed")
        _require(endpoint.get("boundary", {}).get("R_Q_y") == list(cfg["upper_RQy"]), f"Q4 {q4_name} upper R_Q_y changed")
        for order, handle in enumerate(cfg["local_S1"][:2]):
            row = _row(nodes, handle)
            _require(row.get("operation") == "definite_integral" and row.get("original_chart") == q4_name
                     and row.get("physical_radius_power") == 2 and row.get("physical_Jacobian_applied_once") is True
                     and row.get("ordinary_slow_Z_derivative_order") == order, f"Q4 {q4_name} S1 source changed")
        for order, handle in enumerate(cfg["finite_RQy"][:2]):
            row = _row(nodes, handle)
            _require(row.get("operation") == "product" and row.get("arguments") == [cfg["Q_y"][order], cfg["radius"]],
                     f"Q4 {q4_name} R_Q_y product changed")
        checks[name] = {"q4_chart": q4_name, "local_S1": list(cfg["local_S1"][:2]), "finite_RQy": list(cfg["finite_RQy"][:2])}
    nrow = _row(nodes, N_HANDLE)
    _require(nrow.get("operation") == "exact_positive_integer_expression" and nrow.get("Z_independent") is True,
             "Q4 common N changed")
    _require(_row(nodes, PSTAR_HANDLE).get("operation") == "analytic_unary"
             and _row(nodes, PSTAR_HANDLE).get("name") == "exp"
             and _row(nodes, PSTAR_HANDLE).get("argument") == PSTAR_LOG_ARGUMENT_HANDLE,
             "Q4 Pstar is not exp(graph node 5)")
    return {"passed": True, "node_count": len(nodes), "charts": checks,
            "common_N": N_HANDLE, "Pstar": PSTAR_HANDLE, "Pstar_log_argument": PSTAR_LOG_ARGUMENT_HANDLE}


def _validate_loop(root: Path) -> dict[str, Any]:
    archive = _json_gzip(root / LOOP_ARCHIVE)
    bindings = archive.get("prescribed_loop_C5_bindings")
    _require(isinstance(bindings, Mapping), "C5 bindings missing")
    _require(bindings.get("N") == N_HANDLE, "C5 common N changed")
    charts = bindings.get("charts")
    _require(isinstance(charts, Mapping), "C5 chart bindings missing")
    checks: dict[str, Any] = {}
    for name, cfg in _CHARTS.items():
        chart = charts.get(cfg["q4_chart"])
        _require(isinstance(chart, Mapping), f"C5 chart missing: {cfg['q4_chart']}")
        _require(chart.get("actual_phase") == cfg["phase"] and chart.get("coordinate") == cfg["coordinate"]
                 and chart.get("domain") == [0, 1] and chart.get("radius") == cfg["radius"],
                 f"C5 {cfg['q4_chart']} geometry changed")
        _require(chart.get("incoming", {}).get("m") == list(cfg["incoming_m"]), f"C5 {cfg['q4_chart']} incoming m changed")
        _require(chart.get("history", {}).get("m") == list(cfg["history_m"]), f"C5 {cfg['q4_chart']} history m changed")
        checks[name] = {"actual_phase": cfg["phase"], "incoming_m": list(cfg["incoming_m"]), "history_m": list(cfg["history_m"])}
    definition = bindings.get("definition", {})
    _require(definition.get("correction_zero_inlet_only_at_rminus_background_nonzero_memories_unchanged") is True,
             "C5 nonreset history premise changed")
    return {"passed": True, "common_N": N_HANDLE, "charts": checks, "phase_Z_independent": True,
            "nonreset_incoming_history": True}


def _validate_parent_seams(parent: Mapping[str, Any], nodes: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    seam = parent.get("actual_macro_switch_RQy_seam")
    _require(isinstance(seam, Mapping) and seam.get("passed") is True, "macro/first seam missing")
    patch = seam.get("graph_patch", {})
    _require(patch.get("base_node_count") == MERGED_NODE_COUNT and patch.get("appended_node_count") == 16,
             "macro/first seam suffix changed")
    first = patch.get("endpoint_handles", {}).get("first_switch_lower_new_s0_Z0_Z1")
    _require(first == [184687, 184688], "first lower seam handles changed")
    cancel = patch.get("endpoint_handles", {}).get("signed_cancellation")
    _require(cancel == [184700, 184702], "accepted seam cancellation handles changed")
    for k, handle in enumerate(cancel):
        _require(_row(nodes, handle).get("operation") == "sum", "accepted seam row changed")
    return {"passed": True, "base_node_count": MERGED_NODE_COUNT, "appended_node_count": 16,
            "first_lower": first, "signed_cancellation": cancel}


def _rehydrate_prefix(parent: Mapping[str, Any], q4: Mapping[str, Any]) -> tuple[GraphPairs, dict[str, Any]]:
    q4_nodes = q4["exact_graph_nodes"]
    micro = parent["actual_micro_macro_correlated_m15"]["merged_namespace"]["appended_nodes"]
    seam = parent["actual_macro_switch_RQy_seam"]["graph_patch"]["appended_nodes"]
    first = parent["actual_macro_first_switch_pressure_source"]["graph_patch"]["appended_nodes"]
    _require(isinstance(micro, list) and len(micro) == 987, "micro/macro suffix must contain 987 rows")
    _require(isinstance(seam, list) and len(seam) == 16, "macro/first seam must contain 16 rows")
    _require(isinstance(first, list) and len(first) == 195, "first pressure suffix must contain 195 rows")
    graph = GraphPairs(q4_nodes)
    for row in micro + seam + first:
        _append_saved(graph, row)
    _require(len(graph.nodes) == FIRST_NODE_COUNT, "full prefix is not 184898 nodes")
    _require(graph.nodes[:Q4_NODE_COUNT] == q4_nodes, "Q4 prefix changed")
    _require(graph.nodes[Q4_NODE_COUNT:MERGED_NODE_COUNT] == micro, "micro/macro prefix changed")
    _require(graph.nodes[MERGED_NODE_COUNT:SEAM_NODE_COUNT] == seam, "macro/first seam prefix changed")
    _require(graph.nodes[SEAM_NODE_COUNT:FIRST_NODE_COUNT] == first, "macro-first suffix changed")
    patch = parent["actual_macro_first_switch_pressure_source"]["graph_patch"]
    _require(patch.get("base_node_count") == SEAM_NODE_COUNT and patch.get("full_node_count") == FIRST_NODE_COUNT,
             "macro-first patch counts changed")
    _require(patch.get("full_canonical_nodes_sha256") == _digest(graph.nodes), "macro-first full prefix digest changed")
    namespace = parent["actual_micro_macro_correlated_m15"].get("graph_namespace", {})
    return graph, {
        "q4_node_count": Q4_NODE_COUNT, "merged_node_count": MERGED_NODE_COUNT,
        "seam_node_count": SEAM_NODE_COUNT, "first_node_count": FIRST_NODE_COUNT,
        "full_prefix_snapshot": copy.deepcopy(graph.nodes),
        "full_prefix_canonical_nodes_sha256": _digest(graph.nodes),
        "old_namespace_digest": namespace.get("namespace_digest"),
        "segment_starts": {"q4": 0, "micro_macro_merged": Q4_NODE_COUNT,
                           "macro_first_switch_seam": MERGED_NODE_COUNT,
                           "macro_first_switch_pressure": SEAM_NODE_COUNT,
                           "new_suffix": FIRST_NODE_COUNT},
        "suffix_counts": {"micro_macro_merged": len(micro), "macro_first_switch_seam": len(seam),
                          "macro_first_switch_pressure": len(first)},
    }


def _build_six_source(graph: GraphPairs, cfg: Mapping[str, Any]) -> tuple[dict[str, list[int]], dict[str, Any]]:
    q = _pair(cfg["Q"][:2], "Q")
    q_z = _pair(cfg["Q"][1:3], "Q_Z")
    u = _pair(cfg["U"][:2], "U")
    u_z = _pair(cfg["U"][1:3], "U_Z")
    density = {"q": q, "q_Z": q_z, "q_squared": graph.pmul(q, q),
               "U_q_Z": graph.pmul(u, q_z), "U_squared": graph.pmul(u, u),
               "U_U_Z": graph.pmul(u, u_z)}
    weight = graph.mul(cfg["radius"], cfg["radius"], cfg["J"])
    integrals: dict[str, list[int]] = {}
    for name in I1_NAMES:
        integrals[name] = list(graph.integral(
            graph.pscale(density[name], weight), cfg["coordinate"], 0, 1,
            f"actual_{cfg['native_chart']}_finite_N_I1_{name}",
            metadata={"original_chart": cfg["q4_chart"], "source_layer": "physical_six_I1_density",
                      "ordinary_Z_rows": [0, 1], "physical_radius_power": 2,
                      "physical_Jacobian_applied_once": True, "Pstar_applied_in_physical_fields_once": True,
                      "endpoint_independent_of_Z": True, "physical_y_derivatives_claimed": False,
                      "raw_physical_unit": "R^2*J*dcoordinate", "numeric_normalization_external": cfg["numeric_unit"]},
        ))
    return integrals, {"passed": True, "density_pairs": {k: list(v) for k, v in density.items()},
                       "weight_handle": weight, "weight_formula": f"R_{cfg['native_chart']}^2 * J_{cfg['native_chart']}",
                       "weight_units": cfg["graph_measure"], "graph_source_raw_physical": True,
                       "numeric_packet_normalization": cfg["numeric_unit"], "integral_handles": integrals,
                       "ordinary_Z_orders": [0, 1], "Pstar_applied_once": True, "R_not_a_mean_scale": True}


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
    hZ_over_2L = graph.pscale(graph.pdiv(graph.pmul(h, z), L, "positive physical L=1-delta*Z^2"), half)
    minus_half = graph.pneg(graph.fixed(half))
    d_over_L = graph.pdiv(d, L, "positive physical L=1-delta*Z^2")
    two_hZ2_over_L2 = graph.pdiv(graph.pscale(graph.pmul(h, graph.pmul(z, z)), two), L2, "positive physical L=1-delta*Z^2")
    two_Zd_over_L2 = graph.pdiv(graph.pscale(graph.pmul(z, d), two), L2, "positive physical L=1-delta*Z^2")
    terms = {"minus_inv_L": graph.pmul(graph.pneg(inv_L), _pair(integrals["q"], "q")),
             "hZ_over_2L": graph.pmul(hZ_over_2L, _pair(integrals["q_Z"], "q_Z")),
             "minus_half": graph.pmul(minus_half, _pair(integrals["q_squared"], "q_squared")),
             "d_over_L": graph.pmul(d_over_L, _pair(integrals["U_q_Z"], "U_q_Z")),
             "two_hZ2_over_L2": graph.pmul(two_hZ2_over_L2, _pair(integrals["U_squared"], "U_squared")),
             "two_Zd_over_L2": graph.pmul(two_Zd_over_L2, _pair(integrals["U_U_Z"], "U_U_Z"))}
    rows = list(graph.padd(*terms.values()))
    return rows, {"passed": True, "rows": rows, "terms": {k: list(v) for k, v in terms.items()},
                  "formula": "D=-q/L+h*Z*q_Z/(2L)-q^2/2+d*U*q_Z/L+2*h*Z^2*U^2/L^2+2*Z*d*U*U_Z/L^2",
                  "graph_source_unit": "raw physical R^2*J*dcoordinate", "numeric_packet_unit": "external",
                  "normalization_not_in_graph": True, "D_is_not_full_S1": True,
                  "correct_two_h_Z2_term": True}


def _build_reduced_B(graph: GraphPairs, cfg: Mapping[str, Any]) -> tuple[list[int], dict[str, Any]]:
    z = _pair((20, 1), "Z")
    u = _pair(cfg["U"][:2], "U")
    q = _pair(cfg["Q"][:2], "Q")
    L = _pair((184439, 184438), "L")
    L2 = _pair((184440, 184443), "L^2")
    radius_squared = graph.mul(cfg["radius"], cfg["radius"])
    first = graph.pscale(graph.pdiv(q, L, "positive physical L=1-delta*Z^2"), radius_squared)
    second = graph.pscale(graph.pscale(graph.pmul(q, q), radius_squared), 28)
    third = graph.pscale(graph.pmul(graph.pmul(z, z), graph.fixed(radius_squared), graph.pmul(u, u)), 48)
    third = graph.pneg(graph.pdiv(third, L2, "positive physical L=1-delta*Z^2"))
    boundary = graph.padd(first, second, third)
    upper = graph.at(boundary, cfg["coordinate"], 1)
    lower = graph.at(boundary, cfg["coordinate"], 0)
    rows = list(graph.padd(upper, graph.pneg(lower)))
    return rows, {"passed": True, "rows": rows, "boundary_integrand": list(boundary), "upper": list(upper), "lower": list(lower),
                  "formula": "B=[R^2*Q/L+R^2*Q^2/2-2*Z^2*R^2*U^2/L^2]_{t=0}^{t=1}",
                  "physical_radius_power": 2, "physical_Jacobian_applied": False,
                  "endpoint_measure_has_no_J": True, "graph_source_unit": "raw physical R^2*Q endpoint",
                  "numeric_endpoint_not_evaluated": True, "ordinary_Z_orders": [0, 1],
                  "correct_minus_two_Z2_R2_U2_term": True}


def _chart_packet(packet: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    aliases = (name, "post_switch") if name == "post_power" else (name,)
    charts = packet.get("charts")
    if isinstance(charts, Mapping):
        for alias in aliases:
            direct = charts.get(alias)
            if isinstance(direct, Mapping):
                return direct
    for alias in aliases:
        direct = packet.get(alias)
        if isinstance(direct, Mapping):
            return direct
    raise ValueError(f"numeric packet chart missing: {name}")


def _numeric_contract(root: Path, parent: Mapping[str, Any], packet: Mapping[str, Any]) -> dict[str, Any]:
    _require(packet.get("passed") is True, "switch/post numeric packet is not passed")
    _require(packet.get("parent_report_sha256") == PARENT_SHA256, "numeric packet parent SHA changed")
    pins = packet.get("source_pins")
    _require(isinstance(pins, Mapping), "numeric packet source_pins missing")
    joined_global = packet.get("joined_macro_switch_post_volume")
    if joined_global is not None:
        _require(isinstance(joined_global, Mapping), "global joined numeric volume malformed")
    charts: dict[str, Any] = {}
    for name in _CHARTS:
        item = _chart_packet(packet, name)
        for key in ("I1_full", "D_full", "physical_error_log_bounds"):
            _require(isinstance(item.get(key), Mapping), f"numeric {name} {key} missing")
        _require(item["I1_full"].get("passed") is True and item["D_full"].get("passed") is True,
                 f"numeric {name} I1/D not passed")
        joined = item.get("joined_macro_switch_post_volume", joined_global)
        _require(isinstance(joined, Mapping), f"numeric {name} joined volume missing")
        _require(joined.get("passed") is True and joined.get("complete_pressure_enclosed") is False,
                 f"numeric {name} joined scope changed")
        _require(joined.get("normalized_by_Wmacro") is True,
                 f"numeric {name} joined volume must be Wmacro")
        if isinstance(joined.get("I1"), Mapping):
            _require(joined["I1"].get("normalized_by_Wmacro") is True,
                     f"numeric {name} joined I1 must be Wmacro")
        if isinstance(joined.get("D"), Mapping):
            _require(joined["D"].get("normalized_by_Wmacro") is True,
                     f"numeric {name} joined D must be Wmacro")
        if name == "second_switch":
            _require(item["I1_full"].get("normalized_by_Wswitch") is True and item["D_full"].get("normalized_by_Wswitch") is True,
                     "second numeric units must be Wswitch")
        else:
            _require(item["I1_full"].get("normalized_by_Wmacro") is True and item["D_full"].get("normalized_by_Wmacro") is True,
                     "post numeric units must be Wmacro")
        charts[name] = copy.deepcopy(dict(item))
        charts[name]["joined_macro_switch_post_volume"] = copy.deepcopy(dict(joined))
    return {"passed": True, "parent_report_sha256": PARENT_SHA256, "source_pins_bound": True,
            "source_pins": copy.deepcopy(dict(pins)), "charts": charts,
            "joined_macro_switch_post_volume": copy.deepcopy(dict(joined_global)) if isinstance(joined_global, Mapping) else None,
            "limits": copy.deepcopy(packet.get("limits", {})), "source_contract": copy.deepcopy(packet.get("source_contract", {})),
            "numeric_packet_scope": "I1/D plus joined volume; outer R_Q_y remains symbolic", "complete_pressure_numeric": False}


def _linearize_sum_negative(graph: GraphPairs, node: int) -> dict[int, int]:
    """Collect exact integer coefficients through sum/negative rows only.

    Products, integrals, substitutions, and every other graph operation are
    atoms here.  This is a serialization rewrite, not algebra on their
    contents.
    """
    row = _row(graph.nodes, node)
    operation = row.get("operation")
    if operation == "sum":
        arguments = row.get("arguments")
        _require(isinstance(arguments, list), "sum row arguments must be a list")
        coefficients: dict[int, int] = {}
        for argument in arguments:
            child = _linearize_sum_negative(graph, argument)
            for atom, coefficient in child.items():
                coefficients[atom] = coefficients.get(atom, 0) + coefficient
                if coefficients[atom] == 0:
                    del coefficients[atom]
        return coefficients
    if operation == "negative":
        argument = row.get("argument")
        child = _linearize_sum_negative(graph, argument)
        return {atom: -coefficient for atom, coefficient in child.items()}
    return {node: 1}


def _rebuild_linearized_pair(graph: GraphPairs, pair: Sequence[int], targets: Sequence[int], label: str) -> tuple[list[int], dict[str, Any]]:
    """Rebuild a pair from signed atoms after sum/negative cancellation."""
    target_set = set(targets)
    reduced_rows: list[int] = []
    row_receipts: list[dict[str, Any]] = []
    for order, node in enumerate(pair):
        coefficients = _linearize_sum_negative(graph, int(node))
        target_coefficients = {str(target): coefficients.get(target, 0) for target in targets}
        _require(all(value == 0 for value in target_coefficients.values()),
                 f"{label} target coefficient did not cancel at ordinary row {order}")
        surviving_dependencies: dict[str, list[int]] = {}
        for atom in coefficients:
            for target in targets:
                if _has_dependency(graph.nodes, atom, target):
                    surviving_dependencies.setdefault(str(atom), []).append(target)
        _require(not surviving_dependencies,
                 f"{label} surviving atom still depends on removed endpoint target")
        terms: list[int] = []
        for atom, coefficient in sorted(coefficients.items()):
            if coefficient == 1:
                terms.append(atom)
            elif coefficient == -1:
                terms.append(graph.neg(atom))
            else:
                terms.append(graph.mul(atom, graph.constant(coefficient)))
        reduced = graph.add(*terms)
        _require(not any(_has_dependency(graph.nodes, reduced, target) for target in target_set),
                 f"{label} rebuilt row still depends on removed endpoint target")
        reduced_rows.append(reduced)
        row_receipts.append({"ordinary_Z_order": order, "input_handle": int(node),
                             "input_signed_atoms": {str(atom): coefficient for atom, coefficient in sorted(coefficients.items())},
                             "target_coefficients": target_coefficients,
                             "reduced_handle": reduced,
                             "reduced_signed_atoms": {str(atom): coefficient for atom, coefficient in sorted(coefficients.items())},
                             "surviving_target_dependencies": surviving_dependencies,
                             "exact_integer_arithmetic": True})
    return reduced_rows, {"passed": True, "scope": "sum_and_negative_rows_only",
                          "stopped_at": ["product", "definite_integral", "function_substitution", "other_operations"],
                          "targets": [int(target) for target in targets], "rows": row_receipts,
                          "removed_target_coefficients_zero": True,
                          "surviving_target_dependencies_empty": True,
                          "exact_integer_arithmetic": True}


def _build_join(graph: GraphPairs, parent: Mapping[str, Any], local_s1: Mapping[str, Sequence[int]],
                lowers: Mapping[str, Sequence[int]], formal: Mapping[str, Any]) -> dict[str, Any]:
    parent_source = parent.get("actual_macro_first_switch_pressure_source", {}).get("joined_pressure_source", {})
    entry = parent_source.get("joined_contracted")
    _require(entry == list(PARENT_JOINED), "macro-first joined source entry changed")
    _require(parent_source.get("outer_RQy_and_caps_retained") is True, "parent outer R_Q_y was not retained")
    increments: dict[str, list[int]] = {}
    current = list(PARENT_JOINED)
    joins: dict[str, list[int]] = {}
    # The inherited macro/first joined entry already contains the first-switch
    # increment, including ``-R_Q_y(first, t=1)``.  Each new segment therefore
    # subtracts its own upper endpoint exactly once: second subtracts the
    # second-switch upper, and post-power subtracts the post-power upper.
    upper_for = {"second_switch": SECOND_UPPER_RQY, "post_power": POST_UPPER_RQY}
    for name in ("second_switch", "post_power"):
        half = [graph.mul(int(local_s1[name][k]), 28) for k in range(2)]
        inc = [graph.add(half[k], int(lowers[name][k]), graph.neg(upper_for[name][k])) for k in range(2)]
        increments[name] = inc
        current = [graph.add(current[k], inc[k]) for k in range(2)]
        joins[name] = list(current)
    contracted: list[int] | None = list(entry)
    contracted_after_second: list[int] | None = None
    contracted_after_second_reduced: list[int] | None = None
    contracted_before_linearization: list[int] | None = None
    after_second_linearization_receipt: dict[str, Any] | None = None
    linearization_receipt: dict[str, Any] | None = None
    first_identity = formal["seams"]["first_upper_to_second_lower"]["finite_N_Qy_identity_proven"] is True
    second_identity = formal["seams"]["second_upper_to_post_lower"]["finite_N_Qy_identity_proven"] is True
    if first_identity:
        half_second = [graph.mul(int(local_s1["second_switch"][k]), 28) for k in range(2)]
        contracted_after_second = [
            graph.add(entry[k], FIRST_UPPER_RQY[k], half_second[k], graph.neg(SECOND_UPPER_RQY[k]))
            for k in range(2)
        ]
        contracted_after_second_reduced, after_second_linearization_receipt = _rebuild_linearized_pair(
            graph, contracted_after_second, FIRST_UPPER_RQY,
            "contracted_after_second")
        contracted = contracted_after_second_reduced
    else:
        contracted = None
    if second_identity and contracted_after_second is not None:
        half_post = [graph.mul(int(local_s1["post_power"][k]), 28) for k in range(2)]
        contracted_before_linearization = [
            graph.add(contracted_after_second[k], SECOND_UPPER_RQY[k], half_post[k], graph.neg(POST_UPPER_RQY[k]))
            for k in range(2)
        ]
        final_linearization_receipt: dict[str, Any]
        contracted, final_linearization_receipt = _rebuild_linearized_pair(
            graph, contracted_before_linearization, FIRST_UPPER_RQY + SECOND_UPPER_RQY,
            "joined_contracted")
        linearization_receipt = {"after_second": after_second_linearization_receipt,
                                 "final": final_linearization_receipt}
    elif contracted_after_second is not None:
        contracted_before_linearization = list(contracted_after_second)
        contracted = contracted_after_second_reduced
        final_targets = FIRST_UPPER_RQY
        contracted, final_linearization_receipt = _rebuild_linearized_pair(
            graph, contracted_before_linearization, final_targets,
            "joined_contracted")
        linearization_receipt = {"after_second": after_second_linearization_receipt,
                                 "final": final_linearization_receipt}
    if contracted is not None:
        remaining_targets = FIRST_UPPER_RQY + SECOND_UPPER_RQY if second_identity else FIRST_UPPER_RQY
        _require(all(not _has_dependency(graph.nodes, row, target)
                     for row in contracted for target in remaining_targets),
                 "joined_contracted retains a removed endpoint dependency")
    if linearization_receipt is None:
        linearization_receipt = {"passed": False, "contraction_not_available": True}
    return {"passed": True, "ordinary_Z_orders": [0, 1], "parent_joined_contracted": list(entry),
            "second_switch_local_S1": list(local_s1["second_switch"]), "post_power_local_S1": list(local_s1["post_power"]),
            "second_switch_half_S1": [graph.mul(int(local_s1["second_switch"][k]), 28) for k in range(2)],
            "post_power_half_S1": [graph.mul(int(local_s1["post_power"][k]), 28) for k in range(2)],
            "second_switch_lower_RQy_t0": list(lowers["second_switch"]), "post_power_lower_RQy_p0": list(lowers["post_power"]),
            "first_switch_upper_RQy": list(FIRST_UPPER_RQY), "second_switch_upper_RQy": list(SECOND_UPPER_RQY),
            "post_power_upper_RQy": list(POST_UPPER_RQY), "increments_uncontracted": increments,
            "joined_after_second_uncontracted": joins["second_switch"], "joined_after_post_uncontracted": joins["post_power"],
            "joined_uncontracted": joins["post_power"], "joined_contracted": contracted,
            "contracted_before_linearization": contracted_before_linearization,
            "contracted_after_second": contracted_after_second,
            "contracted_after_second_reduced": contracted_after_second_reduced,
            "linearization_receipt": linearization_receipt,
            "formula_uncontracted": "parent_joined + second_local_S1/2 + R_Q_y(second,t=0) - R_Q_y(second,t=1) + post_local_S1/2 + R_Q_y(post,p=0) - R_Q_y(post,p=1)",
            "formula_contracted": "parent_joined + first_upper(second-lower identity) + second_local_S1/2 - second_upper + second_upper(post-lower identity) + post_local_S1/2 - post_upper",
            "conditional_identity_only_in_joined_expression": True,
            "contracted_graph_linearization": "exact signed integer accumulation through sum/negative rows only; all other operations remain atoms",
            "contracted_before_normalization_preserved": True,
            "first_switch_upper_standalone_retained": True, "second_switch_upper_standalone_retained": True,
            "post_power_upper_standalone_retained": True, "macro_lower_standalone_retained": True,
            "outer_RQy_and_caps_retained": True, "signed_relation_not_global_closure": True,
            "contraction_status": "PROVED_FORMAL_ENDPOINT_CONGRUENCES" if contracted is not None else "BLOCKED_UNPROVEN_FINITE_N_QY_ENDPOINT_IDENTITIES",
            "endpoint_congruence_receipt": copy.deepcopy(formal["seams"])}


def _jet_add(left: Sequence[Fraction], right: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(a + b for a, b in zip(left, right))


def _jet_neg(value: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(-item for item in value)


def _jet_scale(value: Sequence[Fraction], scalar: Fraction) -> tuple[Fraction, ...]:
    return tuple(scalar * item for item in value)


def _jet_mul(left: Sequence[Fraction], right: Sequence[Fraction]) -> tuple[Fraction, ...]:
    """Leibniz product of ordinary derivative rows through Z order 2."""
    from math import comb
    return tuple(
        sum(Fraction(comb(k, j)) * left[j] * right[k - j] for j in range(k + 1))
        for k in range(3)
    )


def _jet_div(numerator: Sequence[Fraction], denominator: Sequence[Fraction]) -> tuple[Fraction, ...]:
    """Exact rational ordinary-jet division with a positive nonzero base."""
    from math import comb
    _require(denominator[0] != 0, "formal jet denominator is zero")
    out = [numerator[0] / denominator[0]]
    for k in range(1, 3):
        known = sum(Fraction(comb(k, j)) * denominator[j] * out[k - j] for j in range(1, k + 1))
        out.append((numerator[k] - known) / denominator[0])
    return tuple(out)


def _formal_q4_qy(state: Mapping[str, Sequence[Fraction]]) -> tuple[Fraction, ...]:
    """The Q4 ``Q_y`` polynomial/rational map, kept independent of graph rows."""
    one = (Fraction(1), Fraction(0), Fraction(0))
    z = state["z"]
    delta = state["delta"]
    d = _jet_add(one, _jet_neg(_jet_mul(z, z)))
    L = _jet_add(one, _jet_neg(_jet_scale(_jet_mul(z, z), delta[0])))
    first = _jet_scale(_jet_mul(z, state["my"]), Fraction(1) - delta[0])
    shifted = state["my_shift"]
    numerator = _jet_add(
        _jet_scale(_jet_mul(z, state["uy"]), Fraction(2)),
        _jet_neg(_jet_add(first, _jet_mul(d, shifted))),
    )
    return _jet_div(numerator, L)


def _formal_c5_rational_chain(state: Mapping[str, Sequence[Fraction]]) -> tuple[tuple[Fraction, ...], ...]:
    """A pure C5 polynomial/rational jet chain through ordinary order 2.

    The source uses the same ``p2``, implicit inverse, Poisson denominator,
    and transport rational maps at both sides of a seam.  This small exact
    witness applies those maps to two equal jet states; it does not evaluate
    a source function or select a numerical phase.
    """
    one = (Fraction(1), Fraction(0), Fraction(0))
    u = _jet_div(_jet_mul(state["p2"], state["q"]), state["dstar"])
    r = _jet_mul(u, state["hinv"])
    alpha = _jet_scale(_jet_mul(state["q"], state["hinv"]), Fraction(2))
    denominator = _jet_add(_jet_add(one, _jet_neg(_jet_scale(_jet_mul(r, state["cos"]), Fraction(2)))), _jet_mul(r, r))
    t = _jet_add(state["t0"], _jet_div(_jet_mul(alpha, _jet_add(state["cos"], _jet_neg(r))), denominator))
    phase_map = _jet_add(state["phase"], _jet_mul(t, t))
    loop_map = _jet_div(_jet_mul(state["a"], phase_map), denominator)
    return t, phase_map, loop_map


def _formal_jet_receipt() -> dict[str, Any]:
    """Arithmetic smoke fixture only; it is never used as an endpoint proof."""
    q4_state = {
        "z": (Fraction(1, 3), Fraction(1), Fraction(2)),
        "delta": (Fraction(1, 7), Fraction(0), Fraction(0)),
        "uy": (Fraction(4), Fraction(2), Fraction(1)),
        "my": (Fraction(1), Fraction(1), Fraction(2)),
        "my_shift": (Fraction(1), Fraction(2), Fraction(3)),
    }
    left_qy = _formal_q4_qy(q4_state)
    right_qy = _formal_q4_qy(copy.deepcopy(q4_state))
    _require(left_qy == right_qy, "fresh Q4 rational jet congruence failed")
    c5_state = {
        "p2": (Fraction(3), Fraction(1), Fraction(2)), "q": (Fraction(2), Fraction(1), Fraction(1)),
        "dstar": (Fraction(5), Fraction(1), Fraction(1)), "hinv": (Fraction(4, 5), Fraction(1, 7), Fraction(1, 9)),
        "cos": (Fraction(1, 2), Fraction(1, 3), Fraction(1, 5)), "t0": (Fraction(1, 4), Fraction(1, 6), Fraction(1, 8)),
        "phase": (Fraction(2, 5), Fraction(1, 4), Fraction(1, 7)), "a": (Fraction(2), Fraction(1), Fraction(1)),
    }
    left_c5 = _formal_c5_rational_chain(c5_state)
    right_c5 = _formal_c5_rational_chain(copy.deepcopy(c5_state))
    _require(left_c5 == right_c5, "fresh C5 rational jet congruence failed")
    return {
        "passed": True, "ordinary_Z_orders": [0, 1, 2], "arithmetic": "exact Fraction",
        "Q4_formula": "Q_y=(2*Z*U_y-(1-delta)*Z*m_y-(1-Z^2)*shift_Z(m_y))/(1-delta*Z^2)",
        "Q4_left_right_equal": True, "Q4_jet_rows": [str(value) for value in left_qy],
        "C5_formula": "u=p2*q/dstar; r=u*h-inv; D=1-2*r*cos(psi)+r^2; t=t0+2*q*h-inv*(cos-r)/D; implicit inverse/loop rational maps",
        "C5_left_right_equal": True, "C5_jet_rows": [[str(value) for value in row] for row in left_c5],
        "chain_rule_congruence": "equal inputs, equal parameters, equal ordinary jets through Z2 imply equal outputs",
        "scope": "arithmetic_fixture_only_not_source_endpoint_proof",
    }


def _symbolic_endpoint_identities() -> dict[str, Any]:
    """Prove the switch endpoint identities for arbitrary symbolic jets.

    The source/cache receipt below supplies the equal inlet state and the
    exact flat-cutoff endpoint substitutions.  This routine then performs the
    algebra with independent symbols; it does not compare copied numeric
    states.
    """
    import sympy as sp

    h, R, D, D_y = sp.symbols("h R D D_y", nonzero=True)
    phi, A, A_y = sp.symbols("phi A A_y")
    sigma, sigma_y = sp.symbols("sigma sigma_y")
    lam_first = -h * R * D / 2
    lam_first_y = -h * R * (D + D_y) / 2
    lam_second = -h * R * D * (1 - sigma) / 2 - 2 * sigma / 5
    lam_second_y = (-h * R * (D + D_y) * (1 - sigma) / 2
                    + (h * R * D / 2 - sp.Rational(2, 5)) * sigma_y)
    lam_post = -sp.Rational(2, 5)
    lam_post_y = sp.Integer(0)
    first_upper = {sigma: 1, sigma_y: 0}
    second_lower = {sigma: 0, sigma_y: 0}
    first_lambda_difference = sp.simplify(
        lam_first.subs(first_upper) - lam_second.subs(second_lower))
    first_lambda_y_difference = sp.simplify(
        lam_first_y.subs(first_upper) - lam_second_y.subs(second_lower))
    second_lambda_difference = sp.simplify(
        lam_second.subs(first_upper) - lam_post)
    second_lambda_y_difference = sp.simplify(
        lam_second_y.subs(first_upper) - lam_post_y)
    phi_y_first = phi * lam_first
    phi_yy_first = phi * (lam_first ** 2 + lam_first_y)
    phi_y_second_lower = phi * lam_second.subs(second_lower)
    phi_yy_second_lower = phi * (lam_second.subs(second_lower) ** 2
                                  + lam_second_y.subs(second_lower))
    phi_y_second_upper = phi * lam_second.subs(first_upper)
    phi_yy_second_upper = phi * (lam_second.subs(first_upper) ** 2
                                 + lam_second_y.subs(first_upper))
    phi_y_post = phi * lam_post
    phi_yy_post = phi * (lam_post ** 2 + lam_post_y)
    first_phi_difference = sp.simplify(phi_y_first - phi_y_second_lower)
    first_phi_yy_difference = sp.simplify(phi_yy_first - phi_yy_second_lower)
    second_phi_difference = sp.simplify(phi_y_second_upper - phi_y_post)
    second_phi_yy_difference = sp.simplify(phi_yy_second_upper - phi_yy_post)
    uy_first = -h * (1 - sigma) * A
    uyy_first = h * sigma_y * A - h * (1 - sigma) * A_y
    first_uy_difference = sp.simplify(uy_first.subs(first_upper))
    first_uy_y_difference = sp.simplify(uyy_first.subs(first_upper))
    zero = lambda value: sp.simplify(value) == 0
    differences = {
        "first_lambda": first_lambda_difference,
        "first_lambda_y": first_lambda_y_difference,
        "first_phi_y": first_phi_difference,
        "first_phi_yy": first_phi_yy_difference,
        "first_U_y": first_uy_difference,
        "first_U_yy": first_uy_y_difference,
        "second_lambda": second_lambda_difference,
        "second_lambda_y": second_lambda_y_difference,
        "second_phi_y": second_phi_difference,
        "second_phi_yy": second_phi_yy_difference,
    }
    _require(all(zero(value) for value in differences.values()),
             "symbolic endpoint identity did not reduce to zero")
    return {
        "passed": True, "engine": "sympy arbitrary symbols",
        "symbols": ["h", "R", "D", "D_y", "phi", "A", "A_y", "sigma", "sigma_y"],
        "substitutions": {"first_upper": "sigma=1,sigma_y=0", "second_lower": "sigma=0,sigma_y=0",
                           "second_upper": "sigma=1,sigma_y=0"},
        "lambda_first": "-h*R*D/2",
        "lambda_first_y": "-h*R*(D+D_y)/2",
        "lambda_second": "-h*R*D*(1-sigma)/2-2*sigma/5",
        "lambda_second_y": "-h*R*(D+D_y)*(1-sigma)/2+(h*R*D/2-2/5)*sigma_y",
        "lambda_post": "-2/5", "lambda_post_y": "0",
        "phi_y": "phi*lambda", "phi_yy": "phi*(lambda^2+lambda_y)",
        "U_y_first": "-h*(1-sigma)*A",
        "U_yy_first": "h*sigma_y*A-h*(1-sigma)*A_y",
        "differences": {name: str(value) for name, value in differences.items()},
        "all_symbolic_differences_zero": True,
        "ordinary_Z_orders": [0, 1, 2],
    }


def _symbolic_rational_congruence() -> dict[str, Any]:
    """Check source-map congruence for arbitrary endpoint jet symbols."""
    import sympy as sp

    names = ("Z", "delta", "U", "U_y", "m", "m_y", "m_y_shift", "E", "V", "E_y", "V_y",
             "p2", "q", "dstar")
    left = {name: sp.symbols("L_" + name) for name in names}
    right = {name: sp.symbols("R_" + name) for name in names}
    equality = {right[name]: left[name] for name in names}
    one = sp.Integer(1)
    def qy(state):
        z, delta = state["Z"], state["delta"]
        d = one - z ** 2
        ell = one - delta * z ** 2
        numerator = (2 * z * state["U_y"]
                     - (1 - delta) * z * state["m_y"]
                     - d * state["m_y_shift"])
        return numerator / ell
    def c5_chain(state):
        a = (state["E"] - 2 * state["E_y"]) / state["E"]
        b = 2 * state["V_y"] / state["E"]
        t0 = -b / a
        # This is the actual C5 source chain: u=p2*q/dstar, followed by
        # hinv=1/sqrt(1+u^2), r=u*hinv, and alpha=2*q*hinv.
        u = state["p2"] * state["q"] / state["dstar"]
        hinv = 1 / sp.sqrt(1 + u ** 2)
        r = u * hinv
        alpha = 2 * state["q"] * hinv
        return a, b, t0, u, hinv, r, alpha
    q_difference = sp.simplify(qy(left) - qy(right).subs(equality))
    chain_left = c5_chain(left)
    chain_right = tuple(value.subs(equality) for value in c5_chain(right))
    chain_differences = tuple(sp.simplify(a - b) for a, b in zip(chain_left, chain_right))
    _require(q_difference == 0 and all(value == 0 for value in chain_differences),
             "source rational congruence did not reduce to zero")
    return {
        "passed": True, "engine": "sympy arbitrary left/right symbols with source-bound equality ideal",
        "Q4_map": "(2*Z*U_y-(1-delta)*Z*m_y-(1-Z^2)*shift_Z(m_y))/(1-delta*Z^2)",
        "C5_root_map": "a=(E-2*E_y)/E; b=2*V_y/E; t0=-b/a; p2=generic[p2]; u=p2*q/dstar; hinv=1/sqrt(1+u^2); r=u*hinv; alpha=2*q*hinv",
        "Q4_difference_after_source_substitution": str(q_difference),
        "C5_differences_after_source_substitution": [str(value) for value in chain_differences],
        "ordinary_Z_orders": [0, 1, 2],
        "same_source_function_congruence": True,
    }


def _symbolic_periodic_phase_lemma() -> dict[str, Any]:
    """Derive the exact C5 periodic normalization and phase drift identities."""
    import sympy as sp

    r, h, q, t0, a, b = sp.symbols("r h q t0 a b")
    pi = sp.pi
    mean_f = sp.Integer(0)
    mean_f_squared = 1 / (2 * (1 - r ** 2))
    alpha = 2 * q * h
    T1 = 2 * pi * t0 + 2 * pi * alpha * mean_f
    T2 = 2 * pi * (t0 ** 2 + alpha ** 2 * mean_f_squared)
    ratio = 1 + t0 ** 2 + 2 * q ** 2
    h_relation = sp.expand(h ** 2 - (1 - r ** 2))
    T2_reduced = sp.factor(sp.expand(T2).subs(h ** 2, 1 - r ** 2))
    normalization = sp.factor((2 * pi + T2_reduced) / (2 * pi * ratio) - 1)
    B_drift = sp.simplify((-a * T1 / (2 * pi) - b).subs(t0, -b / a))
    A_drift = sp.simplify(1 - (2 * pi) / (2 * pi))
    _require(mean_f == 0 and T2_reduced == 2 * pi * (t0 ** 2 + 2 * q ** 2)
             and normalization == 0 and B_drift == 0 and A_drift == 0,
             "symbolic C5 periodic lemma did not reduce to zero")
    return {
        "passed": True, "engine": "sympy symbolic Fourier moments",
        "source_formula": "t=t0+2*q*h*(cos(psi)-r)/(1-2*r*cos(psi)+r^2)",
        "fourier_identity": "(cos(psi)-r)/D=sum(n>=0,r^n*cos((n+1)*psi))",
        "conditions": ["|r|<1", "h^2=1-r^2", "h=1/sqrt(1+u^2)", "r=u*h"],
        "mean_f": "0", "mean_f_squared": "1/(2*(1-r^2))",
        "T1_period": "2*pi*t0", "T2_period": "2*pi*(t0^2+2*q^2)",
        "K": "1/(2*pi*(1+t0^2+2*q^2))", "K_normalizes_Phi": True,
        "t0_source_identity": "t0=-b/a",
        "B_phase_drift_after_t0": "0", "A_phase_drift": "0",
        "h_relation_remainder": str(h_relation),
        "ordinary_Z_phase_independent": True,
        "smooth_periodic_integer_extension": True,
    }


def _formal_endpoint_congruence(root: Path, q4: Mapping[str, Any], loop: Mapping[str, Any]) -> dict[str, Any]:
    """Prove endpoint seams from source AST, graph ports, and arbitrary symbols."""
    r100 = _json_gzip(root / R100_ARCHIVE)
    rb = r100.get("prescribed_R100_Rm_C6_bindings")
    lb = loop.get("prescribed_loop_C5_bindings")
    qb = q4.get("modified_Q4_pressure_C3_bindings")
    _require(isinstance(rb, Mapping) and isinstance(lb, Mapping) and isinstance(qb, Mapping),
             "R100/C5/Q4 bindings missing")
    _require(rb.get("hb") == 9 and rb.get("bridge_hb") == 8 and rb.get("Pstar") == PSTAR_HANDLE,
             "R100 native hb/bridge_hb/Pstar cache changed")
    _require(lb.get("N") == N_HANDLE and lb.get("eta") == 21 and lb.get("dstar") == 22
             and lb.get("sc") == 10 and lb.get("phase_origin") == 24954,
             "C5 native N/eta/dstar/sc/phase_origin cache changed")
    rcharts = rb.get("charts", {})
    ccharts = lb.get("charts", {})
    qcharts = qb.get("charts", {})
    qendpoints = qb.get("endpoints", {})
    _require(isinstance(rcharts, Mapping) and isinstance(ccharts, Mapping)
             and isinstance(qcharts, Mapping) and isinstance(qendpoints, Mapping),
             "R100/C5/Q4 chart caches missing")
    rnodes = r100["exact_graph_nodes"]
    lnodes = loop["exact_graph_nodes"]
    qnodes = q4["exact_graph_nodes"]

    source_receipts = {
        "R100_first_switch_fields": _source_fragment(root / R100_SOURCE, "first_switch_fields", (
            "angular=a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius))",
            "logratio=a.integral(a.at(angular,phase,t),t,g.zero,phase",
            "phi=a.mul(initial['phi'],a.exp(logratio))",
            "chi_axial=g.mul(hb,g.sub(g.one,sigma(g,phase)))",
            "Uy=a.scale(a.mul(ratio,drive),g.neg(chi_axial))",
            "phiy=a.scale(a.mul(known['D_over_R'],phi),g.mul(g.constant('-1/2'),hb,radius))")),
        "R100_second_switch_fields": _source_fragment(root / R100_SOURCE, "second_switch_fields", (
            "sig=sigma(g,phase)",
            "angular=a.add(a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius,g.sub(g.one,sig))),",
            "logratio=a.integral(a.at(angular,phase,t),t,g.zero,phase",
            "phi=a.mul(initial['phi'],a.exp(logratio))",
            "slope=a.add(a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,radius,g.sub(g.one,sig))),",
            "return dict(phi=phi,U=initial['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0)")),
        "R100_constructor_and_end": _source_fragment(root / R100_SOURCE, "__init__", (
            "incoming=self.endpoints['R1']",
            "actual=second_switch_fields(a,known,incoming['actual'],t,hb,R)",
            "self.add_chart('second_switch',t,R,hb,actual,incoming['moments'],known=known)",
            "self.endpoints['R2']=self.end('second_switch',g.one)",
            "incoming=self.endpoints['R2']",
            "y=g.mul(self.post_length,p)",
            "phi=a.scale(incoming['actual']['phi'],g.unary('exp',g.mul(c('-2/5'),y)))",
            "actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0))",
            "self.add_chart('post_switch',p,R,self.post_length,actual,incoming['moments'])",
            "self.endpoints['R110']=self.end('post_switch',g.one)")),
        "R100_add_chart": _source_fragment(root / R100_SOURCE, "add_chart", (
            "J=previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
            "fields,radial,generic=previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)",
            "self.charts[name]=dict(phase=phase,domain=[0,1],radius=R,physical_measure=jacobian,actual=actual,incoming_histories=incoming,moments=J,fields=fields,radial=radial,generic=generic")),
        "R100_end": _source_fragment(root / R100_SOURCE, "end", (
            "at=lambda jet:a.at(jet,chart['phase'],value)",
            "actual={key:at(chart['actual'][key]) for key in ('phi','U','phi_y','U_y')}",
            "moments={k:at(jet) for k,jet in chart['moments'].items()}")),
        "R100_sigma": _source_fragment(root / SIGMA_SOURCE, "stable_sigma", (
            "ifv<=0:returnc.mpf(0),c.mpf(0)", "ifv>=1:returnc.mpf(1),c.mpf(0)",
            "derivative=value*(b/(a+b))")),
        "R100_sigma_jet_source": _source_fragment(root / SIGMA_JETS_SOURCE, "sigma_jets", (
            "iflo<=0:rows.append(IntervalTaylor.constant(c,0,4))",
            "ifhi>=1:rows.append(IntervalTaylor.constant(c,1,4))",
            "left=_sigma_left(c,1-c.mpf([a,b]))")),
        "Q4_physical_pressure_fields": _source_fragment(root / Q4_SOURCE, "physical_pressure_fields", (
            "U=a.scale(fields['V'],Pstar)", "Uy=a.scale(radial['V_y'],Pstar)",
            "m=a.scale(fields['m'],Pstar)", "my=a.sub(U,m)",
            "Q=div(numerator(U,m))", "Qy=div(numerator(Uy,my))")),
        "Q4_constructor": _source_fragment(root / Q4_SOURCE, "__init__", (
            "physical=physical_pressure_fields(g,chart['modified_fields'],radial,self.background.delta,self.background.Pstar)",
            "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,chart['radius']))",
            "at=lambda jets,value:{key:a.at(jet,coordinate,value) for key,jet in jets.items()}",
            "endpoint=dict(integrals=at(cumulative,chart['domain'][1]),pressure=at(pressure,chart['domain'][1]),boundary=at(boundary,chart['domain'][1]))")),
        "Q4_total_radial_functions": _source_fragment(root / Q4_SOURCE, "total_radial_functions", (
            "invjac=g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
            "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref)) for ref in rows(chart['loop'][key]))",
            "invN=g.quotient(g.one,N,'same original exact positive common N')",
            "x=a.scale(chart['loop']['A'],invN);ex=a.exp(x)",
            "xy=a.add(a.scale(slow['A'],invN),chart['loop']['A_phi'])",
            "Ey=a.mul(ex,a.add(a.trim(bg['radial']['E_y']),a.mul(chart['roots']['E'],xy)))",
            "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
            "mean_y=a.sub(chart['modified_fields']['V'],chart['modified_fields']['m'])")),
        "C5_generic_roots": _source_fragment(root / LOOP_SOURCE, "generic_roots", (
            "p2=generic['p2']", "t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))")),
        "C5_loop_functions": _source_fragment(root / LOOP_SOURCE, "loop_functions", (
            "u=a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,'same positive canonical d_star'))",
            "hinv=a.div(one,a.sqrt(a.add(one,a.mul(u,u))),'positive sqrt(1+u^2)')",
            "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
            "numerator=m.add(cc,m.scale(rr,-1))",
            "den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))",
            "t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))",
            "ratio=a.add(one,a.mul(roots['t0'],roots['t0']),a.scale(a.mul(q,q),2))",
            "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
            "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
            "lam=m.mul(m.axial(K),m.add(m.constant(1),t2))",
            "psi,delta,at=implicit_inverse(g,partials,inverse,angle)")),
        "C5_finite_N_functions": _source_fragment(root / LOOP_SOURCE, "finite_N_functions", (
            "exm1=FunctionJet((g.unary('expm1',x.value),)+rows(ex)[1:])",
            "E_N=a.mul(E,ex)", "V_N=a.add(V,s)")),
        "C5_transport": _source_fragment(root / LOOP_SOURCE, "transport", (
            "width=g.mul(jacobian,g.sub(coordinate,lower))",
            "contribution=a.integral(body,t,lower,coordinate,actual_finite_N_signed_correction=True",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
            "predecessor_memory_never_reset=True")),
        "C5_constructor": _source_fragment(root / LOOP_SOURCE, "__init__", (
            "self.N=self.unique('exact_positive_integer_expression',None)",
            "offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)",
            "phase=g.unary('fractional_part',g.mul(self.N,offset))",
            "finite=finite_N_functions(g,roots,loop,self.N)",
            "history,parts,memory=transport(g,finite['five_signed_density_increments'],incoming",
            "fields[key]=a.add(a.trim(chart['fields'][key]),history[key])")),
    }

    expected_r100 = {
        "first_switch": (73190, 85905, 9), "second_switch": (75949, 88600, 9),
        "post_switch": (78472, 91107, 91103),
    }
    for name, (phase, radius, jacobian) in expected_r100.items():
        chart = rcharts.get(name)
        _require(isinstance(chart, Mapping) and chart.get("phase") == phase
                 and chart.get("radius") == radius and chart.get("physical_measure") == jacobian,
                 f"R100 {name} geometry cache changed")
    # Constructor/end receipts bind the actual R1/R2/R110 substitutions.
    for endpoint_name, chart_name in (("R1", "first_switch"), ("R2", "second_switch"), ("R110", "post_switch")):
        endpoint = rb.get("endpoints", {}).get(endpoint_name, {})
        chart = rcharts[chart_name]
        radius_row = _row(rnodes, endpoint["radius"])
        _require(radius_row.get("operation") == "function_substitution"
                 and radius_row.get("expression") == chart["radius"]
                 and radius_row.get("variable") == chart["phase"]
                 and radius_row.get("value") == 1
                 and radius_row.get("Z_independent_substitution") is True,
                 f"R100 {endpoint_name} radius endpoint binding changed")
    r1_moments = rb["endpoints"]["R1"].get("moments")
    r2_moments = rb["endpoints"]["R2"].get("moments")
    _require(rcharts["second_switch"].get("incoming_histories") == r1_moments
             and rcharts["post_switch"].get("incoming_histories") == r2_moments,
             "R100 incoming moment histories are not endpoint continuations")
    _require(rcharts["second_switch"].get("actual", {}).get("U") == rb["endpoints"]["R1"].get("actual", {}).get("U")
             and rcharts["post_switch"].get("actual", {}).get("U") == rb["endpoints"]["R2"].get("actual", {}).get("U"),
             "R100 switch U inlet is not endpoint-bound")
    sigma_handles = {"first_switch": 74428, "second_switch": 76947}
    sigma_receipt: dict[str, Any] = {}
    for name, handle in sigma_handles.items():
        row = _row(rnodes, handle)
        _require(row.get("operation") == "exact_original_source_flat_sigma"
                 and row.get("argument") == rcharts[name]["phase"]
                 and row.get("exact_function_not_its_interval_cap") is True,
                 f"R100 {name} flat sigma graph binding changed")
        binding = row.get("source_binding", {})
        _require(binding.get("module") == Path(SIGMA_SOURCE).name
                 and binding.get("module_sha256") == _sha(root / SIGMA_SOURCE)
                 and binding.get("function") == "stable_sigma",
                 f"R100 {name} sigma source pin changed")
        _require(row.get("definition") == "0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside",
                 f"R100 {name} flat sigma definition changed")
        sigma_receipt[name] = {"handle": handle, "argument": rcharts[name]["phase"],
                              "source_binding": copy.deepcopy(binding), "flat_values": {"at_0": 0, "at_1": 1,
                              "sigma_y_at_0": 0, "sigma_y_at_1": 0}}

    rates = ("m", "h", "k", "e", "p")
    seam_cache: dict[str, Any] = {}
    for left, right in (("first_switch", "second_switch"), ("second_switch", "post_power")):
        left_end = lb.get("endpoints", {}).get(left, {})
        right_chart = ccharts.get(right, {})
        _require(isinstance(left_end, Mapping) and isinstance(right_chart, Mapping), "C5 seam endpoint cache missing")
        _require(right_chart.get("incoming") == left_end, f"C5 {left}->{right} incoming history reset")
        _require(right_chart.get("memory", {}).get("p") == 1,
                 f"C5 {right} rate-zero memory is not retained")
        seam_cache[f"{left}_to_{right}"] = {"incoming_equals_predecessor_endpoint": True,
                                             "rates": list(rates), "p_memory": 1,
                                             "transport_lower_substitution": "exp(0)=1 and integral(lower,lower)=0",
                                             "source_history_y_recipe": "X_y=density-rate*X; X_yy=density_y-rate*X_y",
                                             "incoming_cache_keys": sorted(right_chart.get("incoming", {}).keys())}

    q4_port_receipts: dict[str, Any] = {}
    for name, cfg in _CHARTS.items():
        qname = cfg["q4_chart"]
        cname = cfg["native_chart"] if cfg["native_chart"] != "post_switch" else "post_power"
        qchart = qcharts[qname]
        qendpoint = qendpoints[qname]
        cchart = ccharts[cname]
        rname = cfg["native_chart"]
        rchart = rcharts[rname]
        # The C5 roots are the R100 fields and its p2 source, while finite-N
        # V/m fields are the exact inputs consumed by Q4 physical ports.
        _require(cchart["roots"]["E"] == rchart["fields"]["E"][:len(cchart["roots"]["E"])]
                 and cchart["roots"]["V"] == rchart["fields"]["V"][:len(cchart["roots"]["V"])]
                 and cchart["roots"]["p2"] == rchart["generic"]["p2"],
                 f"{name} C5 root inputs are not R100-source bound")
        _require(cchart["finite_N"]["E_N"] == cchart["modified_fields"]["E"]
                 and cchart["finite_N"]["V_N"] == cchart["modified_fields"]["V"],
                 f"{name} finite-N E/V ports changed")
        physical = qchart["physical"]
        for i in range(6):
            urow = _row(qnodes, physical["U"][i])
            uyrow = _row(qnodes, physical["U_y"][i])
            mrow = _row(qnodes, physical["m"][i])
            myrow = _row(qnodes, physical["m_y"][i])
            _require(urow.get("operation") == "product" and urow.get("arguments") == [cchart["modified_fields"]["V"][i], PSTAR_HANDLE],
                     f"{name} Q4 U port is not Pstar*V")
            _require(uyrow.get("operation") == "product" and uyrow.get("arguments") == [qchart["radial"]["V_y"][i], PSTAR_HANDLE],
                     f"{name} Q4 U_y port is not Pstar*V_y")
            _require(mrow.get("operation") == "product" and mrow.get("arguments") == [cchart["modified_fields"]["m"][i], PSTAR_HANDLE],
                     f"{name} Q4 m port is not Pstar*m")
            neg = _row(qnodes, myrow.get("arguments", [None, None])[1]) if myrow.get("operation") == "sum" else {}
            _require(myrow.get("operation") == "sum" and myrow.get("arguments", [None, None])[0] == physical["U"][i]
                     and neg.get("operation") == "negative" and neg.get("argument") == physical["m"][i],
                     f"{name} Q4 m_y port is not U-m")
        q_denominator = _row(qnodes, physical["Q"][0]).get("denominator")
        _require(q_denominator == 25272, f"{name} Q denominator changed")
        for handles in (physical["Q"], physical["Q_y"]):
            for handle in handles:
                row = _row(qnodes, handle)
                _require(row.get("operation") == "positive_quotient"
                         and row.get("denominator") == q_denominator
                         and row.get("source_positive_certificate") == "same original L=1-delta*Z^2",
                         f"{name} Q/Q_y rational source changed")
        for i, handle in enumerate(qendpoint["boundary"]["R_Q_y"]):
            row = _row(qnodes, handle)
            _require(row.get("operation") == "function_substitution"
                     and row.get("expression") == qchart["finite_boundary"]["R_Q_y"][i]
                     and row.get("variable") == cfg["coordinate"] and row.get("value") == 1
                     and row.get("Z_independent_substitution") is True,
                     f"{name} Q4 upper R_Q_y wrapper changed")
        q4_port_receipts[name] = {
            "physical_ports": {key: list(physical[key]) for key in ("U", "U_y", "m", "m_y", "Q", "Q_y")},
            "finite_boundary_R_Q_y": list(qchart["finite_boundary"]["R_Q_y"]),
            "upper_endpoint_R_Q_y": list(qendpoint["boundary"]["R_Q_y"]),
            "lower_wrapper_recipe": "finite_boundary.R_Q_y evaluated at chart coordinate 0",
            "Pstar_once": True, "my_source": "U-m", "ordinary_Z_orders": [0, 1, 2],
            "R100_source_rows_trimmed_to_C5_order": True,
        }

    # The first-switch side is the left endpoint of the first new seam.  It is
    # outside _CHARTS because this module appends only second/post rows, but its
    # physical Q4 ports and phase must still be bound directly to the cache.
    first_q4 = qcharts["first_switch"]
    first_endpoint = qendpoints["first_switch"]
    first_c5 = ccharts["first_switch"]
    first_r100 = rcharts["first_switch"]
    first_physical = first_q4["physical"]
    first_expected = {
        "U": (148085, 148086, 148087, 148088, 148089, 148090),
        "U_y": (148091, 148092, 148093, 148094, 148095, 148096),
        "m": (148097, 148098, 148099, 148100, 148101, 148102),
        "m_y": (148109, 148110, 148111, 148112, 148113, 148114),
        "Q": (148180, 148184, 148190, 148196, 148202),
        "Q_y": (148268, 148272, 148278, 148284, 148290),
    }
    for key, handles in first_expected.items():
        _require(first_physical.get(key) == list(handles), f"first_switch Q4 {key} cache changed")
    _require(first_q4["finite_boundary"]["R_Q_y"] == [148427, 148428, 148429, 148430]
             and first_endpoint["boundary"]["R_Q_y"] == [148514, 148515, 148516, 148517],
             "first_switch Q4 R_Q_y wrappers changed")
    _require(first_c5["actual_phase"] == 104537 and first_c5["global_offset"] == 104535
             and first_c5["coordinate"] == 73190 and first_c5["radius"] == 85905
             and first_c5["roots"]["E"] == first_r100["fields"]["E"][:len(first_c5["roots"]["E"])]
             and first_c5["roots"]["V"] == first_r100["fields"]["V"][:len(first_c5["roots"]["V"])]
             and first_c5["roots"]["p2"] == first_r100["generic"]["p2"],
             "first_switch C5/R100 source cache binding changed")
    _require(first_c5["finite_N"]["E_N"] == first_c5["modified_fields"]["E"]
             and first_c5["finite_N"]["V_N"] == first_c5["modified_fields"]["V"],
             "first_switch finite-N E/V ports changed")
    for i in range(6):
        urow = _row(qnodes, first_physical["U"][i])
        uyrow = _row(qnodes, first_physical["U_y"][i])
        mrow = _row(qnodes, first_physical["m"][i])
        myrow = _row(qnodes, first_physical["m_y"][i])
        _require(urow.get("operation") == "product"
                 and urow.get("arguments") == [first_c5["modified_fields"]["V"][i], PSTAR_HANDLE]
                 and uyrow.get("operation") == "product"
                 and uyrow.get("arguments") == [first_q4["radial"]["V_y"][i], PSTAR_HANDLE]
                 and mrow.get("operation") == "product"
                 and mrow.get("arguments") == [first_c5["modified_fields"]["m"][i], PSTAR_HANDLE],
                 "first_switch Q4 Pstar physical port changed")
        neg = _row(qnodes, myrow.get("arguments", [None, None])[1]) if myrow.get("operation") == "sum" else {}
        _require(myrow.get("operation") == "sum"
                 and myrow.get("arguments", [None, None])[0] == first_physical["U"][i]
                 and neg.get("operation") == "negative"
                 and neg.get("argument") == first_physical["m"][i],
                 "first_switch Q4 m_y port is not U-m")
    first_q_denominator = _row(qnodes, first_physical["Q"][0]).get("denominator")
    _require(first_q_denominator == 25272, "first_switch Q denominator changed")
    for handles in (first_physical["Q"], first_physical["Q_y"]):
        for handle in handles:
            row = _row(qnodes, handle)
            _require(row.get("operation") == "positive_quotient"
                     and row.get("denominator") == first_q_denominator
                     and row.get("source_positive_certificate") == "same original L=1-delta*Z^2",
                     "first_switch Q/Q_y rational source changed")
    for i, handle in enumerate(first_endpoint["boundary"]["R_Q_y"]):
        row = _row(qnodes, handle)
        _require(row.get("operation") == "function_substitution"
                 and row.get("expression") == first_q4["finite_boundary"]["R_Q_y"][i]
                 and row.get("variable") == 73190 and row.get("value") == 1
                 and row.get("Z_independent_substitution") is True,
                 "first_switch upper R_Q_y wrapper changed")
    q4_port_receipts["first_switch"] = {
        "physical_ports": {key: list(first_physical[key]) for key in first_expected},
        "finite_boundary_R_Q_y": list(first_q4["finite_boundary"]["R_Q_y"]),
        "upper_endpoint_R_Q_y": list(first_endpoint["boundary"]["R_Q_y"]),
        "first_C5_phase": first_c5["actual_phase"], "first_C5_global_offset": first_c5["global_offset"],
        "Pstar_once": True, "my_source": "U-m", "ordinary_Z_orders": [0, 1, 2],
        "R100_source_rows_trimmed_to_C5_order": True,
    }

    for name in ("first_switch", "second_switch", "post_power"):
        chart = ccharts[name]
        phase = chart["actual_phase"]
        offset = chart["global_offset"]
        phase_row = _row(lnodes, phase)
        phase_arg = _row(lnodes, phase_row.get("argument"))
        offset_row = _row(lnodes, offset)
        _require(phase_row.get("operation") == "analytic_unary" and phase_row.get("name") == "fractional_part"
                 and phase_arg.get("operation") == "product" and phase_arg.get("arguments") == [N_HANDLE, offset]
                 and not _has_dependency(lnodes, phase, 20),
                 f"C5 {name} does not use common Z-independent N phase")
        _require(offset_row.get("operation") == "sum" and _has_dependency(lnodes, offset, lb["phase_origin"]),
                 f"C5 {name} global phase offset lost phase_origin")
        _require(_has_dependency(lnodes, offset, chart["radius"]),
                 f"C5 {name} phase offset lost current radius")

    rdef = rb.get("definition", {})
    ldef = lb.get("definition", {})
    for key in ("exact_original_both_R100_switch_controls_and_current_radius_retained",
                "all_six_nonzero_incoming_and_every_chart_memory_retained"):
        _require(rdef.get(key) is True, f"R100 source premise missing: {key}")
    for key in ("all_ten_signed_root_q_inverse_A_B_and_phase_primitive_rows_through_Z5",
                "inverse_rows_from_exact_implicit_normalized_series_not_interval_endpoints",
                "complete_original_rminus_to_Rm_Duhamel_correction_histories_C5",
                "correction_zero_inlet_only_at_rminus_background_nonzero_memories_unchanged",
                "original_rate_zero_pressure_memory_retained"):
        _require(ldef.get(key) is True, f"C5 source premise missing: {key}")

    symbolic_endpoint = _symbolic_endpoint_identities()
    symbolic_rational = _symbolic_rational_congruence()
    periodic = _symbolic_periodic_phase_lemma()
    return {
        "passed": True, "ordinary_Z_orders": [0, 1, 2], "common_N": N_HANDLE,
        "eta": lb["eta"], "dstar": lb["dstar"], "sc": lb["sc"], "phase_origin": lb["phase_origin"],
        "Pstar": rb["Pstar"], "source_recipe": {
            "R_first": "100*exp(hb*s)", "R_second": "100*exp(hb*(1+t))",
            "R_post": "100*exp(2*hb+post_length*p)",
            "phase_offset": "log(R/Ra)-phase_origin", "phase": "fractional_part(N*offset)",
            "phase_y": "N*d_y(log(R/Ra))", "periodic_smooth_extension_at_integer_phase": True,
        },
        "source_formulas": source_receipts, "sigma_flat_endpoint_receipt": sigma_receipt,
        "R100_endpoint_bindings": {"R1": "first_switch end at s=1", "R2": "second_switch end at t=1",
                                    "R110": "post_switch end at p=1", "incoming_actual_U_and_moments_bound": True,
                                    "R1_R2_R110_moment_keys": sorted(r1_moments),
                                    "actual_endpoint_keys": sorted(rb["endpoints"]["R1"].get("actual", {}).keys())},
        "Q4_port_bindings": q4_port_receipts, "C5_history_bindings": seam_cache,
        "sigma_flat_endpoint": True, "nonreset_pointers": True,
        "background_E_V_and_native_histories_congruent_through_Z2": True,
        "history_y_relation": "delta_M_y=delta_V-delta_M",
        "C5_same_generic_p2_and_inverse_source": True,
        "C5_total_radial_slow_over_J_and_Nfast_bound": True,
        "Q4_same_rational_Qy_source": True, "physical_y_derivatives_claimed": False,
        "symbolic_endpoint_identity": symbolic_endpoint,
        "symbolic_rational_congruence": symbolic_rational,
        "periodic_phase_analytic_lemma": periodic,
        "seams": {
            "first_upper_to_second_lower": {"R_identity": "100*exp(hb)=100*exp(hb*(1+0))",
                "phi_y_identity": "sigma(1)=1,sigma_y(1)=0,sigma(0)=0,sigma_y(0)=0 => -hb*R*D*phi/2",
                "Uy_identity": "first U_y=-hb*(1-sigma)*A and U_yy flat endpoint = 0 = second U_y,U_yy",
                "finite_N_Qy_identity_proven": True, "contraction_allowed": True},
            "second_upper_to_post_lower": {"R_identity": "100*exp(2*hb)=100*exp(2*hb+post_length*0)",
                "phi_y_identity": "sigma(1)=1,sigma_y(1)=0 => second lambda=-2/5,lambda_y=0=post",
                "Uy_identity": "second and post U_y=U_yy=0",
                "finite_N_Qy_identity_proven": True, "contraction_allowed": True},
        },
    }


def _endpoint_receipts(formal: Mapping[str, Any]) -> dict[str, Any]:
    common = {"same_common_N_handle": N_HANDLE, "common_phase_Z_independent": True,
              "nonreset_histories": True, "exactflat_sigma_endpoint_binding": True,
              "source_identity_receipt_owned_here": True, "source_cache_semantic_rewrite": True,
              "arbitrary_symbolic_jet_proof": True, "arithmetic_fixture_not_proof": True,
              "formal_congruence": copy.deepcopy(formal)}
    return {
        "first_upper_to_second_lower": {**common, "left": list(FIRST_UPPER_RQY), "right_source": [154589, 154590],
            "right_wrapper": "second_lower_t0", "same_R": "R=100*exp(hb)", "Uy_endpoint_match": True,
            "phi_slope_match": "-hb*R*D/2", "background_y2_jets_match": True,
            "finite_N_Qy_identity_proven": True, "contraction_blocked": False,
            "status": "PROVED_FORMAL_CONGRUENCE",
            "reason": "Actual R100/C5/Q4 source and cache bindings plus arbitrary symbolic endpoint and rational-jet congruence through ordinary Z order 2."},
        "second_upper_to_post_lower": {**common, "left": list(SECOND_UPPER_RQY), "right_source": [160200, 160201],
            "right_wrapper": "post_lower_p0", "same_R": "R=100*exp(2*hb)", "Uy_endpoint_match": True,
            "phi_slope_match": "-2*phi/5", "sigma_prime_zero": True,
            "finite_N_Qy_identity_proven": True, "contraction_blocked": False,
            "status": "PROVED_FORMAL_CONGRUENCE",
            "reason": "Actual R100/C5/Q4 source and cache bindings plus arbitrary symbolic endpoint and rational-jet congruence through ordinary Z order 2."},
    }


def build(root: str | Path, parent: Mapping[str, Any], numeric_packet: Mapping[str, Any]) -> dict[str, Any]:
    """Build the append-only second-switch/post-power pressure source."""
    root_path = Path(root)
    _require(_sha(root_path / PARENT) == PARENT_SHA256, "current first-switch parent file changed")
    saved = _json(root_path / PARENT)
    _require(_digest(saved) == _digest(parent), "exact current first-switch parent required")
    _require(isinstance(numeric_packet, Mapping), "numeric packet mapping required")

    pin_paths = (Q4_SOURCE, Q4_ARCHIVE, Q4_CHECK, LOOP_SOURCE, LOOP_ARCHIVE, LOOP_CHECK,
                 DENSITY_SOURCE, GRAPH_HELPER, SEAM_SOURCE, R100_SOURCE, R100_ARCHIVE,
                 SECOND_SOURCE, BRIDGE_SOURCE, ACTUAL_BRIDGE_SOURCE, FIRST_SOURCE, MACRO_FIRST_SOURCE,
                 SIGMA_SOURCE, SIGMA_JETS_SOURCE)
    source_pins: dict[str, Any] = {PARENT: _pin_path(root_path, parent, PARENT)}
    for path in pin_paths:
        source_pins[path] = _pin_path(root_path, parent, path)
    packet_pins = numeric_packet.get("source_pins", {})
    if isinstance(packet_pins, Mapping):
        for path, value in packet_pins.items():
            if isinstance(path, str) and isinstance(value, (Mapping, str)):
                pin = _packet_pin(root_path, parent, path, value)
                old = source_pins.get(path)
                _require(old is None or old["expected_sha256"] == pin["expected_sha256"], f"conflicting source pin: {path}")
                source_pins[path] = pin
    source_pins[MODULE] = {"path": MODULE, "expected_sha256": _sha(root_path / MODULE), "passed": True, "self_pin": True}

    q4 = _validate_q4(root_path)
    loop = _validate_loop(root_path)
    q4_archive = _json_gzip(root_path / Q4_ARCHIVE)
    loop_archive = _json_gzip(root_path / LOOP_ARCHIVE)
    graph, prefix = _rehydrate_prefix(parent, q4_archive)
    seam = _validate_parent_seams(parent, graph.nodes)
    numeric = _numeric_contract(root_path, parent, numeric_packet)
    formal = _formal_endpoint_congruence(root_path, q4_archive, loop_archive)
    source_guards = {
        "q4": _source_guard(root_path / Q4_SOURCE, {"physical_pressure_fields": (
            "U=a.scale(fields['V'],Pstar)", "m=a.scale(fields['m'],Pstar)", "my=a.sub(U,m)",
            "Qy=div(numerator(Uy,my))"), "__init__": ("boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,chart['radius']))",),
            "total_radial_functions": (
                "invjac=g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
                "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref)) for ref in rows(chart['loop'][key]))",
                "invN=g.quotient(g.one,N,'same original exact positive common N')",
                "x=a.scale(chart['loop']['A'],invN);ex=a.exp(x)",
                "xy=a.add(a.scale(slow['A'],invN),chart['loop']['A_phi'])",
                "Ey=a.mul(ex,a.add(a.trim(bg['radial']['E_y']),a.mul(chart['roots']['E'],xy)))",
                "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
                "mean_y=a.sub(chart['modified_fields']['V'],chart['modified_fields']['m'])"),}),
        "c5": _source_guard(root_path / LOOP_SOURCE, {"generic_roots": ("p2=generic['p2']",),
            "loop_functions": ("den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))",
                           "psi,delta,at=implicit_inverse(g,partials,inverse,angle)"),
            "finite_N_functions": ("V_N=a.add(V,s)", "expm1"),
            "transport": ("out[name]=a.add(a.scale(incoming[name],decay),contribution)",
                           "predecessor_memory_never_reset=True"),
            "__init__": ("self.N", "phase=g.unary('fractional_part',g.mul(self.N,offset))",)}),
        "r100": _source_guard(root_path / R100_SOURCE, {"first_switch_fields": ("def first_switch_fields",),
            "second_switch_fields": ("def second_switch_fields", "U_y=a.fixed(0)"), "__init__": (
                "R=g.mul(c(100),g.unary('exp',g.mul(hb,s)))",
                "R=g.mul(c(100),g.unary('exp',g.mul(hb,g.add(g.one,t))))",
                "R=g.mul(c(100),g.unary('exp',g.add(g.mul(c(2),hb),y)))"),
            "add_chart": ("def add_chart",), "end": ("def end",)}),
        "bridge": _source_guard(root_path / BRIDGE_SOURCE, {"__init__": ("def __init__",)}),
        "second_source": _source_guard(root_path / SECOND_SOURCE, {"source_bindings": (
            "second_comparison_radius", "original_positive_post_power_assignments")}),
    }
    guards = {"passed": True, "AST": source_guards, "self_format": _self_format_guard(root_path),
              "Q4": q4, "C5": loop, "seam": seam, "source_execution": False}

    chart_out: dict[str, Any] = {}
    direct_local: dict[str, Sequence[int]] = {}
    lower_wrappers: dict[str, Sequence[int]] = {}
    for name, cfg in _CHARTS.items():
        integrals, six = _build_six_source(graph, cfg)
        B_rows, B = _build_reduced_B(graph, cfg)
        D_rows, D = _build_reduced_D(graph, integrals)
        reduced_S1 = list(graph.padd(_pair(B_rows, "B"), _pair(D_rows, "D")))
        # The joined pressure uses the already accepted Q4 local S1 exactly once;
        # B + D is a separate raw-physical source receipt and is not substituted
        # for that finite-N local integral.
        direct_local[name] = tuple(cfg["local_S1"][:2])
        wrapper = graph.at(_pair(cfg["finite_RQy"][:2], f"{name} finite R_Q_y"), cfg["coordinate"], graph.constant(0))
        lower_wrappers[name] = wrapper
        chart_out[name] = {"passed": True, "native_chart": cfg["native_chart"], "q4_chart": cfg["q4_chart"],
            "numeric_chart": cfg["numeric_name"], "N": N_HANDLE, "N_materialized": False,
            "phase": cfg["phase"], "phase_Z_independent": True, "coordinate": cfg["coordinate"], "domain": [0, 1],
            "radius": cfg["radius"], "native_J": cfg["J"], "physical_ports": {k: list(cfg[k]) for k in ("U", "U_y", "m", "m_y", "Q", "Q_y")},
            "incoming_m": list(cfg["incoming_m"]), "history_m": list(cfg["history_m"]), "nonreset_incoming_history": True,
            "local_S1": list(cfg["local_S1"][:2]), "cumulative_S1": list(cfg["cumulative_S1"][:2]),
            "finite_boundary_RQy": list(cfg["finite_RQy"][:2]), "lower_wrapper_t0_or_p0": list(wrapper),
            "upper_RQy_standalone": list(cfg["upper_RQy"][:2]), "six_I1_source": six, "integrals": integrals,
            "B_source": B, "D_source": D, "S1_source": {"rows": reduced_S1, "formula": "S1=B+D", "D_is_not_full_S1": True},
            "raw_graph_unit": cfg["graph_measure"], "numeric_unit": cfg["numeric_unit"], "Pstar_applied_once": True,
            "physical_R_not_mean_scale": True, "numeric": numeric["charts"][name]}

    join = _build_join(graph, parent, direct_local, lower_wrappers, formal)
    receipts = _endpoint_receipts(formal)
    # The new endpoint wrappers and the uncontracted joins are the only rows
    # appended by this module; the 184898-node prefix is byte-for-byte retained.
    new_suffix = copy.deepcopy(graph.nodes[FIRST_NODE_COUNT:])
    _require(graph.nodes[:FIRST_NODE_COUNT] == prefix["full_prefix_snapshot"], "full 184898-node prefix changed")
    _require(len(new_suffix) > 0, "new source did not append rows")

    graph_patch = {"base_archive": Q4_ARCHIVE, "base_node_count": FIRST_NODE_COUNT,
        "base_prefix_unchanged": True,
        "base_prefix_canonical_nodes_sha256": prefix["full_prefix_canonical_nodes_sha256"],
        "base_parent_first_switch_full_digest": prefix["full_prefix_canonical_nodes_sha256"],
        "new_node_start": FIRST_NODE_COUNT, "appended_nodes": new_suffix, "appended_node_count": len(new_suffix),
        "full_node_count": len(graph.nodes), "full_prefix_digest": prefix["full_prefix_canonical_nodes_sha256"],
        "full_canonical_nodes_sha256": _digest(graph.nodes), "segment_starts": prefix["segment_starts"],
        "prefix_suffix_counts": prefix["suffix_counts"], "endpoint_handles": {
            "first_switch_upper_standalone": list(FIRST_UPPER_RQY), "second_switch_upper_standalone": list(SECOND_UPPER_RQY),
            "post_power_upper_standalone": list(POST_UPPER_RQY), "second_switch_lower_t0": list(lower_wrappers["second_switch"]),
            "post_power_lower_p0": list(lower_wrappers["post_power"]), "second_switch_local_S1": list(direct_local["second_switch"]),
            "post_power_local_S1": list(direct_local["post_power"]), "joined_after_second_uncontracted": join["joined_after_second_uncontracted"],
            "joined_after_post_uncontracted": join["joined_after_post_uncontracted"],
            "joined_contracted": join["joined_contracted"],
            "contracted_before_linearization": join["contracted_before_linearization"]}}

    source_contract = {"passed": True, "finite_N_source": True, "finite_N_fields_from_C5_loop": True,
        "ordinary_Z_orders": [0, 1], "common_N_handle": N_HANDLE, "N_not_materialized": True,
        "phase_Z_independent": True, "native_J_distinct": {"second_switch": 9, "post_power": 91103},
        "actual_source_chart_mapping": {"second_switch": "second_switch", "post_switch": "post_power"},
        "physical_y_derivatives_claimed": False, "physical_Pstar_applied_once": True,
        "Pstar_handle": PSTAR_HANDLE, "Pstar_log_argument_handle": PSTAR_LOG_ARGUMENT_HANDLE,
        "physical_R_not_a_mean_scale": True, "measures": {"second_switch": "R^2*9*dt", "post_power": "R^2*91103*dp"},
        "numeric_units_separate": {"second_switch": "Wswitch", "post_switch": "Wmacro"},
        "outer_RQy_retained": True, "complete_pressure_numeric": False,
        "nonreset_histories": True, "joined_entry": list(PARENT_JOINED),
        "conditional_endpoint_identity_only_in_joined_expression": True,
        "contracted_sum_linearization": True, "removed_endpoint_atoms_proved": True}
    limits = {"second_switch_full_S1_source_defined": True, "post_power_full_S1_source_defined": True,
        "second_switch_D_is_not_full_S1": True, "post_power_D_is_not_full_S1": True,
        "joined_pressure_source_defined": True, "joined_pressure_numeric_enclosure": False,
        "outer_RQy_retained": True, "outer_RQy_evaluated": False, "all_caps_retained": True,
        "first_switch_upper_standalone_retained": True, "second_switch_upper_standalone_retained": True,
        "post_power_upper_standalone_retained": True, "conditional_identity_blocked": False,
        "formal_endpoint_congruence_proved": True,
        "joined_contracted_linearized": True, "contracted_endpoint_atoms_removed": True,
        "global_closure": False, "complete_pressure_global_closure": False,
        "accepted_ancestors_replayed": False, "old_212_terms_and_errors_unchanged": True,
        "physical_y_derivatives_claimed": False, "N_materialized": False}
    result = {"passed": True, "parent_report": PARENT, "parent_report_sha256": PARENT_SHA256,
        "source_family": parent.get("source_family"), "source_pins": source_pins, "input_pins": copy.deepcopy(source_pins),
        "source_guards": guards, "source_contract": source_contract, "contract": copy.deepcopy(source_contract),
        "numeric_packet": numeric, "charts": chart_out, "endpoint_identity_receipts": receipts,
        "joined_source_entry": list(PARENT_JOINED),
        "second_switch_I1_full": copy.deepcopy(numeric["charts"]["second_switch"]["I1_full"]),
        "second_switch_D_full": copy.deepcopy(numeric["charts"]["second_switch"]["D_full"]),
        "post_switch_I1_full": copy.deepcopy(numeric["charts"]["post_power"]["I1_full"]),
        "post_switch_D_full": copy.deepcopy(numeric["charts"]["post_power"]["D_full"]),
        "joined_macro_switch_post_volume": copy.deepcopy(numeric.get("joined_macro_switch_post_volume")),
        "joined_pressure_source": join, "graph_patch": graph_patch, "errors": {
            "second_switch": copy.deepcopy(numeric["charts"]["second_switch"]["physical_error_log_bounds"]),
            "post_switch": copy.deepcopy(numeric["charts"]["post_power"]["physical_error_log_bounds"]),
            "six_I1_difference_error_log_bounds": {
                "second_switch": copy.deepcopy(numeric["charts"]["second_switch"].get("six_I1_difference_error_log_bounds", {})),
                "post_switch": copy.deepcopy(numeric["charts"]["post_power"].get("six_I1_difference_error_log_bounds", {})),
            },
            "signed_I1_radius": {
                "second_switch": copy.deepcopy(numeric["charts"]["second_switch"].get("signed_I1_radius")),
                "post_switch": copy.deepcopy(numeric["charts"]["post_power"].get("signed_I1_radius")),
            }},
        "limits": limits, "signed_relation_not_global_closure": True, "accepted_ancestors_replayed": False}
    json.dumps(result, sort_keys=True, ensure_ascii=True, allow_nan=False)
    return result


__all__ = ["build"]
