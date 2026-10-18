"""Source-only current-N long-reshape/reference pressure extension.

The accepted parent source is rehydrated from its serialized Q4 graph and
saved suffixes.  This module adds the two native charts and their two local
RQy seams.  It binds source formulas and graph ports without importing an
accepted owner, builder, checker, producer, norm, recurrence, quadrature, or
old exact-algebra implementation, and it never materializes a global scale.
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
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_switch_post_finite_N_I1.json"
PARENT_SHA256 = "9d0d957bae129853030bad4cd5c484eb27e8b13b0623601e458fc4e95e93bc2b"
MODULE = "src/openai_ns_reconstruction/original_long_reference_pressure_source.py"

Q4_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
Q4_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
Q4_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions_check.json"
LOOP_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
LOOP_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
LOOP_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check.json"
R100_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
R100_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.json.gz"
R100_CHECK = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions_check.json"
BRIDGE_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
DENSITY_SOURCE = "src/openai_ns_reconstruction/original_finite_N_density_kernel.py"
GRAPH_HELPER = "src/openai_ns_reconstruction/original_signed_graph_pairs.py"
SIGMA_SOURCE = BASE + "lei_ren_part1_paper_compliant_outer_initial.py"
SIGMA_JETS_SOURCE = BASE + "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
LONG_ENDPOINT_SOURCE = BASE + "lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint.py"
LONG_ENDPOINT_ARCHIVE = BASE + "lei_ren_part1_paper_compliant_current_original_long_reshape_endpoint.json"

Q4_NODE_COUNT = 183700
PARENT_NODE_COUNT = 185308
N_HANDLE = 18
PSTAR_HANDLE = 33572
PSTAR_LOG_ARGUMENT_HANDLE = 5
PARENT_JOINED = (185300, 185303)
POST_UPPER_RQY = (160287, 160288)


_CHARTS: dict[str, dict[str, Any]] = {
    "long_reshape": {
        "native_chart": "long_reshape",
        "q4_chart": "long_reshape",
        "phase": 111982,
        "global_offset": 111980,
        "coordinate": 79966,
        "radius": 79970,
        "J": 7,
        "U": (165720, 165721, 165722, 165723, 165724, 165725),
        "U_y": (165726, 165727, 165728, 165729, 165730, 165731),
        "m": (165732, 165733, 165734, 165735, 165736, 165737),
        "m_y": (165744, 165745, 165746, 165747, 165748, 165749),
        "Q": (165815, 165819, 165825, 165831, 165837),
        "Q_y": (165903, 165907, 165913, 165919, 165925),
        "finite_RQy": (166062, 166063, 166064, 166065),
        "local_S1": (166093, 166094, 166095, 166096),
        "cumulative_S1": (166097, 166098, 166099, 166100),
        "incoming_S1": (160271, 160272, 160273, 160274),
        "upper_RQy": (166149, 166150, 166151, 166152),
        "incoming_m": (111948, 111949, 111950, 111951, 111952, 111953),
        "history_m": (114195, 114196, 114197, 114198, 114199, 114200),
        "numeric_unit": "Wlong",
        "graph_measure": "R^2*T*dq",
    },
    "reference": {
        "native_chart": "reference",
        "q4_chart": "reference",
        "phase": 114390,
        "global_offset": 114388,
        "coordinate": 81520,
        "radius": 81524,
        "J": 81519,
        "U": (171440, 171441, 171442, 171443, 171444, 171445),
        "U_y": (171446, 171447, 171448, 171449, 171450, 171451),
        "m": (171452, 171453, 171454, 171455, 171456, 171457),
        "m_y": (171464, 171465, 171466, 171467, 171468, 171469),
        "Q": (171535, 171539, 171545, 171551, 171557),
        "Q_y": (171623, 171627, 171633, 171639, 171645),
        "finite_RQy": (171782, 171783, 171784, 171785),
        "local_S1": (171813, 171814, 171815, 171816),
        "cumulative_S1": (171817, 171818, 171819, 171820),
        "incoming_S1": (166133, 166134, 166135, 166136),
        "upper_RQy": (171869, 171870, 171871, 171872),
        "incoming_m": (114356, 114357, 114358, 114359, 114360, 114361),
        "history_m": (116603, 116604, 116605, 116606, 116607, 116608),
        "numeric_unit": "Wreference",
        "graph_measure": "R^2*reference_length*dq",
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
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("utf-8")
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
    _require(isinstance(row, Mapping) and isinstance(row.get("operation"), str),
             "serialized graph row malformed")
    return graph.node(
        row["operation"],
        **{key: value for key, value in row.items() if key != "operation"},
    )


def _has_dependency(
    nodes: Sequence[Mapping[str, Any]], root: int, target: int
) -> bool:
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
        for key in (
            "argument",
            "numerator",
            "denominator",
            "expression",
            "variable",
            "value",
            "integrand",
            "lower",
            "upper",
        ):
            child = item.get(key)
            if type(child) is int and child >= 0:
                pending.append(child)
        arguments = item.get("arguments")
        if isinstance(arguments, list):
            pending.extend(
                child for child in arguments if type(child) is int and child >= 0
            )
    return False


def _source_guard(
    path: Path, functions: Mapping[str, Sequence[str]]
) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        raw = stream.read()
    tree = ast.parse(raw, filename=str(path))
    found: dict[str, str] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            found[node.name] = ast.get_source_segment(raw, node) or ""
    result: dict[str, Any] = {"path": path.as_posix(), "functions": {}}
    for name, snippets in functions.items():
        body = found.get(name, "")
        compact = "".join(body.split())
        normalized = ["".join(token.split()) for token in snippets]
        matched = bool(body) and all(token in compact for token in normalized)
        result["functions"][name] = {
            "present": bool(body),
            "tokens": list(snippets),
            "matched": matched,
        }
        _require(matched, f"source guard failed: {path.name}:{name}")
    result["passed"] = True
    result["source_executed"] = False
    return result


def _source_fragment(
    path: Path, function: str, snippets: Sequence[str]
) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        raw = stream.read()
    tree = ast.parse(raw, filename=str(path))
    candidates = [
        ast.get_source_segment(raw, node) or ""
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function
    ]
    _require(candidates, f"source function missing: {path.name}:{function}")
    source = max(candidates, key=len)
    compact = "".join(source.split())
    normalized = ["".join(token.split()) for token in snippets]
    _require(
        all(token in compact for token in normalized),
        f"source formula binding failed: {path.name}:{function}",
    )
    return {
        "path": path.as_posix(),
        "function": function,
        "tokens": list(snippets),
        "source_sha256": hashlib.sha256(source.encode("utf-8")).hexdigest(),
        "passed": True,
        "source_executed": False,
    }


def _qualified_fragment(
    path: Path, class_name: str, method_name: str, snippets: Sequence[str]
) -> dict[str, Any]:
    """Bind one method without importing the source owner."""

    with path.open("r", encoding="utf-8", newline="") as stream:
        raw = stream.read()
    tree = ast.parse(raw, filename=str(path))
    candidates: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.ClassDef) or node.name != class_name:
            continue
        for child in node.body:
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)) and child.name == method_name:
                candidates.append(ast.get_source_segment(raw, child) or "")
    _require(candidates, f"source method missing: {path.name}:{class_name}.{method_name}")
    source = max(candidates, key=len)
    compact = "".join(source.split())
    normalized = ["".join(token.split()) for token in snippets]
    _require(
        all(token in compact for token in normalized),
        f"source formula binding failed: {path.name}:{class_name}.{method_name}",
    )
    return {
        "path": path.as_posix(),
        "class": class_name,
        "function": method_name,
        "tokens": list(snippets),
        "source_sha256": hashlib.sha256(source.encode("utf-8")).hexdigest(),
        "passed": True,
        "source_executed": False,
    }


def _self_format_guard(root: Path) -> dict[str, Any]:
    raw = (root / MODULE).read_bytes()
    _require(b"\r" not in raw, "new source must use LF line endings")
    _require(
        raw.endswith(b"\n") and not raw.endswith(b"\n\n"),
        "new source must have one final LF",
    )
    tree = ast.parse(raw.decode("utf-8"), filename=MODULE)
    return {
        "passed": True,
        "LF": True,
        "exactly_one_EOF_LF": True,
        "AST_dump_sha256": _digest(ast.dump(tree, include_attributes=False)),
        "source_executed": False,
    }


def _parent_pin_map(parent: Mapping[str, Any]) -> dict[str, str]:
    table: dict[str, str] = {PARENT: PARENT_SHA256}
    for key in ("source_pins", "input_pins", "input_hashes", "implementation_hashes"):
        value = parent.get(key)
        if isinstance(value, Mapping):
            for path, item in value.items():
                if not isinstance(path, str):
                    continue
                digest = item.get("expected_sha256") if isinstance(item, Mapping) else item
                if isinstance(digest, str):
                    table.setdefault(path, digest)
                    table.setdefault(Path(path).name, digest)
    nested = parent.get("actual_switch_post_pressure_source")
    if isinstance(nested, Mapping):
        for key in ("source_pins", "input_pins"):
            value = nested.get(key)
            if isinstance(value, Mapping):
                for path, item in value.items():
                    digest = item.get("expected_sha256") if isinstance(item, Mapping) else item
                    if isinstance(path, str) and isinstance(digest, str):
                        table.setdefault(path, digest)
                        table.setdefault(Path(path).name, digest)
    return table


def _pin_paths(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    inherited = _parent_pin_map(parent)
    paths = (
        PARENT,
        Q4_SOURCE,
        Q4_ARCHIVE,
        Q4_CHECK,
        LOOP_SOURCE,
        LOOP_ARCHIVE,
        LOOP_CHECK,
        R100_SOURCE,
        R100_ARCHIVE,
        R100_CHECK,
        BRIDGE_SOURCE,
        DENSITY_SOURCE,
        GRAPH_HELPER,
        SIGMA_SOURCE,
        SIGMA_JETS_SOURCE,
        LONG_ENDPOINT_SOURCE,
        LONG_ENDPOINT_ARCHIVE,
    )
    pins: dict[str, Any] = {}
    for path in paths:
        expected = PARENT_SHA256 if path == PARENT else (
            inherited.get(path) or inherited.get(Path(path).name)
        )
        _require(isinstance(expected, str), f"missing transitive source pin: {path}")
        actual = _sha(root / path)
        _require(actual == expected, f"changed source pin: {path}")
        pins[path] = {
            "path": path,
            "expected_sha256": expected,
            "passed": True,
        }
    pins[MODULE] = {
        "path": MODULE,
        "expected_sha256": _sha(root / MODULE),
        "passed": True,
        "self_pin": True,
    }
    return pins


def _bind_background_packet(
    root: Path, parent: Mapping[str, Any], packet: Mapping[str, Any] | None
) -> dict[str, Any]:
    if packet is None:
        return {
            "provided": False,
            "finite_N_corrected": False,
            "source_executed": False,
        }
    _require(isinstance(packet, Mapping), "background_packet must be a mapping")
    for key in (
        "finite_N_corrected",
        "finiteN_corrected",
        "corrected_by_finite_N",
        "numeric_correction_installed",
    ):
        _require(packet.get(key) is not True, "background packet is finite-N corrected")
    inherited = _parent_pin_map(parent)
    packet_pins = packet.get("source_pins", {})
    checked: dict[str, Any] = {}
    if isinstance(packet_pins, Mapping):
        for path, item in packet_pins.items():
            if not isinstance(path, str):
                continue
            expected = item.get("expected_sha256") if isinstance(item, Mapping) else item
            _require(isinstance(expected, str), f"background pin missing: {path}")
            inherited_expected = inherited.get(path) or inherited.get(Path(path).name)
            if inherited_expected is not None:
                _require(expected == inherited_expected, f"background pin changed: {path}")
            actual = _sha(root / path)
            _require(actual == expected, f"changed background source pin: {path}")
            checked[path] = {
                "path": path,
                "expected_sha256": expected,
                "passed": True,
            }
    return {
        "provided": True,
        "source_pins": checked,
        "finite_N_corrected": False,
        "numeric_values_consumed": False,
        "source_executed": False,
    }


def _validate_q4(root: Path) -> dict[str, Any]:
    archive = _json_gzip(root / Q4_ARCHIVE)
    nodes = archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list) and len(nodes) == Q4_NODE_COUNT,
             "Q4 graph node count changed")
    bindings = archive.get("modified_Q4_pressure_C3_bindings")
    _require(isinstance(bindings, Mapping), "Q4 bindings missing")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping),
             "Q4 charts/endpoints missing")
    receipts: dict[str, Any] = {}
    for name, cfg in _CHARTS.items():
        chart = charts.get(name)
        endpoint = endpoints.get(name)
        _require(isinstance(chart, Mapping) and isinstance(endpoint, Mapping),
                 f"Q4 chart missing: {name}")
        physical = chart.get("physical", {})
        _require(isinstance(physical, Mapping), f"Q4 physical ports missing: {name}")
        for field in ("U", "U_y", "m", "m_y", "Q", "Q_y"):
            _require(physical.get(field) == list(cfg[field]),
                     f"Q4 {name} {field} changed")
        _require(
            chart.get("finite_boundary", {}).get("R_Q_y") == list(cfg["finite_RQy"]),
            f"Q4 {name} finite R_Q_y changed",
        )
        _require(
            chart.get("local_integrals", {}).get("S1") == list(cfg["local_S1"]),
            f"Q4 {name} local S1 changed",
        )
        _require(
            chart.get("cumulative_integrals", {}).get("S1")
            == list(cfg["cumulative_S1"]),
            f"Q4 {name} cumulative S1 changed",
        )
        _require(
            chart.get("incoming_integrals", {}).get("S1")
            == list(cfg["incoming_S1"]),
            f"Q4 {name} incoming S1 changed",
        )
        _require(
            endpoint.get("boundary", {}).get("R_Q_y") == list(cfg["upper_RQy"]),
            f"Q4 {name} upper R_Q_y changed",
        )
        for order, handle in enumerate(cfg["local_S1"][:2]):
            row = _row(nodes, handle)
            _require(
                row.get("operation") == "definite_integral"
                and row.get("original_chart") == name
                and row.get("physical_radius_power") == 2
                and row.get("physical_Jacobian_applied_once") is True
                and row.get("ordinary_slow_Z_derivative_order") == order,
                f"Q4 {name} local S1 source changed",
            )
        for order, handle in enumerate(cfg["finite_RQy"][:2]):
            row = _row(nodes, handle)
            _require(
                row.get("operation") == "product"
                and row.get("arguments") == [cfg["Q_y"][order], cfg["radius"]],
                f"Q4 {name} finite R_Q_y product changed",
            )
        receipts[name] = {
            "physical_ports": {
                key: list(cfg[key])
                for key in ("U", "U_y", "m", "m_y", "Q", "Q_y")
            },
            "local_S1": list(cfg["local_S1"][:2]),
            "cumulative_S1": list(cfg["cumulative_S1"][:2]),
            "finite_RQy": list(cfg["finite_RQy"][:2]),
            "upper_RQy": list(cfg["upper_RQy"][:2]),
        }
    nrow = _row(nodes, N_HANDLE)
    _require(
        nrow.get("operation") == "exact_positive_integer_expression"
        and nrow.get("Z_independent") is True,
        "Q4 common N changed",
    )
    pstar = _row(nodes, PSTAR_HANDLE)
    _require(
        pstar.get("operation") == "analytic_unary"
        and pstar.get("name") == "exp"
        and pstar.get("argument") == PSTAR_LOG_ARGUMENT_HANDLE,
        "Q4 Pstar source changed",
    )
    return {
        "passed": True,
        "node_count": len(nodes),
        "charts": receipts,
        "common_N": N_HANDLE,
        "Pstar": PSTAR_HANDLE,
        "Pstar_log_argument": PSTAR_LOG_ARGUMENT_HANDLE,
    }


def _validate_loop(root: Path) -> dict[str, Any]:
    archive = _json_gzip(root / LOOP_ARCHIVE)
    bindings = archive.get("prescribed_loop_C5_bindings")
    _require(isinstance(bindings, Mapping), "C5 bindings missing")
    _require(
        bindings.get("N") == N_HANDLE
        and bindings.get("eta") == 21
        and bindings.get("dstar") == 22
        and bindings.get("sc") == 10
        and bindings.get("phase_origin") == 24954,
        "C5 common controls changed",
    )
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping),
             "C5 charts/endpoints missing")
    nodes = archive.get("exact_graph_nodes")
    _require(isinstance(nodes, list), "C5 graph missing")
    checks: dict[str, Any] = {}
    for name, cfg in _CHARTS.items():
        chart = charts.get(name)
        _require(isinstance(chart, Mapping), f"C5 chart missing: {name}")
        _require(
            chart.get("actual_phase") == cfg["phase"]
            and chart.get("global_offset") == cfg["global_offset"]
            and chart.get("coordinate") == cfg["coordinate"]
            and chart.get("domain") == [0, 1]
            and chart.get("radius") == cfg["radius"]
            and chart.get("physical_measure") == cfg["J"],
            f"C5 {name} geometry changed",
        )
        _require(
            chart.get("incoming", {}).get("m") == list(cfg["incoming_m"])
            and chart.get("history", {}).get("m") == list(cfg["history_m"]),
            f"C5 {name} m history changed",
        )
        for memory_name in ("m", "h", "k", "e", "p"):
            incoming = chart.get("incoming", {}).get(memory_name)
            history = chart.get("history", {}).get(memory_name)
            memory = chart.get("memory", {}).get(memory_name)
            _require(
                isinstance(incoming, list)
                and len(incoming) == 6
                and all(type(handle) is int for handle in incoming)
                and isinstance(history, list)
                and len(history) == 6
                and all(type(handle) is int for handle in history)
                and type(memory) is int,
                f"C5 {name} {memory_name} correction ports changed",
            )
        _require(
            chart.get("finite_N", {}).get("E_N") == chart.get("modified_fields", {}).get("E")
            and chart.get("finite_N", {}).get("V_N") == chart.get("modified_fields", {}).get("V"),
            f"C5 {name} finite corrected fields changed",
        )
        _require(chart.get("memory", {}).get("p") == 1,
                 f"C5 {name} pressure memory was reset")
        phase_row = _row(nodes, cfg["phase"])
        phase_argument = phase_row.get("argument")
        _require(
            phase_row.get("operation") == "analytic_unary"
            and phase_row.get("name") == "fractional_part"
            and type(phase_argument) is int,
            f"C5 {name} phase source changed",
        )
        product = _row(nodes, phase_argument)
        _require(
            product.get("operation") == "product"
            and product.get("arguments") == [N_HANDLE, cfg["global_offset"]],
            f"C5 {name} phase is not common-N/global-offset",
        )
        _require(not _has_dependency(nodes, cfg["global_offset"], 20),
                 f"C5 {name} phase depends on ordinary Z")
        checks[name] = {
            "phase": cfg["phase"],
            "global_offset": cfg["global_offset"],
            "coordinate": cfg["coordinate"],
            "radius": cfg["radius"],
            "native_J": cfg["J"],
            "incoming_m": list(cfg["incoming_m"]),
            "history_m": list(cfg["history_m"]),
            "all_five_correction_memories_retained": True,
            "all_five_incoming_ports": True,
            "all_five_history_ports": True,
            "all_five_memory_ports": True,
            "phase_Z_independent": True,
        }
    _require(
        endpoints.get("post_power", {}).get("m") == charts["long_reshape"]["incoming"]["m"]
        and endpoints.get("long_reshape", {}).get("m")
        == charts["reference"]["incoming"]["m"],
        "C5 long/reference incoming m pointers reset",
    )
    for predecessor, successor in (("post_power", "long_reshape"), ("long_reshape", "reference")):
        for memory_name in ("m", "h", "k", "e", "p"):
            _require(
                endpoints.get(predecessor, {}).get(memory_name)
                == charts[successor].get("incoming", {}).get(memory_name),
                f"C5 {predecessor}->{successor} {memory_name} inlet reset",
            )
    _require(
        charts["long_reshape"].get("incoming") == endpoints.get("post_power")
        and charts["reference"].get("incoming") == endpoints.get("long_reshape"),
        "C5 long/reference full incoming histories are not predecessor endpoints",
    )
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
        "eta": 21,
        "dstar": 22,
        "sc": 10,
        "phase_origin": 24954,
        "charts": checks,
        "native_memories": ["m", "h", "k", "e", "p"],
        "all_five_correction_incoming_history_memory_ports": True,
        "nonreset_incoming_histories": True,
        "phase_Z_independent": True,
        "finite_N_correction_source": "serialized C5 loop only",
    }


def _validate_outer(root: Path) -> dict[str, Any]:
    archive = _json_gzip(root / R100_ARCHIVE)
    bindings = archive.get("prescribed_R100_Rm_C6_bindings")
    _require(isinstance(bindings, Mapping), "R100 bindings missing")
    charts = bindings.get("charts")
    endpoints = bindings.get("endpoints")
    _require(isinstance(charts, Mapping) and isinstance(endpoints, Mapping),
             "R100 charts/endpoints missing")
    expected = {
        "post_switch": (78472, 91107, 91103),
        "long_reshape": (79966, 79970, 7),
        "reference": (81520, 81524, 81519),
    }
    checks: dict[str, Any] = {}
    for name, (phase, radius, measure) in expected.items():
        chart = charts.get(name)
        _require(
            isinstance(chart, Mapping)
            and chart.get("phase") == phase
            and chart.get("radius") == radius
            and chart.get("physical_measure") == measure
            and chart.get("domain") == [0, 1],
            f"R100 {name} geometry changed",
        )
        actual = chart.get("actual")
        radial = chart.get("radial")
        incoming = chart.get("incoming_histories")
        moments = chart.get("moments")
        _require(
            isinstance(actual, Mapping)
            and isinstance(radial, Mapping)
            and isinstance(incoming, Mapping)
            and isinstance(moments, Mapping),
            f"R100 {name} source ports missing",
        )
        for field in ("phi", "U", "phi_y", "U_y"):
            values = actual.get(field)
            _require(
                isinstance(values, list)
                and len(values) == 7
                and all(type(value) is int for value in values),
                f"R100 {name} actual {field} ports changed",
            )
        _require(
            all(value == 0 for value in actual["U_y"])
            and all(value == 0 for value in radial.get("V_y", [])),
            f"R100 {name} actual U ports changed",
        )
        for moment_name in ("H", "M", "K", "A", "B", "C"):
            for port_name, ports in (
                ("incoming", incoming.get(moment_name)),
                ("moments", moments.get(moment_name)),
            ):
                _require(
                    isinstance(ports, list)
                    and len(ports) == 7
                    and all(type(value) is int for value in ports),
                    f"R100 {name} {port_name} {moment_name} ports changed",
                )
        checks[name] = {
            "phase": phase,
            "radius": radius,
            "physical_measure": measure,
            "actual_phi": list(actual["phi"]),
            "actual_U": list(actual["U"]),
            "actual_phi_y": list(actual["phi_y"]),
            "actual_U_y": list(actual["U_y"]),
            "radial_V_y": list(chart["radial"]["V_y"]),
            "incoming_native_moments": {
                key: list(incoming[key]) for key in ("H", "M", "K", "A", "B", "C")
            },
            "native_moments": {
                key: list(moments[key]) for key in ("H", "M", "K", "A", "B", "C")
            },
        }
    _require(
        charts["long_reshape"].get("incoming_histories")
        == endpoints["R110"].get("moments")
        and charts["reference"].get("incoming_histories")
        == endpoints["Rsh"].get("moments"),
        "R100 long/reference histories are not predecessor-bound",
    )
    for chart_name, endpoint_name in (("long_reshape", "R110"), ("reference", "Rsh")):
        incoming = charts[chart_name]["incoming_histories"]
        endpoint_moments = endpoints[endpoint_name]["moments"]
        endpoint_actual = endpoints[endpoint_name].get("actual")
        _require(isinstance(endpoint_actual, Mapping), f"R100 {endpoint_name} actual source missing")
        for field in ("phi", "U", "phi_y", "U_y"):
            values = endpoint_actual.get(field)
            _require(
                isinstance(values, list)
                and len(values) == 7
                and all(type(value) is int for value in values),
                f"R100 {endpoint_name} endpoint {field} changed",
            )
        for moment_name in ("H", "M", "K", "A", "B", "C"):
            _require(
                isinstance(incoming.get(moment_name), list)
                and len(incoming[moment_name]) == 7
                and all(type(handle) is int for handle in incoming[moment_name])
                and incoming[moment_name] == endpoint_moments.get(moment_name)
                and all(
                    isinstance(endpoints[endpoint_name].get("actual", {}).get(field), list)
                    and len(endpoints[endpoint_name]["actual"][field]) == 7
                    for field in ("phi", "U", "phi_y", "U_y")
                ),
                f"R100 {chart_name} six native moment {moment_name} binding changed",
            )
        checks[chart_name]["predecessor_endpoint"] = endpoint_name
        checks[chart_name]["predecessor_actual"] = {
            field: list(endpoint_actual[field])
            for field in ("phi", "U", "phi_y", "U_y")
        }
        checks[chart_name]["predecessor_moments"] = {
            key: list(endpoint_moments[key])
            for key in ("H", "M", "K", "A", "B", "C")
        }
    definition = bindings.get("definition", {})
    for key in (
        "original_B_from_actual_R110_and_anchored_G",
        "all_six_nonzero_incoming_and_every_chart_memory_retained",
        "original_independent_P0_and_physical_Pstar_F0_factors_once",
        "finite_N_modified_fields_correction_histories_and_all17_not_installed",
    ):
        _require(definition.get(key) is True, f"R100 source premise missing: {key}")
    _require(
        bindings.get("Pstar") == PSTAR_HANDLE
        and bindings.get("T") == 7
        and bindings.get("reference_length") == 81519,
        "R100 native T/Pstar/reference length changed",
    )
    geometry = bindings.get("geometry_certificate", {})
    _require(
        geometry.get("actual_native_T_and_logC_match_common_source") is True
        and geometry.get("complete_radius_order")
        == "100<R1<R2<110<Rsh<Rz<Rrestore<Rm",
        "R100 geometry certificate changed",
    )
    return {
        "passed": True,
        "Pstar": PSTAR_HANDLE,
        "T": 7,
        "long_reference": checks,
        "B_definition": "actual_R110_phi, anchored G, log(220)/2, log(1+Z^2)",
        "phi_from_Utheta": "Utheta/(F0*sqrt(2R)); F0 is Z-only",
        "U_constant": True,
        "U_y_zero": True,
        "nonreset_incoming_histories": True,
        "all_six_native_incoming_moments": True,
        "all_six_native_endpoint_moments": True,
        "physical_slow_J_bound_once": True,
        "finite_N_correction_histories_are_C5_only": True,
    }


def _rehydrate_parent(
    parent: Mapping[str, Any], q4: Mapping[str, Any]
) -> tuple[GraphPairs, dict[str, Any]]:
    graph = GraphPairs(q4["exact_graph_nodes"])
    merged = parent["actual_micro_macro_correlated_m15"]["merged_namespace"]
    seam = parent["actual_macro_switch_RQy_seam"]["graph_patch"]
    first = parent["actual_macro_first_switch_pressure_source"]["graph_patch"]
    post = parent["actual_switch_post_pressure_source"]["graph_patch"]
    _require(merged.get("base_node_count") == Q4_NODE_COUNT, "merged Q4 base changed")
    _require(len(merged.get("appended_nodes", [])) == 987, "merged suffix count changed")
    _require(seam.get("base_node_count") == 184687, "accepted seam base changed")
    _require(len(seam.get("appended_nodes", [])) == 16, "accepted seam suffix changed")
    _require(first.get("base_node_count") == 184703, "accepted first base changed")
    _require(len(first.get("appended_nodes", [])) == 195, "accepted first suffix changed")
    _require(post.get("base_node_count") == 184898, "accepted post base changed")
    _require(post.get("full_node_count") == PARENT_NODE_COUNT, "accepted parent node count changed")
    for rows in (
        merged["appended_nodes"],
        seam["appended_nodes"],
        first["appended_nodes"],
        post["appended_nodes"],
    ):
        for row in rows:
            _append_saved(graph, row)
    _require(len(graph.nodes) == PARENT_NODE_COUNT, "parent graph did not rehydrate")
    _require(graph.nodes[:Q4_NODE_COUNT] == q4["exact_graph_nodes"],
             "Q4 graph prefix changed")
    expected_digest = post.get("full_canonical_nodes_sha256")
    _require(isinstance(expected_digest, str) and _digest(graph.nodes) == expected_digest,
             "parent canonical full prefix hash changed")
    _require(
        expected_digest
        == parent["actual_switch_post_pressure_source"]["graph_patch"][
            "full_canonical_nodes_sha256"
        ],
        "parent serialized full digest mismatch",
    )
    return graph, {
        "q4_node_count": Q4_NODE_COUNT,
        "merged_node_count": 184687,
        "seam_node_count": 184703,
        "first_node_count": 184898,
        "parent_node_count": PARENT_NODE_COUNT,
        "segment_starts": {
            "q4": 0,
            "accepted_merged": Q4_NODE_COUNT,
            "accepted_seam": 184687,
            "accepted_first": 184703,
            "accepted_post": 184898,
            "new_suffix": PARENT_NODE_COUNT,
        },
        "suffix_counts": {
            "accepted_merged": 987,
            "accepted_seam": 16,
            "accepted_first": 195,
            "accepted_post": 410,
        },
        "full_prefix_snapshot": copy.deepcopy(graph.nodes),
        "full_prefix_canonical_nodes_sha256": _digest(graph.nodes),
        "serialized_rehydration": True,
    }


def _validate_parent_scope(parent: Mapping[str, Any]) -> dict[str, Any]:
    source = parent.get("actual_switch_post_pressure_source")
    _require(isinstance(source, Mapping), "accepted switch/post source missing")
    limits = source.get("limits")
    _require(
        isinstance(limits, Mapping)
        and limits.get("old_212_terms_and_errors_unchanged") is True,
        "accepted old 212 terms/errors preservation receipt missing",
    )
    return {
        "passed": True,
        "old_212_terms_and_errors_unchanged": True,
        "older_germ_unproved_status_retained": True,
        "older_germ_used_as_proof": False,
        "accepted_numeric_computations_replayed": False,
    }


def _build_six_source(
    graph: GraphPairs, cfg: Mapping[str, Any]
) -> tuple[dict[str, list[int]], dict[str, Any]]:
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
    }


def _build_reduced_D(
    graph: GraphPairs, integrals: Mapping[str, Sequence[int]]
) -> tuple[list[int], dict[str, Any]]:
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
        graph.pdiv(graph.pmul(h, z), L, "positive physical L=1-delta*Z^2"),
        half,
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
        "two_hZ2_over_L2": graph.pmul(
            two_hZ2_over_L2, _pair(integrals["U_squared"], "U_squared")
        ),
        "two_Zd_over_L2": graph.pmul(
            two_Zd_over_L2, _pair(integrals["U_U_Z"], "U_U_Z")
        ),
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
    }


def _build_reduced_B(
    graph: GraphPairs, cfg: Mapping[str, Any]
) -> tuple[list[int], dict[str, Any]]:
    z = _pair((20, 1), "Z")
    u = _pair(cfg["U"][:2], "U")
    q = _pair(cfg["Q"][:2], "Q")
    L = _pair((184439, 184438), "L")
    L2 = _pair((184440, 184443), "L^2")
    radius_squared = graph.mul(cfg["radius"], cfg["radius"])
    first = graph.pscale(
        graph.pdiv(q, L, "positive physical L=1-delta*Z^2"),
        radius_squared,
    )
    second = graph.pscale(
        graph.pscale(graph.pmul(q, q), radius_squared),
        28,
    )
    third = graph.pscale(
        graph.pmul(
            graph.pmul(z, z),
            graph.fixed(radius_squared),
            graph.pmul(u, u),
        ),
        48,
    )
    third = graph.pneg(
        graph.pdiv(third, L2, "positive physical L=1-delta*Z^2")
    )
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
    }


def _build_join(
    graph: GraphPairs,
    local_s1: Mapping[str, Sequence[int]],
    lowers: Mapping[str, Sequence[int]],
    formal: Mapping[str, Any],
) -> dict[str, Any]:
    """Append the signed formula without assuming either new seam identity."""
    _require(formal.get("passed") is False, "new seam proof is still open")
    half = {
        name: [graph.mul(int(value), 28) for value in local_s1[name]]
        for name in _CHARTS
    }
    increments = {
        name: [
            graph.add(half[name][k], int(lowers[name][k]),
                      graph.neg(_CHARTS[name]["upper_RQy"][k]))
            for k in range(2)
        ]
        for name in _CHARTS
    }
    long = [graph.add(PARENT_JOINED[k], increments["long_reshape"][k]) for k in range(2)]
    reference = [graph.add(long[k], increments["reference"][k]) for k in range(2)]
    return {
        "passed": True,
        "ordinary_Z_orders": [0, 1],
        "parent_joined_contracted": list(PARENT_JOINED),
        "long_reshape_local_S1": list(local_s1["long_reshape"]),
        "reference_local_S1": list(local_s1["reference"]),
        "long_reshape_half_S1": half["long_reshape"],
        "reference_half_S1": half["reference"],
        "long_reshape_lower_RQy": list(lowers["long_reshape"]),
        "reference_lower_RQy": list(lowers["reference"]),
        "post_power_upper_RQy": list(POST_UPPER_RQY),
        "long_reshape_upper_RQy": list(_CHARTS["long_reshape"]["upper_RQy"][:2]),
        "reference_upper_RQy": list(_CHARTS["reference"]["upper_RQy"][:2]),
        "increments_uncontracted": increments,
        "joined_after_long_uncontracted": long,
        "joined_after_reference_uncontracted": reference,
        "joined_uncontracted": reference,
        "contracted_after_long_raw": None,
        "contracted_after_long": None,
        "contracted_before_linearization": None,
        "joined_contracted": None,
        "linearization_receipt": {
            "passed": False,
            "contraction_allowed": False,
            "uncontracted_source_preserved": True,
            "standalone_caps_preserved": True,
            "reason": "actual modified physical Uy/Qy seam congruence remains unproved",
        },
        "formula_uncontracted": (
            "parent_joined + long_local_S1/2 + long_R_Q_y(0) - long_R_Q_y(1) "
            "+ reference_local_S1/2 + reference_R_Q_y(0) - reference_R_Q_y(1)"
        ),
        "formula_contracted": None,
        "local_S1_added_once": True,
        "cumulative_S1_not_used": True,
        "outer_RQy_and_caps_retained": True,
        "standalone_endpoint_atoms_retained": True,
        "signed_relation_not_global_closure": True,
        "contraction_status": "BLOCKED_UNPROVEN_SOURCE_CONGRUENCE",
        "contraction_allowed": False,
        "endpoint_congruence_receipt": copy.deepcopy(formal),
    }


def _unproved_endpoint_receipts(
    root: Path,
    q4: Mapping[str, Any],
    loop: Mapping[str, Any],
    outer: Mapping[str, Any],
    outer_check: Mapping[str, Any],
    loop_check: Mapping[str, Any],
) -> dict[str, Any]:
    outer_bindings = outer["prescribed_R100_Rm_C6_bindings"]
    loop_bindings = loop["prescribed_loop_C5_bindings"]
    q4_bindings = q4["modified_Q4_pressure_C3_bindings"]
    source_receipts = {
        "long_reference_constructor": _source_fragment(
            root / R100_SOURCE,
            "__init__",
            (
                "incoming=self.endpoints['R110'];self.B=B=",
                "a.scale(a.trim(bridge_owner.core7['G']),g.neg(self.core.parameters['Lambda']))",
                "a.log(incoming['actual']['phi'],'positive actual R110 phi')",
                "self.logq=logq=a.log(a.add(a.fixed(1),a.mul(z,z)),'positive original 1+Z^2')",
                "a.fixed(g.mul(c('1/2'),g.unary('log',c(220))))",
                "logq)",
                "self.F0=F0=a.trim(bridge_owner.core7['F0'])",
                "q=g.symbol('true_prescribed_long_reshape_fraction')",
                "logUtheta=a.add(a.scale(B,g.sub(g.one,sigma(g,q)))",
                "slope=a.add(a.fixed('-2/5'),a.scale(B,g.neg(g.quotient(sigma_prime(g,q),T",
                "U_y=a.fixed(0)",
                "self.endpoints['Rsh']=self.end('long_reshape',g.one)",
                "logUtheta=a.add(a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
                "U_y=a.fixed(0)",
                "self.endpoints['Rz']=self.end('reference',g.one)",
            ),
        ),
        "phi_from_Utheta": _source_fragment(
            root / R100_SOURCE,
            "phi_from_Utheta",
            (
                "return a.div(Utheta,a.scale(self.F0,g.unary('sqrt',g.mul(g.constant(2),R)))",
                "same positive sqrt(2R)*F0; F0 is Z-only",
            ),
        ),
        "bridge_generic_fields": _source_fragment(
            root / BRIDGE_SOURCE,
            "generic_fields",
            (
                "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)",
                "radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y'])),V_y=a.scale(actual['U_y'],inverse))",
            ),
        ),
        "bridge_densities": _source_fragment(
            root / BRIDGE_SOURCE,
            "densities",
            (
                "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2),",
                "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi))",
            ),
        ),
        "bridge_histories": _source_fragment(
            root / BRIDGE_SOURCE,
            "histories",
            (
                "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase",
                "original_incoming_memory_retained=True",
                "source_measure_applied_once=True",
            ),
        ),
        "q4_physical_ports": _source_fragment(
            root / Q4_SOURCE,
            "physical_pressure_fields",
            (
                "U=a.scale(fields['V'],Pstar)",
                "Uy=a.scale(radial['V_y'],Pstar)",
                "m=a.scale(fields['m'],Pstar)",
                "my=a.sub(U,m)",
                "Q=div(numerator(U,m))",
                "Qy=div(numerator(Uy,my))",
            ),
        ),
        "q4_native_derivative": _qualified_fragment(
            root / Q4_SOURCE,
            "NativeDerivative",
            "__init__",
            (
                "self.cache={phase.node:g.zero};self.core_extensions=[]",
            ),
        ),
        "q4_native_derivative_rules": _qualified_fragment(
            root / Q4_SOURCE,
            "NativeDerivative",
            "derivative",
            (
                "elif op=='sum':value=g.add(*(D(j) for j in row['arguments']))",
                "elif op=='product':",
                "elif op=='positive_quotient':",
                "elif op in ('analytic_unary','exact_original_source_flat_sigma'):",
                "elif op=='function_substitution':",
                "elif op=='definite_integral':",
                "elif op=='substitute_original_inverse_angle':",
                "elif op=='original_monotone_phase_inverse':",
                "value=g.quotient(g.sub(self.derivative(target),at(held)),lam,'original positive phase lambda; exact implicit native derivative')",
                "elif op=='original_lazy_flat_branch':",
            ),
        ),
        "q4_total_radial_functions": _source_fragment(
            root / Q4_SOURCE,
            "total_radial_functions",
            (
                "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref))",
                "xy=a.add(a.scale(slow['A'],invN),chart['loop']['A_phi'])",
                "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
                "total_fast_phase_terms_included=True",
                "physical_y_equals_native_partial_divided_by_Jacobian=True",
            ),
        ),
        "c5_transport": _source_fragment(
            root / LOOP_SOURCE,
            "transport",
            (
                "width=g.mul(jacobian,g.sub(coordinate,lower))",
                "contribution=a.integral(body,t,lower,coordinate,actual_finite_N_signed_correction=True",
                "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
                "physical_Jacobian_applied_once=True",
                "zero_correction_inlet_only_at_original_rminus=True",
                "predecessor_memory_never_reset=True",
            ),
        ),
        "c5_loop_functions": _source_fragment(
            root / LOOP_SOURCE,
            "loop_functions",
            (
                "t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))",
                "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
                "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
                "T2=FunctionJet(tuple(g.mul(c(math.factorial(k)),integrate(t2[k,0])) for k in range(6)))",
                "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
                "inverse=g.node('original_monotone_phase_inverse'",
                "derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction'",
                "partials={(k,0):g.mul(c(Fraction(1,math.factorial(k))),rows(Phi)[k])",
                "chi=a.sub(a.fixed(phase),a.scale(psi,g.quotient(g.one,twopi,'positive 2*pi')))",
                "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),a.scale(roots['b'],phase))",
                "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
                "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
                "Aphi=lazy(",
                "Bphi=lazy(",
            ),
        ),
        "sigma_flat": _source_fragment(
            root / R100_SOURCE,
            "sigma_prime",
            (
                "name='original_flat_sigma_prime'",
                "exact_flat_source_derivative_not_a_bound=True",
            ),
        ),
        "smooth_sigma_source": _source_fragment(
            root / SIGMA_SOURCE,
            "stable_sigma",
            (
                "ifv<=0:returnc.mpf(0),c.mpf(0)",
                "ifv>=1:returnc.mpf(1),c.mpf(0)",
            ),
        ),
        "smooth_sigma_jets_source": _source_fragment(
            root / SIGMA_JETS_SOURCE,
            "sigma_jets",
            (
                "iflo<=0:rows.append(IntervalTaylor.constant(c,0,4))",
                "ifhi>=1:rows.append(IntervalTaylor.constant(c,1,4))",
            ),
        ),
        "anchor_identity": _source_fragment(
            root / LONG_ENDPOINT_SOURCE,
            "source_bindings",
            (
                "'B':'G*(-Lambda)+logarithm(phi)+logarithm(1+z*z)+c.ln(2*R)/2'",
                "phase_zero_function_identity='B-logq-logCstar=-Lambda*G+log(phi110)+log(220)/2=log(Utheta110)'",
                "identity_uses_anchored_F0_definition_not_interval_overlap=True",
            ),
        ),
    }
    for name, receipt in source_receipts.items():
        _require(receipt.get("passed") is True, f"source receipt failed: {name}")

    _require(outer_bindings.get("Pstar") == PSTAR_HANDLE, "R100 Pstar binding changed")
    _require(loop_bindings.get("N") == N_HANDLE, "C5 N binding changed")
    _require(
        q4_bindings["charts"]["long_reshape"]["physical"]["U"]
        == list(_CHARTS["long_reshape"]["U"])
        and q4_bindings["charts"]["reference"]["physical"]["Q_y"]
        == list(_CHARTS["reference"]["Q_y"]),
        "Q4 endpoint port binding changed",
    )
    gap = (
        "Source/cache port bindings do not prove modified physical derivative "
        "congruence. Need actual endpoint C2 controls and all six native/five "
        "correction histories, node-by-node held-phase NativeDerivative rules "
        "with unequal native J, and actual periodic loop derivative identities."
    )
    unproved = {"passed": False, "reason": gap}
    seams = {
        "post_power_upper_to_long_lower": {
            "left": list(POST_UPPER_RQY),
            "right_source": list(_CHARTS["long_reshape"]["finite_RQy"][:2]),
            "finite_N_Qy_identity_proven": False,
            "contraction_allowed": False,
        },
        "long_upper_to_reference_lower": {
            "left": list(_CHARTS["long_reshape"]["upper_RQy"][:2]),
            "right_source": list(_CHARTS["reference"]["finite_RQy"][:2]),
            "finite_N_Qy_identity_proven": False,
            "contraction_allowed": False,
        },
    }
    return {
        "passed": False,
        "status": "SOURCE_BOUND_PORTS_ONLY_SEAM_PROOF_OPEN",
        "ordinary_Z_orders": [0, 1],
        "common_N": N_HANDLE,
        "source_formulas": source_receipts,
        "endpoint_identity": dict(unproved),
        "rational_same_function": dict(unproved),
        "control_naturality": dict(unproved),
        "native_derivative_naturality": dict(unproved),
        "history_two_jet": dict(unproved),
        "periodic_phase_lemma": dict(unproved),
        "smooth_sigma_flatness": {
            "passed": False,
            "source_bound": True,
            "reason": "source token bindings are recorded; no seam theorem is accepted here",
        },
        "q4_physical_port_bindings": {
            name: {
                key: list(q4_bindings["charts"][name]["physical"][key])
                for key in ("U", "U_y", "m", "m_y", "Q", "Q_y")
            }
            for name in _CHARTS
        },
        "c5_corrected_port_bindings": {
            name: {
                "actual_phase": loop_bindings["charts"][name]["actual_phase"],
                "incoming": copy.deepcopy(loop_bindings["charts"][name]["incoming"]),
                "history": copy.deepcopy(loop_bindings["charts"][name]["history"]),
                "memory": copy.deepcopy(loop_bindings["charts"][name]["memory"]),
            }
            for name in _CHARTS
        },
        "port_and_history_bindings_validated": {
            "all_six_native_incoming_moments": outer_check["all_six_native_incoming_moments"],
            "all_six_native_endpoint_moments": outer_check["all_six_native_endpoint_moments"],
            "all_five_correction_incoming_history_memory_ports":
                loop_check["all_five_correction_incoming_history_memory_ports"],
            "bindings_are_not_endpoint_identity_proofs": True,
        },
        "modified_Uy_Qy_closed": False,
        "natural_B_y_closed": False,
        "seams": seams,
        "source_values_selected": False,
        "interval_or_cap_values_used": False,
        "independent_review_rejected_generic_chain_rule_candidate": True,
    }


def build(
    root: str | Path,
    parent: Mapping[str, Any],
    background_packet: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Build the append-only current-N long/reference pressure source."""

    root_path = Path(root)
    _require(isinstance(parent, Mapping), "parent report mapping required")
    _require(_sha(root_path / PARENT) == PARENT_SHA256,
             "current switch/post parent file changed")
    saved_parent = _json(root_path / PARENT)
    _require(_digest(saved_parent) == _digest(parent),
             "exact current switch/post parent required")

    source_pins = _pin_paths(root_path, parent)
    background = _bind_background_packet(root_path, parent, background_packet)
    q4 = _json_gzip(root_path / Q4_ARCHIVE)
    loop = _json_gzip(root_path / LOOP_ARCHIVE)
    outer = _json_gzip(root_path / R100_ARCHIVE)
    q4_check = _validate_q4(root_path)
    loop_check = _validate_loop(root_path)
    outer_check = _validate_outer(root_path)
    graph, prefix = _rehydrate_parent(parent, q4)
    parent_scope = _validate_parent_scope(parent)
    formal = _unproved_endpoint_receipts(
        root_path, q4, loop, outer, outer_check, loop_check
    )

    chart_output: dict[str, Any] = {}
    local_s1: dict[str, Sequence[int]] = {}
    lower_wrappers: dict[str, Sequence[int]] = {}
    for name, cfg in _CHARTS.items():
        integrals, six = _build_six_source(graph, cfg)
        b_rows, b_source = _build_reduced_B(graph, cfg)
        d_rows, d_source = _build_reduced_D(graph, integrals)
        reduced = list(graph.padd(_pair(b_rows, "B"), _pair(d_rows, "D")))
        local_s1[name] = tuple(cfg["local_S1"][:2])
        lower = graph.at(
            _pair(cfg["finite_RQy"][:2], f"{name} finite R_Q_y"),
            cfg["coordinate"],
            graph.constant(0),
        )
        lower_wrappers[name] = lower
        chart_output[name] = {
            "passed": True,
            "native_chart": cfg["native_chart"],
            "q4_chart": cfg["q4_chart"],
            "common_N": N_HANDLE,
            "N_materialized": False,
            "phase": cfg["phase"],
            "global_offset": cfg["global_offset"],
            "phase_Z_independent": True,
            "coordinate": cfg["coordinate"],
            "domain": [0, 1],
            "radius": cfg["radius"],
            "native_J": cfg["J"],
            "physical_ports": {
                key: list(cfg[key])
                for key in ("U", "U_y", "m", "m_y", "Q", "Q_y")
            },
            "incoming_m": list(cfg["incoming_m"]),
            "history_m": list(cfg["history_m"]),
            "native_memories": ["m", "h", "k", "e", "p"],
            "native_history_moments": ["H", "M", "K", "A", "B", "C"],
            "local_S1": list(cfg["local_S1"][:2]),
            "cumulative_S1": list(cfg["cumulative_S1"][:2]),
            "finite_boundary_RQy": list(cfg["finite_RQy"][:2]),
            "lower_wrapper_t0_or_q0": list(lower),
            "upper_RQy_standalone": list(cfg["upper_RQy"][:2]),
            "six_I1_source": six,
            "integrals": integrals,
            "B_source": b_source,
            "D_source": d_source,
            "S1_source": {
                "rows": reduced,
                "formula": "S1=B+D",
                "D_is_not_full_S1": True,
            },
            "raw_graph_unit": cfg["graph_measure"],
            "Pstar_applied_once": True,
            "physical_R_not_mean_scale": True,
            "background_packet_consumed": False,
        }

    join = _build_join(graph, local_s1, lower_wrappers, formal)
    _require(
        graph.nodes[:PARENT_NODE_COUNT] == prefix["full_prefix_snapshot"],
        "full 185308-node parent prefix changed",
    )
    new_suffix = copy.deepcopy(graph.nodes[PARENT_NODE_COUNT:])
    _require(new_suffix, "new source did not append graph rows")

    graph_patch = {
        "base_archive": Q4_ARCHIVE,
        "base_node_count": PARENT_NODE_COUNT,
        "base_prefix_unchanged": True,
        "base_prefix_canonical_nodes_sha256": prefix[
            "full_prefix_canonical_nodes_sha256"
        ],
        "parent_full_canonical_nodes_sha256": prefix[
            "full_prefix_canonical_nodes_sha256"
        ],
        "new_node_start": PARENT_NODE_COUNT,
        "appended_nodes": new_suffix,
        "appended_node_count": len(new_suffix),
        "full_node_count": len(graph.nodes),
        "full_prefix_digest": prefix["full_prefix_canonical_nodes_sha256"],
        "full_canonical_nodes_sha256": _digest(graph.nodes),
        "segment_starts": prefix["segment_starts"],
        "prefix_suffix_counts": prefix["suffix_counts"],
        "endpoint_handles": {
            "parent_joined": list(PARENT_JOINED),
            "post_power_upper_standalone": list(POST_UPPER_RQY),
            "long_reshape_lower_q0": list(lower_wrappers["long_reshape"]),
            "long_reshape_upper_standalone": list(_CHARTS["long_reshape"]["upper_RQy"][:2]),
            "reference_lower_q0": list(lower_wrappers["reference"]),
            "reference_upper_standalone": list(_CHARTS["reference"]["upper_RQy"][:2]),
            "long_reshape_local_S1": list(local_s1["long_reshape"]),
            "reference_local_S1": list(local_s1["reference"]),
            "joined_after_long_uncontracted": join["joined_after_long_uncontracted"],
            "joined_after_reference_uncontracted": join["joined_after_reference_uncontracted"],
            "contracted_after_long_raw": join["contracted_after_long_raw"],
            "contracted_after_long": join["contracted_after_long"],
            "contracted_before_linearization": join["contracted_before_linearization"],
            "joined_contracted": join["joined_contracted"],
        },
    }
    source_contract = {
        "passed": True,
        "current_N_source_only": True,
        "finite_N_correction_histories_from_C5_loop": True,
        "ordinary_Z_orders": [0, 1],
        "common_N_handle": N_HANDLE,
        "N_not_materialized": True,
        "phase_Z_independent": True,
        "native_J_distinct": {"long_reshape": 7, "reference": 81519},
        "actual_source_chart_mapping": {
            "long_reshape": "long_reshape",
            "reference": "reference",
        },
        "physical_y_derivatives_claimed": False,
        "physical_Pstar_applied_once": True,
        "Pstar_handle": PSTAR_HANDLE,
        "Pstar_log_argument_handle": PSTAR_LOG_ARGUMENT_HANDLE,
        "physical_R_not_a_mean_scale": True,
        "physical_slow_J_applied_once": True,
        "Nfast_phase_common_and_Z_independent": True,
        "modified_Uy_Qy_closed_from_actual_source": formal.get(
            "modified_Uy_Qy_closed"
        ) is True,
        "natural_B_y_closed_from_actual_source": formal.get(
            "natural_B_y_closed"
        ) is True,
        "global_scale_exp_materialized": False,
        "all_native_memories_retained": True,
        "all_C5_correction_memories_retained": True,
        "nonreset_histories": True,
        "parent_joined_entry": list(PARENT_JOINED),
        "outer_post_RQy_retained": True,
        "local_S1_added_once": True,
        "cumulative_S1_not_added": True,
        "contraction_allowed": join["contraction_allowed"],
        "contraction_status": join["contraction_status"],
        "complete_pressure_numeric": False,
        "global_five_moment_closure": False,
        "accepted_ancestors_replayed": False,
        "source_execution": False,
    }
    limits = {
        "endpoint_phi_U_C2_proved": formal["endpoint_identity"].get("passed") is True,
        "seams_proved": all(
            formal["seams"][name].get("finite_N_Qy_identity_proven") is True
            for name in ("post_power_upper_to_long_lower", "long_upper_to_reference_lower")
        ),
        "q4_rational_same_function_proved": formal["rational_same_function"].get("passed") is True,
        "smooth_sigma_flatness_bound": formal["smooth_sigma_flatness"].get("passed") is True,
        "periodic_phase_lemma_proved": formal["periodic_phase_lemma"].get("passed") is True,
        "joined_pressure_source_defined": True,
        "joined_pressure_numeric_enclosure": False,
        "outer_RQy_retained": True,
        "outer_RQy_evaluated": False,
        "all_caps_retained": True,
        "contracted_endpoint_atoms_removed": join["joined_contracted"] is not None,
        "raw_contracted_graph_retained": join["contracted_after_long_raw"] is not None,
        "global_closure": False,
        "complete_pressure_global_closure": False,
        "accepted_ancestors_replayed": False,
        "N_materialized": False,
        "background_finite_N_corrected": False,
        "old_212_terms_and_errors_unchanged": True,
        "older_germ_unproved_status_retained": True,
        "older_germ_used_as_proof": False,
    }
    guards = {
        "self_format": _self_format_guard(root_path),
        "Q4": q4_check,
        "C5": loop_check,
        "R100": outer_check,
        "AST": copy.deepcopy(formal["source_formulas"]),
        "source_execution": False,
        "source_values_selected": False,
    }
    result = {
        "passed": True,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": parent.get("source_family"),
        "source_pins": source_pins,
        "input_pins": copy.deepcopy(source_pins),
        "source_guards": guards,
        "parent_scope": parent_scope,
        "source_contract": source_contract,
        "contract": copy.deepcopy(source_contract),
        "background_packet": background,
        "charts": chart_output,
        "endpoint_identity_receipts": formal,
        "joined_source_entry": list(PARENT_JOINED),
        "joined_pressure_source": join,
        "graph_patch": graph_patch,
        "limits": limits,
        "signed_relation_not_global_closure": True,
        "accepted_ancestors_replayed": False,
        "source_execution": False,
    }
    json.dumps(result, sort_keys=True, ensure_ascii=True, allow_nan=False)
    return result


__all__ = ["build"]
