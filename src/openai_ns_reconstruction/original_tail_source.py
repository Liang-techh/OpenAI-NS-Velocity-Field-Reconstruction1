"""Archive-only source correspondence for the two post-:math:`R_c` seams.

This adapter deliberately stops at source functions.  It does not construct an
accepted owner, call a producer/checker, evaluate a field, or turn a saved
enclosure into a value.  The compact and quiet archives use endpoint packets
for their normalized fields.  Their physical ``U``, ``U_y`` and ``m`` rows are
therefore reported as scoped full-chart recipes when an endpoint wrapper is
not serialized; no synthetic graph node is invented.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT


_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_PINNED = {
    "outer_archive": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz",
        "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d",
    ),
    "outer_source": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py",
        "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    ),
    "outer_receipt": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check.json",
        "ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba",
    ),
    "compact_archive": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz",
        "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
    ),
    "compact_source": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py",
        "f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b",
    ),
    "compact_receipt": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check.json",
        "6c933d077f9c151eda246acd4d5f694d32dc8272fc3fb65a031b9b9cf28c8e7b",
    ),
    "quiet_archive": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz",
        "bbc8ae0361bdfbb75255707e7a234eb39f531533e6f39ad718cf262e27ad42c4",
    ),
    "quiet_source": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.py",
        "0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068",
    ),
    "quiet_receipt": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions_check.json",
        "3cb8b2e7c48ee62ce2b3dba9c6b0d68b901010d6a154b18d466d3221e2c8e7c0",
    ),
}

_ARCHIVE_KEYS = {
    "outer": "true_outer_C6_C5_bindings",
    "compact": "true_compact_C5_bindings",
    "quiet": "true_quiet_Rp_C5_C3_bindings",
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _handle(value: Any) -> int | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            return None
    if isinstance(value, Mapping):
        for key in ("node", "handle", "id"):
            if key in value:
                return _handle(value[key])
    return None


def _node(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    handle = _handle(value)
    if handle is None or handle < 0:
        return None
    if isinstance(nodes, Mapping):
        row = nodes.get(handle, nodes.get(str(handle)))
    else:
        try:
            row = nodes[handle]
        except (IndexError, KeyError, TypeError):
            row = None
    return row if isinstance(row, Mapping) else None


def _children(row: Mapping[str, Any]) -> list[int]:
    children: list[int] = []
    for key in (
        "arguments",
        "numerator",
        "denominator",
        "argument",
        "expression",
        "integrand",
        "lower",
        "upper",
        "bounds",
        "body",
    ):
        value = row.get(key)
        values = value if isinstance(value, list) else [value]
        for item in values:
            handle = _handle(item)
            if handle is not None:
                children.append(handle)
    return children


def _walk(nodes: Any, starts: Any, *, limit: int = 100_000) -> tuple[set[int], list[tuple[int, Mapping[str, Any]]]]:
    stack = [_handle(item) for item in (starts if isinstance(starts, list) else [starts])]
    stack = [item for item in stack if item is not None]
    seen: set[int] = set()
    rows: list[tuple[int, Mapping[str, Any]]] = []
    while stack and len(seen) < limit:
        handle = stack.pop()
        if handle in seen:
            continue
        row = _node(nodes, handle)
        if row is None:
            continue
        seen.add(handle)
        rows.append((handle, row))
        stack.extend(child for child in _children(row) if child not in seen)
    return seen, rows


def _reachable(nodes: Any, start: Any, target: Any) -> bool:
    """Return whether a saved graph expression contains a saved child handle."""
    target_handle = _handle(target)
    return target_handle is not None and target_handle in _walk(nodes, start)[0]


def _rows(value: Any, count: int | None = None) -> list[int]:
    if not isinstance(value, list):
        return []
    if count is not None and len(value) < count:
        return []
    output: list[int] = []
    for item in value if count is None else value[:count]:
        handle = _handle(item)
        if handle is None:
            return []
        output.append(handle)
    return output


def _all_handles(nodes: Any, values: Any) -> bool:
    if not isinstance(values, list):
        return False
    return bool(values) and all(_node(nodes, item) is not None for item in values)


def _nonzero_handles(values: Any) -> bool:
    return isinstance(values, list) and bool(values) and all(_handle(item) not in (None, 0) for item in values)


def _nested_equal(left: Any, right: Any) -> bool:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        return set(left) == set(right) and all(_nested_equal(left[key], right[key]) for key in left)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(_nested_equal(a, b) for a, b in zip(left, right))
    return left == right


def _repo_root(root: Any) -> Path:
    if isinstance(root, Path):
        return root
    if isinstance(root, str):
        return Path(root)
    return Path(getattr(root, "root", ROOT))


def _load_serialized(value: Any, repo: Path, name: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"passed": True, "mode": "supplied_mapping", "name": name}
    if not isinstance(value, (str, Path)):
        return {}, {"passed": False, "name": name, "reason": "mapping_or_serialized_path_required"}
    path = Path(value)
    if not path.is_absolute():
        path = repo / path
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError) as exc:
        return {}, {"passed": False, "name": name, "path": str(path), "reason": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {
        "passed": isinstance(data, Mapping),
        "name": name,
        "path": str(path),
    }


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def _lookup_pin(mapping: Mapping[str, Any] | None, relative: str) -> tuple[bool, Any]:
    if not isinstance(mapping, Mapping):
        return False, None
    for key in (relative, Path(relative).name):
        if key not in mapping:
            continue
        value = mapping[key]
        if isinstance(value, Mapping):
            for digest_key in ("sha256", "hash", "digest"):
                if digest_key in value:
                    return True, value[digest_key]
        return True, value
    return False, None


def _file_pin(repo: Path, name: str, hashes: Mapping[str, Any] | None, code_hashes: Mapping[str, Any] | None) -> dict[str, Any]:
    relative, expected = _PINNED[name]
    path = repo / relative
    try:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        return {"path": relative, "expected_sha256": expected, "actual_sha256": None, "passed": False, "error": str(exc)}
    supplied = code_hashes if relative.startswith(("src/", "scripts/")) else hashes
    present, owner_value = _lookup_pin(supplied, relative)
    owner_match = not present or owner_value == expected
    return {
        "path": relative,
        "expected_sha256": expected,
        "actual_sha256": actual,
        "owner_map_entry_present": present,
        "owner_map_entry_matches": owner_match,
        "passed": actual == expected and owner_match,
    }


def _source_statement(text: str, tokens: list[str]) -> dict[str, Any]:
    rows = [{"statement": token, "present": token in text} for token in tokens]
    return {"passed": bool(rows) and all(row["present"] for row in rows), "rows": rows}


def _function_wrapper(nodes: Any, handle: Any, coordinate: int, query: int) -> dict[str, Any]:
    h = _handle(handle)
    row = _node(nodes, h)
    if row is None:
        return {"node": h, "passed": False, "reason": "missing_node"}
    variable = _handle(row.get("variable"))
    value = _handle(row.get("value"))
    operation = row.get("operation")
    z_free = row.get("Z_independent_substitution") is not False
    return {
        "node": h,
        "operation": operation,
        "variable": variable,
        "value": value,
        "Z_independent_substitution": z_free,
        "passed": bool(operation == "function_substitution" and variable == coordinate and value == query and z_free),
    }


def _scope(chart: str, coordinate: int, query: int, *, query_value: str, radius: Any) -> dict[str, Any]:
    return {
        "chart": chart,
        "coordinate_node": coordinate,
        "native_query_node": query,
        "native_query_value": query_value,
        "radius_node": _handle(radius),
        "ordinary_parameter": "Z",
        "ordinary_Z_orders": [0, 1, 2, 3],
    }


def _flat_rows(nodes: Any, starts: Any) -> dict[str, Any]:
    seen, rows = _walk(nodes, starts)
    beta = [(h, row) for h, row in rows if row.get("operation") == "compact_raw_beta"]
    lazy = [(h, row) for h, row in rows if row.get("operation") == "original_lazy_flat_branch"]
    valid_lazy = [
        h
        for h, row in lazy
        if row.get("flat_value") == 0
        and row.get("active_body_not_evaluated_on_flat_branch") is True
        and row.get("original_beta_exact_flat_derivative") is True
    ]
    valid_beta = [h for h, row in beta if isinstance(row.get("exact_definition"), str) and "exact0 elsewhere" in row["exact_definition"]]
    return {
        "reachable_node_count": len(seen),
        "compact_raw_beta_nodes": sorted(valid_beta)[:24],
        "lazy_flat_branch_nodes": sorted(valid_lazy)[:24],
        "raw_beta_count": len(valid_beta),
        "lazy_flat_count": len(valid_lazy),
        "passed": bool(valid_beta or valid_lazy),
        "direct_endpoint_value_evaluated": False,
    }


def _empty_endpoint_integrals(nodes: Any, starts: Any, coordinate: int, endpoint_query: int) -> dict[str, Any]:
    seen, rows = _walk(nodes, starts)
    integrals = []
    for handle, row in rows:
        if row.get("operation") != "definite_integral":
            continue
        lower = _handle(row.get("lower"))
        upper = _handle(row.get("upper"))
        if lower == 1 and upper == coordinate:
            integrals.append({"node": handle, "lower_node": lower, "upper_coordinate_node": upper, "endpoint_query_node": endpoint_query})
    return {
        "rows": integrals[:64],
        "integral_count": len(integrals),
        "coordinate_node": coordinate,
        "endpoint_query_node": endpoint_query,
        "passed": bool(integrals),
        "meaning": "local lower=upper at compact x=1; incoming histories remain separate source functions",
    }


def _prefix_check(nodes: Any, archive_nodes: Any, handles: list[int]) -> dict[str, Any]:
    try:
        current_count = len(nodes)
    except TypeError:
        current_count = 0
    archive_count = len(archive_nodes) if isinstance(archive_nodes, list) else 0
    rows = []
    for handle in sorted(set(handles)):
        current = _node(nodes, handle)
        archived = _node(archive_nodes, handle)
        # The caller may pass the newest quiet graph while each stage archive
        # is an earlier immutable prefix.  A handle introduced by a later
        # stage is therefore outside this stage's prefix and must be skipped,
        # rather than reported as a false mismatch.  Handles that exist in the
        # stage archive remain strict byte-for-byte checks.
        if archived is None:
            continue
        rows.append({"node": handle, "present": current is not None, "matches_archive": current is not None and archived is not None and current == archived})
    return {
        "current_node_count": current_count,
        "archive_node_count": archive_count,
        "rows_checked": len(rows),
        "passed": bool(rows) and all(row["matches_archive"] for row in rows),
        "rows": rows[:80],
    }


def _cross_map(left: Any, right: Any, keys: list[str]) -> dict[str, Any]:
    rows = []
    for key in keys:
        lvalue = left.get(key) if isinstance(left, Mapping) else None
        rvalue = right.get(key) if isinstance(right, Mapping) else None
        rows.append({"field": key, "left": lvalue, "right": rvalue, "passed": _nested_equal(lvalue, rvalue)})
    return {"rows": rows, "passed": bool(rows) and all(row["passed"] for row in rows)}


def _common_parameters(outer: Mapping[str, Any], compact: Mapping[str, Any], quiet: Mapping[str, Any]) -> dict[str, Any]:
    rows = []
    for field in ("N", "mu"):
        values = [outer.get(field), compact.get(field), quiet.get(field)]
        rows.append({"field": field, "values": values, "passed": values[0] == values[1] == values[2] and values[0] is not None})
    for field in ("Pstar", "delta"):
        values = [compact.get(field), quiet.get(field)]
        rows.append({"field": field, "values": values, "passed": values[0] == values[1] and values[0] is not None})
    return {"rows": rows, "passed": bool(rows) and all(row["passed"] for row in rows)}


def _field_rows(
    *,
    seam: str,
    quantity: str,
    left_nodes: list[int],
    right_nodes: list[int],
    left_scope: Mapping[str, Any],
    right_scope: Mapping[str, Any],
    source_identity: Mapping[str, Any],
    orders: int,
    left_recipe: bool = False,
    right_recipe: bool = False,
    left_backing: list[int] | None = None,
    right_backing: list[int] | None = None,
) -> list[dict[str, Any]]:
    rows = []
    for order in range(orders):
        left = left_nodes[order] if order < len(left_nodes) else None
        right = right_nodes[order] if order < len(right_nodes) else None
        rows.append({
            "seam": seam,
            "quantity": quantity,
            "ordinary_Z_order": order,
            "left_node": left,
            "right_node": right,
            "left_endpoint_backing_node": (left_backing or [None] * orders)[order] if left_backing and order < len(left_backing) else None,
            "right_endpoint_backing_node": (right_backing or [None] * orders)[order] if right_backing and order < len(right_backing) else None,
            "left_scope": dict(left_scope),
            "right_scope": dict(right_scope),
            "left_uses_local_endpoint_recipe": left_recipe,
            "right_uses_local_endpoint_recipe": right_recipe,
            # A local full-chart recipe is deliberately not mislabeled as an
            # archived endpoint node.  The actual backing row, where one is
            # serialized, is carried separately above.
            "left_recipe_is_archived_node": False,
            "right_recipe_is_archived_node": False,
            "common_function_identity": dict(source_identity),
            "passed": bool(left is not None and right is not None and source_identity.get("passed") is True),
        })
    return rows


def _leading_rows(
    *,
    seam: str,
    left: Mapping[str, Any],
    right: Mapping[str, Any],
    left_scope: Mapping[str, Any],
    right_scope: Mapping[str, Any],
    source_identity: Mapping[str, Any],
) -> dict[str, Any]:
    rows = []
    for quantity, count in (("V", 3), ("m", 3), ("V_y", 2)):
        lv = _rows(left.get(quantity), count)
        rv = _rows(right.get(quantity), count)
        for order in range(count):
            rows.append({
                "seam": seam,
                "quantity": quantity,
                "ordinary_Z_order": order,
                "left_node": lv[order] if len(lv) > order else None,
                "right_node": rv[order] if len(rv) > order else None,
                "left_scope": dict(left_scope),
                "right_scope": dict(right_scope),
                "source_function_identity": dict(source_identity),
                "independent_of_modified_Qy": True,
                "passed": bool(len(lv) == count and len(rv) == count and source_identity.get("passed") is True),
            })
    return {"rows": rows, "passed": bool(rows) and all(row["passed"] for row in rows), "field_order": ["V", "m", "V_y"]}


def _outer_checks(outer: Mapping[str, Any], nodes: Any, source_text: str, receipt: Mapping[str, Any]) -> dict[str, Any]:
    binding = _mapping(outer.get("true_outer_C6_C5_bindings"))
    endpoint = _mapping(_mapping(binding.get("endpoints")).get("O3_power"))
    chart = _mapping(_mapping(binding.get("charts")).get("O3_power"))
    background = _mapping(_mapping(binding.get("background_endpoints")).get("O3_power"))
    coordinate = _handle(chart.get("coordinate")); query = 48
    endpoint_rows = []
    for group_name, quantity, count in (("physical", "U", 3), ("physical", "U_y", 2), ("physical", "m", 4), ("physical", "Q_y", 2), ("boundary", "R_Q_y", 2)):
        values = _rows(_mapping(endpoint.get(group_name)).get(quantity), count)
        endpoint_rows.append({"group": group_name, "quantity": quantity, "rows": values, "passed": len(values) == count and all(_function_wrapper(nodes, h, coordinate or -1, query)["passed"] for h in values)})
    flat_u = _flat_rows(nodes, _rows(_mapping(endpoint.get("physical")).get("U"), 3))
    flat_uy = _flat_rows(nodes, _rows(_mapping(endpoint.get("physical")).get("U_y"), 2))
    source = _source_statement(source_text, [
        "V=a.fixed(0);Ey=a.scale(E,g.neg(g.add(g.constant('1/2'),g.mul(self.mu,sigma(g,x)))));Vy=a.fixed(0)",
        "self.Rc=self.endpoints['O3_power'];self.Rc_background=self.background_endpoints['O3_power']",
        "physical_radius_power=power",
        "phase=g.unary('fractional_part'",
    ])
    flags = {
        "outer_c6_c5_gate": receipt.get("current_original_true_six_outer_prescribed_C6_generic_C5_functions_installed") is True,
        "outer_finite_N_gate": receipt.get("current_original_true_six_outer_charts_finite_N_C5_functions_installed") is True,
        "outer_modified_Q4_gate": receipt.get("current_original_true_all17_modified_Q4_total_radial_C5_functions_installed") is True,
        "outer_pressure_C3_gate": receipt.get("current_original_true_rminus_to_Rc_signed_pressure_C3_functions_installed") is True,
        "O3_power_endpoint_coordinate": coordinate == 250215,
        "O3_power_endpoint_query": bool(
            _rows(_mapping(endpoint.get("physical")).get("U"), 1)
            and _function_wrapper(nodes, _rows(_mapping(endpoint.get("physical")).get("U"), 1)[0], coordinate or -1, query)["passed"]
        ),
    }
    leading = {
        "V": _rows(_mapping(background.get("fields")).get("V"), 3),
        "m": _rows(_mapping(background.get("fields")).get("m"), 3),
        "V_y": _rows(_mapping(background.get("radial")).get("V_y"), 2),
    }
    leading_ok = all(len(value) == count and _nonzero_handles(value) for value, count in ((leading["V"], 3), (leading["m"], 3), (leading["V_y"], 2)))
    return {
        "binding_key": "true_outer_C6_C5_bindings",
        "endpoint_scope": {"chart": "O3_power", "coordinate_node": coordinate, "query_node": query, "query_value": "2", "radius_node": _handle(endpoint.get("radius"))},
        "endpoint_rows": endpoint_rows,
        "flat_cutoff": {
            "U": flat_u,
            "U_y": flat_uy,
            "required_here": False,
            "passed": True,
            "meaning": "O3_power is the non-flat upstream endpoint; compact x=1 carries the local exact flat cutoff",
            "source_text": source,
        },
        "leading_inputs": {**leading, "passed": leading_ok, "source_handles_are_nonzero": leading_ok},
        "source_receipt_flags": flags,
        "source_text": source,
        "passed": bool(source["passed"] and all(flags.values()) and all(row["passed"] for row in endpoint_rows) and leading_ok),
    }


def _compact_checks(compact: Mapping[str, Any], nodes: Any, source_text: str, receipt: Mapping[str, Any], outer: Mapping[str, Any], outer_checks: Mapping[str, Any]) -> dict[str, Any]:
    binding = _mapping(compact.get("true_compact_C5_bindings"))
    coordinate = _handle(_mapping(binding.get("background")).get("phase"))
    background = _mapping(binding.get("background")); inlet = _mapping(binding.get("inlet")); terminal = _mapping(binding.get("terminal"))
    source = _source_statement(source_text, [
        "self.A=outer.Rc_background['fields']['E']",
        "self.incoming_background=outer.Rc_background['history']",
        "self.incoming_correction=outer.Rc['correction']",
        "self.incoming_integrals=outer.Rc['integrals']",
        "R=g.mul(self.Rc,x)",
        "self.terminal=self.at(g.constant(2));self.inlet=self.at(g.one)",
        "self.terminal_certificates={key:g.node('proved_original_compact_relative_C5_terminal_identity'",
        "def profiles(g,W,N,A,h,x)",
    ])
    outer_binding = _mapping(outer.get("true_outer_C6_C5_bindings"))
    outer_ep = _mapping(_mapping(outer_binding.get("endpoints")).get("O3_power"))
    outer_bg = _mapping(_mapping(outer_binding.get("background_endpoints")).get("O3_power"))
    cross = {
        "Rc_equals_outer_endpoint_radius": binding.get("Rc") == outer_ep.get("radius"),
        "A_equals_outer_background_E": _nested_equal(binding.get("A"), _mapping(outer_bg.get("fields")).get("E")),
        "incoming_background_equals_outer_background_history": _nested_equal(binding.get("incoming_background"), _mapping(outer_bg.get("history"))),
        "incoming_correction_equals_outer_correction": _nested_equal(binding.get("incoming_correction"), _mapping(outer_ep.get("correction"))),
        "incoming_integrals_equals_outer_integrals": _nested_equal(binding.get("incoming_integrals"), _mapping(outer_ep.get("integrals"))),
        "initial_equals_outer_initial": _nested_equal(binding.get("initial"), outer_binding.get("initial")),
    }
    common = _common_parameters(outer_binding, binding, binding)
    endpoint_values = {"x1": 1, "x2": 48}
    scope = {"chart": "true_compact_repair", "coordinate_node": coordinate, "radius_node": _handle(background.get("radius")), "domain": [_handle(v) for v in background.get("domain", [])]}
    inlet_m = _rows(_mapping(inlet.get("modified_fields")).get("m"), 4)
    terminal_m = _rows(_mapping(terminal.get("modified_fields")).get("m"), 4)
    local_empty = _empty_endpoint_integrals(nodes, inlet_m, coordinate or -1, 1)
    flat_start = _rows(_mapping(binding.get("physical")).get("U"), 3) + _rows(_mapping(binding.get("physical")).get("U_y"), 2)
    flat = _flat_rows(nodes, flat_start)
    endpoint_rows = {
        "inlet_modified_m": inlet_m,
        "terminal_modified_m": terminal_m,
        "inlet_background_m": _rows(_mapping(inlet.get("background_fields")).get("m"), 3),
        "terminal_background_m": _rows(_mapping(terminal.get("background_fields")).get("m"), 3),
        "inlet_background_V": _rows(_mapping(inlet.get("background_fields")).get("V"), 3),
        "terminal_background_V": _rows(_mapping(terminal.get("background_fields")).get("V"), 3),
        "inlet_background_V_y": _rows(_mapping(inlet.get("radial")).get("V_y"), 2),
        "terminal_background_V_y": _rows(_mapping(terminal.get("radial")).get("V_y"), 2),
    }
    endpoint_rows_passed = all(_nonzero_handles(value) for value in endpoint_rows.values())
    flags = {
        "compact_targets_gate": receipt.get("current_original_true_all17_normalized_compact_targets_C5_functions_installed") is True,
        "compact_branch_gate": receipt.get("current_original_true_same_branch_compact_controls_C5_functions_installed") is True,
        "compact_modified_C5_gate": receipt.get("current_original_true_Rc_to_2Rc_compact_modified_history_C5_functions_installed") is True,
        "compact_relative_terminal_gate": receipt.get("current_original_true_Rc_to_2Rc_relative_five_C5_zero_identities_installed") is True,
        "compact_pressure_gate": receipt.get("current_original_true_rminus_to_2Rc_signed_pressure_C3_functions_installed") is True,
    }
    return {
        "binding_key": "true_compact_C5_bindings",
        "endpoint_scope": scope,
        "queries": endpoint_values,
        "cross_archive_source_wiring": {"rows": [{"field": key, "passed": value} for key, value in cross.items()], "passed": all(cross.values())},
        "common_parameters": common,
        "endpoint_rows": endpoint_rows,
        "local_flat_cutoff": {"inlet": local_empty, "physical_ports": flat, "source_text": source},
        "source_receipt_flags": flags,
        "source_text": source,
        "passed": bool(source["passed"] and all(flags.values()) and all(cross.values()) and common["passed"] and endpoint_rows_passed and local_empty["passed"] and flat["passed"]),
    }


def _terminal_mass_certificate(compact: Mapping[str, Any], nodes: Any, source_text: str, receipt: Mapping[str, Any]) -> dict[str, Any]:
    binding = _mapping(compact.get("true_compact_C5_bindings"))
    certs = _mapping(binding.get("terminal_certificates")); cert_handle = _handle(certs.get("m")); cert = _node(nodes, cert_handle)
    terminal = _mapping(binding.get("terminal")); terminal_corr = _rows(_mapping(terminal.get("correction")).get("m"), 6)
    rhs = _rows(_mapping(binding.get("terminal_rhs")).get("m"), 6)
    target = _rows(_mapping(binding.get("targets")).get("M"), 6)
    incoming = _rows(_mapping(binding.get("N_scaled_incoming")).get("m"), 6)
    residual = _rows(_mapping(binding.get("residual")).get("M"), 6)
    implicit = _mapping(binding.get("implicit")); recurrence = implicit.get("recurrence")
    W = _mapping(binding.get("W")); B = W.get("B")
    h = binding.get("h")
    m_h = h[0] if isinstance(h, list) and h and isinstance(h[0], list) else []
    source = _source_statement(source_text, [
        "self.terminal_rhs={}",
        "for key,row,amplitude in (('m','M',a.trim(self.A))",
        "self.terminal_certificates={key:g.node('proved_original_compact_relative_C5_terminal_identity'",
        "zero_rows=[g.zero.node]*6,relative_correction_only=True,absolute_leading_histories_and_pressure_not_zeroed=True",
        "original_signed_endpoint_binding=self.bindings['original_partial_history_terminal_identity']",
    ])
    recurrence_rows = recurrence if isinstance(recurrence, list) else []
    recurrence_ok = bool(recurrence_rows) and [row.get("order") for row in recurrence_rows] == [1, 2, 3, 4, 5] and all(row.get("same_full_Jacobian") is True for row in recurrence_rows)
    A_rows = _rows(binding.get("A"), 1)
    A0 = A_rows[0] if A_rows else None
    target_rows = []
    for order, target_handle in enumerate(target):
        target_node = _node(nodes, target_handle)
        numerator = _handle(target_node.get("numerator")) if target_node else None
        denominator = _handle(target_node.get("denominator")) if target_node else None
        target_rows.append({
            "ordinary_Z_order": order,
            "target_node": target_handle,
            "operation": target_node.get("operation") if target_node else None,
            "numerator_node": numerator,
            "denominator_node": denominator,
            "incoming_row_node": incoming[order] if order < len(incoming) else None,
            "previous_target_nodes": target[:order],
            "source_positive_certificate": target_node.get("source_positive_certificate") if target_node else None,
            "numerator_contains_incoming_row": bool(target_node and numerator is not None and _reachable(nodes, numerator, incoming[order])),
            "numerator_contains_previous_targets": bool(target_node and numerator is not None and all(_reachable(nodes, numerator, previous) for previous in target[:order])),
            "passed": bool(
                target_node
                and target_node.get("operation") == "positive_quotient"
                and denominator == A0
                and numerator is not None
                and _node(nodes, numerator) is not None
                and target_node.get("source_positive_certificate") == "original positive Rc amplitude"
                and _reachable(nodes, numerator, incoming[order])
                and all(_reachable(nodes, numerator, previous) for previous in target[:order])
            ),
        })
    residual_rows = []
    h_rows = h if isinstance(h, list) else []
    for order, residual_handle in enumerate(residual):
        residual_node = _node(nodes, residual_handle)
        arguments = _rows(residual_node.get("arguments")) if residual_node else []
        correction_row = h_rows[0][order] if len(h_rows) > 0 and order < len(h_rows[0]) else None
        leading_row = h_rows[1][order] if len(h_rows) > 1 and order < len(h_rows[1]) else None
        residual_rows.append({
            "ordinary_Z_order": order,
            "residual_node": residual_handle,
            "arguments": arguments,
            "target_node": target[order] if order < len(target) else None,
            "correction_control_node": correction_row,
            "leading_control_node": leading_row,
            "passed": bool(
                residual_node
                and residual_node.get("operation") == "sum"
                and arguments
                and arguments[0] == target[order]
                and len(arguments) == 2
                and correction_row is not None
                and leading_row is not None
                and _reachable(nodes, arguments[1], correction_row)
                and _reachable(nodes, arguments[1], leading_row)
            ),
        })
    rhs_rows = []
    common_rhs_scale = None
    for order, rhs_handle in enumerate(rhs):
        rhs_node = _node(nodes, rhs_handle)
        arguments = _rows(rhs_node.get("arguments")) if rhs_node else []
        polynomial = arguments[0] if arguments else None
        scale = arguments[1] if len(arguments) > 1 else None
        if common_rhs_scale is None:
            common_rhs_scale = scale
        rhs_rows.append({
            "ordinary_Z_order": order,
            "terminal_rhs_node": rhs_handle,
            "polynomial_node": polynomial,
            "common_endpoint_scale_node": scale,
            "residual_node": residual[order] if order < len(residual) else None,
            "polynomial_contains_residual": bool(polynomial is not None and order < len(residual) and _reachable(nodes, polynomial, residual[order])),
            "passed": bool(
                rhs_node
                and rhs_node.get("operation") == "product"
                and len(arguments) == 2
                and polynomial is not None
                and scale is not None
                and common_rhs_scale == scale
                and order < len(residual)
                and _reachable(nodes, polynomial, residual[order])
            ),
        })
    m_row_identity = {
        "method": "saved_target_positive_quotient_then_full_implicit_residual_and_terminal_rhs_graph",
        "A0_denominator_node": A0,
        "target_rows": target_rows,
        "residual_rows": residual_rows,
        "terminal_rhs_rows": rhs_rows,
        "common_terminal_rhs_scale_node": common_rhs_scale,
        "passed": bool(
            A0 is not None
            and len(target_rows) == len(residual_rows) == len(rhs_rows) == 6
            and all(row["passed"] for row in target_rows + residual_rows + rhs_rows)
        ),
        "numeric_value_selected": False,
    }
    cert_ok = bool(
        cert is not None
        and cert.get("operation") == "proved_original_compact_relative_C5_terminal_identity"
        and cert.get("source_system") == _handle(binding.get("system"))
        and cert.get("actual_endpoint_rows") == terminal_corr
        and cert.get("equivalent_residual_rows") == rhs
        and cert.get("zero_rows") == [0] * 6
        and cert.get("relative_correction_only") is True
        and cert.get("absolute_leading_histories_and_pressure_not_zeroed") is True
        and recurrence_ok
    )
    transparent = {
        "target_M_rows": target,
        "N_scaled_incoming_m_rows": incoming,
        "residual_M_rows": residual,
        "implicit_system_node": _handle(binding.get("system")),
        "implicit_recurrence_orders": [row.get("order") for row in recurrence_rows],
        "m_control_row": m_h,
        "m_control_row_is_source_graph": _all_handles(nodes, m_h),
        "M_linear_B_row": B[0] if isinstance(B, list) and B else None,
        "M_linear_B_row_is_[1,1,0,0,0]": isinstance(B, list) and bool(B) and B[0] == [1, 1, 0, 0, 0],
        "terminal_rhs_matches_certificate": rhs == (cert.get("equivalent_residual_rows") if cert else None),
        "M_row_identity": m_row_identity,
    }
    transparent["passed"] = bool(
        cert_ok
        and source["passed"]
        and len(target) == len(incoming) == len(residual) == len(terminal_corr) == len(rhs) == 6
        and transparent["m_control_row_is_source_graph"]
        and transparent["M_linear_B_row_is_[1,1,0,0,0]"]
        and transparent["terminal_rhs_matches_certificate"]
        and m_row_identity["passed"]
    )
    return {
        "certificate_node": cert_handle,
        "certificate": {key: cert.get(key) for key in ("operation", "actual_endpoint_rows", "equivalent_residual_rows", "source_system", "zero_rows", "relative_correction_only", "absolute_leading_histories_and_pressure_not_zeroed", "theorem")} if cert else {},
        "transparent_m_identity": transparent,
        "source_text": source,
        "source_receipt_flags": {
            "compact_actual_source": _mapping(_mapping(receipt.get("checks")).get("actual_source")).get("passed") is True,
            "terminal_integral_identity": _mapping(_mapping(receipt.get("checks")).get("actual_inlet_and_endpoint")).get("terminal_integral_identity_uses_pinned_original_exact_bump_theorem") is True,
            "relative_terminal_gate": receipt.get("current_original_true_Rc_to_2Rc_relative_five_C5_zero_identities_installed") is True,
        },
        "ordinary_Z_orders": list(range(6)),
        "relative_correction_only": True,
        "absolute_leading_histories_not_zeroed": bool(cert and cert.get("absolute_leading_histories_and_pressure_not_zeroed") is True),
        "m_relative_zero_source_identity": transparent["passed"],
        "absolute_modified_m_equals_leading_m": transparent["passed"],
        "numeric_value_selected": False,
        "passed": transparent["passed"],
    }


def _quiet_checks(quiet: Mapping[str, Any], nodes: Any, source_text: str, receipt: Mapping[str, Any], compact: Mapping[str, Any], compact_checks: Mapping[str, Any], terminal_mass: Mapping[str, Any]) -> dict[str, Any]:
    binding = _mapping(quiet.get("true_quiet_Rp_C5_C3_bindings")); inlet = _mapping(binding.get("inlet")); incoming = _mapping(binding.get("incoming")); leading_incoming = _mapping(binding.get("leading_incoming")); relative_incoming = _mapping(binding.get("relative_incoming"))
    cb = _mapping(compact.get("true_compact_C5_bindings")); terminal = _mapping(cb.get("terminal"))
    source = _source_statement(source_text, [
        "self.incoming=compact.terminal",
        "self.leading_incoming=compact.terminal['background_fields']",
        "self.relative_incoming=compact.terminal['correction']",
        "self.incoming_integrals=compact.terminal['integrals']",
        "self.incoming_pressure=compact.terminal['pressure']",
        "self.flow=power_flow(g,self.incoming['modified_fields']['E']",
        "self.leading_flow=power_flow(g,self.leading_incoming['E']",
        "self.linear_mass_theorem=g.node('proved_original_quiet_linear_mass_source_identity'",
        "compact_relative_zero=compact.terminal_certificates['m'].node",
        "original_linear_Z_mass_and_nonzero_scalar_K_source_retained=True",
    ])
    cross_rows = []
    for key in ("correction", "modified_fields", "background_fields", "radial", "generic", "integrals", "pressure", "boundary", "radius"):
        lvalue = incoming.get(key); rvalue = terminal.get(key)
        cross_rows.append({"field": key, "passed": _nested_equal(lvalue, rvalue), "left": lvalue, "right": rvalue})
    cross_rows.extend([
        {"field": "leading_incoming", "passed": _nested_equal(leading_incoming, terminal.get("background_fields")), "left": leading_incoming, "right": terminal.get("background_fields")},
        {"field": "relative_incoming", "passed": _nested_equal(relative_incoming, terminal.get("correction")), "left": relative_incoming, "right": terminal.get("correction")},
        {"field": "R2", "passed": binding.get("R2") == terminal.get("radius"), "left": binding.get("R2"), "right": terminal.get("radius")},
        {"field": "initial", "passed": _nested_equal(binding.get("initial"), cb.get("initial")), "left": binding.get("initial"), "right": cb.get("initial")},
    ])
    common = _common_parameters(cb, cb, binding)
    inlet_phys = _mapping(inlet.get("physical")); inlet_fields = _mapping(inlet.get("fields")); inlet_radial = _mapping(inlet.get("radial"))
    field_rows = {
        "U": _rows(inlet_phys.get("U"), 3),
        "U_y": _rows(inlet_phys.get("U_y"), 2),
        "m": _rows(inlet_phys.get("m"), 4),
    }
    quiet_zero_U = bool(field_rows["U"] and field_rows["U_y"] and all(h == 265439 for h in field_rows["U"] + field_rows["U_y"]))
    terminal_source_flat = _flat_rows(nodes, _rows(_mapping(terminal.get("modified_fields")).get("V"), 3) + _rows(_mapping(terminal.get("radial")).get("V_y"), 2))
    source_flags = {
        "quiet_C5_pressure_gate": receipt.get("current_original_true_2Rc_to_Rp_quiet_C5_pressure_C3_functions_installed") is True,
        "quiet_geometry_gate": receipt.get("current_original_true_quiet_K_matches_original_native_source_installed") is True,
        "quiet_y4_gate": receipt.get("current_original_true_quiet_y4_C5_and_direct_Qyy_C4_functions_installed") is True,
        "quiet_terminal_relative_gate": _mapping(_mapping(receipt.get("checks")).get("actual_inlet")).get("compact_flat_velocity_seam_and_nonzero_pressure_history_not_reset") is True,
        "quiet_source_mass_gate": _mapping(_mapping(receipt.get("checks")).get("native_mass_identity")).get("actual_Rh_source_identity_and_same_original_native_mass_FunctionRef_used") is True,
    }
    return {
        "binding_key": "true_quiet_Rp_C5_C3_bindings",
        "endpoint_scope": {"chart": "quiet_2Rc", "coordinate_node": _handle(binding.get("s")), "query_node": 0, "query_value": "0", "radius_node": _handle(binding.get("R2"))},
        "cross_archive_source_wiring": {"rows": cross_rows, "passed": bool(cross_rows) and all(row["passed"] for row in cross_rows)},
        "common_parameters": common,
        "endpoint_rows": field_rows,
        "quiet_zero_velocity": {"physical_U_rows": field_rows["U"], "physical_Uy_rows": field_rows["U_y"], "all_exact_zero_handle": quiet_zero_U, "passed": quiet_zero_U},
        "terminal_flat_cutoff": terminal_source_flat,
        "source_receipt_flags": source_flags,
        "source_text": source,
        "terminal_mass_certificate_consumed": terminal_mass.get("passed") is True,
        "passed": bool(source["passed"] and all(source_flags.values()) and common["passed"] and all(row["passed"] for row in cross_rows) and quiet_zero_U and terminal_source_flat["passed"] and terminal_mass.get("passed") is True),
    }


def build_tail_source(
    nodes: Any,
    outer_archive: Mapping[str, Any] | str | Path,
    compact_archive: Mapping[str, Any] | str | Path,
    quiet_archive: Mapping[str, Any] | str | Path,
    *,
    root: Any = ROOT,
    hashes: Mapping[str, Any] | None = None,
    code_hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Bind the two post-``R_c`` source seams without replaying an owner."""
    repo = _repo_root(root)
    outer, outer_load = _load_serialized(outer_archive, repo, "outer_archive")
    compact, compact_load = _load_serialized(compact_archive, repo, "compact_archive")
    quiet, quiet_load = _load_serialized(quiet_archive, repo, "quiet_archive")
    archives = {"outer": outer, "compact": compact, "quiet": quiet}
    archive_nodes = {name: archive.get("exact_graph_nodes", []) for name, archive in archives.items()}
    pins = {name: _file_pin(repo, name, hashes, code_hashes) for name in _PINNED}
    source_texts = {}
    source_loads = {}
    receipts = {}
    for stage in ("outer", "compact", "quiet"):
        source_name = stage + "_source"; receipt_name = stage + "_receipt"
        source_path = repo / _PINNED[source_name][0]; receipt_path = repo / _PINNED[receipt_name][0]
        try:
            source_texts[stage] = source_path.read_text(encoding="utf-8")
            source_loads[stage] = {"passed": True, "path": _PINNED[source_name][0]}
        except OSError as exc:
            source_texts[stage] = ""; source_loads[stage] = {"passed": False, "error": str(exc)}
        try:
            receipts[stage] = json.loads(receipt_path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            receipts[stage] = {}
    family_rows = []
    for stage, archive in archives.items():
        family_rows.append({"stage": stage, "family": archive.get("source_family"), "matches_expected": archive.get("source_family") == _FAMILY})
    family_passed = bool(family_rows) and all(row["matches_expected"] for row in family_rows)
    selected = [
        250215,
        48,
        258799,
        258805,
        258811,
        251016,
        251017,
        260181,
        260269,
        260275,
        260281,
        262269,
        262275,
        262281,
        262720,
        262769,
        262770,
        262819,
        263344,
        265439,
        265616,
    ]
    prefix_rows = {stage: _prefix_check(nodes, archive_nodes[stage], selected) for stage in archives}
    outer_checks = _outer_checks(outer, archive_nodes["outer"], source_texts["outer"], receipts["outer"])
    compact_checks = _compact_checks(compact, archive_nodes["compact"], source_texts["compact"], receipts["compact"], outer, outer_checks)
    terminal_mass = _terminal_mass_certificate(compact, archive_nodes["compact"], source_texts["compact"], receipts["compact"])
    quiet_checks = _quiet_checks(quiet, archive_nodes["quiet"], source_texts["quiet"], receipts["quiet"], compact, compact_checks, terminal_mass)
    outer_binding = _mapping(outer.get("true_outer_C6_C5_bindings")); compact_binding = _mapping(compact.get("true_compact_C5_bindings")); quiet_binding = _mapping(quiet.get("true_quiet_Rp_C5_C3_bindings"))
    outer_ep = _mapping(_mapping(outer_binding.get("endpoints")).get("O3_power")); outer_bg = _mapping(_mapping(outer_binding.get("background_endpoints")).get("O3_power"))
    compact_bg = _mapping(compact_binding.get("background")); compact_inlet = _mapping(compact_binding.get("inlet")); compact_terminal = _mapping(compact_binding.get("terminal")); quiet_inlet = _mapping(quiet_binding.get("inlet"))
    seam1_identity = {
        "method": "exact_cross_archive_source_transport_and_compact_x1_empty_integrals",
        "passed": bool(compact_checks.get("cross_archive_source_wiring", {}).get("passed") and compact_checks.get("local_flat_cutoff", {}).get("inlet", {}).get("passed") and outer_checks.get("passed")),
        "direct_cross_archive_graph_fingerprint": False,
        "source_text_transport": compact_checks.get("source_text"),
        "incoming_nonzero_correction_preserved": _nonzero_handles(_mapping(compact_binding.get("incoming_correction")).get("m")),
    }
    seam2_identity = {
        "method": "exact_compact_terminal_relative_mass_zero_and_quiet_linear_mass_source",
        "passed": bool(terminal_mass.get("passed") and quiet_checks.get("cross_archive_source_wiring", {}).get("passed")),
        "direct_cross_archive_graph_fingerprint": False,
        "source_text_transport": quiet_checks.get("source_text"),
        "incoming_nonzero_correction_preserved": _nonzero_handles(_mapping(quiet_binding.get("incoming")).get("correction", {}).get("m")),
    }
    outer_scope = _scope("O3_power", 250215, 48, query_value="2", radius=outer_ep.get("radius"))
    compact_x1_scope = _scope("compact_Rc_inlet", 260181, 1, query_value="1", radius=compact_inlet.get("radius"))
    compact_x2_scope = _scope("compact_2Rc_terminal", 260181, 48, query_value="2", radius=compact_terminal.get("radius"))
    quiet_scope = _scope("quiet_2Rc_inlet", 263344, 0, query_value="0", radius=quiet_inlet.get("radius"))
    outer_physical = _mapping(outer_ep.get("physical")); compact_physical = _mapping(compact_binding.get("physical")); quiet_physical = _mapping(_mapping(quiet_inlet).get("physical"))
    outer_leading = {"V": _rows(_mapping(outer_bg.get("fields")).get("V"), 3), "m": _rows(_mapping(outer_bg.get("fields")).get("m"), 3), "V_y": _rows(_mapping(outer_bg.get("radial")).get("V_y"), 2)}
    # Compact's leading flow is the chart background, not either modified
    # endpoint packet.  In particular its radial V_y rows are the exact zero
    # source rows in ``background.radial``.
    compact_background = _mapping(compact_binding.get("background"))
    compact_x1_leading = {"V": _rows(_mapping(compact_background.get("fields")).get("V"), 3), "m": _rows(_mapping(compact_background.get("fields")).get("m"), 3), "V_y": _rows(_mapping(compact_background.get("radial")).get("V_y"), 2)}
    compact_x2_leading = {"V": _rows(_mapping(compact_background.get("fields")).get("V"), 3), "m": _rows(_mapping(compact_background.get("fields")).get("m"), 3), "V_y": _rows(_mapping(compact_background.get("radial")).get("V_y"), 2)}
    # Quiet's original leading flow is serialized separately from the
    # modified inlet fields.  Its first two y-profiles provide V and V_y;
    # leading_fields supplies the original mass history.
    quiet_leading_flow = _mapping(quiet_binding.get("leading_flow"))
    quiet_profiles = quiet_leading_flow.get("profiles_y0_to_y4")
    quiet_profile0 = _mapping(quiet_profiles[0]) if isinstance(quiet_profiles, list) and len(quiet_profiles) >= 1 else {}
    quiet_profile1 = _mapping(quiet_profiles[1]) if isinstance(quiet_profiles, list) and len(quiet_profiles) >= 2 else {}
    quiet_leading = {"V": _rows(quiet_profile0.get("V"), 3), "m": _rows(_mapping(quiet_inlet.get("leading_fields")).get("m"), 3), "V_y": _rows(quiet_profile1.get("V"), 2)}
    leading_identity1 = {"method": "incoming_background_endpoint_transport", "passed": bool(compact_checks.get("cross_archive_source_wiring", {}).get("passed") and compact_checks.get("local_flat_cutoff", {}).get("inlet", {}).get("passed")), "direct_cross_archive_graph_fingerprint": False}
    leading_identity2 = {"method": "quiet_leading_incoming_exact_compact_terminal_background", "passed": bool(quiet_checks.get("cross_archive_source_wiring", {}).get("passed") and quiet_checks.get("source_receipt_flags", {}).get("quiet_source_mass_gate")), "direct_cross_archive_graph_fingerprint": False}
    seams = {}
    seams["power_compact"] = {
        "left_scope": outer_scope,
        "right_scope": compact_x1_scope,
        "modified_U_rows": _field_rows(seam="power_compact", quantity="U", left_nodes=_rows(outer_physical.get("U"), 3), right_nodes=_rows(compact_physical.get("U"), 3), left_scope=outer_scope, right_scope=compact_x1_scope, source_identity=seam1_identity, orders=3, right_recipe=True),
        "modified_Uy_rows": _field_rows(seam="power_compact", quantity="U_y", left_nodes=_rows(outer_physical.get("U_y"), 2), right_nodes=_rows(compact_physical.get("U_y"), 2), left_scope=outer_scope, right_scope=compact_x1_scope, source_identity=seam1_identity, orders=2, right_recipe=True),
        "modified_m_rows": _field_rows(seam="power_compact", quantity="m", left_nodes=_rows(outer_physical.get("m"), 4), right_nodes=_rows(compact_physical.get("m"), 4), left_scope=outer_scope, right_scope=compact_x1_scope, source_identity=seam1_identity, orders=4, right_recipe=True, right_backing=_rows(_mapping(compact_inlet.get("modified_fields")).get("m"), 4)),
        "leading_Qy_inputs": _leading_rows(seam="power_compact", left=outer_leading, right=compact_x1_leading, left_scope=outer_scope, right_scope=compact_x1_scope, source_identity=leading_identity1),
        "actual_endpoint_nodes": {"left_Qy": _rows(outer_physical.get("Q_y"), 2), "left_RQy": _rows(_mapping(outer_ep.get("boundary")).get("R_Q_y"), 2), "right_RQy": _rows(_mapping(compact_inlet.get("boundary")).get("R_Q_y"), 2)},
        "pressure_Qy_RQy_identity_proved_here": False,
    }
    seams["compact_quiet"] = {
        "left_scope": compact_x2_scope,
        "right_scope": quiet_scope,
        "modified_U_rows": _field_rows(seam="compact_quiet", quantity="U", left_nodes=_rows(compact_physical.get("U"), 3), right_nodes=_rows(quiet_physical.get("U"), 3), left_scope=compact_x2_scope, right_scope=quiet_scope, source_identity=seam2_identity, orders=3, left_recipe=True),
        "modified_Uy_rows": _field_rows(seam="compact_quiet", quantity="U_y", left_nodes=_rows(compact_physical.get("U_y"), 2), right_nodes=_rows(quiet_physical.get("U_y"), 2), left_scope=compact_x2_scope, right_scope=quiet_scope, source_identity=seam2_identity, orders=2, left_recipe=True),
        "modified_m_rows": _field_rows(seam="compact_quiet", quantity="m", left_nodes=_rows(compact_physical.get("m"), 4), right_nodes=_rows(quiet_physical.get("m"), 4), left_scope=compact_x2_scope, right_scope=quiet_scope, source_identity=seam2_identity, orders=4, left_recipe=True, left_backing=_rows(_mapping(compact_terminal.get("modified_fields")).get("m"), 4)),
        "leading_Qy_inputs": _leading_rows(seam="compact_quiet", left=compact_x2_leading, right=quiet_leading, left_scope=compact_x2_scope, right_scope=quiet_scope, source_identity=leading_identity2),
        "actual_endpoint_nodes": {"left_RQy": _rows(_mapping(compact_terminal.get("boundary")).get("R_Q_y"), 2), "right_Qy": _rows(quiet_physical.get("Q_y"), 2), "right_RQy": _rows(_mapping(quiet_inlet.get("boundary")).get("R_Q_y"), 2)},
        "pressure_Qy_RQy_identity_proved_here": False,
    }
    for seam in seams.values():
        all_rows = seam["modified_U_rows"] + seam["modified_Uy_rows"] + seam["modified_m_rows"] + seam["leading_Qy_inputs"]["rows"]
        seam["passed"] = bool(all_rows) and all(row["passed"] for row in all_rows) and seam["leading_Qy_inputs"]["passed"]
    source_pins_passed = all(row["passed"] for row in pins.values())
    passed = bool(source_pins_passed and family_passed and all(row["passed"] for row in prefix_rows.values()) and outer_checks["passed"] and compact_checks["passed"] and terminal_mass["passed"] and quiet_checks["passed"] and all(seam["passed"] for seam in seams.values()))
    return {
        "passed": passed,
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "producer_or_checker_executed": False,
        "source_family": {"expected": _FAMILY, "rows": family_rows, "passed": family_passed},
        "source_pins": {"passed": source_pins_passed, "rows": pins},
        "archive_load": {"outer": outer_load, "compact": compact_load, "quiet": quiet_load},
        "source_load": source_loads,
        "source_prefixes": prefix_rows,
        "source_receipts": {"outer": {"all_passed": receipts["outer"].get("all_passed") is True, "explicit_flags_checked": True}, "compact": {"all_passed": receipts["compact"].get("all_passed") is True, "explicit_flags_checked": True}, "quiet": {"all_passed": receipts["quiet"].get("all_passed") is True, "explicit_flags_checked": True}},
        "outer": outer_checks,
        "compact": compact_checks,
        "compact_terminal_mass_zero": terminal_mass,
        "quiet": quiet_checks,
        "seams": seams,
        "terminal_D_zero_inputs": {
            "passed": bool(terminal_mass["passed"] and quiet_checks.get("quiet_zero_velocity", {}).get("passed")),
            "relative_m_exact_source_zero": terminal_mass["m_relative_zero_source_identity"],
            "absolute_modified_m_equals_leading_m": terminal_mass["absolute_modified_m_equals_leading_m"],
            "modified_U_equals_leading_U_at_quiet_inlet": quiet_checks.get("quiet_zero_velocity", {}).get("passed") is True,
            "modified_Uy_equals_leading_Uy_at_quiet_inlet": quiet_checks.get("quiet_zero_velocity", {}).get("passed") is True,
            "pressure_Qy_RQy_recovery_left_to_root": True,
            "numeric_value_selected": False,
        },
        "contract": {
            "modified_fields_bound": ["U", "U_y", "m"],
            "modified_m_ordinary_Z3_included": True,
            "leading_Qy_inputs_bound": ["V", "m", "V_y"],
            "pressure_Qy_RQy_identity_proved_here": False,
            "physical_radius_factor_proved_here": False,
            "ordinary_Z_derivatives_source_C5_or_C6": True,
            "nonzero_histories_retained": True,
            "no_cap_as_function_value": True,
            "no_interval_overlap_alias": True,
        },
        "limits": {
            "pressure_Qy_RQy_recovery": False,
            "signed_numeric_sector_contraction": False,
            "cap_values_as_functions": False,
            "ancestor_constructor_or_checker": False,
        },
        "unproved": {
            "Qy_and_RQy_cross_archive_function_identity": True,
            "pressure_recovery_and_physical_radius_product": True,
            "compact_raw210_numeric_envelope_binding": True,
        },
    }


__all__ = ["build_tail_source"]
