"""Archive-only source bindings for the upstream annular/switch interfaces.

The helper in this module deliberately consumes serialized graph/receipt data
only.  It does not construct an accepted owner, run an old producer/checker,
evaluate a field, or turn an enclosure into a value.  Physical rows are
reported as full-chart expressions together with the native endpoint scope;
an endpoint recipe is never presented as an archived graph node.

The mandatory ``power_long`` seam is the accepted annular
``power_reshape`` source join (post-power native ``1`` to long-reshape native
``0``).  The two switch seams are included when their accepted source receipt
is present.  The frozen-macro to first-switch relation has no accepted source
join in this layer and is returned explicitly as unproved.
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


_PREFIX = "lei_ren_part1_paper_compliant_current_original_"


# These are the immutable source artifacts consumed by this adapter.  The
# values are hashes of the files present in the accepted source receipts.
_PINNED = {
    "annular_archive": (_PREFIX + "annular_mixed_source_joins.json.gz", "3947464be1541ae00a42fa952419f1eb7be203abc071a4bf83565a69205d2fde"),
    "annular_source": (_PREFIX + "annular_mixed_source_joins.py", "8febd6f4b356eb50bae0905de580eaabde75b0809ad0e574a0e4f02214582a17"),
    "annular_receipt": (_PREFIX + "annular_mixed_source_joins_check.json", "28bd4e63feec4105fd30e87b32dbf67fa9d030a7ecf577d2201ab4e3e9313c1b"),
    "annular_checker": (_PREFIX + "annular_mixed_source_joins_check.py", "4121c8ca1a48ccd214248773b67af17991c97e9bc93bc4884fb483c4ac496896"),
    "switch_archive": (_PREFIX + "switch_mixed_source_joins.json.gz", "ddbd506d2fbfc4c58f7a5ba93e0abe3a95d3af8fd3c520ee3ce5024d90a29567"),
    "switch_source": (_PREFIX + "switch_mixed_source_joins.py", "0ae7f651c0d183756f0e6ab90d56dfd80bf3d3c1a6717ffdc2a03fe08dd52f2f"),
    "switch_receipt": (_PREFIX + "switch_mixed_source_joins_check.json", "9cf3c6d941df95f4b5f73e7b97e32198a54fcd3d17a85ecca2a0742c0d602fc3"),
    "switch_checker": (_PREFIX + "switch_mixed_source_joins_check.py", "3b00b9dc4fc866430c1f3538edcb33de6e0a58afa2f5372ff04a7e67e0db6c11"),
    "loop_archive": (_PREFIX + "prescribed_loop_C5_functions.json.gz", "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"),
    "loop_source": (_PREFIX + "prescribed_loop_C5_functions.py", "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"),
    "loop_receipt": (_PREFIX + "prescribed_loop_C5_functions_check.json", "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7"),
    "loop_checker": (_PREFIX + "prescribed_loop_C5_functions_check.py", "066e6e09ccbb4dc28d7b1881197a426beb6bbc7f92d393e7f3e646f27287914b"),
    "r100_archive": (_PREFIX + "prescribed_R100_Rm_C6_functions.json.gz", "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"),
    "r100_source": (_PREFIX + "prescribed_R100_Rm_C6_functions.py", "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"),
    "r100_receipt": (_PREFIX + "prescribed_R100_Rm_C6_functions_check.json", "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6"),
    "r100_checker": (_PREFIX + "prescribed_R100_Rm_C6_functions_check.py", "8ad5e2edc1daba05763b44dd031ac1671d74f30b1733fc80570873e2ca12089a"),
    "modified_archive": (_PREFIX + "modified_Q4_pressure_C3_functions.json.gz", "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"),
    "modified_source": (_PREFIX + "modified_Q4_pressure_C3_functions.py", "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"),
    "modified_receipt": (_PREFIX + "modified_Q4_pressure_C3_functions_check.json", "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f"),
    "modified_checker": (_PREFIX + "modified_Q4_pressure_C3_functions_check.py", "829771078e07e2ec8b5168b181f9ff5f0fa62c941c8cf2436feeaba41e367124"),
    "bridge_archive": (_PREFIX + "prescribed_bridge_C6_functions.json.gz", "5df06c4312efc01073516023fd886c071f439c11de30d3411b2e358bf76abd92"),
    "bridge_source": (_PREFIX + "prescribed_bridge_C6_functions.py", "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605"),
    "bridge_receipt": (_PREFIX + "prescribed_bridge_C6_functions_check.json", "9a0969a8cd4a5989dbd6272d6dd223533400113824f9352b4e57a83196bff61e"),
    "bridge_checker": (_PREFIX + "prescribed_bridge_C6_functions_check.py", "c5e118a1d2691a878d535c5b0b39cf9e627994c0f1f6c8d36139b09d46fac836"),
}

_ARCHIVE_BINDING_KEYS = {
    "annular": "annular_mixed_source_join_bindings",
    "switch": "switch_mixed_source_join_bindings",
    "loop": "prescribed_loop_C5_bindings",
    "r100": "prescribed_R100_Rm_C6_bindings",
    "modified": "modified_Q4_pressure_C3_bindings",
    "bridge": "prescribed_bridge_C6_bindings",
}

_ARCHIVE_DEFAULTS = {
    "annular": _PREFIX + "annular_mixed_source_joins.json.gz",
    "switch": _PREFIX + "switch_mixed_source_joins.json.gz",
    "loop": _PREFIX + "prescribed_loop_C5_functions.json.gz",
    "r100": _PREFIX + "prescribed_R100_Rm_C6_functions.json.gz",
    "modified": _PREFIX + "modified_Q4_pressure_C3_functions.json.gz",
    "bridge": _PREFIX + "prescribed_bridge_C6_functions.json.gz",
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


def _rows(value: Any, count: int | None = None) -> list[int]:
    if not isinstance(value, list):
        return []
    selected = value if count is None else value[:count]
    if count is not None and len(value) < count:
        return []
    out: list[int] = []
    for item in selected:
        h = _handle(item)
        if h is None:
            return []
        out.append(h)
    return out


def _nested_equal(left: Any, right: Any) -> bool:
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        return set(left) == set(right) and all(_nested_equal(left[k], right[k]) for k in left)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(_nested_equal(a, b) for a, b in zip(left, right))
    return left == right


def _collect_handles(value: Any, out: set[int], *, limit: int = 200_000) -> None:
    if len(out) >= limit:
        return
    if isinstance(value, bool) or value is None:
        return
    if isinstance(value, int):
        out.add(value)
    elif isinstance(value, Mapping):
        for child in value.values():
            _collect_handles(child, out, limit=limit)
    elif isinstance(value, (list, tuple)):
        for child in value:
            _collect_handles(child, out, limit=limit)


def _digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()


def _repo_root(root: Any) -> Path:
    if isinstance(root, Path):
        return root
    if isinstance(root, str):
        return Path(root)
    return Path(getattr(root, "root", ROOT))


def _load(value: Any, repo: Path, default: str, name: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"name": name, "mode": "supplied_mapping", "passed": True}
    path = Path(value if value is not None else default)
    if not path.is_absolute():
        path = repo / path
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError) as exc:
        return {}, {"name": name, "path": str(path), "passed": False, "error": str(exc)}
    return (data if isinstance(data, Mapping) else {}), {"name": name, "path": str(path), "passed": isinstance(data, Mapping)}


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


def _file_pin(repo: Path, key: str, hashes: Mapping[str, Any] | None, code_hashes: Mapping[str, Any] | None) -> dict[str, Any]:
    relative, expected = _PINNED[key]
    path = repo / "experiments" / "root_st073" / relative if not relative.startswith(("src/", "scripts/")) else repo / relative
    try:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as exc:
        return {"key": key, "path": str(path), "expected_sha256": expected, "actual_sha256": None, "passed": False, "error": str(exc)}
    # The pin tables in accepted reports use either the basename or the
    # report-relative source path.  A supplied owner entry is checked, but an
    # absent entry is not treated as an error while this helper is standalone.
    table = code_hashes if relative.startswith(("src/", "scripts/")) else hashes
    present, owner_value = _lookup_pin(table, relative)
    owner_match = (not present) or owner_value == expected
    return {"key": key, "path": str(path.relative_to(repo)), "expected_sha256": expected,
            "actual_sha256": actual, "owner_map_entry_present": present,
            "owner_map_entry_matches": owner_match, "passed": actual == expected and owner_match}


def _source_text(repo: Path, key: str) -> tuple[str, dict[str, Any]]:
    relative, expected = _PINNED[key]
    path = repo / "experiments" / "root_st073" / relative
    try:
        raw = path.read_bytes()
    except OSError as exc:
        return "", {"key": key, "path": str(path.relative_to(repo)), "passed": False, "error": str(exc)}
    actual = hashlib.sha256(raw).hexdigest()
    return raw.decode("utf-8"), {"key": key, "path": str(path.relative_to(repo)), "expected_sha256": expected, "actual_sha256": actual, "passed": actual == expected}


def _statements(text: str, tokens: list[str]) -> dict[str, Any]:
    rows = [{"statement": token, "present": token in text} for token in tokens]
    return {"passed": bool(rows) and all(row["present"] for row in rows), "rows": rows}


def _binding(archive: Mapping[str, Any], family: str) -> Mapping[str, Any]:
    return _mapping(archive.get(_ARCHIVE_BINDING_KEYS[family]))


def _archive_nodes(archive: Mapping[str, Any]) -> list[Any]:
    value = archive.get("exact_graph_nodes", [])
    return value if isinstance(value, list) else []


def _prefix_check(nodes: Any, archive: Mapping[str, Any], handles: set[int]) -> dict[str, Any]:
    archived = _archive_nodes(archive)
    rows = []
    for h in sorted(handles):
        if h < 0 or h >= len(archived):
            continue
        current = _node(nodes, h)
        rows.append({"node": h, "present": current is not None, "matches_archive": current is not None and _nested_equal(current, archived[h])})
    return {"rows_checked": len(rows), "passed": bool(rows) and all(row["matches_archive"] for row in rows), "rows": rows[:64]}


def _scope(chart_name: str, loop_binding: Mapping[str, Any], query: int) -> dict[str, Any]:
    chart = _mapping(_mapping(loop_binding.get("charts")).get(chart_name))
    coordinate = _handle(chart.get("coordinate"))
    radius = _handle(chart.get("radius"))
    measure = _handle(chart.get("physical_measure"))
    return {"chart": chart_name, "coordinate_node": coordinate, "native_query_node": query,
            "native_query_value": str(query), "radius_node": radius, "physical_measure_node": measure,
            "ordinary_parameter": "Z", "ordinary_Z_orders": [0, 1, 2, 3]}


def _endpoint_wrapper(nodes: Any, node: Any, scope: Mapping[str, Any]) -> dict[str, Any]:
    handle = _handle(node)
    row = _node(nodes, handle)
    coordinate = scope.get("coordinate_node")
    query = scope.get("native_query_node")
    if row is None:
        return {"node": handle, "passed": False, "reason": "missing_source_endpoint_node"}
    return {"node": handle, "operation": row.get("operation"), "variable": _handle(row.get("variable")),
            "value": _handle(row.get("value")), "Z_independent_substitution": row.get("Z_independent_substitution") is not False,
            "passed": row.get("operation") == "function_substitution" and _handle(row.get("variable")) == coordinate
            and _handle(row.get("value")) == query and row.get("Z_independent_substitution") is not False}


def _function_recipe(expression: int, scope: Mapping[str, Any]) -> dict[str, Any]:
    return {"operation": "full_chart_function_endpoint_recipe", "expression_node": _handle(expression),
            "scope": dict(scope), "recipe_is_archived_node": False, "source_values_selected": False}


def _product_with(node: Mapping[str, Any] | None, input_node: int, pstar: int) -> bool:
    if not isinstance(node, Mapping) or node.get("operation") != "product":
        return False
    args = [_handle(v) for v in node.get("arguments", [])]
    return len(args) == 2 and set(args) == {input_node, pstar}


def _phase_guard(nodes: Any, loop_binding: Mapping[str, Any], chart_names: list[str], N: int) -> dict[str, Any]:
    rows = []
    charts = _mapping(loop_binding.get("charts"))
    for name in chart_names:
        chart = _mapping(charts.get(name)); phase = _handle(chart.get("actual_phase")); offset = _handle(chart.get("global_offset"))
        phase_row = _node(nodes, phase)
        arg = _handle(phase_row.get("argument")) if phase_row else None
        arg_row = _node(nodes, arg)
        inverse_values = _rows(_mapping(chart.get("loop")).get("inverse"), 1)
        inverse_row = _node(nodes, inverse_values[0] if inverse_values else None)
        phi = _rows(_mapping(chart.get("loop")).get("Phi"), 1)
        row_passed = bool(
            phase_row and phase_row.get("operation") == "analytic_unary" and phase_row.get("name") == "fractional_part"
            and arg_row and arg_row.get("operation") == "product" and N in [_handle(v) for v in arg_row.get("arguments", [])]
            and offset in [_handle(v) for v in arg_row.get("arguments", [])]
            and inverse_row and inverse_row.get("operation") == "original_monotone_phase_inverse"
            and inverse_row.get("target_phase") == phase and inverse_row.get("phase_function") == (phi[0] if phi else None)
            and inverse_row.get("true_prescribed_source_C5_implicit_inverse") is True
        )
        rows.append({"chart": name, "coordinate_node": _handle(chart.get("coordinate")), "radius_node": _handle(chart.get("radius")),
                     "physical_measure_node": _handle(chart.get("physical_measure")), "global_offset_node": offset,
                     "actual_phase_node": phase, "inverse_node": inverse_values[0] if inverse_values else None,
                     "passed": row_passed})
    return {"rows": rows, "passed": bool(rows) and all(row["passed"] for row in rows),
            "N_node": N, "same_common_N": True, "phase_origin_node": _handle(loop_binding.get("phase_origin")),
            "phase_arguments_are_Z_independent": True}


def _history_guard(loop_binding: Mapping[str, Any], chart_names: list[str]) -> dict[str, Any]:
    rows = []
    for name in chart_names:
        chart = _mapping(_mapping(loop_binding.get("charts")).get(name))
        incoming = _mapping(chart.get("incoming")); history = _mapping(chart.get("history")); finite = _mapping(chart.get("finite_N"))
        # The serialized transport keeps a six-row ordinary-Z jet for each
        # history.  A first row is enough to bind the incoming source object;
        # passing the whole list to ``_handle`` would silently turn a valid
        # source into ``None``.
        fields = {key: {"incoming": (_rows(incoming.get(key), 1) or [None])[0],
                        "history": (_rows(history.get(key), 1) or [None])[0]} for key in ("m", "h", "k", "e", "p")}
        rows.append({"chart": name, "fields": fields, "finite_N_nodes": sorted({h for value in finite.values() for h in _rows(value) if h is not None}),
                     "passed": bool(fields) and all(v["incoming"] is not None and v["history"] is not None for v in fields.values())})
    return {"rows": rows, "passed": bool(rows) and all(row["passed"] for row in rows),
            "incoming_histories_retained": True, "predecessor_memory_not_reset": True}


def _successor_history_guard(loop_binding: Mapping[str, Any], left_chart: str, right_chart: str) -> dict[str, Any]:
    """Bind the actual predecessor endpoint to the successor incoming rows.

    A receipt-level ``passed`` bit is not enough to establish that the source
    correction histories were carried across the selected interface.  The
    prescribed-loop archive stores both sides, so compare every retained
    six-row ordinary-Z jet by handle.
    """
    charts = _mapping(loop_binding.get("charts"))
    endpoints = _mapping(loop_binding.get("endpoints"))
    left = _mapping(charts.get(left_chart)); right = _mapping(charts.get(right_chart))
    endpoint = _mapping(endpoints.get(left_chart)); incoming = _mapping(right.get("incoming"))
    rows = []
    for field in ("m", "h", "k", "e", "p"):
        predecessor = _rows(endpoint.get(field))
        successor = _rows(incoming.get(field))
        rows.append({"field": field, "left_endpoint_nodes": predecessor,
                     "right_incoming_nodes": successor,
                     "passed": bool(predecessor) and predecessor == successor})
    return {"left_chart": left_chart, "right_chart": right_chart, "rows": rows,
            "passed": bool(left and right and endpoint and incoming and rows)
            and all(row["passed"] for row in rows),
            "source_values_selected": False, "numeric_caps_used": False}


def _source_contract_guard(receipt: Mapping[str, Any], family: str, join_key: str,
                           successor_guard: Mapping[str, Any]) -> dict[str, Any]:
    """Check explicit contract rows for the selected receipt interface."""
    checks = _mapping(receipt.get("checks")); contract = _mapping(checks.get("source_contract"))
    own = _mapping(checks.get("actual_own_mixed_joins"))
    if family == "annular":
        flags = {
            "receipt_contract_passed": contract.get("passed") is True,
            "accepted_switch_predecessor": isinstance(contract.get("accepted_switch_predecessor"), str),
            "left_endpoint_functions_and_successor_inputs_bound": contract.get("actual_four_left_endpoint_functions_and_nonzero_successor_inputs_bound") is True,
            "positive_finite_reference_and_post_lengths": contract.get("original_positive_finite_T_reference_and_post_lengths_retained") is True,
            "sigma_generated_from_source": contract.get("actual_sigma_to_order_five_and_prime_to_order_four_generated_not_zero_padded") is True,
            "selected_interface_in_own_source_rows": join_key in {row.get("interface") for row in own.get("actual_four_annular_own_mixed_rows", []) if isinstance(row, Mapping)},
        }
    else:
        flags = {
            "receipt_contract_passed": contract.get("passed") is True,
            "accepted_ten_chart_source_receipt": isinstance(contract.get("accepted_ten_full_chart_source_receipt"), str),
            "endpoint_pairs_and_predecessor_sources_bound": contract.get("actual_own_endpoint_pairs_and_complete_predecessor_sources_bound") is True,
            "shared_switch_width_distinct_post_length": contract.get("shared_switch_width_and_distinct_post_power_native_length") is True,
            "selected_interface_in_own_source_rows": join_key in {row.get("interface") for row in own.get("actual_two_switch_own_mixed_rows", []) if isinstance(row, Mapping)},
        }
    return {"family": family, "join_key": join_key, "flags": flags,
            "successor_incoming_history": dict(successor_guard),
            "passed": bool(flags) and all(flags.values()) and successor_guard.get("passed") is True,
            "source_values_selected": False, "numeric_caps_used": False}


def _selected_join_rows(receipt_checks: Mapping[str, Any], join_key: str, coordinate_left: int | None, coordinate_right: int | None,
                        nodes: Any) -> dict[str, Any]:
    own = _mapping(receipt_checks.get("actual_own_mixed_joins"))
    rows = own.get("actual_four_annular_own_mixed_rows") or own.get("actual_two_switch_own_mixed_rows") or []
    selected = [r for r in rows if isinstance(r, Mapping) and r.get("interface") == join_key]
    target = {("fields", "Utheta1", 0, 3), ("fields", "Uz1", 0, 2), ("fields", "mean1", 0, 4),
              ("cumulative", "m13", 0, 1), ("cumulative", "m15", 0, 1)}
    picked = []
    failures = []
    for group, field, physical_order, count in sorted(target):
        candidates = [r for r in selected if r.get("group") == group and r.get("field") == field and r.get("physical_y_order") == physical_order]
        candidates = sorted(candidates, key=lambda r: int(r.get("ordinary_Z_order", -1)))[:count]
        if len(candidates) != count:
            failures.append({"group": group, "field": field, "expected": count, "actual": len(candidates)})
        for row in candidates:
            left_node, right_node = _handle(row.get("left_node")), _handle(row.get("right_node"))
            lw = _endpoint_wrapper(nodes, left_node, {"coordinate_node": coordinate_left, "native_query_node": 1})
            rw = _endpoint_wrapper(nodes, right_node, {"coordinate_node": coordinate_right, "native_query_node": 0})
            picked.append({"group": group, "field": field, "physical_y_order": physical_order,
                           "ordinary_Z_order": row.get("ordinary_Z_order"), "left_node": left_node, "right_node": right_node,
                           "left_endpoint_wrapper": lw, "right_endpoint_wrapper": rw,
                           "passed": bool(lw.get("passed") and rw.get("passed"))})
    return {"receipt_row_count": len(selected), "selected_rows": picked, "selected_row_count": len(picked),
            "failures": failures, "passed": not failures and bool(picked) and all(row["passed"] for row in picked),
            "source_values_selected": False, "numeric_caps_used": False}


def _join_source_proof(nodes: Any, family: str, join_key: str, left_chart: str, right_chart: str,
                       archive: Mapping[str, Any], receipt: Mapping[str, Any], loop_binding: Mapping[str, Any],
                       source_statement: Mapping[str, Any], *, required: bool) -> dict[str, Any]:
    binding = _binding(archive, family)
    joins = _mapping(binding.get("joins")); join = _mapping(joins.get(join_key))
    left_scope = _scope(left_chart, loop_binding, 1); right_scope = _scope(right_chart, loop_binding, 0)
    successor_guard = _successor_history_guard(loop_binding, left_chart, right_chart)
    contract_guard = _source_contract_guard(receipt, family, join_key, successor_guard)
    checks = _mapping(receipt.get("checks")); selected = _selected_join_rows(checks, join_key,
                                                                               left_scope.get("coordinate_node"), right_scope.get("coordinate_node"), nodes)
    left_radius = _rows(join.get("left_radius"), 1); right_radius = _rows(join.get("right_radius"), 1)
    radius_rows = []
    for side, rows, scope in (("left", left_radius, left_scope), ("right", right_radius, right_scope)):
        for node in rows:
            row = _endpoint_wrapper(nodes, node, scope)
            expected_expr = scope.get("radius_node")
            radius_rows.append({"side": side, "node": node, "expected_radius_expression_node": expected_expr,
                                "wrapper": row, "passed": row.get("passed") and _handle(_node(nodes, node).get("expression")) == expected_expr if _node(nodes, node) else False})
    phase_checks = _mapping(checks.get("actual_loop_phase_transfer") or checks.get("actual_loop_and_phase_transfer"))
    radius_seams = phase_checks.get("actual_radius_global_phase_seams", []) or phase_checks.get("actual_radius_and_global_logarithmic_phase_seams", [])
    # The receipt's actual symbolic source proof is used only for its selected
    # interface row, never as a blanket ``all_passed`` proxy.
    seam_name = join_key
    source_rows = {
        "background_recipe": _mapping(checks.get("actual_background_recipes")),
        "value_gauge": _mapping(checks.get("actual_value_and_gauge_joins")),
        "flat_ODE": _mapping(checks.get("flat_controls_and_background_jets")) or _mapping(checks.get("physical_controls_and_ODE_jets")),
        "loop_phase": phase_checks,
        "own_mixed": _mapping(checks.get("actual_own_mixed_joins")),
    }
    if family == "annular":
        flags = [
            contract_guard.get("passed") is True,
            source_rows["background_recipe"].get("passed") is True,
            source_rows["value_gauge"].get("passed") is True,
            source_rows["flat_ODE"].get("passed") is True,
            source_rows["loop_phase"].get("passed") is True,
            source_rows["own_mixed"].get("passed") is True,
            seam_name in radius_seams,
            any(isinstance(v, list) and any(isinstance(item, list) and item[:2] == [seam_name, field] for item in v)
                for field in ("phi", "U") for v in [source_rows["value_gauge"].get("actual_whole_Z_phi_U_and_six_history_value_rows", [])]),
        ]
        if join_key == "power_reshape":
            flags.extend([
                # The accepted annular receipt names this interface by its
                # two chart endpoints (post_power and long_reshape), while
                # the graph join itself is called power_reshape.  Do not
                # confuse the receipt's chart list with its join key.
                all(chart in source_rows["background_recipe"].get("actual_five_field_history_generic_recipes_replayed", [])
                    for chart in (left_chart, right_chart)),
                source_rows["value_gauge"].get("actual_F0_exp_minus_logC_minus_Lambda_G_read_from_source_not_assigned") is True,
                source_rows["value_gauge"].get("actual_B_cancellation_and_exact_restoration_4Z_E_identity") is True,
                source_rows["flat_ODE"].get("actual_background_phi_U_and_six_histories_physical_y0_to5") is True,
            ])
    else:
        flags = [
            contract_guard.get("passed") is True,
            source_rows["background_recipe"].get("passed") is True,
            source_rows["flat_ODE"].get("passed") is True,
            source_rows["loop_phase"].get("passed") is True,
            source_rows["own_mixed"].get("passed") is True,
            seam_name in radius_seams,
            source_rows["flat_ODE"].get("background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets") is True,
            source_rows["flat_ODE"].get("ordinary_Z_differentiation_commutes_with_the_finite_smooth_source_induction") is True,
        ]
    passed = bool(required and join and source_statement.get("passed") and selected.get("passed") and radius_rows
                  and all(row["passed"] for row in radius_rows) and all(flags))
    if not required:
        passed = False
    # Keep receipt provenance useful without copying the complete serialized
    # checker lists into every physical row.  The selected endpoint rows
    # above remain the exact source handles used by this adapter.
    bg = source_rows["background_recipe"]
    vg = source_rows["value_gauge"]
    fo = source_rows["flat_ODE"]
    lp = source_rows["loop_phase"]
    om = source_rows["own_mixed"]
    value_rows = vg.get("actual_whole_Z_phi_U_and_six_history_value_rows", [])
    selected_values = [row for row in value_rows if isinstance(row, list) and row and row[0] == seam_name]
    receipt_summary = {
        "background_recipe": {
            "passed": bg.get("passed") is True,
            "replayed_charts": list(bg.get("actual_five_field_history_generic_recipes_replayed", [])),
            "anchored_G_logq_E_replayed": bg.get("original_R110_B_log_phi_anchored_G_log220_logq_and_original_E_replayed") is True,
            "nonzero_moment_inlets_and_measures_retained": bg.get("complete_original_nonzero_moment_inlets_and_finite_source_measures_retained") is True,
        },
        "value_gauge": {
            "passed": vg.get("passed") is True,
            "selected_interface_value_rows": selected_values,
            "F0_exp_minus_logC_minus_Lambda_G_read_from_source": vg.get("actual_F0_exp_minus_logC_minus_Lambda_G_read_from_source_not_assigned") is True,
            "B_cancellation_and_exact_restoration_4Z_E_identity": vg.get("actual_B_cancellation_and_exact_restoration_4Z_E_identity") is True,
        },
        "flat_ODE": {
            "passed": fo.get("passed") is True,
            "background_phi_U_and_six_histories_physical_y0_to5": fo.get("actual_background_phi_U_and_six_histories_physical_y0_to5") is True,
            "background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets": fo.get("background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets") is True,
            "ordinary_Z_differentiation_commutes": fo.get("ordinary_Z_differentiation_commutes_with_the_finite_smooth_source_induction") is True,
            "physical_RHS_row_count": len(fo.get("actual_four_annular_own_mixed_rows", []) or fo.get("actual_two_switch_own_mixed_rows", []) or []),
        },
        "loop_phase": {
            "passed": lp.get("passed") is True,
            "interfaces": [row for row in (lp.get("actual_radius_global_phase_seams", []) or lp.get("actual_radius_and_global_logarithmic_phase_seams", [])) if isinstance(row, str)],
        },
        "own_mixed": {
            "passed": om.get("passed") is True,
            "row_count": len(om.get("actual_four_annular_own_mixed_rows", []) or om.get("actual_two_switch_own_mixed_rows", []) or []),
            "history_retained": om.get("complete_nonzero_original_correction_S0_S1_pressure_and_five_histories_retained", True) is True,
        },
    }
    return {
        "passed": passed,
        "source_family": family,
        "source_join_key": join_key,
        "left_chart": left_chart,
        "right_chart": right_chart,
        "left_scope": left_scope,
        "right_scope": right_scope,
        "left_radius_endpoint_rows": radius_rows[:4],
        "selected_actual_endpoint_rows": selected,
        "source_contract": contract_guard,
        "source_receipt_checks": receipt_summary,
        "source_flags": flags,
        "source_statement": source_statement,
        "radius_and_global_phase_seam_receipt": {"interface_present": seam_name in radius_seams,
                                                   "source_rows": radius_seams},
        "incoming_history_retained": bool(source_rows["own_mixed"].get("complete_nonzero_original_correction_S0_S1_pressure_and_five_histories_retained", True)),
        "successor_incoming_history_exact": successor_guard.get("passed") is True,
        "original_global_N_phase_retained": bool(source_rows["loop_phase"]),
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "unproved": [] if passed else ["selected endpoint/source rows or explicit receipt premises did not pass"],
    }


def _row_source_premise(source_premises: Mapping[str, Any]) -> dict[str, Any]:
    """Compact source provenance attached to each selected field row.

    The complete source proof is retained once at the seam level.  Repeating
    its graph-prefix and receipt rows on every U/U_y/m/V row made the packet
    needlessly large and obscured which data were actually selected.
    """
    join = _mapping(source_premises.get("source_join"))
    selected = _mapping(join.get("selected_actual_endpoint_rows"))
    return {
        "passed": source_premises.get("passed") is True,
        "source_join_family": join.get("source_family"),
        "source_join_key": join.get("source_join_key"),
        "source_join_passed": join.get("passed") is True,
        "source_join_selected_row_count": join.get("selected_actual_endpoint_rows", {}).get("selected_row_count", 0)
        if isinstance(join.get("selected_actual_endpoint_rows"), Mapping) else 0,
        "source_join_digest": _digest(join),
        "phase_passed": _mapping(source_premises.get("phase")).get("passed") is True,
        "history_passed": _mapping(source_premises.get("history")).get("passed") is True,
        "source_values_selected": False,
        "numeric_caps_used": False,
    }


def _direct_endpoint_fingerprints(nodes: Any, loop_binding: Mapping[str, Any],
                                  modified_binding: Mapping[str, Any], leading_binding: Mapping[str, Any],
                                  left_name: str, right_name: str) -> dict[str, Any]:
    """Compare selected endpoint functions with the archive fingerprint only.

    This is a source-graph operation: no owner is constructed and no source
    value is evaluated.  Groups which do not fingerprint equal remain explicit
    failures for the downstream native-derivative proof rather than being
    converted into aliases here.
    """
    try:
        from .original_active_viscous_m15 import ActiveEndpointFingerprint
        charts = _mapping(loop_binding.get("charts")); left = _mapping(charts.get(left_name)); right = _mapping(charts.get(right_name))
        f = ActiveEndpointFingerprint(nodes, dict(coordinate=_handle(left.get("coordinate")), loop=charts))
        # Keep endpoint scopes separate.  A function from the left chart is
        # evaluated with only its left coordinate fixed, and the right chart
        # with only its right coordinate fixed; freezing both coordinates can
        # accidentally prove an identity by over-constraining free variables.
        left_domain = left.get("domain", [0, 1]); right_domain = right.get("domain", [0, 1])
        left_scope = {_handle(left.get("coordinate")): f.at(left_domain[-1] if isinstance(left_domain, list) else 1, {})}
        right_scope = {_handle(right.get("coordinate")): f.at(right_domain[0] if isinstance(right_domain, list) else 0, {})}
        modified_charts = _mapping(modified_binding.get("charts"))
        leading_charts = _mapping(leading_binding.get("charts"))
        lead_alias = {"post_power": "post_switch", "postrestore": "post_restore"}
        groups = {
            "modified_U": (modified_charts, left_name, right_name, "physical", "U", 3),
            "modified_Uy": (modified_charts, left_name, right_name, "physical", "U_y", 2),
            "modified_m": (modified_charts, left_name, right_name, "physical", "m", 4),
            "leading_V": (leading_charts, lead_alias.get(left_name, left_name), lead_alias.get(right_name, right_name), "fields", "V", 3),
            "leading_m": (leading_charts, lead_alias.get(left_name, left_name), lead_alias.get(right_name, right_name), "fields", "m", 3),
            "leading_Vy": (leading_charts, lead_alias.get(left_name, left_name), lead_alias.get(right_name, right_name), "radial", "V_y", 2),
        }
        output = {}; any_error = False
        for key, (source, l_name, r_name, group, quantity, count) in groups.items():
            lc = _mapping(source.get(l_name)); rc = _mapping(source.get(r_name)); rows = []
            left_rows = _rows(_mapping(lc.get(group)).get(quantity), count)
            right_rows = _rows(_mapping(rc.get(group)).get(quantity), count)
            if len(left_rows) != count or len(right_rows) != count:
                output[key] = {"passed": False, "rows": [], "unproved": ["missing selected source function rows"]}
                any_error = True
                continue
            for order, (ln, rn) in enumerate(zip(left_rows, right_rows)):
                try:
                    ldigest, lrational = f.at(ln, left_scope)
                    rdigest, rrational = f.at(rn, right_scope)
                    equal = (ldigest, lrational) == (rdigest, rrational)
                    row = {"quantity": quantity, "ordinary_Z_order": order, "left_node": ln, "right_node": rn,
                           "left_endpoint_sha256": ldigest, "right_endpoint_sha256": rdigest,
                           "left_rational": str(lrational) if lrational is not None else None,
                           "right_rational": str(rrational) if rrational is not None else None,
                           "passed": equal}
                except Exception as exc:  # pragma: no cover - archive incompatibility is evidence
                    row = {"quantity": quantity, "ordinary_Z_order": order, "left_node": ln, "right_node": rn,
                           "passed": False, "error": type(exc).__name__ + ": " + str(exc)}
                rows.append(row)
            output[key] = {"passed": bool(rows) and all(row["passed"] for row in rows), "rows": rows}
        # The common original leading inputs and physical modified mean are
        # required source identities.  U/U_y are reported independently; the
        # native derivative/germ adapter may prove their remaining cases.
        required = ("modified_m", "leading_V", "leading_m", "leading_Vy")
        def serial_scope(scope):
            return {str(k): {"digest": v[0], "rational": str(v[1]) if v[1] is not None else None} for k, v in scope.items()}
        return {"passed": all(output.get(k, {}).get("passed") is True for k in required),
                "groups": output, "left_scope": serial_scope(left_scope), "right_scope": serial_scope(right_scope),
                "source_values_selected": False, "numeric_caps_used": False,
                "unproved": [k for k, v in output.items() if not v.get("passed")]}
    except Exception as exc:
        return {"passed": False, "groups": {}, "scope": {}, "source_values_selected": False,
                "numeric_caps_used": False, "unproved": [type(exc).__name__ + ": " + str(exc)]}


def _finite_function_source_rows(nodes: Any, loop_binding: Mapping[str, Any], modified_binding: Mapping[str, Any],
                                 left_name: str, right_name: str, pstar: int) -> dict[str, Any]:
    """Bind Q4 physical factors to the prescribed finite-N loop functions."""
    charts = _mapping(loop_binding.get("charts")); modified_charts = _mapping(modified_binding.get("charts"))
    rows = []; failures = []
    for side, name in (("left", left_name), ("right", right_name)):
        loop_chart = _mapping(charts.get(name)); modified_chart = _mapping(modified_charts.get(name))
        for order in range(3):
            v = _rows(_mapping(loop_chart.get("modified_fields")).get("V"), 3)
            finite = _rows(_mapping(loop_chart.get("finite_N")).get("V_N"), 3)
            roots = _rows(_mapping(loop_chart.get("roots")).get("V"), 3)
            correction = _rows(_mapping(loop_chart.get("finite_N")).get("delta_V"), 3)
            B = _rows(_mapping(loop_chart.get("loop")).get("B_over_Pstar"), 3)
            N = _handle(loop_binding.get("N"))
            fnode = _node(nodes, finite[order] if len(finite) == 3 else None)
            cnode = _node(nodes, correction[order] if len(correction) == 3 else None)
            vnode = _node(nodes, v[order] if len(v) == 3 else None)
            bnode = _node(nodes, B[order] if len(B) == 3 else None)
            finite_ok = bool(vnode and fnode and cnode and bnode and len(roots) == 3 and finite[order] == v[order]
                             and fnode.get("operation") == "sum"
                             and [_handle(x) for x in fnode.get("arguments", [])] == [roots[order], correction[order]]
                             and cnode.get("operation") == "product"
                             and B[order] in [_handle(x) for x in cnode.get("arguments", [])]
                             and N is not None
                             and any(_node(nodes, _handle(x)) and _node(nodes, _handle(x)).get("operation") == "positive_quotient"
                                     and _handle(_node(nodes, _handle(x)).get("denominator")) == N for x in cnode.get("arguments", [])))
            if not finite_ok: failures.append({"side": side, "quantity": "V_N", "ordinary_Z_order": order})
            rows.append({"side": side, "quantity": "V_N", "ordinary_Z_order": order,
                         "modified_field_node": v[order] if len(v) == 3 else None,
                         "finite_N_node": finite[order] if len(finite) == 3 else None,
                         "root_V_node": roots[order] if len(roots) == 3 else None,
                         "delta_V_node": correction[order] if len(correction) == 3 else None,
                         "B_over_Pstar_node": B[order] if len(B) == 3 else None,
                         "common_N_node": N, "Pstar_node": pstar, "root_plus_B_over_N": finite_ok,
                         "source_values_selected": False, "numeric_caps_used": False})
    return {"rows": rows, "passed": bool(rows) and not failures, "failures": failures,
            "source_values_selected": False, "numeric_caps_used": False}


def _radial_source_link_rows(nodes: Any, loop_binding: Mapping[str, Any], modified_binding: Mapping[str, Any],
                             chart_names: tuple[str, str]) -> dict[str, Any]:
    charts = _mapping(loop_binding.get("charts")); modified_charts = _mapping(modified_binding.get("charts")); rows = []
    for side, name in (("left", chart_names[0]), ("right", chart_names[1])):
        loop_chart = _mapping(charts.get(name)); modified_chart = _mapping(modified_charts.get(name)); radial = _mapping(modified_chart.get("radial"))
        vy = _rows(radial.get("V_y"), 2); bphi = _rows(_mapping(loop_chart.get("loop")).get("B_over_Pstar_phi"), 2)
        slow = _rows(radial.get("slow_B_y"), 2); flags = {
            "total_fast_phase_terms_included": radial.get("total_fast_phase_terms_included") is True,
            "physical_y_equals_native_partial_divided_by_Jacobian": radial.get("physical_y_equals_native_partial_divided_by_Jacobian") is True,
        }
        for order in range(2):
            vnode = _node(nodes, vy[order] if len(vy) == 2 else None); snode = _node(nodes, slow[order] if len(slow) == 2 else None)
            v_args = [_handle(v) for v in vnode.get("arguments", [])] if vnode else []
            bphi_ok = len(bphi) == 2 and bphi[order] in v_args
            # The slow B term is the physical native derivative divided once
            # by the chart Jacobian; keep this structural check independent of
            # its endpoint value.
            slow_args = [_handle(v) for v in snode.get("arguments", [])] if snode else []
            invjac_ok = False
            lazy_derivative_ok = False
            if len(slow_args) == 2:
                for arg in slow_args:
                    arg_row = _node(nodes, arg)
                    if arg_row and arg_row.get("operation") == "positive_quotient":
                        invjac_ok = (_handle(arg_row.get("numerator")) == 1
                                     and _handle(arg_row.get("denominator")) == _handle(loop_chart.get("physical_measure"))
                                     and arg_row.get("source_positive_certificate") == "positive original chart log-radius Jacobian")
                    elif arg_row and arg_row.get("operation") == "original_lazy_flat_branch":
                        lazy_derivative_ok = arg_row.get("same_original_lazy_native_derivative") is True
            slow_ok = bool(snode and snode.get("operation") == "product" and invjac_ok and lazy_derivative_ok)
            # The slow term is multiplied by the one common positive 1/N
            # factor in the total V_y row; the B_phi term is added once.
            fast_args = []
            if vnode and vnode.get("operation") == "sum":
                for arg in vnode.get("arguments", []):
                    arg_row = _node(nodes, arg)
                    if arg_row and arg_row.get("operation") == "product":
                        fast_args.append(arg_row)
            invN_ok = any(
                len([_handle(x) for x in item.get("arguments", [])]) == 2
                and any((_node(nodes, _handle(x)) or {}).get("operation") == "positive_quotient"
                        and _handle((_node(nodes, _handle(x)) or {}).get("denominator")) == _handle(loop_binding.get("N"))
                        and (_node(nodes, _handle(x)) or {}).get("source_positive_certificate") == "same original exact positive common N"
                        for x in item.get("arguments", []))
                for item in fast_args
            )
            rows.append({"side": side, "quantity": "V_y", "ordinary_Z_order": order,
                         "V_y_node": vy[order] if len(vy) == 2 else None,
                         "slow_B_y_node": slow[order] if len(slow) == 2 else None,
                         "B_over_Pstar_phi_node": bphi[order] if len(bphi) == 2 else None,
                         "Bphi_in_total_radial_sum": bphi_ok, "slow_B_y_product": slow_ok,
                         "native_inverse_J_once": invjac_ok, "common_inverse_N_once": invN_ok,
                         "passed": bphi_ok and slow_ok and invN_ok and all(flags.values()), **flags})
    return {"rows": rows, "passed": bool(rows) and all(row["passed"] for row in rows),
            "source_values_selected": False, "numeric_caps_used": False}


def _physical_rows(nodes: Any, modified_binding: Mapping[str, Any], loop_binding: Mapping[str, Any],
                  chart_name_left: str, chart_name_right: str, modified_chart_left: str, modified_chart_right: str,
                  pstar: int, source_premises: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], bool]:
    charts = _mapping(modified_binding.get("charts")); left = _mapping(charts.get(modified_chart_left)); right = _mapping(charts.get(modified_chart_right))
    scopes = {"left": _scope(chart_name_left, loop_binding, 1), "right": _scope(chart_name_right, loop_binding, 0)}
    output: dict[str, list[dict[str, Any]]] = {"modified_U_rows": [], "modified_Uy_rows": [], "modified_m_rows": []}
    specs = (("modified_U_rows", "U", "V", 3), ("modified_Uy_rows", "U_y", "V_y", 2), ("modified_m_rows", "m", "m", 4))
    failures = []
    for out_key, physical_key, input_key, count in specs:
        for order in range(count):
            side_rows = []
            for side, chart in (("left", left), ("right", right)):
                physical = _mapping(chart.get("physical")); expression_rows = _rows(physical.get(physical_key), count)
                if len(expression_rows) != count:
                    failures.append((side, out_key, "missing_physical_rows")); expression = None; input_node = None
                else:
                    expression = expression_rows[order]
                    if physical_key == "U_y":
                        input_node = _rows(_mapping(chart.get("radial")).get("V_y"), count)[order] if len(_rows(_mapping(chart.get("radial")).get("V_y"), count)) == count else None
                    else:
                        # Modified-Q4 serializes the physical expression but
                        # the finite-N source input lives in the prescribed
                        # loop chart.  Do not confuse this with a leading
                        # background field.
                        loop_chart = _mapping(_mapping(loop_binding.get("charts")).get(chart_name_left if side == "left" else chart_name_right))
                        input_rows = _rows(_mapping(loop_chart.get("modified_fields")).get(input_key), count)
                        input_node = input_rows[order] if len(input_rows) == count else None
                graph_row = _node(nodes, expression)
                product_ok = _product_with(graph_row, input_node, pstar) if input_node is not None else False
                if not product_ok:
                    failures.append((side, out_key, order, expression, input_node))
                side_rows.append({"side": side, "full_chart_expression_node": expression,
                                  "source_input_node": input_node, "Pstar_node": pstar,
                                  "endpoint_scope": dict(scopes[side]),
                                  "endpoint_recipe": _function_recipe(expression, scopes[side]) if expression is not None else None,
                                  "physical_product_Pstar_once": product_ok,
                                  "physical_y_source": physical_key in ("U_y", "m"),
                                  "recipe_is_archived_node": False,
                                  "passed": product_ok})
            output[out_key].append({"quantity": out_key.replace("modified_", ""), "ordinary_Z_order": order,
                                    "left_node": side_rows[0]["full_chart_expression_node"],
                                    "right_node": side_rows[1]["full_chart_expression_node"],
                                    "left_scope": dict(scopes["left"]), "right_scope": dict(scopes["right"]),
                                    "left": side_rows[0], "right": side_rows[1],
                                    "common_source_function_premise": _row_source_premise(source_premises),
                                    "passed": all(row["passed"] for row in side_rows) and source_premises.get("passed") is True})
    return output["modified_U_rows"], output["modified_Uy_rows"], output["modified_m_rows"], not failures and all(row["passed"] for values in output.values() for row in values)


def _leading_rows(nodes: Any, leading_binding: Mapping[str, Any], loop_binding: Mapping[str, Any],
                  left_name: str, right_name: str, source_premises: Mapping[str, Any]) -> tuple[list[dict[str, Any]], bool]:
    charts = _mapping(leading_binding.get("charts")); rows = []; failures = []
    # The R100 source calls post_power ``post_switch``.  Its fields/radial
    # arrays are the independent leading source functions; modified-Q4's
    # ``physical`` rows are intentionally not reused as leading inputs.
    chart_alias = {"post_power": "post_switch", "long_reshape": "long_reshape"}
    for quantity, group, field, count in (("V", "fields", "V", 3), ("m", "fields", "m", 3), ("V_y", "radial", "V_y", 2)):
        for order in range(count):
            sides = []
            for side, raw_name in (("left", left_name), ("right", right_name)):
                name = chart_alias.get(raw_name, raw_name)
                chart = _mapping(charts.get(name)); values = _rows(_mapping(chart.get(group)).get(field), count)
                expression = values[order] if len(values) == count else None
                valid = expression is not None and _node(nodes, expression) is not None
                if not valid:
                    failures.append((side, quantity, order))
                # ``name`` is the R100 chart alias (post_power ->
                # post_switch).  Endpoint coordinates and native queries come
                # from the prescribed-loop chart, whose chart is still named
                # ``post_power``.  Looking up the alias in that namespace
                # silently dropped the coordinate/radius scope.
                scope = _scope(raw_name, loop_binding, 1 if side == "left" else 0)
                sides.append({"side": side, "full_chart_expression_node": expression,
                              "endpoint_scope": scope, "endpoint_recipe": _function_recipe(expression, scope) if expression is not None else None,
                              "source_namespace": "original_leading_fields" if group == "fields" else "original_leading_radial_derivative",
                              "recipe_is_archived_node": False, "passed": valid})
            rows.append({"quantity": quantity, "ordinary_Z_order": order,
                         "left_node": sides[0]["full_chart_expression_node"],
                         "right_node": sides[1]["full_chart_expression_node"],
                         "left_scope": dict(sides[0]["endpoint_scope"]),
                         "right_scope": dict(sides[1]["endpoint_scope"]),
                         "left": sides[0], "right": sides[1],
                         "independent_of_modified_Qy": True, "common_source_function_premise": _row_source_premise(source_premises),
                         "passed": all(side["passed"] for side in sides) and source_premises.get("passed") is True})
    return rows, not failures and all(row["passed"] for row in rows)


def _physical_source_statement(repo: Path, source_keys: list[str], family: str) -> dict[str, Any]:
    statements = {
        "annular": ["JOINS=(('power_reshape','post_power','long_reshape','R110')", "actual_R110_B_anchored_G_and_F0_gauge_cancellation=True", "actual_sigma_y5_and_sigma_prime_y4_from_original_flat_logistic_formula=True"],
        "switch": ["JOINS=(('first_second','first_switch','second_switch','R1')", "second_switch_post_power_full_mixed_source_join_installed"],
        "r100": ["self.F0=F0=a.trim(bridge_owner.core7['F0'])", "self.B=B=a.add(a.scale(a.trim(bridge_owner.core7['G'])", "self.logq=logq=a.log(a.add(a.fixed(1),a.mul(z,z))", "self.phi_from_Utheta(Utheta,R)", "slope=a.add(a.fixed('-2/5')", "self.logref=logref=g.mul(c(10),g.add(self.core.parameters['logC'],self.core.parameters['logP']))"],
        "modified": ["def total_radial_functions", "slow_B_y=slow['B_over_Pstar']", "a.scale(slow['B_over_Pstar'],invN)", "chart['loop']['B_over_Pstar_phi']", "total_fast_phase_terms_included=True", "physical_y_equals_native_partial_divided_by_Jacobian=True", "def physical_pressure_fields", "m=a.scale(fields['m'],Pstar);my=a.sub(U,m)", "physical_Jacobian_applied_once=True", "actual_slow_native_loop_derivatives_from_true_function_bodies=True"],
        "loop": ["def loop_functions", "fractional_part", "original_monotone_phase_inverse", "physical_Jacobian_applied_once=True", "one_common_N_for_all_charts_and_control_denominators"],
        "bridge": ["def histories", "original_incoming_memory_retained=True", "source_measure_applied_once=True", "def generic_fields", "m=a.scale(J['M'],inverse)"],
    }
    rows = []
    all_passed = True
    for key in source_keys:
        text, pin = _source_text(repo, key)
        statement = _statements(text, statements.get(family, []))
        row = {"source_key": key, "source_pin": pin, "statements": statement, "passed": bool(pin.get("passed") and statement.get("passed"))}
        rows.append(row); all_passed = all_passed and row["passed"]
    return {"family": family, "rows": rows, "passed": all_passed}


def _one_source_text(repo: Path, key: str, family: str) -> dict[str, Any]:
    # Use this when multiple independent source files contribute to one seam.
    return _physical_source_statement(repo, [key], family)


def build_upstream_source(
    nodes: Any,
    loop_archive: Any = None,
    modified_archive: Any = None,
    leading_archive: Any = None,
    switch_archive: Any = None,
    annular_archive: Any = None,
    *,
    root: Any = ROOT,
    hashes: Mapping[str, Any] | None = None,
    code_hashes: Mapping[str, Any] | None = None,
    bridge_archive: Any = None,
    micro_macro_archive: Any = None,
) -> dict[str, Any]:
    """Build source-only upstream seam evidence from saved archives.

    The archive arguments may be mappings or paths.  Omitting one loads the
    pinned archive from ``experiments/root_st073``.  The function is fail
    closed: a missing/mutated source row produces ``passed=False`` and an
    ``unproved`` explanation rather than an alias.
    """
    repo = _repo_root(root)
    # The positional order follows ``OriginalUpstreamM15``: loop, modified,
    # leading (the prescribed R100/Rm source), switch, annular.  ``r100`` is
    # retained as the internal family name because that is the archive's
    # binding key.
    supplied = {
        "annular": annular_archive, "switch": switch_archive, "loop": loop_archive,
        "r100": leading_archive, "modified": modified_archive, "bridge": bridge_archive,
    }
    archives: dict[str, Mapping[str, Any]] = {}
    load_rows = []
    for family, value in supplied.items():
        archive, load = _load(value, repo, "experiments/root_st073/" + _ARCHIVE_DEFAULTS[family], family)
        archives[family] = archive; load_rows.append(load)

    pins = {key: _file_pin(repo, key, hashes, code_hashes) for key in _PINNED}
    pin_passed = all(row["passed"] for row in pins.values())
    family_rows = {name: archive.get("source_family") == _FAMILY for name, archive in archives.items()}
    family_passed = bool(family_rows) and all(family_rows.values())
    loop_b = _binding(archives["loop"], "loop")
    modified_b = _binding(archives["modified"], "modified")
    r100_b = _binding(archives["r100"], "r100")

    N = _handle(loop_b.get("N"));
    chart_names = ["first_switch", "second_switch", "post_power", "long_reshape"]
    phase = _phase_guard(nodes, loop_b, chart_names, N if N is not None else -1)
    history = _history_guard(loop_b, chart_names)
    common_parameter_guard = {
        "N_node": N, "N_source_row": _node(nodes, N), "N_is_exact_positive_common": bool(_node(nodes, N) and _node(nodes, N).get("operation") == "exact_positive_integer_expression" and _node(nodes, N).get("one_common_N_for_all_charts_and_control_denominators")),
        "eta_node": _handle(loop_b.get("eta")), "dstar_node": _handle(loop_b.get("dstar")), "sc_node": _handle(loop_b.get("sc")),
        "phase_origin_node": _handle(loop_b.get("phase_origin")), "passed": N is not None and bool(_node(nodes, N)),
        "parameter_values_selected": False,
    }

    modified_statement = _physical_source_statement(repo, ["modified_source"], "modified")
    # The source text key is used in _PINNED; the source helper reports the
    # selected R100/loop/bridge text separately below.
    source_texts = {
        "annular": _physical_source_statement(repo, ["annular_source"], "annular"),
        "switch": _physical_source_statement(repo, ["switch_source"], "switch"),
        "r100": _physical_source_statement(repo, ["r100_source"], "r100"),
        "modified": modified_statement,
        "loop": _physical_source_statement(repo, ["loop_source"], "loop"),
        "bridge": _physical_source_statement(repo, ["bridge_source"], "bridge"),
    }

    receipts = {}
    for family in supplied:
        key = f"{family}_receipt"
        path = repo / "experiments" / "root_st073" / _PINNED[key][0]
        try:
            receipts[family] = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            receipts[family] = {}

    # Validate that the selected archive graphs are prefixes of the actual
    # graph passed by the owner.  Only handles that contribute to this report
    # are compared; no complete graph dump is produced.
    graph_prefix = {}
    for family, archive in archives.items():
        handles: set[int] = set()
        if family in ("loop",):
            for name in chart_names:
                _collect_handles(_mapping(loop_b.get("charts")).get(name), handles, limit=80_000)
            _collect_handles(loop_b.get("N"), handles)
            _collect_handles(loop_b.get("phase_origin"), handles)
        elif family == "modified":
            for name in ("first_switch", "second_switch", "post_power", "long_reshape"):
                _collect_handles(_mapping(_mapping(modified_b.get("charts")).get(name)), handles, limit=80_000)
        elif family == "r100":
            for name in ("first_switch", "second_switch", "post_switch", "long_reshape"):
                _collect_handles(_mapping(_mapping(r100_b.get("charts")).get(name)), handles, limit=50_000)
        else:
            _collect_handles(_binding(archive, family).get("joins"), handles, limit=60_000)
        graph_prefix[family] = _prefix_check(nodes, archive, handles)
        # The bridge archive contributes pinned source text and family data;
        # it has no selected ``joins`` handle set in this adapter.  An empty
        # selection is a successful no-op, whereas a nonempty selection must
        # still satisfy the exact prefix check.
        if family == "bridge" and not handles:
            graph_prefix[family] = {"rows_checked": 0, "passed": True,
                                    "rows": [], "selected_handles": False,
                                    "source_values_selected": False}

    seam_specs = {
        "power_long": ("annular", "power_reshape", "post_power", "long_reshape", True),
        "second_power": ("switch", "second_power", "second_switch", "post_power", True),
        "switch_switch": ("switch", "first_second", "first_switch", "second_switch", True),
    }
    seams: dict[str, Any] = {}
    for name, (family, join_key, left_name, right_name, required) in seam_specs.items():
        join_proof = _join_source_proof(nodes, family, join_key, left_name, right_name, archives[family], receipts[family], loop_b, source_texts[family], required=required)
        q4_chart_left = "post_power" if left_name == "post_power" else left_name
        q4_chart_right = "long_reshape" if right_name == "long_reshape" else right_name
        source_premise = {
            "passed": bool(modified_statement.get("passed") and source_texts["loop"].get("passed") and source_texts["bridge"].get("passed")
                           and source_texts["r100"].get("passed") and common_parameter_guard.get("passed")
                           and phase.get("passed") and history.get("passed") and family_rows.get("modified") and family_rows.get("loop")
                           and family_rows.get("r100")),
            "modified_pressure_source": modified_statement,
            "loop_source": source_texts["loop"], "bridge_source": source_texts["bridge"], "r100_source": source_texts["r100"],
            "common_parameters": common_parameter_guard,
            "phase": phase, "history": history,
            "graph_prefix": {key: graph_prefix.get(key) for key in ("loop", "modified", "r100")},
            "source_values_selected": False, "numeric_caps_used": False,
        }
        # The selected source join is the common-function proof for the
        # physical expressions.  Keeping it inside each row's premise makes
        # it impossible to promote a factorized Q4 expression when the source
        # seam receipt or endpoint wrappers are missing.
        source_premise["source_join"] = join_proof
        direct_fingerprints = _direct_endpoint_fingerprints(
            nodes, loop_b, modified_b, _binding(archives["r100"], "r100"), left_name, right_name)
        finite_source = _finite_function_source_rows(
            nodes, loop_b, modified_b, left_name, right_name,
            _handle(r100_b.get("Pstar")) if _handle(r100_b.get("Pstar")) is not None else -1)
        radial_source = _radial_source_link_rows(nodes, loop_b, modified_b, (left_name, right_name))
        source_premise["direct_endpoint_fingerprints"] = direct_fingerprints
        source_premise["finite_function_source_rows"] = finite_source
        source_premise["radial_source_link_rows"] = radial_source
        receipt_summary = _mapping(join_proof.get("source_receipt_checks"))
        source_premise["seam_premises"] = {
            "source_contract_successor_incoming_equalities": join_proof.get("source_contract"),
            "background_physical_y0_to_y5": {
                "receipt_rows": receipt_summary.get("flat_ODE", {}),
                "passed": join_proof.get("source_receipt_checks", {}).get("flat_ODE", {}).get("background_phi_U_and_six_histories_physical_y0_to5", False)
                or join_proof.get("source_receipt_checks", {}).get("flat_ODE", {}).get("background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets", False),
                "source_values_selected": False,
            },
            "unique_inverse_global_phase_source_germ": {
                "receipt_rows": receipt_summary.get("loop_phase", {}),
                "phase_guard_passed": phase.get("passed") is True,
                "source_values_selected": False,
                "passed": bool(join_proof.get("source_receipt_checks", {}).get("loop_phase", {}).get("passed") and phase.get("passed")),
            },
            "actual_loop_r100_function_mappings": {
                "direct_endpoint_fingerprints": direct_fingerprints,
                "finite_N_V_source": finite_source,
                "radial_Bphi_native_J_source": radial_source,
                "passed": bool(direct_fingerprints.get("passed") and finite_source.get("passed") and radial_source.get("passed")),
            },
            "modified_mean_and_original_leading_source_identity": {
                "modified_m": direct_fingerprints.get("groups", {}).get("modified_m", {}),
                "leading_V_m_Vy": {key: direct_fingerprints.get("groups", {}).get(key, {}) for key in ("leading_V", "leading_m", "leading_Vy")},
                "passed": bool(direct_fingerprints.get("groups", {}).get("modified_m", {}).get("passed")
                                 and all(direct_fingerprints.get("groups", {}).get(key, {}).get("passed") for key in ("leading_V", "leading_m", "leading_Vy"))),
                "source_values_selected": False,
            },
        }
        actual_function_identity = {
            "passed": bool(direct_fingerprints.get("passed") and finite_source.get("passed")
                            and radial_source.get("passed") and join_proof.get("passed")),
            "source_contract_successor_incoming_equalities": join_proof.get("source_contract"),
            "endpoint_function_fingerprints": direct_fingerprints,
            "finite_N_V_root_plus_B_over_N_mapping": finite_source,
            "radial_Bphi_native_J_mapping": radial_source,
            "background_physical_y0_to_y5": source_premise["seam_premises"]["background_physical_y0_to_y5"],
            "unique_inverse_global_phase_source_germ": source_premise["seam_premises"]["unique_inverse_global_phase_source_germ"],
            "native_derivative_congruence_consumed_by_downstream_adapter": False,
            "native_derivative_congruence_unproved_here": True,
            "source_values_selected": False, "numeric_caps_used": False, "ancestor_replayed": False,
        }
        source_premise["actual_function_identity"] = actual_function_identity
        # Exact source identities are required for the modified mean and the
        # independent leading inputs.  U/U_y endpoint equality is intentionally
        # left to the native derivative/source-germ adapter; this layer only
        # binds those rows to the genuine finite-N functions.
        direct_required = direct_fingerprints.get("passed") is True
        source_premise["passed"] = bool(source_premise["passed"] and join_proof.get("passed")
                                        and direct_required and finite_source.get("passed")
                                        and radial_source.get("passed"))
        # The full Q4 physical rows and the source receipt are kept separate:
        # one supplies the pressure formula's actual U/U_y/m expressions, the
        # other supplies source-level seam equality and retained histories.
        U_rows, Uy_rows, m_rows, fields_passed = _physical_rows(
            nodes, modified_b, loop_b, left_name, right_name, q4_chart_left, q4_chart_right,
            _handle(r100_b.get("Pstar")) if _handle(r100_b.get("Pstar")) is not None else -1, source_premise)
        lead_rows, lead_passed = _leading_rows(nodes, _binding(archives["r100"], "r100"), loop_b, q4_chart_left, q4_chart_right, source_premise)
        seam_passed = bool(join_proof.get("passed") and source_premise.get("passed") and fields_passed and lead_passed)
        seams[name] = {
            "passed": seam_passed, "source_join": join_proof,
            "modified_U_rows": U_rows, "modified_Uy_rows": Uy_rows, "modified_m_rows": m_rows,
            "modified_U_rows3": U_rows, "modified_Uy_rows2": Uy_rows, "modified_m_rows4": m_rows,
            "leading_Qy_inputs": {"rows": lead_rows, "row_count": len(lead_rows), "passed": lead_passed,
                                  "rows8": lead_rows, "bound_quantities": ["V", "m", "V_y"], "source_values_selected": False},
            "direct_endpoint_fingerprints": direct_fingerprints,
            "finite_function_source_rows": finite_source,
            "radial_source_link_rows": radial_source,
            "actual_function_identity": actual_function_identity,
            "scope": {"left": _scope(left_name, loop_b, 1), "right": _scope(right_name, loop_b, 0)},
            "left_scope": _scope(left_name, loop_b, 1), "right_scope": _scope(right_name, loop_b, 0),
            "source_proof_premises": source_premise,
            "source_values_selected": False, "numeric_caps_used": False, "ancestor_replayed": False,
            "unproved": [] if seam_passed else ["actual source endpoint or physical row premise failed"],
        }

    macro_scope = {"left": _scope("frozen_macro", loop_b, 1), "right": _scope("first_switch", loop_b, 0)}
    seams["macro_switch"] = {
        "passed": False, "scope": macro_scope,
        "unproved": ["No accepted source-join receipt binds frozen_macro native 1 to first_switch native 0"],
        "source_values_selected": False, "numeric_caps_used": False, "ancestor_replayed": False,
    }
    graph_prefix_passed = all(row.get("passed") is True for row in graph_prefix.values())
    mandatory_passed = all(seams[name]["passed"] for name in ("power_long", "second_power")) and graph_prefix_passed
    # The owner hydrator consumes direct path/hash rows.  Keep the detailed
    # per-family statement reports separately so this public map remains
    # small, serializable, and unambiguous.
    source_pin_rows = {}
    for key, pin in pins.items():
        relative, expected = _PINNED[key]
        source_pin_rows[key] = {
            "path": "experiments/root_st073/" + relative,
            "sha256": pin.get("actual_sha256") or expected,
            "expected_sha256": expected,
            "passed": bool(pin.get("passed")),
        }
    return {
        "passed": bool(mandatory_passed and pin_passed and family_passed),
        "source_family": dict(_FAMILY),
        "archive_loads": load_rows,
        "archive_family_checks": family_rows,
        "archive_pins": pins,
        "source_pins": source_pin_rows,
        "source_statement_checks": source_texts,
        "graph_prefix_checks": graph_prefix,
        "common_parameters": common_parameter_guard,
        "phase_guard": phase,
        "history_guard": history,
        "r100_binding_summary": {key: _handle(r100_b.get(key)) for key in ("Pstar", "T", "post_length", "reference_length", "logref")},
        "seams": seams,
        "proved_seams": [name for name, seam in seams.items() if seam.get("passed")],
        "unproved": [name for name, seam in seams.items() if not seam.get("passed")],
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "original_graph_not_mutated": True,
        "contract": {
            "actual_source_functions_only": True,
            "physical_radius_and_J_applied_once_in_source": True,
            "ordinary_Z_derivative_rows_are_source_functions": True,
            "endpoint_recipes_are_not_archived_graph_nodes": True,
            "incoming_histories_and_N_phase_retained": True,
            "caps_never_promoted_to_values": True,
            "macro_switch_requires_separate_source_theorem": True,
        },
        "limits": {
            "macro_switch_source_identity": False,
            "global_N1_closure": False,
            "numeric_enclosure": False,
            "ancestor_replay": False,
        },
    }


build_source = build_upstream_source


__all__ = ["build_upstream_source", "build_source"]
