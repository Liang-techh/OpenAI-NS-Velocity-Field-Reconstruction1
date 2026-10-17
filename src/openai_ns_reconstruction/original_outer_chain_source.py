"""Archive-only source adapter for the six outer Rm seams.

The accepted Rm outer join receipt contains source wiring for every chart in
the postrestore-to-power chain.  This module binds that receipt to the
serialized endpoint ports in the same immutable archive.  It does not run a
producer or checker, evaluate a source value, or turn a cap into a function
identity.  The postrestore-to-Rm seam is retained as an already closed
reference seam; ``seams`` contains only the six remaining joins.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT
from .original_Rm_inlet_source_theorem import (
    _JOIN_ARCHIVE,
    _JOIN_CHECK_SOURCE,
    _JOIN_RECEIPT,
    _JOIN_SOURCE,
    _SOURCE_PINS,
    _digest,
    _handle,
    _mapping,
    _node,
    _pin,
    _repo_root,
)


_BINDING_KEY = "Rm_outer_mixed_source_join_bindings"
_SEAMS = (
    "patch_reference",
    "reference_slope",
    "slope_axial",
    "axial_buffer",
    "buffer_transition",
    "transition_power",
)
_ALL_SEAMS = ("postrestore_patch",) + _SEAMS
_CHARTS = (
    "postrestore",
    "Rm_patch",
    "Rh_reference",
    "O2_slope",
    "O2_axial",
    "O2_buffer",
    "O3_transition",
    "O3_power",
)
_ROOT_GERMS = ("E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa")
_FINITE_CATEGORIES = (
    ("loop", "A"),
    ("loop", "B_over_Pstar"),
    ("loop", "A_phi"),
    ("loop", "B_over_Pstar_phi"),
    ("finite_N", "E_N"),
    ("finite_N", "V_N"),
    ("physical", "U"),
)


def _as_mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _load_serialized(value: Any, repo: Path, *, name: str) -> tuple[Mapping[str, Any], dict[str, Any]]:
    if isinstance(value, Mapping):
        return value, {"passed": True, "mode": "supplied_mapping", "name": name}
    if not isinstance(value, (str, Path)):
        return {}, {"passed": False, "name": name, "reason": "serialized mapping was not supplied"}
    path = Path(value)
    if not path.is_absolute():
        path = repo / path
    try:
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8") as stream:
            data = json.load(stream)
    except (OSError, ValueError) as exc:
        return {}, {"passed": False, "name": name, "reason": "serialized mapping is unreadable", "error": str(exc)}
    return data if isinstance(data, Mapping) else {}, {"passed": isinstance(data, Mapping), "name": name, "path": str(path)}


def _node_count(nodes: Any) -> int:
    try:
        return len(nodes)
    except TypeError:
        return 0


def _node_list(nodes: Any) -> list[Any]:
    if isinstance(nodes, list):
        return nodes
    try:
        return list(nodes)
    except TypeError:
        return []


def _prefix_relation(nodes: Any, archive_nodes: Any) -> dict[str, Any]:
    current = _node_list(nodes)
    archived = _node_list(archive_nodes)
    same = bool(current and current is archived)
    current_prefix_archive = bool(
        same or (len(current) <= len(archived) and archived[: len(current)] == current)
    )
    archive_prefix_current = bool(
        same or (len(archived) <= len(current) and current[: len(archived)] == archived)
    )
    return {
        "passed": bool(current_prefix_archive or archive_prefix_current),
        "current_node_count": len(current),
        "archive_node_count": len(archived),
        "current_prefix_archive": current_prefix_archive,
        "archive_prefix_current": archive_prefix_current,
        "same_serialized_object": same,
    }


def _shape(value: Any) -> tuple[int, int] | None:
    if not isinstance(value, list) or len(value) != 5 or not all(isinstance(row, list) for row in value):
        return None
    widths = {len(row) for row in value}
    if len(widths) != 1 or next(iter(widths), 0) not in (5, 6):
        return None
    return (5, next(iter(widths)))


def _chain_owner_pin_rows(
    repo: Path,
    receipt: Mapping[str, Any],
    pins: Mapping[str, Mapping[str, Any]],
    hashes: Mapping[str, Any] | None,
    code_hashes: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Check owner maps, accepting receipt input only in mapless fixtures.

    The outer join is a newly admitted input family.  A parent report may not
    have copied its three source pins yet; the accepted receipt still binds
    those files.  Such rows remain visible in ``missing_registration`` for the
    owner to merge.  When the caller supplies owner maps, every pin must be
    represented and match; a represented owner hash must always match.
    """
    inherited_hashes = _as_mapping(hashes)
    inherited_code = _as_mapping(code_hashes)
    receipt_inputs = _as_mapping(receipt.get("input_hashes"))
    owner_maps_supplied = hashes is not None or code_hashes is not None

    def lookup(mapping: Mapping[str, Any], relative: str, basename: str) -> Any:
        for key in (relative, basename):
            if key not in mapping:
                continue
            value = mapping[key]
            if isinstance(value, Mapping):
                for digest_key in ("sha256", "hash", "digest"):
                    if digest_key in value:
                        return value[digest_key]
            return value
        return None

    rows: dict[str, Any] = {}
    missing: list[dict[str, Any]] = []
    for name, pin in pins.items():
        relative = str(pin.get("path", ""))
        basename = Path(relative).name
        actual = pin.get("sha256")
        is_code = relative.startswith(("src/", "scripts/"))
        represented = lookup(inherited_code if is_code else inherited_hashes, relative, basename)
        receipt_expected = lookup(receipt_inputs, relative, basename)
        receipt_self = relative == _JOIN_RECEIPT
        receipt_match = (receipt_expected == actual) if receipt_expected is not None else receipt_self
        owner_exact = represented is not None and represented == actual
        receipt_fallback = represented is None and receipt_match and not owner_maps_supplied
        if represented is None:
            missing.append({
                "map": "code_hashes" if is_code else "hashes",
                "path": relative,
                "sha256": actual,
                "receipt_fallback_used": receipt_fallback,
            })
        rows[name] = {
            "path": relative,
            "actual_sha256": actual,
            "receipt_input_hash_present": receipt_expected is not None,
            "receipt_input_hash_matches": receipt_match,
            "owner_map_entry_present": represented is not None,
            "owner_map_entry_matches": owner_exact,
            "receipt_input_fallback_authorized": receipt_fallback,
            "passed": bool(receipt_match and (owner_exact or receipt_fallback)),
        }
    return {
        "passed": bool(rows) and all(row["passed"] for row in rows.values()),
        "rows": rows,
        "missing_registration": missing,
        "registration_required": bool(missing),
        "owner_maps_supplied": owner_maps_supplied,
    }


def _port_row(physical: Mapping[str, Any], quantity: str, physical_y_order: int, count: int) -> list[int]:
    matrix = physical.get(quantity)
    if not isinstance(matrix, list) or len(matrix) <= physical_y_order:
        return []
    row = matrix[physical_y_order]
    if not isinstance(row, list) or len(row) < count:
        return []
    result = []
    for item in row[:count]:
        value = _handle(item)
        if value is None:
            return []
        result.append(value)
    return result


def _children(row: Mapping[str, Any]) -> list[int]:
    result: list[int] = []
    for key in ("arguments", "numerator", "denominator", "argument", "expression", "variable"):
        value = row.get(key)
        if isinstance(value, list):
            result.extend(item for item in (_handle(v) for v in value) if item is not None)
        else:
            item = _handle(value)
            if item is not None:
                result.append(item)
    return result


def _depends_on(nodes: Any, start: Any, target: Any, *, limit: int = 256) -> bool:
    root = _handle(start)
    wanted = _handle(target)
    if root is None or wanted is None:
        return False
    stack = [root]
    seen: set[int] = set()
    while stack and len(seen) < limit:
        current = stack.pop()
        if current == wanted:
            return True
        if current in seen:
            continue
        seen.add(current)
        row = _node(nodes, current)
        if row is not None:
            stack.extend(child for child in _children(row) if child not in seen)
    return False


def _chart_scope(nodes: Any, name: str, chart: Mapping[str, Any]) -> dict[str, Any]:
    physical = _as_mapping(chart.get("physical"))
    coordinate = _handle(chart.get("coordinate"))
    radius = _handle(chart.get("radius"))
    jacobian = _handle(chart.get("jacobian"))
    phase = _handle(chart.get("phase"))
    domain = chart.get("domain")
    domain_values = [_handle(value) for value in domain] if isinstance(domain, list) and len(domain) == 2 else []
    phase_row = _node(nodes, phase)
    phase_ok = bool(
        phase_row is not None
        and phase_row.get("operation") == "analytic_unary"
        and phase_row.get("name") == "fractional_part"
        and _handle(phase_row.get("argument")) is not None
    )
    n_grid = chart.get("N_y")
    n_ok = bool(
        isinstance(n_grid, list)
        and len(n_grid) == 3
        and all(isinstance(row, list) and len(row) == 4 for row in n_grid)
        and all(_node(nodes, item) is not None for row in n_grid for item in row)
    )
    u_shape = _shape(physical.get("U"))
    q_shape = _shape(physical.get("Q"))
    ports = {
        "U": _port_row(physical, "U", 0, 3),
        "U_y": _port_row(physical, "U", 1, 2),
        "Q": _port_row(physical, "Q", 0, 3),
        "Q_y": _port_row(physical, "Q", 1, 2),
    }
    ports_ok = all(len(values) == count for values, count in ((ports["U"], 3), (ports["U_y"], 2), (ports["Q"], 3), (ports["Q_y"], 2)))
    port_nodes_exist = bool(ports_ok and all(_node(nodes, item) is not None for values in ports.values() for item in values))
    handles_ok = all(value is not None and _node(nodes, value) is not None for value in (coordinate, radius, jacobian, phase))
    return {
        "chart": name,
        "coordinate": coordinate,
        "radius": radius,
        "jacobian": jacobian,
        "phase": phase,
        "domain": domain_values,
        "phase_fractional_part": phase_ok,
        "N_y_3x4_source_grid": n_ok,
        "physical_U_shape": u_shape,
        "physical_Q_shape": q_shape,
        "ports": ports,
        "handles_exist": handles_ok,
        "ports_exist": port_nodes_exist,
        "passed": bool(
            handles_ok
            and len(domain_values) == 2
            and phase_ok
            and n_ok
            and u_shape == (5, 6)
            and q_shape == (5, 5)
            and port_nodes_exist
        ),
    }


def _join_scope(nodes: Any, charts: Mapping[str, Any], seam: str, join: Mapping[str, Any]) -> dict[str, Any]:
    left_name = join.get("left_chart")
    right_name = join.get("right_chart")
    left = _as_mapping(join.get("left"))
    right = _as_mapping(join.get("right"))
    left_chart = _as_mapping(charts.get(left_name))
    right_chart = _as_mapping(charts.get(right_name))
    left_query = _handle(left.get("native_query"))
    right_query = _handle(right.get("native_query"))
    left_domain = left_chart.get("domain")
    right_domain = right_chart.get("domain")
    left_domain_values = [_handle(value) for value in left_domain] if isinstance(left_domain, list) else []
    right_domain_values = [_handle(value) for value in right_domain] if isinstance(right_domain, list) else []
    left_physical = _as_mapping(left.get("physical"))
    right_physical = _as_mapping(right.get("physical"))
    shapes_ok = bool(
        _shape(left_physical.get("U")) == (5, 6)
        and _shape(right_physical.get("U")) == (5, 6)
        and _shape(left_physical.get("Q")) == (5, 5)
        and _shape(right_physical.get("Q")) == (5, 5)
    )
    join_nodes_exist = bool(
        shapes_ok
        and all(
            _node(nodes, item) is not None
            for physical in (left_physical, right_physical)
            for quantity in ("U", "Q")
            for row in physical.get(quantity, ())
            for item in row
        )
    )
    handles_ok = bool(
        left_name in charts
        and right_name in charts
        and left_query is not None
        and right_query is not None
        and (left_query in left_domain_values)
        and (right_query in right_domain_values)
        and _node(nodes, left_query) is not None
        and _node(nodes, right_query) is not None
        and isinstance(left.get("radius"), list)
        and isinstance(right.get("radius"), list)
    )
    return {
        "seam": seam,
        "left_chart": left_name,
        "right_chart": right_name,
        "left_native_query": left_query,
        "right_native_query": right_query,
        "left_radius": left.get("radius"),
        "right_radius": right.get("radius"),
        "left_scope": {
            "chart": left_name,
            "native_query": left_query,
            "coordinate": _handle(left_chart.get("coordinate")),
            "radius": _handle(left_chart.get("radius")),
            "jacobian": _handle(left_chart.get("jacobian")),
            "phase": _handle(left_chart.get("phase")),
        },
        "right_scope": {
            "chart": right_name,
            "native_query": right_query,
            "coordinate": _handle(right_chart.get("coordinate")),
            "radius": _handle(right_chart.get("radius")),
            "jacobian": _handle(right_chart.get("jacobian")),
            "phase": _handle(right_chart.get("phase")),
        },
        "join_physical_shapes": {
            "left_U": _shape(left_physical.get("U")),
            "right_U": _shape(right_physical.get("U")),
            "left_Q": _shape(left_physical.get("Q")),
            "right_Q": _shape(right_physical.get("Q")),
        },
        "join_physical_nodes_exist": join_nodes_exist,
        "endpoint_handles_and_shapes": bool(handles_ok and join_nodes_exist),
    }


def _receipt_rows(receipt: Mapping[str, Any]) -> dict[str, Any]:
    checks = _as_mapping(receipt.get("checks"))
    required_sections = (
        "source_contract",
        "original_radius_and_phase",
        "actual_flat_controls",
        "actual_background_germs",
        "actual_memory_and_leading_inlets",
        "actual_own_mixed_joins",
    )
    section_rows = {name: _as_mapping(checks.get(name)).get("passed") is True for name in required_sections}
    memory = _as_mapping(checks.get("actual_memory_and_leading_inlets"))
    finite = _as_mapping(memory.get("actual_finite_N_endpoint_source_germs"))
    calculus = _as_mapping(memory.get("actual_U_source_radial_derivative_calculus"))
    own = _as_mapping(checks.get("actual_own_mixed_joins"))
    radius_phase = _as_mapping(checks.get("original_radius_and_phase"))
    contract = _as_mapping(checks.get("source_contract"))
    root_rows = finite.get("actual_root_whole_function_value_pairs")
    finite_rows = finite.get("actual_inverse_loop_and_finite_N_source_pairs")
    uy_rows = calculus.get("actual_U_y0_to4_C5_source_binding_rows")
    q_rows = own.get("actual_unbound_leading_Q_y0_to4_C4_pairs")
    own_rows = own.get("actual_seven_own_mixed_rows")
    root_map = {(row[0], row[1]) for row in root_rows or () if isinstance(row, list) and len(row) >= 2}
    finite_map = {}
    for row in finite_rows or ():
        if isinstance(row, list) and len(row) >= 6:
            finite_map[(row[0], row[1], row[2], row[3])] = row
    uy_map = {}
    for row in uy_rows or ():
        if isinstance(row, list) and len(row) >= 5:
            uy_map[(row[0], row[1], row[2])] = row
    q_map = {}
    for row in q_rows or ():
        if isinstance(row, list) and len(row) >= 5:
            q_map[(row[0], row[1], row[2])] = row
    per_seam = {}
    expected_finite = set(_FINITE_CATEGORIES)
    for seam in _ALL_SEAMS:
        roots_ok = {(seam, name) for name in _ROOT_GERMS} <= root_map
        finite_ok = all((seam, group, name, order) in finite_map for group, name in expected_finite for order in range(6))
        q_ok = all((seam, physical_y, order) in q_map for physical_y in (0, 1) for order in range(5))
        seam_own = [row for row in own_rows or () if isinstance(row, Mapping) and row.get("interface") == seam]
        per_seam[seam] = {
            "root_germ_pairs": roots_ok,
            "finite_N_and_physical_source_pairs": finite_ok,
            "Q_Qy_wrapper_pairs": q_ok,
            "own_mixed_row_count": len(seam_own),
            "own_mixed_rows": len(seam_own) == 246,
            "passed": bool(roots_ok and finite_ok and q_ok and len(seam_own) == 246),
        }
    chart_names = {row[0] for row in uy_rows or () if isinstance(row, list) and row}
    source_flags = {
        "seven_original_chart_endpoints": contract.get("actual_seven_endpoints_bound_to_original_high_source_charts") is True,
        "common_N_Ra_phase_origin_and_periodic_phase": radius_phase.get("global_N_Ra_phase_origin_and_periodic_phase_unchanged") is True,
        "common_log_radius_phase_interfaces": set(_ALL_SEAMS) <= set(radius_phase.get("actual_radius_and_common_log_radius_phase_interfaces", ())),
        "root_and_physical_germs": finite.get("root_C5_Z_and_physical_y_germs_from_actual_background_IVP_and_original_generic_operator") is True,
        "ordinary_calculus_and_physical_J_inverse_chain_rule": finite.get("original_ordinary_calculus_and_physical_J_inverse_chain_rule_then_supply_U_y0_to4") is True,
        "finite_nonvalue_comparison": finite.get("actual_A_B_inverse_integrals_and_velocity_outputs_not_assigned_common_values") is True,
        "Uy_variable_J_calculus": calculus.get("actual_U_y2_to4_are_derivatives_of_actual_preceding_rows_with_variable_J") is True,
        "unique_inverse_composition": calculus.get("accepted_unique_inverse_composition_and_equal_background_germs_then_identify_high_U_germs") is True,
        "own_variable_radius": own.get("actual_nonconstant_native_radius_germs_not_constant_J_substitutes") is True,
        "own_histories_retained": own.get("nonzero_histories_and_one_plus_delta_own_recovery_and_m15_pressure_sign_retained") is True,
    }
    all_passed = bool(
        receipt.get("all_passed") is True
        and all(section_rows.values())
        and all(row["passed"] for row in per_seam.values())
        and set(_CHARTS) <= chart_names
        and all(source_flags.values())
    )
    return {
        "passed": all_passed,
        "receipt_all_passed": receipt.get("all_passed") is True,
        "sections": section_rows,
        "source_flags": source_flags,
        "per_seam": per_seam,
        "chart_source_rows": sorted(chart_names),
        "root_germ_pair_count": len(root_rows) if isinstance(root_rows, list) else 0,
        "finite_source_pair_count": len(finite_rows) if isinstance(finite_rows, list) else 0,
        "Uy_source_row_count": len(uy_rows) if isinstance(uy_rows, list) else 0,
        "Q_Qy_wrapper_pair_count": len(q_rows) if isinstance(q_rows, list) else 0,
        "own_mixed_row_count": len(own_rows) if isinstance(own_rows, list) else 0,
        "maps": {"root": root_map, "finite": finite_map, "Uy": uy_map, "Q": q_map},
    }


def _seam_source_rows(
    nodes: Any,
    seam: str,
    join: Mapping[str, Any],
    charts: Mapping[str, Any],
    chart_rows: Mapping[str, Mapping[str, Any]],
    receipt_rows: Mapping[str, Any],
) -> dict[str, Any]:
    left_name = join.get("left_chart")
    right_name = join.get("right_chart")
    left_chart = _as_mapping(charts.get(left_name))
    right_chart = _as_mapping(charts.get(right_name))
    left_ports = chart_rows.get(left_name, {}).get("ports", {})
    right_ports = chart_rows.get(right_name, {}).get("ports", {})
    finite_map = receipt_rows["maps"]["finite"]
    uy_map = receipt_rows["maps"]["Uy"]
    q_map = receipt_rows["maps"]["Q"]
    seam_scope = _join_scope(nodes, charts, seam, join)
    u_rows = []
    for order in range(3):
        source = finite_map.get((seam, "physical", "U", order))
        left = left_ports.get("U", [None] * 3)[order] if len(left_ports.get("U", [])) > order else None
        right = right_ports.get("U", [None] * 3)[order] if len(right_ports.get("U", [])) > order else None
        passed = bool(source and source[4] == left and source[5] == right)
        u_rows.append({
            "quantity": "U",
            "ordinary_Z_order": order,
            "left_node": left,
            "right_node": right,
            "left_source_output_node": source[4] if source else None,
            "right_source_output_node": source[5] if source else None,
            "source_binding": "actual_inverse_loop_and_finite_N_source_pairs",
            "passed": passed,
        })
    uy_rows = []
    for order in range(2):
        left_source = uy_map.get((left_name, 1, order))
        right_source = uy_map.get((right_name, 1, order))
        left = left_ports.get("U_y", [None] * 2)[order] if len(left_ports.get("U_y", [])) > order else None
        right = right_ports.get("U_y", [None] * 2)[order] if len(right_ports.get("U_y", [])) > order else None
        passed = bool(
            left_source
            and right_source
            and left_source[4] == left
            and right_source[4] == right
        )
        uy_rows.append({
            "quantity": "U_y",
            "ordinary_Z_order": order,
            "left_node": left,
            "right_node": right,
            "left_source_output_node": left_source[4] if left_source else None,
            "right_source_output_node": right_source[4] if right_source else None,
            "left_native_input_node": left_source[3] if left_source else None,
            "right_native_input_node": right_source[3] if right_source else None,
            "source_binding": "actual_U_y0_to4_C5_source_binding_rows",
            "passed": passed,
        })
    join = _as_mapping(join)
    left_join_physical = _as_mapping(_as_mapping(join.get("left")).get("physical"))
    right_join_physical = _as_mapping(_as_mapping(join.get("right")).get("physical"))
    q_rows = []
    qy_rows = []
    for order in range(3):
        source = q_map.get((seam, 0, order))
        left = _port_row(left_join_physical, "Q", 0, 5)
        right = _port_row(right_join_physical, "Q", 0, 5)
        left_node = left[order] if len(left) > order else None
        right_node = right[order] if len(right) > order else None
        passed = bool(source and source[3] == left_node and source[4] == right_node)
        q_rows.append({
            "quantity": "Q",
            "ordinary_Z_order": order,
            "left_node": left_node,
            "right_node": right_node,
            "left_source_output_node": source[3] if source else None,
            "right_source_output_node": source[4] if source else None,
            "left_endpoint_wrapper": source[3] if source else None,
            "right_endpoint_wrapper": source[4] if source else None,
            "source_binding": "actual_unbound_leading_Q_y0_to4_C4_pairs",
            "passed": passed,
        })
    for order in range(2):
        source = q_map.get((seam, 1, order))
        left = _port_row(left_join_physical, "Q", 1, 5)
        right = _port_row(right_join_physical, "Q", 1, 5)
        left_node = left[order] if len(left) > order else None
        right_node = right[order] if len(right) > order else None
        passed = bool(source and source[3] == left_node and source[4] == right_node)
        qy_rows.append({
            "quantity": "Q_y",
            "ordinary_Z_order": order,
            "left_node": left_node,
            "right_node": right_node,
            "left_source_output_node": source[3] if source else None,
            "right_source_output_node": source[4] if source else None,
            "left_endpoint_wrapper": source[3] if source else None,
            "right_endpoint_wrapper": source[4] if source else None,
            "source_binding": "actual_unbound_leading_Q_y0_to4_C4_pairs",
            "passed": passed,
        })
    source_germs = {
        "root_germs": [name for name in _ROOT_GERMS if (seam, name) in receipt_rows["maps"]["root"]],
        "finite_categories": [
            [group, name, order]
            for group, name in _FINITE_CATEGORIES
            for order in range(3)
            if (seam, group, name, order) in finite_map
        ],
        "finite_endpoint_source_rows": [
            finite_map[(seam, group, name, order)]
            for group, name in _FINITE_CATEGORIES
            for order in range(3)
            if (seam, group, name, order) in finite_map
        ],
    }
    return {
        **seam_scope,
        "U_rows": u_rows,
        "Uy_rows": uy_rows,
        "Q_rows": q_rows,
        "Qy_rows": qy_rows,
        "source_germ_pairs": source_germs,
        "passed": bool(
            seam_scope["endpoint_handles_and_shapes"]
            and all(row["passed"] for row in u_rows + uy_rows + q_rows + qy_rows)
            and len(source_germs["root_germs"]) == len(_ROOT_GERMS)
        ),
    }


def build_outer_chain_source(
    nodes: Any,
    archive: Mapping[str, Any] | str | Path,
    receipt: Mapping[str, Any] | str | Path,
    *,
    root: Any = ROOT,
    hashes: Mapping[str, Any] | None = None,
    code_hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Bind all six remaining outer source joins to actual serialized ports."""
    repo = _repo_root(root)
    archive_data, archive_load = _load_serialized(archive, repo, name="join_archive")
    receipt_data, receipt_load = _load_serialized(receipt, repo, name="join_receipt")
    binding = _as_mapping(archive_data.get(_BINDING_KEY))
    charts = _as_mapping(binding.get("charts"))
    joins = _as_mapping(binding.get("joins"))
    pins = {name: _pin(repo, relative, expected) for name, (relative, expected) in _SOURCE_PINS.items()}
    owner_pins = _chain_owner_pin_rows(repo, receipt_data, pins, hashes, code_hashes)
    prefix = _prefix_relation(nodes, archive_data.get("exact_graph_nodes"))
    family_match = bool(
        _as_mapping(archive_data.get("source_family"))
        and _as_mapping(receipt_data.get("source_family"))
        and archive_data.get("source_family") == receipt_data.get("source_family")
    )
    receipt_rows = _receipt_rows(receipt_data)
    chart_rows = {name: _chart_scope(nodes, name, _as_mapping(charts.get(name))) for name in _CHARTS}
    nonconstant_j = {}
    for name in ("Rm_patch", "O2_axial"):
        row = chart_rows[name]
        nonconstant_j[name] = bool(
            row["jacobian"] is not None
            and row["coordinate"] is not None
            and row["jacobian"] != 1
            and _depends_on(nodes, row["jacobian"], row["coordinate"])
        )
    chart_rows["Rm_patch"]["nonconstant_J_source_dependency"] = nonconstant_j["Rm_patch"]
    chart_rows["O2_axial"]["nonconstant_J_source_dependency"] = nonconstant_j["O2_axial"]
    chart_source_passed = bool(
        all(row["passed"] for row in chart_rows.values())
        and all(nonconstant_j.values())
    )
    seam_rows = {}
    for seam in _SEAMS:
        seam_rows[seam] = _seam_source_rows(
            nodes,
            seam,
            _as_mapping(joins.get(seam)),
            charts,
            chart_rows,
            receipt_rows,
        )
    seam_passed = bool(all(row["passed"] for row in seam_rows.values()))
    source_pins_passed = bool(
        archive_load.get("passed") is True
        and receipt_load.get("passed") is True
        and all(row["passed"] for row in pins.values())
        and owner_pins.get("passed") is True
    )
    common_receipt = _digest({
        "accepted_receipt_sha256": pins["accepted_join_receipt"]["sha256"],
        "source_family": archive_data.get("source_family"),
        "seams": _SEAMS,
        "receipt_checks": {
            "passed": receipt_rows["passed"],
            "sections": receipt_rows["sections"],
            "source_flags": receipt_rows["source_flags"],
            "per_seam": receipt_rows["per_seam"],
            "counts": {
                key: receipt_rows[key]
                for key in (
                    "root_germ_pair_count",
                    "finite_source_pair_count",
                    "Uy_source_row_count",
                    "Q_Qy_wrapper_pair_count",
                    "own_mixed_row_count",
                )
            },
        },
    })
    receipt_public = {key: value for key, value in receipt_rows.items() if key != "maps"}
    passed = bool(
        source_pins_passed
        and prefix["passed"]
        and family_match
        and receipt_rows["passed"]
        and chart_source_passed
        and seam_passed
    )
    u_rows = []
    uy_rows = []
    q_rows = []
    qy_rows = []
    for seam in _SEAMS:
        scope = {
            "left": seam_rows[seam].get("left_scope"),
            "right": seam_rows[seam].get("right_scope"),
        }
        for target, rows in ((u_rows, seam_rows[seam]["U_rows"]), (uy_rows, seam_rows[seam]["Uy_rows"]),
                             (q_rows, seam_rows[seam]["Q_rows"]), (qy_rows, seam_rows[seam]["Qy_rows"])):
            for row in rows:
                row["seam"] = seam
                row["common_response_receipt_sha256"] = common_receipt
                row["endpoint_scope"] = scope
                target.append(row)
    return {
        "passed": passed,
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "seams": seam_rows,
        "rows": seam_rows,
        "U_rows": u_rows,
        "Uy_rows": uy_rows,
        "Q_rows": q_rows,
        "Qy_rows": qy_rows,
        "charts": chart_rows,
        "closed_reference_seam": "postrestore_patch",
        "source_family": {
            "archive": archive_data.get("source_family"),
            "receipt": receipt_data.get("source_family"),
            "matches": family_match,
        },
        "source_prefix": prefix,
        "source_pins": {"passed": source_pins_passed, "rows": pins, "owner": owner_pins},
        "archive_load": archive_load,
        "receipt_load": receipt_load,
        "receipt_checks": receipt_public,
        "common_response_receipt_sha256": common_receipt,
        "limits": {
            "source_values_selected": False,
            "numeric_envelopes_reinterpreted": False,
            "opaque_global_alias": False,
            "postrestore_patch_reproved_here": False,
            "ancestor_constructor_or_checker_rerun": False,
        },
        "unproved": {
            "signed_endpoint_values": True,
            "mean_and_pressure_recovery": True,
            "numerical_cap_contraction": True,
        },
    }


__all__ = ["build_outer_chain_source"]
