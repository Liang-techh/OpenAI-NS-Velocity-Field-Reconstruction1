"""Source-only Q4 native derivative wiring for the upstream seams.

This module inspects the serialized loop, modified-pressure, and leading
source graphs.  It proves the graph-level ``slow_B_y``/``V_y``/``U_y``
calculus used at the actual ``second_switch(1)``, ``post_power(0, 1)``, and
``long_reshape(0)`` scopes.  The source seam certificate is supplied by the
Gibbs/source layer; a generic seam ``passed`` flag is never used as a
substitute for its explicit contract, successor-inlet, physical-y, and
inverse/global-phase premises.

No owner constructor, producer, checker, numerical field evaluation, or error
cap is run here.  Endpoint equality remains a separate source theorem.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT
from .original_upstream_native import _radius_guard as _upstream_radius_guard
from .original_native_derivative_guards import (
    _ast_bundle,
    _contains,
    _handle,
    _int_row,
    _mapping,
    _phase_binding,
    _positive_jacobian_premise,
    _radius_log_jacobian_binding,
    _slow_b_source_binding,
)


_LOOP_KEY = "prescribed_loop_C5_bindings"
_MODIFIED_KEY = "modified_Q4_pressure_C3_bindings"
_LEADING_KEY = "prescribed_R100_Rm_C6_bindings"

_CHARTS = ("second_switch", "post_power", "long_reshape")
_ORDERS = (0, 1)
_EXPECTED = {
    "second_switch": {"coordinate": 75949, "radius": 88600, "J": 9, "leading": "second_switch"},
    "post_power": {"coordinate": 78472, "radius": 91107, "J": 91103, "leading": "post_switch"},
    "long_reshape": {"coordinate": 79966, "radius": 79970, "J": 7, "leading": "long_reshape"},
}
_SEAM_EXPECTED = {
    "second_power": {"left": ("second_switch", 1), "right": ("post_power", 0)},
    "power_long": {"left": ("post_power", 1), "right": ("long_reshape", 0)},
}

# These are source-code pins, rather than checker or owner pins.  The AST
# bundle also reads the bridge and leading source files, so pin those files as
# part of the source calculus certificate.
_SOURCE_PINS = {
    "loop": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
        "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    ),
    "modified": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
        "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    ),
    "leading": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py",
        "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    ),
    "bridge": (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py",
        "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    ),
}


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return Path(__file__).resolve().parents[2]


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping) and isinstance(value.get(key), Mapping):
        return value[key]
    return value if isinstance(value, Mapping) else {}


def _row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    return _int_row(nodes, _handle(value))


def _sequence(value: Any, minimum: int = 2) -> bool:
    return isinstance(value, (list, tuple)) and len(value) >= minimum and all(
        _handle(item) is not None for item in value[:minimum]
    )


def _pin(repo: Path, name: str, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    try:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        actual = None
    return {
        "name": name,
        "path": relative,
        "exists": path.is_file(),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
    }


def _prefix(nodes: Any, archive: Mapping[str, Any]) -> dict[str, Any]:
    archived = archive.get("exact_graph_nodes")
    current = nodes if isinstance(nodes, (list, tuple)) else ()
    ok = bool(isinstance(archived, list) and len(current) >= len(archived))
    if ok:
        try:
            ok = list(current[: len(archived)]) == archived
        except (IndexError, TypeError):
            ok = False
    return {
        "passed": ok,
        "archived_node_count": len(archived) if isinstance(archived, list) else 0,
        "current_node_count": len(current),
        "exact_prefix": ok,
    }


def _same_handles(left: Any, right: Any) -> bool:
    return _handle(left) == _handle(right) and _handle(left) is not None


def _zero_row(nodes: Any, value: Any) -> bool:
    row = _row(nodes, value)
    return bool(row and row.get("operation") == "exact_rational"
                and row.get("numerator") == 0 and row.get("denominator") == 1)


def _direct_product(nodes: Any, value: Any, factors: tuple[Any, ...]) -> dict[str, Any]:
    row = _row(nodes, value)
    expected = tuple(_handle(item) for item in factors)
    args = tuple(_handle(item) for item in row.get("arguments", ())) if row else ()
    return {
        "node": _handle(value),
        "operation": row.get("operation") if row else None,
        "arguments": list(args),
        "expected_factors": list(expected),
        "passed": bool(row and row.get("operation") == "product" and len(args) == len(expected)
                        and all(args.count(item) == expected.count(item) for item in set(expected))),
    }


def _common_n(nodes: Any, loop: Mapping[str, Any]) -> dict[str, Any]:
    n = _handle(loop.get("N"))
    row = _row(nodes, n)
    passed = bool(n == 18 and row and row.get("operation") == "exact_positive_integer_expression"
                  and row.get("exact_positive_integer") is True
                  and row.get("Z_independent") is True
                  and row.get("one_common_N_for_all_charts_and_control_denominators") is True)
    return {
        "passed": passed,
        "node": n,
        "row": row,
        "one_common_N_for_all_charts": passed,
        "Z_independent": bool(row and row.get("Z_independent") is True),
    }


def _inverse_n_node(nodes: Any, values: Any, n: Any) -> tuple[int | None, dict[str, Any]]:
    rows = values if isinstance(values, (list, tuple)) else ()
    candidates = []
    for value in rows:
        row = _row(nodes, value)
        args = tuple(_handle(item) for item in row.get("arguments", ())) if row else ()
        if row and row.get("operation") == "sum":
            # The graph keeps ``slow_B/J`` and ``1/N`` in a nested product,
            # then adds the held-phase B_phi row.  Inspect that product's
            # operands rather than looking only at the outer sum.
            for item in args:
                child = _row(nodes, item)
                child_args = tuple(_handle(arg) for arg in child.get("arguments", ())) if child else ()
                if not (child and child.get("operation") == "product"):
                    continue
                for nested in child_args:
                    candidate = _row(nodes, nested)
                    if (candidate and candidate.get("operation") == "positive_quotient"
                            and candidate.get("numerator") == 1
                            and _same_handles(candidate.get("denominator"), n)
                            and candidate.get("source_positive_certificate") == "same original exact positive common N"):
                        candidates.append(nested)
    unique = sorted(set(candidates))
    return (unique[0] if len(unique) == 1 else None,
            {"candidates": unique, "passed": len(unique) == 1})


def _vy_formula(nodes: Any, vy: Any, slow_b: Any, bphi: Any, n: Any) -> dict[str, Any]:
    slow = slow_b if isinstance(slow_b, (list, tuple)) else ()
    phi = bphi if isinstance(bphi, (list, tuple)) else ()
    values = vy if isinstance(vy, (list, tuple)) else ()
    rows = []
    inv_n, inv_n_evidence = _inverse_n_node(nodes, values, n)
    for order in _ORDERS:
        vy_i = _handle(values[order]) if len(values) > order else None
        slow_i = _handle(slow[order]) if len(slow) > order else None
        phi_i = _handle(phi[order]) if len(phi) > order else None
        row = _row(nodes, vy_i)
        args = tuple(_handle(item) for item in row.get("arguments", ())) if row else ()
        product_rows = []
        for item in args:
            child = _row(nodes, item)
            if child and child.get("operation") == "product":
                child_args = tuple(_handle(arg) for arg in child.get("arguments", ()))
                product_rows.append({
                    "node": item,
                    "arguments": list(child_args),
                    "slow_B_node": slow_i,
                    "inv_N_node": inv_n,
                    "passed": bool(len(child_args) == 2 and set(child_args) == {slow_i, inv_n}),
                })
        slow_term = [item for item in product_rows if item["passed"]]
        phi_count = sum(item == phi_i for item in args)
        row_passed = bool(row and row.get("operation") == "sum" and len(args) == 2
                          and len(slow_term) == 1 and phi_count == 1)
        rows.append({
            "ordinary_Z_order": order,
            "V_y_node": vy_i,
            "slow_B_y_node": slow_i,
            "B_phi_node": phi_i,
            "row": row,
            "arguments": list(args),
            "slow_B_over_N_product": slow_term[0] if slow_term else None,
            "B_phi_direct_count": phi_count,
            "passed": row_passed,
        })
    return {
        "passed": bool(rows) and inv_n_evidence["passed"] and all(item["passed"] for item in rows),
        "rows": rows,
        "inv_N_node": inv_n,
        "inv_N_source": inv_n_evidence,
        "formula": "V_y = leading_V_y + (1/N) d_native(B/Pstar)/J + B_over_Pstar_phi",
        "fast_phase_term_added_once": bool(rows) and all(item["B_phi_direct_count"] == 1 for item in rows),
        "slow_B_divided_by_common_N_once": bool(rows) and all(item["slow_B_over_N_product"] is not None for item in rows),
    }


def _physical_uy(nodes: Any, physical: Mapping[str, Any], vy: Any, pstar: Any) -> dict[str, Any]:
    uy = physical.get("U_y")
    values = uy if isinstance(uy, (list, tuple)) else ()
    vvalues = vy if isinstance(vy, (list, tuple)) else ()
    rows = []
    for order in _ORDERS:
        uy_i = _handle(values[order]) if len(values) > order else None
        vy_i = _handle(vvalues[order]) if len(vvalues) > order else None
        row = _direct_product(nodes, uy_i, (vy_i, pstar))
        row["ordinary_Z_order"] = order
        row["V_y_node"] = vy_i
        row["Pstar_node"] = _handle(pstar)
        rows.append(row)
    return {
        "passed": bool(rows) and all(item["passed"] for item in rows),
        "rows": rows,
        "physical_U_y_equals_Pstar_times_V_y": bool(rows) and all(item["passed"] for item in rows),
        "Pstar_once": bool(rows) and all(item["arguments"].count(_handle(pstar)) == 1 for item in rows),
    }


def _physical_u(nodes: Any, physical: Mapping[str, Any], source_fields: Mapping[str, Any], pstar: Any) -> dict[str, Any]:
    values = physical.get("U")
    inputs = source_fields.get("V")
    rows = []
    for order in range(3):
        u_i = _handle(values[order]) if isinstance(values, (list, tuple)) and len(values) > order else None
        v_i = _handle(inputs[order]) if isinstance(inputs, (list, tuple)) and len(inputs) > order else None
        row = _direct_product(nodes, u_i, (v_i, pstar))
        row["ordinary_Z_order"] = order
        rows.append(row)
    return {"passed": bool(rows) and all(item["passed"] for item in rows), "rows": rows,
            "Pstar_once": bool(rows) and all(item["arguments"].count(_handle(pstar)) == 1 for item in rows)}


def _section_candidates(record: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    values = [record]
    for key in ("source_join", "source_proof_premises", "source_certificate", "source_germ",
                "source_receipt_checks", "germ", "certificate"):
        value = record.get(key)
        if isinstance(value, Mapping):
            values.append(value)
    # The Gibbs/source certificate keeps the dedicated contract, successor,
    # and physical-y premises under source_join.source_receipt_checks.  Those
    # sections remain explicit; this is only a bounded one-level lookup.
    for container in tuple(values):
        for key in ("source_receipt_checks", "source_contract", "successor_incoming_equalities",
                    "background_physical_y0_to5", "unique_inverse_global_phase_source_germ"):
            value = container.get(key)
            if isinstance(value, Mapping):
                values.append(value)
    return values


def _explicit_section(record: Mapping[str, Any], names: tuple[str, ...], flags: tuple[str, ...]) -> dict[str, Any]:
    for container in _section_candidates(record):
        for name in names:
            value = container.get(name)
            if isinstance(value, Mapping):
                flag = next((key for key in flags if value.get(key) is True), None)
                if value.get("passed") is True:
                    return {"passed": True, "name": name, "flag": flag, "evidence": value}
            elif value is True:
                flag = next((key for key in flags if container.get(key) is True), None)
                return {"passed": True, "name": name, "flag": flag, "evidence": container}
    return {"passed": False, "names": list(names), "required_flags": list(flags)}


def _scope(record: Mapping[str, Any], side: str) -> Mapping[str, Any]:
    for container in (record, _mapping(record.get("scope")), _mapping(record.get("endpoint_scope"))):
        value = container.get(side + "_scope")
        if isinstance(value, Mapping):
            return value
        value = _mapping(container.get(side))
        if value:
            return value
    return {}


def _seam_record(seam_premises: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    for container_key in ("seams", "joins", "source_seam_proofs", "source_premises"):
        container = _mapping(seam_premises.get(container_key))
        for key in {
            name,
            "second_switch_post_power" if name == "second_power" else "post_power_long",
            "second_switch_to_post_power" if name == "second_power" else "post_power_to_long",
            "power_reshape" if name == "power_long" else name,
        }:
            value = container.get(key)
            if isinstance(value, Mapping):
                return value
    value = seam_premises.get(name)
    return value if isinstance(value, Mapping) else {}


def _seam_guard(seam_premises: Mapping[str, Any], name: str, loop: Mapping[str, Any], nodes: Any) -> dict[str, Any]:
    record = _seam_record(seam_premises, name)
    expected = _SEAM_EXPECTED[name]
    left = _scope(record, "left")
    right = _scope(record, "right")

    def scoped(value: Mapping[str, Any], expected_chart: str, expected_query: int) -> bool:
        chart = value.get("chart", value.get("chart_name"))
        query = value.get("native_query_node", value.get("native_query", value.get("query_node")))
        coordinate = value.get("coordinate_node", value.get("coordinate"))
        expected_coordinate = _handle(_mapping(loop.get("charts")).get(expected_chart, {}).get("coordinate"))
        return bool(chart == expected_chart and _handle(query) == expected_query
                    and _handle(coordinate) == expected_coordinate)

    contract = _explicit_section(
        record,
        ("source_contract", "actual_source_contract", "contract"),
        ("accepted_successor_source_contract", "actual_source_contract_bound",
         "actual_successor_inputs_bound", "actual_four_left_endpoint_functions_and_nonzero_successor_inputs_bound",
         "actual_first_switch1_and_second_switch1_phi_U_and_six_histories_are_successor_inputs"),
    )
    successor = _explicit_section(
        record,
        ("successor_incoming_equalities", "actual_successor_incoming_equalities",
         "successor_incoming", "incoming_successor_equalities"),
        ("actual_successor_incoming_equalities", "successor_incoming_equalities_passed",
         "actual_successor_inputs_equal_previous_upper", "actual_nonzero_successor_inputs_bound",
         "actual_four_left_endpoint_functions_and_nonzero_successor_inputs_bound",
         "actual_first_switch1_and_second_switch1_phi_U_and_six_histories_are_successor_inputs"),
    )
    background = _explicit_section(
        record,
        ("background_physical_y0_to5", "background_physical_germs", "actual_background_physical_germs",
         "flat_ODE", "background"),
        ("actual_background_phi_U_and_six_histories_physical_y0_to5",
         "background_phi_U_and_six_histories_physical_y0_to5",
         "background_phi_U_and_six_histories_y0_to5_with_arbitrary_nonzero_inlets",
         "actual_background_physical_E_V_and_five_history_y0_to5_IVP_rows"),
    )
    inverse_phase = _explicit_section(
        record,
        ("unique_inverse_global_phase_source_germ", "unique_inverse_global_phase_germ",
         "inverse_global_phase_source_germ", "phase_inverse_source_germ", "loop_phase",
         "actual_loop_phase_transfer", "phase"),
        ("actual_original_root_inverse_phase_J_bindings",
         "actual_original_roots_inverse_phase_and_NJ_bindings",
         "actual_unique_inverse_periodic_phase_germ",
         "actual_radius_and_global_logarithmic_phase_seams",
         "actual_radius_global_phase_seams"),
    )
    scope_ok = scoped(left, expected["left"][0], expected["left"][1]) and scoped(
        right, expected["right"][0], expected["right"][1]
    )
    record_flags = {
        "record_passed": record.get("passed") is True,
        "source_values_selected_false": record.get("source_values_selected", False) is False,
        "numeric_caps_used_false": record.get("numeric_caps_used", False) is False,
        "ancestor_replayed_false": record.get("ancestor_replayed", False) is False,
    }
    passed = bool(scope_ok and all(record_flags.values()) and contract["passed"] and successor["passed"]
                  and background["passed"] and inverse_phase["passed"])
    return {
        "passed": passed,
        "join": name,
        "left_scope": dict(left),
        "right_scope": dict(right),
        "expected_scope": {side: {"chart": chart, "native_query": query}
                           for side, (chart, query) in expected.items()},
        "scope_passed": scope_ok,
        "record_flags": record_flags,
        "source_contract": contract,
        "successor_incoming_equalities": successor,
        "background_physical_y0_to5": background,
        "unique_inverse_global_phase_source_germ": inverse_phase,
        "endpoint_field_seam_equality_claimed": False,
    }


def _chart_guard(nodes: Any, name: str, loop: Mapping[str, Any], modified: Mapping[str, Any],
                 leading: Mapping[str, Any], common_n: Mapping[str, Any], ast_rows: Mapping[str, Any]) -> dict[str, Any]:
    expected = _EXPECTED[name]
    loop_charts = _mapping(loop.get("charts")); chart = _mapping(loop_charts.get(name))
    mod_charts = _mapping(modified.get("charts")); mod_chart = _mapping(mod_charts.get(name))
    leading_charts = _mapping(leading.get("charts")); lead_chart = _mapping(leading_charts.get(expected["leading"]))
    loop_fields = _mapping(chart.get("loop")); radial = _mapping(mod_chart.get("radial")); physical = _mapping(mod_chart.get("physical"))
    coordinate = _handle(chart.get("coordinate")); radius = _handle(chart.get("radius")); jacobian = _handle(chart.get("physical_measure"))
    coord_ok = bool(coordinate == expected["coordinate"] and radius == expected["radius"] and jacobian == expected["J"])
    geometry_aliases = {
        "coordinate_matches_leading_phase": _same_handles(chart.get("coordinate"), lead_chart.get("phase")),
        "radius_matches_leading_radius": _same_handles(chart.get("radius"), lead_chart.get("radius")),
        "J_matches_leading_measure": _same_handles(chart.get("physical_measure"), lead_chart.get("physical_measure")),
    }
    phase = _phase_binding(nodes, chart, common_n.get("node"))
    # The shared derivative helper handles affine log-radius rows, while the
    # upstream source has two exact product variants: second_switch nests
    # ``1+s`` and long_reshape stores ``T*q``.  Reuse the source-only upstream
    # guard for those actual graph shapes instead of accepting a failed
    # generic parser.
    radius_guard = _upstream_radius_guard(nodes, name, chart, loop, common_n.get("node"), loop.get("phase_origin"))
    positive = _positive_jacobian_premise(nodes, chart, leading)
    slow_b = _slow_b_source_binding(nodes, radial.get("slow_B_y"), loop_fields.get("B_over_Pstar"), chart.get("physical_measure"))
    slow_b_seq = radial.get("slow_B_y")
    b_seq = loop_fields.get("B_over_Pstar")
    bphi_seq = loop_fields.get("B_over_Pstar_phi")
    vy_formula = _vy_formula(nodes, radial.get("V_y"), slow_b_seq, bphi_seq, common_n.get("node"))
    leading_vy = _mapping(lead_chart.get("radial")).get("V_y")
    leading_vy_rows = []
    for order in _ORDERS:
        value = _handle(leading_vy[order]) if isinstance(leading_vy, (list, tuple)) and len(leading_vy) > order else None
        leading_vy_rows.append({"ordinary_Z_order": order, "node": value, "exact_zero": _zero_row(nodes, value)})
    leading_zero = bool(leading_vy_rows) and all(row["exact_zero"] for row in leading_vy_rows)
    uy = _physical_uy(nodes, physical, radial.get("V_y"), leading.get("Pstar"))
    u = _physical_u(nodes, physical, _mapping(chart.get("modified_fields")), leading.get("Pstar"))
    source_ast = {
        key: _mapping(ast_rows.get(key))
        for key in ("native_init", "native_seed", "native_substitution", "native_derivative",
                    "native_constructor", "radial_formula", "pressure_formula", "loop_phase", "roots")
    }
    ast_ok = all(item.get("passed") is True for item in source_ast.values())
    rows_ok = bool(_sequence(slow_b_seq) and _sequence(b_seq) and _sequence(bphi_seq)
                   and _sequence(radial.get("V_y")) and _sequence(physical.get("U_y")))
    passed = bool(coord_ok and all(geometry_aliases.values()) and phase.get("passed") is True
                  and radius_guard.get("passed") is True and positive.get("passed") is True
                  and slow_b.get("passed") is True and vy_formula.get("passed") is True
                  and leading_zero and uy.get("passed") is True and ast_ok and rows_ok)
    return {
        "passed": passed,
        "chart": name,
        "coordinate": coordinate,
        "radius": radius,
        "native_J": jacobian,
        "coordinate_geometry": {"passed": coord_ok, "expected": expected, "actual": {"coordinate": coordinate, "radius": radius, "J": jacobian}},
        "geometry_aliases": geometry_aliases,
        "phase": phase,
        "radius_log_jacobian": radius_guard,
        "positive_native_J": positive,
        "ordinary_Z_orders": [0, 1],
        "original_B_over_Pstar_rows": list(b_seq[:2]) if isinstance(b_seq, (list, tuple)) else [],
        "original_B_over_Pstar_phi_rows": list(bphi_seq[:2]) if isinstance(bphi_seq, (list, tuple)) else [],
        "slow_B_source_binding": slow_b,
        "leading_V_y_rows": leading_vy_rows,
        "leading_V_y_exact_zero": leading_zero,
        "V_y_wiring": vy_formula,
        "physical_U_y_wiring": uy,
        "physical_U_wiring": u,
        "source_ast": source_ast,
        "fast_phase_fixed_before_B_phi": ast_rows.get("native_constructor", {}).get("passed") is True,
        "native_J_division_once": slow_b.get("passed") is True and radius_guard.get("passed") is True,
        "total_V_y_keeps_slow_and_B_phi": bool(vy_formula.get("passed")),
        "physical_U_y_has_one_Pstar": uy.get("Pstar_once") is True,
        "source_values_selected": False,
        "numeric_caps_used": False,
    }


def build_upstream_derivative(
    nodes: Any,
    loop_archive: Mapping[str, Any] | None,
    modified_archive: Mapping[str, Any] | None,
    leading_archive: Mapping[str, Any] | None,
    seam_premises: Mapping[str, Any] | None,
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Build the actual upstream Q4 slow-B/native derivative certificate."""
    repo = _repo_root(root)
    loop = _unwrap(loop_archive, _LOOP_KEY)
    modified = _unwrap(modified_archive, _MODIFIED_KEY)
    leading = _unwrap(leading_archive, _LEADING_KEY)
    premises = seam_premises if isinstance(seam_premises, Mapping) else {}
    graph = nodes if isinstance(nodes, (list, tuple, Mapping)) else ()

    pins = {name: _pin(repo, name, relative, expected) for name, (relative, expected) in _SOURCE_PINS.items()}
    source_pin_passed = all(item.get("passed") is True for item in pins.values())
    prefixes = {
        "loop": _prefix(graph, loop_archive if isinstance(loop_archive, Mapping) else {}),
        "modified": _prefix(graph, modified_archive if isinstance(modified_archive, Mapping) else {}),
        "leading": _prefix(graph, leading_archive if isinstance(leading_archive, Mapping) else {}),
    }
    graph_prefix_passed = all(item.get("passed") is True for item in prefixes.values())
    ast_rows = _ast_bundle(repo)
    ast_keys = ("native_init", "native_seed", "native_substitution", "native_derivative",
                "native_constructor", "radial_formula", "pressure_formula", "loop_phase", "roots")
    ast_passed = all(_mapping(ast_rows.get(key)).get("passed") is True for key in ast_keys)
    common_n = _common_n(graph, loop)
    loop_charts = _mapping(loop.get("charts"))
    modified_charts = _mapping(modified.get("charts"))
    leading_charts = _mapping(leading.get("charts"))

    charts = {
        name: _chart_guard(graph, name, loop, modified, leading, common_n, ast_rows)
        for name in _CHARTS
    }
    seams = {
        name: _seam_guard(premises, name, loop, graph)
        for name in ("second_power", "power_long")
    }
    chart_passed = all(item.get("passed") is True for item in charts.values())
    seam_passed = all(item.get("passed") is True for item in seams.values())
    ordinary_z1 = bool(charts) and all(item.get("ordinary_Z_orders") == [0, 1] for item in charts.values())
    physical_r_once = all(item.get("physical_U_y_wiring", {}).get("Pstar_once") is True for item in charts.values())
    required = {
        "source_pins": source_pin_passed,
        "graph_prefixes": graph_prefix_passed,
        "source_ast": ast_passed,
        "common_N": common_n.get("passed") is True,
        "chart_derivative_wiring": chart_passed,
        "seam_source_premises": seam_passed,
        "ordinary_Z1": ordinary_z1,
        "physical_U_y_one_Pstar": physical_r_once,
    }
    unproved = [name for name, value in required.items() if not value]
    return {
        "passed": all(required.values()),
        "method": "actual_upstream_Q4_slow_B_native_derivative_source_algebra",
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "producer_or_checker_executed": False,
        "source_pins": pins,
        "source_hashes": pins,
        "graph_prefixes": prefixes,
        "source_ast": {"passed": ast_passed, "rows": ast_rows, "required_rows": list(ast_keys)},
        "common_N": common_n,
        "charts": charts,
        "chart_guards": charts,
        "seam_premises": seams,
        "source_seam_proofs": seams,
        "checks": {
            "required": required,
            "source_pins": {"passed": source_pin_passed, "rows": pins},
            "graph_prefixes": prefixes,
            "source_ast": {"passed": ast_passed, "rows": ast_rows},
            "common_N": common_n,
            "chart_derivative_wiring": charts,
            "seam_source_premises": seams,
        },
        "premises": {
            "native_derivative_holds_fast_phase_fixed": ast_rows.get("native_constructor", {}).get("passed") is True,
            "chart_specific_J_divided_once": all(item.get("native_J_division_once") is True for item in charts.values()),
            "original_B_over_Pstar_definitions_bound": all(bool(item.get("original_B_over_Pstar_rows")) for item in charts.values()),
            "B_phi_added_once": all(item.get("V_y_wiring", {}).get("fast_phase_term_added_once") is True for item in charts.values()),
            "ordinary_Z1_rows": ordinary_z1,
            "physical_U_y_equals_Pstar_times_leading_Vy_plus_slow_B_over_N_plus_B_phi": physical_r_once and chart_passed,
            "physical_R_once": False,
            "common_source_germ_required_per_join": seam_passed,
            "field_seam_equality_proved": False,
        },
        "scope": {
            "charts": list(_CHARTS),
            "mandatory_endpoint_scopes": {
                "second_switch_native1": {"chart": "second_switch", "native_query": 1},
                "post_power_native0": {"chart": "post_power", "native_query": 0},
                "post_power_native1": {"chart": "post_power", "native_query": 1},
                "long_reshape_native0": {"chart": "long_reshape", "native_query": 0},
            },
            "ordinary_Z_orders": [0, 1],
            "common_N": common_n.get("node"),
            "physical_coordinate": "y=log(R/110)",
        },
        "unproved": unproved + ["field_seam_equality", "numeric_endpoint_values", "cap_as_function_value"],
        "limits": {
            "field_seam_equality": False,
            "numeric_endpoint_values": False,
            "numeric_caps_as_values": False,
            "ancestor_constructor_or_producer_replay": False,
        },
    }


__all__ = ["build_upstream_derivative"]
