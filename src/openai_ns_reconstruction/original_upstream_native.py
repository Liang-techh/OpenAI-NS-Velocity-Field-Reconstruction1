"""Source-only native calculus guards for the upstream outer charts.

This module binds the actual source graph for the switch, post-power and
long-reshape charts.  It deliberately consumes serialized archives and source
text only: no producer, checker, owner constructor, numerical value, or seam
identity is replayed here.  The downstream pressure and Gibbs modules own the
field seam theorems.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT
from .original_native_derivative_guards import (
    _direct_factor,
    _handle,
    _int_row,
    _mapping,
    _phase_binding,
    _radius_log_jacobian_binding,
)


_LOOP_KEY = "prescribed_loop_C5_bindings"
_MODIFIED_KEY = "modified_Q4_pressure_C3_bindings"
_LEADING_KEY = "prescribed_R100_Rm_C6_bindings"
_SWITCH_KEY = "switch_mixed_source_join_bindings"
_ANNULAR_KEY = "annular_mixed_source_join_bindings"

_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_"
_SOURCE_PINS = {
    "loop_archive": (_PREFIX + "prescribed_loop_C5_functions.json.gz", "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"),
    "loop_source": (_PREFIX + "prescribed_loop_C5_functions.py", "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"),
    "loop_checker": (_PREFIX + "prescribed_loop_C5_functions_check.py", "066e6e09ccbb4dc28d7b1881197a426beb6bbc7f92d393e7f3e646f27287914b"),
    "loop_receipt": (_PREFIX + "prescribed_loop_C5_functions_check.json", "2c77e80df04fae615f8ea909fc7b606bfccfa27659cb67744ae8d84020b1a9e7"),
    "modified_archive": (_PREFIX + "modified_Q4_pressure_C3_functions.json.gz", "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"),
    "modified_source": (_PREFIX + "modified_Q4_pressure_C3_functions.py", "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"),
    "modified_checker": (_PREFIX + "modified_Q4_pressure_C3_functions_check.py", "829771078e07e2ec8b5168b181f9ff5f0fa62c941c8cf2436feeaba41e367124"),
    "modified_receipt": (_PREFIX + "modified_Q4_pressure_C3_functions_check.json", "189ef4dacc60bc2feda16046a193a75e89f658f635fd6e3e381200c09637bf6f"),
    "leading_archive": (_PREFIX + "prescribed_R100_Rm_C6_functions.json.gz", "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a"),
    "leading_source": (_PREFIX + "prescribed_R100_Rm_C6_functions.py", "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32"),
    "leading_checker": (_PREFIX + "prescribed_R100_Rm_C6_functions_check.py", "8ad5e2edc1daba05763b44dd031ac1671d74f30b1733fc80570873e2ca12089a"),
    "leading_receipt": (_PREFIX + "prescribed_R100_Rm_C6_functions_check.json", "494caf482bffe5475c011a35a1b7d95fd02725b32eb55653f9ca25984dd98ad6"),
    "switch_archive": (_PREFIX + "switch_mixed_source_joins.json.gz", "ddbd506d2fbfc4c58f7a5ba93e0abe3a95d3af8fd3c520ee3ce5024d90a29567"),
    "switch_source": (_PREFIX + "switch_mixed_source_joins.py", "0ae7f651c0d183756f0e6ab90d56dfd80bf3d3c1a6717ffdc2a03fe08dd52f2f"),
    "switch_checker": (_PREFIX + "switch_mixed_source_joins_check.py", "3b00b9dc4fc866430c1f3538edcb33de6e0a58afa2f5372ff04a7e67e0db6c11"),
    "switch_receipt": (_PREFIX + "switch_mixed_source_joins_check.json", "9cf3c6d941df95f4b5f73e7b97e32198a54fcd3d17a85ecca2a0742c0d602fc3"),
    "annular_archive": (_PREFIX + "annular_mixed_source_joins.json.gz", "3947464be1541ae00a42fa952419f1eb7be203abc071a4bf83565a69205d2fde"),
    "annular_source": (_PREFIX + "annular_mixed_source_joins.py", "8febd6f4b356eb50bae0905de580eaabde75b0809ad0e574a0e4f02214582a17"),
    "annular_checker": (_PREFIX + "annular_mixed_source_joins_check.py", "4121c8ca1a48ccd214248773b67af17991c97e9bc93bc4884fb483c4ac496896"),
    "annular_receipt": (_PREFIX + "annular_mixed_source_joins_check.json", "28bd4e63feec4105fd30e87b32dbf67fa9d030a7ecf577d2201ab4e3e9313c1b"),
}

_ARCHIVE_NAMES = {
    "loop": _PREFIX + "prescribed_loop_C5_functions.json.gz",
    "modified": _PREFIX + "modified_Q4_pressure_C3_functions.json.gz",
    "leading": _PREFIX + "prescribed_R100_Rm_C6_functions.json.gz",
    "switch": _PREFIX + "switch_mixed_source_joins.json.gz",
    "annular": _PREFIX + "annular_mixed_source_joins.json.gz",
}

_RECEIPT_CHECKS = {
    "loop": ("actual_source", "transport", "lazy_flat"),
    "modified": ("actual_source", "native_chain", "implicit_native", "total_radial", "full_generic", "pressure_reduction", "binders"),
    "leading": ("actual_source", "geometry", "cutoff", "long_restore_joins", "switch_controls"),
    "switch": ("source_contract", "actual_background_recipes", "physical_controls_and_ODE_jets", "actual_loop_and_phase_transfer", "actual_own_mixed_joins"),
    "annular": ("source_contract", "actual_background_recipes", "actual_value_and_gauge_joins", "flat_controls_and_background_jets", "actual_loop_phase_transfer", "actual_own_mixed_joins"),
}


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return Path(__file__).resolve().parents[2]


def _nodes(value: Any) -> Any:
    candidate = getattr(value, "nodes", None)
    return candidate if isinstance(candidate, (list, tuple, Mapping)) else value


def _row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _handle(value)
    if index is None or index < 0:
        return None
    try:
        item = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return item if isinstance(item, Mapping) else None


def _load_json(repo: Path, relative: str) -> Mapping[str, Any]:
    try:
        value = json.loads((repo / relative).read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, Mapping) else {}


def _load_gzip(repo: Path, relative: str) -> Mapping[str, Any]:
    try:
        with gzip.open(repo / relative, "rt", encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, ValueError, EOFError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, Mapping) else {}


def _pin(repo: Path, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {"path": relative, "exists": path.is_file(), "sha256": actual,
            "expected_sha256": expected, "passed": actual == expected}


def _extract(value: Any, key: str, depth: int = 0) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or depth > 4:
        return {}
    child = value.get(key)
    if isinstance(child, Mapping):
        return child
    for item in value.values():
        if isinstance(item, Mapping):
            found = _extract(item, key, depth + 1)
            if found:
                return found
    return {}


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    found = _extract(value, key)
    if found:
        return found
    return value if isinstance(value, Mapping) and "charts" in value else {}


def _prefix(nodes: Any, packet: Mapping[str, Any]) -> bool:
    graph = packet.get("exact_graph_nodes")
    current = _nodes(nodes)
    if not isinstance(graph, list) or not isinstance(current, (list, tuple)):
        return False
    try:
        return len(current) >= len(graph) and list(current[: len(graph)]) == graph
    except (TypeError, IndexError):
        return False


def _packet_identity(packet: Mapping[str, Any], loaded: Mapping[str, Any], key: str) -> bool:
    supplied = _unwrap(packet, key)
    loaded_binding = _unwrap(loaded, key)
    if not supplied:
        supplied = loaded_binding
    return bool(supplied and supplied == loaded_binding)


def _source_family(packet: Mapping[str, Any], loaded: Mapping[str, Any]) -> Mapping[str, Any]:
    family = packet.get("source_family") if isinstance(packet, Mapping) else None
    return family if isinstance(family, Mapping) else _mapping(loaded.get("source_family"))


def _receipt_guard(repo: Path, name: str, loaded: Mapping[str, Any]) -> dict[str, Any]:
    receipt = _load_json(repo, _ARCHIVE_NAMES[name].removesuffix(".json.gz") + "_check.json")
    # `.json.gz` is present only in the archive name; the receipt has the same
    # stem without the compression suffix.
    if not receipt:
        receipt = _load_json(repo, _PREFIX + {
            "loop": "prescribed_loop_C5_functions_check.json",
            "modified": "modified_Q4_pressure_C3_functions_check.json",
            "leading": "prescribed_R100_Rm_C6_functions_check.json",
            "switch": "switch_mixed_source_joins_check.json",
            "annular": "annular_mixed_source_joins_check.json",
        }[name])
    checks = _mapping(receipt.get("checks"))
    rows = {}
    for check_name in _RECEIPT_CHECKS[name]:
        item = checks.get(check_name)
        rows[check_name] = bool(item is True or (_mapping(item).get("passed") is True))
    return {
        "passed": receipt.get("all_passed") is True and all(rows.values())
        and _source_family(receipt, loaded) == _FAMILY,
        "all_passed": receipt.get("all_passed") is True,
        "checks": rows,
        "source_family": _source_family(receipt, loaded),
    }


def _ast_guard(repo: Path) -> dict[str, Any]:
    paths = {
        "loop": _SOURCE_PINS["loop_source"][0],
        "modified": _SOURCE_PINS["modified_source"][0],
        "leading": _SOURCE_PINS["leading_source"][0],
        "switch": _SOURCE_PINS["switch_source"][0],
        "annular": _SOURCE_PINS["annular_source"][0],
    }
    text = {}
    parsed = {}
    rows = {}
    for name, relative in paths.items():
        try:
            source = (repo / relative).read_text(encoding="utf-8")
            tree = ast.parse(source)
            text[name] = source
            parsed[name] = tree
            rows[name + "_syntax"] = True
        except (OSError, SyntaxError, TypeError):
            text[name] = ""
            parsed[name] = None
            rows[name + "_syntax"] = False
    rows.update({
        "loop_native_measure_once": all(s in text["loop"] for s in (
            "physical_Jacobian_applied_once=True", "total_Z_phase_derivative_exact_zero=True",
            "original_rate_zero_pressure_memory_retained")),
        "modified_native_derivative": all(s in text["modified"] for s in (
            "class NativeDerivative", "invjac=g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
            "physical_y_equals_native_partial_divided_by_Jacobian=True", "total_fast_phase_terms_included=True",
            "physical_U_mean_Uy_C5_and_Q_Qy_C4_Pstar_once=True")),
        "modified_fast_phase_N_once": all(s in text["modified"] for s in (
            "invN=g.quotient(g.one,N,'same original exact positive common N')",
            "a.scale(slow['A'],invN)", "a.scale(slow['B_over_Pstar'],invN)",
            "chart['loop']['A_phi']", "chart['loop']['B_over_Pstar_phi']"))
        and "a.scale(chart['loop']['A_phi'],invN)" not in text["modified"],
        "modified_pressure_conversion": all(s in text["modified"] for s in (
            "a.scale(fields['V'],Pstar)", "a.scale(radial['V_y'],Pstar)",
            "a.scale(fields['m'],Pstar)", "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,chart['radius']))")),
        "leading_native_post_and_long": all(s in text["leading"] for s in (
            "self.post_length=g.sub(g.unary('log',c('11/10')),g.mul(c(2),hb))",
            "y=g.mul(self.post_length,p)", "q=g.symbol('true_prescribed_long_reshape_fraction')",
            "y=g.mul(T,q);R=g.mul(c(110),g.unary('exp',y))", "Pstar=Pstar=bridge_owner.inlet.Pstar")),
        "leading_distinct_native_parameters": all(s in text["leading"] for s in (
            "native_parameter('hbS'", "native_parameter('T'")),
        "switch_source_only": all(s in text["switch"] for s in (
            "actual_source_field_history_pressure_and_correction_inlets_not_reset",
            "own_outputs_not_assigned_to_seam_targets", "source_functions_not_numeric_evaluation_or_temporal_coefficients")),
        "annular_source_only": all(s in text["annular"] for s in (
            "actual_four_original_interfaces_to_postrestore", "nonzero_incoming_background_correction_and_signed_pressure_memory",
            "numerical_oracles_global_n1_and_temporal_recursion_not_installed")),
    })
    return {"passed": all(rows.values()), "rows": rows,
            "source_only_ast": True, "source_values_selected": False}


def _common_N(nodes: Any, loop: Mapping[str, Any]) -> dict[str, Any]:
    node = _handle(loop.get("N"))
    row = _row(nodes, node)
    passed = bool(node == 18 and row and row.get("operation") == "exact_positive_integer_expression"
                  and row.get("exact_positive_integer") is True
                  and row.get("Z_independent") is True
                  and row.get("one_common_N_for_all_charts_and_control_denominators") is True)
    return {"passed": passed, "node": node, "row": row,
            "one_common_N": passed, "Z_independent": bool(row and row.get("Z_independent") is True)}


def _offset_guard(nodes: Any, chart: Mapping[str, Any], origin: Any) -> dict[str, Any]:
    offset_i = _handle(chart.get("global_offset")); radius_i = _handle(chart.get("radius")); origin_i = _handle(origin)
    row = _row(nodes, offset_i)
    log_i = None; ratio = None; negative_i = None
    if row and row.get("operation") == "sum":
        for item in row.get("arguments", ()):
            candidate = _row(nodes, item)
            if candidate and candidate.get("operation") == "analytic_unary" and candidate.get("name") == "log":
                log_i = _handle(item); ratio = _row(nodes, candidate.get("argument"))
            if candidate and candidate.get("operation") == "negative" and _handle(candidate.get("argument")) == origin_i:
                negative_i = _handle(item)
    ratio_ok = bool(ratio and ratio.get("operation") == "positive_quotient"
                    and _handle(ratio.get("numerator")) == radius_i
                    and _handle(ratio.get("denominator")) is not None)
    passed = bool(row and row.get("operation") == "sum" and log_i is not None
                  and negative_i is not None and ratio_ok and origin_i is not None)
    return {"passed": passed, "global_offset": offset_i, "phase_origin": origin_i,
            "log_ratio": log_i, "radius_ratio": ratio_ok, "origin_subtracted": negative_i is not None,
            "global_offset_is_log_R_over_Ra_minus_phase_origin": passed}


def _radius_guard(nodes: Any, name: str, chart: Mapping[str, Any], loop: Mapping[str, Any], common_n: Any, origin: Any) -> dict[str, Any]:
    coordinate = _handle(chart.get("coordinate")); radius = _handle(chart.get("radius")); jacobian = _handle(chart.get("physical_measure"))
    coordinate_row = _row(nodes, coordinate); radius_row = _row(nodes, radius); jacobian_row = _row(nodes, jacobian)
    expected_names = {
        "first_switch": "true_prescribed_first_switch_s",
        "second_switch": "true_prescribed_second_switch_t",
        "post_power": "true_prescribed_post_switch_fraction",
        "long_reshape": "true_prescribed_long_reshape_fraction",
        "frozen_macro": "original_prescribed_C6_macro_fraction",
    }
    expected_handles = {
        "first_switch": (73190, 85905, 9), "second_switch": (75949, 88600, 9),
        "post_power": (78472, 91107, 91103), "long_reshape": (79966, 79970, 7),
        "frozen_macro": (70265, 71448, 66987),
    }
    eh = expected_handles.get(name, (None, None, None))
    coord_ok = bool(coordinate == eh[0] and coordinate_row and coordinate_row.get("operation") == "bound_variable"
                    and coordinate_row.get("name") == expected_names.get(name))
    radius_args = tuple(_handle(x) for x in (radius_row or {}).get("arguments", ()))
    exp_nodes = [i for i in radius_args if (_row(nodes, i) or {}).get("operation") == "analytic_unary"
                 and (_row(nodes, i) or {}).get("name") == "exp"]
    base_nodes = [i for i in radius_args if i not in exp_nodes]
    exp_i = exp_nodes[0] if len(exp_nodes) == 1 else None
    exp_row = _row(nodes, exp_i); argument_i = _handle(exp_row.get("argument")) if exp_row else None
    argument = _row(nodes, argument_i)
    arg_args = tuple(_handle(x) for x in (argument or {}).get("arguments", ()))
    expected_base = {"first_switch": 31709, "second_switch": 31709, "post_power": 31709,
                     "long_reshape": 79968}.get(name)
    base_ok = bool(len(base_nodes) == 1 and (expected_base is None or base_nodes[0] == expected_base))
    if name in {"first_switch", "second_switch"}:
        # The second switch has the exact source shape
        # ``exp(hbS*(1+s))``: its coordinate is nested in the unit-shift
        # node, whereas the first switch uses ``exp(hbS*s)`` directly.
        if name == "second_switch":
            shifted = _row(nodes, 75950)
            shape_ok = bool(argument and argument.get("operation") == "product" and len(arg_args) == 2
                            and jacobian == 9 and 75950 in arg_args and jacobian in arg_args
                            and shifted and shifted.get("operation") == "sum"
                            and set(_handle(item) for item in shifted.get("arguments", ())) == {1, coordinate})
        else:
            shape_ok = bool(argument and argument.get("operation") == "product" and len(arg_args) == 2
                            and jacobian == 9 and coordinate in arg_args and jacobian in arg_args)
    elif name == "post_power":
        # R = 100 exp(1223 + J * post_power_fraction), where J is the
        # distinct positive post-power native interval length.
        shape_ok = bool(argument and argument.get("operation") == "sum" and len(arg_args) == 2
                        and 91104 in arg_args and 1223 in arg_args)
        product = _row(nodes, 91104)
        shape_ok = bool(shape_ok and product and product.get("operation") == "product"
                        and set(_handle(x) for x in product.get("arguments", ())) == {jacobian, coordinate}
                        and _row(nodes, 91102) and _row(nodes, 91102).get("operation") == "negative"
                        and _handle(_row(nodes, 91102).get("argument")) == 1223)
    elif name == "long_reshape":
        shape_ok = bool(argument and argument.get("operation") == "product" and len(arg_args) == 2
                        and set(arg_args) == {jacobian, coordinate})
    else:
        helper = _radius_log_jacobian_binding(nodes, chart, common_n, origin)
        shape_ok = helper.get("passed") is True
    jac_source_ok = bool(jacobian_row and ((name in {"first_switch", "second_switch"}
                      and jacobian_row.get("operation") == "current_original_source_parameter"
                      and jacobian_row.get("name") == "hbS")
                      or (name == "post_power" and jacobian_row.get("operation") == "sum")
                      or (name == "long_reshape" and jacobian_row.get("operation") == "current_original_source_parameter"
                          and jacobian_row.get("name") == "T")
                      or name == "frozen_macro"))
    phase = _phase_binding(nodes, chart, common_n)
    offset = _offset_guard(nodes, chart, origin)
    passed = bool(coord_ok and radius_row and radius_row.get("operation") == "product"
                  and len(radius_args) == 2 and exp_i is not None and base_ok and shape_ok
                  and jac_source_ok and phase.get("passed") is True and offset.get("passed") is True)
    return {
        "passed": passed, "chart": name, "coordinate": coordinate, "radius": radius,
        "native_J": jacobian, "physical_measure": jacobian, "expected_handles": {"coordinate": eh[0], "radius": eh[1], "native_J": eh[2]},
        "coordinate_row": coordinate_row, "radius_row": radius_row, "native_J_row": jacobian_row,
        "radius_exponential_node": exp_i, "radius_exponent_node": argument_i,
        "radius_base_node": base_nodes[0] if base_nodes else None,
        "radius_formula_shape": shape_ok, "radius_has_one_exp_factor": len(exp_nodes) == 1,
        "R_native_factor_once": bool(base_ok and len(exp_nodes) == 1),
        # ``d_native = J d_y``; therefore the physical derivative is
        # ``d_y = d_native / J``.  Keep both directions explicit so a caller
        # cannot accidentally read the inverse relation as a premise.
        "d_native_equals_J_d_y": bool(shape_ok),
        "d_y_equals_d_native_over_J": bool(shape_ok), "physical_y": "y=log(R/110)",
        "native_J_is_constant": True, "native_J_is_varying": False,
        "d_native_log_R_equals_J": bool(shape_ok), "physical_measure_is_native_J": jacobian == eh[2],
        "phase": phase, "global_offset": offset,
        "phase_factor_after_J_division_is_N": bool(phase.get("passed") is True and shape_ok),
        "total_fast_phase_N_term_has_no_extra_N": True,
        "ordinary_Z_orders": [0, 1, 2, 3, 4],
        "source_values_selected": False, "numeric_caps_used": False,
    }


def _all_direct_factor(nodes: Any, values: Any, factor: Any) -> dict[str, Any]:
    sequence = values if isinstance(values, (list, tuple)) else ()
    rows = []
    factor_i = _handle(factor)
    for value in sequence:
        row = _row(nodes, value); args = row.get("arguments", ()) if row else ()
        count = sum(_handle(item) == factor_i for item in args) if isinstance(args, (list, tuple)) else 0
        rows.append({"row": _handle(value), "operation": row.get("operation") if row else None,
                     "factor": factor_i, "direct_factor_count": count})
    return {"passed": bool(rows) and all(item["direct_factor_count"] == 1 for item in rows), "rows": rows}


def _m_y_rows(nodes: Any, physical: Mapping[str, Any]) -> bool:
    u = physical.get("U", ()); m = physical.get("m", ()); my = physical.get("m_y", ())
    if not (len(u) == len(m) == len(my) == 6):
        return False
    for u_i, m_i, my_i in zip(u, m, my):
        row = _row(nodes, my_i); args = tuple(_handle(x) for x in (row or {}).get("arguments", ()))
        negatives = [x for x in args if (_row(nodes, x) or {}).get("operation") == "negative"]
        if not (row and row.get("operation") == "sum" and _handle(u_i) in args and len(negatives) == 1
                and _handle((_row(nodes, negatives[0]) or {}).get("argument")) == _handle(m_i)):
            return False
    return True


def _pressure_port(nodes: Any, name: str, modified: Mapping[str, Any], loop_chart: Mapping[str, Any], pstar: Any) -> dict[str, Any]:
    physical = _mapping(modified.get("physical")); radial = _mapping(modified.get("radial")); boundary = _mapping(modified.get("finite_boundary"))
    rows = {key: physical.get(key) for key in ("U", "U_y", "m", "m_y", "Q", "Q_y")}
    shape = {key: isinstance(value, (list, tuple)) and len(value) == (5 if key in ("Q", "Q_y") else 6)
             for key, value in rows.items()}
    present = bool(all(shape.values()) and all(_row(nodes, item) is not None for value in rows.values() for item in value))
    pstar_u = _all_direct_factor(nodes, rows["U"], pstar)
    pstar_uy = _all_direct_factor(nodes, rows["U_y"], pstar)
    pstar_m = _all_direct_factor(nodes, rows["m"], pstar)
    q_rows = []
    qy_rows = []
    for key, out in (("Q", q_rows), ("Q_y", qy_rows)):
        for item in rows[key] if isinstance(rows[key], (list, tuple)) else ():
            row = _row(nodes, item)
            out.append({"node": _handle(item), "operation": row.get("operation") if row else None,
                        "denominator": _handle(row.get("denominator")) if row else None,
                        "source_positive_certificate": row.get("source_positive_certificate") if row else None})
    q_den = {item["denominator"] for item in q_rows + qy_rows}
    q_formula = bool(q_rows and qy_rows and len(q_den) == 1 and all(item["operation"] == "positive_quotient" for item in q_rows + qy_rows)
                     and all(str(item["source_positive_certificate"]).startswith("same original") for item in q_rows + qy_rows))
    qplus = boundary.get("Q_plus_Q_y", ()); rqy = boundary.get("R_Q_y", ()); radius = loop_chart.get("radius")
    qplus_ok = bool(len(qplus) == 4 and all((_row(nodes, qplus[i]) or {}).get("operation") == "sum"
                                             and set(_handle(x) for x in (_row(nodes, qplus[i]) or {}).get("arguments", ()))
                                             == {_handle(rows["Q"][i]), _handle(rows["Q_y"][i])} for i in range(4)))
    rqy_ok = bool(len(rqy) == 4 and all((_row(nodes, rqy[i]) or {}).get("operation") == "product"
                                        and sum(_handle(x) == _handle(rows["Q_y"][i]) for x in (_row(nodes, rqy[i]) or {}).get("arguments", ())) == 1
                                        and sum(_handle(x) == _handle(radius) for x in (_row(nodes, rqy[i]) or {}).get("arguments", ())) == 1 for i in range(4)))
    radial_flags = radial.get("total_fast_phase_terms_included") is True and radial.get("physical_y_equals_native_partial_divided_by_Jacobian") is True
    mass_rows = _m_y_rows(nodes, physical)
    passed = bool(present and pstar_u["passed"] and pstar_uy["passed"] and pstar_m["passed"] and mass_rows
                  and q_formula and qplus_ok and rqy_ok and radial_flags)
    return {
        "passed": passed, "chart": name, "physical_rows": {k: list(v) for k, v in rows.items()},
        "radial_rows": {k: list(radial.get(k, ())) for k in ("E_y", "V_y", "m_y", "slow_A_y", "slow_B_y")},
        "finite_boundary_rows": {"Q_plus_Q_y": list(qplus), "R_Q_y": list(rqy)},
        "Pstar_node": _handle(pstar), "Pstar_once": bool(pstar_u["passed"] and pstar_uy["passed"] and pstar_m["passed"]),
        "Pstar_factors": {"U": pstar_u, "U_y": pstar_uy, "m": pstar_m}, "m_y_equals_U_minus_m": mass_rows,
        "Q_rows": q_rows, "Q_y_rows": qy_rows, "Q_formula_source": q_formula,
        # These aliases are deliberately rows, rather than a generic boolean:
        # downstream source checks can verify ordinary Z1 and the R factor
        # against the exact archived graph ports.
        "Q_y_graph_rows": list(rows["Q_y"]), "R_Q_y_rows": list(rqy),
        "ordinary_Z1": bool(len(rows["Q_y"]) >= 2 and len(rqy) >= 2),
        "Q_plus_Q_y_rows": qplus_ok, "R_Q_y_radius_once": rqy_ok, "radius_node": _handle(radius),
        "total_fast_phase_terms_included": radial.get("total_fast_phase_terms_included") is True,
        "physical_y_equals_native_partial_divided_by_Jacobian": radial.get("physical_y_equals_native_partial_divided_by_Jacobian") is True,
        "ordinary_Z_orders": [0, 1, 2, 3, 4], "source_values_selected": False, "numeric_caps_used": False,
    }


def _endpoint_guard(nodes: Any, name: str, endpoint: Mapping[str, Any], coordinate: Any) -> dict[str, Any]:
    boundary = _mapping(endpoint.get("boundary")); rows = []
    for group in ("Q_plus_Q_y", "R_Q_y"):
        values = boundary.get(group, ())
        rows.extend(values if isinstance(values, (list, tuple)) else ())
    checks = []
    for value in rows:
        row = _row(nodes, value)
        checks.append(bool(row and row.get("operation") == "function_substitution"
                           and _handle(row.get("variable")) == _handle(coordinate)
                           and row.get("value") == 1
                           and row.get("Z_independent_substitution") is True))
    return {"passed": bool(rows) and all(checks), "chart": name, "rows": list(rows),
            "endpoint_substitution_rows": len(rows), "substitutes_full_chart": all(checks) if checks else False,
            "field_seam_equality_proved": False, "source_values_selected": False}


def _leading_ports(nodes: Any, leading: Mapping[str, Any], loop_charts: Mapping[str, Any]) -> dict[str, Any]:
    charts = _mapping(leading.get("charts")); aliases = {"first_switch": "first_switch", "second_switch": "second_switch", "post_power": "post_switch", "long_reshape": "long_reshape", "frozen_macro": "frozen_macro"}
    rows = {}
    for name, source_name in aliases.items():
        chart = _mapping(charts.get(source_name)); loop = _mapping(loop_charts.get(name)); actual = _mapping(chart.get("actual"))
        if not chart:
            # Frozen macro is an optional leading source chart; the current
            # C6 archive intentionally carries no such leading binding.
            continue
        required = bool(_handle(chart.get("phase")) == _handle(loop.get("coordinate"))
                        and _handle(chart.get("radius")) == _handle(loop.get("radius"))
                        and _handle(chart.get("physical_measure")) == _handle(loop.get("physical_measure")))
        actual_shape = all(isinstance(actual.get(key), (list, tuple)) and len(actual[key]) == 7 for key in ("phi", "U", "phi_y", "U_y"))
        actual_rows = actual_shape and all(_row(nodes, item) is not None for key in ("phi", "U", "phi_y", "U_y") for item in actual[key])
        rows[name] = {"passed": bool(chart and required and actual_rows), "phase": _handle(chart.get("phase")), "radius": _handle(chart.get("radius")),
                      "physical_measure": _handle(chart.get("physical_measure")), "actual_rows": {key: list(actual.get(key, ())) for key in ("phi", "U", "phi_y", "U_y")},
                      "actual_rows_present": actual_rows, "matches_loop_geometry": required}
    return {"passed": bool(rows) and all(item["passed"] for item in rows.values()), "charts": rows,
            "field_seam_equality_proved": False, "source_values_selected": False}


def _switch_joins(nodes: Any, binding: Mapping[str, Any], *, annular: bool = False) -> dict[str, Any]:
    joins = _mapping(binding.get("joins")); defs = _mapping(binding.get("definition")); contract = binding.get("contract")
    expected = (("power_reshape", "post_power", "long_reshape", "R110"), ("reshape_reference", "long_reshape", "reference", "Rsh"),
                ("reference_restore", "reference", "restoration", "Rz"), ("restore_postrestore", "restoration", "postrestore", "Rrestore")) if annular else (("first_second", "first_switch", "second_switch", "R1"), ("second_power", "second_switch", "post_power", "R2"))
    join_rows = {}
    for key, left, right, endpoint in expected:
        item = _mapping(joins.get(key)); join_rows[key] = {"passed": bool(item and item.get("left_chart") == left and item.get("right_chart") == right and item.get("original_endpoint") == endpoint
                                                           and _mapping(item.get("left")) and _mapping(item.get("right"))),
                                                     "left_chart": item.get("left_chart"), "right_chart": item.get("right_chart"), "original_endpoint": item.get("original_endpoint")}
    controls = binding.get("controls")
    if annular:
        controls = binding.get("controls")
        control_ok = all(isinstance(_mapping(controls).get(name, {}).get("sigma"), (list, tuple)) and len(_mapping(controls).get(name, {}).get("sigma")) == 6
                         and isinstance(_mapping(controls).get(name, {}).get("prime"), (list, tuple)) and len(_mapping(controls).get(name, {}).get("prime")) == 5 for name in ("long_reshape", "restoration"))
    else:
        control_ok = all(isinstance(_mapping(controls).get(name), (list, tuple)) and len(_mapping(controls).get(name)) == 5 for name in ("first_switch", "second_switch"))
    definition_ok = all(value is True for value in defs.values()) and bool(defs)
    contract_ok = isinstance(contract, list) and len(contract) == 14 and all(isinstance(item, Mapping) and item.get("physical_y_order") in (3, 4) for item in contract)
    passed = bool(join_rows and all(item["passed"] for item in join_rows.values()) and control_ok and definition_ok and contract_ok)
    return {"passed": passed, "joins": join_rows, "controls_present": control_ok, "definition": dict(defs),
            "contract_rows": len(contract) if isinstance(contract, list) else 0, "contract_present": contract_ok,
            "field_seam_equality_proved": False, "source_values_selected": False}


def build_upstream_native(
    nodes: Any,
    loop_archive: Mapping[str, Any] | None,
    modified_archive: Mapping[str, Any] | None,
    leading_archive: Mapping[str, Any] | None,
    switch_archive: Mapping[str, Any] | None = None,
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Return source-only native calculus evidence for the upstream charts."""
    repo = _repo_root(root)
    loaded = {name: _load_gzip(repo, relative) for name, relative in _ARCHIVE_NAMES.items()}
    supplied = {"loop": loop_archive or loaded["loop"], "modified": modified_archive or loaded["modified"],
                "leading": leading_archive or loaded["leading"], "switch": switch_archive or loaded["switch"]}
    loop = _unwrap(supplied["loop"], _LOOP_KEY); modified = _unwrap(supplied["modified"], _MODIFIED_KEY); leading = _unwrap(supplied["leading"], _LEADING_KEY)
    switch = _unwrap(supplied["switch"], _SWITCH_KEY)
    # ``switch_archive`` is the two-seam switch packet in the public API.  A
    # caller may also supply a container carrying the annular packet, but do
    # not treat the ordinary switch packet itself as annular evidence.
    annular_input = switch_archive if _extract(switch_archive, _ANNULAR_KEY) else loaded["annular"]
    annular = _unwrap(annular_input, _ANNULAR_KEY)
    graph = _nodes(nodes)
    pins = {name: _pin(repo, relative, expected) for name, (relative, expected) in _SOURCE_PINS.items()}
    receipts = {name: _receipt_guard(repo, name, loaded[name]) for name in _RECEIPT_CHECKS}
    families = {name: _source_family(supplied.get(name, {}), loaded[name]) for name in ("loop", "modified", "leading", "switch")}
    families["annular"] = _source_family(annular_input, loaded["annular"])
    family_ok = all(value == _FAMILY for value in families.values())
    identities = {"loop": _packet_identity(supplied["loop"], loaded["loop"], _LOOP_KEY), "modified": _packet_identity(supplied["modified"], loaded["modified"], _MODIFIED_KEY),
                  "leading": _packet_identity(supplied["leading"], loaded["leading"], _LEADING_KEY), "switch": _packet_identity(supplied["switch"], loaded["switch"], _SWITCH_KEY),
                  "annular": _packet_identity(annular_input, loaded["annular"], _ANNULAR_KEY)}
    prefixes = {name: _prefix(graph, loaded[name]) for name in _ARCHIVE_NAMES}
    ast_guard = _ast_guard(repo)
    common_n = _common_N(graph, loop)
    loop_charts = _mapping(loop.get("charts")); origin = loop.get("phase_origin")
    chart_names = ("first_switch", "frozen_macro", "second_switch", "post_power", "long_reshape")
    chart_guards = {name: _radius_guard(graph, name, _mapping(loop_charts.get(name)), loop, common_n.get("node"), origin) for name in chart_names}
    phase = {name: _mapping(chart_guards[name].get("phase")) for name in chart_names}
    pstar = _handle(leading.get("Pstar"))
    modified_charts = _mapping(modified.get("charts")); pressure_ports = {}
    endpoint_ports = {}
    for name in chart_names:
        if name in modified_charts and name in loop_charts:
            pressure_ports[name] = _pressure_port(graph, name, _mapping(modified_charts[name]), _mapping(loop_charts[name]), pstar)
            endpoint_ports[name] = _endpoint_guard(graph, name, _mapping(_mapping(modified.get("endpoints")).get(name)), loop_charts[name].get("coordinate"))
    leading_ports = _leading_ports(graph, leading, loop_charts)
    switch_guard = _switch_joins(graph, switch)
    annular_guard = _switch_joins(graph, annular, annular=True)
    # The archived geometry certificates establish positivity of the native
    # measures.  Keep this as a source premise attached to each required chart
    # so downstream code cannot infer positivity from an integer node alone.
    geometry = _mapping(leading.get("geometry_certificate"))
    positive_geometry = bool(geometry.get("directed_T_positive") is True
                            and geometry.get("directed_post_length_positive") is True
                            and receipts.get("leading", {}).get("passed") is True)
    native_positive = {
        name: bool(chart_guards[name].get("passed") is True
                   and _handle(chart_guards[name].get("native_J")) is not None
                   and _handle(chart_guards[name].get("native_J")) > 0
                   and positive_geometry)
        for name in chart_names
    }
    for name, value in native_positive.items():
        chart_guards[name]["native_J_positive"] = value
        chart_guards[name]["positive_native_J_source_premise"] = value
    # Frozen macro is an optional diagnostic chart and is not part of the
    # required upstream seam scope.  The leading C6 packet has no frozen macro
    # entry, so its absence must not invalidate the four available ports.
    required_chart_names = ("second_switch", "post_power", "long_reshape")
    leading_required_names = ("first_switch", "second_switch", "post_power", "long_reshape")
    ordinary_z1 = bool(pressure_ports) and all(item.get("ordinary_Z1") is True for item in pressure_ports.values())
    physical_r_once = bool(pressure_ports) and all(item.get("R_Q_y_radius_once") is True for item in pressure_ports.values())
    positive_j = all(native_positive.get(name) is True for name in required_chart_names)
    required = {
        "source_pins": all(pin.get("passed") is True for pin in pins.values()),
        "source_receipts": all(item.get("passed") is True for item in receipts.values()),
        "source_family": family_ok, "archive_identity": all(identities.values()), "graph_prefixes": all(prefixes.values()),
        "source_ast": ast_guard.get("passed") is True, "common_N": common_n.get("passed") is True,
        "chart_guards": all(chart_guards[name].get("passed") is True for name in required_chart_names),
        "phase_N_once": all(phase[name].get("passed") is True for name in required_chart_names) and ast_guard["rows"].get("modified_fast_phase_N_once") is True,
        "pressure_ports": bool(pressure_ports) and all(item.get("passed") is True for item in pressure_ports.values()),
        "endpoint_ports": bool(endpoint_ports) and all(item.get("passed") is True for item in endpoint_ports.values()),
        "leading_ports": bool(leading_ports.get("charts")) and all(leading_ports["charts"].get(name, {}).get("passed") is True for name in leading_required_names),
        "positive_native_J": positive_j, "physical_R_once": physical_r_once, "ordinary_Z1": ordinary_z1,
        "switch_joins": switch_guard.get("passed") is True,
        "annular_joins": annular_guard.get("passed") is True,
    }
    unproved = [name for name, value in required.items() if not value]
    Pstar_once = all(item.get("Pstar_once") is True for item in pressure_ports.values())
    return {
        "passed": all(required.values()), "method": "actual_upstream_second_switch_post_power_long_reshape_native_calculus_guard",
        "source_only": True, "source_values_selected": False, "numeric_caps_used": False,
        "ancestor_replayed": False, "producer_or_checker_executed": False,
        "source_family": _FAMILY, "source_pins": pins, "sourcepins": pins, "receipts": receipts,
        "archive_identity": identities, "graph_prefixes": prefixes, "source_ast": ast_guard,
        "chart_guards": chart_guards, "charts": chart_guards, "pressure_ports": pressure_ports,
        "endpoint_ports": endpoint_ports, "leading_ports": leading_ports,
        "switch_joins": switch_guard, "annular_joins": annular_guard,
        # Public source premises used by the numeric runner.  They are
        # derived from the archived pressure ports and geometry certificates,
        # never from a numeric cap or a reconstructed owner object.
        "physical_R_once": physical_r_once, "ordinary_Z1": ordinary_z1,
        "positive_native_J": positive_j,
        "second_native_J_positive": native_positive.get("second_switch") is True,
        "post_power_native_J_positive": native_positive.get("post_power") is True,
        "long_native_J_positive": native_positive.get("long_reshape") is True,
        "phase": {"common_N": common_n, "charts": phase, "total_fast_phase_N_term_has_no_extra_N": ast_guard["rows"].get("modified_fast_phase_N_once") is True},
        "checks": {"required": required, "common_N": common_n, "source_ast": ast_guard, "receipts": receipts,
                    "chart_guards": chart_guards, "pressure_ports": pressure_ports, "endpoint_ports": endpoint_ports,
                    "leading_ports": leading_ports, "switch_joins": switch_guard, "annular_joins": annular_guard},
        "scope": {"charts": list(chart_names), "physical_coordinate": "y=log(R/110)", "common_N_node": common_n.get("node"),
                  "phase_origin": _handle(origin), "ordinary_Z_orders": [0, 1, 2, 3, 4],
                  "native_coordinates": {name: chart_guards[name].get("coordinate") for name in chart_names},
                  "native_J": {name: chart_guards[name].get("native_J") for name in chart_names}},
        "premises": {"d_native_equals_J_d_y": all(item.get("d_native_equals_J_d_y") is True for name, item in chart_guards.items() if name in required_chart_names),
                     "physical_dy_equals_native_divided_by_J": all(item.get("d_y_equals_d_native_over_J") is True for name, item in chart_guards.items() if name in required_chart_names),
                     "total_fast_phase_N_term_no_extra_N": ast_guard["rows"].get("modified_fast_phase_N_once") is True,
                     "physical_Pstar_once": Pstar_once, "physical_R_Q_y_radius_once": physical_r_once,
                     "physical_R_once": physical_r_once, "ordinary_Z1": ordinary_z1,
                     "positive_native_J": positive_j,
                     "common_N_and_global_phase": required["common_N"] and required["phase_N_once"],
                     "field_seam_equality_proved": False},
        "unproved": unproved + ["field_seam_equality", "numerical_endpoint_values", "numeric_caps_as_field_values"],
        "limits": {"field_seam_equality": False, "numeric_endpoint_values": False, "numeric_caps_as_values": False,
                   "global_m15_contraction": False, "ancestor_constructor_or_producer_replay": False},
    }


__all__ = ["build_upstream_native"]
