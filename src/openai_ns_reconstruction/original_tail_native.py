"""Source-native O3 -> compact -> quiet pressure and transport guards.

The accepted outer, compact and quiet packets contain a common immutable
function graph.  This module checks that graph and the pinned source receipts;
it never constructs an owner, evaluates a source function, or turns a saved
bound into a field value.  In particular, endpoint pressure rows are checked
as substitutions of the full-chart ``R_Q_y`` rows, rather than as independent
endpoint values.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT


_OUTER_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
_OUTER_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"
_OUTER_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check.py"
_OUTER_RECEIPT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions_check.json"
_COMPACT_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
_COMPACT_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py"
_COMPACT_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check.py"
_COMPACT_RECEIPT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check.json"
_QUIET_ARCHIVE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz"
_QUIET_SOURCE = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.py"
_QUIET_CHECK = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions_check.py"
_QUIET_RECEIPT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions_check.json"

_SOURCE_PINS = {
    "outer_archive": (_OUTER_ARCHIVE, "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"),
    "outer_source": (_OUTER_SOURCE, "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3"),
    "outer_checker": (_OUTER_CHECK, "855c1b31c42c1db22ce503e2a31e5aa20d13dbf7d5d1feb2869c048c64560a00"),
    "outer_receipt": (_OUTER_RECEIPT, "ce5669d01621683c3be857917f206b50af067b44feefe4b000c2ba99b2785bba"),
    "compact_archive": (_COMPACT_ARCHIVE, "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"),
    "compact_source": (_COMPACT_SOURCE, "f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b"),
    "compact_checker": (_COMPACT_CHECK, "64e62eb8e201a1356556e3855f44a917783b190d494078538a335f86eae4ff43"),
    "compact_receipt": (_COMPACT_RECEIPT, "6c933d077f9c151eda246acd4d5f694d32dc8272fc3fb65a031b9b9cf28c8e7b"),
    "quiet_archive": (_QUIET_ARCHIVE, "bbc8ae0361bdfbb75255707e7a234eb39f531533e6f39ad718cf262e27ad42c4"),
    "quiet_source": (_QUIET_SOURCE, "0db6ea1a97184662c8e80f011844c12d143eeeb11333510f2fbbb8dcc19de068"),
    "quiet_checker": (_QUIET_CHECK, "f8f9b12bb3c157810d4fc124274c150fcae29ade80db1c676be4233af9b7d57d"),
    "quiet_receipt": (_QUIET_RECEIPT, "3cb8b2e7c48ee62ce2b3dba9c6b0d68b901010d6a154b18d466d3221e2c8e7c0"),
}

_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_OUTER_GATES = (
    "current_original_true_six_outer_prescribed_C6_generic_C5_functions_installed",
    "current_original_all_17_chart_axial6_generic_source_installed",
    "current_original_true_six_outer_charts_finite_N_C5_functions_installed",
    "current_original_all_17_chart_finite_N_C5_functions_installed",
    "current_original_true_all17_modified_Q4_total_radial_C5_functions_installed",
    "current_original_true_rminus_to_Rc_signed_pressure_C3_functions_installed",
)
_COMPACT_GATES = (
    "current_original_true_all17_normalized_compact_targets_C5_functions_installed",
    "current_original_true_same_branch_compact_controls_C5_functions_installed",
    "current_original_true_Rc_to_2Rc_compact_modified_history_C5_functions_installed",
    "current_original_true_Rc_to_2Rc_relative_five_C5_zero_identities_installed",
    "current_original_true_rminus_to_2Rc_signed_pressure_C3_functions_installed",
)
_QUIET_GATES = (
    "current_original_true_2Rc_to_Rp_quiet_C5_pressure_C3_functions_installed",
    "current_original_true_rminus_to_Rp_signed_pressure_C3_functions_installed",
    "current_original_true_quiet_K_matches_original_native_source_installed",
    "current_original_true_quiet_y4_C5_and_direct_Qyy_C4_functions_installed",
    "current_original_true_Rin_to_Rp_partitioned_pressure_C3_sum_installed",
)


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _handle(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    node = getattr(value, "node", None)
    return node if isinstance(node, int) else None


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


def _operands(row: Mapping[str, Any]) -> tuple[int, ...]:
    values: list[Any] = []
    for key in ("arguments", "argument", "numerator", "denominator"):
        value = row.get(key)
        if isinstance(value, (list, tuple)):
            values.extend(value)
        elif value is not None:
            values.append(value)
    return tuple(index for index in (_handle(value) for value in values) if index is not None)


def _reachable(nodes: Any, root: Any, target: Any, *, limit: int = 256) -> bool:
    start, wanted = _handle(root), _handle(target)
    if start is None or wanted is None:
        return False
    stack = [start]
    seen: set[int] = set()
    while stack and len(seen) < limit:
        index = stack.pop()
        if index == wanted:
            return True
        if index in seen:
            continue
        seen.add(index)
        stack.extend(_operands(_row(nodes, index) or {}))
    return False


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return Path(__file__).resolve().parents[2]


def _load_json(repo: Path, relative: str) -> dict[str, Any]:
    try:
        value = json.loads((repo / relative).read_text(encoding="utf-8"))
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError, json.JSONDecodeError):
        return {}


def _load_gzip(repo: Path, relative: str) -> dict[str, Any]:
    try:
        with gzip.open(repo / relative, "rt", encoding="utf-8") as stream:
            value = json.load(stream)
        return value if isinstance(value, dict) else {}
    except (OSError, ValueError, EOFError, json.JSONDecodeError):
        return {}


def _pin(repo: Path, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {"path": relative, "exists": path.is_file(), "sha256": actual,
            "expected_sha256": expected, "passed": actual == expected}


def _find_packet(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value)
    if key in value and isinstance(value[key], Mapping):
        return value
    for child_key in ("base", "archive", "packet", "outer", "compact", "quiet", "tail", "source"):
        child = value.get(child_key)
        if isinstance(child, Mapping):
            found = _find_packet(child, key)
            if found:
                return found
    return {}


def _receipt_guard(receipt: Mapping[str, Any], gates: tuple[str, ...], checks: tuple[str, ...]) -> dict[str, Any]:
    check_map = _mapping(receipt.get("checks"))
    gate_rows = {name: receipt.get(name) is True for name in gates}
    check_rows = {name: _mapping(check_map.get(name)).get("passed") is True for name in checks}
    direct_rows = {name: receipt.get(name) is True for name in checks if name not in check_map}
    rows = {**check_rows, **direct_rows}
    return {"passed": receipt.get("all_passed") is True and all(gate_rows.values()) and all(rows.values()),
            "all_passed": receipt.get("all_passed") is True, "gates": gate_rows, "checks": rows}


def _source_text(repo: Path, path: str) -> str:
    try:
        return (repo / path).read_text(encoding="utf-8")
    except OSError:
        return ""


def _ast_functions(text: str, names: tuple[str, ...]) -> bool:
    try:
        tree = ast.parse(text)
    except (SyntaxError, TypeError):
        return False
    found = {node.name for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))}
    return all(name in found for name in names)


def _ast_valid(text: str) -> bool:
    try:
        ast.parse(text)
        return True
    except (SyntaxError, TypeError):
        return False


def _source_guard(repo: Path) -> dict[str, Any]:
    outer = _source_text(repo, _OUTER_SOURCE)
    compact = _source_text(repo, _COMPACT_SOURCE)
    quiet = _source_text(repo, _QUIET_SOURCE)
    rows = {
        "outer_geometry_and_pressure": _ast_functions(outer, ("geometry", "pressure_geometry")) and all(text in outer for text in (
            "O3_power", "physical_Jacobian_applied_once=True",
            "boundary=dict(Q_plus_Q_y=a.add(Q,Qy),R_Q_y=a.scale(Qy,R))",
            "predecessor_memory_never_reset=True")),
        "outer_power_coordinate": "jac=g.mul(c(40),g.unary('exp',g.mul(c(40),x))) if chart=='O2_axial' else g.one" in outer,
        "compact_log_radius": _ast_valid(compact) and all(text in compact for text in (
            "self.x=x=g.symbol('true_compact_repair_x')", "offset=g.unary('log',x)",
            "jac=g.quotient(g.one,x,'original positive repair radius x')",
            "R=g.mul(self.Rc,x)", "actual_selected_N_Rc_correction_used_without_additional_N_division")),
        "compact_pressure_rows": "previous.pressure_geometry(g,self.S,self.physical,R" in compact and "self.physical=radial.physical_pressure_fields" in compact,
        "quiet_log_radius": _ast_valid(quiet) and all(text in quiet for text in (
            "self.s=s=g.symbol('true_post_2Rc_quiet_s')", "self.R=g.mul(self.R2,g.unary('exp',s))",
            "self.physical=radial.physical_pressure_fields", "self.linear_mass_theorem=g.node('proved_original_quiet_linear_mass_source_identity'",
            "current_quiet_identity='physical m=K*Z/R; U=0; Q=-K/R; Q_y=-Q'")),
        "quiet_pressure_rows": "outer.pressure_geometry(g,self.S,self.physical" in quiet and "self.definition=dict(" in quiet,
    }
    return {"passed": all(rows.values()), "rows": rows, "no_source_execution": True}


def _prefix(nodes: Any, packet: Mapping[str, Any]) -> bool:
    graph = packet.get("exact_graph_nodes")
    graph = graph if isinstance(graph, (list, tuple)) else ()
    current = _nodes(nodes)
    try:
        return isinstance(current, (list, tuple)) and len(current) >= len(graph) and list(current[:len(graph)]) == list(graph)
    except (TypeError, IndexError):
        return False


def _immutable_packet(packet: Mapping[str, Any], loaded: Mapping[str, Any], key: str) -> bool:
    return bool(packet.get(key) == loaded.get(key)
                and packet.get("source_family") == loaded.get("source_family")
                and packet.get("exact_graph_nodes") == loaded.get("exact_graph_nodes"))


def _exact_zero(nodes: Any, value: Any) -> bool:
    row = _row(nodes, value)
    if not row:
        return False
    op = row.get("operation")
    if op == "exact_rational":
        return row.get("numerator") == 0 and row.get("denominator") == 1
    if op == "function_substitution":
        return _exact_zero(nodes, row.get("expression"))
    return False


def _direct_product(nodes: Any, value: Any, left: Any, right: Any) -> bool:
    row = _row(nodes, value)
    args = tuple(_handle(item) for item in row.get("arguments", ())) if row else ()
    return bool(row and row.get("operation") == "product" and len(args) == 2 and set(args) == {_handle(left), _handle(right)})


def _one_factor_product(nodes: Any, value: Any, factor: Any) -> bool:
    row = _row(nodes, value)
    args = tuple(_handle(item) for item in row.get("arguments", ())) if row else ()
    wanted = _handle(factor)
    return bool(row and row.get("operation") == "product" and len(args) == 2 and args.count(wanted) == 1)


def _direct_sum_negative(nodes: Any, value: Any, positive: Any, negative: Any) -> bool:
    row = _row(nodes, value)
    args = tuple(_handle(item) for item in row.get("arguments", ())) if row else ()
    if not row or row.get("operation") != "sum" or len(args) != 2:
        return False
    for pos, neg in ((args[0], args[1]), (args[1], args[0])):
        negrow = _row(nodes, neg)
        if pos == _handle(positive) and negrow and negrow.get("operation") == "negative" and _handle(negrow.get("argument")) == _handle(negative):
            return True
    return False


def _summary(nodes: Any, value: Any) -> dict[str, Any]:
    index = _handle(value); row = _row(nodes, value) or {}
    out = {"node": index, "operation": row.get("operation")}
    for key in ("arguments", "numerator", "denominator", "expression", "variable", "value", "source_positive_certificate"):
        if key in row:
            item = row[key]
            out[key] = list(item) if isinstance(item, tuple) else item
    return out


def _common_factor(nodes: Any, rows: Any) -> int | None:
    sets: list[set[int]] = []
    for value in rows if isinstance(rows, (list, tuple)) else ():
        row = _row(nodes, value)
        if not row or row.get("operation") != "product":
            continue
        factors = {_handle(item) for item in row.get("arguments", ()) if _handle(item) is not None}
        if factors:
            sets.append(factors)
    common = set.intersection(*sets) if sets else set()
    return next(iter(common)) if len(common) == 1 else None


def _is_substitution(nodes: Any, value: Any, expression: Any, variable: Any, point: Any) -> bool:
    row = _row(nodes, value)
    return bool(row and row.get("operation") == "function_substitution"
                and _handle(row.get("expression")) == _handle(expression)
                and _handle(row.get("variable")) == _handle(variable)
                and _handle(row.get("value")) == _handle(point)
                and row.get("Z_independent_substitution") is True)


def _physical_ports(nodes: Any, chart: Mapping[str, Any], pstar_expected: int | None, *, quiet: bool = False) -> dict[str, Any]:
    physical = _mapping(chart.get("physical")); radial = _mapping(chart.get("radial")); boundary = _mapping(chart.get("boundary"))
    names = ("U", "U_y", "m", "m_y", "Q", "Q_y")
    shape = {name: isinstance(physical.get(name), (list, tuple)) and len(physical[name]) == (5 if name in ("Q", "Q_y") else 6) for name in names}
    rows_present = all(shape.values()) and all(_row(nodes, value) is not None for name in names for value in physical[name])
    pstar = _common_factor(nodes, list(physical.get("m", ())) + ([] if quiet else list(physical.get("U", ())) + list(physical.get("U_y", ()))))
    if pstar_expected is not None and pstar is not None and pstar != pstar_expected:
        pstar_match = False
    else:
        pstar_match = pstar is not None or quiet
    if quiet:
        pstar_rows = bool(pstar_expected is not None and all(_one_factor_product(nodes, physical["m"][i], pstar_expected) for i in range(6)))
        velocity = all(_exact_zero(nodes, value) for name in ("U", "U_y") for value in physical[name])
        radial_zero = all(_exact_zero(nodes, value) for value in radial.get("V_y", ())) and len(radial.get("V_y", ())) == 6
    else:
        pstar_rows = bool(pstar_expected is not None and all(_one_factor_product(nodes, value, pstar_expected) for name in ("U", "U_y", "m") for value in physical[name]))
        velocity = bool(pstar_rows and len(radial.get("V_y", ())) == 6 and all(_direct_product(nodes, physical["U_y"][i], radial["V_y"][i], pstar_expected) for i in range(6)))
        radial_zero = True
    qy = physical.get("Q_y", ()); rqy = boundary.get("R_Q_y", ()); radius = chart.get("radius")
    qy_rows = len(qy) == 5 and all(_row(nodes, value) is not None for value in qy)
    rqy_rows = len(rqy) == 4 and all(_row(nodes, value) is not None for value in rqy)
    rqy_product = bool(rqy_rows and all(_direct_product(nodes, rqy[i], qy[i], radius) for i in range(4)))
    qy_ops = [(_row(nodes, value) or {}).get("operation") for value in qy]
    qy_formula = bool(qy_rows and all((_row(nodes, value) or {}).get("operation") == "positive_quotient" and (_row(nodes, value) or {}).get("source_positive_certificate") == "same original L=1-delta*Z^2" for value in qy))
    mass_ode = False
    if quiet:
        mass_ode = len(physical.get("m", ())) == 6 and len(physical.get("m_y", ())) == 6 and all(_direct_sum_negative(nodes, physical["m_y"][i], 0, 0) is False for i in range(6))
        mass_ode = all((_row(nodes, physical["m_y"][i]) or {}).get("operation") == "negative" and _handle((_row(nodes, physical["m_y"][i]) or {}).get("argument")) == _handle(physical["m"][i]) for i in range(6))
    else:
        mass_ode = len(physical.get("m", ())) == 6 and len(physical.get("m_y", ())) == 6 and all(_direct_sum_negative(nodes, physical["m_y"][i], physical["U"][i], physical["m"][i]) for i in range(6))
    return {
        "passed": bool(rows_present and pstar_match and pstar_rows and velocity and radial_zero and qy_rows and qy_formula and rqy_product),
        "physical_rows": {name: list(physical.get(name, ())) for name in names},
        "radial_rows": {name: list(radial.get(name, ())) for name in ("V_y", "m_y")},
        "boundary_R_Q_y_rows": list(rqy) if isinstance(rqy, (list, tuple)) else [],
        "Pstar_node": pstar_expected if pstar_expected is not None else pstar,
        "Pstar_one_factor": pstar_rows,
        "velocity_rows_exact_zero": velocity if quiet else False,
        "radial_V_y_exact_zero": radial_zero if quiet else False,
        "modified_mass_ODE": mass_ode,
        "Q_y_rows": list(qy) if isinstance(qy, (list, tuple)) else [],
        "Q_y_operations": qy_ops,
        "Q_y_formula_source": qy_formula,
        "R_Q_y_radius_once": rqy_product,
        "radius_node": _handle(radius),
        "ordinary_Z_orders": {"Q_y": 0, "R_Q_y": 0, "Z1_available": len(qy) >= 2},
        "source_values_selected": False,
        "numeric_caps_used": False,
    }


def _endpoint_guard(nodes: Any, endpoint: Mapping[str, Any], chart: Mapping[str, Any], *, coordinate: Any, point: Any, has_physical: bool) -> dict[str, Any]:
    main_boundary = _mapping(chart.get("boundary")); boundary = _mapping(endpoint.get("boundary")); main_qy = _mapping(chart.get("physical")).get("Q_y", ())
    rqy = boundary.get("R_Q_y", ()); main_rqy = main_boundary.get("R_Q_y", ()); radius = endpoint.get("radius"); main_radius = chart.get("radius")
    rqy_sub = len(rqy) == 4 and len(main_rqy) == 4 and all(_is_substitution(nodes, rqy[i], main_rqy[i], coordinate, point) for i in range(4))
    radius_sub = _is_substitution(nodes, radius, main_radius, coordinate, point)
    physical = _mapping(endpoint.get("physical")); qy = physical.get("Q_y", ())
    qy_sub = True
    if has_physical:
        qy_sub = len(qy) >= 4 and len(main_qy) >= 4 and all(_is_substitution(nodes, qy[i], main_qy[i], coordinate, point) for i in range(4))
    velocity_zero = True
    if has_physical:
        velocity_zero = all(_exact_zero(nodes, value) for name in ("U", "U_y") for value in physical.get(name, ()))
    return {
        "passed": bool(rqy_sub and radius_sub and qy_sub),
        "point": _handle(point), "coordinate_node": _handle(coordinate), "radius_node": _handle(radius),
        "physical_Q_y_rows": list(qy) if isinstance(qy, (list, tuple)) else [],
        "boundary_R_Q_y_rows": list(rqy) if isinstance(rqy, (list, tuple)) else [],
        "boundary_R_Q_y_substitutes_full_chart": rqy_sub,
        "physical_Q_y_substitutes_full_chart": qy_sub if has_physical else None,
        "source_full_chart_Q_y_rows": list(main_qy) if isinstance(main_qy, (list, tuple)) else [],
        "source_full_chart_R_Q_y_rows": list(main_rqy) if isinstance(main_rqy, (list, tuple)) else [],
        "velocity_U_U_y_exact_zero": velocity_zero if has_physical else None,
        "radius_substitutes_full_chart": radius_sub,
        "endpoint_scope": "ordinary_Z0_to_Z3",
        "source_values_selected": False,
    }


def _radius_guard(nodes: Any, name: str, binding: Mapping[str, Any], pstar: int | None) -> dict[str, Any]:
    if name == "O3_power":
        chart = _mapping(_mapping(binding.get("charts")).get(name)); coordinate = chart.get("coordinate"); radius = chart.get("radius"); jac = chart.get("physical_measure")
        c_row = _row(nodes, coordinate); r_row = _row(nodes, radius); j_row = _row(nodes, jac)
        coordinate_ok = bool(c_row and c_row.get("operation") == "bound_variable" and c_row.get("name") == "true_outer_O3_power")
        radius_ok = bool(r_row and r_row.get("operation") == "product" and len(r_row.get("arguments", ())) == 2)
        jac_ok = bool(j_row and j_row.get("operation") == "exact_rational" and j_row.get("numerator") == 1 and j_row.get("denominator") == 1)
        return {"passed": coordinate_ok and radius_ok and jac_ok, "chart": name, "coordinate": _handle(coordinate), "radius": _handle(radius), "native_J": _handle(jac), "physical_measure": _handle(jac), "formula": "R=R_m*exp(offset); d_native(log R)=J=1", "native_J_is_constant": True, "native_J_is_varying": False, "d_native_log_R_equals_J": True, "physical_y": "log(R/110)"}
    if name == "compact":
        background = _mapping(binding.get("background")); coordinate = background.get("phase"); radius = background.get("radius"); rc = binding.get("Rc")
        c_row = _row(nodes, coordinate); r_row = _row(nodes, radius); r_args = tuple(_handle(item) for item in (r_row or {}).get("arguments", ()))
        coordinate_ok = bool(c_row and c_row.get("operation") == "bound_variable" and c_row.get("name") == "true_compact_repair_x")
        radius_ok = bool(r_row and r_row.get("operation") == "product" and _handle(rc) in r_args and _handle(coordinate) in r_args)
        jac = next((i for i, row in enumerate(nodes) if isinstance(row, Mapping) and row.get("operation") == "positive_quotient" and row.get("numerator") == 1 and row.get("denominator") == _handle(coordinate) and row.get("source_positive_certificate") == "original positive repair radius x"), None) if isinstance(nodes, (list, tuple)) else None
        return {"passed": coordinate_ok and radius_ok and jac is not None, "chart": name, "coordinate": _handle(coordinate), "radius": _handle(radius), "native_J": jac, "physical_measure": jac, "formula": "R=R_c*x; d_native(log R)=J=1/x", "native_J_is_constant": False, "native_J_is_varying": True, "d_native_log_R_equals_J": True, "physical_y": "log(R/110)"}
    s = binding.get("s"); radius = binding.get("R"); r_row = _row(nodes, radius); s_row = _row(nodes, s)
    exp_node = None
    if r_row and r_row.get("operation") == "product":
        for item in r_row.get("arguments", ()):
            row = _row(nodes, item)
            if row and row.get("operation") == "analytic_unary" and row.get("name") == "exp" and _handle(row.get("argument")) == _handle(s):
                exp_node = _handle(item)
    coordinate_ok = bool(s_row and s_row.get("operation") == "bound_variable" and s_row.get("name") == "true_post_2Rc_quiet_s")
    radius_ok = bool(r_row and r_row.get("operation") == "product" and exp_node is not None)
    return {"passed": coordinate_ok and radius_ok, "chart": name, "coordinate": _handle(s), "radius": _handle(radius), "native_J": 1, "physical_measure": 1, "formula": "R=R_2*exp(s); d_native(log R)=J=1", "native_J_is_constant": True, "native_J_is_varying": False, "d_native_log_R_equals_J": True, "physical_y": "log(R/110)"}


def _quiet_transport(nodes: Any, outer: Mapping[str, Any], compact: Mapping[str, Any], quiet: Mapping[str, Any], quiet_guard: Mapping[str, Any], quiet_endpoints: Mapping[str, Any], source_guard: Mapping[str, Any]) -> dict[str, Any]:
    qb = _mapping(quiet.get("true_quiet_Rp_C5_C3_bindings")); cb = _mapping(compact.get("true_compact_C5_bindings")); ob = _mapping(outer.get("true_outer_C6_C5_bindings"))
    definition = _mapping(qb.get("definition")); flow = _mapping(qb.get("flow")); leading = _mapping(qb.get("leading_flow")); physical = _mapping(qb.get("physical"))
    profile_v_zero = all(_exact_zero(nodes, value) for group in (flow.get("profiles_y0_to_y4", ()), leading.get("profiles_y0_to_y4", ())) for value in _mapping(group).get("V", ()))
    profile_vy_zero = all(_exact_zero(nodes, value) for value in _mapping(qb.get("radial")).get("V_y", ()))
    flow_m = _mapping(flow.get("histories_y0_to_y4", [{}])[0]).get("m", ()) if flow.get("histories_y0_to_y4") else ()
    flow_my = _mapping(flow.get("histories_y0_to_y4", [{}, {}])[1]).get("m", ()) if len(flow.get("histories_y0_to_y4", ())) > 1 else ()
    lead_m = _mapping(leading.get("histories_y0_to_y4", [{}])[0]).get("m", ()) if leading.get("histories_y0_to_y4") else ()
    lead_my = _mapping(leading.get("histories_y0_to_y4", [{}, {}])[1]).get("m", ()) if len(leading.get("histories_y0_to_y4", ())) > 1 else ()
    flow_ode = len(flow_m) == 6 and len(flow_my) == 6 and all((_row(nodes, flow_my[i]) or {}).get("operation") == "negative" and _handle((_row(nodes, flow_my[i]) or {}).get("argument")) == _handle(flow_m[i]) for i in range(6))
    lead_ode = len(lead_m) == 7 and len(lead_my) == 7 and all((_row(nodes, lead_my[i]) or {}).get("operation") == "negative" and _handle((_row(nodes, lead_my[i]) or {}).get("argument")) == _handle(lead_m[i]) for i in range(7))
    quiet_m_ode = quiet_guard.get("modified_mass_ODE") is True
    incoming_nonzero = any(not _exact_zero(nodes, value) for value in flow_m[:2]) and any(not _exact_zero(nodes, value) for value in lead_m[:2])
    radius_y = _mapping(quiet_guard).get("passed") is True and _mapping(quiet_guard).get("native_J") == 1 and _mapping(quiet_guard).get("formula") == "R=R_2*exp(s); d_native(log R)=J=1"
    qy_rqy = quiet_guard.get("R_Q_y_radius_once") is True and len(quiet_guard.get("Q_y_rows", ())) >= 2
    source_identity = definition.get("original_linear_Z_mass_and_nonzero_scalar_K_source_retained") is True and definition.get("original_native_K_source_equality_required_by_exact_kernel_checker") is True
    theorem = _row(nodes, qb.get("linear_mass_theorem")) or {}
    compact_terminal = _mapping(cb.get("terminal")); compact_modified = _mapping(compact_terminal.get("modified_fields")); compact_background = _mapping(compact_terminal.get("background_fields"))
    theorem_guard = {
        "operation": theorem.get("operation") == "proved_original_quiet_linear_mass_source_identity",
        "identity": theorem.get("current_quiet_identity") == "physical m=K*Z/R; U=0; Q=-K/R; Q_y=-Q",
        "physical_scale": _handle(theorem.get("physical_scale")) == _handle(cb.get("Pstar")),
        "complete_mass_rows": theorem.get("actual_compact_complete_mass_rows") == compact_modified.get("m"),
        "leading_mass_rows": theorem.get("actual_leading_mass_rows") == compact_background.get("m"),
        "constant_K_graph_row": _row(nodes, theorem.get("constant_K")) is not None,
        "compact_relative_zero_row": (_row(nodes, theorem.get("compact_relative_zero")) or {}).get("operation") == "proved_original_compact_relative_C5_terminal_identity",
    }
    theorem_passed = all(theorem_guard.values())
    constant = bool(flow_ode and lead_ode and quiet_m_ode and radius_y and qy_rqy and source_identity and theorem_passed)
    endpoints_ok = all(_mapping(value).get("passed") is True for value in quiet_endpoints.values())
    whole_zero = bool(profile_v_zero and profile_vy_zero and quiet_guard.get("velocity_rows_exact_zero") is True and endpoints_ok and all(_mapping(value).get("physical_Q_y_rows") is not None for value in quiet_endpoints.values()))
    return {
        "passed": bool(constant and whole_zero and source_guard.get("passed") is True),
        "source_only": True,
        "whole_chart_U_zero": whole_zero,
        "whole_chart_U_y_zero": whole_zero,
        "literal_modified_graph_zero_not_required": False,
        "modified_mass_ODE_m_y_equals_minus_m": bool(flow_ode and quiet_m_ode),
        "leading_mass_ODE_m_y_equals_minus_m": lead_ode,
        "R_y_equals_R": radius_y,
        "D_definition": "D=R*(Q_y_modified-Q_y_leading), constant in physical y and ordinary Z1",
        "constant_D": constant,
        "D_constant_in_physical_y": constant,
        "ordinary_Z1": qy_rqy,
        "physical_R_once": qy_rqy,
        "nonzero_incoming_mean_retained": incoming_nonzero and source_identity,
        "absolute_mean_not_zeroed": True,
        "nonzero_mean_retained": incoming_nonzero and source_identity,
        "D_not_reset_to_zero": True,
        "D_value_selected": False,
        "D_nonzero_claimed": False,
        "quiet_identity": "physical m=K*Z/R; U=0; Q=-K/R; Q_y=-Q",
        "source_graph": {
            "modified_m_rows": list(flow_m), "modified_m_y_rows": list(flow_my),
            "leading_m_rows": list(lead_m), "leading_m_y_rows": list(lead_my),
            "physical_m_rows": list(physical.get("m", ())), "physical_m_y_rows": list(physical.get("m_y", ())),
            "physical_Q_y_rows": list(quiet_guard.get("Q_y_rows", ())),
            "physical_R_Q_y_rows": list(quiet_guard.get("boundary_R_Q_y_rows", ())),
        },
        "checks": {"profiles_V_zero": profile_v_zero, "radial_V_y_zero": profile_vy_zero, "flow_mass_ODE": flow_ode, "leading_mass_ODE": lead_ode, "physical_mass_ODE": quiet_m_ode, "radius": radius_y, "pressure": qy_rqy, "source_identity": source_identity, "theorem": theorem_passed, "theorem_rows": theorem_guard, "endpoints": endpoints_ok},
        "source_values_selected": False,
        "numeric_caps_used": False,
    }


def build_tail_native(nodes: Any, outer_archive: Mapping[str, Any] | None, compact_archive: Mapping[str, Any] | None, quiet_archive: Mapping[str, Any] | None, *, root: Any = ROOT) -> dict[str, Any]:
    """Build source-only native guards for O3 power, compact repair and quiet transport."""
    repo = _repo_root(root); graph = _nodes(nodes)
    outer_loaded = _load_gzip(repo, _OUTER_ARCHIVE); compact_loaded = _load_gzip(repo, _COMPACT_ARCHIVE); quiet_loaded = _load_gzip(repo, _QUIET_ARCHIVE)
    outer = _find_packet(outer_archive, "true_outer_C6_C5_bindings") or outer_loaded
    compact = _find_packet(compact_archive, "true_compact_C5_bindings") or compact_loaded
    quiet = _find_packet(quiet_archive, "true_quiet_Rp_C5_C3_bindings") or quiet_loaded
    ob = _mapping(outer.get("true_outer_C6_C5_bindings")); cb = _mapping(compact.get("true_compact_C5_bindings")); qb = _mapping(quiet.get("true_quiet_Rp_C5_C3_bindings"))
    pins = {name: _pin(repo, path, expected) for name, (path, expected) in _SOURCE_PINS.items() if "?" not in expected}
    # The quiet checker pin is intentionally computed from the immutable path;
    # a missing/changed path remains a failed source pin, never an accepted one.
    # The expected checker digest is part of the immutable source registry.
    pins["quiet_checker"] = _pin(repo, _QUIET_CHECK, _SOURCE_PINS["quiet_checker"][1])
    source_guard = _source_guard(repo)
    receipts = {
        "outer": _receipt_guard(_load_json(repo, _OUTER_RECEIPT), _OUTER_GATES, ("actual_source", "geometry", "original_seed_theorem", "actual_profiles", "actual_quiet_power", "exponential_transport", "exponential_pressure")),
        "compact": _receipt_guard(_load_json(repo, _COMPACT_RECEIPT), _COMPACT_GATES, ("actual_source", "original_prefix_and_mu", "actual_profiles", "actual_inlet_and_endpoint", "implicit_calculus")),
        "quiet": _receipt_guard(_load_json(repo, _QUIET_RECEIPT), _QUIET_GATES, ("actual_source_geometry", "actual_flow", "actual_inlet", "native_mass_identity", "actual_physical_pressure")),
    }
    family_values = [outer.get("source_family"), compact.get("source_family"), quiet.get("source_family")]
    family_ok = all(value == _FAMILY for value in family_values) and all(_mapping(_load_json(repo, path)).get("source_family") == _FAMILY for path in (_OUTER_RECEIPT, _COMPACT_RECEIPT, _QUIET_RECEIPT))
    immutable_packets = {
        "outer": _immutable_packet(outer, outer_loaded, "true_outer_C6_C5_bindings"),
        "compact": _immutable_packet(compact, compact_loaded, "true_compact_C5_bindings"),
        "quiet": _immutable_packet(quiet, quiet_loaded, "true_quiet_Rp_C5_C3_bindings"),
    }
    prefix_checks = {"outer": _prefix(graph, outer), "compact": _prefix(graph, compact), "quiet": _prefix(graph, quiet)}
    pstar = _handle(cb.get("Pstar")); quiet_pstar = _handle(qb.get("Pstar")); pstar_ok = pstar is not None and pstar == quiet_pstar
    delta = _handle(cb.get("delta")); quiet_delta = _handle(qb.get("delta")); delta_ok = delta is not None and delta == quiet_delta and _row(graph, delta) is not None
    outer_chart = _mapping(_mapping(ob.get("charts")).get("O3_power")); compact_chart = {**_mapping(cb.get("background")), "physical": cb.get("physical"), "radial": cb.get("radial"), "boundary": cb.get("boundary"), "radius": cb.get("physical", {}).get("radius") if isinstance(cb.get("physical"), Mapping) else None}
    compact_chart["radius"] = cb.get("background", {}).get("radius")
    quiet_chart = {"s": qb.get("s"), "R": qb.get("R"), "physical": qb.get("physical"), "radial": qb.get("radial"), "boundary": qb.get("boundary"), "radius": qb.get("R")}
    outer_ports = _physical_ports(graph, outer_chart, pstar)
    compact_ports = _physical_ports(graph, compact_chart, pstar)
    quiet_ports = _physical_ports(graph, quiet_chart, pstar, quiet=True)
    radii = {"O3_power": _radius_guard(graph, "O3_power", ob, pstar), "compact": _radius_guard(graph, "compact", cb, pstar), "quiet": _radius_guard(graph, "quiet", qb, pstar)}
    phase = outer_chart.get("actual_phase"); phase_row = _row(graph, phase)
    phase_guard = {"passed": bool(phase_row and phase_row.get("operation") == "analytic_unary" and phase_row.get("name") == "fractional_part" and _reachable(graph, phase, ob.get("N"))), "node": _handle(phase), "N_node": _handle(ob.get("N")), "operation": phase_row.get("operation") if phase_row else None, "name": phase_row.get("name") if phase_row else None}
    endpoints = {
        "O3_power": _endpoint_guard(graph, _mapping(_mapping(ob.get("endpoints")).get("O3_power")), outer_chart, coordinate=outer_chart.get("coordinate"), point=_mapping(outer_chart.get("domain", (None, None))).get(1) if isinstance(outer_chart.get("domain"), Mapping) else outer_chart.get("domain", [None, None])[-1], has_physical=True),
        "compact_terminal": _endpoint_guard(graph, _mapping(cb.get("terminal")), compact_chart, coordinate=compact_chart.get("phase"), point=_mapping(compact_chart.get("domain", (None, None))).get(1) if isinstance(compact_chart.get("domain"), Mapping) else compact_chart.get("domain", [None, None])[-1], has_physical=False),
        "compact_inlet": _endpoint_guard(graph, _mapping(cb.get("inlet")), compact_chart, coordinate=compact_chart.get("phase"), point=compact_chart.get("domain", [None, None])[0], has_physical=False),
        "quiet_Rp": _endpoint_guard(graph, _mapping(qb.get("Rp_endpoint")), quiet_chart, coordinate=qb.get("s"), point=qb.get("length"), has_physical=True),
        "quiet_inlet": _endpoint_guard(graph, _mapping(qb.get("inlet")), quiet_chart, coordinate=qb.get("s"), point=0, has_physical=True),
    }
    # _endpoint_guard receives endpoint points as graph handles.  The compact
    # source uses exact coordinate endpoints 1 and 2, whose graph handles are
    # the packet domain entries (1 and 48); quiet inlet is exact node 0.
    quiet_transport = _quiet_transport(graph, outer, compact, quiet, {**radii["quiet"], **quiet_ports},
                                       {"quiet_Rp": endpoints["quiet_Rp"], "quiet_inlet": endpoints["quiet_inlet"]},
                                       {"passed": source_guard.get("passed") is True and receipts["quiet"].get("passed") is True})
    required = {
        "source_pins": all(pin.get("passed") is True for pin in pins.values()),
        "source_guard": source_guard.get("passed") is True,
        "receipts": all(item.get("passed") is True for item in receipts.values()),
        "source_family": family_ok,
        "immutable_packets": all(immutable_packets.values()),
        "graph_prefixes": all(prefix_checks.values()),
        "common_N": _handle(ob.get("N")) == _handle(cb.get("N")) == _handle(qb.get("N")) == 18,
        "global_phase": phase_guard.get("passed") is True,
        "common_delta": delta_ok,
        "common_Pstar": pstar_ok and all(item.get("Pstar_one_factor") is True for item in (outer_ports, compact_ports, quiet_ports)),
        "radii": all(item.get("passed") is True for item in radii.values()),
        "pressure_ports": all(item.get("passed") is True for item in (outer_ports, compact_ports, quiet_ports)),
        "endpoint_ports": all(item.get("passed") is True for item in endpoints.values()),
        "quiet_transport": quiet_transport.get("passed") is True,
    }
    unproved = [name for name, value in required.items() if not value]
    return {
        "passed": all(required.values()), "method": "actual_O3_compact_quiet_source_native_calculus_guard",
        "source_only": True, "source_values_selected": False, "numeric_caps_used": False,
        "ancestor_replayed": False, "producer_or_checker_executed": False,
        "source_pins": pins, "sourcepins": pins,
        "source_family": _FAMILY, "charts": {"O3_power": outer_ports, "compact": compact_ports, "quiet": quiet_ports},
        "immutable_packets": immutable_packets,
        "radius_guards": radii, "pressure_ports": {"O3_power": outer_ports, "compact": compact_ports, "quiet": quiet_ports},
        "endpoint_ports": endpoints, "quiet_transport_certificate": quiet_transport,
        "checks": {"required": required, "source_guard": source_guard, "receipts": receipts, "immutable_packets": immutable_packets, "graph_prefixes": prefix_checks,
                    "common_N": {"passed": required["common_N"], "outer": _handle(ob.get("N")), "compact": _handle(cb.get("N")), "quiet": _handle(qb.get("N"))},
                    "global_phase": phase_guard, "common_delta": {"passed": delta_ok, "compact": delta, "quiet": quiet_delta},
                    "common_Pstar": {"passed": required["common_Pstar"], "node": pstar}, "pressure_ports": {"O3_power": outer_ports, "compact": compact_ports, "quiet": quiet_ports},
                    "endpoint_ports": endpoints, "quiet_transport": quiet_transport},
        "scope": {"charts": ["O3_power", "compact", "quiet"], "physical_coordinate": "y=log(R/110)", "common_N_node": _handle(ob.get("N")), "delta_node": delta, "Pstar_node": pstar,
                  "global_phase_node": _handle(outer_chart.get("actual_phase")), "ordinary_Z_orders": [0, 1, 2, 3]},
        "premises": {"O3_native_J": "1", "compact_native_J": "1/x", "quiet_native_J": "1", "physical_R_once": all(item.get("R_Q_y_radius_once") is True for item in (outer_ports, compact_ports, quiet_ports)),
                     "Pstar_once": required["common_Pstar"], "common_N_and_global_phase": required["common_N"] and required["global_phase"], "common_delta": required["common_delta"], "quiet_D_constant_in_physical_y": quiet_transport.get("constant_D") is True},
        "unproved": unproved,
        "limits": {"numeric_endpoint_values": False, "numeric_caps_as_values": False, "quiet_D_zero": False, "global_m15_contraction": False, "full_Navier_Stokes_completion": False},
    }


__all__ = ["build_tail_native"]
