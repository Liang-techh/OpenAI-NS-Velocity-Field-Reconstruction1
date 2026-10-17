"""Source-only compact-tail mass identities.

This guard reads the accepted compact C5 graph and its source receipt.  It
does not construct an owner, execute a producer or checker, evaluate a field,
or use a numerical cap.  The three small algebra checks below deliberately
keep the amplitude, incoming correction, common ``N``, target jets and fixed
point jets as opaque symbols.  Only the graph's exact rational/sum/product/
negative/positive-quotient operations are evaluated symbolically.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import inspect
import json
import math
import textwrap
from pathlib import Path
from typing import Any, Mapping

import sympy as _sp

from .original_core_numeric import ROOT


_COMPACT_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
)
_COMPACT_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py"
)
_COMPACT_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check.py"
)
_COMPACT_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check.json"
)
_BUMP_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_repair_band_cone_budget.py"
)

_SOURCE_PINS = {
    "compact_archive": (_COMPACT_ARCHIVE, "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"),
    "compact_source": (_COMPACT_SOURCE, "f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b"),
    "compact_checker": (_COMPACT_CHECK, "64e62eb8e201a1356556e3855f44a917783b190d494078538a335f86eae4ff43"),
    "compact_receipt": (_COMPACT_RECEIPT, "6c933d077f9c151eda246acd4d5f694d32dc8272fc3fb65a031b9b9cf28c8e7b"),
    "original_bump_source": (_BUMP_SOURCE, "4e07cb683841dd90bb14da382da90c83d6639bfa19dd037cd03e72b34fba4c22"),
}

_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

_GATES = (
    "current_original_true_all17_normalized_compact_targets_C5_functions_installed",
    "current_original_true_same_branch_compact_controls_C5_functions_installed",
    "current_original_true_Rc_to_2Rc_compact_modified_history_C5_functions_installed",
    "current_original_true_Rc_to_2Rc_relative_five_C5_zero_identities_installed",
    "current_original_true_rminus_to_2Rc_signed_pressure_C3_functions_installed",
)
_RECEIPT_CHECKS = (
    "actual_source",
    "original_prefix_and_mu",
    "actual_profiles",
    "actual_inlet_and_endpoint",
    "implicit_calculus",
)


class _UnsupportedGraphOperation(Exception):
    """Raised when the bounded proof reaches an operation it was not given."""


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
    if isinstance(value, Mapping) and isinstance(value.get("exact_graph_nodes"), (list, tuple, Mapping)):
        return value["exact_graph_nodes"]
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
    except (OSError, ValueError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _load_gzip(repo: Path, relative: str) -> dict[str, Any]:
    try:
        with gzip.open(repo / relative, "rt", encoding="utf-8") as stream:
            value = json.load(stream)
    except (OSError, ValueError, EOFError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def _pin(repo: Path, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    return {
        "path": relative,
        "exists": path.is_file(),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
    }


def _source_text(repo: Path, relative: str) -> str:
    try:
        return (repo / relative).read_text(encoding="utf-8")
    except OSError:
        return ""


def _ast_valid(text: str) -> bool:
    try:
        ast.parse(text)
    except (SyntaxError, TypeError):
        return False
    return True


def _ast_function_names(text: str) -> set[str]:
    try:
        tree = ast.parse(text)
    except (SyntaxError, TypeError):
        return set()
    return {
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }


def _binding_ast_hash(text: str, function_name: str) -> str | None:
    """Match the source registry's ``ast_binding`` hash convention."""
    try:
        tree = ast.parse(text)
    except (SyntaxError, TypeError):
        return None
    functions = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function_name
    ]
    if len(functions) != 1:
        return None
    segment = ast.get_source_segment(text, functions[0])
    if segment is None:
        return None
    normalized = textwrap.dedent(segment)
    return hashlib.sha256(ast.dump(ast.parse(normalized)).encode()).hexdigest()


def _find_packet(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value)
    if key in value and isinstance(value[key], Mapping):
        return value
    for child_key in ("packet", "archive", "compact", "base", "source"):
        child = value.get(child_key)
        if isinstance(child, Mapping):
            found = _find_packet(child, key)
            if found:
                return found
    return {}


def _compact_bindings(value: Any) -> Mapping[str, Any]:
    packet = _find_packet(value, "true_compact_C5_bindings")
    if packet:
        return _mapping(packet.get("true_compact_C5_bindings"))
    value = _mapping(value)
    return value if "N" in value and "targets" in value else {}


def _packet_graph(value: Any) -> Any:
    packet = _find_packet(value, "true_compact_C5_bindings")
    if packet and isinstance(packet.get("exact_graph_nodes"), (list, tuple, Mapping)):
        return packet["exact_graph_nodes"]
    value = _mapping(value)
    return value.get("exact_graph_nodes", ())


def _receipt_guard(receipt: Mapping[str, Any]) -> dict[str, Any]:
    checks = _mapping(receipt.get("checks"))
    gates = {name: receipt.get(name) is True for name in _GATES}
    check_rows = {
        name: _mapping(checks.get(name)).get("passed") is True
        for name in _RECEIPT_CHECKS
    }
    inlet = _mapping(checks.get("actual_inlet_and_endpoint"))
    theorem = inlet.get("terminal_integral_identity_uses_pinned_original_exact_bump_theorem") is True
    rows = {**gates, **check_rows, "pinned_original_bump_theorem": theorem}
    return {
        "passed": receipt.get("all_passed") is True and all(rows.values()),
        "all_passed": receipt.get("all_passed") is True,
        "gates": gates,
        "checks": check_rows,
        "pinned_original_bump_theorem": theorem,
        "mathematical_fixture_only": inlet.get("mathematical_fixture_only") is True,
    }


def _source_guard(repo: Path, compact: Mapping[str, Any]) -> dict[str, Any]:
    source = _source_text(repo, _COMPACT_SOURCE)
    checker = _source_text(repo, _COMPACT_CHECK)
    bump = _source_text(repo, _BUMP_SOURCE)
    binding = _mapping(_mapping(compact.get("bindings")).get("original_partial_history_terminal_identity"))
    bump_module = binding.get("module") == Path(_BUMP_SOURCE).name
    bump_sha = binding.get("module_sha256") == _SOURCE_PINS["original_bump_source"][1]
    bump_fn = binding.get("function") == "exact_partial_band_graph"
    bump_ast = _binding_ast_hash(bump, "exact_partial_band_graph")
    bump_ast_ok = bump_ast is not None and bump_ast == binding.get("AST_sha256")
    source_rows = {
        "compact_source_ast": _ast_valid(source),
        "compact_checker_ast": _ast_valid(checker),
        "compact_source_functions": {name: name in _ast_function_names(source) for name in (
            "normalized_targets", "implicit_rows", "bump_rows", "density", "profiles",
        )},
        "normalized_target_formula": all(text in source for text in (
            "scaled={k:a.scale(v,N)",
            "d={row:a.div(scaled[k],A if degree==1 else AA",
            "joint=a.sub(scaled['k'],a.mul(A,scaled['m']))",
        )),
        "fixedpoint_source_formula": all(text in source for text in (
            "self.system=g.node('original_compact_unique_Banach_limit_system'",
            "self.h,self.implicit=implicit_rows",
            "self.residual={key:a.add",
            "self.terminal_rhs={}",
            "self.terminal_certificates={key:g.node('proved_original_compact_relative_C5_terminal_identity'",
        )),
        "original_bump_source_ast": _ast_valid(bump),
        "original_bump_function": "exact_partial_band_graph" in _ast_function_names(bump),
        "original_bump_binding_module": bump_module,
        "original_bump_binding_sha256": bump_sha,
        "original_bump_binding_function": bump_fn,
        "original_bump_binding_ast": bump_ast_ok,
    }
    source_rows["passed"] = all(
        value if isinstance(value, bool) else all(value.values())
        for key, value in source_rows.items()
        if key != "passed"
    )
    return {"passed": source_rows["passed"], "rows": source_rows, "no_source_execution": True}


def _prefix(graph: Any, packet: Mapping[str, Any]) -> bool:
    expected = packet.get("exact_graph_nodes")
    if not isinstance(expected, (list, tuple)) or not isinstance(graph, (list, tuple)):
        return False
    try:
        return len(graph) >= len(expected) and list(graph[: len(expected)]) == list(expected)
    except (TypeError, IndexError):
        return False


def _zero_row(nodes: Any, value: Any) -> bool:
    row = _row(nodes, value)
    if not row:
        return False
    if row.get("operation") == "exact_rational":
        try:
            return _sp.Rational(row.get("numerator"), row.get("denominator")) == 0
        except (TypeError, ValueError):
            return False
    if row.get("operation") == "function_substitution":
        return _zero_row(nodes, row.get("expression"))
    return False


class _SymbolicGraph:
    """Small fail-closed evaluator for the four admitted graph operations."""

    _OPS = {"exact_rational", "sum", "product", "negative", "positive_quotient", "analytic_unary"}

    def __init__(self, nodes: Any, opaque: Mapping[int, Any]):
        self.nodes = nodes
        self.opaque = dict(opaque)
        self.cache: dict[int, _sp.Expr] = {}

    def eval(self, value: Any) -> _sp.Expr:
        index = _handle(value)
        if index is None:
            raise _UnsupportedGraphOperation(f"non-node handle {value!r}")
        if index in self.opaque:
            return self.opaque[index]
        if index in self.cache:
            return self.cache[index]
        row = _row(self.nodes, index)
        if row is None:
            raise _UnsupportedGraphOperation(f"missing node {index}")
        operation = row.get("operation")
        if operation not in self._OPS:
            raise _UnsupportedGraphOperation(f"node {index}: {operation!r}")
        if operation == "exact_rational":
            value_expr = _sp.Rational(row.get("numerator"), row.get("denominator"))
        elif operation == "sum":
            value_expr = sum((self.eval(item) for item in row.get("arguments", ())), _sp.Integer(0))
        elif operation == "product":
            value_expr = _sp.Integer(1)
            for item in row.get("arguments", ()):
                value_expr *= self.eval(item)
        elif operation == "negative":
            value_expr = -self.eval(row.get("argument"))
        elif operation == "positive_quotient":
            value_expr = self.eval(row.get("numerator")) / self.eval(row.get("denominator"))
        else:
            name = row.get("name")
            argument = self.eval(row.get("argument"))
            if name == "log":
                value_expr = _sp.log(argument)
            elif name == "exp":
                value_expr = _sp.exp(argument)
            else:
                raise _UnsupportedGraphOperation(f"node {index}: analytic_unary {name!r}")
        value_expr = _sp.powsimp(_sp.cancel(value_expr), force=True)
        self.cache[index] = value_expr
        return value_expr


def _expr_string(value: _sp.Expr) -> str:
    return _sp.sstr(_sp.powsimp(_sp.cancel(value), force=True))


def _symbol(name: str) -> _sp.Symbol:
    return _sp.Symbol(name, commutative=True)


def _ordinary_product(values_left: list[_sp.Expr], values_right: list[_sp.Expr], order: int) -> _sp.Expr:
    return sum(
        _sp.Integer(math.comb(order, j)) * values_left[j] * values_right[order - j]
        for j in range(order + 1)
    )


def _identity_incoming_target(nodes: Any, compact: Mapping[str, Any]) -> dict[str, Any]:
    amplitude = list(compact.get("A", ()))
    incoming = list(_mapping(compact.get("incoming_correction")).get("m", ()))
    scaled = list(_mapping(compact.get("N_scaled_incoming")).get("m", ()))
    targets = list(_mapping(compact.get("targets")).get("M", ()))
    n_node = _handle(compact.get("N"))
    required = len(amplitude) >= 4 and len(incoming) >= 4 and len(scaled) >= 4 and len(targets) >= 4 and n_node is not None
    if not required:
        return {"passed": False, "orders": [], "reason": "missing actual A, incoming m, N or target M rows"}
    opaque = {n_node: _symbol("N")}
    opaque.update({_handle(value): _symbol(f"A{index}") for index, value in enumerate(amplitude) if _handle(value) is not None})
    opaque.update({_handle(value): _symbol(f"m{index}") for index, value in enumerate(incoming) if _handle(value) is not None})
    evaluator = _SymbolicGraph(nodes, opaque)
    rows: list[dict[str, Any]] = []
    try:
        for order in range(4):
            scaled_expr = evaluator.eval(scaled[order])
            incoming_expr = evaluator.eval(incoming[order])
            n_expr = evaluator.eval(n_node)
            target_values = [evaluator.eval(targets[index]) for index in range(order + 1)]
            amplitude_values = [evaluator.eval(amplitude[index]) for index in range(order + 1)]
            product_expr = _ordinary_product(amplitude_values, target_values, order)
            expected_expr = n_expr * incoming_expr
            rows.append({
                "order": order,
                "scaled_incoming_node": _handle(scaled[order]),
                "target_M_node": _handle(targets[order]),
                "passed": _sp.simplify(scaled_expr - expected_expr) == 0 and _sp.simplify(product_expr - expected_expr) == 0,
                "N_times_incoming": _expr_string(expected_expr),
                "scaled_incoming": _expr_string(scaled_expr),
                "ordinary_product_A_M": _expr_string(product_expr),
                "ordinary_product_coefficients": [math.comb(order, j) for j in range(order + 1)],
            })
    except (_UnsupportedGraphOperation, TypeError, ValueError) as exc:
        return {"passed": False, "orders": rows, "reason": str(exc), "bounded_operations_only": True}
    return {
        "passed": all(row["passed"] for row in rows),
        "orders": rows,
        "bounded_operations_only": True,
        "opaque_atoms": {"N": n_node, "A": [_handle(value) for value in amplitude[:4]], "incoming_m": [_handle(value) for value in incoming[:4]]},
        "ordinary_Z_scope": [0, 1, 2, 3],
    }


def _identity_residual(nodes: Any, compact: Mapping[str, Any]) -> dict[str, Any]:
    targets = list(_mapping(compact.get("targets")).get("M", ()))
    residual = list(_mapping(compact.get("residual")).get("M", ()))
    h = compact.get("h")
    n_node = _handle(compact.get("N"))
    b_matrix = _mapping(compact.get("W")).get("B")
    b_row = b_matrix[0] if isinstance(b_matrix, (list, tuple)) and b_matrix else ()
    required = len(targets) >= 4 and len(residual) >= 4 and isinstance(h, list) and len(h) == 5 and all(len(row) >= 4 for row in h) and n_node is not None
    if not required:
        return {"passed": False, "orders": [], "reason": "missing actual residual, target M or h rows"}
    opaque = {n_node: _symbol("N")}
    opaque.update({_handle(value): _symbol(f"T{index}") for index, value in enumerate(targets) if _handle(value) is not None})
    for component, row in enumerate(h):
        opaque.update({_handle(value): _symbol(f"h{component}_{index}") for index, value in enumerate(row) if _handle(value) is not None})
    evaluator = _SymbolicGraph(nodes, opaque)
    rows: list[dict[str, Any]] = []
    try:
        for order in range(4):
            lhs = evaluator.eval(residual[order])
            rhs = evaluator.eval(targets[order]) + evaluator.eval(h[0][order]) + evaluator.eval(h[1][order])
            rows.append({
                "order": order,
                "residual_M_node": _handle(residual[order]),
                "passed": _sp.simplify(lhs - rhs) == 0,
                "residual_M": _expr_string(lhs),
                "target_plus_h0_plus_h1": _expr_string(rhs),
            })
    except (_UnsupportedGraphOperation, TypeError, ValueError) as exc:
        return {"passed": False, "orders": rows, "reason": str(exc), "bounded_operations_only": True}
    b_ok = list(b_row) == [1, 1, 0, 0, 0]
    return {
        "passed": b_ok and all(row["passed"] for row in rows),
        "orders": rows,
        "B_M_row": list(b_row),
        "B_M_is_one_one": b_ok,
        "nonlinear_M": "exact zero from actual residual graph",
        "nonlinear_M_zero": all(row["passed"] for row in rows),
        "bounded_operations_only": True,
        "ordinary_Z_scope": [0, 1, 2, 3],
    }


def _identity_terminal_rhs(nodes: Any, compact: Mapping[str, Any]) -> dict[str, Any]:
    amplitude = list(compact.get("A", ()))
    residual = list(_mapping(compact.get("residual")).get("M", ()))
    rhs_rows = list(_mapping(compact.get("terminal_rhs")).get("m", ()))
    n_node = _handle(compact.get("N"))
    required = len(amplitude) >= 4 and len(residual) >= 4 and len(rhs_rows) >= 4 and n_node is not None
    if not required:
        return {"passed": False, "orders": [], "reason": "missing actual terminal m rows"}
    opaque = {n_node: _symbol("N")}
    opaque.update({_handle(value): _symbol(f"A{index}") for index, value in enumerate(amplitude) if _handle(value) is not None})
    opaque.update({_handle(value): _symbol(f"R{index}") for index, value in enumerate(residual) if _handle(value) is not None})
    evaluator = _SymbolicGraph(nodes, opaque)
    rows: list[dict[str, Any]] = []
    try:
        n_expr = evaluator.eval(n_node)
        amplitude_expr = [evaluator.eval(value) for value in amplitude]
        residual_expr = [evaluator.eval(value) for value in residual]
        for order in range(4):
            lhs = evaluator.eval(rhs_rows[order])
            product_expr = _ordinary_product(amplitude_expr, residual_expr, order)
            expected = product_expr / (2 * n_expr)
            rows.append({
                "order": order,
                "terminal_rhs_m_node": _handle(rhs_rows[order]),
                "passed": _sp.simplify(lhs - expected) == 0,
                "terminal_rhs_m": _expr_string(lhs),
                "one_over_2N_ordinary_product": _expr_string(expected),
                "ordinary_product_coefficients": [math.comb(order, j) for j in range(order + 1)],
            })
    except (_UnsupportedGraphOperation, TypeError, ValueError) as exc:
        return {"passed": False, "orders": rows, "reason": str(exc), "bounded_operations_only": True}
    return {
        "passed": all(row["passed"] for row in rows),
        "orders": rows,
        "exp_minus_log2": "1/2",
        "analytic_constant_identity": "exp(-log(2))=1/2",
        "bounded_operations_only": True,
        "ordinary_Z_scope": [0, 1, 2, 3],
    }


def _identity_fixed_point(nodes: Any, compact: Mapping[str, Any]) -> dict[str, Any]:
    system_handle = _handle(compact.get("system"))
    system = _row(nodes, system_handle)
    implicit = _mapping(compact.get("implicit"))
    jacobian = implicit.get("Jacobian")
    recurrence = implicit.get("recurrence")
    targets = list(_mapping(compact.get("targets")).get("M", ()))
    h = compact.get("h")
    system_ok = bool(system and system.get("operation") == "original_compact_unique_Banach_limit_system")
    b_matrix = _mapping(compact.get("W")).get("B")
    b_row = list(b_matrix[0]) if isinstance(b_matrix, (list, tuple)) and b_matrix else []
    b_ok = b_row == [1, 1, 0, 0, 0]
    values_ok = True
    value_rows: list[dict[str, Any]] = []
    if not isinstance(h, list) or len(h) != 5:
        values_ok = False
    else:
        for component, row in enumerate(h):
            value = _handle(row[0]) if row else None
            item = _row(nodes, value)
            passed = bool(
                item
                and item.get("operation") == "original_compact_unique_Banach_limit_component"
                and item.get("system") == system_handle
                and item.get("component") == component
                and item.get("mathematical_function_not_numeric_interval") is True
            )
            value_rows.append({"component": component, "node": value, "passed": passed})
            values_ok = values_ok and passed
    derivative_rows: list[dict[str, Any]] = []
    rhs_rows: list[dict[str, Any]] = []
    recurrence_ok = isinstance(recurrence, list) and len(recurrence) >= 5
    if recurrence_ok and isinstance(h, list):
        for order in range(1, 6):
            rec = _mapping(recurrence[order - 1])
            solved = list(rec.get("solved_ordinary_row", ()))
            known = list(rec.get("known_rhs", ()))
            derivative_passed = len(solved) == 5 and all(_handle(h[component][order]) == _handle(solved[component]) for component in range(5))
            target_passed = len(targets) > order and bool(known) and _handle(known[0]) == _handle(targets[order])
            full_jacobian = rec.get("same_full_Jacobian") is True
            derivative_rows.append({"order": order, "passed": derivative_passed, "solved_ordinary_row": solved})
            rhs_rows.append({"order": order, "passed": target_passed, "known_rhs_0": _handle(known[0]) if known else None, "target_M": _handle(targets[order]) if len(targets) > order else None, "same_full_Jacobian": full_jacobian})
            recurrence_ok = recurrence_ok and derivative_passed and target_passed and full_jacobian
    branch_transfer = _mapping(compact.get("branch_transfer"))
    c1 = _mapping(branch_transfer.get("original_C1_limit_proof"))
    source_correspondence = _mapping(branch_transfer.get("exact_source_correspondence"))
    unique_branch = c1.get("exact_C1_limit_exists_unique_for_this_repair_only_integer") is True
    jift = "IFT" in str(branch_transfer.get("C5_extension", "")) and "same invertible J" in str(branch_transfer.get("C5_extension", ""))
    no_finite_iterate = "no finite iterate" in str(system.get("definition", "")).lower() and "h*=lim_C1" in str(c1.get("exact_definition", ""))
    no_cap_as_value = all(_mapping(value).get("source_value_or_exponential_not_materialized") is True for value in _mapping(_mapping(c1.get("whole_Z_target_cap")).get("transformed_N_scaled_target_C1_caps")).values())
    source_branch = bool(
        branch_transfer.get("same_original_B_Q_N_and_exact_bump_integrals") is True
        and source_correspondence.get("physical_selected_N_incoming_not_divided_again") is True
        and source_correspondence.get("original_positive_Rc_amplitude_and_joint_mu_target_normalization") is True
    )
    passed = bool(
        system_ok and b_ok and values_ok and recurrence_ok and unique_branch and jift
        and no_finite_iterate and no_cap_as_value and source_branch
    )
    return {
        "passed": passed,
        "system_node": system_handle,
        "system_operation": system.get("operation") if system else None,
        "N_node": _handle(compact.get("N")),
        "W_B_M_row": b_row,
        "W_B_M_row_is_one_one": b_ok,
        "h_value_leaves": value_rows,
        "h_derivative_rows": derivative_rows,
        "known_rhs_target_M_rows": rhs_rows,
        "same_full_Jacobian": all(row.get("same_full_Jacobian") is True for row in rhs_rows),
        "unique_Banach_branch": unique_branch,
        "C5_IFT_extension": jift,
        "zero_residual_from_unique_branch": unique_branch and jift,
        "no_finite_iterate_value": no_finite_iterate,
        "no_numeric_cap_as_value": no_cap_as_value,
        "source_branch_correspondence": source_branch,
    }


def _terminal_identity(nodes: Any, compact: Mapping[str, Any], receipt_guard: Mapping[str, Any]) -> dict[str, Any]:
    terminal = _mapping(compact.get("terminal"))
    correction = _mapping(terminal.get("correction"))
    endpoint_rows = list(correction.get("m", ()))
    rhs_rows = list(_mapping(compact.get("terminal_rhs")).get("m", ()))
    certificate_handle = _handle(_mapping(compact.get("terminal_certificates")).get("m"))
    certificate = _row(nodes, certificate_handle)
    expected_recurrence = _mapping(compact.get("implicit")).get("recurrence")
    binding = _mapping(_mapping(compact.get("bindings")).get("original_partial_history_terminal_identity"))
    cert_endpoint = list(certificate.get("actual_endpoint_rows", ())) if certificate else []
    cert_rhs = list(certificate.get("equivalent_residual_rows", ())) if certificate else []
    zeros = list(certificate.get("zero_rows", ())) if certificate else []
    endpoint_ok = bool(endpoint_rows and cert_endpoint == endpoint_rows and all(_row(nodes, value) is not None for value in endpoint_rows))
    rhs_ok = bool(rhs_rows and cert_rhs == rhs_rows and all(_row(nodes, value) is not None for value in rhs_rows))
    recurrence_ok = bool(certificate and certificate.get("ordinary_implicit_recurrence") == expected_recurrence)
    source_system_ok = bool(certificate and certificate.get("source_system") == _handle(compact.get("system")))
    binding_ok = bool(certificate and certificate.get("original_signed_endpoint_binding") == binding and binding.get("function") == "exact_partial_band_graph")
    zero_ok = len(zeros) >= 4 and all(_zero_row(nodes, value) for value in zeros[:4])
    all_zero_rows = len(zeros) >= 6 and all(_zero_row(nodes, value) for value in zeros[:6])
    relative_only = bool(certificate and certificate.get("relative_correction_only") is True)
    absolute_retained = bool(certificate and certificate.get("absolute_leading_histories_and_pressure_not_zeroed") is True)
    theorem_receipt = receipt_guard.get("pinned_original_bump_theorem") is True
    passed = all((endpoint_ok, rhs_ok, recurrence_ok, source_system_ok, binding_ok, zero_ok, relative_only, absolute_retained, theorem_receipt))
    return {
        "passed": passed,
        "certificate_node": certificate_handle,
        "certificate": {
            "operation": certificate.get("operation") if certificate else None,
            "actual_endpoint_rows": cert_endpoint,
            "equivalent_residual_rows": cert_rhs,
            "source_system": certificate.get("source_system") if certificate else None,
            "zero_rows": zeros,
            "relative_correction_only": relative_only,
            "absolute_leading_histories_and_pressure_not_zeroed": absolute_retained,
            "ordinary_implicit_recurrence": expected_recurrence,
            "original_signed_endpoint_binding": binding,
            "theorem": certificate.get("theorem") if certificate else None,
        },
        "actual_endpoint_rows": [_handle(value) for value in endpoint_rows],
        "equivalent_residual_rows": [_handle(value) for value in rhs_rows],
        "certificate_endpoint_rows_match": endpoint_ok,
        "certificate_rhs_rows_match": rhs_ok,
        "certificate_recurrence_matches": recurrence_ok,
        "certificate_source_system_matches": source_system_ok,
        "original_bump_binding_matches": binding_ok,
        "zero_rows": [_handle(value) for value in zeros],
        "zero_rows_Z0_to_Z3": zero_ok,
        "all_six_zero_rows": all_zero_rows,
        "relative_correction_only": relative_only,
        "absolute_history_not_zeroed": absolute_retained,
        "accepted_source_receipt_theorem": theorem_receipt,
        "terminal_relative_m_zero": zero_ok,
        "ordinary_Z_scope": [0, 1, 2, 3],
        "source_values_selected": False,
    }


def build_tail_mass_identity(nodes: Any, compact_archive: Mapping[str, Any] | None, *, root: Any = ROOT) -> dict[str, Any]:
    """Build the source-only compact-tail mass identity certificate."""
    repo = _repo_root(root)
    graph = _nodes(nodes)
    loaded = _load_gzip(repo, _COMPACT_ARCHIVE)
    packet = _find_packet(compact_archive, "true_compact_C5_bindings") or loaded
    compact = _compact_bindings(packet)
    receipt = _load_json(repo, _COMPACT_RECEIPT)
    pins = {name: _pin(repo, path, expected) for name, (path, expected) in _SOURCE_PINS.items()}
    source_guard = _source_guard(repo, compact)
    receipt_guard = _receipt_guard(receipt)
    family_ok = packet.get("source_family") == _FAMILY and receipt.get("source_family") == _FAMILY
    immutable = bool(
        packet.get("true_compact_C5_bindings") == loaded.get("true_compact_C5_bindings")
        and packet.get("exact_graph_nodes") == loaded.get("exact_graph_nodes")
        and compact
    )
    prefix = _prefix(graph, loaded)
    incoming_target = _identity_incoming_target(graph, compact)
    residual = _identity_residual(graph, compact)
    terminal_rhs = _identity_terminal_rhs(graph, compact)
    fixed_point = _identity_fixed_point(graph, compact)
    terminal = _terminal_identity(graph, compact, receipt_guard)
    required = {
        "source_pins": all(pin.get("passed") is True for pin in pins.values()),
        "source_guard": source_guard.get("passed") is True,
        "receipt": receipt_guard.get("passed") is True,
        "source_family": family_ok,
        "immutable_packet": immutable,
        "graph_prefix": prefix,
        "incoming_target_product": incoming_target.get("passed") is True,
        "residual_M_identity": residual.get("passed") is True,
        "terminal_rhs_identity": terminal_rhs.get("passed") is True,
        "fixed_point_IFT": fixed_point.get("passed") is True,
        "terminal_identity": terminal.get("passed") is True,
    }
    unproved = [name for name, value in required.items() if not value]
    mass_passed = all(required.values())
    certificate = {
        "passed": mass_passed,
        "source_only": True,
        "source_values_selected": False,
        "numeric_value_selected": False,
        "numeric_caps_used": False,
        "incoming_target_product": incoming_target,
        "residual_M_identity": residual,
        "terminal_rhs_identity": terminal_rhs,
        "fixed_point_IFT": fixed_point,
        "terminal_identity": terminal,
        "terminal_relative_m_zero": mass_passed and terminal.get("terminal_relative_m_zero") is True,
        "m_relative_zero_source_identity": mass_passed and terminal.get("terminal_relative_m_zero") is True,
        "absolute_modified_m_equals_leading_m": mass_passed and terminal.get("terminal_relative_m_zero") is True,
        "absolute_history_not_zeroed": mass_passed and terminal.get("absolute_history_not_zeroed") is True,
        "absolute_leading_histories_not_zeroed": mass_passed and terminal.get("absolute_history_not_zeroed") is True,
    }
    return {
        "passed": mass_passed,
        "source_only": True,
        "source_values_selected": False,
        "numeric_value_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "producer_or_checker_executed": False,
        "source_pins": pins,
        "sourcepins": pins,
        "source_family": _FAMILY,
        "mass_identity_certificate": certificate,
        "tail_mass_identity_certificate": certificate,
        "certificate_node": terminal.get("certificate_node"),
        "certificate": terminal.get("certificate", {}),
        "ordinary_Z_orders": [0, 1, 2, 3],
        "relative_correction_only": terminal.get("relative_correction_only") is True,
        "terminal_relative_m_zero": mass_passed and terminal.get("terminal_relative_m_zero") is True,
        "terminal_relative_m_zero_Z0_to_Z3": mass_passed and terminal.get("zero_rows_Z0_to_Z3") is True,
        "terminal_relative_m_zero_rows": terminal.get("zero_rows", [])[:4],
        "m_relative_zero_source_identity": mass_passed and terminal.get("terminal_relative_m_zero") is True,
        "absolute_modified_m_equals_leading_m": mass_passed and terminal.get("terminal_relative_m_zero") is True,
        "absolute_history_not_zeroed": mass_passed and terminal.get("absolute_history_not_zeroed") is True,
        "absolute_leading_histories_not_zeroed": mass_passed and terminal.get("absolute_history_not_zeroed") is True,
        "checks": {
            "required": required,
            "source_guard": source_guard,
            "receipt": receipt_guard,
            "immutable_packet": immutable,
            "graph_prefix": prefix,
            "incoming_target_product": incoming_target,
            "residual_M_identity": residual,
            "terminal_rhs_identity": terminal_rhs,
            "fixed_point_IFT": fixed_point,
            "terminal_identity": terminal,
        },
        "scope": {
            "ordinary_Z_orders": [0, 1, 2, 3],
            "common_N_node": _handle(compact.get("N")),
            "amplitude_nodes": [_handle(value) for value in list(compact.get("A", ()))[:4]],
            "incoming_m_nodes": [_handle(value) for value in list(_mapping(compact.get("incoming_correction")).get("m", ()))[:4]],
            "target_M_nodes": [_handle(value) for value in list(_mapping(compact.get("targets")).get("M", ()))[:4]],
        },
        "unproved": unproved,
        "limits": {
            "numeric_terminal_values": False,
            "numeric_caps_as_values": False,
            "absolute_history_zero_claim": False,
        },
    }


__all__ = ["build_tail_mass_identity"]
