"""Source-only binding for the saved lower ``delta`` and ``logPstar`` records.

The leaf accepts decoded mappings and producer text.  It never opens files,
imports a producer, executes source code, or traverses a source graph.  Its
claim is limited to exact saved-record identity plus an AST check that the
producer setup copied the two native records directly.
"""

from __future__ import annotations

import ast
import copy
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "original_lower_parameter_binding_v1"
EXPECTED_NATIVE_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)
EXPECTED_CORRECTED_SCHEMA = "corrected_reserved_m15_inputs_v1"
DELTA_SOURCE_PATH = ("actual_micro_Q_increment_bounds", "actual_delta")
LOGPSTAR_SOURCE_PATH = (
    "actual_micro_background_I1_numeric",
    "physical_scale",
    "logPstar",
)


class LowerParameterBindingError(ValueError):
    """Raised when a saved parameter binding is not exact."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise LowerParameterBindingError(message)


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    _require(isinstance(value, Mapping), label + " must be a mapping")
    return value


def _at(root: Mapping[str, Any], path: Sequence[str], label: str) -> Any:
    value: Any = root
    for key in path:
        value = _mapping(value, label + " parent").get(key)
        _require(value is not None, label + " is missing at " + key)
    return value


def _binary_tuple(value: Any, label: str) -> tuple[int, int, int, int]:
    _require(isinstance(value, (list, tuple)) and len(value) == 4,
             label + " must be a four-field binary mpf tuple")
    values: list[int] = []
    for index, item in enumerate(value):
        _require(type(item) is int, f"{label}[{index}] must be an integer")
        values.append(item)
    return tuple(values)  # type: ignore[return-value]


def _binary_endpoints(value: Any, label: str) -> dict[str, Any]:
    record = _mapping(value, label)
    _require(record.get("directed") is True, label + " is not directed")
    lower = _binary_tuple(record.get("lower_exact_mpf_tuple"), label + ".lower")
    upper = _binary_tuple(record.get("upper_exact_mpf_tuple"), label + ".upper")
    return {
        "lower_exact_mpf_tuple": list(lower),
        "upper_exact_mpf_tuple": list(upper),
        "directed": True,
    }


def _same_record(left: Any, right: Any, label: str) -> None:
    _binary_endpoints(left, label + ".left")
    _binary_endpoints(right, label + ".right")
    _require(left == right, label + " differs from native saved record")


def _subscript_key(node: ast.Subscript) -> str | None:
    value: Any = node.slice
    if isinstance(value, ast.Index):  # pragma: no cover - old Python AST
        value = value.value
    if isinstance(value, ast.Constant) and type(value.value) is str:
        return value.value
    return None


def _is_deepcopy_path(node: ast.AST, path: Sequence[str]) -> bool:
    """Check exactly ``copy.deepcopy(native[...])`` without evaluating AST."""

    if not isinstance(node, ast.Call) or node.keywords or len(node.args) != 1:
        return False
    function = node.func
    if not (
        isinstance(function, ast.Attribute)
        and function.attr == "deepcopy"
        and isinstance(function.value, ast.Name)
        and function.value.id == "copy"
    ):
        return False
    current: ast.AST = node.args[0]
    for key in reversed(tuple(path)):
        if not isinstance(current, ast.Subscript) or _subscript_key(current) != key:
            return False
        current = current.value
    return isinstance(current, ast.Name) and current.id == "native"


def _static_string(node: ast.AST, assignments: Mapping[str, ast.AST], seen: set[str] | None = None) -> str | None:
    """Resolve only literal string concatenation in the producer AST."""

    seen = set() if seen is None else seen
    if isinstance(node, ast.Constant) and type(node.value) is str:
        return node.value
    if isinstance(node, ast.Name) and node.id not in seen:
        value = assignments.get(node.id)
        if value is not None:
            return _static_string(value, assignments, seen | {node.id})
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        left = _static_string(node.left, assignments, seen)
        right = _static_string(node.right, assignments, seen)
        if left is not None and right is not None:
            return left + right
    return None


def _name(node: ast.AST, expected: str) -> bool:
    return isinstance(node, ast.Name) and node.id == expected


def _is_native_read_pinned(node: ast.AST) -> bool:
    """Check exactly ``native = read_pinned(BASE / NATIVE, NATIVE_HASH)``."""

    if not isinstance(node, ast.Call) or node.keywords or len(node.args) != 2:
        return False
    if not _name(node.func, "read_pinned") or not _name(node.args[1], "NATIVE_HASH"):
        return False
    path = node.args[0]
    return (
        isinstance(path, ast.BinOp)
        and isinstance(path.op, ast.Div)
        and _name(path.left, "BASE")
        and _name(path.right, "NATIVE")
    )


def _source_ast_contract(source_text: str, corrected: Mapping[str, Any]) -> dict[str, Any]:
    _require(isinstance(source_text, str), "producer source must be text")
    try:
        tree = ast.parse(source_text, filename="run_original_corrected_reserved_m15.py")
    except SyntaxError as exc:
        raise LowerParameterBindingError("producer source is not parseable") from exc

    assignments: dict[str, ast.AST] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            value = node.value
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                if isinstance(target, ast.Name):
                    assignments[target.id] = value
    native_hash = _static_string(assignments.get("NATIVE_HASH", ast.Constant(None)), assignments)
    _require(native_hash == EXPECTED_NATIVE_SHA256, "producer NATIVE_HASH changed")
    native_name = _static_string(assignments.get("NATIVE", ast.Constant(None)), assignments)
    _require(native_name is not None and native_name.endswith("native_directed_controls_v2.json"),
             "producer NATIVE report path changed")
    parent_name = _static_string(assignments.get("PARENT", ast.Constant(None)), assignments)
    _require(parent_name is not None and parent_name.endswith("outer_pressure_enclosure.json"),
             "producer PARENT report path changed")
    parent_hash = _static_string(assignments.get("PARENT_HASH", ast.Constant(None)), assignments)
    if "parent_sha256" in corrected:
        _require(parent_hash is not None and corrected.get("parent_sha256") == parent_hash,
                 "producer PARENT_HASH does not match corrected report")

    setup_nodes = [node for node in tree.body
                   if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                   and node.name == "setup"]
    _require(len(setup_nodes) == 1, "producer setup function is missing or duplicated")
    source_call: ast.Call | None = None
    for node in ast.walk(setup_nodes[0]):
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "source"
                   for target in node.targets):
            continue
        if isinstance(node.value, ast.Call) and _name(node.value.func, "dict"):
            source_call = node.value
            break
    _require(source_call is not None, "setup source dict assignment is missing")
    native_read_pinned = False
    for node in ast.walk(setup_nodes[0]):
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "native"
                   for target in node.targets):
            continue
        if _is_native_read_pinned(node.value):
            native_read_pinned = True
            break
    _require(native_read_pinned,
             "setup native is not read_pinned(BASE / NATIVE, NATIVE_HASH)")
    keywords = {item.arg: item.value for item in source_call.keywords if item.arg is not None}
    _require(_name(keywords.get("native_sha256", ast.Constant(None)), "NATIVE_HASH"),
             "setup native_sha256 is not bound to NATIVE_HASH")
    if "parent_sha256" in corrected:
        _require(_name(keywords.get("parent_sha256", ast.Constant(None)), "PARENT_HASH"),
                 "setup parent_sha256 is not bound to PARENT_HASH")
    _require(_is_deepcopy_path(keywords.get("delta", ast.Constant(None)), DELTA_SOURCE_PATH),
             "setup delta is not direct native actual_delta deepcopy")
    _require(_is_deepcopy_path(
        keywords.get("actual_logPstar_parameter", ast.Constant(None)),
        LOGPSTAR_SOURCE_PATH,
    ), "setup logPstar is not direct native physical_scale.logPstar deepcopy")
    return {
        "parsed": True,
        "setup_found": True,
        "source_dict_found": True,
        "native_read_pinned_assignment": True,
        "native_hash_literal": native_hash,
        "native_report_name": native_name,
        "parent_report_name": parent_name,
        "parent_hash_literal": parent_hash,
        "delta_direct_deepcopy": True,
        "actual_logPstar_parameter_direct_deepcopy": True,
        "producer_ast_scope": "setup assignment only; source not executed",
    }


def _family(native: Mapping[str, Any]) -> Any:
    value = native.get("source_family")
    _require(value is not None, "native source_family is missing")
    for section_name in (DELTA_SOURCE_PATH[0], LOGPSTAR_SOURCE_PATH[0]):
        section = _mapping(native.get(section_name), "native." + section_name)
        if "source_family" in section:
            _require(section.get("source_family") == value,
                     "native parameter section source_family changed")
    return value


def _ledger_scope(corrected: Mapping[str, Any], terminal: Mapping[str, Any]) -> dict[str, Any]:
    _require(terminal.get("old_graphs_unchanged") is True,
             "corrected old graph claim changed")
    _require(terminal.get("numerical_ancestor_executed") is False,
             "corrected numerical ancestor scope changed")
    _require(corrected.get("old_212_terms_and_all_errors_unchanged") is True,
             "corrected old ledger preservation claim changed")
    source_result = _mapping(corrected.get("corrected_reserved_m15_source"),
                             "corrected source result")
    _require(source_result.get("source_only") is True,
             "corrected source is not source-only")
    _require(source_result.get("source_executed") is False,
             "corrected source execution scope changed")
    _require(source_result.get("existing_ledger_edited") is False,
             "corrected existing ledger was marked edited")
    return {
        "old_graphs_unchanged": True,
        "old_212_terms_and_all_errors_unchanged": True,
        "source_only": True,
        "source_executed": False,
        "existing_ledger_edited": False,
    }


def _repair_contract(receipt: Mapping[str, Any], family: Any, delta: Any, logpstar: Any) -> dict[str, Any]:
    outer = _mapping(receipt.get("result"), "repair receipt.result")
    result = _mapping(outer.get("result"), "repair receipt.result.result")
    _require(result.get("source_family") == family, "repair receipt source_family changed")
    contract = _mapping(result.get("corrected_terminal_source_contract"),
                        "repair terminal source contract")
    _same_record(contract.get("actual_delta"), delta, "repair receipt actual_delta")
    _same_record(contract.get("actual_logPstar_parameter"), logpstar,
                 "repair receipt actual_logPstar_parameter")
    _require(contract.get("source_only") is True, "repair contract is not source-only")
    _require(result.get("ledger_written") is False, "repair receipt ledger_written changed")
    _require(result.get("old_errors_removed") is False, "repair receipt old_errors_removed changed")
    _require(result.get("source_rows_are_not_point_values") is True,
             "repair receipt promoted source rows to values")
    _require(result.get("endpoint_function_values_selected") is False,
             "repair receipt selected endpoint function values")
    leaf = result.get("leaf")
    if isinstance(leaf, Mapping):
        _require(leaf.get("ordinary_Z1_is_product_rule") is True,
                 "repair receipt ordinary Z1 product rule missing")
        _require(leaf.get("ordinary_Z1_terms_not_merged") is True,
                 "repair receipt ordinary Z1 terms were merged")
    return {
        "path": ["result", "result", "corrected_terminal_source_contract"],
        "source_family_matches": True,
        "actual_delta_matches": True,
        "actual_logPstar_parameter_matches": True,
        "ledger_written": False,
        "old_errors_removed": False,
        "source_rows_are_not_point_values": True,
        "endpoint_function_values_selected": False,
        "ordinary_Z1_flags_checked": isinstance(leaf, Mapping),
    }


def bind(native: Mapping[str, Any], corrected: Mapping[str, Any],
         repair_receipt: Mapping[str, Any], producer_source_text: str) -> dict[str, Any]:
    """Validate exact saved parameter identity without touching real sources."""

    native = _mapping(native, "native_directed_controls_v2")
    corrected = _mapping(corrected, "corrected_reserved_m15_source")
    repair_receipt = _mapping(repair_receipt, "four-cell repair receipt")
    terminal = _mapping(corrected.get("corrected_terminal_source_inputs"),
                        "corrected terminal source inputs")
    _require(terminal.get("schema") == EXPECTED_CORRECTED_SCHEMA,
             "corrected terminal schema changed")
    _require(terminal.get("source_only") is True, "corrected terminal is not source-only")
    _require(terminal.get("native_sha256") == EXPECTED_NATIVE_SHA256,
             "corrected terminal native_sha256 changed")
    _require(terminal.get("numeric_endpoint_cell") == ["0", ".5"],
             "corrected terminal endpoint cell changed")

    family = _family(native)
    native_delta = _at(native, DELTA_SOURCE_PATH, "native delta")
    native_logpstar = _at(native, LOGPSTAR_SOURCE_PATH, "native logPstar")
    corrected_delta = terminal.get("delta")
    corrected_logpstar = terminal.get("actual_logPstar_parameter")
    _same_record(corrected_delta, native_delta, "corrected terminal delta")
    _same_record(corrected_logpstar, native_logpstar,
                 "corrected terminal actual_logPstar_parameter")

    proof = _mapping(terminal.get("physical_endpoint_source_proof"),
                     "corrected endpoint proof")
    _require(proof.get("source_family") == family,
             "corrected terminal source_family changed")
    _require(proof.get("actual_V_y_rows_0_through_3_factorial_converted_once") is True,
             "corrected V_y factorial provenance changed")
    _require(proof.get("common_Pstar_applied_once") is True,
             "corrected common Pstar provenance changed")
    rows = _mapping(terminal.get("physical_endpoint_source_rows"),
                    "corrected endpoint rows")
    _require(isinstance(rows.get("R"), Mapping) and isinstance(rows.get("Pstar"), Mapping),
             "corrected R/Pstar rows are missing")

    ledger = _ledger_scope(corrected, terminal)
    ast_contract = _source_ast_contract(producer_source_text, terminal)
    receipt_contract = _repair_contract(repair_receipt, family,
                                        corrected_delta, corrected_logpstar)
    return {
        "schema": SCHEMA,
        "passed": True,
        "source_family": copy.deepcopy(family),
        "native_sha256": EXPECTED_NATIVE_SHA256,
        "corrected_native_sha256_matches": True,
        "delta": {
            "source_path": list(DELTA_SOURCE_PATH),
            "native_record": copy.deepcopy(native_delta),
            "corrected_record": copy.deepcopy(corrected_delta),
            "binary_endpoints": _binary_endpoints(corrected_delta, "corrected delta"),
            "exact_saved_record_match": True,
        },
        "actual_logPstar_parameter": {
            "source_path": list(LOGPSTAR_SOURCE_PATH),
            "native_record": copy.deepcopy(native_logpstar),
            "corrected_record": copy.deepcopy(corrected_logpstar),
            "binary_endpoints": _binary_endpoints(corrected_logpstar,
                                                   "corrected logPstar"),
            "exact_saved_record_match": True,
        },
        "repair_receipt_contract": receipt_contract,
        "source_ast_contract": ast_contract,
        "ledger_scope": ledger,
        "ordinary_Z1_scope": {
            "source_rows_retained": True,
            "factorial_conversion_checked_once": True,
            "second_derivative_or_second_factorial_conversion": False,
            "product_rule_rederived_here": False,
        },
        "provenance_scope": {
            "archived_parameter_identity_only": True,
            "endpoint_live_owner_identity": False,
            "delta_min_branch_identity": False,
            "Pstar_source_function_identity": False,
            "cell_basis_identity": False,
            "native_corrected_endpoint_identity": False,
            "interval_equality_promoted_to_function_identity": False,
        },
        "source_execution": False,
        "real_source_imported": False,
        "receipt_decoded_here": False,
        "graph_traversed": False,
    }


def build(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Builder-shaped alias for frozen callers."""

    return bind(*args, **kwargs)


__all__ = [
    "SCHEMA",
    "EXPECTED_NATIVE_SHA256",
    "DELTA_SOURCE_PATH",
    "LOGPSTAR_SOURCE_PATH",
    "LowerParameterBindingError",
    "bind",
    "build",
]
