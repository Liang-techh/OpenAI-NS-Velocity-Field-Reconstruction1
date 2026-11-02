"""Conditional source-string AST identity for the two low core matrix leaves.

The runner supplies the source text for the source-matrix leaf and the
matrix-only leaf.  This module deliberately does not import either producer,
read a file, evaluate an expression, or perform numerical algebra.  It parses
the two strings, extracts the actual ``c0_jets``/``b1_jets`` dictionaries and
the seven named intermediate assignments, then compares a small, explicit
AST language after the documented symbol correspondence.

The result is a conditional source-formula receipt.  It is not a numerical
identity, an automatic-differentiation-jet identity, a physical-owner claim,
or a proof of F0 relative-jet equivalence.
"""

from __future__ import annotations

import ast
from typing import Any


SCHEMA = "core_operator_formula_identity_v1"
HIGH_LABEL = "n1_core_low_source_matrix_v2.py"
LOW_LABEL = "n1_core_low_matrix_v3.py"

EXPECTED_C0_SUPPORT = (
    "4,1",
    "5,1",
    "5,2",
    "5,3",
    "5,5",
    "6,1",
    "6,2",
    "6,3",
    "6,4",
    "6,6",
)
EXPECTED_B1_SUPPORT = (
    "5,1",
    "5,2",
    "5,3",
    "6,2",
    "6,3",
    "6,4",
)
EXPECTED_SUPPORT = {
    "C0": EXPECTED_C0_SUPPORT,
    "B1": EXPECTED_B1_SUPPORT,
}
EXPECTED_MATRIX_ITEM_COUNT = len(EXPECTED_C0_SUPPORT) + len(EXPECTED_B1_SUPPORT)

EXPECTED_POWERS = {
    "C0": {
        "4,1": (1, -1),
        "5,1": (0, -1),
        "5,2": (1, -1),
        "5,3": (1, -1),
        "5,5": (0, -1),
        "6,1": (1, -2),
        "6,2": (0, -1),
        "6,3": (0, -1),
        "6,4": (0, -1),
        "6,6": (0, -1),
    },
    "B1": {
        "5,1": (0, -1),
        "5,2": (1, -1),
        "5,3": (1, -1),
        "6,2": (0, -1),
        "6,3": (0, -1),
        "6,4": (0, -1),
    },
}

INTERMEDIATE_NAMES = ("Q", "H", "B", "E", "Cf", "Zf", "Zw")

# The direction is high/source-matrix spelling -> low/matrix spelling.  The
# entries for the sqrt and mpf calls are checked by static definition anchors,
# not by replacing arbitrary calls in the source text.
SYMBOL_CORRESPONDENCE = {
    "z_j": "z_dual",
    "rho_j": "rho_dual",
    "delta_j": "delta_dual",
    "L_j": "L",
    "d_j": "d",
    "U0": "w",
    "U1": "wr",
    "M0": "m",
    "M0_Z": "mz",
    "fz": "fz",
    "context": "c",
    "rho_iv": "rho_value",
    "one": "one",
    "x": "x",
    "context.sqrt(rho_iv)": "x",
    "context.mpf('0.5')": "c.mpf('0.5')",
}
LITERAL_CORRESPONDENCE = {"1": "one"}

ASSUMPTIONS = (
    "Both inputs are source strings supplied by the runner.",
    "Only ast.parse and the controlled expression canonicalizer are used.",
    "The named c0_jets and b1_jets assignments are the producer dictionaries.",
    "The expected supports and F0/Lambda powers are fixed by this module.",
    "The high-to-low symbol correspondence is the explicit mapping in this receipt.",
    "The source literal 1 is identified with the one dual only after both one anchors match.",
    "context.mpf('0.5') and c.mpf('0.5') are identified only after the E anchors match.",
    "context.sqrt(rho_iv) and x are identified only after the static sqrt anchors match.",
)

UNPROVED = (
    "numeric equality is not claimed",
    "automatic differentiation jet equivalence is not claimed",
    "physical owner is not claimed",
    "F0 relative-jet equivalence is not claimed",
)


class SourceIdentityError(ValueError):
    """Reserved for callers that want to classify malformed source input."""


class _NormalizationError(ValueError):
    """An expression contains syntax outside the admitted AST language."""


def _ast_dump(node: ast.AST) -> str:
    """Return the original, location-free AST expression receipt."""

    return ast.dump(node, annotate_fields=True, include_attributes=False)


def _dotted(node: ast.AST) -> str | None:
    """Return a dotted Name/Attribute spelling without evaluating it."""

    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        prefix = _dotted(node.value)
        if prefix is None:
            return None
        return prefix + "." + node.attr
    return None


def _is_name(node: ast.AST, name: str) -> bool:
    return isinstance(node, ast.Name) and node.id == name


def _is_int(node: ast.AST, value: int) -> bool:
    return isinstance(node, ast.Constant) and type(node.value) is int and node.value == value


def _is_string(node: ast.AST, value: str) -> bool:
    return isinstance(node, ast.Constant) and type(node.value) is str and node.value == value


def _is_call(node: ast.AST, function: str, arguments: tuple[bool, ...]) -> bool:
    if not isinstance(node, ast.Call) or _dotted(node.func) != function:
        return False
    if node.keywords or len(node.args) != len(arguments):
        return False
    return all(arguments[index] for index, _ in enumerate(node.args))


def _parse_source(source: Any, label: str) -> tuple[ast.Module | None, dict[str, Any] | None]:
    if not isinstance(source, str):
        return None, {
            "kind": "input_type",
            "label": label,
            "message": "source input must be a string",
        }
    try:
        tree = ast.parse(source, filename=label, mode="exec")
    except SyntaxError as exc:
        return None, {
            "kind": "syntax",
            "label": label,
            "message": exc.msg,
            "line": exc.lineno,
            "offset": exc.offset,
        }
    return tree, None


def _assignments(tree: ast.AST | None, target: str) -> list[ast.Assign]:
    if tree is None:
        return []
    found: list[ast.Assign] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        if _is_name(node.targets[0], target):
            found.append(node)
    return found


def _single_assignment(
    tree: ast.AST | None,
    target: str,
    label: str,
) -> tuple[ast.Assign | None, dict[str, Any] | None]:
    found = _assignments(tree, target)
    if len(found) == 1:
        return found[0], None
    if not found:
        return None, {"kind": "missing_assignment", "path": label + "." + target}
    return None, {
        "kind": "ambiguous_assignment",
        "path": label + "." + target,
        "count": len(found),
    }


def _enclosing_function(tree: ast.Module | None, assignment: ast.Assign | None) -> ast.AST | None:
    """Find the function whose body owns an extracted assignment."""

    if tree is None or assignment is None:
        return None
    functions = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and any(item is assignment for item in ast.walk(node))
    ]
    if len(functions) != 1:
        return None
    return functions[0]


def _assignment_receipt(node: ast.Assign | None) -> dict[str, Any]:
    if node is None:
        return {"present": False, "assignment_ast": None, "expression_ast": None, "line": None}
    return {
        "present": True,
        "assignment_ast": _ast_dump(node),
        "expression_ast": _ast_dump(node.value),
        "line": getattr(node, "lineno", None),
    }


def _call_matches(node: ast.AST, function: str, argument_names: tuple[str, ...]) -> bool:
    if not isinstance(node, ast.Call) or _dotted(node.func) != function:
        return False
    if node.keywords or len(node.args) != len(argument_names):
        return False
    return all(_is_name(argument, name) for argument, name in zip(node.args, argument_names))


def _constant_call_matches(node: ast.AST, function: str, value: str) -> bool:
    if not isinstance(node, ast.Call) or _dotted(node.func) != function:
        return False
    return not node.keywords and len(node.args) == 1 and _is_string(node.args[0], value)


def _one_anchor(node: ast.Assign | None, side: str) -> tuple[bool, str]:
    if node is None:
        return False, "missing one assignment"
    if side == "high":
        value = node.value
        if not isinstance(value, ast.Call) or _dotted(value.func) != "IntervalTaylorJet.constant":
            return False, "one is not an IntervalTaylorJet.constant call"
        if value.keywords or len(value.args) != 3:
            return False, "high one call has the wrong argument shape"
        if not _is_name(value.args[0], "context") or not _is_int(value.args[1], 1):
            return False, "high one call does not anchor context and literal 1"
        if not _is_name(value.args[2], "degree"):
            return False, "high one call does not use degree"
        return True, "high one assignment anchors context, literal 1, and degree"
    value = node.value
    if not isinstance(value, ast.Call) or _dotted(value.func) != "Dual.constant":
        return False, "one is not a Dual.constant call"
    if value.keywords or len(value.args) != 2:
        return False, "low one call has the wrong argument shape"
    if not _is_name(value.args[0], "c") or not _is_int(value.args[1], 1):
        return False, "low one call does not anchor c and literal 1"
    return True, "low one assignment anchors c and literal 1"


def _half_anchor(node: ast.Assign | None, side: str) -> tuple[bool, str]:
    if node is None:
        return False, "missing E assignment"
    expected = "context.mpf" if side == "high" else "c.mpf"
    if _constant_call_matches(node.value, expected, "0.5"):
        return True, expected + "('0.5') is the E constant"
    return False, "E does not contain the anchored " + expected + "('0.5') call"


def _matching_calls(node: ast.AST | None, function: str, arguments: tuple[str, ...]) -> list[ast.Call]:
    if node is None:
        return []
    return [
        item
        for item in ast.walk(node)
        if isinstance(item, ast.Call) and _call_matches(item, function, arguments)
    ]


def _sqrt_anchor(
    high_dicts: dict[str, dict[str, Any]],
    low_dicts: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    assumptions = list(ASSUMPTIONS) + [
        "The high sqrt call is checked in C0 5,5 and C0 6,6, while low x is checked as an assignment."
    ]
    high_calls: dict[str, list[ast.Call]] = {}
    for key in ("5,5", "6,6"):
        entry = high_dicts.get("C0", {}).get("entries", {}).get(key)
        node = entry.get("expression_node") if entry else None
        high_calls[key] = _matching_calls(node, "context.sqrt", ("rho_iv",))
    low_scope = low_dicts.get("C0", {}).get("scope_node")
    low_x, low_issue = _single_assignment(low_scope, "x", LOW_LABEL)
    low_x_match = bool(
        low_x is not None and _call_matches(low_x.value, "c.sqrt", ("rho_value",))
    )
    high_ok = all(len(calls) == 1 for calls in high_calls.values())
    low_ok = low_x_match
    low_rows_use_x = True
    for key in ("5,5", "6,6"):
        entry = low_dicts.get("C0", {}).get("entries", {}).get(key)
        node = entry.get("expression_node") if entry else None
        if node is None or not any(_is_name(item, "x") for item in ast.walk(node)):
            low_rows_use_x = False
    matched = high_ok and low_ok and low_rows_use_x
    differences: list[dict[str, Any]] = []
    for key, calls in high_calls.items():
        if len(calls) != 1:
            differences.append(
                {
                    "kind": "sqrt_anchor",
                    "path": "C0." + key,
                    "high_call_count": len(calls),
                    "expected": "context.sqrt(rho_iv)",
                }
            )
    if not low_ok:
        differences.append(
            {
                "kind": "sqrt_anchor",
                "path": "low.x",
                "low": _assignment_receipt(low_x),
                "expected": "x = c.sqrt(rho_value)",
                "issue": low_issue,
            }
        )
    if not low_rows_use_x:
        differences.append(
            {
                "kind": "sqrt_anchor",
                "path": "low.C0.5,5/6,6",
                "expected": "both rows use x",
            }
        )
    return {
        "matched": matched,
        "normalized": {"high": "x", "low": "x"},
        "high": {
            key: {
                "call_count": len(calls),
                "calls_ast": [_ast_dump(call) for call in calls],
            }
            for key, calls in high_calls.items()
        },
        "low": {
            "x_assignment": _assignment_receipt(low_x),
            "x_call_ast": _ast_dump(low_x.value) if low_x is not None else None,
        },
        "assumptions": assumptions,
        "differences": differences,
    }


def _anchor_receipt(
    high_scope: ast.AST | None,
    low_scope: ast.AST | None,
    sqrt_anchor: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    differences: list[dict[str, Any]] = []
    one_high, one_high_issue = _single_assignment(high_scope, "one", HIGH_LABEL)
    one_low, one_low_issue = _single_assignment(low_scope, "one", LOW_LABEL)
    half_high, half_high_issue = _single_assignment(high_scope, "E", HIGH_LABEL)
    half_low, half_low_issue = _single_assignment(low_scope, "E", LOW_LABEL)
    one_high_ok, one_high_reason = _one_anchor(one_high, "high")
    one_low_ok, one_low_reason = _one_anchor(one_low, "low")
    half_high_ok, half_high_reason = _half_anchor(half_high, "high")
    half_low_ok, half_low_reason = _half_anchor(half_low, "low")
    records = {
        "one_literal": {
            "matched": one_high_ok and one_low_ok,
            "normalized": {"high": "one", "low": "one"},
            "high": _assignment_receipt(one_high),
            "low": _assignment_receipt(one_low),
            "high_reason": one_high_reason,
            "low_reason": one_low_reason,
            "assumptions": list(ASSUMPTIONS),
        },
        "half_mpf": {
            "matched": half_high_ok and half_low_ok,
            "normalized": {"high": "mpf('0.5')", "low": "mpf('0.5')"},
            "high": _assignment_receipt(half_high),
            "low": _assignment_receipt(half_low),
            "high_reason": half_high_reason,
            "low_reason": half_low_reason,
            "assumptions": list(ASSUMPTIONS),
        },
        "sqrt_x": sqrt_anchor,
    }
    for path, issue in (
        ("high.one", one_high_issue),
        ("low.one", one_low_issue),
        ("high.E", half_high_issue),
        ("low.E", half_low_issue),
    ):
        if issue is not None:
            differences.append(issue)
    if not one_high_ok or not one_low_ok:
        differences.append(
            {
                "kind": "anchor_mismatch",
                "path": "one_literal",
                "high_reason": one_high_reason,
                "low_reason": one_low_reason,
            }
        )
    if not half_high_ok or not half_low_ok:
        differences.append(
            {
                "kind": "anchor_mismatch",
                "path": "half_mpf",
                "high_reason": half_high_reason,
                "low_reason": half_low_reason,
            }
        )
    if not sqrt_anchor["matched"]:
        differences.extend(sqrt_anchor["differences"])
    records["matched"] = all(
        records[name]["matched"] for name in ("one_literal", "half_mpf", "sqrt_x")
    )
    records["assumptions"] = list(ASSUMPTIONS)
    records["differences"] = differences
    return records, differences


def _normal(node: ast.AST, side: str) -> tuple[Any, ...]:
    """Canonicalize only the admitted expression grammar.

    The return value is a tuple tree, rather than executable code.  Keeping the
    tuple tree in the receipt makes every normalization step inspectable.
    """

    if isinstance(node, ast.Name):
        name = node.id
        if side == "high":
            name = SYMBOL_CORRESPONDENCE.get(name, name)
        return ("name", name)
    if isinstance(node, ast.Constant):
        if type(node.value) is int:
            if node.value == 1:
                return ("name", LITERAL_CORRESPONDENCE["1"])
            return ("int", str(node.value))
        if type(node.value) is str:
            return ("str", node.value)
        raise _NormalizationError("only integer and string constants are admitted")
    if isinstance(node, ast.UnaryOp):
        if isinstance(node.op, ast.USub):
            operator = "-"
        elif isinstance(node.op, ast.UAdd):
            operator = "+"
        else:
            raise _NormalizationError("unsupported unary operator")
        return ("unary", operator, _normal(node.operand, side))
    if isinstance(node, ast.BinOp):
        operators = {
            ast.Add: "+",
            ast.Sub: "-",
            ast.Mult: "*",
            ast.Div: "/",
        }
        operator = next((symbol for kind, symbol in operators.items() if isinstance(node.op, kind)), None)
        if operator is None:
            raise _NormalizationError("unsupported binary operator")
        return ("bin", operator, _normal(node.left, side), _normal(node.right, side))
    if isinstance(node, ast.Call):
        function = _dotted(node.func)
        if function is None or node.keywords:
            raise _NormalizationError("unsupported call shape")
        if function == "context.sqrt" and len(node.args) == 1 and _is_name(node.args[0], "rho_iv"):
            return ("name", "x")
        if function == "c.sqrt" and len(node.args) == 1 and _is_name(node.args[0], "rho_value"):
            return ("name", "x")
        if function in ("context.mpf", "c.mpf") and len(node.args) == 1:
            return ("call", "mpf", _normal(node.args[0], side))
        raise _NormalizationError("unsupported call: " + (function or "<unknown>"))
    if isinstance(node, ast.Attribute):
        function = _dotted(node)
        if function is None:
            raise _NormalizationError("unsupported attribute")
        return ("attribute", function)
    raise _NormalizationError("unsupported expression node: " + type(node).__name__)


def _render(value: Any) -> str:
    """Render a tuple-tree normalization without invoking a parser/evaluator."""

    if not isinstance(value, tuple) or not value:
        return repr(value)
    tag = value[0]
    if tag == "name":
        return str(value[1])
    if tag == "int":
        return str(value[1])
    if tag == "str":
        return repr(value[1])
    if tag == "attribute":
        return str(value[1])
    if tag == "call":
        return str(value[1]) + "(" + ", ".join(_render(item) for item in value[2:]) + ")"
    if tag == "unary":
        return "(" + str(value[1]) + _render(value[2]) + ")"
    if tag == "bin":
        return "(" + _render(value[2]) + " " + str(value[1]) + " " + _render(value[3]) + ")"
    return repr(value)


def _formula_receipt(node: ast.AST | None, side: str, assumptions: list[str]) -> dict[str, Any]:
    if node is None:
        return {
            "present": False,
            "ast": None,
            "source_ast": None,
            "normalized": None,
            "normalized_value": None,
            "assumptions": assumptions,
            "error": "missing expression",
        }
    original = _ast_dump(node)
    try:
        normalized_value = _normal(node, side)
    except _NormalizationError as exc:
        return {
            "present": True,
            "ast": original,
            "source_ast": original,
            "normalized": None,
            "normalized_value": None,
            "assumptions": assumptions,
            "error": str(exc),
        }
    return {
        "present": True,
        "ast": original,
        "source_ast": original,
        "normalized": _render(normalized_value),
        "normalized_value": normalized_value,
        "assumptions": assumptions,
        "error": None,
    }


def _integer(node: ast.AST | None) -> int | None:
    if isinstance(node, ast.Constant) and type(node.value) is int:
        return node.value
    return None


def _dict_info(tree: ast.Module | None, name: str, label: str) -> dict[str, Any]:
    node, issue = _single_assignment(tree, name, label)
    info: dict[str, Any] = {
        "present": node is not None,
        "assignment": _assignment_receipt(node),
        "scope_node": _enclosing_function(tree, node),
        "entries": {},
        "actual_support": [],
        "errors": [],
    }
    if issue is not None:
        info["errors"].append(issue)
    if node is None:
        return info
    if not isinstance(node.value, ast.Dict):
        info["errors"].append(
            {"kind": "not_dict", "path": label + "." + name, "ast": _ast_dump(node.value)}
        )
        return info
    for key_node, value_node in zip(node.value.keys, node.value.values):
        key = key_node.value if isinstance(key_node, ast.Constant) and type(key_node.value) is str else None
        if key is None:
            info["errors"].append(
                {"kind": "invalid_key", "path": label + "." + name, "key_ast": _ast_dump(key_node)}
            )
            continue
        if key in info["entries"]:
            info["errors"].append(
                {"kind": "duplicate_key", "path": label + "." + name + "." + key}
            )
            continue
        info["actual_support"].append(key)
        entry: dict[str, Any] = {
            "key_ast": _ast_dump(key_node),
            "tuple_ast": _ast_dump(value_node),
            "expression_node": None,
            "f0_power": None,
            "lambda_power": None,
            "f0_power_ast": None,
            "lambda_power_ast": None,
            "entry_error": None,
        }
        if not isinstance(value_node, ast.Tuple) or len(value_node.elts) != 3:
            entry["entry_error"] = "entry must be a three-item tuple"
        else:
            entry["f0_power_ast"] = _ast_dump(value_node.elts[0])
            entry["lambda_power_ast"] = _ast_dump(value_node.elts[1])
            entry["f0_power"] = _integer(value_node.elts[0])
            entry["lambda_power"] = _integer(value_node.elts[1])
            entry["expression_node"] = value_node.elts[2]
            if entry["f0_power"] is None or entry["lambda_power"] is None:
                entry["entry_error"] = "F0/Lambda powers must be integer constants"
        info["entries"][key] = entry
    return info


def _power_receipt(
    high_entry: dict[str, Any] | None,
    low_entry: dict[str, Any] | None,
    expected: tuple[int, int],
    assumptions: list[str],
) -> dict[str, Any]:
    high = None if high_entry is None else (high_entry["f0_power"], high_entry["lambda_power"])
    low = None if low_entry is None else (low_entry["f0_power"], low_entry["lambda_power"])
    return {
        "expected": {"F0_power": expected[0], "Lambda_power": expected[1]},
        "high": {
            "F0_power": high[0] if high is not None else None,
            "Lambda_power": high[1] if high is not None else None,
            "F0_power_ast": high_entry["f0_power_ast"] if high_entry else None,
            "Lambda_power_ast": high_entry["lambda_power_ast"] if high_entry else None,
        },
        "low": {
            "F0_power": low[0] if low is not None else None,
            "Lambda_power": low[1] if low is not None else None,
            "F0_power_ast": low_entry["f0_power_ast"] if low_entry else None,
            "Lambda_power_ast": low_entry["lambda_power_ast"] if low_entry else None,
        },
        "matched": high == expected and low == expected,
        "assumptions": assumptions,
    }


def _term_receipt(
    section: str,
    key: str,
    high_entry: dict[str, Any] | None,
    low_entry: dict[str, Any] | None,
    assumptions: list[str],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    expected = EXPECTED_POWERS[section][key]
    high_formula = _formula_receipt(
        high_entry["expression_node"] if high_entry else None,
        "high",
        assumptions,
    )
    low_formula = _formula_receipt(
        low_entry["expression_node"] if low_entry else None,
        "low",
        assumptions,
    )
    powers = _power_receipt(high_entry, low_entry, expected, assumptions)
    formula_match = (
        high_formula["present"]
        and low_formula["present"]
        and high_formula["error"] is None
        and low_formula["error"] is None
        and high_formula["normalized_value"] == low_formula["normalized_value"]
    )
    matched = powers["matched"] and formula_match
    differences: list[dict[str, Any]] = []
    if high_entry is None or low_entry is None:
        differences.append(
            {
                "kind": "missing_matrix_item",
                "path": section + "." + key,
                "high_present": high_entry is not None,
                "low_present": low_entry is not None,
            }
        )
    if not powers["matched"]:
        differences.append(
            {
                "kind": "power_mismatch",
                "path": section + "." + key,
                "expected": powers["expected"],
                "high": powers["high"],
                "low": powers["low"],
            }
        )
    if not formula_match:
        differences.append(
            {
                "kind": "formula_mismatch",
                "path": section + "." + key,
                "high": high_formula["normalized"],
                "low": low_formula["normalized"],
                "high_error": high_formula["error"],
                "low_error": low_formula["error"],
            }
        )
    return (
        {
            "section": section,
            "key": key,
            "matched": matched,
            "powers": powers,
            "formula": {"high": high_formula, "low": low_formula, "matched": formula_match},
            "high_entry_ast": high_entry["tuple_ast"] if high_entry else None,
            "low_entry_ast": low_entry["tuple_ast"] if low_entry else None,
            "assumptions": assumptions,
            "differences": differences,
        },
        differences,
    )


def _support_receipt(
    section: str,
    high_info: dict[str, Any],
    low_info: dict[str, Any],
    assumptions: list[str],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    expected = list(EXPECTED_SUPPORT[section])
    high_actual = list(high_info.get("actual_support", []))
    low_actual = list(low_info.get("actual_support", []))
    high_missing = [key for key in expected if key not in high_actual]
    low_missing = [key for key in expected if key not in low_actual]
    high_unexpected = [key for key in high_actual if key not in expected]
    low_unexpected = [key for key in low_actual if key not in expected]
    matched = not high_missing and not low_missing and not high_unexpected and not low_unexpected
    differences: list[dict[str, Any]] = []
    if not matched:
        differences.append(
            {
                "kind": "support_mismatch",
                "path": section,
                "expected": expected,
                "high": high_actual,
                "low": low_actual,
                "high_missing": high_missing,
                "low_missing": low_missing,
                "high_unexpected": high_unexpected,
                "low_unexpected": low_unexpected,
            }
        )
    return (
        {
            "expected": expected,
            "high": high_actual,
            "low": low_actual,
            "high_missing": high_missing,
            "low_missing": low_missing,
            "high_unexpected": high_unexpected,
            "low_unexpected": low_unexpected,
            "matched": matched,
            "assumptions": assumptions,
        },
        differences,
    )


def _intermediate_receipt(
    name: str,
    high_node: ast.Assign | None,
    low_node: ast.Assign | None,
    assumptions: list[str],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    high_formula = _formula_receipt(high_node.value if high_node is not None else None, "high", assumptions)
    low_formula = _formula_receipt(low_node.value if low_node is not None else None, "low", assumptions)
    matched = (
        high_formula["present"]
        and low_formula["present"]
        and high_formula["error"] is None
        and low_formula["error"] is None
        and high_formula["normalized_value"] == low_formula["normalized_value"]
    )
    differences: list[dict[str, Any]] = []
    if high_node is None or low_node is None:
        differences.append(
            {
                "kind": "missing_intermediate",
                "path": name,
                "high_present": high_node is not None,
                "low_present": low_node is not None,
            }
        )
    if not matched:
        differences.append(
            {
                "kind": "intermediate_mismatch",
                "path": name,
                "high": high_formula["normalized"],
                "low": low_formula["normalized"],
                "high_error": high_formula["error"],
                "low_error": low_formula["error"],
            }
        )
    return (
        {
            "name": name,
            "matched": matched,
            "high_assignment": _assignment_receipt(high_node),
            "low_assignment": _assignment_receipt(low_node),
            "formula": {"high": high_formula, "low": low_formula, "matched": matched},
            "assumptions": assumptions,
            "differences": differences,
        },
        differences,
    )


def _source_receipt(label: str, tree: ast.Module | None, issue: dict[str, Any] | None) -> dict[str, Any]:
    return {
        "label": label,
        "parsed": tree is not None,
        "top_level_statement_count": len(tree.body) if tree is not None else None,
        "error": issue,
    }


def build(source_high: str, source_low: str) -> dict[str, Any]:
    """Return a conditional source-string AST identity receipt.

    ``source_high`` is the complete text of
    ``n1_core_low_source_matrix_v2.py`` and ``source_low`` is the complete text
    of ``n1_core_low_matrix_v3.py``.  No producer module is imported and no
    expression is evaluated.  A changed support, power, anchor, or formula is
    represented in ``differences`` and yields ``partial`` or ``rejected``.
    """

    high_tree, high_issue = _parse_source(source_high, HIGH_LABEL)
    low_tree, low_issue = _parse_source(source_low, LOW_LABEL)
    differences: list[dict[str, Any]] = []
    if high_issue is not None:
        differences.append(high_issue)
    if low_issue is not None:
        differences.append(low_issue)

    high_dicts = {
        section: _dict_info(high_tree, "c0_jets" if section == "C0" else "b1_jets", HIGH_LABEL)
        for section in ("C0", "B1")
    }
    low_dicts = {
        section: _dict_info(low_tree, "c0_jets" if section == "C0" else "b1_jets", LOW_LABEL)
        for section in ("C0", "B1")
    }
    for section in ("C0", "B1"):
        differences.extend(high_dicts[section]["errors"])
        differences.extend(low_dicts[section]["errors"])

    assumptions = list(ASSUMPTIONS)
    high_scope = high_dicts["C0"].get("scope_node")
    low_scope = low_dicts["C0"].get("scope_node")
    anchors, anchor_differences = _anchor_receipt(
        high_scope,
        low_scope,
        _sqrt_anchor(high_dicts, low_dicts),
    )
    differences.extend(anchor_differences)

    supports: dict[str, Any] = {}
    matrix_items: dict[str, dict[str, Any]] = {}
    for section in ("C0", "B1"):
        support, support_differences = _support_receipt(
            section, high_dicts[section], low_dicts[section], assumptions
        )
        supports[section] = support
        differences.extend(support_differences)
        matrix_items[section] = {}
        for key in EXPECTED_SUPPORT[section]:
            high_entry = high_dicts[section]["entries"].get(key)
            low_entry = low_dicts[section]["entries"].get(key)
            item, item_differences = _term_receipt(
                section, key, high_entry, low_entry, assumptions
            )
            matrix_items[section][key] = item
            differences.extend(item_differences)

    high_intermediate_nodes: dict[str, ast.Assign | None] = {}
    low_intermediate_nodes: dict[str, ast.Assign | None] = {}
    for name in INTERMEDIATE_NAMES:
        high_intermediate_nodes[name], issue = _single_assignment(high_scope, name, HIGH_LABEL)
        if issue is not None:
            differences.append(issue)
        low_intermediate_nodes[name], issue = _single_assignment(low_scope, name, LOW_LABEL)
        if issue is not None:
            differences.append(issue)
    intermediates: dict[str, Any] = {}
    for name in INTERMEDIATE_NAMES:
        record, intermediate_differences = _intermediate_receipt(
            name,
            high_intermediate_nodes[name],
            low_intermediate_nodes[name],
            assumptions,
        )
        intermediates[name] = record
        differences.extend(intermediate_differences)

    checks: list[bool] = [
        high_tree is not None,
        low_tree is not None,
        anchors["matched"],
    ]
    checks.extend(supports[section]["matched"] for section in ("C0", "B1"))
    checks.extend(matrix_items[section][key]["matched"] for section in ("C0", "B1") for key in EXPECTED_SUPPORT[section])
    checks.extend(intermediates[name]["matched"] for name in INTERMEDIATE_NAMES)
    passed_count = sum(1 for check in checks if check)
    passed = bool(checks) and all(checks) and not differences
    if passed:
        status = "proved"
    elif high_tree is not None or low_tree is not None:
        status = "partial"
    else:
        status = "rejected"

    return {
        "schema": SCHEMA,
        "status": status,
        "passed": passed,
        "ok": passed,
        "conditional": True,
        "claim": "source-string AST formula identity under explicit symbol correspondence",
        "source": {
            "high": _source_receipt(HIGH_LABEL, high_tree, high_issue),
            "low": _source_receipt(LOW_LABEL, low_tree, low_issue),
        },
        "symbol_correspondence": dict(SYMBOL_CORRESPONDENCE),
        "literal_correspondence": dict(LITERAL_CORRESPONDENCE),
        "expected_support": {section: list(value) for section, value in EXPECTED_SUPPORT.items()},
        "expected_powers": {
            section: {
                key: {"F0_power": power[0], "Lambda_power": power[1]}
                for key, power in values.items()
            }
            for section, values in EXPECTED_POWERS.items()
        },
        "expected_matrix_item_count": EXPECTED_MATRIX_ITEM_COUNT,
        "anchors": anchors,
        "supports": supports,
        "matrix_items": matrix_items,
        "terms": matrix_items,
        "intermediates": intermediates,
        "checks": {"passed": passed_count, "total": len(checks)},
        "differences": differences,
        "assumptions": assumptions,
        "unproved": list(UNPROVED),
    }


__all__ = [
    "ASSUMPTIONS",
    "EXPECTED_B1_SUPPORT",
    "EXPECTED_C0_SUPPORT",
    "EXPECTED_MATRIX_ITEM_COUNT",
    "EXPECTED_POWERS",
    "EXPECTED_SUPPORT",
    "INTERMEDIATE_NAMES",
    "LITERAL_CORRESPONDENCE",
    "SCHEMA",
    "SYMBOL_CORRESPONDENCE",
    "UNPROVED",
    "SourceIdentityError",
    "build",
]
