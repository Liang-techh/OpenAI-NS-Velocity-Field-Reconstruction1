"""Static O2-to-n1 parameter evidence bridge.

This leaf consumes evidence mappings and source *text* supplied by its caller.
It never opens a path, imports a producer, evaluates a graph, decodes a
receipt, or compares interval endpoints.  Its positive claim is deliberately
small: the supplied O2 branch metadata and the supplied n1 core source use the
same certified expression recipes.  A source-level data-flow check for the
whole-Z bridge is reported separately and may remain ``unresolved``.  None of
these facts proves that two live parameter objects, endpoints, physical bases,
or native ``RQy`` values are the same object.
"""

from __future__ import annotations

import ast
import copy
import hashlib
from collections.abc import Mapping
from typing import Any

SCHEMA = "original_o2_core_bridge_v2"
EXPECTED_BRANCH = "exp(-4logPstar-30)"
CORE_SOURCE_PATH = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
BRIDGE_SOURCE_PATH = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source.py"
O2_SOURCE_PIN_NAME = "lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame.py"
O2_SOURCE_PIN_SHA256 = "dc2378e7cf690c879b48cb0ba0bbccb1fcf00d26b563269ca447fcc03f5c8a2c"
CORE_SOURCE_PIN_NAME = "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py"
CORE_SOURCE_PIN_SHA256 = "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59"

_MISSING = object()


def _copy(value: Any) -> Any:
    try:
        return copy.deepcopy(value)
    except Exception:
        return value


def _sha_text(value: Any) -> str | None:
    return hashlib.sha256(value.encode("utf-8")).hexdigest() if isinstance(value, str) else None


def _path(node: ast.AST) -> tuple[str, ...] | None:
    if isinstance(node, ast.Name):
        return (node.id,)
    if isinstance(node, ast.Attribute):
        prefix = _path(node.value)
        return (*prefix, node.attr) if prefix is not None else None
    return None


def _constant(node: ast.AST, value: Any) -> bool:
    return isinstance(node, ast.Constant) and type(node.value) is type(value) and node.value == value


def _subscript_key(node: ast.AST) -> str | None:
    if not isinstance(node, ast.Subscript):
        return None
    value: ast.AST = node.slice
    if isinstance(value, ast.Index):
        value = value.value
    return value.value if isinstance(value, ast.Constant) and type(value.value) is str else None


def _subscript_path(node: ast.AST) -> tuple[str, ...] | None:
    keys: list[str] = []
    current: ast.AST = node
    while isinstance(current, ast.Subscript):
        key = _subscript_key(current)
        if key is None:
            return None
        keys.append(key)
        current = current.value
    prefix = _path(current)
    return (*prefix, *reversed(keys)) if prefix is not None else None


def _field_path(node: ast.AST) -> tuple[str, ...] | None:
    return _path(node) or _subscript_path(node)


def _call(node: ast.AST, path: tuple[str, ...], argc: int | None = None) -> ast.Call | None:
    if not isinstance(node, ast.Call) or _path(node.func) != path or node.keywords:
        return None
    if argc is not None and len(node.args) != argc:
        return None
    return node


def _assignment_targets(node: ast.AST) -> list[ast.AST]:
    if isinstance(node, ast.Assign):
        return list(node.targets)
    if isinstance(node, ast.AnnAssign):
        return [node.target]
    return []


def _target_name(node: ast.AST) -> str | None:
    return node.id if isinstance(node, ast.Name) else None


def _functions(tree: ast.AST, name: str) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    return [node for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name]


def _assignment_values(function: ast.AST, target: str) -> list[ast.AST]:
    values: list[ast.AST] = []
    for node in ast.walk(function):
        if any(_target_name(item) == target for item in _assignment_targets(node)):
            values.append(node.value)  # type: ignore[attr-defined]
    return values


def _assignment_map(function: ast.AST) -> dict[str, list[ast.AST]]:
    result: dict[str, list[ast.AST]] = {}
    for node in ast.walk(function):
        for item in _assignment_targets(node):
            name = _target_name(item)
            if name is not None:
                result.setdefault(name, []).append(node.value)  # type: ignore[attr-defined]
    return result


def _c_number(node: ast.AST, number: int) -> bool:
    call = _call(node, ("c",), 1)
    return call is not None and _constant(call.args[0], number)


def _expected_logp(node: ast.AST) -> bool:
    outer = _call(node, ("g", "add"), 2)
    if outer is None:
        return False
    unary = _call(outer.args[0], ("g", "unary"), 2)
    return unary is not None and _constant(unary.args[0], "exp") and _c_number(unary.args[1], 40) and _c_number(outer.args[1], 11)


def _expected_delta(node: ast.AST) -> bool:
    outer = _call(node, ("g", "unary"), 2)
    if outer is None or not _constant(outer.args[0], "exp"):
        return False
    neg = _call(outer.args[1], ("g", "neg"), 1)
    add = _call(neg.args[0], ("g", "add"), 2) if neg is not None else None
    mul = _call(add.args[0], ("g", "mul"), 2) if add is not None else None
    return add is not None and mul is not None and _c_number(mul.args[0], 4) and _subscript_path(mul.args[1]) == ("p", "logP") and _c_number(add.args[1], 30)


def _comparison(node: ast.AST, key: str, expected: str, unequal: bool = False) -> bool:
    if not isinstance(node, ast.Compare) or len(node.ops) != 1 or len(node.comparators) != 1:
        return False
    operator = node.ops[0]
    if unequal and not isinstance(operator, ast.NotEq):
        return False
    if not unequal and not isinstance(operator, ast.Eq):
        return False
    return _subscript_path(node.left) == ("p", key) and isinstance(node.comparators[0], ast.Name) and node.comparators[0].id == expected


def _bool_parts(node: ast.AST, operator: type[ast.boolop]) -> list[ast.AST]:
    if isinstance(node, ast.BoolOp) and isinstance(node.op, operator):
        result: list[ast.AST] = []
        for child in node.values:
            result.extend(_bool_parts(child, operator))
        return result
    return [node]


def _positive_guard(node: ast.AST) -> bool:
    parts = _bool_parts(node, ast.And)
    return any(_comparison(item, "logP", "expected_logP") for item in parts) and any(_comparison(item, "delta", "expected_delta") for item in parts)


def _negative_guard(node: ast.AST) -> bool:
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.Not):
        return _positive_guard(node.operand)
    parts = _bool_parts(node, ast.Or)
    return any(_comparison(item, "logP", "expected_logP", True) for item in parts) and any(_comparison(item, "delta", "expected_delta", True) for item in parts)


def _parse_source(source_text: Any, label: str) -> tuple[ast.Module | None, dict[str, Any], list[str]]:
    info = {"path": label, "source_text_supplied": isinstance(source_text, str), "source_executed": False, "sha256": _sha_text(source_text), "parsed": False}
    if not isinstance(source_text, str):
        return None, info, ["source text is not a string"]
    try:
        tree = ast.parse(source_text, filename=label)
    except SyntaxError as exc:
        return None, info, ["source text is not parseable: " + str(exc)]
    info["parsed"] = True
    return tree, info, []


def _core_ast_contract(source_text: Any) -> dict[str, Any]:
    tree, info, errors = _parse_source(source_text, CORE_SOURCE_PATH)
    result: dict[str, Any] = {**info, "function_count": 0, "parameter_identity_found": False, "expected_logP_is_exp_40_plus_11": False, "expected_delta_is_exp_negative_4_logP_minus_30": False, "parameter_binding_guard": False, "guard_uses_expected_logP_and_expected_delta": False, "parameter_identity_ast_contract": False, "errors": list(errors), "source_imported": False, "graph_traversed": False}
    if tree is None:
        result["status"] = "mismatch" if errors else "incomplete"
        return result
    methods = _functions(tree, "parameter_identity")
    result["function_count"] = len(methods)
    if len(methods) != 1:
        result["errors"].append("parameter_identity definition is missing or duplicated")
        result["status"] = "incomplete"
        return result
    method = methods[0]
    result["parameter_identity_found"] = True
    result["ast_sha256"] = hashlib.sha256(ast.dump(method, include_attributes=False).encode("utf-8")).hexdigest()
    logp_values = _assignment_values(method, "expected_logP")
    delta_values = _assignment_values(method, "expected_delta")
    result["expected_logP_assignment_count"] = len(logp_values)
    result["expected_delta_assignment_count"] = len(delta_values)
    if len(logp_values) == 1:
        result["expected_logP_is_exp_40_plus_11"] = _expected_logp(logp_values[0])
    if len(delta_values) == 1:
        result["expected_delta_is_exp_negative_4_logP_minus_30"] = _expected_delta(delta_values[0])
    positive = False
    negative = False
    for node in ast.walk(method):
        if isinstance(node, ast.Assert):
            positive = positive or _positive_guard(node.test)
        elif isinstance(node, ast.If) and any(isinstance(item, ast.Raise) for item in ast.walk(node)):
            negative = negative or _negative_guard(node.test)
    result["parameter_binding_guard"] = positive or negative
    result["guard_uses_expected_logP_and_expected_delta"] = result["parameter_binding_guard"]
    result["parameter_identity_ast_contract"] = bool(result["expected_logP_is_exp_40_plus_11"] and result["expected_delta_is_exp_negative_4_logP_minus_30"] and result["parameter_binding_guard"])
    if not result["expected_logP_is_exp_40_plus_11"]:
        result["errors"].append("expected_logP expression does not match the certified recipe")
    if not result["expected_delta_is_exp_negative_4_logP_minus_30"]:
        result["errors"].append("expected_delta expression does not match the certified recipe")
    if not result["parameter_binding_guard"]:
        result["errors"].append("parameter guard does not bind p['logP'] and p['delta']")
    result["status"] = "proved" if result["parameter_identity_ast_contract"] else "mismatch"
    return result


def _expr_tokens(
    node: ast.AST,
    assignments: Mapping[str, list[ast.AST]],
    seen: set[str] | None = None,
) -> set[str]:
    """Collect source-name tokens through local aliases without evaluating code."""

    seen = set() if seen is None else set(seen)
    tokens: set[str] = set()
    path = _field_path(node)
    if path is not None:
        tokens.update(path)
    # Follow the base name of an attribute/subscript through local aliases;
    # a renamed provider must not be mistaken for a distinct source recipe.
    if isinstance(node, (ast.Attribute, ast.Subscript)):
        tokens.update(_expr_tokens(node.value, assignments, seen))
    if isinstance(node, ast.Name) and node.id not in seen:
        seen.add(node.id)
        for value in assignments.get(node.id, ()):
            tokens.update(_expr_tokens(value, assignments, seen))
    return tokens


def _contains_field_read(
    node: ast.AST,
    field: str,
    required_tokens: set[str],
    assignments: Mapping[str, list[ast.AST]],
) -> bool:
    for item in ast.walk(node):
        path = _field_path(item)
        if path is not None and path[-1] == field:
            if required_tokens.issubset(_expr_tokens(item, assignments)):
                return True
        if not isinstance(item, ast.Call) or item.keywords or len(item.args) < 2:
            continue
        if not any(_constant(arg, field) for arg in item.args[1:]):
            continue
        if required_tokens.issubset(_expr_tokens(item.args[0], assignments)):
            return True
    return False


def _name_reaches(
    node: ast.AST,
    target_names: set[str],
    assignments: Mapping[str, list[ast.AST]],
    seen: set[str] | None = None,
) -> bool:
    seen = set() if seen is None else set(seen)
    if isinstance(node, ast.Name):
        if node.id in target_names:
            return True
        if node.id in seen:
            return False
        seen.add(node.id)
        return any(_name_reaches(value, target_names, assignments, seen) for value in assignments.get(node.id, ()))
    return any(_name_reaches(child, target_names, assignments, seen) for child in ast.iter_child_nodes(node))


def _read_assigned_names(
    field: str,
    required_tokens: set[str],
    assignments: Mapping[str, list[ast.AST]],
) -> set[str]:
    names: set[str] = set()
    for name, values in assignments.items():
        if any(_contains_field_read(value, field, required_tokens, assignments) for value in values):
            names.add(name)
    return names


def _core_constructor_assignment(function: ast.AST) -> list[tuple[str, ast.Call]]:
    found: list[tuple[str, ast.Call]] = []
    for node in ast.walk(function):
        targets = _assignment_targets(node)
        value = getattr(node, "value", None)
        if not targets or not isinstance(value, ast.Call):
            continue
        target_names = [_target_name(target) for target in targets]
        target_paths = [_field_path(target) for target in targets]
        target_is_core = any(name == "core" for name in target_names) or any(
            path is not None and path[-1] == "core" for path in target_paths
        )
        if not target_is_core:
            continue
        constructor_path = _path(value.func)
        keywords = {keyword.arg for keyword in value.keywords if keyword.arg is not None}
        constructor_is_core = constructor_path is not None and (
            constructor_path[-1] == "SimpleNamespace" or "core" in "".join(constructor_path).lower()
        )
        if not (constructor_is_core or {"delta", "logP", "logPstar", "parameter"}.intersection(keywords)):
            continue
        for name in target_names:
            if name == "core":
                found.append((name, value))
        if not any(name == "core" for name in target_names):
            found.append(("self.core", value))
    return found


def _core_binding_flow(
    call: ast.Call,
    delta_names: set[str],
    logp_names: set[str],
    assignments: Mapping[str, list[ast.AST]],
) -> tuple[bool, bool, bool, bool]:
    direct_delta = False
    direct_logp = False
    nested_parameter = False
    for keyword in call.keywords:
        if keyword.arg == "delta":
            direct_delta = _name_reaches(keyword.value, delta_names, assignments) or _contains_field_read(keyword.value, "log_delta", {"parameter_bounds"}, assignments)
        elif keyword.arg in {"logP", "logPstar", "logp"}:
            direct_logp = _name_reaches(keyword.value, logp_names, assignments) or _contains_field_read(keyword.value, "logPstar", {"core_transfer"}, assignments)
        elif keyword.arg == "parameter":
            nested_parameter = (
                (
                    _name_reaches(keyword.value, delta_names, assignments)
                    or _contains_field_read(keyword.value, "log_delta", {"parameter_bounds"}, assignments)
                )
                and (
                    _name_reaches(keyword.value, logp_names, assignments)
                    or _contains_field_read(keyword.value, "logPstar", {"core_transfer"}, assignments)
                )
            )
    # A parameter namespace can be nested in the outer core constructor.
    for item in ast.walk(call):
        if not isinstance(item, ast.Call) or item is call:
            continue
        fields = {keyword.arg for keyword in item.keywords if keyword.arg is not None}
        if "delta" not in fields or not ({"logP", "logPstar", "logp"} & fields):
            continue
        delta_value = next(keyword.value for keyword in item.keywords if keyword.arg == "delta")
        logp_value = next(keyword.value for keyword in item.keywords if keyword.arg in {"logP", "logPstar", "logp"})
        nested_parameter = nested_parameter or (
            (
                _name_reaches(delta_value, delta_names, assignments)
                or _contains_field_read(delta_value, "log_delta", {"parameter_bounds"}, assignments)
            )
            and (
                _name_reaches(logp_value, logp_names, assignments)
                or _contains_field_read(logp_value, "logPstar", {"core_transfer"}, assignments)
            )
        )
    return direct_delta, direct_logp, nested_parameter, (direct_delta or nested_parameter) and (direct_logp or nested_parameter)


def _self_core_connection(function: ast.AST, core_target: str) -> bool:
    for node in ast.walk(function):
        for target in _assignment_targets(node):
            path = _field_path(target)
            if path is None or path[-1] != "core":
                continue
            value = getattr(node, "value", None)
            if isinstance(value, ast.Name) and value.id == core_target:
                return True
    return False


def _bridge_ast_contract(source_text: Any) -> dict[str, Any]:
    tree, info, errors = _parse_source(source_text, BRIDGE_SOURCE_PATH)
    result: dict[str, Any] = {**info, "candidate_function_count": 0, "candidate_function_found": False, "source_logdelta_read": False, "source_log_delta_read": False, "source_logP_read": False, "source_logPstar_read": False, "core_instance_constructed": False, "core_instance_connected": False, "delta_dataflow_to_core": False, "logP_dataflow_to_core": False, "core_parameter_binding_guard": False, "actual_data_flow_proved": False, "status": "incomplete", "errors": list(errors), "source_imported": False, "graph_traversed": False}
    if tree is None:
        result["status"] = "mismatch" if errors else "incomplete"
        return result
    candidates: list[ast.AST] = []
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        assignments = _assignment_map(node)
        delta_read = _contains_field_read(node, "log_delta", {"parameter_bounds"}, assignments)
        logp_read = _contains_field_read(node, "logPstar", {"core_transfer"}, assignments)
        constructors = _core_constructor_assignment(node)
        if delta_read or logp_read or constructors:
            candidates.append(node)
    result["candidate_function_count"] = len(candidates)
    if len(candidates) != 1:
        result["errors"].append("bridge constructor/data-flow function is missing or duplicated")
        result["status"] = "incomplete"
        return result
    function = candidates[0]
    result["candidate_function_found"] = True
    result["function_name"] = function.name
    result["ast_sha256"] = hashlib.sha256(ast.dump(function, include_attributes=False).encode("utf-8")).hexdigest()
    assignments = _assignment_map(function)
    delta_names = _read_assigned_names("log_delta", {"parameter_bounds"}, assignments)
    logp_names = _read_assigned_names("logPstar", {"core_transfer"}, assignments)
    result["source_logdelta_read"] = bool(delta_names) or _contains_field_read(function, "log_delta", {"parameter_bounds"}, assignments)
    result["source_log_delta_read"] = result["source_logdelta_read"]
    result["source_logP_read"] = bool(logp_names) or _contains_field_read(function, "logPstar", {"core_transfer"}, assignments)
    result["source_logPstar_read"] = result["source_logP_read"]
    constructors = _core_constructor_assignment(function)
    result["core_instance_constructed"] = len(constructors) == 1
    if len(constructors) == 1:
        core_target, call = constructors[0]
        direct_delta, direct_logp, nested_parameter, parameter_guard = _core_binding_flow(call, delta_names, logp_names, assignments)
        result["delta_dataflow_to_core"] = direct_delta or nested_parameter
        result["logP_dataflow_to_core"] = direct_logp or nested_parameter
        result["core_parameter_binding_guard"] = parameter_guard
        result["core_instance_connected"] = core_target == "self.core" or _self_core_connection(function, core_target)
    result["actual_data_flow_proved"] = bool(result["source_logdelta_read"] and result["source_logP_read"] and result["core_instance_constructed"] and result["core_instance_connected"] and result["delta_dataflow_to_core"] and result["logP_dataflow_to_core"] and result["core_parameter_binding_guard"])
    if result["actual_data_flow_proved"]:
        result["status"] = "proved"
    else:
        result["status"] = "unresolved"
        result["errors"].append("bridge source does not provide an exact log_delta/logPstar-to-core data flow")
    return result


def _section_candidates(value: Any, required: tuple[str, ...]) -> list[tuple[str, Mapping[str, Any]]]:
    candidates: list[tuple[str, Mapping[str, Any]]] = []
    seen: set[int] = set()

    def visit(path: str, current: Any, depth: int) -> None:
        if depth > 3 or not isinstance(current, Mapping) or id(current) in seen:
            return
        seen.add(id(current))
        if all(key in current for key in required):
            candidates.append((path, current))
        for key in ("result", "check", "leaf", "o2_check", "core_check"):
            child = current.get(key, _MISSING)
            if isinstance(child, Mapping):
                visit(path + "." + key, child, depth + 1)

    visit("root", value, 0)
    return candidates


def _input_hash_pin(container: Mapping[str, Any], name: str, expected: str) -> dict[str, Any]:
    values: list[tuple[str, Any]] = []
    groups: list[tuple[str, Any]] = [("root", container)]
    for key in ("input_hashes", "source_hashes", "hashes"):
        group = container.get(key)
        if isinstance(group, Mapping):
            groups.append((key, group))
    for group_path, group in groups:
        if not isinstance(group, Mapping):
            continue
        for raw_key, raw_value in group.items():
            if not isinstance(raw_key, str):
                continue
            normalized = raw_key.replace("\\", "/")
            if normalized == name or normalized.endswith("/" + name) or normalized.rsplit("/", 1)[-1] == name:
                values.append((group_path + "." + raw_key, raw_value))
    unique = list(dict.fromkeys(value for _, value in values))
    return {
        "name": name,
        "expected_sha256": expected,
        "present": bool(values),
        "paths": [path for path, _ in values],
        "sha256": unique[0] if len(unique) == 1 else None,
        "matches": len(unique) == 1 and unique[0] == expected,
        "conflicting": len(unique) > 1,
        "source_pin_is_provenance_only": True,
    }


def _o2_contract(value: Any) -> dict[str, Any]:
    result: dict[str, Any] = {"present": False, "path": None, "source_definition_binding_passed": False, "selected_branch": None, "selected_branch_expected": EXPECTED_BRANCH, "selected_branch_matches": False, "no_rounded_branch_comparison": None, "no_rounded_branch_comparison_passed": False, "delta_branch_proof_passed": False, "o2_source_pin": None, "passed": False, "status": "incomplete", "errors": []}
    candidates = _section_candidates(value, ("source_definition_binding", "delta_branch_proof"))
    if len(candidates) != 1:
        result["errors"].append("O2 source-definition/branch evidence is missing or ambiguous")
        return result
    path, section = candidates[0]
    result["present"] = True
    result["path"] = path
    result["o2_source_pin"] = _input_hash_pin(section, O2_SOURCE_PIN_NAME, O2_SOURCE_PIN_SHA256)
    binding = section.get("source_definition_binding")
    branch = section.get("delta_branch_proof")
    if not isinstance(binding, Mapping) or not isinstance(branch, Mapping):
        result["errors"].append("O2 evidence sections are not mappings")
        return result
    result["source_definition_binding_passed"] = binding.get("passed") is True
    result["selected_branch"] = branch.get("selected_branch")
    result["selected_branch_matches"] = branch.get("selected_branch") == EXPECTED_BRANCH
    result["no_rounded_branch_comparison"] = branch.get("no_rounded_branch_comparison")
    result["no_rounded_branch_comparison_passed"] = branch.get("no_rounded_branch_comparison") is True
    result["delta_branch_proof_passed"] = branch.get("passed") is True
    if not result["source_definition_binding_passed"]:
        result["errors"].append("O2 source_definition_binding.passed is not true")
    if not result["selected_branch_matches"]:
        result["errors"].append("O2 selected branch differs from the certified exp(-4logPstar-30) branch")
    if not result["no_rounded_branch_comparison_passed"]:
        result["errors"].append("O2 branch comparison is missing or rounded")
    if not result["delta_branch_proof_passed"]:
        result["errors"].append("O2 delta_branch_proof.passed is not true")
    if not result["o2_source_pin"]["present"]:
        result["errors"].append("O2 source parameter-frame source pin is missing")
    elif not result["o2_source_pin"]["matches"]:
        result["errors"].append("O2 source parameter-frame source pin is mismatched")
    result["passed"] = not result["errors"]
    result["status"] = "proved" if result["passed"] else "mismatch"
    return result


def _core_gate_contract(value: Any) -> dict[str, Any]:
    result: dict[str, Any] = {"present": False, "path": None, "all_passed": None, "actual_oracle_passed": None, "identity_passed": None, "exact_graph_recipes_verified": None, "same_exact_recipes": None, "interval_overlap_not_used_as_function_identity": None, "o2_source_pin": None, "core_source_pin": None, "source_pins_passed": False, "passed": False, "status": "incomplete", "errors": []}
    roots: list[tuple[str, Mapping[str, Any]]] = []
    if isinstance(value, Mapping):
        roots.append(("root", value))
        if isinstance(value.get("result"), Mapping):
            roots.append(("root.result", value["result"]))
    candidates: list[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]] = []
    for path, root in roots:
        checks = root.get("checks")
        if not isinstance(checks, Mapping):
            continue
        oracle = checks.get("actual_oracle")
        if not isinstance(oracle, Mapping):
            continue
        identity = oracle.get("exact_parameter_and_native_P0_identity")
        if isinstance(identity, Mapping):
            candidates.append((path + ".checks.actual_oracle.exact_parameter_and_native_P0_identity", root, oracle, identity))
    if len(candidates) != 1:
        result["errors"].append("core exact_parameter_and_native_P0_identity gate is missing or ambiguous")
        return result
    path, root, oracle, identity = candidates[0]
    result["present"] = True
    result["path"] = path
    result["o2_source_pin"] = _input_hash_pin(root, O2_SOURCE_PIN_NAME, O2_SOURCE_PIN_SHA256)
    result["core_source_pin"] = _input_hash_pin(root, CORE_SOURCE_PIN_NAME, CORE_SOURCE_PIN_SHA256)
    result["source_pins_passed"] = bool(result["o2_source_pin"]["matches"] and result["core_source_pin"]["matches"])
    result["all_passed"] = root.get("all_passed")
    result["actual_oracle_passed"] = oracle.get("passed")
    result["identity_passed"] = identity.get("passed")
    result["exact_graph_recipes_verified"] = identity.get("exact_graph_logP_delta_logLambda_epsilon_recipes_verified")
    recipes = identity.get("actual_core_encloses_same_exact_recipes")
    if isinstance(recipes, Mapping):
        result["same_exact_recipes"] = all(recipes.get(key) is True for key in ("logP", "logLambda", "delta", "Lambda", "epsilon"))
    result["interval_overlap_not_used_as_function_identity"] = identity.get("interval_overlap_not_used_as_a_function_identity")
    required = (result["all_passed"] is True, result["actual_oracle_passed"] is True, result["identity_passed"] is True, result["exact_graph_recipes_verified"] is True, result["same_exact_recipes"] is True, result["interval_overlap_not_used_as_function_identity"] is True, result["source_pins_passed"] is True)
    if result["all_passed"] is False or result["actual_oracle_passed"] is False or result["identity_passed"] is False:
        result["errors"].append("core actual-oracle identity gate is explicitly failed")
        result["status"] = "mismatch"
    elif result["interval_overlap_not_used_as_function_identity"] is False:
        result["errors"].append("core gate permits interval overlap as identity")
        result["status"] = "mismatch"
    elif result["o2_source_pin"]["present"] is False or result["core_source_pin"]["present"] is False:
        result["errors"].append("core O2/core source pins are incomplete")
        result["status"] = "incomplete"
    elif result["source_pins_passed"] is False:
        result["errors"].append("core O2/core source pin is mismatched")
        result["status"] = "mismatch"
    elif all(required):
        result["passed"] = True
        result["status"] = "proved"
    else:
        result["errors"].append("core identity gate lacks a required exact recipe field")
        result["status"] = "incomplete"
    return result


def _lineage_leaf(value: Any) -> tuple[str | None, Mapping[str, Any] | None, list[str]]:
    candidates: list[tuple[str, Mapping[str, Any]]] = []
    if isinstance(value, Mapping):
        if (isinstance(value.get("result"), Mapping)
                and isinstance(value["result"].get("result"), Mapping)
                and isinstance(value["result"]["result"].get("leaf"), Mapping)):
            candidates.append(("result.result.leaf", value["result"]["result"]["leaf"]))
        if isinstance(value.get("result"), Mapping) and isinstance(value["result"].get("leaf"), Mapping):
            candidates.append(("result.leaf", value["result"]["leaf"]))
        if isinstance(value.get("leaf"), Mapping):
            candidates.append(("leaf", value["leaf"]))
        if value.get("schema") == "original_core_parameter_lineage_v1":
            candidates.append(("root", value))
    if not candidates:
        return None, None, ["lineage receipt result.leaf is missing"]
    if len(candidates) > 1:
        first = candidates[0][1]
        if any(item[1].get("status") != first.get("status") or item[1].get("passed") != first.get("passed") for item in candidates[1:]):
            return None, None, ["lineage receipt has conflicting leaf statuses"]
    return candidates[0][0], candidates[0][1], []


def _lineage_contract(value: Any) -> dict[str, Any]:
    path, leaf, errors = _lineage_leaf(value)
    result: dict[str, Any] = {"present": leaf is not None, "path": path, "status": leaf.get("status") if leaf is not None else None, "passed_field": leaf.get("passed") if leaf is not None else None, "proved_status_preserved": False, "interval_overlap_not_used_as_identity": True, "live_identity_claimed_by_receipt": False, "status_classification": "incomplete", "errors": list(errors)}
    if leaf is None:
        return result
    result["proved_status_preserved"] = leaf.get("status") == "proved" and leaf.get("passed") is True
    comparisons = leaf.get("comparisons")
    if isinstance(comparisons, Mapping) and comparisons.get("overlap_used_as_identity") is True:
        result["interval_overlap_not_used_as_identity"] = False
        result["errors"].append("lineage receipt promotes interval overlap to identity")
    provenance = leaf.get("provenance_scope")
    if isinstance(provenance, Mapping):
        forbidden = ("cross_run_object_identity_claimed", "live_endpoint_basis_identity", "native_corrected_RQy_identity", "source_function_point_value_identity", "interval_overlap_promoted_to_identity")
        result["live_identity_claimed_by_receipt"] = any(provenance.get(key) is True for key in forbidden)
        if result["live_identity_claimed_by_receipt"]:
            result["errors"].append("lineage receipt contains a forbidden live identity claim")
    if leaf.get("status") == "mismatch":
        result["errors"].append("lineage receipt status is mismatch")
        result["status_classification"] = "mismatch"
    elif leaf.get("status") != "proved" or leaf.get("passed") is not True:
        result["errors"].append("lineage receipt proved status or passed field is missing")
    elif result["errors"]:
        result["status_classification"] = "mismatch"
    else:
        result["status_classification"] = "proved"
    return result


def _dedupe(values: list[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))


def verify(
    o2_check: Mapping[str, Any],
    core_check: Mapping[str, Any],
    lineage_receipt: Mapping[str, Any],
    core_source_text: str,
    bridge_source_text: str,
) -> dict[str, Any]:
    """Return static evidence for the O2-to-core parameter bridge."""

    o2 = _o2_contract(o2_check)
    core = _core_gate_contract(core_check)
    lineage = _lineage_contract(lineage_receipt)
    core_ast = _core_ast_contract(core_source_text)
    bridge_ast = _bridge_ast_contract(bridge_source_text)
    o2_pin_value = o2.get("o2_source_pin", {}).get("sha256") if isinstance(o2.get("o2_source_pin"), Mapping) else None
    core_pin_value = core.get("o2_source_pin", {}).get("sha256") if isinstance(core.get("o2_source_pin"), Mapping) else None
    source_pin_consistency = {
        "o2_report_pin": o2_pin_value,
        "core_report_pin": core_pin_value,
        "same_source_pin": o2_pin_value is not None and o2_pin_value == core_pin_value,
        "source_pin_is_provenance_only": True,
    }
    if o2_pin_value is not None and core_pin_value is not None and o2_pin_value != core_pin_value:
        core["errors"].append("O2 source pin differs between O2 and core reports")
        core["status"] = "mismatch"
        core["passed"] = False

    mismatches = (
        o2["status"] == "mismatch"
        or core["status"] == "mismatch"
        or lineage["status_classification"] == "mismatch"
        or core_ast["status"] == "mismatch"
        or bridge_ast["status"] == "mismatch"
    )
    incomplete = (
        o2["status"] == "incomplete"
        or core["status"] == "incomplete"
        or lineage["status_classification"] == "incomplete"
        or core_ast["status"] == "incomplete"
        or bridge_ast["status"] == "unresolved"
    )
    if mismatches:
        status = "mismatch"
    elif incomplete:
        status = "incomplete"
    else:
        status = "proved"

    recipe_evidence = bool(
        o2["passed"] is True
        and core["passed"] is True
        and core_ast["parameter_identity_ast_contract"] is True
    )
    errors = _dedupe([
        *o2["errors"],
        *core["errors"],
        *lineage["errors"],
        *core_ast["errors"],
        *bridge_ast["errors"],
    ])
    unresolved = {
        "live_parameter_object": True,
        "endpoint": True,
        "native_RQy": True,
        "physical_basis": True,
        "identity": True,
        "cancellation": True,
    }
    proof_scope = {
        "o2_strict_branch_certified": o2["passed"] is True,
        "o2_core_source_pins_consistent": source_pin_consistency["same_source_pin"],
        "core_exact_parameter_and_native_P0_identity_gate_certified": core["passed"] is True,
        "core_source_expression_recipes_exact": core_ast["parameter_identity_ast_contract"] is True,
        "expression_recipes_certified_and_same": recipe_evidence,
        "lineage_status_proved_preserved": lineage["proved_status_preserved"] is True,
        "bridge_source_dataflow_proved": bridge_ast["actual_data_flow_proved"] is True,
        "live_parameter_object_identity": False,
        "actual_live_endpoint_owner_identity": False,
        "live_endpoint_basis_identity": False,
        "native_corrected_RQy_identity": False,
        "physical_basis_identity": False,
        "cancellation_claimed": False,
        "interval_equality_promoted_to_identity": False,
        "parameter_spectral_lineage_promoted_to_identity": False,
        "old_errors_removed": False,
        "old_212_terms_and_all_errors_unchanged": True,
    }
    return {
        "schema": SCHEMA,
        "passed": status == "proved",
        "status": status,
        "errors": errors,
        "o2_branch_proof": o2,
        "o2_check": o2,
        "core_identity_gate": core,
        "exact_parameter_and_native_P0_identity": core,
        "lineage_receipt": lineage,
        "core_source_AST": core_ast,
        "bridge_source_AST": bridge_ast,
        "source_AST": {"core": core_ast, "bridge": bridge_ast},
        "source_pin_consistency": source_pin_consistency,
        "expression_recipe_evidence": {
            "selected_branch": EXPECTED_BRANCH,
            "o2_branch_matches_core_delta_recipe": recipe_evidence,
            "claim_scope": "certified expression recipes only; no live object identity",
        },
        "proof_scope": proof_scope,
        "comparisons": {
            "overlap_used_as_identity": False,
            "interval_equality_promoted_to_identity": False,
            "parameter_spectral_lineage_used_as_identity": False,
        },
        "provenance_scope": {
            "o2_branch_metadata_only": True,
            "core_gate_metadata_only": True,
            "lineage_receipt_status_preserved": lineage["proved_status_preserved"] is True,
            "cross_source_live_parameter_identity": False,
            "live_endpoint_basis_identity": False,
            "native_corrected_RQy_identity": False,
            "physical_basis_identity": False,
            "cancellation_identity": False,
            "interval_overlap_promoted_to_identity": False,
            "rounding_promoted_to_identity": False,
        },
        "unresolved_scope": unresolved,
        "old_212_terms_and_all_errors_unchanged": True,
        "old_errors_removed": False,
        "source_execution": False,
        "real_receipt_opened_here": False,
        "producer_imported": False,
        "source_imported": False,
        "receipt_decoded_here": False,
        "graph_traversed": False,
        "interval_comparison_performed": False,
        "live_identity_claimed": False,
    }


bind = verify


__all__ = ["SCHEMA", "EXPECTED_BRANCH", "verify", "bind"]


