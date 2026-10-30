"""Restricted AST proof for the old/current ``Tw`` recipes.

This module deliberately works on source text and decoded, already-authenticated
opaque references only.  It never imports a historical producer, evaluates a
graph, or interprets a graph-node id as a map/index.  The proof established here
is formula-level; it does not certify that two independently restored graph
objects have the same runtime meaning.
"""
from __future__ import annotations

import ast
from dataclasses import dataclass
from fractions import Fraction
from typing import Any, Mapping, Sequence


SCHEMA = "rp_tw_recipe_proof_v1"


class RecipeProofError(ValueError):
    """Raised when the closed, supported recipe language is not met."""


@dataclass(frozen=True)
class Expr:
    op: str
    args: tuple["Expr", ...] = ()
    value: Fraction | None = None
    real: bool = True


def _const(value: Fraction | int | str) -> Expr:
    return Expr("const", value=Fraction(value))


ZERO = _const(0)
ONE = _const(1)


def _is_zero(value: Expr) -> bool:
    return value.op == "const" and value.value == 0


def _is_one(value: Expr) -> bool:
    return value.op == "const" and value.value == 1


def _positive(value: Expr) -> bool:
    if value.op == "const":
        return value.value is not None and value.value > 0
    if value.op == "exp":
        return True
    if value.op == "mul":
        return all(_positive(a) for a in value.args)
    return False


def _neg(value: Expr) -> Expr:
    if _is_zero(value):
        return ZERO
    if value.op == "const":
        return _const(-value.value)
    if value.op == "neg":
        return value.args[0]
    return Expr("neg", (value,), real=value.real)


def _add(*values: Expr) -> Expr:
    flat: list[Expr] = []
    constant = Fraction(0)
    for value in values:
        if value.op == "add":
            flat.extend(value.args)
        elif value.op == "const":
            constant += value.value
        else:
            flat.append(value)
    if constant:
        flat.append(_const(constant))
    flat = [a for a in flat if not _is_zero(a)]
    if not flat:
        return ZERO
    if len(flat) == 1:
        return flat[0]
    return Expr("add", tuple(sorted(flat, key=expr_text)), real=all(a.real for a in flat))


def _mul(*values: Expr) -> Expr:
    flat: list[Expr] = []
    constant = Fraction(1)
    for value in values:
        if _is_zero(value):
            return ZERO
        if value.op == "mul":
            flat.extend(value.args)
        elif value.op == "const":
            constant *= value.value
        else:
            flat.append(value)
    if constant == 0:
        return ZERO
    if constant != 1 or not flat:
        flat.append(_const(constant))
    flat = [a for a in flat if not _is_one(a)]
    if not flat:
        return ONE
    if len(flat) == 1:
        return flat[0]
    return Expr("mul", tuple(sorted(flat, key=expr_text)), real=all(a.real for a in flat))


def _sub(left: Expr, right: Expr) -> Expr:
    return _add(left, _neg(right))


def _quotient(left: Expr, right: Expr) -> Expr:
    if _is_zero(left):
        return ZERO
    if _is_one(right):
        return left
    return Expr("quotient", (left, right), real=left.real and right.real)


def _exp(value: Expr) -> Expr:
    return Expr("exp", (value,), real=value.real)


def _log(value: Expr) -> Expr:
    if value.op == "exp" and value.real:
        return value.args[0]
    if not _positive(value):
        raise RecipeProofError("log argument has no restricted positive-domain proof")
    return Expr("log", (value,), real=True)


def expr_text(value: Expr) -> str:
    """Stable, non-evaluating rendering used in receipts and comparisons."""
    if value.op == "const":
        if value.value.denominator == 1:
            return str(value.value.numerator)
        return f"{value.value.numerator}/{value.value.denominator}"
    if value.op == "neg":
        return f"-({expr_text(value.args[0])})"
    if value.op == "add":
        return "(" + " + ".join(expr_text(a) for a in value.args) + ")"
    if value.op == "mul":
        return "(" + " * ".join(expr_text(a) for a in value.args) + ")"
    if value.op == "quotient":
        return f"({expr_text(value.args[0])} / {expr_text(value.args[1])})"
    return f"{value.op}({expr_text(value.args[0])})"


def _parse_constant(node: ast.AST) -> Expr:
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = _parse_constant(node.operand)
        return _const(-value.value) if isinstance(node.op, ast.USub) else value
    if not isinstance(node, ast.Constant) or isinstance(node.value, bool):
        raise RecipeProofError("graph constants must be exact integer/rational literals")
    if not isinstance(node.value, (int, str)):
        raise RecipeProofError("graph constants cannot be floating or computed values")
    try:
        return _const(node.value)
    except (TypeError, ValueError, ZeroDivisionError) as exc:
        raise RecipeProofError("invalid exact graph constant") from exc


def _name(node: ast.AST) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "self":
        return "self." + node.attr
    raise RecipeProofError("only local names and self.<field> names are supported")


def _translate(node: ast.AST, env: Mapping[str, Expr]) -> Expr:
    if isinstance(node, ast.Constant):
        # Constants in an expression are accepted only as exact numerics.  A
        # source string is not silently promoted to an unknown symbolic leaf.
        return _parse_constant(node)
    if isinstance(node, (ast.Name, ast.Attribute)):
        key = _name(node)
        if key not in env:
            raise RecipeProofError("unbound source expression: " + key)
        return env[key]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = _translate(node.operand, env)
        return _neg(value) if isinstance(node.op, ast.USub) else value
    if isinstance(node, ast.BinOp):
        left, right = _translate(node.left, env), _translate(node.right, env)
        if isinstance(node.op, ast.Add):
            return _add(left, right)
        if isinstance(node.op, ast.Sub):
            return _sub(left, right)
        if isinstance(node.op, ast.Mult):
            return _mul(left, right)
        if isinstance(node.op, ast.Div):
            return _quotient(left, right)
        raise RecipeProofError("unsupported arithmetic operator in source recipe")
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
        if not isinstance(node.func.value, ast.Name) or node.func.value.id != "g":
            raise RecipeProofError("unknown function node in source recipe")
        op = node.func.attr
        if any(keyword.value is not None for keyword in node.keywords):
            raise RecipeProofError("keyword graph operands are not supported")
        if op == "constant" and len(node.args) == 1:
            return _parse_constant(node.args[0])
        if op == "unary" and len(node.args) == 2:
            if not isinstance(node.args[0], ast.Constant) or not isinstance(node.args[0].value, str):
                raise RecipeProofError("unary operation name must be a literal")
            operand = _translate(node.args[1], env)
            if node.args[0].value == "exp":
                return _exp(operand)
            if node.args[0].value == "log":
                return _log(operand)
            raise RecipeProofError("unsupported unary graph operation")
        if op in {"add", "sub", "mul"} and len(node.args) == 2:
            left, right = (_translate(a, env) for a in node.args)
            return {"add": _add, "sub": _sub, "mul": _mul}[op](left, right)
        if op == "neg" and len(node.args) == 1:
            return _neg(_translate(node.args[0], env))
    raise RecipeProofError("unsupported or unknown source expression node")


def derivative_z(value: Expr) -> Expr:
    """Restricted symbolic first derivative with respect to the ordinary Z."""
    if value.op == "const":
        return ZERO
    if value.op == "z":
        return ONE
    if value.op in {"add"}:
        return _add(*(derivative_z(a) for a in value.args))
    if value.op == "mul":
        terms = []
        for index, argument in enumerate(value.args):
            terms.append(_mul(derivative_z(argument), *(value.args[:index] + value.args[index + 1 :])))
        return _add(*terms)
    if value.op == "neg":
        return _neg(derivative_z(value.args[0]))
    if value.op == "exp":
        return _mul(_exp(value.args[0]), derivative_z(value.args[0]))
    if value.op == "log":
        return _quotient(derivative_z(value.args[0]), value.args[0])
    if value.op == "quotient":
        numerator, denominator = value.args
        return _quotient(_sub(_mul(derivative_z(numerator), denominator),
                              _mul(numerator, derivative_z(denominator))),
                         _mul(denominator, denominator))
    raise RecipeProofError("cannot differentiate an unknown symbolic node")


def _function(tree: ast.Module, name: str) -> ast.FunctionDef | ast.AsyncFunctionDef:
    matches = [node for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
               and node.name == name]
    if len(matches) != 1:
        raise RecipeProofError(f"expected one top-level function {name!r}")
    return matches[0]


def _method(tree: ast.Module, class_name: str, method_name: str) -> ast.FunctionDef:
    classes = [node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == class_name]
    if len(classes) != 1:
        raise RecipeProofError(f"expected one class {class_name!r}")
    methods = [node for node in classes[0].body if isinstance(node, ast.FunctionDef) and node.name == method_name]
    if len(methods) != 1:
        raise RecipeProofError(f"expected one method {class_name}.{method_name}")
    return methods[0]


def _assigned(function: ast.AST, names: Sequence[str], *, fields: bool = False) -> dict[str, ast.AST]:
    found: dict[str, ast.AST] = {}
    wanted = set(names)
    for node in ast.walk(function):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for target in targets:
                key = None
                if isinstance(target, ast.Name) and not fields:
                    key = target.id
                elif (fields and isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name)
                      and target.value.id == "self"):
                    key = "self." + target.attr
                if key in wanted:
                    if key in found:
                        raise RecipeProofError("duplicate recipe assignment: " + key)
                    found[key] = node.value
    missing = wanted - set(found)
    if missing:
        raise RecipeProofError("missing recipe assignment: " + ", ".join(sorted(missing)))
    return found


def _parse(text: str) -> ast.Module:
    try:
        return ast.parse(text)
    except SyntaxError as exc:
        raise RecipeProofError("pinned source is not parseable Python") from exc


def _source_assignment(function: ast.AST, names: Sequence[str], *, fields: bool = False) -> dict[str, tuple[ast.AST, str]]:
    assignments = _assigned(function, names, fields=fields)
    return {key: (value, ast.unparse(value)) for key, value in assignments.items()}


def _extract_old_tw(function: ast.AST, env: Mapping[str, Expr]) -> tuple[Expr, str]:
    # The pinned old frame has a local Tw/Tw_old assignment.  The Dict fallback
    # keeps this extractor explicit if a source revision returns the recipe in
    # one literal mapping, without accepting arbitrary call results.
    assignments = []
    for name in ("Tw", "Tw_old"):
        try:
            assignments.extend(_source_assignment(function, (name,)).values())
        except RecipeProofError:
            pass
    if assignments:
        if len(assignments) != 1:
            raise RecipeProofError("old frame must expose exactly one Tw recipe")
        node, text = assignments[0]
    else:
        candidates = []
        for node in ast.walk(function):
            if isinstance(node, ast.Dict):
                for key, value in zip(node.keys, node.values):
                    if isinstance(key, ast.Constant) and key.value in {"Tw", "Tw_old"}:
                        candidates.append((value, ast.unparse(value)))
        if len(candidates) != 1:
            raise RecipeProofError("missing old Tw recipe")
        node, text = candidates[0]
    return _translate(node, env), text


def _prove_from_ast(old_control: str, old_frame: str, current_outer: str) -> dict[str, Any]:
    control_tree, frame_tree, outer_tree = map(_parse, (old_control, old_frame, current_outer))
    control_fn = _function(control_tree, "build_graph")
    control_assignments = _source_assignment(control_fn, ("P", "mu"))
    p = _translate(control_assignments["P"][0], {})
    old_mu = _translate(control_assignments["mu"][0], {"P": p})
    old_tw_fn = _function(frame_tree, "build")
    mu_binding = _assigned(old_tw_fn, ("mu",))["mu"]
    if ast.dump(mu_binding, include_attributes=False) != ast.dump(ast.parse("params['mu']", mode="eval").body, include_attributes=False):
        raise RecipeProofError("old Tw mu is not read from the frame parameter role")
    old_tw, old_tw_source = _extract_old_tw(old_tw_fn, {"mu": old_mu})
    old_lambda = _log(old_mu)

    outer_fn = _method(outer_tree, "CurrentOriginalTrueOuterC6C5Functions", "__init__")
    outer_assignments = _source_assignment(outer_fn, ("self.logmu", "self.mu", "self.Tw"), fields=True)
    current_lambda = _translate(outer_assignments["self.logmu"][0], {})
    current_mu = _translate(outer_assignments["self.mu"][0], {"self.logmu": current_lambda})
    current_tw = _translate(outer_assignments["self.Tw"][0],
                             {"self.logmu": current_lambda, "self.mu": current_mu})
    # Require the source recipe to actually be the expected -60*log(mu), not a
    # guessed value supplied by a caller.
    expected_old_tw = _mul(_const(-60), old_lambda)
    # ``_translate`` performs only the closed simplifications, including
    # log(exp(lambda)) = lambda when lambda is real.
    old_tw_simplified = old_tw
    if old_tw != expected_old_tw:
        # Source may contain an explicit log(exp(mu)) shape; normalize it once.
        normalized_old_tw = _translate(ast.parse(old_tw_source, mode="eval").body, {"mu": old_mu})
        if normalized_old_tw != expected_old_tw:
            raise RecipeProofError("old Tw is not the restricted -60*log(mu) recipe")
        old_tw_simplified = normalized_old_tw
    else:
        old_tw_simplified = old_tw

    if old_lambda != current_lambda:
        raise RecipeProofError("old and current lambda recipes are not structurally equal")
    if current_mu != _exp(current_lambda):
        raise RecipeProofError("current mu is not exp(current logmu)")
    if old_mu != _exp(old_lambda):
        raise RecipeProofError("old mu is not exp(old lambda)")
    if current_tw != _mul(_const(-60), current_lambda):
        raise RecipeProofError("current Tw is not the restricted -60*logmu recipe")
    old_tw_simplified = _mul(_const(-60), old_lambda)
    current_tw_simplified = _mul(_const(-60), current_lambda)
    old_dz, current_dz = derivative_z(old_tw_simplified), derivative_z(current_tw_simplified)
    if not (_is_zero(old_dz) and _is_zero(current_dz)):
        raise RecipeProofError("ordinary Z first derivative did not reduce to zero")

    return {
        "status": "formula_identity_proved",
        "recipe_scope": {
            "old_control_function": "build_graph",
            "old_frame_function": "build",
            "current_class": "CurrentOriginalTrueOuterC6C5Functions",
            "current_method": "__init__",
            "all_unknown_function_nodes_rejected": True,
            "source_values_evaluated": False,
        },
        "old": {
            "P": expr_text(p),
            "mu": expr_text(old_mu),
            "lambda": expr_text(old_lambda),
            "Tw_source": old_tw_source,
            "Tw": expr_text(old_tw_simplified),
        },
        "current": {
            "logmu_source": outer_assignments["self.logmu"][1],
            "mu_source": outer_assignments["self.mu"][1],
            "Tw_source": outer_assignments["self.Tw"][1],
            "lambda": expr_text(current_lambda),
            "mu": expr_text(current_mu),
            "Tw": expr_text(current_tw_simplified),
        },
        "identities": {
            "lambda_real": old_lambda.real and current_lambda.real,
            "log_exp_identity": "log(exp(lambda)) = lambda for real lambda",
            "restricted_log_exp_identity_proved": True,
            "old_mu_equals_exp_lambda": old_mu == _exp(old_lambda),
            "current_mu_equals_exp_lambda": current_mu == _exp(current_lambda),
            "old_Tw_equals_current_Tw": old_tw_simplified == current_tw_simplified,
            "ordinary_Z1_old": expr_text(old_dz),
            "ordinary_Z1_current": expr_text(current_dz),
            "ordinary_Z1_equals_zero": _is_zero(old_dz) and _is_zero(current_dz),
        },
        "provenance": {
            "ast_recipe_identity": True,
            "old_C1_adapter_mu_value_identity": False,
            "cross_graph_semantics_certified": False,
            "actual_old_current_source_identity": False,
            "physical_identity": False,
            "endpoint_identity": False,
            "error_identity": False,
            "global_identity": False,
            "reason": "AST recipes do not certify restored graph-node provenance or cross-graph semantics.",
        },
    }


def extract_tw_recipes(old_control_source: str, old_frame_source: str,
                       current_outer_source: str) -> dict[str, Any]:
    """Extract and prove the closed recipes from three pinned source texts."""
    return _prove_from_ast(old_control_source, old_frame_source, current_outer_source)


def _literal_path(data: Mapping[str, Any], path: Sequence[str]) -> Any:
    value: Any = data
    for key in path:
        if not isinstance(value, Mapping) or key not in value:
            raise RecipeProofError("protected reference path missing: " + ".".join(path))
        value = value[key]
    return value


def _opaque_equal(left: Any, right: Any) -> bool | None:
    # Node ids are opaque protected references.  They are compared only as
    # values; lists/mappings are not accepted and never become lookup keys.
    if isinstance(left, bool) or isinstance(right, bool):
        return None
    if not isinstance(left, int) or not isinstance(right, int) or left < 0 or right < 0:
        return None
    return left == right


def compare_protected_mu_refs(c3_archive: Mapping[str, Any], control_archive: Mapping[str, Any],
                              f_archive: Mapping[str, Any]) -> dict[str, Any]:
    """Compare the three named protected ``mu`` references without graph access."""
    paths = {
        "c3_limit_controls_C1_adapter_mu": (
            "actual_C3_limit_controls_and_Picard_functions", "exact_C2_limit_functions",
            "exact_C1_limit_adapter", "parameters", "mu"),
        "control_connection_mu": ("exact_current_integral_and_control_graph", "parameters", "mu"),
        "F_snapshot_geometry_mu": ("Rp_native_C5_pulse_source_bindings", "geometry", "mu"),
    }
    values: dict[str, Any] = {}
    missing: dict[str, str] = {}
    for label, path in paths.items():
        try:
            owner = c3_archive if label.startswith("c3_") else control_archive if label.startswith("control_") else f_archive
            values[label] = _literal_path(owner, path)
        except RecipeProofError as exc:
            missing[label] = str(exc)
    comparisons = {}
    labels = tuple(paths)
    if not missing:
        for left, right in ((labels[0], labels[1]), (labels[1], labels[2]), (labels[0], labels[2])):
            equal = _opaque_equal(values[left], values[right])
            comparisons[left + "_vs_" + right] = {
                "status": "equal_opaque_refs" if equal is True else
                          "different_opaque_refs" if equal is False else
                          "unresolved_noninteger_or_missing_reference"
            }
    status = "all_equal_opaque_refs" if not missing and comparisons and all(
        row["status"] == "equal_opaque_refs" for row in comparisons.values()) else (
        "unresolved_protected_reference" if missing or not comparisons else "opaque_refs_differ")
    return {
        "status": status,
        "paths": {label: ".".join(path) for label, path in paths.items()},
        "comparisons": comparisons,
        "values": values,
        "missing": missing,
        "protected_refs_used_as_graph_or_mapping_keys": False,
        "exact_graph_nodes_selected_or_traversed": False,
        "source_values_computed": False,
    }


def prove_tw_recipes_and_compare_refs(old_control_source: str, old_frame_source: str,
                                     current_outer_source: str, c3_archive: Mapping[str, Any],
                                     control_archive: Mapping[str, Any], f_archive: Mapping[str, Any]) -> dict[str, Any]:
    """Single-call public operation used by the one-entry runner."""
    proof = extract_tw_recipes(old_control_source, old_frame_source, current_outer_source)
    proof["protected_mu_refs"] = compare_protected_mu_refs(c3_archive, control_archive, f_archive)
    return proof


__all__ = [
    "SCHEMA", "RecipeProofError", "compare_protected_mu_refs", "derivative_z",
    "extract_tw_recipes", "expr_text", "prove_tw_recipes_and_compare_refs",
]
