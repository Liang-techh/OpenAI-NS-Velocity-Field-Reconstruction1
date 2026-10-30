"""Fail-closed AST audit of the live prescribed-to-outer Rc owner chain.

The caller supplies the seven fixed Python source texts below.  This module
only parses those texts and matches the constructor type guards, owner
assignments, the prescribed ``endpoints['Rm']`` construction, and the C6
checker ``is`` assertion.  It does not import a producer, read an archive or
receipt, inspect a graph, evaluate a value, or perform endpoint mathematics.

The result therefore proves only the live constructor source chain.  It does
not prove archive/receipt identity, cross-graph parameter or unit identity,
old/current geometry equality, endpoint cancellation, or numerical radius
values.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping
from typing import Any


SCHEMA = "rp_rc_owner_chain_audit_v1"
COMPONENT = "rp_rc_owner_chain_audit"


class OwnerChainAuditError(ValueError):
    """Raised when the caller's fixed source manifest is incomplete or changed."""


# Canonical relative paths are intentionally fixed.  A different key, even a
# basename or a Windows separator variant, is rejected instead of silently
# widening the source scope.
REQUIRED_SOURCE_PATHS = (
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions_check.py",
)


def _same(left: ast.AST, right: ast.AST) -> bool:
    return ast.dump(left, include_attributes=False) == ast.dump(right, include_attributes=False)


def _expression(text: str) -> ast.AST:
    return ast.parse(text, mode="eval").body


def _function(tree: ast.Module, name: str, path: str) -> ast.FunctionDef:
    matches = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == name]
    if len(matches) != 1:
        raise OwnerChainAuditError(f"{path}: missing or ambiguous top-level function {name}")
    return matches[0]


def _method(tree: ast.Module, class_name: str, method_name: str, path: str) -> ast.FunctionDef:
    classes = [node for node in tree.body if isinstance(node, ast.ClassDef) and node.name == class_name]
    if len(classes) != 1:
        raise OwnerChainAuditError(f"{path}: missing or ambiguous class {class_name}")
    methods = [
        node for node in classes[0].body
        if isinstance(node, ast.FunctionDef) and node.name == method_name
    ]
    if len(methods) != 1:
        raise OwnerChainAuditError(f"{path}: missing or ambiguous method {class_name}.{method_name}")
    return methods[0]


def _target_text(target: ast.AST) -> str:
    return ast.unparse(target)


def _assignments(function: ast.AST, target: str) -> list[ast.Assign]:
    matches = []
    for node in ast.walk(function):
        if not isinstance(node, ast.Assign):
            continue
        if any(_target_text(candidate) == target for candidate in node.targets):
            matches.append(node)
    return matches


def _assignment(function: ast.AST, target: str, expected: str, path: str) -> ast.Assign:
    matches = _assignments(function, target)
    if len(matches) != 1:
        raise OwnerChainAuditError(
            f"{path}: expected exactly one assignment to {target}, found {len(matches)}"
        )
    node = matches[0]
    if not _same(node.value, _expression(expected)):
        raise OwnerChainAuditError(
            f"{path}:{node.lineno}: assignment {target} does not match {expected}"
        )
    return node


def _qualified_attribute(text: str) -> ast.AST:
    return _expression(text)


def _type_guard(
    function: ast.AST,
    parameter: str,
    expected_type: str,
    path: str,
) -> ast.If:
    """Match the real ``type(x) is not ... or not x.acceptance_loaded`` guard."""

    expected_compare = ast.Compare(
        left=ast.Call(func=ast.Name(id="type", ctx=ast.Load()), args=[ast.Name(id=parameter, ctx=ast.Load())], keywords=[]),
        ops=[ast.IsNot()],
        comparators=[_qualified_attribute(expected_type)],
    )
    expected_loaded = ast.UnaryOp(
        op=ast.Not(),
        operand=ast.Attribute(
            value=ast.Name(id=parameter, ctx=ast.Load()),
            attr="acceptance_loaded",
            ctx=ast.Load(),
        ),
    )
    matches = []
    for node in ast.walk(function):
        if not isinstance(node, ast.If):
            continue
        test = node.test
        if not isinstance(test, ast.BoolOp) or not isinstance(test.op, ast.Or) or len(test.values) != 2:
            continue
        if not (_same(test.values[0], expected_compare) and _same(test.values[1], expected_loaded)):
            continue
        raises_type_error = any(
            isinstance(child, ast.Raise)
            and isinstance(child.exc, ast.Call)
            and isinstance(child.exc.func, ast.Name)
            and child.exc.func.id == "TypeError"
            for statement in node.body
            for child in ast.walk(statement)
        )
        if not raises_type_error:
            raise OwnerChainAuditError(f"{path}:{node.lineno}: type guard does not raise TypeError")
        matches.append(node)
    if len(matches) != 1:
        raise OwnerChainAuditError(
            f"{path}: expected exactly one checked type guard for {parameter}, found {len(matches)}"
        )
    return matches[0]


def _identity_assert(function: ast.AST, expected: str, path: str) -> ast.Assert:
    wanted = _expression(expected)
    matches = [
        node for node in ast.walk(function)
        if isinstance(node, ast.Assert) and _same(node.test, wanted)
    ]
    if len(matches) != 1:
        raise OwnerChainAuditError(
            f"{path}: expected exactly one identity assertion {expected}, found {len(matches)}"
        )
    return matches[0]


def _parse_sources(source_texts: Mapping[str, str]) -> dict[str, ast.Module]:
    if not isinstance(source_texts, Mapping):
        raise OwnerChainAuditError("source_texts must be a fixed path-to-source mapping")
    expected = set(REQUIRED_SOURCE_PATHS)
    actual = set(source_texts)
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    if missing:
        raise OwnerChainAuditError("missing fixed source text: " + ", ".join(missing))
    if extra:
        raise OwnerChainAuditError("unexpected source text path: " + ", ".join(extra))
    trees: dict[str, ast.Module] = {}
    for path in REQUIRED_SOURCE_PATHS:
        text = source_texts[path]
        if not isinstance(text, str):
            raise OwnerChainAuditError(f"{path}: source text must be str")
        try:
            trees[path] = ast.parse(text, filename=path, mode="exec")
        except SyntaxError as exc:
            raise OwnerChainAuditError(f"{path}:{exc.lineno}: syntax error") from exc
    return trees


def prove(source_texts: Mapping[str, str]) -> dict[str, Any]:
    """Prove only the source-level live owner chain, fail-closed.

    The mapping must contain exactly :data:`REQUIRED_SOURCE_PATHS`.  No path
    is opened here; callers provide the already selected source text.
    """

    trees = _parse_sources(source_texts)
    prescribed_path, loop_path, modified_path, c6_path, finite_path, outer_path, c6_check_path = REQUIRED_SOURCE_PATHS

    prescribed = _method(
        trees[prescribed_path],
        "CurrentOriginalPrescribedR100RmC6Functions",
        "__init__",
        prescribed_path,
    )
    loop = _method(
        trees[loop_path],
        "CurrentOriginalPrescribedLoopC5Functions",
        "__init__",
        loop_path,
    )
    modified = _method(
        trees[modified_path],
        "CurrentOriginalModifiedQ4PressureC3Functions",
        "__init__",
        modified_path,
    )
    c6 = _method(
        trees[c6_path],
        "CurrentOriginalTrueRmPatchC6Functions",
        "__init__",
        c6_path,
    )
    finite = _method(
        trees[finite_path],
        "CurrentOriginalTrueRmPatchFiniteNC5Functions",
        "__init__",
        finite_path,
    )
    outer = _method(
        trees[outer_path],
        "CurrentOriginalTrueOuterC6C5Functions",
        "__init__",
        outer_path,
    )
    c6_check = _function(trees[c6_check_path], "actual_source", c6_check_path)

    checks = {
        "prescribed_requires_checked_bridge": {
            "path": prescribed_path,
            "line": _type_guard(
                prescribed,
                "bridge_owner",
                "previous.CurrentOriginalPrescribedBridgeC6Functions",
                prescribed_path,
            ).lineno,
            "type": "previous.CurrentOriginalPrescribedBridgeC6Functions",
        },
        "prescribed_constructs_Rm_endpoint": {
            "path": prescribed_path,
            "line": _assignment(
                prescribed,
                "self.endpoints['Rm']",
                "self.end('post_restore', g.one)",
                prescribed_path,
            ).lineno,
        },
        "loop_requires_prescribed_background": {
            "path": loop_path,
            "line": _type_guard(
                loop,
                "background",
                "previous.CurrentOriginalPrescribedR100RmC6Functions",
                loop_path,
            ).lineno,
        },
        "loop_preserves_background_object": {
            "path": loop_path,
            "line": _assignment(loop, "self.background", "background", loop_path).lineno,
        },
        "modified_requires_loop": {
            "path": modified_path,
            "line": _type_guard(
                modified,
                "loop_owner",
                "previous.CurrentOriginalPrescribedLoopC5Functions",
                modified_path,
            ).lineno,
        },
        "modified_preserves_prescribed_background": {
            "path": modified_path,
            "line": _assignment(
                modified,
                "self.background",
                "loop_owner.background",
                modified_path,
            ).lineno,
        },
        "c6_requires_modified_owner": {
            "path": c6_path,
            "line": _type_guard(
                c6,
                "modified_owner",
                "previous.CurrentOriginalModifiedQ4PressureC3Functions",
                c6_path,
            ).lineno,
        },
        "c6_preserves_modified_background": {
            "path": c6_path,
            "line": _assignment(
                c6,
                "self.background",
                "modified_owner.background",
                c6_path,
            ).lineno,
        },
        "c6_uses_prescribed_Rm_endpoint_object": {
            "path": c6_path,
            "line": _assignment(
                c6,
                "self.inlet",
                "self.background.endpoints['Rm']",
                c6_path,
            ).lineno,
        },
        "finite_requires_C6_owner": {
            "path": finite_path,
            "line": _type_guard(
                finite,
                "leading",
                "previous.CurrentOriginalTrueRmPatchC6Functions",
                finite_path,
            ).lineno,
        },
        "finite_preserves_C6_owner": {
            "path": finite_path,
            "line": _assignment(finite, "self.leading", "leading", finite_path).lineno,
        },
        "outer_requires_finite_owner": {
            "path": outer_path,
            "line": _type_guard(
                outer,
                "patch",
                "previous.CurrentOriginalTrueRmPatchFiniteNC5Functions",
                outer_path,
            ).lineno,
        },
        "outer_preserves_C6_owner": {
            "path": outer_path,
            "line": _assignment(outer, "self.leading", "patch.leading", outer_path).lineno,
        },
        "outer_uses_leading_inlet_radius": {
            "path": outer_path,
            "line": _assignment(
                outer,
                "self.Rm",
                "self.leading.inlet['radius']",
                outer_path,
            ).lineno,
        },
        "c6_checker_asserts_live_identity": {
            "path": c6_check_path,
            "line": _identity_assert(
                c6_check,
                "owner.inlet is owner.background.endpoints['Rm']",
                c6_check_path,
            ).lineno,
        },
    }

    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": "live_constructor_source_chain_proved",
        "source_only": True,
        "checked_source_paths": list(REQUIRED_SOURCE_PATHS),
        "live_constructor_identity_proved": True,
        "constructor_chain": [
            "prescribed.endpoints['Rm']",
            "loop.background",
            "modified.background",
            "true_Rm_patch_C6.background",
            "true_Rm_patch_C6.inlet",
            "finite_N.leading",
            "true_outer.leading",
            "true_outer.leading.inlet['radius']",
        ],
        "checks": checks,
        "archive_receipts_read": False,
        "archive_receipt_identity_proved": False,
        "cross_graph_parameter_identity_proved": False,
        "common_units_proved": False,
        "old_current_geometry_equality_proved": False,
        "endpoint_cancellation_proved": False,
        "numeric_radius_computed": False,
        "limitations": [
            "No JSON or gzip archive/receipt is read or authenticated.",
            "No independently restored object or opaque reference is treated as a live owner identity.",
            "No cross-graph source-parameter or unit identity is claimed.",
            "No old phase_rc/current geometry equality or endpoint cancellation is claimed.",
            "No graph, integral, derivative, numeric radius, or mathematical value is evaluated.",
        ],
    }


audit = prove


__all__ = [
    "COMPONENT",
    "OwnerChainAuditError",
    "REQUIRED_SOURCE_PATHS",
    "SCHEMA",
    "audit",
    "prove",
]
