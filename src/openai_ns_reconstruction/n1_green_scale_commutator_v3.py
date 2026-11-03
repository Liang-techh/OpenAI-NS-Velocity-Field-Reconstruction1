"""Source-bound scale/Green commutator proof for the raw ``n=1`` provider.

This leaf reads source text and checks its AST.  It never imports the raw
provider, constructs a graph, evaluates a kernel, or consumes an interval.
The only asserted operator fact is the following conditional identity.  If
the six raw Green columns are the six literal rows used by ``green_apply``
and the source scale is diagonal with the same factor on each paired output
and input column, then ``S^-1 G S = G``.  The first-chart endpoint and the
fixed radial variables are explicit premises, so the conclusion about
``G_Z`` is a kernel statement and not a claim about the source vector.

The ``B1*S_Z`` term is retained in the transformed first ordinary-Z
operator.  In particular, this module does not replace the raw provider by
the normalized ``Kchi`` resolvent, does not evaluate a residual, and does
not infer a gauge owner from matching tuples or intervals.
"""

from __future__ import annotations

import ast
from collections.abc import Mapping
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable


SCHEMA = "n1_green_scale_commutator_v3"

VOLERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
REGULAR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py"
)
FIRST_RAW_SOURCE = "src/openai_ns_reconstruction/original_n1_first_raw_function.py"

# These three pins are the source pins supplied by the accepted parent cutoff
# report.  Each source is also checked by its literal AST contract; a caller
# label alone is never accepted as provenance.
VOLERRA_SHA256 = "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208"
REGULAR_SHA256 = "b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555"
FIRST_RAW_SHA256 = "f13f29b4d0652e8ed3b3347f19bab7eb1ea0e18fd4be0f50cf6d550a1aeb92c6"

ROWS = ("5", "6", "6", "4", "5", "6")
EXPECTED_FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
EXPECTED_SCALE = (
    "F0(Z)/epsilon",
    "epsilon",
    "epsilon",
    "epsilon*Pstar^2",
    "F0(Z)/epsilon",
    "epsilon",
)
CANONICAL_STATE = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")

_EXPECTED_GREEN_ARGS = ("g", "vectors", "chart", "coordinates", "parameters")
_EXPECTED_GREEN_OWNER = "<module>"
_REQUIRED_GREEN_FLAGS = (
    "actual_raw_regular_n1_Green_body",
    "endpoint_independent_of_Z",
    "complete_history_not_reset_at_core_seam",
)
_EXPECTED_PARENT_HASHES = {
    VOLERRA_SOURCE: VOLERRA_SHA256,
    REGULAR_SOURCE: REGULAR_SHA256,
    FIRST_RAW_SOURCE: FIRST_RAW_SHA256,
}


class SourceContractError(ValueError):
    """Raised when supplied source text cannot establish the AST contract."""


def _dump(node: ast.AST | None) -> str | None:
    if node is None:
        return None
    return ast.dump(node, annotate_fields=True, include_attributes=False)


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse(text: str, label: str) -> ast.Module:
    if not isinstance(text, str):
        raise TypeError(f"{label} must be source text")
    try:
        parsed = ast.parse(text, filename=label, mode="exec")
    except SyntaxError as exc:
        raise SourceContractError(f"invalid AST for {label}: {exc}") from exc
    if not isinstance(parsed, ast.Module):  # pragma: no cover - ast.parse contract
        raise SourceContractError(f"{label} did not produce a module AST")
    return parsed


def _expr(text: str) -> ast.AST:
    try:
        return ast.parse(text, mode="eval").body
    except SyntaxError as exc:  # pragma: no cover - internal expected fragments
        raise SourceContractError(f"internal expression is invalid: {text}") from exc


def _out_index_store_target() -> ast.AST:
    parsed = ast.parse("out[i]=None", mode="exec")
    return parsed.body[0].targets[0]


def _same(node: ast.AST | None, expected: str) -> bool:
    return _dump(node) == _dump(_expr(expected))


def _literal(node: ast.AST | None) -> Any:
    """Decode only constant/tuple/list AST nodes; never evaluate source."""

    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Tuple):
        return tuple(_literal(item) for item in node.elts)
    if isinstance(node, ast.List):
        return [_literal(item) for item in node.elts]
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = _literal(node.operand)
        if isinstance(value, (int, float, complex)) and not isinstance(value, bool):
            return -value if isinstance(node.op, ast.USub) else value
    return None


def _name(node: ast.AST | None) -> str | None:
    return node.id if isinstance(node, ast.Name) else None


def _dotted(node: ast.AST | None) -> str | None:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _dotted(node.value)
        return f"{parent}.{node.attr}" if parent else None
    return None


def _module_functions(tree: ast.Module, name: str) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    return [
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name
    ]


def _all_functions(tree: ast.Module, name: str) -> list[ast.FunctionDef | ast.AsyncFunctionDef]:
    return [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == name
    ]


def _module_assignment(tree: ast.Module, name: str) -> list[ast.AST]:
    values: list[ast.AST] = []
    for node in tree.body:
        if isinstance(node, ast.Assign):
            if any(isinstance(target, ast.Name) and target.id == name for target in node.targets):
                values.append(node.value)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == name:
            if node.value is not None:
                values.append(node.value)
    return values


def _direct_assignments(scope: ast.AST, name: str) -> list[ast.Assign | ast.AnnAssign]:
    body = getattr(scope, "body", ())
    result: list[ast.Assign | ast.AnnAssign] = []
    for node in body:
        if isinstance(node, ast.Assign) and any(
            isinstance(target, ast.Name) and target.id == name for target in node.targets
        ):
            result.append(node)
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.target.id == name:
            result.append(node)
    return result


def _assignment_value(scope: ast.AST, name: str) -> ast.AST | None:
    found = _direct_assignments(scope, name)
    return found[0].value if len(found) == 1 else None


def _plain_args(function: ast.FunctionDef | ast.AsyncFunctionDef) -> tuple[str, ...]:
    args = function.args
    if args.posonlyargs or args.vararg or args.kwonlyargs or args.kwarg or args.defaults or args.kw_defaults:
        return ()
    return tuple(item.arg for item in args.args)


def _call_path(node: ast.AST | None, path: str) -> ast.Call | None:
    if not isinstance(node, ast.Call) or _dotted(node.func) != path:
        return None
    return node


def _subscript_key(node: ast.AST | None, base: str) -> str | None:
    if not isinstance(node, ast.Subscript) or _dotted(node.value) != base:
        return None
    slc = node.slice
    if isinstance(slc, ast.Constant) and isinstance(slc.value, str):
        return slc.value
    return None


def _subscript_name(node: ast.AST | None, base: str, index_name: str) -> bool:
    if not isinstance(node, ast.Subscript) or _dotted(node.value) != base:
        return False
    return _name(node.slice) == index_name


def _vector_row_accesses(node: ast.AST, chart: str) -> list[ast.Subscript]:
    """Return nested ``vectors['chart'][row]`` accesses only.

    The outer subscript's ``value`` is itself a subscript, so a dotted-name
    check on the outer node would miss the actual row access.  Keeping this
    two-level AST check explicit also rejects a hidden literal row or a
    second chart input.
    """

    result: list[ast.Subscript] = []
    for item in ast.walk(node):
        if not isinstance(item, ast.Subscript) or not isinstance(item.value, ast.Subscript):
            continue
        if _subscript_key(item.value, "vectors") == chart and _name(item.slice) == "row":
            result.append(item)
    return result


def _all_vector_row_accesses(node: ast.AST) -> list[ast.Subscript]:
    """Return every nested ``vectors[chart][index]`` access in ``node``."""

    result: list[ast.Subscript] = []
    for item in ast.walk(node):
        if (
            isinstance(item, ast.Subscript)
            and isinstance(item.value, ast.Subscript)
            and _dotted(item.value.value) == "vectors"
        ):
            result.append(item)
    return result


def _comparison(node: ast.AST | None, left: str, op: type[ast.cmpop], right: str) -> bool:
    return (
        isinstance(node, ast.Compare)
        and _name(node.left) == left
        and len(node.ops) == 1
        and isinstance(node.ops[0], op)
        and len(node.comparators) == 1
        and isinstance(node.comparators[0], ast.Constant)
        and node.comparators[0].value == right
    )


def _rows_from_zip(node: ast.AST | None) -> tuple[str, ...] | None:
    if not isinstance(node, ast.Call) or _dotted(node.func) != "zip" or len(node.args) != 2 or node.keywords:
        return None
    rows = _literal(node.args[0])
    if not isinstance(rows, tuple) or not all(isinstance(item, str) for item in rows):
        return None
    if _name(node.args[1]) != "kernels":
        return None
    return tuple(rows)


def _for_target(node: ast.For, expected: tuple[str, ...]) -> bool:
    target = node.target
    if not isinstance(target, ast.Tuple) or len(target.elts) != len(expected):
        return False
    return tuple(_name(item) for item in target.elts) == expected


def _loop_rows(node: ast.For, expected_target: tuple[str, ...]) -> tuple[str, ...] | None:
    if not _for_target(node, expected_target):
        return None
    return _rows_from_zip(node.iter)


def _for_branch(branch: ast.If) -> list[ast.For]:
    return [node for node in branch.body if isinstance(node, ast.For)]


def _find_call(node: ast.AST, path: str) -> list[ast.Call]:
    return [item for item in ast.walk(node) if _call_path(item, path) is not None]  # type: ignore[arg-type]


def _keyword(call: ast.Call, name: str) -> ast.keyword | None:
    for keyword in call.keywords:
        if keyword.arg == name:
            return keyword
    return None


def _keyword_true(call: ast.Call, name: str) -> bool:
    keyword = _keyword(call, name)
    return keyword is not None and isinstance(keyword.value, ast.Constant) and keyword.value.value is True


def _keyword_zero(call: ast.Call, name: str) -> bool:
    keyword = _keyword(call, name)
    return keyword is not None and isinstance(keyword.value, ast.Constant) and keyword.value.value == 0


def _subscript_paths(node: ast.AST, root: str) -> list[ast.Subscript]:
    return [item for item in ast.walk(node) if isinstance(item, ast.Subscript) and _dotted(item.value) == root]


def _function_record(function: ast.FunctionDef | ast.AsyncFunctionDef, owner: str) -> dict[str, Any]:
    return {
        "function": function.name,
        "owner": owner,
        "signature": list(_plain_args(function)),
        "AST_sha256": hashlib.sha256(ast.dump(function, include_attributes=False).encode("utf-8")).hexdigest(),
        "source_executed": False,
    }


def _literal_constants(first_raw_tree: ast.Module, regular_tree: ast.Module) -> dict[str, Any]:
    first_fields_values = _module_assignment(first_raw_tree, "FIELDS")
    first_scale_values = _module_assignment(first_raw_tree, "SCALE_DIAGONAL")
    canonical_values = _module_assignment(first_raw_tree, "CANONICAL_STATE")
    regular_fields_values = _module_assignment(regular_tree, "FIELDS")
    fields_first = _literal(first_fields_values[0]) if len(first_fields_values) == 1 else None
    scale = _literal(first_scale_values[0]) if len(first_scale_values) == 1 else None
    canonical = _literal(canonical_values[0]) if len(canonical_values) == 1 else None
    fields_regular = _literal(regular_fields_values[0]) if len(regular_fields_values) == 1 else None
    checks = {
        "first_raw_FIELDS_module_literal": len(first_fields_values) == 1 and fields_first == EXPECTED_FIELDS,
        "regular_FIELDS_module_literal": len(regular_fields_values) == 1 and fields_regular == EXPECTED_FIELDS,
        "FIELDS_order_agrees_across_sources": fields_first == fields_regular == EXPECTED_FIELDS,
        "first_raw_SCALE_DIAGONAL_module_literal": len(first_scale_values) == 1 and scale == EXPECTED_SCALE,
        "first_raw_CANONICAL_STATE_module_literal": len(canonical_values) == 1 and canonical == CANONICAL_STATE,
    }
    return {
        "checks": checks,
        "fields": tuple(fields_first) if isinstance(fields_first, tuple) else (),
        "regular_fields": tuple(fields_regular) if isinstance(fields_regular, tuple) else (),
        "scale": tuple(scale) if isinstance(scale, tuple) else (),
        "canonical_state": tuple(canonical) if isinstance(canonical, tuple) else (),
        "AST": {
            "first_raw_FIELDS": _dump(first_fields_values[0]) if len(first_fields_values) == 1 else None,
            "first_raw_SCALE_DIAGONAL": _dump(first_scale_values[0]) if len(first_scale_values) == 1 else None,
            "regular_FIELDS": _dump(regular_fields_values[0]) if len(regular_fields_values) == 1 else None,
        },
    }


def _kernel_branch_checks(branch: ast.If, chart: str) -> dict[str, bool]:
    if chart == "core":
        expected = {
            "X": "g.unary('sqrt',r)",
            "sample": "g.mul(r,t,t)",
            "measure": "X",
            "K": "g.mul(c('1/2'),X,t,g.sub(g.one,g.mul(t,t)))",
            "log": "g.mul(X,t,g.unary('log',g.quotient(g.one,t,'positive t inside unit integral')))",
            "kernels": "(K,log,g.neg(K),g.one,g.mul(t,t,t),t)",
        }
        forbidden = {"Z", "z", "v", "vectors", "coordinates", "parameters"}
    else:
        expected = {
            "X": "g.mul(c(2),g.unary('exp',g.mul(c('1/2'),hb,s)))",
            "sample": "g.mul(c(4),t,t)",
            "measure": "c(2)",
            "distance": "g.mul(hb,s)",
            "tail": "g.mul(distance,g.unary('exprel',g.neg(distance)))",
            "K": "g.mul(t,g.add(g.sub(g.one,g.mul(t,t)),g.mul(t,t,tail)))",
            "log": "g.mul(c(2),t,g.add(g.mul(c('1/2'),hb,s),g.unary('log',g.quotient(g.one,t,'positive t inside unit integral'))))",
            "ratio": "g.mul(t,g.unary('exp',g.neg(g.mul(c('1/2'),hb,s))))",
            "kernels": "(K,log,g.neg(K),g.one,g.mul(ratio,ratio,ratio),ratio)",
        }
        forbidden = {"Z", "z", "v", "vectors", "coordinates", "parameters"}
    checks: dict[str, bool] = {}
    for name, fragment in expected.items():
        value = _assignment_value(branch, name)
        checks[f"{chart}_{name}_AST"] = value is not None and _same(value, fragment)
        if name == "kernels" and value is not None:
            names = {item.id for item in ast.walk(value) if isinstance(item, ast.Name)}
            checks[f"{chart}_kernel_has_no_Z_or_endpoint_name"] = not names.intersection(forbidden)
    # Each kernel tuple must be the literal six-column tuple, not a label in a
    # comment or a second hidden assignment in the branch.
    checks[f"{chart}_kernels_assignment_unique"] = len(_direct_assignments(branch, "kernels")) == 1
    return checks


def _integral_flags(loop: ast.For) -> dict[str, bool]:
    calls = [item for item in ast.walk(loop) if isinstance(item, ast.Call) and _dotted(item.func) == "controls.integral"]
    if len(calls) != 1:
        return {"integral_call_unique": False, **{f"integral_{name}": False for name in _REQUIRED_GREEN_FLAGS}}
    call = calls[0]
    checks: dict[str, bool] = {"integral_call_unique": True}
    checks.update({f"integral_{name}": _keyword_true(call, name) for name in _REQUIRED_GREEN_FLAGS})
    checks["integral_zero_endpoint_continuation"] = _keyword_zero(call, "zero_endpoint_integrand_continuation")
    checks["integral_signed_body_not_magnitude_cap"] = _keyword_true(call, "signed_body_not_magnitude_cap")
    order = _keyword(call, "ordinary_slow_Z_derivative_order")
    checks["integral_derivative_order_is_k"] = order is not None and _name(order.value) == "k"
    return checks


def _integral_wrapper_call_unique(loop: ast.For) -> bool:
    return sum(
        1
        for item in ast.walk(loop)
        if isinstance(item, ast.Call) and _dotted(item.func) == "integral"
    ) == 1


def _green_ast(volterra_tree: ast.Module) -> dict[str, Any]:
    top = _module_functions(volterra_tree, "green_apply")
    all_functions = _all_functions(volterra_tree, "green_apply")
    fields_assignments = _module_assignment(volterra_tree, "FIELDS")
    checks: dict[str, bool] = {
        "green_apply_top_level_unique": len(top) == 1 and len(all_functions) == 1,
        "green_apply_owner_is_module": len(top) == 1,
        "green_FIELDS_module_bound_to_regular_FIELDS": (
            len(fields_assignments) == 1
            and _dotted(fields_assignments[0]) == "regular.FIELDS"
        ),
    }
    if not top:
        return {
            "checks": checks,
            "passed": False,
            "function": None,
            "owner": None,
            "loops": (),
            "source_executed": False,
        }
    function = top[0]
    checks["green_apply_signature_exact"] = _plain_args(function) == _EXPECTED_GREEN_ARGS
    checks["green_apply_has_core_first_dispatch"] = any(
        isinstance(item, ast.If) and _comparison(item.test, "chart", ast.Eq, "core") for item in function.body
    )
    dispatch = next(
        (
            item
            for item in function.body
            if isinstance(item, ast.If) and _comparison(item.test, "chart", ast.Eq, "core")
        ),
        None,
    )
    if dispatch is None:
        return {
            "checks": checks,
            "passed": False,
            "function": _function_record(function, _EXPECTED_GREEN_OWNER),
            "owner": _EXPECTED_GREEN_OWNER,
            "loops": (),
            "source_executed": False,
        }
    first_dispatch = dispatch.orelse[0] if len(dispatch.orelse) == 1 and isinstance(dispatch.orelse[0], ast.If) else None
    checks["green_apply_first_dispatch_exact"] = first_dispatch is not None and _comparison(first_dispatch.test, "chart", ast.Eq, "first")
    if first_dispatch is None:
        return {
            "checks": checks,
            "passed": False,
            "function": _function_record(function, _EXPECTED_GREEN_OWNER),
            "owner": _EXPECTED_GREEN_OWNER,
            "loops": (),
            "source_executed": False,
        }

    checks.update(_kernel_branch_checks(dispatch, "core"))
    checks.update(_kernel_branch_checks(first_dispatch, "first"))
    # The core loop is the function-level loop after the chart-kernel
    # dispatch; the first-chart no-reset loop is nested in the first branch.
    # Count every ``for`` in the owner as well, so a hidden third loop cannot
    # satisfy the two visible mapping checks.
    core_loops = [node for node in function.body if isinstance(node, ast.For)]
    first_accumulation_branches = [
        node
        for node in function.body
        if isinstance(node, ast.If) and _comparison(node.test, "chart", ast.Eq, "first")
    ]
    checks["green_first_accumulation_branch_count_one"] = len(first_accumulation_branches) == 1
    first_loops = _for_branch(first_accumulation_branches[0]) if len(first_accumulation_branches) == 1 else []
    checks["green_for_loop_count_two"] = sum(
        isinstance(node, ast.For) for node in ast.walk(function)
    ) == 2
    checks["core_loop_count_one"] = len(core_loops) == 1
    checks["first_loop_count_one"] = len(first_loops) == 1
    loops: list[dict[str, Any]] = []

    out_assignments = _direct_assignments(function, "out")
    checks["out_initialized_once_as_empty_list"] = len(out_assignments) == 1 and _same(out_assignments[0].value, "[]")
    checks["return_is_zip_FIELDS_out"] = any(
        isinstance(node, ast.Return) and _same(node.value, "dict(zip(FIELDS,out))") for node in function.body
    )

    if len(core_loops) == 1:
        loop = core_loops[0]
        row_order = _loop_rows(loop, ("row", "kernel"))
        loop_checks: dict[str, bool] = {
            "target_row_kernel": row_order is not None,
            "rows_literal_exact": row_order == ROWS,
            "iter_is_direct_zip_rows_kernels": row_order is not None,
            "body_statement_count_two": len(loop.body) == 2,
            "body_uses_core_vector_row": False,
            "out_append_accumulation": False,
        }
        if len(loop.body) == 2:
            body_assign, append_expr = loop.body
            expected_body = "collar.at(g,vectors['core'][row],r,sample)"
            loop_checks["body_uses_core_vector_row"] = (
                isinstance(body_assign, ast.Assign)
                and len(body_assign.targets) == 1
                and _name(body_assign.targets[0]) == "body"
                and _same(body_assign.value, expected_body)
            )
            loop_checks["out_append_accumulation"] = (
                isinstance(append_expr, ast.Expr)
                and _same(append_expr.value, "out.append(integral(a.scale(body,g.mul(kernel,measure)),t,g.zero,g.one))")
            )
        loop_checks["no_hidden_cross_row_input"] = (
            len(_vector_row_accesses(loop, "core")) == 1
            and len(_all_vector_row_accesses(loop)) == 1
            and sum(
                isinstance(item, ast.Name) and item.id == "vectors"
                for item in ast.walk(loop)
            ) == 1
        )
        loop_checks["integral_wrapper_call_unique"] = _integral_wrapper_call_unique(loop)
        checks.update({f"core_loop_{key}": value for key, value in loop_checks.items()})
        loops.append({"chart": "core", "rows": list(row_order or ()), "checks": loop_checks, "lineno": loop.lineno})
    else:
        checks.update({f"core_loop_{key}": False for key in ("target_row_kernel", "rows_literal_exact", "iter_is_direct_zip_rows_kernels", "body_statement_count_two", "body_uses_core_vector_row", "out_append_accumulation", "no_hidden_cross_row_input", *[f for f in ("integral_call_unique", *[f"integral_{name}" for name in _REQUIRED_GREEN_FLAGS], "integral_zero_endpoint_continuation", "integral_signed_body_not_magnitude_cap", "integral_derivative_order_is_k")])})

    if len(first_loops) == 1:
        loop = first_loops[0]
        target_ok = (
            isinstance(loop.target, ast.Tuple)
            and len(loop.target.elts) == 2
            and _name(loop.target.elts[0]) == "i"
            and isinstance(loop.target.elts[1], ast.Tuple)
            and tuple(_name(item) for item in loop.target.elts[1].elts) == ("row", "kernel")
        )
        iter_value = loop.iter
        row_order = None
        if isinstance(iter_value, ast.Call) and _dotted(iter_value.func) == "enumerate" and len(iter_value.args) == 1:
            row_order = _rows_from_zip(iter_value.args[0])
        loop_checks = {
            "target_i_row_kernel": target_ok,
            "rows_literal_exact": row_order == ROWS,
            "iter_is_enumerated_direct_zip_rows_kernels": row_order is not None,
            "body_statement_count_two": len(loop.body) == 2,
            "body_uses_first_vector_row": False,
            "out_index_add_accumulation": False,
        }
        if len(loop.body) == 2:
            body_assign, output_assign = loop.body
            expected_body = "collar.at(g,vectors['first'][row],s,v)"
            loop_checks["body_uses_first_vector_row"] = (
                isinstance(body_assign, ast.Assign)
                and len(body_assign.targets) == 1
                and _name(body_assign.targets[0]) == "body"
                and _same(body_assign.value, expected_body)
            )
            loop_checks["out_index_add_accumulation"] = (
                isinstance(output_assign, ast.Assign)
                and len(output_assign.targets) == 1
                and _dump(output_assign.targets[0]) == _dump(_out_index_store_target())
                and _same(output_assign.value, "a.add(out[i],integral(a.scale(body,g.mul(kernel,measure)),v,g.zero,s))")
            )
        loop_checks["no_hidden_cross_row_input"] = (
            len(_vector_row_accesses(loop, "first")) == 1
            and len(_all_vector_row_accesses(loop)) == 1
            and sum(
                isinstance(item, ast.Name) and item.id == "vectors"
                for item in ast.walk(loop)
            ) == 1
        )
        loop_checks["integral_wrapper_call_unique"] = _integral_wrapper_call_unique(loop)
        checks.update({f"first_loop_{key}": value for key, value in loop_checks.items()})
        loops.append({"chart": "first", "rows": list(row_order or ()), "checks": loop_checks, "lineno": loop.lineno})
    else:
        checks.update({f"first_loop_{key}": False for key in ("target_i_row_kernel", "rows_literal_exact", "iter_is_enumerated_direct_zip_rows_kernels", "body_statement_count_two", "body_uses_first_vector_row", "out_index_add_accumulation", "no_hidden_cross_row_input", *[f for f in ("integral_call_unique", *[f"integral_{name}" for name in _REQUIRED_GREEN_FLAGS], "integral_zero_endpoint_continuation", "integral_signed_body_not_magnitude_cap", "integral_derivative_order_is_k")])})

    # The provider's integral helper itself must carry the ordinary-Z order,
    # source-row generator, and caller-owned endpoint arguments.  This is a
    # structural check in addition to the two loop checks above.
    integral_assignments = _direct_assignments(function, "integral")
    integral_lambda_ok = False
    if len(integral_assignments) == 1 and isinstance(integral_assignments[0].value, ast.Lambda):
        lam = integral_assignments[0].value
        integral_lambda_ok = tuple(arg.arg for arg in lam.args.args) == ("q", "var", "lo", "hi")
        calls = [item for item in ast.walk(lam.body) if isinstance(item, ast.Call) and _dotted(item.func) == "controls.integral"]
        integral_lambda_ok = integral_lambda_ok and len(calls) == 1
        if calls:
            control_call = calls[0]
            integral_lambda_ok = integral_lambda_ok and _name(control_call.args[1]) == "row"
            integral_lambda_ok = integral_lambda_ok and _name(control_call.args[3]) == "lo" and _name(control_call.args[4]) == "hi"
            generators = [item for item in ast.walk(lam.body) if isinstance(item, ast.comprehension)]
            if len(generators) != 1 or not isinstance(generators[0].target, ast.Tuple):
                integral_lambda_ok = False
            else:
                target = generators[0].target
                integral_lambda_ok = integral_lambda_ok and len(target.elts) == 2
                integral_lambda_ok = integral_lambda_ok and _name(target.elts[0]) == "k"
                integral_lambda_ok = integral_lambda_ok and _name(target.elts[1]) == "row"
            if generators:
                integral_lambda_ok = integral_lambda_ok and _same(generators[0].iter, "enumerate(c3.rows(q))")
    helper_flags = _integral_flags(integral_assignments[0].value) if integral_assignments and isinstance(integral_assignments[0].value, ast.Lambda) else {}
    checks["integral_helper_source_row_and_endpoints"] = integral_lambda_ok
    checks.update({f"integral_helper_{key}": value for key, value in helper_flags.items()})
    checks["exact_green_flags_present_in_both_integrals"] = all(
        checks.get(f"integral_helper_integral_{name}") is True
        for name in _REQUIRED_GREEN_FLAGS
    )
    checks["exact_rows_present_in_both_for_loops"] = checks.get("core_loop_rows_literal_exact") is True and checks.get("first_loop_rows_literal_exact") is True

    record = _function_record(function, _EXPECTED_GREEN_OWNER)
    return {
        "checks": checks,
        "passed": all(checks.values()),
        "function": record,
        "owner": _EXPECTED_GREEN_OWNER,
        "loops": loops,
        "source_executed": False,
    }


def _scale_records(fields: tuple[str, ...], scale: tuple[str, ...], canonical: tuple[str, ...]) -> dict[str, Any]:
    rows = ROWS
    records: list[dict[str, Any]] = []
    for index, (field, row) in enumerate(zip(fields, rows)):
        input_index = int(row) - 1 if row.isdigit() else -1
        output_scale = scale[index] if index < len(scale) else None
        input_scale = scale[input_index] if 0 <= input_index < len(scale) else None
        input_column = canonical[input_index] if 0 <= input_index < len(canonical) else None
        records.append(
            {
                "output_index": index + 1,
                "output_field": field,
                "input_row": row,
                "input_column": input_column,
                "output_scale": output_scale,
                "input_scale": input_scale,
                "same_scale": output_scale == input_scale,
            }
        )
    return {
        "rows": list(rows),
        "records": records,
        "all_output_input_scales_equal": bool(records) and all(item["same_scale"] for item in records),
    }


_PARENT_FIRST_RAW_FLAGS = (
    "actual_first_raw_source_definition_bound",
    "actual_first_raw_function_bound",
    "actual_first_raw_equation_bound",
    "actual_first_raw_core_entry_bound",
    "actual_first_raw_ordinary_Z0_to_3_bound",
    "actual_first_raw_scale_conjugacy",
    "actual_first_raw_all_word_convergence_bound",
    "source_only",
    "raw_preconditioned_resolvent_factorization",
)
_PARENT_FIRST_RAW_CHECK_FLAGS = (
    "parent",
    "source_pins",
    "graph",
    "equation",
    "scale",
    "core_entry",
    "ordinary_Z",
    "owner_definition",
    "raw_all_word_convergence",
)


def _direct_input_hashes(parent: Mapping[str, Any] | None) -> list[dict[str, str]]:
    """Read only the parent report's direct ``input_hashes`` mapping.

    Inherited/transitive receipts may contain same-named historical files and
    are deliberately not admissible as provenance for this source leaf.
    """

    found: list[dict[str, str]] = []
    child = parent.get("input_hashes") if isinstance(parent, Mapping) else None
    if isinstance(child, Mapping):
        for path, digest in child.items():
            if isinstance(path, str) and isinstance(digest, str):
                found.append({"path": path.replace("\\", "/"), "sha256": digest})
    return found


def _parent_checks(parent: Mapping[str, Any] | None, source_hashes: Mapping[str, str]) -> dict[str, Any]:
    if not isinstance(parent, Mapping):
        return {
            "parent_supplied": False,
            "parent_passed": False,
            "parent_checks_passed": False,
            "parent_actual_first_raw_checks": {key: False for key in _PARENT_FIRST_RAW_FLAGS},
            "parent_input_hashes_present": False,
            "parent_input_hashes_match": False,
            "parent_input_hashes": [],
        }
    pins = _direct_input_hashes(parent)
    checks = parent.get("checks")
    parent_checks_passed = isinstance(checks, Mapping) and checks.get("passed") is True
    first_raw = parent.get("actual_first_raw_n1_function")
    first_raw_checks = {
        key: isinstance(first_raw, Mapping) and first_raw.get(key) is True
        for key in _PARENT_FIRST_RAW_FLAGS
    }
    first_raw_nested_checks = first_raw.get("checks") if isinstance(first_raw, Mapping) else None
    first_raw_checks.update(
        {
            f"actual_first_raw_checks_{key}": (
                isinstance(first_raw_nested_checks, Mapping)
                and first_raw_nested_checks.get(key) is True
            )
            for key in _PARENT_FIRST_RAW_CHECK_FLAGS
        }
    )
    first_raw_checks["actual_first_raw_normalized_Kchi_kernel_substituted_is_false"] = (
        isinstance(first_raw, Mapping)
        and first_raw.get("normalized_Kchi_kernel_substituted") is False
    )
    first_raw_checks["actual_first_raw_status_source_bound_exact_raw_volterra_function"] = (
        isinstance(first_raw, Mapping)
        and first_raw.get("status") == "source_bound_exact_raw_volterra_function"
    )
    per_source: dict[str, bool] = {}
    for relative, expected in _EXPECTED_PARENT_HASHES.items():
        basename = relative.rsplit("/", 1)[-1]
        matches = [
            item["sha256"]
            for item in pins
            if item["path"] == relative or item["path"] == basename
        ]
        # A direct basename is admissible, but two direct entries (including
        # relative-plus-basename duplicates) are rejected rather than merged.
        per_source[relative] = (
            len(matches) == 1
            and matches[0] == expected
            and source_hashes.get(relative) == expected
        )
    parent_input_hashes_present = all(
        any(
            item["path"] == relative
            or item["path"] == relative.rsplit("/", 1)[-1]
            for item in pins
        )
        for relative in _EXPECTED_PARENT_HASHES
    )
    parent_input_hashes_match = all(per_source.values())
    parent_passed = (
        parent_checks_passed
        and all(first_raw_checks.values())
        and parent_input_hashes_present
        and parent_input_hashes_match
    )
    return {
        "parent_supplied": True,
        "parent_passed": parent_passed,
        "parent_checks_passed": parent_checks_passed,
        "parent_actual_first_raw_checks": first_raw_checks,
        "parent_input_hashes_present": parent_input_hashes_present,
        "parent_input_hashes_match": parent_input_hashes_match,
        "per_source": per_source,
        "parent_input_hashes": pins,
    }


def audit_source_texts(
    volterra_text: str,
    regular_text: str,
    first_raw_text: str,
    *,
    parent: Mapping[str, Any] | None = None,
    source_labels: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """Inspect the three supplied source texts and return a conditional proof.

    ``ast.parse`` is the only source interpretation step.  The function is
    intentionally usable by the synthetic fixture, which passes mutated text
    and expects the structural gates to fail closed.
    """

    labels = dict(source_labels or {})
    volterra_label = labels.get("volterra", VOLERRA_SOURCE)
    regular_label = labels.get("regular", REGULAR_SOURCE)
    first_raw_label = labels.get("first_raw", FIRST_RAW_SOURCE)
    volterra_tree = _parse(volterra_text, volterra_label)
    regular_tree = _parse(regular_text, regular_label)
    first_raw_tree = _parse(first_raw_text, first_raw_label)

    green = _green_ast(volterra_tree)
    constants = _literal_constants(first_raw_tree, regular_tree)
    fields = constants["fields"]
    scale = constants["scale"]
    canonical = constants["canonical_state"]
    scale_map = _scale_records(fields, scale, canonical) if fields and scale and canonical else {"rows": list(ROWS), "records": [], "all_output_input_scales_equal": False}

    source_hashes = {
        VOLERRA_SOURCE: _sha256_text(volterra_text),
        REGULAR_SOURCE: _sha256_text(regular_text),
        FIRST_RAW_SOURCE: _sha256_text(first_raw_text),
    }
    parent = dict(parent) if isinstance(parent, Mapping) else parent
    parent_result = _parent_checks(parent, source_hashes)
    checks: dict[str, bool] = dict(constants["checks"])
    checks.update({f"green_{key}": value for key, value in green["checks"].items()})
    checks["source_pin_volterra"] = source_hashes[VOLERRA_SOURCE] == VOLERRA_SHA256
    checks["source_pin_regular"] = source_hashes[REGULAR_SOURCE] == REGULAR_SHA256
    checks["source_pin_first_raw"] = source_hashes[FIRST_RAW_SOURCE] == FIRST_RAW_SHA256
    checks["source_pin_first_raw_AST_checked"] = all(
        constants["checks"].get(key) is True
        for key in (
            "first_raw_FIELDS_module_literal",
            "first_raw_SCALE_DIAGONAL_module_literal",
            "first_raw_CANONICAL_STATE_module_literal",
        )
    )
    checks["scale_mapping_by_literal_input_row"] = scale_map["all_output_input_scales_equal"]
    checks["parent_report_supplied"] = parent_result["parent_supplied"]
    checks["parent_report_passed"] = parent_result["parent_passed"]
    checks["parent_report_checks_passed"] = parent_result["parent_checks_passed"]
    checks.update(
        {
            f"parent_report_{key}": value
            for key, value in parent_result["parent_actual_first_raw_checks"].items()
        }
    )
    checks["parent_report_input_hashes_present"] = parent_result["parent_input_hashes_present"]
    checks["parent_report_input_hashes_match"] = parent_result["parent_input_hashes_match"]
    checks["parent_report_source_scale_AST_compatibility"] = (
        scale_map["all_output_input_scales_equal"]
        and parent_result["parent_actual_first_raw_checks"].get(
            "actual_first_raw_scale_conjugacy"
        ) is True
    )
    checks["A_minus_2chi_E51_scale_entries_equal"] = len(scale) >= 5 and scale[0] == scale[4]
    checks["G_kernel_Z_independence_AST"] = (
        checks.get("green_core_kernel_has_no_Z_or_endpoint_name") is True
        and checks.get("green_first_kernel_has_no_Z_or_endpoint_name") is True
    )
    checks["first_chart_additional_integral_is_no_reset_accumulation"] = checks.get(
        "green_first_loop_out_index_add_accumulation"
    ) is True

    proof = {
        "identity": "S^-1*G*S = G",
        "entrywise_argument": "G[i,j] != 0 only for j=ROWS[i]-1, and SCALE_DIAGONAL[i]=SCALE_DIAGONAL[j]; hence S_i^-1 G[i,j] S_j=G[i,j]",
        "green_kernel_Z_derivative": "G_Z=0 under fixed x, fixed hb, and Z-independent endpoints/phase bounds",
        "A_argument": "conditional algebraic reuse: if A=-2*chi*E51, then S1=S5=F0(Z)/epsilon gives S^-1*A*S=A; the external A support is not reverified by this leaf",
        "conditional_raw_residual_transport": {
            "premise_same_function_owner": "the raw and hatted six-component functions are the same owner",
            "premise_Qhat": "Qhat=S^-1*Q*S, with B1*S_Z retained in the ordinary-Z connection",
            "value_identity": "d_hat=S^-1*d_raw",
            "fixed_x_ordinary_Z_identity": "d_hat_Z=S^-1*d_raw_Z-Gamma*d_hat, Gamma=S^-1*S_Z",
            "source_owner_established_here": False,
        },
        "ordinary_Z_connection": {
            "S_Z": "diag(F0_Z/epsilon,0,0,0,F0_Z/epsilon,0)",
            "Gamma": "S^-1*S_Z=diag(F0_Z/F0,0,0,0,F0_Z/F0,0)",
            "B0_transformed": "S^-1*(B0*S+B1*S_Z)",
            "B1_transformed": "S^-1*B1*S",
            "g_transformed": "S^-1*g",
            "B1_S_Z_retained": True,
            "A_and_B1_source_reverified_here": False,
        },
        "scale_rows": scale_map,
        "parent_raw_scale_binding": {
            "parent_actual_first_raw_scale_conjugacy": parent_result[
                "parent_actual_first_raw_checks"
            ].get("actual_first_raw_scale_conjugacy") is True,
            "parent_raw_preconditioned_resolvent_factorization": parent_result[
                "parent_actual_first_raw_checks"
            ].get("raw_preconditioned_resolvent_factorization") is True,
            "parent_normalized_Kchi_kernel_substituted_is_false": parent_result[
                "parent_actual_first_raw_checks"
            ].get("actual_first_raw_normalized_Kchi_kernel_substituted_is_false") is True,
            "source_AST_scale_rows_match": scale_map["all_output_input_scales_equal"],
            "compatibility_check": checks["parent_report_source_scale_AST_compatibility"],
        },
    }
    premises = {
        "F0_epsilon_Pstar_radial_independence": "F0, epsilon, and Pstar do not depend on the radial input (fixed x); F0 and epsilon are nonzero and Pstar is admissible for S",
        "fixed_x": True,
        "fixed_hb": True,
        "endpoint_independent_of_Z": True,
        "source_rows_are_the_literal_green_rows": list(ROWS),
        "input_field_order_is_literal": list(fields),
        "volterra_FIELDS_bound_to_regular_FIELDS": green["checks"].get(
            "green_FIELDS_module_bound_to_regular_FIELDS"
        ) is True,
        "parent_raw_function_bound": parent_result["parent_passed"],
        "parent_raw_scale_conjugacy_bound": parent_result[
            "parent_actual_first_raw_checks"
        ].get("actual_first_raw_scale_conjugacy") is True,
        "parent_raw_preconditioned_resolvent_factorization_bound": parent_result[
            "parent_actual_first_raw_checks"
        ].get("raw_preconditioned_resolvent_factorization") is True,
        "parent_normalized_Kchi_kernel_not_substituted": parent_result[
            "parent_actual_first_raw_checks"
        ].get("actual_first_raw_normalized_Kchi_kernel_substituted_is_false") is True,
        "source_scale_AST_matches_parent_scale_binding": checks[
            "parent_report_source_scale_AST_compatibility"
        ],
        "no_hidden_cross_row_input": all(
            green["checks"].get(key) is True
            for key in ("core_loop_no_hidden_cross_row_input", "first_loop_no_hidden_cross_row_input")
        ),
    }
    limitations = {
        "complete_Kchi_resolvent_provided": False,
        "raw_numerical_residual_evaluated": False,
        "gauge_owner_proved": False,
        "historical_error_removed": False,
        "historical_error_added": False,
        "interval_or_tuple_overlap_used_as_function_identity": False,
        "source_provider_executed": False,
        "full_B0_B1_operator_conjugacy_beyond_the_proved_connection": False,
        "A_and_B1_support_reverified_here": False,
    }
    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "source_AST": {
            "green_apply": green,
            "first_raw_literals": constants["AST"],
            "regular_FIELDS_literal": constants["AST"]["regular_FIELDS"],
        },
        "function_owner": {
            "green_apply": green.get("owner"),
            "expected": _EXPECTED_GREEN_OWNER,
            "owner_checked_from_AST": True,
        },
        "source_paths": {"volterra": volterra_label, "regular": regular_label, "first_raw": first_raw_label},
        "source_hashes": source_hashes,
        "parent_report": parent_result,
        "fields": list(fields),
        "scale_diagonal": list(scale),
        "rows": list(ROWS),
        "scale_mapping": scale_map,
        "function_level_proof": proof,
        "premises": premises,
        "limitations": limitations,
        "source_executed": False,
        "archive_source_executed": False,
        "old_modules_called": False,
        "numeric_evaluation_performed": False,
    }


def _resolve_root(root: str | Path | None) -> Path:
    starts: list[Path] = []
    if root is not None:
        starts.append(Path(root))
    starts.extend((Path.cwd(), Path(__file__).resolve().parents[2]))
    for start in starts:
        if start.is_file():
            start = start.parent
        for candidate in (start, *start.parents):
            if (candidate / "src").is_dir() and (candidate / "experiments" / "root_st073").is_dir():
                return candidate
    raise FileNotFoundError("st073 repository root is required")


def _load_parent(
    parent: Mapping[str, Any] | str | Path | None,
    base: Path | None = None,
) -> Mapping[str, Any] | None:
    if parent is None:
        return None
    if isinstance(parent, Mapping):
        return parent
    path = Path(parent)
    if not path.is_absolute() and base is not None:
        rooted = base / path
        if rooted.is_file():
            path = rooted
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceContractError(f"parent report is not readable JSON: {path}") from exc


def build_source_bound_proof(
    repo: str | Path | None = None,
    parent: Mapping[str, Any] | str | Path | None = None,
) -> dict[str, Any]:
    """Read the pinned files, parse their ASTs, and return the proof packet."""

    root = _resolve_root(repo)
    paths = {
        "volterra": root / VOLERRA_SOURCE,
        "regular": root / REGULAR_SOURCE,
        "first_raw": root / FIRST_RAW_SOURCE,
    }
    texts: dict[str, str] = {}
    for role, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing source text for {role}: {path}")
        texts[role] = path.read_text(encoding="utf-8")
    parent_mapping = _load_parent(parent, root)
    result = audit_source_texts(
        texts["volterra"],
        texts["regular"],
        texts["first_raw"],
        parent=parent_mapping,
        source_labels={key: str(value) for key, value in paths.items()},
    )
    # Pin file bytes, not the platform-normalized text representation.  The
    # AST proof still consumes text, while provenance hashes remain identical
    # to the parent report even on a newline-translating filesystem.
    byte_hashes = {
        VOLERRA_SOURCE: _sha256_file(paths["volterra"]),
        REGULAR_SOURCE: _sha256_file(paths["regular"]),
        FIRST_RAW_SOURCE: _sha256_file(paths["first_raw"]),
    }
    result["source_hashes"] = byte_hashes
    result["checks"]["source_pin_volterra"] = byte_hashes[VOLERRA_SOURCE] == VOLERRA_SHA256
    result["checks"]["source_pin_regular"] = byte_hashes[REGULAR_SOURCE] == REGULAR_SHA256
    result["checks"]["source_pin_first_raw"] = byte_hashes[FIRST_RAW_SOURCE] == FIRST_RAW_SHA256
    parent_result = _parent_checks(parent_mapping, byte_hashes)
    result["parent_report"] = parent_result
    result["checks"]["parent_report_supplied"] = parent_result["parent_supplied"]
    result["checks"]["parent_report_passed"] = parent_result["parent_passed"]
    result["checks"]["parent_report_checks_passed"] = parent_result["parent_checks_passed"]
    result["checks"].update(
        {
            f"parent_report_{key}": value
            for key, value in parent_result["parent_actual_first_raw_checks"].items()
        }
    )
    result["checks"]["parent_report_input_hashes_present"] = parent_result["parent_input_hashes_present"]
    result["checks"]["parent_report_input_hashes_match"] = parent_result["parent_input_hashes_match"]
    result["passed"] = all(result["checks"].values())
    result["repository_root"] = str(root)
    return result


def prove(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Compatibility wrapper for source-text or repository-root callers."""

    if len(args) >= 3 and all(isinstance(value, str) for value in args[:3]):
        return audit_source_texts(args[0], args[1], args[2], **kwargs)
    return build_source_bound_proof(*args, **kwargs)


build = build_source_bound_proof
run = build_source_bound_proof
audit = audit_source_texts
verify = build_source_bound_proof


__all__ = [
    "SCHEMA",
    "ROWS",
    "EXPECTED_FIELDS",
    "EXPECTED_SCALE",
    "VOLERRA_SOURCE",
    "REGULAR_SOURCE",
    "FIRST_RAW_SOURCE",
    "SourceContractError",
    "audit_source_texts",
    "build_source_bound_proof",
    "build",
    "prove",
    "audit",
    "verify",
    "run",
]



