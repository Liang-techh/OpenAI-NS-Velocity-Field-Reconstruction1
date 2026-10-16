"""Bounded source-semantic guards for the prescribed native derivative.

The accepted C5 graph differentiates with the fast phase held fixed and adds
the admitted phase derivatives afterwards.  This module checks that source
contract and its selected reference/restoration bindings.  It deliberately
does not identify two endpoint derivative graphs, evaluate a cap, replay an
ancestor, or claim the physical inlet seam theorem.
"""

from __future__ import annotations

import hashlib
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping

from .original_reference_germ_guards import (
    _function_evidence,
    _handle_index,
    _node,
    _resolve,
)


_SOURCE_FILES = {
    "modified": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
    "loop": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
    "outer": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py",
    "bridge": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py",
}
_CHARTS = ("reference", "restoration")
_ORDINARY_Z_KEYS = (0, 1)
_REQUIRED_GUARDS = (
    "source_pins",
    "accepted_archives",
    "native_derivative",
    "radial_and_pressure",
    "generic_roots",
    "histories",
    "flat_cutoff",
    "parameters",
    "phase_global_R_Z0",
    "phase_source_recipe",
    "normalized_ODE_and_axial_recipe",
    "leading_germs",
)


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping) and isinstance(value.get(key), Mapping):
        return value[key]
    return _mapping(value)


def _repo_root(root: Any) -> Path:
    if root is not None:
        path = Path(root)
        if path.is_file():
            path = path.parent
        if (path / "experiments").is_dir() or (path / "src").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _strict_pin(root: Path, relative: str, inherited: Any) -> dict[str, Any]:
    """Require the inherited digest; a newly computed digest is not a pin."""
    row = _mapping(inherited)
    canonical = Path(relative).as_posix()
    path = _resolve(root, relative)
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if path is not None else None
    supplied = str(row.get("path", "")).replace("\\", "/")
    expected = row.get("expected_sha256")
    passed = bool(
        supplied == canonical
        and path is not None
        and row.get("exists") is True
        and row.get("pinned") is True
        and isinstance(expected, str)
        and expected == actual
        and row.get("sha256") == actual
    )
    return {
        "passed": passed,
        "path": supplied,
        "canonical_path": canonical,
        "exists": path is not None,
        "sha256": actual,
        "expected_sha256": expected,
        "inherited_pin": row,
    }


def _guard_row(source_guards: Mapping[str, Any], name: str) -> bool:
    checks = _mapping(source_guards.get("checks"))
    row = _mapping(checks.get(name))
    if "passed" in row:
        return row.get("passed") is True
    row = _mapping(source_guards.get(name))
    return row.get("passed") is True


def _handle(value: Any) -> int | None:
    return _handle_index(value)


def _int_row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    return _node(nodes, value)


def _contains(nodes: Any, root: Any, target: Any, *, limit: int = 20000) -> bool:
    """Bounded graph reachability for a selected source handle."""
    wanted = _handle(target)
    start = _handle(root)
    if wanted is None or start is None:
        return False
    stack = [start]
    seen: set[int] = set()
    while stack and len(seen) <= limit:
        index = stack.pop()
        if index in seen:
            continue
        if index == wanted:
            return True
        seen.add(index)
        row = _int_row(nodes, index)
        if row is None:
            continue
        for key in (
            "arguments", "argument", "numerator", "denominator", "integrand",
            "lower", "upper", "expression", "value", "active_body",
            "flat_value", "phase_function", "target_phase", "inverse_angle",
        ):
            value = row.get(key)
            if isinstance(value, (list, tuple)):
                stack.extend(i for i in value if _handle(i) is not None)
            elif _handle(value) is not None:
                stack.append(_handle(value))
    return False


def _poly_clean(value: Mapping[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    return {key: coefficient for key, coefficient in value.items() if coefficient}


def _poly_add(*values: Mapping[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    result: dict[tuple[str, ...], Fraction] = {}
    for value in values:
        for key, coefficient in value.items():
            result[key] = result.get(key, Fraction(0)) + coefficient
    return _poly_clean(result)


def _poly_mul(left: Mapping[tuple[str, ...], Fraction], right: Mapping[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    result: dict[tuple[str, ...], Fraction] = {}
    for lkey, lvalue in left.items():
        for rkey, rvalue in right.items():
            key = tuple(sorted(lkey + rkey))
            result[key] = result.get(key, Fraction(0)) + lvalue * rvalue
    return _poly_clean(result)


def _poly_neg(value: Mapping[tuple[str, ...], Fraction]) -> dict[tuple[str, ...], Fraction]:
    return _poly_clean({key: -coefficient for key, coefficient in value.items()})


def _poly_serial(value: Mapping[tuple[str, ...], Fraction] | None) -> list[list[Any]] | None:
    if value is None:
        return None
    return [[list(key), coefficient.numerator, coefficient.denominator]
            for key, coefficient in sorted(value.items())]


def _native_expr(
    nodes: Any,
    index: Any,
    coordinate: int,
    memo: dict[int, tuple[dict[tuple[str, ...], Fraction], dict[tuple[str, ...], Fraction]] | None],
) -> tuple[dict[tuple[str, ...], Fraction], dict[tuple[str, ...], Fraction]] | None:
    """Small exact polynomial evaluator for the selected radius subgraph.

    Source parameters and other native-independent leaves are opaque symbols;
    only sums, products, negatives, exact rationals and independent quotients
    are normalized.  A coordinate-dependent unissued operation returns None
    instead of being treated as an identity.
    """
    i = _handle(index)
    if i is None:
        return None
    if i in memo:
        return memo[i]
    row = _int_row(nodes, i)
    if row is None:
        memo[i] = None
        return None
    op = row.get("operation")
    result: tuple[dict[tuple[str, ...], Fraction], dict[tuple[str, ...], Fraction]] | None
    if op == "exact_rational":
        try:
            value = Fraction(row["numerator"], row["denominator"])
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            result = None
        else:
            result = ({(): value}, {})
    elif op == "bound_variable":
        expression = {(f"node:{i}",): Fraction(1)}
        result = (expression, {(): Fraction(1)} if i == coordinate else {})
    elif op in {
        "mathematical_pi", "current_original_source_parameter", "original_source_parameter",
        "exact_positive_integer_expression", "original_parameter", "bound_source_parameter",
    }:
        result = ({(f"node:{i}",): Fraction(1)}, {})
    elif op == "sum":
        children = [_native_expr(nodes, child, coordinate, memo) for child in row.get("arguments", ())]
        result = None if any(child is None for child in children) else (
            _poly_add(*(child[0] for child in children if child is not None)),
            _poly_add(*(child[1] for child in children if child is not None)),
        )
    elif op == "product":
        children = [_native_expr(nodes, child, coordinate, memo) for child in row.get("arguments", ())]
        if any(child is None for child in children):
            result = None
        else:
            expression: dict[tuple[str, ...], Fraction] = {(): Fraction(1)}
            derivative: dict[tuple[str, ...], Fraction] = {}
            for position, child in enumerate(children):
                assert child is not None
                # Product rule with the current factor differentiated.
                term = child[1]
                for other_position, other in enumerate(children):
                    if other_position != position and other is not None:
                        term = _poly_mul(term, other[0])
                derivative = _poly_add(derivative, term)
                expression = _poly_mul(expression, child[0])
            result = (expression, derivative)
    elif op == "negative":
        child = _native_expr(nodes, row.get("argument"), coordinate, memo)
        result = None if child is None else (_poly_neg(child[0]), _poly_neg(child[1]))
    elif op == "positive_quotient":
        numerator = _native_expr(nodes, row.get("numerator"), coordinate, memo)
        denominator = _native_expr(nodes, row.get("denominator"), coordinate, memo)
        # No quotient in the selected log-radius affine term is allowed to
        # depend on the native coordinate.  Independent ratios remain opaque.
        if numerator is None or denominator is None or numerator[1] or denominator[1]:
            result = None
        else:
            result = ({(f"node:{i}",): Fraction(1)}, {})
    elif op == "analytic_unary":
        child = _native_expr(nodes, row.get("argument"), coordinate, memo)
        result = None if child is None or child[1] else ({(f"node:{i}",): Fraction(1)}, {})
    else:
        result = None
    memo[i] = result
    return result


def _phase_binding(nodes: Any, chart: Mapping[str, Any], common_n: Any) -> dict[str, Any]:
    phase_handle = _handle(chart.get("actual_phase"))
    offset_handle = _handle(chart.get("global_offset"))
    coordinate = _handle(chart.get("coordinate"))
    phase = _int_row(nodes, phase_handle)
    argument = _handle(phase.get("argument")) if phase else None
    product = _int_row(nodes, argument)
    args = product.get("arguments", ()) if product else ()
    if not isinstance(args, (list, tuple)):
        args = ()
    row_passed = bool(
        phase_handle is not None
        and phase
        and phase.get("operation") == "analytic_unary"
        and phase.get("name") == "fractional_part"
        and product
        and product.get("operation") == "product"
        and len(args) == 2
        and sum(_handle(x) == _handle(common_n) for x in args) == 1
        and sum(_handle(x) == offset_handle for x in args) == 1
        and phase_handle not in {coordinate, _handle(chart.get("phase_center")), _handle(chart.get("phase_origin"))}
    )
    return {
        "passed": row_passed,
        "coordinate": coordinate,
        "global_offset": offset_handle,
        "actual_phase": phase_handle,
        "phase_node": phase,
        "product_node": product,
        "common_N": _handle(common_n),
    }


def _positive_jacobian_premise(
    nodes: Any,
    chart: Mapping[str, Any],
    outer: Mapping[str, Any],
) -> dict[str, Any]:
    jac = _handle(chart.get("physical_measure"))
    row = _int_row(nodes, jac)
    geometry = _mapping(outer.get("geometry_certificate"))
    # The source geometry certificate supplies positivity of the reference
    # length; restoration uses exact J = 1.  The AST guard below binds the
    # same J into every native seed and divides by it once.
    if chart.get("physical_measure") == 1:
        positive = bool(row and row.get("operation") == "exact_rational"
                        and row.get("numerator") == 1 and row.get("denominator") == 1)
        source = "exact_rational_one"
    else:
        positive = bool(
            row
            and geometry.get("directed_reference_gap_positive") is True
            and row.get("operation") in {"sum", "product", "positive_quotient", "analytic_unary", "current_original_source_parameter"}
        )
        source = "outer_geometry_certificate_directed_reference_gap_positive"
    return {
        "passed": positive,
        "handle": jac,
        "node_operation": row.get("operation") if row else None,
        "source": source,
        "geometry_certificate": {
            "directed_reference_gap_positive": geometry.get("directed_reference_gap_positive"),
            "directed_post_length_positive": geometry.get("directed_post_length_positive"),
        },
    }


def _radius_log_jacobian_binding(
    nodes: Any,
    chart: Mapping[str, Any],
    common_n: Any = None,
    phase_origin: Any = None,
) -> dict[str, Any]:
    """Check the source radius map carries the selected native Jacobian.

    The reference radius has ``log(R)`` affine in ``J * coordinate``.  The
    restoration radius has ``log(R)`` affine in its coordinate with ``J=1``.
    This is a structural source check; it does not differentiate an evaluated
    endpoint value.
    """
    coordinate = _handle(chart.get("coordinate"))
    jacobian = _handle(chart.get("physical_measure"))
    if coordinate is None or jacobian is None:
        return {"passed": False, "reason": "missing native coordinate or physical measure"}
    radius = _int_row(nodes, chart.get("radius"))
    if not radius or radius.get("operation") != "product":
        return {"passed": False, "reason": "radius is not a product with an exponential source"}
    exp_rows = []
    for value in radius.get("arguments", ()):
        row = _int_row(nodes, value)
        if row and row.get("operation") == "analytic_unary" and row.get("name") == "exp":
            exp_rows.append((_handle(value), row))
    if len(exp_rows) != 1:
        return {"passed": False, "reason": "radius does not have one exponential factor"}
    exp_handle, exp_row = exp_rows[0]
    log_argument = _handle(exp_row.get("argument"))
    log_row = _int_row(nodes, log_argument)
    if not log_row or log_row.get("operation") != "sum":
        return {"passed": False, "exp": exp_handle, "log_argument": log_argument}
    args = tuple(log_row.get("arguments", ()))
    memo: dict[int, tuple[dict[tuple[str, ...], Fraction], dict[tuple[str, ...], Fraction]] | None] = {}
    jacobian_expr = _native_expr(nodes, jacobian, coordinate, memo)
    log_expr = _native_expr(nodes, log_argument, coordinate, memo)
    if jacobian_expr is None or log_expr is None:
        return {
            "passed": False,
            "radius": _handle(chart.get("radius")),
            "coordinate": coordinate,
            "physical_measure": jacobian,
            "exp": exp_handle,
            "log_argument": log_argument,
            "reason": "selected affine radius subgraph contains an unissued coordinate-dependent operation",
        }
    radius_factor_rows = []
    for value in radius.get("arguments", ()):
        factor = _native_expr(nodes, value, coordinate, memo)
        factor_row = _int_row(nodes, value)
        radius_factor_rows.append({
            "handle": _handle(value),
            "operation": factor_row.get("operation") if factor_row else None,
            "native_derivative_zero": bool(factor is not None and not factor[1]),
        })
    other_radius_factors_independent = all(
        row["native_derivative_zero"] for row in radius_factor_rows
        if row["handle"] != exp_handle
    )
    native_J_constant = not jacobian_expr[1]
    log_derivative_equals_J = log_expr[1] == jacobian_expr[0]

    # Bind global_offset to log(R/Ra)-phase_origin using exact handles.  This
    # supplies the missing native chain factor for phase=N*global_offset.
    offset = _handle(chart.get("global_offset"))
    offset_row = _int_row(nodes, offset)
    log_ratio = None
    ra = None
    origin_term = False
    offset_terms = tuple(offset_row.get('arguments', ())) if offset_row else ()
    matching_log_terms = []
    matching_origin_terms = []
    if offset_row and offset_row.get("operation") == "sum":
        for value in offset_row.get("arguments", ()):
            candidate = _int_row(nodes, value)
            if candidate and candidate.get("operation") == "analytic_unary" and candidate.get("name") == "log":
                log_ratio = _handle(value)
                matching_log_terms.append(_handle(value))
                ratio = _int_row(nodes, candidate.get("argument"))
                if ratio and ratio.get("operation") == "positive_quotient":
                    if _handle(ratio.get("numerator")) == _handle(chart.get("radius")):
                        ra = _handle(ratio.get("denominator"))
            if candidate and candidate.get("operation") == "negative" and _handle(candidate.get("argument")) == _handle(phase_origin):
                origin_term = True
                matching_origin_terms.append(_handle(value))
    ra_expr = _native_expr(nodes, ra, coordinate, memo) if ra is not None else None
    origin_expr = _native_expr(nodes, phase_origin, coordinate, memo) if phase_origin is not None else None
    exact_offset_two_terms = bool(len(offset_terms) == 2 and len(matching_log_terms) == 1 and len(matching_origin_terms) == 1)
    radius_ratio_bound = bool(exact_offset_two_terms and log_ratio is not None and ra_expr is not None and not ra_expr[1])
    phase_origin_independent = bool(origin_term and origin_expr is not None and not origin_expr[1])

    # The fractional part itself is not differentiated.  The source argument
    # is checked exactly as N * global_offset, and the proven radius identity
    # supplies d_native(global_offset)=J, hence the physical factor N after
    # division by J.
    phase = _int_row(nodes, chart.get("actual_phase"))
    phase_argument = _handle(phase.get("argument")) if phase else None
    phase_product = _int_row(nodes, phase_argument)
    phase_args = tuple(phase_product.get("arguments", ())) if phase_product else ()
    n_expr = _native_expr(nodes, common_n, coordinate, memo)
    phase_product_ok = bool(
        phase and phase.get("operation") == "analytic_unary" and phase.get("name") == "fractional_part"
        and phase_product and phase_product.get("operation") == "product"
        and len(phase_args) == 2 and _handle(common_n) in {_handle(value) for value in phase_args}
        and offset in {_handle(value) for value in phase_args}
    )
    physical_phase_factor_N = bool(
        phase_product_ok and n_expr is not None and not n_expr[1]
        and log_derivative_equals_J and radius_ratio_bound and phase_origin_independent
    )
    passed = bool(
        native_J_constant and log_derivative_equals_J and other_radius_factors_independent
        and radius_ratio_bound and phase_origin_independent and physical_phase_factor_N
    )
    return {
        "passed": passed,
        "radius": _handle(chart.get("radius")),
        "coordinate": coordinate,
        "physical_measure": jacobian,
        "exp": exp_handle,
        "log_argument": log_argument,
        "log_argument_args": list(args),
        "radius_factor_rows": radius_factor_rows,
        "other_radius_factors_native_independent": other_radius_factors_independent,
        "canonical_J_expression": _poly_serial(jacobian_expr[0]),
        "canonical_log_R_expression": _poly_serial(log_expr[0]),
        "canonical_log_R_native_derivative": _poly_serial(log_expr[1]),
        "log_R_native_derivative_equals_J": log_derivative_equals_J,
        "native_J_is_constant": native_J_constant,
        "global_offset": offset,
        "global_offset_log_ratio": log_ratio,
        "global_offset_exact_log_ratio_minus_origin": exact_offset_two_terms,
        "Ra": ra,
        "Ra_native_independent": radius_ratio_bound,
        "phase_origin": _handle(phase_origin),
        "phase_origin_native_independent": phase_origin_independent,
        "phase_argument": phase_argument,
        "phase_argument_product": phase_product_ok,
        "physical_phase_factor_after_J_division_is_N": physical_phase_factor_N,
        "source_identity": "exact d_native(log R)=J and d_native(N*global_offset)/J=N" if passed else None,
    }


def _sequence(values: Any, minimum: int = 2) -> bool:
    return isinstance(values, (list, tuple)) and len(values) >= minimum and all(_handle(x) is not None for x in values[:minimum])


def _direct_factor(nodes: Any, row_values: Any, factor: Any) -> dict[str, Any]:
    """Check the direct physical conversion factor, without expanding it."""
    factor_i = _handle(factor)
    values = row_values if isinstance(row_values, (list, tuple)) else ()
    rows = []
    for value in values[:2]:
        node = _int_row(nodes, value)
        args = node.get("arguments", ()) if node else ()
        rows.append({
            "row": _handle(value),
            "operation": node.get("operation") if node else None,
            "factor": factor_i,
            "direct_factor_count": sum(_handle(x) == factor_i for x in args) if isinstance(args, (list, tuple)) else 0,
        })
    return {"passed": bool(rows) and all(row["direct_factor_count"] == 1 for row in rows), "rows": rows}


def _slow_b_source_binding(
    nodes: Any,
    slow_values: Any,
    loop_values: Any,
    jacobian: Any,
) -> dict[str, Any]:
    """Bind each slow-B row to the generated native lazy derivative.

    Native differentiation replaces the original B node by a derivative DAG,
    so reachability of the original B handle is only diagnostic and is never
    used as the identity test.
    """
    rows = []
    slow_seq = slow_values if isinstance(slow_values, (list, tuple)) else ()
    loop_seq = loop_values if isinstance(loop_values, (list, tuple)) else ()
    jacobian_i = _handle(jacobian)
    for order in _ORDINARY_Z_KEYS:
        slow_i = _handle(slow_seq[order]) if len(slow_seq) > order else None
        loop_i = _handle(loop_seq[order]) if len(loop_seq) > order else None
        slow_row = _int_row(nodes, slow_i)
        loop_row = _int_row(nodes, loop_i)
        args = tuple(slow_row.get("arguments", ())) if slow_row else ()
        inverse_row = None
        lazy_row = None
        for value in args:
            candidate = _int_row(nodes, value)
            if not candidate:
                continue
            if candidate.get("operation") == "positive_quotient":
                inverse_row = candidate
            if candidate.get("operation") == "original_lazy_flat_branch":
                lazy_row = candidate
        inverse_ok = bool(
            slow_row and slow_row.get("operation") == "product" and len(args) == 2
            and inverse_row
            and inverse_row.get("numerator") == 1
            and _handle(inverse_row.get("denominator")) == jacobian_i
            and inverse_row.get("source_positive_certificate") == "positive original chart log-radius Jacobian"
        )
        lazy_ok = bool(
            lazy_row
            and lazy_row.get("same_original_lazy_native_derivative") is True
            and lazy_row.get("ordinary_slow_Z_derivative_order") == order
            and lazy_row.get("active_body") is not None
            and lazy_row.get("flat_value") is not None
        )
        source_ok = bool(
            loop_row
            and loop_row.get("operation") == "original_lazy_flat_branch"
            and loop_row.get("same_original_flat_extension_through_order5") is True
            and loop_row.get("ordinary_slow_Z_derivative_order") == order
        )
        rows.append({
            "order": order,
            "slow_B_node": slow_i,
            "loop_B_source_node": loop_i,
            "inverse_J_row": inverse_row,
            "generated_lazy_native_derivative_row": lazy_row,
            "inverse_J_binding": inverse_ok,
            "lazy_native_derivative_binding": lazy_ok,
            "loop_B_source_binding": source_ok,
            "original_B_reachability_diagnostic": _contains(nodes, slow_i, loop_i),
            "original_B_reachability_required": False,
        })
    return {"passed": bool(rows) and all(
        row["inverse_J_binding"] and row["lazy_native_derivative_binding"] and row["loop_B_source_binding"]
        for row in rows
    ), "rows": rows}


def _ast_bundle(root: Path) -> dict[str, Any]:
    paths = {name: _resolve(root, relative) for name, relative in _SOURCE_FILES.items()}

    def ev(name: str, selector: str, statements: Iterable[str]) -> dict[str, Any]:
        path = paths[name]
        if path is None:
            return {"scope": selector, "passed": False, "source": _SOURCE_FILES[name], "statements": []}
        result = _function_evidence(path, selector, statements)
        result["source"] = str(path.relative_to(root)).replace("\\", "/")
        return result

    evidence = {
        "native_init": ev("modified", "NativeDerivative.__init__", (
            "self.cache={phase.node:g.zero}", "self.core_extensions=[]")),
        "native_seed": ev("modified", "NativeDerivative.seed", (
            "self.g.ids([left,right])", "raise ValueError('Conflicting native source derivative')",
            "raise ValueError('A constant cannot have a nonzero radial derivative')",
            "self.cache[left.node]=right")),
        "native_substitution": ev("modified", "NativeDerivative.sub", (
            "self.g.node('function_substitution',expression=body.node,variable=variable.node,value=value.node,Z_independent_substitution=True,exact_native_chain_rule_substitution=True)",)),
        "native_derivative": ev("modified", "NativeDerivative.derivative", (
            "g.quotient(g.sub(g.mul(D(u.node),v),g.mul(u,D(v.node))),g.mul(v,v),row['source_positive_certificate'])",
            "g.mul(outer,dx)", "g.add(direct,chain)",
            "g.add(value,g.mul(self.sub(body,var,upper),self.derivative(upper)),g.neg(g.mul(self.sub(body,var,lower),self.derivative(lower))))",
            "g.add(direct,g.mul(at(self.other(angle).derivative(body)),self.derivative(inverse)))",
            "g.quotient(g.sub(self.derivative(target),at(held)),lam,'original positive phase lambda; exact implicit native derivative')",
            "attrs['active_body']=D(row['active_body']).node", "attrs['flat_value']=D(row['flat_value']).node",
            "attrs['same_original_lazy_native_derivative']=True")),
        "native_constructor": ev("modified", "CurrentOriginalModifiedQ4PressureC3Functions.__init__", (
            "NativeDerivative(g,coordinate,chart['actual_phase'])", "derivative.seed(ref,g.mul(jac,dr))",
            "a.scale(Qy,chart['radius'])")),
        "radial_formula": ev("modified", "total_radial_functions", (
            "g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
            "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref)) for ref in rows(chart['loop'][key]))) for key in ('A','B_over_Pstar')}",
            "g.quotient(g.one,N,'same original exact positive common N')",
            "a.add(a.scale(slow['A'],invN),chart['loop']['A_phi'])",
            "a.mul(ex,a.add(a.trim(bg['radial']['E_y']),a.mul(chart['roots']['E'],xy)))",
            "a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
            "a.sub(chart['modified_fields']['V'],chart['modified_fields']['m'])")),
        "pressure_formula": ev("modified", "physical_pressure_fields", (
            "a.scale(fields['V'],Pstar)", "a.scale(radial['V_y'],Pstar)",
            "a.scale(fields['m'],Pstar)", "a.sub(U,m)",
            "b.div(q,L,'same original L=1-delta*Z^2')",
            "b.sub(b.scale(b.mul(z,b.trim(u)),2),b.add(b.mul(b.sub(one,de),b.mul(z,b.trim(mean))),b.mul(d,b.shift(mean))))",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))")),
        "roots": ev("loop", "generic_roots", (
            "a.sub(E,a.scale(a.trim(radial['E_y']),2))", "a.scale(a.trim(radial['V_y']),2)",
            "a.div(jet,E,'same strictly positive actual prescribed swirl E')",
            "aa,bb=divide(C),divide(B)", "p2=generic['p2']",
            "a.neg(a.div(B,C,'same original a>0 and C=a*E'))")),
        "loop_phase": ev("loop", "loop_functions", (
            "gamma=a.sub(a.fixed(g.mul(c(2),eta)),roots['Delta'])",
            "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
            "u=a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,'same positive canonical d_star'))",
            "a.div(one,a.sqrt(a.add(one,a.mul(u,u))),'positive sqrt(1+u^2)')",
            "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
            "m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))",
            "a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
            "a.add(a.fixed(angle),T2)", "a.mul(K,P)",
            "lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
            "lazy(a.scale(a.mul(roots['E'],a.sub(a.neg(a.scale(a.mul(roots['a'],a.mul(that,invlambda)),g.quotient(g.one,twopi,'positive 2*pi'))),roots['b'])),'1/2'))")),
        "finite_n": ev("loop", "finite_N_functions", (
            "g.quotient(g.one,N,'same exact positive integer N on every source chart')",
            "a.scale(loop['A'],invN)", "a.scale(loop['B_over_Pstar'],invN)",
            "a.mul(roots['E'],exm1)", "a.add(V,s)")),
        "correction_rates": ev("loop", "<module>", (
            "RATES=dict(m=Fraction(1),h=Fraction(3,2),k=Fraction(3,2),e=Fraction(1),p=Fraction(0))",)),
        "finite_N_mean_density": ev("loop", "finite_N_functions", (
            "inc=dict(m=s,h=r,k=a.add(a.mul(V,r),a.mul(E,s),a.mul(r,s)),e=a.sub(a.add(a.scale(a.mul(V,s),2),a.mul(s,s)),a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2'))),p=a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2')))",)),
        "finite_N_transport": ev("loop", "transport", (
            "width=g.mul(jacobian,g.sub(coordinate,lower))",
            "g.unary('exp',g.mul(g.constant(-rate),width))",
            "g.unary('exp',g.mul(g.constant(-rate),jacobian,g.sub(coordinate,t)))",
            "a.scale(a.at(densities[name],coordinate,t),g.mul(jacobian,kernel))",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)")),
        "bridge_histories": ev("bridge", "histories", (
            "g.unary('exp',g.mul(g.constant(-rate),jacobian,phase))",
            "g.unary('exp',g.mul(g.constant(-rate),jacobian,g.sub(phase,t)))",
            "a.scale(a.at(bodies[name],phase,t),g.mul(jacobian,kernel))",
            "a.scale(initial[name],decay)",
            "a.integral(body,t,g.zero,phase,actual_six_history_rate=rate,actual_history_name=name,original_incoming_memory_retained=True,source_measure_applied_once=True,prescribed_bridge_chart=label)")),
        "bridge_fields": ev("bridge", "generic_fields", (
            "g.quotient(g.one,Pstar,'same strictly positive Pstar')",
            "a.scale(F0,g.mul(g.unary('sqrt',g.mul(c(2),radius)),inverse))",
            "a.mul(W,actual['phi'])", "a.scale(actual['U'],inverse)",
            "a.scale(J['M'],inverse)", "a.scale(actual['U_y'],inverse)",
            "a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y']))",
            "a.scale(a.mul(F0,F0),g.mul(radius,inverse,inverse))")),
        "outer_forms": ev("outer", "CurrentOriginalPrescribedR100RmC6Functions.__init__", (
            "a.scale(phi,'-2/5')", "a.add(a.scale(z,4),a.scale(E,g.sub(g.one,sigma(g,t))))",
            "a.scale(E,g.neg(sigma_prime(g,t)))", "self.add_chart('reference',q,R,length,actual,incoming['moments'])",
            "self.add_chart('restoration',t,R,g.one,actual,self.endpoints['Rz']['moments'])")),
        "outer_chart": ev("outer", "CurrentOriginalPrescribedR100RmC6Functions.add_chart", (
            "previous.histories(a,actual['phi'],actual['U'],incoming,phase,jacobian,'true_'+name)",
            "previous.generic_fields(a,actual,J,self.F0,self.P0,self.Pstar,R,self.delta)")),
        "outer_end": ev("outer", "CurrentOriginalPrescribedR100RmC6Functions.end", (
            "self.g.node('function_substitution',expression=chart['radius'].node,variable=chart['phase'].node,value=value.node,Z_independent_substitution=True)",
            "at(chart['actual'][key])", "at(jet)")),
    }
    return evidence


def build_native_derivative_guards(
    nodes: Any,
    l: Mapping[str, Any],
    modified: Mapping[str, Any],
    outer: Mapping[str, Any],
    source_guards: Mapping[str, Any] | None,
    root: Any = None,
) -> dict[str, Any]:
    """Check native-derivative source premises and selected endpoint bindings."""
    guards = _mapping(source_guards)
    loop = _unwrap(l, "prescribed_loop_C5_bindings")
    mb = _unwrap(modified, "modified_Q4_pressure_C3_bindings")
    ob = _unwrap(outer, "prescribed_R100_Rm_C6_bindings")
    repo = _repo_root(root)

    inherited_pins = _mapping(guards.get("source_pins"))
    pins = {
        name: _strict_pin(repo, relative, inherited_pins.get(name))
        for name, relative in _SOURCE_FILES.items()
    }
    source_pin_passed = all(row["passed"] for row in pins.values())
    source_guard_rows = {name: _guard_row(guards, name) for name in _REQUIRED_GUARDS}
    source_guard_passed = guards.get("passed") is True and all(source_guard_rows.values())

    ast_evidence = _ast_bundle(repo)
    ast_passed = all(_mapping(row).get("passed") is True for row in ast_evidence.values())

    loop_charts = _mapping(loop.get("charts"))
    mod_charts = _mapping(mb.get("charts"))
    outer_charts = _mapping(ob.get("charts"))
    chart_rows: dict[str, Any] = {}
    geometry_rows: dict[str, Any] = {}
    ordinary_rows: dict[str, Any] = {}
    trace_rows: dict[str, Any] = {}
    chart_bindings_passed = True
    ordinary_passed = True
    trace_passed = True
    common_n = loop.get("N")
    n_handle = _handle(common_n)
    n_row = _int_row(nodes, n_handle)
    common_n_passed = bool(
        n_handle is not None
        and n_row
        and n_row.get("operation") == "exact_positive_integer_expression"
        and n_row.get("exact_positive_integer") is True
        and n_row.get("Z_independent") is True
        and n_row.get("one_common_N_for_all_charts_and_control_denominators") is True
    )

    for name in _CHARTS:
        lc = _mapping(loop_charts.get(name))
        mc = _mapping(mod_charts.get(name))
        oc = _mapping(outer_charts.get(name))
        phase = _phase_binding(nodes, lc, common_n)
        positive_j = _positive_jacobian_premise(nodes, lc, ob)
        radius_j = _radius_log_jacobian_binding(nodes, lc, common_n, loop.get("phase_origin"))
        required_chart_keys = all(key in lc for key in (
            "coordinate", "domain", "radius", "physical_measure", "global_offset", "actual_phase",
            "roots", "loop", "finite_N", "modified_fields"))
        required_modified = all(key in mc for key in ("radial", "physical", "finite_boundary", "radial_background"))
        required_outer = all(key in oc for key in ("fields", "radial", "actual", "moments"))
        geometry_aliases = {
            "coordinate_is_outer_phase": lc.get("coordinate") == oc.get("phase"),
            "radius_is_outer_radius": lc.get("radius") == oc.get("radius"),
            "physical_measure_is_outer_J": lc.get("physical_measure") == oc.get("physical_measure"),
        }
        chart_ok = bool(required_chart_keys and required_modified and required_outer and phase["passed"]
                        and positive_j["passed"] and radius_j["passed"] and all(geometry_aliases.values()))
        chart_rows[name] = {
            "passed": chart_ok,
            "phase": phase,
            "positive_physical_measure": positive_j,
            "radius_log_jacobian": radius_j,
            "geometry_aliases": geometry_aliases,
            "loop_keys_present": required_chart_keys,
            "modified_keys_present": required_modified,
            "outer_keys_present": required_outer,
        }
        geometry_rows[name] = {
            "coordinate": _handle(lc.get("coordinate")),
            "radius": _handle(lc.get("radius")),
            "physical_measure": _handle(lc.get("physical_measure")),
            "global_offset": _handle(lc.get("global_offset")),
            "actual_phase": _handle(lc.get("actual_phase")),
            "domain": lc.get("domain"),
            "phase_held_by_native_derivative": ast_evidence["native_constructor"].get("passed") is True,
            "physical_phase_factor_after_J_division": phase["passed"] and radius_j["passed"] and common_n_passed,
        }
        expected_sequences = {
            "loop.A": _mapping(lc.get("loop")).get("A"),
            "loop.B_over_Pstar": _mapping(lc.get("loop")).get("B_over_Pstar"),
            "loop.A_phi": _mapping(lc.get("loop")).get("A_phi"),
            "loop.B_over_Pstar_phi": _mapping(lc.get("loop")).get("B_over_Pstar_phi"),
            "finite_N.E_N": _mapping(lc.get("finite_N")).get("E_N"),
            "finite_N.V_N": _mapping(lc.get("finite_N")).get("V_N"),
            "modified.radial.V_y": _mapping(mc.get("radial")).get("V_y"),
            "modified.radial.slow_B_y": _mapping(mc.get("radial")).get("slow_B_y"),
            "modified.physical.U": _mapping(mc.get("physical")).get("U"),
            "modified.physical.U_y": _mapping(mc.get("physical")).get("U_y"),
            "modified.physical.m": _mapping(mc.get("physical")).get("m"),
            "modified.physical.m_y": _mapping(mc.get("physical")).get("m_y"),
            "modified.physical.Q": _mapping(mc.get("physical")).get("Q"),
            "modified.physical.Q_y": _mapping(mc.get("physical")).get("Q_y"),
            "modified.boundary.R_Q_y": _mapping(mc.get("finite_boundary")).get("R_Q_y"),
            "outer.fields.E": _mapping(oc.get("fields")).get("E"),
            "outer.fields.V": _mapping(oc.get("fields")).get("V"),
            "outer.fields.m": _mapping(oc.get("fields")).get("m"),
            "outer.radial.E_y": _mapping(oc.get("radial")).get("E_y"),
            "outer.radial.V_y": _mapping(oc.get("radial")).get("V_y"),
            "outer.actual.U": _mapping(oc.get("actual")).get("U"),
            "outer.actual.U_y": _mapping(oc.get("actual")).get("U_y"),
        }
        seq_rows = {key: {"passed": _sequence(value), "handles": list(value[:2]) if isinstance(value, (list, tuple)) else []}
                    for key, value in expected_sequences.items()}
        seq_ok = all(row["passed"] for row in seq_rows.values())
        ordinary_rows[name] = {"passed": seq_ok, "rows": seq_rows, "ordinary_Z_orders": [0, 1]}
        ordinary_passed = ordinary_passed and seq_ok

        loop_b = _mapping(lc.get("loop")).get("B_over_Pstar")
        loop_bphi = _mapping(lc.get("loop")).get("B_over_Pstar_phi")
        slow_b = _mapping(mc.get("radial")).get("slow_B_y")
        vy = _mapping(mc.get("radial")).get("V_y")
        uy = _mapping(mc.get("physical")).get("U_y")
        rqy = _mapping(mc.get("finite_boundary")).get("R_Q_y")
        jac = lc.get("physical_measure")
        pstar = ob.get("Pstar")
        radius = lc.get("radius")
        slow_b_source = _slow_b_source_binding(nodes, slow_b, loop_b, jac)
        trace_rows[name] = {
            "slow_B_rows": [
                {"order": order, "slow_B": _handle(slow_b[order]) if isinstance(slow_b, (list, tuple)) and len(slow_b) > order else None,
                 "contains_native_J": _contains(nodes, slow_b[order], jac) if isinstance(slow_b, (list, tuple)) and len(slow_b) > order else False,
                 "contains_loop_B": _contains(nodes, slow_b[order], loop_b[order]) if isinstance(slow_b, (list, tuple)) and isinstance(loop_b, (list, tuple)) and len(slow_b) > order and len(loop_b) > order else False}
                for order in _ORDINARY_Z_KEYS
            ],
            "V_y_rows": [
                {"order": order, "contains_slow_B": _contains(nodes, vy[order], slow_b[order]) if isinstance(vy, (list, tuple)) and isinstance(slow_b, (list, tuple)) and len(vy) > order and len(slow_b) > order else False,
                 "contains_B_phi": _contains(nodes, vy[order], loop_bphi[order]) if isinstance(vy, (list, tuple)) and isinstance(loop_bphi, (list, tuple)) and len(vy) > order and len(loop_bphi) > order else False}
                for order in _ORDINARY_Z_KEYS
            ],
            "slow_B_source_binding": slow_b_source,
            "physical_U_y_factor": _direct_factor(nodes, uy, pstar),
            "physical_R_Q_y_factor": _direct_factor(nodes, rqy, radius),
            "one_native_J_division_source": ast_evidence["radial_formula"].get("passed") is True,
            "one_fast_phase_addition_source": ast_evidence["radial_formula"].get("passed") is True,
        }
        local_trace = slow_b_source["passed"]
        local_trace = local_trace and all(row["contains_slow_B"] and row["contains_B_phi"] for row in trace_rows[name]["V_y_rows"])
        local_trace = local_trace and trace_rows[name]["physical_U_y_factor"]["passed"] and trace_rows[name]["physical_R_Q_y_factor"]["passed"]
        trace_rows[name]["passed"] = local_trace
        trace_passed = trace_passed and local_trace

        if not chart_ok:
            chart_bindings_passed = False

    finite_n_rows = {}
    finite_n_passed = common_n_passed
    for name, chart in loop_charts.items():
        finite = _mapping(chart).get("finite_N")
        rows = _mapping(finite).get("A_over_N")
        ok = _sequence(rows) and all(_contains(nodes, row, common_n) for row in rows[:2])
        finite_n_rows[name] = {"passed": ok, "A_over_N": list(rows[:2]) if isinstance(rows, (list, tuple)) else [], "common_N": n_handle}
        finite_n_passed = finite_n_passed and ok

    seed_rows = {}
    seed_passed = True
    for name in _CHARTS:
        mc = _mapping(mod_charts.get(name)); oc = _mapping(outer_charts.get(name))
        mrb = _mapping(mc.get("radial_background")); orad = _mapping(oc.get("radial"))
        e_ok = mrb.get("E") == orad.get("E_y")
        v_ok = mrb.get("V") == orad.get("V_y")
        row = {"passed": e_ok and v_ok, "E_seed_rows": mrb.get("E"), "outer_E_y_rows": orad.get("E_y"),
               "V_seed_rows": mrb.get("V"), "outer_V_y_rows": orad.get("V_y")}
        seed_rows[name] = row
        seed_passed = seed_passed and row["passed"]

    checks = {
        "source_pins": {"passed": source_pin_passed, "rows": pins},
        "inherited_source_guards": {"passed": source_guard_passed, "rows": source_guard_rows,
                                     "source_guards_passed": guards.get("passed") is True},
        "source_ast": {"passed": ast_passed, "rows": ast_evidence},
        "common_positive_finite_N": {"passed": finite_n_passed, "N": n_handle, "N_node": n_row, "charts": finite_n_rows},
        "held_actual_phase_and_geometry": {"passed": chart_bindings_passed, "charts": chart_rows, "handles": geometry_rows},
        "background_radial_seed_binding": {"passed": seed_passed, "charts": seed_rows},
        "ordinary_Z1_rows": {"passed": ordinary_passed, "charts": ordinary_rows},
        "slow_B_to_total_Vy_to_physical_Uy": {"passed": trace_passed, "charts": trace_rows},
    }
    passed = all(row.get("passed") is True for row in checks.values())
    return {
        "passed": passed,
        "method": "pinned_native_derivative_source_semantics_and_selected_bindings",
        "source_only": True,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "checks": checks,
        "source_pins": pins,
        "selected_charts": list(_CHARTS),
        "premises": {
            "native_phase_is_held_then_Aphi_Bphi_added": True,
            "native_background_seed_is_J_times_source_radial_row": True,
            "slow_B_is_native_partial_divided_by_J_then_divided_by_common_N": True,
            "total_Vy_keeps_background_slow_term_and_Bphi": True,
            "physical_Uy_has_one_Pstar_factor": True,
            "physical_RQy_has_one_radius_factor": True,
            "ordinary_Z1_rows_are_retained": True,
        },
        "unproved": {
            "arbitrary_fast_phase_reference_restoration_naturality": True,
            "common_physical_input_two_jet_identity": True,
            "periodic_lift_consumed_by_native_graph": True,
            "endpoint_derivative_graph_equality": True,
            "cap_substitution_or_cap_cancellation": True,
            "numeric_endpoint_value_or_ancestor_constructor": True,
        },
        "limits": {
            "naturality_theorem_proved": False,
            "B_y_alias_authorized": False,
            "Qy_alias_authorized": False,
            "RQy_alias_authorized": False,
            "positive_budget_reinterpretation": False,
        },
    }


__all__ = ["build_native_derivative_guards"]
