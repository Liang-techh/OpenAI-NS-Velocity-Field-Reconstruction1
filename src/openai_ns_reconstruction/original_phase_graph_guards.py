"""Structural guards for the archived fractional-phase loop graph.

The checks in this module inspect graph handles only.  They do not construct a
producer, evaluate a source value, or treat an archive flag or error cap as a
function identity.  In particular, this module validates the graph wiring
needed by a later periodic-lift theorem; it does not prove that theorem.
"""

from __future__ import annotations

from collections import Counter
from typing import Any, Mapping


_REQUIRED_CHARTS = ("reference", "restoration")
_OUTPUTS = ("A", "B_over_Pstar", "A_phi", "B_over_Pstar_phi")


def _map(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _index(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    node = getattr(value, "node", None)
    return node if isinstance(node, int) else None


def _row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _index(value)
    if index is None:
        return None
    try:
        result = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return result if isinstance(result, Mapping) else None


def _children(row: Mapping[str, Any] | None) -> tuple[int, ...]:
    """Return graph edges, excluding lexical variable declarations."""
    if not isinstance(row, Mapping):
        return ()
    operation = row.get("operation")
    keys: tuple[str, ...]
    if operation in {"sum", "product", "logical_or", "logical_and"}:
        keys = ("arguments",)
    elif operation in {
        "negative", "analytic_unary", "exact_original_source_flat_sigma",
        "fractional_part", "regular_even_radial_axis_continuation",
    }:
        keys = ("argument",)
    elif operation == "positive_quotient":
        keys = ("numerator", "denominator")
    elif operation == "original_lazy_flat_branch":
        keys = ("active_body", "Delta", "eta", "flat_value", "flat_predicate",
                "original_collar_coordinate", "original_collar_upper")
    elif operation == "original_monotone_phase_inverse":
        keys = ("phase_function", "target_phase", "angle_lower", "angle_upper", "Delta", "eta")
    elif operation in {"substitute_original_inverse_angle", "substitute_original_native_coordinate",
                       "substitute_free_query_at_native_coordinate"}:
        keys = ("body", "inverse_angle", "coordinate", "value", "lower", "upper")
    elif operation == "function_substitution":
        # ``variable`` is a lexical binder, not a graph dependency of the
        # substituted value.  Its scope is handled explicitly by _z_free.
        keys = ("expression", "value")
    elif operation == "definite_integral":
        keys = ("integrand", "lower", "upper")
    elif operation == "exact_real_comparison":
        keys = ("left", "right")
    else:
        return ()
    result: list[int] = []
    for key in keys:
        value = row.get(key)
        if isinstance(value, (list, tuple)):
            result.extend(item for item in value if _index(item) is not None)
        elif _index(value) is not None:
            result.append(_index(value))
    return tuple(result)


def _graph_indices(nodes: Any):
    if isinstance(nodes, Mapping):
        return tuple(int(key) for key in nodes if isinstance(key, int) or str(key).isdigit())
    return tuple(range(len(nodes))) if isinstance(nodes, (list, tuple)) else ()


def _is_exact(nodes: Any, value: Any, numerator: int, denominator: int = 1) -> bool:
    row = _row(nodes, value)
    return bool(row and row.get("operation") == "exact_rational" and
                row.get("numerator") == numerator and row.get("denominator") == denominator)


def _product(nodes: Any, factors: tuple[int, ...]) -> list[int]:
    wanted = tuple(sorted(factors))
    return [i for i in _graph_indices(nodes)
            if (_row(nodes, i) or {}).get("operation") == "product" and
            tuple(sorted(_children(_row(nodes, i)))) == wanted]


def _flatten_operation(nodes: Any, value: Any, operation: str) -> tuple[int, ...]:
    """Flatten only associative graph nodes of one operation."""
    index = _index(value)
    row = _row(nodes, index)
    if index is None or row is None or row.get("operation") != operation:
        return (index,) if index is not None else ()
    result: list[int] = []
    for child in _children(row):
        result.extend(_flatten_operation(nodes, child, operation))
    return tuple(result)


def _square_of(nodes: Any, value: Any, base: Any) -> bool:
    """Recognize a literal product base*base, allowing associative nesting."""
    index = _index(value)
    base_index = _index(base)
    if index is None or base_index is None:
        return False
    factors = _flatten_operation(nodes, index, "product")
    return Counter(factors) == Counter((base_index, base_index))


def _two_pi(nodes: Any, value: Any) -> bool:
    factors = _flatten_operation(nodes, value, "product")
    return (len(factors) == 2 and
            any(_is_exact(nodes, child, 2) for child in factors) and
            any((_row(nodes, child) or {}).get("operation") == "mathematical_pi"
                for child in factors))


def _reachable(nodes: Any, root: Any) -> set[int]:
    start = _index(root)
    if start is None:
        return set()
    visited: set[int] = set()
    stack = [start]
    while stack:
        index = stack.pop()
        if index in visited:
            continue
        visited.add(index)
        stack.extend(child for child in _children(_row(nodes, index)) if child not in visited)
    return visited


def _z_free(nodes: Any, root: Any) -> dict[str, Any]:
    """Check Z-freeness with lexical function-substitution scopes."""
    visited: set[int] = set()
    unknown: list[dict[str, Any]] = []
    memo: dict[tuple[int, tuple[int, ...]], bool] = {}

    def visit(value: Any, bound: frozenset[int]) -> bool:
        index = _index(value)
        if index is None:
            unknown.append({"node": value, "reason": "noninteger graph handle"})
            return False
        visited.add(index)
        if index in bound:
            return True
        key = (index, tuple(sorted(bound)))
        if key in memo:
            return memo[key]
        row = _row(nodes, index)
        if row is None:
            unknown.append({"node": index, "reason": "missing graph row"})
            memo[key] = False
            return False
        operation = row.get("operation")
        if operation == "bound_variable":
            answer = row.get("name") != "Z"
        elif operation in {"exact_rational", "mathematical_pi"}:
            answer = True
        elif operation == "exact_positive_integer_expression":
            answer = row.get("Z_independent") is True
        elif operation == "current_original_source_parameter":
            # This is a source parameter leaf.  The source graph's explicit
            # Z-independent property is accepted as leaf provenance, while
            # caps, report flags, and numeric values are never consulted.
            answer = row.get("Z_independent") is True
        elif operation == "function_substitution":
            variable = _index(row.get("variable"))
            scope = bound | ({variable} if variable is not None else set())
            answer = visit(row.get("expression"), frozenset(scope)) and visit(row.get("value"), bound)
        elif operation == "definite_integral":
            answer = all(visit(row.get(key), bound) for key in ("integrand", "lower", "upper"))
        else:
            children = _children(row)
            if not children:
                unknown.append({"node": index, "operation": operation,
                                "reason": "opaque or unsupported source leaf"})
                answer = False
            else:
                answer = all(visit(child, bound) for child in children)
        memo[key] = answer
        return answer

    passed = visit(root, frozenset())
    return {
        "passed": bool(passed),
        "root": _index(root),
        "visited_handles": sorted(visited),
        "unknown": unknown[:12],
        "unknown_count": len(unknown),
    }


def _phase_product(nodes: Any, chart: Mapping[str, Any], common_n: Any) -> dict[str, Any]:
    actual = _index(chart.get("actual_phase"))
    offset = _index(chart.get("global_offset"))
    n_handle = _index(common_n)
    phase = _row(nodes, actual)
    argument = _index(phase.get("argument")) if phase else None
    product = _row(nodes, argument)
    args = tuple(_children(product))
    passed = bool(
        actual is not None and phase and phase.get("operation") == "analytic_unary" and
        phase.get("name") == "fractional_part" and argument is not None and
        product and product.get("operation") == "product" and len(args) == 2 and
        n_handle in args and offset in args and n_handle != offset
    )
    return {
        "passed": passed,
        "actual_phase": actual,
        "phase_node": phase,
        "phase_argument": argument,
        "phase_product_nodes": list(args),
        "global_offset": offset,
        "common_N": n_handle,
    }


def _check_q_chain(nodes: Any, chart: Mapping[str, Any], loop: Mapping[str, Any], bindings: Mapping[str, Any]) -> dict[str, Any]:
    roots = _map(chart.get("roots"))
    q_values = loop.get("q")
    q = [_index(value) for value in q_values] if isinstance(q_values, (list, tuple)) else []
    q = [value for value in q if value is not None]
    delta = _index(roots.get("Delta", [None])[0]) if isinstance(roots.get("Delta"), (list, tuple)) else None
    a_root = _index(roots.get("a", [None])[0]) if isinstance(roots.get("a"), (list, tuple)) else None
    eta = _index(bindings.get("eta"))
    q_rows = []
    q_ok = len(q) >= 2 and delta is not None and eta is not None
    for order, q_handle in enumerate(q):
        row = _row(nodes, q_handle)
        local = bool(
            row and row.get("operation") == "original_lazy_flat_branch" and
            _index(row.get("Delta")) == delta and _index(row.get("eta")) == eta and
            _is_exact(nodes, row.get("flat_value"), 0) and
            _index(row.get("flat_predicate")) is not None
        )
        predicate = _row(nodes, row.get("flat_predicate")) if row else None
        local = local and bool(
            predicate and predicate.get("operation") == "exact_real_comparison" and
            predicate.get("operator") == ">=" and _index(predicate.get("left")) == delta and
            _index(predicate.get("right")) == eta
        )
        active = _row(nodes, row.get("active_body")) if row else None
        # Higher ordinary-Z rows are stored as Leibniz sums.  Their immediate
        # children are derivative products, so inspect the full active subtree
        # instead of requiring the order-zero product shape at every order.
        active_reachable = _reachable(nodes, row.get("active_body")) if row else set()
        sigma = [child for child in active_reachable
                 if (_row(nodes, child) or {}).get("operation") == "exact_original_source_flat_sigma"]
        rootsqrt = [child for child in active_reachable
                    if (_row(nodes, child) or {}).get("operation") == "analytic_unary" and
                    (_row(nodes, child) or {}).get("name") == "positive_sqrt"]
        local = local and bool(
            active and active.get("operation") in {"product", "sum"} and
            sigma and rootsqrt and delta in active_reachable and eta in active_reachable
        )
        valid_sigma = [candidate for candidate in sigma
                       if delta in _reachable(nodes, candidate) and
                       eta in _reachable(nodes, candidate)]
        local = local and bool(valid_sigma)
        if rootsqrt:
            root_quotients = [_row(nodes, _row(nodes, root).get("argument"))
                              for root in rootsqrt]
            local = local and any(
                quotient and quotient.get("operation") == "positive_quotient" and
                delta in _reachable(nodes, quotient.get("numerator")) and
                a_root is not None and a_root in _reachable(nodes, quotient.get("denominator"))
                for quotient in root_quotients
            )
        q_rows.append({
            "order": order,
            "handle": q_handle,
            "passed": bool(local),
            "active_body": _index(row.get("active_body")) if row else None,
            "flat_predicate": _index(row.get("flat_predicate")) if row else None,
        })
        q_ok = q_ok and local

    q0 = q[0] if q else None
    p2_values = roots.get("p2")
    p2 = _index(p2_values[0]) if isinstance(p2_values, (list, tuple)) and p2_values else None
    inv_dstar = [i for i in _graph_indices(nodes)
                  if (_row(nodes, i) or {}).get("operation") == "positive_quotient" and
                  _is_exact(nodes, (_row(nodes, i) or {}).get("numerator"), 1) and
                  _index((_row(nodes, i) or {}).get("denominator")) == _index(bindings.get("dstar"))]
    inv_dstar = inv_dstar[0] if inv_dstar else None
    p2q = _product(nodes, tuple(sorted((p2, q0))))[0] if p2 is not None and q0 is not None and _product(nodes, tuple(sorted((p2, q0)))) else None
    u_values = _product(nodes, tuple(sorted((p2q, inv_dstar)))) if p2q is not None and inv_dstar is not None else []
    u = u_values[0] if u_values else None

    # hinv = 1/sqrt(1+u^2), followed by r = u*hinv and alpha = 2*q*hinv.
    u2 = _product(nodes, (u, u)) if u is not None else []
    sums = [i for i in _graph_indices(nodes)
            if (_row(nodes, i) or {}).get("operation") == "sum" and
            u2 and u2[0] in _children(_row(nodes, i)) and
            any(_is_exact(nodes, child, 1) for child in _children(_row(nodes, i)))]
    rootsqrt = [i for i in _graph_indices(nodes)
                if (_row(nodes, i) or {}).get("operation") == "analytic_unary" and
                (_row(nodes, i) or {}).get("name") == "positive_sqrt" and
                _index((_row(nodes, i) or {}).get("argument")) in sums]
    hinv_values = [i for i in _graph_indices(nodes)
                   if (_row(nodes, i) or {}).get("operation") == "positive_quotient" and
                   _is_exact(nodes, (_row(nodes, i) or {}).get("numerator"), 1) and
                   _index((_row(nodes, i) or {}).get("denominator")) in rootsqrt]
    hinv = hinv_values[0] if hinv_values else None
    r_values = _product(nodes, tuple(sorted((u, hinv)))) if u is not None and hinv is not None else []
    qhinv = _product(nodes, tuple(sorted((q0, hinv)))) if q0 is not None and hinv is not None else []
    exact_two = [i for i in _graph_indices(nodes) if _is_exact(nodes, i, 2)]
    alpha_values = [i for inner in qhinv for two in exact_two
                    for i in _product(nodes, tuple(sorted((inner, two))))]

    phi_values = loop.get("Phi")
    phi = _index(phi_values[0]) if isinstance(phi_values, (list, tuple)) and phi_values else None
    phi_row = _row(nodes, phi)
    k_candidates = [child for child in _children(phi_row)
                    if (_row(nodes, child) or {}).get("operation") == "positive_quotient" and
                    _is_exact(nodes, (_row(nodes, child) or {}).get("numerator"), 1)]
    k_handle = k_candidates[0] if k_candidates else None
    k_row = _row(nodes, k_handle)
    k_den = _index(k_row.get("denominator")) if k_row else None
    t0_values = roots.get("t0")
    t0 = _index(t0_values[0]) if isinstance(t0_values, (list, tuple)) and t0_values else None
    ratio_candidates = [child for child in _children(_row(nodes, k_den))
                        if (_row(nodes, child) or {}).get("operation") == "sum" and
                        q0 is not None and q0 in _reachable(nodes, child)] if k_den is not None else []
    twopi = [child for child in _children(_row(nodes, k_den))
             if (_row(nodes, child) or {}).get("operation") == "product" and
             any(_is_exact(nodes, factor, 2) for factor in _children(_row(nodes, child))) and
             any((_row(nodes, factor) or {}).get("operation") == "mathematical_pi" for factor in _children(_row(nodes, child)))] if k_den is not None else []
    k_ok = bool(
        phi_row and phi_row.get("operation") == "product" and k_handle is not None and
        k_den is not None and len(ratio_candidates) == 1 and len(twopi) == 1 and
        (t0 is None or t0 in _reachable(nodes, ratio_candidates[0]))
    )
    passed = bool(q_ok and q0 is not None and p2q is not None and u is not None and
                  hinv is not None and r_values and alpha_values and k_ok)
    return {
        "passed": passed,
        "q_rows": q_rows,
        "q_active_flat_structural_path": bool(q_ok),
        "branch_id": (delta, eta, _index(_row(nodes, q0).get("flat_predicate"))
                      if q0 is not None and _row(nodes, q0) else None),
        "handles": {"q0": q0, "p2q": p2q, "inverse_dstar": inv_dstar, "u": u,
                     "hinv": hinv, "r": r_values[0] if r_values else None,
                     "alpha": alpha_values[0] if alpha_values else None,
                     "K": k_handle, "K_denominator": k_den},
        "K_ratio_and_two_pi_structural_path": k_ok,
        "metadata_flags_not_used_as_function_values": True,
    }


def _check_c0_recipe(
    nodes: Any,
    chart: Mapping[str, Any],
    loop: Mapping[str, Any],
    q_chain: Mapping[str, Any],
    phi: int | None,
    inverse: int | None,
    phase_lambda: int | None,
    angle: int | None,
) -> dict[str, Any]:
    """Check the compact C0 Poisson/Phi recipe exposed by the archive.

    This deliberately checks graph identities and handle equality only.  It
    does not interpret source values or the archived positive-certification
    prose.
    """
    handles = _map(q_chain.get("handles"))
    q0 = _index(handles.get("q0"))
    k = _index(handles.get("K"))
    roots = _map(chart.get("roots"))
    t0_values = roots.get("t0")
    t0 = _index(t0_values[0]) if isinstance(t0_values, (list, tuple)) and t0_values else None
    failures: list[str] = []

    k_row = _row(nodes, k)
    denominator = _index(k_row.get("denominator")) if k_row else None
    denominator_row = _row(nodes, denominator)
    denominator_factors = _children(denominator_row)
    two_pi_factors = [factor for factor in denominator_factors if _two_pi(nodes, factor)]
    sum_factors = [factor for factor in denominator_factors
                   if (_row(nodes, factor) or {}).get("operation") == "sum"]
    sum_terms = _flatten_operation(nodes, sum_factors[0], "sum") if len(sum_factors) == 1 else ()
    exact_two = {i for i in _graph_indices(nodes) if _is_exact(nodes, i, 2)}
    q_square_twice = []
    t_square = []
    for term in sum_terms:
        factors = _flatten_operation(nodes, term, "product")
        if (q0 is not None and factors.count(q0) == 2 and
                len(factors) == 3 and any(factor in exact_two for factor in factors)):
            q_square_twice.append(term)
        if t0 is not None and Counter(factors) == Counter((t0, t0)):
            t_square.append(term)
    t0_is_zero = _is_exact(nodes, t0, 0)
    one_terms = [term for term in sum_terms if _is_exact(nodes, term, 1)]
    expected_sum_terms = 2 if t0_is_zero else 3
    recognized_terms = set(one_terms) | set(q_square_twice) | set(t_square)
    k_recipe = bool(
        k_row and k_row.get("operation") == "positive_quotient" and
        _is_exact(nodes, k_row.get("numerator"), 1) and
        denominator_row and denominator_row.get("operation") == "product" and
        len(denominator_factors) == 2 and len(two_pi_factors) == 1 and len(sum_factors) == 1 and
        len(sum_terms) == expected_sum_terms and len(one_terms) == 1 and
        len(q_square_twice) == 1 and
        ((t0_is_zero and not t_square) or (not t0_is_zero and len(t_square) == 1)) and
        len(recognized_terms) == len(sum_terms)
    )
    if not k_recipe:
        failures.append("K denominator is not the exact 2*pi*(1+t0^2+2*q^2) graph recipe")

    phi_row = _row(nodes, phi)
    phi_factors = _children(phi_row)
    phi_body_candidates = [child for child in phi_factors if child != k]
    phi_body = phi_body_candidates[0] if len(phi_body_candidates) == 1 else None
    phi_terms = _flatten_operation(nodes, phi_body, "sum") if phi_body is not None else ()
    integral_candidates = [term for term in phi_terms
                           if (_row(nodes, term) or {}).get("operation") == "definite_integral"]
    angle_present = angle is not None and angle in phi_terms
    integral = integral_candidates[0] if len(integral_candidates) == 1 else None
    integral_row = _row(nodes, integral)
    integral_upper = _index(integral_row.get("upper")) if integral_row else None
    integral_lower = _index(integral_row.get("lower")) if integral_row else None
    phi_recipe = bool(
        phi_row and phi_row.get("operation") == "product" and len(phi_factors) == 2 and k in phi_factors and
        len(phi_terms) == 2 and len(integral_candidates) == 1 and angle_present and integral_upper == angle and
        _is_exact(nodes, integral_lower, 0) and
        integral_row.get("variable") == (_row(nodes, angle) or {}).get("name")
    )
    if not phi_recipe:
        failures.append("Phi is not the exact K*(angle+integral t^2) graph recipe")

    direction_values = loop.get("direction_at_inverse")
    t1_values = loop.get("T1_at_inverse")
    direction = _index(direction_values[0]) if isinstance(direction_values, (list, tuple)) and direction_values else None
    t1 = _index(t1_values[0]) if isinstance(t1_values, (list, tuple)) and t1_values else None
    direction_row = _row(nodes, direction)
    t1_row = _row(nodes, t1)
    direction_body = _index(direction_row.get("body")) if direction_row else None
    t1_body = _index(t1_row.get("body")) if t1_row else None
    direction_recipe = bool(
        direction_row and direction_row.get("operation") == "substitute_original_inverse_angle" and
        _index(direction_row.get("inverse_angle")) == inverse and direction_body is not None and
        t1_row and t1_row.get("operation") == "substitute_original_inverse_angle" and
        _index(t1_row.get("inverse_angle")) == inverse and
        _row(nodes, t1_body) and (_row(nodes, t1_body) or {}).get("operation") == "definite_integral" and
        _index((_row(nodes, t1_body) or {}).get("integrand")) == direction_body and
        integral is not None and
        bool(_product(nodes, tuple(sorted((direction_body, direction_body)))) and
             _index((_row(nodes, integral) or {}).get("integrand")) in
             _product(nodes, tuple(sorted((direction_body, direction_body)))))
    )
    if not direction_recipe:
        failures.append("T1_at_inverse/direction_at_inverse do not share the Phi t integrand")

    lambda_row = _row(nodes, phase_lambda)
    lambda_body = _index(lambda_row.get("body")) if lambda_row else None
    lambda_body_row = _row(nodes, lambda_body)
    lambda_reach = _reachable(nodes, lambda_body)
    phi_integrand = _index((integral_row or {}).get("integrand")) if integral_row else None
    lambda_recipe = bool(
        lambda_row and lambda_row.get("operation") == "substitute_original_inverse_angle" and
        _index(lambda_row.get("inverse_angle")) == inverse and
        lambda_body_row and lambda_body_row.get("operation") == "product" and
        len(_children(lambda_body_row)) == 2 and k in _children(lambda_body_row) and
        k is not None and k in lambda_reach and direction_body is not None and
        direction_body in lambda_reach and phi_integrand is not None and phi_integrand in lambda_reach
    )
    if not lambda_recipe:
        failures.append("phase_lambda_at_inverse is not bound to the selected K/T1/Phi chain")

    output_t_bindings: dict[str, bool] = {}
    for key in _OUTPUTS:
        values = loop.get(key)
        root = _index(values[0]) if isinstance(values, (list, tuple)) and values else None
        reach = _reachable(nodes, root)
        output_t_bindings[key] = bool(
            direction_body is not None and direction_body in reach and
            all(handle is not None and handle in reach
                for handle in handles.values() if isinstance(handle, int))
        )
    if not all(output_t_bindings.values()):
        failures.append("A/B prefactors are not all bound to the selected Poisson t chain")

    return {
        "passed": not failures,
        "K_recipe": k_recipe,
        "Phi_recipe": phi_recipe,
        "direction_T1_recipe": direction_recipe,
        "phase_lambda_recipe": lambda_recipe,
        "output_t_bindings": output_t_bindings,
        "handles": {"q": q0, "t0": t0, "K": k, "direction_body": direction_body,
                     "T1_body": t1_body, "Phi_integral": integral},
        "failures": failures,
    }


def _check_chart(nodes: Any, chart_name: str, chart: Mapping[str, Any], bindings: Mapping[str, Any]) -> dict[str, Any]:
    loop = _map(chart.get("loop"))
    phase = _phase_product(nodes, chart, bindings.get("N"))
    actual = phase.get("actual_phase")
    inverse_values = loop.get("inverse")
    phi_values = loop.get("Phi")
    lambda_values = loop.get("phase_lambda_at_inverse")
    inverse = _index(inverse_values[0]) if isinstance(inverse_values, (list, tuple)) and inverse_values else None
    phi = _index(phi_values[0]) if isinstance(phi_values, (list, tuple)) and phi_values else None
    lam = _index(lambda_values[0]) if isinstance(lambda_values, (list, tuple)) and lambda_values else None
    inverse_row = _row(nodes, inverse)
    angle = _index(loop.get("angle"))
    angle_row = _row(nodes, angle)
    lambda_row = _row(nodes, lam)
    same_inverse = bool(
        inverse_row and inverse_row.get("operation") == "original_monotone_phase_inverse" and
        _index(inverse_row.get("target_phase")) == actual and
        _index(inverse_row.get("phase_function")) == phi and
        inverse_row.get("phase_endpoints") == [0, 1] and
        angle_row and angle_row.get("operation") == "bound_variable" and
        inverse_row.get("angle_variable") == angle_row.get("name") and
        lambda_row and lambda_row.get("operation") == "substitute_original_inverse_angle" and
        _index(lambda_row.get("inverse_angle")) == inverse and
        lambda_row.get("angle_variable") == angle_row.get("name")
    )
    q_chain = _check_q_chain(nodes, chart, loop, bindings)
    chain_handles = [value for value in _map(q_chain.get("handles")).values()
                     if isinstance(value, int)]
    phi_reach = _reachable(nodes, phi)
    inverse_reach = _reachable(nodes, inverse)
    lambda_reach = _reachable(nodes, lam)
    chain_in_phi = bool(chain_handles) and all(handle in phi_reach for handle in chain_handles)
    chain_in_inverse_substitution = bool(chain_handles) and all(
        handle in inverse_reach and handle in lambda_reach for handle in chain_handles
    )
    chain_binding = {
        "passed": chain_in_phi and chain_in_inverse_substitution,
        "chain_handles": chain_handles,
        "selected_Phi_reachable": chain_in_phi,
        "actual_inverse_reachable": bool(chain_handles) and all(handle in inverse_reach for handle in chain_handles),
        "phase_substitution_reachable": bool(chain_handles) and all(handle in lambda_reach for handle in chain_handles),
    }
    outputs = {}
    branch_ids = set()
    output_graph_ok = True
    for key in _OUTPUTS:
        values = loop.get(key)
        rows = []
        values_ok = isinstance(values, (list, tuple)) and len(values) >= 2
        for order, value in enumerate(values if isinstance(values, (list, tuple)) else ()):
            handle = _index(value)
            row = _row(nodes, handle)
            local = bool(
                row and row.get("operation") == "original_lazy_flat_branch" and
                _is_exact(nodes, row.get("flat_value"), 0) and
                _index(row.get("flat_predicate")) is not None
            )
            if row:
                branch_ids.add((_index(row.get("Delta")), _index(row.get("eta")), _index(row.get("flat_predicate"))))
            if order < 2:
                reachable = _reachable(nodes, handle)
                # A and B are native phase outputs and do not contain the
                # post-inverse substitution node.  The phi-suffixed rows do.
                required = {x for x in (actual, inverse, phi) if x is not None}
                if key.endswith("_phi") and lam is not None:
                    required.add(lam)
                local = local and required <= reachable
            rows.append({"order": order, "handle": handle, "passed": local})
            values_ok = values_ok and local
        outputs[key] = {"passed": values_ok, "rows": rows[:6]}
        output_graph_ok = output_graph_ok and values_ok
    branch_id = q_chain.get("branch_id")
    branch_common = bool(branch_id is not None and branch_ids == {tuple(branch_id)})

    z_handles = {
        "global_offset": _index(chart.get("global_offset")),
        "phase_argument": _index(phase.get("phase_argument")),
        "radius": _index(chart.get("radius")),
        "physical_measure": _index(chart.get("physical_measure")),
        "N": _index(bindings.get("N")),
        "phase_origin": _index(bindings.get("phase_origin")),
    }
    missing_z_handles = [name for name, value in z_handles.items() if value is None]
    z_checks = {name: _z_free(nodes, value) for name, value in z_handles.items() if value is not None}
    z_ok = not missing_z_handles and all(item.get("passed") is True for item in z_checks.values())
    c0_recipe = _check_c0_recipe(nodes, chart, loop, q_chain, phi, inverse, lam, angle)
    passed = bool(
        phase.get("passed") and same_inverse and output_graph_ok and branch_common and
        q_chain["passed"] and chain_binding["passed"] and z_ok and c0_recipe["passed"]
    )
    return {
        "passed": passed,
        "chart": chart_name,
        "phase": phase,
        "inverse_binding": {"passed": same_inverse, "inverse": inverse, "Phi": phi,
                             "phase_lambda_at_inverse": lam, "angle": angle,
                             "angle_variable": angle_row.get("name") if angle_row else None},
        "outputs": outputs,
        "common_lazy_branch_scope": branch_common,
        "q_chain": q_chain,
        "chain_binding": chain_binding,
        "C0_recipe": c0_recipe,
        "missing_required_Z_handles": missing_z_handles,
        "Z_free": z_checks,
        "visited_actual_handles": {name: item["visited_handles"] for name, item in z_checks.items()},
    }


def _diagnose_other_chart(nodes: Any, chart: Mapping[str, Any], bindings: Mapping[str, Any]) -> dict[str, Any]:
    phase = _phase_product(nodes, chart, bindings.get("N"))
    loop = _map(chart.get("loop"))
    inv_values = loop.get("inverse")
    inverse = _index(inv_values[0]) if isinstance(inv_values, (list, tuple)) and inv_values else None
    row = _row(nodes, inverse)
    return {"phase_product": phase.get("passed") is True,
            "inverse_target_matches_phase": bool(row and _index(row.get("target_phase")) == phase.get("actual_phase")),
            "actual_phase": phase.get("actual_phase"), "inverse": inverse}


def build_phase_graph_guards(nodes: Any, l: Mapping[str, Any]) -> dict[str, Any]:
    """Validate reference/restoration phase wiring in an accepted graph."""
    bindings = _map(l.get("prescribed_loop_C5_bindings")) if isinstance(l, Mapping) and "prescribed_loop_C5_bindings" in l else _map(l)
    charts = _map(bindings.get("charts"))
    missing = [name for name in _REQUIRED_CHARTS if not isinstance(charts.get(name), Mapping)]
    if missing:
        return {"passed": False, "method": "archived_graph_structural_phase_binding",
                "charts": {}, "unproved": ["missing charts: " + ", ".join(missing)]}
    rows = {name: _check_chart(nodes, name, charts[name], bindings) for name in _REQUIRED_CHARTS}
    diagnostics = {
        name: _diagnose_other_chart(nodes, chart, bindings)
        for name, chart in charts.items() if name not in _REQUIRED_CHARTS and isinstance(chart, Mapping)
    }
    unproved = [
        "global periodic identity Phi(theta+2*pi)=Phi(theta)+1",
        "endpoint first parameter jets and mixed ordinary-Z/phase jets",
        "analytic positivity/regularity of q=0 branch beyond archived graph wiring",
    ]
    passed = all(row.get("passed") is True for row in rows.values())
    return {
        "passed": passed,
        "method": "archived_graph_structural_phase_binding",
        "source_only": True,
        "numeric_caps_used": False,
        "charts": rows,
        "other_chart_diagnostics": diagnostics,
        "common_structural_chain": passed,
        "unproved": unproved,
        "claims_not_made": ["periodic_C1_lift", "B_y_alias", "RQ_y_alias", "numeric_source_value"],
    }


__all__ = ["build_phase_graph_guards"]
