"""Source-only native-calculus guards for the actual finite-N Rm patch.

The Rm chart is a moving native coordinate: ``R = Rm*x`` and ``J = 1/x``.
This guard intentionally keeps that distinction visible.  It checks the
serialized graph and the pinned producer source without importing a producer,
evaluating a cap, or replaying an accepted ancestor constructor.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

from .original_core_numeric import ROOT
from .original_native_derivative_guards import (
    _contains,
    _direct_factor,
    _function_evidence,
    _handle,
    _int_row,
    _mapping,
    _phase_binding,
)


_LEADING_KEY = "true_Rm_patch_C6_bindings"
_RM_KEY = "Rm_patch_finite_N_C5_bindings"
_LOOP_KEY = "prescribed_loop_C5_bindings"
_ORDERS = tuple(range(6))

_SOURCE_PINS = {
    "leading_patch": (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py",
        "a6f1cb5deaf17e9482aa964aef7867aadbc815e4d703a415e4945cf665f14002",
    ),
    "finite_patch": (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py",
        "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    ),
    "native_derivative": (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py",
        "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    ),
}


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    if isinstance(value, Mapping) and isinstance(value.get(key), Mapping):
        return value[key]
    return value if isinstance(value, Mapping) else {}


def _row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _handle(value)
    if index is None or index < 0:
        return None
    return _int_row(nodes, index)


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return ROOT


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


def _sequence(value: Any, minimum: int = 2) -> bool:
    return isinstance(value, (list, tuple)) and len(value) >= minimum and all(
        _handle(item) is not None for item in value[:minimum]
    )


def _native_independent(nodes: Any, value: Any, coordinate: Any, memo: dict[int, bool] | None = None) -> bool:
    """Conservatively prove that a selected graph expression has no x leaf.

    The general reachability helper intentionally follows every serialized
    integer and therefore reaches unrelated shared suffixes in this archive.
    This smaller evaluator follows operation operands only, which is the
    appropriate test for the substituted Rm factor.
    """
    index = _handle(value)
    coordinate_i = _handle(coordinate)
    if index is None or coordinate_i is None:
        return False
    memo = memo if memo is not None else {}
    if index in memo:
        return memo[index]
    row = _row(nodes, index)
    if row is None:
        memo[index] = False
        return False
    operation = row.get("operation")
    if operation == "bound_variable":
        result = index != coordinate_i
    elif operation in {
        "exact_rational", "mathematical_pi", "current_original_source_parameter",
        "original_source_parameter", "exact_positive_integer_expression",
        "original_parameter", "bound_source_parameter",
    }:
        result = True
    elif operation in {"sum", "product"}:
        args = row.get("arguments")
        result = isinstance(args, (list, tuple)) and bool(args) and all(
            _native_independent(nodes, item, coordinate_i, memo) for item in args
        )
    elif operation == "negative":
        result = _native_independent(nodes, row.get("argument"), coordinate_i, memo)
    elif operation in {"analytic_unary", "exact_original_source_flat_sigma"}:
        result = _native_independent(nodes, row.get("argument"), coordinate_i, memo)
    elif operation == "positive_quotient":
        result = _native_independent(nodes, row.get("numerator"), coordinate_i, memo) and _native_independent(
            nodes, row.get("denominator"), coordinate_i, memo
        )
    elif operation == "function_substitution":
        # The substituted variable is a source endpoint variable, not x in
        # the Rm radius factor.  Its value is an exact endpoint constant.
        result = _native_independent(nodes, row.get("expression"), coordinate_i, memo) and _native_independent(
            nodes, row.get("value"), coordinate_i, memo
        )
    else:
        result = False
    memo[index] = result
    return result


def _direct_local_contains(nodes: Any, root: Any, target: Any, depth: int = 5) -> bool:
    """Bounded operand reachability, avoiding unrelated archive metadata."""
    wanted = _handle(target)
    start = _handle(root)
    if wanted is None or start is None:
        return False
    frontier = [(start, 0)]
    seen: set[int] = set()
    while frontier:
        index, level = frontier.pop()
        if index == wanted:
            return True
        if level >= depth or index in seen:
            continue
        seen.add(index)
        row = _row(nodes, index)
        if row is None:
            continue
        for key in (
            "arguments", "argument", "numerator", "denominator", "expression",
            "active_body", "flat_value", "phase_function", "target_phase",
        ):
            value = row.get(key)
            values = value if isinstance(value, (list, tuple)) else (value,)
            for item in values:
                item_i = _handle(item)
                if item_i is not None:
                    frontier.append((item_i, level + 1))
    return False


def _definition_rows(binding: Mapping[str, Any], names: tuple[str, ...]) -> dict[str, bool]:
    definition = _mapping(binding.get("definition"))
    return {name: definition.get(name) is True for name in names}


def _common_n(nodes: Any, rm: Mapping[str, Any], loop: Mapping[str, Any]) -> dict[str, Any]:
    rm_n = _handle(rm.get("N"))
    loop_n = _handle(loop.get("N"))
    row = _row(nodes, rm_n)
    passed = bool(
        rm_n is not None
        and rm_n == loop_n
        and row
        and row.get("operation") == "exact_positive_integer_expression"
        and row.get("exact_positive_integer") is True
        and row.get("Z_independent") is True
        and row.get("one_common_N_for_all_charts_and_control_denominators") is True
    )
    return {
        "passed": passed,
        "node": rm_n,
        "rm_N": rm_n,
        "loop_N": loop_n,
        "node_row": row,
        "same_rm_and_loop_handle": rm_n == loop_n,
        "Z_independent": bool(row and row.get("Z_independent") is True),
    }


def _radius_geometry(
    nodes: Any,
    rm: Mapping[str, Any],
    loop: Mapping[str, Any],
    common_n: Any,
) -> dict[str, Any]:
    chart = _mapping(rm.get("chart"))
    coordinate = _handle(chart.get("coordinate"))
    radius = _handle(chart.get("radius"))
    jacobian = _handle(chart.get("physical_measure"))
    origin = _handle(rm.get("phase_origin"))
    if origin is None:
        origin = _handle(chart.get("phase_origin"))
    coord_row = _row(nodes, coordinate)
    radius_row = _row(nodes, radius)
    jac_row = _row(nodes, jacobian)
    coordinate_ok = bool(
        coord_row and coord_row.get("operation") == "bound_variable"
        and coord_row.get("name") == "true_Rm_patch_x"
    )

    radius_args = tuple(radius_row.get("arguments", ())) if radius_row else ()
    x_factors = [item for item in radius_args if _handle(item) == coordinate]
    other_factors = [item for item in radius_args if _handle(item) != coordinate]
    base_row = _row(nodes, other_factors[0]) if other_factors else None
    base_value = _handle(base_row.get("value")) if base_row else None
    base_value_row = _row(nodes, base_value)
    radius_factor_ok = bool(
        radius_row and radius_row.get("operation") == "product"
        and len(radius_args) == 2 and len(x_factors) == 1 and len(other_factors) == 1
        and base_row and base_row.get("operation") == "function_substitution"
        and base_row.get("Z_independent_substitution") is True
        and base_value_row and base_value_row.get("operation") == "exact_rational"
        and base_value_row.get("numerator") == 1 and base_value_row.get("denominator") == 1
        and _native_independent(nodes, other_factors[0], coordinate)
        and _row(nodes, other_factors[0]) is not None
    )

    numerator = _handle(jac_row.get("numerator")) if jac_row else None
    denominator = _handle(jac_row.get("denominator")) if jac_row else None
    numerator_row = _row(nodes, numerator)
    jacobian_ok = bool(
        jac_row and jac_row.get("operation") == "positive_quotient"
        and numerator == 1 and numerator_row
        and numerator_row.get("operation") == "exact_rational"
        and numerator_row.get("numerator") == 1
        and numerator_row.get("denominator") == 1
        and denominator == coordinate
        and jac_row.get("source_positive_certificate") == "positive original actual patch x"
    )

    # The offset is checked independently of the Jacobian.  In particular,
    # no constant-J helper is used here.
    offset = _handle(chart.get("global_offset"))
    offset_row = _row(nodes, offset)
    offset_args = tuple(offset_row.get("arguments", ())) if offset_row else ()
    log_nodes = []
    negative_nodes = []
    for item in offset_args:
        item_row = _row(nodes, item)
        if item_row and item_row.get("operation") == "analytic_unary" and item_row.get("name") == "log":
            log_nodes.append((_handle(item), item_row))
        if item_row and item_row.get("operation") == "negative" and _handle(item_row.get("argument")) == origin:
            negative_nodes.append(_handle(item))
    ratio = _row(nodes, log_nodes[0][1].get("argument")) if len(log_nodes) == 1 else None
    ratio_num = _handle(ratio.get("numerator")) if ratio else None
    ratio_den = _handle(ratio.get("denominator")) if ratio else None
    ratio_ok = bool(
        ratio and ratio.get("operation") == "positive_quotient"
        and ratio_num == radius and ratio_den is not None
    )
    # `_native_expr` deliberately treats the moving quotient as opaque.  The
    # denominator is required only to be a native-independent source graph;
    # the explicit x factor in R already supplies the derivative relation.
    ra_independent = False
    if ratio_den is not None:
        from .original_native_derivative_guards import _native_expr

        expression = _native_expr(nodes, ratio_den, coordinate, {})
        ra_independent = expression is not None and not expression[1]
    offset_ok = bool(
        offset_row and offset_row.get("operation") == "sum" and len(offset_args) == 2
        and len(log_nodes) == 1 and len(negative_nodes) == 1 and ratio_ok and ra_independent
        and origin is not None and _handle(loop.get("phase_origin")) == origin
    )

    phase = _phase_binding(nodes, chart, common_n)
    phase_ok = bool(phase.get("passed") is True)
    radius_rows = {
        "coordinate": coordinate,
        "radius": radius,
        "physical_measure": jacobian,
        "global_offset": offset,
        "phase_origin": origin,
        "coordinate_row": coord_row,
        "radius_row": radius_row,
        "jacobian_row": jac_row,
        "radius_base_factor": _handle(other_factors[0]) if other_factors else None,
        "radius_base_factor_row": base_row,
        "radius_base_substitution_at_native_Rm": bool(base_row and base_row.get("operation") == "function_substitution" and base_value == 1),
        "radius_has_exactly_one_native_x_factor": len(x_factors) == 1,
        "radius_base_factor_native_independent": radius_factor_ok,
        "native_J_is_constant": False,
        "native_J_is_varying": jacobian_ok,
        "native_J_denominator_is_x": denominator == coordinate,
        "d_y_equals_x_d_x": bool(radius_factor_ok and jacobian_ok),
        "R_equals_Rm_times_x": radius_factor_ok,
        "J_equals_one_over_x": jacobian_ok,
        "global_offset_is_log_R_over_Ra_minus_phase_origin": offset_ok,
        "Ra_node": ratio_den,
        "Ra_native_independent": ra_independent,
        "phase": phase,
        "phase_factor_after_J_division_is_N": bool(phase_ok and common_n is not None),
    }
    passed = bool(
        coordinate_ok and radius_factor_ok and jacobian_ok and offset_ok and phase_ok
        and common_n is not None
    )
    radius_rows.update({
        "passed": passed,
        "coordinate_is_actual_Rm_patch_x": coordinate_ok,
        "phase_origin_is_loop_origin": _handle(loop.get("phase_origin")) == origin,
    })
    return radius_rows


def _free_seed_binding(nodes: Any, leading: Mapping[str, Any]) -> dict[str, Any]:
    seeds = leading.get("free_source_seeds")
    rows = seeds if isinstance(seeds, (list, tuple)) else ()
    valid_rows = []
    gamma_x = []
    controls = []
    for item in rows:
        item_map = _mapping(item)
        node = _handle(item_map.get("node"))
        variables = item_map.get("variables")
        row = _row(nodes, node)
        valid = bool(node is not None and row is not None and isinstance(variables, list))
        valid_rows.append(valid)
        if valid and "true_Rm_patch_x" in variables:
            gamma_x.append({"node": node, "operation": row.get("operation"), "variables": list(variables)})
        if valid and row.get("operation") in {
            "original_Rm_unique_leading_implicit_control",
            "original_Rm_unique_leading_implicit_system",
        }:
            controls.append({"node": node, "operation": row.get("operation"), "component": row.get("component")})
    chart = _mapping(leading.get("chart"))
    radial = _mapping(chart.get("radial"))
    radial_ports = [
        _handle(item)
        for key in ("E_y", "V_y")
        for item in (radial.get(key) if isinstance(radial.get(key), (list, tuple)) else ())
    ]
    gamma_reachable = bool(gamma_x) and all(
        any(_contains(nodes, port, item["node"]) for port in radial_ports[:2]) for item in gamma_x
    )
    passed = bool(
        rows and all(valid_rows) and len(gamma_x) == 3 and len(controls) >= 6 and gamma_reachable
    )
    return {
        "passed": passed,
        "free_seed_rows": len(rows),
        "all_rows_well_formed": bool(rows) and all(valid_rows),
        "gamma_x_derivative_seed_rows": gamma_x,
        "gamma_x_seed_count": len(gamma_x),
        "gamma_x_seeds_reach_leading_radial_ports": gamma_reachable,
        "leading_control_seed_rows": controls,
        "leading_control_seed_count": len(controls),
        "source_values_selected": False,
    }


def _slow_binding(
    nodes: Any,
    slow_values: Any,
    source_values: Any,
    jacobian: Any,
    label: str,
) -> dict[str, Any]:
    slow = slow_values if isinstance(slow_values, (list, tuple)) else ()
    source = source_values if isinstance(source_values, (list, tuple)) else ()
    jacobian_i = _handle(jacobian)
    rows = []
    for order in _ORDERS:
        slow_i = _handle(slow[order]) if len(slow) > order else None
        source_i = _handle(source[order]) if len(source) > order else None
        slow_row = _row(nodes, slow_i)
        source_row = _row(nodes, source_i)
        args = tuple(slow_row.get("arguments", ())) if slow_row else ()
        inverse = next((_row(nodes, item) for item in args
                        if _row(nodes, item) and _row(nodes, item).get("operation") == "positive_quotient"), None)
        inverse_i = None
        if inverse is not None:
            inverse_i = next((_handle(item) for item in args if _row(nodes, item) is inverse), None)
        lazy = next((_row(nodes, item) for item in args
                     if _row(nodes, item) and _row(nodes, item).get("operation") == "original_lazy_flat_branch"), None)
        inverse_ok = bool(
            slow_row and slow_row.get("operation") == "product" and len(args) == 2
            and inverse and inverse.get("numerator") == 1
            and _handle(inverse.get("denominator")) == jacobian_i
            and inverse.get("source_positive_certificate") == "positive original chart log-radius Jacobian"
        )
        lazy_ok = bool(
            lazy and lazy.get("same_original_lazy_native_derivative") is True
            and lazy.get("ordinary_slow_Z_derivative_order") == order
            and lazy.get("active_body") is not None and lazy.get("flat_value") is not None
        )
        source_ok = bool(
            source_row and source_row.get("operation") == "original_lazy_flat_branch"
            and source_row.get("same_original_flat_extension_through_order5") is True
            and source_row.get("ordinary_slow_Z_derivative_order") == order
        )
        rows.append({
            "order": order,
            "slow_node": slow_i,
            "source_node": source_i,
            "inverse_J_node": inverse_i,
            "inverse_J_binding": inverse_ok,
            "generated_lazy_native_derivative_binding": lazy_ok,
            "source_lazy_binding": source_ok,
            "original_source_reachability_diagnostic": _contains(nodes, slow_i, source_i),
            "original_source_reachability_required": False,
        })
    return {
        "passed": bool(rows) and all(
            row["inverse_J_binding"] and row["generated_lazy_native_derivative_binding"]
            and row["source_lazy_binding"] for row in rows
        ),
        "label": label,
        "rows": rows,
        "all_six_ordinary_orders_checked": len(rows) == 6,
        "division_by_moving_J_once": bool(rows) and all(row["inverse_J_binding"] for row in rows),
    }


def _ast_guards(repo: Path) -> dict[str, Any]:
    leading = repo / _SOURCE_PINS["leading_patch"][0]
    finite = repo / _SOURCE_PINS["finite_patch"][0]
    native = repo / _SOURCE_PINS["native_derivative"][0]
    rows = {
        "native_constructor": _function_evidence(native, "NativeDerivative.__init__", (
            "self.cache={phase.node:g.zero}", "self.core_extensions=[]",
        )),
        "native_seed": _function_evidence(native, "NativeDerivative.seed", (
            "self.g.ids([left,right])", "raise ValueError('Conflicting native source derivative')",
            "raise ValueError('A constant cannot have a nonzero radial derivative')",
            "self.cache[left.node]=right",
        )),
        "native_derivative": _function_evidence(native, "NativeDerivative.derivative", (
            "g.quotient(g.sub(g.mul(D(u.node),v),g.mul(u,D(v.node))),g.mul(v,v),row['source_positive_certificate'])",
            "g.mul(outer,dx)", "g.add(direct,chain)",
            "attrs['active_body']=D(row['active_body']).node",
            "attrs['flat_value']=D(row['flat_value']).node",
            "attrs['same_original_lazy_native_derivative']=True",
        )),
        "background_radial_rows": _function_evidence(native, "background_radial_rows", (
            "result=dict(E=chart['radial']['E_y'],V=chart['radial']['V_y'],m=a.sub(V,f['m']),h=a.sub(E,a.scale(f['h'],'3/2')),k=a.sub(a.mul(E,V),a.scale(f['k'],'3/2')),e=a.sub(a.sub(a.mul(V,V),a.scale(E2,'1/2')),f['e']),p=a.scale(E2,'1/2'),P0=a.fixed(0))",
        )),
        "total_radial_functions": _function_evidence(native, "total_radial_functions", (
            "invjac=g.quotient(g.one,jac,'positive original chart log-radius Jacobian')",
            "slow={key:FunctionJet(tuple(g.mul(invjac,derivative.derivative(ref)) for ref in rows(chart['loop'][key]))) for key in ('A','B_over_Pstar')}",
            "a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
        )),
        "physical_pressure_fields": _function_evidence(native, "physical_pressure_fields", (
            "a.scale(radial['V_y'],Pstar)", "a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
        )),
        "finite_native_derivative": _function_evidence(finite, "native_derivative", (
            "D=radial.NativeDerivative(g,leading.x,phase)",
            "D.shared['free'].update({i:set(v) for i,v in leading.free_source_seeds.items()})",
            "D.seed(ref,g.mul(chart['physical_measure'],dr))",
            "D.seed(gamma[0],gamma[1])", "D.seed(gamma[1],gamma[2])",
        )),
        "leading_patch_fields": _function_evidence(leading, "patch_fields", (
            "Hx=a.add(a.fixed(g.mul(g.constant('1/10'),power(g,x,Fraction(-9,10)))), *(a.scale(h[i+2],gamma[i][1]) for i in range(3)))",
            "radial=dict(E_y=a.scale(a.mul(am,Hx),x),V_y=a.scale(axial_x,g.mul(x,inv)))",
        )),
    }
    source_for = {
        "native_constructor": native,
        "native_seed": native,
        "native_derivative": native,
        "background_radial_rows": native,
        "total_radial_functions": native,
        "physical_pressure_fields": native,
        "finite_native_derivative": finite,
        "leading_patch_fields": leading,
    }
    for name, row in rows.items():
        row["source"] = str(source_for[name].relative_to(repo)).replace("\\", "/")
    return {"passed": bool(rows) and all(row.get("passed") is True for row in rows.values()), "rows": rows}


def _source_definitions(rm: Mapping[str, Any], leading: Mapping[str, Any], loop: Mapping[str, Any]) -> dict[str, Any]:
    leading_names = (
        "actual_E_V_m_h_k_e_p_P0_pressure_and_Ey_Vy_C6",
        "original_x1_inlet_and_Rh_exp1_exit_and_log_Jacobian_1_over_x",
        "full_signed_generic_C5",
        "leading_control_not_finite_N_correction_and_no_history_reset",
        "original_unique_common_box_branch_same_full_Jacobian_through_Z6",
    )
    rm_names = (
        "actual_same_N_Section11_loop_inverse_density_and_history_C5",
        "full_fast_and_slow_total_radial_C5_and_modified_generic_C4",
        "original_gamma_first_second_and_same_implicit_control_native_sources",
        "physical_U_Uy_mean_C5_Q_Qy_C4_and_signed_rminus_to_Rh_pressure_C3",
        "exact_physical_Rm_and_Rm_squared_x_pressure_measures_once",
    )
    loop_names = (
        "one_original_common_N_eta_dstar_and_global_log_R_rminus_phase",
        "complete_original_rminus_to_Rm_Duhamel_correction_histories_C5",
        "correction_zero_inlet_only_at_rminus_background_nonzero_memories_unchanged",
        "original_rate_zero_pressure_memory_retained",
    )
    rows = {
        "leading": _definition_rows(leading, leading_names),
        "finite_Rm": _definition_rows(rm, rm_names),
        "loop": _definition_rows(loop, loop_names),
    }
    return {
        "passed": all(all(group.values()) for group in rows.values()),
        "rows": rows,
    }


def build_Rm_native_guard(
    nodes: Any,
    rm_bindings: Mapping[str, Any],
    leading_bindings: Mapping[str, Any],
    loop_bindings: Mapping[str, Any],
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Return fail-closed source evidence for the actual Rm native calculus.

    ``rm_bindings`` is the finite-N C5 snapshot, ``leading_bindings`` is the
    true leading C6 snapshot, and ``loop_bindings`` is the prescribed loop
    snapshot.  Only their serialized handles are inspected.
    """
    repo = _repo_root(root)
    rm = _unwrap(rm_bindings, _RM_KEY)
    leading = _unwrap(leading_bindings, _LEADING_KEY)
    loop = _unwrap(loop_bindings, _LOOP_KEY)
    chart = _mapping(rm.get("chart"))
    missing = [
        name for name, value in (
            ("rm.chart", chart), ("leading.chart", _mapping(leading.get("chart"))),
            ("loop.charts", _mapping(loop.get("charts"))),
        ) if not value
    ]
    if missing:
        return {
            "passed": False,
            "method": "actual_Rm_C5_moving_native_calculus_source_guard",
            "source_only": True,
            "numeric_caps_used": False,
            "ancestor_replayed": False,
            "unproved": ["missing required serialized binding: " + ", ".join(missing)],
        }

    pins = {name: _pin(repo, relative, expected) for name, (relative, expected) in _SOURCE_PINS.items()}
    definitions = _source_definitions(rm, leading, loop)
    common_n = _common_n(nodes, rm, loop)
    geometry = _radius_geometry(nodes, rm, loop, common_n.get("node"))
    seeds = _free_seed_binding(nodes, leading)
    ast_guards = _ast_guards(repo)

    radial_background = _mapping(rm.get("radial_background"))
    leading_radial = _mapping(_mapping(leading.get("chart")).get("radial"))
    e_seed = radial_background.get("E") == leading_radial.get("E_y") and _sequence(radial_background.get("E"), 7)
    v_seed = radial_background.get("V") == leading_radial.get("V_y") and _sequence(radial_background.get("V"), 7)
    background_seed = {
        "passed": bool(e_seed and v_seed and ast_guards["rows"]["finite_native_derivative"].get("passed")),
        "E_rows": list(radial_background.get("E", ())),
        "V_rows": list(radial_background.get("V", ())),
        "leading_E_y_rows": list(leading_radial.get("E_y", ())),
        "leading_V_y_rows": list(leading_radial.get("V_y", ())),
        "E_seed_is_physical_leading_row": e_seed,
        "V_seed_is_physical_leading_row": v_seed,
        "native_seed_formula": "D.seed(ref, J * physical_row)",
        "J_times_physical_rows_source": ast_guards["rows"]["finite_native_derivative"],
    }

    radial = _mapping(rm.get("radial"))
    loop_chart = _mapping(chart.get("loop"))
    slow_A = _slow_binding(nodes, radial.get("slow_A_y"), loop_chart.get("A"), chart.get("physical_measure"), "A")
    slow_B = _slow_binding(nodes, radial.get("slow_B_y"), loop_chart.get("B_over_Pstar"), chart.get("physical_measure"), "B_over_Pstar")

    v_rows = []
    vy = radial.get("V_y")
    bphi = loop_chart.get("B_over_Pstar_phi")
    for order in (0, 1):
        vy_i = _handle(vy[order]) if isinstance(vy, (list, tuple)) and len(vy) > order else None
        slow_i = _handle(radial.get("slow_B_y", ())[order]) if isinstance(radial.get("slow_B_y"), (list, tuple)) and len(radial.get("slow_B_y")) > order else None
        bphi_i = _handle(bphi[order]) if isinstance(bphi, (list, tuple)) and len(bphi) > order else None
        v_rows.append({
            "ordinary_Z_order": order,
            "V_y_node": vy_i,
            "slow_B_y_node": slow_i,
            "B_phi_node": bphi_i,
            "contains_slow_B_y": _direct_local_contains(nodes, vy_i, slow_i),
            "contains_B_phi": _direct_local_contains(nodes, vy_i, bphi_i),
        })
    v_rows_passed = bool(v_rows) and all(row["contains_slow_B_y"] and row["contains_B_phi"] for row in v_rows)

    physical = _mapping(rm.get("physical"))
    uy = physical.get("U_y")
    u = physical.get("U")
    # Pstar is not a top-level field in this snapshot.  Identify it only as
    # the exact common direct factor in U and U_y, never by evaluating it.
    product_args = []
    for values in (u, uy):
        for item in values[:2] if isinstance(values, (list, tuple)) else ():
            row = _row(nodes, item)
            product_args.append(set(_handle(arg) for arg in row.get("arguments", ())) if row and row.get("operation") == "product" else set())
    pstar_candidates = set.intersection(*product_args) if product_args else set()
    pstar = next(iter(pstar_candidates), None) if len(pstar_candidates) == 1 else None
    uy_factor = _direct_factor(nodes, uy, pstar)
    boundary = _mapping(rm.get("boundary"))
    rqy = boundary.get("R_Q_y")
    rqy_factor = _direct_factor(nodes, rqy, chart.get("radius"))
    physical_ports = {
        "passed": bool(uy_factor.get("passed") and rqy_factor.get("passed")),
        "Pstar_node_inferred_from_U_and_U_y": pstar,
        "U_y_Pstar_factor": uy_factor,
        "R_Q_y_radius_factor": rqy_factor,
        "physical_U_y_has_one_Pstar": uy_factor.get("passed") is True,
        "physical_R_Q_y_has_one_R": rqy_factor.get("passed") is True,
    }

    chart_keys = all(key in chart for key in (
        "coordinate", "domain", "radius", "physical_measure", "global_offset", "actual_phase", "loop",
    ))
    radial_keys = all(key in radial for key in ("E_y", "V_y", "slow_A_y", "slow_B_y"))
    ordinary_rows = {
        "passed": chart_keys and radial_keys and _sequence(vy) and _sequence(uy) and _sequence(rqy),
        "chart_keys_present": chart_keys,
        "radial_keys_present": radial_keys,
        "ordinary_Z_orders": [0, 1],
        "V_y_rows": v_rows,
    }
    checks = {
        "source_pins": {"passed": all(row["passed"] for row in pins.values()), "rows": pins},
        "source_definitions": definitions,
        "source_ast": ast_guards,
        "common_positive_N": common_n,
        "moving_Rm_radius_and_J": geometry,
        "leading_free_controls_and_gamma_x": seeds,
        "background_radial_seed_binding": background_seed,
        "slow_A_native_partial_divided_by_J": slow_A,
        "slow_B_native_partial_divided_by_J": slow_B,
        "total_V_y_slow_and_fast_rows": {"passed": v_rows_passed, "rows": v_rows},
        "physical_U_y_and_R_Q_y_factors": physical_ports,
        "ordinary_Z_rows": ordinary_rows,
    }
    passed = all(_mapping(row).get("passed") is True for row in checks.values())
    return {
        "passed": passed,
        "method": "actual_Rm_C5_moving_native_calculus_source_guard",
        "source_only": True,
        "numeric_caps_used": False,
        "source_values_selected": False,
        "ancestor_replayed": False,
        "checks": checks,
        "source_pins": pins,
        "scope": {
            "chart": "Rm_patch",
            "native_coordinate": geometry.get("coordinate"),
            "native_radius": geometry.get("radius"),
            "native_physical_measure": geometry.get("physical_measure"),
            "common_N": common_n.get("node"),
            "phase_origin": geometry.get("phase_origin"),
            "domain": chart.get("domain"),
        },
        "premises": {
            "R_equals_Rm_times_x": geometry.get("R_equals_Rm_times_x") is True,
            "J_equals_one_over_x": geometry.get("J_equals_one_over_x") is True,
            "d_y_equals_x_d_x": geometry.get("d_y_equals_x_d_x") is True,
            "native_J_is_varying": geometry.get("native_J_is_varying") is True,
            "native_J_is_constant": False,
            "global_offset_log_R_over_Ra_minus_phase_origin": geometry.get("global_offset_is_log_R_over_Ra_minus_phase_origin") is True,
            "phase_uses_common_N": geometry.get("phase", {}).get("passed") is True,
            "native_background_seed_is_J_times_physical_row": background_seed.get("passed") is True,
            "slow_derivative_divides_by_moving_J_once": slow_A.get("passed") is True and slow_B.get("passed") is True,
            "physical_U_y_has_one_Pstar": physical_ports.get("physical_U_y_has_one_Pstar") is True,
            "physical_R_Q_y_has_one_radius": physical_ports.get("physical_R_Q_y_has_one_R") is True,
            "gamma_x_derivatives_seeded": seeds.get("passed") is True,
        },
        "unproved": {
            "numeric_endpoint_value": True,
            "error_cap_as_source_value": True,
            "ancestor_constructor_or_producer_replay": True,
            "endpoint_seam_identity": True,
            "constant_native_J_identity": True,
        },
        "limits": {
            "physical_endpoint_seam": False,
            "global_m15_contraction": False,
            "numeric_source_containment": False,
            "constant_J_assumption": False,
        },
    }


__all__ = ["build_Rm_native_guard"]
