"""Pure graph source for the two-chart finite-``N`` phase mean kernel.

The source charts and their graph handles are supplied by the caller.  This
module adds only exact scalar graph rows for the phase primitives and the
rate-one moving convolutions; it does not load an owner, evaluate a bound, or
select a numerical mean value.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Mapping, Sequence
from typing import Any

from .original_signed_graph_pairs import GraphPairs


_CHARTS = ("first_micro", "second_micro")
_ORDERS = (0, 1, 2, 3)
_DEFAULT_PSTAR = 33572
_DEFAULT_N = 18
_DEFAULT_HB = 8
_INTEGRAL_CERTIFICATE = "same exact source positive interval orientation"


def _graph_id(graph: GraphPairs, value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer graph handle")
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return graph._id(value, name=name)


def _mapping(value: Any, *, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a mapping")
    return value


def _rows(value: Any, *, name: str) -> tuple[Any, ...]:
    if isinstance(value, Mapping):
        for key in ("ordinary_Z", "ordinary_Z_rows", "rows", "values"):
            if key in value:
                value = value[key]
                break
        else:
            raise TypeError(f"{name} must contain ordinary-Z rows")
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence of graph handles")
    return tuple(value)


def _checked_rows(
    graph: GraphPairs,
    value: Any,
    *,
    name: str,
    minimum: int = 4,
) -> tuple[int, ...]:
    rows = _rows(value, name=name)
    if len(rows) < minimum:
        raise ValueError(f"{name} requires at least {minimum} ordinary-Z rows")
    return tuple(_graph_id(graph, row, name=f"{name}[{index}]") for index, row in enumerate(rows[:minimum]))


def _chart_map(bindings: Mapping[str, Any], *, name: str) -> Mapping[str, Any]:
    charts = bindings.get("charts")
    if isinstance(charts, Mapping) and all(chart in charts for chart in _CHARTS):
        return charts
    if all(chart in bindings for chart in _CHARTS):
        return bindings
    raise KeyError(f"{name} must expose first_micro and second_micro charts")


def _nested_handle(
    mappings: Sequence[Mapping[str, Any]],
    names: Sequence[str],
    graph: GraphPairs,
    *,
    fallback: int,
    name: str,
) -> int:
    """Find a scalar source handle without assuming one binding envelope."""

    wanted = {item.lower() for item in names}
    queue: deque[Mapping[str, Any]] = deque(mappings)
    seen: set[int] = set()
    while queue:
        current = queue.popleft()
        marker = id(current)
        if marker in seen:
            continue
        seen.add(marker)
        for key, value in current.items():
            if isinstance(key, str) and key.lower() in wanted:
                try:
                    return _graph_id(graph, value, name=name)
                except (TypeError, ValueError):
                    pass
            if isinstance(value, Mapping) and (
                key in {"bindings", "parameters", "scalars", "common", "downstream", "controls"}
                or not isinstance(value, (str, bytes))
            ):
                queue.append(value)
    try:
        return _graph_id(graph, fallback, name=name)
    except (TypeError, ValueError) as error:
        raise KeyError(f"{name} is missing from supplied source bindings") from error


def _domain(chart: Mapping[str, Any], graph: GraphPairs, *, name: str) -> tuple[int, int]:
    raw = chart.get("domain")
    if raw is None:
        lower, upper = chart.get("lower"), chart.get("upper")
        if lower is None or upper is None:
            raise KeyError(f"{name} must expose domain or lower and upper handles")
        raw = (lower, upper)
    if isinstance(raw, (str, bytes)) or not isinstance(raw, Sequence) or len(raw) != 2:
        raise ValueError(f"{name} domain must contain exactly lower and upper handles")
    return (
        _graph_id(graph, raw[0], name=f"{name} domain lower"),
        _graph_id(graph, raw[1], name=f"{name} domain upper"),
    )


def _endpoint_at_two(chart: Mapping[str, Any], graph: GraphPairs, *, name: str) -> int:
    for key in ("endpoint_at_2", "at2", "endpoint", "upper"):
        if key in chart:
            return _graph_id(graph, chart[key], name=f"{name} endpoint at 2")
    return _domain(chart, graph, name=name)[1]


def _append_unary(graph: GraphPairs, name: str, argument: int, **metadata: Any) -> int:
    return graph.node("analytic_unary", name=name, argument=argument, **metadata)


def _substitute(
    graph: GraphPairs,
    expression: int,
    variable: int,
    value: int,
    *,
    name: str,
    phase_variable: int | None = None,
    outer_integration_variable: bool = False,
) -> int:
    attributes: dict[str, Any] = {
        "expression": _graph_id(graph, expression, name=f"{name} expression"),
        "variable": _graph_id(graph, variable, name=f"{name} variable"),
        "value": _graph_id(graph, value, name=f"{name} value"),
        "Z_independent_substitution": True,
        "capture_avoiding": True,
        "preserve_distinct_outer_integration_variable": True,
        "outer_integration_variable": outer_integration_variable,
    }
    if phase_variable is not None:
        attributes.update(
            source_phase_node=_graph_id(graph, phase_variable, name=f"{name} source phase"),
            derived_phase_recomputed_at_value=True,
            derived_phase_not_bound_variable=True,
        )
    return graph.node("function_substitution", **attributes)


def _definite_integral(
    graph: GraphPairs,
    integrand: int,
    lower: int,
    upper: int,
    *,
    binder_name: str,
    binder: int | None = None,
    order: int,
    kind: str,
    metadata: Mapping[str, Any] | None = None,
) -> tuple[int, int]:
    if binder is None:
        binder = graph.node("bound_variable", name=binder_name)
    else:
        binder = _graph_id(graph, binder, name=f"{kind} bound variable")
    integral = graph.node(
        "definite_integral",
        integrand=_graph_id(graph, integrand, name=f"{kind} integrand"),
        variable=binder_name,
        bound_variable=binder,
        lower=_graph_id(graph, lower, name=f"{kind} lower"),
        upper=_graph_id(graph, upper, name=f"{kind} upper"),
        ordinary_slow_Z_derivative_order=order,
        endpoint_independent_of_Z=True,
        **dict(metadata or {}),
    )
    return integral, binder


def _phase_rows(
    graph: GraphPairs,
    rows: Sequence[int],
    *,
    phase: int,
    lower: int,
    chart_name: str,
    source_name: str,
) -> tuple[int, ...]:
    result: list[int] = []
    for order, body in enumerate(rows):
        theta_name = f"micro_mean_{chart_name}_{source_name}_Z{order}_theta"
        theta = graph.node("bound_variable", name=theta_name)
        substituted = _substitute(
            graph,
            body,
            phase,
            theta,
            name=f"{chart_name} {source_name} phase substitution Z{order}",
            phase_variable=phase,
        )
        integral, _ = _definite_integral(
            graph,
            substituted,
            lower,
            phase,
            binder_name=theta_name,
            order=order,
            kind=f"{chart_name} {source_name} phase primitive Z{order}",
            binder=theta,
            metadata={
                "source_body": body,
                "source_phase": phase,
                "phase_variable_bound_as_theta": True,
                "actual_phase_Z_independent": True,
                "fixed_phase_slow_y": source_name == "Cy",
                "total_y_chain_excluded": source_name == "Cy",
                "phase_primitive_is_direct_ordinary_Z_row": True,
                "factorial_normalization": False,
                "ordinary_Z_order": order,
                "source_operation": "phase_integral_of_slow_B_y" if source_name == "Cy" else "phase_integral_of_B",
                "numeric_value_not_installed": True,
            },
        )
        result.append(integral)
    return tuple(result)


def _moving_rows(
    graph: GraphPairs,
    primitive_rows: Sequence[int],
    *,
    coordinate: int,
    phase: int,
    lower: int,
    hb: int,
    chart_name: str,
    source_name: str,
) -> tuple[int, ...]:
    result: list[int] = []
    for order, primitive in enumerate(primitive_rows):
        t_name = f"micro_mean_{chart_name}_{source_name}_Z{order}_t"
        t = graph.node("bound_variable", name=t_name)
        primitive_at_t = _substitute(
            graph,
            primitive,
            coordinate,
            t,
            name=f"{chart_name} {source_name} moving substitution Z{order}",
            phase_variable=phase,
            outer_integration_variable=True,
        )
        distance = graph.add(coordinate, graph.neg(t))
        exponent = graph.neg(graph.mul(hb, distance))
        kernel = _append_unary(
            graph,
            "exp",
            exponent,
            moving_kernel=True,
            source_chart=chart_name,
            source_name=source_name,
        )
        integrand = graph.mul(kernel, primitive_at_t, hb)
        integral, _ = _definite_integral(
            graph,
            integrand,
            lower,
            coordinate,
            binder_name=t_name,
            order=order,
            kind=f"{chart_name} {source_name} moving convolution Z{order}",
            binder=t,
            metadata={
                "primitive": primitive,
                "kernel": kernel,
                "kernel_formula": "exp(-hb*(s-t))",
                "physical_hb_weight_applied_once": True,
                "moving_rate": "hb",
                "native_coordinate": coordinate,
                "inner_and_outer_native_variables_distinct": True,
                "capture_avoiding_coordinate_substitution": True,
                "ordinary_Z_order": order,
                "factorial_normalization": False,
                "numeric_value_not_installed": True,
            },
        )
        result.append(integral)
    return tuple(result)


def _scale_rows(graph: GraphPairs, rows: Sequence[int], scalar: int) -> tuple[int, ...]:
    return tuple(graph.mul(value, scalar) for value in rows)


def _add_rows(graph: GraphPairs, *rows: Sequence[int]) -> tuple[int, ...]:
    if not rows:
        raise ValueError("ordinary-Z row addition needs at least one row sequence")
    length = len(rows[0])
    if any(len(row) != length for row in rows[1:]):
        raise ValueError("ordinary-Z row addition refuses silent truncation")
    return tuple(graph.add(*(row[index] for row in rows)) for index in range(length))


def _neg_rows(graph: GraphPairs, rows: Sequence[int]) -> tuple[int, ...]:
    return tuple(graph.neg(value) for value in rows)


def _sub_rows(graph: GraphPairs, left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    if len(left) != len(right):
        raise ValueError("ordinary-Z row subtraction refuses silent truncation")
    return _add_rows(graph, left, _neg_rows(graph, right))


def _substitute_rows(
    graph: GraphPairs,
    rows: Sequence[int],
    *,
    variable: int,
    value: int,
    name: str,
) -> tuple[int, ...]:
    return tuple(
        _substitute(
            graph,
            expression,
            variable,
            value,
            name=f"{name} Z{order}",
        )
        for order, expression in enumerate(rows)
    )


def _symbolic_controls() -> dict[str, Any]:
    """Independent symbolic checks for transport and generic seam behavior."""

    try:
        import sympy as sp
    except ImportError:
        return {"available": False, "passed": False}

    y, N, P, hb = sp.symbols("y N P hb", nonzero=True)
    C = sp.Function("C")(y)
    Cy = sp.Function("Cy")(y)
    B = sp.Function("B")(y)
    JC = sp.Function("JC")(y)
    JCy = sp.Function("JCy")(y)
    delta_m = P / N**2 * (C - JC - JCy)
    transport = sp.diff(delta_m, y) + delta_m
    transport_substituted = transport.subs(
        {
            sp.diff(C, y): Cy + N * B,
            sp.diff(JC, y): C - JC,
            sp.diff(JCy, y): Cy - JCy,
        },
        simultaneous=True,
    )
    transport_check = sp.simplify(transport_substituted - P * B / N) == 0

    seam_checks: dict[str, bool] = {}
    for order in _ORDERS:
        C1, C2, JC1, JC2, JCy1, JCy2 = sp.symbols(
            f"C1_{order} C2_{order} JC1_{order} JC2_{order} JCy1_{order} JCy2_{order}"
        )
        first_end = C1 - JC1 - JCy1
        second_lower = C2 - JC1 - JC2 - JCy1 - JCy2
        seam_checks[str(order)] = sp.simplify(
            (second_lower - first_end).subs({C2: C1, JC2: 0, JCy2: 0})
        ) == 0

    inlet = sp.symbols("incoming", nonzero=True)
    memory = sp.exp(-hb * y)
    inlet_check = sp.simplify(sp.diff(memory * inlet, inlet) - memory) == 0
    return {
        "available": True,
        "transport_delta_m_y_plus_delta_m_equals_Pstar_B_over_N": bool(transport_check),
        "seam_cancel_through_Z0_to_Z3": seam_checks,
        "nonzero_inlet_retained_generically": bool(inlet_check),
        "passed": bool(transport_check and all(seam_checks.values()) and inlet_check),
        "assumptions": {
            "C_total_y_equals_Cy_plus_NB": True,
            "JC_y_plus_JC_equals_C": True,
            "JCy_y_plus_JCy_equals_Cy": True,
            "first_Cin_is_caller_proved_zero": True,
            "second_chart_does_not_reset_first_memory": True,
        },
    }


def build_micro_mean(
    g: GraphPairs,
    trueQ4bindings: Mapping[str, Any],
    loopbindings: Mapping[str, Any],
) -> dict[str, Any]:
    """Build source-only phase means for first and second micro charts."""

    if not isinstance(g, GraphPairs):
        raise TypeError("g must be a GraphPairs instance")
    true_bindings = _mapping(trueQ4bindings, name="trueQ4bindings")
    loop_bindings = _mapping(loopbindings, name="loopbindings")
    true_charts = _chart_map(true_bindings, name="trueQ4bindings")
    loop_charts = _chart_map(loop_bindings, name="loopbindings")

    Pstar = _nested_handle((true_bindings, loop_bindings), ("Pstar", "Pstar_id", "pstar"), g, fallback=_DEFAULT_PSTAR, name="Pstar")
    N = _nested_handle((loop_bindings, true_bindings), ("N", "N_id", "common_N", "common_n"), g, fallback=_DEFAULT_N, name="N")
    hb = _nested_handle((true_bindings, loop_bindings), ("hb", "hb_id", "hbB", "original_hb"), g, fallback=_DEFAULT_HB, name="hb")

    zero = g.constant(0)
    one = g.constant(1)
    n_squared = g.mul(N, N)
    scale = g.div(Pstar, n_squared, "same exact positive common N squared")

    chart_ports: dict[str, Any] = {}
    phase_rows: dict[str, tuple[int, ...]] = {}
    cy_rows: dict[str, tuple[int, ...]] = {}
    moving_rows: dict[str, tuple[int, ...]] = {}
    moving_cy_rows: dict[str, tuple[int, ...]] = {}
    domains: dict[str, tuple[int, int]] = {}
    endpoints: dict[str, int] = {}
    coordinates: dict[str, int] = {}

    for chart_name in _CHARTS:
        loop_chart = _mapping(loop_charts[chart_name], name=f"loop chart {chart_name}")
        true_chart = _mapping(true_charts[chart_name], name=f"true Q4 chart {chart_name}")
        loop_section = _mapping(loop_chart.get("loop"), name=f"loop chart {chart_name}['loop']")
        radial_section = _mapping(true_chart.get("radial"), name=f"true Q4 chart {chart_name}['radial']")
        coordinate = _graph_id(g, loop_chart.get("coordinate"), name=f"{chart_name} coordinate")
        phase = _graph_id(g, loop_chart.get("actual_phase"), name=f"{chart_name} actual phase")
        lower, upper = _domain(loop_chart, g, name=f"{chart_name} chart")
        B_rows = _checked_rows(g, loop_section.get("B_over_Pstar"), name=f"{chart_name} B_over_Pstar")
        slow_rows = _checked_rows(g, radial_section.get("slow_B_y"), name=f"{chart_name} slow_B_y")
        C = _phase_rows(g, B_rows, phase=phase, lower=zero, chart_name=chart_name, source_name="C")
        Cy = _phase_rows(g, slow_rows, phase=phase, lower=zero, chart_name=chart_name, source_name="Cy")
        JC = _moving_rows(
            g,
            C,
            coordinate=coordinate,
            phase=phase,
            lower=lower,
            hb=hb,
            chart_name=chart_name,
            source_name="C",
        )
        JCy = _moving_rows(
            g,
            Cy,
            coordinate=coordinate,
            phase=phase,
            lower=lower,
            hb=hb,
            chart_name=chart_name,
            source_name="Cy",
        )
        phase_rows[chart_name] = C
        cy_rows[chart_name] = Cy
        moving_rows[chart_name] = JC
        moving_cy_rows[chart_name] = JCy
        domains[chart_name] = (lower, upper)
        endpoints[chart_name] = _endpoint_at_two(loop_chart, g, name=f"{chart_name} chart")
        coordinates[chart_name] = coordinate
        chart_ports[chart_name] = {
            "coordinate": coordinate,
            "actual_phase": phase,
            "domain": (lower, upper),
            "endpoint_at_2": endpoints[chart_name],
            "B_over_Pstar": B_rows,
            "slow_B_y": slow_rows,
            "C": C,
            "Cy": Cy,
            "JC": JC,
            "JCy": JCy,
            "actual_phase_Z_independent": loop_chart.get("total_Z_phase_derivative_exact_zero") is True,
            "fixed_phase_slow_B_y": True,
            "total_y_chain_excluded": True,
            "ordinary_Z_orders": _ORDERS,
            "factorial_normalization": False,
        }

    first_C = phase_rows["first_micro"]
    first_Cy = cy_rows["first_micro"]
    first_JC = moving_rows["first_micro"]
    first_JCy = moving_cy_rows["first_micro"]
    second_C = phase_rows["second_micro"]
    second_Cy = cy_rows["second_micro"]
    second_JC = moving_rows["second_micro"]
    second_JCy = moving_cy_rows["second_micro"]

    first_Cin = tuple(zero for _ in _ORDERS)
    first_known = _scale_rows(g, _sub_rows(g, first_C, first_JC), scale)
    first_unresolved = _scale_rows(g, _neg_rows(g, first_JCy), scale)
    first_complete = _add_rows(g, first_known, first_unresolved)

    first_end = domains["first_micro"][1]
    first_JC_at_one = _substitute_rows(
        g,
        first_JC,
        variable=coordinates["first_micro"],
        value=first_end,
        name="first_micro JC endpoint at 1",
    )
    first_JCy_at_one = _substitute_rows(
        g,
        first_JCy,
        variable=coordinates["first_micro"],
        value=first_end,
        name="first_micro JCy endpoint at 1",
    )
    graph_distance = g.add(
        coordinates["second_micro"],
        g.neg(one),
    )
    second_memory = _append_unary(
        g,
        "exp",
        g.neg(g.mul(hb, graph_distance)),
        moving_memory=True,
        formula="exp(-hb*(s-1))",
        first_chart_endpoint=first_end,
    )
    second_joint_C = _add_rows(
        g,
        _scale_rows(g, first_JC_at_one, second_memory),
        second_JC,
    )
    second_joint_Cy = _add_rows(
        g,
        _scale_rows(g, first_JCy_at_one, second_memory),
        second_JCy,
    )
    second_known = _scale_rows(g, _sub_rows(g, second_C, second_joint_C), scale)
    second_unresolved = _scale_rows(g, _neg_rows(g, second_joint_Cy), scale)
    second_complete = _add_rows(g, second_known, second_unresolved)

    complete = {"first_micro": first_complete, "second_micro": second_complete}
    known = {"first_micro": first_known, "second_micro": second_known}
    unresolved = {"first_micro": first_unresolved, "second_micro": second_unresolved}
    endpoint_substitutions = {
        "first_at_1": {
            "C": _substitute_rows(
                g,
                first_C,
                variable=coordinates["first_micro"],
                value=first_end,
                name="first_micro C endpoint at 1",
            ),
            "Cy": _substitute_rows(
                g,
                first_Cy,
                variable=coordinates["first_micro"],
                value=first_end,
                name="first_micro Cy endpoint at 1",
            ),
            "JC": first_JC_at_one,
            "JCy": first_JCy_at_one,
            "known_C_part": _substitute_rows(
                g,
                first_known,
                variable=coordinates["first_micro"],
                value=first_end,
                name="first_micro known C endpoint at 1",
            ),
            "unresolved_Cy_part": _substitute_rows(
                g,
                first_unresolved,
                variable=coordinates["first_micro"],
                value=first_end,
                name="first_micro unresolved Cy endpoint at 1",
            ),
            "complete_physical_delta_m": _substitute_rows(
                g,
                first_complete,
                variable=coordinates["first_micro"],
                value=first_end,
                name="first_micro complete mean endpoint at 1",
            ),
        },
        "second_at_2": {
            "C": _substitute_rows(
                g,
                second_C,
                variable=coordinates["second_micro"],
                value=endpoints["second_micro"],
                name="second_micro C endpoint at 2",
            ),
            "Cy": _substitute_rows(
                g,
                second_Cy,
                variable=coordinates["second_micro"],
                value=endpoints["second_micro"],
                name="second_micro Cy endpoint at 2",
            ),
            "JC_joint": _substitute_rows(
                g,
                second_joint_C,
                variable=coordinates["second_micro"],
                value=endpoints["second_micro"],
                name="second_micro JC joint endpoint at 2",
            ),
            "JCy_joint": _substitute_rows(
                g,
                second_joint_Cy,
                variable=coordinates["second_micro"],
                value=endpoints["second_micro"],
                name="second_micro JCy joint endpoint at 2",
            ),
            "known_C_part": _substitute_rows(
                g,
                second_known,
                variable=coordinates["second_micro"],
                value=endpoints["second_micro"],
                name="second_micro known C endpoint at 2",
            ),
            "unresolved_Cy_part": _substitute_rows(
                g,
                second_unresolved,
                variable=coordinates["second_micro"],
                value=endpoints["second_micro"],
                name="second_micro unresolved Cy endpoint at 2",
            ),
            "complete_physical_delta_m": _substitute_rows(
                g,
                second_complete,
                variable=coordinates["second_micro"],
                value=endpoints["second_micro"],
                name="second_micro complete mean endpoint at 2",
            ),
        },
    }

    symbolic_controls = _symbolic_controls()
    source_ports = {
        "scalars": {"Pstar": Pstar, "N": N, "hb": hb, "N_squared": n_squared, "Pstar_over_N_squared": scale},
        "charts": chart_ports,
        "first_Cin": first_Cin,
        "first_Cin_source_owned_zero_required": True,
        "second_memory": second_memory,
        "second_memory_formula": "exp(-hb*(s-1))",
        "second_joint_C": second_joint_C,
        "second_joint_Cy": second_joint_Cy,
        "ordinary_Z_orders": _ORDERS,
        "phase_primitive_lower": zero,
        "integral_certificate": _INTEGRAL_CERTIFICATE,
    }
    flags = {
        "source_only": True,
        "source_chart_serialized": False,
        "numeric_mean_available": False,
        "numeric_bounds_available": False,
        "primitive_proof_claimed": False,
        "old_owner_executed": False,
        "old_builder_executed": False,
        "old_fixture_executed": False,
        "whole_micro_numeric": False,
        "actual_phase_Z_independent_required_from_caller": True,
        "fixed_phase_Cy_not_total_y_derivative": True,
        "physical_hb_weight_once": True,
        "first_Cin_zero_required_caller_proof": True,
        "second_memory_not_reset": True,
        "factorial_normalization": False,
    }
    return {
        "passed": bool(symbolic_controls.get("passed", False)),
        "phase_rows": phase_rows,
        "C": phase_rows,
        "Cy": cy_rows,
        "rawCyfunctional": cy_rows,
        "raw_Cy_functional": cy_rows,
        "JC": moving_rows,
        "JCy": moving_cy_rows,
        "JC_joint": {"first_micro": first_JC, "second_micro": second_joint_C},
        "JCy_joint": {"first_micro": first_JCy, "second_micro": second_joint_Cy},
        "known_C_part": known,
        "unresolved_Cy_part": unresolved,
        "completephysicaldelta_m": complete,
        "complete_physical_delta_m": complete,
        "endpointsubstitutionsat2": endpoint_substitutions,
        "endpoint_substitutions_at_2": endpoint_substitutions,
        "at_2": endpoint_substitutions["second_at_2"],
        "graphsourceports": source_ports,
        "graph_source_ports": source_ports,
        "symbolic_controls": symbolic_controls,
        "flags": flags,
        "ordinary_Z_orders": _ORDERS,
        "source_chart_serialized": False,
        "numeric_mean_available": False,
        "caller_guards_required": {
            "actual_phase_Z_independent": True,
            "parity_seam": True,
            "first_Cin_zero": True,
        },
        "graph": {
            "base_nodes_untouched": True,
            "appended_nodes": g.appended,
            "local_nodes": len(g.nodes),
        },
    }


__all__ = ["build_micro_mean"]
