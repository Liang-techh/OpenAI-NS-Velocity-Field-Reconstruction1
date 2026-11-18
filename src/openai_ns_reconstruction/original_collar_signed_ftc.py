"""Fresh source-only radial FTC leaf for the lower-collar signed source.

The caller supplies actual ordinary-Z C1 source ports.  The leaf never obtains
a field from a cap, midpoint, old numerical row, or an interval overlap.  It
keeps the small ``epsilon**2 * hb * sigma_prime`` factor explicit and builds
``(A-A0)/hb`` from the radial fundamental theorem of calculus, so no two
nearly equal collar fields are evaluated and subtracted.
"""
from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from .original_signed_density_expression import ScalarExpression, SignedDensityExpressionGraph

SCHEMA = "original_collar_signed_ftc_v1"
ORDERS = (0, 1)
REQUIRED = (
    "Z", "L", "hb", "epsilon", "rho", "ratio", "chi", "phi_r_over_phi",
    "psi_r", "psi_rr", "psi_r_at_4", "sigma", "sigma_prime",
)
_BAD = ("cap", "enclosure", "majorant", "bound_only", "magnitude", "upper_bound")


def _expr(graph: SignedDensityExpressionGraph, value: ScalarExpression, label: str) -> ScalarExpression:
    if not isinstance(value, ScalarExpression) or value.graph is not graph:
        raise TypeError(f"{label}: a same-graph source expression is required")
    if type(value.node) is not int or not 0 <= value.node < len(graph.nodes):
        raise ValueError(f"{label}: source node is out of range")
    operation = graph.nodes[value.node]["operation"].lower()
    if any(token in operation for token in _BAD):
        raise ValueError(f"{label}: a bound or cap is not a source port")
    return value


def _pair(graph, value, label):
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise TypeError(f"{label}: a C0/C1 source pair is required")
    return _expr(graph, value[0], label + ".C0"), _expr(graph, value[1], label + ".C1")


def _bound(graph: SignedDensityExpressionGraph, name: str, role: str) -> ScalarExpression:
    return graph.insert({"operation": "bound_variable", "name": name, "role": role,
                         "local_scope": True, "original_graph_nodes_added_or_replaced": False})


def _sub(graph: SignedDensityExpressionGraph, expression: ScalarExpression,
         variable: ScalarExpression, value: ScalarExpression, name: str) -> ScalarExpression:
    return graph.insert({"operation": "function_substitution", "expression": expression.node,
                         "variable": variable.node, "value": value.node, "name": name,
                         "exact_source_occurrence": True,
                         "original_graph_nodes_added_or_replaced": False})


def _add(a, b):
    return a[0] + b[0], a[1] + b[1]


def _neg(a):
    return -a[0], -a[1]


def _mul(a, b):
    return a[0] * b[0], a[1] * b[0] + a[0] * b[1]


def _scale(a, c):
    return c * a[0], c * a[1]


def _sub_pair(a, b):
    return _add(a, _neg(b))


def _quot(a, b):
    return a[0] / b[0], a[1] / b[0] - a[0] * b[1] / b[0] ** 2


def _integral(graph, integrand, variable, upper, *, name: str, limits: Mapping[str, Any]) -> tuple[ScalarExpression, ScalarExpression]:
    rows = []
    for order, source in enumerate(integrand):
        rows.append(graph.insert({
            "operation": "definite_integral", "integrand": source.node,
            "variable": variable.node, "lower": graph.zero.node, "upper": upper.node,
            "variable_name": graph.nodes[variable.node]["name"],
            "name": f"{name}_C{order}", "ordinary_slow_Z_derivative_order": order,
            "endpoint_independent_of_Z": True, "exact_function_integral": True,
            "numerical_value_not_installed": True,
            "limits": dict(limits), "original_graph_nodes_added_or_replaced": False,
        }))
    return tuple(rows)


def _failure(message: str) -> dict[str, Any]:
    return {"schema": SCHEMA, "passed": False, "source_only": True,
            "source_executed": False, "old_row_equality_claimed": False,
            "numeric_value_available": False, "caller_ports_wired": False,
            "actual_source_nodes_appended": False, "limits_recorded": False,
            "incoming_history_retained": False, "gaps": [message]}


def build(*, source_ports: Mapping[str, Any], limits: Mapping[str, Any],
          incoming_history: Any, source_coordinate: ScalarExpression,
          graph: SignedDensityExpressionGraph) -> dict[str, Any]:
    """Build the exact local source factor from caller-owned C1 ports.

    Every named port is a two-entry ``(C0, C1)`` pair from the same graph. ``source_coordinate``
    must be the actual archived collar-coordinate node; no source node is
    appended here. ``limits`` is the caller's exact admitted domain record and
    is copied into every local integral as provenance; the actual FTC interval
    is explicitly ``0 -> s``. Incoming history is required and never reset.
    """
    try:
        if not isinstance(source_ports, Mapping):
            raise TypeError("source_ports mapping is required")
        if not isinstance(limits, Mapping) or "lower" not in limits or "upper" not in limits:
            raise TypeError("exact lower and upper limits are required")
        if incoming_history is None:
            raise TypeError("incoming_history is required and may not be reset")
        g = graph
        values = {name: _pair(g, source_ports[name], name) for name in REQUIRED}
        s = _bound(g, "s", "lower-collar source coordinate")
        t = _bound(g, "t", "radial FTC integration coordinate")
        u = _bound(g, "u", "outer signed-source integration coordinate")

        coordinate = _expr(g, source_coordinate, "source_coordinate")
        local = {}
        for name in ("rho", "ratio", "chi", "phi_r_over_phi", "psi_r", "psi_rr", "sigma", "sigma_prime"):
            local[name] = tuple(_sub(g, item, coordinate, s, f"{name}_at_s_C{order}")
                                for order, item in enumerate(values[name]))
        values.update(local)
        one = (g.one, g.zero)
        two = (g.constant(2), g.zero)
        sixteen = g.constant(16)
        one_minus_hb = _sub_pair(one, values["hb"])
        pref = _quot(_scale(values["Z"], two[0]), values["L"])
        T = _mul(_mul(_mul(values["rho"], values["rho"]), values["ratio"]), values["psi_r"])
        T0 = _scale(values["psi_r_at_4"], sixteen)
        A = _mul(_mul(pref, one_minus_hb), T)
        A0 = _mul(pref, T0)
        rho_s = _mul(values["hb"], values["rho"])
        ratio_s = _mul(_mul(_mul(_mul(values["hb"], values["rho"]),
                                _sub_pair(values["chi"], one)), values["phi_r_over_phi"]), values["ratio"])

        # Radial FTC density: d(rho^2 ratio Psi_r)/ds divided by hb.
        bracket = _add(_add(_scale(values["psi_r"], two[0]),
                            _mul(values["rho"], values["psi_rr"])),
                       _mul(_mul(_mul(values["rho"], _sub_pair(values["chi"], one)),
                                 values["phi_r_over_phi"]), values["psi_r"]))
        radial_density = _mul(_mul(_mul(values["rho"], values["rho"]), values["ratio"]), bracket)
        radial_at_t = tuple(_sub(g, item, s, t, f"radial_density_at_t_C{order}")
                            for order, item in enumerate(radial_density))
        radial_integral = _integral(g, radial_at_t, t, s, name="radial_FTC_0_to_s", limits=limits)
        normalized = _add(_neg(A0), _mul(_mul(one_minus_hb, pref), radial_integral))
        physical_factor = _mul(_mul(_mul(values["epsilon"], values["epsilon"]), values["hb"]), values["sigma_prime"])
        source_factor = _mul(physical_factor, normalized)
        outer_at_u = tuple(_sub(g, item, s, u, f"signed_source_at_u_C{order}")
                           for order, item in enumerate(source_factor))
        primitive = _integral(g, outer_at_u, u, s, name="signed_A_minus_A0_FTC", limits=limits)
        principal = _mul(_mul(_mul(values["epsilon"], values["epsilon"]), A0), values["sigma"])
        full = _add(principal, primitive)
        return {
            "schema": SCHEMA, "passed": True, "source_only": True, "source_executed": False,
            "ordinary_Z_orders": list(ORDERS), "ordinary_Z1_product_rules": True,
            "caller_ports_wired": True, "archive_authentication_performed_by_constructor": False,
            "source_coordinate": coordinate.node,
            "limits": dict(limits), "limits_recorded": True,
            "incoming_history": incoming_history, "incoming_history_retained": True,
            "old_row_equality_claimed": False, "numeric_value_available": False,
            "formulas": {
                "A": "2*Z/L*(1-hb)*rho^2*ratio*Psi_core_rho",
                "A0": "2*Z/L*T0; T0=16*Psi_core_rho(4,Z)",
                "rho_s": "hb*rho",
                "ratio_s": "hb*rho*(chi-1)*(Phi_core_rho/Phi_core)*ratio",
                "normalized_A_minus_A0": "-A0+(1-hb)*(2*Z/L)*integral_0^s rho^2*ratio*(2*PsiR+rho*PsiRR+rho*(chi-1)*(PhiR/Phi)*PsiR)dt",
                "physical_factor": "epsilon^2*hb*sigma_prime",
            },
            "A": {"C0_C1": [x.node for x in A], "T_C0_C1": [x.node for x in T],
                  "rho_s_C0_C1": [x.node for x in rho_s], "ratio_s_C0_C1": [x.node for x in ratio_s]},
            "A0": {"C0_C1": [x.node for x in A0], "T0_C0_C1": [x.node for x in T0]},
            "radial_FTC": {"density_C0_C1": [x.node for x in radial_density],
                           "integral_0_to_s_C0_C1": [x.node for x in radial_integral],
                           "normalized_A_minus_A0_over_hb_C0_C1": [x.node for x in normalized]},
            "source_factor": {"physical_C0_C1": [x.node for x in source_factor],
                              "factor_C0_C1": [x.node for x in physical_factor],
                              "variable": s.node, "epsilon_squared_once": True,
                              "hb_sigma_prime_explicit": True},
            "signed_primitive": {"integral_C0_C1": [x.node for x in primitive],
                                  "principal_epsilon2_A0_sigma_C0_C1": [x.node for x in principal],
                                  "full_C0_C1": [x.node for x in full], "outer_variable": u.node,
                                  "inner_variable": t.node},
            "graph": {"nodes": g.nodes, "packet": g.packet(),
                      "original_graph_nodes_added_or_replaced": False,
                      "actual_source_nodes_appended": False},
            "gaps": [],
        }
    except (TypeError, ValueError, KeyError) as exc:
        return _failure(str(exc))
