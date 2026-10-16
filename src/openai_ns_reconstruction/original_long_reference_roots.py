"""Source-bound root ring for the actual long/reference native seam.

This module proves only the ordinary-Z root identities needed by the
long-reshape native1 to reference native0 interface.  It deliberately does
not call the older reference/restoration certificate: that certificate has a
different chart scope and its endpoint assumptions do not apply here.

The active ``q`` branch is admitted only after the two actual endpoint root
rings and the inherited positive canonical ``eta`` premise have passed.  No
interval endpoint, numerical field value, or error cap is used as a source
function value.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Any, Mapping

from .original_reference_root_Z import _RootZFingerprint, _parameter_receipt


_CHARTS = (
    ("long_reshape", 1, "native1"),
    ("reference", 0, "native0"),
)
_ROOT_KEYS = ("a", "b", "t0", "Delta")
_EXPECTED = {
    "a": Fraction(4, 5),
    "b": Fraction(0),
    "t0": Fraction(0),
    "Delta": Fraction(-6, 5),
}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _handle(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _row(nodes: Any, index: Any) -> Mapping[str, Any] | None:
    i = _handle(index)
    if i is None or i < 0:
        return None
    try:
        value = nodes[i]
    except (IndexError, KeyError, TypeError):
        return None
    return value if isinstance(value, Mapping) else None


def _unwrap_loop(value: Mapping[str, Any]) -> Mapping[str, Any]:
    nested = value.get("prescribed_loop_C5_bindings")
    if isinstance(nested, Mapping):
        return nested
    # An archive's prescribed_loop_C5_bindings itself also has a metadata
    # field named ``bindings``.  Only unwrap that field when the outer object
    # is genuinely a wrapper and the nested value carries the chart map.
    nested = value.get("bindings")
    if "charts" not in value and isinstance(nested, Mapping) and "charts" in nested:
        return nested
    return value


def _unwrap_germ(value: Mapping[str, Any]) -> Mapping[str, Any]:
    for key in (
        "certificate",
        "long_reference_germ_certificate",
        "source_audit_certificate",
        "m15_source_audit",
    ):
        nested = value.get(key)
        if isinstance(nested, Mapping):
            return nested
    return value


def _nested_mapping(value: Mapping[str, Any], keys: tuple[str, ...]) -> Mapping[str, Any]:
    for key in keys:
        candidate = value.get(key)
        if isinstance(candidate, Mapping):
            return candidate
    return {}


def _germ_receipt(germ_packet: Mapping[str, Any]) -> dict[str, Any]:
    germ = _unwrap_germ(_mapping(germ_packet))
    source_only = bool(
        germ.get("source_only") is True
        and germ.get("source_values_selected", False) is False
        and germ.get("numeric_caps_used", False) is False
    )
    evidence = _mapping(germ.get("source_evidence"))
    evidence_groups = evidence.get("groups")
    expected_groups = (
        "long_reference_outer_formula",
        "flat_sigma_prime",
        "bridge_fields_and_radial",
        "duhamel_histories",
        "flat_sigma_source",
        "generic_root_recipe",
        "poisson_loop_recipe",
    )
    evidence_ok = bool(
        evidence.get("passed") is True
        and isinstance(evidence_groups, Mapping)
        and all(_mapping(evidence_groups.get(name)).get("passed") is True for name in expected_groups)
    )
    interface = _mapping(germ.get("common_source_identity_rows"))
    interface_scope = _mapping(interface.get("scope"))
    scope_ok = bool(
        interface.get("passed") is True
        and interface.get("source_values_selected") is False
        and isinstance(interface.get("rows"), list)
        and interface.get("rows")
        and _mapping(interface_scope.get("long_reshape")).get("endpoint") == 1
        and _mapping(interface_scope.get("reference")).get("endpoint") == 0
        and isinstance(_mapping(interface_scope.get("long_reshape")).get("coordinate"), int)
        and isinstance(_mapping(interface_scope.get("reference")).get("coordinate"), int)
    )
    physical = _mapping(germ.get("physical_two_jet_certificate"))
    physical_E_rows = physical.get("E_y_graph_rows")
    physical_V_rows = physical.get("V_y_graph_rows")
    actual_fields_ok = bool(
        isinstance(physical_E_rows, list)
        and physical_E_rows
        and all(_mapping(row).get("passed") is True for row in physical_E_rows)
        and isinstance(physical_V_rows, list)
        and physical_V_rows
        and all(_mapping(row).get("passed") is True for row in physical_V_rows)
    )
    physical_ok = bool(
        physical.get("passed") is True
        and physical.get("source_values_selected", False) is False
        and physical.get("interval_or_cap_values_used", True) is False
        and physical.get("ordinary_Z_orders") == [0, 1, 2]
        and physical.get("E_y") == "E/10"
        and physical.get("E_yy") == "E/100"
        and physical.get("V_y") == "0 at long_reshape native1 and reference native0"
        and physical.get("V_yy") == "0 at the seam"
        and actual_fields_ok
        and _mapping(physical.get("common_Duhamel_second_jets")).get("passed") is True
    )
    return {
        "passed": bool(germ.get("passed") is True and source_only and evidence_ok and scope_ok and physical_ok),
        "source_only": source_only,
        "source_evidence": {
            "passed": evidence_ok,
            "groups": list(expected_groups),
        },
        "common_source_interface": {
            "passed": scope_ok,
            "scope": {
                "long_reshape": dict(_mapping(interface_scope.get("long_reshape")),),
                "reference": dict(_mapping(interface_scope.get("reference")),),
            },
            "rows": len(interface.get("rows", [])) if isinstance(interface.get("rows"), list) else 0,
        },
        "physical_two_jet": {
            "passed": physical_ok,
            "actual_E_y_rows": len(physical_E_rows) if isinstance(physical_E_rows, list) else 0,
            "actual_V_y_rows": len(physical_V_rows) if isinstance(physical_V_rows, list) else 0,
            "actual_fields_checked": actual_fields_ok,
            "endpoint_scope_checked": scope_ok,
            "whole_long_chart_identity_not_claimed": True,
        },
        "numeric_values_or_caps_used": not source_only,
    }


def _eta_receipt(nodes: Any, loop: Mapping[str, Any]) -> dict[str, Any]:
    eta = _handle(loop.get("eta"))
    row = _row(nodes, eta)
    proof = row.get("canonical_parameter_proof") if row else None
    canonical_eta = proof.get("original_canonical_eta_log") if isinstance(proof, Mapping) else None
    passed = bool(
        row
        and row.get("operation") == "current_original_source_parameter"
        and row.get("name") == "eta"
        and row.get("source_attribute") == "eta_log"
        and row.get("Z_independent") is True
        and row.get("quantity_is_exact_original_parameter_not_range_endpoint") is True
        and isinstance(proof, Mapping)
        and isinstance(canonical_eta, Mapping)
        and isinstance(canonical_eta.get("lower_exact_mpf_tuple"), list)
        and isinstance(canonical_eta.get("upper_exact_mpf_tuple"), list)
    )
    return {
        "passed": passed,
        "node": eta,
        "operation": row.get("operation") if row else None,
        "source_attribute": row.get("source_attribute") if row else None,
        "Z_independent": row.get("Z_independent") if row else None,
        "canonical_positive_source": passed,
        "numeric_eta_value_selected": False,
    }


def _port(chart: str, binding: Mapping[str, Any], endpoint_index: int, label: str) -> dict[str, Any]:
    domain = binding.get("domain")
    endpoint = domain[endpoint_index] if isinstance(domain, (list, tuple)) and len(domain) > endpoint_index else None
    coordinate = _handle(binding.get("coordinate"))
    return {
        "chart": chart,
        "native_endpoint": label,
        "coordinate_node": coordinate,
        "endpoint_node": _handle(endpoint),
        "endpoint_scope": {"coordinate_node": coordinate, "endpoint_node": _handle(endpoint)},
        "root_nodes": {
            key: {
                "Z0": _handle(_mapping(binding.get("roots")).get(key, [None, None])[0])
                if isinstance(_mapping(binding.get("roots")).get(key), (list, tuple)) else None,
                "Z1": _handle(_mapping(binding.get("roots")).get(key, [None, None])[1])
                if isinstance(_mapping(binding.get("roots")).get(key), (list, tuple)) else None,
            }
            for key in _ROOT_KEYS
        },
        "q_nodes": {
            "q0": _handle(_mapping(_mapping(binding.get("loop")).get("q")).get("q0"))
            if isinstance(_mapping(binding.get("loop")).get("q"), Mapping) else (
                _handle(_mapping(binding.get("loop")).get("q", [None, None])[0])
                if isinstance(_mapping(binding.get("loop")).get("q"), (list, tuple)) else None
            ),
            "q1": _handle(_mapping(_mapping(binding.get("loop")).get("q")).get("q1"))
            if isinstance(_mapping(binding.get("loop")).get("q"), Mapping) else (
                _handle(_mapping(binding.get("loop")).get("q", [None, None])[1])
                if isinstance(_mapping(binding.get("loop")).get("q"), (list, tuple)) else None
            ),
        },
    }


def _q_nodes(binding: Mapping[str, Any]) -> tuple[int | None, int | None]:
    q = _mapping(binding.get("loop")).get("q")
    if isinstance(q, Mapping):
        return _handle(q.get("q0")), _handle(q.get("q1"))
    if isinstance(q, (list, tuple)):
        return (_handle(q[0]) if len(q) > 0 else None,
                _handle(q[1]) if len(q) > 1 else None)
    return None, None


def _root_row(
    chart: str,
    port: Mapping[str, Any],
    quantity: str,
    order: int,
    node: int | None,
    value: tuple[str, Fraction | None] | None,
    *,
    reason: str | None = None,
) -> dict[str, Any]:
    expected = _EXPECTED[quantity] if order == 0 else Fraction(0)
    exact_value = None if value is None or value[1] is None else [value[1].numerator, value[1].denominator]
    return {
        "chart": chart,
        "native_endpoint": port.get("native_endpoint"),
        "quantity": quantity,
        "ordinary_Z_order": order,
        "original_node": node,
        "expected_value": [expected.numerator, expected.denominator],
        "exact_value": exact_value,
        "ring_digest": value[0] if value is not None else None,
        "exact_ring_row": bool(value is not None and value[1] == expected),
        "source_values_selected": False,
        "source_port": {
            "coordinate_node": port.get("coordinate_node"),
            "endpoint_node": port.get("endpoint_node"),
        },
        **({"unproved_reason": reason} if reason else {}),
    }


def _q_row(
    chart: str,
    port: Mapping[str, Any],
    q0: int | None,
    q1: int | None,
    value0: tuple[str, Fraction | None] | None,
    value1: tuple[str, Fraction | None] | None,
    *,
    reason: str | None = None,
) -> dict[str, Any]:
    exact = bool(value1 is not None and value1[1] == Fraction(0))
    return {
        "chart": chart,
        "native_endpoint": port.get("native_endpoint"),
        "quantity": "q",
        "ordinary_Z_order": 1,
        "original_node": q1,
        "q_value_node": q0,
        "q_value_digest": value0[0] if value0 is not None else None,
        "q_value_exact": None if value0 is None or value0[1] is None else [value0[1].numerator, value0[1].denominator],
        "expected_value": [0, 1],
        "exact_value": None if value1 is None or value1[1] is None else [value1[1].numerator, value1[1].denominator],
        "ring_digest": value1[0] if value1 is not None else None,
        "exact_ring_row": exact,
        "source_values_selected": False,
        "active_branch_evaluated_after_root_ring": reason is None,
        "source_port": {
            "coordinate_node": port.get("coordinate_node"),
            "endpoint_node": port.get("endpoint_node"),
        },
        **({"unproved_reason": reason} if reason else {}),
    }


def build_long_reference_roots(nodes: Any, l: Mapping[str, Any], germ: Mapping[str, Any]) -> dict[str, Any]:
    """Build the exact source certificate for long native1/reference native0 roots.

    The function is intentionally total with respect to missing or failed
    premises: it returns ``passed=False`` and preserves unproved rows instead
    of relabeling an unresolved graph expression as a root identity.
    """
    loop = _unwrap_loop(_mapping(l))
    germ_packet = _unwrap_germ(_mapping(germ))
    charts = _mapping(loop.get("charts"))
    ports: dict[str, dict[str, Any]] = {}
    for name, endpoint_index, label in _CHARTS:
        ports[name] = _port(name, _mapping(charts.get(name)), endpoint_index, label)

    parameter_receipt: dict[str, Any]
    try:
        parameter_receipt = _parameter_receipt(nodes, loop)
    except (ArithmeticError, KeyError, TypeError, ValueError, IndexError) as exc:
        parameter_receipt = {"passed": False, "error": type(exc).__name__ + ": " + str(exc)}
    eta_receipt = _eta_receipt(nodes, loop)
    germ_receipt = _germ_receipt(germ_packet)

    shape_ok = bool(
        len(charts) >= 2
        and all(
            ports[name]["coordinate_node"] is not None
            and ports[name]["endpoint_node"] is not None
            and all(
                isinstance(row, Mapping)
                and row.get("Z0") is not None
                and row.get("Z1") is not None
                for row in ports[name]["root_nodes"].values()
            )
            and ports[name]["q_nodes"]["q0"] is not None
            and ports[name]["q_nodes"]["q1"] is not None
            for name, _, _ in _CHARTS
        )
    )

    q_nodes: dict[int, tuple[int, int]] = {}
    delta_eta: dict[int, int] = {}
    for name, _, _ in _CHARTS:
        chart = _mapping(charts.get(name))
        roots = _mapping(chart.get("roots"))
        delta_values = roots.get("Delta")
        delta = _handle(delta_values[0]) if isinstance(delta_values, (list, tuple)) and delta_values else None
        eta = _handle(loop.get("eta"))
        q0, q1 = ports[name]["q_nodes"]["q0"], ports[name]["q_nodes"]["q1"]
        if delta is not None and eta is not None and q0 is not None and q1 is not None:
            delta_eta[delta] = eta
            q_nodes[q0] = (delta, eta)
            q_nodes[q1] = (delta, eta)

    fingerprint: _RootZFingerprint | None = None
    scopes: dict[str, dict[int, tuple[str, Fraction | None]]] = {}
    root_rows: list[dict[str, Any]] = []
    ring_passed = False
    if shape_ok and len(q_nodes) == 4 and len(delta_eta) == 2:
        try:
            base = ports["long_reshape"]
            fingerprint = _RootZFingerprint(
                nodes,
                {"coordinate": base["coordinate_node"], "loop": charts},
                q_nodes=q_nodes,
                delta_eta=delta_eta,
            )
            shared_scope: dict[int, tuple[str, Fraction | None]] = {}
            for name, _, _ in _CHARTS:
                port = ports[name]
                endpoint_value = fingerprint.at(port["endpoint_node"], {})
                shared_scope[port["coordinate_node"]] = endpoint_value
            for name, _, _ in _CHARTS:
                scopes[name] = dict(shared_scope)
                chart = _mapping(charts.get(name))
                roots = _mapping(chart.get("roots"))
                for quantity in _ROOT_KEYS:
                    values = roots.get(quantity)
                    for order in (0, 1):
                        node = _handle(values[order]) if isinstance(values, (list, tuple)) and len(values) > order else None
                        value = None
                        reason = None
                        if node is None:
                            reason = "missing actual root source node"
                        else:
                            try:
                                value = fingerprint.at(node, shared_scope)
                            except (ArithmeticError, KeyError, TypeError, ValueError, IndexError) as exc:
                                reason = type(exc).__name__ + ": " + str(exc)
                        root_rows.append(_root_row(name, ports[name], quantity, order, node, value, reason=reason))
            ring_passed = len(root_rows) == 16 and all(row["exact_ring_row"] for row in root_rows)
        except (ArithmeticError, KeyError, TypeError, ValueError, IndexError) as exc:
            root_rows = []
            for name, _, _ in _CHARTS:
                for quantity in _ROOT_KEYS:
                    for order in (0, 1):
                        values = ports[name]["root_nodes"].get(quantity, {})
                        root_rows.append(_root_row(
                            name,
                            ports[name],
                            quantity,
                            order,
                            values.get("Z" + str(order)),
                            None,
                            reason=type(exc).__name__ + ": " + str(exc),
                        ))
            ring_passed = False
    else:
        reason = "actual long/reference root source scope incomplete"
        for name, _, _ in _CHARTS:
            for quantity in _ROOT_KEYS:
                for order in (0, 1):
                    values = ports[name]["root_nodes"].get(quantity, {})
                    root_rows.append(_root_row(
                        name,
                        ports[name],
                        quantity,
                        order,
                        values.get("Z" + str(order)),
                        None,
                        reason=reason,
                    ))

    active_gap_passed = bool(
        ring_passed
        and eta_receipt.get("passed") is True
        and germ_receipt.get("passed") is True
    )
    q_rows: list[dict[str, Any]] = []
    q_active_proof: dict[str, Any]
    if active_gap_passed and fingerprint is not None:
        for name, _, _ in _CHARTS:
            q0, q1 = ports[name]["q_nodes"]["q0"], ports[name]["q_nodes"]["q1"]
            try:
                value0 = fingerprint.at(q0, scopes[name])
                value1 = fingerprint.at(q1, scopes[name])
                q_rows.append(_q_row(name, ports[name], q0, q1, value0, value1))
            except (ArithmeticError, KeyError, TypeError, ValueError, IndexError) as exc:
                q_rows.append(_q_row(
                    name, ports[name], q0, q1, None, None,
                    reason=type(exc).__name__ + ": " + str(exc),
                ))
        q_active_proof = {
            "passed": len(q_rows) == 2 and all(row["exact_ring_row"] for row in q_rows),
            "q_nodes": sorted(q_nodes),
            "q_active_rows": sorted(fingerprint.q_active_rows),
            "sigma_one_rows": sorted(fingerprint.sigma_one_rows),
            "sigma_flat_rows": sorted(fingerprint.sigma_flat_rows),
            "active_branch_evaluated_after_root_ring": True,
            "root_ring_passed_before_q": ring_passed,
            "inherited_positive_eta_before_q": eta_receipt.get("passed") is True,
            "exact_gap_identity": "Delta=-6/5<eta; 1-Delta/eta=1+6/(5*eta)>1",
            "sigma_identity": "sigma(1-Delta/eta)=1 and all admitted q flat derivatives are zero",
            "numeric_values_or_caps_used": False,
        }
    else:
        reason = "q active branch withheld until root ring and inherited positivity pass"
        for name, _, _ in _CHARTS:
            q0, q1 = ports[name]["q_nodes"]["q0"], ports[name]["q_nodes"]["q1"]
            q_rows.append(_q_row(name, ports[name], q0, q1, None, None, reason=reason))
        q_active_proof = {
            "passed": False,
            "q_nodes": sorted(q_nodes),
            "q_active_rows": [],
            "sigma_one_rows": [],
            "sigma_flat_rows": [],
            "active_branch_evaluated_after_root_ring": False,
            "root_ring_passed_before_q": ring_passed,
            "inherited_positive_eta_before_q": eta_receipt.get("passed") is True,
            "numeric_values_or_caps_used": False,
            "unproved_reason": reason,
        }

    rows = root_rows + q_rows
    passed = bool(
        len(rows) == 18
        and ring_passed
        and q_active_proof.get("passed") is True
        and parameter_receipt.get("passed") is True
        and eta_receipt.get("passed") is True
        and germ_receipt.get("passed") is True
    )
    return {
        "passed": passed,
        "method": "actual_long_reshape_native1_to_reference_native0_root_fingerprint_with_scoped_q_active_source",
        "scope": {
            name: {
                "chart": name,
                "native_endpoint": label,
                "coordinate_node": ports[name]["coordinate_node"],
                "endpoint_node": ports[name]["endpoint_node"],
                "endpoint_scope": ports[name]["endpoint_scope"],
            }
            for name, _, label in _CHARTS
        },
        "source_only": True,
        "numeric_values_selected": False,
        "numeric_caps_used": False,
        "germ_premise_receipt": germ_receipt,
        "parameter_premise_receipt": parameter_receipt,
        "canonical_eta_positivity_receipt": eta_receipt,
        "active_gap_and_cut_proof": {
            "passed": active_gap_passed,
            "Delta": [-6, 5],
            "eta_positive_canonical": eta_receipt.get("passed") is True,
            "germ_source_premises_passed": germ_receipt.get("passed") is True,
            "identity": "Delta=-6/5<eta; cut=1-Delta/eta=1+6/(5*eta)>1",
            "sigma": "1",
            "q_active_branch_requires_root_ring_and_positivity": True,
        },
        "ports": ports,
        "rows": rows,
        "row_count": len(rows),
        "root_rows": 16,
        "q_rows": 2,
        "root_ring_passed": ring_passed,
        "q_active_proof": q_active_proof,
        "limits": {
            "whole_long_chart_Ey_identity": False,
            "long_q1_flat_jet_only": True,
            "numeric_endpoint_value_or_cap_as_source": False,
            "approved_ancestor_constructor_run": False,
        },
    }


__all__ = ["build_long_reference_roots"]
