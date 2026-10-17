"""Source proof of the upstream loop/root identities at the two real seams.

The loop archive contains separate native chart expressions.  At the
``second_switch(1) = post_power(0)`` and
``post_power(1) = long_reshape(0)`` scopes, the expressions are the same
physical functions only after the actual R 100 source ports and the actual
radius/global-phase endpoint theorem have been admitted.  This module keeps
that admission explicit and bounded.  It never evaluates an owner, producer,
checker, interval enclosure, or numerical fixture.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_active_viscous_m15 import ActiveEndpointFingerprint
from .original_core_numeric import ROOT


_LOOP_KEY = "prescribed_loop_C5_bindings"
_LEADING_KEY = "prescribed_R100_Rm_C6_bindings"
_CHARTS = ("second_switch", "post_power", "long_reshape")
_SEAMS = {
    "second_power": {"left": ("second_switch", 1), "right": ("post_power", 0)},
    "power_long": {"left": ("post_power", 1), "right": ("long_reshape", 0)},
}
_ROOT_KEYS = ("E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa")
_LOOP_KEYS = ("q", "B_over_Pstar", "B_over_Pstar_phi")
_PORTS = (
    ("fields", "E"),
    ("fields", "V"),
    ("radial", "E_y"),
    ("radial", "V_y"),
    ("generic", "theta_linear"),
    ("generic", "theta_quadratic"),
    ("generic", "p2"),
)

_SOURCE_PINS = {
    "loop": (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
        "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    ),
    "leading": (
        "experiments/root_st073/"
        "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py",
        "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    ),
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


def _row(nodes: Any, value: Any) -> Mapping[str, Any] | None:
    index = _handle(value)
    if index is None or index < 0:
        return None
    try:
        result = nodes[index]
    except (IndexError, KeyError, TypeError):
        return None
    return result if isinstance(result, Mapping) else None


def _unwrap(value: Any, key: str) -> Mapping[str, Any]:
    value = _mapping(value)
    nested = value.get(key)
    return nested if isinstance(nested, Mapping) else value


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    ).hexdigest()


def _repo_root(root: Any) -> Path:
    path = Path(root) if root is not None else ROOT
    if path.is_file():
        path = path.parent
    if (path / "src").is_dir() and (path / "experiments").is_dir():
        return path
    return Path(__file__).resolve().parents[2]


def _pin(repo: Path, name: str, relative: str, expected: str) -> dict[str, Any]:
    path = repo / relative
    try:
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        actual = None
    return {
        "name": name,
        "path": relative,
        "exists": path.is_file(),
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
    }


def _prefix(nodes: Any, archive: Mapping[str, Any]) -> dict[str, Any]:
    archived = archive.get("exact_graph_nodes")
    current = nodes if isinstance(nodes, (list, tuple)) else ()
    passed = isinstance(archived, list) and len(current) >= len(archived)
    if passed:
        try:
            passed = list(current[: len(archived)]) == archived
        except (IndexError, TypeError):
            passed = False
    return {
        "passed": bool(passed),
        "archived_node_count": len(archived) if isinstance(archived, list) else 0,
        "current_node_count": len(current),
        "exact_prefix": bool(passed),
    }


def _chart_alias(name: str) -> str:
    return "post_switch" if name == "post_power" else name


def _scope_premise(value: Any, chart: str, endpoint: int) -> dict[str, Any]:
    scope = _mapping(value)
    actual_chart = scope.get("chart", scope.get("chart_name"))
    actual_endpoint = scope.get("native_query_node", scope.get("native_query", scope.get("endpoint")))
    coordinate = _handle(scope.get("coordinate_node", scope.get("coordinate")))
    radius = _handle(scope.get("radius_node", scope.get("radius")))
    passed = bool(
        scope.get("passed", True) is True
        and actual_chart == chart
        and _handle(actual_endpoint) == endpoint
        and coordinate is not None
        and radius is not None
    )
    return {
        "passed": passed,
        "chart": actual_chart,
        "native_query_node": _handle(actual_endpoint),
        "coordinate_node": coordinate,
        "radius_node": radius,
    }


def _seam_record(packet: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    packet = _mapping(packet)
    for container_key in ("seams", "joins", "source_seam_proofs", "source_premises"):
        container = _mapping(packet.get(container_key))
        if isinstance(container.get(name), Mapping):
            return container[name]
    return _mapping(packet.get(name))


def _section(packet: Mapping[str, Any], *keys: str) -> Mapping[str, Any]:
    packet = _mapping(packet)
    for key in keys:
        value = packet.get(key)
        if isinstance(value, Mapping):
            return value
    return {}


def _germ_receipt(germ: Mapping[str, Any], name: str) -> dict[str, Any]:
    """Validate the independently supplied source germ without broad flags."""
    record = _seam_record(germ, name)
    expected = _SEAMS[name]
    left = _scope_premise(record.get("left_scope"), *expected["left"])
    right = _scope_premise(record.get("right_scope"), *expected["right"])
    contract = _section(record, "source_contract")
    successor = _section(record, "successor_incoming_equalities")
    background = _section(
        record,
        "background_physical_y0_to5",
        "background_physical_germ",
        "actual_background_physical_germ",
    )
    inverse = _section(
        record,
        "unique_inverse_global_phase_source_germ",
        "unique_inverse_global_phase",
    )
    common = _section(record, "common_background_operator_inputs")
    native = _section(record, "issued_native_derivative", "native_derivative_congruence")
    native_body = _mapping(germ.get('native_derivative_certificate'))
    native_sha = _digest(native_body) if native_body else None
    native_verified = bool(native and native_body and native.get('passed') is True and
        native_body.get('passed') is True and native.get('certificate_sha256') == native_sha and
        germ.get('native_derivative_certificate_sha256') == native_sha and
        _mapping(_mapping(native_body.get('seam_premises')).get(name)).get('passed') is True)
    native = dict(native, passed=native_verified, actual_certificate_sha256=native_sha,
                  matching_full_derivative_certificate=native_verified)
    gauge = _mapping(common.get("actual_R110_gauge"))
    recipe_charts = common.get("actual_recipe_charts", ())
    root_rows = inverse.get("actual_root_operator_rows", ())
    root_row_map = {
        row.get("chart"): row
        for row in root_rows
        if isinstance(row, Mapping)
    } if isinstance(root_rows, (list, tuple)) else {}
    source_values_selected = any(
        _mapping(value).get("source_values_selected") is True
        for value in (record, background, inverse, common, gauge)
    )
    numeric_values = any(
        _mapping(value).get(key) is True
        for value in (record, background, inverse, common, gauge)
        for key in ("numeric_caps_used", "numeric_values_selected", "interval_or_cap_values_used")
    )
    background_ok = bool(
        background.get("passed") is True
        and background.get("physical_y_orders") == [0, 1, 2, 3, 4, 5]
        and background.get("ordinary_Z_commutes_with_source_induction") is True
        and background.get("nonzero_original_histories_retained") is True
    )
    recipe_ok = bool(
        common.get("passed") is True
        and common.get("source_germ_not_numeric_values") is True
        and common.get("physical_background_orders") == [0, 1, 2, 3, 4, 5]
        and isinstance(recipe_charts, (list, tuple))
        and all(chart in recipe_charts for chart in (expected["left"][0], expected["right"][0]))
    )
    inverse_ok = bool(
        inverse.get("passed") is True
        and inverse.get("actual_selected_radius_global_phase_seam") in (name, "power_reshape" if name == "power_long" else name)
        and root_row_map.get(expected["left"][0], {}).get("passed") is True
        and root_row_map.get(expected["right"][0], {}).get("passed") is True
    )
    gauge_ok = name != "power_long" or bool(
        gauge.get("passed") is True
        and gauge.get("source_values_selected") is False
        and gauge.get("actual_endpoint_wrapper_rows")
    )
    germ_passed = bool(_mapping(germ).get("passed") is True or record.get("passed") is True)
    passed = bool(
        germ_passed
        and record.get("passed") is True
        and left["passed"] and right["passed"]
        and contract.get("passed") is True
        and successor.get("passed") is True
        and background_ok and recipe_ok and inverse_ok and gauge_ok and native_verified
        and not source_values_selected and not numeric_values
    )
    return {
        "passed": passed,
        "source_only": not source_values_selected and not numeric_values,
        "seam": name,
        "left_scope": left,
        "right_scope": right,
        "source_contract": {"passed": contract.get("passed") is True},
        "successor_incoming_equalities": {"passed": successor.get("passed") is True},
        "background_physical_y0_to5": {
            "passed": background_ok,
            "physical_y_orders": background.get("physical_y_orders"),
            "ordinary_Z_commutes_with_source_induction": background.get("ordinary_Z_commutes_with_source_induction") is True,
            "nonzero_original_histories_retained": background.get("nonzero_original_histories_retained") is True,
        },
        "unique_inverse_global_phase_source_germ": {
            "passed": inverse_ok,
            "selected_seam": inverse.get("actual_selected_radius_global_phase_seam"),
            "root_operator_rows": len(root_rows) if isinstance(root_rows, (list, tuple)) else 0,
        },
        "common_background_operator_inputs": {
            "passed": recipe_ok,
            "actual_recipe_receipt_key": common.get("actual_recipe_receipt_key"),
            "actual_recipe_charts": list(recipe_charts) if isinstance(recipe_charts, (list, tuple)) else [],
            "physical_background_orders": common.get("physical_background_orders"),
            "actual_R110_gauge_passed": gauge.get("passed") is True if gauge else None,
            "shared_F0_nodes": common.get("shared_F0_nodes"),
            "Pstar_node": common.get("Pstar_node"),
        },
        "actual_R110_gauge": {
            "passed": gauge_ok,
            "wrapper_rows": len(gauge.get("actual_endpoint_wrapper_rows", ())) if isinstance(gauge, Mapping) else 0,
        },
        "issued_native_derivative": dict(native) if native else {},
        "source_values_selected": source_values_selected,
        "numeric_caps_used": numeric_values,
    }


def _serial_ring(value: Any) -> dict[str, Any]:
    if not isinstance(value, tuple) or len(value) != 2:
        return {"sha256": None, "rational": None}
    digest, rational = value
    return {
        "sha256": digest,
        "rational": None if rational is None else [rational.numerator, rational.denominator],
    }


class ScopedActiveEndpointFingerprint(ActiveEndpointFingerprint):
    """Endpoint ring with exact, side-scoped R 100 source authorizations.

    A source port is replaced by its common physical-function atom only when
    the call carries the exact endpoint scope that was authorized.  Calls at
    another endpoint, from another chart, or with no scope fall back to the
    ordinary graph fingerprint and cannot inherit this identity.
    """

    def __init__(self, nodes: Any, binding: Mapping[str, Any], authorizations: Mapping[int, Any], scope_tokens: Mapping[int, Any]):
        super().__init__(nodes, binding)
        self.authorizations = dict(authorizations)
        self.scope_tokens = dict(scope_tokens)
        self.alias_hits: set[int] = set()

    def at(self, node: Any, scope: Mapping[int, Any] | None = None):
        index = _handle(node)
        authorization = self.authorizations.get(index)
        if authorization is not None and scope is not None:
            candidates = authorization if isinstance(authorization, (list, tuple)) and authorization and isinstance(authorization[0], (list, tuple)) else (authorization,)
            for coordinate, endpoint_token, value in candidates:
                if scope.get(coordinate) == endpoint_token:
                    self.alias_hits.add(index)
                    return value
        return super().at(node, scope)


def _common_port_rows(nodes: Any, loop: Mapping[str, Any], leading: Mapping[str, Any], name: str) -> dict[str, Any]:
    chart = _mapping(_mapping(loop.get("charts")).get(name))
    leading_name = _chart_alias(name)
    source = _mapping(_mapping(leading.get("charts")).get(leading_name))
    rows: list[dict[str, Any]] = []
    passed = bool(chart and source)
    for group, field in _PORTS:
        values = source.get(group, {}).get(field) if isinstance(source.get(group), Mapping) else None
        root_values = chart.get("roots", {}).get(field) if field in _ROOT_KEYS else None
        if field == "p2":
            root_values = _mapping(chart.get("roots")).get("p2")
        valid = isinstance(values, (list, tuple)) and len(values) >= 3 and all(_row(nodes, x) is not None for x in values[:3])
        # The loop source explicitly binds E, V, and p2 to these leading ports.
        if field in ("E", "V", "p2"):
            valid = valid and isinstance(root_values, (list, tuple)) and list(root_values[:6]) == list(values[:6])
        rows.append({
            "group": group,
            "field": field,
            "leading_chart": leading_name,
            "leading_nodes": list(values[:3]) if isinstance(values, (list, tuple)) else [],
            "loop_nodes": list(root_values[:3]) if isinstance(root_values, (list, tuple)) else [],
            "passed": bool(valid),
        })
        passed = passed and valid
    return {"passed": bool(passed), "chart": name, "leading_chart": leading_name, "rows": rows}


def _root_operator_receipt(
    nodes: Any,
    loop: Mapping[str, Any],
    leading: Mapping[str, Any],
    germ: Mapping[str, Any],
    name: str,
) -> dict[str, Any]:
    """Bind each germ row to the exact loop and R 100 source ports."""
    inverse = _mapping(germ.get("unique_inverse_global_phase_source_germ"))
    rows = inverse.get("actual_root_operator_rows", ())
    charts = _mapping(loop.get("charts")); leading_charts = _mapping(leading.get("charts"))
    expected_names = [_SEAMS[name]["left"][0], _SEAMS[name]["right"][0]]
    selected = []
    failures = []
    for chart_name in expected_names:
        chart = _mapping(charts.get(chart_name))
        leading_name = _chart_alias(chart_name)
        leading_chart = _mapping(leading_charts.get(leading_name))
        candidate = next((row for row in rows if isinstance(row, Mapping) and row.get("chart") == chart_name), None)
        actual_phase_row = [chart_name, chart.get("roots"), chart.get("actual_phase"), chart.get("physical_measure")]
        expected_ports = {
            "E": list(_mapping(leading_chart.get("fields")).get("E", ())[:6]),
            "V": list(_mapping(leading_chart.get("fields")).get("V", ())[:6]),
            "p2": list(_mapping(leading_chart.get("generic")).get("p2", ())[:6]),
        }
        actual_ports = _mapping(candidate.get("actual_leading_E_V_p2_nodes")) if candidate else {}
        actual_p2 = candidate.get("leading_p2_nodes") if candidate else None
        passed = bool(
            candidate is not None
            and candidate.get("passed") is True
            and candidate.get("actual_root_phase_J_row") == actual_phase_row
            and candidate.get("coordinate_node") == _handle(chart.get("coordinate"))
            and candidate.get("radius_node") == _handle(chart.get("radius"))
            and actual_ports == {"E": expected_ports["E"], "V": expected_ports["V"]}
            and (list(actual_p2[:6]) == expected_ports["p2"] if isinstance(actual_p2, (list, tuple)) else False)
        )
        row = {
            "chart": chart_name,
            "leading_chart": leading_name,
            "passed": passed,
            "candidate_present": candidate is not None,
            "coordinate_node": _handle(chart.get("coordinate")),
            "radius_node": _handle(chart.get("radius")),
            "actual_root_phase_J_row_bound": candidate.get("actual_root_phase_J_row") if candidate else None,
            "expected_root_phase_J_row_bound": actual_phase_row,
            "actual_leading_E_V_nodes": actual_ports,
            "actual_leading_p2_nodes": actual_p2,
            "expected_leading_E_V_nodes": {key: value for key, value in expected_ports.items() if key in ("E", "V")},
            "expected_leading_p2_nodes": expected_ports["p2"],
            "source_values_selected": False,
        }
        selected.append(row)
        if not passed:
            failures.append({"chart": chart_name, "reason": "germ root/operator row is not bound to the current chart ports"})
    return {
        "passed": bool(len(selected) == 2 and not failures),
        "seam": name,
        "rows": selected,
        "failures": failures,
        "source_values_selected": False,
    }


def _root_operator_row(loop_chart: Mapping[str, Any], leading_chart: Mapping[str, Any], nodes: Any, quantity: str, order: int) -> bool:
    roots = _mapping(loop_chart.get("roots"))
    if quantity in ("E", "V"):
        return _handle((roots.get(quantity) or [None])[order]) == _handle((_mapping(leading_chart.get("fields")).get(quantity) or [None])[order])
    if quantity == "p2":
        return _handle((roots.get(quantity) or [None])[order]) == _handle((_mapping(leading_chart.get("generic")).get("p2") or [None])[order])
    return _row(nodes, (roots.get(quantity) or [None])[order]) is not None


def _scope_tokens(nodes: Any, charts: Mapping[str, Any], left: str, right: str) -> tuple[dict[int, Any], dict[str, Any], dict[str, dict[int, Any]]]:
    left_chart = _mapping(charts.get(left)); right_chart = _mapping(charts.get(right))
    probe = ActiveEndpointFingerprint(nodes, {"coordinate": _handle(left_chart.get("coordinate")), "loop": charts})
    left_value = probe.at(_handle((left_chart.get("domain") or [None, None])[1]), {})
    right_value = probe.at(_handle((right_chart.get("domain") or [None, None])[0]), {})
    left_scope = {_handle(left_chart.get("coordinate")): left_value}
    right_scope = {_handle(right_chart.get("coordinate")): right_value}
    return (
        {**left_scope, **right_scope},
        {"left": _serial_ring(left_value), "right": _serial_ring(right_value)},
        {"left": left_scope, "right": right_scope},
    )


def _make_fingerprint(
    nodes: Any,
    loop: Mapping[str, Any],
    leading: Mapping[str, Any],
    name: str,
    germ: Mapping[str, Any],
) -> tuple[ScopedActiveEndpointFingerprint, dict[int, Any], dict[str, dict[int, Any]], dict[str, Any]]:
    spec = _SEAMS[name]
    left, right = spec["left"][0], spec["right"][0]
    charts = _mapping(loop.get("charts"))
    scope, serialized_scope, endpoint_scopes = _scope_tokens(nodes, charts, left, right)
    authorizations: dict[int, Any] = {}
    authorization_conflicts: list[dict[str, Any]] = []
    authorization_rows: list[dict[str, Any]] = []
    # All actual R 100 primitives are reduced to the same named physical
    # function atom.  The source germ proves that these are common operators,
    # while this class enforces the two endpoint scopes.
    for quantity_group, quantity in _PORTS:
        for order in range(3):
            values: dict[str, int | None] = {}
            for chart_name in (left, right):
                chart = _mapping(_mapping(leading.get("charts")).get(_chart_alias(chart_name)))
                group = _mapping(chart.get(quantity_group))
                node = _handle((group.get(quantity) or [None] * 3)[order])
                values[chart_name] = node
            token = (_digest(["actual_R100_common_primitive", quantity, order]), None)
            # Preserve exact source constants on both sides.  If either
            # chart's port is an exact rational, do not alias the other side
            # to an opaque atom; the ordinary endpoint ring must prove the
            # constant pair directly.
            constant_pair = any(
                _row(nodes, node) is not None and _row(nodes, node).get("operation") == "exact_rational"
                for node in values.values()
            )
            for chart_name, node in values.items():
                coordinate = _handle(_mapping(charts.get(chart_name)).get("coordinate"))
                expected = scope.get(coordinate)
                if constant_pair:
                    continue
                if node is not None and coordinate is not None and expected is not None:
                    # Exact source constants remain exact rationals.  They
                    # must never be replaced by a common opaque atom.
                    source_row = _row(nodes, node)
                    if source_row and source_row.get("operation") == "exact_rational":
                        continue
                    prior = authorizations.get(node)
                    if prior is not None and not isinstance(prior, list):
                        prior = [prior]
                        authorizations[node] = prior
                    if prior is not None and any(item[2] != token for item in prior):
                        authorization_conflicts.append({
                            "node": node,
                            "prior_sha256": prior[0][2][0],
                            "new_sha256": token[0],
                            "quantity": quantity,
                            "ordinary_Z_order": order,
                        })
                    elif prior is None:
                        authorizations[node] = (coordinate, expected, token)
                    elif (coordinate, expected, token) not in prior:
                        prior.append((coordinate, expected, token))
            authorization_rows.append({
                "group": quantity_group,
                "field": quantity,
                "ordinary_Z_order": order,
                "left_node": values[left],
                "right_node": values[right],
                "canonical_sha256": token[0],
                "scope_restricted": True,
                "exact_constant_pair": constant_pair,
            })
    # Radius and phase are admitted separately by the actual radius/global
    # phase seam theorem; they are not accepted as arbitrary endpoint aliases.
    for role in ("radius", "global_offset", "actual_phase"):
        token = (_digest(["actual_radius_global_phase_endpoint", name, role]), None)
        for chart_name in (left, right):
            chart = _mapping(charts.get(chart_name)); coordinate = _handle(chart.get("coordinate"))
            node = _handle(chart.get(role)); expected = scope.get(coordinate)
            if node is not None and coordinate is not None and expected is not None:
                authorizations[node] = (coordinate, expected, token)
        authorization_rows.append({"group": "geometry", "field": role, "canonical_sha256": token[0], "scope_restricted": True})
    fingerprint = ScopedActiveEndpointFingerprint(
        nodes,
        {"coordinate": _handle(_mapping(charts.get(left)).get("coordinate")), "loop": charts},
        authorizations,
        scope,
    )
    return fingerprint, scope, endpoint_scopes, {
        "passed": bool(_mapping(germ).get("passed") is True and not authorization_conflicts),
        "rows": authorization_rows,
        "authorization_conflicts": authorization_conflicts,
        "scope": serialized_scope,
        "authorized_node_count": len(authorizations),
        "canonical_angle_renaming": {
            "source_prefix": "true_loop_C5_angle_",
            "target_name": ActiveEndpointFingerprint.name("true_loop_C5_angle_" + left),
            "arbitrary_common_actual_phase": True,
        },
    }


def _compare(
    fingerprint: ScopedActiveEndpointFingerprint,
    left_scope: Mapping[int, Any],
    right_scope: Mapping[int, Any],
    label: str,
    left: Any,
    right: Any,
    order: int,
) -> dict[str, Any]:
    left_node, right_node = _handle(left), _handle(right)
    row: dict[str, Any] = {
        "label": label,
        "ordinary_Z_order": order,
        "left_node": left_node,
        "right_node": right_node,
        "source_values_selected": False,
    }
    if left_node is None or right_node is None or _row(fingerprint.nodes, left_node) is None or _row(fingerprint.nodes, right_node) is None:
        row.update(passed=False, reason="missing original source graph handle")
        return row
    try:
        # Keep the two native endpoint scopes separate.  A combined map would
        # silently freeze both native variables and could turn a cross-chart
        # expression into an accidental identity.
        left_value = fingerprint.at(left_node, left_scope)
        right_value = fingerprint.at(right_node, right_scope)
    except Exception as exc:  # bounded source proof keeps the concrete failure
        row.update(passed=False, reason=type(exc).__name__ + ": " + str(exc))
        return row
    row.update(
        passed=left_value[0] == right_value[0],
        left_endpoint_sha256=left_value[0],
        right_endpoint_sha256=right_value[0],
        left_ring=_serial_ring(left_value),
        right_ring=_serial_ring(right_value),
    )
    if not row["passed"]:
        row["unproved_reason"] = "actual source operator or endpoint premise did not normalize the two chart expressions"
    return row


def _root_rows(
    fingerprint: ScopedActiveEndpointFingerprint,
    endpoint_scopes: Mapping[str, Mapping[int, Any]],
    loop: Mapping[str, Any],
    name: str,
) -> list[dict[str, Any]]:
    charts = _mapping(loop.get("charts")); left, right = _SEAMS[name]["left"][0], _SEAMS[name]["right"][0]
    output: list[dict[str, Any]] = []
    for quantity in _ROOT_KEYS:
        left_values = _mapping(charts.get(left)).get("roots", {}).get(quantity)
        right_values = _mapping(charts.get(right)).get("roots", {}).get(quantity)
        for order in range(3):
            row = _compare(fingerprint, endpoint_scopes["left"], endpoint_scopes["right"], f"roots.{quantity}[{order}]", left_values[order], right_values[order], order) \
                if isinstance(left_values, (list, tuple)) and isinstance(right_values, (list, tuple)) and len(left_values) > order and len(right_values) > order \
                else {"label": f"roots.{quantity}[{order}]", "ordinary_Z_order": order, "passed": False, "reason": "root ordinary-Z source row missing", "source_values_selected": False}
            row["quantity"] = quantity
            output.append(row)
    return output


def _loop_rows(
    fingerprint: ScopedActiveEndpointFingerprint,
    endpoint_scopes: Mapping[str, Mapping[int, Any]],
    loop: Mapping[str, Any],
    name: str,
) -> list[dict[str, Any]]:
    charts = _mapping(loop.get("charts")); left, right = _SEAMS[name]["left"][0], _SEAMS[name]["right"][0]
    output: list[dict[str, Any]] = []
    for quantity in _LOOP_KEYS:
        left_values = _mapping(charts.get(left)).get("loop", {}).get(quantity)
        right_values = _mapping(charts.get(right)).get("loop", {}).get(quantity)
        for order in range(3):
            row = _compare(fingerprint, endpoint_scopes["left"], endpoint_scopes["right"], f"loop.{quantity}[{order}]", left_values[order], right_values[order], order) \
                if isinstance(left_values, (list, tuple)) and isinstance(right_values, (list, tuple)) and len(left_values) > order and len(right_values) > order \
                else {"label": f"loop.{quantity}[{order}]", "ordinary_Z_order": order, "passed": False, "reason": "loop ordinary-Z source row missing", "source_values_selected": False}
            row["quantity"] = quantity
            output.append(row)
    return output


def _physical_y_congruence(germ: Mapping[str, Any], fingerprint: ScopedActiveEndpointFingerprint, rows: list[dict[str, Any]], name: str) -> dict[str, Any]:
    background = _mapping(germ.get("background_physical_y0_to5"))
    root_ok = all(row.get("passed") is True for row in rows if row.get("quantity") in ("E", "V", "a", "b", "p1", "p2", "t0", "Delta", "kappa"))
    loop_ok = all(row.get("passed") is True for row in rows if row.get("quantity") in _LOOP_KEYS)
    native = _section(germ, "issued_native_derivative", "native_derivative_congruence")
    native_ok = bool(native) and native.get("passed") is True and native.get("matching_full_derivative_certificate") is True
    passed = bool(
        germ.get("passed") is True
        and background.get("passed") is True
        and root_ok and loop_ok and native_ok
    )
    return {
        "passed": passed,
        "seam": name,
        "ordinary_Z_orders": [0, 1],
        "physical_y_orders": [0, 1, 2, 3, 4, 5],
        "same_physical_background_y0_to5": background.get("passed") is True,
        "same_root_operator_and_native_binding": root_ok,
        "same_loop_operator_q_B_Bphi": loop_ok,
        "issued_native_derivative_congruence_required": True,
        "issued_native_derivative_premise": {"present": bool(native), "passed": native_ok},
        "method": "source germ induction plus scoped native root operator; no pressure alias",
        "native_J_applied_once": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
    }


def build_upstream_loop_identity(
    nodes: Any,
    loop_archive: Mapping[str, Any],
    leading_archive: Mapping[str, Any],
    germ_premises: Mapping[str, Any],
    *,
    root: Any = ROOT,
) -> dict[str, Any]:
    """Build the source-only identity certificate for both upstream seams."""
    loop = _unwrap(loop_archive, _LOOP_KEY)
    leading = _unwrap(leading_archive, _LEADING_KEY)
    repo = _repo_root(root)
    pins = {name: _pin(repo, name, path, expected) for name, (path, expected) in _SOURCE_PINS.items()}
    prefixes = {
        "loop": _prefix(nodes, loop_archive),
        "leading": _prefix(nodes, leading_archive),
    }
    chart_bindings = _mapping(loop.get("charts"))
    leading_bindings = _mapping(leading.get("charts"))
    shape_rows = []
    for name in _CHARTS:
        chart = _mapping(chart_bindings.get(name)); leading_chart = _mapping(leading_bindings.get(_chart_alias(name)))
        shape_rows.append({
            "chart": name,
            "leading_chart": _chart_alias(name),
            "loop_coordinate": _handle(chart.get("coordinate")),
            "loop_radius": _handle(chart.get("radius")),
            "loop_domain": list(chart.get("domain", ())) if isinstance(chart.get("domain"), (list, tuple)) else [],
            "leading_present": bool(leading_chart),
            "passed": bool(
                chart and leading_chart and _handle(chart.get("coordinate")) is not None
                and isinstance(chart.get("domain"), (list, tuple)) and len(chart.get("domain")) == 2
            ),
        })
    required_shapes = all(row["passed"] for row in shape_rows)
    germ_rows = {name: _germ_receipt(germ_premises, name) for name in _SEAMS}
    operator_receipts: dict[str, Any] = {}
    for name in _SEAMS:
        operator = _root_operator_receipt(nodes, loop, leading, _seam_record(germ_premises, name), name)
        operator_receipts[name] = operator
        germ_rows[name]["root_operator_binding"] = operator
        germ_rows[name]["passed"] = bool(germ_rows[name]["passed"] and operator["passed"])
    seams: dict[str, Any] = {}
    for name in _SEAMS:
        left, right = _SEAMS[name]["left"][0], _SEAMS[name]["right"][0]
        ports = {chart: _common_port_rows(nodes, loop, leading, chart) for chart in (left, right)}
        port_passed = all(row["passed"] for row in ports.values())
        fingerprint = None
        scope: dict[int, Any] = {}
        auth_evidence: dict[str, Any] = {"passed": False, "rows": []}
        root_rows: list[dict[str, Any]] = []
        loop_rows: list[dict[str, Any]] = []
        error = None
        try:
            fingerprint, scope, endpoint_scopes, auth_evidence = _make_fingerprint(nodes, loop, leading, name, germ_rows[name])
            root_rows = _root_rows(fingerprint, endpoint_scopes, loop, name)
            loop_rows = _loop_rows(fingerprint, endpoint_scopes, loop, name)
        except Exception as exc:
            error = type(exc).__name__ + ": " + str(exc)
        root_passed = bool(root_rows) and len(root_rows) == len(_ROOT_KEYS) * 3 and all(row.get("passed") is True for row in root_rows)
        loop_passed = bool(loop_rows) and len(loop_rows) == len(_LOOP_KEYS) * 3 and all(row.get("passed") is True for row in loop_rows)
        physical = _physical_y_congruence(germ_rows[name], fingerprint, root_rows + loop_rows, name) if fingerprint else {
            "passed": False,
            "seam": name,
            "reason": error or "scoped source fingerprint unavailable",
            "source_values_selected": False,
            "numeric_caps_used": False,
        }
        seams[name] = {
            "passed": bool(
                germ_rows[name]["passed"] and port_passed and auth_evidence.get("passed") is True
                and root_passed and loop_passed and physical.get("passed") is True
            ),
            "left": left,
            "right": right,
            "left_scope": germ_rows[name]["left_scope"],
            "right_scope": germ_rows[name]["right_scope"],
            "germ_premise": germ_rows[name],
            "actual_root_operator_binding": operator_receipts[name],
            "actual_R100_ports": ports,
            "scoped_fingerprint": auth_evidence,
            "root_rows": root_rows,
            "loop_rows": loop_rows,
            "physical_y_derivative_congruence": physical,
            "alias_hits": sorted(fingerprint.alias_hits) if fingerprint else [],
            "endpoint_scope_restricted": True,
            "canonical_angle_renaming": auth_evidence.get("canonical_angle_renaming", {}),
            "error": error,
        }
    source_pin_passed = all(row.get("passed") is True for row in pins.values())
    prefix_passed = all(row.get("passed") is True for row in prefixes.values())
    passed = bool(source_pin_passed and prefix_passed and required_shapes and all(row.get("passed") is True for row in seams.values()))
    return {
        "passed": passed,
        "method": "scoped actual R100 primitive source ports plus exact bounded loop/root endpoint fingerprint",
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "ancestor_replayed": False,
        "producer_or_checker_executed": False,
        "source_pins": pins,
        "source_hashes": pins,
        "graph_prefixes": prefixes,
        "chart_shape": {"passed": required_shapes, "rows": shape_rows},
        "seams": seams,
        "source_seam_proofs": seams,
        "actual_root_operator_bindings": operator_receipts,
        "checks": {
            "source_pins": {"passed": source_pin_passed, "rows": pins},
            "graph_prefixes": prefixes,
            "chart_shape": {"passed": required_shapes, "rows": shape_rows},
            "germ_premises": {name: row["passed"] for name, row in germ_rows.items()},
            "root_operator_binding": {name: row["passed"] for name, row in operator_receipts.items()},
            "actual_R100_ports": {name: all(item["passed"] for item in seams[name]["actual_R100_ports"].values()) for name in seams},
            "all_root_rows": {name: all(item.get("passed") is True for item in seams[name]["root_rows"]) for name in seams},
            "all_loop_rows": {name: all(item.get("passed") is True for item in seams[name]["loop_rows"]) for name in seams},
            "physical_y_derivative_congruence": {name: seams[name]["physical_y_derivative_congruence"].get("passed") is True for name in seams},
        },
        "limits": {
            "endpoint_scope_only": True,
            "common_actual_phase_is_arbitrary_function": True,
            "no_pressure_alias": True,
            "no_numeric_function_values": True,
            "no_cap_as_function_value": True,
            "field_seam_equality_is_source_operator_congruence": True,
        },
    }


__all__ = ["ScopedActiveEndpointFingerprint", "build_upstream_loop_identity"]
