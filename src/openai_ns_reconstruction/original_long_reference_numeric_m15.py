"""Contract the saved long-reshape/reference viscous endpoint seam.

This module is deliberately a packet adapter.  It hydrates no accepted
producer and evaluates no original signed endpoint.  The only new positive
budget is the directed union of the saved long q=0 pressure packet and the
entire saved postrestore sector 155 aggregate.  The source identity is
recorded separately from that error-only budget.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import json
from typing import Any, Mapping

from .original_active_viscous_m15 import ActiveEndpointFingerprint
from .original_core_numeric import PREFIX, ROOT
from .original_micro_viscous_m15 import contract_trace_chain


SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PREVIOUS_REDUCED_COUNT = 198
NEW_REDUCED_COUNT = 197
UNCHANGED_COUNT = 196
LONG_PREFIX_COUNT = 106
LONG_VISCOUS_INDEX = 105
REFERENCE_VISCOUS_INDEX = 123
RESTORATION_VISCOUS_INDEX = 139
POSTRESTORE_VISCOUS_INDEX = 155

LONG_REPORT = PREFIX + "compliant_current_original_long_m15_numeric.json"
OLD_LONG_TERM = "long_finite_N_viscous_error"
OLD_EXTERNAL_TERM = "reference_restoration_postrestore_finite_N_viscous_external_endpoints"
NEW_EXTERNAL_TERM = "long_reference_restore_finite_N_viscous_external_endpoints"

LONG_Q0 = {"kind": "exact_decimal_long_fraction", "value": "0"}
LONG_Q1 = {"kind": "exact_decimal_long_fraction", "value": "1"}


def _canonical(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _digest(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a mapping")
    return value


def _copy_mapping(value: Any, label: str) -> dict[str, Any]:
    return copy.deepcopy(dict(_mapping(value, label)))


def _find(rows: Any, *, label: str, **criteria: Any) -> Mapping[str, Any]:
    for row in rows:
        if isinstance(row, Mapping) and all(row.get(k) == v for k, v in criteria.items()):
            return row
    raise ValueError(f"Saved {label} not found: {criteria}")


def _bindings(value: Any, key: str) -> Mapping[str, Any]:
    """Accept either a raw archive or its already-unwrapped bindings."""
    if isinstance(value, Mapping) and key in value:
        return _mapping(value[key], key)
    return _mapping(value, key)


def _pin_source(owner: Any, relative: str) -> str:
    raw = (ROOT / relative).read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    known = _mapping(getattr(owner, "code_hashes", {}), "owner.code_hashes").get(relative)
    if known != actual:
        raise ValueError(f"Changed accepted source: {relative}")
    return actual


def _ast_statement_set(function: ast.AST) -> set[str]:
    return {ast.dump(node) for node in ast.walk(function)}


def _bind_long_pressure_rows(owner: Any) -> dict[str, Any]:
    """Pin the actual long pressure wrapper's owner and derivative dispatch."""
    source = "src/openai_ns_reconstruction/original_long_pressure_rows.py"
    source_sha = _pin_source(owner, source)
    tree = ast.parse((ROOT / source).read_bytes())
    cls = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "OriginalLongPressureRows"
    )
    init = next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == "__init__")
    statements = [
        "view._covers = self.long._covers",
        "OriginalMacroPressureRows.__init__(self, micro=view)",
        "self.chart, self.derivative_chart = 'long_reshape', 'reshape'",
        "self.loop = json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart']['reshape']",
    ]
    dumps = _ast_statement_set(init)
    for statement in statements:
        if ast.dump(ast.parse(statement).body[0]) not in dumps:
            raise ValueError(f"Changed long pressure owner statement: {statement}")
    return {
        "passed": True,
        "source": source,
        "source_sha256": source_sha,
        "class": "OriginalLongPressureRows",
        "method": "OriginalLongPressureRows.__init__",
        "exact_source_statements": statements,
        "chart": "long_reshape",
        "derivative_chart": "reshape",
        "view_covers_delegate": "self.long._covers",
        "loop_bound_index": "reshape",
        "actual_long_source_owner": True,
        "old_post_power_source_not_reinterpreted": True,
    }


def _long_native_guard(owner: Any) -> dict[str, Any]:
    """Require the new owner's per-chart physical native guard."""
    guards = _mapping(getattr(owner, "native_guards", None), "owner.native_guards")
    charts = _mapping(guards.get("charts"), "owner.native_guards.charts")
    chart = _mapping(charts.get("long_reshape"), "owner.native_guards.charts.long_reshape")
    if guards.get("passed") is not True or chart.get("passed") is not True:
        raise ValueError("Passed long_reshape physical native guard required")
    return {
        "passed": True,
        "chart": "long_reshape",
        "guard": copy.deepcopy(dict(chart)),
        "source_guard_packet": copy.deepcopy(dict(guards)),
        "physical_native_guard_required": True,
        "ordinary_Z_orders": [0, 1],
        "physical_RQy_factor_once": True,
        "physical_Uy_Pstar_factor_once": True,
    }


def _graph_envelope_binding(owner: Any, wrapper: Mapping[str, Any], guard: Mapping[str, Any]) -> dict[str, Any]:
    """Reuse the accepted linear-Qy/source AST certificate at long native q=0."""
    parent = _mapping(owner.parent, "owner.parent")
    prior = _copy_mapping(parent.get("reference0_graph_envelope_binding"), "parent reference0 graph envelope")
    if prior.get("passed") is not True:
        raise ValueError("Accepted reference0 graph-envelope certificate required")
    for key in (
        "exact_numeric_Qy_coefficients",
        "ordinary_Z_product_rule",
        "actual_archived_pressure_formula_certificate",
        "actual_archived_radial_formula_certificate",
    ):
        if key not in prior:
            raise ValueError(f"Accepted universal Qy certificate missing {key}")

    loop = _bindings(owner.l, "prescribed_loop_C5_bindings")
    modified = _bindings(owner.m, "modified_Q4_pressure_C3_bindings")
    charts = _mapping(loop.get("charts"), "loop charts")
    chart = _mapping(charts.get("long_reshape"), "long_reshape loop chart")
    mchart = _mapping(_mapping(modified.get("charts"), "modified charts").get("long_reshape"), "long_reshape modified chart")
    coordinate = chart.get("coordinate")
    endpoint = chart.get("domain", [None, None])[0]
    if coordinate is None or endpoint is None:
        raise ValueError("Long native q=0 endpoint scope required")

    fingerprint = ActiveEndpointFingerprint(owner.nodes, dict(coordinate=coordinate, loop=charts))
    scope = {coordinate: fingerprint.at(endpoint, {})}
    graph_rows = []
    for order in (0, 1):
        qy_node = _mapping(mchart.get("physical"), "long modified physical").get("Q_y", [])[order]
        rqy_node = _mapping(mchart.get("finite_boundary"), "long finite boundary").get("R_Q_y", [])[order]
        qy_hash = fingerprint.at(qy_node, scope)[0]
        rqy_hash = fingerprint.at(rqy_node, scope)[0]
        graph_rows.append(
            {
                "ordinary_Z_order": order,
                "Qy_function_node": qy_node,
                "RQy_function_node": rqy_node,
                "Qy_at_long0_sha256": qy_hash,
                "RQy_at_long0_sha256": rqy_hash,
                "endpoint_scope": {"coordinate_node": coordinate, "endpoint_node": endpoint},
            }
        )

    inherited_keys = (
        "exact_numeric_Qy_coefficients",
        "exact_coefficient_geometry_statements",
        "ordinary_Z_product_rule",
        "actual_archived_pressure_formula_certificate",
        "actual_archived_radial_formula_certificate",
    )
    result = {
        "passed": True,
        "method": "source_bound_linear_Qy_envelope_at_actual_long_native0",
        "inherited_reference0_envelope_sha256": _digest(prior),
        **{key: copy.deepcopy(prior[key]) for key in inherited_keys},
        "endpoint_scope": {"coordinate_node": coordinate, "endpoint_node": endpoint},
        "graph_function_rows": graph_rows,
        "long_pressure_rows_owner": copy.deepcopy(dict(wrapper)),
        "long_native_guard": copy.deepcopy(dict(guard)),
        "generic_nodes_reused_only_under_explicit_scope": True,
        "generic_history_bounds_inherited_not_recertified": True,
        "source_function_values_or_caps_not_selected": True,
        "source_function_not_identified_with_post_power1": True,
        "actual_long_q0_scope_only": True,
    }
    return result


def _long_pressure_report(owner: Any) -> Mapping[str, Any]:
    value = getattr(owner, "long", None)
    if isinstance(value, Mapping):
        report = value
    elif isinstance(getattr(value, "report", None), Mapping):
        report = value.report
    else:
        raise ValueError("Saved long_m15_numeric report required")
    if not isinstance(report.get("pressure_row_samples"), list):
        raise ValueError("Saved long pressure rows required")
    return report


def _own_long_inlet(owner: Any, Z: str, selected_integer: Mapping[str, Any], envelope: Mapping[str, Any]):
    report = _long_pressure_report(owner)
    row = _find(
        report["pressure_row_samples"],
        label="long q=0 pressure row",
        Z=Z,
        coordinate=LONG_Q0,
    )
    if row.get("source_family") != owner.family:
        raise ValueError("Long q=0 pressure source family changed")
    if row.get("selected_integer") != selected_integer:
        raise ValueError("Long q=0 selected N changed")
    cover = _mapping(row.get("source_cover"), "long q=0 source cover")
    if cover.get("generic_derivative_chart") != "reshape":
        raise ValueError("Long q=0 must use generic derivative chart reshape")
    histories = cover.get("cumulative_history_charts")
    if not isinstance(histories, list) or not histories or histories[-1] != "long_reshape":
        raise ValueError("Long history prefix must end at long_reshape")
    contract = _mapping(row.get("contract"), "long q=0 pressure contract")
    if contract.get("original_Pstar_once") is not True:
        raise ValueError("Long q=0 physical Pstar factor must be retained once")
    endpoint = _mapping(row.get("R_Q_y_endpoint"), "long q=0 RQy endpoint")
    errors = _mapping(row.get("physical_pressure_rows"), "long q=0 physical rows").get("Q_y")
    errors = _mapping(errors, "long q=0 Qy rows").get("finite_N_absolute_error_log_bounds")
    endpoint_errors = endpoint.get("finite_N_coefficient_error_log_bounds")
    if errors != endpoint_errors or not isinstance(errors, list) or len(errors) != 2:
        raise ValueError("Long q=0 own Qy C0/Z1 rows required")
    radius = endpoint.get("finite_N_log_physical_factor")
    if radius != _mapping(row.get("geometry"), "long q=0 geometry").get("log_physical_radius"):
        raise ValueError("Long q=0 physical R factor must be sourced once")

    weighted = [
        None if error.get("exact_zero") else owner.decode(error["log_absolute_upper"]) + owner.decode(radius)
        for error in errors
    ]
    physical_rows = _mapping(row.get("physical_pressure_rows"), "long q=0 physical rows")
    own_rows = {name: copy.deepcopy(physical_rows[name].get("finite_N_absolute_error_log_bounds")) for name in ("U", "m", "U_y")}
    slow_fast = copy.deepcopy(row.get("finite_N_axial_radial_error_components"))
    if not isinstance(slow_fast, (dict, list)):
        raise ValueError("Long q=0 slow/fast error components required")
    witness = {
        "source_report": LONG_REPORT,
        "source_packet_sha256": _digest(row),
        "coordinate": copy.deepcopy(row["coordinate"]),
        "selected_integer": copy.deepcopy(selected_integer),
        "source_cover": copy.deepcopy(dict(cover)),
        "actual_own_Qy_error_rows": copy.deepcopy(errors),
        "actual_own_log_radius": copy.deepcopy(radius),
        "actual_linear_error_input_rows": own_rows,
        "own_U_m_Uy_error_rows": copy.deepcopy(own_rows),
        "actual_slow_fast_error_components": slow_fast,
        "own_slow_fast_components": copy.deepcopy(slow_fast),
        "radius_weighted_C0_C1_log_bounds": [owner.record(value) if value is not None else None for value in weighted],
        "actual_long0_graph_envelope_binding": copy.deepcopy(envelope),
        "generic_derivative_chart": "reshape",
        "long_q0_source_packet": True,
        "physical_R_once": True,
        "physical_Pstar_already_in_error": True,
        "extra_epsilon": False,
        "source_function_values_evaluated": False,
        "post_power1_packet_not_reinterpreted": True,
        "reference0_packet_not_used_as_long_q0_cap": True,
    }
    return weighted, witness


def _sample(report: Mapping[str, Any], Z: str, chart: str, coordinate: Mapping[str, Any]) -> Mapping[str, Any]:
    return _find(report.get("samples", []), label=f"{chart} sample", Z=Z, chart=chart, coordinate=coordinate)


def _selected_N_evidence(owner: Any, Z: str, long_q1: Mapping[str, Any]) -> dict[str, Any]:
    reference = _mapping(getattr(owner, "reference", None), "owner.reference")
    restore = _mapping(getattr(owner, "restore", None), "owner.restore")
    ref0 = _sample(reference, Z, "reference", {"kind": "exact_decimal_reference", "value": "0"})
    post1 = _sample(restore, Z, "postrestore", {"kind": "exact_decimal_postrestore", "value": "1"})
    selected = {
        "long_reshape_native1": copy.deepcopy(long_q1.get("actual_selected_integer")),
        "reference_native0": copy.deepcopy(ref0.get("actual_selected_integer")),
        "postrestore_native1": copy.deepcopy(post1.get("actual_selected_integer")),
    }
    if not (selected["long_reshape_native1"] == selected["reference_native0"] == selected["postrestore_native1"]):
        raise ValueError("Saved long/reference/postrestore selected N must agree")
    selected["same_selected_N"] = True
    return selected


def _chart_bindings(owner: Any) -> tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]:
    loop = _bindings(owner.l, "prescribed_loop_C5_bindings")
    modified = _bindings(owner.m, "modified_Q4_pressure_C3_bindings")
    outer = _bindings(owner.outer, "prescribed_R100_Rm_C6_bindings")
    return loop, modified, outer


def _function_ports(owner: Any, charts: tuple[str, ...] = ("long_reshape", "postrestore")) -> list[dict[str, Any]]:
    loop, modified, outer = _chart_bindings(owner)
    loop_charts = _mapping(loop.get("charts"), "loop charts")
    modified_charts = _mapping(modified.get("charts"), "modified charts")
    outer_charts = _mapping(outer.get("charts"), "outer charts")
    ports: list[dict[str, Any]] = []

    for name in charts:
        chart = _mapping(loop_charts.get(name), f"loop chart {name}")
        mchart = _mapping(modified_charts.get(name), f"modified chart {name}")
        outer_name = "post_restore" if name == "postrestore" else name
        bg = _mapping(outer_charts.get(outer_name), f"outer chart {outer_name}")
        radius = chart.get("radius")
        physical_measure = chart.get("physical_measure")
        if radius is None or physical_measure is None:
            raise ValueError(f"Physical radius/measure missing for {name}")

        visited: set[int] = set()

        def z_free(node: int) -> None:
            if node in visited:
                return
            visited.add(node)
            row = owner.nodes[node]
            operation = row["operation"]
            if operation == "bound_variable":
                if node != chart["coordinate"]:
                    raise ValueError(f"Unexpected Z-dependent physical node in {name}")
            elif operation in ("original_source_parameter", "current_original_source_parameter", "exact_positive_integer_expression"):
                if not row.get("Z_independent"):
                    raise ValueError(f"Z-dependent physical parameter in {name}")
            elif operation in ("sum", "product"):
                for argument in row["arguments"]:
                    z_free(argument)
            elif operation in ("negative", "analytic_unary"):
                z_free(row["argument"])
            elif operation == "positive_quotient":
                z_free(row["numerator"])
                z_free(row["denominator"])
            elif operation not in ("exact_rational", "mathematical_pi"):
                raise ValueError(f"Unissued physical-radius operation {operation}")

        z_free(radius)
        z_free(physical_measure)
        finite_boundary = _mapping(mchart.get("finite_boundary"), f"{name} finite boundary")
        rqy = finite_boundary.get("R_Q_y")
        if not isinstance(rqy, (list, tuple)) or len(rqy) < 2:
            raise ValueError(f"Ordinary C0/Z1 RQy rows missing for {name}")
        for node in rqy[:2]:
            row = owner.nodes[node]
            if row.get("operation") != "product" or row.get("arguments", []).count(radius) != 1:
                raise ValueError(f"One explicit physical R factor required for {name}")

        physical = _mapping(mchart.get("physical"), f"{name} physical")
        radial = _mapping(bg.get("radial"), f"{name} leading radial")
        fields = _mapping(bg.get("fields"), f"{name} leading fields")
        for side, endpoint in zip(("lower", "upper"), chart.get("domain", ())):
            ports.append(
                {
                    "chart": name,
                    "side": side,
                    "coordinate_node": chart["coordinate"],
                    "endpoint_node": endpoint,
                    "endpoint_scope": {"coordinate_node": chart["coordinate"], "endpoint_node": endpoint},
                    "physical_radius_node": radius,
                    "native_Jacobian_node": physical_measure,
                    "modified_Qy_nodes": list(physical.get("Q_y", [])[:2]),
                    "modified_RQy_nodes": list(rqy[:2]),
                    "leading_normalized_V_nodes": list(fields.get("V", [])[:3]),
                    "leading_normalized_mean_nodes": list(fields.get("m", [])[:3]),
                    "leading_normalized_Vy_nodes": list(radial.get("V_y", [])[:2]),
                    "definition": "D=R*(Qy_modified-Qy_leading)",
                    "radius_and_measure_Z_free_dependency_nodes": sorted(visited),
                    "physical_R_explicit_once_in_both_Z_rows": True,
                    "ordinary_Z_orders": [0, 1],
                    "signed_coefficient": 1 if side == "lower" else -1,
                    "graph_scope_is_explicit": True,
                }
            )
    return ports


def _seam_bundle(owner: Any, seam: Mapping[str, Any]) -> dict[str, Any]:
    if seam.get("passed") is not True:
        raise ValueError("Actual longreshape1/reference0 RQy seam must pass")
    parent = _mapping(owner.parent, "owner.parent")
    accepted = _mapping(parent.get("source_seam_proofs"), "accepted source seam proofs")
    if accepted.get("passed") is not True:
        raise ValueError("Accepted reference/restoration/post seams required")
    return {
        "passed": True,
        "long_reference": copy.deepcopy(dict(seam)),
        "accepted_reference_restoration": copy.deepcopy(accepted.get("reference_restoration_modified_RQy_rows", [])),
        "accepted_restoration_postrestore": copy.deepcopy(accepted.get("restoration_postrestore_certificate", {})),
        "accepted_source_seam_proofs": copy.deepcopy(dict(accepted)),
        "internal_proofs": {
            "longreshape_native1_reference_native0": True,
            "reference_native1_restoration_native0": True,
            "restoration_native1_postrestore_native0": True,
        },
        "function": "D(longreshape native0)-D(postrestore native1=Rm)",
        "delta_definition": "R*(Qy_modified-Qy_leading)",
        "ordinary_Z_orders": [0, 1],
        "physical_R_once": True,
        "leading_subtraction_retained": True,
        "global_fast_phase_not_reset": True,
        "native_Jacobians_not_equated": True,
        "source_functions_not_error_caps": True,
    }


def _trace_chain() -> dict[str, Any]:
    pieces = [
        {"name": "long_reshape", "lower": "long_native0", "upper": "long_native1"},
        {"name": "reference", "lower": "reference_native0", "upper": "reference_native1"},
        {"name": "restoration", "lower": "restoration_native0", "upper": "restoration_native1"},
        {"name": "postrestore", "lower": "postrestore_native0", "upper": "Rm"},
    ]
    aliases = {
        "long_native1": "reference_native0",
        "reference_native1": "restoration_native0",
        "restoration_native1": "postrestore_native0",
    }
    result = contract_trace_chain(pieces, aliases)
    result["external_endpoints"] = {"lower": "long native0", "upper": "postrestore native1 = Rm"}
    result["cancelled_interfaces_explicit"] = [
        "long native1 = reference native0",
        "reference native1 = restoration native0",
        "restoration native1 = postrestore native0",
    ]
    result["function"] = "D(longreshape native0)-D(postrestore native1=Rm)"
    return result


def _build_source_binding(owner: Any, seam: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    wrapper = _bind_long_pressure_rows(owner)
    native_guard = _long_native_guard(owner)
    graph = _graph_envelope_binding(owner, wrapper, native_guard)
    parent_source = _mapping(_mapping(owner.parent, "owner.parent").get("source_binding"), "parent source binding")
    inherited_keys = (
        "passed",
        "producer_rows",
        "function",
        "physical_radius_power",
        "ordinary_Z_orders",
        "positive_sum_domination",
        "ordinary_Z1_product_rule_retained",
        "physical_Pstar_already_in_Qy_error",
        "slow_B_y_over_N_and_fast_B_phi_both_retained",
        "postrestore_entire_original_155_cap_retained",
        "inherited_pressure_bound_validity_not_independently_recertified",
    )
    source = {
        key: copy.deepcopy(parent_source[key])
        for key in inherited_keys
        if key in parent_source
    }
    source.update(
        {
            "passed": True,
            "function": "D=R*(Qy_modified-Qy_leading); chart contribution=D(lower)-D(upper)",
            "physical_radius_power": 1,
            "ordinary_Z_orders": [0, 1],
            "ordinary_Z1_product_rule_retained": True,
            "physical_Pstar_already_in_Qy_error": True,
            "slow_B_y_over_N_and_fast_B_phi_both_retained": True,
            "separate_inlet_uses_own_longreshape0_pressure_packet": True,
            "longreshape0_source_packet_is_not_post_power1": True,
            "post_power1_not_reinterpreted": True,
            "no_unproved_whole_chart_alias": True,
            "postrestore_entire_original_155_cap_retained": True,
            "inherited_pressure_bound_validity_not_independently_recertified": True,
            "long_pressure_rows_owner": wrapper,
            "long_native_physical_guard": native_guard,
            "actual_long0_graph_envelope_binding": graph,
            "generic_nodes_reused_only_under_explicit_scope": True,
            "generic_history_bounds_inherited_not_recertified": True,
            "new_long_reference_seam_required": True,
            "accepted_reference_restoration_postrestore_seams_required": True,
            "no_assumed_long0_post_power_seam": True,
            "source_function_values_or_caps_not_selected": True,
        }
    )
    rows = list(source.get("producer_rows", []))
    rows.append(
        {
            "source": wrapper["source"],
            "source_sha256": wrapper["source_sha256"],
            "method": wrapper["method"],
            "exact_source_statements": wrapper["exact_source_statements"],
        }
    )
    source["producer_rows"] = rows
    seams = _seam_bundle(owner, seam)
    return source, seams


def _evaluate_sample(owner: Any, source: Mapping[str, Any], seams: Mapping[str, Any], Z: str) -> dict[str, Any]:
    parent = _mapping(owner.parent, "owner.parent")
    parent_sample = _find(parent.get("samples", []), label="active parent sample", Z=Z)
    raw = parent_sample.get("raw_m15_C0_C1_terms")
    previous = parent_sample.get("correlated_reduced_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != RAW_COUNT:
        raise ValueError("Raw 213 m15 state must remain verbatim")
    if not isinstance(previous, list) or len(previous) != PREVIOUS_REDUCED_COUNT:
        raise ValueError("Previous correlated 198 state required")
    long_report = _long_pressure_report(owner)
    long_q1 = _sample(long_report, Z, "long_reshape", LONG_Q1)
    if len(long_q1.get("m15_C0_C1_terms", [])) != LONG_PREFIX_COUNT:
        raise ValueError("Saved long q=1 106-term prefix required")
    if raw[:LONG_PREFIX_COUNT] != long_q1["m15_C0_C1_terms"]:
        raise ValueError("Raw prefix must match long q=1 with exact coordinate kind")
    if raw[LONG_VISCOUS_INDEX].get("name") != OLD_LONG_TERM:
        raise ValueError("Long finite-N viscous sector identity required")
    old_replacements = [term for term in previous if term.get("name") == OLD_EXTERNAL_TERM]
    if len(old_replacements) != 1:
        raise ValueError("Previous reference/restoration/postrestore replacement required")
    if sum(term.get("name") == OLD_LONG_TERM for term in previous) != 1:
        raise ValueError("Previous long finite-N viscous term required")
    selected = _selected_N_evidence(owner, Z, long_q1)
    inlet, inlet_witness = _own_long_inlet(owner, Z, selected["long_reshape_native1"], source["actual_long0_graph_envelope_binding"])

    outgoing = raw[POSTRESTORE_VISCOUS_INDEX]
    if outgoing.get("name") != "postrestore_finite_N_viscous_error":
        raise ValueError("Postrestore finite-N viscous sector identity required")
    outgoing_logs = outgoing.get("log_scale_terms")
    if not isinstance(outgoing_logs, list) or len(outgoing_logs) != 1:
        raise ValueError("Entire original postrestore sector 155 aggregate required")
    if outgoing_logs[0].get("positive_error_bound_atom_only") is not True or outgoing_logs[0].get("source_value_selected") is not False:
        raise ValueError("Postrestore sector 155 must remain an error-only positive aggregate")

    positive_inputs = [value for value in inlet if value is not None]
    positive_inputs.append(owner.decode(outgoing_logs[0]))
    bound = owner.positive_union(positive_inputs)

    replacement = copy.deepcopy(old_replacements[0])
    replacement["name"] = NEW_EXTERNAL_TERM
    replacement["log_scale_terms"] = [owner.record(bound)]
    replacement["actual_function_binding"] = {
        "error_enclosure_only": True,
        "function": "D(longreshape native0)-D(postrestore native1=Rm)",
        "delta_definition": "R*(Qy_modified-Qy_leading)",
        "ordinary_Z1": True,
        "physical_R_once": True,
        "external_endpoints_not_zeroed": True,
        "signed_centers_evaluated": False,
        "no_error_cap_subtraction": True,
        "no_interval_midpoint_or_overlap_closure": True,
        "directed_positive_union_inlet_order_then_full155": True,
        "independent_own_longreshape0_budget": True,
        "full_original_postrestore_155_budget_retained": True,
        "internal_source_seam_proved": True,
        "new_long_reference_seam": copy.deepcopy(seams["long_reference"]),
        "accepted_reference_restoration_postrestore_seams": True,
        "generic_history_bounds_inherited_not_recertified": True,
    }

    reduced: list[dict[str, Any]] = []
    removed: list[str] = []
    for term in previous:
        name = term.get("name")
        if name == OLD_LONG_TERM:
            removed.append(name)
            continue
        if name == OLD_EXTERNAL_TERM:
            removed.append(name)
            reduced.append(copy.deepcopy(replacement))
        else:
            reduced.append(copy.deepcopy(term))
    if set(removed) != {OLD_LONG_TERM, OLD_EXTERNAL_TERM} or len(reduced) != NEW_REDUCED_COUNT:
        raise ValueError("Exact new 197-term view required")
    unchanged_old = [term for term in previous if term.get("name") not in {OLD_LONG_TERM, OLD_EXTERNAL_TERM}]
    unchanged_new = [term for term in reduced if term.get("name") != NEW_EXTERNAL_TERM]
    if len(unchanged_old) != UNCHANGED_COUNT or unchanged_new != unchanged_old:
        raise ValueError("Other 196 correlated sectors must remain verbatim")
    if sum(term.get("name") == NEW_EXTERNAL_TERM for term in reduced) != 1:
        raise ValueError("Exactly one new long/reference replacement required")

    ports = _function_ports(owner)
    old_provenance = [copy.deepcopy(raw[index]) for index in (LONG_VISCOUS_INDEX, REFERENCE_VISCOUS_INDEX, RESTORATION_VISCOUS_INDEX, POSTRESTORE_VISCOUS_INDEX)]
    target = copy.deepcopy(parent_sample.get("independent_target"))
    if target is None:
        target = copy.deepcopy(parent.get("independent_target"))
    if target is None:
        raise ValueError("Independent target must remain available")

    return {
        "Z": Z,
        "source_family": owner.family,
        "raw_m15_C0_C1_terms": copy.deepcopy(raw),
        "correlated_reduced_m15_C0_C1_terms": reduced,
        "original_raw_terms": RAW_COUNT,
        "previous_reduced_terms": PREVIOUS_REDUCED_COUNT,
        "correlated_reduced_terms": NEW_REDUCED_COUNT,
        "unchanged_other_sectors": UNCHANGED_COUNT,
        "signed_trace_chain": _trace_chain(),
        "source_seam_proofs": copy.deepcopy(dict(seams)),
        "actual_signed_function_ports": ports,
        "external_endpoint_budget_witness": {
            "own_longreshape0": inlet_witness,
            "postrestore_full_original_155_aggregate": copy.deepcopy(outgoing),
            "positive_union": owner.record(bound),
            "positive_union_inputs": copy.deepcopy(inlet_witness["radius_weighted_C0_C1_log_bounds"])
            + [copy.deepcopy(outgoing_logs[0])],
            "source_binding": copy.deepcopy(dict(source)),
            "postrestore_upper_separate_budget_claimed": False,
            "long_finite_N_cap_replaced": True,
            "no_cap_subtraction_value_midpoint_or_overlap": True,
        },
        "old_budget_terms_retained_as_provenance": old_provenance,
        "previous_external_replacement_retained_as_provenance": copy.deepcopy(old_replacements[0]),
        "replacement_outer_endpoints": replacement,
        "selected_N_evidence": selected,
        "independent_target": target,
        "largest_remaining_individual_log_caps": owner.largest_caps(reduced, target),
        "contract": {
            "raw_213_verbatim": True,
            "other_196_verbatim": True,
            "independent_target_unchanged": True,
            "reserved_added_this_stage": 0,
            "ordinary_Z1_and_physical_R_once": True,
            "leading_subtraction_retained": True,
            "long_and_postrestore_external_source_functions_retained": True,
            "long_q0_uses_own_pressure_packet": True,
            "post_power1_not_reinterpreted": True,
            "no_assumed_long0_post_power_seam": True,
            "all_other_volume_width_numerical_errors_retained": True,
            "ancestors_not_reexecuted": True,
            "source_bound_signed_trace_contraction": True,
            "error_cap_as_value": False,
        },
        "limits": {
            "long_reference_internal_trace_resolved": True,
            "inherited_pressure_bounds_independently_recertified": False,
            "external_signed_values_resolved": False,
            "independent_numeric_Rp_seam": False,
            "global_five_moment_common_N_closure": False,
            "exact_heat_and_energy": False,
            "admissible_stress_flat_remainder": False,
            "temporal_n_recursion_and_summation": False,
            "mean_oscillatory_correction": False,
            "independent_corrected_Cartesian_NS": False,
        },
    }


def build_long_reference_numeric(owner: Any, seam: Mapping[str, Any]) -> dict[str, Any]:
    """Return the saved long/reference m15 contraction packet.

    ``owner`` is the new source owner assembled by the caller.  It must expose
    only saved reports, graph archives and the new per-chart native guard; this
    function intentionally has no producer or ancestor constructor calls.
    """
    if not isinstance(seam, Mapping):
        raise ValueError("A longreshape1/reference0 seam mapping is required")
    source, seams = _build_source_binding(owner, seam)
    samples = [_evaluate_sample(owner, source, seams, Z) for Z in SAVED_Z]
    return {
        "source_binding": source,
        "samples": samples,
        "source_seam_proofs": seams,
        "contract": {
            "saved_packets_only": True,
            "accepted_ancestor_constructors_not_run": True,
            "cap_is_not_signed_function_value": True,
            "long_q0_graph_scope_explicit": True,
            "post_power1_not_reinterpreted": True,
        },
        "limits": {
            "numeric_Rp_closure": False,
            "global_all_five_closure": False,
            "full_goal_complete": False,
        },
    }


__all__ = ["build_long_reference_numeric"]
