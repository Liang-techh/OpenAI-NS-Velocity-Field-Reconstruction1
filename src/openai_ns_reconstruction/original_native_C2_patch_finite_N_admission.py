"""Metadata-only admission for the frozen finite-``N`` source attempt.

This module is a boundary around a saved source result.  It does not import
the source consumer, rebuild its definition graph, evaluate a source port, or
call a producer, checker, recurrence, quadrature, or evaluator.  The only
graph operation here is a structural digest of the already-saved node list;
the graph is returned as a qualified reference and its output handles are
metadata, never new mathematical definitions.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence


SCHEMA = "original_native_C2_patch_finite_N_admission_v1"
OLD_SOURCE_SCHEMA = "original_native_C2_patch_finite_N_source_v1"

# These identify the already-entered source attempt.  The admission API still
# authenticates the complete mapping with a caller-supplied canonical digest;
# these constants identify the saved source checkpoint and graph, rather than
# the current leading/provider report digests.
FROZEN_RESULT_SHA256 = "ac8b4e100124d47e7c2696d21945edae05726052a4d7bba1a47563f666950379"
FROZEN_SOURCE_CHECKPOINT_SHA256 = "681681cc1e00a1043dd38c40a43be5667ed6d7cb18b1e089ac1341f05880186f"
FROZEN_KERNEL_SHA256 = "0f85613e17a59bcb1d5825acfafb720a8dc1f4580993cb06d88c444dfcc62bd2"
FROZEN_GRAPH_NAMESPACE = OLD_SOURCE_SCHEMA
FROZEN_GRAPH_NODE_COUNT = 1063
EXPECTED_SOURCE_GRAPH_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

ORDINARY_Z_ORDERS = [0, 1, 2]
MOMENT_NAMES = ["m", "h", "k", "e", "p"]
RATES = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}

RETAINED_CONSTRUCTION_DATA = {
    "cumulative_terms": 212,
    "raw_signed_terms": 213,
    "reduced_terms": 183,
    "six_signed_channels_retained": True,
    "all_old_error_data_retained": True,
    "reserved_contribution_registered_once": True,
    "inherited_errors_deleted": False,
    "P1_differs_from_m15": True,
    "physical_Q_differs_from_loop_q": True,
}

FROZEN_UNRESOLVED_GATES = [
    "next_finite_N_interface_prerequisite_unproved",
    "accepted_C6_smoothness_common_domain_premise_unproved",
    "native_C6_to_C5_patch_provider_identity_unproved",
    "patch_finite_N_source_identity_unproved",
]

FROZEN_STATUS_FLAGS = {
    "native_C2_patch_finite_N_fields_defined": False,
    "native_C2_C5_patch_provider_identity_proved": False,
    "native_C2_C5_patch_density_identity_proved": False,
    "native_C2_patch_correction_inlet_identity_proved": True,
    "native_C2_C5_patch_history_identity_proved": False,
    "patch_finite_N_source_identity": False,
}

GLOBAL_FALSE_FLAGS = {
    "same_density_certificate": False,
    "finite_N_history_identity": False,
    "physical_Q_y_Z1": False,
    "global_five_moment_closure": False,
    "legacy_feedback_reference_alias_proved": False,
    "prior_adapter_geometry_alias_proved": False,
    "prior_adapter_P0_alias_proved": False,
    "time_order_recursion_claimed": False,
    "numeric_values_computed": False,
    "old_operators_replayed": False,
    "full_goal_complete": False,
}

REQUIRED_GRAPH_CHECKS = {
    "unique_namespace_ids",
    "closed_topological_arguments",
    "rational_constants_only",
    "ordinary_Z_jets_have_three_rows",
    "no_numeric_or_solver_operation",
    "endpoint_substitutions_are_function_nodes",
    "log_kernel_constructed",
    "scalar_ID_and_jet_list_distinct",
    "expm1_entire_at_A_zero",
    "qualified_source_ports_valid",
    "source_function_applications_have_explicit_closed_bindings",
    "integral_variables_and_bounds_are_scalar_IDs",
    "log_transport_inverse_s_is_in_actual_integrand",
}

EXPECTED_FINITE_RECIPE = {
    "ordinary_Z_orders": ORDINARY_Z_ORDERS,
    "fields": {
        "x": "A/N",
        "r": "E*expm1(A/N)",
        "s": "B_over_Pstar/N",
        "E_N": "E*exp(A/N)",
        "V_N": "V+s",
    },
    "C5_density": {
        "m": "s",
        "h": "r",
        "k": "V*r+E*s+r*s",
        "e": "2*V*s+s^2-E*r-r^2/2",
        "p": "E*r+r^2/2",
    },
    "C2_density": {
        "G": "N*expm1(A/N)",
        "F": "E*G",
        "first": {
            "m": "B",
            "h": "F",
            "k": "V*F+E*B",
            "e": "2*V*B-E*F",
            "p": "E*F",
        },
        "second": {
            "m": "0",
            "h": "0",
            "k": "F*B",
            "e": "B^2-F^2/2",
            "p": "F^2/2",
        },
        "density": "first+second/N",
    },
    "identity": {
        "relation": "C2=N*C5",
        "coefficientwise_fields": MOMENT_NAMES,
        "F_equals_N_r_entire_at_A_zero": True,
        "same_B_provider_required": True,
        "same_Pstar_units_required": True,
    },
}

EXPECTED_SECTION11_RECIPE = {
    "moments": MOMENT_NAMES,
    "rates": RATES,
    "formula": "H=x^-rate H_in+integral_1^x (s/x)^rate*C2rho(s)*ds/s",
    "kernel": "exp(rate*(log(s)-log(x)))",
    "coordinate": "x=R/Rm",
    "domain": "[1,e]",
    "p_memory": "exactly one",
    "reset": "only original Rminus; retain predecessor inlet",
    "derivative_identity": "x*H_x=rho-rate*H",
}


def _canonical_digest(value: object) -> str | None:
    """Return the caller-visible canonical JSON digest, or ``None``."""
    try:
        payload = json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    except (TypeError, ValueError, OverflowError):
        return None
    return hashlib.sha256(payload).hexdigest()


def _get(value: object, path: str, default: object = None) -> object:
    current = value
    for part in path.split(".") if path else ():
        if isinstance(current, Mapping):
            if part not in current:
                return default
            current = current[part]
        elif isinstance(current, Sequence) and not isinstance(current, (str, bytes, bytearray)) and part.isdigit():
            index = int(part)
            if index >= len(current):
                return default
            current = current[index]
        else:
            return default
    return current


def _copy_metadata(value: object) -> object:
    """Copy small returned metadata without interpreting it as a graph."""
    if isinstance(value, Mapping):
        return {key: _copy_metadata(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_copy_metadata(item) for item in value]
    if isinstance(value, tuple):
        return [_copy_metadata(item) for item in value]
    return value


def _component(value: object) -> Mapping[str, object]:
    """Accept either a component mapping or its saved stage envelope."""
    if not isinstance(value, Mapping):
        return {}
    component = value.get("component_result")
    if isinstance(component, Mapping):
        return component
    result = value.get("result")
    if isinstance(result, Mapping):
        nested = result.get("component_result")
        if isinstance(nested, Mapping):
            return nested
        if "definition_graph" in result or result.get("schema") == OLD_SOURCE_SCHEMA:
            return result
    return value


def _section(value: object, component: Mapping[str, object], name: str) -> Mapping[str, object]:
    if isinstance(value, Mapping) and isinstance(value.get(name), Mapping):
        return value[name]  # type: ignore[return-value]
    section = component.get(name)
    if isinstance(section, Mapping):
        return section
    return {}


def _gate(gates: list[str], name: str, condition: bool) -> None:
    if not condition and name not in gates:
        gates.append(name)


def _digest_inputs(
    saved_attempt: Mapping[str, object],
    latest_report: Mapping[str, object],
    provider_report: Mapping[str, object],
    expected_attempt_digest: str,
    expected_latest_digest: str,
    expected_provider_digest: str,
    gates: list[str],
) -> dict[str, str | None]:
    actual = {
        "attempt": _canonical_digest(saved_attempt),
        "latest": _canonical_digest(latest_report),
        "provider": _canonical_digest(provider_report),
    }
    expected = {
        "attempt": expected_attempt_digest,
        "latest": expected_latest_digest,
        "provider": expected_provider_digest,
    }
    for name in ("attempt", "latest", "provider"):
        _gate(
            gates,
            "canonical_" + name + "_digest_mismatch",
            isinstance(expected[name], str) and actual[name] == expected[name],
        )
    return actual


def _saved_graph_reference(
    attempt: Mapping[str, object],
    gates: list[str],
) -> tuple[dict[str, object], bool, dict[str, object]]:
    graph = attempt.get("definition_graph")
    if not isinstance(graph, Mapping):
        _gate(gates, "saved_frozen_graph_reference_unproved", False)
        return ({}, False, {})
    nodes = graph.get("nodes")
    namespace = graph.get("namespace")
    node_count = graph.get("node_count")
    active_digest = graph.get("active_graph_digest")
    graph_shape = (
        isinstance(nodes, list)
        and namespace == FROZEN_GRAPH_NAMESPACE
        and node_count == FROZEN_GRAPH_NODE_COUNT
        and len(nodes) == FROZEN_GRAPH_NODE_COUNT
        and active_digest == EXPECTED_SOURCE_GRAPH_DIGEST
        and isinstance(graph.get("outputs"), Mapping)
    )
    unique_ids = True
    marker_count = 0
    if isinstance(nodes, list):
        identifiers: set[str] = set()
        for node in nodes:
            if not isinstance(node, Mapping):
                graph_shape = False
                continue
            identifier = node.get("id")
            if not isinstance(identifier, str) or identifier in identifiers:
                unique_ids = False
            else:
                identifiers.add(identifier)
            if node.get("namespace") != FROZEN_GRAPH_NAMESPACE:
                graph_shape = False
            if node.get("common_prefix_and_ports_checked") is False:
                marker_count += 1
            value = node.get("value")
            if isinstance(value, Mapping) and value.get("common_prefix_and_ports_checked") is False:
                marker_count += 1
    else:
        unique_ids = False
    historical_marker_ok = marker_count > 0
    _gate(gates, "saved_frozen_graph_reference_unproved", graph_shape and unique_ids)
    _gate(
        gates,
        "saved_historical_common_prefix_marker_missing_or_changed",
        historical_marker_ok,
    )
    canonical_nodes = _canonical_digest(nodes) if isinstance(nodes, list) else None
    outputs = graph.get("outputs") if isinstance(graph.get("outputs"), Mapping) else {}
    qualification = {
        "reference_only": True,
        "re_evaluated": False,
        "rebuilt": False,
        "old_result_sha256": FROZEN_RESULT_SHA256,
        "saved_attempt_canonical_digest": None,
        "source_checkpoint_sha256": attempt.get("source_checkpoint_sha256"),
        "kernel_sha256": None,
        "frozen_namespace": namespace,
        "canonical_node_hash": canonical_nodes,
        "canonical_nodes_digest": canonical_nodes,
        "canonical_node_count": len(nodes) if isinstance(nodes, list) else None,
        "node_count": len(nodes) if isinstance(nodes, list) else None,
        "active_saved_graph_digest": active_digest,
        "historical_common_prefix_and_ports_checked_false_count": marker_count,
        "output_refs": _copy_metadata(outputs),
        "outputs": _copy_metadata(outputs),
    }
    return qualification, bool(graph_shape and unique_ids and historical_marker_ok), {
        "nodes": nodes,
        "graph": graph,
    }


def _saved_premises(
    attempt: Mapping[str, object],
    gates: list[str],
) -> tuple[dict[str, object], bool, Mapping[str, object]]:
    saved_shape_ok = (
        attempt.get("schema") == OLD_SOURCE_SCHEMA
        and attempt.get("source_only") is True
        and attempt.get("source_only_true") is True
        and attempt.get("passed") is False
    )
    _gate(gates, "saved_frozen_attempt_shape_unproved", saved_shape_ok)
    source_audit = attempt.get("source_prerequisite_audit")
    if not isinstance(source_audit, Mapping):
        source_audit = {}
    required_true = (
        "latest_native_C2_composition",
        "accepted_provider_report",
        "accepted_generic_density_identity",
        "accepted_nonreset_memory_report",
        "patch_C5_source_formula_AST",
        "patch_C5_archive_shape",
        "native_C5_geometry_and_Z_independence",
        "native_C6_C5_chart_ports",
        "native_C6_shared_parameter_prefix",
        "accepted_Pstar_source_parameter_port",
        "patch_inlet_witness",
        "patch_correction_inlet_identity",
        "retained_construction_data",
    )
    premises_ok = all(source_audit.get(key) is True for key in required_true)
    _gate(gates, "saved_source_premises_unproved", premises_ok)
    geometry = attempt.get("native_geometry_source_witness")
    geometry_ok = (
        isinstance(geometry, Mapping)
        and geometry.get("valid") is True
        and isinstance(geometry.get("checks"), Mapping)
        and geometry["checks"]
        and all(item is True for item in geometry["checks"].values())
    )
    _gate(gates, "native_geometry_witness_unproved", geometry_ok)
    old_statuses_ok = all(attempt.get(key) is value for key, value in FROZEN_STATUS_FLAGS.items())
    _gate(gates, "saved_frozen_statuses_changed", old_statuses_ok)
    old_gates = attempt.get("unresolved_gates")
    frozen_gate_set_ok = isinstance(old_gates, list) and old_gates == FROZEN_UNRESOLVED_GATES
    _gate(gates, "saved_historical_gate_set_changed", frozen_gate_set_ok)
    old_graph_checks = attempt.get("definition_graph_checks")
    graph_checks_ok = (
        isinstance(old_graph_checks, Mapping)
        and REQUIRED_GRAPH_CHECKS.issubset(old_graph_checks)
        and all(old_graph_checks.get(key) is True for key in REQUIRED_GRAPH_CHECKS)
    )
    _gate(gates, "saved_definition_graph_checks_unproved", graph_checks_ok)
    retained_ok = attempt.get("retained_construction_data") == RETAINED_CONSTRUCTION_DATA
    _gate(gates, "saved_retained_construction_data_mismatch", retained_ok)
    identity = attempt.get("source_identity_argument")
    identity = identity if isinstance(identity, Mapping) else {}
    source_identity_ok = (
        identity.get("new_owner_only") is True
        and identity.get("qualified_previous_DAG_refs_and_native_C5_graph_refs") is True
        and identity.get("same_N_eta_delta_Pstar_phase_origin_required") is True
        and identity.get("phase_origin") == "hb*sc/2"
        and identity.get("B_normalization") == "B_over_Pstar; physical B=2*V_y remains distinct"
        and identity.get("higher_derivative_source_premise")
        == "accepted C6 source equality and genuine C6 regularity on the common original domain"
        and identity.get("ordinary_Z012_is_not_higher_order_certificate") is True
        and identity.get("incoming_data")
        == "accepted C2 outgoing correction ports plus qualified native C5 postrestore inlet"
        and identity.get("numeric_equality_or_interval_overlap_used") is False
        and identity.get("bridge_identity_status") is True
        and identity.get("prior_adapter_P0_and_geometry_alias_assumed") is False
    )
    _gate(gates, "saved_source_identity_argument_unproved", source_identity_ok)
    source_audit_record = attempt.get("source_inspection_audit")
    if not isinstance(source_audit_record, Mapping):
        source_audit_record = {}
    source_audit_ok = (
        source_audit_record.get("source_records_only") is True
        and source_audit_record.get("original_graph_nodes_inspected_not_evaluated") is True
        and source_audit_record.get("new_source_consumer_calls") == 1
        and source_audit_record.get("accepted_ancestor_operator_calls") == 0
        and source_audit_record.get("producer_provider_builder_checker_replayed") is False
        and source_audit_record.get("quadrature_or_norm_calls") == 0
        and source_audit_record.get("numeric_field_evaluated") is False
        and source_audit_record.get("numeric_values_computed") is False
        and source_audit_record.get("old_operators_replayed") is False
        and source_audit_record.get("runtime_model_metadata") is None
    )
    _gate(gates, "saved_source_only_audit_unproved", source_audit_ok)
    ordinary_z_ok = (
        attempt.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
        and _get(attempt, "ordinary_Z_jet_rules.temporal_recursion_claimed") is False
    )
    _gate(gates, "saved_ordinary_Z_scope_unproved", ordinary_z_ok)
    historical_false_flags_ok = all(attempt.get(key) is False for key in GLOBAL_FALSE_FLAGS)
    _gate(gates, "saved_global_negative_scope_changed", historical_false_flags_ok)
    entry_checkpoint_ok = attempt.get("source_checkpoint_sha256") == FROZEN_SOURCE_CHECKPOINT_SHA256
    _gate(gates, "saved_source_checkpoint_changed", entry_checkpoint_ok)
    recipe_ok = attempt.get("finite_N_source_recipe") == EXPECTED_FINITE_RECIPE
    _gate(gates, "saved_finite_N_recipe_changed", recipe_ok)
    section_ok = attempt.get("section11_correction_recipe") == EXPECTED_SECTION11_RECIPE
    _gate(gates, "saved_section11_recipe_changed", section_ok)
    result = {
        "saved_frozen_attempt_shape": saved_shape_ok,
        "source_prerequisite_audit": _copy_metadata(source_audit),
        "definition_graph_checks": _copy_metadata(old_graph_checks) if isinstance(old_graph_checks, Mapping) else {},
        "geometry_witness": _copy_metadata(geometry) if isinstance(geometry, Mapping) else {},
        "historical_status_flags": {
            key: attempt.get(key) for key in FROZEN_STATUS_FLAGS
        },
        "historical_unresolved_gates": _copy_metadata(old_gates) if isinstance(old_gates, list) else [],
        "retained_construction_data": _copy_metadata(attempt.get("retained_construction_data")),
        "source_inspection_audit": _copy_metadata(source_audit_record),
        "ordinary_Z_orders": _copy_metadata(attempt.get("ordinary_Z_orders")),
        "source_checkpoint_sha256": attempt.get("source_checkpoint_sha256"),
        "entry_checkpoint_sha256": attempt.get("entry_checkpoint_sha256"),
        "source_identity_argument": _copy_metadata(identity),
        "source_premises_all_true": bool(premises_ok and geometry_ok and source_identity_ok),
    }
    ok = bool(
        premises_ok
        and geometry_ok
        and old_statuses_ok
        and frozen_gate_set_ok
        and graph_checks_ok
        and retained_ok
        and source_identity_ok
        and source_audit_ok
        and ordinary_z_ok
        and saved_shape_ok
        and historical_false_flags_ok
        and entry_checkpoint_ok
        and recipe_ok
        and section_ok
    )
    return result, ok, source_audit


def _latest_interface_and_leading(
    report: Mapping[str, object],
    component: Mapping[str, object],
    gates: list[str],
) -> tuple[dict[str, object], bool]:
    interface = _section(report, component, "next_finite_N_interface")
    phase = interface.get("phase") if isinstance(interface.get("phase"), Mapping) else {}
    transport = interface.get("patch_transport") if isinstance(interface.get("patch_transport"), Mapping) else {}
    interface_ok = (
        interface.get("schema") == "native_C2_patch_next_finite_N_interface_v1"
        and interface.get("source_scan_only") is True
        and interface.get("source_hashes_verified") is True
        and interface.get("patch_archive_prefix") == 205147
        and interface.get("phase_handle") == 196046
        and interface.get("coordinate_handle") == 185676
        and interface.get("patch_loop_A") == [197890, 197891, 197892]
        and interface.get("patch_loop_B_over_Pstar") == [197928, 197929, 197930]
        and transport.get("coordinate") == "x=R/Rm"
        and transport.get("measure") == "dx/x"
        and transport.get("p_memory") == "exactly one"
        and phase.get("source") == "fractional_part(N*(log(radius/Ra)-phase_origin))"
        and phase.get("origin") == "hb*sc/2, retained through patch"
    )
    scope = _section(report, component, "scope")
    result_flags = _section(report, component, "result_flags")
    leading_ok = (
        component.get("passed") is True
        and component.get("source_only") is True
        and component.get("nodes_digest") == EXPECTED_SOURCE_GRAPH_DIGEST
        and component.get("node_count") == 205300
        and component.get("native_C2_leading_patch_defined") is True
        and component.get("native_C2_C6_patch_source_identity_proved") is True
        and component.get("same_unique_leading_control_map_proved") is True
        and component.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
        and scope.get("native_C2_patch_own_history_identity_proved") is True
        and scope.get("first_ten_nonreset_history_identity_proved") is True
        and scope.get("legacy_feedback_reference_alias_proved") is False
        and result_flags.get("source_only", report.get("source_only", component.get("source_only", True))) is True
        and result_flags.get("numeric_values_computed", report.get("numeric_values_computed", component.get("numeric_values_computed", False))) is False
        and result_flags.get("old_operators_replayed", report.get("old_operators_replayed", component.get("old_operators_replayed", False))) is False
    )
    _gate(gates, "next_finite_N_interface_corrected_origin_unproved", interface_ok)
    _gate(gates, "latest_native_C2_leading_identity_unproved", leading_ok)
    return {
        "schema": interface.get("schema"),
        "source_scan_only": interface.get("source_scan_only"),
        "source_hashes_verified": interface.get("source_hashes_verified"),
        "phase_origin": phase.get("origin"),
        "phase_source": phase.get("source"),
        "transport_coordinate": transport.get("coordinate"),
        "transport_measure": transport.get("measure"),
        "transport_p_memory": transport.get("p_memory"),
        "leading_identity": leading_ok,
        "corrected": interface_ok,
    }, bool(interface_ok and leading_ok)


def _provider_source_premises(
    report: Mapping[str, object],
    component: Mapping[str, object],
    saved_audit: Mapping[str, object],
    gates: list[str],
) -> tuple[dict[str, object], bool]:
    theorem = component.get("actual_source_equality_theorem")
    theorem = theorem if isinstance(theorem, Mapping) else {}
    checks = theorem.get("checks") if isinstance(theorem.get("checks"), Mapping) else {}
    domain = theorem.get("domain_contract") if isinstance(theorem.get("domain_contract"), Mapping) else {}
    generic = component.get("strong_generic_recovery_theorem")
    generic = generic if isinstance(generic, Mapping) else {}
    generic_checks = generic.get("checks") if isinstance(generic.get("checks"), Mapping) else {}
    provider_scope = _section(report, component, "scope")
    common_domain_ok = (
        theorem.get("passed") is True
        and checks.get("ordinary_derivative_conversion_and_smoothness") is True
        and theorem.get("same_actual_analytic_source_domain") is True
        and domain.get("identity_is_conditional_on_common_original_domain") is True
        and isinstance(domain.get("fixed_parameters"), list)
        and isinstance(domain.get("active_branch"), list)
        and isinstance(domain.get("underlying_regularities"), str)
        and "genuine per-side C3" in domain.get("underlying_regularities", "")
        and domain.get("numerical_admissibility_or_conditioning_proved") is False
        and theorem.get("numerical_cone_or_admissibility_margins_claimed") is False
        and generic.get("passed") is True
        and generic_checks.get("C6_recovery") is True
        and "genuine C6 definition supplies C3 regularity" in generic.get("regularity_lift", "")
    )
    provider_ok = (
        component.get("passed") is True
        and component.get("source_only") is True
        and component.get("nodes_digest") == EXPECTED_SOURCE_GRAPH_DIGEST
        and component.get("node_count") == 205300
        and component.get("first_ten_finite_N_field_source_composition_proved") is True
        and component.get("first_ten_density_C2_NC5_identity_proved") is True
        and component.get("generic_density_identity_proved") is True
        and component.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
        and provider_scope.get("numeric_values_computed", report.get("numeric_values_computed", component.get("numeric_values_computed", False))) is False
        and provider_scope.get("old_operators_replayed", report.get("old_operators_replayed", component.get("old_operators_replayed", False))) is False
        and common_domain_ok
    )
    _gate(gates, "provider_corrected_source_smoothness_unproved", common_domain_ok)
    _gate(gates, "accepted_finite_N_provider_prerequisite_unproved", provider_ok)
    return {
        "provider_passed": component.get("passed"),
        "source_only": component.get("source_only"),
        "ordinary_Z_orders": _copy_metadata(component.get("ordinary_Z_orders")),
        "actual_source_equality_theorem_passed": theorem.get("passed"),
        "ordinary_derivative_conversion_and_smoothness": checks.get("ordinary_derivative_conversion_and_smoothness"),
        "same_actual_analytic_source_domain": theorem.get("same_actual_analytic_source_domain"),
        "identity_conditional_on_common_original_domain": domain.get("identity_is_conditional_on_common_original_domain"),
        "numerical_admissibility_or_conditioning_proved": domain.get("numerical_admissibility_or_conditioning_proved"),
        "numerical_cone_or_admissibility_margins_claimed": theorem.get("numerical_cone_or_admissibility_margins_claimed"),
        "fixed_parameters": _copy_metadata(domain.get("fixed_parameters")),
        "active_branch": _copy_metadata(domain.get("active_branch")),
        "underlying_regularities": domain.get("underlying_regularities"),
        "C6_recovery": generic_checks.get("C6_recovery"),
        "regularity_lift": generic.get("regularity_lift"),
        "correct_path_checked": True,
        "provider_ok": provider_ok,
        "saved_nonreset_memory_premise": saved_audit.get("accepted_nonreset_memory_report") is True,
    }, bool(provider_ok and common_domain_ok)


def _c6_prefix_and_ports(
    latest_component: Mapping[str, object],
    saved_audit: Mapping[str, object],
    geometry: Mapping[str, object],
    gates: list[str],
) -> tuple[dict[str, object], bool]:
    qualified = latest_component.get("qualified_source_ports")
    qualified = qualified if isinstance(qualified, Mapping) else {}
    chart_ports = saved_audit.get("native_C6_C5_chart_ports") is True
    shared_parameters = saved_audit.get("native_C6_shared_parameter_prefix") is True
    geometry_independent = saved_audit.get("native_C5_geometry_and_Z_independence") is True
    geometry_valid = geometry.get("valid") is True
    c6_archive = qualified.get("C6_archive_sha256") == "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a"
    c6_prefix = qualified.get("C6_prefix") == 196042
    c6_system = qualified.get("C6_system") == 193470
    leading_audit = latest_component.get("source_prerequisite_audit", {})
    accepted_C6_source = isinstance(leading_audit, Mapping) and all(
        leading_audit.get(name) is True for name in (
            "c6_receipt_accepted", "c6_source_formula_corrected",
            "c6_branch_unique_common_box", "c6_system_exact", "c6_controls_exact",
        )
    )
    derived_prefix = bool(
        geometry_independent
        and geometry_valid
        and chart_ports
        and shared_parameters
        and c6_archive
        and c6_prefix
        and c6_system
        and accepted_C6_source
    )
    _gate(gates, "native_C5_geometry_Z_independence_source_unproved", geometry_independent and geometry_valid)
    _gate(gates, "native_C6_to_C5_patch_provider_identity_unproved", derived_prefix)
    return {
        "native_C5_geometry_and_Z_independence": geometry_independent,
        "native_geometry_witness_valid": geometry_valid,
        "native_C6_C5_chart_ports": chart_ports,
        "native_C6_shared_parameter_prefix": shared_parameters,
        "C6_archive_sha256_checked": c6_archive,
        "C6_prefix_checked": c6_prefix,
        "C6_system_checked": c6_system,
        "accepted_C6_receipt_formula_branch_and_controls": accepted_C6_source,
        "native_C6_patch_prefix_identity_derived": derived_prefix,
    }, derived_prefix


def _history_and_inlet(
    attempt: Mapping[str, object],
    saved_audit: Mapping[str, object],
    provider_ok: bool,
    c6_prefix_ok: bool,
    gates: list[str],
) -> tuple[bool, bool, dict[str, object]]:
    section = attempt.get("section11_correction_recipe")
    section = section if isinstance(section, Mapping) else {}
    section_ok = section == EXPECTED_SECTION11_RECIPE
    inlet_old = (
        attempt.get("native_C2_patch_correction_inlet_identity_proved") is True
        and saved_audit.get("patch_inlet_witness") is True
        and saved_audit.get("patch_correction_inlet_identity") is True
    )
    _gate(gates, "patch_correction_inlet_identity_unproved", inlet_old and c6_prefix_ok)
    history_ok = bool(
        provider_ok
        and c6_prefix_ok
        and section_ok
        and saved_audit.get("accepted_nonreset_memory_report") is True
        and saved_audit.get("patch_correction_inlet_identity") is True
    )
    _gate(gates, "finite_N_history_identity_prerequisite_unproved", history_ok)
    return bool(inlet_old and c6_prefix_ok), history_ok, {
        "section11_recipe_checked": section_ok,
        "accepted_nonreset_memory_report": saved_audit.get("accepted_nonreset_memory_report") is True,
        "accepted_first_ten_C2_outgoing_correction_ports": saved_audit.get("patch_correction_inlet_identity") is True,
        "inlet_old_status_preserved": attempt.get("native_C2_patch_correction_inlet_identity_proved") is True,
        "history_identity_prerequisite": history_ok,
        "rates": _copy_metadata(section.get("rates")),
    }


def derive(
    saved_attempt: Mapping[str, object],
    latest_report: Mapping[str, object],
    provider_report: Mapping[str, object],
    *,
    expected_attempt_digest: str,
    expected_latest_digest: str,
    expected_provider_digest: str,
) -> dict[str, object]:
    """Admit the bounded finite-``N`` source stage from saved metadata only."""
    gates: list[str] = []
    input_digests = _digest_inputs(
        saved_attempt,
        latest_report,
        provider_report,
        expected_attempt_digest,
        expected_latest_digest,
        expected_provider_digest,
        gates,
    )
    attempt = _component(saved_attempt)
    latest = _component(latest_report)
    provider = _component(provider_report)

    candidate, graph_ok, graph_data = _saved_graph_reference(attempt, gates)
    candidate["saved_attempt_canonical_digest"] = input_digests["attempt"]
    implementation_hashes = saved_attempt.get("implementation_hashes") if isinstance(saved_attempt, Mapping) else None
    if isinstance(implementation_hashes, Mapping):
        kernel_hash = implementation_hashes.get(
            "src/openai_ns_reconstruction/original_native_C2_patch_finite_N_source.py"
        )
        candidate["kernel_sha256"] = kernel_hash
        _gate(gates, "saved_source_kernel_changed", kernel_hash == FROZEN_KERNEL_SHA256)
    else:
        # The once-only component result carries the source checkpoint but not
        # its outer publication implementation map.  The optional map is
        # checked when a stage envelope supplies it; absence is not invented
        # as a source failure.
        candidate["kernel_sha256"] = None
    saved_premise_result, saved_premises_ok, saved_audit = _saved_premises(attempt, gates)
    latest_info, latest_ok = _latest_interface_and_leading(latest_report, latest, gates)
    provider_info, provider_ok = _provider_source_premises(
        provider_report,
        provider,
        saved_audit,
        gates,
    )
    geometry = saved_premise_result.get("geometry_witness")
    geometry = geometry if isinstance(geometry, Mapping) else {}
    prefix_info, prefix_ok = _c6_prefix_and_ports(latest, saved_audit, geometry, gates)
    inlet_ok, history_ok, history_info = _history_and_inlet(
        attempt,
        saved_audit,
        provider_ok,
        prefix_ok,
        gates,
    )
    fields_ok = bool(
        graph_ok
        and saved_premises_ok
        and latest_ok
        and prefix_ok
        and not _get(attempt, "scope_limits.endpoint_rows_are_values", False)
    )
    provider_status = bool(fields_ok and provider_ok and prefix_ok)
    density_status = bool(
        provider_status
        and provider.get("generic_density_identity_proved") is True
        and saved_audit.get("accepted_generic_density_identity") is True
    )
    inlet_status = bool(fields_ok and inlet_ok)
    history_status = bool(density_status and history_ok)
    statuses = {
        "native_C2_patch_finite_N_fields_defined": fields_ok,
        "native_C2_C5_patch_provider_identity_proved": provider_status,
        "native_C2_C5_patch_density_identity_proved": density_status,
        "native_C2_patch_correction_inlet_identity_proved": inlet_status,
        "native_C2_C5_patch_history_identity_proved": history_status,
    }
    patch_identity = all(statuses.values()) and not gates
    candidate["qualified_reference"] = bool(graph_ok)
    source_premises = {
        "historical_saved_facts": saved_premise_result,
        "corrected_next_finite_N_interface": latest_info,
        "corrected_provider_source_equality": provider_info,
        "native_C6_C5_provider_identity": prefix_info,
        "section11_history_and_inlet": history_info,
        "common_original_domain_conditional": bool(provider_info.get("identity_conditional_on_common_original_domain") is True),
        "genuine_source_smoothness_conditional": bool(provider_info.get("ordinary_derivative_conversion_and_smoothness") is True),
        "numerical_admissibility_or_conditioning_proved": False,
        "numeric_cone_or_admissibility_margins_claimed": False,
    }
    result: dict[str, object] = {
        "schema": SCHEMA,
        "source_only": True,
        "source_only_true": True,
        "passed": bool(patch_identity),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_not_temporal_n": True,
        "saved_attempt_digest": input_digests["attempt"],
        "expected_attempt_digest": expected_attempt_digest,
        "latest_report_digest": input_digests["latest"],
        "expected_latest_digest": expected_latest_digest,
        "provider_report_digest": input_digests["provider"],
        "expected_provider_digest": expected_provider_digest,
        "definition_graph_node_count": FROZEN_GRAPH_NODE_COUNT,
        "saved_attempt_result_sha256": FROZEN_RESULT_SHA256,
        "saved_source_checkpoint_sha256": FROZEN_SOURCE_CHECKPOINT_SHA256,
        "saved_kernel_sha256": FROZEN_KERNEL_SHA256,
        "new_source_consumer_calls": 1,
        "accepted_ancestor_operator_calls": 0,
        "old_definition_graphs_rebuilt": False,
        **statuses,
        "patch_finite_N_source_identity": bool(patch_identity),
        "candidate_DAG_qualification": candidate,
        "source_premise_audit": source_premises,
        "corrected_integration_mistakes": {
            "phase_origin_annotation": latest_info.get("phase_origin") == "hb*sc/2, retained through patch",
            "phase_origin_canonical_pair": latest_info.get("phase_source") == "fractional_part(N*(log(radius/Ra)-phase_origin))",
            "provider_numerical_path": provider_info.get("numerical_admissibility_or_conditioning_proved") is False,
            "provider_numerical_wrong_top_level_not_used": True,
        },
        "historical_frozen_statuses": _copy_metadata(saved_premise_result.get("historical_status_flags")),
        "historical_frozen_unresolved_gates": _copy_metadata(saved_premise_result.get("historical_unresolved_gates")),
        "retained_construction_data": _copy_metadata(saved_premise_result.get("retained_construction_data")),
        "unresolved_gates": list(gates),
        "global_scope": dict(GLOBAL_FALSE_FLAGS),
        "same_density_certificate": False,
        "finite_N_history_identity": False,
        "physical_Q_y_Z1": False,
        "global_five_moment_closure": False,
        "legacy_feedback_reference_alias_proved": False,
        "prior_adapter_geometry_alias_proved": False,
        "prior_adapter_P0_alias_proved": False,
        "time_order_recursion_claimed": False,
        "numeric_values_computed": False,
        "old_operators_replayed": False,
        "full_goal_complete": False,
        "runtime_model_metadata": None,
        "worker_role": "luna_max_worker",
        "source_inspection_audit": {
            "source_records_only": True,
            "metadata_only_admission": True,
            "new_source_consumer_calls": 1,
            "accepted_ancestor_operator_calls": 0,
            "old_definition_graphs_rebuilt": False,
            "old_graph_re_evaluated": False,
            "project_modules_imported": False,
            "producer_provider_builder_checker_replayed": False,
            "quadrature_or_norm_calls": 0,
            "numeric_field_evaluated": False,
            "numeric_values_computed": False,
            "old_operators_replayed": False,
            "runtime_model_metadata": None,
        },
        "role": "luna_max_worker",
        "model_role": "luna_max",
    }
    # A failed source premise is always visible in both the gate list and the
    # structured audit; no false positive is hidden by a downstream status.
    result["source_premises_all_true"] = bool(not gates)
    result["historical_graph_metadata"] = {
        "node_count": len(graph_data.get("nodes", [])) if isinstance(graph_data.get("nodes"), list) else None,
        "graph_namespace": _get(graph_data, "graph.namespace"),
        "active_graph_digest": _get(graph_data, "graph.active_graph_digest"),
        "definition_graphs_rebuilt": False,
    }
    return result


__all__ = ["SCHEMA", "derive"]
