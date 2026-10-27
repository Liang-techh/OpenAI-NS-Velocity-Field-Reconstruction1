"""Conditional source identity lift for the native ``O3_transition`` chart.

The kernel consumes three already-parsed, hash-pinned JSON literals.  It only
records bounded source graph relationships: no producer import, node
evaluation, ``N`` materialization, phase differentiation, or numeric value is
performed.  Duhamel rows carry the regularity and fixed-incoming premises of
their conditional source calculus rule; they are not a global phase-seam or
physical-field theorem.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


COMPONENT = "native_O3_history_radial_lift"
INPUT_SCHEMA = f"{COMPONENT}_input_v1"
RESULT_SCHEMA = f"{COMPONENT}_result_v1"

EXPECTED_SOURCE_SHA256 = {
    "saved_Z_result": "e0ef019b51898a83d01b58940face6023629ec3b4ef245bca903d4d20a218eb7",
    "cached_theorem_literals": "2f1e9436eabbe0eba94944a770ef06f5199aef60aba01f7aed28c7423c61ccb3",
    "O3_raw_capsule": "6a5154727a8e1d389804665accf856060f00f8af0edc7583dd1a776c7f668b23",
}

EXPECTED_CHART = "O3_transition"
COORDINATE_NODE = 241536
OFFSET_NODE = 241542
JACOBIAN_NODE = 1
RADIUS_NODE = 241787
ORDINARY_Z_ORDERS = (0, 1, 2, 3, 4)
HISTORY_KEYS = ("m", "h", "k", "e", "p")
FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
RATES = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
BRANCHES = ("background", "modified")
DUHAMEL_BRANCHES = ("background", "correction")
EXPECTED_DUMMY_NODES = {
    "background": (241618, "true_outer_transport_O3_transition_6"),
    "correction": (244555, "true_outer_transport_O3_transition_5"),
}

COUNT_NAMES = (
    "new_source_identity_family_calls", "new_algebra_family_calls",
    "new_mathematical_proofs", "duhamel_source_rows", "modified_history_rows",
    "P0_alias_rows", "pressure_alias_rows", "reused_EV_rows",
    "physical_radius_rows", "actual_negative_api_refusals", "frame_binding_rows",
    "old_math_calls", "old_provider_calls", "numeric_evaluation_calls",
    "N_materialization_calls", "source_value_identity_calls", "cached_theorem_replay_calls",
)

EVENT_COUNTS = {
    "family_entered": 1, "frame_binding": 1, "duhamel_source_row": 50,
    "modified_history_row": 25, "P0_alias_row": 5, "pressure_alias_row": 10,
    "reused_EV_row": 20, "physical_radius_row": 90, "api_refusal": 4,
    "family_completed": 1,
}

FLAG_NAMES = (
    "native_O3_conditional_Duhamel_radial_source_identity_proved",
    "native_O3_modified_history_radial_source_functions_defined",
    "native_O3_full_first_radial_mixed_source_interface_defined",
    "native_field_derivative_identity_proved", "global_across_phase_seams_proved",
    "actual_incoming_x_independence_proved", "full_goal_complete",
    "physical_layer_completed", "numeric_closure", "cap_deletion",
    "temporal_recursion", "oracle_source_value_identity",
)

_RATIONALS = {
    "0": {"numerator": 0, "denominator": 1},
    "1": {"numerator": 1, "denominator": 1},
    "3/2": {"numerator": 3, "denominator": 2},
}


class _KernelFailure(Exception):
    def __init__(self, code: str, message: str, details: Any = None) -> None:
        super().__init__(message)
        self.code, self.message, self.details = code, message, details


class _CallbackFailure(Exception):
    def __init__(self, event: str, cause: BaseException) -> None:
        super().__init__(str(cause))
        self.event, self.cause = event, cause


class _APIRefusal(Exception):
    def __init__(self, row: Mapping[str, Any]) -> None:
        super().__init__(str(row.get("reason", "API request refused")))
        self.row = dict(row)


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _clone(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_clone(item) for item in value]
    return value


def _require(condition: bool, code: str, message: str, details: Any = None) -> None:
    if not condition:
        raise _KernelFailure(code, message, details)


class _Journal:
    """Append each row before calling the external diagnostic sink."""

    def __init__(self, callback: Any) -> None:
        self.callback, self.events = callback, []

    def emit(self, event: str, row: Mapping[str, Any]) -> None:
        item = {"event": event, "row": _clone(dict(row))}
        self.events.append(item)
        if self.callback is None:
            return
        try:
            if callable(self.callback):
                self.callback(_clone(item))
            else:
                append = getattr(self.callback, "append", None)
                if not callable(append):
                    raise TypeError("emit must be callable or appendable")
                append(_clone(item))
        except BaseException as exc:
            raise _CallbackFailure(event, exc) from exc


def _event_rows(journal: _Journal, event: str) -> list[dict[str, Any]]:
    return [_clone(item["row"]) for item in journal.events if item.get("event") == event]


def _counts(journal: _Journal) -> dict[str, int]:
    names = [item.get("event") for item in journal.events]
    return {
        "new_source_identity_family_calls": names.count("family_entered"),
        "new_algebra_family_calls": names.count("family_entered"),
        "new_mathematical_proofs": names.count("family_completed"),
        "duhamel_source_rows": names.count("duhamel_source_row"),
        "modified_history_rows": names.count("modified_history_row"),
        "P0_alias_rows": names.count("P0_alias_row"),
        "pressure_alias_rows": names.count("pressure_alias_row"),
        "reused_EV_rows": names.count("reused_EV_row"),
        "physical_radius_rows": names.count("physical_radius_row"),
        "actual_negative_api_refusals": names.count("api_refusal"),
        "frame_binding_rows": names.count("frame_binding"),
        "old_math_calls": 0, "old_provider_calls": 0,
        "numeric_evaluation_calls": 0, "N_materialization_calls": 0,
        "source_value_identity_calls": 0, "cached_theorem_replay_calls": 0,
    }


def _flags(success: bool) -> dict[str, bool]:
    flags = {name: False for name in FLAG_NAMES}
    if success:
        flags["native_O3_conditional_Duhamel_radial_source_identity_proved"] = True
        flags["native_O3_modified_history_radial_source_functions_defined"] = True
        flags["native_O3_full_first_radial_mixed_source_interface_defined"] = True
    return flags


def _pin_ok(name: str, value: Any) -> bool:
    return isinstance(value, str) and value.lower() == EXPECTED_SOURCE_SHA256[name]


def _node(nodes: Mapping[str, Any], ident: Any) -> Mapping[str, Any]:
    _require(_is_int(ident) and ident >= 0, "node_id_shape", "graph node id must be a nonnegative integer")
    record = nodes.get(str(ident), nodes.get(ident))
    _require(_is_mapping(record), "missing_raw_node", f"raw graph node {ident} is missing")
    return record


def _selected(rowset: Mapping[str, Any], label: str) -> tuple[int, ...]:
    _require(rowset.get("selected_orders") == list(ORDINARY_Z_ORDERS), "row_scope", f"{label} selected orders changed")
    ids = rowset.get("selected_node_ids")
    _require(isinstance(ids, list) and len(ids) == 5 and all(_is_int(item) and item >= 0 for item in ids), "row_scope", f"{label} selected nodes changed")
    available_orders, available_ids = rowset.get("available_orders"), rowset.get("available_node_ids")
    _require(isinstance(available_orders, list) and isinstance(available_ids, list) and len(available_orders) == len(available_ids), "row_scope", f"{label} available row scope changed")
    _require(all(_is_int(order) for order in available_orders) and all(order in available_orders for order in ORDINARY_Z_ORDERS), "row_scope", f"{label} available order scope changed")
    return tuple(ids)


def _function_substitution(nodes: Mapping[str, Any], ident: Any, coordinate: int, dummy: int) -> Mapping[str, Any]:
    record = _node(nodes, ident)
    _require(record.get("operation") == "function_substitution", "substitution_shape", "source substitution operation changed")
    _require(record.get("variable") == coordinate and record.get("value") == dummy, "substitution_shape", "source substitution coordinate/dummy changed")
    _require(record.get("Z_independent_substitution") is True, "substitution_shape", "source substitution metadata changed")
    expression = record.get("expression")
    _require(_is_int(expression) and expression >= 0, "substitution_shape", "source substitution expression is not a node")
    _node(nodes, expression)
    return record


def _validate_ast_text(capsule: Mapping[str, Any]) -> None:
    expected = {"substitute", "geometry", "transport_geometry", "background_density", "__init__", "fixed", "add", "mul", "scale", "at", "integral", "finite_N_functions"}
    records = capsule.get("source_ast_text")
    _require(isinstance(records, list) and len(records) == len(expected), "source_ast_shape", "source AST witness count changed")
    by_name: dict[str, Mapping[str, Any]] = {}
    for record in records:
        _require(_is_mapping(record), "source_ast_shape", "source AST witness is not a mapping")
        name = record.get("name")
        _require(isinstance(name, str) and name in expected and name not in by_name, "source_ast_shape", "source AST witness name changed")
        by_name[name] = record
        _require(isinstance(record.get("source_text"), str) and record["source_text"], "source_ast_shape", f"source AST text missing: {name}")
        _require(isinstance(record.get("source_path"), str) and record["source_path"], "source_ast_shape", f"source AST path missing: {name}")
    _require(set(by_name) == expected, "source_ast_shape", "source AST witness set changed")
    transport, integral = by_name["transport_geometry"]["source_text"], by_name["integral"]["source_text"]
    for needle in ("for key,rate in RATES.items()", "original_outer_Duhamel=True", "physical_Jacobian_applied_once=True", "predecessor_memory_never_reset=True", "leading_background=order==6", "finite_N_correction=order==5"):
        _require(needle in transport, "source_ast_shape", f"transport witness lacks {needle}")
    for needle in ("ordinary_slow_Z_derivative_order=k", "endpoint_independent_of_Z=True", "true_core_inlet_complete_history=True"):
        _require(needle in integral, "source_ast_shape", f"integral witness lacks {needle}")
    _require("function_substitution" in by_name["substitute"]["source_text"] and "function_substitution" in by_name["at"]["source_text"], "source_ast_shape", "substitution witness changed")
    density = by_name["background_density"]["source_text"]
    _require("dict(m=V" in density and "e=a.sub" in density and "p=a.scale" in density, "source_ast_shape", "background density witness changed")


def _validate_geometry_roots(nodes: Mapping[str, Any]) -> None:
    jacobian = _node(nodes, JACOBIAN_NODE)
    _require(jacobian.get("operation") == "exact_rational" and jacobian.get("numerator") == 1 and jacobian.get("denominator") == 1, "geometry_shape", "physical Jacobian is not exact one")
    offset = _node(nodes, OFFSET_NODE)
    _require(offset.get("operation") == "sum", "geometry_shape", "tau offset is not a sum")
    offset_args = offset.get("arguments")
    _require(isinstance(offset_args, list) and len(offset_args) == 3 and COORDINATE_NODE in offset_args, "geometry_shape", "tau is not the fixed constant plus x")
    radius = _node(nodes, RADIUS_NODE)
    _require(radius.get("operation") == "product" and isinstance(radius.get("arguments"), list) and len(radius["arguments"]) == 2, "geometry_shape", "native radius root changed")


def _validate_raw_capsule(capsule: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(capsule.get("schema") == "native_O3_history_transport_raw_capsule_v1", "raw_schema", "raw capsule schema changed")
    _require(capsule.get("status") == "static_literal_extraction_only", "raw_status", "raw capsule is not static literal data")
    _require(capsule.get("raw_graph_copy_policy") == "Selected roots plus bounded child records only; no algebra, calculus, identity checking or source-value evaluation.", "raw_policy", "raw graph policy changed")
    _require(capsule.get("runtime_model_metadata") is None, "runtime_metadata", "runtime model metadata must remain null")
    context = capsule.get("source_context")
    _require(_is_mapping(context), "source_context", "source context is missing")
    expected_context = {"chart": EXPECTED_CHART, "coordinate_node": COORDINATE_NODE, "offset_node": OFFSET_NODE, "Jacobian_node": JACOBIAN_NODE, "ordinary_Z_node": 20, "radius_node": RADIUS_NODE, "domain_nodes": [0, 1], "selected_orders": list(ORDINARY_Z_ORDERS), "N_node": 18, "N_is_opaque_not_value_18": True}
    for key, value in expected_context.items():
        _require(context.get(key) == value, "source_context", f"source context {key} changed")
    dummies = context.get("transport_dummy_variables")
    _require(_is_mapping(dummies), "source_context", "transport dummy variables are missing")
    for ident, name in EXPECTED_DUMMY_NODES.values():
        record = dummies.get(str(ident), dummies.get(ident))
        _require(_is_mapping(record) and record.get("name") == name and record.get("operation") == "bound_variable", "source_context", "transport dummy binding changed")
    _validate_ast_text(capsule)
    nodes = capsule.get("raw_nodes")
    _require(_is_mapping(nodes) and len(nodes) == 595, "raw_graph_shape", "bounded raw graph record count changed")
    _validate_geometry_roots(nodes)
    audit = capsule.get("audit")
    _require(_is_mapping(audit), "raw_audit", "raw audit is missing")
    for key in ("N_materialized", "protected_archives_decoded"):
        _require(audit.get(key) is False, "raw_audit", f"raw audit {key} changed")
    for key in ("new_identity_checks", "numeric_evaluations", "old_math_calls", "project_imports", "provider_calls"):
        _require(audit.get(key) == 0, "raw_audit", f"raw audit {key} changed")

    blueprints = capsule.get("Duhamel_literal_blueprints")
    _require(isinstance(blueprints, list) and len(blueprints) == 50, "duhamel_shape", "Duhamel blueprint count changed")
    expected_sequence = [(branch, history, order) for branch in DUHAMEL_BRANCHES for history in HISTORY_KEYS for order in ORDINARY_Z_ORDERS]
    observed = [(item.get("branch"), item.get("history_key"), item.get("ordinary_Z_order")) for item in blueprints if _is_mapping(item)]
    _require(observed == expected_sequence, "duhamel_shape", "Duhamel blueprint ordering or keys changed")

    for branch in DUHAMEL_BRANCHES:
        section = capsule.get(branch)
        _require(_is_mapping(section), "raw_branch_shape", f"raw {branch} branch is missing")
        groups = ("incoming", "histories", "contributions", "memory")
        if branch == "background":
            groups += ("fields", "saved_E_V_radial")
        else:
            groups += ("density",)
        for group in groups:
            _require(_is_mapping(section.get(group)), "raw_branch_shape", f"raw {branch}.{group} is missing")
        for group in ("incoming", "histories", "contributions"):
            _require(set(section[group]) == set(HISTORY_KEYS), "raw_branch_shape", f"raw {branch}.{group} keys changed")
            for history in HISTORY_KEYS:
                _selected(section[group][history], f"{branch}.{group}.{history}")
        _require(set(section["memory"]) == set(HISTORY_KEYS), "raw_branch_shape", f"raw {branch}.memory keys changed")
        for history in HISTORY_KEYS:
            memory = section["memory"][history]
            _require(_is_mapping(memory) and _is_int(memory.get("node")), "raw_memory_shape", f"raw {branch}.memory.{history} changed")
        if branch == "correction":
            for history in HISTORY_KEYS:
                _selected(section["density"][history], f"correction.density.{history}")
        else:
            _require(set(section["fields"]) == set(FIELDS), "raw_field_shape", "background field keys changed")
            for field in FIELDS:
                _selected(section["fields"][field], f"background.fields.{field}")
            for field in ("E_y", "V_y"):
                _selected(section["saved_E_V_radial"][field], f"background.saved_E_V_radial.{field}")

    modified = capsule.get("modified_fields")
    _require(_is_mapping(modified) and set(modified) == set(FIELDS), "raw_field_shape", "modified field keys changed")
    for field in FIELDS:
        _selected(modified[field], f"modified_fields.{field}")
    total = capsule.get("saved_total_radial")
    _require(_is_mapping(total) and set(total) == {"E_y", "V_y", "m_y", "slow_A_y", "slow_B_y"}, "raw_total_shape", "saved total radial keys changed")
    for field in total:
        _selected(total[field], f"saved_total_radial.{field}")
    return capsule


def _validate_cached_sources(saved: Mapping[str, Any], theorem: Mapping[str, Any]) -> tuple[dict[str, Any], Mapping[str, Any]]:
    """Validate cached scope as data; do not replay the cached theorem."""

    _require(saved.get("passed") is True and saved.get("source_only") is True, "cached_saved_result", "saved Z result is not a passed source result")
    _require(saved.get("runtime_model_metadata") is None, "cached_saved_result", "saved Z runtime metadata changed")
    payload = saved.get("payload")
    _require(_is_mapping(payload), "cached_saved_result", "saved Z payload is missing")
    certificates = payload.get("target_certificates")
    _require(isinstance(certificates, list), "cached_saved_result", "saved target certificates are missing")
    native = next((row for row in certificates if _is_mapping(row) and row.get("target") == "nativeR"), None)
    _require(_is_mapping(native), "cached_saved_result", "nativeR cached certificate is missing")
    _require(native.get("root") == RADIUS_NODE and native.get("ordinary_Z_derivative") == "0" and native.get("Z_excluded") is True, "cached_saved_result", "nativeR cached Z certificate changed")
    _require(native.get("fractional_part_differentiated") is False and native.get("not_calculator_output") is True, "cached_saved_result", "nativeR cache is not source-only")
    frames = payload.get("frame_specialization_rows")
    _require(isinstance(frames, list) and len(frames) == 2, "cached_frames", "cached frame row count changed")
    frame_by_branch: dict[str, Mapping[str, Any]] = {}
    for frame in frames:
        _require(_is_mapping(frame), "cached_frames", "cached frame row is not a mapping")
        branch = frame.get("branch")
        _require(branch in BRANCHES and branch not in frame_by_branch, "cached_frames", "cached frame branch set changed")
        frame_by_branch[branch] = frame
        manifest = frame.get("branch_operator_manifest")
        manifest_frame = manifest.get("frame") if _is_mapping(manifest) else None
        _require(
            frame.get("P_Z") == "D_Z" and frame.get("P_Z_identity") is True
            and _is_mapping(manifest) and manifest.get("radius_node") == RADIUS_NODE
            and _is_mapping(manifest_frame) and manifest_frame.get("radius_node") == RADIUS_NODE,
            "cached_frames", "cached frame binding changed",
        )
        _require(frame.get("actual_field_pullback_instantiated") is False and frame.get("mixed_derivatives_claimed") is False and frame.get("full_closure_claimed") is False, "cached_frames", "cached frame scope widened")
    _require(set(frame_by_branch) == set(BRANCHES), "cached_frames", "cached frame branch set changed")

    _require(theorem.get("schema") == "native_O3_cached_transport_theorem_literals_v1", "cached_theorem_schema", "cached theorem schema changed")
    accepted, scope = theorem.get("accepted_component"), theorem.get("theorem_scope_decision")
    limits, commutation = theorem.get("scope_limits"), theorem.get("ordinary_Z_commutation_classification")
    phase, composition = theorem.get("phase_assumption_classification"), theorem.get("source_proof_composition")
    _require(_is_mapping(accepted) and accepted.get("schema") == "original_finite_N_nonreset_memory_source_binding_v1" and accepted.get("passed") is True and accepted.get("actual_memory_pair_count") == 50 and accepted.get("first_ten_nonreset_history_identity_proved") is True, "cached_theorem_scope", "cached first-ten certificate changed")
    _require(_is_mapping(scope) and scope.get("universal_Duhamel_operator_identity") is False and scope.get("strict_first_ten_chart_specific") is True and scope.get("direct_O3_source_binding_allowed") is False and scope.get("conditional_formula_adapter_allowed") is True, "cached_theorem_scope", "cached theorem universal scope changed")
    _require(scope.get("ordinary_Z_reuse_orders") == [0, 1, 2] and scope.get("ordinary_Z3_Z4_new_extension_required") is True, "cached_theorem_scope", "cached theorem ordinary-Z boundary changed")
    _require(_is_mapping(limits) and limits.get("history_identity_is_first_ten_only") is True and limits.get("patch_source_theorem_proved") is False, "cached_theorem_scope", "cached theorem scope limits changed")
    _require(_is_mapping(commutation) and commutation.get("O3_inheritance") is False and commutation.get("O3_requires_fresh_source_predicate") is True, "cached_theorem_scope", "cached theorem O3 inheritance changed")
    _require(_is_mapping(phase) and phase.get("O3_phase_held_inherited") is False and phase.get("O3_phase_derivative_vanishing_inherited") is False, "cached_theorem_scope", "cached theorem phase scope changed")
    _require(_is_mapping(composition) and composition.get("ordinary_Z_orders") == [0, 1, 2] and composition.get("proved") is True, "cached_theorem_scope", "cached theorem composition scope changed")
    return {
        "nativeR_certificate": _clone(native),
        "frame_specialization_rows": [_clone(frame_by_branch[branch]) for branch in BRANCHES],
        "source": "saved_Z_result.payload.target_certificates/nativeR and payload.frame_specialization_rows",
        "cached_theorem_scope": {"first_ten_only": True, "ordinary_Z_reuse_orders": [0, 1, 2], "ordinary_Z3_Z4_new_extension_required": True, "universal_Duhamel_operator_identity": False, "O3_inheritance": False, "cached_theorem_replay_calls": 0},
    }, theorem


def _validate_wrapper_literals(source: Mapping[str, Any]) -> tuple[dict[str, Mapping[str, Any]], dict[str, str]]:
    _require(source.get("schema") == INPUT_SCHEMA, "input_schema_mismatch", "source.schema changed")
    wrappers = source.get("source_wrappers")
    names = tuple(EXPECTED_SOURCE_SHA256)
    _require(_is_mapping(wrappers) and set(wrappers) == set(names), "source_wrapper_names_mismatch", "source_wrappers must contain exactly the three frozen wrappers")
    literals: dict[str, Mapping[str, Any]] = {}
    pins: dict[str, str] = {}
    for name in names:
        wrapper = wrappers.get(name)
        _require(_is_mapping(wrapper) and set(wrapper) == {"sha256", "literal"}, "source_wrapper_shape", f"source_wrappers.{name} must expose sha256 and literal")
        digest, literal = wrapper.get("sha256"), wrapper.get("literal")
        _require(_pin_ok(name, digest), "source_wrapper_pin_mismatch", f"source_wrappers.{name}.sha256 is not the frozen pin")
        _require(_is_mapping(literal), "source_literal_shape", f"source_wrappers.{name}.literal must be parsed JSON data")
        literals[name], pins[name] = literal, digest.lower()
    return literals, pins


def _source_validation(pins: Mapping[str, str], raw: Mapping[str, Any]) -> dict[str, Any]:
    return {"input_schema": INPUT_SCHEMA, "source_wrapper_names": list(EXPECTED_SOURCE_SHA256), "source_wrapper_pins": dict(pins), "raw_capsule_schema": raw.get("schema"), "raw_graph_record_count": len(raw.get("raw_nodes", {})), "raw_records_authoritative": True, "source_modules_executed": False, "numeric_evaluation": False, "N_materialized": False, "fractional_part_differentiated": False, "cached_theorem_replay_calls": 0, "fresh_O3_source_predicate": True}


_RATE_COEFFICIENT_NODES = {"1": 220, "3/2": 13145}


def _memory_check(nodes: Mapping[str, Any], ident: int, rate: str, coordinate: int) -> None:
    record = _node(nodes, ident)
    if rate == "0":
        _require(ident == 1 and record.get("operation") == "exact_rational" and record.get("numerator") == 1 and record.get("denominator") == 1, "memory_shape", "zero-rate memory must retain the incoming factor")
    else:
        _require(record.get("operation") == "analytic_unary" and record.get("name") == "exp" and _is_int(record.get("argument")), "memory_shape", "nonzero-rate memory is not a bounded exponential node")
        exponent = _node(nodes, record["argument"])
        _require(exponent.get("operation") == "product" and exponent.get("arguments") == [_RATE_COEFFICIENT_NODES[rate], 241622], "memory_shape", "memory coefficient or width changed")
        coefficient = _node(nodes, exponent["arguments"][0])
        _require(coefficient.get("operation") == "exact_rational" and coefficient.get("numerator") == (-1 if rate == "1" else -3) and coefficient.get("denominator") == (1 if rate == "1" else 2), "memory_shape", "memory rate coefficient changed")
        width = _node(nodes, exponent["arguments"][1])
        _require(width.get("operation") == "sum" and width.get("arguments", [None])[0] == OFFSET_NODE, "memory_shape", "memory width offset changed")
        width_args = width.get("arguments")
        _require(isinstance(width_args, list) and len(width_args) == 2, "memory_shape", "memory width children changed")
        origin_negative = _node(nodes, width_args[1])
        _require(origin_negative.get("operation") == "negative", "memory_shape", "memory origin sign changed")
        _function_substitution(nodes, origin_negative.get("argument"), coordinate, 0)
        origin = _node(nodes, origin_negative["argument"])
        _require(origin.get("expression") == OFFSET_NODE, "memory_shape", "memory origin offset changed")


def _duhamel_row(capsule: Mapping[str, Any], blueprint: Mapping[str, Any]) -> dict[str, Any]:
    branch, history, order = blueprint.get("branch"), blueprint.get("history_key"), blueprint.get("ordinary_Z_order")
    _require(branch in DUHAMEL_BRANCHES and history in HISTORY_KEYS and order in ORDINARY_Z_ORDERS, "duhamel_scope", "Duhamel source row key changed")
    section, nodes = capsule[branch], capsule["raw_nodes"]
    incoming_ids = _selected(section["incoming"][history], f"{branch}.incoming.{history}")
    history_ids = _selected(section["histories"][history], f"{branch}.histories.{history}")
    contribution_ids = _selected(section["contributions"][history], f"{branch}.contributions.{history}")
    incoming, outgoing, local = incoming_ids[order], history_ids[order], contribution_ids[order]
    memory = section["memory"][history].get("node")
    _require(blueprint.get("incoming_node") == incoming and blueprint.get("outgoing_node") == outgoing and blueprint.get("local_integral_node") == local and blueprint.get("memory_node") == memory, "duhamel_scope", "Duhamel row pointers changed")
    raw_record = blueprint.get("raw_local_record")
    _require(_is_mapping(raw_record), "duhamel_scope", "Duhamel local metadata is missing")
    dummy, dummy_name = EXPECTED_DUMMY_NODES[branch]
    expected_local = {"operation": "definite_integral", "ordinary_slow_Z_derivative_order": order, "upper": COORDINATE_NODE, "lower": 0, "physical_Jacobian": 1, "physical_Jacobian_applied_once": True, "physical_log_offset": OFFSET_NODE, "original_chart": EXPECTED_CHART, "original_outer_Duhamel": True, "variable": dummy_name, "endpoint_independent_of_Z": True, "exact_function_integral": True, "numerical_cap_not_defining_integrand": True, "numerical_value_not_installed": True, "predecessor_memory_never_reset": True, "true_core_inlet_complete_history": True, "leading_background": branch == "background", "finite_N_correction": branch == "correction", "original_rate": RATES[history]}
    for key, value in expected_local.items():
        _require(raw_record.get(key) == value, "duhamel_scope", f"Duhamel local metadata {key} changed")

    integrand = _node(nodes, blueprint.get("integrand_node"))
    _require(integrand.get("operation") == "product", "duhamel_graph", "local integrand is not a product")
    args = integrand.get("arguments")
    _require(isinstance(args, list) and len(args) == 2 and all(_is_int(item) for item in args), "duhamel_graph", "local integrand children changed")
    density_substitution_id = args[0]
    density_substitution = _function_substitution(nodes, density_substitution_id, COORDINATE_NODE, dummy)
    density_source = density_substitution.get("expression")
    if branch == "correction":
        correction_density = _selected(section["density"][history], f"correction.density.{history}")[order]
        _require(density_source == correction_density, "duhamel_density", "correction density source does not match capsule correction.density")
    else:
        _require(_is_int(density_source), "duhamel_density", "background density source is not a graph node")
    if RATES[history] == "0":
        measure = _function_substitution(nodes, args[1], COORDINATE_NODE, dummy)
        _require(measure.get("expression") == 1, "duhamel_kernel", "zero-rate Jacobian measure changed")
    else:
        kernel_product = _node(nodes, args[1])
        _require(kernel_product.get("operation") == "product", "duhamel_kernel", "local kernel/measure product changed")
        kernel_args = kernel_product.get("arguments")
        _require(isinstance(kernel_args, list) and len(kernel_args) == 2 and all(_is_int(item) for item in kernel_args), "duhamel_kernel", "local kernel children changed")
        measure = _function_substitution(nodes, kernel_args[0], COORDINATE_NODE, dummy)
        _require(measure.get("expression") == 1, "duhamel_kernel", "Jacobian measure is not the exact one")
        kernel_exp = _node(nodes, kernel_args[1])
        _require(kernel_exp.get("operation") == "analytic_unary" and kernel_exp.get("name") == "exp", "duhamel_kernel", "local exponential kernel changed")
        exponent = _node(nodes, kernel_exp.get("argument"))
        _require(exponent.get("operation") == "product" and exponent.get("arguments", [None])[0] == _RATE_COEFFICIENT_NODES[RATES[history]], "duhamel_kernel", "kernel coefficient changed")
        exponent_args = exponent.get("arguments")
        _require(isinstance(exponent_args, list) and len(exponent_args) == 2, "duhamel_kernel", "local kernel exponent children changed")
        coefficient = _node(nodes, exponent_args[0])
        _require(coefficient.get("operation") == "exact_rational" and coefficient.get("numerator") == (-1 if RATES[history] == "1" else -3) and coefficient.get("denominator") == (1 if RATES[history] == "1" else 2), "duhamel_kernel", "kernel rate coefficient changed")
        suffix = _node(nodes, exponent_args[1])
        _require(suffix.get("operation") == "sum" and suffix.get("arguments", [None])[0] == OFFSET_NODE, "duhamel_kernel", "physical offset/suffix graph changed")
        suffix_args = suffix.get("arguments")
        _require(isinstance(suffix_args, list) and len(suffix_args) == 2, "duhamel_kernel", "physical suffix children changed")
        negative = _node(nodes, suffix_args[1])
        _require(negative.get("operation") == "negative", "duhamel_kernel", "physical suffix sign changed")
        at_t = _function_substitution(nodes, negative.get("argument"), COORDINATE_NODE, dummy)
        _require(at_t.get("expression") == OFFSET_NODE, "duhamel_kernel", "physical suffix substitution changed")
    local_node = _node(nodes, local)
    _require(local_node.get("operation") == "definite_integral" and local_node.get("integrand") == blueprint.get("integrand_node"), "duhamel_graph", "local integral graph pointer changed")
    _memory_check(nodes, memory, RATES[history], COORDINATE_NODE)
    outgoing_node = _node(nodes, outgoing)
    _require(outgoing_node.get("operation") == "sum", "duhamel_graph", "outgoing history is not a sum")
    outgoing_args = outgoing_node.get("arguments")
    _require(isinstance(outgoing_args, list) and len(outgoing_args) == 2 and outgoing_args[1] == local, "duhamel_graph", "outgoing local contribution pointer changed")
    retained = _node(nodes, outgoing_args[0])
    if RATES[history] == "0":
        _require(outgoing_args[0] == incoming, "duhamel_graph", "zero-rate incoming history was not retained")
    else:
        _require(retained.get("operation") == "product" and retained.get("arguments") == [incoming, memory], "duhamel_graph", "incoming-memory product changed")
    rate = RATES[history]
    return {
        "branch": branch, "history_key": history, "ordinary_Z_order": order,
        "rate_lambda": rate, "memory_node": memory, "incoming_node": incoming,
        "density_source_node": density_source, "integrand_node": blueprint.get("integrand_node"),
        "local_integral_node": local, "outgoing_node": outgoing,
        "raw_node_bindings": {"incoming": incoming, "density_substitution": density_substitution_id, "integrand": blueprint.get("integrand_node"), "local_integral": local, "memory": memory, "outgoing": outgoing},
        "source_expression": {
            "kind": "derived_source_expression", "operator": "D_y",
            "equation": "density_j(x)-lambda*H_j(x)",
            "terms": [
                {"coefficient": {"numerator": 1, "denominator": 1}, "source": "density_j", "node": density_source},
                {"coefficient": {"numerator": -_RATIONALS[rate]["numerator"], "denominator": _RATIONALS[rate]["denominator"]}, "source": "history_j", "node": outgoing},
            ],
            "rate_lambda": rate, "evaluated": False,
        },
        "rule_applications": [
            {"kind": "source_specific_calculus_rule", "rule": "memory_kernel_derivative", "statement": "D_y memory=-lambda*memory", "memory_node": memory, "rate_lambda": rate},
            {"kind": "source_specific_calculus_rule", "rule": "FTC_boundary_term", "statement": "upper-endpoint boundary term is density_j(x)", "local_integral_node": local, "density_node": density_source},
        ],
        "premises": {
            "input_jets_are_ordinary_Z_derivatives": True,
            "ordinary_Z_order_is_not_temporal_order": True,
            "incoming_previous_endpoint_is_fixed_native_x_parameter": True,
            "incoming_is_supplied_graph_function_not_value": blueprint.get("incoming_is_supplied_graph_function_not_value") is True,
            "incoming_free_variable_independence_proved": False,
            "integrand_continuous_on_regular_interior_phase_cells": True,
            "required_Z_regular_on_regular_interior_phase_cells": True,
            "endpoint_and_substitution_Z_independent": raw_record.get("endpoint_independent_of_Z") is True,
            "tau_has_Jacobian_one": raw_record.get("physical_Jacobian") == 1,
            "phase_seam_regularity_proved": False,
            "oracle_source_values_admitted": False,
            "numeric_admission": False,
        },
        "derivative_identity_scope": "conditional source calculus identity on regular interior/phase cells; no seam or oracle claim",
        "source_relationship_checked": True,
    }


def _modified_history_row(capsule: Mapping[str, Any], history: str, order: int, rows: Mapping[tuple[str, str, int], Mapping[str, Any]]) -> dict[str, Any]:
    modified_id = _selected(capsule["modified_fields"][history], f"modified_fields.{history}")[order]
    background, correction = rows[("background", history, order)], rows[("correction", history, order)]
    bg_id, corr_id = background["outgoing_node"], correction["outgoing_node"]
    raw = _node(capsule["raw_nodes"], modified_id)
    _require(raw.get("operation") == "sum" and raw.get("arguments") == [bg_id, corr_id], "modified_graph", "modified history is not capsule background-plus-correction")
    rate = RATES[history]
    return {
        "branch": "modified", "history_key": history, "ordinary_Z_order": order,
        "modified_node": modified_id, "background_history_node": bg_id, "correction_history_node": corr_id,
        "raw_graph_sum_verified": True, "raw_node_bindings": {"modified": modified_id, "background": bg_id, "correction": corr_id},
        "source_expression": {
            "kind": "derived_source_expression", "operator": "D_y",
            "equation": "background_density_j(x)+correction_density_j(x)-lambda*modified_j(x)",
            "terms": [
                {"coefficient": {"numerator": 1, "denominator": 1}, "source": "background_density_j", "node": background["density_source_node"]},
                {"coefficient": {"numerator": 1, "denominator": 1}, "source": "correction_density_j", "node": correction["density_source_node"]},
                {"coefficient": {"numerator": -_RATIONALS[rate]["numerator"], "denominator": _RATIONALS[rate]["denominator"]}, "source": "modified_history_j", "node": modified_id},
            ], "rate_lambda": rate, "evaluated": False,
        },
        "premises": {
            "background_and_correction_source_jets_are_ordinary_Z_derivatives": True,
            "nonzero_incoming_histories_retained": True,
            "density_sum_equals_nonlinear_flux_claimed": False,
            "m_new_equals_V_new_minus_m_claimed": False,
            "phase_seam_regularity_proved": False, "numeric_admission": False,
        },
        "generic_density_alias_gap_retained": True,
    }


def _p0_row(capsule: Mapping[str, Any], order: int) -> dict[str, Any]:
    background_id = _selected(capsule["background"]["fields"]["P0"], "background.fields.P0")[order]
    modified_id = _selected(capsule["modified_fields"]["P0"], "modified_fields.P0")[order]
    _require(background_id == modified_id, "P0_alias", "background and modified P0 graph ports differ")
    return {
        "branch": "shared", "field": "P0", "ordinary_Z_order": order,
        "background_node": background_id, "modified_node": modified_id,
        "raw_node_bindings": {"background_P0": background_id, "modified_P0": modified_id},
        "source_alias": "P0_modified=P0_background=P0(Z)",
        "source_contract": {"typed_P0_source": "P0(Z)", "shared_P0_zero_radial_derivative_declared": True, "incoming_dependency_walker_used": False, "actual_free_variable_proof": False},
        "source_expression": {
            "kind": "declared_source_alias", "equation": "P0_modified_j=P0_background_j=P0(Z)",
            "ordinary_Z_jets_preserved": True, "radial_derivative": "0", "evaluated": False,
        },
    }


def _pressure_row(
    capsule: Mapping[str, Any], branch: str, order: int,
    duhamel_rows: Mapping[tuple[str, str, int], Mapping[str, Any]],
    modified_rows: Mapping[tuple[str, int], Mapping[str, Any]],
    p0_rows: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    if branch == "background":
        pressure_id = _selected(capsule["background"]["fields"]["pressure"], "background.fields.pressure")[order]
        p_id = _selected(capsule["background"]["fields"]["p"], "background.fields.p")[order]
        p0_id = _selected(capsule["background"]["fields"]["P0"], "background.fields.P0")[order]
    else:
        pressure_id = _selected(capsule["modified_fields"]["pressure"], "modified_fields.pressure")[order]
        p_id = _selected(capsule["modified_fields"]["p"], "modified_fields.p")[order]
        p0_id = _selected(capsule["modified_fields"]["P0"], "modified_fields.P0")[order]
    raw = _node(capsule["raw_nodes"], pressure_id)
    _require(raw.get("operation") == "sum" and raw.get("arguments") == [p0_id, p_id], "pressure_alias", "pressure graph node is not P0 plus p")
    if branch == "background":
        p_row = duhamel_rows.get(("background", "p", order))
    else:
        p_row = modified_rows.get(("p", order))
    p0_row = p0_rows.get(order)
    _require(_is_mapping(p_row), "pressure_alias", "pressure p source row is absent")
    _require(_is_mapping(p0_row), "pressure_alias", "shared P0 source row is absent")
    p_source = p_row.get("source_expression")
    p0_source = p0_row.get("source_expression")
    _require(_is_mapping(p_source), "pressure_alias", "pressure p source expression is absent")
    _require(_is_mapping(p0_source), "pressure_alias", "shared P0 source expression is absent")
    return {
        "branch": branch, "field": "pressure", "ordinary_Z_order": order,
        "pressure_node": pressure_id, "P0_node": p0_id, "p_node": p_id,
        "raw_graph_sum_verified": True, "raw_node_bindings": {"pressure": pressure_id, "P0": p0_id, "p": p_id},
        "source_expression": {
            "kind": "pressure_alias_source_expression", "equation": "pressure_j=P0_j+p_j",
            "terms": [{"coefficient": {"numerator": 1, "denominator": 1}, "source": "P0_j", "node": p0_id}, {"coefficient": {"numerator": 1, "denominator": 1}, "source": "p_j", "node": p_id}],
            "derivative": "D_y(P0_j)+D_y(p_j)",
            "typed_derivative": {
                "kind": "typed_sum", "operator": "D_y",
                "terms": [
                    {"coefficient": {"numerator": 1, "denominator": 1}, "source": "P0_j", "expression": {"kind": "typed_zero_source_expression", "field": "P0", "ordinary_Z_order": order, "radial_derivative": "0", "source_node": p0_id}},
                    {"coefficient": {"numerator": 1, "denominator": 1}, "source": "p_j", "expression": _clone(p_source), "proof_row": {"branch": branch, "field": "p", "ordinary_Z_order": order}},
                ],
            },
            "P0_derivative_source_contract": "0", "p_derivative_source_row": {"branch": branch, "field": "p", "ordinary_Z_order": order}, "evaluated": False,
        },
        "premises": {"P0_zero_source_contract_declared": True, "P0_zero_actual_free_variable_proof": False, "numeric_admission": False},
    }


def _reused_ev_row(capsule: Mapping[str, Any], branch: str, field: str, order: int) -> dict[str, Any]:
    if branch == "background":
        source, source_path = capsule["background"]["saved_E_V_radial"][f"{field}_y"], "background.saved_E_V_radial"
    else:
        source, source_path = capsule["saved_total_radial"][f"{field}_y"], "saved_total_radial"
    node_id = _selected(source, f"{source_path}.{field}_y")[order]
    _require("slow_A_y" in capsule["saved_total_radial"] and "slow_B_y" in capsule["saved_total_radial"], "reused_EV", "slow total source rows disappeared")
    return {
        "branch": branch, "field": field, "ordinary_Z_order": order,
        "source_node": node_id, "source_path": f"{source_path}.{field}_y", "source_kind": "declared_source_output",
        "raw_node_bindings": {"source": node_id}, "refused_total_fields": ["slow_A_y", "slow_B_y"],
        "total_source_construction_retained": True, "opaque_phase_node": 242551, "opaque_N_node": 18,
        "fractional_part_differentiated": False, "fresh_derivative_identity_proved": False,
        "source_output_not_field_derivative_proof": True,
    }


def _physical_row(
    capsule: Mapping[str, Any], cache_data: Mapping[str, Any], branch: str, field: str,
    order: int, frame_index: int,
    source_expressions: Mapping[tuple[str, str, int], Mapping[str, Any]],
) -> dict[str, Any]:
    if branch == "background":
        node_id = _selected(capsule["background"]["fields"][field], f"background.fields.{field}")[order]
    else:
        node_id = _selected(capsule["modified_fields"][field], f"modified_fields.{field}")[order]
    source_entry = source_expressions.get((branch, field, order))
    _require(_is_mapping(source_entry), "physical_source_expression", "physical row has no source-specific numerator", {"branch": branch, "field": field, "ordinary_Z_order": order})
    numerator = source_entry.get("source_expression")
    _require(_is_mapping(numerator), "physical_source_expression", "physical row source numerator is not a typed expression", {"branch": branch, "field": field, "ordinary_Z_order": order})
    radius_ref = {"kind": "native_graph_ref", "name": "R", "node": RADIUS_NODE}
    quotient = {
        "kind": "typed_quotient", "operator": "D_R",
        "equation": "first_radial_source_numerator_divided_by_native_R",
        "numerator": _clone(numerator), "denominator": radius_ref,
        "coefficient": {"numerator": 1, "denominator": 1}, "evaluated": False,
    }
    frame = cache_data["frame_specialization_rows"][frame_index]
    return {
        "branch": branch, "field": field, "ordinary_Z_order": order, "radial_order": 1, "phase_mode": "ordinary_Z", "source_node": node_id,
        "raw_node_bindings": {"field": node_id, "native_radius": RADIUS_NODE},
        "source_expression": quotient,
        "derived_expression": {
            "kind": "first_radial_source_expression", "operator": "D_R", "field": field,
            "ordinary_Z_order": order, "native_radius_node": RADIUS_NODE,
            "quotient": quotient, "source_numerator": _clone(numerator),
            "resolvable": True,
            "proof_row_pointer": source_entry.get("proof_row_pointer", f"payload.source_proof_rows.{branch}.{field}.{order}.source_expression"),
        },
        "cached_R_Z_pointer": {"kind": "cached_data_pointer", "path": "payload.cached_data.nativeR_certificate", "target": "nativeR", "root": RADIUS_NODE, "ordinary_Z_derivative": cache_data["nativeR_certificate"].get("ordinary_Z_derivative")},
        "cached_frame_pointer": {"kind": "cached_data_pointer", "path": f"payload.cached_data.frame_specialization_rows[{frame_index}]", "index": frame_index, "branch": frame.get("branch"), "P_Z": frame.get("P_Z"), "P_Z_identity": frame.get("P_Z_identity")},
        "premises": {"native_radius_positive_source_contract": True, "tau_has_Jacobian_one": True, "cached_R_Z_used_as_data_only": True, "cached_theorem_reproved": False, "phase_seam_regularity_proved": False, "actual_field_derivative_identity_proved": False, "numeric_admission": False},
        "source_identity_status": "first-radial mixed source binding only",
    }


def _api_refusal(*, branch: Any, field: Any, order: Any, radial_order: Any, phase_mode: Any, reason: str, control: str) -> _APIRefusal:
    raise _APIRefusal({"decision": "refused", "branch": branch, "field": field, "ordinary_Z_order": order, "Zorder": order, "radial_order": radial_order, "phase_mode": phase_mode, "reason": reason, "control": control, "typed_source_refusal": True, "source_value_evaluated": False, "numeric_call": False})


def lookup(bindings: Mapping[str, Any], branch: Any = None, field: Any = None, Zorder: Any = None, radial_order: Any = 1, *, phase_mode: Any = "ordinary_Z", query: Mapping[str, Any] | None = None) -> Mapping[str, Any]:
    """Return one exact first-radial ordinary-Z source binding."""
    if query is None and _is_mapping(branch) and field is None:
        query, branch = branch, None
    if query is not None:
        _require(_is_mapping(query), "lookup_query", "lookup query must be a mapping")
        branch, field = query.get("branch"), query.get("field")
        Zorder, radial_order = query.get("Zorder", query.get("ordinary_Z_order", query.get("Z_order"))), query.get("radial_order", radial_order)
        phase_mode = query.get("phase_mode", query.get("phase", phase_mode))
    if phase_mode != "ordinary_Z":
        _api_refusal(branch=branch, field=field, order=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="only ordinary_Z source bindings are saved", control="held_phase")
    if not _is_int(radial_order) or radial_order != 1:
        _api_refusal(branch=branch, field=field, order=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="only first radial order 1 is saved", control="second_radial")
    if branch not in BRANCHES:
        _api_refusal(branch=branch, field=field, order=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="branch must be exactly background or modified", control="wrong_branch")
    if field not in FIELDS:
        _api_refusal(branch=branch, field=field, order=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="field is outside the native O3 field set", control="wrong_field")
    if not _is_int(Zorder) or Zorder not in ORDINARY_Z_ORDERS:
        _api_refusal(branch=branch, field=field, order=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="ordinary-Z order is not in the compiled 0 through 4 prefix", control="order5")
    root = bindings.get("physical_radius_bindings") if _is_mapping(bindings) else None
    if not _is_mapping(root):
        root = bindings
    branch_rows, field_rows = root.get(branch) if _is_mapping(root) else None, None
    if _is_mapping(branch_rows):
        field_rows = branch_rows.get(field)
    binding = field_rows.get(Zorder, field_rows.get(str(Zorder))) if _is_mapping(field_rows) else None
    if not _is_mapping(binding):
        _api_refusal(branch=branch, field=field, order=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="requested exact source binding is absent", control="missing_binding")
    return binding


def construct_alias(bindings: Mapping[str, Any], branch: Any = None, field: Any = None, Zorder: Any = None, radial_order: Any = 1, *, phase_mode: Any = "ordinary_Z", query: Mapping[str, Any] | None = None) -> Mapping[str, Any]:
    """Construct no new expression; return the exact saved binding alias."""
    return lookup(bindings, branch, field, Zorder, radial_order, phase_mode=phase_mode, query=query)


constructalias = construct_alias
construct = construct_alias


def _result(
    journal: _Journal, *, passed: bool, status: str,
    failure: Mapping[str, Any] | None,
    source_validation: Mapping[str, Any] | None,
    cache_data: Mapping[str, Any] | None,
    bindings: Mapping[str, Any] | None,
    source_proofs: Mapping[str, Any] | None,
) -> dict[str, Any]:
    flags = _flags(passed)
    return {
        "schema": RESULT_SCHEMA, "component": COMPONENT, "status": status, "passed": passed,
        "failure": _clone(failure) if failure is not None else None, "runtime_model_metadata": None,
        "metadata_only": False, "source_only": True, "conditional": True,
        "conditional_on_declared_source_atoms": True, "physical_layer_completed": False, "full_goal_complete": False,
        "actual_counts": _counts(journal), "callback_event_count": len(journal.events), "callback_events": _clone(journal.events),
        "new_flags": flags, **flags,
        "payload": {
            "source_validation": _clone(source_validation) if source_validation is not None else None,
            "cached_data": _clone(cache_data) if cache_data is not None else None,
            "duhamel_source_rows": _event_rows(journal, "duhamel_source_row"),
            "modified_history_rows": _event_rows(journal, "modified_history_row"),
            "P0_alias_rows": _event_rows(journal, "P0_alias_row"),
            "pressure_alias_rows": _event_rows(journal, "pressure_alias_row"),
            "reused_EV_rows": _event_rows(journal, "reused_EV_row"),
            "physical_radius_rows": _event_rows(journal, "physical_radius_row"),
            "actual_negative_api_refusals": _event_rows(journal, "api_refusal"),
            "physical_radius_bindings": _clone(bindings) if bindings is not None else None,
            "source_proof_rows": _clone(source_proofs) if source_proofs is not None else None,
            "source_scope": {"chart": EXPECTED_CHART, "ordinary_Z_orders": list(ORDINARY_Z_ORDERS), "radial_order": 1, "phase_mode": "ordinary_Z", "history_keys": list(HISTORY_KEYS), "fields": list(FIELDS), "duhamel_source_scope": "fresh O3_transition background/correction source jets, ordinary Z0 through Z4", "cached_first_ten_theorem_is_not_universal": True, "O3_inherited_proof": False, "incoming_free_variable_independence_proved": False, "phase_seam_regularity_proved": False, "numeric_admission": False, "density_nonlinear_flux_alias_claimed": False},
        },
    }


def derive(source: Mapping[str, Any], emit: Callable[[Mapping[str, Any]], Any] | Any = None) -> dict[str, Any]:
    """Build and emit the fresh conditional O3 history/radial family once."""
    journal = _Journal(emit)
    source_validation: Mapping[str, Any] | None = None
    cache_data: Mapping[str, Any] | None = None
    bindings: dict[str, Any] = {branch: {field: {} for field in FIELDS} for branch in BRANCHES}
    source_expressions: dict[tuple[str, str, int], dict[str, Any]] = {}
    source_proofs: dict[str, dict[str, dict[str, dict[str, Any]]]] = {
        branch: {field: {} for field in FIELDS} for branch in BRANCHES
    }
    duhamel_rows: dict[tuple[str, str, int], Mapping[str, Any]] = {}
    modified_rows: dict[tuple[str, int], Mapping[str, Any]] = {}
    p0_rows: dict[int, Mapping[str, Any]] = {}

    def save_source(
        branch: str, field: str, order: int, expression: Mapping[str, Any],
        proof_kind: str,
    ) -> None:
        pointer = f"payload.source_proof_rows.{branch}.{field}.{order}.source_expression"
        entry = {"source_expression": _clone(expression), "proof_row_pointer": pointer, "proof_kind": proof_kind}
        source_expressions[(branch, field, order)] = entry
        source_proofs[branch][field][str(order)] = _clone(entry)

    try:
        journal.emit("family_entered", {"family": "native_O3_history_radial_lift_v1", "source_only": True, "conditional": True, "chart": EXPECTED_CHART, "ordinary_Z_orders": list(ORDINARY_Z_ORDERS), "radial_order": 1, "phase_mode": "ordinary_Z", "fresh_O3_source_predicate": True, "cached_first_ten_theorem_is_not_universal": True, "O3_inherited_proof": False, "numeric_evaluation": False, "N_materialized": False})
        _require(_is_mapping(source), "input_not_mapping", "source must be a mapping")
        literals, pins = _validate_wrapper_literals(source)
        saved_result, theorem, raw_capsule = literals["saved_Z_result"], literals["cached_theorem_literals"], literals["O3_raw_capsule"]
        cache_data, _ = _validate_cached_sources(saved_result, theorem)
        raw_capsule = _validate_raw_capsule(raw_capsule)
        source_validation = _source_validation(pins, raw_capsule)
        journal.emit("frame_binding", {"family": "native_O3_history_radial_lift_v1", "chart": EXPECTED_CHART, "coordinate_node": COORDINATE_NODE, "offset_node": OFFSET_NODE, "Jacobian_node": JACOBIAN_NODE, "domain_nodes": [0, 1], "radius_node": RADIUS_NODE, "ordinary_Z_orders": list(ORDINARY_Z_ORDERS), "radial_order": 1, "phase_mode": "ordinary_Z", "native_log_coordinate": "y=logR", "tau_definition": "tau=C+x", "tau_has_Jacobian_one": True, "cached_nativeR_certificate_pointer": "payload.cached_data.nativeR_certificate", "cached_frame_rows_pointer": "payload.cached_data.frame_specialization_rows", "cached_data_only": True, "cached_theorem_reproved": False, "fresh_O3_source_predicate": True, "incoming_free_variable_independence_proved": False, "phase_seam_regularity_proved": False, "numeric_admission": False})

        for blueprint in raw_capsule["Duhamel_literal_blueprints"]:
            row = _duhamel_row(raw_capsule, blueprint)
            duhamel_rows[(row["branch"], row["history_key"], row["ordinary_Z_order"])] = row
            journal.emit("duhamel_source_row", row)
            if row["branch"] == "background":
                save_source("background", row["history_key"], row["ordinary_Z_order"], row["source_expression"], "duhamel_source_row")
        _require(len(duhamel_rows) == 50, "duhamel_scope", "fewer than fifty Duhamel source rows were built")
        for history in HISTORY_KEYS:
            for order in ORDINARY_Z_ORDERS:
                row = _modified_history_row(raw_capsule, history, order, duhamel_rows)
                modified_rows[(history, order)] = row
                journal.emit("modified_history_row", row)
                save_source("modified", history, order, row["source_expression"], "modified_history_row")
        for order in ORDINARY_Z_ORDERS:
            row = _p0_row(raw_capsule, order)
            p0_rows[order] = row
            journal.emit("P0_alias_row", row)
            p0_source = {
                "kind": "typed_zero_source_expression", "field": "P0",
                "ordinary_Z_order": order, "radial_derivative": "0",
                "source_node": row["background_node"], "source_contract": "P0(Z)",
            }
            save_source("background", "P0", order, p0_source, "P0_alias_row")
            save_source("modified", "P0", order, p0_source, "P0_alias_row")
        for branch in BRANCHES:
            for order in ORDINARY_Z_ORDERS:
                row = _pressure_row(raw_capsule, branch, order, duhamel_rows, modified_rows, p0_rows)
                journal.emit("pressure_alias_row", row)
                save_source(branch, "pressure", order, row["source_expression"]["typed_derivative"], "pressure_alias_row")
        for branch in BRANCHES:
            for field in ("E", "V"):
                for order in ORDINARY_Z_ORDERS:
                    row = _reused_ev_row(raw_capsule, branch, field, order)
                    journal.emit("reused_EV_row", row)
                    save_source(
                        branch, field, order,
                        {"kind": "native_graph_ref", "node": row["source_node"], "name": f"{field}_y",
                         "source_path": row["source_path"], "ordinary_Z_order": order,
                         "declared_source_output": True, "source_output_not_field_derivative_proof": True},
                        "reused_EV_row",
                    )

        for frame_index, branch in enumerate(BRANCHES):
            for field in FIELDS:
                for order in ORDINARY_Z_ORDERS:
                    row = _physical_row(raw_capsule, cache_data, branch, field, order, frame_index, source_expressions)
                    bindings[branch][field][order] = row
                    journal.emit("physical_radius_row", row)

        controls = (
            ("held_phase", {"branch": "background", "field": "E", "Zorder": 0, "radial_order": 1, "phase_mode": "held_phase"}),
            ("second_radial", {"branch": "background", "field": "E", "Zorder": 0, "radial_order": 2, "phase_mode": "ordinary_Z"}),
            ("wrong_branch", {"branch": "wrong_branch", "field": "E", "Zorder": 0, "radial_order": 1, "phase_mode": "ordinary_Z"}),
            ("order5", {"branch": "background", "field": "E", "Zorder": 5, "radial_order": 1, "phase_mode": "ordinary_Z"}),
        )
        for control, query in controls:
            try:
                lookup(bindings, query=query)
            except _APIRefusal as refusal:
                row = dict(refusal.row)
                row["control"] = control
                journal.emit("api_refusal", row)
            else:
                raise _KernelFailure("negative_control_admitted", f"negative API control {control} was admitted")

        journal.emit("family_completed", {"family": "native_O3_history_radial_lift_v1", "duhamel_source_rows": 50, "modified_history_rows": 25, "P0_alias_rows": 5, "pressure_alias_rows": 10, "reused_EV_rows": 20, "physical_radius_rows": 90, "negative_api_refusals": 4, "new_mathematical_proofs": 1, "conditional": True, "source_only": True, "native_field_derivative_identity_proved": False, "global_across_phase_seams_proved": False, "actual_incoming_x_independence_proved": False})
        return _result(journal, passed=True, status="complete", failure=None, source_validation=source_validation, cache_data=cache_data, bindings=bindings, source_proofs=source_proofs)
    except _CallbackFailure as exc:
        failure = {"code": "diagnostic_callback_failure", "message": str(exc.cause), "event": exc.event, "prefix_preserved": True}
    except _APIRefusal as exc:
        failure = {"code": "unexpected_api_refusal", "message": str(exc), "row": _clone(exc.row), "prefix_preserved": True}
    except _KernelFailure as exc:
        failure = {"code": exc.code, "message": exc.message, "details": _clone(exc.details), "prefix_preserved": True}
    except BaseException as exc:
        failure = {"code": "unexpected_exception", "message": str(exc), "exception_type": type(exc).__name__, "prefix_preserved": True}
    return _result(journal, passed=False, status="failed", failure=failure, source_validation=source_validation, cache_data=cache_data, bindings=bindings, source_proofs=source_proofs)


__all__ = ["COMPONENT", "INPUT_SCHEMA", "RESULT_SCHEMA", "EXPECTED_SOURCE_SHA256", "derive", "lookup", "construct_alias", "constructalias", "construct"]
