"""Fresh, inert source binding for the native O3 density/flux identity.

This leaf deliberately works on JSON-shaped source receipts.  It does not
import a producer, materialise ``N``, evaluate a source value, decode an
archive, or call one of the older jet/algebra providers.  The ordinary-Z
rows are typed source expressions under an explicitly declared C4 source-jet
contract; they are not a global phase-seam or numerical field theorem.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


COMPONENT = "native_O3_density_flux_identity"
INPUT_SCHEMA = f"{COMPONENT}_input_v1"
RESULT_SCHEMA = f"{COMPONENT}_result_v1"

EXPECTED_SOURCE_SHA256 = {
    "saved_history_result": "0521dc68caec6e046c0de1b252c518869d20250c824c8e461de1a149d3955eab",
    "density_flux_literals": "ae3a6afbdc8de9c557dff18b3552a6b1c6c1d902decca5298b0b6adbcc216b59",
    "cached_density_scope": "9612240a921b20b006589a7def78cc7900e17c2b6507e7cddb233bfcedba0639",
}

EXPECTED_PARENT_COMMIT = "92b6358946ed026b2a97deb57cd2e43abc984e54"
EXPECTED_PARENT_REPORT_SHA256 = "138f7e61f677294dcadf06e57b903893b61a2bc17afa9cb0d7769490199338fa"
EXPECTED_CHART = "O3_transition"
RADIUS_NODE = 241787
N_NODE = 18
INVERSE_N_NODE = 99411
ORDINARY_Z_ORDERS = (0, 1, 2, 3, 4)
HISTORY_KEYS = ("m", "h", "k", "e", "p")
FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
BRANCHES = ("background", "modified")
RATES = {"m": "1", "h": "3/2", "k": "3/2", "e": "1", "p": "0"}
RATE_RATIONALS = {
    "m": (-1, 1),
    "h": (-3, 2),
    "k": (-3, 2),
    "e": (-1, 1),
    "p": (0, 1),
}
BINOMIALS = (
    (1,),
    (1, 1),
    (1, 2, 1),
    (1, 3, 3, 1),
    (1, 4, 6, 4, 1),
)

COUNT_NAMES = (
    "new_source_identity_family_calls",
    "new_algebra_family_calls",
    "new_mathematical_proofs",
    "scope_binding_rows",
    "analytic_alias_rows",
    "density_flux_rows",
    "derivative_rewrite_rows",
    "saved_m_y_alias_rows",
    "physical_binding_rows",
    "rewritten_physical_bindings",
    "retained_physical_bindings",
    "actual_negative_api_refusals",
    "old_math_calls",
    "old_provider_calls",
    "numeric_evaluation_calls",
    "N_materialization_calls",
    "source_value_reducer_calls",
    "cached_algebra_replay_calls",
)

FLAG_NAMES = (
    "native_O3_exp_expm1_functional_alias_proved",
    "native_O3_density_flux_source_identities_proved",
    "native_O3_saved_m_y_compatibility_proved",
    "native_O3_flux_radial_interface_defined",
    "unconditional_ordinary_Z4_source_admission",
    "native_field_derivative_identity_proved",
    "actual_incoming_x_independence_proved",
    "global_across_phase_seams_proved",
    "oracle_source_value_identity",
    "physical_layer_completed",
    "full_goal_complete",
    "numeric_closure",
    "cap_deletion",
    "temporal_recursion",
)


class _KernelFailure(Exception):
    def __init__(self, code: str, message: str, details: Any = None) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


class _CallbackFailure(Exception):
    def __init__(self, event: str, cause: BaseException) -> None:
        super().__init__(str(cause))
        self.event = event
        self.cause = cause


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
    if isinstance(value, list):
        return [_clone(item) for item in value]
    if isinstance(value, tuple):
        return [_clone(item) for item in value]
    return value


def _require(condition: bool, code: str, message: str, details: Any = None) -> None:
    if not condition:
        raise _KernelFailure(code, message, details)


class _Journal:
    """Record a row before handing it to the durable callback."""

    def __init__(self, callback: Any) -> None:
        self.callback = callback
        self.events: list[dict[str, Any]] = []

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
    rows = journal.events
    names = [item.get("event") for item in rows]
    physical = [item.get("row") for item in rows if item.get("event") == "physical_binding_row"]
    return {
        "new_source_identity_family_calls": names.count("family_entered"),
        "new_algebra_family_calls": names.count("family_entered"),
        "new_mathematical_proofs": names.count("family_completed"),
        "scope_binding_rows": names.count("scope_binding"),
        "analytic_alias_rows": names.count("analytic_alias_row"),
        "density_flux_rows": names.count("density_flux_row"),
        "derivative_rewrite_rows": names.count("derivative_rewrite_row"),
        "saved_m_y_alias_rows": names.count("saved_m_y_alias_row"),
        "physical_binding_rows": names.count("physical_binding_row"),
        "rewritten_physical_bindings": sum(
            1 for row in physical if _is_mapping(row) and row.get("action") == "rewritten_flux"
        ),
        "retained_physical_bindings": sum(
            1 for row in physical if _is_mapping(row) and row.get("action") == "retained_parent"
        ),
        "actual_negative_api_refusals": names.count("api_refusal"),
        "old_math_calls": 0,
        "old_provider_calls": 0,
        "numeric_evaluation_calls": 0,
        "N_materialization_calls": 0,
        "source_value_reducer_calls": 0,
        "cached_algebra_replay_calls": 0,
    }


def _flags(success: bool) -> dict[str, bool]:
    flags = {name: False for name in FLAG_NAMES}
    if success:
        for name in (
            "native_O3_exp_expm1_functional_alias_proved",
            "native_O3_density_flux_source_identities_proved",
            "native_O3_saved_m_y_compatibility_proved",
            "native_O3_flux_radial_interface_defined",
        ):
            flags[name] = True
    return flags


def _pin_ok(name: str, observed: Any) -> bool:
    return isinstance(observed, str) and observed.lower() == EXPECTED_SOURCE_SHA256[name]


def _wrapper_literals(source: Mapping[str, Any]) -> tuple[dict[str, Mapping[str, Any]], dict[str, Any]]:
    _require(source.get("schema") == INPUT_SCHEMA, "input_schema_mismatch", "source schema changed")
    _require(source.get("component") == COMPONENT, "input_component_mismatch", "source component changed")
    wrappers = source.get("source_wrappers")
    _require(_is_mapping(wrappers), "source_wrappers_missing", "source_wrappers must be a mapping")
    _require(set(wrappers) == set(EXPECTED_SOURCE_SHA256), "source_wrapper_names_mismatch", "source wrapper names changed")
    literals: dict[str, Mapping[str, Any]] = {}
    pins: dict[str, str] = {}
    for name, expected in EXPECTED_SOURCE_SHA256.items():
        wrapper = wrappers.get(name)
        _require(_is_mapping(wrapper), "source_wrapper_shape", f"source_wrappers.{name} must be a mapping")
        _require(set(wrapper) == {"sha256", "literal"}, "source_wrapper_shape", f"source_wrappers.{name} keys changed")
        observed = wrapper.get("sha256")
        _require(_pin_ok(name, observed), "source_wrapper_pin_mismatch", f"source_wrappers.{name} SHA changed")
        literal = wrapper.get("literal")
        _require(_is_mapping(literal), "source_literal_shape", f"source_wrappers.{name} literal must be a mapping")
        literals[name] = literal
        pins[name] = str(observed).lower()
    policy = source.get("policy")
    _require(_is_mapping(policy), "policy_shape", "policy must be a mapping")
    _require(policy.get("metadata_only") is False, "metadata_policy_mismatch", "metadata-only policy changed")
    _require(policy.get("source_only") is True, "source_only_policy_mismatch", "source-only policy changed")
    _require(policy.get("full_goal_complete") is False, "goal_policy_mismatch", "full-goal policy changed")
    _require(policy.get("numeric_evaluation") is False, "numeric_policy_mismatch", "numeric policy changed")
    if source.get("parent_commit") is not None:
        _require(source.get("parent_commit") == EXPECTED_PARENT_COMMIT, "parent_commit_mismatch", "parent commit changed")
    if source.get("parent_head") is not None:
        _require(source.get("parent_head") == EXPECTED_PARENT_COMMIT, "parent_head_mismatch", "parent head changed")
    report = source.get("parent_report")
    if report is not None:
        _require(_is_mapping(report), "parent_report_shape", "parent report must be a mapping")
        _require(report.get("sha256") == EXPECTED_PARENT_REPORT_SHA256, "parent_report_pin_mismatch", "parent report SHA changed")
    return literals, {
        "input_schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "source_wrapper_names": list(EXPECTED_SOURCE_SHA256),
        "source_wrapper_pins": pins,
        "source_modules_imported": False,
        "project_modules_imported": False,
        "archive_opened_or_decoded": False,
        "source_value_evaluation": False,
        "numeric_evaluation": False,
        "N_materialized": False,
        "cached_algebra_replay_calls": 0,
    }


def _selected(rowset: Any, label: str) -> list[int]:
    _require(_is_mapping(rowset), "row_scope", f"{label} is missing")
    _require(rowset.get("selected_orders") == list(ORDINARY_Z_ORDERS), "row_scope", f"{label} selected orders changed")
    selected = rowset.get("selected_node_ids")
    _require(
        isinstance(selected, list)
        and len(selected) == 5
        and all(_is_int(node) and node >= 0 for node in selected),
        "row_scope",
        f"{label} selected node IDs changed",
    )
    _require(rowset.get("rows_are_genuine") is True, "row_scope", f"{label} rows are not genuine")
    _require(rowset.get("padding_rows_added") is False, "row_scope", f"{label} has padding rows")
    rows = rowset.get("rows")
    _require(isinstance(rows, list) and len(rows) >= 5, "row_scope", f"{label} rows are missing")
    for index, row in enumerate(rows[:5]):
        _require(_is_mapping(row), "row_scope", f"{label} row is not a mapping")
        _require(row.get("ordinary_Z_order") == index and row.get("node") == selected[index], "row_scope", f"{label} row port changed")
    return [int(node) for node in selected]


def _node(heads: Mapping[str, Any], ident: Any, label: str) -> Mapping[str, Any]:
    _require(_is_int(ident) and ident >= 0, "graph_node_shape", f"{label} node is not a nonnegative integer")
    record = heads.get(str(ident), heads.get(ident))
    _require(_is_mapping(record), "graph_node_missing", f"{label} graph node {ident} is missing")
    _require(record.get("node") == ident, "graph_node_shape", f"{label} graph node identity changed")
    return record


def _record_children(record: Mapping[str, Any], label: str) -> list[int]:
    children = record.get("ordered_children")
    _require(
        children is None or (isinstance(children, list) and all(_is_int(child) and child >= 0 for child in children)),
        "graph_children_shape",
        f"{label} ordered children changed",
    )
    return list(children or [])


def _source_children(
    record: Mapping[str, Any], constructor: Mapping[str, Any], ident: int, label: str
) -> list[int]:
    """Read direct children from the selected head or its constructor receipt."""
    children = record.get("ordered_children")
    if children is None:
        ordered = constructor.get("ordered_child_lists")
        _require(_is_mapping(ordered), "graph_children_shape", "constructor child receipt is missing")
        children = ordered.get(str(ident), ordered.get(ident))
    _require(
        children is None or (isinstance(children, list) and all(_is_int(child) and child >= 0 for child in children)),
        "graph_children_shape",
        f"{label} ordered children changed",
    )
    return list(children or [])


def _row_id_map(rows: Mapping[str, Any], label: str) -> dict[str, list[int]]:
    _require(_is_mapping(rows), "row_map_shape", f"{label} is missing")
    return {str(key): _selected(value, f"{label}.{key}") for key, value in rows.items()}


def _validate_ast_and_graph(raw: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    _require(raw.get("schema") == "native_O3_density_flux_raw_literals_v1", "raw_schema", "density-flux packet schema changed")
    _require(raw.get("packet_revision") == 2, "raw_revision", "the v2 density-flux packet is required")
    _require(raw.get("status") == "static_literal_extraction_only", "raw_status", "density-flux packet status changed")
    _require(raw.get("packet_type") == "native_O3_density_flux_raw_literals", "raw_packet_type", "density-flux packet type changed")
    scope = raw.get("scope")
    _require(_is_mapping(scope), "raw_scope", "density-flux packet scope missing")
    _require(scope.get("chart") == EXPECTED_CHART, "raw_scope", "density-flux chart changed")
    _require(scope.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS), "raw_scope", "density-flux order scope changed")
    _require(scope.get("selected_rows_are_genuine") is True, "raw_scope", "selected rows are not genuine")
    _require(scope.get("padding_rows_added") is False, "raw_scope", "padding rows were added")

    audit = raw.get("audit")
    _require(_is_mapping(audit), "raw_audit", "density-flux audit missing")
    for key in (
        "J_arithmetic_performed",
        "N_materialized",
        "algebra_called",
        "calculus_executed",
        "checkers_called",
        "derivative_polynomials_constructed",
        "differentiators_called",
        "flux_graph_constructed",
        "graph_reduction_or_normalization",
        "math_identity_functions_called",
        "project_consumers_called",
        "project_helpers_imported",
        "quadrature_or_recurrence_executed",
        "source_modules_imported",
        "source_proof_executed",
        "source_value_evaluation",
        "forbidden_continuation_decoded",
        "forbidden_huge_Rm_loop_density_decoded",
    ):
        _require(audit.get(key) in (False, 0), "raw_audit", f"forbidden packet action was recorded: {key}")
    _require(audit.get("identity_checks") == 0, "raw_audit", "packet identity checks were executed")
    _require(audit.get("raw_archive_authenticated_before_decode") is True, "raw_audit", "archive authentication receipt changed")

    literal_rows = raw.get("literal_rows")
    _require(_is_mapping(literal_rows), "literal_rows_shape", "literal rows missing")
    maps: dict[str, dict[str, list[int]]] = {}
    selected_total = 0
    for group, rows in literal_rows.items():
        maps[group] = _row_id_map(rows, f"literal_rows.{group}")
        selected_total += sum(len(ids) for ids in maps[group].values())
    _require(selected_total == 120, "literal_rows_scope", "selected literal row count changed")
    required_groups = {
        "background_fields",
        "chart_roots",
        "finite_N",
        "finite_N_constructor_rows",
        "five_signed_density_increments",
        "loop",
        "modified_fields",
        "saved_total_radial",
    }
    _require(set(maps) == required_groups, "literal_rows_scope", "literal row groups changed")

    ast_data = raw.get("source_ast")
    _require(_is_mapping(ast_data), "source_ast_shape", "source AST witnesses missing")
    expected_ast_keys = {
        "BridgeAlgebra.exp",
        "JetAlgebra.add",
        "JetAlgebra.mul",
        "JetAlgebra.scale",
        "background_density",
        "finite_N_functions",
    }
    _require(set(ast_data) == expected_ast_keys, "source_ast_shape", "source AST witness keys changed")
    for name in expected_ast_keys:
        _require(_is_mapping(ast_data[name]), "source_ast_shape", f"source AST witness changed: {name}")
    finite_text = str(ast_data["finite_N_functions"].get("source_text", ""))
    for marker in (
        "ex=a.exp(x)",
        "expm1=FunctionJet",
        "r=a.mul(roots['E'],expm1)",
        "s=a.scale(loop['B_over_Pstar'],invN)",
        "E_N=a.mul(E,ex)",
        "V_N=a.add(V,s)",
        "five_signed_density_increments=inc",
    ):
        _require(marker in finite_text, "source_ast_shape", f"finite_N_functions witness lacks {marker}")
    density_text = str(ast_data["background_density"].get("source_text", ""))
    for marker in ("dict(m=V", "a.mul(E,V)", "a.scale(EE,'1/2')"):
        _require(marker in density_text, "source_ast_shape", f"background_density witness lacks {marker}")
    for member in ("add", "mul", "scale"):
        text = str(ast_data[f"JetAlgebra.{member}"].get("source_text", ""))
        _require("math.comb" in text or member in ("add", "scale"), "source_ast_shape", f"JetAlgebra.{member} witness changed")
    exp_text = str(ast_data["BridgeAlgebra.exp"].get("source_text", ""))
    _require("unary('exp'" in exp_text and "math.comb" in exp_text, "source_ast_shape", "BridgeAlgebra.exp witness changed")

    graph = raw.get("source_graph")
    _require(_is_mapping(graph), "source_graph_shape", "source graph missing")
    _require(graph.get("chart") == EXPECTED_CHART, "source_graph_shape", "source graph chart changed")
    _require(graph.get("head_count") == 84, "source_graph_shape", "source graph head count changed")
    constructor = graph.get("constructor_skeleton")
    _require(_is_mapping(constructor), "source_graph_shape", "constructor skeleton missing")
    _require(constructor.get("root_count") == 50, "source_graph_shape", "constructor root count changed")
    _require(constructor.get("unique_record_count") == 189, "source_graph_shape", "constructor record count changed")
    _require(constructor.get("no_derivative_polynomials_constructed") is True, "source_graph_shape", "derivative polynomial construction changed")
    _require(constructor.get("no_flux_graph_constructed") is True, "source_graph_shape", "flux graph construction changed")
    heads = graph.get("heads")
    _require(_is_mapping(heads), "source_graph_shape", "source graph heads missing")
    supplement = raw.get("static_child_record_supplement")
    _require(
        _is_mapping(supplement)
        and supplement.get("new_identity_checks") == 0
        and supplement.get("raw_nodes_added") == [249091, 249092, 249093, 249094, 249095],
        "source_graph_shape",
        "v2 negative-child supplement changed",
    )
    zero = graph.get("zero_node_provenance")
    _require(_is_mapping(zero) and zero.get("node") == 0 and zero.get("genuine_background_zero") is True and zero.get("padding") is False, "source_graph_shape", "zero-node provenance changed")

    opaque = raw.get("opaque_source_parameters")
    _require(_is_mapping(opaque), "opaque_source_shape", "opaque source parameters missing")
    n_data, inv_data = opaque.get("N"), opaque.get("inverse_N_graph_port")
    _require(_is_mapping(n_data) and n_data.get("node") == N_NODE and n_data.get("materialized") is False and n_data.get("arithmetic_performed") is False, "opaque_source_shape", "opaque N contract changed")
    _require(_is_mapping(inv_data) and inv_data.get("node") == INVERSE_N_NODE and inv_data.get("materialized") is False and inv_data.get("arithmetic_performed") is False, "opaque_source_shape", "inverse-N contract changed")

    # Authenticate only the selected source ports and their direct constructor relationships.
    background = maps["background_fields"]
    chart = maps["chart_roots"]
    finite = maps["finite_N"]
    constructor_rows = maps["finite_N_constructor_rows"]
    increments = maps["five_signed_density_increments"]
    loop = maps["loop"]
    modified = maps["modified_fields"]
    saved_radial = maps["saved_total_radial"]
    _require(background["E"] == chart["E"], "source_port_shape", "background/chart E ports changed")
    _require(background["V"] == chart["V"], "source_port_shape", "background/chart V ports changed")
    _require(finite["A_over_N"] == [244325, 244326, 244327, 244328, 244329], "source_port_shape", "A/N source ports changed")
    _require(constructor_rows["exp_A_over_N"] == [244331, 244332, 244335, 244339, 244344], "source_port_shape", "exp source ports changed")
    _require(constructor_rows["expm1_A_over_N"] == [244351, 244332, 244335, 244339, 244344], "source_port_shape", "expm1 source ports changed")
    _require(opaque.get("N", {}).get("node") == N_NODE, "source_port_shape", "N node changed")
    for index, node in enumerate(finite["A_over_N"]):
        record = _node(heads, node, f"A_over_N[{index}]")
        _require(record.get("operation") == "product", "source_port_shape", "A/N is not an actual product source port")
        children = _source_children(record, constructor, node, f"A_over_N[{index}]")
        _require(children == [loop["A"][index], INVERSE_N_NODE], "source_port_shape", "A/N constructor relationship changed")
    for index, node in enumerate(constructor_rows["exp_A_over_N"]):
        record = _node(heads, node, f"exp_A_over_N[{index}]")
        expm_record = _node(heads, constructor_rows["expm1_A_over_N"][index], f"expm1_A_over_N[{index}]")
        if index == 0:
            _require(record.get("operation") == "analytic_unary" and record.get("raw_record", {}).get("name") == "exp", "source_port_shape", "exp value creator changed")
            _require(expm_record.get("operation") == "analytic_unary" and expm_record.get("raw_record", {}).get("name") == "expm1", "source_port_shape", "expm1 value creator changed")
            _require(_source_children(record, constructor, node, f"exp_A_over_N[{index}]") == [finite["A_over_N"][0]], "source_port_shape", "exp argument changed")
            _require(_source_children(expm_record, constructor, constructor_rows["expm1_A_over_N"][index], f"expm1_A_over_N[{index}]") == [finite["A_over_N"][0]], "source_port_shape", "expm1 argument changed")
        else:
            _require(constructor_rows["expm1_A_over_N"][index] == node, "source_port_shape", "higher exp/expm1 source rows are not shared")
            _require(record.get("operation") in ("product", "sum"), "source_port_shape", "higher exp source row operation changed")
    for key in ("E_N", "V_N", "delta_E", "delta_V"):
        for index, node in enumerate(finite[key]):
            record = _node(heads, node, f"finite_N.{key}[{index}]")
            _require(record.get("operation") in ("product", "sum"), "source_port_shape", f"finite_N.{key} source row operation changed")
    for key, ids in increments.items():
        for index, node in enumerate(ids):
            record = _node(heads, node, f"increment.{key}[{index}]")
            _require(record.get("operation") in ("product", "sum"), "source_port_shape", f"increment.{key} source row operation changed")
    for index, node in enumerate(saved_radial["m_y"]):
        record = _node(heads, node, f"saved_total_radial.m_y[{index}]")
        children = _record_children(record, f"saved_total_radial.m_y[{index}]")
        _require(record.get("operation") == "sum" and len(children) == 2, "source_port_shape", "saved m_y is not a strict binary sum")
        _require(children[0] == finite["V_N"][index] and children[1] != children[0], "source_port_shape", "saved m_y child relationship changed")
    aliases = graph.get("aliases_from_actual_handles")
    _require(isinstance(aliases, list) and len(aliases) == 5, "source_alias_shape", "actual source alias receipts changed")
    _require(all(_is_mapping(alias) and alias.get("identity_test_performed") is False for alias in aliases), "source_alias_shape", "old alias identity checks were executed")

    return {
        "raw_packet_schema": raw.get("schema"),
        "raw_status": raw.get("status"),
        "chart": EXPECTED_CHART,
        "selected_literal_rows": selected_total,
        "source_graph_head_count": graph.get("head_count"),
        "source_graph_constructor_record_count": constructor.get("unique_record_count"),
        "selected_source_ports_only": True,
        "source_ast_witness_names": {
            "BridgeAlgebra.exp": ["exp"],
            "JetAlgebra.add": ["add"],
            "JetAlgebra.mul": ["mul"],
            "JetAlgebra.scale": ["scale"],
            "background_density": ["background_density"],
            "finite_N_functions": ["finite_N_functions"],
        },
        "handles": {
            "background_E": background["E"],
            "background_V": background["V"],
            "chart_E": chart["E"],
            "chart_V": chart["V"],
            "A_over_N": finite["A_over_N"],
            "E_N": finite["E_N"],
            "V_N": finite["V_N"],
            "delta_E": finite["delta_E"],
            "delta_V": finite["delta_V"],
            "exp_A_over_N": constructor_rows["exp_A_over_N"],
            "expm1_A_over_N": constructor_rows["expm1_A_over_N"],
            "increments": increments,
            "modified_fields": modified,
            "saved_m_y": saved_radial["m_y"],
            "loop": loop,
            "N": N_NODE,
            "inverse_N": INVERSE_N_NODE,
        },
        "ordinary_C4_source_jet_semantics": True,
        "global_ordinary_jet_oracle": False,
        "source_modules_executed": False,
        "archive_decoded": False,
        "N_materialized": False,
        "numeric_evaluation": False,
        "old_algebra_replayed": False,
        "cached_generic_proof_reused_as_data": True,
    }, {"heads": heads, "literal_rows": literal_rows, "ast": ast_data}


def _validate_generic_scope(generic: Mapping[str, Any]) -> dict[str, Any]:
    _require(generic.get("schema") == "native_O3_generic_density_certificate_scope_v1", "generic_schema", "generic density scope schema changed")
    _require(generic.get("status") == "static_bounded_audit", "generic_status", "generic density scope status changed")
    decision = generic.get("decision")
    _require(_is_mapping(decision), "generic_scope", "generic density decision missing")
    for key in ("abstract_polynomial_schema_exists", "abstract_polynomial_schema_is_exact", "accepted_generic_identity_witness", "new_source_binding_leaf_required"):
        _require(decision.get(key) is True, "generic_scope", f"generic density decision changed: {key}")
    _require(decision.get("accepted_universal_ordinary_jet_theorem_through_Z4") is False, "generic_scope", "old universal jet theorem was accepted")
    abstract = generic.get("abstract_density_schema")
    _require(_is_mapping(abstract), "generic_scope", "abstract density schema missing")
    expected_base = {"m": "V", "h": "E", "k": "E*V", "e": "V^2-E^2/2", "p": "E^2/2"}
    expected_inc = {"m": "s", "h": "r", "k": "V*r+E*s+r*s", "e": "2*V*s+s^2-E*r-r^2/2", "p": "E*r+r^2/2"}
    _require(abstract.get("base_density_functions") == expected_base, "generic_scope", "abstract base density schema changed")
    _require(abstract.get("flux_or_density_increment") == expected_inc, "generic_scope", "abstract increment schema changed")
    _require(abstract.get("field_alias_contract") == {"E_N": "E+r", "V_N": "V+s"}, "generic_scope", "abstract alias contract changed")
    accepted = generic.get("accepted_generic_proof")
    _require(_is_mapping(accepted), "generic_scope", "accepted generic proof data missing")
    _require(_is_mapping(accepted.get("exact_identity")), "generic_scope", "generic exact identity data missing")
    _require(_is_mapping(accepted.get("identity_check_witness")), "generic_scope", "generic identity witness missing")
    _require(accepted["identity_check_witness"].get("ordinary_Z_orders") == [0, 1, 2], "generic_scope", "generic proof order scope changed")
    _require(generic.get("static_policy", {}).get("source_imports_or_calls") == 0, "generic_scope", "generic provider was called")
    _require(generic.get("static_policy", {}).get("math_reexecuted") is False, "generic_scope", "generic proof was reexecuted")
    _require(generic.get("static_policy", {}).get("numeric_values_computed") is False, "generic_scope", "generic values were computed")
    return {
        "schema": generic.get("schema"),
        "status": generic.get("status"),
        "abstract_density_schema": _clone(abstract),
        "accepted_generic_proof": _clone(accepted),
        "accepted_orders": [0, 1, 2],
        "ordinary_Z3_Z4_fresh_extension_required": True,
        "O3_source_alias_E_plus_r_equals_E_N": False,
        "O3_flux_E_N_V_N_certificate": False,
        "consumed_as_data_only": True,
        "recomputed": False,
        "numeric_evaluation": False,
    }


def _parent_payload(saved: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    _require(saved.get("schema") is not None, "parent_shape", "saved history result schema missing")
    result = saved.get("result")
    _require(_is_mapping(result), "parent_shape", "saved history nested result missing")
    _require(result.get("passed") is True, "parent_shape", "saved history result is not passed")
    _require(result.get("source_only") is True and result.get("full_goal_complete") is False, "parent_shape", "saved history scope changed")
    _require(result.get("conditional") is True, "parent_shape", "saved history conditional scope changed")
    payload = result.get("payload")
    _require(_is_mapping(payload), "parent_shape", "saved history payload missing")
    duhamel = payload.get("duhamel_source_rows")
    modified = payload.get("modified_history_rows")
    _require(isinstance(duhamel, list) and len(duhamel) == 50, "parent_rows", "saved Duhamel row count changed")
    _require(isinstance(modified, list) and len(modified) == 25, "parent_rows", "saved modified history row count changed")
    bindings_key = "bindings" if _is_mapping(payload.get("bindings")) else "physical_radius_bindings"
    bindings = payload.get(bindings_key)
    _require(_is_mapping(bindings), "parent_bindings", "saved physical binding map missing")
    for branch in BRANCHES:
        branch_rows = bindings.get(branch)
        _require(_is_mapping(branch_rows), "parent_bindings", f"saved branch missing: {branch}")
        _require(set(branch_rows) == set(FIELDS), "parent_bindings", f"saved field set changed: {branch}")
        for field in FIELDS:
            field_rows = branch_rows.get(field)
            _require(_is_mapping(field_rows), "parent_bindings", f"saved field rows missing: {branch}.{field}")
            for order in ORDINARY_Z_ORDERS:
                row = field_rows.get(str(order), field_rows.get(order))
                _require(_is_mapping(row), "parent_bindings", f"saved binding missing: {branch}.{field}.{order}")
                _require(row.get("branch") == branch and row.get("field") == field and row.get("ordinary_Z_order", order) == order, "parent_bindings", "saved binding identity changed")
                _require(row.get("radial_order") == 1 and row.get("phase_mode") == "ordinary_Z", "parent_bindings", "saved binding scope changed")
    by_duhamel: dict[tuple[str, str, int], Mapping[str, Any]] = {}
    for row in duhamel:
        _require(_is_mapping(row), "parent_rows", "Duhamel row is not a mapping")
        branch, history, order = row.get("branch"), row.get("history_key"), row.get("ordinary_Z_order")
        _require(branch in ("background", "correction") and history in HISTORY_KEYS and _is_int(order) and order in ORDINARY_Z_ORDERS, "parent_rows", "Duhamel row identity changed")
        _require((branch, history, order) not in by_duhamel, "parent_rows", "duplicate Duhamel row")
        _require(_is_int(row.get("density_source_node")) and _is_int(row.get("outgoing_node")), "parent_rows", "Duhamel source node missing")
        _require(_is_mapping(row.get("source_expression")) and _is_mapping(row.get("raw_node_bindings")), "parent_rows", "Duhamel source expression missing")
        _require(row.get("source_relationship_checked") is True, "parent_rows", "Duhamel source relationship was not checked")
        by_duhamel[(branch, history, order)] = row
    _require(len(by_duhamel) == 50, "parent_rows", "Duhamel row keys incomplete")
    by_modified: dict[tuple[str, int], Mapping[str, Any]] = {}
    for row in modified:
        _require(_is_mapping(row), "parent_rows", "modified history row is not a mapping")
        history, order = row.get("history_key"), row.get("ordinary_Z_order")
        _require(row.get("branch") == "modified" and history in HISTORY_KEYS and _is_int(order) and order in ORDINARY_Z_ORDERS, "parent_rows", "modified history row identity changed")
        _require((history, order) not in by_modified, "parent_rows", "duplicate modified history row")
        _require(row.get("raw_graph_sum_verified") is True, "parent_rows", "modified history sum receipt changed")
        premises = row.get("premises")
        _require(_is_mapping(premises) and premises.get("nonzero_incoming_histories_retained") is True, "parent_rows", "nonzero incoming history premise changed")
        source = row.get("source_expression")
        _require(_is_mapping(source) and source.get("equation") == "background_density_j(x)+correction_density_j(x)-lambda*modified_j(x)", "parent_rows", "modified Duhamel source equation changed")
        _require(source.get("rate_lambda") == RATES[history], "parent_rows", "modified Duhamel rate changed")
        terms = source.get("terms")
        _require(isinstance(terms, list) and len(terms) == 3, "parent_rows", "modified Duhamel term count changed")
        bg = by_duhamel[("background", history, order)]
        corr = by_duhamel[("correction", history, order)]
        _require(terms[0].get("node") == bg.get("density_source_node"), "parent_rows", "background density source node changed")
        _require(terms[1].get("node") == corr.get("density_source_node"), "parent_rows", "correction density source node changed")
        _require(terms[2].get("node") == row.get("modified_node"), "parent_rows", "modified history source node changed")
        by_modified[(history, order)] = row
    _require(len(by_modified) == 25, "parent_rows", "modified history keys incomplete")
    cached = payload.get("cached_data")
    _require(_is_mapping(cached), "parent_cache", "cached radial frame data missing")
    native_r = cached.get("nativeR_certificate")
    _require(_is_mapping(native_r) and native_r.get("root") == RADIUS_NODE and native_r.get("Z_excluded") is True, "parent_cache", "native R certificate changed")
    frames = cached.get("frame_specialization_rows")
    _require(isinstance(frames, list) and len(frames) == 2, "parent_cache", "cached native R frame rows changed")
    for index, branch in enumerate(BRANCHES):
        frame = frames[index]
        _require(_is_mapping(frame) and frame.get("branch") == branch and frame.get("P_Z") == "D_Z" and frame.get("P_Z_identity") is True, "parent_cache", "cached frame specialization changed")
        manifest = frame.get("branch_operator_manifest")
        _require(_is_mapping(manifest), "parent_cache", "cached branch frame manifest missing")
        frame_data = manifest.get("frame")
        _require(_is_mapping(frame_data) and frame_data.get("radius_node") == RADIUS_NODE, "parent_cache", "cached frame radius changed")
    return {
        "payload": payload,
        "bindings_key": bindings_key,
        "bindings": bindings,
        "duhamel": by_duhamel,
        "modified": by_modified,
        "cached": cached,
    }, {
        "saved_result_schema": saved.get("schema"),
        "nested_result_schema": result.get("schema"),
        "payload_bindings_key": bindings_key,
        "duhamel_rows": 50,
        "modified_history_rows": 25,
        "physical_binding_rows": 90,
        "native_radius_node": RADIUS_NODE,
        "cached_native_R_retained": True,
        "source_unit": "native_R",
        "physical_time_map_instantiated": False,
    }


def _graph_ref(name: str, node: int, *, source_path: str, ordinary_Z_order: int | None = None) -> dict[str, Any]:
    value: dict[str, Any] = {"kind": "native_graph_ref", "name": name, "node": node, "source_path": source_path}
    if ordinary_Z_order is not None:
        value["ordinary_Z_order"] = ordinary_Z_order
    return value


def _rational(numerator: int, denominator: int = 1, *, source: str | None = None) -> dict[str, Any]:
    value: dict[str, Any] = {"kind": "rational", "numerator": numerator, "denominator": denominator}
    if source is not None:
        value["source"] = source
    return value


def _product(*arguments: Mapping[str, Any]) -> dict[str, Any]:
    return {"kind": "product", "arguments": [_clone(argument) for argument in arguments]}


def _sum(*arguments: Mapping[str, Any]) -> Mapping[str, Any]:
    return {"kind": "sum", "arguments": [_clone(argument) for argument in arguments]}


def _field_ref(handles: Mapping[str, Any], field: str, order: int, *, source_path: str) -> dict[str, Any]:
    return _graph_ref(field, int(handles[field][order]), source_path=source_path, ordinary_Z_order=order)


def _validate_expression_kinds(value: Any) -> None:
    """Reject accidental evaluator/normalizer nodes in fresh source DAGs."""
    if isinstance(value, Mapping):
        if "kind" in value:
            _require(value.get("kind") in {"native_graph_ref", "rational", "product", "sum", "typed_quotient"}, "expression_shape", "fresh expression contains an unapproved node kind")
        for child in value.values():
            _validate_expression_kinds(child)
    elif isinstance(value, list):
        for child in value:
            _validate_expression_kinds(child)


def _alias_row(handles: Mapping[str, Any], order: int) -> dict[str, Any]:
    e_node = int(handles["background_E"][order])
    r_node = int(handles["delta_E"][order])
    en_node = int(handles["E_N"][order])
    x_node = int(handles["A_over_N"][order])
    exp_node = int(handles["exp_A_over_N"][order])
    expm1_node = int(handles["expm1_A_over_N"][order])
    e_ref = _graph_ref("E", e_node, source_path=f"literal_rows.background_fields.E[{order}]", ordinary_Z_order=order)
    r_ref = _graph_ref("r", r_node, source_path=f"literal_rows.finite_N.delta_E[{order}]", ordinary_Z_order=order)
    en_ref = _graph_ref("E_N", en_node, source_path=f"literal_rows.finite_N.E_N[{order}]", ordinary_Z_order=order)
    x_ref = _graph_ref("A_over_N", x_node, source_path=f"literal_rows.finite_N.A_over_N[{order}]", ordinary_Z_order=order)
    exp_ref = _graph_ref("exp(A_over_N)", exp_node, source_path=f"literal_rows.finite_N_constructor_rows.exp_A_over_N[{order}]", ordinary_Z_order=order)
    expm1_ref = _graph_ref("expm1(A_over_N)", expm1_node, source_path=f"literal_rows.finite_N_constructor_rows.expm1_A_over_N[{order}]", ordinary_Z_order=order)
    rhs = _sum(e_ref, r_ref)
    _validate_expression_kinds(rhs)
    return {
        "ordinary_Z_order": order,
        "source_inputs": {"E": e_ref, "r": r_ref, "A_over_N": x_ref, "exp": exp_ref, "expm1": expm1_ref, "E_N": en_ref},
        "analytic_value_law": {
            "equation": "expm1(a)=exp(a)-1",
            "valid_at_a_zero": True,
            "division_by_A": False,
            "value_zero_creator": "analytic_expm1",
            "value_zero_source_refs": {"exp": exp_ref, "expm1": expm1_ref, "argument": x_ref},
            "blind_exp_expm1_node_alias": False,
        },
        "higher_row_source_law": {
            "same_exp_derivative_source_ref": order > 0,
            "expm1_source_node_is_exp_source_node": order > 0,
            "source_factory_rule": "value row uses analytic_expm1; higher ordinary rows use full exp derivative source",
        },
        "functional_identity": {
            "equation": "E_N=E+r",
            "lhs": en_ref,
            "rhs": rhs,
            "source_factory": "finite_N_functions",
            "same_E_and_A_over_N_source_namespace": True,
            "constant_N_node": N_NODE,
            "inverse_N_node": INVERSE_N_NODE,
            "A_zero_branch_included": True,
            "functional_not_numeric": True,
        },
        "source_expression": rhs,
        "ordinary_jet_lift": {
            "source_jet_semantics": "genuine ordinary C4 source jet",
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "current_order": order,
            "Leibniz_binomial_factorial_semantics": True,
            "binomial_table": [list(row) for row in BINOMIALS],
            "constant_N_rows_through_order4": True,
            "global_jet_oracle": False,
            "phase_mode": "ordinary_Z",
            "phase_seam_regularity_proved": False,
            "incoming_x_independence_proved": False,
        },
        "evaluated": False,
        "source_identity_status": "fresh O3 functional alias under conditional ordinary C4 source-jet semantics",
    }


def _flux_expression(handles: Mapping[str, Any], field: str, order: int) -> Mapping[str, Any]:
    en = lambda index: _field_ref(handles, "E_N", index, source_path=f"literal_rows.finite_N.E_N[{index}]")
    vn = lambda index: _field_ref(handles, "V_N", index, source_path=f"literal_rows.finite_N.V_N[{index}]")
    if field == "m":
        expression = vn(order)
    elif field == "h":
        expression = en(order)
    elif field == "k":
        terms = [_product(_rational(BINOMIALS[order][index], source="ordinary_Leibniz_binomial"), en(index), vn(order - index)) for index in range(order + 1)]
        expression = terms[0] if len(terms) == 1 else _sum(*terms)
    elif field == "e":
        vv = [_product(_rational(BINOMIALS[order][index], source="ordinary_Leibniz_binomial"), vn(index), vn(order - index)) for index in range(order + 1)]
        ee = [_product(_rational(BINOMIALS[order][index], source="ordinary_Leibniz_binomial"), en(index), en(order - index)) for index in range(order + 1)]
        vv_expr = vv[0] if len(vv) == 1 else _sum(*vv)
        ee_expr = ee[0] if len(ee) == 1 else _sum(*ee)
        expression = _sum(vv_expr, _product(_rational(-1, 2, source="base_density_minus_one_half"), ee_expr))
    elif field == "p":
        ee = [_product(_rational(BINOMIALS[order][index], source="ordinary_Leibniz_binomial"), en(index), en(order - index)) for index in range(order + 1)]
        ee_expr = ee[0] if len(ee) == 1 else _sum(*ee)
        expression = _product(_rational(1, 2, source="base_density_one_half"), ee_expr)
    else:
        raise _KernelFailure("flux_field", f"unsupported density flux field: {field}")
    _validate_expression_kinds(expression)
    return expression


def _density_row(
    raw_handles: Mapping[str, Any],
    generic: Mapping[str, Any],
    duhamel: Mapping[tuple[str, str, int], Mapping[str, Any]],
    history: str,
    order: int,
) -> dict[str, Any]:
    bg = duhamel[("background", history, order)]
    increment_node = int(raw_handles["increments"][history][order])
    bg_ref = _graph_ref(
        f"background_density_{history}",
        int(bg["density_source_node"]),
        source_path=f"saved_history_result.result.payload.duhamel_source_rows[background,{history},{order}].density_source_node",
        ordinary_Z_order=order,
    )
    inc_ref = _graph_ref(
        f"finite_N_increment_{history}",
        increment_node,
        source_path=f"literal_rows.five_signed_density_increments.{history}[{order}]",
        ordinary_Z_order=order,
    )
    base_density = _sum(bg_ref, inc_ref)
    flux = _flux_expression(raw_handles, history, order)
    _validate_expression_kinds(base_density)
    return {
        "history_key": history,
        "field": history,
        "ordinary_Z_order": order,
        "base_density": base_density,
        "source_density_identity": {
            "equation": "saved_background_density_j+finite_N_increment_j=D_Z^j F(E_N,V_N)",
            "background_density_source": bg_ref,
            "finite_N_increment_source": inc_ref,
            "pointwise_flux": flux,
            "generic_schema_data_only": True,
            "abstract_polynomial_schema": _clone(generic["abstract_density_schema"]),
        },
        "flux_expression": flux,
        "source_expression": flux,
        "ordinary_jet_semantics": {
            "source_jet_semantics": "genuine ordinary C4 source jet",
            "ordinary_Z_order": order,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "binomial_table": [list(row) for row in BINOMIALS],
            "Leibniz_rule_applied": True,
            "factorial_normalization_preserved": True,
            "constant_N_node": N_NODE,
            "N_materialized": False,
            "global_jet_oracle": False,
            "phase_seam_regularity_proved": False,
        },
        "evaluated": False,
        "source_identity_status": "fresh O3 density-flux source identity",
    }


def _parent_binding(parent_bindings: Mapping[str, Any], branch: str, field: str, order: int) -> Mapping[str, Any]:
    branch_rows = parent_bindings.get(branch)
    field_rows = branch_rows.get(field) if _is_mapping(branch_rows) else None
    if not _is_mapping(field_rows):
        raise _KernelFailure("parent_bindings", f"parent binding map missing {branch}.{field}")
    row = field_rows.get(str(order), field_rows.get(order))
    if not _is_mapping(row):
        raise _KernelFailure("parent_bindings", f"parent binding map missing {branch}.{field}.{order}")
    return row


def _rewrite_row(
    parent_row: Mapping[str, Any],
    density: Mapping[str, Any],
    modified: Mapping[str, Any],
    history: str,
    order: int,
) -> dict[str, Any]:
    lambda_num, lambda_den = RATE_RATIONALS[history]
    history_ref = _graph_ref(
        f"modified_{history}_history",
        int(modified["modified_node"]),
        source_path=f"saved_history_result.result.payload.modified_history_rows[{history},{order}].modified_node",
        ordinary_Z_order=order,
    )
    rate_term = _product(_rational(lambda_num, lambda_den, source="negative_Duhamel_rate"), history_ref)
    numerator = _sum(density["flux_expression"], rate_term)
    _validate_expression_kinds(numerator)
    prior_expression = parent_row.get("source_expression")
    return {
        "history_key": history,
        "ordinary_Z_order": order,
        "prior_conditional_source_row": _clone(parent_row),
        "prior_source_expression": _clone(prior_expression),
        "fresh_flux_expression": _clone(density["flux_expression"]),
        "same_actual_modified_history_source": history_ref,
        "rewritten_source_expression": {
            "equation": "F(E_N,V_N)_j-lambda*H_j",
            "source_operator": "D_y",
            "rate_lambda": RATES[history],
            "numerator": numerator,
            "evaluated": False,
        },
        "source_expression": {
            "equation": "F(E_N,V_N)_j-lambda*H_j",
            "source_operator": "D_y",
            "rate_lambda": RATES[history],
            "numerator": _clone(numerator),
            "evaluated": False,
        },
        "premises": {
            "background_and_correction_source_jets_are_ordinary_Z_derivatives": True,
            "genuine_ordinary_C4_source_jet": True,
            "nonzero_incoming_histories_retained": True,
            "all_parent_regular_interior_phase_cell_premises_retained": True,
            "phase_seam_regularity_proved": False,
            "incoming_x_independence_proved": False,
            "numeric_admission": False,
            "history_reset": False,
            "factor_Pstar_or_N_applied": False,
        },
        "evaluated": False,
        "source_identity_status": "fresh flux rewrite of current conditional Duhamel source row",
    }


def _saved_m_y_row(
    handles: Mapping[str, Any],
    heads: Mapping[str, Any],
    parent_bindings: Mapping[str, Any],
    modified: Mapping[tuple[str, int], Mapping[str, Any]],
    order: int,
) -> dict[str, Any]:
    node = int(handles["saved_m_y"][order])
    head = _node(heads, node, f"saved_total_radial.m_y[{order}]")
    children = _record_children(head, f"saved_total_radial.m_y[{order}]")
    vn_node = int(handles["V_N"][order])
    _require(head.get("operation") == "sum" and len(children) == 2, "saved_m_y_shape", "saved m_y source is not a strict binary sum")
    _require(children[0] == vn_node and children[1] != vn_node, "saved_m_y_shape", "saved m_y source children changed")
    parent_m = _parent_binding(parent_bindings, "modified", "m", order)
    modified_m_node = int(modified[("m", order)]["modified_node"])
    negative = _node(heads, children[1], f"saved_total_radial.m_y[{order}].second_child")
    negative_children = _record_children(negative, f"saved_total_radial.m_y[{order}].second_child")
    negative_raw = negative.get("raw_record")
    _require(
        negative.get("operation") == "negative"
        and negative_children == [modified_m_node]
        and _is_mapping(negative_raw)
        and negative_raw.get("operation") == "negative"
        and negative_raw.get("argument") == modified_m_node,
        "saved_m_y_shape",
        "saved m_y negative child record does not bind the actual modified m node",
    )
    m_ref = _graph_ref("modified_m", modified_m_node, source_path=f"saved_history_result.result.payload.modified_history_rows[m,{order}].modified_node", ordinary_Z_order=order)
    vn_ref = _graph_ref("V_N", vn_node, source_path=f"literal_rows.finite_N.V_N[{order}]", ordinary_Z_order=order)
    fresh_alias = _sum(vn_ref, _product(_rational(-1, source="saved_m_y_negative_child"), m_ref))
    _validate_expression_kinds(fresh_alias)
    return {
        "ordinary_Z_order": order,
        "saved_m_y": _graph_ref("saved_m_y", node, source_path=f"literal_rows.saved_total_radial.m_y[{order}]", ordinary_Z_order=order),
        "saved_graph_structure": {
            "operation": "sum",
            "ordered_children": children,
            "strict_sum_verified": True,
            "positive_child": _graph_ref("V_N", children[0], source_path=f"literal_rows.finite_N.V_N[{order}]", ordinary_Z_order=order),
            "negative_child_node": children[1],
            "negative_child_record": {
                "node": children[1],
                "operation": negative.get("operation"),
                "ordered_children": negative_children,
                "raw_record_operation": negative_raw.get("operation"),
                "raw_record_argument": negative_raw.get("argument"),
                "raw_record_sha256": negative.get("raw_record_sha256"),
            },
            "negative_child_declared_by_saved_source": True,
            "actual_negative_child_checked": True,
            "child_operation_not_reconstructed": False,
        },
        "fresh_rewritten_flux": _clone(vn_ref),
        "fresh_alias_expression": fresh_alias,
        "source_expression": fresh_alias,
        "modified_history_source": m_ref,
        "parent_modified_binding_source": _clone(parent_m.get("source_expression")),
        "premises": {
            "strict_sum_and_negative_child_source_verified": True,
            "actual_negative_child_checked": True,
            "actual_saved_radial_source_used": True,
            "fresh_rewritten_mflux_is_V_N": True,
            "factor_Pstar_or_N_applied": False,
            "numeric_admission": False,
            "phase_seam_regularity_proved": False,
        },
        "evaluated": False,
    }


def _radius_ref() -> dict[str, Any]:
    return _graph_ref("R", RADIUS_NODE, source_path="saved_history_result.result.payload.cached_data.nativeR_certificate.root")


def _quotient_row(
    parent_row: Mapping[str, Any],
    field: str,
    order: int,
    numerator: Mapping[str, Any],
    *,
    pressure: bool = False,
    p0_node: int | None = None,
) -> dict[str, Any]:
    radius = _radius_ref()
    quotient = {
        "kind": "typed_quotient",
        "operator": "D_R",
        "equation": "D_R_source_numerator_divided_by_native_R",
        "numerator": _clone(numerator),
        "denominator": radius,
        "coefficient": _rational(1),
        "evaluated": False,
    }
    _validate_expression_kinds(quotient)
    row = _clone(dict(parent_row))
    row["source_expression"] = quotient
    row["derived_expression"] = {
        "field": field,
        "operator": "D_R",
        "ordinary_Z_order": order,
        "native_radius_node": RADIUS_NODE,
        "source_expression": _clone(quotient),
        "resolvable": True,
        "evaluated": False,
    }
    row["source_identity_status"] = "fresh O3 density-flux first-radial source binding"
    row["fresh_source_binding"] = True
    row["source_unit"] = "native_R"
    row["physical_time_map_instantiated"] = False
    row["parent_binding_premises_retained"] = _clone(parent_row.get("premises"))
    row["parent_binding_metadata_retained"] = True
    row["premises"] = {
        "native_radius_positive_source_contract": True,
        "cached_R_Z_used_as_data_only": True,
        "cached_native_R_frame_Z_invariance_retained": True,
        "cached_theorem_reproved": False,
        "genuine_ordinary_C4_source_jet": True,
        "actual_field_derivative_identity_proved": False,
        "phase_seam_regularity_proved": False,
        "numeric_admission": False,
    }
    if pressure:
        row["pressure_source_contract"] = {
            "p_derivative_expression": _clone(numerator.get("arguments", [numerator])[-1] if _is_mapping(numerator) else numerator),
            "cached_declared_P0_radial_zero": True,
            "P0_source_node": p0_node,
            "P0_ordinary_Z_value_claimed": False,
            "ordinary_Z0_claimed": False,
        }
    return row


def _api_refusal(*, branch: Any, field: Any, order: Any, radial_order: Any, phase_mode: Any, reason: str, control: str) -> _APIRefusal:
    raise _APIRefusal(
        {
            "decision": "refused",
            "branch": branch,
            "field": field,
            "j": order,
            "ordinary_Z_order": order,
            "radial_order": radial_order,
            "phase_mode": phase_mode,
            "reason": reason,
            "control": control,
            "typed_source_refusal": True,
            "source_value_evaluated": False,
            "numeric_call": False,
        }
    )


def lookup(
    bindings: Mapping[str, Any],
    branch: Any = None,
    field: Any = None,
    j: Any = None,
    radial_order: Any = 1,
    phase_mode: Any = "ordinary_Z",
    *,
    query: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    """Look up one public first-radial ordinary-Z source binding."""
    if query is None and _is_mapping(branch) and field is None:
        query, branch = branch, None
    if query is not None:
        _require(_is_mapping(query), "lookup_query", "lookup query must be a mapping")
        branch = query.get("branch")
        field = query.get("field")
        j = query.get("j", query.get("ordinary_Z_order", query.get("Zorder", query.get("Z_order"))))
        radial_order = query.get("radial_order", radial_order)
        phase_mode = query.get("phase_mode", query.get("phase", phase_mode))
    if phase_mode != "ordinary_Z":
        _api_refusal(branch=branch, field=field, order=j, radial_order=radial_order, phase_mode=phase_mode, reason="only ordinary_Z source bindings are saved", control="held_phase")
    if not _is_int(radial_order) or radial_order != 1:
        _api_refusal(branch=branch, field=field, order=j, radial_order=radial_order, phase_mode=phase_mode, reason="only first radial order 1 is saved", control="second_radial")
    if branch not in BRANCHES:
        _api_refusal(branch=branch, field=field, order=j, radial_order=radial_order, phase_mode=phase_mode, reason="branch must be exactly background or modified", control="wrong_branch")
    if field not in FIELDS:
        _api_refusal(branch=branch, field=field, order=j, radial_order=radial_order, phase_mode=phase_mode, reason="field is outside the native O3 field set", control="wrong_field")
    if not _is_int(j) or j not in ORDINARY_Z_ORDERS:
        _api_refusal(branch=branch, field=field, order=j, radial_order=radial_order, phase_mode=phase_mode, reason="ordinary-Z order is not in the compiled 0 through 4 prefix", control="order5")
    root: Any = bindings
    if _is_mapping(bindings.get("bindings")):
        root = bindings.get("bindings")
    elif _is_mapping(bindings.get("physical_bindings")):
        root = bindings.get("physical_bindings")
    branch_rows = root.get(branch) if _is_mapping(root) else None
    field_rows = branch_rows.get(field) if _is_mapping(branch_rows) else None
    value = field_rows.get(str(j), field_rows.get(j)) if _is_mapping(field_rows) else None
    if not _is_mapping(value):
        _api_refusal(branch=branch, field=field, order=j, radial_order=radial_order, phase_mode=phase_mode, reason="requested exact source binding is absent", control="missing_binding")
    return value


def construct_alias(
    bindings: Mapping[str, Any],
    branch: Any = None,
    field: Any = None,
    j: Any = None,
    radial_order: Any = 1,
    phase_mode: Any = "ordinary_Z",
    *,
    query: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    """Return an exact saved binding; no new expression is constructed."""
    return lookup(bindings, branch, field, j, radial_order, phase_mode, query=query)


constructalias = construct_alias
construct = construct_alias


def _result(
    journal: _Journal,
    *,
    passed: bool,
    status: str,
    failure: Mapping[str, Any] | None,
    source_validation: Mapping[str, Any] | None,
    generic_data: Mapping[str, Any] | None,
    parent_data: Mapping[str, Any] | None,
    alias_rows: list[Mapping[str, Any]],
    density_rows: list[Mapping[str, Any]],
    rewrite_rows: list[Mapping[str, Any]],
    saved_rows: list[Mapping[str, Any]],
    final_bindings: Mapping[str, Any] | None,
) -> dict[str, Any]:
    flags = _flags(passed)
    counts = _counts(journal)
    payload: dict[str, Any] = {
        "source_validation": _clone(source_validation) if source_validation is not None else None,
        "generic_density_scope": _clone(generic_data) if generic_data is not None else None,
        "parent_source_data": _clone({key: value for key, value in parent_data.items() if key not in ("cached_data", "source_proof_rows")}) if parent_data is not None else None,
        "cached_data": _clone(parent_data.get("cached_data")) if parent_data is not None else None,
        "source_proof_rows": _clone(parent_data.get("source_proof_rows")) if parent_data is not None else None,
        "source_proof_rows_scope": "retained parent namespace; fresh flux rewrites are recorded separately",
        "analytic_alias_rows": _clone(alias_rows),
        "density_flux_rows": _clone(density_rows),
        "derivative_rewrite_rows": _clone(rewrite_rows),
        "saved_m_y_alias_rows": _clone(saved_rows),
        "physical_binding_rows": _event_rows(journal, "physical_binding_row"),
        "actual_negative_api_refusals": _event_rows(journal, "api_refusal"),
        "bindings": _clone(final_bindings) if final_bindings is not None else None,
        "physical_bindings": _clone(final_bindings) if final_bindings is not None else None,
        "source_scope": {
            "chart": EXPECTED_CHART,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "radial_order": 1,
            "phase_mode": "ordinary_Z",
            "history_keys": list(HISTORY_KEYS),
            "fields": list(FIELDS),
            "source_jet_semantics": "genuine ordinary C4 source jet",
            "conditional": True,
            "pointwise_polynomial_schema_consumed_as_data": True,
            "ordinary_Z3_Z4_fresh_extension": True,
            "unconditional_ordinary_Z4_source_admission": False,
            "incoming_x_independence_proved": False,
            "phase_seam_regularity_proved": False,
            "numeric_admission": False,
            "source_unit": "native_R",
            "physical_time_map_instantiated": False,
        },
    }
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "status": status,
        "passed": passed,
        "failure": _clone(failure) if failure is not None else None,
        "runtime_model_metadata": None,
        "metadata_only": False,
        "source_only": True,
        "conditional": True,
        "conditional_on_declared_source_atoms": True,
        "physical_layer_completed": False,
        "full_goal_complete": False,
        "numeric_evaluation": False,
        "proof_count": 1 if passed else 0,
        "actual_counts": counts,
        "callback_event_count": len(journal.events),
        "callback_events": _clone(journal.events),
        "new_flags": flags,
        **flags,
        "payload": payload,
    }


def derive(source: Mapping[str, Any], emit: Callable[[Mapping[str, Any]], Any] | Any = None) -> dict[str, Any]:
    """Build the fresh O3 density-flux source family exactly once."""
    journal = _Journal(emit)
    source_validation: Mapping[str, Any] | None = None
    generic_data: Mapping[str, Any] | None = None
    parent_data: Mapping[str, Any] | None = None
    alias_rows: list[Mapping[str, Any]] = []
    density_rows: list[Mapping[str, Any]] = []
    rewrite_rows: list[Mapping[str, Any]] = []
    saved_rows: list[Mapping[str, Any]] = []
    final_bindings: dict[str, Any] | None = None
    try:
        # This is intentionally the first callback, before any source identity validation.
        journal.emit(
            "family_entered",
            {
                "family": "native_O3_density_flux_identity_v1",
                "source_only": True,
                "conditional": True,
                "chart": EXPECTED_CHART,
                "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
                "radial_order": 1,
                "phase_mode": "ordinary_Z",
                "source_jet_semantics": "genuine ordinary C4 source jet",
                "fresh_O3_source_predicate": True,
                "pointwise_polynomial_schema_consumed_as_data": True,
                "ordinary_Z3_Z4_fresh_extension": True,
                "unconditional_ordinary_Z4_source_admission": False,
                "numeric_evaluation": False,
                "N_materialized": False,
            },
        )
        literals, source_validation = _wrapper_literals(source)
        raw = literals["density_flux_literals"]
        generic_data = _validate_generic_scope(literals["cached_density_scope"])
        raw_validation, raw_parts = _validate_ast_and_graph(raw)
        source_validation = {**dict(source_validation), **raw_validation, "wrapper_pins": dict(source_validation["source_wrapper_pins"])}
        parent_parts, parent_validation = _parent_payload(literals["saved_history_result"])
        parent_proof_rows = parent_parts["payload"].get("source_proof_rows")
        _require(_is_mapping(parent_proof_rows), "parent_proof_namespace", "retained parent proof-row namespace is missing")
        parent_data = {**parent_validation, "cached_data": parent_parts["cached"], "source_proof_rows": parent_proof_rows}
        source_validation = {**dict(source_validation), "parent": parent_validation}
        handles = raw_parts["literal_rows"]
        source_handles = raw_validation["handles"]
        heads = raw_parts["heads"]

        journal.emit(
            "scope_binding",
            {
                "family": "native_O3_density_flux_identity_v1",
                "chart": EXPECTED_CHART,
                "source_wrappers": list(EXPECTED_SOURCE_SHA256),
                "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
                "history_keys": list(HISTORY_KEYS),
                "fields": list(FIELDS),
                "radial_order": 1,
                "phase_mode": "ordinary_Z",
                "source_jet_semantics": "genuine ordinary C4 source jet",
                "selected_literal_rows": 120,
                "source_graph_heads": 84,
                "source_graph_constructor_records": 189,
                "opaque_N_node": N_NODE,
                "inverse_N_node": INVERSE_N_NODE,
                "native_radius_node": RADIUS_NODE,
                "source_unit": "native_R",
                "physical_time_map_instantiated": False,
                "conditional": True,
                "global_jet_oracle": False,
                "numeric_admission": False,
            },
        )

        # Fresh functional E_N=E+r aliases.  The exp/expm1 value row is kept
        # distinct from the shared higher derivative source rows.
        for order in ORDINARY_Z_ORDERS:
            row = _alias_row(source_handles, order)
            alias_rows.append(row)
            journal.emit("analytic_alias_row", row)
        _require(len(alias_rows) == 5, "alias_scope", "five analytic alias rows were not built")

        # Fresh pointwise density/flux source rows, lifted with the fixed
        # ordinary-Z binomial table; no old jet builder is called.
        for history in HISTORY_KEYS:
            for order in ORDINARY_Z_ORDERS:
                row = _density_row(source_handles, generic_data, parent_parts["duhamel"], history, order)
                density_rows.append(row)
                journal.emit("density_flux_row", row)
        _require(len(density_rows) == 25, "density_scope", "twenty-five density-flux rows were not built")
        density_by_key = {(str(row["history_key"]), int(row["ordinary_Z_order"])): row for row in density_rows}

        # Rewrite only the current conditional modified Duhamel rows; retain
        # all incoming pointers and the parent's regular interior premises.
        for history in HISTORY_KEYS:
            for order in ORDINARY_Z_ORDERS:
                prior = parent_parts["modified"][(history, order)]
                row = _rewrite_row(prior, density_by_key[(history, order)], prior, history, order)
                rewrite_rows.append(row)
                journal.emit("derivative_rewrite_row", row)
        _require(len(rewrite_rows) == 25, "rewrite_scope", "twenty-five derivative rewrite rows were not built")
        rewrite_by_key = {(str(row["history_key"]), int(row["ordinary_Z_order"])): row for row in rewrite_rows}

        # Verify the saved total radial m_y source before recording the fresh
        # VN - modified-m compatibility expression.
        for order in ORDINARY_Z_ORDERS:
            row = _saved_m_y_row(source_handles, heads, parent_parts["bindings"], parent_parts["modified"], order)
            saved_rows.append(row)
            journal.emit("saved_m_y_alias_row", row)
        _require(len(saved_rows) == 5, "saved_m_y_scope", "five saved m_y alias rows were not built")

        final_bindings = _clone(parent_parts["bindings"])
        rewritten_keys: set[tuple[str, str, int]] = set()
        for branch in BRANCHES:
            for field in FIELDS:
                for order in ORDINARY_Z_ORDERS:
                    parent_row = _parent_binding(parent_parts["bindings"], branch, field, order)
                    rewritten = branch == "modified" and field in HISTORY_KEYS
                    rewritten_pressure = branch == "modified" and field == "pressure"
                    if rewritten:
                        rewrite = rewrite_by_key[(field, order)]
                        p_row = _rewrite_row(parent_row, density_by_key[(field, order)], parent_parts["modified"][(field, order)], field, order)
                        qrow = _quotient_row(parent_row, field, order, rewrite["rewritten_source_expression"]["numerator"])
                        qrow["rewritten_derivative_source"] = _clone(rewrite)
                        row = qrow
                        rewritten_keys.add((branch, field, order))
                        action = "rewritten_flux"
                    elif rewritten_pressure:
                        p_density = density_by_key[("p", order)]
                        p_flux = p_density["flux_expression"]
                        p0_binding = _parent_binding(parent_parts["bindings"], "modified", "P0", order)
                        p0_node = p0_binding.get("raw_node_bindings", {}).get("field") if _is_mapping(p0_binding.get("raw_node_bindings")) else None
                        p0_zero = _rational(0, source="cached_declared_P0_radial_zero")
                        pressure_numerator = _sum(p0_zero, p_flux)
                        _validate_expression_kinds(pressure_numerator)
                        row = _quotient_row(parent_row, "pressure", order, pressure_numerator, pressure=True, p0_node=p0_node if _is_int(p0_node) else None)
                        row["rewritten_pressure_source"] = {
                            "p_derivative_expression": _clone(p_flux),
                            "cached_declared_P0_radial_zero": p0_zero,
                            "P0_ordinary_Z_value_claimed": False,
                            "ordinary_Z0_claimed": False,
                            "source_jet_semantics": "genuine ordinary C4 source jet",
                        }
                        rewritten_keys.add((branch, field, order))
                        action = "rewritten_flux"
                    else:
                        row = _clone(parent_row)
                        action = "retained_parent"
                    field_rows = final_bindings[branch][field]
                    key = str(order) if str(order) in field_rows else order
                    field_rows[key] = _clone(row)
                    journal.emit("physical_binding_row", {"action": action, "binding": _clone(row)})
        _require(len(rewritten_keys) == 30, "physical_scope", "exactly thirty physical rows were not rewritten")
        _require(len(journal.events) >= 0, "journal_shape", "journal shape changed")

        controls = (
            ("held_phase", {"branch": "background", "field": "m", "j": 0, "radial_order": 1, "phase_mode": "held_phase"}),
            ("second_radial", {"branch": "background", "field": "m", "j": 0, "radial_order": 2, "phase_mode": "ordinary_Z"}),
            ("wrong_branch", {"branch": "wrong_branch", "field": "m", "j": 0, "radial_order": 1, "phase_mode": "ordinary_Z"}),
            ("order5", {"branch": "background", "field": "m", "j": 5, "radial_order": 1, "phase_mode": "ordinary_Z"}),
        )
        for control, query in controls:
            try:
                lookup(final_bindings, query=query)
            except _APIRefusal as refusal:
                row = dict(refusal.row)
                row["control"] = control
                journal.emit("api_refusal", row)
            else:
                raise _KernelFailure("negative_control_admitted", f"negative API control {control} was admitted")

        journal.emit(
            "family_completed",
            {
                "family": "native_O3_density_flux_identity_v1",
                "analytic_alias_rows": 5,
                "density_flux_rows": 25,
                "derivative_rewrite_rows": 25,
                "saved_m_y_alias_rows": 5,
                "physical_binding_rows": 90,
                "rewritten_physical_bindings": 30,
                "retained_physical_bindings": 60,
                "negative_api_refusals": 4,
                "new_mathematical_proofs": 1,
                "conditional": True,
                "source_only": True,
                "native_field_derivative_identity_proved": False,
                "global_across_phase_seams_proved": False,
                "actual_incoming_x_independence_proved": False,
                "unconditional_ordinary_Z4_source_admission": False,
            },
        )
        return _result(
            journal,
            passed=True,
            status="complete",
            failure=None,
            source_validation=source_validation,
            generic_data=generic_data,
            parent_data=parent_data,
            alias_rows=alias_rows,
            density_rows=density_rows,
            rewrite_rows=rewrite_rows,
            saved_rows=saved_rows,
            final_bindings=final_bindings,
        )
    except _CallbackFailure as exc:
        failure = {"code": "diagnostic_callback_failure", "message": str(exc.cause), "event": exc.event, "prefix_preserved": True}
    except _APIRefusal as exc:
        failure = {"code": "unexpected_api_refusal", "message": str(exc), "row": _clone(exc.row), "prefix_preserved": True}
    except _KernelFailure as exc:
        failure = {"code": exc.code, "message": exc.message, "details": _clone(exc.details), "prefix_preserved": True}
    except BaseException as exc:
        failure = {"code": "unexpected_exception", "message": str(exc), "exception_type": type(exc).__name__, "prefix_preserved": True}
    return _result(
        journal,
        passed=False,
        status="failed",
        failure=failure,
        source_validation=source_validation,
        generic_data=generic_data,
        parent_data=parent_data,
        alias_rows=alias_rows,
        density_rows=density_rows,
        rewrite_rows=rewrite_rows,
        saved_rows=saved_rows,
        final_bindings=final_bindings,
    )


__all__ = [
    "COMPONENT",
    "INPUT_SCHEMA",
    "RESULT_SCHEMA",
    "EXPECTED_SOURCE_SHA256",
    "derive",
    "lookup",
    "construct_alias",
    "constructalias",
    "construct",
]
