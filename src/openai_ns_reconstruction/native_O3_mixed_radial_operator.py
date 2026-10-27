"""Inert source-operator kernel for the native ``O3_transition`` chart.

The kernel consumes already-parsed, hash-pinned literal packets.  It builds
typed *source expressions* for the rows that are actually present in those
packets; it does not evaluate a field, differentiate a phase, decode an
archive, materialize ``N``, or call an older provider/checker.  A compiled row
therefore records a quotient of two native graph references and the source
premises under which that quotient may later be consumed.

The saved ``R_Z = 0`` fact is carried as cached data only.  In particular, the
two flags for an actual field derivative identity and complete mixed field
coverage remain false.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


COMPONENT = "native_O3_mixed_radial_operator"
INPUT_SCHEMA = "native_O3_mixed_radial_operator_input_v1"
RESULT_SCHEMA = "original_native_O3_mixed_radial_operator"

EXPECTED_SOURCE_SHA256 = {
    "saved_Z_result": "e0ef019b51898a83d01b58940face6023629ec3b4ef245bca903d4d20a218eb7",
    "background_literals": "2748333d5dc25c5a34b23235555dc7e285f52c05e5b1ae2995e50aa2281bca14",
    "modified_literals": "d1834f1ed518b8a4525edd110fc3fc9d551db6504090d7668db3dad027cc3906",
}
EXPECTED_PARENT_COMMIT = "72dccb761638c0e14001bf2e4fc76c3f6dff7762"
EXPECTED_PARENT_REPORT_SHA256 = "0587c2c5ec692f6090f02511b56adbb2a79238836081c802d1211eb3348b9d3c"
EXPECTED_ARCHIVE_SHA256 = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
EXPECTED_CHART = "O3_transition"
EXPECTED_COORDINATE_NODE = 241536
EXPECTED_RADIUS_NODE = 241787

FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
MODIFIED_FIELDS = ("E", "V", "m")
ORDINARY_Z_ORDERS = (0, 1, 2, 3, 4)
MISSING_MODIFIED_FIELDS = ("h", "k", "e", "p", "P0", "pressure")

EXPECTED_BACKGROUND_NODES = {
    "E": (241563, 241564, 241565, 241566, 241567),
    "V": (0, 0, 0, 0, 0),
    "m": (244746, 244747, 244748, 244749, 244750),
    "h": (244767, 244768, 244769, 244770, 244771),
    "k": (244781, 244782, 244783, 244784, 244785),
    "e": (244795, 244796, 244797, 244798, 244799),
    "p": (241604, 241605, 241606, 241607, 241608),
    "P0": (0, 0, 0, 0, 0),
    "pressure": (241604, 241605, 241606, 241607, 241608),
}
EXPECTED_MODIFIED_NODES = {
    "E": (249053, 249056, 249060, 249065, 249071),
    "V": (249085, 249086, 249087, 249088, 249089),
    "m": (249097, 249098, 249099, 249100, 249101),
}

COUNT_NAMES = (
    "new_algebra_family_calls",
    "new_compiled_operator_bindings",
    "unavailable_source_rows",
    "actual_negative_api_refusals",
    "old_math_calls",
    "old_provider_calls",
    "numeric_evaluation_calls",
    "N_materialization_calls",
    "source_value_identity_calls",
    "new_mathematical_proofs",
)

FLAG_NAMES = (
    "native_O3_first_radial_mixed_source_operator_compiled",
    "native_O3_total_phase_source_rows_retained",
    "native_field_derivative_identity_proved",
    "full_mixed_field_rows_complete",
    "native_Rm_Z_zero_proved",
    "actual_native_field_pullback_instantiated",
    "fractional_part_differentiated",
    "N_materialized",
    "physical_layer_completed",
    "full_goal_complete",
    "numeric_closure",
    "temporal_recursion",
)


class _KernelFailure(Exception):
    """Controlled validation or source-shape failure."""

    def __init__(self, code: str, message: str, details: Any = None) -> None:
        super().__init__(message)
        self.code = code
        self.message = message
        self.details = details


class _CallbackFailure(Exception):
    """A diagnostic callback failed after an event was journalled."""

    def __init__(self, event: str, cause: BaseException) -> None:
        super().__init__(str(cause))
        self.event = event
        self.cause = cause


class _APIRefusal(Exception):
    """A requested lookup/construct operation is outside the saved interface."""

    def __init__(self, row: Mapping[str, Any]) -> None:
        super().__init__(str(row.get("reason", "API request refused")))
        self.row = dict(row)


def _is_mapping(value: Any) -> bool:
    return isinstance(value, Mapping)


def _is_int(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _clone(value: Any) -> Any:
    """Clone only JSON-shaped data; opaque source values are never interpreted."""

    if isinstance(value, Mapping):
        return {key: _clone(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_clone(item) for item in value]
    if isinstance(value, tuple):
        return [_clone(item) for item in value]
    return value


def _new_counts() -> dict[str, int]:
    return {name: 0 for name in COUNT_NAMES}


def _flags(success: bool) -> dict[str, bool]:
    flags = {name: False for name in FLAG_NAMES}
    if success:
        flags["native_O3_first_radial_mixed_source_operator_compiled"] = True
        flags["native_O3_total_phase_source_rows_retained"] = True
    return flags


def _pin_ok(name: str, observed: Any) -> bool:
    if not isinstance(observed, str) or len(observed) != 64:
        return False
    if any(char not in "0123456789abcdefABCDEF" for char in observed):
        return False
    expected = EXPECTED_SOURCE_SHA256[name]
    return expected is None or observed.lower() == expected


def _require(condition: bool, code: str, message: str, details: Any = None) -> None:
    if not condition:
        raise _KernelFailure(code, message, details)


def _wrapper_literals(source: Mapping[str, Any]) -> tuple[dict[str, Mapping[str, Any]], dict[str, Any]]:
    _require(source.get("schema") == INPUT_SCHEMA, "input_schema_mismatch", "source.schema changed")
    _require(source.get("component") == COMPONENT, "input_component_mismatch", "source.component changed")
    wrappers = source.get("source_wrappers")
    names = tuple(EXPECTED_SOURCE_SHA256)
    _require(_is_mapping(wrappers), "source_wrappers_missing", "source_wrappers must be a mapping")
    _require(set(wrappers) == set(names), "source_wrapper_names_mismatch", "source_wrappers must contain exactly the three frozen wrappers")
    literals: dict[str, Mapping[str, Any]] = {}
    fingerprints: dict[str, str] = {}
    for name in names:
        wrapper = wrappers.get(name)
        _require(_is_mapping(wrapper), "source_wrapper_shape", f"source_wrappers.{name} must be a mapping")
        _require(set(wrapper) == {"sha256", "literal"}, "source_wrapper_shape", f"source_wrappers.{name} must expose sha256 and literal")
        observed = wrapper.get("sha256")
        _require(_pin_ok(name, observed), "source_wrapper_pin_mismatch", f"source_wrappers.{name}.sha256 is not the frozen pin")
        literal = wrapper.get("literal")
        _require(_is_mapping(literal), "source_literal_shape", f"source_wrappers.{name}.literal must be parsed JSON data")
        literals[name] = literal
        fingerprints[name] = str(observed).lower()

    parent_commit = source.get("parent_commit")
    if parent_commit is not None:
        _require(parent_commit == EXPECTED_PARENT_COMMIT, "parent_commit_mismatch", "parent_commit is not the frozen source head")
    archive_sha = source.get("source_archive_sha256")
    if archive_sha is not None:
        _require(archive_sha == EXPECTED_ARCHIVE_SHA256, "archive_pin_mismatch", "source archive pin changed")
    parent_report = source.get("parent_report")
    if parent_report is not None:
        _require(_is_mapping(parent_report), "parent_report_shape", "parent_report must be hash-only parsed data")
        _require(parent_report.get("sha256") == EXPECTED_PARENT_REPORT_SHA256, "parent_report_pin_mismatch", "parent report pin changed")
    policy = source.get("policy")
    if policy is not None:
        _require(_is_mapping(policy), "policy_shape", "policy must be parsed data")
        _require(policy.get("source_only") is True, "source_only_policy_mismatch", "source-only policy changed")
        _require(policy.get("full_goal_complete") is False, "goal_policy_mismatch", "full-goal policy changed")
        _require(policy.get("numeric_evaluation") is False, "numeric_policy_mismatch", "numeric evaluation policy changed")
    parent_certificate = source.get("parent_certificate")
    if parent_certificate is not None:
        _require(_is_mapping(parent_certificate), "parent_certificate_shape", "parent_certificate must be parsed data")
    return literals, {
        "input_schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "source_wrapper_names": list(names),
        "source_wrapper_pins": fingerprints,
        "parent_commit": parent_commit,
        "source_archive_sha256": archive_sha,
        "parent_report_sha256": parent_report.get("sha256") if _is_mapping(parent_report) else None,
        "parent_certificate_present": parent_certificate is not None,
        "source_modules_imported": False,
        "project_modules_imported": False,
        "archive_opened_or_decoded": False,
        "source_value_evaluation": False,
        "fractional_part_differentiated": False,
        "N_materialized": False,
    }


def _row_by_order(rowset: Any, order: int) -> Mapping[str, Any] | None:
    if not _is_mapping(rowset):
        return None
    rows = rowset.get("rows")
    if isinstance(rows, (list, tuple)):
        for row in rows:
            if _is_mapping(row) and row.get("order") == order:
                return row
        return None
    if _is_mapping(rows):
        row = rows.get(order)
        if row is None:
            row = rows.get(str(order))
        return row if _is_mapping(row) else None
    row = rowset.get(order)
    if row is None:
        row = rowset.get(str(order))
    return row if _is_mapping(row) else None


def _validate_rowset(
    rowset: Any,
    *,
    field: str,
    branch: str,
    order_bound: int,
    source_status: str,
    total_row: bool,
    expected_nodes: tuple[int, ...] | None = None,
) -> dict[int, dict[str, Any]]:
    _require(_is_mapping(rowset), "rowset_missing", f"{branch}.{field} rowset is missing")
    available_count = rowset.get("available_row_count")
    _require(_is_int(available_count) and available_count >= order_bound, "rowset_order_bound", f"{branch}.{field} does not expose the genuine order bound")
    available_orders = rowset.get("available_orders")
    if available_orders is not None:
        _require(list(available_orders[:5]) == list(ORDINARY_Z_ORDERS), "rowset_orders_mismatch", f"{branch}.{field} selected orders changed")
    _require(rowset.get("selected_orders", list(ORDINARY_Z_ORDERS)) == list(ORDINARY_Z_ORDERS), "rowset_selected_orders_mismatch", f"{branch}.{field} must select ordinary Z orders 0 through 4")
    _require(rowset.get("selected_row_count", 5) == 5, "rowset_selected_count_mismatch", f"{branch}.{field} selected row count changed")
    _require(rowset.get("rows_are_genuine") is not False, "rowset_not_genuine", f"{branch}.{field} rows must be genuine source rows")
    _require(rowset.get("padding_rows_added") is not True, "rowset_padding_rejected", f"{branch}.{field} cannot use padding rows")
    selected: dict[int, dict[str, Any]] = {}
    for order in ORDINARY_Z_ORDERS:
        item = _row_by_order(rowset, order)
        _require(item is not None, "source_row_missing", f"{branch}.{field} order {order} row is missing")
        node = item.get("node")
        _require(_is_int(node) and node >= 0, "source_node_malformed", f"{branch}.{field} order {order} node is not a graph reference")
        if expected_nodes is not None:
            _require(node == expected_nodes[order], "source_node_pin_mismatch", f"{branch}.{field} order {order} node changed")
        source_order = item.get("source_order", order)
        _require(source_order == order, "source_order_mismatch", f"{branch}.{field} source order changed")
        selected[order] = {
            "order": order,
            "source_order": source_order,
            "node": node,
            "archive_json_path": item.get("archive_json_path", rowset.get("archive_json_path")),
            "literal_json_path": item.get("literal_json_path", rowset.get("literal_json_path")),
            "source_status": source_status,
            "branch": branch,
            "field": field,
            "order_bound": order_bound,
            "rows_are_genuine": True,
            "padding_added": False,
            "total_row": total_row,
        }
    return selected


def _background_rowsets(packet: Mapping[str, Any]) -> tuple[dict[str, dict[int, dict[str, Any]]], Mapping[str, Any]]:
    section = packet.get("declared_background_radial_rows")
    _require(_is_mapping(section), "background_rows_section_missing", "declared_background_radial_rows is required")
    rows = section.get("rows")
    _require(_is_mapping(rows), "background_rows_missing", "declared_background_radial_rows.rows is required")
    fields = section.get("field_keys")
    _require(list(fields) == list(FIELDS), "background_field_order_mismatch", "background field order changed")
    source_status = str(section.get("source_status", "declared source-equation rows only; derivative identity unproved"))
    result: dict[str, dict[int, dict[str, Any]]] = {}
    for field in FIELDS:
        result[field] = _validate_rowset(
            rows.get(field),
            field=field,
            branch="background",
            order_bound=7,
            source_status=source_status,
            total_row=False,
            expected_nodes=EXPECTED_BACKGROUND_NODES[field],
        )
    _require(section.get("selected_orders_are_genuine") is not False, "background_orders_not_genuine", "background selected orders are not genuine")
    _require(section.get("padding_added") is not True, "background_padding_rejected", "background packet reports padding")
    return result, section


def _modified_total_section(packet: Mapping[str, Any]) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    candidates = (
        "modified_total_radial_rows",
        "saved_total_radial_rows",
        "actual_saved_total_radial_rows",
        "total_radial_rows",
        "total_log_radial_physical_y_output",
        "saved_total_radial_y_output",
    )
    for name in candidates:
        section = packet.get(name)
        if not _is_mapping(section):
            continue
        rows = section.get("rows")
        if not _is_mapping(rows):
            rows = section.get("derivative_rows")
        if _is_mapping(rows) and any(key in rows for key in ("E_y", "V_y", "m_y")):
            return section, rows
    # A few producers place the section under a branch object.  Only exact
    # total-y keys are admitted; slow_A_y/slow_B_y are intentionally ignored.
    for container_name in ("modified", "modified_branch", "saved_modified"):
        container = packet.get(container_name)
        if not _is_mapping(container):
            continue
        for name in candidates:
            section = container.get(name)
            if not _is_mapping(section):
                continue
            rows = section.get("rows") or section.get("derivative_rows")
            if _is_mapping(rows) and any(key in rows for key in ("E_y", "V_y", "m_y")):
                return section, rows
    # The frozen modified packet intentionally stores only compact root arrays
    # (rather than a fabricated background-like rowset).  Convert those root
    # arrays into the same local rowset shape without touching any graph
    # children or arithmetic values.
    compact = packet.get("available_orders_and_root_counts")
    if _is_mapping(compact):
        modified_chart = compact.get("modified_chart")
        radial_total = modified_chart.get("radial_total") if _is_mapping(modified_chart) else None
        if _is_mapping(radial_total) and any(key in radial_total for key in ("E_y", "V_y", "m_y")):
            rows: dict[str, dict[str, Any]] = {}
            for key in ("E_y", "V_y", "m_y"):
                entry = radial_total.get(key)
                if not _is_mapping(entry):
                    continue
                roots = entry.get("roots")
                orders = entry.get("orders")
                if not isinstance(roots, (list, tuple)):
                    continue
                archive_json_path = f"$.true_outer_C6_C5_bindings.charts.O3_transition.radial.{key}"
                literal_json_path = f"$.available_orders_and_root_counts.modified_chart.radial_total.{key}"
                rows[key] = {
                    "archive_json_path": archive_json_path,
                    "literal_json_path": literal_json_path,
                    "available_row_count": entry.get("count", len(roots)),
                    "available_orders": list(orders) if isinstance(orders, (list, tuple)) else list(range(len(roots))),
                    "selected_orders": list(ORDINARY_Z_ORDERS),
                    "selected_row_count": len(ORDINARY_Z_ORDERS),
                    "rows_are_genuine": True,
                    "padding_rows_added": False,
                    "rows": [
                        {
                            "order": order,
                            "source_order": order,
                            "node": roots[order],
                            "archive_json_path": archive_json_path,
                            "literal_json_path": literal_json_path,
                        }
                        for order in ORDINARY_Z_ORDERS
                        if order < len(roots)
                    ],
                }
            if all(key in rows for key in ("E_y", "V_y", "m_y")):
                return {
                    "source_status": "actual saved total helper output; total E_y/V_y/m_y rows only",
                    "archive_json_path": "$.true_outer_C6_C5_bindings.charts.O3_transition.radial",
                    "literal_json_path": "$.available_orders_and_root_counts.modified_chart.radial_total",
                    "rows": rows,
                    "metadata_values": {
                        "total_fast_phase_terms_included": packet["held_phase_and_fast_response_literal_roots"]["radial_flags"]["total_fast_phase_terms_included"],
                        "physical_y_equals_native_partial_divided_by_Jacobian": packet["held_phase_and_fast_response_literal_roots"]["radial_flags"]["physical_y_equals_native_partial_divided_by_Jacobian"],
                    },
                }, rows
    raise _KernelFailure("modified_total_rows_section_missing", "modified packet has no saved total E_y/V_y/m_y row section")


def _modified_rowsets(packet: Mapping[str, Any]) -> tuple[dict[str, dict[int, dict[str, Any]]], Mapping[str, Any], dict[str, Any]]:
    binding_scope = packet.get("binding_scope")
    if _is_mapping(binding_scope):
        if binding_scope.get("native_chart") is not None:
            _require(binding_scope.get("native_chart") == EXPECTED_CHART, "modified_chart_mismatch", "modified packet chart changed")
        if binding_scope.get("chart_coordinate_node") is not None:
            _require(binding_scope.get("chart_coordinate_node") == EXPECTED_COORDINATE_NODE, "modified_coordinate_mismatch", "modified packet coordinate node changed")
    section, rows = _modified_total_section(packet)
    source_status = str(section.get("source_status", "actual saved total helper output"))
    held_phase = packet.get("held_phase_and_fast_response_literal_roots")
    radial_flags = held_phase.get("radial_flags") if _is_mapping(held_phase) else None
    _require(_is_mapping(radial_flags), "modified_fastphase_metadata_missing", "modified packet radial_flags are required")
    total_fast_phase_terms_included = radial_flags.get("total_fast_phase_terms_included")
    physical_y_equals_native_partial_divided_by_Jacobian = radial_flags.get("physical_y_equals_native_partial_divided_by_Jacobian")
    _require(total_fast_phase_terms_included is True, "modified_fastphase_terms_missing", "saved total rows must record included fast-phase terms")
    _require(physical_y_equals_native_partial_divided_by_Jacobian is True, "modified_jacobian_metadata_missing", "saved total rows must record the native/Jacobian relation")
    selected: dict[str, dict[int, dict[str, Any]]] = {}
    for field in MODIFIED_FIELDS:
        key = f"{field}_y"
        selected[field] = _validate_rowset(
            rows.get(key),
            field=field,
            branch="modified",
            order_bound=6,
            source_status=source_status,
            total_row=True,
            expected_nodes=EXPECTED_MODIFIED_NODES[field],
        )
    metadata_values = section.get("metadata_values")
    if not _is_mapping(metadata_values):
        metadata_values = section
    fastphase = {
        "total_fast_phase_terms_included": total_fast_phase_terms_included,
        "physical_y_equals_native_partial_divided_by_Jacobian": physical_y_equals_native_partial_divided_by_Jacobian,
        "slow_A_y_B_y_are_not_total_E_or_V": True,
        "slow_A_y/B_y_alone_not_accepted_as_total_E_or_V": True,
        "total_E_y_and_V_y_already_include_A_phi_B_phi": True,
        "total_already_includes_A_phi/B_phi_and_N_once": True,
        "N_included_once_in_saved_total_helper": True,
        "extra_N_factor_generated": False,
        "fractional_part_differentiated": False,
        "source_status": source_status,
    }
    return selected, section, fastphase


def _radius_ref(packet: Mapping[str, Any]) -> dict[str, Any]:
    chart = packet.get("saved_chart_definition")
    radius = chart.get("radius") if _is_mapping(chart) else None
    if not _is_mapping(radius):
        radius = packet.get("radius")
    if _is_mapping(radius) and radius.get("node") is not None:
        _require(radius.get("node") == EXPECTED_RADIUS_NODE, "radius_node_mismatch", "saved radius node changed")
        path = radius.get("archive_json_path")
    else:
        path = None
    return {
        "kind": "native_graph_ref",
        "role": "native_radius",
        "node": EXPECTED_RADIUS_NODE,
        "archive_json_path": path,
        "chart": EXPECTED_CHART,
        "positiveR_premise": "regular native radius R > 0 (source premise; no value evaluated)",
        "coordinate_frame": EXPECTED_CHART,
    }


def _cached_RZ_ref(saved_result: Mapping[str, Any]) -> dict[str, Any]:
    _require(saved_result.get("passed") is True, "saved_Z_result_not_passed", "saved Z result must retain passed=true")
    new_flags = saved_result.get("new_flags")
    _require(_is_mapping(new_flags), "saved_Z_result_flags_missing", "saved Z result new_flags are required")
    _require(
        new_flags.get("native_fixed_radius_Z_frame_specialized") is True,
        "saved_Z_result_radius_frame_missing",
        "saved Z result must retain native fixed-radius Z frame specialization",
    )
    payload = saved_result.get("payload") if _is_mapping(saved_result) else None
    target_certificates = payload.get("target_certificates") if _is_mapping(payload) else None
    _require(isinstance(target_certificates, (list, tuple)), "saved_Z_result_certificates_missing", "saved Z target certificates are required")
    native_R_certificate = None
    for certificate in target_certificates:
        if _is_mapping(certificate) and certificate.get("target") == "nativeR" and certificate.get("root") == EXPECTED_RADIUS_NODE:
            native_R_certificate = dict(certificate)
            break
    _require(native_R_certificate is not None, "saved_Z_result_nativeR_certificate_missing", "saved nativeR target certificate is required")
    frame_rows = payload.get("frame_specialization_rows") if _is_mapping(payload) else None
    _require(isinstance(frame_rows, (list, tuple)), "saved_Z_result_frame_rows_missing", "saved frame specialization rows are required")
    selected_frame_rows = []
    for row in frame_rows:
        if _is_mapping(row) and row.get("branch") in ("background", "modified"):
            selected_frame_rows.append(dict(row))
    _require(
        len(selected_frame_rows) == 2 and {row.get("branch") for row in selected_frame_rows} == {"background", "modified"},
        "saved_Z_result_frame_branches_missing",
        "saved frame specialization rows must retain both native branches",
    )
    return {
        "value": "0",
        "identity": "R_Z=0",
        "status": "cached_data_not_reproved",
        "source_result_sha256": EXPECTED_SOURCE_SHA256["saved_Z_result"],
        "saved_result_passed": True,
        "saved_native_fixed_radius_Z_frame_specialized": True,
        "saved_native_Rm_Z_zero_flag": saved_result.get("native_Rm_Z_zero_proved"),
        "native_R_target_certificate": native_R_certificate,
        "selected_frame_specialization_rows": selected_frame_rows,
        "selected_frame_branches": [row.get("branch") for row in selected_frame_rows],
        "used_as_premise_only": True,
        "new_mathematical_proof": False,
    }


def _binding_row(
    source_row: Mapping[str, Any],
    radius_ref: Mapping[str, Any],
    cached_RZ: Mapping[str, Any],
    *,
    branch: str,
    field: str,
    order: int,
    fastphase: Mapping[str, Any] | None,
) -> dict[str, Any]:
    numerator = {
        "kind": "native_graph_ref",
        "role": "source_total_y_row",
        "source_origin": "actual_saved_total_helper_output" if branch == "modified" else "declared_background_radial_source_equation",
        "node": source_row["node"],
        "archive_json_path": source_row.get("archive_json_path"),
        "branch": branch,
        "field": field,
        "ordinary_Z_order": order,
        "coordinate_frame": EXPECTED_CHART,
    }
    expression = {
        "kind": "typed_quotient",
        "operator": "native_source_quotient",
        "operator_semantics": "supplied_total_y_source_row / native_graph_ref(radius241787)",
        "source_prerequisite": "supplied_total_y_source_row",
        "conditional_source_equation": True,
        "not_an_actual_derivative_identity": True,
        "numerator": numerator,
        "denominator": dict(radius_ref),
        "typed": True,
        "source_only": True,
        "not_evaluated": True,
        "text": f"native_graph_ref({source_row['node']}) / native_graph_ref({EXPECTED_RADIUS_NODE})",
    }
    row = {
        "branch": branch,
        "chart": EXPECTED_CHART,
        "coordinate_frame": EXPECTED_CHART,
        "field": field,
        "ordinary_Z": True,
        "ordinary_Z_order": order,
        "Zorder": order,
        "radial_order": 1,
        "phase_mode": "ordinary_Z_total",
        "total_vs_held_phase": "total",
        "held_phase": False,
        "source_definition_status": source_row.get("source_status"),
        "source_row": numerator,
        "radius": dict(radius_ref),
        "expression": expression,
        "compiled_expression": expression["text"],
        "supplied_total_y_source_row": True,
        "operator_semantics": "supplied_total_y_source_row divided by native_graph_ref(radius241787)",
        "conditional_source_equation": True,
        "not_an_actual_derivative_identity": True,
        "cached_R_Z": {
            "identity": cached_RZ["identity"],
            "value": cached_RZ["value"],
            "status": "cached_data_not_reproved",
            "source_result_sha256": cached_RZ["source_result_sha256"],
            "saved_target_certificate_path": "$.payload.target_certificates[target=nativeR,root=241787]",
            "full_cached_record_in_this_result": "$.payload.cached_R_Z",
            "used_as_premise_only": True,
            "new_mathematical_proof": False,
        },
        "regular_positiveR_premise": "R > 0",
        "genuine_order_bound": source_row.get("order_bound"),
        "source_rows_are_genuine": True,
        "padding_added": False,
        "factorial_reapplied": False,
        "source_value_evaluated": False,
        "derivative_identity_proved": False,
        "native_field_derivative_identity_proved": False,
        "full_mixed_field_rows_complete": False,
        "physical_time_or_velocity_claimed": False,
        "cartesian_residual_or_norm_claimed": False,
        "log_measure_as_3D_volume_claimed": False,
    }
    if branch == "background" and field == "P0" and source_row.get("node") == 0:
        row["source_zero"] = True
        row["zero_provenance"] = "explicit native source zero graph; not padding or substitution"
    else:
        row["source_zero"] = False
    if fastphase is not None:
        row["fastphase_provenance"] = dict(fastphase)
    return row


def _query_parts(
    branch: Any,
    field: Any,
    Zorder: Any,
    radial_order: Any,
    phase_mode: Any,
    query: Mapping[str, Any] | None,
) -> tuple[Any, Any, Any, Any, Any]:
    if query is not None:
        branch = query.get("branch")
        field = query.get("field")
        if "Zorder" in query:
            Zorder = query.get("Zorder")
        elif "ordinary_Z_order" in query:
            Zorder = query.get("ordinary_Z_order")
        elif "Z_order" in query:
            Zorder = query.get("Z_order")
        radial_order = query.get("radial_order", radial_order)
        phase_mode = query.get("phase_mode", query.get("phase", phase_mode))
        if query.get("held_phase") is True:
            phase_mode = "held_phase"
        if query.get("total_vs_held_phase") not in (None, "total"):
            phase_mode = query.get("total_vs_held_phase")
    return branch, field, Zorder, radial_order, phase_mode


def _refusal(*, branch: Any, field: Any, Zorder: Any, radial_order: Any, phase_mode: Any, reason: str, control: str | None = None) -> _APIRefusal:
    row: dict[str, Any] = {
        "decision": "refused",
        "branch": branch,
        "field": field,
        "ordinary_Z_order": Zorder,
        "Zorder": Zorder,
        "radial_order": radial_order,
        "phase_mode": phase_mode,
        "total_vs_held_phase": "held_phase" if phase_mode != "ordinary_Z" and phase_mode != "ordinary_Z_total" else "total",
        "reason": reason,
        "typed_source_refusal": True,
        "source_value_evaluated": False,
        "numeric_call": False,
    }
    if control is not None:
        row["control"] = control
    return _APIRefusal(row)


def lookup(
    bindings: Mapping[str, Any],
    branch: Any = None,
    field: Any = None,
    Zorder: Any = None,
    radial_order: Any = 1,
    *,
    phase_mode: Any = "ordinary_Z",
    query: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    """Return one existing binding, rejecting every non-exact request."""

    if query is None and _is_mapping(branch) and field is None:
        query = branch
        branch = None
    branch, field, Zorder, radial_order, phase_mode = _query_parts(
        branch, field, Zorder, radial_order, phase_mode, query
    )
    if phase_mode not in ("ordinary_Z", "ordinary_Z_total"):
        raise _refusal(branch=branch, field=field, Zorder=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="held_phase requests are not total ordinary-Z bindings", control="held_phase")
    if not _is_int(radial_order) or radial_order != 1:
        raise _refusal(branch=branch, field=field, Zorder=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="only first radial order 1 is saved", control="second_radial")
    if branch not in ("background", "modified"):
        raise _refusal(branch=branch, field=field, Zorder=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="branch must be exactly background or modified", control="wrong_branch")
    if field not in FIELDS:
        raise _refusal(branch=branch, field=field, Zorder=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="field is outside the native O3 field set")
    if not _is_int(Zorder) or Zorder not in ORDINARY_Z_ORDERS:
        raise _refusal(branch=branch, field=field, Zorder=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason="ordinary-Z order is not in the compiled 0 through 4 prefix", control="order5" if Zorder == 5 else None)
    branch_bindings = bindings.get(branch) if _is_mapping(bindings) else None
    field_bindings = branch_bindings.get(field) if _is_mapping(branch_bindings) else None
    binding = None
    if _is_mapping(field_bindings):
        binding = field_bindings.get(Zorder)
        if binding is None:
            binding = field_bindings.get(str(Zorder))
    if not _is_mapping(binding):
        control = "missing_modified_field" if branch == "modified" and field in MISSING_MODIFIED_FIELDS else None
        reason = "modified total-y source row is unavailable; no background substitution or padding is permitted" if control else "requested exact source binding is absent"
        raise _refusal(branch=branch, field=field, Zorder=Zorder, radial_order=radial_order, phase_mode=phase_mode, reason=reason, control=control)
    return binding


def construct(
    bindings: Mapping[str, Any],
    branch: Any = None,
    field: Any = None,
    Zorder: Any = None,
    radial_order: Any = 1,
    *,
    phase_mode: Any = "ordinary_Z",
    query: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    """Construct the saved typed quotient for one exact binding query."""

    return lookup(
        bindings,
        branch,
        field,
        Zorder,
        radial_order,
        phase_mode=phase_mode,
        query=query,
    )


class _Journal:
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
                return
            append = getattr(self.callback, "append", None)
            if callable(append):
                append(_clone(item))
                return
            raise TypeError("emit must be callable or appendable")
        except BaseException as exc:
            raise _CallbackFailure(event, exc) from exc


def _event_rows(journal: _Journal, event: str) -> list[dict[str, Any]]:
    return [
        _clone(item["row"])
        for item in journal.events
        if item.get("event") == event and _is_mapping(item.get("row"))
    ]


def _counts(journal: _Journal) -> dict[str, int]:
    event_names = [item.get("event") for item in journal.events]
    return {
        "new_algebra_family_calls": event_names.count("family_entered"),
        "new_compiled_operator_bindings": event_names.count("compiled_row"),
        "unavailable_source_rows": event_names.count("unavailable_source"),
        "actual_negative_api_refusals": event_names.count("api_refusal"),
        "old_math_calls": 0,
        "old_provider_calls": 0,
        "numeric_evaluation_calls": 0,
        "N_materialization_calls": 0,
        "source_value_identity_calls": 0,
        "new_mathematical_proofs": 0,
    }


def _result(
    *,
    journal: _Journal,
    passed: bool,
    status: str,
    failure: Mapping[str, Any] | None,
    source_validation: Mapping[str, Any] | None,
    parent_certificate: Mapping[str, Any] | None,
    bindings: Mapping[str, Any] | None,
    radius_ref: Mapping[str, Any] | None,
    cached_RZ: Mapping[str, Any] | None,
    background_section: Mapping[str, Any] | None,
    modified_section: Mapping[str, Any] | None,
    fastphase: Mapping[str, Any] | None,
) -> dict[str, Any]:
    flags = _flags(passed)
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
        "actual_counts": _counts(journal),
        "callback_event_count": len(journal.events),
        "callback_events": _clone(journal.events),
        "new_flags": flags,
        **flags,
        "payload": {
            "source_validation": _clone(source_validation) if source_validation is not None else None,
            "parent_certificate": _clone(parent_certificate) if parent_certificate is not None else {
                "scope": "typed declared-fixed-native-parameter semantics only",
                "cached_data_only": True,
                "new_proof": False,
            },
            "compiled_operator_bindings": _event_rows(journal, "compiled_row"),
            "unavailable_source_rows": _event_rows(journal, "unavailable_source"),
            "actual_negative_api_refusals": _event_rows(journal, "api_refusal"),
            "bindings": _clone(bindings) if bindings is not None else None,
            "radius": _clone(radius_ref) if radius_ref is not None else None,
            "cached_R_Z": _clone(cached_RZ) if cached_RZ is not None else None,
            "background_source_definition": _clone(background_section) if background_section is not None else None,
            "modified_source_definition": _clone(modified_section) if modified_section is not None else None,
            "modified_fastphase_provenance": _clone(fastphase) if fastphase is not None else None,
            "source_operator_scope": {
                "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
                "background_fields": list(FIELDS),
                "modified_fields": list(MODIFIED_FIELDS),
                "missing_modified_fields": list(MISSING_MODIFIED_FIELDS),
                "background_order_bound": 7,
                "modified_order_bound": 6,
                "radial_order": 1,
                "held_phase": False,
                "fractional_part_differentiated": False,
                "N_materialized": False,
                "source_value_identity_calls": 0,
                "derivative_identity_proved": False,
                "full_mixed_field_rows_complete": False,
                "physical_time_velocity_cartesian_claims": False,
                "log_measure_as_3D_volume": False,
            },
        },
    }


def derive(source: dict[str, Any], emit: Callable[[dict[str, Any]], Any] | Any = None) -> dict[str, Any]:
    """Compile saved first-radial/ordinary-Z source quotient expressions once.

    ``emit`` receives one ``{"event": ..., "row": ...}`` mapping per
    journalled event.  A callback failure stops later work and returns the
    complete event prefix with all new flags false.
    """

    journal = _Journal(emit)
    source_validation: Mapping[str, Any] | None = None
    parent_certificate: Mapping[str, Any] | None = source.get("parent_certificate") if _is_mapping(source) else None
    bindings: dict[str, Any] | None = None
    radius_ref: Mapping[str, Any] | None = None
    cached_RZ: Mapping[str, Any] | None = None
    background_section: Mapping[str, Any] | None = None
    modified_section: Mapping[str, Any] | None = None
    fastphase: Mapping[str, Any] | None = None
    try:
        journal.emit(
            "family_entered",
            {
                "family": "native_O3_mixed_radial_source_operator_v1",
                "source_only": True,
                "conditional": True,
                "coordinate_frame": EXPECTED_CHART,
                "cached_R_Z": "0",
                "regular_positiveR_premise": True,
                "background_order_bound": 7,
                "modified_order_bound": 6,
                "radial_order": 1,
                "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
                "native_field_derivative_identity_proved": False,
                "full_mixed_field_rows_complete": False,
                "fractional_part_differentiated": False,
                "N_materialized": False,
            },
        )
        _require(_is_mapping(source), "input_not_mapping", "source must be a mapping")
        literals, source_validation = _wrapper_literals(source)
        saved_result = literals["saved_Z_result"]
        background_packet = literals["background_literals"]
        modified_packet = literals["modified_literals"]
        background_rows, background_section = _background_rowsets(background_packet)
        modified_rows, modified_section, fastphase = _modified_rowsets(modified_packet)
        radius_ref = _radius_ref(background_packet)
        cached_RZ = _cached_RZ_ref(saved_result)

        bindings = {"background": {}, "modified": {}}
        for field in FIELDS:
            bindings["background"][field] = {}
            for order in ORDINARY_Z_ORDERS:
                row = _binding_row(
                    background_rows[field][order], radius_ref, cached_RZ,
                    branch="background", field=field, order=order, fastphase=None,
                )
                bindings["background"][field][order] = row
                journal.emit("compiled_row", construct(bindings, "background", field, order, 1))
        for field in MODIFIED_FIELDS:
            bindings["modified"][field] = {}
            for order in ORDINARY_Z_ORDERS:
                row = _binding_row(
                    modified_rows[field][order], radius_ref, cached_RZ,
                    branch="modified", field=field, order=order, fastphase=fastphase,
                )
                bindings["modified"][field][order] = row
                journal.emit("compiled_row", construct(bindings, "modified", field, order, 1))

        for field in MISSING_MODIFIED_FIELDS:
            try:
                lookup(bindings, "modified", field, 0, 1)
            except _APIRefusal as refusal:
                unavailable = dict(refusal.row)
                unavailable.update(
                    {
                        "decision": "unavailable",
                        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
                        "source_definition_status": "no saved total helper output for this field",
                        "reason": "modified total-y source rows are absent",
                        "background_substitution": False,
                        "padding_added": False,
                        "factorial_reapplied": False,
                        "native_field_derivative_identity_proved": False,
                        "full_mixed_field_rows_complete": False,
                    }
                )
                journal.emit("unavailable_source", unavailable)
            else:
                raise _KernelFailure("missing_source_admitted", f"missing modified field {field} was admitted")

        controls = (
            ("held_phase", {"branch": "background", "field": "E", "Zorder": 0, "radial_order": 1, "phase_mode": "held_phase"}),
            ("second_radial", {"branch": "background", "field": "E", "Zorder": 0, "radial_order": 2, "phase_mode": "ordinary_Z"}),
            ("wrong_branch", {"branch": "wrongbranch", "field": "E", "Zorder": 0, "radial_order": 1, "phase_mode": "ordinary_Z"}),
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

        journal.emit(
            "family_completed",
            {
                "family": "native_O3_mixed_radial_source_operator_v1",
                "compiled_rows": 60,
                "unavailable_source_rows": 6,
                "negative_api_refusals": 4,
                "new_mathematical_proofs": 0,
                "source_only": True,
                "conditional": True,
                "native_field_derivative_identity_proved": False,
                "full_mixed_field_rows_complete": False,
            },
        )
        return _result(
            journal=journal,
            passed=True,
            status="complete",
            failure=None,
            source_validation=source_validation,
            parent_certificate=parent_certificate,
            bindings=bindings,
            radius_ref=radius_ref,
            cached_RZ=cached_RZ,
            background_section=background_section,
            modified_section=modified_section,
            fastphase=fastphase,
        )
    except _CallbackFailure as exc:
        return _result(
            journal=journal,
            passed=False,
            status="failed",
            failure={"code": "diagnostic_callback_failure", "message": str(exc.cause), "event": exc.event, "prefix_preserved": True},
            source_validation=source_validation,
            parent_certificate=parent_certificate,
            bindings=bindings,
            radius_ref=radius_ref,
            cached_RZ=cached_RZ,
            background_section=background_section,
            modified_section=modified_section,
            fastphase=fastphase,
        )
    except _APIRefusal as exc:
        return _result(
            journal=journal,
            passed=False,
            status="failed",
            failure={"code": "unexpected_api_refusal", "message": str(exc), "row": _clone(exc.row), "prefix_preserved": True},
            source_validation=source_validation,
            parent_certificate=parent_certificate,
            bindings=bindings,
            radius_ref=radius_ref,
            cached_RZ=cached_RZ,
            background_section=background_section,
            modified_section=modified_section,
            fastphase=fastphase,
        )
    except _KernelFailure as exc:
        return _result(
            journal=journal,
            passed=False,
            status="failed",
            failure={"code": exc.code, "message": exc.message, "details": _clone(exc.details), "prefix_preserved": True},
            source_validation=source_validation,
            parent_certificate=parent_certificate,
            bindings=bindings,
            radius_ref=radius_ref,
            cached_RZ=cached_RZ,
            background_section=background_section,
            modified_section=modified_section,
            fastphase=fastphase,
        )
    except BaseException as exc:
        return _result(
            journal=journal,
            passed=False,
            status="failed",
            failure={"code": "unexpected_exception", "message": str(exc), "exception_type": type(exc).__name__, "prefix_preserved": True},
            source_validation=source_validation,
            parent_certificate=parent_certificate,
            bindings=bindings,
            radius_ref=radius_ref,
            cached_RZ=cached_RZ,
            background_section=background_section,
            modified_section=modified_section,
            fastphase=fastphase,
        )


__all__ = [
    "COMPONENT",
    "INPUT_SCHEMA",
    "RESULT_SCHEMA",
    "EXPECTED_SOURCE_SHA256",
    "derive",
    "lookup",
    "construct",
]
