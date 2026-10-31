"""Bounded, record-only join for the saved Q4/core inlet moments.

The caller authenticates the four dictionaries passed to :func:`bridge`.  This
module only compares rows already present in those dictionaries.  In
particular, it does not import a producer, walk a source implementation,
evaluate a graph node, perform an integral, or turn a directed enclosure into
a selected value.

The graph snapshot stores a Jet as a list of integer node handles.  A handle
is treated as an index local to its snapshot; matching rows are evidence at
the record layer and are deliberately not promoted to a dynamic binder or
source-function identity.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


SCHEMA = "original_q4_saved_moment_join_v1"
MOMENT_ORDERS = {
    "M": (0, 1, 2),
    "A": (0, 1, 2),
    "B": (0, 1, 2, 3),
    "C": (0, 1, 2),
}
MOMENT_NAMES = tuple(MOMENT_ORDERS)
PHYSICAL_NAMES = tuple(name for name in MOMENT_NAMES)
FAMILY_FIELDS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)
RHO_NAME = "original_n1_true_core_rho"
MOMENT_VARIABLE_NAME = "original_core_inlet_moment_t"
FACTORIALS = {0: 1, 1: 1, 2: 2, 3: 6}


class _Unresolved(ValueError):
    """Internal fail-closed signal; all public failures become JSON records."""


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise _Unresolved(label + " must be a mapping")
    return value


def _integer(value: Any, label: str) -> int:
    if type(value) is not int:
        raise _Unresolved(label + " must be an integer")
    return value


def _node(nodes: Any, handle: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(nodes, list):
        raise _Unresolved("graph node table missing")
    if type(handle) is not int or handle < 0 or handle >= len(nodes):
        raise _Unresolved(label + " is not a valid graph handle")
    row = nodes[handle]
    if not isinstance(row, Mapping):
        raise _Unresolved(label + " row is not a mapping")
    return row


def _lookup(container: Any, key: Any, label: str) -> Any:
    """Read a JSON mapping/list slot without coercing arbitrary values."""

    if isinstance(container, Mapping):
        if key in container:
            return container[key]
        string_key = str(key)
        if string_key in container:
            return container[string_key]
        raise _Unresolved(label + " is missing")
    if isinstance(container, list) and type(key) is int and 0 <= key < len(container):
        return container[key]
    raise _Unresolved(label + " is missing")


def _family(record: Any, label: str) -> dict[str, Any]:
    record = _mapping(record, label)
    value = record.get("source_family")
    if value is None and all(field in record for field in FAMILY_FIELDS):
        value = record
    value = _mapping(value, label + ".source_family")
    result: dict[str, Any] = {}
    for field in FAMILY_FIELDS:
        if field not in value:
            raise _Unresolved(label + ".source_family." + field + " is missing")
        result[field] = value[field]
    return result


def _nodes(record: Any, label: str) -> list[Any]:
    record = _mapping(record, label)
    nodes = record.get("exact_graph_nodes")
    if not isinstance(nodes, list):
        raise _Unresolved(label + ".exact_graph_nodes is missing")
    return nodes


def _is_exact_rational(row: Any, numerator: int, denominator: int = 1) -> bool:
    if not isinstance(row, Mapping) or row.get("operation") != "exact_rational":
        return False
    return row.get("numerator") == numerator and row.get("denominator") == denominator


def _endpoint_raw(row: Any, side: str) -> Any:
    """Extract an archive endpoint without decoding or rounding it."""

    if isinstance(row, Mapping):
        for key in (
            side + "_exact_mpf_tuple",
            side + "_exact_tuple",
            side + "_endpoint",
            side,
        ):
            if key in row:
                return row[key]
        nested = row.get("interval")
        if nested is not None:
            return _endpoint_raw(nested, side)
        return None
    if isinstance(row, (list, tuple)) and len(row) == 2:
        # An interval is a pair of endpoints.  A two-integer signed dyadic is
        # itself an endpoint and must not be split.
        if all(type(item) is int for item in row):
            return list(row)
        return row[0 if side == "lower" else 1]
    return None


def _error_upper_raw(row: Any) -> Any:
    if isinstance(row, Mapping):
        for key in ("upper_exact_mpf_tuple", "error_upper", "upper_exact_tuple", "upper"):
            if key in row:
                value = row[key]
                if key in ("error_upper",) and isinstance(value, Mapping):
                    return _error_upper_raw(value)
                return value
        nested = row.get("interval")
        if nested is not None:
            return _error_upper_raw(nested)
        return None
    if isinstance(row, (list, tuple)) and len(row) == 2 and not all(type(item) is int for item in row):
        return row[1]
    return None


def _listify(value: Any) -> Any:
    if isinstance(value, tuple):
        return [_listify(item) for item in value]
    if isinstance(value, list):
        return [_listify(item) for item in value]
    if isinstance(value, Mapping):
        return {key: _listify(item) for key, item in value.items()}
    return value


def _endpoint_integer(value: Any) -> int | None:
    """Recognise only small exact integer metadata used for rho/Z guards."""

    if type(value) is int:
        return value
    if isinstance(value, str):
        if value in {"-1", "1", "4"}:
            return int(value)
        return None
    if isinstance(value, Mapping):
        if value.get("operation") == "exact_rational":
            numerator, denominator = value.get("numerator"), value.get("denominator")
            if type(numerator) is int and type(denominator) is int and denominator:
                if numerator % denominator == 0:
                    return numerator // denominator
            return None
        for key in ("value", "exact_value", "integer"):
            if key in value:
                return _endpoint_integer(value[key])
        return None
    if not isinstance(value, (list, tuple)):
        return None
    if len(value) == 4 and all(type(item) is int for item in value):
        sign, mantissa, exponent, bitcount = value
        if sign not in (0, 1) or mantissa < 0 or bitcount < 0:
            return None
        if mantissa == 0:
            return 0 if sign == 0 else None
        if bitcount != mantissa.bit_length():
            return None
        if exponent >= 0:
            # This helper is only a guard for the exact metadata values
            # -1, 1, and 4.  Do not expand a hostile astronomical exponent.
            if exponent > 2:
                return None
            magnitude = mantissa << exponent
            return -magnitude if sign else magnitude
        if -exponent > bitcount:
            return None
        divisor = 1 << (-exponent)
        if mantissa % divisor:
            return None
        magnitude = mantissa // divisor
        return -magnitude if sign else magnitude
    if len(value) == 2 and all(type(item) is int for item in value):
        signed_mantissa, exponent = value
        if exponent >= 0:
            if exponent > 2:
                return None
            return signed_mantissa << exponent
        if -exponent > abs(signed_mantissa).bit_length():
            return None
        divisor = 1 << (-exponent)
        if signed_mantissa % divisor:
            return None
        return signed_mantissa // divisor
    return None


def _metadata_endpoints(row: Any, label: str) -> tuple[Any, Any]:
    if isinstance(row, Mapping):
        lower = _endpoint_raw(row, "lower")
        upper = _endpoint_raw(row, "upper")
        if lower is not None and upper is not None:
            return lower, upper
        # A few JSON snapshots use ``bounds`` around the same two raw tuples.
        bounds = row.get("bounds")
        if bounds is not None:
            return _metadata_endpoints(bounds, label + ".bounds")
    if isinstance(row, (list, tuple)) and len(row) == 2 and not all(type(item) is int for item in row):
        return row[0], row[1]
    raise _Unresolved(label + " directed endpoints are missing")


def _require_metadata(view: Mapping[str, Any], checks: dict[str, bool]) -> None:
    rho = view.get("rho")
    z = view.get("Z")
    rho_lower, rho_upper = _metadata_endpoints(rho, "inlet.rho")
    z_lower, z_upper = _metadata_endpoints(z, "inlet.Z")
    checks["saved_inlet_rho_lower_is_4"] = _endpoint_integer(rho_lower) == 4
    checks["saved_inlet_rho_upper_is_4"] = _endpoint_integer(rho_upper) == 4
    checks["saved_inlet_Z_lower_is_minus_1"] = _endpoint_integer(z_lower) == -1
    checks["saved_inlet_Z_upper_is_1"] = _endpoint_integer(z_upper) == 1
    if not all(checks[key] for key in (
        "saved_inlet_rho_lower_is_4",
        "saved_inlet_rho_upper_is_4",
        "saved_inlet_Z_lower_is_minus_1",
        "saved_inlet_Z_upper_is_1",
    )):
        raise _Unresolved("saved inlet rho/Z metadata is not rho=4, Z=[-1,1]")


def _view(interior: Any) -> tuple[Mapping[str, Any], dict[str, Any]]:
    interior = _mapping(interior, "core_interior")
    views = interior.get("current_core_interior_views", interior.get("views"))
    views = _mapping(views, "core_interior.views")
    inlet = _mapping(views.get("inlet"), "core_interior.views.inlet")
    return inlet, _family(interior, "core_interior")


def _grid(view: Mapping[str, Any], interior: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = view.get(name)
    if value is None:
        value = interior.get(name)
    return _mapping(value, "saved inlet " + name)


def _saved_rows(
    view: Mapping[str, Any],
    interior: Mapping[str, Any],
    name: str,
    order: int,
) -> dict[str, Any]:
    key = "rho0_Z" + str(order)
    finite_grid = _grid(view, interior, "finite_integrated_axial_Taylor_coefficient_grids")
    tail_grid = _grid(view, interior, "directed_integrated_axial_Taylor_tail_bounds")
    ordinary_grid = _grid(view, interior, "ordinary_rho4_axial6_moment_grids")
    finite_by_name = _mapping(finite_grid.get(name), "saved finite " + name)
    tail_by_name = _mapping(tail_grid.get(name), "saved tail " + name)
    ordinary_by_name = _mapping(ordinary_grid.get(name), "saved ordinary " + name)
    finite = _lookup(finite_by_name, key, "saved finite " + name + "." + key)
    tail = _lookup(tail_by_name, key, "saved tail " + name + "." + key)
    ordinary = _lookup(ordinary_by_name, key, "saved ordinary " + name + "." + key)
    if finite is None or tail is None or ordinary is None:
        raise _Unresolved("saved inlet row is null for " + name + "." + key)
    if _endpoint_raw(finite, "lower") is None or _endpoint_raw(finite, "upper") is None:
        raise _Unresolved("saved finite row has no raw endpoints for " + name + "." + key)
    if _error_upper_raw(tail) is None:
        raise _Unresolved("saved tail row has no upper endpoint for " + name + "." + key)
    return {
        "finite_integrated_axial_Taylor_coefficient_row": finite,
        "directed_integrated_axial_Taylor_tail_bound_row": tail,
        "ordinary_rho4_axial6_moment_enclosure_row": ordinary,
        "factorial": FACTORIALS[order],
        "range_is_value": False,
    }


def _prior_item(item: Any) -> tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any] | None]:
    item = _mapping(item, "prior inspection")
    inspection = item.get("inspection", item)
    inspection = _mapping(inspection, "prior inspection.result")
    evidence = _mapping(inspection.get("evidence"), "prior inspection.evidence")
    scope_value = inspection.get("scope", {})
    scope = _mapping(scope_value, "prior inspection.scope")
    frontier_value = item.get("frontier")
    frontier = _mapping(frontier_value, "prior inspection.frontier") if frontier_value is not None else None
    return item, inspection, evidence, frontier


def _check_prior_scope(scope: Mapping[str, Any], label: str) -> None:
    # These flags are optional across the historical inspection receipts.  If
    # present, a true value would overstate what the prior inspection proved.
    for key in (
        "capture_safe_substitution_proved",
        "capture_avoidance_proved",
        "numeric_function_value_computed",
        "numeric_field_value_computed",
        "old_errors_removed",
        "terminal_closure_proved",
    ):
        if key in scope and scope[key] is not False:
            raise _Unresolved(label + " prior scope flag is not false: " + key)
    if scope.get("opaque_frontier_present") is True or scope.get("opaque_scope") is True:
        raise _Unresolved(label + " has an opaque prior scope")


def _physical_view(physical_result: Any, family: Mapping[str, Any]) -> Mapping[str, Any]:
    physical_result = _mapping(physical_result, "physical_result")
    if physical_result.get("status") not in (None, "returned"):
        raise _Unresolved("physical_result is not a returned v2 receipt")
    # The v2 receipt records the family under family_evidence, rather than a
    # top-level source_family field.  Do not invent a physical-family field.
    family_evidence = _mapping(physical_result.get("family_evidence"), "physical_result.family_evidence")
    saved_family = _mapping(family_evidence.get("saved_core_family"),
                            "physical_result.family_evidence.saved_core_family")
    for field in ("all_three_Q4_core_family_records_match", "fixed_point_family_matches"):
        if field in family_evidence and family_evidence[field] is not True:
            raise _Unresolved("physical family evidence flag is not true: " + field)
    physical_family = {field: saved_family.get(field) for field in FAMILY_FIELDS}
    if any(value is None for value in physical_family.values()):
        raise _Unresolved("physical_result.family_evidence.saved_core_family is incomplete")
    if physical_family != dict(family):
        raise _Unresolved("physical_result source family mismatch")
    views = _mapping(physical_result.get("views"), "physical_result.views")
    return _mapping(views.get("inlet"), "physical_result.views.inlet")


def _physical_rows(
    physical_result: Any,
    saved: dict[str, dict[str, dict[str, Any]]],
    family: Mapping[str, Any],
) -> dict[str, dict[str, dict[str, Any]]]:
    inlet = _physical_view(physical_result, family)
    moments = _mapping(inlet.get("input_moments"), "physical_result.views.inlet.input_moments")
    output: dict[str, dict[str, dict[str, Any]]] = {}
    for name in PHYSICAL_NAMES:
        raw_name = moments.get(name, moments.get(name.lower()))
        raw_name = _mapping(raw_name, "physical input moment " + name)
        output[name] = {}
        for order in (0, 1):
            zone_key = "z" + str(order)
            zone = _mapping(raw_name.get(zone_key), "physical input moment " + name + "." + zone_key)
            finite_actual = zone.get("finite")
            error_actual = zone.get("error_upper")
            if finite_actual is None or error_actual is None:
                raise _Unresolved("physical input row is missing finite/error for " + name + "." + zone_key)
            saved_row = saved[name][str(order)]
            finite_saved = saved_row["finite_integrated_axial_Taylor_coefficient_row"]
            tail_saved = saved_row["directed_integrated_axial_Taylor_tail_bound_row"]
            expected_finite = [
                _listify(_endpoint_raw(finite_saved, "lower")),
                _listify(_endpoint_raw(finite_saved, "upper")),
            ]
            expected_error = _listify(_error_upper_raw(tail_saved))
            actual_finite = _listify(finite_actual)
            if isinstance(finite_actual, Mapping):
                actual_finite = [
                    _listify(_endpoint_raw(finite_actual, "lower")),
                    _listify(_endpoint_raw(finite_actual, "upper")),
                ]
            actual_error = _listify(error_actual)
            if isinstance(error_actual, Mapping):
                actual_error = _listify(_error_upper_raw(error_actual))
            if expected_finite[0] is None or expected_finite[1] is None or expected_error is None:
                raise _Unresolved("saved physical comparison endpoint is missing")
            if actual_finite != expected_finite:
                raise _Unresolved("physical finite endpoints differ for " + name + "." + zone_key)
            if actual_error != expected_error:
                raise _Unresolved("physical tail upper endpoint differs for " + name + "." + zone_key)
            output[name][zone_key] = {
                "v2_finite": finite_actual,
                "v2_error_upper": error_actual,
                "saved_finite_lower_exact_mpf_tuple": expected_finite[0],
                "saved_finite_upper_exact_mpf_tuple": expected_finite[1],
                "saved_tail_upper_exact_mpf_tuple": expected_error,
                "exact_raw_tuple_match": True,
                "range_is_value": False,
            }
    return output


def bridge(
    core: Any,
    q4: Any,
    inspections: Any,
    interior: Any,
    physical_result: Any = None,
) -> dict[str, Any]:
    """Join thirteen already-authenticated source rows and saved ranges.

    The optional fifth argument is intentionally fail-closed when omitted:
    this revision's physical-scale contract requires the v2 inlet receipt.
    """

    checks: dict[str, bool] = {}
    errors: list[str] = []
    evidence: dict[str, Any] = {}
    scope: dict[str, Any] = {
        "capture_safe_substitution_proved": False,
        "capture_avoidance_proved": False,
        "dynamic_binder_identity_proved": False,
        "source_evaluator_invoked": False,
        "numeric_function_value_computed": False,
        "numeric_field_value_computed": False,
        "old_errors_removed": False,
        "terminal_closure_proved": False,
        "q4_equalities_and_intervals_only_conditional": True,
        "saved_enclosures_are_selected_values": False,
        "factorial_applied_to_saved_rows": False,
        "matching_rows_are_record_layer_evidence_only": True,
    }

    def require(label: str, condition: bool, message: str | None = None) -> None:
        checks[label] = bool(condition)
        if not condition:
            raise _Unresolved(message or label)

    try:
        core_record = _mapping(core, "core")
        q4_record = _mapping(q4, "q4")
        core_family = _family(core_record, "core")
        q4_family = _family(q4_record, "q4")
        require("q4_core_source_family_matches", q4_family == core_family, "q4/core source family differs")
        inlet_view, interior_family = _view(interior)
        require("core_interior_source_family_matches", interior_family == core_family, "core interior source family differs")
        _require_metadata(inlet_view, checks)

        q4_nodes = _nodes(q4_record, "q4")
        core_nodes = _nodes(core_record, "core")
        core_bindings = _mapping(core_record.get("core_inlet_axial6_bindings"), "core.core_inlet_axial6_bindings")
        moments = _mapping(core_bindings.get("moments"), "core moments")
        inlet_moments = _mapping(core_bindings.get("inlet_moments"), "core inlet_moments")
        require("selected_moment_families_are_MABC", all(name in moments and name in inlet_moments for name in MOMENT_NAMES))

        if not isinstance(inspections, list):
            raise _Unresolved("prior inspections must be a list")
        require("prior_inspection_count_is_13", len(inspections) == 13)

        by_integral: dict[int, tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any] | None]] = {}
        for item in inspections:
            parsed = _prior_item(item)
            evidence_row = parsed[2]
            integral_handle = _integer(evidence_row.get("integral_node"), "prior integral_node")
            if integral_handle in by_integral:
                raise _Unresolved("duplicate prior integral_node")
            by_integral[integral_handle] = parsed
        require("prior_integral_nodes_are_unique", len(by_integral) == 13)

        joined: list[dict[str, Any]] = []
        saved_rows: dict[str, dict[str, dict[str, Any]]] = {}
        for name in MOMENT_NAMES:
            saved_rows[name] = {}
            for order in MOMENT_ORDERS[name]:
                slot = _lookup(moments, name, "core moments " + name)
                moment_handle = _integer(_lookup(slot, order, "core moments " + name + "." + str(order)),
                                         "core moments " + name + "." + str(order))
                inlet_slot = _lookup(inlet_moments, name, "core inlet_moments " + name)
                inlet_handle = _integer(_lookup(inlet_slot, order, "core inlet_moments " + name + "." + str(order)),
                                        "core inlet_moments " + name + "." + str(order))
                q4_moment = _node(q4_nodes, moment_handle, "q4 moment " + name + "." + str(order))
                core_moment = _node(core_nodes, moment_handle, "core moment " + name + "." + str(order))
                q4_inlet = _node(q4_nodes, inlet_handle, "q4 inlet " + name + "." + str(order))
                core_inlet = _node(core_nodes, inlet_handle, "core inlet " + name + "." + str(order))
                require("row_match_moment_" + name + "_Z" + str(order), q4_moment == core_moment,
                        "q4/core moment row differs for " + name + "." + str(order))
                require("row_match_inlet_" + name + "_Z" + str(order), q4_inlet == core_inlet,
                        "q4/core inlet substitution row differs for " + name + "." + str(order))
                require("outer_substitution_operation_" + name + "_Z" + str(order),
                        q4_inlet.get("operation") == "function_substitution",
                        "outer substitution operation differs for " + name + "." + str(order))
                require("outer_expression_" + name + "_Z" + str(order),
                        q4_inlet.get("expression") == moment_handle,
                        "outer substitution expression differs for " + name + "." + str(order))
                variable_handle = q4_inlet.get("variable")
                require("outer_variable_handle_" + name + "_Z" + str(order),
                        type(variable_handle) is int,
                        "outer substitution variable is not a graph handle for " + name + "." + str(order))
                variable_row = _node(q4_nodes, variable_handle, "outer substitution variable")
                require("outer_variable_" + name + "_Z" + str(order),
                        variable_row == RHO_NAME or (
                            isinstance(variable_row, Mapping)
                            and variable_row.get("operation") == "bound_variable"
                            and variable_row.get("name") == RHO_NAME
                        ),
                        "outer substitution variable differs for " + name + "." + str(order))
                value_handle = q4_inlet.get("value")
                require("outer_value_handle_" + name + "_Z" + str(order),
                        type(value_handle) is int,
                        "outer substitution value is not a graph handle for " + name + "." + str(order))
                value_row = _node(q4_nodes, value_handle, "outer substitution value")
                require("outer_value_" + name + "_Z" + str(order), _is_exact_rational(value_row, 4),
                        "outer substitution value differs for " + name + "." + str(order))
                require("moment_operation_" + name + "_Z" + str(order),
                        q4_moment.get("operation") == "definite_integral",
                        "moment operation differs for " + name + "." + str(order))
                require("moment_name_" + name + "_Z" + str(order),
                        q4_moment.get("actual_core_normalized_moment") == name,
                        "moment name differs for " + name + "." + str(order))
                require("moment_order_" + name + "_Z" + str(order),
                        q4_moment.get("ordinary_slow_Z_derivative_order") == order,
                        "moment ordinary order differs for " + name + "." + str(order))
                require("moment_variable_" + name + "_Z" + str(order),
                        q4_moment.get("variable") == MOMENT_VARIABLE_NAME,
                        "moment variable differs for " + name + "." + str(order))
                require("moment_unit_lower_" + name + "_Z" + str(order),
                        _is_exact_rational(_node(q4_nodes, q4_moment.get("lower"), "moment lower"), 0),
                        "moment lower endpoint differs for " + name + "." + str(order))
                require("moment_unit_upper_" + name + "_Z" + str(order),
                        _is_exact_rational(_node(q4_nodes, q4_moment.get("upper"), "moment upper"), 1),
                        "moment upper endpoint differs for " + name + "." + str(order))

                prior = by_integral.get(moment_handle)
                if prior is None:
                    raise _Unresolved("prior inspection missing integral " + str(moment_handle))
                item, inspection, prior_evidence, frontier = prior
                require("prior_observed_" + name + "_Z" + str(order), inspection.get("status") == "observed",
                        "prior inspection is not observed for " + name + "." + str(order))
                _check_prior_scope(_mapping(inspection.get("scope", {}), "prior scope"),
                                   name + ".Z" + str(order))
                require("prior_integrand_link_" + name + "_Z" + str(order),
                        prior_evidence.get("integral_node") == moment_handle
                        and prior_evidence.get("integrand_node") == q4_moment.get("integrand"),
                        "prior integral/integrand handles differ for " + name + "." + str(order))
                require("prior_integral_source_row_" + name + "_Z" + str(order),
                        prior_evidence.get("integral_source_row") == q4_moment,
                        "prior integral source row differs for " + name + "." + str(order))
                if frontier is not None:
                    require("prior_frontier_link_" + name + "_Z" + str(order),
                            frontier.get("node_id") == moment_handle
                            and frontier.get("integrand_node_id", frontier.get("integrand_node")) == q4_moment.get("integrand"),
                            "prior frontier differs for " + name + "." + str(order))
                if "binder_name" in prior_evidence:
                    require("prior_binder_name_" + name + "_Z" + str(order),
                            prior_evidence.get("binder_name") == MOMENT_VARIABLE_NAME,
                            "prior binder name differs for " + name + "." + str(order))
                if "unit_interval_bounds_verified" in prior_evidence:
                    require("prior_unit_interval_" + name + "_Z" + str(order),
                            prior_evidence.get("unit_interval_bounds_verified") is True,
                            "prior unit interval record differs for " + name + "." + str(order))

                saved = _saved_rows(inlet_view, _mapping(interior, "core interior"), name, order)
                saved_rows[name][str(order)] = saved
                joined.append({
                    "name": name,
                    "ordinary_Z_order": order,
                    "factorial": FACTORIALS[order],
                    "q4": {
                        "moment_node": moment_handle,
                        "inlet_node": inlet_handle,
                        "moment_row": q4_moment,
                        "inlet_substitution_row": q4_inlet,
                    },
                    "core": {
                        "moment_node": moment_handle,
                        "inlet_node": inlet_handle,
                        "moment_row": core_moment,
                        "inlet_substitution_row": core_inlet,
                    },
                    "prior_observation": {
                        "integral_node": prior_evidence.get("integral_node"),
                        "integrand_node": prior_evidence.get("integrand_node"),
                        "integral_source_row": prior_evidence.get("integral_source_row"),
                        "status": inspection.get("status"),
                    },
                    "saved_inlet_rows": saved,
                })

        require("all_13_moment_rows_joined", len(joined) == 13)
        require("all_prior_integrals_consumed", set(by_integral) == {
            item["prior_observation"]["integral_node"] for item in joined
        })
        physical = _physical_rows(physical_result, saved_rows, core_family)
        require("physical_v2_inlet_rows_joined", len(physical) == 4)
        scope["factorial_applied_to_saved_rows"] = False
        evidence = {
            "source_family": core_family,
            "ordinary_Z_orders": [0, 1, 2, 3],
            "moment_rows": joined,
            "physical_v2_inlet_rows": physical,
            "factorial_conversion": {
                name: {str(order): FACTORIALS[order] for order in orders}
                for name, orders in MOMENT_ORDERS.items()
            },
            "saved_rows_are_retained_raw": True,
            "finite_taylor_rows_are_not_values": True,
            "tail_rows_are_not_values": True,
            "ordinary_enclosure_rows_are_not_values": True,
        }
        checks["no_H_or_K_claimed_in_13_rows"] = all(item["name"] in MOMENT_NAMES for item in joined)
        require("no_H_or_K_claimed_in_13_rows", checks["no_H_or_K_claimed_in_13_rows"])
        return {
            "schema": SCHEMA,
            "status": "matched",
            "passed": True,
            "checks": checks,
            "evidence": evidence,
            "scope": scope,
            "errors": errors,
        }
    except _Unresolved as exc:
        errors.append(str(exc))
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        errors.append(type(exc).__name__ + ": " + str(exc))
    return {
        "schema": SCHEMA,
        "status": "unresolved",
        "passed": False,
        "checks": checks,
        "evidence": evidence,
        "scope": scope,
        "errors": errors,
    }


__all__ = ["SCHEMA", "MOMENT_ORDERS", "FACTORIALS", "bridge"]
