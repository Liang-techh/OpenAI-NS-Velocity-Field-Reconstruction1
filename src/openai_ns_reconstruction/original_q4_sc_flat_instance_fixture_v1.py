"""Fresh synthetic fixture for joining a new Q4 ``sc`` instance to flat Q4.

The fixture models only the JSON-like receipts consumed by the new join
adapter.  It does not import a producer, decode an archive, walk a graph, or
evaluate a field.  Every case builds new dictionaries, including its pinned
input-hash metadata, before applying a case-local mutation.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


SCHEMA = "original_q4_sc_flat_instance_fixture_v1"

# These are opaque raw-source identities.  They are the two source hashes
# which the existing receipts expose for the Q4 and prescribed-loop inputs.
Q4_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
)
LOOP_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.json.gz"
)
Q4_SHA = "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"
LOOP_SHA = "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"

# Small, opaque synthetic handles.  Their values have no graph meaning.
SC_NODE = 17
SC_OVER_TWO_NODE = 18
COORDINATE_NODE = 29
PHYSICAL_U_Y_ROWS = (901, 902, 903)

# 3 * 2**(-6) is a canonical positive dyadic below 1/8.  The half tuple is
# carried as a separate exact binding; no power-of-two calculation is used.
SC_EXACT_MPF_TUPLE = (0, 3, -6, 2)
SC_OVER_TWO_EXACT_MPF_TUPLE = (0, 3, -7, 2)
# The endpoint radial coordinate carried by a source port is the exact sc/2
# dyadic, not a separately invented coordinate value.
COORDINATE_EXACT_MPF_TUPLE = SC_OVER_TWO_EXACT_MPF_TUPLE
ORDINARY_Z_ORDERS = (0, 1, 2)


def _graph_id() -> dict[str, str]:
    """Return a fresh graph-identity mapping for one synthetic case."""

    return {
        "q4": Q4_SHA,
        "loop": LOOP_SHA,
        "q4_input_hash": Q4_SHA,
        "loop_input_hash": LOOP_SHA,
    }


def _instance_scope(graph_id: Mapping[str, str]) -> dict[str, Any]:
    """Return the deliberately narrow scope of the new explicit instance."""

    return {
        "newly_chosen_binding": True,
        "existing_archive_assignment_proved": False,
        "corrected_owner_identity_proved": False,
        "old_errors_removed": False,
        "flat_conclusion_proved_here": False,
        "binding_materialized": True,
        "graph_id": dict(graph_id),
        "prior_q4_flat_conclusion_used": False,
        "q4_flat_conclusion_used_only_conditionally": True,
        "numeric_values_computed": False,
        "source_graph_walked": False,
    }


def _instance_bundle() -> dict[str, Any]:
    """Build one fresh synthetic explicit-instance result bundle."""

    graph_id = _graph_id()
    scope = _instance_scope(graph_id)
    sc_tuple = list(SC_EXACT_MPF_TUPLE)
    half_tuple = list(SC_OVER_TWO_EXACT_MPF_TUPLE)
    selected_source_parameter = {
        "name": "sc",
        "node_id": SC_NODE,
        "sc_node_id": SC_NODE,
        "exact_mpf_tuple": list(sc_tuple),
        "sc_exact_mpf_tuple": list(sc_tuple),
        "graph_id": dict(graph_id),
        "selected_from_singleton_endpoint": True,
        "old_same_source_claimed": False,
    }
    half_record = {
        "name": "sc_over_two",
        "node_id": SC_OVER_TWO_NODE,
        "exact_mpf_tuple": list(half_tuple),
        "sc_over_two_exact_mpf_tuple": list(half_tuple),
        "power_of_two_materialized": False,
    }
    binding = {
        "binding_kind": "new_explicit_q4_sc_parameter_instance",
        "graph_id": dict(graph_id),
        "node_id": SC_NODE,
        "sc_node_id": SC_NODE,
        "sc_exact_mpf_tuple": list(sc_tuple),
        "selected_source_parameter": selected_source_parameter,
        "derived_dyadics": {
            "sc_over_two": half_record,
            "sc_half": dict(half_record),
            "comparison": {
                "sc_over_two_lt_three_sc_over_four": True,
                "three_sc_over_four_lt_sc": True,
                "power_of_two_materialized": False,
            },
        },
        "coordinate_node": COORDINATE_NODE,
        "coordinate_exact_mpf_tuple": list(COORDINATE_EXACT_MPF_TUPLE),
        "scope": dict(scope),
    }
    instance_result = {
        "status": "verified",
        "passed": True,
        "binding": binding,
        "selected_source_parameter": selected_source_parameter,
        "derived_dyadics": {
            "sc_over_two": dict(half_record),
            "sc_half": dict(half_record),
            "comparison": {
                "sc_over_two_lt_three_sc_over_four": True,
                "three_sc_over_four_lt_sc": True,
                "power_of_two_materialized": False,
            },
        },
        "scope": scope,
        "checks": {
            "new_binding_materialized": True,
            "sc_positive_proved_from_canonical_tuple": True,
            "derived_sc_over_two_is_canonical": True,
            "graph_identifier_explicit": True,
        },
    }
    return {
        "result": {
            "status": "returned",
            "instance_result": instance_result,
            "input_hashes": {
                Q4_ARCHIVE: Q4_SHA,
                LOOP_ARCHIVE: LOOP_SHA,
                "experiments/root_st073/synthetic_sc_instance_receipt.pin": "1" * 64,
                "experiments/root_st073/synthetic_q4_sc_instance_source.py": "2" * 64,
            },
            "receipt_pins": {
                "instance_receipt": {
                    "status": "verified",
                    "source_sha256": "3" * 64,
                    "input_hashes": {Q4_ARCHIVE: Q4_SHA, LOOP_ARCHIVE: LOOP_SHA},
                }
            },
            "raw_inputs_authenticated_before_decode": True,
            "source_code_executed": False,
            "source_graph_walked": False,
            "source_values_evaluated": False,
        }
    }


def _flat_scope() -> dict[str, Any]:
    return {
        "Q4_Uy_correction_zero_at_endpoint_conditional_on_sc_positive": True,
        "sc_positive_independently_certified": False,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "corrected_source_identity_proved": False,
        "history_function_equality_claimed": False,
        "old_errors_removed": False,
        "numeric_values_computed": False,
    }


def _flat_bundle() -> dict[str, Any]:
    """Build one fresh synthetic prior-flat result bundle."""

    source_hashes = {"q4": Q4_SHA, "loop": LOOP_SHA}
    coordinate_tuple = list(COORDINATE_EXACT_MPF_TUPLE)
    selected_tuple = list(SC_EXACT_MPF_TUPLE)
    half_tuple = list(SC_OVER_TWO_EXACT_MPF_TUPLE)
    evidence = {
        "physical_U_y_rows": list(PHYSICAL_U_Y_ROWS),
        "background_V_y_rows": [611, 612, 613],
        "Pstar_node_from_prior_receipt": 707,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "coordinate": COORDINATE_NODE,
        "coordinate_node": COORDINATE_NODE,
        "coordinate_exact_mpf_tuple": coordinate_tuple,
        "selected_sc_node": SC_NODE,
        "selected_sc_exact_mpf_tuple": selected_tuple,
        "sc_over_two_node": SC_OVER_TWO_NODE,
        "sc_over_two_exact_mpf_tuple": half_tuple,
        "source_hashes": dict(source_hashes),
    }
    archive_result = {
        "status": "verified",
        "passed": True,
        "checks": {
            "prior_flat_is_conditional": True,
            "prior_physical_receipt_binds_actual_Uy_rows": True,
            "keeps_three_ordinary_rows": True,
        },
        "evidence": evidence,
        "source_hashes": dict(source_hashes),
        "scope": _flat_scope(),
    }
    return {
        "result": {
            "status": "returned",
            "archive_result": archive_result,
            "input_hashes": {
                Q4_ARCHIVE: Q4_SHA,
                LOOP_ARCHIVE: LOOP_SHA,
                "experiments/root_st073/synthetic_flat_receipt.pin": "4" * 64,
                "experiments/root_st073/synthetic_physical_rows_receipt.pin": "5" * 64,
            },
            "receipt_pins": {
                "flat_receipt": {
                    "status": "verified",
                    "source_sha256": "6" * 64,
                    "input_hashes": {Q4_ARCHIVE: Q4_SHA, LOOP_ARCHIVE: LOOP_SHA},
                }
            },
            "raw_inputs_authenticated_before_decode": True,
            "source_code_executed": False,
            "source_graph_walked": False,
            "source_values_evaluated": False,
        }
    }


def _new_pair() -> tuple[dict[str, Any], dict[str, Any]]:
    """Create independent instance and flat dictionaries for one case."""

    return _instance_bundle(), _flat_bundle()


def _safe_get(value: Any, *keys: str) -> Any:
    for key in keys:
        if not isinstance(value, Mapping) or key not in value:
            return None
        value = value[key]
    return value


def _result_rejected(value: Any) -> bool:
    return (
        isinstance(value, Mapping)
        and value.get("passed") is False
        and value.get("status") in ("mismatch", "unresolved")
    )


def _positive_scope_is_narrow(value: Any) -> bool:
    scope = _safe_get(value, "scope")
    if not isinstance(scope, Mapping):
        return False
    new_binding = (
        scope.get("new_instance_only") is True
        or scope.get("newly_chosen_binding") is True
    )
    prohibited_claims = (
        "existing_archive_assignment_proved",
        "corrected_owner_identity_proved",
        "old_errors_removed",
        "flat_conclusion_proved_here",
        "cross_owner_identity_proved",
        "global_closure_proved",
    )
    return new_binding and all(
        key not in scope or scope.get(key) is False for key in prohibited_claims
    )


def _tuple_is(value: Any, expected: tuple[int, int, int, int]) -> bool:
    return value == list(expected) or value == expected


def _selected_sc_binding(row: Mapping[str, Any]) -> bool:
    direct = (
        row.get("selected_sc_node") == SC_NODE
        and _tuple_is(row.get("selected_sc_exact_mpf_tuple"), SC_EXACT_MPF_TUPLE)
    )
    aliases = (
        row.get("sc_node_id") == SC_NODE
        and _tuple_is(row.get("sc_exact_mpf_tuple"), SC_EXACT_MPF_TUPLE)
    )
    parameter = (
        row.get("parameter_node") == SC_NODE
        and _tuple_is(row.get("parameter_exact_mpf_tuple"), SC_EXACT_MPF_TUPLE)
    )
    nested = row.get("parameter_binding")
    if isinstance(nested, Mapping):
        selected = nested.get("selected_sc")
        if isinstance(selected, Mapping):
            direct = direct or (
                selected.get("node_id") == SC_NODE
                and _tuple_is(selected.get("exact_mpf_tuple"), SC_EXACT_MPF_TUPLE)
            )
    return direct or aliases or parameter


def _sc_over_two_binding(row: Mapping[str, Any]) -> bool:
    direct_pairs = (
        ("sc_over_two_node", "sc_over_two_exact_mpf_tuple"),
        ("sc_half_node", "sc_half_exact_mpf_tuple"),
        ("parameter_half_node", "parameter_half_exact_mpf_tuple"),
    )
    for node_key, tuple_key in direct_pairs:
        if row.get(node_key) == SC_OVER_TWO_NODE and _tuple_is(
            row.get(tuple_key), SC_OVER_TWO_EXACT_MPF_TUPLE
        ):
            return True
    nested = row.get("parameter_binding")
    if isinstance(nested, Mapping):
        half = nested.get("sc_over_two") or nested.get("sc_half")
        if isinstance(half, Mapping) and (
            half.get("node_id") == SC_OVER_TWO_NODE
            and _tuple_is(half.get("exact_mpf_tuple"), SC_OVER_TWO_EXACT_MPF_TUPLE)
        ):
            return True
    # A port may designate sc/2 by the two required parameter fields and carry
    # selected sc in the explicit fields checked above.
    return (
        row.get("parameter_node") == SC_OVER_TWO_NODE
        and _tuple_is(
            row.get("parameter_exact_mpf_tuple"), SC_OVER_TWO_EXACT_MPF_TUPLE
        )
        and (
            row.get("selected_sc_node") == SC_NODE
            or row.get("sc_node_id") == SC_NODE
        )
    )


def _coordinate_is_sc_over_two(row: Mapping[str, Any]) -> bool:
    """The canonical join port carries sc/2 in its coordinate tuple."""

    return row.get("coordinate_node") == COORDINATE_NODE and _tuple_is(
        row.get("coordinate_exact_mpf_tuple"), SC_OVER_TWO_EXACT_MPF_TUPLE
    )


def _port_row_is_bound(row: Any, order: int) -> bool:
    if not isinstance(row, Mapping):
        return False
    base = (
        row.get("archive_sha256") == Q4_SHA
        and row.get("expression_node") == PHYSICAL_U_Y_ROWS[order]
        and row.get("ordinary_Z_order") == order
        and row.get("coordinate_node") == COORDINATE_NODE
        and _tuple_is(
            row.get("coordinate_exact_mpf_tuple"), COORDINATE_EXACT_MPF_TUPLE
        )
    )
    return (
        base
        and _selected_sc_binding(row)
        and (_sc_over_two_binding(row) or _coordinate_is_sc_over_two(row))
    )


def _ports_are_three_rows(value: Any) -> bool:
    rows = _safe_get(value, "source_ports", "U_y")
    return isinstance(rows, list) and len(rows) == 3 and all(
        _port_row_is_bound(row, order) for order, row in enumerate(rows)
    )


def _diagnostic(case: str, value: Any) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        return {"case": case, "passed": False, "status": None}
    return {
        "case": case,
        "passed": value.get("passed") is True,
        "status": value.get("status"),
    }


def run_fixture(join: Callable[..., Any]) -> dict[str, Any]:
    """Run the positive join and independent fail-closed admission cases."""

    checks: dict[str, bool] = {}
    diagnostics: list[dict[str, Any]] = []
    if not callable(join):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"join_callable": False},
            "diagnostics": [{"case": "input", "passed": False}],
            "synthetic_only": True,
        }
    checks["join_callable"] = True

    instance_bundle, flat_bundle = _new_pair()
    try:
        positive = join(instance_bundle, flat_bundle)
        checks["positive_verified"] = (
            isinstance(positive, Mapping)
            and positive.get("status") == "verified"
            and positive.get("passed") is True
        )
        checks["positive_has_checks_mapping"] = isinstance(
            _safe_get(positive, "checks"), Mapping
        )
        checks["positive_scope_is_new_instance_only"] = _positive_scope_is_narrow(
            positive
        )
        checks["positive_source_ports_U_y_three_rows"] = _ports_are_three_rows(
            positive
        )
        diagnostics.append(_diagnostic("positive", positive))
    except Exception as exc:
        checks["positive_verified"] = False
        diagnostics.append(
            {"case": "positive", "passed": False, "exception": type(exc).__name__}
        )

    def reject(case: str, mutate: Callable[[dict[str, Any], dict[str, Any]], None]) -> None:
        case_instance, case_flat = _new_pair()
        mutate(case_instance, case_flat)
        try:
            result = join(case_instance, case_flat)
            checks["reject_" + case] = _result_rejected(result)
            diagnostics.append(_diagnostic(case, result))
        except Exception as exc:
            checks["reject_" + case] = False
            diagnostics.append(
                {"case": case, "passed": False, "exception": type(exc).__name__}
            )

    reject(
        "wrong_owner",
        lambda instance, _flat: instance["result"]["instance_result"]["binding"][
            "graph_id"
        ].update({"q4": "f" * 64, "q4_input_hash": "f" * 64}),
    )
    reject(
        "wrong_hash",
        lambda _instance, flat: flat["result"]["input_hashes"].update(
            {Q4_ARCHIVE: "e" * 64}
        ),
    )
    reject(
        "non_positive_sc",
        lambda instance, _flat: instance["result"]["instance_result"]["binding"].update(
            {"sc_exact_mpf_tuple": [1, 3, -6, 2]}
        ),
    )
    reject(
        "noncanonical_tuple",
        lambda instance, _flat: instance["result"]["instance_result"]["binding"].update(
            {"sc_exact_mpf_tuple": [0, 6, -7, 3]}
        ),
    )
    reject(
        "old_archive_assignment_claimed",
        lambda instance, _flat: instance["result"]["instance_result"]["scope"].update(
            {"existing_archive_assignment_proved": True}
        ),
    )
    reject(
        "flat_scope_not_established",
        lambda _instance, flat: flat["result"]["archive_result"]["scope"].update(
            {"Q4_Uy_correction_zero_at_endpoint_conditional_on_sc_positive": False}
        ),
    )
    reject(
        "modified_ordinary_Z_order",
        lambda _instance, flat: (
            flat["result"]["archive_result"]["scope"].update(
                {"ordinary_Z_orders": [0, 1, 3]}
            ),
            flat["result"]["archive_result"]["evidence"].update(
                {"ordinary_Z_orders": [0, 1, 3]}
            ),
        ),
    )
    reject(
        "old_errors_removed",
        lambda instance, _flat: instance["result"]["instance_result"]["scope"].update(
            {"old_errors_removed": True}
        ),
    )
    reject(
        "missing_source_rows",
        lambda _instance, flat: flat["result"]["archive_result"]["evidence"].update(
            {"physical_U_y_rows": [PHYSICAL_U_Y_ROWS[0], PHYSICAL_U_Y_ROWS[1]]}
        ),
    )

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(value is True for value in checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "test_count": len(checks),
        "synthetic_only": True,
        "source_code_executed": False,
        "source_graph_walked": False,
        "source_values_evaluated": False,
        "cross_owner_identity_proved": False,
        "old_archive_assignment_proved": False,
        "global_closure_proved": False,
    }


def _valid_port() -> dict[str, Any]:
    """Build one fresh port for the optional port-evaluator fixture."""

    return {
        "schema": "q4_physical_Uy_endpoint_source_port_v1",
        "archive_sha256": Q4_SHA,
        "expression_node": PHYSICAL_U_Y_ROWS[0],
        "ordinary_Z_order": 0,
        "parameter_node": SC_NODE,
        "parameter_exact_mpf_tuple": list(SC_EXACT_MPF_TUPLE),
        "selected_sc_node": SC_NODE,
        "selected_sc_exact_mpf_tuple": list(SC_EXACT_MPF_TUPLE),
        "sc_over_two_node": SC_OVER_TWO_NODE,
        "sc_over_two_exact_mpf_tuple": list(SC_OVER_TWO_EXACT_MPF_TUPLE),
        "coordinate_node": COORDINATE_NODE,
        "coordinate_exact_mpf_tuple": list(COORDINATE_EXACT_MPF_TUPLE),
    }


def run_port_fixture(evaluate_port: Callable[..., Any]) -> dict[str, Any]:
    """Check one-call port forwarding and fail-closed port admission."""

    checks: dict[str, bool] = {}
    diagnostics: list[dict[str, Any]] = []
    if not callable(evaluate_port):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"evaluate_port_callable": False},
            "diagnostics": [{"case": "input", "passed": False}],
            "synthetic_only": True,
        }
    checks["evaluate_port_callable"] = True

    port = _valid_port()
    calls: list[Any] = []
    sentinel = {"resolved": "synthetic-port-result"}

    def resolver(value: Any) -> Any:
        calls.append(value)
        return sentinel

    try:
        result = evaluate_port(port, resolver)
        checks["positive_resolver_called_once"] = (
            len(calls) == 1 and isinstance(calls[0], Mapping) and calls[0] == port
        )
        checks["positive_returns_resolver_result"] = result is sentinel
        diagnostics.append(
            {"case": "positive", "passed": result is sentinel, "calls": len(calls)}
        )
    except Exception as exc:
        checks["positive_resolver_called_once"] = False
        checks["positive_returns_resolver_result"] = False
        diagnostics.append(
            {"case": "positive", "passed": False, "exception": type(exc).__name__}
        )

    def reject_port(case: str, mutate: Callable[[dict[str, Any]], None]) -> None:
        case_port = _valid_port()
        mutate(case_port)
        case_calls: list[Any] = []

        def case_resolver(value: Any) -> Any:
            case_calls.append(value)
            return sentinel

        try:
            rejected = evaluate_port(case_port, case_resolver)
            checks["reject_" + case + "_without_resolver"] = (
                len(case_calls) == 0 and rejected is not sentinel
            )
            diagnostics.append(
                {
                    "case": case,
                    "passed": checks["reject_" + case + "_without_resolver"],
                    "calls": len(case_calls),
                }
            )
        except Exception as exc:
            checks["reject_" + case + "_without_resolver"] = len(case_calls) == 0
            diagnostics.append(
                {
                    "case": case,
                    "passed": checks["reject_" + case + "_without_resolver"],
                    "calls": len(case_calls),
                    "exception": type(exc).__name__,
                }
            )

    reject_port("missing_expression_node", lambda value: value.pop("expression_node"))
    reject_port(
        "damaged_parameter_tuple",
        lambda value: value.update({"parameter_exact_mpf_tuple": [0, 6, -7, 3]}),
    )

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(value is True for value in checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "test_count": len(checks),
        "synthetic_only": True,
        "source_code_executed": False,
        "source_graph_walked": False,
        "source_values_evaluated": False,
    }


__all__ = ["SCHEMA", "run_fixture", "run_port_fixture"]
