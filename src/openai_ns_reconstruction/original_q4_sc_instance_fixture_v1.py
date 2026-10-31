"""Fresh synthetic positive and negative cases for the Q4 ``sc`` adapter.

The fixture has no dependency on an archive, graph producer, or historical
receipt.  Each case constructs new dictionaries before its own mutation, so
one case cannot contaminate another.
"""
from __future__ import annotations

from typing import Any, Callable


SCHEMA = "original_q4_sc_instance_fixture_v1"


def _positive() -> tuple[dict[str, Any], dict[str, Any], int, str]:
    """Create one independent synthetic receipt and endpoint."""

    q4_hash = "1" * 64
    loop_hash = "2" * 64
    sc_node = 17
    sc_hash = "a" * 64
    receipt = {
        "result": {
            "archive_result": {
                "status": "verified",
                "passed": True,
                "checks": {"archive_family_metadata_matches": True},
                "evidence": {
                    "sc": sc_node,
                    "selected_shared_node_sha256": {"sc": sc_hash},
                },
                "scope": {
                    "archived_selected_source_nodes_only": True,
                    "full_graph_prefix_equality_claimed": False,
                    "old_errors_removed": False,
                },
            },
            "input_hashes": {
                "synthetic_q4_archive.json.gz": q4_hash,
                "synthetic_loop_archive.json.gz": loop_hash,
            },
        }
    }
    # 3 * 2**(-6) = 3/64; bitcount + exponent is -4 <= -3.
    endpoint = {
        "lower_exact_mpf_tuple": [0, 3, -6, 2],
        "upper_exact_mpf_tuple": [0, 3, -6, 2],
    }
    return receipt, endpoint, sc_node, sc_hash


def _safe_get(value: Any, *keys: str) -> Any:
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            return None
        value = value[key]
    return value


def _scope_flags(value: Any) -> bool:
    scope = _safe_get(value, "scope")
    return (
        isinstance(scope, dict)
        and scope.get("newly_chosen_binding") is True
        and scope.get("existing_archive_assignment_proved") is False
        and scope.get("corrected_owner_identity_proved") is False
        and scope.get("old_errors_removed") is False
    )


def run_fixture(bind: Callable[..., Any]) -> dict[str, Any]:
    """Run the fresh synthetic admission matrix against ``bind``."""

    checks: dict[str, bool] = {}
    diagnostics: dict[str, Any] = {}
    if not callable(bind):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"bind_callable": False},
            "diagnostics": {"input": "bind is not callable"},
            "synthetic_only": True,
        }
    checks["bind_callable"] = True

    receipt, endpoint, sc_node, sc_hash = _positive()
    try:
        positive = bind(receipt, endpoint, sc_node, sc_hash)
        checks["positive_passed"] = isinstance(positive, dict) and positive.get("passed") is True
        checks["positive_verified"] = (
            isinstance(positive, dict) and positive.get("status") == "verified"
        )
        checks["positive_graph_identifier"] = (
            _safe_get(positive, "binding", "graph_id") == {
                "q4": "1" * 64,
                "loop": "2" * 64,
                "q4_input_hash": "1" * 64,
                "loop_input_hash": "2" * 64,
            }
        )
        checks["positive_sc_node_identity"] = (
            _safe_get(positive, "binding", "sc_node_id") == sc_node
        )
        checks["positive_exact_sc_tuple"] = (
            _safe_get(positive, "binding", "sc_exact_mpf_tuple") == [0, 3, -6, 2]
        )
        checks["positive_source_family_is_scoped"] = (
            _safe_get(positive, "binding", "source_family",
                      "receipt_scope_present") is True
            and _safe_get(positive, "binding", "source_family",
                          "source_family_value_exposed") is False
        )
        checks["positive_sc_half_explicit"] = (
            _safe_get(positive, "derived_dyadics", "sc_over_two", "exact_mpf_tuple")
            == [0, 3, -7, 2]
        )
        checks["positive_three_sc_quarters_explicit"] = (
            _safe_get(positive, "derived_dyadics", "three_sc_over_four", "exact_mpf_tuple")
            == [0, 9, -8, 4]
        )
        checks["positive_comparison_proved"] = (
            _safe_get(positive, "derived_dyadics", "comparison",
                      "sc_over_two_lt_three_sc_over_four") is True
            and _safe_get(positive, "derived_dyadics", "comparison",
                          "three_sc_over_four_lt_sc") is True
        )
        checks["positive_no_power_expansion"] = (
            _safe_get(positive, "derived_dyadics", "comparison",
                      "power_of_two_materialized") is False
        )
        checks["positive_scope_flags"] = _scope_flags(positive)
        checks["positive_conditional_flat_requires_external_receipt"] = (
            _safe_get(positive, "q4_flat_conclusion", "used") is False
            and _safe_get(positive, "q4_flat_conclusion",
                          "requires_external_authenticated_receipt") is True
            and _safe_get(positive, "q4_flat_conclusion",
                          "numeric_field_evaluation") is False
        )
        checks["positive_no_midpoint_claim"] = (
            _safe_get(positive, "binding", "endpoint",
                      "analytic_value_equal_interval_midpoint_claimed") is False
        )
        checks["positive_no_source_evaluation"] = (
            positive.get("source_values_evaluated") is False
            and positive.get("source_graph_walked") is False
        )
        diagnostics["positive"] = {
            "passed": checks["positive_passed"],
            "status": positive.get("status") if isinstance(positive, dict) else None,
        }
    except Exception as exc:
        checks["positive_no_exception"] = False
        diagnostics["positive_exception"] = type(exc).__name__ + ": " + str(exc)

    def reject(case: str, mutate: Callable[[dict[str, Any], dict[str, Any]], None],
               *, expected_hash: str | None = None) -> None:
        # _positive() is called for every case: no case mutates a shared base.
        case_receipt, case_endpoint, case_id, case_hash = _positive()
        mutate(case_receipt, case_endpoint)
        supplied_hash = case_hash if expected_hash is None else expected_hash
        try:
            result = bind(case_receipt, case_endpoint, case_id, supplied_hash)
            checks["reject_" + case] = (
                isinstance(result, dict) and result.get("passed") is False
            )
            checks["reject_" + case + "_no_binding"] = (
                isinstance(result, dict) and result.get("binding") is None
            )
            diagnostics[case] = {
                "passed": checks["reject_" + case],
                "status": result.get("status") if isinstance(result, dict) else None,
            }
        except Exception as exc:
            checks["reject_" + case] = False
            checks["reject_" + case + "_no_binding"] = False
            diagnostics[case + "_exception"] = type(exc).__name__ + ": " + str(exc)

    reject(
        "zero",
        lambda _receipt, endpoint: endpoint.update({
            "lower_exact_mpf_tuple": [0, 0, 0, 0],
            "upper_exact_mpf_tuple": [0, 0, 0, 0],
        }),
    )
    reject(
        "negative",
        lambda _receipt, endpoint: endpoint.update({
            "lower_exact_mpf_tuple": [1, 3, -6, 2],
            "upper_exact_mpf_tuple": [1, 3, -6, 2],
        }),
    )
    reject(
        "non_singleton",
        lambda _receipt, endpoint: endpoint.update({
            "lower_exact_mpf_tuple": [0, 3, -6, 2],
            "upper_exact_mpf_tuple": [0, 1, -4, 1],
        }),
    )
    reject(
        "wrong_bitcount",
        lambda _receipt, endpoint: endpoint.update({
            "lower_exact_mpf_tuple": [0, 3, -6, 3],
            "upper_exact_mpf_tuple": [0, 3, -6, 3],
        }),
    )
    reject(
        "bool_input",
        lambda _receipt, endpoint: endpoint.update({
            "lower_exact_mpf_tuple": [False, 3, -6, 2],
            "upper_exact_mpf_tuple": [False, 3, -6, 2],
        }),
    )
    reject(
        "non_int_exponent",
        lambda _receipt, endpoint: endpoint.update({
            "lower_exact_mpf_tuple": [0, 3, -6.0, 2],
            "upper_exact_mpf_tuple": [0, 3, -6.0, 2],
        }),
    )
    reject(
        "wrong_receipt_status",
        lambda receipt, _endpoint: receipt["result"]["archive_result"].update({
            "status": "mismatch",
        }),
    )
    reject(
        "wrong_receipt_passed_flag",
        lambda receipt, _endpoint: receipt["result"]["archive_result"].update({
            "passed": False,
        }),
    )
    reject(
        "sc_node_tampered",
        lambda receipt, _endpoint: receipt["result"]["archive_result"]["evidence"].update({
            "sc": 18,
        }),
    )
    reject(
        "sc_hash_tampered",
        lambda receipt, _endpoint: receipt["result"]["archive_result"]["evidence"].update({
            "selected_shared_node_sha256": {"sc": "b" * 64},
        }),
    )
    reject(
        "wrong_expected_hash",
        lambda _receipt, _endpoint: None,
        expected_hash="c" * 64,
    )
    reject(
        "missing_scope",
        lambda receipt, _endpoint: receipt["result"]["archive_result"].pop(
            "scope"
        ),
    )
    reject(
        "missing_q4_hash",
        lambda receipt, _endpoint: receipt["result"]["input_hashes"].pop(
            "synthetic_q4_archive.json.gz"
        ),
    )
    reject(
        "not_below_one_eighth",
        lambda _receipt, endpoint: endpoint.update({
            "lower_exact_mpf_tuple": [0, 1, -2, 1],
            "upper_exact_mpf_tuple": [0, 1, -2, 1],
        }),
    )
    reject(
        "noncanonical_even_mantissa",
        lambda _receipt, endpoint: endpoint.update({
            "lower_exact_mpf_tuple": [0, 6, -7, 3],
            "upper_exact_mpf_tuple": [0, 6, -7, 3],
        }),
    )

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "synthetic_only": True,
        "source_code_executed": False,
        "source_graph_walked": False,
        "source_values_evaluated": False,
        "old_archive_data_reused": False,
        "fresh_dictionary_per_case": True,
    }


__all__ = ["SCHEMA", "run_fixture"]
