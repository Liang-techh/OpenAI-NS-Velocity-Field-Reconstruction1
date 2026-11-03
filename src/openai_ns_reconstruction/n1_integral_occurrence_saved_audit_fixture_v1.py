"""Fresh synthetic checks for the saved occurrence-result validator.

The fixture creates no graph and does not import or call any earlier leaf.
Its seven scalar records use new IDs and short synthetic paths only to test
saved-record alignment, rejection boundaries, and the non-identity claim.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


SCHEMA = "n1_integral_occurrence_saved_audit_fixture_v1"


def _saved_failure() -> dict[str, Any]:
    catalog: list[dict[str, Any]] = []
    leaf_records: list[dict[str, Any]] = []
    for index in range(7):
        integral = 7000 + index
        integrand = 7100 + index
        candidate = 7200 + (index % 2)
        variable = "synthetic_scalar_v"
        catalog.append(
            {
                "integral_node_id": integral,
                "unique_node_id": candidate,
                "variable_name": variable,
                "binding_status": "unique_candidate",
                "candidate_node_ids": [candidate],
                "original_record": {
                    "integrand": integrand,
                    "operation": "definite_integral",
                    "variable": variable,
                },
            }
        )
        leaf_records.append(
            {
                "integral_node_id": integral,
                "integrand_node_id": integrand,
                "candidate_node_id": candidate,
                "candidate_node_ids": [candidate],
                "variable_name": variable,
                "status": "found_but_unresolved",
                "occurrences": [
                    {
                        "node_id": candidate,
                        "path": [
                            {"node_id": integrand, "via": "integrand_root"},
                            {"node_id": candidate, "via": "synthetic.scalar"},
                        ],
                        "target_id_unshadowed": True,
                    }
                ],
                "opaque_frontier": [
                    {
                        "node_id": 7300 + index,
                        "reason": "synthetic_frontier",
                        "path": [
                            {"node_id": integrand, "via": "integrand_root"},
                            {"node_id": 7300 + index, "via": "opaque"},
                        ],
                    }
                ],
                "candidate_actual_integral_binder": False,
                "complete_identity_proved": False,
                "numeric_field_value_computed": False,
                "integrals_evaluated": False,
            }
        )
    return {
        "schema": "n1_integral_scoped_occurrence_v1_runner_result_v1",
        "component": "n1_integral_scoped_occurrence_v1",
        "status": "failed",
        "passed": False,
        "error": "NameError: name '_is_int' is not defined",
        "fixture": {
            "schema": "n1_integral_scoped_occurrence_fixture_v1",
            "passed": True,
            "synthetic_only": True,
        },
        "progress": {
            "fixture": {
                "schema": "n1_integral_scoped_occurrence_fixture_v1",
                "passed": True,
                "synthetic_only": True,
            },
            "leaf_calls": 1,
            "integral_records": catalog,
            "leaf_raw_output": {
                "schema": "n1_integral_scoped_occurrence_v1_result",
                "status": "unresolved",
                "catalog_record_count": 7,
                "records": leaf_records,
                "candidate_actual_integral_binder_proved": False,
                "complete_identity_proved": False,
                "numeric_field_value_computed": False,
                "integrals_evaluated": False,
            },
        },
        "candidate_actual_integral_binder_proved": False,
        "complete_identity_proved": False,
        "numeric_field_value_computed": False,
        "integrals_evaluated": False,
    }


def _copy_saved() -> dict[str, Any]:
    """Build a fresh value instead of sharing or mutating a prior case."""
    return _saved_failure()


def _audit(module: Any, saved: Any) -> Any:
    function = getattr(module, "audit_saved_result", None)
    if not callable(function):
        return None
    return function(saved)


def _status(result: Any) -> Any:
    return result.get("status") if isinstance(result, Mapping) else None


def _diagnostic(
    name: str,
    result: Any,
    checks: Mapping[str, bool],
    exception: str | None,
) -> dict[str, Any]:
    item: dict[str, Any] = {
        "case": name,
        "status": _status(result),
        "passed": all(value is True for value in checks.values()),
        "failed_checks": [
            key for key, value in checks.items() if value is not True
        ],
    }
    if exception is not None:
        item["exception"] = exception
    return item


def run_fixture(module: Any) -> dict[str, Any]:
    """Run saved-result boundary checks against the supplied new module."""

    function = getattr(module, "audit_saved_result", None)
    if not callable(function):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"audit_callable": False},
            "diagnostics": [{"case": "input", "passed": False}],
            "case_count": 0,
            "test_count": 1,
            "synthetic_only": True,
            "graph_read": False,
            "saved_metadata_only": True,
            "graph_or_identity_authenticated": False,
            "numeric_values_computed": False,
        }

    checks: dict[str, bool] = {"audit_callable": True}
    diagnostics: list[dict[str, Any]] = []

    cases = (
        "valid_saved_failure",
        "runner_schema_mismatch",
        "component_family_mismatch",
        "fixture_schema_mismatch",
        "progress_fixture_schema_mismatch",
        "leaf_schema_mismatch",
        "variable_name_mismatch",
        "catalog_variable_name_mismatch",
        "original_variable_mismatch",
        "integral_id_mismatch",
        "integrand_value_mismatch",
        "candidate_ids_mismatch",
        "forged_identity_claim",
        "forged_numeric_claim",
        "wrong_failure_error",
        "unsupported_record_status",
    )
    for name in cases:
        saved = _copy_saved()
        if name == "runner_schema_mismatch":
            saved["schema"] = "n1_other_runner_result_v1"
        elif name == "component_family_mismatch":
            saved["component"] = "n1_other_component_v1"
        elif name == "fixture_schema_mismatch":
            saved["fixture"]["schema"] = "n1_other_fixture_v1"
        elif name == "progress_fixture_schema_mismatch":
            saved["progress"]["fixture"]["schema"] = "n1_other_fixture_v1"
        elif name == "leaf_schema_mismatch":
            saved["progress"]["leaf_raw_output"][
                "schema"
            ] = "n1_other_leaf_result_v1"
        elif name == "variable_name_mismatch":
            saved["progress"]["leaf_raw_output"]["records"][0][
                "variable_name"
            ] = "wrong_variable"
        elif name == "catalog_variable_name_mismatch":
            saved["progress"]["integral_records"][0][
                "variable_name"
            ] = "wrong_variable"
        elif name == "original_variable_mismatch":
            saved["progress"]["integral_records"][0]["original_record"][
                "variable"
            ] = "wrong_variable"
        elif name == "integral_id_mismatch":
            saved["progress"]["leaf_raw_output"]["records"][3][
                "integral_node_id"
            ] += 1
        elif name == "integrand_value_mismatch":
            saved["progress"]["integral_records"][4]["original_record"][
                "integrand"
            ] += 1
        elif name == "candidate_ids_mismatch":
            saved["progress"]["integral_records"][2][
                "candidate_node_ids"
            ] = [99999]
        elif name == "forged_identity_claim":
            saved["progress"]["leaf_raw_output"]["records"][0][
                "candidate_actual_integral_binder"
            ] = True
        elif name == "forged_numeric_claim":
            saved["progress"]["leaf_raw_output"]["records"][0][
                "numeric_field_value_computed"
            ] = True
        elif name == "wrong_failure_error":
            saved["error"] = "ValueError: unrelated failure"
        elif name == "unsupported_record_status":
            saved["progress"]["leaf_raw_output"]["records"][1][
                "status"
            ] = "complete"

        result: Any = None
        exception: str | None = None
        try:
            result = _audit(module, saved)
        except Exception as exc:
            exception = type(exc).__name__

        prefix = name + "_"
        if name == "valid_saved_failure":
            records = result.get("saved_leaf_records") if isinstance(result, Mapping) else None
            checks[prefix + "accepted"] = (
                _status(result) == "validated_saved_failure"
            )
            checks[prefix + "seven_preserved"] = (
                isinstance(records, list)
                and len(records) == 7
                and records == saved["progress"]["leaf_raw_output"]["records"]
            )
            checks[prefix + "not_identity"] = (
                isinstance(result, Mapping)
                and result.get("structural_occurrence_is_not_identity") is True
                and result.get("candidate_actual_integral_binder_proved") is False
                and result.get("complete_identity_proved") is False
            )
            checks[prefix + "not_revalidated"] = (
                isinstance(result, Mapping)
                and result.get("paths_independently_revalidated") is False
                and result.get("scopes_independently_revalidated") is False
                and result.get("graph_traversed") is False
            )
        else:
            checks[prefix + "rejected"] = _status(result) == "rejected"
            checks[prefix + "no_identity"] = (
                isinstance(result, Mapping)
                and result.get("candidate_actual_integral_binder_proved") is False
                and result.get("complete_identity_proved") is False
            )
        diagnostics.append(
            _diagnostic(
                name,
                result,
                {
                    key: value
                    for key, value in checks.items()
                    if key.startswith(prefix)
                },
                exception,
            )
        )

    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(value is True for value in checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "case_count": len(cases),
        "test_count": len(checks),
        "synthetic_only": True,
        "graph_read": False,
        "graph_traversed": False,
        "source_traversal_reexecuted": False,
        "saved_metadata_only": True,
        "graph_or_identity_authenticated": False,
        "numeric_values_computed": False,
        "integrals_evaluated": False,
        "binder_identity_proved": False,
        "complete_identity_proved": False,
    }


__all__ = ["SCHEMA", "run_fixture"]
