"""Validate a saved failed occurrence-leaf result without rerunning it.

This module consumes only a caller-supplied saved result mapping.  It does not
read a graph, call the earlier occurrence leaf, inspect expression paths, or
revalidate scope.  The seven leaf records are preserved as prior structural
evidence.  Structural occurrence is explicitly not binder identity.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any


SCHEMA = "n1_integral_occurrence_saved_audit_v1"
RESULT_SCHEMA = SCHEMA + "_result"
RUNNER_SCHEMA = "n1_integral_scoped_occurrence_v1_runner_result_v1"
COMPONENT = "n1_integral_scoped_occurrence_v1"
FIXTURE_SCHEMA = "n1_integral_scoped_occurrence_fixture_v1"
LEAF_SCHEMA = "n1_integral_scoped_occurrence_v1_result"
EXPECTED_COUNT = 7
LEAF_STATUS = "unresolved"
RECORD_STATUSES = frozenset(
    {
        "unresolved",
        "found_structural_only",
        "found_but_unresolved",
        "not_found_complete",
    }
)

_FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "candidate_actual_integral_binder",
        "candidate_actual_integral_binder_proved",
        "candidate_actual_integral_binder_verified",
        "binder_identity",
        "binder_identity_proved",
        "binder_identity_verified",
        "complete_identity_proved",
        "complete_identity_verified",
        "same_binder_proved",
        "same_binder_verified",
        "numeric_proof",
        "certified",
        "integral_value_computed",
        "actual_integral_value_computed",
        "source_value_computed",
        "source_values_computed",
        "numeric_values_computed",
        "numeric_value_computed",
        "numeric_values",
        "numeric_evaluated",
        "numeric_values_evaluated",
        "actual_point_values_computed",
        "numeric_equivalence_certified",
        "integral_value_evaluated",
        "integral_values_evaluated",
        "integral_evaluated",
        "integral_values",
        "m15_claimed",
        "kchi_claimed",
    }
)
_FORBIDDEN_FALSE_ONLY_KEYS = _FORBIDDEN_TRUE_KEYS | frozenset(
    {
        "numeric_field_value_computed",
        "source_value_evaluated",
        "integrals_evaluated",
        "binder_mappings_installed",
        "actual_graph_evaluate_called",
        "graph_evaluated",
        "source_values_evaluated",
    }
)


def _is_int(value: Any) -> bool:
    return type(value) is int


def _is_nonnegative_int(value: Any) -> bool:
    return _is_int(value) and value >= 0


def _fixture_errors(value: Any, path: str) -> list[str]:
    if not isinstance(value, Mapping):
        return [path + " is not an object"]
    errors: list[str] = []
    if value.get("schema") != FIXTURE_SCHEMA:
        errors.append(path + ".schema mismatch")
    if value.get("passed") is not True:
        errors.append(path + ".passed is not true")
    return errors


def _claim_errors(value: Any, path: str = "saved_result") -> list[str]:
    """Reject true proof/value claims without interpreting structural paths."""

    errors: list[str] = []
    if isinstance(value, Mapping):
        for raw_key, item in value.items():
            key = str(raw_key).strip().lower().replace("-", "_")
            item_path = path + "." + str(raw_key)
            if key in _FORBIDDEN_FALSE_ONLY_KEYS:
                allowed_false = item is False or item is None or (
                    type(item) is int and item == 0
                )
                if not allowed_false:
                    errors.append("forbidden claim at " + item_path)
            errors.extend(_claim_errors(item, item_path))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            errors.extend(_claim_errors(item, path + "[" + str(index) + "]"))
    return errors


def _failure_error(saved: Mapping[str, Any]) -> str | None:
    value = saved.get("error")
    if (
        isinstance(value, str)
        and "NameError" in value
        and "name '_is_int' is not defined" in value
    ):
        return None
    return "saved error is not the expected _is_int NameError"


def _expected_record(
    row: Any, index: int
) -> tuple[int, int, int, str] | str:
    if not isinstance(row, Mapping):
        return "progress.integral_records[" + str(index) + "] is not an object"
    integral = row.get("integral_node_id")
    candidate = row.get("unique_node_id")
    original = row.get("original_record")
    integrand = (
        original.get("integrand")
        if isinstance(original, Mapping)
        else None
    )
    if not _is_nonnegative_int(integral):
        return "saved integral record has invalid nonnegative integral_node_id at " + str(index)
    if not _is_nonnegative_int(integrand):
        return "saved integral record has invalid nonnegative original_record.integrand at " + str(index)
    if not _is_nonnegative_int(candidate):
        return "saved integral record has invalid nonnegative unique_node_id at " + str(index)
    if row.get("binding_status") != "unique_candidate":
        return "saved integral record binding_status is not unique_candidate at " + str(index)
    candidate_ids = row.get("candidate_node_ids")
    if (
        type(candidate_ids) is not list
        or len(candidate_ids) != 1
        or not _is_nonnegative_int(candidate_ids[0])
        or candidate_ids[0] != candidate
    ):
        return "saved integral record candidate_node_ids mismatch at " + str(index)
    if not isinstance(original, Mapping):
        return "saved integral record original_record is not an object at " + str(index)
    if original.get("operation") != "definite_integral":
        return "saved integral record operation is not definite_integral at " + str(index)
    variable = row.get("variable_name")
    original_variable = original.get("variable")
    if not isinstance(variable, str) or not variable:
        return "saved integral record variable_name is invalid at " + str(index)
    if original_variable != variable:
        return "saved integral record original variable mismatch at " + str(index)
    return integral, integrand, candidate, variable


def _reject(
    errors: list[str],
    *,
    saved_records: Any = None,
    catalog_records: Any = None,
) -> dict[str, Any]:
    return {
        "schema": RESULT_SCHEMA,
        "status": "rejected",
        "errors": list(errors),
        "saved_leaf_records": saved_records if isinstance(saved_records, list) else [],
        "saved_catalog_records": (
            catalog_records if isinstance(catalog_records, list) else []
        ),
        "source_traversal_output_preserved": False,
        "saved_metadata_only": True,
        "graph_or_identity_authenticated": False,
        "source_traversal_reexecuted": False,
        "graph_read": False,
        "graph_traversed": False,
        "paths_independently_revalidated": False,
        "scopes_independently_revalidated": False,
        "candidate_actual_integral_binder_proved": False,
        "binder_identity_proved": False,
        "complete_identity_proved": False,
        "numeric_field_value_computed": False,
        "integrals_evaluated": False,
        "binder_mappings_installed": False,
    }


def audit_saved_result(saved_result: Any) -> dict[str, Any]:
    """Validate one saved failed run and preserve its raw structural records."""

    if not isinstance(saved_result, Mapping):
        return _reject(["saved result must be a mapping"])

    errors: list[str] = []
    if saved_result.get("schema") != RUNNER_SCHEMA:
        errors.append("saved schema mismatch")
    if saved_result.get("component") != COMPONENT:
        errors.append("saved component/family mismatch")
    if saved_result.get("status") != "failed":
        errors.append("saved status is not failed")
    if saved_result.get("passed") is not False:
        errors.append("saved passed is not false")
    failure_error = _failure_error(saved_result)
    if failure_error is not None:
        errors.append(failure_error)
    errors.extend(_fixture_errors(saved_result.get("fixture"), "fixture"))

    progress = saved_result.get("progress")
    if not isinstance(progress, Mapping):
        return _reject(errors + ["saved progress is missing"])
    errors.extend(_fixture_errors(progress.get("fixture"), "progress.fixture"))
    if progress.get("leaf_calls") != 1:
        errors.append("progress.leaf_calls is not one")

    leaf_raw = progress.get("leaf_raw_output")
    if not isinstance(leaf_raw, Mapping):
        return _reject(errors + ["progress.leaf_raw_output is missing"])
    if leaf_raw.get("schema") != LEAF_SCHEMA:
        errors.append("saved leaf_raw_output schema mismatch")
    if leaf_raw.get("status") != LEAF_STATUS:
        errors.append("saved leaf_raw_output status is not unresolved")
    if leaf_raw.get("catalog_record_count") != EXPECTED_COUNT:
        errors.append("saved leaf catalog count is not seven")

    catalog_records = progress.get("integral_records")
    leaf_records = leaf_raw.get("records")
    if not isinstance(catalog_records, list) or len(catalog_records) != EXPECTED_COUNT:
        errors.append("progress.integral_records does not contain seven records")
        catalog_records = catalog_records if isinstance(catalog_records, list) else []
    if not isinstance(leaf_records, list) or len(leaf_records) != EXPECTED_COUNT:
        errors.append("saved leaf records do not contain seven records")
        leaf_records = leaf_records if isinstance(leaf_records, list) else []

    if len(catalog_records) == EXPECTED_COUNT and len(leaf_records) == EXPECTED_COUNT:
        seen_integrals: set[int] = set()
        for index, (catalog_row, leaf_row) in enumerate(
            zip(catalog_records, leaf_records)
        ):
            expected = _expected_record(catalog_row, index)
            if isinstance(expected, str):
                errors.append(expected)
                continue
            integral, integrand, candidate, variable = expected
            if integral in seen_integrals:
                errors.append("duplicate integral_node_id at " + str(index))
            seen_integrals.add(integral)
            if not isinstance(leaf_row, Mapping):
                errors.append("leaf record " + str(index) + " is not an object")
                continue
            if leaf_row.get("status") not in RECORD_STATUSES:
                errors.append("leaf record " + str(index) + " has unsupported status")
            if leaf_row.get("integral_node_id") != integral:
                errors.append("leaf record " + str(index) + " integral_node_id mismatch")
            if leaf_row.get("integrand_node_id") != integrand:
                errors.append("leaf record " + str(index) + " integrand_node_id mismatch")
            if leaf_row.get("candidate_node_id") != candidate:
                errors.append("leaf record " + str(index) + " candidate_node_id mismatch")
            if not _is_nonnegative_int(leaf_row.get("integral_node_id")):
                errors.append("leaf record " + str(index) + " integral_node_id is not nonnegative int")
            if not _is_nonnegative_int(leaf_row.get("integrand_node_id")):
                errors.append("leaf record " + str(index) + " integrand_node_id is not nonnegative int")
            if not _is_nonnegative_int(leaf_row.get("candidate_node_id")):
                errors.append("leaf record " + str(index) + " candidate_node_id is not nonnegative int")
            if leaf_row.get("variable_name") != variable:
                errors.append("leaf record " + str(index) + " variable_name mismatch")
            candidate_ids = leaf_row.get("candidate_node_ids")
            if candidate_ids is not None:
                if (
                    type(candidate_ids) is not list
                    or len(candidate_ids) != 1
                    or not _is_nonnegative_int(candidate_ids[0])
                    or candidate_ids[0] != candidate
                ):
                    errors.append("leaf record " + str(index) + " candidate_node_ids mismatch")

    claim_errors = _claim_errors(saved_result)
    errors.extend(claim_errors)
    if errors:
        return _reject(
            errors,
            saved_records=leaf_records,
            catalog_records=catalog_records,
        )

    # The records are deliberately retained without inspecting occurrence
    # paths, frontier reasons, or lexical tokens.  This is saved evidence from
    # the earlier traversal, not an independent traversal or proof.
    return {
        "schema": RESULT_SCHEMA,
        "status": "validated_saved_failure",
        "saved_schema": RUNNER_SCHEMA,
        "component": COMPONENT,
        "saved_status": "failed",
        "saved_error": saved_result.get("error"),
        "fixture_passed": True,
        "fixture_schema": FIXTURE_SCHEMA,
        "leaf_schema": LEAF_SCHEMA,
        "saved_metadata_only": True,
        "graph_or_identity_authenticated": False,
        "leaf_calls": 1,
        "catalog_record_count": EXPECTED_COUNT,
        "saved_leaf_records": leaf_records,
        "saved_catalog_records": catalog_records,
        "source_traversal_output_preserved": True,
        "source_traversal_reexecuted": False,
        "graph_read": False,
        "graph_traversed": False,
        "paths_independently_revalidated": False,
        "scopes_independently_revalidated": False,
        "structural_occurrence_is_not_identity": True,
        "candidate_actual_integral_binder_proved": False,
        "binder_identity_proved": False,
        "complete_identity_proved": False,
        "numeric_field_value_computed": False,
        "source_value_evaluated": False,
        "integrals_evaluated": False,
        "binder_mappings_installed": False,
        "claims": {
            "saved_records_only": True,
            "structural_occurrence_is_not_identity": True,
            "independent_paths_revalidated": False,
            "complete_proof": False,
            "candidate_actual_integral_binder": False,
            "numeric_values": False,
            "integral_values": False,
        },
    }


__all__ = [
    "SCHEMA",
    "RESULT_SCHEMA",
    "RUNNER_SCHEMA",
    "COMPONENT",
    "FIXTURE_SCHEMA",
    "LEAF_SCHEMA",
    "EXPECTED_COUNT",
    "LEAF_STATUS",
    "RECORD_STATUSES",
    "audit_saved_result",
]
