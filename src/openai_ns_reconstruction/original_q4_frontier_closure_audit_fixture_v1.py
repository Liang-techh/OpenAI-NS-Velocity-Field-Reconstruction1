"""Fresh synthetic checks for the Q4 frontier-closure audit.

All rows and evidence in this module are synthetic.  The fixture does not
import or load the original source; it is inert until ``run_fixture`` is
called with a supplied audit callable (or module exposing ``audit``).
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


SCHEMA = "original_q4_frontier_closure_audit_fixture_v1"
TARGET_NODE_IDS = (
    35559, 35674, 35675, 35676, 35677, 35678, 35679, 35680, 35558, 20
)
PROTECTED_NODE_IDS = frozenset(
    {0, 1, 17, 33572, 25220, 263344, 263348, 263349, 262918, 205147, 263346, 33982}
)
_VALID_STATUSES = frozenset({"conditionally_provable", "unresolved", "mismatch"})


def _fresh_rows() -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    for node_id in TARGET_NODE_IDS:
        edges = [] if node_id in {35558, 20} else [35558, 20]
        rows[node_id] = {
            "node_id": node_id,
            "metadata": {
                "operation": "synthetic_q4_frontier_component",
                "quantity": "synthetic_Phi",
                "synthetic_row_id": f"fresh-{node_id}",
            },
            "declared_edges": [{"node_id": edge, "role": "synthetic_child"} for edge in edges],
            "declared_edges_complete": True,
            "scope_aware": True,
            "source_binding_identity": {
                "binder_node_id": 35558,
                "binder_name": "synthetic_q4_rho",
                "binding_relation": "synthetic_occurrence_to_binder",
                "identity_proved": True,
                # Equality of this field must not be promoted to proof.
                "provider_refs": ("synthetic-provider",),
            },
            "occurrence_evidence": {
                "complete": True,
                "scope_aware": True,
                "authenticated": True,
                "binding_identity_proved": True,
            },
            "free_variable_evidence": {
                "complete": True,
                "scope_aware": True,
                "authenticated": True,
                "free_variable_set_complete": True,
                "target_binder_absent": True,
            },
            "source_authentication": {
                "authenticated": True,
                "basis": "synthetic-caller-claim-only",
            },
        }
    return rows


def _adapter(node_id: int, row: Any) -> Mapping[str, Any]:
    if not isinstance(row, Mapping) or row.get("node_id") != node_id:
        raise ValueError("synthetic row identity mismatch")
    return row


def _status(result: Any) -> str | None:
    return result.get("status") if isinstance(result, Mapping) and isinstance(result.get("status"), str) else None


def _scope(result: Any) -> Mapping[str, Any] | None:
    value = result.get("scope") if isinstance(result, Mapping) else None
    return value if isinstance(value, Mapping) else None


def _counts(result: Any) -> Mapping[str, Any] | None:
    value = result.get("counts") if isinstance(result, Mapping) else None
    return value if isinstance(value, Mapping) else None


def _check_common(result: Any) -> dict[str, bool]:
    scope = _scope(result)
    counts = _counts(result)
    return {
        "status_enum": _status(result) in _VALID_STATUSES,
        "scope_shape": scope is not None,
        "target_count": counts is not None and counts.get("target_ids") == len(TARGET_NODE_IDS),
        "source_auth_not_proved": scope is not None and scope.get("source_authentication_proved") is False,
        "graph_not_traversed": (
            scope is not None
            and scope.get("source_graph_traversed") is False
            and scope.get("declared_edges_traversed") is False
        ),
        "numeric_values_not_computed": (
            scope is not None and scope.get("numeric_field_value_computed") is False
        ),
        "provider_text_not_proof": (
            scope is not None and scope.get("provider_refs_or_text_used_as_proof") is False
        ),
    }


def _diagnostic(name: str, result: Any, checks: Mapping[str, bool], exc: str | None) -> dict[str, Any]:
    item = {
        "case": name,
        "status": _status(result),
        "passed": all(value is True for value in checks.values()),
        "failed_checks": [key for key, value in checks.items() if value is not True],
    }
    if exc is not None:
        item["exception"] = exc
    return item


def _case(audit: Callable[..., Any], name: str) -> tuple[Any, dict[str, bool], str | None]:
    rows = _fresh_rows()
    kwargs: dict[str, Any] = {"adapter": _adapter}
    expected = "conditionally_provable"
    if name == "positive_conditional":
        pass
    elif name == "raw_schema_default_unresolved":
        # Deliberately use a source-shaped but noncanonical row.  Without an
        # adapter the leaf must not guess that rho/Z fields are edges or that
        # provider metadata proves lexical closure.
        rows = {
            node_id: {
                "operation": "synthetic_unknown_source_operation",
                "rho": 35558,
                "Z_variable": 20,
                "provider_refs": ("same-provider",),
                "text": "same text",
            }
            for node_id in TARGET_NODE_IDS
        }
        kwargs = {}
        expected = "unresolved"
    elif name == "missing_row":
        rows.pop(35678)
        expected = "unresolved"
    elif name == "external_edge":
        rows[35559]["declared_edges"].append({"node_id": 999001, "role": "opaque_external"})
        expected = "unresolved"
    elif name == "protected_edge":
        rows[35674]["declared_edges"].append({"node_id": 0, "role": "opaque_protected"})
        expected = "unresolved"
    elif name == "provider_text_only":
        for row in rows.values():
            row["source_binding_identity"] = {"provider_refs": ("same",), "text": "same"}
            row["occurrence_evidence"] = None
            row["free_variable_evidence"] = None
        expected = "unresolved"
    elif name == "incomplete_free_variables":
        rows[35559]["free_variable_evidence"]["complete"] = False
        expected = "unresolved"
    elif name == "wrong_adapter_id":
        def wrong_adapter(node_id: int, row: Any) -> Mapping[str, Any]:
            value = dict(_adapter(node_id, row))
            value["node_id"] = node_id + 1
            return value
        kwargs = {"adapter": wrong_adapter}
        expected = "mismatch"
    elif name == "extra_protected_mapping_key":
        rows[0] = {"node_id": 0, "opaque": True}
    else:
        raise ValueError("unknown synthetic frontier case: " + name)

    result: Any = None
    exc: str | None = None
    try:
        result = audit(rows, **kwargs)
    except Exception as error:  # A raised error is not a returned contract.
        exc = type(error).__name__
    checks = _check_common(result)
    checks["expected_status"] = _status(result) == expected
    counts = _counts(result)
    if name == "positive_conditional":
        checks["all_rows_provable"] = counts is not None and counts.get("rows_conditionally_provable") == 10
        scope = _scope(result)
        checks["conditional_scope"] = (
            scope is not None
            and scope.get("closure_provable") is True
            and scope.get("conditional_on_caller_source_authentication") is True
        )
    elif name == "raw_schema_default_unresolved":
        checks["adapter_required"] = counts is not None and counts.get("rows_conditionally_provable") == 0
    elif name == "missing_row":
        checks["missing_counted"] = counts is not None and counts.get("rows_missing") == 1
    elif name in {"external_edge", "protected_edge"}:
        frontier = result.get("frontier") if isinstance(result, Mapping) else None
        checks["opaque_frontier"] = isinstance(frontier, list) and bool(frontier)
        checks["edges_not_traversed"] = all(
            item.get("traversed") is False for item in frontier if isinstance(item, Mapping)
        ) if isinstance(frontier, list) else False
    elif name == "provider_text_only":
        checks["same_strings_not_proof"] = _status(result) == "unresolved"
    elif name == "incomplete_free_variables":
        checks["free_variable_gap_unresolved"] = _status(result) == "unresolved"
    elif name == "wrong_adapter_id":
        checks["wrong_id_rejected"] = _status(result) == "mismatch"
    elif name == "extra_protected_mapping_key":
        access = result.get("access") if isinstance(result, Mapping) else None
        checks["extra_key_not_indexed"] = (
            isinstance(access, Mapping)
            and 0 not in list(access.get("indexed_ids", []))
            and set(access.get("indexed_ids", [])) == set(TARGET_NODE_IDS)
        )
    if exc is not None:
        checks["no_exception"] = False
    return result, checks, exc


def run_fixture(leaf_or_audit: Any) -> dict[str, Any]:
    """Run synthetic positive and adversarial checks against a supplied leaf."""

    audit = leaf_or_audit if callable(leaf_or_audit) else getattr(leaf_or_audit, "audit", None)
    if not callable(audit):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"audit_callable": False},
            "diagnostics": [{"case": "input", "passed": False}],
            "case_count": 0,
            "test_count": 1,
            "synthetic_only": True,
            "source_code_executed": False,
            "source_graph_walked": False,
        }

    cases = (
        "positive_conditional",
        "raw_schema_default_unresolved",
        "missing_row",
        "external_edge",
        "protected_edge",
        "provider_text_only",
        "incomplete_free_variables",
        "wrong_adapter_id",
        "extra_protected_mapping_key",
    )
    checks: dict[str, bool] = {"audit_callable": True}
    diagnostics: list[dict[str, Any]] = []
    for name in cases:
        result, case_checks, exc = _case(audit, name)
        checks.update({f"{name}_{key}": value for key, value in case_checks.items()})
        diagnostics.append(_diagnostic(name, result, case_checks, exc))
    return {
        "schema": SCHEMA,
        "passed": bool(checks) and all(value is True for value in checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "case_count": len(cases),
        "test_count": len(checks),
        "synthetic_only": True,
        "source_code_executed": False,
        "source_graph_walked": False,
        "numeric_values_computed": False,
        "scope": {
            "numeric_field_value_computed": False,
            "source_authentication_proved": False,
            "complete_identity_proved": False,
            "provider_refs_or_text_used_as_proof": False,
            "declared_edges_traversed": False,
        },
    }


run = run_fixture

__all__ = ["SCHEMA", "TARGET_NODE_IDS", "PROTECTED_NODE_IDS", "run_fixture", "run"]
