"""Independent synthetic fixture for the native Q4 endpoint consumer.

The fixture creates tiny local ``FunctionTransportGraph``, ``FunctionRef`` and
``FunctionJet`` objects.  It does not import a producer or rebuild an archive.
"""

from __future__ import annotations

from typing import Any


SCHEMA = "original_q4_lower_substitution_fixture_v1"


class FunctionRef:
    """Synthetic equivalent of native ``FunctionRef(graph, node)``."""

    def __init__(self, graph: Any, node: int):
        self.graph = graph
        self.node = node


class FunctionTransportGraph:
    def __init__(self):
        self.nodes: list[dict[str, Any]] = []

    def insert(self, row: dict[str, Any]) -> FunctionRef:
        node = len(self.nodes)
        self.nodes.append(dict(row))
        return FunctionRef(self, node)


class FunctionJet:
    def __init__(self, *rows: FunctionRef):
        self.rows = tuple(rows)


class _Owner:
    pass


def _node(graph: FunctionTransportGraph, operation: str, **fields: Any) -> FunctionRef:
    row = {"operation": operation}
    row.update(fields)
    return graph.insert(row)


def _make_case(case: str | None = None) -> tuple[Any, FunctionTransportGraph, FunctionTransportGraph]:
    """Build a fresh Q4-like live object for one positive/negative case."""

    graph = FunctionTransportGraph()
    other_graph = FunctionTransportGraph()
    sc = _node(graph, "current_original_source_parameter", name="sc",
               Z_independent=True, source_family="synthetic-native-family")
    half = _node(graph, "exact_rational", numerator=1, denominator=2)
    lower = _node(graph, "product", arguments=[sc.node, half.node])
    upper = _node(graph, "exact_rational", numerator=1, denominator=1)
    coordinate = _node(graph, "bound_variable", name="first_micro_coordinate")
    boundary_body0 = _node(graph, "product", arguments=[coordinate.node, sc.node])
    boundary_body1 = _node(graph, "sum", arguments=[boundary_body0.node, half.node])
    boundary_rows = (FunctionRef(graph, boundary_body0.node),
                     FunctionRef(graph, boundary_body1.node))
    source_chart: dict[str, Any] = {
        "coordinate": FunctionRef(graph, coordinate.node),
        "domain": [FunctionRef(graph, lower.node), FunctionRef(graph, upper.node)],
    }
    loop_first = source_chart
    first_micro: dict[str, Any] = {
        "source_chart": source_chart,
        "finite_boundary": {"R_Q_y": FunctionJet(*boundary_rows)},
    }
    loop = _Owner()
    loop.g = graph
    loop.sc = FunctionRef(graph, sc.node)
    loop.charts = {"first_micro": loop_first}
    q4 = _Owner()
    q4.g = graph
    q4.family = "synthetic-native-family"
    q4.loop = loop
    q4.charts = {"first_micro": first_micro}

    substitution_rows: list[FunctionRef] = []
    for boundary_ref in boundary_rows:
        substitution_rows.append(graph.insert({
            "operation": "function_substitution",
            "expression": boundary_ref.node,
            "variable": coordinate.node,
            "value": lower.node,
            "Z_independent_substitution": True,
        }))

    if case == "wrong_variable":
        substitution_rows[0] = graph.insert({
            "operation": "function_substitution", "expression": boundary_rows[0].node,
            "variable": upper.node, "value": lower.node,
            "Z_independent_substitution": True,
        })
    elif case == "wrong_value":
        substitution_rows[0] = graph.insert({
            "operation": "function_substitution", "expression": boundary_rows[0].node,
            "variable": coordinate.node, "value": upper.node,
            "Z_independent_substitution": True,
        })
    elif case == "wrong_body":
        substitution_rows[0] = graph.insert({
            "operation": "function_substitution", "expression": boundary_rows[1].node,
            "variable": coordinate.node, "value": lower.node,
            "Z_independent_substitution": True,
        })
    elif case == "z_dependent":
        substitution_rows[0] = graph.insert({
            "operation": "function_substitution", "expression": boundary_rows[0].node,
            "variable": coordinate.node, "value": lower.node,
            "Z_independent_substitution": False,
        })
    elif case == "wrong_graph":
        other_graph.nodes.extend(dict(row) for row in graph.nodes)
        substitution_rows[0] = FunctionRef(other_graph, substitution_rows[0].node)
    elif case == "same_numeric_different_graph":
        other_graph.nodes.extend(dict(row) for row in graph.nodes)
        substitution_rows[1] = FunctionRef(other_graph, substitution_rows[1].node)
    elif case == "source_chart_identity":
        first_micro["source_chart"] = {
            "coordinate": FunctionRef(graph, coordinate.node),
            "domain": [FunctionRef(graph, lower.node), FunctionRef(graph, upper.node)],
        }
    elif case == "source_node_out_of_range":
        substitution_rows[0] = FunctionRef(graph, len(graph.nodes) + 17)
    elif case == "sc_not_independent":
        graph.nodes[sc.node]["Z_independent"] = False
    elif case == "sc_metadata_missing":
        graph.nodes[sc.node].pop("Z_independent", None)
    elif case == "sc_family_mismatch":
        graph.nodes[sc.node]["source_family"] = "other-family"
    elif case == "wrong_lower_product":
        wrong_lower = _node(graph, "product", arguments=[half.node, sc.node])
        source_chart["domain"][0] = FunctionRef(graph, wrong_lower.node)
    elif case == "wrong_half":
        wrong_half = _node(graph, "exact_rational", numerator=1, denominator=3)
        wrong_lower = _node(graph, "product", arguments=[sc.node, wrong_half.node])
        source_chart["domain"][0] = FunctionRef(graph, wrong_lower.node)

    q4.initial = {"R_Q_y": FunctionJet(*substitution_rows)}
    return q4, graph, other_graph


def _json_safe(value: Any, seen: set[int] | None = None) -> bool:
    if isinstance(value, (str, int, float, bool)) or value is None:
        return True
    if seen is None:
        seen = set()
    identity = id(value)
    if identity in seen:
        return False
    seen.add(identity)
    if isinstance(value, dict):
        return all(isinstance(key, str) and _json_safe(child, seen)
                   for key, child in value.items())
    if isinstance(value, (list, tuple)):
        return all(_json_safe(child, seen) for child in value)
    return False


def _digest_shape(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 64:
        return False
    return all(char in "0123456789abcdef" for char in value)


def _failed_checks(result: Any) -> list[str]:
    if not isinstance(result, dict):
        return ["result_not_mapping"]
    return [name for name, passed in result.get("checks", {}).items() if passed is not True]


def _case_diagnostic(case: str, result: Any) -> dict[str, Any]:
    if not isinstance(result, dict):
        return {"case": case, "passed": False, "reason": "verifier did not return metadata"}
    return {
        "case": case,
        "passed": result.get("passed") is True,
        "status": result.get("status"),
        "failed_checks": _failed_checks(result)[:8],
    }


def _archive_fixture(consumer: Any) -> tuple[bool, dict[str, bool], dict[str, Any]]:
    keys = (
        "actual_common_prefix_and_native_initial_flux_ports",
        "actual_native_RQy_product_and_Z_independent_pullback",
        "original_loop_flat_branch_is_radial_not_a_Z_sample",
    )
    report = {
        "component_result": {
            "ordinary_Z_orders": [0, 1],
            "source_owned_Z_cell": ["0", "1/2"],
            "source_prerequisite_audit": {key: True for key in keys},
        }
    }
    bad_report = {
        "component_result": {
            "ordinary_Z_orders": [0, 1],
            "source_owned_Z_cell": ["0", "1/2"],
            "source_prerequisite_audit": {
                **{key: True for key in keys}, keys[1]: False,
            },
        }
    }
    good = consumer(report)
    bad = consumer(bad_report)
    checks = {
        "archive_positive_passed": isinstance(good, dict) and good.get("passed") is True,
        "archive_positive_scope_preserved": (
            isinstance(good, dict)
            and good.get("scope", {}).get("ordinary_Z_orders") == [0, 1]
            and good.get("scope", {}).get("source_owned_Z_cell") == ["0", "1/2"]
        ),
        "archive_negative_rejected": isinstance(bad, dict) and bad.get("passed") is False,
        "archive_no_live_identity_claim": (
            isinstance(good, dict)
            and good.get("live_q4_identity_proved") is False
            and good.get("corrected_source_identity_proved") is False
            and good.get("cross_owner_identity_claimed") is False
        ),
    }
    diagnostics = {
        "positive": _case_diagnostic("archived_positive", good),
        "negative": _case_diagnostic("archived_negative", bad),
    }
    return all(checks.values()), checks, diagnostics


def run_fixture(verify: Any, consume_archived_evidence: Any = None) -> dict[str, Any]:
    """Run positive and fail-closed synthetic live-Q4 checks."""

    checks: dict[str, bool] = {}
    diagnostics: list[dict[str, Any]] = []
    if not callable(verify):
        return {
            "schema": SCHEMA,
            "passed": False,
            "checks": {"verify_callable": False},
            "diagnostics": [{"case": "input", "passed": False,
                             "reason": "verify is not callable"}],
        }
    checks["verify_callable"] = True

    q4, _, _ = _make_case()
    positive = verify(q4)
    checks["positive_passed"] = isinstance(positive, dict) and positive.get("passed") is True
    checks["positive_verified_status"] = (
        isinstance(positive, dict) and positive.get("status") == "verified"
    )
    checks["positive_two_ordinary_rows"] = (
        isinstance(positive, dict)
        and positive.get("ordinary_Z_orders") == [0, 1]
        and positive.get("checks", {}).get("both_ordinary_Z_rows_native_substitution_exact") is True
    )
    checks["positive_Z0_Z1_matched"] = (
        isinstance(positive, dict)
        and set(positive.get("matched", {})) == {"Z0", "Z1"}
    )
    checks["positive_exact_source_literal_metadata"] = (
        isinstance(positive, dict)
        and positive.get("matched", {}).get("Z0", {}).get("expression") is not None
        and positive.get("matched", {}).get("Z1", {}).get("expression") is not None
        and positive.get("matched", {}).get("Z0", {}).get("variable") is not None
        and positive.get("matched", {}).get("Z0", {}).get("value") is not None
    )
    checks["positive_hashes_are_runtime_sha256"] = (
        isinstance(positive, dict)
        and all(
            _digest_shape(positive.get("matched", {}).get(label, {}).get("initial_node_sha256"))
            and _digest_shape(positive.get("matched", {}).get(label, {}).get("boundary_node_sha256"))
            for label in ("Z0", "Z1")
        )
        and positive.get("hash_checks") == {"Z0": True, "Z1": True}
    )
    checks["positive_sc_over_two_definition"] = (
        isinstance(positive, dict)
        and positive.get("checks", {}).get("first_micro_lower_is_sc_over_two_exact") is True
        and positive.get("source_ids", {}).get("sc") is not None
        and positive.get("source_ids", {}).get("half") is not None
    )
    checks["positive_metadata_is_json_safe"] = _json_safe(positive)
    checks["positive_no_numeric_evaluation"] = (
        isinstance(positive, dict)
        and positive.get("numeric_values_computed") is False
        and positive.get("graph_traversed") is False
    )
    checks["positive_no_cross_owner_claim"] = (
        isinstance(positive, dict)
        and positive.get("cross_owner_identity_claimed") is False
        and positive.get("corrected_source_identity_proved") is False
    )
    diagnostics.append(_case_diagnostic("positive", positive))

    negative_cases = (
        ("wrong_variable", "wrong_variable_literal_rejected"),
        ("wrong_value", "wrong_value_literal_rejected"),
        ("wrong_body", "wrong_body_literal_rejected"),
        ("wrong_graph", "wrong_graph_rejected"),
        ("z_dependent", "Z_dependence_rejected"),
        ("source_chart_identity", "source_chart_identity_rejected"),
        ("source_node_out_of_range", "source_node_out_of_range_rejected"),
        ("same_numeric_different_graph", "same_numeric_different_graph_rejected"),
        ("sc_not_independent", "sc_Z_dependence_rejected"),
        ("sc_metadata_missing", "missing_sc_Z_metadata_unresolved"),
        ("sc_family_mismatch", "sc_family_mismatch_rejected"),
        ("wrong_lower_product", "wrong_sc_product_rejected"),
        ("wrong_half", "wrong_half_rejected"),
    )
    for case, check_name in negative_cases:
        bad = verify(_make_case(case)[0])
        checks[check_name] = isinstance(bad, dict) and bad.get("passed") is False
        checks[case + "_status_unresolved"] = (
            isinstance(bad, dict) and bad.get("status") == "unresolved"
        )
        checks[case + "_no_cross_owner_claim"] = (
            isinstance(bad, dict)
            and bad.get("cross_owner_identity_claimed") is False
            and bad.get("corrected_source_identity_proved") is False
        )
        diagnostics.append(_case_diagnostic(case, bad))

    if consume_archived_evidence is None:
        checks["archive_consumer_optional"] = True
    elif callable(consume_archived_evidence):
        archive_passed, archive_checks, archive_diagnostics = _archive_fixture(
            consume_archived_evidence
        )
        checks.update(archive_checks)
        checks["archive_fixture_passed"] = archive_passed
        diagnostics.append({"case": "archived_metadata", **archive_diagnostics})
    else:
        checks["archive_consumer_optional"] = False
        diagnostics.append({"case": "archived_metadata", "passed": False,
                            "reason": "archive consumer is not callable"})

    return {
        "schema": SCHEMA,
        "passed": all(checks.values()),
        "checks": checks,
        "diagnostics": diagnostics,
        "synthetic_only": True,
        "source_code_executed": False,
        "source_graph_walked": False,
        "source_values_evaluated": False,
        "old_owner_constructed": False,
        "old_helper_called": False,
        "cross_owner_identity_claimed": False,
    }


run = run_fixture

__all__ = [
    "SCHEMA", "FunctionRef", "FunctionTransportGraph", "FunctionJet",
    "run_fixture", "run",
]

