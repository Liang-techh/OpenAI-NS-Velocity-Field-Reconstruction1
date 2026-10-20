"""Bind the actual finite-``N`` nonreset memory rows.

This is a small source/data consumer for the frozen first-ten history attempt.
It authenticates one saved component result and the two graph archives, then
checks the fifty symbolic memory factors.  It deliberately does not import a
project module, execute a producer or operator, reduce a field, or recompute
the finite-``N`` history.
"""

from __future__ import annotations

import base64
import copy
import gzip
import hashlib
import json
from collections.abc import Mapping, Sequence
from fractions import Fraction
from typing import Any


SCHEMA = "original_finite_N_nonreset_memory_source_binding_v1"
OLD_SCHEMA = "original_finite_N_nonreset_history_source_binding_v1"
OLD_STAGE_SCHEMA = "original_finite_N_nonreset_history_source_binding_stage_v1"

EXPECTED_NODE_COUNT = 205300
EXPECTED_NODES_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"

PREFIX = "lei_ren_part1_paper_compliant_current_original_"
PREVIOUS_REPORT = PREFIX + "finite_N_nonreset_history_source_binding.json"
PREVIOUS_REPORT_SHA256 = "f42277a1c5001573e90ae7f7a6915ad0da6d0008a7229cd6befda867246114b5"
CONTROL_ARCHIVE = PREFIX + "whole_Z_all_N_control_connection.json.gz"
CONTROL_ARCHIVE_SHA256 = "588f411e9122b9a76e6a42efc2b5754137f6af66a9cae42f8d6a7c2c59920de8"
C5_ARCHIVE = PREFIX + "prescribed_loop_C5_functions.json.gz"
C5_ARCHIVE_SHA256 = "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"

OLD_HISTORY_SOURCE = "src/openai_ns_reconstruction/original_finite_N_nonreset_history_source_binding.py"
OLD_HISTORY_SOURCE_SHA256 = "ae2b95d3effa1a1ed43ce9e1004c7ffa5be6f16394430dc534b2f58049ce4d22"
OLD_HISTORY_RUNNER = "scripts/run_original_finite_N_nonreset_history_source_binding.py"
OLD_HISTORY_RUNNER_SHA256 = "34d96d41421dc9eee267af6671e1737d5e364a5dbe2031afac1128e18c7cdc3c"
LOADER_SOURCE = "scripts/diagnose_original_Rm_physical_qy0_access.py"
LOADER_SOURCE_SHA256 = "765667cb242c4dbb21fb7ef5ed030b6d4d6889b48136cbff5cf8475fd8967fdb"
PARENT_REPORT = PREFIX + "Rm_physical_Qy_Z_calculus.json"
PARENT_REPORT_SHA256 = "8903b59267bae01d5182bf2dcc84ed3c02c44eb8c7d3dd9d7b14af67ac092746"
PROVIDER_REPORT = PREFIX + "finite_N_provider_source_binding.json"
PROVIDER_REPORT_SHA256 = "1f23b6188abffb6d99cf00fb1274a58a2fba1100a326334a1a798d0b920f6c41"

CHECKPOINT_PREFIX = "tmp/original_finite_N_nonreset_history_source_binding/new_history_binding_"
CHECKPOINT_PATHS = (
    CHECKPOINT_PREFIX + "source.json",
    CHECKPOINT_PREFIX + "entered.json",
    CHECKPOINT_PREFIX + "result.json",
)

RECORD_HASHES = {
    PREVIOUS_REPORT: PREVIOUS_REPORT_SHA256,
    CONTROL_ARCHIVE: CONTROL_ARCHIVE_SHA256,
    C5_ARCHIVE: C5_ARCHIVE_SHA256,
}

FIRST_TEN_CHARTS = (
    "first_micro",
    "second_micro",
    "frozen_macro",
    "first_switch",
    "second_switch",
    "post_power",
    "long_reshape",
    "reference",
    "restoration",
    "postrestore",
)
RATES = ("m", "h", "k", "e", "p")
RATE_VALUES = {
    "m": Fraction(1),
    "h": Fraction(3, 2),
    "k": Fraction(3, 2),
    "e": Fraction(1),
    "p": Fraction(0),
}
ORDINARY_Z_ORDERS = (0, 1, 2)

SOURCE_FAMILY = {
    "actual_five_defect_family_sha256": "3983d0ddb33fa85e6ab152ef7e29960f8b95aca3e1e86f1bda0882b39d825894",
    "implicit_source_sha256": "5aec111986d745459eb2e2fece291f1dc3fa7bf986529494df802c2aa2daceae",
    "datum_enclosure_sha256": "dd4040ee3ded75a56f65407fdd1cebee6aa6f559cce9d1bb57e2146a8a436010",
}

OLD_MEMORY_GATES = tuple(
    "actual_own_rate_memory_unproved:" + chart for chart in FIRST_TEN_CHARTS
)
OLD_NEGATIVE_GATES = OLD_MEMORY_GATES + ("first_ten_nonreset_history_identity_unproved",)


def _sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def canonical_digest(value: Any) -> str:
    """Hash JSON graph data without interpreting any graph node."""

    return _sha(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        ).encode("utf-8")
    )


def _gate(gates: list[str], name: str) -> None:
    if name not in gates:
        gates.append(name)


def _valid_index(graph: Any, index: Any) -> bool:
    return (
        isinstance(graph, Sequence)
        and not isinstance(graph, (str, bytes, bytearray))
        and isinstance(index, int)
        and not isinstance(index, bool)
        and 0 <= index < len(graph)
    )


def _product_signature(
    graph: Sequence[Any], index: Any, *, depth: int = 0
) -> tuple[Fraction, tuple[int, ...]] | None:
    """Normalize products and exact rational signs, leaving other nodes opaque.

    A graph node is never evaluated.  Only ``exact_rational``, ``negative``
    and ``product`` are opened; every other mapping is an opaque factor handle.
    """

    if depth > 256 or not _valid_index(graph, index):
        return None
    node = graph[index]
    if not isinstance(node, Mapping):
        return None
    operation = node.get("operation")
    if operation == "exact_rational":
        numerator = node.get("numerator")
        denominator = node.get("denominator")
        if (
            type(numerator) is not int
            or type(denominator) is not int
            or denominator == 0
        ):
            return None
        try:
            return Fraction(numerator, denominator), ()
        except (TypeError, ValueError, ZeroDivisionError):
            return None
    if operation == "negative":
        value = _product_signature(graph, node.get("argument"), depth=depth + 1)
        if value is None:
            return None
        return -value[0], value[1]
    if operation == "product":
        arguments = node.get("arguments")
        if not isinstance(arguments, list):
            return None
        coefficient = Fraction(1)
        factors: list[int] = []
        for argument in arguments:
            value = _product_signature(graph, argument, depth=depth + 1)
            if value is None:
                return None
            coefficient *= value[0]
            factors.extend(value[1])
        return coefficient, tuple(sorted(factors))
    # Analytic functions, sums, parameters, and every other graph node are
    # opaque source-port handles in this bounded matcher.
    return Fraction(1), (index,)


def _difference_matches(
    graph: Sequence[Any], index: Any, left: Any, right: Any
) -> bool:
    """Recognize ``left - right`` from source graph structure only."""

    right_form = _product_signature(graph, right)
    if right_form == (Fraction(0), ()):
        return index == left
    if not _valid_index(graph, index):
        return False
    node = graph[index]
    if not isinstance(node, Mapping) or node.get("operation") != "sum":
        return False
    arguments = node.get("arguments")
    if not isinstance(arguments, list) or len(arguments) != 2:
        return False
    return bool(
        arguments[0] == left
        and right_form is not None
        and _product_signature(graph, arguments[1])
        == (-right_form[0], right_form[1])
    )


def _signature_json(
    value: tuple[Fraction, tuple[int, ...]] | None,
) -> dict[str, Any] | None:
    if value is None:
        return None
    coefficient, factors = value
    return {"coefficient": str(coefficient), "factors": list(factors)}


def _json(text: str | None, gates: list[str], name: str) -> Mapping[str, Any]:
    if not isinstance(text, str):
        _gate(gates, f"source_record_text_missing:{name}")
        return {}
    try:
        value = json.loads(text)
    except (TypeError, ValueError, UnicodeError):
        _gate(gates, f"json_invalid:{name}")
        return {}
    if not isinstance(value, Mapping):
        _gate(gates, f"json_not_mapping:{name}")
        return {}
    return value


def _record_text(
    source_records: Mapping[str, Any], name: str, gates: list[str]
) -> str | None:
    record = source_records.get(name)
    if not isinstance(record, Mapping):
        _gate(gates, f"source_record_missing:{name}")
        return None
    text = record.get("source_text")
    if (
        record.get("path") != "experiments/root_st073/" + name
        or record.get("sha256") != RECORD_HASHES[name]
        or not isinstance(text, str)
        or _sha(text.encode("utf-8")) != RECORD_HASHES[name]
    ):
        _gate(gates, f"source_record_contract_mismatch:{name}")
        return None
    return text


def _archive_text(
    source_records: Mapping[str, Any], name: str, gates: list[str]
) -> str | None:
    """Authenticate compressed bytes before accepting any decoded text."""

    record = source_records.get(name)
    if not isinstance(record, Mapping):
        _gate(gates, f"source_record_missing:{name}")
        return None
    text = record.get("source_text")
    encoded = record.get("raw_bytes_base64")
    if (
        record.get("path") != "experiments/root_st073/" + name
        or record.get("sha256") != RECORD_HASHES[name]
        or record.get("encoding") != "gzip"
        or not isinstance(text, str)
        or not isinstance(encoded, str)
    ):
        _gate(gates, f"archive_record_contract_mismatch:{name}")
        return None
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (TypeError, ValueError):
        _gate(gates, f"raw_archive_base64_invalid:{name}")
        return None
    if _sha(raw) != RECORD_HASHES[name] or not raw.startswith(b"\x1f\x8b"):
        _gate(gates, f"raw_archive_hash_mismatch:{name}")
        return None
    try:
        decoded = gzip.decompress(raw)
        decoded_text = decoded.decode("utf-8")
    except (OSError, EOFError, UnicodeError):
        _gate(gates, f"raw_archive_decode_error:{name}")
        return None
    if _sha(decoded) != record.get("decoded_sha256") or decoded_text != text:
        _gate(gates, f"decoded_archive_mismatch:{name}")
        return None
    return text


def _all_true(value: Any) -> bool:
    """Require nonempty structural evidence with only true boolean leaves."""

    if isinstance(value, Mapping):
        return bool(value) and all(_all_true(item) for item in value.values())
    if isinstance(value, list):
        return bool(value) and all(_all_true(item) for item in value)
    if isinstance(value, bool):
        return value
    return False


def _prior_checkpoint_provenance(report: Mapping[str, Any]) -> bool:
    hashes = report.get("checkpoint_hashes")
    checkpoints = report.get("checkpoints")
    if (
        not isinstance(hashes, Mapping)
        or not isinstance(checkpoints, list)
        or report.get("checkpoint_count") != len(CHECKPOINT_PATHS)
        or tuple(hashes) != CHECKPOINT_PATHS
        or len(checkpoints) != len(CHECKPOINT_PATHS)
    ):
        return False
    entries: dict[str, Mapping[str, Any]] = {}
    for item in checkpoints:
        if not isinstance(item, Mapping) or not isinstance(item.get("path"), str):
            return False
        entries[item["path"]] = item
    if tuple(entries) != CHECKPOINT_PATHS:
        return False
    for path in CHECKPOINT_PATHS:
        expected = hashes.get(path)
        item = entries.get(path)
        if not isinstance(expected, str) or not isinstance(item, Mapping):
            return False
        if item.get("sha256") != expected or not isinstance(item.get("bytes"), int):
            return False
        encoded = item.get("base64_bytes")
        if not isinstance(encoded, str):
            return False
        try:
            raw = base64.b64decode(encoded, validate=True)
        except (TypeError, ValueError):
            return False
        if len(raw) != item["bytes"] or _sha(raw) != expected:
            return False
    return True


def _prior_nonmemory_checks(
    report: Mapping[str, Any], component: Mapping[str, Any]
) -> tuple[bool, dict[str, bool]]:
    """Validate the saved negative attempt without replaying any old helper."""

    checks: dict[str, bool] = {}
    checks["outer_schema"] = report.get("schema") == OLD_STAGE_SCHEMA
    checks["outer_source_only"] = report.get("source_only") is True
    checks["outer_scope_flags"] = all(
        report.get(name) is expected
        for name, expected in (
            ("numeric_values_computed", False),
            ("old_operators_replayed", False),
            ("full_goal_complete", False),
            ("source_only_validation", True),
            ("saved_component_results_used", True),
            ("no_target_or_helper_execution_by_worker", True),
        )
    )
    checks["component_schema"] = component.get("schema") == OLD_SCHEMA
    checks["component_negative"] = (
        component.get("passed") is False
        and component.get("failed_attempt_preserved") is not True
    )
    checks["component_graph"] = (
        component.get("nodes_digest") == EXPECTED_NODES_DIGEST
        and component.get("expected_nodes_digest") == EXPECTED_NODES_DIGEST
        and component.get("node_count") == EXPECTED_NODE_COUNT
    )
    checks["component_nonmemory_flags"] = all(
        component.get(name) is expected
        for name, expected in (
            ("source_only", True),
            ("ordinary_Z_not_temporal_n", True),
            ("first_ten_density_C2_NC5_identity_proved", True),
            ("first_ten_finite_N_field_source_composition_proved", True),
            ("generic_density_identity_reused", True),
            ("generic_density_identity_proved", True),
            ("numeric_values_computed", False),
            ("old_operators_replayed", False),
            ("mathematical_operator_replayed", False),
            ("finite_N_history_identity", False),
            ("same_density_certificate", False),
            ("physical_Q_y_Z1", False),
            ("global_five_moment_closure", False),
            ("patch_source_composition_proved", False),
            ("full_goal_complete", False),
            ("time_order_recursion_claimed", False),
            ("timeorder_recursion_claimed", False),
        )
    ) and component.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS)
    checks["negative_gate_set"] = (
        isinstance(component.get("unresolved_gates"), list)
        and len(component["unresolved_gates"]) == len(OLD_NEGATIVE_GATES)
        and set(component["unresolved_gates"]) == set(OLD_NEGATIVE_GATES)
    )
    saved_ast = component.get("source_ast_witness")
    checks["source_ast_passed"] = isinstance(saved_ast, Mapping) and bool(saved_ast) and _all_true(saved_ast)
    saved_interiors = component.get("first_ten_interior_source_specs")
    checks["source_interior_definitions_reused"] = isinstance(saved_interiors, Mapping) and all(
        isinstance(saved_interiors.get(chart), Mapping) and saved_interiors[chart].get("passed") is True
        for chart in FIRST_TEN_CHARTS)

    partitions = component.get("partition_coverage")
    checks["native_partitions_passed"] = (
        isinstance(partitions, Mapping)
        and all(
            isinstance(partitions.get(chart), Mapping)
            and partitions[chart].get("passed") is True
            for chart in FIRST_TEN_CHARTS
        )
    )
    jacobians = component.get("physical_jacobian_audit")
    checks["jacobians_passed"] = (
        isinstance(jacobians, Mapping)
        and all(
            isinstance(jacobians.get(chart), Mapping)
            and jacobians[chart].get("passed") is True
            for chart in FIRST_TEN_CHARTS
        )
    )
    history = component.get("first_ten_history_ports")
    checks["source_history_passed"] = (
        isinstance(history, Mapping)
        and all(
            isinstance(history.get(chart), Mapping)
            and isinstance(history[chart].get("checks"), Mapping)
            and history[chart]["checks"].get("chart_passed") is True
            and all(flag is True for flag in history[chart]["checks"].values())
            for chart in FIRST_TEN_CHARTS
        )
    )

    endpoint = component.get("endpoint_ports")
    chart_ports = endpoint.get("chart_ports") if isinstance(endpoint, Mapping) else None
    endpoint_rows_ok = True
    if not isinstance(endpoint, Mapping) or endpoint.get("endpoint_rows_are_values") is not False:
        endpoint_rows_ok = False
    if not isinstance(chart_ports, Mapping):
        endpoint_rows_ok = False
    else:
        for chart in FIRST_TEN_CHARTS:
            row = chart_ports.get(chart)
            substitution = row.get("C5_endpoint_substitution_nodes_match_all_six_orders") if isinstance(row, Mapping) else None
            ports = row.get("endpoint_ports") if isinstance(row, Mapping) else None
            endpoint_rows_ok = endpoint_rows_ok and (
                isinstance(row, Mapping)
                and row.get("endpoint_rows_are_values") is False
                and row.get("C2_incoming_is_previous_outgoing") is True
                and row.get("C5_incoming_is_previous_endpoint_or_original_zero") is True
                and isinstance(substitution, Mapping)
                and all(substitution.get(rate) is True for rate in RATES)
                and row.get("becomes_next_incoming") is True
                and row.get("passed") is True
                and isinstance(ports, Mapping)
                and all(
                    isinstance(ports.get(rate), list)
                    and len(ports[rate]) == 6
                    for rate in RATES
                )
            )
    checks["endpoint_ports_passed"] = endpoint_rows_ok

    geometry = component.get("source_geometry_and_jacobians")
    geometry_charts = geometry.get("charts") if isinstance(geometry, Mapping) else None
    checks["source_geometry_reused"] = (
        isinstance(geometry, Mapping)
        and geometry.get("positive_widths_numerically_admitted") is False
        and isinstance(geometry_charts, Mapping)
        and all(
            isinstance(geometry_charts.get(chart), Mapping)
            and geometry_charts[chart].get("Jacobian_is_source_derivative") is True
            and geometry_charts[chart].get("rate_kernel_uses_same_offset") is True
            for chart in FIRST_TEN_CHARTS
        )
    )
    domain = component.get("domain_evidence")
    checks["conditional_domain_reused"] = (
        isinstance(domain, Mapping)
        and domain.get("source_only") is True
        and domain.get("conditional_original_common_admissible_domain") is True
        and domain.get("numerical_domain_admission_proved") is False
        and domain.get("endpoint_rows_are_values") is False
        and domain.get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS)
    )

    audit = component.get("source_inspection_audit")
    checks["prior_static_audit"] = (
        isinstance(audit, Mapping)
        and audit.get("source_records_only") is True
        and audit.get("nodes_hashed_only") is True
        and audit.get("new_source_consumer_calls") == 1
        and audit.get("accepted_ancestor_operator_calls") == 0
        and audit.get("accepted_mathematical_operator_calls") == 0
        and audit.get("old_mathematical_operator_calls") == 0
        and audit.get("history_operator_calls") == 0
        and audit.get("quadrature_or_norm_calls") == 0
        and audit.get("numeric_field_evaluated") is False
    )

    pins = report.get("source_pins")
    input_hashes = report.get("input_hashes")
    implementation_hashes = report.get("implementation_hashes")
    component_hashes = report.get("component_implementation_hashes")
    checks["source_pins"] = (
        isinstance(pins, Mapping)
        and isinstance(input_hashes, Mapping)
        and isinstance(implementation_hashes, Mapping)
        and isinstance(component_hashes, Mapping)
        and pins.get(OLD_HISTORY_SOURCE) == OLD_HISTORY_SOURCE_SHA256
        and pins.get(OLD_HISTORY_RUNNER) == OLD_HISTORY_RUNNER_SHA256
        and pins.get(LOADER_SOURCE) == LOADER_SOURCE_SHA256
        and implementation_hashes.get(OLD_HISTORY_SOURCE) == OLD_HISTORY_SOURCE_SHA256
        and implementation_hashes.get(OLD_HISTORY_RUNNER) == OLD_HISTORY_RUNNER_SHA256
        and component_hashes.get(OLD_HISTORY_SOURCE) == OLD_HISTORY_SOURCE_SHA256
        and component_hashes.get(OLD_HISTORY_RUNNER) == OLD_HISTORY_RUNNER_SHA256
        and component_hashes.get(LOADER_SOURCE) == LOADER_SOURCE_SHA256
        and all(input_hashes.get(name) == digest for name, digest in (
            (CONTROL_ARCHIVE, CONTROL_ARCHIVE_SHA256),
            (C5_ARCHIVE, C5_ARCHIVE_SHA256),
        ))
    )
    checks["checkpoint_provenance"] = _prior_checkpoint_provenance(report)
    checks["parent_provenance"] = (
        report.get("parent_report") == PARENT_REPORT
        and report.get("parent_report_sha256") == PARENT_REPORT_SHA256
        and report.get("preserved_parent_report_sha256") == PARENT_REPORT_SHA256
        and report.get("accepted_source_stage_report") == PROVIDER_REPORT
        and report.get("accepted_source_stage_sha256") == PROVIDER_REPORT_SHA256
    )
    checks["saved_scope"] = (
        isinstance(report.get("scope"), Mapping)
        and report["scope"].get("ordinary_Z_orders") == list(ORDINARY_Z_ORDERS)
        and report["scope"].get("ordinary_Z_not_temporal_n") is True
        and report["scope"].get("same_density_certificate") is False
        and report["scope"].get("finite_N_history_identity") is False
        and report["scope"].get("physical_Q_y_Z1") is False
        and report["scope"].get("global_five_moment_closure") is False
        and report["scope"].get("full_goal_complete") is False
    )
    checks["saved_result_flags"] = (
        isinstance(report.get("result_flags"), Mapping)
        and report["result_flags"].get("passed") is False
        and report["result_flags"].get("first_ten_nonreset_history_identity_proved") is False
        and report["result_flags"].get("first_ten_history_endpoint_identity_proved") is False
        and report["result_flags"].get("patch_source_composition_proved") is False
        and report["result_flags"].get("full_goal_complete") is False
    )
    checks["all_prior_nonmemory_premises"] = all(checks.values())
    return checks["all_prior_nonmemory_premises"], checks


def _exp_node(graph: Sequence[Any], index: Any) -> dict[str, Any] | None:
    if not _valid_index(graph, index) or not isinstance(graph[index], Mapping):
        return None
    node = graph[index]
    return {
        "index": index,
        "operation": node.get("operation"),
        "name": node.get("name"),
        "argument": node.get("argument"),
    }


def _is_exp_node(graph: Sequence[Any], index: Any) -> bool:
    node = _exp_node(graph, index)
    return bool(
        node is not None
        and node.get("operation") == "analytic_unary"
        and node.get("name") == "exp"
        and _valid_index(graph, node.get("argument"))
    )


def _memory_rate_audit(
    control: Mapping[str, Any], c5: Mapping[str, Any], gates: list[str]
) -> dict[str, Any]:
    """Check all fifty analytic memory rows without evaluating a graph node."""

    c2_graph = control.get("exact_function_graph_nodes")
    c5_graph = c5.get("exact_graph_nodes")
    c2_binding = control.get("exact_current_integral_and_control_graph")
    c5_binding = c5.get("prescribed_loop_C5_bindings")
    control_rows = c2_binding.get("windows") if isinstance(c2_binding, Mapping) else None
    charts = c5_binding.get("charts") if isinstance(c5_binding, Mapping) else None
    rows = {
        row.get("chart"): row
        for row in control_rows
        if isinstance(row, Mapping) and isinstance(row.get("chart"), str)
    } if isinstance(control_rows, list) else {}

    out: dict[str, Any] = {}
    pair_count = 0
    for chart in FIRST_TEN_CHARTS:
        old = rows.get(chart, {})
        current = charts.get(chart, {}) if isinstance(charts, Mapping) else {}
        width_index = old.get("width") if isinstance(old, Mapping) else None
        width_form = _product_signature(c2_graph, width_index) if isinstance(c2_graph, Sequence) else None
        control_memory = old.get("memory") if isinstance(old, Mapping) else None
        c5_memory = current.get("memory") if isinstance(current, Mapping) else None
        checks: dict[str, Any] = {}
        for field in RATES:
            rate = RATE_VALUES[field]
            c2_index = control_memory.get(field) if isinstance(control_memory, Mapping) else None
            c5_index = c5_memory.get(field) if isinstance(c5_memory, Mapping) else None
            c2_argument = None
            c2_exponent = None
            c2_expected = None
            c2_ok = False
            if rate == 0:
                c2_ok = _product_signature(c2_graph, c2_index) == (Fraction(1), ())
            elif isinstance(c2_graph, Sequence) and _is_exp_node(c2_graph, c2_index):
                c2_node = _exp_node(c2_graph, c2_index)
                c2_argument = c2_node.get("argument") if c2_node else None
                c2_exponent = _product_signature(c2_graph, c2_argument)
                if width_form is not None:
                    c2_expected = (-rate * width_form[0], width_form[1])
                    c2_ok = c2_exponent == c2_expected

            c5_argument = None
            c5_exponent = None
            c5_expected = None
            delta_index = None
            delta_form = None
            measure_form = None
            c5_ok = False
            if rate == 0:
                c5_ok = _product_signature(c5_graph, c5_index) == (Fraction(1), ())
            elif isinstance(c5_graph, Sequence) and _is_exp_node(c5_graph, c5_index):
                c5_node = _exp_node(c5_graph, c5_index)
                c5_argument = c5_node.get("argument") if c5_node else None
                c5_exponent = _product_signature(c5_graph, c5_argument)
                domain = current.get("domain") if isinstance(current, Mapping) else None
                coordinate = current.get("coordinate") if isinstance(current, Mapping) else None
                if isinstance(domain, list) and len(domain) == 2:
                    lower = domain[0]
                    # Find the symbolic difference port from the exponent's
                    # factors.  It is checked structurally, never evaluated.
                    candidates = list(c5_exponent[1]) if c5_exponent is not None else []
                    for candidate in candidates:
                        if _difference_matches(c5_graph, candidate, coordinate, lower):
                            delta_index = candidate
                            break
                    if delta_index is None and _difference_matches(c5_graph, coordinate, coordinate, lower):
                        delta_index = coordinate
                    if delta_index is not None:
                        delta_form = _product_signature(c5_graph, delta_index)
                    measure_form = _product_signature(
                        c5_graph,
                        current.get("physical_measure") if isinstance(current, Mapping) else None,
                    )
                    if measure_form is not None and delta_form is not None:
                        c5_expected = (
                            -rate * measure_form[0] * delta_form[0],
                            tuple(sorted(measure_form[1] + delta_form[1])),
                        )
                        c5_ok = c5_exponent == c5_expected

            pair_count += 1
            checks[field] = {
                "rate": str(rate),
                "C2_archive": CONTROL_ARCHIVE,
                "C5_archive": C5_ARCHIVE,
                "C2_memory_node": c2_index,
                "C2_width_node": width_index,
                "C2_exponential_node_schema": _exp_node(c2_graph, c2_index),
                "C2_exponent_signature": _signature_json(c2_exponent),
                "C2_expected_exponent_signature": _signature_json(c2_expected),
                "C5_memory_node": c5_index,
                "C5_exponential_node_schema": _exp_node(c5_graph, c5_index),
                "C5_native_difference_node": delta_index,
                "C5_physical_measure_node": current.get("physical_measure") if isinstance(current, Mapping) else None,
                "C5_exponent_signature": _signature_json(c5_exponent),
                "C5_expected_exponent_signature": _signature_json(c5_expected),
                "rate_zero_literal_identity_memory": rate == 0,
                "C2_actual_exponent_matches": bool(c2_ok),
                "C5_actual_exponent_matches": bool(c5_ok),
                "passed": bool(c2_ok and c5_ok),
            }
        row_passed = len(checks) == len(RATES) and all(
            row.get("passed") is True for row in checks.values()
        )
        if not row_passed:
            _gate(gates, "actual_own_rate_memory_unproved:" + chart)
        out[chart] = {
            "checks": checks,
            "passed": row_passed,
        }
    return out


def derive(
    nodes: Sequence[Any],
    source_records: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Compose the saved non-memory proof with a corrected memory audit."""

    gates: list[str] = []
    try:
        actual_digest = canonical_digest(nodes)
    except (TypeError, ValueError, OverflowError):
        actual_digest = ""
        _gate(gates, "nodes_digest_error")
    valid_nodes = isinstance(nodes, Sequence) and not isinstance(nodes, (str, bytes, bytearray))
    if expected_nodes_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "expected_nodes_digest_mismatch")
    if actual_digest != EXPECTED_NODES_DIGEST:
        _gate(gates, "nodes_digest_mismatch")
    if not valid_nodes or len(nodes) != EXPECTED_NODE_COUNT:
        _gate(gates, "node_count_mismatch")
    if not isinstance(source_records, Mapping):
        _gate(gates, "source_records_not_mapping")
        records: Mapping[str, Any] = {}
    else:
        records = source_records

    texts: dict[str, str] = {}
    previous_text = _record_text(records, PREVIOUS_REPORT, gates)
    control_text = _archive_text(records, CONTROL_ARCHIVE, gates)
    c5_text = _archive_text(records, C5_ARCHIVE, gates)
    if previous_text is not None:
        texts[PREVIOUS_REPORT] = previous_text
    if control_text is not None:
        texts[CONTROL_ARCHIVE] = control_text
    if c5_text is not None:
        texts[C5_ARCHIVE] = c5_text

    previous_report = _json(texts.get(PREVIOUS_REPORT), gates, PREVIOUS_REPORT)
    previous_component = previous_report.get("component_result")
    if not isinstance(previous_component, Mapping):
        previous_component = {}
        _gate(gates, "saved_prior_component_missing")
    prior_ok, prior_checks = _prior_nonmemory_checks(previous_report, previous_component)
    if not prior_ok:
        _gate(gates, "saved_nonmemory_prerequisites_unproved")

    control = _json(texts.get(CONTROL_ARCHIVE), gates, CONTROL_ARCHIVE)
    c5 = _json(texts.get(C5_ARCHIVE), gates, C5_ARCHIVE)
    memories = _memory_rate_audit(control, c5, gates)
    memory_ok = (
        isinstance(memories, Mapping)
        and len(memories) == len(FIRST_TEN_CHARTS)
        and all(
            isinstance(memories.get(chart), Mapping)
            and isinstance(memories[chart].get("checks"), Mapping)
            and len(memories[chart]["checks"]) == len(RATES)
            and memories[chart].get("passed") is True
            and all(
                isinstance(memories[chart]["checks"].get(field), Mapping)
                and memories[chart]["checks"][field].get("rate") == str(RATE_VALUES[field])
                and memories[chart]["checks"][field].get("C2_actual_exponent_matches") is True
                and memories[chart]["checks"][field].get("C5_actual_exponent_matches") is True
                and memories[chart]["checks"][field].get("passed") is True
                for field in RATES
            )
            for chart in FIRST_TEN_CHARTS
        )
    )
    if not memory_ok:
        _gate(gates, "actual_memory_pairs_unproved")

    overall_passed = bool(
        valid_nodes
        and len(nodes) == EXPECTED_NODE_COUNT
        and actual_digest == EXPECTED_NODES_DIGEST
        and expected_nodes_digest == EXPECTED_NODES_DIGEST
        and prior_ok
        and memory_ok
        and not gates
    )

    if isinstance(previous_component, Mapping):
        result: dict[str, Any] = copy.deepcopy(dict(previous_component))
        # Capture the exact old negative component before changing any outer
        # flags.  This is saved JSON data, not a re-entry into the old leaf.
        result["prior_component_result"] = copy.deepcopy(dict(previous_component))
    else:
        result = {}
    result.update(
        {
            "schema": SCHEMA,
            "source_only": True,
            "passed": overall_passed,
            "nodes_digest": actual_digest,
            "expected_nodes_digest": EXPECTED_NODES_DIGEST,
            "node_count": len(nodes) if valid_nodes else None,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "ordinary_Z_not_temporal_n": True,
            "time_order_recursion_claimed": False,
            "timeorder_recursion_claimed": False,
            "actual_own_rate_memory_audit": memories,
            "actual_memory_pair_count": sum(len(row.get("checks", {})) for row in memories.values()),
            "coveredcharts": list(FIRST_TEN_CHARTS) if overall_passed else [],
            "covered_charts": list(FIRST_TEN_CHARTS) if overall_passed else [],
            "first_ten_nonreset_history_identity_proved": overall_passed,
            "first_ten_history_endpoint_identity_proved": overall_passed,
            "record_hashes": dict(RECORD_HASHES),
            "unresolved_gates": sorted(set(gates)),
            "saved_nonmemory_prerequisites_reused": prior_ok,
            "generic_density_identity_reused": prior_ok,
            "numeric_values_computed": False,
            "old_operators_replayed": False,
            "mathematical_operator_replayed": False,
            "finite_N_history_identity": False,
            "same_density_certificate": False,
            "physical_Q_y_Z1": False,
            "global_five_moment_closure": False,
            "patch_source_composition_proved": False,
            "full_goal_complete": False,
            "runtime_model_metadata": None,
        }
    )
    endpoint_ports = result.get("endpoint_ports")
    if isinstance(endpoint_ports, Mapping):
        endpoint_ports = copy.deepcopy(dict(endpoint_ports))
        endpoint_ports["a_at_endpoint_becomes_next_incoming"] = overall_passed
        result["endpoint_ports"] = endpoint_ports

    history_proof = result.get("history_identity_proof")
    if not isinstance(history_proof, Mapping):
        history_proof = {}
    history_proof = copy.deepcopy(dict(history_proof))
    history_proof["proved"] = overall_passed
    history_proof["endpoint_composite"] = overall_passed
    history_proof["endpoint_composite_identity_proved"] = overall_passed
    history_proof["source_proof_composition"] = {
        "saved_rho_C2_equals_N_rho_C5": prior_ok,
        "same_rates_kernels_J_and_incoming_reused": prior_ok,
        "actual_memory_exponents_newly_bound": memory_ok,
        "proved": overall_passed,
        "conclusion": (
            "H_C2=N*H_C5 for ordinary Z0/1/2"
            if overall_passed else "unproved: source premises or memory binding failed"
        ),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "original_domain_conditional": prior_ok,
        "numeric_values_computed": False,
    }
    result["history_identity_proof"] = history_proof
    result["source_proof_composition"] = copy.deepcopy(history_proof["source_proof_composition"])

    result["source_inspection_audit"] = {
        "source_records_only": True,
        "static_source_AST_graph_and_saved_JSON_data_inspection_only": True,
        "nodes_hashed_only": True,
        "new_source_consumer_calls": 1,
        "actual_memory_pair_count": sum(len(row.get("checks", {})) for row in memories.values()),
        "prior_history_or_density_consumer_reentered": False,
        "accepted_ancestor_operator_calls": 0,
        "accepted_mathematical_operator_calls": 0,
        "old_mathematical_operator_calls": 0,
        "history_operator_calls": 0,
        "quadrature_or_norm_calls": 0,
        "numeric_field_evaluated": False,
        "old_operators_replayed": False,
        "numeric_values_computed": False,
        "full_goal_complete": False,
        "time_order_recursion_claimed": False,
        "timeorder_recursion_claimed": False,
    }

    source_pins = previous_report.get("source_pins")
    input_hashes = previous_report.get("input_hashes")
    implementation_hashes = previous_report.get("implementation_hashes")
    component_hashes = previous_report.get("component_implementation_hashes")
    checkpoint_hashes = previous_report.get("checkpoint_hashes")
    result["parent_provenance"] = {
        "raw_report": PREVIOUS_REPORT,
        "raw_report_sha256": PREVIOUS_REPORT_SHA256,
        "old_history_source": {"path": OLD_HISTORY_SOURCE, "sha256": OLD_HISTORY_SOURCE_SHA256},
        "old_history_runner": {"path": OLD_HISTORY_RUNNER, "sha256": OLD_HISTORY_RUNNER_SHA256},
        "loader": {"path": LOADER_SOURCE, "sha256": LOADER_SOURCE_SHA256},
        "parent_report": {"path": PARENT_REPORT, "sha256": PARENT_REPORT_SHA256},
        "accepted_source_stage_report": {"path": PROVIDER_REPORT, "sha256": PROVIDER_REPORT_SHA256},
        "source_pins": copy.deepcopy(dict(source_pins)) if isinstance(source_pins, Mapping) else {},
        "input_hashes": copy.deepcopy(dict(input_hashes)) if isinstance(input_hashes, Mapping) else {},
        "implementation_hashes": copy.deepcopy(dict(implementation_hashes)) if isinstance(implementation_hashes, Mapping) else {},
        "component_implementation_hashes": copy.deepcopy(dict(component_hashes)) if isinstance(component_hashes, Mapping) else {},
        "checkpoint_hashes": copy.deepcopy(dict(checkpoint_hashes)) if isinstance(checkpoint_hashes, Mapping) else {},
        "checkpoint_paths": list(CHECKPOINT_PATHS),
        "checkpoint_count": previous_report.get("checkpoint_count"),
        "validated": prior_checks.get("checkpoint_provenance") is True and prior_ok,
    }
    result["prior_report"] = {
        "path": PREVIOUS_REPORT,
        "sha256": PREVIOUS_REPORT_SHA256,
        "schema": previous_report.get("schema"),
        "component_schema": previous_component.get("schema") if isinstance(previous_component, Mapping) else None,
        "negative_gates": copy.deepcopy(previous_component.get("unresolved_gates", []))
        if isinstance(previous_component, Mapping) else [],
        "validation": prior_checks,
    }
    result["scope_limits"] = copy.deepcopy(result.get("scope_limits", {}))
    if isinstance(result["scope_limits"], dict):
        result["scope_limits"]["global_five_moment_closure"] = False
        result["scope_limits"]["full_goal_complete"] = False
        result["scope_limits"]["patch_source_theorem_proved"] = False
        result["scope_limits"]["history_identity_is_first_ten_only"] = True

    return result


__all__ = [
    "EXPECTED_NODE_COUNT",
    "EXPECTED_NODES_DIGEST",
    "RECORD_HASHES",
    "SCHEMA",
    "canonical_digest",
    "derive",
]
