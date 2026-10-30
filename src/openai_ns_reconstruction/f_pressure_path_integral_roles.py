"""Bounded source-role audit for the two saved F ``pressure.m15`` paths.

This module is deliberately self contained.  It does not import any previous
producer, graph class, checker, or runner.  The one-shot consumer authenticates
every manifest byte before decoding the saved prior receipt or the F archive.
After decoding it performs only the following finite reads from the F graph:

* integral rows ``279383`` and ``279479``;
* their integrand rows ``279382`` and ``279478``;
* the shared product row ``34721``; and
* at most one unique, non-protected sibling of ``34721``.

The radius target ``33982`` and all protected references are opaque terminals.
No graph traversal, BFS, source-graph evaluation, quadrature, or numerical
selection is performed.  Integer node references in the result are wrapped
with a namespace and are never treated as values.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import json
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any


SCHEMA = "f_pressure_path_integral_roles_v1"
COMPONENT = SCHEMA
PRIOR_COMPONENT = "f_pressure_radius_path_v1"
PRIOR_WRAPPER_SCHEMA = PRIOR_COMPONENT + "_result_v1"
PRIOR_INNER_SCHEMA = PRIOR_COMPONENT

F_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
F_ARCHIVE_SHA256 = "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
PRIOR_RESULT = Path("tmp/f_pressure_radius_path_v1/result.json")
PRIOR_RESULT_SHA256 = "82674aa39957d87cdbb067e9332a31bceea36cae85e397cf6e708c6dc1949847"

# This is the saved name used as ``pulse5.NAME`` by the producer.  Keeping the
# name literal avoids importing the producer module in a source-only consumer.
PULSE5_NAME = (
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions.json.gz"
)
PULSE5_NAMESPACE = "pulse5_C3_snapshot"
F_GRAPH_NAMESPACE = "F_exact_graph"

RADIUS_TARGET = 33982
SHARED_PRODUCT = 34721
DIRECT_NODE_IDS = frozenset({279383, 279479, 279382, 279478, SHARED_PRODUCT})
EXPECTED_ROOTS = {
    0: {"root": 279426, "integral": 279383, "integrand": 279382},
    1: {"root": 279542, "integral": 279479, "integrand": 279478},
}

PROTECTED_NODE_IDS = frozenset(
    {
        0,
        1,
        17,
        33572,
        25220,
        263344,
        263348,
        263349,
        262918,
        205147,
        263346,
    }
)

CHARTS = ("main", "gap", "end")
I1_CHANNELS = (
    "I1_q",
    "I1_q_Z",
    "I1_q_squared",
    "I1_U_q_Z",
    "I1_U_squared",
    "I1_U_U_Z",
)
ORDERS = (0, 1)

FORMAL_INTEGRANDS = {
    "I1_q": "Q",
    "I1_q_Z": "Q_Z",
    "I1_q_squared": "Q^2",
    "I1_U_q_Z": "U*Q_Z",
    "I1_U_squared": "U^2",
    "I1_U_U_Z": "U*U_Z",
}
FORMAL_COEFFICIENTS = {
    "I1_q": "-1/(2*L)",
    "I1_q_Z": "(1-delta)*Z/(4*L)",
    "I1_q_squared": "-1/4",
    "I1_U_q_Z": "(1-Z^2)/(2*L)",
    "I1_U_squared": "(1-delta)*Z^2/L^2",
    "I1_U_U_Z": "Z*(1-Z^2)/L^2",
}

_OPERAND_FIELDS = frozenset(
    {
        "integrand",
        "lower",
        "upper",
        "arguments",
        "argument",
        "numerator",
        "denominator",
        "left",
        "right",
        "body",
        "coordinate",
        "Z_variable",
        "free_query",
        "query_coordinate",
        "expression",
        "variable",
        "value",
    }
)


class FPressurePathIntegralRolesInputError(ValueError):
    """The frozen manifest, prior receipt, or selected rows are malformed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FPressurePathIntegralRolesInputError(message)


def _is_node_id(value: Any) -> bool:
    return type(value) is int and value >= 0


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canonical_json(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")


def _row_sha(row: Mapping[str, Any]) -> str:
    return _sha(_canonical_json(row))


def _opaque_ref(
    node: Any, namespace: str, role: str | None = None
) -> dict[str, Any]:
    _require(_is_node_id(node), "opaque node reference is malformed")
    result: dict[str, Any] = {
        "namespace": namespace,
        "node": node,
        "node_id_is_not_a_numeric_value": True,
    }
    if role is not None:
        result["source_role"] = role
    return result


def _path_name(value: Any) -> str:
    _require(isinstance(value, str) and value, "manifest path missing")
    return value.replace("\\", "/")


def _digest(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64, label + " SHA malformed")
    result = value.lower()
    _require(
        all(character in "0123456789abcdef" for character in result),
        label + " SHA malformed",
    )
    return result


def _safe(root: Path, name: str) -> Path:
    relative = Path(name)
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        "manifest path escapes root",
    )
    result = (root / relative).resolve()
    _require(result.is_relative_to(root.resolve()), "manifest path escapes root")
    return result


def _canonical_manifest_name(raw_name: str) -> str:
    """Accept only relative names or absolute aliases of the two inputs."""

    name = raw_name.replace("\\", "/")
    if name.startswith("/") or ":/" in name:
        aliases = (F_ARCHIVE.as_posix(), PRIOR_RESULT.as_posix())
        matches = [item for item in aliases if name.endswith("/" + item)]
        _require(len(matches) == 1, "absolute manifest path is not a required input")
        return matches[0]
    return name


def _manifest_hashes(
    source: Mapping[str, Any], hashes: Mapping[str, Any] | None
) -> dict[str, str]:
    entries: dict[str, str] = {}

    def add(name_value: Any, digest_value: Any, label: str) -> None:
        name = _path_name(name_value)
        digest = _digest(digest_value, label)
        _require(
            name not in entries or entries[name] == digest,
            "manifest hash conflict: " + name,
        )
        entries[name] = digest

    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), "manifest files entry malformed")
            add(item.get("path"), item.get("sha256"), f"manifest.files[{index}]")

    for field in ("sha256", "input_hashes", "hashes"):
        table = source.get(field)
        if isinstance(table, Mapping):
            for name, digest in table.items():
                add(name, digest, "manifest." + field)

    if isinstance(hashes, Mapping):
        for name, digest in hashes.items():
            add(name, digest, "caller hashes")

    _require(entries, "caller frozen manifest has no hashes")
    return entries


def _find_hash(entries: Mapping[str, str], name: str) -> str:
    if name in entries:
        return entries[name]
    normalized = name.replace("\\", "/")
    candidates = [
        value
        for key, value in entries.items()
        if key.replace("\\", "/").endswith("/" + normalized)
    ]
    _require(len(candidates) == 1, "required manifest hash missing: " + name)
    return candidates[0]


def _authenticate_all_raw(
    root: Path, source: Mapping[str, Any], hashes: Mapping[str, Any] | None
) -> tuple[dict[str, bytes], dict[str, str], int]:
    """Read and hash every manifest entry; do not decode anything yet."""

    entries = _manifest_hashes(source, hashes)
    required_expected = {
        F_ARCHIVE.as_posix(): _find_hash(entries, F_ARCHIVE.as_posix()),
        PRIOR_RESULT.as_posix(): _find_hash(entries, PRIOR_RESULT.as_posix()),
    }
    _require(
        required_expected[F_ARCHIVE.as_posix()] == F_ARCHIVE_SHA256,
        "required F archive SHA changed",
    )
    _require(
        required_expected[PRIOR_RESULT.as_posix()] == PRIOR_RESULT_SHA256,
        "required prior result SHA changed",
    )

    raw_by_path: dict[str, bytes] = {}
    actual_by_path: dict[str, str] = {}
    for raw_name, expected in entries.items():
        canonical = _canonical_manifest_name(raw_name)
        data = _safe(root, canonical).read_bytes()
        observed = _sha(data)
        _require(observed == expected, "raw SHA changed: " + canonical)
        if canonical in actual_by_path:
            _require(
                actual_by_path[canonical] == observed,
                "manifest aliases disagree: " + canonical,
            )
        raw_by_path[canonical] = data
        actual_by_path[canonical] = observed

    _require(F_ARCHIVE.as_posix() in raw_by_path, "authenticated F archive bytes missing")
    _require(PRIOR_RESULT.as_posix() in raw_by_path, "authenticated prior bytes missing")
    return raw_by_path, actual_by_path, len(entries)


def _decode_json(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise FPressurePathIntegralRolesInputError(label + " is not JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return value


def _decode_archive(raw: bytes) -> Mapping[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise FPressurePathIntegralRolesInputError("F archive is not valid gzip") from exc
    return _decode_json(decoded, "F archive")


def _unwrap_prior_result(prior: Mapping[str, Any]) -> Mapping[str, Any]:
    """Validate the saved prior-result envelope before using its inner leaf."""

    _require(
        prior.get("schema") == PRIOR_WRAPPER_SCHEMA,
        "prior result wrapper schema mismatch",
    )
    _require(
        prior.get("component") == PRIOR_COMPONENT,
        "prior result wrapper component mismatch",
    )
    _require(
        prior.get("status") == "returned",
        "prior result wrapper status is not returned",
    )
    _require(prior.get("error") is None, "prior result wrapper contains an error")

    inner = prior.get("f_pressure_radius_path")
    _require(isinstance(inner, Mapping), "prior result inner path is missing")
    _require(
        inner.get("schema") == PRIOR_INNER_SCHEMA,
        "prior result inner schema mismatch",
    )
    _require(
        inner.get("component") == PRIOR_COMPONENT,
        "prior result inner component mismatch",
    )
    _require(
        inner.get("status") == "bounded_F_pressure_radius_path_partial",
        "prior result inner status is outside the saved partial-hit scope",
    )
    _require(inner.get("existence") is True, "prior result inner has no structural hit")
    _require(inner.get("root_count") == 2, "prior result inner root count mismatch")
    _require(
        inner.get("pressure_m15_row_count") == 4,
        "prior result inner pressure row scope mismatch",
    )
    _require(inner.get("traversed_pressure_rows") == [0, 1], "prior result row scope mismatch")
    _require(inner.get("blocked_root_count") == 2, "prior result protected-boundary scope mismatch")
    _require(inner.get("blocked_roots") == [0, 1], "prior result blocked-root scope mismatch")

    scope = inner.get("scope")
    _require(isinstance(scope, Mapping), "prior result inner scope is missing")
    for name in (
        "source_values_computed",
        "source_values_selected",
        "numeric_values_computed",
        "numeric_use_admissible",
        "integrals_evaluated",
        "old_producer_replayed",
        "ancestor_replayed",
        "fresh_graph_generated",
        "prior_20_root_query_reentered",
        "P_F_G_compared",
        "target_row_dereferenced",
        "protected_root_indexed",
    ):
        _require(scope.get(name) is False, "prior result scope guard is not false: " + name)
    _require(
        scope.get("F_archive_exact_graph_traversed") is True,
        "prior result scope has no bounded F traversal marker",
    )

    archive_scope = inner.get("archive_scope")
    _require(isinstance(archive_scope, Mapping), "prior result archive scope is missing")
    for name in (
        "protected_nodes_indexed",
        "target_node_dereferenced",
        "prior_20_root_query_reentered",
    ):
        _require(
            archive_scope.get(name) is False,
            "prior result archive scope guard is not false: " + name,
        )
    _require(
        archive_scope.get("F_archive_exact_graph_traversed") is True,
        "prior result archive scope has no exact-graph marker",
    )
    _require(
        archive_scope.get("only_pressure_m15_Z0_Z1_roots_traversed") is True,
        "prior result archive root scope mismatch",
    )
    return inner


def _step_node(step: Any, label: str) -> int:
    _require(isinstance(step, Mapping), label + " path step malformed")
    node_ref = step.get("node_ref")
    _require(isinstance(node_ref, Mapping), label + " path node_ref missing")
    node = node_ref.get("node")
    _require(_is_node_id(node), label + " path node reference malformed")
    _require(
        node_ref.get("namespace") == F_GRAPH_NAMESPACE,
        label + " path namespace mismatch",
    )
    _require(
        node_ref.get("node_id_is_not_a_numeric_value") is True,
        label + " path opaque marker missing",
    )
    return node


def _validate_prior(prior: Mapping[str, Any]) -> dict[int, dict[str, Any]]:
    prior = _unwrap_prior_result(prior)
    _require(
        prior.get("schema") == PRIOR_INNER_SCHEMA,
        "prior result schema mismatch",
    )
    _require(prior.get("existence") is True, "prior result has no structural existence")
    _require(
        prior.get("traversed_pressure_rows") == [0, 1],
        "prior result does not have only pressure Z0/Z1 rows",
    )
    scope = prior.get("scope")
    _require(
        isinstance(scope, Mapping)
        and scope.get("prior_20_root_query_reentered") is False,
        "prior result reentered the retired 20-root query",
    )
    roots = prior.get("roots")
    _require(isinstance(roots, list) and len(roots) == 2, "prior result root scope mismatch")

    result: dict[int, dict[str, Any]] = {}
    for item in roots:
        _require(isinstance(item, Mapping), "prior root record malformed")
        order = item.get("row_order")
        _require(order in EXPECTED_ROOTS and order not in result, "prior root order mismatch")
        expected = EXPECTED_ROOTS[order]

        root_ref = item.get("root_ref")
        _require(isinstance(root_ref, Mapping), "prior root_ref missing")
        _require(root_ref.get("node") == expected["root"], "prior pressure root mismatch")
        _require(
            root_ref.get("namespace") == F_GRAPH_NAMESPACE
            and root_ref.get("node_id_is_not_a_numeric_value") is True,
            "prior pressure root opacity mismatch",
        )
        _require(
            root_ref.get("source_role")
            == "$.Rp_native_C5_pulse_source_bindings.pressure.m15[%d]" % order,
            "prior pressure root source role mismatch",
        )
        _require(item.get("existence") is True, "prior pressure root has no hit")
        _require(item.get("stopped_on_target") is True, "prior path did not stop on target")
        _require(item.get("target_row_dereferenced") is False, "prior target row was dereferenced")

        target_ref = item.get("target_ref")
        _require(isinstance(target_ref, Mapping), "prior target_ref missing")
        _require(target_ref.get("node") == RADIUS_TARGET, "prior radius target mismatch")
        _require(
            target_ref.get("namespace") == F_GRAPH_NAMESPACE
            and target_ref.get("node_id_is_not_a_numeric_value") is True,
            "prior target opacity mismatch",
        )

        path = item.get("path")
        _require(isinstance(path, list) and path, "prior witness path missing")
        ids = [_step_node(step, "prior") for step in path]
        _require(ids[0] == expected["root"], "prior witness path root mismatch")
        _require(ids[-1] == RADIUS_TARGET, "prior witness path target mismatch")
        _require(ids.count(RADIUS_TARGET) == 1, "prior witness has repeated radius target")
        _require(
            not any(node in PROTECTED_NODE_IDS for node in ids),
            "prior witness crosses a protected terminal",
        )

        try:
            integral_pos = ids.index(expected["integral"])
            integrand_pos = ids.index(expected["integrand"], integral_pos + 1)
            product_pos = ids.index(SHARED_PRODUCT, integrand_pos + 1)
            target_pos = ids.index(RADIUS_TARGET, product_pos + 1)
        except ValueError as exc:
            raise FPressurePathIntegralRolesInputError(
                "prior witness lacks the exact integral-to-radius chain"
            ) from exc
        _require(
            (integral_pos + 1 == integrand_pos)
            and (integrand_pos + 1 == product_pos)
            and (product_pos + 1 == target_pos),
            "prior witness chain is not the exact adjacent operand path",
        )
        integral_step = path[integral_pos]
        integrand_step = path[integrand_pos]
        product_step = path[product_pos]
        target_step = path[target_pos]
        _require(
            integral_step.get("operator") == "definite_integral",
            "prior integral operator mismatch",
        )
        _require(
            integrand_step.get("operand_label") == "integrand",
            "prior integral-to-integrand edge mismatch",
        )
        _require(
            integrand_step.get("operator") == "product",
            "prior integrand operator mismatch",
        )
        _require(
            product_step.get("operator") == "product",
            "prior shared-product operator mismatch",
        )
        _require(
            target_step.get("dereferenced") is False
            and target_step.get("operator") is None,
            "prior radius target is not terminal",
        )

        result[order] = {
            "row_order": order,
            "root_id": expected["root"],
            "integral_id": expected["integral"],
            "integrand_id": expected["integrand"],
            "shared_product_id": SHARED_PRODUCT,
            "radius_target_id": RADIUS_TARGET,
            "path_depth": len(path) - 1,
            "path_node_count": len(path),
            "path_chain_validated": True,
        }
    _require(set(result) == {0, 1}, "prior result does not contain both rows")
    return result


def _node_row(
    nodes: Any,
    node: int,
    read_ids: set[int],
    cache: dict[int, Mapping[str, Any]],
    *,
    allow_sibling: bool = False,
) -> Mapping[str, Any]:
    """Read one explicitly allowed row, never an operand discovered later."""

    _require(_is_node_id(node), "direct row reference malformed")
    _require(node != RADIUS_TARGET, "radius target row is protected from indexing")
    _require(node not in PROTECTED_NODE_IDS, "protected row is protected from indexing")
    _require(
        allow_sibling or node in DIRECT_NODE_IDS,
        "row is outside the fixed direct-index scope",
    )
    if node in cache:
        return cache[node]
    if isinstance(nodes, Mapping):
        _require(node in nodes, "direct graph row missing")
        row = nodes[node]
    else:
        _require(isinstance(nodes, list), "F exact graph nodes is not a list")
        _require(node < len(nodes), "direct graph row is out of bounds")
        row = nodes[node]
    _require(isinstance(row, Mapping), "direct graph row is malformed")
    read_ids.add(node)
    cache[node] = row
    return row


def _typed_edges(row: Mapping[str, Any], operation: Any) -> list[dict[str, Any]]:
    fields = {
        "definite_integral": ("integrand", "lower", "upper"),
        "product": ("arguments",),
    }.get(operation, ())
    edges: list[dict[str, Any]] = []
    for field in fields:
        if field not in row or row.get(field) is None:
            continue
        value = row.get(field)
        values: Iterable[Any]
        if isinstance(value, list):
            values = value
        else:
            values = (value,)
        for index, operand in enumerate(values):
            label = f"{field}[{index}]" if isinstance(value, list) else field
            if _is_node_id(operand):
                edges.append(
                    {
                        "field": label,
                        "ref": _opaque_ref(operand, F_GRAPH_NAMESPACE),
                        "typed_edge": True,
                    }
                )
            else:
                edges.append(
                    {
                        "field": label,
                        "typed_edge": False,
                        "malformed_operand_kind": type(operand).__name__,
                    }
                )
    return edges


def _metadata_fields(row: Mapping[str, Any]) -> dict[str, Any]:
    fields = (
        "measure",
        "variable",
        "endpoint_independent_of_Z",
        "physical_Pstar_not_reapplied",
        "ordinary_Z_derivative_order",
        "ordinary_slow_Z_derivative_order",
        "physical_radius_power",
        "physical_Jacobian_applied_once",
        "original_complete_pulse_chart",
        "original_pulse_own_history_retained",
        "auxiliary_pressure_q_equals_own_Q",
        "no_finite_history_or_chart_cutoff",
        "signed_integrand_not_its_magnitude_cap",
    )
    result: dict[str, Any] = {}
    for field in fields:
        if field not in row:
            continue
        value = row[field]
        if field == "variable" and _is_node_id(value):
            result[field] = _opaque_ref(value, F_GRAPH_NAMESPACE, "integral variable")
        else:
            result[field] = value
    return result


def _normalize_measure(value: Any) -> str:
    return "".join(str(value).split())


def _measure_candidate(chart: str, measure: Any) -> dict[str, Any] | None:
    normalized = _normalize_measure(measure)
    if chart in ("main", "gap") and normalized == "(R^2/mu)*dxi":
        return {
            "formal_measure": "R^2/mu dxi",
            "coordinate": "xi",
            "candidate_only": True,
            "source_measure_matches_chart": True,
        }
    if chart == "end" and normalized == "R^2*ds":
        return {
            "formal_measure": "R^2 ds",
            "coordinate": "s",
            "candidate_only": True,
            "source_measure_matches_chart": True,
        }
    return None


def _variable_matches_chart(chart: str, variable: Any) -> bool:
    if not isinstance(variable, str):
        return False
    normalized = variable.replace(" ", "")
    if chart == "main":
        return normalized.endswith("main_xi")
    if chart == "gap":
        return normalized.endswith("gap_xi")
    if chart == "end":
        return normalized.endswith("end_s")
    return False


def _integral_metadata(
    row: Mapping[str, Any], chart: str, channel: str, order: int, node: int
) -> dict[str, Any]:
    ordinary = row.get("ordinary_Z_derivative_order")
    if ordinary is None:
        ordinary = row.get("ordinary_slow_Z_derivative_order")
    declared_chart = row.get("original_complete_pulse_chart")
    if declared_chart is None:
        declared_chart = row.get("chart")
    measure = row.get("measure")
    candidate = _measure_candidate(chart, measure)
    reasons: list[str] = []
    if row.get("operation") != "definite_integral":
        reasons.append("operation_not_definite_integral")
    if declared_chart != chart:
        reasons.append("metadata_chart_mismatch")
    if measure is None:
        reasons.append("measure_missing")
    elif candidate is None:
        reasons.append("measure_chart_mismatch")
    if not _variable_matches_chart(chart, row.get("variable")):
        reasons.append("variable_chart_mismatch")
    if ordinary != order:
        reasons.append("ordinary_Z_order_mismatch")
    if row.get("endpoint_independent_of_Z") is not True:
        reasons.append("endpoint_independence_only_declaration_missing_or_false")
    if row.get("physical_Pstar_not_reapplied") is not True:
        reasons.append("physical_Pstar_not_reapplied_missing_or_false")

    lower = row.get("lower")
    upper = row.get("upper")
    result: dict[str, Any] = {
        "node_ref": _opaque_ref(node, F_GRAPH_NAMESPACE, f"F integral row {node}"),
        "operation": row.get("operation"),
        "chart": chart,
        "channel": channel,
        "formal_integrand": FORMAL_INTEGRANDS[channel],
        "measure": measure,
        "variable": row.get("variable"),
        "lower": _opaque_ref(lower, F_GRAPH_NAMESPACE, "integral lower"),
        "upper": _opaque_ref(upper, F_GRAPH_NAMESPACE, "integral upper"),
        "endpoint_independent_of_Z_declared": row.get("endpoint_independent_of_Z"),
        "physical_Pstar_not_reapplied_declared": row.get(
            "physical_Pstar_not_reapplied"
        ),
        "ordinary_Z_order": ordinary,
        "row_metadata": _metadata_fields(row),
        "typed_edges": _typed_edges(row, row.get("operation")),
        "row_sha256": _row_sha(row),
        "formal_measure_candidate": candidate,
        "metadata_consistent": not reasons,
        "unresolved_reasons": reasons,
        "actual_ordinary_Z_proved": False,
    }
    return result


def _integral_row_record(
    row: Mapping[str, Any], node: int, order: int
) -> dict[str, Any]:
    """Record the selected row without assigning it a chart or channel."""

    declared_chart = row.get("original_complete_pulse_chart")
    if declared_chart is None:
        declared_chart = row.get("chart")
    variable = row.get("variable")
    if _is_node_id(variable):
        variable_output: Any = _opaque_ref(variable, F_GRAPH_NAMESPACE, "integral variable")
    else:
        variable_output = variable
    lower = row.get("lower")
    upper = row.get("upper")
    ordinary = row.get("ordinary_Z_derivative_order")
    if ordinary is None:
        ordinary = row.get("ordinary_slow_Z_derivative_order")
    return {
        "row_order": order,
        "node_ref": _opaque_ref(node, F_GRAPH_NAMESPACE, "selected integral row"),
        "operation": row.get("operation"),
        "chart": declared_chart,
        "measure": row.get("measure"),
        "variable": variable_output,
        "lower": _opaque_ref(lower, F_GRAPH_NAMESPACE, "integral lower"),
        "upper": _opaque_ref(upper, F_GRAPH_NAMESPACE, "integral upper"),
        "endpoint_independent_of_Z_declared": row.get("endpoint_independent_of_Z"),
        "physical_Pstar_not_reapplied_declared": row.get(
            "physical_Pstar_not_reapplied"
        ),
        "ordinary_Z_order": ordinary,
        "row_metadata": _metadata_fields(row),
        "typed_edges": _typed_edges(row, row.get("operation")),
        "row_sha256": _row_sha(row),
    }


def _direct_integral_rows(
    nodes: Any,
    prior_info: Mapping[int, Mapping[str, Any]],
    read_ids: set[int],
    cache: dict[int, Mapping[str, Any]],
) -> tuple[dict[int, dict[str, Any]], dict[int, dict[str, Any]], dict[str, Any]]:
    integral_rows: dict[int, dict[str, Any]] = {}
    integrand_rows: dict[int, dict[str, Any]] = {}
    product_info: dict[str, Any] = {
        "node_ref": _opaque_ref(SHARED_PRODUCT, F_GRAPH_NAMESPACE, "shared product"),
        "integrand_edges": {},
        "row": None,
        "row_sha256": None,
        "typed_edges": [],
        "radius_target_edges": [],
        "sibling": None,
        "unresolved_reasons": [],
    }

    for order in ORDERS:
        expected = prior_info[order]
        integral_id = expected["integral_id"]
        integrand_id = expected["integrand_id"]
        integral_row = _node_row(nodes, integral_id, read_ids, cache)
        integrand_row = _node_row(nodes, integrand_id, read_ids, cache)
        integral_rows[order] = {
            "node_id": integral_id,
            "row": integral_row,
            "expected_integrand_id": integrand_id,
        }
        integrand_rows[order] = {
            "node_id": integrand_id,
            "row": integrand_row,
        }
        product_info["integrand_edges"][str(order)] = {
            "integral_ref": _opaque_ref(integral_id, F_GRAPH_NAMESPACE),
            "integrand_ref": _opaque_ref(integrand_id, F_GRAPH_NAMESPACE),
            "integral_integrand_edge_matches": integral_row.get("integrand")
            == integrand_id,
            "integral_typed_edges": _typed_edges(
                integral_row, integral_row.get("operation")
            ),
            "integrand_typed_edges": _typed_edges(
                integrand_row, integrand_row.get("operation")
            ),
            "integral_row_sha256": _row_sha(integral_row),
            "integrand_row_sha256": _row_sha(integrand_row),
            "unresolved_reasons": [],
        }
        if integral_row.get("integrand") != integrand_id:
            product_info["integrand_edges"][str(order)]["unresolved_reasons"].append(
                "integral_integrand_edge_mismatch"
            )
        if integrand_row.get("operation") != "product":
            product_info["integrand_edges"][str(order)]["unresolved_reasons"].append(
                "integrand_not_product"
            )
        arguments = integrand_row.get("arguments")
        if not isinstance(arguments, list) or arguments.count(SHARED_PRODUCT) != 1:
            product_info["integrand_edges"][str(order)]["unresolved_reasons"].append(
                "integrand_shared_product_edge_mismatch"
            )
        product_info["integrand_edges"][str(order)]["shared_product_edge"] = (
            arguments.count(SHARED_PRODUCT) == 1 if isinstance(arguments, list) else False
        )

    product_row = _node_row(nodes, SHARED_PRODUCT, read_ids, cache)
    product_info["row"] = {
        "operation": product_row.get("operation"),
        "row_metadata": _metadata_fields(product_row),
    }
    product_info["row_sha256"] = _row_sha(product_row)
    product_info["typed_edges"] = _typed_edges(product_row, product_row.get("operation"))
    arguments = product_row.get("arguments")
    if product_row.get("operation") != "product":
        product_info["unresolved_reasons"].append("shared_product_not_product")
        arguments = []
    if not isinstance(arguments, list):
        product_info["unresolved_reasons"].append("shared_product_arguments_missing")
        arguments = []
    radius_count = arguments.count(RADIUS_TARGET)
    product_info["radius_target_edges"] = [
        _opaque_ref(RADIUS_TARGET, F_GRAPH_NAMESPACE, "saved radius target")
        for item in arguments
        if item == RADIUS_TARGET
    ]
    if radius_count != 1:
        product_info["unresolved_reasons"].append("shared_product_radius_edge_mismatch")

    siblings = [item for item in arguments if item != RADIUS_TARGET]
    if len(siblings) == 1 and _is_node_id(siblings[0]):
        sibling = siblings[0]
        if sibling in PROTECTED_NODE_IDS:
            product_info["sibling"] = {
                "status": "protected_sibling_not_indexed",
                "ref": _opaque_ref(sibling, F_GRAPH_NAMESPACE, "protected sibling"),
            }
            product_info["unresolved_reasons"].append("protected_sibling")
        else:
            sibling_row = _node_row(
                nodes, sibling, read_ids, cache, allow_sibling=True
            )
            product_info["sibling"] = {
                "status": "indexed_for_jacobian_metadata_only",
                "ref": _opaque_ref(sibling, F_GRAPH_NAMESPACE, "jacobian sibling"),
                "operation": sibling_row.get("operation"),
                "row_metadata": _metadata_fields(sibling_row),
                "typed_edges": _typed_edges(sibling_row, sibling_row.get("operation")),
                "row_sha256": _row_sha(sibling_row),
                "jacobian_identity_proved": False,
            }
    elif len(siblings) == 0:
        product_info["unresolved_reasons"].append("jacobian_sibling_missing")
        product_info["sibling"] = {"status": "unresolved_no_sibling"}
    elif len(siblings) > 1:
        product_info["unresolved_reasons"].append("jacobian_sibling_ambiguous")
        product_info["sibling"] = {
            "status": "unresolved_multiple_siblings",
            "candidate_refs": [
                _opaque_ref(item, F_GRAPH_NAMESPACE, "unindexed sibling")
                if _is_node_id(item)
                else {"typed_edge": False, "kind": type(item).__name__}
                for item in siblings
            ],
        }
    else:
        product_info["unresolved_reasons"].append("jacobian_sibling_malformed")
        product_info["sibling"] = {"status": "unresolved_malformed_sibling"}

    return integral_rows, integrand_rows, product_info


def _mapping_records(
    bindings: Mapping[str, Any]
) -> tuple[list[dict[str, Any]], list[str]]:
    inputs = bindings.get("inputs")
    _require(isinstance(inputs, Mapping), "F source inputs mapping missing")
    pulse = inputs.get(PULSE5_NAME)
    _require(isinstance(pulse, Mapping), "saved pulse5 input snapshot missing")
    integrals = pulse.get("integrals")
    _require(isinstance(integrals, Mapping), "saved pulse5 integral snapshot missing")
    rewritten = bindings.get("rewritten_nodes")
    _require(isinstance(rewritten, Mapping), "saved rewritten_nodes map missing")

    records: list[dict[str, Any]] = []
    unresolved: list[str] = []
    for chart in CHARTS:
        chart_rows = integrals.get(chart)
        if not isinstance(chart_rows, Mapping):
            unresolved.append("missing_chart:" + chart)
            for channel in I1_CHANNELS:
                for order in ORDERS:
                    records.append(
                        {
                            "chart": chart,
                            "channel": channel,
                            "order": order,
                            "source_role": f"inputs[{PULSE5_NAME}].integrals[{chart}][{channel}][{order}]",
                            "status": "unresolved_missing_chart",
                        }
                    )
            continue
        for channel in I1_CHANNELS:
            rows = chart_rows.get(channel)
            if not isinstance(rows, list) or len(rows) < 2:
                unresolved.append(f"missing_rows:{chart}:{channel}")
                for order in ORDERS:
                    records.append(
                        {
                            "chart": chart,
                            "channel": channel,
                            "order": order,
                            "source_role": f"inputs[{PULSE5_NAME}].integrals[{chart}][{channel}][{order}]",
                            "status": "unresolved_missing_row",
                        }
                    )
                continue
            for order in ORDERS:
                source_id = rows[order]
                role = f"inputs[{PULSE5_NAME}].integrals[{chart}][{channel}][{order}]"
                record: dict[str, Any] = {
                    "chart": chart,
                    "channel": channel,
                    "order": order,
                    "source_role": role,
                    "formal_integrand": FORMAL_INTEGRANDS[channel],
                    "formal_m15_coefficient": FORMAL_COEFFICIENTS[channel],
                }
                if not _is_node_id(source_id):
                    record["status"] = "unresolved_malformed_source_ref"
                    record["source_ref_kind"] = type(source_id).__name__
                    unresolved.append("malformed_source_ref:" + role)
                    records.append(record)
                    continue
                record["source_ref"] = _opaque_ref(source_id, PULSE5_NAMESPACE, role)
                if source_id in PROTECTED_NODE_IDS:
                    record["status"] = "unresolved_protected_source_ref"
                    unresolved.append("protected_source_ref:" + role)
                    records.append(record)
                    continue
                key = str(source_id)
                if key not in rewritten:
                    record["status"] = "unresolved_missing_rewrite"
                    unresolved.append("missing_rewrite:" + role)
                    records.append(record)
                    continue
                target_id = rewritten[key]
                if not _is_node_id(target_id):
                    record["status"] = "unresolved_malformed_rewrite"
                    record["target_ref_kind"] = type(target_id).__name__
                    unresolved.append("malformed_rewrite:" + role)
                    records.append(record)
                    continue
                record["target_ref"] = _opaque_ref(target_id, F_GRAPH_NAMESPACE, role)
                record["status"] = (
                    "mapped_to_protected_target"
                    if target_id in PROTECTED_NODE_IDS
                    else "mapped"
                )
                if target_id in PROTECTED_NODE_IDS:
                    unresolved.append("protected_rewrite_target:" + role)
                records.append(record)
    _require(len(records) == 36, "bounded pulse5 lookup did not produce 36 rows")
    return records, unresolved


def _target_matches(records: list[dict[str, Any]], target: int) -> list[dict[str, Any]]:
    return [
        item
        for item in records
        if isinstance(item.get("target_ref"), Mapping)
        and item["target_ref"].get("node") == target
    ]


def _ordinary_z1_candidate(
    records: list[dict[str, Any]], chart: str, channel: str
) -> dict[str, Any]:
    pair = [
        item
        for item in records
        if item.get("chart") == chart and item.get("channel") == channel
    ]
    zeros = [item for item in pair if item.get("order") == 0 and item.get("status") == "mapped"]
    ones = [item for item in pair if item.get("order") == 1 and item.get("status") == "mapped"]
    candidate: dict[str, Any] = {
        "formula": "c*I_Z + c_Z*I",
        "coefficient_c": FORMAL_COEFFICIENTS[channel],
        "coefficient_c_Z": "unresolved",
        "base_integral_role": None,
        "integral_Z_role": None,
        "candidate_only": True,
        "actual_ordinary_Z_proved": False,
        "capture_free_substitution_proved": False,
        "sign_assumed": False,
        "unresolved_reasons": [],
    }
    if len(zeros) == 1:
        candidate["base_integral_role"] = zeros[0]["source_role"]
    else:
        candidate["unresolved_reasons"].append("base_integral_zero_or_multiple")
    if len(ones) == 1:
        candidate["integral_Z_role"] = ones[0]["source_role"]
    else:
        candidate["unresolved_reasons"].append("integral_Z_zero_or_multiple")
    return candidate


def _assemble_from_packets(
    prior: Mapping[str, Any], packet: Mapping[str, Any]
) -> dict[str, Any]:
    prior_info = _validate_prior(prior)
    bindings = packet.get("Rp_native_C5_pulse_source_bindings")
    _require(isinstance(bindings, Mapping), "F source bindings missing")
    nodes = packet.get("exact_graph_nodes")
    _require(isinstance(nodes, (list, Mapping)), "F exact graph nodes missing")

    records, map_unresolved = _mapping_records(bindings)
    read_ids: set[int] = set()
    cache: dict[int, Mapping[str, Any]] = {}
    integral_rows, integrand_rows, product_info = _direct_integral_rows(
        nodes, prior_info, read_ids, cache
    )

    bindings_out: list[dict[str, Any]] = []
    semantic_unresolved = list(map_unresolved)
    for order in ORDERS:
        expected = prior_info[order]
        integral_id = expected["integral_id"]
        matches = _target_matches(records, integral_id)
        binding: dict[str, Any] = {
            "root_order": order,
            "root_ref": _opaque_ref(expected["root_id"], F_GRAPH_NAMESPACE),
            "integral_ref": _opaque_ref(integral_id, F_GRAPH_NAMESPACE),
            "integrand_ref": _opaque_ref(expected["integrand_id"], F_GRAPH_NAMESPACE),
            "candidate_source_roles": [
                {
                    "chart": item.get("chart"),
                    "channel": item.get("channel"),
                    "order": item.get("order"),
                    "source_role": item.get("source_role"),
                    "source_ref": item.get("source_ref"),
                    "target_ref": item.get("target_ref"),
                    "status": item.get("status"),
                }
                for item in matches
            ],
            "unique_chart_channel": len(matches) == 1,
            "unresolved_reasons": [],
        }
        if len(matches) != 1:
            binding["unresolved_reasons"].append(
                "integral_target_zero_or_multiple_source_roles"
            )
            semantic_unresolved.append(
                f"integral_target_zero_or_multiple:{order}:{integral_id}"
            )
        else:
            match = matches[0]
            chart = match["chart"]
            channel = match["channel"]
            binding["chart"] = chart
            binding["channel"] = channel
            binding["source_role"] = match["source_role"]
            binding["formal_integrand"] = FORMAL_INTEGRANDS[channel]
            binding["formal_m15_coefficient"] = FORMAL_COEFFICIENTS[channel]
            binding["integral_metadata"] = _integral_metadata(
                integral_rows[order]["row"], chart, channel, match["order"], integral_id
            )
            if match["order"] != order:
                binding["unresolved_reasons"].append("source_order_target_mismatch")
            if binding["integral_metadata"]["unresolved_reasons"]:
                binding["unresolved_reasons"].extend(
                    binding["integral_metadata"]["unresolved_reasons"]
                )
            if binding["unresolved_reasons"]:
                semantic_unresolved.extend(
                    f"integral_metadata:{order}:{reason}"
                    for reason in binding["unresolved_reasons"]
                )
        bindings_out.append(binding)

    common: dict[str, Any] = {
        "same_chart_channel": False,
        "chart_channel_comparison": {
            "performed": False,
            "Z0": None,
            "Z1": None,
            "same": False,
            "comparison_only": True,
        },
        "metadata_comparison": {
            "performed": False,
            "comparison_only": True,
            "measure_field_equal": None,
            "variable_field_equal": None,
            "endpoint_declaration_equal": None,
            "physical_Pstar_declaration_equal": None,
            "measure_owner_proved": False,
            "variable_owner_proved": False,
            "endpoint_owner_proved": False,
            "common_measure_proved": False,
            "unresolved_reasons": [],
        },
        "ordinary_Z1_candidate": None,
        "unresolved_reasons": [],
    }
    if len(bindings_out) == 2 and all(item.get("unique_chart_channel") for item in bindings_out):
        left, right = bindings_out
        same = left.get("chart") == right.get("chart") and left.get("channel") == right.get("channel")
        common["same_chart_channel"] = same
        common["chart_channel_comparison"] = {
            "performed": True,
            "Z0": {"chart": left.get("chart"), "channel": left.get("channel")},
            "Z1": {"chart": right.get("chart"), "channel": right.get("channel")},
            "same": same,
            "comparison_only": True,
        }
        left_metadata = left.get("integral_metadata")
        right_metadata = right.get("integral_metadata")
        if isinstance(left_metadata, Mapping) and isinstance(right_metadata, Mapping):
            common["metadata_comparison"] = {
                "performed": True,
                "comparison_only": True,
                "measure_field_equal": left_metadata.get("measure")
                == right_metadata.get("measure"),
                "variable_field_equal": left_metadata.get("variable")
                == right_metadata.get("variable"),
                "endpoint_declaration_equal": left_metadata.get(
                    "endpoint_independent_of_Z_declared"
                )
                == right_metadata.get("endpoint_independent_of_Z_declared"),
                "physical_Pstar_declaration_equal": left_metadata.get(
                    "physical_Pstar_not_reapplied_declared"
                )
                == right_metadata.get("physical_Pstar_not_reapplied_declared"),
                "measure_owner_proved": False,
                "variable_owner_proved": False,
                "endpoint_owner_proved": False,
                "common_measure_proved": False,
                "unresolved_reasons": [],
            }
        else:
            common["metadata_comparison"]["unresolved_reasons"].append(
                "integral_metadata_missing_for_comparison"
            )
        if same:
            common["ordinary_Z1_candidate"] = _ordinary_z1_candidate(
                records, left["chart"], left["channel"]
            )
            if common["ordinary_Z1_candidate"]["unresolved_reasons"]:
                common["unresolved_reasons"].extend(
                    common["ordinary_Z1_candidate"]["unresolved_reasons"]
                )
        else:
            common["unresolved_reasons"].append("Z0_Z1_chart_channel_mismatch")
    else:
        common["unresolved_reasons"].append("Z0_Z1_nonunique_source_roles")
        common["metadata_comparison"]["unresolved_reasons"].append(
            "metadata_comparison_requires_two_unique_roles"
        )
    semantic_unresolved.extend(
        "common_Z0_Z1:" + reason for reason in common["unresolved_reasons"]
    )
    semantic_unresolved.extend(
        "product:" + reason for reason in product_info["unresolved_reasons"]
    )

    direct_rows = [
        _opaque_ref(node, F_GRAPH_NAMESPACE, "directly indexed row")
        for node in sorted(read_ids)
    ]
    allowed_read_ids = set(DIRECT_NODE_IDS)
    sibling_record = product_info.get("sibling")
    if isinstance(sibling_record, Mapping):
        sibling_ref = sibling_record.get("ref")
        if isinstance(sibling_ref, Mapping) and _is_node_id(sibling_ref.get("node")):
            allowed_read_ids.add(sibling_ref["node"])
    _require(read_ids.issubset(allowed_read_ids), "direct row scope widened")
    status = (
        "bounded_F_integral_roles_found"
        if not semantic_unresolved
        else "bounded_F_integral_roles_partial"
    )
    return {
        "schema": SCHEMA,
        "component": COMPONENT,
        "status": status,
        "source_only": True,
        "archive_whole_json_decoded": True,
        "bounded_graph_direct_reads": True,
        "unresolved_matches_preserved": True,
        "source_values_computed": False,
        "source_values_selected": False,
        "numeric_values_computed": False,
        "numeric_use_admissible": False,
        "integrals_evaluated": False,
        "quadrature_evaluated": False,
        "source_graph_traversed": False,
        "old_source_graph_traversed": False,
        "actual_source_values": False,
        "graph_traversal": False,
        "protected_nodes_indexed": False,
        "target_node_dereferenced": False,
        "prior_scope": {
            "only_pressure_m15_Z0_Z1": True,
            "roots": [
                {
                    "row_order": prior_info[order]["row_order"],
                    "root_ref": _opaque_ref(prior_info[order]["root_id"], F_GRAPH_NAMESPACE),
                    "integral_ref": _opaque_ref(prior_info[order]["integral_id"], F_GRAPH_NAMESPACE),
                    "integrand_ref": _opaque_ref(prior_info[order]["integrand_id"], F_GRAPH_NAMESPACE),
                    "shared_product_ref": _opaque_ref(SHARED_PRODUCT, F_GRAPH_NAMESPACE),
                    "radius_target_ref": _opaque_ref(RADIUS_TARGET, F_GRAPH_NAMESPACE),
                    "path_chain_validated": True,
                    "path_depth": prior_info[order]["path_depth"],
                }
                for order in ORDERS
            ],
            "prior_20_root_query_reentered": False,
        },
        "bounded_lookup": {
            "pulse5_name": PULSE5_NAME,
            "charts": list(CHARTS),
            "channels": list(I1_CHANNELS),
            "ordinary_Z_orders": list(ORDERS),
            "source_ref_count": len(records),
            "source_refs": records,
            "rewritten_nodes_direct_lookup_only": True,
            "source_graph_nodes_read": False,
            "source_ids_are_not_values": True,
        },
        "integral_role_mapping_count": len(records),
        "integral_role_mappings": records,
        "integral_bindings": bindings_out,
        "integral_rows": {
            str(order): _integral_row_record(
                integral_rows[order]["row"], integral_rows[order]["node_id"], order
            )
            for order in ORDERS
        },
        "integrand_rows": {
            str(order): {
                "node_ref": _opaque_ref(integrand_rows[order]["node_id"], F_GRAPH_NAMESPACE),
                "operation": integrand_rows[order]["row"].get("operation"),
                "typed_edges": _typed_edges(
                    integrand_rows[order]["row"], integrand_rows[order]["row"].get("operation")
                ),
                "row_sha256": _row_sha(integrand_rows[order]["row"]),
            }
            for order in ORDERS
        },
        "shared_product": product_info,
        "common_Z0_Z1_source": common,
        "direct_index_scope": {
            "fixed_rows": [
                _opaque_ref(node, F_GRAPH_NAMESPACE, "fixed direct row")
                for node in sorted(DIRECT_NODE_IDS)
            ],
            "indexed_rows": direct_rows,
            "radius_target_indexed": False,
            "protected_rows_indexed": False,
            "at_most_one_sibling_indexed": True,
            "graph_edges_followed": False,
            "BFS_rerun": False,
            "ancestor_replayed": False,
        },
        "direct_graph_row_ids": sorted(read_ids),
        "source_gates": {
            "capture_free_substitution_proved": False,
            "source_identity_proved": False,
            "actual_ordinary_Z_transport_proved": False,
            "endpoint_transport_proved": False,
            "common_measure_proved": False,
            "error_transport_proved": False,
            "m15_closure_proved": False,
            "numeric_use_admissible": False,
        },
        "physical_limits": {
            "R2_over_mu_vs_R2_selected_from_metadata_only": True,
            "measure_from_interval_estimate": False,
            "unknown_factors_assigned_positive_sign": False,
            "source_values_computed": False,
            "integrals_evaluated": False,
        },
        "semantic_unresolved": sorted(set(semantic_unresolved)),
    }


def build_integral_roles(
    root: str | Path = ".",
    source: Mapping[str, Any] | None = None,
    hashes: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Authenticate and perform the fixed five-row source-role audit once."""

    _require(isinstance(source, Mapping), "caller frozen manifest is required")
    repo = Path(root).resolve()
    raw, actual, manifest_count = _authenticate_all_raw(repo, source, hashes)

    # Decoding begins only after every manifest raw SHA has passed.
    prior = _decode_json(raw[PRIOR_RESULT.as_posix()], "prior result")
    archive = _decode_archive(raw[F_ARCHIVE.as_posix()])
    result = _assemble_from_packets(prior, archive)
    result["raw_hash_authentication"] = {
        "all_manifest_raw_hashes_checked_before_decode": True,
        "manifest_entry_count": manifest_count,
        "prior_result": {
            "path": PRIOR_RESULT.as_posix(),
            "expected_sha256": PRIOR_RESULT_SHA256,
            "actual_sha256": actual[PRIOR_RESULT.as_posix()],
            "decoded_after_raw_hash": True,
        },
        "F_archive": {
            "path": F_ARCHIVE.as_posix(),
            "expected_sha256": F_ARCHIVE_SHA256,
            "actual_sha256": actual[F_ARCHIVE.as_posix()],
            "decompressed_once_after_raw_hash": True,
        },
    }
    return result


def _synthetic_prior() -> dict[str, Any]:
    def step(node: int, operation: str | None, via: str | None, dereferenced: bool) -> dict[str, Any]:
        return {
            "node_ref": _opaque_ref(node, F_GRAPH_NAMESPACE),
            "operator": operation,
            "operand_label": via,
            "depth": 0,
            "dereferenced": dereferenced,
        }

    roots = []
    for order, values in EXPECTED_ROOTS.items():
        path = [
            step(values["root"], "sum", "<root>", True),
            step(values["integral"], "definite_integral", "arguments[0]", True),
            step(values["integrand"], "product", "integrand", True),
            step(SHARED_PRODUCT, "product", "arguments[0]", True),
            step(RADIUS_TARGET, None, "arguments[0]", False),
        ]
        for depth, item in enumerate(path):
            item["depth"] = depth
        roots.append(
            {
                "row_order": order,
                "root_ref": _opaque_ref(
                    values["root"],
                    F_GRAPH_NAMESPACE,
                    "$.Rp_native_C5_pulse_source_bindings.pressure.m15[%d]" % order,
                ),
                "target_ref": _opaque_ref(RADIUS_TARGET, F_GRAPH_NAMESPACE),
                "existence": True,
                "stopped_on_target": True,
                "target_row_dereferenced": False,
                "path": path,
            }
        )
    inner = {
        "schema": PRIOR_INNER_SCHEMA,
        "component": PRIOR_COMPONENT,
        "status": "bounded_F_pressure_radius_path_partial",
        "existence": True,
        "root_count": 2,
        "pressure_m15_row_count": 4,
        "traversed_pressure_rows": [0, 1],
        "blocked_root_count": 2,
        "blocked_roots": [0, 1],
        "roots": roots,
        "archive_scope": {
            "F_archive_exact_graph_traversed": True,
            "protected_nodes_indexed": False,
            "target_node_dereferenced": False,
            "only_pressure_m15_Z0_Z1_roots_traversed": True,
            "prior_20_root_query_reentered": False,
        },
        "scope": {
            "F_archive_exact_graph_traversed": True,
            "source_values_computed": False,
            "source_values_selected": False,
            "numeric_values_computed": False,
            "numeric_use_admissible": False,
            "integrals_evaluated": False,
            "old_producer_replayed": False,
            "ancestor_replayed": False,
            "fresh_graph_generated": False,
            "prior_20_root_query_reentered": False,
            "P_F_G_compared": False,
            "target_row_dereferenced": False,
            "protected_root_indexed": False,
        },
    }
    return {
        "schema": PRIOR_WRAPPER_SCHEMA,
        "component": PRIOR_COMPONENT,
        "phase": PRIOR_COMPONENT + "_pending",
        "status": "returned",
        "error": None,
        "f_pressure_radius_path": inner,
    }


def _synthetic_packet(
    *,
    metadata_mismatch: bool = False,
    zero_mapping: bool = False,
    multiple_mapping: bool = False,
    protected_sibling: bool = False,
    moving_endpoint: bool = False,
) -> dict[str, Any]:
    nodes: dict[int, dict[str, Any]] = {
        279383: {
            "operation": "definite_integral",
            "integrand": 279382,
            "lower": 5,
            "upper": 6,
            "variable": "original_own_pulse_main_xi",
            "measure": "(R^2/mu)*dxi",
            "original_complete_pulse_chart": "gap" if metadata_mismatch else "main",
            "endpoint_independent_of_Z": not moving_endpoint,
            "physical_Pstar_not_reapplied": True,
            "ordinary_slow_Z_derivative_order": 0,
            "physical_radius_power": 2,
            "physical_Jacobian_applied_once": True,
        },
        279479: {
            "operation": "definite_integral",
            "integrand": 279478,
            "lower": 7,
            "upper": 8,
            "variable": "original_own_pulse_main_xi",
            "measure": "(R^2/mu)*dxi",
            "original_complete_pulse_chart": "main",
            "endpoint_independent_of_Z": True,
            "physical_Pstar_not_reapplied": True,
            "ordinary_slow_Z_derivative_order": 1,
            "physical_radius_power": 2,
            "physical_Jacobian_applied_once": True,
        },
        279382: {"operation": "product", "arguments": [34721, 500]},
        279478: {"operation": "product", "arguments": [34721, 501]},
        34721: {
            "operation": "product",
            "arguments": [RADIUS_TARGET, 400 if not protected_sibling else 17],
        },
        400: {"operation": "exact_rational"},
        5: {"operation": "exact_rational"},
        6: {"operation": "exact_rational"},
        7: {"operation": "exact_rational"},
        8: {"operation": "exact_rational"},
    }

    integral_data: dict[str, dict[str, list[int]]] = {}
    rewritten: dict[str, int] = {}
    next_source = 100
    next_target = 50000
    for chart in CHARTS:
        integral_data[chart] = {}
        for channel in I1_CHANNELS:
            rows: list[int] = []
            for order in range(4):
                source_id = next_source
                next_source += 1
                rows.append(source_id)
                if chart == "main" and channel == "I1_q" and order == 0:
                    rewritten[str(source_id)] = 279383
                elif chart == "main" and channel == "I1_q" and order == 1:
                    rewritten[str(source_id)] = 279479
                else:
                    rewritten[str(source_id)] = next_target
                    next_target += 1
            integral_data[chart][channel] = rows

    if zero_mapping:
        source_id = integral_data["main"]["I1_q"][1]
        rewritten.pop(str(source_id), None)
    if multiple_mapping:
        source_id = integral_data["gap"]["I1_q"][0]
        rewritten[str(source_id)] = 279383

    return {
        "exact_graph_nodes": nodes,
        "Rp_native_C5_pulse_source_bindings": {
            "inputs": {PULSE5_NAME: {"integrals": integral_data}},
            "rewritten_nodes": rewritten,
        },
    }


def verify_synthetic_fixture() -> dict[str, Any]:
    """Exercise metadata-only guards without filesystem or graph imports."""

    baseline = _assemble_from_packets(_synthetic_prior(), _synthetic_packet())
    metadata_mismatch = _assemble_from_packets(
        _synthetic_prior(), _synthetic_packet(metadata_mismatch=True)
    )
    zero_mapping = _assemble_from_packets(
        _synthetic_prior(), _synthetic_packet(zero_mapping=True)
    )
    multiple_mapping = _assemble_from_packets(
        _synthetic_prior(), _synthetic_packet(multiple_mapping=True)
    )
    protected_sibling = _assemble_from_packets(
        _synthetic_prior(), _synthetic_packet(protected_sibling=True)
    )
    moving_endpoint = _assemble_from_packets(
        _synthetic_prior(), _synthetic_packet(moving_endpoint=True)
    )

    candidate = baseline["common_Z0_Z1_source"]["ordinary_Z1_candidate"]
    checks = {
        "baseline_unique_same_chart_channel": baseline["common_Z0_Z1_source"][
            "same_chart_channel"
        ],
        "metadata_mismatch_is_unresolved": any(
            "metadata_chart_mismatch" in item
            for item in metadata_mismatch["semantic_unresolved"]
        ),
        "zero_mapping_is_unresolved": any(
            "missing_rewrite" in item for item in zero_mapping["semantic_unresolved"]
        ),
        "multiple_mapping_is_unresolved": any(
            "zero_or_multiple" in item
            for item in multiple_mapping["semantic_unresolved"]
        ),
        "protected_sibling_not_indexed": protected_sibling["shared_product"]["sibling"][
            "status"
        ]
        == "protected_sibling_not_indexed",
        "moving_endpoint_is_unresolved": any(
            "endpoint_independence" in item
            for item in moving_endpoint["semantic_unresolved"]
        ),
        "complete_c_z_guard": (
            isinstance(candidate, Mapping)
            and candidate.get("formula") == "c*I_Z + c_Z*I"
            and candidate.get("coefficient_c_Z") == "unresolved"
            and candidate.get("candidate_only") is True
            and candidate.get("actual_ordinary_Z_proved") is False
            and candidate.get("sign_assumed") is False
        ),
        "radius_and_protected_rows_not_indexed": (
            baseline["direct_index_scope"]["radius_target_indexed"] is False
            and baseline["direct_index_scope"]["protected_rows_indexed"] is False
        ),
    }
    result = {
        "schema": SCHEMA + "_synthetic_fixture",
        "passed": all(checks.values()),
        "source_values_computed": False,
        "source_values_selected": False,
        "actual_source_values": False,
        "numeric_values_computed": False,
        "numeric_use_admissible": False,
        "integrals_evaluated": False,
        "quadrature_evaluated": False,
        "source_graph_traversed": False,
        "old_source_graph_traversed": False,
        "graph_traversal": False,
        "protected_nodes_indexed": False,
        "target_node_dereferenced": False,
        "source_identity_proved": False,
        "capture_free_substitution_proved": False,
        "ordinary_Z_transport_proved": False,
        "endpoint_transport_proved": False,
        "common_measure_proved": False,
        "error_transport_proved": False,
        "m15_closure_proved": False,
        "checks": checks,
        "json_safe": False,
    }
    try:
        json.dumps(result, sort_keys=True, separators=(",", ":"))
        result["json_safe"] = True
    except (TypeError, ValueError):
        result["json_safe"] = False
    result["passed"] = result["passed"] and result["json_safe"]
    return result


__all__ = [
    "FPressurePathIntegralRolesInputError",
    "build_integral_roles",
    "verify_synthetic_fixture",
]
