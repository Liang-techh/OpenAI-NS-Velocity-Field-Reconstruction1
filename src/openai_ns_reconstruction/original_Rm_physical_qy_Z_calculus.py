"""Source provenance for the ordinary ``Z`` derivative of physical ``Q_y``.

This consumer is intentionally narrower than the existing ``Q_y[Z0]`` value
seam.  It does not create a derivative graph, evaluate a field, differentiate
an alias, or infer a paired ``Z2`` value.  It checks the saved arithmetic DAG
and uses only the six actual source-jet edges admitted by the primitive
receipt.  A receipt failure leaves the arithmetic result conditional and
cannot produce a physical ``Q_y[Z1]`` seam.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
import hashlib
import json
from typing import Any


SCHEMA = "original_Rm_physical_qy_Z_derivative_provenance_v1"
PRIMITIVE_SCHEMA = "original_Rm_physical_primitive_Z_jets_v1"
QY0_SCHEMA = "original_Rm_physical_qy0_value_seam_v1"
COMMON_PHASE_SCHEMA = "actual_Rm_common_phase_source_v1"
SIDES = ("postrestore", "patch")
Z_NODE = 20
PSTAR_NODE = 33572
COMMON_DELTA_NODE = 25220
COMMON_DELTA_ONE_NODE = 25222
COMMON_D_NODE = 25251
COMMON_D1_NODE = 25247
COMMON_L_NODE = 25272
COMMON_L1_NODE = 25275
GRAPH_NODE_COUNT = 205300
DERIVATIVE_PREFIX = 205147
EXPECTED_GRAPH_DIGEST = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
EXPECTED_FINITE_PREFIX_DIGEST = "c7ab7da8cc5c51cd782a5794e3076a6cf021a45235aae1c75bd246b01491b639"
EXPECTED_QY0_PACKET_DIGEST = "5bad9b551ebe3d48e382fa2602c2967ed7ca51661911692396d693aa24c8b0d7"
EXPECTED_COMMON_PHASE_PACKET_DIGEST = "05d77b68c94251784fbf919205ef59404e1b6124ee4d7599a89485c13155b769"
POST_COORDINATE = 84528
PATCH_COORDINATE = 185676
EXPECTED_ENDPOINT = 1
PRODUCER_MODULE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py"
)
PRODUCER_MODULE_SHA256 = "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7"


EXPECTED_CONTEXTS = {
    "postrestore": {"coordinate": POST_COORDINATE, "endpoint": EXPECTED_ENDPOINT},
    "patch": {"coordinate": PATCH_COORDINATE, "endpoint": EXPECTED_ENDPOINT},
}

EXPECTED_SOURCE_EDGES: dict[str, dict[str, dict[int, int]]] = {
    "postrestore": {
        "U_y": {0: 183274, 1: 183275},
        "physical_mean_y": {0: 183292, 1: 183293, 2: 183294},
    },
    "patch": {
        "U_y": {0: 204323, 1: 204324},
        "physical_mean_y": {0: 204341, 1: 204342, 2: 204343},
    },
}
EXPECTED_SOURCE_DERIVATIVES = {
    (side, field, order): EXPECTED_SOURCE_EDGES[side][field][order + 1]
    for side in SIDES
    for field, orders in (("U_y", (0,)), ("physical_mean_y", (0, 1)))
    for order in orders
}

EXPECTED_QY0 = {"postrestore": 183451, "patch": 204500}
EXPECTED_QY1 = {"postrestore": 183455, "patch": 204504}


class PhysicalQyDerivativeProvenanceError(ValueError):
    """Malformed graph input; ordinary unresolved gates are returned."""


def canonical_digest(value: Any) -> str:
    """Return the canonical digest used by the saved source graph packets."""

    raw = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _digest(value: Any, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise PhysicalQyDerivativeProvenanceError(f"{label} must be a SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise PhysicalQyDerivativeProvenanceError(f"{label} must be a SHA-256 digest") from exc
    return value.lower()


def _contains_true_flag(value: Any, key: str) -> bool:
    """Find a contradictory true flag at any saved-packet nesting level."""

    if isinstance(value, Mapping):
        if value.get(key) is True:
            return True
        return any(_contains_true_flag(child, key) for child in value.values())
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return any(_contains_true_flag(child, key) for child in value)
    return False


def _mapping(value: Any, label: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise PhysicalQyDerivativeProvenanceError(f"{label} must be a mapping")
    return value


def _row(nodes: Sequence[Any], index: int, label: str = "node") -> Mapping[str, Any]:
    if not isinstance(index, int) or isinstance(index, bool) or index < 0 or index >= len(nodes):
        raise PhysicalQyDerivativeProvenanceError(f"{label} {index!r} is outside the graph")
    value = nodes[index]
    if not isinstance(value, Mapping):
        raise PhysicalQyDerivativeProvenanceError(f"{label} {index} is not a mapping")
    return value


def _operation(row: Mapping[str, Any], index: int) -> str:
    value = row.get("operation", row.get("op"))
    if not isinstance(value, str) or not value:
        raise PhysicalQyDerivativeProvenanceError(f"node {index} has no operation")
    return value


def _expected_row(nodes: Sequence[Any], index: int, expected: Mapping[str, Any], errors: list[dict[str, Any]], label: str) -> None:
    try:
        actual = _row(nodes, index)
    except PhysicalQyDerivativeProvenanceError as exc:
        errors.append({"gate": label, "node": index, "error": str(exc)})
        return
    if dict(actual) != dict(expected):
        errors.append({"gate": label, "node": index, "expected": dict(expected), "actual": dict(actual)})


def _source_rows() -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {
        Z_NODE: {"operation": "bound_variable", "name": "Z"},
        25218: {"operation": "sum", "arguments": [14, 4194]},
        25219: {"operation": "negative", "argument": 25218},
        COMMON_DELTA_NODE: {"operation": "analytic_unary", "name": "exp", "argument": 25219},
        25221: {"operation": "negative", "argument": COMMON_DELTA_NODE},
        COMMON_DELTA_ONE_NODE: {"operation": "sum", "arguments": [1, 25221]},
        25246: {"operation": "sum", "arguments": [Z_NODE, Z_NODE]},
        COMMON_D1_NODE: {"operation": "negative", "argument": 25246},
        25249: {"operation": "product", "arguments": [Z_NODE, Z_NODE]},
        25250: {"operation": "negative", "argument": 25249},
        COMMON_D_NODE: {"operation": "sum", "arguments": [1, 25250]},
        25270: {"operation": "product", "arguments": [COMMON_DELTA_NODE, 25249]},
        25271: {"operation": "negative", "argument": 25270},
        COMMON_L_NODE: {"operation": "sum", "arguments": [1, 25271]},
        25274: {"operation": "product", "arguments": [COMMON_DELTA_NODE, 25246]},
        COMMON_L1_NODE: {"operation": "negative", "argument": 25274},
        PSTAR_NODE: {"operation": "analytic_unary", "name": "exp", "argument": 5},
    }
    return rows


def _physical_rows() -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {}
    charts = {
        "postrestore": {
            "U": ([183268, 183269, 183270], [121571, 121572, 121352]),
            "U_y": ([183274, 183275, 183276], [182757, 182758, 182759]),
            "m": ([183280, 183281, 183282], [121720, 121721, 121722]),
            "neg_m": ([183286, 183287, 183288], [183280, 183281, 183282]),
            "mean_y": ([183292, 183293, 183294], [183268, 183269, 183270]),
        },
        "patch": {
            "U": ([204317, 204318, 204319], [198423, 198424, 198425]),
            "U_y": ([204323, 204324, 204325], [203806, 203807, 203808]),
            "m": ([204329, 204330, 204331], [198591, 198592, 198593]),
            "neg_m": ([204335, 204336, 204337], [204329, 204330, 204331]),
            "mean_y": ([204341, 204342, 204343], [204317, 204318, 204319]),
        },
    }
    for chart in charts.values():
        for target, source in zip(*chart["U"]):
            rows[target] = {"operation": "product", "arguments": [source, PSTAR_NODE]}
        for target, source in zip(*chart["U_y"]):
            rows[target] = {"operation": "product", "arguments": [source, PSTAR_NODE]}
        for target, source in zip(*chart["m"]):
            rows[target] = {"operation": "product", "arguments": [source, PSTAR_NODE]}
        for target, source in zip(*chart["neg_m"]):
            rows[target] = {"operation": "negative", "argument": source}
        for target, source in zip(*chart["mean_y"]):
            rows[target] = {"operation": "sum", "arguments": [source, target - 6]}
    return rows


def _qy_rows() -> dict[int, dict[str, Any]]:
    rows: dict[int, dict[str, Any]] = {
        # The source formula and quotient recurrence are checked by exact
        # descriptor equality below.  These rows are the only reachable DAG
        # arithmetic needed by the outer ordinary-Z closure.
        183386: {"operation": "product", "arguments": [20, 183274]},
        183387: {"operation": "product", "arguments": [20, 183275]},
        183388: {"operation": "sum", "arguments": [183387, 183274]},
        183398: {"operation": "product", "arguments": [183386, 48]},
        183399: {"operation": "product", "arguments": [183388, 48]},
        183403: {"operation": "product", "arguments": [20, 183292]},
        183404: {"operation": "product", "arguments": [20, 183293]},
        183405: {"operation": "sum", "arguments": [183404, 183292]},
        183415: {"operation": "product", "arguments": [25222, 183403]},
        183416: {"operation": "product", "arguments": [25222, 183405]},
        183420: {"operation": "product", "arguments": [25251, 183293]},
        183421: {"operation": "product", "arguments": [25251, 183294]},
        183422: {"operation": "product", "arguments": [25247, 183293]},
        183423: {"operation": "sum", "arguments": [183421, 183422]},
        183436: {"operation": "sum", "arguments": [183415, 183420]},
        183437: {"operation": "sum", "arguments": [183416, 183423]},
        183441: {"operation": "negative", "argument": 183436},
        183442: {"operation": "negative", "argument": 183437},
        183446: {"operation": "sum", "arguments": [183398, 183441]},
        183447: {"operation": "sum", "arguments": [183399, 183442]},
        183451: {"operation": "positive_quotient", "numerator": 183446, "denominator": 25272, "source_positive_certificate": "same original L=1-delta*Z^2"},
        183452: {"operation": "product", "arguments": [25275, 183451]},
        183453: {"operation": "negative", "argument": 183452},
        183454: {"operation": "sum", "arguments": [183447, 183453]},
        183455: {"operation": "positive_quotient", "numerator": 183454, "denominator": 25272, "source_positive_certificate": "same original L=1-delta*Z^2"},
        204435: {"operation": "product", "arguments": [20, 204323]},
        204436: {"operation": "product", "arguments": [20, 204324]},
        204437: {"operation": "sum", "arguments": [204436, 204323]},
        204447: {"operation": "product", "arguments": [204435, 48]},
        204448: {"operation": "product", "arguments": [204437, 48]},
        204452: {"operation": "product", "arguments": [20, 204341]},
        204453: {"operation": "product", "arguments": [20, 204342]},
        204454: {"operation": "sum", "arguments": [204453, 204341]},
        204464: {"operation": "product", "arguments": [25222, 204452]},
        204465: {"operation": "product", "arguments": [25222, 204454]},
        204469: {"operation": "product", "arguments": [25251, 204342]},
        204470: {"operation": "product", "arguments": [25251, 204343]},
        204471: {"operation": "product", "arguments": [25247, 204342]},
        204472: {"operation": "sum", "arguments": [204470, 204471]},
        204485: {"operation": "sum", "arguments": [204464, 204469]},
        204486: {"operation": "sum", "arguments": [204465, 204472]},
        204490: {"operation": "negative", "argument": 204485},
        204491: {"operation": "negative", "argument": 204486},
        204495: {"operation": "sum", "arguments": [204447, 204490]},
        204496: {"operation": "sum", "arguments": [204448, 204491]},
        204500: {"operation": "positive_quotient", "numerator": 204495, "denominator": 25272, "source_positive_certificate": "same original L=1-delta*Z^2"},
        204501: {"operation": "product", "arguments": [25275, 204500]},
        204502: {"operation": "negative", "argument": 204501},
        204503: {"operation": "sum", "arguments": [204496, 204502]},
        204504: {"operation": "positive_quotient", "numerator": 204503, "denominator": 25272, "source_positive_certificate": "same original L=1-delta*Z^2"},
    }
    return rows


def _constant_closure(nodes: Sequence[Any], root: int) -> tuple[bool, list[dict[str, Any]], set[int]]:
    """Check that a source constant is Z-free by recursively inspecting rows."""

    errors: list[dict[str, Any]] = []
    seen: set[int] = set()

    def visit(index: int) -> None:
        if index in seen:
            return
        seen.add(index)
        try:
            row = _row(nodes, index)
            op = _operation(row, index)
        except PhysicalQyDerivativeProvenanceError as exc:
            errors.append({"gate": "constant_source_closure", "node": index, "error": str(exc)})
            return
        if op == "exact_rational":
            numerator = row.get("numerator", row.get("value"))
            denominator = row.get("denominator", 1)
            if isinstance(numerator, bool) or not isinstance(numerator, int) or isinstance(denominator, bool) or not isinstance(denominator, int) or denominator == 0:
                errors.append({"gate": "constant_source_closure", "node": index, "error": "malformed exact rational"})
            return
        if op == "exact_positive_integer_expression":
            if row.get("Z_independent") is not True or row.get("exact_positive_integer") is not True:
                errors.append({"gate": "constant_source_closure", "node": index, "error": "positive integer source is not explicitly Z independent"})
            return
        if op in {"sum", "product"}:
            children = row.get("arguments")
        elif op in {"negative", "analytic_unary"}:
            children = [row.get("argument")]
        elif op in {"positive_quotient", "quotient"}:
            children = [row.get("numerator"), row.get("denominator")]
            if op == "positive_quotient" and not row.get("source_positive_certificate"):
                errors.append({"gate": "constant_source_closure", "node": index, "error": "positive quotient lacks source certificate"})
        elif op == "function_substitution":
            if row.get("Z_independent_substitution") is not True:
                errors.append({"gate": "constant_source_closure", "node": index, "error": "unguarded constant substitution"})
            children = [row.get("expression"), row.get("variable"), row.get("value")]
        else:
            errors.append({"gate": "constant_source_closure", "node": index, "error": f"unhandled source operation {op}"})
            return
        if not isinstance(children, Sequence) or isinstance(children, (str, bytes, bytearray)):
            errors.append({"gate": "constant_source_closure", "node": index, "error": "source children are not a sequence"})
            return
        for child in children:
            if not isinstance(child, int) or isinstance(child, bool) or child < 0 or child >= len(nodes) or child >= index:
                errors.append({"gate": "constant_source_closure", "node": index, "error": f"invalid backward child {child!r}"})
                continue
            visit(child)

    visit(root)
    return not errors, errors, seen


def _validate_graph(nodes: Sequence[Any], expected_nodes_digest: str) -> dict[str, Any]:
    errors: list[dict[str, Any]] = []
    try:
        expected = _digest(expected_nodes_digest, "expected_nodes_digest")
        actual = canonical_digest(nodes)
    except (TypeError, ValueError, PhysicalQyDerivativeProvenanceError) as exc:
        return {
            "passed": False,
            "expected_nodes_digest": str(expected_nodes_digest).lower(),
            "nodes_digest": None,
            "graph_node_count": len(nodes) if isinstance(nodes, Sequence) else None,
            "prefix_digest": None,
            "errors": [{"gate": "graph_digest", "error": str(exc)}],
        }
    if actual != expected:
        errors.append({"gate": "graph_digest", "expected": expected, "actual": actual})
    if expected != EXPECTED_GRAPH_DIGEST:
        errors.append({"gate": "fixed_graph_digest", "expected": EXPECTED_GRAPH_DIGEST, "actual": expected})
    if len(nodes) != GRAPH_NODE_COUNT:
        errors.append({"gate": "graph_node_count", "expected": GRAPH_NODE_COUNT, "actual": len(nodes)})
    prefix = canonical_digest(nodes[:DERIVATIVE_PREFIX])
    if prefix != EXPECTED_FINITE_PREFIX_DIGEST:
        errors.append({"gate": "finite_prefix_digest", "expected": EXPECTED_FINITE_PREFIX_DIGEST, "actual": prefix})
    return {
        "passed": not errors,
        "expected_nodes_digest": expected,
        "nodes_digest": actual,
        "graph_node_count": len(nodes),
        "prefix_digest": prefix,
        "errors": errors,
    }


def _validate_descriptors(nodes: Sequence[Any]) -> dict[str, Any]:
    errors: list[dict[str, Any]] = []
    for index, expected in _source_rows().items():
        _expected_row(nodes, index, expected, errors, "source_descriptor")
    physical = _physical_rows()
    for index, expected in physical.items():
        _expected_row(nodes, index, expected, errors, "physical_producer_descriptor")
    for index, expected in _qy_rows().items():
        _expected_row(nodes, index, expected, errors, "Qy_arithmetic_descriptor")
    closures: dict[str, Any] = {}
    for name, root in {
        "delta": COMMON_DELTA_NODE,
        "one_minus_delta": COMMON_DELTA_ONE_NODE,
        "Pstar": PSTAR_NODE,
    }.items():
        passed, closure_errors, seen = _constant_closure(nodes, root)
        closures[name] = {"root": root, "passed": passed, "used_nodes": sorted(seen), "errors": closure_errors}
        errors.extend(closure_errors)
    return {"passed": not errors, "errors": errors, "constant_closures": closures}


def _context_from_row(row: Mapping[str, Any], side: str, contexts: Mapping[str, Any]) -> tuple[int | None, int | None]:
    raw = row.get("context")
    if isinstance(raw, str):
        raw = contexts.get(raw)
    if raw is None:
        raw = contexts.get(side)
    if not isinstance(raw, Mapping):
        return None, None
    return raw.get("coordinate"), raw.get("endpoint")


def _receipt_rows(receipt: Mapping[str, Any]) -> list[dict[str, Any]]:
    raw = receipt.get("rows")
    if isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
        return [dict(item) for item in raw if isinstance(item, Mapping)]
    if not isinstance(raw, Mapping):
        return []
    result: list[dict[str, Any]] = []
    for side, side_data in raw.items():
        if not isinstance(side_data, Mapping):
            continue
        for field, field_data in side_data.items():
            if isinstance(field_data, Mapping):
                if any(key in field_data for key in ("source_handle", "ordinary_order", "order")):
                    row = dict(field_data)
                    row.setdefault("side", side)
                    row.setdefault("field", field)
                    result.append(row)
                    continue
                field_data = [field_data]
            if isinstance(field_data, Sequence) and not isinstance(field_data, (str, bytes, bytearray)):
                for order, item in enumerate(field_data):
                    if not isinstance(item, Mapping):
                        continue
                    row = dict(item)
                    row.setdefault("side", side)
                    row.setdefault("field", field)
                    row.setdefault("ordinary_order", order)
                    result.append(row)
    return result


def _validate_primitive_receipt(
    receipt: Any,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    errors: list[dict[str, Any]] = []
    rows_out: list[dict[str, Any]] = []
    if not isinstance(receipt, Mapping):
        return {"passed": False, "rows": rows_out, "errors": [{"gate": "primitive_receipt", "error": "receipt is missing or not a mapping"}]}
    if receipt.get("schema") != PRIMITIVE_SCHEMA:
        errors.append({"gate": "primitive_schema", "expected": PRIMITIVE_SCHEMA, "actual": receipt.get("schema")})
    if receipt.get("source_only") is not True:
        errors.append({"gate": "primitive_source_only", "error": "source_only must be true"})
    if receipt.get("Z_variable") != Z_NODE:
        errors.append({"gate": "primitive_Z_variable", "expected": Z_NODE, "actual": receipt.get("Z_variable")})
    try:
        receipt_digest = _digest(receipt.get("nodes_digest"), "primitive_receipt.nodes_digest")
        if receipt_digest != expected_nodes_digest:
            errors.append({"gate": "primitive_nodes_digest", "expected": expected_nodes_digest, "actual": receipt_digest})
    except PhysicalQyDerivativeProvenanceError as exc:
        errors.append({"gate": "primitive_nodes_digest", "error": str(exc)})
    contexts = receipt.get("contexts")
    if not isinstance(contexts, Mapping):
        errors.append({"gate": "primitive_contexts", "error": "contexts are missing"})
        contexts = {}
    else:
        for side, expected in EXPECTED_CONTEXTS.items():
            actual = contexts.get(side)
            if not isinstance(actual, Mapping) or actual.get("coordinate") != expected["coordinate"] or actual.get("endpoint") != expected["endpoint"]:
                errors.append({"gate": "primitive_context", "side": side, "expected": expected, "actual": actual})
            if isinstance(actual, Mapping) and (actual.get("changed_Z") is True or actual.get("changed_native_environment") is True):
                errors.append({"gate": "primitive_context_changed", "side": side})
    expected_keys = {
        (side, field, order)
        for side in SIDES
        for field, orders in (("U_y", (0,)), ("physical_mean_y", (0, 1)))
        for order in orders
    }
    seen: set[tuple[str, str, int]] = set()
    for row in _receipt_rows(receipt):
        side = row.get("side")
        field = row.get("field")
        order = row.get("ordinary_order", row.get("order"))
        key = (side, field, order)
        if key in seen:
            errors.append({"gate": "primitive_duplicate_row", "row": key})
            continue
        seen.add(key)
        expected_source = EXPECTED_SOURCE_EDGES.get(side, {}).get(field, {}).get(order)
        expected_derivative = EXPECTED_SOURCE_DERIVATIVES.get(key)
        row_status: dict[str, Any] = {
            "side": side,
            "field": field,
            "ordinary_order": order,
            "source_handle": row.get("source_handle"),
            "derivative_handle": row.get("derivative_handle"),
            "expected_source_handle": expected_source,
            "expected_derivative_handle": expected_derivative,
            "source_derivative_proved": row.get("source_derivative_proved") is True,
            "regularity_proved": row.get("regularity_proved") is True,
            "context": row.get("context", side),
        }
        rows_out.append(row_status)
        if key not in expected_keys:
            errors.append({"gate": "primitive_unexpected_row", "row": key})
            continue
        if row.get("source_handle") != expected_source or row.get("derivative_handle") != expected_derivative:
            errors.append({"gate": "primitive_handle_binding", "row": key, "expected_source": expected_source, "expected_derivative": expected_derivative})
        if row.get("source_derivative_proved") is not True:
            errors.append({"gate": "primitive_source_derivative", "row": key})
        if row.get("regularity_proved") is not True:
            errors.append({"gate": "primitive_regularity", "row": key})
        coordinate, endpoint = _context_from_row(row, side, contexts)
        expected_context = EXPECTED_CONTEXTS.get(side)
        if coordinate != (expected_context or {}).get("coordinate") or endpoint != (expected_context or {}).get("endpoint"):
            errors.append({"gate": "primitive_row_context", "row": key, "coordinate": coordinate, "endpoint": endpoint})
        if row.get("changed_Z") is True or row.get("changed_native_environment") is True or row.get("differentiated_opaque_alias") is True:
            errors.append({"gate": "primitive_alias_or_context", "row": key})
        if any("alias" in str(k).lower() for k in row):
            errors.append({"gate": "primitive_alias_field", "row": key})
    missing = sorted(expected_keys - seen, key=repr)
    for key in missing:
        errors.append({"gate": "primitive_missing_row", "row": key})
    if set(seen) != expected_keys:
        errors.append({"gate": "primitive_row_set", "expected": sorted(expected_keys, key=repr), "actual": sorted(seen, key=repr)})
    if receipt.get("passed") is not True:
        errors.append({"gate": "primitive_receipt_passed", "error": "receipt.passed is not true"})
    if receipt.get("differentiated_opaque_alias") is True:
        errors.append({"gate": "primitive_opaque_alias", "error": "opaque alias differentiation is disallowed"})
    return {
        "passed": not errors,
        "schema": receipt.get("schema"),
        "source_only": receipt.get("source_only") is True,
        "Z_variable": receipt.get("Z_variable"),
        "nodes_digest": receipt.get("nodes_digest"),
        "contexts": {side: dict(contexts.get(side, {})) if isinstance(contexts.get(side), Mapping) else None for side in SIDES},
        "rows": rows_out,
        "errors": errors,
        "input_digest": canonical_digest(receipt),
    }


def _validate_qy0_receipt(receipt: Any, expected_nodes_digest: str) -> dict[str, Any]:
    errors: list[dict[str, Any]] = []
    if not isinstance(receipt, Mapping):
        return {"passed": False, "all_Z_identity_proved": False, "errors": [{"gate": "qy0_receipt", "error": "receipt is missing or not a mapping"}]}
    if receipt.get("schema") != QY0_SCHEMA:
        errors.append({"gate": "qy0_schema", "expected": QY0_SCHEMA, "actual": receipt.get("schema")})
    if receipt.get("passed") is not True or receipt.get("source_only") is not True:
        errors.append({"gate": "qy0_passed_source_only"})
    if receipt.get("differentiated_opaque_alias") is True:
        errors.append({"gate": "qy0_opaque_alias"})
    try:
        packet_digest = canonical_digest(receipt)
    except (TypeError, ValueError) as exc:
        packet_digest = None
        errors.append({"gate": "qy0_packet_digest", "error": str(exc)})
    if packet_digest != EXPECTED_QY0_PACKET_DIGEST:
        errors.append({"gate": "qy0_packet_digest", "expected": EXPECTED_QY0_PACKET_DIGEST, "actual": packet_digest})
    for flag in (
        "numeric_values_computed",
        "pressure_quadrature_executed",
        "old_mathematical_operator_replayed",
        "accepted_graph_replayed",
        "differentiated_opaque_alias",
        "full_RQy_germ_Z1_proved",
        "physical_RQy_Z1_proved",
    ):
        if _contains_true_flag(receipt, flag):
            errors.append({"gate": "qy0_contradictory_flag", "flag": flag})
    try:
        if _digest(receipt.get("nodes_digest"), "qy0_receipt.nodes_digest") != expected_nodes_digest:
            errors.append({"gate": "qy0_nodes_digest"})
    except PhysicalQyDerivativeProvenanceError as exc:
        errors.append({"gate": "qy0_nodes_digest", "error": str(exc)})
    formula = receipt.get("physical_formula")
    if not isinstance(formula, Mapping) or formula.get("mean_y_source") != "U-mean" or formula.get("radial_mean_y_used") is not False:
        errors.append({"gate": "qy0_mean_y_source"})
    alias_policy = receipt.get("alias_policy")
    if not isinstance(alias_policy, Mapping) or alias_policy.get("source_bound_endpoint_alias_only") is not True or alias_policy.get("ordinary_Z1_alias_is_independent") is not True or alias_policy.get("ordinary_Z1_was_not_differentiated") is not True:
        errors.append({"gate": "qy0_alias_policy"})
    pairs = receipt.get("canonical_value_pairs", receipt.get("value_pairs"))
    pair = None
    if isinstance(pairs, Sequence) and not isinstance(pairs, (str, bytes, bytearray)):
        for candidate in pairs:
            if isinstance(candidate, Mapping) and candidate.get("field") == "Q_y" and candidate.get("ordinary_order") == 0:
                pair = candidate
                break
    if pair is None:
        errors.append({"gate": "qy0_pair", "error": "Q_y ordinary order 0 pair is missing"})
    else:
        if pair.get("postrestore_target_handle") != EXPECTED_QY0["postrestore"] or pair.get("patch_target_handle") != EXPECTED_QY0["patch"]:
            errors.append({"gate": "qy0_pair_targets"})
        if pair.get("equal") is not True or pair.get("status") != "proved":
            errors.append({"gate": "qy0_pair_proved"})
        post = pair.get("postrestore")
        patch = pair.get("patch")
        if not isinstance(post, Mapping) or not isinstance(patch, Mapping) or post.get("srepr") != patch.get("srepr") or "Z" not in str(post.get("srepr", "")):
            errors.append({"gate": "qy0_all_Z_witness"})
    ports = receipt.get("physical_ports")
    for side, expected in EXPECTED_CONTEXTS.items():
        chart = ports.get(side, {}).get("chart") if isinstance(ports, Mapping) and isinstance(ports.get(side), Mapping) else None
        if not isinstance(chart, Mapping) or chart.get("coordinate") != expected["coordinate"] or chart.get("endpoint") != expected["endpoint"]:
            errors.append({"gate": "qy0_context", "side": side})
    return {
        "passed": not errors,
        "all_Z_identity_proved": not errors,
        "schema": receipt.get("schema"),
        "errors": errors,
        "input_digest": packet_digest,
        "external_mean_premise": receipt.get("provenance", {}).get("modified_mean_external_premise") if isinstance(receipt.get("provenance"), Mapping) else None,
    }


def _validate_common_phase(receipt: Any) -> dict[str, Any]:
    errors: list[dict[str, Any]] = []
    if not isinstance(receipt, Mapping):
        return {"passed": False, "errors": [{"gate": "common_phase", "error": "receipt is missing or not a mapping"}]}
    try:
        packet_digest = canonical_digest(receipt)
    except (TypeError, ValueError) as exc:
        packet_digest = None
        errors.append({"gate": "common_phase_packet_digest", "error": str(exc)})
    if packet_digest != EXPECTED_COMMON_PHASE_PACKET_DIGEST:
        errors.append({"gate": "common_phase_packet_digest", "expected": EXPECTED_COMMON_PHASE_PACKET_DIGEST, "actual": packet_digest})
    if receipt.get("schema") != COMMON_PHASE_SCHEMA:
        errors.append({"gate": "common_phase_schema"})
    for key in ("passed", "source_only", "actual_common_phase_proved", "actual_phase_Z_independence_proved"):
        if receipt.get(key) is not True:
            errors.append({"gate": "common_phase_flag", "key": key})
    for flag in (
        "fractional_part_differentiated",
        "numerical_phase_value_computed",
        "inverse_conditioning_proved",
        "native_loop_derivative_proved",
        "density_seam_proved",
        "full_goal_complete",
    ):
        if _contains_true_flag(receipt, flag):
            errors.append({"gate": "common_phase_contradictory_flag", "flag": flag})
    positive_endpoint_guard = receipt.get("coordinates") == [POST_COORDINATE, PATCH_COORDINATE] and receipt.get("endpoints") == [1, 1]
    if not positive_endpoint_guard:
        errors.append({"gate": "common_phase_context"})
    phase_rows = receipt.get("phase_rows")
    if not isinstance(phase_rows, Sequence) or isinstance(phase_rows, (str, bytes, bytearray)) or len(phase_rows) != 2:
        positive_endpoint_guard = False
        errors.append({"gate": "common_phase_positive_endpoint_guard", "error": "two phase rows are required"})
    else:
        expected_rows = (
            {"coordinate": POST_COORDINATE, "endpoint": EXPECTED_ENDPOINT},
            {"coordinate": PATCH_COORDINATE, "endpoint": EXPECTED_ENDPOINT},
        )
        for offset, (row, expected) in enumerate(zip(phase_rows, expected_rows)):
            if not isinstance(row, Mapping) or row.get("coordinate") != expected["coordinate"] or row.get("endpoint") != expected["endpoint"]:
                positive_endpoint_guard = False
                errors.append({"gate": "common_phase_positive_endpoint_guard", "row": offset, "expected": expected, "actual": row})
    return {"passed": not errors, "positive_endpoint_guard_proved": positive_endpoint_guard, "actual_Z_independence_proved": receipt.get("actual_phase_Z_independence_proved") is True, "errors": errors, "input_digest": packet_digest}


@dataclass(frozen=True)
class _Poly:
    terms: tuple[tuple[tuple[int, ...], Fraction], ...]


def _poly_from(terms: Mapping[tuple[int, ...], Fraction]) -> _Poly:
    clean = [(tuple(monomial), coefficient) for monomial, coefficient in terms.items() if coefficient]
    clean.sort(key=lambda item: item[0])
    return _Poly(tuple(clean))


def _poly_const(value: Fraction | int = 0) -> _Poly:
    fraction = value if isinstance(value, Fraction) else Fraction(value, 1)
    return _poly_from({(): fraction})


def _poly_var(handle: int) -> _Poly:
    return _poly_from({(handle,): Fraction(1, 1)})


def _poly_add(left: _Poly, right: _Poly) -> _Poly:
    terms: dict[tuple[int, ...], Fraction] = {}
    for monomial, coefficient in left.terms + right.terms:
        terms[monomial] = terms.get(monomial, Fraction(0, 1)) + coefficient
    return _poly_from(terms)


def _poly_neg(value: _Poly) -> _Poly:
    return _poly_from({monomial: -coefficient for monomial, coefficient in value.terms})


def _poly_mul(left: _Poly, right: _Poly) -> _Poly:
    terms: dict[tuple[int, ...], Fraction] = {}
    for left_monomial, left_coefficient in left.terms:
        for right_monomial, right_coefficient in right.terms:
            monomial = tuple(sorted(left_monomial + right_monomial))
            terms[monomial] = terms.get(monomial, Fraction(0, 1)) + left_coefficient * right_coefficient
    return _poly_from(terms)


@dataclass(frozen=True)
class _Rat:
    numerator: _Poly
    denominator: _Poly


def _rat_const(value: Fraction | int = 0) -> _Rat:
    return _Rat(_poly_const(value), _poly_const(1))


def _rat_atom(handle: int) -> _Rat:
    return _Rat(_poly_var(handle), _poly_const(1))


def _rat_add(left: _Rat, right: _Rat) -> _Rat:
    return _Rat(_poly_add(_poly_mul(left.numerator, right.denominator), _poly_mul(right.numerator, left.denominator)), _poly_mul(left.denominator, right.denominator))


def _rat_neg(value: _Rat) -> _Rat:
    return _Rat(_poly_neg(value.numerator), value.denominator)


def _rat_sub(left: _Rat, right: _Rat) -> _Rat:
    return _rat_add(left, _rat_neg(right))


def _rat_mul(left: _Rat, right: _Rat) -> _Rat:
    return _Rat(_poly_mul(left.numerator, right.numerator), _poly_mul(left.denominator, right.denominator))


def _rat_div(left: _Rat, right: _Rat) -> _Rat:
    if not right.numerator.terms:
        raise PhysicalQyDerivativeProvenanceError("symbolic quotient has a zero numerator denominator")
    return _Rat(_poly_mul(left.numerator, right.denominator), _poly_mul(left.denominator, right.numerator))


def _rat_equal(left: _Rat, right: _Rat) -> bool:
    return _poly_mul(left.numerator, right.denominator) == _poly_mul(right.numerator, left.denominator)


@dataclass(frozen=True)
class _Dual:
    value: _Rat
    derivative: _Rat


class _ArithmeticEvaluator:
    """Exact rational-polynomial dual evaluator over an existing graph."""

    def __init__(self, nodes: Sequence[Any], source_derivatives: Mapping[int, int]) -> None:
        self.nodes = nodes
        self.source_derivatives = dict(source_derivatives)
        self.used_nodes: set[int] = set()
        self.operation_counts: Counter[str] = Counter()

    def _record(self, index: int, operation: str) -> None:
        self.used_nodes.add(index)
        self.operation_counts[operation] += 1

    def _row_op(self, index: int) -> tuple[Mapping[str, Any], str]:
        row = _row(self.nodes, index)
        operation = _operation(row, index)
        self._record(index, operation)
        return row, operation

    def _children(self, row: Mapping[str, Any], key: str, index: int) -> list[int]:
        raw = row.get(key)
        if not isinstance(raw, Sequence) or isinstance(raw, (str, bytes, bytearray)):
            raise PhysicalQyDerivativeProvenanceError(f"node {index}.{key} is not a child sequence")
        children = list(raw)
        if not all(isinstance(child, int) and not isinstance(child, bool) for child in children):
            raise PhysicalQyDerivativeProvenanceError(f"node {index}.{key} has a non-integer child")
        return children

    def _rational(self, row: Mapping[str, Any], index: int) -> _Rat:
        numerator = row.get("numerator", row.get("value"))
        denominator = row.get("denominator", 1)
        if isinstance(numerator, bool) or not isinstance(numerator, int) or isinstance(denominator, bool) or not isinstance(denominator, int) or denominator == 0:
            raise PhysicalQyDerivativeProvenanceError(f"node {index} is not an exact rational")
        return _rat_const(Fraction(numerator, denominator))

    def value(self, index: int) -> _Rat:
        if index in self.source_derivatives or index in EXPECTED_SOURCE_EDGES["postrestore"]["U_y"].values() or index in EXPECTED_SOURCE_EDGES["patch"]["U_y"].values() or index in {183292, 183293, 183294, 204341, 204342, 204343}:
            self.used_nodes.add(index)
            return _rat_atom(index)
        if index == Z_NODE:
            # ``value`` traverses the actual Q_y[Z1] source row.  Z is a
            # source atom there; only ``dual`` assigns it derivative one.
            self.used_nodes.add(index)
            return _rat_atom(index)
        if index in {COMMON_DELTA_ONE_NODE, COMMON_D_NODE, COMMON_D1_NODE, COMMON_L_NODE, COMMON_L1_NODE, COMMON_DELTA_NODE, PSTAR_NODE}:
            self.used_nodes.add(index)
            return _rat_atom(index)
        row, operation = self._row_op(index)
        if operation == "exact_rational":
            return self._rational(row, index)
        if operation == "sum":
            result = _rat_const(0)
            for child in self._children(row, "arguments", index):
                result = _rat_add(result, self.value(child))
            return result
        if operation == "product":
            result = _rat_const(1)
            for child in self._children(row, "arguments", index):
                result = _rat_mul(result, self.value(child))
            return result
        if operation == "negative":
            return _rat_neg(self.value(row["argument"]))
        if operation in {"positive_quotient", "quotient"}:
            if operation == "positive_quotient" and not row.get("source_positive_certificate"):
                raise PhysicalQyDerivativeProvenanceError(f"node {index} positive quotient lacks source certificate")
            return _rat_div(self.value(row["numerator"]), self.value(row["denominator"]))
        raise PhysicalQyDerivativeProvenanceError(f"unsupported arithmetic operation {operation} at node {index}")

    def dual(self, index: int) -> _Dual:
        if index in self.source_derivatives:
            self.used_nodes.add(index)
            return _Dual(_rat_atom(index), _rat_atom(self.source_derivatives[index]))
        if index == Z_NODE:
            self.used_nodes.add(index)
            return _Dual(_rat_atom(index), _rat_const(1))
        if index == COMMON_DELTA_ONE_NODE or index == COMMON_DELTA_NODE or index == PSTAR_NODE:
            self.used_nodes.add(index)
            return _Dual(_rat_atom(index), _rat_const(0))
        if index == COMMON_D_NODE:
            self.used_nodes.add(index)
            return _Dual(_rat_atom(index), _rat_atom(COMMON_D1_NODE))
        if index == COMMON_L_NODE:
            self.used_nodes.add(index)
            return _Dual(_rat_atom(index), _rat_atom(COMMON_L1_NODE))
        if index in {COMMON_D1_NODE, COMMON_L1_NODE}:
            self.used_nodes.add(index)
            return _Dual(_rat_atom(index), _rat_const(0))
        row, operation = self._row_op(index)
        if operation == "exact_rational":
            return _Dual(self._rational(row, index), _rat_const(0))
        if operation == "sum":
            result = _Dual(_rat_const(0), _rat_const(0))
            for child in self._children(row, "arguments", index):
                result = _dual_add(result, self.dual(child))
            return result
        if operation == "product":
            result = _Dual(_rat_const(1), _rat_const(0))
            for child in self._children(row, "arguments", index):
                result = _dual_mul(result, self.dual(child))
            return result
        if operation == "negative":
            return _dual_neg(self.dual(row["argument"]))
        if operation in {"positive_quotient", "quotient"}:
            if operation == "positive_quotient" and not row.get("source_positive_certificate"):
                raise PhysicalQyDerivativeProvenanceError(f"node {index} positive quotient lacks source certificate")
            return _dual_div(self.dual(row["numerator"]), self.dual(row["denominator"]))
        raise PhysicalQyDerivativeProvenanceError(f"unsupported dual operation {operation} at node {index}")


def _dual_add(left: _Dual, right: _Dual) -> _Dual:
    return _Dual(_rat_add(left.value, right.value), _rat_add(left.derivative, right.derivative))


def _dual_mul(left: _Dual, right: _Dual) -> _Dual:
    return _Dual(_rat_mul(left.value, right.value), _rat_add(_rat_mul(left.derivative, right.value), _rat_mul(left.value, right.derivative)))


def _dual_neg(value: _Dual) -> _Dual:
    return _Dual(_rat_neg(value.value), _rat_neg(value.derivative))


def _dual_div(left: _Dual, right: _Dual) -> _Dual:
    value = _rat_div(left.value, right.value)
    numerator = _rat_sub(_rat_mul(left.derivative, right.value), _rat_mul(left.value, right.derivative))
    denominator = _rat_mul(right.value, right.value)
    derivative = _rat_div(numerator, denominator)
    return _Dual(value, derivative)


def _rat_json(value: _Rat) -> dict[str, Any]:
    def poly_json(poly: _Poly) -> list[dict[str, Any]]:
        return [
            {"monomial": list(monomial), "coefficient": [coefficient.numerator, coefficient.denominator]}
            for monomial, coefficient in poly.terms
        ]

    return {"numerator": poly_json(value.numerator), "denominator": poly_json(value.denominator)}


def _side_derivative(
    nodes: Sequence[Any],
    side: str,
) -> dict[str, Any]:
    source = EXPECTED_SOURCE_EDGES[side]
    derivatives = {
        source["U_y"][0]: source["U_y"][1],
        source["physical_mean_y"][0]: source["physical_mean_y"][1],
        source["physical_mean_y"][1]: source["physical_mean_y"][2],
    }
    evaluator = _ArithmeticEvaluator(nodes, derivatives)
    try:
        actual = evaluator.dual(EXPECTED_QY0[side])
        target = evaluator.value(EXPECTED_QY1[side])
        equal = _rat_equal(actual.derivative, target)
        return {
            "side": side,
            "Q_y_Z0_handle": EXPECTED_QY0[side],
            "Q_y_Z1_handle": EXPECTED_QY1[side],
            "source_edges": {field: {str(order): node for order, node in rows.items()} for field, rows in source.items()},
            "source_derivative_edges": {str(key): value for key, value in derivatives.items()},
            "D_Z_Q_y_Z0": _rat_json(actual.derivative),
            "actual_Q_y_Z1": _rat_json(target),
            "derivative_equals_actual": equal,
            "status": "proved" if equal else "failed",
            "used_node_witness": sorted(evaluator.used_nodes),
            "used_node_digest": canonical_digest(sorted(evaluator.used_nodes)),
            "operation_counts": dict(sorted(evaluator.operation_counts.items())),
            "error": None,
        }
    except (KeyError, TypeError, ValueError, PhysicalQyDerivativeProvenanceError) as exc:
        return {
            "side": side,
            "Q_y_Z0_handle": EXPECTED_QY0[side],
            "Q_y_Z1_handle": EXPECTED_QY1[side],
            "source_edges": {field: {str(order): node for order, node in rows.items()} for field, rows in source.items()},
            "source_derivative_edges": {str(key): value for key, value in derivatives.items()},
            "D_Z_Q_y_Z0": None,
            "actual_Q_y_Z1": None,
            "derivative_equals_actual": False,
            "status": "unresolved",
            "used_node_witness": sorted(evaluator.used_nodes),
            "used_node_digest": canonical_digest(sorted(evaluator.used_nodes)),
            "operation_counts": dict(sorted(evaluator.operation_counts.items())),
            "error": str(exc),
        }


def derive(
    nodes: Sequence[Mapping[str, Any]],
    qy0_value_seam: Mapping[str, Any],
    primitive_Z_source_receipt: Mapping[str, Any] | None,
    common_phase_source: Mapping[str, Any],
    *,
    expected_nodes_digest: str,
) -> dict[str, Any]:
    """Check the conditional source-to-``D_Z`` closure for physical ``Q_y``.

    The six primitive receipt rows are the only source derivatives admitted:
    ``U_y[0] -> U_y[1]`` and ``mean_y[0] -> mean_y[1] -> mean_y[2]`` on each
    chart.  The arithmetic evaluator remains conditional when that receipt is
    absent or open, but no VALUE alias receives a derivative.
    """

    graph = _validate_graph(nodes, expected_nodes_digest)
    graph_nodes_ok = graph["passed"]
    descriptors = _validate_descriptors(nodes) if isinstance(nodes, Sequence) else {"passed": False, "errors": [{"gate": "nodes", "error": "nodes is not a sequence"}], "constant_closures": {}}
    primitive = _validate_primitive_receipt(primitive_Z_source_receipt, EXPECTED_GRAPH_DIGEST)
    qy0 = _validate_qy0_receipt(qy0_value_seam, EXPECTED_GRAPH_DIGEST)
    common = _validate_common_phase(common_phase_source)
    side_results = {side: _side_derivative(nodes, side) if graph_nodes_ok and descriptors["passed"] else {"side": side, "status": "unresolved", "derivative_equals_actual": False, "error": "graph or source descriptor gate failed", "Q_y_Z0_handle": EXPECTED_QY0[side], "Q_y_Z1_handle": EXPECTED_QY1[side], "source_edges": EXPECTED_SOURCE_EDGES[side]} for side in SIDES}
    side_derivatives_proved = all(result.get("derivative_equals_actual") is True for result in side_results.values())
    outer_calculus_proved = bool(graph_nodes_ok and descriptors["passed"] and side_derivatives_proved)
    primitive_contracts_passed = primitive.get("passed") is True
    qy0_identity_proved = qy0.get("all_Z_identity_proved") is True
    common_phase_passed = common.get("passed") is True
    passed = bool(outer_calculus_proved and primitive_contracts_passed and qy0_identity_proved and common_phase_passed)
    unresolved: list[dict[str, Any]] = []
    unresolved.extend(graph.get("errors", []))
    unresolved.extend(descriptors.get("errors", []))
    unresolved.extend(primitive.get("errors", []))
    unresolved.extend(qy0.get("errors", []))
    unresolved.extend(common.get("errors", []))
    return {
        "schema": SCHEMA,
        "passed": passed,
        "source_only": True,
        "actual_source_graph_inspected": True,
        "actual_saved_graph_run": False,
        "numeric_values_computed": False,
        "numeric_field_evaluated": False,
        "pressure_quadrature_executed": False,
        "old_mathematical_operator_replayed": False,
        "old_operators_replayed": False,
        "accepted_graph_replayed": False,
        "differentiated_opaque_alias": False,
        "new_value_admission": False,
        "expected_nodes_digest": graph.get("expected_nodes_digest"),
        "nodes_digest": graph.get("nodes_digest"),
        "graphdigest": graph.get("nodes_digest"),
        "graph_node_count": graph.get("graph_node_count"),
        "derivative_source_prefix_count": DERIVATIVE_PREFIX,
        "derivative_source_prefix_digest": graph.get("prefix_digest"),
        "outer_calculus_proved": outer_calculus_proved,
        "ordinary_Q_y_Z0_identity_admitted": qy0_identity_proved,
        "physical_Qy_Z1_value_seam_proved": passed,
        "ordinary_Q_y_Z1_value_congruence": passed,
        "C1_regularity_proved": primitive_contracts_passed,
        "C1_regular_proved": primitive_contracts_passed,
        "primitive_Z_source_receipt_passed": primitive_contracts_passed,
        "primitive_source_contracts_passed": primitive_contracts_passed,
        "common_phase_positive_endpoint_guard_proved": common.get("positive_endpoint_guard_proved") is True,
        "q_y_z1": {
            "field": "Q_y",
            "ordinary_order": 1,
            "postrestore_target_handle": EXPECTED_QY1["postrestore"],
            "patch_target_handle": EXPECTED_QY1["patch"],
            "status": "proved" if passed else ("conditional" if outer_calculus_proved else "open"),
            "equal": passed,
            "proved": passed,
            "derivative_provenance": "actual per-side D_Z(Q_y[0]) equals actual Q_y[1]",
            "cross_side_identity_rule": "differentiate admitted all-Z Q_y[0] identity under C1 regularity",
            "independent_Z2_value_proof": False,
            "reason": None if passed else "primitive source, C1, all-Z Q_y[0], or graph gate remains open",
        },
        "per_side": side_results,
        "primitive_source_contract": primitive,
        "qy0_value_seam": qy0,
        "common_phase_source": common,
        "source_descriptor_validation": descriptors,
        "conditional_obligations": {
            "primitive_schema": PRIMITIVE_SCHEMA,
            "required_rows": [
                {"side": side, "field": field, "ordinary_order": order, "source_handle": EXPECTED_SOURCE_EDGES[side][field][order], "derivative_handle": EXPECTED_SOURCE_DERIVATIVES.get((side, field, order))}
                for side in SIDES
                for field, orders in (("U_y", (0,)), ("physical_mean_y", (0, 1)))
                for order in orders
            ],
            "same_endpoint_native_context": EXPECTED_CONTEXTS,
            "Z_variable": Z_NODE,
            "no_opaque_alias_differentiation": True,
            "no_new_value_admission": True,
            "strict_constant_source_closure": descriptors.get("constant_closures", {}),
            "primitive_contracts_genuine_passed": primitive_contracts_passed,
            "C1_regularity_genuine_passed": primitive_contracts_passed,
            "admitted_all_Z_Qy0_identity": qy0_identity_proved,
            "common_phase_actual_Z_independence": common_phase_passed,
        },
        "actual_used_node_witness": {side: {"nodes": result.get("used_node_witness", []), "digest": result.get("used_node_digest"), "operations": result.get("operation_counts", {})} for side, result in side_results.items()},
        "input_data_digests": {
            "primitive_Z_source_receipt": primitive.get("input_digest"),
            "qy0_value_seam": qy0.get("input_digest"),
            "common_phase_source": common.get("input_digest"),
            "nodes": graph.get("nodes_digest"),
        },
        "provenance": {
            "producer_module": PRODUCER_MODULE,
            "producer_module_sha256": PRODUCER_MODULE_SHA256,
            "formula": "Q_y=div(2*Z*trim(U_y)-((1-delta)*Z*trim(U-mean)+(1-Z^2)*shift(U-mean)))",
            "quotient_rule": "D_Z(N/L)=(D_Z(N)*L-N*D_Z(L))/L^2; saved Q_y[1] is checked against this dual result",
            "positive_quotient_certificate": "same original L=1-delta*Z^2",
        },
        "unresolved_gates": unresolved,
        "full_RQy_germ_Z1_proved": False,
        "full_RQy_germ_proved": False,
        "physical_RQy_Z1_proved": False,
        "full_physical_join": False,
        "global_five_moment_closure": False,
        "Cartesian_NS": False,
        "full_goal_complete": False,
    }


__all__ = [
    "COMMON_D1_NODE",
    "COMMON_D_NODE",
    "COMMON_DELTA_ONE_NODE",
    "COMMON_L1_NODE",
    "COMMON_L_NODE",
    "EXPECTED_GRAPH_DIGEST",
    "EXPECTED_QY0",
    "EXPECTED_QY1",
    "EXPECTED_SOURCE_EDGES",
    "PhysicalQyDerivativeProvenanceError",
    "SCHEMA",
    "canonical_digest",
    "derive",
]
