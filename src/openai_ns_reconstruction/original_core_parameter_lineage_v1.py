"""Static lineage audit for the saved native core parameters.

This leaf accepts already decoded mappings and two producer source texts.  It
does not open a report, import a producer, execute source, or replay an
ancestor.  The claim is deliberately narrower than a live-object identity:
the native core ``delta`` and ``logPstar`` records, their corresponding native
actual records, and the two records carried by the previous parameter leaf
must have identical directed binary endpoints.  Display strings and other
interval metadata are retained and reported, but are not used as identity.
"""

from __future__ import annotations

import ast
import copy
import hashlib
from collections.abc import Mapping, Sequence
from typing import Any


SCHEMA = "original_core_parameter_lineage_v1"
PREVIOUS_SCHEMA = "original_lower_parameter_binding_v1"

# These are the immutable paths and hashes exposed by the published previous
# lower-parameter runner.  They authenticate provenance records already
# supplied by the caller; this leaf never opens any of the files.
PREVIOUS_NATIVE_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_directed_controls_v2.json"
)
PREVIOUS_NATIVE_SHA256 = (
    "10abe68b98fdf87709a7a7c9f27103523701e746c357ae62386038700052490a"
)
PREVIOUS_CORRECTED_PATH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_corrected_reserved_m15_source.json"
)
PREVIOUS_CORRECTED_SHA256 = (
    "fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5"
)
PREVIOUS_PRODUCER_PATH = "scripts/run_original_corrected_reserved_m15.py"
PREVIOUS_PRODUCER_SHA256 = (
    "9b41d60ab277656734cc98aa3272123a1ce383da4337ba85982dd41142350d6d"
)

Q_SOURCE_PATH = "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"
BACKGROUND_SOURCE_PATH = (
    "src/openai_ns_reconstruction/original_micro_background_I1_numeric.py"
)

CORE_PARAMETER_PATH = (
    "actual_core_signed_direction",
    "parameters",
)
NATIVE_DELTA_PATH = (
    "actual_micro_Q_increment_bounds",
    "actual_delta",
)
NATIVE_LOGPSTAR_PATH = (
    "actual_micro_background_I1_numeric",
    "physical_scale",
    "logPstar",
)


_MISSING = object()
_BINARY_KEYS = {
    "directed",
    "lower_exact_mpf_tuple",
    "upper_exact_mpf_tuple",
}


def _sha_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _copy(value: Any) -> Any:
    return copy.deepcopy(value)


def _path(value: Any, keys: Sequence[str]) -> Any:
    current = value
    for key in keys:
        if not isinstance(current, Mapping) or key not in current:
            return _MISSING
        current = current[key]
    return current


def _path_text(keys: Sequence[str]) -> str:
    return ".".join(keys)


def _tuple4(value: Any) -> tuple[int, int, int, int] | None:
    if not isinstance(value, (list, tuple)) or len(value) != 4:
        return None
    if any(type(item) is not int for item in value):
        return None
    return tuple(value)  # type: ignore[return-value]


def _record_info(value: Any, label: str) -> dict[str, Any]:
    """Decode only the exact binary fields, preserving the raw record."""

    info: dict[str, Any] = {
        "label": label,
        "raw": _copy(value) if value is not _MISSING else None,
        "present": value is not _MISSING,
        "valid": False,
        "binary": None,
        "errors": [],
    }
    if value is _MISSING:
        info["errors"].append("record is missing")
        return info
    if not isinstance(value, Mapping):
        info["errors"].append("record is not a mapping")
        return info
    if value.get("directed") is not True:
        info["errors"].append("directed flag is not true")
    lower = _tuple4(value.get("lower_exact_mpf_tuple"))
    upper = _tuple4(value.get("upper_exact_mpf_tuple"))
    if lower is None:
        info["errors"].append("lower_exact_mpf_tuple is not four integers")
    if upper is None:
        info["errors"].append("upper_exact_mpf_tuple is not four integers")
    if not info["errors"]:
        info["valid"] = True
        info["binary"] = {
            "directed": True,
            "lower_exact_mpf_tuple": list(lower),
            "upper_exact_mpf_tuple": list(upper),
        }
    return info


def _metadata_differences(left: Mapping[str, Any], right: Mapping[str, Any]) -> list[str]:
    differences: list[str] = []
    for key in sorted(set(left) | set(right)):
        if key in _BINARY_KEYS:
            continue
        if left.get(key, _MISSING) != right.get(key, _MISSING):
            differences.append(str(key))
    return differences


def _compare_records(label: str, left: Any, right: Any) -> dict[str, Any]:
    """Compare exact directed binary tuples, never interval overlap."""

    left_info = _record_info(left, label + ".left")
    right_info = _record_info(right, label + ".right")
    mismatches: list[dict[str, Any]] = []
    if not left_info["valid"] or not right_info["valid"]:
        mismatches.append({
            "kind": "invalid_record",
            "left_errors": list(left_info["errors"]),
            "right_errors": list(right_info["errors"]),
        })
    else:
        left_binary = left_info["binary"]
        right_binary = right_info["binary"]
        for key in ("directed", "lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
            if left_binary[key] != right_binary[key]:
                mismatches.append({
                    "kind": "binary_tuple_mismatch",
                    "field": key,
                    "left": _copy(left_binary[key]),
                    "right": _copy(right_binary[key]),
                })

    metadata: list[str] = []
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        metadata = _metadata_differences(left, right)
    exact = not mismatches
    return {
        "label": label,
        "exact_binary_tuple_equal": exact,
        "binary_tuple_comparison": "exact",
        "left": _copy(left) if left is not _MISSING else None,
        "right": _copy(right) if right is not _MISSING else None,
        "mismatches": mismatches,
        "nonbinary_metadata_differences": metadata,
        "metadata_differences_do_not_fail_exact_binary_identity": True,
        "rounding_or_serialization_difference_possible": bool(mismatches),
        "overlap_used_as_identity": False,
    }


def _unpack_native(report: Any) -> tuple[Any, str, list[str]]:
    errors: list[str] = []
    if not isinstance(report, Mapping):
        return _MISSING, "", ["native report is not a mapping"]
    if "actual_core_signed_direction" in report:
        return report, "", errors
    nested = report.get("actual_native_directed_controls_v2")
    if isinstance(nested, Mapping) and "actual_core_signed_direction" in nested:
        return nested, "actual_native_directed_controls_v2", errors
    errors.append("native report lacks actual_core_signed_direction")
    return _MISSING, "", errors


def _unpack_previous(receipt: Any) -> tuple[Any, str, Any, list[str]]:
    errors: list[str] = []
    if not isinstance(receipt, Mapping):
        return _MISSING, "", _MISSING, ["previous parameter receipt is not a mapping"]
    if "delta" in receipt and "actual_logPstar_parameter" in receipt:
        return receipt, "", receipt.get("source_provenance", _MISSING), errors
    result = receipt.get("result")
    if isinstance(result, Mapping):
        leaf = result.get("leaf")
        if isinstance(leaf, Mapping):
            return leaf, "result.leaf", result.get("source_provenance", _MISSING), errors
        if "delta" in result and "actual_logPstar_parameter" in result:
            return result, "result", result.get("source_provenance", _MISSING), errors
        nested = result.get("result")
        if isinstance(nested, Mapping):
            leaf = nested.get("leaf")
            if isinstance(leaf, Mapping):
                return leaf, "result.result.leaf", nested.get("source_provenance", _MISSING), errors
    errors.append("previous receipt lacks the published leaf result")
    return _MISSING, "", _MISSING, errors


def _subscript_key(node: ast.Subscript) -> str | None:
    value: Any = node.slice
    if isinstance(value, ast.Index):  # pragma: no cover - old Python AST
        value = value.value
    if isinstance(value, ast.Constant) and type(value.value) is str:
        return value.value
    return None


def _subscript_path(node: ast.AST) -> tuple[str, ...] | None:
    keys: list[str] = []
    current: ast.AST = node
    while isinstance(current, ast.Subscript):
        key = _subscript_key(current)
        if key is None:
            return None
        keys.append(key)
        current = current.value
    if isinstance(current, ast.Name):
        return (current.id, *reversed(keys))
    return None


def _is_read_parameter(node: ast.AST, parameter: str) -> bool:
    if not isinstance(node, ast.Call) or node.keywords or len(node.args) != 2:
        return False
    if not isinstance(node.func, ast.Name) or node.func.id != "read":
        return False
    if not isinstance(node.args[0], ast.Name) or node.args[0].id != "c":
        return False
    return _subscript_path(node.args[1]) == (
        "parent",
        "actual_core_signed_direction",
        "parameters",
        parameter,
    )


def _is_name(node: ast.AST | None, name: str) -> bool:
    return isinstance(node, ast.Name) and node.id == name


def _function(tree: ast.AST, name: str) -> tuple[ast.AST | None, int]:
    nodes = [
        node for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == name
    ]
    return (nodes[0], len(nodes)) if len(nodes) == 1 else (None, len(nodes))


def _assignment_read(function: ast.AST, target: str, parameter: str) -> bool:
    for node in ast.walk(function):
        if isinstance(node, ast.Assign):
            targets = node.targets
        elif isinstance(node, ast.AnnAssign):
            targets = [node.target]
        else:
            continue
        if any(_is_name(item, target) for item in targets):
            if _is_read_parameter(node.value, parameter):
                return True
    return False


def _dict_record_binding(function: ast.AST, key: str, variable: str) -> bool:
    for node in ast.walk(function):
        if isinstance(node, ast.Dict):
            for item_key, item_value in zip(node.keys, node.values):
                if not (isinstance(item_key, ast.Constant)
                        and item_key.value == key):
                    continue
                if (isinstance(item_value, ast.Call)
                        and _is_name(item_value.func, "record")
                        and len(item_value.args) == 1
                        and not item_value.keywords
                        and _is_name(item_value.args[0], variable)):
                    return True
        if isinstance(node, ast.Call) and _is_name(node.func, "dict"):
            for keyword in node.keywords:
                if keyword.arg != key:
                    continue
                value = keyword.value
                if (isinstance(value, ast.Call)
                        and _is_name(value.func, "record")
                        and len(value.args) == 1
                        and not value.keywords
                        and _is_name(value.args[0], variable)):
                    return True
    return False


def _dict_name_binding(function: ast.AST, key: str, variable: str) -> bool:
    for node in ast.walk(function):
        if isinstance(node, ast.Dict):
            for item_key, item_value in zip(node.keys, node.values):
                if (isinstance(item_key, ast.Constant)
                        and item_key.value == key
                        and _is_name(item_value, variable)):
                    return True
        if isinstance(node, ast.Call) and _is_name(node.func, "dict"):
            for keyword in node.keywords:
                if keyword.arg == key and _is_name(keyword.value, variable):
                    return True
    return False


def _build_call(function: ast.AST, name: str) -> bool:
    for node in ast.walk(function):
        if not isinstance(node, ast.Call) or not _is_name(node.func, name):
            continue
        if node.keywords or len(node.args) < 2:
            continue
        if _is_name(node.args[0], "c") and _is_name(node.args[1], "parent"):
            return True
    return False


def _ast_contract(source_text: Any, path: str, kind: str) -> dict[str, Any]:
    result: dict[str, Any] = {
        "path": path,
        "source_executed": False,
        "source_text_supplied": isinstance(source_text, str),
        "sha256": _sha_text(source_text) if isinstance(source_text, str) else None,
        "parsed": False,
        "function_found": False,
        "ast_sha256": None,
        "direct_parameter_read": False,
        "direct_record": False,
        "build_calls_helper": False,
        "build_saves_helper_result": False,
        "errors": [],
    }
    if not isinstance(source_text, str):
        result["errors"].append("source text is not a string")
        return result
    try:
        tree = ast.parse(source_text, filename=path)
    except SyntaxError as exc:
        result["errors"].append("source text is not parseable: " + str(exc))
        return result
    result["parsed"] = True
    function_name = "build" if kind == "q" else "_physical_log_scale"
    function, count = _function(tree, function_name)
    result["function_count"] = count
    if function is None:
        result["errors"].append(function_name + " definition is missing or duplicated")
        return result
    result["function_found"] = True
    result["ast_sha256"] = hashlib.sha256(
        ast.dump(function, include_attributes=False).encode("utf-8")
    ).hexdigest()
    if kind == "q":
        result["direct_parameter_read"] = _assignment_read(
            function, "de", "delta"
        )
        result["direct_record"] = _dict_record_binding(
            function, "actual_delta", "de"
        )
        build, build_count = _function(tree, "build")
        result["build_count"] = build_count
        if build is not None:
            # The selected function is build, so this is intentionally
            # repeated as an explicit producer-boundary check.
            result["direct_parameter_read"] = _assignment_read(
                build, "de", "delta"
            )
            result["direct_record"] = _dict_record_binding(
                build, "actual_delta", "de"
            )
    else:
        result["direct_parameter_read"] = _assignment_read(
            function, "logP", "logPstar"
        )
        result["direct_record"] = _dict_record_binding(
            function, "logPstar", "logP"
        )
        build, build_count = _function(tree, "build")
        result["build_count"] = build_count
        if build is not None:
            result["build_calls_helper"] = _build_call(
                build, "_physical_log_scale"
            )
            result["build_saves_helper_result"] = _dict_name_binding(
                build, "physical_scale", "scale"
            )
    if not result["direct_parameter_read"]:
        result["errors"].append("direct read of the required core parameter is missing")
    if not result["direct_record"]:
        result["errors"].append("direct record of the required actual field is missing")
    if kind == "background":
        if not result["build_calls_helper"]:
            result["errors"].append("build does not call _physical_log_scale(c, parent, ...)")
        if not result["build_saves_helper_result"]:
            result["errors"].append("build does not save physical_scale=scale")
    result["passed"] = not result["errors"]
    return result


def _declared_hashes(container: Any, path: str) -> list[str]:
    """Collect explicit hashes without searching arbitrary source text."""

    if not isinstance(container, Mapping):
        return []
    aliases = (path, path.split("/")[-1])
    found: list[str] = []
    for group_name in ("implementation_hashes", "source_pins", "source_provenance"):
        group = container.get(group_name)
        if not isinstance(group, Mapping):
            continue
        for alias in aliases:
            item = group.get(alias, _MISSING)
            if isinstance(item, str) and len(item) == 64:
                found.append(item.lower())
            elif isinstance(item, Mapping):
                for key in ("sha256", "expected_sha256", "actual_sha256"):
                    value = item.get(key)
                    if isinstance(value, str) and len(value) == 64:
                        found.append(value.lower())
    return list(dict.fromkeys(found))


def _hash_contract(
    source_text: Any,
    path: str,
    containers: Sequence[Any],
) -> dict[str, Any]:
    computed = _sha_text(source_text) if isinstance(source_text, str) else None
    declared: list[str] = []
    for container in containers:
        declared.extend(_declared_hashes(container, path))
    declared = list(dict.fromkeys(declared))
    matches = bool(computed is not None and declared and all(item == computed for item in declared))
    return {
        "path": path,
        "computed_sha256": computed,
        "declared_sha256": declared,
        "declared": bool(declared),
        "matches": matches,
        "missing_hash_evidence": not bool(declared),
        "source_executed": False,
    }


def _family_contract(native: Any, previous: Any) -> dict[str, Any]:
    values: dict[str, Any] = {
        "native": _path(native, ("source_family",)),
        "native_core_packet": _path(
            native, ("actual_core_signed_direction", "source_family")
        ),
        "native_delta_packet": _path(
            native, ("actual_micro_Q_increment_bounds", "source_family")
        ),
        "native_background_packet": _path(
            native, ("actual_micro_background_I1_numeric", "source_family")
        ),
        "previous_parameter_leaf": _path(previous, ("source_family",)),
    }
    missing = [key for key, value in values.items() if value is _MISSING]
    mismatches: list[str] = []
    present = [value for value in values.values() if value is not _MISSING]
    if present:
        first = present[0]
        mismatches = [key for key, value in values.items()
                      if value is not _MISSING and value != first]
    passed = not missing and not mismatches
    return {
        "passed": passed,
        "values": {key: None if value is _MISSING else _copy(value)
                   for key, value in values.items()},
        "missing": missing,
        "mismatches": mismatches,
        "comparison": "exact serialized source_family equality",
    }


def _previous_provenance(value: Any) -> dict[str, Any]:
    expected = {
        "native": (PREVIOUS_NATIVE_PATH, PREVIOUS_NATIVE_SHA256),
        "corrected": (PREVIOUS_CORRECTED_PATH, PREVIOUS_CORRECTED_SHA256),
        "producer": (PREVIOUS_PRODUCER_PATH, PREVIOUS_PRODUCER_SHA256),
    }
    missing: list[str] = []
    mismatches: list[dict[str, Any]] = []
    checks: dict[str, Any] = {}
    if not isinstance(value, Mapping):
        return {
            "present": False,
            "raw": None,
            "passed": False,
            "missing": ["source_provenance"],
            "mismatches": [],
            "checks": checks,
            "source_provenance_scope": "previous receipt sibling, not leaf",
        }
    for name, (path, digest) in expected.items():
        item = value.get(name, _MISSING)
        if not isinstance(item, Mapping):
            missing.append(name)
            checks[name] = {"present": False, "path": path, "sha256": digest}
            continue
        actual_path = item.get("path")
        actual_digest = item.get("sha256")
        passed = actual_path == path and actual_digest == digest
        checks[name] = {
            "present": True,
            "path": actual_path,
            "sha256": actual_digest,
            "expected_path": path,
            "expected_sha256": digest,
            "exact": passed,
        }
        if not passed:
            mismatches.append({
                "name": name,
                "path": actual_path,
                "sha256": actual_digest,
                "expected_path": path,
                "expected_sha256": digest,
            })
    return {
        "present": True,
        "raw": _copy(value),
        "passed": not missing and not mismatches,
        "missing": missing,
        "mismatches": mismatches,
        "checks": checks,
        "source_provenance_scope": "previous receipt sibling, not leaf",
    }


def _record_at(root: Any, path: Sequence[str], label: str) -> dict[str, Any]:
    value = _path(root, path)
    return {
        "path": list(path),
        "label": label,
        "record": _copy(value) if value is not _MISSING else None,
        "present": value is not _MISSING,
        "info": _record_info(value, label),
    }


def _comparison_table(
    native_core_delta: Any,
    native_actual_delta: Any,
    native_core_logpstar: Any,
    native_actual_logpstar: Any,
    previous: Any,
) -> dict[str, Any]:
    previous_delta = _path(previous, ("delta",))
    previous_logpstar = _path(previous, ("actual_logPstar_parameter",))
    previous_delta_native = (
        _path(previous_delta, ("native_record",))
        if previous_delta is not _MISSING else _MISSING
    )
    previous_delta_corrected = (
        _path(previous_delta, ("corrected_record",))
        if previous_delta is not _MISSING else _MISSING
    )
    previous_log_native = (
        _path(previous_logpstar, ("native_record",))
        if previous_logpstar is not _MISSING else _MISSING
    )
    previous_log_corrected = (
        _path(previous_logpstar, ("corrected_record",))
        if previous_logpstar is not _MISSING else _MISSING
    )
    delta_pairs = {
        "native_core_vs_native_actual": (native_core_delta, native_actual_delta),
        "native_core_vs_previous_native": (native_core_delta, previous_delta_native),
        "native_core_vs_previous_corrected": (native_core_delta, previous_delta_corrected),
        "native_actual_vs_previous_native": (native_actual_delta, previous_delta_native),
        "native_actual_vs_previous_corrected": (native_actual_delta, previous_delta_corrected),
        "previous_native_vs_previous_corrected": (
            previous_delta_native, previous_delta_corrected
        ),
    }
    log_pairs = {
        "native_core_vs_native_actual": (native_core_logpstar, native_actual_logpstar),
        "native_core_vs_previous_native": (native_core_logpstar, previous_log_native),
        "native_core_vs_previous_corrected": (native_core_logpstar, previous_log_corrected),
        "native_actual_vs_previous_native": (native_actual_logpstar, previous_log_native),
        "native_actual_vs_previous_corrected": (native_actual_logpstar, previous_log_corrected),
        "previous_native_vs_previous_corrected": (
            previous_log_native, previous_log_corrected
        ),
    }
    delta = {
        key: _compare_records("delta." + key, left, right)
        for key, (left, right) in delta_pairs.items()
    }
    logpstar = {
        key: _compare_records("logPstar." + key, left, right)
        for key, (left, right) in log_pairs.items()
    }
    all_comparisons = [*delta.values(), *logpstar.values()]
    mismatches = [item for item in all_comparisons
                  if item["exact_binary_tuple_equal"] is False]
    return {
        "delta": delta,
        "logPstar": logpstar,
        "all_exact_binary_tuple_equal": not mismatches,
        "mismatch_count": len(mismatches),
        "mismatches": mismatches,
        "overlap_used_as_identity": False,
    }


def verify(
    native_report: Mapping[str, Any],
    previous_parameter_receipt: Mapping[str, Any],
    q_source_text: str,
    background_source_text: str,
) -> dict[str, Any]:
    """Audit direct saved parameter lineage without executing any producer."""

    native, native_prefix, native_errors = _unpack_native(native_report)
    previous, previous_prefix, previous_provenance_value, previous_errors = _unpack_previous(
        previous_parameter_receipt
    )
    errors = [*native_errors, *previous_errors]

    if native is _MISSING:
        core_delta = core_logpstar = actual_delta = actual_logpstar = _MISSING
    else:
        core_parameters = _path(native, CORE_PARAMETER_PATH)
        core_delta = _path(core_parameters, ("delta",))
        core_logpstar = _path(core_parameters, ("logPstar",))
        actual_delta = _path(native, NATIVE_DELTA_PATH)
        actual_logpstar = _path(native, NATIVE_LOGPSTAR_PATH)
        for label, value in (
            ("native core delta", core_delta),
            ("native core logPstar", core_logpstar),
            ("native actual delta", actual_delta),
            ("native actual logPstar", actual_logpstar),
        ):
            if value is _MISSING:
                errors.append(label + " is missing")

    if previous is _MISSING:
        previous_for_comparison: Any = {}
    else:
        previous_for_comparison = previous
        if _path(previous, ("delta", "native_record")) is _MISSING:
            errors.append("previous delta.native_record is missing")
        if _path(previous, ("delta", "corrected_record")) is _MISSING:
            errors.append("previous delta.corrected_record is missing")
        if _path(previous, ("actual_logPstar_parameter", "native_record")) is _MISSING:
            errors.append("previous actual_logPstar_parameter.native_record is missing")
        if _path(previous, ("actual_logPstar_parameter", "corrected_record")) is _MISSING:
            errors.append("previous actual_logPstar_parameter.corrected_record is missing")

    ast_q = _ast_contract(q_source_text, Q_SOURCE_PATH, "q")
    ast_background = _ast_contract(
        background_source_text, BACKGROUND_SOURCE_PATH, "background"
    )
    ast_passed = ast_q["passed"] and ast_background["passed"]

    containers: list[Any] = []
    # Keep the wrapper available when a caller passes a report whose native
    # packet is nested under ``actual_native_directed_controls_v2``.
    containers.append(native_report)
    if native is not _MISSING:
        containers.append(native)
        containers.append(_path(native, NATIVE_DELTA_PATH[:1]))
        containers.append(_path(native, NATIVE_LOGPSTAR_PATH[:2]))
    if previous is not _MISSING:
        containers.append(previous)
    hash_q = _hash_contract(q_source_text, Q_SOURCE_PATH, containers)
    hash_background = _hash_contract(
        background_source_text, BACKGROUND_SOURCE_PATH, containers
    )
    hash_passed = hash_q["matches"] and hash_background["matches"]
    if not hash_q["matches"]:
        errors.append("Q producer source hash evidence is missing or mismatched")
    if not hash_background["matches"]:
        errors.append("background producer source hash evidence is missing or mismatched")

    family = _family_contract(native, previous)
    if not family["passed"]:
        errors.append("source_family evidence is missing or mismatched")
    provenance = _previous_provenance(previous_provenance_value)
    if not provenance["passed"]:
        errors.append("previous source_provenance is missing")

    comparisons = _comparison_table(
        core_delta,
        actual_delta,
        core_logpstar,
        actual_logpstar,
        previous_for_comparison,
    )
    if not comparisons["all_exact_binary_tuple_equal"]:
        errors.append("one or more directed binary tuples differ")

    structural = {
        "native_scope_path": native_prefix,
        "previous_scope_path": previous_prefix,
        "core_and_actual_in_same_native_packet": native is not _MISSING,
        "parameter_paths": {
            "core_delta": list((*CORE_PARAMETER_PATH, "delta")),
            "core_logPstar": list((*CORE_PARAMETER_PATH, "logPstar")),
            "actual_delta": list(NATIVE_DELTA_PATH),
            "actual_logPstar": list(NATIVE_LOGPSTAR_PATH),
        },
        "previous_leaf_schema": (
            previous.get("schema") if isinstance(previous, Mapping) else None
        ),
        "previous_leaf_schema_expected": PREVIOUS_SCHEMA,
    }
    if structural["previous_leaf_schema"] not in (None, PREVIOUS_SCHEMA):
        errors.append("previous parameter leaf schema changed")
        structural["previous_leaf_schema_match"] = False
    else:
        structural["previous_leaf_schema_match"] = True

    if not errors and ast_passed and hash_passed and family["passed"] \
            and provenance["passed"] and comparisons["all_exact_binary_tuple_equal"]:
        status = "proved"
    elif comparisons["mismatch_count"]:
        status = "mismatch"
    else:
        status = "incomplete"

    return {
        "schema": SCHEMA,
        "passed": status == "proved",
        "status": status,
        "errors": list(dict.fromkeys(errors)),
        "source_execution": False,
        "real_receipt_opened_here": False,
        "producer_imported": False,
        "graph_traversed": False,
        "native_scope": structural,
        "source_family": family,
        "source_hashes": {
            "Q_increment_producer": hash_q,
            "background_scale_producer": hash_background,
        },
        "source_AST": {
            "Q_increment_producer": ast_q,
            "background_scale_producer": ast_background,
        },
        "previous_provenance": provenance,
        "comparisons": comparisons,
        "provenance_scope": {
            "archived_core_parameter_lineage_only": True,
            "same_native_packet_core_to_actual_checked": structural[
                "core_and_actual_in_same_native_packet"
            ],
            "source_family_checked": family["passed"],
            "source_hashes_checked": hash_passed,
            "direct_read_record_lineage_checked": ast_passed,
            "cross_run_object_identity_claimed": False,
            "delta_branch_identity_checked": False,
            "live_endpoint_basis_identity": False,
            "native_corrected_RQy_identity": False,
            "source_function_point_value_identity": False,
            "interval_overlap_promoted_to_identity": False,
            "rounding_difference_promoted_to_identity": False,
        },
    }


consume = verify
bind = verify


__all__ = ["verify", "consume", "bind", "SCHEMA"]
