"""Source-only inspection of saved pulse endpoint rewrites.

The archive and acceptance receipt are authenticated before parsing.  Only the
four explicitly named endpoint roles (``gap.upper``, ``end.lower``,
``main.upper``, and ``gap.lower``) are read.  For each role this consumer
checks only ``R_q_y`` and ``q_plus_q_y`` rows 0 and 1, for sixteen possible
opaque-reference lookups in total.

Protected references are never used as mapping keys or dereferenced.  A
rewrite-map hit is evidence about a saved mapping entry only; it does not
prove endpoint identity, common source functions, derivatives, physical
radius, or numerical values.  Same-seam mapped-target equality is reported
only when both non-protected source keys have mappings.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


SCHEMA = "f_internal_endpoint_saved_rewrite_v1"
DEFAULT_FRAME_ARCHIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source.json.gz"
)
DEFAULT_FRAME_ARCHIVE_SHA256 = (
    "442df70e5ab4051066b241e943c8d0531067d057e774a04dcc341dadbd703ebe"
)
DEFAULT_RECEIPT = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source_check.json"
)
DEFAULT_RECEIPT_SHA256 = (
    "320f4623859c4ad64a4b5c6c4912e43e42651b58d8359b334b96b5b0d6d38d15"
)

ORDERS = (0, 1)
FIELDS = ("R_q_y", "q_plus_q_y")
ENDPOINT_ROLES = (
    ("main.upper", "gap.lower"),
    ("gap.lower", "main.upper"),
    ("gap.upper", "end.lower"),
    ("end.lower", "gap.upper"),
)
SEAM_PAIRS = (("main.upper", "gap.lower"), ("gap.upper", "end.lower"))

# Protected opaque references may be compared or reported, but never used as
# a rewrite-map key or graph index.
PROTECTED_SOURCE_REFS = frozenset(
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
        33982,
    }
)


class CaptureInputError(ValueError):
    """A requested saved input was not authenticated or well formed."""


class _UnresolvedCapture(Exception):
    """Internal fail-closed signal for a missing direct source role."""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CaptureInputError(message)


def _path(root: Path, relative: str | Path) -> Path:
    path = Path(relative)
    return path if path.is_absolute() else root / path


def _raw_sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _read_json_after_raw_auth(
    root: Path, relative: str | Path, expected_sha256: str
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = _path(root, relative)
    _require(path.is_file(), f"missing saved input: {Path(relative).as_posix()}")
    raw = path.read_bytes()
    actual = _raw_sha256(raw)
    _require(actual == expected_sha256, f"raw SHA changed: {Path(relative).as_posix()}")
    value = json.loads(raw.decode("utf-8"))
    _require(isinstance(value, Mapping), f"saved input is not an object: {Path(relative).as_posix()}")
    return value, {
        "path": Path(relative).as_posix(),
        "raw_sha256": actual,
        "raw_hash_authenticated_before_parse": True,
        "gzip_decoded": False,
    }


def _read_gzip_json_after_raw_auth(
    root: Path, relative: str | Path, expected_sha256: str
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = _path(root, relative)
    _require(path.is_file(), f"missing raw gzip input: {Path(relative).as_posix()}")
    raw = path.read_bytes()
    actual = _raw_sha256(raw)
    _require(actual == expected_sha256, f"raw gzip SHA changed: {Path(relative).as_posix()}")
    decoded = gzip.decompress(raw)
    value = json.loads(decoded.decode("utf-8"))
    _require(isinstance(value, Mapping), f"decoded gzip input is not an object: {Path(relative).as_posix()}")
    return value, {
        "path": Path(relative).as_posix(),
        "raw_sha256": actual,
        "raw_hash_authenticated_before_parse": True,
        "gzip_decoded": True,
    }


def _direct(value: Any, path: Sequence[str], label: str) -> Any:
    """Read a literal mapping path without interpreting graph nodes."""

    current = value
    for key in path:
        if not isinstance(current, Mapping) or key not in current:
            raise _UnresolvedCapture(f"missing_direct_field:{label}")
        current = current[key]
    return current


def _family_evidence(
    archive: Mapping[str, Any],
    receipt: Mapping[str, Any],
    frame_archive_path: str | Path,
    frame_archive_sha256: str,
) -> dict[str, Any]:
    archive_family = archive.get("source_family")
    receipt_family = receipt.get("source_family")
    archive_visible = isinstance(archive_family, Mapping) and bool(archive_family)
    receipt_visible = isinstance(receipt_family, Mapping) and bool(receipt_family)
    families_match = archive_visible and receipt_visible and archive_family == receipt_family
    input_hashes = receipt.get("input_hashes")
    archive_key = Path(frame_archive_path).name
    archive_pin = input_hashes.get(archive_key) if isinstance(input_hashes, Mapping) else None
    archive_pin_matches = archive_pin == frame_archive_sha256
    all_passed = receipt.get("all_passed") is True
    supported = all_passed and families_match and archive_pin_matches
    return {
        "status": "proved_receipt_acceptance" if supported else "unresolved_receipt_acceptance",
        "supported": supported,
        "all_passed": all_passed,
        "archive_source_family_visible": archive_visible,
        "receipt_source_family_visible": receipt_visible,
        "source_family_matches": families_match,
        "archive_pin_visible": archive_pin is not None,
        "archive_pin_matches": archive_pin_matches,
        "archive_source_family": dict(archive_family) if archive_visible else None,
        "receipt_source_family": dict(receipt_family) if receipt_visible else None,
        "receipt_fields_used": ["all_passed", "source_family", "input_hashes"],
    }


def _base_flags() -> dict[str, Any]:
    return {
        "source_value_selected": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "numeric_values_computed": False,
        "numeric_difference_evaluated": False,
        "numeric_error_computed": False,
        "ordinary_Z3_semantics_proved": False,
        "m3_ordinary_Z3_semantics_proved": False,
        "F_source_identity_proved": False,
        "F_chart_finite_endpoints_proved": False,
        "endpoint_identity_proved": False,
        "same_source_proved": False,
        "source_function_identity_proved": False,
        "rewrite_map_coverage_proved": False,
        "physical_radius_proved": False,
        "global_closure_proved": False,
        "full_goal_complete": False,
        "error_removed": False,
        "old_error_removed": False,
        "old_errors_removed": False,
        "ledger_212_213_183_changed": False,
    }


def _empty_position(
    source_role: str,
    target_role: str | None,
    field: str,
    order: int,
    reason: str,
) -> dict[str, Any]:
    return {
        "status": "unresolved",
        "source_role": source_role,
        "target_role": target_role,
        "field": field,
        "ordinary_Z_order": order,
        "reason": reason,
        "source_value_selected": False,
        "node_id_is_not_a_numeric_value": True,
        "rewrite_map_lookup_attempted": False,
        "same_source_proved": False,
    }


def _empty_checks(reason: str) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    for source_role, target_role in ENDPOINT_ROLES:
        role: dict[str, Any] = {
            "source_role": source_role,
            "target_role": target_role,
            "fields": {},
            "source_values_selected": False,
        }
        for field in FIELDS:
            role["fields"][field] = {
                str(order): _empty_position(
                    source_role,
                    target_role,
                    field,
                    order,
                    reason,
                )
                for order in ORDERS
            }
        checks[source_role] = role
    return checks


def _is_integer_ref(value: Any) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def _inspect_ref(
    ref: Any,
    *,
    source_role: str,
    target_role: str | None,
    field: str,
    order: int,
    rewritten_nodes: Mapping[str, Any],
) -> dict[str, Any]:
    """Check the source ref before attempting the string-key map lookup."""

    position: dict[str, Any] = {
        "source_role": source_role,
        "target_role": target_role,
        "field": field,
        "ordinary_Z_order": order,
        "source_value_selected": False,
        "node_id_is_not_a_numeric_value": True,
        "same_source_proved": False,
    }
    if not _is_integer_ref(ref):
        position.update(
            {
                "status": "unresolved",
                "reason": "saved endpoint row is not an opaque integer reference",
                "rewrite_map_lookup_attempted": False,
            }
        )
        return position
    position["source_node_ref"] = ref
    position["source_ref_protected"] = ref in PROTECTED_SOURCE_REFS
    if ref in PROTECTED_SOURCE_REFS:
        position.update(
            {
                "status": "protected",
                "reason": "protected source reference is not used as a map or graph index",
                "rewrite_map_lookup_attempted": False,
            }
        )
        return position

    map_key = str(ref)
    position["rewrite_map_key"] = map_key
    position["rewrite_map_lookup_attempted"] = True
    if map_key not in rewritten_nodes:
        position.update(
            {
                "status": "not_recorded",
                "reason": "old expression reference is absent from rewritten_nodes",
            }
        )
        return position

    rewritten_ref = rewritten_nodes[map_key]
    position["rewrite_map_recorded"] = True
    if not _is_integer_ref(rewritten_ref):
        position.update(
            {
                "status": "unresolved",
                "reason": "rewritten_nodes target is not an opaque integer reference",
            }
        )
        return position
    position["rewritten_source_ref"] = rewritten_ref
    position["rewritten_ref_protected"] = rewritten_ref in PROTECTED_SOURCE_REFS
    if rewritten_ref in PROTECTED_SOURCE_REFS:
        position.update(
            {
                "status": "protected",
                "reason": "rewritten target is protected and is not dereferenced",
            }
        )
        return position
    position.update(
        {
            "status": "hit",
            "reason": "rewrite map entry observed; target remains opaque",
            "rewrite_map_hit": True,
        }
    )
    return position


def _inspect_role(
    endpoints: Mapping[str, Any],
    rewritten_nodes: Mapping[str, Any],
    *,
    source_role: str,
    target_role: str | None,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "source_role": source_role,
        "fields": {},
        "source_values_selected": False,
    }
    try:
        chart, edge = source_role.split(".")
        chart_endpoints = _direct(endpoints, (chart,), f"endpoints.{chart}")
        if not isinstance(chart_endpoints, list) or len(chart_endpoints) != 2:
            raise _UnresolvedCapture(f"endpoints.{chart}:expected_lower_upper_pair")
        if edge not in ("lower", "upper"):
            raise _UnresolvedCapture(f"endpoints.{source_role}:unknown_edge")
        role_object = chart_endpoints[0 if edge == "lower" else 1]
    except _UnresolvedCapture as exc:
        for field in FIELDS:
            result["fields"][field] = {
                str(order): _empty_position(source_role, target_role, field, order, exc.reason)
                for order in ORDERS
            }
        return result
    if not isinstance(role_object, Mapping):
        for field in FIELDS:
            result["fields"][field] = {
                str(order): _empty_position(
                    source_role,
                    target_role,
                    field,
                    order,
                    "endpoint role is not an object",
                )
                for order in ORDERS
            }
        return result
    for field in FIELDS:
        try:
            raw_rows = _direct(
                role_object,
                (field,),
                f"endpoints.{source_role}.{field}",
            )
        except _UnresolvedCapture as exc:
            raw_rows = None
            missing_reason = exc.reason
        else:
            missing_reason = "endpoint field does not contain ordinary rows 0 and 1"
        if not isinstance(raw_rows, list) or len(raw_rows) < len(ORDERS):
            result["fields"][field] = {
                str(order): _empty_position(
                    source_role,
                    target_role,
                    field,
                    order,
                    missing_reason,
                )
                for order in ORDERS
            }
            continue
        result["fields"][field] = {
            str(order): _inspect_ref(
                raw_rows[order],
                source_role=source_role,
                target_role=target_role,
                field=field,
                order=order,
                rewritten_nodes=rewritten_nodes,
            )
            for order in ORDERS
        }
    return result


def _seam_comparisons(endpoint_checks: Mapping[str, Any]) -> dict[str, Any]:
    comparisons: dict[str, Any] = {}
    for left_role, right_role in SEAM_PAIRS:
        seam_name = left_role + "__" + right_role
        seam_fields: dict[str, Any] = {}
        left = endpoint_checks.get(left_role, {})
        right = endpoint_checks.get(right_role, {})
        for field in FIELDS:
            order_results: dict[str, Any] = {}
            left_rows = left.get("fields", {}).get(field, {}) if isinstance(left, Mapping) else {}
            right_rows = right.get("fields", {}).get(field, {}) if isinstance(right, Mapping) else {}
            for order in ORDERS:
                lrow = left_rows.get(str(order), {}) if isinstance(left_rows, Mapping) else {}
                rrow = right_rows.get(str(order), {}) if isinstance(right_rows, Mapping) else {}
                result: dict[str, Any] = {
                    "left_role": left_role,
                    "right_role": right_role,
                    "field": field,
                    "ordinary_Z_order": order,
                    "status": "unresolved",
                    "same_source_proved": False,
                    "source_value_selected": False,
                    "equality_attempted": False,
                }
                if isinstance(lrow, Mapping) and isinstance(rrow, Mapping):
                    lkey_ok = lrow.get("status") == "hit" and lrow.get("source_ref_protected") is False
                    rkey_ok = rrow.get("status") == "hit" and rrow.get("source_ref_protected") is False
                    if lkey_ok and rkey_ok:
                        lref = lrow.get("rewritten_source_ref")
                        rref = rrow.get("rewritten_source_ref")
                        result["equality_attempted"] = True
                        result["left_rewritten_source_ref"] = lref
                        result["right_rewritten_source_ref"] = rref
                        result["status"] = "equal" if lref == rref else "different"
                order_results[str(order)] = result
            seam_fields[field] = order_results
        comparisons[seam_name] = {
            "left_role": left_role,
            "right_role": right_role,
            "fields": seam_fields,
            "source_values_selected": False,
            "endpoint_identity_proved": False,
        }
    return comparisons


def _unresolved_result(
    *,
    reason: str,
    inputs: Mapping[str, Any] | None = None,
    receipt: Mapping[str, Any] | None = None,
    pulse5_archive_name: str | None = None,
) -> dict[str, Any]:
    checks = _empty_checks(reason)
    return {
        "schema": SCHEMA,
        "status": "unresolved_endpoint_rewrite",
        "reason": reason,
        "pulse5_archive_name": pulse5_archive_name,
        "admissible_for_numeric_use": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "inputs": dict(inputs or {}),
        "receipt_evidence": dict(receipt or {}),
        "endpoint_checks": checks,
        "seam_comparisons": _seam_comparisons(checks),
        "rewrite_map_coverage_observed": False,
        "flags": _base_flags(),
        "scope_audit": {
            "entire_archive_json_decoded": bool(
                isinstance(inputs, Mapping)
                and inputs.get("frame_archive_decoded_only_after_raw_hash") is True
            ),
            "exact_graph_nodes_explicitly_selected": False,
            "exact_graph_nodes_accessed": False,
            "graph_accessed": False,
            "protected_references_used_only_for_equality": True,
            "protected_references_used_as_graph_or_mapping_indices": False,
            "source_values_selected": False,
            "numeric_values_computed": False,
        },
    }


def recover(
    root: str | Path,
    *,
    frame_archive_path: str | Path = DEFAULT_FRAME_ARCHIVE,
    frame_archive_sha256: str = DEFAULT_FRAME_ARCHIVE_SHA256,
    receipt_path: str | Path = DEFAULT_RECEIPT,
    receipt_sha256: str = DEFAULT_RECEIPT_SHA256,
    pulse5_archive_name: str | None = None,
) -> dict[str, Any]:
    """Inspect four endpoint roles and their sixteen saved row mappings."""

    repo_root = Path(root)
    inputs: dict[str, Any] = {
        "frame_archive_path": Path(frame_archive_path).as_posix(),
        "acceptance_receipt_path": Path(receipt_path).as_posix(),
        "frame_archive_decoded_only_after_raw_hash": False,
        "entire_archive_json_decoded": False,
        "pulse5_archive_name": pulse5_archive_name,
    }
    try:
        archive, archive_meta = _read_gzip_json_after_raw_auth(
            repo_root, frame_archive_path, frame_archive_sha256
        )
        receipt, receipt_meta = _read_json_after_raw_auth(
            repo_root, receipt_path, receipt_sha256
        )
    except (CaptureInputError, OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        return _unresolved_result(reason=str(exc), inputs=inputs, pulse5_archive_name=pulse5_archive_name)

    inputs.update(
        {
            "frame_archive": archive_meta,
            "acceptance_receipt": receipt_meta,
            "frame_archive_decoded_only_after_raw_hash": True,
            "entire_archive_json_decoded": True,
        }
    )
    family = _family_evidence(archive, receipt, frame_archive_path, frame_archive_sha256)
    if not family["supported"]:
        return _unresolved_result(
            reason="archive_receipt_source_family_or_acceptance_unresolved",
            inputs=inputs,
            receipt=family,
            pulse5_archive_name=pulse5_archive_name,
        )
    if not isinstance(pulse5_archive_name, str) or not pulse5_archive_name:
        return _unresolved_result(
            reason="pulse5_archive_name_missing_or_uncertain",
            inputs=inputs,
            receipt=family,
            pulse5_archive_name=pulse5_archive_name,
        )
    name_path = Path(pulse5_archive_name)
    if name_path.name != pulse5_archive_name or ".." in name_path.parts:
        return _unresolved_result(
            reason="pulse5_archive_name_must_be_a_direct_archive_key",
            inputs=inputs,
            receipt=family,
            pulse5_archive_name=pulse5_archive_name,
        )

    try:
        bindings = _direct(archive, ("Rp_native_C5_pulse_source_bindings",), "archive.Rp_native_C5_pulse_source_bindings")
        archive_inputs = _direct(bindings, ("inputs",), "snapshot.inputs")
        pulse5_snapshot = _direct(archive_inputs, (pulse5_archive_name,), f"snapshot.inputs[{pulse5_archive_name}]")
        endpoints = _direct(pulse5_snapshot, ("endpoints",), f"snapshot.inputs[{pulse5_archive_name}].endpoints")
        rewritten_nodes = _direct(bindings, ("rewritten_nodes",), "snapshot.rewritten_nodes")
        if not isinstance(endpoints, Mapping):
            raise _UnresolvedCapture("snapshot pulse5 endpoints are not an object")
        if not isinstance(rewritten_nodes, Mapping):
            raise _UnresolvedCapture("snapshot.rewritten_nodes is not an object")
        endpoint_checks = {
            source_role: _inspect_role(
                endpoints,
                rewritten_nodes,
                source_role=source_role,
                target_role=target_role,
            )
            for source_role, target_role in ENDPOINT_ROLES
        }
    except _UnresolvedCapture as exc:
        return _unresolved_result(
            reason=exc.reason,
            inputs=inputs,
            receipt=family,
            pulse5_archive_name=pulse5_archive_name,
        )

    positions = [
        position
        for role in endpoint_checks.values()
        for field in role["fields"].values()
        for position in field.values()
    ]
    hits = sum(position.get("status") == "hit" for position in positions)
    protected = sum(position.get("status") == "protected" for position in positions)
    not_recorded = sum(position.get("status") == "not_recorded" for position in positions)
    seam_comparisons = _seam_comparisons(endpoint_checks)
    equalities = sum(
        item.get("status") == "equal"
        for seam in seam_comparisons.values()
        for field in seam["fields"].values()
        for item in field.values()
    )
    differences = sum(
        item.get("status") == "different"
        for seam in seam_comparisons.values()
        for field in seam["fields"].values()
        for item in field.values()
    )
    return {
        "schema": SCHEMA,
        "status": "returned",
        "pulse5_archive_name": pulse5_archive_name,
        "admissible_for_numeric_use": False,
        "source_values_selected": False,
        "source_values_computed": False,
        "inputs": inputs,
        "receipt_evidence": family,
        "endpoint_checks": endpoint_checks,
        "seam_comparisons": seam_comparisons,
        "rewrite_map_coverage_observed": hits == len(positions),
        "rewrite_map_hit_count": hits,
        "rewrite_map_not_recorded_count": not_recorded,
        "protected_unresolved_count": protected,
        "seam_equal_count": equalities,
        "seam_different_count": differences,
        "flags": _base_flags(),
        "scope_audit": {
            "selected_direct_source_roles": [
                f"snapshot.inputs[{pulse5_archive_name}].endpoints.{role}"
                for role, _target_role in ENDPOINT_ROLES
            ],
            "selected_fields": list(FIELDS),
            "selected_ordinary_Z_orders": list(ORDERS),
            "possible_rewrite_map_lookups": len(ENDPOINT_ROLES) * len(FIELDS) * len(ORDERS),
            "rewritten_nodes_map_read": True,
            "entire_archive_json_decoded": True,
            "exact_graph_nodes_explicitly_selected": False,
            "exact_graph_nodes_accessed": False,
            "graph_accessed": False,
            "protected_references_used_only_for_equality": True,
            "protected_references_used_as_graph_or_mapping_indices": False,
            "source_values_selected": False,
            "numeric_values_computed": False,
            "derivatives_evaluated": False,
            "endpoint_identity_proved": False,
            "same_source_proved": False,
            "error_terms_removed": False,
        },
    }


build_f_internal_endpoint_saved_rewrite = recover


__all__ = [
    "CaptureInputError",
    "DEFAULT_FRAME_ARCHIVE",
    "DEFAULT_FRAME_ARCHIVE_SHA256",
    "DEFAULT_RECEIPT",
    "DEFAULT_RECEIPT_SHA256",
    "PROTECTED_SOURCE_REFS",
    "SCHEMA",
    "build_f_internal_endpoint_saved_rewrite",
    "recover",
]
