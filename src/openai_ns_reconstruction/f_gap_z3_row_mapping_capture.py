"""Bounded, source-only capture check for the fourth F-source jet row.

This consumer is deliberately narrower than the older frame proof.  It reads
one saved proof result, one authenticated Rp archive, and its authenticated
receipt.  It compares the proof's saved role rows with the corresponding
snapshot rows and then checks only row index ``3`` (the fourth saved row):

* ``m1`` uses ``input_substitution``;
* ``ap`` uses ``rewritten_nodes``.

Node identifiers are opaque references.  The authenticated gzip is decoded as
the complete archive JSON, but this module never selects or indexes the
``exact_graph_nodes`` field, follows children, evaluates a node, imports an
older source helper, or computes a numerical or mathematical quantity.  In
particular, an exact index-3 mapping is not promoted to an ordinary-Z
derivative claim.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "f_gap_z3_row_mapping_capture_v1"
ROW_INDEX = 3

DEFAULT_SNAPSHOT = Path("tmp/frame_main_gap_source_proof_v1/result.json")
DEFAULT_SNAPSHOT_SHA256 = (
    "af33465e04cd8901060227f4ab75f072efdd5c4d2fd59893a7ff9add27404144"
)
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
DEFAULT_SOURCE_INPUT_KEY = (
    "lei_ren_part1_paper_compliant_current_original_"
    "n1_pulse_pressure_C3_functions.json.gz"
)

# These identifiers are protected by the caller's boundary.  They may only
# be compared as opaque references.  A nonconstant protected identifier is a
# fail-closed unresolved case; it is never used as a mapping or graph key.
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
OPAQUE_CONSTANT_REFS = frozenset({0, 1})


class CaptureInputError(ValueError):
    """A requested input could not be authenticated or parsed."""


class _UnresolvedCapture(Exception):
    """Internal fail-closed signal for absent or protected evidence."""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CaptureInputError(message)


def _raw_sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _path(root: Path, relative: str | Path) -> Path:
    return Path(relative) if Path(relative).is_absolute() else root / Path(relative)


def _read_json_after_raw_auth(
    root: Path, relative: str | Path, expected_sha256: str
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    path = _path(root, relative)
    _require(path.is_file(), f"missing saved input: {Path(relative).as_posix()}")
    raw = path.read_bytes()
    actual = _raw_sha256(raw)
    _require(
        actual == expected_sha256,
        f"raw SHA changed: {Path(relative).as_posix()}",
    )
    value = json.loads(raw.decode("utf-8"))
    _require(
        isinstance(value, Mapping),
        f"saved input is not an object: {Path(relative).as_posix()}",
    )
    return value, {
        "path": Path(relative).as_posix(),
        "raw_sha256": actual,
        "raw_hash_authenticated_before_parse": True,
        "gzip_decoded": False,
    }


def _read_gzip_json_after_raw_auth(
    root: Path, relative: str | Path, expected_sha256: str
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    """Authenticate the complete gzip byte stream before decompression."""

    path = _path(root, relative)
    _require(path.is_file(), f"missing raw gzip input: {Path(relative).as_posix()}")
    raw = path.read_bytes()
    actual = _raw_sha256(raw)
    _require(
        actual == expected_sha256,
        f"raw gzip SHA changed: {Path(relative).as_posix()}",
    )
    decoded = gzip.decompress(raw)
    value = json.loads(decoded.decode("utf-8"))
    _require(
        isinstance(value, Mapping),
        f"decoded gzip input is not an object: {Path(relative).as_posix()}",
    )
    return value, {
        "path": Path(relative).as_posix(),
        "raw_sha256": actual,
        "raw_hash_authenticated_before_parse": True,
        "gzip_decoded": True,
    }


def _direct(value: Any, path: Sequence[str], label: str) -> Any:
    current = value
    for key in path:
        if not isinstance(current, Mapping) or key not in current:
            raise _UnresolvedCapture(f"missing direct field: {label}")
        current = current[key]
    return current


def _as_opaque_ref(value: Any, label: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise _UnresolvedCapture(f"{label}:source_ref_is_not_an_opaque_integer")
    return value


def _guard_ref(value: Any, label: str) -> int:
    """Validate a direct reference without treating it as a numeric value."""

    ref = _as_opaque_ref(value, label)
    if ref in PROTECTED_SOURCE_REFS and ref not in OPAQUE_CONSTANT_REFS:
        raise _UnresolvedCapture(f"{label}:protected_source_ref_unresolved")
    return ref


def _guard_index3(value: Any, label: str) -> int:
    if not isinstance(value, list) or len(value) <= ROW_INDEX:
        raise _UnresolvedCapture(f"{label}:fewer_than_four_saved_rows")
    # Deliberately do not inspect any row other than index 3.  Other rows may
    # contain protected references that are outside this bounded proof.
    return _guard_ref(value[ROW_INDEX], f"{label}[{ROW_INDEX}]")


def _prior_payload(
    snapshot: Mapping[str, Any],
    frame_relative: str | Path,
    frame_sha256: str,
) -> tuple[Mapping[str, Any], dict[str, Any]]:
    """Validate the saved proof wrapper before reading its source roles."""

    if snapshot.get("status") != "returned":
        raise _UnresolvedCapture("prior_result:wrapper_status_not_returned")
    nested = snapshot.get("result")
    if not isinstance(nested, Mapping):
        raise _UnresolvedCapture("prior_result:wrapper_result_missing")
    if nested.get("schema") != "frame_main_gap_source_proof_v1":
        raise _UnresolvedCapture("prior_result:payload_schema_mismatch")
    if nested.get("status") != "conditional_algebra_bound_to_saved_source_roles":
        raise _UnresolvedCapture("prior_result:payload_status_mismatch")
    source_raw_hashes = nested.get("source_raw_hashes")
    if not isinstance(source_raw_hashes, Mapping):
        raise _UnresolvedCapture("prior_result:source_raw_hashes_missing")
    archive_key = Path(frame_relative).name
    if source_raw_hashes.get(archive_key) != frame_sha256:
        raise _UnresolvedCapture("prior_result:RP_archive_hash_mismatch")
    return nested, {
        "status": "proved_prior_wrapper_and_RP_hash",
        "supported": True,
        "wrapper_status": snapshot.get("status"),
        "payload_schema": nested.get("schema"),
        "payload_status": nested.get("status"),
        "RP_archive_key": archive_key,
        "RP_archive_hash": source_raw_hashes.get(archive_key),
    }


def _prior_role_ref(
    payload: Mapping[str, Any], role: str, label: str
) -> int:
    roles = payload.get("source_roles")
    if not isinstance(roles, Mapping):
        raise _UnresolvedCapture("prior_result:source_roles_missing")
    return _guard_index3(roles.get(role), f"prior.source_roles.{label}")


def _compare_prior_rows(
    prior: int, current: int, label: str
) -> dict[str, Any]:
    """Compare only the bounded opaque index-3 row before rewrite lookup."""

    if prior != current:
        return {
            "status": "unresolved_prior_role_mismatch",
            "supported": False,
            "label": label,
            "index": ROW_INDEX,
        }
    return {
        "status": "proved_prior_role_matches_snapshot",
        "supported": True,
        "label": label,
        "index": ROW_INDEX,
    }


def _receipt_evidence(
    receipt: Mapping[str, Any],
    frame_relative: str | Path,
    frame_sha256: str,
    frame_source_family: Any,
) -> dict[str, Any]:
    """Use only visible receipt fields for source-family and acceptance."""

    source_family = receipt.get("source_family")
    accepted = receipt.get("all_passed") is True
    family_visible = isinstance(source_family, Mapping) and bool(source_family)
    frame_family_visible = isinstance(frame_source_family, Mapping) and bool(
        frame_source_family
    )
    family_matches_frame = (
        family_visible
        and frame_family_visible
        and source_family == frame_source_family
    )

    input_hashes = receipt.get("input_hashes")
    archive_pin = None
    archive_pin_matches = False
    if isinstance(input_hashes, Mapping):
        archive_pin = input_hashes.get(Path(frame_relative).name)
        archive_pin_matches = archive_pin == frame_sha256

    supported = (
        accepted
        and family_visible
        and frame_family_visible
        and family_matches_frame
        and archive_pin_matches
    )
    return {
        "status": "proved_receipt_acceptance" if supported else "unresolved_receipt_acceptance",
        "supported": supported,
        "all_passed": accepted,
        "source_family_visible": family_visible,
        "source_family_kind": "mapping" if family_visible else None,
        "source_family": dict(source_family) if family_visible else None,
        "frame_source_family_visible": frame_family_visible,
        "frame_source_family_kind": "mapping" if frame_family_visible else None,
        "source_family_matches_frame": family_matches_frame,
        "archive_pin_visible": archive_pin is not None,
        "archive_pin_matches": archive_pin_matches,
        "receipt_fields_used": ["all_passed", "source_family", "input_hashes"],
    }


def _mapping_evidence(
    *,
    kind: str,
    old_ref: int,
    current_ref: int,
    mapping: Any,
) -> dict[str, Any]:
    """Check one opaque row without dereferencing a graph node."""

    old = old_ref
    current = current_ref
    result: dict[str, Any] = {
        "kind": kind,
        "index": ROW_INDEX,
        "old_source_ref": old,
        "current_source_ref": current,
        "source_values_selected": False,
        "graph_accessed": False,
        "graph_accessed_definition": (
            "no graph node selected/indexed/traversed; complete archive JSON parsing is not graph access"
        ),
        "exact_mapping_supported": False,
    }

    # 0 and 1 are opaque global source references, never mapping keys.
    if old in OPAQUE_CONSTANT_REFS:
        result["status"] = (
            "proved_opaque_constant_reference_equal"
            if old == current
            else "unresolved_opaque_constant_reference_changed"
        )
        result["exact_mapping_supported"] = old == current
        return result

    # _guard_index3 already rejects these, but keep this local guard next to the
    # only index-3 mapping boundary so a future caller cannot bypass it.
    if old in PROTECTED_SOURCE_REFS or current in PROTECTED_SOURCE_REFS:
        result["status"] = "unresolved_protected_source_ref"
        return result
    if not isinstance(mapping, Mapping):
        result["status"] = "unresolved_mapping_owner_missing"
        return result

    key = str(old)
    if key not in mapping:
        result["status"] = "unresolved_exact_mapping_missing"
        return result
    mapped = mapping[key]
    try:
        mapped_ref = _guard_ref(mapped, f"{kind}.mapping[{key}]")
    except _UnresolvedCapture as exc:
        result["status"] = exc.reason
        return result
    result["mapped_source_ref"] = mapped_ref
    result["status"] = (
        "proved_exact_mapping"
        if mapped_ref == current
        else "unresolved_exact_mapping_target_mismatch"
    )
    result["exact_mapping_supported"] = mapped_ref == current
    return result


def _unresolved_result(
    *,
    reason: str,
    inputs: Mapping[str, Any] | None = None,
    receipt: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": "unresolved_bounded_capture_evidence",
        "reason": reason,
        "admissible_for_numeric_use": False,
        "exact_mapping_supported": False,
        "index": ROW_INDEX,
        "inputs": dict(inputs or {}),
        "receipt_evidence": dict(receipt or {}),
        "flags": {
            "m3_ordinary_Z3_semantics_proved": False,
            "ordinary_Z3_semantics_proved": False,
            "ordinary_Z_transport_proved": False,
            "F_source_identity_proved": False,
            "physical_radius_proved": False,
            "numeric_values_computed": False,
            "numeric_difference_evaluated": False,
            "numeric_error_computed": False,
            "global_closure_proved": False,
            "full_goal_complete": False,
            "source_values_selected": False,
            "graph_accessed": False,
            "graph_accessed_definition": (
                "no graph node selected/indexed/traversed; complete archive JSON parsing is not graph access"
            ),
            "old_helper_imported": False,
            "old_constructor_called": False,
        },
    }


def recover(
    root: str | Path,
    *,
    snapshot_path: str | Path = DEFAULT_SNAPSHOT,
    snapshot_sha256: str = DEFAULT_SNAPSHOT_SHA256,
    frame_archive_path: str | Path = DEFAULT_FRAME_ARCHIVE,
    frame_archive_sha256: str = DEFAULT_FRAME_ARCHIVE_SHA256,
    receipt_path: str | Path = DEFAULT_RECEIPT,
    receipt_sha256: str = DEFAULT_RECEIPT_SHA256,
    source_input_key: str = DEFAULT_SOURCE_INPUT_KEY,
) -> dict[str, Any]:
    """Return the bounded index-3 capture result.

    The three byte streams are authenticated before parsing; gzip bytes are
    decompressed only after their complete raw SHA-256 matches.  Missing
    source roles, protected references, or absent rewrite evidence remain
    unresolved rather than being inferred.
    """

    repo_root = Path(root)
    snapshot, snapshot_receipt = _read_json_after_raw_auth(
        repo_root, snapshot_path, snapshot_sha256
    )
    frame_archive, frame_receipt = _read_gzip_json_after_raw_auth(
        repo_root, frame_archive_path, frame_archive_sha256
    )
    acceptance_receipt, acceptance_receipt_meta = _read_json_after_raw_auth(
        repo_root, receipt_path, receipt_sha256
    )
    input_meta = {
        "snapshot": snapshot_receipt,
        "frame_archive": frame_receipt,
        "acceptance_receipt": acceptance_receipt_meta,
        "frame_archive_decoded_only_after_raw_hash": True,
    }

    try:
        prior, prior_wrapper_evidence = _prior_payload(
            snapshot, frame_archive_path, frame_archive_sha256
        )
        bindings = _direct(
            frame_archive,
            ("Rp_native_C5_pulse_source_bindings",),
            "frame.Rp_native_C5_pulse_source_bindings",
        )
        frame_source_family = frame_archive.get("source_family")
        inputs = _direct(bindings, ("inputs",), "frame.inputs")
        source_input = _direct(inputs, (source_input_key,), f"frame.inputs.{source_input_key}")
        old_native = _guard_index3(
            _direct(source_input, ("native", "m1"), "frame.inputs.source.native.m1"),
            "snapshot.original_native_m1",
        )
        current_native = _guard_index3(
            _direct(bindings, ("current_frame", "m1"), "frame.current_frame.m1"),
            "snapshot.current_native_m1",
        )
        old_ap = _guard_index3(
            _direct(source_input, ("selected_pulse", "ap"), "frame.inputs.source.selected_pulse.ap"),
            "snapshot.original_selected_ap",
        )
        current_ap = _guard_index3(
            _direct(bindings, ("selected_pulse", "ap"), "frame.selected_pulse.ap"),
            "snapshot.current_selected_ap",
        )

        prior_rows = {
            "original_native_m1": _prior_role_ref(
                prior, "original_native_m1", "original_native_m1"
            ),
            "current_native_m1": _prior_role_ref(
                prior, "current_native_m1", "current_native_m1"
            ),
            "original_selected_ap": _prior_role_ref(
                prior, "original_selected_ap", "original_selected_ap"
            ),
            "current_selected_ap": _prior_role_ref(
                prior, "current_selected_ap", "current_selected_ap"
            ),
        }
        prior_checks = {
            "original_native_m1": _compare_prior_rows(
                prior_rows["original_native_m1"], old_native, "original_native_m1"
            ),
            "current_native_m1": _compare_prior_rows(
                prior_rows["current_native_m1"], current_native, "current_native_m1"
            ),
            "original_selected_ap": _compare_prior_rows(
                prior_rows["original_selected_ap"], old_ap, "original_selected_ap"
            ),
            "current_selected_ap": _compare_prior_rows(
                prior_rows["current_selected_ap"], current_ap, "current_selected_ap"
            ),
        }
        prior_rows_match = all(
            check["supported"] for check in prior_checks.values()
        )
        if not prior_rows_match:
            return _unresolved_result(
                reason="prior_saved_role_does_not_match_archive_snapshot",
                inputs=input_meta,
            ) | {
                "prior_wrapper_evidence": prior_wrapper_evidence,
                "prior_role_checks": prior_checks,
            }

        receipt_evidence = _receipt_evidence(
            acceptance_receipt,
            frame_archive_path,
            frame_archive_sha256,
            frame_source_family,
        )
        if not receipt_evidence["supported"]:
            return _unresolved_result(
                reason="visible_receipt_does_not_prove_source_acceptance",
                inputs=input_meta,
                receipt=receipt_evidence,
            ) | {
                "prior_wrapper_evidence": prior_wrapper_evidence,
                "prior_role_checks": prior_checks,
            }

        input_substitution = _direct(
            bindings, ("input_substitution",), "frame.input_substitution"
        )
        rewritten_nodes = _direct(
            bindings, ("rewritten_nodes",), "frame.rewritten_nodes"
        )
        native_mapping = _mapping_evidence(
            kind="m1",
            old_ref=old_native,
            current_ref=current_native,
            mapping=input_substitution,
        )
        ap_mapping = _mapping_evidence(
            kind="ap",
            old_ref=old_ap,
            current_ref=current_ap,
            mapping=rewritten_nodes,
        )
        exact_supported = bool(
            native_mapping["exact_mapping_supported"]
            and ap_mapping["exact_mapping_supported"]
        )
        return {
            "schema": SCHEMA,
            "status": (
                "bounded_index3_exact_mapping_proved"
                if exact_supported
                else "unresolved_bounded_index3_mapping"
            ),
            "admissible_for_numeric_use": False,
            "exact_mapping_supported": exact_supported,
            "index": ROW_INDEX,
            "inputs": input_meta,
            "prior_wrapper_evidence": prior_wrapper_evidence,
            "receipt_evidence": receipt_evidence,
            "prior_role_checks": prior_checks,
            "mapping_evidence": {"m1": native_mapping, "ap": ap_mapping},
            "source_values_selected": False,
            "flags": {
                "m3_ordinary_Z3_semantics_proved": False,
                "ordinary_Z3_semantics_proved": False,
                "ordinary_Z_transport_proved": False,
                "F_source_identity_proved": False,
                "physical_radius_proved": False,
                "numeric_values_computed": False,
                "numeric_difference_evaluated": False,
                "numeric_error_computed": False,
                "global_closure_proved": False,
                "full_goal_complete": False,
                "source_values_selected": False,
                "graph_accessed": False,
                "graph_accessed_definition": (
                    "no graph node selected/indexed/traversed; complete archive JSON parsing is not graph access"
                ),
                "old_helper_imported": False,
                "old_constructor_called": False,
            },
            "scope_audit": {
                "snapshot_roles_read": True,
                "current_frame_index3_read": True,
                "selected_pulse_index3_read": True,
                "other_rows_read": False,
                "input_substitution_index3_only": True,
                "rewritten_nodes_index3_only": True,
                "entire_archive_json_decoded": True,
                "exact_graph_nodes_explicitly_selected": False,
                "exact_graph_nodes_accessed": False,
                "exact_graph_nodes_access_definition": (
                    "false means no exact_graph_nodes node was selected/indexed/traversed; the complete JSON was parsed"
                ),
                "children_traversed": False,
                "source_constructors_called": 0,
                "mathematical_evaluations": 0,
                "numeric_values_computed": 0,
            },
        }
    except _UnresolvedCapture as exc:
        return _unresolved_result(reason=exc.reason, inputs=input_meta)


build_f_gap_z3_row_mapping_capture = recover


__all__ = [
    "CaptureInputError",
    "DEFAULT_FRAME_ARCHIVE",
    "DEFAULT_FRAME_ARCHIVE_SHA256",
    "DEFAULT_RECEIPT",
    "DEFAULT_RECEIPT_SHA256",
    "DEFAULT_SNAPSHOT",
    "DEFAULT_SNAPSHOT_SHA256",
    "DEFAULT_SOURCE_INPUT_KEY",
    "PROTECTED_SOURCE_REFS",
    "SCHEMA",
    "build_f_gap_z3_row_mapping_capture",
    "recover",
]
