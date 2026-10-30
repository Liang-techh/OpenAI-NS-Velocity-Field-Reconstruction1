"""Bounded source-only identity capture for the accepted quiet ``Rp`` owner.

The consumer authenticates the complete Rp-native gzip and its acceptance
receipt before parsing either object.  It then reads only three direct saved
references from ``Rp_native_C5_pulse_source_bindings``:

* ``geometry.Rp``;
* ``inputs[quiet archive name].Rp``; and
* ``inputs[quiet archive name].Rp_endpoint.radius``.

The references are opaque source handles.  Equality, including equality of a
protected handle, is allowed; no handle is used as a graph or mapping index.
The complete archive JSON is parsed honestly, but ``exact_graph_nodes`` is
never selected, indexed, or traversed.  This leaf does not evaluate source
values or prove an F-chart endpoint, a physical radius, ordinary-Z semantics,
closure, or error reduction.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SCHEMA = "f_quiet_rp_owner_link_v1"
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
DEFAULT_QUIET_ARCHIVE_NAME = (
    "lei_ren_part1_paper_compliant_current_original_"
    "true_quiet_Rp_C5_pressure_C3_functions.json.gz"
)

# These are protected opaque source references from the surrounding handoff.
# They may be compared for equality, but none may be used as a graph or map
# index by this leaf.
PROTECTED_SOURCE_REFS = frozenset(
    {
        0,
        1,
        17,
        25220,
        33572,
        205147,
        262918,
        263344,
        263346,
        263348,
        263349,
        33982,
    }
)


class CaptureInputError(ValueError):
    """A requested saved input was not authenticated or well formed."""


class _UnresolvedCapture(Exception):
    """Internal fail-closed signal for missing or mismatched evidence."""

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
    """Authenticate all gzip bytes before decompression and JSON parsing."""

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
    """Read only a named mapping path; never interpret a value as a node id."""

    current = value
    for key in path:
        if not isinstance(current, Mapping) or key not in current:
            raise _UnresolvedCapture(f"missing_direct_field:{label}")
        current = current[key]
    return current


def _opaque_ref(value: Any, label: str) -> tuple[int, bool]:
    if not isinstance(value, int) or isinstance(value, bool):
        raise _UnresolvedCapture(f"{label}:not_an_opaque_integer_reference")
    return value, value in PROTECTED_SOURCE_REFS


def _family_evidence(
    archive: Mapping[str, Any],
    receipt: Mapping[str, Any],
    frame_archive_path: str | Path,
    frame_archive_sha256: str,
) -> dict[str, Any]:
    frame_family = archive.get("source_family")
    receipt_family = receipt.get("source_family")
    frame_visible = isinstance(frame_family, Mapping) and bool(frame_family)
    receipt_visible = isinstance(receipt_family, Mapping) and bool(receipt_family)
    families_match = frame_visible and receipt_visible and frame_family == receipt_family

    input_hashes = receipt.get("input_hashes")
    archive_key = Path(frame_archive_path).name
    archive_pin = input_hashes.get(archive_key) if isinstance(input_hashes, Mapping) else None
    archive_pin_matches = archive_pin == frame_archive_sha256
    accepted = receipt.get("all_passed") is True
    supported = accepted and frame_visible and receipt_visible and families_match and archive_pin_matches
    return {
        "status": "proved_receipt_acceptance" if supported else "unresolved_receipt_acceptance",
        "supported": supported,
        "all_passed": accepted,
        "archive_source_family_visible": frame_visible,
        "receipt_source_family_visible": receipt_visible,
        "source_family_matches": families_match,
        "archive_pin_visible": archive_pin is not None,
        "archive_pin_matches": archive_pin_matches,
        "archive_source_family": dict(frame_family) if frame_visible else None,
        "receipt_source_family": dict(receipt_family) if receipt_visible else None,
        "receipt_fields_used": ["all_passed", "source_family", "input_hashes"],
    }


def _reference_evidence(
    geometry_rp: tuple[int, bool],
    quiet_rp: tuple[int, bool],
    endpoint_radius: tuple[int, bool],
) -> dict[str, Any]:
    geometry_ref, geometry_protected = geometry_rp
    quiet_ref, quiet_protected = quiet_rp
    endpoint_ref, endpoint_protected = endpoint_radius
    geometry_equals_quiet = geometry_ref == quiet_ref
    quiet_equals_endpoint = quiet_ref == endpoint_ref
    geometry_equals_endpoint = geometry_ref == endpoint_ref
    same_reference = geometry_equals_quiet and quiet_equals_endpoint and geometry_equals_endpoint
    return {
        "geometry_Rp": {
            "opaque_reference": geometry_ref,
            "protected": geometry_protected,
        },
        "quiet_input_Rp": {
            "opaque_reference": quiet_ref,
            "protected": quiet_protected,
        },
        "quiet_input_Rp_endpoint_radius": {
            "opaque_reference": endpoint_ref,
            "protected": endpoint_protected,
        },
        "comparisons": {
            "geometry_Rp_equals_quiet_input_Rp": geometry_equals_quiet,
            "quiet_input_Rp_equals_endpoint_radius": quiet_equals_endpoint,
            "geometry_Rp_equals_endpoint_radius": geometry_equals_endpoint,
            "all_three_are_same_opaque_reference": same_reference,
        },
        "protected_reference_equality_only": True,
        "protected_references_used_as_graph_or_mapping_indices": False,
    }


def _base_flags() -> dict[str, Any]:
    return {
        "ordinary_Z3_semantics_proved": False,
        "m3_ordinary_Z3_semantics_proved": False,
        "ordinary_Z_transport_proved": False,
        "F_source_identity_proved": False,
        "F_chart_finite_endpoints_proved": False,
        "physical_radius_proved": False,
        "physical_radius_numeric_proved": False,
        "source_value_selected": False,
        "source_values_selected": False,
        "numeric_values_computed": False,
        "numeric_difference_evaluated": False,
        "numeric_error_computed": False,
        "global_closure_proved": False,
        "full_goal_complete": False,
        "error_removed": False,
        "old_error_removed": False,
    }


def _unresolved_result(
    *,
    reason: str,
    inputs: Mapping[str, Any] | None = None,
    receipt: Mapping[str, Any] | None = None,
    references: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "status": "unresolved_quiet_rp_owner_link",
        "reason": reason,
        "same_namespace_saved_reference_identity": False,
        "admissible_for_numeric_use": False,
        "inputs": dict(inputs or {}),
        "receipt_evidence": dict(receipt or {}),
        "reference_evidence": dict(references or {}),
        "flags": _base_flags(),
        "scope_audit": {
            "entire_archive_json_decoded": True,
            "exact_graph_nodes_explicitly_selected": False,
            "exact_graph_nodes_accessed": False,
            "graph_accessed": False,
            "graph_accessed_definition": (
                "complete JSON parsing occurred; no exact_graph_nodes node was selected, indexed, or traversed"
            ),
            "protected_references_used_only_for_equality": True,
            "source_values_selected": False,
        },
    }


def recover(
    root: str | Path,
    *,
    frame_archive_path: str | Path = DEFAULT_FRAME_ARCHIVE,
    frame_archive_sha256: str = DEFAULT_FRAME_ARCHIVE_SHA256,
    receipt_path: str | Path = DEFAULT_RECEIPT,
    receipt_sha256: str = DEFAULT_RECEIPT_SHA256,
    quiet_archive_name: str = DEFAULT_QUIET_ARCHIVE_NAME,
) -> dict[str, Any]:
    """Capture only the three saved quiet-Rp references after raw auth."""

    _require(isinstance(quiet_archive_name, str) and quiet_archive_name, "quiet archive name is missing")
    _require(
        Path(quiet_archive_name).name == quiet_archive_name
        and ".." not in Path(quiet_archive_name).parts,
        "quiet archive name must be a direct archive key",
    )
    repo_root = Path(root)
    archive, archive_meta = _read_gzip_json_after_raw_auth(
        repo_root, frame_archive_path, frame_archive_sha256
    )
    receipt, receipt_meta = _read_json_after_raw_auth(
        repo_root, receipt_path, receipt_sha256
    )
    inputs = {
        "frame_archive": archive_meta,
        "acceptance_receipt": receipt_meta,
        "frame_archive_decoded_only_after_raw_hash": True,
        "entire_archive_json_decoded": True,
        "quiet_archive_name": quiet_archive_name,
    }
    family = _family_evidence(archive, receipt, frame_archive_path, frame_archive_sha256)
    if not family["supported"]:
        return _unresolved_result(
            reason="archive_receipt_source_family_or_acceptance_unresolved",
            inputs=inputs,
            receipt=family,
        )

    try:
        bindings = _direct(
            archive,
            ("Rp_native_C5_pulse_source_bindings",),
            "archive.Rp_native_C5_pulse_source_bindings",
        )
        geometry = _direct(bindings, ("geometry",), "snapshot.geometry")
        archive_inputs = _direct(bindings, ("inputs",), "snapshot.inputs")
        quiet_input = _direct(
            archive_inputs,
            (quiet_archive_name,),
            f"snapshot.inputs[{quiet_archive_name}]",
        )
        geometry_rp = _opaque_ref(
            _direct(geometry, ("Rp",), "snapshot.geometry.Rp"),
            "snapshot.geometry.Rp",
        )
        quiet_rp = _opaque_ref(
            _direct(quiet_input, ("Rp",), f"snapshot.inputs[{quiet_archive_name}].Rp"),
            f"snapshot.inputs[{quiet_archive_name}].Rp",
        )
        endpoint_radius = _opaque_ref(
            _direct(
                quiet_input,
                ("Rp_endpoint", "radius"),
                f"snapshot.inputs[{quiet_archive_name}].Rp_endpoint.radius",
            ),
            f"snapshot.inputs[{quiet_archive_name}].Rp_endpoint.radius",
        )
        references = _reference_evidence(geometry_rp, quiet_rp, endpoint_radius)
        same_reference = references["comparisons"]["all_three_are_same_opaque_reference"]
        return {
            "schema": SCHEMA,
            "status": (
                "bounded_quiet_rp_geometry_endpoint_reference_identity_proved"
                if same_reference
                else "unresolved_quiet_rp_geometry_endpoint_reference_mismatch"
            ),
            "same_namespace_saved_reference_identity": same_reference,
            "admissible_for_numeric_use": False,
            "inputs": inputs,
            "receipt_evidence": family,
            "reference_evidence": references,
            "flags": _base_flags(),
            "scope_audit": {
                "selected_direct_fields": [
                    "Rp_native_C5_pulse_source_bindings.geometry.Rp",
                    f"Rp_native_C5_pulse_source_bindings.inputs[{quiet_archive_name}].Rp",
                    f"Rp_native_C5_pulse_source_bindings.inputs[{quiet_archive_name}].Rp_endpoint.radius",
                ],
                "entire_archive_json_decoded": True,
                "exact_graph_nodes_explicitly_selected": False,
                "exact_graph_nodes_accessed": False,
                "graph_accessed": False,
                "graph_accessed_definition": (
                    "complete JSON parsing occurred; no exact_graph_nodes node was selected, indexed, or traversed"
                ),
                "protected_references_used_only_for_equality": True,
                "protected_references_used_as_graph_or_mapping_indices": False,
                "source_values_selected": False,
                "numeric_values_computed": False,
            },
        }
    except _UnresolvedCapture as exc:
        return _unresolved_result(
            reason=exc.reason,
            inputs=inputs,
            receipt=family,
        )


build_f_quiet_rp_owner_link = recover


__all__ = [
    "CaptureInputError",
    "DEFAULT_FRAME_ARCHIVE",
    "DEFAULT_FRAME_ARCHIVE_SHA256",
    "DEFAULT_RECEIPT",
    "DEFAULT_RECEIPT_SHA256",
    "DEFAULT_QUIET_ARCHIVE_NAME",
    "PROTECTED_SOURCE_REFS",
    "SCHEMA",
    "build_f_quiet_rp_owner_link",
    "recover",
]
