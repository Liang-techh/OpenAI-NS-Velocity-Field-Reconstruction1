"""Fail-closed admission for the two canonical Gamma source archives.

This module is intentionally a provenance gate, not a source evaluator.  It
has one archive-admission entry point, :func:`admit_canonical_gamma_sources`,
whose paths, raw SHA-256 pins, receipt pins, and JSON pointers are all fixed in
this file.  The historical ``.py`` files are hash-only inputs: they are never
imported or executed.  The two gzip archives are decoded only after every one
of the six canonical files has passed its raw-byte hash check and both small
receipt documents have passed their flags and archive cross-pins.

The returned records preserve direct ``function_bindings`` metadata only.
An integer or list in a saved role is opaque and is not interpreted as a
function value or a proven node identity. This leaf makes no P/F ownership, common-source, numerical,
certified-bound, or mathematical-value claim.

``synthetic_negative_fixture`` is pure in-memory validation.  It does not read
the repository, decode gzip, import a historical source, or call the admission
entry point; it exercises the rejection gates against deliberately malformed
metadata.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Any, Mapping


SCHEMA = "gamma_source_admission_v1"

# These are the only accepted source paths.  They are the paths recorded by
# the frozen Gamma weighted-energy checkpoint and the preceding runner/result.
GAMMA_STEM = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions"
)
COMPLETE_FUTURE_STEM = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions"
)

GAMMA_PY = GAMMA_STEM + ".py"
GAMMA_CHECK_JSON = GAMMA_STEM + "_check.json"
GAMMA_ARCHIVE = GAMMA_STEM + ".json.gz"
COMPLETE_FUTURE_PY = COMPLETE_FUTURE_STEM + ".py"
COMPLETE_FUTURE_CHECK_JSON = COMPLETE_FUTURE_STEM + "_check.json"
COMPLETE_FUTURE_ARCHIVE = COMPLETE_FUTURE_STEM + ".json.gz"

GAMMA_ARCHIVE_SHA256 = "dbc3dd702930b85599e7967b5e14d934f3c326a731cd1b9528990158ff74c612"
GAMMA_PY_SHA256 = "38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01"
GAMMA_CHECK_SHA256 = "ee6fc38317768a68b5775f6224eeb250f6383ee3f732d83c68b30316c720f927"
COMPLETE_FUTURE_ARCHIVE_SHA256 = (
    "b63c1cdae91fa6b083a0bc7582c28bd5f9664225f223639b28eea8060ab7185b"
)
COMPLETE_FUTURE_PY_SHA256 = (
    "ace67570df36883b7a7faefae49f95e92fd8f4ce3a80126cdec46b6b778c67ac"
)
COMPLETE_FUTURE_CHECK_SHA256 = (
    "fd7fbf5ffd6305198c9f7ed4b14f2fa2d6e660d2f414a0debd78fc140812860f"
)

# The order is deliberate: _authenticate_raw_files hashes all six byte
# strings first and compares all six pins before any gzip decoder is called.
CANONICAL_FILE_PINS: tuple[tuple[str, str], ...] = (
    (GAMMA_PY, GAMMA_PY_SHA256),
    (GAMMA_CHECK_JSON, GAMMA_CHECK_SHA256),
    (GAMMA_ARCHIVE, GAMMA_ARCHIVE_SHA256),
    (COMPLETE_FUTURE_PY, COMPLETE_FUTURE_PY_SHA256),
    (COMPLETE_FUTURE_CHECK_JSON, COMPLETE_FUTURE_CHECK_SHA256),
    (COMPLETE_FUTURE_ARCHIVE, COMPLETE_FUTURE_ARCHIVE_SHA256),
)

# These direct pointers come from the previous runner/result and the new
# handoff.  No caller can replace them with a nearby node or a pseudo-path.
# The returned ``reference`` is metadata only; it is never called or treated
# as a numerical/function value.
GAMMA_POINTERS: Mapping[str, str] = MappingProxyType({
    "source_S": "/function_bindings/schedule/S",
    "delta": "/function_bindings/schedule/delta",
    "epsilon": "/function_bindings/schedule/epsilon",
    "waiting": "/function_bindings/schedule/waiting",
    "logRp": "/function_bindings/schedule/logRp",
    "logRtail": "/function_bindings/schedule/logRtail",
    # The saved log-Rp bridge independently locates this direct child port.
    "native_binding": "/function_bindings/native_binding/imported_logRp",
    "Dhat": "/function_bindings/functions/Dhat",
    "collar_energy": "/function_bindings/pieces/collar/energy",
    "outer_energy": "/function_bindings/pieces/exterior/energy",
    "actual_scaled_deficit": "/function_bindings/pieces/actual_scaled_deficit",
    "preheat_bracket": "/function_bindings/pieces/preheat_bracket",
    "a": "/function_bindings/pieces/a",
})

COMPLETE_FUTURE_POINTERS: Mapping[str, str] = MappingProxyType({
    "flatten": "/function_bindings/pieces/flatten",
    "power": "/function_bindings/pieces/power",
    "angular_change": "/function_bindings/pieces/angular_change",
    "steep_waiting_preheat_tail": "/function_bindings/pieces/steep_waiting_preheat_tail",
    "Gamma_heat_difference": "/function_bindings/pieces/Gamma_heat_difference",
    "Nrel": "/function_bindings/pieces/Nrel",
    "tail_multiplier": "/function_bindings/pieces/tail_multiplier",
    "energy_complete": "/function_bindings/energy/complete",
    "energy_weighted": "/function_bindings/energy/weighted",
    "selected_ap": "/function_bindings/selected_pulse/ap",
    "selected_C1": "/function_bindings/selected_pulse/C1",
    "selected_C2": "/function_bindings/selected_pulse/C2",
    "selector_target": "/function_bindings/pulse_equations/target",
    "selector_A2": "/function_bindings/pulse_equations/A2",
    "selector_A1": "/function_bindings/pulse_equations/A1",
    "selector_A0": "/function_bindings/pulse_equations/A0",
    "selector_discriminant": "/function_bindings/pulse_equations/discriminant",
})

# The handoff asks for actual A/pre/radius ownership, but neither the previous
# runner nor the static handoff supplied a fixed JSON pointer for those roles.
# Keep that limitation explicit instead of inventing a nearby path.
UNRESOLVED_STATIC_ROLE_PATHS: Mapping[str, tuple[str, ...]] = MappingProxyType({
    "gamma": ("P_F_A_owners", "P_F_pre_owners", "F_physical_radius"),
    "complete_future": (),
})

RECEIPT_FLAGS: Mapping[str, tuple[str, ...]] = MappingProxyType({
    "gamma": (
        "current_original_n1_exact_Gamma_three_future_C1_functions_installed",
    ),
    "complete_future": (
        "current_original_n1_complete_future_energy_C1_function_installed",
        "current_original_n1_actual_small_angular_C1_functions_installed",
        "current_original_n1_actual_selected_pulse_coefficient_C1_functions_installed",
        "original_pulse_positive_branch_proved",
    ),
})


class GammaSourceAdmissionError(RuntimeError):
    """Raised whenever canonical source admission cannot be established."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise GammaSourceAdmissionError(message)


def _sha256_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _canonical_relative(raw: Any, label: str) -> str:
    """Validate a slash-separated repository-relative path without aliases."""

    _require(isinstance(raw, str) and bool(raw), f"{label} path is not text")
    _require("\\" not in raw, f"{label} path uses a non-canonical separator")
    path = PurePosixPath(raw)
    _require(not path.is_absolute(), f"{label} path is absolute")
    _require(".." not in path.parts and "." not in path.parts, f"{label} path is pseudo")
    _require(path.as_posix() == raw, f"{label} path is not canonical")
    return raw


def _canonical_file(root: Path, relative: str, label: str) -> Path:
    relative = _canonical_relative(relative, label)
    path = root.joinpath(*relative.split("/"))
    try:
        path.resolve().relative_to(root)
    except ValueError as exc:
        raise GammaSourceAdmissionError(f"{label} path escapes workspace") from exc
    return path


def _workspace_root(value: Any) -> Path:
    try:
        root = Path(value).resolve()
    except (TypeError, ValueError, OSError) as exc:
        raise GammaSourceAdmissionError("workspace root is not a usable path") from exc
    _require(root.is_dir(), "workspace root is not a directory")
    return root


def _check_hashes(expected: Mapping[str, str], actual: Mapping[str, str]) -> None:
    """Compare a complete hash set only after all actual hashes are available."""

    for relative, digest in expected.items():
        _require(actual.get(relative) == digest, f"raw SHA-256 mismatch: {relative}")


def _authenticate_raw_files(root: Path) -> tuple[dict[str, bytes], dict[str, str]]:
    """Read and hash all six raw files, then compare every pin.

    The returned gzip bytes are the exact bytes that were hashed.  Keeping them
    in memory prevents a later file replacement from changing what gets
    decoded and makes the raw-before-decompression ordering explicit.
    """

    raw_bytes: dict[str, bytes] = {}
    actual: dict[str, str] = {}
    expected = {relative: digest for relative, digest in CANONICAL_FILE_PINS}
    for relative, _ in CANONICAL_FILE_PINS:
        path = _canonical_file(root, relative, "canonical input")
        _require(path.is_file(), f"canonical input is missing: {relative}")
        try:
            raw = path.read_bytes()
        except OSError as exc:
            raise GammaSourceAdmissionError(f"canonical input cannot be read: {relative}") from exc
        raw_bytes[relative] = raw
        actual[relative] = _sha256_bytes(raw)

    # No decoder or JSON parser has been called for an archive at this point.
    _check_hashes(expected, actual)
    return raw_bytes, actual


def _json_document(raw: bytes, label: str) -> Mapping[str, Any]:
    try:
        document = json.loads(raw)
    except (TypeError, ValueError, UnicodeDecodeError) as exc:
        raise GammaSourceAdmissionError(f"{label} is not valid JSON") from exc
    _require(isinstance(document, Mapping), f"{label} is not a JSON object")
    return document


def _required_bool(document: Mapping[str, Any], key: str, label: str) -> bool:
    _require(key in document, f"{label} is missing required field: {key}")
    value = document[key]
    _require(type(value) is bool, f"{label}.{key} is not boolean")
    _require(value is True, f"{label}.{key} is not true")
    return True


def _archive_basename(relative: str) -> str:
    return PurePosixPath(relative).name


def _validate_receipt(
    group: str,
    document: Mapping[str, Any],
    archive_relative: str,
    archive_sha256: str,
) -> dict[str, Any]:
    """Check only the small receipt fields needed for this source gate."""

    _required_bool(document, "all_passed", f"{group} receipt")
    for flag in RECEIPT_FLAGS[group]:
        _required_bool(document, flag, f"{group} receipt")
    input_hashes = document.get("input_hashes")
    _require(isinstance(input_hashes, Mapping), f"{group} receipt input_hashes is missing")
    archive_name = _archive_basename(archive_relative)
    _require(archive_name in input_hashes, f"{group} receipt lacks archive cross-pin")
    cross_pin = input_hashes[archive_name]
    _require(cross_pin == archive_sha256, f"{group} receipt archive cross-pin mismatch")
    return {
        "all_passed": True,
        "required_flags": {flag: True for flag in RECEIPT_FLAGS[group]},
        "archive_sha256_pointer": f"/input_hashes/{archive_name}",
        "archive_sha256": cross_pin,
        "input_hashes_list_materialized": False,
    }


def _validate_complete_gamma_cross_pin(
    complete_document: Mapping[str, Any], gamma_archive_sha256: str
) -> dict[str, Any]:
    input_hashes = complete_document.get("input_hashes")
    _require(isinstance(input_hashes, Mapping), "complete future receipt input_hashes is missing")
    gamma_name = _archive_basename(GAMMA_ARCHIVE)
    _require(gamma_name in input_hashes, "complete future receipt lacks Gamma archive cross-pin")
    cross_pin = input_hashes[gamma_name]
    _require(cross_pin == gamma_archive_sha256, "complete future receipt Gamma archive cross-pin mismatch")
    return {
        "json_pointer": f"/input_hashes/{gamma_name}",
        "reference": cross_pin,
        "reference_is_value": False,
        "matches_gamma_archive": True,
    }


def _json_pointer(document: Any, pointer: str) -> Any:
    """Resolve one fixed RFC 6901 pointer, without reference/graph traversal."""

    _require(isinstance(pointer, str) and pointer.startswith("/"), "fixed JSON pointer is malformed")
    current = document
    for raw_token in pointer.split("/")[1:]:
        token = raw_token.replace("~1", "/").replace("~0", "~")
        if isinstance(current, Mapping):
            _require(token in current, f"fixed JSON pointer component is absent: {token}")
            current = current[token]
        elif isinstance(current, list):
            _require(token != "" and token.isdigit(), "fixed JSON pointer list component is invalid")
            index = int(token)
            _require(index < len(current), f"fixed JSON pointer list index is absent: {token}")
            current = current[index]
        else:
            raise GammaSourceAdmissionError("fixed JSON pointer descends through a scalar")
    return current


def _read_direct_roles(
    group: str, document: Mapping[str, Any], pointers: Mapping[str, str], archive_relative: str
) -> dict[str, Any]:
    roles: dict[str, Any] = {}
    missing: list[str] = []
    for role, pointer in pointers.items():
        _require(pointer.startswith("/function_bindings/"), f"{group}.{role} leaves function_bindings")
        try:
            reference = _json_pointer(document, pointer)
        except GammaSourceAdmissionError as exc:
            missing.append(role)
            roles[role] = {
                "path": archive_relative,
                "json_pointer": pointer,
                "available": False,
                "missing_reason": str(exc),
                "metadata_interpretation": "opaque_uninterpreted",
            }
        else:
            roles[role] = {
                "path": archive_relative,
                "json_pointer": pointer,
                "reference": reference,
                "available": True,
                "metadata_interpretation": "opaque_uninterpreted",
            }
    _require(not missing, f"{group} archive has missing fixed roles: {', '.join(missing)}")
    return roles


def _decode_and_read_roles(
    raw_archive: bytes, group: str, archive_relative: str, pointers: Mapping[str, str],
    receipt: Mapping[str, Any],
) -> dict[str, Any]:
    """Decode one already-authenticated gzip and read only fixed pointers."""

    try:
        decoded = gzip.decompress(raw_archive)
    except (OSError, EOFError, ValueError) as exc:
        raise GammaSourceAdmissionError(f"{group} archive gzip decoding failed") from exc
    document = _json_document(decoded, f"{group} archive")
    _require("source_family" in document and "source_family" in receipt,
             f"{group} source family is missing")
    _require(document["source_family"] == receipt["source_family"],
             f"{group} archive/receipt source family differs")
    roles = _read_direct_roles(group, document, pointers, archive_relative)
    return {
        "path": archive_relative,
        "encoding": "gzip_json",
        "roles": roles,
        "source_family": document["source_family"],
        "own_receipt_family_matches": True,
        "raw_authenticated_before_decompression": True,
        "source_graph_traversed": False,
        "direct_json_pointer_only": True,
        "role_metadata_is_opaque": True,
    }


def admit_canonical_gamma_sources(workspace_root: str | Path) -> dict[str, Any]:
    """Admit only the pinned Gamma and complete-future archive role metadata.

    The function raises :class:`GammaSourceAdmissionError` on every failed
    gate.  It never accepts caller-provided paths, hashes, or pointers.  A
    successful result still contains metadata-only references and explicitly
    leaves P/F ownership and all source-value claims false.
    """

    root = _workspace_root(workspace_root)
    raw_bytes, actual_hashes = _authenticate_raw_files(root)

    # Receipt JSON is read only after all six raw hashes have passed.  No
    # archive bytes are decoded until both receipt documents and their pins pass.
    gamma_receipt_doc = _json_document(raw_bytes[GAMMA_CHECK_JSON], "Gamma receipt")
    complete_receipt_doc = _json_document(
        raw_bytes[COMPLETE_FUTURE_CHECK_JSON], "complete future receipt"
    )
    gamma_receipt = _validate_receipt(
        "gamma", gamma_receipt_doc, GAMMA_ARCHIVE, GAMMA_ARCHIVE_SHA256
    )
    complete_receipt = _validate_receipt(
        "complete_future",
        complete_receipt_doc,
        COMPLETE_FUTURE_ARCHIVE,
        COMPLETE_FUTURE_ARCHIVE_SHA256,
    )
    complete_receipt["gamma_archive_cross_pin"] = _validate_complete_gamma_cross_pin(
        complete_receipt_doc, GAMMA_ARCHIVE_SHA256
    )

    # Only now may gzip decoding occur, and the decoder receives the exact raw
    # bytes whose hashes were checked above.
    gamma_archive = _decode_and_read_roles(
        raw_bytes[GAMMA_ARCHIVE], "gamma", GAMMA_ARCHIVE, GAMMA_POINTERS, gamma_receipt_doc
    )
    complete_archive = _decode_and_read_roles(
        raw_bytes[COMPLETE_FUTURE_ARCHIVE],
        "complete_future",
        COMPLETE_FUTURE_ARCHIVE,
        COMPLETE_FUTURE_POINTERS,
        complete_receipt_doc,
    )

    return {
        "schema": SCHEMA,
        "admitted": True,
        "canonical_paths_are_fixed": True,
        "canonical_files": {
            relative: {
                "expected_sha256": expected,
                "actual_sha256": actual_hashes[relative],
                "raw_sha256_verified": True,
                "read_mode": "raw_bytes_hash_only" if not relative.endswith(".json.gz") else "raw_bytes_hash_then_decode",
            }
            for relative, expected in CANONICAL_FILE_PINS
        },
        "receipts": {
            "gamma": {"path": GAMMA_CHECK_JSON, **gamma_receipt},
            "complete_future": {"path": COMPLETE_FUTURE_CHECK_JSON, **complete_receipt},
        },
        "archives": {"gamma": gamma_archive, "complete_future": complete_archive},
        "unresolved_static_role_paths": {
            group: list(names) for group, names in UNRESOLVED_STATIC_ROLE_PATHS.items()
        },
        "source_graph_traversed": False,
        "direct_json_pointer_only": True,
        "role_metadata_is_opaque": True,
        "source_values_computed": False,
        "source_value_computation": False,
        "certified": False,
        "P_F_owner_correspondence_proved": False,
        "common_a_pre_A_proved": False,
        "full_future_P_F_proved": False,
    }


def _expect_rejection(label: str, action: Any, checks: dict[str, bool]) -> None:
    try:
        action()
    except GammaSourceAdmissionError:
        checks[label] = True
    else:
        checks[label] = False


def synthetic_negative_fixture() -> dict[str, Any]:
    """Run pure in-memory negative controls for every admission gate.

    This fixture deliberately never calls the public admission interface and
    never reads a filesystem path.  It is therefore not evidence about the
    canonical archives themselves.
    """

    checks: dict[str, bool] = {}
    fake_sha = "a" * 64
    other_sha = "b" * 64

    _expect_rejection(
        "pseudo_path_rejected",
        lambda: _canonical_relative("../experiments/root_st073/Gamma.json.gz", "synthetic"),
        checks,
    )
    _expect_rejection(
        "raw_hash_mismatch_rejected",
        lambda: _check_hashes({"synthetic.json.gz": fake_sha}, {"synthetic.json.gz": other_sha}),
        checks,
    )
    _expect_rejection(
        "receipt_missing_rejected",
        lambda: _validate_receipt("gamma", {}, GAMMA_ARCHIVE, GAMMA_ARCHIVE_SHA256),
        checks,
    )
    _expect_rejection(
        "receipt_false_all_passed_rejected",
        lambda: _validate_receipt(
            "gamma",
            {
                "all_passed": False,
                RECEIPT_FLAGS["gamma"][0]: True,
                "input_hashes": {_archive_basename(GAMMA_ARCHIVE): GAMMA_ARCHIVE_SHA256},
            },
            GAMMA_ARCHIVE,
            GAMMA_ARCHIVE_SHA256,
        ),
        checks,
    )
    _expect_rejection(
        "receipt_false_required_flag_rejected",
        lambda: _validate_receipt(
            "gamma",
            {
                "all_passed": True,
                RECEIPT_FLAGS["gamma"][0]: False,
                "input_hashes": {_archive_basename(GAMMA_ARCHIVE): GAMMA_ARCHIVE_SHA256},
            },
            GAMMA_ARCHIVE,
            GAMMA_ARCHIVE_SHA256,
        ),
        checks,
    )
    complete_with_wrong_gamma_pin = {
        "all_passed": True,
        **{flag: True for flag in RECEIPT_FLAGS["complete_future"]},
        "input_hashes": {
            _archive_basename(COMPLETE_FUTURE_ARCHIVE): COMPLETE_FUTURE_ARCHIVE_SHA256,
            _archive_basename(GAMMA_ARCHIVE): other_sha,
        },
    }
    _expect_rejection(
        "complete_receipt_gamma_cross_pin_rejected",
        lambda: (
            _validate_receipt(
                "complete_future",
                complete_with_wrong_gamma_pin,
                COMPLETE_FUTURE_ARCHIVE,
                COMPLETE_FUTURE_ARCHIVE_SHA256,
            ),
            _validate_complete_gamma_cross_pin(complete_with_wrong_gamma_pin, GAMMA_ARCHIVE_SHA256),
        ),
        checks,
    )
    _expect_rejection(
        "missing_role_rejected",
        lambda: _read_direct_roles(
            "gamma",
            {"function_bindings": {"schedule": {"S": 123}}},
            {"source_S": GAMMA_POINTERS["source_S"], "delta": GAMMA_POINTERS["delta"]},
            GAMMA_ARCHIVE,
        ),
        checks,
    )

    return {
        "schema": SCHEMA,
        "fixture": "synthetic_negative_archive_admission",
        "checks": checks,
        "checks_passed": bool(checks) and all(checks.values()),
        "filesystem_read": False,
        "gzip_decoded": False,
        "historical_source_executed": False,
        "source_values_computed": False,
        "certified": False,
        "P_F_owner_correspondence_proved": False,
        "full_future_P_F_proved": False,
    }


__all__ = [
    "SCHEMA",
    "GAMMA_POINTERS",
    "COMPLETE_FUTURE_POINTERS",
    "CANONICAL_FILE_PINS",
    "UNRESOLVED_STATIC_ROLE_PATHS",
    "RECEIPT_FLAGS",
    "GammaSourceAdmissionError",
    "admit_canonical_gamma_sources",
    "synthetic_negative_fixture",
]
