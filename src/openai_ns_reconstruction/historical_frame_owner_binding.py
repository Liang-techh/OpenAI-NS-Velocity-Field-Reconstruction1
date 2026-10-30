"""Authenticate and bind the saved historical ``logRp`` frame owner.

This leaf is intentionally a narrow receipt consumer.  It reads only the
saved receipt metadata and the six named fields in the three authenticated
archives.  Every raw input is hashed before any gzip payload is decoded.  No
graph, source function, old constructor, or pressure report is imported or
traversed here.

The returned packet records the historical frame node and Gamma's imported
node as saved references.  Equality below is equality of those opaque saved
references only; it is not a comparison of graph-node values or a claim
about ordinary-``Z`` independence, current source-function equivalence, or
endpoint closure.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


SCHEMA = "historical_frame_owner_binding_v1"

PREFIX = "experiments/root_st073/lei_ren_part1_paper_compliant_"
PRE_RP_ARCHIVE = PREFIX + "current_original_n1_preRp_pressure_source.json.gz"
GAMMA_ARCHIVE = PREFIX + "current_original_n1_Gamma_future_functions.json.gz"
FRAME_ARCHIVE = PREFIX + "current_original_C3_Rp_native_frame.json.gz"
PRE_RP_CHECK = PREFIX + "current_original_n1_preRp_pressure_source_check.json"
GAMMA_CHECK = PREFIX + "current_original_n1_Gamma_future_functions_check.json"
FRAME_CHECK = PREFIX + "current_original_C3_Rp_native_frame_check.json"

ARCHIVES = {
    "preRp": PRE_RP_ARCHIVE,
    "Gamma": GAMMA_ARCHIVE,
    "frame": FRAME_ARCHIVE,
}
CHECKS = {
    "preRp": PRE_RP_CHECK,
    "Gamma": GAMMA_CHECK,
    "frame": FRAME_CHECK,
}

# These values are a fallback for direct library use.  The one-shot runner
# always supplies the caller-owned frozen source.json manifest, and the
# manifest values take precedence there.
KNOWN_RAW_SHA256 = {
    PRE_RP_ARCHIVE: "9b2eedd03d4544acb8fff9c175875c8f35848e2992a82106e0a01c2b8e331119",
    GAMMA_ARCHIVE: "dbc3dd702930b85599e7967b5e14d934f3c326a731cd1b9528990158ff74c612",
    FRAME_ARCHIVE: "859b59d422104d184092bcf60f8afa95267335ce87e8b41f2c271c18130e9a28",
    PRE_RP_CHECK: "ea874353136a936d2e637354ec604bfd7e0f77af8b1fcbabfac2f185dc04e541",
    GAMMA_CHECK: "ee6fc38317768a68b5775f6224eeb250f6383ee3f732d83c68b30316c720f927",
    FRAME_CHECK: "d34d0e4f90b131e68bbde16caad33029ff89b75174741ad02d03d0a6e135d9e8",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha256(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _read_json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw)
    except Exception as exc:  # pragma: no cover - malformed frozen input
        raise ValueError(f"{label} is not valid JSON") from exc
    _require(isinstance(value, dict), f"{label} is not a JSON object")
    return value


def _manifest_hashes(source: Mapping[str, Any], root: Path) -> dict[str, str]:
    """Extract the file hash map from the caller-owned source checkpoint.

    Checkpoints in this repository use either ``sha256``/``input_hashes``
    maps or a ``files`` list.  Supporting both keeps the leaf independent of
    checkpoint formatting while still requiring the six exact inputs below.
    """

    result: dict[str, str] = {}

    def add(name: Any, digest: Any, where: str) -> None:
        _require(isinstance(name, str), f"{where} has a non-text path")
        path = Path(name)
        if path.is_absolute():
            try:
                name = path.resolve().relative_to(root.resolve()).as_posix()
            except ValueError as exc:
                raise ValueError(f"{where} path is outside the repository") from exc
        else:
            _require(".." not in path.parts, f"{where} path escapes the repository")
            name = path.as_posix()
        _require(isinstance(digest, str), f"{where} has a non-text SHA-256")
        normalized = digest.lower()
        _require(
            len(normalized) == 64
            and all(char in "0123456789abcdef" for char in normalized),
            f"{where} has an invalid SHA-256",
        )
        prior = result.get(name)
        _require(prior is None or prior == normalized, f"conflicting SHA-256 for {name}")
        result[name] = normalized

    for key in ("sha256", "input_hashes", "input_sha256", "code_sha256"):
        value = source.get(key)
        if isinstance(value, Mapping):
            for name, digest in value.items():
                add(name, digest, f"source.{key}[{name!r}]")

    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            add(item.get("path"), item.get("sha256"), f"source.files[{index}]")

    return result


def _expected_hashes(
    source: Mapping[str, Any] | None, root: Path
) -> dict[str, str]:
    expected = dict(KNOWN_RAW_SHA256)
    if source is None:
        return expected
    manifest = _manifest_hashes(source, root)
    for relative in KNOWN_RAW_SHA256:
        _require(relative in manifest, f"source checkpoint does not pin {relative}")
        expected[relative] = manifest[relative]
    return expected


def _read_and_authenticate_raw(
    root: Path, relative: str, expected_sha256: str
) -> tuple[bytes, str]:
    path = root / relative
    _require(path.is_file(), f"missing frozen input: {relative}")
    raw = path.read_bytes()
    actual = _sha256(raw)
    _require(actual == expected_sha256, f"raw SHA-256 mismatch: {relative}")
    return raw, actual


def _saved_node(value: Any, label: str) -> int:
    _require(isinstance(value, int) and not isinstance(value, bool), f"{label} is not an integer node")
    return value


def _saved_field(container: Mapping[str, Any], key: str, label: str) -> Any:
    _require(key in container, f"missing saved field: {label}.{key}")
    return container[key]


def recover(
    root: str | Path,
    source_checkpoint: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a strict saved-reference receipt for the historical log-radius frame.

    ``source_checkpoint`` is the already-loaded frozen ``source.json`` passed
    by the runner.  It is consumed only for the six raw input hashes; no
    archive is decoded until all six hashes have passed.
    """

    root_path = Path(root)
    expected = _expected_hashes(source_checkpoint, root_path)
    required = (
        ("preRp", PRE_RP_ARCHIVE),
        ("Gamma", GAMMA_ARCHIVE),
        ("frame", FRAME_ARCHIVE),
        ("preRp", PRE_RP_CHECK),
        ("Gamma", GAMMA_CHECK),
        ("frame", FRAME_CHECK),
    )

    # Authenticate every raw byte first.  In particular, no gzip.decompress
    # call is reachable until this complete loop has succeeded.
    raw_inputs: dict[str, bytes] = {}
    raw_sha256: dict[str, str] = {}
    for _owner, relative in required:
        raw, actual = _read_and_authenticate_raw(root_path, relative, expected[relative])
        raw_inputs[relative] = raw
        raw_sha256[relative] = actual

    receipt_raw = {
        owner: _read_json(raw_inputs[path], path)
        for owner, path in CHECKS.items()
    }
    receipt_input_hashes: dict[str, dict[str, Any]] = {}
    receipt_source_family: dict[str, Any] = {}
    receipt_all_passed: dict[str, bool] = {}
    archive_for_owner = {
        "preRp": PRE_RP_ARCHIVE,
        "Gamma": GAMMA_ARCHIVE,
        "frame": FRAME_ARCHIVE,
    }

    for owner, receipt in receipt_raw.items():
        _require(receipt.get("all_passed") is True, f"{owner} receipt did not pass")
        family = receipt.get("source_family")
        _require(family is not None, f"{owner} receipt has no source_family")
        hashes = receipt.get("input_hashes")
        _require(isinstance(hashes, Mapping), f"{owner} receipt has no input_hashes map")
        own_archive = archive_for_owner[owner]
        _require(
            hashes.get(Path(own_archive).name) == expected[own_archive],
            f"{owner} receipt does not pin its own archive SHA",
        )
        receipt_all_passed[owner] = True
        receipt_source_family[owner] = family
        receipt_input_hashes[owner] = {
            key: hashes.get(Path(key).name)
            for key in (own_archive, FRAME_ARCHIVE)
        }

    _require(
        receipt_source_family["preRp"]
        == receipt_source_family["Gamma"]
        == receipt_source_family["frame"],
        "receipt source_family metadata differs",
    )

    pre_frame_pin = receipt_input_hashes["preRp"].get(FRAME_ARCHIVE)
    gamma_frame_pin = receipt_input_hashes["Gamma"].get(FRAME_ARCHIVE)
    _require(pre_frame_pin is not None, "preRp receipt has no frame archive pin")
    _require(gamma_frame_pin is not None, "Gamma receipt has no frame archive pin")
    _require(pre_frame_pin == gamma_frame_pin, "preRp and Gamma frame pins differ")
    _require(pre_frame_pin == expected[FRAME_ARCHIVE], "receipt frame pin differs from raw frame SHA")

    # Only now decode the three gzip archives.  Their payloads are accessed by
    # the explicitly named paths below; no graph or nested node walk occurs.
    archives = {
        owner: _read_json(gzip.decompress(raw_inputs[path]), path)
        for owner, path in ARCHIVES.items()
    }
    _require(all(archives[owner].get("source_family") == receipt_source_family[owner]
                 for owner in ARCHIVES), "archive and receipt source families differ")
    pre_geometry = _saved_field(
        _saved_field(archives["preRp"], "source_data", "preRp"),
        "geometry",
        "preRp.source_data",
    )
    _require(isinstance(pre_geometry, Mapping), "preRp source_data.geometry is not an object")
    pre_original = _saved_node(
        _saved_field(pre_geometry, "original_logRp_node", "preRp.source_data.geometry"),
        "preRp.source_data.geometry.original_logRp_node",
    )
    accepted_frame_sha = _saved_field(
        pre_geometry,
        "accepted_native_frame_sha256",
        "preRp.source_data.geometry",
    )
    _require(isinstance(accepted_frame_sha, str), "accepted_native_frame_sha256 is not text")
    _require(
        accepted_frame_sha.lower() == expected[FRAME_ARCHIVE],
        "preRp accepted_native_frame_sha256 does not identify the frame archive",
    )

    frame_functions = _saved_field(
        archives["frame"],
        "actual_Rp_native_frame_functions",
        "frame",
    )
    _require(isinstance(frame_functions, Mapping), "frame actual_Rp_native_frame_functions is not an object")
    frame_original = _saved_node(
        _saved_field(frame_functions, "original_logRp", "frame.actual_Rp_native_frame_functions"),
        "frame.actual_Rp_native_frame_functions.original_logRp",
    )

    gamma_bindings = _saved_field(archives["Gamma"], "function_bindings", "Gamma")
    _require(isinstance(gamma_bindings, Mapping), "Gamma function_bindings is not an object")
    native_binding = _saved_field(gamma_bindings, "native_binding", "Gamma.function_bindings")
    schedule = _saved_field(gamma_bindings, "schedule", "Gamma.function_bindings")
    _require(isinstance(native_binding, Mapping), "Gamma native_binding is not an object")
    _require(isinstance(schedule, Mapping), "Gamma schedule is not an object")
    _require(native_binding.get("frame_archive") == Path(FRAME_ARCHIVE).name,
             "Gamma native binding names a different frame archive")
    gamma_original = _saved_node(
        _saved_field(native_binding, "original_logRp", "Gamma.function_bindings.native_binding"),
        "Gamma.function_bindings.native_binding.original_logRp",
    )
    gamma_imported = _saved_node(
        _saved_field(native_binding, "imported_logRp", "Gamma.function_bindings.native_binding"),
        "Gamma.function_bindings.native_binding.imported_logRp",
    )
    gamma_schedule = _saved_node(
        _saved_field(schedule, "logRp", "Gamma.function_bindings.schedule"),
        "Gamma.function_bindings.schedule.logRp",
    )

    checks = {
        "all_three_receipts_all_passed": all(receipt_all_passed.values()),
        "receipt_source_family_equal": (
            receipt_source_family["preRp"]
            == receipt_source_family["Gamma"]
            == receipt_source_family["frame"]
        ),
        "archive_and_receipt_source_families_equal": True,
        "Gamma_native_binding_names_frame_archive": True,
        "receipt_self_archive_sha256": {
            owner: receipt_input_hashes[owner][archive_for_owner[owner]]
            == expected[archive_for_owner[owner]]
            for owner in archive_for_owner
        },
        "preRp_and_Gamma_frame_input_sha256_equal": pre_frame_pin == gamma_frame_pin,
        "preRp_frame_input_sha256_matches_raw_frame": pre_frame_pin == expected[FRAME_ARCHIVE],
        "accepted_native_frame_sha256_matches_raw_frame": (
            accepted_frame_sha.lower() == expected[FRAME_ARCHIVE]
        ),
        "preRp_original_logRp_matches_frame_original_logRp": pre_original == frame_original,
        "Gamma_original_logRp_matches_frame_original_logRp": gamma_original == frame_original,
        "historical_frame_original_logRp_shared": (
            pre_original == frame_original == gamma_original
        ),
        "Gamma_imported_logRp_matches_schedule_logRp": gamma_imported == gamma_schedule,
    }
    checks["all_required_checks_passed"] = bool(
        checks["all_three_receipts_all_passed"]
        and checks["receipt_source_family_equal"]
        and all(checks["receipt_self_archive_sha256"].values())
        and checks["preRp_and_Gamma_frame_input_sha256_equal"]
        and checks["preRp_frame_input_sha256_matches_raw_frame"]
        and checks["accepted_native_frame_sha256_matches_raw_frame"]
        and checks["historical_frame_original_logRp_shared"]
        and checks["Gamma_imported_logRp_matches_schedule_logRp"]
    )
    _require(checks["all_required_checks_passed"], "historical frame owner checks failed")

    return {
        "schema": SCHEMA,
        "status": "saved_historical_same_frame_receipt",
        "source_only": True,
        "archive_paths": dict(ARCHIVES),
        "receipt_paths": dict(CHECKS),
        "archive_raw_sha256": {
            owner: raw_sha256[path] for owner, path in ARCHIVES.items()
        },
        "receipt_raw_sha256": {
            owner: raw_sha256[path] for owner, path in CHECKS.items()
        },
        "receipts": {
            owner: {
                "all_passed": receipt_all_passed[owner],
                "source_family": receipt_source_family[owner],
                "input_hashes": receipt_input_hashes[owner],
                "self_archive": archive_for_owner[owner],
                "self_archive_sha256": receipt_input_hashes[owner][archive_for_owner[owner]],
                "frame_input_sha256": receipt_input_hashes[owner].get(FRAME_ARCHIVE),
            }
            for owner in ("preRp", "Gamma", "frame")
        },
        "saved_numbers": {
            "preRp_source_data_geometry_original_logRp_node": pre_original,
            "frame_actual_Rp_native_frame_functions_original_logRp": frame_original,
            "Gamma_function_bindings_native_binding_original_logRp": gamma_original,
            "Gamma_function_bindings_native_binding_imported_logRp": gamma_imported,
            "Gamma_function_bindings_schedule_logRp": gamma_schedule,
        },
        "saved_frame_sha256": {
            "preRp_source_data_geometry_accepted_native_frame_sha256": accepted_frame_sha,
            "raw_frame_archive_sha256": expected[FRAME_ARCHIVE],
            "preRp_receipt_frame_input_sha256": pre_frame_pin,
            "Gamma_receipt_frame_input_sha256": gamma_frame_pin,
        },
        "checks": checks,
        "binding": {
            "historical_frame_owner": "frame.actual_Rp_native_frame_functions.original_logRp",
            "preRp_geometry_original_logRp": pre_original,
            "frame_original_logRp": frame_original,
            "Gamma_native_original_logRp": gamma_original,
            "Gamma_imported_logRp": gamma_imported,
            "Gamma_schedule_logRp": gamma_schedule,
            "same_historical_frame_original_node": True,
            "same_Gamma_imported_schedule_node": True,
        },
        "nonclaims": {
            "ordinary_Z_independence": False,
            "current_source_function_equivalence": False,
            "endpoint_closure": False,
            "cross_namespace_node_value_comparison": False,
            "graph_traversal": False,
            "source_values_computed": False,
            "old_methods_called": False,
            "pressure_report_traversed": False,
        },
        "runtime_model_metadata": None,
    }


__all__ = ["recover", "SCHEMA"]
