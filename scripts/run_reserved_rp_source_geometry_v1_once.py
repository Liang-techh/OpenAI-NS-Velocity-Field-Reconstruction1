"""Run the authenticated reserved/quiet/native source-geometry audit once.

The caller freezes ``tmp/reserved_rp_source_geometry_v1/source.json`` before
invoking this file.  Every listed file is hashed first, and an exclusive entry
receipt is written before any gzip member is decoded or the new geometry
consumer is imported.  Historical producer/checker modules are never
imported.  The checker source is parsed as static evidence; the saved checker
receipt and the four saved quiet archive rows provide the quiet-length
positive witness.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import importlib.util
import json
import re
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
NAME = "reserved_rp_source_geometry_v1"
DIRECTORY = ROOT / "tmp" / NAME
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entry.json"
RESULT = DIRECTORY / "result.json"

GEOMETRY = "src/openai_ns_reconstruction/reserved_rp_source_geometry_v1.py"
RUNNER = "scripts/run_reserved_rp_source_geometry_v1_once.py"
RESERVED_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator.json.gz"
)
QUIET_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp.json.gz"
)
NATIVE_GRAPH = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.json.gz"
)
RESERVED_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator.py"
)
QUIET_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp.py"
)
NATIVE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.py"
)
C3_CHECK_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp_check.py"
)
C3_CHECK_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_C3_power_to_Rp_check.json"
)

RESERVED_GRAPH_SHA256 = (
    "98fd45a48d72656ec91f19a5564577baed0dbb43fe85ce8c842336856dc17a17"
)
QUIET_GRAPH_SHA256 = (
    "cb2a69ca971965a6e487bd0fa818204216747094f37819a2cc3e42a3af79911a"
)
NATIVE_GRAPH_SHA256 = (
    "859b59d422104d184092bcf60f8afa95267335ce87e8b41f2c271c18130e9a28"
)
C3_CHECK_SOURCE_SHA256 = (
    "d73cdd6bca22788b3ed4a3c97df594d5e99255c6ae2421d513151988a776b500"
)
C3_CHECK_RECEIPT_SHA256 = (
    "49eb13920a5a57acfecb237b9e606ce80afb58726be8423273940fbfc1bbf362"
)
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

FIXED_GRAPH_PINS = {
    RESERVED_GRAPH: RESERVED_GRAPH_SHA256,
    QUIET_GRAPH: QUIET_GRAPH_SHA256,
    NATIVE_GRAPH: NATIVE_GRAPH_SHA256,
}
REQUIRED_INPUTS = {
    GEOMETRY,
    RUNNER,
    RESERVED_GRAPH,
    QUIET_GRAPH,
    NATIVE_GRAPH,
    RESERVED_SOURCE,
    QUIET_SOURCE,
    NATIVE_SOURCE,
    C3_CHECK_SOURCE,
    C3_CHECK_RECEIPT,
}


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry input contract is not satisfied."""


def _require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_constant(token: str) -> Any:
    raise ValueError("non-finite JSON value is not admitted: " + token)


def _load_json(raw: bytes, label: str) -> Any:
    try:
        return json.loads(raw.decode("utf-8"), parse_constant=_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RunnerInputError(label + " is not finite UTF-8 JSON") from exc


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def _repository_path(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(), label + " is missing")
    return key, path


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    rows = source.get("files")
    _require(isinstance(rows, list), "source.files is missing")
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        label = f"source.files[{index}]"
        _require(isinstance(row, Mapping), label + " is not an object")
        key, path = _repository_path(row.get("path"), label)
        _require(key not in hashes, "duplicate source file: " + key)
        expected = row.get("raw_sha256", row.get("sha256"))
        if row.get("raw_sha256") is not None and row.get("sha256") is not None:
            _require(
                str(row["raw_sha256"]).lower() == str(row["sha256"]).lower(),
                label + " SHA fields disagree",
            )
        _require(
            isinstance(expected, str) and SHA256_RE.fullmatch(expected) is not None,
            label + " has invalid SHA-256",
        )
        actual = _sha_file(path)
        _require(actual.lower() == expected.lower(), "raw SHA-256 mismatch: " + key)
        if key.endswith(".gz"):
            _require(
                row.get("raw_sha256_before_decode") is True
                or row.get("raw_bytes_authenticated_before_decode") is True,
                "gzip input is not marked raw-authenticated: " + key,
            )
        hashes[key] = actual
    _require(REQUIRED_INPUTS.issubset(hashes), "source manifest is missing required frozen inputs")
    for path, expected in FIXED_GRAPH_PINS.items():
        _require(hashes[path].lower() == expected, "fixed graph raw SHA pin mismatch: " + path)
    _require(
        hashes[C3_CHECK_SOURCE].lower() == C3_CHECK_SOURCE_SHA256,
        "C3 checker source raw SHA pin mismatch",
    )
    _require(
        hashes[C3_CHECK_RECEIPT].lower() == C3_CHECK_RECEIPT_SHA256,
        "C3 checker receipt raw SHA pin mismatch",
    )
    return hashes


def _preflight() -> tuple[dict[str, Any], dict[str, str], str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    _require(SOURCE.is_file(), "source manifest is missing")
    source_raw = SOURCE.read_bytes()
    source = _load_json(source_raw, "source manifest")
    _require(isinstance(source, Mapping), "source manifest must be an object")
    _require(source.get("component") in (None, NAME), "source manifest component mismatch")
    _require(source.get("one_entry_only") in (None, True), "source manifest one-entry guard")
    _require(source.get("do_not_reenter") in (None, True), "source manifest do-not-reenter guard")
    hashes = _manifest_hashes(source)
    return dict(source), hashes, _sha(source_raw)


def _write_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(value)
    try:
        with path.open("xb") as stream:
            stream.write(raw)
    except FileExistsError as exc:
        raise RunnerInputError("one-entry file already exists: " + str(path)) from exc
    return _sha(raw)


def _recheck_raw(hashes: Mapping[str, str], key: str) -> bytes:
    _require(key in hashes, "unhashed input: " + key)
    raw = (ROOT / key).read_bytes()
    _require(_sha(raw).lower() == hashes[key].lower(), "input changed after entry: " + key)
    return raw


def _import_geometry() -> Any:
    spec = importlib.util.spec_from_file_location(NAME + "_consumer", ROOT / GEOMETRY)
    _require(spec is not None and spec.loader is not None, "geometry module load failed")
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[NAME + "_consumer"] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _line_hits(text: str, needle: str) -> list[int]:
    return [line_no for line_no, line in enumerate(text.splitlines(), 1) if needle in line]


def _static_quiet_length_evidence(
    checker_source_raw: bytes,
    checker_receipt_raw: bytes,
    quiet_archive: Mapping[str, Any],
) -> dict[str, Any]:
    """Read the pinned checker and receipt without importing or executing them."""

    try:
        checker_text = checker_source_raw.decode("utf-8")
        ast.parse(checker_text, filename=C3_CHECK_SOURCE, mode="exec")
    except (UnicodeDecodeError, SyntaxError) as exc:
        raise RunnerInputError("C3 checker static AST parse failed") from exc

    normalized = "".join(checker_text.split())
    required_source_fragments = {
        "four_cell_loop": "forrowinraw['actual_four_Z_C3_power_ranges']",
        "quiet_length_interval_read": "length=current.current.packets.interval(c,row['actual_positive_quiet_length'])",
        "strict_lower_positive_guard": "current.ep(length)[0]>0",
        "four_cell_count_guard": "len(raw['actual_four_Z_C3_power_ranges'])==4",
    }
    for label, fragment in required_source_fragments.items():
        _require(fragment in normalized, "pinned C3 checker guard missing: " + label)

    receipt = _load_json(checker_receipt_raw, C3_CHECK_RECEIPT)
    _require(isinstance(receipt, Mapping), "C3 checker receipt is not an object")
    _require(receipt.get("all_passed") is True, "C3 checker receipt all_passed is not True")
    geometry_receipt = receipt.get("independent_native_frame_and_geometry")
    _require(isinstance(geometry_receipt, Mapping), "C3 checker geometry receipt is missing")
    _require(
        geometry_receipt.get("exact_Rc_2Rc_Rp_and_quiet_length_identity") is True,
        "C3 checker quiet-length identity receipt is not True",
    )
    directed_receipt = receipt.get("actual_directed_ranges")
    _require(isinstance(directed_receipt, Mapping), "C3 checker directed-range receipt is missing")
    _require(
        directed_receipt.get("actual_four_original_axial_cells") == 4,
        "C3 checker receipt does not record four original axial cells",
    )

    rows = quiet_archive.get("actual_four_Z_C3_power_ranges")
    _require(isinstance(rows, list) and len(rows) == 4, "quiet archive does not contain four Z cells")
    cells: list[dict[str, Any]] = []
    seen_cells: set[tuple[str, ...]] = set()
    for index, row in enumerate(rows):
        _require(isinstance(row, Mapping), f"quiet archive Z row {index} is not an object")
        cell = row.get("exact_Z_cell")
        _require(
            isinstance(cell, list) and len(cell) == 2 and all(isinstance(item, str) for item in cell),
            f"quiet archive Z row {index} exact_Z_cell is not a two-coordinate cell",
        )
        cell_key = tuple(cell)
        _require(cell_key not in seen_cells, "quiet archive Z cells are not distinct")
        seen_cells.add(cell_key)
        length = row.get("actual_positive_quiet_length")
        _require(isinstance(length, Mapping), f"quiet archive Z row {index} quiet length is missing")
        lower_tuple = length.get("lower_exact_mpf_tuple")
        _require(
            isinstance(lower_tuple, list) and len(lower_tuple) == 4,
            f"quiet archive Z row {index} exact mpf lower tuple is missing",
        )
        sign, mantissa, exponent, precision = lower_tuple
        _require(type(sign) is int and sign == 0, f"quiet archive Z row {index} lower mpf sign is not positive")
        _require(
            type(mantissa) is int and mantissa > 0,
            f"quiet archive Z row {index} lower mpf mantissa is not strictly positive",
        )
        _require(type(exponent) is int and type(precision) is int, f"quiet archive Z row {index} mpf tuple is malformed")
        _require(
            row.get("positive_swirl_denominators_are_actual_source_functions") is True,
            f"quiet archive Z row {index} source-function positivity flag is not True",
        )
        _require(
            row.get("whole_quiet_domain_and_terminal_Rp_covered") is True,
            f"quiet archive Z row {index} quiet-domain flag is not True",
        )
        cells.append(
            {
                "exact_Z_cell": list(cell),
                "actual_positive_quiet_length_lower_exact_mpf_tuple": list(lower_tuple),
                "lower_sign": sign,
                "lower_mantissa_strictly_positive": True,
            }
        )

    return {
        "checked": True,
        "expression": "Tw-2-log(2)>0",
        "source": (
            "pinned C3_power_to_Rp_check.py static per-cell lower guard plus "
            "pinned C3_power_to_Rp_check.json all_passed receipt and four "
            "authenticated quiet archive rows"
        ),
        "checker_source_sha256": C3_CHECK_SOURCE_SHA256,
        "checker_receipt_sha256": C3_CHECK_RECEIPT_SHA256,
        "checker_source_static_guards": {
            key: _line_hits(checker_text, fragment)
            for key, fragment in (
                ("four_cell_loop", "actual_four_Z_C3_power_ranges"),
                ("quiet_length_interval_read", "actual_positive_quiet_length"),
                ("strict_lower_positive_guard", "current.ep(length)[0]>0"),
            )
        },
        "checker_receipt": {
            "all_passed": True,
            "exact_Rc_2Rc_Rp_and_quiet_length_identity": True,
            "actual_four_original_axial_cells": 4,
        },
        "four_Z_cells": cells,
    }


def _source_bindings(
    reserved_archive: Mapping[str, Any],
    quiet_archive: Mapping[str, Any],
    native_archive: Mapping[str, Any],
    evidence: Mapping[str, Any],
) -> dict[str, Any]:
    reserved_binding = reserved_archive.get("actual_background_binding")
    quiet_functions = quiet_archive.get("actual_C3_repaired_exit_power_and_pulse_functions")
    native_functions = native_archive.get("actual_Rp_native_frame_functions")
    _require(isinstance(reserved_binding, Mapping), "reserved background binding is missing")
    _require(isinstance(quiet_functions, Mapping), "quiet C3 function binding is missing")
    _require(isinstance(native_functions, Mapping), "native frame function binding is missing")
    parameters = quiet_functions.get("parameters")
    _require(isinstance(parameters, Mapping), "quiet C3 parameter binding is missing")
    source_family = quiet_archive.get("source_family")
    _require(isinstance(source_family, Mapping), "quiet source family is missing")
    return {
        "parameters": dict(parameters),
        "source_family": dict(source_family),
        "reserved": dict(reserved_binding),
        "quiet_c3": dict(quiet_functions),
        "native_frame": dict(native_functions),
        "quiet_positive_length": dict(evidence),
    }


def _alter_graph_node(
    archive: Mapping[str, Any], node: int, **changes: Any
) -> dict[str, Any]:
    rows = archive.get("exact_graph_nodes")
    _require(isinstance(rows, list), "cannot alter archive without exact_graph_nodes")
    _require(0 <= node < len(rows), "negative-case node is outside exact_graph_nodes")
    altered = dict(archive)
    altered_rows = list(rows)
    _require(isinstance(altered_rows[node], Mapping), "negative-case graph row is not an object")
    altered_rows[node] = dict(altered_rows[node], **changes)
    altered["exact_graph_nodes"] = altered_rows
    return altered


def _rejection_case(
    name: str,
    audit: Any,
    *,
    reserved: Any,
    quiet: Any,
    native: Any,
    source_bindings: Mapping[str, Any],
    support_certificate: Mapping[str, Any],
    quiet_length_evidence: Any,
) -> dict[str, Any]:
    try:
        audit(
            reserved=reserved,
            quiet_c3=quiet,
            native_frame=native,
            source_bindings=source_bindings,
            support_certificate=support_certificate,
            quiet_length_evidence=quiet_length_evidence,
        )
    except BaseException as exc:
        _require(
            type(exc).__name__ == "SourceGeometryError",
            name + " raised an unexpected exception type",
        )
        return {
            "name": name,
            "status": "passed",
            "rejected": True,
            "exception_type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    raise RunnerInputError(name + " was accepted after a post-prefix dependency change")


def _consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    keys = (
        GEOMETRY,
        RESERVED_GRAPH,
        QUIET_GRAPH,
        NATIVE_GRAPH,
        RESERVED_SOURCE,
        QUIET_SOURCE,
        NATIVE_SOURCE,
        C3_CHECK_SOURCE,
        C3_CHECK_RECEIPT,
    )
    raw = {key: _recheck_raw(hashes, key) for key in keys}
    progress["raw_recheck_complete"] = True

    # Only the new pure consumer is imported.  Historical source files remain
    # bytes for the static provenance record and are never imported.
    geometry = _import_geometry()
    progress["geometry_consumer_imported"] = True

    reserved_archive = _load_json(gzip.decompress(raw[RESERVED_GRAPH]), RESERVED_GRAPH)
    quiet_archive = _load_json(gzip.decompress(raw[QUIET_GRAPH]), QUIET_GRAPH)
    native_archive = _load_json(gzip.decompress(raw[NATIVE_GRAPH]), NATIVE_GRAPH)
    for label, archive in (
        ("reserved", reserved_archive),
        ("quiet_c3", quiet_archive),
        ("native_frame", native_archive),
    ):
        _require(isinstance(archive, Mapping), label + " graph archive is not an object")
    progress["gzip_decoded_after_entry"] = True

    evidence = _static_quiet_length_evidence(
        raw[C3_CHECK_SOURCE],
        raw[C3_CHECK_RECEIPT],
        quiet_archive,
    )
    progress["quiet_length_evidence_checked"] = True
    source_bindings = _source_bindings(
        reserved_archive,
        quiet_archive,
        native_archive,
        evidence,
    )
    support_certificate = reserved_archive.get("support_certificate")
    _require(isinstance(support_certificate, Mapping), "reserved support certificate is missing")

    reserved_graph = geometry.AuthenticatedGraph(
        role="reserved", archive=reserved_archive, raw_sha256=hashes[RESERVED_GRAPH]
    )
    quiet_graph = geometry.AuthenticatedGraph(
        role="quiet_c3", archive=quiet_archive, raw_sha256=hashes[QUIET_GRAPH]
    )
    native_graph = geometry.AuthenticatedGraph(
        role="native_frame", archive=native_archive, raw_sha256=hashes[NATIVE_GRAPH]
    )
    report = geometry.audit_reserved_rp_source_geometry(
        reserved=reserved_graph,
        quiet_c3=quiet_graph,
        native_frame=native_graph,
        source_bindings=source_bindings,
        support_certificate=support_certificate,
        quiet_length_evidence=geometry.QuietLengthEvidence(
            checked=True,
            expression="Tw-2-log(2)>0",
            source=evidence["source"],
        ),
    )
    progress["primary_geometry_audit_calls"] = 1
    _require(report.source_recipe_full_support is True, "source-recipe support did not close")

    negative_cases = [
        _rejection_case(
            "post_prefix_Rc_dependency_change_rejected",
            geometry.audit_reserved_rp_source_geometry,
            reserved=reserved_graph,
            quiet=geometry.AuthenticatedGraph(
                role="quiet_c3",
                archive=_alter_graph_node(quiet_archive, 25494, argument=4892),
                raw_sha256=hashes[QUIET_GRAPH],
            ),
            native=native_graph,
            source_bindings=source_bindings,
            support_certificate=support_certificate,
            quiet_length_evidence=evidence,
        ),
        _rejection_case(
            "post_prefix_Tw_dependency_change_rejected",
            geometry.audit_reserved_rp_source_geometry,
            reserved=reserved_graph,
            quiet=geometry.AuthenticatedGraph(
                role="quiet_c3",
                archive=_alter_graph_node(quiet_archive, 25007, arguments=[25005, 25005]),
                raw_sha256=hashes[QUIET_GRAPH],
            ),
            native=native_graph,
            source_bindings=source_bindings,
            support_certificate=support_certificate,
            quiet_length_evidence=evidence,
        ),
    ]
    progress["negative_geometry_audit_calls"] = len(negative_cases)
    _require(all(case.get("status") == "passed" for case in negative_cases), "negative graph cases did not reject")

    return {
        "geometry_report": report.report(),
        "quiet_length_evidence": evidence,
        "negative_cases": negative_cases,
        "audit_calls": {
            "primary_successful_audit": 1,
            "negative_rejection_audits": len(negative_cases),
        },
        "raw_sha256_pins": {
            "reserved_graph": hashes[RESERVED_GRAPH],
            "quiet_c3_graph": hashes[QUIET_GRAPH],
            "native_frame_graph": hashes[NATIVE_GRAPH],
            "C3_checker_source": hashes[C3_CHECK_SOURCE],
            "C3_checker_receipt": hashes[C3_CHECK_RECEIPT],
        },
        "claims": {
            "source_recipe_full_support": True,
            "live_owner_proved": False,
            "full_c1_matrix_proved": False,
            "physical_m15_proved": False,
            "historical_producer_imported": False,
            "historical_fixture_called": False,
            "old_math_recomputed": False,
        },
    }


def _error_record(exc: BaseException) -> dict[str, str]:
    return {
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
    }


def main() -> int:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    if ENTRY.exists() or RESULT.exists():
        print(json.dumps({"component": NAME, "status": "rejected", "error": "component already entered"}))
        return 1

    source_sha: str | None = None
    hashes: dict[str, str] = {}
    preflight_error: dict[str, str] | None = None
    try:
        _, hashes, source_sha = _preflight()
    except BaseException as exc:
        preflight_error = _error_record(exc)

    entry_value: dict[str, Any] = {
        "component": NAME,
        "source_sha256": source_sha,
        "input_hashes": hashes,
        "entered": True,
        "preflight_passed": preflight_error is None,
        "raw_sha256_verified_before_import_or_decode": preflight_error is None,
        "gzip_decode_before_entry": False,
        "ancestor_producer_called": False,
        "historical_producer_imported": False,
        "historical_fixture_called": False,
        "manifest_created_by_runner": False,
    }
    if preflight_error is not None:
        entry_value["error"] = preflight_error
    try:
        entry_sha = _write_exclusive(ENTRY, entry_value)
    except BaseException as exc:
        print(json.dumps({"component": NAME, "status": "failed", "error": _error_record(exc)}))
        return 1

    progress: dict[str, Any] = {
        "raw_recheck_complete": False,
        "geometry_consumer_imported": False,
        "gzip_decoded_after_entry": False,
        "quiet_length_evidence_checked": False,
        "primary_geometry_audit_calls": 0,
        "negative_geometry_audit_calls": 0,
    }
    outcome: dict[str, Any] = {
        "component": NAME,
        "status": "failed",
        "entry_sha256": entry_sha,
        "source_sha256": source_sha,
        "input_hashes": hashes,
        "ancestor_producer_called": False,
        "historical_producer_imported": False,
        "historical_fixture_called": False,
        "raw_sha256_verified_before_import_or_decode": preflight_error is None,
        "progress": progress,
        "output": None,
    }
    if preflight_error is not None:
        outcome["error"] = preflight_error
    else:
        try:
            output = _consume(hashes, progress)
            outcome["output"] = output
            outcome["status"] = "passed"
        except BaseException as exc:
            outcome["error"] = _error_record(exc)
        finally:
            outcome["progress"] = progress
    try:
        result_sha = _write_exclusive(RESULT, outcome)
    except BaseException as exc:
        print(json.dumps({"component": NAME, "status": "failed", "entry_sha256": entry_sha, "error": _error_record(exc)}))
        return 1
    print(
        json.dumps(
            {
                "component": NAME,
                "status": outcome["status"],
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "error": (outcome.get("error") or {}).get("message"),
            },
            sort_keys=True,
        )
    )
    return 0 if outcome["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
