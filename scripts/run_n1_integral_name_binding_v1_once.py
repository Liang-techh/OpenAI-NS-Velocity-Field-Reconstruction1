"""Run the n1 name-binding candidate boundary exactly once.

This runner is intentionally a transport boundary, not an evaluator.  It
freezes the new leaf, its synthetic fixture, this runner, the two previously
saved receipt archives, and the five actual source files.  It authenticates
all raw bytes and the pending source checkpoint, writes the exclusive entry
record, and only then imports the new leaf and fixture.  The synthetic
fixture runs once.  Only after the fixture passes are the two pinned gzip
receipts decoded and the source files decoded as text.

The only post-fixture leaf call is::

    leaf.analyze(saved_result, source_texts)

The leaf receives one bounded adapter object made from the decoded scoped
catalog result plus the saved occurrence records, and a five-alias source-text
mapping.  This runner does not decode or traverse a graph, execute an old
producer, evaluate a numeric value or integral, install a binder, or claim an
owner, identity, m15, or parent.  A successful output is only a conditional
source-name-binding candidate.
"""

from __future__ import annotations

import builtins
import gzip
import hashlib
import importlib.util
import json
import math
import os
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "n1_integral_name_binding_v1"
PHASE = COMPONENT + "_pending"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/n1_integral_name_binding_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_integral_name_binding_fixture_v1.py"
RUNNER_REL = "scripts/run_n1_integral_name_binding_v1_once.py"

PRIOR_SCOPED_REL = (
    "experiments/root_st073/n1_scoped_graph_evaluator_v1_receipt.json.gz"
)
PRIOR_OCCURRENCE_REL = (
    "experiments/root_st073/n1_integral_occurrence_saved_audit_v1_receipt.json.gz"
)

SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_Rc_function_transport.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_Rc_functional_controls.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.py",
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_C3_functions.py",
)

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    RUNNER_REL,
    PRIOR_SCOPED_REL,
    PRIOR_OCCURRENCE_REL,
    *SOURCE_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)
PRIOR_REL = (PRIOR_SCOPED_REL, PRIOR_OCCURRENCE_REL)

# These are the accepted raw-byte pins for the already-saved inputs.  The
# leaf, fixture, and runner hashes are supplied by the source manifest after
# those new files have been frozen by the parent.
PRIOR_SHA256 = {
    PRIOR_SCOPED_REL: "2cf73004ec11b38777182cdb2ce6c96d82e699331a868c9056ff88a19c10678a",
    PRIOR_OCCURRENCE_REL: "b2bcf11ead7e6e5991225e394436b2bae34a22191c854dbf3ca70d5f3efffdfb",
}
SOURCE_SHA256 = {
    SOURCE_REL[0]: "57863c9c3e5746f770ece127cb233193d0e022735dd69994b8ca06c9a584d9d5",
    SOURCE_REL[1]: "cb6f11d4a62d59a45f05d0961236d5c6236bd544dca8b18499f68b2d8d98d009",
    SOURCE_REL[2]: "d71ff19a68cfed8ace9a6857ffeca860ecee90d8d89b35b2292a7310e180fb59",
    SOURCE_REL[3]: "d00dfe285735844eabc985f99d94190caa1cfd187f668eb14a2686055cf175fa",
    SOURCE_REL[4]: "c2b80f895f03c9a4157b056e572f564843620049aaced34061efd8eb357a699e",
}

SOURCE_TEXT_ALIASES = {
    SOURCE_REL[0]: "graph",
    SOURCE_REL[1]: "controls",
    SOURCE_REL[2]: "core",
    SOURCE_REL[3]: "lower_collar",
    SOURCE_REL[4]: "c3_derivative",
}

LEAF_SCHEMA = COMPONENT
LEAF_RESULT_SCHEMA = COMPONENT + "_result"
FIXTURE_SCHEMA = "n1_integral_name_binding_fixture_v1"
SCOPED_COMPONENT = "n1_scoped_graph_evaluator_v1"
SCOPED_RESULT_SCHEMA = SCOPED_COMPONENT + "_runner_result_v1"
OCCURRENCE_COMPONENT = "n1_integral_occurrence_saved_audit_v1"
OCCURRENCE_RESULT_SCHEMA = OCCURRENCE_COMPONENT + "_runner_result_v1"
OCCURRENCE_LEAF_COMPONENT = "n1_integral_scoped_occurrence_v1"
OCCURRENCE_LEAF_SCHEMA = OCCURRENCE_COMPONENT + "_result"
SHA256_HEX = frozenset("0123456789abcdefABCDEF")

_LOADED_NAMES = (
    "_n1_integral_name_binding_v1_runner_leaf",
    "_n1_integral_name_binding_v1_runner_fixture",
)
_MISSING = object()


class RunnerInputError(ValueError):
    """Raised for a fail-closed input, provenance, or scope violation."""


def _require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_safe(value: Any, seen: set[int] | None = None) -> Any:
    """Serialize an ordinary result without invoking project operations."""

    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"type": "opaque_bytes", "length": len(value), "sha256": _sha_bytes(value)}

    active = seen if seen is not None else set()
    identity = id(value)
    if identity in active:
        return {"type": "opaque_cycle", "repr": repr(value)}
    active.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): _json_safe(item, active) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(item, active) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((_json_safe(item, active) for item in value), key=repr)
        return {"type": "opaque_value", "class": type(value).__name__, "repr": repr(value)}
    finally:
        active.discard(identity)


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            _json_safe(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _write_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(value)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
    descriptor = os.open(str(path), flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
    finally:
        if descriptor != -1:
            os.close(descriptor)
    return _sha_bytes(raw)


def _decode_object(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _sha_text(value: Any, label: str) -> str:
    _require(
        isinstance(value, str)
        and len(value) == 64
        and all(character in SHA256_HEX for character in value),
        "invalid SHA-256 for " + label,
    )
    return value.lower()


def _repo_file(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " path must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " path is missing or escapes repository: " + key,
    )
    return key, path


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    files = source.get("files")
    _require(isinstance(files, list), "source.files list is missing")
    hashes: dict[str, str] = {}
    for index, row in enumerate(files):
        _require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_file(row.get("path"), f"source.files[{index}]")
        _require(key in FROZEN_SET, "manifest contains out-of-scope file: " + key)
        _require(key not in hashes, "duplicate manifest path: " + key)
        expected = _sha_text(row.get("sha256"), key)
        actual = _sha_file(path)
        _require(actual == expected, "raw SHA mismatch: " + key)
        if key in PRIOR_REL:
            _require(
                row.get("raw_sha256_before_decode") is True,
                "prior gzip is not marked raw-pinned before decode: " + key,
            )
        if key in SOURCE_REL:
            _require(
                row.get("source_text_supplied") is True,
                "actual source is not marked for text supply: " + key,
            )
        hashes[key] = actual

    _require(
        frozenset(hashes) == FROZEN_SET,
        "manifest must contain exactly the ten frozen inputs",
    )
    for key, expected in PRIOR_SHA256.items():
        _require(hashes[key] == expected, "accepted prior gzip raw SHA mismatch: " + key)
    for key, expected in SOURCE_SHA256.items():
        _require(hashes[key] == expected, "accepted actual source SHA mismatch: " + key)
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _validate_contract(source: Mapping[str, Any]) -> None:
    contract = source.get("input_contract")
    _require(isinstance(contract, Mapping), "input_contract must be an object")
    expected = {
        "fixture_api": "run_fixture(leaf)",
        "leaf_api": "analyze(saved_result, source_texts)",
        "prior_receipt_paths": list(PRIOR_REL),
        "source_paths": list(SOURCE_REL),
        "source_text_aliases": dict(SOURCE_TEXT_ALIASES),
        "scoped_receipt_component": SCOPED_COMPONENT,
        "scoped_receipt_schema": SCOPED_RESULT_SCHEMA,
        "scoped_receipt_status": "returned",
        "occurrence_receipt_component": OCCURRENCE_COMPONENT,
        "occurrence_receipt_schema": OCCURRENCE_RESULT_SCHEMA,
        "occurrence_receipt_status": "returned",
        "occurrence_leaf_component": OCCURRENCE_LEAF_COMPONENT,
        "occurrence_leaf_schema": OCCURRENCE_LEAF_SCHEMA,
        "occurrence_leaf_status": "validated_saved_failure",
        "record_id_field": "integral_node_id",
        "record_id_count": 7,
        "record_ids_one_to_one": True,
        "decode_after_fixture": True,
        "authenticate_raw_gzip_before_decode": True,
        "no_graph_decode": True,
        "no_graph_traversal": True,
        "no_graph_evaluation": True,
        "no_numeric_evaluation": True,
        "no_old_producer_import": True,
        "no_binder_install": True,
        "conditional_source_name_binding_candidate_only": True,
        "source_ast_consumer_inside_leaf": True,
        "leaf_invocations_after_fixture": 1,
    }
    for key, value in expected.items():
        _require(key in contract, "input_contract key is missing: " + key)
        _require(contract.get(key) == value, key + " contract mismatch")
    for key in ("prior_scoped_sha256", "prior_occurrence_sha256"):
        expected_sha = PRIOR_SHA256[PRIOR_REL[0] if key == "prior_scoped_sha256" else PRIOR_REL[1]]
        if key in contract:
            _require(contract.get(key) == expected_sha, key + " contract mismatch")
    if "actual_source_sha256" in contract:
        _require(
            contract.get("actual_source_sha256") == dict(SOURCE_SHA256),
            "actual_source_sha256 contract mismatch",
        )


def _brief_record(brief: Mapping[str, Any]) -> Mapping[str, Any]:
    if brief.get("component") == COMPONENT:
        return brief
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        return nested
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping) and continuation.get("component") == COMPONENT:
        return continuation
    return brief


def _checkpoint_candidate(
    brief: Mapping[str, Any], record: Mapping[str, Any]
) -> Mapping[str, Any] | None:
    candidates: list[Any] = [record.get("source_checkpoint"), brief.get("source_checkpoint")]
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        candidates.append(nested.get("source_checkpoint"))
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping):
        candidates.append(continuation.get("source_checkpoint"))
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            return candidate
    return None


def _checkpoint_path(raw_path: Any) -> Path:
    _require(isinstance(raw_path, str) and raw_path, "pending source checkpoint path is missing")
    path_value = Path(raw_path.replace("\\", "/"))
    if path_value.is_absolute():
        return path_value.resolve()
    _require(".." not in path_value.parts, "pending source checkpoint escapes repository")
    return (ROOT / path_value).resolve()


def _validate_optional_preimage(checkpoint: Mapping[str, Any]) -> None:
    path_value = checkpoint.get("preimage_path")
    expected = checkpoint.get("preimage_sha256")
    if path_value is None and expected is None:
        return
    _require(
        isinstance(path_value, str) and isinstance(expected, str),
        "pending preimage pointer is incomplete",
    )
    path = Path(path_value).resolve()
    _require(path.is_file(), "pending preimage is missing")
    _require(
        _sha_text(expected, "pending preimage") == _sha_file(path),
        "pending preimage SHA mismatch",
    )


def _validate_pending_brief(source_sha: str) -> tuple[str, str]:
    _require(BRIEF.is_file(), "pending resume brief is missing")
    brief_raw = BRIEF.read_bytes()
    brief_sha = _sha_bytes(brief_raw)
    brief = _decode_object(brief_raw, "pending resume brief")
    record = _brief_record(brief)
    _require(record.get("component") == COMPONENT, "pending brief component mismatch")
    _require(record.get("phase") == PHASE, "pending brief phase mismatch")
    _require(record.get("one_entry_only") is True, "pending brief one-entry guard missing")
    _require(record.get("do_not_reenter") is True, "pending brief re-entry guard missing")
    if "goal_active" in record:
        _require(record.get("goal_active") is True, "pending brief goal is not active")
    if "full_goal_complete" in record:
        _require(
            record.get("full_goal_complete") is False,
            "pending brief marks the full goal complete",
        )
    checkpoint = _checkpoint_candidate(brief, record)
    _require(checkpoint is not None, "pending brief source_checkpoint is missing")
    checkpoint_path = _checkpoint_path(checkpoint.get("path"))
    _require(
        checkpoint_path == SOURCE.resolve(),
        "pending brief points to another source checkpoint",
    )
    _require(
        _sha_text(checkpoint.get("sha256"), "pending source checkpoint") == source_sha,
        "pending brief source checkpoint SHA mismatch",
    )
    _validate_optional_preimage(checkpoint)
    return brief_sha, SOURCE.relative_to(ROOT).as_posix()


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered or has a result")
    _require(SOURCE.is_file(), "source manifest is missing")
    source_raw = SOURCE.read_bytes()
    source_sha = _sha_bytes(source_raw)
    source = _decode_object(source_raw, "source manifest")
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    _require(source.get("phase") == PHASE, "source manifest phase mismatch")
    _require(source.get("one_entry_only") is True, "source manifest one-entry guard missing")
    _require(source.get("do_not_reenter") is True, "source manifest re-entry guard missing")
    _validate_contract(source)
    hashes = _manifest_hashes(source)
    brief_sha, checkpoint_key = _validate_pending_brief(source_sha)
    return source, hashes, source_sha, brief_sha, checkpoint_key


def _reauthenticate_after_entry(
    hashes: Mapping[str, str], source_sha: str, brief_sha: str
) -> None:
    """Recheck manifest, brief, and all frozen raw inputs before new import."""

    source_raw = SOURCE.read_bytes()
    _require(_sha_bytes(source_raw) == source_sha, "source manifest changed after exclusive entry")
    source = _decode_object(source_raw, "source manifest after entry")
    _require(source.get("component") == COMPONENT, "post-entry source component mismatch")
    _require(source.get("phase") == PHASE, "post-entry source phase mismatch")
    _require(source.get("one_entry_only") is True, "post-entry source one-entry guard missing")
    _require(source.get("do_not_reenter") is True, "post-entry source re-entry guard missing")
    _validate_contract(source)
    current_hashes = _manifest_hashes(source)
    _require(current_hashes == dict(hashes), "frozen input manifest changed after exclusive entry")

    brief_raw = BRIEF.read_bytes()
    _require(_sha_bytes(brief_raw) == brief_sha, "pending brief changed after exclusive entry")
    current_brief_sha, _ = _validate_pending_brief(source_sha)
    _require(current_brief_sha == brief_sha, "pending brief checkpoint changed after exclusive entry")

    for relative in FROZEN_INPUTS:
        _require(
            _sha_file(ROOT / relative) == hashes[relative],
            "frozen input changed before new import: " + relative,
        )


def _read_frozen_raw(hashes: Mapping[str, str]) -> dict[str, bytes]:
    """Read raw bytes and authenticate them without decoding the prior gzip."""

    raw: dict[str, bytes] = {}
    for relative in FROZEN_INPUTS:
        data = (ROOT / relative).read_bytes()
        _require(
            _sha_bytes(data) == hashes[relative],
            "frozen input changed after fixture: " + relative,
        )
        raw[relative] = data
    return raw


def _forbidden_project_modules() -> list[str]:
    """Reject historical package/source modules already present in the process."""

    found: list[str] = []
    package_root = (ROOT / "src" / "openai_ns_reconstruction").resolve()
    experiments_root = (ROOT / "experiments").resolve()
    for name, module in list(sys.modules.items()):
        if name in _LOADED_NAMES:
            continue
        if name == "openai_ns_reconstruction" or name.startswith("openai_ns_reconstruction."):
            found.append(name)
            continue
        module_file = getattr(module, "__file__", None)
        if not isinstance(module_file, str):
            continue
        try:
            resolved = Path(module_file).resolve()
        except OSError:
            continue
        if resolved.is_relative_to(package_root) or resolved.is_relative_to(experiments_root):
            found.append(name)
    return sorted(set(found))


def _guarded_import(original: Callable[..., Any]) -> Callable[..., Any]:
    def checked(
        name: str,
        globals_: Any = None,
        locals_: Any = None,
        fromlist: tuple[str, ...] | None = None,
        level: int = 0,
    ) -> Any:
        resolved = name
        if level:
            package = globals_.get("__package__", "") if isinstance(globals_, Mapping) else ""
            try:
                resolved = importlib.util.resolve_name(name, package)
            except (ImportError, AttributeError):
                resolved = name
        if resolved == "openai_ns_reconstruction" or resolved.startswith("openai_ns_reconstruction."):
            if resolved not in _LOADED_NAMES:
                raise RunnerInputError("historical package import blocked: " + resolved)
        if resolved == "experiments" or resolved.startswith("experiments."):
            raise RunnerInputError("actual source/producer import blocked: " + resolved)
        return original(name, globals_, locals_, fromlist, level)

    return checked


def _cleanup_new_modules() -> None:
    for name in _LOADED_NAMES:
        sys.modules.pop(name, None)


def _load_module(relative: str, module_name: str) -> ModuleType:
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(module_name, path)
    _require(spec is not None and spec.loader is not None, "new module spec is unavailable: " + relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _load_new_modules() -> tuple[ModuleType, ModuleType]:
    _require(not _forbidden_project_modules(), "historical package or source module already loaded")
    original_import = builtins.__import__
    builtins.__import__ = _guarded_import(original_import)
    try:
        leaf = _load_module(LEAF_REL, _LOADED_NAMES[0])
        fixture = _load_module(FIXTURE_REL, _LOADED_NAMES[1])
        _require(not _forbidden_project_modules(), "historical module loaded by new files")
        return leaf, fixture
    except BaseException:
        _cleanup_new_modules()
        raise
    finally:
        builtins.__import__ = original_import


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not a mapping")
    result = dict(value)
    _require(result.get("schema") == FIXTURE_SCHEMA, "fixture schema mismatch")
    _require(result.get("passed") is True, "synthetic fixture did not pass")
    _require(result.get("synthetic_only") is True, "fixture is not synthetic-only")
    checks = result.get("checks")
    _require(
        isinstance(checks, Mapping) and bool(checks),
        "fixture checks must be a non-empty mapping",
    )
    _require(all(item is True for item in checks.values()), "synthetic fixture contains a failed check")
    for key in (
        "graph_read",
        "graph_decoded",
        "graph_traversed",
        "source_graph_imported",
        "source_graph_executed",
        "old_producer_imported",
        "old_leaf_called",
        "numeric_values_computed",
        "actual_integral_value_computed",
        "integrals_evaluated",
        "binder_mapping_installed",
        "binder_identity_proved",
        "candidate_actual_integral_binder",
        "candidate_actual_integral_binder_proved",
        "actual_graph_owner_identity",
        "source_callback_certified",
        "source_function_equality_proved",
        "actual_owner_proved",
        "identity_proved",
        "parent_identity_proved",
        "m15_claimed",
        "parent_claimed",
    ):
        if key in result:
            _require(result[key] in (False, 0, None), "fixture reported forbidden activity: " + key)
    if "leaf_calls" in result:
        _require(result.get("leaf_calls") == 0, "fixture called the saved-record leaf entry")
    return _json_safe(result)


def _decode_gzip_json(raw: bytes, label: str) -> dict[str, Any]:
    """Decode one authenticated saved receipt, never a graph archive."""

    try:
        decoded = gzip.decompress(raw)
        value = json.loads(decoded)
    except (OSError, UnicodeDecodeError, gzip.BadGzipFile, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid authenticated gzip JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _receipt_envelope_candidates(value: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """Unwrap only the known receipt envelope layers.

    This is deliberately not a recursive walk.  It never follows graph
    fields, instruction operands, or arbitrary nested objects.
    """

    candidates: list[Mapping[str, Any]] = [value]
    first = value.get("result")
    if isinstance(first, Mapping):
        candidates.append(first)
        second = first.get("result")
        if isinstance(second, Mapping):
            candidates.append(second)
        output = first.get("output")
        if isinstance(output, Mapping):
            candidates.append(output)
    output = value.get("output")
    if isinstance(output, Mapping):
        candidates.append(output)
    # Published audit receipts keep the saved rows under the runner result's
    # leaf_output field: root.result.leaf_output.saved_leaf_records.  Add
    # this one named layer only; do not recursively search receipt contents.
    for candidate in tuple(candidates):
        leaf_output = candidate.get("leaf_output")
        if isinstance(leaf_output, Mapping):
            candidates.append(leaf_output)
    return candidates


def _scoped_saved_result(
    receipt: Mapping[str, Any], label: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Extract the saved seven-row scoped catalog without graph inspection."""

    for candidate in _receipt_envelope_candidates(receipt):
        catalog_value: Any = None
        catalog_key: str | None = None
        for key in ("integral_catalog", "catalog"):
            value = candidate.get(key)
            if isinstance(value, Mapping) and isinstance(value.get("records"), list):
                catalog_value = dict(value)
                catalog_key = key
                break
        if catalog_value is None:
            for key in ("integral_candidates", "integral_records"):
                value = candidate.get(key)
                if isinstance(value, list):
                    catalog_value = {"records": list(value)}
                    catalog_key = key
                    break
        if catalog_value is None:
            progress = candidate.get("progress")
            if isinstance(progress, Mapping):
                for key in ("integral_catalog", "catalog", "integral_candidates", "integral_records"):
                    value = progress.get(key)
                    if isinstance(value, Mapping) and isinstance(value.get("records"), list):
                        catalog_value = dict(value)
                        catalog_key = "progress." + key
                        break
                    if isinstance(value, list):
                        catalog_value = {"records": list(value)}
                        catalog_key = "progress." + key
                        break
        if catalog_value is None:
            continue
        _require(
            candidate.get("component") == SCOPED_COMPONENT,
            label + " component mismatch",
        )
        _require(
            candidate.get("schema") == SCOPED_RESULT_SCHEMA,
            label + " result schema mismatch",
        )
        _require(candidate.get("status") == "returned", label + " status is not returned")
        _require(candidate.get("passed") is True, label + " passed marker is not true")
        records = catalog_value.get("records")
        _require(
            catalog_value.get("schema") == SCOPED_COMPONENT + "_integral_catalog_v1",
            label + " catalog schema mismatch",
        )
        _require(
            catalog_value.get("integral_record_count") == 7,
            label + " catalog integral_record_count is not seven",
        )
        _require(isinstance(records, list) and len(records) == 7, label + " must contain seven catalog records")
        integral_ids: list[int] = []
        for index, row in enumerate(records):
            _require(isinstance(row, Mapping), f"{label} catalog record {index} is not an object")
            integral_id = row.get("integral_node_id")
            _require(
                type(integral_id) is int and integral_id >= 0,
                f"{label} catalog record {index} integral_node_id is malformed",
            )
            integral_ids.append(integral_id)
        _require(
            len(set(integral_ids)) == 7,
            label + " catalog integral_node_id values are not unique",
        )
        saved_result = dict(candidate)
        saved_result["integral_catalog"] = catalog_value
        saved_result["integral_candidates"] = list(records)
        return saved_result, {
            "catalog_key": catalog_key,
            "catalog_record_count": len(records),
            "integral_node_ids": integral_ids,
            "top_level_component": candidate.get("component"),
            "top_level_schema": candidate.get("schema"),
        }
    raise RunnerInputError(label + " has no bounded seven-row integral catalog")


def _occurrence_rows(
    receipt: Mapping[str, Any], label: str
) -> tuple[list[Mapping[str, Any]], list[Mapping[str, Any]], dict[str, Any]]:
    """Extract saved occurrence rows from the known audit envelope only."""

    for candidate in _receipt_envelope_candidates(receipt):
        leaf_candidate: Mapping[str, Any] | None = None
        leaf_output = candidate.get("leaf_output")
        if isinstance(leaf_output, Mapping):
            _require(
                candidate.get("component") == OCCURRENCE_COMPONENT,
                label + " runner component mismatch",
            )
            _require(
                candidate.get("schema") == OCCURRENCE_RESULT_SCHEMA,
                label + " runner result schema mismatch",
            )
            _require(candidate.get("status") == "returned", label + " runner status is not returned")
            _require(candidate.get("passed") is True, label + " runner passed marker is not true")
            leaf_candidate = leaf_output
        elif candidate.get("component") == OCCURRENCE_LEAF_COMPONENT:
            # This direct form is retained only for a saved leaf result; the
            # published receipt normally reaches it through root.result.
            leaf_candidate = candidate
        if leaf_candidate is None:
            continue
        _require(
            leaf_candidate.get("component") == OCCURRENCE_LEAF_COMPONENT,
            label + " leaf component mismatch",
        )
        _require(
            leaf_candidate.get("schema") == OCCURRENCE_LEAF_SCHEMA,
            label + " leaf result schema mismatch",
        )
        _require(
            leaf_candidate.get("status") == "validated_saved_failure",
            label + " leaf status mismatch",
        )
        _require(
            leaf_candidate.get("saved_status") == "failed",
            label + " leaf saved_status mismatch",
        )
        _require(
            leaf_candidate.get("fixture_passed") is True,
            label + " leaf fixture marker is not true",
        )
        _require(leaf_candidate.get("leaf_calls") == 1, label + " leaf call count is not one")
        _require(
            leaf_candidate.get("catalog_record_count") == 7,
            label + " leaf catalog record count is not seven",
        )
        rows = leaf_candidate.get("saved_leaf_records")
        catalog_rows = leaf_candidate.get("saved_catalog_records")
        _require(
            isinstance(rows, list) and len(rows) == 7,
            label + " saved_leaf_records must contain seven rows",
        )
        _require(
            isinstance(catalog_rows, list) and len(catalog_rows) == 7,
            label + " saved_catalog_records must contain seven rows",
        )
        checked: list[Mapping[str, Any]] = []
        checked_catalog: list[Mapping[str, Any]] = []
        occurrence_ids: list[int] = []
        catalog_ids: list[int] = []
        for index, row in enumerate(rows):
            _require(isinstance(row, Mapping), f"{label} saved occurrence record {index} is not an object")
            integral_id = row.get("integral_node_id")
            _require(
                type(integral_id) is int and integral_id >= 0,
                f"{label} saved occurrence record {index} integral_node_id is malformed",
            )
            checked.append(row)
            occurrence_ids.append(integral_id)
        for index, row in enumerate(catalog_rows):
            _require(isinstance(row, Mapping), f"{label} saved catalog record {index} is not an object")
            integral_id = row.get("integral_node_id")
            _require(
                type(integral_id) is int and integral_id >= 0,
                f"{label} saved catalog record {index} integral_node_id is malformed",
            )
            checked_catalog.append(row)
            catalog_ids.append(integral_id)
        _require(len(set(occurrence_ids)) == 7, label + " saved occurrence IDs are not unique")
        _require(len(set(catalog_ids)) == 7, label + " saved catalog IDs are not unique")
        _require(occurrence_ids == catalog_ids, label + " saved occurrence/catalog IDs do not correspond")
        return checked, checked_catalog, {
            "occurrence_key": "leaf_output.saved_leaf_records",
            "occurrence_record_count": len(checked),
            "catalog_record_count": len(checked_catalog),
            "integral_node_ids": occurrence_ids,
            "catalog_integral_node_ids": catalog_ids,
            "top_level_component": candidate.get("component"),
            "top_level_schema": candidate.get("schema"),
            "leaf_component": leaf_candidate.get("component"),
            "leaf_schema": leaf_candidate.get("schema"),
        }
    raise RunnerInputError(label + " has no bounded seven-row saved occurrence result")


def _build_saved_result(
    scoped_receipt: Mapping[str, Any],
    occurrence_receipt: Mapping[str, Any],
    scoped_path: str,
    occurrence_path: str,
    scoped_sha256: str,
    occurrence_sha256: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Join the two prior receipt results into the leaf's real input shape."""

    saved_result, scoped_meta = _scoped_saved_result(scoped_receipt, scoped_path)
    occurrence_records, occurrence_catalog_records, occurrence_meta = _occurrence_rows(
        occurrence_receipt, occurrence_path
    )
    scoped_ids = scoped_meta.get("integral_node_ids")
    occurrence_ids = occurrence_meta.get("integral_node_ids")
    occurrence_catalog_ids = occurrence_meta.get("catalog_integral_node_ids")
    _require(
        isinstance(scoped_ids, list)
        and isinstance(occurrence_ids, list)
        and isinstance(occurrence_catalog_ids, list),
        "prior receipt ID metadata is incomplete",
    )
    _require(
        scoped_ids == occurrence_ids == occurrence_catalog_ids,
        "scoped catalog and saved occurrence IDs do not correspond one-to-one",
    )
    # The leaf's occurrence adapter accepts a holder with ``records``.  This
    # is a bounded saved-record join, not a graph or witness traversal.
    saved_occurrence_output = {"records": list(occurrence_records)}
    saved_result["occurrence_output"] = saved_occurrence_output
    saved_result["saved_occurrence_output"] = saved_occurrence_output
    saved_result["saved_occurrence_counts_preserved"] = True
    meta = {
        "scoped_receipt": {
            "path": scoped_path,
            "raw_sha256": scoped_sha256,
            **scoped_meta,
        },
        "occurrence_receipt": {
            "path": occurrence_path,
            "raw_sha256": occurrence_sha256,
            **occurrence_meta,
        },
        "joined_catalog_record_count": 7,
        "joined_occurrence_record_count": 7,
        "joined_integral_node_ids": list(scoped_ids),
        "occurrence_catalog_rows_supplied": len(occurrence_catalog_records),
        "component_status_schema_authenticated": True,
        "graph_decoded": False,
        "graph_traversed": False,
    }
    return saved_result, meta


def _receipt_summary(value: Mapping[str, Any], raw_sha256: str, path: str) -> dict[str, Any]:
    """Record only top-level receipt metadata; do not walk saved records."""

    return {
        "path": path,
        "raw_sha256": raw_sha256,
        "decoded": True,
        "top_level_keys": sorted(str(key) for key in value),
        "top_level_component": value.get("component"),
        "top_level_schema": value.get("schema"),
    }


_FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "actual_owner",
        "actual_owner_name",
        "actual_owner_proved",
        "actual_owner_verified",
        "owner_identity",
        "owner_identity_proved",
        "owner_identity_verified",
        "source_identity",
        "source_identity_proved",
        "source_identity_verified",
        "actual_graph_owner_identity",
        "candidate_actual_integral_binder",
        "candidate_actual_integral_binder_proved",
        "source_callback_certified",
        "source_function_equality_proved",
        "identity_proved",
        "identity_verified",
        "same_binder_proved",
        "binder_identity",
        "binder_identity_proved",
        "binder_identity_verified",
        "binder_mapping_installed",
        "binder_mappings_installed",
        "m15",
        "m15_claimed",
        "parent",
        "parent_claimed",
        "actual_integral_value",
        "actual_integral_value_computed",
        "integral_value_computed",
        "numeric_values_computed",
        "numeric_proof",
        "source_value_evaluated",
        "source_values_evaluated",
        "integrals_evaluated",
        "graph_decoded",
        "graph_traversed",
        "graph_evaluated",
        "old_producer_imported",
        "certified",
        "complete_identity_proved",
        "parent_identity",
        "parent_identity_proved",
        "parent_owner_identity",
    }
)


def _false_marker(value: Any) -> bool:
    return value is False or value is None or (type(value) is int and value == 0)


def _reject_forbidden_claims(value: Any, path: str = "leaf_output") -> None:
    if isinstance(value, Mapping):
        for raw_key, item in value.items():
            key = str(raw_key).strip().lower().replace("-", "_")
            item_path = path + "." + str(raw_key)
            if key in _FORBIDDEN_TRUE_KEYS:
                _require(_false_marker(item), "forbidden claim at " + item_path)
            _reject_forbidden_claims(item, item_path)
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _reject_forbidden_claims(item, path + "[" + str(index) + "]")


def _validate_leaf_output(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "name-binding leaf result is not a mapping")
    result = dict(value)
    _require(result.get("schema") == LEAF_RESULT_SCHEMA, "name-binding leaf schema mismatch")
    _require(
        result.get("status") in {
            "conditional_candidates",
            "partial_conditional_candidates",
            "unresolved",
        },
        "name-binding leaf returned an unsupported status",
    )
    _require(
        type(result.get("source_name_binding_candidate")) is bool,
        "name-binding source candidate marker is missing",
    )
    records = result.get("records")
    _require(isinstance(records, list) and len(records) == 7, "name-binding leaf must return seven records")
    for index, row in enumerate(records):
        _require(isinstance(row, Mapping), f"name-binding leaf record {index} is not an object")
        for key in (
            "actual_graph_owner_identity",
            "candidate_actual_integral_binder",
            "candidate_actual_integral_binder_proved",
            "binder_identity_proved",
            "binder_mapping_installed",
            "source_callback_certified",
            "source_function_equality_proved",
            "m15_claimed",
            "parent_identity",
            "parent_identity_proved",
            "parent_owner_identity",
        ):
            if key in row:
                _require(_false_marker(row.get(key)), f"per-item forbidden claim at records[{index}].{key}")
        if "parent_identity_unproved" in row:
            _require(
                row.get("parent_identity_unproved") is True,
                f"per-item parent boundary is not conservative at records[{index}]",
            )
    _reject_forbidden_claims(result)
    return _json_safe(result)


def _call_leaf_once(
    leaf: ModuleType,
    saved_result: Mapping[str, Any],
    source_texts: Mapping[str, str],
    progress: dict[str, Any],
) -> tuple[str, Any]:
    api_name = "analyze"
    function = getattr(leaf, api_name, None)
    _require(callable(function), "new leaf analyze is missing or not callable")
    progress["leaf_calls"] = 1
    try:
        return api_name, function(saved_result, source_texts)
    except BaseException:
        progress["leaf_traceback"] = traceback.format_exc()
        raise


def _consume(
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    progress: dict[str, Any],
) -> dict[str, Any]:
    # No historical package/source module may be imported before this point.
    _reauthenticate_after_entry(hashes, source_sha, brief_sha)
    leaf, fixture = _load_new_modules()
    progress["new_modules_imported"] = True

    run_fixture = getattr(fixture, "run_fixture", None)
    _require(callable(run_fixture), "fixture.run_fixture is missing or not callable")
    progress["fixture_calls"] = 1
    try:
        fixture_raw = run_fixture(leaf)
    except BaseException:
        progress["fixture_traceback"] = traceback.format_exc()
        raise
    progress["fixture_raw"] = _json_safe(fixture_raw)
    fixture_result = _validate_fixture(fixture_raw)
    progress["fixture"] = fixture_result
    _require(not _forbidden_project_modules(), "historical module loaded during fixture")

    # This is the second raw authentication pass.  It occurs after the
    # fixture and before either prior gzip is decompressed.
    raw = _read_frozen_raw(hashes)
    progress["raw_gzip_authenticated_before_decode"] = True

    scoped_catalog_receipt = _decode_gzip_json(raw[PRIOR_SCOPED_REL], PRIOR_SCOPED_REL)
    progress["scoped_catalog_decoded"] = True
    occurrence_audit_receipt = _decode_gzip_json(raw[PRIOR_OCCURRENCE_REL], PRIOR_OCCURRENCE_REL)
    progress["occurrence_audit_decoded"] = True
    saved_result, saved_result_meta = _build_saved_result(
        scoped_catalog_receipt,
        occurrence_audit_receipt,
        PRIOR_SCOPED_REL,
        PRIOR_OCCURRENCE_REL,
        hashes[PRIOR_SCOPED_REL],
        hashes[PRIOR_OCCURRENCE_REL],
    )
    progress["saved_result"] = saved_result_meta

    # Source files are transported as text only.  No AST parse or source
    # lookup occurs here; the new leaf owns its one AST-consumer call.
    source_texts: dict[str, str] = {}
    for relative in SOURCE_REL:
        try:
            source_texts[SOURCE_TEXT_ALIASES[relative]] = raw[relative].decode("utf-8")
        except UnicodeDecodeError as exc:
            raise RunnerInputError("actual source is not UTF-8 text: " + relative) from exc
    source_hashes = {relative: hashes[relative] for relative in SOURCE_REL}
    progress["source_texts_supplied"] = True
    progress["source_text_paths"] = list(SOURCE_REL)
    progress["source_text_aliases"] = dict(SOURCE_TEXT_ALIASES)
    progress["source_text_hashes"] = dict(source_hashes)
    progress["source_ast_consumer_in_runner"] = False

    progress["leaf_calls"] = 0
    api_name, leaf_raw = _call_leaf_once(
        leaf, saved_result, source_texts, progress
    )
    progress["leaf_api"] = api_name
    progress["leaf_raw_output"] = _json_safe(leaf_raw)
    leaf_result = _validate_leaf_output(leaf_raw)
    progress["leaf_output"] = leaf_result

    return {
        "status": "returned",
        "passed": True,
        "conditional_source_name_binding_candidate": bool(
            leaf_result.get("source_name_binding_candidate") is True
        ),
        "source_name_binding_candidate": bool(
            leaf_result.get("source_name_binding_candidate") is True
        ),
        "fixture": fixture_result,
        "saved_result": progress["saved_result"],
        "source_text_paths": list(SOURCE_REL),
        "source_text_aliases": dict(SOURCE_TEXT_ALIASES),
        "source_text_hashes": dict(source_hashes),
        "leaf_api": api_name,
        "leaf_output": leaf_result,
        "fixture_calls": 1,
        "leaf_calls": 1,
        "raw_gzip_authenticated_before_decode": True,
        "prior_receipts_decoded_after_fixture_only": True,
        "source_ast_consumer_in_runner": False,
        "graph_decoded": False,
        "graph_traversed": False,
        "graph_evaluated": False,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "old_producer_imported": False,
        "numeric_values_computed": False,
        "actual_integral_value_computed": False,
        "integrals_evaluated": False,
        "binder_mapping_installed": False,
        "binder_identity_proved": False,
        "all_non_source_identity_claims_false": True,
        "per_item_parent_identity_proved": False,
        "per_item_parent_identity_unproved": True,
        "actual_owner_proved": False,
        "source_identity_verified": False,
        "m15_claimed": False,
        "parent_claimed": False,
        "certified": False,
        "traceback": None,
    }


def _base_result(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    checkpoint_key: str,
    entry_sha: str,
    progress: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "schema": COMPONENT + "_runner_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "source_checkpoint_path": SOURCE.relative_to(ROOT).as_posix(),
        "pending_brief_path": BRIEF.relative_to(ROOT.parent.parent).as_posix(),
        "pending_brief_sha256": brief_sha,
        "pending_brief_checkpoint_path": checkpoint_key,
        "source_checkpoint_sha256": source_sha,
        "input_hashes": dict(hashes),
        "raw_file_sha256": dict(hashes),
        "one_entry_only": True,
        "do_not_reenter": True,
        "brief_modified": False,
        "raw_inputs_authenticated_before_import": True,
        "raw_gzip_authenticated_before_decode": True,
        "prior_receipts_decoded_after_fixture_only": True,
        "graph_decoded": False,
        "graph_traversed": False,
        "graph_evaluated": False,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "old_producer_imported": False,
        "numeric_values_computed": False,
        "actual_integral_value_computed": False,
        "integrals_evaluated": False,
        "binder_mapping_installed": False,
        "binder_identity_proved": False,
        "all_non_source_identity_claims_false": True,
        "per_item_parent_identity_proved": False,
        "per_item_parent_identity_unproved": True,
        "actual_owner_proved": False,
        "source_identity_verified": False,
        "m15_claimed": False,
        "parent_claimed": False,
        "certified": False,
        "fixture_calls": int(progress.get("fixture_calls", 0)),
        "leaf_calls": int(progress.get("leaf_calls", 0)),
        "progress": _json_safe(progress),
        "source_manifest": _json_safe(source),
        "entry_sha256": entry_sha,
    }


def main() -> int:
    try:
        source, hashes, source_sha, brief_sha, checkpoint_key = _preflight()
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_runner_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed_preflight",
            "passed": False,
            "entry_written": False,
            "result_written": False,
            "error": type(exc).__name__ + ": " + str(exc),
            "traceback": traceback.format_exc(),
        }
        print(json.dumps(output, ensure_ascii=False, sort_keys=True))
        return 1

    progress: dict[str, Any] = {
        "fixture_calls": 0,
        "leaf_calls": 0,
        "new_modules_imported": False,
        "fixture": None,
        "fixture_raw": None,
        "saved_result": None,
        "scoped_catalog_decoded": False,
        "occurrence_audit_decoded": False,
        "source_texts_supplied": False,
        "source_text_paths": [],
        "source_text_hashes": {},
        "source_ast_consumer_in_runner": False,
        "raw_gzip_authenticated_before_decode": False,
        "prior_receipts_decoded_after_fixture_only": False,
        "leaf_output": None,
        "leaf_raw_output": None,
        "graph_decoded": False,
        "graph_traversed": False,
        "graph_evaluated": False,
        "old_producer_imported": False,
        "numeric_values_computed": False,
        "actual_integral_value_computed": False,
        "integrals_evaluated": False,
        "binder_mapping_installed": False,
        "binder_identity_proved": False,
        "actual_owner_proved": False,
        "source_identity_verified": False,
        "m15_claimed": False,
        "parent_claimed": False,
    }

    try:
        entry_sha = _write_exclusive(
            ENTRY,
            {
                "schema": COMPONENT + "_entered_v1",
                "component": COMPONENT,
                "phase": PHASE,
                "entered": True,
                "one_entry_only": True,
                "do_not_reenter": True,
                "source_checkpoint_sha256": source_sha,
                "pending_brief_sha256": brief_sha,
                "pending_brief_checkpoint_path": checkpoint_key,
                "input_hashes": dict(hashes),
                "raw_inputs_authenticated_before_import": True,
                "raw_gzip_authenticated_before_decode": True,
                "prior_receipts_decoded_after_fixture_only": True,
                "graph_decoded": False,
                "graph_traversed": False,
                "graph_evaluated": False,
                "old_producer_imported": False,
                "numeric_values_computed": False,
                "binder_mapping_installed": False,
                "fixture_calls": 0,
                "leaf_calls": 0,
            },
        )
    except BaseException as exc:
        output = {
            "schema": COMPONENT + "_runner_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed_entry",
            "passed": False,
            "entry_written": False,
            "result_written": False,
            "error": type(exc).__name__ + ": " + str(exc),
            "traceback": traceback.format_exc(),
        }
        print(json.dumps(output, ensure_ascii=False, sort_keys=True))
        return 1

    status = "returned"
    try:
        try:
            output = _consume(hashes, source_sha, brief_sha, progress)
            output.update({"status": "returned", "passed": True, "traceback": None})
        except BaseException as exc:
            status = "failed"
            error_traceback = traceback.format_exc()
            progress["traceback"] = error_traceback
            output = {
                "status": status,
                "passed": False,
                "conditional_source_name_binding_candidate": False,
                "source_name_binding_candidate": False,
                "fixture": progress.get("fixture"),
                "saved_result": progress.get("saved_result"),
                "source_text_paths": progress.get("source_text_paths"),
                "source_text_aliases": progress.get("source_text_aliases"),
                "source_text_hashes": progress.get("source_text_hashes"),
                "leaf_api": progress.get("leaf_api"),
                "leaf_output": progress.get("leaf_output"),
                "leaf_raw_output": progress.get("leaf_raw_output"),
                "error": type(exc).__name__ + ": " + str(exc),
                "traceback": error_traceback,
                "raw_gzip_authenticated_before_decode": bool(
                    progress.get("raw_gzip_authenticated_before_decode")
                ),
                "prior_receipts_decoded_after_fixture_only": bool(
                    progress.get("prior_receipts_decoded_after_fixture_only")
                ),
                "graph_decoded": False,
                "graph_traversed": False,
                "graph_evaluated": False,
                "old_producer_imported": False,
                "numeric_values_computed": False,
                "actual_integral_value_computed": False,
                "integrals_evaluated": False,
                "binder_mapping_installed": False,
                "binder_identity_proved": False,
                "actual_owner_proved": False,
                "source_identity_verified": False,
                "m15_claimed": False,
                "parent_claimed": False,
                "certified": False,
            }
    finally:
        _cleanup_new_modules()

    output.update(_base_result(source, hashes, source_sha, brief_sha, checkpoint_key, entry_sha, progress))
    try:
        result_sha = _write_exclusive(RESULT, output)
    except BaseException as exc:
        print(
            json.dumps(
                {
                    "component": COMPONENT,
                    "phase": PHASE,
                    "status": "failed_result_write",
                    "entry_sha256": entry_sha,
                    "error": type(exc).__name__ + ": " + str(exc),
                    "traceback": traceback.format_exc(),
                },
                ensure_ascii=False,
                sort_keys=True,
            )
        )
        return 1

    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": status,
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "source_checkpoint_sha256": source_sha,
                "pending_brief_sha256": brief_sha,
                "input_hashes": dict(hashes),
                "fixture_calls": int(progress.get("fixture_calls", 0)),
                "leaf_calls": int(progress.get("leaf_calls", 0)),
                "scoped_catalog_decoded": bool(progress.get("scoped_catalog_decoded")),
                "occurrence_audit_decoded": bool(progress.get("occurrence_audit_decoded")),
                "raw_gzip_authenticated_before_decode": bool(
                    progress.get("raw_gzip_authenticated_before_decode")
                ),
                "graph_decoded": False,
                "graph_traversed": False,
                "numeric_values_computed": False,
                "binder_mapping_installed": False,
                "actual_owner_proved": False,
                "source_identity_verified": False,
                "m15_claimed": False,
                "parent_claimed": False,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
