"""Audit one already-saved failed occurrence result exactly once.

This runner authenticates the new leaf, its synthetic fixture, the previously
saved plain result, and the prior failure archive as raw bytes.  The failure
archive is raw-only and is never decompressed.  No graph, graph receipt, old
leaf, producer, witness traversal, numeric evaluation, or binder proof is
loaded or called.  The fixture is called once; only after it passes is the
saved JSON decoded and passed once to ``audit_saved_result``.

The audit output is retained as saved structural evidence.  It is not an
independent verification of the earlier graph traversal or of binder identity.
"""

from __future__ import annotations

import builtins
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

COMPONENT = "n1_integral_occurrence_saved_audit_v1"
PHASE = COMPONENT + "_pending"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/n1_integral_occurrence_saved_audit_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_integral_occurrence_saved_audit_fixture_v1.py"
RUNNER_REL = "scripts/run_n1_integral_occurrence_saved_audit_v1_once.py"
SAVED_RESULT_REL = "tmp/n1_integral_scoped_occurrence_v1/result.json"
FAILURE_ARCHIVE_REL = (
    "experiments/root_st073/n1_integral_scoped_occurrence_v1_failure.json.gz"
)

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    RUNNER_REL,
    SAVED_RESULT_REL,
    FAILURE_ARCHIVE_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)

SAVED_RESULT_SHA256 = (
    "67a47626c8dfdd3208848a5908072bf0d1d4389865b1184c3a5ec2dbfe9253ac"
)
FAILURE_ARCHIVE_SHA256 = (
    "c1308e07c2d13de25aa0f602b10eef1362c59f6a6b514ff3a821f888cd5ce7f6"
)

LEAF_SCHEMA = "n1_integral_occurrence_saved_audit_v1"
LEAF_RESULT_SCHEMA = LEAF_SCHEMA + "_result"
FIXTURE_SCHEMA = "n1_integral_occurrence_saved_audit_fixture_v1"

SHA256_HEX = frozenset("0123456789abcdefABCDEF")
_LOADED_NAMES = (
    "_n1_integral_occurrence_saved_audit_v1_runner_leaf",
    "_n1_integral_occurrence_saved_audit_v1_runner_fixture",
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
    """Keep complete JSON-safe output without evaluating opaque values."""

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
        return {
            "type": "opaque_value",
            "class": type(value).__name__,
            "repr": repr(value),
        }
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
        if key == SAVED_RESULT_REL:
            _require(
                row.get("raw_sha256_before_decode") is True,
                "saved result is not marked raw-pinned before decode",
            )
        if key == FAILURE_ARCHIVE_REL:
            _require(
                row.get("raw_sha256_before_decode") is True,
                "failure archive is not marked raw-pinned before decode",
            )
            _require(
                row.get("raw_sha256_only") is True,
                "failure archive is not marked raw-only",
            )
            for forbidden_key in ("decode", "decompressed", "decoded"):
                if forbidden_key in row:
                    _require(
                        row.get(forbidden_key) in (False, None, 0),
                        "failure archive decode marker is enabled: " + forbidden_key,
                    )
        hashes[key] = actual
    _require(
        frozenset(hashes) == FROZEN_SET,
        "manifest must contain exactly the five frozen inputs",
    )
    _require(
        hashes[SAVED_RESULT_REL] == SAVED_RESULT_SHA256,
        "saved result raw SHA is not the accepted hard pin",
    )
    _require(
        hashes[FAILURE_ARCHIVE_REL] == FAILURE_ARCHIVE_SHA256,
        "failure archive raw SHA is not the accepted hard pin",
    )
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _validate_contract(source: Mapping[str, Any]) -> None:
    contract = source.get("input_contract")
    _require(isinstance(contract, Mapping), "input_contract must be an object")
    expected = {
        "fixture_api": "run_fixture(leaf)",
        "leaf_api": "audit_saved_result(saved_result)",
        "graph_path": None,
        "prior_path": None,
        "saved_result_path": SAVED_RESULT_REL,
        "failure_archive_path": FAILURE_ARCHIVE_REL,
        "decode_after_fixture": True,
        "saved_result_decode_after_fixture": True,
        "failure_archive_raw_only": True,
        "no_graph_read": True,
        "no_graph_import": True,
        "no_old_leaf_call": True,
        "no_old_producer_import": True,
        "no_witness_reanalysis": True,
        "saved_evidence_not_independent_verification": True,
        "no_actual_graph_evaluation": True,
    }
    for key, value in expected.items():
        _require(key in contract, "input_contract key is missing: " + key)
        _require(contract.get(key) == value, key + " contract mismatch")
    for key, value in (
        ("saved_result_sha256", SAVED_RESULT_SHA256),
        ("failure_archive_sha256", FAILURE_ARCHIVE_SHA256),
    ):
        if key in contract:
            _require(contract.get(key) == value, key + " contract mismatch")


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
    _require(
        isinstance(raw_path, str) and raw_path,
        "pending source checkpoint path is missing",
    )
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
    """Recheck source, brief, and every frozen raw input before import."""

    source_raw = SOURCE.read_bytes()
    _require(
        _sha_bytes(source_raw) == source_sha,
        "source manifest changed after exclusive entry",
    )
    source = _decode_object(source_raw, "source manifest after entry")
    _require(source.get("component") == COMPONENT, "post-entry source component mismatch")
    _require(source.get("phase") == PHASE, "post-entry source phase mismatch")
    _require(
        source.get("one_entry_only") is True,
        "post-entry source one-entry guard missing",
    )
    _require(
        source.get("do_not_reenter") is True,
        "post-entry source re-entry guard missing",
    )
    _validate_contract(source)
    current_hashes = _manifest_hashes(source)
    _require(
        current_hashes == dict(hashes),
        "frozen input manifest changed after exclusive entry",
    )

    brief_raw = BRIEF.read_bytes()
    _require(
        _sha_bytes(brief_raw) == brief_sha,
        "pending brief changed after exclusive entry",
    )
    current_brief_sha, _ = _validate_pending_brief(source_sha)
    _require(
        current_brief_sha == brief_sha,
        "pending brief checkpoint changed after exclusive entry",
    )


def _rehash_failure_archive(expected: str) -> None:
    """Authenticate gzip bytes only; intentionally never decode or decompress."""

    path = ROOT / FAILURE_ARCHIVE_REL
    _require(
        _sha_file(path) == expected == FAILURE_ARCHIVE_SHA256,
        "failure archive raw SHA changed or is not the accepted hard pin",
    )


def _read_saved_result_after_fixture(expected: str) -> tuple[dict[str, Any], str]:
    """Read and decode the plain saved JSON only after fixture success."""

    path = ROOT / SAVED_RESULT_REL
    raw = path.read_bytes()
    actual = _sha_bytes(raw)
    _require(
        actual == expected == SAVED_RESULT_SHA256,
        "saved result raw SHA changed or is not the accepted hard pin",
    )
    return _decode_object(raw, SAVED_RESULT_REL), actual


def _forbidden_project_modules() -> list[str]:
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
        return original(name, globals_, locals_, fromlist, level)

    return checked


def _cleanup_new_modules() -> None:
    for name in _LOADED_NAMES:
        sys.modules.pop(name, None)


def _load_module(relative: str, module_name: str) -> ModuleType:
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(module_name, path)
    _require(
        spec is not None and spec.loader is not None,
        "new module spec is unavailable: " + relative,
    )
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
        _require(
            not _forbidden_project_modules(),
            "historical module loaded by new files",
        )
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
    _require(
        all(item is True for item in checks.values()),
        "synthetic fixture contains a failed check",
    )
    for key in (
        "graph_read",
        "graph_traversed",
        "source_traversal_reexecuted",
        "numeric_values_computed",
        "integrals_evaluated",
        "binder_identity_proved",
        "complete_identity_proved",
        "actual_graph_evaluate_called",
        "actual_integral_value_computed",
        "source_graph_imported",
        "source_graph_executed",
        "binder_mapping_installed",
        "cancellation_claimed",
    ):
        if key in result:
            _require(
                result[key] in (False, 0, None),
                "fixture reported forbidden activity: " + key,
            )
    return _json_safe(result)


_FORBIDDEN_TRUE_KEYS = frozenset(
    {
        "candidate_actual_integral_binder",
        "candidate_actual_integral_binder_proved",
        "candidate_actual_integral_binder_verified",
        "binder_identity",
        "binder_identity_proved",
        "binder_identity_verified",
        "complete_identity_proved",
        "complete_identity_verified",
        "same_binder_proved",
        "same_binder_verified",
        "numeric_proof",
        "certified",
        "integral_value_computed",
        "actual_integral_value_computed",
        "source_value_computed",
        "source_values_computed",
        "numeric_values_computed",
        "actual_point_values_computed",
        "numeric_equivalence_certified",
        "m15_claimed",
        "kchi_claimed",
        "cancellation",
        "cancellation_claimed",
    }
)
_FORBIDDEN_FALSE_ONLY_KEYS = _FORBIDDEN_TRUE_KEYS | frozenset(
    {
        "numeric_field_value_computed",
        "source_value_evaluated",
        "integrals_evaluated",
        "binder_mappings_installed",
        "actual_graph_evaluate_called",
        "graph_evaluated",
        "source_values_evaluated",
        "saved_evidence_independent_verification",
        "saved_evidence_reanalyzed",
        "witness_reanalyzed",
    }
)


def _reject_forbidden_claims(value: Any, path: str = "audit_result") -> None:
    if isinstance(value, Mapping):
        for raw_key, item in value.items():
            key = str(raw_key).strip().lower().replace("-", "_")
            item_path = path + "." + str(raw_key)
            if key in _FORBIDDEN_FALSE_ONLY_KEYS:
                allowed_false = item is False or item is None or (
                    type(item) is int and item == 0
                )
                if not allowed_false:
                    raise RunnerInputError("forbidden claim at " + item_path)
            _reject_forbidden_claims(item, item_path)
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _reject_forbidden_claims(item, path + "[" + str(index) + "]")


def _validate_audit_result(value: Any) -> dict[str, Any]:
    """Validate the saved-only audit envelope without reanalyzing witnesses."""

    _require(isinstance(value, Mapping), "saved audit result is not a mapping")
    result = dict(value)
    _require(result.get("schema") == LEAF_RESULT_SCHEMA, "saved audit schema mismatch")
    _require(
        result.get("status") == "validated_saved_failure",
        "saved audit did not validate the pinned failed result",
    )
    _require(result.get("saved_status") == "failed", "saved audit status mismatch")
    _require(result.get("fixture_passed") is True, "saved audit fixture marker mismatch")
    _require(result.get("leaf_calls") == 1, "saved audit leaf call count is not one")
    _require(
        result.get("catalog_record_count") == 7,
        "saved audit catalog record count is not seven",
    )
    saved_records = result.get("saved_leaf_records")
    catalog_records = result.get("saved_catalog_records")
    _require(
        isinstance(saved_records, list) and len(saved_records) == 7,
        "saved audit does not preserve seven leaf records",
    )
    _require(
        isinstance(catalog_records, list) and len(catalog_records) == 7,
        "saved audit does not preserve seven catalog records",
    )
    for key in (
        "source_traversal_reexecuted",
        "graph_read",
        "graph_traversed",
        "paths_independently_revalidated",
        "scopes_independently_revalidated",
        "candidate_actual_integral_binder_proved",
        "binder_identity_proved",
        "complete_identity_proved",
        "numeric_field_value_computed",
        "source_value_evaluated",
        "integrals_evaluated",
        "binder_mappings_installed",
    ):
        _require(result.get(key) is False, "saved audit boundary is enabled: " + key)
    _require(
        result.get("structural_occurrence_is_not_identity") is True,
        "saved audit did not preserve the non-identity boundary",
    )
    claims = result.get("claims")
    _require(isinstance(claims, Mapping), "saved audit claims are missing")
    for key in (
        "saved_records_only",
        "structural_occurrence_is_not_identity",
    ):
        _require(claims.get(key) is True, "saved audit claim is missing: " + key)
    for key in (
        "independent_paths_revalidated",
        "complete_proof",
        "candidate_actual_integral_binder",
        "numeric_values",
        "integral_values",
    ):
        _require(claims.get(key) is False, "saved audit claim is enabled: " + key)
    _reject_forbidden_claims(result)
    return _json_safe(result)


def _call_audit_once(leaf: ModuleType, saved_result: Mapping[str, Any]) -> tuple[str, Any]:
    name = "audit_saved_result"
    function = getattr(leaf, name, None)
    _require(callable(function), "new leaf audit_saved_result is missing or not callable")
    return name, function(saved_result)


def _consume(
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    progress: dict[str, Any],
) -> dict[str, Any]:
    # This reauthentication is before importing either new file.  It hashes
    # every frozen input, but it does not decode the saved JSON or archive.
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

    # The failure archive remains raw-only.  This check deliberately performs
    # no gzip import, decompression, JSON decode, or interpretation.
    _rehash_failure_archive(hashes[FAILURE_ARCHIVE_REL])
    progress["failure_archive_raw_sha256"] = hashes[FAILURE_ARCHIVE_REL]
    progress["failure_archive_raw_only"] = True
    progress["failure_archive_decoded"] = False
    progress["failure_archive_decompressed"] = False

    # This is the sole saved-result read/decode, reachable only after fixture
    # success.  It is prior evidence, not a new graph or witness traversal.
    progress["saved_result_raw_sha256"] = hashes[SAVED_RESULT_REL]
    saved_result, saved_sha = _read_saved_result_after_fixture(hashes[SAVED_RESULT_REL])
    progress["saved_result_raw_sha256"] = saved_sha
    progress["saved_result_decoded"] = True
    progress["saved_result"] = _json_safe(saved_result)
    progress["saved_evidence_independent_verification"] = False
    progress["witness_reanalyzed"] = False

    progress["leaf_calls"] = 1
    try:
        api_name, audit_raw = _call_audit_once(leaf, saved_result)
    except BaseException:
        progress["leaf_traceback"] = traceback.format_exc()
        raise
    progress["leaf_api"] = api_name
    progress["leaf_raw_output"] = _json_safe(audit_raw)
    audit_result = _validate_audit_result(audit_raw)
    progress["leaf_output"] = audit_result
    return {
        "status": "returned",
        "passed": True,
        "fixture": fixture_result,
        "fixture_raw": _json_safe(fixture_raw),
        "leaf_api": api_name,
        "leaf_output": audit_result,
        "leaf_raw_output": _json_safe(audit_raw),
        "saved_result": _json_safe(saved_result),
        "saved_result_raw_sha256": saved_sha,
        "failure_archive_path": FAILURE_ARCHIVE_REL,
        "failure_archive_raw_sha256": hashes[FAILURE_ARCHIVE_REL],
        "failure_archive_raw_only": True,
        "failure_archive_decoded": False,
        "failure_archive_decompressed": False,
        "fixture_calls": 1,
        "leaf_calls": 1,
        "saved_record_count": 7,
        "graph_read": False,
        "graph_decoded": False,
        "graph_traversed": False,
        "graph_edges_traversed": False,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "old_leaf_called": False,
        "old_producer_imported": False,
        "witness_reanalyzed": False,
        "saved_evidence_independent_verification": False,
        "actual_graph_evaluate_called": False,
        "actual_integral_value_computed": False,
        "source_value_evaluated": False,
        "integrals_evaluated": False,
        "numeric_values_computed": False,
        "numeric_proof": False,
        "binder_mapping_installed": False,
        "candidate_actual_integral_binder_proved": False,
        "binder_identity_proved": False,
        "complete_identity_proved": False,
        "cancellation_claimed": False,
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
        "saved_result_decoded_after_fixture_only": True,
        "failure_archive_raw_only": True,
        "failure_archive_decoded": False,
        "failure_archive_decompressed": False,
        "graph_read": False,
        "graph_decoded": False,
        "graph_traversed": False,
        "graph_edges_traversed": False,
        "source_graph_imported": False,
        "source_graph_executed": False,
        "old_leaf_called": False,
        "old_producer_imported": False,
        "witness_reanalyzed": False,
        "saved_evidence_independent_verification": False,
        "actual_graph_evaluate_called": False,
        "actual_integral_value_computed": False,
        "source_value_evaluated": False,
        "integrals_evaluated": False,
        "numeric_values_computed": False,
        "numeric_proof": False,
        "binder_mapping_installed": False,
        "candidate_actual_integral_binder_proved": False,
        "binder_identity_proved": False,
        "complete_identity_proved": False,
        "cancellation_claimed": False,
        "certified": False,
        "fixture_calls": int(progress.get("fixture_calls", 0)),
        "leaf_calls": int(progress.get("leaf_calls", 0)),
        "saved_result_raw_sha256": progress.get("saved_result_raw_sha256"),
        "failure_archive_raw_sha256": hashes[FAILURE_ARCHIVE_REL],
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
        "saved_result_raw_sha256": None,
        "saved_result_decoded": False,
        "failure_archive_raw_sha256": hashes[FAILURE_ARCHIVE_REL],
        "failure_archive_raw_only": True,
        "failure_archive_decoded": False,
        "failure_archive_decompressed": False,
        "leaf_output": None,
        "leaf_raw_output": None,
        "graph_read": False,
        "graph_decoded": False,
        "graph_traversed": False,
        "graph_edges_traversed": False,
        "saved_evidence_independent_verification": False,
        "witness_reanalyzed": False,
        "actual_graph_evaluate_called": False,
        "actual_integral_value_computed": False,
        "numeric_values_computed": False,
        "binder_identity_proved": False,
        "complete_identity_proved": False,
        "cancellation_claimed": False,
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
                "saved_result_decoded_after_fixture_only": True,
                "failure_archive_raw_only": True,
                "failure_archive_decoded": False,
                "failure_archive_decompressed": False,
                "graph_read": False,
                "graph_traversed": False,
                "old_leaf_called": False,
                "old_producer_imported": False,
                "witness_reanalyzed": False,
                "saved_evidence_independent_verification": False,
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
                "fixture": progress.get("fixture"),
                "fixture_raw": progress.get("fixture_raw"),
                "saved_result": progress.get("saved_result"),
                "saved_result_raw_sha256": progress.get("saved_result_raw_sha256"),
                "leaf_api": progress.get("leaf_api"),
                "leaf_output": progress.get("leaf_output"),
                "leaf_raw_output": progress.get("leaf_raw_output"),
                "error": type(exc).__name__ + ": " + str(exc),
                "traceback": error_traceback,
                "failure_archive_path": FAILURE_ARCHIVE_REL,
                "failure_archive_raw_sha256": hashes[FAILURE_ARCHIVE_REL],
                "failure_archive_raw_only": True,
                "failure_archive_decoded": False,
                "failure_archive_decompressed": False,
                "graph_read": False,
                "graph_decoded": False,
                "graph_traversed": False,
                "graph_edges_traversed": False,
                "old_leaf_called": False,
                "old_producer_imported": False,
                "witness_reanalyzed": False,
                "saved_evidence_independent_verification": False,
                "actual_graph_evaluate_called": False,
                "actual_integral_value_computed": False,
                "numeric_values_computed": False,
                "binder_mapping_installed": False,
                "binder_identity_proved": False,
                "complete_identity_proved": False,
                "cancellation_claimed": False,
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
                "saved_result_decoded": bool(progress.get("saved_result_decoded")),
                "failure_archive_raw_only": True,
                "failure_archive_decoded": False,
                "graph_read": False,
                "graph_traversed": False,
                "saved_evidence_independent_verification": False,
                "witness_reanalyzed": False,
            },
            ensure_ascii=False,
            sort_keys=True,
        )
    )
    return 0 if status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
