"""One-shot provenance-only bridge for the saved v3 green-scale result.

The runner does not re-run v3 AST mathematics and never imports the v3 leaf.
The v1/v2/v3 failure archives are raw-byte anchors only: they are hashed and
never decoded or decompressed.  The new fixture receives only fresh synthetic
mappings.  Only after that fixture succeeds are the saved v3 result and the
decoded cutoff parent passed to the new provenance leaf once.
"""

from __future__ import annotations

import hashlib
import importlib.machinery
import importlib.util
import json
import math
import os
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from types import ModuleType
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
COMPONENT = "n1_green_scale_provenance_bridge_v1"
PHASE = COMPONENT + "_pending"
PREVIOUS_COMPONENT = "n1_green_scale_commutator_v3"
PREVIOUS_PHASE = PREVIOUS_COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_MODULE = "openai_ns_reconstruction.n1_green_scale_provenance_bridge_v1"
FIXTURE_MODULE = "openai_ns_reconstruction.n1_green_scale_provenance_bridge_fixture_v1"
LEAF_REL = "src/openai_ns_reconstruction/n1_green_scale_provenance_bridge_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_green_scale_provenance_bridge_fixture_v1.py"
RUNNER_REL = "scripts/run_n1_green_scale_provenance_bridge_v1_once.py"
V3_RESULT_REL = "tmp/n1_green_scale_commutator_v3/result.json"
PARENT_CUTOFF_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
)
V1_FAILURE_REL = "experiments/root_st073/n1_green_scale_commutator_v1_failure.json.gz"
V2_FAILURE_REL = "experiments/root_st073/n1_green_scale_commutator_v2_failure.json.gz"
V3_FAILURE_REL = "experiments/root_st073/n1_green_scale_commutator_v3_failure.json.gz"
COMMON_SOURCE_DOC_REL = "docs/CURRENT_N1_COMMON_SOURCE_BRIDGE_2026_11_02.md"
RAW_BRIDGE_DOC_REL = "docs/CURRENT_N1_RAW_KCHI_RESIDUAL_BRIDGE_2026_11_02.md"

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    RUNNER_REL,
    V3_RESULT_REL,
    PARENT_CUTOFF_REL,
    V1_FAILURE_REL,
    V2_FAILURE_REL,
    V3_FAILURE_REL,
    COMMON_SOURCE_DOC_REL,
    RAW_BRIDGE_DOC_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)
FAILURE_ARCHIVES = (V1_FAILURE_REL, V2_FAILURE_REL, V3_FAILURE_REL)
PROVENANCE_DOCUMENTS = (COMMON_SOURCE_DOC_REL, RAW_BRIDGE_DOC_REL)
SHA256_HEX = frozenset("0123456789abcdefABCDEF")

PINNED_INPUTS = {
    V3_RESULT_REL: "87bc982493ca1f3dcce0183e5eeb368c439b1ddf57ec295296895fbe5363db95",
    PARENT_CUTOFF_REL: "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856",
    V1_FAILURE_REL: "d0ccd6e1ffd646eaa702568e9cae3d930acc440fa8b977fde4eaba48c068edf0",
    V2_FAILURE_REL: "719a57707fa5b478bcd4e3ff61905ff50f1fdb9b084983f388444db8be449bc4",
    V3_FAILURE_REL: "3b07192d99b21d8127775145f444e2d64767f0b2987a39cbbd5a32b47253393a",
    COMMON_SOURCE_DOC_REL: "9c1e3c5107bcf8baaa00a9d8e99a59fe3d2dee5414eb02d245d0af8ff4621593",
    RAW_BRIDGE_DOC_REL: "b9cab8da096a995336be8cdd7f8cf88b305337f51ade2aabaa5e0d2075cefdfb",
}

LEAF_API = "audit_saved_mappings(v3_result, parent_report)"
FIXTURE_API = "run_fixture(auditor)"
LEAF_CANDIDATES = (
    "audit_saved_mappings",
    "audit_provenance",
    "audit_saved_result",
    "audit_provenance_metadata",
    "bridge_provenance",
)
FORBIDDEN_MODULE_BASES = frozenset(
    {
        "original_n1_first_raw_function",
        "n1_kchi_first_correction_v1",
        "n1_raw_kchi_residual_bridge_v1",
        "n1_raw_kchi_residual_bridge_fixture_v1",
        "n1_green_scale_commutator_v1",
        "n1_green_scale_commutator_fixture_v1",
        "n1_green_scale_commutator_v2",
        "n1_green_scale_commutator_fixture_v2",
        "n1_green_scale_commutator_v3",
        "n1_green_scale_commutator_fixture_v3",
    }
)
FORBIDDEN_FILE_NAMES = frozenset(name + ".py" for name in FORBIDDEN_MODULE_BASES)


class RunnerInputError(ValueError):
    """Raised when the frozen provenance bridge contract is not satisfied."""


def _require(condition: bool, message: str) -> None:
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


def _forbidden_module_names() -> list[str]:
    found: list[str] = []
    for name, module in sys.modules.items():
        base = name.rsplit(".", 1)[-1]
        module_file = getattr(module, "__file__", None)
        file_name = Path(module_file).name if isinstance(module_file, str) else ""
        if base in FORBIDDEN_MODULE_BASES or file_name in FORBIDDEN_FILE_NAMES:
            found.append(name)
    return sorted(set(found))


def _json_safe(value: Any, seen: set[int] | None = None) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"type": "bytes", "bytes": len(value), "sha256": _sha_bytes(value)}
    active = seen if seen is not None else set()
    identity = id(value)
    if identity in active:
        return {"type": "cycle", "repr": repr(value)}
    active.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): _json_safe(item, active) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(item, active) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((_json_safe(item, active) for item in value), key=repr)
        return repr(value)
    finally:
        active.discard(identity)


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(_json_safe(value), ensure_ascii=False, sort_keys=True,
                   separators=(",", ":"), allow_nan=False).encode("utf-8")
        + b"\n"
    )


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


def _load_object(path: Path, label: str) -> dict[str, Any]:
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise RunnerInputError(label + " is unreadable") from exc
    return _decode_object(raw, label)


def _repo_file(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    normalized = raw_path.replace("\\", "/")
    relative = Path(normalized)
    _require(not relative.is_absolute() and ".." not in relative.parts,
             label + " path must be repository-relative")
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(),
             label + " path is missing or escapes repository: " + normalized)
    return key, path


def _sha_text(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64
             and all(character in SHA256_HEX for character in value),
             "invalid SHA-256 for " + label)
    return value.lower()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    files = source.get("files")
    _require(isinstance(files, list), "provenance source manifest files list is missing")
    hashes: dict[str, str] = {}
    for index, row in enumerate(files):
        _require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        key, path = _repo_file(row.get("path"), f"source.files[{index}]")
        _require(key in FROZEN_SET, "manifest contains out-of-scope file: " + key)
        _require(key not in hashes, "duplicate manifest path: " + key)
        expected = _sha_text(row.get("sha256"), key)
        actual = _sha_file(path)
        _require(actual == expected, "raw SHA mismatch: " + key)
        if key in PINNED_INPUTS:
            _require(actual == PINNED_INPUTS[key], "pinned input changed: " + key)
        if key in FAILURE_ARCHIVES:
            marker = row.get("raw_sha256_before_decode",
                            row.get("authenticate_raw_before_decode"))
            _require(marker is True, key + " is not marked raw-only")
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET,
             "manifest must freeze exactly the ten provenance inputs")
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _record_for(brief: Mapping[str, Any], component: str) -> Mapping[str, Any]:
    if brief.get("component") == component:
        return brief
    nested = brief.get(component)
    if isinstance(nested, Mapping):
        return nested
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping) and continuation.get("component") == component:
        return continuation
    return brief


def _brief_checkpoint(brief: Mapping[str, Any]) -> Mapping[str, Any]:
    checkpoint: Any = brief.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        nested = brief.get(COMPONENT)
        if isinstance(nested, Mapping):
            checkpoint = nested.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        continuation = brief.get("continuation")
        if isinstance(continuation, Mapping):
            checkpoint = continuation.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending provenance source checkpoint missing")
    return checkpoint


def _resolve_external(raw_path: Any, label: str) -> Path:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    path = Path(raw_path)
    if not path.is_absolute():
        path = ROOT / path
    path = path.resolve()
    _require(path.is_file(), label + " file is missing")
    return path


def _contract_path(contract: Mapping[str, Any], names: tuple[str, ...], label: str) -> str:
    for name in names:
        if isinstance(contract.get(name), str):
            return contract[name]
    raise RunnerInputError(label + " contract path is missing")


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "provenance component already entered")
    _require(BRIEF.is_file(), "resume brief is missing")
    _require(SOURCE.is_file(), "provenance source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    brief = _decode_object(brief_raw, "resume brief")
    record = _record_for(brief, COMPONENT)
    _require(record.get("component") == COMPONENT, "pending provenance brief component mismatch")
    _require(record.get("phase") == PHASE, "pending provenance brief phase mismatch")
    _require(record.get("one_entry_only") is True, "pending provenance one-entry guard missing")
    _require(record.get("do_not_reenter") is True, "pending provenance reentry guard missing")

    checkpoint = _brief_checkpoint(brief)
    checkpoint_path = _resolve_external(checkpoint.get("path"), "provenance source checkpoint")
    _require(checkpoint_path.resolve() == SOURCE.resolve(), "provenance brief points elsewhere")
    source_sha = _sha_file(SOURCE)
    _require(_sha_text(checkpoint.get("sha256"), "provenance checkpoint") == source_sha,
             "provenance source checkpoint SHA mismatch")
    preimage_path = _resolve_external(checkpoint.get("preimage_path"), "provenance brief preimage")
    preimage_raw = preimage_path.read_bytes()
    _require(_sha_text(checkpoint.get("preimage_sha256"), "provenance preimage")
             == _sha_bytes(preimage_raw), "provenance preimage SHA mismatch")
    previous = _record_for(_decode_object(preimage_raw, "provenance brief preimage"),
                           PREVIOUS_COMPONENT)
    _require(previous.get("component") == PREVIOUS_COMPONENT,
             "provenance preimage is not v3")
    _require(previous.get("phase") == PREVIOUS_PHASE,
             "provenance preimage is not v3 pending")

    source = _load_object(SOURCE, "provenance source manifest")
    _require(source.get("component") == COMPONENT, "provenance source component mismatch")
    _require(source.get("one_entry_only") is True, "provenance source one-entry guard missing")
    if "phase" in source:
        _require(source.get("phase") == PHASE, "provenance source phase mismatch")
    _require(source.get("previous_component") == PREVIOUS_COMPONENT,
             "provenance source previous component mismatch")
    contract = source.get("input_contract")
    _require(isinstance(contract, Mapping), "provenance input contract missing")
    _require(contract.get("leaf_api") == LEAF_API, "provenance leaf API contract mismatch")
    _require(contract.get("fixture_api") == FIXTURE_API, "provenance fixture API contract mismatch")
    _require(_contract_path(contract, ("saved_v3_result", "v3_result", "saved_result"),
                            "saved v3 result") == V3_RESULT_REL,
             "saved v3 result contract mismatch")
    _require(contract.get("parent_report") == PARENT_CUTOFF_REL,
             "provenance parent report contract mismatch")
    archive_contract = contract.get("failure_archives")
    if isinstance(archive_contract, (list, tuple)):
        archive_paths = {
            item if isinstance(item, str) else item.get("path")
            for item in archive_contract
            if isinstance(item, (str, Mapping))
        }
        _require(archive_paths == set(FAILURE_ARCHIVES),
                 "provenance failure archive contract mismatch")
    else:
        _require(contract.get("v1_failure_archive") == V1_FAILURE_REL,
                 "v1 failure archive contract mismatch")
        _require(contract.get("v2_failure_archive") == V2_FAILURE_REL,
                 "v2 failure archive contract mismatch")
        _require(contract.get("v3_failure_archive") == V3_FAILURE_REL,
                 "v3 failure archive contract mismatch")
    documents = contract.get("provenance_documents", contract.get("documents"))
    if isinstance(documents, (list, tuple)):
        document_paths = {
            item if isinstance(item, str) else item.get("path")
            for item in documents
            if isinstance(item, (str, Mapping))
        }
        _require(document_paths == set(PROVENANCE_DOCUMENTS),
                 "provenance document contract mismatch")
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, _sha_bytes(brief_raw), _sha_bytes(preimage_raw)


def _load_new_modules() -> tuple[Any, Any]:
    _require(not _forbidden_module_names(), "historical or v3 module already loaded")
    package_name = "openai_ns_reconstruction"
    package_dir = SRC / package_name
    previous_package = sys.modules.get(package_name)
    previous_leaf = sys.modules.get(LEAF_MODULE)
    previous_fixture = sys.modules.get(FIXTURE_MODULE)
    namespace = ModuleType(package_name)
    namespace.__path__ = [str(package_dir)]
    namespace.__package__ = package_name
    namespace.__spec__ = importlib.machinery.ModuleSpec(package_name, loader=None, is_package=True)
    sys.modules[package_name] = namespace
    old_path = list(sys.path)
    source_root = str(SRC)
    if source_root not in sys.path:
        sys.path.insert(0, source_root)
    try:
        leaf_spec = importlib.util.spec_from_file_location(LEAF_MODULE, ROOT / LEAF_REL)
        fixture_spec = importlib.util.spec_from_file_location(FIXTURE_MODULE, ROOT / FIXTURE_REL)
        _require(leaf_spec is not None and leaf_spec.loader is not None,
                 "provenance leaf spec unavailable")
        _require(fixture_spec is not None and fixture_spec.loader is not None,
                 "provenance fixture spec unavailable")
        leaf = importlib.util.module_from_spec(leaf_spec)
        sys.modules[LEAF_MODULE] = leaf
        leaf_spec.loader.exec_module(leaf)
        fixture = importlib.util.module_from_spec(fixture_spec)
        sys.modules[FIXTURE_MODULE] = fixture
        fixture_spec.loader.exec_module(fixture)
    finally:
        sys.path[:] = old_path
        if previous_package is None:
            sys.modules.pop(package_name, None)
        else:
            sys.modules[package_name] = previous_package
        if previous_leaf is None:
            sys.modules.pop(LEAF_MODULE, None)
        else:
            sys.modules[LEAF_MODULE] = previous_leaf
        if previous_fixture is None:
            sys.modules.pop(FIXTURE_MODULE, None)
        else:
            sys.modules[FIXTURE_MODULE] = previous_fixture
    _require(not _forbidden_module_names(), "new provenance import loaded forbidden module")
    return leaf, fixture


def _resolve_auditor(leaf: Any) -> tuple[Any, str]:
    for name in LEAF_CANDIDATES:
        candidate = getattr(leaf, name, None)
        if callable(candidate):
            return candidate, name
    raise RunnerInputError("provenance leaf exports no supported audit callable")


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "provenance fixture result is not an object")
    result = dict(value)
    _require(result.get("schema") == "n1_green_scale_provenance_bridge_fixture_v1",
             "provenance fixture schema mismatch")
    _require(result.get("passed") is True, "provenance synthetic fixture did not pass")
    checks = result.get("checks")
    _require(isinstance(checks, Mapping) and bool(checks), "provenance fixture checks missing")
    _require(all(item is True for item in checks.values()),
             "provenance fixture has a failed check")
    records = result.get("synthetic_records", result.get("cases", result.get("mappings")))
    _require(isinstance(records, Mapping) and bool(records),
             "provenance fixture has no synthetic mapping records")
    _require(len(records) == 7, "provenance fixture must expose seven new synthetic mappings")
    for flag in (
        "historical_modules_called", "historical_producer_imported", "old_modules_called",
        "v3_leaf_imported", "v3_result_consumed", "source_ast_reexecuted",
        "source_texts_consumed", "source_text_consumed", "source_AST_recomputed",
        "numeric_evaluation_performed", "v3_actual_leaf_called_here",
        "v1_or_v2_fixture_called",
        "v1_failure_decompressed", "v2_failure_decompressed", "v3_failure_decompressed",
        "v1_or_v2_synthetic_cases_reentered",
    ):
        if flag in result:
            _require(result[flag] is False, "provenance fixture guard failed: " + flag)
    if "synthetic_only" in result:
        _require(result["synthetic_only"] is True, "provenance fixture is not synthetic-only")
    return _json_safe(result)


def _validate_saved_v3(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "saved v3 result is not an object")
    result = dict(value)
    _require(result.get("schema") == "n1_green_scale_commutator_v3_result_v1",
             "saved v3 result schema mismatch")
    _require(result.get("component") == PREVIOUS_COMPONENT,
             "saved v3 result component mismatch")
    _require(result.get("actual_leaf_calls") == 1,
             "saved v3 result must record one actual v3 leaf call")
    _require(result.get("source_only") is True, "saved v3 result source-only guard missing")
    _require(result.get("actual_source_metadata_only") is False,
             "saved v3 metadata-only flag changed")
    actual = result.get("actual_result")
    _require(isinstance(actual, Mapping), "saved v3 actual leaf result is missing")
    _require(actual.get("schema") == "n1_green_scale_commutator_v3",
             "saved v3 actual result schema mismatch")
    checks = actual.get("checks")
    _require(isinstance(checks, Mapping) and bool(checks),
             "saved v3 actual checks are missing")
    false_keys = {str(key) for key, item in checks.items() if item is False}
    _require(false_keys <= {
        "parent_report_input_hashes_present",
        "parent_report_input_hashes_match",
        "parent_report_passed",
    }, "saved v3 result has a non-provenance source proof failure")
    for key in ("source_executed", "archive_source_executed", "old_modules_called"):
        if key in actual:
            _require(actual[key] is False, "saved v3 source-only guard changed: " + key)
    return _json_safe(result)


def _archive_scope(hashes: Mapping[str, str]) -> dict[str, Any]:
    return {
        str(index + 1): {
            "path": path, "sha256": hashes[path],
            "raw_sha_authenticated_before_entry": True,
            "raw_only": True, "decoded": False, "decompressed": False,
        }
        for index, path in enumerate(FAILURE_ARCHIVES)
    }


def _base_result(source: Mapping[str, Any], hashes: Mapping[str, str], source_sha: str,
                 brief_sha: str, preimage_sha: str, entry_sha: str,
                 calls: list[str], invocations: int) -> dict[str, Any]:
    return {
        "schema": f"{COMPONENT}_result_v1", "component": COMPONENT, "phase": PHASE,
        "source_checkpoint_path": "tmp/" + COMPONENT + "/source.json",
        "source_checkpoint_sha256": source_sha, "source_manifest_sha256": source_sha,
        "brief_path": str(BRIEF), "brief_sha256": brief_sha,
        "brief_preimage_sha256": preimage_sha, "entry_sha256": entry_sha,
        "input_hashes": dict(hashes), "raw_file_sha256": dict(hashes),
        "source_manifest": dict(source), "consumer_calls": list(calls),
        "verification_invocations": invocations, "actual_leaf_calls": 0,
        "source_only": True, "metadata_only": True,
        "actual_source_metadata_only": True, "v3_source_ast_reexecuted": False,
        "fixture_completed_before_saved_result_decode": False,
        "failure_archives_raw_only": True, "v3_result_decoded": False,
        "parent_report_decoded": False,
        "complete_kchi_proved": False, "ledger_modified": False,
        "full_goal_complete": False, "failure_archives": _archive_scope(hashes),
        "provenance_documents": {
            path: {"path": path, "sha256": hashes[path]}
            for path in PROVENANCE_DOCUMENTS
        },
        "runtime_model_metadata": None,
    }


def _consume(source: Mapping[str, Any], hashes: Mapping[str, str], source_sha: str,
             brief_sha: str, preimage_sha: str, entry_sha: str,
             progress: dict[str, Any]) -> dict[str, Any]:
    raw: dict[str, bytes] = {}
    for path in FROZEN_INPUTS:
        data = (ROOT / path).read_bytes()
        _require(_sha_bytes(data) == hashes[path], "input changed after provenance entry: " + path)
        raw[path] = data

    leaf, fixture = _load_new_modules()
    auditor, auditor_name = _resolve_auditor(leaf)
    progress["leaf_callable"] = auditor_name
    progress["historical_modules_imported"] = _forbidden_module_names()
    run_fixture = getattr(fixture, "run_fixture", None)
    _require(callable(run_fixture), "provenance fixture must export run_fixture(auditor)")
    progress["verification_invocations"] = 1
    progress["consumer_calls"] = ["fixture.run_fixture(audit_saved_mappings)"]
    # The fixture owns all mappings here; the saved v3 result and parent are
    # intentionally still opaque until this synthetic gate passes.
    fixture_output = run_fixture(auditor)
    progress["fixture_raw"] = _json_safe(fixture_output)
    fixture_result = _validate_fixture(fixture_output)
    progress["fixture"] = fixture_result
    progress["fixture_completed_before_saved_result_decode"] = True
    progress["failure_archives_raw_only"] = True
    _require(not _forbidden_module_names(), "fixture imported a forbidden module")

    saved_v3_result = _validate_saved_v3(_decode_object(raw[V3_RESULT_REL], "saved v3 result"))
    progress["saved_v3_result"] = saved_v3_result
    progress["v3_result_decoded"] = True
    parent_report = _decode_object(raw[PARENT_CUTOFF_REL], "parent cutoff report")
    progress["parent_report_decoded"] = True
    progress["saved_v3_result_reference"] = {
        "path": V3_RESULT_REL, "sha256": hashes[V3_RESULT_REL],
        "bytes": len(raw[V3_RESULT_REL]), "decoded_after_fixture": True,
    }
    progress["parent_report"] = {
        "path": PARENT_CUTOFF_REL, "sha256": hashes[PARENT_CUTOFF_REL],
        "bytes": len(raw[PARENT_CUTOFF_REL]), "schema": parent_report.get("schema"),
        "decoded_after_fixture": True,
    }
    progress["failure_archives"] = _archive_scope(hashes)

    progress["actual_leaf_calls"] = 1
    progress["consumer_calls"].append("leaf.audit_saved_mappings(v3_result, parent_report)")
    actual_output = auditor(saved_v3_result, parent_report)
    progress["verification_invocations"] = 2
    encoded_actual = _json_safe(actual_output)
    progress["actual_result"] = encoded_actual
    _require(isinstance(actual_output, Mapping), "provenance leaf result is not an object")
    _require(isinstance(encoded_actual, Mapping), "encoded provenance result is not an object")
    _require(not _forbidden_module_names(), "provenance leaf imported a forbidden module")
    _require(encoded_actual.get("schema") == "n1_green_scale_provenance_bridge_v1",
             "provenance leaf schema mismatch")
    _require(encoded_actual.get("passed") is True, "provenance bridge did not pass")
    checks = encoded_actual.get("checks")
    _require(isinstance(checks, Mapping) and bool(checks),
             "provenance bridge checks are missing")
    _require(all(item is True for item in checks.values()),
             "provenance bridge has a failed check")
    for flag in (
        "v3_source_ast_reexecuted", "v3_leaf_imported", "source_ast_rerun",
        "source_AST_recomputed", "source_text_consumed",
        "numeric_evaluation_performed", "source_provider_executed",
        "historical_producer_imported", "v1_failure_decompressed",
        "v2_failure_decompressed", "v3_failure_decompressed", "ledger_modified",
        "full_goal_complete", "complete_kchi_proved",
    ):
        if flag in encoded_actual:
            _require(encoded_actual[flag] is not True,
                     "provenance leaf made an out-of-scope claim: " + flag)
    if "metadata_only" in encoded_actual:
        _require(encoded_actual["metadata_only"] is True,
                 "provenance leaf is not metadata-only")
    return {
        "status": "returned", "passed": True, "fixture": fixture_result,
        "saved_v3_result": saved_v3_result,
        "saved_v3_result_reference": progress["saved_v3_result_reference"],
        "parent_report": progress["parent_report"], "actual_result": encoded_actual,
        "failure_archives": progress["failure_archives"], "actual_leaf_calls": 1,
        "ancestor_math_calls": 0, "source_only": True, "metadata_only": True,
        "actual_source_metadata_only": True, "v3_source_ast_reexecuted": False,
        "fixture_completed_before_saved_result_decode": True,
        "failure_archives_raw_only": True, "v3_result_decoded": True,
        "parent_report_decoded": True,
        "complete_kchi_proved": False, "ledger_modified": False,
        "full_goal_complete": False,
    }


def main() -> int:
    source, hashes, source_sha, brief_sha, preimage_sha = _preflight()
    entry_sha = _write_exclusive(ENTRY, {
        "schema": f"{COMPONENT}_entered_v1", "component": COMPONENT, "phase": PHASE,
        "entered": True, "source_checkpoint_sha256": source_sha,
        "source_manifest_sha256": source_sha, "brief_preimage_sha256": preimage_sha,
        "input_hashes": hashes, "raw_sha_authenticated_before_entry": True,
        "v1_failure_archive_decoded": False, "v1_failure_archive_decompressed": False,
        "v2_failure_archive_decoded": False, "v2_failure_archive_decompressed": False,
        "v3_failure_archive_decoded": False, "v3_failure_archive_decompressed": False,
        "v3_result_decoded": False, "parent_report_decoded": False,
        "consumer_calls_started": 0, "actual_leaf_calls": 0,
        "ancestor_math_calls": 0, "runtime_model_metadata": None,
    })
    progress: dict[str, Any] = {
        "fixture": None, "actual_result": None, "saved_v3_result": None,
        "actual_leaf_calls": 0, "consumer_calls": [], "verification_invocations": 0,
    }
    try:
        nested = _consume(source, hashes, source_sha, brief_sha, preimage_sha, entry_sha, progress)
        calls = list(progress.get("consumer_calls", ()))
        output = _base_result(source, hashes, source_sha, brief_sha, preimage_sha,
                              entry_sha, calls, int(progress.get("verification_invocations", 0)))
        output.update(nested)
        output["consumer_calls"] = calls
        output["verification_invocations"] = int(progress.get("verification_invocations", 0))
        output["actual_leaf_calls"] = int(progress.get("actual_leaf_calls", 0))
        output["error"] = None
        output["traceback"] = None
    except BaseException as exc:
        calls = list(progress.get("consumer_calls", ()))
        output = _base_result(source, hashes, source_sha, brief_sha, preimage_sha,
                              entry_sha, calls, int(progress.get("verification_invocations", 0)))
        output.update({
            "status": "failed", "passed": False,
            "fixture": progress.get("fixture") if progress.get("fixture") is not None
            else progress.get("fixture_raw"),
            "saved_v3_result": progress.get("saved_v3_result"),
            "actual_result": progress.get("actual_result"),
            "progress": _json_safe(progress),
            "actual_leaf_calls": int(progress.get("actual_leaf_calls", 0)),
            "error": {"type": type(exc).__name__, "message": str(exc)},
            "traceback": traceback.format_exc(),
        })
    result_sha = _write_exclusive(RESULT, output)
    print(json.dumps({
        "component": COMPONENT, "phase": PHASE, "status": output["status"],
        "result_sha256": result_sha, "consumer_calls": output["consumer_calls"],
        "verification_invocations": output["verification_invocations"],
        "actual_leaf_calls": output["actual_leaf_calls"],
    }, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
