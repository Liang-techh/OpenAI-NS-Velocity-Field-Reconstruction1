"""One-shot v2 source-bound entry for the N1 green-scale commutator.

This is a fresh leaf/fixture entry after the v1 synthetic-fixture failure.
All manifest inputs are authenticated by raw SHA before the exclusive entry.
The preserved v1 failure gzip is hashed only; it is never decompressed or
decoded.  The v2 fixture owns its three synthetic negative cases.  After the
fixture, the v2 leaf receives the three authenticated source texts and the
decoded parent cutoff report exactly once through ``audit_source_texts``.
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
COMPONENT = "n1_green_scale_commutator_v2"
PHASE = COMPONENT + "_pending"
PREVIOUS_COMPONENT = "n1_green_scale_commutator_v1"
PREVIOUS_PHASE = PREVIOUS_COMPONENT + "_pending"

BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"
PREVIOUS_SOURCE = ROOT / "tmp" / PREVIOUS_COMPONENT / "source.json"

LEAF_MODULE = "openai_ns_reconstruction.n1_green_scale_commutator_v2"
FIXTURE_MODULE = "openai_ns_reconstruction.n1_green_scale_commutator_fixture_v2"

LEAF_REL = "src/openai_ns_reconstruction/n1_green_scale_commutator_v2.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_green_scale_commutator_fixture_v2.py"
RUNNER_REL = "scripts/run_n1_green_scale_commutator_v2_once.py"
INNER_VOLTERRA_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
REGULAR_SERIES_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_regular_series_function.py"
)
ORIGINAL_RAW_REL = "src/openai_ns_reconstruction/original_n1_first_raw_function.py"
PARENT_CUTOFF_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cutoff_source_binding.json"
)
BRIDGE_DOC_REL = "docs/CURRENT_N1_RAW_KCHI_RESIDUAL_BRIDGE_2026_11_02.md"
V1_FAILURE_REL = "experiments/root_st073/n1_green_scale_commutator_v1_failure.json.gz"

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    RUNNER_REL,
    INNER_VOLTERRA_REL,
    REGULAR_SERIES_REL,
    ORIGINAL_RAW_REL,
    PARENT_CUTOFF_REL,
    BRIDGE_DOC_REL,
    V1_FAILURE_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)
SOURCE_TEXT_INPUTS = (INNER_VOLTERRA_REL, REGULAR_SERIES_REL, ORIGINAL_RAW_REL)
PROVENANCE_INPUTS = (BRIDGE_DOC_REL, V1_FAILURE_REL)
SOURCE_TEXT_ROLES = ("volterra", "regular", "first_raw")
LEAF_API = (
    "audit_source_texts(volterra_text, regular_text, first_raw_text, "
    "*, parent, source_labels)"
)
FIXTURE_API = "run_fixture(auditor)"
SHA256_HEX = frozenset("0123456789abcdefABCDEF")

PINNED_INPUTS = {
    INNER_VOLTERRA_REL: "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208",
    REGULAR_SERIES_REL: "b64b7cae8529cfc5caa158fe06ac5660ec1d2217fc4d845d02002dac3a195555",
    ORIGINAL_RAW_REL: "f13f29b4d0652e8ed3b3347f19bab7eb1ea0e18fd4be0f50cf6d550a1aeb92c6",
    PARENT_CUTOFF_REL: "4c36c2dc2e709d7ca740a689af9c9e16802569f300dfe9997851b9bb72142856",
    BRIDGE_DOC_REL: "b9cab8da096a995336be8cdd7f8cf88b305337f51ade2aabaa5e0d2075cefdfb",
    V1_FAILURE_REL: "d0ccd6e1ffd646eaa702568e9cae3d930acc440fa8b977fde4eaba48c068edf0",
}

FORBIDDEN_MODULE_BASES = frozenset(
    {
        "original_n1_first_raw_function",
        "n1_kchi_first_correction_v1",
        "n1_raw_kchi_residual_bridge_v1",
        "n1_raw_kchi_residual_bridge_fixture_v1",
        "n1_green_scale_commutator_v1",
        "n1_green_scale_commutator_fixture_v1",
        "1_green_scale_commutator_v1",
    }
)
FORBIDDEN_FILE_NAMES = frozenset(name + ".py" for name in FORBIDDEN_MODULE_BASES)


class RunnerInputError(ValueError):
    """Raised when the frozen v2 one-entry contract is not satisfied."""


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
        json.dumps(
            _json_safe(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
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
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " path must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " path is missing or escapes repository: " + normalized,
    )
    return key, path


def _sha_text(value: Any, label: str) -> str:
    _require(
        isinstance(value, str)
        and len(value) == 64
        and all(character in SHA256_HEX for character in value),
        "invalid SHA-256 for " + label,
    )
    return value.lower()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    entries: list[tuple[Any, Any, Mapping[str, Any] | None, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            entries.append((item.get("path"), item.get("sha256"), item, f"source.files[{index}]"))
    else:
        for field in ("raw_file_sha256", "input_sha256", "input_hashes", "sha256"):
            table = source.get(field)
            if isinstance(table, Mapping):
                entries.extend((path, digest, None, f"source.{field}") for path, digest in table.items())
                break
    _require(entries, "source manifest contains no raw-file SHA entries")

    hashes: dict[str, str] = {}
    for raw_path, expected, row, label in entries:
        key, path = _repo_file(raw_path, label)
        _require(key in FROZEN_SET, "manifest contains out-of-scope file: " + key)
        digest = _sha_text(expected, key)
        _require(key not in hashes, "duplicate manifest path: " + key)
        actual = _sha_file(path)
        _require(actual == digest, "raw SHA mismatch: " + key)
        if key in PINNED_INPUTS:
            _require(actual == PINNED_INPUTS[key], "pinned input changed: " + key)
        if key == V1_FAILURE_REL and row is not None:
            marker = row.get("raw_sha256_before_decode", row.get("authenticate_raw_before_decode"))
            _require(marker is True, "v1 failure archive is not marked raw-only")
        hashes[key] = actual

    _require(frozenset(hashes) == FROZEN_SET, "manifest must freeze exactly the nine required files")
    return {key: hashes[key] for key in FROZEN_INPUTS}


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
    _require(isinstance(checkpoint, Mapping), "pending v2 source checkpoint missing")
    return checkpoint


def _resolve_external(raw_path: Any, label: str) -> Path:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    path = Path(raw_path)
    if not path.is_absolute():
        path = ROOT / path
    path = path.resolve()
    _require(path.is_file(), label + " file is missing")
    return path


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "v2 component already entered")
    _require(BRIEF.is_file(), "resume brief is missing")
    _require(SOURCE.is_file(), "v2 source manifest is missing")

    brief_raw = BRIEF.read_bytes()
    brief = _decode_object(brief_raw, "resume brief")
    record = _brief_record(brief)
    _require(record.get("component") == COMPONENT, "pending v2 brief component mismatch")
    _require(record.get("phase") == PHASE, "pending v2 brief phase mismatch")
    _require(record.get("one_entry_only") is True, "pending v2 brief one-entry guard missing")
    _require(record.get("do_not_reenter") is True, "pending v2 brief reentry guard missing")

    checkpoint = _brief_checkpoint(brief)
    raw_path = checkpoint.get("path")
    checkpoint_path = _resolve_external(raw_path, "v2 source checkpoint")
    _require(checkpoint_path.resolve() == SOURCE.resolve(), "v2 brief points to another source")
    source_sha = _sha_file(SOURCE)
    _require(_sha_text(checkpoint.get("sha256"), "v2 source checkpoint") == source_sha,
             "v2 source checkpoint SHA mismatch")

    preimage_path = _resolve_external(checkpoint.get("preimage_path"), "v2 brief preimage")
    preimage_raw = preimage_path.read_bytes()
    _require(
        _sha_text(checkpoint.get("preimage_sha256"), "v2 brief preimage") == _sha_bytes(preimage_raw),
        "v2 brief preimage SHA mismatch",
    )
    preimage = _decode_object(preimage_raw, "v2 brief preimage")
    previous_record = _brief_record(preimage)
    _require(previous_record.get("component") == PREVIOUS_COMPONENT,
             "v2 brief preimage is not the v1 brief")
    _require(previous_record.get("phase") == PREVIOUS_PHASE,
             "v2 brief preimage is not v1 pending")

    source = _load_object(SOURCE, "v2 source manifest")
    _require(source.get("component") == COMPONENT, "v2 source manifest component mismatch")
    _require(source.get("one_entry_only") is True, "v2 source manifest one-entry guard missing")
    if "phase" in source:
        _require(source.get("phase") == PHASE, "v2 source manifest phase mismatch")
    _require(source.get("previous_component") == PREVIOUS_COMPONENT,
             "v2 source manifest previous component mismatch")

    previous_checkpoint = source.get("previous_source_checkpoint")
    if isinstance(previous_checkpoint, Mapping):
        previous_path = _resolve_external(previous_checkpoint.get("path"), "previous source checkpoint")
        _require(previous_path.resolve() == PREVIOUS_SOURCE.resolve(),
                 "v2 previous source checkpoint path mismatch")
        _require(
            _sha_text(previous_checkpoint.get("sha256"), "previous source checkpoint")
            == _sha_file(PREVIOUS_SOURCE),
            "v2 previous source checkpoint SHA mismatch",
        )

    contract = source.get("input_contract")
    _require(isinstance(contract, Mapping), "v2 input contract missing")
    _require(contract.get("leaf_api") == LEAF_API, "v2 leaf API contract mismatch")
    _require(contract.get("fixture_api") == FIXTURE_API, "v2 fixture API contract mismatch")
    _require(tuple(contract.get("source_texts", ())) in (SOURCE_TEXT_ROLES, SOURCE_TEXT_INPUTS),
             "v2 source-text contract mismatch")
    _require(contract.get("parent_report") == PARENT_CUTOFF_REL,
             "v2 parent report contract mismatch")
    _require(contract.get("v1_failure_archive") == V1_FAILURE_REL,
             "v2 failure archive contract mismatch")

    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, _sha_bytes(brief_raw), _sha_bytes(preimage_raw)


def _load_new_modules() -> tuple[Any, Any]:
    preexisting = _forbidden_module_names()
    _require(not preexisting, "historical or v1 module already loaded")

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
        _require(leaf_spec is not None and leaf_spec.loader is not None, "v2 leaf spec unavailable")
        _require(fixture_spec is not None and fixture_spec.loader is not None, "v2 fixture spec unavailable")
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
    _require(not _forbidden_module_names(), "v2 import loaded historical or v1 module")
    return leaf, fixture


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "v2 fixture result is not an object")
    result = dict(value)
    _require(result.get("schema") == "n1_green_scale_commutator_fixture_v2",
             "v2 fixture schema mismatch")
    _require(result.get("passed") is True, "v2 synthetic fixture did not pass")
    checks = result.get("checks")
    _require(isinstance(checks, Mapping) and bool(checks), "v2 fixture checks are missing")
    _require(all(item is True for item in checks.values()), "v2 fixture has a failed check")
    negatives = result.get("negative_cases", result.get("cases"))
    if isinstance(negatives, Mapping):
        negative_count = len(negatives)
    elif isinstance(negatives, (list, tuple)):
        negative_count = len(negatives)
    else:
        negative_count = result.get("negative_case_count")
    _require(negative_count == 3, "v2 fixture must expose three negative cases")
    for flag in (
        "historical_modules_called",
        "historical_producer_imported",
        "old_modules_called",
        "implementation_imported",
        "provider_executed",
        "numerical_residual_evaluated",
        "mpmath_used",
        "v1_cases_reused",
        "reused_v1_synthetic_cases",
        "v1_fixture_called",
        "v1_canonical_case_called",
        "v1_old_mutation_cases_called",
    ):
        if flag in result:
            _require(result[flag] is False, "v2 fixture guard failed: " + flag)
    return _json_safe(result)


def _actual_scope(hashes: Mapping[str, str]) -> dict[str, Any]:
    return {
        "source_only": True,
        "actual_source_metadata_only": False,
        "source_text_inputs": [
            {"path": path, "sha256": hashes[path], "passed_verbatim": True}
            for path in SOURCE_TEXT_INPUTS
        ],
        "parent_report_input": {
            "path": PARENT_CUTOFF_REL,
            "sha256": hashes[PARENT_CUTOFF_REL],
            "raw_sha_authenticated_before_entry": True,
            "decoded_after_entry_only": True,
        },
        "v1_failure_archive": {
            "path": V1_FAILURE_REL,
            "sha256": hashes[V1_FAILURE_REL],
            "raw_sha_authenticated_before_entry": True,
            "decoded": False,
            "decompressed": False,
        },
        "authenticated_bridge_doc": {"path": BRIDGE_DOC_REL, "sha256": hashes[BRIDGE_DOC_REL]},
        "complete_kchi_proved": False,
        "ledger_modified": False,
        "full_goal_complete": False,
        "ancestor_math_calls": 0,
        "historical_producer_imported": bool(_forbidden_module_names()),
    }


def _base_result(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    preimage_sha: str,
    entry_sha: str,
    calls: list[str],
    verification_invocations: int,
) -> dict[str, Any]:
    return {
        "schema": f"{COMPONENT}_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "source_checkpoint_path": "tmp/" + COMPONENT + "/source.json",
        "source_checkpoint_sha256": source_sha,
        "source_manifest_sha256": source_sha,
        "brief_path": str(BRIEF),
        "brief_sha256": brief_sha,
        "brief_preimage_sha256": preimage_sha,
        "entry_sha256": entry_sha,
        "input_hashes": dict(hashes),
        "raw_file_sha256": dict(hashes),
        "source_manifest": dict(source),
        "consumer_calls": list(calls),
        "verification_invocations": verification_invocations,
        "ancestor_math_calls": 0,
        "historical_producer_imported": bool(_forbidden_module_names()),
        "source_only": True,
        "actual_source_metadata_only": False,
        "complete_kchi_proved": False,
        "ledger_modified": False,
        "full_goal_complete": False,
        "actual_source_scope": _actual_scope(hashes),
        "runtime_model_metadata": None,
    }


def _consume(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    preimage_sha: str,
    entry_sha: str,
    progress: dict[str, Any],
) -> dict[str, Any]:
    raw: dict[str, bytes] = {}
    for path in FROZEN_INPUTS:
        data = (ROOT / path).read_bytes()
        _require(_sha_bytes(data) == hashes[path], "input changed after v2 entry: " + path)
        raw[path] = data

    leaf, fixture = _load_new_modules()
    progress["historical_modules_imported"] = _forbidden_module_names()

    source_texts: dict[str, str] = {}
    for role, path in zip(SOURCE_TEXT_ROLES, SOURCE_TEXT_INPUTS):
        try:
            source_texts[role] = raw[path].decode("utf-8")
        except UnicodeDecodeError as exc:
            raise RunnerInputError("source input is not UTF-8: " + path) from exc
    try:
        parent_report_value = json.loads(raw[PARENT_CUTOFF_REL].decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError("parent cutoff report is not valid UTF-8 JSON") from exc
    _require(isinstance(parent_report_value, Mapping), "parent cutoff report is not an object")
    parent_report = dict(parent_report_value)

    auditor = getattr(leaf, "audit_source_texts", None)
    _require(callable(auditor), "v2 leaf must export audit_source_texts")
    run_fixture = getattr(fixture, "run_fixture", None)
    _require(callable(run_fixture), "v2 fixture must export run_fixture(auditor)")
    progress["verification_invocations"] = 1
    progress["consumer_calls"] = ["fixture.run_fixture(audit_source_texts)"]
    fixture_output = run_fixture(auditor)
    progress["fixture_raw"] = _json_safe(fixture_output)
    fixture_result = _validate_fixture(fixture_output)
    progress["fixture"] = fixture_result
    _require(not _forbidden_module_names(), "v2 fixture imported a historical or v1 module")
    progress["fixture_inputs"] = {
        "fixture_owned_synthetic_texts": True,
        "fixture_source_texts_supplied_by_runner": False,
        "fixture_file_path": FIXTURE_REL,
        "fixture_file_sha256": hashes[FIXTURE_REL],
        "negative_case_count": 3,
        "actual_source_texts_consumed": False,
    }

    progress["source_texts"] = {
        "roles": list(SOURCE_TEXT_ROLES),
        "paths": list(SOURCE_TEXT_INPUTS),
        "sha256": {path: hashes[path] for path in SOURCE_TEXT_INPUTS},
        "bytes": {path: len(raw[path]) for path in SOURCE_TEXT_INPUTS},
        "passed_verbatim": True,
        "decoded_after_entry": True,
    }
    progress["parent_report"] = {
        "path": PARENT_CUTOFF_REL,
        "sha256": hashes[PARENT_CUTOFF_REL],
        "bytes": len(raw[PARENT_CUTOFF_REL]),
        "decoded_after_entry": True,
    }
    progress["v1_failure_archive"] = {
        "path": V1_FAILURE_REL,
        "sha256": hashes[V1_FAILURE_REL],
        "bytes": len(raw[V1_FAILURE_REL]),
        "raw_only": True,
        "decoded": False,
        "decompressed": False,
    }

    source_labels = {
        "volterra": INNER_VOLTERRA_REL,
        "regular": REGULAR_SERIES_REL,
        "first_raw": ORIGINAL_RAW_REL,
    }
    progress["actual_leaf_calls"] = 1
    progress["consumer_calls"].append(
        "leaf.audit_source_texts(volterra_text, regular_text, first_raw_text, parent=..., source_labels=...)"
    )
    result = auditor(
        source_texts["volterra"],
        source_texts["regular"],
        source_texts["first_raw"],
        parent=parent_report,
        source_labels=source_labels,
    )
    progress["verification_invocations"] = 2
    encoded_result = _json_safe(result)
    progress["actual_result"] = encoded_result
    _require(isinstance(result, Mapping), "v2 leaf result is not an object")
    _require(isinstance(encoded_result, Mapping), "v2 encoded leaf result is not an object")
    _require(not _forbidden_module_names(), "v2 leaf imported a historical or v1 module")
    _require(encoded_result.get("passed") is True, "v2 actual source proof did not pass")
    checks = encoded_result.get("checks")
    _require(isinstance(checks, Mapping) and bool(checks), "v2 actual source checks are missing")
    _require(all(item is True for item in checks.values()), "v2 actual source proof has a failed check")
    for flag in (
        "complete_kchi_proved",
        "full_kchi_proved",
        "complete_Kchi_resolvent_provided",
        "full_Kchi_proved",
        "ledger_modified",
        "full_goal_complete",
    ):
        _require(encoded_result.get(flag) is not True, "v2 leaf attempted an out-of-scope claim: " + flag)

    return {
        "status": "returned",
        "passed": True,
        "fixture": fixture_result,
        "leaf_result": encoded_result,
        "source_texts": progress["source_texts"],
        "fixture_inputs": progress["fixture_inputs"],
        "parent_report": progress["parent_report"],
        "v1_failure_archive": progress["v1_failure_archive"],
        "actual_leaf_calls": 1,
        "ancestor_math_calls": 0,
        "historical_producer_imported": False,
        "source_only": True,
        "actual_source_metadata_only": False,
        "actual_source_scope": _actual_scope(hashes),
    }


def main() -> int:
    source, hashes, source_sha, brief_sha, preimage_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTRY,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "source_manifest_sha256": source_sha,
            "brief_preimage_sha256": preimage_sha,
            "input_hashes": hashes,
            "raw_sha_authenticated_before_entry": True,
            "v1_failure_archive_decoded": False,
            "v1_failure_archive_decompressed": False,
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "historical_producer_imported": False,
            "runtime_model_metadata": None,
        },
    )

    progress: dict[str, Any] = {
        "fixture": None,
        "actual_result": None,
        "actual_leaf_calls": 0,
        "consumer_calls": [],
        "verification_invocations": 0,
    }
    try:
        nested = _consume(source, hashes, source_sha, brief_sha, preimage_sha, entry_sha, progress)
        calls = list(progress.get("consumer_calls", ()))
        output = _base_result(source, hashes, source_sha, brief_sha, preimage_sha, entry_sha,
                              calls, int(progress.get("verification_invocations", 0)))
        output.update(nested)
        output["consumer_calls"] = calls
        output["verification_invocations"] = int(progress.get("verification_invocations", 0))
        output["error"] = None
        output["traceback"] = None
    except BaseException as exc:
        calls = list(progress.get("consumer_calls", ()))
        output = _base_result(source, hashes, source_sha, brief_sha, preimage_sha, entry_sha,
                              calls, int(progress.get("verification_invocations", 0)))
        output.update(
            {
                "status": "failed",
                "passed": False,
                "fixture": progress.get("fixture") if progress.get("fixture") is not None
                else progress.get("fixture_raw"),
                "actual_result": progress.get("actual_result"),
                "progress": _json_safe(progress),
                "actual_leaf_calls": int(progress.get("actual_leaf_calls", 0)),
                "error": {"type": type(exc).__name__, "message": str(exc)},
                "traceback": traceback.format_exc(),
            }
        )

    result_sha = _write_exclusive(RESULT, output)
    print(json.dumps({
        "component": COMPONENT,
        "phase": PHASE,
        "status": output["status"],
        "result_sha256": result_sha,
        "consumer_calls": output["consumer_calls"],
        "verification_invocations": output["verification_invocations"],
    }, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
