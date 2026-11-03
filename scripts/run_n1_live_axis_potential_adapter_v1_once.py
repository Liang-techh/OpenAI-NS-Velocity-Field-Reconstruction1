"""One-shot runner for the fresh live-axis potential adapter fixture.

This runner authenticates only the manifest bytes, writes its entry marker before
loading the new files, and calls exactly ``fixture.run_fixture(adapter_module)``.
The weighted propagator and raw Volterra source files are hash-only anchors; they
are never imported, decoded, or otherwise inspected by this runner.
"""

from __future__ import annotations

import builtins
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
PACKAGE = "openai_ns_reconstruction"
PACKAGE_DIR = SRC / PACKAGE
COMPONENT = "n1_live_axis_potential_adapter_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_MODULE = PACKAGE + ".n1_live_axis_potential_adapter_v1"
FIXTURE_MODULE = PACKAGE + ".n1_live_axis_potential_adapter_fixture_v1"
LEAF_REL = "src/openai_ns_reconstruction/n1_live_axis_potential_adapter_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_live_axis_potential_adapter_fixture_v1.py"
RUNNER_REL = "scripts/run_n1_live_axis_potential_adapter_v1_once.py"
WEIGHT_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_weighted_propagator.py"
)
VOLTERRA_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
FROZEN_INPUTS = (LEAF_REL, FIXTURE_REL, RUNNER_REL, WEIGHT_REL, VOLTERRA_REL)
FROZEN_SET = frozenset(FROZEN_INPUTS)
FIXTURE_SCHEMA = "n1_live_axis_potential_adapter_fixture_v1"
SHA256_HEX = frozenset("0123456789abcdefABCDEF")
ALLOWED_MODULES = frozenset({PACKAGE, LEAF_MODULE, FIXTURE_MODULE})


class RunnerInputError(ValueError):
    """Raised when the one-shot contract is not satisfied."""


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


def _json_safe(value: Any, seen: set[int] | None = None) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"type": "bytes", "length": len(value), "sha256": _sha_bytes(value)}
    active = seen if seen is not None else set()
    identity = id(value)
    if identity in active:
        return {"type": "cycle", "repr": repr(value)}
    active.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(k): _json_safe(v, active) for k, v in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(v, active) for v in value]
        if isinstance(value, (set, frozenset)):
            return sorted((_json_safe(v, active) for v in value), key=repr)
        return repr(value)
    finally:
        active.discard(identity)


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(_json_safe(value), ensure_ascii=False, sort_keys=True,
                       separators=(",", ":"), allow_nan=False).encode("utf-8") + b"\n")


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
        return _decode_object(path.read_bytes(), label)
    except OSError as exc:
        raise RunnerInputError(label + " is unreadable") from exc


def _sha_text(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64 and
             all(c in SHA256_HEX for c in value), "invalid SHA-256 for " + label)
    return value.lower()


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
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET,
             "manifest must contain exactly the five frozen inputs")
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _record_for(brief: Mapping[str, Any]) -> Mapping[str, Any]:
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
    candidates: list[Any] = [brief.get("source_checkpoint"),
                             brief.get(COMPONENT + "_source_checkpoint")]
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        candidates.append(nested.get("source_checkpoint"))
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping):
        candidates.append(continuation.get("source_checkpoint"))
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            return candidate
    raise RunnerInputError("pending source checkpoint is missing")


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(),
             "component already entered or has a result")
    _require(BRIEF.is_file(), "resume brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    brief = _decode_object(brief_raw, "resume brief")
    record = _record_for(brief)
    _require(record.get("component") == COMPONENT, "pending brief component mismatch")
    _require(record.get("phase") == PHASE, "pending brief phase mismatch")
    _require(record.get("one_entry_only") is True, "pending one-entry guard missing")
    _require(record.get("do_not_reenter") is True, "pending reentry guard missing")
    checkpoint = _brief_checkpoint(brief)
    checkpoint_path = _repo_file(checkpoint.get("path"), "source checkpoint")[1]
    _require(checkpoint_path.resolve() == SOURCE.resolve(),
             "pending brief points to a different source checkpoint")
    source_sha = _sha_file(SOURCE)
    _require(_sha_text(checkpoint.get("sha256"), "source checkpoint") == source_sha,
             "source checkpoint SHA mismatch")
    source = _load_object(SOURCE, "source manifest")
    _require(source.get("component") == COMPONENT, "source component mismatch")
    _require(source.get("one_entry_only") is True, "source one-entry guard missing")
    if "phase" in source:
        _require(source.get("phase") == PHASE, "source phase mismatch")
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, _sha_bytes(brief_raw)


def _forbidden_modules() -> list[str]:
    found: list[str] = []
    package_root = PACKAGE_DIR.resolve()
    experiments_root = (ROOT / "experiments").resolve()
    for name, module in list(sys.modules.items()):
        if name in ALLOWED_MODULES:
            continue
        if name == PACKAGE or name.startswith(PACKAGE + "."):
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


def _guarded_import(original: Any):
    def checked(name: str, globals_: Any = None, locals_: Any = None,
                fromlist: tuple[str, ...] | None = None, level: int = 0):
        resolved = name
        if level:
            package = globals_.get("__package__", "") if isinstance(globals_, Mapping) else ""
            try:
                resolved = importlib.util.resolve_name(name, package)
            except (ImportError, AttributeError):
                resolved = name
        if resolved == PACKAGE:
            for item in fromlist or ():
                if item != "*" and PACKAGE + "." + item not in ALLOWED_MODULES:
                    raise RunnerInputError("historical package import blocked: " + item)
        elif resolved.startswith(PACKAGE + ".") and resolved not in ALLOWED_MODULES:
            raise RunnerInputError("historical package import blocked: " + resolved)
        return original(name, globals_, locals_, fromlist, level)
    return checked


def _cleanup_modules() -> None:
    for name in list(sys.modules):
        if name == PACKAGE or name.startswith(PACKAGE + "."):
            sys.modules.pop(name, None)


def _load_new_modules() -> tuple[ModuleType, ModuleType]:
    _require(not _forbidden_modules(), "historical package or source module already loaded")
    namespace = ModuleType(PACKAGE)
    namespace.__path__ = [str(PACKAGE_DIR)]
    namespace.__package__ = PACKAGE
    namespace.__spec__ = importlib.machinery.ModuleSpec(PACKAGE, loader=None, is_package=True)
    sys.modules[PACKAGE] = namespace
    old_path = list(sys.path)
    original_import = builtins.__import__
    builtins.__import__ = _guarded_import(original_import)
    try:
        if str(SRC) not in sys.path:
            sys.path.insert(0, str(SRC))
        loaded: list[ModuleType] = []
        for module_name, relative in ((LEAF_MODULE, LEAF_REL), (FIXTURE_MODULE, FIXTURE_REL)):
            path = ROOT / Path(relative)
            spec = importlib.util.spec_from_file_location(module_name, path)
            _require(spec is not None and spec.loader is not None,
                     "new module spec is unavailable: " + relative)
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            spec.loader.exec_module(module)
            loaded.append(module)
        _require(not _forbidden_modules(), "historical module loaded by new files")
        return loaded[0], loaded[1]
    except BaseException:
        _cleanup_modules()
        raise
    finally:
        builtins.__import__ = original_import
        sys.path[:] = old_path


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not a mapping")
    result = dict(value)
    if "schema" in result:
        _require(result.get("schema") == FIXTURE_SCHEMA, "fixture schema mismatch")
    _require(result.get("passed") is True, "synthetic fixture did not pass")
    checks = result.get("checks")
    if checks is not None:
        _require(isinstance(checks, Mapping) and bool(checks),
                 "fixture checks must be a non-empty mapping")
        _require(all(item is True for item in checks.values()),
                 "synthetic fixture contains a failed check")
    for key in ("historical_module_imported", "historical_modules_imported",
                "old_code_called", "old_receipt_decoded", "old_receipt_decompressed",
                "weighted_propagator_imported", "raw_volterra_imported"):
        if key in result:
            _require(result[key] in (False, 0, None), "fixture reported forbidden activity: " + key)
    return _json_safe(result)


def _consume(source: Mapping[str, Any], hashes: Mapping[str, str],
             progress: dict[str, Any]) -> dict[str, Any]:
    after = _manifest_hashes(source)
    _require(dict(after) == dict(hashes), "manifest input changed after entry")
    progress["raw_sha_reauthenticated_after_entry"] = True
    adapter: ModuleType | None = None
    try:
        adapter, fixture = _load_new_modules()
        progress["new_modules_imported"] = True
        run_fixture = getattr(fixture, "run_fixture", None)
        _require(callable(run_fixture), "fixture.run_fixture is missing or not callable")
        progress["fixture_calls"] = 1
        progress["fixture_call"] = "fixture.run_fixture(adapter_module)"
        raw = run_fixture(adapter)
        progress["fixture_raw"] = _json_safe(raw)
        checked = _validate_fixture(raw)
        _require(not _forbidden_modules(), "historical module loaded during fixture")
        return {
            "status": "returned",
            "passed": True,
            "fixture": checked,
            "raw_hashes_after_entry": dict(after),
            "source_anchors_read_only": True,
            "source_anchors_executed": False,
            "runner_adapter_calls": 0,
        }
    finally:
        _cleanup_modules()


def _base_result(source: Mapping[str, Any], hashes: Mapping[str, str], source_sha: str,
                 brief_sha: str, entry_sha: str, progress: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema": COMPONENT + "_runner_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "source_checkpoint_path": SOURCE.relative_to(ROOT).as_posix(),
        "source_checkpoint_sha256": source_sha,
        "source_manifest_sha256": source_sha,
        "brief_path": BRIEF.as_posix(),
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": dict(hashes),
        "raw_file_sha256": dict(hashes),
        "source_manifest": _json_safe(source),
        "one_entry_only": True,
        "do_not_reenter": True,
        "brief_modified": False,
        "fixture_calls": int(progress.get("fixture_calls", 0)),
        "runner_adapter_calls": 0,
        "progress": _json_safe(progress),
    }


def main() -> int:
    try:
        source, hashes, source_sha, brief_sha = _preflight()
    except BaseException as exc:
        print(json.dumps({"component": COMPONENT, "status": "preflight_failed",
                          "error": str(exc), "traceback": traceback.format_exc()},
                         ensure_ascii=False, sort_keys=True))
        return 2

    entry = {
        "schema": COMPONENT + "_entered_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "entered": True,
        "one_entry_only": True,
        "do_not_reenter": True,
        "source_checkpoint_sha256": source_sha,
        "source_manifest_sha256": source_sha,
        "brief_sha256": brief_sha,
        "input_hashes": dict(hashes),
        "raw_sha_authenticated_before_entry": True,
        "source_anchors_read_only": True,
        "source_anchors_executed": False,
        "fixture_calls": 0,
        "runner_adapter_calls": 0,
    }
    try:
        entry_sha = _write_exclusive(ENTRY, entry)
    except BaseException as exc:
        print(json.dumps({"component": COMPONENT, "status": "entry_failed",
                          "error": str(exc), "traceback": traceback.format_exc()},
                         ensure_ascii=False, sort_keys=True))
        return 2

    progress: dict[str, Any] = {"fixture_calls": 0, "runner_adapter_calls": 0}
    try:
        nested = _consume(source, hashes, progress)
        output = _base_result(source, hashes, source_sha, brief_sha, entry_sha, progress)
        output.update(nested)
        output.update({"error": None, "traceback": None})
        status = "returned"
    except BaseException as exc:
        output = _base_result(source, hashes, source_sha, brief_sha, entry_sha, progress)
        output.update({"status": "failed", "passed": False,
                       "error": type(exc).__name__ + ": " + str(exc),
                       "traceback": traceback.format_exc(),
                       "source_anchors_read_only": True,
                       "source_anchors_executed": False})
        status = "failed"
    try:
        result_sha = _write_exclusive(RESULT, output)
    except BaseException as exc:
        print(json.dumps({"component": COMPONENT, "status": "result_write_failed",
                          "entry_sha256": entry_sha, "error": str(exc),
                          "traceback": traceback.format_exc()},
                         ensure_ascii=False, sort_keys=True))
        return 2
    print(json.dumps({"component": COMPONENT, "phase": PHASE, "status": status,
                      "entry_sha256": entry_sha, "result_sha256": result_sha,
                      "source_checkpoint_sha256": source_sha,
                      "brief_sha256": brief_sha, "input_hashes": dict(hashes),
                      "fixture_calls": int(progress.get("fixture_calls", 0)),
                      "runner_adapter_calls": 0}, ensure_ascii=False, sort_keys=True))
    return 0 if status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
