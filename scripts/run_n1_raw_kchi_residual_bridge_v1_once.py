"""One-shot source binding receipt for the N1 raw/Kchi residual bridge.

Only the new fully-qualified leaf and fixture are imported.  The fixture is
called once and is the only mathematical consumer in this runner: it checks
the exact-Fraction synthetic operator identity.  The actual-source section is
metadata-only.  It binds the frozen document/source files by raw SHA and
records the four prerequisites that remain open; it does not evaluate a real
residual or chi function and it does not call a historical producer.
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
COMPONENT = "n1_raw_kchi_residual_bridge_v1"
PHASE = COMPONENT + "_pending"

BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_MODULE = "openai_ns_reconstruction.n1_raw_kchi_residual_bridge_v1"
FIXTURE_MODULE = "openai_ns_reconstruction.n1_raw_kchi_residual_bridge_fixture_v1"

LEAF_REL = "src/openai_ns_reconstruction/n1_raw_kchi_residual_bridge_v1.py"
FIXTURE_REL = (
    "src/openai_ns_reconstruction/"
    "n1_raw_kchi_residual_bridge_fixture_v1.py"
)
RUNNER_REL = "scripts/run_n1_raw_kchi_residual_bridge_v1_once.py"
DOC_FEEDBACK_REL = "docs/CURRENT_N1_CORE_K0_FEEDBACK_M15_2026_10_30.md"
DOC_COMMON_REL = "docs/CURRENT_N1_COMMON_SOURCE_BRIDGE_2026_11_02.md"
ORIGINAL_RAW_REL = "src/openai_ns_reconstruction/original_n1_first_raw_function.py"
KCHI_FIRST_REL = "src/openai_ns_reconstruction/n1_kchi_first_correction_v1.py"

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    RUNNER_REL,
    DOC_FEEDBACK_REL,
    DOC_COMMON_REL,
    ORIGINAL_RAW_REL,
    KCHI_FIRST_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)
SOURCE_REFERENCE_INPUTS = (
    DOC_FEEDBACK_REL,
    DOC_COMMON_REL,
    ORIGINAL_RAW_REL,
    KCHI_FIRST_REL,
)
SHA256_HEX = frozenset("0123456789abcdefABCDEF")

UNMET_PREREQUISITES = (
    {
        "key": "source_owner",
        "satisfied": False,
        "status": "unmet",
        "requirement": "one shared live owner for the real residual and chi functions",
        "reason": "no shared live owner is available in this source-only round",
    },
    {
        "key": "resolvent",
        "satisfied": False,
        "status": "unmet",
        "requirement": "a real-data resolvent/domain and invertibility certificate",
        "reason": "no real resolvent was evaluated or certified",
    },
    {
        "key": "same_space_unit",
        "satisfied": False,
        "status": "unmet",
        "requirement": "residual and chi in the same space and unit convention",
        "reason": "same-space and unit compatibility was not established",
    },
    {
        "key": "derivative",
        "satisfied": False,
        "status": "unmet",
        "requirement": "the real-data derivative/domain premise for the comparison",
        "reason": "no real derivative hypothesis or identity was checked",
    },
)


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


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
    """Make fixture output finite and JSON serializable without re-evaluation."""

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
            items = [_json_safe(item, active) for item in value]
            return sorted(items, key=repr)
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
    """Persist one receipt with OS-level exclusive creation."""

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


def _load_object(path: Path, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_bytes())
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not valid JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


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
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(),
             label + " path is missing or escapes repository: " + normalized)
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
    entries: list[tuple[Any, Any, str]] = []
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            entries.append((item.get("path"), item.get("sha256"), f"source.files[{index}]"))
    else:
        for field in ("raw_file_sha256", "input_sha256", "input_hashes", "sha256"):
            table = source.get(field)
            if isinstance(table, Mapping):
                entries.extend(
                    (raw_path, expected, f"source.{field}")
                    for raw_path, expected in table.items()
                )
                break
    _require(entries, "source manifest contains no raw-file SHA entries")

    hashes: dict[str, str] = {}
    for raw_path, expected, label in entries:
        key, path = _repo_file(raw_path, label)
        _require(key in FROZEN_SET, "manifest contains out-of-scope file: " + key)
        digest = _sha_text(expected, key)
        _require(key not in hashes, "duplicate manifest path: " + key)
        actual = _sha_file(path)
        _require(actual.lower() == digest, "raw SHA mismatch: " + key)
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET,
             "manifest must freeze exactly the seven required files")
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _checkpoint_sha(brief: Mapping[str, Any]) -> str:
    checkpoint: Any = brief.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        checkpoint = brief.get(COMPONENT + "_source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        continuation = brief.get(COMPONENT)
        if isinstance(continuation, Mapping):
            checkpoint = continuation.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        continuation = brief.get("continuation")
        if isinstance(continuation, Mapping):
            checkpoint = continuation.get("source_checkpoint")
    _require(isinstance(checkpoint, Mapping), "pending brief source checkpoint missing")
    raw_path = checkpoint.get("path")
    _require(isinstance(raw_path, str), "pending brief source checkpoint path missing")
    path = Path(raw_path)
    if not path.is_absolute():
        path = ROOT / path
    _require(path.resolve() == SOURCE.resolve(), "pending brief points to another source")
    expected = checkpoint.get("sha256", checkpoint.get("manifest_sha256"))
    actual = _sha_file(SOURCE)
    _require(_sha_text(expected, "pending brief source checkpoint") == actual,
             "pending brief source checkpoint SHA mismatch")
    return actual


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


def _preflight() -> tuple[dict[str, Any], dict[str, str], str, str]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    _require(BRIEF.is_file(), "short continuation brief is missing")
    _require(SOURCE.is_file(), "source manifest is missing")
    brief = _load_object(BRIEF, "short continuation brief")
    record = _brief_record(brief)
    _require(record.get("component") == COMPONENT, "pending brief component mismatch")
    _require(record.get("phase") == PHASE, "pending brief phase mismatch")
    source = _load_object(SOURCE, "source manifest")
    _require(source.get("component") == COMPONENT, "source manifest component mismatch")
    if "phase" in source:
        _require(source.get("phase") == PHASE, "source manifest phase mismatch")
    source_sha = _checkpoint_sha(brief)
    hashes = _manifest_hashes(source)
    return source, hashes, source_sha, _sha_file(BRIEF)


def _load_new_modules() -> tuple[Any, Any, list[str]]:
    before = set(sys.modules)
    source_root = str(SRC)
    package_name = "openai_ns_reconstruction"
    package_dir = SRC / package_name
    previous_package = sys.modules.get(package_name)
    namespace = ModuleType(package_name)
    namespace.__path__ = [str(package_dir)]
    namespace.__package__ = package_name
    namespace.__spec__ = importlib.machinery.ModuleSpec(
        package_name, loader=None, is_package=True
    )
    sys.modules[package_name] = namespace
    old_path = list(sys.path)
    if source_root not in sys.path:
        sys.path.insert(0, source_root)
    try:
        # Use a namespace package so importing these two fully-qualified new
        # modules does not execute openai_ns_reconstruction/__init__.py or any
        # unrelated ancestor module.  The fixture's relative import then
        # reuses the exact leaf module loaded immediately above.
        leaf_spec = importlib.util.spec_from_file_location(LEAF_MODULE, ROOT / LEAF_REL)
        fixture_spec = importlib.util.spec_from_file_location(FIXTURE_MODULE, ROOT / FIXTURE_REL)
        _require(leaf_spec is not None and leaf_spec.loader is not None,
                 "new leaf module spec is unavailable")
        _require(fixture_spec is not None and fixture_spec.loader is not None,
                 "new fixture module spec is unavailable")
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
    forbidden = []
    for name in sorted(set(sys.modules) - before):
        base = name.rsplit(".", 1)[-1]
        if base in {"original_n1_first_raw_function", "n1_kchi_first_correction_v1"}:
            forbidden.append(name)
    _require(not forbidden, "new leaf/fixture imported a historical producer")
    return leaf, fixture, forbidden


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not a JSON object")
    result = dict(value)
    _require(result.get("passed") is True, "synthetic fixture did not pass")
    checks = result.get("checks")
    _require(isinstance(checks, Mapping) and bool(checks), "fixture checks are missing")
    _require(all(item is True for item in checks.values()), "fixture has a failed check")
    _require(result.get("arithmetic") == "exact Fraction only",
             "fixture did not declare exact Fraction arithmetic")
    _require(result.get("historical_modules_called") is False,
             "fixture called a historical module")
    _require(result.get("mpmath_used") is False, "fixture used mpmath")
    _json_bytes(result)
    return result


def _actual_source_scope(hashes: Mapping[str, str]) -> dict[str, Any]:
    return {
        "source_only": True,
        "binding_kind": "saved_source_document_and_source_semantics_only",
        "source_references": [
            {"path": path, "sha256": hashes[path], "raw_file_reference": True}
            for path in SOURCE_REFERENCE_INPUTS
        ],
        "source_values_computed": False,
        "numeric_values_computed": False,
        "real_residual_or_chi_evaluated": False,
        "function_identity_checked": False,
        "actual_source_identity_proved": False,
        "unsatisfied_prerequisites": list(UNMET_PREREQUISITES),
    }


def _base_result(
    source: Mapping[str, Any],
    hashes: Mapping[str, str],
    source_sha: str,
    brief_sha: str,
    entry_sha: str,
    calls: list[str],
    invocation_count: int,
    actual_scope: Mapping[str, Any],
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
        "entry_sha256": entry_sha,
        "input_hashes": dict(hashes),
        "raw_file_sha256": dict(hashes),
        "source_manifest": dict(source),
        "consumer_calls": list(calls),
        "verification_invocations": invocation_count,
        "ancestor_math_calls": 0,
        "historical_producer_imported": False,
        "metadata_only": False,
        "actual_source_metadata_only": True,
        "source_only": True,
        "function_identity_checked": False,
        "synthetic_exact_fraction_operator_identity": False,
        "real_residual_or_chi_evaluated": False,
        "actual_leaf_called": False,
        "actual_source_scope": dict(actual_scope),
        "full_goal_complete": False,
        "runtime_model_metadata": None,
    }


def main() -> int:
    source, hashes, source_sha, brief_sha = _preflight()
    entry_sha = _write_exclusive(
        ENTRY,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "source_manifest_sha256": source_sha,
            "input_hashes": hashes,
            "raw_file_sha256": hashes,
            "consumer_calls_started": 0,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
        },
    )

    calls: list[str] = []
    invocation_count = 0
    fixture_result: dict[str, Any] | None = None
    fixture_output: Any = None
    actual_scope = _actual_source_scope(hashes)
    try:
        _leaf, fixture, _forbidden = _load_new_modules()
        run_fixture = getattr(fixture, "run_fixture", None)
        _require(callable(run_fixture), "fixture has no run_fixture()")
        calls.append("fixture.run_fixture()")
        invocation_count += 1
        fixture_output = run_fixture()
        fixture_result = _validate_fixture(fixture_output)
        output = _base_result(
            source, hashes, source_sha, brief_sha, entry_sha,
            calls, invocation_count, actual_scope,
        )
        output.update(
            {
                "status": "returned",
                "passed": True,
                "synthetic_fixture_passed": True,
                "synthetic_exact_fraction_operator_identity": True,
                "fixture": fixture_result,
                "error": None,
                "traceback": None,
            }
        )
    except Exception as exc:
        output = _base_result(
            source, hashes, source_sha, brief_sha, entry_sha,
            calls, invocation_count, actual_scope,
        )
        output.update(
            {
                "status": "failed",
                "passed": False,
                "synthetic_fixture_passed": fixture_result is not None,
                "fixture": fixture_result if fixture_result is not None else _json_safe(fixture_output),
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
        "consumer_calls": calls,
        "verification_invocations": invocation_count,
    }, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
