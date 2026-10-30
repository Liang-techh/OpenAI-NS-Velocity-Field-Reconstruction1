"""One-shot, non-reentrant runner for the pulse log-radius transport fixture.

The runner is intentionally separate from the mathematical module.  It first
authenticates the pending brief and its source checkpoint, writes an exclusive
entry receipt, then imports and calls the new fixture exactly once.  A result
artifact is written for both success and failure, including the complete error
traceback on failure.  It does not import or execute any historical fixture.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from pathlib import Path
from typing import Any, Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "tmp" / "pulse_log_radius_transport_v1"
SOURCE = LEAF / "source.json"
CAPSULE = LEAF / "capsule.json"
ENTRY = LEAF / "entered.json"
RESULT = LEAF / "result.json"
CODE = ROOT / "src" / "openai_ns_reconstruction" / "pulse_log_radius_transport.py"
COMPONENT = "pulse_log_radius_transport_v1"
PHASE = "pulse_log_radius_transport_pending"
HEX64 = re.compile(r"^[0-9a-fA-F]{64}$")


def sha256(path: Path) -> str:
    """Hash a file in bounded chunks."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def save_new(path: Path, value: Any) -> str:
    """Write exactly once and return the hash of the persisted bytes."""

    raw = _json_bytes(value)
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _resolve(raw_path: Any) -> Path:
    if not isinstance(raw_path, str):
        raise RuntimeError(f"hash record path is not text: {raw_path!r}")
    path = Path(raw_path)
    if not path.is_absolute():
        path = ROOT / path
    return path


def _check_file_hash(path: Path, expected: Any, label: str) -> str:
    if not isinstance(expected, str) or HEX64.fullmatch(expected) is None:
        raise RuntimeError(f"{label} has no valid SHA-256")
    actual = sha256(path)
    if actual.lower() != expected.lower():
        raise RuntimeError(f"{label} hash mismatch: {path}")
    return actual


def _record_path_hashes(value: Any) -> Iterable[tuple[str, str, str]]:
    """Yield path/hash records from the supported checkpoint shapes.

    Existing checkpoint writers use either ``files: [{path, sha256}]`` or a
    ``sha256: {path: digest}`` map.  A small recursive reader also accepts a
    nested ``artifacts``/``hashes`` record without trusting unrelated scalar
    metadata.
    """

    if isinstance(value, Mapping):
        path_value = value.get("path")
        digest_value = value.get("sha256")
        if isinstance(path_value, str) and isinstance(digest_value, str):
            if HEX64.fullmatch(digest_value):
                yield (path_value, digest_value, "path/sha256")
        for key in ("sha256", "hashes", "files", "artifacts", "scoped_file_hashes"):
            nested = value.get(key)
            if isinstance(nested, Mapping):
                for name, digest in nested.items():
                    if isinstance(name, str) and isinstance(digest, str) and HEX64.fullmatch(digest):
                        yield (name, digest, f"{key}[{name}]")
            elif isinstance(nested, list):
                for item in nested:
                    yield from _record_path_hashes(item)
        # Check other nested mappings only when they explicitly look like a
        # checkpoint; this avoids traversing historical brief inventories.
        for key, nested in value.items():
            if key in {"sha256", "hashes", "files", "artifacts", "scoped_file_hashes"}:
                continue
            if isinstance(nested, Mapping) and ("path" in nested or "sha256" in nested):
                yield from _record_path_hashes(nested)
    elif isinstance(value, list):
        for item in value:
            yield from _record_path_hashes(item)


def _deduplicate_records(records: Iterable[tuple[str, str, str]]) -> list[tuple[str, str, str]]:
    seen: set[tuple[str, str]] = set()
    result: list[tuple[str, str, str]] = []
    for path, digest, label in records:
        key = (path, digest.lower())
        if key not in seen:
            seen.add(key)
            result.append((path, digest, label))
    return result


def _brief_checkpoint_record(brief: Mapping[str, Any]) -> tuple[Path, str]:
    """Find the target checkpoint pointer and its expected hash in the brief."""

    names = (
        f"{COMPONENT}_source_checkpoint",
        "pulse_log_radius_transport_source_checkpoint",
        "source_checkpoint",
    )
    for name in names:
        value = brief.get(name)
        if isinstance(value, Mapping):
            path_value = value.get("path")
            digest_value = value.get("sha256")
            if isinstance(path_value, str) and isinstance(digest_value, str):
                return _resolve(path_value), digest_value
        if isinstance(value, str):
            digest_value = brief.get(name + "_sha256")
            if isinstance(digest_value, str):
                return _resolve(value), digest_value
    # Some pending briefs place the component record under ``authority`` or
    # ``pending``.  Only records naming this component are considered.
    stack: list[Any] = [brief]
    while stack:
        current = stack.pop()
        if isinstance(current, Mapping):
            component_text = " ".join(str(current.get(k, "")) for k in ("component", "name", "key"))
            if COMPONENT in component_text or "pulse_log_radius_transport" in component_text:
                for key in ("source_checkpoint", "checkpoint"):
                    value = current.get(key)
                    if isinstance(value, Mapping):
                        path_value = value.get("path")
                        digest_value = value.get("sha256")
                        if isinstance(path_value, str) and isinstance(digest_value, str):
                            return _resolve(path_value), digest_value
            stack.extend(current.values())
        elif isinstance(current, list):
            stack.extend(current)
    raise RuntimeError("pending brief has no pulse log-radius source checkpoint hash")


def _preflight() -> tuple[dict[str, Any], dict[str, Any], str, str]:
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError("already entered; never reenter")
    if not BRIEF.is_file() or not SOURCE.is_file() or not CAPSULE.is_file():
        raise RuntimeError("brief, source checkpoint, and capsule must already exist")
    brief = json.loads(BRIEF.read_bytes())
    source = json.loads(SOURCE.read_bytes())
    if not isinstance(brief, Mapping) or not isinstance(source, Mapping):
        raise RuntimeError("brief/source checkpoint must be JSON objects")
    if brief.get("phase") != PHASE or brief.get("component") != COMPONENT:
        raise RuntimeError("brief is not pending for this component")
    for flag in ("goal_active", "full_goal_complete", "do_not_reenter"):
        if flag in brief:
            expected = {"goal_active": True, "full_goal_complete": False, "do_not_reenter": True}[flag]
            if brief.get(flag) is not expected:
                raise RuntimeError(f"brief guard {flag!r} failed")

    checkpoint_path, checkpoint_expected = _brief_checkpoint_record(brief)
    checkpoint_actual = _check_file_hash(checkpoint_path, checkpoint_expected, "brief checkpoint")
    if checkpoint_path.resolve() != SOURCE.resolve():
        raise RuntimeError("brief checkpoint does not point to this component source.json")

    records = _deduplicate_records(_record_path_hashes(source))
    if not records:
        raise RuntimeError("source checkpoint has no file hashes")
    for raw_path, expected, label in records:
        _check_file_hash(_resolve(raw_path), expected, label)
    # A source checkpoint is expected to authenticate the capsule as well, but
    # tolerate manifests that use only the source files map while still making
    # the capsule hash explicit in the entry receipt.
    capsule_actual = sha256(CAPSULE)
    return dict(brief), dict(source), checkpoint_actual, capsule_actual


def main() -> int:
    brief, source, checkpoint_sha, capsule_sha = _preflight()
    source_sha = sha256(SOURCE)
    receipt = save_new(
        ENTRY,
        {
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "brief_sha256": sha256(BRIEF),
            "source_sha256": source_sha,
            "checkpoint_sha256": checkpoint_sha,
            "capsule_sha256": capsule_sha,
            "consumer_calls_started": 0,
            "fixture_calls": 0,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
        },
    )
    try:
        spec = importlib.util.spec_from_file_location(
            "_pulse_log_radius_transport_v1", CODE
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("new fixture load failed")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        # Exactly one mathematical fixture invocation, after the exclusive
        # entry receipt has been persisted.
        value = module.verify_synthetic_fixture()
        output: dict[str, Any] = {
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "brief_sha256": sha256(BRIEF),
            "source_sha256": sha256(SOURCE),
            "entry_sha256": receipt,
            "checkpoint_sha256": sha256(SOURCE),
            "capsule_sha256": sha256(CAPSULE),
            "fixture_calls": 1,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
            "error": None,
            "result": value,
        }
    except Exception as exc:  # noqa: BLE001 - complete error artifact is required
        output = {
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "brief_sha256": sha256(BRIEF),
            "source_sha256": sha256(SOURCE),
            "entry_sha256": receipt,
            "checkpoint_sha256": sha256(SOURCE),
            "capsule_sha256": sha256(CAPSULE),
            "fixture_calls": 1,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }
    result_sha = save_new(RESULT, output)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PHASE,
                "status": output["status"],
                "result_sha256": result_sha,
                "fixture_calls": 1,
                "ancestor_math_calls": 0,
                "runtime_model_metadata": None,
            },
            sort_keys=True,
        )
    )
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
