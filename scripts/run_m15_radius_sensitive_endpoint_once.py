"""Single-entry runner for the caller-owned radius-sensitive m15 leaf.

The pending brief and frozen source checkpoint are supplied by the root task.
This script never creates or edits either one.  After validating every frozen
input hash it writes ``entered.json`` exclusively, imports only the new leaf,
invokes its synthetic fixture once, and saves the complete returned packet or
the complete failure traceback in ``result.json``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "tmp/m15_radius_sensitive_endpoint_v1"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "src/openai_ns_reconstruction/m15_radius_sensitive_endpoint.py"
COMPONENT = "m15_radius_sensitive_endpoint_v1"
PENDING_PHASE = "m15_radius_sensitive_endpoint_pending"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_default(value: Any) -> Any:
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "_mpf_"):
        return str(value)
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def _write_exclusive(path: Path, value: Mapping[str, Any]) -> str:
    raw = (
        json.dumps(value, sort_keys=True, separators=(",", ":"),
                   ensure_ascii=False, default=_json_default) + "\n"
    ).encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _hash_mapping(source: Mapping[str, Any]) -> dict[str, str]:
    hashes = source.get("input_hashes")
    if not isinstance(hashes, Mapping) or not hashes:
        hashes = source.get("sha256")
    if not isinstance(hashes, Mapping) or not hashes:
        raise RuntimeError("source checkpoint has no input hash mapping")
    result: dict[str, str] = {}
    for name, expected in hashes.items():
        relative = str(name)
        path = Path(relative)
        if path.is_absolute() or ".." in path.parts:
            raise RuntimeError("checkpoint input path is not repository-relative: " + relative)
        digest = str(expected).lower()
        if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise RuntimeError("invalid SHA-256 in source checkpoint: " + relative)
        result[relative] = digest
    return result


def _preflight() -> tuple[dict[str, Any], dict[str, Any], str]:
    if not BRIEF.is_file():
        raise RuntimeError("pending brief is missing: " + str(BRIEF))
    brief = json.loads(BRIEF.read_bytes())
    if brief.get("phase") != PENDING_PHASE or brief.get("component") != COMPONENT:
        raise RuntimeError("wrong pending brief phase or component")
    source_path = DIRECTORY / "source.json"
    if not source_path.is_file():
        raise RuntimeError("source checkpoint is missing: " + str(source_path))
    if brief.get("source_checkpoint", {}).get("sha256") != _sha(source_path):
        raise RuntimeError("pending brief checkpoint hash mismatch")
    source = json.loads(source_path.read_bytes())
    if source.get("component") != COMPONENT:
        raise RuntimeError("source checkpoint component mismatch")
    for relative, expected in _hash_mapping(source).items():
        path = ROOT / relative
        if not path.is_file() or _sha(path) != expected:
            raise RuntimeError("frozen input changed: " + relative)
    return brief, source, _sha(source_path)


def main() -> int:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    entered_path, result_path = DIRECTORY / "entered.json", DIRECTORY / "result.json"
    if entered_path.exists() or result_path.exists():
        raise RuntimeError("m15 radius-sensitive endpoint already entered; never re-enter")

    _brief, source, source_sha256 = _preflight()
    entry_sha256 = _write_exclusive(
        entered_path,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "entered": True,
            "source_sha256": source_sha256,
            "input_hashes": dict(_hash_mapping(source)),
            "verification_function": "verify_synthetic_fixture",
            "verification_invocations": 0,
            "source_values_computed": False,
            "source_equality_proved": False,
            "runtime_model_metadata": None,
        },
    )

    invocation_count = 0
    try:
        spec = importlib.util.spec_from_file_location(
            "m15_radius_sensitive_endpoint_new_leaf", LEAF
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("new leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        invocation_count = 1
        returned = leaf.verify_synthetic_fixture()
        if returned.get("passed") is not True:
            raise RuntimeError("synthetic fixture did not pass")
        result: dict[str, Any] = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "status": "returned",
            "source_sha256": source_sha256,
            "entry_sha256": entry_sha256,
            "verification_function": "verify_synthetic_fixture",
            "verification_invocations": invocation_count,
            "consumer_calls": invocation_count,
            "source_values_computed": False,
            "source_equality_proved": False,
            "runtime_model_metadata": None,
            "result": returned,
        }
    except Exception:
        result = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "status": "failed",
            "source_sha256": source_sha256,
            "entry_sha256": entry_sha256,
            "verification_function": "verify_synthetic_fixture",
            "verification_invocations": invocation_count,
            "consumer_calls": invocation_count,
            "source_values_computed": False,
            "source_equality_proved": False,
            "runtime_model_metadata": None,
            "traceback": traceback.format_exc(),
        }

    result_sha256 = _write_exclusive(result_path, result)
    print(json.dumps({"component": COMPONENT, "phase": PENDING_PHASE,
                      "status": result["status"], "result_sha256": result_sha256,
                      "verification_invocations": invocation_count}, sort_keys=True))
    return 0 if result["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())


