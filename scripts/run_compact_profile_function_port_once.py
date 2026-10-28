"""Single-entry runner for the compact profile numerical port.

The runner is intentionally inert until the root task provides its pending
brief and source checkpoint.  It validates every checkpoint hash and creates
an exclusive ``entered.json`` receipt before importing the leaf.  The leaf's
synthetic verification function is invoked exactly once; source-value
computation remains explicitly false.  A complete result receipt is written
exclusively even when validation or the one invocation fails.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "tmp/compact_profile_function_port_v1"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "src/openai_ns_reconstruction/compact_profile_function_port.py"
COMPONENT = "compact_profile_function_port_v1"
PENDING_PHASE = "compact_profile_function_port_pending"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_default(value: Any) -> Any:
    """Serialize mpmath values without converting huge exponents to float."""

    if hasattr(value, "_mpf_"):
        return str(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def _write_exclusive(path: Path, value: Mapping[str, Any]) -> str:
    raw = (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            default=_json_default,
        )
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _checkpoint_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    hashes = source.get("sha256")
    if not isinstance(hashes, Mapping) or not hashes:
        # A few older staged checkpoints used input_hashes; accepting the
        # equivalent shape keeps validation strict without requiring a new
        # checkpoint format.
        hashes = source.get("input_hashes")
    if not isinstance(hashes, Mapping) or not hashes:
        raise RuntimeError("source checkpoint has no input hash mapping")
    result: dict[str, str] = {}
    for name, expected in hashes.items():
        relative = str(name)
        digest = str(expected)
        if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest.lower()):
            raise RuntimeError("invalid SHA-256 in source checkpoint: " + relative)
        result[relative] = digest.lower()
    return result


def _preflight() -> tuple[dict[str, Any], dict[str, Any], str]:
    """Validate pending state and all source-checkpoint inputs pre-import."""

    if not BRIEF.is_file():
        raise RuntimeError("pending brief is missing: " + str(BRIEF))
    state = json.loads(BRIEF.read_bytes())
    if state.get("phase") != PENDING_PHASE or state.get("component") != COMPONENT:
        raise RuntimeError("wrong pending state for compact profile port")
    source_path = DIRECTORY / "source.json"
    if not source_path.is_file():
        raise RuntimeError("source checkpoint is missing: " + str(source_path))
    source = json.loads(source_path.read_bytes())
    if source.get("component") != COMPONENT:
        raise RuntimeError("source checkpoint component mismatch")
    hashes = _checkpoint_hashes(source)
    for name, expected in hashes.items():
        path = ROOT / name
        if not path.is_file() or _sha(path) != expected:
            raise RuntimeError("frozen input changed: " + name)
    return state, source, _sha(source_path)


def main() -> int:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    entered_path = DIRECTORY / "entered.json"
    result_path = DIRECTORY / "result.json"
    # Both receipts are one-shot artifacts.  Do this before any import and do
    # not overwrite either artifact on a repeated invocation.
    if entered_path.exists() or result_path.exists():
        raise RuntimeError("compact profile port already entered; never re-enter")

    source_sha256: str | None = None
    entry_sha256: str | None = None
    invocation_count = 0
    try:
        state, source, source_sha256 = _preflight()
        entry_sha256 = _write_exclusive(
            entered_path,
            {
                "component": COMPONENT,
                "phase": PENDING_PHASE,
                "entered": True,
                "source_sha256": source_sha256,
                "verification_function": "verify_synthetic_fixture",
                "verification_invocations": 0,
                "source_values_computed": False,
                "source_value_computation": False,
                "numeric_evaluation": True,
                "certified": False,
                "runtime_model_metadata": None,
            },
        )
        spec = importlib.util.spec_from_file_location(
            "compact_profile_function_port_new_leaf", LEAF
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        # Exactly one call into the leaf verification function.  It owns all
        # synthetic inputs, including the interior polynomial and endpoints.
        invocation_count = 1
        returned = leaf.verify_synthetic_fixture()
        if returned.get("source_value_computation") is not False:
            raise RuntimeError("leaf reported source-value computation")
        if returned.get("checks_passed") is not True:
            raise RuntimeError("synthetic verification did not pass")
        result: dict[str, Any] = {
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "status": "returned",
            "source_sha256": source_sha256,
            "entry_sha256": entry_sha256,
            "verification_function": "verify_synthetic_fixture",
            "verification_invocations": invocation_count,
            "consumer_calls": invocation_count,
            "source_values_computed": False,
            "source_value_computation": False,
            "numeric_evaluation": True,
            "certified": False,
            "runtime_model_metadata": None,
            "result": returned,
        }
    except Exception:
        result = {
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "status": "failed",
            "source_sha256": source_sha256,
            "entry_sha256": entry_sha256,
            "verification_function": "verify_synthetic_fixture",
            "verification_invocations": invocation_count,
            "consumer_calls": invocation_count,
            "source_values_computed": False,
            "source_value_computation": False,
            "numeric_evaluation": True,
            "certified": False,
            "runtime_model_metadata": None,
            "traceback": traceback.format_exc(),
        }

    result_sha256 = _write_exclusive(result_path, result)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PENDING_PHASE,
                "status": result["status"],
                "result_sha256": result_sha256,
                "verification_invocations": invocation_count,
                "source_values_computed": False,
                "source_value_computation": False,
                "numeric_evaluation": True,
                "certified": False,
            },
            sort_keys=True,
        )
    )
    return 0 if result["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
