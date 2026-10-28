"""Single-entry runner for the compact scaled-profile numerical fixture.

The runner does not create a pending brief or a source checkpoint.  It first
validates the existing pending state and every input hash declared by the
checkpoint, writes an exclusive ``entered.json`` receipt, and only then loads
the new leaf.  The leaf fixture is invoked exactly once.  A complete
``result.json`` receipt records either the full returned value or the complete
exception traceback.  Re-entry never overwrites either receipt.
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
DIRECTORY = ROOT / "tmp/compact_scaled_profiles_v1"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "src/openai_ns_reconstruction/compact_scaled_profiles.py"
COMPONENT = "compact_scaled_profiles_v1"
PENDING_PHASE = "compact_scaled_profiles_pending"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_default(value: Any) -> Any:
    """Serialize mpmath values without converting them to binary floats."""

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


def _manifest_path(raw_path: Any) -> Path:
    if not isinstance(raw_path, str) or not raw_path:
        raise RuntimeError("checkpoint hash entry has no path")
    path = Path(raw_path)
    return path if path.is_absolute() else ROOT / path


def _add_hash(
    result: dict[str, str], raw_path: Any, raw_digest: Any
) -> None:
    if not isinstance(raw_path, str) or not raw_path:
        raise RuntimeError("checkpoint hash entry has no path")
    digest = str(raw_digest).lower()
    if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
        raise RuntimeError("invalid SHA-256 in source checkpoint: " + raw_path)
    previous = result.get(raw_path)
    if previous is not None and previous != digest:
        raise RuntimeError("conflicting SHA-256 entries for: " + raw_path)
    result[raw_path] = digest


def _checkpoint_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Collect every supported input-hash declaration from one checkpoint."""

    result: dict[str, str] = {}
    found = False
    for key in ("sha256", "input_hashes", "input_sha256"):
        mapping = source.get(key)
        if mapping is None:
            continue
        if not isinstance(mapping, Mapping) or not mapping:
            raise RuntimeError(f"source checkpoint has invalid {key} mapping")
        found = True
        for name, expected in mapping.items():
            _add_hash(result, name, expected)
    files = source.get("files")
    if files is not None:
        if not isinstance(files, list) or not files:
            raise RuntimeError("source checkpoint has invalid files manifest")
        found = True
        for item in files:
            if not isinstance(item, Mapping):
                raise RuntimeError("invalid source checkpoint file entry")
            _add_hash(result, item.get("path"), item.get("sha256"))
    if not found or not result:
        raise RuntimeError("source checkpoint has no input hash mapping")
    return result


def _require_brief_checkpoint_hash(
    brief: Mapping[str, Any], source_path: Path
) -> None:
    """If the brief pins the checkpoint file itself, verify that pin too."""

    checkpoint = brief.get("source_checkpoint")
    if not isinstance(checkpoint, Mapping):
        return
    expected = checkpoint.get("sha256")
    if isinstance(expected, str):
        normalized = expected.lower()
        if len(normalized) != 64 or any(
            character not in "0123456789abcdef" for character in normalized
        ):
            raise RuntimeError("pending brief has invalid source checkpoint SHA")
        if _sha(source_path) != normalized:
            raise RuntimeError("source checkpoint SHA mismatch in pending brief")


def _preflight() -> tuple[dict[str, Any], dict[str, Any], str]:
    """Validate pending state and all source-checkpoint input hashes."""

    if not BRIEF.is_file():
        raise RuntimeError("pending brief is missing: " + str(BRIEF))
    state = json.loads(BRIEF.read_bytes())
    if state.get("phase") != PENDING_PHASE:
        raise RuntimeError("wrong pending phase for compact scaled profiles")
    if state.get("component") != COMPONENT:
        raise RuntimeError("wrong pending component for compact scaled profiles")

    source_path = DIRECTORY / "source.json"
    if not source_path.is_file():
        raise RuntimeError("source checkpoint is missing: " + str(source_path))
    source = json.loads(source_path.read_bytes())
    if source.get("component") != COMPONENT:
        raise RuntimeError("source checkpoint component mismatch")
    _require_brief_checkpoint_hash(state, source_path)
    hashes = _checkpoint_hashes(source)
    for name, expected in hashes.items():
        path = _manifest_path(name)
        if not path.is_file():
            raise RuntimeError("frozen input is missing: " + name)
        if _sha(path) != expected:
            raise RuntimeError("frozen input changed: " + name)
    return state, source, _sha(source_path)


def main() -> int:
    # The runner may create only its receipt directory; it never creates or
    # edits the pending brief or source checkpoint.
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    entered_path = DIRECTORY / "entered.json"
    result_path = DIRECTORY / "result.json"
    if entered_path.exists() or result_path.exists():
        raise RuntimeError("compact scaled profiles already entered; never re-enter")

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
                "source_admissibility": "not_verified",
                "numeric_N_materialized": False,
                "runtime_model_metadata": None,
            },
        )
        spec = importlib.util.spec_from_file_location(
            "compact_scaled_profiles_new_leaf", LEAF
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        # Exactly one call into the new leaf fixture.  It performs only
        # synthetic caller-data checks and never loads the old source leaf.
        invocation_count = 1
        returned = leaf.verify_synthetic_fixture()
        if returned.get("source_value_computation") is not False:
            raise RuntimeError("leaf reported source-value computation")
        if returned.get("numeric_N_materialized") is not False:
            raise RuntimeError("leaf materialized a numerical N")
        if returned.get("source_admissibility") != "not_verified":
            raise RuntimeError("leaf changed source-admissibility status")
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
            "source_admissibility": "not_verified",
            "numeric_N_materialized": False,
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
            "source_admissibility": "not_verified",
            "numeric_N_materialized": False,
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
                "source_admissibility": "not_verified",
                "numeric_N_materialized": False,
                "certified": False,
            },
            sort_keys=True,
        )
    )
    return 0 if result["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
