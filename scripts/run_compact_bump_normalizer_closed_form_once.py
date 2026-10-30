"""Enter and invoke the closed-form compact-bump leaf exactly once.

All pending-state and frozen-input checks happen before the exclusive
``entered.json`` write.  The new leaf is then imported and its fixture is
called exactly once.  The runner never calls a numerical integrator or an
ancestor normalizer, and it writes a complete result receipt (including a
traceback on leaf failure) without modifying the pending brief or source
checkpoint.
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
DIRECTORY = ROOT / "tmp/compact_bump_normalizer_closed_form_v1"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "src/openai_ns_reconstruction/compact_bump_normalizer_closed_form.py"
COMPONENT = "compact_bump_normalizer_closed_form_v1"
PENDING_PHASE = "compact_bump_normalizer_closed_form_pending"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_default(value: Any) -> Any:
    """Keep high-precision mpmath values as decimal strings in receipts."""

    if hasattr(value, "_mpf_"):
        return str(value)
    if hasattr(value, "_mpc_"):
        return str(value)
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def _exclusive_json(path: Path, value: Mapping[str, Any]) -> str:
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


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _manifest_hashes(checkpoint: Mapping[str, Any]) -> dict[str, str]:
    """Read either common frozen-input manifest shape, without weakening it."""

    hashes = checkpoint.get("sha256")
    if not isinstance(hashes, Mapping) or not hashes:
        hashes = checkpoint.get("input_hashes")
    if isinstance(hashes, Mapping) and hashes:
        result = {}
        for name, expected in hashes.items():
            relative = str(name)
            digest = str(expected).lower()
            _require(
                len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest),
                "invalid SHA-256 in source checkpoint: " + relative,
            )
            result[relative] = digest
        return result

    files = checkpoint.get("files")
    if isinstance(files, list) and files:
        result = {}
        for item in files:
            _require(isinstance(item, Mapping), "invalid frozen-file hash entry")
            relative = str(item.get("path"))
            digest = str(item.get("sha256", "")).lower()
            _require(relative not in {"", "None"}, "frozen-file path is missing")
            _require(
                len(digest) == 64
                and all(character in "0123456789abcdef" for character in digest),
                "invalid SHA-256 in frozen-file manifest: " + relative,
            )
            result[relative] = digest
        return result
    raise RuntimeError("checkpoint has no frozen-input hash mapping")


def _manifest_path(raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else ROOT / path


def _validate_checkpoint_hashes(
    source: Mapping[str, Any], capsule: Mapping[str, Any] | None
) -> tuple[str, str | None, int]:
    """Validate every declared frozen hash before writing ``entered``."""

    _require(source.get("component") == COMPONENT, "source checkpoint component mismatch")
    source_hashes = _manifest_hashes(source)
    for name, expected in source_hashes.items():
        path = _manifest_path(name)
        _require(path.is_file(), "missing frozen input: " + name)
        _require(_sha(path) == expected, "frozen input changed: " + name)

    capsule_sha256: str | None = None
    validated_inputs = len(source_hashes)
    declared_capsule_sha = source.get("capsule_sha256")
    if declared_capsule_sha is not None:
        capsule_path = DIRECTORY / "capsule.json"
        _require(capsule_path.is_file(), "declared capsule checkpoint is missing")
        capsule_sha256 = _sha(capsule_path)
        _require(
            capsule_sha256 == str(declared_capsule_sha).lower(),
            "capsule checkpoint SHA mismatch",
        )

    if capsule is not None:
        _require(capsule.get("component") == COMPONENT, "capsule component mismatch")
        capsule_hashes = capsule.get("input_sha256")
        if not isinstance(capsule_hashes, Mapping) or not capsule_hashes:
            capsule_hashes = capsule.get("sha256")
        if isinstance(capsule_hashes, Mapping) and capsule_hashes:
            for name, expected in capsule_hashes.items():
                relative = str(name)
                digest = str(expected).lower()
                _require(
                    len(digest) == 64
                    and all(character in "0123456789abcdef" for character in digest),
                    "invalid SHA-256 in capsule input map: " + relative,
                )
                path = _manifest_path(relative)
                _require(path.is_file(), "missing capsule input: " + relative)
                _require(_sha(path) == digest, "capsule input changed: " + relative)
                validated_inputs += 1
    return _sha(DIRECTORY / "source.json"), capsule_sha256, validated_inputs


def _preflight() -> tuple[dict[str, Any], dict[str, Any], dict[str, Any] | None, str, str | None, int]:
    """Validate pending phase/component and all frozen hashes pre-entry."""

    _require(BRIEF.is_file(), "pending brief is missing: " + str(BRIEF))
    state = json.loads(BRIEF.read_bytes())
    _require(
        state.get("phase") == PENDING_PHASE,
        "wrong pending phase for compact bump closed form",
    )
    _require(state.get("component") == COMPONENT, "wrong pending component")
    if "goal_active" in state:
        _require(state.get("goal_active") is True, "pending goal is not active")
    if "full_goal_complete" in state:
        _require(
            state.get("full_goal_complete") is False,
            "pending brief already marks the full goal complete",
        )
    if "do_not_reenter" in state:
        _require(state.get("do_not_reenter") is True, "pending brief lacks no-reentry guard")

    source_path = DIRECTORY / "source.json"
    _require(source_path.is_file(), "source checkpoint is missing: " + str(source_path))
    _require(state.get("source_checkpoint", {}).get("sha256") == _sha(source_path), "pending checkpoint hash mismatch")
    source = json.loads(source_path.read_bytes())
    capsule_path = DIRECTORY / "capsule.json"
    capsule = json.loads(capsule_path.read_bytes()) if capsule_path.is_file() else None
    source_sha256, capsule_sha256, input_count = _validate_checkpoint_hashes(source, capsule)
    return state, source, capsule, source_sha256, capsule_sha256, input_count


def main() -> int:
    DIRECTORY.mkdir(parents=True, exist_ok=True)
    entered_path = DIRECTORY / "entered.json"
    result_path = DIRECTORY / "result.json"
    if entered_path.exists() or result_path.exists():
        raise RuntimeError(
            "compact bump normalizer closed form already entered; never re-enter"
        )

    state, source, capsule, source_sha256, capsule_sha256, input_count = _preflight()
    entry_sha256 = _exclusive_json(
        entered_path,
        {
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "entered": True,
            "source_sha256": source_sha256,
            "capsule_sha256": capsule_sha256,
            "validated_frozen_input_count": input_count,
            "verification_function": "verify_special_function_fixture",
            "verification_invocations": 0,
            "consumer_calls": 0,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_value_computation": False,
            "numeric_evaluation": True,
            "certified": False,
            "runtime_model_metadata": None,
        },
    )

    invocation_count = 0
    try:
        spec = importlib.util.spec_from_file_location(
            "compact_bump_normalizer_closed_form_new_leaf", LEAF
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        invocation_count = 1
        returned = leaf.verify_special_function_fixture(ROOT)
        _require(returned.get("schema") == COMPONENT, "leaf schema mismatch")
        _require(returned.get("checks_passed") is True, "special-function fixture did not pass")
        _require(returned.get("certified") is False, "leaf claimed a certified value")
        _require(
            returned.get("source_provider_graph_same_proof") is False,
            "leaf claimed source-provider graph same proof",
        )
        result: dict[str, Any] = {
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "status": "returned",
            "source_sha256": source_sha256,
            "capsule_sha256": capsule_sha256,
            "entry_sha256": entry_sha256,
            "validated_frozen_input_count": input_count,
            "verification_function": "verify_special_function_fixture",
            "verification_invocations": invocation_count,
            "consumer_calls": invocation_count,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_value_computation": False,
            "numeric_evaluation": True,
            "certified": False,
            "runtime_model_metadata": None,
            "pending_state_component": state.get("component"),
            "result": returned,
        }
    except Exception:
        result = {
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "status": "failed",
            "source_sha256": source_sha256,
            "capsule_sha256": capsule_sha256,
            "entry_sha256": entry_sha256,
            "validated_frozen_input_count": input_count,
            "verification_function": "verify_special_function_fixture",
            "verification_invocations": invocation_count,
            "consumer_calls": invocation_count,
            "ancestor_math_calls": 0,
            "source_values_computed": False,
            "source_value_computation": False,
            "numeric_evaluation": True,
            "certified": False,
            "runtime_model_metadata": None,
            "traceback": traceback.format_exc(),
        }

    result_sha256 = _exclusive_json(result_path, result)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "phase": PENDING_PHASE,
                "status": result["status"],
                "result_sha256": result_sha256,
                "verification_invocations": invocation_count,
                "consumer_calls": invocation_count,
                "ancestor_math_calls": 0,
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
