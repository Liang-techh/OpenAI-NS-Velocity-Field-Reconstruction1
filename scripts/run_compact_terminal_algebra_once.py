"""Single-entry runner for the compact terminal algebra leaf.

The pending state and every checkpoint hash are validated before loading the
leaf.  ``entered.json`` and the completed ``result.json`` are exclusive files:
an existing receipt is never overwritten and the leaf cannot be re-entered.
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
DIRECTORY = ROOT / "tmp/compact_terminal_algebra_v1"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "src/openai_ns_reconstruction/compact_terminal_algebra.py"
COMPONENT = "compact_terminal_algebra_v1"
PENDING_PHASE = "compact_terminal_algebra_pending"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _exclusive_json(path: Path, value: Mapping[str, Any]) -> str:
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _checkpoint_files(source: Mapping[str, Any]) -> dict[str, str]:
    hashes = source.get("sha256")
    if not isinstance(hashes, Mapping) or not hashes:
        raise RuntimeError("source checkpoint has no sha256 mapping")
    return {str(name): str(expected) for name, expected in hashes.items()}


def main() -> int:
    source_path = DIRECTORY / "source.json"
    entered_path = DIRECTORY / "entered.json"
    result_path = DIRECTORY / "result.json"
    if entered_path.exists() or result_path.exists():
        raise RuntimeError("compact terminal algebra already entered; never re-enter")
    state = json.loads(BRIEF.read_bytes())
    if state.get("phase") != PENDING_PHASE or state.get("component") != COMPONENT:
        raise RuntimeError("wrong pending state")
    source = json.loads(source_path.read_bytes())
    if source.get("component") != COMPONENT:
        raise RuntimeError("source checkpoint component mismatch")
    checkpoints = _checkpoint_files(source)
    for name, expected in checkpoints.items():
        path = ROOT / name
        if not path.is_file() or _sha(path) != expected:
            raise RuntimeError("frozen input changed: " + name)
    entry_sha = _exclusive_json(
        entered_path,
        {
            "component": COMPONENT,
            "phase": PENDING_PHASE,
            "entered": True,
            "source_sha256": _sha(source_path),
            "consumer_calls": 0,
            "runtime_model_metadata": None,
        },
    )
    try:
        spec = importlib.util.spec_from_file_location("compact_terminal_algebra_new_leaf", LEAF)
        if spec is None or spec.loader is None:
            raise RuntimeError("leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        returned = leaf.recover(ROOT)
        result: dict[str, Any] = {
            "component": COMPONENT,
            "status": "returned",
            "source_sha256": _sha(source_path),
            "entry_sha256": entry_sha,
            "consumer_calls": 1,
            "runtime_model_metadata": None,
            "result": returned,
        }
    except Exception:
        result = {
            "component": COMPONENT,
            "status": "failed",
            "source_sha256": _sha(source_path),
            "entry_sha256": entry_sha,
            "consumer_calls": 1,
            "runtime_model_metadata": None,
            "traceback": traceback.format_exc(),
        }
    result_sha = _exclusive_json(result_path, result)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "status": result["status"],
                "result_sha256": result_sha,
                "consumer_calls": 1,
                "runtime_model_metadata": None,
            },
            sort_keys=True,
        )
    )
    return 0 if result["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
