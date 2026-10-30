"""Run the authenticated historical-frame owner receipt exactly once.

The root task owns and prepares ``tmp/historical_frame_owner_binding_v1`` and
its frozen ``source.json`` checkpoint.  This runner never creates that
checkpoint or a cache.  It verifies the pending brief, source checkpoint,
all six saved inputs, and both code files before writing an exclusive entry
receipt.  Only then is the new leaf imported and called once.  A complete
success or failure packet is persisted before anything is printed.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import re
import sys
import traceback
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "tmp" / "historical_frame_owner_binding_v1"
SOURCE = DIRECTORY / "source.json"
ENTERED = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "src" / "openai_ns_reconstruction" / "historical_frame_owner_binding.py"
COMPONENT = "historical_frame_owner_binding_v1"
PHASE = "prepared_historical_frame_owner_binding"
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

REQUIRED_INPUTS = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source.json.gz",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.json.gz",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame.json.gz",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source_check.json",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions_check.json",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame_check.json",
}
CODE_PATHS = {
    "src/openai_ns_reconstruction/historical_frame_owner_binding.py",
    "scripts/run_historical_frame_owner_binding_once.py",
}


def sha256(path: Path) -> str:
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


def _write_exclusive(path: Path, value: Any) -> str:
    raw = _json_bytes(value)
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _relative_path(raw: Any, label: str) -> str:
    _require(isinstance(raw, str), f"{label} path is not text")
    path = Path(raw)
    if path.is_absolute():
        try:
            return path.resolve().relative_to(ROOT.resolve()).as_posix()
        except ValueError as exc:
            raise RuntimeError(f"{label} path is outside the repository") from exc
    _require(".." not in path.parts, f"{label} path escapes the repository")
    return path.as_posix()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Read supported source checkpoint hash layouts without broad traversal."""

    result: dict[str, str] = {}

    def add(raw_path: Any, raw_digest: Any, label: str) -> None:
        relative = _relative_path(raw_path, label)
        _require(isinstance(raw_digest, str) and SHA256_RE.fullmatch(raw_digest) is not None,
                 f"{label} has no valid SHA-256")
        digest = raw_digest.lower()
        prior = result.get(relative)
        _require(prior is None or prior == digest, f"conflicting SHA-256 for {relative}")
        result[relative] = digest

    for key in ("sha256", "input_hashes", "input_sha256", "code_sha256"):
        mapping = source.get(key)
        if isinstance(mapping, Mapping):
            for raw_path, raw_digest in mapping.items():
                add(raw_path, raw_digest, f"source.{key}[{raw_path!r}]")

    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            add(item.get("path"), item.get("sha256"), f"source.files[{index}]")

    return result


def _brief_checkpoint_hash(brief: Mapping[str, Any]) -> tuple[str, str | None]:
    """Return the brief's source checkpoint path and hash."""

    records: list[Mapping[str, Any]] = []
    for key in (f"{COMPONENT}_source_checkpoint", "source_checkpoint"):
        value = brief.get(key)
        if isinstance(value, Mapping):
            records.append(value)
    for record in records:
        raw_path = record.get("path")
        raw_sha = record.get("sha256")
        if isinstance(raw_path, str):
            return raw_path, raw_sha if isinstance(raw_sha, str) else None

    raw_path = brief.get("source_checkpoint_path")
    raw_sha = brief.get("source_checkpoint_sha256")
    _require(isinstance(raw_path, str), "pending brief has no source checkpoint path")
    return raw_path, raw_sha if isinstance(raw_sha, str) else None


def _preflight() -> tuple[dict[str, Any], dict[str, Any], str, dict[str, str]]:
    _require(not ENTERED.exists() and not RESULT.exists(), "already entered; never re-enter")
    _require(BRIEF.is_file(), f"pending brief is missing: {BRIEF}")
    _require(SOURCE.is_file(), f"source checkpoint is missing: {SOURCE}")
    brief = json.loads(BRIEF.read_bytes())
    source = json.loads(SOURCE.read_bytes())
    _require(isinstance(brief, Mapping), "pending brief is not a JSON object")
    _require(isinstance(source, Mapping), "source checkpoint is not a JSON object")
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    _require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False
             and brief.get("do_not_reenter") is True, "pending goal guards mismatch")
    _require(source.get("component") == COMPONENT, "source checkpoint component mismatch")

    raw_checkpoint_path, checkpoint_expected = _brief_checkpoint_hash(brief)
    checkpoint_relative = _relative_path(raw_checkpoint_path, "brief source checkpoint")
    _require((ROOT / checkpoint_relative).resolve() == SOURCE.resolve(),
             "pending brief points at a different source checkpoint")
    checkpoint_actual = sha256(SOURCE)
    _require(isinstance(checkpoint_expected, str) and SHA256_RE.fullmatch(checkpoint_expected) is not None,
             "pending brief source checkpoint SHA-256 is invalid")
    _require(checkpoint_actual == checkpoint_expected.lower(), "source checkpoint SHA-256 mismatch")

    hashes = _manifest_hashes(source)
    _require(REQUIRED_INPUTS.issubset(hashes), "source checkpoint does not pin all six saved inputs")
    _require(CODE_PATHS.issubset(hashes), "source checkpoint does not pin both code files")
    for relative, expected in hashes.items():
        path = ROOT / relative
        _require(path.is_file(), f"frozen input is missing: {relative}")
        _require(sha256(path) == expected, f"frozen input changed: {relative}")
    return dict(brief), dict(source), checkpoint_actual, hashes


def main() -> int:
    brief, source, source_sha256, input_hashes = _preflight()
    entered_sha256 = _write_exclusive(
        ENTERED,
        {
            "schema": f"{COMPONENT}_entered_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "brief_sha256": sha256(BRIEF),
            "source_checkpoint_sha256": source_sha256,
            "input_hashes": dict(input_hashes),
            "consumer_calls": 0,
            "runtime_model_metadata": None,
        },
    )

    consumer_calls = 0
    try:
        spec = importlib.util.spec_from_file_location(
            "historical_frame_owner_binding_new_leaf", LEAF
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("new leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        consumer_calls = 1
        value = leaf.recover(ROOT, source)
        output: dict[str, Any] = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "brief_sha256": sha256(BRIEF),
            "source_checkpoint_sha256": source_sha256,
            "entry_sha256": entered_sha256,
            "input_hashes": dict(input_hashes),
            "consumer_calls": consumer_calls,
            "runtime_model_metadata": None,
            "error": None,
            "result": value,
        }
    except Exception as exc:  # noqa: BLE001 - complete failure artifact is required
        output = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "brief_sha256": sha256(BRIEF),
            "source_checkpoint_sha256": source_sha256,
            "entry_sha256": entered_sha256,
            "input_hashes": dict(input_hashes),
            "consumer_calls": consumer_calls,
            "runtime_model_metadata": None,
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
            "traceback": traceback.format_exc(),
        }

    result_sha256 = _write_exclusive(RESULT, output)
    print(json.dumps({
        "component": COMPONENT,
        "phase": PHASE,
        "status": output["status"],
        "result_sha256": result_sha256,
        "consumer_calls": consumer_calls,
        "runtime_model_metadata": None,
    }, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
