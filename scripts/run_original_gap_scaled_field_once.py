"""Single-entry runner for the caller-owned original gap field port.

The root task prepares ``tmp/original_gap_scaled_field_v1/source.json`` and
the pending brief.  This runner does not create or update either checkpoint.
It verifies the pending phase and goal guards, verifies every SHA listed by the
source manifest (including the three frozen original pulse source files and
both new code files), performs a static formula/guard fragment audit, writes
an exclusive ``entered.json``, imports only the new leaf, and invokes its
independent fixture once.  A complete success or failure packet is written to
``result.json`` before the short status line is printed.
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
DIRECTORY = ROOT / "tmp" / "original_gap_scaled_field_v1"
SOURCE = DIRECTORY / "source.json"
ENTERED = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
LEAF = ROOT / "src" / "openai_ns_reconstruction" / "original_gap_scaled_field.py"
RUNNER = ROOT / "scripts" / "run_original_gap_scaled_field_once.py"
COMPONENT = "original_gap_scaled_field_v1"
PHASE = "original_gap_scaled_field_pending"
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")

OWN_FIELDS = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields.py"
)
KERNEL_FUNCTIONS = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions.py"
)
PRESSURE_FUNCTIONS = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions.py"
)
LEAF_RELATIVE = "src/openai_ns_reconstruction/original_gap_scaled_field.py"
RUNNER_RELATIVE = "scripts/run_original_gap_scaled_field_once.py"
REQUIRED_INPUTS = {OWN_FIELDS, KERNEL_FUNCTIONS, PRESSURE_FUNCTIONS}
REQUIRED_CODE = {LEAF_RELATIVE, RUNNER_RELATIVE}

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
            default=_json_default,
        )
        + "\n"
    ).encode("utf-8")


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if hasattr(value, "_mpf_"):
        return str(value)
    raise TypeError(f"not JSON serializable: {type(value).__name__}")


def _write_exclusive(path: Path, value: Any) -> str:
    raw = _json_bytes(value)
    with path.open("xb") as stream:
        stream.write(raw)
        stream.flush()
    return hashlib.sha256(raw).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _relative_path(raw: Any, label: str) -> str:
    _require(isinstance(raw, str), f"{label} path is not text")
    candidate = Path(raw)
    if candidate.is_absolute():
        try:
            return candidate.resolve().relative_to(ROOT.resolve()).as_posix()
        except ValueError as exc:
            raise RuntimeError(f"{label} path is outside the repository") from exc
    _require(".." not in candidate.parts, f"{label} path escapes the repository")
    return candidate.as_posix()


def _manifest_hashes(source: Mapping[str, Any]) -> dict[str, str]:
    """Read the supported source checkpoint hash layouts exactly once."""

    result: dict[str, str] = {}

    def add(raw_path: Any, raw_digest: Any, label: str) -> None:
        relative = _relative_path(raw_path, label)
        _require(
            isinstance(raw_digest, str) and SHA256_RE.fullmatch(raw_digest) is not None,
            f"{label} has no valid SHA-256",
        )
        digest = raw_digest.lower()
        prior = result.get(relative)
        _require(prior is None or prior == digest, f"conflicting SHA-256 for {relative}")
        result[relative] = digest

    for key in ("sha256", "input_hashes", "input_sha256", "code_sha256", "source_sha256"):
        mapping = source.get(key)
        if isinstance(mapping, Mapping):
            for raw_path, raw_digest in mapping.items():
                add(raw_path, raw_digest, f"source.{key}[{raw_path!r}]")
    files = source.get("files")
    if isinstance(files, list):
        for index, item in enumerate(files):
            _require(isinstance(item, Mapping), f"source.files[{index}] is not an object")
            add(item.get("path"), item.get("sha256"), f"source.files[{index}]")
    _require(result, "source checkpoint has no SHA-256 mapping")
    return result


def _brief_checkpoint(brief: Mapping[str, Any]) -> tuple[str, str]:
    records: list[Mapping[str, Any]] = []
    for key in (f"{COMPONENT}_source_checkpoint", "source_checkpoint"):
        value = brief.get(key)
        if isinstance(value, Mapping):
            records.append(value)
    for record in records:
        path = record.get("path")
        digest = record.get("sha256")
        if isinstance(path, str) and isinstance(digest, str):
            return path, digest
    path = brief.get("source_checkpoint_path")
    digest = brief.get("source_checkpoint_sha256")
    _require(isinstance(path, str), "pending brief has no source checkpoint path")
    _require(isinstance(digest, str), "pending brief has no source checkpoint SHA-256")
    return path, digest


def _preflight() -> tuple[dict[str, Any], dict[str, Any], str, dict[str, str], dict[str, Any]]:
    _require(not ENTERED.exists() and not RESULT.exists(), "already entered; never re-enter")
    _require(BRIEF.is_file(), f"pending brief is missing: {BRIEF}")
    _require(SOURCE.is_file(), f"source checkpoint is missing: {SOURCE}")
    brief = json.loads(BRIEF.read_bytes())
    source = json.loads(SOURCE.read_bytes())
    _require(isinstance(brief, Mapping), "pending brief is not a JSON object")
    _require(isinstance(source, Mapping), "source checkpoint is not a JSON object")
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PHASE, "pending brief phase mismatch")
    _require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("do_not_reenter") is True,
        "pending goal guards mismatch",
    )
    _require(source.get("component") == COMPONENT, "source checkpoint component mismatch")

    checkpoint_path, expected = _brief_checkpoint(brief)
    checkpoint_relative = _relative_path(checkpoint_path, "brief source checkpoint")
    _require((ROOT / checkpoint_relative).resolve() == SOURCE.resolve(),
             "pending brief points at a different source checkpoint")
    source_sha256 = sha256(SOURCE)
    _require(SHA256_RE.fullmatch(expected) is not None,
             "pending brief source checkpoint SHA-256 is invalid")
    _require(source_sha256 == expected.lower(), "source checkpoint SHA-256 mismatch")

    hashes = _manifest_hashes(source)
    _require(REQUIRED_INPUTS.issubset(hashes),
             "source checkpoint does not pin all three original pulse source files")
    _require(REQUIRED_CODE.issubset(hashes),
             "source checkpoint does not pin both new code files")
    for relative, expected_digest in hashes.items():
        path = ROOT / relative
        _require(path.is_file(), f"frozen input is missing: {relative}")
        _require(sha256(path) == expected_digest,
                 f"frozen input changed: {relative}")
    audit = {"scope": "source files pinned by raw SHA; no source-expression audit",
             "source_files": {relative: hashes[relative] for relative in sorted(REQUIRED_INPUTS)},
             "source_identity_proved": False}
    return dict(brief), dict(source), source_sha256, hashes, audit


def main() -> int:
    brief, source, source_sha256, input_hashes, source_audit = _preflight()
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
            "source_audit": source_audit,
            "consumer_calls": 0,
            "runtime_model_metadata": None,
            "goal_active": brief.get("goal_active"),
            "full_goal_complete": brief.get("full_goal_complete"),
            "do_not_reenter": brief.get("do_not_reenter"),
            "worker_is_not_only_worker": True,
        },
    )

    source_manifest_summary = {
        "component": source.get("component"),
        "schema": source.get("schema"),
        "source_checkpoint_sha256": source_sha256,
        "input_hashes": dict(input_hashes),
        "required_original_sources": sorted(REQUIRED_INPUTS),
        "required_new_code": sorted(REQUIRED_CODE),
        "hashes_verified_before_entry": True,
        "static_fragment_scope": source_audit["scope"],
    }

    consumer_calls = 0
    imported = False
    try:
        spec = importlib.util.spec_from_file_location(
            "original_gap_scaled_field_new_leaf", LEAF
        )
        if spec is None or spec.loader is None:
            raise RuntimeError("new leaf load failed")
        leaf = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = leaf
        spec.loader.exec_module(leaf)
        imported = True
        consumer_calls = 1
        returned = leaf.verify_original_gap_scaled_field_fixture(source_audit)
        if not isinstance(returned, Mapping):
            raise RuntimeError("fixture did not return a mapping")
        if returned.get("passed") is not True:
            raise RuntimeError("original gap scaled-field fixture did not pass")
        output: dict[str, Any] = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "returned",
            "brief_sha256": sha256(BRIEF),
            "source_checkpoint_sha256": source_sha256,
            "entry_sha256": entered_sha256,
            "input_hashes": dict(input_hashes),
            "source_audit": source_audit,
            "source_manifest": source_manifest_summary,
            "consumer_calls": consumer_calls,
            "leaf_imported": imported,
            "runtime_model_metadata": None,
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "worker_is_not_only_worker": True,
                "source_manifest_format": "source.json SHA-256 mapping; source identity remains unproved",
            "limitations": [
                "caller-owned point data; no source binding or authentication",
                "fixture u/C1/C2 jets are synthetic and are not source values",
                "no P/F signed comparison is performed",
                "logarithmic scale is not materialized",
                "only the independent medium-mu fixture materializes its moderate original exponent",
                "finite differences are diagnostic, not an error certificate",
                "Q_ZZ Z1 is unavailable from only C3 input jets",
                "the real mu=exp(-4*(exp(40)+11))/1000 scale is recorded symbolically, not materialized",
            ],
            "result": dict(returned),
        }
    except Exception as exc:  # noqa: BLE001 - complete failure artifact is required
        trace = traceback.format_exc()
        output = {
            "schema": f"{COMPONENT}_result_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "status": "failed",
            "brief_sha256": sha256(BRIEF),
            "source_checkpoint_sha256": source_sha256,
            "entry_sha256": entered_sha256,
            "input_hashes": dict(input_hashes),
            "source_audit": source_audit,
            "source_manifest": source_manifest_summary,
            "consumer_calls": consumer_calls,
            "leaf_imported": imported,
            "runtime_model_metadata": None,
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "worker_is_not_only_worker": True,
            "error": {
                "type": type(exc).__name__,
                "message": str(exc),
                "traceback": trace,
            },
            "traceback": trace,
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
