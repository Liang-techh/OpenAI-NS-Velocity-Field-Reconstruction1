"""Durable metadata driver for the native O3 physical binding leaf.

The prepare action authenticates frozen bytes, JSON literals, and Python ASTs
only.  The run action is the sole source consumer: it persists entry before a
single import/derive call, journals every two-argument callback immediately,
and saves the complete result before validation.  This file contains no
mathematical implementation.
"""

from __future__ import annotations

import argparse
import ast
from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
from typing import Any


REPO = Path(
    r"C:/Users/z5242/Documents/Codex/2026-09-13/"
    r"before-doing-any-project-work-verify/work/st073-transition-next"
)
WORK = Path(
    r"C:/Users/z5242/Documents/Codex/2026-10-16/"
    r"https-github-com-liang-techh-openai/work"
)

COMPONENT = "native_O3_total_radius_velocity_pressure_binding"
INPUT_SCHEMA = f"{COMPONENT}_input_v1"
MANIFEST_SCHEMA = f"{COMPONENT}_manifest_v1"
SOURCE_SCHEMA = f"{COMPONENT}_source_v1"
ENTERED_SCHEMA = f"{COMPONENT}_entered_v1"
RESULT_SCHEMA = f"{COMPONENT}_result_checkpoint_v1"
JOURNAL_SCHEMA = f"{COMPONENT}_journal_v1"

PARENT_REPORT = REPO / (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_O3_finite_N_ast_flux_semantics.json"
)
PARENT_REPORT_SHA256 = "962a4953bde3a9f026d84d36e23eea3bb678924a3dda0201c6e513016397dd3d"
PARENT_REPORT_INTERNAL_PARENT_COMMIT = "4f049dfef1f70cbee9b09953cb882a09aa81fbd5"
CURRENT_HEAD_COMMIT = "32c9a65aa8cc17f893a94052644c4a70a728eab8"
PARENT_COMMIT = CURRENT_HEAD_COMMIT
PARENT_RESULT = REPO / "tmp/native_O3_finite_N_ast_flux_semantics/result.json"
PARENT_RESULT_SHA256 = "d4cf1b39195a608e32197c58d572f86f0d974ba940ef2ff577295fdfc95b092d"
PARENT_SOURCE_PIN_COUNT = 2325
PARENT_CHECKPOINT_COUNT = 288
PARENT_FAILED_ARTIFACT_COUNT = 190
PARENT_PROTECTED_FILE_COUNT = 16
PARENT_DIRTY_PATH_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
PUBLISHED_CHECKPOINT_COUNT = 292
FAILED_COUNT_IF_NEW_FAILURE = 197
EXPECTED_SOURCE_UNION_COUNT = 2330
SOURCE_ARCHIVE_SHA256 = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"

MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
KERNEL_RELATIVE = Path(
    "src/openai_ns_reconstruction/"
    "native_O3_total_radius_velocity_pressure_binding.py"
)
RUNNER_RELATIVE = Path("scripts/run_native_O3_total_radius_velocity_pressure_binding.py")
KERNEL_PATH = REPO / KERNEL_RELATIVE
RUNNER_PATH = REPO / RUNNER_RELATIVE
INPUT_PATH = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF_PATH = WORK / f"{COMPONENT}_brief_2026_10_27.json"
CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_PATH = CHECKPOINT_DIR / "source.json"
ENTERED_PATH = CHECKPOINT_DIR / "entered.json"
JOURNAL_PATH = CHECKPOINT_DIR / "journal.jsonl"
RESULT_PATH = CHECKPOINT_DIR / "result.json"

WRAPPER_PATHS = {
    "saved_flux_result": PARENT_RESULT,
    "source_literals": WORK / "native_O3_density_flux_raw_literals_2026_10_27_v2.json",
    "physical_source_capsule": WORK / "native_O3_physical_source_capsule_2026_10_27_v2.json",
}
WRAPPER_NAMES = tuple(WRAPPER_PATHS)
WRAPPER_SHA256 = {
    "saved_flux_result": PARENT_RESULT_SHA256,
    "source_literals": "ae3a6afbdc8de9c557dff18b3552a6b1c6c1d902decca5298b0b6adbcc216b59",
    "physical_source_capsule": "90ce8bbacd274cf7ea022afd976799944b60709fef69162a6736a891fcef7788",
}

EXECUTION_RULES = [
    "Never rerun any entered ancestor mathematical/helper family, including failed leaves, old builders/providers/normalizers/checkers, source-value reducers, quadrature or recurrence. Static source/AST/literal JSON/hash inspection is allowed.",
    "Before a NEW mathematical/source-expression consumer: freeze code and input, persist source checkpoint and pending brief, persist entry before invocation, invoke once, save complete result before review. Preserve a failure and repair in a separate frozen leaf; never reenter.",
    "Authenticate raw gzip SHA BEFORE decoding. Never decode docs/continuation/RESUME_FIRST_SWITCH_2026_10_17.json or the huge Rm_loop_density_identity.json; raw SHA only.",
    "Never restart stopped session 17273. Retain nonzero incoming histories and ordinary Z1. Bounds, midpoints and interval overlap are not function values; remove old numerical errors only with a same-source packet theorem.",
    "Use current tree and GitHub; preserve all dirty and untracked files. Commit/push only exact accepted scoped files. Goal remains ACTIVE / INCOMPLETE.",
    "Children: luna_max_worker, gpt-5.6-luna/max fixed role, fork_turns:none, complete brief. Runtime model metadata unavailable means null; never claim it was verified.",
    "A worker function that checks mathematical source identities counts as a source consumer even when it consumes static receipts or is called a pure self-check. Workers return unexecuted code; root alone freezes code and inputs, records entry, invokes once and saves the complete result before review. AST parsing and metadata-only audits remain permitted. Do not execute or reconstruct the rejected memory draft result.",
]


def require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def valid_sha(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(
        character in "0123456789abcdefABCDEF" for character in value
    )


def wire(path: str | Path) -> str:
    return Path(path).as_posix()


def encode(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        + b"\n"
    )


def read_bytes(path: Path, expected: str | None = None) -> bytes:
    require(path.is_file(), f"missing artifact: {wire(path)}")
    raw = path.read_bytes()
    if expected is not None:
        require(sha_bytes(raw) == expected.lower(), f"raw SHA mismatch: {wire(path)}")
    return raw


def read_json(path: Path, expected: str | None = None) -> tuple[bytes, dict[str, Any]]:
    raw = read_bytes(path, expected)
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict), f"JSON object required: {wire(path)}")
    return raw, value


def write_new(path: Path, value: bytes | Any) -> bytes:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("xb") as handle:
        handle.write(raw)
        handle.flush()
        os.fsync(handle.fileno())
    return raw


def write_replace(path: Path, value: bytes | Any) -> bytes:
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("wb") as handle:
        handle.write(raw)
        handle.flush()
        os.fsync(handle.fileno())
    return raw


def _candidates(value: str | Path) -> list[Path]:
    path = Path(value)
    if path.is_absolute():
        return [path]
    return [
        Path.cwd() / path,
        REPO / path,
        WORK / path,
        REPO / "experiments/root_st073" / path,
    ]


def resolve_artifact(value: str | Path, expected: str | None = None) -> Path:
    candidates = _candidates(value)
    if expected is not None:
        expected = expected.lower()
        for candidate in candidates:
            if candidate.is_file() and sha_bytes(candidate.read_bytes()) == expected:
                return candidate
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    raise ValueError(f"missing artifact: {value}")


def artifact(path: Path, raw: bytes, *, hash_only: bool = False) -> dict[str, Any]:
    value: dict[str, Any] = {
        "path": wire(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "raw_bytes_exact": True,
        "raw_readonly": True,
    }
    if hash_only:
        value.update({"hash_only": True, "literal_embedded": False})
    return value


def authenticate_spec(
    spec: Mapping[str, Any],
    label: str,
    *,
    expected_path: Path | None = None,
) -> tuple[Path, bytes, str]:
    path_value = spec.get("path", spec.get("artifact_path"))
    digest = spec.get("sha256", spec.get("raw_sha256"))
    require(isinstance(path_value, str) and path_value, f"{label} path missing")
    require(valid_sha(digest), f"{label} SHA malformed")
    path = resolve_artifact(path_value, digest)
    if expected_path is not None:
        require(path.resolve() == expected_path.resolve(), f"{label} path is not frozen")
    raw = read_bytes(path, digest)
    return path, raw, digest.lower()


def parse_ast(path: Path, raw: bytes, label: str) -> dict[str, Any]:
    try:
        tree = ast.parse(raw.decode("utf-8"), filename=str(path), mode="exec")
    except (SyntaxError, UnicodeDecodeError) as exc:
        raise ValueError(f"{label} AST parse failed: {exc}") from exc
    definitions = [
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
    ]
    return {
        "path": wire(path),
        "sha256": sha_bytes(raw),
        "ast_parse": "ok",
        "node_count": sum(1 for _ in ast.walk(tree)),
        "top_level_definitions": definitions,
    }


def git_head() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def load_active_brief() -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    raw, brief = read_json(MAIN_BRIEF)
    require(len(raw) < 100_000, "active main brief exceeds 100 KB")
    require(brief.get("execution_rules") == EXECUTION_RULES, "global execution rules changed")
    require(brief.get("execution_rules_preserved_exact") is True, "execution rules are not sealed")
    require(brief.get("goal_active") is True, "main brief is not active")
    require(brief.get("full_goal_complete") is False, "main brief is already complete")
    accepted = brief.get("last_accepted_mathematical_state")
    require(isinstance(accepted, Mapping), "last accepted mathematical state missing")
    require(
        accepted.get("component") == "native_O3_finite_N_ast_flux_semantics",
        "accepted finite-N state changed",
    )
    require(
        accepted.get("commit") == CURRENT_HEAD_COMMIT
        and accepted.get("head") == CURRENT_HEAD_COMMIT,
        "accepted finite-N commit changed",
    )
    require(
        accepted.get("report_sha256") == PARENT_REPORT_SHA256,
        "accepted finite-N report changed",
    )
    require(
        accepted.get("source_only") is True
        and accepted.get("metadata_only") is False
        and accepted.get("full_goal_complete") is False,
        "accepted finite-N scope changed",
    )
    baseline = brief.get("accepted_progress")
    require(
        isinstance(baseline, Mapping)
        and baseline.get("commit") == "8f36f19e3c524b424a9286c46cc8e0b855f75307",
        "physical baseline changed",
    )
    context = {
        "execution_rules": list(brief["execution_rules"]),
        "accepted_mathematical_state": dict(accepted),
        "physical_baseline": dict(baseline),
        "goal_active": brief["goal_active"],
        "full_goal_complete": brief["full_goal_complete"],
    }
    return raw, brief, context


def brief_preimage(raw: bytes) -> tuple[Path, bytes]:
    path = WORK / f"{COMPONENT}_main_resume_brief_preimage_{sha_bytes(raw)[:12]}.json"
    require(not path.exists(), "main brief preimage already exists")
    return path, write_new(path, raw)


def update_brief(raw: bytes, brief: Mapping[str, Any], patch: Mapping[str, Any]) -> dict[str, Any]:
    require(brief.get("execution_rules") == EXECUTION_RULES, "execution rules changed before brief update")
    updated = dict(brief)
    for key, value in patch.items():
        if key == "in_progress" and isinstance(updated.get(key), dict) and isinstance(value, Mapping):
            merged = dict(updated[key])
            merged.update(dict(value))
            updated[key] = merged
        else:
            updated[key] = value
    updated["execution_rules"] = list(EXECUTION_RULES)
    updated["execution_rules_preserved_exact"] = True
    new_raw = encode(updated)
    require(len(new_raw) < 100_000, "updated main brief exceeds 100 KB")
    write_replace(MAIN_BRIEF, new_raw)
    return {
        "path": wire(MAIN_BRIEF),
        "sha256": sha_bytes(new_raw),
        "bytes": len(new_raw),
        "preimage_sha256": sha_bytes(raw),
    }


def map_hashes(value: Any, label: str) -> dict[str, str]:
    require(isinstance(value, Mapping), f"{label} must be a mapping")
    result: dict[str, str] = {}
    for path, digest in value.items():
        require(isinstance(path, str) and path, f"{label} path malformed")
        require(valid_sha(digest), f"{label} SHA malformed: {path}")
        result[path] = digest.lower()
    return result


def authenticate_hash_map(value: Any, label: str) -> dict[str, str]:
    mapping = map_hashes(value, label)
    for path_value, digest in mapping.items():
        read_bytes(resolve_artifact(path_value, digest), digest)
    return mapping


def load_parent() -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    # Preserve the parent as a hash-only pointer; never recursively embed its
    # approximately 19 MB report or callback/result payload.
    report_raw, report = read_json(PARENT_REPORT, PARENT_REPORT_SHA256)
    require(
        report.get("component") == "native_O3_finite_N_ast_flux_semantics",
        "parent component changed",
    )
    require(report.get("passed") is True and report.get("attempt_status") == "passed", "parent status changed")
    require(
        report.get("source_only") is True
        and report.get("metadata_only") is False
        and report.get("full_goal_complete") is False,
        "parent scope changed",
    )
    require(report.get("conditional") is True, "parent conditional seal changed")
    require(
        report.get("math_source_identity_accepted") is True
        and report.get("mathematical_proof_accepted") is True,
        "parent accepted proof state changed",
    )
    require(report.get("physical_layer_completed") is False, "parent physical state changed")
    require(report.get("parent_commit") == PARENT_REPORT_INTERNAL_PARENT_COMMIT, "parent lineage changed")
    maps = {
        "source_pins": authenticate_hash_map(report.get("source_pins"), "parent source pins"),
        "protected_file_hashes": authenticate_hash_map(
            report.get("protected_file_hashes"), "parent protected files"
        ),
        "legacy_dirty_path_hashes": authenticate_hash_map(
            report.get("legacy_dirty_path_hashes"), "parent dirty files"
        ),
        "preserved_failed_artifact_hashes_union": authenticate_hash_map(
            report.get("preserved_failed_artifact_hashes_union"), "parent failed files"
        ),
    }
    require(len(maps["source_pins"]) == PARENT_SOURCE_PIN_COUNT, "parent source-pin count changed")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT, "parent source_pin_count changed")
    require(report.get("source_union_count") == PARENT_SOURCE_PIN_COUNT, "parent source union changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT, "parent checkpoint count changed")
    require(
        len(maps["preserved_failed_artifact_hashes_union"]) == PARENT_FAILED_ARTIFACT_COUNT,
        "parent failed artifact count changed",
    )
    require(
        report.get("preserved_failed_artifact_count") == PARENT_FAILED_ARTIFACT_COUNT,
        "parent failed artifact field changed",
    )
    require(
        len(maps["protected_file_hashes"]) == PARENT_PROTECTED_FILE_COUNT,
        "parent protected count changed",
    )
    require(
        len(maps["legacy_dirty_path_hashes"]) == PARENT_DIRTY_PATH_COUNT,
        "parent dirty count changed",
    )
    false_flags = report.get("parent_false_flags")
    false_values = report.get("parent_false_flag_values")
    require(
        isinstance(false_flags, list) and len(false_flags) == PARENT_FALSE_FLAG_COUNT,
        "parent false-flag ledger changed",
    )
    require(
        false_values == {name: False for name in false_flags},
        "parent false-flag values changed",
    )
    new_flags = report.get("new_flags")
    require(
        isinstance(new_flags, Mapping)
        and all(isinstance(value, bool) for value in new_flags.values()),
        "parent new flags malformed",
    )
    retained: dict[str, Any] = {}
    for key in (
        "preserved_error_cap_metadata",
        "retained_construction_data",
        "retained_integral_term_inventory",
        "retained_cap_provenance",
    ):
        require(key in report, f"parent retained ledger missing: {key}")
        retained[key] = report[key]
    parent_counts = {
        "source_pins": len(maps["source_pins"]),
        "checkpoints": int(report["checkpoint_count"]),
        "failed": len(maps["preserved_failed_artifact_hashes_union"]),
        "protected": len(maps["protected_file_hashes"]),
        "dirty": len(maps["legacy_dirty_path_hashes"]),
        "false_flags": len(false_flags),
        "parent_source_pin_count": report.get("parent_source_pin_count"),
        "parent_checkpoint_count": report.get("parent_checkpoint_count"),
        "parent_failed_artifact_count": report.get("parent_failed_artifact_count"),
        "parent_protected_file_count": report.get("parent_protected_file_count"),
        "parent_dirty_path_count": report.get("parent_dirty_path_count"),
        "actual_counts": report.get("actual_counts", {}),
    }
    return report_raw, report, {
        **maps,
        **retained,
        "parent_false_flags": list(false_flags),
        "parent_false_flag_values": dict(false_values),
        "parent_new_flags": dict(new_flags),
        "parent_counts": parent_counts,
        "parent_expected_counts": report.get("expected_counts", {}),
        "parent_expected_new_flags": report.get("expected_new_flags", {}),
        "parent_report": artifact(PARENT_REPORT, report_raw, hash_only=True),
    }


def manifest_value(manifest: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in manifest:
            return manifest[name]
    return None


def manifest_artifact(manifest: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = manifest_value(manifest, name, f"{name}_artifact")
    if not isinstance(value, Mapping):
        for container_name in ("artifacts", "code2"):
            container = manifest.get(container_name)
            if isinstance(container, Mapping) and isinstance(container.get(name), Mapping):
                value = container[name]
                break
    require(isinstance(value, Mapping), f"manifest {name} artifact missing")
    return value


def int_mapping(value: Any, label: str) -> dict[str, int]:
    require(isinstance(value, Mapping), f"{label} missing")
    result: dict[str, int] = {}
    for name, count in value.items():
        require(isinstance(name, str) and name, f"{label} name malformed")
        require(
            isinstance(count, int) and not isinstance(count, bool) and count >= 0,
            f"{label} count malformed: {name}",
        )
        result[name] = count
    return result


def expected_event_counts(manifest: Mapping[str, Any]) -> dict[str, int]:
    policy = manifest_value(manifest, "api_policy", "policy")
    value = manifest_value(manifest, "callback_event_counts", "expected_callback_event_counts")
    if value is None and isinstance(policy, Mapping):
        value = manifest_value(policy, "callback_event_counts", "event_counts", "expected_event_counts")
    return int_mapping(value, "API callback event counts")


def normalize_count_event_map(
    manifest: Mapping[str, Any],
    callback_counts: Mapping[str, int],
    expected_counts: Mapping[str, int],
) -> dict[str, tuple[str, ...]]:
    policy = manifest_value(manifest, "api_policy", "policy")
    value = manifest_value(manifest, "count_event_map")
    if value is None and isinstance(policy, Mapping):
        value = manifest_value(policy, "count_event_map")
    require(isinstance(value, Mapping), "manifest count_event_map missing")
    event_names = set(callback_counts)
    count_names = set(expected_counts)
    raw_map = dict(value)
    canonical: dict[str, tuple[str, ...]] = {}
    if set(raw_map).issubset(event_names):
        for event, targets in raw_map.items():
            if targets is None:
                canonical[str(event)] = ()
            elif isinstance(targets, str):
                canonical[str(event)] = (targets,)
            elif isinstance(targets, (list, tuple)):
                canonical[str(event)] = tuple(str(target) for target in targets)
            else:
                raise ValueError(f"count_event_map targets malformed: {event}")
    else:
        require(
            set(raw_map).issubset(count_names),
            "count_event_map keys are not callback events or expected counts",
        )
        for count_name, events in raw_map.items():
            if isinstance(events, str):
                events = (events,)
            require(
                isinstance(events, (list, tuple)),
                f"count_event_map events malformed: {count_name}",
            )
            for event in events:
                event_name = str(event)
                require(event_name in event_names, f"count_event_map event unknown: {event_name}")
                canonical[event_name] = (*canonical.get(event_name, ()), str(count_name))
    for event, targets in canonical.items():
        require(event in event_names, f"count_event_map event unknown: {event}")
        require(
            all(target in count_names for target in targets),
            f"count_event_map count unknown: {event}",
        )
    for event in event_names:
        canonical.setdefault(event, ())
    return canonical


def normalize_action_rules(
    manifest: Mapping[str, Any],
    expected_counts: Mapping[str, int],
    callback_counts: Mapping[str, int],
) -> list[dict[str, Any]]:
    # This is an optional generic manifest adapter only.  There is deliberately
    # no old binding-family fallback or retained/rewritten classification.
    policy = manifest_value(manifest, "api_policy", "policy")
    value = manifest_value(manifest, "count_action_rules")
    if value is None and isinstance(policy, Mapping):
        value = manifest_value(policy, "count_action_rules")
    if value is None:
        return []
    require(isinstance(value, (list, tuple)), "manifest count_action_rules must be a list")
    result: list[dict[str, Any]] = []
    for index, item in enumerate(value):
        require(isinstance(item, Mapping), f"count_action_rules[{index}] malformed")
        counter = item.get("counter")
        event = item.get("event")
        row_field = item.get("row_field")
        require(
            isinstance(counter, str) and counter in expected_counts,
            f"count_action_rules[{index}] counter unknown",
        )
        require(isinstance(event, str) and event in callback_counts, f"count_action_rules[{index}] event unknown")
        require(row_field is None or isinstance(row_field, str), f"count_action_rules[{index}] row_field malformed")
        require("equals" in item, f"count_action_rules[{index}] equals missing")
        result.append(
            {
                "counter": counter,
                "event": event,
                "row_field": row_field,
                "equals": item["equals"],
            }
        )
    return result


def load_manifest(path_value: str) -> tuple[Path, bytes, dict[str, Any]]:
    path = resolve_artifact(path_value)
    raw, manifest = read_json(path)
    require(manifest.get("schema") == MANIFEST_SCHEMA, "manifest schema changed")
    require(manifest.get("component") == COMPONENT, "manifest component changed")
    require(
        manifest.get("source_archive_sha256", SOURCE_ARCHIVE_SHA256) == SOURCE_ARCHIVE_SHA256,
        "source archive pin changed",
    )
    return path, raw, manifest


def authenticate_manifest(
    manifest_path: Path,
    manifest_raw: bytes,
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    expected_head = manifest_value(manifest, "expected_head", "current_head", "head", "expected_HEAD")
    require(expected_head == CURRENT_HEAD_COMMIT, "manifest expected HEAD changed")
    require(git_head() == CURRENT_HEAD_COMMIT, "repository HEAD changed")
    parent_pin = manifest_value(manifest, "parent_report", "parent_report_pin")
    require(isinstance(parent_pin, Mapping), "manifest parent report pin missing")
    parent_path, parent_raw, parent_sha = authenticate_spec(
        parent_pin,
        "manifest parent report",
        expected_path=PARENT_REPORT,
    )
    require(parent_sha == PARENT_REPORT_SHA256, "manifest parent report SHA changed")
    kernel_path, kernel_raw, kernel_sha = authenticate_spec(
        manifest_artifact(manifest, "kernel"),
        "manifest kernel",
        expected_path=KERNEL_PATH,
    )
    runner_path, runner_raw, runner_sha = authenticate_spec(
        manifest_artifact(manifest, "runner"),
        "manifest runner",
        expected_path=RUNNER_PATH,
    )
    kernel_ast = parse_ast(kernel_path, kernel_raw, "kernel")
    require("derive" in kernel_ast["top_level_definitions"], "frozen kernel has no derive API")
    runner_ast = parse_ast(runner_path, runner_raw, "runner")
    wrappers = manifest.get("source_wrappers")
    require(
        isinstance(wrappers, Mapping) and set(wrappers) == set(WRAPPER_NAMES),
        "manifest source wrapper names changed",
    )
    wrapper_data: dict[str, dict[str, Any]] = {}
    wrapper_artifacts: dict[str, dict[str, Any]] = {}
    for name in WRAPPER_NAMES:
        spec = wrappers[name]
        require(isinstance(spec, Mapping), f"wrapper {name} is not an artifact")
        path, raw, digest = authenticate_spec(
            spec,
            f"wrapper {name}",
            expected_path=WRAPPER_PATHS[name],
        )
        require(digest == WRAPPER_SHA256[name], f"wrapper {name} SHA changed")
        try:
            literal = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise ValueError(f"wrapper {name} is not JSON") from exc
        require(isinstance(literal, dict), f"wrapper {name} JSON object required")
        wrapper_data[name] = {"sha256": digest, "literal": literal}
        wrapper_artifacts[name] = artifact(path, raw)
    expected_counts = int_mapping(
        manifest_value(manifest, "expected_counts"),
        "manifest expected_counts",
    )
    expected_flags_value = manifest_value(
        manifest,
        "expected_new_flags",
        "newflags",
        "expected_newflags",
        "new_flags",
    )
    require(isinstance(expected_flags_value, Mapping), "manifest newflags missing")
    expected_flags = {str(name): flag for name, flag in expected_flags_value.items()}
    require(
        all(isinstance(flag, bool) for flag in expected_flags.values()),
        "manifest new flag value malformed",
    )
    policy = manifest_value(manifest, "api_policy", "policy")
    require(isinstance(policy, Mapping), "manifest API policy missing")
    require(
        policy.get("source_only") is True
        and policy.get("metadata_only") is False
        and policy.get("full_goal_complete") is False
        and policy.get("numeric_evaluation") is False,
        "source completion policy changed",
    )
    callback_counts = expected_event_counts(manifest)
    count_event_map = normalize_count_event_map(
        manifest,
        callback_counts,
        expected_counts,
    )
    action_rules = normalize_action_rules(manifest, expected_counts, callback_counts)
    nexttask = manifest_value(manifest, "nexttask", "next_task")
    if nexttask is None:
        nexttask = "source_native_O3_physical_spatial_jet_and_volume_norm_route"
    require(isinstance(nexttask, str) and nexttask, "manifest nexttask malformed")
    kernel_result_schema = manifest_value(manifest, "kernel_result_schema", "result_schema")
    if kernel_result_schema is None:
        kernel_result_schema = f"{COMPONENT}_result_v1"
    require(isinstance(kernel_result_schema, str) and kernel_result_schema, "manifest result schema malformed")
    return {
        "manifest": {
            "path": wire(manifest_path),
            "sha256": sha_bytes(manifest_raw),
            "bytes": len(manifest_raw),
        },
        "kernel": {**kernel_ast, "path": wire(kernel_path), "sha256": kernel_sha},
        "runner": {**runner_ast, "path": wire(runner_path), "sha256": runner_sha},
        "parent_report": artifact(parent_path, parent_raw, hash_only=True),
        "source_wrappers": wrapper_data,
        "source_wrapper_artifacts": wrapper_artifacts,
        "expected_counts": expected_counts,
        "expected_new_flags": expected_flags,
        "api_policy": dict(policy),
        "callback_event_counts": callback_counts,
        "count_event_map": count_event_map,
        "count_action_rules": action_rules,
        "nexttask": nexttask,
        "kernel_result_schema": kernel_result_schema,
        "proof_scope": manifest.get("proof_scope"),
        "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
    }


def merge_source_hashes(parent: Mapping[str, str], additions: Mapping[str, str]) -> dict[str, str]:
    result = dict(parent)
    for path, digest in additions.items():
        require(valid_sha(digest), f"new source pin SHA malformed: {path}")
        if path in result:
            require(result[path].lower() == digest.lower(), f"conflicting duplicate source pin: {path}")
        else:
            result[path] = digest.lower()
    return result


def no_base64_keys(value: Any) -> None:
    if isinstance(value, Mapping):
        for key, child in value.items():
            require("base64" not in str(key).lower(), "input contains a base64 field")
            no_base64_keys(child)
    elif isinstance(value, list):
        for child in value:
            no_base64_keys(child)


def build_input(
    parent: Mapping[str, Any],
    checked: Mapping[str, Any],
    context: Mapping[str, Any],
) -> dict[str, Any]:
    value = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_commit": CURRENT_HEAD_COMMIT,
        "parent_head": CURRENT_HEAD_COMMIT,
        "parent_report_internal_parent_commit": PARENT_REPORT_INTERNAL_PARENT_COMMIT,
        "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
        "source_wrappers": checked["source_wrappers"],
        "proof_scope": checked.get("proof_scope"),
        "policy": {
            "metadata_only": False,
            "source_only": True,
            "full_goal_complete": False,
            "numeric_evaluation": False,
        },
        "frozen_context": {
            "expected_parent_commit": CURRENT_HEAD_COMMIT,
            "parent_head": CURRENT_HEAD_COMMIT,
        },
        "parent_report": {
            "path": wire(PARENT_REPORT),
            "sha256": PARENT_REPORT_SHA256,
            "hash_only": True,
        },
        "saved_parent_counts": parent["parent_counts"],
        "execution_rules": context["execution_rules"],
        "accepted_mathematical_state": context["accepted_mathematical_state"],
        "physical_baseline": context["physical_baseline"],
        "goal_active": True,
        "full_goal_complete": False,
    }
    no_base64_keys(value)
    return value


def build_source(
    parent: Mapping[str, Any],
    checked: Mapping[str, Any],
    input_raw: bytes,
    context: Mapping[str, Any],
) -> dict[str, Any]:
    code2 = {"kernel": checked["kernel"], "runner": checked["runner"]}
    additions = {
        checked["kernel"]["path"]: checked["kernel"]["sha256"],
        checked["runner"]["path"]: checked["runner"]["sha256"],
        checked["manifest"]["path"]: checked["manifest"]["sha256"],
        wire(INPUT_PATH): sha_bytes(input_raw),
        checked["source_wrapper_artifacts"]["physical_source_capsule"]["path"]: checked[
            "source_wrapper_artifacts"
        ]["physical_source_capsule"]["sha256"],
    }
    source_union = merge_source_hashes(parent["source_pins"], additions)
    require(
        len(source_union) == EXPECTED_SOURCE_UNION_COUNT,
        "source union count is not frozen at 2330",
    )
    return {
        "schema": SOURCE_SCHEMA,
        "component": COMPONENT,
        "parent_component": "native_O3_finite_N_ast_flux_semantics",
        "parent_commit": CURRENT_HEAD_COMMIT,
        "parent_head": CURRENT_HEAD_COMMIT,
        "parent_report_internal_parent_commit": PARENT_REPORT_INTERNAL_PARENT_COMMIT,
        "current_head_commit": CURRENT_HEAD_COMMIT,
        "parent_report": parent["parent_report"],
        "parent_counts": parent["parent_counts"],
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_ARTIFACT_COUNT,
        "expected_failed_count_if_new_failure": FAILED_COUNT_IF_NEW_FAILURE,
        "parent_protected_file_count": PARENT_PROTECTED_FILE_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_PATH_COUNT,
        "parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "expected_published_checkpoint_count": PUBLISHED_CHECKPOINT_COUNT,
        "parent_expected_counts": parent["parent_expected_counts"],
        "parent_expected_new_flags": parent["parent_expected_new_flags"],
        "parent_new_flags": parent["parent_new_flags"],
        "source_pins": dict(parent["source_pins"]),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "protected_file_hashes": dict(parent["protected_file_hashes"]),
        "legacy_dirty_path_hashes": dict(parent["legacy_dirty_path_hashes"]),
        "preserved_failed_artifact_hashes_union": dict(parent["preserved_failed_artifact_hashes_union"]),
        "parent_false_flags": list(parent["parent_false_flags"]),
        "parent_false_flag_values": dict(parent["parent_false_flag_values"]),
        "preserved_error_cap_metadata": parent["preserved_error_cap_metadata"],
        "retained_construction_data": parent["retained_construction_data"],
        "retained_integral_term_inventory": parent["retained_integral_term_inventory"],
        "retained_cap_provenance": parent["retained_cap_provenance"],
        "code2": code2,
        "manifest": checked["manifest"],
        "input_spec": {
            "path": wire(INPUT_PATH),
            "sha256": sha_bytes(input_raw),
            "bytes": len(input_raw),
        },
        "input": json.loads(input_raw.decode("utf-8")),
        "source_wrappers": {
            name: {
                **checked["source_wrapper_artifacts"][name],
                "literal": checked["source_wrappers"][name]["literal"],
            }
            for name in checked["source_wrapper_artifacts"]
        },
        "expected_counts": dict(checked["expected_counts"]),
        "expected_new_flags": dict(checked["expected_new_flags"]),
        "callback_event_counts": dict(checked["callback_event_counts"]),
        "count_event_map": {
            name: list(targets) for name, targets in checked["count_event_map"].items()
        },
        "count_action_rules": list(checked["count_action_rules"]),
        "callback_event_count": sum(checked["callback_event_counts"].values()),
        "journal_event_count": sum(checked["callback_event_counts"].values()) + 2,
        "nexttask": checked["nexttask"],
        "kernel_result_schema": checked["kernel_result_schema"],
        "proof_scope": checked.get("proof_scope"),
        "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
        "execution_rules": context["execution_rules"],
        "accepted_mathematical_state": context["accepted_mathematical_state"],
        "physical_baseline": context["physical_baseline"],
        "goal_active": True,
        "metadata_only": False,
        "source_only": True,
        "full_goal_complete": False,
        "numeric_evaluation": False,
        "runtime_model_metadata": None,
        "no_project_imports_before_run": True,
    }


def prepare(manifest_value_path: str) -> dict[str, Any]:
    require(not INPUT_PATH.exists(), "input artifact already exists")
    require(not SOURCE_PATH.exists(), "source checkpoint already exists")
    require(not BRIEF_PATH.exists(), "pending brief already exists")
    require(
        not ENTERED_PATH.exists() and not JOURNAL_PATH.exists() and not RESULT_PATH.exists(),
        "attempt already entered",
    )
    manifest_path, manifest_raw, manifest = load_manifest(manifest_value_path)
    checked = authenticate_manifest(manifest_path, manifest_raw, manifest)
    _, _, parent = load_parent()
    current_raw, current_brief, context = load_active_brief()
    input_value = build_input(parent, checked, context)
    input_raw = write_new(INPUT_PATH, input_value)
    preimage_path, preimage_raw = brief_preimage(current_raw)
    source = build_source(parent, checked, input_raw, context)
    source["main_brief_preimage"] = {
        "path": wire(preimage_path),
        "sha256": sha_bytes(preimage_raw),
        "bytes": len(preimage_raw),
        "raw_bytes_exact": True,
        "raw_readonly": True,
    }
    source_raw = write_new(SOURCE_PATH, source)
    pending = {
        "schema": f"{COMPONENT}_brief_v1",
        "component": COMPONENT,
        "phase": "prepared",
        "current_component": COMPONENT,
        "source_checkpoint": {"path": wire(SOURCE_PATH), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire(INPUT_PATH), "sha256": sha_bytes(input_raw)},
        "code2": source["code2"],
        "manifest": source["manifest"],
        "source_wrappers": source["source_wrappers"],
        "expected_counts": source["expected_counts"],
        "expected_new_flags": source["expected_new_flags"],
        "callback_event_counts": source["callback_event_counts"],
        "count_event_map": source["count_event_map"],
        "count_action_rules": source["count_action_rules"],
        "proof_scope": source["proof_scope"],
        "main_brief_preimage": source["main_brief_preimage"],
        "journal": {"path": wire(JOURNAL_PATH), "schema": JOURNAL_SCHEMA},
        "metadata_only": False,
        "source_only": True,
        "full_goal_complete": False,
        "runtime_model_metadata": None,
    }
    pending_raw = write_new(BRIEF_PATH, pending)
    main_brief = update_brief(
        current_raw,
        current_brief,
        {
            "phase": "prepared",
            "current_component": COMPONENT,
            "in_progress": {
                COMPONENT: {
                    "source": wire(SOURCE_PATH),
                    "source_sha256": sha_bytes(source_raw),
                    "input": wire(INPUT_PATH),
                    "input_sha256": sha_bytes(input_raw),
                    "code2": source["code2"],
                    "manifest": source["manifest"],
                    "metadata_only": False,
                    "source_only": True,
                    "full_goal_complete": False,
                }
            },
            "goal_active": True,
            "full_goal_complete": False,
            "next_task": checked["nexttask"],
            "next_executable_task": checked["nexttask"],
        },
    )
    return {
        "status": "prepared",
        "component": COMPONENT,
        "source_checkpoint": {"path": wire(SOURCE_PATH), "sha256": sha_bytes(source_raw)},
        "input": {"path": wire(INPUT_PATH), "sha256": sha_bytes(input_raw)},
        "pending_brief": {"path": wire(BRIEF_PATH), "sha256": sha_bytes(pending_raw)},
        "main_brief": main_brief,
        "source_union_count": source["source_union_count"],
        "expected_counts": source["expected_counts"],
        "expected_new_flags": source["expected_new_flags"],
        "callback_event_counts": source["callback_event_counts"],
        "count_action_rules": source["count_action_rules"],
        "callback_event_count": source["callback_event_count"],
        "journal_event_count": source["journal_event_count"],
        "metadata_only": False,
        "source_only": True,
        "runtime_model_metadata": None,
        "no_run_performed": True,
    }


def load_source_for_run() -> tuple[bytes, dict[str, Any]]:
    source_raw, source = read_json(SOURCE_PATH)
    require(
        source.get("schema") == SOURCE_SCHEMA and source.get("component") == COMPONENT,
        "source checkpoint schema changed",
    )
    require(
        source.get("parent_commit") == CURRENT_HEAD_COMMIT
        and source.get("parent_head") == CURRENT_HEAD_COMMIT
        and source.get("current_head_commit") == CURRENT_HEAD_COMMIT,
        "source parent HEAD changed",
    )
    require(
        source.get("parent_report_internal_parent_commit") == PARENT_REPORT_INTERNAL_PARENT_COMMIT,
        "source report lineage changed",
    )
    require(source.get("source_archive_sha256") == SOURCE_ARCHIVE_SHA256, "source archive pin changed")
    preimage = source.get("main_brief_preimage")
    require(isinstance(preimage, Mapping), "main brief preimage pointer missing")
    read_bytes(resolve_artifact(preimage["path"], preimage["sha256"]), preimage["sha256"])
    input_path = resolve_artifact(source["input_spec"]["path"], source["input_spec"]["sha256"])
    require(input_path.resolve() == INPUT_PATH.resolve(), "frozen input path changed")
    input_raw, input_value = read_json(input_path, source["input_spec"]["sha256"])
    require(input_value == source.get("input"), "frozen input literal changed")
    manifest_path = resolve_artifact(source["manifest"]["path"], source["manifest"]["sha256"])
    manifest_raw, manifest = read_json(manifest_path, source["manifest"]["sha256"])
    checked = authenticate_manifest(manifest_path, manifest_raw, manifest)
    require(checked["kernel"]["sha256"] == source["code2"]["kernel"]["sha256"], "frozen kernel changed")
    require(checked["runner"]["sha256"] == source["code2"]["runner"]["sha256"], "frozen runner changed")
    require(checked["expected_counts"] == source["expected_counts"], "frozen expected counts changed")
    require(checked["expected_new_flags"] == source["expected_new_flags"], "frozen expected flags changed")
    require(
        checked["callback_event_counts"] == source["callback_event_counts"],
        "frozen callback counts changed",
    )
    require(
        {name: list(targets) for name, targets in checked["count_event_map"].items()}
        == source["count_event_map"],
        "frozen count_event_map changed",
    )
    require(
        checked["count_action_rules"] == source.get("count_action_rules", []),
        "frozen count_action_rules changed",
    )
    require(set(source.get("source_wrappers", {})) == set(WRAPPER_NAMES), "frozen wrapper names changed")
    for name, wrapper in source["source_wrappers"].items():
        require(name in WRAPPER_PATHS, f"unknown frozen wrapper: {name}")
        path, raw, digest = authenticate_spec(
            wrapper,
            f"frozen wrapper {name}",
            expected_path=WRAPPER_PATHS[name],
        )
        require(digest == wrapper["sha256"], f"frozen wrapper SHA changed: {name}")
        require(digest == WRAPPER_SHA256[name], f"frozen wrapper SHA changed: {name}")
        require(
            json.loads(raw.decode("utf-8")) == wrapper.get("literal"),
            f"frozen wrapper literal changed: {name}",
        )
    for key in (
        "source_pins",
        "protected_file_hashes",
        "legacy_dirty_path_hashes",
        "preserved_failed_artifact_hashes_union",
    ):
        authenticate_hash_map(source[key], f"frozen {key}")
    require(source.get("source_union_count") == EXPECTED_SOURCE_UNION_COUNT, "source union count changed")
    require(
        len(source.get("source_union_hashes", {})) == EXPECTED_SOURCE_UNION_COUNT,
        "source union ledger changed",
    )
    for path_value, digest in source["source_union_hashes"].items():
        read_bytes(resolve_artifact(path_value, digest), digest)
    return source_raw, source


def append_journal(row: Mapping[str, Any]) -> None:
    raw = encode(dict(row))
    with JOURNAL_PATH.open("ab") as handle:
        handle.write(raw)
        handle.flush()
        os.fsync(handle.fileno())


def journal_rows() -> list[dict[str, Any]]:
    if not JOURNAL_PATH.exists():
        return []
    rows: list[dict[str, Any]] = []
    for line in JOURNAL_PATH.read_bytes().splitlines():
        if line:
            value = json.loads(line.decode("utf-8"))
            require(isinstance(value, dict), "journal rows must be objects")
            rows.append(value)
    return rows


def import_kernel(path: Path) -> Any:
    require(path.resolve() != Path(__file__).resolve(), "kernel path aliases metadata runner")
    spec = importlib.util.spec_from_file_location(
        "_native_O3_total_radius_velocity_pressure_binding_entry",
        path,
    )
    require(spec is not None and spec.loader is not None, "kernel import spec unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    require(sys.modules.get(spec.name) is module, "kernel was not registered in sys.modules")
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "frozen kernel has no derive API")
    return module


def callback_rows_callback(source: Mapping[str, Any]) -> tuple[Any, list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    allowed = set(source["callback_event_counts"])

    def emit(*args: Any) -> None:
        require(
            len(args) == 2
            and isinstance(args[0], str)
            and isinstance(args[1], Mapping),
            "emit expects event name plus mapping",
        )
        row = {"event": args[0], "row": dict(args[1])}
        require(row["event"] in allowed, f"unexpected callback event: {row['event']}")
        append_journal(row)
        rows.append(row)

    return emit, rows


def _field_value(value: Any, path: str | None) -> Any:
    if path is None or path == "":
        return value
    if path.startswith("row."):
        path = path[4:]
    current = value
    for part in path.split("."):
        if not isinstance(current, Mapping):
            return None
        current = current.get(part)
    return current


def count_event_prefix(
    source: Mapping[str, Any],
    rows: list[Mapping[str, Any]],
    *,
    failed: bool = False,
) -> dict[str, int]:
    expected = source.get("expected_counts")
    event_map = source.get("count_event_map")
    require(isinstance(expected, Mapping) and isinstance(event_map, Mapping), "frozen count metadata missing")
    actual = {str(name): 0 for name in expected}
    for record in rows:
        event = record.get("event")
        targets = event_map.get(event, ())
        if isinstance(targets, str):
            targets = (targets,)
        if isinstance(targets, (list, tuple)):
            for target in targets:
                if str(target) in actual:
                    actual[str(target)] += 1
    action_rules = source.get("count_action_rules") or []
    for rule in action_rules:
        counter = rule["counter"]
        event = rule["event"]
        row_field = rule.get("row_field")
        equals = rule.get("equals")
        actual[counter] = sum(
            1
            for record in rows
            if record.get("event") == event
            and _field_value(record.get("row"), row_field) == equals
        )
    if failed and "failed_family_calls" in actual:
        actual["failed_family_calls"] = int(
            any(record.get("event") == "family_entered" for record in rows)
        )
    if failed and "new_mathematical_proofs" in actual:
        actual["new_mathematical_proofs"] = 0
    return actual


def false_flags(source: Mapping[str, Any]) -> dict[str, bool]:
    expected = source.get("expected_new_flags")
    require(isinstance(expected, Mapping), "frozen expected flags missing")
    return {str(name): False for name in expected}


def kernel_failure_result(
    source: Mapping[str, Any],
    error: BaseException,
    rows: list[dict[str, Any]],
    *,
    imported: bool,
    derive_entry_calls: int,
    derive_calls: int,
) -> dict[str, Any]:
    return {
        "schema": source.get("kernel_result_schema", f"{COMPONENT}_result_v1"),
        "component": COMPONENT,
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "full_goal_complete": False,
        "runtime_model_metadata": None,
        "failure": {
            "kind": "kernel_import_or_derive_exception",
            "error_type": type(error).__name__,
            "error": repr(error),
            "traceback": traceback.format_exc(),
            "kernel_imported": imported,
            "derive_entry_calls": derive_entry_calls,
            "derive_calls": derive_calls,
            "callback_prefix": list(rows),
        },
        "proof_count": 0,
        "new_flags": false_flags(source),
        "actual_counts": count_event_prefix(source, rows, failed=True),
        "callback_events": list(rows),
    }


def callback_summary(source: Mapping[str, Any], rows: list[Mapping[str, Any]]) -> dict[str, Any]:
    expected = source.get("callback_event_counts")
    require(isinstance(expected, Mapping), "frozen callback counts missing")
    counts = {str(name): 0 for name in expected}
    for row in rows:
        name = row.get("event")
        require(name in counts, f"unexpected callback event: {name}")
        counts[str(name)] += 1
    return {
        "callback_event_count": len(rows),
        "callback_event_counts": counts,
        "journal_event_count": len(journal_rows()),
        "journal_event_count_expected": sum(expected.values()) + 2,
    }


def callback_counts_prefix(source: Mapping[str, Any], rows: list[Mapping[str, Any]]) -> dict[str, int]:
    expected = source.get("callback_event_counts")
    require(isinstance(expected, Mapping), "frozen callback counts missing")
    counts = {str(name): 0 for name in expected}
    for row in rows:
        name = row.get("event")
        require(name in counts, f"unexpected callback event: {name}")
        counts[str(name)] += 1
    return counts


def wrap_result(
    source: Mapping[str, Any],
    kernel_result: Mapping[str, Any],
    *,
    imported: bool,
    derive_entry_calls: int,
    derive_calls: int,
    rows: list[dict[str, Any]],
    manager_failure: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    nested = dict(kernel_result)
    kernel_failed = nested.get("passed") is not True
    failure_record: Mapping[str, Any] | None = manager_failure
    if failure_record is None and kernel_failed:
        failure_record = {
            "kind": "kernel_reported_failure",
            "kernel_failure": nested.get("failure"),
            "kernel_prefix_preserved": True,
        }
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "failed" if kernel_failed else "passed",
        "passed": not kernel_failed,
        "metadata_only": False,
        "source_only": True,
        "full_goal_complete": False,
        "numeric_evaluation": False,
        "result": nested,
        "actual_counts": count_event_prefix(source, rows, failed=kernel_failed),
        "new_flags": false_flags(source) if kernel_failed else nested.get("new_flags", {}),
        "proof_count": 0 if kernel_failed else nested.get("proof_count", 0),
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_imported": imported,
        "callback_events": list(rows),
        "callback_event_count": len(rows),
        "callback_event_counts": callback_counts_prefix(source, rows),
        "failure": dict(failure_record) if failure_record is not None else None,
        "runtime_model_metadata": None,
        "source_checkpoint_sha256": sha_bytes(SOURCE_PATH.read_bytes()),
    }


def validate_saved(
    source: Mapping[str, Any],
    saved: Mapping[str, Any],
    rows: list[dict[str, Any]],
    *,
    imported: bool,
    derive_entry_calls: int,
    derive_calls: int,
) -> dict[str, Any]:
    violations: list[str] = []
    if saved.get("schema") != RESULT_SCHEMA:
        violations.append("result schema")
    if saved.get("component") != COMPONENT:
        violations.append("result component")
    if saved.get("metadata_only") is not False or saved.get("source_only") is not True:
        violations.append("scope")
    if saved.get("runtime_model_metadata") is not None:
        violations.append("runtime model metadata")
    try:
        summary = callback_summary(source, rows)
    except (TypeError, ValueError) as exc:
        summary = {
            "callback_event_counts": {},
            "journal_event_count": len(journal_rows()),
            "journal_event_count_expected": sum(source.get("callback_event_counts", {}).values()) + 2,
        }
        violations.append(f"callback journal: {exc}")
    if saved.get("callback_event_count") != len(rows):
        violations.append("callback event count")
    if saved.get("callback_event_counts") != summary["callback_event_counts"]:
        violations.append("callback event counts prefix")
    failed = saved.get("passed") is not True
    if saved.get("actual_counts") != count_event_prefix(source, rows, failed=failed):
        violations.append("actual prefix counts")
    nested = saved.get("result")
    if not isinstance(nested, Mapping):
        violations.append("nested result")
    if failed:
        if saved.get("new_flags") != false_flags(source) or saved.get("proof_count") != 0:
            violations.append("failed flags or proof count")
        if not isinstance(nested, Mapping) or nested.get("passed") is not False:
            violations.append("nested failure result")
        if not isinstance(nested, Mapping) or not isinstance(nested.get("failure"), Mapping):
            violations.append("nested kernel failure")
        if not isinstance(saved.get("failure"), Mapping):
            violations.append("manager or kernel failure")
        expected_manager_events = 2 if imported else 1
        if len(journal_rows()) != expected_manager_events + len(rows):
            violations.append("failed journal prefix length")
    else:
        if saved.get("actual_counts") != source.get("expected_counts"):
            violations.append("actual count values")
        if saved.get("new_flags") != source.get("expected_new_flags"):
            violations.append("new flag values")
        if (
            saved.get("derive_entry_calls") != 1
            or saved.get("derive_calls") != 1
            or saved.get("kernel_imported") is not True
        ):
            violations.append("entry counts")
        if summary["callback_event_counts"] != source.get("callback_event_counts"):
            violations.append("callback event counts")
        if summary["journal_event_count"] != summary["journal_event_count_expected"]:
            violations.append("journal event count")
        if nested.get("passed") is not True:
            violations.append("nested success result")
    return {
        "schema": f"{COMPONENT}_validation_v1",
        "component": COMPONENT,
        "passed": not violations and not failed,
        "attempt_failed": failed,
        "actual_counts": saved.get("actual_counts"),
        "new_flags": saved.get("new_flags"),
        "callback_event_count": len(rows),
        "journal_event_count": len(journal_rows()),
        "violations": violations,
    }


def run() -> dict[str, Any]:
    source_raw, source = load_source_for_run()
    require(not ENTERED_PATH.exists(), "entered checkpoint already exists; reentry forbidden")
    require(
        not JOURNAL_PATH.exists() and not RESULT_PATH.exists(),
        "run checkpoints already exist; reentry forbidden",
    )
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "kernel": source["code2"]["kernel"],
        "runner": source["code2"]["runner"],
        "import_attempts": 1,
        "derive_entry_calls": 0,
        "derive_calls": 0,
        "sys_modules_registration_before_exec": True,
        "runtime_model_metadata": None,
        "entered_at_utc": now(),
        "do_not_reenter": True,
    }
    entered_raw = write_new(ENTERED_PATH, entered)
    write_new(JOURNAL_PATH, b"")
    append_journal(
        {
            "schema": JOURNAL_SCHEMA,
            "event": "manager_entry_before_import",
            "source_checkpoint_sha256": sha_bytes(source_raw),
            "entered_checkpoint_sha256": sha_bytes(entered_raw),
            "import_attempts": 1,
            "derive_entry_calls": 0,
            "derive_calls": 0,
        }
    )
    emit, rows = callback_rows_callback(source)
    imported = False
    derive_entry_calls = 0
    derive_calls = 0
    kernel_result: Mapping[str, Any] | None = None
    manager_failure: Mapping[str, Any] | None = None
    try:
        module = import_kernel(KERNEL_PATH)
        imported = True
        derive_entry_calls = 1
        append_journal(
            {
                "schema": JOURNAL_SCHEMA,
                "event": "manager_derive_entry",
                "import_attempts": 1,
                "derive_entry_calls": 1,
                "derive_calls": 0,
                "sys_modules_registered": True,
            }
        )
        derive_calls = 1
        returned = module.derive(source["input"], emit)
        require(isinstance(returned, Mapping), "kernel derive returned a non-mapping")
        kernel_result = dict(returned)
    except BaseException as exc:
        kernel_result = kernel_failure_result(
            source,
            exc,
            rows,
            imported=imported,
            derive_entry_calls=derive_entry_calls,
            derive_calls=derive_calls,
        )
        manager_failure = {
            "kind": "manager_import_or_derive_failure",
            "error_type": type(exc).__name__,
            "error": repr(exc),
            "traceback": traceback.format_exc(),
            "kernel_imported": imported,
            "derive_entry_calls": derive_entry_calls,
            "derive_calls": derive_calls,
            "callback_prefix": list(rows),
        }
    wrapped = wrap_result(
        source,
        kernel_result or {},
        imported=imported,
        derive_entry_calls=derive_entry_calls,
        derive_calls=derive_calls,
        rows=rows,
        manager_failure=manager_failure,
    )
    # The complete result is durable before any validation or review.
    result_raw = write_new(RESULT_PATH, wrapped)
    saved_raw, saved = read_json(RESULT_PATH, sha_bytes(result_raw))
    validation = validate_saved(
        source,
        saved,
        rows,
        imported=imported,
        derive_entry_calls=derive_entry_calls,
        derive_calls=derive_calls,
    )
    return {
        "status": "run_complete",
        "component": COMPONENT,
        "attempt_status": saved.get("attempt_status"),
        "passed": bool(validation["passed"]),
        "result": {
            "path": wire(RESULT_PATH),
            "sha256": sha_bytes(saved_raw),
            "schema": saved.get("schema"),
        },
        "actual_counts": saved.get("actual_counts"),
        "new_flags": saved.get("new_flags"),
        "proof_count": saved.get("proof_count"),
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_imported": imported,
        "callback_event_count": len(rows),
        "journal_event_count": len(journal_rows()),
        "validation": validation,
        "failure": saved.get("failure"),
        "runtime_model_metadata": None,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--manifest", required=True)
    subparsers.add_parser("run")
    args = parser.parse_args(argv)
    try:
        receipt = prepare(args.manifest) if args.command == "prepare" else run()
    except Exception as exc:
        print(json.dumps({"status": "error", "error": f"{type(exc).__name__}: {exc}"}, sort_keys=True))
        return 2
    print(json.dumps(receipt, sort_keys=True))
    return 0 if receipt.get("passed", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
