"""Durable, one-shot runner for the native mixed radial operator leaf.

The prepare action authenticates only raw artifacts, JSON literals, and Python
ASTs.  It never imports project code.  The run action is intentionally the
only place where the frozen new kernel is imported, and it records durable
entry and callback bytes before doing any metadata validation.
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

COMPONENT = "native_O3_mixed_radial_operator"
INPUT_SCHEMA = f"{COMPONENT}_input_v1"
MANIFEST_SCHEMA = f"{COMPONENT}_manifest_v1"
SOURCE_SCHEMA = f"{COMPONENT}_source_v1"
ENTERED_SCHEMA = f"{COMPONENT}_entered_v1"
RESULT_SCHEMA = f"{COMPONENT}_result_checkpoint_v1"
JOURNAL_SCHEMA = f"{COMPONENT}_journal_v1"

PARENT_COMMIT = "72dccb761638c0e14001bf2e4fc76c3f6dff7762"
PARENT_REPORT = REPO / (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_parameter_Z_independence.json"
)
PARENT_REPORT_SHA256 = "0587c2c5ec692f6090f02511b56adbb2a79238836081c802d1211eb3348b9d3c"
PARENT_RESULT = REPO / "tmp/native_parameter_Z_independence/result.json"
PARENT_RESULT_SHA256 = "e0ef019b51898a83d01b58940face6023629ec3b4ef245bca903d4d20a218eb7"
PARENT_SOURCE_PIN_COUNT = 2287
PARENT_CHECKPOINT_COUNT = 272
PARENT_FAILED_COUNT = 183
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
SOURCE_ARCHIVE_SHA256 = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"

MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
KERNEL_RELATIVE = Path("src/openai_ns_reconstruction/native_O3_mixed_radial_operator.py")
RUNNER_RELATIVE = Path("scripts/run_native_O3_mixed_radial_operator.py")
KERNEL_PATH = REPO / KERNEL_RELATIVE
RUNNER_PATH = REPO / RUNNER_RELATIVE
INPUT_PATH = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF_PATH = WORK / f"{COMPONENT}_brief_2026_10_27.json"
CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_PATH = CHECKPOINT_DIR / "source.json"
ENTERED_PATH = CHECKPOINT_DIR / "entered.json"
JOURNAL_PATH = CHECKPOINT_DIR / "journal.jsonl"
RESULT_PATH = CHECKPOINT_DIR / "result.json"

EXPECTED_WRAPPERS = {
    "saved_Z_result": {
        "path": PARENT_RESULT,
        "sha256": PARENT_RESULT_SHA256,
    },
    "background_literals": {
        "path": WORK / "native_O3_background_mixed_source_literals_2026_10_27.json",
        "sha256": "2748333d5dc25c5a34b23235555dc7e285f52c05e5b1ae2995e50aa2281bca14",
    },
    "modified_literals": {
        "path": WORK / "native_O3_modified_mixed_source_literals_2026_10_27.json",
        "sha256": "d1834f1ed518b8a4525edd110fc3fc9d551db6504090d7668db3dad027cc3906",
    },
}

# These are the seven global rules.  They are compared byte-for-byte as JSON
# strings with the active resume brief and copied back unchanged on update.
EXECUTION_RULES = [
    "Never rerun any entered ancestor mathematical/helper family, including failed leaves, old builders/providers/normalizers/checkers, source-value reducers, quadrature or recurrence. Static source/AST/literal JSON/hash inspection is allowed.",
    "Before a NEW mathematical/source-expression consumer: freeze code and input, persist source checkpoint and pending brief, persist entry before invocation, invoke once, save complete result before review. Preserve a failure and repair in a separate frozen leaf; never reenter.",
    "Authenticate raw gzip SHA BEFORE decoding. Never decode docs/continuation/RESUME_FIRST_SWITCH_2026_10_17.json or the huge Rm_loop_density_identity.json; raw SHA only.",
    "Never restart stopped session 17273. Retain nonzero incoming histories and ordinary Z1. Bounds, midpoints and interval overlap are not function values; remove old numerical errors only with a same-source packet theorem.",
    "Use current tree and GitHub; preserve all dirty and untracked files. Commit/push only exact accepted scoped files. Goal remains ACTIVE / INCOMPLETE.",
    "Children: luna_max_worker, gpt-5.6-luna/max fixed role, fork_turns:none, complete brief. Runtime model metadata unavailable means null; never claim it was verified.",
    "A worker function that checks mathematical source identities counts as a source consumer even when it consumes static receipts or is called a pure self-check. Workers return unexecuted code; root alone freezes code and inputs, records entry, invokes once and saves the complete result before review. AST parsing and metadata-only audits remain permitted. Do not execute or reconstruct the rejected memory draft result.",
]

DEFAULT_EVENT_COUNTS = {
    "compiled_row": 60,
    "unavailable_source": 6,
    "api_refusal": 4,
    "family_entered": 1,
    "family_completed": 1,
}


def require(condition: Any, message: str) -> None:
    if not condition:
        raise ValueError(message)


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def valid_sha(value: Any) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(character in "0123456789abcdefABCDEF" for character in value)
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
        REPO / "experiments/root_st073" / path,
        WORK / path,
        WORK / "experiments/root_st073" / path,
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
    value = {
        "path": wire(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "raw_bytes_exact": True,
        "raw_readonly": True,
    }
    if hash_only:
        value["hash_only"] = True
        value["literal_embedded"] = False
    return value


def authenticate_spec(spec: Mapping[str, Any], label: str, *, expected_path: Path | None = None) -> tuple[Path, bytes, str]:
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
    return {"path": wire(path), "sha256": sha_bytes(raw), "ast_parse": "ok", "node_count": sum(1 for _ in ast.walk(tree))}


def git_head() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def load_active_brief() -> tuple[bytes, dict[str, Any]]:
    raw, brief = read_json(MAIN_BRIEF)
    require(brief.get("execution_rules") == EXECUTION_RULES, "global execution rules changed")
    require(brief.get("goal_active") is True, "main brief is not active")
    require(brief.get("full_goal_complete") is False, "main brief is already complete")
    return raw, brief


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
    return {"path": wire(MAIN_BRIEF), "sha256": sha_bytes(new_raw), "bytes": len(new_raw), "preimage_sha256": sha_bytes(raw)}


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
        path = resolve_artifact(path_value, digest)
        read_bytes(path, digest)
    return mapping


def load_parent() -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    report_raw, report = read_json(PARENT_REPORT, PARENT_REPORT_SHA256)
    require(report.get("component") == "native_parameter_Z_independence", "parent component changed")
    require(report.get("passed") is True, "parent report is not passed")
    require(report.get("source_only") is True and report.get("full_goal_complete") is False, "parent scope changed")
    maps = {
        "source_pins": authenticate_hash_map(report.get("source_pins"), "parent source pins"),
        "protected_file_hashes": authenticate_hash_map(report.get("protected_file_hashes"), "parent protected files"),
        "legacy_dirty_path_hashes": authenticate_hash_map(report.get("legacy_dirty_path_hashes"), "parent dirty files"),
        "preserved_failed_artifact_hashes_union": authenticate_hash_map(report.get("preserved_failed_artifact_hashes_union"), "parent failed files"),
    }
    require(len(maps["source_pins"]) == PARENT_SOURCE_PIN_COUNT, "parent source-pin count changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT, "parent checkpoint count changed")
    require(len(maps["preserved_failed_artifact_hashes_union"]) == PARENT_FAILED_COUNT, "parent failed count changed")
    require(len(maps["protected_file_hashes"]) == PARENT_PROTECTED_COUNT, "parent protected count changed")
    require(len(maps["legacy_dirty_path_hashes"]) == PARENT_DIRTY_COUNT, "parent dirty count changed")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT, "parent source_pin_count field changed")
    require(report.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT, "parent preserved_failed_artifact_count field changed")
    require(report.get("expected_parent_false_flag_count") == PARENT_FALSE_FLAG_COUNT, "parent false-flag count field changed")
    false_flags = report.get("parent_false_flags")
    false_values = report.get("parent_false_flag_values")
    require(isinstance(false_flags, list) and len(false_flags) == PARENT_FALSE_FLAG_COUNT, "parent false-flag count changed")
    require(false_values == {name: False for name in false_flags}, "parent false-flag values changed")
    retained = {}
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
        "checkpoints": int(report.get("checkpoint_count")),
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
    parent_pointer = artifact(PARENT_REPORT, report_raw, hash_only=True)
    return report_raw, report, {
        **maps,
        **retained,
        "parent_false_flags": list(false_flags),
        "parent_false_flag_values": dict(false_values),
        "parent_counts": parent_counts,
        "parent_expected_counts": report.get("expected_counts", {}),
        "parent_expected_new_flags": report.get("expected_new_flags", {}),
        "parent_report": parent_pointer,
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


def expected_event_counts(manifest: Mapping[str, Any]) -> dict[str, int]:
    policy = manifest_value(manifest, "api_policy", "policy")
    value = policy.get("callback_event_counts") if isinstance(policy, Mapping) else None
    if value is None and isinstance(policy, Mapping):
        value = policy.get(
            "event_counts",
            policy.get("expected_event_counts", policy.get("expected_callback_event_counts")),
        )
    if value is None:
        value = DEFAULT_EVENT_COUNTS
    require(isinstance(value, Mapping), "API callback event counts missing")
    result = {str(name): int(count) for name, count in value.items()}
    require(all(count >= 0 for count in result.values()), "negative callback event count")
    require(sum(result.values()) == 72, "kernel callback event total must be 72")
    return result


def load_manifest(path_value: str) -> tuple[Path, bytes, dict[str, Any]]:
    path = resolve_artifact(path_value)
    raw, manifest = read_json(path)
    require(manifest.get("schema") == MANIFEST_SCHEMA, "manifest schema changed")
    require(manifest.get("component") == COMPONENT, "manifest component changed")
    require(manifest.get("source_archive_sha256", SOURCE_ARCHIVE_SHA256) == SOURCE_ARCHIVE_SHA256, "source archive pin changed")
    return path, raw, manifest


def authenticate_manifest(manifest_path: Path, manifest_raw: bytes, manifest: Mapping[str, Any]) -> dict[str, Any]:
    expected_head = manifest_value(manifest, "expected_head", "head", "expected_HEAD")
    require(expected_head == PARENT_COMMIT, "manifest expected HEAD changed")
    require(git_head() == PARENT_COMMIT, "repository HEAD changed")
    parent_pin = manifest_value(manifest, "parent_report", "parent_report_pin")
    require(isinstance(parent_pin, Mapping), "manifest parent report pin missing")
    parent_path, parent_raw, parent_sha = authenticate_spec(parent_pin, "manifest parent report", expected_path=PARENT_REPORT)
    require(parent_sha == PARENT_REPORT_SHA256, "manifest parent report SHA changed")
    kernel_path, kernel_raw, kernel_sha = authenticate_spec(
        manifest_artifact(manifest, "kernel"), "manifest kernel", expected_path=KERNEL_PATH
    )
    runner_path, runner_raw, runner_sha = authenticate_spec(
        manifest_artifact(manifest, "runner"), "manifest runner", expected_path=RUNNER_PATH
    )
    kernel_ast = parse_ast(kernel_path, kernel_raw, "kernel")
    runner_ast = parse_ast(runner_path, runner_raw, "runner")
    wrappers = manifest.get("source_wrappers")
    require(isinstance(wrappers, Mapping), "manifest source_wrappers missing")
    require(set(wrappers) == set(EXPECTED_WRAPPERS), "manifest wrapper names changed")
    wrapper_data: dict[str, dict[str, Any]] = {}
    wrapper_artifacts: dict[str, dict[str, Any]] = {}
    for name, expected in EXPECTED_WRAPPERS.items():
        spec = wrappers[name]
        require(isinstance(spec, Mapping), f"wrapper {name} is not an artifact")
        expected_digest = expected["sha256"]
        require(spec.get("sha256") == expected_digest, f"wrapper {name} SHA contract changed")
        path, raw, digest = authenticate_spec(spec, f"wrapper {name}", expected_path=expected["path"])
        literal = json.loads(raw.decode("utf-8"))
        wrapper_data[name] = {"sha256": digest, "literal": literal}
        wrapper_artifacts[name] = artifact(path, raw)
    expected_counts = manifest_value(manifest, "expected_counts")
    expected_flags = manifest_value(manifest, "expected_new_flags", "newflags", "expected_newflags", "new_flags")
    require(isinstance(expected_counts, Mapping), "manifest expected_counts missing")
    require(isinstance(expected_flags, Mapping), "manifest newflags missing")
    policy = manifest_value(manifest, "api_policy", "policy")
    require(isinstance(policy, Mapping), "manifest API policy missing")
    require(policy.get("source_only") is True, "source_only policy changed")
    require(policy.get("full_goal_complete") is False, "full_goal_complete policy changed")
    require(policy.get("numeric_evaluation") is False, "numeric_evaluation policy changed")
    callback_counts = expected_event_counts(manifest)
    nexttask = manifest_value(manifest, "nexttask", "next_task")
    require(isinstance(nexttask, str) and nexttask, "manifest nexttask missing")
    return {
        "manifest": {"path": wire(manifest_path), "sha256": sha_bytes(manifest_raw), "bytes": len(manifest_raw)},
        "kernel": {**kernel_ast, "path": wire(kernel_path), "sha256": kernel_sha},
        "runner": {**runner_ast, "path": wire(runner_path), "sha256": runner_sha},
        "parent_report": artifact(parent_path, parent_raw, hash_only=True),
        "source_wrappers": wrapper_data,
        "source_wrapper_artifacts": wrapper_artifacts,
        "expected_counts": dict(expected_counts),
        "expected_new_flags": dict(expected_flags),
        "api_policy": dict(policy),
        "callback_event_counts": callback_counts,
        "nexttask": nexttask,
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


def build_input(parent: Mapping[str, Any], checked: Mapping[str, Any]) -> dict[str, Any]:
    value = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
        "source_wrappers": checked["source_wrappers"],
        "policy": {
            "source_only": True,
            "full_goal_complete": False,
            "numeric_evaluation": False,
        },
        "parent_report": {"path": wire(PARENT_REPORT), "sha256": PARENT_REPORT_SHA256, "hash_only": True},
        "saved_parent_counts": parent["parent_counts"],
    }
    no_base64_keys(value)
    return value


def build_source(parent: Mapping[str, Any], checked: Mapping[str, Any], input_raw: bytes, manifest_path: Path) -> dict[str, Any]:
    code2 = {"kernel": checked["kernel"], "runner": checked["runner"]}
    additions = {
        checked["kernel"]["path"]: checked["kernel"]["sha256"],
        checked["runner"]["path"]: checked["runner"]["sha256"],
        checked["manifest"]["path"]: checked["manifest"]["sha256"],
        wire(INPUT_PATH): sha_bytes(input_raw),
    }
    additions.update({item["path"]: item["sha256"] for item in checked["source_wrapper_artifacts"].values()})
    source_union = merge_source_hashes(parent["source_pins"], additions)
    return {
        "schema": SOURCE_SCHEMA,
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_report": parent["parent_report"],
        "parent_counts": parent["parent_counts"],
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "expected_failed_count_if_new_failure": PARENT_FAILED_COUNT + 7,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "parent_expected_counts": parent["parent_expected_counts"],
        "parent_expected_new_flags": parent["parent_expected_new_flags"],
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
        "input_spec": {"path": wire(INPUT_PATH), "sha256": sha_bytes(input_raw), "bytes": len(input_raw)},
        "input": json.loads(input_raw.decode("utf-8")),
        "source_wrappers": {
            name: {
                **checked["source_wrapper_artifacts"][name],
                "literal": checked["source_wrappers"][name]["literal"],
            }
            for name in checked["source_wrapper_artifacts"]
        },
        "expected_counts": checked["expected_counts"],
        "expected_new_flags": checked["expected_new_flags"],
        "newflags": checked["expected_new_flags"],
        "api_policy": checked["api_policy"],
        "callback_event_counts": checked["callback_event_counts"],
        "nexttask": checked["nexttask"],
        "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
        "policy": {"source_only": True, "full_goal_complete": False, "numeric_evaluation": False},
        "runtime_model_metadata": None,
        "no_project_imports_before_run": True,
        "numeric_evaluation": False,
    }


def prepare(manifest_value_path: str) -> dict[str, Any]:
    require(not INPUT_PATH.exists(), "input artifact already exists")
    require(not SOURCE_PATH.exists(), "source checkpoint already exists")
    require(not BRIEF_PATH.exists(), "pending brief already exists")
    require(not ENTERED_PATH.exists() and not JOURNAL_PATH.exists() and not RESULT_PATH.exists(), "attempt already entered")
    manifest_path, manifest_raw, manifest = load_manifest(manifest_value_path)
    checked = authenticate_manifest(manifest_path, manifest_raw, manifest)
    _, _, parent = load_parent()
    current_raw, current_brief = load_active_brief()
    input_value = build_input(parent, checked)
    input_raw = write_new(INPUT_PATH, input_value)
    preimage_path, preimage_raw = brief_preimage(current_raw)
    source = build_source(parent, checked, input_raw, manifest_path)
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
        "source_wrappers": checked["source_wrapper_artifacts"],
        "main_brief_preimage": {"path": wire(preimage_path), "sha256": sha_bytes(preimage_raw), "bytes": len(preimage_raw)},
        "journal": {"path": wire(JOURNAL_PATH), "schema": JOURNAL_SCHEMA},
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
                    "source_only": True,
                    "full_goal_complete": False,
                }
            },
            "goal_active": True,
            "full_goal_complete": False,
            "next_task": checked["nexttask"],
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
        "runtime_model_metadata": None,
    }


def load_source_for_run() -> tuple[bytes, dict[str, Any]]:
    source_raw, source = read_json(SOURCE_PATH)
    require(source.get("schema") == SOURCE_SCHEMA and source.get("component") == COMPONENT, "source checkpoint schema changed")
    require(source.get("parent_commit") == PARENT_COMMIT, "source parent commit changed")
    input_path = resolve_artifact(source["input_spec"]["path"], source["input_spec"]["sha256"])
    require(input_path.resolve() == INPUT_PATH.resolve(), "frozen input path changed")
    input_raw, input_value = read_json(input_path, source["input_spec"]["sha256"])
    require(input_value == source.get("input"), "frozen input literal changed")
    manifest_path = resolve_artifact(source["manifest"]["path"], source["manifest"]["sha256"])
    manifest_raw, manifest = read_json(manifest_path, source["manifest"]["sha256"])
    checked = authenticate_manifest(manifest_path, manifest_raw, manifest)
    require(checked["kernel"]["sha256"] == source["code2"]["kernel"]["sha256"], "frozen kernel changed")
    require(checked["runner"]["sha256"] == source["code2"]["runner"]["sha256"], "frozen runner changed")
    require(sha_bytes(input_raw) == source["input_spec"]["sha256"], "frozen input SHA changed")
    for key in (
        "source_pins",
        "protected_file_hashes",
        "legacy_dirty_path_hashes",
        "preserved_failed_artifact_hashes_union",
    ):
        authenticate_hash_map(source[key], f"frozen {key}")
    for path_value, digest in source["source_union_hashes"].items():
        path = resolve_artifact(path_value, digest)
        read_bytes(path, digest)
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
    spec = importlib.util.spec_from_file_location("_native_O3_mixed_radial_operator_entry", path)
    require(spec is not None and spec.loader is not None, "kernel import spec unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    require(sys.modules.get(spec.name) is module, "kernel was not registered in sys.modules")
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "frozen kernel has no derive API")
    return module


def callback_rows_callback() -> tuple[Any, list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []

    def emit(*args: Any) -> None:
        if len(args) == 1 and isinstance(args[0], Mapping):
            row = dict(args[0])
        elif len(args) == 2 and isinstance(args[0], str) and isinstance(args[1], Mapping):
            row = {"event": args[0], "payload": dict(args[1])}
        else:
            raise TypeError("emit expects event mapping or event name plus mapping")
        rows.append(row)
        append_journal(row)

    return emit, rows


def failure_result(
    source: Mapping[str, Any],
    error: BaseException,
    rows: list[dict[str, Any]],
    imported: bool,
    derive_entry_calls: int,
    derive_calls: int,
    partial: Mapping[str, Any] | None,
) -> dict[str, Any]:
    event_to_count = {
        "family_entered": "new_algebra_family_calls",
        "compiled_row": "new_compiled_operator_bindings",
        "unavailable_source": "unavailable_source_rows",
        "api_refusal": "actual_negative_api_refusals",
    }
    actual_counts = {
        str(name): 0
        for name in (source.get("expected_counts") if isinstance(source.get("expected_counts"), Mapping) else {})
    }
    for row in rows:
        count_name = event_to_count.get(row.get("event"))
        if count_name is not None:
            actual_counts[count_name] = actual_counts.get(count_name, 0) + 1
    for count_name in event_to_count.values():
        actual_counts.setdefault(count_name, 0)
    expected_flags = source.get("expected_new_flags")
    new_flags = {
        str(name): False
        for name in (expected_flags if isinstance(expected_flags, Mapping) else {})
    }
    partial_result = dict(partial) if isinstance(partial, Mapping) else None
    failure = {
        "kind": "manager_import_or_derive_failure",
        "error_type": type(error).__name__,
        "error": repr(error),
        "traceback": traceback.format_exc(),
        "kernel_imported": imported,
        "stopped": True,
        "callback_prefix": list(rows),
    }
    result = partial_result or {}
    result.update({
        "schema": result.get("schema", source.get("expected_result_schema", f"{COMPONENT}_result_v1")),
        "component": COMPONENT,
        "passed": False,
        "actual_counts": actual_counts,
        "new_flags": new_flags,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_imported": imported,
        "partial_result": partial_result,
        "payload": {
            "partial_result": partial_result,
            "recovered_callback_prefix": list(rows),
        },
        "failure": failure,
        "runtime_model_metadata": None,
    })
    return result


def wrap_result(source: Mapping[str, Any], kernel_result: Mapping[str, Any], *, imported: bool, derive_entry_calls: int, derive_calls: int, rows: list[dict[str, Any]], failure: Mapping[str, Any] | None = None) -> dict[str, Any]:
    passed = kernel_result.get("passed") is True and failure is None
    wrapper: dict[str, Any] = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "passed" if passed else "failed",
        "passed": passed,
        "metadata_only": False,
        "source_only": True,
        "full_goal_complete": False,
        "numeric_evaluation": False,
        "result": dict(kernel_result),
        "payload": dict(kernel_result),
        "actual_counts": kernel_result.get("actual_counts", {}),
        "new_flags": kernel_result.get("new_flags", {}),
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "import_attempts": 1,
        "kernel_imported": imported,
        "callback_events": list(rows),
        "callback_event_count": len(rows),
        "failure": dict(failure) if failure is not None else None,
        "runtime_model_metadata": None,
        "source_checkpoint_sha256": sha_bytes(SOURCE_PATH.read_bytes()),
    }
    return wrapper


def validate_saved(source: Mapping[str, Any], saved: Mapping[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    violations: list[str] = []
    if saved.get("schema") != RESULT_SCHEMA:
        violations.append("result schema")
    if saved.get("component") != COMPONENT:
        violations.append("result component")
    if saved.get("runtime_model_metadata") is not None:
        violations.append("runtime model metadata")
    if saved.get("passed") is True:
        expected_counts = source.get("expected_counts", {})
        actual_counts = saved.get("actual_counts")
        if not isinstance(actual_counts, Mapping) or set(actual_counts) != set(expected_counts):
            violations.append("actual count keys")
        elif any(actual_counts.get(key) != expected for key, expected in expected_counts.items()):
            violations.append("actual count values")
        expected_flags = source.get("expected_new_flags", {})
        actual_flags = saved.get("new_flags")
        if not isinstance(actual_flags, Mapping) or set(actual_flags) != set(expected_flags):
            violations.append("new flag keys")
        elif any(not isinstance(actual_flags.get(key), bool) for key in expected_flags):
            violations.append("new flag value types")
        if saved.get("derive_entry_calls") != 1 or saved.get("derive_calls") != 1 or saved.get("kernel_imported") is not True:
            violations.append("entry counts")
        event_counts: dict[str, int] = {}
        for row in rows:
            name = row.get("event")
            event_counts[name] = event_counts.get(name, 0) + 1
        if event_counts != source.get("callback_event_counts"):
            violations.append("callback event counts")
        if len(journal_rows()) != 74:
            violations.append("journal event count")
    return {
        "schema": f"{COMPONENT}_validation_v1",
        "component": COMPONENT,
        "passed": not violations and saved.get("passed") is True,
        "actual_counts": saved.get("actual_counts"),
        "callback_event_count": len(rows),
        "journal_event_count": len(journal_rows()),
        "violations": violations,
    }


def run() -> dict[str, Any]:
    source_raw, source = load_source_for_run()
    require(not ENTERED_PATH.exists(), "entered checkpoint already exists; reentry forbidden")
    require(not JOURNAL_PATH.exists() and not RESULT_PATH.exists(), "run checkpoints already exist; reentry forbidden")
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
    append_journal({
        "schema": JOURNAL_SCHEMA,
        "event": "manager_entry_before_import",
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_bytes(entered_raw),
        "import_attempts": 1,
        "derive_entry_calls": 0,
        "derive_calls": 0,
    })
    emit, rows = callback_rows_callback()
    imported = False
    derive_entry_calls = 0
    derive_calls = 0
    kernel_result: Mapping[str, Any] | None = None
    failure: Mapping[str, Any] | None = None
    try:
        module = import_kernel(KERNEL_PATH)
        imported = True
        derive_entry_calls = 1
        append_journal({
            "schema": JOURNAL_SCHEMA,
            "event": "manager_derive_entry",
            "import_attempts": 1,
            "derive_entry_calls": 1,
            "derive_calls": 0,
            "sys_modules_registered": True,
        })
        derive_calls = 1
        returned = module.derive(source["input"], emit)
        require(isinstance(returned, Mapping), "kernel derive returned a non-mapping")
        kernel_result = dict(returned)
    except BaseException as exc:
        partial = kernel_result if isinstance(kernel_result, Mapping) else None
        kernel_result = failure_result(source, exc, rows, imported, derive_entry_calls, derive_calls, partial)
        failure = kernel_result["failure"]
    wrapped = wrap_result(
        source,
        kernel_result or {},
        imported=imported,
        derive_entry_calls=derive_entry_calls,
        derive_calls=derive_calls,
        rows=rows,
        failure=failure,
    )
    # This is the durability boundary: no result validation occurs before it.
    result_raw = write_new(RESULT_PATH, wrapped)
    saved_raw, saved = read_json(RESULT_PATH, sha_bytes(result_raw))
    validation = validate_saved(source, saved, rows)
    return {
        "status": "run_complete",
        "component": COMPONENT,
        "attempt_status": saved.get("attempt_status"),
        "passed": bool(validation["passed"]),
        "result": {"path": wire(RESULT_PATH), "sha256": sha_bytes(saved_raw), "schema": saved.get("schema")},
        "actual_counts": saved.get("actual_counts"),
        "new_flags": saved.get("new_flags"),
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
        "kernel_imported": imported,
        "callback_event_count": len(rows),
        "journal_event_count": len(journal_rows()),
        "validation": validation,
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
