"""Durable runner for the native O3 Cartesian spatial-jet source leaf.

``prepare`` authenticates only frozen bytes, JSON literals, hashes and ASTs.
``run`` is the sole source consumer: it records entry before one import and
one ``derive`` call, journals each strict two-argument callback row before
retaining it, and saves the complete nested result before validation.

The runner deliberately contains no mathematical implementation.  Capsule
bytes and all expected counts/flags are supplied by the frozen manifest; in
particular, the capsule SHA is never a source-code placeholder.
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

COMPONENT = "native_O3_Cartesian_spatial_jet"
INPUT_SCHEMA = f"{COMPONENT}_input_v1"
MANIFEST_SCHEMA = f"{COMPONENT}_manifest_v1"
SOURCE_SCHEMA = f"{COMPONENT}_source_v1"
ENTERED_SCHEMA = f"{COMPONENT}_entered_v1"
RESULT_SCHEMA = f"{COMPONENT}_result_checkpoint_v1"
JOURNAL_SCHEMA = f"{COMPONENT}_journal_v1"

# The radial leaf is the immediate accepted parent.  Its publication was
# made on 3a4e..., then the current checkout advanced to d2bbc... for this
# new leaf.  These values are intentionally not derived from mutable files.
PARENT_REPORT = REPO / (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_native_O3_radial_variation.json"
)
PARENT_REPORT_SHA256 = "ccafb3840072efb5de6f525bbdeb17cc62b489c8c9fb170ecf30279c6578812c"
PARENT_REPORT_INTERNAL_PARENT_COMMIT = "3a4e6849db35b1a757022d3e9afc9e058721ab50"
CURRENT_HEAD_COMMIT = "d2bbc985429f343834ba30a7676dbd20fd6eb25a"
PARENT_COMMIT = CURRENT_HEAD_COMMIT
PARENT_SOURCE_PIN_COUNT = 2351
PARENT_CHECKPOINT_COUNT = 300
PARENT_FAILED_ARTIFACT_COUNT = 197
PARENT_PROTECTED_FILE_COUNT = 16
PARENT_DIRTY_PATH_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
PUBLISHED_CHECKPOINT_COUNT = 304
FAILED_COUNT_IF_NEW_FAILURE = 204
EXPECTED_SOURCE_UNION_COUNT = 2356
SOURCE_ARCHIVE_SHA256 = "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d"
SOURCE_ARCHIVE_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
)

MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
MAIN_BRIEF_SHA256 = "e9f7de2ae4e1e185f4d3ce70498e102fb241c3b6d18979072dd4bcdfa635e45b"
MAIN_BRIEF_BYTES = 92824
KERNEL_RELATIVE = Path("src/openai_ns_reconstruction/native_O3_Cartesian_spatial_jet.py")
RUNNER_RELATIVE = Path("scripts/run_native_O3_Cartesian_spatial_jet.py")
KERNEL_PATH = REPO / KERNEL_RELATIVE
RUNNER_PATH = REPO / RUNNER_RELATIVE
INPUT_PATH = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF_PATH = WORK / f"{COMPONENT}_brief_2026_10_27.json"
CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_PATH = CHECKPOINT_DIR / "source.json"
ENTERED_PATH = CHECKPOINT_DIR / "entered.json"
JOURNAL_PATH = CHECKPOINT_DIR / "journal.jsonl"
RESULT_PATH = CHECKPOINT_DIR / "result.json"

# This is the only new source wrapper.  Its SHA is authenticated from the
# manifest after the root creates and freezes the capsule.
CAPSULE_PATH = WORK / "native_O3_Cartesian_spatial_source_capsule_2026_10_27.json"
WRAPPER_NAMES = ("capsule",)

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
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdefABCDEF" for c in value)


def wire(path: str | Path) -> str:
    return Path(path).as_posix()


def encode(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def encode_compact(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def encode_brief(value: Any) -> bytes:
    raw = encode(value)
    if len(raw) < 100_000:
        return raw
    return encode_compact(value)


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
    temporary = path.with_name(path.name + ".Cartesian_spatial_jet_replace")
    with temporary.open("xb") as handle:
        handle.write(raw)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    return raw


def candidates(value: str | Path) -> list[Path]:
    path = Path(value)
    if path.is_absolute():
        return [path]
    return list(dict.fromkeys((Path.cwd() / path, REPO / path, WORK / path, REPO / "experiments/root_st073" / path)))


def resolve_artifact(value: str | Path, expected: str | None = None) -> Path:
    pool = candidates(value)
    if expected:
        for path in pool:
            if path.is_file() and sha_bytes(path.read_bytes()) == expected.lower():
                return path
    for path in pool:
        if path.is_file():
            return path
    raise FileNotFoundError(str(value))


def artifact(path: str | Path, raw: bytes | None = None, *, hash_only: bool = False) -> dict[str, Any]:
    path = Path(path)
    raw = path.read_bytes() if raw is None else raw
    item = {"path": wire(path), "sha256": sha_bytes(raw), "bytes": len(raw), "raw_bytes_exact": True, "raw_readonly": True}
    if hash_only:
        item["hash_only"] = True
    else:
        item["raw_bytes_base64"] = __import__("base64").b64encode(raw).decode("ascii")
    return item


def authenticate_spec(spec: Mapping[str, Any], label: str, expected_path: Path | None = None) -> tuple[Path, bytes, str]:
    require(isinstance(spec, Mapping), f"{label} spec missing")
    digest = spec.get("sha256")
    require(valid_sha(digest), f"{label} SHA missing")
    path = resolve_artifact(spec.get("path"), digest)
    raw = read_bytes(path, digest)
    if expected_path is not None:
        require(path.resolve() == expected_path.resolve(), f"{label} path changed")
    return path, raw, sha_bytes(raw)


def parse_ast(path: Path, raw: bytes, label: str) -> dict[str, Any]:
    try:
        tree = ast.parse(raw.decode("utf-8"), filename=str(path), mode="exec")
    except (SyntaxError, UnicodeDecodeError) as exc:
        raise ValueError(f"{label} AST parse failed: {exc}") from exc
    return {
        "path": wire(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "ast_parse": "ok",
        "node_count": sum(1 for _ in ast.walk(tree)),
        "top_level_definitions": [node.name for node in tree.body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))],
    }


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip()


def load_active_brief() -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    raw, brief = read_json(MAIN_BRIEF)
    require(len(raw) == MAIN_BRIEF_BYTES and sha_bytes(raw) == MAIN_BRIEF_SHA256, "active brief raw pin changed")
    require(brief.get("execution_rules") == EXECUTION_RULES and brief.get("execution_rules_preserved_exact") is True, "execution rules changed")
    require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False, "goal state changed")
    accepted = brief.get("last_accepted_mathematical_state")
    require(isinstance(accepted, Mapping) and accepted.get("component") == "native_O3_radial_variation", "accepted state changed")
    require(accepted.get("commit") == CURRENT_HEAD_COMMIT and accepted.get("head") == CURRENT_HEAD_COMMIT, "accepted HEAD changed")
    require(accepted.get("report_sha256") == PARENT_REPORT_SHA256, "accepted parent report changed")
    baseline = brief.get("accepted_progress")
    require(isinstance(baseline, Mapping) and baseline.get("commit") == "8f36f19e3c524b424a9286c46cc8e0b855f75307", "physical baseline changed")
    return raw, brief, {"accepted": dict(accepted), "baseline": dict(baseline), "execution_rules": list(EXECUTION_RULES)}


def update_brief(raw: bytes, brief: Mapping[str, Any], patch: Mapping[str, Any]) -> dict[str, Any]:
    require(brief.get("execution_rules") == EXECUTION_RULES, "execution rules changed before brief update")
    updated = dict(brief)
    for key, value in patch.items():
        if key == "in_progress" and isinstance(updated.get(key), dict) and isinstance(value, Mapping):
            merged = dict(updated[key]); merged.update(dict(value)); updated[key] = merged
        else:
            updated[key] = value
    updated["execution_rules"] = list(EXECUTION_RULES)
    updated["execution_rules_preserved_exact"] = True
    new_raw = encode_brief(updated)
    require(len(new_raw) < 100_000, "updated brief exceeds 100 KB")
    write_replace(MAIN_BRIEF, new_raw)
    return {"path": wire(MAIN_BRIEF), "sha256": sha_bytes(new_raw), "bytes": len(new_raw), "preimage_sha256": sha_bytes(raw)}


def hash_map(value: Any, label: str) -> dict[str, str]:
    require(isinstance(value, Mapping), f"{label} must be a mapping")
    out: dict[str, str] = {}
    for path, digest in value.items():
        require(isinstance(path, str) and valid_sha(digest), f"{label} malformed: {path}")
        out[path] = digest.lower()
    return out


def authenticate_hash_map(value: Any, label: str) -> dict[str, str]:
    mapping = hash_map(value, label)
    for path, digest in mapping.items():
        read_bytes(resolve_artifact(path, digest), digest)
    return mapping


def load_parent() -> tuple[bytes, dict[str, Any], dict[str, Any]]:
    raw, report = read_json(PARENT_REPORT, PARENT_REPORT_SHA256)
    require(report.get("component") == "native_O3_radial_variation", "parent component changed")
    require(report.get("passed") is True and report.get("attempt_status") == "passed", "parent status changed")
    require(report.get("source_only") is True and report.get("metadata_only") is False and report.get("full_goal_complete") is False, "parent scope changed")
    require(report.get("publication_metadata_only") is True and report.get("physical_layer_completed") is False, "parent publication scope changed")
    require(report.get("parent_commit") == PARENT_REPORT_INTERNAL_PARENT_COMMIT and report.get("prepublication_base_head") == PARENT_REPORT_INTERNAL_PARENT_COMMIT, "parent lineage changed")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT and report.get("source_union_count") == PARENT_SOURCE_PIN_COUNT, "parent source count changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT and report.get("preserved_failed_artifact_count") == PARENT_FAILED_ARTIFACT_COUNT, "parent checkpoint/failure count changed")
    require(report.get("parent_protected_file_count") == PARENT_PROTECTED_FILE_COUNT and report.get("parent_dirty_path_count") == PARENT_DIRTY_PATH_COUNT, "parent caps changed")
    flags = report.get("parent_false_flags")
    require(isinstance(flags, list) and len(flags) == PARENT_FALSE_FLAG_COUNT and report.get("parent_false_flag_values") == {name: False for name in flags}, "parent false flags changed")
    maps = {key: authenticate_hash_map(report.get(key), f"parent {key}") for key in ("source_pins", "protected_file_hashes", "legacy_dirty_path_hashes", "preserved_failed_artifact_hashes_union")}
    require(len(maps["source_pins"]) == PARENT_SOURCE_PIN_COUNT and len(maps["protected_file_hashes"]) == PARENT_PROTECTED_FILE_COUNT and len(maps["legacy_dirty_path_hashes"]) == PARENT_DIRTY_PATH_COUNT and len(maps["preserved_failed_artifact_hashes_union"]) == PARENT_FAILED_ARTIFACT_COUNT, "parent map count changed")
    retained = {key: report[key] for key in ("preserved_error_cap_metadata", "retained_construction_data", "retained_integral_term_inventory", "retained_cap_provenance")}
    context = {**maps, **retained, "parent_false_flags": list(flags), "parent_false_flag_values": dict(report["parent_false_flag_values"]), "parent_new_flags": dict(report.get("new_flags", {})), "parent_counts": {"source_pins": PARENT_SOURCE_PIN_COUNT, "checkpoints": PARENT_CHECKPOINT_COUNT, "failed": PARENT_FAILED_ARTIFACT_COUNT, "protected": PARENT_PROTECTED_FILE_COUNT, "dirty": PARENT_DIRTY_PATH_COUNT, "false_flags": PARENT_FALSE_FLAG_COUNT, "actual_counts": report.get("actual_counts", {})}, "parent_expected_counts": report.get("expected_counts", {}), "parent_expected_new_flags": report.get("expected_new_flags", {}), "parent_report": artifact(PARENT_REPORT, raw, hash_only=True)}
    return raw, report, context


def manifest_value(manifest: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in manifest:
            return manifest[name]
    return None


def load_manifest(value: str | Path) -> tuple[Path, bytes, dict[str, Any]]:
    path = resolve_artifact(value)
    raw, manifest = read_json(path)
    require(manifest.get("schema") == MANIFEST_SCHEMA and manifest.get("component") == COMPONENT, "manifest identity changed")
    require(manifest.get("source_archive_sha256", SOURCE_ARCHIVE_SHA256) == SOURCE_ARCHIVE_SHA256, "source archive SHA changed")
    return path, raw, manifest


def int_mapping(value: Any, label: str) -> dict[str, int]:
    require(isinstance(value, Mapping) and all(isinstance(k, str) and isinstance(v, int) and not isinstance(v, bool) and v >= 0 for k, v in value.items()), f"{label} malformed")
    return {str(k): int(v) for k, v in value.items()}


def map_events(manifest: Mapping[str, Any], callback_counts: Mapping[str, int]) -> dict[str, tuple[str, ...]]:
    raw = manifest_value(manifest, "count_event_map")
    require(isinstance(raw, Mapping), "count_event_map missing")
    names = set(callback_counts); out: dict[str, tuple[str, ...]] = {}
    for event, targets in raw.items():
        require(str(event) in names, f"unknown callback event: {event}")
        if isinstance(targets, str):
            targets = (targets,)
        require(isinstance(targets, (list, tuple)), f"event map malformed: {event}")
        out[str(event)] = tuple(str(v) for v in targets)
    for event in names:
        out.setdefault(event, ())
    return out


def capsule_status(capsule: Mapping[str, Any]) -> Mapping[str, Any]:
    status = capsule.get("status")
    # The frozen capsule uses the compact literal status string; accept the
    # mapping form as well so the manifest can preserve either exact schema
    # representation without changing the boundary checks.
    if isinstance(status, Mapping):
        require(status.get("static_literal_extraction_only") is True, "capsule static extraction status changed")
        view: Mapping[str, Any] = status
    else:
        require(status == "static_literal_extraction_only", "capsule static extraction status changed")
        view = capsule
    require(capsule.get("metadata_only") is True and capsule.get("mathematical_functions_called") == 0, "capsule metadata boundary changed")
    require(capsule.get("recipe_definitions_unexecuted") is True, "capsule recipe execution boundary changed")
    return view


def authenticate_manifest(path: Path, raw: bytes, manifest: Mapping[str, Any]) -> dict[str, Any]:
    require(manifest_value(manifest, "expected_head", "current_head", "head") == CURRENT_HEAD_COMMIT, "manifest expected HEAD changed")
    require(git_head() == CURRENT_HEAD_COMMIT, "repository HEAD changed")
    parent_spec = manifest_value(manifest, "parent_report", "parent_report_pin")
    _, _, parent_sha = authenticate_spec(parent_spec, "manifest parent report", PARENT_REPORT)
    require(parent_sha == PARENT_REPORT_SHA256, "manifest parent report SHA changed")
    artifacts = manifest.get("artifacts", manifest.get("code2", manifest))
    require(isinstance(artifacts, Mapping), "manifest artifacts missing")
    kernel_path, kernel_raw, kernel_sha = authenticate_spec(artifacts.get("kernel"), "manifest kernel", KERNEL_PATH)
    runner_path, runner_raw, runner_sha = authenticate_spec(artifacts.get("runner"), "manifest runner", RUNNER_PATH)
    kernel_ast = parse_ast(kernel_path, kernel_raw, "kernel")
    runner_ast = parse_ast(runner_path, runner_raw, "runner")
    require("derive" in kernel_ast["top_level_definitions"], "kernel derive API missing")
    wrappers = manifest.get("source_wrappers")
    if isinstance(wrappers, Mapping) and "capsule" not in wrappers and "source_capsule" in wrappers:
        wrappers = {"capsule": wrappers["source_capsule"]}
    require(isinstance(wrappers, Mapping) and set(wrappers) == set(WRAPPER_NAMES), "manifest wrappers changed")
    wrapper_data: dict[str, dict[str, Any]] = {}; wrapper_artifacts: dict[str, dict[str, Any]] = {}
    for name in WRAPPER_NAMES:
        path_value, wrapper_raw, digest = authenticate_spec(wrappers[name], f"wrapper {name}", CAPSULE_PATH)
        literal = json.loads(wrapper_raw.decode("utf-8"))
        require(isinstance(literal, dict), f"wrapper {name} is not an object")
        require(literal.get("schema") == "native_O3_Cartesian_spatial_source_capsule_v1", "capsule schema changed")
        capsule_status(literal)
        wrapper_data[name] = {"path": wire(path_value), "sha256": digest, "bytes": len(wrapper_raw), "literal": literal}
        wrapper_artifacts[name] = artifact(path_value, wrapper_raw, hash_only=True)
    expected_counts = int_mapping(manifest_value(manifest, "expected_counts"), "expected counts")
    flags = manifest_value(manifest, "expected_new_flags", "new_flags", "expected_flags")
    require(isinstance(flags, Mapping) and all(isinstance(v, bool) for v in flags.values()), "expected flags malformed")
    expected_flags = {str(k): bool(v) for k, v in flags.items()}
    callback_counts = int_mapping(manifest_value(manifest, "callback_event_counts", "expected_callback_event_counts"), "callback counts")
    event_map = map_events(manifest, callback_counts)
    policy = manifest_value(manifest, "api_policy", "policy")
    require(isinstance(policy, Mapping) and policy.get("source_only") is True and policy.get("metadata_only") is False and policy.get("full_goal_complete") is False and policy.get("numeric_evaluation") is False, "manifest policy changed")
    nexttask = manifest_value(manifest, "nexttask", "next_task")
    require(isinstance(nexttask, str) and nexttask and nexttask != COMPONENT, "manifest next task missing")
    schema = manifest_value(manifest, "kernel_result_schema", "result_schema") or f"{COMPONENT}_kernel_result_v1"
    next_plan = manifest_value(manifest, "next_unexecuted_plan", "next_plan", "future_unexecuted_plan", "next_independent_method_plan")
    if next_plan is not None:
        require(isinstance(next_plan, Mapping) and valid_sha(next_plan.get("sha256")), "future plan pointer malformed")
        plan_name = str(next_plan.get("path", ""))
        require("Cartesian_spatial_jet_and_volume_integrands_plan" not in Path(plan_name).name, "current completed plan cannot be next plan")
    require(sha_bytes(read_bytes(REPO / SOURCE_ARCHIVE_RELATIVE)) == SOURCE_ARCHIVE_SHA256, "raw source archive SHA changed")
    return {"manifest": {"path": wire(path), "sha256": sha_bytes(raw), "bytes": len(raw)}, "kernel": {**kernel_ast, "path": wire(kernel_path), "sha256": kernel_sha}, "runner": {**runner_ast, "path": wire(runner_path), "sha256": runner_sha}, "source_wrappers": wrapper_data, "source_wrapper_artifacts": wrapper_artifacts, "expected_counts": expected_counts, "expected_new_flags": expected_flags, "callback_event_counts": callback_counts, "count_event_map": event_map, "api_policy": dict(policy), "nexttask": nexttask, "next_plan": dict(next_plan) if isinstance(next_plan, Mapping) else None, "kernel_result_schema": schema, "proof_scope": manifest.get("proof_scope")}


def merge_hashes(parent: Mapping[str, str], additions: Mapping[str, str]) -> dict[str, str]:
    out = dict(parent)
    for path, digest in additions.items():
        require(path not in out or out[path].lower() == digest.lower(), f"conflicting source pin: {path}")
        out[path] = digest.lower()
    return out


def build_input(checked: Mapping[str, Any], context: Mapping[str, Any]) -> dict[str, Any]:
    capsule = checked["source_wrappers"]["capsule"]
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_commit": CURRENT_HEAD_COMMIT,
        "parent_report_internal_parent_commit": PARENT_REPORT_INTERNAL_PARENT_COMMIT,
        "capsule": capsule["literal"],
        "source_wrappers": {"capsule": {key: capsule[key] for key in ("path", "sha256", "bytes")}},
        "policy": {"metadata_only": False, "source_only": True, "full_goal_complete": False, "numeric_evaluation": False},
        "execution_rules": context["execution_rules"],
        "accepted_mathematical_state": context["accepted"],
        "physical_baseline": context["baseline"],
        "goal_active": True,
        "full_goal_complete": False,
        "next_task": checked["nexttask"],
        "next_plan": checked.get("next_plan"),
        "runtime_model_metadata": None,
    }


def build_source(parent: Mapping[str, Any], checked: Mapping[str, Any], input_raw: bytes, input_value: Mapping[str, Any], context: Mapping[str, Any], manifest_path: Path) -> dict[str, Any]:
    code2 = {"kernel": checked["kernel"], "runner": checked["runner"]}
    capsule = checked["source_wrappers"]["capsule"]
    additions = {code2["kernel"]["path"]: code2["kernel"]["sha256"], code2["runner"]["path"]: code2["runner"]["sha256"], checked["manifest"]["path"]: checked["manifest"]["sha256"], wire(INPUT_PATH): sha_bytes(input_raw), capsule["path"]: capsule["sha256"]}
    union = merge_hashes(parent["source_pins"], additions)
    require(len(union) == EXPECTED_SOURCE_UNION_COUNT, "source union count must be 2,356")
    return {
        "schema": SOURCE_SCHEMA,
        "component": COMPONENT,
        "parent_component": parent.get("component"),
        "parent_commit": CURRENT_HEAD_COMMIT,
        "parent_head": CURRENT_HEAD_COMMIT,
        "current_head_commit": CURRENT_HEAD_COMMIT,
        "parent_report_internal_parent_commit": PARENT_REPORT_INTERNAL_PARENT_COMMIT,
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
        "source_pins": dict(parent["source_pins"]),
        "source_union_hashes": union,
        "source_union_count": len(union),
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
        "input": dict(input_value),
        "source_wrappers": checked["source_wrappers"],
        "expected_counts": dict(checked["expected_counts"]),
        "expected_new_flags": dict(checked["expected_new_flags"]),
        "callback_event_counts": dict(checked["callback_event_counts"]),
        "count_event_map": {key: list(value) for key, value in checked["count_event_map"].items()},
        "callback_event_count": sum(checked["callback_event_counts"].values()),
        "journal_event_count": sum(checked["callback_event_counts"].values()) + 2,
        "nexttask": checked["nexttask"],
        "next_plan": checked.get("next_plan"),
        "kernel_result_schema": checked["kernel_result_schema"],
        "proof_scope": checked.get("proof_scope"),
        "source_archive_sha256": SOURCE_ARCHIVE_SHA256,
        "execution_rules": context["execution_rules"],
        "accepted_mathematical_state": context["accepted"],
        "physical_baseline": context["baseline"],
        "goal_active": True,
        "metadata_only": False,
        "source_only": True,
        "full_goal_complete": False,
        "numeric_evaluation": False,
        "runtime_model_metadata": None,
        "no_project_imports_before_run": True,
    }


def prepare(manifest_value_path: str) -> dict[str, Any]:
    for path in (INPUT_PATH, SOURCE_PATH, BRIEF_PATH, ENTERED_PATH, JOURNAL_PATH, RESULT_PATH):
        require(not path.exists(), f"artifact already exists: {wire(path)}")
    manifest_path, manifest_raw, manifest = load_manifest(manifest_value_path)
    checked = authenticate_manifest(manifest_path, manifest_raw, manifest)
    _, _, parent = load_parent()
    current_raw, current_brief, context = load_active_brief()
    input_value = build_input(checked, context)
    input_raw = write_new(INPUT_PATH, input_value)
    preimage_path = WORK / f"{COMPONENT}_main_resume_brief_preimage_{sha_bytes(current_raw)[:12]}.json"
    require(not preimage_path.exists(), "main brief preimage already exists")
    preimage_raw = write_new(preimage_path, current_raw)
    source = build_source(parent, checked, input_raw, input_value, context, manifest_path)
    source["main_brief_preimage"] = {"path": wire(preimage_path), "sha256": sha_bytes(preimage_raw), "bytes": len(preimage_raw), "raw_bytes_exact": True, "raw_readonly": True}
    source_raw = write_new(SOURCE_PATH, source)
    pending = {"schema": f"{COMPONENT}_brief_v1", "component": COMPONENT, "phase": "prepared", "current_component": COMPONENT, "source_checkpoint": {"path": wire(SOURCE_PATH), "sha256": sha_bytes(source_raw)}, "input_spec": source["input_spec"], "code2": source["code2"], "manifest": source["manifest"], "source_wrappers": source["source_wrappers"], "expected_counts": source["expected_counts"], "expected_new_flags": source["expected_new_flags"], "callback_event_counts": source["callback_event_counts"], "count_event_map": source["count_event_map"], "main_brief_preimage": source["main_brief_preimage"], "journal": {"path": wire(JOURNAL_PATH), "schema": JOURNAL_SCHEMA, "format": "jsonl_compact_one_object_per_line"}, "metadata_only": False, "source_only": True, "full_goal_complete": False, "runtime_model_metadata": None}
    pending_raw = write_new(BRIEF_PATH, pending)
    main_brief = update_brief(current_raw, current_brief, {"phase": "prepared", "current_component": COMPONENT, "in_progress": {COMPONENT: {"source": wire(SOURCE_PATH), "source_sha256": sha_bytes(source_raw), "input": wire(INPUT_PATH), "input_sha256": sha_bytes(input_raw), "code2": source["code2"], "manifest": source["manifest"], "metadata_only": False, "source_only": True, "full_goal_complete": False}}, "goal_active": True, "full_goal_complete": False, "next_task": source["nexttask"], "next_executable_task": source["nexttask"]})
    return {"status": "prepared", "component": COMPONENT, "source_checkpoint": {"path": wire(SOURCE_PATH), "sha256": sha_bytes(source_raw)}, "input": {"path": wire(INPUT_PATH), "sha256": sha_bytes(input_raw)}, "pending_brief": {"path": wire(BRIEF_PATH), "sha256": sha_bytes(pending_raw)}, "main_brief": main_brief, "source_union_count": source["source_union_count"], "expected_counts": source["expected_counts"], "expected_new_flags": source["expected_new_flags"], "callback_event_counts": source["callback_event_counts"], "callback_event_count": source["callback_event_count"], "journal_event_count": source["journal_event_count"], "metadata_only": False, "source_only": True, "runtime_model_metadata": None, "no_run_performed": True}


def load_source_for_run() -> tuple[bytes, dict[str, Any]]:
    source_raw, source = read_json(SOURCE_PATH)
    require(source.get("schema") == SOURCE_SCHEMA and source.get("component") == COMPONENT, "source checkpoint changed")
    require(source.get("parent_commit") == CURRENT_HEAD_COMMIT and source.get("parent_head") == CURRENT_HEAD_COMMIT and source.get("current_head_commit") == CURRENT_HEAD_COMMIT, "source HEAD changed")
    require(source.get("parent_report_internal_parent_commit") == PARENT_REPORT_INTERNAL_PARENT_COMMIT, "source lineage changed")
    require(source.get("source_archive_sha256") == SOURCE_ARCHIVE_SHA256, "source archive changed")
    preimage = source.get("main_brief_preimage"); require(isinstance(preimage, Mapping), "brief preimage missing"); read_bytes(resolve_artifact(preimage["path"], preimage["sha256"]), preimage["sha256"])
    input_path = resolve_artifact(source["input_spec"]["path"], source["input_spec"]["sha256"]); require(input_path.resolve() == INPUT_PATH.resolve(), "input path changed")
    _, input_value = read_json(input_path, source["input_spec"]["sha256"]); require(input_value == source.get("input"), "input literal changed")
    manifest_path = resolve_artifact(source["manifest"]["path"], source["manifest"]["sha256"]); manifest_raw, manifest = read_json(manifest_path, source["manifest"]["sha256"])
    checked = authenticate_manifest(manifest_path, manifest_raw, manifest)
    require(checked["kernel"]["sha256"] == source["code2"]["kernel"]["sha256"] and checked["runner"]["sha256"] == source["code2"]["runner"]["sha256"], "code2 changed")
    require(checked["expected_counts"] == source["expected_counts"] and checked["expected_new_flags"] == source["expected_new_flags"] and checked["callback_event_counts"] == source["callback_event_counts"], "manifest contract changed")
    require({key: list(value) for key, value in checked["count_event_map"].items()} == source["count_event_map"], "event map changed")
    require(checked.get("nexttask") == source.get("nexttask") and checked.get("next_plan") == source.get("next_plan"), "next task contract changed")
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})) == EXPECTED_SOURCE_UNION_COUNT, "source union count changed")
    for path, digest in source["source_union_hashes"].items():
        read_bytes(resolve_artifact(path, digest), digest)
    for name, wrapper in source["source_wrappers"].items():
        path = resolve_artifact(wrapper["path"], wrapper["sha256"]); raw = read_bytes(path, wrapper["sha256"]); require(json.loads(raw.decode("utf-8")) == wrapper.get("literal"), f"wrapper {name} changed"); capsule_status(wrapper["literal"])
    return source_raw, source


def append_journal(row: Mapping[str, Any]) -> None:
    # Compact JSONL is deliberately one complete object per line.  The old
    # radial runner's pretty encoder produced invalid line-oriented journal
    # validation after the result had already been saved; this leaf does not.
    raw = encode_compact(dict(row))
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
            value = json.loads(line.decode("utf-8")); require(isinstance(value, dict), "journal row is not an object"); rows.append(value)
    return rows


def import_kernel(path: Path) -> Any:
    require(path.resolve() != Path(__file__).resolve(), "kernel aliases runner")
    spec = importlib.util.spec_from_file_location("_native_O3_Cartesian_spatial_jet_entry", path)
    require(spec is not None and spec.loader is not None, "kernel import spec unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    require(sys.modules.get(spec.name) is module, "kernel module registration failed")
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "derive API missing")
    return module


def callback_rows_callback() -> tuple[Any, list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []

    def emit(*args: Any) -> None:
        require(len(args) == 2 and isinstance(args[0], str) and isinstance(args[1], Mapping), "callback must be exactly (event,row)")
        row = {"event": args[0], "row": dict(args[1])}
        append_journal(row)
        rows.append(row)

    return emit, rows


def count_prefix(source: Mapping[str, Any], rows: list[Mapping[str, Any]], failed: bool = False) -> dict[str, int]:
    actual = {str(key): 0 for key in source["expected_counts"]}
    for row in rows:
        for target in source["count_event_map"].get(row.get("event"), ()):
            require(target in actual, f"unknown count target: {target}")
            actual[target] += 1
    if failed and "failed_family_calls" in actual:
        actual["failed_family_calls"] = 1
    return actual


def false_flags(source: Mapping[str, Any]) -> dict[str, bool]:
    return {str(name): False for name in source["expected_new_flags"]}


def kernel_failure(source: Mapping[str, Any], exc: BaseException, rows: list[dict[str, Any]], imported: bool, entry_calls: int, calls: int) -> dict[str, Any]:
    failure = {"kind": "manager_import_or_derive_failure", "error_type": type(exc).__name__, "error": repr(exc), "traceback": traceback.format_exc(), "kernel_imported": imported, "derive_entry_calls": entry_calls, "derive_calls": calls, "callback_prefix": list(rows), "stopped": True}
    return {"schema": source.get("kernel_result_schema", f"{COMPONENT}_kernel_result_v1"), "component": COMPONENT, "passed": False, "metadata_only": False, "source_only": True, "full_goal_complete": False, "numeric_evaluation": False, "actual_counts": count_prefix(source, rows, True), "new_flags": false_flags(source), "proof_count": 0, "derive_entry_calls": entry_calls, "derive_calls": calls, "kernel_imported": imported, "failure": failure, "runtime_model_metadata": None}


def wrap_result(source: Mapping[str, Any], kernel_result: Mapping[str, Any], imported: bool, entry_calls: int, calls: int, rows: list[dict[str, Any]], manager_failure: Mapping[str, Any] | None) -> dict[str, Any]:
    passed = kernel_result.get("passed") is True and manager_failure is None
    return {"schema": RESULT_SCHEMA, "component": COMPONENT, "attempt_status": "passed" if passed else "failed", "passed": passed, "metadata_only": False, "source_only": True, "full_goal_complete": False, "numeric_evaluation": False, "result": dict(kernel_result), "actual_counts": kernel_result.get("actual_counts", count_prefix(source, rows, not passed)), "new_flags": kernel_result.get("new_flags", false_flags(source)), "proof_count": kernel_result.get("proof_count", 0), "derive_entry_calls": entry_calls, "derive_calls": calls, "kernel_imported": imported, "callback_events": list(rows), "callback_event_count": len(rows), "failure": dict(manager_failure) if manager_failure else kernel_result.get("failure"), "runtime_model_metadata": None, "source_checkpoint_sha256": sha_bytes(SOURCE_PATH.read_bytes())}


def validate_saved(source: Mapping[str, Any], saved: Mapping[str, Any], rows: list[dict[str, Any]], imported: bool, entry_calls: int, calls: int) -> dict[str, Any]:
    violations: list[str] = []
    if saved.get("schema") != RESULT_SCHEMA: violations.append("result schema")
    if saved.get("component") != COMPONENT: violations.append("result component")
    if saved.get("runtime_model_metadata") is not None: violations.append("runtime model metadata")
    allowed = set(source["callback_event_counts"])
    if any(row.get("event") not in allowed for row in rows): violations.append("callback event name")
    callback_counts = {name: sum(1 for row in rows if row.get("event") == name) for name in source["callback_event_counts"]}
    if saved.get("passed") is True:
        if saved.get("actual_counts") != source["expected_counts"]: violations.append("actual counts")
        if saved.get("new_flags") != source["expected_new_flags"]: violations.append("new flags")
        if saved.get("proof_count") != 1: violations.append("proof count")
        if not (imported and entry_calls == 1 and calls == 1): violations.append("lifecycle counts")
        if callback_counts != source["callback_event_counts"]: violations.append("callback counts")
        if len(journal_rows()) != source["journal_event_count"]: violations.append("journal count")
    else:
        if saved.get("actual_counts") != count_prefix(source, rows, True): violations.append("failed prefix counts")
        if saved.get("new_flags") != false_flags(source): violations.append("failed flags")
        if saved.get("proof_count") != 0: violations.append("failed proof count")
        nested = saved.get("result")
        if not isinstance(nested, Mapping) or nested.get("passed") is not False or (not isinstance(nested.get("failure"), Mapping) and not isinstance(saved.get("failure"), Mapping)): violations.append("nested failure")
        if len(journal_rows()) != 1 + int(imported) + len(rows): violations.append("failed journal prefix")
    return {"schema": f"{COMPONENT}_validation_v1", "component": COMPONENT, "passed": not violations and saved.get("passed") is True, "attempt_failed": saved.get("passed") is not True, "actual_counts": saved.get("actual_counts"), "new_flags": saved.get("new_flags"), "proof_count": saved.get("proof_count"), "callback_event_count": len(rows), "journal_event_count": len(journal_rows()), "violations": violations}


def run() -> dict[str, Any]:
    source_raw, source = load_source_for_run()
    require(not ENTERED_PATH.exists() and not JOURNAL_PATH.exists() and not RESULT_PATH.exists(), "run checkpoints already exist; reentry forbidden")
    entered = {"schema": ENTERED_SCHEMA, "component": COMPONENT, "source_checkpoint_sha256": sha_bytes(source_raw), "kernel": source["code2"]["kernel"], "runner": source["code2"]["runner"], "import_attempts": 1, "derive_entry_calls": 0, "derive_calls": 0, "sys_modules_registration_before_exec": True, "runtime_model_metadata": None, "entered_at_utc": now(), "do_not_reenter": True}
    entered_raw = write_new(ENTERED_PATH, entered)
    write_new(JOURNAL_PATH, b"")
    append_journal({"schema": JOURNAL_SCHEMA, "event": "manager_entry_before_import", "source_checkpoint_sha256": sha_bytes(source_raw), "entered_checkpoint_sha256": sha_bytes(entered_raw), "import_attempts": 1, "derive_entry_calls": 0, "derive_calls": 0})
    emit, rows = callback_rows_callback(); imported = False; entry_calls = 0; calls = 0; manager_failure = None
    kernel_result: Mapping[str, Any]
    try:
        module = import_kernel(KERNEL_PATH); imported = True; entry_calls = 1
        append_journal({"schema": JOURNAL_SCHEMA, "event": "manager_derive_entry", "import_attempts": 1, "derive_entry_calls": 1, "derive_calls": 0, "sys_modules_registered": True})
        calls = 1
        returned = module.derive(source["input"]["capsule"], callback=emit)
        require(isinstance(returned, Mapping), "derive returned a non-mapping")
        kernel_result = dict(returned)
    except BaseException as exc:
        kernel_result = kernel_failure(source, exc, rows, imported, entry_calls, calls); manager_failure = kernel_result["failure"]
    result_value = wrap_result(source, kernel_result, imported, entry_calls, calls, rows, manager_failure)
    result_raw = write_new(RESULT_PATH, result_value)
    saved_raw, saved = read_json(RESULT_PATH, sha_bytes(result_raw))
    validation = validate_saved(source, saved, rows, imported, entry_calls, calls)
    return {"status": "run_complete", "component": COMPONENT, "attempt_status": saved.get("attempt_status"), "passed": bool(validation["passed"]), "result": {"path": wire(RESULT_PATH), "sha256": sha_bytes(saved_raw), "schema": saved.get("schema")}, "actual_counts": saved.get("actual_counts"), "new_flags": saved.get("new_flags"), "proof_count": saved.get("proof_count"), "derive_entry_calls": entry_calls, "derive_calls": calls, "kernel_imported": imported, "callback_event_count": len(rows), "journal_event_count": len(journal_rows()), "validation": validation, "failure": saved.get("failure"), "runtime_model_metadata": None}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare"); prep.add_argument("--manifest", required=True)
    sub.add_parser("run")
    args = parser.parse_args(argv)
    try:
        receipt = prepare(args.manifest) if args.command == "prepare" else run()
    except Exception as exc:
        print(json.dumps({"status": "error", "error": f"{type(exc).__name__}: {exc}"}, sort_keys=True)); return 2
    print(json.dumps(receipt, sort_keys=True)); return 0 if receipt.get("passed", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
