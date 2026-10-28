"""Durable, unexecuted runner for the finite gauge-jet repair leaf.

The runner is a lifecycle manager, not a source-expression implementation.
``prepare`` authenticates a root-frozen manifest and records source/input
checkpoints.  A later root-only ``run`` persists entry, imports the canonical
kernel once, calls ``derive(data, emit)`` once, saves the complete returned
mapping before any payload validation, and never re-enters a failed leaf.
"""

from __future__ import annotations

import argparse
import ast
from collections.abc import Mapping
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import importlib.util
import json
import os
from pathlib import Path
import sys
import traceback
from typing import Any


REPO = Path(r"C:/Users/z5242/Documents/Codex/2026-09-13/before-doing-any-project-work-verify/work/st073-transition-next")
WORK = Path(r"C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work")
COMPONENT = "native_O3_eta_n1_finite_gauge_jet_repair"
INPUT_SCHEMA = f"{COMPONENT}_input_v1"
MANIFEST_SCHEMA = f"{COMPONENT}_manifest_v1"
SOURCE_SCHEMA = f"{COMPONENT}_source_v1"
ENTERED_SCHEMA = f"{COMPONENT}_entered_v1"
RESULT_SCHEMA = f"{COMPONENT}_result_checkpoint_v1"
JOURNAL_SCHEMA = f"{COMPONENT}_journal_v1"
CAPSULE_SCHEMA = f"{COMPONENT}_source_capsule_v1"

PARENT_REPORT = REPO / "experiments/root_st073" / "lei_ren_part1_paper_compliant_current_native_O3_eta_n1_core_source.json"
PARENT_REPORT_SHA256 = "330dccdef9ca61030370fa1d7873e75962c6d0f54ed9b16934d9e202327f3694"
PARENT_REPORT_BYTES = 57982348
PARENT_COMPONENT = "native_O3_eta_n1_core_source"
PARENT_INTERNAL_PARENT_COMMIT = "119d81ff2e2303b90457dbec436859b5cd922944"
PARENT_GRANDPARENT_COMMIT = "cddf394af87d3508f51718690d01a0ee4c4d349f"
CURRENT_HEAD_COMMIT = "e25ed616b2f93b53d47322b4117cf67bcd79041a"
PARENT_SOURCE_PIN_COUNT = 2468
PARENT_CHECKPOINT_COUNT = 352
PARENT_FAILED_ARTIFACT_COUNT = 211
PARENT_PROTECTED_FILE_COUNT = 16
PARENT_DIRTY_PATH_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
FROZEN_SOURCE_UNION_COUNT = 2473
PUBLISHED_SOURCE_UNION_COUNT = 2477
PUBLISHED_CHECKPOINT_COUNT = 356
FAILED_COUNT_IF_NEW_FAILURE = 218
PHYSICAL_BASELINE = "8f36f19e3c524b424a9286c46cc8e0b855f75307"

MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
MAIN_BRIEF_SHA256 = "cfd46de62152e57f63bcf9203e64d94c34a09104f5b92525b93852b88ed206fb"
MAIN_BRIEF_BYTES = 87219

KERNEL_RELATIVE = Path("src/openai_ns_reconstruction/native_O3_eta_n1_finite_gauge_jet_repair.py")
RUNNER_RELATIVE = Path("scripts/run_native_O3_eta_n1_finite_gauge_jet_repair.py")
KERNEL_PATH = REPO / KERNEL_RELATIVE
RUNNER_PATH = REPO / RUNNER_RELATIVE
INPUT_PATH = WORK / f"{COMPONENT}_input_2026_10_28.json"
MANIFEST_PATH = WORK / f"{COMPONENT}_manifest_2026_10_28.json"
BRIEF_PATH = WORK / f"{COMPONENT}_brief_2026_10_28.json"
CAPSULE_PATH = WORK / f"{COMPONENT}_source_capsule_2026_10_28.json"
CAPSULE_QA_PATH = WORK / f"{COMPONENT}_source_capsule_QA_2026_10_28.json"
CAPSULE_RAW_SHA256 = "de3973fbf0815dd2ba77a19d35e5dd1eb7eafe6ca6a08bf45b866a41d36e62bf"
CAPSULE_CANONICAL_SHA256 = "741fa4f436ed1a4ff13ca9c78ddf8fa2313b64ea6b7565d787b8cd74c4d02ba8"
CAPSULE_BYTES = 3145068
CAPSULE_QA_SHA256 = "bc3fac1b492561b58e7d8cae7c0fef07dadc83d7257a331fdfd2b2aa444bccbb"
CAPSULE_PROMOTION_RECEIPT_PATH = WORK / f"{COMPONENT}_capsule_v2_promotion_receipt_2026_10_28.json"
CAPSULE_PROMOTION_RECEIPT_SHA256 = "e98153fb00a6da3897112724eaac7851c973533c8404f1aa929b7823bda90f1e"
KERNEL_RAW_SHA256 = "7e1f77a602b13ef829850c6ca4a7e59febe3a36aa024f169730de7b9afe631bc"
KERNEL_BYTES = 56419
INTERVAL_PRECISION = 280
CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_PATH = CHECKPOINT_DIR / "source.json"
ENTERED_PATH = CHECKPOINT_DIR / "entered.json"
JOURNAL_PATH = CHECKPOINT_DIR / "journal.jsonl"
RESULT_PATH = CHECKPOINT_DIR / "result.json"
VALIDATION_PATH = WORK / f"{COMPONENT}_validation_receipt_2026_10_28.json"
BRIEF_SOFT_LIMIT = 97000

RHO_GRID = ("4",)
UNIFORM_RHO_GRID = ("0", "4")
ETA_GRID = ("-1", "-1/2", "-1/4", "0", "1/4", "1/2", "1")
UNIFORM_ETA_GRID = ("-1", "1")
PAYLOAD_RHO_LABELS = ("4", "uniform_rho")
PAYLOAD_ETA_LABELS = (*ETA_GRID, "uniform_eta")
NEW_MATH_FAMILIES = (
    "positive_denominator_gap",
    "rational_log_gauge_jet_ports",
    "angular_coefficient_contraction",
    "finite_angular_rhs_binding",
)
NEXT_SUCCESS_TASK = "source_native_O3_eta_n1_core_response_and_m15_binding"
PRE_ENTRY_TASK = "source_native_O3_eta_n1_finite_gauge_jet_repair"
GLOBAL_FALSE_KEYS = (
    "source_ftc_zero_identity",
    "source_C5_installation_claimed",
    "new_loop_domain_admission_claimed",
    "admissible_cone_admission_claimed",
    "cap_to_source_identity_proved",
    "quantitative_packet_error_reduction_proved",
    "numeric_values_computed",
    "numeric_containment_claimed",
    "physical_layer_completed",
    "global_five_moment_closure",
    "time_order_recursion_claimed",
    "old_operators_replayed",
    "old_pointwise_S_recovered",
    "old_mass_relation_proved",
    "old212_adapter_identity",
    "old212_quantclosure",
    "native_old_density_identity",
    "native_old_density_closure",
    "old_errors_removed",
    "old_error_terms_removed",
    "full_goal_complete",
    "actual_C3_C5_incoming_source_identity_proved",
    "actual_C3_C5_control_source_identity_proved",
    "actual_C3_C5_profile_identity_proved",
    "actual_C3_cap_to_C5_source_bound_proved",
    "compact_volume_cap_same_source_bound_proved",
)
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


def runtime_environment_metadata() -> dict[str, Any]:
    try:
        mpmath_version = importlib.metadata.version("mpmath")
    except importlib.metadata.PackageNotFoundError:
        mpmath_version = None
    return {"python_version": sys.version, "mpmath_version": mpmath_version, "interval_precision": INTERVAL_PRECISION, "runtime_model_metadata": None}


def sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def valid_sha(value: Any) -> bool:
    return isinstance(value, str) and len(value) == 64 and all(char in "0123456789abcdefABCDEF" for char in value)


def wire(path: str | Path) -> str:
    return Path(path).as_posix()


def encode(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def encode_compact(value: Any) -> bytes:
    return (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")


def write_new(path: Path, value: bytes | Any) -> bytes:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def write_replace(path: Path, value: bytes | Any) -> bytes:
    raw = value if isinstance(value, bytes) else encode(value)
    temporary = path.with_name(path.name + f".{COMPONENT}.replace")
    if temporary.exists():
        temporary.unlink()
    with temporary.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)
    return raw


def preserve_preimage(path: Path, raw: bytes, label: str) -> dict[str, Any]:
    target = WORK / f"{COMPONENT}_{label}_preimage_{sha_bytes(raw)[:16]}{path.suffix}"
    if target.exists():
        require(target.read_bytes() == raw, "preimage collision: " + wire(target))
    else:
        write_new(target, raw)
    return artifact(target, raw)


def candidates(value: str | Path) -> list[Path]:
    if isinstance(value, Mapping):
        value = value.get("path", "")
    path = Path(value)
    if path.is_absolute():
        return [path]
    return list(dict.fromkeys((Path.cwd() / path, WORK / path, REPO / path, REPO / "experiments/root_st073" / path)))


def resolve_artifact(value: str | Path, expected: str | None = None) -> Path:
    pool = candidates(value)
    if expected is not None:
        require(valid_sha(expected), f"expected SHA malformed for {value}")
        for path in pool:
            if path.is_file() and sha_bytes(path.read_bytes()) == expected.lower():
                return path
    for path in pool:
        if path.is_file():
            return path
    raise FileNotFoundError(str(value))


def read_bytes(value: str | Path, expected: str | None = None) -> bytes:
    path = resolve_artifact(value, expected)
    raw = path.read_bytes()
    if expected is not None:
        require(sha_bytes(raw) == expected.lower(), f"raw SHA changed: {wire(path)}")
    return raw


def read_json(value: str | Path, expected: str | None = None) -> tuple[Path, bytes, dict[str, Any]]:
    path = resolve_artifact(value, expected)
    raw = read_bytes(path, expected)
    parsed = json.loads(raw.decode("utf-8"))
    require(isinstance(parsed, dict), f"JSON object required: {wire(path)}")
    return path, raw, parsed


def artifact(path: str | Path, raw: bytes | None = None, *, expected: str | None = None, hash_only: bool = True) -> dict[str, Any]:
    resolved = resolve_artifact(path, expected)
    data = resolved.read_bytes() if raw is None else raw
    if expected is not None:
        require(sha_bytes(data) == expected.lower(), f"artifact SHA changed: {wire(resolved)}")
    item = {"path": wire(resolved), "sha256": sha_bytes(data), "bytes": len(data), "raw_bytes_exact": True, "raw_readonly": True}
    if hash_only:
        item["hash_only"] = True
    return item


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


def top_level_literals(raw: bytes, wanted: set[str], label: str) -> dict[str, Any]:
    tree = ast.parse(raw.decode("utf-8"), filename=label, mode="exec")
    values: dict[str, Any] = {}
    for node in tree.body:
        targets: list[ast.expr] = []
        if isinstance(node, ast.Assign):
            targets = [target for target in node.targets if isinstance(target, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target]
        for target in targets:
            if target.id not in wanted:
                continue
            try:
                values[target.id] = ast.literal_eval(node.value)
            except (ValueError, TypeError, SyntaxError) as exc:
                raise ValueError(f"{label} literal {target.id} is not a direct literal") from exc
    return values


def ast_kernel_contract(path: Path, raw: bytes) -> dict[str, Any]:
    # EXPECTED_COUNTS/EXPECTED_FLAGS and an API mapping containing SCHEMA are
    # derived expressions; only direct literal assignments are inspected.
    values = top_level_literals(raw, {"SUCCESS_COUNTS", "SUCCESS_FLAGS", "CALLBACK_EVENT_COUNTS", "COUNT_EVENT_MAP", "COMPONENT", "INTERVAL_DPS"}, "kernel")
    for name in ("SUCCESS_COUNTS", "SUCCESS_FLAGS", "CALLBACK_EVENT_COUNTS", "COUNT_EVENT_MAP"):
        require(isinstance(values.get(name), Mapping), f"kernel literal missing {name}")
    require(values.get("COMPONENT") == COMPONENT, "kernel component literal changed")
    if "INTERVAL_DPS" in values:
        require(values["INTERVAL_DPS"] == INTERVAL_PRECISION, "kernel interval precision changed")
    defs = parse_ast(path, raw, "kernel")["top_level_definitions"]
    require("derive" in defs, "kernel derive API missing")
    counts = dict(values["SUCCESS_COUNTS"])
    callbacks = dict(values["CALLBACK_EVENT_COUNTS"])
    require(counts.get("source_family_calls") == 1 and counts.get("new_math_family_calls") == 4 and counts.get("derive_calls") == 1 and counts.get("numeric_evaluation_calls") == 1 and counts.get("new_mathematical_proofs") == 1 and counts.get("rhs_rows_evaluated") == 48, "kernel count literals changed")
    require(sum(int(value) for value in callbacks.values()) == 11, "kernel callback literal count changed")
    return {"SUCCESS_COUNTS": counts, "SUCCESS_FLAGS": dict(values["SUCCESS_FLAGS"]), "CALLBACK_EVENT_COUNTS": callbacks, "COUNT_EVENT_MAP": dict(values["COUNT_EVENT_MAP"]), "SCHEMA": f"{COMPONENT}_kernel_result_v1", "CAPSULE_SCHEMA": CAPSULE_SCHEMA, "component_literal": values["COMPONENT"]}


def authenticate_spec(spec: Mapping[str, Any], label: str, expected_path: Path | None = None) -> tuple[Path, bytes, str]:
    require(isinstance(spec, Mapping), f"{label} spec missing")
    expected = spec.get("sha256")
    require(valid_sha(expected), f"{label} SHA missing")
    path = resolve_artifact(spec.get("path"), expected)
    raw = read_bytes(path, expected)
    if expected_path is not None:
        require(path.resolve() == expected_path.resolve(), f"{label} path changed")
    if "bytes" in spec:
        require(type(spec["bytes"]) is int and spec["bytes"] == len(raw), f"{label} byte count changed")
    return path, raw, sha_bytes(raw)


def capsule_status(capsule: Mapping[str, Any]) -> dict[str, Any]:
    require(capsule.get("schema") == CAPSULE_SCHEMA and capsule.get("component") == COMPONENT, "capsule identity changed")
    for name in ("parameters", "source_refs", "shape_contract", "proof_packet", "gauge_data", "mixed_phi", "parent_rhs_rows", "exact_view_registration", "source_atoms", "grid", "parent_rows_full_result", "repair_input_policy"):
        require(name in capsule, f"capsule field missing: {name}")
    require(any(name in capsule for name in ("retained_native_packet", "retainedpacket", "retained_packet", "retained_native_source_packet", "preservation")), "capsule retained packet field missing")
    require(isinstance(capsule["parameters"], Mapping) and isinstance(capsule["source_refs"], Mapping) and isinstance(capsule["shape_contract"], Mapping) and isinstance(capsule["proof_packet"], Mapping), "capsule contract fields malformed")
    for name, source_ref in capsule["source_refs"].items():
        if isinstance(source_ref, Mapping) and source_ref.get("path"):
            require(source_ref.get("path") != CAPSULE_PATH.as_posix(), f"capsule source_refs self-hash: {name}")
    shape = capsule["shape_contract"]
    require(shape.get("eta_count") == 8 and shape.get("rho_count") == 2 and shape.get("parent_rhs_row_count") == 48, "capsule shape contract changed")
    metadata = capsule.get("metadata")
    require(isinstance(metadata, Mapping) and metadata.get("metadata_only") is True and metadata.get("new_math_calls") == 0 and metadata.get("projection_arithmetic") is False and metadata.get("runtime_model_metadata") is None, "capsule metadata scope changed")
    return {"schema": capsule.get("schema"), "component": capsule.get("component"), "metadata_only": True, "math_recomputed": False}


def manifest_value(manifest: Mapping[str, Any], *names: str) -> Any:
    for name in names:
        if name in manifest:
            return manifest[name]
    raise ValueError("manifest field missing: " + "/".join(names))


def load_manifest(value: str | Path) -> tuple[Path, bytes, dict[str, Any]]:
    path, raw, manifest = read_json(value)
    require(manifest.get("schema") == MANIFEST_SCHEMA and manifest.get("component") == COMPONENT, "manifest identity changed")
    return path, raw, manifest


def load_parent() -> dict[str, Any]:
    _, raw, report = read_json(PARENT_REPORT, PARENT_REPORT_SHA256)
    require(len(raw) == PARENT_REPORT_BYTES, "parent report byte count changed")
    require(report.get("component") == PARENT_COMPONENT and report.get("passed") is True and report.get("attempt_status") == "passed", "parent report status changed")
    require(report.get("parent_commit") == PARENT_INTERNAL_PARENT_COMMIT, "parent internal commit changed")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT and report.get("source_union_count") == PARENT_SOURCE_PIN_COUNT and report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT and report.get("preserved_failed_artifact_count") == PARENT_FAILED_ARTIFACT_COUNT, "parent counts changed")
    require(isinstance(report.get("protected_file_hashes"), Mapping) and len(report["protected_file_hashes"]) == PARENT_PROTECTED_FILE_COUNT, "parent protected map changed")
    require(isinstance(report.get("legacy_dirty_path_hashes"), Mapping) and len(report["legacy_dirty_path_hashes"]) == PARENT_DIRTY_PATH_COUNT, "parent dirty map changed")
    require(isinstance(report.get("source_pins"), Mapping) and len(report["source_pins"]) == PARENT_SOURCE_PIN_COUNT, "parent source pin map changed")
    source_pins = {str(path): str(expected).lower() for path, expected in report["source_pins"].items()}
    for path, expected in source_pins.items():
        require(valid_sha(expected), "parent source pin SHA malformed")
        read_bytes(path, expected)
    for label in ("protected_file_hashes", "legacy_dirty_path_hashes"):
        for path, expected in report[label].items():
            require(valid_sha(expected), f"parent {label} contains malformed SHA")
            read_bytes(path, expected)
    require(report.get("finite_rhs_packet_admitted") is False and report.get("regular_n1_response_ready") is False, "parent finite admission status changed")
    return {"component": report.get("component"), "report": artifact(PARENT_REPORT, raw, expected=PARENT_REPORT_SHA256), "sha256": PARENT_REPORT_SHA256, "bytes": len(raw), "parent_commit": report.get("parent_commit"), "protected_file_count": len(report["protected_file_hashes"]), "dirty_path_count": len(report["legacy_dirty_path_hashes"]), "failed_artifact_count": report.get("preserved_failed_artifact_count"), "source_pin_count": report.get("source_pin_count"), "checkpoint_count": report.get("checkpoint_count"), "global_false_count": report.get("global_false_count", PARENT_FALSE_FLAG_COUNT), "finite_rhs_packet_admitted": False, "regular_n1_response_ready": False, "source_pins": source_pins, "protected_file_hashes": dict(report["protected_file_hashes"]), "legacy_dirty_path_hashes": dict(report["legacy_dirty_path_hashes"]), "preserved_failed_artifact_hashes_union": dict(report.get("preserved_failed_artifact_hashes_union", {})), "global_false_flags": dict(report.get("global_false_flags", {})), "parent_false_flags": list(report.get("parent_false_flags", [])), "parent_false_flag_values": dict(report.get("parent_false_flag_values", {}))}


def build_source(manifest: Mapping[str, Any], manifest_path: Path, manifest_raw: bytes, input_path: Path, input_raw: bytes, input_value: Mapping[str, Any], kernel_ast: Mapping[str, Any], runner_ast: Mapping[str, Any], parent: Mapping[str, Any], capsule_art: Mapping[str, Any], kernel_art: Mapping[str, Any], runner_art: Mapping[str, Any], main_brief_preimage: Mapping[str, Any], environment_metadata: Mapping[str, Any]) -> dict[str, Any]:
    source_union = dict(parent["source_pins"])
    for item in (kernel_art, runner_art, artifact(input_path, input_raw), artifact(manifest_path, manifest_raw), capsule_art):
        path = str(item["path"])
        sha = str(item["sha256"]).lower()
        if path in source_union:
            require(source_union[path] == sha, "source union path has conflicting SHA: " + path)
        source_union[path] = sha
    require(len(source_union) == FROZEN_SOURCE_UNION_COUNT, "frozen source union count changed")
    source = {
        "schema": SOURCE_SCHEMA,
        "component": COMPONENT,
        "metadata_only": False,
        "source_only": True,
        "numeric_evaluation": True,
        "numeric_values_computed": False,
        "current_head": CURRENT_HEAD_COMMIT,
        "parent_commit": CURRENT_HEAD_COMMIT,
        "parent_report": dict(parent["report"]),
        "parent_component": PARENT_COMPONENT,
        "physical_baseline": PHYSICAL_BASELINE,
        "last_accepted_mathematical_state": {"commit": CURRENT_HEAD_COMMIT, "component": PARENT_COMPONENT, "conditional": True, "source_only": True, "metadata_only": False, "finite_rhs_packet_admitted": False, "regular_n1_response_ready": False},
        "source_checkpoint_input": artifact(input_path, input_raw),
        "source_checkpoint_manifest": artifact(manifest_path, manifest_raw),
        "source_wrappers": {"capsule": dict(capsule_art)},
        "main_brief_preimage": dict(main_brief_preimage),
        "environment_metadata": dict(environment_metadata),
        "input": dict(input_value),
        "kernel": dict(kernel_art),
        "runner": dict(runner_art),
        "code2": {"kernel": dict(kernel_ast), "runner": dict(runner_ast)},
        "expected_counts": dict(manifest["expected_counts"]),
        "expected_new_flags": dict(manifest["expected_new_flags"]),
        "callback_event_counts": dict(manifest["callback_event_counts"]),
        "count_event_map": dict(manifest["count_event_map"]),
        "kernel_result_schema": manifest["kernel_result_schema"],
        "capsule_schema": CAPSULE_SCHEMA,
        "new_math_families": list(NEW_MATH_FAMILIES),
        "execution_rules": list(EXECUTION_RULES),
        "execution_rules_preserved_exact": True,
        "global_false_keys": list(GLOBAL_FALSE_KEYS),
        "next_task": PRE_ENTRY_TASK,
        "next_task_on_success": NEXT_SUCCESS_TASK,
        "continuation_requires_checkpoint_inspection": True,
        "finite_rhs_packet_admitted": False,
        "finite_g5_uniform_C0_C1_admitted": False,
        "regular_n1_response_ready": False,
        "old_error_terms_retained": True,
        "do_not_reenter": True,
    }
    source["source_pins"] = dict(source_union)
    source["source_union_hashes"] = dict(source_union)
    source["parent_source_pins"] = dict(parent["source_pins"])
    source["source_pin_count"] = FROZEN_SOURCE_UNION_COUNT
    source["source_union_count"] = FROZEN_SOURCE_UNION_COUNT
    return source


def prepare(manifest_value_path: str) -> dict[str, Any]:
    manifest_path, manifest_raw, manifest = load_manifest(manifest_value_path)
    parent = load_parent()
    brief_path, brief_raw, brief = read_json(MAIN_BRIEF, MAIN_BRIEF_SHA256)
    require(brief_path.resolve() == MAIN_BRIEF.resolve() and len(brief_raw) == MAIN_BRIEF_BYTES, "main brief identity changed")
    require(brief.get("execution_rules") == EXECUTION_RULES and brief.get("execution_rules_preserved_exact") is True, "main brief execution rules changed")
    require(not SOURCE_PATH.exists() and not BRIEF_PATH.exists() and not ENTERED_PATH.exists() and not RESULT_PATH.exists(), "source or pending run checkpoint already exists; freeze is single-entry")
    main_brief_preimage = preserve_preimage(MAIN_BRIEF, brief_raw, "main_brief")
    environment_metadata = runtime_environment_metadata()
    input_path, input_raw, input_value = read_json(manifest_value(manifest, "input"))
    require(input_value.get("schema") == INPUT_SCHEMA and input_value.get("component") == COMPONENT, "input identity changed")
    kernel_path, kernel_raw, _ = authenticate_spec(manifest["artifacts"]["kernel"], "kernel", KERNEL_PATH)
    runner_path, runner_raw, _ = authenticate_spec(manifest["artifacts"]["runner"], "runner", RUNNER_PATH)
    require(len(kernel_raw) == KERNEL_BYTES and sha_bytes(kernel_raw) == KERNEL_RAW_SHA256, "kernel final raw pin changed")
    capsule_path, capsule_raw, _ = authenticate_spec(manifest["source_wrappers"]["capsule"], "capsule", CAPSULE_PATH)
    require(len(capsule_raw) == CAPSULE_BYTES and sha_bytes(capsule_raw) == CAPSULE_RAW_SHA256, "capsule raw pin changed")
    capsule = json.loads(capsule_raw.decode("utf-8"))
    require(isinstance(capsule, Mapping), "capsule JSON object missing")
    capsule_status(capsule)
    capsule_ref = manifest["source_wrappers"]["capsule"]
    canonical_capsule = json.dumps(capsule, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    require(capsule_ref.get("canonical_sha256") == sha_bytes(canonical_capsule) and capsule_ref.get("canonical_bytes") == len(canonical_capsule), "capsule canonical hash changed")
    require(capsule_ref.get("canonical_sha256") == CAPSULE_CANONICAL_SHA256, "capsule canonical pin changed")
    qa_path, qa_raw, qa = read_json(CAPSULE_QA_PATH, CAPSULE_QA_SHA256)
    require(isinstance(qa, Mapping) and qa.get("metadata_only") is True and qa.get("new_math_calls") == 0, "capsule QA metadata scope changed")
    qa_ref = capsule_ref.get("qa")
    require(isinstance(qa_ref, Mapping) and qa_ref.get("path") == qa_path.as_posix() and qa_ref.get("sha256") == CAPSULE_QA_SHA256 and qa_ref.get("bytes") == len(qa_raw), "capsule QA pointer changed")
    _, promotion_raw, promotion = read_json(CAPSULE_PROMOTION_RECEIPT_PATH, CAPSULE_PROMOTION_RECEIPT_SHA256)
    require(promotion.get("metadata_only") is True and promotion.get("new_math_calls") == 0 and promotion.get("source_entry_performed") is False, "capsule promotion receipt scope changed")
    kernel_ast = ast_kernel_contract(kernel_path, kernel_raw)
    runner_ast = parse_ast(runner_path, runner_raw, "runner")
    require(runner_ast["path"] == wire(RUNNER_PATH), "runner canonical path changed")
    require(tuple(manifest.get("new_math_families", ())) == NEW_MATH_FAMILIES, "new family contract changed")
    require(list(manifest.get("execution_rules", ())) == EXECUTION_RULES and manifest.get("execution_rules_preserved_exact") is True, "execution rules changed")
    require(dict(manifest["expected_counts"]) == kernel_ast["SUCCESS_COUNTS"], "manifest/kernel count contract changed")
    require(dict(manifest["expected_new_flags"]) == kernel_ast["SUCCESS_FLAGS"], "manifest/kernel flag contract changed")
    require(dict(manifest["callback_event_counts"]) == kernel_ast["CALLBACK_EVENT_COUNTS"], "manifest/kernel callback contract changed")
    capsule_source_art = artifact(capsule_path, capsule_raw)
    capsule_source_art.update({"canonical_sha256": sha_bytes(canonical_capsule), "canonical_bytes": len(canonical_capsule), "qa": artifact(qa_path, qa_raw), "promotion_receipt": artifact(CAPSULE_PROMOTION_RECEIPT_PATH, promotion_raw)})
    source = build_source(manifest, manifest_path, manifest_raw, input_path, input_raw, input_value, kernel_ast, runner_ast, parent, capsule_source_art, artifact(kernel_path, kernel_raw), artifact(runner_path, runner_raw), main_brief_preimage, environment_metadata)
    source_raw = write_new(SOURCE_PATH, source)
    pending = {
        "schema": f"{COMPONENT}_pending_brief_v1",
        "component": COMPONENT,
        "status": "prepared_unexecuted",
        "metadata_only": True,
        "new_math_calls": 0,
        "source_checkpoint": artifact(SOURCE_PATH, source_raw),
        "manifest": artifact(manifest_path, manifest_raw),
        "input": artifact(input_path, input_raw),
        "parent_report": dict(parent["report"]),
        "main_brief_preimage": dict(main_brief_preimage),
        "environment_metadata": dict(environment_metadata),
        "main_brief": artifact(MAIN_BRIEF, brief_raw),
        "source_union_count": FROZEN_SOURCE_UNION_COUNT,
        "next_task": PRE_ENTRY_TASK,
        "next_executable_task": PRE_ENTRY_TASK,
        "next_task_on_success": NEXT_SUCCESS_TASK,
        "continuation_requires_checkpoint_inspection": True,
        "finite_rhs_packet_admitted": False,
        "regular_n1_response_ready": False,
        "pre_entry": True,
    }
    pending_raw = write_new(BRIEF_PATH, pending)
    prepared = json.loads(json.dumps(brief, ensure_ascii=False))
    prepared.update({
        "component": COMPONENT,
        "current_component": COMPONENT,
        "phase": "prepared_unexecuted_finite_gauge_jet_repair",
        "publication_status": "prepared_unexecuted",
        "source_checkpoint": artifact(SOURCE_PATH, source_raw),
        "pending_brief": artifact(BRIEF_PATH, pending_raw),
        "main_brief_preimage": dict(main_brief_preimage),
        "next_task": PRE_ENTRY_TASK,
        "next_executable_task": PRE_ENTRY_TASK,
        "next_task_on_success": NEXT_SUCCESS_TASK,
        "continuation_requires_checkpoint_inspection": True,
        "next_unexecuted_plan": {"component": PRE_ENTRY_TASK, "executed": False, "metadata_only": True, "new_math_calls": 0},
        "do_not_reenter": True,
        "environment_metadata": dict(environment_metadata),
        "goal_active": True,
        "full_goal_complete": False,
        "finite_rhs_packet_admitted": False,
        "finite_g5_uniform_C0_C1_admitted": False,
        "regular_n1_response_ready": False,
        "actual_n1_solution_materialized": False,
        "regular_n1_binding_installed": False,
        "source_pin_count": FROZEN_SOURCE_UNION_COUNT,
        "source_union_count": FROZEN_SOURCE_UNION_COUNT,
        "checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "preserved_failed_artifact_count": PARENT_FAILED_ARTIFACT_COUNT,
        "protected_file_hashes": dict(parent["protected_file_hashes"]),
        "legacy_dirty_path_hashes": dict(parent["legacy_dirty_path_hashes"]),
        "preserved_failed_artifact_hashes_union": dict(parent.get("preserved_failed_artifact_hashes_union", {})),
        "global_false_count": PARENT_FALSE_FLAG_COUNT,
        "global_false_flags": dict(parent.get("global_false_flags", {})),
        "last_accepted_mathematical_state": {"commit": CURRENT_HEAD_COMMIT, "component": PARENT_COMPONENT, "conditional": True, "source_only": True, "metadata_only": False, "finite_rhs_packet_admitted": False, "regular_n1_response_ready": False},
    })
    prepared_raw = encode(prepared)
    if len(prepared_raw) >= BRIEF_SOFT_LIMIT:
        prepared_raw = encode_compact(prepared)
    require(len(prepared_raw) < BRIEF_SOFT_LIMIT, "prepared main brief exceeds 97 KB")
    write_replace(MAIN_BRIEF, prepared_raw)
    pending["main_brief"] = artifact(MAIN_BRIEF, prepared_raw)
    pending["prepared_main_brief"] = artifact(MAIN_BRIEF, prepared_raw)
    pending_raw = write_replace(BRIEF_PATH, pending)
    return {"schema": f"{COMPONENT}_prepare_receipt_v1", "component": COMPONENT, "status": "prepared_unexecuted", "metadata_only": True, "new_math_calls": 0, "source_checkpoint": artifact(SOURCE_PATH, source_raw), "pending_brief": artifact(BRIEF_PATH, pending_raw), "main_brief": artifact(MAIN_BRIEF, prepared_raw), "main_brief_preimage": dict(main_brief_preimage), "kernel": artifact(kernel_path, kernel_raw), "runner": artifact(runner_path, runner_raw), "capsule": artifact(capsule_path, capsule_raw), "parent_report": dict(parent["report"]), "input": artifact(input_path, input_raw), "source_union_count": FROZEN_SOURCE_UNION_COUNT, "ast_only": True, "kernel_executed": False, "runner_executed": False, "runtime_model_metadata": None}


def load_source_for_run() -> tuple[bytes, dict[str, Any]]:
    source_path, raw, source = read_json(SOURCE_PATH)
    require(source_path.resolve() == SOURCE_PATH.resolve(), "source checkpoint path changed")
    require(source.get("schema") == SOURCE_SCHEMA and source.get("component") == COMPONENT, "source checkpoint identity changed")
    require(source.get("execution_rules") == EXECUTION_RULES and source.get("execution_rules_preserved_exact") is True, "source execution rules changed")
    require(source.get("finite_rhs_packet_admitted") is False and source.get("regular_n1_response_ready") is False, "source admission scope changed before run")
    environment_metadata = source.get("environment_metadata")
    require(isinstance(environment_metadata, Mapping) and environment_metadata.get("interval_precision") == INTERVAL_PRECISION and environment_metadata.get("runtime_model_metadata") is None and isinstance(environment_metadata.get("python_version"), str), "source environment metadata changed")
    union = source.get("source_union_hashes")
    require(isinstance(union, Mapping) and len(union) == FROZEN_SOURCE_UNION_COUNT, "source union map changed")
    for path, expected in union.items():
        require(valid_sha(expected), "source union SHA malformed")
        read_bytes(path, expected)
    require(source.get("source_pins") == dict(union) and source.get("source_pin_count") == FROZEN_SOURCE_UNION_COUNT and source.get("source_union_count") == FROZEN_SOURCE_UNION_COUNT, "source pin aliases changed")
    parent_source_pins = source.get("parent_source_pins")
    require(isinstance(parent_source_pins, Mapping) and len(parent_source_pins) == PARENT_SOURCE_PIN_COUNT and all(union.get(path) == expected for path, expected in parent_source_pins.items()), "parent source pin projection changed")
    for name, expected_path in (("kernel", KERNEL_PATH), ("runner", RUNNER_PATH)):
        spec = source.get(name)
        require(isinstance(spec, Mapping) and spec.get("path") == wire(expected_path) and valid_sha(spec.get("sha256")), f"source {name} pointer changed")
        code_raw = read_bytes(spec["path"], spec["sha256"])
        if name == "kernel":
            require(len(code_raw) == KERNEL_BYTES and sha_bytes(code_raw) == KERNEL_RAW_SHA256, "source kernel final raw pin changed")
    input_spec = source.get("source_checkpoint_input")
    manifest_spec = source.get("source_checkpoint_manifest")
    require(isinstance(input_spec, Mapping) and isinstance(manifest_spec, Mapping), "source input/manifest pointers missing")
    _, input_raw, input_value = read_json(input_spec["path"], input_spec["sha256"])
    _, manifest_raw, manifest = read_json(manifest_spec["path"], manifest_spec["sha256"])
    require(input_value.get("schema") == INPUT_SCHEMA and input_value.get("component") == COMPONENT, "source input identity changed")
    require(manifest.get("schema") == MANIFEST_SCHEMA and manifest.get("component") == COMPONENT, "source manifest identity changed")
    capsule_spec = source.get("source_wrappers", {}).get("capsule") if isinstance(source.get("source_wrappers"), Mapping) else None
    require(isinstance(capsule_spec, Mapping) and capsule_spec.get("path") == wire(CAPSULE_PATH) and capsule_spec.get("sha256") == CAPSULE_RAW_SHA256 and capsule_spec.get("bytes") == CAPSULE_BYTES, "source capsule pointer changed")
    _, capsule_raw, capsule = read_json(capsule_spec["path"], CAPSULE_RAW_SHA256)
    capsule_status(capsule)
    canonical_capsule = json.dumps(capsule, ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    require(sha_bytes(canonical_capsule) == CAPSULE_CANONICAL_SHA256 and capsule_spec.get("canonical_sha256") == CAPSULE_CANONICAL_SHA256, "source capsule canonical pin changed")
    read_bytes(CAPSULE_QA_PATH, CAPSULE_QA_SHA256)
    read_bytes(CAPSULE_PROMOTION_RECEIPT_PATH, CAPSULE_PROMOTION_RECEIPT_SHA256)
    preimage = source.get("main_brief_preimage")
    require(isinstance(preimage, Mapping) and preimage.get("sha256") == MAIN_BRIEF_SHA256 and preimage.get("bytes") == MAIN_BRIEF_BYTES, "main brief preimage pointer changed")
    read_bytes(preimage["path"], MAIN_BRIEF_SHA256)
    require(len(input_raw) == int(input_spec.get("bytes", len(input_raw))) and len(manifest_raw) == int(manifest_spec.get("bytes", len(manifest_raw))) and len(capsule_raw) == CAPSULE_BYTES, "source checkpoint bytes changed")
    brief_path, brief_raw, brief = read_json(MAIN_BRIEF)
    require(brief_path.resolve() == MAIN_BRIEF.resolve() and brief.get("component") == COMPONENT and brief.get("phase") == "prepared_unexecuted_finite_gauge_jet_repair", "prepared main brief changed before entry")
    require(brief.get("source_checkpoint", {}).get("sha256") == sha_bytes(raw) and brief.get("execution_rules") == EXECUTION_RULES and brief.get("execution_rules_preserved_exact") is True and brief.get("do_not_reenter") is True and brief.get("goal_active") is True and brief.get("full_goal_complete") is False, "prepared main brief safety flags changed")
    brief_preimage = brief.get("main_brief_preimage")
    require(isinstance(brief_preimage, Mapping) and brief_preimage.get("sha256") == preimage.get("sha256") and brief_preimage.get("bytes") == preimage.get("bytes"), "prepared main brief preimage alias changed")
    return raw, source


def append_journal(row: Mapping[str, Any]) -> None:
    with JOURNAL_PATH.open("ab") as stream:
        stream.write(encode_compact(dict(row)))
        stream.write(b"\n")
        stream.flush()
        os.fsync(stream.fileno())


def journal_rows() -> list[dict[str, Any]]:
    if not JOURNAL_PATH.exists():
        return []
    rows: list[dict[str, Any]] = []
    for line in JOURNAL_PATH.read_bytes().splitlines():
        if line:
            value = json.loads(line.decode("utf-8"))
            require(isinstance(value, dict), "journal row is not an object")
            rows.append(value)
    return rows


def import_kernel(path: Path) -> Any:
    require(path.resolve() != Path(__file__).resolve(), "kernel aliases runner")
    spec = importlib.util.spec_from_file_location(f"_{COMPONENT}_entry", path)
    require(spec is not None and spec.loader is not None, "kernel import spec unavailable")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "derive API missing")
    return module


def callback_rows_callback() -> tuple[Any, list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []

    def emit(*args: Any) -> None:
        require(len(args) == 2 and isinstance(args[0], str) and isinstance(args[1], Mapping), "emit must be exactly (event,row)")
        row = {"event": args[0], "row": dict(args[1])}
        append_journal({"schema": JOURNAL_SCHEMA, "event": "callback", **row})
        rows.append(row)

    return emit, rows


def count_prefix(source: Mapping[str, Any], rows: list[Mapping[str, Any]]) -> dict[str, int]:
    actual = {str(key): 0 for key in source["expected_counts"]}
    for row in rows:
        event = row.get("event")
        for target, events in source["count_event_map"].items():
            require(target in actual, f"unknown count target: {target}")
            if event in events:
                actual[target] += 1
    return actual


def false_flags(source: Mapping[str, Any]) -> dict[str, bool]:
    return {str(name): False for name in source["expected_new_flags"]}


def exact_nonfinite(value: Any) -> bool:
    if isinstance(value, str):
        return value in {"-inf", "+inf", "inf", "nan", "NaN", "Infinity", "-Infinity", "+Infinity"}
    if isinstance(value, Mapping):
        return any(exact_nonfinite(key) or exact_nonfinite(item) for key, item in value.items())
    if isinstance(value, (list, tuple)):
        return any(exact_nonfinite(item) for item in value)
    return False


def record_bundle(value: Any, label: str, *, allow_none: bool = False) -> None:
    if value is None and allow_none:
        return
    require(isinstance(value, Mapping), f"{label} bundle missing")
    require("actual" in value and "error" in value, f"{label} bundle lacks actual/error")
    candidate = value.get("candidate")
    if candidate is not None:
        require(isinstance(candidate, Mapping) or isinstance(candidate, (list, tuple)), f"{label} candidate malformed")


def record_interval(value: Any, label: str, *, allow_none: bool = False) -> None:
    if value is None and allow_none:
        return
    require(isinstance(value, Mapping) and "lower" in value and "upper" in value, f"{label} interval record missing")


def validate_payload(value: Mapping[str, Any]) -> None:
    payload = value.get("payload") if isinstance(value.get("payload"), Mapping) else value
    require(isinstance(payload, Mapping), "kernel payload missing")
    require(value.get("finite_rhs_packet_admitted") is True and value.get("finite_g5_uniform_C0_C1_admitted") is True, "kernel finite-admission flags missing")
    proof = payload.get("proof_packet")
    require(isinstance(proof, Mapping), "proof packet missing")
    require(proof.get("schema") == COMPONENT + "_proof_v1" and proof.get("proof_count") == 1 and proof.get("conditional_local_packet") is True, "proof packet identity changed")
    gap = proof.get("positive_denominator_gap")
    require(isinstance(gap, Mapping) and gap.get("sigma_squared_lower_strictly_positive") is True and gap.get("H_squared_nonnegative_by_sign_branch") is True and gap.get("D_lower_ge_sigma_lower") is True, "positive denominator gap proof not established")
    certificate = proof.get("angular_coefficient_certificate")
    require(isinstance(certificate, Mapping) and certificate.get("passed") is True and certificate.get("schema") == "finite_gauge_angular_coefficient_normal_form_v1" and certificate.get("C11_corrected") == "-4*d*z*rho/L^2" and certificate.get("C11_has_no_L_minus_3_term") is True, "normal-form certificate changed")
    definitions = payload.get("definitions")
    require(isinstance(definitions, Mapping) and isinstance(definitions.get("denominator"), Mapping) and definitions["denominator"].get("gap") == "lower(D)>=lower(sigma^2)>0", "denominator normal form missing")
    gauge = payload.get("gauge_table")
    require(isinstance(gauge, (list, Mapping)) and len(gauge) == 8, "gauge table must contain eight eta rows")
    rhs = payload.get("rhs_rows", payload.get("rhs_packet", {}).get("rows") if isinstance(payload.get("rhs_packet"), Mapping) else None)
    require(isinstance(rhs, list) and len(rhs) == 48, "rhs table must contain 48 rows")
    expected = {(eta, rho, quantity) for eta in PAYLOAD_ETA_LABELS for rho in PAYLOAD_RHO_LABELS for quantity in ("g4", "g5", "g6")}
    actual_selectors: set[tuple[str, str, str]] = set()
    uniform_g5_rows = 0
    new_g5_c0_rows = 0
    new_g5_c1_rows = 0
    for index, row in enumerate(rhs):
        require(isinstance(row, Mapping), f"rhs row {index} is not a mapping")
        selector = (str(row.get("eta")), str(row.get("rho")), str(row.get("quantity")))
        require(selector in expected, f"rhs selector changed at row {index}")
        actual_selectors.add(selector)
        require(row.get("atom") in ("one", "F0_anchor"), f"rhs atom changed at row {index}")
        for name in ("C0", "C1"):
            record_bundle(row.get(name), f"rhs {name} row {index}")
        for name in ("actual_C0", "actual_C1", "candidate_C0", "candidate_C1", "error_C0", "error_C1"):
            if name in row:
                record_interval(row.get(name), f"rhs {name} row {index}", allow_none=name.startswith("candidate_"))
        sectors = row.get("sectors")
        require(isinstance(sectors, (list, Mapping)) and len(sectors) > 0, f"rhs sectors missing at row {index}")
        sector_values = sectors.values() if isinstance(sectors, Mapping) else sectors
        for sector_index, sector in enumerate(sector_values):
            require(isinstance(sector, Mapping), f"rhs sector {index}/{sector_index} malformed")
            actual = sector.get("actual")
            require(isinstance(actual, (list, tuple)) and len(actual) == 2, f"rhs sector actual pair missing at {index}/{sector_index}")
            for endpoint_index, endpoint in enumerate(actual):
                record_interval(endpoint, f"rhs sector actual {index}/{sector_index}/{endpoint_index}")
            candidate = sector.get("candidate")
            if candidate is not None:
                require(isinstance(candidate, (list, tuple)) and len(candidate) == 2, f"rhs sector candidate pair malformed at {index}/{sector_index}")
                for endpoint_index, endpoint in enumerate(candidate):
                    record_interval(endpoint, f"rhs sector candidate {index}/{sector_index}/{endpoint_index}")
        if selector[2] == "g5":
            require(row.get("F0_anchor_formal", False) is True or row.get("atom") == "F0_anchor", f"g5 formal F0 status missing at row {index}")
            require(row.get("new_g5_C1") is True, f"g5 C1 repair flag missing at row {index}")
            new_g5_c1_rows += int(row.get("new_g5_C1") is True)
            new_g5_c0_rows += int(row.get("new_g5_C0") is True)
            if selector[0] == "uniform_eta":
                uniform_g5_rows += 1
    require(actual_selectors == expected, "rhs selectors are not exactly 8 eta by 2 rho by 3 quantities")
    summary = payload.get("summary")
    require(isinstance(summary, Mapping) and summary.get("eta_count") == 8 and summary.get("rho_count") == 2 and summary.get("rhs_row_count") == 48 and summary.get("copied_g4_g6_rows") == 32 and summary.get("copied_point_g5_C0_rows") == 14 and summary.get("new_g5_C0_rows") == 2 and summary.get("new_g5_C1_rows") == 16 and summary.get("new_ell2_uniform_rows") == 1 and summary.get("new_ell3_rows") == 8, "finite packet row-count summary changed")
    require(uniform_g5_rows == 2 and new_g5_c0_rows == 2 and new_g5_c1_rows == 16, "repaired g5 row partition changed")
    require(proof.get("finite_rhs_packet_admitted") is True and summary.get("finite_rhs_packet_admitted") is True and proof.get("regular_n1_response_ready") is False and summary.get("regular_n1_response_ready") is False, "finite admission or regular-response flags changed")
    require(proof.get("old_212_error_recertification") is False and proof.get("global_closure") is False and proof.get("no_old_provider_calls") is True and proof.get("no_checker_or_selftest") is True and proof.get("no_absolute_F0_exp") is True, "old/global boundary changed")
    require(not exact_nonfinite(payload), "saved finite payload contains an exact nonfinite scalar string")


def serializable(value: Any) -> Any:
    try:
        json.dumps(value, ensure_ascii=False, allow_nan=False)
        return value
    except (TypeError, ValueError):
        return repr(value)


_NO_RETURN = object()


def kernel_failure(source: Mapping[str, Any], exc: BaseException, rows: list[dict[str, Any]], imported: bool, entry_calls: int, calls: int, returned_value: Any = _NO_RETURN) -> dict[str, Any]:
    failure = {"kind": "manager_import_or_derive_failure", "failure_stage": "derive" if imported else "import", "error_type": type(exc).__name__, "error": repr(exc), "traceback": traceback.format_exc(), "kernel_imported": imported, "derive_entry_calls": entry_calls, "derive_calls": calls, "callback_prefix": list(rows), "callback_prefix_count": len(rows), "callback_prefix_counts": count_prefix(source, rows), "stopped": True}
    if returned_value is not _NO_RETURN:
        failure["returned_value_before_mapping_guard"] = serializable(returned_value)
    return {"schema": source["kernel_result_schema"], "component": COMPONENT, "passed": False, "failure": failure}


def wrap_result(source: Mapping[str, Any], kernel_result: Mapping[str, Any], imported: bool, entry_calls: int, calls: int, rows: list[dict[str, Any]], manager_failure: Mapping[str, Any] | None) -> dict[str, Any]:
    passed = kernel_result.get("passed") is True and manager_failure is None and kernel_result.get("finite_rhs_packet_admitted") is True and kernel_result.get("finite_g5_uniform_C0_C1_admitted") is True
    nested = dict(kernel_result)
    return {"schema": RESULT_SCHEMA, "component": COMPONENT, "attempt_status": "passed" if passed else "failed", "passed": passed, "metadata_only": False, "source_only": True, "numeric_evaluation": True, "numeric_values_computed": passed, "current_leaf_numeric_evaluation": True, "finite_rhs_packet_admitted": bool(passed and kernel_result.get("finite_rhs_packet_admitted") is True), "finite_g5_uniform_C0_C1_admitted": bool(passed and kernel_result.get("finite_g5_uniform_C0_C1_admitted") is True), "regular_n1_solution_view_installed": False, "regular_n1_binding_installed": False, "actual_n1_solution_materialized": False, "global_false_flags": {name: False for name in GLOBAL_FALSE_KEYS}, "result": nested, "actual_counts": dict(kernel_result.get("actual_counts", count_prefix(source, rows))) if passed else count_prefix(source, rows), "new_flags": dict(kernel_result.get("new_flags", false_flags(source))) if passed else false_flags(source), "proof_count": int(kernel_result.get("proof_count", 0)) if passed else 0, "derive_entry_calls": entry_calls, "derive_calls": calls, "kernel_imported": imported, "callback_events": list(rows), "callback_event_count": len(rows), "journal_event_count": len(journal_rows()), "callback_prefix_counts": count_prefix(source, rows), "failure": dict(manager_failure) if manager_failure else nested.get("failure"), "mathematical_proof_accepted": passed, "source_operator_implementation_accepted": passed, "math_source_identity_accepted": passed, "runtime_model_metadata": None, "source_checkpoint_sha256": sha_bytes(SOURCE_PATH.read_bytes())}


def validate_saved(source: Mapping[str, Any], saved: Mapping[str, Any], rows: list[dict[str, Any]], imported: bool, entry_calls: int, calls: int) -> dict[str, Any]:
    violations: list[str] = []
    if saved.get("schema") != RESULT_SCHEMA: violations.append("result schema")
    if saved.get("component") != COMPONENT: violations.append("result component")
    if saved.get("runtime_model_metadata") is not None: violations.append("runtime metadata")
    allowed = set(source["callback_event_counts"])
    if any(row.get("event") not in allowed for row in rows): violations.append("callback event name")
    callback_counts = {name: sum(1 for row in rows if row.get("event") == name) for name in source["callback_event_counts"]}
    if saved.get("passed") is True:
        if saved.get("actual_counts") != source["expected_counts"]: violations.append("actual counts")
        if saved.get("new_flags") != source["expected_new_flags"]: violations.append("new flags")
        if saved.get("proof_count") != source["expected_counts"].get("new_mathematical_proofs", 1): violations.append("proof count")
        if not (imported and entry_calls == 1 and calls == 1): violations.append("lifecycle")
        if callback_counts != source["callback_event_counts"]: violations.append("callback counts")
        if saved.get("journal_event_count") != len(journal_rows()): violations.append("journal count")
        try:
            validate_payload(saved["result"])
        except Exception as exc:
            violations.append(f"payload: {exc}")
    else:
        if saved.get("proof_count") != 0: violations.append("failed proof count")
        if saved.get("journal_event_count") != len(journal_rows()): violations.append("failed journal count")
        if len(journal_rows()) != 1 + int(imported) + len(rows): violations.append("failed journal prefix")
    return {"schema": f"{COMPONENT}_validation_v1", "component": COMPONENT, "passed": not violations and saved.get("passed") is True, "attempt_failed": saved.get("passed") is not True, "actual_counts": saved.get("actual_counts"), "new_flags": saved.get("new_flags"), "proof_count": saved.get("proof_count"), "callback_event_count": len(rows), "journal_event_count": len(journal_rows()), "violations": violations}


def run() -> dict[str, Any]:
    source_raw, source = load_source_for_run()
    require(not ENTERED_PATH.exists() and not JOURNAL_PATH.exists() and not RESULT_PATH.exists() and not VALIDATION_PATH.exists(), "run checkpoints already exist; reentry forbidden")
    entered = {"schema": ENTERED_SCHEMA, "component": COMPONENT, "source_checkpoint_sha256": sha_bytes(source_raw), "kernel": source["kernel"], "runner": source["runner"], "environment_metadata": dict(source.get("environment_metadata", {})), "import_attempts": 1, "derive_entry_calls": 0, "derive_calls": 0, "sys_modules_registration_before_exec": True, "metadata_only": False, "source_only": True, "numeric_evaluation": True, "numeric_values_computed": False, "finite_rhs_packet_admitted": False, "finite_g5_uniform_C0_C1_admitted": False, "regular_n1_solution_view_installed": False, "runtime_model_metadata": None, "entered_at_utc": now(), "do_not_reenter": True}
    entered_raw = write_new(ENTERED_PATH, entered)
    write_new(JOURNAL_PATH, b"")
    append_journal({"schema": JOURNAL_SCHEMA, "event": "manager_entry_before_import", "source_checkpoint_sha256": sha_bytes(source_raw), "entered_checkpoint_sha256": sha_bytes(entered_raw), "import_attempts": 1, "derive_entry_calls": 0, "derive_calls": 0})
    emit, rows = callback_rows_callback()
    imported = False; entry_calls = 0; calls = 0; manager_failure: Mapping[str, Any] | None = None; returned_value: Any = _NO_RETURN
    try:
        module = import_kernel(KERNEL_PATH)
        imported = True; entry_calls = 1
        append_journal({"schema": JOURNAL_SCHEMA, "event": "manager_derive_entry", "import_attempts": 1, "derive_entry_calls": 1, "derive_calls": 0, "sys_modules_registered": True})
        capsule_ref = source["source_wrappers"]["capsule"]
        _, _, capsule = read_json(capsule_ref["path"], capsule_ref["sha256"])
        capsule_status(capsule)
        data = dict(source["input"])
        data["capsule"] = capsule
        calls = 1
        returned_value = module.derive(data, emit)
        require(isinstance(returned_value, Mapping), "derive returned a non-mapping")
        kernel_result = dict(returned_value)
    except BaseException as exc:
        kernel_result = kernel_failure(source, exc, rows, imported, entry_calls, calls, returned_value)
        manager_failure = kernel_result["failure"]
    result_value = wrap_result(source, kernel_result, imported, entry_calls, calls, rows, manager_failure)
    result_raw = write_new(RESULT_PATH, result_value)
    _, saved_raw, saved = read_json(RESULT_PATH, sha_bytes(result_raw))
    validation = validate_saved(source, saved, rows, imported, entry_calls, calls)
    validation_receipt = {"schema": f"{COMPONENT}_validation_receipt_v1", "component": COMPONENT, "metadata_only": True, "new_math_calls": 0, "source_checkpoint_sha256": sha_bytes(source_raw), "result_sha256": sha_bytes(saved_raw), "result_bytes": len(saved_raw), "passed": bool(validation["passed"]), "attempt_failed": bool(validation["attempt_failed"]), "violations": list(validation["violations"]), "actual_counts": validation.get("actual_counts"), "new_flags": validation.get("new_flags"), "proof_count": validation.get("proof_count"), "callback_event_count": validation.get("callback_event_count"), "journal_event_count": validation.get("journal_event_count"), "kernel_finite_rhs_packet_admitted": saved.get("finite_rhs_packet_admitted") is True, "kernel_finite_g5_uniform_C0_C1_admitted": saved.get("finite_g5_uniform_C0_C1_admitted") is True, "finite_rhs_packet_admitted": bool(validation["passed"]), "finite_g5_uniform_C0_C1_admitted": bool(validation["passed"]), "regular_n1_response_ready": False, "runtime_model_metadata": None}
    validation_raw = write_new(VALIDATION_PATH, validation_receipt)
    return {"status": "run_complete", "component": COMPONENT, "attempt_status": saved.get("attempt_status"), "passed": bool(validation["passed"]), "metadata_only": False, "numeric_evaluation": True, "numeric_values_computed": bool(saved.get("passed")), "finite_rhs_packet_admitted": bool(validation["passed"]), "finite_g5_uniform_C0_C1_admitted": bool(validation["passed"]), "regular_n1_solution_view_installed": False, "actual_n1_solution_materialized": False, "environment_metadata": dict(source.get("environment_metadata", {})), "result": {"path": wire(RESULT_PATH), "sha256": sha_bytes(saved_raw), "bytes": len(saved_raw), "schema": saved.get("schema")}, "validation_receipt": {"path": wire(VALIDATION_PATH), "sha256": sha_bytes(validation_raw), "bytes": len(validation_raw), "schema": validation_receipt.get("schema")}, "actual_counts": saved.get("actual_counts"), "new_flags": saved.get("new_flags"), "proof_count": saved.get("proof_count"), "derive_entry_calls": entry_calls, "derive_calls": calls, "kernel_imported": imported, "callback_event_count": len(rows), "journal_event_count": len(journal_rows()), "validation": validation, "failure": saved.get("failure"), "runtime_model_metadata": None}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    prep.add_argument("--manifest", default=str(MANIFEST_PATH))
    sub.add_parser("run")
    args = parser.parse_args(argv)
    try:
        receipt = prepare(args.manifest) if args.command == "prepare" else run()
    except Exception as exc:
        print(json.dumps({"status": "error", "component": COMPONENT, "error": f"{type(exc).__name__}: {exc}", "metadata_only": False, "runtime_model_metadata": None}, sort_keys=True))
        return 2
    print(json.dumps(receipt, ensure_ascii=False, sort_keys=True))
    return 0 if receipt.get("passed", True) else 1


if __name__ == "__main__":
    raise SystemExit(main())
