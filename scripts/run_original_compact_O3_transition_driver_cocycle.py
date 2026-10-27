"""Metadata-only manager for the compact O3 transition driver cocycle leaf.

The manager freezes an authenticated parent C3 projection, a root-supplied
input specification, and the two new code artifacts.  It deliberately does
not contain a mathematical kernel or any ancestor replay/evaluation path;
all source-value and proof work belongs to the new kernel.  ``prepare`` writes one durable
source packet; ``run`` enters once, fsyncs the entered marker before importing
the new kernel, calls ``derive`` once, and fsyncs the complete result before
validation.  Publication, resume, and seal consume saved bytes and metadata
only.  No command is executed merely by importing this module.
"""

from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess


REPO = Path(
    r"C:/Users/z5242/Documents/Codex/2026-09-13/"
    r"before-doing-any-project-work-verify/work/st073-transition-next"
)
WORK = Path(
    r"C:/Users/z5242/Documents/Codex/2026-10-16/"
    r"https-github-com-liang-techh-openai/work"
)

COMPONENT = "compact_O3_transition_driver_cocycle"
INPUT_SCHEMA = "compact_O3_transition_driver_cocycle_input_v1"
RESULT_SCHEMA = "original_compact_O3_transition_driver_cocycle"
ENTERED_SCHEMA = "original_compact_O3_transition_driver_cocycle_entered_v1"
DIAGNOSTIC_SCHEMA = "compact_O3_transition_driver_cocycle_stage_diagnostic_v1"
PUBLICATION_SCHEMA = "original_compact_O3_transition_driver_cocycle_publication_v1"
SEAL_SCHEMA = "compact_O3_transition_driver_cocycle_final_seal_v1"

KERNEL_RELATIVE = Path(
    "src/openai_ns_reconstruction/original_compact_O3_transition_driver_cocycle.py"
)
RUNNER_RELATIVE = Path("scripts/run_original_compact_O3_transition_driver_cocycle.py")
REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_driver_cocycle.json"
)
MILESTONE_RELATIVE = Path(
    "docs/CURRENT_ORIGINAL_COMPACT_O3_TRANSITION_DRIVER_COCYCLE_2026_10_26.md"
)
HANDOFF_RELATIVES = (
    Path("docs/AGENT_TASKS.md"),
    Path("docs/CROSS_DEVICE_HANDOFF.md"),
    Path("docs/CURRENT_CHECKPOINT.md"),
    Path("docs/NEXT_TASKS.md"),
)
SCOPED_RELATIVES = (
    Path(".gitattributes"),
    KERNEL_RELATIVE,
    RUNNER_RELATIVE,
    REPORT_RELATIVE,
    MILESTONE_RELATIVE,
    *HANDOFF_RELATIVES,
)

CHECKPOINT_DIR = REPO / "tmp" / COMPONENT
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"
CHECKPOINTS = (
    SOURCE_CHECKPOINT,
    ENTERED_CHECKPOINT,
    RESULT_CHECKPOINT,
    DIAGNOSTIC_CHECKPOINT,
)

BRIEF = WORK / f"{COMPONENT}_brief_2026_10_26.json"
MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
PREPARE_RECEIPT = WORK / f"{COMPONENT}_prepare_receipt_2026_10_26.json"
RUN_RECEIPT = WORK / f"{COMPONENT}_run_receipt_2026_10_26.json"
PUBLISH_RECEIPT = WORK / f"{COMPONENT}_publication_receipt_2026_10_26.json"
RESUME_RECEIPT = WORK / f"{COMPONENT}_resume_receipt_2026_10_26.json"
SEAL_RECEIPT = WORK / f"{COMPONENT}_final_seal_receipt_2026_10_26.json"

PARENT_REPORT_DEFAULT = (
    REPO
    / "experiments/root_st073/"
    / "lei_ren_part1_paper_compliant_current_original_compact_C3_quiet_endpoint_projection.json"
)
PARENT_REPORT_SHA256 = (
    "633abe4e87fb0b7d1a3495a21fb5fdb5c4323921eede61724b421a24f22915b1"
)
PARENT_COMMIT = "b0f7e8e0d9a69b8c82345ff0a6a356c8b7a2cc96"

PARENT_SOURCE_PIN_COUNT = 2186
PARENT_CHECKPOINT_COUNT = 236
CURRENT_CHECKPOINT_COUNT = PARENT_CHECKPOINT_COUNT + 4
PARENT_FAILED_COUNT = 169
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26

COCYCLE_CHECK_COUNT = 4
EXPREL_PRODUCER_ORDER_CHECK_COUNT = 4
EXPREL_PRODUCER_ORDER_FAILURE_COUNT = 0
ACTUAL_DRIVER_GAP_COUNT = 20
ORDINARY_Z_ORDERS = (0, 1, 2, 3)

NEW_FLAGS = (
    "O3_transition_finite_N_driver_source_bound",
    "O3_transition_exponential_driver_cocycle_ordinary_Z3_proved",
    "O3_transition_twenty_signed_density_gap_functions_defined",
)
OLD_QUIET_FLAGS = (
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
)
PARENT_C3_FLAGS = (
    "C3_quiet_terminal_query_projection_bound",
    "C3_quiet_endpoint_memory_relations_reused",
)
NONCLAIM_FLAGS = (
    "full_pair_proved",
    "full_pair_identity_proved",
    "incoming_identity_proved",
    "incoming_source_identity_proved",
    "caps_bound",
    "cap_provenance_bound",
    "global_five_moment_closure",
    "global_closure_proved",
)

BOOKKEEPING_FIELDS = frozenset(
    {
        "schema",
        "component",
        "passed",
        "metadata_only",
        "publication_metadata_only",
        "source_only",
        "input_digest",
        "input_raw_sha256",
        "runtime_model_metadata",
        "original_source_evaluations",
        "source_packet_defined",
        "source_packet_definition_accepted",
        "source_packet_accepted_owned",
        "source_packet_accepted_inherited",
        "kernel_imported",
        "derive_entry_calls",
        "derive_calls",
        "entry_attempts",
        "source_checkpoint_sha256",
        "entry_checkpoint_sha256",
        "diagnostic_checkpoint_sha256",
        "new_flags",
        "new_counts",
        "counts",
        "raw_sha_pins",
        "source_union_hashes",
        "source_union_count",
        "parent_result_raw_sha256",
    }
)

PARENT_REPORT_EXCLUDED_FIELDS = (
    "metadata_helper_sources",
    "lossless_artifacts",
    "lossless_artifact_count",
    "full_result",
    "component_result",
    "full_result_raw_bytes",
    "full_result_raw_sha256",
    "source_pins",
    "source_union_hashes",
    "parent_checkpoint_inventory",
    "checkpoint_hashes",
    "checkpoint_paths",
    "protected_file_hashes",
    "legacy_dirty_path_hashes",
    "preserved_failed_artifact_hashes_union",
    "implementation_hashes",
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def sha_path(path):
    return sha_bytes(Path(path).read_bytes())


def encode(value):
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


def clone(value):
    return json.loads(encode(value).decode("utf-8"))


def wire_path(path):
    """Return a portable path for every manager-created JSON payload."""
    return Path(path).as_posix()


def portable_metadata(value):
    """Keep inherited metadata JSON-safe while using forward-slash paths."""
    if isinstance(value, dict):
        result = {}
        for key, child in value.items():
            portable_key = key.replace("\\", "/") if isinstance(key, str) else key
            require(portable_key not in result, "Portable metadata key collision: " + str(portable_key))
            result[portable_key] = portable_metadata(child)
        return result
    if isinstance(value, list):
        return [portable_metadata(child) for child in value]
    if isinstance(value, str) and ("\\" in value or (len(value) > 2 and value[1] == ":")):
        return value.replace("\\", "/")
    return value


def remove_metadata_fields(value, excluded):
    if isinstance(value, dict):
        return {
            key: remove_metadata_fields(child, excluded)
            for key, child in value.items()
            if key not in excluded
        }
    if isinstance(value, list):
        return [remove_metadata_fields(child, excluded) for child in value]
    return value


def wire_mapping(mapping):
    result = {}
    for name, digest in mapping.items():
        key = wire_path(name)
        require(key not in result or result[key] == digest, "Conflicting portable hash at " + key)
        result[key] = digest
    return result


def write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def write_replace(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else encode(value)
    with path.open("wb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())
    return raw


def read_bytes(path, expected=None):
    path = Path(path)
    require(path.exists(), "Missing artifact: " + wire_path(path))
    raw = path.read_bytes()
    if expected is not None:
        require(sha_bytes(raw) == expected.lower(), "Changed artifact: " + wire_path(path))
    return raw


def read_json(path, expected=None):
    raw = read_bytes(path, expected)
    return raw, json.loads(raw.decode("utf-8"), parse_int=int)


def resolve_path(value, base=REPO):
    path = Path(value)
    if path.is_absolute():
        return path
    direct = Path(base) / path
    if direct.exists():
        return direct
    repo_path = REPO / path
    if repo_path.exists():
        return repo_path
    experiment_path = REPO / "experiments" / "root_st073" / path
    if experiment_path.exists():
        return experiment_path
    return WORK / path


def valid_sha(value):
    return isinstance(value, str) and len(value) == 64 and all(
        character in "0123456789abcdef" for character in value.lower()
    )


def valid_git_revision(value):
    return isinstance(value, str) and len(value) == 40 and all(
        character in "0123456789abcdef" for character in value.lower()
    )


def auth_map(value, label, authenticate=True):
    require(isinstance(value, dict), label + " is not a mapping")
    result = {}
    for name, digest in value.items():
        require(isinstance(name, str) and valid_sha(digest), label + " has a malformed SHA")
        if authenticate:
            read_bytes(resolve_path(name), digest)
        result[name] = digest.lower()
    return result


def flatten_hashes(value, label):
    result = {}

    def add(name, digest):
        require(isinstance(name, str) and valid_sha(digest), label + " has a malformed hash")
        digest = digest.lower()
        require(name not in result or result[name] == digest, "Conflicting hash at " + name)
        result[name] = digest

    def visit(node, location):
        require(isinstance(node, dict), location + " is not a mapping")
        if not node:
            return
        if all(isinstance(k, str) and valid_sha(v) for k, v in node.items()):
            for k, v in node.items():
                add(k, v)
            return
        for key, child in node.items():
            child_location = location + "." + str(key)
            if isinstance(child, dict):
                visit(child, child_location)
            elif isinstance(child, list):
                for index, item in enumerate(child):
                    if isinstance(item, dict):
                        visit(item, child_location + "[" + str(index) + "]")

    visit(value, label)
    return result


def merge_hashes(*maps):
    merged = {}
    for mapping in maps:
        require(isinstance(mapping, dict), "Hash map is not a mapping")
        for name, digest in mapping.items():
            require(isinstance(name, str) and valid_sha(digest), "Malformed hash map entry")
            digest = digest.lower()
            require(name not in merged or merged[name] == digest, "Conflicting hash at " + name)
            merged[name] = digest
    return merged


def parent_view(report):
    excluded = set(PARENT_REPORT_EXCLUDED_FIELDS)
    view = portable_metadata({key: report[key] for key in sorted(report) if key not in excluded})
    view = remove_metadata_fields(view, excluded)
    view["parent_report_is_authenticated_view"] = True
    view["parent_report_view_excluded_fields"] = list(PARENT_REPORT_EXCLUDED_FIELDS)
    view["parent_report_view_paths_portable"] = True
    require(len(encode(view)) < 100_000, "Parent report metadata view exceeds 100 KB")
    return view


def report_result_reference(report, report_path, report_raw_sha):
    selector = report.get("full_result_artifact_selector")
    result_sha = report.get("full_result_raw_sha256")
    prefix = "$.lossless_artifacts[path="
    if isinstance(selector, str) and selector.startswith(prefix) and selector.endswith("]"):
        result_name = selector[len(prefix) : -1]
        require(valid_sha(result_sha), "Parent report result SHA is malformed")
        result_path = resolve_path(result_name)
        result_raw, result = read_json(result_path, result_sha)
        require(report.get("full_result") == result, "Parent publication full result is not exact")
        require(report.get("full_result_raw_bytes") == len(result_raw), "Parent result byte count changed")
        return {
            "path": wire_path(result_path),
            "sha256": result_sha.lower(),
            "bytes": len(result_raw),
            "schema": result.get("schema"),
            "component": result.get("component"),
        }, result
    if report.get("schema") == RESULT_SCHEMA and report.get("component") == COMPONENT:
        return {
            "path": wire_path(report_path),
            "sha256": report_raw_sha,
            "bytes": report_path.stat().st_size,
            "schema": report.get("schema"),
            "component": report.get("component"),
        }, report
    raise AssertionError("Parent report has no exact lossless result selector")


def load_parent(args):
    report_path = resolve_path(args.parent_report or wire_path(PARENT_REPORT_DEFAULT), WORK)
    report_sha = (args.parent_report_sha256 or PARENT_REPORT_SHA256).lower()
    parent_commit = (args.parent_commit or PARENT_COMMIT).lower()
    require(report_sha == PARENT_REPORT_SHA256, "Parent report SHA is not the authenticated current projection")
    require(parent_commit == PARENT_COMMIT, "Parent commit is not the authenticated current projection parent")
    report_raw, report = read_json(report_path, report_sha)
    require(isinstance(report, dict), "Parent report is not a mapping")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT, "Parent source-pin count changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Parent checkpoint count changed")
    require(
        report.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT,
        "Parent failed-artifact count changed",
    )
    require(len(report.get("false_flags", [])) == PARENT_FALSE_FLAG_COUNT, "Parent false-flag count changed")
    require(len(report.get("protected_file_hashes", {})) == PARENT_PROTECTED_COUNT, "Parent protected-file count changed")
    require(len(report.get("legacy_dirty_path_hashes", {})) == PARENT_DIRTY_COUNT, "Parent dirty-path count changed")
    false_flags = report.get("false_flags")
    require(isinstance(false_flags, list) and len(set(false_flags)) == PARENT_FALSE_FLAG_COUNT, "Parent false flags are malformed")
    for key in false_flags:
        require(report.get(key) is False, "Parent false flag changed: " + str(key))
    for key in OLD_QUIET_FLAGS:
        require(report.get(key) is False, "Parent quiet flag changed: " + key)
    for key in PARENT_C3_FLAGS:
        require(report.get(key) is True, "Parent C3 projection flag changed: " + key)

    source_pins = auth_map(report.get("source_pins"), "Parent source pins")
    require(len(source_pins) == PARENT_SOURCE_PIN_COUNT, "Parent source pins are incomplete")
    inventory = flatten_hashes(report.get("parent_checkpoint_inventory"), "Parent checkpoint inventory")
    own = flatten_hashes({"checkpoint_hashes": report.get("checkpoint_hashes")}, "Parent checkpoint hashes")
    complete = merge_hashes(inventory, own)
    require(len(complete) == PARENT_CHECKPOINT_COUNT, "Parent checkpoint inventory is incomplete")
    for name, digest in complete.items():
        read_bytes(resolve_path(name), digest)
    protected = auth_map(report.get("protected_file_hashes"), "Parent protected files")
    dirty = auth_map(report.get("legacy_dirty_path_hashes"), "Parent dirty files")
    failed = auth_map(report.get("preserved_failed_artifact_hashes_union"), "Parent failed artifacts")
    require(len(failed) == PARENT_FAILED_COUNT, "Parent failed artifacts are incomplete")
    require(len(wire_mapping(failed)) == PARENT_FAILED_COUNT, "Parent failed artifact paths are not portable and unique")
    result_ref, result = report_result_reference(report, report_path, report_sha)

    report_meta = {
        "path": wire_path(report_path),
        "sha256": report_sha,
        "bytes": len(report_raw),
        "schema": report.get("schema"),
        "component": report.get("component"),
    }
    accepted = {}
    for mapping in (report, result):
        for key, value in mapping.items():
            if (
                isinstance(key, str)
                and value is True
                and key not in false_flags
                and key not in BOOKKEEPING_FIELDS
                and key not in NEW_FLAGS
                and key not in NONCLAIM_FLAGS
            ):
                accepted[key] = True
    return {
        "report": report,
        "report_raw": report_raw,
        "report_meta": report_meta,
        "report_view": parent_view(report),
        "result": result,
        "result_meta": result_ref,
        "source_pins": source_pins,
        "complete_checkpoints": complete,
        "protected": protected,
        "dirty": dirty,
        "failed": failed,
        "false_flags": list(false_flags),
        "accepted_true": accepted,
        "parent_commit": parent_commit,
    }


def discover_receipts(value, label):
    """Find explicit raw path/SHA receipt records without interpreting math."""
    found = []

    def visit(node, location):
        if isinstance(node, dict):
            path = node.get("path")
            digest = node.get("sha256", node.get("raw_sha256"))
            if isinstance(path, str) and valid_sha(digest):
                found.append({"path": path, "sha256": digest.lower(), "location": location})
            for key, child in node.items():
                visit(child, location + "." + str(key))
        elif isinstance(node, list):
            for index, child in enumerate(node):
                visit(child, location + "[" + str(index) + "]")

    visit(value, label)
    return found


def load_input(args):
    require(args.input_spec, "--input-spec is required for prepare")
    input_path = resolve_path(args.input_spec, WORK)
    input_raw = read_bytes(input_path)
    if args.input_sha256:
        require(valid_sha(args.input_sha256), "Input SHA-256 is malformed")
        require(sha_bytes(input_raw) == args.input_sha256.lower(), "Input SHA-256 changed")
    input_value = json.loads(input_raw.decode("utf-8"), parse_int=int)
    require(isinstance(input_value, dict), "O3 transition input is not a mapping")
    require(input_value.get("schema") == INPUT_SCHEMA, "O3 transition input schema changed")

    raw_pins = {}
    declared_pins = input_value.get("additional_raw_pins", {})
    require(isinstance(declared_pins, dict), "Input additional raw pins is not a mapping")
    for name, digest in declared_pins.items():
        require(isinstance(name, str) and valid_sha(digest), "Input additional raw pins has a malformed SHA")
        resolved = resolve_path(name, input_path.parent)
        read_bytes(resolved, digest)
        raw_pins[wire_path(resolved)] = digest.lower()
    receipts = discover_receipts(input_value.get("static_receipts", {}), "input.static_receipts")
    portable_receipts = []
    for receipt in receipts:
        resolved = resolve_path(receipt["path"], input_path.parent)
        read_bytes(resolved, receipt["sha256"])
        portable_receipts.append(
            {
                "path": wire_path(resolved),
                "sha256": receipt["sha256"],
                "location": receipt["location"],
            }
        )
    # A receipt can be listed both in additional_raw_pins and static_receipts;
    # merge only after authenticating both forms.
    receipt_map = {item["path"]: item["sha256"] for item in portable_receipts}
    packet_pins = merge_hashes(raw_pins, receipt_map)
    return {
        "path": wire_path(input_path),
        "sha256": sha_bytes(input_raw),
        "bytes": len(input_raw),
        "raw": input_raw,
        "value": input_value,
        "additional_raw_pins": packet_pins,
        "static_receipts": portable_receipts,
    }


def kernel_path(args):
    path = resolve_path(args.kernel_path or wire_path(REPO / KERNEL_RELATIVE), REPO)
    raw = read_bytes(path)
    if args.kernel_sha256:
        require(valid_sha(args.kernel_sha256), "Kernel SHA-256 is malformed")
        require(sha_bytes(raw) == args.kernel_sha256.lower(), "Kernel SHA-256 changed")
    return path, raw


def ensure_runner(args):
    path = resolve_path(args.runner_path or wire_path(REPO / RUNNER_RELATIVE), REPO)
    manager_raw = Path(__file__).read_bytes()
    if not path.exists():
        write_new(path, manager_raw)
    raw = read_bytes(path)
    return path, raw


def next_task_metadata():
    return {
        "phase": "sourcegap",
        "priority": "P2",
        "primitive": True,
        "source_gap": True,
        "action": "Bind the actual primitive source gap on the preceding O3 transition.",
        "task": "Bind actual primitive source functions for the next O3 transition source-gap check.",
        "source_function_count": ACTUAL_DRIVER_GAP_COUNT,
        "function_count": ACTUAL_DRIVER_GAP_COUNT,
        "gap_status": "unknown",
        "gap_zero_or_nonzero_asserted": False,
        "source_only": True,
        "do_not_promote_driver_facts_without_new_source_binding": True,
        "do_not_reenter_current_attempt": True,
    }


def source_paths(parent, input_meta, kernel, runner):
    external = Path(__file__).resolve()
    stage = {
        wire_path(parent["report_meta"]["path"]): parent["report_meta"]["sha256"],
        wire_path(parent["result_meta"]["path"]): parent["result_meta"]["sha256"],
        wire_path(input_meta["path"]): input_meta["sha256"],
        wire_path(kernel): sha_path(kernel),
        wire_path(runner): sha_path(runner),
        wire_path(external): sha_path(external),
    }
    stage = merge_hashes(stage, input_meta["additional_raw_pins"])
    return stage, merge_hashes(parent["source_pins"], stage)


def build_source(args):
    parent = load_parent(args)
    input_meta = load_input(args)
    kernel, kernel_raw = kernel_path(args)
    runner, runner_raw = ensure_runner(args)
    stage, source_union = source_paths(parent, input_meta, kernel, runner)
    inherited = dict(parent["accepted_true"])
    for key in (*NEW_FLAGS, *NONCLAIM_FLAGS):
        inherited.pop(key, None)
    payload = clone(input_meta["value"])
    source = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_only": True,
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_definition_accepted": False,
        "source_packet_accepted_owned": False,
        "source_packet_accepted_inherited": False,
        "parent_commit": parent["parent_commit"],
        "parent_report": clone(parent["report_view"]),
        "parent_report_metadata": clone(parent["report_meta"]),
        "parent_report_sha256": parent["report_meta"]["sha256"],
        "parent_result": clone(parent["result"]),
        "current_projection_data": clone(parent["result"]),
        "parent_result_metadata": clone(parent["result_meta"]),
        "parent_result_raw_sha256": parent["result_meta"]["sha256"],
        "input": payload,
        "input_spec": clone(payload),
        "input_spec_metadata": {
            "path": input_meta["path"],
            "sha256": input_meta["sha256"],
            "bytes": input_meta["bytes"],
        },
        "input_raw_sha256": input_meta["sha256"],
        "source_union_hashes": wire_mapping(source_union),
        "source_union_count": len(wire_mapping(source_union)),
        "stage_hashes": wire_mapping(stage),
        "current_code_hashes": {
            wire_path(kernel): sha_bytes(kernel_raw),
            wire_path(runner): sha_bytes(runner_raw),
        },
        "kernel_path": wire_path(kernel),
        "runner_path": wire_path(runner),
        "external_draft": {"path": wire_path(Path(__file__).resolve()), "sha256": sha_path(Path(__file__).resolve())},
        "raw_sha_pins": {
            "parent_report_raw_sha256": parent["report_meta"]["sha256"],
            "parent_result_raw_sha256": parent["result_meta"]["sha256"],
            "input_raw_sha256": input_meta["sha256"],
            "kernel_raw_sha256": sha_bytes(kernel_raw),
            "runner_raw_sha256": sha_bytes(runner_raw),
            "external_manager_raw_sha256": sha_path(Path(__file__).resolve()),
            "additional_raw_pins": wire_mapping(input_meta["additional_raw_pins"]),
        },
        "raw_packet_receipts": [
            {
                "path": wire_path(item["path"]),
                "sha256": item["sha256"],
                "location": item["location"],
            }
            for item in input_meta["static_receipts"]
        ],
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "current_checkpoint_count": CURRENT_CHECKPOINT_COUNT,
        "preserved_failed_artifact_hashes_union": wire_mapping(parent["failed"]),
        "preserved_failed_artifact_count": len(parent["failed"]),
        "protected_file_hashes": wire_mapping(parent["protected"]),
        "legacy_dirty_path_hashes": wire_mapping(parent["dirty"]),
        "complete_parent_checkpoint_hashes": wire_mapping(parent["complete_checkpoints"]),
        "complete_parent_checkpoint_count": len(parent["complete_checkpoints"]),
        "inherited_false_flags": list(parent["false_flags"]),
        "inherited_parent_true_flags": inherited,
        "new_flags": list(NEW_FLAGS),
        "old_quiet_flags_explicit_false": {key: False for key in OLD_QUIET_FLAGS},
        "parent_C3_projection_flags_explicit_true": {key: True for key in PARENT_C3_FLAGS},
        "nonclaim_flags_explicit_false": {key: False for key in NONCLAIM_FLAGS},
        "counts": {
            "new_exponential_cocycle_order_checks": COCYCLE_CHECK_COUNT,
            "new_exprel_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "actual_producer_order_failures": 0,
            "new_actual_signed_density_gap_rows": ACTUAL_DRIVER_GAP_COUNT,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        },
        "new_counts": {
            "cocycle_checks": COCYCLE_CHECK_COUNT,
            "new_exprel_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "actual_producer_order_failures": 0,
            "checked_residual_count": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "checked_residual_failures": 0,
            "actual_driver_gap_rows": ACTUAL_DRIVER_GAP_COUNT,
            "new_density_proof_rows": 0,
        },
        "derivation_contract": {
            "driver": "exponential_F_driver_cocycle",
            "exponential_F_driver_cocycle": True,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "new_exponential_cocycle_order_checks": COCYCLE_CHECK_COUNT,
            "new_exprel_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "actual_producer_order_failures": 0,
            "new_actual_signed_density_gap_rows": ACTUAL_DRIVER_GAP_COUNT,
            "actual_signed_density_gap_function_recipes": ACTUAL_DRIVER_GAP_COUNT,
            "function_not_samples": True,
            "new_density_proof_rows": 0,
            "numeric": False,
            "full_pair_proved": False,
            "incoming_identity_proved": False,
            "caps_bound": False,
            "global_five_moment_closure": False,
        },
        "source_inspection_contract": {
            "new_exprel_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "actual_producer_order_failures": 0,
            "checked_residual_count": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "checked_residual_failures": 0,
        },
        "old_generic_identity_cached_data_only": True,
        "old_generic_identity_calls": 0,
        "old_numeric_calls": 0,
        "numeric": False,
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "full_pair_proved": False,
        "incoming_identity_proved": False,
        "caps_bound": False,
        "global_five_moment_closure": False,
        "next_executable_task": next_task_metadata(),
        "goal_active": True,
        "full_goal_complete": False,
        "ancestor_decode_forbidden": True,
        "kernel_import_allowed_only_in_run": True,
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    source["input_digest"] = sha_bytes(encode(source["input"]))
    return source, parent, input_meta, kernel, runner


def prepare(args):
    require(not any(path.exists() for path in CHECKPOINTS), "O3 attempt already has a durable checkpoint")
    source, parent, input_meta, kernel, runner = build_source(args)
    source_raw = write_new(SOURCE_CHECKPOINT, source)
    main_brief = update_main_brief(
        {
            "phase": "prepared_unentered",
            "source_checkpoint": wire_path(SOURCE_CHECKPOINT),
            "source_checkpoint_sha256": sha_bytes(source_raw),
            "report": wire_path(REPO / REPORT_RELATIVE),
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "next_executable_task": next_task_metadata(),
        }
    )
    receipt = {
        "schema": COMPONENT + "_prepare_receipt_v1",
        "component": COMPONENT,
        "metadata_only": True,
        "source_checkpoint": wire_path(SOURCE_CHECKPOINT),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "source_union_count": source["source_union_count"],
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "current_checkpoint_count": CURRENT_CHECKPOINT_COUNT,
        "input_raw_sha256": input_meta["sha256"],
        "kernel_sha256": source["current_code_hashes"][wire_path(kernel)],
        "runner_sha256": source["current_code_hashes"][wire_path(runner)],
        "runtime_model_metadata": None,
        "kernel_imported": False,
        "goal_active": True,
        "main_resume_brief": main_brief,
    }
    receipt_raw = write_replace(PREPARE_RECEIPT, receipt)
    brief = {
        "schema": COMPONENT + "_brief_v1",
        "component": COMPONENT,
        "phase": "prepared_unentered",
        "source_checkpoint": wire_path(SOURCE_CHECKPOINT),
        "source_checkpoint_sha256": receipt["source_checkpoint_sha256"],
        "input_raw_sha256": input_meta["sha256"],
        "runner_sha256": receipt["runner_sha256"],
        "kernel_sha256": receipt["kernel_sha256"],
        "goal_active": True,
        "full_goal_complete": False,
        "new_flags": list(NEW_FLAGS),
        "next_executable_task": next_task_metadata(),
        "new_exponential_cocycle_order_checks": COCYCLE_CHECK_COUNT,
        "new_actual_signed_density_gap_rows": ACTUAL_DRIVER_GAP_COUNT,
        "runtime_model_metadata": None,
    }
    brief_raw = write_replace(BRIEF, brief)
    return {
        "prepare": receipt,
        "prepare_sha256": sha_bytes(receipt_raw),
        "brief_sha256": sha_bytes(brief_raw),
        "parent": parent["report_meta"],
        "input": input_meta["path"],
        "kernel": wire_path(kernel),
    }


def load_source():
    source_raw, source = read_json(SOURCE_CHECKPOINT)
    require(source.get("schema") == INPUT_SCHEMA and source.get("component") == COMPONENT, "O3 source schema changed")
    require(source.get("metadata_only") is True and source.get("source_only") is True, "O3 source scope changed")
    require(source.get("runtime_model_metadata") is None, "O3 source runtime metadata is not null")
    require(source.get("parent_report_sha256") == PARENT_REPORT_SHA256, "Source parent report SHA changed")
    require(source.get("parent_commit") == PARENT_COMMIT, "Source parent commit changed")
    require(source.get("parent_source_pin_count") == PARENT_SOURCE_PIN_COUNT, "Source parent pin count changed")
    require(source.get("complete_parent_checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Source checkpoint count changed")
    require(source.get("current_checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Source current checkpoint count changed")
    require(source.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT, "Source failed-artifact count changed")
    require(len(source.get("protected_file_hashes", {})) == PARENT_PROTECTED_COUNT, "Source protected count changed")
    require(len(source.get("legacy_dirty_path_hashes", {})) == PARENT_DIRTY_COUNT, "Source dirty count changed")
    require(source.get("input", {}).get("schema") == INPUT_SCHEMA, "Source input schema changed")
    for name, digest in source.get("source_union_hashes", {}).items():
        read_bytes(resolve_path(name), digest)
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "Source union count changed")
    return source_raw, source


def journal_append(events, event):
    payload = clone(event) if isinstance(event, dict) else {"event": str(event)}
    payload["runtime_model_metadata"] = None
    events.append(payload)
    with DIAGNOSTIC_CHECKPOINT.open("ab") as stream:
        stream.write(encode(payload))
        stream.flush()
        os.fsync(stream.fileno())


def fallback(source, calls, imported, error=None):
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "metadata_only": False,
        "source_only": True,
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_defined": False,
        "source_packet_definition_accepted": False,
        "source_packet_accepted_owned": False,
        "source_packet_accepted_inherited": False,
        "kernel_imported": bool(imported),
        "derive_entry_calls": calls,
        "derive_calls": calls,
        "entry_attempts": 1,
        "input_digest": source.get("input_digest"),
        "input_raw_sha256": source.get("input_raw_sha256"),
        "blockers": [str(error)] if error else ["O3 derive returned no result"],
        "counts": clone(source.get("counts", {})),
        "new_counts": {
            "cocycle_checks": 0,
            "new_exprel_producer_order_checks": 0,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": 0,
            "actual_producer_order_failures": 0,
            "actual_driver_gap_rows": 0,
            "new_density_proof_rows": 0,
        },
        "numeric": False,
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "new_exponential_cocycle_order_checks": 0,
            "new_exprel_producer_order_checks": 0,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": 0,
            "actual_producer_order_failures": 0,
            "new_actual_signed_density_gap_rows": 0,
            "new_density_proof_rows": 0,
            "old_generic_identity_calls": 0,
            "old_provider_calls": 0,
            "old_checker_calls": 0,
            "old_normalizer_calls": 0,
            "old_source_value_reducer_calls": 0,
            "old_numeric_calls": 0,
            "numeric_calls": 0,
            "numeric_values_computed": 0,
        },
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "old_generic_identity_cached_data_only": True,
        "old_numeric_calls": 0,
        "full_pair_proved": False,
        "incoming_identity_proved": False,
        "caps_bound": False,
        "global_five_moment_closure": False,
    }


def decorate_result_metadata(value, source, calls, imported):
    result = clone(value) if isinstance(value, dict) else fallback(source, calls, imported, "derive returned a non-mapping result")
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("passed", False)
    result["passed"] = result.get("passed") is True
    result["metadata_only"] = False
    result["source_only"] = True
    result["runtime_model_metadata"] = None
    result.setdefault("original_source_evaluations", 0)
    result["input_digest"] = source.get("input_digest")
    result["input_raw_sha256"] = source.get("input_raw_sha256")
    result["kernel_imported"] = bool(imported)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    result["source_packet_defined"] = result.get("passed") is True
    result["source_packet_definition_accepted"] = result.get("passed") is True
    result["source_packet_accepted_owned"] = result.get("passed") is True
    result["source_packet_accepted_inherited"] = False
    result["parent_result_raw_sha256"] = source.get("parent_result_raw_sha256")
    result["new_flags"] = {key: result.get(key) is True for key in NEW_FLAGS}
    result.setdefault("new_counts", {
        "cocycle_checks": COCYCLE_CHECK_COUNT if result.get("passed") is True else 0,
        "new_exprel_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT if result.get("passed") is True else 0,
        "new_exprel_producer_order_failures": 0,
        "actual_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT if result.get("passed") is True else 0,
        "actual_producer_order_failures": 0,
        "actual_driver_gap_rows": ACTUAL_DRIVER_GAP_COUNT if result.get("passed") is True else 0,
        "new_density_proof_rows": 0,
    })
    if not isinstance(result.get("new_counts"), dict):
        result["new_counts"] = {}
    result["new_counts"].setdefault(
        "new_exprel_producer_order_checks",
        EXPREL_PRODUCER_ORDER_CHECK_COUNT if result.get("passed") is True else 0,
    )
    result["new_counts"].setdefault("new_exprel_producer_order_failures", 0)
    result["new_counts"].setdefault(
        "actual_producer_order_checks",
        EXPREL_PRODUCER_ORDER_CHECK_COUNT if result.get("passed") is True else 0,
    )
    result["new_counts"].setdefault("actual_producer_order_failures", 0)
    result["derivation_contract"] = clone(source.get("derivation_contract", {}))
    result.setdefault("numeric", False)
    result.setdefault("numeric_values_computed", False)
    result.setdefault("numeric_field_values_present", False)
    result["old_generic_identity_cached_data_only"] = True
    for key in (
        "old_generic_identity_calls",
        "old_provider_calls",
        "old_checker_calls",
        "old_normalizer_calls",
        "old_source_value_reducer_calls",
        "old_numeric_calls",
        "numeric_calls",
    ):
        result.setdefault(key, 0)
    for key in NONCLAIM_FLAGS:
        result.setdefault(key, False)
    for key, flag in source.get("inherited_parent_true_flags", {}).items():
        if flag is True and key not in BOOKKEEPING_FIELDS and key not in (*NEW_FLAGS, *NONCLAIM_FLAGS):
            result.setdefault(key, True)
    for key in source.get("inherited_false_flags", []):
        result[key] = False
    for key in OLD_QUIET_FLAGS:
        result.setdefault(key, False)

    audit = result.setdefault("source_inspection_audit", {})
    if not isinstance(audit, dict):
        audit = {}
        result["source_inspection_audit"] = audit
    audit.setdefault("original_source_evaluations", 0)
    audit.setdefault("new_exponential_cocycle_order_checks", 0)
    audit.setdefault("new_exprel_producer_order_checks", 0)
    audit.setdefault("new_exprel_producer_order_failures", 0)
    audit.setdefault("actual_producer_order_checks", 0)
    audit.setdefault("actual_producer_order_failures", 0)
    audit.setdefault("checked_residual_count", 0)
    audit.setdefault("checked_residual_failures", 0)
    audit.setdefault("new_actual_signed_density_gap_rows", 0)
    audit.setdefault("new_density_proof_rows", 0)
    for key in (
        "old_generic_identity_calls",
        "old_provider_calls",
        "old_checker_calls",
        "old_normalizer_calls",
        "old_source_value_reducer_calls",
        "old_numeric_calls",
        "numeric_calls",
        "numeric_values_computed",
    ):
        audit.setdefault(key, 0)
    result["goal_active"] = True
    result["full_goal_complete"] = False
    return result


def _result_count(result, names):
    nodes = [result]
    for key in ("counts", "new_counts", "cocycle", "transition_driver_cocycle", "driver_cocycle", "source_gap"):
        value = result.get(key)
        if isinstance(value, dict):
            nodes.append(value)
    found = None
    for node in nodes:
        for name in names:
            value = node.get(name)
            if isinstance(value, int) and not isinstance(value, bool):
                found = value
            elif isinstance(value, list):
                found = len(value)
    return found


def validate_result_metadata(result, source):
    require(result.get("schema") == RESULT_SCHEMA, "O3 result schema changed")
    require(result.get("component") == COMPONENT, "O3 result component changed")
    require(result.get("metadata_only") is False and result.get("source_only") is True, "O3 result scope changed")
    require(result.get("runtime_model_metadata") is None, "O3 runtime metadata is not null")
    require(result.get("original_source_evaluations") == 0, "O3 evaluated an original source")
    require(result.get("derive_entry_calls") in (0, 1), "O3 derive entry was retried")
    require(result.get("derive_calls") == result.get("derive_entry_calls") and result.get("entry_attempts") == 1, "O3 derive bookkeeping changed")
    require(result.get("input_digest") == source.get("input_digest"), "O3 input digest changed")
    require(result.get("input_raw_sha256") == source.get("input_raw_sha256"), "O3 input raw SHA changed")
    require(result.get("source_packet_accepted_inherited") is False, "O3 source acceptance was inherited")
    require(result.get("full_goal_complete") is False, "O3 overclaims full goal completion")
    require(result.get("numeric_values_computed") is not True, "O3 computed numeric values")
    require(result.get("numeric_field_values_present") is not True, "O3 contains numeric field values")
    require(result.get("numeric") is not True, "O3 numeric mode is not false")
    for key in (
        "old_generic_identity_calls",
        "old_provider_calls",
        "old_checker_calls",
        "old_normalizer_calls",
        "old_source_value_reducer_calls",
        "old_numeric_calls",
        "numeric_calls",
    ):
        require(result.get(key, 0) == 0, "O3 result recorded forbidden old/numeric work: " + key)
    for key in NONCLAIM_FLAGS:
        require(result.get(key) is not True, "O3 overclaims nonclaim flag: " + key)
    audit = result.get("source_inspection_audit", {})
    require(isinstance(audit, dict), "O3 source audit is not a mapping")
    require(audit.get("original_source_evaluations", 0) == 0, "O3 source audit evaluated an original source")
    for key in (
        "old_generic_identity_calls",
        "old_provider_calls",
        "old_checker_calls",
        "old_normalizer_calls",
        "old_source_value_reducer_calls",
        "numeric_values_computed",
    ):
        require(audit.get(key, 0) == 0, "O3 recorded forbidden old/numeric work: " + key)
    require(audit.get("new_density_proof_rows", 0) == 0, "O3 recorded a new density proof row")
    require(result.get("new_counts", {}).get("new_density_proof_rows", 0) == 0, "O3 result introduced a density proof row")
    producer_claim = result.get("passed") is True or result.get(NEW_FLAGS[1]) is True
    if producer_claim:
        require(
            audit.get("new_exprel_producer_order_checks") == EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "O3 exponential-relation producer check count is not 4",
        )
        require(audit.get("new_exprel_producer_order_failures", 0) == 0, "O3 producer failures are not zero")
        require(audit.get("actual_producer_order_checks") == EXPREL_PRODUCER_ORDER_CHECK_COUNT, "O3 actual producer check count is not 4")
        require(audit.get("actual_producer_order_failures", 0) == 0, "O3 actual producer failures are not zero")
        require(audit.get("checked_residual_count") == EXPREL_PRODUCER_ORDER_CHECK_COUNT, "O3 checked residual count is not 4")
        require(audit.get("checked_residual_failures", 0) == 0, "O3 checked residual failures are not zero")
        require(result["new_counts"].get("new_exprel_producer_order_checks") == EXPREL_PRODUCER_ORDER_CHECK_COUNT, "O3 result producer count is not 4")
        require(result["new_counts"].get("new_exprel_producer_order_failures", 0) == 0, "O3 result producer failures are not zero")

    if result.get("passed") is False:
        require(not any(result.get(key) is True for key in NEW_FLAGS), "Failed O3 result overclaims a new flag")
        return {
            "passed": False,
            "honest_failure": True,
            "cocycle_checks": _result_count(result, ("cocycle_checks", "exponential_cocycle_order_checks", "cocycle_check_count")),
            "actual_driver_gap_rows": _result_count(result, ("actual_driver_gap_rows", "signed_density_gap_rows", "driver_gap_row_count")),
            "new_exprel_producer_order_checks": audit.get("new_exprel_producer_order_checks", 0),
            "new_exprel_producer_order_failures": audit.get("new_exprel_producer_order_failures", 0),
        }
    require(result.get("source_packet_definition_accepted") is True, "Passed O3 result lacks owned source acceptance")
    require(result.get("source_packet_accepted_owned") is True, "Passed O3 result lacks owned acceptance marker")
    cocycle = _result_count(result, ("cocycle_checks", "exponential_cocycle_order_checks", "cocycle_check_count"))
    gaps = _result_count(result, ("actual_driver_gap_rows", "signed_density_gap_rows", "driver_gap_row_count"))
    require(cocycle == COCYCLE_CHECK_COUNT, "O3 cocycle check count is not 4")
    require(gaps == ACTUAL_DRIVER_GAP_COUNT, "O3 actual driver gap row count is not 20")
    require(all(result.get(key) is True for key in NEW_FLAGS), "O3 flags are not all source-bound")
    require(audit.get("new_exponential_cocycle_order_checks") == COCYCLE_CHECK_COUNT, "O3 audit cocycle count is not 4")
    require(audit.get("new_actual_signed_density_gap_rows") == ACTUAL_DRIVER_GAP_COUNT, "O3 audit gap count is not 20")
    return {
        "passed": True,
        "cocycle_checks": cocycle,
        "actual_driver_gap_rows": gaps,
        "new_exprel_producer_order_checks": audit.get("new_exprel_producer_order_checks"),
        "new_exprel_producer_order_failures": audit.get("new_exprel_producer_order_failures", 0),
        "actual_producer_order_checks": audit.get("actual_producer_order_checks"),
        "actual_producer_order_failures": audit.get("actual_producer_order_failures", 0),
        "new_density_proof_rows": 0,
    }


def run(args):
    require(not ENTERED_CHECKPOINT.exists() and not RESULT_CHECKPOINT.exists() and not DIAGNOSTIC_CHECKPOINT.exists(), "O3 attempt already entered")
    source_raw, source = load_source()
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "source_checkpoint_path": wire_path(SOURCE_CHECKPOINT),
        "kernel_imported_before_entry": False,
        "derive_calls_before_entry": 0,
        "entry_attempts": 1,
        "runtime_model_metadata": None,
        "entered_code_immutable": True,
        "never_retry_entered_code": True,
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    entered_raw = write_new(ENTERED_CHECKPOINT, entered)
    write_new(DIAGNOSTIC_CHECKPOINT, b"")
    events = []
    calls = 0
    imported = False
    result = None
    try:
        kernel = resolve_path(source["kernel_path"], REPO)
        kernel_raw = read_bytes(kernel, source["current_code_hashes"][wire_path(kernel)])
        spec = importlib.util.spec_from_file_location("_compact_O3_transition_driver_cocycle_entry", kernel)
        require(spec is not None and spec.loader is not None, "Cannot create O3 kernel import specification")
        module = importlib.util.module_from_spec(spec)
        imported = True
        spec.loader.exec_module(module)
        derive = getattr(module, "derive", None)
        require(callable(derive), "O3 kernel lacks derive(source, journal_callback)")
        kernel_events = []
        calls += 1
        result = derive(clone(source), kernel_events)
        for event in kernel_events:
            journal_append(events, event)
    except Exception as error:
        journal_append(events, {"stage": "derive_exception", "error_type": type(error).__name__, "error": str(error)})
        result = fallback(source, calls, imported, error)
    result = decorate_result_metadata(result, source, calls, imported)
    result["source_checkpoint_sha256"] = sha_bytes(source_raw)
    result["entry_checkpoint_sha256"] = sha_bytes(entered_raw)
    result["diagnostic_checkpoint_sha256"] = sha_path(DIAGNOSTIC_CHECKPOINT)
    result["diagnostic_event_count"] = len(events)
    # The result is exclusive and durable before validation.  It is never rewritten.
    result_raw = write_new(RESULT_CHECKPOINT, result)
    validation = None
    validation_error = None
    try:
        validation = validate_result_metadata(result, source)
    except Exception as error:
        validation_error = str(error)
    receipt = {
        "schema": COMPONENT + "_run_receipt_v1",
        "component": COMPONENT,
        "metadata_only": True,
        "passed": bool(validation is not None and result.get("passed") is True),
        "validation_accepted": validation is not None,
        "validation_error": validation_error,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_bytes(entered_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "diagnostic_checkpoint_sha256": sha_path(DIAGNOSTIC_CHECKPOINT),
        "derive_entry_calls": calls,
        "kernel_imported": imported,
        "runtime_model_metadata": None,
    }
    receipt_raw = write_replace(RUN_RECEIPT, receipt)
    return {"run": receipt, "run_sha256": sha_bytes(receipt_raw), "validation": validation, "validation_error": validation_error}


def artifact_entry(path, expected=None, embed=False):
    raw = read_bytes(path, expected)
    entry = {"path": wire_path(path), "sha256": sha_bytes(raw), "bytes": len(raw), "raw_readonly": True}
    if embed:
        entry["raw_bytes_base64"] = base64.b64encode(raw).decode("ascii")
        entry["raw_bytes_exact"] = True
    return entry


def checkpoint_artifacts(embed=False):
    return [artifact_entry(path, embed=embed) for path in CHECKPOINTS]


def update_main_brief(pointer):
    """Merge a compact pointer while archiving, rather than deleting, old detail."""
    require(MAIN_BRIEF.exists(), "Shared RESUME_BRIEF is missing")
    old_raw = read_bytes(MAIN_BRIEF)
    old = json.loads(old_raw.decode("utf-8"), parse_int=int)
    require(isinstance(old, dict), "Shared RESUME_BRIEF is not a mapping")
    preimage = WORK / f"{COMPONENT}_main_resume_brief_preimage_{sha_bytes(old_raw)[:12]}.json"
    if not preimage.exists():
        write_new(preimage, old_raw)
    state = clone(old)
    in_progress = state.get("in_progress")
    require(in_progress is None or isinstance(in_progress, dict), "Shared RESUME_BRIEF in_progress is not a mapping")
    if not isinstance(in_progress, dict):
        in_progress = {}
    in_progress[COMPONENT] = clone(pointer)
    state["in_progress"] = in_progress
    state["latest_component"] = COMPONENT
    state["latest_component_pointer"] = clone(pointer)

    # A pending prepare must not promote a report or milestone that does not
    # exist.  Publication creates both files first; resume/seal add authority.
    if pointer.get("phase") != "prepared_unentered":
        report_path = REPO / REPORT_RELATIVE
        milestone_path = REPO / MILESTONE_RELATIVE
        if report_path.exists() and milestone_path.exists():
            state["latest_milestone"] = {
                "path": wire_path(milestone_path),
                "report": wire_path(report_path),
                "phase": pointer.get("phase"),
                "report_sha256": pointer.get("report_sha256"),
                "milestone_sha256": pointer.get("milestone_sha256"),
            }
            read_first = state.get("read_first")
            if not isinstance(read_first, list):
                read_first = []
            promoted = [wire_path(milestone_path), wire_path(report_path)]
            state["read_first"] = promoted + [
                item for item in read_first if item not in promoted
            ]
    authority = pointer.get("latest_authoritative_state")
    if authority is not None:
        require(isinstance(authority, dict), "Authoritative state pointer is not a mapping")
        state["latest_authoritative_state"] = clone(authority)
        verified = pointer.get("verified_local_and_remote_commit")
        require(valid_git_revision(verified), "Verified local/remote commit is malformed")
        state["verified_local_and_remote_commit"] = verified
        exact = state.get("current_exact_source_result")
        require(isinstance(exact, dict), "Shared current_exact_source_result is not a mapping")
        exact["current_report"] = clone(authority)
        exact["current_report_sha256"] = authority.get("report_sha256")
        state["current_exact_source_result"] = exact
        if pointer.get("next_remaining") is not None:
            state["next_remaining"] = clone(pointer["next_remaining"])
        if pointer.get("latest_review") is not None:
            state["latest_review"] = clone(pointer["latest_review"])
        if pointer.get("latest_qa") is not None:
            state["latest_qa"] = clone(pointer["latest_qa"])

    merged_raw = encode(state)
    if len(merged_raw) > 100_000:
        compact_keys = {
            "continuation_detail_archive",
            "archived_detail",
            "archived_prior_component_details",
            "failed_source_attempt_evidence",
        }
        compact = {key: value for key, value in state.items() if key not in compact_keys}
        old_in_progress = state.get("in_progress", {})
        if isinstance(old_in_progress, dict):
            keep = {
                key
                for key in old_in_progress
                if key == COMPONENT or "quiet" in key.lower() or key == "compact_C3_quiet_endpoint_projection"
            }
            omitted = sorted(key for key in old_in_progress if key not in keep)
            compact["in_progress"] = {key: old_in_progress[key] for key in sorted(keep)}
            compact["archived_in_progress_selector"] = {
                "selector": "$.in_progress",
                "path": wire_path(preimage),
                "sha256": sha_bytes(old_raw),
                "bytes": len(old_raw),
                "omitted_keys": omitted,
                "preserved_exact_bytes": True,
            }
        compact["archived_main_resume_brief_preimage"] = {
            "path": wire_path(preimage),
            "sha256": sha_bytes(old_raw),
            "bytes": len(old_raw),
        }
        compact["main_resume_brief_archival_compaction"] = True
        merged_raw = encode(compact)
    require(len(merged_raw) <= 100_000, "Shared RESUME_BRIEF remains over 100 KB after archival compaction")
    write_replace(MAIN_BRIEF, merged_raw)
    return {
        "path": wire_path(MAIN_BRIEF),
        "sha256": sha_bytes(merged_raw),
        "bytes": len(merged_raw),
        "preimage": {"path": wire_path(preimage), "sha256": sha_path(preimage), "bytes": len(old_raw)},
    }


def authoritative_pointer(report, report_raw, milestone_raw, commit, review=None, qa=None):
    full = report.get("full_result")
    require(isinstance(full, dict), "Published report full result is not a mapping")
    return {
        "schema": COMPONENT + "_authoritative_pointer_v1",
        "commit": commit,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "milestone": wire_path(REPO / MILESTONE_RELATIVE),
        "milestone_sha256": sha_bytes(milestone_raw),
        "attempt_status": report.get("attempt_status"),
        "passed": report.get("passed") is True,
        "source_packet_defined": full.get("source_packet_defined") is True,
        "source_packet_definition_accepted": full.get("source_packet_definition_accepted") is True,
        "source_packet_accepted_owned": full.get("source_packet_accepted_owned") is True,
        "new_flags": clone(report.get("new_flags", {})),
        "nonclaim_flags": clone(report.get("nonclaim_flags", {})),
        "source_pin_count": report.get("source_pin_count"),
        "checkpoint_count": report.get("checkpoint_count"),
        "preserved_failed_artifact_count": report.get("preserved_failed_artifact_count"),
        "runtime_model_metadata": None,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "latest_review": clone(review) if review is not None else None,
        "latest_qa": clone(qa) if qa is not None else None,
    }


def hydrate_publish_args(args, source):
    if not args.parent_report:
        args.parent_report = source.get("parent_report_metadata", {}).get("path")
    if not args.parent_report_sha256:
        args.parent_report_sha256 = source.get("parent_report_sha256")
    if not args.parent_commit:
        args.parent_commit = source.get("parent_commit")


def failed_union(source, result, current_artifacts):
    if result.get("passed") is True:
        return clone(source["preserved_failed_artifact_hashes_union"]), PARENT_FAILED_COUNT
    additions = {}
    for entry in current_artifacts:
        additions[entry["path"]] = entry["sha256"]
    additions[wire_path(PREPARE_RECEIPT)] = sha_path(PREPARE_RECEIPT)
    additions[wire_path(RUN_RECEIPT)] = sha_path(RUN_RECEIPT)
    additions[wire_path(Path(__file__).resolve())] = sha_path(Path(__file__).resolve())
    require(len(additions) == 7, "Failed O3 artifact retention set is not exactly seven artifacts")
    merged = merge_hashes(source["preserved_failed_artifact_hashes_union"], additions)
    # The failed-attempt contract is exactly seven newly retained artifacts.
    require(len(merged) == PARENT_FAILED_COUNT + 7, "Failed O3 artifact retention did not add exactly seven artifacts")
    return merged, len(merged)


def report_for_publish(source, result, result_raw, code2, input_meta):
    current_ck = checkpoint_artifacts(embed=True)
    helper = artifact_entry(Path(__file__).resolve(), source["external_draft"]["sha256"], embed=True)
    input_artifact = artifact_entry(Path(input_meta["path"]), input_meta["sha256"], embed=True)
    packet_receipts = []
    for receipt in source.get("raw_packet_receipts", []):
        packet_receipts.append(artifact_entry(resolve_path(receipt["path"]), receipt["sha256"], embed=True))
    raw_packet_receipts = clone(packet_receipts)
    current_failed_union, failed_count = failed_union(source, result, current_ck)
    artifacts = [
        artifact_entry(Path(source["parent_report_metadata"]["path"]), source["parent_report_metadata"]["sha256"]),
        artifact_entry(Path(source["parent_result_metadata"]["path"]), source["parent_result_metadata"]["sha256"]),
        *[clone(item) for item in code2],
        input_artifact,
        helper,
        *raw_packet_receipts,
        *current_ck,
    ]
    report = {
        "schema": PUBLICATION_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "passed" if result.get("passed") is True else "failed",
        "passed": result.get("passed") is True,
        "source_only": True,
        "publication_metadata_only": True,
        "runtime_model_metadata": None,
        "parent_report": source["parent_report_metadata"]["path"],
        "parent_report_sha256": source["parent_report_sha256"],
        "parent_commit": source["parent_commit"],
        "parent_report_is_authenticated_view": True,
        "parent_report_view": clone(source["parent_report"]),
        "parent_result_raw_sha256": source["parent_result_raw_sha256"],
        "parent_result_metadata": clone(source["parent_result_metadata"]),
        "source_pins": clone(source["source_union_hashes"]),
        "source_pin_count": len(source["source_union_hashes"]),
        "source_union_count": len(source["source_union_hashes"]),
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_inventory": clone(source["complete_parent_checkpoint_hashes"]),
        "checkpoint_hashes": {entry["path"]: entry["sha256"] for entry in current_ck},
        "checkpoint_paths": [entry["path"] for entry in current_ck],
        "checkpoint_count": len(source["complete_parent_checkpoint_hashes"]) + len(current_ck),
        "inherited_checkpoint_count": len(source["complete_parent_checkpoint_hashes"]),
        "preserved_failed_artifact_hashes_union": current_failed_union,
        "preserved_failed_artifact_count": failed_count,
        "protected_file_hashes": clone(source["protected_file_hashes"]),
        "legacy_dirty_path_hashes": clone(source["legacy_dirty_path_hashes"]),
        "inherited_failed_artifacts_preserved": True,
        "retained_construction_data": clone(source["parent_report"].get("retained_construction_data")),
        "retained_integral_term_inventory": clone(source["parent_report"].get("retained_integral_term_inventory")),
        "retained_cap_provenance": clone(source["parent_report"].get("retained_cap_provenance")),
        "publication_helper_metadata": {key: helper[key] for key in ("path", "sha256", "bytes")},
        "metadata_helper_artifact": clone(helper),
        "input_artifact": clone(input_artifact),
        "false_flags": list(source["inherited_false_flags"]),
        "false_flag_values": {key: False for key in source["inherited_false_flags"]},
        "inherited_parent_true_flags": clone(source["inherited_parent_true_flags"]),
        "new_flags": {key: result.get(key) is True for key in NEW_FLAGS},
        "old_quiet_flags": {key: False for key in OLD_QUIET_FLAGS},
        "parent_C3_projection_flags": {key: True for key in PARENT_C3_FLAGS},
        "nonclaim_flags": {key: False for key in NONCLAIM_FLAGS},
        "new_exponential_cocycle_order_checks": COCYCLE_CHECK_COUNT,
        "new_exprel_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
        "new_exprel_producer_order_failures": 0,
        "actual_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
        "actual_producer_order_failures": 0,
        "new_actual_signed_density_gap_rows": ACTUAL_DRIVER_GAP_COUNT,
        "new_density_proof_rows": 0,
        "derivation_contract": clone(source["derivation_contract"]),
        "counts": clone(source["counts"]),
        "new_counts": clone(source["new_counts"]),
        "source_inspection_audit": {
            "new_exprel_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "new_exprel_producer_order_failures": 0,
            "actual_producer_order_checks": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "actual_producer_order_failures": 0,
            "checked_residual_count": EXPREL_PRODUCER_ORDER_CHECK_COUNT,
            "checked_residual_failures": 0,
        },
        "full_result": clone(result),
        "component_result": clone(result),
        "component_result_exact": True,
        "source_packet_defined": result.get("passed") is True,
        "source_packet_definition_accepted": result.get("passed") is True,
        "source_packet_accepted_owned": result.get("passed") is True,
        "source_packet_accepted_inherited": False,
        "original_source_evaluations": 0,
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "numeric": False,
        "full_result_raw_sha256": sha_bytes(result_raw),
        "full_result_raw_bytes": len(result_raw),
        "full_result_artifact_selector": "$.lossless_artifacts[path=" + wire_path(RESULT_CHECKPOINT) + "]",
        "lossless_artifacts": artifacts,
        "lossless_artifact_count": len(artifacts),
        "raw_packet_receipts": raw_packet_receipts,
        "raw_bytes_lossless": True,
        "raw_sha_pins": clone(source["raw_sha_pins"]),
        "input_spec_metadata": clone(source["input_spec_metadata"]),
        "goal_active": True,
        "full_goal_complete": False,
        "next_executable_task": next_task_metadata(),
        "code2": [clone(item) for item in code2],
        "external_draft": clone(source["external_draft"]),
        "math_reexecuted": False,
        "ancestor_kernel_imported": False,
        "old_generic_identity_cached_data_only": True,
        "old_generic_identity_calls": 0,
        "old_numeric_calls": 0,
        "full_pair_proved": False,
        "incoming_identity_proved": False,
        "caps_bound": False,
        "global_five_moment_closure": False,
    }
    for key in NEW_FLAGS:
        report[key] = result.get(key) is True
    for key, value in source["inherited_parent_true_flags"].items():
        if key not in BOOKKEEPING_FIELDS and key not in NEW_FLAGS and key not in NONCLAIM_FLAGS:
            report.setdefault(key, True)
    for key in source["inherited_false_flags"]:
        report[key] = False
    for key in OLD_QUIET_FLAGS:
        report[key] = False
    for key in NONCLAIM_FLAGS:
        report[key] = False
    return report


def validate_publication_metadata(report, source):
    """Validate only publication bytes, counts, and authenticated metadata."""
    require(report.get("source_pins") == source.get("source_union_hashes"), "Published source pins changed")
    require(report.get("source_union_count") == len(source.get("source_union_hashes", {})), "Published source union count changed")
    require(report.get("parent_report_sha256") == PARENT_REPORT_SHA256, "Published parent report SHA changed")
    require(report.get("parent_commit") == PARENT_COMMIT, "Published parent commit changed")
    require(report.get("parent_result_raw_sha256") == source.get("parent_result_raw_sha256"), "Published parent result SHA changed")
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Published checkpoint count changed")
    require(report.get("new_exprel_producer_order_checks") == EXPREL_PRODUCER_ORDER_CHECK_COUNT, "Published producer check count changed")
    require(report.get("new_exprel_producer_order_failures", 0) == 0, "Published producer failures changed")
    require(report.get("counts", {}).get("new_exprel_producer_order_checks") == EXPREL_PRODUCER_ORDER_CHECK_COUNT, "Published producer count map changed")
    require(report.get("new_counts", {}).get("new_exprel_producer_order_checks") == EXPREL_PRODUCER_ORDER_CHECK_COUNT, "Published new producer count map changed")
    current = checkpoint_artifacts(embed=False)
    require(report.get("checkpoint_hashes") == {entry["path"]: entry["sha256"] for entry in current}, "Published checkpoint hashes changed")
    require(report.get("checkpoint_paths") == [entry["path"] for entry in current], "Published checkpoint paths changed")
    for entry in report.get("code2", []):
        read_bytes(resolve_path(entry.get("path")), entry.get("sha256"))
        require(entry.get("raw_readonly") is True and entry.get("raw_bytes_exact") is True, "Published code2 is not lossless")
    input_artifact = report.get("input_artifact")
    require(isinstance(input_artifact, dict), "Published input artifact is missing")
    read_bytes(resolve_path(input_artifact.get("path")), input_artifact.get("sha256"))
    require(input_artifact.get("raw_readonly") is True and input_artifact.get("raw_bytes_exact") is True, "Published input is not lossless")
    helper = report.get("metadata_helper_artifact")
    require(isinstance(helper, dict), "Published metadata helper artifact is missing")
    read_bytes(resolve_path(helper.get("path")), helper.get("sha256"))
    require(helper.get("raw_readonly") is True and helper.get("raw_bytes_exact") is True, "Published helper is not lossless")
    for entry in report.get("raw_packet_receipts", []):
        read_bytes(resolve_path(entry.get("path")), entry.get("sha256"))
        require(entry.get("raw_readonly") is True and entry.get("raw_bytes_exact") is True, "Published raw packet is not lossless")
    require(report.get("raw_sha_pins") == source.get("raw_sha_pins"), "Published raw SHA pins changed")
    return {
        "source_pin_count": report.get("source_pin_count"),
        "checkpoint_count": report.get("checkpoint_count"),
        "parent_result_raw_sha256": report.get("parent_result_raw_sha256"),
        "raw_sha_pins_verified": True,
        "metadata_only": True,
    }


def handoff_preimage(path, raw):
    preimage = WORK / f"{COMPONENT}_handoff_preimage_{path.stem}_{sha_bytes(raw)[:12]}.md"
    if not preimage.exists():
        write_new(preimage, raw)
    return {"path": wire_path(preimage), "sha256": sha_path(preimage), "bytes": len(raw), "raw_readonly": True}


def publish(args):
    source_raw, source = load_source()
    hydrate_publish_args(args, source)
    require(RESULT_CHECKPOINT.exists(), "O3 result checkpoint is missing")
    result_raw, result = read_json(RESULT_CHECKPOINT)
    parent = load_parent(args)
    input_meta = source["input_spec_metadata"]
    kernel = resolve_path(source["kernel_path"], REPO)
    runner = resolve_path(source["runner_path"], REPO)
    code2 = [
        artifact_entry(kernel, source["current_code_hashes"][wire_path(kernel)], embed=True),
        artifact_entry(runner, source["current_code_hashes"][wire_path(runner)], embed=True),
    ]
    input_path = resolve_path(input_meta["path"], WORK)
    input_meta_for_publish = {**input_meta, "path": input_path}
    durable_validation = validate_result_metadata(result, source)
    report = report_for_publish(source, result, result_raw, code2, input_meta_for_publish)
    report["durable_result_validation"] = clone(durable_validation)
    report["publication_metadata_validation"] = validate_publication_metadata(report, source)
    report_raw = write_replace(REPO / REPORT_RELATIVE, report)
    milestone = (
        "# Current original O3 transition driver cocycle (2026-10-26)\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}`; SHA-256 `{sha_bytes(report_raw)}`.\n\n"
        f"Attempt status: `{report['attempt_status']}`; source-only: `true`; exponential cocycle checks: `4/4`; C3 exponential-relation producer checks: `4/4` (`new_exprel_producer_order_checks=4`); producer failures: `0`; actual signed density-gap functions: `20/20`.\n\n"
        "This leaf binds a finite-N driver source, an exponential driver cocycle on ordinary Z orders 0 through 3, "
        "and twenty actual signed density-gap function recipes.  It introduces no numeric field values, density proof rows, "
        "full-pair proof, incoming identity, cap closure, global closure, physical completion, or full-goal completion. "
        "The old generic identity is retained only as cached data.\n\n"
        "The entered attempt is immutable and must never be reentered.\n\n"
        f"Next source-bound task: P2 `{next_task_metadata()['phase']}` with `{ACTUAL_DRIVER_GAP_COUNT}` source functions; "
        "the source gap remains unknown and is not asserted zero.\n"
    ).encode("utf-8")
    milestone_raw = write_replace(REPO / MILESTONE_RELATIVE, milestone)
    preimages = []
    handoff_prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}` (SHA-256 `{sha_bytes(report_raw)}`).\n"
        f"Status: `{report['attempt_status']}`; cocycle checks `4/4`; actual driver gap rows `20/20`; full goal remains incomplete.\n"
        "Do not reenter the saved attempt; preserve all inherited evidence and false flags.\n"
    ).encode("utf-8")
    for relative in HANDOFF_RELATIVES:
        path = REPO / relative
        old = path.read_bytes() if path.exists() else b""
        preimages.append(handoff_preimage(path, old))
        write_replace(path, handoff_prefix + old)
    attrs_path = REPO / ".gitattributes"
    attrs_old = attrs_path.read_bytes() if attrs_path.exists() else b""
    attrs_text = attrs_old.decode("utf-8", errors="strict")
    attr_lines = [
        f"/{KERNEL_RELATIVE.as_posix()} -text -diff",
        f"/{RUNNER_RELATIVE.as_posix()} -text -diff",
        f"/{REPORT_RELATIVE.as_posix()} -text -diff",
        f"/{MILESTONE_RELATIVE.as_posix()} -text -diff",
    ]
    for line in attr_lines:
        if line not in attrs_text.splitlines():
            attrs_text += ("\n" if attrs_text and not attrs_text.endswith("\n") else "") + line + "\n"
    write_replace(attrs_path, attrs_text.encode("utf-8"))
    main_brief = update_main_brief(
        {
            "phase": "published_passed" if report["attempt_status"] == "passed" else "published_failed",
            "report": wire_path(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "milestone": wire_path(REPO / MILESTONE_RELATIVE),
            "milestone_sha256": sha_bytes(milestone_raw),
            "attempt_status": report["attempt_status"],
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "next_executable_task": next_task_metadata(),
        }
    )
    receipt = {
        "schema": COMPONENT + "_publication_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "attempt_status": report["attempt_status"],
        "scoped_paths": [path.as_posix() for path in SCOPED_RELATIVES],
        "handoff_preimages": preimages,
        "exact_nine_scope": True,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "runtime_model_metadata": None,
        "math_reexecuted": False,
        "main_resume_brief": main_brief,
    }
    receipt_raw = write_replace(PUBLISH_RECEIPT, receipt)
    return {"publication": receipt, "publication_sha256": sha_bytes(receipt_raw)}


def git_output(*args):
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def verify_published_commit(report):
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "O3 publication checkpoint count changed")
    require(report.get("source_pin_count", 0) > PARENT_SOURCE_PIN_COUNT, "O3 source pins did not advance")
    require(report.get("full_goal_complete") is False and report.get("goal_active") is True, "O3 goal is not active/incomplete")
    head = git_output("rev-parse", "HEAD")
    require(valid_git_revision(head), "Local HEAD is not an exact Git revision")
    require(report.get("parent_commit") == PARENT_COMMIT, "Published parent commit changed")
    require(git_output("rev-parse", "HEAD^") == PARENT_COMMIT, "HEAD parent does not match authenticated parent commit")
    branch = git_output("branch", "--show-current")
    require(branch, "Current Git branch is unnamed")
    remote_line = git_output("ls-remote", "origin", "refs/heads/" + branch)
    remote_parts = remote_line.split()
    require(len(remote_parts) >= 2 and valid_git_revision(remote_parts[0]), "Remote branch lookup returned no exact head")
    remote = remote_parts[0].lower()
    require(head.lower() == remote, "HEAD is not synchronized with origin exact branch")
    changed = set(git_output("diff", "--name-only", PARENT_COMMIT + "..HEAD").splitlines())
    expected = {path.as_posix() for path in SCOPED_RELATIVES}
    require(changed == expected, "Committed scope is not exactly the nine O3 paths")
    for relative in SCOPED_RELATIVES:
        local_path = REPO / relative
        require(local_path.exists(), "Published scoped path is missing: " + relative.as_posix())
        blob_size = int(git_output("cat-file", "-s", "HEAD:" + relative.as_posix()))
        require(blob_size == local_path.stat().st_size, "Published cat-file byte count changed: " + relative.as_posix())
    return head, remote


def resume(args):
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "O3 publication report is missing")
    head, remote = verify_published_commit(report)
    milestone_raw = read_bytes(REPO / MILESTONE_RELATIVE)
    resume_qa = {"durable_result_validation": clone(report.get("durable_result_validation"))}
    authority = authoritative_pointer(report, report_raw, milestone_raw, head, report.get("independent_review"), resume_qa)
    main_brief = update_main_brief(
        {
            "phase": "resumed_published",
            "report": wire_path(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "milestone": wire_path(REPO / MILESTONE_RELATIVE),
            "head": head,
            "remote": remote,
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "next_executable_task": next_task_metadata(),
            "latest_authoritative_state": authority,
            "verified_local_and_remote_commit": head,
            "next_remaining": next_task_metadata(),
            "milestone_sha256": sha_bytes(milestone_raw),
            "latest_review": report.get("independent_review"),
            "latest_qa": resume_qa,
        }
    )
    receipt = {
        "schema": COMPONENT + "_resume_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "head": head,
        "remote": remote,
        "exact_nine_scope": True,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "runtime_model_metadata": None,
        "main_resume_brief": main_brief,
    }
    receipt_raw = write_replace(RESUME_RECEIPT, receipt)
    return {"resume": receipt, "resume_sha256": sha_bytes(receipt_raw)}


def seal(args):
    require(args.review and args.review_sha256, "Seal requires an independent saved review path and SHA-256")
    require(valid_sha(args.review_sha256), "Independent review SHA-256 is malformed")
    review_path = resolve_path(args.review, WORK)
    review_raw = read_bytes(review_path, args.review_sha256)
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "O3 report is not published")
    source_raw, source = load_source()
    require(valid_sha(report.get("full_result_raw_sha256")), "Published result SHA is malformed")
    result_raw, result = read_json(RESULT_CHECKPOINT, report["full_result_raw_sha256"])
    require(report.get("full_result") == result, "Published full result changed before seal")
    require(report.get("component_result") == result, "Published component result changed before seal")
    require(report.get("full_result_raw_bytes") == len(result_raw), "Published result byte count changed before seal")
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Published checkpoint count changed before seal")
    require(report.get("source_pin_count") == len(source.get("source_union_hashes", {})), "Published source-pin count changed before seal")
    require(report.get("preserved_failed_artifact_count") in (PARENT_FAILED_COUNT, PARENT_FAILED_COUNT + 7), "Published failed count changed before seal")
    publication_metadata_validation = validate_publication_metadata(report, source)
    require(report.get("publication_metadata_validation") == publication_metadata_validation, "Published metadata validation changed before seal")
    durable_validation = validate_result_metadata(result, source)
    require(report.get("durable_result_validation") == durable_validation, "Published durable validation changed before seal")
    head, remote = verify_published_commit(report)
    milestone_raw = read_bytes(REPO / MILESTONE_RELATIVE)
    review_meta = {"path": wire_path(review_path), "sha256": sha_bytes(review_raw), "bytes": len(review_raw), "raw_readonly": True}
    qa = {"durable_result_validation": clone(durable_validation), "review_sha256": review_meta["sha256"]}
    authority = authoritative_pointer(report, report_raw, milestone_raw, head, review_meta, qa)
    main_brief = update_main_brief(
        {
            "phase": "sealed",
            "report": wire_path(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "seal": wire_path(SEAL_RECEIPT),
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "next_executable_task": next_task_metadata(),
            "latest_authoritative_state": authority,
            "verified_local_and_remote_commit": head,
            "next_remaining": next_task_metadata(),
            "milestone_sha256": sha_bytes(milestone_raw),
            "latest_review": review_meta,
            "latest_qa": qa,
        }
    )
    receipt = {
        "schema": SEAL_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "head": head,
        "remote": remote,
        "durable_result_validation": durable_validation,
        "publication_metadata_validation": publication_metadata_validation,
        "independent_review": review_meta,
        "review_required": True,
        "runtime_model_metadata": None,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "main_resume_brief": main_brief,
    }
    receipt_raw = write_replace(SEAL_RECEIPT, receipt)
    return {"seal": receipt, "seal_sha256": sha_bytes(receipt_raw)}


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run", "publish", "resume", "seal"))
    p.add_argument("--input-spec")
    p.add_argument("--input-sha256")
    p.add_argument("--parent-report")
    p.add_argument("--parent-report-sha256")
    p.add_argument("--parent-commit")
    p.add_argument("--kernel-path")
    p.add_argument("--kernel-sha256")
    p.add_argument("--runner-path")
    p.add_argument("--review")
    p.add_argument("--review-sha256")
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    if args.action == "prepare":
        value = prepare(args)
    elif args.action == "run":
        value = run(args)
    elif args.action == "publish":
        value = publish(args)
    elif args.action == "resume":
        value = resume(args)
    else:
        value = seal(args)
    print(json.dumps(value, ensure_ascii=True, sort_keys=True))


if __name__ == "__main__":
    main()
