"""Metadata-only manager for the compact O3 transition primitive gap leaf.

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

COMPONENT = "compact_O3_transition_primitive_gap"
INPUT_SCHEMA = "compact_O3_transition_primitive_gap_input_v1"
RESULT_SCHEMA = "original_compact_O3_transition_primitive_gap"
ENTERED_SCHEMA = "original_compact_O3_transition_primitive_gap_entered_v1"
DIAGNOSTIC_SCHEMA = "compact_O3_transition_primitive_gap_stage_diagnostic_v1"
PUBLICATION_SCHEMA = "original_compact_O3_transition_primitive_gap_publication_v1"
SEAL_SCHEMA = "compact_O3_transition_primitive_gap_final_seal_v1"

KERNEL_RELATIVE = Path(
    "src/openai_ns_reconstruction/original_compact_O3_transition_primitive_gap.py"
)
RUNNER_RELATIVE = Path("scripts/run_original_compact_O3_transition_primitive_gap.py")
REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_primitive_gap.json"
)
MILESTONE_RELATIVE = Path(
    "docs/CURRENT_ORIGINAL_COMPACT_O3_TRANSITION_PRIMITIVE_GAP_2026_10_27.md"
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

BRIEF = WORK / f"{COMPONENT}_brief_2026_10_27.json"
MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
PREPARE_RECEIPT = WORK / f"{COMPONENT}_prepare_receipt_2026_10_27.json"
RUN_RECEIPT = WORK / f"{COMPONENT}_run_receipt_2026_10_27.json"
PUBLISH_RECEIPT = WORK / f"{COMPONENT}_publication_receipt_2026_10_27.json"
RESUME_RECEIPT = WORK / f"{COMPONENT}_resume_receipt_2026_10_27.json"
SEAL_RECEIPT = WORK / f"{COMPONENT}_final_seal_receipt_2026_10_27.json"

PARENT_REPORT_DEFAULT = (
    REPO
    / "experiments/root_st073/"
    / "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_driver_cocycle.json"
)
PARENT_REPORT_SHA256 = (
    "16fab70f7599c01a70159a344aa2386055f41f6af26a3145a216e132ed6fa5d2"
)
PARENT_COMMIT = "098c46878efd98c84ec2ef438c9d7f1596f529d4"
PARENT_REPORT_PARENT_COMMIT = "b0f7e8e0d9a69b8c82345ff0a6a356c8b7a2cc96"

PARENT_SOURCE_PIN_COUNT = 2207
PARENT_CHECKPOINT_COUNT = 240
CURRENT_CHECKPOINT_COUNT = PARENT_CHECKPOINT_COUNT + 4
PARENT_FAILED_COUNT = 169
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26

ACTUAL_PRIMITIVE_GAP_COUNT = 12
ROOT_SHAPE_ROW_COUNT = 4
ROOT_SHAPE_CHECKS_PER_ROW = 4
ROOT_SHAPE_SUBRESIDUAL_COUNT = 16
ROOT_SHAPE_ORDER_FAILURE_COUNT = 0
SIGMA_IDENTITY_CHECK_COUNT = 2
PRIMITIVE_COCYCLE_ORDER_CHECK_COUNT = 4
PRIMITIVE_COCYCLE_SUBRESIDUAL_COUNT = 8
PRIMITIVE_COCYCLE_RESIDUALS_PER_ROW = 2
PRIMITIVE_GAP_FAMILY_COUNTS = {"E": 4, "A": 4, "B": 4}
ORDINARY_Z_ORDERS = (0, 1, 2, 3)

NEW_FLAGS = (
    # Provisional names; the root may finalize these names before prepare.
    "O3_transition_primitive_sources_bound",
    "O3_transition_root_shape_correlations_proved",
    "O3_transition_signed_primitive_gap_functions_defined",
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
    "O3_transition_E_zero_claimed",
    "O3_transition_A_zero_claimed",
    "O3_transition_B_zero_claimed",
    "O3_transition_P2_zero_claimed",
    "O3_transition_phase_equality_proved",
    "O3_transition_signed_gap_zero_claimed",
    "numeric_values_computed",
    "numeric_containment_claimed",
    "physical_layer_completed",
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
    # Lossless payloads remain preserved in source/report artifacts, but are
    # omitted from the compact authenticated parent view so the main brief
    # and metadata packet stay portable and below 100 KB.
    "code2",
    "input_artifact",
    "metadata_helper_artifact",
    "raw_packet_receipts",
    "raw_sha_pins",
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
    require(
        report.get("component") == "compact_O3_transition_driver_cocycle",
        "Parent report component is not the authenticated preceding O3 transition",
    )
    require(
        report.get("parent_commit") == PARENT_REPORT_PARENT_COMMIT,
        "Parent report internal parent commit changed",
    )
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
        "parent_report_parent_commit": report.get("parent_commit"),
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
    # The external manager never manufactures or edits the runner.  Prepare
    # must freeze the root-owned runner bytes exactly as supplied.
    require(path.exists(), "Primitive-gap runner is missing; root must supply it before prepare")
    raw = read_bytes(path)
    return path, raw


def next_task_metadata():
    return {
        "schema": "O3_transition_primitive_source_gap_next_v1",
        "phase": "O3_transition_primitive_source_gap",
        "priority": "P2",
        "primitive": True,
        "source_gap": True,
        "action": "Bind actual E/A/B primitive source gaps and root-shape correlations on the preceding O3 transition.",
        "task": "Bind actual E/A/B primitive source functions for the next O3 transition source-gap check.",
        "primitive_families": ["E", "A", "B"],
        "family_function_counts": clone(PRIMITIVE_GAP_FAMILY_COUNTS),
        "source_function_count": ACTUAL_PRIMITIVE_GAP_COUNT,
        "signed_primitive_gap_function_count": ACTUAL_PRIMITIVE_GAP_COUNT,
        "function_count": ACTUAL_PRIMITIVE_GAP_COUNT,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "root_shape_row_count": ROOT_SHAPE_ROW_COUNT,
        "root_shape_checks_per_row": ROOT_SHAPE_CHECKS_PER_ROW,
        "root_shape_subresidual_count": ROOT_SHAPE_SUBRESIDUAL_COUNT,
        "root_shape_identity_fields": ["a", "b", "t0", "Delta"],
        "sigma_identity_families": [
            "positive_interior_logistic_reflection",
            "source_endpoint_binding",
        ],
        "sigma_identity_checks": SIGMA_IDENTITY_CHECK_COUNT,
        "primitive_cocycle_order_checks": PRIMITIVE_COCYCLE_ORDER_CHECK_COUNT,
        "primitive_cocycle_subresiduals": PRIMITIVE_COCYCLE_SUBRESIDUAL_COUNT,
        "primitive_cocycle_residuals_per_row": PRIMITIVE_COCYCLE_RESIDUALS_PER_ROW,
        "root_shape_sigma_same_exact_definition": True,
        "root_shape_sigma_interval_evaluated_value_equality": False,
        "gap_status": "unknown",
        "gap_zero_or_nonzero_asserted": False,
        "E_zero_claimed": False,
        "A_zero_claimed": False,
        "B_zero_claimed": False,
        "E_gap_zero_claimed": False,
        "A_gap_zero_claimed": False,
        "B_gap_zero_claimed": False,
        "P2_zero_claimed": False,
        "phase_equality_claimed": False,
        "phase_equality_proved": False,
        "signed_gap_zero_claimed": False,
        "gap_zero_claimed": False,
        "primitive_gap_zero_claimed": False,
        "source_only": True,
        "do_not_promote_primitive_facts_without_new_source_binding": True,
        "do_not_reenter_current_attempt": True,
    }


def count_contract(payload):
    """Return the fixed canonical O3 leaf counts without inspecting source values."""
    defaults = {
        "primitive_gap_functions": ACTUAL_PRIMITIVE_GAP_COUNT,
        "E_gap_functions": PRIMITIVE_GAP_FAMILY_COUNTS["E"],
        "A_gap_functions": PRIMITIVE_GAP_FAMILY_COUNTS["A"],
        "B_gap_functions": PRIMITIVE_GAP_FAMILY_COUNTS["B"],
        "root_shape_rows": ROOT_SHAPE_ROW_COUNT,
        "root_shape_checks_per_row": ROOT_SHAPE_CHECKS_PER_ROW,
        "root_shape_subresiduals": ROOT_SHAPE_SUBRESIDUAL_COUNT,
        "sigma_identity_checks": SIGMA_IDENTITY_CHECK_COUNT,
        "primitive_cocycle_order_checks": PRIMITIVE_COCYCLE_ORDER_CHECK_COUNT,
        "primitive_cocycle_subresiduals": PRIMITIVE_COCYCLE_SUBRESIDUAL_COUNT,
        "primitive_cocycle_residuals_per_row": PRIMITIVE_COCYCLE_RESIDUALS_PER_ROW,
    }
    declared = payload.get("count_contract", {})
    require(isinstance(declared, dict), "Input count_contract is not a mapping")
    require(set(declared).issubset(defaults), "Input count_contract has unsupported fields")
    result = dict(defaults)
    for key, default in defaults.items():
        if key in declared:
            require(declared[key] == default, "Input count_contract is not the canonical O3 leaf contract for " + key)
    require(
        result["primitive_gap_functions"] == result["E_gap_functions"] + result["A_gap_functions"] + result["B_gap_functions"],
        "Input primitive gap count does not equal E/A/B family counts",
    )
    require(
        result["root_shape_subresiduals"] == result["root_shape_rows"] * result["root_shape_checks_per_row"],
        "Input root-shape subresidual count does not match rows times checks",
    )
    require(
        result["primitive_cocycle_subresiduals"] == result["primitive_cocycle_order_checks"] * result["primitive_cocycle_residuals_per_row"],
        "Input primitive cocycle residual count does not match rows times residuals",
    )
    return result


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
    contract = count_contract(payload)
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
        "parent_head": PARENT_COMMIT,
        "parent_report_parent_commit": parent["parent_report_parent_commit"],
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
        "count_contract": clone(contract),
        "count_contract_parameterized": False,
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
        "inherited_code2": clone(parent["report"].get("code2", [])),
        "preserved_error_cap_metadata": {
            key: clone(parent["report"].get(key))
            for key in (
                "errors",
                "error_terms",
                "error_caps",
                "caps",
                "cap_provenance",
                "retained_cap_provenance",
                "retained_construction_data",
                "retained_integral_term_inventory",
                "inherited_failed_artifacts_preserved",
            )
            if key in parent["report"]
        },
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
            "new_root_shape_order_checks": ROOT_SHAPE_ROW_COUNT,
            "new_root_shape_order_failures": ROOT_SHAPE_ORDER_FAILURE_COUNT,
            "new_actual_primitive_gap_rows": contract["primitive_gap_functions"],
            "actual_primitive_gap_functions": contract["primitive_gap_functions"],
            "actual_E_gap_functions": contract["E_gap_functions"],
            "actual_A_gap_functions": contract["A_gap_functions"],
            "actual_B_gap_functions": contract["B_gap_functions"],
            "root_shape_ordinary_Z_rows": contract["root_shape_rows"],
            "root_shape_checks_per_row": contract["root_shape_checks_per_row"],
            "root_shape_subresiduals": contract["root_shape_subresiduals"],
            "new_sigma_identity_checks": contract["sigma_identity_checks"],
            "new_primitive_cocycle_order_checks": contract["primitive_cocycle_order_checks"],
            "new_primitive_cocycle_subresiduals": contract["primitive_cocycle_subresiduals"],
            "primitive_cocycle_residuals_per_row": contract["primitive_cocycle_residuals_per_row"],
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        },
        "new_counts": {
            "new_root_shape_order_checks": contract["root_shape_rows"],
            "new_root_shape_order_failures": ROOT_SHAPE_ORDER_FAILURE_COUNT,
            "new_actual_primitive_gap_rows": contract["primitive_gap_functions"],
            "actual_primitive_gap_functions": contract["primitive_gap_functions"],
            "actual_E_gap_functions": contract["E_gap_functions"],
            "actual_A_gap_functions": contract["A_gap_functions"],
            "actual_B_gap_functions": contract["B_gap_functions"],
            "root_shape_ordinary_Z_rows": contract["root_shape_rows"],
            "root_shape_checks_per_row": contract["root_shape_checks_per_row"],
            "root_shape_subresiduals": contract["root_shape_subresiduals"],
            "new_sigma_identity_checks": contract["sigma_identity_checks"],
            "new_primitive_cocycle_order_checks": contract["primitive_cocycle_order_checks"],
            "new_primitive_cocycle_subresiduals": contract["primitive_cocycle_subresiduals"],
            "primitive_cocycle_residuals_per_row": contract["primitive_cocycle_residuals_per_row"],
            "new_density_proof_rows": 0,
        },
        "derivation_contract": {
            "primitive_gap": "signed_E_A_B_primitive_gap_functions",
            "root_shape_correlations": True,
            "cached_ancestor_sections": {
                "exponential_F_driver_cocycle": {"cached_only": True, "used_for_new_proof": False},
                "exprel_producer_order_checks": {"cached_only": True, "used_for_new_proof": False},
            },
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "new_root_shape_order_checks": contract["root_shape_rows"],
            "new_root_shape_order_failures": ROOT_SHAPE_ORDER_FAILURE_COUNT,
            "new_actual_primitive_gap_rows": contract["primitive_gap_functions"],
            "actual_primitive_gap_functions": contract["primitive_gap_functions"],
            "actual_E_gap_functions": contract["E_gap_functions"],
            "actual_A_gap_functions": contract["A_gap_functions"],
            "actual_B_gap_functions": contract["B_gap_functions"],
            "root_shape_ordinary_Z_rows": contract["root_shape_rows"],
            "root_shape_checks_per_row": contract["root_shape_checks_per_row"],
            "root_shape_subresiduals": contract["root_shape_subresiduals"],
            "new_sigma_identity_checks": contract["sigma_identity_checks"],
            "new_primitive_cocycle_order_checks": contract["primitive_cocycle_order_checks"],
            "new_primitive_cocycle_subresiduals": contract["primitive_cocycle_subresiduals"],
            "primitive_cocycle_residuals_per_row": contract["primitive_cocycle_residuals_per_row"],
            "function_not_samples": True,
            "new_density_proof_rows": 0,
            "numeric": False,
            "full_pair_proved": False,
            "incoming_identity_proved": False,
            "caps_bound": False,
            "global_five_moment_closure": False,
        },
        "source_inspection_contract": {
            "new_root_shape_order_checks": contract["root_shape_rows"],
            "new_root_shape_order_failures": ROOT_SHAPE_ORDER_FAILURE_COUNT,
            "new_actual_primitive_gap_rows": contract["primitive_gap_functions"],
            "actual_primitive_gap_functions": contract["primitive_gap_functions"],
            "actual_E_gap_functions": contract["E_gap_functions"],
            "actual_A_gap_functions": contract["A_gap_functions"],
            "actual_B_gap_functions": contract["B_gap_functions"],
            "root_shape_ordinary_Z_rows": contract["root_shape_rows"],
            "root_shape_checks_per_row": contract["root_shape_checks_per_row"],
            "root_shape_subresiduals": contract["root_shape_subresiduals"],
            "new_sigma_identity_checks": contract["sigma_identity_checks"],
            "new_primitive_cocycle_order_checks": contract["primitive_cocycle_order_checks"],
            "new_primitive_cocycle_subresiduals": contract["primitive_cocycle_subresiduals"],
            "primitive_cocycle_residuals_per_row": contract["primitive_cocycle_residuals_per_row"],
        },
        "sigma_definition_contract": {
            "root_shape_sigma_same_exact_definition": True,
            "root_shape_sigma_interval_evaluated_value_equality": False,
            "crossmultiply_identity_families": [
                "positive_interior_logistic_reflection",
                "source_endpoint_binding",
            ],
        },
        "old_generic_identity_cached_data_only": True,
        "old_generic_identity_calls": 0,
        "old_numeric_calls": 0,
        "old_math_calls": 0,
        "old_source_calls": 0,
        "old_provider_calls": 0,
        "old_checker_calls": 0,
        "old_normalizer_calls": 0,
        "old_source_value_reducer_calls": 0,
        "old_operator_calls": 0,
        "old_source_consumer_calls": 0,
        "source_value_reducer_calls": 0,
        "numeric_field_evaluations": 0,
        "numeric_calls": 0,
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
        "E_zero_claimed": False,
        "A_zero_claimed": False,
        "B_zero_claimed": False,
        "E_gap_zero_claimed": False,
        "A_gap_zero_claimed": False,
        "B_gap_zero_claimed": False,
        "P2_zero_claimed": False,
        "phase_equality_claimed": False,
        "phase_equality_proved": False,
        "numerical": False,
        "global_closure": False,
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
        "new_root_shape_order_checks": source["count_contract"]["root_shape_rows"],
        "new_root_shape_order_failures": ROOT_SHAPE_ORDER_FAILURE_COUNT,
        "new_actual_primitive_gap_rows": source["count_contract"]["primitive_gap_functions"],
        "actual_primitive_gap_functions": source["count_contract"]["primitive_gap_functions"],
        "root_shape_ordinary_Z_rows": source["count_contract"]["root_shape_rows"],
        "root_shape_subresiduals": source["count_contract"]["root_shape_subresiduals"],
        "new_sigma_identity_checks": source["count_contract"]["sigma_identity_checks"],
        "new_primitive_cocycle_order_checks": source["count_contract"]["primitive_cocycle_order_checks"],
        "new_primitive_cocycle_subresiduals": source["count_contract"]["primitive_cocycle_subresiduals"],
        "primitive_cocycle_residuals_per_row": source["count_contract"]["primitive_cocycle_residuals_per_row"],
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
    _prepare_raw, prepare = read_json(PREPARE_RECEIPT)
    require(prepare.get("schema") == COMPONENT + "_prepare_receipt_v1", "O3 prepare receipt schema changed")
    require(prepare.get("component") == COMPONENT, "O3 prepare receipt component changed")
    require(prepare.get("metadata_only") is True, "O3 prepare receipt scope changed")
    require(prepare.get("source_checkpoint") == wire_path(SOURCE_CHECKPOINT), "O3 prepare source path changed")
    require(prepare.get("source_checkpoint_sha256") == sha_bytes(source_raw), "O3 source SHA is not authenticated by the durable prepare receipt")
    require(source.get("schema") == INPUT_SCHEMA and source.get("component") == COMPONENT, "O3 source schema changed")
    require(source.get("metadata_only") is True and source.get("source_only") is True, "O3 source scope changed")
    require(source.get("runtime_model_metadata") is None, "O3 source runtime metadata is not null")
    require(source.get("parent_report_sha256") == PARENT_REPORT_SHA256, "Source parent report SHA changed")
    require(source.get("parent_commit") == PARENT_COMMIT, "Source parent commit changed")
    require(source.get("parent_head") == PARENT_COMMIT, "Source parent HEAD changed")
    require(source.get("parent_report_parent_commit") == PARENT_REPORT_PARENT_COMMIT, "Source parent report parent commit changed")
    require(source.get("parent_source_pin_count") == PARENT_SOURCE_PIN_COUNT, "Source parent pin count changed")
    require(source.get("complete_parent_checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Source checkpoint count changed")
    require(source.get("current_checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Source current checkpoint count changed")
    require(source.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT, "Source failed-artifact count changed")
    require(len(source.get("protected_file_hashes", {})) == PARENT_PROTECTED_COUNT, "Source protected count changed")
    require(len(source.get("legacy_dirty_path_hashes", {})) == PARENT_DIRTY_COUNT, "Source dirty count changed")
    require(source.get("input", {}).get("schema") == INPUT_SCHEMA, "Source input schema changed")
    require(source.get("count_contract") == count_contract(source.get("input", {})), "Source count contract changed")
    require(source.get("new_flags") == list(NEW_FLAGS), "Source provisional new flags changed before freeze")
    input_spec = source.get("input_spec_metadata", {})
    require(input_spec.get("sha256") == source.get("input_raw_sha256"), "Source input pin changed")
    require(prepare.get("input_raw_sha256") == input_spec.get("sha256"), "Prepare input pin does not match source")
    read_bytes(resolve_path(input_spec.get("path"), WORK), input_spec.get("sha256"))
    code_hashes = source.get("current_code_hashes", {})
    kernel_path = resolve_path(source.get("kernel_path"), REPO)
    runner_path = resolve_path(source.get("runner_path"), REPO)
    kernel_key = wire_path(kernel_path)
    runner_key = wire_path(runner_path)
    require(prepare.get("kernel_sha256") == code_hashes.get(kernel_key), "Prepare kernel pin does not match source")
    require(prepare.get("runner_sha256") == code_hashes.get(runner_key), "Prepare runner pin does not match source")
    read_bytes(kernel_path, prepare.get("kernel_sha256"))
    read_bytes(runner_path, prepare.get("runner_sha256"))
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
        "new_algebra_family_calls": 0,
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
            "new_algebra_family_calls": 0,
            "new_root_shape_order_checks": 0,
            "new_root_shape_order_failures": 0,
            "new_actual_primitive_gap_rows": 0,
            "actual_primitive_gap_functions": 0,
            "actual_E_gap_functions": 0,
            "actual_A_gap_functions": 0,
            "actual_B_gap_functions": 0,
            "root_shape_ordinary_Z_rows": 0,
            "root_shape_checks_per_row": 0,
            "root_shape_subresiduals": 0,
            "new_sigma_identity_checks": 0,
            "new_primitive_cocycle_order_checks": 0,
            "new_primitive_cocycle_subresiduals": 0,
            "new_root_shape_order_checks": 0,
            "new_actual_primitive_gap_rows": 0,
            "new_density_proof_rows": 0,
        },
        "numeric": False,
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "new_root_shape_order_checks": 0,
            "new_root_shape_order_failures": 0,
            "new_actual_primitive_gap_rows": 0,
            "actual_primitive_gap_functions": 0,
            "actual_E_gap_functions": 0,
            "actual_A_gap_functions": 0,
            "actual_B_gap_functions": 0,
            "root_shape_ordinary_Z_rows": 0,
            "root_shape_checks_per_row": 0,
            "root_shape_subresiduals": 0,
            "new_sigma_identity_checks": 0,
            "new_primitive_cocycle_order_checks": 0,
            "new_primitive_cocycle_subresiduals": 0,
            "new_algebra_family_calls": 0,
            "new_density_proof_rows": 0,
            "old_generic_identity_calls": 0,
            "old_provider_calls": 0,
            "old_checker_calls": 0,
            "old_normalizer_calls": 0,
            "old_source_value_reducer_calls": 0,
            "old_numeric_calls": 0,
            "numeric_calls": 0,
            "old_math_calls": 0,
            "old_source_calls": 0,
            "numeric_values_computed": 0,
        },
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "old_generic_identity_cached_data_only": True,
        "old_numeric_calls": 0,
        "old_math_calls": 0,
        "old_source_calls": 0,
        "old_provider_calls": 0,
        "old_checker_calls": 0,
        "old_normalizer_calls": 0,
        "old_source_value_reducer_calls": 0,
        "numeric_calls": 0,
        "full_pair_proved": False,
        "incoming_identity_proved": False,
        "caps_bound": False,
        "global_five_moment_closure": False,
        "E_zero_claimed": False,
        "A_zero_claimed": False,
        "B_zero_claimed": False,
        "P2_zero_claimed": False,
        "phase_equality_claimed": False,
        "numerical": False,
        "global_closure": False,
        "signed_gap_zero_claimed": False,
        "gap_zero_claimed": False,
        "primitive_gap_zero_claimed": False,
    }


def decorate_result_metadata(value, source, calls, imported):
    result = clone(value) if isinstance(value, dict) else fallback(source, calls, imported, "derive returned a non-mapping result")
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("passed", False)
    result["passed"] = result.get("passed") is True
    result["new_algebra_family_calls"] = 1 if result.get("passed") is True else 0
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
    supplied_flags = result.get("new_flags", {})
    if not isinstance(supplied_flags, dict):
        supplied_flags = {}
    for key in NEW_FLAGS:
        result[key] = result.get(key) is True or supplied_flags.get(key) is True
    result["new_flags"] = {key: result.get(key) is True for key in NEW_FLAGS}
    # Missing source-bound counts remain zero; the manager never fabricates
    # primitive or root-shape facts from the expected contract.
    result.setdefault("new_counts", {})
    if not isinstance(result.get("new_counts"), dict):
        result["new_counts"] = {}
    result["new_counts"]["new_algebra_family_calls"] = result["new_algebra_family_calls"]
    result["new_counts"].setdefault("new_root_shape_order_checks", 0)
    result["new_counts"].setdefault("new_root_shape_order_failures", 0)
    result["new_counts"].setdefault("new_actual_primitive_gap_rows", 0)
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
        "old_operator_calls",
        "old_source_consumer_calls",
        "source_value_reducer_calls",
        "numeric_field_evaluations",
        "old_numeric_calls",
        "numeric_calls",
        "old_math_calls",
        "old_source_calls",
    ):
        result.setdefault(key, 0)
    for key in NONCLAIM_FLAGS:
        result.setdefault(key, False)
    for key in ("E_zero_claimed", "A_zero_claimed", "B_zero_claimed", "E_gap_zero_claimed", "A_gap_zero_claimed", "B_gap_zero_claimed", "P2_zero_claimed", "phase_equality_claimed", "phase_equality_proved", "signed_gap_zero_claimed", "gap_zero_claimed", "primitive_gap_zero_claimed"):
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
    audit["new_algebra_family_calls"] = result["new_algebra_family_calls"]
    audit.setdefault("new_root_shape_order_checks", result.get("root_shape_ordinary_Z_rows", 0))
    audit.setdefault("new_root_shape_order_failures", 0)
    audit.setdefault("new_actual_primitive_gap_rows", result.get("actual_primitive_gap_functions", 0))
    audit.setdefault("new_root_shape_order_checks", result.get("root_shape_ordinary_Z_rows", 0))
    audit.setdefault("actual_primitive_gap_functions", result.get("actual_primitive_gap_functions", 0))
    audit.setdefault("actual_E_gap_functions", 0)
    audit.setdefault("actual_A_gap_functions", 0)
    audit.setdefault("actual_B_gap_functions", 0)
    audit.setdefault("root_shape_ordinary_Z_rows", result.get("root_shape_ordinary_Z_rows", 0))
    audit.setdefault("root_shape_checks_per_row", 0)
    audit.setdefault("root_shape_subresiduals", result.get("root_shape_subresiduals", 0))
    audit.setdefault("new_sigma_identity_checks", result.get("new_sigma_identity_checks", 0))
    audit.setdefault("new_primitive_cocycle_order_checks", result.get("new_primitive_cocycle_order_checks", 0))
    audit.setdefault("new_primitive_cocycle_subresiduals", result.get("new_primitive_cocycle_subresiduals", 0))
    audit.setdefault("new_density_proof_rows", 0)
    for key in (
        "old_generic_identity_calls",
        "old_provider_calls",
        "old_checker_calls",
        "old_normalizer_calls",
        "old_source_value_reducer_calls",
        "old_operator_calls",
        "old_source_consumer_calls",
        "source_value_reducer_calls",
        "numeric_field_evaluations",
        "old_numeric_calls",
        "numeric_calls",
        "old_math_calls",
        "old_source_calls",
        "numeric_values_computed",
    ):
        audit.setdefault(key, 0)
    result["goal_active"] = True
    result["full_goal_complete"] = False
    return result


def _result_count(result, names):
    nodes = [result]
    for key in (
        "counts",
        "new_counts",
        "cocycle",
        "transition_driver_cocycle",
        "driver_cocycle",
        "primitive_gap",
        "transition_primitive_gap",
        "root_shape",
        "source_gap",
        "source_inspection_audit",
        "audit",
    ):
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
    """Validate durable primitive-gap metadata without recomputing mathematics."""
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
    require(result.get("numerical") is not True, "O3 numerical mode is not false")
    require(result.get("global_closure") is not True, "O3 global closure is not false")
    expected_algebra_calls = 1 if result.get("passed") is True else 0
    require(result.get("new_algebra_family_calls") == expected_algebra_calls, "O3 algebra-family call count is not truthful")

    forbidden_calls = (
        "old_generic_identity_calls",
        "old_provider_calls",
        "old_checker_calls",
        "old_normalizer_calls",
        "old_source_value_reducer_calls",
        "old_operator_calls",
        "old_source_consumer_calls",
        "source_value_reducer_calls",
        "numeric_field_evaluations",
        "old_numeric_calls",
        "numeric_calls",
        "old_math_calls",
        "old_source_calls",
    )
    for key in forbidden_calls:
        require(result.get(key, 0) == 0, "O3 result recorded forbidden old/math/numeric work: " + key)
    for key in NONCLAIM_FLAGS:
        require(result.get(key) is not True, "O3 overclaims nonclaim flag: " + key)
    for key in ("E_zero_claimed", "A_zero_claimed", "B_zero_claimed", "E_gap_zero_claimed", "A_gap_zero_claimed", "B_gap_zero_claimed", "P2_zero_claimed", "phase_equality_claimed", "phase_equality_proved", "signed_gap_zero_claimed", "gap_zero_claimed", "primitive_gap_zero_claimed"):
        require(result.get(key) is not True, "O3 overclaims forbidden primitive-gap claim: " + key)

    audit = result.get("source_inspection_audit", {})
    require(isinstance(audit, dict), "O3 source audit is not a mapping")
    require(audit.get("original_source_evaluations", 0) == 0, "O3 source audit evaluated an original source")
    require(audit.get("new_algebra_family_calls") == expected_algebra_calls, "O3 audit algebra-family call count changed")
    for key in forbidden_calls + ("numeric_values_computed",):
        require(audit.get(key, 0) in (0, False), "O3 audit recorded forbidden old/math/numeric work: " + key)
    require(audit.get("new_density_proof_rows", 0) == 0, "O3 recorded a new density proof row")
    require(result.get("new_counts", {}).get("new_density_proof_rows", 0) == 0, "O3 result introduced a density proof row")

    producer_claim = result.get("passed") is True or result.get(NEW_FLAGS[1]) is True
    if producer_claim:
        require(audit.get("new_root_shape_order_checks") == ROOT_SHAPE_ROW_COUNT, "O3 root-shape check rows are not 4")
        require(audit.get("new_root_shape_order_failures", 0) == ROOT_SHAPE_ORDER_FAILURE_COUNT, "O3 root-shape failures are not zero")

    contract = source.get("count_contract", {})
    primitive_count = _result_count(
        result,
        (
            "actual_primitive_gap_functions",
            "signed_primitive_gap_function_count",
            "primitive_gap_function_count",
            "actual_EAB_gap_functions",
            "new_actual_primitive_gap_rows",
        ),
    )
    root_rows = _result_count(
        result,
        ("root_shape_ordinary_Z_rows", "new_root_shape_order_checks", "root_shape_rows", "root_shape_row_count", "ordinary_Z_root_shape_rows"),
    )
    root_subresiduals = _result_count(
        result,
        ("root_shape_subresiduals", "new_root_shape_subresiduals", "root_shape_subresidual_count", "root_shape_identity_residuals"),
    )
    sigma_checks = _result_count(
        result,
        ("new_sigma_identity_checks", "sigma_identity_checks", "sigma_identity_check_count"),
    )
    primitive_cocycle_checks = _result_count(
        result,
        ("new_primitive_cocycle_order_checks", "primitive_cocycle_order_checks", "primitive_cocycle_check_count"),
    )
    primitive_cocycle_subresiduals = _result_count(
        result,
        ("new_primitive_cocycle_subresiduals", "primitive_cocycle_subresiduals", "primitive_cocycle_residual_count"),
    )
    if result.get("passed") is False:
        require(not any(result.get(key) is True for key in NEW_FLAGS), "Failed O3 result overclaims a new flag")
        return {
            "passed": False,
            "honest_failure": True,
            "actual_primitive_gap_functions": primitive_count,
            "root_shape_ordinary_Z_rows": root_rows,
            "root_shape_subresiduals": root_subresiduals,
            "new_sigma_identity_checks": sigma_checks,
            "new_primitive_cocycle_order_checks": primitive_cocycle_checks,
            "new_primitive_cocycle_subresiduals": primitive_cocycle_subresiduals,
            "new_density_proof_rows": 0,
            "forbidden_work_calls": 0,
        }

    require(result.get("source_packet_definition_accepted") is True, "Passed O3 result lacks owned source acceptance")
    require(result.get("source_packet_accepted_owned") is True, "Passed O3 result lacks owned acceptance marker")
    require(primitive_count == contract.get("primitive_gap_functions"), "O3 primitive gap function count changed")
    for family in ("E", "A", "B"):
        family_count = _result_count(result, ("actual_" + family + "_gap_functions", family + "_gap_function_count"))
        require(family_count == contract.get(family + "_gap_functions"), "O3 " + family + " primitive gap function count changed")
    require(root_rows == contract.get("root_shape_rows"), "O3 root-shape ordinary-Z row count changed")
    require(root_subresiduals == contract.get("root_shape_subresiduals"), "O3 root-shape subresidual count changed")
    require(sigma_checks == contract.get("sigma_identity_checks"), "O3 sigma identity check count changed")
    require(primitive_cocycle_checks == contract.get("primitive_cocycle_order_checks"), "O3 primitive A/B cocycle order count changed")
    require(primitive_cocycle_subresiduals == contract.get("primitive_cocycle_subresiduals"), "O3 primitive A/B cocycle residual count changed")
    sigma_contract = source.get("sigma_definition_contract", {})
    sigma_exact = result.get("root_shape_sigma_same_exact_definition")
    if sigma_exact is None:
        sigma_exact = result.get("sigma_same_exact_definition")
    require(sigma_exact is True, "O3 root-shape sigma definition is not exact")
    interval_equality = result.get("root_shape_sigma_interval_evaluated_value_equality")
    if interval_equality is None:
        interval_equality = result.get("sigma_interval_evaluated_value_equality", False)
    require(interval_equality is not True, "O3 root-shape sigma used interval value equality")
    require(sigma_contract.get("root_shape_sigma_same_exact_definition") is True, "O3 sigma source contract changed")
    identity_fields = result.get("root_shape_identity_fields", ["a", "b", "t0", "Delta"])
    require(identity_fields == ["a", "b", "t0", "Delta"], "O3 root-shape identity field contract changed")
    sigma_families = result.get("sigma_identity_families", sigma_contract.get("crossmultiply_identity_families"))
    require(sigma_families == sigma_contract.get("crossmultiply_identity_families"), "O3 sigma identity family contract changed")
    require(all(result.get(key) is True for key in NEW_FLAGS), "O3 flags are not all source-bound")
    audit_root_rows = audit.get("new_root_shape_order_checks")
    require(audit_root_rows == contract.get("root_shape_rows"), "O3 ordinary-Z root-shape row count changed")
    require(audit.get("new_root_shape_order_failures", 0) == ROOT_SHAPE_ORDER_FAILURE_COUNT, "O3 root-shape failure count changed")
    audit_gap = audit.get("new_actual_primitive_gap_rows")
    require(audit_gap == contract.get("primitive_gap_functions"), "O3 audit primitive gap count changed")
    require(audit.get("root_shape_ordinary_Z_rows", root_rows) == contract.get("root_shape_rows"), "O3 audit root-shape row count changed")
    require(audit.get("root_shape_subresiduals", root_subresiduals) == contract.get("root_shape_subresiduals"), "O3 audit root-shape residual count changed")
    return {
        "passed": True,
        "actual_primitive_gap_functions": primitive_count,
        "actual_E_gap_functions": _result_count(result, ("actual_E_gap_functions", "E_gap_function_count")),
        "actual_A_gap_functions": _result_count(result, ("actual_A_gap_functions", "A_gap_function_count")),
        "actual_B_gap_functions": _result_count(result, ("actual_B_gap_functions", "B_gap_function_count")),
        "root_shape_ordinary_Z_rows": root_rows,
        "root_shape_checks_per_row": contract.get("root_shape_checks_per_row"),
        "root_shape_subresiduals": root_subresiduals,
        "new_sigma_identity_checks": sigma_checks,
        "new_primitive_cocycle_order_checks": primitive_cocycle_checks,
        "new_primitive_cocycle_subresiduals": primitive_cocycle_subresiduals,
        "primitive_cocycle_residuals_per_row": contract.get("primitive_cocycle_residuals_per_row"),
        "new_root_shape_order_checks": audit.get("new_root_shape_order_checks"),
        "new_root_shape_order_failures": audit.get("new_root_shape_order_failures", 0),
        "new_actual_primitive_gap_rows": audit.get("new_actual_primitive_gap_rows"),
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
        spec = importlib.util.spec_from_file_location("_compact_O3_transition_primitive_gap_entry", kernel)
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


def validate_result_checkpoint_links(source_raw, result_raw, result):
    """Authenticate every saved result-to-checkpoint digest before publication."""
    entered_raw = read_bytes(ENTERED_CHECKPOINT)
    diagnostic_raw = read_bytes(DIAGNOSTIC_CHECKPOINT)
    require(result.get("source_checkpoint_sha256") == sha_bytes(source_raw), "Saved result source checkpoint SHA changed")
    require(result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw), "Saved result entered checkpoint SHA changed")
    require(result.get("diagnostic_checkpoint_sha256") == sha_bytes(diagnostic_raw), "Saved result diagnostic checkpoint SHA changed")
    return {
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entry_checkpoint_sha256": sha_bytes(entered_raw),
        "diagnostic_checkpoint_sha256": sha_bytes(diagnostic_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
    }


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

    # Repair stale aliases in the shared brief.  The canonical AGENTS path,
    # new milestone, and new report are always recorded in read_first; a
    # last_turn progress record is changed only when the pointer explicitly
    # marks that progress as accepted.
    canonical_read_first = [
        wire_path(REPO / "AGENTS.md"),
        wire_path(REPO / MILESTONE_RELATIVE),
        wire_path(REPO / REPORT_RELATIVE),
    ]
    read_first = state.get("read_first")
    if not isinstance(read_first, list):
        read_first = []
    state["read_first"] = canonical_read_first + [
        item for item in read_first if item not in canonical_read_first
    ]
    state["read_first_alias_repairs"] = {
        "AGENTS.md": wire_path(REPO / "AGENTS.md"),
        "newdoc": wire_path(REPO / MILESTONE_RELATIVE),
        "report": wire_path(REPO / REPORT_RELATIVE),
    }
    if pointer.get("accepted_progress") is True and isinstance(pointer.get("last_turn"), dict):
        state["last_turn"] = clone(pointer["last_turn"])

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
            state["read_first"] = canonical_read_first + [
                item for item in read_first if item not in canonical_read_first
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
        "actual_primitive_gap_functions": report.get("actual_primitive_gap_functions"),
        "root_shape_ordinary_Z_rows": report.get("root_shape_ordinary_Z_rows"),
        "root_shape_subresiduals": report.get("root_shape_subresiduals"),
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
        "parent_head": source["parent_head"],
        "parent_report_parent_commit": source["parent_report_parent_commit"],
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
        "preserved_error_cap_metadata": clone(source.get("preserved_error_cap_metadata", {})),
        "inherited_code2": clone(source.get("inherited_code2", [])),
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
        "new_algebra_family_calls": result.get("new_algebra_family_calls", 0),
        "new_root_shape_order_checks": source["count_contract"]["root_shape_rows"],
        "new_root_shape_order_failures": ROOT_SHAPE_ORDER_FAILURE_COUNT,
        "new_actual_primitive_gap_rows": source["count_contract"]["primitive_gap_functions"],
        "actual_primitive_gap_functions": source["count_contract"]["primitive_gap_functions"],
        "actual_E_gap_functions": source["count_contract"]["E_gap_functions"],
        "actual_A_gap_functions": source["count_contract"]["A_gap_functions"],
        "actual_B_gap_functions": source["count_contract"]["B_gap_functions"],
        "root_shape_ordinary_Z_rows": source["count_contract"]["root_shape_rows"],
        "root_shape_checks_per_row": source["count_contract"]["root_shape_checks_per_row"],
        "root_shape_subresiduals": source["count_contract"]["root_shape_subresiduals"],
        "new_sigma_identity_checks": source["count_contract"]["sigma_identity_checks"],
        "new_primitive_cocycle_order_checks": source["count_contract"]["primitive_cocycle_order_checks"],
        "new_primitive_cocycle_subresiduals": source["count_contract"]["primitive_cocycle_subresiduals"],
        "primitive_cocycle_residuals_per_row": source["count_contract"]["primitive_cocycle_residuals_per_row"],
        "root_shape_sigma_same_exact_definition": True,
        "root_shape_sigma_interval_evaluated_value_equality": False,
        "root_shape_identity_fields": ["a", "b", "t0", "Delta"],
        "sigma_identity_families": [
            "positive_interior_logistic_reflection",
            "source_endpoint_binding",
        ],
        "new_density_proof_rows": 0,
        "derivation_contract": clone(source["derivation_contract"]),
        "counts": {**clone(source["counts"]), "new_algebra_family_calls": result.get("new_algebra_family_calls", 0)},
        "new_counts": {**clone(source["new_counts"]), "new_algebra_family_calls": result.get("new_algebra_family_calls", 0)},
        "source_inspection_audit": {
            "new_algebra_family_calls": result.get("new_algebra_family_calls", 0),
            "new_root_shape_order_checks": source["count_contract"]["root_shape_rows"],
            "new_root_shape_order_failures": ROOT_SHAPE_ORDER_FAILURE_COUNT,
            "new_actual_primitive_gap_rows": source["count_contract"]["primitive_gap_functions"],
            "actual_primitive_gap_functions": source["count_contract"]["primitive_gap_functions"],
            "actual_E_gap_functions": source["count_contract"]["E_gap_functions"],
            "actual_A_gap_functions": source["count_contract"]["A_gap_functions"],
            "actual_B_gap_functions": source["count_contract"]["B_gap_functions"],
            "root_shape_ordinary_Z_rows": source["count_contract"]["root_shape_rows"],
            "root_shape_checks_per_row": source["count_contract"]["root_shape_checks_per_row"],
            "root_shape_subresiduals": source["count_contract"]["root_shape_subresiduals"],
            "new_sigma_identity_checks": source["count_contract"]["sigma_identity_checks"],
            "new_primitive_cocycle_order_checks": source["count_contract"]["primitive_cocycle_order_checks"],
            "new_primitive_cocycle_subresiduals": source["count_contract"]["primitive_cocycle_subresiduals"],
            "primitive_cocycle_residuals_per_row": source["count_contract"]["primitive_cocycle_residuals_per_row"],
            "new_root_shape_order_checks": source["count_contract"]["root_shape_rows"],
            "new_actual_primitive_gap_rows": source["count_contract"]["primitive_gap_functions"],
            "new_density_proof_rows": 0,
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
        "old_operator_calls": 0,
        "old_source_consumer_calls": 0,
        "source_value_reducer_calls": 0,
        "numeric_field_evaluations": 0,
        "old_provider_calls": 0,
        "old_checker_calls": 0,
        "old_normalizer_calls": 0,
        "old_source_value_reducer_calls": 0,
        "numeric_calls": 0,
        "old_math_calls": 0,
        "old_source_calls": 0,
        "numeric": False,
        "numerical": False,
        "E_zero_claimed": False,
        "A_zero_claimed": False,
        "B_zero_claimed": False,
        "E_gap_zero_claimed": False,
        "A_gap_zero_claimed": False,
        "B_gap_zero_claimed": False,
        "P2_zero_claimed": False,
        "phase_equality_claimed": False,
        "phase_equality_proved": False,
        "signed_gap_zero_claimed": False,
        "gap_zero_claimed": False,
        "primitive_gap_zero_claimed": False,
        "global_closure": False,
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
    require(report.get("parent_head") == PARENT_COMMIT, "Published parent HEAD changed")
    require(report.get("parent_report_parent_commit") == PARENT_REPORT_PARENT_COMMIT, "Published parent report parent commit changed")
    require(report.get("parent_result_raw_sha256") == source.get("parent_result_raw_sha256"), "Published parent result SHA changed")
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Published checkpoint count changed")
    contract = source.get("count_contract", {})
    expected_algebra_calls = 1 if report.get("passed") is True else 0
    require(report.get("new_algebra_family_calls") == expected_algebra_calls, "Published algebra-family call count changed")
    require(report.get("counts", {}).get("new_algebra_family_calls") == expected_algebra_calls, "Published algebra-family count map changed")
    require(report.get("new_counts", {}).get("new_algebra_family_calls") == expected_algebra_calls, "Published new algebra-family count map changed")
    require(report.get("actual_primitive_gap_functions") == contract.get("primitive_gap_functions"), "Published primitive gap count changed")
    require(report.get("new_actual_primitive_gap_rows") == contract.get("primitive_gap_functions"), "Published primitive gap row count changed")
    require(report.get("new_root_shape_order_checks") == contract.get("root_shape_rows"), "Published root-shape order count changed")
    require(report.get("new_root_shape_order_failures", 0) == ROOT_SHAPE_ORDER_FAILURE_COUNT, "Published root-shape order failures changed")
    require(report.get("counts", {}).get("new_root_shape_order_checks") == contract.get("root_shape_rows"), "Published root-shape count map changed")
    require(report.get("new_counts", {}).get("new_root_shape_order_checks") == contract.get("root_shape_rows"), "Published new root-shape count map changed")
    for family in ("E", "A", "B"):
        require(report.get("actual_" + family + "_gap_functions") == contract.get(family + "_gap_functions"), "Published " + family + " primitive gap count changed")
    require(report.get("root_shape_ordinary_Z_rows") == contract.get("root_shape_rows"), "Published root-shape row count changed")
    require(report.get("root_shape_subresiduals") == contract.get("root_shape_subresiduals"), "Published root-shape subresidual count changed")
    require(report.get("new_sigma_identity_checks") == contract.get("sigma_identity_checks"), "Published sigma identity count changed")
    require(report.get("new_primitive_cocycle_order_checks") == contract.get("primitive_cocycle_order_checks"), "Published primitive cocycle order count changed")
    require(report.get("new_primitive_cocycle_subresiduals") == contract.get("primitive_cocycle_subresiduals"), "Published primitive cocycle residual count changed")
    require(report.get("primitive_cocycle_residuals_per_row") == contract.get("primitive_cocycle_residuals_per_row"), "Published primitive cocycle residuals-per-row changed")
    require(report.get("root_shape_sigma_same_exact_definition") is True, "Published root-shape sigma definition changed")
    require(report.get("root_shape_sigma_interval_evaluated_value_equality") is not True, "Published root-shape sigma interval equality changed")
    require(
        report.get("sigma_identity_families") == ["positive_interior_logistic_reflection", "source_endpoint_binding"],
        "Published sigma identity families changed",
    )
    audit = report.get("source_inspection_audit", {})
    require(audit.get("new_algebra_family_calls") == expected_algebra_calls, "Published algebra-family audit count changed")
    require(audit.get("new_sigma_identity_checks") == contract.get("sigma_identity_checks"), "Published sigma audit count changed")
    require(audit.get("new_primitive_cocycle_order_checks") == contract.get("primitive_cocycle_order_checks"), "Published primitive cocycle audit count changed")
    require(audit.get("new_root_shape_order_checks") == contract.get("root_shape_rows"), "Published root-shape audit count changed")
    require(audit.get("new_actual_primitive_gap_rows") == contract.get("primitive_gap_functions"), "Published primitive gap audit count changed")
    require(audit.get("new_density_proof_rows") == 0, "Published density audit count changed")
    for key in ("E_zero_claimed", "A_zero_claimed", "B_zero_claimed", "E_gap_zero_claimed", "A_gap_zero_claimed", "B_gap_zero_claimed", "P2_zero_claimed", "phase_equality_claimed", "phase_equality_proved", "signed_gap_zero_claimed", "gap_zero_claimed", "primitive_gap_zero_claimed"):
        require(report.get(key) is not True, "Published forbidden zero/phase claim changed: " + key)
    require(report.get("new_density_proof_rows") == 0, "Published density proof row count changed")
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
    checkpoint_links = validate_result_checkpoint_links(source_raw, result_raw, result)
    durable_validation = validate_result_metadata(result, source)
    report = report_for_publish(source, result, result_raw, code2, input_meta_for_publish)
    report["result_checkpoint_links"] = checkpoint_links
    report["durable_result_validation"] = clone(durable_validation)
    report["publication_metadata_validation"] = validate_publication_metadata(report, source)
    report_raw = write_replace(REPO / REPORT_RELATIVE, report)
    milestone = (
        "# Current original O3 transition primitive gap (2026-10-27)\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}`; SHA-256 `{sha_bytes(report_raw)}`.\n\n"
        f"Attempt status: `{report['attempt_status']}`; source-only: `true`; sigma cross-multiply checks: `{SIGMA_IDENTITY_CHECK_COUNT}/{SIGMA_IDENTITY_CHECK_COUNT}`; primitive A/B cocycle rows: `{PRIMITIVE_COCYCLE_ORDER_CHECK_COUNT}/{PRIMITIVE_COCYCLE_ORDER_CHECK_COUNT}`; root-shape ordinary-Z rows: `{ROOT_SHAPE_ROW_COUNT}/{ROOT_SHAPE_ROW_COUNT}`; root-shape identities: `{ROOT_SHAPE_SUBRESIDUAL_COUNT}`; actual signed E/A/B primitive gap functions: `{ACTUAL_PRIMITIVE_GAP_COUNT}/{ACTUAL_PRIMITIVE_GAP_COUNT}`.\n\n"
        "This leaf binds two sigma cross-multiply identities, four primitive A/B cocycle rows, actual E/A/B primitive source gaps, four ordinary-Z root-shape rows, and sixteen root-shape subresiduals. "
        "It introduces no numeric field values, density proof rows, E/A/B zero claim, P2 zero claim, phase equality, "
        "full-pair proof, incoming identity, cap closure, global closure, physical completion, or full-goal completion. "
        "The old generic identity is retained only as cached data.\n\n"
        "The entered attempt is immutable and must never be reentered.\n\n"
        f"Next source-bound task: P2 `{next_task_metadata()['phase']}` with `{ACTUAL_PRIMITIVE_GAP_COUNT}` primitive source functions; "
        "the source gap remains unknown and is not asserted zero.\n"
    ).encode("utf-8")
    milestone_raw = write_replace(REPO / MILESTONE_RELATIVE, milestone)
    preimages = []
    handoff_prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}` (SHA-256 `{sha_bytes(report_raw)}`).\n"
        f"Status: `{report['attempt_status']}`; sigma checks `{SIGMA_IDENTITY_CHECK_COUNT}/{SIGMA_IDENTITY_CHECK_COUNT}`; primitive cocycle rows `{PRIMITIVE_COCYCLE_ORDER_CHECK_COUNT}/{PRIMITIVE_COCYCLE_ORDER_CHECK_COUNT}`; root-shape rows `{ROOT_SHAPE_ROW_COUNT}/{ROOT_SHAPE_ROW_COUNT}`; primitive gap functions `{ACTUAL_PRIMITIVE_GAP_COUNT}/{ACTUAL_PRIMITIVE_GAP_COUNT}`; full goal remains incomplete.\n"
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
        f"/{KERNEL_RELATIVE.as_posix()} -text",
        f"/{RUNNER_RELATIVE.as_posix()} -text",
        f"/{REPORT_RELATIVE.as_posix()} -text",
        f"/{MILESTONE_RELATIVE.as_posix()} -text",
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
            "accepted_progress": report["attempt_status"] == "passed",
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
    require(report.get("parent_head") == PARENT_COMMIT, "Published parent HEAD changed")
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
        committed_blob = git_output("rev-parse", "HEAD:" + relative.as_posix()).lower()
        local_blob = git_output("hash-object", wire_path(local_path)).lower()
        require(valid_sha(committed_blob) and committed_blob == local_blob, "Published bytes differ from committed scope: " + relative.as_posix())
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
            "accepted_progress": report.get("passed") is True,
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
    checkpoint_links = validate_result_checkpoint_links(source_raw, result_raw, result)
    require(report.get("result_checkpoint_links") == checkpoint_links, "Published result checkpoint links changed before seal")
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
            "accepted_progress": True,
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





