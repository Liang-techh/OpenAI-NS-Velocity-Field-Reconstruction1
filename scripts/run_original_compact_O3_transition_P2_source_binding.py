"""Metadata-only manager for the independent O3 transition P2 source binding.

This file is intentionally a workflow draft.  It has no mathematical proof
kernel and does not replay the failed primitive, sigma, root-shape, cocycle,
source-provider, checker, normalizer, reducer, or numeric paths.  The root
supplies the input, kernel, and runner bytes before ``prepare``.  ``prepare``
freezes those bytes and an authenticated failed-parent projection.  ``run``
enters once, fsyncs the entered marker before importing only the new kernel,
invokes ``derive`` once, journals each callback immediately, and fsyncs the
result before metadata validation.  The result is never rewritten.  Publish,
resume, and seal consume saved bytes and metadata only.

The planned source-port contract is kept separate from actual durable result
counts.  A failed attempt therefore reports only rows saved by the kernel,
never the planned eight rows.  New mathematical proofs, new algebra-family
calls, numeric calls, source-provider calls, and old-math calls are always
zero for this binding-only leaf.
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

COMPONENT = "compact_O3_transition_P2_source_binding"
INPUT_SCHEMA = "compact_O3_transition_P2_source_binding_input_v1"
RESULT_SCHEMA = "original_compact_O3_transition_P2_source_binding"
ENTERED_SCHEMA = "original_compact_O3_transition_P2_source_binding_entered_v1"
DIAGNOSTIC_SCHEMA = "compact_O3_transition_P2_source_binding_stage_diagnostic_v1"
PUBLICATION_SCHEMA = "original_compact_O3_transition_P2_source_binding_publication_v1"
SEAL_SCHEMA = "compact_O3_transition_P2_source_binding_final_seal_v1"

KERNEL_RELATIVE = Path(
    "src/openai_ns_reconstruction/"
    "original_compact_O3_transition_P2_source_binding.py"
)
RUNNER_RELATIVE = Path(
    "scripts/run_original_compact_O3_transition_P2_source_binding.py"
)
REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_P2_source_binding.json"
)
MILESTONE_RELATIVE = Path(
    "docs/CURRENT_ORIGINAL_COMPACT_O3_TRANSITION_P2_SOURCE_BINDING_2026_10_27.md"
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
# Match the repository's established binary/source artifact convention: the
# four science artifacts receive ``-text``; the other five exact-scope files
# remain ordinary text handoffs or attribute metadata.
SCOPE_ATTRIBUTE_RELATIVES = (
    KERNEL_RELATIVE,
    RUNNER_RELATIVE,
    REPORT_RELATIVE,
    MILESTONE_RELATIVE,
)
SCOPE_ATTRIBUTE_LINES = tuple(
    f"/{relative.as_posix()} -text" for relative in SCOPE_ATTRIBUTE_RELATIVES
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

PARENT_REPORT_DEFAULT = REPO / (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_primitive_gap.json"
)
PARENT_REPORT_SHA256 = (
    "3261bb535e3987c7e74c02cef84871834f6fa2bdf894b0bc0c907f08095d1cd7"
)
PARENT_RESULT_SHA256 = (
    "a257f372f3291ecc921ab8fe0d42609e375663b6e8a0f67b713a154b3b0f7525"
)
# This is the current parent HEAD; the parent report itself records the prior
# commit 098c4687... in its own parent_commit/parent_head fields.
PARENT_HEAD = "45bad2f3d7bd2e98ab3cee137f64a3f5a390c289"
# Compatibility alias used by the earlier metadata managers; this component
# treats it as the authenticated parent HEAD, not the report's internal
# parent_commit field.
PARENT_COMMIT = PARENT_HEAD
PARENT_REPORT_COMMIT = "098c46878efd98c84ec2ef438c9d7f1596f529d4"
PARENT_REPORT_PARENT_COMMIT = "b0f7e8e0d9a69b8c82345ff0a6a356c8b7a2cc96"

PARENT_SOURCE_PIN_COUNT = 2219
PARENT_CHECKPOINT_COUNT = 244
CURRENT_CHECKPOINT_COUNT = PARENT_CHECKPOINT_COUNT + 4
PARENT_FAILED_COUNT = 176
FAILED_COUNT_IF_NEW_FAILURE = PARENT_FAILED_COUNT + 7
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26

PLANNING_CONTRACT = {
    "source_port_rows": 8,
    "native_source_port_rows": 4,
    "C3_source_port_rows": 4,
    "projection_rows": 8,
    "projection_orders": [0, 1, 2, 3],
    "pressure_route_sides": 2,
    "new_mathematical_proofs": 0,
    "new_algebra_family_calls": 0,
    "numeric_calls": 0,
    "source_provider_calls": 0,
    "old_math_calls": 0,
}
ACTUAL_COUNT_FIELDS = (
    "source_port_rows",
    "native_source_port_rows",
    "C3_source_port_rows",
    "projection_rows",
    "pressure_route_sides",
    "new_mathematical_proofs",
    "new_algebra_family_calls",
    "numeric_calls",
    "source_provider_calls",
    "old_math_calls",
)
ZERO_ACTUAL_COUNTS = {key: 0 for key in ACTUAL_COUNT_FIELDS}

NEW_FLAGS = (
    "O3_transition_P2_source_ports_bound",
    "O3_transition_P2_pressure_restoration_route_bound",
    "O3_transition_P2_ordinary_Z_source_projection_bound",
)
EXPECTED_NATIVE_P2_ROOT_IDS = (242304, 242308, 242314, 242321)
EXPECTED_C3_P2_ROOT_IDS = (4598, 4599, 9992, 21729)
EXPECTED_PROJECTION_ORDERS = (0, 1, 2, 3)
# These names are inspected only in saved result/audit structures.  Any
# nonzero value is a protocol violation for this source-binding-only leaf.
FORBIDDEN_WORK_FIELD_ALIASES = {
    "new_mathematical_proofs": ("new_mathematical_proofs", "actual_new_mathematical_proofs"),
    "new_algebra_family_calls": ("new_algebra_family_calls", "actual_new_algebra_family_calls"),
    "numeric_calls": ("numeric_calls", "actual_numeric_calls"),
    "source_provider_calls": ("source_provider_calls", "actual_source_provider_calls"),
    "old_math_calls": ("old_math_calls", "actual_old_math_calls"),
    "provider_calls": ("provider_calls", "actual_provider_calls"),
    "checker_calls": ("checker_calls", "actual_checker_calls"),
    "normalizer_calls": ("normalizer_calls", "actual_normalizer_calls"),
    "reducer_calls": ("reducer_calls", "actual_reducer_calls"),
    "quadrature_calls": ("quadrature_calls", "actual_quadrature_calls"),
    "recurrence_calls": ("recurrence_calls", "actual_recurrence_calls"),
    "field_calls": ("field_calls", "actual_field_calls"),
    "field_evaluations": ("field_evaluations", "actual_field_evaluations"),
    "integral_calls": ("integral_calls", "actual_integral_calls"),
    "integral_evaluations": ("integral_evaluations", "actual_integral_evaluations"),
    "proof_calls": ("proof_calls", "actual_proof_calls"),
    "math_calls": ("math_calls", "actual_math_calls"),
}
OLD_QUIET_FLAGS = (
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
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
        "kernel_imported",
        "derive_entry_calls",
        "derive_calls",
        "entry_attempts",
        "source_checkpoint_sha256",
        "entry_checkpoint_sha256",
        "diagnostic_checkpoint_sha256",
        "new_flags",
        "actual_counts",
        "planning_contract",
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
    "parent_report_view",
    "parent_report",
    "parent_checkpoint_inventory",
    "checkpoint_hashes",
    "checkpoint_paths",
    "protected_file_hashes",
    "legacy_dirty_path_hashes",
    "preserved_failed_artifact_hashes_union",
    "implementation_hashes",
    "code2",
    "input_artifact",
    "metadata_helper_artifact",
    "metadata_repair_artifact",
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
    return Path(path).as_posix()


def portable_metadata(value):
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
        require(sha_bytes(raw) == str(expected).lower(), "Changed artifact: " + wire_path(path))
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
    # Git object/revision IDs are SHA-1 here: 40 hexadecimal characters.
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
        require(isinstance(name, str) and valid_sha(digest), label + " has a malformed SHA")
        digest = digest.lower()
        require(name not in result or result[name] == digest, "Conflicting hash at " + name)
        result[name] = digest

    def visit(node, location):
        require(isinstance(node, dict), location + " is not a mapping")
        if not node:
            return
        if all(isinstance(key, str) and valid_sha(value) for key, value in node.items()):
            for key, value in node.items():
                add(key, value)
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
    # All original bytes remain in the committed parent report; this is a pointer view.
    keys = ("schema", "component", "attempt_status", "passed", "source_only",
            "source_pin_count", "source_union_count", "checkpoint_count",
            "preserved_failed_artifact_count", "new_flags", "actual_counts",
            "planned_count_contract", "full_goal_complete", "global_five_moment_closure",
            "retained_construction_data", "retained_integral_term_inventory")
    view = {key: clone(report[key]) for key in keys if key in report}
    view["metadata_only"] = True
    view["complete_parent_report_raw_sha256"] = PARENT_REPORT_SHA256
    view["selection"] = "bounded metadata only; full flags, caps, checkpoint and failed-artifact maps are retained separately from the authenticated parent"
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
    if report.get("full_result_raw_sha256") == report_raw_sha:
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
    parent_commit = (args.parent_commit or PARENT_HEAD).lower()
    require(report_sha == PARENT_REPORT_SHA256, "Parent report SHA is not the authenticated failed projection")
    require(parent_commit == PARENT_HEAD, "Parent HEAD is not the authenticated failed projection parent")
    report_raw, report = read_json(report_path, report_sha)
    require(isinstance(report, dict), "Parent report is not a mapping")
    require(report.get("schema") == "original_compact_O3_transition_primitive_gap_publication_v1", "Parent report schema changed")
    require(report.get("component") == "compact_O3_transition_primitive_gap", "Parent report component changed")
    require(report.get("attempt_status") == "failed" and report.get("passed") is False, "Parent attempt is not the saved failure")
    require(report.get("source_only") is True, "Parent report is not source-only")
    require(report.get("parent_commit") == PARENT_REPORT_COMMIT, "Parent report internal commit changed")
    require(report.get("parent_head") == PARENT_REPORT_COMMIT, "Parent report internal HEAD changed")
    require(report.get("parent_report_parent_commit") == PARENT_REPORT_PARENT_COMMIT, "Parent report ancestry changed")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT, "Parent source-pin count changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Parent checkpoint count changed")
    require(report.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT, "Parent failed count changed")
    require(len(report.get("false_flags", [])) == PARENT_FALSE_FLAG_COUNT, "Parent false-flag count changed")
    require(len(report.get("protected_file_hashes", {})) == PARENT_PROTECTED_COUNT, "Parent protected-file count changed")
    require(len(report.get("legacy_dirty_path_hashes", {})) == PARENT_DIRTY_COUNT, "Parent dirty-path count changed")
    false_flags = report.get("false_flags")
    require(isinstance(false_flags, list) and len(set(false_flags)) == PARENT_FALSE_FLAG_COUNT, "Parent false flags are malformed")
    for key in false_flags:
        require(report.get(key) is False, "Parent false flag changed: " + str(key))

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
    require(result_ref["sha256"] == PARENT_RESULT_SHA256, "Parent result SHA changed")
    require(result.get("component") == "compact_O3_transition_primitive_gap", "Parent result component changed")
    require(result.get("passed") is False, "Parent saved result is not the failed result")
    require(report.get("failure_actual_counts") == report.get("actual_counts"), "Parent failure counts are not actual counts")
    require(all(value == 0 for value in report.get("failure_actual_counts", {}).values()), "Parent failure actual counts are nonzero")

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
    preserved_keys = (
        "retained_construction_data",
        "retained_integral_term_inventory",
        "retained_cap_provenance",
        "preserved_error_cap_metadata",
        "inherited_failed_artifacts_preserved",
        "old_errors_removed",
        "old_error_terms_removed",
    )
    preserved = {key: clone(report[key]) for key in preserved_keys if key in report}
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
        "parent_head": parent_commit,
        "parent_report_commit": report.get("parent_commit"),
        "parent_report_parent_commit": report.get("parent_report_parent_commit"),
        "preserved": preserved,
    }


def discover_receipts(value, label):
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
    require(isinstance(input_value, dict), "P2 input is not a mapping")
    require(input_value.get("schema") == INPUT_SCHEMA, "P2 input schema changed")
    # The source packet owns the frozen planning contract.  The current input
    # is a source-only extraction packet and therefore may omit that metadata;
    # if it supplies a contract, authenticate only the known fields without
    # requiring the packet to duplicate manager-owned planning metadata.
    declared = input_value.get("planning_contract")
    if declared is not None:
        require(isinstance(declared, dict), "Input planning contract is not a mapping")
        require(set(declared).issubset(PLANNING_CONTRACT), "Input planning contract has unknown fields")
        for key, expected in declared.items():
            require(expected == PLANNING_CONTRACT[key], "Input planning contract changed: " + key)

    # ``additional_raw_pins`` is source metadata in the current packet (root
    # IDs, formulas, audit booleans, and similar values), not a path-to-SHA
    # map.  Only an explicitly named raw_sha_pins map, or a legacy map whose
    # values are all SHA-256 strings, is authenticated as extra file input.
    additional = input_value.get("additional_raw_pins", {})
    require(isinstance(additional, dict), "Input additional raw pins is not a mapping")
    declared_pins = input_value.get("raw_sha_pins", {})
    require(isinstance(declared_pins, dict), "Input raw SHA pins is not a mapping")
    if not declared_pins and additional and all(valid_sha(digest) for digest in additional.values()):
        declared_pins = additional

    raw_pins = {}
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
            {"path": wire_path(resolved), "sha256": receipt["sha256"], "location": receipt["location"]}
        )
    packet_pins = merge_hashes(raw_pins, {item["path"]: item["sha256"] for item in portable_receipts})
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
    require(path == REPO / KERNEL_RELATIVE, "P2 kernel path is outside the frozen source path")
    raw = read_bytes(path)
    if args.kernel_sha256:
        require(valid_sha(args.kernel_sha256), "Kernel SHA-256 is malformed")
        require(sha_bytes(raw) == args.kernel_sha256.lower(), "Kernel SHA-256 changed")
    return path, raw


def ensure_runner(args):
    path = resolve_path(args.runner_path or wire_path(REPO / RUNNER_RELATIVE), REPO)
    require(path == REPO / RUNNER_RELATIVE, "P2 runner path is outside the frozen source path")
    require(path.exists(), "P2 runner is missing; root must supply it before prepare")
    raw = read_bytes(path)
    if args.runner_sha256:
        require(valid_sha(args.runner_sha256), "Runner SHA-256 is malformed")
        require(sha_bytes(raw) == args.runner_sha256.lower(), "Runner SHA-256 changed")
    return path, raw


def next_task_metadata(passed=None):
    if passed is True:
        return {
            "schema": "O3_transition_P2_pressure_gap_operator_next_v1",
            "phase": "O3_transition_P2_pressure_gap_operator",
            "priority": "P2",
            "source_only": True,
            "independent_of_source_binding": True,
            "fresh_operator_required": True,
            "action": "Apply a fresh independent P2 pressure-gap operator to the accepted source-binding data.",
            "task": "Derive the physical pressure restoration gap only after this source-binding result is accepted.",
            "do_not_reenter_current_attempt": True,
            "do_not_replay_failed_primitive": True,
        }
    if passed is False:
        return {
            "schema": "O3_transition_P2_source_binding_failure_inspection_next_v1",
            "phase": "O3_transition_P2_source_binding_failure_inspection",
            "priority": "P2",
            "source_only": True,
            "cached_binding_prefix_only": True,
            "action": "Inspect the cached P2 source-binding prefix and diagnose the saved failure without reentry.",
            "task": "Use only durably saved binding data; do not retry or enter the failed attempt.",
            "do_not_reenter_current_attempt": True,
            "do_not_retry_current_attempt": True,
        }
    return {
        "schema": "O3_transition_P2_source_binding_next_v1",
        "phase": "O3_transition_P2_source_binding",
        "priority": "P2",
        "source_only": True,
        "action": "Bind native and rebuilt C3 P2 source ports before any pressure-gap algebra.",
        "task": "Freeze source ports, pressure-restoration route sides, and ordinary-Z projection orders 0 through 3.",
        "do_not_reenter_current_attempt": True,
        "planned_contract": clone(PLANNING_CONTRACT),
    }


def source_paths(parent, input_meta, kernel, runner):
    # The external manager and post-repair artifacts are metadata helpers, not
    # science pins.  Their hashes are retained separately in raw_sha_pins.
    stage = {
        wire_path(parent["report_meta"]["path"]): parent["report_meta"]["sha256"],
        wire_path(parent["result_meta"]["path"]): parent["result_meta"]["sha256"],
        wire_path(input_meta["path"]): input_meta["sha256"],
        wire_path(kernel): sha_path(kernel),
        wire_path(runner): sha_path(runner),
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
    external = Path(__file__).resolve()
    source = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_only": True,
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "parent_head": parent["parent_head"],
        "parent_commit": parent["parent_head"],
        "parent_report_commit": parent["parent_report_commit"],
        "parent_report_parent_commit": parent["parent_report_parent_commit"],
        "parent_report": clone(parent["report_view"]),
        "parent_report_metadata": clone(parent["report_meta"]),
        "parent_report_sha256": parent["report_meta"]["sha256"],
        "parent_result": clone(parent["result"]),
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
        "planning_contract": clone(PLANNING_CONTRACT),
        "planned_counts_are_not_results": True,
        "new_flags": list(NEW_FLAGS),
        "actual_counts_at_prepare": None,
        "actual_projection_orders_at_prepare": [],
        "source_union_hashes": wire_mapping(source_union),
        "source_union_count": len(source_union),
        "stage_hashes": wire_mapping(stage),
        "current_code_hashes": {
            wire_path(kernel): sha_bytes(kernel_raw),
            wire_path(runner): sha_bytes(runner_raw),
        },
        "kernel_path": wire_path(kernel),
        "runner_path": wire_path(runner),
        "metadata_helper_artifact": {
            "path": wire_path(external),
            "sha256": sha_path(external),
            "bytes": external.stat().st_size,
            "raw_readonly": True,
            "outside_science_pins": True,
        },
        "raw_sha_pins": {
            "parent_report_raw_sha256": parent["report_meta"]["sha256"],
            "parent_result_raw_sha256": parent["result_meta"]["sha256"],
            "input_raw_sha256": input_meta["sha256"],
            "kernel_raw_sha256": sha_bytes(kernel_raw),
            "runner_raw_sha256": sha_bytes(runner_raw),
            "external_manager_raw_sha256": sha_path(external),
            "additional_raw_pins": wire_mapping(input_meta["additional_raw_pins"]),
        },
        "raw_packet_receipts": [clone(item) for item in input_meta["static_receipts"]],
        "inherited_code2": clone(parent["report"].get("code2", [])),
        "preserved_error_cap_metadata": clone(parent["preserved"]),
        "retained_construction_data": clone(parent["preserved"].get("retained_construction_data")),
        "retained_integral_term_inventory": clone(parent["preserved"].get("retained_integral_term_inventory")),
        "retained_cap_provenance": clone(parent["preserved"].get("retained_cap_provenance")),
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
        "old_quiet_flags_explicit_false": {key: False for key in OLD_QUIET_FLAGS},
        "nonclaim_flags_explicit_false": {key: False for key in NONCLAIM_FLAGS},
        "new_mathematical_proofs": 0,
        "new_algebra_family_calls": 0,
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "goal_active": True,
        "full_goal_complete": False,
        "next_executable_task": next_task_metadata(),
        "ancestor_decode_forbidden": True,
        "kernel_import_allowed_only_in_run": True,
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    source["input_digest"] = sha_bytes(encode(source["input"]))
    return source, parent, input_meta, kernel, runner


def validate_metadata_brief(path, raw):
    require(len(raw) <= 100_000, "Metadata brief exceeds 100 KB: " + wire_path(path))
    for marker in (b"raw_bytes_base64", b'"raw_text"', b'"raw_bytes"'):
        require(marker not in raw, "Metadata brief contains embedded bytes: " + wire_path(path))
    return True


def validate_scope_attribute_lines():
    path = REPO / ".gitattributes"
    raw = read_bytes(path)
    text = raw.decode("utf-8", errors="strict")
    lines = text.splitlines()
    for line in SCOPE_ATTRIBUTE_LINES:
        require(lines.count(line) == 1, "P2 .gitattributes scope line is not exact: " + line)
        require(line.endswith(" -text"), "P2 scope attribute is not -text: " + line)
    return list(SCOPE_ATTRIBUTE_LINES)


def update_main_brief(pointer):
    # Reports retain raw bytes; the shared brief receives SHA pointers only.
    def compact(value):
        if isinstance(value, dict):
            return {
                key: compact(child)
                for key, child in value.items()
                if key not in {"raw_bytes_base64", "raw_text", "raw_bytes"}
            }
        if isinstance(value, list):
            return [compact(child) for child in value]
        return value

    pointer = compact(pointer)
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
    canonical_read_first = [
        wire_path(REPO / "AGENTS.md"),
        wire_path(REPO / MILESTONE_RELATIVE),
        wire_path(REPO / REPORT_RELATIVE),
    ]
    # Keep only the current instruction, milestone, and report aliases in the
    # active read-first list; older aliases remain recoverable in the preimage.
    state["read_first"] = canonical_read_first
    state["read_first_alias_repairs"] = {
        "AGENTS.md": wire_path(REPO / "AGENTS.md"),
        "newdoc": wire_path(REPO / MILESTONE_RELATIVE),
        "report": wire_path(REPO / REPORT_RELATIVE),
    }
    if pointer.get("accepted_progress") is True and isinstance(pointer.get("last_turn"), dict):
        state["last_turn"] = clone(pointer["last_turn"])
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
        for key in ("next_remaining", "latest_review", "latest_qa"):
            if pointer.get(key) is not None:
                state[key] = clone(pointer[key])
    merged_raw = encode(state)
    if len(merged_raw) > 100_000:
        compact_keys = {
            "continuation_detail_archive",
            "archived_detail",
            "archived_prior_component_details",
            "failed_source_attempt_evidence",
        }
        compact_state = {key: value for key, value in state.items() if key not in compact_keys}
        old_in_progress = state.get("in_progress", {})
        if isinstance(old_in_progress, dict):
            keep = {
                key for key in old_in_progress
                if key == COMPONENT or "quiet" in key.lower()
            }
            omitted = sorted(key for key in old_in_progress if key not in keep)
            compact_state["in_progress"] = {key: old_in_progress[key] for key in sorted(keep)}
            compact_state["archived_in_progress_selector"] = {
                "selector": "$.in_progress",
                "path": wire_path(preimage),
                "sha256": sha_bytes(old_raw),
                "bytes": len(old_raw),
                "omitted_keys": omitted,
                "preserved_exact_bytes": True,
            }
        compact_state["archived_main_resume_brief_preimage"] = {
            "path": wire_path(preimage),
            "sha256": sha_bytes(old_raw),
            "bytes": len(old_raw),
        }
        compact_state["main_resume_brief_archival_compaction"] = True
        merged_raw = encode(compact_state)
    require(len(merged_raw) <= 100_000, "Shared RESUME_BRIEF remains over 100 KB after archival compaction")
    validate_metadata_brief(MAIN_BRIEF, merged_raw)
    write_replace(MAIN_BRIEF, merged_raw)
    return {
        "path": wire_path(MAIN_BRIEF),
        "sha256": sha_bytes(merged_raw),
        "bytes": len(merged_raw),
        "preimage": {"path": wire_path(preimage), "sha256": sha_path(preimage), "bytes": len(old_raw)},
    }


def prepare(args):
    require(not any(path.exists() for path in CHECKPOINTS), "P2 attempt already has a durable checkpoint")
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
            "planning_contract": clone(PLANNING_CONTRACT),
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
        "planning_contract": clone(PLANNING_CONTRACT),
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
        "planning_contract": clone(PLANNING_CONTRACT),
        "actual_counts": None,
        "next_executable_task": next_task_metadata(),
        "runtime_model_metadata": None,
    }
    brief_raw = write_replace(BRIEF, brief)
    validate_metadata_brief(BRIEF, brief_raw)
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
    _prepare_raw, prepare_receipt = read_json(PREPARE_RECEIPT)
    require(prepare_receipt.get("schema") == COMPONENT + "_prepare_receipt_v1", "P2 prepare receipt schema changed")
    require(prepare_receipt.get("component") == COMPONENT, "P2 prepare receipt component changed")
    require(prepare_receipt.get("source_checkpoint") == wire_path(SOURCE_CHECKPOINT), "P2 prepare source path changed")
    require(prepare_receipt.get("source_checkpoint_sha256") == sha_bytes(source_raw), "P2 source SHA is not authenticated")
    require(source.get("schema") == INPUT_SCHEMA and source.get("component") == COMPONENT, "P2 source schema changed")
    require(source.get("metadata_only") is True and source.get("source_only") is True, "P2 source scope changed")
    require(source.get("runtime_model_metadata") is None, "P2 source runtime metadata is not null")
    require(source.get("parent_report_sha256") == PARENT_REPORT_SHA256, "P2 parent report SHA changed")
    require(source.get("parent_result_raw_sha256") == PARENT_RESULT_SHA256, "P2 parent result SHA changed")
    require(source.get("parent_commit") == PARENT_HEAD and source.get("parent_head") == PARENT_HEAD, "P2 parent HEAD changed")
    require(source.get("parent_report_commit") == PARENT_REPORT_COMMIT, "P2 parent report commit changed")
    require(source.get("parent_report_parent_commit") == PARENT_REPORT_PARENT_COMMIT, "P2 parent ancestry changed")
    require(source.get("parent_source_pin_count") == PARENT_SOURCE_PIN_COUNT, "P2 parent pin count changed")
    require(source.get("complete_parent_checkpoint_count") == PARENT_CHECKPOINT_COUNT, "P2 parent checkpoint count changed")
    require(source.get("current_checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "P2 current checkpoint count changed")
    require(source.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT, "P2 parent failed count changed")
    require(len(source.get("protected_file_hashes", {})) == PARENT_PROTECTED_COUNT, "P2 protected count changed")
    require(len(source.get("legacy_dirty_path_hashes", {})) == PARENT_DIRTY_COUNT, "P2 dirty count changed")
    require(source.get("input", {}).get("schema") == INPUT_SCHEMA, "P2 input schema changed")
    require(source.get("planning_contract") == PLANNING_CONTRACT, "P2 planning contract changed")
    require(source.get("new_flags") == list(NEW_FLAGS), "P2 new flags changed before run")
    input_spec = source.get("input_spec_metadata", {})
    require(input_spec.get("sha256") == source.get("input_raw_sha256"), "P2 input pin changed")
    require(prepare_receipt.get("input_raw_sha256") == input_spec.get("sha256"), "P2 input pin does not match source")
    read_bytes(resolve_path(input_spec.get("path"), WORK), input_spec.get("sha256"))
    kernel_path_value = resolve_path(source.get("kernel_path"), REPO)
    runner_path_value = resolve_path(source.get("runner_path"), REPO)
    require(kernel_path_value == REPO / KERNEL_RELATIVE, "P2 kernel path changed")
    require(runner_path_value == REPO / RUNNER_RELATIVE, "P2 runner path changed")
    code_hashes = source.get("current_code_hashes", {})
    require(prepare_receipt.get("kernel_sha256") == code_hashes.get(wire_path(kernel_path_value)), "P2 kernel pin changed")
    require(prepare_receipt.get("runner_sha256") == code_hashes.get(wire_path(runner_path_value)), "P2 runner pin changed")
    read_bytes(kernel_path_value, prepare_receipt.get("kernel_sha256"))
    read_bytes(runner_path_value, prepare_receipt.get("runner_sha256"))
    for name, digest in source.get("source_union_hashes", {}).items():
        read_bytes(resolve_path(name), digest)
    helper = source.get("metadata_helper_artifact", {})
    read_bytes(resolve_path(helper.get("path")), helper.get("sha256"))
    require(helper.get("outside_science_pins") is True, "P2 metadata helper entered science pins")
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "P2 source union count changed")
    return source_raw, source


def journal_append(events, event):
    payload = clone(event) if isinstance(event, dict) else {"event": str(event)}
    payload["schema"] = DIAGNOSTIC_SCHEMA
    payload["runtime_model_metadata"] = None
    events.append(payload)
    with DIAGNOSTIC_CHECKPOINT.open("ab") as stream:
        stream.write(encode(payload))
        stream.flush()
        os.fsync(stream.fileno())


def _int_value(mapping, aliases):
    if not isinstance(mapping, dict):
        return 0
    for key in aliases:
        value = mapping.get(key)
        if isinstance(value, int) and not isinstance(value, bool) and value >= 0:
            return value
    return 0


def _count_value(mapping, aliases):
    """Return (present, value) without treating a missing key as zero evidence."""
    if not isinstance(mapping, dict):
        return False, 0
    for key in aliases:
        if key not in mapping:
            continue
        value = mapping.get(key)
        if isinstance(value, int) and not isinstance(value, bool) and value >= 0:
            return True, value
        return True, 0
    return False, 0


def _durable_shape_counts(value):
    """Count only explicit returned binding rows, never the pass bit."""
    if not isinstance(value, dict):
        return {"source_port_rows": 0, "native_source_port_rows": 0, "C3_source_port_rows": 0,
                "projection_rows": 0, "pressure_route_sides": 0}
    references = value.get("typed_source_references")
    native = references.get("native") if isinstance(references, dict) else None
    c3 = references.get("rebuilt_C3") if isinstance(references, dict) else None

    def side_rows(side):
        if not isinstance(side, dict):
            return 0
        handles = side.get("literal_handles")
        roots = handles.get("root_ids") if isinstance(handles, dict) else None
        if isinstance(roots, list):
            return len(roots)
        roots = side.get("root_ids")
        return len(roots) if isinstance(roots, list) else 0

    native_rows = side_rows(native)
    c3_rows = side_rows(c3)
    views = value.get("whole_function_pressure_numerator_views")
    if not isinstance(views, dict):
        views = value.get("whole_function_views")
    native_view = views.get("native") if isinstance(views, dict) else None
    c3_view = views.get("rebuilt_C3") if isinstance(views, dict) else None

    def projection_rows(side):
        if not isinstance(side, dict):
            return 0
        for key in ("ordinary_Z_projection", "ordinary_Z"):
            rows = side.get(key)
            if isinstance(rows, list):
                return len(rows)
        quotient_projection = side.get("p2_quotient_projection")
        if isinstance(quotient_projection, dict) and isinstance(quotient_projection.get("rows"), list):
            return len(quotient_projection["rows"])
        raw_views = side.get("raw_W_derivative_views")
        if isinstance(raw_views, list):
            return len(raw_views)
        return 0

    projection = value.get("projection")
    projection_native = projection.get("native") if isinstance(projection, dict) else None
    projection_c3 = projection.get("rebuilt_C3") if isinstance(projection, dict) else None
    projected = projection_rows(projection_native) + projection_rows(projection_c3)
    if projected == 0:
        projected = projection_rows(native) + projection_rows(c3)
    if projected == 0:
        projected = projection_rows(native_view) + projection_rows(c3_view)
    route_sides = 0
    for side in (native_view, c3_view):
        if isinstance(side, dict) and isinstance(side.get("pressure_route"), dict):
            route_sides += 1
    return {
        "source_port_rows": native_rows + c3_rows,
        "native_source_port_rows": native_rows,
        "C3_source_port_rows": c3_rows,
        "projection_rows": projected,
        "pressure_route_sides": route_sides,
    }


def _actual_counts(value):
    if not isinstance(value, dict):
        return clone(ZERO_ACTUAL_COUNTS)
    durable = _durable_shape_counts(value)
    aliases = {
        "source_port_rows": ("source_port_rows", "actual_source_port_rows"),
        "native_source_port_rows": ("native_source_port_rows", "actual_native_source_port_rows"),
        "C3_source_port_rows": ("C3_source_port_rows", "actual_C3_source_port_rows"),
        "projection_rows": ("projection_rows", "actual_projection_rows", "ordinary_Z_projection_rows"),
        "pressure_route_sides": (
            "pressure_route_sides",
            "actual_pressure_route_sides",
            "physical_pressure_route_sides",
            "actual_physical_pressure_route_sides",
        ),
        "new_mathematical_proofs": ("new_mathematical_proofs", "actual_new_mathematical_proofs"),
        "new_algebra_family_calls": ("new_algebra_family_calls", "actual_new_algebra_family_calls"),
        "numeric_calls": ("numeric_calls", "actual_numeric_calls"),
        "source_provider_calls": ("source_provider_calls", "actual_source_provider_calls"),
        "old_math_calls": ("old_math_calls", "actual_old_math_calls"),
    }
    candidates = []
    for key in ("actual_counts", "audit_counts", "source_binding_audit_counts", "source_binding_audit"):
        candidate = value.get(key)
        if isinstance(candidate, dict):
            candidates.append(candidate)
            if isinstance(candidate.get("actual_counts"), dict):
                candidates.append(candidate["actual_counts"])
    candidates.append(value)
    result = {}
    for key in ACTUAL_COUNT_FIELDS:
        found = False
        count = 0
        for candidate in candidates:
            # Top-level canonical names are accepted only inside an explicit
            # audit/count map; a top-level planned field is not evidence.
            local_aliases = aliases[key] if candidate is not value else tuple(
                alias for alias in aliases[key] if alias.startswith("actual_")
            )
            present, count = _count_value(candidate, local_aliases)
            if present:
                found = True
                break
        result[key] = count if found else durable[key]
    return result


def _projection_orders(value):
    if not isinstance(value, dict):
        return []
    orders = value.get("actual_projection_orders")
    if orders is None:
        for key in ("audit_counts", "source_binding_audit_counts", "source_binding_audit", "actual_counts"):
            candidate = value.get(key)
            if isinstance(candidate, dict) and isinstance(candidate.get("projection_orders"), list):
                orders = candidate["projection_orders"]
                break
    if orders is None and isinstance(value.get("ordinary_Z_orders"), list):
        orders = value["ordinary_Z_orders"]
    if not isinstance(orders, list):
        return []
    return [item for item in orders if isinstance(item, int) and not isinstance(item, bool)]


def _forbidden_work_evidence(value):
    """Find nonzero work counters anywhere in the saved result/audit tree."""
    aliases = {
        alias: label
        for label, names in FORBIDDEN_WORK_FIELD_ALIASES.items()
        for alias in names
    }
    evidence = {}

    def visit(node, location):
        if isinstance(node, dict):
            for key, child in node.items():
                child_location = location + "." + str(key)
                label = aliases.get(key)
                if label is not None and isinstance(child, (int, float, bool)) and child != 0:
                    evidence[child_location] = {"field": label, "value": clone(child)}
                visit(child, child_location)
        elif isinstance(node, list):
            for index, child in enumerate(node):
                visit(child, location + "[" + str(index) + "]")

    visit(value, "result")
    return evidence


def _add_protocol_violation(result, name, detail):
    violations = result.get("protocol_violations")
    violations = clone(violations) if isinstance(violations, dict) else {}
    violations[name] = clone(detail)
    result["protocol_violations"] = violations


def _typed_source_binding_facts(value):
    """Validate counts from the saved typed source rows, not count aliases."""
    facts = {
        "valid": False,
        "source_port_rows": 0,
        "native_source_port_rows": 0,
        "C3_source_port_rows": 0,
        "projection_rows": 0,
        "projection_orders": [],
        "pressure_route_sides": 0,
        "reasons": [],
    }
    if not isinstance(value, dict):
        facts["reasons"].append("result is not a mapping")
        return facts
    references = value.get("typed_source_references")
    if not isinstance(references, dict):
        facts["reasons"].append("typed_source_references is missing")
        return facts

    def side_facts(name, expected_ids, node_key):
        side = references.get(name)
        if not isinstance(side, dict):
            facts["reasons"].append("typed source side missing: " + name)
            return []
        handles = side.get("literal_handles")
        if not isinstance(handles, dict):
            facts["reasons"].append("typed literal handles missing: " + name)
            return []
        root_ids = handles.get("root_ids")
        if root_ids != list(expected_ids):
            facts["reasons"].append("typed P2 root IDs changed: " + name)
        rows = side.get("ordinary_Z_projection")
        if not isinstance(rows, list) or len(rows) != len(expected_ids):
            facts["reasons"].append("typed ordinary-Z projection rows incomplete: " + name)
            return []
        orders = []
        for index, (row, expected_id) in enumerate(zip(rows, expected_ids)):
            if not isinstance(row, dict):
                facts["reasons"].append("typed ordinary-Z row is not a mapping: " + name)
                continue
            if row.get("order") != index:
                facts["reasons"].append("typed ordinary-Z order changed: " + name)
            if row.get(node_key) != expected_id:
                facts["reasons"].append("typed P2 node ID changed: " + name)
            if row.get("dependency_order") != index + 1:
                facts["reasons"].append("typed dependency order changed: " + name)
            orders.append(row.get("order"))
        return orders

    native_orders = side_facts("native", EXPECTED_NATIVE_P2_ROOT_IDS, "root_node")
    c3_orders = side_facts("rebuilt_C3", EXPECTED_C3_P2_ROOT_IDS, "recipe_node")
    if native_orders == list(EXPECTED_PROJECTION_ORDERS) and c3_orders == list(EXPECTED_PROJECTION_ORDERS):
        facts["projection_orders"] = list(EXPECTED_PROJECTION_ORDERS)
    else:
        facts["reasons"].append("typed projection orders are not exactly 0 through 3 on both sides")
    # Do not use count aliases to establish completeness; these values are
    # derived only from the saved typed side rows inspected above.
    facts["native_source_port_rows"] = len(native_orders)
    facts["C3_source_port_rows"] = len(c3_orders)
    facts["source_port_rows"] = facts["native_source_port_rows"] + facts["C3_source_port_rows"]
    facts["projection_rows"] = len(native_orders) + len(c3_orders)

    views = value.get("whole_function_pressure_numerator_views")
    if not isinstance(views, dict):
        views = value.get("whole_function_views")
    if isinstance(views, dict):
        for name in ("native", "rebuilt_C3"):
            side = views.get(name)
            route = side.get("pressure_route") if isinstance(side, dict) else None
            if not isinstance(route, dict) or not route:
                facts["reasons"].append("typed pressure route missing: " + name)
            else:
                facts["pressure_route_sides"] += 1
    else:
        facts["reasons"].append("typed whole-function views missing")

    facts["valid"] = (
        not facts["reasons"]
        and facts["source_port_rows"] == 8
        and facts["native_source_port_rows"] == 4
        and facts["C3_source_port_rows"] == 4
        and facts["projection_rows"] == 8
        and facts["projection_orders"] == list(EXPECTED_PROJECTION_ORDERS)
        and facts["pressure_route_sides"] == 2
    )
    return facts


def fallback(source, calls, imported, error=None):
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "runtime_model_metadata": None,
        "kernel_imported": bool(imported),
        "derive_entry_calls": calls,
        "derive_calls": calls,
        "entry_attempts": 1,
        "input_digest": source.get("input_digest"),
        "input_raw_sha256": source.get("input_raw_sha256"),
        "blockers": [str(error)] if error else ["P2 derive returned no result"],
        "actual_counts": clone(ZERO_ACTUAL_COUNTS),
        "actual_projection_orders": [],
        "new_flags": {key: False for key in NEW_FLAGS},
        "new_mathematical_proofs": 0,
        "new_algebra_family_calls": 0,
        "numeric_calls": 0,
        "source_provider_calls": 0,
        "old_math_calls": 0,
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "full_goal_complete": False,
    }


def decorate_result_metadata(value, source, calls, imported):
    result = clone(value) if isinstance(value, dict) else fallback(source, calls, imported, "derive returned a non-mapping result")
    kernel_result_schema = result.get("schema")
    if kernel_result_schema is not None:
        result["kernel_result_schema"] = kernel_result_schema
    kernel_flags = result.get("new_flags")
    if not isinstance(kernel_flags, dict):
        kernel_flags = result.get("flags")
    if not isinstance(kernel_flags, dict):
        kernel_flags = {}
    kernel_passed = result.get("passed")
    if not isinstance(kernel_passed, bool):
        kernel_passed = False
        _add_protocol_violation(result, "missing_explicit_passed", {"required": True})
    result["passed"] = kernel_passed
    forbidden_evidence = _forbidden_work_evidence(result)
    if forbidden_evidence:
        _add_protocol_violation(
            result,
            "forbidden_nonzero_work_counts",
            {
                "evidence": forbidden_evidence,
                "binding_kernel_must_not_claim_math_or_numeric_work": True,
            },
        )
        # Preserve every source/audit count verbatim.  The durable result is
        # intentionally rejected after its fsync; no count is normalized to 0.
        result["passed"] = False
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result["schema"] = RESULT_SCHEMA
    result["component"] = COMPONENT
    result["passed"] = result.get("passed") is True and not result.get("protocol_violations")
    result["source_only"] = True
    result["metadata_only"] = False
    result["runtime_model_metadata"] = None
    result["input_digest"] = source.get("input_digest")
    result["input_raw_sha256"] = source.get("input_raw_sha256")
    result["kernel_imported"] = bool(imported)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    actual = _actual_counts(result)
    typed_facts = _typed_source_binding_facts(result)
    # Valid results receive the manager's canonical count view.  When the
    # kernel/audit already contains forbidden work, retain its raw maps and
    # counters verbatim so the post-fsync validator can reject them honestly.
    if "actual_counts" not in result or not forbidden_evidence:
        result["actual_counts"] = actual
    if not forbidden_evidence:
        for key in ("source_port_rows", "native_source_port_rows", "C3_source_port_rows", "projection_rows", "pressure_route_sides"):
            result[key] = actual[key]
    result["actual_projection_orders"] = (
        list(typed_facts["projection_orders"])
        if typed_facts["valid"]
        else _projection_orders(result)
    )
    result["new_flags"] = {
        key: result.get(key) is True or kernel_flags.get(key) is True
        for key in NEW_FLAGS
    }
    result["flags"] = clone(result["new_flags"])
    for key in NEW_FLAGS:
        result[key] = result["new_flags"][key]
    if result["passed"] and calls != 1:
        _add_protocol_violation(
            result,
            "successful_derive_call_count",
            {"derive_entry_calls": calls, "derive_calls": calls, "required": 1},
        )
        result["passed"] = False
    # A kernel may return a nominal pass with incomplete typed source facts.
    # Convert that protocol mismatch to an honest failed saved result while
    # retaining the rows and orders it actually returned.
    if result["passed"]:
        if not typed_facts["valid"]:
            _add_protocol_violation(
                result,
                "incomplete_typed_source_facts",
                {"reasons": list(typed_facts["reasons"])},
            )
        complete_binding = (
            result["new_flags"] == {key: True for key in NEW_FLAGS}
            and typed_facts["valid"]
            and typed_facts["source_port_rows"] == PLANNING_CONTRACT["source_port_rows"]
            and typed_facts["native_source_port_rows"] == PLANNING_CONTRACT["native_source_port_rows"]
            and typed_facts["C3_source_port_rows"] == PLANNING_CONTRACT["C3_source_port_rows"]
            and typed_facts["projection_rows"] == PLANNING_CONTRACT["projection_rows"]
            and typed_facts["pressure_route_sides"] == PLANNING_CONTRACT["pressure_route_sides"]
            and typed_facts["projection_orders"] == PLANNING_CONTRACT["projection_orders"]
            and actual["source_port_rows"] == typed_facts["source_port_rows"]
            and actual["native_source_port_rows"] == typed_facts["native_source_port_rows"]
            and actual["C3_source_port_rows"] == typed_facts["C3_source_port_rows"]
            and actual["projection_rows"] == typed_facts["projection_rows"]
            and actual["pressure_route_sides"] == typed_facts["pressure_route_sides"]
        )
        if not complete_binding:
            _add_protocol_violation(
                result,
                "incomplete_typed_source_binding_contract",
                {"facts": typed_facts},
            )
            result["passed"] = False
    if not result["passed"]:
        result["new_flags"] = {key: False for key in NEW_FLAGS}
        result["flags"] = clone(result["new_flags"])
        for key in NEW_FLAGS:
            result[key] = False
    if "audit_counts" not in result:
        result["audit_counts"] = clone(actual)
    result["source_packet_defined"] = result["passed"]
    result["source_packet_definition_accepted"] = result["passed"]
    result["source_packet_accepted_owned"] = result["passed"]
    result["source_packet_accepted_inherited"] = False
    result["parent_result_raw_sha256"] = source.get("parent_result_raw_sha256")
    result["planning_contract"] = clone(source["planning_contract"])
    result["planned_counts_are_not_results"] = True
    for key in ("full_goal_complete", "numeric_values_computed", "numeric_field_values_present"):
        if result.get(key) is True:
            _add_protocol_violation(result, "forbidden_completion_or_numeric_flag", {key: True})
        else:
            result.setdefault(key, False)
    result["runtime_model_metadata"] = None
    for key, flag in source.get("inherited_parent_true_flags", {}).items():
        if flag is True and key not in BOOKKEEPING_FIELDS and key not in (*NEW_FLAGS, *NONCLAIM_FLAGS):
            result.setdefault(key, True)
    for key in source.get("inherited_false_flags", []):
        result[key] = False
    for key in OLD_QUIET_FLAGS:
        result.setdefault(key, False)
    for key in NONCLAIM_FLAGS:
        if result.get(key) is True:
            _add_protocol_violation(result, "forbidden_nonclaim_flag", {key: True})
        else:
            result.setdefault(key, False)
    if result.get("protocol_violations"):
        result["passed"] = False
        result["new_flags"] = {key: False for key in NEW_FLAGS}
        result["flags"] = clone(result["new_flags"])
        for key in NEW_FLAGS:
            result[key] = False
        result["source_packet_defined"] = False
        result["source_packet_definition_accepted"] = False
        result["source_packet_accepted_owned"] = False
    return result


def validate_actual_counts(actual, passed):
    require(isinstance(actual, dict), "P2 actual counts are not a mapping")
    require(set(actual) == set(ACTUAL_COUNT_FIELDS), "P2 actual count fields changed")
    for key, value in actual.items():
        require(isinstance(value, int) and not isinstance(value, bool) and value >= 0, "P2 actual count is malformed: " + key)
    require(actual["source_port_rows"] <= PLANNING_CONTRACT["source_port_rows"], "P2 source rows exceed the plan")
    require(actual["native_source_port_rows"] <= PLANNING_CONTRACT["native_source_port_rows"], "P2 native rows exceed the plan")
    require(actual["C3_source_port_rows"] <= PLANNING_CONTRACT["C3_source_port_rows"], "P2 C3 rows exceed the plan")
    require(actual["projection_rows"] <= PLANNING_CONTRACT["projection_rows"], "P2 projection rows exceed the plan")
    require(actual["pressure_route_sides"] <= PLANNING_CONTRACT["pressure_route_sides"], "P2 route sides exceed the plan")
    for key in ("new_mathematical_proofs", "new_algebra_family_calls", "numeric_calls", "source_provider_calls", "old_math_calls"):
        require(actual[key] == 0, "P2 binding recorded forbidden work: " + key)
    require(
        actual["source_port_rows"] == actual["native_source_port_rows"] + actual["C3_source_port_rows"],
        "P2 source row total does not equal native plus C3 rows",
    )
    if passed:
        require(actual["source_port_rows"] == PLANNING_CONTRACT["source_port_rows"], "P2 source rows are not fully bound")
        require(actual["native_source_port_rows"] == PLANNING_CONTRACT["native_source_port_rows"], "P2 native source rows are incomplete")
        require(actual["C3_source_port_rows"] == PLANNING_CONTRACT["C3_source_port_rows"], "P2 C3 source rows are incomplete")
        require(actual["projection_rows"] == PLANNING_CONTRACT["projection_rows"], "P2 projection rows are incomplete")
        require(actual["pressure_route_sides"] == PLANNING_CONTRACT["pressure_route_sides"], "P2 route sides are incomplete")


def validate_result_metadata(result, source):
    require(result.get("schema") == RESULT_SCHEMA, "P2 result schema changed")
    require(result.get("component") == COMPONENT, "P2 result component changed")
    require(isinstance(result.get("passed"), bool), "P2 result passed field is not explicit")
    require(result.get("metadata_only") is False and result.get("source_only") is True, "P2 result scope changed")
    require(result.get("runtime_model_metadata") is None, "P2 runtime metadata is not null")
    require(result.get("derive_entry_calls") in (0, 1), "P2 derive entry was retried")
    require(result.get("derive_calls") == result.get("derive_entry_calls") and result.get("entry_attempts") == 1, "P2 derive bookkeeping changed")
    require(not _forbidden_work_evidence(result), "P2 saved result contains nonzero forbidden work counts")
    require(not result.get("protocol_violations"), "P2 result contains a protocol violation")
    if result.get("passed") is True:
        require(result.get("derive_entry_calls") == 1 and result.get("derive_calls") == 1, "Successful P2 derive did not run exactly once")
    require(result.get("input_digest") == source.get("input_digest"), "P2 input digest changed")
    require(result.get("input_raw_sha256") == source.get("input_raw_sha256"), "P2 input raw SHA changed")
    require(result.get("source_packet_accepted_inherited") is False, "P2 source acceptance was inherited")
    require(result.get("full_goal_complete") is False, "P2 overclaims full goal completion")
    require(result.get("numeric_values_computed") is not True and result.get("numeric_field_values_present") is not True, "P2 contains numeric values")
    require(result.get("new_mathematical_proofs") == 0, "P2 new mathematical proofs are not zero")
    require(result.get("new_algebra_family_calls") == 0, "P2 algebra-family calls are not zero")
    require(result.get("numeric_calls") == 0 and result.get("source_provider_calls") == 0 and result.get("old_math_calls") == 0, "P2 forbidden calls are nonzero")
    require(result.get("flags") == result.get("new_flags"), "P2 source-binding flags changed")
    actual = result.get("actual_counts")
    validate_actual_counts(actual, result.get("passed") is True)
    for key in ("source_port_rows", "native_source_port_rows", "C3_source_port_rows", "projection_rows", "pressure_route_sides"):
        require(result.get(key) == actual[key], "P2 actual count alias changed: " + key)
    require(result.get("audit_counts") == actual, "P2 audit counts changed")
    orders = result.get("actual_projection_orders")
    require(isinstance(orders, list) and len(set(orders)) == len(orders), "P2 projection orders are malformed")
    require(orders == sorted(orders), "P2 projection orders are not ordered")
    require(all(order in PLANNING_CONTRACT["projection_orders"] for order in orders), "P2 projection order is outside the plan")
    if result.get("passed") is False:
        require(not any(result.get(key) is True for key in NEW_FLAGS), "Failed P2 result overclaims a new flag")
        require(not result.get("source_packet_defined"), "Failed P2 result defines an accepted source packet")
        return {
            "passed": False,
            "honest_failure": True,
            "actual_counts": clone(actual),
            "actual_projection_orders": list(orders),
            "forbidden_work_calls": 0,
        }
    require(all(result.get(key) is True for key in NEW_FLAGS), "P2 flags are not all source-bound")
    require(result.get("source_packet_definition_accepted") is True and result.get("source_packet_accepted_owned") is True, "P2 acceptance is not owned")
    typed_facts = _typed_source_binding_facts(result)
    require(typed_facts["valid"], "P2 typed source facts are incomplete")
    require(orders == typed_facts["projection_orders"], "P2 projection orders are not from typed source facts")
    require(actual["source_port_rows"] == typed_facts["source_port_rows"], "P2 source rows are not from typed source facts")
    require(actual["native_source_port_rows"] == typed_facts["native_source_port_rows"], "P2 native rows are not from typed source facts")
    require(actual["C3_source_port_rows"] == typed_facts["C3_source_port_rows"], "P2 C3 rows are not from typed source facts")
    require(actual["projection_rows"] == typed_facts["projection_rows"], "P2 projection rows are not from typed source facts")
    require(actual["pressure_route_sides"] == typed_facts["pressure_route_sides"], "P2 pressure routes are not from typed source facts")
    require(orders == PLANNING_CONTRACT["projection_orders"], "P2 projection order plan is incomplete")
    return {
        "passed": True,
        "actual_counts": clone(actual),
        "actual_projection_orders": list(orders),
        "forbidden_work_calls": 0,
    }


def run(args):
    require(not ENTERED_CHECKPOINT.exists() and not RESULT_CHECKPOINT.exists() and not DIAGNOSTIC_CHECKPOINT.exists(), "P2 attempt already entered")
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
        read_bytes(kernel, source["current_code_hashes"][wire_path(kernel)])
        spec = importlib.util.spec_from_file_location("_compact_O3_transition_P2_source_binding_entry", kernel)
        require(spec is not None and spec.loader is not None, "Cannot create P2 kernel import specification")
        module = importlib.util.module_from_spec(spec)
        imported = True
        spec.loader.exec_module(module)
        derive = getattr(module, "derive", None)
        require(callable(derive), "P2 kernel lacks derive(source, journal_callback)")

        def callback(event):
            journal_append(events, event)

        calls = 1
        result = derive(clone(source), callback)
    except Exception as error:
        journal_append(events, {"stage": "derive_exception", "error_type": type(error).__name__, "error": str(error)})
        result = fallback(source, calls, imported, error)
    result = decorate_result_metadata(result, source, calls, imported)
    result["source_checkpoint_sha256"] = sha_bytes(source_raw)
    result["entry_checkpoint_sha256"] = sha_bytes(entered_raw)
    result["diagnostic_checkpoint_sha256"] = sha_path(DIAGNOSTIC_CHECKPOINT)
    result["diagnostic_event_count"] = len(events)
    # Exclusive durable result: validation and repair metadata never rewrite it.
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
        "actual_counts": clone(result["actual_counts"]),
        "actual_projection_orders": list(result["actual_projection_orders"]),
        "derive_entry_calls": calls,
        "kernel_imported": imported,
        "protocol_violations": clone(result.get("protocol_violations", {})),
        "protocol_violation_receipt": bool(result.get("protocol_violations")),
        "runtime_model_metadata": None,
    }
    receipt_raw = write_replace(RUN_RECEIPT, receipt)
    return {"run": receipt, "run_sha256": sha_bytes(receipt_raw), "validation": validation, "validation_error": validation_error}


def validate_result_checkpoint_links(source_raw, result_raw, result):
    entered_raw = read_bytes(ENTERED_CHECKPOINT)
    diagnostic_raw = read_bytes(DIAGNOSTIC_CHECKPOINT)
    require(result.get("source_checkpoint_sha256") == sha_bytes(source_raw), "P2 saved result source SHA changed")
    require(result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw), "P2 saved result entry SHA changed")
    require(result.get("diagnostic_checkpoint_sha256") == sha_bytes(diagnostic_raw), "P2 saved result diagnostic SHA changed")
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


def failed_union(source, result, current_artifacts):
    if result.get("passed") is True:
        return clone(source["preserved_failed_artifact_hashes_union"]), PARENT_FAILED_COUNT
    additions = {entry["path"]: entry["sha256"] for entry in current_artifacts}
    additions[wire_path(PREPARE_RECEIPT)] = sha_path(PREPARE_RECEIPT)
    additions[wire_path(RUN_RECEIPT)] = sha_path(RUN_RECEIPT)
    additions[wire_path(Path(__file__).resolve())] = sha_path(Path(__file__).resolve())
    require(len(additions) == 7, "P2 failed retention set is not exactly seven artifacts")
    merged = merge_hashes(source["preserved_failed_artifact_hashes_union"], additions)
    require(len(merged) == FAILED_COUNT_IF_NEW_FAILURE, "P2 failed retention did not add exactly seven artifacts")
    return merged, len(merged)


def report_for_publish(source, result, result_raw, code2, input_meta):
    current_ck = checkpoint_artifacts(embed=True)
    helper = artifact_entry(Path(__file__).resolve(), source["metadata_helper_artifact"]["sha256"], embed=True)
    helper["outside_science_pins"] = True
    input_artifact = artifact_entry(Path(input_meta["path"]), input_meta["sha256"], embed=True)
    packet_receipts = [
        artifact_entry(resolve_path(item["path"]), item["sha256"], embed=True)
        for item in source.get("raw_packet_receipts", [])
    ]
    current_failed_union, failed_count = failed_union(source, result, current_ck)
    artifacts = [
        artifact_entry(Path(source["parent_report_metadata"]["path"]), source["parent_report_metadata"]["sha256"]),
        artifact_entry(Path(source["parent_result_metadata"]["path"]), source["parent_result_metadata"]["sha256"]),
        *[clone(item) for item in code2],
        input_artifact,
        helper,
        *packet_receipts,
        *current_ck,
    ]
    actual = clone(result["actual_counts"])
    passed = result.get("passed") is True
    report = {
        "schema": PUBLICATION_SCHEMA,
        "component": COMPONENT,
        "attempt_status": "passed" if passed else "failed",
        "passed": passed,
        "source_only": True,
        "publication_metadata_only": True,
        "runtime_model_metadata": None,
        "parent_report": source["parent_report_metadata"]["path"],
        "parent_report_sha256": source["parent_report_sha256"],
        "parent_commit": source["parent_commit"],
        "parent_head": source["parent_head"],
        "parent_report_commit": source["parent_report_commit"],
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
        "preserved_error_cap_metadata": clone(source["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(source["retained_construction_data"]),
        "retained_integral_term_inventory": clone(source["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(source["retained_cap_provenance"]),
        "inherited_code2": clone(source.get("inherited_code2", [])),
        "planning_contract": clone(source["planning_contract"]),
        "planned_counts_are_not_results": True,
        "actual_counts": actual,
        "failure_actual_counts": clone(actual) if not passed else None,
        "source_port_rows": actual["source_port_rows"],
        "native_source_port_rows": actual["native_source_port_rows"],
        "C3_source_port_rows": actual["C3_source_port_rows"],
        "projection_rows": actual["projection_rows"],
        "pressure_route_sides": actual["pressure_route_sides"],
        "actual_projection_orders": list(result["actual_projection_orders"]),
        "new_flags": {key: result.get(key) is True for key in NEW_FLAGS},
        "new_mathematical_proofs": result.get("new_mathematical_proofs", actual["new_mathematical_proofs"]),
        "new_algebra_family_calls": result.get("new_algebra_family_calls", actual["new_algebra_family_calls"]),
        "numeric_calls": result.get("numeric_calls", actual["numeric_calls"]),
        "source_provider_calls": result.get("source_provider_calls", actual["source_provider_calls"]),
        "old_math_calls": result.get("old_math_calls", actual["old_math_calls"]),
        "protocol_violations": clone(result.get("protocol_violations", {})),
        "numeric_values_computed": result.get("numeric_values_computed", False),
        "numeric_field_values_present": result.get("numeric_field_values_present", False),
        "source_packet_defined": passed,
        "source_packet_definition_accepted": passed,
        "source_packet_accepted_owned": passed,
        "source_packet_accepted_inherited": False,
        "original_source_evaluations": 0,
        "full_result": clone(result),
        "component_result": clone(result),
        "component_result_exact": True,
        "full_result_raw_sha256": sha_bytes(result_raw),
        "full_result_raw_bytes": len(result_raw),
        "full_result_artifact_selector": "$.lossless_artifacts[path=" + wire_path(RESULT_CHECKPOINT) + "]",
        "lossless_artifacts": artifacts,
        "lossless_artifact_count": len(artifacts),
        "raw_packet_receipts": packet_receipts,
        "raw_bytes_lossless": True,
        "raw_sha_pins": clone(source["raw_sha_pins"]),
        "input_spec_metadata": clone(source["input_spec_metadata"]),
        "metadata_helper_artifact": clone(helper),
        "input_artifact": clone(input_artifact),
        "code2": [clone(item) for item in code2],
        "goal_active": True,
        "full_goal_complete": False,
        "next_executable_task": next_task_metadata(passed),
        "external_draft": clone(source["metadata_helper_artifact"]),
        "math_reexecuted": False,
        "ancestor_kernel_imported": False,
        "scope_attribute_lines": list(SCOPE_ATTRIBUTE_LINES),
        "exact_nine_scope": True,
    }
    for key, value in source["inherited_parent_true_flags"].items():
        if key not in BOOKKEEPING_FIELDS and key not in NEW_FLAGS and key not in NONCLAIM_FLAGS:
            report.setdefault(key, True)
    for key in source["inherited_false_flags"]:
        report[key] = False
    for key in OLD_QUIET_FLAGS:
        report[key] = False
    for key in NONCLAIM_FLAGS:
        report[key] = False
    for key in NEW_FLAGS:
        report[key] = result.get(key) is True
    return report


def validate_embedded_artifact(entry):
    raw = read_bytes(resolve_path(entry.get("path")), entry.get("sha256"))
    require(entry.get("raw_readonly") is True, "Lossless artifact is not read-only")
    require(entry.get("bytes") == len(raw), "Lossless artifact byte count changed")
    if "raw_bytes_base64" in entry:
        require(entry.get("raw_bytes_exact") is True, "Embedded artifact is not marked exact")
        require(base64.b64decode(entry["raw_bytes_base64"].encode("ascii")) == raw, "Embedded artifact bytes changed")


def validate_publication_metadata(report, source):
    require(report.get("source_pins") == source.get("source_union_hashes"), "Published source pins changed")
    require(report.get("source_union_count") == len(source.get("source_union_hashes", {})), "Published source union count changed")
    require(report.get("parent_report_sha256") == PARENT_REPORT_SHA256, "Published parent report SHA changed")
    require(report.get("parent_result_raw_sha256") == PARENT_RESULT_SHA256, "Published parent result SHA changed")
    require(report.get("parent_commit") == PARENT_HEAD and report.get("parent_head") == PARENT_HEAD, "Published parent HEAD changed")
    require(report.get("parent_report_commit") == PARENT_REPORT_COMMIT, "Published parent report commit changed")
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Published checkpoint count changed")
    require(report.get("planning_contract") == PLANNING_CONTRACT, "Published planning contract changed")
    require(not report.get("protocol_violations"), "Published report contains a protocol violation")
    actual = report.get("actual_counts")
    validate_actual_counts(actual, report.get("passed") is True)
    require(report.get("failure_actual_counts") == (None if report.get("passed") is True else actual), "Published failure counts changed")
    for key in ("source_port_rows", "native_source_port_rows", "C3_source_port_rows", "projection_rows", "pressure_route_sides"):
        require(report.get(key) == actual[key], "Published actual count alias changed: " + key)
    require(report.get("actual_projection_orders") == sorted(set(report.get("actual_projection_orders", []))), "Published projection orders changed")
    require(report.get("new_mathematical_proofs") == 0 and report.get("new_algebra_family_calls") == 0, "Published math counts changed")
    require(report.get("numeric_calls") == 0 and report.get("source_provider_calls") == 0 and report.get("old_math_calls") == 0, "Published forbidden counts changed")
    require(report.get("new_flags") == {key: report.get(key) is True for key in NEW_FLAGS}, "Published flags changed")
    if report.get("passed") is True:
        typed_facts = _typed_source_binding_facts(report.get("full_result"))
        require(typed_facts["valid"], "Published typed source facts changed")
        require(actual["source_port_rows"] == typed_facts["source_port_rows"], "Published source rows changed from typed facts")
        require(actual["projection_rows"] == typed_facts["projection_rows"], "Published projection rows changed from typed facts")
        require(actual["pressure_route_sides"] == typed_facts["pressure_route_sides"], "Published pressure routes changed from typed facts")
    require(report.get("source_packet_defined") is (report.get("passed") is True), "Published packet definition changed")
    require(report.get("source_packet_accepted_inherited") is False, "Published inherited acceptance changed")
    require(report.get("raw_sha_pins") == source.get("raw_sha_pins"), "Published raw SHA pins changed")
    require(report.get("scoped_paths") == [path.as_posix() for path in SCOPED_RELATIVES], "Published nine-path scope changed")
    require(report.get("exact_nine_scope") is True, "Published nine-path marker changed")
    require(report.get("scope_attribute_lines") == list(SCOPE_ATTRIBUTE_LINES), "Published scope attributes changed")
    current = checkpoint_artifacts(embed=False)
    require(report.get("checkpoint_hashes") == {entry["path"]: entry["sha256"] for entry in current}, "Published checkpoint hashes changed")
    require(report.get("checkpoint_paths") == [entry["path"] for entry in current], "Published checkpoint paths changed")
    for entry in report.get("lossless_artifacts", []):
        validate_embedded_artifact(entry)
    helper = report.get("metadata_helper_artifact")
    require(isinstance(helper, dict) and helper.get("outside_science_pins") is True, "Published helper scope changed")
    require(helper.get("path") not in report.get("source_pins", {}), "Metadata helper entered source pins")
    return {
        "source_pin_count": report.get("source_pin_count"),
        "checkpoint_count": report.get("checkpoint_count"),
        "actual_counts": clone(actual),
        "actual_projection_orders": list(report.get("actual_projection_orders", [])),
        "raw_sha_pins_verified": True,
        "metadata_only": True,
    }


def handoff_preimage(path, raw):
    preimage = WORK / f"{COMPONENT}_handoff_preimage_{path.stem}_{sha_bytes(raw)[:12]}.md"
    if not preimage.exists():
        write_new(preimage, raw)
    return {"path": wire_path(preimage), "sha256": sha_path(preimage), "bytes": len(raw), "raw_readonly": True}


def hydrate_publish_args(args, source):
    if not args.parent_report:
        args.parent_report = source.get("parent_report_metadata", {}).get("path")
    if not args.parent_report_sha256:
        args.parent_report_sha256 = source.get("parent_report_sha256")
    if not args.parent_commit:
        args.parent_commit = source.get("parent_commit")


def publish(args):
    source_raw, source = load_source()
    hydrate_publish_args(args, source)
    require(RESULT_CHECKPOINT.exists(), "P2 result checkpoint is missing")
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
    report["scoped_paths"] = [path.as_posix() for path in SCOPED_RELATIVES]
    report["publication_metadata_validation"] = validate_publication_metadata(report, source)
    report_raw = write_replace(REPO / REPORT_RELATIVE, report)
    actual = result["actual_counts"]
    milestone = (
        "# Current original O3 transition P2 source binding (2026-10-27)\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}`; SHA-256 `{sha_bytes(report_raw)}`.\n\n"
        f"Attempt status: `{report['attempt_status']}`; source-only: `true`; actual source-port rows: `{actual['source_port_rows']}` (planned `{PLANNING_CONTRACT['source_port_rows']}`); actual native/C3 rows: `{actual['native_source_port_rows']}/{actual['C3_source_port_rows']}` (planned `4/4`); actual projection rows: `{actual['projection_rows']}` (planned `{PLANNING_CONTRACT['projection_rows']}`); actual projection orders: `{report['actual_projection_orders']}` (planned `{PLANNING_CONTRACT['projection_orders']}`); actual pressure-route sides: `{actual['pressure_route_sides']}` (planned `{PLANNING_CONTRACT['pressure_route_sides']}`).\n\n"
        "This leaf binds source metadata only. New mathematical proofs, new algebra-family calls, numeric calls, source-provider calls, and old-math calls are zero. A failed attempt reports only durable saved rows; planned rows are not results. No P2 zero claim, phase equality, cap closure, global closure, physical completion, or full-goal completion is asserted.\n\n"
        "The entered attempt is immutable and must never be reentered.\n\n"
        f"Next task: `{next_task_metadata(report['passed'])['phase']}`.\n"
    ).encode("utf-8")
    milestone_raw = write_replace(REPO / MILESTONE_RELATIVE, milestone)
    handoff_prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}` (SHA-256 `{sha_bytes(report_raw)}`).\n"
        f"Status: `{report['attempt_status']}`; actual source-port rows `{actual['source_port_rows']}` / planned `{PLANNING_CONTRACT['source_port_rows']}`; projection orders `{report['actual_projection_orders']}` / planned `{PLANNING_CONTRACT['projection_orders']}`.\n"
        "Do not reenter this attempt. After a pass, run a fresh independent P2 pressure-gap operator; after a failure, inspect only the cached binding prefix.\n"
    ).encode("utf-8")
    preimages = []
    for relative in HANDOFF_RELATIVES:
        path = REPO / relative
        old = path.read_bytes() if path.exists() else b""
        preimages.append(handoff_preimage(path, old))
        write_replace(path, handoff_prefix + old)
    attrs_path = REPO / ".gitattributes"
    attrs_old = attrs_path.read_bytes() if attrs_path.exists() else b""
    attrs_text = attrs_old.decode("utf-8", errors="strict")
    for line in SCOPE_ATTRIBUTE_LINES:
        if line not in attrs_text.splitlines():
            attrs_text += ("\n" if attrs_text and not attrs_text.endswith("\n") else "") + line + "\n"
    write_replace(attrs_path, attrs_text.encode("utf-8"))
    validate_scope_attribute_lines()
    main_brief = update_main_brief(
        {
            "phase": "published_passed" if report["passed"] else "published_failed",
            "report": wire_path(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "milestone": wire_path(REPO / MILESTONE_RELATIVE),
            "milestone_sha256": sha_bytes(milestone_raw),
            "attempt_status": report["attempt_status"],
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "next_executable_task": next_task_metadata(report["passed"]),
            "actual_counts": clone(actual),
            "actual_projection_orders": list(report["actual_projection_orders"]),
            "accepted_progress": report["passed"],
        }
    )
    receipt = {
        "schema": COMPONENT + "_publication_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "attempt_status": report["attempt_status"],
        "scoped_paths": [path.as_posix() for path in SCOPED_RELATIVES],
        "scope_attribute_lines": list(SCOPE_ATTRIBUTE_LINES),
        "exact_nine_scope": True,
        "handoff_preimages": preimages,
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
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "P2 publication checkpoint count changed")
    require(report.get("source_pin_count", 0) > PARENT_SOURCE_PIN_COUNT, "P2 source pins did not advance")
    require(report.get("full_goal_complete") is False and report.get("goal_active") is True, "P2 goal state changed")
    head = git_output("rev-parse", "HEAD")
    require(valid_git_revision(head), "Local HEAD is not an exact Git SHA-1 revision")
    require(report.get("parent_commit") == PARENT_HEAD and report.get("parent_head") == PARENT_HEAD, "P2 published parent HEAD changed")
    require(git_output("rev-parse", "HEAD^") == PARENT_HEAD, "HEAD parent does not match authenticated P2 parent")
    branch = git_output("branch", "--show-current")
    require(branch, "Current Git branch is unnamed")
    remote_line = git_output("ls-remote", "origin", "refs/heads/" + branch)
    remote_parts = remote_line.split()
    require(len(remote_parts) >= 2 and valid_git_revision(remote_parts[0]), "Remote branch lookup returned no exact Git SHA-1")
    remote = remote_parts[0].lower()
    require(head.lower() == remote, "HEAD is not synchronized with origin exact branch")
    changed = set(git_output("diff", "--name-only", PARENT_HEAD + "..HEAD").splitlines())
    expected = {path.as_posix() for path in SCOPED_RELATIVES}
    require(changed == expected, "Committed scope is not exactly the nine P2 paths")
    for relative in SCOPED_RELATIVES:
        local_path = REPO / relative
        require(local_path.exists(), "Published scoped path is missing: " + relative.as_posix())
        blob_size = int(git_output("cat-file", "-s", "HEAD:" + relative.as_posix()))
        require(blob_size == local_path.stat().st_size, "Published cat-file byte count changed: " + relative.as_posix())
        committed_blob = git_output("rev-parse", "HEAD:" + relative.as_posix()).lower()
        local_blob = git_output("hash-object", wire_path(local_path)).lower()
        require(valid_git_revision(committed_blob) and committed_blob == local_blob, "Published bytes differ from committed scope: " + relative.as_posix())
    validate_scope_attribute_lines()
    return head, remote


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
        "new_flags": clone(report.get("new_flags", {})),
        "actual_counts": clone(report.get("actual_counts", {})),
        "actual_projection_orders": list(report.get("actual_projection_orders", [])),
        "source_pin_count": report.get("source_pin_count"),
        "checkpoint_count": report.get("checkpoint_count"),
        "preserved_failed_artifact_count": report.get("preserved_failed_artifact_count"),
        "runtime_model_metadata": None,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "latest_review": clone(review) if review is not None else None,
        "latest_qa": clone(qa) if qa is not None else None,
    }


def resume(args):
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "P2 publication report is missing")
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
            "next_executable_task": clone(report["next_executable_task"]),
            "latest_authoritative_state": authority,
            "verified_local_and_remote_commit": head,
            "next_remaining": clone(report["next_executable_task"]),
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
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "P2 report is not published")
    source_raw, source = load_source()
    require(valid_sha(report.get("full_result_raw_sha256")), "Published result SHA is malformed")
    result_raw, result = read_json(RESULT_CHECKPOINT, report["full_result_raw_sha256"])
    checkpoint_links = validate_result_checkpoint_links(source_raw, result_raw, result)
    require(report.get("result_checkpoint_links") == checkpoint_links, "Published result checkpoint links changed before seal")
    require(report.get("full_result") == result and report.get("component_result") == result, "Published full result changed before seal")
    require(report.get("full_result_raw_bytes") == len(result_raw), "Published result byte count changed before seal")
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Published checkpoint count changed before seal")
    require(report.get("preserved_failed_artifact_count") in (PARENT_FAILED_COUNT, FAILED_COUNT_IF_NEW_FAILURE), "Published failed count changed before seal")
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
            "next_executable_task": clone(report["next_executable_task"]),
            "latest_authoritative_state": authority,
            "verified_local_and_remote_commit": head,
            "next_remaining": clone(report["next_executable_task"]),
            "milestone_sha256": sha_bytes(milestone_raw),
            "latest_review": review_meta,
            "latest_qa": qa,
            "accepted_progress": report.get("passed") is True,
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
    parser_value = argparse.ArgumentParser(description=__doc__)
    parser_value.add_argument("action", choices=("prepare", "run", "publish", "resume", "seal"))
    parser_value.add_argument("--input-spec")
    parser_value.add_argument("--input-sha256")
    parser_value.add_argument("--parent-report")
    parser_value.add_argument("--parent-report-sha256")
    parser_value.add_argument("--parent-commit")
    parser_value.add_argument("--kernel-path")
    parser_value.add_argument("--kernel-sha256")
    parser_value.add_argument("--runner-path")
    parser_value.add_argument("--runner-sha256")
    parser_value.add_argument("--review")
    parser_value.add_argument("--review-sha256")
    return parser_value


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
