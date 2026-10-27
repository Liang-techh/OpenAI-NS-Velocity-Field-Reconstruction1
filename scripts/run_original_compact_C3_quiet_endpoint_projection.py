"""Metadata-only manager for the compact C3 quiet endpoint projection leaf.

This file is deliberately a durable orchestration scaffold.  ``prepare``
authenticates the supplied parent publication, the saved failed quiet-power
attempt, the old accepted parent result, the current input pin, and the two
code artifacts.  It writes one source packet.  ``run`` enters once, fsyncs the
entered marker before importing only the new projection kernel, and invokes
its ``derive`` entry point once.  ``publish``, ``resume``, and ``seal`` only
consume saved metadata and raw hashes.  No ancestor source, provider, proof,
checker, normalizer, quadrature, recurrence, density, or pair function is
imported or replayed by this manager.

The fresh kernel is intentionally supplied by the C3 source worker.  Its
source path and hash are frozen at ``prepare``; absence of that source is a
hard preparation error rather than an invitation to use an ancestor module.
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
import sys
import traceback


REPO = Path(
    r"C:/Users/z5242/Documents/Codex/2026-09-13/"
    r"before-doing-any-project-work-verify/work/st073-transition-next"
)
WORK = Path(
    r"C:/Users/z5242/Documents/Codex/2026-10-16/"
    r"https-github-com-liang-techh-openai/work"
)

COMPONENT = "compact_C3_quiet_endpoint_projection"
INPUT_SCHEMA = "original_compact_C3_quiet_endpoint_projection_inputs_v1"
RESULT_SCHEMA = "original_compact_C3_quiet_endpoint_projection"
ENTERED_SCHEMA = "original_compact_C3_quiet_endpoint_projection_entered_v1"
DIAGNOSTIC_SCHEMA = "compact_C3_quiet_endpoint_projection_stage_diagnostic_v1"
PUBLICATION_SCHEMA = "original_compact_C3_quiet_endpoint_projection_publication_v1"
SEAL_SCHEMA = "compact_C3_quiet_endpoint_projection_final_seal_v1"

KERNEL_RELATIVE = Path(
    "src/openai_ns_reconstruction/original_compact_C3_quiet_endpoint_projection.py"
)
RUNNER_RELATIVE = Path("scripts/run_original_compact_C3_quiet_endpoint_projection.py")
REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_C3_quiet_endpoint_projection.json"
)
MILESTONE_RELATIVE = Path(
    "docs/CURRENT_ORIGINAL_COMPACT_C3_QUIET_ENDPOINT_PROJECTION_2026_10_26.md"
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

# These are authenticated raw artifacts from the immediately preceding quiet
# attempt.  They are references/data only; the manager never imports the old
# kernel or replays the failed attempt.
OLD_ACCEPTED_RESULT_SHA256 = (
    "b91da155898c38a9ca4b26595983a38aa04b2965ae3db0a908fadcea97b34241"
)
FAILED_RESULT_SHA256 = (
    "e9b1711e64408b1d016d2a0e73adc55deef8db64e70c76acbc9f0b111b97790d"
)
FAILED_SOURCE_SHA256 = (
    "0bf78d19e496c44bfaffd61b4c070963f5b5aca2bb5ff22bbdf1d7afd46024f1"
)
QUIET_INPUT_SHA256 = (
    "6783a68cf00b575c717d64e6a157b87f9678ffde09b298cb9add34ebe7e4fd78"
)
NEXT_CONTRACT_PATH = WORK / "compact_nonquiet_predecessor_static_contract_2026_10_26.json"
NEXT_CONTRACT_SHA256 = (
    "6e34296a2f79c4020f3279f13488e7532cb5c6f7814621ca25b1465a600a9297"
)
FROZEN_KERNEL_SHA256 = (
    "baf61439c1587e481029d8b25d0822b7e209a482d4d72f47fb712cc2ccc44383"
)

PARENT_SOURCE_PIN_COUNT = 2179
PARENT_CHECKPOINT_COUNT = 232
CURRENT_CHECKPOINT_COUNT = 236
PARENT_FAILED_COUNT = 169
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26
TERMINAL_QUERY_COUNT = 20
MEMORY_QUERY_COUNT = 6
TRANSPORT_RECIPE_COUNT = 20

NEW_FLAGS = (
    "C3_quiet_terminal_query_projection_bound",
    "C3_quiet_endpoint_memory_relations_reused",
)
QUIET_FLAGS = (
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
)
FORBIDDEN_TRUE_FLAGS = (
    "global_five_moment_closure",
    "full_goal_complete",
    "numeric_values_computed",
    "native_old_density_closure",
    "native_old_density_identity",
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
        "runtime_model_metadata",
        "original_source_evaluations",
        "source_packet_defined",
        "source_packet_definition_accepted",
        "kernel_imported",
        "derive_entry_calls",
        "derive_calls",
        "entry_attempts",
        "source_checkpoint_sha256",
        "entry_checkpoint_sha256",
        "diagnostic_checkpoint_sha256",
    }
)

PARENT_REPORT_EXCLUDED_FIELDS = (
    "metadata_helper_sources",
    "lossless_artifacts",
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
    require(path.exists(), "Missing artifact: " + str(path))
    raw = path.read_bytes()
    if expected is not None:
        require(sha_bytes(raw) == expected.lower(), "Changed artifact: " + str(path))
    return raw


def read_json(path, expected=None):
    raw = read_bytes(path, expected)
    return raw, json.loads(raw.decode("utf-8"), parse_int=int)


def resolve_path(value, base=REPO):
    path = Path(value)
    if path.is_absolute():
        return path
    direct = base / path
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
    view = {
        key: report[key]
        for key in sorted(report)
        if key not in set(PARENT_REPORT_EXCLUDED_FIELDS)
    }
    view["parent_report_is_authenticated_view"] = True
    view["parent_report_view_excluded_fields"] = list(PARENT_REPORT_EXCLUDED_FIELDS)
    require(len(encode(view)) < 100_000, "Parent report metadata view exceeds 100 KB")
    return view


def report_result_reference(report, report_path, report_raw_sha):
    selector = report.get("full_result_artifact_selector")
    result_sha = report.get("full_result_raw_sha256")
    if isinstance(selector, str) and selector.startswith("$.lossless_artifacts[path=") and selector.endswith("]"):
        result_name = selector[len("$.lossless_artifacts[path=") : -1]
        require(valid_sha(result_sha), "Current report result SHA is malformed")
        result_path = resolve_path(result_name)
        result_raw, result = read_json(result_path, result_sha)
        require(report.get("full_result") == result, "Publication report full result is not exact")
        require(report.get("full_result_raw_bytes") == len(result_raw), "Publication report result bytes changed")
        return {
            "path": str(result_path),
            "sha256": result_sha.lower(),
            "bytes": len(result_raw),
            "schema": result.get("schema"),
            "component": result.get("component"),
        }, result
    if report.get("schema") == RESULT_SCHEMA and report.get("component") == COMPONENT:
        return {
            "path": str(report_path),
            "sha256": report_raw_sha,
            "bytes": report_path.stat().st_size,
            "schema": report.get("schema"),
            "component": report.get("component"),
        }, report
    raise AssertionError("Current failed report has no exact lossless result selector")


def load_parent(args):
    require(args.parent_report and args.parent_report_sha256 and args.parent_commit, "Parent report, SHA, and commit are required")
    report_path = resolve_path(args.parent_report, WORK)
    report_raw, report = read_json(report_path, args.parent_report_sha256)
    require(isinstance(report, dict), "Parent report is not a mapping")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT, "Parent source-pin count changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Parent checkpoint count changed")
    require(len(report.get("false_flags", [])) == PARENT_FALSE_FLAG_COUNT, "Parent false-flag count changed")
    require(report.get("preserved_failed_artifact_count", len(report.get("preserved_failed_artifact_hashes_union", {}))) == PARENT_FAILED_COUNT, "Parent failed-artifact count changed")
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
    result_ref, result = report_result_reference(report, report_path, args.parent_report_sha256.lower())
    report_meta = {
        "path": str(report_path),
        "sha256": args.parent_report_sha256.lower(),
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
                and key not in QUIET_FLAGS
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
        "parent_commit": args.parent_commit,
    }


def locate_fixed(repo, relative, expected):
    path = resolve_path(relative, repo)
    return path, read_bytes(path, expected)


def load_failed(args):
    failed_report_path = resolve_path(args.current_failed_report, WORK) if args.current_failed_report else None
    if failed_report_path is None:
        failed_report_path = REPO / "tmp" / "original_compact_quiet_power_history_pair" / "result.json"
        failed_report_sha = FAILED_RESULT_SHA256
    else:
        failed_report_sha = args.current_failed_report_sha256
        require(valid_sha(failed_report_sha), "Current failed report SHA is required")
    raw, report_or_result = read_json(failed_report_path, failed_report_sha)
    result_ref, result = report_result_reference(report_or_result, failed_report_path, failed_report_sha.lower())
    require(result.get("passed") is False, "Current failed attempt is not failed")
    require(result.get("runtime_model_metadata") is None, "Current failed runtime metadata is not null")
    require(result_ref["sha256"] == FAILED_RESULT_SHA256, "Current failed result SHA changed")
    for key in QUIET_FLAGS:
        require(result.get(key) is False, "Current failed quiet flag changed: " + key)
    source_path = REPO / "tmp" / "original_compact_quiet_power_history_pair" / "source.json"
    source_raw, source = read_json(source_path, FAILED_SOURCE_SHA256)
    input_path = WORK / "compact_quiet_power_history_pair_input_2026_10_26.json"
    input_raw = read_bytes(input_path, QUIET_INPUT_SHA256)
    old_parent_path = REPO / "tmp" / "original_compact_signed_history_defect_transport" / "result.json"
    old_parent_raw = read_bytes(old_parent_path, OLD_ACCEPTED_RESULT_SHA256)
    quiet_input = json.loads(input_raw.decode("utf-8"), parse_int=int)
    frozen_input_slice = quiet_input.get("C3_literal_slice")
    old_parent_result = quiet_input.get("parent_result")
    require(isinstance(frozen_input_slice, dict), "Frozen C3 literal input slice is missing")
    require(
        frozen_input_slice.get("archive_raw_sha256")
        == "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193",
        "Frozen C3 input archive SHA changed",
    )
    require(
        isinstance(old_parent_result, dict)
        and old_parent_result.get("_raw_sha256") == OLD_ACCEPTED_RESULT_SHA256,
        "Old accepted parent result is not the pinned result",
    )
    embedded_parent = source.get("parent_result")
    require(isinstance(embedded_parent, dict), "Failed source checkpoint lacks embedded parent result")
    require(source.get("parent_result_sha256") == OLD_ACCEPTED_RESULT_SHA256, "Failed source parent-result SHA is not the fixed accepted result")
    parent_without_marker = dict(old_parent_result)
    parent_without_marker.pop("_raw_sha256", None)
    require(embedded_parent == parent_without_marker, "Failed source embedded parent result differs from the fixed accepted result")
    source_hashes = source.get("source_union_hashes", {})
    input_names = {str(input_path), input_path.as_posix(), str(input_path).replace("\\", "/")}
    require(any(source_hashes.get(name) == QUIET_INPUT_SHA256 for name in input_names), "Failed source does not cross-bind the fixed quiet input SHA")
    terminal = result.get("C3_source_binding", {}).get("terminal_row_checks", [])
    memory = result.get("C3_source_binding", {}).get("memory_row_checks", [])
    require(isinstance(terminal, list) and len(terminal) == TERMINAL_QUERY_COUNT, "Saved C3 terminal checks are not 20 rows")
    require(isinstance(memory, list) and len(memory) == MEMORY_QUERY_COUNT, "Saved C3 memory checks are not 6 rows")
    require(all(row.get("passed") is True for row in terminal), "Saved C3 terminal rows are not all passed")
    require(all(row.get("passed") is True for row in memory), "Saved C3 memory rows are not all passed")
    parent_result = result.get("parent_result_retained_exact")
    terminal_ports = parent_result.get("C3_rebuilt_history_source_ports", {}).get("terminal_ports") if isinstance(parent_result, dict) else None
    require(isinstance(terminal_ports, dict), "Old accepted parent C3 terminal ports are missing")
    require(sum(len(v) for v in terminal_ports.values() if isinstance(v, list)) == TERMINAL_QUERY_COUNT, "Old accepted C3 terminal query count changed")
    failed_meta = {
        "path": str(failed_report_path),
        "sha256": failed_report_sha.lower(),
        "bytes": len(raw),
        "schema": report_or_result.get("schema"),
        "component": report_or_result.get("component"),
    }
    return {
        "report": report_or_result,
        "report_meta": failed_meta,
        "result": result,
        "result_meta": result_ref,
        "source_meta": {"path": str(source_path), "sha256": FAILED_SOURCE_SHA256, "bytes": len(source_raw), "schema": source.get("schema"), "component": source.get("component")},
        "input_meta": {"path": str(input_path), "sha256": QUIET_INPUT_SHA256, "bytes": len(input_raw)},
        "old_parent_result_meta": {
            "path": str(old_parent_path),
            "sha256": OLD_ACCEPTED_RESULT_SHA256,
            "bytes": len(old_parent_raw),
        },
        "terminal_rows": clone(terminal),
        "memory_rows": clone(memory),
        "terminal_ports": clone(terminal_ports),
        "old_parent_result": clone(old_parent_result),
        "frozen_input_slice": clone(frozen_input_slice),
    }


def kernel_path(args):
    path = resolve_path(args.kernel_path or str(REPO / KERNEL_RELATIVE), REPO)
    read_bytes(path, FROZEN_KERNEL_SHA256)
    return path


def ensure_runner(args):
    path = resolve_path(args.runner_path or str(REPO / RUNNER_RELATIVE), REPO)
    if not path.exists():
        write_new(path, Path(__file__).read_bytes())
    raw = read_bytes(path)
    return path, raw


def source_paths(parent, failed, args, kernel, runner):
    external = Path(__file__).resolve()
    read_bytes(NEXT_CONTRACT_PATH, NEXT_CONTRACT_SHA256)
    stage = {
        parent["report_meta"]["path"]: parent["report_meta"]["sha256"],
        parent["result_meta"]["path"]: parent["result_meta"]["sha256"],
        failed["report_meta"]["path"]: failed["report_meta"]["sha256"],
        failed["result_meta"]["path"]: failed["result_meta"]["sha256"],
        failed["old_parent_result_meta"]["path"]: OLD_ACCEPTED_RESULT_SHA256,
        failed["source_meta"]["path"]: FAILED_SOURCE_SHA256,
        failed["input_meta"]["path"]: QUIET_INPUT_SHA256,
        str(kernel): sha_path(kernel),
        str(runner): sha_path(runner),
        str(external): sha_path(external),
        str(NEXT_CONTRACT_PATH): NEXT_CONTRACT_SHA256,
    }
    return stage, merge_hashes(parent["source_pins"], stage)


def projection_input(parent, failed):
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "projection_scope": {
            "field": "C3_quiet_endpoint",
            "typed_terminal_query_count": TERMINAL_QUERY_COUNT,
            "transport_recipe_count": TRANSPORT_RECIPE_COUNT,
            "memory_relation_count": MEMORY_QUERY_COUNT,
            "ordinary_Z_orders": [0, 1, 2, 3],
            "coordinate_interval": ["0", "2"],
            "ordinary_Z_domain": ["-1", "1"],
            "outgoing_selector_claimed": False,
            "p_out_equals_in_alias": True,
            "quiet_density_claimed": False,
            "native_claimed": False,
            "pair_proof_claimed": False,
        },
        "parent_result": clone(failed["old_parent_result"]),
        "failed_result": clone(failed["result"]),
        "frozen_input_slice": clone(failed["frozen_input_slice"]),
        "parent_result_raw_sha256": OLD_ACCEPTED_RESULT_SHA256,
        "failed_result_raw_sha256": FAILED_RESULT_SHA256,
        "frozen_input_slice_archive_sha256": failed["frozen_input_slice"].get("archive_raw_sha256"),
        "parent_report_metadata_only": clone(parent["report_view"]),
        "old_accepted_parent_result_raw_sha256": OLD_ACCEPTED_RESULT_SHA256,
        "failed_quiet_result_raw_sha256": FAILED_RESULT_SHA256,
        "failed_quiet_source_checkpoint_raw_sha256": FAILED_SOURCE_SHA256,
        "failed_quiet_input_raw_sha256": QUIET_INPUT_SHA256,
        "metadata_only": True,
        "source_only": True,
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
    }


def build_source(args):
    parent = load_parent(args)
    failed = load_failed(args)
    kernel = kernel_path(args)
    runner, runner_raw = ensure_runner(args)
    stage, source_union = source_paths(parent, failed, args, kernel, runner)
    payload = projection_input(parent, failed)
    false_flags = parent["false_flags"]
    inherited = dict(parent["accepted_true"])
    for key in (*QUIET_FLAGS, *NEW_FLAGS):
        inherited.pop(key, None)
    source = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_only": True,
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "parent_commit": parent["parent_commit"],
        "parent_report": parent["report_view"],
        "parent_report_metadata": parent["report_meta"],
        "parent_report_sha256": parent["report_meta"]["sha256"],
        "parent_result_metadata": failed["old_parent_result_meta"],
        "parent_result_raw_sha256": OLD_ACCEPTED_RESULT_SHA256,
        "failed_result_metadata": failed["result_meta"],
        "current_failed_report_metadata": failed["report_meta"],
        "current_failed_result_metadata": failed["result_meta"],
        "current_failed_source_checkpoint_metadata": failed["source_meta"],
        "current_failed_input_metadata": failed["input_meta"],
        "old_accepted_parent_result_raw_sha256": OLD_ACCEPTED_RESULT_SHA256,
        "failed_quiet_result_raw_sha256": FAILED_RESULT_SHA256,
        "failed_quiet_source_checkpoint_raw_sha256": FAILED_SOURCE_SHA256,
        "failed_quiet_input_raw_sha256": QUIET_INPUT_SHA256,
        "input": payload,
        # The fresh kernel receives these exact data objects directly as well
        # as under ``input``; they are JSON metadata, never encoded blobs.
        "parent_result": clone(payload["parent_result"]),
        "failed_result": clone(payload["failed_result"]),
        "frozen_input_slice": clone(payload["frozen_input_slice"]),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "stage_hashes": stage,
        "current_code_hashes": {
            str(kernel): sha_path(kernel),
            str(runner): sha_bytes(runner_raw),
        },
        "kernel_path": str(kernel),
        "runner_path": str(runner),
        "external_draft": {"path": str(Path(__file__).resolve()), "sha256": sha_path(Path(__file__).resolve())},
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "current_checkpoint_count": CURRENT_CHECKPOINT_COUNT,
        "preserved_failed_artifact_hashes_union": parent["failed"],
        "preserved_failed_artifact_count": len(parent["failed"]),
        "protected_file_hashes": parent["protected"],
        "legacy_dirty_path_hashes": parent["dirty"],
        "complete_parent_checkpoint_hashes": parent["complete_checkpoints"],
        "complete_parent_checkpoint_count": len(parent["complete_checkpoints"]),
        "inherited_false_flags": false_flags,
        "inherited_parent_true_flags": inherited,
        "new_flags": list(NEW_FLAGS),
        "quiet_flags_explicit_false": {key: False for key in QUIET_FLAGS},
        "projection_counts": {
            "typed_terminal_query_count": TERMINAL_QUERY_COUNT,
            "transport_recipe_count": TRANSPORT_RECIPE_COUNT,
            "memory_relation_count": MEMORY_QUERY_COUNT,
        },
        "next_executable_task": next_task_metadata(),
        "next_contract_metadata": {
            "path": str(NEXT_CONTRACT_PATH),
            "sha256": NEXT_CONTRACT_SHA256,
            "is_proof": False,
            "source_trace_only": True,
        },
        "goal_active": True,
        "full_goal_complete": False,
        "ancestor_decode_forbidden": True,
        "kernel_import_allowed_only_in_run": True,
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    source["input_digest"] = sha_bytes(encode(source["input"]))
    return source, parent, failed, kernel, runner


def prepare(args):
    require(not any(path.exists() for path in CHECKPOINTS), "Projection checkpoint directory already contains an attempt")
    source, parent, failed, kernel, runner = build_source(args)
    raw = write_new(SOURCE_CHECKPOINT, source)
    main_brief = update_main_brief(
        {
            "phase": "prepared_unentered",
            "source_checkpoint": str(SOURCE_CHECKPOINT),
            "source_checkpoint_sha256": sha_bytes(raw),
            "report": str(REPO / REPORT_RELATIVE),
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
        "source_checkpoint": str(SOURCE_CHECKPOINT),
        "source_checkpoint_sha256": sha_bytes(raw),
        "source_union_count": source["source_union_count"],
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "current_checkpoint_count": CURRENT_CHECKPOINT_COUNT,
        "typed_terminal_query_count": TERMINAL_QUERY_COUNT,
        "transport_recipe_count": TRANSPORT_RECIPE_COUNT,
        "memory_relation_count": MEMORY_QUERY_COUNT,
        "kernel_sha256": source["current_code_hashes"][str(kernel)],
        "runner_sha256": source["current_code_hashes"][str(runner)],
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
        "source_checkpoint": str(SOURCE_CHECKPOINT),
        "source_checkpoint_sha256": receipt["source_checkpoint_sha256"],
        "runner_sha256": receipt["runner_sha256"],
        "kernel_sha256": receipt["kernel_sha256"],
        "goal_active": True,
        "full_goal_complete": False,
        "new_flags": list(NEW_FLAGS),
        "next_executable_task": next_task_metadata(),
        "typed_terminal_query_count": TERMINAL_QUERY_COUNT,
        "transport_recipe_count": TRANSPORT_RECIPE_COUNT,
        "memory_relation_count": MEMORY_QUERY_COUNT,
        "runtime_model_metadata": None,
    }
    brief_raw = write_replace(BRIEF, brief)
    return {"prepare": receipt, "prepare_sha256": sha_bytes(receipt_raw), "brief_sha256": sha_bytes(brief_raw), "parent": parent["report_meta"], "failed": failed["result_meta"], "kernel": str(kernel)}


def load_source():
    source_raw, source = read_json(SOURCE_CHECKPOINT)
    require(source.get("schema") == INPUT_SCHEMA and source.get("component") == COMPONENT, "Projection source checkpoint schema changed")
    require(source.get("metadata_only") is True and source.get("runtime_model_metadata") is None, "Projection source checkpoint scope changed")
    require(source.get("parent_result_raw_sha256") == OLD_ACCEPTED_RESULT_SHA256, "Source parent-result SHA changed")
    require(source.get("failed_quiet_result_raw_sha256") == FAILED_RESULT_SHA256, "Source failed-result SHA changed")
    require(source.get("failed_quiet_source_checkpoint_raw_sha256") == FAILED_SOURCE_SHA256, "Source failed-source SHA changed")
    require(source.get("failed_quiet_input_raw_sha256") == QUIET_INPUT_SHA256, "Source quiet-input SHA changed")
    for name, digest in source.get("source_union_hashes", {}).items():
        read_bytes(resolve_path(name), digest)
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "Projection source union changed")
    require(source.get("complete_parent_checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Inherited checkpoint count changed")
    require(source.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT, "Inherited failed-artifact count changed")
    require(len(source.get("protected_file_hashes", {})) == PARENT_PROTECTED_COUNT, "Inherited protected count changed")
    require(len(source.get("legacy_dirty_path_hashes", {})) == PARENT_DIRTY_COUNT, "Inherited dirty count changed")
    return source_raw, source


def journal_append(events, event):
    payload = clone(event) if isinstance(event, dict) else {"event": str(event)}
    payload.setdefault("runtime_model_metadata", None)
    events.append(payload)
    with DIAGNOSTIC_CHECKPOINT.open("ab") as stream:
        stream.write(encode(payload))
        stream.flush()
        os.fsync(stream.fileno())


def fallback(source, calls, imported, error=None):
    result = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "metadata_only": False,
        "source_only": True,
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_defined": False,
        "kernel_imported": bool(imported),
        "derive_entry_calls": calls,
        "derive_calls": calls,
        "entry_attempts": 1,
        "input_digest": source.get("input_digest"),
        "blockers": [str(error)] if error else ["projection derive returned no result"],
        "projection_scope": clone(source.get("input", {}).get("projection_scope", {})),
        "projection_counts": {
            "typed_terminal_query_count": 0,
            "transport_recipe_count": 0,
            "memory_relation_count": 0,
        },
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "new_C3_terminal_query_rows": 0,
            "new_C3_endpoint_memory_rows": 0,
            "new_C3_quiet_density_rows": 0,
            "new_native_quiet_density_rows": 0,
            "new_common_kernel_rows": 0,
            "new_signed_quiet_transport_rows": 0,
        },
    }
    return result


def normalise_result(value, source, calls, imported):
    result = clone(value) if isinstance(value, dict) else fallback(source, calls, imported, "derive returned a non-mapping result")
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("passed", False)
    result["metadata_only"] = False
    result["source_only"] = True
    result["runtime_model_metadata"] = None
    result.setdefault("original_source_evaluations", 0)
    result["input_digest"] = source.get("input_digest")
    result["kernel_imported"] = bool(imported)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    result["source_packet_defined"] = bool(result.get("passed") is True)
    result["source_packet_definition_accepted"] = bool(result.get("passed") is True)
    result["parent_result_raw_sha256"] = OLD_ACCEPTED_RESULT_SHA256
    result["failed_quiet_result_raw_sha256"] = FAILED_RESULT_SHA256
    result["new_flags"] = list(NEW_FLAGS)
    result.setdefault("memory_relation_count", MEMORY_QUERY_COUNT)
    result.setdefault("projection_counts", {
        "typed_terminal_query_count": result.get("typed_terminal_query_count", 0),
        "transport_recipe_count": result.get("transport_recipe_count", 0),
        "memory_relation_count": MEMORY_QUERY_COUNT,
    })
    audit = result.setdefault("source_inspection_audit", {})
    if not isinstance(audit, dict):
        audit = {}
        result["source_inspection_audit"] = audit
    audit.setdefault("original_source_evaluations", 0)
    audit.setdefault("new_C3_terminal_query_rows", result.get("typed_terminal_query_count", 0))
    audit.setdefault("new_C3_endpoint_memory_rows", MEMORY_QUERY_COUNT)
    for key in (
        "new_C3_quiet_density_rows",
        "new_native_quiet_density_rows",
        "new_common_kernel_rows",
        "new_signed_quiet_transport_rows",
    ):
        audit.setdefault(key, 0)
    for key, value in source.get("inherited_parent_true_flags", {}).items():
        if value is True and key not in BOOKKEEPING_FIELDS and key not in (*QUIET_FLAGS, *NEW_FLAGS):
            result.setdefault(key, True)
    for key in source.get("inherited_false_flags", []):
        result[key] = False
    for key in QUIET_FLAGS:
        result[key] = False
    for key in FORBIDDEN_TRUE_FLAGS:
        result[key] = False
    result["goal_active"] = True
    result["full_goal_complete"] = False
    return result


def _find_projection_counts(result):
    nodes = [result]
    for key in ("projection", "C3_quiet_endpoint_projection", "endpoint_projection", "projection_certificate"):
        value = result.get(key)
        if isinstance(value, dict):
            nodes.append(value)
    terminal = None
    transport = None
    memory = None
    for node in nodes:
        for key in ("projected_terminal_count", "typed_terminal_query_count", "terminal_query_count", "terminal_row_count", "terminal_checks"):
            value = node.get(key)
            if isinstance(value, int) and not isinstance(value, bool):
                terminal = value
            elif isinstance(value, list):
                terminal = len(value)
        for key in ("transport_recipe_count", "projected_transport_recipe_count", "transport_recipes"):
            value = node.get(key)
            if isinstance(value, int) and not isinstance(value, bool):
                transport = value
            elif isinstance(value, list):
                transport = len(value)
        for key in ("projected_memory_count", "memory_relation_count", "memory_row_count", "memory_checks"):
            value = node.get(key)
            if isinstance(value, int) and not isinstance(value, bool):
                memory = value
            elif isinstance(value, list):
                memory = len(value)
    return terminal, transport, memory


def validate_result(result, source):
    require(result.get("schema") == RESULT_SCHEMA, "Projection result schema changed")
    require(result.get("component") == COMPONENT, "Projection result component changed")
    require(result.get("metadata_only") is False and result.get("source_only") is True, "Projection result scope changed")
    require(result.get("runtime_model_metadata") is None, "Projection runtime metadata is not null")
    require(result.get("original_source_evaluations") == 0, "Projection evaluated an original source")
    require(result.get("derive_entry_calls") == 1 and result.get("derive_calls") == 1 and result.get("entry_attempts") == 1, "Projection derive bookkeeping changed")
    require(result.get("input_digest") == source.get("input_digest"), "Projection input digest changed")
    if result.get("passed") is True:
        require(result.get("source_packet_definition_accepted") is True, "Passed projection lacks owned source-packet acceptance")
    else:
        require(result.get("source_packet_definition_accepted") is False, "Failed projection owns an accepted source-packet flag")
    require(result.get("full_goal_complete") is False, "Projection overclaims full goal completion")
    require(result.get("numeric_values_computed") is not True, "Projection computed numeric values")
    for key in QUIET_FLAGS:
        require(result.get(key) is False, "Quiet density/pair flag was not kept false: " + key)
    for key in FORBIDDEN_TRUE_FLAGS:
        require(result.get(key) is not True, "Projection overclaims forbidden flag: " + key)
    audit = result.get("source_inspection_audit", {})
    require(isinstance(audit, dict), "Projection source audit is not a mapping")
    for key in ("new_C3_quiet_density_rows", "new_native_quiet_density_rows", "new_common_kernel_rows", "new_signed_quiet_transport_rows"):
        require(audit.get(key, 0) == 0, "Projection recorded forbidden quiet computation: " + key)
    terminal, transport, memory = _find_projection_counts(result)
    if result.get("passed") is False:
        require(result.get(NEW_FLAGS[0]) is not True and result.get(NEW_FLAGS[1]) is not True, "Failed projection overclaims a new bound")
        return {
            "passed": False,
            "honest_failure": True,
            "terminal_count": terminal,
            "transport_recipe_count": transport,
            "memory_count": memory,
        }
    require(terminal == TERMINAL_QUERY_COUNT, "Projection terminal count is not 20")
    require(transport == TRANSPORT_RECIPE_COUNT, "Projection transport-recipe count is not 20")
    require(memory == MEMORY_QUERY_COUNT, "Projection memory count is not 6")
    require(result.get(NEW_FLAGS[0]) is True and result.get(NEW_FLAGS[1]) is True, "Projection flags are not both source-bound")
    outgoing = result.get("outgoing_ports")
    if outgoing is not None:
        require(isinstance(outgoing, dict) and all(value in ([], None) for value in outgoing.values()), "Projection introduced an outgoing selector")
    return {"passed": True, "terminal_count": terminal, "transport_recipe_count": transport, "memory_count": memory}


def run(args):
    require(not ENTERED_CHECKPOINT.exists() and not RESULT_CHECKPOINT.exists() and not DIAGNOSTIC_CHECKPOINT.exists(), "Projection attempt already entered")
    source_raw, source = load_source()
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "source_checkpoint_path": str(SOURCE_CHECKPOINT),
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
        kernel_raw = read_bytes(kernel, source["current_code_hashes"][str(kernel)])
        spec = importlib.util.spec_from_file_location("_compact_C3_quiet_endpoint_projection_entry", kernel)
        require(spec is not None and spec.loader is not None, "Cannot create projection kernel import specification")
        module = importlib.util.module_from_spec(spec)
        imported = True
        spec.loader.exec_module(module)
        derive = getattr(module, "derive", None)
        require(callable(derive), "Projection kernel lacks derive(source, journal_callback)")

        kernel_events = []
        calls += 1
        result = derive(clone(source), kernel_events)
        for event in kernel_events:
            journal_append(events, event)
    except Exception as error:
        journal_append(events, {"stage": "derive_exception", "error_type": type(error).__name__, "error": str(error)})
        result = fallback(source, calls, imported, error)
    result = normalise_result(result, source, calls, imported)
    result["source_checkpoint_sha256"] = sha_bytes(source_raw)
    result["entry_checkpoint_sha256"] = sha_bytes(entered_raw)
    result["diagnostic_checkpoint_sha256"] = sha_path(DIAGNOSTIC_CHECKPOINT)
    result["diagnostic_event_count"] = len(events)
    # Durable result first; validation is intentionally after this fsync.
    result_raw = write_new(RESULT_CHECKPOINT, result)
    validation = None
    validation_error = None
    try:
        validation = validate_result(result, source)
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
    entry = {"path": str(path), "sha256": sha_bytes(raw), "bytes": len(raw)}
    if embed:
        entry["raw_bytes_base64"] = base64.b64encode(raw).decode("ascii")
        entry["raw_bytes_exact"] = True
    return entry


def checkpoint_artifacts(embed=False):
    return [artifact_entry(path, embed=embed) for path in CHECKPOINTS]


def update_main_brief(pointer):
    """Merge a compact component pointer into the shared resume brief."""
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
    # Preparation must not make the shared brief point at a milestone or
    # report that does not exist yet. Publication can point at the newly
    # written files; resume/seal additionally promote a verified authority.
    if pointer.get("phase") != "prepared_unentered":
        milestone = {
            "path": str(REPO / MILESTONE_RELATIVE),
            "report": str(REPO / REPORT_RELATIVE),
            "phase": pointer.get("phase"),
        }
        if pointer.get("milestone_sha256") is not None:
            milestone["sha256"] = pointer["milestone_sha256"]
        if pointer.get("report_sha256") is not None:
            milestone["report_sha256"] = pointer["report_sha256"]
        state["latest_milestone"] = milestone
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
            keep = {key for key in (COMPONENT, "compact_quiet_power_history_pair") if key in old_in_progress}
            omitted = sorted(key for key in old_in_progress if key not in keep)
            compact["in_progress"] = {key: old_in_progress[key] for key in sorted(keep)}
            compact["archived_in_progress_selector"] = {
                "selector": "$.in_progress",
                "path": str(preimage),
                "sha256": sha_bytes(old_raw),
                "bytes": len(old_raw),
                "omitted_keys": omitted,
                "preserved_exact_bytes": True,
            }
        compact["archived_main_resume_brief_preimage"] = {
            "path": str(preimage),
            "sha256": sha_bytes(old_raw),
            "bytes": len(old_raw),
        }
        compact["main_resume_brief_archival_compaction"] = True
        merged_raw = encode(compact)
    require(len(merged_raw) <= 100_000, "Shared RESUME_BRIEF remains over 100 KB after archival compaction")
    write_replace(MAIN_BRIEF, merged_raw)
    return {
        "path": str(MAIN_BRIEF),
        "sha256": sha_bytes(merged_raw),
        "bytes": len(merged_raw),
        "preimage": {"path": str(preimage), "sha256": sha_path(preimage), "bytes": len(old_raw)},
    }


def next_task_metadata():
    return {
        "action": (
            "Bind and evaluate actual signed driver recipes on the preceding "
            "O3 transition; retain the unknown inlet and all caps."
        ),
        "contract": str(NEXT_CONTRACT_PATH),
        "contract_sha256": NEXT_CONTRACT_SHA256,
        "contract_is_proof": False,
        "contract_is_source_trace_only": True,
        "do_not_promote_driver_facts_without_new_source_binding": True,
        "do_not_reenter_current_attempt": True,
    }


def authoritative_pointer(report, report_raw, milestone_raw, commit):
    full = report.get("full_result")
    require(isinstance(full, dict), "Published report full result is not a mapping")
    return {
        "schema": COMPONENT + "_authoritative_pointer_v1",
        "commit": commit,
        "report": str(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "milestone": str(REPO / MILESTONE_RELATIVE),
        "milestone_sha256": sha_bytes(milestone_raw),
        "attempt_status": report.get("attempt_status"),
        "passed": report.get("passed") is True,
        "source_packet_defined": full.get("source_packet_defined") is True,
        "source_packet_definition_accepted": full.get("source_packet_definition_accepted") is True,
        "new_flags": clone(report.get("new_flags", {})),
        "quiet_flags": clone(report.get("quiet_flags", {})),
        "source_pin_count": report.get("source_pin_count"),
        "checkpoint_count": report.get("checkpoint_count"),
        "preserved_failed_artifact_count": report.get("preserved_failed_artifact_count"),
        "runtime_model_metadata": None,
        "full_goal_complete": False,
        "do_not_reenter": True,
    }


def hydrate_publish_args(args, source):
    parent_meta = source.get("parent_report_metadata", {})
    failed_meta = source.get("current_failed_report_metadata", {})
    if not args.parent_report and isinstance(parent_meta, dict):
        args.parent_report = parent_meta.get("path")
    if not args.parent_report_sha256 and isinstance(parent_meta, dict):
        args.parent_report_sha256 = parent_meta.get("sha256")
    if not args.parent_commit:
        args.parent_commit = source.get("parent_commit")
    if not args.current_failed_report and isinstance(failed_meta, dict):
        args.current_failed_report = failed_meta.get("path")
    if not args.current_failed_report_sha256 and isinstance(failed_meta, dict):
        args.current_failed_report_sha256 = failed_meta.get("sha256")


def report_for_publish(source, result, result_raw, parent, failed, code2, external, input_meta):
    current_ck = checkpoint_artifacts(embed=True)
    artifacts = [
        artifact_entry(Path(source["parent_report_metadata"]["path"]), source["parent_report_metadata"]["sha256"]),
        artifact_entry(Path(failed["old_parent_result_meta"]["path"]), failed["old_parent_result_meta"]["sha256"]),
        artifact_entry(Path(failed["result_meta"]["path"]), failed["result_meta"]["sha256"]),
        artifact_entry(Path(failed["source_meta"]["path"]), failed["source_meta"]["sha256"]),
        artifact_entry(Path(input_meta["path"]), input_meta["sha256"], embed=True),
        artifact_entry(Path(code2[0]["path"]), code2[0]["sha256"], embed=True),
        artifact_entry(Path(code2[1]["path"]), code2[1]["sha256"], embed=True),
        artifact_entry(external["path"], external["sha256"], embed=True),
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
        "parent_result_raw_sha256": OLD_ACCEPTED_RESULT_SHA256,
        "parent_result_metadata": clone(source["parent_result_metadata"]),
        "current_failed_result_raw_sha256": FAILED_RESULT_SHA256,
        "current_failed_result_metadata": clone(source["failed_result_metadata"]),
        "current_failed_source_checkpoint_raw_sha256": FAILED_SOURCE_SHA256,
        "current_failed_input_raw_sha256": QUIET_INPUT_SHA256,
        "source_pins": clone(source["source_union_hashes"]),
        "source_pin_count": len(source["source_union_hashes"]),
        "source_union_count": len(source["source_union_hashes"]),
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_inventory": clone(source["complete_parent_checkpoint_hashes"]),
        "checkpoint_hashes": {entry["path"]: entry["sha256"] for entry in current_ck},
        "checkpoint_paths": [entry["path"] for entry in current_ck],
        "checkpoint_count": len(source["complete_parent_checkpoint_hashes"]) + len(current_ck),
        "inherited_checkpoint_count": len(source["complete_parent_checkpoint_hashes"]),
        "preserved_failed_artifact_hashes_union": clone(source["preserved_failed_artifact_hashes_union"]),
        "preserved_failed_artifact_count": PARENT_FAILED_COUNT,
        "protected_file_hashes": clone(source["protected_file_hashes"]),
        "legacy_dirty_path_hashes": clone(source["legacy_dirty_path_hashes"]),
        "inherited_failed_artifacts_preserved": True,
        "false_flags": list(source["inherited_false_flags"]),
        "false_flag_values": {key: False for key in source["inherited_false_flags"]},
        "inherited_parent_true_flags": clone(source["inherited_parent_true_flags"]),
        "new_flags": {key: result.get(key) is True for key in NEW_FLAGS},
        "quiet_flags": {key: False for key in QUIET_FLAGS},
        "typed_terminal_query_count": TERMINAL_QUERY_COUNT,
        "transport_recipe_count": TRANSPORT_RECIPE_COUNT,
        "memory_relation_count": MEMORY_QUERY_COUNT,
        "projection_counts": {
            "typed_terminal_query_count": TERMINAL_QUERY_COUNT,
            "transport_recipe_count": TRANSPORT_RECIPE_COUNT,
            "memory_relation_count": MEMORY_QUERY_COUNT,
        },
        "full_result": clone(result),
        "component_result": clone(result),
        "component_result_exact": True,
        "source_packet_defined": result.get("passed") is True,
        "source_packet_definition_accepted": result.get("passed") is True,
        "original_source_evaluations": 0,
        "numeric_values_computed": False,
        "full_result_raw_sha256": sha_bytes(result_raw),
        "full_result_raw_bytes": len(result_raw),
        "full_result_artifact_selector": "$.lossless_artifacts[path=" + str(RESULT_CHECKPOINT) + "]",
        "lossless_artifacts": artifacts,
        "lossless_artifact_count": len(artifacts),
        "raw_bytes_lossless": True,
        "goal_active": True,
        "full_goal_complete": False,
        "next_executable_task": next_task_metadata(),
        "next_contract_metadata": {
            "path": str(NEXT_CONTRACT_PATH),
            "sha256": NEXT_CONTRACT_SHA256,
            "is_proof": False,
            "source_trace_only": True,
        },
        "code2": [clone(code2[0]), clone(code2[1])],
        "external_draft": clone(external),
        "math_reexecuted": False,
        "ancestor_kernel_imported": False,
    }
    for key in NEW_FLAGS:
        report[key] = result.get(key) is True
    for key, value in source["inherited_parent_true_flags"].items():
        if key not in BOOKKEEPING_FIELDS and key not in NEW_FLAGS and key not in QUIET_FLAGS:
            report.setdefault(key, True)
    for key in source["inherited_false_flags"]:
        report[key] = False
    for key in QUIET_FLAGS:
        report[key] = False
    for key in FORBIDDEN_TRUE_FLAGS:
        report[key] = False
    return report


def handoff_preimage(path, raw):
    preimage = WORK / f"{COMPONENT}_handoff_preimage_{path.stem}_{sha_bytes(raw)[:12]}.md"
    if not preimage.exists():
        write_new(preimage, raw)
    return {"path": str(preimage), "sha256": sha_path(preimage), "bytes": len(raw)}


def publish(args):
    source_raw, source = load_source()
    hydrate_publish_args(args, source)
    require(RESULT_CHECKPOINT.exists(), "Projection result checkpoint is missing")
    result_raw, result = read_json(RESULT_CHECKPOINT)
    parent = load_parent(args)
    failed = load_failed(args)
    kernel = resolve_path(source["kernel_path"], REPO)
    runner = resolve_path(source["runner_path"], REPO)
    code2 = [artifact_entry(kernel, source["current_code_hashes"][str(kernel)]), artifact_entry(runner, source["current_code_hashes"][str(runner)])]
    external = {"path": Path(source["external_draft"]["path"]), "sha256": source["external_draft"]["sha256"]}
    artifact_entry(external["path"], external["sha256"])
    input_meta = source["current_failed_input_metadata"]
    durable_validation = validate_result(result, source)
    report = report_for_publish(source, result, result_raw, parent, failed, code2, external, input_meta)
    report["durable_result_validation"] = clone(durable_validation)
    report_raw = write_replace(REPO / REPORT_RELATIVE, report)
    milestone = (
        f"# Current original C3 quiet endpoint projection (2026-10-26)\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}`; SHA-256 `{sha_bytes(report_raw)}`.\n\n"
        f"Attempt status: `{report['attempt_status']}`; source-only: `true`; typed terminal queries: `20/20`; transport recipes: `20/20`; cached endpoint memory relations: `6/6`.\n\n"
        "This leaf projects explicit parent C3 terminal queries onto cached passed rows. "
        "It claims no quiet density, native source, paired-history proof, outgoing selector, numeric value, physical completion, or full-goal completion. "
        "All inherited accepted flags, false flags, failed artifacts, caps, and checkpoint maps remain retained.\n\n"
        "The entered attempt is immutable and must never be reentered.\n\n"
        f"Next source-bound task: bind/evaluate actual signed driver recipes on the preceding O3 transition; "
        f"unproved metadata contract `{NEXT_CONTRACT_PATH}` (SHA-256 `{NEXT_CONTRACT_SHA256}`) is only a source-trace pointer.\n"
    ).encode("utf-8")
    milestone_raw = write_replace(REPO / MILESTONE_RELATIVE, milestone)
    preimages = []
    handoff_prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}` (SHA-256 `{sha_bytes(report_raw)}`).\n"
        f"Status: `{report['attempt_status']}`; terminal projection `20/20`; transport recipes `20/20`; memory reuse `6/6`; full goal remains incomplete.\n"
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
        f"/{KERNEL_RELATIVE.as_posix()} -diff",
        f"/{RUNNER_RELATIVE.as_posix()} -diff",
        f"/{REPORT_RELATIVE.as_posix()} -diff",
    ]
    for line in attr_lines:
        if line not in attrs_text.splitlines():
            attrs_text += ("\n" if attrs_text and not attrs_text.endswith("\n") else "") + line + "\n"
    write_replace(attrs_path, attrs_text.encode("utf-8"))
    main_brief = update_main_brief(
        {
            "phase": "published_passed" if report["attempt_status"] == "passed" else "published_failed",
            "report": str(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "milestone": str(REPO / MILESTONE_RELATIVE),
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
        "report": str(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "attempt_status": report["attempt_status"],
        "scoped_paths": [str(path.as_posix()) for path in SCOPED_RELATIVES],
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
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Projection publication checkpoint count changed")
    require(report.get("source_pin_count", 0) > PARENT_SOURCE_PIN_COUNT, "Projection source pins did not advance")
    require(report.get("full_goal_complete") is False and report.get("goal_active") is True, "Projection goal is not active/incomplete")
    head = git_output("rev-parse", "HEAD")
    require(valid_git_revision(head), "Local HEAD is not an exact Git revision")
    parent_commit = report.get("parent_commit")
    require(isinstance(parent_commit, str) and valid_git_revision(parent_commit), "Published parent commit is malformed")
    require(git_output("rev-parse", "HEAD^") == parent_commit, "HEAD parent does not match the published parent commit")
    branch = git_output("branch", "--show-current")
    require(branch, "Current Git branch is unnamed")
    remote_line = git_output("ls-remote", "origin", "refs/heads/" + branch)
    remote_parts = remote_line.split()
    require(len(remote_parts) >= 2 and valid_git_revision(remote_parts[0]), "Remote branch lookup returned no exact head")
    remote = remote_parts[0].lower()
    require(head.lower() == remote, "HEAD is not synchronized with origin's exact branch")
    changed = set(git_output("diff", "--name-only", parent_commit + "..HEAD").splitlines())
    expected = {path.as_posix() for path in SCOPED_RELATIVES}
    require(changed == expected, "Committed scope is not exactly the nine projection paths")
    for relative in SCOPED_RELATIVES:
        local_path = REPO / relative
        require(local_path.exists(), "Published scoped path is missing: " + relative.as_posix())
        blob_size = int(git_output("cat-file", "-s", "HEAD:" + relative.as_posix()))
        require(blob_size == local_path.stat().st_size, "Published cat-file byte count changed: " + relative.as_posix())
    return head, remote


def resume(args):
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "Projection publication report is missing")
    head, remote = verify_published_commit(report)
    milestone_raw = read_bytes(REPO / MILESTONE_RELATIVE)
    authority = authoritative_pointer(report, report_raw, milestone_raw, head)
    main_brief = update_main_brief(
        {
            "phase": "resumed_published",
            "report": str(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "milestone": str(REPO / MILESTONE_RELATIVE),
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
        }
    )
    receipt = {
        "schema": COMPONENT + "_resume_receipt_v1",
        "component": COMPONENT,
        "report": str(REPO / REPORT_RELATIVE),
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
    require(args.review and args.review_sha256, "Seal requires an independent review path and SHA-256")
    review_path = resolve_path(args.review, WORK)
    review_raw = read_bytes(review_path, args.review_sha256)
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "Projection report is not published")
    source_raw, source = load_source()
    require(valid_sha(report.get("full_result_raw_sha256")), "Published full-result SHA is malformed")
    result_raw, result = read_json(RESULT_CHECKPOINT, report["full_result_raw_sha256"])
    require(report.get("full_result") == result, "Published full result changed before seal")
    require(report.get("full_result_raw_bytes") == len(result_raw), "Published full-result byte count changed before seal")
    require(report.get("checkpoint_count") == CURRENT_CHECKPOINT_COUNT, "Published checkpoint count changed before seal")
    require(report.get("source_pin_count") == len(source.get("source_union_hashes", {})), "Published source-pin count changed before seal")
    durable_validation = validate_result(result, source)
    require(report.get("durable_result_validation") == durable_validation, "Published durable-result validation changed before seal")
    head, remote = verify_published_commit(report)
    milestone_raw = read_bytes(REPO / MILESTONE_RELATIVE)
    authority = authoritative_pointer(report, report_raw, milestone_raw, head)
    main_brief = update_main_brief(
        {
            "phase": "sealed",
            "report": str(REPO / REPORT_RELATIVE),
            "report_sha256": sha_bytes(report_raw),
            "seal": str(SEAL_RECEIPT),
            "goal_active": True,
            "full_goal_complete": False,
            "do_not_reenter": True,
            "next_executable_task": next_task_metadata(),
            "latest_authoritative_state": authority,
            "verified_local_and_remote_commit": head,
            "next_remaining": next_task_metadata(),
            "milestone_sha256": sha_bytes(milestone_raw),
        }
    )
    receipt = {
        "schema": SEAL_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "report": str(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "head": head,
        "remote": remote,
        "durable_result_validation": durable_validation,
        "independent_review": {"path": str(review_path), "sha256": sha_bytes(review_raw), "bytes": len(review_raw)},
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
    p.add_argument("--parent-report")
    p.add_argument("--parent-report-sha256")
    p.add_argument("--parent-commit")
    p.add_argument("--current-failed-report")
    p.add_argument("--current-failed-report-sha256")
    p.add_argument("--kernel-path")
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
