"""Standalone metadata workflow for the P2 pressure-gap leaf.

This file is a workflow draft, not the pressure-gap kernel.  The only
scientific entry point it knows about is the future repository kernel's
``derive(source, diagnostic)`` function.  ``prepare`` authenticates the
saved source-binding publication and freezes a small input specification;
``run`` writes the entered marker before importing that kernel, calls it once,
and saves its raw result before doing any validation.  Publication, resume,
and seal operate on saved bytes and metadata.  Importing this module performs
no filesystem, Git, kernel, or runner action.
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

COMPONENT = "compact_O3_transition_P2_pressure_gap"
INPUT_SCHEMA = "compact_O3_transition_P2_pressure_gap_input_v1"
RESULT_SCHEMA = "original_compact_O3_transition_P2_pressure_gap"
ENTERED_SCHEMA = "original_compact_O3_transition_P2_pressure_gap_entered_v1"
DIAGNOSTIC_SCHEMA = "compact_O3_transition_P2_pressure_gap_stage_diagnostic_v1"
PUBLICATION_SCHEMA = "original_compact_O3_transition_P2_pressure_gap_publication_v1"
SEAL_SCHEMA = "compact_O3_transition_P2_pressure_gap_final_seal_v1"

KERNEL_RELATIVE = Path(
    "src/openai_ns_reconstruction/"
    "original_compact_O3_transition_P2_pressure_gap.py"
)
RUNNER_RELATIVE = Path(
    "scripts/run_original_compact_O3_transition_P2_pressure_gap.py"
)
REPORT_RELATIVE = Path(
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_P2_pressure_gap.json"
)
MILESTONE_RELATIVE = Path(
    "docs/CURRENT_ORIGINAL_COMPACT_O3_TRANSITION_P2_PRESSURE_GAP_2026_10_27.md"
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

INPUT_SPEC = WORK / f"{COMPONENT}_input_2026_10_27.json"
BRIEF = WORK / f"{COMPONENT}_brief_2026_10_27.json"
MAIN_BRIEF = REPO.parents[1] / "RESUME_BRIEF_2026_10_16.json"
PREPARE_RECEIPT = WORK / f"{COMPONENT}_prepare_receipt_2026_10_27.json"
RUN_RECEIPT = WORK / f"{COMPONENT}_run_receipt_2026_10_27.json"
PUBLISH_RECEIPT = WORK / f"{COMPONENT}_publication_receipt_2026_10_27.json"
RESUME_RECEIPT = WORK / f"{COMPONENT}_resume_receipt_2026_10_27.json"
SEAL_RECEIPT = WORK / f"{COMPONENT}_final_seal_receipt_2026_10_27.json"

# This stage is based on the sealed source-binding publication, while that
# publication itself retains the preceding primitive-gap parent anchors.
PARENT_REPORT_DEFAULT = (
    REPO
    / "experiments/root_st073/"
    / "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_P2_source_binding.json"
)
PARENT_REPORT_SHA256 = (
    "ac6768f0a842b8a0c3dba0fb8c24721bdce7c3c1435591c9bfdf6709db9dde3f"
)
LEGACY_PARENT_REPORT_DEFAULT = (
    REPO
    / "experiments/root_st073/"
    / "lei_ren_part1_paper_compliant_current_original_compact_O3_transition_primitive_gap.json"
)
LEGACY_PARENT_REPORT_SHA256 = (
    "3261bb535e3987c7e74c02cef84871834f6fa2bdf894b0bc0c907f08095d1cd7"
)
PARENT_COMMIT = "b7687e8f187ab5e54552d9a3d28201d69a38476c"
PARENT_BINDING_REPORT_PARENT_HEAD = "45bad2f3d7bd2e98ab3cee137f64a3f5a390c289"
PARENT_BINDING_REPORT_COMMIT = "098c46878efd98c84ec2ef438c9d7f1596f529d4"
PARENT_BINDING_REPORT_PARENT_COMMIT = "b0f7e8e0d9a69b8c82345ff0a6a356c8b7a2cc96"
PARENT_BINDING_RESULT_SHA256 = (
    "007ffc9391381a6d62cfccb763dc6a4ce2b430cf51d6209013b576f2f36fc521"
)
PARENT_BINDING_SOURCE_CHECKPOINT_SHA256 = (
    "79a781ea0fa4a6509e907abb0f6b8c03f01f06ae1700bc9123e1f8064bf26a41"
)
PARENT_SOURCE_PIN_COUNT = 2228
PARENT_CHECKPOINT_COUNT = 248
CURRENT_CHECKPOINT_COUNT = PARENT_CHECKPOINT_COUNT + 4
PARENT_FAILED_COUNT = 176
FAILED_COUNT_IF_NEW_FAILURE = PARENT_FAILED_COUNT + 7
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2
PARENT_FALSE_FLAG_COUNT = 26

EXPECTED_COUNTS = {
    "new_mathematical_proofs": 1,
    "new_algebra_family_calls": 1,
    "pressure_identity_orders": 4,
    "gap_function_rows": 4,
    "side_quotient_premise_rows": 8,
    "numeric_calls": 0,
    "source_provider_calls": 0,
    "old_math_calls": 0,
}
COUNT_KEYS = tuple(EXPECTED_COUNTS)
NEW_FLAGS = (
    "O3_transition_P2_signed_jet_gap_functions_defined",
    "O3_transition_P2_pressure_quotient_correlation_ordinary_Z3_proved",
)
EXPECTED_PROJECTION_ORDERS = (0, 1, 2, 3)
NONCLAIM_FLAGS = (
    "time_order_recursion_claimed",
    "global_closure_proved",
    "O3_transition_P2_zero_claimed",
    "physical_layer_completed",
    "numeric_values_computed",
    "numeric_containment_claimed",
    "O3_transition_phase_equality_proved",
    "O3_transition_signed_gap_zero_claimed",
    "O3_transition_P2_pressure_gap_zero_claimed",
    "full_goal_complete",
)
PRESERVED_FALSE_FIELD_NAMES = (
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
    "old_errors_removed",
    "old_error_terms_removed",
    "global_five_moment_closure",
    "global_closure_proved",
    "O3_transition_P2_zero_claimed",
    "O3_transition_phase_equality_proved",
    "time_order_recursion_claimed",
    "numeric_values_computed",
    "numeric_containment_claimed",
    "physical_layer_completed",
)
FORBIDDEN_COUNT_TOKENS = (
    "numeric",
    "source_provider",
    "provider",
    "provider_calls",
    "checker",
    "checker_calls",
    "normalizer",
    "normalizer_calls",
    "reducer",
    "reducer_calls",
    "quadrature",
    "quadrature_calls",
    "recurrence",
    "recurrence_calls",
    "field",
    "field_calls",
    "field_evaluations",
    "integral",
    "integral_calls",
    "old_math",
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


def portable(value):
    if isinstance(value, dict):
        result = {}
        for key, child in value.items():
            key = key.replace("\\", "/") if isinstance(key, str) else key
            require(key not in result, "Portable metadata key collision: " + str(key))
            result[key] = portable(child)
        return result
    if isinstance(value, list):
        return [portable(child) for child in value]
    if isinstance(value, str) and ("\\" in value or (len(value) > 2 and value[1] == ":")):
        return value.replace("\\", "/")
    return value


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
        char in "0123456789abcdef" for char in value.lower()
    )


def valid_git_revision(value):
    return isinstance(value, str) and len(value) == 40 and all(
        char in "0123456789abcdef" for char in value.lower()
    )


def auth_map(value, label, authenticate=True):
    require(isinstance(value, dict), label + " is not a mapping")
    result = {}
    for name, digest in value.items():
        require(isinstance(name, str) and valid_sha(digest), label + " has a malformed SHA")
        if authenticate:
            read_bytes(resolve_path(name), digest)
        result[wire_path(name)] = digest.lower()
    return result


def flatten_hashes(value, label):
    result = {}

    def add(name, digest):
        require(isinstance(name, str) and valid_sha(digest), label + " has a malformed hash")
        name = wire_path(name)
        digest = digest.lower()
        require(name not in result or result[name] == digest, "Conflicting hash at " + name)
        result[name] = digest

    def visit(node, location):
        require(isinstance(node, dict), location + " is not a mapping")
        if not node:
            return
        if all(isinstance(k, str) and valid_sha(v) for k, v in node.items()):
            for key, value in node.items():
                add(key, value)
            return
        for key, child in node.items():
            if isinstance(child, dict):
                visit(child, location + "." + str(key))
            elif isinstance(child, list):
                for index, item in enumerate(child):
                    if isinstance(item, dict):
                        visit(item, location + "." + str(key) + "[" + str(index) + "]")

    visit(value, label)
    return result


def merge_hashes(*maps):
    merged = {}
    for mapping in maps:
        require(isinstance(mapping, dict), "Hash map is not a mapping")
        for name, digest in mapping.items():
            name = wire_path(name)
            require(valid_sha(digest), "Malformed hash map entry: " + name)
            digest = digest.lower()
            require(name not in merged or merged[name] == digest, "Conflicting hash at " + name)
            merged[name] = digest
    return merged


def bool_fields(*mappings):
    result = {}
    for mapping in mappings:
        if not isinstance(mapping, dict):
            continue
        for key, value in mapping.items():
            if isinstance(key, str) and isinstance(value, bool):
                result[key] = value
    return dict(sorted(result.items()))


def preserved_false_fields(mapping):
    return {
        key: mapping[key]
        for key in PRESERVED_FALSE_FIELD_NAMES
        if mapping.get(key) is False
    }


def bounded_parent_view(report):
    keys = (
        "schema", "component", "attempt_status", "passed", "source_only",
        "parent_commit", "parent_head", "parent_report_commit", "parent_report_parent_commit",
        "source_pin_count", "checkpoint_count", "preserved_failed_artifact_count",
        "actual_counts", "actual_projection_orders", "new_flags", "false_flags", "false_flag_values",
        "full_goal_complete", "source_packet_definition_accepted",
    )
    view = portable({key: report[key] for key in keys if key in report})
    view["parent_report_is_authenticated_view"] = True
    view["complete_parent_report_raw_sha256"] = PARENT_REPORT_SHA256
    view["omitted_fields_saved_in_complete_parent_report"] = sorted(set(report) - set(keys))
    require(len(encode(view)) < 100_000, "Bounded parent view exceeds 100 KB")
    return view


def result_reference(report, report_path):
    selector = report.get("full_result_artifact_selector")
    digest = report.get("full_result_raw_sha256")
    prefix = "$.lossless_artifacts[path="
    require(isinstance(selector, str) and selector.startswith(prefix) and selector.endswith("]"), "Binding result selector missing")
    result_name = selector[len(prefix) : -1]
    require(valid_sha(digest), "Binding result SHA is malformed")
    result_path = resolve_path(result_name)
    result_raw, result = read_json(result_path, digest)
    require(report.get("full_result") == result, "Binding full result is not exact")
    require(report.get("full_result_raw_bytes") == len(result_raw), "Binding result byte count changed")
    return {
        "path": wire_path(result_path),
        "sha256": digest.lower(),
        "bytes": len(result_raw),
        "schema": result.get("schema"),
        "component": result.get("component"),
    }, result


def load_parent(args):
    report_path = resolve_path(args.parent_report or wire_path(PARENT_REPORT_DEFAULT), WORK)
    report_sha = (args.parent_report_sha256 or PARENT_REPORT_SHA256).lower()
    require(report_sha == PARENT_REPORT_SHA256, "Parent binding report SHA is not authenticated")
    report_raw, report = read_json(report_path, report_sha)
    require(report.get("schema") == "original_compact_O3_transition_P2_source_binding_publication_v1", "Wrong binding report schema")
    require(report.get("component") == "compact_O3_transition_P2_source_binding", "Wrong binding report component")
    require(report.get("passed") is True and report.get("attempt_status") == "passed", "Binding publication did not pass")
    require(report.get("source_only") is True, "Binding publication is not source-only")
    require(report.get("parent_commit") == PARENT_BINDING_REPORT_PARENT_HEAD, "Binding parent commit changed")
    require(report.get("parent_head") == PARENT_BINDING_REPORT_PARENT_HEAD, "Binding parent head changed")
    require(report.get("parent_report_commit") == PARENT_BINDING_REPORT_COMMIT, "Binding report commit changed")
    require(report.get("parent_report_parent_commit") == PARENT_BINDING_REPORT_PARENT_COMMIT, "Binding report parent changed")
    require(report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT, "Binding source-pin count changed")
    require(report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Binding checkpoint count changed")
    require(report.get("preserved_failed_artifact_count") == PARENT_FAILED_COUNT, "Binding failed-artifact count changed")
    source_pins = auth_map(report.get("source_pins"), "Binding source pins")
    require(len(source_pins) == PARENT_SOURCE_PIN_COUNT, "Binding source pins are incomplete")
    inventory = flatten_hashes(report.get("parent_checkpoint_inventory"), "Binding checkpoint inventory")
    checkpoint_hashes = flatten_hashes({"checkpoint_hashes": report.get("checkpoint_hashes")}, "Binding checkpoint hashes")
    checkpoints = merge_hashes(inventory, checkpoint_hashes)
    require(len(checkpoints) == PARENT_CHECKPOINT_COUNT, "Binding checkpoints are incomplete")
    for name, digest in checkpoints.items():
        read_bytes(resolve_path(name), digest)
    protected = auth_map(report.get("protected_file_hashes"), "Binding protected files")
    require(len(protected) == PARENT_PROTECTED_COUNT, "Binding protected-file count changed")
    dirty = auth_map(report.get("legacy_dirty_path_hashes"), "Binding dirty files")
    require(len(dirty) == PARENT_DIRTY_COUNT, "Binding dirty-file count changed")
    failed = auth_map(report.get("preserved_failed_artifact_hashes_union"), "Binding failed artifacts")
    require(len(failed) == PARENT_FAILED_COUNT, "Binding failed-artifact map is incomplete")
    result_meta, result = result_reference(report, report_path)
    require(result_meta["sha256"] == PARENT_BINDING_RESULT_SHA256, "Binding result SHA changed")
    require(result.get("schema") == "original_compact_O3_transition_P2_source_binding", "Binding result schema changed")
    require(result.get("component") == "compact_O3_transition_P2_source_binding", "Binding result component changed")
    require(result.get("passed") is True and result.get("source_only") is True, "Binding result is not accepted source data")
    binding_flags = result.get("new_flags")
    require(isinstance(binding_flags, dict), "Binding new flags are missing")
    require(all(binding_flags.get(name) is True for name in (
        "O3_transition_P2_source_ports_bound",
        "O3_transition_P2_pressure_restoration_route_bound",
        "O3_transition_P2_ordinary_Z_source_projection_bound",
    )), "Binding source flags are not all true")
    legacy_raw, legacy = read_json(LEGACY_PARENT_REPORT_DEFAULT, LEGACY_PARENT_REPORT_SHA256)
    require(legacy.get("schema") == "original_compact_O3_transition_primitive_gap_publication_v1", "Legacy ledger schema changed")
    legacy_flags = legacy.get("false_flags")
    require(isinstance(legacy_flags, list) and len(legacy_flags) == PARENT_FALSE_FLAG_COUNT, "Inherited false-flag ledger is not exact 26")
    require(len(set(legacy_flags)) == PARENT_FALSE_FLAG_COUNT, "Inherited false-flag ledger has duplicates")
    legacy_values = legacy.get("false_flag_values")
    require(isinstance(legacy_values, dict), "Inherited false-flag values are missing")
    for key in legacy_flags:
        require(legacy_values.get(key) is False, "Inherited false flag changed: " + str(key))
    parent_bool = bool_fields(report, result, legacy)
    return {
        "report": report,
        "report_raw": report_raw,
        "report_meta": {
            "path": wire_path(report_path),
            "sha256": report_sha,
            "bytes": len(report_raw),
            "schema": report.get("schema"),
            "component": report.get("component"),
        },
        "report_view": bounded_parent_view(report),
        "result": result,
        "result_meta": result_meta,
        "source_pins": source_pins,
        "checkpoints": checkpoints,
        "protected": protected,
        "dirty": dirty,
        "failed": failed,
        "boolean_fields": parent_bool,
        "legacy_boolean_fields": bool_fields(legacy),
        "false_flags": list(legacy_flags),
        "false_flag_values": {key: False for key in legacy_flags},
        "legacy_report_meta": {
            "path": wire_path(LEGACY_PARENT_REPORT_DEFAULT),
            "sha256": LEGACY_PARENT_REPORT_SHA256,
            "bytes": len(legacy_raw),
            "schema": legacy.get("schema"),
            "component": legacy.get("component"),
        },
    }


def artifact_entry(path, digest=None, embed=False, outside_science_pins=False):
    path = Path(path)
    raw = read_bytes(path, digest)
    item = {
        "path": wire_path(path),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "raw_readonly": True,
        "raw_bytes_exact": True,
    }
    if outside_science_pins:
        item["outside_science_pins"] = True
    if embed:
        item["raw_bytes_base64"] = base64.b64encode(raw).decode("ascii")
    return item


def checkpoint_artifacts():
    result = []
    for path in CHECKPOINTS:
        if path.exists():
            raw = path.read_bytes()
            result.append({"path": wire_path(path), "sha256": sha_bytes(raw), "bytes": len(raw)})
    return result


def checkpoint_link_map(source_raw, result_raw):
    return {
        "source": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "entered": {"path": wire_path(ENTERED_CHECKPOINT), "sha256": sha_path(ENTERED_CHECKPOINT)},
        "result": {"path": wire_path(RESULT_CHECKPOINT), "sha256": sha_bytes(result_raw)},
        "diagnostics": {"path": wire_path(DIAGNOSTIC_CHECKPOINT), "sha256": sha_path(DIAGNOSTIC_CHECKPOINT)},
    }


def append_diagnostic(event):
    raw = encode(portable(event))
    DIAGNOSTIC_CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    with DIAGNOSTIC_CHECKPOINT.open("ab") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())


def load_source():
    raw, source = read_json(SOURCE_CHECKPOINT)
    require(source.get("schema") == INPUT_SCHEMA, "Source packet schema changed")
    require(source.get("component") == COMPONENT, "Source packet component changed")
    require(source.get("parent_commit") == PARENT_COMMIT, "Source packet parent commit changed")
    require(source.get("binding_result_raw_sha256") == PARENT_BINDING_RESULT_SHA256, "Source binding SHA changed")
    require(source.get("binding_source_checkpoint_sha256") == PARENT_BINDING_SOURCE_CHECKPOINT_SHA256, "Binding source checkpoint SHA changed")
    require(source.get("binding_report_raw_sha256") == PARENT_REPORT_SHA256, "Binding report SHA changed")
    require(isinstance(source.get("parent_false_flags"), list) and len(source["parent_false_flags"]) == PARENT_FALSE_FLAG_COUNT, "Parent false-flag ledger changed")
    require(source.get("parent_false_flag_values") == {key: False for key in source["parent_false_flags"]}, "Parent false-flag values changed")
    require(isinstance(source.get("binding_result"), dict), "Source binding result is missing")
    require(sha_bytes(encode(source["binding_result"])) == source["binding_result_json_sha256"], "Source binding bytes changed")
    require(source.get("source_union_count") == len(source.get("source_union_hashes", {})), "Source union count changed")
    for name, digest in source.get("source_union_hashes", {}).items():
        read_bytes(resolve_path(name), digest)
    return raw, source


def build_input_spec(parent, review_path=None, review_sha=None):
    review = None
    if review_path:
        require(review_sha and valid_sha(review_sha), "Review SHA is required when review path is supplied")
        review_file = resolve_path(review_path, WORK)
        review_raw = read_bytes(review_file, review_sha)
        review_data = json.loads(review_raw.decode("utf-8"))
        require(review_data.get("passed") is True, "Saved source-binding review did not pass")
        review = {
            "path": wire_path(review_file),
            "sha256": sha_bytes(review_raw),
            "bytes": len(review_raw),
            "schema": review_data.get("schema"),
            "component": review_data.get("component"),
            "passed": review_data.get("passed"),
            "outside_science_pins": True,
        }
    value = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "source_binding_report": parent["report_meta"],
        "source_binding_result": parent["result_meta"],
        "binding_result_raw_sha256": PARENT_BINDING_RESULT_SHA256,
        "parent_result_raw_sha256": PARENT_BINDING_RESULT_SHA256,
        "binding_source_checkpoint_sha256": PARENT_BINDING_SOURCE_CHECKPOINT_SHA256,
        "binding_report_raw_sha256": parent["report_meta"]["sha256"],
        "binding_actual_counts": clone(parent["result"].get("actual_counts")),
        "binding_projection_orders": clone(parent["result"].get("actual_projection_orders")),
        "binding_result": parent["result"],
        "binding_new_flags": parent["result"].get("new_flags"),
        "binding_actual_counts": parent["result"].get("actual_counts"),
        "binding_projection_orders": parent["result"].get("actual_projection_orders"),
        "review": review,
        "metadata_only_assembly": True,
        "science_reexecuted": False,
    }
    raw = write_new(INPUT_SPEC, value)
    return raw, value


def build_source(parent, input_raw, input_value, kernel_path, runner_path):
    kernel_raw = read_bytes(kernel_path)
    runner_raw = read_bytes(runner_path)
    additions = {
        wire_path(parent["report_meta"]["path"]): parent["report_meta"]["sha256"],
        wire_path(parent["legacy_report_meta"]["path"]): parent["legacy_report_meta"]["sha256"],
        wire_path(parent["result_meta"]["path"]): parent["result_meta"]["sha256"],
        wire_path(INPUT_SPEC): sha_bytes(input_raw),
        wire_path(kernel_path): sha_bytes(kernel_raw),
        wire_path(runner_path): sha_bytes(runner_raw),
    }
    binding_source_paths = [
        name for name, digest in parent["checkpoints"].items()
        if digest == PARENT_BINDING_SOURCE_CHECKPOINT_SHA256
    ]
    require(len(binding_source_paths) == 1, "Binding source checkpoint path is not uniquely pinned")
    additions[binding_source_paths[0]] = PARENT_BINDING_SOURCE_CHECKPOINT_SHA256
    source_union = merge_hashes(parent["source_pins"], additions)
    source = {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_report": parent["report_meta"],
        "parent_passed": True,
        "parent_result_data_only": True,
        "parent_report_view": parent["report_view"],
        "parent_report_boolean_fields": parent["boolean_fields"],
        "legacy_parent_report": parent["legacy_report_meta"],
        "legacy_parent_boolean_fields": parent["legacy_boolean_fields"],
        "preserved_false_fields": preserved_false_fields(parent["boolean_fields"]),
        "parent_false_flags": parent["false_flags"],
        "parent_false_flag_values": parent["false_flag_values"],
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "expected_failed_count_if_new_failure": FAILED_COUNT_IF_NEW_FAILURE,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "parent_source_pins": parent["source_pins"],
        "parent_checkpoint_hashes": parent["checkpoints"],
        "parent_protected_hashes": parent["protected"],
        "parent_dirty_hashes": parent["dirty"],
        "parent_failed_hashes": parent["failed"],
        "binding_result": clone(parent["result"]),
        "preserved_error_cap_metadata": clone(parent["report"]["preserved_error_cap_metadata"]),
        "retained_construction_data": clone(parent["report"]["retained_construction_data"]),
        "retained_integral_term_inventory": clone(parent["report"]["retained_integral_term_inventory"]),
        "retained_cap_provenance": clone(parent["report"]["retained_cap_provenance"]),
        "binding_result_raw_sha256": PARENT_BINDING_RESULT_SHA256,
        "binding_source_checkpoint_sha256": PARENT_BINDING_SOURCE_CHECKPOINT_SHA256,
        "binding_report_raw_sha256": parent["report_meta"]["sha256"],
        "binding_result_json_sha256": sha_bytes(encode(parent["result"])),
        "input_spec_metadata": {
            "path": wire_path(INPUT_SPEC),
            "sha256": sha_bytes(input_raw),
            "bytes": len(input_raw),
            "raw_readonly": True,
        },
        "input_raw_sha256": sha_bytes(input_raw),
        "kernel_path": wire_path(kernel_path),
        "kernel_sha256": sha_bytes(kernel_raw),
        "runner_path": wire_path(runner_path),
        "runner_sha256": sha_bytes(runner_raw),
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": list(NEW_FLAGS),
        "expected_projection_orders": list(EXPECTED_PROJECTION_ORDERS),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "raw_sha_pins": source_union,
        "source_packet_defined": True,
        "source_packet_definition_accepted": True,
        "source_packet_accepted_owned": True,
        "source_packet_accepted_inherited": False,
        "source_only": True,
        "metadata_only": False,
        "runtime_model_metadata": None,
        "math_reexecuted": False,
        "no_runner_execution": True,
    }
    return source


def prepare(args):
    require(not SOURCE_CHECKPOINT.exists(), "Source packet already exists; do not reenter")
    require(not INPUT_SPEC.exists(), "Input specification already exists; do not rewrite")
    require(all(not path.exists() for path in CHECKPOINTS), "A pressure-gap checkpoint already exists; do not reenter")
    parent = load_parent(args)
    kernel_path = resolve_path(args.kernel_path or wire_path(REPO / KERNEL_RELATIVE), REPO)
    runner_path = resolve_path(args.runner_path or wire_path(REPO / RUNNER_RELATIVE), REPO)
    require(kernel_path.resolve() == (REPO / KERNEL_RELATIVE).resolve(), "Kernel path is outside the frozen pressure-gap component")
    require(runner_path.resolve() == (REPO / RUNNER_RELATIVE).resolve(), "Runner path is outside the frozen pressure-gap component")
    require(kernel_path != Path(__file__).resolve(), "Metadata helper cannot be kernel")
    require(kernel_path.exists() and runner_path.exists(), "Kernel and runner must be supplied before prepare")
    if args.kernel_sha256:
        require(sha_path(kernel_path) == args.kernel_sha256.lower(), "Kernel SHA does not match supplied value")
    if args.runner_sha256:
        require(sha_path(runner_path) == args.runner_sha256.lower(), "Runner SHA does not match supplied value")
    input_raw, input_value = build_input_spec(parent, args.review, args.review_sha256)
    source = build_source(parent, input_raw, input_value, kernel_path, runner_path)
    source_raw = write_new(SOURCE_CHECKPOINT, source)
    diagnostic_raw = write_new(DIAGNOSTIC_CHECKPOINT, b"")
    brief = {
        "schema": COMPONENT + "_brief_v1",
        "component": COMPONENT,
        "phase": "prepared",
        "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_raw)},
        "diagnostic_checkpoint": {"path": wire_path(DIAGNOSTIC_CHECKPOINT), "sha256": sha_bytes(diagnostic_raw)},
        "parent_report": parent["report_meta"],
        "kernel": {"path": wire_path(kernel_path), "sha256": source["kernel_sha256"]},
        "runner": {"path": wire_path(runner_path), "sha256": source["runner_sha256"]},
        "metadata_only": True,
        "no_runner_execution": True,
        "do_not_reenter": True,
    }
    require(len(encode(brief)) < 100_000, "Stage brief exceeds 100 KB")
    brief_raw = write_new(BRIEF, brief)
    pending_parent = {
        "schema": COMPONENT + "_pending_parent_authority_v1",
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_report": parent["report_meta"],
        "parent_result": parent["result_meta"],
        "parent_passed": True,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "full_goal_complete": False,
        "do_not_reenter": True,
    }
    main_brief_meta = update_main_brief({
        "phase": "prepared_unentered",
        "component": COMPONENT,
        "source_checkpoint": {"path": wire_path(SOURCE_CHECKPOINT), "sha256": sha_bytes(source_raw)},
        "input_spec": {"path": wire_path(INPUT_SPEC), "sha256": sha_bytes(input_raw)},
        "diagnostic_checkpoint": {"path": wire_path(DIAGNOSTIC_CHECKPOINT), "sha256": sha_bytes(diagnostic_raw)},
        "code2": [
            {"path": wire_path(kernel_path), "sha256": source["kernel_sha256"], "bytes": kernel_path.stat().st_size},
            {"path": wire_path(runner_path), "sha256": source["runner_sha256"], "bytes": runner_path.stat().st_size},
        ],
        "source_union_count": source["source_union_count"],
        "parent_report": parent["report_meta"],
        "parent_result": parent["result_meta"],
        "pending_authoritative_parent": pending_parent,
        "no_math_imported": True,
        "no_math_called": True,
        "runner_executed": False,
        "metadata_only": True,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "do_not_reenter_after_entry": True,
        "next_executable_task": {
            "phase": "single_pressure_gap_derive_after_prepare",
            "source": "frozen_source_packet_only",
            "reentry": "forbidden",
        },
    })
    receipt = {
        "schema": COMPONENT + "_prepare_receipt_v1",
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "input_spec_sha256": sha_bytes(input_raw),
        "diagnostic_checkpoint_sha256": sha_bytes(diagnostic_raw),
        "source_union_count": source["source_union_count"],
        "source_union_hashes": source["source_union_hashes"],
        "metadata_only": True,
        "math_reexecuted": False,
        "runner_executed": False,
        "brief_sha256": sha_bytes(brief_raw),
        "main_resume_brief": main_brief_meta,
        "pending_authoritative_parent": pending_parent,
        "runtime_model_metadata": None,
    }
    receipt_raw = write_new(PREPARE_RECEIPT, receipt)
    return {"prepare": receipt, "prepare_sha256": sha_bytes(receipt_raw)}


def import_kernel(path):
    path = Path(path)
    require(path.resolve() != Path(__file__).resolve(), "Kernel path aliases metadata helper")
    spec = importlib.util.spec_from_file_location("_compact_o3_p2_pressure_gap_kernel", path)
    require(spec is not None and spec.loader is not None, "Cannot load pressure-gap kernel spec")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(callable(getattr(module, "derive", None)), "Pressure-gap kernel has no derive API")
    return module


def failure_result(error):
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "attempt_status": "failed",
        "derive_entry_calls": 1,
        "derive_calls": 0,
        "actual_counts": None,
        "new_flags": {},
        "protocol_violations": [],
        "failure_reason": "derive_exception",
        "failure_error_type": type(error).__name__,
        "failure_error": repr(error),
        "source_only": True,
        "metadata_only": False,
        "runtime_model_metadata": None,
    }


def run(args):
    source_raw, source = load_source()
    require(not ENTERED_CHECKPOINT.exists(), "Entered marker already exists; do not reenter")
    require(not RESULT_CHECKPOINT.exists(), "Result already exists; do not reenter")
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "kernel_path": source["kernel_path"],
        "kernel_sha256": source["kernel_sha256"],
        "derive_entry_calls": 1,
        "entry_attempts": 1,
        "kernel_imported": False,
        "runtime_model_metadata": None,
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
        "do_not_reenter": True,
    }
    entered_raw = write_new(ENTERED_CHECKPOINT, entered)
    append_diagnostic({
        "schema": DIAGNOSTIC_SCHEMA,
        "event": "entered_durable",
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_bytes(entered_raw),
    })
    result = None
    error = None
    kernel_imported = False
    derive_invoked = False
    try:
        module = import_kernel(resolve_path(source["kernel_path"], REPO))
        kernel_imported = True
        append_diagnostic({"schema": DIAGNOSTIC_SCHEMA, "event": "kernel_imported", "kernel_sha256": source["kernel_sha256"]})
        derive_invoked = True
        result = module.derive(clone(source), append_diagnostic)
        require(isinstance(result, dict), "derive returned no mapping")
    except Exception as exc:  # durable failure evidence is written below
        error = exc
        result = failure_result(exc)
        append_diagnostic({
            "schema": DIAGNOSTIC_SCHEMA,
            "event": "derive_failed",
            "error_type": type(exc).__name__,
            "error": repr(exc),
        })
    derive_calls = 1 if derive_invoked else 0
    # Keep the kernel object byte-for-byte under component_result while adding
    # manager-owned call accounting to the durable envelope.  The envelope is
    # what is validated and linked to the entered/journal checkpoints; it never
    # rewrites any kernel count or flag.
    kernel_result = clone(result)
    durable_result = clone(result)
    durable_result["component_result"] = kernel_result
    durable_result["execution_metadata"] = {
        "derive_entry_calls": 1,
        "derive_calls": derive_calls,
        "kernel_imported": kernel_imported,
        "new_math_entry": derive_calls,
        "old_math_reexecuted": False,
        "runner_executed": False,
    }
    durable_result["derive_entry_calls"] = 1
    durable_result["derive_calls"] = derive_calls
    result_raw = write_new(RESULT_CHECKPOINT, durable_result)
    validation = validate_result(durable_result, source, derive_entry_calls=1, derive_calls=derive_calls)
    receipt = {
        "schema": COMPONENT + "_run_receipt_v1",
        "component": COMPONENT,
        "entered_checkpoint_sha256": sha_bytes(entered_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "result_checkpoint_links": checkpoint_link_map(source_raw, result_raw),
        "diagnostic_checkpoint_sha256": sha_path(DIAGNOSTIC_CHECKPOINT),
        "derive_entry_calls": 1,
        "derive_calls": derive_calls,
        "kernel_imported": kernel_imported,
        "durable_result_saved_before_validation": True,
        "validation": validation,
        "protocol_violation_receipt": validation.get("protocol_violation_receipt"),
        "attempt_status": "passed" if validation["passed"] else "failed",
        "metadata_only": True,
        "component_result_preserved": True,
        "new_math_entry": derive_calls,
        "old_math_reexecuted": False,
        "math_reexecuted": False,
        "runner_executed": False,
        "runtime_model_metadata": None,
        "do_not_reenter": True,
    }
    receipt_raw = write_replace(RUN_RECEIPT, receipt)
    return {"run": receipt, "run_sha256": sha_bytes(receipt_raw)}


def recursive_forbidden_counts(node, location="result"):
    found = []
    if isinstance(node, dict):
        for key, value in node.items():
            key_text = str(key).lower()
            if isinstance(value, (int, float)) and not isinstance(value, bool) and value != 0:
                if any(token in key_text for token in FORBIDDEN_COUNT_TOKENS):
                    found.append({"path": location + "." + str(key), "value": value})
            if isinstance(value, bool) and value and any(token in key_text for token in FORBIDDEN_COUNT_TOKENS):
                found.append({"path": location + "." + str(key), "value": True})
            found.extend(recursive_forbidden_counts(value, location + "." + str(key)))
    elif isinstance(node, list):
        for index, value in enumerate(node):
            found.extend(recursive_forbidden_counts(value, location + "[" + str(index) + "]"))
    return found


def validate_result(result, source, derive_entry_calls=None, derive_calls=None):
    require(isinstance(result, dict), "Pressure-gap result is not a mapping")
    require(result.get("schema") == RESULT_SCHEMA, "Pressure-gap result schema changed")
    require(result.get("component") == COMPONENT, "Pressure-gap result component changed")
    require(isinstance(result.get("passed"), bool), "Pressure-gap result must contain explicit passed bool")
    violations = []
    raw_protocol = result.get("protocol_violations")
    if raw_protocol:
        require(isinstance(raw_protocol, list), "protocol_violations is not a list")
        violations.extend(clone(raw_protocol))
    violations.extend(recursive_forbidden_counts(result))
    for flag in NONCLAIM_FLAGS:
        if result.get(flag) is True:
            violations.append({"flag": flag, "reason": "nonclaim flag is true"})
    if result.get("attempt_status") is not None:
        expected_status = "passed" if result["passed"] else "failed"
        if result.get("attempt_status") != expected_status:
            violations.append({"field": "attempt_status", "reason": "status disagrees with explicit passed", "value": result.get("attempt_status")})
    if result.get("source_only") is not True:
        violations.append({"field": "source_only", "reason": "pressure-gap result must remain source-only"})
    if result.get("metadata_only") is not False:
        violations.append({"field": "metadata_only", "reason": "pressure-gap result must be non-metadata"})
    counts = result.get("actual_counts")
    if counts is not None:
        if not isinstance(counts, dict) or set(counts) != set(COUNT_KEYS):
            violations.append({"field": "actual_counts", "reason": "known count key set changed"})
        else:
            for key in COUNT_KEYS:
                if not isinstance(counts[key], int) or isinstance(counts[key], bool) or counts[key] < 0:
                    violations.append({"field": "actual_counts." + key, "reason": "malformed actual count", "value": counts[key]})
    elif result["passed"]:
        violations.append({"field": "actual_counts", "reason": "passed result has no durable actual counts"})
    flags = result.get("new_flags")
    if flags is not None:
        if not isinstance(flags, dict) or set(flags) != set(NEW_FLAGS):
            violations.append({"field": "new_flags", "reason": "pressure-gap new flag set changed"})
        elif not all(isinstance(flags[name], bool) for name in NEW_FLAGS):
            violations.append({"field": "new_flags", "reason": "pressure-gap new flags are not booleans"})
    elif result["passed"]:
        violations.append({"field": "new_flags", "reason": "passed result has no new flags"})
    if result.get("passed") is True:
        if counts != EXPECTED_COUNTS:
            violations.append({"field": "actual_counts", "reason": "passed result counts do not match expected contract", "value": clone(counts)})
        if flags is None or not isinstance(flags, dict) or not all(flags.get(name) is True for name in NEW_FLAGS):
            violations.append({"field": "new_flags", "reason": "passed result flags are incomplete"})
        if derive_entry_calls != 1 or derive_calls != 1:
            violations.append({"field": "derive_calls", "reason": "successful derive entry/call count is not one", "entry": derive_entry_calls, "calls": derive_calls})
    if result.get("derive_entry_calls") is not None and result.get("derive_entry_calls") != derive_entry_calls:
        violations.append({"field": "derive_entry_calls", "reason": "kernel result disagrees with manager call tracker"})
    if result.get("derive_calls") is not None and result.get("derive_calls") != derive_calls:
        violations.append({"field": "derive_calls", "reason": "kernel result disagrees with manager call tracker"})
    if isinstance(result.get("flags"), dict) and result.get("flags") != flags:
        violations.append({"field": "flags", "reason": "flags and new_flags differ"})
    if isinstance(result.get("audit_counts"), dict) and result.get("audit_counts") != counts:
        violations.append({"field": "audit_counts", "reason": "audit_counts and actual_counts differ"})
    if isinstance(counts, dict):
        for key in COUNT_KEYS:
            if key in result and result[key] != counts[key]:
                violations.append({"field": key, "reason": "top-level count disagrees with actual_counts"})
    binding_pins = result.get("binding_pins")
    if binding_pins is not None:
        expected_pins = {
            "binding_result_raw_sha256": PARENT_BINDING_RESULT_SHA256,
            "binding_source_checkpoint_sha256": PARENT_BINDING_SOURCE_CHECKPOINT_SHA256,
            "binding_report_raw_sha256": PARENT_REPORT_SHA256,
        }
        if binding_pins != expected_pins:
            violations.append({"field": "binding_pins", "reason": "binding input pins changed"})
    if result.get("projection_orders") is not None:
        if result["projection_orders"] != list(EXPECTED_PROJECTION_ORDERS):
            violations.append({"field": "projection_orders", "reason": "projection order contract changed"})
    if result.get("ordinary_Z_orders") is not None:
        if result["ordinary_Z_orders"] != list(EXPECTED_PROJECTION_ORDERS):
            violations.append({"field": "ordinary_Z_orders", "reason": "ordinary-Z order contract changed"})
    passed = result.get("passed") is True and not violations
    protocol_receipt = None
    if violations:
        protocol_receipt = {
            "schema": COMPONENT + "_protocol_violation_receipt_v1",
            "component": COMPONENT,
            "passed": False,
            "raw_result_preserved": True,
            "violations": clone(violations),
            "actual_counts": clone(counts),
            "forbidden_counts_preserved": True,
        }
    return {
        "schema": COMPONENT + "_result_validation_v1",
        "component": COMPONENT,
        "passed": passed,
        "explicit_passed": result["passed"],
        "attempt_status": "passed" if passed else "failed",
        "actual_counts": clone(counts),
        "new_flags": clone(flags),
        "protocol_violations": clone(violations),
        "protocol_violation_receipt": protocol_receipt,
        "raw_result_preserved": True,
        "numeric_source_old_math_counts_rewritten": False,
        "derive_entry_calls": derive_entry_calls,
        "derive_calls": derive_calls,
    }


def validate_scope_attributes():
    path = REPO / ".gitattributes"
    raw = read_bytes(path)
    lines = raw.decode("utf-8").splitlines()
    expected = scope_attribute_lines()
    require(all(lines.count(line) == 1 for line in expected), "Scope attribute lines are not exact once")
    return expected


def scope_attribute_lines():
    return [
        f"/{KERNEL_RELATIVE.as_posix()} -text",
        f"/{RUNNER_RELATIVE.as_posix()} -text",
        f"/{REPORT_RELATIVE.as_posix()} -text",
        f"/{MILESTONE_RELATIVE.as_posix()} -text",
    ]


def failed_artifacts(source_raw, result_raw, validation):
    if validation["passed"]:
        return {}
    result = {
        wire_path(SOURCE_CHECKPOINT): sha_bytes(source_raw),
        wire_path(ENTERED_CHECKPOINT): sha_path(ENTERED_CHECKPOINT),
        wire_path(RESULT_CHECKPOINT): sha_bytes(result_raw),
        wire_path(DIAGNOSTIC_CHECKPOINT): sha_path(DIAGNOSTIC_CHECKPOINT),
        wire_path(PREPARE_RECEIPT): sha_path(PREPARE_RECEIPT),
        wire_path(RUN_RECEIPT): sha_path(RUN_RECEIPT),
        wire_path(Path(__file__).resolve()): sha_path(Path(__file__).resolve()),
    }
    require(len(result) == 7, "New failure ledger must add exactly seven artifacts")
    return result


def next_task_metadata(passed):
    if passed:
        return {
            "phase": "independent_signed_full_W_sector_correlation",
            "source": "accepted_pressure_gap_result",
            "reentry": "forbidden",
            "retained_mismatches": ["geometry", "history", "pressure"],
            "claims": "no_global_closure_or_physical_completion",
        }
    return {
        "phase": "inspect_cached_pressure_gap_prefix",
        "source": "durable_saved_result_and_journal_only",
        "reentry": "forbidden",
        "claims": "no_global_closure_or_physical_completion",
    }


def report_for_publish(source, result_raw, result, validation, code2, input_artifact):
    source_raw = SOURCE_CHECKPOINT.read_bytes()
    current = checkpoint_artifacts()
    require(len(current) == 4, "Expected exactly four current checkpoint artifacts")
    current_hashes = {item["path"]: item["sha256"] for item in current}
    new_failed = failed_artifacts(source_raw, result_raw, validation)
    failed = merge_hashes(source["parent_failed_hashes"], new_failed)
    report = {
        "schema": PUBLICATION_SCHEMA,
        "component": COMPONENT,
        "attempt_status": validation["attempt_status"],
        "passed": validation["passed"],
        "explicit_result_passed": result.get("passed"),
        "metadata_only": True,
        "publication_metadata_only": True,
        "source_only": False,
        "parent_commit": PARENT_COMMIT,
        "parent_head": PARENT_COMMIT,
        "parent_report": source["parent_report"],
        "parent_report_sha256": source["parent_report"]["sha256"],
        "parent_passed": True,
        "parent_result_data_only": True,
        "legacy_parent_report": source["legacy_parent_report"],
        "parent_binding_result_raw_sha256": PARENT_BINDING_RESULT_SHA256,
        "parent_result_raw_sha256": PARENT_BINDING_RESULT_SHA256,
        "parent_source_pin_count": PARENT_SOURCE_PIN_COUNT,
        "parent_checkpoint_count": PARENT_CHECKPOINT_COUNT,
        "parent_failed_artifact_count": PARENT_FAILED_COUNT,
        "parent_protected_file_count": PARENT_PROTECTED_COUNT,
        "parent_dirty_path_count": PARENT_DIRTY_COUNT,
        "source_pin_count": source["source_union_count"],
        "source_union_count": source["source_union_count"],
        "source_pins": source["source_union_hashes"],
        "source_union_hashes": source["source_union_hashes"],
        "raw_sha_pins": source["raw_sha_pins"],
        "checkpoint_count": CURRENT_CHECKPOINT_COUNT,
        "checkpoint_paths": [item["path"] for item in current],
        "checkpoint_hashes": current_hashes,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "entered_checkpoint_sha256": sha_path(ENTERED_CHECKPOINT),
        "diagnostic_checkpoint_sha256": sha_path(DIAGNOSTIC_CHECKPOINT),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "result_checkpoint_links": checkpoint_link_map(source_raw, result_raw),
        "parent_checkpoint_hashes": source["parent_checkpoint_hashes"],
        "protected_file_hashes": source["parent_protected_hashes"],
        "legacy_dirty_path_hashes": source["parent_dirty_hashes"],
        "preserved_failed_artifact_hashes_union": failed,
        "preserved_failed_artifact_count": len(failed),
        "inherited_failed_artifacts_preserved": True,
        "parent_report_view": source["parent_report_view"],
        "parent_report_boolean_fields": source["parent_report_boolean_fields"],
        "legacy_parent_boolean_fields": source["legacy_parent_boolean_fields"],
        "preserved_false_fields": source["preserved_false_fields"],
        "parent_false_flags": source["parent_false_flags"],
        "parent_false_flag_values": source["parent_false_flag_values"],
        "expected_parent_false_flag_count": PARENT_FALSE_FLAG_COUNT,
        "binding_result": source["binding_result"],
        "preserved_error_cap_metadata": source["preserved_error_cap_metadata"],
        "retained_construction_data": source["retained_construction_data"],
        "retained_integral_term_inventory": source["retained_integral_term_inventory"],
        "retained_cap_provenance": source["retained_cap_provenance"],
        "binding_result_raw_sha256": PARENT_BINDING_RESULT_SHA256,
        "input_artifact": input_artifact,
        "input_raw_sha256": input_artifact["sha256"],
        "code2": code2,
        "current_code_hashes": {
            source["kernel_path"]: source["kernel_sha256"],
            source["runner_path"]: source["runner_sha256"],
        },
        "metadata_helper_artifact": artifact_entry(Path(__file__).resolve(), outside_science_pins=True),
        "external_draft": artifact_entry(Path(__file__).resolve(), outside_science_pins=True),
        "metadata_repair_not_code2": True,
        "full_result": clone(result),
        "component_result": clone(result.get("component_result", result)),
        "component_result_exact": True,
        "full_result_raw_bytes": len(result_raw),
        "full_result_raw_sha256": sha_bytes(result_raw),
        "raw_bytes_lossless": True,
        "result_checkpoint": {"path": wire_path(RESULT_CHECKPOINT), "sha256": sha_bytes(result_raw)},
        "actual_counts": clone(validation["actual_counts"]),
        "new_flags": clone(validation["new_flags"]),
        "protocol_violations": clone(validation["protocol_violations"]),
        "durable_result_validation": clone(validation),
        "execution_metadata": clone(result.get("execution_metadata")),
        "raw_result_preserved": True,
        "runtime_model_metadata": None,
        "new_math_entry": validation["derive_calls"],
        "old_math_reexecuted": False,
        "math_reexecuted": False,
        "runner_executed": False,
        "no_runner_execution": True,
        "full_goal_complete": False,
        "goal_active": True,
        "do_not_reenter": True,
        "next_executable_task": next_task_metadata(validation["passed"]),
        "exact_nine_scope": True,
        "scoped_paths": [path.as_posix() for path in SCOPED_RELATIVES],
        "scope_attribute_lines": scope_attribute_lines(),
        "input_spec_schema": INPUT_SCHEMA,
        "expected_counts": clone(EXPECTED_COUNTS),
        "expected_new_flags": list(NEW_FLAGS),
        "expected_projection_orders": list(EXPECTED_PROJECTION_ORDERS),
        "new_algebra_family_calls": (validation["actual_counts"] or {}).get("new_algebra_family_calls"),
        "pressure_identity_orders": (validation["actual_counts"] or {}).get("pressure_identity_orders"),
        "gap_function_rows": (validation["actual_counts"] or {}).get("gap_function_rows"),
        "side_quotient_premise_rows": (validation["actual_counts"] or {}).get("side_quotient_premise_rows"),
        "numeric_calls": (validation["actual_counts"] or {}).get("numeric_calls"),
        "source_provider_calls": (validation["actual_counts"] or {}).get("source_provider_calls"),
        "old_math_calls": (validation["actual_counts"] or {}).get("old_math_calls"),
    }
    return report


def update_main_brief(patch):
    require(MAIN_BRIEF.exists(), "Main resume brief is missing")
    old_raw, old = read_json(MAIN_BRIEF)
    require(len(old_raw) < 100_000, "Existing main resume brief exceeds 100 KB")
    preimage = WORK / f"{COMPONENT}_main_resume_brief_preimage_{sha_bytes(old_raw)[:12]}.json"
    if not preimage.exists():
        write_new(preimage, old_raw)
    result = clone(old)
    in_progress = result.get("in_progress")
    if not isinstance(in_progress, dict):
        in_progress = {}
    in_progress[COMPONENT] = portable(patch)
    result["in_progress"] = in_progress
    result["latest_component"] = COMPONENT
    result["latest_component_pointer"] = portable(patch)
    canonical = [
        wire_path(REPO / "AGENTS.md"),
        wire_path(REPO / MILESTONE_RELATIVE),
        wire_path(REPO / REPORT_RELATIVE),
    ]
    old_read_first = result.get("read_first")
    if not isinstance(old_read_first, list):
        old_read_first = []
    report_path = REPO / REPORT_RELATIVE
    milestone_path = REPO / MILESTONE_RELATIVE
    published = report_path.exists() and milestone_path.exists()
    if published:
        result["read_first"] = canonical + [item for item in old_read_first if item not in canonical]
        result["read_first_alias_repairs"] = {
            "AGENTS.md": canonical[0],
            "newdoc": canonical[1],
            "report": canonical[2],
        }
    else:
        # During prepare, retain the accepted parent's read_first list.  The
        # pending stage is exposed through in_progress/latest_component_pointer
        # and must not promote an unpublished report or milestone.
        result["read_first"] = clone(old_read_first)
    authority = patch.get("latest_authoritative_state")
    if authority is not None:
        require(isinstance(authority, dict), "Authoritative state is not a mapping")
        result["latest_authoritative_state"] = portable(authority)
        result["verified_local_and_remote_commit"] = authority.get("verified_local_and_remote_commit")
        result["latest_review"] = portable(patch.get("latest_review"))
        result["latest_qa"] = portable(patch.get("latest_qa"))
        result["next_remaining"] = portable(patch.get("next_remaining"))
    if published:
        result["latest_milestone"] = {
            "path": wire_path(milestone_path),
            "report": wire_path(report_path),
            "phase": patch.get("phase"),
            "report_sha256": patch.get("report_sha256"),
            "milestone_sha256": patch.get("milestone_sha256"),
        }
    raw = encode(result)
    if len(raw) > 100_000:
        compact_keys = {"continuation_detail_archive", "archived_detail", "archived_prior_component_details", "failed_source_attempt_evidence"}
        compact = {key: value for key, value in result.items() if key not in compact_keys}
        compact["archived_main_resume_brief_preimage"] = {
            "path": wire_path(preimage),
            "sha256": sha_bytes(old_raw),
            "bytes": len(old_raw),
            "preserved_exact_bytes": True,
        }
        compact["main_resume_brief_archival_compaction"] = True
        raw = encode(compact)
    require(len(raw) < 100_000, "Updated main resume brief exceeds 100 KB")
    write_replace(MAIN_BRIEF, raw)
    brief_meta = {
        "path": wire_path(MAIN_BRIEF),
        "sha256": sha_bytes(raw),
        "bytes": len(raw),
        "preimage": {"path": wire_path(preimage), "sha256": sha_bytes(old_raw), "bytes": len(old_raw)},
    }
    return brief_meta


def authoritative_state(report, report_raw, validation, head, review=None):
    return {
        "schema": COMPONENT + "_authoritative_state_v1",
        "component": COMPONENT,
        "commit": head,
        "head": head,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "result_sha256": report.get("full_result_raw_sha256"),
        "passed": report.get("passed") is True,
        "attempt_status": report.get("attempt_status"),
        "actual_counts": clone(report.get("actual_counts")),
        "new_flags": clone(report.get("new_flags")),
        "source_pin_count": report.get("source_pin_count"),
        "checkpoint_count": report.get("checkpoint_count"),
        "preserved_failed_artifact_count": report.get("preserved_failed_artifact_count"),
        "full_goal_complete": False,
        "read_first": [
            wire_path(REPO / "AGENTS.md"),
            wire_path(REPO / MILESTONE_RELATIVE),
            wire_path(REPO / REPORT_RELATIVE),
        ],
        "latest_review": clone(review) if review is not None else report.get("independent_review"),
        "latest_qa": clone(validation),
        "next_remaining": clone(report.get("next_executable_task")),
        "verified_local_and_remote_commit": head,
        "do_not_reenter": True,
    }


def publish(args):
    source_raw, source = load_source()
    require(ENTERED_CHECKPOINT.exists() and RESULT_CHECKPOINT.exists(), "Run checkpoints are missing")
    result_raw, result = read_json(RESULT_CHECKPOINT)
    require(RUN_RECEIPT.exists(), "Run receipt is missing")
    _, run_receipt = read_json(RUN_RECEIPT)
    validation = validate_result(
        result,
        source,
        derive_entry_calls=run_receipt.get("derive_entry_calls"),
        derive_calls=run_receipt.get("derive_calls"),
    )
    kernel = resolve_path(source["kernel_path"], REPO)
    runner = resolve_path(source["runner_path"], REPO)
    code2 = [artifact_entry(kernel, source["kernel_sha256"], embed=True), artifact_entry(runner, source["runner_sha256"], embed=True)]
    input_artifact = artifact_entry(resolve_path(source["input_spec_metadata"]["path"], WORK), source["input_spec_metadata"]["sha256"], embed=True)
    report = report_for_publish(source, result_raw, result, validation, code2, input_artifact)
    report_raw = write_replace(REPO / REPORT_RELATIVE, report)
    milestone = (
        "# Current original O3 transition P2 pressure gap (2026-10-27)\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}`; SHA-256 `{sha_bytes(report_raw)}`.\n\n"
        f"Attempt status: `{report['attempt_status']}`; explicit result passed: `{report['explicit_result_passed']}`; "
        f"actual counts: `{report['actual_counts']}`; new flags: `{report['new_flags']}`.\n\n"
        "This leaf consumes the accepted source-binding packet and records a formal P2 pressure-coefficient gap, "
        "not a restored physical-pressure (P0+p) gap, density-p2/P1 gap, or m15 result. The conditional source "
        "identities are En*G=dW-dE*Pc and G=Pn-Pc, with Pn=Wn/En and Pc=Wc/Ec, four ordinary-Z orders, "
        "eight side quotient premises, and En_0/Ec_0 nonzero domain assumptions. W is retained as an actual "
        "typed full-source producer and is not recovered from E times P2; same-Z pairing is formal only.\n\n"
        "Inherited 212 cumulative / 213 raw / 183 reduced terms, six signed channels, one reserved contribution, "
        "17 K terms, two outer B terms, and 16 J terms remain retained metadata, together with the 2B/16J "
        "cap/error provenance. The ordinary-Z recurrence is "
        "G_j=(dW_j-sum C(j,k)dE_k Pc_(j-k)-sum(k=1..j)C(j,k)En_k G_(j-k))/En_0; each row checks both "
        "the product identity and the recurrence-versus-independent-witness induction. No P2 value, pressure closure, global closure, "
        "temporal recursion, phase equality, numeric containment, signed-gap zero, or physical completion is claimed.\n\n"
        "Next concrete task: an independent signed full-W sector correlation using the retained actual geometry, "
        "history, and pressure mismatches; do not retry the failed primitive path. The saved entered attempt is "
        "immutable and must not be reentered.\n"
    ).encode("utf-8")
    milestone_raw = write_replace(REPO / MILESTONE_RELATIVE, milestone)
    handoff_prefix = (
        f"\n\n## {COMPONENT} publication\n\n"
        f"Report: `{REPORT_RELATIVE.as_posix()}` (SHA-256 `{sha_bytes(report_raw)}`).\n"
        f"Status: `{report['attempt_status']}`; actual counts `{report['actual_counts']}`; do not reenter.\n"
    ).encode("utf-8")
    preimages = []
    for relative in HANDOFF_RELATIVES:
        path = REPO / relative
        old = path.read_bytes() if path.exists() else b""
        preimage = WORK / f"{COMPONENT}_handoff_preimage_{path.stem}_{sha_bytes(old)[:12]}.md"
        if not preimage.exists():
            write_new(preimage, old)
        preimages.append({"path": wire_path(preimage), "sha256": sha_path(preimage), "bytes": len(old)})
        write_replace(path, handoff_prefix + old)
    attrs_path = REPO / ".gitattributes"
    attrs_old = attrs_path.read_bytes() if attrs_path.exists() else b""
    attrs_text = attrs_old.decode("utf-8")
    for line in (
        f"/{KERNEL_RELATIVE.as_posix()} -text",
        f"/{RUNNER_RELATIVE.as_posix()} -text",
        f"/{REPORT_RELATIVE.as_posix()} -text",
        f"/{MILESTONE_RELATIVE.as_posix()} -text",
    ):
        if line not in attrs_text.splitlines():
            attrs_text += ("\n" if attrs_text and not attrs_text.endswith("\n") else "") + line + "\n"
    write_replace(attrs_path, attrs_text.encode("utf-8"))
    brief = update_main_brief({
        "phase": "published_passed" if report["passed"] else "published_failed",
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "milestone": wire_path(REPO / MILESTONE_RELATIVE),
        "milestone_sha256": sha_bytes(milestone_raw),
        "attempt_status": report["attempt_status"],
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "next_executable_task": report["next_executable_task"],
    })
    receipt = {
        "schema": COMPONENT + "_publication_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "attempt_status": report["attempt_status"],
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "handoff_preimages": preimages,
        "exact_nine_scope": True,
        "runtime_model_metadata": None,
        "main_resume_brief": brief,
    }
    receipt_raw = write_replace(PUBLISH_RECEIPT, receipt)
    return {"publication": receipt, "publication_sha256": sha_bytes(receipt_raw)}


def git_output(*args):
    return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()


def verify_published_commit(report):
    head = git_output("rev-parse", "HEAD")
    require(valid_git_revision(head), "Local Git head is not a 40-hex revision")
    require(git_output("rev-parse", "HEAD^").lower() == PARENT_COMMIT.lower(), "Published head parent changed")
    branch = git_output("branch", "--show-current")
    require(branch, "Published branch is unnamed")
    remote_parts = git_output("ls-remote", "origin", "refs/heads/" + branch).split()
    require(remote_parts and valid_git_revision(remote_parts[0]), "Remote branch head is not a 40-hex revision")
    require(remote_parts[0].lower() == head.lower(), "Published head is not synchronized with origin")
    changed = set(git_output("diff", "--name-only", PARENT_COMMIT + "..HEAD").splitlines())
    expected = {path.as_posix() for path in SCOPED_RELATIVES}
    require(changed == expected, "Published committed scope is not exactly nine paths")
    validate_scope_attributes()
    for relative in SCOPED_RELATIVES:
        local_path = REPO / relative
        require(local_path.exists(), "Published scoped path is missing: " + relative.as_posix())
        blob_size = int(git_output("cat-file", "-s", "HEAD:" + relative.as_posix()))
        require(blob_size == local_path.stat().st_size, "Published blob byte count changed: " + relative.as_posix())
        committed_blob = git_output("rev-parse", "HEAD:" + relative.as_posix()).lower()
        worktree_blob = git_output("hash-object", wire_path(local_path)).lower()
        require(valid_git_revision(committed_blob) and committed_blob == worktree_blob, "Published bytes differ from worktree: " + relative.as_posix())
    require(report.get("parent_commit") == PARENT_COMMIT, "Published parent commit metadata changed")
    require(report.get("parent_head") == PARENT_COMMIT, "Published parent head metadata changed")
    return head


def resume(args):
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "Publication report is missing")
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT, report.get("full_result_raw_sha256"))
    require(report.get("full_result") == result, "Published result changed")
    stored_validation = report.get("durable_result_validation")
    require(isinstance(stored_validation, dict), "Published validation metadata is missing")
    validation = validate_result(
        result,
        source,
        derive_entry_calls=stored_validation.get("derive_entry_calls"),
        derive_calls=stored_validation.get("derive_calls"),
    )
    require(report.get("durable_result_validation") == validation, "Published validation changed")
    require(report.get("component_result") == result.get("component_result", result), "Published kernel result changed")
    require(report.get("full_result_raw_bytes") == len(result_raw), "Published result byte count changed")
    require(report.get("source_checkpoint_sha256") in (None, sha_bytes(source_raw)), "Published source link changed")
    require(report.get("result_checkpoint_links") == checkpoint_link_map(source_raw, result_raw), "Published checkpoint links changed")
    head = verify_published_commit(report)
    authority = authoritative_state(report, report_raw, validation, head)
    brief = update_main_brief({
        "phase": "resumed_published",
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "head": head,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "next_executable_task": report["next_executable_task"],
        "latest_authoritative_state": authority,
        "read_first": authority["read_first"],
        "latest_review": authority["latest_review"],
        "latest_qa": authority["latest_qa"],
        "next_remaining": authority["next_remaining"],
        "verified_local_and_remote_commit": head,
    })
    receipt = {
        "schema": COMPONENT + "_resume_receipt_v1",
        "component": COMPONENT,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "head": head,
        "do_not_reenter": True,
        "main_resume_brief": brief,
    }
    receipt_raw = write_replace(RESUME_RECEIPT, receipt)
    return {"resume": receipt, "resume_sha256": sha_bytes(receipt_raw)}


def seal(args):
    require(args.review and args.review_sha256, "Seal requires an independent review path and SHA")
    require(valid_sha(args.review_sha256), "Independent review SHA is malformed")
    review_path = resolve_path(args.review, WORK)
    review_raw = read_bytes(review_path, args.review_sha256.lower())
    review = json.loads(review_raw.decode("utf-8"))
    require(review.get("passed") is True, "Independent review did not pass")
    report_raw, report = read_json(REPO / REPORT_RELATIVE)
    require(report.get("schema") == PUBLICATION_SCHEMA and report.get("component") == COMPONENT, "Publication report is missing")
    source_raw, source = load_source()
    result_raw, result = read_json(RESULT_CHECKPOINT, report.get("full_result_raw_sha256"))
    stored_validation = report.get("durable_result_validation")
    require(isinstance(stored_validation, dict), "Published validation metadata is missing before seal")
    validation = validate_result(
        result,
        source,
        derive_entry_calls=stored_validation.get("derive_entry_calls"),
        derive_calls=stored_validation.get("derive_calls"),
    )
    require(report.get("full_result") == result, "Published result changed before seal")
    require(report.get("component_result") == result.get("component_result", result), "Published kernel result changed before seal")
    require(report.get("result_checkpoint_links") == checkpoint_link_map(source_raw, result_raw), "Published checkpoint links changed before seal")
    require(report.get("durable_result_validation") == validation, "Published validation changed before seal")
    head = verify_published_commit(report)
    review_meta = {"path": wire_path(review_path), "sha256": sha_bytes(review_raw), "bytes": len(review_raw), "raw_readonly": True, "outside_science_pins": True}
    authority = authoritative_state(report, report_raw, validation, head, review_meta)
    brief = update_main_brief({
        "phase": "sealed",
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "seal": wire_path(SEAL_RECEIPT),
        "head": head,
        "goal_active": True,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "next_executable_task": report["next_executable_task"],
        "latest_authoritative_state": authority,
        "read_first": authority["read_first"],
        "latest_review": review_meta,
        "latest_qa": authority["latest_qa"],
        "next_remaining": authority["next_remaining"],
        "verified_local_and_remote_commit": head,
    })
    receipt = {
        "schema": SEAL_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "report": wire_path(REPO / REPORT_RELATIVE),
        "report_sha256": sha_bytes(report_raw),
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "result_checkpoint_sha256": sha_bytes(result_raw),
        "durable_result_validation": validation,
        "independent_review": review_meta,
        "head": head,
        "full_goal_complete": False,
        "do_not_reenter": True,
        "main_resume_brief": brief,
    }
    receipt_raw = write_replace(SEAL_RECEIPT, receipt)
    return {"seal": receipt, "seal_sha256": sha_bytes(receipt_raw)}


def parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "publish", "resume", "seal"))
    parser.add_argument("--parent-report")
    parser.add_argument("--parent-report-sha256")
    parser.add_argument("--kernel-path")
    parser.add_argument("--kernel-sha256")
    parser.add_argument("--runner-path")
    parser.add_argument("--runner-sha256")
    parser.add_argument("--review")
    parser.add_argument("--review-sha256")
    return parser


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


