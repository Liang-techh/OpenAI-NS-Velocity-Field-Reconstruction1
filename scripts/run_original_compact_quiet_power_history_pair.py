"""Unexecuted metadata runner for the compact quiet-power history-pair leaf.

This file freezes the authenticated parent, input specification, new-kernel
paths, and durable checkpoint protocol.  A later explicit ``run`` action may
import only ``derive(source, journal_callback)`` from the new kernel once.
This scaffold does not import ancestor mathematics, providers, proofs,
checkers, normalizers, quadratures, or recurrence code, and it is intentionally
unexecuted in the preparation turn.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import traceback


ROOT = Path(
    "C:/Users/z5242/Documents/Codex/2026-09-13/"
    "before-doing-any-project-work-verify/work/st073-transition-next"
)
BRIEF = ROOT.parents[1] / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "compact_quiet_power_history_pair"
RUNNER = "scripts/run_original_compact_quiet_power_history_pair.py"
KERNEL = "src/openai_ns_reconstruction/original_compact_quiet_power_history_pair.py"
PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_signed_history_defect_transport.json"
)
PARENT_REPORT_SHA256 = (
    "dc367908bd44991ef6f4a973d688b4bd27bb0cab43464e0e396b3abb1201c861"
)
PARENT_RESULT_SHA256 = (
    "b91da155898c38a9ca4b26595983a38aa04b2965ae3db0a908fadcea97b34241"
)
# The next repository commit is expected to be directly based on the current
# signed-history commit.  The parent report itself records 870f9548 as its
# parent commit.
PARENT_COMMIT = "edf2b3fed5f44fb255fe424b2dddd393417087a4"
PARENT_REPORT_PARENT_COMMIT = "870f95489b93c91f8f3dabdba0ddf283ffe6c927"

INPUT_SCHEMA = "original_compact_quiet_power_history_pair_inputs_v1"
INPUT_SPEC_SCHEMA = "compact_quiet_power_history_pair_input_v1"
ENTERED_SCHEMA = "original_compact_quiet_power_history_pair_entered_v1"
RESULT_SCHEMA = "original_compact_quiet_power_history_pair"
DIAGNOSTIC_SCHEMA = "compact_quiet_power_history_pair_stage_diagnostic_v1"
INHERITED_FALSE_FLAG_COUNT = 26
PARENT_SOURCE_PIN_COUNT = 2153
PARENT_INHERITED_CHECKPOINT_COUNT = 224
PARENT_CHECKPOINT_COUNT = 228
PARENT_FAILED_COUNT = 162
PARENT_PROTECTED_COUNT = 16
PARENT_DIRTY_COUNT = 2

NEW_FLAGS = (
    "quiet_power_native_density_zero_source_proved",
    "quiet_power_C3_density_zero_source_proved",
    "quiet_power_common_history_kernel_source_bound",
    "quiet_power_signed_history_gap_transport_proved",
)

# These are source-bound scope descriptors, not numerical values or proofs.
# They make the intended leaf contract explicit without importing an ancestor
# source consumer.
QUIET_POWER_SOURCE_SCOPE = {
    "field": "quiet_O3_power",
    "ordinary_Z_orders": [0, 1, 2, 3],
    "coordinate_interval": ["0", "2"],
    "ordinary_Z_domain": ["-1", "1"],
    "ordinary_Z_function_identity_not_samples": True,
    "zero_local_densities_only": True,
    "zero_density_whole_interval": True,
    "zero_density_provenance": "source_bound_flat_branch_and_admission_certificate",
    "existing_C3_partial_history_claimed": False,
}
QUIET_POWER_KERNEL_SCOPE = {
    "same_relative_log_offset": True,
    "same_physical_measure": True,
    "absolute_Rm_Rc_cross_archive_identity_claimed": False,
    "source_integrand_kernel": "exp(-lambda*(2-s))",
    "memory_kernel": "exp(-2*lambda)",
}
QUIET_POWER_TRANSPORT_SCOPE = {
    "formula": "deltaT_out=exp(-2*lambda)*deltaT_in",
    "endpoint_coordinate": "2",
    "ordinary_Z_orders": [0, 1, 2, 3],
    "unknown_entry_gap_retained": True,
    "predecessor_inlet_retained": True,
    "predecessor_inlet_zeroed": False,
}

CHECKPOINT_DIR = ROOT / "tmp/original_compact_quiet_power_history_pair"
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

# These fields are owned by the current attempt.  In particular, parent
# result booleans named passed, source_packet_defined, kernel_imported, or
# metadata_only must never be copied over current fallback state.
CURRENT_RESULT_FIELDS = frozenset(
    {
        "schema",
        "component",
        "passed",
        "metadata_only",
        "source_only",
        "input_digest",
        "runtime_model_metadata",
        "original_source_evaluations",
        "source_packet_defined",
        "kernel_imported",
        "derive_entry_calls",
        "derive_calls",
        "entry_attempts",
    }
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def sha_path(path):
    return sha_bytes(Path(path).read_bytes())


def _encode(value):
    return (
        json.dumps(
            value,
            ensure_ascii=True,
            sort_keys=True,
            allow_nan=False,
            separators=(",", ":"),
        ).encode("utf-8")
        + b"\n"
    )


def _clone(value):
    return json.loads(_encode(value).decode("utf-8"))


def _digest(value):
    return sha_bytes(
        json.dumps(
            value,
            ensure_ascii=True,
            sort_keys=True,
            allow_nan=False,
            separators=(",", ":"),
        ).encode("utf-8")
    )


def _write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = value if isinstance(value, bytes) else _encode(value)
    with path.open("xb") as stream:
        stream.write(raw)
        stream.flush()
        os.fsync(stream.fileno())


def _read_json(path):
    raw = Path(path).read_bytes()
    return raw, json.loads(raw.decode("utf-8"))


def _resolve(name):
    path = Path(name)
    if path.is_absolute():
        return path
    direct = ROOT / path
    if direct.exists():
        return direct
    return ROOT / "experiments/root_st073" / path


def _raw_auth(name, expected):
    require(
        isinstance(expected, str)
        and len(expected) == 64
        and all(character in "0123456789abcdef" for character in expected),
        "Malformed SHA-256 for " + str(name),
    )
    path = _resolve(name)
    require(path.exists(), "Missing authenticated artifact: " + str(name))
    raw = path.read_bytes()
    require(sha_bytes(raw) == expected, "Changed authenticated artifact: " + str(name))
    return raw


def _is_sha(value):
    return isinstance(value, str) and len(value) == 64 and all(
        character in "0123456789abcdef" for character in value
    )


def _hash_map(value, label):
    require(isinstance(value, dict), label + " is not a mapping")
    result = {}
    for name, expected in value.items():
        require(isinstance(name, str) and _is_sha(expected), label + " has a malformed SHA")
        _raw_auth(name, expected)
        result[name] = expected
    return result


def _flatten_hashes(value, label):
    """Flatten nested raw-hash maps without decoding referenced files."""

    result = {}

    def add(name, expected):
        require(isinstance(name, str) and _is_sha(expected), label + " has a malformed SHA")
        require(name not in result or result[name] == expected, "Conflicting hash at " + name)
        result[name] = expected

    def visit(node, path):
        require(isinstance(node, dict), path + " is not a mapping")
        if not node:
            return
        if all(isinstance(key, str) and _is_sha(child) for key, child in node.items()):
            for key, child in node.items():
                add(key, child)
            return
        for key, child in node.items():
            child_path = path + "." + str(key)
            if isinstance(child, dict):
                visit(child, child_path)
            elif isinstance(child, list):
                for index, item in enumerate(child):
                    if isinstance(item, dict):
                        visit(item, child_path + "[" + str(index) + "]")

    visit(value, label)
    return result


def _merge_hashes(*maps):
    merged = {}
    for mapping in maps:
        require(isinstance(mapping, dict), "Hash map is not a mapping")
        for name, expected in mapping.items():
            require(isinstance(name, str) and _is_sha(expected), "Malformed hash map entry")
            require(name not in merged or merged[name] == expected, "Conflicting hash at " + name)
            merged[name] = expected
    return merged


def _parent_result_reference(report):
    selector = report.get("full_result_artifact_selector")
    prefix = "$.lossless_artifacts[path="
    require(
        isinstance(selector, str) and selector.startswith(prefix) and selector.endswith("]"),
        "Parent result selector is malformed",
    )
    path = selector[len(prefix) : -1]
    expected = report.get("full_result_raw_sha256")
    require(path and _is_sha(expected), "Parent result selector or SHA is malformed")
    return path, expected


def _parent_report_view(report):
    excluded = [
        "lossless_artifacts",
        "full_result",
        "component_result",
        "full_result_raw_bytes",
        "full_result_raw_sha256",
    ]
    view = {key: _clone(value) for key, value in report.items() if key not in excluded}
    view["parent_report_is_authenticated_view"] = True
    view["parent_report_view_excluded_fields"] = list(excluded)
    return view


def _inherited_true_flags(report, result, false_flags):
    inherited = {}
    for mapping in (report, result):
        for key, value in mapping.items():
            if isinstance(value, bool) and value and key not in false_flags:
                inherited[key] = True
    return inherited


def _parent_scope():
    report_raw = _raw_auth(PARENT_REPORT, PARENT_REPORT_SHA256)
    report = json.loads(report_raw.decode("utf-8"))
    require(
        report.get("schema")
        == "original_compact_signed_history_defect_transport_publication_v1"
        and report.get("component") == "compact_signed_history_defect_transport"
        and report.get("parent_commit") == PARENT_REPORT_PARENT_COMMIT
        and report.get("attempt_status") == "passed"
        and report.get("runtime_model_metadata") is None
        and report.get("source_pin_count") == PARENT_SOURCE_PIN_COUNT
        and report.get("checkpoint_count") == PARENT_CHECKPOINT_COUNT
        and report.get("inherited_checkpoint_count") == PARENT_INHERITED_CHECKPOINT_COUNT,
        "Signed-history parent report scope changed",
    )
    result_path, result_sha = _parent_result_reference(report)
    require(result_sha == PARENT_RESULT_SHA256, "Signed-history parent result SHA changed")
    result_raw = _raw_auth(result_path, result_sha)
    result = json.loads(result_raw.decode("utf-8"))
    require(
        report.get("full_result_raw_bytes") == len(result_raw)
        and report.get("full_result") == result
        and report.get("component_result") == result,
        "Signed-history parent result is not exact",
    )
    require(
        result.get("schema") == "original_compact_signed_history_defect_transport"
        and result.get("component") == "compact_signed_history_defect_transport"
        and result.get("passed") is True
        and result.get("metadata_only") is False
        and result.get("source_only") is True
        and result.get("runtime_model_metadata") is None
        and result.get("original_source_evaluations") == 0,
        "Signed-history parent result scope changed",
    )
    for key in (
        "actual_native_C3_amplitude_identity_proved",
        "signed_history_source_gap_functions_defined",
        "native_N_scaling_applied_once_bound",
        "five_signed_density_defect_identities_proved",
        "paired_history_defect_transport_proved",
    ):
        require(report.get(key) is True and result.get(key) is True, "Parent flag lost: " + key)
    false_flags = report.get("false_flags")
    require(
        isinstance(false_flags, list)
        and len(false_flags) == INHERITED_FALSE_FLAG_COUNT
        and len(set(false_flags)) == INHERITED_FALSE_FLAG_COUNT,
        "Signed-history parent false-flag inventory changed",
    )
    for name in false_flags:
        require(report.get(name) is False and result.get(name) is False, "Parent false flag changed: " + name)
    source_pins = _hash_map(report.get("source_pins"), "Signed-history parent source pins")
    require(len(source_pins) == PARENT_SOURCE_PIN_COUNT, "Signed-history parent source-pin count changed")
    inventory = _flatten_hashes(report.get("parent_checkpoint_inventory"), "Parent checkpoint inventory")
    own = _flatten_hashes({"checkpoint_hashes": report.get("checkpoint_hashes")}, "Parent checkpoint hashes")
    complete = _merge_hashes(inventory, own)
    require(
        len(inventory) == PARENT_INHERITED_CHECKPOINT_COUNT
        and len(own) == 4
        and len(complete) == PARENT_CHECKPOINT_COUNT,
        "Signed-history parent checkpoint maps changed",
    )
    for name, expected in complete.items():
        _raw_auth(name, expected)
    protected = _hash_map(report.get("protected_file_hashes"), "Parent protected files")
    dirty = _hash_map(report.get("legacy_dirty_path_hashes"), "Parent dirty files")
    failed = _hash_map(report.get("preserved_failed_artifact_hashes_union"), "Parent failed artifacts")
    require(
        len(protected) == PARENT_PROTECTED_COUNT
        and len(dirty) == PARENT_DIRTY_COUNT
        and len(failed) == PARENT_FAILED_COUNT,
        "Signed-history parent protected, dirty, or failed maps changed",
    )
    return {
        "report_raw": report_raw,
        "report": report,
        "report_view": _parent_report_view(report),
        "result_raw": result_raw,
        "result": result,
        "result_path": result_path,
        "result_sha256": result_sha,
        "source_pins": source_pins,
        "inventory": inventory,
        "own": own,
        "complete": complete,
        "protected": protected,
        "dirty": dirty,
        "failed": failed,
        "false_flags": list(false_flags),
        "inherited_true": _inherited_true_flags(report, result, false_flags),
    }


def _code_hashes():
    result = {}
    for name in (RUNNER, KERNEL):
        path = _resolve(name)
        require(path.exists(), "New quiet-power code blob is missing: " + name)
        result[name] = sha_path(path)
    return result


def _scope_bundle(spec):
    source_scope = spec.get("quiet_power_source_scope", QUIET_POWER_SOURCE_SCOPE)
    kernel_scope = spec.get("quiet_power_common_kernel_scope", QUIET_POWER_KERNEL_SCOPE)
    transport_scope = spec.get("quiet_power_transport_scope", QUIET_POWER_TRANSPORT_SCOPE)
    require(source_scope == QUIET_POWER_SOURCE_SCOPE, "Quiet-power source scope changed")
    require(kernel_scope == QUIET_POWER_KERNEL_SCOPE, "Quiet-power kernel scope changed")
    require(transport_scope == QUIET_POWER_TRANSPORT_SCOPE, "Quiet-power transport scope changed")
    return {
        "source": _clone(source_scope),
        "common_kernel": _clone(kernel_scope),
        "history_transport": _clone(transport_scope),
    }


def _input_spec(path):
    path = Path(path)
    require(path.is_absolute(), "--input-spec must be an absolute path")
    require(path.exists(), "Quiet-power input spec is missing: " + str(path))
    raw = path.read_bytes()
    value = json.loads(raw.decode("utf-8"))
    require(
        isinstance(value, dict)
        and value.get("schema") == INPUT_SPEC_SCHEMA
        and value.get("component") == COMPONENT,
        "Quiet-power input-spec schema changed",
    )
    require(isinstance(value.get("additional_raw_pins", {}), dict), "Input additional_raw_pins is not a map")
    require(isinstance(value.get("static_receipts", []), list), "Input static_receipts is not a list")
    scope = _scope_bundle(value)
    flat = value.get("flat_admission_metadata")
    if flat is not None:
        require(isinstance(flat, dict), "flat_admission_metadata is not a mapping")
        require(flat.get("zero_local_densities_only") is True, "Flat admission metadata overclaims density scope")
        provenance = flat.get("zero_density_provenance", QUIET_POWER_SOURCE_SCOPE["zero_density_provenance"])
        require(provenance == QUIET_POWER_SOURCE_SCOPE["zero_density_provenance"], "Flat admission density provenance changed")
    return path, raw, sha_bytes(raw), value, scope


def _input_pins(spec):
    additional = _hash_map(spec.get("additional_raw_pins", {}), "Input additional raw pins")
    receipts = {}
    for row in spec.get("static_receipts", []):
        require(isinstance(row, dict) and isinstance(row.get("path"), str), "Malformed input static receipt")
        expected = row.get("sha256")
        require(_is_sha(expected), "Malformed input static receipt SHA")
        _raw_auth(row["path"], expected)
        receipts[row["path"]] = expected
    return _merge_hashes(additional, receipts)


def _build_packet(input_spec_path):
    parent = _parent_scope()
    path, raw, spec_sha, spec, scope = _input_spec(input_spec_path)
    code = _code_hashes()
    extra = _input_pins(spec)
    stage = _merge_hashes(
        {
            PARENT_REPORT: PARENT_REPORT_SHA256,
            parent["result_path"]: parent["result_sha256"],
            RUNNER: code[RUNNER],
            KERNEL: code[KERNEL],
            str(path): spec_sha,
        },
        extra,
    )
    source_union = _merge_hashes(parent["source_pins"], stage)
    for name, expected in source_union.items():
        _raw_auth(name, expected)
    digest = _digest(
        {
            "parent_report_sha256": PARENT_REPORT_SHA256,
            "parent_result_sha256": parent["result_sha256"],
            "current_code_hashes": code,
            "input_spec_sha256": spec_sha,
            "stage_hashes": stage,
        }
    )
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "parent_report": parent["report_view"],
        "parent_report_is_authenticated_view": True,
        "parent_report_view_excluded_fields": parent["report_view"][
            "parent_report_view_excluded_fields"
        ],
        "parent_result": parent["result"],
        "parent_report_metadata": {
            "path": PARENT_REPORT,
            "raw_sha256": PARENT_REPORT_SHA256,
            "raw_bytes": len(parent["report_raw"]),
            "schema": parent["report"].get("schema"),
            "component": parent["report"].get("component"),
        },
        "parent_result_metadata": {
            "path": parent["result_path"],
            "raw_sha256": parent["result_sha256"],
            "raw_bytes": len(parent["result_raw"]),
            "schema": parent["result"].get("schema"),
            "component": parent["result"].get("component"),
        },
        "input": spec,
        "input_spec": {"path": str(path), "raw_sha256": spec_sha, "raw_bytes": len(raw)},
        "quiet_power_scope": scope,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_sha256": parent["result_sha256"],
        "parent_result_path": parent["result_path"],
        "parent_commit": PARENT_COMMIT,
        "input_digest": digest,
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "stage_hashes": stage,
        "current_code_hashes": dict(code),
        "implementation_hashes": dict(code),
        "parent_checkpoint_inventory": dict(parent["inventory"]),
        "checkpoint_hashes": dict(parent["own"]),
        "complete_parent_checkpoint_hashes": dict(parent["complete"]),
        "complete_parent_checkpoint_count": len(parent["complete"]),
        "preserved_source_pin_count": len(parent["source_pins"]),
        "protected_file_hashes": dict(parent["protected"]),
        "legacy_dirty_path_hashes": dict(parent["dirty"]),
        "preserved_failed_artifact_hashes_union": dict(parent["failed"]),
        "preserved_failed_artifact_count": len(parent["failed"]),
        "inherited_false_flags": list(parent["false_flags"]),
        "inherited_parent_true_flags": dict(parent["inherited_true"]),
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "raw_only_inherited_sources": True,
        "ancestor_decode_forbidden": True,
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def prepare(input_spec_path):
    require(not any(path.exists() for path in CHECKPOINTS), "Quiet-power checkpoint directory already contains an attempt")
    packet = _build_packet(input_spec_path)
    _write_new(SOURCE_CHECKPOINT, packet)
    return {
        "passed": True,
        "metadata_only": True,
        "component": COMPONENT,
        "source": SOURCE_CHECKPOINT.relative_to(ROOT).as_posix(),
        "source_sha256": sha_path(SOURCE_CHECKPOINT),
        "input_digest": packet["input_digest"],
        "source_union_count": packet["source_union_count"],
        "runtime_model_metadata": None,
    }


def _verify_packet(packet, parent):
    require(
        packet.get("schema") == INPUT_SCHEMA
        and packet.get("component") == COMPONENT
        and packet.get("metadata_only") is True
        and packet.get("runtime_model_metadata") is None
        and packet.get("original_source_evaluations") == 0,
        "Frozen quiet-power source packet scope changed",
    )
    require(
        packet.get("parent_report_sha256") == PARENT_REPORT_SHA256
        and packet.get("parent_result_sha256") == parent["result_sha256"]
        and packet.get("parent_result_path") == parent["result_path"]
        and packet.get("parent_commit") == PARENT_COMMIT,
        "Frozen quiet-power parent metadata changed",
    )
    require(packet.get("parent_report") == parent["report_view"] and packet.get("parent_result") == parent["result"], "Frozen parent view changed")
    code = packet.get("current_code_hashes")
    require(
        isinstance(code, dict)
        and set(code) == {RUNNER, KERNEL}
        and code == packet.get("implementation_hashes")
        and code == _code_hashes(),
        "Frozen quiet-power code changed",
    )
    spec_meta = packet.get("input_spec")
    require(isinstance(spec_meta, dict) and isinstance(spec_meta.get("path"), str) and Path(spec_meta["path"]).is_absolute() and _is_sha(spec_meta.get("raw_sha256")), "Frozen input spec metadata changed")
    spec_raw = _raw_auth(spec_meta["path"], spec_meta["raw_sha256"])
    require(len(spec_raw) == spec_meta.get("raw_bytes"), "Frozen input-spec byte count changed")
    spec = json.loads(spec_raw.decode("utf-8"))
    require(spec == packet.get("input") and spec.get("schema") == INPUT_SPEC_SCHEMA and packet.get("quiet_power_scope") == _scope_bundle(spec), "Frozen input payload or quiet scope changed")
    extra = _input_pins(spec)
    expected_stage = _merge_hashes(
        {
            PARENT_REPORT: PARENT_REPORT_SHA256,
            parent["result_path"]: parent["result_sha256"],
            RUNNER: code[RUNNER],
            KERNEL: code[KERNEL],
            spec_meta["path"]: spec_meta["raw_sha256"],
        },
        extra,
    )
    require(packet.get("stage_hashes") == expected_stage, "Frozen stage hash manifest changed")
    expected_union = _merge_hashes(parent["source_pins"], expected_stage)
    require(packet.get("source_union_hashes") == expected_union and packet.get("source_union_count") == len(expected_union) and len(expected_union) > PARENT_SOURCE_PIN_COUNT, "Frozen source union changed")
    for name, expected in expected_union.items():
        _raw_auth(name, expected)
    require(packet.get("complete_parent_checkpoint_hashes") == parent["complete"] and packet.get("complete_parent_checkpoint_count") == PARENT_CHECKPOINT_COUNT, "Frozen parent checkpoint map changed")
    require(packet.get("protected_file_hashes") == parent["protected"] and packet.get("legacy_dirty_path_hashes") == parent["dirty"], "Frozen protected or dirty maps changed")
    require(packet.get("preserved_failed_artifact_hashes_union") == parent["failed"] and packet.get("preserved_failed_artifact_count") == len(parent["failed"]) and packet.get("inherited_false_flags") == parent["false_flags"], "Frozen failure or false-flag maps changed")
    require(packet.get("inherited_parent_true_flags") == parent["inherited_true"], "Frozen inherited accepted flags changed")
    expected_digest = _digest(
        {
            "parent_report_sha256": PARENT_REPORT_SHA256,
            "parent_result_sha256": parent["result_sha256"],
            "current_code_hashes": code,
            "input_spec_sha256": spec_meta["raw_sha256"],
            "stage_hashes": expected_stage,
        }
    )
    require(packet.get("input_digest") == expected_digest, "Frozen input digest changed")


def _load_packet():
    source_raw, packet = _read_json(SOURCE_CHECKPOINT)
    parent = _parent_scope()
    _verify_packet(packet, parent)
    return source_raw, packet, parent


def _brief_marker(packet, source_sha):
    require(BRIEF.exists(), "Resume brief is missing")
    brief_raw, brief = _read_json(BRIEF)
    marker = brief.get("in_progress", {}).get(COMPONENT)
    source_rel = SOURCE_CHECKPOINT.relative_to(ROOT).as_posix()
    require(isinstance(marker, dict), "Resume brief has no quiet-power marker")
    require(
        marker.get("phase") == "prepared_unentered"
        and marker.get("source_checkpoint_sha256") == source_sha
        and marker.get("source_checkpoint_path") in (source_rel, str(SOURCE_CHECKPOINT))
        and marker.get("input_digest") == packet["input_digest"]
        and marker.get("kernel_sha256") == packet["current_code_hashes"][KERNEL]
        and marker.get("runner_sha256") == packet["current_code_hashes"][RUNNER],
        "Resume brief does not point to the prepared quiet-power packet",
    )
    return sha_bytes(brief_raw)


def _append_event(path, events, packet, stage, payload):
    event = {
        "schema": DIAGNOSTIC_SCHEMA,
        "component": COMPONENT,
        "index": len(events),
        "stage": stage,
        "input_digest": packet["input_digest"],
        "payload": _clone(payload),
    }
    with path.open("ab") as stream:
        stream.write(_encode(event))
        stream.flush()
        os.fsync(stream.fileno())
    events.append(event)


def _source_for_derive(packet, parent):
    inherited = dict(packet.get("inherited_parent_true_flags", {}))
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "input": _clone(packet["input"]),
        "input_spec": _clone(packet["input_spec"]),
        "quiet_power_scope": _clone(packet["quiet_power_scope"]),
        "parent_report": _clone(parent["report_view"]),
        "parent_report_is_authenticated_view": True,
        "parent_result": _clone(parent["result"]),
        "parent_report_metadata": _clone(packet["parent_report_metadata"]),
        "parent_result_metadata": _clone(packet["parent_result_metadata"]),
        "parent_report_sha256": packet["parent_report_sha256"],
        "parent_result_sha256": packet["parent_result_sha256"],
        "input_digest": packet["input_digest"],
        "source_union_hashes": _clone(packet["source_union_hashes"]),
        "complete_parent_checkpoint_hashes": _clone(packet["complete_parent_checkpoint_hashes"]),
        "preserved_failed_artifact_hashes_union": _clone(packet["preserved_failed_artifact_hashes_union"]),
        "inherited_false_flags": _clone(packet["inherited_false_flags"]),
        "inherited_parent_flags": inherited,
        "inherited_parent_true_flags": _clone(inherited),
        "retained_parent_data": {
            "retained_construction_data": _clone(parent["result"].get("retained_construction_data")),
            "retained_integral_term_inventory": _clone(parent["result"].get("retained_integral_term_inventory")),
            "retained_cap_provenance": _clone(parent["result"].get("retained_cap_provenance")),
        },
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
    }


def _apply_inherited_flags(result, packet):
    for key, value in packet.get("inherited_parent_true_flags", {}).items():
        if key not in CURRENT_RESULT_FIELDS:
            result.setdefault(key, value)
    return result


def _fallback(packet, error=None, calls=0, imported=False):
    result = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "metadata_only": False,
        "source_only": True,
        "input_digest": packet["input_digest"],
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "derive_entry_calls": calls,
        "derive_calls": calls,
        "entry_attempts": 1,
        "kernel_imported": imported,
        "source_packet_defined": False,
        **{key: False for key in NEW_FLAGS},
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "numeric_field_evaluations": 0,
            "old_source_consumer_calls": 0,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "old_source_normalizer_calls": 0,
            "quadratures": 0,
            "temporal_recurrence_calls": 0,
        },
    }
    _apply_inherited_flags(result, packet)
    # Reassert all current-attempt fields after applying inherited metadata.
    result.update(
        {
            "schema": RESULT_SCHEMA,
            "component": COMPONENT,
            "passed": False,
            "metadata_only": False,
            "source_only": True,
            "input_digest": packet["input_digest"],
            "runtime_model_metadata": None,
            "original_source_evaluations": 0,
            "derive_entry_calls": calls,
            "derive_calls": calls,
            "entry_attempts": 1,
            "kernel_imported": imported,
            "source_packet_defined": False,
        }
    )
    if error is not None:
        result.update(
            {
                "caught_exception": True,
                "error_type": type(error).__name__,
                "error": str(error),
                "traceback": traceback.format_exc(),
            }
        )
    for name in packet.get("inherited_false_flags", []):
        result[name] = False
    return result


def _normalise_result(value, packet, calls, imported):
    result = _clone(value) if isinstance(value, dict) else _fallback(packet, TypeError("derive returned a non-mapping result"), calls, imported)
    _apply_inherited_flags(result, packet)
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("passed", False)
    result.setdefault("metadata_only", False)
    result.setdefault("source_only", True)
    result.setdefault("input_digest", packet["input_digest"])
    result.setdefault("runtime_model_metadata", None)
    result.setdefault("original_source_evaluations", 0)
    result.setdefault("source_packet_defined", False)
    for key in NEW_FLAGS:
        result.setdefault(key, False)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    result["kernel_imported"] = imported
    return result


def _validate_quiet_result_scope(result):
    scope = result.get("quiet_power_scope")
    if scope is not None:
        require(isinstance(scope, dict), "Quiet-power result scope is not a mapping")
        require(scope.get("source") == QUIET_POWER_SOURCE_SCOPE, "Quiet-power result source scope changed")
        require(scope.get("common_kernel") == QUIET_POWER_KERNEL_SCOPE, "Quiet-power result kernel scope changed")
        require(scope.get("history_transport") == QUIET_POWER_TRANSPORT_SCOPE, "Quiet-power result transport scope changed")
    certificate = result.get("quiet_power_history_pair_certificate")
    consequence = result.get("quiet_power_history_pair_consequence")
    if result.get("passed") is True:
        require(isinstance(certificate, dict), "Passed quiet-power result lacks certificate metadata")
        require(isinstance(consequence, dict), "Passed quiet-power result lacks consequence metadata")


def _validate_result(result, packet):
    require(
        isinstance(result, dict)
        and result.get("schema") == RESULT_SCHEMA
        and result.get("component") == COMPONENT
        and result.get("metadata_only") is False
        and result.get("source_only") is True
        and isinstance(result.get("passed"), bool),
        "Quiet-power result schema or scope changed",
    )
    require(
        result.get("input_digest") == packet["input_digest"]
        and result.get("runtime_model_metadata") is None
        and result.get("original_source_evaluations") == 0
        and result.get("entry_attempts") == 1,
        "Quiet-power execution scope changed",
    )
    calls = result.get("derive_entry_calls")
    require(isinstance(calls, int) and not isinstance(calls, bool) and calls in (0, 1) and result.get("derive_calls") == calls and isinstance(result.get("kernel_imported"), bool), "Quiet-power derive bookkeeping changed")
    if result.get("passed") is True:
        require(calls == 1 and result.get("kernel_imported") is True and result.get("source_packet_defined") is True, "Passed quiet-power result requires one imported derive and a source packet")
        for key in NEW_FLAGS:
            require(result.get(key) is True, "Passed quiet-power result lacks narrow flag: " + key)
    for key in NEW_FLAGS:
        require(isinstance(result.get(key), bool), "Quiet-power narrow flag is not boolean: " + key)
    forbidden = (
        "actual_C3_C5_incoming_source_identity_proved",
        "actual_C3_C5_control_source_identity_proved",
        "actual_C3_C5_profile_identity_proved",
        "actual_C3_cap_to_C5_source_bound_proved",
        "compact_volume_cap_same_source_bound_proved",
        "global_five_moment_closure",
        "numeric_values_computed",
        "full_goal_complete",
    )
    for key in forbidden:
        if key in result:
            require(result[key] is False, "Forbidden completion flag is true: " + key)
    for key in packet.get("inherited_false_flags", []):
        require(result.get(key) is False, "Inherited false flag changed: " + key)
    for key, value in packet.get("inherited_parent_true_flags", {}).items():
        if key not in CURRENT_RESULT_FIELDS:
            require(result.get(key) is value, "Inherited accepted flag changed: " + key)
    audit = result.get("source_inspection_audit")
    if result.get("passed") is True:
        require(isinstance(audit, dict), "Passed quiet-power result lacks source inspection audit")
        for key in (
            "new_native_quiet_density_rows",
            "new_C3_quiet_density_rows",
            "new_common_kernel_rows",
            "new_signed_quiet_transport_rows",
        ):
            require(audit.get(key) == 20 and isinstance(audit.get(key), int) and not isinstance(audit.get(key), bool), "Passed quiet-power audit count is not 20: " + key)
    if isinstance(audit, dict):
        for key in (
            "new_native_quiet_density_rows",
            "new_C3_quiet_density_rows",
            "new_common_kernel_rows",
            "new_signed_quiet_transport_rows",
        ):
            if key in audit:
                require(isinstance(audit[key], int) and not isinstance(audit[key], bool) and audit[key] >= 0, "Malformed quiet-power audit count: " + key)
        for key in (
            "original_source_evaluations",
            "numeric_field_evaluations",
            "old_source_consumer_calls",
            "old_proof_functions_called",
            "old_operator_calls",
            "old_source_normalizer_calls",
            "quadratures",
            "temporal_recurrence_calls",
        ):
            if key in audit and audit[key] is not None:
                require(audit[key] == 0, "Forbidden source call recorded: " + key)
    _validate_quiet_result_scope(result)


def _load_saved(packet):
    require(ENTERED_CHECKPOINT.exists() and RESULT_CHECKPOINT.exists() and DIAGNOSTIC_CHECKPOINT.exists(), "Saved quiet-power attempt is incomplete")
    entered_raw, entered = _read_json(ENTERED_CHECKPOINT)
    result_raw, result = _read_json(RESULT_CHECKPOINT)
    journal_raw = DIAGNOSTIC_CHECKPOINT.read_bytes()
    source_sha = sha_path(SOURCE_CHECKPOINT)
    require(entered.get("schema") == ENTERED_SCHEMA and entered.get("component") == COMPONENT and entered.get("source_checkpoint_sha256") == source_sha and entered.get("input_digest") == packet["input_digest"] and entered.get("kernel_imported_before_entry") is False and entered.get("derive_calls_before_entry") == 0 and entered.get("entry_attempts") == 1, "Entered quiet-power checkpoint changed")
    require(result.get("source_checkpoint_sha256") == source_sha and result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw) and result.get("diagnostic_checkpoint_sha256") == sha_bytes(journal_raw), "Saved quiet-power provenance changed")
    _validate_result(result, packet)
    return entered_raw, result_raw, result, journal_raw


def run():
    require(SOURCE_CHECKPOINT.exists() and not ENTERED_CHECKPOINT.exists() and not RESULT_CHECKPOINT.exists() and not DIAGNOSTIC_CHECKPOINT.exists(), "Frozen quiet-power source is absent or attempt already entered")
    source_raw, packet, parent = _load_packet()
    brief_sha = _brief_marker(packet, sha_bytes(source_raw))
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "source_checkpoint_path": SOURCE_CHECKPOINT.relative_to(ROOT).as_posix(),
        "input_digest": packet["input_digest"],
        "current_code_hashes": dict(packet["current_code_hashes"]),
        "brief_sha256_before_entry": brief_sha,
        "kernel_imported_before_entry": False,
        "derive_calls_before_entry": 0,
        "entry_attempts": 1,
        "entered_code_immutable": True,
        "never_retry_entered_code": True,
        "runtime_model_metadata": None,
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    # Durable entry and empty journal are committed before the sole import.
    _write_new(ENTERED_CHECKPOINT, entered)
    _write_new(DIAGNOSTIC_CHECKPOINT, b"")
    events = []

    def journal_callback(event=None, *args, **fields):
        if isinstance(event, dict):
            payload = dict(event)
            payload.update(fields)
        elif event is None:
            payload = dict(fields)
        else:
            payload = {"event": event, **fields}
        if args:
            payload["positional_args"] = list(args)
        _append_event(DIAGNOSTIC_CHECKPOINT, events, packet, "kernel_callback", payload)

    calls = 0
    imported = False
    try:
        source = _source_for_derive(packet, parent)
        spec = importlib.util.spec_from_file_location("_original_compact_quiet_power_history_pair_entry", _resolve(KERNEL))
        require(spec is not None and spec.loader is not None, "Cannot create quiet-power kernel import specification")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        imported = True
        derive = getattr(module, "derive", None)
        require(callable(derive), "New quiet-power kernel lacks derive(source, journal_callback)")
        calls = 1
        value = derive(source, journal_callback)
        result = _normalise_result(value, packet, calls, imported)
        _append_event(DIAGNOSTIC_CHECKPOINT, events, packet, "derive_returned", {"derive_entry_calls": calls, "passed": result.get("passed"), "runtime_model_metadata": None})
    except Exception as error:
        result = _fallback(packet, error, calls, imported)
        _append_event(DIAGNOSTIC_CHECKPOINT, events, packet, "derive_exception", {"error_type": type(error).__name__, "error": str(error), "derive_entry_calls": calls})
    result["source_checkpoint_sha256"] = sha_path(SOURCE_CHECKPOINT)
    result["entry_checkpoint_sha256"] = sha_path(ENTERED_CHECKPOINT)
    result["diagnostic_checkpoint_sha256"] = sha_path(DIAGNOSTIC_CHECKPOINT)
    result["diagnostic_event_count"] = len(events)
    result["diagnostic_summary"] = {"event_count": len(events), "journal_sha256": result["diagnostic_checkpoint_sha256"]}
    result["kernel_imported"] = imported
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    result["runtime_model_metadata"] = None
    # The result is durable before any validation can reject it.
    _write_new(RESULT_CHECKPOINT, result)
    _validate_result(result, packet)
    return {"passed": result.get("passed") is True, "metadata_only": True, "component": COMPONENT, "result": RESULT_CHECKPOINT.relative_to(ROOT).as_posix(), "result_sha256": sha_path(RESULT_CHECKPOINT), "derive_entry_calls": result.get("derive_entry_calls"), "runtime_model_metadata": None}


def inspect():
    source_raw, packet, _ = _load_packet()
    entered_raw, result_raw, result, journal_raw = _load_saved(packet)
    return {"passed": result.get("passed") is True, "metadata_only": True, "component": COMPONENT, "source_sha256": sha_bytes(source_raw), "entered_sha256": sha_bytes(entered_raw), "result_sha256": sha_bytes(result_raw), "diagnostics_sha256": sha_bytes(journal_raw), "derive_entry_calls": result.get("derive_entry_calls"), "runtime_model_metadata": None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "inspect"))
    parser.add_argument("--input-spec", help="Absolute JSON input spec for prepare")
    args = parser.parse_args()
    if args.action == "prepare":
        require(args.input_spec is not None, "prepare requires --input-spec")
        value = prepare(args.input_spec)
    else:
        require(args.input_spec is None, "--input-spec is only valid with prepare")
        value = globals()[args.action]()
    print(json.dumps(value, ensure_ascii=True, allow_nan=False))


if __name__ == "__main__":
    main()


