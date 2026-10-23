"""Bounded durable runner for the cached C2 volume-source leaf.

``prepare`` authenticates the saved parent report/result, one selected signed
pullback report, inherited raw pins, checkpoints, and the two new code blobs,
then writes one frozen source packet.  ``run`` enters once before importing
only the new kernel and invokes ``derive(source, journal_callback)`` once.
``inspect`` reads saved metadata only.  No native archive, legacy artifact,
old checker, normalizer, provider, quadrature, or proof implementation is
decoded or executed here.
"""

from __future__ import annotations

from datetime import datetime, timezone
import argparse
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
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "compact_cached_C2_volume_source"
RUNNER = "scripts/run_original_compact_cached_C2_volume_source.py"
KERNEL = "src/openai_ns_reconstruction/original_compact_cached_C2_volume_source.py"
PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_volume_error_source_packet.json"
)
PARENT_REPORT_SHA256 = (
    "b554f056f73b1f038d9c12b162e95d24b50d1e7c90173364496e09fd2b74029f"
)
PARENT_RESULT = "tmp/original_compact_volume_error_source_packet/result.json"
PARENT_RESULT_SHA256 = (
    "809ae44b56ee4428da44708630a72a38f0538102671b9ab05d5f404ac4cffb5d"
)
CURRENT_PARENT_COMMIT = "0ce549d64c97f9d36f047eac57a6fa6359c4bf12"
SIGNED_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_signed_pullback_difference.json"
)
SIGNED_REPORT_SHA256 = (
    "4a44b291afb4a0de8f82c652b0e99b4cd6e02eece8c86f44e2e94c76c5e4b9f7"
)

INPUT_SCHEMA = "original_compact_cached_C2_volume_source_inputs_v1"
ENTERED_SCHEMA = "original_compact_cached_C2_volume_source_entered_v1"
RESULT_SCHEMA = "original_compact_cached_C2_volume_source"
DIAGNOSTIC_SCHEMA = "compact_cached_C2_volume_source_stage_diagnostic_v1"

CHECKPOINT_DIR = ROOT / "tmp/original_compact_cached_C2_volume_source"
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"
CHECKPOINTS = {
    "source": SOURCE_CHECKPOINT,
    "entered": ENTERED_CHECKPOINT,
    "result": RESULT_CHECKPOINT,
    "diagnostics": DIAGNOSTIC_CHECKPOINT,
}

RAW_ONLY_MARKERS = ("resume_first_switch", "rm_loop_density_identity")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def sha_path(path):
    return sha_bytes(Path(path).read_bytes())


def _encode(value):
    return (
        json.dumps(value, ensure_ascii=True, sort_keys=True,
                   allow_nan=False, separators=(",", ":")).encode("utf-8")
        + b"\n"
    )


def _clone(value):
    return json.loads(_encode(value).decode("utf-8"))


def _digest(value):
    return sha_bytes(json.dumps(value, ensure_ascii=True, sort_keys=True,
                                allow_nan=False, separators=(",", ":"))
                     .encode("utf-8"))


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
    require(isinstance(expected, str) and len(expected) == 64,
            "Malformed SHA-256 for " + str(name))
    path = _resolve(name)
    require(path.exists(), "Missing authenticated artifact: " + str(name))
    raw = path.read_bytes()
    require(sha_bytes(raw) == expected,
            "Changed authenticated artifact: " + str(name))
    return raw


def _raw_name_forbidden(name):
    normalized = str(name).replace("-", "_").lower()
    return any(marker in normalized for marker in RAW_ONLY_MARKERS)


def _is_sha(value):
    return isinstance(value, str) and len(value) == 64


def _hash_map(value, label, authenticate=True):
    require(isinstance(value, dict), label + " is not a mapping")
    result = {}
    for name, expected in value.items():
        require(isinstance(name, str) and _is_sha(expected),
                label + " contains a malformed SHA-256")
        if authenticate:
            _raw_auth(name, expected)
        result[name] = expected
    return result


def _flatten_hashes(value, label):
    """Flatten flat or nested checkpoint maps, merging identical repeats."""
    flattened = {}

    def add(name, expected):
        require(isinstance(name, str) and _is_sha(expected),
                label + " contains a malformed SHA-256")
        if name in flattened:
            require(flattened[name] == expected,
                    "Conflicting checkpoint hash collision at " + name)
        else:
            flattened[name] = expected

    def visit(node, path):
        require(isinstance(node, dict), path + " is not a mapping")
        if not node:
            return
        if all(isinstance(key, str) and _is_sha(child)
               for key, child in node.items()):
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
    return flattened


def _merge_hashes(*maps):
    merged = {}
    for mapping in maps:
        for name, expected in mapping.items():
            if name in merged:
                require(merged[name] == expected,
                        "Conflicting authenticated hash at " + name)
            else:
                merged[name] = expected
    return merged


def _parent_scope():
    report_raw = _raw_auth(PARENT_REPORT, PARENT_REPORT_SHA256)
    result_raw = _raw_auth(PARENT_RESULT, PARENT_RESULT_SHA256)
    report = json.loads(report_raw.decode("utf-8"))
    result = json.loads(result_raw.decode("utf-8"))
    require(
        report.get("schema")
        == "original_compact_volume_error_source_packet_publication_v1"
        and report.get("component") == "compact_volume_error_source_packet"
        and report.get("runtime_model_metadata") is None
        and report.get("source_pin_count") == 1910
        and report.get("checkpoint_count") == 184,
        "Parent volume-error report scope changed",
    )
    require(
        report.get("full_result_raw_sha256") == PARENT_RESULT_SHA256
        and report.get("full_result_raw_bytes") == len(result_raw)
        and report.get("full_result") == result,
        "Parent full result does not equal the saved result checkpoint",
    )
    require(
        result.get("schema") == "original_compact_volume_error_source_packet"
        and result.get("component") == "compact_volume_error_source_packet"
        and result.get("runtime_model_metadata") is None
        and result.get("source_ftc_transport_defined") is True
        and result.get("source_packet_defined") is False
        and result.get("cap_to_source_identity_proved") is False
        and result.get("quantitative_packet_error_reduction_proved") is False,
        "Parent result scope changed",
    )
    source_pins = _hash_map(report.get("source_pins"),
                            "Parent source pins")
    require(len(source_pins) == report.get("source_pin_count") == 1910,
            "Parent source-pin count changed")
    inventory = _flatten_hashes(report.get("parent_checkpoint_inventory"),
                                "Parent checkpoint inventory")
    own = _flatten_hashes({"checkpoint_hashes": report.get("checkpoint_hashes")},
                          "Parent report checkpoint hashes")
    immediate_value = report.get("parent_checkpoint_hashes")
    immediate = {} if immediate_value is None else _hash_map(
        immediate_value, "Parent immediate checkpoint hashes"
    )
    require(len(inventory) == 180 and len(own) == 4,
            "Parent checkpoint map sizes changed")
    complete = _merge_hashes(inventory, immediate, own)
    require(len(complete) == 184,
            "Complete inherited checkpoint count changed")
    for name, expected in complete.items():
        _raw_auth(name, expected)
    failed = _hash_map(report.get("preserved_failed_artifact_hashes_union"),
                       "Parent failed-artifact union")
    protected = _hash_map(report.get("protected_file_hashes"),
                          "Parent protected files")
    dirty = _hash_map(report.get("legacy_dirty_path_hashes"),
                      "Parent dirty files")
    require(len(protected) == 16 and len(dirty) == 2,
            "Parent protected/dirty mapping changed")
    return {
        "report_raw": report_raw,
        "report": report,
        "result_raw": result_raw,
        "result": result,
        "source_pins": source_pins,
        "inventory": inventory,
        "own": own,
        "immediate": immediate,
        "complete": complete,
        "failed": failed,
        "protected": protected,
        "dirty": dirty,
    }


def _code_hashes():
    result = {}
    for name in (RUNNER, KERNEL):
        path = _resolve(name)
        require(path.exists(), "New code blob is missing: " + name)
        result[name] = sha_path(path)
    return result


def _signed_record():
    raw = _raw_auth(SIGNED_REPORT, SIGNED_REPORT_SHA256)
    report = json.loads(raw.decode("utf-8"))
    require(report.get("schema") == "original_signed_pullback_difference_stage_v1"
            and report.get("component") == "signed_pullback_difference"
            and isinstance(report.get("component_result"), dict),
            "Signed pullback source record shape changed")
    component = report["component_result"]
    require(component.get("schema") == "original_signed_pullback_difference"
            and component.get("component") == "signed_pullback_difference"
            and component.get("passed") is True,
            "Signed pullback component result is not the saved accepted record")
    return raw, report, component


def _build_packet():
    parent = _parent_scope()
    code = _code_hashes()
    # Authenticate the selected JSON source record but defer its decode until
    # run, after the new entry checkpoint is durable.
    signed_raw = _raw_auth(SIGNED_REPORT, SIGNED_REPORT_SHA256)
    stage = {
        PARENT_REPORT: PARENT_REPORT_SHA256,
        PARENT_RESULT: PARENT_RESULT_SHA256,
        SIGNED_REPORT: SIGNED_REPORT_SHA256,
        **code,
    }
    source_union = _merge_hashes(parent["source_pins"], stage)
    for name, expected in source_union.items():
        _raw_auth(name, expected)
    for name, expected in parent["protected"].items():
        _raw_auth(name, expected)
    for name, expected in parent["dirty"].items():
        _raw_auth(name, expected)
    input_digest = _digest({
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_sha256": PARENT_RESULT_SHA256,
        "signed_report_sha256": SIGNED_REPORT_SHA256,
        "current_code_hashes": code,
    })
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_packet_definition_only": True,
        "cap_to_source_identity_proved": False,
        "quantitative_packet_error_reduction_proved": False,
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result": PARENT_RESULT,
        "parent_result_sha256": PARENT_RESULT_SHA256,
        "signed_report": SIGNED_REPORT,
        "signed_report_sha256": SIGNED_REPORT_SHA256,
        "signed_component_selector": "$.component_result",
        "parent_commit": CURRENT_PARENT_COMMIT,
        "input_digest": input_digest,
        "inherited_source_pins": dict(parent["source_pins"]),
        "inherited_source_pin_count": len(parent["source_pins"]),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "stage_hashes": stage,
        "current_code_hashes": code,
        "implementation_hashes": dict(code),
        "parent_checkpoint_inventory": dict(parent["inventory"]),
        "parent_checkpoint_hashes": dict(parent["immediate"]),
        "report_checkpoint_hashes": dict(parent["own"]),
        "complete_parent_checkpoint_hashes": dict(parent["complete"]),
        "parent_checkpoint_inventory_count": len(parent["inventory"]),
        "report_checkpoint_hashes_count": len(parent["own"]),
        "complete_parent_checkpoint_count": len(parent["complete"]),
        "parent_checkpoint_count": len(parent["complete"]),
        "preserved_failed_artifact_hashes_union": dict(parent["failed"]),
        "preserved_failed_artifact_count": len(parent["failed"]),
        "protected_file_hashes": dict(parent["protected"]),
        "legacy_dirty_path_hashes": dict(parent["dirty"]),
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "raw_only_inherited_sources": True,
        "history_decode_forbidden": True,
        "selected_signed_record_raw_bytes": len(signed_raw),
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def prepare():
    require(not any(path.exists() for path in CHECKPOINTS.values()),
            "Cached C2 checkpoint directory already contains an attempt")
    packet = _build_packet()
    _write_new(SOURCE_CHECKPOINT, packet)
    return {
        "passed": True,
        "metadata_only": True,
        "component": COMPONENT,
        "source": str(SOURCE_CHECKPOINT.relative_to(ROOT).as_posix()),
        "source_sha256": sha_path(SOURCE_CHECKPOINT),
        "input_digest": packet["input_digest"],
        "source_union_count": packet["source_union_count"],
        "runtime_model_metadata": None,
    }


def _verify_packet(packet, parent):
    require(packet.get("schema") == INPUT_SCHEMA
            and packet.get("component") == COMPONENT
            and packet.get("metadata_only") is True
            and packet.get("source_packet_definition_only") is True
            and packet.get("cap_to_source_identity_proved") is False
            and packet.get("quantitative_packet_error_reduction_proved") is False
            and packet.get("runtime_model_metadata") is None
            and packet.get("original_source_evaluations") == 0,
            "Frozen cached C2 source packet scope changed")
    require(packet.get("parent_report_sha256") == PARENT_REPORT_SHA256
            and packet.get("parent_result_sha256") == PARENT_RESULT_SHA256
            and packet.get("signed_report_sha256") == SIGNED_REPORT_SHA256
            and packet.get("parent_commit") == CURRENT_PARENT_COMMIT,
            "Frozen parent/source hashes changed")
    code = packet.get("current_code_hashes")
    require(code == packet.get("implementation_hashes")
            and isinstance(code, dict)
            and set(code) == {RUNNER, KERNEL},
            "Frozen cached C2 code map changed")
    live = _code_hashes()
    require(code == live, "New code changed after prepare")
    expected_stage = {
        PARENT_REPORT: PARENT_REPORT_SHA256,
        PARENT_RESULT: PARENT_RESULT_SHA256,
        SIGNED_REPORT: SIGNED_REPORT_SHA256,
        **code,
    }
    require(packet.get("stage_hashes") == expected_stage,
            "Frozen cached C2 stage manifest changed")
    pins = packet.get("source_union_hashes")
    require(isinstance(pins, dict)
            and packet.get("source_union_count") == len(pins),
            "Frozen source union changed")
    expected_pins = _merge_hashes(parent["source_pins"], expected_stage)
    require(pins == expected_pins,
            "Frozen cached C2 source union does not match its stages")
    for name, expected in pins.items():
        _raw_auth(name, expected)
    require(packet.get("inherited_source_pins") == parent["source_pins"],
            "Frozen inherited source pins changed")
    require(packet.get("parent_checkpoint_inventory") == parent["inventory"]
            and packet.get("parent_checkpoint_hashes") == parent["immediate"]
            and packet.get("report_checkpoint_hashes") == parent["own"]
            and packet.get("complete_parent_checkpoint_hashes") == parent["complete"]
            and packet.get("complete_parent_checkpoint_count") == 184,
            "Frozen checkpoint maps changed")
    require(packet.get("protected_file_hashes") == parent["protected"]
            and packet.get("legacy_dirty_path_hashes") == parent["dirty"],
            "Frozen protected or dirty maps changed")


def _load_packet(decode_signed=False):
    source_raw, packet = _read_json(SOURCE_CHECKPOINT)
    parent = _parent_scope()
    _verify_packet(packet, parent)
    signed = None
    if decode_signed:
        signed_raw, signed_report, signed_component = _signed_record()
        require(sha_bytes(signed_raw) == packet["signed_report_sha256"],
                "Signed source record changed after preparation")
        signed = {
            "raw": signed_raw,
            "report": signed_report,
            "component": signed_component,
        }
    else:
        _raw_auth(SIGNED_REPORT, packet["signed_report_sha256"])
    return source_raw, packet, parent, signed


def _brief_marker(packet, source_sha):
    require(BRIEF.exists(), "Resume brief is missing")
    brief_raw, brief = _read_json(BRIEF)
    node = brief.get("in_progress", {}).get(COMPONENT)
    require(isinstance(node, dict), "Resume brief has no cached C2 marker")
    source_rel = SOURCE_CHECKPOINT.relative_to(ROOT).as_posix()
    require(node.get("phase") == "prepared_unentered"
            and node.get("source_checkpoint_sha256") == source_sha
            and node.get("source_checkpoint_path") in (source_rel,
                                                          str(SOURCE_CHECKPOINT))
            and node.get("kernel_sha256") == packet["current_code_hashes"][KERNEL]
            and node.get("runner_sha256") == packet["current_code_hashes"][RUNNER],
            "Resume brief does not point to the prepared cached C2 packet")
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


def _source_for_derive(packet, parent, signed):
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "parent_report": _clone(parent["report"]),
        "parent_result": _clone(parent["result"]),
        "signed_component": _clone(signed["component"]),
        "signed_report_metadata": {
            "path": SIGNED_REPORT,
            "raw_sha256": packet["signed_report_sha256"],
            "selector": packet["signed_component_selector"],
            "raw_authenticated_before_decode": True,
        },
        "parent_report_sha256": packet["parent_report_sha256"],
        "parent_result_sha256": packet["parent_result_sha256"],
        "input_digest": packet["input_digest"],
        "source_union_hashes": _clone(packet["source_union_hashes"]),
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_definition_only": True,
        "cap_to_source_identity_proved": False,
        "quantitative_packet_error_reduction_proved": False,
    }


def _fallback(packet, error=None, calls=0, imported=False):
    result = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "metadata_only": False,
        "input_digest": packet["input_digest"],
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_definition_only": True,
        "source_ftc_transport_defined": True,
        "cap_to_source_identity_proved": False,
        "quantitative_packet_error_reduction_proved": False,
        "derive_entry_calls": calls,
        "derive_calls": calls,
        "entry_attempts": 1,
        "kernel_imported": imported,
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "old_source_consumer_calls": 0,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "old_source_normalizer_calls": 0,
            "runtime_model_metadata": None,
        },
    }
    if error is not None:
        result.update({
            "caught_exception": True,
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": traceback.format_exc(),
        })
    return result


def _normalise_result(value, packet, calls, imported):
    if isinstance(value, dict):
        result = _clone(value)
    else:
        result = _fallback(packet,
                           TypeError("derive returned a non-mapping result"),
                           calls, imported)
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("passed", False)
    result.setdefault("metadata_only", False)
    result.setdefault("input_digest", packet["input_digest"])
    result.setdefault("runtime_model_metadata", None)
    result.setdefault("original_source_evaluations", 0)
    result.setdefault("source_packet_definition_only", True)
    result.setdefault("source_ftc_transport_defined", True)
    result.setdefault("cap_to_source_identity_proved", False)
    result.setdefault("quantitative_packet_error_reduction_proved", False)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    result["kernel_imported"] = imported
    return result


def _validate_result(result, packet):
    require(isinstance(result, dict)
            and result.get("schema") == RESULT_SCHEMA
            and result.get("component") == COMPONENT
            and result.get("metadata_only") is False,
            "Cached C2 result schema or metadata scope changed")
    require(result.get("input_digest") == packet["input_digest"]
            and result.get("runtime_model_metadata") is None
            and result.get("original_source_evaluations") == 0
            and result.get("source_packet_definition_only") is True
            and result.get("cap_to_source_identity_proved") is False
            and result.get("quantitative_packet_error_reduction_proved") is False
            and result.get("entry_attempts") == 1,
            "Cached C2 execution scope changed")
    calls = result.get("derive_entry_calls")
    require(isinstance(calls, int) and not isinstance(calls, bool)
            and calls in (0, 1)
            and result.get("derive_calls") == calls
            and isinstance(result.get("kernel_imported"), bool),
            "Cached C2 derive bookkeeping changed")
    if result.get("passed") is True:
        require(calls == 1 and result.get("kernel_imported") is True,
                "Passed result requires one successful derive call")
    audit = result.get("source_inspection_audit")
    if isinstance(audit, dict):
        for key in ("original_source_evaluations", "old_source_consumer_calls",
                    "old_proof_functions_called", "old_operator_calls",
                    "old_source_normalizer_calls"):
            if key in audit and audit[key] is not None:
                require(audit[key] == 0,
                        "Forbidden source consumer/operator call recorded")


def _load_saved(packet):
    require(ENTERED_CHECKPOINT.exists() and RESULT_CHECKPOINT.exists()
            and DIAGNOSTIC_CHECKPOINT.exists(),
            "Saved cached C2 attempt is incomplete")
    entered_raw, entered = _read_json(ENTERED_CHECKPOINT)
    result_raw, result = _read_json(RESULT_CHECKPOINT)
    journal_raw = DIAGNOSTIC_CHECKPOINT.read_bytes()
    source_sha = sha_path(SOURCE_CHECKPOINT)
    require(entered.get("schema") == ENTERED_SCHEMA
            and entered.get("component") == COMPONENT
            and entered.get("source_checkpoint_sha256") == source_sha
            and entered.get("kernel_imported_before_entry") is False
            and entered.get("derive_calls_before_entry") == 0
            and entered.get("entry_attempts") == 1,
            "Entered cached C2 checkpoint changed")
    require(result.get("source_checkpoint_sha256") == source_sha
            and result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw)
            and result.get("diagnostic_checkpoint_sha256") == sha_bytes(journal_raw),
            "Saved cached C2 result provenance changed")
    _validate_result(result, packet)
    return {
        "entered_raw": entered_raw,
        "entered": entered,
        "result_raw": result_raw,
        "result": result,
        "journal_raw": journal_raw,
    }


def run():
    require(SOURCE_CHECKPOINT.exists()
            and not ENTERED_CHECKPOINT.exists()
            and not RESULT_CHECKPOINT.exists()
            and not DIAGNOSTIC_CHECKPOINT.exists(),
            "Frozen source is absent or cached C2 attempt already entered")
    source_raw, packet, parent, signed = _load_packet(decode_signed=True)
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
        _append_event(DIAGNOSTIC_CHECKPOINT, events, packet,
                      "kernel_callback", payload)

    calls = 0
    imported = False
    try:
        source = _source_for_derive(packet, parent, signed)
        spec = importlib.util.spec_from_file_location(
            "_original_compact_cached_C2_volume_source_entry",
            _resolve(KERNEL),
        )
        require(spec is not None and spec.loader is not None,
                "Cannot create cached C2 kernel import specification")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        imported = True
        derive = getattr(module, "derive", None)
        require(callable(derive),
                "New cached C2 kernel lacks derive(source, journal_callback)")
        calls = 1
        value = derive(source, journal_callback)
        result = _normalise_result(value, packet, calls, imported)
        _append_event(DIAGNOSTIC_CHECKPOINT, events, packet, "derive_returned", {
            "derive_entry_calls": calls,
            "passed": result.get("passed"),
            "runtime_model_metadata": None,
        })
    except Exception as error:
        result = _fallback(packet, error, calls, imported)
        _append_event(DIAGNOSTIC_CHECKPOINT, events, packet, "derive_exception", {
            "error_type": type(error).__name__,
            "error": str(error),
            "derive_entry_calls": calls,
        })
    result["source_checkpoint_sha256"] = sha_path(SOURCE_CHECKPOINT)
    result["entry_checkpoint_sha256"] = sha_path(ENTERED_CHECKPOINT)
    result["diagnostic_checkpoint_sha256"] = sha_path(DIAGNOSTIC_CHECKPOINT)
    result["diagnostic_event_count"] = len(events)
    result["kernel_imported"] = imported
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    result["runtime_model_metadata"] = None
    _write_new(RESULT_CHECKPOINT, result)
    # The complete result and journal are durable before metadata validation.
    _validate_result(result, packet)
    return {
        "passed": result.get("passed") is True,
        "metadata_only": True,
        "component": COMPONENT,
        "result": RESULT_CHECKPOINT.relative_to(ROOT).as_posix(),
        "result_sha256": sha_path(RESULT_CHECKPOINT),
        "derive_entry_calls": result.get("derive_entry_calls"),
        "runtime_model_metadata": None,
    }


def inspect():
    source_raw, packet, _, _ = _load_packet(decode_signed=False)
    saved = _load_saved(packet)
    return {
        "passed": saved["result"].get("passed") is True,
        "metadata_only": True,
        "component": COMPONENT,
        "source_sha256": sha_bytes(source_raw),
        "entered_sha256": sha_bytes(saved["entered_raw"]),
        "result_sha256": sha_bytes(saved["result_raw"]),
        "diagnostics_sha256": sha_bytes(saved["journal_raw"]),
        "derive_entry_calls": saved["result"].get("derive_entry_calls"),
        "runtime_model_metadata": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "inspect"))
    args = parser.parse_args()
    value = globals()[args.action]()
    print(json.dumps(value, ensure_ascii=True, allow_nan=False))


if __name__ == "__main__":
    main()
