"""Bounded durable runner for the compact signed-history defect transport leaf.

The runner freezes authenticated metadata and one external input specification,
then (when explicitly invoked later) enters exactly once before importing only
the new kernel and calling ``derive(source, journal_callback)`` once.  It does
not decode or call ancestor source consumers, proofs, normalizers, providers,
quadratures, or recurrence code.  This draft is intentionally unexecuted.
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
BRIEF = ROOT.parents[1] / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "compact_signed_history_defect_transport"
RUNNER = "scripts/run_original_compact_signed_history_defect_transport.py"
KERNEL = "src/openai_ns_reconstruction/original_compact_signed_history_defect_transport.py"
PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_C3_outer_amplitude_identity.json"
)
PARENT_REPORT_SHA256 = (
    "987e1aba8d33a22b325e1f0386a00e522778d296b65f98ac7d826da228292975"
)
PARENT_RESULT_SHA256 = (
    "b34968d776189f7c259fe43105123a997ac4eff57c474239c40cc565e5ed62ed"
)
PARENT_COMMIT = "870f95489b93c91f8f3dabdba0ddf283ffe6c927"
PARENT_REPORT_PARENT_COMMIT = "bb95f192f244ff4f1e07068cf19191eaf296be76"

INPUT_SCHEMA = "original_compact_signed_history_defect_transport_inputs_v1"
INPUT_SPEC_SCHEMA = "compact_signed_history_defect_transport_input_v1"
ENTERED_SCHEMA = "original_compact_signed_history_defect_transport_entered_v1"
RESULT_SCHEMA = "original_compact_signed_history_defect_transport"
DIAGNOSTIC_SCHEMA = "compact_signed_history_defect_transport_stage_diagnostic_v1"
INHERITED_FALSE_FLAG_COUNT = 26

CHECKPOINT_DIR = ROOT / "tmp/original_compact_signed_history_defect_transport"
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


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def sha_path(path):
    return sha_bytes(Path(path).read_bytes())


def _encode(value):
    return (
        json.dumps(value, ensure_ascii=True, sort_keys=True, allow_nan=False,
                   separators=(",", ":")).encode("utf-8") + b"\n"
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


def _is_sha(value):
    return isinstance(value, str) and len(value) == 64 and all(
        character in "0123456789abcdef" for character in value
    )


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
    """Flatten nested checkpoint maps, allowing equal repeated hashes."""
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


def _parent_result_reference(report):
    selector = report.get("full_result_artifact_selector")
    prefix = "$.lossless_artifacts[path="
    require(isinstance(selector, str) and selector.startswith(prefix)
            and selector.endswith("]"), "Parent result selector is malformed")
    path = selector[len(prefix):-1]
    expected = report.get("full_result_raw_sha256")
    require(path and _is_sha(expected), "Parent result selector or SHA is malformed")
    return path, expected


def _parent_report_view(report):
    excluded = [
        "lossless_artifacts", "full_result", "component_result",
        "full_result_raw_bytes", "full_result_raw_sha256",
    ]
    view = {key: _clone(value) for key, value in report.items()
            if key not in excluded}
    view["parent_report_is_authenticated_view"] = True
    view["parent_report_view_excluded_fields"] = list(excluded)
    return view


def _parent_scope():
    report_raw = _raw_auth(PARENT_REPORT, PARENT_REPORT_SHA256)
    report = json.loads(report_raw.decode("utf-8"))
    require(
        report.get("schema")
        == "original_compact_C3_outer_amplitude_identity_publication_v1"
        and report.get("component") == "compact_C3_outer_amplitude_identity"
        and report.get("runtime_model_metadata") is None
        and report.get("source_pin_count") == 2126
        and report.get("checkpoint_count") == 224
        and report.get("parent_commit") == PARENT_REPORT_PARENT_COMMIT
        and report.get("attempt_status") == "passed",
        "Parent C3 outer-amplitude report scope changed",
    )
    result_path, result_sha = _parent_result_reference(report)
    require(result_sha == PARENT_RESULT_SHA256,
            "Parent result SHA changed")
    result_raw = _raw_auth(result_path, result_sha)
    result = json.loads(result_raw.decode("utf-8"))
    require(
        report.get("full_result_raw_bytes") == len(result_raw)
        and report.get("full_result") == result
        and report.get("component_result") == result,
        "Parent full result does not equal the saved result checkpoint",
    )
    require(
        result.get("schema") == "original_compact_C3_outer_amplitude_identity"
        and result.get("component") == "compact_C3_outer_amplitude_identity"
        and result.get("passed") is True
        and result.get("metadata_only") is False
        and result.get("runtime_model_metadata") is None
        and result.get("original_source_evaluations") == 0,
        "Parent C3 outer-amplitude result scope changed",
    )
    false_flags = report.get("false_flags")
    require(
        isinstance(false_flags, list)
        and len(false_flags) == INHERITED_FALSE_FLAG_COUNT
        and len(set(false_flags)) == INHERITED_FALSE_FLAG_COUNT
        and all(isinstance(name, str) and name for name in false_flags),
        "Parent inherited false-flag list changed",
    )
    for name in false_flags:
        require(report.get(name) is False and result.get(name) is False,
                "Parent inherited false flag is not false: " + name)
    source_pins = _hash_map(report.get("source_pins"), "Parent source pins")
    require(len(source_pins) == report.get("source_pin_count") == 2126,
            "Parent source-pin count changed")
    inventory = _flatten_hashes(report.get("parent_checkpoint_inventory"),
                                "Parent checkpoint inventory")
    own = _flatten_hashes({"checkpoint_hashes": report.get("checkpoint_hashes")},
                          "Parent report checkpoint hashes")
    require(len(inventory) == report.get("inherited_checkpoint_count") == 220,
            "Parent inherited checkpoint count changed")
    require(len(own) == 4, "Parent report checkpoint map changed")
    complete = _merge_hashes(inventory, own)
    require(len(complete) == report.get("checkpoint_count") == 224,
            "Complete parent checkpoint count changed")
    for name, expected in complete.items():
        _raw_auth(name, expected)
    protected = _hash_map(report.get("protected_file_hashes"),
                          "Parent protected files")
    dirty = _hash_map(report.get("legacy_dirty_path_hashes"),
                      "Parent dirty files")
    failed = _hash_map(
        report.get("preserved_failed_artifact_hashes_union"),
        "Parent failed-artifact union",
    )
    require(len(protected) == 16 and len(dirty) == 2 and len(failed) == 162,
            "Parent protected, dirty, or failed mapping changed")
    return {
        "report_raw": report_raw, "report": report,
        "report_view": _parent_report_view(report),
        "result_raw": result_raw, "result": result,
        "result_path": result_path, "result_sha256": result_sha,
        "source_pins": source_pins, "inventory": inventory, "own": own,
        "complete": complete, "protected": protected, "dirty": dirty,
        "failed": failed, "false_flags": list(false_flags),
    }


def _code_hashes():
    result = {}
    for name in (RUNNER, KERNEL):
        path = _resolve(name)
        require(path.exists(), "New code blob is missing: " + name)
        result[name] = sha_path(path)
    return result


def _input_spec(path):
    path = Path(path)
    require(path.is_absolute(), "--input-spec must be an absolute path")
    require(path.exists(), "Input spec is missing: " + str(path))
    raw = path.read_bytes()
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict)
            and value.get("schema") == INPUT_SPEC_SCHEMA
            and value.get("component") == COMPONENT,
            "Signed-history input-spec schema changed")
    return path, raw, sha_bytes(raw), value


def _input_pins(spec):
    additional = _hash_map(spec.get("additional_raw_pins", {}),
                           "Input additional raw pins")
    receipts = {}
    for row in spec.get("static_receipts", []):
        require(isinstance(row, dict) and isinstance(row.get("path"), str),
                "Malformed input static receipt")
        expected = row.get("sha256")
        require(_is_sha(expected), "Malformed input static receipt SHA")
        receipts[row["path"]] = expected
        _raw_auth(row["path"], expected)
    return _merge_hashes(additional, receipts)


def _build_packet(input_spec_path):
    parent = _parent_scope()
    inherited_true = {}
    for mapping in (parent["report"], parent["result"]):
        for key, value in mapping.items():
            if isinstance(value, bool) and value and key not in parent["false_flags"]:
                inherited_true[key] = True
    code = _code_hashes()
    spec_path, spec_raw, spec_sha, spec = _input_spec(input_spec_path)
    extra = _input_pins(spec)
    stage = _merge_hashes(
        {PARENT_REPORT: PARENT_REPORT_SHA256,
         parent["result_path"]: parent["result_sha256"],
         RUNNER: code[RUNNER], KERNEL: code[KERNEL],
         str(spec_path): spec_sha},
        extra,
    )
    source_union = _merge_hashes(parent["source_pins"], stage)
    for name, expected in source_union.items():
        _raw_auth(name, expected)
    digest = _digest({
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_sha256": parent["result_sha256"],
        "current_code_hashes": code,
        "input_spec_sha256": spec_sha,
        "stage_hashes": stage,
    })
    return {
        "schema": INPUT_SCHEMA, "component": COMPONENT, "metadata_only": True,
        "parent_report": parent["report_view"],
        "parent_report_is_authenticated_view": True,
        "parent_report_view_excluded_fields": parent["report_view"][
            "parent_report_view_excluded_fields"],
        "parent_result": parent["result"],
        "parent_report_metadata": {"path": PARENT_REPORT,
            "raw_sha256": PARENT_REPORT_SHA256, "raw_bytes": len(parent["report_raw"]),
            "schema": parent["report"].get("schema"),
            "component": parent["report"].get("component")},
        "parent_result_metadata": {"path": parent["result_path"],
            "raw_sha256": parent["result_sha256"], "raw_bytes": len(parent["result_raw"]),
            "schema": parent["result"].get("schema"),
            "component": parent["result"].get("component")},
        "input": spec,
        "input_spec": {"path": str(spec_path), "raw_sha256": spec_sha,
                        "raw_bytes": len(spec_raw)},
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
        "inherited_parent_true_flags": dict(inherited_true),
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "raw_only_inherited_sources": True,
        "ancestor_decode_forbidden": True,
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def prepare(input_spec_path):
    require(not any(path.exists() for path in CHECKPOINTS),
            "Signed-history checkpoint directory already contains an attempt")
    packet = _build_packet(input_spec_path)
    _write_new(SOURCE_CHECKPOINT, packet)
    return {"passed": True, "metadata_only": True, "component": COMPONENT,
            "source": SOURCE_CHECKPOINT.relative_to(ROOT).as_posix(),
            "source_sha256": sha_path(SOURCE_CHECKPOINT),
            "input_digest": packet["input_digest"],
            "source_union_count": packet["source_union_count"],
            "runtime_model_metadata": None}


def _verify_packet(packet, parent):
    require(packet.get("schema") == INPUT_SCHEMA
            and packet.get("component") == COMPONENT
            and packet.get("metadata_only") is True
            and packet.get("runtime_model_metadata") is None
            and packet.get("original_source_evaluations") == 0,
            "Frozen signed-history source packet scope changed")
    require(packet.get("parent_report_sha256") == PARENT_REPORT_SHA256
            and packet.get("parent_result_sha256") == parent["result_sha256"]
            and packet.get("parent_result_path") == parent["result_path"]
            and packet.get("parent_commit") == PARENT_COMMIT,
            "Frozen parent metadata changed")
    require(packet.get("parent_report") == parent["report_view"]
            and packet.get("parent_result") == parent["result"],
            "Frozen authenticated parent view/result changed")
    code = packet.get("current_code_hashes")
    require(code == packet.get("implementation_hashes")
            and isinstance(code, dict) and set(code) == {RUNNER, KERNEL}
            and code == _code_hashes(), "Frozen new code changed")
    spec_meta = packet.get("input_spec")
    require(isinstance(spec_meta, dict) and isinstance(spec_meta.get("path"), str)
            and Path(spec_meta["path"]).is_absolute()
            and _is_sha(spec_meta.get("raw_sha256")), "Frozen input spec changed")
    spec_raw = _raw_auth(spec_meta["path"], spec_meta["raw_sha256"])
    require(len(spec_raw) == spec_meta.get("raw_bytes"),
            "Frozen input-spec byte count changed")
    spec = json.loads(spec_raw.decode("utf-8"))
    require(spec == packet.get("input") and spec.get("schema") == INPUT_SPEC_SCHEMA,
            "Frozen input payload changed")
    extra = _input_pins(spec)
    expected_stage = _merge_hashes(
        {PARENT_REPORT: PARENT_REPORT_SHA256,
         parent["result_path"]: parent["result_sha256"],
         RUNNER: code[RUNNER], KERNEL: code[KERNEL],
         spec_meta["path"]: spec_meta["raw_sha256"]}, extra)
    require(packet.get("stage_hashes") == expected_stage,
            "Frozen stage manifest changed")
    expected_union = _merge_hashes(parent["source_pins"], expected_stage)
    require(packet.get("source_union_hashes") == expected_union
            and packet.get("source_union_count") == len(expected_union),
            "Frozen source union changed")
    for name, expected in expected_union.items():
        _raw_auth(name, expected)
    require(packet.get("parent_checkpoint_inventory") == parent["inventory"]
            and packet.get("checkpoint_hashes") == parent["own"]
            and packet.get("complete_parent_checkpoint_hashes") == parent["complete"]
            and packet.get("complete_parent_checkpoint_count") == 224,
            "Frozen checkpoint maps changed")
    require(packet.get("protected_file_hashes") == parent["protected"]
            and packet.get("legacy_dirty_path_hashes") == parent["dirty"],
            "Frozen protected or dirty maps changed")
    require(
        packet.get("preserved_failed_artifact_hashes_union") == parent["failed"]
        and packet.get("preserved_failed_artifact_count") == len(parent["failed"])
        and packet.get("inherited_false_flags") == parent["false_flags"],
        "Frozen failed-artifact union or inherited false flags changed",
    )
    inherited_true = {}
    for mapping in (parent["report"], parent["result"]):
        for key, value in mapping.items():
            if isinstance(value, bool) and value and key not in parent["false_flags"]:
                inherited_true[key] = True
    require(packet.get("inherited_parent_true_flags") == inherited_true,
            "Frozen inherited true flags changed")
    expected_digest = _digest({
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result_sha256": parent["result_sha256"],
        "current_code_hashes": code,
        "input_spec_sha256": spec_meta["raw_sha256"],
        "stage_hashes": expected_stage,
    })
    require(packet.get("input_digest") == expected_digest,
            "Frozen input digest changed")


def _load_packet():
    source_raw, packet = _read_json(SOURCE_CHECKPOINT)
    parent = _parent_scope()
    _verify_packet(packet, parent)
    return source_raw, packet, parent


def _brief_marker(packet, source_sha):
    require(BRIEF.exists(), "Resume brief is missing")
    brief_raw, brief = _read_json(BRIEF)
    marker = brief.get("in_progress", {}).get(COMPONENT)
    require(isinstance(marker, dict), "Resume brief has no signed-history marker")
    source_rel = SOURCE_CHECKPOINT.relative_to(ROOT).as_posix()
    require(marker.get("phase") == "prepared_unentered"
            and marker.get("source_checkpoint_sha256") == source_sha
            and marker.get("source_checkpoint_path") in (source_rel, str(SOURCE_CHECKPOINT))
            and marker.get("input_digest") == packet["input_digest"]
            and marker.get("kernel_sha256") == packet["current_code_hashes"][KERNEL]
            and marker.get("runner_sha256") == packet["current_code_hashes"][RUNNER],
            "Resume brief does not point to the prepared signed-history packet")
    return sha_bytes(brief_raw)


def _append_event(path, events, packet, stage, payload):
    event = {"schema": DIAGNOSTIC_SCHEMA, "component": COMPONENT,
             "index": len(events), "stage": stage,
             "input_digest": packet["input_digest"], "payload": _clone(payload)}
    with path.open("ab") as stream:
        stream.write(_encode(event)); stream.flush(); os.fsync(stream.fileno())
    events.append(event)


def _source_for_derive(packet, parent):
    inherited = {}
    for mapping in (parent["report"], parent["result"]):
        for key, value in mapping.items():
            if isinstance(value, bool) and value and key not in parent["false_flags"]:
                inherited[key] = value
    return {
        "schema": INPUT_SCHEMA, "component": COMPONENT,
        "input": _clone(packet["input"]), "input_spec": _clone(packet["input_spec"]),
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
        "preserved_failed_artifact_hashes_union": _clone(
            packet["preserved_failed_artifact_hashes_union"]
        ),
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


def _fallback(packet, error=None, calls=0, imported=False):
    result = {
        "schema": RESULT_SCHEMA, "component": COMPONENT, "passed": False,
        "metadata_only": False, "source_only": True,
        "input_digest": packet["input_digest"], "runtime_model_metadata": None,
        "original_source_evaluations": 0, "derive_entry_calls": calls,
        "derive_calls": calls, "entry_attempts": 1, "kernel_imported": imported,
        "source_packet_defined": False,
        "signed_history_source_gap_functions_defined": False,
        "native_N_scaling_applied_once_bound": False,
        "five_signed_density_defect_identities_proved": False,
        "paired_history_defect_transport_proved": False,
        "source_inspection_audit": {
            "original_source_evaluations": 0, "numeric_field_evaluations": 0,
            "old_source_consumer_calls": 0, "old_proof_functions_called": 0,
            "old_operator_calls": 0, "old_source_normalizer_calls": 0,
        },
    }
    if error is not None:
        result.update({"caught_exception": True, "error_type": type(error).__name__,
                      "error": str(error), "traceback": traceback.format_exc()})
    for name in packet.get("inherited_false_flags", []):
        result[name] = False
    for name, value in packet.get("inherited_parent_true_flags", {}).items():
        result[name] = value
    return result


def _normalise_result(value, packet, calls, imported):
    result = _clone(value) if isinstance(value, dict) else _fallback(
        packet, TypeError("derive returned a non-mapping result"), calls, imported)
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("passed", False)
    result.setdefault("metadata_only", False)
    result.setdefault("source_only", True)
    result.setdefault("input_digest", packet["input_digest"])
    result.setdefault("runtime_model_metadata", None)
    result.setdefault("original_source_evaluations", 0)
    for name, value in packet.get("inherited_parent_true_flags", {}).items():
        result.setdefault(name, value)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    result["kernel_imported"] = imported
    return result


def _validate_result(result, packet):
    require(isinstance(result, dict) and result.get("schema") == RESULT_SCHEMA
            and result.get("component") == COMPONENT
            and result.get("metadata_only") is False
            and isinstance(result.get("passed"), bool),
            "Signed-history result schema or scope changed")
    require(result.get("input_digest") == packet["input_digest"]
            and result.get("runtime_model_metadata") is None
            and result.get("original_source_evaluations") == 0
            and result.get("entry_attempts") == 1,
            "Signed-history execution scope changed")
    calls = result.get("derive_entry_calls")
    require(isinstance(calls, int) and not isinstance(calls, bool)
            and calls in (0, 1) and result.get("derive_calls") == calls
            and isinstance(result.get("kernel_imported"), bool),
            "Signed-history derive bookkeeping changed")
    if result.get("passed") is True:
        require(calls == 1 and result.get("kernel_imported") is True,
                "Passed signed-history result requires one successful derive call")
    for key in (
        "signed_history_source_gap_functions_defined",
        "native_N_scaling_applied_once_bound",
        "five_signed_density_defect_identities_proved",
        "paired_history_defect_transport_proved",
    ):
        if key in result:
            require(isinstance(result[key], bool),
                    "Signed-history narrow flag is not boolean: " + key)
    forbidden = (
        "actual_C3_C5_incoming_source_identity_proved",
        "actual_C3_C5_control_source_identity_proved",
        "actual_C3_C5_profile_identity_proved",
        "actual_C3_cap_to_C5_source_bound_proved",
        "compact_volume_cap_same_source_bound_proved",
        "global_five_moment_closure", "global_cap_proved",
        "numeric_values_computed", "full_goal_complete",
    )
    for key in forbidden:
        if key in result:
            require(result[key] is False, "Forbidden completion flag is true: " + key)
    for key in packet.get("inherited_false_flags", []):
        require(result.get(key) is False,
                "Inherited false flag became true or disappeared: " + key)
    for key, value in packet.get("inherited_parent_true_flags", {}).items():
        require(result.get(key) is value,
                "Inherited true flag changed: " + key)
    if result.get("passed") is True:
        audit = result.get("source_inspection_audit")
        require(isinstance(audit, dict),
                "Passed signed-history result lacks source inspection audit")
        for key in (
            "new_signed_history_source_rows",
            "new_signed_density_order_checks",
            "new_paired_transport_order_checks",
            "new_energy_pressure_local_correlation_checks",
        ):
            require(isinstance(audit.get(key), int)
                    and not isinstance(audit.get(key), bool)
                    and audit[key] >= 0,
                    "Passed signed-history result lacks audit key: " + key)
    audit = result.get("source_inspection_audit")
    if isinstance(audit, dict):
        for key in (
            "new_signed_history_source_rows",
            "new_signed_density_order_checks",
            "new_paired_transport_order_checks",
            "new_energy_pressure_local_correlation_checks",
        ):
            if key in audit:
                require(isinstance(audit[key], int) and not isinstance(audit[key], bool)
                        and audit[key] >= 0,
                        "Signed-history audit count is malformed: " + key)
        for key in ("original_source_evaluations", "numeric_field_evaluations",
                    "old_source_consumer_calls", "old_proof_functions_called",
                    "old_operator_calls", "old_source_normalizer_calls"):
            if key in audit and audit[key] is not None:
                require(audit[key] == 0, "Forbidden source call recorded")


def _load_saved(packet):
    require(ENTERED_CHECKPOINT.exists() and RESULT_CHECKPOINT.exists()
            and DIAGNOSTIC_CHECKPOINT.exists(), "Saved attempt is incomplete")
    entered_raw, entered = _read_json(ENTERED_CHECKPOINT)
    result_raw, result = _read_json(RESULT_CHECKPOINT)
    journal_raw = DIAGNOSTIC_CHECKPOINT.read_bytes()
    source_sha = sha_path(SOURCE_CHECKPOINT)
    require(entered.get("schema") == ENTERED_SCHEMA
            and entered.get("component") == COMPONENT
            and entered.get("source_checkpoint_sha256") == source_sha
            and entered.get("input_digest") == packet["input_digest"]
            and entered.get("kernel_imported_before_entry") is False
            and entered.get("derive_calls_before_entry") == 0
            and entered.get("entry_attempts") == 1,
            "Entered signed-history checkpoint changed")
    require(result.get("source_checkpoint_sha256") == source_sha
            and result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw)
            and result.get("diagnostic_checkpoint_sha256") == sha_bytes(journal_raw),
            "Saved result provenance changed")
    _validate_result(result, packet)
    return entered_raw, result_raw, result, journal_raw


def run():
    require(SOURCE_CHECKPOINT.exists()
            and not ENTERED_CHECKPOINT.exists() and not RESULT_CHECKPOINT.exists()
            and not DIAGNOSTIC_CHECKPOINT.exists(),
            "Frozen source is absent or attempt already entered")
    source_raw, packet, parent = _load_packet()
    brief_sha = _brief_marker(packet, sha_bytes(source_raw))
    entered = {
        "schema": ENTERED_SCHEMA, "component": COMPONENT,
        "source_checkpoint_sha256": sha_bytes(source_raw),
        "source_checkpoint_path": SOURCE_CHECKPOINT.relative_to(ROOT).as_posix(),
        "input_digest": packet["input_digest"],
        "current_code_hashes": dict(packet["current_code_hashes"]),
        "brief_sha256_before_entry": brief_sha,
        "kernel_imported_before_entry": False, "derive_calls_before_entry": 0,
        "entry_attempts": 1, "entered_code_immutable": True,
        "never_retry_entered_code": True, "runtime_model_metadata": None,
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write_new(ENTERED_CHECKPOINT, entered)
    _write_new(DIAGNOSTIC_CHECKPOINT, b"")
    events = []

    def journal_callback(event=None, *args, **fields):
        if isinstance(event, dict):
            payload = dict(event); payload.update(fields)
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
        source = _source_for_derive(packet, parent)
        spec = importlib.util.spec_from_file_location(
            "_original_compact_signed_history_defect_transport_entry",
            _resolve(KERNEL),
        )
        require(spec is not None and spec.loader is not None,
                "Cannot create new-kernel import specification")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        imported = True
        derive = getattr(module, "derive", None)
        require(callable(derive), "New kernel lacks derive(source, journal_callback)")
        calls = 1
        value = derive(source, journal_callback)
        result = _normalise_result(value, packet, calls, imported)
        _append_event(DIAGNOSTIC_CHECKPOINT, events, packet, "derive_returned", {
            "derive_entry_calls": calls, "passed": result.get("passed"),
            "runtime_model_metadata": None,
        })
    except Exception as error:
        result = _fallback(packet, error, calls, imported)
        _append_event(DIAGNOSTIC_CHECKPOINT, events, packet, "derive_exception", {
            "error_type": type(error).__name__, "error": str(error),
            "derive_entry_calls": calls,
        })
    result["source_checkpoint_sha256"] = sha_path(SOURCE_CHECKPOINT)
    result["entry_checkpoint_sha256"] = sha_path(ENTERED_CHECKPOINT)
    result["diagnostic_checkpoint_sha256"] = sha_path(DIAGNOSTIC_CHECKPOINT)
    result["diagnostic_event_count"] = len(events)
    result["diagnostic_summary"] = {"event_count": len(events),
                                    "journal_sha256": result["diagnostic_checkpoint_sha256"]}
    result["kernel_imported"] = imported
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["entry_attempts"] = 1
    result["runtime_model_metadata"] = None
    _write_new(RESULT_CHECKPOINT, result)
    _validate_result(result, packet)
    return {"passed": result.get("passed") is True, "metadata_only": True,
            "component": COMPONENT,
            "result": RESULT_CHECKPOINT.relative_to(ROOT).as_posix(),
            "result_sha256": sha_path(RESULT_CHECKPOINT),
            "derive_entry_calls": result.get("derive_entry_calls"),
            "runtime_model_metadata": None}


def inspect():
    source_raw, packet, _ = _load_packet()
    entered_raw, result_raw, result, journal_raw = _load_saved(packet)
    return {"passed": result.get("passed") is True, "metadata_only": True,
            "component": COMPONENT, "source_sha256": sha_bytes(source_raw),
            "entered_sha256": sha_bytes(entered_raw),
            "result_sha256": sha_bytes(result_raw),
            "diagnostics_sha256": sha_bytes(journal_raw),
            "derive_entry_calls": result.get("derive_entry_calls"),
            "runtime_model_metadata": None}


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
