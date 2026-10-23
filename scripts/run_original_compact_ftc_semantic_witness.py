"""Durable metadata runner for the compact FTC semantic-witness leaf.

This file is a checkpoint harness only.  It authenticates the saved parent
report/result, one producer source file, a caller-supplied audit receipt, and
the two new repository code blobs.  ``prepare`` writes one frozen source
packet, ``run`` enters at most once and calls only the new semantic-witness
module, and ``inspect`` reads saved checkpoints and result metadata.

The new module receives ``derive(source, journal_callback)`` exactly once.
``source`` is a mapping with ``records`` for the parent result and producer
source, the decoded audit receipt, and their authenticated metadata.  The
callback accepts one optional dictionary or event string plus keyword fields;
each event is appended and fsynced to ``diagnostics.jsonl``.

No old checker, normalizer, quadrature, provider, proof operator, source
consumer, or forbidden history artifact is imported or executed here.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import os
import sys
import traceback


ROOT = Path(
    "C:/Users/z5242/Documents/Codex/2026-09-13/"
    "before-doing-any-project-work-verify/work/st073-transition-next"
)
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "compact_ftc_semantic_witness"
RUNNER = "scripts/run_original_compact_ftc_semantic_witness.py"
KERNEL = "src/openai_ns_reconstruction/original_compact_ftc_semantic_witness.py"

PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_ftc_port_binding.json"
)
PARENT_REPORT_SHA256 = (
    "cc15d8dadfc3cd269f13a64f4ed48e5451be7ae4405af372ddcab95bb1ae2c70"
)
PARENT_RESULT = "tmp/original_compact_ftc_port_binding/result.json"
PARENT_RESULT_SHA256 = (
    "8b4fcfbf6bb44f84b033c415f303762f716e9eeebf813a342bcc1d5c643f24a8"
)
CURRENT_PARENT_COMMIT = "3eee4e118de2d07cadc8acc378a94c0f517c9265"
PRODUCER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py"
)
PRODUCER_SOURCE_SHA256 = (
    "661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e"
)

INPUT_NAMES = (PARENT_RESULT, PRODUCER_SOURCE)
CODE2_PATHS = (RUNNER, KERNEL)
RAW_ONLY_MARKERS = ("resume_first_switch", "rm_loop_density_identity")

INPUT_SCHEMA = "original_compact_ftc_semantic_witness_inputs_v1"
ENTERED_SCHEMA = "original_compact_ftc_semantic_witness_entered_v1"
RESULT_SCHEMA = "original_compact_ftc_semantic_witness"
DIAGNOSTIC_SCHEMA = "compact_ftc_semantic_witness_stage_diagnostic_v1"

CHECKPOINT_DIR = ROOT / "tmp/original_compact_ftc_semantic_witness"
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"
CHECKPOINT_NAMES = {
    "source": SOURCE_CHECKPOINT,
    "entered": ENTERED_CHECKPOINT,
    "result": RESULT_CHECKPOINT,
    "diagnostics": DIAGNOSTIC_CHECKPOINT,
}

EXPECTED_NEW_PROOF_COUNT = 5
EXPECTED_CACHED_WITNESS_COUNT = 54
EXPECTED_OLD_SOURCE_CALLS = 0
PREPARED_PHASES = ("prepared_not_entered", "prepared_unentered")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def sha_path(path):
    return sha_bytes(Path(path).read_bytes())


def _encode(value):
    return (json.dumps(value, ensure_ascii=True, sort_keys=True,
                       allow_nan=False).encode("utf-8") + b"\n")


def _clone(value):
    return json.loads(_encode(value).decode("utf-8"))


def _digest(value):
    return sha_bytes(json.dumps(value, ensure_ascii=True, sort_keys=True,
                                separators=(",", ":"), allow_nan=False)
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


def _stored_path(path):
    path = Path(path)
    if path.is_absolute():
        return str(path)
    return path.as_posix()


def _raw_auth(name, expected):
    path = _resolve(name)
    require(path.exists(), "Missing authenticated artifact: " + str(name))
    raw = path.read_bytes()
    require(sha_bytes(raw) == expected, "Changed authenticated artifact: " + str(name))
    return raw


def _raw_name_forbidden(name):
    normalized = str(name).replace("-", "_").lower()
    return any(marker in normalized for marker in RAW_ONLY_MARKERS)


def _hash_map(value, label, authenticate=True):
    require(isinstance(value, dict), label + " is not a mapping")
    result = {}
    for name, expected in value.items():
        require(isinstance(name, str) and isinstance(expected, str)
                and len(expected) == 64,
                label + " contains a malformed SHA-256")
        if authenticate:
            _raw_auth(name, expected)
        result[name] = expected
    return result


def _flatten_hash_maps(value, label):
    """Flatten mapping keys ending in ``hashes`` and reject collisions."""
    flattened = {}

    def visit(node, path):
        if isinstance(node, list):
            for index, child in enumerate(node):
                if isinstance(child, (dict, list)):
                    visit(child, path + "[" + str(index) + "]")
            return
        require(isinstance(node, dict), path + " is not a mapping")
        for key, child in node.items():
            child_path = path + "." + str(key)
            if isinstance(key, str) and key.endswith("hashes"):
                require(isinstance(child, dict), child_path + " is not a hash map")
                for name, expected in child.items():
                    require(isinstance(name, str) and isinstance(expected, str)
                            and len(expected) == 64,
                            child_path + " contains a malformed SHA-256")
                    if name in flattened:
                        require(flattened[name] == expected,
                                "Conflicting hash collision at " + name)
                    else:
                        flattened[name] = expected
            elif isinstance(child, (dict, list)):
                visit(child, child_path)

    visit(value, label)
    return flattened


def _parent_scope():
    report_raw = _raw_auth(PARENT_REPORT, PARENT_REPORT_SHA256)
    result_raw = _raw_auth(PARENT_RESULT, PARENT_RESULT_SHA256)
    report = json.loads(report_raw.decode("utf-8"))
    parent_result = json.loads(result_raw.decode("utf-8"))
    require(report.get("schema") == "original_compact_ftc_port_binding_stage_v1"
            and report.get("component") == "compact_ftc_port_binding"
            and report.get("source_only") is True
            and report.get("full_goal_complete") is False
            and report.get("metadata_only") is True
            and report.get("runtime_model_metadata") is None,
            "Parent report scope changed")
    require(report.get("full_result_raw_sha256") == PARENT_RESULT_SHA256
            and report.get("full_result_raw_bytes") == len(result_raw)
            and report.get("full_result") == parent_result,
            "Parent result does not equal parent report.full_result")
    require(parent_result.get("schema") == "original_compact_ftc_port_binding"
            and parent_result.get("component") == "compact_ftc_port_binding"
            and parent_result.get("runtime_model_metadata") is None,
            "Parent result shape changed")
    source_pins = _hash_map(report.get("source_pins"), "Parent source pins")
    require(report.get("source_pin_count") == len(source_pins),
            "Parent source-pin count changed")
    inventory = report.get("parent_checkpoint_inventory")
    require(isinstance(inventory, dict), "Parent checkpoint inventory missing")
    inventory_hashes = _flatten_hash_maps(inventory, "parent_checkpoint_inventory")
    own_hashes = _flatten_hash_maps(
        {"checkpoint_hashes": report.get("checkpoint_hashes")},
        "parent_report_checkpoint_hashes",
    )
    parent_hashes = _hash_map(
        report.get("parent_checkpoint_hashes"),
        "Parent immediate checkpoint hashes",
    )
    require(len(inventory_hashes) == 168
            and len(parent_hashes) == 4
            and len(own_hashes) == 4,
            "Parent checkpoint map sizes changed")
    complete_hashes = {}
    for mapping in (inventory_hashes, parent_hashes, own_hashes):
        for name, expected in mapping.items():
            if name in complete_hashes:
                require(complete_hashes[name] == expected,
                        "Conflicting parent checkpoint hash collision at " + name)
            else:
                complete_hashes[name] = expected
    require(len(complete_hashes) == report.get("checkpoint_count") == 176,
            "Complete parent checkpoint count changed")
    for name, expected in complete_hashes.items():
        _raw_auth(name, expected)
    require(inventory.get("raw_sha_only") is True
            and inventory.get("history_decode_forbidden") is True,
            "Parent raw-only checkpoint policy changed")
    failed_union = _hash_map(
        report.get("preserved_failed_artifact_hashes_union"),
        "Parent failed-artifact union",
    )
    protected = _hash_map(report.get("protected_file_hashes"),
                          "Parent protected files")
    dirty = _hash_map(report.get("legacy_dirty_path_hashes"),
                      "Parent dirty files")
    require(len(protected) == 16 and len(dirty) == 2,
            "Parent protected/dirty mapping changed")
    retained = report.get("retained_construction_data")
    require(isinstance(retained, dict)
            and retained.get("all_old_error_data_retained") is True,
            "Parent retained construction data changed")
    return {
        "report_raw": report_raw,
        "report": report,
        "result_raw": result_raw,
        "result": parent_result,
        "source_pins": source_pins,
        "inventory": _clone(inventory),
        "inventory_hashes": inventory_hashes,
        "parent_hashes": parent_hashes,
        "own_hashes": own_hashes,
        "complete_hashes": complete_hashes,
        "failed_union": failed_union,
        "protected": protected,
        "dirty": dirty,
    }


def _code_hashes():
    result = {}
    for name in CODE2_PATHS:
        path = _resolve(name)
        require(path.exists(), "New code blob is missing: " + name)
        result[name] = sha_path(path)
    return result


def _decode_input(path, expected, kind):
    raw = _raw_auth(path, expected)
    decoded = raw.decode("utf-8")
    value = json.loads(decoded) if kind == "json" else decoded
    return {
        "path": _stored_path(path),
        "raw_sha256": sha_bytes(raw),
        "decoded_sha256": sha_bytes(decoded.encode("utf-8")),
        "raw_authenticated_before_decode": True,
        "encoding": kind,
        "source_text": decoded,
        "value": value,
    }


def _audit_path(path):
    candidate = Path(path).expanduser()
    require(candidate.is_absolute(), "--audit must be an absolute path")
    require(candidate.exists(), "Audit receipt is missing")
    return candidate


def _build_packet(audit_path):
    parent = _parent_scope()
    producer_raw = _raw_auth(PRODUCER_SOURCE, PRODUCER_SOURCE_SHA256)
    audit = _audit_path(audit_path)
    audit_raw = audit.read_bytes()
    audit_sha = sha_bytes(audit_raw)
    audit_value = json.loads(audit_raw.decode("utf-8"))
    require(isinstance(audit_value, dict), "Audit receipt is not a JSON mapping")
    code = _code_hashes()
    stage_hashes = {
        PARENT_REPORT: PARENT_REPORT_SHA256,
        PARENT_RESULT: PARENT_RESULT_SHA256,
        PRODUCER_SOURCE: PRODUCER_SOURCE_SHA256,
        str(audit): audit_sha,
        **code,
    }
    source_union = dict(parent["source_pins"])
    for name, expected in stage_hashes.items():
        previous = source_union.get(name)
        require(previous in (None, expected),
                "Stage artifact collides with an inherited source pin: " + name)
        source_union[name] = expected
    for name, expected in source_union.items():
        _raw_auth(name, expected)
    for name, expected in parent["protected"].items():
        _raw_auth(name, expected)
    for name, expected in parent["dirty"].items():
        _raw_auth(name, expected)
    input_meta = {
        PARENT_RESULT: {
            "path": PARENT_RESULT,
            "raw_sha256": PARENT_RESULT_SHA256,
            "encoding": "json",
        },
        PRODUCER_SOURCE: {
            "path": PRODUCER_SOURCE,
            "raw_sha256": PRODUCER_SOURCE_SHA256,
            "encoding": "utf-8",
        },
    }
    audit_meta = {
        "path": str(audit),
        "raw_sha256": audit_sha,
        "decoded_sha256": sha_bytes(audit_raw),
        "encoding": "json",
        "raw_authenticated_before_decode": True,
    }
    input_digest = _digest({
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "inputs": input_meta,
        "audit": audit_meta,
    })
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_only": True,
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result": PARENT_RESULT,
        "parent_result_sha256": PARENT_RESULT_SHA256,
        "producer_source": PRODUCER_SOURCE,
        "producer_source_sha256": PRODUCER_SOURCE_SHA256,
        "parent_commit": CURRENT_PARENT_COMMIT,
        "input_names": list(INPUT_NAMES),
        "input_count": len(INPUT_NAMES),
        "inputs": input_meta,
        "audit": audit_meta,
        "audit_value_schema": audit_value.get("schema"),
        "input_digest": input_digest,
        "inherited_source_pins": dict(parent["source_pins"]),
        "inherited_source_pin_count": len(parent["source_pins"]),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "stage_hashes": stage_hashes,
        "current_code_hashes": code,
        "implementation_hashes": code,
        "parent_checkpoint_inventory": parent["inventory"],
        "parent_checkpoint_hashes": parent["parent_hashes"],
        "report_checkpoint_hashes": parent["own_hashes"],
        "complete_parent_checkpoint_hashes": parent["complete_hashes"],
        "parent_checkpoint_inventory_count": len(parent["inventory_hashes"]),
        "parent_checkpoint_hashes_count": len(parent["parent_hashes"]),
        "report_checkpoint_hashes_count": len(parent["own_hashes"]),
        "complete_parent_checkpoint_count": len(parent["complete_hashes"]),
        "parent_inherited_checkpoint_count": parent["report"].get(
            "inherited_checkpoint_count"),
        "parent_checkpoint_count": parent["report"].get("checkpoint_count"),
        "preserved_failed_artifact_hashes_union": dict(parent["failed_union"]),
        "preserved_failed_artifact_count": len(parent["failed_union"]),
        "new_checkpoint_count": 4,
        "final_checkpoint_count": (parent["report"].get("checkpoint_count", 0)
                                    + 4),
        "protected_file_hashes": dict(parent["protected"]),
        "legacy_dirty_path_hashes": dict(parent["dirty"]),
        "retained_construction_data": _clone(
            parent["report"].get("retained_construction_data", {})),
        "retained_integral_term_inventory": _clone(
            parent["report"].get("retained_integral_term_inventory", {})),
        "expected_new_proof_count": EXPECTED_NEW_PROOF_COUNT,
        "expected_cached_witness_count": EXPECTED_CACHED_WITNESS_COUNT,
        "expected_old_source_consumer_calls": EXPECTED_OLD_SOURCE_CALLS,
        "raw_only_inherited_sources": True,
        "history_decode_forbidden": True,
        "runtime_model_metadata": None,
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def _paths():
    return dict(CHECKPOINT_NAMES)


def prepare(audit_path):
    paths = _paths()
    require(not any(path.exists() for path in paths.values()),
            "Semantic-witness checkpoint directory already contains an attempt")
    packet = _build_packet(audit_path)
    _write_new(paths["source"], packet)
    return {
        "passed": True,
        "metadata_only": True,
        "component": COMPONENT,
        "source": paths["source"].relative_to(ROOT).as_posix(),
        "source_sha256": sha_path(paths["source"]),
        "input_digest": packet["input_digest"],
        "input_count": packet["input_count"],
        "source_union_count": packet["source_union_count"],
        "runtime_model_metadata": None,
    }


def _load_packet(decode_inputs=False):
    source_raw, packet = _read_json(SOURCE_CHECKPOINT)
    require(packet.get("schema") == INPUT_SCHEMA
            and packet.get("component") == COMPONENT
            and packet.get("metadata_only") is True
            and packet.get("input_names") == list(INPUT_NAMES)
            and packet.get("input_count") == len(INPUT_NAMES)
            and packet.get("runtime_model_metadata") is None,
            "Frozen semantic-witness source packet changed")
    code = packet.get("current_code_hashes")
    require(code == packet.get("implementation_hashes")
            and isinstance(code, dict)
            and set(code) == set(CODE2_PATHS),
            "Frozen semantic-witness code map changed")
    for name, expected in code.items():
        require(sha_path(_resolve(name)) == expected,
                "New code changed after prepare: " + name)
    source_union = packet.get("source_union_hashes")
    require(isinstance(source_union, dict)
            and packet.get("source_union_count") == len(source_union),
            "Frozen source union changed")
    for name, expected in source_union.items():
        require(not _raw_name_forbidden(name) or _resolve(name).exists(),
                "Forbidden history artifact is unavailable: " + name)
        _raw_auth(name, expected)
    for name, expected in packet.get("protected_file_hashes", {}).items():
        _raw_auth(name, expected)
    for name, expected in packet.get("legacy_dirty_path_hashes", {}).items():
        _raw_auth(name, expected)
    parent = _parent_scope()
    require(packet.get("parent_report_sha256") == PARENT_REPORT_SHA256
            and packet.get("parent_result_sha256") == PARENT_RESULT_SHA256
            and packet.get("producer_source_sha256") == PRODUCER_SOURCE_SHA256
            and packet.get("parent_commit") == CURRENT_PARENT_COMMIT,
            "Frozen parent/source hashes changed")
    require(packet.get("inherited_source_pins") == parent["source_pins"]
            and packet.get("parent_checkpoint_inventory") == parent["inventory"]
            and packet.get("parent_checkpoint_hashes") == parent["parent_hashes"]
            and packet.get("report_checkpoint_hashes") == parent["own_hashes"]
            and packet.get("complete_parent_checkpoint_hashes")
            == parent["complete_hashes"]
            and packet.get("complete_parent_checkpoint_count") == 176
            and packet.get("protected_file_hashes") == parent["protected"]
            and packet.get("legacy_dirty_path_hashes") == parent["dirty"],
            "Frozen inherited metadata maps changed")
    records = None
    audit = None
    if decode_inputs:
        result_record = _decode_input(
            PARENT_RESULT, PARENT_RESULT_SHA256, "json")
        producer_record = _decode_input(
            PRODUCER_SOURCE, PRODUCER_SOURCE_SHA256, "utf-8")
        audit_path = packet["audit"]["path"]
        audit = _decode_input(audit_path, packet["audit"]["raw_sha256"], "json")
        require(result_record["value"] == parent["result"],
                "Parent result changed after packet preparation")
        require(packet["inputs"][PARENT_RESULT]["raw_sha256"]
                == result_record["raw_sha256"]
                and packet["inputs"][PRODUCER_SOURCE]["raw_sha256"]
                == producer_record["raw_sha256"]
                and packet["audit"]["raw_sha256"] == audit["raw_sha256"],
                "Selected input bytes changed after preparation")
        records = {
            PARENT_RESULT: result_record,
            PRODUCER_SOURCE: producer_record,
        }
    return source_raw, packet, records, audit


def _brief_prepared_source(packet, source_sha):
    require(BRIEF.exists(), "Resume brief is missing")
    brief_raw, brief = _read_json(BRIEF)
    node = brief.get("in_progress", {}).get(COMPONENT)
    require(isinstance(node, dict),
            "Resume brief has no semantic-witness prepared marker")
    source_rel = SOURCE_CHECKPOINT.relative_to(ROOT).as_posix()
    require(node.get("phase") in PREPARED_PHASES
            and node.get("state", "prepared_not_entered")
            == "prepared_not_entered"
            and node.get("source_checkpoint_sha256") == source_sha
            and node.get("source_checkpoint_path") in (source_rel,
                                                          str(SOURCE_CHECKPOINT))
            and node.get("kernel_sha256") == packet["current_code_hashes"][KERNEL]
            and node.get("runner_sha256") == packet["current_code_hashes"][RUNNER],
            "Resume brief does not point to the prepared semantic-witness source")
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


def _source_for_derive(packet, records, audit):
    metadata = {
        name: {key: value for key, value in record.items()
               if key not in ("value", "source_text")}
        for name, record in records.items()
    }
    values = {
        PARENT_RESULT: records[PARENT_RESULT]["value"],
        PRODUCER_SOURCE: records[PRODUCER_SOURCE]["source_text"],
    }
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "input_digest": packet["input_digest"],
        "records": values,
        "record_metadata": metadata,
        "audit": audit["value"],
        "audit_metadata": {
            key: value for key, value in audit.items()
            if key not in ("value", "source_text")
        },
        "parent_report_sha256": packet["parent_report_sha256"],
        "parent_result_sha256": packet["parent_result_sha256"],
        "producer_source_sha256": packet["producer_source_sha256"],
        "expected_new_proof_count": packet["expected_new_proof_count"],
        "expected_cached_witness_count": packet["expected_cached_witness_count"],
        "expected_old_source_consumer_calls": (
            packet["expected_old_source_consumer_calls"]
        ),
    }


def _fallback(packet, error=None, calls=0, import_succeeded=False):
    result = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "input_digest": packet["input_digest"],
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "cached_formal_proof_calls": 0,
        "derive_entry_calls": calls,
        "derive_calls": calls,
        "entry_attempts": 1,
        "kernel_imported": import_succeeded,
        "source_functions_called": 0,
        "fresh_source_proof_normalizations": 0,
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "old_source_consumer_calls": EXPECTED_OLD_SOURCE_CALLS,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "old_source_normalizer_calls": 0,
            "cached_formal_proof_calls": 0,
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


def _normalise_result(value, packet, calls, import_succeeded):
    if isinstance(value, dict):
        result = _clone(value)
    else:
        result = _fallback(packet,
                           TypeError("derive returned a non-mapping result"),
                           calls=calls, import_succeeded=import_succeeded)
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("passed", False)
    result.setdefault("source_only", True)
    result.setdefault("metadata_only", False)
    result.setdefault("input_digest", packet["input_digest"])
    result.setdefault("runtime_model_metadata", None)
    result.setdefault("original_source_evaluations", 0)
    result.setdefault("cached_formal_proof_calls", 0)
    result.setdefault("source_functions_called", 0)
    result.setdefault("entry_attempts", 1)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["kernel_imported"] = import_succeeded
    return result


def _candidate_count(result, names):
    for name in names:
        value = result.get(name)
        if value is not None:
            return value
    return None


def _validate_result(result, packet):
    require(isinstance(result, dict)
            and result.get("schema") == RESULT_SCHEMA
            and result.get("component") == COMPONENT
            and result.get("metadata_only") is False,
            "Semantic-witness result schema or metadata scope changed")
    require(result.get("input_digest") == packet["input_digest"]
            and result.get("runtime_model_metadata") is None
            and result.get("original_source_evaluations") == 0
            and result.get("cached_formal_proof_calls") == 0
            and result.get("entry_attempts") == 1,
            "Semantic-witness execution scope changed")
    calls = result.get("derive_entry_calls")
    require(isinstance(calls, int) and not isinstance(calls, bool)
            and calls in (0, 1)
            and result.get("derive_calls") == calls
            and isinstance(result.get("kernel_imported"), bool),
            "Semantic-witness derive bookkeeping changed")
    if result.get("passed") is True:
        require(calls == 1 and result.get("kernel_imported") is True,
                "A passed result requires one successful derive call")
    audit = result.get("source_inspection_audit")
    if isinstance(audit, dict):
        require(audit.get("old_source_consumer_calls",
                          EXPECTED_OLD_SOURCE_CALLS)
                == EXPECTED_OLD_SOURCE_CALLS
                and audit.get("old_proof_functions_called", 0) == 0
                and audit.get("old_operator_calls", 0) == 0
                and audit.get("old_source_normalizer_calls", 0) == 0,
                "Old source consumer/operator calls were recorded")
    fresh = _candidate_count(
        result, ("new_proof_count", "fresh_proof_count",
                 "fresh_source_proof_count", "witness_proof_count"))
    if fresh is not None:
        require(fresh == packet["expected_new_proof_count"],
                "New semantic-witness proof count changed")
    cached = _candidate_count(
        result, ("cached_witness_count", "cached_proof_count",
                 "inherited_witness_count"))
    if cached is not None:
        require(cached == packet["expected_cached_witness_count"],
                "Cached semantic-witness count changed")


def _load_saved(packet):
    paths = _paths()
    require(paths["entered"].exists() and paths["result"].exists()
            and paths["diagnostics"].exists(),
            "Saved semantic-witness attempt is incomplete")
    entered_raw, entered = _read_json(paths["entered"])
    result_raw, result = _read_json(paths["result"])
    journal_raw = paths["diagnostics"].read_bytes()
    source_sha = sha_path(paths["source"])
    require(entered.get("schema") == ENTERED_SCHEMA
            and entered.get("component") == COMPONENT
            and entered.get("source_checkpoint_sha256") == source_sha
            and entered.get("kernel_imported_before_entry") is False
            and entered.get("derive_calls_before_entry") == 0
            and entered.get("entry_attempts") == 1,
            "Entered semantic-witness checkpoint changed")
    require(result.get("source_checkpoint_sha256") == source_sha
            and result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw)
            and result.get("diagnostic_checkpoint_sha256")
            == sha_bytes(journal_raw),
            "Saved semantic-witness result provenance changed")
    _validate_result(result, packet)
    return {
        "entered_raw": entered_raw,
        "entered": entered,
        "result_raw": result_raw,
        "result": result,
        "journal_raw": journal_raw,
    }


def run():
    paths = _paths()
    require(paths["source"].exists()
            and not paths["entered"].exists()
            and not paths["result"].exists()
            and not paths["diagnostics"].exists(),
            "Frozen source is absent or semantic-witness attempt already entered")
    source_raw, packet, records, audit = _load_packet(decode_inputs=True)
    brief_sha = _brief_prepared_source(packet, sha_bytes(source_raw))
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
    _write_new(paths["entered"], entered)
    _write_new(paths["diagnostics"], b"")
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
        _append_event(paths["diagnostics"], events, packet,
                      "kernel_callback", payload)

    calls = 0
    import_succeeded = False
    try:
        source = _source_for_derive(packet, records, audit)
        spec = importlib.util.spec_from_file_location(
            "_original_compact_ftc_semantic_witness_entry", _resolve(KERNEL)
        )
        require(spec is not None and spec.loader is not None,
                "Cannot create new semantic-witness import specification")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        import_succeeded = True
        derive = getattr(module, "derive", None)
        require(callable(derive),
                "New semantic-witness module lacks derive(source, journal_callback)")
        calls = 1
        value = derive(source, journal_callback)
        result = _normalise_result(value, packet, calls, import_succeeded)
        _append_event(paths["diagnostics"], events, packet, "derive_returned", {
            "derive_entry_calls": calls,
            "passed": result.get("passed"),
            "runtime_model_metadata": None,
        })
    except Exception as error:
        result = _fallback(packet, error, calls=calls,
                           import_succeeded=import_succeeded)
        _append_event(paths["diagnostics"], events, packet, "derive_exception", {
            "error_type": type(error).__name__,
            "error": str(error),
            "derive_entry_calls": calls,
        })
    result["source_checkpoint_sha256"] = sha_path(paths["source"])
    result["entry_checkpoint_sha256"] = sha_path(paths["entered"])
    result["diagnostic_checkpoint_sha256"] = sha_path(paths["diagnostics"])
    result["diagnostic_event_count"] = len(events)
    result["kernel_imported"] = import_succeeded
    result["entry_attempts"] = 1
    result["runtime_model_metadata"] = None
    _write_new(paths["result"], result)
    # The complete result and journal are durable before this metadata check.
    _validate_result(result, packet)
    return {
        "passed": result.get("passed") is True,
        "metadata_only": True,
        "component": COMPONENT,
        "result": paths["result"].relative_to(ROOT).as_posix(),
        "result_sha256": sha_path(paths["result"]),
        "attempt_status": "passed" if result.get("passed") else "failed",
        "derive_entry_calls": result.get("derive_entry_calls"),
        "cached_formal_proof_calls": result.get("cached_formal_proof_calls"),
        "runtime_model_metadata": None,
    }


def inspect():
    source_raw, packet, _, _ = _load_packet(decode_inputs=False)
    saved = _load_saved(packet)
    return {
        "passed": saved["result"].get("passed") is True,
        "metadata_only": True,
        "component": COMPONENT,
        "source_sha256": sha_bytes(source_raw),
        "entered_sha256": sha_bytes(saved["entered_raw"]),
        "result_sha256": sha_bytes(saved["result_raw"]),
        "diagnostics_sha256": sha_bytes(saved["journal_raw"]),
        "attempt_status": ("passed" if saved["result"].get("passed")
                            else "failed"),
        "derive_entry_calls": saved["result"].get("derive_entry_calls"),
        "cached_formal_proof_calls": saved["result"].get(
            "cached_formal_proof_calls"),
        "runtime_model_metadata": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "inspect"))
    parser.add_argument("--audit", help="Absolute JSON audit receipt for prepare")
    args = parser.parse_args()
    if args.action == "prepare":
        require(args.audit is not None, "prepare requires --audit ABSOLUTE_PATH")
        value = prepare(args.audit)
    else:
        require(args.audit is None, "--audit is valid only with prepare")
        value = globals()[args.action]()
    print(json.dumps(value, ensure_ascii=True, allow_nan=False))


if __name__ == "__main__":
    main()
