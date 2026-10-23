"""Durable metadata runner for the compact volume-error source packet leaf.

This harness freezes a caller-supplied input manifest and the authenticated
parent semantic-witness packet.  ``prepare`` writes one source checkpoint
without decoding selected source records.  ``run`` requires the external
prepared brief marker, enters exactly once, imports only the new volume-error
module, and calls ``derive(source, journal_callback)`` at most once.  The
complete result and journal are fsynced before result metadata is inspected.

The runner contains no source consumer, checker, normalizer, quadrature,
provider, recurrence, or proof implementation.  It never evaluates the old
source functions and treats this source packet as definition/provenance data,
not as quantitative error reduction.
"""

from __future__ import annotations

from datetime import datetime, timezone
import argparse
import gzip
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

COMPONENT = "compact_volume_error_source_packet"
RUNNER = "scripts/run_original_compact_volume_error_source_packet.py"
KERNEL = "src/openai_ns_reconstruction/original_compact_volume_error_source_packet.py"

PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_ftc_semantic_witness.json"
)
PARENT_REPORT_SHA256 = (
    "e03410d1d4575377402ce3887b863a0610694b8917346090ffcc7d3b946166c9"
)
PARENT_RESULT = "tmp/original_compact_ftc_semantic_witness/result.json"
PARENT_RESULT_SHA256 = (
    "9a574b9c4e353174c9f934fb7bf23fd1570aac27f7bf24ba07352b45d99793f1"
)
# The report records the prior parent commit.  The source packet is prepared
# from the current parent checkout supplied for this leaf.
CURRENT_PARENT_COMMIT = "6685ba640edb9d7e72c2466712eb157609ec5ae5"

INPUT_SCHEMA = "original_compact_volume_error_source_packet_inputs_v1"
ENTERED_SCHEMA = "original_compact_volume_error_source_packet_entered_v1"
RESULT_SCHEMA = "original_compact_volume_error_source_packet"
DIAGNOSTIC_SCHEMA = "compact_volume_error_source_packet_stage_diagnostic_v1"

CHECKPOINT_DIR = ROOT / "tmp/original_compact_volume_error_source_packet"
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

ALLOWED_ENCODINGS = ("json", "utf8", "gzip_json")
PREPARED_PHASES = ("prepared_unentered", "prepared_not_entered")
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
            separators=(",", ":"),
            allow_nan=False,
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


def _stored_path(path):
    path = Path(path)
    if path.is_absolute():
        return str(path)
    return path.as_posix()


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


def _flatten_checkpoint_hashes(value, label):
    """Accept either a flat path map or nested ``*hashes`` maps.

    The current semantic-witness parent stores its 176 checkpoint hashes as a
    flat path map.  A small recursive fallback keeps this harness compatible
    with the nested shape used by earlier parent reports, while merging an
    identical path/hash pair and rejecting a conflicting pair.
    """

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


def _merge_hash_maps(*mappings):
    merged = {}
    for mapping in mappings:
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
    parent_result = json.loads(result_raw.decode("utf-8"))
    require(
        report.get("schema")
        == "original_compact_ftc_semantic_witness_publication_v1"
        and report.get("component") == "compact_ftc_semantic_witness"
        and report.get("runtime_model_metadata") is None,
        "Parent semantic-witness report scope changed",
    )
    require(
        report.get("full_result_raw_sha256") == PARENT_RESULT_SHA256
        and report.get("full_result_raw_bytes") == len(result_raw)
        and report.get("full_result") == parent_result,
        "Parent result does not equal parent report.full_result",
    )
    require(
        parent_result.get("schema") == "original_compact_ftc_semantic_witness"
        and parent_result.get("component") == "compact_ftc_semantic_witness"
        and parent_result.get("runtime_model_metadata") is None,
        "Parent result shape changed",
    )

    source_pins = _hash_map(report.get("source_pins"), "Parent source pins")
    require(report.get("source_pin_count") == len(source_pins) == 1897,
            "Parent source-pin count changed")

    inventory = report.get("parent_checkpoint_inventory")
    require(isinstance(inventory, dict), "Parent checkpoint inventory missing")
    inventory_hashes = _flatten_checkpoint_hashes(
        inventory, "parent_checkpoint_inventory"
    )
    own_hashes = _flatten_checkpoint_hashes(
        {"checkpoint_hashes": report.get("checkpoint_hashes")},
        "parent_report_checkpoint_hashes",
    )
    immediate_value = report.get("parent_checkpoint_hashes")
    immediate_hashes = (
        {}
        if immediate_value is None
        else _hash_map(immediate_value, "Parent immediate checkpoint hashes")
    )
    require(len(inventory_hashes) == 176 and len(own_hashes) == 4,
            "Parent checkpoint map sizes changed")
    complete_hashes = _merge_hash_maps(
        inventory_hashes, immediate_hashes, own_hashes
    )
    require(len(complete_hashes) == report.get("checkpoint_count") == 180,
            "Complete parent checkpoint count changed")
    for name, expected in complete_hashes.items():
        _raw_auth(name, expected)

    failed_union = _hash_map(
        report.get("preserved_failed_artifact_hashes_union"),
        "Parent failed-artifact union",
    )
    protected = _hash_map(
        report.get("protected_file_hashes"), "Parent protected files"
    )
    dirty = _hash_map(
        report.get("legacy_dirty_path_hashes"), "Parent dirty files"
    )
    require(len(protected) == 16 and len(dirty) == 2,
            "Parent protected/dirty mapping changed")
    return {
        "report_raw": report_raw,
        "report": report,
        "result_raw": result_raw,
        "result": parent_result,
        "source_pins": source_pins,
        "inventory": _clone(inventory),
        "inventory_hashes": inventory_hashes,
        "parent_hashes": immediate_hashes,
        "own_hashes": own_hashes,
        "complete_hashes": complete_hashes,
        "failed_union": failed_union,
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


def _load_input_spec(path):
    candidate = Path(path).expanduser()
    require(candidate.is_absolute(), "--inputs must be an absolute path")
    require(candidate.exists(), "Input specification is missing")
    raw = candidate.read_bytes()
    raw_sha = sha_bytes(raw)
    value = json.loads(raw.decode("utf-8"))
    require(isinstance(value, dict) and value,
            "Input specification must be a nonempty JSON mapping")
    records = {}
    for record_path, entry in value.items():
        require(isinstance(record_path, str) and record_path,
                "Input specification contains a malformed path")
        require(isinstance(entry, dict),
                "Input specification entry is not a mapping: " + record_path)
        expected = entry.get("raw_sha256")
        encoding = entry.get("encoding")
        decoded_expected = entry.get("decoded_sha256")
        require(_is_sha(expected),
                "Input specification has malformed raw_sha256: " + record_path)
        require(encoding in ALLOWED_ENCODINGS,
                "Unsupported input encoding: " + str(encoding))
        if decoded_expected is not None:
            require(_is_sha(decoded_expected),
                    "Malformed decoded_sha256: " + record_path)
        if encoding == "gzip_json":
            require(decoded_expected is not None,
                    "gzip_json requires decoded_sha256: " + record_path)
        _raw_auth(record_path, expected)
        records[record_path] = {
            "path": record_path,
            "raw_sha256": expected,
            "encoding": encoding,
            "decoded_sha256": decoded_expected,
            "raw_authenticated_before_decode": True,
        }
    return {
        "path": str(candidate),
        "raw_sha256": raw_sha,
        "records": records,
        "raw_bytes": len(raw),
    }


def _build_packet(inputs_path):
    parent = _parent_scope()
    spec = _load_input_spec(inputs_path)
    code = _code_hashes()

    stage_hashes = _merge_hash_maps(
        {
            PARENT_REPORT: PARENT_REPORT_SHA256,
            PARENT_RESULT: PARENT_RESULT_SHA256,
            spec["path"]: spec["raw_sha256"],
            **code,
        },
        {
            name: metadata["raw_sha256"]
            for name, metadata in spec["records"].items()
        },
    )
    source_union = _merge_hash_maps(parent["source_pins"], stage_hashes)
    for name, expected in source_union.items():
        _raw_auth(name, expected)
    for name, expected in parent["protected"].items():
        _raw_auth(name, expected)
    for name, expected in parent["dirty"].items():
        _raw_auth(name, expected)

    input_digest = _digest(
        {
            "parent_report_sha256": PARENT_REPORT_SHA256,
            "parent_result_sha256": PARENT_RESULT_SHA256,
            "input_spec_path": spec["path"],
            "input_spec_sha256": spec["raw_sha256"],
            "records": spec["records"],
            "current_code_hashes": code,
        }
    )
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_packet_definition_only": True,
        "quantitative_packet_error_reduction_proved": False,
        "source_only": True,
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "parent_result": PARENT_RESULT,
        "parent_result_sha256": PARENT_RESULT_SHA256,
        "parent_commit": CURRENT_PARENT_COMMIT,
        "input_spec_path": spec["path"],
        "input_spec_raw_sha256": spec["raw_sha256"],
        "input_spec_raw_bytes": spec["raw_bytes"],
        "input_spec": _clone(spec["records"]),
        "input_count": len(spec["records"]),
        "input_digest": input_digest,
        "inherited_source_pins": dict(parent["source_pins"]),
        "inherited_source_pin_count": len(parent["source_pins"]),
        "source_union_hashes": source_union,
        "source_union_count": len(source_union),
        "stage_hashes": stage_hashes,
        "current_code_hashes": code,
        "implementation_hashes": dict(code),
        "parent_checkpoint_inventory": parent["inventory"],
        "parent_checkpoint_hashes": parent["parent_hashes"],
        "report_checkpoint_hashes": parent["own_hashes"],
        "complete_parent_checkpoint_hashes": parent["complete_hashes"],
        "parent_checkpoint_inventory_count": len(parent["inventory_hashes"]),
        "parent_checkpoint_hashes_count": len(parent["parent_hashes"]),
        "report_checkpoint_hashes_count": len(parent["own_hashes"]),
        "complete_parent_checkpoint_count": len(parent["complete_hashes"]),
        "parent_checkpoint_count": parent["report"].get("checkpoint_count"),
        "preserved_failed_artifact_hashes_union": dict(parent["failed_union"]),
        "preserved_failed_artifact_count": len(parent["failed_union"]),
        "new_checkpoint_count": 4,
        "final_checkpoint_count": (
            parent["report"].get("checkpoint_count", 0) + 4
        ),
        "protected_file_hashes": dict(parent["protected"]),
        "legacy_dirty_path_hashes": dict(parent["dirty"]),
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "raw_only_inherited_sources": True,
        "history_decode_forbidden": True,
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def _paths():
    return dict(CHECKPOINT_NAMES)


def prepare(inputs_path):
    paths = _paths()
    require(
        not any(path.exists() for path in paths.values()),
        "Volume-error source-packet checkpoint directory already contains an attempt",
    )
    packet = _build_packet(inputs_path)
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


def _verify_packet(packet):
    require(
        packet.get("schema") == INPUT_SCHEMA
        and packet.get("component") == COMPONENT
        and packet.get("metadata_only") is True
        and packet.get("source_packet_definition_only") is True
        and packet.get("quantitative_packet_error_reduction_proved") is False
        and packet.get("runtime_model_metadata") is None
        and packet.get("original_source_evaluations") == 0,
        "Frozen volume-error source packet scope changed",
    )
    require(packet.get("input_count") == len(packet.get("input_spec", {})),
            "Frozen input specification count changed")
    code = packet.get("current_code_hashes")
    require(
        code == packet.get("implementation_hashes")
        and isinstance(code, dict)
        and set(code) == {RUNNER, KERNEL},
        "Frozen volume-error code map changed",
    )
    for name, expected in code.items():
        require(sha_path(_resolve(name)) == expected,
                "New code changed after prepare: " + name)
    source_union = packet.get("source_union_hashes")
    require(
        isinstance(source_union, dict)
        and packet.get("source_union_count") == len(source_union),
        "Frozen source union changed",
    )
    for name, expected in source_union.items():
        _raw_auth(name, expected)
    for name, expected in packet.get("protected_file_hashes", {}).items():
        _raw_auth(name, expected)
    for name, expected in packet.get("legacy_dirty_path_hashes", {}).items():
        _raw_auth(name, expected)


def _load_packet(decode_inputs=False):
    source_raw, packet = _read_json(SOURCE_CHECKPOINT)
    _verify_packet(packet)
    parent = _parent_scope()
    require(
        packet.get("parent_report_sha256") == PARENT_REPORT_SHA256
        and packet.get("parent_result_sha256") == PARENT_RESULT_SHA256
        and packet.get("parent_commit") == CURRENT_PARENT_COMMIT,
        "Frozen parent hashes or commit changed",
    )
    require(
        packet.get("inherited_source_pins") == parent["source_pins"]
        and packet.get("parent_checkpoint_inventory") == parent["inventory"]
        and packet.get("parent_checkpoint_hashes") == parent["parent_hashes"]
        and packet.get("report_checkpoint_hashes") == parent["own_hashes"]
        and packet.get("complete_parent_checkpoint_hashes")
        == parent["complete_hashes"]
        and packet.get("complete_parent_checkpoint_count") == 180
        and packet.get("protected_file_hashes") == parent["protected"]
        and packet.get("legacy_dirty_path_hashes") == parent["dirty"],
        "Frozen inherited metadata maps changed",
    )

    spec = _load_input_spec(packet["input_spec_path"])
    require(
        spec["raw_sha256"] == packet["input_spec_raw_sha256"]
        and spec["records"] == packet["input_spec"],
        "Frozen input specification changed",
    )
    expected_stage_hashes = _merge_hash_maps(
        {
            PARENT_REPORT: PARENT_REPORT_SHA256,
            PARENT_RESULT: PARENT_RESULT_SHA256,
            spec["path"]: spec["raw_sha256"],
            **packet["current_code_hashes"],
        },
        {
            name: metadata["raw_sha256"]
            for name, metadata in spec["records"].items()
        },
    )
    require(packet.get("stage_hashes") == expected_stage_hashes,
            "Frozen stage hash manifest changed")
    require(
        packet.get("source_union_hashes")
        == _merge_hash_maps(parent["source_pins"], expected_stage_hashes),
        "Frozen source union does not match inherited pins and stage inputs",
    )
    records = None
    if decode_inputs:
        records = {}
        for record_path, metadata in packet["input_spec"].items():
            records[record_path] = _decode_record(metadata)
    return source_raw, packet, parent, records


def _decode_record(metadata):
    name = metadata["path"]
    require(not _raw_name_forbidden(name),
            "Forbidden history artifact cannot be decoded: " + name)
    raw = _raw_auth(name, metadata["raw_sha256"])
    encoding = metadata["encoding"]
    if encoding == "gzip_json":
        decoded_raw = gzip.decompress(raw)
        text = decoded_raw.decode("utf-8")
        value = json.loads(text)
    elif encoding == "json":
        decoded_raw = raw
        text = decoded_raw.decode("utf-8")
        value = json.loads(text)
    elif encoding == "utf8":
        decoded_raw = raw
        text = decoded_raw.decode("utf-8")
        value = text
    else:
        raise AssertionError("Unsupported frozen input encoding: " + str(encoding))
    decoded_sha = sha_bytes(decoded_raw)
    expected_decoded = metadata.get("decoded_sha256")
    if expected_decoded is not None:
        require(decoded_sha == expected_decoded,
                "Changed decoded input artifact: " + name)
    return {
        "path": name,
        "raw_sha256": metadata["raw_sha256"],
        "encoding": encoding,
        "decoded_sha256": decoded_sha,
        "raw_authenticated_before_decode": True,
        "value": value,
    }


def _brief_prepared_source(packet, source_sha):
    require(BRIEF.exists(), "Resume brief is missing")
    brief_raw, brief = _read_json(BRIEF)
    node = brief.get("in_progress", {}).get(COMPONENT)
    require(isinstance(node, dict),
            "Resume brief has no volume-error prepared marker")
    source_rel = SOURCE_CHECKPOINT.relative_to(ROOT).as_posix()
    require(
        node.get("phase") in PREPARED_PHASES
        and node.get("source_checkpoint_sha256") == source_sha
        and node.get("source_checkpoint_path")
        in (source_rel, str(SOURCE_CHECKPOINT))
        and node.get("kernel_sha256") == packet["current_code_hashes"][KERNEL]
        and node.get("runner_sha256") == packet["current_code_hashes"][RUNNER],
        "Resume brief does not point to the prepared volume-error source packet",
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


def _source_for_derive(packet, parent, records):
    values = {name: record["value"] for name, record in records.items()}
    metadata = {
        name: {
            key: value
            for key, value in record.items()
            if key != "value"
        }
        for name, record in records.items()
    }
    return {
        "schema": INPUT_SCHEMA,
        "component": COMPONENT,
        "records": values,
        "record_metadata": metadata,
        "parent_report": _clone(parent["report"]),
        "parent_result": _clone(parent["result"]),
        "parent_report_sha256": packet["parent_report_sha256"],
        "parent_result_sha256": packet["parent_result_sha256"],
        "input_digest": packet["input_digest"],
        "input_spec": _clone(packet["input_spec"]),
        "input_spec_path": packet["input_spec_path"],
        "input_spec_sha256": packet["input_spec_raw_sha256"],
        "runtime_model_metadata": None,
        "original_source_evaluations": 0,
        "source_packet_definition_only": True,
        "quantitative_packet_error_reduction_proved": False,
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
        "source_packet_definition_only": True,
        "quantitative_packet_error_reduction_proved": False,
        "derive_entry_calls": calls,
        "derive_calls": calls,
        "entry_attempts": 1,
        "kernel_imported": import_succeeded,
        "source_functions_called": 0,
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
        result.update(
            {
                "caught_exception": True,
                "error_type": type(error).__name__,
                "error": str(error),
                "traceback": traceback.format_exc(),
            }
        )
    return result


def _normalise_result(value, packet, calls, import_succeeded):
    if isinstance(value, dict):
        result = _clone(value)
    else:
        result = _fallback(
            packet,
            TypeError("derive returned a non-mapping result"),
            calls=calls,
            import_succeeded=import_succeeded,
        )
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("passed", False)
    result.setdefault("source_only", True)
    result.setdefault("metadata_only", False)
    result.setdefault("input_digest", packet["input_digest"])
    result.setdefault("runtime_model_metadata", None)
    result.setdefault("original_source_evaluations", 0)
    result.setdefault("source_packet_definition_only", True)
    result.setdefault("quantitative_packet_error_reduction_proved", False)
    result.setdefault("entry_attempts", 1)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["kernel_imported"] = import_succeeded
    return result


def _validate_result(result, packet):
    require(
        isinstance(result, dict)
        and result.get("schema") == RESULT_SCHEMA
        and result.get("component") == COMPONENT
        and result.get("metadata_only") is False,
        "Volume-error result schema or metadata scope changed",
    )
    require(
        result.get("input_digest") == packet["input_digest"]
        and result.get("runtime_model_metadata") is None
        and result.get("original_source_evaluations") == 0
        and result.get("entry_attempts") == 1,
        "Volume-error execution scope changed",
    )
    if "source_packet_definition_only" in result:
        require(result["source_packet_definition_only"] is True,
                "Source-packet definition scope changed")
    if "quantitative_packet_error_reduction_proved" in result:
        require(result["quantitative_packet_error_reduction_proved"] is False,
                "Source packet was mislabeled as quantitative reduction")
    calls = result.get("derive_entry_calls")
    require(
        isinstance(calls, int)
        and not isinstance(calls, bool)
        and calls in (0, 1)
        and result.get("derive_calls") == calls
        and isinstance(result.get("kernel_imported"), bool),
        "Volume-error derive bookkeeping changed",
    )
    if result.get("passed") is True:
        require(calls == 1 and result.get("kernel_imported") is True,
                "A passed result requires one successful derive call")
    audit = result.get("source_inspection_audit")
    if isinstance(audit, dict):
        for key in (
            "original_source_evaluations",
            "old_source_consumer_calls",
            "old_proof_functions_called",
            "old_operator_calls",
            "old_source_normalizer_calls",
        ):
            if key in audit and audit[key] is not None:
                require(audit[key] == 0,
                        "Forbidden source consumer/operator call was recorded")


def _load_saved(packet):
    paths = _paths()
    require(
        paths["entered"].exists()
        and paths["result"].exists()
        and paths["diagnostics"].exists(),
        "Saved volume-error attempt is incomplete",
    )
    entered_raw, entered = _read_json(paths["entered"])
    result_raw, result = _read_json(paths["result"])
    journal_raw = paths["diagnostics"].read_bytes()
    source_sha = sha_path(paths["source"])
    require(
        entered.get("schema") == ENTERED_SCHEMA
        and entered.get("component") == COMPONENT
        and entered.get("source_checkpoint_sha256") == source_sha
        and entered.get("kernel_imported_before_entry") is False
        and entered.get("derive_calls_before_entry") == 0
        and entered.get("entry_attempts") == 1,
        "Entered volume-error checkpoint changed",
    )
    require(
        result.get("source_checkpoint_sha256") == source_sha
        and result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw)
        and result.get("diagnostic_checkpoint_sha256") == sha_bytes(journal_raw),
        "Saved volume-error result provenance changed",
    )
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
    require(
        paths["source"].exists()
        and not paths["entered"].exists()
        and not paths["result"].exists()
        and not paths["diagnostics"].exists(),
        "Frozen source is absent or volume-error attempt already entered",
    )
    source_raw, packet, parent, records = _load_packet(decode_inputs=True)
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
        _append_event(
            paths["diagnostics"], events, packet, "kernel_callback", payload
        )

    calls = 0
    import_succeeded = False
    try:
        source = _source_for_derive(packet, parent, records)
        spec = importlib.util.spec_from_file_location(
            "_original_compact_volume_error_source_packet_entry",
            _resolve(KERNEL),
        )
        require(spec is not None and spec.loader is not None,
                "Cannot create new volume-error import specification")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        import_succeeded = True
        derive = getattr(module, "derive", None)
        require(callable(derive),
                "New volume-error module lacks derive(source, journal_callback)")
        calls = 1
        value = derive(source, journal_callback)
        result = _normalise_result(value, packet, calls, import_succeeded)
        _append_event(paths["diagnostics"], events, packet, "derive_returned", {
            "derive_entry_calls": calls,
            "passed": result.get("passed"),
            "runtime_model_metadata": None,
        })
    except Exception as error:
        result = _fallback(
            packet,
            error,
            calls=calls,
            import_succeeded=import_succeeded,
        )
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
    # The complete result and journal are durable before metadata validation.
    _validate_result(result, packet)
    return {
        "passed": result.get("passed") is True,
        "metadata_only": True,
        "component": COMPONENT,
        "result": paths["result"].relative_to(ROOT).as_posix(),
        "result_sha256": sha_path(paths["result"]),
        "attempt_status": "passed" if result.get("passed") else "failed",
        "derive_entry_calls": result.get("derive_entry_calls"),
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
        "attempt_status": (
            "passed" if saved["result"].get("passed") else "failed"
        ),
        "derive_entry_calls": saved["result"].get("derive_entry_calls"),
        "runtime_model_metadata": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "inspect"))
    parser.add_argument(
        "--inputs", help="Absolute JSON input specification for prepare"
    )
    args = parser.parse_args()
    if args.action == "prepare":
        require(args.inputs is not None,
                "prepare requires --inputs ABSOLUTE_PATH")
        value = prepare(args.inputs)
    else:
        require(args.inputs is None, "--inputs is valid only with prepare")
        value = globals()[args.action]()
    print(json.dumps(value, ensure_ascii=True, allow_nan=False))


if __name__ == "__main__":
    main()
