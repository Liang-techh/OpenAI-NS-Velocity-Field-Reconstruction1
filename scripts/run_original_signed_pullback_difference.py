"""Run the signed pullback-difference control once, with durable checkpoints.

``prepare`` authenticates the current parent report, selected compact source
records, and raw checkpoint pins before writing a frozen input manifest.
``run`` reloads those four records, writes the entry checkpoint before
importing the new kernel, calls its source-bound ``derive`` entry exactly
once, and saves the complete result and diagnostic journal before validation.
``inspect`` only reads saved metadata; it never imports the kernel or reruns
an entered attempt.

The source-bound graph references authenticated serialized records and uses
exact rational-jet arithmetic.  It does not execute original source functions,
produce numerical physical fields, or make a physical completion claim.
"""

from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import sys
import time
import traceback


ROOT = Path(__file__).resolve().parents[1]
RUNNER = "scripts/run_original_signed_pullback_difference.py"
KERNEL = "src/openai_ns_reconstruction/original_signed_pullback_difference.py"
PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_cross2_B_v2.json"
)
PARENT_SHA256 = (
    "342a685c1515329e6c04827f60e2e18d8cfef4d6bbf49087cf062e3fb3c78b67"
)
PARENT_COMMIT = "128d5818eae7c2f394b591b6da2fe7368827c532"

COMPACT_C5_RECORD = (
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
)
N1_PRERP_RECORD = (
    "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz"
)
NATIVE_COMPACT_RECORD = (
    "lei_ren_part1_paper_compliant_current_original_native_compact_pressure_polarization.json"
)
RECORD_NAMES = (
    PARENT_REPORT,
    COMPACT_C5_RECORD,
    N1_PRERP_RECORD,
    NATIVE_COMPACT_RECORD,
)
RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    COMPACT_C5_RECORD: "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
    N1_PRERP_RECORD: "3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce",
    NATIVE_COMPACT_RECORD: "be73a703c4d02c35f9470ffd7a8e694b4be3af0df39124ad6e44eebcc0390610",
}
DECODED_RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    COMPACT_C5_RECORD: "b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510",
    N1_PRERP_RECORD: "348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64",
    NATIVE_COMPACT_RECORD: "be73a703c4d02c35f9470ffd7a8e694b4be3af0df39124ad6e44eebcc0390610",
}

COMPONENT = "signed_pullback_difference"
INPUT_SCHEMA = "original_signed_pullback_difference_inputs_v1"
ENTERED_SCHEMA = "original_signed_pullback_difference_entered_v1"
RESULT_SCHEMA = "original_signed_pullback_difference"
DIAGNOSTIC_SCHEMA = "signed_pullback_difference_stage_diagnostic_v1"

CHECKPOINT_DIR = ROOT / "tmp/original_signed_pullback_difference"
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"

INHERITED_SOURCE_COUNT = 1882
FINAL_SOURCE_COUNT = 1885
AUTHENTICATED_OLD_CHECKPOINT_COUNT = 152
ACCEPTED_PARENT_CURRENT_CHECKPOINT_COUNT = 4
FAILED_PARENT_CURRENT_CHECKPOINT_COUNT = 4
PARENT_REPORT_CHECKPOINT_COUNT = 4
INHERITED_CHECKPOINT_COUNT = 164
NEW_CHECKPOINT_COUNT = 4
FINAL_CHECKPOINT_COUNT = 168
PRESERVED_FAILED_ARTIFACT_COUNT = 38

ACCEPTED_NATIVE_B_NODES = list(range(520, 536))
ORDINARY_Z_ORDERS = [0, 1]
OWNED_Z_CELL = ["0", "1/2"]

FALSE_FLAGS = (
    "old_errors_removed",
    "old_error_terms_removed",
    "quantitative_packet_error_reduction_proved",
    "native_old_density_identity",
    "global_five_moment_closure",
    "physical_layer_completed",
    "time_order_recursion_claimed",
    "full_goal_complete",
    "numeric_values_computed",
    "numeric_containment_claimed",
    "source_C5_installation_claimed",
    "new_loop_domain_admission_claimed",
    "admissible_cone_admission_claimed",
    "old_operators_replayed",
    "old212_adapter_identity",
    "old212_quantclosure",
)

# These names are used only as raw-name policy markers.  Matching artifacts
# are hashed as bytes and are never decoded by this runner.
RAW_ONLY_MARKERS = ("resume_first_switch", "rm_loop_density_identity")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def sha_path(path):
    return sha_bytes(Path(path).read_bytes())


def _json_default(value):
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    raise TypeError(f"Object of type {type(value).__name__} is not JSON serializable")


def _encode(value):
    return (
        json.dumps(
            value,
            ensure_ascii=True,
            allow_nan=False,
            sort_keys=True,
            default=_json_default,
        ).encode("utf-8")
        + b"\n"
    )


def _clone(value):
    return json.loads(_encode(value).decode("utf-8"))


def digest(value):
    return sha_bytes(
        json.dumps(
            value,
            ensure_ascii=True,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
            default=_json_default,
        ).encode("utf-8")
    )


def _write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(_encode(value))
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
    path = _resolve(name)
    require(path.exists(), "Missing authenticated raw artifact: " + str(name))
    raw = path.read_bytes()
    require(sha_bytes(raw) == expected, "Changed authenticated raw artifact: " + str(name))
    return raw


def _read_record(name, expected):
    """Authenticate one selected raw record before decoding its payload."""
    path = _resolve(name)
    raw = _raw_auth(name, expected)
    decoded = gzip.decompress(raw) if name.endswith(".gz") else raw
    return {
        "source_text": decoded.decode("utf-8"),
        "raw_sha256": sha_bytes(raw),
        "decoded_sha256": sha_bytes(decoded),
        "raw_authenticated_before_decode": True,
        "is_gzip": name.endswith(".gz"),
        "source_path": path.relative_to(ROOT).as_posix(),
    }


def _read_selected_records():
    records = {
        name: _read_record(name, RECORD_HASHES[name])
        for name in RECORD_NAMES
    }
    return records


def _record_metadata(records):
    for name in RECORD_NAMES:
        require(records[name]["decoded_sha256"] == DECODED_RECORD_HASHES[name],
                "Changed decoded selected-record hash: " + name)
    return {
        "record_names": list(RECORD_NAMES),
        "record_hashes": {
            name: records[name]["raw_sha256"] for name in RECORD_NAMES
        },
        "decoded_record_hashes": {
            name: records[name]["decoded_sha256"] for name in RECORD_NAMES
        },
        "record_compression": {
            name: records[name]["is_gzip"] for name in RECORD_NAMES
        },
        "record_source_paths": {
            name: records[name]["source_path"] for name in RECORD_NAMES
        },
        "record_count": len(RECORD_NAMES),
        "input_digest": digest(records),
    }


def _hash_map(value, count, label):
    require(isinstance(value, dict) and len(value) == count, label + " count changed")
    for name, expected in value.items():
        require(
            isinstance(name, str)
            and isinstance(expected, str)
            and len(expected) == 64,
            label + " contains a malformed SHA256",
        )
        _raw_auth(name, expected)
    return dict(value)


def _raw_only_name(name):
    normalized = name.replace("-", "_").lower()
    return any(marker in normalized for marker in RAW_ONLY_MARKERS)


def _parent_scope():
    # Hash the report bytes before decoding its JSON envelope.
    parent_path = _resolve(PARENT_REPORT)
    parent_raw = parent_path.read_bytes()
    require(sha_bytes(parent_raw) == PARENT_SHA256, "Changed authenticated parent report")
    parent = json.loads(parent_raw.decode("utf-8"))
    require(
        parent.get("schema") == "original_native_cross2_B_v2_stage_v1"
        and parent.get("component") == "cross2_B_source_v2"
        and parent.get("source_only") is True
        and parent.get("accepted_identity") is True
        and parent.get("full_goal_complete") is False
        and parent.get("parent_commit") == "2b70546657ff89d445111d8aef0f4a24f669037c",
        "Parent report scope changed",
    )

    source_pins = _hash_map(parent.get("source_pins"), INHERITED_SOURCE_COUNT,
                            "Parent source pins")
    failed_union = _hash_map(
        parent.get("preserved_failed_artifact_hashes_union"),
        PRESERVED_FAILED_ARTIFACT_COUNT,
        "Parent failed-artifact union",
    )
    inventory = parent.get("old_checkpoint_inventory")
    require(isinstance(inventory, dict), "Parent checkpoint inventory missing")
    authenticated = _hash_map(
        inventory.get("authenticated_checkpoint_hashes"),
        AUTHENTICATED_OLD_CHECKPOINT_COUNT,
        "Authenticated old checkpoints",
    )
    accepted_current = _hash_map(
        inventory.get("accepted_parent_current_checkpoint_hashes"),
        ACCEPTED_PARENT_CURRENT_CHECKPOINT_COUNT,
        "Accepted parent-current checkpoints",
    )
    failed_current = _hash_map(
        inventory.get("failed_parent_current_checkpoint_hashes"),
        FAILED_PARENT_CURRENT_CHECKPOINT_COUNT,
        "Failed parent-current checkpoints",
    )
    report_current = _hash_map(
        parent.get("checkpoint_hashes"),
        PARENT_REPORT_CHECKPOINT_COUNT,
        "Parent report checkpoints",
    )
    checkpoint_union = {
        **authenticated,
        **accepted_current,
        **failed_current,
        **report_current,
    }
    require(len(checkpoint_union) == INHERITED_CHECKPOINT_COUNT,
            "Inherited checkpoint union is not 164 distinct entries")
    require(
        inventory.get("raw_sha_only") is True
        and inventory.get("history_decode_forbidden") is True,
        "Parent checkpoint raw-only policy changed",
    )
    accepted_nodes = parent.get("accepted_native_B_nodes")
    require(accepted_nodes == ACCEPTED_NATIVE_B_NODES,
            "Parent accepted native B node set changed")

    # Every inherited pin and checkpoint is authenticated as bytes.  No source
    # record, forbidden history, or old checkpoint is JSON-decoded here.
    return {
        "raw": parent_raw,
        "report": parent,
        "source_pins": source_pins,
        "failed_union": failed_union,
        "checkpoint_inventory": {
            "authenticated_checkpoint_hashes": authenticated,
            "accepted_parent_current_checkpoint_hashes": accepted_current,
            "failed_parent_current_checkpoint_hashes": failed_current,
            "parent_report_checkpoint_hashes": report_current,
            "authenticated_old_count": AUTHENTICATED_OLD_CHECKPOINT_COUNT,
            "accepted_parent_current_count": ACCEPTED_PARENT_CURRENT_CHECKPOINT_COUNT,
            "failed_parent_current_count": FAILED_PARENT_CURRENT_CHECKPOINT_COUNT,
            "parent_report_count": PARENT_REPORT_CHECKPOINT_COUNT,
            "inherited_count": INHERITED_CHECKPOINT_COUNT,
            "raw_sha_only": True,
            "history_decode_forbidden": True,
        },
    }


def _paths():
    return {
        "source": SOURCE_CHECKPOINT,
        "entered": ENTERED_CHECKPOINT,
        "result": RESULT_CHECKPOINT,
        "diagnostics": DIAGNOSTIC_CHECKPOINT,
    }


def _build_packet():
    parent = _parent_scope()
    records = _read_selected_records()
    record_metadata = _record_metadata(records)
    for name in RECORD_NAMES[1:]:
        require(
            parent["source_pins"].get(name) == RECORD_HASHES[name],
            "Selected compact source is not inherited by the parent: " + name,
        )
    require((ROOT / KERNEL).exists(), "New signed pullback kernel is missing")
    current_code_hashes = {
        RUNNER: sha_path(ROOT / RUNNER),
        KERNEL: sha_path(ROOT / KERNEL),
    }
    final_union = dict(parent["source_pins"])
    require(PARENT_REPORT not in final_union and RUNNER not in final_union
            and KERNEL not in final_union,
            "New publication paths already occur in inherited source pins")
    final_union[PARENT_REPORT] = PARENT_SHA256
    final_union.update(current_code_hashes)
    require(len(final_union) == FINAL_SOURCE_COUNT,
            "Final source pin manifest must contain 1885 entries")

    manifest_digest = digest({
        "inherited_source_hashes": parent["source_pins"],
        "parent_report_hashes": {PARENT_REPORT: PARENT_SHA256},
        "current_code_hashes": current_code_hashes,
        "final_union_hashes": final_union,
        "selected_record_metadata": record_metadata,
    })
    return {
        "schema": INPUT_SCHEMA,
        "result_schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_only": True,
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_SHA256,
        "parent_commit": PARENT_COMMIT,
        "input_digest": record_metadata["input_digest"],
        "input_manifest_sha256": manifest_digest,
        "inherited_source_hashes": dict(parent["source_pins"]),
        "inherited_source_pin_count": INHERITED_SOURCE_COUNT,
        "source_pins": final_union,
        "source_pin_count": FINAL_SOURCE_COUNT,
        "input_source_pin_count": FINAL_SOURCE_COUNT,
        "source_pin_manifest": {
            "inherited_source_hashes": dict(parent["source_pins"]),
            "inherited_source_pin_count": INHERITED_SOURCE_COUNT,
            "parent_report_hashes": {PARENT_REPORT: PARENT_SHA256},
            "current_code_hashes": dict(current_code_hashes),
            "final_union_hashes": final_union,
            "final_union_count": FINAL_SOURCE_COUNT,
            "input_manifest_sha256": manifest_digest,
            "new_pins": [PARENT_REPORT, RUNNER, KERNEL],
            "new_pin_count": 3,
            "raw_only_inherited_sources": True,
            "history_decode_forbidden": True,
        },
        "current_code_hashes": current_code_hashes,
        "implementation_hashes": current_code_hashes,
        "record_names": list(RECORD_NAMES),
        "record_hashes": dict(record_metadata["record_hashes"]),
        "decoded_record_hashes": dict(record_metadata["decoded_record_hashes"]),
        "record_compression": dict(record_metadata["record_compression"]),
        "record_source_paths": dict(record_metadata["record_source_paths"]),
        "record_count": len(RECORD_NAMES),
        "record_input_digest": record_metadata["input_digest"],
        "selected_source_records": list(RECORD_NAMES[1:]),
        "source_refs_opaque": True,
        "source_record_dag_replay": False,
        "parent_checkpoint_inventory": parent["checkpoint_inventory"],
        "inherited_checkpoint_count": INHERITED_CHECKPOINT_COUNT,
        "new_checkpoint_count": NEW_CHECKPOINT_COUNT,
        "final_checkpoint_count": FINAL_CHECKPOINT_COUNT,
        "preserved_failed_artifact_hashes": dict(parent["failed_union"]),
        "preserved_failed_artifact_count": PRESERVED_FAILED_ARTIFACT_COUNT,
        "accepted_native_B_nodes": list(ACCEPTED_NATIVE_B_NODES),
        "remaining_native_B_nodes": [],
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "ordinary_Z_not_temporal_n": True,
        "owned_Z_cell": list(OWNED_Z_CELL),
        "synthetic_rational_jet_controls": True,
        "original_source_evaluations": 0,
        "native_old_density_identity": False,
        "false_flags": list(FALSE_FLAGS),
        "runtime_model_metadata": None,
        "kernel_imported": False,
        "derive_controls_calls": 0,
        "source_functions_called": 0,
        "source_inspection_audit": {
            "synthetic_rational_arithmetic_only": True,
            "original_source_evaluations": 0,
            "source_records_decoded": len(RECORD_NAMES),
            "source_refs_opaque": True,
            "dag_walks": 0,
            "old_source_consumer_calls": 0,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "numerical_evaluations": 0,
            "quadratures": 0,
            "recurrence_calls": 0,
            "runtime_model_metadata": None,
        },
        "raw_only_records": [],
        "checkpoint_policy": (
            "Authenticate inherited bytes only; do not decode forbidden history "
            "or replay accepted ancestors."
        ),
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def prepare():
    paths = _paths()
    require(not any(path.exists() for path in paths.values()),
            "Signed pullback checkpoint directory already contains an attempt")
    packet = _build_packet()
    _write_new(paths["source"], packet)
    return {
        "passed": True,
        "metadata_only": True,
        "component": COMPONENT,
        "source": paths["source"].relative_to(ROOT).as_posix(),
        "source_sha256": sha_path(paths["source"]),
        "input_digest": packet["input_digest"],
        "source_pin_count": FINAL_SOURCE_COUNT,
        "inherited_source_pin_count": INHERITED_SOURCE_COUNT,
        "checkpoint_count": FINAL_CHECKPOINT_COUNT,
        "runtime_model_metadata": None,
    }


def _load_packet():
    raw, packet = _read_json(SOURCE_CHECKPOINT)
    require(packet.get("schema") == INPUT_SCHEMA
            and packet.get("component") == COMPONENT,
            "Frozen input manifest schema changed")
    require(packet.get("source_pin_count") == FINAL_SOURCE_COUNT
            and packet.get("inherited_source_pin_count") == INHERITED_SOURCE_COUNT
            and packet.get("final_checkpoint_count") == FINAL_CHECKPOINT_COUNT,
            "Frozen input manifest counts changed")
    require(packet.get("parent_report") == PARENT_REPORT
            and packet.get("parent_report_sha256") == PARENT_SHA256
            and packet.get("parent_commit") == PARENT_COMMIT,
            "Frozen parent pointer changed")
    code = packet.get("current_code_hashes")
    require(code == packet.get("implementation_hashes")
            and isinstance(code, dict)
            and set(code) == {RUNNER, KERNEL},
            "Frozen implementation hash map changed")
    for name, expected in code.items():
        require(sha_path(ROOT / name) == expected,
                "Implementation changed after prepare: " + name)
    pins = packet.get("source_pins")
    require(isinstance(pins, dict) and len(pins) == FINAL_SOURCE_COUNT,
            "Frozen final source union changed")
    require(packet.get("source_pin_manifest", {}).get("final_union_hashes") == pins,
            "Frozen source pin manifest changed")
    require(packet.get("input_manifest_sha256")
            == packet.get("source_pin_manifest", {}).get("input_manifest_sha256"),
            "Frozen input manifest digest changed")
    # Authenticate all final pins as bytes.  This intentionally never decodes
    # inherited source records or the forbidden history artifacts.
    for name, expected in pins.items():
        _raw_auth(name, expected)
    require(packet.get("record_names") == list(RECORD_NAMES)
            and packet.get("selected_source_records") == list(RECORD_NAMES[1:])
            and packet.get("record_count") == len(RECORD_NAMES),
            "Frozen selected-record inventory changed")
    require(packet.get("record_hashes") == RECORD_HASHES,
            "Frozen selected-record raw hashes changed")
    records = _read_selected_records()
    metadata = _record_metadata(records)
    require(packet.get("decoded_record_hashes") == DECODED_RECORD_HASHES
            and packet.get("decoded_record_hashes") == metadata["decoded_record_hashes"]
            and packet.get("record_compression") == metadata["record_compression"]
            and packet.get("record_source_paths") == metadata["record_source_paths"]
            and packet.get("record_input_digest") == metadata["input_digest"]
            and packet.get("input_digest") == metadata["input_digest"],
            "Selected-record bytes or canonical input digest changed")
    return raw, packet, records


def _append_event(path, events, packet, stage, payload):
    frozen = _clone(payload)
    event = {
        "schema": DIAGNOSTIC_SCHEMA,
        "component": COMPONENT,
        "index": len(events),
        "stage": stage,
        "input_digest": packet["input_digest"],
        "payload": frozen,
    }
    with path.open("ab") as stream:
        stream.write(_encode(event))
        stream.flush()
        os.fsync(stream.fileno())
    events.append(event)


def _empty_file(path):
    with path.open("xb") as stream:
        stream.flush()
        os.fsync(stream.fileno())


def _fallback(packet, error=None):
    caught = error is not None
    return {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "input_digest": packet["input_digest"],
        "record_count": len(RECORD_NAMES),
        "record_input_digest": packet["input_digest"],
        "accepted_native_B_nodes": list(ACCEPTED_NATIVE_B_NODES),
        "remaining_native_B_nodes": [],
        "native_old_density_identity": False,
        "old212_adapter_identity": False,
        "old212_quantclosure": False,
        "source_bound_signed_compact_difference_defined": False,
        "source_record_dag_replay": False,
        "original_source_evaluations": 0,
        "synthetic_rational_jet_controls": True,
        "derive_entry_calls": 0,
        "derive_controls_calls": 0,
        "runtime_model_metadata": None,
        "runner_caught_exception": caught,
        "runtime_error": None if error is None else str(error),
        "error_type": None if error is None else type(error).__name__,
        "traceback": "" if error is None else "".join(traceback.format_exception(error)),
    }


def _result_with_harness_metadata(value, packet, caught, calls, elapsed,
                                  source_sha, entered_sha, journal_sha, events):
    graph_built = (
        not caught
        and isinstance(value, dict)
        and isinstance(value.get("signed_compact_graph"), dict)
        and isinstance(value.get("kernel_difference_projection"), dict)
    )
    if isinstance(value, dict):
        result = _clone(value)
    else:
        result = _fallback(packet, TypeError("derive returned a non-mapping result"))
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("source_only", True)
    result.setdefault("metadata_only", False)
    result.setdefault("input_digest", packet["input_digest"])
    result.setdefault("record_count", len(RECORD_NAMES))
    result.setdefault("record_input_digest", packet["input_digest"])
    result.setdefault("record_names", list(RECORD_NAMES))
    result.setdefault("record_hashes", _clone(packet["record_hashes"]))
    result.setdefault("decoded_record_hashes", _clone(packet["decoded_record_hashes"]))
    result.setdefault("source_refs_opaque", True)
    result.setdefault("source_record_dag_replay", False)
    result.setdefault("accepted_native_B_nodes", list(ACCEPTED_NATIVE_B_NODES))
    result.setdefault("remaining_native_B_nodes", [])
    result.setdefault("native_old_density_identity", False)
    result.setdefault("old212_adapter_identity", False)
    result.setdefault("old212_quantclosure", False)
    result.setdefault("source_bound_signed_compact_difference_defined", graph_built)
    result.setdefault("original_source_evaluations", 0)
    result.setdefault("synthetic_rational_jet_controls", True)
    result.setdefault("passed", not caught)
    result["derive_entry_calls"] = calls
    result.setdefault("derive_controls_calls", 1 if calls and not caught else 0)
    for flag in FALSE_FLAGS:
        result.setdefault(flag, False)
    result.setdefault("retained_construction_data", {
        "cumulative_terms": 212,
        "raw_signed_terms": 213,
        "reduced_terms": 183,
        "six_signed_channels_retained": True,
        "all_old_error_data_retained": True,
        "reserved_contribution_registered_once": True,
        "inherited_errors_deleted": False,
    })
    terms = result.get("retained_integral_term_inventory")
    if not isinstance(terms, dict):
        terms = {}
        result["retained_integral_term_inventory"] = terms
    terms.setdefault("K_terms", 17)
    terms.setdefault("outer_B_terms", 2)
    terms.setdefault("J_terms", 16)
    terms.setdefault("K_terms_retained", True)
    terms.setdefault("outer_B_terms_retained", True)
    terms.setdefault("J_terms_retained", True)
    terms.setdefault("reserved_contribution_registered_once", True)
    result.setdefault("runtime_model_metadata", None)
    audit = result.get("source_inspection_audit")
    if not isinstance(audit, dict):
        audit = {}
        result["source_inspection_audit"] = audit
    audit.setdefault("synthetic_rational_arithmetic_only", True)
    audit.setdefault("original_source_evaluations", 0)
    audit.setdefault("source_records_decoded", len(RECORD_NAMES))
    audit.setdefault("source_refs_opaque", True)
    audit.setdefault("source_record_dag_replay", False)
    audit.setdefault("dag_walks", 0)
    audit.setdefault("old_source_consumer_calls", 0)
    audit.setdefault("old_proof_functions_called", 0)
    audit.setdefault("old_operator_calls", 0)
    audit.setdefault("numerical_evaluations", 0)
    audit.setdefault("quadratures", 0)
    audit.setdefault("recurrence_calls", 0)
    audit.setdefault("runtime_model_metadata", result.get("runtime_model_metadata"))
    result.update({
        "execution_seconds": elapsed,
        "source_checkpoint_sha256": source_sha,
        "entry_checkpoint_sha256": entered_sha,
        "diagnostic_checkpoint_path": DIAGNOSTIC_CHECKPOINT.relative_to(ROOT).as_posix(),
        "diagnostic_checkpoint_sha256": journal_sha,
        "diagnostic_events": list(events),
        "diagnostic_events_saved_before_validation": len(events),
        "full_result_saved_before_validation": True,
        "entered_code_immutable": True,
        "never_retry_entered_code": True,
        "runner_caught_exception": caught,
        "runner_full_goal_complete": False,
        "kernel_imported_before_entry": False,
        "derive_entry_calls": result.get("derive_entry_calls", calls),
        "derive_controls_calls": result.get("derive_controls_calls", 0),
        "runtime_model_metadata": result.get("runtime_model_metadata"),
    })
    return result


def _validate_result(result, packet):
    require(isinstance(result, dict), "Saved result is not a mapping")
    require(result.get("schema") == RESULT_SCHEMA
            and result.get("component") == COMPONENT
            and result.get("source_only") is True
            and result.get("metadata_only") is False
            and result.get("input_digest") == packet["input_digest"],
            "Saved result schema or input digest changed")
    require(result.get("record_count") == len(RECORD_NAMES)
            and result.get("record_input_digest") == packet["input_digest"]
            and result.get("record_names") == list(RECORD_NAMES)
            and result.get("record_hashes") == packet["record_hashes"]
            and result.get("decoded_record_hashes") == packet["decoded_record_hashes"]
            and result.get("source_refs_opaque") is True
            and result.get("source_record_dag_replay") is False,
            "Saved selected-record metadata changed")
    require(result.get("full_result_saved_before_validation") is True
            and result.get("runner_full_goal_complete") is False,
            "Result durability metadata changed")
    require(result.get("accepted_native_B_nodes") == ACCEPTED_NATIVE_B_NODES
            and result.get("remaining_native_B_nodes") == [],
            "Accepted native B scope changed")
    require(result.get("native_old_density_identity") is False
            and result.get("original_source_evaluations") == 0
            and result.get("synthetic_rational_jet_controls") is True,
            "Synthetic/source distinction changed")
    require(result.get("old_errors_removed") is False
            and result.get("old_error_terms_removed") is False,
            "Old error-removal flags changed")
    for flag in FALSE_FLAGS:
        require(result.get(flag) is False, "False gate was raised: " + flag)
    if result.get("passed") is True:
        require(result.get("source_bound_signed_compact_difference_defined") is True
                and result.get("derive_entry_calls") == 1
                and result.get("derive_controls_calls") == 1,
                "Source-bound derive did not complete exactly once")
    audit = result.get("source_inspection_audit")
    require(isinstance(audit, dict)
            and audit.get("original_source_evaluations") == 0
            and audit.get("synthetic_rational_arithmetic_only") is True
            and audit.get("source_records_decoded") == len(RECORD_NAMES)
            and audit.get("source_refs_opaque") is True
            and audit.get("source_record_dag_replay") is False
            and audit.get("dag_walks") == 0
            and audit.get("old_source_consumer_calls") == 0
            and audit.get("old_proof_functions_called") == 0
            and audit.get("old_operator_calls") == 0
            and audit.get("numerical_evaluations") == 0
            and audit.get("quadratures") == 0
            and audit.get("recurrence_calls") == 0,
            "Source inspection audit changed")

    ledger = result.get("retained_construction_data")
    require(isinstance(ledger, dict)
            and ledger.get("cumulative_terms") == 212
            and ledger.get("raw_signed_terms") == 213
            and ledger.get("reduced_terms") == 183
            and ledger.get("six_signed_channels_retained") is True
            and ledger.get("all_old_error_data_retained") is True
            and ledger.get("reserved_contribution_registered_once") is True
            and ledger.get("inherited_errors_deleted") is False,
            "Retained 212/213/183 construction ledger changed")
    terms = result.get("retained_integral_term_inventory")
    require(isinstance(terms, dict)
            and terms.get("K_terms") == 17
            and terms.get("outer_B_terms") == 2
            and terms.get("J_terms") == 16
            and terms.get("K_terms_retained") is True
            and terms.get("outer_B_terms_retained") is True
            and terms.get("J_terms_retained") is True
            and terms.get("reserved_contribution_registered_once") is True,
            "Retained integral-term inventory changed")


def run():
    paths = _paths()
    require(paths["source"].exists() and not paths["entered"].exists()
            and not paths["result"].exists() and not paths["diagnostics"].exists(),
            "Frozen source is absent or already entered")
    source_raw, packet, records = _load_packet()
    # The entry checkpoint is durable before importing any new kernel code.
    _write_new(paths["entered"], {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_sha256": sha_bytes(source_raw),
        "input_digest": packet["input_digest"],
        "implementation_hashes": packet["current_code_hashes"],
        "kernel_sha256": packet["current_code_hashes"][KERNEL],
        "inherited_source_pin_count": INHERITED_SOURCE_COUNT,
        "record_names": list(RECORD_NAMES),
        "record_input_digest": packet["input_digest"],
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_imported_before_entry": False,
        "derive_calls_before_entry": 0,
        "entered_code_immutable": True,
        "never_retry_entered_code": True,
        "runtime_model_metadata": None,
    })
    _empty_file(paths["diagnostics"])
    events = []
    _append_event(paths["diagnostics"], events, packet, "entry_durable", {
        "kernel_imported_before_entry": False,
        "derive_calls_before_entry": 0,
        "record_count": len(RECORD_NAMES),
    })

    started = time.perf_counter()
    calls = 0
    caught = False
    value = None
    try:
        def persist_diagnostic(stage, payload=None):
            # The binding fragment emits one diagnostic item; accept the
            # conventional two-argument form as well for a stable kernel API.
            if payload is None and isinstance(stage, dict):
                item = stage
                stage = item.get("stage", "diagnostic")
                payload = item.get("payload", item)
            _append_event(paths["diagnostics"], events, packet,
                          "kernel." + str(stage), payload)

        _append_event(paths["diagnostics"], events, packet, "kernel_import_begin", {
            "kernel": KERNEL,
            "kernel_sha256": packet["current_code_hashes"][KERNEL],
        })
        spec = importlib.util.spec_from_file_location(
            "new_signed_pullback_difference", ROOT / KERNEL
        )
        require(spec is not None and spec.loader is not None,
                "Unable to construct new kernel import spec")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        derive = getattr(module, "derive", None)
        require(callable(derive), "New kernel does not expose source-bound derive")
        calls = 1
        value = derive(records, diagnostic_callback=persist_diagnostic)
        _append_event(paths["diagnostics"], events, packet, "derive_returned", {
            "returned_mapping": isinstance(value, dict),
            "returned_keys": sorted(value) if isinstance(value, dict) else [],
            "derive_entry_calls": calls,
            "derive_controls_calls": value.get("derive_controls_calls")
            if isinstance(value, dict) else None,
        })
    except Exception as error:
        caught = True
        value = _fallback(packet, error)
        _append_event(paths["diagnostics"], events, packet, "exception_caught", {
            "error_type": type(error).__name__,
            "error": str(error),
            "derive_entry_calls": calls,
        })

    elapsed = time.perf_counter() - started
    journal_sha = sha_path(paths["diagnostics"])
    result = _result_with_harness_metadata(
        value,
        packet,
        caught,
        calls,
        elapsed,
        sha_bytes(source_raw),
        sha_path(paths["entered"]),
        journal_sha,
        events,
    )
    # Save the complete result before checking any acceptance or ledger field.
    _write_new(paths["result"], result)
    _validate_result(result, packet)
    return {
        "passed": True,
        "metadata_only": True,
        "component": COMPONENT,
        "result": paths["result"].relative_to(ROOT).as_posix(),
        "result_sha256": sha_path(paths["result"]),
        "diagnostic_sha256": journal_sha,
        "derive_entry_calls": result.get("derive_entry_calls"),
        "derive_controls_calls": result.get("derive_controls_calls"),
        "runtime_model_metadata": result.get("runtime_model_metadata"),
    }


def inspect():
    paths = _paths()
    require(all(path.exists() for path in paths.values()),
            "Saved source, entry, result, and diagnostics are required")
    source_raw, packet = _read_json(paths["source"])
    entered_raw, entered = _read_json(paths["entered"])
    result_raw, result = _read_json(paths["result"])
    journal_raw = paths["diagnostics"].read_bytes()
    events = []
    if journal_raw:
        for line in journal_raw.splitlines():
            events.append(json.loads(line.decode("utf-8")))
    require(packet.get("schema") == INPUT_SCHEMA
            and packet.get("source_pin_count") == FINAL_SOURCE_COUNT
            and packet.get("record_names") == list(RECORD_NAMES)
            and packet.get("record_count") == len(RECORD_NAMES)
            and packet.get("record_input_digest") == packet.get("input_digest"),
            "Saved input metadata changed")
    require(entered.get("source_sha256") == sha_bytes(source_raw)
            and entered.get("kernel_imported_before_entry") is False
            and entered.get("derive_calls_before_entry") == 0
            and entered.get("record_names") == list(RECORD_NAMES)
            and entered.get("record_input_digest") == packet.get("input_digest"),
            "Saved entry metadata changed")
    require(result.get("source_checkpoint_sha256") == sha_bytes(source_raw)
            and result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw)
            and result.get("diagnostic_checkpoint_sha256") == sha_bytes(journal_raw)
            and result.get("diagnostic_events") == events
            and result.get("full_result_saved_before_validation") is True,
            "Saved result/journal metadata changed")
    return {
        "passed": True,
        "metadata_only": True,
        "component": COMPONENT,
        "source_sha256": sha_bytes(source_raw),
        "entered_sha256": sha_bytes(entered_raw),
        "result_sha256": sha_bytes(result_raw),
        "diagnostic_sha256": sha_bytes(journal_raw),
        "diagnostic_events": len(events),
        "derive_entry_calls": result.get("derive_entry_calls"),
        "derive_controls_calls": result.get("derive_controls_calls"),
        "runner_caught_exception": result.get("runner_caught_exception"),
        "runtime_model_metadata": result.get("runtime_model_metadata"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "inspect"))
    args = parser.parse_args()
    print(json.dumps(globals()[args.action](), ensure_ascii=True, allow_nan=False))


if __name__ == "__main__":
    main()
