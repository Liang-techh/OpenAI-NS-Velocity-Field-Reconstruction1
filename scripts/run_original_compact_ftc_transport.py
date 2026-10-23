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


ROOT = Path(
    "C:/Users/z5242/Documents/Codex/2026-09-13/"
    "before-doing-any-project-work-verify/work/st073-transition-next"
)
RUNNER = "scripts/run_original_compact_ftc_transport.py"
KERNEL = "src/openai_ns_reconstruction/original_compact_ftc_transport.py"
PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_signed_pullback_difference.json"
)
PARENT_SHA256 = (
    "4a44b291afb4a0de8f82c652b0e99b4cd6e02eece8c86f44e2e94c76c5e4b9f7"
)
PARENT_COMMIT = "41d2a5d70707bdaf62dd10e9662a62d9d2a4dfbf"
PARENT_REPORT_BASE_COMMIT = "128d5818eae7c2f394b591b6da2fe7368827c532"

COMPACT_C5_RECORD = (
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
)
N1_PRERP_RECORD = (
    "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz"
)
RECORD_NAMES = (
    PARENT_REPORT,
    COMPACT_C5_RECORD,
    N1_PRERP_RECORD,
)
RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    COMPACT_C5_RECORD: "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269",
    N1_PRERP_RECORD: "3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce",
}
DECODED_RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    COMPACT_C5_RECORD: "b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510",
    N1_PRERP_RECORD: "348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64",
}

COMPONENT = "compact_ftc_transport"
INPUT_SCHEMA = "original_compact_ftc_transport_inputs_v1"
ENTERED_SCHEMA = "original_compact_ftc_transport_entered_v1"
RESULT_SCHEMA = "original_compact_ftc_transport"
DIAGNOSTIC_SCHEMA = "compact_ftc_transport_stage_diagnostic_v1"

CHECKPOINT_DIR = ROOT / "tmp/original_compact_ftc_transport"
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"

INHERITED_SOURCE_COUNT = 1885
FINAL_SOURCE_COUNT = 1888
AUTHENTICATED_OLD_CHECKPOINT_COUNT = 152
ACCEPTED_PARENT_CURRENT_CHECKPOINT_COUNT = 4
FAILED_PARENT_CURRENT_CHECKPOINT_COUNT = 4
PARENT_REPORT_CHECKPOINT_COUNT = 4
PARENT_INVENTORY_COUNT = 164
INHERITED_CHECKPOINT_COUNT = 168
NEW_CHECKPOINT_COUNT = 4
FINAL_CHECKPOINT_COUNT = 172
PRESERVED_FAILED_ARTIFACT_COUNT = 38
PARENT_LOSSLESS_CHECKPOINT_PATHS = {
    "source_checkpoint": "tmp/original_signed_pullback_difference/source.json",
    "entry_checkpoint": "tmp/original_signed_pullback_difference/entered.json",
    "result_checkpoint": "tmp/original_signed_pullback_difference/result.json",
    "diagnostic_checkpoint": "tmp/original_signed_pullback_difference/diagnostics.jsonl",
}
FINITE_ALLOWED_OPERATIONS = [
    "exact_rational",
    "negative",
    "positive_quotient",
    "product",
    "sum",
]
FINITE_WHITELIST_CONTEXT_COUNT = 5
FINITE_WHITELIST_HARD_CAP = 1024

ACCEPTED_NATIVE_B_NODES = list(range(520, 536))
FALSE_FLAGS = (
    "old_errors_removed",
    "old_error_terms_removed",
    "quantitative_packet_error_reduction_proved",
    "native_old_density_identity",
    "native_old_density_closure",
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
    "source_ftc_zero_identity",
    "old_pointwise_S_recovered",
    "old_mass_relation_proved",
)
SOURCE_FTC_GATE_FLAGS = (
    "source_bound_native_C1_FTC_proved",
    "source_bound_old_C1_IBP_proved",
    "source_bound_signed_transport_matches_saved_expression_C1_proved",
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


def _validate_finite_whitelist_contexts(contexts, label):
    require(isinstance(contexts, list)
            and len(contexts) == FINITE_WHITELIST_CONTEXT_COUNT,
            label + " must contain five finite-whitelist contexts")
    labels = set()
    for context in contexts:
        require(isinstance(context, dict), label + " context is not a mapping")
        context_label = context.get("label")
        require(isinstance(context_label, str) and context_label
                and context_label not in labels,
                label + " context labels changed")
        labels.add(context_label)
        row_ids = context.get("observed_allowed_row_ids")
        require(isinstance(row_ids, list)
                and all(isinstance(row_id, int) and not isinstance(row_id, bool)
                        and row_id >= 0 for row_id in row_ids)
                and row_ids == sorted(set(row_ids)),
                label + " context row IDs changed")
        require(context.get("observed_allowed_row_count") == len(row_ids)
                and len(row_ids) <= FINITE_WHITELIST_HARD_CAP
                and context.get("hard_cap") == FINITE_WHITELIST_HARD_CAP,
                label + " context row count or hard cap changed")
        require(context.get("allowed_operations") == FINITE_ALLOWED_OPERATIONS
                and context.get("unsupported_operations_rejected") is True,
                label + " context operation whitelist changed")
        archive_sha = context.get("archive_decoded_sha256")
        require(isinstance(archive_sha, str)
                and archive_sha in set(DECODED_RECORD_HASHES.values()),
                label + " context archive hash changed")
        observed_json = json.dumps(row_ids, separators=(",", ":"))
        require(context.get("observed_allowed_row_ids_sha256")
                == sha_bytes(observed_json.encode("utf-8")),
                label + " context row-ID hash changed")


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


def _parent_lossless_checkpoint_map(parent):
    artifacts = parent.get("lossless_artifacts")
    require(isinstance(artifacts, list), "Parent lossless artifacts missing")
    by_role = {}
    for item in artifacts:
        role = item.get("role") if isinstance(item, dict) else None
        if role not in PARENT_LOSSLESS_CHECKPOINT_PATHS:
            continue
        require(role not in by_role, "Duplicate parent lossless checkpoint role")
        path = item.get("path")
        expected = item.get("sha256")
        require(path == PARENT_LOSSLESS_CHECKPOINT_PATHS[role]
                and isinstance(expected, str) and len(expected) == 64,
                "Parent lossless checkpoint artifact changed")
        by_role[role] = (path, expected)
    require(set(by_role) == set(PARENT_LOSSLESS_CHECKPOINT_PATHS),
            "Parent lossless checkpoint rows changed")
    result = {path: expected for path, expected in by_role.values()}
    require(len(result) == PARENT_REPORT_CHECKPOINT_COUNT,
            "Parent lossless checkpoint count changed")
    for name, expected in result.items():
        _raw_auth(name, expected)
    return result


def _parent_scope():
    parent_raw = _raw_auth(PARENT_REPORT, PARENT_SHA256)
    parent = json.loads(parent_raw.decode("utf-8"))
    require(
        parent.get("schema") == "original_signed_pullback_difference_stage_v1"
        and parent.get("component") == "signed_pullback_difference"
        and parent.get("source_only") is True
        and parent.get("accepted_identity") is False
        and parent.get("full_goal_complete") is False
        and parent.get("parent_commit") == PARENT_REPORT_BASE_COMMIT,
        "Parent report scope changed",
    )

    source_pins = _hash_map(parent.get("source_pins"), INHERITED_SOURCE_COUNT,
                            "Parent source pins")
    failed_union = _hash_map(
        parent.get("preserved_failed_artifact_hashes_union"),
        PRESERVED_FAILED_ARTIFACT_COUNT,
        "Parent failed-artifact union",
    )
    inventory = parent.get("parent_checkpoint_inventory")
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
    inventory_parent_report = _hash_map(
        inventory.get("parent_report_checkpoint_hashes"),
        PARENT_REPORT_CHECKPOINT_COUNT,
        "Parent-report inventory checkpoints",
    )
    lossless_current = _parent_lossless_checkpoint_map(parent)
    require(
        inventory.get("inherited_count") == PARENT_INVENTORY_COUNT
        and inventory.get("raw_sha_only") is True
        and inventory.get("history_decode_forbidden") is True,
        "Parent checkpoint raw-only policy changed",
    )
    checkpoint_union = {
        **authenticated,
        **accepted_current,
        **failed_current,
        **inventory_parent_report,
    }
    require(len(checkpoint_union) == PARENT_INVENTORY_COUNT,
            "Parent checkpoint inventory is not 164 distinct entries")
    require(len({**checkpoint_union, **lossless_current}) == INHERITED_CHECKPOINT_COUNT,
            "Inherited checkpoint union is not 168 distinct entries")
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
            "parent_report_checkpoint_hashes": inventory_parent_report,
            "parent_lossless_checkpoint_hashes": lossless_current,
            "authenticated_old_count": AUTHENTICATED_OLD_CHECKPOINT_COUNT,
            "accepted_parent_current_count": ACCEPTED_PARENT_CURRENT_CHECKPOINT_COUNT,
            "failed_parent_current_count": FAILED_PARENT_CURRENT_CHECKPOINT_COUNT,
            "parent_report_count": PARENT_REPORT_CHECKPOINT_COUNT,
            "parent_lossless_count": PARENT_REPORT_CHECKPOINT_COUNT,
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
    require((ROOT / KERNEL).exists(), "New compact FTC transport kernel is missing")
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
            "Final source pin manifest must contain 1888 entries")

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
        "source_refs_opaque": False,
        "primitive_profile_leaves_opaque": True,
        "source_record_dag_replay": False,
        "source_record_dag_replay_scope": (
            "False means no old consumer or whole-archive replay; selected "
            "archive rows may receive bounded symbolic expansion."
        ),
        "parent_checkpoint_inventory": parent["checkpoint_inventory"],
        "inherited_checkpoint_count": INHERITED_CHECKPOINT_COUNT,
        "new_checkpoint_count": NEW_CHECKPOINT_COUNT,
        "final_checkpoint_count": FINAL_CHECKPOINT_COUNT,
        "preserved_failed_artifact_hashes": dict(parent["failed_union"]),
        "preserved_failed_artifact_count": PRESERVED_FAILED_ARTIFACT_COUNT,
        "accepted_native_B_nodes": list(ACCEPTED_NATIVE_B_NODES),
        "remaining_native_B_nodes": [],
        "retained_construction_data": {
            "cumulative_terms": 212,
            "raw_signed_terms": 213,
            "reduced_terms": 183,
            "six_signed_channels_retained": True,
            "all_old_error_data_retained": True,
            "reserved_contribution_registered_once": True,
            "inherited_errors_deleted": False,
        },
        "retained_integral_term_inventory": {
            "K_terms": 17,
            "outer_B_terms": 2,
            "J_terms": 16,
            "K_terms_retained": True,
            "outer_B_terms_retained": True,
            "J_terms_retained": True,
            "reserved_contribution_registered_once": True,
        },
        "synthetic_rational_jet_controls": False,
        "original_source_evaluations": 0,
        "source_ftc_transport_defined": False,
        "native_old_density_identity": False,
        "old212_adapter_identity": False,
        "old212_quantclosure": False,
        "false_flags": list(FALSE_FLAGS),
        "runtime_model_metadata": None,
        "kernel_imported": False,
        "derive_entry_calls": 0,
        "source_functions_called": 0,
        "source_inspection_audit": {
            "synthetic_rational_arithmetic_only": False,
            "original_source_evaluations": 0,
            "source_records_decoded": len(RECORD_NAMES),
            "source_refs_opaque": False,
            "primitive_profile_leaves_opaque": True,
            "whole_archive_DAG_walks": 0,
            "new_source_consumer_calls": 0,
            "publication_source_consumer_calls": 0,
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
            "Authenticate inherited bytes only; never decode forbidden history "
            "or replay accepted ancestors."
        ),
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def prepare():
    paths = _paths()
    require(not any(path.exists() for path in paths.values()),
            "Compact FTC checkpoint directory already contains an attempt")
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
        "source_ftc_transport_defined": False,
        "source_refs_opaque": False,
        "primitive_profile_leaves_opaque": True,
        "source_record_dag_replay": False,
        "original_source_evaluations": 0,
        "finite_whitelist_enforced": False,
        "finite_whitelist_contexts": [],
        "synthetic_rational_jet_controls": False,
        "symbolic_source_row_evaluations": 0,
        "bounded_source_row_symbolic_expansions": 0,
        "new_symbolic_normalizations": 0,
        "derive_entry_calls": 0,
        "runtime_model_metadata": None,
        "runner_caught_exception": caught,
        "runtime_error": None if error is None else str(error),
        "error_type": None if error is None else type(error).__name__,
        "traceback": "" if error is None else "".join(traceback.format_exception(error)),
    }


def _result_with_harness_metadata(value, packet, caught, calls, elapsed,
                                  source_sha, entered_sha, journal_sha, events):
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
    result["source_refs_opaque"] = False
    result.setdefault("primitive_profile_leaves_opaque", True)
    result.setdefault("source_record_dag_replay", False)
    result.setdefault(
        "source_record_dag_replay_scope",
        "False means no old consumer or whole-archive replay; selected archive "
        "rows may receive bounded symbolic expansion.",
    )
    result.setdefault("accepted_native_B_nodes", list(ACCEPTED_NATIVE_B_NODES))
    result.setdefault("remaining_native_B_nodes", [])
    for gate in SOURCE_FTC_GATE_FLAGS:
        result.setdefault(gate, False)
    result.setdefault("native_old_density_identity", False)
    result.setdefault("old212_adapter_identity", False)
    result.setdefault("old212_quantclosure", False)
    result.setdefault("finite_whitelist_enforced", False)
    result.setdefault("finite_whitelist_contexts", [])
    result.setdefault(
        "source_ftc_transport_defined",
        not caught and isinstance(value, dict)
        and value.get("source_ftc_transport_defined") is True,
    )
    result["original_source_evaluations"] = 0
    result.setdefault("symbolic_source_row_evaluations", 0)
    result.setdefault("bounded_source_row_symbolic_expansions", 0)
    result.setdefault("new_symbolic_normalizations", 0)
    result["synthetic_rational_jet_controls"] = False
    result.setdefault("passed", not caught)
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    if not caught:
        result["source_ftc_transport_defined"] = all(
            result.get(gate) is True for gate in SOURCE_FTC_GATE_FLAGS
        )
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
    result["runtime_model_metadata"] = None
    audit = result.get("source_inspection_audit")
    if not isinstance(audit, dict):
        audit = {}
        result["source_inspection_audit"] = audit
    audit["synthetic_rational_arithmetic_only"] = False
    audit.setdefault("original_source_evaluations", 0)
    audit.setdefault("source_records_decoded", len(RECORD_NAMES))
    audit["source_refs_opaque"] = False
    audit.setdefault("primitive_profile_leaves_opaque", True)
    audit.setdefault("source_record_dag_replay", False)
    audit.setdefault("whole_archive_DAG_walks", 0)
    audit["new_source_consumer_calls"] = calls
    audit.setdefault("publication_source_consumer_calls", 0)
    audit.setdefault("old_source_consumer_calls", 0)
    audit.setdefault("old_proof_functions_called", 0)
    audit.setdefault("old_operator_calls", 0)
    audit.setdefault("numerical_evaluations", 0)
    audit.setdefault("quadratures", 0)
    audit.setdefault("recurrence_calls", 0)
    audit.setdefault("old_source_normalizer_calls", 0)
    audit["symbolic_source_row_evaluations"] = result.get(
        "symbolic_source_row_evaluations", 0)
    audit["bounded_source_row_symbolic_expansions"] = result.get(
        "bounded_source_row_symbolic_expansions", 0)
    audit["new_symbolic_normalizations"] = result.get(
        "new_symbolic_normalizations", 0)
    audit["finite_whitelist_enforced"] = result.get(
        "finite_whitelist_enforced") is True
    audit["finite_whitelist_contexts"] = _clone(
        result.get("finite_whitelist_contexts", []))
    audit["runtime_model_metadata"] = None
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
        "derive_calls": calls,
        "runtime_model_metadata": None,
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
            and result.get("source_refs_opaque") is False
            and result.get("primitive_profile_leaves_opaque") is True
            and result.get("source_record_dag_replay") is False,
            "Saved selected-record metadata changed")
    require(result.get("full_result_saved_before_validation") is True
            and result.get("runner_full_goal_complete") is False
            and result.get("kernel_imported_before_entry") is False,
            "Result durability metadata changed")
    require(result.get("accepted_native_B_nodes") == ACCEPTED_NATIVE_B_NODES
            and result.get("remaining_native_B_nodes") == [],
            "Accepted native B scope changed")
    require(result.get("native_old_density_identity") is False
            and result.get("original_source_evaluations") == 0
            and result.get("runtime_model_metadata") is None
            and isinstance(result.get("finite_whitelist_enforced"), bool)
            and isinstance(result.get("finite_whitelist_contexts"), list)
            and result.get("source_record_dag_replay") is False,
            "Source FTC/native-old distinction changed")
    require(result.get("old_errors_removed") is False
            and result.get("old_error_terms_removed") is False,
            "Old error-removal flags changed")
    for flag in FALSE_FLAGS:
        require(result.get(flag) is False, "False gate was raised: " + flag)
    if result.get("passed") is True:
        require(result.get("source_ftc_transport_defined") is True
                and all(result.get(gate) is True
                        for gate in SOURCE_FTC_GATE_FLAGS)
                and result.get("derive_entry_calls") == 1
                and result.get("derive_calls") == 1,
                "Source FTC derive did not complete exactly once")
    else:
        require(result.get("source_ftc_transport_defined") is False,
                "Failed FTC result cannot claim source transport")
    audit = result.get("source_inspection_audit")
    require(isinstance(audit, dict)
            and audit.get("original_source_evaluations") == 0
            and audit.get("synthetic_rational_arithmetic_only") is False
            and audit.get("source_records_decoded") == len(RECORD_NAMES)
            and audit.get("source_refs_opaque") is False
            and audit.get("primitive_profile_leaves_opaque") is True
            and audit.get("source_record_dag_replay") is False
            and audit.get("whole_archive_DAG_walks") == 0
            and audit.get("new_source_consumer_calls")
                == result.get("derive_entry_calls")
            and audit.get("publication_source_consumer_calls") == 0
            and audit.get("old_source_consumer_calls") == 0
            and audit.get("old_proof_functions_called") == 0
            and audit.get("old_operator_calls") == 0
            and audit.get("old_source_normalizer_calls") == 0
            and audit.get("numerical_evaluations") == 0
            and audit.get("quadratures") == 0
            and audit.get("recurrence_calls") == 0
            and isinstance(audit.get("finite_whitelist_enforced"), bool)
            and isinstance(audit.get("finite_whitelist_contexts"), list)
            and audit.get("runtime_model_metadata") is None,
            "Source inspection audit changed")
    if result.get("passed") is True:
        require(result.get("finite_whitelist_enforced") is True
                and audit.get("finite_whitelist_enforced") is True
                and audit.get("finite_whitelist_contexts")
                    == result.get("finite_whitelist_contexts"),
                "Finite whitelist proof gate is missing")
        _validate_finite_whitelist_contexts(
            result.get("finite_whitelist_contexts"), "Result")
        _validate_finite_whitelist_contexts(
            audit.get("finite_whitelist_contexts"), "Source audit")
        require(isinstance(result.get("bounded_source_row_symbolic_expansions"), int)
                and result.get("bounded_source_row_symbolic_expansions") > 0
                and isinstance(result.get("new_symbolic_normalizations"), int)
                and result.get("new_symbolic_normalizations") > 0
                and audit.get("bounded_source_row_symbolic_expansions")
                    == result.get("bounded_source_row_symbolic_expansions")
                and audit.get("new_symbolic_normalizations")
                    == result.get("new_symbolic_normalizations"),
                "Bounded symbolic source-row audit counts are missing")
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
            "new_original_compact_ftc_transport", ROOT / KERNEL
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
            "derive_calls": calls,
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
        "derive_calls": result.get("derive_calls"),
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
        "derive_calls": result.get("derive_calls"),
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
