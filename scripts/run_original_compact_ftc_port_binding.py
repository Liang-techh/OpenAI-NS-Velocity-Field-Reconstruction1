"""Run the compact FTC source-port binding leaf with durable checkpoints.

The runner authenticates the inherited report, source pins, and checkpoint
inventories as raw bytes before decoding the four selected source records.
``run`` writes the entry checkpoint before importing the new kernel, calls
``derive`` exactly once, and durably saves the complete result and diagnostic
journal before validation.  ``inspect`` reads saved metadata only.

The runner never imports an old checker, source provider, normalizer, proof
operator, or forbidden history artifact.  The new kernel owns the bounded
``PortCodec`` source recipes and the fresh symbolic proof work.
"""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
import argparse
import gzip
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
RUNNER = "scripts/run_original_compact_ftc_port_binding.py"
KERNEL = "src/openai_ns_reconstruction/original_compact_ftc_port_binding.py"

PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_compact_ftc_transport.json"
)
PARENT_SHA256 = (
    "c9857568a2e553c1963a94f670be4f6687a4b9c296dd2d10746b323aba64e033"
)
PARENT_COMMIT = "9534d2c460b4baf5d4d28759b6385e734e2a9e35"
PARENT_REPORT_PARENT_COMMIT = "41d2a5d70707bdaf62dd10e9662a62d9d2a4dfbf"

SIGNED_GRAPH_RECORD = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_signed_pullback_difference.json"
)
COMPACT_C5_RECORD = (
    "lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz"
)
N1_PRERP_RECORD = (
    "lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz"
)
RECORD_NAMES = (
    PARENT_REPORT,
    SIGNED_GRAPH_RECORD,
    COMPACT_C5_RECORD,
    N1_PRERP_RECORD,
)
RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    SIGNED_GRAPH_RECORD: (
        "4a44b291afb4a0de8f82c652b0e99b4cd6e02eece8c86f44e2e94c76c5e4b9f7"
    ),
    COMPACT_C5_RECORD: (
        "3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269"
    ),
    N1_PRERP_RECORD: (
        "3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce"
    ),
}
DECODED_RECORD_HASHES = {
    PARENT_REPORT: PARENT_SHA256,
    SIGNED_GRAPH_RECORD: (
        "4a44b291afb4a0de8f82c652b0e99b4cd6e02eece8c86f44e2e94c76c5e4b9f7"
    ),
    COMPACT_C5_RECORD: (
        "b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510"
    ),
    N1_PRERP_RECORD: (
        "348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64"
    ),
}

COMPONENT = "compact_ftc_port_binding"
INPUT_SCHEMA = "original_compact_ftc_port_binding_inputs_v1"
ENTERED_SCHEMA = "original_compact_ftc_port_binding_entered_v1"
RESULT_SCHEMA = "original_compact_ftc_port_binding"
DIAGNOSTIC_SCHEMA = "compact_ftc_port_binding_stage_diagnostic_v1"
CHECKPOINT_DIR = ROOT / "tmp/original_compact_ftc_port_binding"
SOURCE_CHECKPOINT = CHECKPOINT_DIR / "source.json"
ENTERED_CHECKPOINT = CHECKPOINT_DIR / "entered.json"
RESULT_CHECKPOINT = CHECKPOINT_DIR / "result.json"
DIAGNOSTIC_CHECKPOINT = CHECKPOINT_DIR / "diagnostics.jsonl"

INHERITED_SOURCE_COUNT = 1888
FINAL_SOURCE_COUNT = 1891
AUTHENTICATED_OLD_CHECKPOINT_COUNT = 152
ACCEPTED_PARENT_CURRENT_CHECKPOINT_COUNT = 4
FAILED_PARENT_CURRENT_CHECKPOINT_COUNT = 4
PARENT_REPORT_CHECKPOINT_COUNT = 4
PARENT_LOSSLESS_CHECKPOINT_COUNT = 4
PARENT_INVENTORY_COUNT = 168
INHERITED_CHECKPOINT_COUNT = 172
NEW_CHECKPOINT_COUNT = 4
FINAL_CHECKPOINT_COUNT = 176
PRESERVED_FAILED_ARTIFACT_COUNT = 44
FAILED_ATTEMPT_ARTIFACT_COUNT = 50

CHECKPOINT_MAP_KEYS = (
    ("authenticated_checkpoint_hashes", AUTHENTICATED_OLD_CHECKPOINT_COUNT),
    ("accepted_parent_current_checkpoint_hashes",
     ACCEPTED_PARENT_CURRENT_CHECKPOINT_COUNT),
    ("failed_parent_current_checkpoint_hashes",
     FAILED_PARENT_CURRENT_CHECKPOINT_COUNT),
    ("parent_report_checkpoint_hashes", PARENT_REPORT_CHECKPOINT_COUNT),
    ("parent_lossless_checkpoint_hashes", PARENT_LOSSLESS_CHECKPOINT_COUNT),
)

SOURCE_FTC_GATE_FLAGS = (
    "source_bound_native_C1_FTC_proved",
    "source_bound_old_C1_IBP_proved",
    "source_bound_signed_transport_matches_saved_expression_C1_proved",
)
RAW_ONLY_MARKERS = ("resume_first_switch", "rm_loop_density_identity")
FALSE_FLAGS = (
    "native_old_density_closure",
    "source_ftc_zero_identity",
    "native_old_density_identity",
    "old212_adapter_identity",
    "old212_quantclosure",
    "numeric_containment_claimed",
    "source_C5_installation_claimed",
    "new_loop_domain_admission_claimed",
    "admissible_cone_admission_claimed",
    "quantitative_packet_error_reduction_proved",
    "physical_layer_completed",
    "global_five_moment_closure",
    "time_order_recursion_claimed",
    "old_operators_replayed",
    "old_pointwise_S_recovered",
    "old_mass_relation_proved",
    "old_errors_removed",
    "old_error_terms_removed",
    "full_goal_complete",
    "numeric_values_computed",
)
ACCEPTED_NATIVE_B_NODES = list(range(520, 536))


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


def digest(value):
    return sha_bytes(json.dumps(value, ensure_ascii=True, sort_keys=True,
                                separators=(",", ":"), allow_nan=False)
                     .encode("utf-8"))


def _write_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(value if isinstance(value, bytes) else _encode(value))
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
    require(sha_bytes(raw) == expected, "Changed raw artifact: " + str(name))
    return raw


def _raw_only_name(name):
    normalized = str(name).replace("-", "_").lower()
    return any(marker in normalized for marker in RAW_ONLY_MARKERS)


def _read_record(name, expected):
    require(not _raw_only_name(name),
            "Forbidden history artifact cannot be decoded: " + str(name))
    raw = _raw_auth(name, expected)
    decoded = gzip.decompress(raw) if name.endswith(".gz") else raw
    return {
        "source_text": decoded.decode("utf-8"),
        "raw_sha256": sha_bytes(raw),
        "decoded_sha256": sha_bytes(decoded),
        "raw_authenticated_before_decode": True,
        "is_gzip": name.endswith(".gz"),
        "source_path": _resolve(name).relative_to(ROOT).as_posix(),
    }


def _read_selected_records():
    return {name: _read_record(name, RECORD_HASHES[name])
            for name in RECORD_NAMES}


def _record_metadata(records):
    for name in RECORD_NAMES:
        require(records[name]["raw_authenticated_before_decode"] is True,
                "Selected record was not raw-authenticated before decode")
        require(records[name]["decoded_sha256"] == DECODED_RECORD_HASHES[name],
                "Changed decoded selected record: " + name)
    return {
        "record_names": list(RECORD_NAMES),
        "record_hashes": {name: records[name]["raw_sha256"]
                          for name in RECORD_NAMES},
        "decoded_record_hashes": {name: records[name]["decoded_sha256"]
                                   for name in RECORD_NAMES},
        "record_compression": {name: records[name]["is_gzip"]
                               for name in RECORD_NAMES},
        "record_source_paths": {name: records[name]["source_path"]
                                for name in RECORD_NAMES},
        "record_count": len(RECORD_NAMES),
        "input_digest": digest(records),
    }


def _hash_map(value, count, label, authenticate=True):
    require(isinstance(value, dict) and len(value) == count,
            label + " count changed")
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
    """Flatten only mapping keys ending in ``hashes`` and detect collisions."""
    flattened = {}

    def visit(node, path):
        if isinstance(node, list):
            for index, item in enumerate(node):
                if isinstance(item, (dict, list)):
                    visit(item, path + "[" + str(index) + "]")
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
                    require(name not in flattened,
                            "Checkpoint hash collision at " + name)
                    flattened[name] = expected
            elif isinstance(child, dict):
                visit(child, child_path)
            elif isinstance(child, list):
                visit(child, child_path)

    visit(value, label)
    return flattened


def _checkpoint_scope(parent):
    inventory = parent.get("parent_checkpoint_inventory")
    require(isinstance(inventory, dict), "Parent checkpoint inventory missing")
    for key, count in CHECKPOINT_MAP_KEYS:
        _hash_map(inventory.get(key), count, "Parent " + key,
                  authenticate=False)
    inventory_hashes = _flatten_hash_maps(inventory, "parent_checkpoint_inventory")
    require(len(inventory_hashes) == PARENT_INVENTORY_COUNT,
            "Parent checkpoint inventory union count changed")
    own_hashes = _flatten_hash_maps(
        {"checkpoint_hashes": parent.get("checkpoint_hashes")},
        "parent_report_checkpoint_hashes",
    )
    require(len(own_hashes) == PARENT_REPORT_CHECKPOINT_COUNT,
            "Parent report checkpoint count changed")
    require(not set(inventory_hashes).intersection(own_hashes),
            "Parent inventory and report checkpoint hashes collide")
    for name, expected in {**inventory_hashes, **own_hashes}.items():
        _raw_auth(name, expected)
    require(parent.get("checkpoint_count") == INHERITED_CHECKPOINT_COUNT,
            "Parent checkpoint count changed")
    require(inventory.get("inherited_count") == PARENT_INVENTORY_COUNT
            and inventory.get("raw_sha_only") is True
            and inventory.get("history_decode_forbidden") is True,
            "Parent checkpoint raw-only policy changed")
    return {
        "inventory": _clone(inventory),
        "inventory_hashes": inventory_hashes,
        "own_hashes": own_hashes,
        "inherited_count": INHERITED_CHECKPOINT_COUNT,
    }


def _parent_scope():
    parent_raw = _raw_auth(PARENT_REPORT, PARENT_SHA256)
    parent = json.loads(parent_raw.decode("utf-8"))
    require(parent.get("schema") == "original_compact_ftc_transport_stage_v1"
            and parent.get("component") == "compact_ftc_transport"
            and parent.get("source_only") is True
            and parent.get("accepted_identity") is False
            and parent.get("full_goal_complete") is False
            and parent.get("parent_commit") == PARENT_REPORT_PARENT_COMMIT
            and parent.get("source_pin_count") == INHERITED_SOURCE_COUNT
            and parent.get("checkpoint_count") == INHERITED_CHECKPOINT_COUNT
            and parent.get("preserved_failed_artifact_count")
            == PRESERVED_FAILED_ARTIFACT_COUNT,
            "Parent compact FTC transport report scope changed")
    source_pins = _hash_map(parent.get("source_pins"), INHERITED_SOURCE_COUNT,
                            "Parent source pins")
    failed_union = _hash_map(
        parent.get("preserved_failed_artifact_hashes_union"),
        PRESERVED_FAILED_ARTIFACT_COUNT,
        "Parent failed-artifact union",
    )
    checkpoint = _checkpoint_scope(parent)
    formal = parent.get("formal_algebra")
    require(parent.get("formal_algebra_cached") is True
            and isinstance(formal, dict)
            and formal.get("passed") is True
            and formal.get("fresh_normalizations") == 32
            and formal.get("negative_control_count") == 9
            and formal.get("sympy_version") == "1.14.0",
            "Parent cached formal proof summary changed")
    require(parent.get("retained_construction_data", {}).get(
        "all_old_error_data_retained") is True,
        "Parent old error data is not retained")
    require(parent.get("native_old_density_identity") is False
            and parent.get("old212_adapter_identity") is False
            and parent.get("old212_quantclosure") is False,
            "Parent identity flags changed")
    return {
        "raw": parent_raw,
        "report": parent,
        "source_pins": source_pins,
        "failed_union": failed_union,
        "checkpoint": checkpoint,
        "formal": _clone(formal),
    }


def _build_packet():
    parent = _parent_scope()
    require((ROOT / KERNEL).exists(), "New compact FTC port-binding kernel missing")
    code = {RUNNER: sha_path(ROOT / RUNNER), KERNEL: sha_path(ROOT / KERNEL)}
    new_pins = {PARENT_REPORT: PARENT_SHA256, **code}
    final_pins = dict(parent["source_pins"])
    for name, expected in new_pins.items():
        require(name not in final_pins, "New source pin already inherited: " + name)
        final_pins[name] = expected
    require(len(final_pins) == FINAL_SOURCE_COUNT,
            "Final source pin manifest must contain 1891 entries")
    # Re-authenticate the complete source union before decoding any selected
    # record, including the new parent report and both code blobs.
    for name, expected in final_pins.items():
        _raw_auth(name, expected)
    records = _read_selected_records()
    metadata = _record_metadata(records)
    for name in RECORD_NAMES[1:]:
        require(parent["source_pins"].get(name) == RECORD_HASHES[name],
                "Selected inherited source changed: " + name)
    manifest_digest = digest({
        "inherited_source_hashes": parent["source_pins"],
        "new_source_hashes": new_pins,
        "final_union_hashes": final_pins,
        "selected_record_metadata": metadata,
    })
    formal = {
        "passed": True,
        "fresh_normalizations": 32,
        "negative_control_count": 9,
        "sympy_version": "1.14.0",
        "cached_formal_proof_calls": 0,
        "rerun": False,
        "source": "parent_report.formal_algebra",
    }
    return {
        "schema": INPUT_SCHEMA,
        "result_schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "metadata_only": True,
        "source_only": True,
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_SHA256,
        "parent_commit": PARENT_COMMIT,
        "parent_report_parent_commit": PARENT_REPORT_PARENT_COMMIT,
        "input_digest": metadata["input_digest"],
        "input_manifest_sha256": manifest_digest,
        "inherited_source_hashes": dict(parent["source_pins"]),
        "inherited_source_pin_count": INHERITED_SOURCE_COUNT,
        "source_pins": final_pins,
        "source_pin_count": FINAL_SOURCE_COUNT,
        "source_pin_manifest": {
            "inherited_source_hashes": dict(parent["source_pins"]),
            "inherited_source_pin_count": INHERITED_SOURCE_COUNT,
            "new_source_hashes": new_pins,
            "current_code_hashes": dict(code),
            "final_union_hashes": final_pins,
            "final_union_count": FINAL_SOURCE_COUNT,
            "input_manifest_sha256": manifest_digest,
            "new_pins": list(new_pins),
            "new_pin_count": len(new_pins),
            "raw_only_inherited_sources": True,
            "history_decode_forbidden": True,
        },
        "current_code_hashes": code,
        "implementation_hashes": code,
        "record_names": list(RECORD_NAMES),
        "record_hashes": dict(metadata["record_hashes"]),
        "decoded_record_hashes": dict(metadata["decoded_record_hashes"]),
        "record_compression": dict(metadata["record_compression"]),
        "record_source_paths": dict(metadata["record_source_paths"]),
        "record_count": len(RECORD_NAMES),
        "record_input_digest": metadata["input_digest"],
        "selected_source_records": list(RECORD_NAMES),
        "selected_record_raw_authenticated_before_decode": True,
        "source_refs_opaque": False,
        "primitive_profile_leaves_opaque": True,
        "source_record_dag_replay": False,
        "source_record_dag_replay_scope": (
            "False means no old consumer or whole-archive replay; the new "
            "kernel may perform bounded symbolic expansion of selected rows."
        ),
        "parent_checkpoint_inventory": parent["checkpoint"]["inventory"],
        "parent_checkpoint_hashes": parent["checkpoint"]["own_hashes"],
        "inherited_checkpoint_count": INHERITED_CHECKPOINT_COUNT,
        "new_checkpoint_count": NEW_CHECKPOINT_COUNT,
        "final_checkpoint_count": FINAL_CHECKPOINT_COUNT,
        "preserved_failed_artifact_hashes": dict(parent["failed_union"]),
        "preserved_failed_artifact_count": PRESERVED_FAILED_ARTIFACT_COUNT,
        "accepted_native_B_nodes": list(ACCEPTED_NATIVE_B_NODES),
        "remaining_native_B_nodes": [],
        "retained_construction_data": _clone(
            parent["report"]["retained_construction_data"]),
        "retained_integral_term_inventory": _clone(
            parent["report"].get("retained_integral_term_inventory", {})),
        "cached_formal_algebra": formal,
        "cached_formal_proof_calls": 0,
        "source_port_codec_policy": {
            "codec": "PortCodec",
            "compile_algorithm": (
                "breadth-first fixed-root closure to seeds; iterative Kahn "
                "topological compilation with ascending ready IDs"
            ),
            "recipe_manifests_are_kernel_owned": True,
            "per_recipe_slice_cap": 1024,
            "physical_factor_rule": (
                "exactly one authenticated Pstar factor and one expected "
                "integer field child; preserve raw argument order"
            ),
            "cached_formal_proof_calls": 0,
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
            "original_source_evaluations": 0,
            "source_records_decoded": len(RECORD_NAMES),
            "source_refs_opaque": False,
            "primitive_profile_leaves_opaque": True,
            "source_record_dag_replay": False,
            "whole_archive_DAG_walks": 0,
            "new_source_consumer_calls": 0,
            "publication_source_consumer_calls": 0,
            "old_source_consumer_calls": 0,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "old_source_normalizer_calls": 0,
            "numerical_evaluations": 0,
            "quadratures": 0,
            "recurrence_calls": 0,
            "cached_formal_proof_calls": 0,
            "runtime_model_metadata": None,
        },
        "raw_only_records": [],
        "checkpoint_policy": (
            "Authenticate inherited checkpoint bytes only; never decode "
            "forbidden history or replay accepted ancestors."
        ),
        "prepared_at_utc": datetime.now(timezone.utc).isoformat(),
    }


def prepare():
    paths = _paths()
    require(not any(path.exists() for path in paths.values()),
            "Port-binding checkpoint directory already contains an attempt")
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


def _paths():
    return {
        "source": SOURCE_CHECKPOINT,
        "entered": ENTERED_CHECKPOINT,
        "result": RESULT_CHECKPOINT,
        "diagnostics": DIAGNOSTIC_CHECKPOINT,
    }


def _auth_packet_checkpoint_scope(packet):
    inventory = packet.get("parent_checkpoint_inventory")
    inventory_hashes = _flatten_hash_maps(inventory, "packet.parent_checkpoint_inventory")
    require(len(inventory_hashes) == PARENT_INVENTORY_COUNT,
            "Frozen parent checkpoint inventory changed")
    own = _hash_map(packet.get("parent_checkpoint_hashes"),
                    PARENT_REPORT_CHECKPOINT_COUNT,
                    "Frozen parent report checkpoints", authenticate=False)
    require(not set(inventory_hashes).intersection(own),
            "Frozen checkpoint maps collide")
    for name, expected in {**inventory_hashes, **own}.items():
        _raw_auth(name, expected)


def _load_packet(decode_records):
    raw, packet = _read_json(SOURCE_CHECKPOINT)
    require(packet.get("schema") == INPUT_SCHEMA
            and packet.get("component") == COMPONENT,
            "Frozen input manifest schema changed")
    require(packet.get("source_pin_count") == FINAL_SOURCE_COUNT
            and packet.get("final_checkpoint_count") == FINAL_CHECKPOINT_COUNT
            and packet.get("parent_report") == PARENT_REPORT
            and packet.get("parent_report_sha256") == PARENT_SHA256
            and packet.get("parent_commit") == PARENT_COMMIT,
            "Frozen input manifest scope changed")
    code = packet.get("current_code_hashes")
    require(code == packet.get("implementation_hashes")
            and isinstance(code, dict) and set(code) == {RUNNER, KERNEL},
            "Frozen implementation hash map changed")
    for name, expected in code.items():
        require(sha_path(ROOT / name) == expected,
                "Implementation changed after prepare: " + name)
    pins = packet.get("source_pins")
    require(isinstance(pins, dict) and len(pins) == FINAL_SOURCE_COUNT,
            "Frozen final source union changed")
    require(packet.get("source_pin_manifest", {}).get("final_union_hashes") == pins,
            "Frozen source pin manifest changed")
    for name, expected in pins.items():
        _raw_auth(name, expected)
    _auth_packet_checkpoint_scope(packet)
    require(packet.get("record_names") == list(RECORD_NAMES)
            and packet.get("record_count") == len(RECORD_NAMES)
            and packet.get("record_hashes") == RECORD_HASHES,
            "Frozen selected-record inventory changed")
    records = None
    if decode_records:
        records = _read_selected_records()
        metadata = _record_metadata(records)
        require(packet.get("decoded_record_hashes") == metadata["decoded_record_hashes"]
                and packet.get("record_input_digest") == metadata["input_digest"]
                and packet.get("input_digest") == metadata["input_digest"],
                "Selected record bytes or canonical input digest changed")
    return raw, packet, records


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


def _empty_file(path):
    with path.open("xb") as stream:
        stream.flush()
        os.fsync(stream.fileno())


def _formal_summary(packet):
    return _clone(packet["cached_formal_algebra"])


def _fallback(packet, error=None, calls=0, import_succeeded=False):
    result = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "passed": False,
        "source_only": True,
        "metadata_only": False,
        "input_digest": packet["input_digest"],
        "record_count": len(RECORD_NAMES),
        "record_input_digest": packet["input_digest"],
        "record_names": list(RECORD_NAMES),
        "record_hashes": dict(packet["record_hashes"]),
        "decoded_record_hashes": dict(packet["decoded_record_hashes"]),
        "source_ftc_transport_defined": False,
        "native_old_density_identity": False,
        "old212_adapter_identity": False,
        "old212_quantclosure": False,
        "full_goal_complete": False,
        "original_source_evaluations": 0,
        # A failed entry must not turn inherited parent metadata into a
        # validated result.  Keep the inherited summary under a distinct
        # field so publication can report the distinction truthfully.
        "cached_formal_algebra": None,
        "inherited_cached_formal_algebra": _formal_summary(packet),
        "cached_formal_proof_calls": 0,
        "fresh_source_proof_normalizations": 0,
        "recipe_manifests": {},
        "runtime_model_metadata": None,
        "derive_entry_calls": calls,
        "derive_calls": calls,
        "entry_attempts": 1,
        "kernel_imported": import_succeeded,
        "source_functions_called": 0,
        "retained_construction_data": _clone(
            packet["retained_construction_data"]),
        "source_inspection_audit": {
            "original_source_evaluations": 0,
            "whole_archive_DAG_walks": 0,
            "source_record_dag_replay": False,
            "new_source_consumer_calls": calls,
            "old_source_consumer_calls": 0,
            "old_proof_functions_called": 0,
            "old_operator_calls": 0,
            "old_source_normalizer_calls": 0,
            "cached_formal_proof_calls": 0,
            "runtime_model_metadata": None,
        },
    }
    for flag in SOURCE_FTC_GATE_FLAGS:
        result[flag] = False
    for flag in FALSE_FLAGS:
        result[flag] = False
    if error is not None:
        result["caught_exception"] = True
        result["error_type"] = type(error).__name__
        result["error"] = str(error)
        result["traceback"] = traceback.format_exc()
    return result


def _normalise_result(value, packet, calls, import_succeeded=True):
    result = dict(value) if isinstance(value, dict) else _fallback(
        packet, TypeError("derive returned a non-mapping result"),
        calls=calls, import_succeeded=import_succeeded)
    result.setdefault("schema", RESULT_SCHEMA)
    result.setdefault("component", COMPONENT)
    result.setdefault("source_only", True)
    result.setdefault("metadata_only", False)
    result.setdefault("input_digest", packet["input_digest"])
    result.setdefault("record_count", len(RECORD_NAMES))
    result.setdefault("record_input_digest", packet["input_digest"])
    result.setdefault("record_names", list(RECORD_NAMES))
    result.setdefault("record_hashes", dict(packet["record_hashes"]))
    result.setdefault("decoded_record_hashes", dict(packet["decoded_record_hashes"]))
    result.setdefault("runtime_model_metadata", None)
    result.setdefault("original_source_evaluations", 0)
    result.setdefault("cached_formal_algebra", None)
    result.setdefault("inherited_cached_formal_algebra", _formal_summary(packet))
    result.setdefault("cached_formal_proof_calls", 0)
    result.setdefault(
        "fresh_source_proof_normalizations",
        result.get("fresh_source_normalizations", 0),
    )
    result.setdefault("recipe_manifests", {})
    result.setdefault("source_functions_called", 0)
    result.setdefault("entry_attempts", 1)
    result.setdefault("retained_construction_data",
                      _clone(packet["retained_construction_data"]))
    result["derive_entry_calls"] = calls
    result["derive_calls"] = calls
    result["kernel_imported"] = import_succeeded
    for flag in SOURCE_FTC_GATE_FLAGS:
        result.setdefault(flag, False)
    result.setdefault("source_ftc_transport_defined", False)
    audit = result.setdefault("source_inspection_audit", {})
    for key, default in {
        "original_source_evaluations": 0,
        "whole_archive_DAG_walks": 0,
        "source_record_dag_replay": False,
        "new_source_consumer_calls": calls,
        "old_source_consumer_calls": 0,
        "old_proof_functions_called": 0,
        "old_operator_calls": 0,
        "old_source_normalizer_calls": 0,
        "cached_formal_proof_calls": 0,
        "runtime_model_metadata": None,
    }.items():
        audit.setdefault(key, default)
    return result


def _manifest_items(value):
    if isinstance(value, dict):
        if any(key in value for key in
               ("slice_cap", "cap", "closed_slice_cap")):
            yield value
        for child in value.values():
            yield from _manifest_items(child)
    elif isinstance(value, list):
        for child in value:
            yield from _manifest_items(child)


def _formal_counts(formal):
    require(isinstance(formal, dict),
            "Cached formal algebra summary is missing")
    fresh = formal.get(
        "fresh_normalizations",
        formal.get("fresh_normalizations_in_cached_attempt"),
    )
    negative = formal.get(
        "negative_control_count",
        formal.get("negative_controls_in_cached_attempt"),
    )
    return fresh, negative


def _validate_recipe_manifests(result, require_nonempty=True):
    manifests = [result.get(key) for key in (
        "source_recipe_result", "source_slices", "recipe_manifests",
        "codec_manifests", "source_recipe_manifests", "recipe_manifest",
    ) if result.get(key) is not None]
    if not manifests:
        require(not require_nonempty,
                "Fresh codec recipe manifests are missing")
        return
    items = list(_manifest_items(manifests))
    if not items:
        require(not require_nonempty,
                "Fresh codec recipe manifests contain no bounded slices")
        return
    for item in items:
        cap = item.get("slice_cap", item.get("closed_slice_cap",
                                             item.get("cap")))
        require(isinstance(cap, int) and not isinstance(cap, bool)
                and 0 < cap <= 1024,
                "Codec recipe slice cap exceeds 1024")
        for key in ("compiled_row_ids", "closed_slice_ids", "observed_row_ids",
                    "roots", "seeds"):
            values = item.get(key)
            if isinstance(values, list):
                require(len(values) <= cap,
                        "Codec recipe bounded row list exceeds its cap")


def _validate_result(result, packet):
    require(result.get("schema") == RESULT_SCHEMA
            and result.get("component") == COMPONENT
            and result.get("source_only") is True
            and result.get("metadata_only") is False,
            "Port-binding result schema or metadata scope changed")
    require(result.get("input_digest") == packet["input_digest"]
            and result.get("record_input_digest") == packet["input_digest"]
            and result.get("record_count") == len(RECORD_NAMES)
            and result.get("record_names") == list(RECORD_NAMES)
            and result.get("record_hashes") == packet["record_hashes"]
            and result.get("decoded_record_hashes")
            == packet["decoded_record_hashes"],
            "Port-binding result input provenance changed")
    derive_entry_calls = result.get("derive_entry_calls")
    derive_calls = result.get("derive_calls")
    kernel_imported = result.get("kernel_imported")
    require(result.get("runtime_model_metadata") is None
            and result.get("original_source_evaluations") == 0
            and result.get("cached_formal_proof_calls") == 0
            and isinstance(derive_entry_calls, int)
            and not isinstance(derive_entry_calls, bool)
            and derive_entry_calls in (0, 1)
            and derive_calls == derive_entry_calls
            and isinstance(kernel_imported, bool)
            and result.get("entry_attempts") == 1
            and result.get("source_functions_called") == 0,
            "Port-binding execution scope changed")
    if result.get("passed") is True:
        require(derive_entry_calls == 1 and kernel_imported is True,
                "A passed result requires one imported derive call")
    elif derive_entry_calls == 1:
        require(kernel_imported is True,
                "A recorded derive call requires a successful import")
    formal = result.get("cached_formal_algebra")
    if result.get("passed") is True or formal is not None:
        formal_fresh, formal_negative = _formal_counts(formal)
        require(isinstance(formal, dict)
                and formal.get("passed") is True
                and formal_fresh == 32
                and formal_negative == 9
                and formal.get("sympy_version") in (None, "1.14.0")
                and formal.get("cached_formal_proof_calls") == 0,
                "Cached formal proof summary changed")
    else:
        inherited = result.get("inherited_cached_formal_algebra")
        inherited_fresh, inherited_negative = _formal_counts(inherited)
        require(isinstance(inherited, dict)
                and inherited.get("passed") is True
                and inherited_fresh == 32
                and inherited_negative == 9
                and inherited.get("sympy_version") in (None, "1.14.0")
                and inherited.get("cached_formal_proof_calls") == 0,
                "Inherited formal proof metadata changed")
    fresh = result.get("fresh_source_normalizations",
                       result.get("fresh_source_proof_normalizations"))
    if result.get("passed") is True:
        require(isinstance(fresh, int) and not isinstance(fresh, bool)
                and fresh > 0,
                "Fresh source proof normalization count is missing")
    elif fresh is not None:
        require(isinstance(fresh, int) and not isinstance(fresh, bool)
                and fresh >= 0,
                "Failed result has an invalid fresh normalization count")
    controls = result.get("codec_controls")
    if result.get("passed") is True or controls is not None:
        require(isinstance(controls, dict)
                and controls.get("count") == 12
                and isinstance(controls.get("controls"), list)
                and len(controls["controls"]) == 12,
                "PortCodec control manifest changed")
        require(isinstance(controls.get("passed"), bool),
                "PortCodec control result is missing a boolean status")
        if result.get("passed") is True:
            require(controls.get("passed") is True,
                    "A passed result requires passed PortCodec controls")
    manifest_present = any(
        result.get(key) not in (None, {}, []) for key in (
            "source_recipe_result", "source_slices", "recipe_manifests",
            "codec_manifests", "source_recipe_manifests", "recipe_manifest",
        )
    )
    if result.get("passed") is True or manifest_present:
        _validate_recipe_manifests(
            result, require_nonempty=result.get("passed") is True)
    gates = [result.get(flag) is True for flag in SOURCE_FTC_GATE_FLAGS]
    require(result.get("source_ftc_transport_defined") is all(gates),
            "Source FTC transport flag does not match its three gates")
    require(result.get("passed") is all(gates),
            "Overall result does not match the three source FTC gates")
    retained = result.get("retained_construction_data")
    require(isinstance(retained, dict)
            and retained.get("all_old_error_data_retained") is True,
            "Old construction data was not retained")
    for flag in FALSE_FLAGS:
        require(result.get(flag) is False,
                "Forbidden completion flag became true: " + flag)
    audit = result.get("source_inspection_audit")
    require(isinstance(audit, dict)
            and audit.get("original_source_evaluations") == 0
            and audit.get("whole_archive_DAG_walks") == 0
            and audit.get("source_record_dag_replay") is False
            and audit.get("old_source_consumer_calls") == 0
            and audit.get("old_proof_functions_called") == 0
            and audit.get("old_operator_calls") == 0
            and audit.get("old_source_normalizer_calls") == 0
            and audit.get("cached_formal_proof_calls") == 0,
            "Source inspection audit scope changed")


def _load_attempt(parent):
    source_raw, packet, _ = _load_packet(decode_records=False)
    paths = _paths()
    require(paths["entered"].exists() and paths["result"].exists()
            and paths["diagnostics"].exists(),
            "Entered port-binding attempt is incomplete")
    entered_raw, entered = _read_json(paths["entered"])
    result_raw, result = _read_json(paths["result"])
    journal_raw = paths["diagnostics"].read_bytes()
    require(entered.get("schema") == ENTERED_SCHEMA
            and entered.get("component") == COMPONENT
            and entered.get("source_sha256") == sha_bytes(source_raw)
            and entered.get("kernel_imported_before_entry") is False
            and entered.get("derive_calls_before_entry") == 0
            and entered.get("entry_attempts") == 1
            and entered.get("record_names") == list(RECORD_NAMES),
            "Entered port-binding checkpoint changed")
    require(result.get("source_checkpoint_sha256") == sha_bytes(source_raw)
            and result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw)
            and result.get("diagnostic_checkpoint_sha256") == sha_bytes(journal_raw)
            and result.get("result_schema") == RESULT_SCHEMA,
            "Port-binding result checkpoint provenance changed")
    _validate_result(result, packet)
    artifacts = {
        "runner_source": (RUNNER, (ROOT / RUNNER).read_bytes()),
        "kernel_source": (KERNEL, (ROOT / KERNEL).read_bytes()),
        "source_checkpoint": (paths["source"].relative_to(ROOT).as_posix(), source_raw),
        "entry_checkpoint": (paths["entered"].relative_to(ROOT).as_posix(), entered_raw),
        "result_checkpoint": (paths["result"].relative_to(ROOT).as_posix(), result_raw),
        "diagnostic_checkpoint": (paths["diagnostics"].relative_to(ROOT).as_posix(), journal_raw),
    }
    artifact_rows = []
    for role, (name, raw) in artifacts.items():
        artifact_rows.append({"role": role, "path": name,
                              "sha256": sha_bytes(raw), "bytes": len(raw),
                              "raw_authenticated": True})
    failed_union = dict(parent["failed_union"])
    if result.get("passed") is False:
        for row in artifact_rows:
            failed_union[row["path"]] = row["sha256"]
    require(len(failed_union) == (PRESERVED_FAILED_ARTIFACT_COUNT
                                  if result.get("passed")
                                  else FAILED_ATTEMPT_ARTIFACT_COUNT),
            "Failed-artifact union count changed")
    return {
        "source_raw": source_raw,
        "packet": packet,
        "entered_raw": entered_raw,
        "entered": entered,
        "result_raw": result_raw,
        "result": result,
        "journal_raw": journal_raw,
        "artifacts": artifact_rows,
        "checkpoint_hashes": {
            paths["source"].relative_to(ROOT).as_posix(): sha_bytes(source_raw),
            paths["entered"].relative_to(ROOT).as_posix(): sha_bytes(entered_raw),
            paths["result"].relative_to(ROOT).as_posix(): sha_bytes(result_raw),
            paths["diagnostics"].relative_to(ROOT).as_posix(): sha_bytes(journal_raw),
        },
        "failed_union": failed_union,
        "attempt_status": "passed" if result.get("passed") else "failed",
    }


def run():
    paths = _paths()
    require(paths["source"].exists() and not paths["entered"].exists()
            and not paths["result"].exists()
            and not paths["diagnostics"].exists(),
            "Frozen source is absent or the attempt was already entered")
    _, packet, records = _load_packet(decode_records=True)
    entered = {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_sha256": sha_path(paths["source"]),
        "input_digest": packet["input_digest"],
        "record_names": list(RECORD_NAMES),
        "record_input_digest": packet["record_input_digest"],
        "current_code_hashes": dict(packet["current_code_hashes"]),
        "kernel_imported_before_entry": False,
        "derive_calls_before_entry": 0,
        "entry_attempts": 1,
        "entered_code_immutable": True,
        "never_retry_entered_code": True,
        "runtime_model_metadata": None,
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    _write_new(paths["entered"], entered)
    _empty_file(paths["diagnostics"])
    events = []

    def persist_diagnostic(stage, payload):
        _append_event(paths["diagnostics"], events, packet, stage, payload)

    calls = 0
    import_succeeded = False
    try:
        spec = importlib.util.spec_from_file_location(
            "_original_compact_ftc_port_binding_entry", ROOT / KERNEL
        )
        require(spec is not None and spec.loader is not None,
                "Cannot create a new kernel import specification")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        import_succeeded = True
        derive = getattr(module, "derive", None)
        require(callable(derive), "New kernel lacks derive(records, diagnostic_callback)")
        calls = 1
        value = derive(records, diagnostic_callback=persist_diagnostic)
        result = _normalise_result(value, packet, calls,
                                   import_succeeded=import_succeeded)
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
    # The complete result and diagnostic journal are durable before validation.
    _validate_result(result, packet)
    return {
        "passed": result.get("passed") is True,
        "metadata_only": True,
        "component": COMPONENT,
        "result": paths["result"].relative_to(ROOT).as_posix(),
        "result_sha256": sha_path(paths["result"]),
        "attempt_status": "passed" if result.get("passed") else "failed",
        "source_ftc_transport_defined": result.get("source_ftc_transport_defined"),
        "derive_entry_calls": result.get("derive_entry_calls"),
        "cached_formal_proof_calls": result.get("cached_formal_proof_calls"),
        "runtime_model_metadata": None,
    }


def inspect():
    paths = _paths()
    source_raw, packet, _ = _load_packet(decode_records=False)
    entered_raw, entered = _read_json(paths["entered"])
    result_raw, result = _read_json(paths["result"])
    journal_raw = paths["diagnostics"].read_bytes()
    require(entered.get("source_sha256") == sha_bytes(source_raw)
            and entered.get("kernel_imported_before_entry") is False
            and entered.get("derive_calls_before_entry") == 0
            and entered.get("entry_attempts") == 1,
            "Saved entered checkpoint changed")
    require(result.get("source_checkpoint_sha256") == sha_bytes(source_raw)
            and result.get("entry_checkpoint_sha256") == sha_bytes(entered_raw)
            and result.get("diagnostic_checkpoint_sha256") == sha_bytes(journal_raw),
            "Saved result checkpoint changed")
    _validate_result(result, packet)
    return {
        "passed": result.get("passed") is True,
        "metadata_only": True,
        "component": COMPONENT,
        "source_sha256": sha_bytes(source_raw),
        "entered_sha256": sha_bytes(entered_raw),
        "result_sha256": sha_bytes(result_raw),
        "diagnostics_sha256": sha_bytes(journal_raw),
        "attempt_status": "passed" if result.get("passed") else "failed",
        "source_ftc_transport_defined": result.get("source_ftc_transport_defined"),
        "derive_entry_calls": result.get("derive_entry_calls"),
        "cached_formal_proof_calls": result.get("cached_formal_proof_calls"),
        "runtime_model_metadata": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "inspect"))
    args = parser.parse_args()
    print(json.dumps(globals()[args.action](), ensure_ascii=True,
                     allow_nan=False))


if __name__ == "__main__":
    main()
