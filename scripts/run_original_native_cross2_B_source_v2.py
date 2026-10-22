"""Freeze and consume one fresh cross2-B v2 source packet exactly once.

The prepare action authenticates the frozen failure report, all inherited source
pins, preserved failed artifacts, and old checkpoints before decoding the nine
selected source records.  The cross action records entry before importing only
the new v2 kernel, calls its single derive entry point once, and durably saves
the complete result and append-only diagnostics before validation.  Inspecting
a saved attempt reads metadata only.
"""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import gzip
import hashlib
import json
import os
import time
import traceback


ROOT = Path(__file__).resolve().parents[1]
RUNNER = "scripts/run_original_native_cross2_B_source_v2.py"
KERNEL = "src/openai_ns_reconstruction/original_native_cross2_B_source_v2.py"

PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_"
    "cross2_B_source_failure.json"
)
PARENT_SHA = "19d09adbfac2b236bcb1ff01712f6e6f1db20fe94c041284cb847e5f5916e60e"
PARENT_FAILURE_SCHEMA = "original_native_cross2_B_source_failure_stage_v1"
PARENT_FAILURE_COMPONENT = "cross2_B_source_v1"
PARENT_COMMIT = "8e08b206fd68f3cba97662d1bdaf3da02277b867"
PARENT_INTERNAL_COMMIT = "89d6f6d46fdcc6e83e9910e43bc544e606a00897"
PARENT_ANCESTOR_COMMIT = "89d6f6d46fdcc6e83e9910e43bc544e606a00897"
REPOSITORY_HEAD = "2b70546657ff89d445111d8aef0f4a24f669037c"

PREFIX = "lei_ren_part1_paper_compliant_current_original_"
C6_ARCHIVE = PREFIX + "prescribed_R100_Rm_C6_functions.json.gz"
LOOP_ARCHIVE = PREFIX + "prescribed_loop_C5_functions.json.gz"
PATCH_ARCHIVE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
OUTER_ARCHIVE = PREFIX + "true_outer_C6_C5_functions.json.gz"
C6_CREATOR = PREFIX + "prescribed_R100_Rm_C6_functions.py"
BRIDGE_CREATOR = PREFIX + "prescribed_bridge_C6_functions.py"
PATCH_CREATOR = PREFIX + "true_Rm_patch_finite_N_C5_functions.py"
OUTER_CREATOR = PREFIX + "true_outer_C6_C5_functions.py"

RECORD_RAW_HASHES = {
    PARENT: PARENT_SHA,
    C6_ARCHIVE: "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
    LOOP_ARCHIVE: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    PATCH_ARCHIVE: "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad",
    OUTER_ARCHIVE: "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d",
    C6_CREATOR: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    BRIDGE_CREATOR: "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    PATCH_CREATOR: "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    OUTER_CREATOR: "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
}
RECORD_DECODED_HASHES = {
    PARENT: PARENT_SHA,
    C6_ARCHIVE: "18d06627b935a46de5eca87c286a9a51c5a2d5ce14f11ae0c1c06a7c6308221a",
    LOOP_ARCHIVE: "f4825741f568c001f0f4658dfc3ce8eec3dc0466db976c29f2f4b609b8ceadc2",
    PATCH_ARCHIVE: "c1026e05015584a7ce29bd7e8da95b17d97669523f1d063ac5b5c800559d9b69",
    OUTER_ARCHIVE: "53b1bbf15c8c2c5a8400d7e22ef940c6a502397df8cf26bb23e51ff3869e9139",
    C6_CREATOR: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    BRIDGE_CREATOR: "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    PATCH_CREATOR: "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    OUTER_CREATOR: "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
}

CK = ROOT / "tmp/original_native_cross2_B_source_v2"
COMPONENT = "cross2_B_source_v2"
INPUT_SCHEMA = "native_cross2_B_source_inputs_v2"
RESULT_SCHEMA = "original_native_cross2_B_source_v2"
DIAGNOSTIC_SCHEMA = "native_cross2_B_source_stage_diagnostic_v2"
ENTERED_SCHEMA = "native_cross2_B_source_entered_v2"

INHERITED_SOURCE_COUNT = 1879
FINAL_SOURCE_COUNT = 1882
FAILED_ARTIFACT_COUNT = 38
AUTHENTICATED_OLD_CHECKPOINT_COUNT = 152
PARENT_CURRENT_CHECKPOINT_COUNT = 4
PARENT_INHERITED_CHECKPOINT_COUNT = 156
FAILURE_CURRENT_CHECKPOINT_COUNT = 4
INHERITED_CHECKPOINT_COUNT = 160
NEW_CHECKPOINT_COUNT = 4
FINAL_CHECKPOINT_COUNT = 164

TARGET_NODES = [529, 530]
PARENT_ACCEPTED_NODES = [520, 521, 522, 523, 524, 525, 526, 527, 528, 531, 532, 533, 534, 535]
ACCEPTED_NODES = list(range(520, 536))
ORDINARY_Z_ORDERS = [0, 1]
OWNED_Z_CELL = ["0", "1/2"]

EXPECTED_RECORD_NAMES = [
    PARENT, C6_ARCHIVE, LOOP_ARCHIVE, PATCH_ARCHIVE, OUTER_ARCHIVE,
    C6_CREATOR, BRIDGE_CREATOR, PATCH_CREATOR, OUTER_CREATOR,
]
EXPECTED_SOURCE_TEXT_ALIASES = {
    "c6": C6_CREATOR,
    "bridge": BRIDGE_CREATOR,
    "patch": PATCH_CREATOR,
    "outer": OUTER_CREATOR,
}

STAGE_FLAGS = [
    "actual_batch_source_ports_bound",
    "all_target_endpoint_identities_proved",
    "weighted_B_C1_source_composition_proved",
    "native_IN_batch_source_reduction_proved",
    "all_old_terms_retained",
    "accepted_ancestors_not_reexecuted",
]
FALSE_FLAGS = [
    "quantitative_packet_error_reduction_proved",
    "numeric_values_computed",
    "source_C5_installation_claimed",
    "numeric_containment_claimed",
    "new_loop_domain_admission_claimed",
    "admissible_cone_admission_claimed",
    "physical_layer_completed",
    "old_error_terms_removed",
    "global_five_moment_closure",
    "time_order_recursion_claimed",
    "full_goal_complete",
    "old_operators_replayed",
]

RAW_ONLY_MARKERS = (
    "resume_first_switch_2026_10_17",
    "first_switch_history",
    "true_first_switch_history",
    "rm_loop_density_identity",
)
FORBIDDEN_HISTORY_PATHS = (
    "docs/continuation/RESUME_FIRST_SWITCH_2026_10_17.json",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_loop_density_identity.json",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_Rm_loop_density_identity.json.gz",
)

def require(condition, message):
    if not condition:
        raise AssertionError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("utf-8")).hexdigest()


def source_path(name):
    candidate = ROOT / name
    return candidate if candidate.exists() else ROOT / "experiments/root_st073" / name


def _raw_only(name):
    normalized = name.replace("-", "_").lower()
    return any(marker in normalized for marker in RAW_ONLY_MARKERS)


def _resolve_literal(node, assignments):
    """Resolve only the small literal grammar used in the assembled header."""
    if isinstance(node, ast.Constant):
        return node.value
    if isinstance(node, ast.Name):
        require(node.id in assignments, "Unknown literal name: " + node.id)
        return _resolve_literal(assignments[node.id], assignments)
    if isinstance(node, (ast.List, ast.Tuple, ast.Set)):
        values = [_resolve_literal(item, assignments) for item in node.elts]
        return tuple(values) if isinstance(node, ast.Tuple) else values
    if isinstance(node, ast.Dict):
        return {
            _resolve_literal(key, assignments): _resolve_literal(value, assignments)
            for key, value in zip(node.keys, node.values)
        }
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Add):
        return _resolve_literal(node.left, assignments) + _resolve_literal(node.right, assignments)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        value = _resolve_literal(node.operand, assignments)
        return -value if isinstance(node.op, ast.USub) else value
    raise ValueError("Nonliteral kernel configuration: " + ast.dump(node))


def literal_configuration():
    """Read top-level literals without importing or executing the kernel."""
    tree = ast.parse((ROOT / KERNEL).read_text(encoding="utf-8"), filename=KERNEL)
    assignments = {}
    for statement in tree.body:
        if isinstance(statement, ast.Assign):
            for target in statement.targets:
                if isinstance(target, ast.Name):
                    assignments[target.id] = statement.value
        elif isinstance(statement, ast.AnnAssign) and isinstance(statement.target, ast.Name):
            assignments[statement.target.id] = statement.value
    names = (
        "SCHEMA", "RECORD_NAMES", "RECORD_HASHES", "DECODED_RECORD_HASHES",
        "SOURCE_TEXT_ALIASES", "TARGET_NODES", "STAGE_FLAGS", "FALSE_FLAGS",
    )
    return {name: _resolve_literal(assignments[name], assignments) for name in names}


def _read_record(name, expected, decoded_expected=None, decode=True):
    """Authenticate raw bytes before optional gzip or text decoding."""
    path = source_path(name)
    raw = path.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    require(actual == expected, "Changed raw source: " + name)
    record = {
        "raw_sha256": actual,
        "raw_authenticated_before_decode": True,
        "is_gzip": name.endswith(".gz"),
        "source_path": path.relative_to(ROOT).as_posix(),
        "raw_only": not decode,
    }
    if not decode:
        record["decoded_sha256"] = None
        return record
    decoded = gzip.decompress(raw) if name.endswith(".gz") else raw
    text = decoded.decode("utf-8")
    decoded_sha = hashlib.sha256(decoded).hexdigest()
    require(decoded_expected is None or decoded_sha == decoded_expected,
            "Changed decoded source: " + name)
    record.update(source_text=text, decoded_sha256=decoded_sha)
    return record


def _read_raw_hash(name, expected):
    """Authenticate a preserved artifact without decoding its bytes."""
    actual = hashlib.sha256(source_path(name).read_bytes()).hexdigest()
    require(actual == expected, "Changed preserved raw artifact: " + name)
    return actual


def _write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())


def paths():
    return tuple(CK / ("new_cross2_B_v2_" + kind + ".json")
                 for kind in ("source", "entered", "result"))


def diagnostic_path():
    return CK / "new_cross2_B_v2_diagnostics.jsonl"




def _failure_scope(failure):
    require(
        failure.get("schema") == PARENT_FAILURE_SCHEMA
        and failure.get("component") == PARENT_FAILURE_COMPONENT
        and failure.get("failure_frozen") is True
        and failure.get("source_only") is True
        and failure.get("do_not_reenter") is True
        and failure.get("do_not_reenter_any_copy") is True,
        "Frozen v1 failure report scope changed",
    )
    require(
        failure.get("accepted_identity") is False
        and failure.get("accepted_identity_flag") is False
        and failure.get("full_goal_complete") is False
        and failure.get("failure_attempt_count") == 1,
        "Failure report admission state changed",
    )
    require(
        failure.get("source_pin_count") == INHERITED_SOURCE_COUNT
        and failure.get("input_source_pin_count") == INHERITED_SOURCE_COUNT
        and failure.get("accepted_native_B_nodes") == PARENT_ACCEPTED_NODES
        and failure.get("remaining_B_identity_nodes") == TARGET_NODES
        and failure.get("remaining_B_identity_count") == len(TARGET_NODES)
        and failure.get("all_two_remaining_B_identities_open") is True,
        "Failure report source scope changed",
    )
    pins = failure.get("source_pins")
    require(
        isinstance(pins, dict)
        and len(pins) == INHERITED_SOURCE_COUNT
        and all(isinstance(name, str) and isinstance(value, str) and len(value) == 64
                for name, value in pins.items()),
        "Frozen inherited source pin map changed",
    )
    manifest = failure.get("source_pin_manifest", {})
    require(
        isinstance(manifest, dict)
        and manifest.get("inherited_source_pin_count") == INHERITED_SOURCE_COUNT - 3
        and manifest.get("new_pin_count") == 3
        and manifest.get("final_union_count") == INHERITED_SOURCE_COUNT
        and isinstance(manifest.get("final_union_hashes"), dict)
        and manifest.get("final_union_hashes") == pins,
        "Frozen source pin manifest changed",
    )
    component = failure.get("component_result", {})
    require(
        isinstance(component, dict)
        and component.get("schema") == "original_native_cross2_B_source_v1"
        and component.get("component") == PARENT_FAILURE_COMPONENT
        and component.get("passed") is False
        and component.get("source_only") is True
        and component.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
        and component.get("ordinary_Z_not_temporal_n") is True
        and component.get("owned_Z_cell") == OWNED_Z_CELL
        and component.get("numeric_values_computed") is False
        and component.get("old_operators_replayed") is False
        and component.get("remaining_B_identity_nodes") == TARGET_NODES
        and component.get("remaining_B_identity_count") == len(TARGET_NODES),
        "Frozen failed component scope changed",
    )
    require(
        failure.get("component_result_exact") is True
        and failure.get("component_result_is_projection") is False
        and failure.get("full_failed_result_exact_saved_json") is True
        and isinstance(failure.get("full_failed_result"), dict),
        "Frozen full failure result is not exact",
    )
    stage = component.get("stage_flags")
    require(
        isinstance(stage, dict)
        and stage.get("actual_batch_source_ports_bound") is False
        and stage.get("all_target_endpoint_identities_proved") is False
        and stage.get("weighted_B_C1_source_composition_proved") is False
        and stage.get("native_IN_batch_source_reduction_proved") is False
        and stage.get("all_old_terms_retained") is True
        and stage.get("accepted_ancestors_not_reexecuted") is True,
        "Frozen failed stage flags changed",
    )
    for flag in FALSE_FLAGS:
        require(component.get(flag) is False, "Frozen failed false flag changed: " + flag)
    retained = component.get("retained_construction_data", {})
    require(
        retained.get("cumulative_terms") == 212
        and retained.get("raw_signed_terms") == 213
        and retained.get("reduced_terms") == 183
        and retained.get("six_signed_channels_retained") is True
        and retained.get("all_old_error_data_retained") is True
        and retained.get("reserved_contribution_registered_once") is True
        and retained.get("inherited_errors_deleted") is False,
        "Frozen construction ledger changed",
    )
    integrals = component.get("retained_integral_term_inventory", {})
    require(
        integrals.get("K_terms") == 17
        and integrals.get("outer_B_terms") == 2
        and integrals.get("K_terms_retained") is True
        and integrals.get("outer_B_terms_retained") is True,
        "Frozen integral-term inventory changed",
    )
    errors = component.get("old_errors_frozen", {})
    require(
        errors.get("cumulative_terms") == 212
        and errors.get("raw_signed_terms") == 213
        and errors.get("reduced_terms") == 183
        and errors.get("all_old_packet_errors_retained") is True,
        "Frozen error ledger changed",
    )
    audit = failure.get("source_inspection_audit", {})
    require(
        isinstance(audit, dict)
        and audit.get("new_source_consumer_calls") == 1
        and audit.get("accepted_source_consumer_calls") == 0
        and audit.get("accepted_ancestor_operator_calls") == 0
        and audit.get("old_proof_functions_called") == 0
        and audit.get("numerical_evaluations") == 0
        and audit.get("quadratures") == 0
        and audit.get("recurrence_calls") == 0
        and audit.get("runtime_model_metadata") is None,
        "Frozen failure source audit changed",
    )
    failed = failure.get("preserved_failed_artifact_hashes_union")
    require(
        isinstance(failed, dict)
        and len(failed) == FAILED_ARTIFACT_COUNT
        and failure.get("preserved_failed_artifact_count_union") == FAILED_ARTIFACT_COUNT
        and all(isinstance(name, str) and isinstance(value, str) and len(value) == 64
                for name, value in failed.items()),
        "Frozen failed-artifact union changed",
    )
    return pins, failed


def _failure_checkpoint_inventory(failure):
    old = failure.get("old_checkpoint_inventory")
    current = failure.get("checkpoint_hashes")
    require(isinstance(old, dict) and isinstance(current, dict),
            "Frozen checkpoint inventory missing")
    authenticated = old.get("authenticated_checkpoint_hashes")
    parent_current = old.get("parent_current_checkpoint_hashes")
    require(
        isinstance(authenticated, dict)
        and len(authenticated) == AUTHENTICATED_OLD_CHECKPOINT_COUNT
        and isinstance(parent_current, dict)
        and len(parent_current) == PARENT_CURRENT_CHECKPOINT_COUNT
        and len(current) == FAILURE_CURRENT_CHECKPOINT_COUNT,
        "Frozen checkpoint counts changed",
    )
    require(
        old.get("authenticated_old_count") == PARENT_INHERITED_CHECKPOINT_COUNT
        and old.get("authenticated_total_count") == PARENT_INHERITED_CHECKPOINT_COUNT
        and old.get("parent_current_count") == PARENT_CURRENT_CHECKPOINT_COUNT
        and old.get("distinct_parent_count") == PARENT_INHERITED_CHECKPOINT_COUNT
        and old.get("raw_sha_only") is True
        and old.get("history_decode_forbidden") is True,
        "Frozen checkpoint metadata changed",
    )
    combined = {**authenticated, **parent_current, **current}
    require(
        len(combined) == INHERITED_CHECKPOINT_COUNT
        and all(isinstance(name, str) and isinstance(value, str) and len(value) == 64
                for name, value in combined.items()),
        "Frozen checkpoint hash union changed",
    )
    return {
        "authenticated_checkpoint_hashes": dict(authenticated),
        "parent_current_checkpoint_hashes": dict(parent_current),
        "failure_current_checkpoint_hashes": dict(current),
        "authenticated_old_count": AUTHENTICATED_OLD_CHECKPOINT_COUNT,
        "parent_current_count": PARENT_CURRENT_CHECKPOINT_COUNT,
        "failure_current_count": FAILURE_CURRENT_CHECKPOINT_COUNT,
        "inherited_count": INHERITED_CHECKPOINT_COUNT,
        "raw_sha_only": True,
        "history_decode_forbidden": True,
    }


def _static_scope(config):
    require(config["SCHEMA"] == RESULT_SCHEMA, "Kernel schema changed")
    names = list(config["RECORD_NAMES"])
    hashes = dict(config["RECORD_HASHES"])
    decoded = dict(config["DECODED_RECORD_HASHES"])
    require(names == EXPECTED_RECORD_NAMES
            and names == list(dict.fromkeys(names)),
            "Kernel record inventory changed")
    require(set(hashes) == set(names)
            and set(decoded) == set(names)
            and hashes == RECORD_RAW_HASHES
            and decoded == RECORD_DECODED_HASHES,
            "Kernel record hashes changed")
    require(all(isinstance(value, str) and len(value) == 64
                for value in hashes.values())
            and all(isinstance(value, str) and len(value) == 64
                    for value in decoded.values()),
            "Kernel record hashes are malformed")
    require(config["TARGET_NODES"] == TARGET_NODES,
            "Kernel target nodes changed")
    stage = list(config["STAGE_FLAGS"])
    false = list(config["FALSE_FLAGS"])
    require(stage == STAGE_FLAGS and false == FALSE_FLAGS,
            "Kernel flag ledger changed")
    aliases = dict(config["SOURCE_TEXT_ALIASES"])
    require(aliases == EXPECTED_SOURCE_TEXT_ALIASES,
            "Kernel creator aliases changed")
    return names, hashes, decoded, aliases, stage, false


def _prepare_packet():
    config = literal_configuration()
    names, expected_hashes, decoded_expected, aliases, stage, false = _static_scope(config)
    parent_raw = source_path(PARENT).read_bytes()
    require(hashlib.sha256(parent_raw).hexdigest() == PARENT_SHA,
            "Changed frozen cross2 failure report")
    failure = json.loads(parent_raw)
    inherited, failed = _failure_scope(failure)
    require(PARENT not in inherited and RUNNER not in inherited and KERNEL not in inherited,
            "Fresh v2 publication is already in inherited pins")

    for name, expected in inherited.items():
        _read_raw_hash(name, expected)
    for name, expected in failed.items():
        _read_raw_hash(name, expected)
    checkpoint_inventory = _failure_checkpoint_inventory(failure)
    for group in (
        "authenticated_checkpoint_hashes",
        "parent_current_checkpoint_hashes",
        "failure_current_checkpoint_hashes",
    ):
        for name, expected in checkpoint_inventory[group].items():
            _read_raw_hash(name, expected)

    current_code_hashes = {
        RUNNER: sha(ROOT / RUNNER),
        KERNEL: sha(ROOT / KERNEL),
    }
    final_union = dict(inherited)
    final_union[PARENT] = PARENT_SHA
    final_union.update(current_code_hashes)
    require(len(final_union) == FINAL_SOURCE_COUNT,
            "Final source union must be 1882 pins")
    require(names == EXPECTED_RECORD_NAMES,
            "Fresh v2 record ordering changed")
    require(set(final_union) >= set(names),
            "Selected record escapes authenticated source union")
    for name, expected in expected_hashes.items():
        require(final_union.get(name) == expected,
                "Selected record hash changed: " + name)

    records = {}
    for name in names:
        records[name] = _read_record(
            name, expected_hashes[name], decoded_expected[name],
            decode=not _raw_only(name))
    input_digest = digest(records)
    record_decoded = {name: item["decoded_sha256"] for name, item in records.items()}
    raw_only = [name for name, item in records.items() if item.get("raw_only")]
    static_receipt = {
        "metadata_only": True,
        "actions_executed": False,
        "AST_receipt_source_calls": 0,
        "source_functions_executed": 0,
        "kernel_imported": False,
        "derive_calls": 0,
        "old_proof_calls": 0,
        "old_operator_calls": 0,
        "numeric_values_computed": False,
        "runtime_model_metadata": None,
        "requested_worker_model": "gpt-5.6-luna",
        "requested_reasoning_effort": "max",
    }
    return dict(
        schema=INPUT_SCHEMA,
        result_schema=RESULT_SCHEMA,
        component=COMPONENT,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA,
        parent_failure_schema=PARENT_FAILURE_SCHEMA,
        parent_commit=PARENT_COMMIT,
        parent_internal_commit=PARENT_INTERNAL_COMMIT,
        parent_ancestor_commit=PARENT_ANCESTOR_COMMIT,
        repository_head=REPOSITORY_HEAD,
        input_digest=input_digest,
        record_count=len(names),
        record_names=names,
        record_hashes=expected_hashes,
        record_raw_sha256={name: item["raw_sha256"] for name, item in records.items()},
        record_decoded_sha256=record_decoded,
        decoded_record_hashes=record_decoded,
        record_source_paths={name: item["source_path"] for name, item in records.items()},
        raw_only_records=raw_only,
        inherited_source_hashes=dict(inherited),
        inherited_source_pin_count=INHERITED_SOURCE_COUNT,
        current_code_hashes=current_code_hashes,
        implementation_hashes=current_code_hashes,
        source_pin_manifest={
            "inherited_source_hashes": dict(inherited),
            "inherited_source_pin_count": INHERITED_SOURCE_COUNT,
            "parent_report_hashes": {PARENT: PARENT_SHA},
            "parent_commit": PARENT_COMMIT,
            "current_code_hashes": dict(current_code_hashes),
            "record_raw_sha256": {
                name: item["raw_sha256"] for name, item in records.items()},
            "record_decoded_sha256": dict(record_decoded),
            "input_digest": input_digest,
            "new_pins": [PARENT, RUNNER, KERNEL],
            "new_pin_count": 3,
            "final_union_hashes": final_union,
            "final_union_count": FINAL_SOURCE_COUNT,
            "parent_report_separate_publication": True,
        },
        preserved_failed_artifact_hashes=dict(failed),
        preserved_failed_artifact_count=FAILED_ARTIFACT_COUNT,
        parent_checkpoint_inventory=checkpoint_inventory,
        new_checkpoint_count=NEW_CHECKPOINT_COUNT,
        final_checkpoint_count=FINAL_CHECKPOINT_COUNT,
        schema_expected=config["SCHEMA"],
        stage_flags=stage,
        false_flags=false,
        target_nodes=TARGET_NODES,
        accepted_nodes_if_pass=ACCEPTED_NODES,
        accepted_nodes_if_failed=PARENT_ACCEPTED_NODES,
        remaining_nodes_if_failed=TARGET_NODES,
        ordinary_Z_orders=ORDINARY_Z_ORDERS,
        ordinary_Z_not_temporal_n=True,
        owned_Z_cell=OWNED_Z_CELL,
        physical_scope="m through Z2 on original owned cell 0 and 1/2",
        Q_y_identity_claim=False,
        binding_fix_count=3,
        endpoint_case_nodes=TARGET_NODES,
        source_only=True,
        ancestor_operators_replayed=False,
        old_producer_function_calls=0,
        old_operator_calls=0,
        old_proof_functions_called=0,
        numeric_values_computed=False,
        full_goal_complete=False,
        retained_construction_data={
            "cumulative_terms": 212,
            "raw_signed_terms": 213,
            "reduced_terms": 183,
            "six_signed_channels_retained": True,
            "all_old_error_data_retained": True,
            "reserved_contribution_registered_once": True,
            "inherited_errors_deleted": False,
        },
        retained_integral_term_inventory={
            "K_terms": 17, "outer_B_terms": 2,
            "K_terms_retained": True, "outer_B_terms_retained": True,
        },
        old_errors_frozen={
            "cumulative_terms": 212, "raw_signed_terms": 213,
            "reduced_terms": 183, "all_old_packet_errors_retained": True,
        },
        static_receipt=static_receipt,
        diagnostics_policy="append-only JSONL; snapshot, flush, and fsync before callback return",
        diagnostics_checkpoint_additional_to_source_entered_result=True,
        source_functions_executed=0,
        runtime_model_metadata=None,
        requested_worker_model="gpt-5.6-luna",
        requested_reasoning_effort="max",
        raw_authentication_policy={
            "all_inherited_source_pins_before_decode": True,
            "all_failed_artifacts_before_decode": True,
            "all_old_checkpoints_before_decode": True,
            "forbidden_history_raw_only": True,
            "all_gzip_raw_sha_before_decode": True,
            "selected_record_count": len(names),
        },
        forbidden_history_paths=list(FORBIDDEN_HISTORY_PATHS),
        no_cached_record_decode=True,
    )


def prepare_cross():
    source, entered, result = paths()
    journal = diagnostic_path()
    require(not any(path.exists() for path in (source, entered, result, journal)),
            "Cross2 v2 checkpoints already exist; never rerun")
    packet = _prepare_packet()
    _write_new(source, packet)
    return {
        "passed": True,
        "prepared_only": True,
        "ready_hash": sha(source),
        "source_sha256": sha(source),
        "code_sha256": dict(packet["current_code_hashes"]),
        "input_digest": packet["input_digest"],
        "source_pin_count": packet["source_pin_manifest"]["final_union_count"],
        "record_count": len(packet["record_names"]),
        "api": ["prepare_cross", "cross", "inspect_cross"],
        "static_receipt": packet["static_receipt"],
    }


def _load_packet(source):
    packet = json.loads(source.read_bytes())
    require(
        packet.get("schema") == INPUT_SCHEMA
        and packet.get("result_schema") == RESULT_SCHEMA
        and packet.get("component") == COMPONENT,
        "Frozen v2 source packet schema changed",
    )
    manifest = packet.get("source_pin_manifest", {})
    require(
        manifest.get("final_union_count") == FINAL_SOURCE_COUNT
        and packet.get("record_count") == len(EXPECTED_RECORD_NAMES)
        and packet.get("record_names") == EXPECTED_RECORD_NAMES
        and packet.get("inherited_source_pin_count") == INHERITED_SOURCE_COUNT,
        "Frozen v2 source union or record count changed",
    )
    require(
        packet.get("new_checkpoint_count") == NEW_CHECKPOINT_COUNT
        and packet.get("final_checkpoint_count") == FINAL_CHECKPOINT_COUNT
        and packet.get("parent_report_sha256") == PARENT_SHA,
        "Frozen v2 checkpoint or parent hash changed",
    )
    inherited = packet.get("inherited_source_hashes", {})
    require(
        isinstance(inherited, dict)
        and len(inherited) == INHERITED_SOURCE_COUNT
        and inherited == manifest.get("inherited_source_hashes"),
        "Frozen inherited source map changed",
    )
    require(
        set(packet.get("current_code_hashes", {})) == {RUNNER, KERNEL}
        and packet.get("implementation_hashes") == packet["current_code_hashes"],
        "Frozen v2 implementation inventory changed",
    )
    inventory = packet.get("parent_checkpoint_inventory", {})
    require(
        isinstance(inventory, dict)
        and inventory.get("inherited_count") == INHERITED_CHECKPOINT_COUNT
        and len(inventory.get("authenticated_checkpoint_hashes", {})) == AUTHENTICATED_OLD_CHECKPOINT_COUNT
        and len(inventory.get("parent_current_checkpoint_hashes", {})) == PARENT_CURRENT_CHECKPOINT_COUNT
        and len(inventory.get("failure_current_checkpoint_hashes", {})) == FAILURE_CURRENT_CHECKPOINT_COUNT,
        "Frozen v2 checkpoint inventory changed",
    )
    for name, expected in packet["current_code_hashes"].items():
        require(sha(ROOT / name) == expected,
                "Current implementation changed after prepare: " + name)
    return packet


def _load_frozen_records(packet):
    parent_path = source_path(packet["parent_report"])
    parent_raw = parent_path.read_bytes()
    require(hashlib.sha256(parent_raw).hexdigest() == packet["parent_report_sha256"],
            "Changed frozen parent failure report")
    failure = json.loads(parent_raw)
    inherited, failed = _failure_scope(failure)
    require(inherited == packet["inherited_source_hashes"],
            "Packet inherited pins do not match frozen failure report")
    require(len(inherited) == INHERITED_SOURCE_COUNT,
            "Inherited source pin count changed")
    for name, expected in inherited.items():
        _read_raw_hash(name, expected)
    require(failed == packet["preserved_failed_artifact_hashes"],
            "Packet failed-artifact union does not match failure report")
    for name, expected in failed.items():
        _read_raw_hash(name, expected)
    inventory = packet["parent_checkpoint_inventory"]
    for group in (
        "authenticated_checkpoint_hashes",
        "parent_current_checkpoint_hashes",
        "failure_current_checkpoint_hashes",
    ):
        for name, expected in inventory[group].items():
            _read_raw_hash(name, expected)
    raw_only = set(packet.get("raw_only_records", []))
    require(
        raw_only == {name for name in packet["record_names"] if _raw_only(name)}
        and not any(name in packet["record_names"] for name in FORBIDDEN_HISTORY_PATHS),
        "Frozen raw-only source policy changed",
    )
    records = {}
    for name in packet["record_names"]:
        records[name] = _read_record(
            name, packet["record_hashes"][name],
            packet["record_decoded_sha256"][name],
            decode=name not in raw_only,
        )
    require(digest(records) == packet["input_digest"],
            "Frozen v2 input digest changed")
    return records


def _validate_ledger(result):
    retained = result.get("retained_construction_data")
    require(
        isinstance(retained, dict)
        and retained.get("cumulative_terms") == 212
        and retained.get("raw_signed_terms") == 213
        and retained.get("reduced_terms") == 183
        and retained.get("six_signed_channels_retained") is True
        and retained.get("all_old_error_data_retained") is True
        and retained.get("reserved_contribution_registered_once") is True
        and retained.get("inherited_errors_deleted") is False,
        "Construction ledger changed",
    )
    integrals = result.get("retained_integral_term_inventory")
    require(
        isinstance(integrals, dict)
        and integrals.get("K_terms") == 17
        and integrals.get("outer_B_terms") == 2
        and integrals.get("K_terms_retained") is True
        and integrals.get("outer_B_terms_retained") is True,
        "Integral-term inventory changed",
    )
    errors = result.get("old_errors_frozen")
    require(
        isinstance(errors, dict)
        and errors.get("cumulative_terms") == 212
        and errors.get("raw_signed_terms") == 213
        and errors.get("reduced_terms") == 183
        and errors.get("all_old_packet_errors_retained") is True,
        "Frozen error ledger changed",
    )


def _validate_result(result, packet):
    require(
        isinstance(result, dict)
        and result.get("schema") == RESULT_SCHEMA
        and result.get("component") == COMPONENT
        and result.get("source_only") is True
        and type(result.get("passed")) is bool,
        "Kernel result schema or component changed",
    )
    require(
        result.get("input_digest") == packet["input_digest"]
        and result.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
        and result.get("ordinary_Z_not_temporal_n") is True
        and result.get("owned_Z_cell") == OWNED_Z_CELL
        and result.get("target_B_identity_nodes") == TARGET_NODES,
        "Kernel result scope changed",
    )
    for flag in packet["false_flags"]:
        require(result.get(flag) is False, "Unadmitted scope changed: " + flag)
    stage_values = [result.get(flag) for flag in packet["stage_flags"]]
    require(all(type(value) is bool for value in stage_values),
            "Kernel stage flags are not booleans")
    declared = result.get("stage_flags")
    if declared is not None:
        require(
            isinstance(declared, dict)
            and all(declared.get(flag) is result.get(flag)
                    for flag in packet["stage_flags"]),
            "Kernel declared stage flags changed",
        )
    all_stages = all(stage_values)
    require((not all_stages) or result["passed"],
            "All declared kernel stages cannot be true on a failed result")
    if result["passed"]:
        require(
            all_stages
            and result.get("newly_proved_B_identity_nodes") == TARGET_NODES
            and result.get("accepted_native_B_nodes") == ACCEPTED_NODES
            and result.get("remaining_B_identity_nodes") == []
            and result.get("remaining_B_identity_count") == 0
            and result.get("all_native_B_jumps_zero_proved") is True,
            "Cross2 v2 acceptance set changed",
        )
    else:
        require(
            result.get("accepted_native_B_nodes") == PARENT_ACCEPTED_NODES
            and result.get("remaining_B_identity_nodes") == TARGET_NODES
            and result.get("remaining_B_identity_count") == len(TARGET_NODES)
            and result.get("all_native_B_jumps_zero_proved") is False,
            "Failed v2 result must keep both target identities open",
        )
    proofs = result.get("component_proofs")
    require(isinstance(proofs, dict), "Component proof metadata missing")
    fresh = proofs.get("fresh_cross")
    require(isinstance(fresh, dict)
            and type(fresh.get("passed")) is bool
            and isinstance(fresh.get("per_node"), dict)
            and set(fresh["per_node"]) == {str(node) for node in TARGET_NODES},
            "Fresh cross proof metadata changed")
    require(
        fresh.get("passed") is result["passed"]
        and all(isinstance(fresh["per_node"][str(node)], dict)
                for node in TARGET_NODES),
        "Fresh cross proof scope changed",
    )
    if "additional_source_contracts" in fresh:
        require(isinstance(fresh["additional_source_contracts"], dict),
                "Fresh source contracts changed")
    if "cached_v3_source_evidence" in proofs:
        require(isinstance(proofs["cached_v3_source_evidence"], dict),
                "Cached evidence shape changed")
    for key in (
        "old_proof_functions_called", "old_operator_calls",
        "old_operator_function_calls", "source_producer_function_calls",
        "old_producer_function_calls", "new_cross_calls",
        "new_source_consumer_calls", "numerical_evaluations",
        "quadratures", "recurrence_calls",
    ):
        if key in result:
            require(
                result[key] == 0
                or (key in {"new_cross_calls", "new_source_consumer_calls"}
                    and result[key] == 1),
                "Unexpected call count: " + key,
            )
    require(
        result.get("old_operators_replayed", False) is False
        and result.get("numeric_values_computed") is False
        and result.get("runtime_model_metadata") is None,
        "Unadmitted replay, numeric evaluation, or model metadata",
    )
    _validate_ledger(result)
    audit = result.get("source_inspection_audit")
    if isinstance(audit, dict):
        require(
            audit.get("new_source_consumer_calls") == 1
            and audit.get("accepted_source_consumer_calls") == 0
            and audit.get("accepted_ancestor_operator_calls") == 0
            and audit.get("old_proof_functions_called", 0) == 0
            and audit.get("runtime_model_metadata") is None,
            "Source consumer audit changed",
        )
    require(
        result.get("native_IN_reduction_is_source_only") is True
        and result.get("archived_all_16_J_terms_retained") is True
        and result.get("source_jump_zero_does_not_delete_quantitative_errors") is True,
        "Inherited source construction scope changed",
    )


def _empty_journal(path):
    with path.open("xb") as stream:
        stream.flush()
        os.fsync(stream.fileno())


def cross():
    source, entered, target = paths()
    journal = diagnostic_path()
    require(
        source.exists() and not entered.exists()
        and not target.exists() and not journal.exists(),
        "Frozen v2 packet is absent or already entered",
    )
    packet = _load_packet(source)
    _write_new(entered, {
        "schema": ENTERED_SCHEMA,
        "component": COMPONENT,
        "source_sha256": sha(source),
        "input_digest": packet["input_digest"],
        "implementation_hashes": packet["current_code_hashes"],
        "kernel_sha256": packet["current_code_hashes"][KERNEL],
        "inherited_source_pin_count": INHERITED_SOURCE_COUNT,
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_imported_before_entry": False,
        "derive_calls_before_entry": 0,
        "entered_code_immutable": True,
        "never_retry_entered_code": True,
    })
    _empty_journal(journal)
    events = []

    def persist_diagnostic(stage, payload):
        require(isinstance(stage, str) and bool(stage), "Invalid diagnostic stage")
        frozen = json.loads(json.dumps(
            payload, sort_keys=True, ensure_ascii=True, allow_nan=False))
        event = {
            "schema": DIAGNOSTIC_SCHEMA,
            "index": len(events),
            "stage": stage,
            "input_digest": packet["input_digest"],
            "payload": frozen,
        }
        raw = (json.dumps(event, sort_keys=True, ensure_ascii=True,
                          allow_nan=False) + "\n").encode("utf-8")
        with journal.open("ab") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        events.append(event)

    started = time.perf_counter()
    caught = False
    records = None
    try:
        records = _load_frozen_records(packet)
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "new_native_cross2_B_source_v2", ROOT / KERNEL)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        result = module.derive(
            records,
            expected_input_digest=packet["input_digest"],
            authenticated_inherited_source_hashes=packet["inherited_source_hashes"],
            diagnostic_callback=persist_diagnostic,
        )
    except Exception as error:
        caught = True
        result = {
            "schema": RESULT_SCHEMA,
            "component": COMPONENT,
            "passed": False,
            "source_only": True,
            "failed_attempt_preserved": True,
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": "".join(traceback.format_exception(error)),
            "input_digest": packet["input_digest"],
            "ordinary_Z_orders": ORDINARY_Z_ORDERS,
            "ordinary_Z_not_temporal_n": True,
            "owned_Z_cell": OWNED_Z_CELL,
            "target_B_identity_nodes": TARGET_NODES,
            "accepted_native_B_nodes": PARENT_ACCEPTED_NODES,
            "remaining_B_identity_nodes": TARGET_NODES,
            "remaining_B_identity_count": len(TARGET_NODES),
            "all_native_B_jumps_zero_proved": False,
            "runtime_model_metadata": None,
        }
    if not isinstance(result, dict):
        caught = True
        result = {
            "schema": RESULT_SCHEMA,
            "component": COMPONENT,
            "passed": False,
            "source_only": True,
            "failed_attempt_preserved": True,
            "error_type": "TypeError",
            "error": "derive returned a non-mapping result",
            "traceback": "",
            "input_digest": packet["input_digest"],
            "target_B_identity_nodes": TARGET_NODES,
            "accepted_native_B_nodes": PARENT_ACCEPTED_NODES,
            "remaining_B_identity_nodes": TARGET_NODES,
            "remaining_B_identity_count": len(TARGET_NODES),
            "all_native_B_jumps_zero_proved": False,
            "runtime_model_metadata": None,
        }
    result.update(
        execution_seconds=time.perf_counter() - started,
        source_checkpoint_sha256=sha(source),
        entry_checkpoint_sha256=sha(entered),
        runner_caught_exception=caught,
        runner_full_goal_complete=False,
        runtime_model_metadata=None,
        authenticated_inherited_source_pin_count=INHERITED_SOURCE_COUNT,
        authenticated_inherited_source_hashes_sha256=digest(
            packet["inherited_source_hashes"]),
        stage_diagnostics=list(events),
        diagnostic_events=list(events),
        diagnostic_checkpoint_path=journal.relative_to(ROOT).as_posix(),
        diagnostic_checkpoint_sha256=sha(journal),
        diagnostic_events_saved_before_validation=len(events),
        full_result_saved_before_validation=True,
        entered_code_immutable=True,
        never_retry_entered_code=True,
    )
    _write_new(target, result)
    if records is not None:
        require(digest(records) == packet["input_digest"],
                "Input mutated during v2 source consumer")
    if not caught:
        _validate_result(result, packet)
    return {
        "passed": result["passed"],
        "ready_hash": sha(target),
        "result_sha256": sha(target),
        "code_sha256": dict(packet["current_code_hashes"]),
        "seconds": result["execution_seconds"],
        "error": result.get("error"),
        "diagnostic_events": len(events),
        "api": ["prepare_cross", "cross", "inspect_cross"],
    }


def inspect_cross():
    """Inspect saved metadata only; never decode selected inputs or import code."""
    source, entered, target = paths()
    journal = diagnostic_path()
    require(all(path.exists() for path in (source, entered, target, journal)),
            "Saved v2 checkpoints are incomplete")
    packet = json.loads(source.read_bytes())
    entry = json.loads(entered.read_bytes())
    result = json.loads(target.read_bytes())
    require(
        packet.get("schema") == INPUT_SCHEMA
        and packet.get("result_schema") == RESULT_SCHEMA
        and packet.get("component") == COMPONENT
        and packet.get("source_pin_manifest", {}).get("final_union_count")
        == FINAL_SOURCE_COUNT
        and packet.get("record_names") == EXPECTED_RECORD_NAMES
        and packet.get("record_count") == len(EXPECTED_RECORD_NAMES),
        "Saved v2 packet metadata changed",
    )
    inherited = packet.get("inherited_source_hashes", {})
    require(isinstance(inherited, dict)
            and len(inherited) == INHERITED_SOURCE_COUNT,
            "Saved inherited source map changed")
    for name, expected in packet.get("current_code_hashes", {}).items():
        require(sha(ROOT / name) == expected,
                "Saved implementation changed: " + name)
    require(
        entry.get("schema") == ENTERED_SCHEMA
        and entry.get("source_sha256") == sha(source)
        and entry.get("input_digest") == packet.get("input_digest")
        and entry.get("implementation_hashes") == packet.get("current_code_hashes")
        and entry.get("kernel_imported_before_entry") is False
        and entry.get("derive_calls_before_entry") == 0
        and entry.get("entered_code_immutable") is True,
        "Saved entry metadata changed",
    )
    require(
        result.get("source_checkpoint_sha256") == sha(source)
        and result.get("entry_checkpoint_sha256") == sha(entered)
        and result.get("input_digest") == packet.get("input_digest")
        and result.get("schema") == RESULT_SCHEMA
        and result.get("component") == COMPONENT
        and result.get("source_only") is True
        and type(result.get("passed")) is bool
        and result.get("diagnostic_checkpoint_path")
        == journal.relative_to(ROOT).as_posix()
        and result.get("diagnostic_checkpoint_sha256") == sha(journal)
        and result.get("full_result_saved_before_validation") is True
        and result.get("authenticated_inherited_source_pin_count")
        == INHERITED_SOURCE_COUNT,
        "Saved result checkpoint links changed",
    )
    lines = journal.read_text(encoding="utf-8").splitlines()
    events = [json.loads(line) for line in lines]
    require(
        all(event.get("schema") == DIAGNOSTIC_SCHEMA
            and event.get("input_digest") == packet.get("input_digest")
            and event.get("index") == index
            for index, event in enumerate(events)),
        "Saved diagnostic journal metadata changed",
    )
    require(
        result.get("stage_diagnostics") == events
        and result.get("diagnostic_events") == events
        and result.get("diagnostic_events_saved_before_validation") == len(events),
        "Saved diagnostic/result snapshots changed",
    )
    if result.get("runner_caught_exception"):
        require(
            result.get("passed") is False
            and result.get("failed_attempt_preserved") is True
            and isinstance(result.get("traceback"), str),
            "Saved exception envelope changed",
        )
    return {
        "passed": True,
        "saved_only": True,
        "ready_hash": sha(target),
        "component_passed": result.get("passed"),
        "full_goal_complete": False,
        "code_sha256": packet.get("current_code_hashes"),
        "static_receipt": packet.get("static_receipt"),
        "diagnostic_events": len(events),
        "api": ["prepare_cross", "cross", "inspect_cross"],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare_cross", "cross", "inspect_cross"))
    action = parser.parse_args().action
    output = prepare_cross() if action == "prepare_cross" else (
        cross() if action == "cross" else inspect_cross())
    print(json.dumps(output, ensure_ascii=True, allow_nan=False))
