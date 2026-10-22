"""Freeze and consume one cross2-B source packet exactly once.

The prepare action reads only literal kernel configuration and authenticated
metadata.  The cross action writes entry state before importing the assembled
kernel, invokes its single ``derive`` entry point once, and writes the complete
result before validation.  Inspecting a saved attempt checks metadata and the
append-only diagnostic journal only.
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
RUNNER = "scripts/run_original_native_cross2_B_source_v1.py"
KERNEL = "src/openai_ns_reconstruction/original_native_cross2_B_source_v1.py"
PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_inner_B_v1.json"
)
PARENT_SHA = "fbc4aefecf54530183f9da696da38dc98b4201793c50f4d3f0be3afe8682647b"
PARENT_COMMIT = "8e08b206fd68f3cba97662d1bdaf3da02277b867"
PARENT_INTERNAL_COMMIT = "89d6f6d46fdcc6e83e9910e43bc544e606a00897"
PARENT_ANCESTOR_COMMIT = "89d6f6d46fdcc6e83e9910e43bc544e606a00897"
PARENT_INNER_RUNNER = "scripts/run_original_native_inner_B_source_v1.py"
PARENT_INNER_RUNNER_SHA = "a93b790b97672aaa2369884228b5d39a2c680ed793269efd2b1120a2dc3a258f"
PARENT_INNER_KERNEL = "src/openai_ns_reconstruction/original_native_inner_B_source_v1.py"
PARENT_INNER_KERNEL_SHA = "09d9bc12832dad13535dbb1a51127bf25e9a6dd19f5e5d4799107b467e66ec8d"
PRIOR_J3 = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_"
    "J3_endpoint_transport_v3.json"
)
PRIOR_J3_SHA = "c798cc171350236483e4e5ee4e35ba2bcc2dc4d6b714a491d926eb149175f7d8"
WEIGHTED = (
    "lei_ren_part1_paper_compliant_current_original_native_"
    "active_S1_weighted_ibp.json"
)
BACKGROUND = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_"
    "R100_Rm_C6_functions.json.gz"
)
LOOP = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_"
    "loop_C5_functions.json.gz"
)
MOD = (
    "lei_ren_part1_paper_compliant_current_original_modified_"
    "Q4_pressure_C3_functions.json.gz"
)
PATCH = (
    "lei_ren_part1_paper_compliant_current_original_true_"
    "Rm_patch_finite_N_C5_functions.json.gz"
)
OUTER = (
    "lei_ren_part1_paper_compliant_current_original_true_"
    "outer_C6_C5_functions.json.gz"
)
C6_CREATOR = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_"
    "R100_Rm_C6_functions.py"
)
BRIDGE_CREATOR = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_"
    "bridge_C6_functions.py"
)
LOOP_CREATOR = (
    "lei_ren_part1_paper_compliant_current_original_prescribed_"
    "loop_C5_functions.py"
)
MOD_CREATOR = (
    "lei_ren_part1_paper_compliant_current_original_modified_"
    "Q4_pressure_C3_functions.py"
)
PATCH_CREATOR = (
    "lei_ren_part1_paper_compliant_current_original_true_"
    "Rm_patch_finite_N_C5_functions.py"
)
OUTER_CREATOR = (
    "lei_ren_part1_paper_compliant_current_original_true_"
    "outer_C6_C5_functions.py"
)
GAMMA_CREATOR = (
    "lei_ren_part1_paper_compliant_current_original_n1_"
    "Gamma_future_functions.py"
)

RECORD_RAW_HASHES = {
    PRIOR_J3: PRIOR_J3_SHA,
    WEIGHTED: "aab72bbc7a882e9bf089a43dd3407babe085f1c0ec61583578d20a5e7c9bec28",
    BACKGROUND: "5e5b565dd8d2abc41e6f3df0d9e17000e524d493e58a5778dfcf217c6e42f44a",
    LOOP: "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa",
    MOD: "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a",
    PATCH: "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad",
    OUTER: "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d",
    C6_CREATOR: "70dd6d6022f68d601ac4f0fbaf38614bbaf215baab272f134f1c5a159af65d32",
    BRIDGE_CREATOR: "8ab2c26092db06453ed5a6047d8905e37fdad17eaf9c4947048498551425c605",
    LOOP_CREATOR: "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    MOD_CREATOR: "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    PATCH_CREATOR: "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    OUTER_CREATOR: "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    GAMMA_CREATOR: "38d3a8e5abd964c42b135410485f6ad58dfe9f7656035b4b5c7b1ed48c2ebb01",
}

CK = ROOT / "tmp/original_native_cross2_B_source_v1"
COMPONENT = "cross2_B_source_v1"
INPUT_SCHEMA = "native_cross2_B_source_inputs_v1"
RESULT_SCHEMA = "original_native_cross2_B_source_v1"
DIAGNOSTIC_SCHEMA = "native_cross2_B_source_stage_diagnostic_v1"
INHERITED_SOURCE_COUNT = 1876
FINAL_SOURCE_COUNT = 1879
FAILED_ARTIFACT_COUNT = 32
PARENT_AUTHENTICATED_CHECKPOINT_COUNT = 152
INHERITED_CHECKPOINT_COUNT = 156
PARENT_CURRENT_CHECKPOINT_COUNT = 4
NEW_CHECKPOINT_COUNT = 4
FINAL_CHECKPOINT_COUNT = 160
TARGET_NODES = [529, 530]
PARENT_ACCEPTED_NODES = [520, 521, 522, 523, 524, 525, 526, 527, 528, 531, 532, 533, 534, 535]
ACCEPTED_NODES = PARENT_ACCEPTED_NODES + TARGET_NODES
ORDINARY_Z_ORDERS = [0, 1]
OWNED_Z_CELL = ["0", "1/2"]

EXPECTED_RECORD_NAMES = [
    PARENT, PRIOR_J3, WEIGHTED, BACKGROUND, LOOP, MOD, PATCH, OUTER,
    C6_CREATOR, BRIDGE_CREATOR, LOOP_CREATOR, MOD_CREATOR, PATCH_CREATOR,
    OUTER_CREATOR, GAMMA_CREATOR,
]

# Keep these names aligned with the accepted outer-B source report so a later
# Publisher can merge the cross2 result without inventing a second flag ledger.
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
    return tuple(CK / ("new_cross2_B_v1_" + kind + ".json")
                 for kind in ("source", "entered", "result"))


def diagnostic_path():
    return CK / "new_cross2_B_v1_diagnostics.jsonl"


def _parent_scope(parent):
    require(parent.get("schema") == "original_native_inner_B_v1_stage_v1",
            "Inner parent schema changed")
    require(parent.get("component") == "inner_B_source_v1"
            and parent.get("source_only") is True
            and parent.get("accepted_identity") is True
            and parent.get("full_goal_complete") is False,
            "Inner parent acceptance scope changed")
    expected_code = {
        PARENT_INNER_RUNNER: PARENT_INNER_RUNNER_SHA,
        PARENT_INNER_KERNEL: PARENT_INNER_KERNEL_SHA,
    }
    require(parent.get("current_code_hashes") == expected_code
            and parent.get("implementation_hashes") == expected_code,
            "Frozen inner implementation inventory changed")
    pins = parent.get("source_pins")
    require(isinstance(pins, dict)
            and len(pins) == INHERITED_SOURCE_COUNT
            and parent.get("source_pin_count") == INHERITED_SOURCE_COUNT
            and parent.get("input_source_pin_count") == INHERITED_SOURCE_COUNT,
            "Inner parent source inventory changed")
    require(all(isinstance(name, str) and isinstance(value, str)
                and len(value) == 64 for name, value in pins.items()),
            "Malformed inner source hash")
    require(pins.get(PARENT_INNER_RUNNER) == PARENT_INNER_RUNNER_SHA
            and pins.get(PARENT_INNER_KERNEL) == PARENT_INNER_KERNEL_SHA,
            "Inner code pins are not authenticated")
    require(parent.get("accepted_native_B_nodes") == PARENT_ACCEPTED_NODES
            and parent.get("remaining_B_identity_nodes") == TARGET_NODES,
            "Inner accepted or remaining nodes changed")
    component = parent.get("component_result", {})
    require(component.get("schema") == "original_native_inner_B_source_v1"
            and component.get("component") == "inner_B_source_v1"
            and component.get("passed") is True
            and component.get("source_only") is True
            and component.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
            and component.get("ordinary_Z_not_temporal_n") is True
            and component.get("owned_Z_cell") == OWNED_Z_CELL
            and component.get("numeric_values_computed") is False
            and component.get("old_operators_replayed") is False,
            "Inner component certificate changed")
    require(parent.get("component_result_exact") is True
            and parent.get("component_result_is_projection") is False
            and parent.get("full_result_exact_saved_json") is True
            and parent.get("full_result") == parent.get("component_result"),
            "Inner full result is not the exact saved result")
    require(parent.get("raw_bytes_lossless") is True,
            "Inner lossless artifact marker changed")
    checkpoints = parent.get("checkpoints")
    require(isinstance(checkpoints, list) and len(checkpoints) == 4,
            "Inner lossless checkpoint list changed")
    for item in checkpoints:
        require(isinstance(item, dict)
                and item.get("encoding") == "base64"
                and isinstance(item.get("base64_bytes"), str)
                and isinstance(item.get("bytes"), int)
                and isinstance(item.get("sha256"), str)
                and len(item["sha256"]) == 64,
                "Lossless checkpoint metadata missing")
    result_checkpoint = checkpoints[2]
    require(parent.get("component_result_raw_sha256") == result_checkpoint.get("sha256")
            and parent.get("component_result_raw_bytes") == result_checkpoint.get("bytes"),
            "Saved full-result raw hash or byte count changed")
    execution_artifact = parent.get("execution_artifact", {})
    require(execution_artifact.get("metadata_only") is True
            and execution_artifact.get("full_result_exact") is True
            and execution_artifact.get("lossless_checkpoint_count") == 4,
            "Saved execution artifact is not lossless")
    audit = component.get("source_inspection_audit", {})
    require(audit.get("new_source_consumer_calls") == 1
            and audit.get("accepted_source_consumer_calls") == 0
            and audit.get("old_proof_functions_called") == 0
            and audit.get("runtime_model_metadata") is None,
            "Saved inner raw-result audit changed")
    require(parent.get("accepted_source_consumer_calls") == 1
            and parent.get("execution_seconds") is not None
            and parent.get("diagnostic_event_count") == 8,
            "Inner report acceptance audit changed")
    require(component.get("execution_seconds") == parent.get("execution_seconds")
            and component.get("diagnostic_events_saved_before_validation") == 8,
            "Inner component execution metadata changed")
    for audit_record in (parent.get("source_inspection_audit", {}), audit):
        require(audit_record.get("old_proof_functions_called", 0) == 0
                and audit_record.get("old_operator_calls", 0) == 0
                and audit_record.get("quadratures", 0) == 0
                and audit_record.get("recurrence_calls", 0) == 0,
                "Inner saved audit contains an old or numeric call")
    failed = parent.get("preserved_failed_artifact_hashes_union")
    require(isinstance(failed, dict) and len(failed) == FAILED_ARTIFACT_COUNT
            and parent.get("preserved_failed_artifact_count_union") == FAILED_ARTIFACT_COUNT,
            "Inner failed-artifact union changed")
    require(parent.get("parent_commit") == PARENT_ANCESTOR_COMMIT,
            "Inner parent ancestry changed")
    return pins, failed


def _parent_checkpoint_inventory(parent):
    """Authenticate the parent canonical 152 plus current 4 checkpoints."""
    inventory = parent.get("old_checkpoint_inventory")
    require(isinstance(inventory, dict), "Inner checkpoint inventory missing")
    inherited = inventory.get("authenticated_checkpoint_hashes")
    inherited_old = inventory.get("inherited_parent_old_checkpoint_hashes")
    parent_current = inventory.get("parent_current_checkpoint_hashes")
    current = parent.get("checkpoint_hashes")
    require(isinstance(inherited, dict) and len(inherited) == PARENT_AUTHENTICATED_CHECKPOINT_COUNT,
            "Authenticated parent checkpoint count changed")
    require(isinstance(inherited_old, dict) and len(inherited_old) == 148
            and isinstance(parent_current, dict) and len(parent_current) == 4
            and set(inherited_old) | set(parent_current) == set(inherited),
            "Canonical parent checkpoint topology changed")
    require(isinstance(current, dict) and len(current) == PARENT_CURRENT_CHECKPOINT_COUNT
            and set(inherited).isdisjoint(current)
            and len(set(inherited) | set(current)) == INHERITED_CHECKPOINT_COUNT,
            "Current parent checkpoint set changed")
    require(inventory.get("authenticated_parent_old_count") == 148
            and inventory.get("authenticated_total_count") == PARENT_AUTHENTICATED_CHECKPOINT_COUNT
            and inventory.get("parent_current_count") == 4
            and inventory.get("distinct_parent_count") == PARENT_AUTHENTICATED_CHECKPOINT_COUNT
            and inventory.get("raw_sha_only") is True
            and inventory.get("history_decode_forbidden") is True,
            "Canonical parent checkpoint counters changed")
    historical = inventory.get("historical_parent_inventory")
    require(isinstance(historical, dict), "Historical parent inventory missing")
    combined = {**inherited, **current}
    require(all(isinstance(name, str) and isinstance(value, str)
                and len(value) == 64 for name, value in combined.items()),
            "Malformed parent checkpoint hash")
    for name, expected in combined.items():
        _read_raw_hash(name, expected)
    return {
        "authenticated_checkpoint_hashes": dict(inherited),
        "parent_current_checkpoint_hashes": dict(current),
        "authenticated_parent_old_count": PARENT_AUTHENTICATED_CHECKPOINT_COUNT,
        "authenticated_total_count": PARENT_AUTHENTICATED_CHECKPOINT_COUNT,
        "parent_current_count": PARENT_CURRENT_CHECKPOINT_COUNT,
        "distinct_parent_count": INHERITED_CHECKPOINT_COUNT,
        "historical_parent_inventory": json.loads(json.dumps(
            historical, sort_keys=True, ensure_ascii=True, allow_nan=False)),
        "raw_sha_only": True,
        "history_decode_forbidden": True,
    }


def _static_scope(config):
    require(config["SCHEMA"] == RESULT_SCHEMA, "Kernel schema changed")
    names = list(config["RECORD_NAMES"])
    hashes = dict(config["RECORD_HASHES"])
    decoded = dict(config["DECODED_RECORD_HASHES"])
    require(names == list(dict.fromkeys(names))
            and set(names) == set(EXPECTED_RECORD_NAMES),
            "Kernel record inventory changed")
    require(set(hashes) == set(names) and set(decoded) == set(names),
            "Kernel record maps do not match")
    require(all(isinstance(value, str) and len(value) == 64 for value in hashes.values()),
            "Kernel raw hashes are not SHA-256 literals")
    require(all(value is None or (isinstance(value, str) and len(value) == 64)
                for value in decoded.values()), "Kernel decoded hashes are malformed")
    require(config["TARGET_NODES"] == TARGET_NODES,
            "Kernel target nodes changed")
    stage = list(config["STAGE_FLAGS"])
    false = list(config["FALSE_FLAGS"])
    require(stage and len(stage) == len(set(stage))
            and len(false) == len(set(false))
            and set(stage).isdisjoint(false), "Kernel flag ledger changed")
    aliases = dict(config["SOURCE_TEXT_ALIASES"])
    require(set(aliases) == {"c6", "bridge", "loop", "mod", "patch", "outer", "gamma"}
            and set(aliases.values()) <= set(names),
            "Kernel creator aliases changed")
    require(aliases["c6"] == C6_CREATOR and aliases["bridge"] == BRIDGE_CREATOR
            and aliases["loop"] == LOOP_CREATOR and aliases["mod"] == MOD_CREATOR
            and aliases["patch"] == PATCH_CREATOR and aliases["outer"] == OUTER_CREATOR
            and aliases["gamma"] == GAMMA_CREATOR,
            "Kernel creator aliases are not the authenticated records")
    return names, hashes, decoded, aliases, stage, false


def _prepare_packet():
    config = literal_configuration()
    names, expected_hashes, decoded_expected, aliases, stage, false = _static_scope(config)
    # Authenticate the parent report as raw bytes before using its metadata to
    # enumerate preserved artifacts and checkpoints.  No selected source
    # record is decoded until all inherited, failed, and checkpoint bytes have
    # completed their raw hash pass below.
    parent_raw = source_path(PARENT).read_bytes()
    require(hashlib.sha256(parent_raw).hexdigest() == PARENT_SHA,
            "Changed inner parent report")
    parent = json.loads(parent_raw)
    inherited, failed = _parent_scope(parent)
    require(PARENT not in inherited and RUNNER not in inherited and KERNEL not in inherited,
            "Parent source union already contains cross2 publication")

    # Raw inventory authentication precedes every optional JSON or gzip decode.
    for name, expected in inherited.items():
        _read_raw_hash(name, expected)
    for name, expected in failed.items():
        _read_raw_hash(name, expected)
    checkpoint_inventory = _parent_checkpoint_inventory(parent)

    current_code_hashes = {RUNNER: sha(ROOT / RUNNER), KERNEL: sha(ROOT / KERNEL)}
    final_union = dict(inherited)
    final_union[PARENT] = PARENT_SHA
    final_union.update(current_code_hashes)
    require(len(final_union) == FINAL_SOURCE_COUNT,
            "Final source union must be parent report plus code2")
    require(names == EXPECTED_RECORD_NAMES,
            "Cross2 kernel record ordering changed")
    require(set(final_union) >= set(names),
            "Kernel record escapes the authenticated source union")
    for name, expected in expected_hashes.items():
        require(expected == final_union[name], "Kernel record hash changed: " + name)

    records = {}
    for name in names:
        records[name] = _read_record(
            name, expected_hashes[name], decoded_expected[name], decode=not _raw_only(name))
    input_digest = digest(records)
    record_decoded = {name: item["decoded_sha256"] for name, item in records.items()}
    raw_only = [name for name, item in records.items() if item.get("raw_only")]
    static_receipt = {
        "ast_configuration_read": True,
        "AST_receipt_source_calls": 0,
        "source_helper_hash": None,
        "kernel_imported": False,
        "derive_calls": 0,
        "source_consumer_calls": 0,
        "old_producer_calls": 0,
        "old_operator_calls": 0,
        "old_proof_calls": 0,
        "provider_calls": 0,
        "numeric_values_computed": False,
        "runtime_model_metadata": None,
        "actions_executed": False,
    }
    return dict(
        schema=INPUT_SCHEMA,
        result_schema=RESULT_SCHEMA,
        component=COMPONENT,
        parent_report=PARENT,
        parent_report_sha256=PARENT_SHA,
        parent_commit=PARENT_COMMIT,
        parent_internal_commit=PARENT_INTERNAL_COMMIT,
        parent_ancestor_commit=PARENT_ANCESTOR_COMMIT,
        prior_J3_report=PRIOR_J3,
        prior_J3_report_sha256=PRIOR_J3_SHA,
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
        static_receipt=static_receipt,
        diagnostics_policy="append-only JSONL; snapshot, flush, and fsync before callback return",
        diagnostics_checkpoint_additional_to_source_entered_result=True,
        source_functions_executed=0,
        runtime_model_metadata=None,
        raw_authentication_policy={
            "all_inherited_source_pins_before_decode": True,
            "all_failed_artifacts_before_decode": True,
            "all_parent_checkpoints_before_decode": True,
            "forbidden_history_raw_only": True,
            "all_gzip_raw_sha_before_decode": True,
        },
    )


def prepare_cross():
    source, entered, result = paths()
    journal = diagnostic_path()
    require(not any(path.exists() for path in (source, entered, result, journal)),
            "Cross2 checkpoints already exist; never rerun")
    packet = _prepare_packet()
    _write_new(source, packet)
    return {
        "passed": True,
        "prepared_only": True,
        "source_sha256": sha(source),
        "code_sha256": dict(packet["current_code_hashes"]),
        "input_digest": packet["input_digest"],
        "source_pin_count": packet["source_pin_manifest"]["final_union_count"],
        "record_count": len(packet["record_names"]),
        "static_receipt": packet["static_receipt"],
    }


def _load_packet(source):
    packet = json.loads(source.read_bytes())
    require(packet.get("schema") == INPUT_SCHEMA
            and packet.get("result_schema") == RESULT_SCHEMA
            and packet.get("component") == COMPONENT,
            "Frozen source packet schema changed")
    manifest = packet.get("source_pin_manifest", {})
    require(manifest.get("final_union_count") == FINAL_SOURCE_COUNT
            and packet.get("record_count") == len(EXPECTED_RECORD_NAMES)
            and packet.get("record_names") == EXPECTED_RECORD_NAMES,
            "Frozen source union or record count changed")
    require(packet.get("new_checkpoint_count") == NEW_CHECKPOINT_COUNT
            and packet.get("final_checkpoint_count") == FINAL_CHECKPOINT_COUNT,
            "Frozen checkpoint count changed")
    require(packet.get("parent_report_sha256") == PARENT_SHA
            and packet.get("prior_J3_report_sha256") == PRIOR_J3_SHA,
            "Frozen parent hashes changed")
    require(set(packet.get("current_code_hashes", {})) == {RUNNER, KERNEL}
            and packet.get("implementation_hashes") == packet["current_code_hashes"],
            "Frozen implementation inventory changed")
    inventory = packet.get("parent_checkpoint_inventory", {})
    require(isinstance(inventory, dict)
            and len(inventory.get("authenticated_checkpoint_hashes", {})) == PARENT_AUTHENTICATED_CHECKPOINT_COUNT
            and len(inventory.get("parent_current_checkpoint_hashes", {})) == 4
            and inventory.get("distinct_parent_count") == INHERITED_CHECKPOINT_COUNT,
            "Frozen parent checkpoint inventory changed")
    for name, expected in packet["current_code_hashes"].items():
        require(sha(ROOT / name) == expected,
                "Current implementation changed after prepare: " + name)
    return packet


def _load_frozen_records(packet):
    for name, expected in packet["inherited_source_hashes"].items():
        _read_raw_hash(name, expected)
    for name, expected in packet["preserved_failed_artifact_hashes"].items():
        _read_raw_hash(name, expected)
    inventory = packet["parent_checkpoint_inventory"]
    for group in ("authenticated_checkpoint_hashes", "parent_current_checkpoint_hashes"):
        for name, expected in inventory[group].items():
            _read_raw_hash(name, expected)
    raw_only = set(packet.get("raw_only_records", []))
    require(raw_only == {name for name in packet["record_names"] if _raw_only(name)},
            "Frozen raw-only source policy changed")
    records = {}
    for name in packet["record_names"]:
        records[name] = _read_record(
            name, packet["record_hashes"][name],
            packet["record_decoded_sha256"][name], decode=name not in raw_only)
        require(records[name]["raw_sha256"] == packet["record_raw_sha256"][name],
                "Frozen raw record changed: " + name)
        require(records[name]["decoded_sha256"] == packet["record_decoded_sha256"][name],
                "Frozen decoded record changed: " + name)
    require(digest(records) == packet["input_digest"], "Frozen input digest changed")
    return records


def _validate_result(result, packet):
    require(isinstance(result, dict)
            and result.get("schema") == RESULT_SCHEMA
            and result.get("component") == COMPONENT
            and result.get("source_only") is True
            and type(result.get("passed")) is bool,
            "Kernel result schema or component changed")
    require(result.get("input_digest") == packet["input_digest"],
            "Kernel result input digest changed")
    require(result.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
            and result.get("ordinary_Z_not_temporal_n") is True
            and result.get("owned_Z_cell") == OWNED_Z_CELL,
            "Kernel ordinary-Z scope changed")
    require(result.get("target_B_identity_nodes") == TARGET_NODES,
            "Kernel target scope changed")
    for flag in packet["false_flags"]:
        require(result.get(flag) is False, "Unadmitted scope changed: " + flag)
    for optional in ("Q_y_identity_claim", "Q_y_identity_proved"):
        if optional in result:
            require(result[optional] is False, "Q_y equality was claimed")
    stage_values = [result.get(flag) for flag in packet["stage_flags"]]
    require(all(type(value) is bool for value in stage_values),
            "Kernel stage flags are not booleans")
    declared = result.get("stage_flags")
    if declared is not None:
        require(isinstance(declared, dict)
                and all(declared.get(flag) is result.get(flag)
                        for flag in packet["stage_flags"]),
                "Kernel declared stage flags changed")
    all_stages = all(stage_values)
    require((not all_stages) or result["passed"],
            "All declared kernel stages cannot be true on a failed result")
    gates = result.get("unresolved_gates")
    require(isinstance(gates, list), "Kernel unresolved gates missing")
    if result["passed"]:
        require(all_stages and not gates, "Passed result lacks declared stage completion")
        require(result.get("newly_proved_B_identity_nodes") == TARGET_NODES
                and result.get("accepted_native_B_nodes") == ACCEPTED_NODES
                and result.get("remaining_B_identity_nodes") == []
                and result.get("remaining_B_identity_count") == 0
                and result.get("all_native_B_jumps_zero_proved") is True,
                "Cross2 acceptance set changed")
        scope = result.get("B_identity_scope", {})
        per_node = scope.get("per_node") if isinstance(scope, dict) else None
        require(isinstance(per_node, dict)
                and set(per_node) == {str(node) for node in TARGET_NODES},
                "Cross2 per-node scope changed")
        for node in TARGET_NODES:
            evidence = per_node[str(node)]
            require(isinstance(evidence, dict)
                    and all(evidence.get(key) is True for key in (
                        "passed", "R_C1_proved", "U_C1_proved", "m_Z2_proved",
                        "Q_C1_proved", "B_C1_proved",
                        "native_weighted_jump_zero_proved"))
                    and evidence.get("source_only") is True
                    and evidence.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS,
                    "Cross2 per-node proof flags changed")
            if "Q_y_identity_claim" in evidence:
                require(evidence["Q_y_identity_claim"] is False,
                        "Cross2 per-node Q_y claim changed")
    else:
        if "all_native_B_jumps_zero_proved" in result:
            require(result["all_native_B_jumps_zero_proved"] is False,
                    "Failed result admitted all native B jumps")
        require(not all_stages or bool(gates), "Failed result lacks a failed stage or gate")
    for key in ("old_proof_functions_called", "old_operator_calls",
                "old_operator_function_calls", "source_producer_function_calls",
                "old_producer_function_calls", "new_cross_calls",
                "numerical_evaluations", "quadratures", "recurrence_calls"):
        if key in result:
            require(result[key] == 0 or (key == "new_cross_calls" and result[key] == 1),
                    "Unexpected call count: " + key)
    require(result.get("old_operators_replayed", False) is False
            and result.get("numeric_values_computed") is False,
            "Unadmitted replay or numeric evaluation")
    require(result.get("native_IN_reduction_is_source_only") is True
            and result.get("archived_all_16_J_terms_retained") is True
            and result.get("source_jump_zero_does_not_delete_quantitative_errors") is True,
            "Inherited construction ledger was not preserved")
    retained = result.get("retained_construction_data")
    require(isinstance(retained, dict)
            and retained.get("cumulative_terms") == 212
            and retained.get("raw_signed_terms") == 213
            and retained.get("reduced_terms") == 183
            and retained.get("six_signed_channels_retained") is True
            and retained.get("all_old_error_data_retained") is True
            and retained.get("reserved_contribution_registered_once") is True
            and retained.get("inherited_errors_deleted") is False,
            "Inherited term ledger changed")
    integrals = result.get("retained_integral_term_inventory")
    require(isinstance(integrals, dict)
            and integrals.get("K_terms") == 17
            and integrals.get("outer_B_terms") == 2
            and integrals.get("K_terms_retained") is True
            and integrals.get("outer_B_terms_retained") is True,
            "Integral-term inventory changed")
    frozen_errors = result.get("old_errors_frozen")
    require(isinstance(frozen_errors, dict)
            and frozen_errors.get("cumulative_terms") == 212
            and frozen_errors.get("raw_signed_terms") == 213
            and frozen_errors.get("reduced_terms") == 183
            and frozen_errors.get("all_old_packet_errors_retained") is True,
            "Frozen error ledger changed")
    audit = result.get("source_inspection_audit")
    if isinstance(audit, dict):
        require(audit.get("new_source_consumer_calls") == 1
                and audit.get("accepted_source_consumer_calls") == 0
                and audit.get("accepted_ancestor_operator_calls") == 0
                and audit.get("old_proof_functions_called", 0) == 0
                and audit.get("runtime_model_metadata") is None,
                "Cross2 source consumer audit changed")
    require(result.get("runtime_model_metadata") is None,
            "Runtime model metadata must remain null")
    remaining = result.get("remaining_construction_gates")
    require(isinstance(remaining, list) and bool(remaining),
            "Open construction gates were dropped")
    if "binding_fix_count" in result:
        require(result["binding_fix_count"] == 3, "Binding-fix count changed")


def _empty_journal(path):
    with path.open("xb") as stream:
        stream.flush()
        os.fsync(stream.fileno())


def cross():
    source, entered, target = paths()
    journal = diagnostic_path()
    require(source.exists() and not entered.exists() and not target.exists()
            and not journal.exists(), "Frozen cross2 packet is absent or already entered")
    packet = _load_packet(source)
    _write_new(entered, {
        "schema": "native_cross2_B_source_entered_v1",
        "component": COMPONENT,
        "source_sha256": sha(source),
        "input_digest": packet["input_digest"],
        "implementation_hashes": packet["current_code_hashes"],
        "kernel_sha256": packet["current_code_hashes"][KERNEL],
        "entered_at_utc": datetime.now(timezone.utc).isoformat(),
        "kernel_imported_before_entry": False,
        "derive_calls_before_entry": 0,
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
            "new_native_cross2_B_source_v1", ROOT / KERNEL)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        result = module.derive(
            records, expected_input_digest=packet["input_digest"],
            diagnostic_callback=persist_diagnostic)
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
            "all_native_B_jumps_zero_proved": False,
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
            "all_native_B_jumps_zero_proved": False,
        }
    result.update(
        execution_seconds=time.perf_counter() - started,
        source_checkpoint_sha256=sha(source),
        entry_checkpoint_sha256=sha(entered),
        runner_caught_exception=caught,
        runner_full_goal_complete=False,
        runtime_model_metadata=None,
        stage_diagnostics=list(events),
        diagnostic_events=list(events),
        diagnostic_checkpoint_path=journal.relative_to(ROOT).as_posix(),
        diagnostic_checkpoint_sha256=sha(journal),
        diagnostic_events_saved_before_validation=len(events),
        full_result_saved_before_validation=True,
    )
    # The result is durable before any validation can reject it; it is never
    # overwritten by this runner.
    _write_new(target, result)
    if records is not None:
        require(digest(records) == packet["input_digest"],
                "Input mutated during cross2 source consumer")
    if not caught:
        _validate_result(result, packet)
    return {
        "passed": result["passed"],
        "result_sha256": sha(target),
        "code_sha256": dict(packet["current_code_hashes"]),
        "seconds": result["execution_seconds"],
        "error": result.get("error"),
        "diagnostic_events": len(events),
    }


def inspect_cross():
    """Inspect saved metadata only; never decode inputs or import the kernel."""
    source, entered, target = paths()
    journal = diagnostic_path()
    require(all(path.exists() for path in (source, entered, target, journal)),
            "Saved cross2 checkpoints are incomplete")
    packet = json.loads(source.read_bytes())
    entry = json.loads(entered.read_bytes())
    result = json.loads(target.read_bytes())
    require(packet.get("schema") == INPUT_SCHEMA
            and packet.get("result_schema") == RESULT_SCHEMA
            and packet.get("component") == COMPONENT,
            "Saved packet metadata changed")
    require(packet.get("source_pin_manifest", {}).get("final_union_count") == FINAL_SOURCE_COUNT
            and packet.get("record_names") == EXPECTED_RECORD_NAMES
            and packet.get("record_count") == len(EXPECTED_RECORD_NAMES),
            "Saved packet inventory changed")
    for name, expected in packet.get("current_code_hashes", {}).items():
        require(sha(ROOT / name) == expected, "Saved implementation changed: " + name)
    require(entry.get("schema") == "native_cross2_B_source_entered_v1"
            and entry.get("source_sha256") == sha(source)
            and entry.get("input_digest") == packet.get("input_digest")
            and entry.get("implementation_hashes") == packet.get("current_code_hashes")
            and entry.get("kernel_imported_before_entry") is False
            and entry.get("derive_calls_before_entry") == 0,
            "Saved entry metadata changed")
    require(result.get("source_checkpoint_sha256") == sha(source)
            and result.get("entry_checkpoint_sha256") == sha(entered)
            and result.get("input_digest") == packet.get("input_digest")
            and result.get("schema") == RESULT_SCHEMA
            and result.get("component") == COMPONENT
            and result.get("source_only") is True
            and type(result.get("passed")) is bool
            and result.get("diagnostic_checkpoint_path") ==
            journal.relative_to(ROOT).as_posix()
            and result.get("diagnostic_checkpoint_sha256") == sha(journal)
            and result.get("full_result_saved_before_validation") is True,
            "Saved result checkpoint links changed")
    lines = journal.read_text(encoding="utf-8").splitlines()
    events = [json.loads(line) for line in lines]
    require(all(event.get("schema") == DIAGNOSTIC_SCHEMA
                and event.get("input_digest") == packet.get("input_digest")
                and event.get("index") == index
                for index, event in enumerate(events)),
            "Saved diagnostic journal metadata changed")
    require(result.get("stage_diagnostics") == events
            and result.get("diagnostic_events") == events
            and result.get("diagnostic_events_saved_before_validation") == len(events),
            "Saved diagnostic/result snapshots changed")
    if result.get("runner_caught_exception"):
        require(result.get("passed") is False
                and result.get("failed_attempt_preserved") is True
                and isinstance(result.get("traceback"), str),
                "Saved exception envelope changed")
    return {
        "passed": True,
        "saved_only": True,
        "component_passed": result.get("passed"),
        "full_goal_complete": False,
        "code_sha256": packet.get("current_code_hashes"),
        "static_receipt": packet.get("static_receipt"),
        "diagnostic_events": len(events),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare_cross", "cross", "inspect_cross"))
    action = parser.parse_args().action
    output = prepare_cross() if action == "prepare_cross" else (
        cross() if action == "cross" else inspect_cross())
    print(json.dumps(output, ensure_ascii=True, allow_nan=False))
