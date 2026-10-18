"""Checkpoint a paired constructor proof without polynomial expansion.

The older canonical expansion is an unfinished, unaccepted attempt. This
runner neither restarts it nor invokes an accepted mathematical operator.
"""
from __future__ import annotations

import argparse
import ast
import base64
import copy
from datetime import datetime, timezone
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
PARENT = PREFIX + "Rm_patch_native_provenance.json"
PARENT_HASH = "e649b4f1556c3c04eadbf9469d79e14795d6be2da8a8db52669a9e5fdebc05ea"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
FINITE_HASH = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
PORTS = "experiments/root_st073/original_Rm_loop_ports_admission.json"
PORTS_HASH = "0331c5fc37bd3ec6c0d386c38b07d9863ede622e14f1c6dc65bf8cdf675a16b9"
GRAPH_HASH = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
RUNNER = "scripts/run_original_Rm_slow_native_structural_values.py"
KERNEL = "src/openai_ns_reconstruction/original_Rm_slow_native_structural_values.py"
FROZEN_LIBRARY = "src/openai_ns_reconstruction/original_Rm_slow_native_value_seam.py"
OLD_RUNNER = "scripts/run_original_Rm_slow_native_value_seam.py"
OLD_CHECKPOINTS = ROOT / "tmp/original_Rm_slow_native_value_seam"
OLD_SOURCE = OLD_CHECKPOINTS / "new_source.json"
OLD_ENTRY = OLD_CHECKPOINTS / "new_operator_entered.json"
OLD_RESULT = OLD_CHECKPOINTS / "new_slow_native_value_seam.json"
OLD_SOURCE_HASH = "cb5f3b5c61d7c5175a11847ca66c1f20c57b4623dc91aac45ee429459e1c409a"
OLD_ENTRY_HASH = "32ac21b523efab8d58f25d18c9a2de363e1babb776cb9883a8955ba5837dd14d"
BUNDLE = ROOT / "experiments/root_st073/original_Rm_slow_canonical_expansion_attempt.json"
BUNDLE_HASH = "53327015c5779c388c6becc18d6727ac90441d3b4bca44840dbda8e502d3ed66"
DISPOSITION = ROOT / "experiments/root_st073/original_Rm_slow_canonical_expansion_disposition.json"
CHECKPOINTS = ROOT / "tmp/original_Rm_slow_native_structural_values"
SOURCE = CHECKPOINTS / "new_source.json"
ENTRY = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_structural_values.json"
REPORT = BASE / (PREFIX + "Rm_slow_native_structural_values.json")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_for(name):
    path = ROOT / name if name.startswith(("src/", "scripts/", "experiments/", "docs/")) else BASE / name
    path = path.resolve()
    require(path.is_relative_to(ROOT.resolve()), "source pin escapes repository: " + name)
    return path


def read_saved(path):
    return json.loads(path.read_bytes())


def read_pinned(path, expected):
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, "source pin changed: " + str(path))
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def write_new(path, value):
    require(not path.exists(), "artifact already exists: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True)
        stream.write("\n")


def preflight():
    for name in (RUNNER, KERNEL):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + name)
        ast.parse(raw, filename=name)


def append_saved_delta(nodes, packet, base_key, full_key):
    require(digest(nodes) == packet[base_key], "saved delta base changed")
    require(len(nodes) == packet.get("base_node_count", packet.get("base_graph_node_count")), "saved delta count changed")
    result = nodes + copy.deepcopy(packet["appended_nodes"])
    require(len(result) == packet["full_node_count"] and digest(result) == packet[full_key], "saved delta changed")
    return result


def load_saved_graph(parent):
    archive = read_pinned(BASE / FINITE, FINITE_HASH)
    require(archive["source_family"] == parent["source_family"], "source family changed")
    nodes = archive["exact_graph_nodes"]
    require(len(nodes) == 205147, "actual finite source count changed")
    for packet, base_key, full_key in (
        (parent["actual_Rm_native_radial_endpoint_source"], "base_digest", "full_nodes_digest"),
        (parent["actual_Rm_second_radial_endpoint_functions"], "base_digest", "full_nodes_digest"),
        (parent["actual_Rm_history_native_source"]["graph_delta"], "base_nodes_sha256", "full_canonical_nodes_sha256"),
        (parent["actual_Rm_history_seam_identity"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256"),
        (parent["actual_Rm_modified_moment_germ"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256"),
    ):
        nodes = append_saved_delta(nodes, packet, base_key, full_key)
    require(len(nodes) == 205300 and digest(nodes) == GRAPH_HASH, "current saved graph changed")
    return nodes


def paired_plan(left, right, common):
    require(common["passed"] is True and common["source_only"] is True
        and common["actual_common_phase_proved"] is True, "common phase is not admitted")
    require(common["coordinates"] == [84528, 185676] and common["endpoints"] == [1, 1], "common phase endpoints changed")
    require(common["phase_nodes"] == [119416, 196046], "common phase handles changed")
    require(common["radius_lemma_digest"] == "a719562a3503be0d7b4672633adb88f7f16ed587638d58c0d5e027d9bc406b24", "radius lemma changed")
    reduction = common["endpoint_radius_reduction"]
    require(reduction["left"]["closed_node"] == 85743 and reduction["right"]["closed_node"] == 187963
        and reduction["right"]["body_node"] == 187214 and reduction["common_radius_source_node"] == 85743, "actual endpoint radius ports changed")
    pairs = []
    for field in ("A", "B"):
        for order in (0, 1):
            row = dict(field=field, ordinary_order=order)
            for side, provenance in (("postrestore", left), ("patch", right)):
                matches = [r for r in provenance["rows"] if r["field"] == field and r["ordinary_order"] == order]
                require(len(matches) == 1, "actual target row missing or duplicated")
                selected = matches[0]
                row[side] = {key:selected[key] for key in ("source_handle", "derivative_handle", "physical_target_handle")}
            pairs.append(row)
    return dict(schema="actual_Rm_paired_native_value_targets_v1",
        native_coordinates=dict(zip(("postrestore", "patch"), common["coordinates"])),
        native_endpoints=dict(zip(("postrestore", "patch"), common["endpoints"])),
        held_phase=dict(zip(("postrestore", "patch"), common["phase_nodes"])), targets=pairs)


def setup():
    preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "inherited source pin changed: " + name)
    bundle = read_pinned(BUNDLE, BUNDLE_HASH)
    original_bytes = {}
    for item in bundle["files"]:
        raw = base64.b64decode(item["base64_bytes"], validate=True)
        require(len(raw) == item["bytes"] and hashlib.sha256(raw).hexdigest() == item["sha256"], "old attempt snapshot changed")
        original_bytes[item["repository_path"]] = raw
        actual_path = ROOT / item["repository_path"]
        if actual_path.exists():
            require(sha(actual_path) == item["sha256"], "original attempt source bytes changed")
    source_bytes = original_bytes[OLD_SOURCE.relative_to(ROOT).as_posix()]
    entry_bytes = original_bytes[OLD_ENTRY.relative_to(ROOT).as_posix()]
    require(hashlib.sha256(source_bytes).hexdigest() == OLD_SOURCE_HASH
        and hashlib.sha256(entry_bytes).hexdigest() == OLD_ENTRY_HASH, "original input/entry archive changed")
    old_source = json.loads(source_bytes)
    require(old_source["parent_sha256"] == PARENT_HASH and old_source["nodes_digest"] == GRAPH_HASH, "original snapshot belongs to another source")
    require(hashlib.sha256(original_bytes[FROZEN_LIBRARY]).hexdigest() == old_source["implementation_hashes"][FROZEN_LIBRARY], "frozen metadata library snapshot changed")
    nodes = load_saved_graph(parent)
    ports = read_pinned(ROOT / PORTS, PORTS_HASH)
    left = parent["actual_Rm_slow_native_provenance"]["sides"]["postrestore"]
    right = parent["actual_Rm_patch_native_provenance"]
    require(left["passed"] is True and right["passed"] is True, "accepted actual derivative provenance is missing")
    plan = paired_plan(left, right, ports["common_phase_source"])
    inputs = dict(generic_packet=parent["actual_Rm_generic_native_identity"], loop_ports_admission=ports,
        left_provenance=left, right_provenance=right, paired_target_plan=plan, expected_nodes_digest=GRAPH_HASH)
    packet = dict(schema="actual_Rm_paired_constructor_VALUE_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, graph_node_count=len(nodes), graph_digest=GRAPH_HASH,
        paired_target_plan=plan, paired_target_plan_digest=digest(plan), input_digest=digest(inputs),
        input_data_digests={key:digest(inputs[key]) for key in ("generic_packet", "loop_ports_admission", "left_provenance", "right_provenance", "paired_target_plan")},
        implementation_hashes={name:sha(ROOT / name) for name in (RUNNER, KERNEL)},
        bounded_loop_ports=PORTS, bounded_loop_ports_sha256=PORTS_HASH,
        original_attempt_bundle=BUNDLE.relative_to(ROOT).as_posix(), original_attempt_bundle_sha256=sha(BUNDLE),
        proof_method="exact_constructor_congruence", polynomial_expansion_or_GCD=False,
        accepted_math_operators_replayed=False, derivatives_generated=False,
        numeric_field_evaluated=False, full_goal_complete=False)
    return parent, nodes, inputs, packet


def prepare_source():
    *_, packet = setup()
    require(not any(path.exists() for path in (SOURCE, ENTRY, RESULT)), "new paired source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, sha256=sha(SOURCE), nodes=packet["graph_node_count"],
        paired_target_plan_digest=packet["paired_target_plan_digest"])


def validate_kernel_result(result, packet):
    require(isinstance(result, dict) and type(result.get("passed")) is bool, "malformed kernel result")
    require(result.get("schema") == "original_Rm_slow_native_value_seam_v1"
        and result.get("proof_method") == "exact_constructor_congruence", "kernel proof contract changed")
    require(result.get("source_only") is True and result.get("numeric_values_evaluated") is False
        and result.get("new_native_derivative_generated") is False
        and result.get("no_derivatives_generated") is True and result.get("old_operators_replayed") is False,
        "kernel execution flags violate the source-only contract")
    require(result.get("static_checks", {}).get("derivatives_generated") is False,
        "kernel derivative-generation flag missing or changed")
    for key in ("accepted_mathematical_operator_replayed", "numeric_field_evaluated", "full_goal_complete", "derivatives_generated"):
        require(key not in result or result[key] is False, "kernel returned forbidden action flag: " + key)
    require(result.get("nodes_digest") == result.get("expected_nodes_digest") == GRAPH_HASH,
        "kernel graph binding changed")
    require(result.get("input_data_digests") == packet["input_data_digests"], "kernel individual input digests changed")
    result_plan = result.get("paired_target_plan", {})
    require(all(result_plan.get(key) == value for key, value in packet["paired_target_plan"].items())
        and result_plan.get("digest") == packet["paired_target_plan_digest"], "kernel target plan binding changed")
    rows = result.get("canonical_value_pairs", [])
    require(len(rows) == 4, "kernel target coverage is incomplete")
    expected = {(row["field"],row["ordinary_order"]):row for row in packet["paired_target_plan"]["targets"]}
    seen = set()
    for row in rows:
        key = (row["field"],row["ordinary_order"])
        require(key in expected and key not in seen, "kernel target row changed or duplicated")
        seen.add(key)
        require(all(row.get(side+"_target_handle") == expected[key][side]["physical_target_handle"]
            for side in ("postrestore", "patch")), "kernel actual target handles changed")
        require((row.get("status"),row.get("equal")) in (("proved",True),("unresolved",False)), "invalid row admission status")
    admitted = not result.get("unresolved_gates") and all(row["status"] == "proved" for row in rows)
    require(result["passed"] == admitted, "incomplete structural proof marked passed")


def validate_saved_result(result, packet):
    if result.get("failed_attempt_preserved") is True:
        require(result.get("passed") is False and isinstance(result.get("error"), str), "malformed failed attempt")
    else:
        validate_kernel_result(result, packet)
    entry = read_saved(ENTRY)
    require(read_saved(SOURCE) == packet, "saved paired source changed")
    require(entry["source_sha256"] == result["source_checkpoint_sha256"] == sha(SOURCE), "new source link changed")
    require(entry["kernel_sha256"] == packet["implementation_hashes"][KERNEL], "new kernel link changed")
    require(entry["input_digest"] == packet["input_digest"]
        and entry["paired_target_plan_digest"] == packet["paired_target_plan_digest"], "new paired input changed")
    require(result["entry_checkpoint_sha256"] == sha(ENTRY), "new entry link changed")


def structural_values():
    _, nodes, inputs, packet = setup()
    require(read_saved(SOURCE) == packet, "new paired source changed")
    require(not ENTRY.exists() and not RESULT.exists(), "new paired operator already entered")
    if OLD_RESULT.exists():
        require(read_saved(OLD_RESULT)["passed"] is not True,
            "original value operator returned a successful result; inspect it instead of repeating a completed proof")
    sys.path.insert(0, str(ROOT / "src"))
    spec = importlib.util.spec_from_file_location("new_Rm_slow_native_structural_values", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTRY, dict(source_sha256=sha(SOURCE), kernel_sha256=packet["implementation_hashes"][KERNEL],
        input_digest=packet["input_digest"], paired_target_plan_digest=packet["paired_target_plan_digest"],
        entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started = time.perf_counter()
    result = None
    try:
        result = module.derive(nodes, **inputs)
        require(digest(nodes) == GRAPH_HASH and digest(inputs) == packet["input_digest"], "new operator mutated admitted input data")
        validate_kernel_result(result, packet)
    except Exception as error:
        result = dict(passed=False, error_type=type(error).__name__, error=str(error),
            failed_attempt_preserved=True, rejected_return=result)
    result.update(execution_seconds=time.perf_counter() - started, source_checkpoint_sha256=sha(SOURCE),
        entry_checkpoint_sha256=sha(ENTRY), actual_source_graph_inspected=True,
        runner_full_goal_complete=False)
    write_new(RESULT, result)
    return dict(passed=result["passed"], sha256=sha(RESULT), seconds=result["execution_seconds"],
        error=result.get("error"), unresolved_gate_count=len(result.get("unresolved_gates", [])),
        unresolved_gates=[{key:gate[key] for key in ("gate", "field", "ordinary_order", "side", "error") if key in gate}
            for gate in result.get("unresolved_gates", [])])


def inspect_result():
    _, _, _, packet = setup()
    result = read_saved(RESULT)
    validate_saved_result(result, packet)
    return dict(passed=True, serialized_only=True, mathematical_slow_seam_passed=result["passed"],
        result_sha256=sha(RESULT), source_sha256=sha(SOURCE), input_digest=packet["input_digest"],
        accepted_or_new_operator_replayed=False, full_goal_complete=False)


def build_report(parent, packet, result):
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_slow_native_paired_constructor_values", parent_report=PARENT,
        parent_report_sha256=PARENT_HASH, actual_Rm_slow_native_value_seam_inputs=packet,
        actual_Rm_slow_native_value_seam=result,
        slow_native_value_seam_checkpoint_hashes={path.name:sha(path) for path in (SOURCE, ENTRY, RESULT)},
        original_unaccepted_canonical_expansion=dict(attempt_bundle=BUNDLE.relative_to(ROOT).as_posix(),
            attempt_bundle_sha256=sha(BUNDLE), mathematical_result_transferred=False,
            disposition=read_saved(DISPOSITION) if DISPOSITION.exists() else
                dict(state="separately managed unfinished attempt; no terminal state asserted")))
    report["input_hashes"].update({PARENT:PARENT_HASH, PORTS:PORTS_HASH,
        BUNDLE.relative_to(ROOT).as_posix():sha(BUNDLE)})
    if DISPOSITION.exists():
        report["input_hashes"][DISPOSITION.relative_to(ROOT).as_posix()] = sha(DISPOSITION)
    report["implementation_hashes"].update(packet["implementation_hashes"])
    report["execution_seconds"]["new_actual_paired_constructor_values"] = result["execution_seconds"]
    report["scope"].update(slow_native_loop_seam_proved=result["passed"], modified_E_V_native_seam_proved=False,
        current_Rm_RQy_germ_identity_proved=False, full_goal_complete=False)
    report["current_turn_slow_native_value_seam_execution_audit"] = dict(
        new_paired_constructor_operator_once=True, paired_source_contract_saved_before_entry=True,
        accepted_derivative_provenance_consumed_as_data=True, accepted_math_operators_replayed=False,
        polynomial_expansion_or_GCD=False, old_unfinished_value_operator_restarted=False,
        original_attempt_bytes_preserved=True, saved_before_inspection=True,
        worker_fixed_role="luna_max_worker", worker_fork_turns="none", runtime_worker_model_metadata=None,
        numerical_field_or_quadrature_or_norm_computed=False, full_goal_complete=False)
    return report


def write_report():
    parent, _, _, packet = setup()
    result = read_saved(RESULT)
    validate_saved_result(result, packet)
    write_new(REPORT, build_report(parent, packet, result))
    return dict(passed=True, mathematical_slow_seam_passed=result["passed"],
        bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, _, _, packet = setup()
    report, result = read_saved(REPORT), read_saved(RESULT)
    validate_saved_result(result, packet)
    require(report == build_report(parent, packet, result), "serialized report content or scope changed")
    for name, expected in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "saved source pin changed: " + name)
    return dict(passed=True, serialized_only=True, mathematical_slow_seam_passed=result["passed"],
        bytes=REPORT.stat().st_size, sha256=sha(REPORT), accepted_or_new_operator_replayed=False,
        full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "structural-values", "inspect-result", "write-report", "inspect-report"):
        group.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        answer = dict(passed=True, AST_LF_EOF_only=True)
    else:
        answer = globals()[next(name for name, value in vars(args).items() if value)]()
    print(json.dumps(answer))
