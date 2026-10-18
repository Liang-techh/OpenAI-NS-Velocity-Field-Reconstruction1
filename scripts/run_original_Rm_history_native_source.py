"""Install actual Rm Duhamel endpoint/native-y sources without density boxing."""
from __future__ import annotations

import argparse
import ast
import copy
import importlib.util
import json
import sys
import time

from run_original_Rm_second_radial_identity import (
    ROOT, BASE, PREFIX, require, sha, digest, path_for, read_pinned, read_saved, write_new,
)

PARENT = PREFIX + "Rm_second_radial_identity.json"
PARENT_HASH = "f221e697094a0b290b0674b11bcca534c688e0982ca775ffe923a42c839afab9"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
INLET_RECEIPT = PREFIX + "true_Rm_patch_finite_N_C5_functions_check.json"
KERNEL = "src/openai_ns_reconstruction/original_Rm_history_native_source.py"
KERNEL_HASH = "76ed0808975ad4ea83c1ba572f986d9d2e911f88d029600bac25769eba14a4d2"
FIXTURE = "experiments/root_st073/original_Rm_history_native_source_fixture.json"
FIXTURE_HASH = "bcc62dc012579573aa75d19c86bb1aea8aef543939fc17f58f1a9ffd9056cf09"
RUNNER = "scripts/run_original_Rm_history_native_source.py"
REPORT = BASE / (PREFIX + "Rm_history_native_source.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_history_native_source"
SOURCE = CHECKPOINTS / "new_source.json"
ENTERED = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_history_native_source.json"


def preflight():
    for relative in (RUNNER, KERNEL):
        raw = (ROOT / relative).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: "+relative)
        ast.parse(raw, filename=relative)
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "history adapter changed")
    fixture = read_pinned(ROOT / FIXTURE, FIXTURE_HASH)
    require(fixture.get("passed") is True or fixture.get("toy_harness", {}).get("passed") is True,
            "new history fixture not admitted")


def setup():
    preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: "+name)
    finite = read_pinned(BASE / FINITE, parent["input_hashes"][FINITE])
    receipt = read_pinned(BASE / INLET_RECEIPT, parent["input_hashes"][INLET_RECEIPT])
    nodes = finite["exact_graph_nodes"]
    first, second = parent["actual_Rm_native_radial_endpoint_source"], parent["actual_Rm_second_radial_endpoint_functions"]
    require(digest(nodes) == first["base_digest"], "base source changed")
    nodes = nodes + copy.deepcopy(first["appended_nodes"])
    require(digest(nodes) == first["full_nodes_digest"] == second["base_digest"], "first endpoint delta changed")
    nodes = nodes + copy.deepcopy(second["appended_nodes"])
    require(len(nodes) == 205160 and digest(nodes) == second["full_nodes_digest"], "second endpoint delta changed")
    fb = finite["Rm_patch_finite_N_C5_bindings"]
    packet = dict(schema="actual_Rm_history_native_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, graph_node_count=len(nodes), graph_digest=digest(nodes),
        finite_sha256=parent["input_hashes"][FINITE], old_receipt_sha256=parent["input_hashes"][INLET_RECEIPT],
        old_receipt_provenance_only=True, no_old_toy_field_identity_admitted=True,
        coordinate=fb["chart"]["coordinate"], physical_measure=fb["chart"]["physical_measure"],
        history=copy.deepcopy(fb["chart"]["history"]), incoming=copy.deepcopy(fb["chart"]["incoming"]),
        ordinary_Z_orders=[0, 1], actual_new_endpoint_and_FTC_source_required=True,
        all_old_212_terms_and_errors_preserved=True)
    return parent, nodes, fb, receipt, packet


def prepare_source():
    *_, packet = setup()
    require(not any(path.exists() for path in (SOURCE, ENTERED, RESULT)), "source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE), graph_node_count=packet["graph_node_count"])


def compute_new():
    parent, nodes, fb, receipt, packet = setup()
    require(read_saved(SOURCE) == packet and not ENTERED.exists() and not RESULT.exists(), "source changed or entered; no replay")
    spec = importlib.util.spec_from_file_location("new_Rm_history_native_source", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTERED, dict(source_sha256=sha(SOURCE), kernel_sha256=KERNEL_HASH))
    started = time.perf_counter()
    result = module.derive(nodes, {"Rm_patch_finite_N_C5_bindings":fb}, receipt,
                           expected_nodes_digest=packet["graph_digest"], ordinary_orders=(0, 1))
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
        entry_checkpoint_sha256=sha(ENTERED), actual_saved_graph_run=True, actual_saved_graph_loaded=True, source_executed=True)
    result["scope"]["no_saved_graph_execution"] = False
    result["scope"]["accepted_ancestor_operator_executed"] = False
    write_new(RESULT, result)
    return dict(passed=result.get("passed"), result_sha256=sha(RESULT),
                new_graph_nodes=result["graph_delta"]["appended_node_count"])


def write_report():
    parent, nodes, fb, receipt, expected = setup()
    packet, result = read_saved(SOURCE), read_saved(RESULT)
    require(packet == expected and result["passed"] is True, "new history source not admitted")
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_history_native_source", parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        actual_Rm_history_native_inputs=packet, actual_Rm_history_native_source=result,
        history_native_checkpoint_hashes={path.name:sha(path) for path in (SOURCE, ENTERED, RESULT)})
    report["execution_seconds"]["new_actual_history_native_source"] = result["execution_seconds"]
    report["input_hashes"].update({PARENT:PARENT_HASH, FIXTURE:FIXTURE_HASH})
    report["implementation_hashes"].update({KERNEL:KERNEL_HASH, RUNNER:sha(ROOT / RUNNER)})
    report["scope"].update(actual_current_Rm_history_endpoint_identity=True, actual_current_Rm_history_native_source=True,
        remaining_density_seam_unproved=True, modified_finite_N_germ=False, current_Rm_RQy_germ_identity_proved=False,
        current_Rm_to_Rp_join=False, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, fb, receipt, expected = setup()
    report, result = read_saved(REPORT), read_saved(RESULT)
    require(report["actual_Rm_history_native_inputs"] == read_saved(SOURCE) == expected
            and report["actual_Rm_history_native_source"] == result, "saved history source changed")
    require(report["history_native_checkpoint_hashes"] == {path.name:sha(path) for path in (SOURCE, ENTERED, RESULT)}
            and result["source_checkpoint_sha256"] == sha(SOURCE) and result["entry_checkpoint_sha256"] == sha(ENTERED),
            "history checkpoint link changed")
    delta = result["graph_delta"]
    require(delta["base_graph_node_count"] == len(nodes) and delta["base_nodes_sha256"] == digest(nodes)
            and digest(nodes+delta["appended_nodes"]) == delta["full_canonical_nodes_sha256"], "history graph delta changed")
    for key in ("actual_Rm_second_radial_identity", "actual_Rm_second_radial_endpoint_functions", "actual_Rm_leading_seed_admission",
                "current_Rm_pressure_inlet_identity", "actual_Rm_absolute_m15", "corrected_own_swirl_provider", "corrected_reserved_point_numeric"):
        require(report[key] == parent[key], "inherited packet changed: "+key)
    for name, expected_hash in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected_hash, "report pin changed: "+name)
    require(result["remaining_density_seam_unproved"] is True and report["scope"]["full_goal_complete"] is False,
            "invalid density/full germ claim")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
        new_native_history_sources=len(result["rows"]), old_mathematics_executed=False, new_operator_replayed=False,
        density_seam_proved=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "compute-new", "write-report", "inspect-report"):
        group.add_argument("--"+name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        print(json.dumps(dict(passed=True, source_preflight_only=True)))
    else:
        action = prepare_source if args.prepare_source else compute_new if args.compute_new else write_report if args.write_report else inspect_report
        print(json.dumps(action(), indent=2))
