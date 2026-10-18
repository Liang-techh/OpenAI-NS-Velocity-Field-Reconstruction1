"""Checkpoint the corrected physical Q_y Z0 consumer only.

The current diagnostic's loader reads saved data and authenticates archives;
neither the accepted radial nor slow mathematical consumer is called.
"""
from __future__ import annotations

import argparse
import ast
import copy
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PARENT = "lei_ren_part1_paper_compliant_current_original_Rm_modified_physical_values.json"
PARENT_HASH = "5c26b7bb3b70e874a5790a5f8424233678b5a964200b8ffd88dce102b22110f1"
GRAPH_HASH = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
INPUT_HASH = "081b95e2df85e4f9573d31277db95e009e181a08ab6d84f0a02703b407ea1195"
LOADER = "scripts/diagnose_original_Rm_physical_qy0_access.py"
LOADER_HASH = "765667cb242c4dbb21fb7ef5ed030b6d4d6889b48136cbff5cf8475fd8967fdb"
DIAGNOSTIC = "experiments/root_st073/original_Rm_physical_Qy_Z0_traceback_attempt.json"
DIAGNOSTIC_HASH = "dec29dbaea2dd4498b841c74744085eefd552003f95ccaa3a87c6094e25770c1"
RUNNER = "scripts/run_original_Rm_physical_qy0_access_values.py"
KERNEL = "src/openai_ns_reconstruction/original_Rm_physical_qy0_access_values.py"
KERNEL_HASH = "b4196776aabc9dae09ef83dfbeb5e274441f02897592fa33be8757677760f1c3"  # New internal moment-map access correction.
CHECKPOINTS = ROOT / "tmp/original_Rm_physical_qy0_access_values"
SOURCE = CHECKPOINTS / "new_source.json"
ENTRY = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_physical_values.json"
REPORT = BASE / "lei_ren_part1_paper_compliant_current_original_Rm_physical_Qy_Z0_access_values.json"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_pinned(path, expected):
    raw = path.read_bytes()
    require(isinstance(expected, str) and hashlib.sha256(raw).hexdigest() == expected, "source pin changed: " + str(path))
    return json.loads(raw)


def write_new(path, value):
    require(not path.exists(), "artifact already exists: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True)
        stream.write("\n")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def setup():
    require(sha(ROOT / LOADER) == LOADER_HASH and sha(ROOT / KERNEL) == KERNEL_HASH, "new code is not pinned")
    for name in (RUNNER, KERNEL, LOADER):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF changed")
        ast.parse(raw, filename=name)
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    diagnostic = read_pinned(ROOT / DIAGNOSTIC, DIAGNOSTIC_HASH)
    require(diagnostic["exception_captured"] is True and diagnostic["error_type"] == "KeyError"
        and diagnostic["error"] == "'postrestore'" and diagnostic["input_unchanged"] is True
        and diagnostic["graph_unchanged"] is True, "preserved diagnostic does not show the access failure")
    fault = diagnostic["frames"][-1]
    require(fault["line"] == 978 and fault["function"] == "<listcomp>" and fault["side"] == "postrestore"
        and "'sides'" in fault["mapping_key_summaries"]["moment_row"]["keys"], "first failed access changed")
    loader = load_module("new_Rm_saved_input_loader", ROOT / LOADER)
    nodes, inputs, _ = loader.load_inputs()
    require(digest(inputs) == INPUT_HASH and digest(nodes) == GRAPH_HASH, "corrected attempt changed physical inputs")
    packet = dict(schema="Rm_corrected_physical_Qy_Z0_access_inputs_v1", parent_report=PARENT,
        parent_sha256=PARENT_HASH, graph_digest=GRAPH_HASH, input_digest=INPUT_HASH,
        diagnostic=DIAGNOSTIC, diagnostic_sha256=DIAGNOSTIC_HASH,
        implementation_hashes={name:sha(ROOT / name) for name in (RUNNER, KERNEL, LOADER)},
        correction="four internal moment_row side accesses use the validated sides map",
        accepted_radial_and_slow_consumed_as_data=True, accepted_math_operators_replayed=False,
        physical_Qy_Z1_admitted=False, full_goal_complete=False)
    return parent, nodes, inputs, packet


def prepare():
    *_, packet = setup()
    require(not any(p.exists() for p in (SOURCE, ENTRY, RESULT)), "corrected attempt already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, sha256=sha(SOURCE), input_digest=INPUT_HASH)


def validate_kernel_result(result):
    require(isinstance(result, dict) and type(result.get("passed")) is bool, "malformed new result")
    require(result.get("schema") == "original_Rm_physical_qy0_value_seam_v1" and result.get("source_only") is True,
        "corrected kernel schema/source contract changed")
    for key in ("numeric_values_computed", "pressure_quadrature_executed", "old_mathematical_operator_replayed",
        "accepted_graph_replayed", "differentiated_opaque_alias", "full_RQy_germ_Z1_proved", "physical_RQy_Z1_proved"):
        require(result.get(key) is False, "forbidden or missing kernel action/scope flag: " + key)
    require(result.get("expected_nodes_digest") == result.get("nodes_digest") == GRAPH_HASH, "kernel graph binding changed")
    rows = result.get("canonical_value_pairs", [])
    require(len(rows) == 1 and rows[0].get("field") == "Q_y" and rows[0].get("ordinary_order") == 0
        and rows[0].get("postrestore_target_handle") == 183451 and rows[0].get("patch_target_handle") == 204500,
        "corrected result targets changed")
    require(rows[0].get("equal") == result["passed"] and rows[0].get("numeric_value") is False,
        "physical admission flag inconsistent")
    require(not result["passed"] or rows[0].get("status") == "proved", "physical row was not proved")
    z1 = result.get("q_y_z1", {})
    require(z1.get("ordinary_order") == 1 and z1.get("status") == "open" and z1.get("equal") is False
        and z1.get("proved") is False and z1.get("postrestore_target_handle") == 183455
        and z1.get("patch_target_handle") == 204504, "open ordinary-Z1 scope changed")
    require(result.get("input_data_digests", {}).get("actual_source_bindings") == digest(result["physical_ports"]),
        "checked physical port binding changed")


def run():
    _, nodes, inputs, packet = setup()
    require(json.loads(SOURCE.read_bytes()) == packet, "corrected source changed")
    require(not ENTRY.exists() and not RESULT.exists(), "corrected consumer already entered")
    module = load_module("new_corrected_Rm_physical_Qy_Z0", ROOT / KERNEL)
    write_new(ENTRY, dict(source_sha256=sha(SOURCE), kernel_sha256=KERNEL_HASH,
        input_digest=INPUT_HASH, entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started = time.perf_counter()
    result = None
    try:
        result = module.derive(nodes, **inputs)
        require(digest(nodes) == GRAPH_HASH and digest(inputs) == INPUT_HASH, "new operator mutated immutable input")
        validate_kernel_result(result)
        for key in ("bounded_loop_ports_admission", "modified_moment_germ", "modified_radial_value_seam"):
            require(result["input_data_digests"][key] == digest(inputs[key]), "new packet digest changed: " + key)
    except Exception as error:
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__,
            error=str(error), traceback="".join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
        entry_checkpoint_sha256=sha(ENTRY), runner_full_goal_complete=False)
    write_new(RESULT, result)
    return dict(passed=result["passed"], sha256=sha(RESULT), seconds=result["execution_seconds"], error=result.get("error"))


def validate_saved(result, packet):
    require(json.loads(SOURCE.read_bytes()) == packet, "saved corrected source changed")
    entry = json.loads(ENTRY.read_bytes())
    require(entry["source_sha256"] == result["source_checkpoint_sha256"] == sha(SOURCE)
        and entry["kernel_sha256"] == KERNEL_HASH and entry["input_digest"] == INPUT_HASH
        and result["entry_checkpoint_sha256"] == sha(ENTRY), "saved corrected entry/result binding changed")
    if result.get("failed_attempt_preserved") is True:
        require(result["passed"] is False and isinstance(result.get("traceback"), str), "malformed preserved failure")
    else:
        validate_kernel_result(result)


def build_report(parent, packet, result):
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_corrected_physical_Qy_Z0_access_values", parent_report=PARENT,
        parent_report_sha256=PARENT_HASH, actual_Rm_corrected_physical_Qy_Z0_inputs=packet,
        actual_Rm_physical_Qy_Z0_value_seam=result,
        corrected_physical_Qy_Z0_checkpoint_hashes={p.name:sha(p) for p in (SOURCE, ENTRY, RESULT)},
        original_physical_Qy_Z0_traceback_attempt=dict(path=DIAGNOSTIC, sha256=DIAGNOSTIC_HASH,
            first_fault_line=978, prior_unaccepted_physical_reentry_for_traceback_once=True))
    report["input_hashes"].update({PARENT:PARENT_HASH, DIAGNOSTIC:DIAGNOSTIC_HASH})
    report["implementation_hashes"].update(packet["implementation_hashes"])
    report["execution_seconds"]["new_corrected_physical_Qy_Z0"] = result["execution_seconds"]
    report["scope"].update(physical_Q_y_Z0_value_seam_proved=result["passed"],
        physical_Q_y_Z1_value_seam_proved=False, current_Rm_RQy_germ_identity_proved=False, full_goal_complete=False)
    report["current_turn_corrected_physical_Qy_Z0_execution_audit"] = dict(new_corrected_consumer_calls=1,
        failed_unaccepted_consumer_reentries_for_traceback=1, accepted_radial_or_slow_replayed=False,
        source_and_entry_saved_before_operator=True, saved_before_inspection=True,
        numerical_field_or_quadrature_or_norm_computed=False, full_goal_complete=False)
    return report


def write_report():
    parent, _, _, packet = setup()
    result = json.loads(RESULT.read_bytes())
    validate_saved(result, packet)
    write_new(REPORT, build_report(parent, packet, result))
    return dict(passed=True, physical_Qy_Z0_passed=result["passed"], bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect():
    parent, _, _, packet = setup()
    result = json.loads(RESULT.read_bytes())
    validate_saved(result, packet)
    report = json.loads(REPORT.read_bytes())
    require(report == build_report(parent, packet, result), "serialized corrected report changed")
    return dict(passed=True, serialized_only=True, physical_Qy_Z0_passed=result["passed"],
        report_sha256=sha(REPORT), report_bytes=REPORT.stat().st_size, accepted_operators_replayed=False,
        full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "write_report", "inspect"))
    print(json.dumps(globals()[parser.parse_args().action]()))
