"""Checkpoint new physical ordinary-Z source and derivative obligations.

Existing source functions and proofs are read as data. Candidate jet ports
never become derivative premises without a passed primitive source receipt.
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
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
PARENT = PREFIX + "Rm_physical_Qy_Z0_access_values.json"
PARENT_HASH = "4b8a2c9333939d330d04e4d39d9372644102f008d7a533836f5f14a4bca4a6a9"
GRAPH_HASH = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
LOADER = "scripts/diagnose_original_Rm_physical_qy0_access.py"
LOADER_HASH = "765667cb242c4dbb21fb7ef5ed030b6d4d6889b48136cbff5cf8475fd8967fdb"
RUNNER = "scripts/run_original_Rm_physical_Qy_Z_provenance.py"
KERNELS = dict(primitive="src/openai_ns_reconstruction/original_Rm_physical_primitive_Z_jets.py",
    outer="src/openai_ns_reconstruction/original_Rm_physical_qy_Z_derivative_provenance.py")
RECIPE_HASHES = {
    PREFIX + "modified_Q4_pressure_C3_functions.py": "7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7",
    PREFIX + "true_Rm_patch_finite_N_C5_functions.py": "c4445ef8acfb33d29cf9022862f8b6396e6aa2f38bda7467409b5378a5deaf15",
    PREFIX + "core_inlet_axial6_generic_source.py": "29738d46c6eb5897f0fbeefaec741d0eda4b97baed57b43ce9a7890a7dc03890",
    PREFIX + "Rm_outer_mixed_source_joins_check.json": "14d1cc6b641435c60a7e3662c996a260c07b328df5b8f486c4344f3d7eeb4867",
}
CHECKPOINTS = ROOT / "tmp/original_Rm_physical_Qy_Z_provenance"
SOURCE = CHECKPOINTS / "new_source.json"
RESULTS = {key:CHECKPOINTS / ("new_"+key+"_result.json") for key in KERNELS}
ENTRIES = {key:CHECKPOINTS / ("new_"+key+"_entered.json") for key in KERNELS}
REPORT = BASE / (PREFIX + "Rm_physical_Qy_Z_derivative_provenance.json")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()).hexdigest()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_pinned(path, expected):
    raw = path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, "source pin changed: " + str(path))
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
    require(sha(ROOT / LOADER) == LOADER_HASH, "saved-data loader changed")
    for name in (RUNNER, *KERNELS.values()):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF changed")
        ast.parse(raw, filename=name)
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    require(parent["actual_Rm_physical_Qy_Z0_value_seam"]["passed"] is True, "accepted Qy0 function identity missing")
    require(parent["scope"]["physical_Q_y_Z1_value_seam_proved"] is False, "Qy1 already admitted; do not repeat")
    loader = load_module("new_Qy_Z_saved_input_loader", ROOT / LOADER)
    nodes, old_inputs, _ = loader.load_inputs()
    require(digest(nodes) == GRAPH_HASH and len(nodes) == 205300, "saved graph changed")
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(loader.path_for(name)) == expected, "inherited pin changed: " + name)
    recipes = {}
    for name, expected in RECIPE_HASHES.items():
        raw = (BASE / name).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == expected, "producer source changed")
        recipes[name] = dict(path="experiments/root_st073/"+name, sha256=expected, source_text=raw.decode("utf-8"))
    primitive = dict(actual_source_bindings=old_inputs["actual_source_bindings"], source_recipe_records=recipes,
        expected_nodes_digest=GRAPH_HASH)
    outer = dict(qy0_value_seam=parent["actual_Rm_physical_Qy_Z0_value_seam"],
        common_phase_source=old_inputs["bounded_loop_ports_admission"]["common_phase_source"], expected_nodes_digest=GRAPH_HASH)
    packet = dict(schema="actual_Rm_physical_Qy_Z_provenance_inputs_v1", parent_report=PARENT,
        parent_sha256=PARENT_HASH, graph_digest=GRAPH_HASH, graph_node_count=len(nodes),
        input_digests=dict(primitive=digest(primitive), outer_without_primitive=digest(outer)),
        recipe_hashes=RECIPE_HASHES, implementation_hashes={name:sha(ROOT / name) for name in (RUNNER, LOADER, *KERNELS.values())},
        accepted_operators_replayed=False, numerical_field_evaluated=False,
        candidate_source_jets_are_not_admitted_by_array_order=True, full_goal_complete=False)
    return parent, nodes, dict(primitive=primitive, outer=outer), packet


def prepare():
    *_, packet = setup()
    require(not any(p.exists() for p in (SOURCE, *ENTRIES.values(), *RESULTS.values())), "new attempt already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, sha256=sha(SOURCE), input_digests=packet["input_digests"])


def actual_inputs(key, inputs):
    values = dict(inputs[key])
    dependencies = {}
    if key == "outer":
        receipt = json.loads(RESULTS["primitive"].read_bytes())
        require(receipt.get("schema") == "original_Rm_physical_primitive_Z_jets_v1", "primitive source receipt missing")
        validate_saved("primitive", receipt, inputs, json.loads(SOURCE.read_bytes()))
        # False primitive admission is allowed only for a conditional outer
        # calculation. The result contract below forbids promoting its lift.
        values["primitive_Z_source_receipt"] = receipt
        dependencies[RESULTS["primitive"].name] = sha(RESULTS["primitive"])
    return values, dependencies


def validate_result(key, result, values):
    require(isinstance(result, dict) and type(result.get("passed")) is bool, "malformed new result")
    require(result.get("source_only") is True and result.get("nodes_digest") == GRAPH_HASH, "source/graph result binding changed")
    if key == "primitive":
        require(result.get("schema") == "original_Rm_physical_primitive_Z_jets_v1" and result.get("Z_variable") == 20,
            "primitive receipt schema/context changed")
        rows = result.get("rows", [])
        require(len(rows) == 6, "primitive receipt needs all six source obligations")
        admitted = all(r.get("source_derivative_proved") is True and r.get("regularity_proved") is True for r in rows)
        require(not result["passed"] or admitted, "candidate primitive source rows promoted without proof")
    else:
        require(result.get("schema") == "original_Rm_physical_qy_Z_derivative_provenance_v1", "outer schema changed")
        require(type(result.get("outer_calculus_proved")) is bool
            and result.get("physical_Qy_Z1_value_seam_proved") == result["passed"], "outer admission inconsistent")
        require(result.get("full_RQy_germ_proved") is False and result.get("old_operators_replayed") is False
            and result.get("numeric_values_computed") is False, "outer execution/scope changed")
        require(not result["passed"] or (values["primitive_Z_source_receipt"]["passed"] is True
            and result["outer_calculus_proved"] is True and not result.get("unresolved_gates")),
            "conditional outer proof promoted without primitive source provenance")


def run_component(key):
    _, nodes, inputs, packet = setup()
    require(json.loads(SOURCE.read_bytes()) == packet, "frozen foundation changed")
    require(not ENTRIES[key].exists() and not RESULTS[key].exists(), "new component already entered")
    values, dependencies = actual_inputs(key, inputs)
    source_path = CHECKPOINTS / ("new_"+key+"_source.json")
    source = dict(foundation_sha256=sha(SOURCE), component=key, input_digest=digest(values),
        dependency_hashes=dependencies, graph_digest=GRAPH_HASH, kernel_sha256=packet["implementation_hashes"][KERNELS[key]])
    write_new(source_path, source)
    module = load_module("new_Rm_Qy_Z_"+key, ROOT / KERNELS[key])
    write_new(ENTRIES[key], dict(source_sha256=sha(source_path), kernel_sha256=source["kernel_sha256"],
        input_digest=source["input_digest"], entered_at_utc=datetime.now(timezone.utc).isoformat()))
    started = time.perf_counter()
    result = None
    try:
        result = module.derive(nodes, **values)
        require(digest(nodes) == GRAPH_HASH and digest(values) == source["input_digest"], "new operator mutated source input")
        validate_result(key, result, values)
    except Exception as error:
        result = dict(passed=False, failed_attempt_preserved=True, error_type=type(error).__name__, error=str(error),
            traceback="".join(traceback.format_exception(error)), rejected_return=result)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(source_path),
        entry_checkpoint_sha256=sha(ENTRIES[key]), runner_full_goal_complete=False)
    write_new(RESULTS[key], result)
    return dict(component=key, passed=result["passed"], sha256=sha(RESULTS[key]), seconds=result["execution_seconds"],
        outer_calculus_proved=result.get("outer_calculus_proved"), error=result.get("error"),
        unresolved_gate_count=len(result.get("unresolved_gates", [])))


def validate_saved(key, result, inputs, packet):
    values, dependencies = actual_inputs(key, inputs)
    source_path = CHECKPOINTS / ("new_"+key+"_source.json")
    source, entry = json.loads(source_path.read_bytes()), json.loads(ENTRIES[key].read_bytes())
    require(source["foundation_sha256"] == sha(SOURCE) and source["input_digest"] == digest(values)
        and source["dependency_hashes"] == dependencies and source["graph_digest"] == GRAPH_HASH
        and source["kernel_sha256"] == packet["implementation_hashes"][KERNELS[key]], "saved component source changed")
    require(entry["source_sha256"] == result["source_checkpoint_sha256"] == sha(source_path)
        and entry["input_digest"] == source["input_digest"] and entry["kernel_sha256"] == source["kernel_sha256"]
        and result["entry_checkpoint_sha256"] == sha(ENTRIES[key]), "saved component entry/result binding changed")
    if result.get("failed_attempt_preserved") is True:
        require(result["passed"] is False and isinstance(result.get("traceback"), str), "malformed failed attempt")
    else:
        validate_result(key, result, values)


def build_report(parent, packet, results):
    report = copy.deepcopy(parent)
    outer = results["outer"]
    report.update(stage="actual_Rm_physical_Qy_ordinary_Z_derivative_provenance", parent_report=PARENT,
        parent_report_sha256=PARENT_HASH, actual_Rm_physical_Qy_Z_provenance_inputs=packet,
        actual_Rm_physical_primitive_Z_jets=results["primitive"],
        actual_Rm_physical_Qy_Z_derivative_provenance=outer,
        physical_Qy_Z_provenance_checkpoint_hashes={p.name:sha(p) for p in sorted(CHECKPOINTS.glob("*.json"))})
    report["input_hashes"].update({PARENT:PARENT_HASH, **RECIPE_HASHES})
    report["implementation_hashes"].update(packet["implementation_hashes"])
    for key, result in results.items():
        report["execution_seconds"]["new_physical_Qy_Z_"+key] = result["execution_seconds"]
    report["scope"].update(physical_Q_y_Z1_value_seam_proved=outer["passed"],
        current_Rm_RQy_germ_identity_proved=False, full_goal_complete=False)
    report["current_turn_physical_Qy_Z_execution_audit"] = dict(primitive_source_operator_calls=1,
        outer_calculus_operator_calls=1, accepted_mathematical_operators_replayed=False,
        numerical_field_or_quadrature_or_norm_computed=False, source_and_entry_saved_before_each_call=True,
        conditional_calculus_is_not_unconditional_Qy1_admission=True, full_goal_complete=False)
    return report


def saved_setup():
    parent, _, inputs, packet = setup()
    require(json.loads(SOURCE.read_bytes()) == packet, "saved foundation changed")
    results = {key:json.loads(path.read_bytes()) for key,path in RESULTS.items()}
    for key, result in results.items():
        validate_saved(key, result, inputs, packet)
    return parent, packet, results


def write_report():
    parent, packet, results = saved_setup()
    write_new(REPORT, build_report(parent, packet, results))
    return dict(saved=True, sha256=sha(REPORT), bytes=REPORT.stat().st_size,
        primitive_sources_passed=results["primitive"]["passed"], Qy_Z1_passed=results["outer"]["passed"],
        outer_calculus_proved=results["outer"].get("outer_calculus_proved"))


def inspect():
    parent, packet, results = saved_setup()
    require(json.loads(REPORT.read_bytes()) == build_report(parent, packet, results), "serialized report content/scope changed")
    return dict(passed=True, serialized_only=True, sha256=sha(REPORT), bytes=REPORT.stat().st_size,
        primitive_sources_passed=results["primitive"]["passed"], Qy_Z1_passed=results["outer"]["passed"],
        outer_calculus_proved=results["outer"].get("outer_calculus_proved"), accepted_operators_replayed=False,
        full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "primitive", "outer", "write_report", "inspect"))
    action = parser.parse_args().action
    print(json.dumps(run_component(action) if action in KERNELS else globals()[action]()))
