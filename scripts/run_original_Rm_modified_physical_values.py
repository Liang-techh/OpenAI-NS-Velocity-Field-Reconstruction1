"""Checkpoint new modified radial and physical Q_y Z0 value consumers.

Saved graph nodes and proof packets are loaded as immutable data. No accepted
producer, derivative traversal, recurrence, integral or value reducer is run.
"""
from __future__ import annotations

import argparse
import ast
import copy
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
PARENT = PREFIX + "Rm_slow_native_value_seam.json"
PARENT_HASH = None  # Filled only after the existing live result is admitted.
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
FINITE_HASH = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
MODIFIED = PREFIX + "modified_Q4_pressure_C3_functions.json.gz"
MODIFIED_HASH = "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"
PORTS = "experiments/root_st073/original_Rm_loop_ports_admission.json"
PORTS_HASH = "0331c5fc37bd3ec6c0d386c38b07d9863ede622e14f1c6dc65bf8cdf675a16b9"
GRAPH_HASH = "a2dcc782ba3033ff6ca5f54a66632ac6ef007bb0c5b9960f10567d44ef2059a3"
RUNNER = "scripts/run_original_Rm_modified_physical_values.py"
KERNELS = {
    "radial": "src/openai_ns_reconstruction/original_Rm_modified_radial_value_seam.py",
    "physical": "src/openai_ns_reconstruction/original_Rm_physical_qy0_value_seam.py",
}
CHECKPOINTS = ROOT / "tmp/original_Rm_modified_physical_values"
SOURCE = CHECKPOINTS / "new_source.json"
RESULTS = {key: CHECKPOINTS / ("new_" + key + "_values.json") for key in KERNELS}
ENTRIES = {key: CHECKPOINTS / ("new_" + key + "_entered.json") for key in KERNELS}
REPORT = BASE / (PREFIX + "Rm_modified_physical_values.json")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode()
    return hashlib.sha256(encoded).hexdigest()


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def path_for(name):
    path = ROOT / name if name.startswith(("src/", "scripts/", "experiments/", "docs/")) else BASE / name
    resolved = path.resolve()
    require(resolved.is_relative_to(ROOT.resolve()), "source pin escapes repository: " + name)
    return resolved


def read_saved(path):
    return json.loads(path.read_bytes())


def read_pinned(path, expected):
    raw = path.read_bytes()
    require(isinstance(expected, str) and hashlib.sha256(raw).hexdigest() == expected, "source pin changed: " + str(path))
    # Compressed bytes are authenticated before decompression.
    return json.loads(gzip.decompress(raw) if path.suffix == ".gz" else raw)


def write_new(path, value):
    require(not path.exists(), "checkpoint already exists: " + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, ensure_ascii=True)
        stream.write("\n")


def preflight():
    for name in (RUNNER, *KERNELS.values()):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + name)
        ast.parse(raw, filename=name)


def append_saved_delta(nodes, packet, base_hash_key, full_hash_key):
    require(digest(nodes) == packet[base_hash_key], "saved delta base changed")
    base_count = packet.get("base_node_count", packet.get("base_graph_node_count"))
    require(len(nodes) == base_count, "saved delta base count changed")
    result = nodes + copy.deepcopy(packet["appended_nodes"])
    require(len(result) == packet["full_node_count"] and digest(result) == packet[full_hash_key], "saved graph delta changed")
    return result


def load_saved_source(parent):
    finite = read_pinned(BASE / FINITE, FINITE_HASH)
    modified = read_pinned(BASE / MODIFIED, MODIFIED_HASH)
    require(finite["source_family"] == modified["source_family"] == parent["source_family"], "source family changed")
    nodes = finite["exact_graph_nodes"]
    require(len(nodes) == 205147 and nodes[:183700] == modified["exact_graph_nodes"], "actual original prefixes changed")
    for packet, base_key, full_key in (
        (parent["actual_Rm_native_radial_endpoint_source"], "base_digest", "full_nodes_digest"),
        (parent["actual_Rm_second_radial_endpoint_functions"], "base_digest", "full_nodes_digest"),
        (parent["actual_Rm_history_native_source"]["graph_delta"], "base_nodes_sha256", "full_canonical_nodes_sha256"),
        (parent["actual_Rm_history_seam_identity"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256"),
        (parent["actual_Rm_modified_moment_germ"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256"),
    ):
        nodes = append_saved_delta(nodes, packet, base_key, full_key)
    require(len(nodes) == 205300 and digest(nodes) == GRAPH_HASH, "current saved graph changed")
    bindings = {
        "postrestore": modified["modified_Q4_pressure_C3_bindings"]["charts"]["postrestore"],
        "patch": finite["Rm_patch_finite_N_C5_bindings"],
    }
    return nodes, bindings


def setup():
    preflight()
    require(PARENT_HASH is not None, "the already running slow-target VALUE result has not yet been admitted")
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "inherited source pin changed: " + name)
    require(parent["actual_Rm_slow_native_value_seam"]["passed"] is True, "slow-target VALUE seam must pass first")
    nodes, bindings = load_saved_source(parent)
    ports = read_pinned(ROOT / PORTS, PORTS_HASH)
    implementation = {name: sha(ROOT / name) for name in (RUNNER, *KERNELS.values())}
    # Kernel-specific input tables are finalized before this source is frozen.
    inputs = component_inputs(parent, ports, bindings)
    packet = dict(schema="actual_Rm_modified_physical_value_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, graph_node_count=len(nodes), graph_digest=GRAPH_HASH,
        input_digests={key: digest(value) for key, value in inputs.items()},
        source_bindings_digest=digest(bindings), implementation_hashes=implementation,
        archive_hashes={FINITE: FINITE_HASH, MODIFIED: MODIFIED_HASH, PORTS: PORTS_HASH},
        accepted_math_operators_replayed=False, derivatives_generated=False, numerical_field_evaluated=False,
        physical_Q_y_Z1_admitted=False, full_goal_complete=False)
    return parent, nodes, inputs, packet


def component_inputs(parent, ports, bindings):
    return {
        "radial": dict(generic_packet=parent["actual_Rm_generic_native_identity"],
            loop_ports_admission=ports, slow_value_seam=parent["actual_Rm_slow_native_value_seam"],
            raw_native_identity=parent["actual_Rm_raw_native_jet_identity"],
            second_radial_identity=parent["actual_Rm_second_radial_identity"],
            expected_nodes_digest=GRAPH_HASH),
        "physical": dict(bounded_loop_ports_admission=ports,
            modified_moment_germ=parent["actual_Rm_modified_moment_germ"],
            actual_source_bindings=bindings, expected_nodes_digest=GRAPH_HASH),
    }


def prepare_source():
    *_, packet = setup()
    require(not any(path.exists() for path in (SOURCE, *ENTRIES.values(), *RESULTS.values())), "source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, sha256=sha(SOURCE), nodes=packet["graph_node_count"])


def actual_inputs(key, inputs):
    value = dict(inputs[key])
    dependencies = {}
    if key == "physical":
        radial = read_saved(RESULTS["radial"])
        require(radial["passed"] is True, "new total E_y/V_y values must pass before physical Q_y Z0")
        require(radial["source_checkpoint_sha256"] == sha(CHECKPOINTS / "new_radial_source.json"), "radial source link changed")
        require(radial["entry_checkpoint_sha256"] == sha(ENTRIES["radial"]), "radial entry link changed")
        value["modified_radial_value_seam"] = radial
        dependencies[RESULTS["radial"].name] = sha(RESULTS["radial"])
    return value, dependencies


def run_component(key):
    _, nodes, inputs, packet = setup()
    require(read_saved(SOURCE) == packet, "new frozen foundation changed")
    require(not ENTRIES[key].exists() and not RESULTS[key].exists(), "new component already entered")
    values, dependencies = actual_inputs(key, inputs)
    component_source = CHECKPOINTS / ("new_" + key + "_source.json")
    source_packet = dict(foundation_sha256=sha(SOURCE), component=key, input_digest=digest(values),
        dependency_checkpoint_hashes=dependencies, graph_digest=GRAPH_HASH,
        kernel_sha256=packet["implementation_hashes"][KERNELS[key]],
        accepted_math_operators_replayed=False, numerical_field_evaluated=False)
    write_new(component_source, source_packet)
    sys.path.insert(0, str(ROOT / "src"))
    spec = importlib.util.spec_from_file_location("new_Rm_" + key + "_value_consumer", ROOT / KERNELS[key])
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTRIES[key], dict(component_source_sha256=sha(component_source),
        kernel_sha256=source_packet["kernel_sha256"], input_digest=source_packet["input_digest"]))
    started = time.perf_counter()
    try:
        result = module.derive(nodes, **values)
    except Exception as error:
        result = dict(passed=False, error_type=type(error).__name__, error=str(error), failed_attempt_preserved=True)
    result.update(execution_seconds=time.perf_counter() - started,
        source_checkpoint_sha256=sha(component_source), entry_checkpoint_sha256=sha(ENTRIES[key]),
        actual_source_graph_inspected=True, accepted_math_operators_replayed=False,
        numeric_field_evaluated=False, full_goal_complete=False)
    write_new(RESULTS[key], result)
    return dict(component=key, passed=result["passed"], seconds=result["execution_seconds"],
        result_sha256=sha(RESULTS[key]), error=result.get("error"),
        unresolved_gate_count=len(result.get("unresolved_gates", [])))


def checkpoint_hashes():
    paths = [SOURCE]
    for key in KERNELS:
        paths.extend((CHECKPOINTS / ("new_" + key + "_source.json"), ENTRIES[key], RESULTS[key]))
    return {path.name: sha(path) for path in paths if path.exists()}


def write_report():
    parent, _, _, packet = setup()
    require(read_saved(SOURCE) == packet, "new source checkpoint changed")
    radial = read_saved(RESULTS["radial"])
    if RESULTS["physical"].exists():
        physical = read_saved(RESULTS["physical"])
    else:
        require(radial["passed"] is False and not ENTRIES["physical"].exists(), "missing physical component result")
        physical = dict(passed=False, not_run=True, reason="total radial VALUE seam did not pass")
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_modified_radial_and_physical_Qy_Z0_values",
        parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        actual_Rm_modified_physical_value_inputs=packet,
        actual_Rm_modified_radial_value_seam=radial, actual_Rm_physical_Qy_Z0_value_seam=physical,
        modified_physical_value_checkpoint_hashes=checkpoint_hashes())
    report["input_hashes"].update({PARENT: PARENT_HASH, FINITE: FINITE_HASH, MODIFIED: MODIFIED_HASH, PORTS: PORTS_HASH})
    report["implementation_hashes"].update(packet["implementation_hashes"])
    report["execution_seconds"]["new_actual_modified_radial_values"] = radial["execution_seconds"]
    if not physical.get("not_run"):
        report["execution_seconds"]["new_actual_physical_Qy_Z0_values"] = physical["execution_seconds"]
    report["scope"].update(modified_E_V_native_seam_proved=radial["passed"],
        physical_Q_y_Z0_value_seam_proved=physical["passed"], physical_Q_y_Z1_value_seam_proved=False,
        current_Rm_RQy_germ_identity_proved=False, full_goal_complete=False)
    report["current_turn_modified_physical_values_execution_audit"] = dict(
        new_radial_value_operator_once=True, new_physical_Qy_Z0_value_operator_once=not physical.get("not_run", False),
        accepted_packets_consumed_as_data=True, accepted_math_operators_replayed=False,
        derivatives_generated=False, numerical_field_or_quadrature_or_norm_computed=False,
        source_and_entry_saved_before_each_new_operator=True, worker_fixed_role="luna_max_worker",
        worker_fork_turns="none", runtime_worker_model_metadata=None, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=True, radial_passed=radial["passed"], physical_Qy_Z0_passed=physical["passed"],
        bytes=REPORT.stat().st_size, sha256=sha(REPORT), physical_Qy_Z1_admitted=False)


def inspect_report():
    parent, _, inputs, packet = setup()
    report = read_saved(REPORT)
    require(report["actual_Rm_modified_physical_value_inputs"] == read_saved(SOURCE) == packet, "frozen source packet changed")
    require(report["modified_physical_value_checkpoint_hashes"] == checkpoint_hashes(), "saved checkpoint bytes changed")
    for key, report_key in (("radial", "actual_Rm_modified_radial_value_seam"), ("physical", "actual_Rm_physical_Qy_Z0_value_seam")):
        if not RESULTS[key].exists():
            require(key == "physical" and report[report_key].get("not_run") is True, "unexpected missing result")
            continue
        result = read_saved(RESULTS[key])
        source_path = CHECKPOINTS / ("new_" + key + "_source.json")
        source, entry = read_saved(source_path), read_saved(ENTRIES[key])
        values, dependencies = actual_inputs(key, inputs)
        require(source["foundation_sha256"] == sha(SOURCE) and source["input_digest"] == digest(values), "component input data changed")
        require(source["dependency_checkpoint_hashes"] == dependencies, "component dependency changed")
        require(entry["component_source_sha256"] == result["source_checkpoint_sha256"] == sha(source_path), "component source link changed")
        require(result["entry_checkpoint_sha256"] == sha(ENTRIES[key]) and report[report_key] == result, "component saved result changed")
    changed_keys = {"stage", "parent_report", "parent_report_sha256", "input_hashes", "implementation_hashes", "execution_seconds", "scope"}
    for key, value in parent.items():
        if key not in changed_keys:
            require(report[key] == value, "inherited packet changed: " + key)
    for name, expected in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "saved source pin changed: " + name)
    require(not report["scope"]["physical_Q_y_Z1_value_seam_proved"]
        and not report["scope"]["current_Rm_RQy_germ_identity_proved"]
        and not report["scope"]["full_goal_complete"], "open physical/global gate changed")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
        radial_passed=report["actual_Rm_modified_radial_value_seam"]["passed"],
        physical_Qy_Z0_passed=report["actual_Rm_physical_Qy_Z0_value_seam"]["passed"],
        accepted_or_new_operator_replayed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "radial-values", "physical-qy0", "write-report", "inspect-report"):
        group.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        answer = dict(passed=True, AST_LF_EOF_only=True)
    elif args.radial_values:
        answer = run_component("radial")
    elif args.physical_qy0:
        answer = run_component("physical")
    else:
        answer = globals()[next(name for name, value in vars(args).items() if value)]()
    print(json.dumps(answer))
