"""Checkpoint new actual Duhamel and modified-moment Rm seam consumers."""
from __future__ import annotations

import argparse
import ast
import copy
import importlib.util
import json
import sys
import time

from run_original_Rm_second_radial_identity import ROOT, BASE, PREFIX, require, sha, digest, path_for, read_pinned, read_saved, write_new

PARENT = PREFIX + "Rm_loop_density_identity.json.gz"
PARENT_HASH = "29c4fec78d3c622b36bb10867792b3d4df1b62f332cf3bc93cd5ed684706fe29"
BASE_REPORT = PREFIX + "Rm_generic_native_identity.json"
BASE_HASH = "71db4b9121e591878b4047025ed6d98e2a9a78d3659d82430ecfc294c0cafbf2"
PORTS = "experiments/root_st073/original_Rm_loop_ports_admission.json"
PORTS_HASH = "0331c5fc37bd3ec6c0d386c38b07d9863ede622e14f1c6dc65bf8cdf675a16b9"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
LOOP = PREFIX + "prescribed_loop_C5_functions.json.gz"
EXTRACTOR = "scripts/extract_original_Rm_loop_ports.py"
HISTORY_KERNEL = "src/openai_ns_reconstruction/original_Rm_history_seam_identity.py"
MOMENT_KERNEL = "src/openai_ns_reconstruction/original_Rm_modified_moment_germ.py"
RUNNER = "scripts/run_original_Rm_history_moment_seam.py"
REPORT = BASE / (PREFIX + "Rm_history_moment_seam.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_history_moment_seam"
SOURCE = CHECKPOINTS / "new_source.json"
HISTORY_ENTRY = CHECKPOINTS / "new_history_operator_entered.json"
HISTORY_RESULT = CHECKPOINTS / "new_history_seam_identity.json"
MOMENT_ENTRY = CHECKPOINTS / "new_moment_operator_entered.json"
MOMENT_RESULT = CHECKPOINTS / "new_modified_moment_germ.json"


def load_module(relative, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def preflight():
    for name in (RUNNER, HISTORY_KERNEL, MOMENT_KERNEL, EXTRACTOR):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + name)
        ast.parse(raw, filename=name)


def setup():
    preflight()
    ports = read_pinned(ROOT / PORTS, PORTS_HASH)
    require(ports["source_archive_sha256"] == PARENT_HASH and sha(BASE / PARENT) == PARENT_HASH
            and ports["extraction_implementation_sha256"] == sha(ROOT / EXTRACTOR), "bounded saved admission changed")
    for name, expected in {**ports["input_hashes"], **ports["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "inherited source pin changed: " + name)
    base = read_pinned(BASE / BASE_REPORT, BASE_HASH)
    for key, expected in ports["inherited_packet_digests"].items():
        require(digest(base[key]) == expected, "bounded projection changed inherited packet: " + key)
    finite = read_pinned(BASE / FINITE, ports["input_hashes"][FINITE])
    loop = read_pinned(BASE / LOOP, ports["input_hashes"][LOOP])
    nodes = finite["exact_graph_nodes"]
    require(len(nodes) == 205147 and nodes[:121834] == loop["exact_graph_nodes"], "saved exact prefix changed")
    first, second = base["actual_Rm_native_radial_endpoint_source"], base["actual_Rm_second_radial_endpoint_functions"]
    history = base["actual_Rm_history_native_source"]["graph_delta"]
    require(digest(nodes) == first["base_digest"], "finite graph digest changed")
    nodes = nodes + copy.deepcopy(first["appended_nodes"])
    require(digest(nodes) == first["full_nodes_digest"] == second["base_digest"], "first endpoint delta changed")
    nodes = nodes + copy.deepcopy(second["appended_nodes"])
    require(digest(nodes) == second["full_nodes_digest"] == history["base_nodes_sha256"], "second endpoint delta changed")
    nodes = nodes + copy.deepcopy(history["appended_nodes"])
    require(len(nodes) == 205200 and digest(nodes) == history["full_canonical_nodes_sha256"] == ports["nodes_digest"], "source205200 changed")
    lb = loop["prescribed_loop_C5_bindings"]
    left = copy.deepcopy(lb["charts"]["postrestore"])
    left["outlet"] = dict(Rm_correction=copy.deepcopy(lb["Rm_correction"]))
    right = finite["Rm_patch_finite_N_C5_bindings"]["chart"]
    packet = dict(schema="actual_Rm_history_moment_seam_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, inherited_payload_report=BASE_REPORT,
        inherited_payload_sha256=BASE_HASH, bounded_loop_ports=PORTS, bounded_loop_ports_sha256=PORTS_HASH,
        graph_node_count=len(nodes), graph_digest=digest(nodes), left_chart_digest=digest(left), right_chart_digest=digest(right),
        patch_history_source_digest=digest(base["actual_Rm_history_native_source"]),
        seed_admission_digest=digest(base["actual_Rm_leading_seed_admission"]), raw_native_digest=digest(base["actual_Rm_raw_native_jet_identity"]),
        ordinary_Z_orders=[0, 1], kernel_sha256={name:sha(ROOT / name) for name in (RUNNER, HISTORY_KERNEL, MOMENT_KERNEL, EXTRACTOR)},
        verbose_old_result_not_copied=True, accepted_mathematical_operator_replayed=False,
        all_old_212_terms_and_errors_preserved=True, full_goal_complete=False)
    return base, ports, nodes, left, right, packet


def prepare_source():
    *_, packet = setup()
    require(not any(p.exists() for p in (SOURCE, HISTORY_ENTRY, HISTORY_RESULT, MOMENT_ENTRY, MOMENT_RESULT)), "source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE), graph_node_count=packet["graph_node_count"])


def history_seam():
    base, ports, nodes, left, right, packet = setup()
    require(read_saved(SOURCE) == packet and not HISTORY_ENTRY.exists() and not HISTORY_RESULT.exists(), "history source changed or already entered")
    module = load_module(HISTORY_KERNEL, "new_actual_Rm_history_seam")
    write_new(HISTORY_ENTRY, dict(source_sha256=sha(SOURCE), kernel_sha256=packet["kernel_sha256"][HISTORY_KERNEL]))
    started = time.perf_counter()
    try:
        result = module.derive(nodes, left, right, base["actual_Rm_history_native_source"], ports,
                               expected_nodes_digest=packet["graph_digest"])
    except Exception as error:
        result = dict(passed=False, error_type=type(error).__name__, error=str(error), failed_attempt_preserved=True)
    result.update(execution_seconds=time.perf_counter()-started, actual_source_graph_inspected=True,
                  source_checkpoint_sha256=sha(SOURCE), entry_checkpoint_sha256=sha(HISTORY_ENTRY))
    write_new(HISTORY_RESULT, result)
    return dict(passed=result["passed"], sha256=sha(HISTORY_RESULT), seconds=result["execution_seconds"], error=result.get("error"))


def moment_germ():
    base, ports, nodes, left, right, packet = setup()
    history = read_saved(HISTORY_RESULT)
    require(read_saved(SOURCE) == packet and history["passed"] is True
            and not MOMENT_ENTRY.exists() and not MOMENT_RESULT.exists(), "moment source changed or already entered")
    delta = history["graph_delta"]
    require(delta["base_nodes_sha256"] == packet["graph_digest"], "history source base changed")
    nodes = nodes + copy.deepcopy(delta["appended_nodes"])
    require(digest(nodes) == delta["full_nodes_sha256"], "history source definitions changed")
    module = load_module(MOMENT_KERNEL, "new_actual_Rm_modified_moment_germ")
    write_new(MOMENT_ENTRY, dict(source_sha256=sha(SOURCE), history_checkpoint_sha256=sha(HISTORY_RESULT),
                               kernel_sha256=packet["kernel_sha256"][MOMENT_KERNEL]))
    started = time.perf_counter()
    try:
        result = module.derive(nodes, left, right, base["actual_Rm_leading_seed_admission"],
            base["actual_Rm_leading_seed_inputs"], base["actual_Rm_raw_native_jet_identity"], history,
            expected_nodes_digest=digest(nodes))
    except Exception as error:
        result = dict(passed=False, error_type=type(error).__name__, error=str(error), failed_attempt_preserved=True)
    result.update(execution_seconds=time.perf_counter()-started, actual_source_graph_inspected=True,
                  source_checkpoint_sha256=sha(SOURCE), entry_checkpoint_sha256=sha(MOMENT_ENTRY))
    write_new(MOMENT_RESULT, result)
    return dict(passed=result["passed"], sha256=sha(MOMENT_RESULT), seconds=result["execution_seconds"], error=result.get("error"))


def write_report():
    base, ports, nodes, left, right, packet = setup()
    history, moment = read_saved(HISTORY_RESULT), read_saved(MOMENT_RESULT)
    require(read_saved(SOURCE) == packet and history["passed"] is True and moment["passed"] is True, "new source seam not admitted")
    report = copy.deepcopy(base)
    bounded = {key:value for key, value in ports.items() if key not in ("input_hashes", "implementation_hashes")}
    report.update(stage="actual_Rm_history_modified_moment_seam", parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        inherited_payload_report=BASE_REPORT, inherited_payload_sha256=BASE_HASH,
        bounded_Rm_loop_source_admission=bounded, actual_Rm_history_moment_inputs=packet,
        actual_Rm_history_seam_identity=history, actual_Rm_modified_moment_germ=moment,
        history_moment_checkpoint_hashes={p.name:sha(p) for p in (SOURCE, HISTORY_ENTRY, HISTORY_RESULT, MOMENT_ENTRY, MOMENT_RESULT)})
    report["input_hashes"] = dict(ports["input_hashes"])
    report["input_hashes"].update({PARENT:PARENT_HASH, PORTS:PORTS_HASH})
    report["implementation_hashes"] = dict(ports["implementation_hashes"])
    report["implementation_hashes"].update(packet["kernel_sha256"])
    report["execution_seconds"].update(new_actual_history_seam=history["execution_seconds"], new_actual_modified_moment_germ=moment["execution_seconds"])
    report["scope"].update(actual_current_Rm_common_phase_proved=True, actual_current_Rm_signed_density_seam_proved=True,
        remaining_density_seam_unproved=False, modified_five_moment_value_native_seam_Z01_proved=True,
        slow_native_loop_seam_proved=False, modified_E_V_native_seam_proved=False,
        modified_finite_N_germ=False, current_Rm_RQy_germ_identity_proved=False,
        current_Rm_to_Rp_join=False, full_goal_complete=False)
    report["current_turn_history_moment_execution_audit"] = dict(new_actual_components_once=[HISTORY_KERNEL, MOMENT_KERNEL],
        accepted_math_operators_replayed=False, saved_before_inspection=True, bounded_saved_projection_used=True,
        runtime_worker_model_metadata=None, workers_fixed_role="luna_max_worker", worker_fork_turns="none",
        numerical_field_or_quadrature_or_norm_computed=False, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    base, ports, nodes, left, right, packet = setup()
    report, history, moment = read_saved(REPORT), read_saved(HISTORY_RESULT), read_saved(MOMENT_RESULT)
    require(report["actual_Rm_history_moment_inputs"] == read_saved(SOURCE) == packet
            and report["actual_Rm_history_seam_identity"] == history
            and report["actual_Rm_modified_moment_germ"] == moment, "saved new packet changed")
    require(report["history_moment_checkpoint_hashes"] == {p.name:sha(p) for p in (SOURCE, HISTORY_ENTRY, HISTORY_RESULT, MOMENT_ENTRY, MOMENT_RESULT)}, "checkpoint bytes changed")
    require(len(history["rows"]) == len(moment["rows"]) == 10 and history["passed"] is True and moment["passed"] is True,
            "new source proof rows incomplete")
    for result, entry in ((history, HISTORY_ENTRY), (moment, MOMENT_ENTRY)):
        require(result["source_checkpoint_sha256"] == sha(SOURCE) and result["entry_checkpoint_sha256"] == sha(entry), "execution checkpoint link changed")
        delta = result["graph_delta"]
        require(digest(nodes) == delta["base_nodes_sha256"], "graph delta base changed")
        nodes = nodes + copy.deepcopy(delta["appended_nodes"])
        require(digest(nodes) == delta["full_nodes_sha256"], "graph delta definition changed")
    for key in ("actual_Rm_generic_native_identity", "actual_Rm_history_native_source", "actual_Rm_absolute_m15",
                "corrected_own_swirl_provider", "corrected_reserved_point_numeric", "current_Rm_pressure_inlet_identity"):
        require(report[key] == base[key], "inherited packet changed: " + key)
    for name, expected in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "report source pin changed: " + name)
    require(report["scope"]["full_goal_complete"] is False and report["scope"]["current_Rm_RQy_germ_identity_proved"] is False, "open construction gate changed")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
        new_history_pairs=len(history["rows"]), new_modified_moment_pairs=len(moment["rows"]), graph_node_count=len(nodes),
        accepted_or_new_operator_replayed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "history-seam", "moment-germ", "write-report", "inspect-report"):
        group.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        answer = dict(passed=True, AST_LF_EOF_only=True)
    else:
        answer = globals()[next(name for name, value in vars(args).items() if value)]()
    print(json.dumps(answer))
