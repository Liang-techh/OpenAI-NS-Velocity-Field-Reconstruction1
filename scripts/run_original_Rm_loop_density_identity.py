"""Checkpoint new actual common-phase, loop and signed-density seam operators."""
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

PARENT = PREFIX + "Rm_generic_native_identity.json"
PARENT_HASH = "71db4b9121e591878b4047025ed6d98e2a9a78d3659d82430ecfc294c0cafbf2"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
LOOP = PREFIX + "prescribed_loop_C5_functions.json.gz"
LIFT = PREFIX + "reference_phase_lift.json"
LIFT_HASH = "e1640aaa18b99fb376ebdf5e97b9cfa333b6eab567da24020d11d20b0bffdecf"
LOOP_SOURCE = "experiments/root_st073/" + PREFIX + "prescribed_loop_C5_functions.py"
LOOP_SOURCE_HASH = "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d"
PHASE_KERNEL = "src/openai_ns_reconstruction/original_Rm_common_phase_source.py"
KERNEL = "src/openai_ns_reconstruction/original_Rm_loop_density_identity.py"
RUNNER = "scripts/run_original_Rm_loop_density_identity.py"
REPORT = BASE / (PREFIX + "Rm_loop_density_identity.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_loop_density_identity"
SOURCE = CHECKPOINTS / "new_source.json"
PHASE_ENTRY = CHECKPOINTS / "new_phase_operator_entered.json"
PHASE_RESULT = CHECKPOINTS / "new_common_phase_source.json"
ENTRY = CHECKPOINTS / "new_loop_operator_entered.json"
RESULT = CHECKPOINTS / "new_loop_density_identity.json"


def module_at(relative, name):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def preflight():
    for relative in (RUNNER, PHASE_KERNEL, KERNEL):
        raw = (ROOT / relative).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + relative)
        ast.parse(raw, filename=relative)


def setup():
    preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: " + name)
    finite = read_pinned(BASE / FINITE, parent["input_hashes"][FINITE])
    loop = read_pinned(BASE / LOOP, parent["input_hashes"][LOOP])
    lift = read_pinned(BASE / LIFT, LIFT_HASH)
    require(sha(ROOT / LOOP_SOURCE) == LOOP_SOURCE_HASH, "loop theorem recipe changed")
    nodes = finite["exact_graph_nodes"]
    require(nodes[:121834] == loop["exact_graph_nodes"] and len(nodes) == 205147, "saved source prefix changed")
    first = parent["actual_Rm_native_radial_endpoint_source"]
    second = parent["actual_Rm_second_radial_endpoint_functions"]
    history = parent["actual_Rm_history_native_source"]["graph_delta"]
    require(digest(nodes) == first["base_digest"], "finite base changed")
    nodes = nodes + copy.deepcopy(first["appended_nodes"])
    require(digest(nodes) == first["full_nodes_digest"] == second["base_digest"], "first delta changed")
    nodes = nodes + copy.deepcopy(second["appended_nodes"])
    require(digest(nodes) == second["full_nodes_digest"] == history["base_nodes_sha256"], "second delta changed")
    nodes = nodes + copy.deepcopy(history["appended_nodes"])
    require(digest(nodes) == history["full_canonical_nodes_sha256"] and len(nodes) == 205200, "history delta changed")
    left = loop["prescribed_loop_C5_bindings"]["charts"]["postrestore"]
    right = finite["Rm_patch_finite_N_C5_bindings"]["chart"]
    generic = parent["actual_Rm_generic_native_identity"]
    require(generic["passed"] is True and not generic["unresolved_gates"], "leading generic source not admitted")
    packet = dict(schema="actual_Rm_loop_density_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, graph_node_count=len(nodes), graph_digest=digest(nodes),
        inherited_generic_digest=digest(generic), inherited_history_source_digest=digest(parent["actual_Rm_history_native_source"]),
        left_chart_digest=digest(left), right_chart_digest=digest(right), ordinary_Z_orders=[0, 1],
        loop_archive_sha256=parent["input_hashes"][LOOP], finite_archive_sha256=parent["input_hashes"][FINITE],
        saved_periodic_theorem_sha256=LIFT_HASH, same_loop_recipe_sha256=LOOP_SOURCE_HASH,
        runner_sha256=sha(ROOT / RUNNER),
        kernel_sha256=sha(ROOT / KERNEL), phase_kernel_sha256=sha(ROOT / PHASE_KERNEL),
        all_212_terms_and_errors_unchanged=True, no_accepted_operator_replayed=True)
    return parent, nodes, left, right, lift, packet


def prepare_source():
    *_, packet = setup()
    require(not any(p.exists() for p in (SOURCE, PHASE_ENTRY, PHASE_RESULT, ENTRY, RESULT)), "source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE), graph_node_count=packet["graph_node_count"])


def common_phase():
    parent, nodes, left, right, lift, packet = setup()
    require(read_saved(SOURCE) == packet and not PHASE_ENTRY.exists() and not PHASE_RESULT.exists(), "phase source changed or entered")
    module = module_at(PHASE_KERNEL, "new_Rm_common_phase_source")
    write_new(PHASE_ENTRY, dict(source_sha256=sha(SOURCE), kernel_sha256=packet["phase_kernel_sha256"]))
    started = time.perf_counter()
    try:
        result = module.prove(nodes, left, right, parent["actual_Rm_generic_native_inputs"]["radius_source_lemma"],
                              lift, parent["actual_Rm_leading_seed_admission"], parent["actual_Rm_generic_native_identity"],
                              loop_source_sha256=LOOP_SOURCE_HASH)
    except Exception as error:
        result = dict(passed=False, error_type=type(error).__name__, error=str(error),
                      failed_attempt_preserved=True, no_retry_authorized_by_runner=True)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
                  entry_checkpoint_sha256=sha(PHASE_ENTRY), actual_saved_graph_run=True)
    write_new(PHASE_RESULT, result)
    return dict(passed=result["passed"], sha256=sha(PHASE_RESULT), execution_seconds=result["execution_seconds"])


def loop_density():
    parent, nodes, left, right, lift, packet = setup()
    phase = read_saved(PHASE_RESULT)
    require(read_saved(SOURCE) == packet and phase["passed"] is True and not ENTRY.exists() and not RESULT.exists(),
            "loop source changed or entered")
    require(phase["same_loop_source_sha256"] == LOOP_SOURCE_HASH
            and phase["source_checkpoint_sha256"] == sha(SOURCE)
            and phase["entry_checkpoint_sha256"] == sha(PHASE_ENTRY)
            and phase["source_graph_digest"] == packet["graph_digest"], "phase proof link changed")
    module = module_at(KERNEL, "new_Rm_loop_density_identity")
    write_new(ENTRY, dict(source_sha256=sha(SOURCE), kernel_sha256=packet["kernel_sha256"], phase_sha256=sha(PHASE_RESULT)))
    started = time.perf_counter()
    try:
        result = module.derive(nodes, left, right, parent["actual_Rm_generic_native_identity"], phase,
                               expected_nodes_digest=packet["graph_digest"])
    except Exception as error:
        result = dict(passed=False, error_type=type(error).__name__, error=str(error),
                      failed_attempt_preserved=True, no_retry_authorized_by_runner=True)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
                  entry_checkpoint_sha256=sha(ENTRY), actual_saved_graph_run=True)
    write_new(RESULT, result)
    return dict(passed=result.get("passed"), sha256=sha(RESULT), execution_seconds=result["execution_seconds"],
                unresolved=result.get("unresolved_gates"))


def write_report():
    parent, nodes, left, right, lift, packet = setup()
    phase, result = read_saved(PHASE_RESULT), read_saved(RESULT)
    require(read_saved(SOURCE) == packet, "source checkpoint changed")
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_loop_density_identity", parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        actual_Rm_loop_density_inputs=packet, actual_Rm_common_phase_source=phase, actual_Rm_loop_density_identity=result,
        loop_density_checkpoint_hashes={p.name:sha(p) for p in (SOURCE, PHASE_ENTRY, PHASE_RESULT, ENTRY, RESULT)})
    report["input_hashes"].update({PARENT:PARENT_HASH, LIFT:LIFT_HASH})
    report["implementation_hashes"].update({p:sha(ROOT / p) for p in (RUNNER, PHASE_KERNEL, KERNEL)})
    report["execution_seconds"].update(new_actual_common_phase=phase["execution_seconds"], new_actual_loop_density=result["execution_seconds"])
    report["scope"].update(actual_current_Rm_common_phase_proved=phase["passed"],
        actual_current_Rm_signed_density_seam_proved=result["passed"], remaining_density_seam_unproved=not result["passed"],
        modified_finite_N_germ=False, current_Rm_RQy_germ_identity_proved=False,
        current_Rm_to_Rp_join=False, full_goal_complete=False)
    report["current_turn_loop_execution_audit"] = dict(actual_new_components_executed_once=[PHASE_KERNEL, KERNEL],
        accepted_ancestor_operators_replayed=False, saved_before_inspection=True,
        workers_fixed_role="luna_max_worker", worker_fork_turns="none", runtime_worker_model_metadata=None,
        numeric_field_or_quadrature_or_norm_computed=False, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=result["passed"], bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, left, right, lift, packet = setup()
    report, phase, result = read_saved(REPORT), read_saved(PHASE_RESULT), read_saved(RESULT)
    require(report["actual_Rm_loop_density_inputs"] == read_saved(SOURCE) == packet
            and report["actual_Rm_common_phase_source"] == phase
            and report["actual_Rm_loop_density_identity"] == result, "new source packet changed")
    require(report["loop_density_checkpoint_hashes"] == {p.name:sha(p) for p in (SOURCE, PHASE_ENTRY, PHASE_RESULT, ENTRY, RESULT)},
            "checkpoint link changed")
    require(phase["same_loop_source_sha256"] == LOOP_SOURCE_HASH
            and phase["source_checkpoint_sha256"] == sha(SOURCE)
            and phase["entry_checkpoint_sha256"] == sha(PHASE_ENTRY)
            and phase["source_graph_digest"] == packet["graph_digest"]
            and result["source_checkpoint_sha256"] == sha(SOURCE)
            and result["entry_checkpoint_sha256"] == sha(ENTRY), "execution proof link changed")
    for key in ("actual_Rm_generic_native_identity", "actual_Rm_history_native_source", "actual_Rm_absolute_m15",
                "corrected_own_swirl_provider", "corrected_reserved_point_numeric", "current_Rm_pressure_inlet_identity"):
        require(report[key] == parent[key], "parent packet changed: " + key)
    for name, expected in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "report pin changed: " + name)
    require(phase["passed"] is True and report["scope"]["full_goal_complete"] is False, "phase/goal status changed")
    return dict(passed=True, source_identity_passed=result["passed"], serialized_only=True,
        bytes=REPORT.stat().st_size, sha256=sha(REPORT), new_operator_replayed=False,
        old_mathematics_executed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "common-phase", "loop-density", "write-report", "inspect-report"):
        group.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        answer = dict(passed=True, AST_LF_EOF_only=True)
    else:
        answer = globals()[next(name for name, value in vars(args).items() if value)]()
    print(json.dumps(answer))
