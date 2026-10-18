"""Run the new control-bound patch-only native provenance operator once."""
from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import importlib.util
import json
import sys
import time

from run_original_Rm_slow_native_provenance import ROOT, BASE, PREFIX, require, sha, digest, read_pinned, read_saved, path_for, write_new, append_delta

PARENT = PREFIX + "Rm_slow_native_provenance.json"
PARENT_HASH = "7ba0883d9a6c0af55f7ab3c566da85f1287d30e70e25a493827e159f36a6caab"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
FINITE_HASH = "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"
LEADING = PREFIX + "true_Rm_patch_C6_functions.json.gz"
LEADING_HASH = "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a"
BRANCH = PREFIX + "whole_Z_Rm_centered_inverse_check.json"
BRANCH_HASH = "29781f93710407553ca958a192b326f7c9bae6ebf4f269794116c452837b0702"
RUNNER = "scripts/run_original_Rm_patch_native_provenance.py"
KERNEL = "src/openai_ns_reconstruction/original_Rm_patch_native_provenance.py"
REPORT = BASE / (PREFIX + "Rm_patch_native_provenance.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_patch_native_provenance_v2"
SOURCE = CHECKPOINTS / "new_source.json"
ENTRY = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_patch_native_provenance.json"
GAMMA = ((185701, 185702, 185703), (185728, 185729, 185730), (185755, 185756, 185757))
FAILED_EVIDENCE = "experiments/root_st073/original_Rm_patch_native_provenance_attempt1.json"


def preflight():
    for name in (RUNNER, KERNEL):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + name)
        ast.parse(raw, filename=name)


def setup():
    preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "inherited pin changed: " + name)
    finite = read_pinned(BASE / FINITE, FINITE_HASH)
    leading = read_pinned(BASE / LEADING, LEADING_HASH)
    branch = read_pinned(BASE / BRANCH, BRANCH_HASH)
    family = finite["source_family"]
    require(family == leading["source_family"] == branch["source_family"] and branch["all_passed"] is True
            and branch["same_common_unit_box_proves_unique_whole_Z_smooth_leading_family"] is True
            and branch["residual_zero_containment_not_used_as_functional_closure_proof"] is True,
            "original unique Z-family premise changed")
    nodes = finite["exact_graph_nodes"]
    require(len(nodes) == 205147 and nodes[:196042] == leading["exact_graph_nodes"], "actual source prefix changed")
    prefix_digest = digest(nodes)
    nodes = append_delta(nodes, parent["actual_Rm_native_radial_endpoint_source"], "base_digest", "full_nodes_digest")
    nodes = append_delta(nodes, parent["actual_Rm_second_radial_endpoint_functions"], "base_digest", "full_nodes_digest")
    nodes = append_delta(nodes, parent["actual_Rm_history_native_source"]["graph_delta"], "base_nodes_sha256", "full_canonical_nodes_sha256")
    nodes = append_delta(nodes, parent["actual_Rm_history_seam_identity"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256")
    nodes = append_delta(nodes, parent["actual_Rm_modified_moment_germ"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256")
    require(len(nodes) == 205300 and digest(nodes) == parent["actual_Rm_slow_native_provenance_inputs"]["graph_digest"], "current source changed")
    bindings = leading["true_Rm_patch_C6_bindings"]
    f = finite["Rm_patch_finite_N_C5_bindings"]
    loop = copy.deepcopy(f["chart"])
    loop.update(radial=copy.deepcopy(f["radial"]), invJ=198869)
    require(f["radial_background"] == parent["actual_Rm_raw_native_jet_inputs"]["patch_radial_background"], "accepted raw seeds changed")
    left = parent["actual_Rm_slow_native_provenance"]["sides"]["postrestore"]
    require(left["passed"] is True and len(left["rows"]) == 4, "left admitted source data missing")
    prior_patch = parent["actual_Rm_slow_native_provenance"]["sides"]["patch"]
    require(prior_patch["passed"] is False and prior_patch["visited_source_D_pairs"] == [], "prior patch unexpectedly admitted")
    require(all(p["source_handle"] < 183700 and p["derivative_handle"] < 183700
                for p in left["visited_source_D_pairs"]), "left provenance unexpectedly depends on later patch gamma")
    contracts = {}
    for name, class_name, function_name in (
        ("experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_C6_functions.py", "ExactBumpWeights", "gamma"),
        ("experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py", None, "native_derivative"),
    ):
        tree = ast.parse((ROOT / name).read_bytes())
        owner = tree if class_name is None else next(q for q in tree.body if isinstance(q, ast.ClassDef) and q.name == class_name)
        function = next(q for q in owner.body if isinstance(q, ast.FunctionDef) and q.name == function_name)
        contracts[name] = dict(sha256=sha(ROOT / name), function=function_name, class_name=class_name,
            AST_sha256=hashlib.sha256(ast.dump(function, include_attributes=False).encode()).hexdigest())
    inputs = dict(leading_bindings=bindings, loop_chart=loop, radial_background=f["radial_background"],
        accepted_left_provenance=left, expected_nodes_digest=digest(nodes), derivative_prefix_digest=prefix_digest,
        source_family=family, gamma_seed_rows=GAMMA, branch_receipt_sha256=BRANCH_HASH)
    packet = dict(schema="Rm_control_bound_patch_native_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, nodes_digest=digest(nodes), node_count=len(nodes),
        derivative_source_prefix_count=205147, derivative_source_prefix_digest=prefix_digest,
        input_digest=digest(inputs), accepted_left_provenance_digest=digest(left),
        archive_hashes={FINITE:FINITE_HASH, LEADING:LEADING_HASH, BRANCH:BRANCH_HASH},
        implementation_hashes={name:sha(ROOT / name) for name in (RUNNER, KERNEL)},
        original_normalized_gamma_contracts=contracts,
        failed_attempt_evidence=dict(path=FAILED_EVIDENCE, sha256=sha(ROOT / FAILED_EVIDENCE)),
        new_scope="actual control native-independence contract and patch-only provenance",
        accepted_left_traversal_replayed=False, ordinary_Z_orders=[0, 1], control_support_rows=35,
        gamma_seed_scope="actual normalized gamma products with 40/1600/64000 scales; supersedes untraversed prior patch raw-beta seeds",
        old_accepted_mathematical_operator_replayed=False, full_goal_complete=False)
    return parent, nodes, inputs, packet


def prepare_source():
    *_, packet = setup()
    require(not any(p.exists() for p in (SOURCE, ENTRY, RESULT)), "new source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, sha256=sha(SOURCE), node_count=packet["node_count"])


def patch_provenance():
    parent, nodes, inputs, packet = setup()
    require(read_saved(SOURCE) == packet and not ENTRY.exists() and not RESULT.exists(), "input changed or operator already entered")
    sys.path.insert(0, str(ROOT / "src"))
    spec = importlib.util.spec_from_file_location("new_Rm_patch_native_provenance", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTRY, dict(source_sha256=sha(SOURCE), kernel_sha256=packet["implementation_hashes"][KERNEL]))
    started = time.perf_counter()
    try:
        result = module.derive(nodes, **inputs)
    except Exception as error:
        result = dict(passed=False, error_type=type(error).__name__, error=str(error), failed_attempt_preserved=True)
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
        entry_checkpoint_sha256=sha(ENTRY), actual_source_graph_inspected=True, actual_saved_graph_run=False)
    write_new(RESULT, result)
    return dict(passed=result["passed"], sha256=sha(RESULT), seconds=result["execution_seconds"],
        error=result.get("error"), unresolved_gates=result.get("unresolved_gates", []))


def write_report():
    parent, nodes, inputs, packet = setup()
    result = read_saved(RESULT)
    require(read_saved(SOURCE) == packet, "source checkpoint changed")
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_control_bound_patch_native_provenance", parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        actual_Rm_patch_native_provenance_inputs=packet, actual_Rm_patch_native_provenance=result,
        patch_native_checkpoint_hashes={p.name:sha(p) for p in (SOURCE, ENTRY, RESULT)})
    report["input_hashes"][PARENT] = PARENT_HASH
    report["input_hashes"][FAILED_EVIDENCE] = packet["failed_attempt_evidence"]["sha256"]
    report["input_hashes"].update(packet["archive_hashes"])
    report["implementation_hashes"].update(packet["implementation_hashes"])
    report["execution_seconds"]["new_control_bound_patch_native_provenance"] = result["execution_seconds"]
    report["scope"].update(postrestore_slow_native_provenance_proved=True,
        actual_slow_native_provenance_checked=result["passed"],
        patch_slow_native_provenance_proved=result["passed"], full_slow_native_provenance_passed=result["passed"],
        slow_native_loop_seam_proved=False, modified_E_V_native_seam_proved=False,
        current_Rm_RQy_germ_identity_proved=False, full_goal_complete=False)
    report["current_turn_patch_native_execution_audit"] = dict(new_patch_operator_once=True,
        accepted_left_result_consumed_as_data=True, accepted_left_traversal_replayed=False,
        accepted_ancestors_or_fixtures_or_providers_replayed=False, saved_before_inspection=True,
        worker_fixed_role="luna_max_worker", worker_fork_turns="none", runtime_worker_model_metadata=None,
        numerical_field_or_quadrature_or_norm_computed=False)
    write_new(REPORT, report)
    return dict(passed=True, mathematical_patch_provenance_passed=result["passed"], bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, inputs, packet = setup()
    report, result = read_saved(REPORT), read_saved(RESULT)
    require(report["actual_Rm_patch_native_provenance_inputs"] == read_saved(SOURCE) == packet
            and report["actual_Rm_patch_native_provenance"] == result, "saved new packet changed")
    require(report["patch_native_checkpoint_hashes"] == {p.name:sha(p) for p in (SOURCE, ENTRY, RESULT)}, "checkpoint bytes changed")
    require(result["source_checkpoint_sha256"] == sha(SOURCE) and result["entry_checkpoint_sha256"] == sha(ENTRY), "entry links changed")
    for key in ("actual_Rm_slow_native_provenance", "actual_Rm_generic_native_identity", "actual_Rm_history_seam_identity",
                "actual_Rm_modified_moment_germ", "actual_Rm_absolute_m15", "corrected_own_swirl_provider",
                "corrected_reserved_point_numeric", "current_Rm_pressure_inlet_identity"):
        require(report[key] == parent[key], "inherited data changed: " + key)
    for name, expected in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "saved source pin changed: " + name)
    require(not report["scope"]["slow_native_loop_seam_proved"] and not report["scope"]["current_Rm_RQy_germ_identity_proved"]
            and not report["scope"]["full_goal_complete"], "open physical gate changed")
    return dict(passed=True, serialized_only=True, mathematical_patch_provenance_passed=result["passed"],
        patch_rows=result.get("rows", []), unresolved_gates=result.get("unresolved_gates", []),
        bytes=REPORT.stat().st_size, sha256=sha(REPORT), accepted_or_new_operator_replayed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "patch-provenance", "write-report", "inspect-report"):
        group.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        answer = dict(passed=True, AST_LF_EOF_only=True)
    else:
        answer = globals()[next(name for name, value in vars(args).items() if value)]()
    print(json.dumps(answer))
