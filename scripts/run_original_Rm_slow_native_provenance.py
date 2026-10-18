"""Checkpoint a new read-only binding of saved slow native derivative roots."""
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
PARENT = PREFIX + "Rm_history_moment_seam.json"
PARENT_HASH = "22734d97f79fd16cfa9e07e05df824b550e4049d040aac9168791e60f0760d00"
ARCHIVES = {
    "loop": ("prescribed_loop_C5", "8419fdc338dc746099c7cad119ea61a8e928652319d5caf5314f8a730c533caa"),
    "modified": ("modified_Q4_pressure_C3", "b96316728406fddb9d94e866042529b7d84b7146167cba221ed5a88b427f147a"),
    "leading": ("true_Rm_patch_C6", "4352ebd8862b5c92e3e4026728183014063baf556802a88421cdbd517ff1346a"),
    "finite": ("true_Rm_patch_finite_N_C5", "8f3bc343e9dbffbdda465ce94d6183735184645a70e291cb900d3cde726769ad"),
}
RUNNER = "scripts/run_original_Rm_slow_native_provenance.py"
KERNEL = "src/openai_ns_reconstruction/original_Rm_native_derivative_provenance.py"
REPORT = BASE / (PREFIX + "Rm_slow_native_provenance.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_slow_native_provenance_v2"
SOURCE = CHECKPOINTS / "new_source.json"
ENTRY = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_native_derivative_provenance.json"
GAMMA = ((185679, 185699, 185700), (185706, 185726, 185727), (185733, 185753, 185754))
FAILED_EVIDENCE = "experiments/root_st073/original_Rm_slow_native_provenance_attempt1.json"


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def read_pinned(path, expected):
    raw = Path(path).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, "changed pinned input: " + str(path))
    return json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)


def read_saved(path):
    return json.loads(path.read_bytes())


def path_for(name):
    path = Path(name)
    require(not path.is_absolute(), "absolute input path")
    if ".." in path.parts:
        require(name == "../../src/openai_ns_reconstruction/schedule_pressure.py", "unexpected traversal")
        return (BASE / path).resolve()
    return ROOT / path if path.parts[0] in ("experiments", "src", "scripts", "tmp") else BASE / path.name


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write((json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode())


def preflight():
    for name in (RUNNER, KERNEL):
        raw = (ROOT / name).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + name)
        ast.parse(raw, filename=name)


def append_delta(nodes, delta, base_key, full_key):
    require(digest(nodes) == delta[base_key], "saved delta base changed")
    nodes = nodes + copy.deepcopy(delta["appended_nodes"])
    require(digest(nodes) == delta[full_key], "saved delta content changed")
    return nodes


def setup():
    preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "inherited pin changed: " + name)
    archive_paths = {key:PREFIX + stem + "_functions.json.gz" for key, (stem, _) in ARCHIVES.items()}
    archives = {key:read_pinned(BASE / archive_paths[key], expected) for key, (_, expected) in ARCHIVES.items()}
    nodes = archives["finite"]["exact_graph_nodes"]
    require(len(nodes) == 205147, "finite source count changed")
    derivative_prefix_digest = digest(nodes)
    for key in ("loop", "modified", "leading"):
        prefix = archives[key]["exact_graph_nodes"]
        require(nodes[:len(prefix)] == prefix, "actual saved prefix changed: " + key)
    nodes = append_delta(nodes, parent["actual_Rm_native_radial_endpoint_source"], "base_digest", "full_nodes_digest")
    nodes = append_delta(nodes, parent["actual_Rm_second_radial_endpoint_functions"], "base_digest", "full_nodes_digest")
    nodes = append_delta(nodes, parent["actual_Rm_history_native_source"]["graph_delta"], "base_nodes_sha256", "full_canonical_nodes_sha256")
    nodes = append_delta(nodes, parent["actual_Rm_history_seam_identity"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256")
    nodes = append_delta(nodes, parent["actual_Rm_modified_moment_germ"]["graph_delta"], "base_nodes_sha256", "full_nodes_sha256")
    require(len(nodes) == 205300, "latest source count changed")
    left = copy.deepcopy(archives["loop"]["prescribed_loop_C5_bindings"]["charts"]["postrestore"])
    q4 = archives["modified"]["modified_Q4_pressure_C3_bindings"]["charts"]["postrestore"]
    left["radial"] = q4["radial"]
    left["invJ"] = 172066
    left_source = copy.deepcopy(parent["actual_Rm_leading_seed_inputs"]["postrestore"])
    left_source["physical_measure"] = left["physical_measure"]
    left_source["source_family"] = archives["modified"]["source_family"]
    left_source["gamma_seed_rows"] = []
    patch = archives["finite"]["Rm_patch_finite_N_C5_bindings"]
    right = copy.deepcopy(patch["chart"])
    right["radial"] = patch["radial"]
    right["invJ"] = 198869
    right_source = copy.deepcopy(archives["leading"]["true_Rm_patch_C6_bindings"]["chart"])
    right_source["source_family"] = archives["leading"]["source_family"]
    right_source["gamma_seed_rows"] = GAMMA
    require(q4["radial_background"] == parent["actual_Rm_raw_native_jet_inputs"]["modified_radial_background"]
            and patch["radial_background"] == parent["actual_Rm_raw_native_jet_inputs"]["patch_radial_background"],
            "actual raw native source seed handles changed")
    # These are producer-contract seeds, not new analytic derivative proofs.
    # Authenticate their exact saved topology; the final certificate keeps
    # this premise separate from the independently checked chain-rule DAG.
    for raw, first, second in GAMMA:
        beta = nodes[raw]
        require(beta["operation"] == "compact_raw_beta", "gamma source changed")
        for ref in (first, second):
            row = nodes[ref]
            require(row["operation"] == "original_lazy_flat_branch" and row["Delta"] == beta["argument"]
                    and row["flat_value"] == 0 and row["original_beta_exact_flat_derivative"] is True, "gamma derivative contract changed")
    source_contract = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py"
    source_tree = ast.parse((ROOT / source_contract).read_bytes())
    contract_function = next(n for n in source_tree.body if isinstance(n, ast.FunctionDef) and n.name == "native_derivative")
    contract_ast = hashlib.sha256(ast.dump(contract_function, include_attributes=False).encode()).hexdigest()
    inputs = {
        "postrestore": dict(source_chart=left_source, loop_chart=left, radial_background=q4["radial_background"], gamma_seed_rows=[]),
        "patch": dict(source_chart=right_source, loop_chart=right, radial_background=patch["radial_background"], gamma_seed_rows=GAMMA),
    }
    packet = dict(schema="actual_Rm_slow_native_provenance_inputs_v2", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, graph_node_count=len(nodes), graph_digest=digest(nodes),
        chart_input_digests={side:digest(value) for side, value in inputs.items()},
        archive_hashes={archive_paths[key]:expected for key, (_, expected) in ARCHIVES.items()},
        implementation_hashes={name:sha(ROOT / name) for name in (RUNNER, KERNEL)},
        original_patch_native_derivative_contract=dict(path=source_contract, sha256=sha(ROOT / source_contract),
            function="native_derivative", AST_sha256=contract_ast, gamma_seed_rows=GAMMA,
            interpretation="saved producer derivative seed contract; not a new proof of beta derivative formulas"),
        accepted_raw_native_identity_digest=digest(parent["actual_Rm_raw_native_jet_identity"]),
        derivative_source_prefix_count=205147, derivative_source_prefix_digest=derivative_prefix_digest,
        failed_attempt_evidence=dict(path=FAILED_EVIDENCE, sha256=sha(ROOT / FAILED_EVIDENCE)),
        ordinary_Z_orders=[0, 1], raw_seed_ordinary_Z_orders=[0, 1, 2],
        mathematical_scope="existing actual slow target provenance only; no cross-chart derivative equality",
        accepted_mathematical_operator_replayed=False, full_goal_complete=False)
    return parent, nodes, inputs, packet


def prepare_source():
    *_, packet = setup()
    require(not any(p.exists() for p in (SOURCE, ENTRY, RESULT)), "source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, sha256=sha(SOURCE), graph_node_count=packet["graph_node_count"])


def provenance():
    parent, nodes, inputs, packet = setup()
    require(read_saved(SOURCE) == json.loads(json.dumps(packet)) and not ENTRY.exists() and not RESULT.exists(), "source changed or new operator already entered")
    spec = importlib.util.spec_from_file_location("new_Rm_native_derivative_provenance", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTRY, dict(source_sha256=sha(SOURCE), kernel_sha256=packet["implementation_hashes"][KERNEL]))
    started = time.perf_counter()
    sides = {}
    for side, supplied in inputs.items():
        try:
            sides[side] = module.verify(nodes, **supplied, expected_nodes_digest=packet["graph_digest"], side=side,
                derivative_source_prefix_count=packet["derivative_source_prefix_count"],
                expected_derivative_prefix_digest=packet["derivative_source_prefix_digest"])
        except Exception as error:
            sides[side] = dict(passed=False, error_type=type(error).__name__, error=str(error), failed_attempt_preserved=True)
    result = dict(schema="actual_Rm_slow_native_provenance_result_v1", passed=all(s.get("passed") is True for s in sides.values()),
        sides=sides, source_only=True, actual_source_graph_inspected=True, actual_saved_graph_run=False,
        derivative_seed_scope=packet["original_patch_native_derivative_contract"],
        cross_chart_native_seam_proved=False, numeric_field_evaluated=False, execution_seconds=time.perf_counter()-started,
        source_checkpoint_sha256=sha(SOURCE), entry_checkpoint_sha256=sha(ENTRY), full_goal_complete=False)
    write_new(RESULT, result)
    return dict(passed=result["passed"], sha256=sha(RESULT), seconds=result["execution_seconds"],
        sides={side:dict(passed=value.get("passed"), error=value.get("error"), unresolved_count=len(value.get("unresolved", []))) for side, value in sides.items()})


def write_report():
    parent, nodes, inputs, packet = setup()
    source, result = read_saved(SOURCE), read_saved(RESULT)
    require(source == json.loads(json.dumps(packet)), "source checkpoint changed")
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_slow_native_derivative_provenance", parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        actual_Rm_slow_native_provenance_inputs=source, actual_Rm_slow_native_provenance=result,
        slow_native_provenance_checkpoint_hashes={p.name:sha(p) for p in (SOURCE, ENTRY, RESULT)})
    report["input_hashes"][PARENT] = PARENT_HASH
    report["input_hashes"][FAILED_EVIDENCE] = packet["failed_attempt_evidence"]["sha256"]
    report["input_hashes"].update(packet["archive_hashes"])
    report["implementation_hashes"].update(packet["implementation_hashes"])
    report["execution_seconds"]["new_actual_slow_native_provenance"] = result["execution_seconds"]
    report["scope"].update(actual_slow_native_provenance_checked=result["passed"],
        slow_native_loop_seam_proved=False, modified_E_V_native_seam_proved=False,
        current_Rm_RQy_germ_identity_proved=False, full_goal_complete=False)
    report["current_turn_slow_native_execution_audit"] = dict(new_actual_component_once=KERNEL,
        accepted_mathematical_operators_replayed=False, saved_before_inspection=True,
        runtime_worker_model_metadata=None, worker_fixed_role="luna_max_worker", worker_fork_turns="none",
        worker_reported_new_kernel_toy_verify_calls=4, worker_reported_new_gamma_seed_toy_calls=3,
        worker_toy_checkpoint_available=False, worker_execution_deviated_from_read_only_brief=True,
        numerical_field_or_quadrature_or_norm_computed=False)
    write_new(REPORT, report)
    return dict(passed=True, mathematical_provenance_passed=result["passed"], bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, inputs, packet = setup()
    report, result = read_saved(REPORT), read_saved(RESULT)
    require(report["actual_Rm_slow_native_provenance_inputs"] == read_saved(SOURCE) == json.loads(json.dumps(packet))
            and report["actual_Rm_slow_native_provenance"] == result, "saved packet changed")
    require(report["slow_native_provenance_checkpoint_hashes"] == {p.name:sha(p) for p in (SOURCE, ENTRY, RESULT)}, "checkpoint bytes changed")
    for key in ("actual_Rm_generic_native_identity", "actual_Rm_history_seam_identity", "actual_Rm_modified_moment_germ",
                "actual_Rm_absolute_m15", "corrected_own_swirl_provider", "corrected_reserved_point_numeric", "current_Rm_pressure_inlet_identity"):
        require(report[key] == parent[key], "inherited packet changed: " + key)
    for name, expected in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "saved pin changed: " + name)
    require(not report["scope"]["slow_native_loop_seam_proved"] and not report["scope"]["current_Rm_RQy_germ_identity_proved"]
            and not report["scope"]["full_goal_complete"], "open construction scope changed")
    return dict(passed=True, serialized_only=True, mathematical_provenance_passed=result["passed"],
        sides={side:dict(passed=v.get("passed"), proved_rows=sum(r.get("status") == "proved" for r in v.get("rows", [])),
                         unresolved=v.get("unresolved", []), error=v.get("error")) for side, v in result["sides"].items()},
        bytes=REPORT.stat().st_size, sha256=sha(REPORT), accepted_or_new_operator_replayed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "provenance", "write-report", "inspect-report"):
        group.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        answer = dict(passed=True, AST_LF_EOF_only=True)
    else:
        answer = globals()[next(name for name, value in vars(args).items() if value)]()
    print(json.dumps(answer))
