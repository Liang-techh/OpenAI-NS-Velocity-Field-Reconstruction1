"""Run a new actual endpoint second-radial source consumer exactly once."""
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
PARENT = PREFIX + "Rm_raw_native_jet_identity.json"
PARENT_HASH = "00abde807d6ce5e693cf8eef9233b0fdb135687353edfc6bb6a7472148763317"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
KERNEL = "src/openai_ns_reconstruction/original_Rm_second_radial_identity.py"
KERNEL_HASH = "8611b01445fafbe2aae4d89f87941d6c4827c1965b6cfd91e395b76e7214908b"
FIXTURE = "experiments/root_st073/original_Rm_second_radial_identity_fixture.json"
FIXTURE_HASH = "4f86bd110c22dae6abbae4c8ecb1623223d2f8d4a352284af3ec8b6d6e1c3934"
RUNNER = "scripts/run_original_Rm_second_radial_identity.py"
REPORT = BASE / (PREFIX + "Rm_second_radial_identity.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_second_radial_identity"
SOURCE = CHECKPOINTS / "new_source.json"
ENTERED = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_second_radial_identity.json"
ENDPOINT = CHECKPOINTS / "new_second_radial_endpoint_functions.json"


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                   ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def path_for(name):
    path = Path(name)
    require(not path.is_absolute(), "absolute input path")
    if ".." in path.parts:
        require(name == "../../src/openai_ns_reconstruction/schedule_pressure.py", "unexpected traversal")
        return (BASE / path).resolve()
    return ROOT / path if path.parts[0] in ("experiments", "src", "scripts", "tmp") else BASE / path.name


def read_pinned(path, expected):
    raw = Path(path).read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, "changed pinned input: " + str(path))
    return json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)


def read_saved(path):
    return json.loads(path.read_bytes())


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write((json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode())


def preflight():
    for relative in (RUNNER, KERNEL):
        raw = (ROOT / relative).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + relative)
        ast.parse(raw, filename=relative)
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "new derivative consumer changed")
    fixture = read_pinned(ROOT / FIXTURE, FIXTURE_HASH)
    require(fixture.get("passed") is True or fixture.get("toy_harness", {}).get("passed") is True,
            "new derivative fixture not admitted")


def setup():
    preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: " + name)
    finite = read_pinned(BASE / FINITE, parent["input_hashes"][FINITE])
    nodes = finite["exact_graph_nodes"]
    radial = parent["actual_Rm_native_radial_endpoint_source"]
    require(len(nodes) == 205147 and digest(nodes) == radial["base_digest"], "base source changed")
    nodes = nodes + copy.deepcopy(radial["appended_nodes"])
    require(len(nodes) == 205153 and digest(nodes) == radial["full_nodes_digest"], "saved endpoint source changed")
    seed = parent["actual_Rm_leading_seed_admission"]
    raw = parent["actual_Rm_raw_native_jet_identity"]
    require(seed["passed"] is True and raw["passed"] is True and not raw["RQy_identity_proved"],
            "leading seed/raw native premise not admitted")
    inputs = parent["actual_Rm_raw_native_jet_inputs"]
    left, right = inputs["modified_radial_background"], inputs["patch_radial_background"]
    require(left["E"][:3] == [85146, 85149, 85153] and left["V"][:3] == [0, 0, 0]
            and right["E"][:3] == [195164, 195165, 195166] and right["V"][:3] == [195171, 195172, 195173],
            "actual first-radial input ports changed")
    pairs = []
    for row in seed["fields"]:
        if row["field"] not in ("E_y", "V_y"):
            continue
        require(row["status"] == "proved", "radial source pair not proved")
        for port, variable, body in ((row["inlet_node"], 84528, left["E" if row["field"] == "E_y" else "V"][row["order"]]),
                (row["patch_node"], 185676, right["E" if row["field"] == "E_y" else "V"][row["order"]])):
            endpoint = nodes[port]
            require(endpoint["operation"] == "function_substitution" and endpoint["expression"] == body
                    and endpoint["variable"] == variable and endpoint["value"] == 1
                    and endpoint["Z_independent_substitution"] is True, "radial pair scope changed")
        pairs.append(dict(field=row["field"], order=row["order"], left_node=row["inlet_node"],
            right_node=row["patch_node"], proof_ref=dict(parent_sha256=PARENT_HASH, seed_admission_digest=digest(seed))))
    require(len(pairs) == 6, "missing radial source pairs")
    controls = parent["actual_Rm_leading_seed_normalization"]["source_parameter"]["leading_control_rows"]
    ports = dict(postrestore=dict(native_variable=84528, endpoint_node=1, E_y=left["E"][:3], V_y=left["V"][:3]),
        patch=dict(native_variable=185676, endpoint_node=1, E_y=right["E"][:3], V_y=right["V"][:3]),
        frozen_native_constant_nodes=copy.deepcopy(controls),
        saved_flat_premise=copy.deepcopy(parent["actual_Rm_leading_seed_inputs"]["saved_flat_certificate"]),
        accepted_Ey_Vy_source_pairs=pairs, ordinary_Z_orders=[0, 1, 2],
        physical_native_endpoint_Jacobians_exact_one=True,
        physical_directions=dict(postrestore="D_y p=1", patch="D_y x=x; x=1"))
    packet = dict(schema="actual_Rm_second_radial_inputs_v1", parent_report=PARENT,
        parent_sha256=PARENT_HASH, source_only=True, graph_node_count=len(nodes), graph_digest=digest(nodes),
        accepted_seed_digest=digest(seed), accepted_raw_jet_digest=digest(raw), ports=ports,
        accepted_ancestors_replayed=False, all_old_terms_and_errors_preserved=True)
    return parent, nodes, packet


def prepare_source():
    parent, nodes, packet = setup()
    require(not any(path.exists() for path in (SOURCE, ENTERED, RESULT, ENDPOINT)), "source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE), graph_node_count=len(nodes))


def compute_new():
    parent, nodes, packet = setup()
    require(read_saved(SOURCE) == packet and not ENTERED.exists() and not RESULT.exists(), "source changed or entered; no replay")
    spec = importlib.util.spec_from_file_location("new_Rm_second_radial", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTERED, dict(source_sha256=sha(SOURCE), kernel_sha256=KERNEL_HASH))
    started = time.perf_counter()
    result = module.prove(nodes, packet["ports"], expected_nodes_digest=packet["graph_digest"],
                          accepted_seed_digest=packet["accepted_seed_digest"])
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
                  entry_checkpoint_sha256=sha(ENTERED), actual_saved_graph_run=True, actual_saved_graph_loaded=True)
    result["symbolic_runtime"] = dict(sympy_version=module._sp.__version__)
    result["scope"]["no_saved_graph_execution"] = False
    result["scope"]["accepted_ancestor_operator_executed"] = False
    write_new(RESULT, result)
    return dict(passed=result.get("passed"), result_sha256=sha(RESULT), unresolved_gates=result.get("unresolved_gates"))


def define_endpoint_functions():
    parent, nodes, packet = setup()
    result = read_saved(RESULT)
    require(result["passed"] is True and result["source_checkpoint_sha256"] == sha(SOURCE), "second-radial identity not admitted")
    require(not ENDPOINT.exists(), "endpoint definitions already saved")
    started = time.perf_counter()
    delta = [dict(operation="exact_rational", numerator="1", denominator="10")]
    tenth = len(nodes)
    bindings = {side: dict(E_yy=[], V_yy=[0, 0, 0]) for side in ("postrestore", "patch")}
    for side, key in (("postrestore", "left_node"), ("patch", "right_node")):
        for pair in packet["ports"]["accepted_Ey_Vy_source_pairs"]:
            if pair["field"] != "E_y":
                continue
            bindings[side]["E_yy"].append(len(nodes)+len(delta))
            delta.append(dict(operation="product", arguments=[tenth, pair[key]],
                certified_actual_endpoint_native_second_derivative=True, ordinary_slow_Z_derivative_order=pair["order"],
                actual_derivative_identity_sha256=sha(RESULT), no_extra_physical_scale=True))
    value = dict(passed=True, source_only=True, base_node_count=len(nodes), base_digest=packet["graph_digest"],
        appended_nodes=delta, endpoint_bindings=bindings, full_node_count=len(nodes)+len(delta),
        full_nodes_digest=digest(nodes+delta), identity_checkpoint_sha256=sha(RESULT),
        numerical_values_computed=False, execution_seconds=time.perf_counter()-started)
    write_new(ENDPOINT, value)
    return dict(passed=True, sha256=sha(ENDPOINT), full_node_count=value["full_node_count"])


def write_report():
    parent, nodes, expected = setup()
    packet, result, endpoint = read_saved(SOURCE), read_saved(RESULT), read_saved(ENDPOINT)
    require(packet == expected and result["passed"] is True and endpoint["passed"] is True, "new component not admitted")
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_leading_second_radial_identity", parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        actual_Rm_second_radial_inputs=packet, actual_Rm_second_radial_identity=result,
        actual_Rm_second_radial_endpoint_functions=endpoint,
        second_radial_checkpoint_hashes={path.name:sha(path) for path in (SOURCE, ENTERED, RESULT, ENDPOINT)})
    report["execution_seconds"].update(new_actual_second_radial_identity=result["execution_seconds"],
                                       new_second_radial_endpoint_functions=endpoint["execution_seconds"])
    report["input_hashes"].update({PARENT:PARENT_HASH, FIXTURE:FIXTURE_HASH})
    report["implementation_hashes"].update({KERNEL:KERNEL_HASH, RUNNER:sha(ROOT / RUNNER)})
    report["scope"].update(actual_leading_second_radial_identity=True, modified_finite_N_germ=False,
        current_Rm_RQy_germ_identity_proved=False, current_Rm_to_Rp_join=False, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, expected = setup()
    report, result, endpoint = read_saved(REPORT), read_saved(RESULT), read_saved(ENDPOINT)
    require(report["actual_Rm_second_radial_inputs"] == read_saved(SOURCE) == expected
            and report["actual_Rm_second_radial_identity"] == result
            and report["actual_Rm_second_radial_endpoint_functions"] == endpoint, "saved packets changed")
    require(report["second_radial_checkpoint_hashes"] == {path.name:sha(path) for path in (SOURCE, ENTERED, RESULT, ENDPOINT)}
            and endpoint["identity_checkpoint_sha256"] == sha(RESULT)
            and digest(nodes+endpoint["appended_nodes"]) == endpoint["full_nodes_digest"], "new source/proof link changed")
    for key in ("actual_Rm_leading_seed_admission", "actual_Rm_raw_native_jet_identity", "current_Rm_pressure_inlet_identity",
                "actual_Rm_absolute_m15", "corrected_own_swirl_provider", "corrected_reserved_point_numeric"):
        require(report[key] == parent[key], "inherited packet changed: " + key)
    for name, expected_hash in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected_hash, "report pin changed: " + name)
    require(report["scope"]["full_goal_complete"] is False and report["current_Rm_RQy_germ_identity_proved"] is False,
            "invalid physical seam/completion claim")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
        old_mathematics_executed=False, new_operator_replayed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "compute-new", "define-endpoint-functions", "write-report", "inspect-report"):
        group.add_argument("--"+name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        print(json.dumps(dict(passed=True, source_preflight_only=True)))
    else:
        action = (prepare_source if args.prepare_source else compute_new if args.compute_new
                  else define_endpoint_functions if args.define_endpoint_functions else write_report if args.write_report else inspect_report)
        print(json.dumps(action(), indent=2))
