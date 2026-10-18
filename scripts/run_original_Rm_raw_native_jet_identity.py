"""Consume a saved Rm seed proof and reduce only the new raw native-y jet."""
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
PARENT = PREFIX + "Rm_leading_seed_normalization.json"
PARENT_HASH = "794f5e4412c099512d1012b8de8a6d65959417c30a3c881b1fc59c200abc2d83"
Q4 = PREFIX + "modified_Q4_pressure_C3_functions.json.gz"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
KERNEL = "src/openai_ns_reconstruction/original_Rm_raw_native_jet_identity.py"
KERNEL_HASH = "cd4b9d76ad3f7424e457e6986d373e87a94fd357e2037d42590073f4b1e867b5"
FIXTURE = "experiments/root_st073/original_Rm_raw_native_jet_identity_fixture.json"
FIXTURE_HASH = "dee5654d5fa8cdf492064395f29a1b2d0d15900c531815726f2c85128cb1ad98"
RUNNER = "scripts/run_original_Rm_raw_native_jet_identity.py"
REPORT = BASE / (PREFIX + "Rm_raw_native_jet_identity.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_raw_native_jet_identity"
SOURCE = CHECKPOINTS / "new_source.json"
ENTERED = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_raw_native_jet.json"


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


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write((json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode())


def read_saved(path):
    return json.loads(path.read_bytes())


def preflight():
    for relative in (RUNNER, KERNEL):
        raw = (ROOT / relative).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + relative)
        ast.parse(raw, filename=relative)
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "unreviewed raw-jet reducer")
    fixture = read_pinned(ROOT / FIXTURE, FIXTURE_HASH)
    require(fixture.get("passed") is True or fixture.get("toy_harness", {}).get("passed") is True
            or fixture.get("harness_checkpoint", {}).get("result", {}).get("passed") is True,
            "fixture did not pass")


def setup():
    preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: " + name)
    q4, finite = (read_pinned(BASE / name, parent["input_hashes"][name]) for name in (Q4, FINITE))
    nodes = finite["exact_graph_nodes"]
    radial = parent["actual_Rm_native_radial_endpoint_source"]
    require(len(nodes) == 205147 and nodes[:183700] == q4["exact_graph_nodes"]
            and digest(nodes) == radial["base_digest"], "actual graph prefix changed")
    nodes = nodes + copy.deepcopy(radial["appended_nodes"])
    require(len(nodes) == 205153 and digest(nodes) == radial["full_nodes_digest"], "saved endpoint delta changed")
    seed = parent["actual_Rm_leading_seed_normalization"]
    admission = parent["actual_Rm_leading_seed_admission"]
    require(admission["passed"] is True and seed["expected_nodes_digest"] == digest(nodes)
            and seed["scope"]["nodes_input_unmutated"] is True
            and seed["accepted_m_reproved"] is False, "seed premise not admitted")
    source = parent["actual_Rm_leading_seed_inputs"]
    require(source["postrestore"]["phase"] == 84528 and source["postrestore"]["physical_measure"] == 1,
            "postrestore native source metadata changed")
    aliases = []
    for index, row in enumerate(admission["fields"]):
        require(row["status"] == ("external_premise" if row["field"] == "m" else
                "proved_by_new_pressure_composition" if row["field"] == "pressure" else "proved"), "unproved seed row")
        left, right = nodes[row["inlet_node"]], nodes[row["patch_node"]]
        for endpoint, variable in ((left, 84528), (right, 185676)):
            require(endpoint["operation"] == "function_substitution"
                    and endpoint["variable"] == variable and endpoint["value"] == 1
                    and endpoint["Z_independent_substitution"] is True, "seed endpoint scope changed")
        if row["field"] in ("E_y", "V_y"):
            require(row["inlet_node"] == radial["radial_endpoint_roots"][row["field"]][row["order"]]
                    and left["actual_current_Rm_postrestore_radial_source"] is True
                    and left["ordinary_slow_Z_derivative_order"] == row["order"]
                    and left["physical_native_Jacobian_exact_one"] is True,
                    "new native radial source alias changed")
        if row["field"] == "m":
            mean = source["saved_mean_certificate"]["rows"][row["order"]]
            require((left["expression"], right["expression"])
                    == (mean["left_node"], mean["right_node"]) and mean["passed"] is True and mean["ordinary_Z_order"] == row["order"],
                    "saved mean does not prove this endpoint pair")
        aliases.append(dict(field=row["field"], order=row["order"], left_node=row["inlet_node"],
            right_node=row["patch_node"], proof_ref=dict(parent_sha256=PARENT_HASH,
                seed_proof_digest=digest(admission), field_row=index, admission_digest=digest(admission),
                saved_mean_report_sha256=source["saved_mean_report_sha256"] if row["field"] == "m" else None)))
    require(len(aliases) == 33, "incomplete seed aliases")
    mb, fb = q4["modified_Q4_pressure_C3_bindings"], finite["Rm_patch_finite_N_C5_bindings"]
    packet = dict(schema="actual_Rm_raw_native_jet_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, seed_proof_digest=digest(admission),
        nodes_digest=digest(nodes), node_count=len(nodes), seed_aliases=aliases,
        modified_radial_background=copy.deepcopy(mb["charts"]["postrestore"]["radial_background"]),
        patch_radial_background=copy.deepcopy(fb["radial_background"]),
        ordinary_Z_orders=[0, 1, 2], physical_native_endpoint_Jacobians_exact_one=True,
        leading_raw_native_jet_only=True, finite_N_modified_germ_proved=False,
        old_reducers_providers_or_fixtures_replayed=False)
    return parent, nodes, mb, fb, packet


def prepare_source():
    *_, packet = setup()
    require(not any(path.exists() for path in (SOURCE, ENTERED, RESULT)), "source already prepared")
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE), seed_aliases=len(packet["seed_aliases"]))


def compute_new():
    parent, nodes, mb, fb, packet = setup()
    require(read_saved(SOURCE) == packet and not ENTERED.exists() and not RESULT.exists(), "source changed or operator entered; no replay")
    spec = importlib.util.spec_from_file_location("new_Rm_raw_native_jet", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTERED, dict(kernel_sha256=KERNEL_HASH, source_sha256=sha(SOURCE), expression_reduction_may_now_run=True))
    started = time.perf_counter()
    result = module.prove(nodes, {"modified_Q4_pressure_C3_bindings": mb},
        {"Rm_patch_finite_N_C5_bindings": fb}, packet["seed_aliases"],
        expected_nodes_digest=packet["nodes_digest"], seed_proof_digest=packet["seed_proof_digest"])
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
                  entry_checkpoint_sha256=sha(ENTERED), actual_saved_graph_run=True, actual_saved_graph_loaded=True)
    write_new(RESULT, result)
    return dict(passed=result["passed"], result_sha256=sha(RESULT), unresolved_gates=result["unresolved_gates"])


def write_report():
    parent, nodes, mb, fb, expected = setup()
    packet, result = read_saved(SOURCE), read_saved(RESULT)
    require(packet == expected, "saved native source changed")
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_leading_raw_native_jet_identity", parent_report=PARENT,
        parent_report_sha256=PARENT_HASH, actual_Rm_raw_native_jet_inputs=packet,
        actual_Rm_raw_native_jet_identity=result,
        raw_native_jet_checkpoint_hashes={path.name:sha(path) for path in (SOURCE, ENTERED, RESULT)})
    report["execution_seconds"]["new_actual_Rm_raw_native_jet_identity"] = result["execution_seconds"]
    report["input_hashes"].update({PARENT:PARENT_HASH, FIXTURE:FIXTURE_HASH})
    report["implementation_hashes"].update({KERNEL:KERNEL_HASH, RUNNER:sha(ROOT / RUNNER)})
    report["scope"].update(current_Rm_RQy_germ_identity_proved=False, current_Rm_to_Rp_join=False,
        leading_raw_native_jet_only=True, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, mb, fb, expected = setup()
    report, result = read_saved(REPORT), read_saved(RESULT)
    require(report["actual_Rm_raw_native_jet_inputs"] == read_saved(SOURCE) == expected
            and report["actual_Rm_raw_native_jet_identity"] == result, "saved inputs/result changed")
    require(report["raw_native_jet_checkpoint_hashes"] == {path.name:sha(path) for path in (SOURCE, ENTERED, RESULT)}
            and result["source_checkpoint_sha256"] == sha(SOURCE)
            and result["entry_checkpoint_sha256"] == sha(ENTERED), "saved checkpoint links changed")
    for key in ("actual_Rm_leading_seed_normalization", "actual_Rm_leading_seed_admission",
                "actual_Rm_seed_pressure_composition", "current_Rm_pressure_inlet_identity",
                "actual_Rm_absolute_m15", "corrected_own_swirl_provider", "corrected_reserved_point_numeric"):
        require(report[key] == parent[key], "inherited packet changed: " + key)
    for name, expected_hash in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected_hash, "report pin changed: " + name)
    require(report["scope"]["full_goal_complete"] is False and result["RQy_identity_proved"] is False,
            "invalid germ/completion claim")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
                raw_native_rows_proved=sum(row["status"] == "proved" for row in result["row_forms"]),
                actual_operator_replayed=False, accepted_ancestors_replayed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "compute-new", "write-report", "inspect-report"):
        group.add_argument("--" + name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        print(json.dumps(dict(passed=True, source_preflight_only=True)))
    else:
        action = (prepare_source if args.prepare_source else compute_new if args.compute_new
                  else write_report if args.write_report else inspect_report)
        print(json.dumps(action(), indent=2))
