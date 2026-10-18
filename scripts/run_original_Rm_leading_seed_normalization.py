"""Prepare saved source inputs and run a new Rm leading seed normalizer once.

Source preparation, report packaging and inspection are data-only. Accepted
owners, mean proofs, producers and numerical providers are never imported.
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
PARENT = PREFIX + "Rm_pressure_inlet_identity.json"
PARENT_HASH = "2a95b27e9851c225592c11e4481e3f8f4c92af317972d6521e45b29f1145e18b"
C6 = PREFIX + "prescribed_R100_Rm_C6_functions.json.gz"
LEADING = PREFIX + "true_Rm_patch_C6_functions.json.gz"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
MEAN = PREFIX + "Rm_inlet_m15_numeric.json"
KERNEL = "src/openai_ns_reconstruction/original_Rm_leading_seed_normalization.py"
KERNEL_HASH = "5ea9cc6f03989d5d466fa246b6f0e93cbc29de89fdf353be9d15119de3120b53"
FIXTURE = "experiments/root_st073/original_Rm_leading_seed_normalization_fixture.json"
FIXTURE_HASH = "c66948fb5e6001a41233b8c400dddb1fdbcf6c9e12ec443c2d07e495678bf86d"
RUNNER = "scripts/run_original_Rm_leading_seed_normalization.py"
REPORT = BASE / (PREFIX + "Rm_leading_seed_normalization.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_leading_seed_normalization"
SOURCE = CHECKPOINTS / "new_source.json"
RADIAL = CHECKPOINTS / "new_radial_endpoint.json"
ENTERED = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_normalization.json"
PRESSURE_KERNEL = "src/openai_ns_reconstruction/original_Rm_seed_pressure_composition.py"
PRESSURE_KERNEL_HASH = "13a42608cee955c1e4a16e253de2959e718c7c28e056c991496804d80ae89c2b"
PRESSURE_FIXTURE = "experiments/root_st073/original_Rm_seed_pressure_composition_fixture.json"
PRESSURE_FIXTURE_HASH = "290b2e55840f19ead2b4c0d3170bafc05f6e6e8fea4b3a117d49d27f1412a0ee"
PRESSURE_ENTERED = CHECKPOINTS / "new_pressure_composition_entered.json"
PRESSURE_RESULT = CHECKPOINTS / "new_pressure_composition.json"


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                   ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


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


def source_preflight():
    for relative in (RUNNER, KERNEL):
        raw = (ROOT / relative).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: " + relative)
        ast.parse(raw, filename=relative)
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "unreviewed new normalizer")
    fixture = read_pinned(ROOT / FIXTURE, FIXTURE_HASH)
    require(fixture.get("passed") is True or fixture.get("toy_harness", {}).get("passed") is True,
            "new fixture witness not admitted")


def setup():
    source_preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: " + name)
    c6, leading, finite, old = (read_pinned(BASE / name, parent["input_hashes"][name]) for name in (C6, LEADING, FINITE, MEAN))
    nodes = finite["exact_graph_nodes"]
    require(len(c6["exact_graph_nodes"]) == 97295 and len(leading["exact_graph_nodes"]) == 196042 and len(nodes) == 205147
            and nodes[:196042] == leading["exact_graph_nodes"]
            and nodes[:97295] == c6["exact_graph_nodes"], "prescribed/leading prefix changed")
    require(c6["source_family"] == leading["source_family"] == finite["source_family"] == old["source_family"] == parent["source_family"],
            "source family changed")
    cb, lb = c6["prescribed_R100_Rm_C6_bindings"], leading["true_Rm_patch_C6_bindings"]
    require(lb["inlet"] == cb["endpoints"]["Rm"], "patch uses another prescribed inlet")
    require(lb["inlet"]["fields"]["E"][:3] == [85803, 85804, 85805]
            and lb["am"][:3] == [183702, 183705, 183710], "amplitude ports changed")
    require(lb["inlet"]["fields"]["m"][:3] == [97211, 97212, 97213]
            and lb["chart"]["fields"]["m"][:3] == [195000, 195001, 195002], "saved mean scope changed")
    post = cb["charts"]["post_restore"]
    require(post["phase"] == 84528 and post["physical_measure"] == 1
            and post["radial"]["E_y"][:3] == [85146, 85149, 85153]
            and post["radial"]["V_y"][:3] == [0, 0, 0], "actual native radial source changed")
    require(cb["F0"][:3] == [35657, 35658, 35661], "F0 ordinary source ports changed")
    require(nodes[33572]["operation"] == "analytic_unary" and nodes[33572]["name"] == "exp"
            and nodes[33572]["argument"] == 5 and nodes[5]["operation"] == "sum"
            and nodes[5]["arguments"] == [3, 4], "frozen logP source changed")
    for order, body in enumerate((32693, 32694, 32695)):
        recipe = nodes[body]
        require(recipe["operation"] == "current_original_leading_function_recipe"
                and recipe["quantity"] == "P0" and recipe["Z_variable"] == 20
                and recipe["coordinate_independent"] is True
                and recipe["phase_independent_leading_source"] is True
                and recipe["defining_quantity_not_a_range_value"] is True
                and recipe["exact_source_op_P0_function_not_enclosure_endpoint"] is True,
                "P0 axial source recipe changed")
        for port, variable in ((cb["endpoints"]["Rm"]["fields"]["P0"][order], 84528),
                               (lb["endpoints"]["Rm"]["fields"]["P0"][order], 185676)):
            wrapper = nodes[port]
            require(wrapper["operation"] == "function_substitution" and wrapper["expression"] == body
                    and wrapper["variable"] == variable and wrapper["value"] == 1
                    and wrapper["Z_independent_substitution"] is True, "P0 wrapper source changed")
    mean = old["leading_mean_certificate"]
    flat = old["compact_flat_and_history_certificate"]
    require(mean["passed"] is True and mean["ordinary_Z_orders"] == [0, 1, 2]
            and mean["source_only"] is True and mean["numeric_caps_used"] is False
            and mean["global_opaque_alias_used"] is False and flat["passed"] is True,
            "saved mean/flat source premise not admitted")
    require([(row["left_node"], row["right_node"]) for row in mean["rows"]]
            == [(96436, 195000), (96437, 195001), (96438, 195002)]
            and all(row["passed"] is True for row in mean["rows"]), "saved mean roots changed")
    for order, row in enumerate(mean["rows"]):
        require(row["ordinary_Z_order"] == order, "saved mean order changed")
        for port, variable, body in ((lb["inlet"]["fields"]["m"][order], 84528, row["left_node"]),
                (lb["endpoints"]["Rm"]["fields"]["m"][order], 185676, row["right_node"])):
            endpoint = nodes[port]
            require(endpoint["operation"] == "function_substitution" and endpoint["expression"] == body
                    and endpoint["variable"] == variable and endpoint["value"] == 1
                    and endpoint["Z_independent_substitution"] is True, "saved mean endpoint scope changed")
    source = dict(schema="Rm_leading_seed_normalization_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, source_family=copy.deepcopy(parent["source_family"]),
        graph_node_count=len(nodes), graph_digest=digest(nodes), shared_prescribed_prefix_digest=digest(nodes[:97295]),
        c6_sha256=parent["input_hashes"][C6], leading_sha256=parent["input_hashes"][LEADING], finite_sha256=parent["input_hashes"][FINITE],
        saved_mean_report_sha256=parent["input_hashes"][MEAN], saved_mean_certificate=copy.deepcopy(mean),
        saved_flat_certificate=copy.deepcopy(flat), saved_mean_consumed_not_reexecuted=True,
        leading_inlet=copy.deepcopy(lb["inlet"]), leading_endpoint=copy.deepcopy(lb["endpoints"]["Rm"]),
        leading_am=copy.deepcopy(lb["am"]), postrestore=copy.deepcopy(cb["charts"]["post_restore"]),
        ordinary_Z_orders=[0, 1, 2], all_old_errors_and_212_terms_unchanged=True,
        current_RQy_defect_packet_digest=digest(parent["current_Rm_pressure_inlet_identity"]),
        RQy_current_and_patch_not_yet_aliased=True, accepted_ancestors_replayed=False)
    return parent, nodes, cb, lb, source


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write((json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode())


def read_checkpoint(path):
    return json.loads(path.read_bytes())


def prepare_source():
    parent, nodes, cb, lb, source = setup()
    require(not any(path.exists() for path in (SOURCE, RADIAL, ENTERED, RESULT)), "checkpoints already exist")
    write_new(SOURCE, source)
    started = time.perf_counter()
    delta, roots = [], {}
    for name in ("E_y", "V_y"):
        roots[name] = []
        for order in range(3):
            roots[name].append(len(nodes) + len(delta))
            delta.append(dict(operation="function_substitution",
                expression=source["postrestore"]["radial"][name][order],
                variable=84528, value=1, Z_independent_substitution=True,
                actual_current_Rm_postrestore_radial_source=True,
                ordinary_slow_Z_derivative_order=order, physical_native_Jacobian_exact_one=True))
    radial = dict(passed=True, source_only=True, base_node_count=len(nodes), base_digest=source["graph_digest"],
        appended_nodes=delta, radial_endpoint_roots=roots, full_node_count=len(nodes) + len(delta),
        full_nodes_digest=digest(nodes + delta), native_postrestore_coordinate=84528, endpoint_value=1,
        physical_native_Jacobian_exact_one=True, additional_Pstar_or_factorial=False,
        numerical_value_computed=False, execution_seconds=time.perf_counter()-started,
        source_checkpoint_sha256=sha(SOURCE))
    write_new(RADIAL, radial)
    return dict(passed=True, source_sha256=sha(SOURCE), radial_sha256=sha(RADIAL),
                new_source_nodes=len(delta), full_graph_node_count=radial["full_node_count"])


def compute_new():
    parent, nodes, cb, lb, expected = setup()
    source = read_checkpoint(SOURCE)
    require(source == expected and not ENTERED.exists() and not RESULT.exists(), "source mismatch or operator entered; do not replay")
    radial = read_checkpoint(RADIAL)
    require(radial["base_digest"] == source["graph_digest"] and radial["source_checkpoint_sha256"] == sha(SOURCE),
            "radial source checkpoint changed")
    nodes = nodes + copy.deepcopy(radial["appended_nodes"])
    require(digest(nodes) == radial["full_nodes_digest"], "new native radial graph changed")
    cb = copy.deepcopy(cb)
    cb["endpoints"]["Rm"]["radial"] = copy.deepcopy(radial["radial_endpoint_roots"])
    spec = importlib.util.spec_from_file_location("new_Rm_leading_seed_normalization", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTERED, dict(source_sha256=sha(SOURCE), kernel_sha256=KERNEL_HASH, expression_reduction_may_now_run=True))
    started = time.perf_counter()
    result = module.prove(nodes, {"true_Rm_patch_C6_bindings": lb, "prescribed_R100_Rm_C6_bindings": cb},
        expected_nodes_digest=radial["full_nodes_digest"], ordinary_orders=(0, 1, 2),
        saved_mean_premise=dict(report_sha256=source["saved_mean_report_sha256"],
                                mean=source["saved_mean_certificate"], flat=source["saved_flat_certificate"]))
    result["execution_seconds"] = time.perf_counter() - started
    result["actual_saved_graph_run"] = True
    result["actual_saved_graph_loaded"] = True
    result["source_checkpoint_sha256"] = sha(SOURCE)
    result["operator_entry_sha256"] = sha(ENTERED)
    result["native_radial_checkpoint_sha256"] = sha(RADIAL)
    write_new(RESULT, result)
    return dict(passed=result.get("passed"), source_sha256=sha(SOURCE), result_sha256=sha(RESULT))


def compute_pressure_composition():
    parent, nodes, cb, lb, expected = setup()
    require(read_checkpoint(SOURCE) == expected and not PRESSURE_ENTERED.exists()
            and not PRESSURE_RESULT.exists(), "pressure component already entered; no replay")
    raw = (ROOT / PRESSURE_KERNEL).read_bytes()
    require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "pressure source LF/EOF")
    ast.parse(raw, filename=PRESSURE_KERNEL)
    require(sha(ROOT / PRESSURE_KERNEL) == PRESSURE_KERNEL_HASH, "pressure component source changed")
    fixture = read_pinned(ROOT / PRESSURE_FIXTURE, PRESSURE_FIXTURE_HASH)
    require(fixture.get("passed") is True or fixture.get("toy_harness", {}).get("passed") is True,
            "pressure component fixture not admitted")
    seed = read_checkpoint(RESULT)
    require(sha(RESULT) == "f0dcd69cc9c6e1ecf332787b49fc8d12fb26e3d934e94b646a43396f4111bd50",
            "partial seed checkpoint changed")
    radial = read_checkpoint(RADIAL)
    nodes = nodes + copy.deepcopy(radial["appended_nodes"])
    spec = importlib.util.spec_from_file_location("new_Rm_seed_pressure_composition", ROOT / PRESSURE_KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(PRESSURE_ENTERED, dict(kernel_sha256=PRESSURE_KERNEL_HASH, seed_checkpoint_sha256=sha(RESULT)))
    started = time.perf_counter()
    result = module.prove(nodes, seed, expected_nodes_digest=radial["full_nodes_digest"])
    result.update(execution_seconds=time.perf_counter()-started, seed_checkpoint_sha256=sha(RESULT),
                  entry_checkpoint_sha256=sha(PRESSURE_ENTERED), actual_saved_graph_run=True, actual_saved_graph_loaded=True)
    write_new(PRESSURE_RESULT, result)
    return dict(passed=result["passed"], sha256=sha(PRESSURE_RESULT))


def write_report():
    parent, nodes, cb, lb, expected = setup()
    source, result = read_checkpoint(SOURCE), read_checkpoint(RESULT)
    radial = read_checkpoint(RADIAL)
    pressure = read_checkpoint(PRESSURE_RESULT)
    require(source == expected, "source checkpoint changed")
    require(pressure["passed"] is True, "pressure composition not proved")
    rows = copy.deepcopy(result["fields"])
    for row in rows:
        if row["field"] == "pressure":
            require(row["status"] == "unresolved", "partial checkpoint pressure row changed")
            row.pop("error", None)
            row.pop("gate", None)
            row["status"] = "proved_by_new_pressure_composition"
            row["composition_checkpoint_sha256"] = sha(PRESSURE_RESULT)
        else:
            require(row["status"] == ("external_premise" if row["field"] == "m" else "proved"), "unproved seed row")
    admission = dict(passed=True, ordinary_Z_orders=[0, 1, 2], fields=rows,
        partial_seed_checkpoint_sha256=sha(RESULT), pressure_composition_checkpoint_sha256=sha(PRESSURE_RESULT),
        original_partial_result_preserved=True, old_seed_operator_replayed=False)
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_leading_seed_normalization",
        parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        actual_Rm_leading_seed_inputs=source, actual_Rm_leading_seed_normalization=result,
        actual_Rm_native_radial_endpoint_source=radial,
        actual_Rm_seed_pressure_composition=pressure, actual_Rm_leading_seed_admission=admission,
        execution_seconds=dict(new_native_radial_endpoint_source=radial["execution_seconds"],
                               new_actual_Rm_leading_seed_normalization=result["execution_seconds"],
                               new_pressure_composition=pressure["execution_seconds"]),
        checkpoint_hashes={path.name: sha(path) for path in (SOURCE, RADIAL, ENTERED, RESULT, PRESSURE_ENTERED, PRESSURE_RESULT)},
        current_Rm_RQy_germ_identity_proved=False, current_Rm_to_Rp_join=False,
        new_velocity_numeric_values_computed=False, old_212_terms_and_all_errors_unchanged=True)
    report["input_hashes"] = dict(parent["input_hashes"], **{PARENT: PARENT_HASH, FIXTURE: FIXTURE_HASH,
                                                           PRESSURE_FIXTURE: PRESSURE_FIXTURE_HASH})
    report["implementation_hashes"] = dict(parent["implementation_hashes"], **{KERNEL: KERNEL_HASH,
        PRESSURE_KERNEL: PRESSURE_KERNEL_HASH, RUNNER: sha(ROOT / RUNNER)})
    report["scope"] = dict(parent["scope"], current_Rm_RQy_germ_identity_proved=False,
        current_Rm_to_Rp_join=False, new_velocity_numeric_values_computed=False,
        global_five_moment_closure=False, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, cb, lb, expected = setup()
    report = read_checkpoint(REPORT)
    source, result = read_checkpoint(SOURCE), read_checkpoint(RESULT)
    radial = read_checkpoint(RADIAL)
    require(source == expected and report["actual_Rm_leading_seed_inputs"] == source
            and report["actual_Rm_leading_seed_normalization"] == result, "saved source/result changed")
    require(report["actual_Rm_native_radial_endpoint_source"] == radial
            and report["checkpoint_hashes"] == {path.name: sha(path) for path in (SOURCE, RADIAL, ENTERED, RESULT, PRESSURE_ENTERED, PRESSURE_RESULT)}
            and result["source_checkpoint_sha256"] == sha(SOURCE)
            and result["operator_entry_sha256"] == sha(ENTERED)
            and result["native_radial_checkpoint_sha256"] == sha(RADIAL), "checkpoint link changed")
    require(report["actual_Rm_seed_pressure_composition"] == read_checkpoint(PRESSURE_RESULT)
            and report["actual_Rm_leading_seed_admission"]["passed"] is True,
            "new pressure composition/admission changed")
    require(radial["base_node_count"] == len(nodes) == 205147 and len(radial["appended_nodes"]) == 6
            and digest(nodes + radial["appended_nodes"]) == radial["full_nodes_digest"], "native radial graph changed")
    pressure = read_checkpoint(PRESSURE_RESULT)
    require(pressure["seed_checkpoint_sha256"] == sha(RESULT)
            and pressure["entry_checkpoint_sha256"] == sha(PRESSURE_ENTERED)
            and pressure["passed"] is True and len(pressure["fields"]) == 3,
            "pressure component link changed")
    for admitted, original in zip(report["actual_Rm_leading_seed_admission"]["fields"], result["fields"]):
        require(all(admitted[key] == original[key] for key in ("field", "order", "patch_node", "inlet_node")),
                "admission source pair changed")
        if original["field"] == "pressure":
            proved = pressure["fields"][original["order"]]
            require(proved["status"] == "proved" and all(proved[key] == original[key]
                    for key in ("field", "order", "patch_node", "inlet_node"))
                    and admitted["composition_checkpoint_sha256"] == sha(PRESSURE_RESULT),
                    "pressure admission scope changed")
        else:
            require(admitted == original, "already proved seed row changed")
    for offset, row in enumerate(radial["appended_nodes"]):
        name, order = ("E_y", "V_y")[offset // 3], offset % 3
        require(row["operation"] == "function_substitution" and row["variable"] == 84528 and row["value"] == 1
                and row["Z_independent_substitution"] is True
                and row["expression"] == source["postrestore"]["radial"][name][order]
                and row["ordinary_slow_Z_derivative_order"] == order, "native radial source operand changed")
    for key in ("current_Rm_pressure_inlet_identity", "actual_Rm_absolute_m15", "corrected_own_swirl_provider",
                "corrected_reserved_point_numeric"):
        require(report[key] == parent[key], "inherited packet changed: " + key)
    for name, expected_hash in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected_hash, "report pin changed: " + name)
    require(report["scope"]["full_goal_complete"] is False
            and report["current_Rm_RQy_germ_identity_proved"] is False, "invalid physical join/completion claim")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
        source_checkpoint_sha256=sha(SOURCE), result_checkpoint_sha256=sha(RESULT),
        old_mathematics_executed=False, new_actual_operator_replayed=False,
        current_Rm_RQy_germ_identity_proved=False, current_Rm_to_Rp_join=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--preflight", action="store_true")
    group.add_argument("--prepare-source", action="store_true")
    group.add_argument("--compute-new", action="store_true")
    group.add_argument("--compute-pressure-composition", action="store_true")
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        source_preflight()
        print(json.dumps(dict(passed=True, source_preflight_only=True, old_mathematics_executed=False)))
    else:
        action = (prepare_source if args.prepare_source else compute_new if args.compute_new
                  else compute_pressure_composition if args.compute_pressure_composition
                  else write_report if args.write_report else inspect_report)
        print(json.dumps(action(), indent=2))
