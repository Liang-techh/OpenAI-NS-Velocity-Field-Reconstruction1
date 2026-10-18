"""New source-only Rm pressure inlet and patch-increment consumer.

The actual operator runs only with --compute-new. Later report actions
consume checkpoints and never replay the endpoint reduction or ancestors.
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
PARENT = PREFIX + "corrected_own_swirl_provider.json"
PARENT_HASH = "bf4ab120edd8dc02925977dc37614215a2be00a83b8393d1d5763c9ad50cba42"
MODIFIED = PREFIX + "modified_Q4_pressure_C3_functions.json.gz"
LEADING = PREFIX + "true_Rm_patch_C6_functions.json.gz"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
KERNEL = "src/openai_ns_reconstruction/original_Rm_pressure_inlet_identity.py"
KERNEL_HASH = "ef554c0e4b6d3c6e20c90b1165984dfe0347004f51d64e7bfc31d54ceb8aeda8"
FIXTURE = "experiments/root_st073/original_Rm_pressure_inlet_identity_fixture.json"
FIXTURE_HASH = "23a04df253ef6fc5ecd7e961899874d158b5311af78f797f3c71118eefc85277"
RUNNER = "scripts/run_original_Rm_pressure_inlet_identity.py"
REPORT = BASE / (PREFIX + "Rm_pressure_inlet_identity.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_pressure_inlet_identity"
SOURCE = CHECKPOINTS / "new_source.json"
IDENTITY = CHECKPOINTS / "new_identity.json"
IMPORT_FAILURE = CHECKPOINTS / "new_import_failure.json"
LAYOUT_FAILURE = CHECKPOINTS / "new_layout_failure.json"
REDUCTION_ENTERED = CHECKPOINTS / "new_reduction_entered.json"


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
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "unreviewed new kernel")
    fixture = read_pinned(ROOT / FIXTURE, FIXTURE_HASH)
    witness = fixture.get("toy_harness", {})
    require(witness.get("passed") is True and witness.get("source_sha256") == KERNEL_HASH,
            "new fixture witness not admitted")


def setup():
    source_preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: " + name)
    modified, leading, finite = (read_pinned(BASE / name, parent["input_hashes"][name])
                                 for name in (MODIFIED, LEADING, FINITE))
    nodes = finite["exact_graph_nodes"]
    require(len(modified["exact_graph_nodes"]) == 183700 and len(leading["exact_graph_nodes"]) == 196042
            and len(nodes) == 205147, "saved graph counts changed")
    require(nodes[:183700] == modified["exact_graph_nodes"]
            and nodes[:196042] == leading["exact_graph_nodes"], "shared source prefixes differ")
    family = parent["corrected_own_swirl_source_inputs"]["source_family"]
    require(modified["source_family"] == leading["source_family"] == finite["source_family"] == family,
            "source family changed")
    current = parent["source_namespace"]
    require(current["external_Rm_handles"] == [183696, 183697]
            and parent["actual_restoration_pressure_adoption"]["adopted_pressure_join"]["external_Rm_retained"] == [183696, 183697],
            "current pressure namespace changed")
    mb, fb = modified["modified_Q4_pressure_C3_bindings"], finite["Rm_patch_finite_N_C5_bindings"]
    require(mb["Rm"]["boundary"]["R_Q_y"][:2] == [183696, 183697]
            and fb["endpoints"]["Rm"]["boundary"]["R_Q_y"][:2] == [204932, 204933], "boundary ports changed")
    require(mb["Rm"]["integrals"]["S1"][:2] == [183680, 183681]
            and mb["Rm"]["pressure"]["m15"][:2] == [183688, 183689]
            and mb["initial"]["R_Q_y"][:2] == [129567, 129568], "current pressure operands changed")
    require(fb["endpoints"]["Rm"]["pressure"]["m15"][:2] == [204924, 204925]
            and fb["endpoints"]["Rm"]["integrals"]["S1"][:2] == [204916, 204917]
            and fb["endpoints"]["Rh"]["pressure"]["m15"][:2] == [205135, 205136]
            and fb["endpoints"]["Rh"]["integrals"]["S1"][:2] == [205127, 205128]
            and fb["endpoints"]["Rh"]["boundary"]["R_Q_y"][:2] == [205143, 205144], "patch pressure operands changed")
    require(fb["incoming_integrals"] == mb["Rm"]["integrals"] and fb["initial"] == mb["initial"],
            "the patch reset its signed pressure memory or initial boundary")
    source = dict(schema="Rm_pressure_inlet_source_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, source_family=copy.deepcopy(family),
        modified_sha256=parent["input_hashes"][MODIFIED], leading_sha256=parent["input_hashes"][LEADING],
        finite_sha256=parent["input_hashes"][FINITE], graph_node_count=len(nodes), graph_digest=digest(nodes),
        modified_prefix_digest=digest(nodes[:183700]), leading_prefix_digest=digest(nodes[:196042]),
        exact_shared_prefixes=True, ordinary_Z_orders=[0, 1],
        current_pressure_namespace=copy.deepcopy(current),
        current_Rm_absolute_packet_digest=digest(parent["actual_Rm_absolute_m15"]),
        current_Rm_numeric_domain=copy.deepcopy(parent["axial_domain"]),
        current_boundary=copy.deepcopy(mb["Rm"]["boundary"]), current_pressure=copy.deepcopy(mb["Rm"]["pressure"]),
        current_integrals=copy.deepcopy(mb["Rm"]["integrals"]), initial=copy.deepcopy(mb["initial"]),
        patch_coordinate=fb["chart"]["coordinate"], patch_incoming_integrals=copy.deepcopy(fb["incoming_integrals"]),
        patch_local_integrals=copy.deepcopy(fb["local_integrals"]), patch_endpoints=copy.deepcopy(fb["endpoints"]),
        pressure_function_is_not_a_selected_error_cap=True, old_table_rows_not_added=True,
        all_212_forward_terms_and_errors_unchanged=True, own_swirl_registration_count=1,
        accepted_ancestors_replayed=False)
    return parent, nodes, mb, fb, source


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write((json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode())


def read_checkpoint(path):
    return json.loads(path.read_bytes())


def compute_new():
    parent, nodes, modified, finite, source = setup()
    require(not IDENTITY.exists() and not REDUCTION_ENTERED.exists(), "reduction was entered; do not replay")
    if SOURCE.exists():
        failure = read_checkpoint(IMPORT_FAILURE)
        require(failure.get("failed_during_module_import") is True
                and failure.get("new_actual_operator_entered") is False
                and failure.get("mathematics_executed") is False
                and failure.get("source_checkpoint_sha256") == sha(SOURCE)
                and read_checkpoint(SOURCE) == source,
                "a source checkpoint may resume only after the documented pre-math import failure")
        if LAYOUT_FAILURE.exists():
            layout_failure = read_checkpoint(LAYOUT_FAILURE)
            require(layout_failure.get("failed_before_reducer_construction") is True
                    and layout_failure.get("expression_reduction_executed") is False
                    and layout_failure.get("source_checkpoint_sha256") == sha(SOURCE),
                    "layout failure did not occur before expression reduction")
    else:
        write_new(SOURCE, source)
    spec = importlib.util.spec_from_file_location("new_Rm_pressure_identity", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(REDUCTION_ENTERED, dict(source_sha256=sha(SOURCE), kernel_sha256=KERNEL_HASH,
        explicit_archive_binding_wrappers=True, expression_reduction_may_now_run=True))
    started = time.perf_counter()
    result = module.prove(nodes, {"modified_Q4_pressure_C3_bindings": modified},
                          {"Rm_patch_finite_N_C5_bindings": finite},
                          expected_nodes_digest=source["graph_digest"], common_prefix_count=183700)
    result["execution_seconds"] = time.perf_counter() - started
    result["source_checkpoint_sha256"] = sha(SOURCE)
    result["mathematical_entry_receipt_sha256"] = sha(REDUCTION_ENTERED)
    write_new(IDENTITY, result)
    return dict(passed=True, source_sha256=sha(SOURCE), identity_sha256=sha(IDENTITY))


def write_report():
    parent, nodes, modified, finite, expected = setup()
    source, identity = read_checkpoint(SOURCE), read_checkpoint(IDENTITY)
    require(source == expected and identity["passed"] is True, "new checkpoint mismatch")
    report = copy.deepcopy(parent)
    report.update(stage="current_Rm_actual_pressure_inlet_and_relative_patch_source",
        source_family=copy.deepcopy(source["source_family"]),
        parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        current_Rm_pressure_inlet_source_inputs=source, current_Rm_pressure_inlet_identity=identity,
        execution_seconds=dict(new_actual_Rm_pressure_inlet_identity=identity["execution_seconds"]),
        checkpoint_hashes={path.name: sha(path) for path in (SOURCE, IDENTITY)},
        pre_math_import_failure=read_checkpoint(IMPORT_FAILURE) if IMPORT_FAILURE.exists() else None,
        pre_reduction_layout_failure=read_checkpoint(LAYOUT_FAILURE) if LAYOUT_FAILURE.exists() else None,
        actual_reduction_entry=read_checkpoint(REDUCTION_ENTERED),
        old_212_terms_and_all_errors_unchanged=True, reserved_n1_own_swirl_registration_count=1,
        current_Rm_RQy_germ_identity_proved=False, current_Rm_to_Rp_join=False,
        current_patch_increment_numerically_evaluated=False)
    report["input_hashes"] = dict(parent["input_hashes"], **{PARENT: PARENT_HASH, FIXTURE: FIXTURE_HASH})
    report["implementation_hashes"] = dict(parent["implementation_hashes"], **{KERNEL: KERNEL_HASH, RUNNER: sha(ROOT / RUNNER)})
    report["scope"] = dict(parent["scope"], current_actual_Rm_pressure_inlet_identity=True,
        current_Rm_RQy_germ_identity_proved=False, current_Rm_to_Rp_join=False,
        current_patch_increment_numerically_evaluated=False, global_five_moment_closure=False, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, modified, finite, expected = setup()
    report = read_checkpoint(REPORT)
    source, identity = read_checkpoint(SOURCE), read_checkpoint(IDENTITY)
    require(source == expected and report["current_Rm_pressure_inlet_source_inputs"] == source, "source changed")
    require(report["current_Rm_pressure_inlet_identity"] == identity
            and identity["source_checkpoint_sha256"] == sha(SOURCE), "identity checkpoint changed")
    require(identity["mathematical_entry_receipt_sha256"] == sha(REDUCTION_ENTERED)
            and report["actual_reduction_entry"] == read_checkpoint(REDUCTION_ENTERED), "reduction entry receipt changed")
    require(report["pre_math_import_failure"] == (read_checkpoint(IMPORT_FAILURE) if IMPORT_FAILURE.exists() else None)
            and report["pre_reduction_layout_failure"] == (read_checkpoint(LAYOUT_FAILURE) if LAYOUT_FAILURE.exists() else None),
            "pre-reduction failure record changed")
    require(report["checkpoint_hashes"] == {path.name: sha(path) for path in (SOURCE, IDENTITY)}, "checkpoint hashes changed")
    for key in ("actual_Rm_absolute_m15", "actual_H_m15_Rm_transport", "corrected_own_swirl_provider",
                "corrected_own_swirl_rebinding", "corrected_reserved_point_numeric"):
        require(report[key] == parent[key], "inherited packet changed: " + key)
    for name, expected_hash in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected_hash, "report pin changed: " + name)
    require(report["scope"]["full_goal_complete"] is False
            and report["current_Rm_RQy_germ_identity_proved"] is False, "invalid seam/completion claim")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
        source_checkpoint_sha256=sha(SOURCE), identity_checkpoint_sha256=sha(IDENTITY),
        old_mathematics_executed=False, new_actual_operator_replayed=False,
        current_Rm_RQy_germ_identity_proved=False, current_Rm_to_Rp_join=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--preflight", action="store_true")
    group.add_argument("--compute-new", action="store_true")
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        source_preflight()
        print(json.dumps(dict(passed=True, source_preflight_only=True, old_mathematics_executed=False)))
    else:
        action = compute_new if args.compute_new else write_report if args.write_report else inspect_report
        print(json.dumps(action(), indent=2))
