"""Compare actual leading loop roots and first native jets, ordinary Z0/Z1."""
from __future__ import annotations

import argparse
import ast
import copy
import importlib.util
from pathlib import Path
import sys
import time
import json

# These helpers only authenticate/read/write saved bytes. No operator runs
# when importing the preceding runner.
from run_original_Rm_second_radial_identity import (
    ROOT, BASE, PREFIX, require, sha, digest, path_for, read_pinned, read_saved, write_new,
)

PARENT = PREFIX + "Rm_second_radial_identity.json"
PARENT_HASH = "f221e697094a0b290b0674b11bcca534c688e0982ca775ffe923a42c839afab9"
FINITE = PREFIX + "true_Rm_patch_finite_N_C5_functions.json.gz"
Q4 = PREFIX + "modified_Q4_pressure_C3_functions.json.gz"
LOOP = PREFIX + "prescribed_loop_C5_functions.json.gz"
KERNEL = "src/openai_ns_reconstruction/original_Rm_generic_native_identity.py"
KERNEL_HASH = "91dd7abc32fe362538626fd2aacf870655608bbeaa3e820dc67a984ea6dff193"
FIXTURE = "experiments/root_st073/original_Rm_generic_native_identity_fixture.json"
FIXTURE_HASH = "95505acface4fbbccac52a72ad5bca743ed7e15912915726a2fc72f4c0b5cfe8"
RUNNER = "scripts/run_original_Rm_generic_native_identity.py"
REPORT = BASE / (PREFIX + "Rm_generic_native_identity.json")
CHECKPOINTS = ROOT / "tmp/original_Rm_generic_native_identity"
SOURCE = CHECKPOINTS / "new_source.json"
ENTERED = CHECKPOINTS / "new_operator_entered.json"
RESULT = CHECKPOINTS / "new_generic_native_identity.json"
RADIUS = CHECKPOINTS / "new_radius_source_layout.json"
DIAGNOSTIC = "experiments/root_st073/original_Rm_radius_probe_diagnostic.json"
DIAGNOSTIC_HASH = "fc946e9aad9caa5056e84fbc6c4dc9d6253fa89f76175f9fa0df4e55bb3ada97"
HISTORY_REPORT = PREFIX + "Rm_history_native_source.json"
HISTORY_REPORT_HASH = "c3c10f2c622d6079805723ae2b1e168f36621e45812ca89676e7816e8d0acd1a"


def saved_history_report(parent, nodes):
    history = read_pinned(BASE / HISTORY_REPORT, HISTORY_REPORT_HASH)
    require(history["parent_report"] == PARENT and history["parent_report_sha256"] == PARENT_HASH,
            "independent history source has another parent")
    source = history["actual_Rm_history_native_source"]
    require(source["passed"] is True and source["graph_delta"]["base_nodes_sha256"] == digest(nodes)
            and source["remaining_density_seam_unproved"] is True, "history source not admitted on this prefix")
    for key in ("actual_Rm_second_radial_identity", "actual_Rm_second_radial_endpoint_functions",
                "current_Rm_pressure_inlet_identity", "actual_Rm_absolute_m15", "corrected_own_swirl_provider"):
        require(history[key] == parent[key], "independent history parent packet changed: "+key)
    return history


def radius_source_lemma(nodes):
    """Record an analytic lemma from exact source constructors, without AD.

    R(p)=110*exp(L-7+p), where node L precedes p in the pinned backward DAG.
    The chain rule gives R_p=R. For x=R/Rm, R=Rm*x and D_y=x*d/dx.
    Closed endpoint substitutions have derivative zero. No symbolic radius
    helper, value approximation or captured diagnostic is used as a premise.
    """
    require(nodes[84533] == dict(operation="product", arguments=[79968, 84532]), "post radius constructor changed")
    require(nodes[79968] == dict(operation="exact_rational", numerator=110, denominator=1), "post radius scale changed")
    require(nodes[84532] == dict(operation="analytic_unary", name="exp", argument=84531), "post radius exponential changed")
    require(nodes[84531] == dict(operation="sum", arguments=[81514, 84530])
            and nodes[84530] == dict(operation="sum", arguments=[84529, 84528])
            and nodes[84529] == dict(operation="exact_rational", numerator=-7, denominator=1), "post radius affine argument changed")
    require(81514 < 84528 and nodes[84528]["operation"] == "bound_variable", "native independence guard changed")
    for endpoint, body, variable in ((85743, 84533, 84528), (187963, 187214, 185676)):
        row = nodes[endpoint]
        require(row["operation"] == "function_substitution" and row["expression"] == body
                and row["variable"] == variable and row["value"] == 1
                and row["Z_independent_substitution"] is True, "radius endpoint definition changed")
    require(nodes[187214] == dict(operation="product", arguments=[85743, 185676])
            and nodes[187215]["operation"] == "positive_quotient" and nodes[187215]["numerator"] == 1
            and nodes[187215]["denominator"] == 185676, "patch radius or Jacobian changed")
    ids = (79968, 81514, 84528, 84529, 84530, 84531, 84532, 84533, 85743, 185676, 187214, 187215, 187963)
    return dict(schema="Rm_analytical_radius_source_layout_lemma_v1", verified=True, status="proved",
        proof_kind="analytical_source_layout_lemma", source_only=True,
        exact_source_descriptors={str(index):copy.deepcopy(nodes[index]) for index in ids},
        postrestore_formula="R(p)=110*exp(L-7+p); L=node81514 is independent of p=node84528 by backward DAG",
        patch_formula="R(x)=Rm*x; D_y=x*d/dx; closed Rm has native derivative zero",
        proof_rule="elementary chain and product rules on the displayed exact source definitions",
        automatic_radius_differentiation_replayed=False, diagnostic_probe_used_as_premise=False,
        proof=dict(postrestore=dict(body_node=84533, endpoint_node=85743, coordinate_node=84528, identity_text="D_y R - R = 0"),
                   patch=dict(body_node=187214, endpoint_node=187963, coordinate_node=185676, physical_measure_node=187215,
                              identity_text="R_patch = Rm*x, J = 1/x")),
        proof_ref=dict(kind="analytical_source_definition_lemma", parent_report_sha256=PARENT_HASH,
                       exact_descriptors_digest=digest({str(index):nodes[index] for index in ids})))


def preflight():
    for relative in (RUNNER, KERNEL):
        raw = (ROOT / relative).read_bytes()
        require(b"\r" not in raw and raw.endswith(b"\n") and not raw.endswith(b"\n\n"), "LF/EOF: "+relative)
        ast.parse(raw, filename=relative)
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "new generic consumer changed")
    fixture = read_pinned(ROOT / FIXTURE, FIXTURE_HASH)
    require(fixture.get("source_only") is True and isinstance(fixture.get("nodes"), list)
            and fixture.get("root_data_only_execution_audit", {}).get("fixture_not_used_as_source_identity_premise") is True,
            "preserved generic input fixture changed")


def setup():
    preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: "+name)
    finite, q4, loop = (read_pinned(BASE / name, parent["input_hashes"][name]) for name in (FINITE, Q4, LOOP))
    nodes = finite["exact_graph_nodes"]
    first = parent["actual_Rm_native_radial_endpoint_source"]
    second = parent["actual_Rm_second_radial_endpoint_functions"]
    require(len(nodes) == 205147 and nodes[:183700] == q4["exact_graph_nodes"]
            and nodes[:121834] == loop["exact_graph_nodes"]
            and digest(nodes) == first["base_digest"], "saved source prefix changed")
    nodes = nodes + copy.deepcopy(first["appended_nodes"])
    require(digest(nodes) == first["full_nodes_digest"] == second["base_digest"], "first endpoint delta changed")
    nodes = nodes + copy.deepcopy(second["appended_nodes"])
    require(len(nodes) == 205160 and digest(nodes) == second["full_nodes_digest"], "second endpoint delta changed")
    seed = parent["actual_Rm_leading_seed_admission"]
    raw = parent["actual_Rm_raw_native_jet_identity"]
    proof = parent["actual_Rm_second_radial_identity"]
    require(seed["passed"] is True and raw["passed"] is True and proof["passed"] is True,
            "actual leading source premises incomplete")
    mb, fb = q4["modified_Q4_pressure_C3_bindings"], finite["Rm_patch_finite_N_C5_bindings"]
    lemma = radius_source_lemma(nodes)
    mb, fb = copy.deepcopy(mb), copy.deepcopy(fb)
    leading_post = loop["prescribed_loop_C5_bindings"]["charts"]["postrestore"]
    for key in ("coordinate", "radius", "physical_measure", "roots"):
        mb["charts"]["postrestore"][key] = copy.deepcopy(leading_post[key])
    mb["radius_source"] = copy.deepcopy(lemma)
    fb["radius_source"] = copy.deepcopy(lemma)
    require(sha(ROOT / DIAGNOSTIC) == DIAGNOSTIC_HASH, "original probe diagnostic changed")
    packet = dict(schema="actual_Rm_generic_native_inputs_v1", source_only=True,
        parent_report=PARENT, parent_sha256=PARENT_HASH, graph_node_count=len(nodes), graph_digest=digest(nodes),
        ordinary_Z_orders=[0, 1], leading_roots_only=True, modified_p2_excluded=True,
        postrestore_leading_roots_archive=LOOP, postrestore_leading_roots_archive_sha256=parent["input_hashes"][LOOP],
        left_roots=copy.deepcopy(mb["charts"]["postrestore"]["roots"]), right_roots=copy.deepcopy(fb["chart"]["roots"]),
        seed_admission_digest=digest(seed), raw_native_jet_digest=digest(raw), second_radial_identity_digest=digest(proof),
        radius_source_lemma=lemma, radius_source_lemma_digest=digest(lemma),
        auxiliary_probe_diagnostic=dict(path=DIAGNOSTIC, sha256=DIAGNOSTIC_HASH, admitted_as_premise=False,
            successful_helper_calls=2, duplicate_helper_calls=1, execution_module_hash_available=False,
            checkpoint_saved_from_original_stdout_after_inspection=True),
        no_old_operator_replayed=True, all_212_terms_and_errors_unchanged=True)
    return parent, nodes, mb, fb, packet


def prepare_source():
    *_, packet = setup()
    require(not any(path.exists() for path in (SOURCE, ENTERED, RESULT, RADIUS)), "source already prepared")
    write_new(RADIUS, packet["radius_source_lemma"])
    write_new(SOURCE, packet)
    return dict(passed=True, source_sha256=sha(SOURCE), graph_node_count=packet["graph_node_count"])


def compute_new():
    parent, nodes, mb, fb, packet = setup()
    require(read_saved(SOURCE) == packet and not ENTERED.exists() and not RESULT.exists(), "source changed or entered; no replay")
    spec = importlib.util.spec_from_file_location("new_Rm_generic_native", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    write_new(ENTERED, dict(source_sha256=sha(SOURCE), kernel_sha256=KERNEL_HASH))
    started = time.perf_counter()
    result = module.prove(nodes, {"modified_Q4_pressure_C3_bindings":mb}, {"Rm_patch_finite_N_C5_bindings":fb},
        parent["actual_Rm_leading_seed_admission"], parent["actual_Rm_raw_native_jet_identity"],
        parent["actual_Rm_second_radial_identity"], expected_nodes_digest=packet["graph_digest"])
    result.update(execution_seconds=time.perf_counter()-started, source_checkpoint_sha256=sha(SOURCE),
        entry_checkpoint_sha256=sha(ENTERED), actual_saved_graph_run=True, actual_saved_graph_loaded=True)
    write_new(RESULT, result)
    return dict(passed=result.get("passed"), result_sha256=sha(RESULT), unresolved_gates=result.get("unresolved_gates"))


def write_report():
    parent, nodes, mb, fb, expected = setup()
    packet, result = read_saved(SOURCE), read_saved(RESULT)
    require(packet == expected, "saved source changed")
    history = saved_history_report(parent, nodes)
    report = copy.deepcopy(parent)
    report.update(stage="actual_Rm_leading_generic_native_identity", parent_report=PARENT,
        parent_report_sha256=PARENT_HASH, actual_Rm_generic_native_inputs=packet,
        actual_Rm_generic_native_identity=result,
        generic_native_checkpoint_hashes={path.name:sha(path) for path in (SOURCE, ENTERED, RESULT, RADIUS)})
    report["execution_seconds"]["new_actual_generic_native_identity"] = result["execution_seconds"]
    report["input_hashes"].update({PARENT:PARENT_HASH, FIXTURE:FIXTURE_HASH, DIAGNOSTIC:DIAGNOSTIC_HASH})
    report["implementation_hashes"].update({KERNEL:KERNEL_HASH, RUNNER:sha(ROOT / RUNNER)})
    for kind in ("input_hashes", "implementation_hashes"):
        for name, value in history[kind].items():
            require(report[kind].get(name, value) == value, "independent history pin conflict: "+name)
            report[kind][name] = value
    report["input_hashes"][HISTORY_REPORT] = HISTORY_REPORT_HASH
    for key in ("actual_Rm_history_native_inputs", "actual_Rm_history_native_source", "history_native_checkpoint_hashes"):
        report[key] = copy.deepcopy(history[key])
    report["execution_seconds"]["new_actual_history_native_source"] = history["execution_seconds"]["new_actual_history_native_source"]
    report["independent_history_component"] = dict(report=HISTORY_REPORT, sha256=HISTORY_REPORT_HASH,
        common_source_prefix_node_count=len(nodes), common_source_prefix_digest=digest(nodes),
        remaining_density_seam_unproved=True)
    report["current_turn_execution_audit"] = dict(primary_current_source_operators_run_once_each=True,
        primary_components=["second_radial", "leading_generic_native", "actual_history_FTC_source"],
        accepted_ancestor_builders_checkers_providers_quadratures_replayed=False,
        auxiliary_radius_helper_successful_calls=2, auxiliary_radius_duplicate_calls=1,
        auxiliary_probe_execution_module_hash_captured=False, auxiliary_probe_saved_after_stdout_inspection=True,
        auxiliary_probe_is_diagnostic_only=True, radius_premise="analytical_source_layout_lemma",
        generic_positive_fixture_execution_receipt_available=False,
        old_fixture_results_not_rebound_to_new_source_hashes=True,
        raw_generic_no_saved_graph_execution_flag_is_template_metadata=True,
        actual_current_generic_source_traversal_executed=True,
        full_goal_complete=False)
    report["scope"].update(actual_leading_generic_native_identity=result["passed"], generic_ordinary_Z_orders=[0, 1],
        modified_finite_N_germ=False, current_Rm_RQy_germ_identity_proved=False, current_Rm_to_Rp_join=False,
        global_five_moment_closure=False, full_goal_complete=False)
    report["scope"].update(actual_current_Rm_history_endpoint_identity=True, actual_current_Rm_history_native_source=True,
                           remaining_density_seam_unproved=True)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, mb, fb, expected = setup()
    report, result = read_saved(REPORT), read_saved(RESULT)
    history = saved_history_report(parent, nodes)
    require(report["actual_Rm_generic_native_inputs"] == read_saved(SOURCE) == expected
            and report["actual_Rm_generic_native_identity"] == result, "saved generic component changed")
    require(read_saved(RADIUS) == expected["radius_source_lemma"], "analytic radius source lemma changed")
    for key in ("actual_Rm_history_native_inputs", "actual_Rm_history_native_source", "history_native_checkpoint_hashes"):
        require(report[key] == history[key], "saved independent history component changed: "+key)
    require(report["generic_native_checkpoint_hashes"] == {path.name:sha(path) for path in (SOURCE, ENTERED, RESULT, RADIUS)}
            and result["source_checkpoint_sha256"] == sha(SOURCE) and result["entry_checkpoint_sha256"] == sha(ENTERED),
            "generic checkpoint link changed")
    for key in ("actual_Rm_second_radial_identity", "actual_Rm_second_radial_endpoint_functions", "actual_Rm_leading_seed_admission",
                "current_Rm_pressure_inlet_identity", "actual_Rm_absolute_m15", "corrected_own_swirl_provider", "corrected_reserved_point_numeric"):
        require(report[key] == parent[key], "inherited packet changed: "+key)
    for name, expected_hash in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected_hash, "report pin changed: "+name)
    require(report["scope"]["full_goal_complete"] is False and report["current_Rm_RQy_germ_identity_proved"] is False,
            "invalid physical germ/completion claim")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
        old_mathematics_executed=False, new_operator_replayed=False, full_goal_complete=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    for name in ("preflight", "prepare-source", "compute-new", "write-report", "inspect-report"):
        group.add_argument("--"+name, action="store_true")
    args = parser.parse_args()
    if args.preflight:
        preflight()
        print(json.dumps(dict(passed=True, source_preflight_only=True)))
    else:
        action = prepare_source if args.prepare_source else compute_new if args.compute_new else write_report if args.write_report else inspect_report
        print(json.dumps(action(), indent=2))
