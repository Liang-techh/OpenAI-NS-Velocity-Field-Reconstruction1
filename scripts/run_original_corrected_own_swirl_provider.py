"""Install a corrected C1 own-swirl function source and cached point provider.

Only --compute-new traverses the new dependency replacement. Packaging and
inspection read checkpoints. No accepted numerical constructor is imported.
"""
from __future__ import annotations

import argparse
import ast
import copy
from fractions import Fraction
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PREFIX = "lei_ren_part1_paper_compliant_current_original_"
PARENT = PREFIX + "corrected_reserved_m15_source.json"
PARENT_HASH = "fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5"
FULL = PREFIX + "n1_full_defect_linear_functions.json.gz"
PULSE_BRIDGE = PREFIX + "n1_pulse_numeric_bridge.json.gz"
FULL_BRIDGE = PREFIX + "n1_full_five_numeric_bridge.json.gz"
LEDGER = PREFIX + "macro_endpoint_m15_numeric.json"
REPAIRED = PREFIX + "repaired_Rp_m15_numeric.json"
KERNEL = "src/openai_ns_reconstruction/original_corrected_own_swirl_provider.py"
KERNEL_HASH = "dab9e9975a6e381b286b8b05b4850400d74540ae1c71815162162ce88dc92e6b"
FIXTURE = "experiments/root_st073/original_corrected_own_swirl_provider_fixture.json"
FIXTURE_HASH = "91b61cc485259c6c100ad6b9b5681628e0be194dc9f16a68f57c08126953a448"
RUNNER = "scripts/run_original_corrected_own_swirl_provider.py"
REPORT = BASE / (PREFIX + "corrected_own_swirl_provider.json")
CHECKPOINTS = ROOT / "tmp/original_corrected_own_swirl_provider"
SOURCE = CHECKPOINTS / "new_source.json"
REBOUND = CHECKPOINTS / "new_rebound.json"
PROVIDER = CHECKPOINTS / "new_provider.json"
ORDERS = ("m11", "m12", "m13", "m14", "m15")


def require(value, message):
    if not value:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                   ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
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
        tree = ast.parse(raw, filename=relative)
        names = [n.name for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))]
        require(len(names) == len(set(names)), "duplicate definitions: " + relative)
    require(sha(ROOT / KERNEL) == KERNEL_HASH, "unreviewed new kernel")
    fixture = read_pinned(ROOT / FIXTURE, FIXTURE_HASH)
    require(fixture["sourcechecks"]["module_sha256"] == KERNEL_HASH
            and all(fixture["checks"].values()), "new fixture witness changed")


def ledger_receipt(ledger, repaired, expected_hash):
    """Read saved ledger positions; do not evaluate, sum or reduce any term."""
    require(ledger["source_family"] == repaired["source_family"], "ledger source family changed")
    receipts = []
    require(len(ledger["samples"]) == len(repaired["samples"]) == 3, "saved ledger rows changed")
    for index, (sample, old) in enumerate(zip(ledger["samples"], repaired["samples"])):
        require(sample["Z"] == old["Z"], "saved Z order changed")
        raw = sample["raw_m15_C0_C1_terms"]
        reduced = sample["correlated_reduced_m15_C0_C1_terms"]
        require(len(raw) == 213 and len(reduced) == 183, "wrong saved ledger sizes")
        require(raw == old["m15_C0_C1_terms"], "raw saved ledger was changed by endpoint reduction")
        require(raw[212]["name"] == reduced[182]["name"] == "reserved_n1_own_swirl_once",
                "wrong reserved slot")
        require(raw[212] == reduced[182], "reserved source changed during saved reduction")
        require(sum(term["name"] == "reserved_n1_own_swirl_once" for term in raw) == 1,
                "reserved source not registered exactly once")
        receipts.append(dict(saved_Z=sample["Z"], sample_index=index,
            raw_forward_212_digest=digest(raw[:212]), raw_all_213_digest=digest(raw),
            reduced_all_183_digest=digest(reduced), six_structural_channels_digest=digest(raw[:6]),
            six_structural_channel_names=[term["name"] for term in raw[:6]],
            old_reserved_digest=digest(raw[212]),
            raw_forward_ref=f"samples[{index}].raw_m15_C0_C1_terms[0:212]",
            raw_reserved_ref=f"samples[{index}].raw_m15_C0_C1_terms[212]",
            reduced_reserved_ref=f"samples[{index}].correlated_reduced_m15_C0_C1_terms[182]",
            old_reserved_record=copy.deepcopy(raw[212]),
            all_212_forward_and_error_sectors_retained=True,
            old_numeric_rows_are_provenance_only=True))
    return dict(parent_report=LEDGER, parent_sha256=expected_hash,
                raw_count=213, forward_count=212, reduced_count=183,
                raw_reserved_index=212, reduced_reserved_index=182,
                six_channels_are_five_signed_plus_core_swirl_error_only=True,
                samples=receipts, numerical_ledger_not_reevaluated=True)


def setup():
    source_preflight()
    parent = read_pinned(BASE / PARENT, PARENT_HASH)
    for name, expected in {**parent["input_hashes"], **parent["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected, "changed inherited pin: " + name)
    full = read_pinned(BASE / FULL, parent["input_hashes"][FULL])
    pulse = read_pinned(BASE / PULSE_BRIDGE, parent["input_hashes"][PULSE_BRIDGE])
    bridge = read_pinned(BASE / FULL_BRIDGE, parent["input_hashes"][FULL_BRIDGE])
    ledger = read_pinned(BASE / LEDGER, parent["input_hashes"][LEDGER])
    repaired = read_pinned(BASE / REPAIRED, parent["input_hashes"][REPAIRED])
    nodes = copy.deepcopy(full["exact_graph_nodes"])
    rows = pulse["normalization_graph_nodes"]
    require(len(nodes) == 53398 and pulse["normalization_graph_span"] == [35558, 53438], "wrong normalization span")
    require(nodes[35558:] == rows[:17840] and len(rows) == 17880, "normalization overlap changed")
    nodes.extend(copy.deepcopy(rows[17840:]))
    require(bridge["normalization_graph_span"] == [53438, 53444]
            and len(bridge["normalization_graph_nodes"]) == 6, "wrong full normalization span")
    nodes.extend(copy.deepcopy(bridge["normalization_graph_nodes"]))
    corrected = parent["corrected_reserved_m15_source"]
    normalized = parent["corrected_reserved_normalized_source"]
    require(len(nodes) == 53444, "wrong saved complete prefix")
    nodes.extend(copy.deepcopy(corrected["delta_nodes"]))
    require(digest(nodes) == corrected["graph_patch"]["full_canonical_nodes_sha256"], "corrected saved graph changed")
    nodes.extend(copy.deepcopy(normalized["appended_nodes"]))
    require(len(nodes) == 53500 and digest(nodes) == normalized["full_nodes_digest"], "normalized saved graph changed")
    require(normalized["common_unit_node"] == 53427 and normalized["unit_Z1_exact_zero"] is True,
            "common positive frozen unit changed")
    old = full["full_defect_linear_bindings"]
    require(old["moments"]["m15"] == [53175, 53176], "old full m15 roots changed")
    require(corrected["physical_C1_new_roots"]["corrected_full_m15"] == [53474, 53475]
            and corrected["physical_C1_new_roots"]["corrected_preRp_m15"] == [53472, 53473], "wrong corrected sources")
    require(normalized["normalized_C1_roots"]["corrected_reserved_m15"] == [53492, 53493]
            and normalized["normalized_C1_roots"]["independent_target_m15"] == [53494, 53495], "wrong normalized sources")
    selected = {key: copy.deepcopy(old[key]) for key in
                ("coefficients", "linear_response_residuals", "repair_velocity")}
    selected["response_matrix"] = [copy.deepcopy(old["response_matrix"][name]) for name in ORDERS]
    source = dict(schema="corrected_own_swirl_provider_inputs_v1", source_only=True,
        parent_sha256=PARENT_HASH, graph_node_count=len(nodes), graph_digest=digest(nodes),
        graph_reassembled_from_saved_rows_only=True, full_graph_sha256=parent["input_hashes"][FULL],
        old_bindings=selected, response_matrix_row_order=list(ORDERS),
        replacements={"53175": 53474, "53176": 53475},
        old_physical_moments=copy.deepcopy(old["moments"]), old_segments=copy.deepcopy(old["segments"]),
        corrected_physical_roots=copy.deepcopy(corrected["physical_C1_new_roots"]),
        corrected_normalized_roots=copy.deepcopy(normalized["normalized_C1_roots"]),
        unit_receipt=copy.deepcopy(parent["corrected_terminal_source_inputs"]["unit_receipt"]),
        source_family=copy.deepcopy(full["source_family"]),
        ledger=ledger_receipt(ledger, repaired, parent["input_hashes"][LEDGER]),
        cached_corrected_reserved=copy.deepcopy(parent["corrected_reserved_point_numeric"]["rawCorrectedReserved"]),
        cached_point_source_ref="parent.corrected_reserved_point_numeric.rawCorrectedReserved",
        cached_point_exact_Z="371/1000", accepted_ancestors_replayed=False)
    return parent, nodes, source


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write((json.dumps(value, indent=2, ensure_ascii=True, allow_nan=False) + "\n").encode())


def read_checkpoint(path):
    return json.loads(path.read_bytes())


def load_kernel():
    spec = importlib.util.spec_from_file_location("new_corrected_own_swirl", ROOT / KERNEL)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def install_provider(module, source, rebound):
    """Copy existing source roots and numeric tuples; perform no new quadrature."""
    require(rebound["passed"] is True, "new rebound checkpoint failed")
    updated = rebound["updated_bindings"]
    old = source["old_bindings"]
    for name in ("alpha1", "alpha2"):
        require(updated["coefficients"][name] == old["coefficients"][name], "axial coefficient changed")
    for name in ("beta1", "beta2", "beta3"):
        require(updated["coefficients"][name] != old["coefficients"][name], "angular coefficient was not rebound")
    require(updated["response_matrix"] == old["response_matrix"], "response matrix changed")
    for name in ("Uz1", "physical_radius", "log_R_over_Rc"):
        require(updated["repair_velocity"][name] == old["repair_velocity"][name], "independent field port changed")
    for name in ("Utheta1", "F1"):
        require(updated["repair_velocity"][name] != old["repair_velocity"][name], "swirl field was not rebound")
    moments = copy.deepcopy(source["old_physical_moments"])
    moments["m15"] = source["corrected_physical_roots"]["corrected_full_m15"]
    segments = copy.deepcopy(source["old_segments"])
    segments["preRp"]["m15"] = source["corrected_physical_roots"]["corrected_preRp_m15"]
    bindings = dict(field_source_bindings=dict(repair_velocity=updated["repair_velocity"], moments=moments,
                    segments=segments, graph_digest=rebound["new_graph_digest"], ordinary_Z_orders=[0, 1]),
        coefficient_source_bindings=dict(coefficients=updated["coefficients"],
            response_matrix=dict(zip(ORDERS, updated["response_matrix"])),
            linear_response_residuals=updated["linear_response_residuals"], ordinary_Z_orders=[0, 1]))
    leaf = dict(name="reserved_n1_own_swirl_once", physical_measure="-R*theta0*Utheta1 ds",
        physical_C1_roots=source["corrected_physical_roots"]["corrected_reserved_m15"],
        normalized_C1_roots=source["corrected_normalized_roots"]["corrected_reserved_m15"],
        physical_full_C1_roots=moments["m15"],
        normalized_full_C1_roots=source["corrected_normalized_roots"]["corrected_full_m15"],
        unit_node=53427, unit_Z1_exact_zero=True, ordinary_Z_orders=[0, 1], full_support_x=2,
        exact_angular_response=[0, 0, 1], source_identity="reserved_corrected=-full_corrected",
        coefficient_and_field_bindings=bindings, own_swirl_registration_count=1,
        current_forward_germ_join_admitted=False, legacy_terminal_handles_are_provenance_only=True)
    def cached_record(row):
        return dict(source_enclosure_exact_mpf_tuples=[row["lower_exact_mpf_tuple"], row["upper_exact_mpf_tuple"]],
                    directed=True, bounds_are_function_values=False, cached_parent_enclosure=True)
    packet = dict(Z="371/1000", native_compact_x="2", source_family=source["source_family"],
        coefficient_C0_C1=[cached_record(row) for row in source["cached_corrected_reserved"]],
        angular_inverse_response_exact=[0, 0, 1], original_source_bindings=leaf,
        contract=dict(ordinary_Z1_retained=True, fixed_reference_A_Rv_squared_unit_once=True,
            full_support_response_is_exact_matrix_row_identity=True,
            independent_complete_pulse_not_used_to_define_this_response=True,
            numeric_packet_copied_without_recomputation=True),
        limits=dict(point_only=True, exact_Z="371/1000", numerical_Utheta1_or_F1=False,
                    current_forward_germ_join=False, all_Z_functional_closure=False))
    provider = module.CorrectedOwnSwirlProvider(bindings, packet)
    output = provider.evaluate(Fraction(371, 1000), 2)
    require(output == packet and provider.field_source_bindings() == bindings["field_source_bindings"]
            and provider.coefficient_source_bindings() == bindings["coefficient_source_bindings"], "data-only provider changed packet")
    indices = [index for index in range(183) if index not in (38, 39, 40)]
    overlay = [dict(reduced_index=index, source_override=leaf) if index == 182 else
               dict(reduced_index=index, unchanged_saved_source_reference=True) for index in indices]
    return dict(passed=True, source_only_fields=True, binding_packet=bindings, response_source_leaf=leaf,
        cached_point_packet=output, cached_point_equal_to_parent_exact_tuples=True,
        raw_ledger_override=dict(index=212, replacement=leaf, preceding_212_unchanged=True),
        reduced_ledger_override=dict(index=182, replacement=leaf, all_other_182_unchanged=True),
        joint_source_overlay=overlay, joint_source_count=180, unchanged_joint_source_count=179,
        source_slot_replaced_once=True, extra_own_swirl_added=False, old_numeric_rows_transferred=False,
        numerical_ledger_not_reevaluated=True, legacy_default_provider_chain_reconfigured=False,
        full_forward_chain_uses_new_provider=False, numerical_Utheta1_F1_evaluated=False,
        repaired_own_M_V_P_C3_rebuilt=False, current_Rm_to_Rp_join=False, global_closure=False,
        source_checkpoint_sha256=sha(SOURCE), rebound_checkpoint_sha256=sha(REBOUND))


def compute_new():
    parent, nodes, source = setup()
    require(not any(path.exists() for path in (SOURCE, REBOUND, PROVIDER)), "checkpoints exist; do not replay")
    write_new(SOURCE, source)
    module = load_kernel()
    started = time.perf_counter()
    result = module.rebind(nodes, source["old_bindings"], {53175: 53474, 53176: 53475},
                           expected_nodes_digest=source["graph_digest"])
    result["execution_seconds"] = time.perf_counter() - started
    result["source_checkpoint_sha256"] = sha(SOURCE)
    write_new(REBOUND, result)
    result = read_checkpoint(REBOUND)
    started = time.perf_counter()
    provider = install_provider(module, source, result)
    provider["execution_seconds"] = time.perf_counter() - started
    write_new(PROVIDER, provider)
    return dict(passed=True, source_sha256=sha(SOURCE), rebound_sha256=sha(REBOUND), provider_sha256=sha(PROVIDER))


def write_report():
    parent, nodes, expected = setup()
    source, rebound, provider = (read_checkpoint(path) for path in (SOURCE, REBOUND, PROVIDER))
    require(source == expected and rebound["passed"] and provider["passed"], "new checkpoint mismatch")
    report = copy.deepcopy(parent)
    report.update(stage="corrected_own_swirl_C1_function_and_cached_point_provider",
        parent_report=PARENT, parent_report_sha256=PARENT_HASH,
        corrected_own_swirl_source_inputs=source, corrected_own_swirl_rebinding=rebound,
        corrected_own_swirl_provider=provider,
        execution_seconds=dict(new_C1_dependency_rebinding=rebound["execution_seconds"],
                               data_only_provider_installation=provider["execution_seconds"]),
        checkpoint_hashes={path.name: sha(path) for path in (SOURCE, REBOUND, PROVIDER)},
        corrected_C1_field_source_provider_installed=True,
        corrected_numeric_m15_response_provider_installed_at_cached_point=True,
        legacy_provider_files_unchanged=True, full_forward_chain_uses_new_provider=False,
        existing_velocity_provider_rebound=False, old_212_terms_and_all_errors_unchanged=True,
        reserved_n1_own_swirl_registration_count=1, extra_reserved_term_added=False)
    report["input_hashes"] = dict(parent["input_hashes"], **{PARENT: PARENT_HASH, FIXTURE: FIXTURE_HASH})
    report["implementation_hashes"] = dict(parent["implementation_hashes"], **{KERNEL: KERNEL_HASH, RUNNER: sha(ROOT / RUNNER)})
    report["scope"] = dict(parent["scope"], corrected_C1_field_source_provider=True,
        corrected_numeric_m15_cached_point_provider=True, full_forward_chain_uses_new_provider=False,
        numerical_Utheta1_F1_evaluated=False, repaired_own_M_V_P_C3_rebuilt=False,
        current_Rm_to_Rp_join=False, independent_terminal_function_match=False,
        global_five_moment_closure=False, full_goal_complete=False)
    write_new(REPORT, report)
    return dict(passed=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT))


def inspect_report():
    parent, nodes, expected = setup()
    report = read_checkpoint(REPORT)
    source, rebound, provider = (read_checkpoint(path) for path in (SOURCE, REBOUND, PROVIDER))
    require(source == expected and report["corrected_own_swirl_source_inputs"] == source, "source changed")
    require(report["corrected_own_swirl_rebinding"] == rebound and report["corrected_own_swirl_provider"] == provider,
            "report differs from checkpoint")
    require(report["checkpoint_hashes"] == {path.name: sha(path) for path in (SOURCE, REBOUND, PROVIDER)}, "checkpoint link changed")
    require(rebound["source_checkpoint_sha256"] == sha(SOURCE)
            and provider["source_checkpoint_sha256"] == sha(SOURCE)
            and provider["rebound_checkpoint_sha256"] == sha(REBOUND), "component link changed")
    require(digest(nodes + rebound["delta_nodes"]) == rebound["new_graph_digest"], "new graph digest changed")
    require(len(nodes) + len(rebound["delta_nodes"]) == rebound["new_graph_node_count"], "new graph count changed")
    require(provider["source_slot_replaced_once"] and not provider["extra_own_swirl_added"], "wrong reserved count")
    require([row["reduced_index"] for row in provider["joint_source_overlay"]]
            == [index for index in range(183) if index not in (38, 39, 40)], "joint indices changed")
    require(sum("source_override" in row for row in provider["joint_source_overlay"]) == 1, "multiple source overrides")
    for row, cached in zip(provider["cached_point_packet"]["coefficient_C0_C1"], source["cached_corrected_reserved"]):
        require(row["source_enclosure_exact_mpf_tuples"] ==
                [cached["lower_exact_mpf_tuple"], cached["upper_exact_mpf_tuple"]], "cached tuples changed")
    for key in ("corrected_reserved_m15_source", "corrected_reserved_normalized_source", "corrected_reserved_point_numeric",
                "actual_Rm_absolute_m15", "actual_H_m15_Rm_transport"):
        require(report[key] == parent[key], "inherited packet changed: " + key)
    for name, expected_hash in {**report["input_hashes"], **report["implementation_hashes"]}.items():
        require(sha(path_for(name)) == expected_hash, "report pin changed: " + name)
    require(report["scope"]["full_goal_complete"] is False, "invalid completion claim")
    return dict(passed=True, serialized_only=True, bytes=REPORT.stat().st_size, sha256=sha(REPORT),
        new_graph_node_count=rebound["new_graph_node_count"], appended_dependency_nodes=len(rebound["delta_nodes"]),
        all_212_forward_terms_retained=True, reserved_source_overrides=1,
        ordinary_Z_orders=[0, 1], unchanged_alpha_Uz_response=True,
        old_mathematics_executed=False, actual_rebinding_replayed=False, current_Rm_to_Rp_join=False)


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
