"""Produce the signed macro source primitive; inspect the saved packet without replay."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction.original_macro_weighted_m15 import (  # noqa: E402
    MASS,
    PARENT,
    PARENT_SHA256,
    OriginalMacroWeightedM15,
    digest,
    load,
    sha,
)
from openai_ns_reconstruction.original_weighted_leaf_m15 import (  # noqa: E402
    materialize_view,
)


DIRECTORY = ROOT / "experiments/root_st073"
NAME = "lei_ren_part1_paper_compliant_current_original_macro_weighted_m15_numeric.json"
CODE = (
    "src/openai_ns_reconstruction/original_macro_weighted_kernel.py",
    "src/openai_ns_reconstruction/original_macro_weighted_m15.py",
    "scripts/check_original_macro_weighted_kernel.py",
    "scripts/run_original_macro_weighted_m15.py",
)
EXPECTED_Z = ("-.3", "0", ".3")
EXPECTED_RATES = [2, 2, 1, 3, 2, 1]


def require(value, message):
    if not value:
        raise AssertionError(message)


def parent():
    require(sha(DIRECTORY / PARENT) == PARENT_SHA256, "Pinned weighted parent bytes changed")
    report = load(DIRECTORY / PARENT)
    require(report["checks"]["passed"] is True, "Passed weighted parent required")
    return report


def _assert_unselected(value, path="source"):
    """Reject every serialized approximation promoted to an exact source value."""
    if isinstance(value, dict):
        if "approximation_is_exact_source_value" in value:
            require(value["approximation_is_exact_source_value"] is False, f"Promoted approximation: {path}")
        if "source_enclosure_is_not_a_selected_value" in value:
            require(value["source_enclosure_is_not_a_selected_value"] is True, f"Selected enclosure: {path}")
        for key, item in value.items():
            _assert_unselected(item, f"{path}.{key}")
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _assert_unselected(item, f"{path}[{index}]")


def _positive_kernel(row):
    require(isinstance(row, dict), "Malformed radial kernel row")
    _assert_unselected(row, "kernel")
    endpoints = row.get("source_enclosure_exact_mpf_tuples")
    require(isinstance(endpoints, list) and len(endpoints) == 2, "Kernel enclosure missing")
    values = [mp.make_mpf(tuple(endpoint)) for endpoint in endpoints]
    require(all(value > 0 for value in values), "Radial kernel is not strictly positive")
    require(mp.mpf(row["approximation"]) > 0, "Radial kernel point is not positive")


def _source_rows(report):
    rows = report.get("macro_weighted_source_rows")
    require(isinstance(rows, list) and [row.get("Z") for row in rows] == list(EXPECTED_Z), "Source sample order changed")
    for row in rows:
        require(row["ordinary_source_orders"] == [0, 1, 2], "Ordinary source Z2 rows missing")
        require(row["ordinary_coefficient_orders"] == [0, 1], "Ordinary C0/C1 rows missing")
        require(len(row["weighted_response_rows"]) == 3, "Weighted ordinary row count changed")
        require(len(row["weighted_C0_C1"]) == 2, "Weighted C0/C1 row count changed")
        require(len(row["native_lower_endpoint_C0_C1"]) == 2, "Native endpoint C0/C1 row count changed")
        require(len(row["kernel_rows"]) == 6, "Kernel row count changed")
        require(row["mode_rates"] == EXPECTED_RATES, "Mode rate order changed")
        require(row["kernel_positive_weights_and_L_Z_independent"] is True, "Kernel positivity certificate missing")
        require(row["rho_inside_dy_and_external_R0_have_distinct_roles"] is True, "Radius roles changed")
        require(row["uniform_macro_angular_enclosure_not_frozen_endpoint_value"] is True, "Uniform enclosure semantics changed")
        require(row["native_endpoint_is_macro0_endpoint"] is False, "Unproved native/macro endpoint identity selected")
        require(row["native_endpoint_used_in_ledger"] is False, "Native endpoint entered the ledger")
        require(row["old_mass_or_core_computed"] is False, "Old mass/core recomputation selected")
        require(row["physical_integral_scale"] == "R100*epsilon*hbB", "Integral physical scale changed")
        require(row["physical_lower_scale"] == "4*epsilon^2*hbB", "Native endpoint physical scale changed")
        for key in ("weighted_response_rows", "weighted_C0_C1", "native_lower_endpoint_C0_C1"):
            _assert_unselected(row[key], f"{row['Z']}.{key}")
        for kernel in row["kernel_rows"]:
            _positive_kernel(kernel)
    return rows


def _check_source_guards(report, parent_report):
    guards = report["source_guards"]
    require(guards["passed"] is True, "Macro source guards failed")
    require(all(value["passed"] is True for value in guards["AST"].values()), "Macro AST source guard failed")
    require(all(value is True for value in guards["checks"].values()), "Macro source guard check failed")
    require(guards["inherited_actual_FTC_and_seam_digest"] == digest(parent_report["source_guards"]), "Inherited FTC witness changed")
    lower = load(DIRECTORY / MASS)["transported_lower_dy"]["source_bindings"]
    require(guards["inherited_lower_source_digest"] == digest(lower), "Inherited lower source witness changed")
    require(guards["uniform_prefix_proof_count"] == 18, "Uniform prefix proof count changed")
    require(guards["native_lower_is_macro_source_port_proved"] is False, "Native lower source port claimed")
    require(guards["native_lower_and_macro0_value_continuity_proved"] is False, "Native/macro value continuity claimed")
    require(guards["coupled_raw23_raw24_raw44_replacement_applied"] is False, "Unproved coupled replacement applied")
    require(guards["source_caps_or_interval_overlap_used_as_function_identity"] is False, "Interval overlap promoted to identity")


def _ledger_samples(parent_report):
    packet = parent_report["inherited_m15_packet"]
    ledger_name = packet["report"]
    ledger_path = DIRECTORY / ledger_name
    require(sha(ledger_path) == packet["report_sha256"], "Pinned original ledger bytes changed")
    require(parent_report["input_hashes"].get(ledger_name) == packet["report_sha256"], "Ledger input pin changed")
    ledger = load(ledger_path)
    by_z = {sample["Z"]: sample for sample in ledger["samples"]}
    require(tuple(by_z) == EXPECTED_Z, "Original ledger sample order changed")
    return packet, by_z, ledger["source_seam_proofs"]["joins"]["micro_macro"]


def _check_view(view, parent_view, sample):
    raw = sample["raw_m15_C0_C1_terms"]
    reduced = sample["correlated_reduced_m15_C0_C1_terms"]
    require(len(raw) == 213 and len(reduced) == 183, "Original archive size changed")
    require(view == parent_view, "Signed 183 view changed from parent")
    require(view["parent_report"] == "lei_ren_part1_paper_compliant_current_original_macro_endpoint_m15_numeric.json", "View parent changed")
    require(view["raw213_digest"] == digest(raw), "Raw 213 archive digest changed")
    require(view["archived_reduced183_digest"] == digest(reduced), "Reduced 183 archive digest changed")
    require(view["parent_reduced_terms_digest"] == digest(reduced), "Parent reduced archive digest changed")
    require(view["view_term_count"] == len(view["entries"]) == 183, "Signed view term count changed")
    require(view["unchanged_parent_references"] == 181, "Signed view reference count changed")
    require([entry["parent_reduced_index"] for entry in view["entries"]] == list(range(183)), "Parent view indices changed")

    replacements = [entry for entry in view["entries"] if entry["operation"] != "unchanged_parent_reference"]
    require([(entry["parent_reduced_index"], entry["raw_index"]) for entry in replacements] == [(41, 46), (42, 47)], "Signed replacement scope changed")
    require(view["entries"][19] == parent_view["entries"][19] and view["entries"][20] == parent_view["entries"][20], "Original physical atoms 19/20 changed")
    require(view["entries"][39] == parent_view["entries"][39], "Original index 39 changed")
    require(view["entries"][41] == parent_view["entries"][41] and view["entries"][42] == parent_view["entries"][42], "Signed 41/42 parent entries changed")
    require(all(view["entries"][index] == parent_view["entries"][index] for index in list(range(183)) if index not in (41, 42)), "Unchanged original references changed")

    expanded = materialize_view(reduced, view["entries"])
    require(digest(expanded) == view["expanded_signed_view_digest"], "Signed view expansion changed")
    require(expanded[:6] == reduced[:6], "Six inner/lower terms changed")
    require(view["six_inner_lower_digest"] == digest(raw[:6]), "Six inner/lower witness changed")
    require(sum(term["name"] == "reserved_n1_own_swirl_once" for term in expanded) == 1, "Reserved term count changed")
    require(view["independent_target"] == sample["independent_target"], "Independent target changed")
    require(view["independent_target_digest"] == digest(sample["independent_target"]), "Independent target digest changed")
    binding = view["independent_target"]["actual_function_binding"]
    require(binding["sign"] == -1 and binding["physical_measure"] == "-R*pressure_density*dR", "Independent target sign/measure changed")
    require(view["selected_N_evidence_digest"] == digest(sample["selected_N_evidence"]), "Selected N witness changed")
    require(view["largest_unchanged_caps"] == sample["largest_remaining_individual_log_caps"], "Original error caps changed")
    require(view["closure_or_small_residual_claim"] is False, "Closure claim introduced")

    # Raw 23/24/43/44/45 and every other non-replaced archive term remain exact.
    for index in (19, 20, 23, 24, 39, 43, 44, 45):
        require(expanded[index] == reduced[index], f"Original error term {index} changed")
    _assert_unselected(expanded, f"view.{view['Z']}")
    return expanded


def check_packet(report, parent_report):
    require(report["parent_report"] == PARENT and report["parent_report_sha256"] == PARENT_SHA256, "Parent identity changed")
    require(report["source_family"] == parent_report["source_family"], "Source family changed")
    _check_source_guards(report, parent_report)
    require(report["inherited_m15_packet"] == parent_report["inherited_m15_packet"], "Inherited packet changed")
    require(report["inherited_m15_packet_digest"] == parent_report["inherited_m15_packet_digest"], "Inherited packet digest changed")
    require(report["inherited_signed183_views_digest"] == digest(parent_report["signed_m15_query_views"]), "Inherited signed view digest changed")
    require(report["source_primitive_only"] is True and report["signed183_views_unchanged"] is True, "Source-only scope changed")
    coupling = report["endpoint_coupling"]
    require(coupling["applied"] is False and coupling["raw23_raw24_raw44_retained"] is True, "Endpoint coupling scope changed")
    require(coupling["native_lower_equals_macro0_Uy_proved"] is False and coupling["native_lower_macro0_value_and_jump_identity_proved"] is False, "Unproved endpoint identity selected")

    rows = _source_rows(report)
    require(len(rows) == 3, "Source row count changed")
    packet, samples, seam = _ledger_samples(parent_report)
    partial = report["native_macro0_partial_seam_evidence"]
    require(partial["saved_physical_seam_passed"] is True and partial["saved_physical_seam_digest"] == digest(seam), "Physical seam witness changed")
    require(partial["left_scope"] == seam["left_scope"] and partial["right_scope"] == seam["right_scope"], "Native/macro source scope changed")
    require(partial["common_Qy_operator"] == seam["leading_Qy_common"], "Common leading operator changed")
    require(partial["explicit_normalized_Psi_Mpsi_ledger_theorem_bound_here"] is False and partial["raw23_raw24_signed_endpoint_ledger_binding_bound_here"] is False, "Unproved ledger seam selected")
    require(partial["source_values_selected"] is False, "Endpoint value selected")
    views = report["signed_m15_query_views"]
    parent_views = parent_report["signed_m15_query_views"]
    require(views == parent_views and len(views) == 3, "Accepted signed 183 views changed")
    for view, parent_view in zip(views, parent_views):
        expanded = _check_view(view, parent_view, samples[view["Z"]])
        require(len(expanded) == 183, "Expanded signed view size changed")
    require(packet["ledger_hydrated_by_reference"] is True and packet["new_source_interfaces_do_not_change_error_budget"] is True, "Inherited witness flags changed")

    contract = report["contract"]
    require(contract["accepted_ancestor_or_fixture_replayed"] is False, "Accepted ancestor replayed")
    require(contract["mass_core_or_old_integrals_recomputed"] is False, "Old mass/core recomputed")
    require(contract["raw213_reduced183_archives_unchanged"] is True and contract["signed_raw46_raw47_unchanged"] is True, "Original archives changed")
    require(contract["six_inner_lower_reserved_once_target_and_N_unchanged"] is True, "Original witnesses changed")
    require(contract["positive_caps_midpoints_or_interval_overlap_used_as_values"] is False, "Positive cap promoted to value")
    require(contract["hbB_hbS_or_absolute_swirl_amplitude_materialized"] is False, "Width/amplitude materialized")
    require(contract["raw23_raw24_raw43_raw44_raw45_and_all_other_errors_retained"] is True, "Original error budget changed")
    require(all(value is False for value in report["limits"].values()), "Global completion claim introduced")
    return {
        "passed": True,
        "source_rows": 3,
        "ordinary_Z_orders": [0, 1, 2],
        "ordinary_C_orders": [0, 1],
        "positive_kernel_rows_per_source": 6,
        "mode_rates": EXPECTED_RATES,
        "signed_views": 3,
        "view_terms": 183,
        "unchanged_terms_per_view": 181,
        "replaced_parent_slots": [41, 42],
        "raw213_reduced183_archives_unchanged": True,
        "native_endpoint_coupling_applied": False,
        "global_terminal_closure": False,
    }


def pins(report):
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        for name, expected in report[key].items():
            require(sha(base / name) == expected, "Changed pin: " + name)
            count += 1
    return count


def run():
    started = time.perf_counter()
    from check_original_macro_weighted_kernel import run_fixture

    fixture = run_fixture()
    require(fixture["passed"] is True, "Independent macro kernel fixture failed")
    owner = OriginalMacroWeightedM15()
    report = owner.build()
    code = dict(owner.implementation_hashes)
    code.update({name: sha(ROOT / name) for name in CODE})
    report.update(input_hashes=owner.input_hashes, implementation_hashes=code)
    report["checks"] = {
        "passed": True,
        "independent_new_kernel_fixture": fixture,
        "packet_invariants": check_packet(report, owner.parent),
        "no_accepted_ancestor_replay": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path):
    report = load(path)
    require(report["checks"]["passed"] is True and report["checks"]["independent_new_kernel_fixture"]["passed"] is True, "Saved report failed")
    return {
        "passed": True,
        "pins": pins(report),
        "packet_invariants": check_packet(report, parent()),
        "serialized_only": True,
        "owner_or_fixture_or_ancestor_run": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    if args.write_report:
        report = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8", newline="\n")
        temporary.replace(path)
    else:
        print(json.dumps(inspect(path), indent=2))
