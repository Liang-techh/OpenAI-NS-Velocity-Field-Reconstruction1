"""Produce or inspect signed native macro background I1 and the next RQy seam.

Only the two new adapters execute during production.  Inspection reads saved
packets, graph rows and pinned bytes.  Finite-N macro error and complete
pressure remain separate, unfinished targets.
"""
from __future__ import annotations

import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction import original_macro_background_I1_numeric as numeric
from openai_ns_reconstruction import original_macro_switch_RQy_seam as seam

PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_micro_macro_correlated_m15.json"
PARENT_SHA256 = "8354af18754496450fb41ff0335e2c82f503c0d64852496da11ee1546fb23577"
NAME = "lei_ren_part1_paper_compliant_current_original_macro_background_I1.json"
NUMERIC_KEY = "actual_macro_background_I1_numeric"
SEAM_KEY = "actual_macro_switch_RQy_seam"
CODE = (
    "src/openai_ns_reconstruction/original_macro_background_I1_numeric.py",
    "src/openai_ns_reconstruction/original_macro_switch_RQy_seam.py",
    "scripts/run_original_macro_background_I1.py",
)
REVIEWED_NUMERIC_SHA256 = "0602ab6b769e30c748517cf588d7675e82a1e1611c4795b950fb37165b5f8560"
REVIEWED_SEAM_SHA256 = "d61e83b29b7af16e5e9e1951f8926f172870ba7a408c618a75fde93b0d9704b4"


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    raw = Path(path).read_bytes()
    return json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)


def parent():
    require(sha(BASE / PARENT_NAME) == PARENT_SHA256, "Accepted parent bytes changed")
    result = load(BASE / PARENT_NAME)
    require(result["checks"]["passed"] is True, "Accepted parent did not pass")
    return result


def check_pins(report):
    count = 0
    for folder, group in ((BASE, report["input_hashes"]), (ROOT, report["implementation_hashes"])):
        for name, expected in group.items():
            require(sha(folder / name) == expected, "Source bytes changed: " + name)
            count += 1
    return count


def merge_pins(accepted, *packets):
    inputs, implementations = dict(accepted["input_hashes"]), dict(accepted["implementation_hashes"])
    inputs[PARENT_NAME] = PARENT_SHA256
    for packet in packets:
        for key, item in packet["source_pins"].items():
            path = item["path"] if isinstance(item, dict) else key
            expected = item["expected_sha256"] if isinstance(item, dict) else item
            path = Path(path)
            if len(path.parts) == 1:
                actual, name, target = BASE / path, path.name, inputs
            elif path.as_posix().startswith("experiments/root_st073/"):
                actual, name, target = ROOT / path, path.name, inputs
            else:
                actual, name, target = ROOT / path, path.as_posix(), implementations
            require(sha(actual) == expected, "New packet source bytes changed: " + str(path))
            require(target.get(name, expected) == expected, "Conflicting source pin: " + name)
            target[name] = expected
    implementations.update({name: sha(ROOT / name) for name in CODE})
    return inputs, implementations


def check_records(value):
    """Check exact binary record ordering without expanding giant powers."""
    from mpmath.libmp import mpf_cmp
    count = 0
    if isinstance(value, dict):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            lo, hi = value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"]
            require(all(isinstance(t, list) and len(t) == 4 and all(type(x) is int for x in t)
                        and t[0] in (0, 1) and t[1] >= 0 and t[3] >= 0 for t in (lo, hi)),
                    "Malformed or nonfinite directed interval")
            require(mpf_cmp(tuple(lo), tuple(hi)) <= 0 and value.get("directed") is True,
                    "Reversed or nondirected interval")
            count += 1
        count += sum(check_records(item) for item in value.values())
    elif isinstance(value, list):
        count += sum(check_records(item) for item in value)
    return count


def references(row):
    fields = {
        "sum": ("arguments",), "product": ("arguments",),
        "negative": ("argument",), "analytic_unary": ("argument",),
        "positive_quotient": ("numerator", "denominator"),
        "function_substitution": ("expression", "variable", "value"),
        "definite_integral": ("integrand", "variable", "lower", "upper"),
    }
    result = []
    for key in fields.get(row.get("operation"), ()):
        item = row[key]
        result.extend(item if isinstance(item, list) else [item])
    return [item for item in result if type(item) is int]


def validate_seam(packet, accepted):
    require(packet["passed"] is True and packet["parent_report_sha256"] == PARENT_SHA256,
            "Wrong macro/switch seam parent")
    patch = packet["graph_patch"]
    old = accepted["actual_micro_macro_correlated_m15"]["merged_namespace"]
    archive = load(ROOT / old["base_archive"])
    nodes = archive["exact_graph_nodes"] + old["appended_nodes"]
    require(len(nodes) == old["merged_node_count"] == patch["base_node_count"] == 184687,
            "Seam does not extend the complete accepted namespace")
    require(patch["base_parent_report"] == PARENT_NAME
            and patch["base_parent_report_sha256"] == PARENT_SHA256
            and patch["base_parent_merged_node_count"] == 184687
            and patch["base_parent_merged_namespace_digest"]
            == accepted["actual_micro_macro_correlated_m15"]["graph_namespace"]["namespace_digest"],
            "Seam base parent namespace binding changed")
    require(patch["base_prefix_unchanged"] is True
            and len(patch["appended_nodes"]) == patch["appended_node_count"],
            "Seam suffix serialization changed")
    for row in patch["appended_nodes"]:
        require(all(0 <= ref < len(nodes) for ref in references(row)), "Invalid new graph reference")
        nodes.append(row)
    ports = patch["endpoint_handles"]
    require(ports["frozen_macro_upper_used_Z0_Z1"] == [141114, 141115],
            "Macro outer source was replaced")
    require(len(packet["proof_rows"]) == 2 and packet["symbolic_controls"]["passed"] is True,
            "Missing ordinary C1 source identities")
    for order, proof in enumerate(packet["proof_rows"]):
        require(proof["ordinary_Z_order"] == order and proof["source_identity"]["passed"] is True,
                "Seam conditional source row failed")
        left, right = proof["frozen_macro_upper_m1"], proof["first_switch_lower_s0"]
        require(left == 141114 + order and right == ports["first_switch_lower_new_s0_Z0_Z1"][order],
                "Seam proof endpoint binding changed")
        row = nodes[proof["signed_boundary_term"]]
        require(row["operation"] == "sum" and len(row["arguments"]) == 2
                and nodes[row["arguments"][0]] == {"operation": "negative", "argument": left}
                and row["arguments"][1] == right, "Wrong opposite FTC trace signs")
    contract = packet["contract"]
    require(contract["ordinary_Z_orders"] == [0, 1]
            and contract["frozen_macro_native_J"] == 66987 and contract["first_switch_native_J"] == 9
            and all(contract[key] is True for key in (
                "native_J_remain_distinct", "same_phase_source_N", "nonreset_M_retained",
                "first_switch_upper_standalone_retained", "macro_outer_upper_standalone_retained",
                "identity_only_in_adjacent_joined_expression", "old_caps_and_parent_graph_retained",
                "no_old_owner_checker_producer_executed"))
            and contract["source_values_selected"] is False, "Seam scope changed")
    return dict(passed=True, serialized_nodes=len(nodes), ordinary_Z_orders=[0, 1],
                physical_FTC_signs=True, native_J_distinct=True, source_only=True)


def validate_numeric(packet, accepted):
    require(packet["passed"] is True, "Macro background numeric adapter failed")
    json.dumps(packet, allow_nan=False)
    require(packet["parent_report_sha256"] == PARENT_SHA256
            and packet["source_family"] == accepted["source_family"]
            and packet["source_contract"]["passed"] is True
            and packet["actual_native_pre_finite_N_macro_I1"] is True
            and packet["signed_macro_I1_enclosures"] is True
            and packet["normalized_by_Wmacro"] is True
            and packet["accepted_ancestors_replayed"] is False, "Macro source or unit scope changed")
    require(packet["domain"]["m"] == ["0", "1"]
            and packet["domain"]["Z"] == ["37/100", "93/250"]
            and packet["domain"]["I1_ordinary_Z_orders"] == [0, 1]
            and packet["domain"]["pre_finite_N_native_background"] is True,
            "Macro ordinary derivative domain changed")
    integrals = packet["macro_I1"]
    require(integrals == packet["I1"] and integrals["D_rows"] == packet["D"]
            and integrals["normalized_by_Wmacro"] is True
            and integrals["positive_geometry_multiplied_once"] is True
            and integrals["ordinary_Z_orders"] == [0, 1]
            and integrals["finite_N_correction_added"] is False, "Macro I1 measure or scope changed")
    names = {"q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z"}
    require(set(integrals["six_integrals"]) == names, "Incomplete six-I1 source")
    for rows in list(integrals["six_integrals"].values()) + [integrals["D_rows"]]:
        require(isinstance(rows, list) and len(rows) == 2
                and [row["ordinary_Z_order"] for row in rows] == [0, 1], "Incomplete ordinary C1 result")
    fields = packet["actual_macro_fields"]
    require(len(fields["U_b"]) == len(fields["M_b"]) == 4 and len(fields["Q_b"]) == 3
            and fields["M_is_nonreset_hull"] is True
            and fields["native_M_not_comparison_M"] is True,
            "Native M3 or Q2 source is missing")
    comparison = packet["comparison"]
    require(comparison["decay_is_Z_independent"] is True and comparison["rows_are_convex_hulls"] is True
            and set(comparison["history_rows"]) == {"H", "K", "M", "A", "B", "C"}
            and all(len(rows) >= 5 for rows in comparison["history_rows"].values()),
            "Macro differentiated drive requires source histories through Z4")
    limits = packet["limits"]
    require(all(limits[key] is False for key in (
        "finite_N_correction_added", "finite_N_error_radius_available", "full_finite_N_enclosed",
        "outer_RQy_deleted", "complete_pressure_global_closure", "global_five_moment_closure")),
        "Background result was promoted to full finite-N or global closure")
    target = packet["target_full_finite_N_graph"]
    original = accepted["actual_micro_macro_correlated_m15"]["actual_macro_six_I1_source"]
    require(target["numerically_evaluated"] is False and target["retained_for_future_comparison"] is True
            and target["target_source_handles"] == {key: original[key] for key in (
                "physical_fields", "source", "integrals", "physical_measure")},
            "Full finite-N graph was substituted or mislabelled")
    count = check_records(packet)
    require(count > 40, "Missing macro signed source and error records")
    return dict(passed=True, directed_binary_records=count, scope="pre-finite-N native macro background")


def validate_report(report, accepted):
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes", "implementation_hashes",
               "contract", "checks", "execution_seconds", "execution_audit"}
    for key, value in accepted.items():
        if key not in mutable:
            require(report.get(key) == value, "Accepted payload changed: " + key)
    require(report["parent_report"] == PARENT_NAME and report["parent_report_sha256"] == PARENT_SHA256,
            "Wrong accepted predecessor")
    require(report["execution_audit"]["current_turn_accepted_ancestors_replayed"] is False
            and report["execution_audit"]["inherited_parent_execution_audit"] == accepted["execution_audit"],
            "Current or inherited execution audit changed")
    return dict(numeric=validate_numeric(report[NUMERIC_KEY], accepted), seam=validate_seam(report[SEAM_KEY], accepted),
                accepted_math_payload_preserved=True)


def run():
    started = time.perf_counter()
    require(REVIEWED_NUMERIC_SHA256 != "REVIEW_PENDING" and REVIEWED_SEAM_SHA256 != "REVIEW_PENDING",
            "Independent final mathematical review is pending")
    require(sha(ROOT / CODE[0]) == REVIEWED_NUMERIC_SHA256
            and sha(ROOT / CODE[1]) == REVIEWED_SEAM_SHA256, "Reviewed source bytes changed")
    accepted = parent()
    check_pins(accepted)
    numbers, endpoint = numeric.build(ROOT, accepted), seam.build(ROOT, accepted)
    require(accepted == parent(), "New adapters mutated accepted data")
    report = copy.deepcopy(accepted)
    inputs, implementations = merge_pins(accepted, numbers, endpoint)
    report.update(stage="native_macro_background_signed_I1_C1_and_macro_switch_RQy_source",
                  parent_report=PARENT_NAME, parent_report_sha256=PARENT_SHA256,
                  input_hashes=inputs, implementation_hashes=implementations)
    report[NUMERIC_KEY], report[SEAM_KEY] = numbers, endpoint
    report["execution_audit"] = dict(current_turn_accepted_ancestors_replayed=False,
        current_turn_accepted_producer_owner_norm_recurrence_quadrature_replayed=False,
        unaccepted_report_interface_repairs=[
            "Q4 archive lookup corrected from graph.nodes to exact_graph_nodes before any report write",
            "Raw interval Wmacro metadata changed to JSON integer 10000; exact numeric record and arithmetic unchanged",
        ],
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    report["contract"] = dict(native_macro_background_I1_C1=True, macro_finite_N_I1_enclosed=False,
        macro_switch_RQy_conditional_source_identity=True, full_pressure_numeric=False,
        terminal_five_moment_closure=False, full_goal_complete=False)
    report["checks"] = dict(passed=True, packet_invariants=validate_report(report, accepted),
        independent_new_math_review_completed=True, independent_runtime_worker_model_metadata_available=False,
        reviewed_numeric_sha256=REVIEWED_NUMERIC_SHA256, reviewed_seam_sha256=REVIEWED_SEAM_SHA256)
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path):
    report, accepted = load(path), parent()
    require(report["checks"]["passed"] is True
            and report["checks"]["independent_new_math_review_completed"] is True
            and report["checks"]["reviewed_numeric_sha256"] == REVIEWED_NUMERIC_SHA256
            == report["implementation_hashes"][CODE[0]]
            and report["checks"]["reviewed_seam_sha256"] == REVIEWED_SEAM_SHA256
            == report["implementation_hashes"][CODE[1]]
            and sha(ROOT / CODE[2]) == report["implementation_hashes"][CODE[2]],
            "Saved report review or runner source pin changed")
    return dict(passed=True, pins=check_pins(report), packet_invariants=validate_report(report, accepted),
                serialized_only=True, accepted_builders_producers_or_checkers_replayed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = BASE / NAME
    if args.write_report:
        result = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        temporary.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=result["execution_seconds"],
                              bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
