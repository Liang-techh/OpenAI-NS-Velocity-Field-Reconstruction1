"""Produce or inspect current-N macro six-I1 signed enclosures.

Production executes only this stage's new logarithmic input, loop, and
physical-error adapters. Inspection reads saved packets and pinned bytes.
The inherited endpoint sources and cumulative ledger remain unchanged.
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

from openai_ns_reconstruction import original_macro_loop_Z3_inputs as inputs
from openai_ns_reconstruction import original_macro_loop_Z3_bounds as loop
from openai_ns_reconstruction import original_macro_finite_N_I1_numeric as numeric

PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_macro_background_I1.json"
PARENT_SHA256 = "5a77d3d8b9b294516cc0be05ed83667e2710d64ff8427c7e917b8aed3379aede"
NAME = "lei_ren_part1_paper_compliant_current_original_macro_finite_N_I1.json"
INPUT_KEY = "actual_macro_loop_Z3_inputs"
LOOP_KEY = "actual_macro_loop_Z3_bounds"
NUMERIC_KEY = "actual_macro_finite_N_I1_numeric"
CODE = (
    "src/openai_ns_reconstruction/original_log_Z_phase_bounds.py",
    "src/openai_ns_reconstruction/original_macro_loop_Z3_inputs.py",
    "src/openai_ns_reconstruction/original_macro_loop_Z3_bounds.py",
    "src/openai_ns_reconstruction/original_macro_finite_N_I1_numeric.py",
    "src/openai_ns_reconstruction/original_macro_native_Z4_extension.py",
    "scripts/run_original_macro_finite_N_I1.py",
)
REVIEWED = {
    CODE[0]: "21b442cd8f0044f4652ed7d93d0e54b17a380913925540cb79b17324936fee59",
    CODE[1]: "7117bb8b545e62795d6a99a541dbd7d024ebdc50b734eff902bd48ce8d24b2bd",
    CODE[2]: "e05c37e18ab6874b030780e3e6839d499658d63277b6b6e1db6a1e467293c997",
    CODE[3]: "deedb650c9260822cf556ec1dbf434ee4906803d29b711b39d15ee5e79f4a8ba",
    CODE[4]: "ac0c29264fb791c52559035151a8cffc0653bb7189902d88b7e5de09575485ce",
}


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
    saved, implementations = dict(accepted["input_hashes"]), dict(accepted["implementation_hashes"])
    saved[PARENT_NAME] = PARENT_SHA256
    for packet in packets:
        for key, item in packet["source_pins"].items():
            path = Path(item["path"] if isinstance(item, dict) else key)
            expected = item["expected_sha256"] if isinstance(item, dict) else item
            if len(path.parts) == 1:
                actual, name, target = BASE / path, path.name, saved
            elif path.as_posix().startswith("experiments/root_st073/"):
                actual, name, target = ROOT / path, path.name, saved
            else:
                actual, name, target = ROOT / path, path.as_posix(), implementations
            require(sha(actual) == expected, "New packet source bytes changed: " + str(path))
            require(target.get(name, expected) == expected, "Conflicting source pin: " + name)
            target[name] = expected
    implementations.update({name: sha(ROOT / name) for name in CODE})
    return saved, implementations


def check_records(value):
    """Check directed binary ordering without expanding giant powers."""
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


def validate_inputs(packet, accepted):
    from mpmath.libmp import mpf_cmp
    proof = packet["proofs"]
    require(all(proof[key] is True for key in (
        "kappa_less_than_one_cut_one", "eta_positive_le_one", "a_positive",
        "actual_N18", "actual_N18_is_graph_handle", "actual_phase_Z_independent",
        "real_macro_source_parameters")), "Current macro source premise changed")
    require(mpf_cmp(tuple(packet["kappa_log_upper"]["upper_exact_mpf_tuple"]), (0, 0, 0, 0)) < 0,
            "Saved actual macro kappa<1 proof failed")
    extension = packet["native_Z4_extension"]
    require(extension["passed"] is True and extension["parent_report_sha256"] == PARENT_SHA256,
            "Missing source-related macro native Z4 extension")
    fourth = extension["comparison_fifth_radius_proof"]
    phi_floor = fourth["comparison_phi_floor"]
    require(mpf_cmp(tuple(phi_floor["lower_exact_mpf_tuple"]), (0, 0, 0, 0)) > 0
            and mpf_cmp(tuple(fourth["comparison_phi"][0]["lower_exact_mpf_tuple"]),
                    tuple(phi_floor["lower_exact_mpf_tuple"])) >= 0,
            "New reciprocal-phi cap uses a floor above the accepted comparison floor")
    identity = fourth["L_identity"]
    require(identity["saved_L_encloses_identity"] is True
            and fourth["conditioning"]["reciprocal_rows_le_two"] is True,
            "Missing actual L polynomial or positive reciprocal certificate")
    for saved, computed in zip(identity["saved_rows"], identity["computed_rows"]):
        require(mpf_cmp(tuple(saved["lower_exact_mpf_tuple"]), tuple(computed["lower_exact_mpf_tuple"])) <= 0
                and mpf_cmp(tuple(computed["upper_exact_mpf_tuple"]), tuple(saved["upper_exact_mpf_tuple"])) <= 0,
                "Saved L does not enclose the current polynomial source")
    require(mpf_cmp(tuple(extension["F2_fifth_proof"]["log_upper"]["upper_exact_mpf_tuple"]),
                    (0, 0, 0, 0)) < 0, "Saved F0-squared fifth derivative proof failed")
    for name in ("angular_error_log", "transport_error_log"):
        require(mpf_cmp(tuple(fourth[name]["upper_exact_mpf_tuple"]),
                        tuple(fourth["angular_error_threshold"]["lower_exact_mpf_tuple"])) < 0,
                "Saved macro fourth-row radius failed: " + name)
    fields = extension["native_fields"]
    for name in ("phi", "U", "M", "A", "B", "C"):
        require(len(fields[name]) == 5
                and [row["ordinary_Z_order"] for row in fields[name]] == [0, 1, 2, 3, 4],
                "Missing native macro Z4 field: " + name)
    original = accepted["actual_macro_background_I1_numeric"]["actual_macro_fields"]
    for name, old_name in (("phi", "phi_b"), ("U", "U_b"), ("M", "M_b")):
        for new, old in zip(fields[name][:4], original[old_name]):
            require(all(new["value"][side] == old["value"][side] for side in (
                "lower_exact_mpf_tuple", "upper_exact_mpf_tuple")),
                "Accepted macro native prefix changed: " + name)
    return dict(passed=True, source_branch="kappa<1, cut=1", native_macro_Z4=True,
                accepted_native_prefix_retained=True)


def validate_loop(packet):
    proof = packet["inverse_proof"]
    require(proof["passed"] is True and proof["Z_only"] is True
            and proof["physical_y_derivatives"] is False
            and all(proof[key] is True for key in (
                "ordinary_Z_derivatives", "phase_terms_divided_by_Z_and_angle_factorials",
                "returned_Z_factorials_restored", "delta_zero_is_exact")),
            "Implicit inverse derivative certificate changed")
    require(packet["ordinary_Z_orders"] == [0, 1, 2, 3]
            and packet["physical_y_derivatives_claimed"] is False
            and packet["universal_loop_controls_consumed"] is False
            and packet["universal_B_over_E_value_cap_tightened"] is False,
            "Wrong phase/derivative or generic C0 tightening scope")
    return dict(passed=True, ordinary_Z_orders=[0, 1, 2, 3], fixed_actual_phase=True,
                generic_C0_tightening_used=False)


def validate_numeric(packet, accepted):
    from mpmath.libmp import mpf_cmp
    require(packet["passed"] is True and packet["actual_current_N18"] is True
            and packet["accepted_ancestors_replayed"] is False,
            "Macro finite-N numeric scope changed")
    contract = packet["source_contract"]
    require(contract["passed"] is True and contract["ordinary_Z_orders"] == [0, 1, 2, 3]
            and contract["physical_y_derivatives_claimed"] is False
            and contract["common_N_handle"] == 18
            and all(contract[key] is True for key in (
                "phase_Z_independent", "N_not_materialized", "macro_incoming_nonreset",
                "macro_delta_U_has_no_incoming", "micro_delta_U_not_added_as_macro_inlet",
                "physical_Pstar_applied_once", "q_recovery_linear")),
            "Wrong current-N source, derivative, or physical scaling")
    fields = packet["physical_error_log_bounds"]
    for name, length in (("delta_U", 4), ("delta_M", 4), ("delta_Q", 3)):
        require(len(fields[name]) == length
                and [row["ordinary_Z_order"] for row in fields[name]] == list(range(length)),
                "Missing physical error derivative: " + name)
    I1, D = packet["I1_full"], packet["D_full"]
    require(I1["passed"] is True and D["passed"] is True
            and set(I1["six_integrals"]) == set(numeric.I1_NAMES)
            and I1["signed_background_unchanged"] is True
            and I1["error_radius_is_not_a_selected_value"] is True
            and D["computed_from_full_I1_once"] is True
            and D["separate_D_error_not_added"] is True
            and packet["all_normalized_I1_errors_below_2_power_minus1024"] is True,
            "Six-I1 finite-N radius was not added exactly once")
    for rows in list(I1["six_integrals"].values()) + [D["D_rows"]]:
        require(len(rows) == 2 and [row["ordinary_Z_order"] for row in rows] == [0, 1],
                "Missing ordinary C1 signed result")
    radius = packet["signed_I1_radius"]
    require(radius["lower_exact_mpf_tuple"] == radius["upper_exact_mpf_tuple"] == [0, 1, -1024, 1]
            and I1["symmetric_radius_added_once"] == radius, "Wrong finite-N I1 radius")
    threshold = tuple(packet["I1_error_threshold"]["lower_exact_mpf_tuple"])
    errors = packet["six_I1_difference_error_log_bounds"]
    require(set(errors) == set(numeric.I1_NAMES), "Incomplete saved six-I1 error certificate")
    for name, rows in errors.items():
        require(len(rows) == 2 and [row["ordinary_Z_order"] for row in rows] == [0, 1],
                "Missing saved error row: " + name)
        require(all(mpf_cmp(tuple(row["log_absolute_upper_bound"]["upper_exact_mpf_tuple"]), threshold) < 0
                    for row in rows), "Saved finite-N error exceeds radius: " + name)
    limits = packet["limits"]
    require(limits["finite_N_six_I1_enclosed"] is True
            and limits["outer_RQy_retained"] is True
            and all(limits[key] is False for key in (
                "outer_RQy_evaluated", "complete_pressure_global_closure",
                "global_five_moment_closure", "physical_y_derivatives_claimed")),
            "Finite-N I1 result was promoted to complete pressure or closure")
    target = accepted["actual_micro_macro_correlated_m15"]["actual_macro_six_I1_source"]
    require(contract["source_graph_handles"]["target_from_parent"] == target,
            "Full finite-N target source changed")
    return dict(passed=True, scope="current-N macro six-I1 ordinary C1", directed_records=check_records(packet))


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
    counts = {}
    for key in (INPUT_KEY, LOOP_KEY, NUMERIC_KEY):
        packet = report[key]
        require(packet["passed"] is True and packet["parent_report_sha256"] == PARENT_SHA256,
                "New packet belongs to a different source: " + key)
        counts[key] = check_records(packet)
        require(counts[key] > 0, "Missing directed numeric records: " + key)
    return dict(inputs=validate_inputs(report[INPUT_KEY], accepted), loop=validate_loop(report[LOOP_KEY]),
                numeric=validate_numeric(report[NUMERIC_KEY], accepted), directed_records=counts,
                accepted_math_payload_preserved=True)


def run():
    started = time.perf_counter()
    for path, expected in REVIEWED.items():
        require(expected != "REVIEW_PENDING", "Independent final mathematical review is pending: " + path)
        require(sha(ROOT / path) == expected, "Reviewed source bytes changed: " + path)
    accepted = parent()
    check_pins(accepted)
    source = inputs.build(ROOT, accepted)
    bounds = loop.build(ROOT, accepted, source)
    result = numeric.build(ROOT, accepted, bounds)
    require(accepted == parent(), "New adapters mutated accepted data")
    report = copy.deepcopy(accepted)
    saved, implementations = merge_pins(accepted, source, bounds, result)
    report.update(stage="current_N_macro_six_I1_signed_C1_with_once_added_finite_N_radius",
                  parent_report=PARENT_NAME, parent_report_sha256=PARENT_SHA256,
                  input_hashes=saved, implementation_hashes=implementations)
    report[INPUT_KEY], report[LOOP_KEY], report[NUMERIC_KEY] = source, bounds, result
    report["execution_audit"] = dict(current_turn_accepted_ancestors_replayed=False,
        current_turn_accepted_producer_owner_norm_recurrence_quadrature_replayed=False,
        rejected_unaccepted_input_draft_sha256="6d8c472a55be822504c68be652bef781c7c6e5b0ec534fa1d09320afb0bf06f3",
        pre_acceptance_review_repairs=[
            "Native macro Z4 extension separated from micro inlet; derivative and pressure source omissions repaired",
            "D ordinary C1 products use both Leibniz terms; no repeated D radius",
            "C5 actual AST guards replace a local recipe guard; older generic C0 tightening removed",
        ],
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    report["contract"] = dict(macro_finite_N_six_I1_C1_enclosed=True, macro_complete_pressure_numeric=False,
        macro_switch_RQy_conditional_source_identity=True, terminal_five_moment_closure=False,
        full_goal_complete=False)
    report["checks"] = dict(passed=True, packet_invariants=validate_report(report, accepted),
        independent_new_math_review_completed=True, independent_runtime_worker_model_metadata_available=False,
        reviewed_implementation_sha256=REVIEWED)
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    json.dumps(report, allow_nan=False)
    return report


def inspect(path):
    report, accepted = load(path), parent()
    require(report["checks"]["passed"] is True
            and report["checks"]["independent_new_math_review_completed"] is True
            and report["checks"]["reviewed_implementation_sha256"] == REVIEWED,
            "Saved report review binding changed")
    for path, expected in REVIEWED.items():
        require(report["implementation_hashes"][path] == expected == sha(ROOT / path),
                "Saved reviewed source changed: " + path)
    require(sha(ROOT / CODE[-1]) == report["implementation_hashes"][CODE[-1]], "Runner source pin changed")
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
        temporary.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                             encoding="utf-8", newline="\n")
        temporary.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=result["execution_seconds"],
                              bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
