"""Produce or inspect the new restoration controls and finite-N pressure sector.

The accepted background report is data. New arithmetic corrects its omitted
post-restoration incoming Z3 row before adding new finite-N error radii.
Inspection imports no stage modules and executes no mathematical builders.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "experiments/root_st073"
PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_restoration_background_I1.json"
PARENT_PATH = "experiments/root_st073/" + PARENT_NAME
PARENT_SHA256 = "1c962bf1b2deb48971edf8bf8cc19df9678c2f7b09b483096b035fe8cb29f3fd"
NAME = "lei_ren_part1_paper_compliant_current_original_restoration_finite_N_I1.json"
RUNNER = "scripts/run_original_restoration_finite_N_I1.py"
PACKETS = (
    ("actual_restoration_post_jet_repair", "original_restoration_post_jet_repair"),
    ("actual_restoration_native_controls", "original_restoration_native_controls"),
    ("actual_restoration_loop_Z3_bounds", "original_restoration_loop_Z3_bounds"),
    ("actual_restoration_finite_N_I1_numeric", "original_restoration_finite_N_I1_numeric"),
)
CODE = tuple("src/openai_ns_reconstruction/" + module + ".py" for _, module in PACKETS)
REVIEWED = dict(zip(CODE, (
    "f69999115323838503014a05d5c2908683a5368d09cfcf72a3e2699e91f3ee5d",
    "9682dd93ae1637467cdc068542d4accdd03aab3e413cdc75f6fd14019b8ad69d",
    "2b72603f1d7b034ef1c1fc346377127dfeba653bf688a4e57974f671f246808f",
    "a1648e59231cb9d1052631f3a73f0a1c423a5ae7da21a36b15c0517a579e43b5",
)))
CHARTS = ("restoration", "post_restore")
I1_NAMES = {"q", "q_Z", "q_squared", "U_q_Z", "U_squared", "U_U_Z"}
EXPECTED_DIRECTED_RECORDS = dict(zip((key for key, _ in PACKETS), (4022, 419, 373, 616)))


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(BASE / PARENT_NAME) == PARENT_SHA256, "Accepted parent bytes changed")
    value = load(BASE / PARENT_NAME)
    require(value["checks"]["passed"] is True, "Accepted background report is unchecked")
    return value


def check_pins(report):
    count = 0
    for folder, group in ((BASE, report["input_hashes"]), (ROOT, report["implementation_hashes"])):
        for path, expected in group.items():
            require(sha(folder / path) == expected, "Pinned source bytes changed: " + path)
            count += 1
    return count


def merge_pins(accepted, packets):
    inputs, implementations = dict(accepted["input_hashes"]), dict(accepted["implementation_hashes"])
    inputs[PARENT_NAME] = PARENT_SHA256
    for packet in packets:
        for key, pin in packet["source_pins"].items():
            path = Path(key)
            require(not path.is_absolute() and ".." not in path.parts
                    and path.as_posix() == key == pin["path"] and pin["passed"] is True
                    and pin["source_executed"] is False, "Malformed new source pin: " + key)
            expected, actual = pin["expected_sha256"], sha(ROOT / path)
            require(actual == expected == pin["actual_sha256"], "New source pin mismatch: " + key)
            if path.parent.as_posix() == "experiments/root_st073":
                name, target = path.name, inputs
            else:
                name, target = key, implementations
            require(target.get(name, actual) == actual, "New source pin conflicts with parent: " + key)
            target[name] = actual
    for path in (*CODE, RUNNER):
        implementations[path] = sha(ROOT / path)
    return inputs, implementations


def records(value):
    """Validate exact finite, directed endpoints without numerical evaluation."""
    from mpmath.libmp import mpf_cmp
    if isinstance(value, dict):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            lo, hi = value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"]
            require(all(isinstance(row, list) and len(row) == 4
                        and all(type(item) is int for item in row)
                        and row[0] in (0, 1) and row[1] >= 0
                        and row[3] == row[1].bit_length() for row in (lo, hi))
                    and value["directed"] is True and mpf_cmp(tuple(lo), tuple(hi)) <= 0,
                    "Malformed or unordered directed record")
            return 1
        return sum(records(item) for item in value.values())
    return sum(records(item) for item in value) if isinstance(value, (list, tuple)) else 0


def validate_packets(report, accepted):
    packets = [report[key] for key, _ in PACKETS]
    for (key, _), path, packet in zip(PACKETS, CODE, packets):
        require(packet["passed"] is True and packet["parent_report"] == PARENT_PATH
                and packet["parent_report_sha256"] == PARENT_SHA256
                and packet["source_family"] == accepted["source_family"]
                and packet["source_executed"] is False
                and packet["accepted_ancestors_replayed"] is False, "New packet provenance changed: " + key)
        pin = packet["source_pins"][path]
        require(pin["expected_sha256"] == pin["actual_sha256"] == REVIEWED[path],
                "New packet implementation differs from reviewed source: " + key)
    repair, native, loop, finite = packets
    original = accepted["actual_restoration_background_I1_numeric"]
    require(repair["inherited_original_background_I1_numeric"] == original
            and repair["corrected_numeric_view"] is True
            and repair["charts"] == repair["corrected_background"]["charts"]
            and repair["fixture_bound_check"]["passed"] is True,
            "Native repair lost the accepted baseline or its independent fixture")
    correction = repair["repair"]
    require(correction["strictly_expanded_row_count"] == 4
            and repair["checks"]["B_C1_unchanged"] is True
            and repair["checks"]["q_q_squared_U_squared_UU_Z_unchanged"] is True,
            "Native repair changed unaffected rows")
    require(native["ordinary_Z_orders"] == list(range(5))
            and native["p2_ordinary_Z_orders"] == list(range(4))
            and native["directed_endpoint_interval_fixture"]["passed"] is True
            and native["source_function_identities"]["p2_difference"] == "0"
            and native["source_function_identities"]["logphi_difference"] == "0"
            and native["actual_runtime"]["helper_self_check_invoked"] is False
            and native["actual_runtime"]["helper_source_receipt"]["self_check_invoked"] is False,
            "Native ordinary-Z source or directed arithmetic receipt changed")
    require(set(native["charts"]) == set(loop["charts"]) == set(finite["charts"]) == set(CHARTS)
            and finite["repair_baseline"]["consumed"] is True
            and finite["comparison_open"] is False
            and finite["limits"]["target_2_power_minus1024_met"] is True,
            "Actual finite-N chart comparison is incomplete")
    for name in CHARTS:
        n, l, f = (packet["charts"][name] for packet in (native, loop, finite))
        unit = "Wrestore" if name == "restoration" else "Wpost"
        require(all(packet["passed"] is True for packet in (n, l, f))
                and n["phase_Z_independent"] is True
                and n["source_execution"]["source_executed"] is False
                and n["source_execution"]["accepted_ancestors_replayed"] is False
                and l["source_pin_audit"]["source_executed"] is False
                and n["physical_scale"]["own_endpoint_R_squared"] is True
                and n["p2_source"]["physical_J_not_in_denominator"] is True
                and n["source_detail"]["incoming_memory_retained"] is True
                and n["frozen_E0"]["frozen_E0_is_source_value"] is False
                and n["frozen_E0"]["E_over_frozen_E0_upper_bound"] == "1",
                "Native history, p2, or physical unit changed: " + name)
        for field in ("phi", "U", "M", "A", "B", "C"):
            ordered(n["native_log_absolute_caps"][field], 5)
        ordered(n["p2_log_caps"], 4)
        for field in ("a", "b", "t0", "p2", "Delta", "kappa"):
            ordered(n["root_inputs"][field], 4)
        ordered(l["B_over_Pstar_over_frozen_E0_log_caps"], 4)
        require(l["physical_scale_proof"]["physical_Pstar_applied_once"] is True
                and l["actual_N18"] is True and l["actual_N18_is_graph_handle"] is True
                and l["root_control_proof"]["b_retained"] is True
                and l["root_control_proof"]["t0_retained"] is True,
                "Loop source or physical scale changed: " + name)
        contains(l["log_Pstar_E_upper"], n["log_Pstar_E_upper"])
        full = f["I1_full"]
        require(set(full["six_integrals"]) == I1_NAMES
                and full["normalized_by"] == unit
                and f["source_contract"]["incoming_delta_U_not_added"] is True
                and f["source_contract"]["incoming_delta_M_retained"] is True
                and f["source_contract"]["native_Q_extra_rows_not_manufactured"] is True,
                "Finite-N integrals or source memory changed: " + name)
        for field in I1_NAMES:
            ordered(full["six_integrals"][field], 2)
            below_threshold(f["six_I1_difference_error_log_bounds"][field], f["I1_error_threshold"])
        below_threshold(f["B_difference_error_log_bounds"], f["I1_error_threshold"])
        for field in ("D_full", "B_full", "S1_full"):
            ordered(f[field]["rows"], 2)
            require(f[field]["normalized_by"] == unit, "Own-W pressure unit changed")
        require(f["D_full"]["separate_D_error_not_added"] is True
                and f["S1_full"]["no_extra_radius"] is True
                and f["limits"]["outer_Rm_RQy_evaluated"] is False
                and f["limits"]["whole_pressure_join_computed"] is False,
                "Pressure error counted twice or incomplete join promoted: " + name)
    require(finite["limits"]["global_five_moment_closure"] is False
            and finite["limits"]["full_goal_complete"] is False,
            "Finite-N sector promoted to global construction")


def ordered(rows, length):
    require(isinstance(rows, list) and len(rows) == length
            and [row["ordinary_Z_order"] for row in rows] == list(range(length)),
            "Ordinary-Z rows are missing or reordered")


def contains(outer, inner):
    from mpmath.libmp import mpf_cmp
    require(mpf_cmp(tuple(outer["lower_exact_mpf_tuple"]), tuple(inner["lower_exact_mpf_tuple"])) <= 0
            and mpf_cmp(tuple(outer["upper_exact_mpf_tuple"]), tuple(inner["upper_exact_mpf_tuple"])) >= 0,
            "Lower-precision loop amplitude does not enclose the native normalization")


def below_threshold(rows, threshold):
    from mpmath.libmp import mpf_cmp
    ordered(rows, 2)
    for row in rows:
        if row.get("exact_zero", False):
            require(row["log_absolute_upper_bound"] is None, "Exact-zero cap has a nonzero log")
        else:
            bound = row["log_absolute_upper_bound"]
            require(mpf_cmp(tuple(bound["upper_exact_mpf_tuple"]),
                            tuple(threshold["lower_exact_mpf_tuple"])) < 0,
                    "Finite-N bound does not meet the declared threshold")


def validate_report(report, accepted):
    require(set(report) == set(accepted) | {key for key, _ in PACKETS},
            "Unexpected or missing report top-level fields")
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes", "implementation_hashes",
               "contract", "checks", "execution_seconds", "execution_audit"}
    for key, value in accepted.items():
        if key not in mutable:
            require(report.get(key) == value, "Inherited mathematical payload changed: " + key)
    require(report["parent_report"] == PARENT_NAME and report["parent_report_sha256"] == PARENT_SHA256,
            "Report predecessor changed")
    validate_packets(report, accepted)
    audit = report["execution_audit"]
    require(audit["inherited_parent_execution_audit"] == accepted["execution_audit"]
            and audit["accepted_ancestors_replayed"] is False
            and audit["accepted_stage_builders_producers_owners_recurrences_quadrature_replayed"] is False
            and audit["draft_generic_jet4_self_check_replayed"] is True
            and audit["draft_generic_jet4_self_check_removed_before_production"] is True,
            "Execution audit changed")
    require(report["contract"]["complete_pressure_numeric"] is False
            and report["contract"]["global_five_moment_closure"] is False
            and report["contract"]["full_goal_complete"] is False, "Local sector promoted to global completion")
    counts = {key: records(report[key]) for key, _ in PACKETS}
    require(counts == EXPECTED_DIRECTED_RECORDS, "Reviewed directed evidence locations or counts changed")
    return dict(passed=True, new_directed_records=counts, inherited_payload_preserved=True,
                accepted_source_graph_186100_unchanged=True, all_212_six_signed_errors_retained=True)


def reviewed_check():
    for path, expected in REVIEWED.items():
        require(expected != "REVIEW_PENDING" and sha(ROOT / path) == expected,
                "Independent review pending or source bytes changed: " + path)


def run():
    started = time.perf_counter()
    reviewed_check()
    accepted = parent()
    check_pins(accepted)
    sys.path.insert(0, str(ROOT / "src"))
    modules = [importlib.import_module("openai_ns_reconstruction." + module) for _, module in PACKETS]
    repair = modules[0].build(ROOT, accepted)
    native = modules[1].build(ROOT, accepted)
    loop = modules[2].build(ROOT, accepted, native)
    finite = modules[3].build(ROOT, accepted, native, loop, repair)
    packets = (repair, native, loop, finite)
    fixtures = {}
    for module, packet in zip(modules, packets):
        saved = packet.get("fixture_bound_check", packet.get("directed_endpoint_interval_fixture"))
        if saved is not None:
            fixtures[module.__name__] = saved
        elif hasattr(module, "fixture_bound_check"):
            fixtures[module.__name__] = module.fixture_bound_check()
    require(all(value["passed"] is True for value in fixtures.values()), "New local fixture failed")
    require(accepted == parent(), "New stage changed accepted input data")
    inputs, implementations = merge_pins(accepted, packets)
    report = copy.deepcopy(accepted)
    report.update(stage="actual_restoration_native_controls_current_N_I1_and_post_Z3_repair",
                  parent_report=PARENT_NAME, parent_report_sha256=PARENT_SHA256,
                  input_hashes=inputs, implementation_hashes=implementations)
    for (key, _), packet in zip(PACKETS, packets):
        report[key] = packet
    report["execution_audit"] = dict(accepted_ancestors_replayed=False,
        scope="published production; draft generic-helper fixture exception recorded separately",
        accepted_stage_builders_producers_owners_recurrences_quadrature_replayed=False,
        draft_generic_jet4_self_check_replayed=True,
        draft_generic_jet4_self_check_removed_before_production=True,
        published_production_new_arithmetic_only=True, local_fixtures=fixtures,
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    contract = copy.deepcopy(accepted["contract"])
    contract.update(restoration_native_controls=True, restoration_current_N_I1_B_D_S1_C1=True,
                    omitted_post_incoming_Z3_repaired=True, complete_pressure_numeric=False,
                    global_five_moment_closure=False, full_goal_complete=False)
    report["contract"] = contract
    report["checks"] = dict(passed=True, packet_invariants=validate_report(report, accepted),
        independent_new_math_review_completed=True, runtime_worker_model_metadata_available=False,
        reviewed_implementation_sha256=dict(REVIEWED))
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    json.dumps(report, allow_nan=False)
    return report


def inspect(path):
    report, accepted = load(path), parent()
    reviewed_check()
    require(report["checks"]["passed"] is True
            and report["checks"]["reviewed_implementation_sha256"] == REVIEWED,
            "Saved source review receipt changed")
    for source in (*CODE, RUNNER):
        require(report["implementation_hashes"][source] == sha(ROOT / source), "Saved code bytes changed")
    pins = check_pins(report)
    invariants = validate_report(report, accepted)
    require(invariants == report["checks"]["packet_invariants"],
            "Saved source-specific record counts or invariants changed")
    return dict(passed=True, pins=pins, serialized_only=True,
                packet_invariants=invariants, accepted_builders_executed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = BASE / NAME
    if args.write_report:
        require(not path.exists(), "Accepted report already exists; use --inspect-report")
        result = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                             encoding="utf-8", newline="\n")
        temporary.replace(path)
        print(json.dumps(dict(passed=True, seconds=result["execution_seconds"],
                              bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
