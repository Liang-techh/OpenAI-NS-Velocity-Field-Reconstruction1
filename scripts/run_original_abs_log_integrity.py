"""Publish an arithmetic integrity correction without replaying ancestor stages.

Historical parent results are retained as data, not accepted as current
certification. Inspection is serialized-only and imports no stage module.
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
PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_restoration_finite_N_I1.json"
PARENT_PATH = "experiments/root_st073/" + PARENT_NAME
PARENT_SHA256 = "fbede8e8c948a2025c65671bce6265af84e57fcb25d4a7e497bf456769d647e4"
NAME = "lei_ren_part1_paper_compliant_current_original_abs_log_integrity.json"
RUNNER = "scripts/run_original_abs_log_integrity.py"
KEY = "actual_abs_log_integrity"
PRIMITIVE = "src/openai_ns_reconstruction/original_directed_abs_log.py"
AUDIT = "src/openai_ns_reconstruction/original_abs_log_integrity.py"
REVIEWED = {PRIMITIVE: "96126cf50406cfc87c8d9af854a42ac9b4ed9a0222b55587d5c913bfcd9ec767",
            AUDIT: "485551576307cbe7be545fda443521d0d1df532576dcf780b4a38ad991387fd5"}
EXPECTED_RECORDS = 8


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            result.update(block)
    return result.hexdigest()


def load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(BASE/PARENT_NAME) == PARENT_SHA256, "Historical parent bytes changed")
    return load(BASE/PARENT_NAME)


def check_pins(report):
    count = 0
    for folder, pins in ((BASE, report["input_hashes"]), (ROOT, report["implementation_hashes"])):
        for path, expected in pins.items():
            require(sha(folder/path) == expected, "Pinned bytes changed: "+path)
            count += 1
    return count


def reviewed_check():
    for path, expected in REVIEWED.items():
        require(expected != "REVIEW_PENDING" and sha(ROOT/path) == expected,
                "Independent source review pending or bytes changed: "+path)


def merge_pins(accepted, packet):
    inputs, implementations = dict(accepted["input_hashes"]), dict(accepted["implementation_hashes"])
    inputs[PARENT_NAME] = PARENT_SHA256
    for key, pin in packet["source_pins"].items():
        path = Path(key)
        require(not path.is_absolute() and ".." not in path.parts
                and path.as_posix() == key == pin["path"] and pin["passed"] is True
                and pin["source_executed"] is False, "Noncanonical source pin: "+key)
        actual = sha(ROOT/path)
        require(actual == pin["expected_sha256"] == pin["actual_sha256"], "Source pin changed: "+key)
        name, target = (path.name, inputs) if path.parent.as_posix() == "experiments/root_st073" \
            else (key, implementations)
        require(target.get(name, actual) == actual, "Source pin conflicts with historical bytes")
        target[name] = actual
    for path in (*REVIEWED, RUNNER):
        implementations[path] = sha(ROOT/path)
    return inputs, implementations


def records(value):
    from mpmath.libmp import mpf_cmp
    if isinstance(value, dict):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            lo, hi = value["lower_exact_mpf_tuple"], value["upper_exact_mpf_tuple"]
            require(all(isinstance(row, list) and len(row) == 4
                        and all(type(item) is int for item in row)
                        and row[0] in (0,1) and row[1] >= 0
                        and row[3] == row[1].bit_length() for row in (lo,hi))
                    and value["directed"] is True and mpf_cmp(tuple(lo),tuple(hi)) <= 0,
                    "Malformed exact directed record")
            return 1
        return sum(records(item) for item in value.values())
    return sum(records(item) for item in value) if isinstance(value,(list,tuple)) else 0


def validate(report, accepted):
    require(set(report) == set(accepted)|{KEY}, "Unexpected report fields")
    mutable = {"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes",
               "contract","checks","execution_seconds","execution_audit"}
    for key,value in accepted.items():
        if key not in mutable:
            require(report.get(key) == value, "Historical mathematical data changed: "+key)
    require(report["parent_report"] == PARENT_NAME and report["parent_report_sha256"] == PARENT_SHA256,
            "Historical predecessor changed")
    packet = report[KEY]
    require(packet["passed"] is True and packet["parent_report"] == PARENT_PATH
            and packet["parent_report_sha256"] == PARENT_SHA256
            and packet["source_family"] == accepted["source_family"], "Audit source binding changed")
    for path,expected in REVIEWED.items():
        pin = packet["source_pins"][path]
        require(pin["expected_sha256"] == pin["actual_sha256"] == expected, "Reviewed audit bytes changed")
    from mpmath.libmp import mpf_cmp
    require(packet["source_executed"] is False and packet["accepted_ancestors_replayed"] is False
            and packet["passed_semantics"] == "integrity_audit_valid_not_legacy_bounds_valid"
            and packet["inherited_log_cap_certification_verified"] is False
            and packet["revalidation_required"] is True
            and packet["old_caps_and_data_preserved"] is True
            and packet["no_new_Rm_certification"] is True
            and packet["parent_scope"]["historical_pass_flag_not_consumed_as_certification"] is True,
            "Integrity audit scope or inherited certification status changed")
    witness = packet["legacy_abs_log_witness"]
    require(witness["passed"] is True and witness["witness_count"] == 2
            and witness["global_mp_precision_restored"] is True
            and witness["old_helper_fixture_not_called"] is True
            and witness["old_helper_called_only_on_new_witnesses"] is True
            and [row["rational"] for row in witness["witnesses"]] == ["7/10","1/3"],
            "Fresh underbound witness scope changed")
    for row in witness["witnesses"]:
        old = row["legacy_log_cap"]["upper_exact_mpf_tuple"]
        reference = row["rational_reference_log_cap"]["lower_exact_mpf_tuple"]
        directed = row["directed_reference_log_cap"]["lower_exact_mpf_tuple"]
        require(old == row["legacy_upper_exact_mpf_tuple"]
                and reference == row["rational_reference_lower_exact_mpf_tuple"]
                and directed == row["directed_lower_exact_mpf_tuple"]
                and row["scalar_precision_bits"] == 53 and row["interval_precision_bits"] == 1600
                and mpf_cmp(tuple(old),tuple(reference)) < 0
                and mpf_cmp(tuple(old),tuple(directed)) < 0,
                "Legacy cap no longer strictly underbounds the rational reference")
    receipt = packet["self_hash_fixture_receipt"]
    fixture,gate = receipt["fixture"],receipt["gate_fixture"]
    require(receipt["self_sha256"] == REVIEWED[AUDIT]
            and fixture["passed"] is True and fixture["rational_cases"] == 18
            and fixture["scalar_precisions"] == [37,53,113]
            and fixture["scalar_abs_called"] is False
            and fixture["accepted_self_checks_replayed"] is False
            and gate["passed"] is True and gate["historical_passed_flag_rejected"] is True
            and gate["schema_only_forged_receipt_rejected"] is True
            and gate["explicit_new_receipt_accepted"] is False
            and gate["supported_verified_revalidation_methods"] == []
            and gate["no_revalidation_receipt_admitted_yet"] is True,
            "Primitive fixture or no-bypass certification gate changed")
    scan = packet["direct_old_abs_log_static_scan"]
    require(scan["passed"] is True and scan["pinned_python_path_count"] == 431
            and scan["pinned_python_paths_sha_checked"] == 431
            and scan["affected_direct_module_count"] == 2
            and scan["direct_old_abs_log_call_count"] == 23
            and scan["transitive_calls_not_scanned"] is True
            and scan["dynamic_aliases_not_scanned"] is True,
            "Bounded source impact scan changed")
    expected_calls = {"src/openai_ns_reconstruction/original_long_reference_native_controls.py":17,
                      "src/openai_ns_reconstruction/original_restoration_native_controls.py":6}
    require({row["path"]:row["direct_call_count"] for row in scan["affected_direct_modules"]} == expected_calls,
            "Direct source impact changed")
    for row in scan["affected_direct_modules"]:
        require(row["expected_sha256"] == row["actual_sha256"]
                == accepted["implementation_hashes"][row["path"]]
                and len(row["call_lines"]) == row["direct_call_count"], "Direct consumer evidence changed")
    contract,audit = report["contract"],report["execution_audit"]
    require(contract["inherited_log_caps_require_revalidation"] is True
            and contract["inherited_log_caps_currently_certified"] is False
            and contract["finite_N_error_thresholds_currently_certified"] is False
            and contract["new_Rm_pipeline_accepted"] is False
            and contract["complete_pressure_numeric"] is False
            and contract["global_five_moment_closure"] is False
            and contract["full_goal_complete"] is False
            and audit["accepted_ancestors_replayed"] is False
            and audit["historical_parent_passed_flag_consumed_as_certification"] is False
            and audit["inherited_parent_execution_audit"] == accepted["execution_audit"],
            "Historical results promoted to currently verified construction")
    count = records(packet)
    require(count == EXPECTED_RECORDS, "Directed witness record locations/count changed")
    return dict(passed=True,scope="arithmetic integrity only",directed_witness_records=count,
                inherited_payload_preserved=True,old_212_six_signed_errors_retained=True,
                inherited_log_caps_require_revalidation=True,new_Rm_pipeline_accepted=False)


def run():
    started = time.perf_counter()
    reviewed_check()
    accepted = parent()
    check_pins(accepted)
    sys.path.insert(0,str(ROOT/"src"))
    module = importlib.import_module("openai_ns_reconstruction.original_abs_log_integrity")
    packet = module.build(ROOT,accepted)
    require(accepted == parent(), "Integrity audit changed historical input")
    inputs,implementations = merge_pins(accepted,packet)
    report = copy.deepcopy(accepted)
    report.update(stage="directed_absolute_log_primitive_and_inherited_certification_integrity",
                  parent_report=PARENT_NAME,parent_report_sha256=PARENT_SHA256,
                  input_hashes=inputs,implementation_hashes=implementations)
    report[KEY] = packet
    report["execution_audit"] = dict(accepted_ancestors_replayed=False,
        historical_parent_passed_flag_consumed_as_certification=False,
        new_rational_witnesses_and_source_scan_only=True,
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    contract = copy.deepcopy(accepted["contract"])
    contract.update(inherited_log_caps_require_revalidation=True,inherited_log_caps_currently_certified=False,
                    finite_N_error_thresholds_currently_certified=False,
                    new_Rm_pipeline_accepted=False,complete_pressure_numeric=False,
                    global_five_moment_closure=False,full_goal_complete=False)
    report["contract"] = contract
    report["checks"] = dict(passed=True,scope="new arithmetic primitive and integrity audit only",
        historical_numeric_claims_revalidated=False,reviewed_implementation_sha256=dict(REVIEWED),
        packet_invariants=validate(report,accepted),runtime_worker_model_metadata_available=False)
    report["execution_seconds"] = round(time.perf_counter()-started,3)
    json.dumps(report,allow_nan=False)
    return report


def inspect(path):
    report,accepted = load(path),parent()
    reviewed_check()
    require(report["checks"]["passed"] is True
            and report["checks"]["historical_numeric_claims_revalidated"] is False
            and report["checks"]["reviewed_implementation_sha256"] == REVIEWED,
            "Saved audit scope changed")
    pins = check_pins(report)
    invariants = validate(report,accepted)
    require(invariants == report["checks"]["packet_invariants"], "Saved audit evidence changed")
    return dict(passed=True,scope="integrity audit only; historical log caps require revalidation",
                pins=pins,serialized_only=True,packet_invariants=invariants,accepted_builders_executed=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report",action="store_true")
    group.add_argument("--inspect-report",action="store_true")
    args = parser.parse_args()
    path = BASE/NAME
    if args.write_report:
        require(not path.exists(), "Report exists; use --inspect-report")
        result = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                             encoding="utf-8",newline="\n")
        temporary.replace(path)
        print(json.dumps(dict(passed=True,seconds=result["execution_seconds"],bytes=path.stat().st_size,sha256=sha(path))))
    else:
        print(json.dumps(inspect(path),indent=2))
