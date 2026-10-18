"""Publish scoped directed repairs while retaining historical construction data.

This stage repairs named numerical blocks, not entire inherited certificates.
Inspection reads serialized evidence only; it never imports a stage builder.
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
PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_abs_log_integrity.json"
PARENT_PATH = "experiments/root_st073/" + PARENT_NAME
PARENT_SHA256 = "dd4cb5741cd7e04c112d2968bf9a81aaceecfac3351a7ca93475c5c9aae646dd"
NAME = "lei_ren_part1_paper_compliant_current_original_native_control_directed_repair.json"
RUNNER = "scripts/run_original_native_control_directed_repair.py"
KEY = "actual_native_control_directed_repair"
REVIEWED = {
    "src/openai_ns_reconstruction/original_directed_log_reciprocal.py": "5f948b11559c0c7533f3b6f98a210a18b06f3a959ae96be31d3ef139fcad72df",
    "src/openai_ns_reconstruction/original_native_control_directed_repair.py": "3f13a01fd6c4d85ca526970fedd8baffa0245176336e5dfe3696ad92dbe4dfb4",
}


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
    require(sha(BASE / PARENT_NAME) == PARENT_SHA256, "Corrective parent bytes changed")
    return load(BASE / PARENT_NAME)


def check_pins(report):
    count = 0
    for folder, pins in ((BASE, report["input_hashes"]),
                         (ROOT, report["implementation_hashes"])):
        for path, expected in pins.items():
            require(sha(folder / path) == expected, "Pinned bytes changed: " + path)
            count += 1
    return count


def reviewed_check():
    for path, expected in REVIEWED.items():
        require(expected != "REVIEW_PENDING" and sha(ROOT / path) == expected,
                "Independent source review pending or bytes changed: " + path)


def merge_pins(accepted, packet):
    inputs = dict(accepted["input_hashes"])
    implementations = dict(accepted["implementation_hashes"])
    inputs[PARENT_NAME] = PARENT_SHA256
    for key, pin in packet["source_pins"].items():
        path = Path(key)
        require(not path.is_absolute() and ".." not in path.parts
                and path.as_posix() == key == pin["path"]
                and pin["passed"] is True and pin["source_executed"] is False,
                "Noncanonical source pin: " + key)
        actual = sha(ROOT / path)
        require(actual == pin["expected_sha256"] == pin["actual_sha256"],
                "Source pin changed: " + key)
        name, target = (path.name, inputs) if path.parent.as_posix() == "experiments/root_st073" \
            else (key, implementations)
        require(target.get(name, actual) == actual, "Historical source pin conflicts")
        target[name] = actual
    for path in (*REVIEWED, RUNNER):
        implementations[path] = sha(ROOT / path)
    return inputs, implementations


def finite_tuple(value):
    require(isinstance(value, (list, tuple)) and len(value) == 4
            and all(type(part) is int for part in value)
            and value[0] in (0, 1) and value[1] >= 0
            and value[3] == value[1].bit_length()
            and (tuple(value) == (0, 0, 0, 0) if value[1] == 0 else value[1] % 2 == 1),
            "Malformed finite canonical MPF tuple")
    return tuple(value)


def directed_record(value):
    from mpmath.libmp import mpf_cmp
    require(isinstance(value, dict) and value.get("directed") is True,
            "Canonical directed interval record required")
    lo = finite_tuple(value["lower_exact_mpf_tuple"])
    hi = finite_tuple(value["upper_exact_mpf_tuple"])
    require(mpf_cmp(lo, hi) <= 0, "Reversed directed interval record")


def record_count(value):
    if isinstance(value, dict):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            directed_record(value)
            return 1
        return sum(record_count(item) for item in value.values())
    return sum(record_count(item) for item in value) if isinstance(value, (list, tuple)) else 0


def validate_packet(packet, accepted):
    from mpmath.libmp import mpf_cmp
    require(packet["passed"] is True
            and packet["schema"] == "original_native_control_directed_repair_v1"
            and packet["passed_semantics"] ==
            "local_directed_abs_input_inverse_L_and_restoration_inverse_Utheta_valid; "
            "inherited_and_downstream_caps_not_recertified"
            and packet["source_family"] == accepted["source_family"]
            and packet["parent_scope"]["report"] == PARENT_PATH
            and packet["parent_scope"]["report_sha256"] == PARENT_SHA256
            and packet["parent_scope"]["exact_parent_mapping_checked"] is True
            and packet["accepted_ancestors_replayed"] is False
            and packet["source_executed"] is False
            and packet["full_goal_complete"] is False,
            "Corrective source binding/scope changed")
    charts = packet["charts"]
    require(set(charts) == {"long_reference:long_reshape", "long_reference:reference",
                            "restoration:restoration", "restoration:post_restore"},
            "Corrective chart set changed")
    shared = charts["long_reference:long_reshape"]["directed_L_block"]
    for key, chart in charts.items():
        family, name = key.split(":")
        require(chart["family"] == family and chart["chart"] == name
                and chart["ordinary_Z_orders"] == list(range(5))
                and chart["downstream_recertification"] is False
                and chart["directed_L_block"] == shared,
                "Chart dependency or local scope changed")
        old_section = "actual_" + ("long_reference" if family == "long_reference"
                                    else "restoration") + "_native_controls"
        old = accepted[old_section]["charts"][name]
        require(chart["historical_inverse_L_cap"]["rows"] == old["p2_source"]["inverse_L"],
                "Historical inverse-L evidence was replaced")
        for leaf in chart["signed_leaf_checkpoints"].values():
            require(leaf["ordinary_Z_orders"] == list(range(5))
                    and len(leaf["rows"]) == 5 and leaf["signed_value_selected"] is False,
                    "Signed leaf rows changed")
            for k, row in enumerate(leaf["rows"]):
                directed_record(row["signed_interval"])
                if row["directed_abs_log_upper"] is not None:
                    directed_record(row["directed_abs_log_upper"])
                require(row["ordinary_Z_order"] == k
                        and row["exact_zero"] == (row["directed_abs_log_upper"] is None),
                        "Signed leaf zero/log status changed")
        if family == "restoration":
            block = chart["directed_inverse_Utheta_block"]
            require(block["formula"] == "1/Utheta=(1+Z^2)*exp(-logP-offset/10)"
                    and block["offset_domain"] == (["-8", "-7"] if name == "restoration"
                                                   else ["-7", "-6"])
                    and block["ordinary_Z_orders"] == list(range(5))
                    and block["exact_zero_orders"] == [3, 4]
                    and block["physical_y_derivatives_claimed"] is False
                    and block["inherited_relative_caps_used"] is False
                    and block["stored_logarithms_exponentiated"] is False
                    and block["complete_p2_recertified"] is False
                    and block["historical_inverse_Utheta_caps"] == old["p2_source"]["inverse_Utheta"],
                    "Exact polynomial inverse-Utheta scope changed")
            for field in ("logP_signed_interval", "offset_signed_interval",
                          "log_amplitude_interval", "inverse_log_amplitude_upper"):
                directed_record(block[field])
            require(len(block["polynomial_signed_rows"]) == 5
                    and len(block["polynomial_abs_log_rows"]) == 5,
                    "Polynomial source row count changed")
            for row in block["polynomial_signed_rows"]:
                directed_record(row)
            for row in block["polynomial_abs_log_rows"]:
                if row is not None:
                    directed_record(row)
            lower = tuple(block["log_amplitude_interval"]["lower_exact_mpf_tuple"])
            flipped = (1-lower[0], *lower[1:]) if lower[1] else (0, *lower[1:])
            inv = block["inverse_log_amplitude_upper"]
            require(tuple(inv["lower_exact_mpf_tuple"]) == flipped
                    and tuple(inv["upper_exact_mpf_tuple"]) == flipped,
                    "Inverse-Utheta amplitude lost exact directed endpoint")
            rows = block["directed_inverse_Utheta_rows"]
            require(len(rows) == 5 and [r["ordinary_Z_order"] for r in rows] == list(range(5))
                    and [r["exact_zero"] for r in rows] == [False, False, False, True, True]
                    and all((r["log_absolute_upper_bound"] is None) == r["exact_zero"] for r in rows),
                    "Polynomial inverse-Utheta derivative rows changed")
            for row in rows:
                if row["log_absolute_upper_bound"] is not None:
                    directed_record(row["log_absolute_upper_bound"])
    require(shared["ordinary_Z_orders"] == list(range(5))
            and shared["signed_value_selected"] is False
            and shared["inherited_log_caps_used_as_signed_inputs"] is False,
            "L block scope changed")
    directed_record(shared["z_signed_interval"])
    directed_record(shared["delta_signed_interval"])
    require(len(shared["signed_L_rows"]) == 5, "Signed L row count changed")
    for k, row in enumerate(shared["signed_L_rows"]):
        require(row["ordinary_Z_order"] == k and row["exact_zero"] == (k >= 3),
                "Signed L derivative order changed")
        directed_record(row["signed_interval"])
    for field in ("directed_abs_log_upper_rows", "directed_inverse_L_rows"):
        rows = shared[field]
        zeros = [False, False, False, True, True] if field == "directed_abs_log_upper_rows" else [False]*5
        require(len(rows) == 5 and [r["ordinary_Z_order"] for r in rows] == list(range(5))
                and [r["exact_zero"] for r in rows] == zeros
                and all((r["log_absolute_upper_bound"] is None) == r["exact_zero"] for r in rows),
                "L reciprocal derivative rows changed")
        for row in rows:
            if row["log_absolute_upper_bound"] is not None:
                directed_record(row["log_absolute_upper_bound"])
    fixture = packet["fixture_bound_check"]
    require(fixture["passed"] is True and fixture["rational_cases"] == 18
            and fixture["scalar_precisions"] == [37, 53, 113]
            and fixture["interval_precision_bits"] == 1600
            and fixture["nonzero_reciprocal_rows"] == 5
            and fixture["scalar_abs_called"] is False
            and fixture["scalar_endpoint_negation"] is False,
            "New native-block fixture changed")
    reciprocal = packet["reciprocal_fixture"]
    require(reciprocal["passed"] is True and reciprocal["checked_cases"] == 6
            and len(reciprocal["fixture_rows"]) == 6
            and reciprocal["huge_rational_legacy_underbound_cases"] == 3
            and reciprocal["global_scalar_precision_restored"] is True
            and reciprocal["old_helper_called_only_on_fresh_inputs"] is True,
            "Fresh reciprocal witness scope changed")
    underbounds = 0
    for row in reciprocal["fixture_rows"]:
        pair = row["certified_lower_log_exact_tuples"]
        require(len(pair) == 2, "Malformed certified fixture interval")
        lo, hi = (finite_tuple(item) for item in pair)
        require(mpf_cmp(lo, hi) <= 0, "Reversed certified fixture interval")
        flipped = (1-lo[0], *lo[1:]) if lo[1] else (0, *lo[1:])
        corrected = row["corrected_inverse_exact_tuples"]
        require(len(corrected) == 2
                and finite_tuple(row["expected_inverse_lower_exact_tuple"]) == flipped
                and all(finite_tuple(v) == flipped for v in corrected),
                "Reciprocal tuple witness changed")
        strict = mpf_cmp(finite_tuple(row["legacy_inverse_upper_exact_tuple"]), flipped) < 0
        require(strict == row["legacy_upper_lt_expected_lower"], "False reciprocal underbound witness")
        underbounds += strict
    require(underbounds == reciprocal["old_legacy_underbound_cases"], "Witness count changed")
    limits = packet["limits"]
    require(limits["restoration_post_inverse_Utheta_repaired"] is True
            and all(limits[k] is True for k in ("local_blocks_only",
                    "long_reference_inverse_Utheta_not_repaired", "native_full_p2_not_rebuilt",
                    "long_reference_full_caps_not_recertified", "restoration_full_caps_not_recertified",
                    "mixed_loop_and_Rm_not_recertified", "finite_N_not_recertified"))
            and limits["full_goal_complete"] is False,
            "Local factor repair promoted to a full certificate")
    require(packet["primitive_api"]["scalar_abs_endpoint_rounding"] is False
            and packet["primitive_api"]["signed_function_value_selected"] is False
            and packet["reciprocal_helper_api"]["stored_logarithms_exponentiated"] is False
            and packet["reciprocal_helper_api"]["lower_log_endpoint_tuple_flipped_in_new_helper"] is True,
            "Primitive execution scope changed")
    preservation = packet["preservation"]
    require(all(preservation[k] is True for k in (
                "inherited_historical_inverse_L_caps_retained",
                "inherited_log_caps_not_treated_as_signed_values",
                "signed_interval_data_not_deleted_or_midpointed"))
            and all(preservation[k] is False for k in (
                "old_native_modules_edited", "accepted_builders_or_producers_called",
                "accepted_checkers_or_self_checks_called", "accepted_recurrences_or_quadratures_called")),
            "Preservation/execution scope changed")
    checks, audit = packet["checks"], packet["execution_audit"]
    require(checks["passed"] is True
            and all(checks[k] is True for k in ("source_pins_checked", "source_guards_checked",
                    "signed_U_M_leaf_checkpoints_built", "directed_L_rows_built",
                    "directed_inverse_L_rows_built", "directed_restoration_post_inverse_Utheta_rows_built",
                    "inverse_Utheta_formula_guarded"))
            and all(checks[k] is False for k in ("inherited_caps_recertified",
                    "downstream_finite_N_recertified", "downstream_Rm_recertified", "source_executed"))
            and all(audit[k] is False for k in ("source_executed", "accepted_ancestors_replayed",
                    "accepted_builders_or_owners_replayed", "accepted_checkers_replayed",
                    "accepted_recurrences_replayed", "accepted_quadratures_replayed"))
            and audit["legacy_generic_reciprocal_called_only_on_six_fresh_inputs"] is True,
            "Corrective execution audit changed")


def validate(report, accepted):
    require(set(report) == set(accepted) | {KEY}, "Unexpected report fields")
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes",
               "implementation_hashes", "contract", "checks", "execution_seconds",
               "execution_audit"}
    for key, value in accepted.items():
        if key not in mutable:
            require(report.get(key) == value, "Historical mathematical data changed: " + key)
    require(report["parent_report"] == PARENT_NAME
            and report["parent_report_sha256"] == PARENT_SHA256,
            "Corrective predecessor changed")
    validate_packet(report[KEY], accepted)
    require(report["checks"]["passed"] is True
            and report["checks"]["scope"] == "named_directed_corrective_blocks_only"
            and report["checks"]["inherited_caps_revalidated"] is False,
            "Report check scope changed")
    count = record_count(report[KEY])
    require(count > 0 and report["checks"]["directed_records"] == count,
            "Directed corrective evidence count changed")
    contract = report["contract"]
    require(contract["inherited_log_cap_certification_verified"] is False
            and contract["inherited_log_caps_require_revalidation"] is True
            and contract["new_Rm_pipeline_accepted"] is False
            and contract["global_five_moment_closure"] is False
            and contract["full_goal_complete"] is False,
            "Scoped repair promoted to a complete certificate")
    inputs, implementations = merge_pins(accepted, report[KEY])
    require(report["input_hashes"] == inputs and report["implementation_hashes"] == implementations,
            "Corrective report pin set changed")
    return count


def produce(accepted):
    reviewed_check()
    before = check_pins(accepted)
    started = time.monotonic()
    sys.path.insert(0, str(ROOT / "src"))
    module = importlib.import_module("openai_ns_reconstruction.original_native_control_directed_repair")
    packet = module.build(ROOT, accepted)
    report = copy.deepcopy(accepted)
    report[KEY] = packet
    report["stage"] = "current_original_native_control_directed_repair"
    report["parent_report"] = PARENT_NAME
    report["parent_report_sha256"] = PARENT_SHA256
    report["input_hashes"], report["implementation_hashes"] = merge_pins(accepted, packet)
    report["contract"].update(
        inherited_log_cap_certification_verified=False,
        inherited_log_caps_require_revalidation=True,
        new_Rm_pipeline_accepted=False,
        global_five_moment_closure=False,
        full_goal_complete=False,
    )
    report["checks"] = dict(passed=True, scope="named_directed_corrective_blocks_only",
                            inherited_caps_revalidated=False,
                            directed_records=record_count(packet),
                            accepted_ancestors_replayed=False)
    report["execution_audit"] = dict(
        new_builder="original_native_control_directed_repair.build",
        accepted_construction_builders_called=False,
        accepted_generic_self_checks_called=False,
        historical_mathematical_data_preserved=True,
        old_Rm_drafts_executed=False,
        input_pin_count=before,
    )
    report["execution_seconds"] = time.monotonic() - started
    validate(report, accepted)
    check_pins(report)
    return report


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-report", action="store_true")
    mode.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    accepted = parent()
    path = BASE / NAME
    if args.write_report:
        require(not path.exists(), "Corrective report already exists; production replay refused")
        report = produce(accepted)
        path.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                        encoding="utf-8", newline="\n")
    else:
        reviewed_check()
        report = load(path)
    records = validate(report, accepted)
    pins = check_pins(report)
    print(json.dumps(dict(passed=True, report=str(path), bytes=path.stat().st_size,
                          sha256=sha(path), pins=pins, directed_records=records,
                          scope="named_directed_corrective_blocks_only",
                          inherited_caps_revalidated=False,
                          accepted_ancestors_replayed=False), ensure_ascii=False))


if __name__ == "__main__":
    main()
