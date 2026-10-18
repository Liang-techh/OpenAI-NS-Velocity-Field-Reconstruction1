"""Save or inspect same-source six-I1 C1 defect and pressure-error caps.

This runner preserves the accepted physical U/Q increment packet and checks
the new same-source defect arithmetic, physical R^2*hb integration weights,
and the open outer R_Q_y gate.  It does not select signed fields or execute
an accepted producer, owner, norm, recurrence, fixture, or quadrature.
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "experiments/root_st073"
sys.path.insert(0, str(ROOT / "src"))

from openai_ns_reconstruction import original_micro_I1_defect_bounds as current


build = current.build
PARENT = current.PARENT
PARENT_SHA256 = current.PARENT_SHA256
DEFECT = current.DEFECT
CORRELATED = current.CORRELATED
BRIDGE = current.BRIDGE
LOOP = current.LOOP
LOOP_ARCHIVE = LOOP[:-3] + ".json.gz"
Q_HELPER = "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"
NAMES = tuple(current.NAMES)
require = current.require
sha = current.sha

REVIEWED_MODULE_SHA256 = "cc6df5d1694784e0c41c83c7a9f6771c38fe43db84104b7f217e7af74f6b6cb1"

NAME = "lei_ren_part1_paper_compliant_current_original_micro_I1_defect_bounds.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_I1_defect_bounds.py",
    "scripts/run_original_micro_I1_defect_bounds.py",
)


def _load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted Q-increment parent changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True, "Checked Q-increment parent required")
    return report


def pins(report):
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        for name, expected in report[key].items():
            require(sha(base / name) == expected, "Changed source pin: " + name)
            count += 1
    return count


def _passed(value, label):
    require(isinstance(value, dict) and value.get("passed") is True, label)


def _directed(value, label):
    require(isinstance(value, dict) and value.get("directed") is True,
            label + " directed record missing")


def _rows(rows, count, label):
    require(isinstance(rows, list) and len(rows) == count,
            label + " row count changed")
    orders = []
    for row in rows:
        require(isinstance(row, dict), label + " row is not an object")
        orders.append(row.get("ordinary_Z_order"))
        _directed(row.get("log_absolute_upper_bound"), label + " log cap")
    require(orders == list(range(count)), label + " ordinary-Z order changed")


def _records(rows, count, label):
    require(isinstance(rows, list) and len(rows) == count,
            label + " record count changed")
    for row in rows:
        _directed(row, label + " directed coefficient")


def _ratio_rows(rows, count, label):
    require(isinstance(rows, list) and len(rows) == count,
            label + " ratio row count changed")
    orders = []
    for row in rows:
        require(isinstance(row, dict), label + " ratio row is not an object")
        orders.append(row.get("ordinary_Z_order"))
        require("log_cap_over_selected_logC" in row
                and "log_absolute_upper_bound" not in row,
                label + " explicit ratio field missing")
        _directed(row["log_cap_over_selected_logC"], label + " ratio")
    require(orders == list(range(count)), label + " ratio ordinary-Z order changed")


def _ast(receipt, label, path, function, snippets):
    _passed(receipt, label + " source AST receipt missing")
    require(receipt.get("path") == path
            and receipt.get("function") == function
            and receipt.get("source_executed") is False,
            label + " source AST binding changed")
    required = receipt.get("required_text", ())
    if isinstance(required, dict):
        required = required.keys()
    require(set(required) == set(snippets),
            label + " source AST fragments changed")
    require(isinstance(receipt.get("AST_sha256"), str),
            label + " source AST digest missing")


def _source_contract(packet):
    contract = packet.get("source_contract")
    _passed(contract, "I1 source contract missing")
    require(
        contract.get("original_physical_I1_measure") == "R^2*hb*ds exactly once"
        and contract.get("accepted_free_Z_internal_B_and_RQy_seam_consumed") is True
        and contract.get("first_lower_finite_N_error_exact_zero_by_source_flat_collar") is True
        and contract.get("same_native_background_on_both_micro_charts") is True
        and contract.get("outer_RQy_retained") is True
        and contract.get("old_212_and_all_errors_unchanged") is True
        and contract.get("source_graph_functions_not_cap_values") is True
        and contract.get("accepted_seam_checker_replayed") is False,
        "I1 source contract changed",
    )
    handles = contract.get("six_signed_defect_handles")
    require(isinstance(handles, dict)
            and set(handles) == {"first_micro", "second_micro"},
            "Six signed-defect chart handles changed")
    for chart in ("first_micro", "second_micro"):
        values = handles[chart]
        require(isinstance(values, dict) and set(values) == set(NAMES),
                "Six signed-defect handle fields changed: " + chart)
        for name in NAMES:
            row = values[name]
            require(isinstance(row, list) and len(row) == 2
                    and all(isinstance(value, int)
                            and not isinstance(value, bool) and value >= 0
                            for value in row),
                    "C1 signed-defect handles changed: " + chart + "." + name)


def _background(packet):
    bg = packet.get("background_native_norm_provider")
    _passed(bg, "Native background norm provider missing")
    for key, count in (
        ("physical_native_U", 3),
        ("physical_native_mean", 4),
        ("physical_native_Q", 3),
    ):
        _rows(bg.get(key), count, "Native " + key)
    proof = bg.get("Q_recovery")
    _passed(proof, "Native background Q recovery proof missing")
    for key in ("positive_L", "L_Z", "L_ZZ"):
        _directed(proof.get(key), "Native Q recovery " + key)
    _rows(proof.get("numerator_log_caps"), 3, "Native Q recovery numerator")
    require(proof.get("ordinary_derivatives") is True
            and proof.get("no_top_row_padding") is True,
            "Native Q recovery controls changed")
    require(
        bg.get("inputs_are_saved_actual_native_bounds") is True
        and bg.get("saved_U_M_already_physical") is True
        and bg.get("inverse_Pstar_in_C6_restored_once_by_Q4") is True
        and bg.get("no_extra_Pstar_on_saved_U_M") is True
        and bg.get("native_mean_is_M_not_R_times_M") is True
        and bg.get("native_not_identified_with_comparison") is True
        and bg.get("new_background_Q_norm_computation") is True
        and bg.get("accepted_U_M_producer_not_replayed") is True,
        "Native U/M physical-scale provenance changed",
    )


def _parent_provenance(accepted):
    source = accepted.get("actual_micro_b_p2_mixed_bounds")
    require(isinstance(source, dict) and source.get("passed") is True,
            "Accepted native b/p2 parent missing")
    native = source.get("native_angular_and_axial_provider")
    require(
        isinstance(native, dict)
        and native.get("passed") is True
        and native.get("ordinary_Z0_Z4") is True
        and native.get("physical_y0_y1") is True
        and native.get("original_core_U_inlet_and_all_corrections_retained") is True
        and native.get("native_and_comparison_never_identified") is True
        and native.get("physical_dy_hb_ds_once") is True,
        "Saved native physical source proof changed",
    )
    increment = accepted.get("actual_micro_Q_increment_bounds")
    require(isinstance(increment, dict) and increment.get("passed") is True,
            "Accepted Q-increment parent missing")
    radius = increment.get("radius_provenance")
    require(
        isinstance(radius, dict)
        and radius.get("passed") is True
        and radius.get("R_equals_epsilon_rho") is True
        and radius.get("epsilon_Lambda_equals_one") is True
        and radius.get("physical_R_Z_independent") is True
        and radius.get("actual_micro_radius_handle") == 68644,
        "Accepted physical R provenance changed",
    )


def _weight_provenance(packet):
    weight = packet.get("weight_provenance")
    require(isinstance(weight, dict), "Physical I1 weight provenance missing")
    require(
        set(weight) == {
            "passed",
            "physical_measure_handle",
            "physical_measure_is_original_hb_not_integer_eight",
            "loghb_exact_tuples_equal_accepted_source_basis",
            "same_actual_radius_handle",
            "native_chart_widths_at_most_one",
            "actual_second_upper_number",
            "source_first_lower_is_sc_half",
            "positive_sc_upper_receipt_consumed",
            "Z_independent_measure_and_native_endpoints",
            "source_ledger_sha256",
            "no_geometry_owner_or_hb_producer_executed",
        },
        "Physical I1 weight provenance fields changed",
    )
    require(
        weight.get("passed") is True
        and weight.get("physical_measure_handle") == 8
        and weight.get("physical_measure_is_original_hb_not_integer_eight") is True
        and weight.get("loghb_exact_tuples_equal_accepted_source_basis") is True
        and weight.get("same_actual_radius_handle") == 68644
        and weight.get("native_chart_widths_at_most_one") is True
        and weight.get("actual_second_upper_number") == 2
        and weight.get("source_first_lower_is_sc_half") is True
        and weight.get("positive_sc_upper_receipt_consumed") is True
        and weight.get("Z_independent_measure_and_native_endpoints") is True
        and isinstance(weight.get("source_ledger_sha256"), str)
        and weight.get("no_geometry_owner_or_hb_producer_executed") is True,
        "Physical I1 weight provenance changed",
    )


def _coefficients(packet):
    coeff = packet.get("coefficient_ordinary_Z0_Z1")
    require(isinstance(coeff, dict)
            and set(coeff) == {"volume", "boundary"},
            "Signed C1 coefficient packet changed")
    require(set(coeff["volume"]) == set(NAMES),
            "Six volume coefficient fields changed")
    require(set(coeff["boundary"]) == {"q", "q_squared", "U_squared"},
            "Three boundary coefficient fields changed")
    for key, values in coeff["volume"].items():
        _records(values, 2, "Volume coefficient " + key)
    for key, values in coeff["boundary"].items():
        _records(values, 2, "Boundary coefficient " + key)


def _chart(chart, name):
    require(isinstance(chart, dict), "Actual " + name + " I1 chart missing")
    require(set(chart) == {"six_defect_log_caps", "six_local_I1_log_caps"},
            name + " I1 chart fields changed")
    for field in ("six_defect_log_caps", "six_local_I1_log_caps"):
        caps = chart.get(field)
        require(isinstance(caps, dict) and set(caps) == set(NAMES),
                name + " " + field + " fields changed")
        for key in NAMES:
            _rows(caps[key], 2, name + " " + field + "." + key)


def _pressure(packet):
    fields = {
        "reduced_volume_D",
        "outer_reduced_boundary_B",
        "half_S1_pressure_without_RQy",
        "retained_outer_RQy",
        "complete_micro_pressure_error",
    }
    caps = packet.get("pressure_error_log_caps")
    require(isinstance(caps, dict) and set(caps) == fields,
            "Five pressure error cap fields changed")
    for key in fields:
        _rows(caps[key], 2, "Pressure error " + key)
    ratios = packet.get("pressure_error_log_cap_ratios")
    require(isinstance(ratios, dict) and set(ratios) == fields,
            "Five pressure error ratio fields changed")
    for key in fields:
        _ratio_rows(ratios[key], 2, "Pressure error " + key)
    strict = packet.get("strict_log_envelopes")
    _passed(strict, "Strict I1 pressure envelope missing")
    require(strict.get("six_defects") == -2600
            and strict.get("joint_I1_and_reduced_volume_and_half_S1") == -2700
            and strict.get("reference") == "selected positive logC",
            "Strict I1 pressure envelope changed")
    require(packet.get("physical_I1_error_C1_below_2_power_minus1024") is True,
            "Physical I1 threshold missing")
    small = packet.get("full_micro_pressure_error_small")
    require(isinstance(small, list) and len(small) == 2,
            "Full pressure threshold rows missing")
    require([row.get("ordinary_Z_order") for row in small] == [0, 1]
            and all(row.get("below_2_power_minus1024") is False for row in small),
            "Outer R_Q_y gate must remain false at both C1 orders")


def _new_controls(packet):
    controls = packet.get("new_six_defect_controls")
    _passed(controls, "New six-defect controls missing")
    require(
        controls.get("six_polynomial_differences_and_ordinary_Z1") is True
        and controls.get("square_and_cross_terms_retained") is True
        and controls.get("ordinary_product_binomials") is True
        and controls.get("no_extra_Z_on_U_QZ_coefficient") is True
        and controls.get("U_squared_coefficient")
        == "2*(1-delta)*Z^2/L^2"
        and controls.get("U_UZ_coefficient")
        == "2*Z*(1-Z^2)/L^2"
        and controls.get("boundary_U_squared_coefficient") == "-2*Z^2/L^2"
        and controls.get("actual_native_upper_is_two_not_graph_handle_48") is True
        and controls.get("Pstar_in_physical_fields_not_added_to_I1_weight") is True,
        "Six-defect coefficient controls changed",
    )


def _packet_fields(report, accepted):
    require(report["parent_report"] == Path(PARENT).name
            and report["parent_report_sha256"] == PARENT_SHA256,
            "Same Q-increment predecessor required")
    mutable = {
        "stage",
        "parent_report",
        "parent_report_sha256",
        "input_hashes",
        "implementation_hashes",
        "contract",
        "checks",
        "execution_seconds",
    }
    for key, value in accepted.items():
        if key not in mutable:
            require(report[key] == value,
                    "Inherited signed packet or error changed: " + key)

    packet = report.get("actual_micro_I1_defect_bounds")
    _passed(packet, "Actual I1 defect packet missing")
    require(packet.get("source_family") == accepted["source_family"],
            "Source family changed")
    _source_contract(packet)
    _parent_provenance(accepted)
    _background(packet)
    _weight_provenance(packet)
    _new_controls(packet)
    _coefficients(packet)

    ast = packet.get("source_AST")
    require(isinstance(ast, dict)
            and set(ast) == {"defects", "native_scale", "shared_background", "joined"},
            "I1 source AST field set changed")
    _ast(ast["defects"], "Signed defects", DEFECT, "build_defect", (
        "du = _jet_scale(g, delta_V[:3], Pstar)",
        "dm = _jet_scale(g, history_m[:4], Pstar)",
        "U0 = _jet_sub(g, physical_U[:3], du)",
        "Q0 = _jet_sub(g, physical_Q[:3], dq)",
    ))
    _ast(ast["native_scale"], "Native C6 scale", BRIDGE + ".py", "generic_fields", (
        "inverse=g.quotient(g.one,Pstar,'same strictly positive Pstar')",
        "V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)",
    ))
    _ast(ast["shared_background"], "Shared loop background", LOOP, "__init__", (
        "'first_micro':bridge.micro,'second_micro':bridge.micro",
        "fields={key:a.trim(jet) for key,jet in chart['fields'].items()}",
        "for key in RATES:fields[key]=a.add(a.trim(chart['fields'][key]),history[key])",
    ))
    _ast(ast["joined"], "Joined six-I1 reduction", CORRELATED, "build", (
        'radius=chart["radius"];weight=fixed(g.mul(radius,radius,chart["physical_measure"]))',
        'S1=add(B,D);increment=add(half(S1),neg(boundary("R_Q_y")))',
        'neg(half(integrals["q_squared"])),divide(mul(d,integrals["U_q_Z"]))',
        'g.pscale(divide(mul(h,z,z,integrals["U_squared"]),L2),two)',
        'g.pscale(divide(mul(z,d,integrals["U_U_Z"]),L2),two)',
    ))

    require(packet.get("domain")
            == accepted["actual_micro_Q_increment_bounds"]["domain"],
            "I1 domain changed")
    charts = packet.get("charts")
    require(isinstance(charts, dict) and set(charts) == {"first", "second"},
            "Both I1 charts required")
    _chart(charts["first"], "first")
    _chart(charts["second"], "second")
    joint = packet.get("six_joint_I1_log_caps")
    require(isinstance(joint, dict) and set(joint) == set(NAMES),
            "Joined six-I1 fields changed")
    for key in NAMES:
        _rows(joint[key], 2, "Joined I1 " + key)

    _pressure(packet)
    require(packet.get("physical_R_squared_and_hb_applied_once") is True
            and packet.get("signed_function_values_evaluated") is False
            and packet.get("signed_six_integrals_evaluated") is False
            and packet.get("outer_RQy_certified_small") is False
            and packet.get("old_signed_sources_and_all_error_caps_unchanged") is True
            and packet.get("accepted_ancestors_replayed") is False
            and packet.get("global_five_moment_closure") is False
            and packet.get("full_goal_complete") is False,
            "I1 completion gates changed")

    source_pins = packet.get("source_pins")
    require(isinstance(source_pins, dict), "Current I1 source pins missing")
    expected_pins = (
        (PARENT, PARENT_SHA256),
        (DEFECT, accepted["implementation_hashes"][DEFECT]),
        (CORRELATED, accepted["implementation_hashes"][CORRELATED]),
        (BRIDGE + ".py", accepted["input_hashes"][Path(BRIDGE + ".py").name]),
        (LOOP, accepted["input_hashes"][Path(LOOP).name]),
        (LOOP_ARCHIVE, accepted["input_hashes"][Path(LOOP_ARCHIVE).name]),
        (Q_HELPER, accepted["implementation_hashes"][Q_HELPER]),
    )
    require(set(source_pins) == {path for path, _ in expected_pins},
            "I1 source pin set changed")
    for path, expected in expected_pins:
        row = source_pins.get(path)
        require(isinstance(row, dict)
                and row.get("path") == path
                and row.get("expected_sha256") == expected,
                "Reviewed I1 source pin missing: " + path)

    return dict(
        passed=True,
        actual_micro_six_I1_C1_defect_and_pressure_log_caps=True,
        native_background_U3_mean4_Q3=True,
        two_micro_six_defect_and_local_I1_C1=True,
        joined_six_I1_C1=True,
        strict_defect_I1_D_halfS1_envelopes=True,
        physical_I1_threshold_passed=True,
        full_pressure_outer_RQy_gate_retained_false=True,
        signed_values_integrals_globalclosure_open=True,
        old_signed_sources_and_all_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(REVIEWED_MODULE_SHA256 != "REVIEW_PENDING",
            "Independent review hash pending")
    require(sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256,
            "Only the independently reviewed I1 module may produce this report")
    accepted = parent()
    pins(accepted)
    packet = build(ROOT, accepted)
    require(accepted == parent(), "New enclosure mutated its accepted parent")

    inputs = dict(accepted["input_hashes"])
    implementation = dict(accepted["implementation_hashes"])
    for row in packet["source_pins"].values():
        path = row["path"]
        is_input = path.startswith("experiments/root_st073/")
        name = Path(path).name if is_input else path
        target = inputs if is_input else implementation
        require(sha(ROOT / path) == row["expected_sha256"],
                "Changed new input: " + path)
        require(target.get(name, row["expected_sha256"]) == row["expected_sha256"],
                "Conflicting source pin: " + path)
        target[name] = row["expected_sha256"]

    implementation.update({path: sha(ROOT / path) for path in CODE})
    report = copy.deepcopy(accepted)
    report.update(
        stage="actual_micro_six_I1_C1_defect_bounds",
        parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_I1_defect_bounds=packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
    )
    report["contract"] = dict(
        actual_micro_six_I1_C1_defect_and_pressure_log_caps=True,
        native_background_U3_mean4_Q3=True,
        two_micro_six_defect_and_local_I1_C1=True,
        joined_six_I1_C1=True,
        strict_defect_I1_D_halfS1_envelopes=True,
        physical_I1_threshold_passed=True,
        full_pressure_outer_RQy_gate_retained_false=True,
        signed_values_integrals_globalclosure_open=True,
        old_signed_sources_and_all_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Actual same-source six signed finite-N defect C1 pairs, native "
            "physical U/M/Q background rows, R^2*hb local and joined I1 caps, "
            "signed C1 reduced coefficients, and D/half-S1 pressure caps; "
            "the outer R_Q_y threshold, signed values, full pressure closure, "
            "and global five-moment closure remain open"
        ),
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independent_runtime_worker_model_metadata_available=False,
        accepted_ancestors_replayed=False,
    )
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path):
    report = _load(path)
    accepted = parent()
    checks = report["checks"]
    require(
        checks["passed"]
        and checks.get("independent_new_math_review_completed") is True
        and checks.get("independent_review_completed") is True
        and checks["independently_reviewed_module_sha256"]
        == report["implementation_hashes"][CODE[0]]
        == REVIEWED_MODULE_SHA256,
        "Same independently reviewed I1 implementation required",
    )
    return dict(
        passed=True,
        pins=pins(report),
        packet_invariants=_packet_fields(report, accepted),
        serialized_only=True,
        accepted_owner_norm_recurrence_or_quadrature_replayed=False,
    )


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
        temporary.write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        temporary.replace(path)
        print(json.dumps(dict(
            passed=True,
            execution_seconds=report["execution_seconds"],
            bytes=path.stat().st_size,
            sha256=sha(path),
        )))
    else:
        print(json.dumps(inspect(path), indent=2))
