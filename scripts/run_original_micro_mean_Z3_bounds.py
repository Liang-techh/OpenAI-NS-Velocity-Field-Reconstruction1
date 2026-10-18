"""Save or inspect source-owned ordinary-Z3 phase-mean logarithmic caps.

This runner preserves the accepted fixed-phase inverse packet and records the
same-source C/Cy primitives, moving histories, inherited second-chart memory,
and complete physical mean caps.  It leaves velocity finite-N errors, QZZ,
full C5 containment, and global closure open.
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

# Keep the source module as one object so this checker follows its explicit
# source-provenance contract without replaying an accepted producer.
from openai_ns_reconstruction import original_micro_mean_Z3_bounds as current


build = current.build
PARENT = current.PARENT
PARENT_SHA256 = current.PARENT_SHA256
KERNEL = current.KERNEL
KERNEL_SHA256 = current.KERNEL_SHA256
require = current.require
sha = current.sha

# The parent replaces this placeholder after independent review freezes the
# source module.  A report cannot be written while that review is pending.
REVIEWED_MODULE_SHA256 = "f7892b4134fa81ecaf3b4348c41a58c43d70ef9605ab681dcf16b8d9fd42bf54"

NAME = "lei_ren_part1_paper_compliant_current_original_micro_mean_Z3_bounds.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_mean_Z3_bounds.py",
    "scripts/run_original_micro_mean_Z3_bounds.py",
    KERNEL,
)


def _load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted phase-inverse parent changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True, "Checked phase-inverse parent required")
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


def _directed(record, label):
    require(isinstance(record, dict) and record.get("directed") is True,
            label + " directed record missing")


def _rows(rows, label):
    """Check four ordinary-Z rows without replaying any numerical cap."""
    require(isinstance(rows, list) and len(rows) == 4,
            label + " must have four ordinary-Z rows")
    orders = []
    for row in rows:
        require(isinstance(row, dict), label + " row is not an object")
        orders.append(row.get("ordinary_Z_order"))
        _directed(row.get("log_absolute_upper_bound"), label + " log cap")
    require(orders == [0, 1, 2, 3], label + " must cover ordinary Z0..Z3")


def _ratio_rows(rows, label):
    """Check the explicit selected-logC ratio field."""
    require(isinstance(rows, list) and len(rows) == 4,
            label + " must have four ratio rows")
    orders = []
    for row in rows:
        require(isinstance(row, dict), label + " ratio row is not an object")
        orders.append(row.get("ordinary_Z_order"))
        require("log_cap_over_selected_logC" in row
                and "log_absolute_upper_bound" not in row,
                label + " ratio field changed")
        _directed(row["log_cap_over_selected_logC"], label + " ratio")
    require(orders == [0, 1, 2, 3], label + " must cover ordinary Z0..Z3")


def _ast(receipt, label, function, snippets):
    _passed(receipt, label + " source AST receipt missing")
    require(receipt.get("path") == KERNEL
            and receipt.get("function") == function
            and receipt.get("source_executed") is False,
            label + " source AST binding changed")
    required = receipt.get("required_text", ())
    if isinstance(required, dict):
        required = required.keys()
    require(set(required) == set(snippets), label + " source AST fragments changed")
    require(isinstance(receipt.get("AST_sha256"), str),
            label + " source AST digest missing")


def _source_contract(packet):
    contract = packet.get("source_contract")
    _passed(contract, "Source mean contract missing")
    for key in (
        "accepted_source_proof_consumed",
        "exact_transport_identity_consumed",
        "ordinary_Z_seam_identity_consumed",
        "first_inlet_zero_only_at_source_owned_left_collar",
        "second_history_not_reset",
        "actual_phase_Z_independent",
        "ordinary_Z_derivatives_have_no_factorials",
        "source_phase_values_not_chosen",
        "source_or_checker_not_executed",
    ):
        require(contract.get(key) is True, "Source mean contract changed: " + key)
    require(contract.get("second_memory_formula") == "exp(-hb*(s-1))",
            "Joined second-chart memory formula changed")
    handles = contract.get("source_graph_handles")
    require(isinstance(handles, dict), "Source graph handles missing")
    fields = {"C", "Cy", "JC", "JC_joint", "JCy", "JCy_joint",
              "complete_physical_delta_m"}
    require(set(handles) == fields, "Source graph handle field set changed")
    for field, charts in handles.items():
        require(isinstance(charts, dict)
                and set(charts) == {"first_micro", "second_micro"},
                "Source graph chart handles changed: " + field)
        for chart, rows in charts.items():
            require(isinstance(rows, (list, tuple)) and len(rows) == 4,
                    "Source graph Z3 handles missing: " + field + "." + chart)
            require(all(isinstance(value, int) and not isinstance(value, bool)
                        and value >= 0 for value in rows),
                    "Source graph handles are not graph IDs: " + field + "." + chart)


def _history_controls(packet):
    controls = packet.get("new_history_controls")
    _passed(controls, "New joined-history controls missing")
    require(
        controls.get("positive_kernel_mass") == "1-exp(-hb*length)"
        and controls.get("first_native_width_at_most_one") is True
        and controls.get("second_native_width_at_most_one") is True
        and controls.get("first_mass_upper") == "hb"
        and controls.get("joined_mass_upper") == "2*hb"
        and controls.get("physical_hb_applied_once") is True
        and controls.get("second_memory_times_first_kernel_equals_joined_kernel") is True
        and controls.get("mean_formula") ==
        "delta_m=Pstar/N^2*(C-JC_joint-JCy_joint)"
        and controls.get("total_y_identity") ==
        "delta_m_total_y+delta_m=Pstar*B_over_Pstar/N"
        and controls.get("norm_bound_does_not_drop_Cy_or_incoming_memory") is True
        and controls.get("ordinary_Z_commutes_with_phase_and_moving_integrals") is True
        and controls.get("periodic_phase_cuts") ==
        "C(0)=C(1)=0 and derivatives of this parameter identity; Cy=C_slow_y"
        and controls.get("threshold_is_local_mean_error_not_full_NS_residual") is True,
        "Joined-history controls changed",
    )


def _chart(chart, name):
    _passed(chart, "Actual " + name + " mean chart missing")
    require(chart.get("ordinary_Z_orders") == [0, 1, 2, 3],
            name + " ordinary Z0..Z3 required")
    require(chart.get("source_chart") == name + "_micro"
            and chart.get("C_and_Cy_not_reset") is True
            and chart.get("second_incoming_history_retained") is (name == "second"),
            name + " source-chart history controls changed")
    fields = {
        "physical_Pstar_C",
        "physical_Pstar_Cy",
        "physical_Pstar_JC_joint",
        "physical_Pstar_JCy_joint",
        "known_C_mean_part",
        "retained_Cy_mean_part",
        "physical_complete_delta_m",
        "physical_total_y_delta_m",
    }
    caps = chart.get("log_caps")
    require(isinstance(caps, dict) and set(caps) == fields,
            name + " eight mean log-cap fields changed")
    for key, value in caps.items():
        _rows(value, name + " " + key)
    source_caps = chart.get("source_physical_B_log_caps")
    require(isinstance(source_caps, dict) and set(source_caps) == {"B", "slow_By"},
            name + " source physical B/slow-By caps missing")
    for key, value in source_caps.items():
        _rows(value, name + " source " + key)
    require(chart.get("primitive_caps_inherited_from_B_by_length_one_not_function_equality") is True,
            name + " primitive cap inheritance control missing")
    ratios = chart.get("log_cap_ratios")
    require(isinstance(ratios, dict) and set(ratios) == fields,
            name + " eight ratio fields changed")
    for key, value in ratios.items():
        _ratio_rows(value, name + " " + key)
    require(chart.get("ratios_are_logcaps_divided_by_selected_positive_logC") is True,
            name + " selected-logC ratio control missing")
    require(chart.get("complete_physical_mean_Z0_Z3_below_2_power_minus1024") is True
            and chart.get("complete_physical_total_y_mean_Z0_Z3_below_2_power_minus1024") is True,
            name + " local physical mean thresholds missing")
    strict = chart.get("strict_log_envelope")
    _passed(strict, name + " strict logarithmic mean envelope missing")
    require(strict.get("mean_multiplier") == -4500
            and strict.get("total_y_mean_multiplier") == -600
            and strict.get("reference") == "selected positive logC",
            name + " strict logarithmic mean envelope changed")


def _packet_fields(report, accepted):
    require(
        report["parent_report"] == Path(PARENT).name
        and report["parent_report_sha256"] == PARENT_SHA256,
        "Same signed phase-inverse predecessor required",
    )
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

    packet = report.get("actual_micro_mean_Z3_bounds")
    _passed(packet, "Actual mean Z3 packet missing")
    require(packet.get("source_family") == accepted["source_family"],
            "Source family changed")

    _source_contract(packet)
    _history_controls(packet)

    ast = packet.get("source_AST")
    require(isinstance(ast, dict) and set(ast) == {"phase", "moving", "joined"},
            "Mean source AST field set changed")
    _ast(ast["phase"], "Phase primitive", "_phase_rows", (
        '"fixed_phase_slow_y": source_name == "Cy"',
        '"total_y_chain_excluded": source_name == "Cy"',
        '"phase_primitive_is_direct_ordinary_Z_row": True',
        '"source_operation": "phase_integral_of_slow_B_y" if source_name == "Cy" else "phase_integral_of_B"',
    ))
    _ast(ast["moving"], "Moving convolution", "_moving_rows", (
        "exponent = graph.neg(graph.mul(hb, distance))",
        "integrand = graph.mul(kernel, primitive_at_t, hb)",
        '"physical_hb_weight_applied_once": True',
    ))
    _ast(ast["joined"], "Joined mean", "build_micro_mean", (
        'scale = g.div(Pstar, n_squared, "same exact positive common N squared")',
        'formula="exp(-hb*(s-1))"',
        "_scale_rows(g, first_JC_at_one, second_memory)",
        "_scale_rows(g, first_JCy_at_one, second_memory)",
        "second_complete = _add_rows(g, second_known, second_unresolved)",
    ))

    require(packet.get("domain") == accepted["actual_micro_phase_inverse_bounds"]["domain"]
            and packet.get("ordinary_Z_orders") == [0, 1, 2, 3],
            "Mean Z3 domain changed")
    charts = packet.get("charts")
    require(isinstance(charts, dict) and set(charts) == {"first", "second"},
            "Both mean charts required")
    _chart(charts["first"], "first")
    _chart(charts["second"], "second")

    scales = packet.get("physical_scale_logs")
    require(isinstance(scales, dict)
            and set(scales) == {"Pstar_E_upper", "hb", "N", "selected_logC"},
            "Physical Pstar/E/hb/N/logC scales changed")
    for key, value in scales.items():
        _directed(value, "Physical scale " + key)

    amplitude = packet.get("amplitude_provenance")
    _passed(amplitude, "Actual amplitude provenance missing")
    require(
        amplitude.get("actual_positive_E0_is_the_same_native_E") is True
        and amplitude.get("frozen_E0_at_every_source_point_below_saved_uniform_E_upper") is True
        and amplitude.get("domain") == {"s": ["0", "2"], "Z": ["37/100", "93/250"]}
        and amplitude.get("bound")
        == "Pstar*E=F0*sqrt(2R)*phi <= sqrt(2*rho_upper/Lambda)*Phi_core_upper"
        and amplitude.get("F0_positive_and_at_most_one") is True
        and amplitude.get("native_phi_positive_and_nonincreasing") is True
        and amplitude.get("actual_core_argument_domain_receipt_consumed") is True
        and amplitude.get("saved_global_log_bound_consumed") is True,
        "Actual amplitude provenance controls changed",
    )

    for key in (
        "actual_Cy_Z3_norm_provider",
        "actual_complete_mean_Z3_norm_provider",
        "actual_complete_total_y_mean_Z3_norm_provider",
        "primitive_phase_integration_length_bound_one",
        "physical_Pstar_applied_once_and_common_N_squared_retained",
    ):
        require(packet.get(key) is True, "Mean provider control missing: " + key)
    for key in (
        "finite_N_velocity_derivative_errors",
        "numerical_QZZ",
        "full_C5_numeric_containment",
        "signed_mean_function_values_evaluated",
        "global_five_moment_closure",
        "full_goal_complete",
    ):
        require(packet.get(key) is False, "Unproved completion gate changed: " + key)
    require(packet.get("zero_is_symmetric_error_box_center_not_selected_value") is True
            and packet.get("old_signed_sources_and_all_error_caps_unchanged") is True
            and packet.get("accepted_ancestors_replayed") is False,
            "Mean source ancestry controls changed")

    source_pins = packet.get("source_pins")
    require(isinstance(source_pins, dict), "Current mean source pins missing")
    for path, expected in ((PARENT, PARENT_SHA256), (KERNEL, KERNEL_SHA256)):
        row = source_pins.get(path)
        require(isinstance(row, dict) and row.get("path") == path
                and row.get("expected_sha256") == expected,
                "Reviewed mean source pin missing: " + path)

    return dict(
        passed=True,
        actual_Cy_Z3_and_complete_physical_mean_Z3_log_caps=True,
        ordinary_Z0_Z3_two_charts=True,
        joined_history_and_second_memory_retained=True,
        physical_Pstar_E_hb_N_scales_pinned=True,
        finite_N_velocity_QZZ_fullC5_globalclosure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(REVIEWED_MODULE_SHA256 != "REVIEW_PENDING",
            "Independent review hash pending")
    require(sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256,
            "Only the independently reviewed mean module may produce this report")
    require(sha(ROOT / KERNEL) == KERNEL_SHA256,
            "Changed reviewed mean kernel")
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
        require(sha(ROOT / path) == row["expected_sha256"], "Changed new input: " + path)
        require(target.get(name, row["expected_sha256"]) == row["expected_sha256"],
                "Conflicting source pin: " + path)
        target[name] = row["expected_sha256"]

    implementation.update({path: sha(ROOT / path) for path in CODE})
    report = copy.deepcopy(accepted)
    report.update(
        stage="actual_micro_Cy_and_complete_physical_mean_Z3_bounds",
        parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_mean_Z3_bounds=packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
    )
    report["contract"] = dict(
        actual_Cy_Z3_and_complete_physical_mean_Z3_log_caps=True,
        ordinary_Z0_to_Z3_two_charts=True,
        joined_history_and_second_memory_retained=True,
        physical_Pstar_E_hb_N_scales_pinned=True,
        finite_N_velocity_QZZ_fullC5_globalclosure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Actual same-source C/Cy phase primitives and joined moving histories "
            "through ordinary Z3 on both micro charts, inherited second-chart memory, "
            "physical Pstar/E/hb/N scale records, and complete physical mean/total-y "
            "log caps; velocity finite-N errors, QZZ, C5 containment, and global "
            "closure remain open"
        ),
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independently_reviewed_kernel_sha256=sha(ROOT / KERNEL),
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
        == REVIEWED_MODULE_SHA256
        and checks["independently_reviewed_kernel_sha256"]
        == report["implementation_hashes"][KERNEL],
        "Same independently reviewed implementation required",
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
