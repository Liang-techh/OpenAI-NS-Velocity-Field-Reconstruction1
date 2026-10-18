"""Save or inspect actual ordinary-source b/t0/p2 mixed logarithmic bounds.

This runner preserves the accepted E/a packet and records the partial
same-source b/t0/p2 provider.  It deliberately leaves q, physical C_y Z3,
mean Z3, QZZ, and the global five-moment closure open.
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
# source-provenance contract without replaying any accepted producer.
from openai_ns_reconstruction import original_micro_b_p2_mixed_bounds as current


build = current.build
PARENT = current.PARENT
PARENT_SHA256 = current.PARENT_SHA256
Z5 = current.Z5
PRESSURE = current.PRESSURE
require = current.require
sha = current.sha

# The parent replaces this placeholder after its independent review freezes
# the source module.  A report cannot be written while the review is pending.
REVIEWED_MODULE_SHA256 = "2aaebf895d7d3a1b01a209055016ab6350c0c37a74aeda2182eaf9386918aa78"

NAME = "lei_ren_part1_paper_compliant_current_original_micro_b_p2_mixed_bounds.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_b_p2_mixed_bounds.py",
    "scripts/run_original_micro_b_p2_mixed_bounds.py",
    Z5,
    PRESSURE,
)


def _load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted E/a parent changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True, "Checked E/a parent required")
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


def _row_orders(rows, orders, label):
    require(isinstance(rows, list), label + " rows missing")
    require(
        [row.get("ordinary_Z_order") for row in rows] == list(orders),
        label + " orders missing",
    )


def _packet_fields(report, accepted):
    require(
        report["parent_report"] == Path(PARENT).name
        and report["parent_report_sha256"] == PARENT_SHA256,
        "Same signed E/a predecessor required",
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
            require(report[key] == value, "Inherited signed packet or error changed: " + key)

    packet = report["actual_micro_b_p2_mixed_bounds"]
    _passed(packet, "Actual b/t0/p2 packet missing")
    require(packet["source_family"] == accepted["source_family"], "Source family changed")
    require(packet.get("ordinary_derivatives_and_original_inlets") is True,
            "Original ordinary derivatives and inlets required")

    extension = packet.get("source_extension")
    _passed(extension, "Current source extension receipt missing")
    require(
        extension.get("ordinary_Z_orders") == list(range(6))
        and extension.get("direction_ordinary_Z_orders") == list(range(5))
        and extension.get("old_source_Z0_Z4_consumed_without_recomputation") is True
        and extension.get("old_Dbar_Z0_Z3_consumed_byte_exact") is True
        and extension.get("new_Dbar_Z4_and_yZ4_only") is True
        and extension.get("ordinary_derivatives_no_factorial_discount") is True
        and extension.get("no_zero_padding") is True
        and extension.get("no_accepted_writer_or_producer_replayed") is True,
        "Actual source-extension controls missing",
    )

    inlet = packet.get("complete_inlet_integral_bound")
    _passed(inlet, "Complete uniform inlet integral provider missing")
    require(inlet.get("same_complete_inlets_passed_to_native_and_comparison") is True,
            "Both paths must retain the same complete inlet functions")
    require({(v["scaled_radial_order"],v["ordinary_Z_order"])
             for v in extension["new_core_rows"]["new_rows"]} == {(0,5),(1,5)},
            "Only two new core Xh rows are admitted")
    pressure = packet.get("whole_tube_pressure")
    _passed(pressure, "Whole-tube pressure provider missing")
    require(
        pressure.get("source_executed") is False
        and pressure.get("accepted_ancestors_replayed") is False
        and pressure.get("ordinary_Z_orders") == [0, 1, 2, 3, 4, 5]
        and pressure.get("normalized_units") == "P0/Pstar^2"
        and pressure.get("numeric_cap_not_promoted_to_value") is True,
        "Whole-tube pressure provenance controls missing",
    )
    physical_scale = pressure.get("physical_scale")
    require(
        isinstance(physical_scale, dict)
        and physical_scale.get("Pstar_squared_not_applied") is True,
        "Pressure physical Pstar scaling must remain explicit",
    )

    _passed(pressure.get("ordinary_jet_controls"), "Exact pressure ordinary-jet controls missing")
    _passed(packet.get("reduced_source_algebra"), "Actual source algebra controls missing")
    amplitude = packet.get("whole_tube_amplitude")
    _passed(amplitude, "Whole-tube amplitude provider missing")
    require(
        amplitude.get("inverse_F0_not_materialized") is True
        and amplitude.get("F0_squared_ordinary_Z0_Z4_coarse_absolute_cap_one") is True
        and amplitude.get("selected_logC_not_replaced_by_minimum_guard") is True,
        "Whole-tube amplitude controls missing",
    )

    native = packet.get("native_angular_and_axial_provider")
    _passed(native, "Native angular/axial provider missing")
    require(
        native.get("ordinary_Z0_Z4") is True
        and native.get("physical_y0_y1") is True
        and native.get("original_core_U_inlet_and_all_corrections_retained") is True
        and native.get("native_and_comparison_never_identified") is True
        and native.get("physical_dy_hb_ds_once") is True,
        "Native history controls missing",
    )

    caps = packet.get("actual_b_t0_p2_mixed_bounds")
    _passed(caps, "Actual b/t0/p2 mixed provider missing")
    require(
        caps.get("ordinary_Z_orders") == [0, 1, 2, 3]
        and caps.get("physical_y_orders") == [0, 1]
        and caps.get("complete_b_t0_p2_mixed_provider") is True,
        "Actual b/t0/p2 ordinary and physical orders required",
    )
    for key in (
        "b_Z_first",
        "b_Z_second",
        "b_yZ_first",
        "b_yZ_second",
        "t0_Z_both",
        "t0_yZ_both",
        "p2_Z_both",
        "p2_yZ_both",
    ):
        _row_orders(caps.get(key), range(4), "b/t0/p2 " + key)
    require(
        caps.get("log_caps_are_not_signed_values") is True
        and caps.get("F0_and_hb_and_N_not_materialized") is True
        and caps.get("b_exact_phi_cancellation") is True
        and caps.get("t0_exact_chi_cancellation") is True
        and caps.get("p2_uses_original_native_histories") is True
        and caps.get("p2_not_identified_with_b") is True
        and caps.get("p2_physical_sqrt_R_over_two_and_pressure_Pstar_squared_retained") is True
        and caps.get("q_mixed_provider_complete") is False
        and caps.get("numerical_Cy_Z3") is False
        and caps.get("numerical_QZZ") is False,
        "b/t0/p2 logarithmic-cap controls missing",
    )

    _passed(caps.get("coarse_logC_envelopes"), "Coarse actual mixed log envelopes missing")
    # These remain intentionally false: this packet is only the requested
    # partial mixed provider.
    for key in (
        "complete_mixed_root_provider",
        "q_mixed_provider_complete",
        "numerical_Cy_Z3",
        "numerical_mean_Z3",
        "numerical_QZZ",
        "global_five_moment_closure",
        "full_goal_complete",
    ):
        require(packet.get(key) is False, "Unproved completion gate changed: " + key)
    require(packet.get("accepted_ancestors_replayed") is False,
            "Accepted ancestor replayed")
    return dict(
        passed=True,
        actual_b_t0_p2_mixed_bounds=True,
        ordinary_Z0_Z3=True,
        physical_y0_y1=True,
        complete_b_t0_p2_mixed_provider=True,
        native_history_provider_passed=True,
        whole_pressure_provider_passed=True,
        whole_amplitude_provider_passed=True,
        q_mixed_provider_open=True,
        Cy_Z3_open=True,
        mean_Z3_open=True,
        QZZ_open=True,
        global_closure_open=True,
        inherited_signed_sources_and_errors_unchanged=True,
        accepted_ancestors_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(REVIEWED_MODULE_SHA256 != "REVIEW_PENDING",
            "Independent review hash pending")
    require(
        sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256,
        "Only the independently reviewed b/t0/p2 module may produce this report",
    )
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
        require(
            target.get(name, row["expected_sha256"]) == row["expected_sha256"],
            "Conflicting source pin: " + path,
        )
        target[name] = row["expected_sha256"]

    implementation.update({path: sha(ROOT / path) for path in CODE})
    report = copy.deepcopy(accepted)
    report.update(
        stage="actual_micro_b_p2_mixed_source_jet_bounds",
        parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_b_p2_mixed_bounds=packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
    )
    report["contract"] = dict(
        actual_b_t0_p2_mixed_log_caps=True,
        ordinary_Z0_to_Z3_and_physical_y0_to_y1=True,
        native_history_and_whole_pressure_amplitude_providers=True,
        q_Cy_Z3_mean_Z3_QZZ_global_closure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Actual partial b/t0/p2 mixed logarithmic caps using native histories, "
            "whole pressure, and amplitude providers; q, Cy_Z3, mean_Z3, QZZ, "
            "and global closure remain open"
        ),
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independently_reviewed_model_sha256=sha(ROOT / Z5),
        independently_reviewed_helper_sha256=sha(ROOT / Z5),
        independently_reviewed_pressure_sha256=sha(ROOT / PRESSURE),
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
        and checks["independently_reviewed_model_sha256"]
        == report["implementation_hashes"][Z5]
        and checks["independently_reviewed_pressure_sha256"]
        == report["implementation_hashes"][PRESSURE],
        "Same independently reviewed implementation required",
    )
    require(all(report["implementation_hashes"][path] == expected
                for path,expected in current.REVIEWED_HELPERS.items()),
            "Same independently reviewed source helpers required")
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
        print(
            json.dumps(
                dict(
                    passed=True,
                    execution_seconds=report["execution_seconds"],
                    bytes=path.stat().st_size,
                    sha256=sha(path),
                )
            )
        )
    else:
        print(json.dumps(inspect(path), indent=2))
