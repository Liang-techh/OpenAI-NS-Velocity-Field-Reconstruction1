"""Save or inspect same-source physical U/Q finite-N increment caps.

This runner preserves the accepted complete mean packet and records the
source-linear Q recovery, total-y velocity rows, and R-weighted Q_y rows.  It
does not select signed fields and leaves full QZZ/C5/global closure open.
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
from openai_ns_reconstruction import original_micro_Q_increment_bounds as current


build = current.build
PARENT = current.PARENT
PARENT_SHA256 = current.PARENT_SHA256
RECOVERY = current.RECOVERY
Q4 = current.Q4
LOOP = current.LOOP
CORE_CHECK = current.CORE_CHECK
require = current.require
sha = current.sha

# The parent replaces this placeholder after independent review freezes the
# source module.  A report cannot be written while that review is pending.
REVIEWED_MODULE_SHA256 = "db80b104dd0dc6060e6f10cc9d65a5aa75b447eab11fae38177119affc1da0fb"

NAME = "lei_ren_part1_paper_compliant_current_original_micro_Q_increment_bounds.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py",
    "scripts/run_original_micro_Q_increment_bounds.py",
)


def _load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted mean-Z3 parent changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True, "Checked mean-Z3 parent required")
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
    require(set(required) == set(snippets), label + " source AST fragments changed")
    require(isinstance(receipt.get("AST_sha256"), str),
            label + " source AST digest missing")


def _source_contract(packet):
    contract = packet.get("source_contract")
    _passed(contract, "Q increment source contract missing")
    for key in (
        "accepted_Q_recovery_controls_consumed",
        "actual_physical_U_and_mean",
        "same_source_Q_replacement_not_extra_ledger_contribution",
        "phase_Z_independent_and_total_y_phase_derivative_N",
        "original_mean_transport_rate_one_AST_bound",
        "old_source_function_and_all_error_packets_retained",
        "accepted_recovery_controls_or_source_producers_replayed",
    ):
        expected = False if key.endswith("replayed") else True
        require(contract.get(key) is expected,
                "Q increment source contract changed: " + key)
    require(contract.get("mean_replacement_identity")
            == "same delta_m_y+delta_m=Pstar*B0/N, same zero left inlet and nonreset seam; linear ODE uniqueness"
            and contract.get("full_recovered_Q_equals_original_by_same_U_and_mean_function_identity") is True,
            "Q4 replacement identity controls changed")
    handles = contract.get("source_graph_handles")
    require(isinstance(handles, dict)
            and set(handles) == {"first_micro", "second_micro"},
            "Q recovery source graph chart handles changed")
    for chart, values in handles.items():
        require(isinstance(values, dict)
                and set(values) == {
                    "physical_mean_ordinary_Z0_Z3",
                    "physical_Q_ordinary_Z0_Z2",
                    "original_Q_source_nodes",
                }, "Q recovery source handle fields changed: " + chart)
        require(len(values["physical_mean_ordinary_Z0_Z3"]) == 4
                and len(values["physical_Q_ordinary_Z0_Z2"]) == 3
                and len(values["original_Q_source_nodes"]) == 3,
                "Q recovery source handle orders changed: " + chart)
        for field, rows in values.items():
            require(all(isinstance(value, int) and not isinstance(value, bool)
                        and value >= 0 for value in rows),
                    "Q recovery source handles are not graph IDs: " + chart + "." + field)


def _linear_controls(packet):
    controls = packet.get("new_linear_Q_controls")
    _passed(controls, "New linear Q controls missing")
    require(
        controls.get("operator_linearity_exact") is True
        and controls.get("nonunit_delta_control") == "2/7"
        and controls.get("ordinary_quotient_rows_checked") == [0, 1, 2]
        and controls.get("ordinary_Z_factorials_not_inserted") is True
        and controls.get("Q_error_recovered_before_absolute_norm") is True
        and controls.get("mean_y_identity") == "delta_m_y=delta_U-delta_m"
        and controls.get("total_y_velocity")
        == "delta_U_y=Pstar*(slow_B_y/N+B_over_Pstar_phi)"
        and controls.get("no_extra_N_on_B_phi") is True
        and controls.get("physical_Pstar_applied_once") is True
        and controls.get("source_linear_error_not_subtraction_of_large_interval_boxes") is True,
        "Linear Q controls changed",
    )


def _recovery(proof, count, label):
    _passed(proof, label + " proof missing")
    for key in ("positive_L", "L_Z", "L_ZZ"):
        _directed(proof.get(key), label + " " + key)
    _rows(proof.get("numerator_log_caps"), count, label + " numerator")
    require(proof.get("ordinary_derivatives") is True
            and proof.get("no_top_row_padding") is True,
            label + " quotient controls changed")


def _chart(chart, name):
    _passed(chart, "Actual " + name + " U/Q chart missing")
    require(chart.get("source_chart") == name + "_micro",
            name + " source chart changed")
    fields = {
        "physical_delta_U": 4,
        "physical_delta_m": 4,
        "physical_delta_m_y": 4,
        "physical_delta_Q": 3,
        "physical_delta_U_y": 2,
        "physical_delta_Q_y": 2,
        "physical_R_delta_Q_y": 2,
    }
    caps = chart.get("log_caps")
    require(isinstance(caps, dict) and set(caps) == set(fields),
            name + " seven U/Q log-cap fields changed")
    for key, count in fields.items():
        _rows(caps[key], count, name + " " + key)
    ratios = chart.get("log_cap_ratios")
    require(isinstance(ratios, dict) and set(ratios) == set(fields),
            name + " seven U/Q ratio fields changed")
    for key, count in fields.items():
        _ratio_rows(ratios[key], count, name + " " + key)

    _recovery(chart.get("Q_recovery"), 3, name + " Q")
    _recovery(chart.get("Q_y_recovery"), 2, name + " Q_y")
    require(chart.get("physical_delta_U_Z0_Z3_below_2_power_minus1024") is True
            and chart.get("physical_delta_Q_Z0_Z2_below_2_power_minus1024") is True,
            name + " local U/Q thresholds missing")
    strict = chart.get("strict_log_envelopes")
    _passed(strict, name + " strict U/Q envelope missing")
    require(strict.get("U_multiplier") == -600
            and strict.get("Q_multiplier") == -2600
            and strict.get("reference") == "selected positive logC",
            name + " strict U/Q envelope changed")
    outer = chart.get("outer_RQy_smallness")
    require(isinstance(outer, list) and len(outer) == 2,
            name + " outer RQy rows missing")
    _orders = []
    for row in outer:
        require(isinstance(row, dict)
                and isinstance(row.get("below_2_power_minus1024"), bool),
                name + " outer RQy boolean control changed")
        _orders.append(row.get("ordinary_Z_order"))
    require(_orders == [0, 1], name + " outer RQy ordinary-Z order changed")
    require(chart.get("non_small_bounds_retained") is True
            and chart.get("total_y_fast_Bphi_term_not_divided_by_N") is True
            and chart.get("physical_R_applied_once") is True,
            name + " R/Qy scale controls changed")


def _packet_fields(report, accepted):
    require(
        report["parent_report"] == Path(PARENT).name
        and report["parent_report_sha256"] == PARENT_SHA256,
        "Same signed mean-Z3 predecessor required",
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

    packet = report.get("actual_micro_Q_increment_bounds")
    _passed(packet, "Actual U/Q increment packet missing")
    require(packet.get("source_family") == accepted["source_family"],
            "Source family changed")
    _source_contract(packet)
    _linear_controls(packet)

    ast = packet.get("source_AST")
    require(isinstance(ast, dict)
            and set(ast) == {
                "recovery", "physical", "total_y", "density", "history",
                "original_fields", "recovered_fields",
            },
            "Q increment source AST field set changed")
    _ast(ast["recovery"], "Q recovery", RECOVERY, "recover_Q", (
        "n0=sub(g.mul(two,z,U[0]),g.add(g.mul(h,z,m[0]),g.mul(d,m[1])))",
        "n2=g.add(g.mul(g.constant(4),U[1]),g.mul(two,z,U[2]),g.mul(two,de,m[1]),",
        "q2=g.div(sub(n2,g.add(g.mul(two,L1,q1),g.mul(L2,q0))),L,cert)",
    ))
    _ast(ast["physical"], "Physical pressure", Q4, "physical_pressure_fields", (
        "U=a.scale(fields['V'],Pstar)",
        "my=a.sub(U,m)",
        "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))",
    ))
    _ast(ast["total_y"], "Total-y velocity", Q4, "total_radial_functions", (
        "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",
    ))
    _ast(ast["density"], "Finite-N density", LOOP, "finite_N_functions", (
        "s=a.scale(loop['B_over_Pstar'],invN)",
        "inc=dict(m=s,",
        "V_N=a.add(V,s)",
    ))
    _ast(ast["history"], "Transport history", LOOP, "transport", (
        "for name,rate in RATES.items():",
        "predecessor_memory_never_reset=True",
        "out[name]=a.add(a.scale(incoming[name],decay),contribution)",
    ))
    _ast(ast["original_fields"], "Original finite-N fields", LOOP, "__init__", (
        "incoming={key:a.fixed(0) for key in RATES}",
        "history,parts,memory=transport(g,finite['five_signed_density_increments'],incoming,",
        "for key in RATES:fields[key]=a.add(a.trim(chart['fields'][key]),history[key])",
        "incoming=self.endpoints[name]",
    ))
    _ast(ast["recovered_fields"], "Recovered physical fields", RECOVERY, "build", (
        'mean=[g.add(g.mul(33572,bg[k]),kernel["complete_physical_delta_m"][name][k]) for k in range(4)]',
        "q=recover_Q(g,U,mean)",
    ))

    require(packet.get("domain") == accepted["actual_micro_mean_Z3_bounds"]["domain"],
            "U/Q increment domain changed")
    charts = packet.get("charts")
    require(isinstance(charts, dict) and set(charts) == {"first", "second"},
            "Both U/Q increment charts required")
    _chart(charts["first"], "first")
    _chart(charts["second"], "second")

    for key in ("actual_delta", "positive_L", "physical_R_upper_log"):
        _directed(packet.get(key), "Actual " + key)
    radius = packet.get("radius_provenance")
    _passed(radius, "Actual radius provenance missing")
    require(
        radius.get("R_equals_epsilon_rho") is True
        and radius.get("epsilon_Lambda_equals_one") is True
        and radius.get("exact_logLambda") == "4*logPstar+1000"
        and radius.get("actual_micro_radius_handle") == 68644
        and radius.get("Q4_RQy_C1_handles") == {
            "first_micro": [129559, 129560],
            "second_micro": [133660, 133661],
        }
        and radius.get("physical_R_applied_once") is True
        and radius.get("physical_R_Z_independent") is True
        and radius.get("same_loop_Q4_radius_graph_prefix") is True
        and radius.get("actual_rho_source_AST_consumed") is True
        and radius.get("core_recipe_receipt_consumed_not_replayed") is True
        and radius.get("source_producers_executed") is False,
        "Actual radius provenance controls changed",
    )
    require(packet.get("actual_physical_delta_U_Z3_error_provider") is True
            and packet.get("actual_physical_delta_Q_Z2_error_provider") is True
            and packet.get("actual_physical_delta_Q_y_Z1_error_provider") is True
            and packet.get("actual_physical_R_delta_Q_y_Z1_error_provider") is True,
            "U/Q increment providers missing")
    for key in (
        "signed_function_values_evaluated",
        "all_outer_RQy_errors_certified_small",
        "six_I1_defect_bounds",
        "complete_full_Q_Z2_numeric_values",
        "full_C5_numeric_containment",
        "global_five_moment_closure",
        "full_goal_complete",
    ):
        require(packet.get(key) is False, "Unproved U/Q completion gate changed: " + key)
    require(packet.get("zero_is_error_box_center_not_selected_value") is True
            and packet.get("old_signed_sources_and_all_error_caps_unchanged") is True
            and packet.get("accepted_ancestors_replayed") is False,
            "U/Q source ancestry controls changed")

    source_pins = packet.get("source_pins")
    require(isinstance(source_pins, dict), "Current U/Q source pins missing")
    for path, expected in ((PARENT, PARENT_SHA256),
                           (RECOVERY, accepted["implementation_hashes"][RECOVERY]),
                           (Q4, accepted["input_hashes"][Path(Q4).name]),
                           (LOOP, accepted["input_hashes"][Path(LOOP).name]),
                           (CORE_CHECK, accepted["input_hashes"][Path(CORE_CHECK).name]),
                           (LOOP[:-3] + ".json.gz",
                            accepted["input_hashes"][Path(LOOP[:-3] + ".json.gz").name]),
                           (Q4[:-3] + ".json.gz",
                            accepted["input_hashes"][Path(Q4[:-3] + ".json.gz").name]),
                           (Q4[:-3] + "_check.json",
                            accepted["input_hashes"][Path(Q4[:-3] + "_check.json").name])):
        row = source_pins.get(path)
        require(isinstance(row, dict) and row.get("path") == path
                and row.get("expected_sha256") == expected,
                "Reviewed U/Q source pin missing: " + path)

    return dict(
        passed=True,
        actual_physical_U_Q_increment_log_caps=True,
        ordinary_U_m_Z0_Z3_Q_Z0_Z2_Qy_Z0_Z1=True,
        positive_L_R_and_Q_recovery=True,
        outer_RQy_smallness_gate_retained_false=True,
        signed_values_full_Q_C5_globalclosure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(REVIEWED_MODULE_SHA256 != "REVIEW_PENDING",
            "Independent review hash pending")
    require(sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256,
            "Only the independently reviewed U/Q module may produce this report")
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
        stage="actual_micro_physical_U_Q_increment_bounds",
        parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_Q_increment_bounds=packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
    )
    report["contract"] = dict(
        actual_physical_U_Q_increment_log_caps=True,
        ordinary_U_m_Z0_Z3_Q_Z0_Z2_Qy_Z0_Z1=True,
        positive_L_R_and_Q_recovery=True,
        outer_RQy_smallness_gate_retained_false=True,
        signed_values_full_Q_C5_globalclosure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Actual same-source physical U/m/m_y caps through ordinary Z3, Q and "
            "Q_y recovery through their available source rows, R-weighted Q_y Z0/Z1, "
            "positive delta/L/R controls, and total-y B_phi scaling; signed values, "
            "outer RQ_y certification, full QZZ, C5 containment, and global closure "
            "remain open"
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
