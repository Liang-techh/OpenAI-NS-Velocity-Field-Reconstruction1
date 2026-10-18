"""Save or inspect same-source fixed-phase inverse mixed logarithmic bounds.

This runner preserves the accepted q and scaled-Poisson packet and records
the active fixed-phase inverse/loop provider.  It leaves total-y finite-N
errors, C_y Z3, mean Z3, QZZ, full C5 containment, and global closure open.
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
from openai_ns_reconstruction import original_micro_phase_inverse_bounds as current


build = current.build
PARENT = current.PARENT
PARENT_SHA256 = current.PARENT_SHA256
HELPER = current.HELPER
HELPER_SHA256 = current.HELPER_SHA256
LOOP = current.LOOP
require = current.require
sha = current.sha

# The parent replaces this placeholder after independent review freezes the
# source module.  A report cannot be written while that review is pending.
REVIEWED_MODULE_SHA256 = "4123926c522cd212f25f3bd5f754b07975f6203da718dbe5ea79e021920ba4ed"

NAME = "lei_ren_part1_paper_compliant_current_original_micro_phase_inverse_bounds.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_phase_inverse_bounds.py",
    "scripts/run_original_micro_phase_inverse_bounds.py",
    HELPER,
)


def _load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted q parent changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True, "Checked q parent required")
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


def _mixed_rows(rows, label):
    """Check the eight (physical y, ordinary Z) rows without value replay."""
    require(isinstance(rows, list) and len(rows) == 8, label + " must have eight rows")
    expected = {(i, k) for i in (0, 1) for k in range(4)}
    actual = set()
    for row in rows:
        require(isinstance(row, dict), label + " row is not an object")
        actual.add((row.get("physical_y_order"), row.get("ordinary_Z_order")))
        require(isinstance(row.get("exact_zero"), bool), label + " zero flag missing")
        if row["exact_zero"]:
            require(row.get("log_absolute_upper_bound") is None,
                    label + " exact-zero row has a log record")
        else:
            record = row.get("log_absolute_upper_bound")
            require(isinstance(record, dict) and record.get("directed") is True,
                    label + " directed log cap missing")
    require(actual == expected, label + " must cover physical y0..y1 and ordinary Z0..Z3")


def _phase_rows(rows, label):
    """Check all 24 (physical y, ordinary Z, free-angle) rows."""
    require(isinstance(rows, list) and len(rows) == 24, label + " must have 24 rows")
    expected = {
        (i, k, j)
        for total in range(5)
        for i in (0, 1)
        for k in range(4)
        for j in range(5)
        if i + k + j == total
    }
    actual = set()
    for row in rows:
        require(isinstance(row, dict), label + " row is not an object")
        key = (
            row.get("physical_y_order"),
            row.get("ordinary_Z_order"),
            row.get("free_angle_order"),
        )
        actual.add(key)
        require(isinstance(row.get("exact_zero"), bool), label + " zero flag missing")
        if row["exact_zero"]:
            require(row.get("log_absolute_upper_bound") is None,
                    label + " exact-zero row has a log record")
        else:
            record = row.get("log_absolute_upper_bound")
            require(isinstance(record, dict) and record.get("directed") is True,
                    label + " directed log cap missing")
    require(actual == expected,
            label + " must cover all physical y/Z/free-angle indices through total order four")


def _chart(chart, name):
    _passed(chart, "Actual " + name + " phase-inverse chart missing")
    require(chart.get("ordinary_Z_orders") == [0, 1, 2, 3]
            and chart.get("physical_y_orders") == [0, 1],
            name + " mixed orders changed")
    require(chart.get("active_branch_only") is True
            and chart.get("free_angle") == ["0", "2*pi"]
            and chart.get("target_phase") == ["0", "1"],
            name + " active domain controls missing")
    for key in ("positive_D_log_lower", "positive_lambda_log_lower"):
        record = chart.get(key)
        require(isinstance(record, dict) and record.get("directed") is True,
                name + " " + key + " directed lower bound missing")

    free = chart.get("free_angle_jets")
    require(isinstance(free, dict), name + " free-angle jets missing")
    require(set(free) == {"D", "t", "T1", "T2", "Phi", "lambda_value"},
            name + " free-angle field set changed")
    for key, rows in free.items():
        _phase_rows(rows, name + " free-angle " + key)

    mixed = chart.get("mixed_log_caps")
    require(isinstance(mixed, dict), name + " mixed logarithmic caps missing")
    require(set(mixed) == {
        "inverse",
        "T1_at_inverse",
        "direction_at_inverse",
        "lambda_at_inverse",
        "A",
        "B_over_Pstar_over_frozen_E0",
        "A_phi",
        "B_over_Pstar_phi_over_frozen_E0",
    }, name + " mixed field set changed")
    for key, rows in mixed.items():
        _mixed_rows(rows, name + " mixed " + key)

    strict = chart.get("strict_log_envelopes")
    require(
        isinstance(strict, dict)
        and strict.get("passed") is True
        and strict.get("all_eight_loop_inverse_fields") == 9000
        and strict.get("A_and_B_over_Pstar_over_frozen_E0_Z0_Z3") == 5600
        and strict.get("reference") == "selected positive logC",
        name + " strict logarithmic envelope controls changed",
    )

    _passed(chart.get("inverse_recurrence"), name + " inverse recurrence proof missing")
    recurrence = chart["inverse_recurrence"]
    require(
        recurrence.get("ordinary_derivatives") is True
        and recurrence.get("normalized_increment_is_taylor_coefficient") is True
        and recurrence.get("coefficient_product_has_no_binomial_factors") is True
        and recurrence.get("phase_terms_divided_by_parameter_and_angle_factorials") is True
        and recurrence.get("returned_parameter_factorials_restored") is True
        and recurrence.get("leading_beta_zero_angle_one_excluded") is True
        and recurrence.get("leading_denominator_log_lower_bound") is True
        and recurrence.get("inverse_value_bound") == "|psi|<=2*pi"
        and recurrence.get("inverse_value_is_not_selected") is True
        and recurrence.get("parameter_rows") == 8,
        name + " inverse recurrence controls changed",
    )
    require(chart.get("E_product_rule_retained") is True
            and chart.get("E_reference_is_point_frozen_not_function_division") is True
            and chart.get("bound_signs_are_not_actual_source_signs") is True,
            name + " E product-rule controls changed")
    require(chart.get("exact_source_formulas") == {
        "A": "a/2*(phase-psi/(2*pi))",
        "B_over_Pstar": "E/2*(-a*T1hat/(2*pi)-b*phase)",
        "A_phi": "a/2*(1-1/(2*pi*lambda_hat))",
        "B_over_Pstar_phi": "E/2*(-a*t_hat/(2*pi*lambda_hat)-b)",
    }, name + " exact source formulas changed")


def _packet_fields(report, accepted):
    require(
        report["parent_report"] == Path(PARENT).name
        and report["parent_report_sha256"] == PARENT_SHA256,
        "Same signed q predecessor required",
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

    packet = report["actual_micro_phase_inverse_bounds"]
    _passed(packet, "Actual phase-inverse packet missing")
    require(packet["source_family"] == accepted["source_family"], "Source family changed")

    ast = packet.get("source_AST")
    _passed(ast, "Source AST receipt missing")
    require(ast.get("function") == "loop_functions"
            and ast.get("path") == LOOP + ".py"
            and ast.get("source_executed") is False,
            "Current loop source AST receipt changed")
    require(set(ast.get("required_text", ())) == {
        "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
        "den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))",
        "t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))",
        "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
        "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
        "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
        "lam=m.mul(m.axial(K),m.add(m.constant(1),t2))",
        "phase_endpoints=[0,1]",
        "derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction'",
        "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
    }, "Phase source AST guard fragments changed")

    root_ast = packet.get("root_source_AST")
    _passed(root_ast, "Generic-root source AST receipt missing")
    require(
        root_ast.get("function") == "generic_roots"
        and root_ast.get("path") == LOOP + ".py"
        and root_ast.get("source_executed") is False,
        "Generic-root source AST receipt changed",
    )
    require(
        set(root_ast.get("required_text", ())) == {
            "aa,bb=divide(C),divide(B)",
            "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
        },
        "Generic-root AST guard fragments changed",
    )

    endpoint = packet.get("endpoint_normalization")
    _passed(endpoint, "Endpoint normalization receipt missing")
    require(
        endpoint.get("cosine_generating_identity_exact") is True
        and endpoint.get("Fourier_series")
        == "H=(cos(psi)-r)/D=sum(n>=1,r^(n-1)*cos(n*psi))"
        and endpoint.get("convergence")
        == "absolute for |r|<1, including r=0 by continuation; accepted 1-r^2>0"
        and endpoint.get("mean_H") == "0"
        and endpoint.get("mean_H_squared") == "1/(2*(1-r^2))"
        and endpoint.get("alpha") == "2*q*hinv"
        and endpoint.get("gap") == "1-r^2=hinv^2"
        and endpoint.get("mean_t") == "t0"
        and endpoint.get("mean_t_squared") == "t0^2+2*q^2"
        and endpoint.get("T1_end") == "2*pi*t0"
        and endpoint.get("T2_end") == "2*pi*(t0^2+2*q^2)"
        and endpoint.get("Phi_zero") == "0"
        and endpoint.get("Phi_twopi") == "1"
        and endpoint.get("Phi_strictly_increasing") is True
        and endpoint.get("inverse_exists_uniquely_for_every_target_phase") is True
        and endpoint.get("exact_t0_identity") == "t0=-b/a"
        and endpoint.get("A_B_both_phase_endpoints_exact_zero") is True
        and endpoint.get("identities_hold_as_parameter_functions_not_at_sampled_points") is True
        and endpoint.get("normalized_ordinary_derivative_conversion_not_used_in_endpoint_proof") is True,
        "Endpoint normalization identities changed",
    )

    controls = packet.get("new_implicit_jet_controls")
    _passed(controls, "New implicit phase-jet controls missing")
    require(controls.get("ordinary_derivatives") is True
            and controls.get("normalized_taylor_coefficients") is True
            and controls.get("no_source_execution") is True,
            "Implicit phase-jet controls changed")

    charts = packet.get("charts")
    require(isinstance(charts, dict) and set(charts) == {"first", "second"},
            "Both phase-inverse charts required")
    _chart(charts["first"], "first")
    _chart(charts["second"], "second")

    require(packet.get("domain") == accepted["actual_micro_q_mixed_bounds"]["domain"]
            and packet.get("ordinary_Z_orders") == [0, 1, 2, 3]
            and packet.get("physical_y_orders") == [0, 1],
            "Phase-inverse domain orders changed")
    require(packet.get("derivative_scope") ==
            "physical slow y and ordinary Z at fixed target phase; not total y"
            and packet.get("total_y_chain_rule") ==
            "d_total_y = d_slow_y + N*d_phase; phase_Z=0",
            "Fixed-phase derivative scope changed")
    require(packet.get("B_scale_convention") ==
            "source B_over_Pstar=E*M/2; physical B=Pstar*B_over_Pstar",
            "B physical-scale convention changed")

    flat = packet.get("flat_branch")
    _passed(flat, "Flat branch proof missing")
    require(
        flat.get("apply_before_active_formula") is True
        and flat.get("q_A_B_Aphi_Bphi_exact_zero") is True
        and flat.get("inverse") == "psi=2*pi*phase"
        and flat.get("slow_inverse_y_Z_derivatives_exact_zero") is True
        and flat.get("phase_inverse_phi_derivative") == "2*pi"
        and flat.get("raw_flat_a_b_not_bounded_by_active_caps") is True
        and flat.get("natural_collar_receipt_consumed") is True
        and flat.get("shear_flat_receipt_consumed") is True,
        "Flat branch controls changed",
    )

    require(packet.get("phase_or_inverse_mixed_provider_complete") is True
            and packet.get("loop_AB_fixed_phase_mixed_provider_complete") is True,
            "Fixed-phase inverse/loop provider required")
    for key in (
        "finite_N_derivative_errors",
        "numerical_Cy_Z3",
        "numerical_mean_Z3",
        "numerical_QZZ",
        "full_C5_numeric_containment",
        "global_five_moment_closure",
        "full_goal_complete",
    ):
        require(packet.get(key) is False, "Unproved completion gate changed: " + key)
    require(packet.get("inherited_signed_sources_and_all_errors_unchanged") is True
            and packet.get("accepted_ancestors_replayed") is False,
            "Accepted source ancestry controls changed")

    source_pins = packet.get("source_pins")
    require(isinstance(source_pins, dict), "Current phase-inverse source pins missing")
    for path, expected in ((PARENT, PARENT_SHA256), (HELPER, HELPER_SHA256)):
        row = source_pins.get(path)
        require(isinstance(row, dict) and row.get("path") == path
                and row.get("expected_sha256") == expected,
                "Reviewed phase-inverse source pin missing: " + path)
    loop_path = LOOP + ".py"
    row = source_pins.get(loop_path)
    require(isinstance(row, dict) and row.get("path") == loop_path,
            "Accepted C5 loop source pin missing")

    return dict(
        passed=True,
        actual_fixed_phase_inverse_and_loop_mixed_log_caps=True,
        ordinary_Z0_Z3_and_physical_y0_y1=True,
        free_angle_order_0_to_4_mixed_jets=True,
        positive_D_lambda_and_inverse_recurrence=True,
        natural_collar_shear_flat_glue_consumed=True,
        finite_N_Cy_Z3_mean_Z3_QZZ_fullC5_globalclosure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(REVIEWED_MODULE_SHA256 != "REVIEW_PENDING",
            "Independent review hash pending")
    require(sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256,
            "Only the independently reviewed phase-inverse module may produce this report")
    require(sha(ROOT / HELPER) == HELPER_SHA256,
            "Changed reviewed phase-jet helper")
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
        stage="actual_micro_fixed_phase_inverse_loop_mixed_bounds",
        parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_phase_inverse_bounds=packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
    )
    report["contract"] = dict(
        actual_fixed_phase_inverse_and_loop_mixed_log_caps=True,
        ordinary_Z0_to_Z3_and_physical_y0_to_y1=True,
        free_angle_order_0_to_4_mixed_jets=True,
        positive_D_lambda_and_inverse_recurrence=True,
        natural_collar_shear_flat_glue_consumed=True,
        finite_N_Cy_Z3_mean_Z3_QZZ_fullC5_globalclosure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Actual same-source fixed-phase inverse and loop A/B mixed logarithmic "
            "caps through ordinary Z3/physical y1, free-angle jets through total "
            "order four, positive D/lambda conditioning, inverse recurrence, and "
            "flat/collar glue; total-y finite-N, Cy_Z3, mean_Z3, QZZ, C5 containment, "
            "and global closure remain open"
        ),
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independently_reviewed_helper_sha256=sha(ROOT / HELPER),
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
        and checks["independently_reviewed_helper_sha256"]
        == report["implementation_hashes"][HELPER],
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
