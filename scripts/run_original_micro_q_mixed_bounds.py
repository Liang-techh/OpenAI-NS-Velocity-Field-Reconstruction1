"""Save or inspect actual q and scaled-Poisson mixed logarithmic bounds.

This runner preserves the accepted b/t0/p2 packet and records the current
same-source q and scaled-Poisson input provider.  It deliberately leaves
phase or inverse derivatives, finite-N errors, C_y Z3, mean Z3, QZZ, full C5
containment, and the global five-moment closure open.
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
from openai_ns_reconstruction import original_micro_q_mixed_bounds as current


build = current.build
PARENT = current.PARENT
PARENT_SHA256 = current.PARENT_SHA256
LOG_HELPER = current.LOG_HELPER
SIGMA_HELPER = current.SIGMA_HELPER
HELPER_HASHES = dict(current.HELPER_HASHES)
require = current.require
sha = current.sha

# The parent replaces this placeholder after its independent review freezes
# the source module.  A report cannot be written while the review is pending.
REVIEWED_MODULE_SHA256 = "73b70df8ed52e05fb45c1a95789ed38b75f70f83bfe9f002737accbcd8b9645d"

NAME = "lei_ren_part1_paper_compliant_current_original_micro_q_mixed_bounds.json"
CODE = (
    "src/openai_ns_reconstruction/original_micro_q_mixed_bounds.py",
    "scripts/run_original_micro_q_mixed_bounds.py",
    LOG_HELPER,
    SIGMA_HELPER,
)


def _load(path):
    return json.loads(Path(path).read_bytes())


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted b/t0/p2 parent changed")
    report = _load(ROOT / PARENT)
    require(report["checks"]["passed"] is True, "Checked b/t0/p2 parent required")
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


def _rows(rows, label):
    """Check the eight (physical y, ordinary Z) rows without value replay."""
    require(isinstance(rows, list) and len(rows) == 8, label + " must have eight rows")
    expected = {(i, k) for i in (0, 1) for k in range(4)}
    actual = set()
    for row in rows:
        require(isinstance(row, dict), label + " row is not an object")
        key = (row.get("physical_y_order"), row.get("ordinary_Z_order"))
        actual.add(key)
        require(isinstance(row.get("exact_zero"), bool), label + " zero flag missing")
        if row["exact_zero"]:
            require(row.get("log_absolute_upper_bound") is None,
                    label + " exact-zero row has a log record")
        else:
            record = row.get("log_absolute_upper_bound")
            require(isinstance(record, dict) and record.get("directed") is True,
                    label + " directed log cap missing")
    require(actual == expected, label + " must cover physical y0..y1 and ordinary Z0..Z3")


def _chart(chart, name):
    _passed(chart, "Actual " + name + " chart missing")
    require(chart.get("ordinary_Z_orders") == [0, 1, 2, 3],
            name + " ordinary Z0..Z3 required")
    require(chart.get("physical_y_orders") == [0, 1],
            name + " physical y0..y1 required")
    for key in (
        "q_mixed_log_caps",
        "u_mixed_log_caps",
        "hinv_mixed_log_caps",
        "r_mixed_log_caps",
        "alpha_mixed_log_caps",
        "K_mixed_log_caps",
        "kappa_mixed_log_caps",
        "gamma_mixed_log_caps",
        "cut_argument_mixed_log_caps",
        "sigma_cut_mixed_log_caps",
        "qhat_mixed_log_caps",
    ):
        _rows(chart.get(key), name + " " + key)

    envelope = chart.get("q_coarse_strict_log_envelope")
    _passed(envelope, name + " q coarse logarithmic envelope missing")
    require(envelope.get("multiplier") == 2048
            and envelope.get("reference") == "selected positive logC",
            name + " q envelope reference changed")

    conditioning = chart.get("scaled_Poisson_conditioning")
    _passed(conditioning, name + " scaled Poisson conditioning missing")
    require(
        conditioning.get("exact_gap_identity") == "1-r^2=1/(1+u^2)>0"
        and conditioning.get("exact_D") == "1-2r*cos(psi)+r^2"
        and conditioning.get("D_lower_identity") ==
        "D>=(1-|r|)^2>=(1-r^2)^2/4"
        and conditioning.get("lambda_lower_identity") ==
        "lambda=K*(1+t^2)>=K>0"
        and conditioning.get("free_angle") == ["0", "2*pi"]
        and conditioning.get("r_not_rounded_to_one") is True
        and conditioning.get("q_u_r_hinv_not_materialized") is True
        and conditioning.get("inverse_angle_derivatives_enclosed") is False,
        name + " positive D/lambda source identities missing",
    )
    for key in (
        "one_minus_r_squared_log_lower",
        "full_free_angle_D_log_lower",
        "full_free_angle_lambda_log_lower",
    ):
        require(isinstance(conditioning.get(key), dict)
                and conditioning[key].get("directed") is True,
                name + " " + key + " directed lower bound missing")
    require(chart.get("caps_not_signed_values") is True
            and chart.get("sqrt_not_evaluated_on_flat_branch") is True,
            name + " cap and flat-branch controls missing")


def _packet_fields(report, accepted):
    require(
        report["parent_report"] == Path(PARENT).name
        and report["parent_report_sha256"] == PARENT_SHA256,
        "Same signed b/t0/p2 predecessor required",
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

    packet = report["actual_micro_q_mixed_bounds"]
    _passed(packet, "Actual q mixed packet missing")
    require(packet["source_family"] == accepted["source_family"], "Source family changed")

    for key, label in (
        ("parameter_envelopes", "Parameter envelope proof missing"),
        ("ordinary_log_mixed_controls", "Ordinary log-mixed controls missing"),
        ("scalar_sigma_derivative_proof", "Scalar sigma derivative proof missing"),
        ("scalar_scaled_Poisson_derivative_proof", "Scalar scaled-Poisson proof missing"),
        ("direct_source_provenance", "Direct source provenance missing"),
        ("natural_collar_flat_proof", "Natural collar/flat proof missing"),
        ("shear_flat_glue", "Shear/flat glue proof missing"),
    ):
        _passed(packet.get(key), label)

    parameters = packet["parameter_envelopes"]
    require(
        parameters.get("approximate_ratios_not_exact_equalities") is True
        and parameters.get("eta_hb_N_not_materialized") is True
        and parameters.get("C_denotes_selected_positive_logC") is True,
        "Selected parameter logarithmic controls missing",
    )
    ratios = parameters.get("log_ratios_to_selected_logC")
    require(isinstance(ratios, dict) and {"eta", "hb", "N"}.issubset(ratios),
            "Selected eta/hb/N log ratios missing")

    mixed = packet["ordinary_log_mixed_controls"]
    require(mixed.get("ordinary_derivatives") is True
            and mixed.get("log_absolute_only") is True
            and mixed.get("zero_sentinel_is_exact") is True
            and mixed.get("no_underlying_exp_or_signed_value") is True,
            "Ordinary log-mixed controls changed")
    scalar = packet["scalar_scaled_Poisson_derivative_proof"]
    require(scalar.get("globally_real_scalar_caps") is True
            and scalar.get("no_sampling_or_near_unit_rounding") is True,
            "Scalar Poisson cap controls changed")

    provenance = packet["direct_source_provenance"]
    require(provenance.get("source_producers_or_interface_checker_replayed") is False
            and provenance.get("alpha_one_core_history_functional_identity_consumed") is True
            and provenance.get("original_sc_is_Z_independent_exact_parameter") is True,
            "Direct source provenance controls changed")

    require(provenance.get("actual_active_Delta_less_eta_gate_consumed") is True
            and provenance.get("actual_a_and_E_positive_gate_consumed") is True
            and provenance.get("original_lazy_flat_checked_rows") == 144
            and provenance.get("active_inverse_and_sqrt_not_evaluated_on_flat") is True,
            "Actual active/flat source admission gates missing")
    collar = packet["natural_collar_flat_proof"]
    require(
        collar.get("natural_flat_neighborhood") == ["0", "1/512"]
        and collar.get("one_sided_at_s_zero") is True
        and collar.get("forced_collar_strictly_inside_natural_neighborhood") is True
        and collar.get("sc_positive_less_than_2_power_minus10") is True
        and collar.get("comparison_alpha_exact_one") is True
        and collar.get("q_exact_zero_on_natural_neighborhood") is True
        and collar.get("physical_y_and_Z_jets_through_total_order_four_zero") is True
        and collar.get("original_hard_collar_branch_not_smoothed_or_replaced") is True,
        "Natural collar and flat controls changed",
    )
    glue = packet["shear_flat_glue"]
    require(
        glue.get("flat_predicate") == "Delta>=eta or original first collar"
        and glue.get("source_q_and_all_mixed_derivatives_exact_zero_on_flat") is True
        and glue.get("gamma_at_shear_boundary_at_least_eta") is True
        and glue.get("sigma_flat_derivatives_through_order_four") is True
        and glue.get("active_formula_uniform_bounds_extend_through_flat_boundary") is True,
        "Shear flat glue controls changed",
    )

    ast_receipt = packet.get("source_AST")
    require(isinstance(ast_receipt, dict), "Source AST receipt missing")
    for key in ("q", "roots", "comparison"):
        _passed(ast_receipt.get(key), "Source AST " + key + " receipt missing")

    charts = packet.get("charts")
    require(isinstance(charts, dict), "Actual q charts missing")
    _chart(charts.get("first"), "first")
    _chart(charts.get("second"), "second")

    require(packet.get("ordinary_Z_orders") == [0, 1, 2, 3]
            and packet.get("physical_y_orders") == [0, 1],
            "Actual q domain orders changed")
    require(packet.get("q_mixed_provider_complete") is True
            and packet.get("actual_scaled_Poisson_inputs_mixed_provider") is True
            and packet.get("actual_free_angle_positive_D_and_lambda_log_bounds") is True,
            "Actual q and scaled-Poisson providers required")

    # These are intentionally false: this packet is the requested partial
    # q provider and does not claim the remaining finite-N or C5 closure.
    for key in (
        "phase_or_inverse_mixed_provider_complete",
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
    require(isinstance(source_pins, dict), "Current source pins missing")
    for path in (LOG_HELPER, SIGMA_HELPER):
        row = source_pins.get(path)
        require(isinstance(row, dict) and row.get("path") == path
                and row.get("expected_sha256") == HELPER_HASHES[path],
                "Reviewed helper source pin missing: " + path)

    return dict(
        passed=True,
        actual_q_and_scaled_Poisson_mixed_log_caps=True,
        ordinary_Z0_Z3_and_physical_y0_y1=True,
        positive_D_and_lambda_conditioning=True,
        natural_collar_and_shear_flat_glue=True,
        phase_or_inverse_finite_N_Cy_Z3_mean_Z3_QZZ_fullC5_globalclosure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )


def run():
    started = time.perf_counter()
    require(REVIEWED_MODULE_SHA256 != "REVIEW_PENDING",
            "Independent review hash pending")
    require(sha(ROOT / CODE[0]) == REVIEWED_MODULE_SHA256,
            "Only the independently reviewed q module may produce this report")
    for path, expected in HELPER_HASHES.items():
        require(sha(ROOT / path) == expected, "Changed reviewed helper: " + path)
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
        stage="actual_micro_q_mixed_source_jet_bounds",
        parent_report=Path(PARENT).name,
        parent_report_sha256=PARENT_SHA256,
        actual_micro_q_mixed_bounds=packet,
        input_hashes=inputs,
        implementation_hashes=implementation,
    )
    report["contract"] = dict(
        actual_q_and_scaled_Poisson_mixed_log_caps=True,
        ordinary_Z0_to_Z3_and_physical_y0_to_y1=True,
        positive_D_and_lambda_conditioning=True,
        natural_collar_and_shear_flat_glue=True,
        phase_or_inverse_finite_N_Cy_Z3_mean_Z3_QZZ_fullC5_globalclosure_open=True,
        old_signed_sources_and_error_caps_retained=True,
        accepted_ancestors_replayed=False,
    )
    report["checks"] = dict(
        passed=True,
        packet_invariants=_packet_fields(report, accepted),
        independent_new_math_review_completed=True,
        independent_review_completed=True,
        independent_review_scope=(
            "Actual same-source q and scaled-Poisson mixed logarithmic caps through "
            "ordinary Z3/physical y1; positive D/lambda conditioning and natural "
            "collar/shear flat glue; phase/inverse, finite-N, Cy_Z3, mean_Z3, QZZ, "
            "C5 containment, and global closure remain open"
        ),
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independently_reviewed_log_helper_sha256=sha(ROOT / LOG_HELPER),
        independently_reviewed_sigma_helper_sha256=sha(ROOT / SIGMA_HELPER),
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
        and checks["independently_reviewed_log_helper_sha256"]
        == report["implementation_hashes"][LOG_HELPER]
        and checks["independently_reviewed_sigma_helper_sha256"]
        == report["implementation_hashes"][SIGMA_HELPER],
        "Same independently reviewed implementation required",
    )
    for path_key, expected in HELPER_HASHES.items():
        require(report["implementation_hashes"][path_key] == expected,
                "Same independently reviewed helper required: " + path_key)
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
