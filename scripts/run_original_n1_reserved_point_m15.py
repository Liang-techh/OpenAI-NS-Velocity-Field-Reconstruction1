"""Consume the saved signed inner oracle in a correlated reserved m15 point view."""
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
from openai_ns_reconstruction.original_finite_N_signed_source import digest, load, require, sha
from openai_ns_reconstruction.original_n1_reserved_point_m15 import build, PARENT, PARENT_SHA256, FULL, PULSE, DYADIC_BITS
NAME = "lei_ren_part1_paper_compliant_current_original_n1_reserved_point_m15.json"
CODE = ("src/openai_ns_reconstruction/original_n1_reserved_point_m15.py",
        "scripts/run_original_n1_reserved_point_m15.py")


def parent():
    require(sha(ROOT / PARENT) == PARENT_SHA256, "Accepted actual inner parent changed")
    p = load(ROOT / PARENT)
    require(p["checks"]["passed"] is True, "Passed source parent required")
    return p


def pins(r):
    count = 0
    for base, key in ((DIRECTORY, "input_hashes"), (ROOT, "implementation_hashes")):
        for n, h in r[key].items():
            require(sha(base / n) == h, "Changed source pin: " + n)
            count += 1
    return count


def merge_pins(p, a):
    inputs, impl = dict(p["input_hashes"]), dict(p["implementation_hashes"])
    for row in a["source_pins"].values():
        n, h = row["path"], row["expected_sha256"]
        require(sha(ROOT / n) == h, "Changed new source: " + n)
        if n.startswith("experiments/root_st073/"):
            n = n[len("experiments/root_st073/"):]; target = inputs
        else:
            target = impl
        require(target.get(n, h) == h, "Conflicting source pin: " + n)
        target[n] = h
    return inputs, impl


def packet_check(r, p):
    require(r["parent_report"] == Path(PARENT).name and r["parent_report_sha256"] == PARENT_SHA256,
            "Different actual inner point parent")
    mutable = {"stage", "parent_report", "parent_report_sha256", "input_hashes", "implementation_hashes",
               "contract", "checks", "execution_seconds"}
    for k, v in p.items():
        if k not in mutable:
            require(r[k] == v, "Inherited source/error packet changed: " + k)
    a = r["actual_reserved_point_m15"]
    require(a["passed"] is True and a["parent_mapping_sha256"] == digest(p) and
            a["source_family"] == p["source_family"] and a["ordinary_Z_orders"] == [0, 1],
            "Original source provenance changed")
    true_flags = ("common_point_tube_contained_in_inner_tube", "independent_pulse_target_enclosure_admitted",
        "reserved_response_enclosure_admitted", "pulse_function_canceled_before_boxing", "independent_target_consumed_once",
        "inner_evaluation_replaced_once_in_new_point_view", "lower_preRp_saved_errors_preserved",
        "ordinary_Z1_is_physical_derivative_divided_by_fixed_unit", "F0_connection_already_in_source_coefficients",
        "inner_common_unit_cover_is_not_a_function_value", "old_six_signed_contributions_retained", "one_point_only",
        "P1_is_distinct_from_m15", "reserved_response_full_support_only")
    false_flags = ("old_ledger_or_error_caps_changed", "discrepancy_zero_identity_proved", "independent_terminal_matching_closed",
        "all_Z_functional_coverage", "global_five_moment_closure", "full_goal_complete", "source_executed", "accepted_ancestors_replayed",
        "forward_Rp_discrepancy_enclosed", "forward_six_source_function_identity")
    require(all(a[k] is True for k in true_flags) and all(a[k] is False for k in false_flags),
            "Source comparison scope or consumed-once invariant failed")
    require(a["exact_algebra"]["passed"] is True and all(a["exact_algebra"]["checks"].values()) and
            a["local_inner_unit"] == "epsilon^2*Pstar^2" and a["shared_unit"]["physical"] == "reference*A*Rv^2" and
            a["shared_unit"]["U_Z_exact_zero"] is True, "Common physical unit or source algebra failed")
    require(a["exact_core_pressure_scaling_checked"] == p["core_pressure_source_identity"]["scaling"] and
            a["local_dominant_source_exponents_checked"] == [0, 0, 1, -2, 0] and
            a["source_admission"]["cached_joint_full_source_hash_matches"] is True and
            a["source_admission"]["positive_reference_and_energy_epsilon_applied_once"] is True,
            "Exact core epsilon/pressure scaling or same full-source provenance missing")
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction.original_n1_selected_word_numeric import interval
    c = MPIntervalContext(); c.prec = 1800
    full_point = load(ROOT / FULL)["samples"]["point"]
    pulse_point = load(ROOT / PULSE)["samples"]["point"]
    require(a["selected_Z_interval"] == full_point["Z"] == pulse_point["Z"],
            "Serialized comparison must use the same admitted point tube")
    ratio = interval(c, a["local_to_common_unit_log_ratio"])
    for k in ("0", "1"):
        row = a["ordinary_Z"][k]
        part = full_point["pressure_component_enclosures"]["m15"]
        order = int(k)
        require(row["lower_common_unit_enclosure"] == part["lower"]["C1_enclosures"][order] and
                row["preRp_common_unit_enclosure"] == part["preRp"][order] and
                row["complete_pulse_common_unit_enclosure"] == pulse_point["normalized_true_pulse_C1_enclosures"]["m15"][order],
                "A saved actual source/error enclosure was changed")
        local = interval(c, row["actual_inner_local_enclosure"])
        expected_local = p["actual_signed_C_inner_m15"]["ordinary_Z"][k]["actual_C_inner_m15"]["normalized"]["normalized_integral_enclosure"]
        require(local._mpi_ == interval(c, expected_local)._mpi_, "The actual signed inner oracle was changed")
        log_abs = ratio + c.ln(abs(local))
        require(log_abs._mpi_ == interval(c, row["inner_common_unit_log_absolute_enclosure"])._mpi_ and
                log_abs.b < -DYADIC_BITS * c.ln(2), "Serialized physical scale cover failed")
        inner = interval(c, row["inner_common_unit_enclosure"])
        require(inner._mpi_ == c.mpf([-(c.mpf(2)**-DYADIC_BITS), 0])._mpi_, "Inner cover is not the admitted one-sided box")
        target = -interval(c, row["complete_pulse_common_unit_enclosure"])
        diff = -(inner + interval(c, row["lower_common_unit_enclosure"]) + interval(c, row["preRp_common_unit_enclosure"]))
        require(target._mpi_ == interval(c, row["independent_negative_complete_pulse_target"])._mpi_ and
                diff._mpi_ == interval(c, row["correlated_reserved_minus_target"])._mpi_ and
                (target + diff)._mpi_ == interval(c, row["reserved_full_support_response"])._mpi_ and
                (abs(diff) / target)._mpi_ == interval(c, row["relative_correlated_discrepancy_enclosure"])._mpi_,
                "Serialized target, response and correlated discrepancy do not recompose")
        require(interval(c, row["actual_inner_local_enclosure"]).b < 0 and row["inner_strict_sign"] == "negative" and
                interval(c, row["independent_negative_complete_pulse_target"]).a > 0 and
                interval(c, row["reserved_full_support_response"]).a > 0 and
                interval(c, row["relative_correlated_discrepancy_enclosure"]).b < c.mpf("3e-431") and
                row["cap_or_midpoint_selected_as_value"] is False and row["interval_overlap_used_as_closure"] is False,
                "Independent signed target/response comparison failed")
    require(all(v is False for v in r["limits"].values()), "Moment point comparison promoted to full physical closure")
    return dict(passed=True, actual_inner_oracle_consumed_once_in_existing_joint_source=True,
        independent_negative_complete_pulse_m15_and_reserved_response_enclosed=True,
        correlated_difference_evaluated_before_pulse_boxes=True, ordinary_Z0_Z1_and_fixed_physical_units=True,
        all_inherited_six_signed_terms_and_error_caps_unchanged=True,
        all_Z_or_independent_terminal_matching=False, accepted_ancestors_replayed=False)


def run():
    started = time.perf_counter()
    p = parent(); pins(p)
    a = build(ROOT, p); inputs, impl = merge_pins(p, a)
    r = copy.deepcopy(p)
    r.update(stage="actual_reserved_minus_negative_complete_pulse_m15_at_saved_Z",
        parent_report=Path(PARENT).name, parent_report_sha256=PARENT_SHA256,
        actual_reserved_point_m15=a, input_hashes=inputs,
        implementation_hashes=dict(impl, **{n: sha(ROOT / n) for n in CODE}))
    r["contract"] = dict(actual_inner_point_oracle_replaces_only_existing_inner_source=True,
        independent_pulse_target_and_reserved_response_enclosed=True,
        same_pulse_canceled_as_function_before_boxing=True, target_consumed_once=True,
        lower_preRp_source_errors_retained=True, all_old_error_caps_and_six_signed_terms_retained=True,
        ordinary_Z1_and_physical_m15_units_applied_once=True,
        all_Z_or_global_five_moment_closure=False, accepted_ancestors_replayed=False)
    r["checks"] = dict(passed=True, packet_invariants=packet_check(r, p),
        independent_new_math_review_completed=True,
        independent_review_scope="Saved-point inner scale conversion and reserved/independent complete-pulse m15 C1 source comparison",
        independent_review_resolution="Point source identity, exact epsilon/Pstar scale gates, pulse admission and cached source provenance accepted",
        independently_reviewed_module_sha256=sha(ROOT / CODE[0]),
        independent_runtime_worker_model_metadata_available=False,
        focused_exact_algebra_and_actual_source_arithmetic=True, accepted_ancestors_replayed=False)
    r["execution_seconds"] = round(time.perf_counter() - started, 3)
    return r


def inspect(path):
    r = load(path)
    require(r["checks"]["passed"] is True and r["checks"]["independent_new_math_review_completed"] is True and
            r["checks"]["independently_reviewed_module_sha256"] == r["implementation_hashes"][CODE[0]],
            "Same independently reviewed passed implementation required")
    return dict(passed=True, pins=pins(r), packet_invariants=packet_check(r, parent()),
        serialized_only=True, owner_fixture_recurrence_norm_or_quadrature_replayed=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report", action="store_true")
    group.add_argument("--inspect-report", action="store_true")
    args = ap.parse_args(); path = DIRECTORY / NAME
    if args.write_report:
        r = run(); tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(r, indent=2, ensure_ascii=False) + "\n", encoding="utf-8", newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True, execution_seconds=r["execution_seconds"], bytes=path.stat().st_size, sha256=sha(path))))
    else:
        print(json.dumps(inspect(path), indent=2))
