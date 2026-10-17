"""Enclose signed actual lower m15 and join its shared q1 pressure seam."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "experiments/root_st073"
sys.path.insert(0,str(ROOT/"src"))
from openai_ns_reconstruction.original_finite_N_signed_source import digest, load, require, sha
from openai_ns_reconstruction.original_n1_signed_lower_m15 import build, PARENT, PARENT_SHA256
NAME = "lei_ren_part1_paper_compliant_current_original_n1_signed_lower_m15.json"
CODE = ("src/openai_ns_reconstruction/original_n1_signed_lower_m15.py", "scripts/run_original_n1_signed_lower_m15.py")


def parent():
    require(sha(ROOT/PARENT)==PARENT_SHA256,"Accepted reserved-point parent changed")
    p=load(ROOT/PARENT);require(p["checks"]["passed"] is True,"Passed source parent required")
    return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,"input_hashes"),(ROOT,"implementation_hashes")):
        for n,h in r[key].items():
            require(sha(base/n)==h,"Changed source pin: "+n);count+=1
    return count


def merge_pins(p,a):
    inputs,impl=dict(p["input_hashes"]),dict(p["implementation_hashes"])
    for row in a["source_pins"].values():
        n,h=row["path"],row["expected_sha256"]
        require(sha(ROOT/n)==h,"Changed new lower source: "+n)
        if n.startswith("experiments/root_st073/"):
            n=n[len("experiments/root_st073/"):];target=inputs
        else:target=impl
        require(target.get(n,h)==h,"Conflicting source pin: "+n);target[n]=h
    return inputs,impl


def packet_check(r,p):
    require(r["parent_report"]==Path(PARENT).name and r["parent_report_sha256"]==PARENT_SHA256,"Different lower parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes",
             "contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source/error packet changed: "+k)
    a=r["actual_signed_lower_m15"]
    require(a["passed"] is True and a["parent_mapping_sha256"]==digest(p) and a["source_family"]==p["source_family"] and
            a["ordinary_Z_orders"]==[0,1] and a["lower_source_handles"]==[36592,36593] and
            a["inner_source_handles"]==[53161,53162],"Actual lower source/derivative scope changed")
    true_flags=("raw_velocity_cutoff_identically_zero_on_lower","current_native_first_window_covers_q0_to_q2",
        "same_core_inlet_source_used_for_both_lower_endpoints","independent_endpoint_Psi_uncertainty_subtraction_not_used",
        "shared_Qaxis_Psi_and_mean_constants_canceled_before_boxing","actual_signed_lower_value_enclosed",
        "actual_inner_plus_lower_value_enclosed","own_mean_and_pressure_histories_retained",
        "distinct_nonzero_flat_endpoint_atoms_retained","new_join_uses_direct_zero_to_2a_primitive",
        "shared_q1_pressure_endpoint_canceled_as_function","independent_q1_pressure_model_error_box_not_added_to_new_join",
        "old_six_signed_terms_and_all_error_caps_unchanged","P1_is_distinct_from_m15","one_point_only",
        "accepted_C_inner_contains_identical_q1_pressure_primitive","lower_direct_differential_source_bound",
        "lower_sc_and_hb_Z_independent","both_prescribed_flat_atoms_ordinary_Z_independent",
        "reserved_comparison_fixed_units_checked")
    false_flags=("old_finite_prefix_equals_actual_raw_function","old_six_numerical_source_join_closed",
        "forward_Rp_discrepancy_enclosed","all_Z_functional_coverage","global_five_moment_closure","full_goal_complete",
        "accepted_ancestors_replayed","source_executed")
    require(all(a[k] is True for k in true_flags) and all(a[k] is False for k in false_flags),"Lower/prefix scope contract failed")
    require(a["exact_algebra"]["passed"] is True and all(a["exact_algebra"]["checks"].values()),"Pressure seam algebra failed")
    for k in ("0","1"):
        proof=a["accepted_q1_pressure_core_and_velocity_assembly"][k]
        require(proof["passed"] is True and
            proof["q1_pressure_digest"]==digest(p["actual_correlated_leading_pressure_increment"]["ordinary_Z"][k]["m15_leading_N1p_local"]) and
            proof["first_collar_digest"]==digest(p["actual_complete_first_collar_m15_increment"]["ordinary_Z"][k]) and
            proof["actual_core_and_inner_digest"]==digest(p["actual_signed_C_inner_m15"]["ordinary_Z"][k]),
            "Joined q1 seam does not bind the accepted original pressure/velocity/core rows")
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction.original_n1_selected_word_numeric import interval
    from openai_ns_reconstruction.original_n1_velocity_m15_increment import saved_polynomial
    from openai_ns_reconstruction.original_n1_selected_pressure_primitive import add,export
    from openai_ns_reconstruction.original_n1_signed_lower_m15 import extend
    c=MPIntervalContext();c.prec=1800
    for k in ("0","1"):
        row=a["ordinary_Z"][k]
        lower=row["actual_lower_m15"];joined=row["actual_inner_plus_lower_m15"]
        for value in (lower,joined):
            cert=value["normalized"]
            require(cert["strict_sign"]=="negative" and interval(c,cert["normalized_integral_enclosure"]).b<0 and
                    interval(c,cert["relative_interval_width_upper"]).b<c.mpf("1e-27") and
                    cert["caps_not_selected_as_values"] is True,"Actual signed lower/prefix result failed")
        ls=add(*(saved_polynomial(c,row[n]) for n in ("regular_viscous","flat_start","flat_end","axis_inertial")))
        le=add(*(saved_polynomial(c,row[n]) for n in ("paired_pressure_error_radius","lower_inertial_error_radius")))
        require(export(ls)==lower["signed_approximation"] and export(le,error=True)==lower["error_radius"],
                "Serialized lower pressure source/error does not recompose")
        core=p["actual_signed_C_inner_m15"]["ordinary_Z"][k]["actual_core_m15"]
        js=add(extend(saved_polynomial(c,core["signed_approximation"])),
            saved_polynomial(c,row["retained_cut_velocity_signed"]),
            saved_polynomial(c,row["direct_zero_to_2a_leading_pressure"]["signed_approximation"]))
        je=add(extend(saved_polynomial(c,core["error_radius"])),
            saved_polynomial(c,row["retained_cut_velocity_error"]),
            saved_polynomial(c,row["direct_zero_to_2a_leading_pressure"]["error_radius"]))
        require(export(js)==joined["signed_approximation"] and export(je,error=True)==joined["error_radius"],
                "Joined core/cut-velocity/direct pressure source does not recompose")
        point=row["point_reserved_comparison"];old=p["actual_reserved_point_m15"]["ordinary_Z"][k]
        require(point["preRp_common_unit_enclosure"]==old["preRp_common_unit_enclosure"] and
            point["independent_negative_complete_pulse_target"]==old["independent_negative_complete_pulse_target"] and
            point["original_pulse_consumed_once"] is True and point["prefix_cover_is_not_function_value"] is True,
            "Independent target/preRp errors or cover semantics changed")
        diff=-(interval(c,point["joined_prefix_common_unit_enclosure"])+interval(c,point["preRp_common_unit_enclosure"]))
        target=interval(c,point["independent_negative_complete_pulse_target"])
        require(diff._mpi_==interval(c,point["correlated_reserved_minus_target"])._mpi_ and
            (target+diff)._mpi_==interval(c,point["reserved_response"])._mpi_ and
            interval(c,point["relative_correlated_discrepancy"]).b<c.mpf("2e-432"),"Point source comparison failed")
    require(all(v is False for v in r["limits"].values()),"Lower source result promoted to completed physical layers")
    return dict(passed=True,actual_signed_lower_Z0_Z1_enclosed=True,
        shared_q1_pressure_endpoint_canceled_before_boxing_in_new_join=True,
        core_cut_velocity_and_direct_zero_to_2a_pressure_each_once=True,
        inherited_six_signed_terms_and_all_error_caps_unchanged=True,
        all_Z_or_forward_Rp_matching=False,accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p);inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="actual_signed_lower_m15_and_correlated_q1_pressure_seam",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,actual_signed_lower_m15=a,
        input_hashes=inputs,implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(actual_lower_m15_value_and_ordinary_Z1_enclosed=True,
        native_q1_to_q2_covered_by_accepted_first_window=True,
        independent_nonzero_prescribed_shear_flat_atoms_retained=True,
        shared_q1_pressure_endpoint_canceled_by_direct_FTC=True,
        all_old_error_caps_and_six_signed_terms_retained=True,
        all_Z_or_forward_Rp_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),independent_new_math_review_completed=True,
        independent_review_scope="Actual lower leading-N1p m15, same-inlet endpoint correlation and direct q0..2 pressure join at saved Z",
        independent_review_resolution="Lower pressure factors and errors, identical accepted q1 source assembly, direct joined FTC and constant physical units accepted",
        independently_reviewed_module_sha256=sha(ROOT/CODE[0]),independent_runtime_worker_model_metadata_available=False,
        focused_exact_algebra_and_actual_source_arithmetic=True,accepted_ancestors_replayed=False)
    r["execution_seconds"]=round(time.perf_counter()-started,3);return r


def inspect(path):
    r=load(path)
    require(r["checks"]["passed"] is True and r["checks"]["independent_new_math_review_completed"] is True and
        r["checks"]["independently_reviewed_module_sha256"]==r["implementation_hashes"][CODE[0]],
        "Same independently reviewed passed implementation required")
    return dict(passed=True,pins=pins(r),packet_invariants=packet_check(r,parent()),serialized_only=True,
                owner_fixture_recurrence_norm_or_quadrature_replayed=False)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report",action="store_true");group.add_argument("--inspect-report",action="store_true")
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8",newline="\n");tmp.replace(path)
        print(json.dumps(dict(passed=True,execution_seconds=r["execution_seconds"],bytes=path.stat().st_size,sha256=sha(path))))
    else:print(json.dumps(inspect(path),indent=2))
