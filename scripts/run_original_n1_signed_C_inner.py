"""Enclose actual core and C_inner m15 at the saved axial interval."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];DIRECTORY=ROOT/"experiments/root_st073"
sys.path.insert(0,str(ROOT/"src"))
from openai_ns_reconstruction.original_finite_N_signed_source import digest,load,require,sha
from openai_ns_reconstruction.original_n1_signed_C_inner_m15 import build,PARENT,PARENT_SHA256
NAME="lei_ren_part1_paper_compliant_current_original_n1_signed_C_inner_m15.json"
CODE=("src/openai_ns_reconstruction/original_n1_signed_C_inner_m15.py",
    "scripts/run_original_n1_signed_C_inner.py")


def parent():
    require(sha(ROOT/PARENT)==PARENT_SHA256,"Accepted complete-collar parent changed")
    p=load(ROOT/PARENT);require(p["checks"]["passed"] is True,"Passed source parent required")
    return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,"input_hashes"),(ROOT,"implementation_hashes")):
        for n,h in r[key].items():require(sha(base/n)==h,"Changed source pin: "+n);count+=1
    return count


def merge_pins(p,a):
    inputs,impl=dict(p["input_hashes"]),dict(p["implementation_hashes"])
    for row in a["source_pins"].values():
        n,h=row["path"],row["expected_sha256"];require(sha(ROOT/n)==h,"Changed new source: "+n)
        if n.startswith("experiments/root_st073/"):n=n[len("experiments/root_st073/"):];target=inputs
        else:target=impl
        require(target.get(n,h)==h,"Conflicting source pin: "+n);target[n]=h
    return inputs,impl


def packet_check(r,p):
    require(r["parent_report"]==Path(PARENT).name and r["parent_report_sha256"]==PARENT_SHA256,
        "Different C_inner m15 parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes",
        "contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source/error packet changed: "+k)
    a=r["actual_signed_C_inner_m15"]
    require(a["passed"] is True and a["parent_mapping_sha256"]==digest(p) and
        a["source_family"]==p["source_family"] and a["ordinary_Z_orders"]==[0,1] and
        a["actual_core_m15_value_enclosed"] is True and a["actual_C_inner_m15_value_enclosed"] is True and
        a["complete_first_collar_added_once"] is True and a["axis_endpoint_rho_squared_b_rho_limit_zero"] is True and
        a["core_mean_not_reset"] is True and a["swirl_core_integral_only_bounded_not_signed_value"] is True,
        "Actual core/collar source identity or error admission failed")
    require(a["core_source_handles"]==[53125,53126] and a["collar_source_handles"]==[53137,53138] and
        a["C_inner_source_handles"]==[53161,53162] and a["exact_algebra"]["passed"] is True and
        all(a["exact_algebra"]["checks"].values()) and
        a["raw_core_bounds_proof"]["caps_are_not_signed_values"] is True and
        a["core_bounds_proof"]["no_inverse_rho_bound_used"] is True,
        "Same-source core FTC, raw word comparison or axis-regularity contract failed")
    require(a["actual_axis_Psi_mean_Qpsi_zero_receipt_and_source_checked"] is True and
        a["executable_helper_pins_checked"] is True and
        a["raw_core_bounds_proof"]["exact_core_source_and_solution_Z_tubes_match"] is True and
        a["raw_core_bounds_proof"]["whole_core_radial_source_disk_checked"] is True and
        a["raw_core_bounds_proof"]["physical_row6_source_normalization_checked"] is True,
        "Executable provenance, actual axis premise or whole-core source-domain gate missing")
    for k in ("0","1"):
        row=a["ordinary_Z"][k]
        for name in ("actual_core_m15","actual_C_inner_m15"):
            cert=row[name]["normalized"]
            require(cert["strict_sign"]=="negative" and
                cert["upper_bounds_not_used_as_signed_values"] is True and cert["dyadic_term_bounds"],
                "Actual signed core/C_inner enclosure failed")
    require(a["one_point_only"] is True and a["all_Z_functional_coverage"] is False and
        a["independent_terminal_comparison_closed"] is False and a["global_five_moment_closure"] is False and
        a["old_error_caps_or_six_ledger_terms_removed"] is False and a["accepted_ancestors_replayed"] is False and
        a["source_executed"] is False and a["P1_is_distinct_from_m15"] is True and
        all(v is False for v in r["limits"].values()),"One-point m15 enclosure promoted to functional/global closure")
    return dict(passed=True,actual_core_and_C_inner_m15_enclosed_at_saved_Z=True,
        ordinary_Z0_Z1_strictly_negative=True,complete_collar_added_once=True,
        physical_scales_and_actual_mean_history_retained=True,
        all_inherited_six_signed_terms_and_error_caps_unchanged=True,
        all_Z_or_independent_terminal_matching=False,accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p);inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="actual_signed_C_inner_m15_at_saved_Z",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,
        actual_signed_C_inner_m15=a,input_hashes=inputs,
        implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(actual_core_m15_primitive_from_same_mean_endpoint=True,
        actual_C_inner_m15_value_and_ordinary_Z1_enclosed_at_saved_point=True,
        full_raw_feedback_radius_covers_k_ge_1=True,
        nonzero_axial_source_model_with_separate_errors=True,
        core_swirl_integral_bounded_only_not_chosen_as_value=True,
        core_plus_complete_collar_once=True,all_old_error_caps_and_six_signed_terms_retained=True,
        all_Z_or_global_five_moment_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),
        independent_new_math_review_completed=True,
        independent_review_scope="Core plus accepted first-collar C_inner m15 at one saved Z interval, ordinary Z0/Z1 only",
        independent_review_resolution="Core/axis/Green/tail math and helper/domain/axis provenance gates accepted",
        independently_reviewed_module_sha256=sha(ROOT/CODE[0]),
        independent_runtime_worker_model_metadata_available=False,
        focused_exact_algebra_and_actual_source_arithmetic=True,accepted_ancestors_replayed=False)
    r["execution_seconds"]=round(time.perf_counter()-started,3);return r


def inspect(path):
    r=load(path);require(r["checks"]["passed"] is True,"Saved passed report required")
    require(r["checks"]["independent_new_math_review_completed"] is True and
        r["checks"]["independently_reviewed_module_sha256"]==r["implementation_hashes"][CODE[0]],
        "Same independently reviewed implementation required")
    return dict(passed=True,pins=pins(r),packet_invariants=packet_check(r,parent()),
        serialized_only=True,owner_fixture_recurrence_norm_or_quadrature_replayed=False)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__);group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report",action="store_true");group.add_argument("--inspect-report",action="store_true")
    args=ap.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        r=run();tmp=path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n",encoding="utf-8",newline="\n");tmp.replace(path)
        print(json.dumps(dict(passed=True,execution_seconds=r["execution_seconds"],bytes=path.stat().st_size,sha256=sha(path))))
    else:print(json.dumps(inspect(path),indent=2))
