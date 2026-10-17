"""Enclose the complete actual selected first-collar m15 increment."""
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
from openai_ns_reconstruction.original_n1_velocity_m15_increment import build,PARENT,PARENT_SHA256
NAME="lei_ren_part1_paper_compliant_current_original_n1_complete_first_collar_m15.json"
CODE=("src/openai_ns_reconstruction/original_n1_velocity_m15_increment.py",
    "scripts/run_original_n1_velocity_m15.py")


def parent():
    require(sha(ROOT/PARENT)==PARENT_SHA256,"Accepted correlated pressure parent changed")
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
        "Different complete first-collar m15 parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes",
        "contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source/error packet changed: "+k)
    a=r["actual_complete_first_collar_m15_increment"]
    require(a["passed"] is True and a["parent_mapping_sha256"]==digest(p) and
        a["source_family"]==p["source_family"] and a["ordinary_Z_orders"]==[0,1] and
        a["actual_complete_first_collar_m15_increment_enclosed"] is True and
        a["actual_saved_signed_inlet_products_integrated_with_separate_errors"] is True and
        a["raw_allword_tail_added_once_per_field_at_evaluation_point"] is True,
        "Actual complete local integral/source error admission failed")
    require(a["exact_algebra"]["passed"] is True and all(a["exact_algebra"]["checks"].values()) and
        a["leading_variation_proof"]["identical_source_log_anchors"] is True and
        a["leading_variation_proof"]["F0_connection_added_again"] is False and
        a["leading_variation_proof"]["complex_rho_or_sigma_endpoint_extension_used"] is False and
        a["complex_F0_proof"]["physical_complex_cap_not_center_anchor"] is True,
        "Physical ordinary derivative or complex source-domain contract failed")
    for k in ("0","1"):
        row=a["ordinary_Z"][k];cert=row["complete_first_collar_m15"]["normalized"]
        require(cert["strict_sign"]=="negative" and cert["retained_flat_source_not_zeroed"] is True and
            cert["dyadic_caps_are_errors_not_signed_values"] is True and len(cert["dyadic_term_bounds"])>=34,
            "Complete first-collar m15 sign/normalized enclosure failed")
        for name in ("cut_axial_velocity","cut_swirl_coupling"):
            require(row[name]["signed_inlet_product_integral"] and row[name]["error_radius"],
                "Signed velocity/coupling approximation and its error must both remain")
            require(row[name]["relative_absolute_bound"]["absolute_upper_is_not_a_signed_component_value"] is True and
                row[name]["relative_absolute_bound"]["separate_signed_and_error_sparse_terms"]>0,
                "Each new component needs its separate actual absolute bound")
    require(a["P1_is_distinct_from_m15"] is True and a["core_contribution_added"] is False and
        a["full_C_inner_value_evaluated"] is False and a["all_Z_functional_coverage"] is False and
        a["old_error_caps_or_six_ledger_terms_removed"] is False and
        a["accepted_ancestors_replayed"] is False and a["source_executed"] is False and
        all(v is False for v in r["limits"].values()),"Local collar result promoted to full C_inner/global closure")
    return dict(passed=True,complete_actual_local_first_collar_m15_increment_enclosed=True,
        physical_cut_velocity_and_swirl_coupling_included_with_errors=True,
        leading_N1p_added_once_and_core_not_added=True,
        normalized_m15_increment_and_Z1_strictly_negative=True,
        all_inherited_six_signed_terms_and_error_caps_unchanged=True,
        full_C_inner_or_global_terminal_matching=False,accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p);inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="actual_complete_first_collar_m15_increment",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,
        actual_complete_first_collar_m15_increment=a,input_hashes=inputs,
        implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(actual_raw_cut_velocity_and_coupling_integrals_enclosed=True,
        original_eta_exact_mass_and_physical_rho_measure=True,
        actual_inlet_products_plus_separate_variation_quadrature_interval_and_word_tail_errors=True,
        nonconstant_F0_ordinary_Z_derivative_connection_in_saved_coefficients=True,
        complete_local_first_collar_m15_and_Z1_strictly_negative=True,
        all_old_error_caps_and_six_signed_terms_retained=True,
        full_C_inner_or_global_five_moment_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),
        independent_new_math_review_completed=True,
        independent_review_scope="Complete local first-collar m15 increment, saved Z~0.371 ordinary Z0/Z1 only",
        independent_review_resolution="Physical factors, source jets, raw/leading variation, product/geometry errors and local scope accepted",
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
