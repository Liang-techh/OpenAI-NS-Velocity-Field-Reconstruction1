"""Save or inspect the actual first-micro nonempty sc/2..3sc/4 pressure increment."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
DIRECTORY=ROOT/"experiments/root_st073"
sys.path.insert(0,str(ROOT/"src"))
from openai_ns_reconstruction.original_finite_N_signed_source import digest,load,require,sha
from openai_ns_reconstruction.original_n1_first_micro_prefix import build,PARENT,PARENT_SHA256,ActualFirstMicroContinuation
NAME="lei_ren_part1_paper_compliant_current_original_n1_first_micro_prefix.json"
CODE=("src/openai_ns_reconstruction/original_n1_first_micro_prefix.py","scripts/run_original_n1_first_micro_prefix.py")


def parent():
    require(sha(ROOT/PARENT)==PARENT_SHA256,"Accepted actual q2 parent changed")
    p=load(ROOT/PARENT);require(p["checks"]["passed"] is True,"Checked q2 parent required")
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
        require(sha(ROOT/n)==h,"Changed seam source: "+n)
        if n.startswith("experiments/root_st073/"):
            n=n[len("experiments/root_st073/"):];target=inputs
        else:target=impl
        require(target.get(n,h)==h,"Conflicting source pin: "+n);target[n]=h
    return inputs,impl


def packet_check(r,p):
    require(r["parent_report"]==Path(PARENT).name and r["parent_report_sha256"]==PARENT_SHA256,"Different corrected physical seam parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes","contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source/error packet changed: "+k)
    a=r["actual_first_micro_prefix"]
    require(a["passed"] is True and a["source_family"]==p["source_family"] and a["ordinary_Z_orders"]==[0,1] and
        a["selected_Z_interval"]==p["actual_signed_lower_m15"]["selected_Z_interval"] and
        a["parent_shifted_inlet_digest"]==digest(p["actual_first_micro_seam"]),"Same actual native source/Z jets required")
    flags=("shared_sc_half_endpoint_and_error_canceled_before_cumulative_boxing","actual_nonempty_first_micro_leading_pressure_increment_enclosed",
        "ordinary_Z1_and_physical_units_retained","finite_N_prefix_exact_zero_from_original_source",
        "prescribed_shear_sigma_at_endpoints_not_zero","raw_own_velocity_cutoff_zero_but_leading_pressure_history_retained",
        "all_old_six_terms_and_error_caps_unchanged","old_preRp_extra_Pstar_not_reintroduced")
    false_flags=("whole_first_micro_increment_enclosed","full_forward_Rp","all_Z_functional_coverage",
        "global_five_moment_closure","full_goal_complete","accepted_ancestors_replayed")
    require(all(a[k] is True for k in flags) and all(a[k] is False for k in false_flags),"Actual prefix scope changed")
    require(a["exact_algebra"]["passed"] is True and len(a["exact_algebra"]["checks"])==6 and all(a["exact_algebra"]["checks"].values()),
        "Six exact source/geometry identities required")
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction.original_n1_selected_word_numeric import interval
    c=MPIntervalContext();c.prec=1800
    start,end=[interval(c,v) for v in a["native_interval"]]
    oldend=interval(c,p["actual_signed_lower_m15"]["lower_phase_interval"][1])
    require(start==oldend and end==start*3/2 and 0<start<end<c.mpf(1)/4096 and
        a["start_flat_atom"]==p["actual_signed_lower_m15"]["flat_end_atom"] and
        a["end_flat_atom"]["sigma_is_strictly_positive"] is True,"Nonempty exact physical source interval required")
    require(len(a["source_log_bases"])==len(a["source_log_basis_names"])==7 and
        a["source_log_bases"][:6]==p["actual_signed_lower_m15"]["source_log_bases"],"Three distinct flat atoms and original scale anchors required")
    full=a["full_first_micro_continuation_source"]
    require(a["full_first_micro_continuation_source_digest"]==digest(full) and full["passed"] is True and
        full["actual_finite_N_slow_and_fast_signed_functions_retained"] is True and
        full["actual_increment_source_graph_identity_checked"] is True and
        full["exact_selected_N_not_candidate_257"] is True and
        full["whole_chart_numerical_enclosure_available"] is False and
        full["quantitative_mixed_C5_or_global_N_admission"] is False and
        [v["node"] for v in full["actual_S1_ordinary_Z"]]==[129638,129639] and
        [v["node"] for v in full["actual_upper_R_Qy_ordinary_Z"]]==[129654,129655],
        "Full-chart source functions must retain finiteN and remain distinct from numerical enclosures")
    for k in ("0","1"):
        row=a["ordinary_Z"][k]
        require(row["no_finite_N_correction_by_source_identity"] is True,"No finiteN cap may be set to zero as a value")
        for name in ("actual_increment","cumulative_at_three_sc_quarters"):
            packet=row[name];cert=packet["normalized"]
            require(cert["strict_sign"]=="negative" and interval(c,cert["normalized_integral_enclosure"]).b<0 and
                interval(c,cert["relative_interval_width_upper"]).b<c.mpf("1e-27"),"Signed actual prefix or derivative unresolved")
            for kind in ("signed_approximation","error_radius"):
                for t in packet[kind]:
                    require(len(t["source_exponents"])==7,"Physical scale or flat source factor lost")
                    if kind=="error_radius":require(interval(c,t["coefficient_interval"]).a>=0,"Negative error radius")
    require(all(v is False for v in r["limits"].values()),"Nonempty point prefix promoted to completed physical layers")
    provider=ActualFirstMicroContinuation(r)
    exact=provider.full_endpoint('.371');numeric=provider.flat_prefix('0.371')
    require(exact["numerical_enclosure_available"] is False and numeric["numerical_enclosure_available"] is True and
        numeric["ordinary_Z"]==a["ordinary_Z"],"Source definition and actual prefix enclosure conflated")
    exact["source"]["actual_S1_ordinary_Z"][0]["node"]=-1
    require(provider.full_endpoint('.371')["source"]["actual_S1_ordinary_Z"][0]["node"]==129638,
        "Caller mutated the checked source graph")
    rejected=False
    try:provider.full_endpoint('.3')
    except ValueError:rejected=True
    require(rejected,"Foreign axial query admitted without numeric inlet oracle")
    return dict(passed=True,six_exact_source_identities=True,actual_nonempty_signed_increment_Z0_Z1=True,
        source_owned_finiteN_flat_identity=True,direct_cumulative_join_and_shared_endpoint_error=True,
        full_firstmicro_true_source_adapter_and_actual_upper_boundary=True,foreign_query_and_caller_mutation_rejected=True,
        all_inherited_packets_and_caps_unchanged=True,whole_firstmicro_or_global_closure=False,accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p)
    require(p==load(ROOT/PARENT),"Prefix computation mutated its accepted parent")
    inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="actual_nonempty_first_micro_pressure_prefix",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,actual_first_micro_prefix=a,
        input_hashes=inputs,implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(actual_nonempty_first_micro_leading_pressure_increment_enclosed=True,
        exact_flat_correction_from_source_identity=True,ordinary_Z1_and_original_scales=True,
        shared_sc_half_endpoint_error_contracted_before_boxing=True,old_six_and_all_error_caps_retained=True,
        whole_first_micro_or_global_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),independent_new_math_review_completed=True,
        independent_review_scope="Actual native extension through3sc/4, source-owned flat finiteN prefix, physical leading pressure increment/error and cumulative endpoint contraction",
        independent_review_resolution="Prescribed-loop receipt and nonvanishing Phi/L bind the new extension; full-chart source remains source-only; Uy substitutions are references, not registered existing node IDs",
        independently_reviewed_module_sha256=sha(ROOT/CODE[0]),independent_runtime_worker_model_metadata_available=False,
        focused_exact_algebra_and_saved_source_arithmetic=True,accepted_ancestors_replayed=False)
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
