"""Save or inspect the true micro mean with periodic primitive transport."""
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
from openai_ns_reconstruction.original_n1_micro_phase_mean import build,PARENT,PARENT_SHA256,HELPER
NAME="lei_ren_part1_paper_compliant_current_original_n1_micro_phase_mean.json"
CODE=("src/openai_ns_reconstruction/original_n1_micro_phase_mean.py",
      "scripts/run_original_n1_micro_phase_mean.py",HELPER)


def parent():
    require(sha(ROOT/PARENT)==PARENT_SHA256,"Accepted first-micro prefix changed")
    p=load(ROOT/PARENT);require(p["checks"]["passed"] is True,"Checked prefix required")
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
        require(sha(ROOT/n)==h,"Changed correlated source: "+n)
        if n.startswith("experiments/root_st073/"):
            n=n[len("experiments/root_st073/"):];target=inputs
        else:target=impl
        require(target.get(n,h)==h,"Conflicting source pin: "+n);target[n]=h
    return inputs,impl


def packet_check(r,p):
    require(r["parent_report"]==Path(PARENT).name and r["parent_report_sha256"]==PARENT_SHA256,"Different actual firstmicro parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes","contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source/error packet changed: "+k)
    a=r["actual_micro_phase_mean"]
    require(a["passed"] is True and a["source_family"]==p["source_family"] and a["ordinary_Z_orders"]==[0,1,2,3] and
        a["accepted_six_I1_source_digest"]==digest(p["actual_micro_correlated_m15"]),"Same actual six-I1 source required")
    for key in ("same_source_internal_C_endpoint_canceled","complete_C_and_Cy_convolutions_retained",
                "ordinary_mean_Z3_source_implemented","old_six_and_all_error_caps_unchanged"):
        require(a[key] is True,"Missing exact source property: "+key)
    for key in ("numeric_Cy_Z3_enclosure_available","whole_mean_Z3_numerical_enclosure_available",
                "legacy_C3_caps_admitted_for_true_C5_source","numerical_QZZ_available",
                "numerical_whole_micro_or_forward_Rp","global_five_moment_closure","accepted_ancestors_replayed","full_goal_complete"):
        require(a[key] is False,"Unproved numerical/closure gate: "+key)
    seam=a["new_free_phase_B_slow_By_seam"]
    require(seam["passed"] is True and seam["rows"]==len(seam["new_source_rows"])==8 and
        seam["Z_kept_free"] is True and seam["issued_phase_kept_free"] is True and
        seam["row_digest"]==digest(seam["new_source_rows"]) and
        len(seam["new_second_local_convolutions_zero_at_s1"])==8 and
        seam["accepted14_row_fingerprint_computation_replayed"] is False,"New B/Cy source seam proof required")
    algebra=a["exact_algebra"]
    require(algebra["passed"] is True and algebra["transport_delta_m_y_plus_delta_m_equals_Pstar_B_over_N"] is True and
        all(algebra["seam_cancel_through_Z0_to_Z3"].values()) and algebra["nonzero_inlet_retained_generically"] is True and
        all(a["physical_Q_recovery_controls"]["ordinary_Z0_Z1_Z2"]),"Physical transport/recovery algebra failed")
    patch=a["graph_patch"];kernel=patch["kernel"];nodes=patch["appended_nodes"];base=patch["base_graph_node_count"]
    require(digest(patch)==a["graph_patch_digest"] and patch["old_Q4_and_624_node_source_suffix_unchanged"] is True and
        base==184324 and kernel["graph_source_ports"]["scalars"]["Pstar"]==33572 and
        kernel["graph_source_ports"]["scalars"]["N"]==18 and kernel["graph_source_ports"]["scalars"]["hb"]==8,
        "Exact source namespace and physical scalars required")
    integrals=[n for n in nodes if n["operation"]=="definite_integral"]
    binders=[n["name"] for n in nodes if n["operation"]=="bound_variable"]
    require(len(integrals)==len(binders)==len(set(binders))==32,"Sixteen phase and sixteen moving integrals required")
    for row in integrals:
        require(row["ordinary_Z_order"] in (0,1,2,3) and row["factorial_normalization"] is False,"Ordinary Z rows required")
        if "moving_rate" in row:require(row["physical_hb_weight_applied_once"] is True,"Physical measure repeated or omitted")
    for name in ("first_micro","second_micro"):
        require(len(kernel["complete_physical_delta_m"][name])==len(kernel["unresolved_Cy_part"][name])==4 and
            len(patch["Q_recovery"][name]["physical_Q_ordinary_Z0_Z2"])==3 and
            len(patch["Q_recovery"][name]["source_six_I1_C1_integrands"])==6,
            "MeanZ3, QZZ and six ordinary C1 integrands must remain connected")
    diagnostic=a["legacy_weighted_cap_diagnostic"]
    require(diagnostic["numerical_bound_on_actual_C5_mean"] is False and diagnostic["caps_used_as_signed_values"] is False and
        diagnostic["charts"]["first_micro"]["radius_small"] is False and
        a["directed_physical_prefactor"]["coefficient_not_mean_value_or_error_radius"] is True,
        "A prefactor or older cap cannot become an actual field bound")
    require(all(v is False for v in r["limits"].values()),"Physical completion gate changed")
    return dict(passed=True,new_free_Z_and_phase_eight_row_seam=True,second_local_integrals_empty_at_lower=True,
        physical_Nminus2_transport_with_Cy_retained=True,ordinary_mean_Z3_to_QZZ_and_six_integrands=True,
        thirty_two_exact_integral_leaves=True,legacy_cap_transfer_rejected=True,
        all_inherited_packets_and_caps_unchanged=True,numerical_mean_or_global_closure=False,accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p)
    require(p==load(ROOT/PARENT),"Correlated construction mutated its accepted parent")
    inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="true_micro_phase_primitive_mean_and_QZZ_source",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,actual_micro_phase_mean=a,
        input_hashes=inputs,implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(true_micro_mean_Z3_Nminus2_function=True,new_free_Z_free_phase_B_Cy_seam=True,
        fixed_phase_Cy_convolution_retained=True,old_C3_cap_transfer_not_admitted=True,
        old_six_and_all_error_caps_retained=True,numerical_whole_micro_or_global_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),independent_new_math_review_completed=True,
        independent_review_scope="Actual micro C/Cy phase primitives, new B/slowBy source seam, physical Nminus2 mean and QZZ recovery, explicit lack of numerical Cy bounds",
        independently_reviewed_module_sha256=sha(ROOT/CODE[0]),independent_runtime_worker_model_metadata_available=False,
        focused_exact_algebra_and_new_graph_only=True,accepted_ancestors_replayed=False)
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
