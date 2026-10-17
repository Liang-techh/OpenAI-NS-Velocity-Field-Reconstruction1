"""Save or inspect true first/second-micro pressure in the six-I1 basis."""
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
from openai_ns_reconstruction.original_n1_micro_correlated_m15 import build,PARENT,PARENT_SHA256,HELPER,DEFECT
NAME="lei_ren_part1_paper_compliant_current_original_n1_micro_correlated_m15.json"
CODE=("src/openai_ns_reconstruction/original_n1_micro_correlated_m15.py",
      "scripts/run_original_n1_micro_correlated_m15.py",HELPER,DEFECT)


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
    a=r["actual_micro_correlated_m15"]
    require(a["passed"] is True and a["source_family"]==p["source_family"] and a["ordinary_Z_orders"]==[0,1] and
        a["shifted_inlet_source_digest"]==digest(p["actual_first_micro_seam"]),"Same actual source and shifted inlet required")
    flags=("actual_full_Q4_six_I1_source_implemented","bulk_Qy_and_Uy_eliminated_by_exact_IBP",
        "full_slow_and_fast_finiteN_functions_retained_at_outer_endpoint",
        "internal_RQy_and_reduced_B_seams_cancel_as_functions",
        "same_source_error_function_shared_before_any_future_boxing","old_six_numeric_terms_and_all_error_caps_unchanged")
    false_flags=("accepted_seam_fingerprints_replayed","numerical_signed_joint_integral_evaluated",
        "numerical_full_micro_pair_enclosure","all_Z_five_moment_closure","full_forward_Rp","full_goal_complete","accepted_ancestors_replayed")
    require(all(a[k] is True for k in flags) and all(a[k] is False for k in false_flags),"Source-only scope changed")
    require(a["exact_algebra"]["passed"] is True and len(a["exact_algebra"]["checks"])==6 and
        all(a["exact_algebra"]["checks"].values()),"Six exact pressure/source/seam identities required")
    seam=a["seam_certificate"]
    require(seam["passed"] is True and len(seam["source_identity_rows"])==14 and
        seam["source_graph_prefix_and_all_14_field_ports_checked"] is True and
        seam["free_Z_not_numeric_Z_substitution"] is True and seam["accepted_fingerprint_computation_replayed"] is False,
        "Accepted actual free-Z seam must be consumed without replay")
    patch=a["graph_patch"]
    require(digest(patch)==a["graph_patch_digest"] and patch["old_graph_unchanged"] is True and
        patch["exactly_one_representation_consumed"] is True and patch["output_formula"]=="M15(s=2)=H+shifted_inlet_tail",
        "Append-only exact pressure representation required")
    require(set(patch["local_I1"])=={"first_micro","second_micro"} and
        set(patch["joint_I1"])==set(a["six_I1_names"]) and len(a["six_I1_names"])==6,"Six integrands per actual chart required")
    base=patch["base_graph_node_count"];nodes=patch["appended_nodes"];integrals=[];binders=[]
    for index,row in enumerate(nodes,base):
        if row["operation"]=="bound_variable":binders.append(row["name"])
        if row["operation"]=="definite_integral":
            require(row["physical_radius_power"]==2 and row["physical_Jacobian_applied_once"] is True and
                row["no_total_y_derivative_in_volume"] is True and row["endpoint_independent_of_Z"] is True and
                row["ordinary_slow_Z_derivative_order"] in (0,1),"Wrong physical measure or derivative order")
            integrals.append(index)
    require(len(integrals)==24 and len(binders)==len(set(binders))==24 and
        set(integrals)=={n for chart in patch["local_I1"].values() for pair in chart.values() for n in pair},
        "Exactly24 new signed I1 leaves and unique lexical binders required")
    for defect in patch["finite_N_signed_integrand_defects"].values():
        require(defect["passed"] is True and defect["symbolic_checks"]["passed"] is True and
            all(defect["exact_symbolic_expansion_checks"].values()) and all(defect["ordinary_Z1_checks"].values()) and
            len(defect["U0"])==len(defect["Q0"])==3 and set(defect["defect_integrands"])==set(a["six_I1_names"]) and
            defect["mean_Z3_input_ref"]==defect["dm"][3] and defect["whole_micro_numeric"] is False and
            defect["flags"]["no_Z3_baseline_claim"] is True and
            defect["uniform_in_N_bound"] is False and defect["cap_as_value"] is False,
            "Same-source raw six-integrand differences and QZZ/meanZZZ gap required")
    for name in ("actual_micro_pair_increment","shifted_inlet_tail","reduced_S1","reduced_boundary","reduced_volume"):
        require(len(patch[name])==2 and all(type(i) is int and base<=i<base+len(nodes) for i in patch[name]),"Unregistered ordinary pair")
    require(all(v is False for v in r["limits"].values()),"Source representation promoted to completed physical layer")
    return dict(passed=True,six_exact_pressure_identities=True,twenty_four_true_I1_leaves=True,
        exact_same_function_internal_seam=True,ordinary_Z1_and_physical_measure=True,
        signed_raw_six_integrand_defect_and_mean_Z3_dependency=True,
        true_outer_RQy_retained=True,all_inherited_packets_and_caps_unchanged=True,
        whole_micro_pair_or_global_closure=False,accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p)
    require(p==load(ROOT/PARENT),"Correlated construction mutated its accepted parent")
    inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="true_micro_pair_correlated_six_I1_source",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,actual_micro_correlated_m15=a,
        input_hashes=inputs,implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(true_full_Q4_six_I1_source=True,bulk_total_y_derivatives_eliminated=True,
        accepted_free_Z_seam_consumed_without_replay=True,outer_slow_fast_boundary_retained=True,
        old_six_and_all_error_caps_retained=True,numerical_whole_micro_or_global_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),independent_new_math_review_completed=True,
        independent_review_scope="True Q4 six-I1 pressure IBP, ordinaryZ1 graph arithmetic, saved free-Z endpoint identity and retained physical outer RQy",
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
