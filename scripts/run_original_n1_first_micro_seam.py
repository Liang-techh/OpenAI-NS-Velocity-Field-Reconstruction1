"""Save or inspect the actual first-micro sc/2 viscous seam contraction."""
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
from openai_ns_reconstruction.original_n1_first_micro_seam import build,PARENT,PARENT_SHA256
NAME="lei_ren_part1_paper_compliant_current_original_n1_first_micro_seam.json"
CODE=("src/openai_ns_reconstruction/original_n1_first_micro_seam.py","scripts/run_original_n1_first_micro_seam.py")


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
    require(r["parent_report"]==Path(PARENT).name and r["parent_report_sha256"]==PARENT_SHA256,"Different actual q2 parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes","contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source/error packet changed: "+k)
    a=r["actual_first_micro_seam"]
    require(a["passed"] is True and a["source_family"]==p["source_family"] and
        a["ordinary_Z_orders"]==[0,1] and a["selected_Z_interval"]==p["actual_signed_lower_m15"]["selected_Z_interval"] and
        a["exact_phase"]==p["actual_signed_lower_m15"]["lower_phase_interval"][1] and
        a["accepted_q2_inlet_digest"]==p["actual_correlated_q2_inlet"]["shared_q2_inlet_digest"],
        "Same actual sc/2 source and ordinary Z derivative required")
    flags=("same_error_function_at_shared_endpoint","endpoint_contraction_before_interval_boxing",
        "lower_leading_Vy_kept_as_actual_source","leading_Qaxis_and_own_mean_not_reset",
        "old_six_numeric_terms_and_all_error_caps_unchanged","new_representation_only","one_saved_Z_tube_only")
    false_flags=("downstream_actual_S1_evaluated","first_micro_upper_R_Qy_evaluated","full_forward_Rp",
        "all_Z_functional_closure","global_five_moment_closure","full_goal_complete","accepted_ancestors_replayed")
    require(all(a[k] is True for k in flags) and all(a[k] is False for k in false_flags),"Seam scope or error policy changed")
    require(a["exact_algebra"]["passed"] is True and len(a["exact_algebra"]["checks"])==6 and
        all(a["exact_algebra"]["checks"].values()) and a["lower_R_Qy_source_nodes"]==[30442,30443] and
        a["lower_actual_Vy_source_nodes"]==[30363,30364],"Six seam identities and true finite endpoints required")
    require(a["archived_preRp_endpoint_has_extra_Pstar"] is True and
        a["archived_preRp_m15_is_not_corrected_physical_source"] is True and
        a["actual_graph_scale_defect_checks"]["passed"] is True and
        len(a["actual_graph_scale_defect_checks"]["checks"])==6 and
        all(a["actual_graph_scale_defect_checks"]["checks"].values()) and
        a["corrected_physical_endpoint_graph_digest"]==digest(a["corrected_physical_endpoint_graph"]),
        "Physical V_y normalization and real graph repair must be checked, not inferred")
    patch=a["corrected_physical_endpoint_graph"]
    require(patch["exact_replacements"]=={"30365":30363,"30366":30364} and
        patch["old_graph_unchanged"] is True and len(patch["appended_nodes"])==30,
        "Correct only the extra scale on the actual lower V_y rows")
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction.original_n1_selected_word_numeric import interval,endpoint
    c=MPIntervalContext();c.prec=1800
    require(interval(c,a["positive_inertia_radius_ratio_upper"]).a>0,"Positive inherited-density geometry extension required")
    for k in ("0","1"):
        row=a["ordinary_Z"][k]
        require(row["strict_sign_claimed"] is False and row["caps_not_selected_as_function_values"] is True and
            row["core_viscous_and_collar_derivative_endpoint_not_boxed"] is True and
            row["common_unit"]=="reference*A*Rv^2","Only shifted-source model plus radius admitted")
        for kind in ("signed_approximation","error_radius"):
            for t in row[kind]:
                require(len(t["source_exponents"])==6,"Source basis truncated")
                if kind=="error_radius":require(interval(c,t["coefficient_interval"]).a>=0,"Negative error radius")
        logs=row["logarithmic_triangle_terms"]
        upper=max(endpoint(c,interval(c,t["log_absolute_upper"])) for t in logs)+c.ln(len(logs))
        require(upper._mpi_==interval(c,row["absolute_common_unit_log_upper"])._mpi_,"Common-unit bound is not the saved triangle upper")
    require(all(v is False for v in r["limits"].values()),"Source seam promoted to completed physical layers")
    return dict(passed=True,six_exact_chain_rule_and_FTC_identities=True,six_real_graph_scale_defect_identities=True,
        archived_extra_Pstar_fixed_in_new_graph=True,ordinary_Z0_Z1_retained=True,
        original_pressure_graph_and_actual_Vy_source_bound=True,shared_endpoint_contracted_before_boxing=True,
        all_inherited_packets_and_error_caps_unchanged=True,first_micro_increment_or_global_closure=False,
        accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p)
    require(p==load(ROOT/PARENT),"Seam computation mutated its accepted parent")
    inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="actual_first_micro_sc_half_viscous_seam",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,actual_first_micro_seam=a,
        input_hashes=inputs,implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(actual_q2_shifted_inlet_evaluated=True,source_owned_lower_viscous_endpoint_contracted=True,
        actual_raw_V_y_is_physical=True,extra_Pstar_removed_from_new_endpoint_graph=True,
        same_actual_leading_field_and_error_function=True,nonzero_leading_Vy_and_histories_retained=True,
        old_six_and_all_error_caps_retained=True,downstream_or_global_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),independent_new_math_review_completed=True,
        independent_review_scope="Raw physical Vy normalization, 30-node corrected graph, native sc/2 field, flat finite-N inlet, RQy chain rule, pressure FTC and shifted inlet model/error",
        independent_review_resolution="Raw Vy=Uy=epsilon*rho*Psi_rho; Qy=epsilon*rho*b_rho; archived extra Pstar corrected only in appended graph; old results not equated",
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
