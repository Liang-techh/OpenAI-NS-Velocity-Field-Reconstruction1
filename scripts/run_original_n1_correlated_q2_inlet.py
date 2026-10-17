"""Enclose signed actual A-A0 and supply a correlated q2 inlet."""
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
from openai_ns_reconstruction.original_n1_correlated_q2_inlet import build, PARENT, PARENT_SHA256, ActualQ2Inlet, SIX_WEIGHTS, SIX_NAMES, REDUCED_WEIGHTS, reduce_weights
NAME = "lei_ren_part1_paper_compliant_current_original_n1_correlated_q2_inlet.json"
CODE = ("src/openai_ns_reconstruction/original_n1_correlated_q2_inlet.py", "scripts/run_original_n1_correlated_q2_inlet.py")


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
    require(r["parent_report"]==Path(PARENT).name and r["parent_report_sha256"]==PARENT_SHA256,"Different actual q2 parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes",
             "contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source/error packet changed: "+k)
    a=r["actual_correlated_q2_inlet"]
    require(a["passed"] is True and a["parent_mapping_sha256"]==digest(p) and a["source_family"]==p["source_family"] and
        a["ordinary_Z_orders"]==[0,1] and a["signed_remainder_normalization"]=="hb*epsilon^2*Pstar^2*sigma(sc/2)",
        "Original physical q2/remainder source required")
    true_flags=("signed_A_minus_A0_value_and_ordinary_Z1_enclosed","exact_six_source_weights_reduce_to_accepted_q2",
        "shared_remainder_error_canceled_before_aggregate_boxing","original_nonzero_sigma_not_replaced_by_cap",
        "old_six_and_all_error_caps_unchanged","one_point_only")
    false_flags=("old_finite_numeric_prefix_equals_exact_function","downstream_first_micro_increment_evaluated",
        "forward_Rp_discrepancy_enclosed","all_Z_functional_coverage","global_five_moment_closure","full_goal_complete",
        "accepted_ancestors_replayed","source_executed")
    require(all(a[k] is True for k in true_flags) and all(a[k] is False for k in false_flags),"Point q2 source scope failed")
    require(a["exact_algebra"]["passed"] is True and all(a["exact_algebra"]["checks"].values()) and
        a["actual_primitive_binding"]["packet_digest"]==digest(p["canonical_collar_signed_primitive"]) and
        a["actual_primitive_binding"]["radial_Psi_order"]==1,"Actual radial1 signed primitive/source algebra changed")
    packet=a["shared_q2_inlet"]
    require(a["shared_q2_inlet_digest"]==digest(packet) and
        packet["exact_phase_tuple"]==p["actual_signed_lower_m15"]["lower_phase_interval"][1]["lower_exact_mpf_tuple"] and
        packet["q_scale_sc_quarter_tuple"]==p["actual_signed_lower_m15"]["lower_phase_interval"][0]["lower_exact_mpf_tuple"] and
        packet["source_leaves_digest"]==digest({k:v["source_leaves"] for k,v in a["ordinary_Z"].items()}) and
        packet["component_four_exact_pressure_partition_digest"]==digest(p["actual_collar_pressure_partition"]) and
        [v["name"] for v in packet["components"]]==list(SIX_NAMES) and
        [v["source_weights"] for v in packet["components"]]==list(SIX_WEIGHTS) and
        [v["error_leaf_weights"] for v in packet["components"]]==list(SIX_WEIGHTS) and
        reduce_weights([v["source_weights"] for v in packet["components"]])==REDUCED_WEIGHTS and
        packet["source_weights_reduce_before_boxing"] is True and
        packet["do_not_sum_independent_component_error_radii"] is True and
        packet["legacy_six_entry_format"] is False,"Shared signed source/error leaves must cancel before boxing")
    from mpmath.ctx_iv import MPIntervalContext
    from openai_ns_reconstruction.original_n1_selected_word_numeric import interval
    c=MPIntervalContext();c.prec=1800
    for k in ("0","1"):
        row=a["ordinary_Z"][k];cert=row["actual_signed_A_minus_A0"]["normalized"]
        require(cert["strict_sign"]=="negative" and interval(c,cert["normalized_integral_enclosure"]).b<0 and
            interval(c,cert["relative_interval_width_upper"]).b<c.mpf("1e-27") and
            cert["caps_not_selected_as_values"] is True and
            row["aggregate_evaluated_after_source_weight_reduction"] is True and
            row["remainder_radius_subtracted_as_value"] is False and
            row["accepted_joined_q2_enclosure"]==p["actual_signed_lower_m15"]["ordinary_Z"][k]["actual_inner_plus_lower_m15"],
            "Actual signed remainder or accepted aggregate changed")
        for leaf in row["source_leaves"].values():
            count=leaf["source_log_basis_count"]
            require(count in (5,6) and all(len(t["source_exponents"])==count
                for key in ("signed_approximation","error_radius") for t in leaf[key]),"Mixed source basis silently truncated")
    provider=ActualQ2Inlet(r)
    first=provider.first_micro_inlet(".371");pre=provider.preRp_inner_lower("0.371")
    for packet in (first,pre):
        packet.pop("consumer");packet.pop("Z")
    require(first==pre,"Downstream consumers received different q2 sources")
    require(first["source_leaves_by_ordinary_Z"]=={k:v["source_leaves"] for k,v in a["ordinary_Z"].items()},
        "Downstream shared source/error payload changed")
    foreign_rejected=False
    try:provider.first_micro_inlet(".3")
    except ValueError:foreign_rejected=True
    require(foreign_rejected,"Foreign point admitted without a source oracle")
    first["components"][0]["name"]="mutated caller copy"
    first["source_leaves_by_ordinary_Z"]["0"]["remainder"]["actual_function_binding"]["name"]="mutated caller leaf"
    require(provider.first_micro_inlet(".371")["components"][0]["name"]==SIX_NAMES[0],"Caller mutated the shared inlet")
    require(provider.first_micro_inlet(".371")["source_leaves_by_ordinary_Z"]["0"]["remainder"]["actual_function_binding"]["name"]=="remainder",
        "Caller mutated the shared source/error leaf")
    require(all(v is False for v in r["limits"].values()),"Source inlet result promoted to completed physical layers")
    return dict(passed=True,actual_signed_A_minus_A0_Z0_Z1_enclosed=True,shared_six_source_and_error_leaves_reduced=True,
        both_new_consumers_share_actual_q2_inlet=True,foreign_query_and_caller_mutation_rejected=True,
        all_inherited_six_numeric_terms_and_error_caps_unchanged=True,downstream_or_full_forward_closure=False,
        accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p)
    require(p==load(ROOT/PARENT),"New inlet computation mutated its accepted parent")
    inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="actual_signed_A_minus_A0_and_correlated_q2_inlet",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,actual_correlated_q2_inlet=a,
        input_hashes=inputs,implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(actual_signed_A_minus_A0_integral_and_ordinary_Z1_enclosed=True,
        principal_and_remainder_are_shared_function_leaves=True,source_weights_reduce_before_aggregate_boxing=True,
        actual_q2_inlet_shared_by_new_first_micro_and_preRp_consumers=True,old_numeric_six_format_distinct=True,
        old_six_and_all_error_caps_retained=True,full_forward_or_global_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),independent_new_math_review_completed=True,
        independent_review_scope="Actual signed A-A0 at sc/2, positive d-sigma error bound and shared q2 source/error-leaf reduction",
        independent_review_resolution="Phase tuple is S=sc/2; sc/4 is a separate q scale; provider verifies both plus shared physical source/error leaf digests",
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
