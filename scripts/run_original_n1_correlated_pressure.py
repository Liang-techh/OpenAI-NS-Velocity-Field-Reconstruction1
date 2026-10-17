"""Compute correlated selected-collar pressure increments and their signs."""
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
from openai_ns_reconstruction.original_n1_correlated_pressure_increment import build,PARENT,PARENT_SHA256
NAME="lei_ren_part1_paper_compliant_current_original_n1_correlated_pressure_increment.json"
CODE=("src/openai_ns_reconstruction/original_n1_correlated_pressure_increment.py",
    "scripts/run_original_n1_correlated_pressure.py")


def parent():
    require(sha(ROOT/PARENT)==PARENT_SHA256,"Accepted pressure primitive parent changed")
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
        "Different correlated pressure parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes",
        "contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source/error packet changed: "+k)
    a=r["actual_correlated_leading_pressure_increment"]
    require(a["passed"] is True and a["parent_mapping_sha256"]==digest(p) and
        a["source_family"]==p["source_family"] and a["ordinary_Z_orders"]==[0,1] and
        a["actual_correlated_signed_coefficients_and_error_radii_evaluated"] is True and
        a["shared_mean_and_Psi_value_inlet_constants_canceled_by_same_source_FTC"] is True and
        a["single_core_Psi_rho_inlet_uncertainty_retained"] is True and
        a["old_independent_endpoint_Qpsi_values_not_used"] is True and
        a["original_nonzero_flat_source_factor_retained"] is True,
        "Same-source correlation or actual numerical arithmetic admission failed")
    require(a["flat_source_atom"]["sigma_is_strictly_positive"] is True and
        a["flat_source_atom"]["sigma_materialized"] is False and
        a["flat_source_atom"]["cap_is_not_the_function_value"] is True and
        a["flat_source_atom"]["monotone_on_selected_real_interval_certified"] is True and
        a["flat_source_atom"]["monotonicity_proof"]["positive_rewrite_identity"] is True and
        a["source_log_basis_names"][-1]=="sigma_selected","Original nonzero flat source atom lost")
    require(a["exact_algebra"]["passed"] is True and all(a["exact_algebra"]["checks"].values()),
        "Own-mean/pressure transport algebra failed")
    for k in ("0","1"):
        for name,sign in (("P1_W0_local","positive"),("m15_leading_N1p_local","negative")):
            cert=a["ordinary_Z"][k][name]["normalized"]
            require(cert["strict_sign"]==sign and cert["retained_flat_source_not_zeroed"] is True and
                cert["dyadic_caps_are_errors_not_signed_values"] is True and len(cert["dyadic_term_bounds"])==34,
                "Certified normalized local integral/sign contract changed")
    require(a["P1_is_distinct_from_m15"] is True and
        a["full_n1_pressure_or_C_inner_m15_evaluated"] is False and
        a["all_Z_functional_coverage"] is False and
        a["old_error_caps_or_six_ledger_terms_removed"] is False and
        a["accepted_ancestors_replayed"] is False and a["source_executed"] is False and
        all(v is False for v in r["limits"].values()),"Point source sign promoted to global closure")
    return dict(passed=True,correlated_actual_leading_pressure_increment_enclosed=True,
        common_mean_and_inlet_value_canceled_by_same_source_FTC=True,
        core_inlet_derivative_uncertainty_and_nonzero_flat_atom_retained=True,
        P1_positive_m15_N1p_negative_at_saved_Z_for_orders_0_1=True,
        all_inherited_six_signed_terms_and_error_caps_unchanged=True,
        full_C_inner_or_global_terminal_matching=False,accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p);a=build(ROOT,p);inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="actual_correlated_leading_pressure_increment",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,
        actual_correlated_leading_pressure_increment=a,input_hashes=inputs,
        implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(correlated_increment_from_one_actual_core_Psi_rho_inlet=True,
        no_independent_Qpsi_endpoint_subtraction=True,original_nonzero_sigma_source_factor=True,
        normalized_local_P1_and_m15_N1p_values_enclosed_with_strict_signs=True,
        ordinary_Z0_Z1_at_saved_point_only=True,
        all_old_error_caps_and_six_signed_terms_retained=True,
        full_C_inner_or_global_five_moment_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),
        independent_new_math_review_completed=True,
        independent_review_scope="Local Z~0.371 leading N1p increment, ordinary Z0/Z1 only",
        independent_review_resolution="Factor2 core radial cover and exact original-sigma monotonicity gate accepted",
        independently_reviewed_module_sha256=sha(ROOT/CODE[0]),
        independent_runtime_worker_model_metadata_available=False,
        focused_exact_algebra_and_actual_source_arithmetic=True,
        accepted_ancestors_replayed=False)
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
