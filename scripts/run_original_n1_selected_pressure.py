"""Compute signed selected-collar leading-pressure primitives from saved sources."""
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
from openai_ns_reconstruction.original_n1_selected_pressure_primitive import build,PARENT,PARENT_SHA256

NAME="lei_ren_part1_paper_compliant_current_original_n1_selected_pressure_primitive.json"
CODE=("src/openai_ns_reconstruction/original_n1_selected_pressure_primitive.py",
    "scripts/run_original_n1_selected_pressure.py")


def parent():
    require(sha(ROOT/PARENT)==PARENT_SHA256,"Accepted selected-word parent changed")
    p=load(ROOT/PARENT);require(p["checks"]["passed"] is True,"Passed source parent required")
    return p


def pins(r):
    count=0
    for base,key in ((DIRECTORY,"input_hashes"),(ROOT,"implementation_hashes")):
        for name,h in r[key].items():
            require(sha(base/name)==h,"Changed source pin: "+name);count+=1
    return count


def merge_pins(p,a):
    inputs,impl=dict(p["input_hashes"]),dict(p["implementation_hashes"])
    for row in a["source_pins"].values():
        n,h=row["path"],row["expected_sha256"]
        require(sha(ROOT/n)==h,"Changed new source: "+n)
        if n.startswith("experiments/root_st073/"):
            n=n[len("experiments/root_st073/"):];target=inputs
        else:target=impl
        require(target.get(n,h)==h,"Conflicting source pin: "+n);target[n]=h
    return inputs,impl


def packet_check(r,p):
    require(r["parent_report"]==Path(PARENT).name and r["parent_report_sha256"]==PARENT_SHA256,
        "Different pressure primitive parent")
    mutable={"stage","parent_report","parent_report_sha256","input_hashes","implementation_hashes",
        "contract","checks","execution_seconds"}
    for k,v in p.items():
        if k not in mutable:require(r[k]==v,"Inherited source or error changed: "+k)
    a=r["actual_selected_leading_pressure_primitive"]
    require(a["passed"] is True and a["source_family"]==p["source_family"] and
        a["parent_mapping_sha256"]==digest(p) and
        a["actual_signed_endpoint_and_axis_approximation_evaluated"] is True and
        a["actual_leading_N1p_integral_enclosed"] is True and
        a["full_W0_pressure_increment_enclosed"] is True and
        a["physical_N1p_no_Pstar_normalization_used"] is True and
        a["radial_constant_Qaxis_canceled_as_same_function_before_boxing"] is True and
        a["saved_endpoint_source_uncertainty_retained"] is True and
        a["geometry_difference_factored_before_endpoint_boxing"] is True and
        a["ordinary_Z_orders"]==[0,1] and
        a["source_executed"] is False and a["accepted_ancestors_replayed"] is False,
        "Actual signed primitive/source arithmetic admission failed")
    require(a["exact_algebra"]["passed"] is True and
        len(a["exact_algebra"]["checks"])==9 and all(a["exact_algebra"]["checks"].values()),
        "Physical FTC, epsilon split or geometry algebra failed")
    require(a["P1_is_distinct_from_m15"] is True and a["m15_component_is_only_leading_N1p"] is True and
        a["full_C_inner_m15_evaluated"] is False and a["old_error_caps_removed"] is False and
        a["six_term_ledger_modified"] is False and a["all_Z_functional_coverage"] is False and
        a["endpoint_sign_or_useful_relative_accuracy_certified"] is False and
        all(v is False for v in r["limits"].values()),"Local primitive scope was promoted")
    return dict(passed=True,inherited_source_graphs_six_signed_terms_and_all_errors_unchanged=True,
        signed_pressure_approximation_and_inertial_radius_separate=True,
        P1_and_m15_have_distinct_sign_and_epsilon_scaling=True,
        same_source_Qaxis_and_geometry_correlation_used_before_boxing=True,
        full_C_inner_or_global_closure=False,accepted_ancestors_replayed=False)


def run():
    started=time.perf_counter();p=parent();pins(p)
    a=build(ROOT,p);inputs,impl=merge_pins(p,a)
    r=copy.deepcopy(p);r.update(stage="actual_selected_leading_pressure_signed_primitive",
        parent_report=Path(PARENT).name,parent_report_sha256=PARENT_SHA256,
        actual_selected_leading_pressure_primitive=a,input_hashes=inputs,
        implementation_hashes=dict(impl,**{n:sha(ROOT/n) for n in CODE}))
    r["contract"]=dict(actual_signed_leading_pressure_component_enclosed=True,
        P1_W0_is_not_complete_n1_pressure_coefficient=True,
        m15_result_is_leading_N1p_component_only=True,
        exact_Qaxis_cancellation_and_factored_geometry_before_boxing=True,
        remaining_inertial_volume_error_is_not_a_signed_value=True,
        saved_endpoint_errors_and_six_signed_terms_retained=True,
        full_C_inner_m15_or_independent_terminal_matching=False,
        global_five_moment_closure=False,accepted_ancestors_replayed=False)
    r["checks"]=dict(passed=True,packet_invariants=packet_check(r,p),
        independent_new_math_review_completed=False,
        focused_exact_algebra_and_actual_source_arithmetic=True,
        accepted_ancestors_replayed=False)
    r["execution_seconds"]=round(time.perf_counter()-started,3);return r


def inspect(path):
    r=load(path);require(r["checks"]["passed"] is True,"Saved passed report required")
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
