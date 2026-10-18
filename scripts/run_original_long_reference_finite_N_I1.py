"""Produce/inspect actual long/reference native controls and finite-N I1.

Reads the accepted native background report without replaying its builders.
Only new native-control, loop and finite-N stage arithmetic executes.
"""
from __future__ import annotations
import argparse
import copy
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/"experiments/root_st073"
PARENT_NAME="lei_ren_part1_paper_compliant_current_original_long_reference_background_I1.json"
PARENT_SHA256="7fce21dc6e763f42b99ef7fcf7416decaf848a5124cda99741222dba9619eb1b"
NAME="lei_ren_part1_paper_compliant_current_original_long_reference_finite_N_I1.json"
NATIVE="actual_long_reference_native_controls"
LOOP="actual_long_reference_loop_Z3_bounds"
FINITE="actual_long_reference_finite_N_I1_numeric"
CODE=("src/openai_ns_reconstruction/original_native_log_jet4.py",
      "src/openai_ns_reconstruction/original_long_reference_native_controls.py",
      "src/openai_ns_reconstruction/original_long_reference_loop_Z3_bounds.py",
      "src/openai_ns_reconstruction/original_long_reference_finite_N_I1_numeric.py",
      "scripts/run_original_long_reference_finite_N_I1.py")
REVIEWED={CODE[0]:"ba11bb60af654e58271a735560aeae2f1f9320ef3bfb8cee0966e9b7e50e0df4",
          CODE[1]:"b15f020def2316b4e91a38dea0faa6b6a04cc35fa50f1cf6f0b78fa4d2235094",
          CODE[2]:"405f61729f10bd545480fe1446cdf47b11ccdf0cef9e670a08e3d83a1aa6e85e",
          CODE[3]:"dd32288b88e7ccc168ad13cc6d74b6671832bafa40dae35c698ffa3d32a93176"}


def require(value,message):
    if not value:
        raise ValueError(message)


def sha(path):
    d=hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda:stream.read(1<<20),b""):
            d.update(block)
    return d.hexdigest()


def load(path):
    raw=Path(path).read_bytes()
    return json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)


def check_pins(report):
    count=0
    for folder,group in ((BASE,report["input_hashes"]),(ROOT,report["implementation_hashes"])):
        for name,expected in group.items():
            require(isinstance(expected,str) and sha(folder/name)==expected,
                    "Pinned source bytes changed: "+name)
            count+=1
    return count


def merge_pins(accepted,*packets):
    inputs=dict(accepted["input_hashes"])
    implementations=dict(accepted["implementation_hashes"])
    inputs[PARENT_NAME]=PARENT_SHA256
    for packet in packets:
        require(packet["passed"] is True and packet["parent_report_sha256"]==PARENT_SHA256,
                "New packet belongs to another predecessor")
        for key,item in packet["source_pins"].items():
            path=Path(item.get("path",key) if isinstance(item,dict) else key)
            expected=item.get("expected_sha256") if isinstance(item,dict) else item
            if len(path.parts)==1:
                actual=BASE/path
                name,target=path.name,inputs
            elif path.as_posix().startswith("experiments/root_st073/"):
                actual=ROOT/path
                name,target=path.name,inputs
            else:
                actual=ROOT/path
                name,target=path.as_posix(),implementations
            require(sha(actual)==expected and target.get(name,expected)==expected,
                    "New source pin mismatch: "+path.as_posix())
            target[name]=expected
    for path in CODE:
        implementations[path]=sha(ROOT/path)
    return inputs,implementations


def check_records(value):
    from mpmath.libmp import mpf_cmp
    count=0
    if isinstance(value,dict):
        if "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value:
            lo,hi=value["lower_exact_mpf_tuple"],value["upper_exact_mpf_tuple"]
            require(all(isinstance(t,list) and len(t)==4 and all(type(n) is int for n in t)
                        and t[0] in (0,1) and t[1]>=0 and t[3]>=0 for t in (lo,hi)),
                    "Malformed or nonfinite numerical record")
            require(value.get("directed") is True and mpf_cmp(tuple(lo),tuple(hi))<=0,
                    "New numerical record must be directed and ordered")
            count+=1
        count+=sum(check_records(v) for v in value.values())
    elif isinstance(value,(list,tuple)):
        count+=sum(check_records(v) for v in value)
    return count

def parent():
    require(sha(BASE/PARENT_NAME)==PARENT_SHA256,"Accepted native-background bytes changed")
    result=load(BASE/PARENT_NAME)
    require(result["checks"]["passed"] is True,"Accepted native-background checks failed")
    return result


def validate_native(packet,accepted):
    require(packet["passed"] is True and packet["parent_report_sha256"]==PARENT_SHA256,
            "Wrong native-control predecessor")
    require(packet["source_family"]==accepted["source_family"]
            and packet["actual_anchor_binding"]["anchored_G_not_Gbar"] is True
            and packet["actual_anchor_binding"]["old_reference_join_family_not_consumed"] is True,
            "Actual anchored family was replaced")
    require(packet["source_function_identities"]["passed"] is True
            and packet["source_function_identities"]["p2_difference"]=="0"
            and packet["source_function_identities"]["logphi_difference"]=="0"
            and packet["slope_proof"]["passed"] is True,
            "Native source formula/slope identity missing")
    for name in ("long_reshape","reference"):
        chart=packet["charts"][name]
        require(chart["passed"] is True and chart["p2_source"]["full_pressure_and_ABC_retained"] is True
                and chart["p2_source"]["physical_J_not_in_denominator"] is True,
                "Native full-pressure controls are incomplete")
        require(set(chart["native_log_absolute_caps"])=={"phi","U","M","A","B","C"},
                "Native field controls incomplete")
        require(chart["limits"]["native_H_K_numeric"] is False
                and chart["limits"]["signed_phi_ABC_values_evaluated"] is False
                and chart["limits"]["seam_proof"] is False,
                "Native bounds were promoted to values or a seam theorem")
    require(packet["limits"]["modified_endpoint_seams_proved"] is False
            and packet["limits"]["full_goal_complete"] is False,
            "New native controls were promoted to global closure")


def validate_loop(packet,native):
    require(packet["passed"] is True and packet["parent_report_sha256"]==PARENT_SHA256,
            "Wrong loop predecessor")
    require(packet["physical_y_derivatives_claimed"] is False,
            "Fixed-phase bounds were promoted to physical-y bounds")
    for name in ("long_reshape","reference"):
        chart=packet["charts"][name]
        require(chart["passed"] is True and chart["logN"]==native["charts"][name]["logN"]
                and chart["log_Pstar_E_upper"]==native["charts"][name]["log_Pstar_E_upper"],
                "Native loop scale changed")
        rows=chart["B_over_Pstar_over_frozen_E0_log_caps"]
        require([v["ordinary_Z_order"] for v in rows]==[0,1,2,3],
                "Loop must retain all four ordinary Z rows")


def validate_finite(packet,accepted):
    from mpmath.libmp import mpf_cmp
    require(packet["passed"] is True and packet["parent_report_sha256"]==PARENT_SHA256,
            "Wrong finite-N predecessor")
    require(packet["all_normalized_I1_errors_below_2_power_minus1024"] is True,
            "Current-N I1 error comparison failed")
    require(packet["limits"]["complete_pressure_global_closure"] is False
            and packet["limits"]["global_five_moment_closure"] is False
            and packet["limits"]["outer_RQy_evaluated"] is False,
            "Finite-N I1 was promoted to complete pressure/global closure")
    require(packet["inherited_lemma_reference"]==accepted[
                "actual_long_reference_background_I1_numeric"]["inherited_lemma_reference"]
            and packet["pressure_identity_reexecuted"] is False,
            "Accepted pressure lemma was replayed/replaced")
    names={"q","q_Z","q_squared","U_q_Z","U_squared","U_U_Z"}
    for name,unit in (("long_reshape","Wlong"),("reference","Wref")):
        chart=packet["charts"][name]
        require(chart["passed"] is True and chart["ordinary_Z_orders"]==[0,1]
                and chart["geometry"]["normalization"]==unit
                and chart["geometry"]["own_weight_used"] is True
                and chart["geometry"]["no_extra_R_weight"] is True,
                "Wrong current-N physical measure")
        radius=chart["signed_I1_radius"]
        require(radius["lower_exact_mpf_tuple"]==radius["upper_exact_mpf_tuple"]==[0,1,-1024,1],
                "Wrong one-time finite-N radius")
        require(set(chart["I1_full"]["six_integrals"])==names
                and chart["I1_full"]["integrated_Q_used"] is False
                and chart["I1_full"]["symmetric_radius_added_once"]==radius
                and chart["B_full"]["symmetric_radius_added_once"]==radius
                and chart["D_full"]["separate_D_error_not_added"] is True
                and chart["S1_full"]["separate_S1_error_not_added"] is True,
                "Finite-N radius repeated or Q roles mixed")
        rows=[row for key in names for row in chart["six_I1_difference_error_log_bounds"][key]]
        rows+=chart["B_difference_error_log_bounds"]
        limit=tuple(chart["I1_error_threshold"]["lower_exact_mpf_tuple"])
        for row in rows:
            if not row.get("exact_zero",False):
                raw=row["log_absolute_upper_bound"]
                require(mpf_cmp(tuple(raw["upper_exact_mpf_tuple"]),limit)<0,
                        "Saved finite-N defect does not fit its claimed radius")
        require(chart["limits"]["seams_proved"] is False
                and chart["limits"]["outer_RQy_evaluated"] is False,
                "Unproved RQy endpoint lost")


def validate_report(report,accepted):
    mutable={"stage","parent_report","parent_report_sha256","input_hashes",
             "implementation_hashes","contract","checks","execution_seconds","execution_audit"}
    for key,value in accepted.items():
        if key not in mutable:
            require(report.get(key)==value,"Accepted mathematical payload changed: "+key)
    require(report["execution_audit"]["inherited_parent_execution_audit"]==accepted["execution_audit"]
            and report["execution_audit"]["current_turn_accepted_ancestors_replayed"] is False,
            "Execution audit was lost/reset")
    require(report["actual_long_reference_pressure_source"]["joined_pressure_source"]["joined_contracted"] is None,
            "Unproved source endpoint was removed")
    counts={key:check_records(report[key]) for key in (NATIVE,LOOP,FINITE)}
    require(all(v>0 for v in counts.values()),"Missing new directed numerical records")
    validate_native(report[NATIVE],accepted)
    validate_loop(report[LOOP],report[NATIVE])
    validate_finite(report[FINITE],accepted)
    return dict(passed=True,new_directed_records=counts,
        accepted_math_payload_preserved=True,actual_native_controls=True,
        full_pressure_numeric=False,new_modified_seams_proved=False)


def run():
    started=time.perf_counter()
    for path,expected in REVIEWED.items():
        require(expected!="REVIEW_PENDING" and sha(ROOT/path)==expected,
                "Independent source review pending or bytes changed: "+path)
    accepted=parent()
    check_pins(accepted)
    sys.path.insert(0,str(ROOT/"src"))
    from openai_ns_reconstruction import original_long_reference_native_controls as native
    from openai_ns_reconstruction import original_long_reference_loop_Z3_bounds as loop
    from openai_ns_reconstruction import original_long_reference_finite_N_I1_numeric as finite
    n=native.build(ROOT,accepted)
    b=loop.build(ROOT,accepted,n)
    f=finite.build(ROOT,accepted,n,b)
    require(accepted==parent(),"New stage changed its accepted input")
    inputs,implementations=merge_pins(accepted,n,b,f)
    report=copy.deepcopy(accepted)
    report.update(stage="actual_native_long_reference_controls_and_current_N_I1_B_D_S1",
        parent_report=PARENT_NAME,parent_report_sha256=PARENT_SHA256,
        input_hashes=inputs,implementation_hashes=implementations)
    report[NATIVE],report[LOOP],report[FINITE]=n,b,f
    report["execution_audit"]=dict(current_turn_accepted_ancestors_replayed=False,
        new_stage_arithmetic_only=True,
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    report["contract"]=dict(actual_native_long_reference_controls=True,
        long_reference_current_N_I1_B_D_S1_C1=True,
        full_pressure_numeric=False,new_modified_endpoint_seams_proved=False,
        global_five_moment_closure=False,full_goal_complete=False)
    report["checks"]=dict(passed=True,packet_invariants=validate_report(report,accepted),
        independent_new_math_review_completed=True,
        independent_runtime_worker_model_metadata_available=False,
        reviewed_implementation_sha256=REVIEWED)
    report["execution_seconds"]=round(time.perf_counter()-started,3)
    json.dumps(report,allow_nan=False)
    return report


def inspect(path):
    report,accepted=load(path),parent()
    require(report["checks"]["passed"] is True
            and report["checks"]["reviewed_implementation_sha256"]==REVIEWED,
            "Saved source reviews changed")
    for path in CODE:
        require(report["implementation_hashes"][path]==sha(ROOT/path),"Saved code changed: "+path)
    return dict(passed=True,pins=check_pins(report),packet_invariants=validate_report(report,accepted),
        serialized_only=True,accepted_builders_or_checkers_replayed=False)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write-report",action="store_true")
    group.add_argument("--inspect-report",action="store_true")
    args=parser.parse_args()
    path=BASE/NAME
    if args.write_report:
        report=run()
        tmp=path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8",newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True,execution_seconds=report["execution_seconds"],bytes=path.stat().st_size,sha256=sha(path))))
    else:
        print(json.dumps(inspect(path),indent=2))
