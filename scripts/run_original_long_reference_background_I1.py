"""Produce or inspect the native long/reference background pressure sector.

Only the two new stage builders run. Accepted data and graph patches are read
without replaying accepted arithmetic, constructors, norms or quadrature.
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
PARENT_NAME="lei_ren_part1_paper_compliant_current_original_switch_post_finite_N_I1.json"
PARENT_SHA256="9d0d957bae129853030bad4cd5c484eb27e8b13b0623601e458fc4e95e93bc2b"
NAME="lei_ren_part1_paper_compliant_current_original_long_reference_background_I1.json"
BG_KEY="actual_long_reference_background_I1_numeric"
SOURCE_KEY="actual_long_reference_pressure_source"
CHARTS=("long_reshape","reference")
CODE=("src/openai_ns_reconstruction/original_long_reference_background_I1_numeric.py",
      "src/openai_ns_reconstruction/original_long_reference_pressure_source.py",
      "scripts/run_original_long_reference_background_I1.py")
REVIEWED={CODE[0]:"ff73a7ccdadf1f9cbeec200700284ceaacf5fbdd26f371c7b7ae1162b9c2f5b9",
          CODE[1]:"7bd5e325dcb6534bfb8ac7d880948f13d52777f3b8397bda032f74385804cc37"}

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

def parent():
    path=BASE/PARENT_NAME
    require(sha(path)==PARENT_SHA256,"Accepted switch/post parent bytes changed")
    result=load(path)
    require(result["checks"]["passed"] is True,"Accepted parent checks failed")
    return result

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

def orders(rows,label,count):
    require(isinstance(rows,list) and len(rows)==count
            and [r.get("ordinary_Z_order") for r in rows]==list(range(count)),
            "Incomplete ordinary Z rows: "+label)

def same_values(left,right):
    return len(left)==len(right) and all(
        a["value"][side]==b["value"][side] for a,b in zip(left,right)
        for side in ("lower_exact_mpf_tuple","upper_exact_mpf_tuple"))

def validate_background(packet,accepted):
    from mpmath.libmp import mpf_cmp
    require(packet["passed"] is True and packet["parent_report_sha256"]==PARENT_SHA256
            and set(packet["charts"])==set(CHARTS),"Wrong background source parent/charts")
    contract=packet["source_contract"]
    require(contract["M_not_reset_to_U"] is True
            and contract["reference_incoming_is_long_endpoint"] is True
            and contract["C2_cancelled_before_interval_boxing"] is True
            and contract["U_squared_coefficient"]=="2*h*Z^2/Lz^2",
            "Native mean memory or pressure algebra changed")
    inherited=packet["inherited_lemma_reference"]
    require(inherited["passed"] is True and inherited["symbolic_remainder"]=="0"
            and inherited["Q4_source_pin_verified"] is True
            and inherited["reexecuted"] is False,"B+D pressure lemma changed")
    radius=packet["tail_control"]["r"]
    require(radius["lower_exact_mpf_tuple"]==[0,1,-1024,1]
            and radius["upper_exact_mpf_tuple"]==[0,1,-1024,1]
            and packet["tail_control"]["lower_zero_is_an_interval_enclosure_not_a_selected_value"] is True,
            "Tail radius was selected as a field value")
    logarithmic_cap=packet["tail_control"]["log_r"]["lower_exact_mpf_tuple"]
    for name in CHARTS:
        chart=packet["charts"][name]
        require(chart["passed"] is True and chart["normalization"]==
                ("Wlong" if name=="long_reshape" else "Wref"),"Wrong endpoint physical scale")
        for key in ("U_constant_rows","incoming_M","endpoint_M","whole_segment_M_hull"):
            orders(chart[key],name+" "+key,5)
        proof=chart["tail_proof"]
        require(all(proof["inequalities"].values()) and proof["tails_not_set_to_zero"] is True
                and proof["exp_minus_L_nonzero"] is True
                and mpf_cmp(tuple(proof["log_L_minus_2L"]["upper_exact_mpf_tuple"]),
                            tuple(logarithmic_cap))<0,"Unproved normalized tail cap")
        pointwise=proof["pointwise_decay_interval"]
        endpoint=proof["endpoint_decay_interval"]
        require(pointwise["lower_exact_mpf_tuple"]==[0,0,0,0]
                and pointwise["upper_exact_mpf_tuple"]==[0,1,0,1]
                and endpoint["lower_exact_mpf_tuple"]==[0,0,0,0]
                and endpoint["upper_exact_mpf_tuple"]==[0,1,-1024,1],
                "Whole-segment and endpoint decay enclosures were confused")
        a=chart["arithmetic"]
        for key in ("Q_pointwise","Q_endpoint","integrated_q","alpha","beta"):
            orders(a[key],name+" "+key,3)
        require(a["source"]["Q_and_integrated_Iq_distinct"] is True
                and a["source"]["pointwise_decay_interval"]=="[0,1]"
                and a["source"]["endpoint_decay_interval"]=="[0,r]"
                and "two_h_Z2_over_L2" in a["coefficients"],"Pointwise Q or D coefficient changed")
        require(set(a["six_integrals"])==
                {"q","q_Z","q_squared","U_q_Z","U_squared","U_U_Z"},"Incomplete six I1")
        for key,rows in a["six_integrals"].items():
            orders(rows,name+" "+key,2)
        for key in ("D","B_over_W","S1_over_W"):
            orders(a[key],name+" "+key,2)
        require(a["algebra"]["normalized_formula"]=="B/W=2*C0*J0+C1*J1"
                and a["algebra"]["C2_cancelled_before_interval_boxing"] is True
                and a["algebra"]["S1_formula"]=="S1=B+D", "Boundary correlation changed")
    old=accepted["actual_switch_post_background_I1_numeric"]["charts"]["post_switch"]["native_fields"]
    long,reference=(packet["charts"][key] for key in CHARTS)
    require(same_values(long["U_constant_rows"],old["U"])
            and same_values(reference["U_constant_rows"],old["U"])
            and same_values(long["incoming_M"],old["M"])
            and same_values(reference["incoming_M"],long["endpoint_M"]),"Source inlet was replaced/reset")
    joined=packet["joined_long_reference_background"]
    require(joined["normalization"]=="Wref=Rz^2"
            and joined["source_proof"]["ratio_positive_nonzero"] is True
            and joined["source_proof"]["no_selected_zero"] is True
            and joined["source_proof"]["no_independent_radius_added"] is True,
            "Joined background scale was selected or counted twice")
    require(packet["radii"]["radius_not_materialized"] is True
            and packet["limits"]["current_N_not_enclosed"] is True
            and packet["limits"]["finite_N_loop_not_evaluated"] is True
            and packet["limits"]["full_pressure_global_closure"] is False,
            "Background enclosures were promoted to finite-N closure")
    return dict(passed=True,ordinary_C1=True,native_mean_memory_retained=True,
                logarithmic_endpoint_normalizations=True,nonzero_tail_enclosures=True,
                background_only=True,current_N_numeric=False)

def validate_source(packet,accepted):
    require(packet["passed"] is True and packet["parent_report_sha256"]==PARENT_SHA256,
            "Wrong pressure-source parent")
    old=accepted["actual_switch_post_pressure_source"]["graph_patch"]
    graph=packet["graph_patch"]
    require(graph["base_node_count"]==graph["new_node_start"]==185308
            and graph["base_prefix_unchanged"] is True
            and graph["base_prefix_canonical_nodes_sha256"]==old["full_canonical_nodes_sha256"]
            and graph["appended_node_count"]==len(graph["appended_nodes"])
            and graph["full_node_count"]==185308+len(graph["appended_nodes"]),
            "Accepted graph prefix changed")
    limits=packet["limits"]
    require(limits["outer_RQy_retained"] is True and limits["all_caps_retained"] is True
            and limits["outer_RQy_evaluated"] is False
            and limits["complete_pressure_global_closure"] is False
            and limits["background_finite_N_corrected"] is False,
            "Source extension was promoted to a numerical pressure solution")
    contract=packet["source_contract"]
    require(contract["current_N_source_only"] is True and contract["common_N_handle"]==18
            and contract["local_S1_added_once"] is True
            and contract["cumulative_S1_not_added"] is True
            and contract["nonreset_histories"] is True,
            "Current-N source, local integral or memory scope changed")
    joined=packet["joined_pressure_source"]
    require(isinstance(joined["joined_uncontracted"],list)
            and joined["outer_RQy_and_caps_retained"] is True,
            "Uncontracted endpoint evidence was lost")
    if limits["seams_proved"]:
        proof=packet["endpoint_identity_receipts"]
        require(proof["passed"] is True and proof["modified_Uy_Qy_closed"] is True
                and proof["native_derivative_naturality"]["passed"] is True,
                "New modified endpoint contraction lacks a source proof")
        for name in ("post_power_upper_to_long_lower","long_upper_to_reference_lower"):
            require(proof["seams"][name]["finite_N_Qy_identity_proven"] is True,
                    "Unproved endpoint contracted")
        linear=joined["linearization_receipt"]["long_to_reference"]
        require(linear["passed"] is True and linear["scope"]=="sum_and_negative_rows_only"
                and linear["targets"]==[160287,160288,166149,166150]
                and linear["removed_target_coefficients_zero"] is True
                and linear["surviving_target_dependencies_empty"] is True,
                "Internal endpoint targets survived signed normalization")
        orders(linear["rows"],"signed normalization",2)
        for k,row in enumerate(linear["rows"]):
            require(all(value==0 for value in row["target_coefficients"].values())
                    and row["reduced_handle"]==joined["joined_contracted"][k],
                    "Endpoint integer coefficients did not cancel")
    else:
        proof=packet["endpoint_identity_receipts"]
        require(joined["joined_contracted"] is None
                and joined["contracted_after_long_raw"] is None
                and joined["contracted_after_long"] is None
                and joined["contracted_before_linearization"] is None
                and joined["contraction_allowed"] is False
                and limits["contracted_endpoint_atoms_removed"] is False
                and proof["passed"] is False
                and proof["modified_Uy_Qy_closed"] is False
                and proof["natural_B_y_closed"] is False,
                "Unproved source endpoint was removed or advertised as proved")
        for name in ("post_power_upper_to_long_lower","long_upper_to_reference_lower"):
            require(proof["seams"][name]["finite_N_Qy_identity_proven"] is False
                    and proof["seams"][name]["contraction_allowed"] is False,
                    "Unaccepted modified endpoint theorem was retained")
    require(limits["old_212_terms_and_errors_unchanged"] is True
            and limits["older_germ_unproved_status_retained"] is True
            and limits["older_germ_used_as_proof"] is False,
            "Old source terms or unproved status were lost")
    return dict(passed=True,source_prefix_unchanged=True,
                current_N_source_only=True,new_internal_seams_proved=limits["seams_proved"],
                numerical_pressure_complete=False)

def validate_report(report,accepted):
    mutable={"stage","parent_report","parent_report_sha256","input_hashes",
             "implementation_hashes","contract","checks","execution_seconds","execution_audit"}
    for key,value in accepted.items():
        if key not in mutable:
            require(report.get(key)==value,"Accepted mathematical payload changed: "+key)
    require(report["parent_report"]==PARENT_NAME and report["parent_report_sha256"]==PARENT_SHA256,
            "Wrong accepted predecessor")
    require(report["execution_audit"]["inherited_parent_execution_audit"]==accepted["execution_audit"]
            and report["execution_audit"]["current_turn_accepted_ancestors_replayed"] is False,
            "Execution audit was reset")
    counts={key:check_records(report[key]) for key in (BG_KEY,SOURCE_KEY)}
    require(counts[BG_KEY]>0,"No new directed background enclosures")
    return dict(background=validate_background(report[BG_KEY],accepted),
                pressure_source=validate_source(report[SOURCE_KEY],accepted),
                new_directed_records=counts,accepted_math_payload_preserved=True)

def run():
    started=time.perf_counter()
    for path,expected in REVIEWED.items():
        require(expected!="REVIEW_PENDING" and sha(ROOT/path)==expected,
                "Independent source review pending or bytes changed: "+path)
    accepted=parent()
    check_pins(accepted)
    sys.path.insert(0,str(ROOT/"src"))
    from openai_ns_reconstruction import original_long_reference_background_I1_numeric as bg
    from openai_ns_reconstruction import original_long_reference_pressure_source as source
    background=bg.build(ROOT,accepted)
    pressure=source.build(ROOT,accepted,background)
    require(accepted==parent(),"New stage changed its accepted input")
    inputs,implementations=merge_pins(accepted,background,pressure)
    report=copy.deepcopy(accepted)
    report.update(stage="native_long_reference_background_I1_B_D_S1_and_current_N_pressure_source",
        parent_report=PARENT_NAME,parent_report_sha256=PARENT_SHA256,
        input_hashes=inputs,implementation_hashes=implementations)
    report[BG_KEY],report[SOURCE_KEY]=background,pressure
    report["execution_audit"]=dict(current_turn_accepted_ancestors_replayed=False,
        new_stage_arithmetic_only=True,
        inherited_parent_execution_audit=copy.deepcopy(accepted["execution_audit"]))
    report["contract"]=dict(long_reference_background_pressure_sector_C1=True,
        same_current_N_pressure_source_extension=True,
        long_reference_finite_N_numeric=False,complete_pressure_numeric=False,
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
        require(report["implementation_hashes"][path]==sha(ROOT/path),
                "Saved implementation bytes changed: "+path)
    return dict(passed=True,pins=check_pins(report),
        packet_invariants=validate_report(report,accepted),serialized_only=True,
        accepted_builders_or_checkers_replayed=False)

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
        tmp.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",
                       encoding="utf-8",newline="\n")
        tmp.replace(path)
        print(json.dumps(dict(passed=True,execution_seconds=report["execution_seconds"],
                              bytes=path.stat().st_size,sha256=sha(path))))
    else:
        print(json.dumps(inspect(path),indent=2))
