"""Produce/inspect new signed source C1 without replaying saved mass or ancestors."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
DIRECTORY=ROOT/"experiments/root_st073"
sys.path.insert(0,str(ROOT/"src"))
NAME="lei_ren_part1_paper_compliant_current_original_signed_C1_source_numeric.json"
PARENT="lei_ren_part1_paper_compliant_current_original_source_jet_mass_numeric.json"
PARENT_SHA256="c3b2e2a974ebca402e2d2762ab2c508150fdd72f80aa2882152e944bdad84efe"
CODE=("src/openai_ns_reconstruction/original_signed_C1_source.py",
      "src/openai_ns_reconstruction/original_modal_source_point.py",
      "src/openai_ns_reconstruction/original_factored_source_jet.py",
      "scripts/check_original_signed_C1_source.py","scripts/run_original_signed_C1_source.py")


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()


def require(condition,message):
    if not condition:raise AssertionError(message)


def load_parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,"Accepted mass/source parent bytes changed")
    result=json.loads((DIRECTORY/PARENT).read_bytes())
    require(result["checks"]["passed"] is True,"Parent not passed")
    return result


def check_packet(report,parent):
    require(report["parent_report"]==PARENT and report["parent_report_sha256"]==PARENT_SHA256,"Wrong source parent")
    require(report["source_family"]==parent["source_family"],"Source family changed")
    require(report["inherited_m15_packet"]==parent["inherited_m15_packet"] and
            report["inherited_m15_packet_digest"]==parent["inherited_m15_packet_digest"],"m15 ledger witnesses changed")
    require(report["inherited_mass_packet_digest"]==digest(parent["interior_mass"]),"Saved mass changed")
    source=report["source_jet"]
    require(source["passed"] is True and len(source["rows"])==3,"New same-generator source interface failed")
    require(source["source_bindings"]["passed"] is True and all(
        row["passed"] is True for row in source["source_bindings"]["AST"].values()),"Specific source AST proof failed")
    require(all(row["passed"] is True and row["bound_to_accepted_parent"] is True
                for row in source["source_pins"].values()),"New source input not bound to accepted parent")
    require(all(row["passed"] is True and row["shared_Phi0_cancellation_is_algebraic"] is True
                for row in source["comparison_q"]["correlated_shared_phi_proof"]),"Shared-Phi source correlation failed")
    require(isinstance(source["source_scope"],str) and isinstance(source["full_source_scope"],bool),"Explicit source scope required")
    for row in source["rows"]:
        require(row["ordinary_Z_orders"]==[0,1,2] and len(row["source_rows"])==3,"Missing ordinary B_ZZ")
        require(all(isinstance(value.get("approximation"),str)
                    and isinstance(value.get("source_enclosure_exact_mpf_tuples"),list)
                    and value.get("approximation_is_exact_source_value") is False
                    for value in row["source_rows"]),"Signed source approximation/enclosure separation failed")
    integral=report["signed_source_integral"]
    require(integral["passed"] is True and integral["row_count"]==19 and len(integral["rows"])==19,"New source C1 interface failed")
    for row,original in zip(integral["rows"],parent["signed_source_integral_C0"]["rows"]):
        require(row["Z"]==original["Z"] and row["phase"]==original["phase"]
                and row["physical_scale"]==original["physical_scale"],"Source phase or physical scale changed")
        require(row["source_scope"]==source["source_scope"] and row["ordinary_Z_orders"]==[0,1]
                and len(row["main_component_C0_C1"])==2,"Source C1 scope/order changed")
        require(row["mass_Z_independent"] is True and row["weighted_raw47_remainder_retained"] is True
                and row["full_weighted_primitive_evaluated"] is False,"Weighted primitive promoted or raw47 removed")
    guard=report["rejected_draft_source_origin"]
    require(guard["passed"] is False and not guard["rows"] and guard["source_origin_guard"]["passed"] is False,
            "Mismatched R100-origin draft was accepted")
    require(guard["contract"]["ordinary_B_ZZ_recovered"] is False,"Invalid draft B_ZZ claim survived")
    origin=guard["source_origin_guard"]
    require(origin["blocked_before_compute_source"] is True and origin["family_triple_check"]["passed"] is True
            and origin["source_input_bindings"]["passed"] is True
            and all(row["passed"] is True for row in origin["ast_bindings"].values()),"Rejected draft provenance guard failed")
    require(report["contract"]["accepted_ancestor_replayed"] is False
            and report["contract"]["saved_mass_recomputed"] is False
            and report["contract"]["raw213_reduced183_and_target_retained"] is True,"Preservation contract changed")
    require(all(value is False for value in report["limits"].values()),"Global completion claim introduced")
    return {"passed":True,"new_source_jet_queries":3,"raw46_component_C0_C1_queries":19,
            "source_scope":source["source_scope"],"full_source_scope":source["full_source_scope"],
            "invalid_R100_macro_alias_rejected":True,"raw213_reduced183_and_mass_exact_witnesses_retained":True,
            "full_weighted_source_or_global_five_closure":False}


def check_pins(report):
    count=0
    for base,key in ((DIRECTORY,"input_hashes"),(ROOT,"implementation_hashes")):
        require(bool(report.get(key)),"Missing source pins")
        for name,expected in report[key].items():
            require(sha(base/name)==expected,"Changed pin: "+name);count+=1
    return count


def run():
    started=time.perf_counter()
    from check_original_signed_C1_source import run_fixture
    from openai_ns_reconstruction.original_signed_C1_source import OriginalSignedC1Source
    fixture=run_fixture()
    owner=OriginalSignedC1Source();packet=owner.build()
    code=dict(owner.implementation_hashes);code.update({name:sha(ROOT/name) for name in CODE})
    report=dict(packet,input_hashes=owner.input_hashes,implementation_hashes=code)
    checks=check_packet(report,owner.parent)
    report["checks"]={"passed":True,"packet_invariants":checks,"independent_fixture":fixture,
                      "source_bound_new_operators_only":True,"no_accepted_ancestor_replay":True}
    report["execution_seconds"]=round(time.perf_counter()-started,3)
    return report


def inspect(path):
    report=json.loads(path.read_bytes());require(report["checks"]["passed"] is True,"Saved report not passed")
    pins=check_pins(report);checks=check_packet(report,load_parent())
    fixture=report["checks"]["independent_fixture"]
    require(fixture["passed"] is True and all(row["passed"] is True for row in fixture["controls"]),"Saved fixture failed")
    return {"passed":True,"pins":pins,"packet_invariants":checks,"serialized_only":True,
            "owner_fixture_or_ancestor_run":False}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);modes=parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-report",action="store_true");modes.add_argument("--inspect-report",action="store_true")
    args=parser.parse_args();path=DIRECTORY/NAME
    if args.write_report:
        result=run();temporary=path.with_suffix('.json.tmp')
        temporary.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8',newline='\n')
        temporary.replace(path)
    else:print(json.dumps(inspect(path),indent=2))
