"""Write or inspect the compact source-bound leading-viscous adapter report."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT/"experiments/root_st073"
sys.path.insert(0,str(ROOT/"src"))

PREFIX = "lei_ren_part1_paper_"
NAME = PREFIX+"compliant_current_original_leading_viscous_source_numeric.json"
PARENT = PREFIX+"compliant_current_original_macro_endpoint_m15_numeric.json"
PARENT_SHA256 = "16eb888da56fad8e5696b5fd8400b827ff7634e2c8f524110267cfcc6f525b80"
FIXTURE = "scripts/check_original_leading_viscous_m15.py"
RUNNER = "scripts/run_original_leading_viscous_m15.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()


def require(condition,message):
    if not condition:
        raise AssertionError(message)


def load_parent():
    path = DIRECTORY/PARENT
    require(sha(path)==PARENT_SHA256,"Accepted parent bytes changed")
    packet = json.loads(path.read_bytes())
    require(packet["checks"]["passed"] is True,"Parent is not passed")
    return packet


def check_pins(report):
    count = 0
    for base,key in ((DIRECTORY,"input_hashes"),(ROOT,"implementation_hashes")):
        require(isinstance(report.get(key),dict) and report[key],"Missing "+key)
        for name,expected in report[key].items():
            require(sha(base/name)==expected,"Changed pin "+name)
            count += 1
    return count


def check_packet(report,parent):
    require(report.get("parent_report")==PARENT and report.get("parent_report_sha256")==PARENT_SHA256,
            "Wrong source parent")
    require(report.get("source_family")==parent["source_family"],"Source family changed")
    inherited = report["inherited_m15_packet"]
    require(inherited["passed"] is True and inherited["report_sha256"]==PARENT_SHA256,
            "Inherited m15 packet reference missing")
    require(inherited["ledger_hydrated_by_reference"] is True,"Ledger representation changed")
    samples = inherited["samples"]
    require(len(samples)==len(parent["samples"])==3,"Three inherited axial packets required")
    for row,original in zip(samples,parent["samples"]):
        require(row["Z"]==original["Z"] and row["raw_count"]==213 and row["reduced_count"]==183,
                "Inherited ledger count or axial query changed")
        for key,value in (("raw_sha256",original["raw_m15_C0_C1_terms"]),
                          ("reduced_sha256",original["correlated_reduced_m15_C0_C1_terms"]),
                          ("signed_inner_lower_six_sha256",original["raw_m15_C0_C1_terms"][:6]),
                          ("independent_target_sha256",original["independent_target"])):
            require(row[key]==digest(value),"Changed inherited "+key)
        require(row["reserved_count"]==1 and row["terms_replaced_or_removed"]==0,
                "Source adapter changed the active ledger")
    summaries = {}
    for name in ("macro_source","switch_source"):
        packet = report[name]
        require(packet.get("source_bindings_passed") is True
                and packet.get("endpoint_interface_passed") is True,
                name+" advertised endpoint interface failed")
        require(all(row.get("passed") is True for row in packet["family_checks"].values()),
                name+" source family failed")
        require(isinstance(packet.get("source_bindings"),(dict,list)),name+" source binding missing")
        rows = packet.get("point_remainder_rows")
        require(isinstance(rows,list) and rows,name+" new coefficient rows missing")
        require(isinstance(packet.get("blocked_rows"),list),name+" missing-input records required")
        summaries[name] = {"point_remainder_row_count":len(rows),
                           "blocked_row_count":len(packet["blocked_rows"])}
    ftc = report["independent_leading_FTC"]
    require(ftc.get("passed") is True and ftc.get("old_error_caps_cancelled_or_removed") is False,
            "Independent leading FTC or error preservation failed")
    require(ftc.get("physical_R_once") is True and ftc.get("physical_Pstar_once") is True
            and ftc.get("ordinary_Z1") is True and ftc.get("first1_H_not_zeroed") is True,
            "Independent leading source scale/history failed")
    require(len(ftc["endpoint_recipes"])==4 and all(
        recipe.get("modified_Qy_alias") is False and recipe.get("recipe_is_archived_graph_node") is False
        for recipe in ftc["endpoint_recipes"].values()),"Leading pressure port namespace failed")
    require(report["numeric_integral_scope_complete"] is False,"Unsupported numeric integral closure")
    switch = report["switch_source"]
    require(switch.get("integral_C0_interface_passed") is True
            and switch["limits"]["integral_C0_rows_compared"]==3
            and switch["limits"]["integral_C0_rows_comparison_passed"]==3,
            "Phase-one signed source integral C0 comparisons failed")
    require(all(value is True for value in report["contract"].values()),"Adapter contract changed")
    require(all(value is False for value in report["limits"].values()),"A global completion claim was introduced")
    return {"passed":True,"inherited_raw213_reduced183_exact":True,
            "target_six_signed_and_reserved_once_retained":True,"adapters":summaries,
            "global_completion_claim":False}


def run():
    started = time.perf_counter()
    from check_original_leading_viscous_m15 import run_fixture
    from openai_ns_reconstruction.original_leading_viscous_m15 import OriginalLeadingViscousM15
    fixture = run_fixture()
    require(fixture["passed"] is True,"Independent operator fixture failed")
    owner = OriginalLeadingViscousM15()
    packet = owner.build()
    code = dict(owner.implementation_hashes)
    for relative in (FIXTURE,RUNNER):
        code[relative] = sha(ROOT/relative)
    report = dict(packet,input_hashes=owner.input_hashes,implementation_hashes=code)
    checks = check_packet(report,owner.parent)
    report["checks"] = {"passed":True,"packet_invariants":checks,"independent_fixture":fixture,
                       "source_bound_new_operators_only":True,"no_ancestor_replay":True}
    report["execution_seconds"] = round(time.perf_counter()-started,3)
    return report


def inspect(path):
    report = json.loads(path.read_bytes())
    require(report["checks"]["passed"] is True,"Saved report is not passed")
    pins = check_pins(report)
    checks = check_packet(report,load_parent())
    fixture = report["checks"]["independent_fixture"]
    require(fixture["passed"] is True and all(row["passed"] is True for row in fixture["controls"]),
            "Saved fixture evidence failed")
    return {"passed":True,"pins":pins,"packet_invariants":checks,
            "serialized_only":True,"owner_fixture_or_ancestor_run":False}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-report",action="store_true")
    modes.add_argument("--inspect-report",action="store_true")
    args = parser.parse_args()
    path = DIRECTORY/NAME
    if args.write_report:
        result = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8",newline="\n")
        temporary.replace(path)
    else:
        print(json.dumps(inspect(path),indent=2))
