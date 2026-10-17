"""Produce or inspect only the new source-jet / bounded mass stage."""
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
NAME = "lei_ren_part1_paper_compliant_current_original_source_jet_mass_numeric.json"
PARENT = "lei_ren_part1_paper_compliant_current_original_leading_viscous_source_numeric.json"
PARENT_SHA256 = "d1aa3644ea8158fe78dccb210a0a75a4c81ad257b8504758cb5797845b8003d4"
CODE = (
    "src/openai_ns_reconstruction/original_source_jet_mass.py",
    "src/openai_ns_reconstruction/original_source_jet_kernel.py",
    "src/openai_ns_reconstruction/original_interior_mass_source.py",
    "src/openai_ns_reconstruction/original_lower_dy_source.py",
    "scripts/check_original_source_jet_kernel.py",
    "scripts/run_original_source_jet_mass.py",
)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=False).encode()).hexdigest()


def require(condition,message):
    if not condition:
        raise AssertionError(message)


def load_parent():
    require(sha(DIRECTORY/PARENT)==PARENT_SHA256,"Accepted leading source report bytes changed")
    result = json.loads((DIRECTORY/PARENT).read_bytes())
    require(result["checks"]["passed"] is True,"Accepted source parent not passed")
    return result


def check_packet(report,parent):
    require(report["parent_report"]==PARENT and report["parent_report_sha256"]==PARENT_SHA256,"Wrong source parent")
    require(report["source_family"]==parent["source_family"],"Source family changed")
    inherited = report["inherited_m15_packet"]
    require(inherited==parent["inherited_m15_packet"] and report["inherited_m15_packet_digest"]==digest(inherited),
            "Accepted m15 ledger witness changed")
    require(all(row["raw_count"]==213 and row["reduced_count"]==183 and row["reserved_count"]==1
                and row["terms_replaced_or_removed"]==0 for row in inherited["samples"]),"Active m15 ledger changed")
    mass = report["interior_mass"]
    require(mass["passed"] is True and len(mass["rows"])==7,"New bounded mass interface failed")
    require(mass["source"]["signed_source_integral"] is True,"Mass source identity changed")
    for row in mass["rows"]:
        require(row["mass_Z_independent"] is True,"Axial-independent pulse mass changed")
        require(isinstance(row.get("source_enclosure_exact_mpf_tuples"),list),"Exact mass enclosure missing")
        require(row["convergence"]["passed"] is True,"128/256 mass consistency failed")
        if row["phase"] not in ("0","1"):
            require(row["approximation_is_exact_source_value"] is False and row["source_value_not_selected"] is True,
                    "Interior mass approximation was promoted to a function value")
    c0 = report["signed_source_integral_C0"]
    require(c0["passed"] is True and c0["row_count"]==19,"Signed raw46 main-component C0 interface failed")
    require(c0["signed_full_primitive_evaluated"] is False and c0["old_error_caps_removed"] is False,
            "Main component was promoted to full weighted primitive")
    original = {(row["Z"],row["coordinate"]["value"]):row for row in parent["switch_source"]["point_remainder_rows"]
                if row.get("name")=="first_switch_viscous_source_integral"}
    from decimal import Decimal
    for row in c0["rows"]:
        match = next(v for (z,s),v in original.items() if z==row["Z"] and Decimal(s)==Decimal(row["phase"]))
        require(row["physical_scale"]==match["physical_scale"],"Raw46 physical scale changed")
        require(row["main_component_C0"]["approximation_is_exact_source_value"] is False,"Raw46 approximation promoted")
        require(row["raw47_weighted_remainder_retained"] is True,"Raw47 removed")
    require(report["contract"]["accepted_ancestor_replayed"] is False,"Ancestor replay introduced")
    require(report["contract"]["positive_caps_selected_as_values"] is False,"Positive cap selected as source value")
    require(report["contract"]["interval_midpoints_claimed_as_actual_values"] is False,"Midpoint promoted")
    lower = report["transported_lower_dy"]
    require(lower["passed"] is True and len(lower["rows"])==3 and not lower["blocked_rows"],"Transported lower source failed")
    require(all(row["ordinary_Z_orders"]==[0,1,2] and len(row["result"])==3
                and row["transport"]["macro_J66987_applied"] is False for row in lower["rows"]),
            "Transported source convention or geometry changed")
    require(lower["source_bindings"]["coefficient_convention"]["passed"] is True and all(
        row["coefficient_convention"]["factorial_conversion_applied"] is True
        and row["actual"]["inherited_conservative_envelope"] is True
        and row["comparison"]["positive_caps_only_widen_enclosure"] is True for row in lower["rows"]),
        "Taylor convention or conservative source envelope changed")
    require(all(set(row["family"])==set(parent["source_family"]) for row in lower["family_checks"].values()),
            "Source-family witness must not duplicate a source report")
    require(report["lower_local_operators"]["passed"] is True and
            report["lower_local_operators"]["raw24_positive_bound_replaced_or_removed"] is False,
            "Lower source operators replaced a positive cap")
    require(all(value is False for value in report["limits"].values()),"Global completion claim introduced")
    return {"passed":True,"raw213_reduced183_exact_witnesses_retained":True,"mass_phase_count":7,
            "raw46_C0_main_component_count":19,"transported_lower_Z0_Z1_Z2_count":3,
            "old_source_errors_and_target_retained":True,
            "source_function_equality_from_interval_overlap":False,"global_completion_claim":False}


def check_pins(report):
    count = 0
    for base,key in ((DIRECTORY,"input_hashes"),(ROOT,"implementation_hashes")):
        require(bool(report.get(key)),"Missing source pins")
        for name,expected in report[key].items():
            require(sha(base/name)==expected,"Changed pin: "+name)
            count += 1
    return count


def run():
    started = time.perf_counter()
    from check_original_source_jet_kernel import run_fixture
    from openai_ns_reconstruction.original_source_jet_mass import OriginalSourceJetMass
    fixture = run_fixture()
    owner = OriginalSourceJetMass()
    packet = owner.build()
    code = dict(owner.implementation_hashes)
    code.update({relative:sha(ROOT/relative) for relative in CODE})
    report = dict(packet,input_hashes=owner.input_hashes,implementation_hashes=code)
    checks = check_packet(report,owner.parent)
    report["checks"] = {"passed":True,"packet_invariants":checks,"independent_fixture":fixture,
                       "source_bound_new_operators_only":True,"no_accepted_ancestor_replay":True}
    report["execution_seconds"] = round(time.perf_counter()-started,3)
    return report


def inspect(path):
    report = json.loads(path.read_bytes())
    require(report["checks"]["passed"] is True,"Saved report not passed")
    count = check_pins(report)
    checks = check_packet(report,load_parent())
    fixture = report["checks"]["independent_fixture"]
    require(fixture["passed"] is True and all(row["passed"] is True for row in fixture["controls"]),"Saved fixture evidence failed")
    return {"passed":True,"pins":count,"packet_invariants":checks,"serialized_only":True,
            "owner_fixture_or_ancestor_run":False}


if __name__=="__main__":
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
