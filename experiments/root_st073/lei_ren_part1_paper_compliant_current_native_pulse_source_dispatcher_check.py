"""Focused current native pulse ownership and complete chart coverage."""
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import (
    CurrentNativePulseSourceDispatcher, PULSE_CHARTS, CHARTS,
    UNIFORM, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"current_native_pulse_source_dispatcher.json"


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentNativePulseSourceDispatcher(require_checked=False)
        manifest=encode(pack(field.manifest()))
        extra=("whole_current_native_pulse_charts","whole_current_gap_coordinate_overlap")
        if {k:v for k,v in raw.items() if k not in extra}!=manifest:
            raise ValueError("Current native chart sources, graph or coverage changed")
        if (raw["current_downstream_chart_owner_count"]!=20
                or tuple(raw["ordered_current_chart_registry"])!=CHARTS
                or not raw["current_Rp_external_pulse_join_certified"]
                or not raw["all_current_pulse_charts_callable"]
                or not raw["all_current_pulse_charts_installed"]
                or not raw["current_native_pulse_source_ownership_proved"]
                or raw["current_native_pulse_source_ownership_certified"]
                or not all(raw["current_native_pulse_method_identity"].values())
                or not all(raw["retained_original_functional_pulse_identities"].values())
                or raw[UNIFORM] or raw["full_pulse_C4_installed"]):
            raise ValueError("Same native chart ownership and limited interface scope required")
        coverage=raw["current_native_pulse_coverage"]
        if (not coverage["original_six_chart_domain_union_fully_enclosed"]
                or not coverage["exact_gap_boundary_covered_by_main_overlap"]
                or not coverage["actual_exact_source_endpoint_not_replaced"]):
            raise ValueError("Original full gap domain requires coordinate overlap")
        c=MPIntervalContext();c.dps=240
        keys={"y"+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)}
        packets=raw["whole_current_native_pulse_charts"]
        if tuple(packets)!=PULSE_CHARTS:raise ValueError("Native whole chart packet omitted")
        cases={**packets,"pulse_gap_overlap":raw["whole_current_gap_coordinate_overlap"]}
        counts={}
        for chart,packet in cases.items():
            if (field.provider("pulse_gap" if chart=="pulse_gap_overlap" else chart) is not field.native_pulse
                    or packet["actual_five_defect_family_sha256"]!=field.family
                    or packet["implicit_source_sha256"]!=field.source
                    or packet["datum_enclosure_sha256"]!=field.datum_sha
                    or not packet["same_single_current_native_pulse_object_used"]
                    or not packet["current_Rp_external_pulse_join_certified"]
                    or packet["current_native_pulse_source_ownership_certified"]
                    or packet[UNIFORM] or any(packet[k] for k in SCOPES)):
                raise ValueError("Current native packet source/owner/scope differs")
            source=packet["source_packet"];coord=source["coordinate"]
            if chart=="pulse_entrance":key="entrance_t"
            elif chart in ("pulse_gap_end","pulse_end"):key="offset_from_Rv"
            else:key="xi"
            domain=(field.coverage()["supplemental_gap_main_box"] if chart=="pulse_gap_overlap"
                else field.coverage()["chart_boxes"][chart])
            if encode(pack(field.native_pulse.ctx.mpf(domain)))!=coord[key]:
                raise ValueError("Generated packet does not cover declared whole coordinate box: "+chart)
            Zlo,Zhi=endpoints(read_interval(c,source["Z"]))
            if Zlo>-1 or Zhi<1:raise ValueError("Whole axial domain omitted")
            groups=packet["physical_mixed_grids"]
            if set(groups)!={"physical_mixed_derivatives_total_order_le4"}:
                raise ValueError("Ordinary native mixed-four grid required")
            grid=groups["physical_mixed_derivatives_total_order_le4"];count=0
            if len(grid)!=4:raise ValueError("Native velocity/pressure component omitted")
            for rows in grid.values():
                if set(rows)!=keys:raise ValueError("Total-order-four native row omitted")
                for value in rows.values():
                    lo,hi=endpoints(read_interval(c,value))
                    if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                        raise ArithmeticError("Nonfinite current native chart enclosure")
                    count+=1
            if count!=60:raise ValueError("Native source packet row scope differs")
            counts[chart]=count
        # These native methods accept broader intervals; the dispatcher must
        # enforce its advertised subdomains before calling them.
        for chart,value in (("pulse_main",".019"),("pulse_exit","9.9"),("pulse_gap","12.0001")):
            try:field.evaluate(chart,0,value)
            except ValueError:pass
            else:raise ValueError("Native route subdomain guard did not reject "+chart)
        if (any(raw[k] for k in SCOPES)
                or raw["all_profile_source_charts_callable"]
                or raw["full_current_core_to_heat_physical_assembly"]
                or raw["quantitative_flat_velocity_interface_bound_ledger_available"]):
            raise ValueError("Native chart coverage promoted unfinished global scope")
        result=dict(actual_five_defect_family_sha256=field.family,
            implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
            current_input_hashes_checked=len(raw["input_hashes"]),
            current_downstream_chart_owners_checked=20,
            current_native_pulse_chart_names=list(PULSE_CHARTS),
            same_native_pulse_object_for_every_chart_checked=True,
            current_Rp_external_pulse_receipt_consumed=True,
            retained_original_functional_pulse_identities_checked=True,
            whole_native_source_rows_checked=counts,total_whole_native_source_rows_checked=sum(counts.values()),
            original_six_chart_domain_union_fully_enclosed=True,
            reciprocal_boundary_covered_by_legal_coordinate_overlap=True,
            actual_source_endpoint_not_replaced=True,
            narrower_native_route_guards_checked=True,
            current_native_pulse_source_ownership_certified=True,
            all_current_pulse_charts_installed=True,all_current_pulse_charts_callable=True,
            **{UNIFORM:False},full_pulse_C4_installed=False,
            quantitative_flat_velocity_interface_bound_ledger_available=False,
            full_current_core_to_heat_physical_assembly=False,
            **dict.fromkeys(SCOPES,False),all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),
                Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current native pulse chain: 20 owners, six same-object charts, 420 whole-domain rows",flush=True)
    return result


if __name__=="__main__":
    run()
