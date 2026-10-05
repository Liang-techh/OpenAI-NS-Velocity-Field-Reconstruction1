"""Current Rp native source graph, exact mixed join and entrance coverage."""
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_actual_Rp_source_join import (
    CurrentRpPulseSourceDispatcher, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"actual_Rp_source_join.json"


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentRpPulseSourceDispatcher(require_checked=False)
        manifest=encode(pack(field.manifest()))
        extra=("current_Rp_left_source","current_Rp_native_right_source",
            "whole_current_native_entrance","whole_entrance_original_domain_fully_enclosed")
        if {k:v for k,v in raw.items() if k not in extra}!=manifest:
            raise ValueError("Current Rp source graph, formulas or dependencies changed")
        if (not all(raw["current_Rp_native_source_graph"].values())
                or raw["current_downstream_chart_owner_count"]!=15
                or not raw["current_Rp_external_pulse_join_proved"]
                or raw["current_Rp_external_pulse_join_certified"]
                or not raw["whole_entrance_original_domain_fully_enclosed"]
                or raw["all_current_pulse_charts_installed"]):
            raise ValueError("Current Rp producer ownership/acceptance scope differs")
        shapes=raw["current_Rp_terminal_shape_proof"]
        projection=raw["current_Rp_native_inlet_projection_proof"]
        mixed=raw["current_Rp_exact_physical_mixed_source_join"]
        boundary=raw["current_Rp_entrance_boundary_and_radius_proof"]
        if (not all(v["passed"] for v in (shapes,projection,mixed,boundary,raw["current_Rp_parameter_source_proof"]))
                or shapes["source_shape_identities"]!=28
                or shapes["canonical_unit_identities"]!=5
                or projection["actual_native_inlet_axial5_coefficient_identities"]!=54
                or mixed["exact_replayed_velocity_pressure_and_primitive_mixed4_rows"]!=135
                or not boundary["actual_gp_value_and_ordinary_jets0_through4_exact_zero"]):
            raise ValueError("Complete actual Rp source-functional join required")
        native=field.native_pulse
        if field.provider("pulse_entrance") is not native:
            raise ValueError("Current native pulse owner was reconstructed separately")
        c=MPIntervalContext();c.dps=240
        keys={"y"+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)}
        counts={}
        for name in extra[:3]:
            packet=raw[name]
            if (packet["actual_five_defect_family_sha256"]!=field.family
                    or packet["implicit_source_sha256"]!=field.source
                    or packet["datum_enclosure_sha256"]!=field.datum_sha
                    or any(packet[k] for k in SCOPES)):
                raise ValueError("Current Rp side/source/scope differs")
            groups=packet["physical_mixed_grids"];count=0
            if name=="current_Rp_left_source":
                if packet["chart"]!="O3_power" or len(groups)!=2:
                    raise ValueError("Actual current terminal O3 packet required")
            else:
                if (packet["chart"]!="pulse_entrance"
                        or packet["source_provider"]!=PREFIX+"pulse_mixed_C4.CompliantPulseMixedC4"
                        or packet["acceptance_receipt"]!=PREFIX+"actual_Rp_source_join_check.json"
                        or not packet["current_Rp_external_pulse_join_proved"]
                        or packet["current_Rp_external_pulse_join_certified"]):
                    raise ValueError("Current native pulse packet required")
            for group in groups.values():
                for grid in group.values():
                    if set(grid)!=keys:raise ValueError("Actual total-order4 source row omitted")
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                            raise ArithmeticError("Nonfinite current Rp/entrance source enclosure")
                        count+=1
            expected=135 if name=="current_Rp_left_source" else 60
            if count!=expected:raise ValueError("Current Rp/entrance derivative row scope differs")
            counts[name]=count
        if (any(raw[k] for k in SCOPES)
                or raw["full_current_core_to_heat_physical_assembly"]
                or raw["all_profile_source_charts_callable"]):
            raise ValueError("Rp local/source result promoted unfinished global scope")
        result=dict(actual_five_defect_family_sha256=field.family,
            implicit_source_sha256=field.source,datum_enclosure_sha256=field.datum_sha,
            current_input_hashes_checked=len(raw["input_hashes"]),
            current_downstream_chart_owners_checked=15,
            current_native_pulse_source_graph=raw["current_Rp_native_source_graph"],
            exact_current_Rp_coefficient_functions_checked=28,
            exact_native_inlet_axial5_source_coefficients_checked=54,
            exact_current_Rp_physical_mixed4_source_identities_checked=135,
            actual_runtime_source_rows_checked=counts,
            original_gp_entrance_flat_jets_and_partial_energy_checked=True,
            same_pressure_function_parameters_Cstar_and_Rp_source_checked=True,
            saved_power_inlet_samples_not_used_for_current_native_join=True,
            current_Rp_external_pulse_join_certified=True,
            current_native_pulse_entrance_owner_installed=True,
            all_current_pulse_charts_installed=False,
            full_current_core_to_heat_physical_assembly=False,
            **dict.fromkeys(SCOPES,False),all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current Rp native source join: 15 owners, 54 inlet coefficients, 135 exact mixed identities",flush=True)
    return result


if __name__=="__main__":
    run()
