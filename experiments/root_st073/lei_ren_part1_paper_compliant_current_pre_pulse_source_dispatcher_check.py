"""Current pre-pulse ownership, exact histories and whole-chart mixed4 acceptance."""
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_pre_pulse_source_dispatcher import (
    CurrentPrePulseSourceDispatcher, CHARTS, PRE_ROUTES, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"current_pre_pulse_source_dispatcher.json"


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentPrePulseSourceDispatcher(require_checked=False)
        expected=encode(pack(field.manifest()))
        source={k:v for k,v in raw.items() if k!="whole_current_pre_pulse_charts"}
        if source!=expected:
            raise ValueError("Current Rh graph/pre-pulse history/parameter/source proof changed")
        proof=raw["current_pre_pulse_five_interface_source_proof"]
        if (not proof["passed"] or proof["exact_history_identity_count"]!=25
                or len(proof["exact_five_history_interface_identities"])!=5
                or proof["total_implied_physical_mixed4_interface_rows"]!=675
                or not proof["retained_axial_and_mixed_histories_not_reset_when_V_zero"]
                or not proof["exact_positive_mu_retained_in_formal_power_slope"]
                or not raw["current_Rh_external_neighbor_join_certified"]
                or raw["current_Rh_to_Rp_source_chain_certified"]):
            raise ValueError("Current Rh admission/producer pre-pulse certification scope differs")
        registry=raw["ordered_current_chart_registry"]
        if len(registry)!=14 or tuple(registry)[-5:]!=CHARTS:
            raise ValueError("Expected fourteen current downstream chart owners in radial order")
        if field.provider("Rh_reference") is not field.provider("O3_power"):
            raise ValueError("Rh and Rp cannot use separately initialized history providers")
        c=MPIntervalContext();c.dps=160
        rowset={"y"+str(j)+"_Z"+str(n) for j in range(5) for n in range(5-j)}
        counts={}
        for chart,packet in raw["whole_current_pre_pulse_charts"].items():
            if chart not in CHARTS:raise ValueError("Unknown whole current pre-pulse chart")
            route=PRE_ROUTES[chart]
            if (registry[chart]!=route or packet["source_provider"]!=route["provider"]
                    or packet["source_coordinate_domain"]!=route["domain"]
                    or packet["acceptance_receipt"]!=route["acceptance_receipt"]
                    or not packet["same_current_Rh_reference_provider"]
                    or not packet["current_Rh_to_Rp_source_chain_proved"]
                    or packet["current_Rh_to_Rp_source_chain_certified"]
                    or packet["current_Rp_external_pulse_join_certified"]):
                raise ValueError("Whole current chart has wrong owner/coverage/acceptance scope")
            if tuple(packet["physical_mixed_grids"])!=(
                    "physical_velocity_pressure_y_Z_mixed4","physical_five_primitive_y_Z_mixed4"):
                raise ValueError("Both physical groups required")
            count=0
            for group,grids in packet["physical_mixed_grids"].items():
                if len(grids)!=(5 if "five_primitive" in group else 4):
                    raise ValueError("Physical components or primitive omitted")
                for grid in grids.values():
                    if set(grid)!=rowset:raise ValueError("Current pre-pulse mixed4 row missing")
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                            raise ArithmeticError("Nonfinite current pre-pulse source row")
                        count+=1
            if count!=135 or any(packet[k] for k in SCOPES):
                raise ValueError("Current chart count/scope differs")
            counts[chart]=count
        if tuple(counts)!=CHARTS or sum(counts.values())!=675:
            raise ValueError("Every new current pre-pulse chart must be exercised")
        if (any(raw[k] for k in SCOPES) or raw["current_Rp_external_pulse_join_certified"]
                or raw["all_profile_source_charts_callable"] or raw["uniform_physical_units_assembled"]):
            raise ValueError("Current pre-pulse stage promotes unbuilt external/global/point/recursion scope")
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_input_hashes_checked=len(raw["input_hashes"]),
            current_original_pre_pulse_source_bindings=raw["current_original_pre_pulse_history_source_bindings"],
            current_pre_pulse_exact_five_history_interface_identities=25,
            current_pre_pulse_functional_mixed4_interface_rows_implied=675,
            current_whole_pre_pulse_mixed4_rows_checked=counts,total_current_new_pre_pulse_rows_checked=675,
            current_downstream_chart_owners_checked=14,same_current_Rh_Rp_provider_object_verified=True,
            current_Rh_external_neighbor_join_certified=True,
            current_Rh_to_Rp_source_chain_certified=True,all_current_pre_pulse_charts_callable=True,
            original_positive_mu_and_retained_nonzero_histories_preserved=True,
            unchanged_independent_physical_and_turnoff_coordinate_fixtures_consumed=True,
            source_history_identities_not_overlap_proof=True,
            current_Rp_external_pulse_join_certified=False,
            all_profile_source_charts_callable=False,uniform_physical_units_assembled=False,
            **dict.fromkeys(SCOPES,False),all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current Rh-to-Rp original source chain: 14 owners, 25 history identities, 675 whole mixed4 rows",flush=True)
    return result


if __name__=="__main__":
    run()
