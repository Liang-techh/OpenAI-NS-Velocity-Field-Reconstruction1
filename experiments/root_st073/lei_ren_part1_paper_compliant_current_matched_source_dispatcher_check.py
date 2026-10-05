"""Current eight-chart ownership, source composition and whole-domain rows."""
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_matched_source_dispatcher import (
    CurrentMatchedSourceDispatcher, CURRENT_ROUTES, Rm_functional_identity,
    GATES, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_source_dispatcher import ROUTES
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import canonical_source, encode_parent
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"current_matched_source_dispatcher.json"


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentMatchedSourceDispatcher(require_checked=False)
        live=field.manifest();c=field.anchor.ctx
        if (raw["ordered_current_chart_registry"]!=live["ordered_current_chart_registry"]
                or raw["datum_enclosure_sha256"]!=field.datum_sha
                or raw["current_route_specific_acceptance_gates"]!=encode_parent(GATES)
                or raw["provider_graph_identity"]!=field.provider_graph_identity()
                or not all(field.provider_graph_identity().values())):
            raise ValueError("Current chart owner/common graph/source units differ")
        if raw["current_external_Rm_functional_join"]!=Rm_functional_identity():
            raise ValueError("Current external Rm source/coordinate identity changed")
        count=0;by_chart={}
        packets=raw["whole_current_chart_evaluations"]
        if tuple(packets)!=tuple(CURRENT_ROUTES) or len(packets)!=8:
            raise ValueError("All eight current charts must be exercised in original radial order")
        for chart,spec in CURRENT_ROUTES.items():
            packet=packets[chart];provider=field.provider(chart)
            stem,cls,method,description,domain,receipt,extra=spec
            if (type(provider).__name__!=cls
                    or packet["source_provider"]!=PREFIX+stem+"."+cls
                    or packet["acceptance_receipt"]!=PREFIX+receipt+".json"
                    or packet["source_coverage_coordinate"]!=description
                    or packet["source_coordinate_domain"]!=domain
                    or packet["actual_five_defect_family_sha256"]!=field.family
                    or packet["implicit_source_sha256"]!=field.source
                    or packet["datum_enclosure_sha256"]!=field.datum_sha
                    or not packet["same_current_provider_graph"]):
                raise ValueError("Current per-chart provider/acceptance/coordinate mismatch")
            if (spec[2:5]!=ROUTES[chart][2:5] or spec[6]!=ROUTES[chart][6]):
                raise ValueError("Original chart method/domain/order must be preserved")
            rows=0
            for group,grids in packet["physical_mixed_grids"].items():
                prefix="s" if "_phase_" in group else "x" if "_x_" in group else "y"
                expected={prefix+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)}
                if len(grids)!=(5 if "five_primitive" in group else 4):
                    raise ValueError("Current physical component/primitive omitted")
                for grid in grids.values():
                    if set(grid)!=expected:
                        raise ValueError("Current chart mixed4 row set incomplete")
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                            raise ArithmeticError("Nonfinite current whole-domain field row")
                        rows+=1
            if rows!=(195 if chart=="actual_patch" else 135):
                raise ValueError("Current chart mixed4 group count changed")
            by_chart[chart]=rows;count+=rows
            if any(packet[k] for k in SCOPES):
                raise ValueError("Current source dispatcher promoted unfinished physical/global scope")
        proof=raw["current_external_Rm_functional_join"]
        if (not proof["passed"] or len(proof["exact_physical_function_identities"])!=9
                or proof["total_implied_physical_mixed4_rows"]!=135
                or not proof["fixed_normalization_multipliers_not_differentiated_in_Z"]):
            raise ValueError("Current external Rm functional proof scope incomplete")
        # The common graph makes source-general accepted local identities
        # applicable to the actual neighboring chart owners. Grids with
        # different fixed physical units are not compared as raw boxes.
        if (not field.receipts["actual_long_reshape_mixed_C4_check"]["R110_postswitch_power_reshape_local_mixed4_join_certified"]
                or not field.receipts["actual_Rsh_source_join_check"]["current_Rsh_source_functional_join_certified"]
                or not field.receipts["actual_reference_restore_mixed_C4_check"]["current_original_restore_end_exact_4Z"]
                or not field.receipts["actual_feedback_patch_mixed_C4_check"]["current_Rm_source_transport_and_P0_retained"]):
            raise ValueError("Current neighboring source interfaces lack accepted owners")
        if (raw["current_Rh_external_neighbor_join_certified"] or raw["all_profile_source_charts_callable"]
                or raw["uniform_physical_units_assembled"] or any(raw[k] for k in SCOPES)):
            raise ValueError("Unbuilt Rh/core-to-heat/Cartesian/global/recursion scope promoted")
        try:
            field.provider("core")
        except ValueError:
            pass
        else:
            raise ValueError("Eight-chart dispatcher cannot silently load a legacy core owner")
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_input_hashes_checked=len(raw["input_hashes"]),
            current_chart_owners_checked=8,current_provider_graph_identity=field.provider_graph_identity(),
            current_route_specific_gates_checked=GATES,
            whole_current_chart_mixed4_rows_checked=by_chart,total_current_physical_mixed4_rows_checked=count,
            current_external_Rm_exact_physical_source_identities=9,
            current_external_Rm_implied_mixed4_rows=135,
            current_R110_Rsh_Rz_restore_exit_Rm_functional_mixed4_joins_composed=True,
            same_actual_source_graph_P0_and_cache_objects_used=True,
            source_joins_not_proved_by_enclosure_overlap=True,
            current_Rh_outgoing_unique_implicit_reference_closure_available=True,
            current_Rh_external_neighbor_join_certified=False,
            all_current_matched_charts_exercised=True,all_current_matched_charts_callable=True,
            all_profile_source_charts_callable=False,uniform_physical_units_assembled=False,
            **dict.fromkeys(SCOPES,False),all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS one current matched source chain: 8 owners,1140 mixed rows,current external Rm join",flush=True)
    return result


if __name__=="__main__":
    run()
