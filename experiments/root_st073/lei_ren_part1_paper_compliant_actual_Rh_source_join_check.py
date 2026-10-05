"""Focused acceptance for the current external Rh reference interface."""
import json
from pathlib import Path

from lei_ren_part1_paper_compliant_actual_Rh_source_join import (
    build, RH_ROUTE, SCOPES, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"actual_Rh_source_join.json"


def run():
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=build()
    if raw!=expected:
        raise ValueError("Current Rh defining source, provider, proof or reference rows changed")
    pressure=raw["current_Rh_pressure_defining_function_proof"]
    proof=raw["current_external_Rh_functional_join"]
    if (not pressure["passed"] or pressure["compliant_epsilon"]!=".001"
            or not pressure["exact_same_inherited_normalized_jets_callable"]
            or not pressure["same_fourteen_implicit_stage_sources"]
            or not pressure["common_flatten_function_retained_not_replaced_by_zero"]
            or pressure["exact_first_six_Taylor_defining_coefficients_projected"]!=6
            or pressure["q_recurrence_symbolic_rows_checked"]!=6):
        raise ValueError("Current Rh P0 defining-function/first-six projection required")
    rowset={"y"+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)}
    rows=proof["physical_mixed4_rows_implied_by_exact_function_identities"]
    if (not proof["passed"] or len(rows)!=9 or any(set(v)!=rowset for v in rows.values())
            or proof["total_implied_physical_mixed4_rows"]!=135
            or not proof["current_unique_implicit_map_full_weight_closure_used"]
            or not proof["normalization_multipliers_frozen_at_basepoint_not_differentiated_in_Z"]):
        raise ValueError("All nine physical source functions/mixed4 rows required")
    registry=raw["ordered_current_chart_registry"]
    if (len(registry)!=9 or tuple(registry)[-1]!="Rh_reference"
            or registry["Rh_reference"]!=RH_ROUTE
            or raw["current_Rh_external_neighbor_join_certified"]
            or not raw["current_Rh_external_neighbor_source_join_proved"]
            or not raw["native_eight_chart_Rh_false_flag_preserved"]):
        raise ValueError("Current external Rh needs its own ninth owner/receipt")
    native=json.loads((HERE/(PREFIX+"current_matched_source_dispatcher.json")).read_bytes())
    _verify_hashes(native)
    if native["current_Rh_external_neighbor_join_certified"]:
        raise ValueError("Historical eight-chart external Rh scope must remain false")
    packet=raw["whole_current_Rh_reference_chart"]
    if (packet["source_provider"]!=RH_ROUTE["provider"]
            or packet["acceptance_receipt"]!=RH_ROUTE["acceptance_receipt"]
            or packet["source_coordinate_domain"]!=[-5,0]
            or packet["current_Rh_external_neighbor_join_certified"]
            or not packet["current_Rh_external_neighbor_source_join_proved"]):
        raise ValueError("Current whole Rh reference owner/coverage differs")
    import mpmath as mp
    c=mp.iv;c.dps=160;count=0
    for group,grids in packet["physical_mixed_grids"].items():
        if len(grids)!=(5 if "five_primitive" in group else 4):
            raise ValueError("Whole current reference components/primitive missing")
        for grid in grids.values():
            if set(grid)!=rowset:raise ValueError("Whole current reference mixed4 row missing")
            for value in grid.values():
                lo,hi=endpoints(read_interval(c,value))
                if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                    raise ArithmeticError("Nonfinite whole current Rh reference source row")
                count+=1
    if count!=135:raise ValueError("Expected 135 whole reference mixed4 rows")
    if (any(raw[k] for k in SCOPES) or any(packet[k] for k in SCOPES)
            or raw["all_profile_source_charts_callable"] or raw["uniform_physical_units_assembled"]):
        raise ValueError("External Rh admission cannot promote global/point/recursion scope")
    result=dict(actual_five_defect_family_sha256=raw["actual_five_defect_family_sha256"],
        implicit_source_sha256=raw["implicit_source_sha256"],
        datum_enclosure_sha256=raw["datum_enclosure_sha256"],
        current_input_hashes_checked=len(raw["input_hashes"]),
        current_Rh_pressure_source_AST_bindings=pressure["source_AST_bindings"],
        compliant_pressure_epsilon=".001",same_inherited_P0_defining_function_verified=True,
        exact_first_six_pressure_Taylor_source_coefficients_projected=6,
        q_Taylor_recurrence_symbolic_rows_checked=6,
        first_six_enclosure_prefix_diagnostic_rows_checked=18,
        physical_functions_with_exact_Rh_source_identity=9,
        physical_mixed4_rows_implied_by_exact_source_identities=135,
        whole_current_Rh_reference_mixed4_rows_checked=count,
        current_downstream_chart_owner_count=9,current_Rh_reference_owner_installed=True,
        current_Rh_unique_implicit_map_full_support_closure_retained=True,
        current_Rh_external_neighbor_join_certified=True,
        native_eight_chart_Rh_false_flag_preserved=True,
        pressure_enclosures_not_selected_as_point_values=True,
        source_equality_proof_not_enclosure_overlap=True,
        all_profile_source_charts_callable=False,uniform_physical_units_assembled=False,
        **dict.fromkeys(SCOPES,False),all_passed=True,
        input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current external Rh: shared P0 function, 6 projected source coefficients, 9 physical identities, 135 mixed4 rows",flush=True)
    return result


if __name__=="__main__":
    run()
