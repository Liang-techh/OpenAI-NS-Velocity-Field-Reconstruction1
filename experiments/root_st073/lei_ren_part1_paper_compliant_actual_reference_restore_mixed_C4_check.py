"""Current source/row acceptance for correlated-E reference/restoration.

Unchanged independent physical/restoration fixtures are reused by hashes.
The current Rsh full mixed4 join remains a separate open interface gate.
"""
import json
import math
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import (
    CompliantActualReferenceRestoreMixedC4, current_E_source_bindings, current_centered_core,
    accepted, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import (
    canonical_source, encode_parent)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import norm
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"actual_reference_restore_mixed_C4.json"
GROUPS=("physical_velocity_pressure_y_Z_mixed4","physical_five_primitive_y_Z_mixed4")


def exact(left,right,label):
    if canonical_source(left)!=canonical_source(right):
        raise ValueError("Current source transfer differs: "+label)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        provider=CompliantActualReferenceRestoreMixedC4();c=provider.ctx
        family,source=provider.family,provider.source
        name=PREFIX+"reference_restore_profiles_check.json"
        profiles=accepted(name,family,source,"actual_axial_restoration_axial5_enclosures_available")
        generic=provider.original_check
        if (not profiles["restoration_fixture"]["passed"] or not profiles["structural_identities"]["passed"]
                or not generic["structural_checks"]["passed"] or not generic["independent_physical_fixture"]["passed"]):
            raise ValueError("Current unchanged generic fixture receipts required")
        if (raw["current_centered_E_source_bindings"]!=current_E_source_bindings(provider)
                or raw["datum_enclosure_sha256"]!=provider.reference.core.datum.datum_sha):
            raise ValueError("Current E source or datum binding differs")
        graph=raw["shared_exact_axial_source"]
        if graph!=provider.shared_axial_source or graph!=provider.long_mixed.shared_axial_source:
            raise ValueError("Current signed axial source graph changed")
        if graph["E_V110_minus_4Z"]["args"]!=graph["V110"]["args"][1:]:
            raise ValueError("Centered E lost common exact source references")
        inp=provider.reference.inputs([-1,1])
        parent=provider.long_mixed.reshape.evaluate([-1,1],phase=1)
        exact(raw["current_actual_Rsh_parent"],encode_parent(parent),"actual Rsh parent")
        exact(encode_parent(inp["parent"]),encode_parent(parent),"reference current parent")
        exact(raw["current_centered_E_component_enclosures"],
              encode_parent(inp["current_E_component_enclosures"]),"current E component source boxes")
        exact(raw["current_centered_E_axial5"],encode_parent(list(inp["E"].coefficients)),"current E enclosure")
        core_packet=current_centered_core(provider.long_mixed.history,[-1,1])
        if core_packet["baseline_rows_checked"]!=6 or not core_packet["defining_field_not_changed"]:
            raise ValueError("Current fresh affine baseline removal is not source-bound")
        incoming=provider.long_mixed.history.upstream.packet([-1,1],1,"macro")
        switch=provider.long_mixed.history.switch.inputs([-1,1])
        parts=inp["current_E_component_enclosures"]
        core_E=core_packet["E"]
        exact(encode_parent(parts["core"]),encode_parent(list(core_E.coefficients)),"centered core source")
        exact(encode_parent(parts["bridge"]),encode_parent(incoming["actual_delta_V_axial5"]),"actual bridge increment")
        exact(encode_parent(parts["first_switch"]),encode_parent(switch["velocity_increment"]),"actual first-switch increment")
        if len(parts["raw_sum_before_C2_theorem"])!=6:
            raise ValueError("Current E axial5 coefficient count lost")
        ledger=provider.reference.reshape.switch.bridge.records["K1_ledger"]
        gate=provider.reference.reshape.switch.bridge.records["global_exit_certificate"]
        budget=read_interval(c,gate["rho_bridge_C2_upper"])
        total=sum((norm(c,parts["bridge"][k]+parts["first_switch"][k])*math.factorial(k)
                   for k in range(3)),c.mpf(0))
        if endpoints(total)[1]>endpoints(budget)[0]:
            raise ValueError("Current actual increment C2 budget not proved")
        if (not gate["actual_whole_axis_Ra_R110_relaxed_cone_analytically_certified"]
                or not gate["short_switch_moments_included_by_exact_velocity_averaging"]
                or "includes both short-switch C2 velocity increments" not in gate["proof"]["axial_budget"]):
            raise ValueError("Original source theorem does not cover current route")
        rho=read_interval(c,ledger["rho_core_C2_bound"])+budget
        exact(encode_parent(inp["rho"]),encode_parent(rho),"same source C2 bound")
        for k in range(3):
            source_cover=c.mpf([-endpoints(rho)[1],endpoints(rho)[1]])/math.factorial(k)
            if k==0:
                source_cover+=provider.reference.core.j
            lo,hi=endpoints(inp["E"][k]);boundlo,boundhi=endpoints(source_cover)
            if lo<boundlo or hi>boundhi:
                raise ValueError("Final E source theorem envelope lost")
        keys=("whole_reference","whole_restoration","whole_postrestore","actual_Rsh_exit",
              "actual_Rz_reference_side","actual_restore_inlet","actual_restore_exit",
              "actual_postrestore_inlet","actual_Rm_exit")
        packets=[raw[k] for k in keys]+raw["interior_packets"]
        expected={"y"+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)}
        counts={group:0 for group in GROUPS}
        for packet in packets:
            for group in GROUPS:
                if len(packet[group])!=(4 if group==GROUPS[0] else 5):
                    raise ValueError("Physical source field/primitive group omitted")
                for grid in packet[group].values():
                    if set(grid)!=expected:
                        raise ValueError("Current restoration mixed4 row domain incomplete")
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                            raise ArithmeticError("Nonfinite current reference/restore row")
                        counts[group]+=1
            packet_inp=provider.reference.inputs(read_interval(c,packet["Z"]))
            exact(packet["actual_inherited_axial5_packet"]["pressure_axis_axial5_coefficients"],
                  encode_parent(list(packet_inp["original_axis_pressure"].coefficients)),"canonical P0")
            if any(packet[k] for k in ("Rsh_reshape_join_certified","full_inner_interfaces_certified",
                "full_cartesian_vector_derivatives_certified","admissible_stress_lift_constructed","temporal_recursion")):
                raise ValueError("Unbuilt Rsh/global scope promoted")
        for a,b in (("actual_Rz_reference_side","actual_restore_inlet"),
                    ("actual_restore_exit","actual_postrestore_inlet")):
            for group in GROUPS:
                exact(raw[a][group],raw[b][group],"current exact two-sided "+a)
        exact(raw["actual_Rsh_exit"]["actual_inherited_axial5_packet"]["actual_centered_moment_axial5_coefficients"],
              encode_parent({k:list(row.coefficients) for k,row in inp["source_centered_Rsh"].items()}),
              "current six centered Rsh histories")
        for key in ("actual_restore_exit","actual_postrestore_inlet","actual_Rm_exit"):
            rows=raw[key]["actual_inherited_axial5_packet"]["Uz_actual_axial5_coefficients"]
            true=list((inp["z"]*4).coefficients)
            exact(rows,encode_parent(true),"exact terminal 4Z")
            for rowkey,value in raw[key][GROUPS[0]]["Uz"].items():
                k,n=[int(v[1:]) for v in rowkey.split("_")]
                target=(inp["z"]*4)[n]*math.factorial(n) if k==0 else c.mpf(0)
                exact(value,encode_parent(target),"terminal Uz mixed derivative")
        defects=raw["current_unpatched_five_defects_at_Rm"]["actual_normalized_five_defect_axial5_coefficients"]
        if len(defects)!=5 or any(len(row)!=6 for row in defects):
            raise ValueError("Current five prepatch defect axial5 shape lost")
        for row in defects:
            for value in row:
                lo,hi=endpoints(read_interval(c,value))
                if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                    raise ArithmeticError("Nonfinite current prepatch defect")
        for key in ("current_Rsh_source_functional_join_certified","current_actual_moment_patch_installed",
                    "full_implicit_leading_inputs_recomputed","actual_point_moment_history_recovered",
                    "global_completed_tensor_admissibility","full_inner_interfaces_certified",
                    "full_cartesian_vector_derivatives_certified","admissible_stress_lift_constructed","temporal_recursion"):
            if raw[key]:
                raise ValueError("Unbuilt current scope promoted: "+key)
        result=dict(actual_five_defect_family_sha256=family,implicit_source_sha256=source,
            datum_enclosure_sha256=provider.reference.core.datum.datum_sha,
            current_input_hashes_checked=len(raw["input_hashes"]),
            original_profile_fixture_receipt_consumed=name,
            original_mixed_fixture_receipt_consumed=PREFIX+"reference_restore_mixed_C4_check.json",
            current_centered_E_source_bindings=raw["current_centered_E_source_bindings"],
            current_core_source_ordinary_rows_checked=6,
            current_actual_increment_C2_budget_verified=True,
            current_E_source_C2_theorem_verified_before_intersection=True,
            current_Rsh_parent_and_six_centered_histories_retained=True,
            current_mixed_bounds_checked=counts,
            current_exact_Rz_and_restore_exit_two_sided_rows_checked=270,
            current_original_restore_end_exact_4Z=True,
            current_unpatched_five_defect_coefficients_checked=30,
            current_correlated_E_installed=True,
            current_actual_reference_restore_mixed4_available=True,
            current_Rsh_source_functional_join_certified=False,
            current_actual_moment_patch_installed=False,
            full_implicit_leading_inputs_recomputed=False,
            actual_point_moment_history_recovered=False,
            global_completed_tensor_admissibility=False,
            full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    print("PASS current correlated E, reference/restoration mixed4 and exact final 4Z",flush=True)
    return result


if __name__=="__main__":
    run()
