"""Focused current implicit-source patch mixed4 and Rm/Rh interface checks."""
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_actual_feedback_patch_mixed_C4 import (
    CompliantActualFeedbackPatchMixedC4, current_mixed_source_bindings, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import canonical_source, encode_parent
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"actual_feedback_patch_mixed_C4.json"
GROUPS=("physical_velocity_pressure_x_Z_mixed4","physical_velocity_pressure_y_Z_mixed4",
        "physical_five_primitive_x_Z_mixed4")


def exact(a,b,label):
    if canonical_source(a)!=canonical_source(b):
        raise ValueError("Current patch mixed source differs: "+label)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        provider=CompliantActualFeedbackPatchMixedC4();c=provider.ctx
        if (raw["datum_enclosure_sha256"]!=provider.patch.core.datum.datum_sha
                or raw["current_patch_mixed_source_bindings"]!=current_mixed_source_bindings()):
            raise ValueError("Current patch source/datum/units changed")
        generic=provider.generic_check;joins=generic["exact_coordinate_and_join_checks"]
        if (not joins["actual_Rm_join_uses_open_zero_correction_neighborhood"]
                or not joins["actual_Rh_join_uses_full_support_and_same_implicit_functional_closure"]):
            raise ValueError("Original functional patch endpoint proof missing")
        current_name=PREFIX+"actual_feedback_moment_patch_check.json"
        current=json.loads((HERE/current_name).read_bytes());_verify_hashes(current)
        if not (current["all_passed"] and current["current_actual_patch_implicit_axial5_recomputed"]
                and current["current_five_functional_terminal_identities_connected"]
                and current["current_Rm_source_parent_and_P0_retained"]):
            raise ValueError("Current endpoint joins need current accepted implicit source family")
        packets=[raw[k] for k in ("whole_actual_patch","actual_Rm_inlet","actual_Rh_exit")]
        packets+=raw["source_bump_edge_packets"]+raw["interior_packets"]
        counts={g:0 for g in GROUPS}
        for packet in packets:
            for group in GROUPS:
                grids=packet[group]
                if len(grids)!=(5 if group==GROUPS[2] else 4):
                    raise ValueError("Current physical velocity/primitive omitted")
                prefix="y" if group==GROUPS[1] else "x"
                expected={prefix+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)}
                for grid in grids.values():
                    if set(grid)!=expected:
                        raise ValueError("Current patch mixed4 grid incomplete")
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                            raise ArithmeticError("Nonfinite current physical patch row")
                        counts[group]+=1
            if (not packet["same_actual_coefficient_family_and_P0_retained"]
                    or not packet["ordinary_derivatives_not_Taylor_radial_coefficients"]
                    or not packet["radial_prefactors_differentiated_before_mixed_grid"]):
                raise ValueError("Current patch history or normalization lost")
        for key,x in (("actual_Rm_inlet",c.mpf(1)),("actual_Rh_exit",c.exp(1))):
            live=provider.evaluate(x,[-1,1])
            exact(raw[key]["actual_inherited_patch_packet"],encode_parent(live["actual_inherited_patch_packet"]),key+" current parent")
            for group in GROUPS:
                exact(raw[key][group],encode_parent(live[group]),key+" "+group)
        for edge,packet in zip((49,51,59,61,69,71),raw["source_bump_edge_packets"]):
            if endpoints(read_interval(c,packet["x"]))!=endpoints(c.mpf(edge)/40):
                raise ValueError("Original compact support edge changed")
            for rows in packet["actual_gamma_x_derivatives"]:
                for value in rows:
                    lo,hi=endpoints(read_interval(c,value))
                    if not lo<=0<=hi:
                        raise ArithmeticError("Beta support-edge enclosure excludes its source zero")
        # Exact rational support edges and the accepted flat endpoint limits
        # prove source zeros. Rounded edge boxes can enclose a neighborhood;
        # their widths are not a counterexample or a functional proof.
        if endpoints(c.exp(1))[0]<=mp.mpf(71)/40:
            raise ArithmeticError("Original Rh beyond-support ordering changed")
        rh=raw["actual_Rh_exit"]["actual_inherited_patch_packet"]
        if not rh["terminal_refinement_is_functional_identity_not_moment_reset"]:
            raise ValueError("Current Rh closure cannot reset source histories")
        zeros=0
        for row in rh["actual_normalized_five_defects"]:
            for value in row["coefficients"]:
                if endpoints(read_interval(c,value))!=(mp.mpf(0),mp.mpf(0)):
                    raise ArithmeticError("Current five terminal functions not exactly closed")
                zeros+=1
        for key in ("full_implicit_leading_inputs_recomputed","actual_point_moment_history_recovered",
            "global_completed_tensor_admissibility","full_inner_interfaces_certified",
            "full_cartesian_vector_derivatives_certified","admissible_stress_lift_constructed","temporal_recursion"):
            if raw[key]:
                raise ValueError("Unbuilt global/point/recursion scope promoted: "+key)
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            datum_enclosure_sha256=provider.patch.core.datum.datum_sha,
            current_input_hashes_checked=len(raw["input_hashes"]),
            current_patch_mixed_source_bindings=raw["current_patch_mixed_source_bindings"],
            current_physical_mixed_bounds_checked=counts,
            current_Rm_Rh_endpoint_provider_rows_replayed=390,
            exact_coordinate_and_support_fixture_consumed=PREFIX+"actual_patch_mixed_C4_check.json",
            unchanged_independent_physical_fixture=generic["independent_physical_fixture"],
            original_flat_beta_support_edges_checked=6,
            current_five_terminal_functional_zero_rows_checked=zeros,
            current_Rm_source_transport_and_P0_retained=True,
            current_Rh_same_unique_implicit_function_and_full_support_retained=True,
            current_actual_moment_patch_installed=True,current_actual_patch_mixed4_available=True,
            current_actual_patch_implicit_axial5_recomputed=True,
            current_Rm_and_Rh_functional_mixed4_joins_certified=True,
            current_five_functional_terminal_identities_connected=True,
            source_joins_not_proved_by_interval_overlap=True,
            full_implicit_leading_inputs_recomputed=False,actual_point_moment_history_recovered=False,
            global_completed_tensor_admissibility=False,full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current implicit-source patch physical mixed4 and functional Rm/Rh joins",flush=True)
    return result


if __name__=="__main__":
    run()
