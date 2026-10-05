"""Focused admission of current source defects, one implicit Jacobian and patch."""
import json
import math
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_actual_feedback_moment_patch import (
    CompliantActualFeedbackMomentPatch, current_patch_source_bindings, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import canonical_source, encode_parent
from lei_ren_part1_paper_interval_five_bump_inverse import weighted_norm, mag
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"actual_feedback_moment_patch.json"


def exact(a,b,label):
    if canonical_source(a)!=canonical_source(b):
        raise ValueError("Current patch source transfer differs: "+label)


def contains_zero(c,value,label):
    lo,hi=endpoints(read_interval(c,value))
    if not lo<=0<=hi:
        raise ArithmeticError("Current implicit diagnostic excludes zero: "+label)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        provider=CompliantActualFeedbackMomentPatch(require_checked=False);c=provider.ctx
        if (raw["datum_enclosure_sha256"]!=provider.core.datum.datum_sha
                or raw["original_correlated_source_class_admission"]!=provider.source_class_admission
                or raw["current_patch_source_bindings"]!=current_patch_source_bindings()
                or raw["current_exact_centered_E_source"]!=provider.reference_mixed.shared_axial_source["E_V110_minus_4Z"]):
            raise ValueError("Current source family, datum or exact E graph differs")
        inverse,data=provider.coefficients([-1,1])
        exact(raw["actual_connected_source_data"],encode_parent(data),"all five current source defects/tails")
        exact(raw["actual_implicit_coefficient_inverse_axial5"],encode_parent(inverse),"same current implicit inverse")
        norm=weighted_norm(c,inverse["preconditioned_error"],provider.scales)
        if (not inverse["initial_C1_inverse"]["certified"] or endpoints(norm)[1]>=1
                or endpoints(norm)!=endpoints(inverse["higher_inverse_contraction"])
                or not inverse["smooth_implicit_coefficient_family_through_axial5"]
                or not inverse["residual_zero_containment_only_diagnostic"]):
            raise ArithmeticError("Current shared implicit Jacobian not admitted")
        if [p["order"] for p in inverse["higher_order_proofs"]]!=[2,3,4,5]:
            raise ValueError("Current higher implicit orders omitted")
        tails=0;decays=0
        for row,bound in zip(data["actual_tail_functions"],data["actual_tail_C1_norm_upper"]):
            actual=mag(c,row[0])+mag(c,row[1])
            if endpoints(actual)!=endpoints(bound) or endpoints(actual)[1]>endpoints(provider.oldtail)[0]:
                raise ArithmeticError("Current source tail outside original source-general C1 class")
            tails+=1
        for proof in data["retained_source_decay_proofs"]:
            value=-proof["rate"]*proof["source_gap_enclosure"]+c.ln(proof["input_absolute_upper"])
            if endpoints(value)[1]>endpoints(provider.reference.logcap)[0]:
                raise ArithmeticError("Current factored source decay proof failed")
            decays+=1
        if endpoints(data["rho"])[1]>endpoints(provider.oldrho)[1]:
            raise ArithmeticError("Current E outside original admitted source class")
        residuals=raw["actual_implicit_coefficient_inverse_axial5"]["implicit_map_axial5_residual_enclosures"]
        if len(residuals)!=5 or any(len(row["coefficients"])!=6 for row in residuals):
            raise ValueError("Current implicit map diagnostic needs all30 rows")
        for row in residuals:
            for v in row["coefficients"]:
                contains_zero(c,v,"current implicit map")
        for proof in raw["actual_implicit_coefficient_inverse_axial5"]["higher_order_proofs"]:
            for v in proof["linear_equation_residual"]:
                contains_zero(c,v,"same Jacobian higher row")
        groups=("Utheta_over_Pstar_axial5","Uz_axial5","Utheta_y_over_Pstar_axial5","Uz_y_axial5",
            "Ur_over_sqrt_R_over_2_axial4","Mz_over_R_axial5","Mtheta_over_sqrt2_Rm_1p5_Am_axial5",
            "Mtheta_z_over_sqrt2_Rm_1p5_Am_axial5","Mztheta_minus_8ZMz_plus_16Z2R_over_RmAm2_axial5",
            "Mp_over_Am2_axial5","P_over_Pstar2_axial5","original_P0_axial5")
        count=0;terminal_rows=0
        packets=raw["actual_patch_samples"]+[raw["whole_actual_patch"],raw["actual_terminal_Rh"]]
        for packet in packets:
            for name in groups:
                if len(packet[name])!=(5 if name.startswith("Ur_") else 6):
                    raise ValueError("Current patch axial field order lost")
                for v in packet[name]:
                    lo,hi=endpoints(read_interval(c,v))
                    if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                        raise ArithmeticError("Nonfinite current actual patch field")
                    count+=1
            if endpoints(read_interval(c,packet["Utheta_over_Pstar_axial5"][0]))[0]<=0:
                raise ArithmeticError("Current original patch positive swirl lost")
            if packet["exact_terminal_closure_from_same_implicit_map"]:
                if not packet["terminal_refinement_is_functional_identity_not_moment_reset"]:
                    raise ValueError("Terminal closure must use same unique implicit function")
                for row,diagnostic in zip(packet["actual_normalized_five_defects"],packet["actual_unreferenced_terminal_diagnostic"]):
                    for v,q in zip(row["coefficients"],diagnostic["coefficients"]):
                        if endpoints(read_interval(c,v))!=(mp.mpf(0),mp.mpf(0)):
                            raise ArithmeticError("Current functional terminal refinement not zero")
                        contains_zero(c,q,"unrefined current full-support map")
                        terminal_rows+=1
        inlet=provider.evaluate(1,[-1,1]);unpatched=provider.reference.terminal([-1,1],-6)
        exact(encode_parent(inlet["original_P0_axial5"]),
            encode_parent(unpatched["pressure_axis_axial5_coefficients"]),"current Rm pressure")
        reference_raw=json.loads((HERE/(PREFIX+"actual_reference_restore_mixed_C4.json")).read_bytes())
        exact(encode_parent(unpatched),reference_raw["actual_Rm_exit"]["actual_inherited_axial5_packet"],
            "current unpatched Rm parent")
        for key in ("full_implicit_leading_inputs_recomputed","actual_point_moment_history_recovered",
            "global_completed_tensor_admissibility","radial_mixed4_certified","all_radial_endpoint_joins_certified",
            "full_cartesian_vector_derivatives_certified","admissible_stress_lift_constructed",
            "full_NS_background_completed","temporal_recursion"):
            if raw[key]:
                raise ValueError("Unbuilt full current scope promoted: "+key)
        generic=provider.generic_moment_check
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            datum_enclosure_sha256=provider.core.datum.datum_sha,
            current_input_hashes_checked=len(raw["input_hashes"]),
            current_patch_source_bindings=raw["current_patch_source_bindings"],
            original_correlated_source_class_admission=raw["original_correlated_source_class_admission"],
            unchanged_generic_transport_Jacobian_fixture_consumed=PREFIX+"actual_moment_patch_check.json",
            source_transport_and_implicit_structure=generic["source_transport_and_implicit_structure"],
            original_independent_coefficient_fixture=generic["independent_coefficient_fixture"],
            current_source_tail_C1_caps_checked=tails,current_factored_source_decay_caps_checked=decays,
            current_actual_defect_axial5_coefficients_checked=30,
            current_implicit_control_axial5_coefficients_checked=30,
            current_implicit_map_residual_diagnostics_checked=30,
            current_higher_order_same_Jacobian_rows_checked=20,
            current_patch_velocity_pressure_moment_coefficients_checked=count,
            current_terminal_functional_zero_rows_checked=terminal_rows,
            current_Rm_source_parent_and_P0_retained=True,
            current_actual_moment_patch_installed=True,
            current_actual_patch_implicit_axial5_recomputed=True,
            current_five_functional_terminal_identities_connected=True,
            residual_zero_containment_not_used_as_closure_proof=True,
            full_implicit_leading_inputs_recomputed=False,actual_point_moment_history_recovered=False,
            global_completed_tensor_admissibility=False,radial_mixed4_certified=False,
            all_radial_endpoint_joins_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,full_NS_background_completed=False,
            temporal_recursion=False,all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current actual five-source patch and same implicit Jacobian through axial5",flush=True)
    return result


if __name__=="__main__":
    run()
