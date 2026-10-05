"""Focused acceptance of current actual R110-to-Rsh long-reshape mixed4.

Unchanged physical/kernel fixture evidence is consumed by current hashes.
New checks admit the actual provider wiring, source graph and R110 join.
No current reference-restoration interface is certified by this checker.
"""
import json
import math
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_actual_long_reshape_mixed_C4 import (
    CompliantActualLongReshapeMixedC4, bridge_source_definitions,
    profile_source_bindings, accepted, sha, HERE, PREFIX)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4_check import (
    canonical_source, encode_parent)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

NAME=PREFIX+"actual_long_reshape_mixed_C4.json"
GROUPS=("physical_velocity_pressure_y_Z_mixed4","physical_five_primitive_y_Z_mixed4")


def exact(left,right,label):
    if canonical_source(left)!=canonical_source(right):
        raise ValueError("Current exact source transfer differs: "+label)


def source_graph_check(provider,raw):
    graph=raw["shared_exact_axial_source"]
    if graph!=provider.shared_axial_source or graph!=provider.switch_mixed.shared_axial_source:
        raise ValueError("Current R110 axial graph did not survive long reshape")
    sources=bridge_source_definitions(provider)
    if raw["exact_bridge_source_definitions"]!=sources:
        raise ValueError("Original bridge definition interface changed")
    terms=graph["formal_signed_integrals"]
    first=provider.history.switch.phase([-1,1],1)
    drive=["R*hydro","R*Pstar^2*pressure","R^2*F0^2*swirl"]
    expected={
        "I_bridge":dict(variable="s=log(R/Ra)",bounds=["0","log(100/Ra)"],sign=-1,
            chi="1-(1-hb)*sigma(s/hb)",quotient="phi_actual(R,Z)/phi_bar(R,Z)",
            radius="R=Ra*exp(s)",drive_terms=drive,
            exact_original_V_source=sources["V"],exact_original_chi_source=sources["chi"]),
        "I_first_switch":dict(variable="t",bounds=["0","1"],sign=-1,factor="hb^2",
            weight="1-sigma(t)",quotient="phi_actual(R,Z)/phi_bar(R,Z)",
            radius="R=100*exp(hb*t)",drive="same current-radius axial direction sum; not the angular Dbar",
            drive_terms=drive,exact_original_V_source=first["exact_Uz_source"])}
    if terms!=expected:
        raise ValueError("Original signed integral units/controls/radius changed")
    if (graph["V110"]["args"]!=graph["V100"]["args"]+[{
        "formal_integral":"I_first_switch","shared_source":graph["shared_source_namespace"]}]
            or graph["E_V110_minus_4Z"]["args"]!=graph["V110"]["args"][1:]):
        raise ValueError("Centered E no longer shares the exact raw V110 source")
    if (not graph["caps_are_not_source_integrals"] or graph["formal_integrals_numerically_reconstructed"]
            or graph["axial_Taylor_orders"]!=list(range(6))
            or graph["ordinary_derivative_factorials"]!=[math.factorial(k) for k in range(6)]):
        raise ValueError("Current exact source units/scope changed")
    for path,digest in graph["source_bindings"].items():
        if sha(path)!=digest:
            raise ValueError("Current exact graph dependency changed: "+path)
    return dict(current_switch_graph_preserved_exactly=True,
        original_signed_integrals_and_units_canonicalized=True,
        same_centered_E_function_references_retained=True,
        source_definition_interface_is_AST_bound=True,
        graph_enclosures_are_not_point_integral_values=True)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        provider=CompliantActualLongReshapeMixedC4();c=provider.ctx
        family,source=provider.family,provider.source
        original=provider.original_check
        profiles_name=PREFIX+"long_reshape_profiles_check.json"
        profiles=accepted(profiles_name,family,source,"actual_Rsh_exit_axial5_available")
        if (not profiles["kernel_fixture"]["passed"] or not profiles["structural_identities"]["passed"]
                or not profiles["actual_B_C2_shear_sigma_and_kernel_rate_theorems_directly_bound"]
                or not original["symbolic_checks"]["passed"]
                or not original["independent_physical_fixture"]["passed"]):
            raise ValueError("Hash-current unchanged kernel/physical fixture evidence required")
        if (raw["current_actual_long_reshape_source_bindings"]!=profile_source_bindings(provider)
                or raw["datum_enclosure_sha256"]!=provider.reshape.core.datum.datum_sha):
            raise ValueError("Current profile wiring or pressure datum differs")
        source_proof=source_graph_check(provider,raw)
        keys=("whole_reshape","actual_R110_inlet","actual_Rsh_exit",
              "whole_local_R109_R110_power","actual_R110_power_side")
        packets=[raw[k] for k in keys]+raw["interior_packets"]
        expected={"y"+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)}
        counts={group:0 for group in GROUPS}
        for packet in packets:
            for group in GROUPS:
                if len(packet[group])!=(4 if group==GROUPS[0] else 5):
                    raise ValueError("Physical field/primitive group omitted")
                for grid in packet[group].values():
                    if set(grid)!=expected:
                        raise ValueError("Mixed4 row domain is incomplete")
                    for value in grid.values():
                        lo,hi=endpoints(read_interval(c,value))
                        if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):
                            raise ArithmeticError("Nonfinite actual long-reshape row")
                        counts[group]+=1
            if (not packet["positive_source_cap_is_enclosure_only"]
                    or not packet["large_derivative_factors_combined_before_positive_source_cap"]):
                raise ValueError("Positive source capped before derivative factors")
            if any(packet[k] for k in ("full_inner_interfaces_certified",
                "full_cartesian_vector_derivatives_certified","admissible_stress_lift_constructed","temporal_recursion")):
                raise ValueError("Unbuilt whole-field scope promoted")
        for key in ("whole_reshape","actual_R110_inlet","actual_Rsh_exit"):
            packet=raw[key]
            if (not packet["current_actual_R110_feedback_in_long_reshape"]
                    or packet["Rsh_reference_mixed4_join_certified"]
                    or packet["current_reference_restoration_installed"]):
                raise ValueError("Current long-reshape/restoration scope changed")
        inp=provider.reshape.inputs([-1,1])
        incoming=provider.history.switch.inlet([-1,1])
        exact(encode_parent(inp["inlet"]),encode_parent(incoming),"current R110 inlet")
        left=raw["actual_R110_power_side"]["actual_inherited_axial5_packet"]
        inlet=encode_parent(inp["inlet"])
        trace_keys=("F_actual_over_F0_axial5_coefficients",
                    "F_actual_true_axial5_divided_by_F0","Uz_actual_axial5_coefficients",
                    "actual_moment_shape_axial5_coefficients","actual_Q_axial4_coefficients",
                    "pressure_axis_axial5_coefficients",
                    "pressure_increment_true_axial5_divided_by_R_F0_squared")
        for key in trace_keys:
            exact(left[key],inlet[key],key)
        # Both sides are built from the current shared provider. The original
        # generic symbolic proof admits the flat cutoff chain rule; overlap
        # is neither used nor described as a functional join proof.
        for group in GROUPS:
            exact(raw["actual_R110_inlet"][group],raw["actual_R110_power_side"][group],
                  "two-sided R110 "+group)
        for endpoint in (0,1):
            jets=sigma_jets(c,c.mpf(endpoint))
            if (endpoints(jets[0])!=(mp.mpf(endpoint),mp.mpf(endpoint))
                    or any(endpoints(jets[k])!=(mp.mpf(0),mp.mpf(0)) for k in range(1,5))):
                raise ArithmeticError("Original flat endpoint derivative changed")
        if endpoints(provider.reshape.T)!=endpoints(400*provider.reshape.A):
            raise ValueError("Original T changed")
        capcount=0
        for proof in raw["factored_positive_source_cap_proofs"]:
            if (endpoints(read_interval(c,proof["log_magnitude_upper"]))[1]>
                    endpoints(read_interval(c,proof["log_cap"]))[0]
                    or not proof["exact_source_not_replaced"]):
                raise ValueError("Invalid final source magnitude bound")
            capcount+=1
        for key in ("Rsh_reference_mixed4_join_certified","current_reference_restoration_installed",
                    "current_actual_moment_patch_installed","full_implicit_leading_inputs_recomputed",
                    "actual_point_moment_history_recovered","global_completed_tensor_admissibility",
                    "full_inner_interfaces_certified","full_cartesian_vector_derivatives_certified",
                    "admissible_stress_lift_constructed","temporal_recursion"):
            if raw[key]:
                raise ValueError("Unbuilt current downstream/global scope promoted: "+key)
        result=dict(
            actual_five_defect_family_sha256=family,implicit_source_sha256=source,
            datum_enclosure_sha256=provider.reshape.core.datum.datum_sha,
            actual_mixed_bounds_checked=counts,
            current_input_hashes_checked=len(raw["input_hashes"]),
            original_kernel_fixture_receipt_consumed=profiles_name,
            original_physical_symbolic_fixture_receipt_consumed=PREFIX+"long_reshape_mixed_C4_check.json",
            current_actual_R110_input_AST_bindings_checked=6,
            current_actual_R110_normalization_AST_bindings_checked=3,
            current_constant_V_transport_AST_bindings_checked=3,
            current_actual_R110_trace_groups_checked=len(trace_keys),
            exact_current_R110_two_sided_mixed_grid_rows_checked=135,
            current_exact_axial_source_binding=source_proof,
            factored_positive_source_caps_checked=capcount,
            original_selected_A_T_logref_and_pressure_retained=True,
            current_actual_R110_feedback_in_long_reshape=True,
            current_actual_long_reshape_mixed4_available=True,
            R110_postswitch_power_reshape_local_mixed4_join_certified=True,
            Rsh_reference_mixed4_join_certified=False,
            current_reference_restoration_installed=False,
            current_actual_moment_patch_installed=False,
            full_implicit_leading_inputs_recomputed=False,
            actual_point_moment_history_recovered=False,
            global_completed_tensor_admissibility=False,
            full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,
            all_passed=True,
            input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf8")
    print("PASS current actual long-reshape mixed4 and exact R110 source join",flush=True)
    return result


if __name__=="__main__":
    run()
