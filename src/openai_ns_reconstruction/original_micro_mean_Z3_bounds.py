"""Source-owned C/Cy and complete physical micro mean Z0..Z3 norm bounds.

The accepted signed graphs remain the functions.  This module attaches
directed log-absolute caps to their ordinary derivative rows, including both
moving histories and the inherited memory on the second original chart.
"""
from __future__ import annotations

from pathlib import Path

from .original_core_signed_direction import require, sha, read, record, source_guard
from .original_micro_direction_C0 import load, PREFIX
from .original_log_mixed_jets import _logsum

PARENT = PREFIX + "compliant_current_original_micro_phase_inverse_bounds.json"
PARENT_SHA256 = "0551d3b8acbcc4578cf97591a299e679f3e589efd4ce6d8e8e66b765e2cbabdc"
KERNEL = "src/openai_ns_reconstruction/original_micro_phase_mean_kernel.py"
KERNEL_SHA256 = "3c22ebb63788ab1b8dd50dd8fa5923a67f579ec464ef00d6a9ec81eea7de1181"


def rows(jet):
    return [dict(ordinary_Z_order=k,log_absolute_upper_bound=record(value)) for k,value in enumerate(jet)]


def saved_B_rows(c, chart, y_order):
    source=chart["mixed_log_caps"]["B_over_Pstar_over_frozen_E0"]
    out={row["ordinary_Z_order"]:read(c,row["log_absolute_upper_bound"])
         for row in source if row["physical_y_order"]==y_order and not row["exact_zero"]}
    require(set(out)==set(range(4)),"Actual ordinary B/slow-By rows missing")
    require(chart["E_product_rule_retained"] and chart["E_reference_is_point_frozen_not_function_division"],
            "Actual E product rule and frozen reference required")
    return [out[k] for k in range(4)]


def source_contract(root,parent):
    old=parent["actual_micro_phase_mean"]; proof=old["source_proof"]
    require(old["passed"] and old["ordinary_mean_Z3_source_implemented"]
            and old["complete_C_and_Cy_convolutions_retained"]
            and old["same_source_internal_C_endpoint_canceled"],"Complete source mean packet required")
    for key in ("C_phi_equals_B","C_0_equals_C_1_zero","C_y_is_fixed_phase_primitive_of_slow_B_y",
                "source_owned_first_collar_makes_Cin_zero","only_original_rminus_mean_inlet_zero",
                "first_second_C_and_Cy_seam_equal_as_free_Z_functions",
                "ordinary_Z3_from_same_C5_functions_not_differentiated_caps",
                "phase_total_y_derivative_equals_N","phase_integer_cuts_use_periodic_extension",
                "actual_huge_N_phase_not_chosen_or_materialized"):
        require(proof[key],"Missing accepted source identity: "+key)
    algebra=old["exact_algebra"]
    require(algebra["passed"] and algebra["transport_delta_m_y_plus_delta_m_equals_Pstar_B_over_N"]
            and all(algebra["seam_cancel_through_Z0_to_Z3"].values()),"Source mean transport/seam identity required")
    ports=old["graph_patch"]["kernel"]["graph_source_ports"]
    require(ports["second_memory_formula"]=="exp(-hb*(s-1))"
            and ports["first_Cin_source_owned_zero_required"],"Original joined history memory required")
    for chart in ("first_micro","second_micro"):
        port=ports["charts"][chart]
        require(port["actual_phase_Z_independent"] and port["fixed_phase_slow_B_y"]
                and port["total_y_chain_excluded"] and port["ordinary_Z_orders"]==[0,1,2,3]
                and not port["factorial_normalization"],"Actual source primitive derivative port differs")
    require(parent["actual_micro_q_mixed_bounds"]["direct_source_provenance"]
            ["original_sc_is_Z_independent_exact_parameter"],"Fixed original chart endpoints required")
    kernel=old["graph_patch"]["kernel"]
    fields=("C","Cy","JC","JCy","JC_joint","JCy_joint","complete_physical_delta_m")
    handles={key:kernel[key] for key in fields}
    for key in fields:
        require(all(len(value)==4 for value in handles[key].values()),"Source graph Z3 rows missing: "+key)
    return dict(passed=True,accepted_source_proof_consumed=True,
        exact_transport_identity_consumed=True,ordinary_Z_seam_identity_consumed=True,
        source_graph_handles=handles,second_memory_formula=ports["second_memory_formula"],
        first_inlet_zero_only_at_source_owned_left_collar=True,second_history_not_reset=True,
        actual_phase_Z_independent=True,ordinary_Z_derivatives_have_no_factorials=True,
        source_phase_values_not_chosen=True,source_or_checker_not_executed=True),dict(
        phase=source_guard(root,KERNEL,"_phase_rows",(
            '"fixed_phase_slow_y": source_name == "Cy"',
            '"total_y_chain_excluded": source_name == "Cy"',
            '"phase_primitive_is_direct_ordinary_Z_row": True',
            '"source_operation": "phase_integral_of_slow_B_y" if source_name == "Cy" else "phase_integral_of_B"')),
        moving=source_guard(root,KERNEL,"_moving_rows",(
            'exponent = graph.neg(graph.mul(hb, distance))',
            'integrand = graph.mul(kernel, primitive_at_t, hb)',
            '"physical_hb_weight_applied_once": True')),
        joined=source_guard(root,KERNEL,"build_micro_mean",(
            'scale = g.div(Pstar, n_squared, "same exact positive common N squared")',
            'formula="exp(-hb*(s-1))"',
            '_scale_rows(g, first_JC_at_one, second_memory)',
            '_scale_rows(g, first_JCy_at_one, second_memory)',
            'second_complete = _add_rows(g, second_known, second_unresolved)')))


def history_controls():
    """New exact controls for joined convolution and physical transport."""
    import sympy as s
    h,u,v=s.symbols("h u v",positive=True)
    mass=s.integrate(h*s.exp(-h*(u-v)),(v,0,u))
    require(s.simplify(mass-(1-s.exp(-h*u)))==0,"Moving mass identity failed")
    # Source first width <=1 and second width <=1.  The exponential is <=1;
    # hence first mass <=hb and joined mass <=2hb, including incoming memory.
    N,P,C,Cy,B,JC,JCy=s.symbols("N P C Cy B JC JCy",nonzero=True)
    dm=P/N**2*(C-JC-JCy)
    total_y=P/N**2*(Cy+N*B-(C-JC)-(Cy-JCy))
    require(s.simplify(total_y+dm-P*B/N)==0,"Complete physical transport identity failed")
    require(s.simplify(s.exp(-h*(u-1))*s.exp(-h*(1-v))-s.exp(-h*(u-v)))==0,
            "Second history is not the joined kernel")
    return dict(passed=True,positive_kernel_mass="1-exp(-hb*length)",
        first_native_width_at_most_one=True,second_native_width_at_most_one=True,
        first_mass_upper="hb",joined_mass_upper="2*hb",physical_hb_applied_once=True,
        second_memory_times_first_kernel_equals_joined_kernel=True,
        mean_formula="delta_m=Pstar/N^2*(C-JC_joint-JCy_joint)",
        total_y_identity="delta_m_total_y+delta_m=Pstar*B_over_Pstar/N",
        norm_bound_does_not_drop_Cy_or_incoming_memory=True,
        ordinary_Z_commutes_with_phase_and_moving_integrals=True,
        reason="phase_Z=0; hb,N,sc and native chart endpoints independent of Z; saved smooth source jets",
        periodic_phase_cuts="C(0)=C(1)=0 and derivatives of this parameter identity; Cy=C_slow_y",
        threshold_is_local_mean_error_not_full_NS_residual=True)


def build(root,parent):
    root=Path(root)
    require(sha(root/PARENT)==PARENT_SHA256 and load(root/PARENT)==parent,"Accepted fixed-phase provider required")
    require(sha(root/KERNEL)==KERNEL_SHA256,"Source-owned mean kernel changed")
    phase=parent["actual_micro_phase_inverse_bounds"]
    require(phase["passed"] and phase["loop_AB_fixed_phase_mixed_provider_complete"],"Actual B and slow By source norms required")
    contract,ast=source_contract(root,parent)
    from mpmath.ctx_iv import MPIntervalContext
    c=MPIntervalContext();c.prec=600
    old=parent["actual_micro_direction_C0"]["actual_C0_errors"]
    require(old["passed"] and old["phase_free_uniform_bounds"] and old["F0_positive_and_at_most_one"],
            "Global positive actual E amplitude bound required")
    direction=parent["actual_micro_direction_C0"]["comparison_direction"]
    require(direction["passed"] and direction["domain"]==dict(s=["0","2"],Z=["37/100","93/250"])
            and direction["actual_native_phi_positive_and_nonincreasing"]
            and parent["actual_micro_direction_C0"]["both_micro_actual_a_and_E_positive"]
            and parent["actual_core_signed_direction"]["micro_core_argument_domain"]["passed"],
            "Saved global amplitude cap must cover the actual frozen E0 on both micro charts")
    logPE=read(c,old["log_Pstar_E_upper"])
    logs=parent["actual_micro_q_mixed_bounds"]["parameter_envelopes"]["actual_logs"]
    C,loghb,logN=[read(c,logs[key]) for key in ("C","hb","N")]
    require(C.a>0 and loghb.b<0 and logN.a>0,"Positive actual scales required")
    primitive={}
    for chart in ("first","second"):
        current=phase["charts"][chart]
        # Uniform |phase|<=1: each primitive row is bounded by its complete
        # source B or slow By row, with the saved positive Pstar*E floor-free cap.
        physical_B=[logPE+value for value in saved_B_rows(c,current,0)]
        physical_slow_By=[logPE+value for value in saved_B_rows(c,current,1)]
        primitive[chart]=dict(B=physical_B,slow_By=physical_slow_By,
                              C=list(physical_B),Cy=list(physical_slow_By))
    histories={}
    for chart in ("first","second"):
        histories[chart]={}
        for key in ("C","Cy"):
            first=[loghb+value for value in primitive["first"][key]]
            histories[chart][key]=first if chart=="first" else [
                _logsum(c,[first[k],loghb+primitive["second"][key][k]]) for k in range(4)]
    charts={}
    for chart in ("first","second"):
        pc,py=primitive[chart]["C"],primitive[chart]["Cy"]
        jc,jy=histories[chart]["C"],histories[chart]["Cy"]
        known=[_logsum(c,[pc[k],jc[k]])-2*logN for k in range(4)]
        unresolved=[jy[k]-2*logN for k in range(4)]
        mean=[_logsum(c,[pc[k],jc[k],jy[k]])-2*logN for k in range(4)]
        mean_y=[_logsum(c,[primitive[chart]["B"][k]-logN,mean[k]]) for k in range(4)]
        fields=dict(physical_Pstar_C=pc,physical_Pstar_Cy=py,
                    physical_Pstar_JC_joint=jc,physical_Pstar_JCy_joint=jy,
                    known_C_mean_part=known,retained_Cy_mean_part=unresolved,
                    physical_complete_delta_m=mean,physical_total_y_delta_m=mean_y)
        require(all(value.b<(-4500*C).a for value in mean),"Complete mean Z3 coarse log margin failed")
        require(all(value.b<(-600*C).a for value in mean_y),"Complete mean total-y Z3 coarse log margin failed")
        threshold=-1024*c.ln(2)
        require(all(value.b<threshold.a for value in mean+mean_y),"Local physical mean error radius failed")
        charts[chart]=dict(passed=True,ordinary_Z_orders=[0,1,2,3],
            source_physical_B_log_caps={key:rows(primitive[chart][key]) for key in ("B","slow_By")},
            primitive_caps_inherited_from_B_by_length_one_not_function_equality=True,
            log_caps={key:rows(value) for key,value in fields.items()},
            log_cap_ratios={key:[dict(ordinary_Z_order=k,log_cap_over_selected_logC=record(value/C))
                                for k,value in enumerate(values)] for key,values in fields.items()},
            ratios_are_logcaps_divided_by_selected_positive_logC=True,
            complete_physical_mean_Z0_Z3_below_2_power_minus1024=True,
            complete_physical_total_y_mean_Z0_Z3_below_2_power_minus1024=True,
            strict_log_envelope=dict(mean_multiplier=-4500,total_y_mean_multiplier=-600,reference="selected positive logC",passed=True),
            source_chart=chart+"_micro",C_and_Cy_not_reset=True,
            second_incoming_history_retained=(chart=="second"))
    return dict(passed=True,source_family=parent["source_family"],
        source_pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256),
                     KERNEL:dict(path=KERNEL,expected_sha256=KERNEL_SHA256)},
        source_contract=contract,source_AST=ast,new_history_controls=history_controls(),
        domain=phase["domain"],ordinary_Z_orders=[0,1,2,3],charts=charts,
        physical_scale_logs=dict(Pstar_E_upper=record(logPE),hb=record(loghb),N=record(logN),selected_logC=record(C)),
        amplitude_provenance=dict(passed=True,actual_positive_E0_is_the_same_native_E=True,
            frozen_E0_at_every_source_point_below_saved_uniform_E_upper=True,
            domain=dict(s=["0","2"],Z=["37/100","93/250"]),
            bound="Pstar*E=F0*sqrt(2R)*phi <= sqrt(2*rho_upper/Lambda)*Phi_core_upper",
            F0_positive_and_at_most_one=True,native_phi_positive_and_nonincreasing=True,
            actual_core_argument_domain_receipt_consumed=True,saved_global_log_bound_consumed=True),
        actual_Cy_Z3_norm_provider=True,actual_complete_mean_Z3_norm_provider=True,
        actual_complete_total_y_mean_Z3_norm_provider=True,
        primitive_phase_integration_length_bound_one=True,
        physical_Pstar_applied_once_and_common_N_squared_retained=True,
        finite_N_velocity_derivative_errors=False,numerical_QZZ=False,full_C5_numeric_containment=False,
        signed_mean_function_values_evaluated=False,zero_is_symmetric_error_box_center_not_selected_value=True,
        old_signed_sources_and_all_error_caps_unchanged=True,accepted_ancestors_replayed=False,
        global_five_moment_closure=False,full_goal_complete=False)
