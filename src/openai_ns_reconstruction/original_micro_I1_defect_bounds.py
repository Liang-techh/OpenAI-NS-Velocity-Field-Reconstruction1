"""Same-source six I1 C1 defects and joined physical pressure error bounds."""
from __future__ import annotations

from pathlib import Path
import math

from .original_core_signed_direction import require,sha,read,record,source_guard
from .original_micro_direction_C0 import load,PREFIX,BRIDGE
from .original_log_mixed_jets import _logsum
from .original_ordinary_jet_bounds import jet_mul,jet_div
from .original_micro_Q_increment_bounds import recover

PARENT=PREFIX+"compliant_current_original_micro_Q_increment_bounds.json"
PARENT_SHA256="9c41dd5f3bf0dccb3e65582bc2452ede22bfe50ce223dadc57c6eccdf40e823c"
DEFECT="src/openai_ns_reconstruction/original_micro_finiteN_volume_defect.py"
CORRELATED="src/openai_ns_reconstruction/original_n1_micro_correlated_m15.py"
LOOP=PREFIX+"compliant_current_original_prescribed_loop_C5_functions.py"
NAMES=("q","q_Z","q_squared","U_q_Z","U_squared","U_U_Z")


def pack(values):
    return [dict(ordinary_Z_order=k,log_absolute_upper_bound=record(v)) for k,v in enumerate(values)]


def saved(c,rows):
    require([row["ordinary_Z_order"] for row in rows]==list(range(len(rows))),"Ordinary derivative rows required")
    return [read(c,row["log_absolute_upper_bound"]) for row in rows]


def add(c,*jets):
    require(len({len(jet) for jet in jets})==1,"Equal ordinary jet lengths required")
    return [_logsum(c,[jet[k] for jet in jets if jet[k] is not None]) for k in range(len(jets[0]))]


def mul(c,left,right):
    require(len(left)==len(right),"Equal ordinary jet lengths required")
    return [_logsum(c,[left[j]+right[k-j]+c.ln(math.comb(k,j))
                      for j in range(k+1) if left[j] is not None and right[k-j] is not None])
            for k in range(len(left))]


def scale(jet,logfactor):
    return [None if v is None else v+logfactor for v in jet]


def coefficient_logs(c,jet):
    return [None if c.fabs(v).b==0 else c.ln(c.fabs(v).b) for v in jet]


def coefficients(c,Z,de):
    """Actual signed C1 coefficient functions, differentiated before norms."""
    one=[c.mpf(1),c.mpf(0)];z=[Z,c.mpf(1)];h=[1-de,c.mpf(0)]
    z2=jet_mul(c,z,z);d=[one[k]-z2[k] for k in range(2)]
    L=[one[k]-de*z2[k] for k in range(2)]
    require(L[0].a>0,"Actual L floor required")
    L2=jet_mul(c,L,L)
    times=lambda jet,value:[value*v for v in jet]
    volume=dict(q=times(jet_div(c,one,L),-1),q_Z=times(jet_div(c,jet_mul(c,h,z),L),c.mpf(1)/2),
        q_squared=[-c.mpf(1)/2,c.mpf(0)],U_q_Z=jet_div(c,d,L),
        U_squared=times(jet_div(c,jet_mul(c,h,z2),L2),2),
        U_U_Z=times(jet_div(c,jet_mul(c,z,d),L2),2))
    boundary=dict(q=jet_div(c,one,L),q_squared=[c.mpf(1)/2,c.mpf(0)],
                  U_squared=times(jet_div(c,z2,L2),-2))
    return volume,boundary


def six_defects(c,U,Q,dU,dQ):
    require(len(U)==len(Q)==len(dU)==len(dQ)==3,"C1 product defects require source rows through order two")
    u,q,du,dq=U[:2],Q[:2],dU[:2],dQ[:2]
    uz,qz,duz,dqz=U[1:3],Q[1:3],dU[1:3],dQ[1:3]
    return dict(q=dq,q_Z=dqz,
        q_squared=add(c,scale(mul(c,q,dq),c.ln(2)),mul(c,dq,dq)),
        U_q_Z=add(c,mul(c,u,dqz),mul(c,du,qz),mul(c,du,dqz)),
        U_squared=add(c,scale(mul(c,u,du),c.ln(2)),mul(c,du,du)),
        U_U_Z=add(c,mul(c,u,duz),mul(c,du,uz),mul(c,du,duz)))


def background(c,parent):
    """New background-Q norm arithmetic on saved native U/M intervals."""
    source=parent["actual_micro_b_p2_mixed_bounds"]
    native=source["native_angular_and_axial_provider"]
    require(source["passed"] and native["passed"] and native["ordinary_Z0_Z4"]
            and native["physical_y0_y1"] and native["original_core_U_inlet_and_all_corrections_retained"]
            and native["native_and_comparison_never_identified"] and native["physical_dy_hb_ds_once"]
            and source["complete_inlet_integral_bound"]["same_complete_inlets_passed_to_native_and_comparison"],
            "Complete original native U/M input bounds required")
    def physical(rows,count):
        result=[]
        for row in rows[:count]:
            value=c.fabs(read(c,row)).b
            require(value>0,"Nonzero saved native absolute norm required")
            result.append(c.ln(value))
        require(len(result)==count,"Missing saved native derivative")
        return result
    U=physical(native["native_U"],3)
    M=physical(source["native_complete_history_rows"]["M"],4)
    de=read(c,parent["actual_core_signed_direction"]["parameters"]["delta"])
    Z=c.mpf([c.mpf(37)/100,c.mpf(93)/250])
    Q,proof=recover(c,U,M,Z,de,2)
    return U,Q,dict(passed=True,physical_native_U=pack(U),physical_native_mean=pack(M),
        physical_native_Q=pack(Q),Q_recovery=proof,
        inputs_are_saved_actual_native_bounds=True,saved_U_M_already_physical=True,
        inverse_Pstar_in_C6_restored_once_by_Q4=True,no_extra_Pstar_on_saved_U_M=True,
        native_mean_is_M_not_R_times_M=True,native_not_identified_with_comparison=True,
        new_background_Q_norm_computation=True,accepted_U_M_producer_not_replayed=True)


def weight_provenance(root,parent):
    path=LOOP[:-3]+".json.gz";expected=parent["input_hashes"][Path(path).name]
    require(sha(root/path)==expected,"Actual I1 physical-measure archive changed")
    graph=load(root/path);nodes=graph["exact_graph_nodes"]
    require(parent["actual_micro_phase_mean"]["graph_patch"]["kernel"]["graph_source_ports"]["scalars"]["hb"]==8
            and nodes[8]["name"]=="hbB" and nodes[8]["Z_independent"]
            and nodes[8]["quantity_is_exact_original_parameter_not_range_endpoint"],
            "Actual source hb handle is required, not the numeric value eight")
    ledger=nodes[8]["current_source_receipt"]
    require(parent["input_hashes"][ledger]==nodes[8]["current_source_receipt_sha256"]
            and sha(root/"experiments/root_st073"/ledger)==nodes[8]["current_source_receipt_sha256"],
            "Actual hb source ledger receipt must remain pinned")
    a=parent["actual_micro_mean_Z3_bounds"]["physical_scale_logs"]["hb"]
    basis=parent["actual_signed_lower_m15"]
    require(basis["source_log_basis_names"][0]=="hb"
            and all(a[key]==basis["source_log_bases"][0][key] for key in ("lower_exact_mpf_tuple","upper_exact_mpf_tuple")),
            "Integration loghb must be the accepted physical-measure logarithm")
    for name in ("first_micro","second_micro"):
        chart=graph["prescribed_loop_C5_bindings"]["charts"][name]
        require(chart["physical_measure"]==8 and chart["radius"]==68644,"Same actual R² hb I1 weight required")
    charts=graph["prescribed_loop_C5_bindings"]["charts"]
    first,second=charts["first_micro"],charts["second_micro"]
    lower=nodes[first["domain"][0]]
    require(lower["operation"]=="product" and len(lower["arguments"])==2 and 10 in lower["arguments"],
            "Original first lower must be sc/2")
    half=next(k for k in lower["arguments"] if k!=10)
    require(nodes[half]==dict(operation="exact_rational",numerator=1,denominator=2)
            and first["domain"][1]==second["domain"][0]==1 and second["domain"][1]==48
            and nodes[10]["Z_independent"]
            and parent["actual_micro_q_mixed_bounds"]["natural_collar_flat_proof"]["sc_positive_less_than_2_power_minus10"],
            "Original positive sc and chart widths <=1 required")
    require(nodes[48]==dict(operation="exact_rational",numerator=2,denominator=1),"Source endpoint handle 48 must denote two")
    return dict(passed=True,physical_measure_handle=8,physical_measure_is_original_hb_not_integer_eight=True,
        loghb_exact_tuples_equal_accepted_source_basis=True,same_actual_radius_handle=68644,
        native_chart_widths_at_most_one=True,actual_second_upper_number=2,
        source_first_lower_is_sc_half=True,positive_sc_upper_receipt_consumed=True,
        Z_independent_measure_and_native_endpoints=True,source_ledger_sha256=nodes[8]["current_source_receipt_sha256"],
        no_geometry_owner_or_hb_producer_executed=True),{path:dict(path=path,expected_sha256=expected)}


def source_contract(root,parent):
    source=parent["actual_micro_correlated_m15"]
    require(source["passed"] and source["actual_full_Q4_six_I1_source_implemented"]
            and source["physical_I1_measure"]=="R^2*hb*ds exactly once"
            and source["internal_RQy_and_reduced_B_seams_cancel_as_functions"],"Exact same-source joined I1 representation required")
    require(source["seam_certificate"]["source_graph_prefix_and_all_14_field_ports_checked"]
            and source["seam_certificate"]["free_Z_not_numeric_Z_substitution"]
            and not source["seam_certificate"]["accepted_fingerprint_computation_replayed"],
            "Accepted free-Z same-source endpoint identity required")
    mean=parent["actual_micro_phase_mean"]
    require(mean["source_proof"]["source_owned_first_collar_makes_Cin_zero"]
            and parent["actual_micro_phase_inverse_bounds"]["flat_branch"]["q_A_B_Aphi_Bphi_exact_zero"],
            "Source-owned flat lower endpoint required for the zero finite-N error")
    require(parent["actual_micro_Q_increment_bounds"]["source_contract"]
            ["full_recovered_Q_equals_original_by_same_U_and_mean_function_identity"],"Source Q function identity required")
    handles={chart:{name:source["graph_patch"]["finite_N_signed_integrand_defects"][chart][name] for name in NAMES}
             for chart in ("first_micro","second_micro")}
    return dict(passed=True,six_signed_defect_handles=handles,
        original_physical_I1_measure="R^2*hb*ds exactly once",
        accepted_free_Z_internal_B_and_RQy_seam_consumed=True,
        first_lower_finite_N_error_exact_zero_by_source_flat_collar=True,
        same_native_background_on_both_micro_charts=True,
        outer_RQy_retained=True,old_212_and_all_errors_unchanged=True,
        source_graph_functions_not_cap_values=True,accepted_seam_checker_replayed=False),dict(
        defects=source_guard(root,DEFECT,"build_defect",(
            "du = _jet_scale(g, delta_V[:3], Pstar)","dm = _jet_scale(g, history_m[:4], Pstar)",
            "U0 = _jet_sub(g, physical_U[:3], du)","Q0 = _jet_sub(g, physical_Q[:3], dq)")),
        native_scale=source_guard(root,BRIDGE+".py","generic_fields",(
            "inverse=g.quotient(g.one,Pstar,'same strictly positive Pstar')",
            "V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse)")),
        shared_background=source_guard(root,LOOP,"__init__",(
            "'first_micro':bridge.micro,'second_micro':bridge.micro",
            "fields={key:a.trim(jet) for key,jet in chart['fields'].items()}",
            "for key in RATES:fields[key]=a.add(a.trim(chart['fields'][key]),history[key])")),
        joined=source_guard(root,CORRELATED,"build",(
            'radius=chart["radius"];weight=fixed(g.mul(radius,radius,chart["physical_measure"]))',
            'S1=add(B,D);increment=add(half(S1),neg(boundary("R_Q_y")))',
            'neg(half(integrals["q_squared"])),divide(mul(d,integrals["U_q_Z"]))',
            'g.pscale(divide(mul(h,z,z,integrals["U_squared"]),L2),two)',
            'g.pscale(divide(mul(z,d,integrals["U_U_Z"]),L2),two)')))


def new_controls():
    import sympy as s
    z,de=s.symbols("z de");U,Q,dU,dQ=[s.Function(name)(z) for name in ("U","Q","dU","dQ")]
    terms=lambda u,q:dict(q=q,q_Z=s.diff(q,z),q_squared=q*q,U_q_Z=u*s.diff(q,z),U_squared=u*u,U_U_Z=u*s.diff(u,z))
    base,full=terms(U,Q),terms(U+dU,Q+dQ)
    expansion=dict(q=dQ,q_Z=s.diff(dQ,z),q_squared=2*Q*dQ+dQ*dQ,
        U_q_Z=U*s.diff(dQ,z)+dU*s.diff(Q,z)+dU*s.diff(dQ,z),
        U_squared=2*U*dU+dU*dU,U_U_Z=U*s.diff(dU,z)+dU*s.diff(U,z)+dU*s.diff(dU,z))
    require(all(s.expand(s.diff(full[key]-base[key]-expansion[key],z,k))==0 for key in NAMES for k in (0,1)),
            "Exact six C1 defect expansion failed")
    return dict(passed=True,six_polynomial_differences_and_ordinary_Z1=True,
        square_and_cross_terms_retained=True,ordinary_product_binomials=True,
        no_extra_Z_on_U_QZ_coefficient=True,U_squared_coefficient="2*(1-delta)*Z^2/L^2",
        U_UZ_coefficient="2*Z*(1-Z^2)/L^2",boundary_U_squared_coefficient="-2*Z^2/L^2",
        actual_native_upper_is_two_not_graph_handle_48=True,
        Pstar_in_physical_fields_not_added_to_I1_weight=True)


def build(root,parent):
    root=Path(root)
    require(sha(root/PARENT)==PARENT_SHA256 and load(root/PARENT)==parent,"Accepted Q-increment parent required")
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256)}
    for path in (DEFECT,CORRELATED,BRIDGE+".py",LOOP,"src/openai_ns_reconstruction/original_micro_Q_increment_bounds.py"):
        expected=parent["input_hashes"][Path(path).name] if path.startswith("experiments/") else parent["implementation_hashes"][path]
        require(sha(root/path)==expected,"Same-source I1 input changed")
        pins[path]=dict(path=path,expected_sha256=expected)
    contract,ast=source_contract(root,parent)
    weight,weight_pins=weight_provenance(root,parent);pins.update(weight_pins)
    from mpmath.ctx_iv import MPIntervalContext
    c=MPIntervalContext();c.prec=600
    U,Q,bg=background(c,parent)
    parameter=parent["actual_micro_mean_Z3_bounds"]["physical_scale_logs"]
    C,loghb=[read(c,parameter[key]) for key in ("selected_logC","hb")]
    inc=parent["actual_micro_Q_increment_bounds"]
    radius=inc["radius_provenance"]
    require(inc["passed"] and radius["passed"] and radius["physical_R_Z_independent"]
            and radius["R_equals_epsilon_rho"] and radius["epsilon_Lambda_equals_one"]
            and radius["actual_micro_radius_handle"]==68644,"Same actual physical R cap required")
    logR=read(c,inc["physical_R_upper_log"])
    de=read(c,inc["actual_delta"]);Z=c.mpf([c.mpf(37)/100,c.mpf(93)/250])
    volume,boundary=coefficients(c,Z,de)
    defects={};local={}
    for chart in ("first","second"):
        source=inc["charts"][chart]["log_caps"]
        dU=saved(c,source["physical_delta_U"])[:3];dQ=saved(c,source["physical_delta_Q"])
        defects[chart]=six_defects(c,U,Q,dU,dQ)
        local[chart]={key:scale(value,loghb+2*logR) for key,value in defects[chart].items()}
    joint={key:add(c,local["first"][key],local["second"][key]) for key in NAMES}
    D=add(c,*[mul(c,coefficient_logs(c,volume[key]),joint[key]) for key in NAMES])
    B=scale(add(c,*[mul(c,coefficient_logs(c,value),defects["second"][key]) for key,value in boundary.items()]),2*logR)
    partial=scale(add(c,B,D),-c.ln(2))
    outer=saved(c,inc["charts"]["second"]["log_caps"]["physical_R_delta_Q_y"])
    complete=add(c,partial,outer)
    require(all(v.b<(-2600*C).a for chart in defects.values() for values in chart.values() for v in values),"Six C1 defects envelope failed")
    require(all(v.b<(-2700*C).a for values in joint.values() for v in values)
            and all(v.b<(-2700*C).a for v in D+partial),"Joined I1 pressure volume envelope failed")
    threshold=-1024*c.ln(2)
    require(all(v.b<threshold.a for values in joint.values() for v in values)
            and all(v.b<threshold.a for v in B+D+partial),"Local I1 pressure error radius failed")
    fields=dict(reduced_volume_D=D,outer_reduced_boundary_B=B,half_S1_pressure_without_RQy=partial,
                retained_outer_RQy=outer,complete_micro_pressure_error=complete)
    return dict(passed=True,source_family=parent["source_family"],source_pins=pins,source_contract=contract,source_AST=ast,
        domain=inc["domain"],background_native_norm_provider=bg,weight_provenance=weight,new_six_defect_controls=new_controls(),
        coefficient_ordinary_Z0_Z1=dict(volume={key:[record(v) for v in values] for key,values in volume.items()},
            boundary={key:[record(v) for v in values] for key,values in boundary.items()}),
        charts={chart:dict(six_defect_log_caps={key:pack(value) for key,value in defects[chart].items()},
                          six_local_I1_log_caps={key:pack(value) for key,value in local[chart].items()}) for chart in defects},
        six_joint_I1_log_caps={key:pack(value) for key,value in joint.items()},
        pressure_error_log_caps={key:pack(value) for key,value in fields.items()},
        pressure_error_log_cap_ratios={key:[dict(ordinary_Z_order=k,log_cap_over_selected_logC=record(v/C)) for k,v in enumerate(value)] for key,value in fields.items()},
        strict_log_envelopes=dict(six_defects=-2600,joint_I1_and_reduced_volume_and_half_S1=-2700,reference="selected positive logC",passed=True),
        physical_I1_error_C1_below_2_power_minus1024=True,
        full_micro_pressure_error_small=[dict(ordinary_Z_order=k,below_2_power_minus1024=bool(v.b<threshold.a)) for k,v in enumerate(complete)],
        physical_R_squared_and_hb_applied_once=True,signed_function_values_evaluated=False,
        signed_six_integrals_evaluated=False,outer_RQy_certified_small=False,
        old_signed_sources_and_all_error_caps_unchanged=True,accepted_ancestors_replayed=False,
        global_five_moment_closure=False,full_goal_complete=False)
