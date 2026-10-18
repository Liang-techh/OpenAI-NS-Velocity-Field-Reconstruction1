"""Physical finite-N U/Q increment norms from the same complete mean.

Q is recovered by its linear source operator before taking absolute caps.
The total-y velocity chain includes B_phi without a further N factor.
Neither huge signed fields nor an unknown issued phase are selected.
"""
from __future__ import annotations

from pathlib import Path
import ast

from .original_core_signed_direction import require,sha,read,record,source_guard
from .original_micro_direction_C0 import load,PREFIX
from .original_log_mixed_jets import _logsum

PARENT=PREFIX+"compliant_current_original_micro_mean_Z3_bounds.json"
PARENT_SHA256="557195be5893c17e9398371017b3b4b7a52b6a0eae473b6befc19cafe84b8aef"
RECOVERY="src/openai_ns_reconstruction/original_n1_micro_phase_mean.py"
Q4=PREFIX+"compliant_current_original_modified_Q4_pressure_C3_functions.py"
LOOP=PREFIX+"compliant_current_original_prescribed_loop_C5_functions.py"
CORE_CHECK=PREFIX+"compliant_current_original_n1_core_fixed_point_functions_check.json"


def pack(values):
    return [dict(ordinary_Z_order=k,log_absolute_upper_bound=record(v)) for k,v in enumerate(values)]


def saved(c,rows):
    require([r["ordinary_Z_order"] for r in rows]==list(range(len(rows))),"Ordered ordinary derivative rows required")
    return [read(c,r["log_absolute_upper_bound"]) for r in rows]


def linear(c,terms):
    """Absolute linear combination; source signs are retained by the AST."""
    values=[]
    for coefficient,logcap in terms:
        absolute=c.fabs(coefficient).b
        if absolute>0 and logcap is not None:values.append(c.ln(absolute)+logcap)
    return _logsum(c,values)


def recover(c,U,m,Z,delta,order):
    """Ordinary Q derivative caps; no invented top derivative or padding."""
    require(order in (1,2) and len(U)==order+1 and len(m)==order+2,"Exact source derivative order required")
    h,d,L=1-delta,1-Z*Z,1-delta*Z*Z
    L1,L2=-2*delta*Z,-2*delta
    require(L.a>0,"Actual L=1-delta*Z^2 must be positive")
    n=[linear(c,[(2*Z,U[0]),(-h*Z,m[0]),(-d,m[1])]),
       linear(c,[(c.mpf(2),U[0]),(2*Z,U[1]),(-h,m[0]),((1+delta)*Z,m[1]),(-d,m[2])])]
    if order==2:n.append(linear(c,[(c.mpf(4),U[1]),(2*Z,U[2]),(2*delta,m[1]),((3+delta)*Z,m[2]),(-d,m[3])]))
    inv=-c.ln(L.a)
    q=[n[0]+inv]
    q.append(linear(c,[(c.mpf(1),n[1]),(-L1,q[0])])+inv)
    if order==2:q.append(linear(c,[(c.mpf(1),n[2]),(-2*L1,q[1]),(-L2,q[0])])+inv)
    return q,dict(passed=True,positive_L=record(L),L_Z=record(L1),L_ZZ=record(L2),
                  numerator_log_caps=pack(n),ordinary_derivatives=True,no_top_row_padding=True)


def new_controls():
    """Independent linearity and nonunit-denominator derivative controls."""
    import sympy as s
    z,de=s.symbols("z de")
    U,M,A,B=[s.Function(name)(z) for name in ("U","M","A","B")]
    op=lambda u,m:(2*z*u-(1-de)*z*m-(1-z*z)*s.diff(m,z))/(1-de*z*z)
    require(s.cancel(op(U+A,M+B)-op(U,M)-op(A,B))==0,"Q source linearity failed")
    uu=1-3*z+2*z*z;mm=2+z-z*z+3*z**3
    L=1-s.Rational(2,7)*z*z
    numerator=2*z*uu-(1-s.Rational(2,7))*z*mm-(1-z*z)*s.diff(mm,z)
    n=[numerator,s.diff(numerator,z),s.diff(numerator,z,2)]
    q=[n[0]/L]
    q.append((n[1]-s.diff(L,z)*q[0])/L)
    q.append((n[2]-2*s.diff(L,z)*q[1]-s.diff(L,z,2)*q[0])/L)
    require(all(s.cancel(value-s.diff(numerator/L,z,k))==0 for k,value in enumerate(q)),
            "Nonunit Q quotient recurrence differs from direct derivative")
    return dict(passed=True,operator_linearity_exact=True,
        nonunit_delta_control="2/7",ordinary_quotient_rows_checked=[0,1,2],
        ordinary_Z_factorials_not_inserted=True,
        Q_error_recovered_before_absolute_norm=True,
        mean_y_identity="delta_m_y=delta_U-delta_m",
        total_y_velocity="delta_U_y=Pstar*(slow_B_y/N+B_over_Pstar_phi)",
        no_extra_N_on_B_phi=True,physical_Pstar_applied_once=True,
        source_linear_error_not_subtraction_of_large_interval_boxes=True)


def provenance(root,parent):
    old=parent["actual_micro_phase_mean"]
    require(old["passed"] and old["exact_algebra"]["transport_delta_m_y_plus_delta_m_equals_Pstar_B_over_N"],
            "Actual complete mean source transport required")
    control=old["physical_Q_recovery_controls"]
    require(control["passed"] and all(control["ordinary_Z0_Z1_Z2"]) and control["mean_Z3_retained"]
            and control["no_factorial_normalization"],"Accepted same-source Q recurrence required")
    handles={}
    for chart in ("first_micro","second_micro"):
        item=old["graph_patch"]["Q_recovery"][chart]
        require(item["same_source_replacement_not_extra_additive_contribution"],"Q recovery cannot be an extra contribution")
        handles[chart]={key:item[key] for key in ("physical_mean_ordinary_Z0_Z3","physical_Q_ordinary_Z0_Z2","original_Q_source_nodes")}
    require(old["source_proof"]["phase_total_y_derivative_equals_N"]
            and old["source_proof"]["actual_huge_N_phase_not_chosen_or_materialized"],"Same unchanged actual phase chain required")
    require(old["source_proof"]["source_owned_first_collar_makes_Cin_zero"]
            and old["source_proof"]["first_second_C_and_Cy_seam_equal_as_free_Z_functions"]
            and old["exact_algebra"]["assumptions"]["second_chart_does_not_reset_first_memory"],
            "Original mean inlet and joined memory are required for function equality")
    tree=ast.parse((root/LOOP).read_text(encoding="utf-8"))
    assignments=[node for node in tree.body if isinstance(node,ast.Assign)
                 and any(isinstance(t,ast.Name) and t.id=="RATES" for t in node.targets)]
    require(len(assignments)==1,"Original density transport rates changed")
    rate=next(value.value for value in assignments[0].value.keywords if value.arg=="m")
    require(ast.dump(rate)==ast.dump(ast.parse("Fraction(1)",mode="eval").body),"Physical mean transport rate must be one")
    receipts=dict(recovery=source_guard(root,RECOVERY,"recover_Q",(
        "n0=sub(g.mul(two,z,U[0]),g.add(g.mul(h,z,m[0]),g.mul(d,m[1])))",
        "n2=g.add(g.mul(g.constant(4),U[1]),g.mul(two,z,U[2]),g.mul(two,de,m[1]),",
        "q2=g.div(sub(n2,g.add(g.mul(two,L1,q1),g.mul(L2,q0))),L,cert)")),
        physical=source_guard(root,Q4,"physical_pressure_fields",(
            "U=a.scale(fields['V'],Pstar)","my=a.sub(U,m)",
            "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))")),
        total_y=source_guard(root,Q4,"total_radial_functions",(
            "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",)),
        density=source_guard(root,LOOP,"finite_N_functions",(
            "s=a.scale(loop['B_over_Pstar'],invN)","inc=dict(m=s,",
            "V_N=a.add(V,s)")),
        history=source_guard(root,LOOP,"transport",(
            "for name,rate in RATES.items():", "predecessor_memory_never_reset=True",
            "out[name]=a.add(a.scale(incoming[name],decay),contribution)")),
        original_fields=source_guard(root,LOOP,"__init__",(
            "incoming={key:a.fixed(0) for key in RATES}",
            "history,parts,memory=transport(g,finite['five_signed_density_increments'],incoming,",
            "for key in RATES:fields[key]=a.add(a.trim(chart['fields'][key]),history[key])",
            "incoming=self.endpoints[name]")),
        recovered_fields=source_guard(root,RECOVERY,"build",(
            'mean=[g.add(g.mul(33572,bg[k]),kernel["complete_physical_delta_m"][name][k]) for k in range(4)]',
            "q=recover_Q(g,U,mean)")))
    return dict(passed=True,source_graph_handles=handles,
        accepted_Q_recovery_controls_consumed=True,actual_physical_U_and_mean=True,
        same_source_Q_replacement_not_extra_ledger_contribution=True,
        phase_Z_independent_and_total_y_phase_derivative_N=True,
        original_mean_transport_rate_one_AST_bound=True,
        mean_replacement_identity="same delta_m_y+delta_m=Pstar*B0/N, same zero left inlet and nonreset seam; linear ODE uniqueness",
        full_recovered_Q_equals_original_by_same_U_and_mean_function_identity=True,
        old_source_function_and_all_error_packets_retained=True,
        accepted_recovery_controls_or_source_producers_replayed=False),receipts


def radius_provenance(root,parent):
    """Read exact graph recipes; never execute the geometry/core owners."""
    paths=(CORE_CHECK,LOOP[:-3]+".json.gz",Q4[:-3]+".json.gz",Q4[:-3]+"_check.json")
    pins={}
    for path in paths:
        expected=parent["input_hashes"][Path(path).name]
        require(sha(root/path)==expected,"Actual radius source archive changed")
        pins[path]=dict(path=path,expected_sha256=expected)
    recipe=load(root/CORE_CHECK)["checks"]["actual_oracle"]["exact_parameter_and_native_P0_identity"]
    require(recipe["exact_graph_logP_delta_logLambda_epsilon_recipes_verified"],"Exact epsilon/Lambda source recipe required")
    source=parent["actual_core_signed_direction"]["source_AST"]["actual_micro_core_argument"]
    require(source["passed"] and source["source_executed"] is False
            and "rho=lambda q:g.mul(c(4),g.unary('exp',g.mul(hb,q)))" in source["required_text"],
            "Saved actual rho(s) source AST required")
    loop,q4=load(root/paths[1]),load(root/paths[2])
    ln,qn=loop["exact_graph_nodes"],q4["exact_graph_nodes"]
    expected={65604:dict(operation="product",arguments=[13,35565]),
              68642:dict(operation="product",arguments=[8,66988]),
              68643:dict(operation="analytic_unary",name="exp",argument=68642),
              68644:dict(operation="product",arguments=[65604,68643])}
    require(all(ln[k]==qn[k]==value for k,value in expected.items()),"Same physical R graph prefix required")
    require(ln[8]["Z_independent"] and ln[20]==dict(operation="bound_variable",name="Z")
            and ln[66988]==dict(operation="bound_variable",name="original_prescribed_C6_micro_s"),
            "Actual native s/hb must be independent of ordinary Z")
    require(ln[13]==dict(operation="exact_rational",numerator=4,denominator=1),"Actual Ra=4*epsilon required")
    require(ln[35565]==dict(operation="analytic_unary",name="exp",argument=35564)
            and ln[35564]==dict(operation="negative",argument=35562)
            and ln[35562]==dict(operation="sum",arguments=[14,643])
            and ln[14]==dict(operation="product",arguments=[13,5])
            and ln[643]==dict(operation="exact_rational",numerator=1000,denominator=1),
            "Exact epsilon=exp[-(4logPstar+1000)] graph recipe required")
    endpoints={}
    for name in ("first_micro","second_micro"):
        chart=loop["prescribed_loop_C5_bindings"]["charts"][name]
        require(chart["radius"]==68644 and chart["physical_measure"]==8,"Same actual micro radius/Jacobian required")
        physical=q4["modified_Q4_pressure_C3_bindings"]["charts"][name]
        boundary=physical["finite_boundary"]["R_Q_y"][:2]
        require(all(qn[ref]["operation"]=="product" and sorted(qn[ref]["arguments"])==sorted([68644,q])
                    for ref,q in zip(boundary,physical["physical"]["Q_y"][:2])),"Physical R must multiply Q_y once")
        endpoints[name]=boundary
    check=load(root/paths[3])
    require(check["all_passed"] and check["checks"]["actual_source"]["Pstar_once_and_physical_R_Rsquared_measures_once"],
            "Accepted Q4 physical-radius receipt required")
    return dict(passed=True,R_equals_epsilon_rho=True,epsilon_Lambda_equals_one=True,
        exact_logLambda="4*logPstar+1000",actual_micro_radius_handle=68644,
        Q4_RQy_C1_handles=endpoints,physical_R_applied_once=True,
        physical_R_Z_independent=True,
        same_loop_Q4_radius_graph_prefix=True,actual_rho_source_AST_consumed=True,
        core_recipe_receipt_consumed_not_replayed=True,source_producers_executed=False),pins


def build(root,parent):
    root=Path(root)
    require(sha(root/PARENT)==PARENT_SHA256 and load(root/PARENT)==parent,"Accepted complete mean Z3 parent required")
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256)}
    for path in (RECOVERY,Q4,LOOP):
        key=Path(path).name if path.startswith("experiments/") else path
        group="input_hashes" if path.startswith("experiments/") else "implementation_hashes"
        expected=parent[group][key]
        require(sha(root/path)==expected,"Actual Q recovery source changed")
        pins[path]=dict(path=path,expected_sha256=expected)
    contract,ast=provenance(root,parent)
    radius,radius_pins=radius_provenance(root,parent);pins.update(radius_pins)
    from mpmath.ctx_iv import MPIntervalContext
    c=MPIntervalContext();c.prec=600
    mean=parent["actual_micro_mean_Z3_bounds"]
    require(mean["passed"] and mean["actual_complete_mean_Z3_norm_provider"]
            and mean["amplitude_provenance"]["passed"],"Complete physical mean and same E cap required")
    logs=mean["physical_scale_logs"]
    C,logN,logPE=[read(c,logs[key]) for key in ("selected_logC","N","Pstar_E_upper")]
    de=read(c,parent["actual_core_signed_direction"]["parameters"]["delta"])
    Z=c.mpf([c.mpf(37)/100,c.mpf(93)/250])
    require(de.a>0 and de.b<1 and C.a>0 and logN.a>0,"Actual positive parameter enclosures required")
    logP=read(c,parent["actual_core_signed_direction"]["parameters"]["logPstar"])
    domain=parent["actual_core_signed_direction"]["micro_core_argument_domain"]
    require(domain["passed"] and domain["native_s_interval"]==["0","2"],"Actual R geometry domain required")
    logR=c.ln(read(c,domain["core_argument_upper"]).b)-(4*logP+1000)
    out={}
    for chart in ("first","second"):
        source=mean["charts"][chart]
        physical_B=saved(c,source["source_physical_B_log_caps"]["B"])
        physical_By=saved(c,source["source_physical_B_log_caps"]["slow_By"])
        U=[v-logN for v in physical_B]
        m=saved(c,source["log_caps"]["physical_complete_delta_m"])
        my=saved(c,source["log_caps"]["physical_total_y_delta_m"])
        q,proof=recover(c,U[:3],m,Z,de,2)
        phase=parent["actual_micro_phase_inverse_bounds"]["charts"][chart]
        ph={r["ordinary_Z_order"]:read(c,r["log_absolute_upper_bound"])
            for r in phase["mixed_log_caps"]["B_over_Pstar_phi_over_frozen_E0"] if r["physical_y_order"]==0 and not r["exact_zero"]}
        require(set(ph)==set(range(4)) and phase["E_product_rule_retained"],"Actual Bphi ordinary rows required")
        Uy=[_logsum(c,[physical_By[k]-logN,logPE+ph[k]]) for k in range(2)]
        qy,yproof=recover(c,Uy,my[:3],Z,de,1)
        rqy=[value+logR for value in qy]
        fields=dict(physical_delta_U=U,physical_delta_m=m,physical_delta_m_y=my,
                    physical_delta_Q=q,physical_delta_U_y=Uy,physical_delta_Q_y=qy,
                    physical_R_delta_Q_y=rqy)
        threshold=-1024*c.ln(2)
        require(all(v.b<(-600*C).a for v in U) and all(v.b<(-2600*C).a for v in q),"Physical U/Q increment strict log bound failed")
        require(all(value.b<threshold.a for value in U+q),"Local U/Q increment radius failed")
        out[chart]=dict(passed=True,source_chart=chart+"_micro",Q_recovery=proof,Q_y_recovery=yproof,
            log_caps={key:pack(values) for key,values in fields.items()},
            log_cap_ratios={key:[dict(ordinary_Z_order=k,log_cap_over_selected_logC=record(v/C)) for k,v in enumerate(values)] for key,values in fields.items()},
            physical_delta_U_Z0_Z3_below_2_power_minus1024=True,
            physical_delta_Q_Z0_Z2_below_2_power_minus1024=True,
            strict_log_envelopes=dict(U_multiplier=-600,Q_multiplier=-2600,reference="selected positive logC",passed=True),
            outer_RQy_smallness=[dict(ordinary_Z_order=k,below_2_power_minus1024=bool(value.b<threshold.a)) for k,value in enumerate(rqy)],
            non_small_bounds_retained=True,total_y_fast_Bphi_term_not_divided_by_N=True,
            physical_R_applied_once=True)
    return dict(passed=True,source_family=parent["source_family"],source_pins=pins,
        source_contract=contract,source_AST=ast,radius_provenance=radius,new_linear_Q_controls=new_controls(),
        domain=mean["domain"],charts=out,actual_delta=record(de),positive_L=record(1-de*Z*Z),
        physical_R_upper_log=record(logR),
        actual_physical_delta_U_Z3_error_provider=True,actual_physical_delta_Q_Z2_error_provider=True,
        actual_physical_delta_Q_y_Z1_error_provider=True,actual_physical_R_delta_Q_y_Z1_error_provider=True,
        signed_function_values_evaluated=False,zero_is_error_box_center_not_selected_value=True,
        all_outer_RQy_errors_certified_small=all(row["below_2_power_minus1024"] for chart in out.values() for row in chart["outer_RQy_smallness"]),six_I1_defect_bounds=False,
        complete_full_Q_Z2_numeric_values=False,full_C5_numeric_containment=False,
        global_five_moment_closure=False,full_goal_complete=False,
        old_signed_sources_and_all_error_caps_unchanged=True,accepted_ancestors_replayed=False)
