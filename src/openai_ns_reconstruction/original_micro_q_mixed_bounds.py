"""Actual micro q and scaled Poisson input mixed logarithmic bounds."""
from __future__ import annotations
import math
from pathlib import Path
from .original_core_signed_direction import require,sha,read,record,source_guard
from .original_micro_direction_C0 import load,PREFIX,LOOP,BRIDGE,INLET

PARENT=PREFIX+"compliant_current_original_micro_b_p2_mixed_bounds.json"
PARENT_SHA256="fc28f40dcf58e3a5123725361baca36005902f0733080c2fff6852293a703e5d"
LOG_HELPER="src/openai_ns_reconstruction/original_log_mixed_jets.py"
SIGMA_HELPER="src/openai_ns_reconstruction/original_flat_sigma_derivative_bounds.py"
INTERFACE=PREFIX+"compliant_current_core_first_interface"
INTERIOR=PREFIX+"compliant_current_core_interior_moments"
HELPER_HASHES={LOG_HELPER:"ceeebf88a02207b02fedfdff727bc2ae1d09f4a52ebf27b8a84ff8876719361a",
              SIGMA_HELPER:"c528a552f8c9942165c799aa292656afe1b9233b49aa26c7874c2576ede2b64d"}


def pack(jet):
    return [dict(physical_y_order=i,ordinary_Z_order=k,exact_zero=(v is None),
                 log_absolute_upper_bound=None if v is None else record(v)) for (i,k),v in jet.items()]


def parameters(c,parent):
    old=parent["actual_micro_direction_C0"]["actual_C0_errors"]
    C=read(c,parent["actual_micro_b_p2_mixed_bounds"]["whole_tube_amplitude"]["logC"])
    out=dict(C=C,eta=read(c,old["log_eta"]),hb=read(c,parent["actual_signed_lower_m15"]["source_log_bases"][0]),
             dstar=read(c,old["log_dstar"]),N=read(c,old["log_common_N"]))
    require(old["eta_positive_and_below_one"] and old["dstar_strict_positive"] and C.a>0,
            "Actual positive canonical parameters required")
    require((-245*C).b<out["eta"].a and out["eta"].b<(-243*C).a and
            (-101*C).b<out["hb"].a and out["hb"].b<(-99*C).a and out["N"].a>(6131*C).b,
            "Saved actual parameter log envelopes differ")
    return out,dict(passed=True,actual_logs={k:record(v) for k,v in out.items()},
        log_ratios_to_selected_logC={k:record(out[k]/C) for k in ("eta","hb","N")},
        approximate_ratios_not_exact_equalities=True,
        coarse_strict_envelopes="-245C<logeta<-243C; -101C<loghb<-99C; logN>6131C",
        C_denotes_selected_positive_logC=True,eta_hb_N_not_materialized=True)


def source_inputs(c,parent,chart,engine,logs):
    """Source derivative bounds; the active branch only tightens value caps."""
    from .original_log_mixed_jets import LogMixed
    Ea=parent["actual_micro_Ea_mixed_bounds"]["actual_Ea_mixed_bounds"]
    cmp=parent["actual_micro_direction_C0"]["comparison_direction"]["bounds"]
    D=read(c,cmp["R_D_comparison"])
    logfloor=logs["hb"]+c.ln(D.a)
    logaUpper=c.ln(3)
    if chart=="second":logaUpper=c.mpf(min(logaUpper.b,(logs["hb"]+c.ln(D.b)).b))
    aNorm={(0,0):c.mpf(0)}
    for k in range(1,4):aNorm[0,k]=c.ln(read(c,Ea["a_Z_over_a"][k]["absolute_upper_enclosure"]).b)
    for k in range(4):
        aNorm[1,k]=read(c,Ea["a_yZ_log_absolute_upper_"+chart][k]["log_absolute_upper_bound"])-logfloor
    normalized=LogMixed.from_mapping(c,aNorm)
    inv=engine.scale(engine.reciprocal(normalized,c.mpf(0)),-logfloor)
    actualA={}
    for i,k in engine.indices:
        if (i,k)==(0,0):actualA[i,k]=logaUpper
        elif i==0:
            cap=read(c,Ea["a_Z_absolute_factored_chi_R"][k]["absolute_upper_enclosure"])
            actualA[i,k]=c.ln(cap.b)+(logs["hb"] if chart=="second" else 0)
        else:actualA[i,k]=read(c,Ea["a_yZ_log_absolute_upper_"+chart][k]["log_absolute_upper_bound"])
    root=parent["actual_micro_b_p2_mixed_bounds"]["actual_b_t0_p2_mixed_bounds"]
    def channel(value_key,y_key):
        result={}
        for k in range(4):
            result[0,k]=read(c,root[value_key][k]["log_absolute_upper_bound"])
            result[1,k]=read(c,root[y_key][k]["log_absolute_upper_bound"])
        return LogMixed.from_mapping(c,result)
    b=channel("b_Z_"+chart,"b_yZ_"+chart)
    bd=b.as_dict();bd[0,0]=c.mpf(min(c.ln(2).b,bd[0,0].b));b=LogMixed.from_mapping(c,bd)
    return dict(a=LogMixed.from_mapping(c,actualA),inverse_a=inv,b=b,
        p2=channel("p2_Z_both","p2_yZ_both"),t0=channel("t0_Z_both","t0_yZ_both"),
        loga_lower=logfloor,loga_upper=logaUpper),dict(passed=True,
        a_positive_lower_log=record(logfloor),active_a_upper_log=record(logaUpper),
        normalized_a_derivative_ratios=pack(normalized),
        reciprocal_uses_point_frozen_reference_a0_not_global_a_over_a=True,
        inverse_a_scaled_by_actual_positive_floor=True,active_value_caps="a<=3, |b|<=2, kappa<=3",
        higher_derivatives_consumed_from_source_not_differentiated_value_caps=True)


def scalar_poisson_controls():
    """Global scalar derivative bounds for h(u) and r(u), valid on all R."""
    import sympy as s
    x=s.symbols("x",real=True);out={}
    for name,P in (("hinv",s.Integer(1)),("r",x)):
        rows=[];caps=[]
        for k in range(5):
            poly=s.Poly(s.expand(P),x)
            cap=sum(abs(v) for v in poly.all_coeffs())
            require(poly.degree()<=2*k+1,"Scalar denominator does not dominate each monomial")
            exact=(1+x*x)**s.Rational(-1,2) if name=="hinv" else x*(1+x*x)**s.Rational(-1,2)
            require(s.simplify(s.diff(exact,x,k)-P/(1+x*x)**s.Rational(2*k+1,2))==0,
                    "Scaled Poisson scalar derivative differs")
            rows.append(str(P));caps.append(int(cap))
            P=s.expand((1+x*x)*s.diff(P,x)-(2*k+1)*x*P)
        out[name]=dict(derivative_upper_bounds=caps,numerator_polynomials=rows)
    return dict(passed=True,functions=out,
        proof="D^k f=P_k(u)/(1+u^2)^(k+1/2); each monomial degree<=2k+1 so absolute quotient<=sum absolute coefficients",
        globally_real_scalar_caps=True,no_sampling_or_near_unit_rounding=True)


def active_chart(c,parent,chart,engine,logs,sigma,scalar):
    from .original_log_mixed_jets import LogMixed
    roots,rootProof=source_inputs(c,parent,chart,engine,logs)
    kappa=engine.add(roots["a"],engine.mul(engine.mul(roots["b"],roots["b"]),roots["inverse_a"]))
    kd=kappa.as_dict();kd[0,0]=c.ln(3);kappa=LogMixed.from_mapping(c,kd)
    gd=kappa.as_dict();gd[0,0]=c.ln(4);gamma=LogMixed.from_mapping(c,gd)
    cd={key:(c.ln(3)-logs["eta"] if key==(0,0) else value-logs["eta"]) for key,value in kappa.items()}
    cutArgument=LogMixed.from_mapping(c,cd)
    cut=engine.compose_scalar(cutArgument,[c.ln(c.mpf(n)) for n in sigma["derivative_upper_bounds"]])
    quotient=engine.scale(engine.mul(gamma,roots["inverse_a"]),-c.ln(2))
    rootLower=logs["eta"]-c.ln(2)-roots["loga_upper"]
    rootUpper=c.ln(2)-roots["loga_lower"]
    qhat=engine.sqrt(quotient,rootLower,rootUpper)
    q=engine.mul(cut,qhat)
    u=engine.scale(engine.mul(roots["p2"],q),-logs["dstar"])
    h=engine.compose_scalar(u,[c.ln(c.mpf(n)) for n in scalar["functions"]["hinv"]["derivative_upper_bounds"]])
    r=engine.compose_scalar(u,[c.ln(c.mpf(n)) for n in scalar["functions"]["r"]["derivative_upper_bounds"]])
    alpha=engine.scale(engine.mul(q,h),c.ln(2))
    ratio=engine.add(engine.constant(c.mpf(0)),engine.mul(roots["t0"],roots["t0"]),engine.scale(engine.mul(q,q),c.ln(2)))
    K=engine.scale(engine.reciprocal(ratio,c.mpf(0)),-c.ln(2*c.pi))
    # These are lower bounds on positive quantities, not log-absolute caps.
    gapLower=-c.mpf(max(c.mpf(0),(2*u[0,0]).b))-c.ln(2)
    Dlower=2*gapLower-c.ln(4)
    lambdaLower=-c.ln(2*c.pi)-ratio[0,0]
    coarse=2048
    require(all(v is not None and v.b<(coarse*logs["C"]).a for _,v in q.items()),
            "Actual q mixed bounds exceed the stated coarse log envelope")
    require(gapLower.a > (-105*logs["C"]).b and Dlower.a > (-211*logs["C"]).b
            and lambdaLower.a > (-101*logs["C"]).b,"Actual free-angle conditioning log floor failed")
    require(all(v is not None and v.b<(2048*logs["C"]).a
                for jet in (u,h,r,alpha,K) for _,v in jet.items()),
            "Scaled Poisson mixed input log envelope failed")
    return dict(passed=True,chart=chart,ordinary_Z_orders=[0,1,2,3],physical_y_orders=[0,1],
        root_source_inputs=rootProof,kappa_mixed_log_caps=pack(kappa),gamma_mixed_log_caps=pack(gamma),
        cut_argument_mixed_log_caps=pack(cutArgument),sigma_cut_mixed_log_caps=pack(cut),
        qhat_mixed_log_caps=pack(qhat),q_mixed_log_caps=pack(q),
        u_mixed_log_caps=pack(u),hinv_mixed_log_caps=pack(h),r_mixed_log_caps=pack(r),
        alpha_mixed_log_caps=pack(alpha),K_mixed_log_caps=pack(K),
        active_domain="Delta<eta, a>0; gamma=2eta-Delta>eta",
        square_root_log_lower=record(rootLower),square_root_log_upper=record(rootUpper),
        q_coarse_strict_log_envelope=dict(multiplier=coarse,reference="selected positive logC",passed=True),
        scaled_Poisson_conditioning=dict(passed=True,
            one_minus_r_squared_log_lower=record(gapLower),
            full_free_angle_D_log_lower=record(Dlower),
            full_free_angle_lambda_log_lower=record(lambdaLower),
            exact_gap_identity="1-r^2=1/(1+u^2)>0",exact_D="1-2r*cos(psi)+r^2",
            D_lower_identity="D>=(1-|r|)^2>=(1-r^2)^2/4",
            lambda_lower_identity="lambda=K*(1+t^2)>=K>0",
            free_angle=["0","2*pi"],r_not_rounded_to_one=True,
            q_u_r_hinv_not_materialized=True,inverse_angle_derivatives_enclosed=False,
            coarse_strict_log_lower_multipliers=dict(one_minus_r_squared=-105,D=-211,lambda_value=-101),
            selected_positive_logC_is_reference=True),
        scaled_Poisson_inputs_common_strict_log_upper_multiplier=2048,
        caps_not_signed_values=True,sqrt_not_evaluated_on_flat_branch=True)


def collar_flat(c,root,parent):
    core=parent["actual_core_signed_direction"]
    D=read(c,core["cases"]["core_tube"]["Dbar_ordinary_Z"][0])
    domain=core["micro_core_argument_domain"]
    require(domain["passed"] and domain["native_s_interval"]==["0","2"] and
            read(c,domain["core_argument_upper"]).b<c.mpf(41)/10,
            "Actual rho(s) must lie inside the signed core receipt")
    native=parent["actual_native_collar_transport"]["domain"]
    require(native["sc_positive_and_less_than_2_power_minus10"],"Actual sc source bound required")
    require(c.mpf(3)/4/1024<c.mpf(1)/512 and c.mpf(999)/1000*D.a>3,
            "Natural collar shear margin failed")
    require(4-512**2<c.ln(c.mpf(1)/1000).a,"Exact sigma tail bound too weak")
    identity=parent["actual_native_core_angular_direction"]["angular_direction_identity"]
    require(identity["symbolic_proof"]["Dbar_identity_exact"] and identity["nonvanishing_domain"]["passed"],
            "Accepted core angular identity required")
    require(parent["actual_micro_direction_C0"]["actual_C0_errors"]["eta_positive_and_below_one"],"eta<1 missing")
    return dict(passed=True,natural_flat_neighborhood=["0","1/512"],one_sided_at_s_zero=True,
        forced_collar_strictly_inside_natural_neighborhood=True,sc_positive_less_than_2_power_minus10=True,
        comparison_alpha_exact_one=True,
        comparison_core_identity="FTC for barphi/barU; complete core moment inlets; same linear ODEs and uniqueness",
        same_core_angular_direction_identity_consumed=True,
        core_Dbar_positive_interval=record(D),chi_lower="999/1000",strict_a_lower=record(c.mpf(999)/1000*D.a),
        actual_rho_argument_domain=["4","41/10"],full_micro_argument_receipt_consumed=True,
        sigma_tail_inequality="sigma(s)<=exp(4-1/s^2) on 0<s<=1/2; at0 sigma=0 exactly",
        sigma_tail_log_upper=4-512**2,strict_Delta_lower=1,eta_strictly_below_one=True,
        q_exact_zero_on_natural_neighborhood=True,physical_y_and_Z_jets_through_total_order_four_zero=True,
        original_hard_collar_branch_not_smoothed_or_replaced=True)


def build(root,parent):
    root=Path(root)
    require(sha(root/PARENT)==PARENT_SHA256 and load(root/PARENT)==parent,"Accepted actual b/t0/p2 parent required")
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256)}
    for path,h in HELPER_HASHES.items():
        require(sha(root/path)==h,"Changed reviewed new log/sigma helper")
        pins[path]=dict(path=path,expected_sha256=h)
    for base in (LOOP,BRIDGE,INLET,INTERIOR):
        for suffix in (".py",".json.gz","_check.json"):
            path=base+suffix;h=parent["input_hashes"][Path(path).name]
            require(sha(root/path)==h,"Changed actual loop source")
            pins[path]=dict(path=path,expected_sha256=h)
    for suffix in (".py","_check.json"):
        path=INTERFACE+suffix;h=parent["input_hashes"][Path(path).name]
        require(sha(root/path)==h,"Changed accepted first-interface source")
        pins[path]=dict(path=path,expected_sha256=h)
    from mpmath.ctx_iv import MPIntervalContext
    from .original_log_mixed_jets import LogMixed,symbolic_controls
    from .original_flat_sigma_derivative_bounds import build as sigma_build
    c=MPIntervalContext();c.prec=600;engine=LogMixed(c)
    logs,param=parameters(c,parent);sigma=sigma_build(root,parent);pins.update(sigma["source_pins"])
    scalar=scalar_poisson_controls()
    return dict(passed=True,source_family=parent["source_family"],source_pins=pins,
        parameter_envelopes=param,ordinary_log_mixed_controls=symbolic_controls(),
        scalar_sigma_derivative_proof=sigma,scalar_scaled_Poisson_derivative_proof=scalar,
        charts={chart:active_chart(c,parent,chart,engine,logs,sigma,scalar) for chart in ("first","second")},
        natural_collar_flat_proof=collar_flat(c,root,parent),source_AST=source_guards(root),
        direct_source_provenance=source_provenance(root,parent),
        shear_flat_glue=dict(passed=True,flat_predicate="Delta>=eta or original first collar",
            source_q_and_all_mixed_derivatives_exact_zero_on_flat=True,
            gamma_at_shear_boundary_at_least_eta=True,sigma_flat_derivatives_through_order_four=True,
            active_formula_uniform_bounds_extend_through_flat_boundary=True),
        domain=dict(s=[0,2],Z=["37/100","93/250"],physical_y="log(R/Ra)=hb*s",
            first_original_chart=["sc/2","1"],second_original_chart=["1","2"],
            source_owned_zero_left_extension=["0","sc/2"]),
        ordinary_Z_orders=[0,1,2,3],physical_y_orders=[0,1],q_mixed_provider_complete=True,
        actual_scaled_Poisson_inputs_mixed_provider=True,
        actual_free_angle_positive_D_and_lambda_log_bounds=True,
        phase_or_inverse_mixed_provider_complete=False,finite_N_derivative_errors=False,
        numerical_Cy_Z3=False,numerical_mean_Z3=False,numerical_QZZ=False,
        full_C5_numeric_containment=False,global_five_moment_closure=False,full_goal_complete=False,
        inherited_signed_sources_and_all_errors_unchanged=True,accepted_ancestors_replayed=False)


def source_guards(root):
    return dict(q=source_guard(root,LOOP+".py","loop_functions",(
        "gamma=a.sub(a.fixed(g.mul(c(2),eta)),roots['Delta'])",
        "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,",
        "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),",
        "collar=g.mul(c('3/4'),sc) if chart=='first_micro' else None",
        "u=a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,",
        "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)")),
        roots=source_guard(root,LOOP+".py","generic_roots",(
        "kap=a.add(aa,a.div(a.mul(bb,bb),aa,",)),
        comparison=source_guard(root,BRIDGE+".py","__init__",(
        "alpha=g.sub(g.one,sigma(g.sub(u,g.one)))",
        "self.initial7={name:extend(d,jet,7) for name,jet in inlet_owner.inlet_moments.items()}")))


def source_provenance(root,parent):
    """Read checked archives and consume exact functional identities only."""
    gates=parent["actual_micro_direction_C0"]["actual_root_domain_gates"]
    require(gates["passed"] and gates["active_implies_Delta_less_eta"] and
            parent["actual_micro_direction_C0"]["both_micro_actual_a_and_E_positive"],
            "Actual positive roots and active Delta<eta gate required")
    lazy=load(root/(LOOP+"_check.json"))["checks"]["lazy_flat"]
    require(lazy["passed"] and lazy["lazy_original_collar_and_shear_flat_rows"]==144 and
            lazy["active_inverse_and_singular_sqrt_not_evaluated_on_flat_branch"],
            "Original checked lazy flat functions required")
    interface=load(root/(INTERFACE+"_check.json"));interior=load(root/(INTERIOR+"_check.json"))
    q=interface["phase0_ODE_and_exact_positive_width_pullback"]
    require(interface["all_passed"] and q["passed"] and q["first_comparison_alpha_equals_one"] and
            q["equality_by_functions_and_ODE_uniqueness_not_interval_overlap"] and
            q["original_full_three_scale_controls_AST_bound"] and q["exact_mixed4_pullback_identities"]==15,
            "Accepted alpha-one functional source identity required")
    require(interior["all_passed"] and interior["current_core_interior_source_theorem"]["passed"] and
            interior["current_core_interior_source_theorem"]["source_functions_not_chosen_from_enclosures"],
            "Actual normalized core primitives required")
    bridge,inlet,loop=[load(root/(path+".json.gz")) for path in (BRIDGE,INLET,LOOP)]
    for path in (BRIDGE,INLET,LOOP):
        check=load(root/(path+"_check.json"))
        require(check["all_passed"] and check["source_family"]==parent["source_family"],"Same checked C6/C5 source required")
    bb,ib,lb=bridge["prescribed_bridge_C6_bindings"],inlet["core_inlet_axial6_bindings"],loop["prescribed_loop_C5_bindings"]
    require(all(bb["initial7"][name][:7]==ib["inlet_moments"][name] for name in ib["inlet_moments"]) and
            all(bb["core7"][name][:7]==ib["core_functions"][name] for name in ib["core_functions"]),
            "Exact common core and complete inlets required")
    handles=parent["actual_micro_direction_C0"]["actual_micro_source_handles"]
    require(bb["comparison_micro"]==handles["comparison"] and bb["micro"]==handles["native"] and
            all(lb["charts"][key+"_micro"]==handles[key] for key in ("first","second")),
            "Actual source handles differ from current original functions")
    node=loop["exact_graph_nodes"][10]
    require(node["name"]=="sc" and node["Z_independent"] and node["quantity_is_exact_original_parameter_not_range_endpoint"] and
            node["source_family"]==parent["source_family"],"Original fixed sc must be the source parameter")
    return dict(passed=True,actual_active_Delta_less_eta_gate_consumed=True,
        actual_a_and_E_positive_gate_consumed=True,original_lazy_flat_checked_rows=144,
        active_inverse_and_sqrt_not_evaluated_on_flat=True,
        initial7_prefix_equals_complete_core_inlets=True,core7_prefix_equals_current_core=True,
        alpha_one_core_history_functional_identity_consumed=True,
        accepted_first_interface_receipt_only_not_replayed=True,
        alpha_one_history_uniqueness_proof=q["same_moment_history_join_proof"][:3],
        comparison_core_identification_scope="0<=s<=1/512 only; alpha=1 by actual source",
        no_second_comparison_core_identification=True,
        first_q_source_handles=handles["first"]["loop"]["q"],second_q_source_handles=handles["second"]["loop"]["q"],
        first_original_flat_predicate=handles["first"]["loop"]["flat_predicate"],
        second_original_flat_predicate=handles["second"]["loop"]["flat_predicate"],
        actual_eta_node=21,actual_sc_node=10,original_sc_is_Z_independent_exact_parameter=True,
        source_producers_or_interface_checker_replayed=False)
