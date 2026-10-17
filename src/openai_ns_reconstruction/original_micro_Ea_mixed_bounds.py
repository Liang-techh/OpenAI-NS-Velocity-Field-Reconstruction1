"""Actual micro E/a mixed bounds from ordinary source jets and Xh norms."""
from __future__ import annotations
import math
from pathlib import Path
from .original_core_signed_direction import require, sha, digest, read, record, source_guard
from .original_micro_direction_C0 import load, BASE, PREFIX, BRIDGE, INLET, LOOP, CORE, UNIFORM, SIGMA, RATES_SOURCE

PARENT = PREFIX+"compliant_current_original_micro_direction_C0.json"
PARENT_SHA256 = "45e41b113fa861934497c28cf5e2cd29a02e14874723d20fd13ab5ddcf6cef52"
HELPER = "src/openai_ns_reconstruction/original_ordinary_jet_bounds.py"
HELPER_SHA256 = "52b724d646218603de1dfc03568b06d5c29e7931b46e8c291d522bddce0658d2"
TRANSFER = PREFIX+"compliant_core_transfer.json"
LINEAR = PREFIX+"shared_linear_resolvent.json"
TUBE = PREFIX+"shared_analytic_tube.json"


def analytic_rows(c, uniform, major, linear, tube):
    """Consume the old rows; compute only the three missing mixed rows."""
    existing = {(r["scaled_radial_order"],r["axial_order"]):r for r in uniform["derivative_bounds"]}
    h = read(c,tube["Xh_parameter"]); radius = c.mpf(41)/10
    correction = read(c,major["scaled_map_size_upper"])
    PhiNorm = read(c,linear["Phi_model_Xh_norm_upper"])+correction
    PsiNorm = read(c,major["fresh_Psi_model_Xh_norm_upper"])+correction
    eps = read(c,major["epsilon"])
    rows, fresh = {}, []
    for i in (0,1):
        for k in range(5):
            if (i,k) in existing:
                rows[i,k] = {field:read(c,existing[i,k][field]).b for field in
                    ("Phi_absolute_bound","Uz_minus_U0_absolute_bound")}
            else:
                factor = c.mpf(math.factorial(i+k))/((k+1)**2*h**k*20**i)/(1-radius/20)**(i+k+1)
                rows[i,k] = dict(Phi_absolute_bound=(PhiNorm*factor).b,
                                  Uz_minus_U0_absolute_bound=(eps*PsiNorm*factor).b)
                fresh.append(dict(scaled_radial_order=i,ordinary_Z_order=k,embedding=record(factor),
                    Phi_absolute_bound=record(PhiNorm*factor),
                    Uz_minus_U0_absolute_bound=record(eps*PsiNorm*factor)))
    require(len(fresh)==3 and {(r["scaled_radial_order"],r["ordinary_Z_order"]) for r in fresh} == {(0,4),(1,3),(1,4)},
            "Only the required new all-coefficient bounds should be computed")
    proof = dict(passed=True,new_rows=fresh,old_rows_consumed_without_recomputation=True,
        norm_definition=tube["norm_definition"],
        embedding="(i+k)!/[20^i*h^k*(k+1)^2*(1-rho_max/20)^(i+k+1)]",
        correction_is_uniform_self_map_bound_on_unit_ball=True,
        fixed_point_argument="X*=X0+N(X*), X* in unit ball; sup_ball ||N(X)|| <= scaled_map_size_upper",
        no_new_coefficient_or_temporal_recursion=True,ordinary_derivatives_no_factorial_discount=True,
        rho_domain=["0","41/10"],Z_domain=["-1","1"])
    return rows, proof


def source_jets(c, parent, uniform, rows):
    from .original_ordinary_jet_bounds import jet_add,jet_mul,jet_div,jet_exp,jet_log,absolute_majorant,bound_product
    ctx=c; zero=c.mpf(0); const=lambda x:[c.mpf(x)]+[zero]*4
    sym=lambda x:c.mpf([-1,1])*c.fabs(x).b
    cap=lambda key,i,k:rows[i,k][key]
    old=parent["actual_micro_direction_C0"]["comparison_direction"]["bounds"]
    core=parent["actual_core_signed_direction"]
    ph=read(c,core["cases"]["core_tube"]["Phi"]);pr=read(c,core["cases"]["core_tube"]["Phi_rho"])
    z=[c.mpf([c.mpf(37)/100,c.mpf(93)/250]),c.mpf(1),zero,zero,zero]
    j=read(c,core["parameters"]["j"]);de=read(c,core["parameters"]["delta"])
    eps=read(c,old["epsilon"]);rho=c.mpf([4,c.mpf(41)/10]); length=2*c.mpf(1)/1000
    Phi=[ph]+[sym(cap("Phi_absolute_bound",0,k)) for k in range(1,5)]
    Pr=[pr]+[sym(cap("Phi_absolute_bound",1,k)) for k in range(1,5)]
    U=[4*z[0]+j+sym(cap("Uz_minus_U0_absolute_bound",0,0)),
       4+sym(cap("Uz_minus_U0_absolute_bound",0,1))]+[sym(cap("Uz_minus_U0_absolute_bound",0,k)) for k in range(2,5)]
    Ur=[sym(cap("Uz_minus_U0_absolute_bound",1,k)) for k in range(5)]
    coreLog=jet_log(c,Phi);ratio=jet_div(c,Pr,Phi)
    barLog=[c.ln(read(c,old["Phi_comparison"]))]+[
        coreLog[k]+sym(length*rho.b*c.fabs(ratio[k]).b) for k in range(1,5)]
    barPhi=jet_exp(c,barLog,base=read(c,old["Phi_comparison"]))
    barU=[U[k]+sym(length*rho.b*c.fabs(Ur[k]).b) for k in range(5)]
    barU[0]=read(c,old["U_comparison"])
    globalPhi=[read(c,uniform["actual_Phi_upper"]).b]+[cap("Phi_absolute_bound",0,k) for k in range(1,5)]
    productIn=bound_product(c,globalPhi,absolute_majorant(c,U))
    productBar=jet_mul(c,barPhi,barU)
    H=[read(c,old["JH"])]+[sym(max(globalPhi[k].b,c.fabs(barPhi[k]).b)) for k in range(1,5)]
    K=[read(c,old["JK"])]+[sym(max(productIn[k].b,c.fabs(productBar[k]).b)) for k in range(1,5)]
    M=[read(c,old["M"])]+[sym(max(c.fabs(U[k]).b,c.fabs(barU[k]).b)) for k in range(1,5)]
    d=jet_add(c,const(1),[-q for q in jet_mul(c,z,z)])
    L=jet_add(c,const(1),[-de*q for q in jet_mul(c,z,z)])
    H0=jet_add(c,[(1-de)*q/2 for q in z],jet_mul(c,d,jet_add(c,[4*q for q in z],const(j))))
    denominator=jet_add(c,jet_mul(c,H0,H0),const((j/500)**2))
    g=jet_div(c,jet_mul(c,L,H0),denominator)
    alpha=c.mpf([0,1])
    barPhiY=[alpha*rho*q for q in jet_mul(c,barPhi,ratio)]
    barUY=[alpha*rho*q for q in Ur]
    HY=[2*(barPhi[k]-H[k]) for k in range(5)]
    KY=[2*(productBar[k]-K[k]) for k in range(5)]
    MY=[barU[k]-M[k] for k in range(5)]
    values=dict(core_Phi=Phi,core_Phi_rho=Pr,core_U=U,core_U_rho=Ur,core_logPhi=coreLog,
        comparison_logPhi=barLog,comparison_Phi=barPhi,comparison_U=barU,
        JH=H,JK=K,M=M,JH_y=HY,JK_y=KY,M_y=MY,
        comparison_Phi_y=barPhiY,comparison_U_y=barUY,gaxis=g,L=L,d=d,Z=z)
    return values,dict(epsilon=eps,rho=rho,delta=de), dict(passed=True,
        ordinary_Z_orders=list(range(5)),radius_and_alpha_Z_independent=True,
        complete_inlet_derivatives_bounded_from_original_normalized_integrals=True,
        integral_measure_bound="0<=y<=2hb<2/1000",actual_U0_Z1_equals_four=True,
        U_minus_U0_rows_already_include_epsilon=True,second_comparison_not_replaced_by_core=True,
        positive_Duhamel_kernel_bounds_each_differentiated_source=True,
        no_derivative_of_C0_cap=True,source_y_ODEs="JH_y=2(barphi-JH); JK_y=2(barphi*barU-JK); M_y=barU-M")


def direction_jets(c, values, scalar, parent):
    from .original_ordinary_jet_bounds import jet_add,jet_mul,jet_reciprocal,jet_derivative
    zero=[c.mpf(0)]*4
    const=lambda x:([c.mpf(x)]+zero[1:],zero)
    def val(name):
        if name in ("Z","d","L","gaxis"):
            return values[name][:4],zero
        require(name+"_y" in values,"Missing physical-y source jet: "+name)
        return values[name][:4],values[name+"_y"][:4]
    add=lambda a,b:(jet_add(c,a[0],b[0]),jet_add(c,a[1],b[1]))
    neg=lambda a:([-x for x in a[0]],[-x for x in a[1]])
    def mul(a,b):return (jet_mul(c,a[0],b[0]),jet_add(c,jet_mul(c,a[1],b[0]),jet_mul(c,a[0],b[1])))
    def inv(a):
        q=jet_reciprocal(c,a[0]);return q,[-x for x in jet_mul(c,jet_mul(c,q,q),a[1])]
    scale=lambda a,x:mul(a,const(x))
    derivative=lambda name:(jet_derivative(c,values[name]),jet_derivative(c,values[name+"_y"]))
    z,d,L,g=[val(k) for k in ("Z","d","L","gaxis")]
    H,K,M=[val(k) for k in ("JH","JK","M")]
    HZ,KZ,MZ=[derivative(k) for k in ("JH","JK","M")]
    ph=val("comparison_Phi");de=scalar["delta"];eps=scalar["epsilon"]
    A=add(add(scale(H,1-de/2),neg(scale(mul(z,HZ),(1-de)/2))),
          add(neg(mul(d,KZ)),scale(mul(z,K),2*de-1)))
    B=add(scale(mul(z,H),(1-de)/2),mul(d,K))
    W=add(const(1),neg(add(scale(mul(z,M),1-de),mul(d,MZ))))
    T=add(scale(mul(mul(g,B),inv(ph)),c.mpf(1)/2),
          scale(add(neg(W),scale(mul(A,inv(ph)),c.mpf(1)/2)),eps))
    rho=( [scalar["rho"]]+zero[1:], [scalar["rho"]]+zero[1:] )
    result=mul(mul(rho,inv(L)),T)
    result[0][0]=read(c,parent["actual_micro_direction_C0"]["comparison_direction"]["bounds"]["R_D_comparison"])
    return result


def root_Ea_bounds(c, values, directions, parent, major):
    from .original_ordinary_jet_bounds import exp_bound,bound_product
    D,Dy=directions
    ab=lambda v:c.fabs(v).b
    normD,normY=[ab(v) for v in D],[ab(v) for v in Dy]
    hbcap=c.mpf(1)/1000;lam=read(c,major["Lambda"])
    nativeLog=[c.mpf(0)]+[ab(values["core_logPhi"][k])+hbcap*normD[k] for k in range(1,4)]
    logE=[c.mpf(0)]+[lam*ab(values["gaxis"][k-1])+nativeLog[k] for k in range(1,4)]
    relativeE=exp_bound(c,logE,base_upper=c.mpf(1))
    oneMinusA=[1+normD[0]]+normD[1:]
    relativeEy=[v/2 for v in bound_product(c,relativeE,oneMinusA)]
    floor=D[0].a
    relativeA=[c.mpf(1)]+[v/floor for v in normD[1:]]
    loghb=read(c,parent["actual_signed_lower_m15"]["source_log_bases"][0])
    first,second=[],[]
    for dk,yk in zip(normD,normY):
        dx,dy=max(c.mpf(1),dk),max(c.mpf(1),yk)
        left=c.ln(32*dx)-loghb;right=c.ln(dy)
        first.append(c.mpf(max(left.b,right.b))+c.ln(2))
        second.append(loghb+c.ln(dy))
    def pack(jet,start=0,key="absolute_upper_enclosure"):
        return [dict(ordinary_Z_order=k,**{key:record(v)}) for k,v in enumerate(jet,start)]
    return dict(passed=True,ordinary_Z_orders=[0,1,2,3],physical_y_orders=[0,1],
        E_Z_over_E=pack(relativeE),E_yZ_over_E=pack(relativeEy),a_Z_over_a=pack(relativeA),
        a_Z_absolute_factored_chi_R=pack(normD),
        a_yZ_log_absolute_upper_first=pack(first,key="log_absolute_upper_bound"),
        a_yZ_log_absolute_upper_second=pack(second,key="log_absolute_upper_bound"),
        a_y_first_bound="32*abs(Dbar_Zk)/hb+abs(Dbar_yZk)",
        a_y_second_bound="hb*abs(Dbar_yZk); chi_y=0 on s in [1,2]",
        logE_ordinary_Z_derivative_absolute_bounds=pack(logE[1:],start=1),
        relative_exp_base_one_is_ratio_not_selected_E_value=True,
        source_identity_E_y="E_y=(1-a)*E/2",source_identity_a="a=chi_R*Dbar",
        ordinary_Z_does_not_differentiate_chi_R=True,physical_y_inverse_hb_retained=True,
        log_absolute_caps_not_signed_derivative_values=True,F0_hb_N_not_materialized=True,
        actual_pre_finite_N_background_Ea=True,finite_N_modified_E_N_not_bounded=True,
        b_t0_q_p2_bounds_complete=False)


def sigma_derivative_proof():
    import sympy as s
    u=s.symbols("u",positive=True)
    require(s.simplify(s.diff(u**3*s.exp(-u*u),u)-u*u*(3-2*u*u)*s.exp(-u*u))==0,
            "Flat sigma derivative envelope proof differs")
    return dict(passed=True,global_absolute_derivative_upper=32,
        proof="For 0<x<=1/2, sigma<=exp(4-1/x^2), (1-x)^-3<=8. With u=1/x>=2, u^3 exp(-u^2) and exp(-u^2) decrease, so sigma'<=2 exp(4)*(8+8)exp(-4)=32. Reflection covers the other half; flat ends have derivative0.",
        original_sigma_not_its_sample_or_cap=True)


def build(root,parent):
    root=Path(root)
    require(sha(root/PARENT)==PARENT_SHA256 and load(root/PARENT)==parent,"Accepted two-micro C0 source required")
    require(sha(root/HELPER)==HELPER_SHA256,"Changed reviewed ordinary-jet helper")
    pins={PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256),HELPER:dict(path=HELPER,expected_sha256=HELPER_SHA256)}
    paths=[TRANSFER,TRANSFER[:-5]+".py",LINEAR,TUBE,UNIFORM,UNIFORM[:-5]+".py",CORE,SIGMA,RATES_SOURCE]
    paths += [base+suffix for base in (BRIDGE,INLET,LOOP) for suffix in (".py",".json.gz","_check.json")]
    for path in paths:
        expected=parent["input_hashes"][Path(path).name]
        require(sha(root/path)==expected,"Changed same-source analytic input: "+path)
        pins[path]=dict(path=path,expected_sha256=expected)
    major,linear,tube,uniform=[load(root/p) for p in (TRANSFER,LINEAR,TUBE,UNIFORM)]
    require(linear["input_sha256"]==sha(root/TUBE) and linear["analytic_linear_inverse_bound_certified"],
            "Linear model norm must use this exact analytic tube")
    require(major["ball_radius"]==1 and len(major["terms"])==20 and major["all_paper_8_50_terms_included"] and
            major["contraction_proved"] and major["bounds_for_new_implicit_continuous_source"] and
            major["analytic_core_family_sha256"]==tube["analytic_core_family_sha256"]==uniform["analytic_core_family_sha256"] and
            uniform["implicit_source_sha256"]==parent["source_family"]["implicit_source_sha256"],"One actual analytic fixed point required")
    from mpmath.ctx_iv import MPIntervalContext
    from .original_ordinary_jet_bounds import symbolic_controls
    c=MPIntervalContext();c.prec=600
    rows,analytic=analytic_rows(c,uniform,major,linear,tube)
    values,scalar,state=source_jets(c,parent,uniform,rows)
    directions=direction_jets(c,values,scalar,parent)
    Ea=root_Ea_bounds(c,values,directions,parent,major)
    return dict(passed=True,source_family=parent["source_family"],source_pins=pins,
        source_AST=source_guards(root),source_provenance=source_provenance(root,parent),
        ordinary_jet_controls=symbolic_controls(),analytic_derivative_provider=analytic,
        comparison_source_jet_proof=state,sigma_derivative_bound=sigma_derivative_proof(),
        comparison_source_Z0_Z4_enclosures={key:[record(v) for v in jet] for key,jet in values.items()},
        comparison_Dbar_Z0_Z3=[record(v) for v in directions[0]],
        comparison_Dbar_yZ0_Z3=[record(v) for v in directions[1]],actual_Ea_mixed_bounds=Ea,
        domain=dict(s=[0,2],Z=["37/100","93/250"],rho=["4","41/10"]),
        highest_actual_mixed_core_order=5,actual_Ea_root_bound_provider=True,complete_mixed_root_provider=False,
        numerical_Cy_Z3=False,numerical_mean_Z3=False,numerical_QZZ=False,global_five_moment_closure=False,
        inherited_signed_sources_and_all_errors_unchanged=True,accepted_ancestors_replayed=False,full_goal_complete=False)


def source_guards(root):
    from .original_micro_direction_C0 import source_guards as accepted_guards
    # Only AST reads; no accepted arithmetic, graph producer or seam proof is executed.
    old=accepted_guards(root)
    selected={k:old[k] for k in ("comparison_densities","histories","directions","core_inlets","core_gauge","native_E","native_phi","sigma","loop_roots")}
    return dict(**selected,self_map=source_guard(root,TRANSFER[:-5]+".py","run",(
        "Bphi=get(linear,'Phi_model_Xh_norm_upper')+1","Bpsi=Psi_model+1",
        "value=c*Bphi**p*Bpsi**q","ball_center='paper X0=(Phi0,Psi0)',ball_radius=1")),
        norm_semantics=source_guard(root,UNIFORM[:-5]+".py","run",(
        "correction_norm = get(major,'scaled_map_size_upper')",
        "phi_norm = get(linear,'Phi_model_Xh_norm_upper')+correction_norm",
        "psi_norm = get(major,'fresh_Psi_model_Xh_norm_upper')+correction_norm",
        "Uz_minus_U0_absolute_bound=eps*psi_norm*factor")),embedding=source_guard(root,UNIFORM[:-5]+".py","embedding",(
        "math.factorial(i+k)","(k+1)**2 * h**k * 20**i","(1-r/20)**(i+k+1)")),
        comparison=source_guard(root,BRIDGE+".py","__init__",(
        "alpha=g.sub(g.one,sigma(g.sub(u,g.one)))", "rho=lambda q:g.mul(c(4),g.unary('exp',g.mul(hb,q)))",
        "self.initial7={name:extend(d,jet,7) for name,jet in inlet_owner.inlet_moments.items()}",
        "chi=g.add(g.one,g.mul(g.sub(hb,g.one),sigma(s)))",
        "R=g.mul(Ra,g.unary('exp',g.mul(hb,s)))")),
        native=source_guard(root,BRIDGE+".py","actual_fields",(
        "angular=a.scale(known['D_over_R'],g.mul(chi,radius,jacobian,g.constant('-1/2')))",)),
        physical_y=source_guard(root,BRIDGE+".py","generic_fields",(
        "radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y']))",)))


def source_provenance(root,parent):
    """Read frozen DAGs/receipts; bind these bounds to their original ports."""
    bridge,inlet,loop=[load(root/(base+".json.gz")) for base in (BRIDGE,INLET,LOOP)]
    for base in (BRIDGE,INLET,LOOP):
        receipt=load(root/(base+"_check.json"))
        require(receipt["all_passed"] and receipt["source_family"]==parent["source_family"],
                "Current C6/C5 source receipt mismatch")
    bn,inn,ln=[r["exact_graph_nodes"] for r in (bridge,inlet,loop)]
    require(ln[:len(bn)]==bn and bn[:len(inn)]==inn,"Actual DAG prefixes differ")
    bb,lb=bridge["prescribed_bridge_C6_bindings"],loop["prescribed_loop_C5_bindings"]
    handles=parent["actual_micro_direction_C0"]["actual_micro_source_handles"]
    require(bb["comparison_micro"]==handles["comparison"] and bb["micro"]==handles["native"] and
            lb["charts"]["first_micro"]==handles["first"] and lb["charts"]["second_micro"]==handles["second"],
            "New source bounds must bind the accepted current comparison/native/root functions")
    gates=parent["actual_micro_direction_C0"]["actual_root_domain_gates"]
    require(gates["passed"] and gates["s"]==[0,2] and gates["hb_strictly_positive_and_below_one"],
            "Accepted actual geometry gate missing")
    require(bb["micro"]["domain"]==[0,2] and bb["micro"]["physical_measure"]==8 and
            handles["first"]["roots"]==handles["second"]["roots"],"Common physical source roots required")
    return dict(passed=True,current_checked_C6_C5_source=True,exact_graph_prefixes_equal=True,
        parent_handle_identity=True,comparison=handles["comparison"],native=handles["native"],
        first_roots=handles["first"]["roots"],second_roots=handles["second"]["roots"],
        actual_root_Ea_not_independent_surrogate=True,physical_y_measure_node=8,
        comparison_physical_y_is_hb_inverse_s_derivative=True,
        original_nonzero_inlets_retained=True,source_producers_replayed=False)
