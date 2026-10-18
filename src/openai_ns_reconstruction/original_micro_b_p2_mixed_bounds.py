"""Actual micro b/t0/p2: source-bound logarithmic mixed derivative caps."""
from __future__ import annotations
import math
from pathlib import Path
from .original_core_signed_direction import require, sha, read, record, source_guard
from .original_micro_direction_C0 import load, PREFIX, BRIDGE, INLET, CORE, CORE_CHECK, NORMS
from .original_micro_Ea_mixed_bounds import TRANSFER, TUBE
from .original_ordinary_jet_bounds import jet_add,jet_mul,jet_reciprocal,jet_div,jet_exp,jet_log,exp_bound,bound_product

PARENT=PREFIX+"compliant_current_original_micro_Ea_mixed_bounds.json"
PARENT_SHA256="a48e8ba9de2f2865598d78b6ba7cfcc94f95ec4b94cdda45b7acb2c168413576"
Z5="src/openai_ns_reconstruction/original_micro_source_Z5.py"
PRESSURE="src/openai_ns_reconstruction/original_pressure_axial_tube.py"
REVIEWED_HELPERS={
    Z5:"666a59c454de29d4be48daef3ed8adffad455c3de98bbd64afad40a2332992af",
    PRESSURE:"bd5b6ccc52a26d35ed0b89c98180c6c8e81cb8afbb14976f6eefa52f4c0e8685"}


def symmetric(c,v):
    return c.mpf([-1,1])*c.fabs(v).b


def absolute(c,jet):
    return [c.fabs(v).b for v in jet]


class Dual:
    """Ordinary Z jets paired with physical-y jets; never pad missing rows."""
    def __init__(self,c,order):self.c,self.n=c,order+1
    def fixed(self,v):
        z=[self.c.mpf(0)]*self.n
        return [self.c.mpf(v)]+z[1:],z
    def pair(self,v,y):
        require(len(v)>=self.n and len(y)>=self.n,"Missing actual mixed source row")
        return list(v[:self.n]),list(y[:self.n])
    def add(self,*args):
        out=self.fixed(0)
        for x in args:out=jet_add(self.c,out[0],x[0]),jet_add(self.c,out[1],x[1])
        return out
    def neg(self,x):return [-v for v in x[0]],[-v for v in x[1]]
    def mul(self,x,y):
        return jet_mul(self.c,x[0],y[0]),jet_add(self.c,jet_mul(self.c,x[1],y[0]),jet_mul(self.c,x[0],y[1]))
    def scale(self,x,v):return self.mul(x,self.fixed(v))
    def inv(self,x):
        q=jet_reciprocal(self.c,x[0])
        return q,[-v for v in jet_mul(self.c,jet_mul(self.c,q,q),x[1])]
    def div(self,x,y):return self.mul(x,self.inv(y))
    def source(self,values,name):return self.pair(values[name],values[name+"_y"])
    def shift(self,values,name):return self.pair(values[name][1:],values[name+"_y"][1:])


def amplitude(c,root,parent,values,major,tube,norms):
    old=parent["actual_micro_direction_C0"]["actual_C0_errors"]
    require(old["exact_anchored_G_nonnegative"] and old["F0_positive_and_at_most_one"],"Accepted anchored amplitude sign required")
    receipt=load(root/CORE_CHECK)
    scalar=receipt["checks"]["actual_oracle"]["exact_parameter_and_native_P0_identity"]["actual_source_scalar_and_unique_anchor_identity"]
    require(scalar["logC_matches_original_selected_source"] and scalar["anchor_leaf_uses_this_exact_source"],"Actual selected logC and anchor required")
    logC=read(c,norms["selected_logCstar"]);lam=read(c,major["Lambda"])
    G=read(c,tube["G_modulus_upper"])
    logF=c.mpf([(-logC-lam*G).a,(-logC).b])
    logs=[c.mpf(0)]+[lam*values["gaxis"][k-1] for k in range(1,5)]
    invRelative=exp_bound(c,logs,base_upper=c.mpf(1))
    f2Relative=exp_bound(c,[2*v for v in logs],base_upper=c.mpf(1))
    logcaps=[-2*logC+c.ln(c.fabs(v).b) for v in f2Relative]
    require(all(v.b<0 for v in logcaps),"Actual F0ÃƒÆ’Ã†â€™Ãƒâ€ Ã¢â‚¬â„¢ÃƒÆ’Ã‚Â¢ÃƒÂ¢Ã¢â‚¬Å¡Ã‚Â¬Ãƒâ€¦Ã‚Â¡ÃƒÆ’Ã†â€™ÃƒÂ¢Ã¢â€šÂ¬Ã…Â¡ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â² derivatives do not satisfy the coarse unit cap")
    F2=[c.mpf([0,1])]+[c.mpf([-1,1]) for _ in range(4)]
    logP=read(c,parent["actual_core_signed_direction"]["parameters"]["logPstar"])
    eps=read(c,major["epsilon"]);rho=c.mpf([4,c.mpf(41)/10]);R=eps*rho
    require(tube["common_complex_axis_poles_excluded"] and tube["tube_definition"]=="union of disks |Z-x|<=eta over real x in [-1,1]",
            "Saved whole-axis anchored G modulus required")
    return dict(logF=logF,logC=logC,inv_relative=invRelative,F2=F2,P2=c.exp(2*logP),rho=rho,R=R,epsilon=eps),dict(
        passed=True,Z=["37/100","93/250"],log_F0=record(logF),logC=record(logC),
        G_modulus_upper=record(G),log_Pstar=record(logP),
        exact_log_F0="-logC-Lambda*G; 0<=G<=Gbar on this tube",
        F0_squared_derivative_log_caps=[record(v) for v in logcaps],
        F0_squared_ordinary_Z0_Z4_coarse_absolute_cap_one=True,
        inverse_F0_relative_ordinary_Z0_Z4=[record(v) for v in invRelative],
        inverse_F0_not_materialized=True,selected_logC_not_replaced_by_minimum_guard=True,
        DAG_exp_argument_5_is_a_node_handle_not_literal_exp5=True)


def complete_inlet_norms(c,values,globalPhi):
    """Bounds on the COMPLETE normalized integrals, uniform in the inlet t."""
    # core_U contains U0(Z)+GLOBAL Xh epsilon*Psi bounds for 0<=rho<=4.1.
    # These are not point-only values at the inlet rho=4.
    Uglobal=absolute(c,values["core_U"])
    Phi2=bound_product(c,globalPhi,globalPhi)
    return dict(M=Uglobal,A=bound_product(c,Uglobal,Uglobal),
                B=[v/2 for v in Phi2],C=Phi2)


def add_histories(c,values,phi_name,U_name,inlets,names=("A","B","C")):
    """Original complete inlet plus positive kernels, with each exact rate."""
    phi,U=values[phi_name],values[U_name]
    phi2=jet_mul(c,phi,phi);U2=jet_mul(c,U,U)
    for name in names:
        rate=2 if name=="B" else 1
        source=U2 if name=="A" else phi2
        caps=[max(c.fabs(v).b,c.fabs(w).b/rate) for v,w in zip(inlets[name],source)]
        jet=[c.mpf([0,caps[0].b])]+[symmetric(c,v) for v in caps[1:]]
        values[name]=jet
        values[name+"_y"]=[source[k]-rate*jet[k] for k in range(len(jet))]
    return values


def drive(c,values,p0,scalars,order=4):
    """R*Q is the exact comparison full drive, after R*Lambda=rho."""
    a=Dual(c,order);zero=[c.mpf(0)]*(order+2)
    vals=dict(values)
    for name in ("Z","d","L","gaxis"):vals[name+"_y"]=zero
    z,d,L,g=[a.source(vals,k) for k in ("Z","d","L","gaxis")]
    U,M,A,B,C=[a.source(vals,k) for k in ("comparison_U","M","A","B","C")]
    Mz,Az,Bz,Cz=[a.shift(vals,k) for k in ("M","A","B","C")]
    de=read(c,scalars["delta_record"])
    W=a.add(a.fixed(1),a.neg(a.add(a.scale(a.mul(z,M),1-de),a.mul(d,Mz))))
    hydro=a.div(a.add(a.neg(a.mul(W,U)),a.scale(a.add(M,a.neg(a.mul(z,Mz))),(1-de)/2),
                         a.scale(a.mul(z,A),2*de),a.neg(a.mul(d,Az))),a.scale(L,2))
    p=a.pair(p0,zero);pz=a.pair(p0[1:],zero)
    Pr=a.add(a.scale(a.mul(z,p),2*(1+de)),a.neg(a.mul(d,pz)))
    pressure=a.div(Pr,a.scale(L,2))
    S0=a.add(a.scale(a.mul(z,B),-2*de),a.mul(d,Bz),a.scale(a.mul(z,C),2*(1+de)),a.neg(a.mul(d,Cz)))
    rho=a.pair([scalars["rho"]]+zero,[scalars["rho"]]+zero)
    R=a.pair([scalars["R"]]+zero,[scalars["R"]]+zero)
    F2=a.pair(scalars["F2"],zero)
    swirl=a.mul(F2,a.add(a.div(a.mul(R,S0),a.scale(L,2)),
        a.neg(a.div(a.mul(a.mul(a.mul(rho,g),d),a.add(B,a.neg(C))),L))))
    return a.add(hydro,a.scale(pressure,scalars["P2"]),swirl)


def native(c,values,D,Q,scalars,inlet_M):
    """Native angular equation first, then prescribed axial integral and history."""
    kmax=4;phCore=values["core_Phi"][:5];coreLog=jet_log(c,phCore)
    hbcap=c.mpf(1)/1000;Dnorm=absolute(c,D[0]);floor=phCore[0].a/2
    require(hbcap*Dnorm[0] < c.ln(2).a,"Native phi floor proof unresolved")
    logphi=[c.ln(c.mpf([floor,phCore[0].b]))]+[coreLog[k]+symmetric(c,hbcap*Dnorm[k]) for k in range(1,5)]
    phi=jet_exp(c,logphi,base=c.mpf([floor,phCore[0].b]))
    phiy=[-v/2 for v in jet_mul(c,D[0],phi)]
    # chi_R is Z-independent and in [0,1]; do not differentiate its cap.
    phiy=[c.mpf([0,1])*v for v in phiy]
    ratio=jet_div(c,phi,values["comparison_Phi"][:5])
    Uycap=[scalars["R"].b*v for v in bound_product(c,ratio,Q[0])]
    Uycaps=absolute(c,Uycap)
    errors=[2*hbcap*v for v in Uycaps]
    require(all(v.b<1 for v in errors),"Native prescribed U integral bound unresolved")
    U=[values["core_U"][k]+symmetric(c,errors[k]) for k in range(5)]
    out=dict(core_U=values["core_U"][:5],phi=phi,phi_y=phiy,U=U,U_y=[symmetric(c,v) for v in Uycaps])
    # Native histories are never identified with comparison histories.
    Mmax=[max(c.fabs(v).b,c.fabs(w).b) for v,w in zip(inlet_M[:5],U)]
    out["M"]=[c.mpf([-Mmax[0].b,Mmax[0].b])]+[symmetric(c,v) for v in Mmax[1:]]
    out["M_y"]=[U[k]-out["M"][k] for k in range(5)]
    for name in ("Z","d","L","gaxis"):
        out[name]=values[name][:5];out[name+"_y"]=[c.mpf(0)]*5
    return out,dict(passed=True,native_phi= [record(v) for v in phi],
        native_U=[record(v) for v in U],native_U_y_absolute=[record(v) for v in Uycaps],
        native_U_minus_core_inlet_absolute=[record(v) for v in errors],
        native_phi_positive_lower=record(c.mpf(floor)),
        exact_source="phi_y=-chi*Dbar*phi/2; U_y=-chi*(phi/barphi)*R*Q",
        ordinary_Z0_Z4=True,physical_y0_y1=True,original_core_U_inlet_and_all_corrections_retained=True,
        native_and_comparison_never_identified=True,physical_dy_hb_ds_once=True)


def p2_numerator(c,values,p0,scalars):
    a=Dual(c,3);zero=[c.mpf(0)]*5
    z,d,L,g=[a.source(values,k) for k in ("Z","d","L","gaxis")]
    U,M,A,B,C=[a.source(values,k) for k in ("U","M","A","B","C")]
    Mz,Az,Bz,Cz=[a.shift(values,k) for k in ("M","A","B","C")]
    de=read(c,scalars["delta_record"])
    T=a.add(a.scale(a.mul(z,M),1-de),a.mul(d,Mz))
    lin=a.add(a.neg(U),a.scale(a.add(M,a.neg(a.mul(z,Mz))),(1-de)/2))
    p,pz=a.pair(p0,zero),a.pair(p0[1:],zero)
    Pr=a.add(a.scale(a.mul(z,p),2*(1+de)),a.neg(a.mul(d,pz)))
    S0=a.add(a.scale(a.mul(z,B),-2*de),a.mul(d,Bz),a.scale(a.mul(z,C),2*(1+de)),a.neg(a.mul(d,Cz)))
    rho=a.pair([scalars["rho"]]+zero,[scalars["rho"]]+zero)
    R=a.pair([scalars["R"]]+zero,[scalars["R"]]+zero)
    F2=a.pair(scalars["F2"],zero)
    centrifugal=a.mul(F2,a.add(a.mul(R,S0),a.scale(a.mul(a.mul(a.mul(rho,g),d),a.add(B,a.neg(C))),-2)))
    N=a.add(lin,a.mul(U,T),a.scale(a.mul(z,A),2*de),a.neg(a.mul(d,Az)),a.scale(Pr,scalars["P2"]),centrifugal)
    return a.div(N,a.mul(a.source(values,"phi"),L))


def log_upper(c,v):
    cap=c.fabs(v).b
    require(cap>0,"Positive norm envelope needed for logarithmic cap")
    return c.ln(cap)


def logged_products(c,jet,beta,offset):
    norm=bound_product(c,jet,beta[:len(jet)])
    return [offset+log_upper(c,v) for v in norm]


def logadd(c,left,right):
    return c.mpf(max(left.b,right.b))+c.ln(2)


def root_caps(c,values,D,Q,nat,p0,scalars,parent):
    a=Dual(c,3);zero=[c.mpf(0)]*5
    C=a.div(a.pair(Q[0],Q[1]),a.source(values,"comparison_Phi"))
    djet=a.pair(D[0],D[1]);TC=a.div(C,djet)
    PC=p2_numerator(c,nat,p0,scalars)
    inv=scalars["inv_relative"][:4]
    base=-scalars["logF"].a+c.ln(2*scalars["R"].b)/2
    pbase=-scalars["logF"].a+c.ln(scalars["R"].b/2)/2
    loghb=read(c,parent["actual_signed_lower_m15"]["source_log_bases"][0])
    bZ=logged_products(c,C[0],inv,base)
    bslow=jet_add(c,C[1],[v/2 for v in C[0]])
    bYregular=logged_products(c,bslow,inv,base)
    bYfirst=[logadd(c,z+c.ln(32)-loghb,y) for z,y in zip(bZ,bYregular)]
    bYsecond=[v+loghb for v in bYregular]
    tZ=logged_products(c,TC[0],inv,base)
    tY=logged_products(c,jet_add(c,TC[1],[v/2 for v in TC[0]]),inv,base)
    pZ=logged_products(c,PC[0],inv,pbase)
    pY=logged_products(c,jet_add(c,PC[1],[v/2 for v in PC[0]]),inv,pbase)
    coarse=dict(b_Z_first=(bZ,2),b_yZ_first=(bYfirst,102),
                t0_Z_both=(tZ,2),t0_yZ_both=(tY,2),p2_Z_both=(pZ,2),p2_yZ_both=(pY,2),
                b_Z_second=([v+loghb for v in bZ],-98),b_yZ_second=(bYsecond,-98))
    require(all(all(v.b < (factor*scalars["logC"]).a for v in jet) for jet,factor in coarse.values()),
            "Actual mixed logarithmic caps fail the stated coarse logC envelopes")
    def rows(jet):return [dict(ordinary_Z_order=k,log_absolute_upper_bound=record(v)) for k,v in enumerate(jet)]
    return dict(passed=True,ordinary_Z_orders=[0,1,2,3],physical_y_orders=[0,1],
        b_Z_first=rows(bZ),b_Z_second=rows([v+loghb for v in bZ]),b_yZ_first=rows(bYfirst),b_yZ_second=rows(bYsecond),
        t0_Z_both=rows(tZ),t0_yZ_both=rows(tY),p2_Z_both=rows(pZ),p2_yZ_both=rows(pY),
        log_caps_are_not_signed_values=True,F0_and_hb_and_N_not_materialized=True,
        b_exact_phi_cancellation=True,t0_exact_chi_cancellation=True,
        p2_uses_original_native_histories=True,p2_not_identified_with_b=True,
        p2_physical_sqrt_R_over_two_and_pressure_Pstar_squared_retained=True,
        coarse_logC_envelopes=dict(passed=True,logC=record(scalars["logC"]),
            strict_log_absolute_upper_multipliers={name:factor for name,(_,factor) in coarse.items()},
            each_multiplier_applies_to_ordinary_Z0_Z3=True,not_actual_signed_values=True),
        complete_b_t0_p2_mixed_provider=True,q_mixed_provider_complete=False,
        actual_finite_N_derivative_errors=False,numerical_Cy_Z3=False,numerical_QZZ=False)



def algebra_controls():
    """Exact checks of the two source reductions, before any boxing."""
    import sympy as s
    R,F,P,phi,barphi,chi,L=s.symbols("R F P phi barphi chi L",positive=True)
    rho,lam,g,d,de,z,U,M,A,B,C,Mz,Az,Bz,Cz,p,pz,W=s.symbols("rho lam g d de z U M A B C Mz Az Bz Cz p pz W")
    f=-lam*g
    S0=-2*de*z*B+d*Bz+2*(1+de)*z*C-d*Cz
    hydro=(-W*U+(1-de)*(M-z*Mz)/2+2*de*z*A-d*Az)/(2*L)
    Pr=2*(1+de)*z*p-d*pz
    pressure=Pr/(2*L)
    rawQ=hydro+P**2*pressure+R*F**2*(S0+2*f*d*(B-C))/(2*L)
    Q=hydro+P**2*pressure+F**2*(R*S0/(2*L)-rho*g*d*(B-C)/L)
    E=F*s.sqrt(2*R)*phi/P
    Uy=-chi*(phi/barphi)*R*Q
    b=2*Uy/(P*E); a=chi*s.symbols("Dbar",positive=True)
    bExpected=-s.sqrt(2*R)*chi*Q/(F*barphi)
    tExpected=s.sqrt(2*R)*Q/(F*barphi*s.symbols("Dbar",positive=True))
    transport=((1-de)*z*M+d*Mz)/P
    e=(A-R*F**2*B)/P**2;ez=(Az-R*F**2*(Bz+2*f*B))/P**2
    pressureField=p+R*F**2*C/P**2
    pressureZ=pz+R*F**2*(Cz+2*f*C)/P**2
    linear=(-U/P+(1-de)*(M/P-z*Mz/P)/2)/L
    quadratic=((U/P)*transport+2*de*z*e-d*ez+2*(1+de)*z*pressureField-d*pressureZ)/L
    rawP2=R*(linear+P*quadratic)/E
    N=-U+(1-de)*(M-z*Mz)/2+U*((1-de)*z*M+d*Mz)+2*de*z*A-d*Az+P**2*Pr+F**2*(R*S0-2*rho*g*d*(B-C))
    reducedP2=s.sqrt(R/2)*N/(F*phi*L)
    checks=dict(comparison_R_Lambda_cancellation=s.simplify(rawQ.subs(lam,rho/R)-Q)==0,
        actual_b_phi_and_Pstar_cancel=s.simplify(b-bExpected)==0,
        actual_t0_chi_cancel=s.simplify(-b/a-tExpected)==0,
        native_generic_p2_reduction=s.simplify(rawP2.subs(lam,rho/R)-reducedP2)==0)
    require(all(checks.values()),"Actual comparison/generic scale identity failed")
    return dict(passed=True,checks=checks,signed_sources_combined_before_interval_boxing=True,
        generic_p2_distinct_from_comparison_pressure_drive=True)


def build(root,parent):
    root=Path(root)
    require(sha(root/PARENT)==PARENT_SHA256 and load(root/PARENT)==parent,"Accepted E/a parent required")
    for path,h in REVIEWED_HELPERS.items():
        require(sha(root/path)==h,"Changed reviewed new source helper: "+path)
    from mpmath.ctx_iv import MPIntervalContext
    from .original_micro_source_Z5 import extend
    from .original_pressure_axial_tube import build as pressure_build
    c=MPIntervalContext();c.prec=600
    values,D,extension,pins=extend(c,root,parent)
    p0,pressure,pressurePins=pressure_build(c,root,parent,order=5);pins.update(pressurePins)
    for path in (PARENT,Z5,PRESSURE,TRANSFER,TUBE,NORMS,CORE_CHECK,BRIDGE+".py",INLET+".py",CORE):
        expected=PARENT_SHA256 if path==PARENT else sha(root/path)
        if path.startswith("experiments/root_st073/") and path!=PARENT:
            require(expected==parent["input_hashes"][Path(path).name],"Changed accepted analytic/source input")
        pins[path]=dict(path=path,expected_sha256=expected)
    major,tube,norms=[load(root/path) for path in (TRANSFER,TUBE,NORMS)]
    scalars,amp=amplitude(c,root,parent,values,major,tube,norms)
    scalars["delta_record"]=parent["actual_core_signed_direction"]["parameters"]["delta"]
    globalPhi=[read(c,load(root/(PREFIX+"compliant_core_uniform_bounds.json"))["actual_Phi_upper"]).b]+absolute(c,values["core_Phi"][1:])
    inlets=complete_inlet_norms(c,values,globalPhi)
    add_histories(c,values,"comparison_Phi","comparison_U",inlets)
    Q=drive(c,values,p0,scalars)
    nat,nativeProof=native(c,values,D,Q,scalars,inlets["M"])
    add_histories(c,nat,"phi","U",inlets)
    caps=root_caps(c,values,D,Q,nat,p0,scalars,parent)
    return dict(passed=True,source_family=parent["source_family"],source_pins=pins,
        source_extension=extension,whole_tube_pressure=pressure,whole_tube_amplitude=amp,
        reduced_source_algebra=algebra_controls(),
        native_angular_and_axial_provider=nativeProof,actual_b_t0_p2_mixed_bounds=caps,
        domain=dict(s=[0,2],Z=["37/100","93/250"],physical_y="log(R/Ra)=hb*s"),
        comparison_complete_ABC_history_rows={k:[record(v) for v in values[k]] for k in ("A","B","C","A_y","B_y","C_y")},
        native_complete_history_rows={k:[record(v) for v in nat[k]] for k in ("M","A","B","C","M_y","A_y","B_y","C_y")},
        comparison_reduced_drive_Z0_Z4=[record(v) for v in Q[0]],
        comparison_reduced_drive_yZ0_Z4=[record(v) for v in Q[1]],
        source_AST=source_guards(root),ordinary_derivatives_and_original_inlets=True,
        complete_inlet_integral_bound=dict(passed=True,radial_domain=["0","41/10"],
            inlet_substitution="rho=4*t, 0<=t<=1; independent of Z",
            U_input="4Z+j plus global Xh epsilon*Psi errors; NOT a point-only U enclosure",
            Phi_input="global Phi0 ceiling; global Xh ordinary derivatives at orders1..5",
            inlet_M="integral_0^1 U(4t,Z) dt; bound by global U ordinary norms",
            inlet_A="integral_0^1 U(4t,Z)^2 dt; binomial product of global U norms",
            inlet_B="integral_0^1 t*Phi(4t,Z)^2 dt; global product norm times 1/2",
            inlet_C="integral_0^1 Phi(4t,Z)^2 dt; global Phi product norm",
            exact_inlet_weights=dict(M=1,A=1,B="1/2",C=1),
            complete_inlet_norm_rows={name:[record(v) for v in row] for name,row in inlets.items()},
            same_complete_inlets_passed_to_native_and_comparison=True,
            global_U_norms=[record(v) for v in absolute(c,values["core_U"])],
            global_Phi_norms=[record(v) for v in globalPhi],
            history_rule="max(initial_derivative_cap, differentiated_source_cap/rate)",
            incoming_functions_retained_not_replaced_by_caps=True),
        exact_reduced_functions=dict(
            b="-sqrt(2R)*chi_R*Q/(F0*barphi)",
            t0="sqrt(2R)*Q/(F0*barphi*Dbar)",
            p2="sqrt(R/2)*N/(F0*phi*L)",
            Q="hydro+Pstar^2*pressure+F0^2*(R*S0/(2L)-rho*gaxis*d*(JB-JC)/L), comparison histories",
            N="-U+(1-delta)*(JM-Z*JM_Z)/2+U*((1-delta)*Z*JM+d*JM_Z)+2delta*Z*JA-d*JA_Z+Pstar^2*Pr+F0^2*(R*S0-2rho*gaxis*d*(JB-JC)), native histories",
            S0="-2delta*Z*JB+d*JB_Z+2(1+delta)*Z*JC-d*JC_Z",
            Pr="2(1+delta)*Z*p0-d*p0_Z"),
        actual_source_handles=dict(
            first_roots=parent["actual_micro_direction_C0"]["actual_micro_source_handles"]["first"]["roots"],
            second_roots=parent["actual_micro_direction_C0"]["actual_micro_source_handles"]["second"]["roots"],
            native=parent["actual_micro_direction_C0"]["actual_micro_source_handles"]["native"]),
        physical_y_jets_derived_from_exact_defining_ODEs=True,
        no_invented_p2_y_graph_handle=True,
        complete_mixed_root_provider=False,q_mixed_provider_complete=False,
        numerical_Cy_Z3=False,numerical_mean_Z3=False,numerical_QZZ=False,
        global_five_moment_closure=False,full_goal_complete=False,
        inherited_signed_sources_and_all_errors_unchanged=True,accepted_ancestors_replayed=False)


def source_guards(root):
    return dict(complete_inlets=source_guard(root,INLET+".py","__init__",(
        "sample=g.mul(rho,t)","A=a.mul(u,u),B=a.scale(a.mul(phi,phi),t),C=a.mul(phi,phi)",
        "self.moment_functions={key:a.integral(body,t,g.zero,g.one,actual_core_normalized_moment=key,",
        "self.inlet_moments={key:a.at(jet,rho,c(4)) for key,jet in self.moment_functions.items()}")),
        global_U_caps=source_guard(root,"src/openai_ns_reconstruction/original_micro_Ea_mixed_bounds.py","source_jets",(
        '4*z[0]+j+sym(cap("Uz_minus_U0_absolute_bound",0,0))',
        '4+sym(cap("Uz_minus_U0_absolute_bound",0,1))')),
        complete_histories=source_guard(root,BRIDGE+".py","densities",(
        "A=a.mul(U,U),B=a.mul(phi,phi),C=a.mul(phi,phi)")),
        comparison_drives=source_guard(root,BRIDGE+".py","comparison_directions",(
        "a.scale(a.mul(a.add(one,de),a.mul(z,a.trim(P0))),2)",
        "swirl=a.div(dressed_swirl,a.trim(F2)")),
        actual_native=source_guard(root,BRIDGE+".py","actual_fields",(
        "Vy=a.scale(a.mul(ratio,full_drive),g.neg(chi))",
        "actual_prescribed_axial_equation=True,all_hydro_pressure_swirl_sectors_retained=True")),
        generic_p2=source_guard(root,INLET+".py","generic_functions",(
        "p2=a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,",
        "a.scale(a.mul(a.add(one,de),a.mul(z,pressure)),2)")))
