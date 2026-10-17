"""Correlated leading-pressure increments from the actual mean FTC.

The shared inlet is one function. Psi_rho is transported from that inlet,
with the original nonzero prescribed-shear sigma retained as a fifth factor.
No independently boxed endpoint Qpsi values are subtracted here.
"""
from __future__ import annotations
from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
import sympy as sp

from .original_finite_N_signed_source import ast_guard,digest,load,require,sha
from .original_n1_selected_word_numeric import interval,endpoint,record
from .original_n1_selected_pressure_primitive import add,scale,shift,product,poly,export,jadd,jmul,jscale,jdiv

BASE="experiments/root_st073/"
PREFIX=BASE+"lei_ren_part1_paper_compliant_current_original_n1_"
PARENT=PREFIX+"selected_pressure_primitive.json"
PARENT_SHA256="0a8131af5c9a331affbe3f18db0a8f525876f80666e205fbaeda0d2b8fe86a21"
MEAN=PREFIX+"leading_mean_epsilon_factor.json.gz"
MEAN_CHECK=PREFIX+"leading_mean_epsilon_factor_check.json"
MEAN_SOURCE=PREFIX+"leading_mean_epsilon_factor.py"
DOMAIN=PREFIX+"inner_analytic_domain.json.gz"
MAJOR=BASE+"lei_ren_part1_paper_compliant_core_transfer.json"
MAJOR_CHECK=BASE+"lei_ren_part1_paper_compliant_core_transfer_check.json"
PHYSICAL_SOURCE=BASE+"lei_ren_part1_paper_compliant_core_physical_field.py"
FLAT_SOURCE=BASE+"lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
NATIVE_SOURCE="src/openai_ns_reconstruction/original_native_collar_transport.py"
ZERO=(0,0,0,0,0)
H=(1,0,0,0,0)
SIGMA=(0,0,0,0,1)
PRECISION_BITS=1800
SMALL_RATIO_BITS=1024


def lift(c,row,order=0,upper=False):
    return {k+(0,):v for k,v in poly(c,row,order,positive_upper=upper).items()}


def hydrate(c,rows):
    return {tuple(t["source_exponents"])+(0,):interval(c,t["coefficient_interval"]) for t in rows}


def scalar(c,v):
    return {ZERO:c.mpf(v)}


def jet_linear(coeff,p):
    return [scale(p[0],coeff[0]),add(scale(p[0],coeff[1]),scale(p[1],coeff[0]))]


def model_jets(c,sample,delta):
    z=[interval(c,sample["Z"]),c.mpf(1)]
    one=[c.mpf(1),c.mpf(0)]
    d=jadd(one,jscale(jmul(z,z),-1));L=jadd(one,jscale(jmul(z,z),-delta))
    p=[lift(c,sample["leading_Psi_rho_rows"][1],i) for i in (0,1,2)]
    A=jdiv(jscale(z,(3+delta)/2),L);B=jdiv(jscale(d,-c.mpf(1)/2),L)
    K=jet_linear(A,p[:2]);other=jet_linear(B,p[1:3])
    K=[add(a,b) for a,b in zip(K,other)]
    V=jet_linear(jdiv(jscale(z,2),L),p[:2])
    return K,V


def flat_atom(c,s):
    require(s>0 and s<c.mpf(1)/4096,"Selected positive first phase required")
    odds=1/(1-s)**2-1/s**2
    cap=s*c.exp(-1000)
    require(odds.b<(c.log(s)-1000).a,"Flat logistic cap must follow from exact logit")
    kappa=c.make_mpf(((1/(1+cap))._mpi_[0],c.mpf(1)._mpi_[1]))
    log_bounds=c.make_mpf(((odds-cap)._mpi_[0],odds._mpi_[1]))
    x=sp.symbols("x");y,E=sp.symbols("y E",positive=True)
    L1=2/(1-x)**3+2/x**3
    derivative_identity=sp.simplify(sp.diff(1/(1-x)**2-1/x**2,x)-L1)==0
    # x=y/(1+y) parametrizes (0,1); E=exp(odds)>0. Both the
    # original logistic factor and its logit derivative are strictly positive.
    positive_derivative=2*(1+y)**3*(1+y**(-3))*E/(1+E)**2
    positive_rewrite=sp.cancel(L1.subs(x,y/(1+y))*E/(1+E)**2-positive_derivative)==0
    require(derivative_identity and positive_rewrite and positive_derivative.is_positive is True,
        "Original sigma must be monotone before bounding its integrals")
    return dict(name="sigma_selected",argument=record(s),
        exact_function="sigma(s)=exp(odds)/(1+exp(odds)), odds=(1-s)^(-2)-s^(-2)",
        odds=record(odds),kappa=record(kappa),log_sigma_bounds=record(log_bounds),
        positive_cap=record(cap),cap_is_not_the_function_value=True,
        sigma_is_strictly_positive=True,sigma_materialized=False,
        monotone_on_selected_real_interval_certified=True,
        monotonicity_proof=dict(logit_derivative_identity=True,positive_rewrite_identity=True,
            derivative="sigma*(1-sigma)*(2/(1-x)^3+2/x^3)",
            positivity_substitution="x=y/(1+y), y>0; E=exp(odds)>0",
            transformed_derivative=str(positive_derivative),strictly_positive=True,
            flat_left_endpoint="sigma(0)=0 from pinned source branch",
            consequence="0<=sigma(t)<=sigma(s) for 0<=t<=s"),
        ordinary_Z_derivative_is_zero=True,
        distinct_from_own_velocity_cutoff_eta_at_q1=True)


def correlated_error(c,s,Cd,Cw,mean,domain,major):
    proof=mean["uniform_complex_bound_proof"]
    eta=endpoint(c,interval(c,proof["original_outer_eta"]),0)
    Ac=endpoint(c,interval(c,proof["actual_core_logarithmic_derivative_modulus_upper"]))
    ratio=interval(c,proof["actual_Phi_over_core_modulus_upper"])
    require(ratio.a>=1,"Saved uniform Psi_rho bound must dominate the core bound")
    B1=lift(c,mean["uniform_complex_factored_bounds"]["Psi"][1],upper=True)
    require(all(k[0]==0 for k in B1),"Core first-derivative comparison bound must have no hb pole")
    size=endpoint(c,interval(c,major["scaled_map_size_upper"]))
    q=endpoint(c,interval(c,major["scaled_map_Lipschitz_upper"]))
    require(size>0 and 0<=q<1,"Accepted nonlinear core contraction required")
    # Whole-core nonlinear value <=4*C, rho Cauchy radius1/2 gives
    # ||Psi_core,rhorho|| <= 2!*4*C/(1/2)^2 =32*C.
    C=size/(1-q);B2=scalar(c,32*C)
    rho=endpoint(c,interval(c,domain["actual_inner_geometry"]["real_rho_upper_bound"]))
    L=endpoint(c,interval(c,domain["denominator_bounds"]["L_modulus_lower"]),0)
    require(L>0,"Positive common complex denominator required")
    M=s*rho*Ac/2
    require(M.b<c.mpf(1)/2,"Complex ratio exponential must stay below the factor2 cover")
    # R(0)=1, R_s=hb*rho*(chi-1)*Acore*R:
    # |R-1|<=2*hb*s*rho_max*Ac*sigma(s).
    # Use |R|<=2 on the core radial change and |R-1| on pA.
    T=add(scale(B2,2*Cd),shift(scale(B1,2*s*rho*Ac),SIGMA))
    end_error=shift(T,H)
    # Epsi=DeltaPsi-pA*DeltaRho; EJ=int rho*Psi_rho-pA*int rho.
    Epsi=add(shift(scale(B1,Cd),(1,0,0,0,1)),shift(scale(T,Cd),(2,0,0,0,0)))
    EJ=add(shift(scale(B1,Cw),(1,0,0,0,1)),shift(scale(T,Cw),(2,0,0,0,0)))
    loss,gap=eta/8,eta/16
    z=1+eta;d=1+z*z
    # k=1-delta and 1+delta both have modulus <=1+delta_upper<2.
    common=2*z+d/loss
    P=scale(add(scale(end_error,2*z*rho),scale(Epsi,common)),1/L)
    J=scale(add(scale(end_error,2*z*rho*rho),scale(EJ,common)),1/L)
    return shift(P,(0,0,0,-1,0)),shift(J,(0,0,0,-2,0)),gap,dict(
        core_correction_bound=record(C),core_Psi_rho2_bound=record(32*C),
        ratio_exponent_cover=record(M),core_log_derivative_upper=record(Ac),
        source_ratio="R(0)=1; R_s=hb*rho*(chi-1)*(Phi_core_rho/Phi_core)*R",
        ratio_error="|R-1|<=2*hb*s*rho_max*Ac*sigma_selected",
        derivative_error="|Psi_rho-chi*pA|<=hb*(2*Cd*B2+2*s*rho_max*Ac*B1*sigma_selected)",
        Epsi="DeltaPsi-pA*DeltaRho; <=hb*Cd*sigma*B1+hb^2*Cd*T",
        EJ="int rho*Psi_rho-pA*int rho; <=hb*Cw*sigma*B1+hb^2*Cw*T",
        first_Z_loss=record(loss),whole_error_Z1_gap=record(gap),
        no_independent_inlet_error_subtraction=True)


def normalized_enclosure(c,signed,errors,bases,s,eps_power):
    dominant=(1,0,1,-eps_power,0)
    require(dominant in signed,"Expected physical pressure-sector term missing")
    D=signed[dominant]/s
    require(D.a>0 or D.b<0,"Dominant coefficient must have a strict sign")
    threshold=-SMALL_RATIO_BITS*c.log(2)
    tiny=c.mpf(1)/(1<<SMALL_RATIO_BITS)
    radius=c.mpf(0);proof=[]
    terms=[("other_signed",k,v) for k,v in signed.items() if k!=dominant]
    terms += [("error_radius",k,v) for k,v in errors.items()]
    for kind,key,value in terms:
        bound=endpoint(c,abs(value))
        if bound==0:continue
        logratio=c.log(bound)-c.log(s)
        diff=tuple(a-b for a,b in zip(key,dominant))
        # Identical source factors cancel BEFORE boxing; sigma is a true
        # retained source atom, not its positive vanishing cap.
        for power,base in zip(diff,bases):
            if power:logratio+=power*base
        require(logratio.b<=threshold.a,"A normalized term is not within the certified dyadic tail")
        radius+=tiny
        proof.append(dict(kind=kind,source_exponents=list(key),factor_difference=list(diff),
            log_absolute_normalized_upper=record(logratio),dyadic_upper="2^-"+str(SMALL_RATIO_BITS)))
    result=c.make_mpf(((D-radius)._mpi_[0],(D+radius)._mpi_[1]))
    strict="positive" if result.a>0 else "negative" if result.b<0 else "undetermined"
    require(strict!="undetermined","Computed local leading pressure sign unresolved")
    relative_width=(result.b-result.a)/min(abs(result.a),abs(result.b))
    return dict(normalizing_positive_factor="hb*Pstar^2*epsilon^"+str(eps_power)+"*s",
        dominant_source_exponents=list(dominant),dominant_coefficient_div_phase=record(D),
        all_other_terms_absolute_upper=record(radius),dyadic_term_bounds=proof,
        normalized_integral_enclosure=record(result),strict_sign=strict,
        relative_interval_width_upper=record(relative_width),
        retained_flat_source_not_zeroed=True,dyadic_caps_are_errors_not_signed_values=True)


def exact_algebra():
    r,z=sp.symbols("rho Z",real=True);de=sp.symbols("delta")
    f,m=sp.Function("Psi")(r,z),sp.Function("m")(r,z)
    L,d=1-de*z*z,1-z*z
    b=(2*z*f-(1-de)*z*m-d*sp.diff(m,z))/L
    meanr=(f-m)/r
    br=sp.diff(b,r).subs(sp.diff(m,r),meanr).subs(
        sp.Derivative(m,r,z),sp.diff(meanr,z))
    br=br.subs(sp.Derivative(m,z,r),sp.diff(meanr,z))
    P=(2*z*r*sp.diff(f,r)+(1+de)*z*f-d*sp.diff(f,z))/L
    checks={"mean_removed_from_P1_primitive":sp.simplify((b+r*br-P).doit())==0}
    primitive=r*(f-m)
    for k in (0,1):
        defect=sp.diff(primitive,r)-r*sp.diff(f,r)
        defect=defect.subs(sp.diff(m,r),meanr)
        checks["mean_removed_from_m15_difference_Z"+str(k)]=sp.simplify(sp.diff(defect,z,k))==0
    p,pz,dr,jw,rB,sg,h=sp.symbols("p pz DeltaR Jw rhoB Sigma hb")
    K=((3+de)*z*p-d*pz)/(2*L)
    regP=(2*z*dr*p+(1+de)*z*dr*p-d*dr*pz)/L
    regM=-(2*z*(2*jw)*p-(1-de)*z*jw*p-d*jw*pz)/L
    checks["P1_regular_factor_two"]=sp.simplify(regP-2*K*dr)==0
    checks["m15_regular_factor_minus_two"]=sp.simplify(regM+2*K*jw)==0
    wb=p*(1-(1-h)*sg);distance=rB-4;weight=(rB*rB-16)/2
    Pmodel=(2*z*(rB*wb-4*p)+(1+de)*z*p*distance-d*pz*distance)/L
    Mmodel=-(2*z*(rB*rB*wb-16*p)-(1-de)*z*p*weight-d*pz*weight)/L
    Pflat=-2*z*rB*(1-h)*sg*p/L;Mflat=2*z*rB*rB*(1-h)*sg*p/L
    checks["P1_flat_term_sign"]=sp.cancel(sp.expand(Pmodel-2*K*distance-Pflat))==0
    checks["m15_flat_term_sign"]=sp.cancel(sp.expand(Mmodel+2*K*weight-Mflat))==0
    require(all(checks.values()),"Correlated own-mean primitive algebra failed")
    return dict(passed=True,checks=checks)


def build(root,parent=None):
    root=Path(root);require(sha(root/PARENT)==PARENT_SHA256,"Accepted pressure parent changed")
    saved=load(root/PARENT)
    if parent is None:parent=saved
    require(digest(parent)==digest(saved) and parent["checks"]["passed"] is True,"Same passed pressure parent required")
    paths=(PARENT,MEAN,MEAN_CHECK,MEAN_SOURCE,DOMAIN,MAJOR,MAJOR_CHECK,PHYSICAL_SOURCE,FLAT_SOURCE,NATIVE_SOURCE)
    pins={}
    for path in paths:
        expected=(PARENT_SHA256 if path==PARENT else parent["input_hashes"][path[len(BASE):]]
            if path.startswith(BASE) else parent["implementation_hashes"][path])
        require(sha(root/path)==expected,"Changed correlated pressure source: "+path)
        pins[path]=dict(path=path,expected_sha256=expected,passed=True)
    data,check,domain,major,mc=load(root/MEAN),load(root/MEAN_CHECK),load(root/DOMAIN),load(root/MAJOR),load(root/MAJOR_CHECK)
    require(check["all_passed"] is True and mc["all_passed"] is True and major["contraction_proved"] is True and
        major["implicit_source_sha256"]==parent["source_family"]["implicit_source_sha256"] and
        major["datum_enclosure_sha256"]==parent["source_family"]["datum_enclosure_sha256"] and
        data["source_family"]==domain["source_family"]==parent["source_family"] and
        check["input_hashes"][Path(MEAN).name]==pins[MEAN]["expected_sha256"],
        "Same actual core contraction and checked mean source required")
    native=parent["actual_native_collar_transport"]
    require(native["passed"] and native["actual_quotient_equals_ratio"] is True and
        native["native_defining_leading_function_equals_canonical_collar"] is True and
        native["first_window_only"] is True and
        native["domain"]["ordinary_Z_independent_endpoints"] is True and
        "0<hb<1/2" in native["domain"]["rho_bound_proof"],
        "Actual native/core comparison ratio identity required")
    old=parent["actual_selected_leading_pressure_primitive"]
    A,B=data["samples"]["actual_first_inlet"],data["samples"]["selected_actual_inner_exit"]
    require(A["Z"]==B["Z"]==old["selected_Z_interval"] and
        A["own_mean_identity"]==B["own_mean_identity"]==
        "rho*meanPsi=4*meanPsi_core(4)+integral_4^rho Psi_actual drho" and
        old["physical_N1p_no_Pstar_normalization_used"] is True,
        "Same own-mean source and physical N1p required")
    guards=dict(
        mean=ast_guard(root,MEAN_SOURCE,"actual_rows",(
            "mean1=op.truncate((psi0-mean0)/rho,4)","mean2=op.truncate((psi1-2*mean1)/rho,3)")),
        core_bounds=ast_guard(root,MEAN_SOURCE,"complex_mean_bounds",(
            "correction=domain.upper(c,4*physical.correction)",
            "core_psi=[slope*r+lift(correction),slope+lift(2*correction),lift(8*correction)]",
            "psi1=core_psi[1]*(ratio*(1+cap))")),
        correction=ast_guard(root,PHYSICAL_SOURCE,"__init__",(
            "q=read(major,'scaled_map_Lipschitz_upper')",
            "self.correction=read(major,'scaled_map_size_upper')/(1-q)")),
        native=ast_guard(root,NATIVE_SOURCE,"_algebra",(
            "ar=(chi-1)*h*r*sp.diff(p,r)/p*a",
            "psi_s=a*chi*h*r*sp.diff(v,r)","mean_s=h*(Psi-mean)")),
        flat=ast_guard(root,FLAT_SOURCE,"_sigma_left",(
            "odds=1/(1-x)**2-1/x**2","e=positive_exp(c,odds);value=e/(1+e)",
            "L1=2/(1-x)**3+2/x**3","q=value*(1-value)",
            "derivatives=[value,q*L1", "if hi<=0:return IntervalTaylor.constant(c,0,4)")),
    )
    c=MPIntervalContext();c.prec=PRECISION_BITS
    s=interval(c,old["selected_phase"]);delta=interval(c,old["delta_enclosure"])
    require(s==interval(c,B["geometry"]["phase"]) and delta>0 and delta<c.mpf("1e-200"),
        "Same exact phase and delta required")
    Cd=interval(c,old["radius_factors"]["rhoB_minus4_over_hb"])
    Cw=interval(c,old["radius_factors"]["integral_rho_drho_over_hb"])
    atom=flat_atom(c,s)
    pError,mError,gap,error_proof=correlated_error(c,s,Cd,Cw,B,domain,major)
    require(gap==interval(c,old["radius_factors"]["complex_S1_S2_derivative_gap"]),
        "Same source solution Cauchy gap required")
    K,V=model_jets(c,A,delta)
    one_minus_h=add(scalar(c,1),shift(scalar(c,-1),H))
    rhoB=add(scalar(c,4),shift(scalar(c,Cd),H))
    rhoB2=add(scalar(c,16),shift(scalar(c,2*Cw),H))
    bases=[interval(c,v) for v in old["source_log_bases"]]+[interval(c,atom["log_sigma_bounds"])]
    rows={}
    for k in (0,1):
        P_reg=shift(scale(K[k],2*Cd),(1,0,0,-1,0))
        M_reg=shift(scale(K[k],-2*Cw),(1,0,0,-2,0))
        P_flat=shift(scale(product(product(one_minus_h,rhoB),V[k]),-1),(0,0,0,-1,1))
        M_flat=shift(product(product(one_minus_h,rhoB2),V[k]),(0,0,0,-2,1))
        old_row=old["ordinary_Z"][str(k)]
        P_axis=hydrate(c,old_row["P1_W0_local"]["axis_inertial"])
        M_axis=hydrate(c,old_row["m15_leading_N1p_local"]["axis_inertial"])
        P_signed=add(P_reg,P_flat,P_axis);M_signed=add(M_reg,M_flat,M_axis)
        factor=1 if k==0 else 1/gap
        P_errors=add(scale(pError,factor),hydrate(c,old_row["P1_W0_local"]["remaining_inertial_error_radius"]))
        M_errors=add(scale(mError,factor),hydrate(c,old_row["m15_leading_N1p_local"]["remaining_inertial_error_radius"]))
        require(all(len(key)==5 and key[1]==0 for p in (P_signed,M_signed,P_errors,M_errors) for key in p),
            "Z-dependent F0 factors cannot be differentiated as constant scales")
        rows[str(k)]=dict(
            P1_W0_local=dict(regular_viscous=export(P_reg),flat_viscous=export(P_flat),
                axis_inertial=export(P_axis),signed_approximation=export(P_signed),
                correlated_error_radius=export(P_errors,error=True),
                normalized=normalized_enclosure(c,P_signed,P_errors,bases,s,1)),
            m15_leading_N1p_local=dict(regular_viscous=export(M_reg),flat_viscous=export(M_flat),
                axis_inertial=export(M_axis),signed_approximation=export(M_signed),
                correlated_error_radius=export(M_errors,error=True),
                normalized=normalized_enclosure(c,M_signed,M_errors,bases,s,2)))
    algebra=exact_algebra()
    return dict(passed=True,source_family=parent["source_family"],parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,parent_mapping_sha256=digest(parent),
        source_pins=pins,source_AST=guards,exact_algebra=algebra,ordinary_Z_orders=[0,1],
        selected_phase=old["selected_phase"],selected_Z_interval=A["Z"],ordinary_Z=rows,
        source_log_basis_names=old["source_log_basis_names"]+["sigma_selected"],
        first_four_source_log_bases=old["source_log_bases"],flat_source_atom=atom,
        error_proof=error_proof,
        shared_mean_and_Psi_value_inlet_constants_canceled_by_same_source_FTC=True,
        single_core_Psi_rho_inlet_uncertainty_retained=True,
        independent_endpoint_uncertainty_subtraction_eliminated_in_new_representation=True,
        actual_selected_leading_N1p_increment_sign_certified=True,
        actual_correlated_signed_coefficients_and_error_radii_evaluated=True,
        old_independent_endpoint_Qpsi_values_not_used=True,
        raw_first_keep_s_half_values_not_used=True,
        original_nonzero_flat_source_factor_retained=True,
        P1_is_distinct_from_m15=True,full_n1_pressure_or_C_inner_m15_evaluated=False,
        all_Z_functional_coverage=False,global_five_moment_closure=False,
        old_error_caps_or_six_ledger_terms_removed=False,
        accepted_ancestors_replayed=False,source_executed=False,source_only=False,
        precision_bits=PRECISION_BITS)
