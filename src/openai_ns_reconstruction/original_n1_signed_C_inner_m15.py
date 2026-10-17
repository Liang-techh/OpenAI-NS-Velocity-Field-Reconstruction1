"""Actual core plus first-collar m15 at the saved axial interval.

The signed core pressure primitive is evaluated from one actual endpoint.
Other source integrals have explicit signed approximations and/or symmetric
error bounds. Norm caps are never substituted for signed function values.
"""
from __future__ import annotations

from pathlib import Path
from fractions import Fraction
import hashlib
import json
from mpmath.ctx_iv import MPIntervalContext
import sympy as sp

from .original_finite_N_signed_source import ast_guard,digest,load,require,sha
from .original_n1_selected_word_numeric import interval,endpoint,record
from .original_n1_selected_pressure_primitive import add,scale,shift,product,export,axis_jet,norm_polynomials
from .original_n1_correlated_pressure_increment import lift,scalar
from .original_n1_velocity_m15_increment import saved_polynomial

BASE="experiments/root_st073/"
PREFIX=BASE+"lei_ren_part1_paper_compliant_current_original_n1_"
PARENT=PREFIX+"complete_first_collar_m15.json"
PARENT_SHA256="333d9afeb4d16542b922aecb41101875eac12c3256f980bce1a13f1574a167bd"
MEAN=PREFIX+"leading_mean_epsilon_factor.json.gz"
MEAN_CHECK=PREFIX+"leading_mean_epsilon_factor_check.json"
MEAN_SOURCE=PREFIX+"leading_mean_epsilon_factor.py"
DOMAIN=PREFIX+"inner_analytic_domain.json.gz"
MAJOR=BASE+"lei_ren_part1_paper_compliant_core_transfer.json"
PHYSICAL_SOURCE=BASE+"lei_ren_part1_paper_compliant_core_physical_field.py"
RAW_SOURCE=PREFIX+"inner_volterra_functions.py"
COUPLED_SOURCE=PREFIX+"coupled_first_step.py"
ZERO=(0,0,0,0,0)
EPS=(0,0,0,-1,0)
PRECISION_BITS=1800
HELPERS=tuple("src/openai_ns_reconstruction/"+name+".py" for name in (
    "original_finite_N_signed_source","original_n1_selected_word_numeric",
    "original_n1_selected_pressure_primitive","original_n1_correlated_pressure_increment",
    "original_n1_velocity_m15_increment"))


def absolute(c,p):
    return {k:endpoint(c,abs(v)) for k,v in p.items()}


def axis_Nz(c,z,j,de):
    u=4*z+j;d=1-z*z;L=1-de*z*z
    N=4-(8+4*de)*z*z-(1+de)*j*z
    Nz=-2*(8+4*de)*z-(1+de)*j
    Nzz=-2*(8+4*de)
    Lz=-2*de*z;Lzz=-2*de
    A=N/L;Az=(Nz-Lz*A)/L;Azz=(Nzz-Lzz*A-2*Lz*Az)/L
    value=(2*z*A-d*Az)/L
    derivative=(2*A+4*z*Az-d*Azz-Lz*value)/L
    return [value,derivative],[u,c.mpf(4)]


def core_bounds(c,parent,mean,domain,delta,j):
    sample=mean["samples"]["selected_actual_inner_exit"]
    proof=sample["uniform_complex_bound_proof"]
    eta=endpoint(c,interval(c,proof["original_outer_eta"]),0)
    require(interval(c,domain["common_tubes"]["outer_Z_radius"])==eta and
        interval(c,proof["Psi_and_meanPsi_bound_tube_radius"])==eta/2,
        "Core and admitted native Psi bound must share the larger Z tube")
    require(interval(c,proof["actual_Phi_over_core_modulus_upper"]).a>=1,
        "Actual native Psi_rho cover must dominate the core comparison bound")
    B1=lift(c,sample["uniform_complex_factored_bounds"]["Psi"][1],upper=True)
    require(all(k[0]==k[1]==k[4]==0 for k in B1),"Core derivative cover must not omit width/amplitude factors")
    z=1+eta;d=1+z*z;de=endpoint(c,delta)
    L=endpoint(c,interval(c,domain["denominator_bounds"]["L_modulus_lower"]),0)
    require(L>0,"Nonzero core axial denominator required")
    u=4*z+endpoint(c,abs(j));a=((1+de)*z*u+4*d)/L
    # Axis Psi=0, so real-rho FTC and the averaging identity give
    # |Psi|<=4B1, |m|<=2B1, |m_rho|<=B1/2, without division by rho.
    qloss=eta/4;loss=eta/8;gap=eta/16
    b=scale(B1,(8*z+2*(1+de)*z+2*d/qloss)/L)
    br=scale(B1,(2*z+(1+de)*z/2+d/(2*qloss))/L)
    norms=dict(z=z,d=d,L=L,rho=c.mpf(4),a=a,uv=u,b=b,br=br,v=scale(B1,4),loss=loss,gap=gap)
    S1,S2=norm_polynomials(c,norms,delta)
    A=((1+de)*z*u+4*d)/L
    Az=((1+de)*u+4*(3+de)*z)/L+A*(2*de*z)/L
    Nza=(2*z*A+d*Az)/L
    return norms,S1,S2,Nza,dict(original_outer_eta=record(eta),core_radial_interval=[0,4],
        core_Psi_rho_cover=export(B1,error=True),Psi_and_mean_source_Z_radius=record(eta/2),
        mean_Psi_value_and_radial_derivative_bounds="|Psi|<=4B1, |m|<=2B1, |m_rho|<=B1/2 by axis FTC",
        mean_identity="m(rho)=integral_0^1 Psi(t*rho) dt; m_rho=integral_0^1 t*Psi_rho(t*rho) dt",
        Qpsi_source_Z_radius=record(eta/4),first_Z_loss=record(loss),whole_error_Z1_gap=record(gap),
        no_inverse_rho_bound_used=True,native_bound_used_only_as_core_derivative_upper=True)


def raw_core_bounds(c,parent,domain):
    source=parent["n1_complex_actual_source_bounds"]
    require(source["passed"] is True and source["actual_source_bounds_proved"] is True and
        source["source_family_match"] is True and
        source["domain"]["actual_collar_used_as_core_extension"] is False,
        "A proved genuine core source bound is required")
    sd=source["domain"];tubes=domain["common_tubes"]
    require(interval(c,sd["outer_Z_radius"])==interval(c,tubes["outer_Z_radius"]) and
        interval(c,sd["Gg_Z_radius"])==interval(c,tubes["Gg_Z_radius"]) and
        interval(c,sd["solution_Z_radius"])==interval(c,tubes["solution_Z_radius"]),
        "Core source and solution bounds must use the exact same complex Z tubes")
    require(interval(c,tubes["Gg_Z_radius"])==interval(c,tubes["outer_Z_radius"])/4 and
        interval(c,tubes["solution_Z_radius"])==interval(c,tubes["outer_Z_radius"])/8,
        "The raw value/error bounds must cover the eta/16 derivative disk")
    rad=Fraction(sd["rho_bound_radius"]);target=Fraction(sd["rho_target_radius"])
    require(rad>=4 and target==4 and sd["rho_source_disk"]=="|rho| < 5" and
        source["raw_ast_binding"]["passed"] is True and source["gates"]["source_normalization"] is True and
        source["weighted_A0_A1"]["normalized_source"]=="(g4/x,g5,g6)",
        "Row6 must be the physical raw source on a disk containing the whole core")
    G6=endpoint(c,interval(c,source["source_bounds"]["g6"]))
    F=endpoint(c,interval(c,source["n0_fields"]["F"]))
    tail=parent["actual_selected_raw_word_tail"]
    ar=tail["arithmetic"]
    require(tail["passed"] is True and tail["tail_is_a_directed_majorant_not_a_field_value"] is True and
        ar["N"]==1 and ar["positive_Mbase_and_beta_below_one"] is True and
        interval(c,tail["raw_tail_enclosures"]["0"]["cauchy_gap"])==interval(c,tubes["outer_Z_radius"])/16,
        "Saved uniform raw preconditioned majorants required")
    M=endpoint(c,interval(c,ar["Mbase_upper"]))
    beta1=endpoint(c,interval(c,ar["beta_N_upper"]))
    beta0=2*beta1
    require(G6>0 and F>0 and M>0 and beta0>=0 and beta0<1,
        "Whole-feedback beta0 comparison must be admitted")
    feedback=M*beta0/(1-beta0);full=M/(1-beta0)
    return G6,F,feedback,full,dict(row6_physical_source_upper=record(G6),
        physical_leading_F_complex_upper=record(F),Mbase_upper=record(M),saved_beta1_upper=record(beta1),
        beta0_comparison_upper=record(beta0),all_feedback_upper=record(feedback),full_raw_norm_upper=record(full),
        comparison="for k>=1: Cp/k+C1/loss <= Cp+C1/loss <= 2*(Cp/2+C1/loss)",
        all_feedback_formula="Mbase*beta0/(1-beta0), beta0=2*saved_beta1",
        raw_Uz1_initial_Green_bound="|Wbase_Uz1(rho)|<=rho*G6/4",
        g6_pressure_feedback_already_in_source_cap=True,
        exact_core_source_and_solution_Z_tubes_match=True,whole_core_radial_source_disk_checked=True,
        physical_row6_source_normalization_checked=True,
        caps_are_not_signed_values=True,accepted_norm_producer_replayed=False)


def normalized_enclosure(c,signed,errors,bases):
    dominant=(0,0,1,-2,0)
    require(dominant in signed,"Physical core pressure source sector missing")
    D=signed[dominant];threshold=-1024*c.log(2);tiny=c.mpf(1)/(1<<1024)
    radius=c.mpf(0);proof=[]
    terms=[("other_signed",k,v) for k,v in signed.items() if k!=dominant]
    terms += [("error_radius",k,v) for k,v in errors.items()]
    for kind,key,value in terms:
        upper=endpoint(c,abs(value))
        if upper==0:continue
        logratio=c.log(upper)
        diff=tuple(a-b for a,b in zip(key,dominant))
        for power,base in zip(diff,bases):
            if power:logratio+=power*base
        require(logratio.b<=threshold.a,"A core/C_inner normalized source term exceeds the dyadic cover")
        radius+=tiny
        proof.append(dict(kind=kind,source_exponents=list(key),factor_difference=list(diff),
            log_absolute_normalized_upper=record(logratio),dyadic_upper="2^-1024"))
    result=c.make_mpf(((D-radius)._mpi_[0],(D+radius)._mpi_[1]))
    sign="negative" if result.b<0 else "positive" if result.a>0 else "undetermined"
    require(sign!="undetermined","Actual core/C_inner sign is unresolved")
    width=(result.b-result.a)/min(abs(result.a),abs(result.b))
    return dict(normalizing_positive_factor="Pstar^2*epsilon^2",normalized_integral_enclosure=record(result),
        dominant_source_exponents=list(dominant),dominant_coefficient=record(D),strict_sign=sign,
        relative_interval_width_upper=record(width),all_other_terms_absolute_upper=record(radius),
        dyadic_term_bounds=proof,upper_bounds_not_used_as_signed_values=True)


def exact_algebra():
    r,z,de,j,e=sp.symbols("rho Z delta j epsilon",real=True)
    b=sp.Function("b")(r,z);S=sp.Function("S")(r,z)
    Np=-S/2+2*sp.diff(b,r)+r*sp.diff(b,r,2)
    defect=-e**2*r*Np-sp.diff(-e**2*r*r*sp.diff(b,r),r)-e**2*r*S/2
    checks={"core_Np_primitive_Z"+str(k):sp.expand(sp.diff(defect,z,k))==0 for k in (0,1)}
    u=4*z+j;L=1-de*z*z;d=1-z*z
    A=((-1-de)*z*u+d*4)/L
    nz=(2*z*A-d*sp.diff(A,z))/L
    N=4-(8+4*de)*z*z-(1+de)*j*z
    A0=N/L;A1=(sp.diff(N,z)-sp.diff(L,z)*A0)/L
    A2=(sp.diff(N,z,2)-sp.diff(L,z,2)*A0-2*sp.diff(L,z)*A1)/L
    value=(2*z*A0-d*A1)/L
    derivative=(2*A0+4*z*A1-d*A2-sp.diff(L,z)*value)/L
    checks["axis_N1z_value_formula"]=sp.cancel(value-nz)==0
    checks["axis_N1z_ordinary_Z1_formula"]=sp.cancel(derivative-sp.diff(nz,z))==0
    qa=((1+de)*z*u-4*d)/L
    checks["axis_N1z_matches_actual_Qaxis_transport"]=sp.cancel(nz-(-2*z*qa+d*sp.diff(qa,z))/L)==0
    checks["axis_velocity_integral_model"]=sp.cancel(sp.integrate(2*e*u*(e*nz*r/2),(r,0,4))-8*e**2*u*nz)==0
    checks["axis_inertial_weight"]=sp.integrate(r/2,(r,0,4))==4
    t=sp.symbols("t",positive=True)
    checks["core_axial_Green_constant_source_projection"]=sp.integrate(4*t*sp.log(1/t),(t,0,1))==1
    C0,C1,loss=sp.symbols("Cp C1 loss",positive=True)
    checks["beta0_from_beta1_comparison"]=sp.simplify(2*(C0/2+C1/loss)-(C0+C1/loss))==C1/loss
    require(all(checks.values()),"Core m15 exact source algebra failed")
    return dict(passed=True,checks=checks)


def build(root,parent=None):
    root=Path(root);raw_parent=(root/PARENT).read_bytes()
    if hashlib.sha256(raw_parent).hexdigest()!=PARENT_SHA256:
        raise ValueError("Accepted complete-collar parent changed")
    saved=json.loads(raw_parent)
    # Verify executable arithmetic helpers using stdlib hashing and the
    # hard-pinned parent before trusting their numerical/gating functions.
    for path in HELPERS:
        if hashlib.sha256((root/path).read_bytes()).hexdigest()!=saved["implementation_hashes"][path]:
            raise ValueError("Changed executable core/C_inner helper: "+path)
    if parent is None:parent=saved
    require(digest(parent)==digest(saved) and parent["checks"]["passed"] is True,"Same passed parent required")
    paths=(PARENT,MEAN,MEAN_CHECK,MEAN_SOURCE,DOMAIN,MAJOR,PHYSICAL_SOURCE,RAW_SOURCE,COUPLED_SOURCE)+HELPERS
    pins={}
    for path in paths:
        expected=(PARENT_SHA256 if path==PARENT else parent["input_hashes"][path[len(BASE):]]
            if path.startswith(BASE) else parent["implementation_hashes"][path])
        require(sha(root/path)==expected,"Changed core/C_inner source: "+path)
        pins[path]=dict(path=path,expected_sha256=expected,passed=True)
    mean,check,domain,major=[load(root/path) for path in (MEAN,MEAN_CHECK,DOMAIN,MAJOR)]
    require(check["all_passed"] is True and
        check["input_hashes"][Path(MEAN).name]==pins[MEAN]["expected_sha256"] and
        check["checks"]["original_axis_and_identical_core_actual_first_inlet_factors"] is True and
        mean["source_family"]==domain["source_family"]==parent["source_family"] and
        major["implicit_source_sha256"]==parent["source_family"]["implicit_source_sha256"] and
        major["datum_enclosure_sha256"]==parent["source_family"]["datum_enclosure_sha256"],
        "Checked actual axis/core/own mean and pressure family required")
    c=MPIntervalContext();c.prec=PRECISION_BITS
    A=mean["samples"]["actual_first_inlet"];O=mean["samples"]["core_axis"]
    require(A["Z"]==O["Z"]==parent["actual_complete_first_collar_m15_increment"]["selected_Z_interval"],
        "Core and collar must use the identical saved axial interval")
    require(O["geometry"]["chart"]=="core" and interval(c,O["geometry"]["rho"])==0,
        "Canonical analytic core axis required for the zero inlet premise")
    for k in (0,1):
        for key in ("leading_Psi_rho_rows","leading_own_meanPsi_rho_rows","leading_Qpsi_rho_rows"):
            axis=lift(c,O[key][0],k)
            require(all(value==0 for value in axis.values()),"Actual regular core Psi/mean/Qpsi must vanish at the axis")
    delta=interval(c,parent["actual_selected_leading_pressure_primitive"]["delta_enclosure"])
    j=interval(c,major["required_j"]);z=interval(c,A["Z"])
    require(delta>0 and delta<1 and z.a>-1 and z.b<1,"Actual core parameters required")
    norms,S1,S2,Nza,core_proof=core_bounds(c,parent,mean,domain,delta,j)
    G6,BF,feedback,Mfull,raw_proof=raw_core_bounds(c,parent,domain)
    nz,u=axis_Nz(c,z,j,delta);Saxis=axis_jet(c,A,delta)
    B1=lift(c,mean["samples"]["selected_actual_inner_exit"]["uniform_complex_factored_bounds"]["Psi"][1],upper=True)
    Ubound=add(scalar(c,norms["uv"]),shift(scale(B1,4),EPS))
    raw_error=add(scalar(c,G6+feedback),shift(scalar(c,2*Nza),EPS))
    leading_difference=shift(scale(B1,4),EPS)
    axial_error=shift(scale(add(product(Ubound,raw_error),
        product(leading_difference,shift(scalar(c,2*Nza),EPS))),8),EPS)
    swirl_error=shift(scalar(c,16*BF*Mfull),(0,0,0,-2,0))
    pressure=parent["actual_correlated_leading_pressure_increment"]
    bases=[interval(c,v) for v in pressure["first_four_source_log_bases"]]+[
        interval(c,pressure["flat_source_atom"]["log_sigma_bounds"])]
    rows={}
    for k in (0,1):
        visc=shift(scale(lift(c,A["leading_Qpsi_rho_rows"][1],k),-16),(0,0,0,-2,0))
        inertia=shift(scalar(c,4*Saxis[k]),(0,0,0,-2,0))
        value=u[0]*nz[0] if k==0 else u[1]*nz[0]+u[0]*nz[1]
        axial=shift(scalar(c,8*value),(0,0,0,-2,0))
        factor=1 if k==0 else 1/norms["gap"]
        inertial_error=add(shift(scale(S1,4*factor),(0,0,0,-3,0)),
            shift(scale(S2,4*factor),(0,0,0,-4,0)))
        errors=add(inertial_error,scale(axial_error,factor),scale(swirl_error,factor))
        signed=add(visc,inertia,axial)
        collar=parent["actual_complete_first_collar_m15_increment"]["ordinary_Z"][str(k)]["complete_first_collar_m15"]
        full_signed=add(signed,saved_polynomial(c,collar["signed_approximation"]))
        full_error=add(errors,saved_polynomial(c,collar["error_radius"]))
        rows[str(k)]=dict(core_viscous_endpoint=export(visc),core_axis_inertia=export(inertia),
            core_axial_axis_source_approximation=export(axial),core_inertial_error=export(inertial_error,error=True),
            core_axial_error=export(scale(axial_error,factor),error=True),
            core_swirl_coupling_unsigned_error=export(scale(swirl_error,factor),error=True),
            actual_core_m15=dict(signed_approximation=export(signed),error_radius=export(errors,error=True),
                normalized=normalized_enclosure(c,signed,errors,bases)),
            actual_C_inner_m15=dict(signed_approximation=export(full_signed),error_radius=export(full_error,error=True),
                normalized=normalized_enclosure(c,full_signed,full_error,bases)))
    guards=dict(
        core_physical_axis=ast_guard(root,PHYSICAL_SOURCE,"axis_inputs",("u=4*z+self.j",)),
        regular_core_axis=ast_guard(root,PHYSICAL_SOURCE,"normalized_jets",(
            "if axis and i==0:",
            "phi[index]=c.mpf(1 if k==0 else 0);psi[index]=average[index]=c.mpf(0)")),
        raw_pressure=ast_guard(root,RAW_SOURCE,"raw_system",(
            "Np=a.add(a.scale('-1/2',S),a.scale(parameters['Lambda'],a.add(a.scale(2,c1(Q[1])),a.scale(r,c1(Q[2])))))",
            "Nz=g.neg(axial_second(g,U,r,de,g.sub(c(-1),de)))")),
        raw_core_Green=ast_guard(root,RAW_SOURCE,"green_apply",(
            "sample=g.mul(r,t,t);measure=X",
            "log=g.mul(X,t,g.unary('log',g.quotient(g.one,t,'positive t inside unit integral')))")),
        raw_feedback_tail=ast_guard(root,COUPLED_SOURCE,"preconditioned_tail",(
            "beta=domain.upper(c,c.e*a*S*(Cp/(N+1)+C1/loss))",
            "arbitrary_k_term_formula='Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!'")))
    return dict(passed=True,source_family=parent["source_family"],parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,parent_mapping_sha256=digest(parent),source_pins=pins,
        source_AST=guards,exact_algebra=exact_algebra(),ordinary_Z_orders=[0,1],ordinary_Z=rows,
        scope="One saved Z interval near0.371; ordinary derivatives at fixed native rho and endpoints",
        one_point_only=True,
        selected_Z_interval=A["Z"],precision_bits=PRECISION_BITS,core_bounds_proof=core_proof,
        raw_core_bounds_proof=raw_proof,source_log_basis_names=pressure["source_log_basis_names"],
        core_identity="Mcore=-16*eps^2*b_rho(4)+eps^2/2*int_0^4 rho*S+2eps*intUlead*Uz1-2eps^2*int rho*Flead*F1",
        signed_axial_model="8*epsilon^2*Uaxis*N1z_axis; N1z_axis=-Z_-2 Z_(-1-delta) Uaxis",
        axis_endpoint_rho_squared_b_rho_limit_zero=True,core_mean_not_reset=True,
        actual_axis_Psi_mean_Qpsi_zero_receipt_and_source_checked=True,executable_helper_pins_checked=True,
        actual_core_m15_value_enclosed=True,actual_C_inner_m15_value_enclosed=True,
        complete_first_collar_added_once=True,core_source_handles=[53125,53126],collar_source_handles=[53137,53138],
        C_inner_source_handles=[53161,53162],swirl_core_integral_only_bounded_not_signed_value=True,
        old_error_caps_or_six_ledger_terms_removed=False,all_Z_functional_coverage=False,
        independent_terminal_comparison_closed=False,global_five_moment_closure=False,
        source_executed=False,accepted_ancestors_replayed=False,P1_is_distinct_from_m15=True)
