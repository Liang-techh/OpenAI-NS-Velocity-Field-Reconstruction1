"""Signed first-collar cut-velocity m15 integrals with separate source errors.

The approximation is an actual saved inlet product, integrated against the
exact original eta mass. It is never a midpoint of a positive norm bound.
All field variation, omitted words and radial geometry remain error radii.
"""
from __future__ import annotations

from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
import sympy as sp

from .original_finite_N_signed_source import ast_guard,digest,load,require,sha
from .original_n1_selected_word_numeric import interval,endpoint,record
from .original_n1_selected_pressure_primitive import add,scale,shift,product,export
from .original_n1_correlated_pressure_increment import lift,hydrate,scalar,normalized_enclosure

BASE="experiments/root_st073/"
PREFIX=BASE+"lei_ren_part1_paper_compliant_current_original_n1_"
PARENT=PREFIX+"correlated_pressure_increment.json"
PARENT_SHA256="895e772d41b0f0756648aa66c08c788e2cd93c096b50dfb7cd62d56a9499e3ee"
ZERO=(0,0,0,0,0)
H=(1,0,0,0,0)
PRECISION_BITS=1800
MEAN=PREFIX+"leading_mean_epsilon_factor.json.gz"
MEAN_CHECK=PREFIX+"leading_mean_epsilon_factor_check.json"
MEAN_SOURCE=PREFIX+"leading_mean_epsilon_factor.py"
OPERATOR=PREFIX+"inner_operator.json.gz"
OPERATOR_CHECK=PREFIX+"inner_operator_check.json"
OPERATOR_SOURCE=PREFIX+"inner_operator.py"
DOMAIN=PREFIX+"inner_analytic_domain.json.gz"
DOMAIN_CHECK=PREFIX+"inner_analytic_domain_check.json"
DOMAIN_SOURCE=PREFIX+"inner_analytic_domain.py"


def complex_F0_cap(c,root,parent,pins,domain):
    check=load(root/DOMAIN_CHECK)
    require(check["all_passed"] is True and
        check["input_hashes"][Path(DOMAIN).name]==pins[DOMAIN]["expected_sha256"],
        "Accepted complex-domain receipt must bind the saved archive")
    source=parent["n1_complex_actual_source_bounds"]
    embedding=source["xh_embedding"]
    require(source["passed"] is True and source["actual_source_bounds_proved"] is True and
        all(embedding[k] is True for k in ("proved","F0_complex_axis_binding",
            "nested_complex_Z_tubes","outer_tube_within_Xh_Z_disk")),
        "Actual physical F0 bound must be admitted on the common complex Z tube")
    amplitude=domain["original_axis_amplitude"]
    require(amplitude["source_formula"]=="exp(-selected_logCstar-Lambda*G(Z))",
        "The uniform F0 cap cannot be replaced by a center F0_anchor amplitude")
    cap=endpoint(c,interval(c,amplitude["complex_modulus_upper"]))
    require(cap>0,"Positive physical complex F0 upper bound required")
    return cap,dict(archive=DOMAIN,checked_receipt=DOMAIN_CHECK,
        source_formula=amplitude["source_formula"],complex_modulus_upper=record(cap),
        parent_actual_complex_source_bounds_passed=True,common_complex_Z_tubes_admitted=True,
        physical_complex_cap_not_center_anchor=True,cap_used_only_as_error_bound=True)


def source_jet(c,row,order):
    # This row already represents the physical ordinary derivative. The
    # F0 connection is in its coefficients, not a missing derivative of
    # the fixed anchor log. Never add that connection a second time.
    result={}
    for term in row["terms"]:
        key=tuple(term["source_exponents"])
        require(len(key)==4,"Four original fixed source factors required")
        result=add(result,{key+(0,):interval(c,term["ordinary_Z_derivatives"][order])})
    return result


def leading_inlets(c,parent,operator,mean,domain,F0_complex_cap):
    A=operator["samples"]["first_inlet"]
    W=parent["actual_selected_first_word_enclosure"]
    MA=mean["samples"]["actual_first_inlet"]
    MB=mean["samples"]["selected_actual_inner_exit"]
    require(A["Z"]==MA["Z"]==W["selected_Z_interval"] and
        A["amplitude_anchor_fixed_during_jet_algebra"] is True and
        A["ordinary_Z_amplitude_derivatives_in_coefficients"] is True and
        A["outputs_are_enclosures_of_original_functions_not_selected_point_values"] is True,
        "Same actual leading inlet with complete physical Z derivatives required")
    require(A["source_log_basis_names"]==W["source_log_basis_names"] and
        A["source_log_bases"]==W["source_log_bases"],
        "Identical fixed F0/source anchors required before multiplying leading and raw jets")
    require(interval(c,A["geometry"]["phase"])==0 and
        interval(c,A["geometry"]["rho_enclosure_only"])==4,
        "Actual core-to-first inlet must be rho4, phase0")
    fields=A["ordinary_leading_rho_rows"]
    F=[source_jet(c,fields["F"][0],k) for k in (0,1)]
    U=[lift(c,MA["leading_factored_Uz_rho_rows"][0],k) for k in (0,1)]
    old=parent["actual_selected_leading_pressure_primitive"]
    Cd=interval(c,old["radius_factors"]["rhoB_minus4_over_hb"])
    Cw=interval(c,old["radius_factors"]["integral_rho_drho_over_hb"])
    proof=MB["uniform_complex_bound_proof"]
    eta=endpoint(c,interval(c,proof["original_outer_eta"]),0)
    gap=eta/16
    require(interval(c,domain["common_tubes"]["outer_Z_radius"])==eta,
        "Leading F/Psi source bounds must share the same admitted outer Z tube")
    Ac=endpoint(c,interval(c,proof["actual_core_logarithmic_derivative_modulus_upper"]))
    floor=endpoint(c,interval(c,domain["denominator_bounds"]["Phi_core_modulus_lower"]),0)
    phirow=domain["actual_collar_radial_bounds"]["Phi"][0]["0"]
    require(phirow["exact_width_power"]==0 and interval(c,phirow["log_width_scale"])==0,
        "Actual leading Phi value bound cannot omit an inverse-width factor")
    phi=endpoint(c,interval(c,phirow["coefficient"]))
    require(floor>0 and phi>0 and F0_complex_cap>0 and gap>0,
        "Common Z tube needs an actual nonzero core Phi and positive physical F cover")
    require(interval(c,proof["Psi_and_meanPsi_bound_tube_radius"])==eta/2,
        "Actual native Psi_rho complex Z tube required")
    actual_Psi_r=lift(c,MB["uniform_complex_factored_bounds"]["Psi"][1],upper=True)
    EF=shift(scalar(c,Cd*Ac*phi*F0_complex_cap),H)
    EU=shift(scale(actual_Psi_r,Cd),(1,0,0,-1,0))
    errors=dict(F=[EF,scale(EF,1/gap)],U=[EU,scale(EU,1/gap)])
    return dict(F=F,U=U),errors,Cd,Cw,dict(
        actual_signed_leading_inlet_hydrated=True,archive=OPERATOR,sample="first_inlet",
        physical_F_definition="F_lead=F0*Phi_actual",
        F0_Z_connection_in_saved_leading_and_raw_coefficients=True,
        F0_connection_added_again=False,identical_source_log_anchors=True,
        common_complex_Z_source_radius=record(eta/2),whole_difference_Cauchy_gap=record(gap),
        Phi_core_modulus_lower=record(floor),actual_Phi_modulus_upper=record(phi),
        F0_complex_cap=record(F0_complex_cap),core_log_derivative_cover=record(Ac),
        leading_F_variation="|DeltaF|<=hb*Cd*Ac*Phi_bound*F0_complex_cap",
        leading_U_variation="|DeltaU|<=epsilon*hb*Cd*actual_Psi_rho_bound",
        Z1="Cauchy only in Z on whole difference, gap eta/16; rho/s integration stays real",
        real_chi_bound="0<=sigma<=1 and 0<hb<1/2 imply hb<=chi<=1",
        complex_rho_or_sigma_endpoint_extension_used=False)


def raw_inlet_and_variation(c,parent):
    words=parent["actual_selected_first_word_enclosure"]
    tail=parent["actual_selected_raw_word_tail"]
    require(words["passed"] is True and words["ordinary_Z_orders"]==[0,1] and
        words["raw_allword_tail_kept_separate"] is True and
        tail["passed"] is True and tail["tail_is_a_directed_majorant_not_a_field_value"] is True,
        "Actual finite prefix and separate uniform raw all-word tail required")
    rows=words["selected_coordinate_rows"]
    A=next(r for r in rows if r["q"]=="0")
    B=next(r for r in rows if r["q"]=="1")
    s=interval(c,B["s_exact_selected_enclosure"])
    ymax=endpoint(c,interval(c,B["y_upper"]))
    require(s.a>0 and ymax==s/2 and B["s_definition"]=="sc*q/4",
        "Selected homogeneous displacement must be sc/8, not source-kernel dx")
    inlets,errors={},{}
    for field,dx in (("F1","d_x_F1"),("Uz1","d_x_Uz1")):
        inlets[field]=[];errors[field]=[]
        for k in ("0","1"):
            ar,br=A["ordinary_Z"][k],B["ordinary_Z"][k]
            inlet=hydrate(c,ar["prefix"][field])
            derivative=hydrate(c,ar["prefix"][dx])
            local=endpoint(c,interval(c,br["local_raw_error_radii"]["prefix"][field]))
            omitted=endpoint(c,interval(c,tail["raw_tail_enclosures"][k]["tail_upper"]))
            require(local>=0 and omitted>=0 and ar["allword_tail_included"] is False,
                "Local and omitted-word radii must be separately nonnegative")
            inlets[field].append(inlet)
            errors[field].append(add(scale(absolute(c,derivative),ymax),scalar(c,local+omitted)))
    return inlets,errors,dict(selected_phase=record(s),homogeneous_y_upper=record(ymax),
        uniformity="q1 local_radii monotone in q; y in [0,sc/8] over the whole selected window",
        raw_variation="ymax*abs(prefix_d_x_field_A)+local_prefix_radius(q1)+raw_allword_tail",
        one_uniform_tail_added_at_evaluation_point=True,
        exact_raw_inlet_not_replaced_by_finite_prefix=True,
        source_interval_and_core_quadrature_errors_retained_in_prefix=True)


def absolute(c,p):
    return {k:endpoint(c,abs(v)) for k,v in p.items()}


def pjet_product(a,b):
    return [product(a[0],b[0]),add(product(a[1],b[0]),product(a[0],b[1]))]


def pjet_add(*jets):
    return [add(*(j[k] for j in jets)) for k in (0,1)]


def product_error(c,a,b,ea,eb):
    # a,b are actual signed source approximations; ea,eb are symmetric radii
    # to the true fields, for values and ordinary Z derivatives separately.
    aa,bb=[[absolute(c,p) for p in j] for j in (a,b)]
    return pjet_add(pjet_product(aa,eb),pjet_product(bb,ea),pjet_product(ea,eb))


def integrate_product(c,a,b,ea,eb,s,Cd,Cw,radial_power,eps_power,sign):
    require(radial_power in (1,2) and eps_power in (1,2) and sign in (-1,1),
        "Only the two physical cut-velocity m15 terms are admitted")
    require(s.a>0 and Cd.a>0 and Cw.a>0,
        "Positive same-source phase and directed geometry factors required")
    mass=5*s/8
    model=pjet_product(a,b)
    err=product_error(c,a,b,ea,eb)
    aa,bb=[[absolute(c,p) for p in j] for j in (a,b)]
    uniform=pjet_add(pjet_product(aa,bb),err)
    # rhoB-4=hb*Cd, rhoB^2-16=2*hb*Cw. Geometry is a
    # source-factor error, not independently boxed endpoint subtraction.
    delta=Cd if radial_power==1 else 2*Cw
    weight=4**radial_power
    exponents=(1,0,0,-eps_power,0)
    signed=[shift(scale(p,sign*2*weight*mass),exponents) for p in model]
    radii=[shift(scale(add(scale(e,weight),shift(scale(u,delta),H)),2*mass),exponents)
        for e,u in zip(err,uniform)]
    return signed,radii


def exact_algebra():
    x,a,b,ea,eb=sp.symbols("x a b ea eb")
    checks={"bilinear_error_identity":sp.expand((a+ea)*(b+eb)-a*b-a*eb-b*ea-ea*eb)==0}
    S=sp.symbols("S",positive=True)
    checks["axial_inlet_eta_mass_factor"]=2*4*sp.Rational(5,8)*S==5*S
    checks["swirl_inlet_eta_mass_factor"]=-2*16*sp.Rational(5,8)*S==-20*S
    az,bz,eaz,ebz=sp.symbols("az bz eaz ebz")
    first=az*eb+a*ebz+bz*ea+b*eaz+eaz*eb+ea*ebz
    direct=sp.diff((a+x*az+ea+x*eaz)*(b+x*bz+eb+x*ebz)-(a+x*az)*(b+x*bz),x).subs(x,0)
    checks["ordinary_Z_product_error_identity"]=sp.expand(direct-first)==0
    h=sp.symbols("hb",positive=True)
    rho=4*sp.exp(h*S)
    cd=4*(sp.exp(h*S)-1)/h
    cw=8*(sp.exp(2*h*S)-1)/h
    checks["rho_geometry_identity"]=sp.simplify(rho-4-h*cd)==0
    checks["rho_squared_geometry_identity"]=sp.simplify(rho*rho-16-2*h*cw)==0
    t=sp.symbols("t",real=True)
    checks["positive_expm1_by_real_FTC"]=(sp.diff(sp.exp(t)-1,t)==sp.exp(t) and
        (sp.exp(t)-1).subs(t,0)==0 and sp.exp(t).is_positive is True)
    require(all(checks.values()),"Cut-velocity integral product algebra failed")
    return dict(passed=True,checks=checks)


def saved_polynomial(c,rows):
    return {tuple(t["source_exponents"]):interval(c,t["coefficient_interval"]) for t in rows}


def relative_component_bound(c,signed,errors,bases,s):
    # Bound each new component separately before it is merged with pressure.
    # This is an absolute upper bound, never the component's signed value.
    dominant=(1,0,1,-2,0);threshold=-1024*c.log(2)
    count=0;maximum=None
    for p in (signed,errors):
        for key,value in p.items():
            bound=endpoint(c,abs(value))
            if bound==0:continue
            logratio=c.log(bound)-c.log(s)
            for power,base in zip((a-b for a,b in zip(key,dominant)),bases):
                if power:logratio+=power*base
            require(logratio.b<=threshold.a,"New component error exceeds admitted relative dyadic bound")
            top=endpoint(c,logratio)
            maximum=top if maximum is None or top.b>maximum.b else maximum
            count+=1
    require(count>0,"A physical component must retain actual source/error terms")
    return dict(normalizing_positive_factor="hb*Pstar^2*epsilon^2*s",
        separate_signed_and_error_sparse_terms=count,each_term_absolute_upper="2^-1024",
        normalized_absolute_upper=record(c.mpf(count)/(1<<1024)),
        maximum_log_absolute_normalized_term_upper=record(maximum),
        absolute_upper_is_not_a_signed_component_value=True)


def build(root,parent=None):
    root=Path(root);require(sha(root/PARENT)==PARENT_SHA256,"Accepted correlated pressure parent changed")
    saved=load(root/PARENT)
    if parent is None:parent=saved
    require(digest(parent)==digest(saved) and parent["checks"]["passed"] is True,
        "Same passed correlated pressure parent required")
    paths=(PARENT,MEAN,MEAN_CHECK,MEAN_SOURCE,OPERATOR,OPERATOR_CHECK,OPERATOR_SOURCE,
        DOMAIN,DOMAIN_CHECK,DOMAIN_SOURCE)
    pins={}
    for path in paths:
        expected=PARENT_SHA256 if path==PARENT else parent["input_hashes"][path[len(BASE):]]
        require(sha(root/path)==expected,"Changed velocity/coupling source: "+path)
        pins[path]=dict(path=path,expected_sha256=expected,passed=True)
    mean,mc,op,oc,domain=[load(root/path) for path in (MEAN,MEAN_CHECK,OPERATOR,OPERATOR_CHECK,DOMAIN)]
    require(mc["all_passed"] is True and oc["all_passed"] is True and
        mc["input_hashes"][Path(MEAN).name]==pins[MEAN]["expected_sha256"] and
        oc["input_hashes"][Path(OPERATOR).name]==pins[OPERATOR]["expected_sha256"] and
        mean["source_family"]==op["source_family"]==domain["source_family"]==parent["source_family"],
        "Checked actual leading/raw sources of one family required")
    require(parent["actual_cutoff_complete_support_measure"]["passed"] is True and
        parent["actual_cutoff_complete_support_measure"]["exact_unweighted_q_mass"]=="1/2+(1/4)*(1/2)=5/8" and
        parent["actual_cutoff_complete_support_measure"]["measure_is_Z_independent"] is True and
        parent["actual_n1_cutoff_history_identity"]["passed"] is True,
        "Original own velocity cutoff and exact source mass required")
    c=MPIntervalContext();c.prec=PRECISION_BITS
    # Filled by the source-bound complex-amplitude gate below; never infer a
    # physical F bound from its center value or a complex-rho continuation.
    Fcap,Fproof=complex_F0_cap(c,root,parent,pins,domain)
    lead,le,cd,cw,lp=leading_inlets(c,parent,op,mean,domain,Fcap)
    raw,re,rp=raw_inlet_and_variation(c,parent)
    s=interval(c,rp["selected_phase"])
    require(s==interval(c,parent["actual_correlated_leading_pressure_increment"]["selected_phase"]),
        "Same exact full first-collar endpoint required")
    Us,Ue=integrate_product(c,raw["Uz1"],lead["U"],re["Uz1"],le["U"],s,cd,cw,1,1,1)
    Fs,Fe=integrate_product(c,raw["F1"],lead["F"],re["F1"],le["F"],s,cd,cw,2,2,-1)
    pressure=parent["actual_correlated_leading_pressure_increment"]
    bases=[interval(c,x) for x in pressure["first_four_source_log_bases"]]+[
        interval(c,pressure["flat_source_atom"]["log_sigma_bounds"])]
    rows={}
    for k in (0,1):
        P=pressure["ordinary_Z"][str(k)]["m15_leading_N1p_local"]
        ps=saved_polynomial(c,P["signed_approximation"])
        pe=saved_polynomial(c,P["correlated_error_radius"])
        signed=add(ps,Us[k],Fs[k]);errors=add(pe,Ue[k],Fe[k])
        rows[str(k)]=dict(
            cut_axial_velocity=dict(signed_inlet_product_integral=export(Us[k]),error_radius=export(Ue[k],error=True),
                relative_absolute_bound=relative_component_bound(c,Us[k],Ue[k],bases,s)),
            cut_swirl_coupling=dict(signed_inlet_product_integral=export(Fs[k]),error_radius=export(Fe[k],error=True),
                relative_absolute_bound=relative_component_bound(c,Fs[k],Fe[k],bases,s)),
            complete_first_collar_m15=dict(signed_approximation=export(signed),
                error_radius=export(errors,error=True),normalized=normalized_enclosure(c,signed,errors,bases,s,2)))
    guards=dict(
        complex_physical_F_cover=ast_guard(root,DOMAIN_SOURCE,"hierarchy_bounds",(
            "f=[v*Fmax for v in phi]")),
        original_complex_F0_cap=ast_guard(root,DOMAIN_SOURCE,"domain",(
            "Fmax=upper(c,c.exp(-2*core.logLambda-1000))",
            "original_axis_amplitude=dict(complex_modulus_upper=Fmax")),
        F0_connection=ast_guard(root,OPERATOR_SOURCE,"first_rows",(
            "ratioF=Jet(c,[inputs['F0_ratios'][k]/math.factorial(k) for k in range(6)])",
            "frows=[algebra.shift(value*ratioF,(0,1,0,0)) for value in phirows]")),
        physical_derivative_export=ast_guard(root,OPERATOR_SOURCE,"export",(
            "ordinary_Z_derivatives=[math.factorial(k)*v for k,v in enumerate(jet.coefficients)]")),
        actual_Psi_complex_bound=ast_guard(root,MEAN_SOURCE,"complex_mean_bounds",(
            "psi1=core_psi[1]*(ratio*(1+cap))",
            "Psi_and_meanPsi_bound_tube_radius=eta/2")))
    return dict(passed=True,source_family=parent["source_family"],parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,parent_mapping_sha256=digest(parent),source_pins=pins,
        source_AST=guards,exact_algebra=exact_algebra(),raw_variation_proof=rp,leading_variation_proof=lp,
        complex_F0_proof=Fproof,ordinary_Z_orders=[0,1],ordinary_Z=rows,
        selected_phase=rp["selected_phase"],selected_Z_interval=pressure["selected_Z_interval"],
        source_log_basis_names=pressure["source_log_basis_names"],precision_bits=PRECISION_BITS,
        physical_density="2*epsilon*Ulead*eta*Uz1raw-epsilon^2*rho*(2*Flead*eta*F1raw+N1p)",
        actual_complete_first_collar_m15_increment_enclosed=True,
        actual_saved_signed_inlet_products_integrated_with_separate_errors=True,
        original_nonzero_leading_sigma_and_own_eta_distinct=True,
        raw_allword_tail_added_once_per_field_at_evaluation_point=True,
        core_contribution_added=False,full_C_inner_value_evaluated=False,
        all_Z_functional_coverage=False,old_error_caps_or_six_ledger_terms_removed=False,
        P1_is_distinct_from_m15=True,source_executed=False,accepted_ancestors_replayed=False,
        global_five_moment_closure=False)
