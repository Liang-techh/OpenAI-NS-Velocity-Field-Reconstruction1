"""Signed selected-collar leading-pressure primitives with factored remainders.

Only saved leading endpoint/complex bounds are hydrated. The physical source
is N1p=-S/2+Lambda*(2*Q_rho+rho*Q_rhorho), Q=Qaxis+epsilon*Qpsi.
The viscous part is evaluated by exact FTC before endpoint interval boxing.
"""
from __future__ import annotations

from pathlib import Path
from mpmath.ctx_iv import MPIntervalContext
import sympy as sp

from .original_finite_N_signed_source import ast_guard, digest, load, require, sha
from .original_n1_selected_word_numeric import interval, endpoint, record

BASE = "experiments/root_st073/"
PREFIX = BASE + "lei_ren_part1_paper_compliant_current_original_n1_"
PARENT = PREFIX + "selected_words.json"
PARENT_SHA256 = "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e"
MEAN = PREFIX + "leading_mean_epsilon_factor.json.gz"
MEAN_CHECK = PREFIX + "leading_mean_epsilon_factor_check.json"
MEAN_SOURCE = PREFIX + "leading_mean_epsilon_factor.py"
DOMAIN = PREFIX + "inner_analytic_domain.json.gz"
GRAPH = PREFIX + "lower_collar_functions.json.gz"
RAW_SOURCE = PREFIX + "inner_volterra_functions.py"
OP_SOURCE = PREFIX + "inner_operator.py"
CUT_SOURCE = PREFIX + "inner_cutoff_functions.py"
ZERO = (0, 0, 0, 0)
PRECISION_BITS = 1800


def add(*polys):
    out = {}
    for p in polys:
        for k, v in p.items():
            out[k] = out[k] + v if k in out else v
    return out


def scale(p, a):
    return {k: v*a for k, v in p.items()}


def shift(p, delta):
    return {tuple(a+b for a,b in zip(k,delta)):v for k,v in p.items()}


def product(p, q):
    out = {}
    for a, av in p.items():
        for b, bv in q.items():
            k = tuple(x+y for x,y in zip(a,b))
            out[k] = out[k]+av*bv if k in out else av*bv
    return out


def poly(c, row, order=0, positive_upper=False):
    out = {}
    for term in row["terms"]:
        key = tuple(term["source_exponents"])
        require(len(key)==4 and key[1]==0, "Only Z-independent scale sectors are admitted")
        v = interval(c, term["ordinary_Z_derivatives"][order])
        if positive_upper:
            v = endpoint(c,v)
            require(v>=0, "Positive complex norm required")
        out = add(out,{key:v})
    return out


def export(p, error=False):
    result = []
    for k,v in sorted(p.items()):
        if error:
            require(v.b>=0, "Nonnegative error radius required")
            bits = ((0,0,0,0), v._mpi_[1])
            value = v.ctx.make_mpf(bits)
        else:
            value = v
        result.append(dict(source_exponents=list(k),coefficient_interval=record(value)))
    return result


def jadd(*jets):
    return [sum(v[k] for v in jets) for k in (0,1)]


def jmul(a,b):
    return [a[0]*b[0],a[1]*b[0]+a[0]*b[1]]


def jscale(a,s):
    return [s*v for v in a]


def jdiv(a,b):
    return [a[0]/b[0],(a[1]*b[0]-a[0]*b[1])/(b[0]*b[0])]


def axis_jet(c,sample,delta):
    qa=sample["leading_Qaxis"]
    require(len(qa["terms"])==1 and tuple(qa["terms"][0]["source_exponents"])==ZERO,
        "Exact radial-constant Qaxis sector required")
    av=[interval(c,q) for q in qa["terms"][0]["ordinary_Z_derivatives"][:3]]
    up=poly(c,sample["leading_factored_Uz_rho_rows"][0],0)
    uz=poly(c,sample["leading_factored_Uz_rho_rows"][0],1)
    require(ZERO in up and ZERO in uz, "Saved actual Uaxis and derivative required")
    z=[interval(c,sample["Z"]),c.mpf(1)]
    one=[c.mpf(1),c.mpf(0)]
    d=jadd(one,jscale(jmul(z,z),-1))
    L=jadd(one,jscale(jmul(z,z),-delta))
    a,az,u=av[:2],av[1:3],[up[ZERO],uz[ZERO]]
    time=jdiv(jadd(a,jscale(jmul(z,az),(1-delta)/2)),L)
    transport=jdiv(jmul(u,jadd(jscale(jmul(z,a),-2),jmul(d,az))),L)
    return jadd(time,jscale(jmul(a,a),c.mpf(1)/2),transport)


def graph_scalar(c,nodes,index,visited):
    require(len(visited)<64, "Parameter subgraph exceeds bounded scalar scope")
    visited.add(index);n=nodes[index];op=n["operation"]
    if op=="exact_rational":
        return c.mpf(n["numerator"])/n["denominator"]
    if op=="negative":
        return -graph_scalar(c,nodes,n["argument"],visited)
    if op in ("sum","product"):
        values=[graph_scalar(c,nodes,i,visited) for i in n["arguments"]]
        out=c.mpf(0 if op=="sum" else 1)
        for v in values:
            out=out+v if op=="sum" else out*v
        return out
    if op=="analytic_unary" and n["name"]=="exp":
        x=graph_scalar(c,nodes,n["argument"],visited)
        require(abs(x)<c.mpf("1e20"), "Refuse materializing microscopic width exponentials")
        return c.exp(x)
    raise ValueError("Unsupported exact parameter operation: "+str(op))


def inertial_norms(c,sample,domain):
    u=sample["uniform_complex_factored_bounds"]
    proof=sample["uniform_complex_bound_proof"]
    eta=endpoint(c,interval(c,proof["original_outer_eta"]),0)
    require(interval(c,proof["Qpsi_and_Q_bound_tube_radius"])==eta/4,
        "Same-source complex Qpsi tube required")
    require(interval(c,proof["Psi_and_meanPsi_bound_tube_radius"])==eta/2 and
        proof["norm_bounds_not_selected_comparison_field"] is True and
        proof["all_original_pressure_atoms_retained"] is True,
        "Actual whole-source complex Psi/Qpsi bounds required")
    loss,gap=eta/8,eta/16
    z=1+eta;d=1+z*z
    L=endpoint(c,interval(c,domain["denominator_bounds"]["L_modulus_lower"]),0)
    require(L>0,"Positive common complex denominator bound required")
    # The actual delta is Z independent; its graph value is passed separately.
    a=poly(c,u["Q"][0],positive_upper=True)[ZERO]
    uv=poly(c,u["Uz"][0],positive_upper=True)[ZERO]
    b=poly(c,u["Qpsi"][0],positive_upper=True)
    br=poly(c,u["Qpsi"][1],positive_upper=True)
    v=poly(c,u["Psi"][0],positive_upper=True)
    rho=endpoint(c,interval(c,domain["actual_inner_geometry"]["real_rho_upper_bound"]))
    return dict(eta=eta,loss=loss,gap=gap,z=z,d=d,L=L,a=a,uv=uv,b=b,br=br,v=v,rho=rho)


def norm_polynomials(c,n,delta):
    z,d,L,r=n["z"],n["d"],n["L"],n["rho"]
    a,u,b,br,v=n["a"],n["uv"],n["b"],n["br"],n["v"]
    az=a/n["loss"];bz=scale(b,1/n["loss"])
    Baxis=(2*z*a+d*az)/L
    Bb=scale(add(scale(b,2*z),scale(bz,d),scale(br,2*z*r)),1/L)
    time=scale(add(b,scale(bz,z*(1+endpoint(c,delta))/2),scale(br,r)),1/L)
    S1=add(time,scale(b,a),scale(br,r*a),scale(Bb,u),scale(v,Baxis))
    S2=add(scale(product(b,b),c.mpf(1)/2),scale(product(b,br),r),product(v,Bb))
    return S1,S2


def exact_algebra():
    r,z=sp.symbols("rho Z",real=True)
    e=sp.symbols("epsilon",positive=True);de=sp.symbols("delta",real=True)
    a,u=sp.Function("a")(z),sp.Function("u")(z)
    b,v=sp.Function("b")(r,z),sp.Function("v")(r,z)
    L,d=1-de*z*z,1-z*z
    A=lambda q:(q+z*(1-de)*sp.diff(q,z)/2+r*sp.diff(q,r))/L
    B=lambda q:(-2*z*q+d*sp.diff(q,z)-2*z*r*sp.diff(q,r))/L
    S=lambda q,w:A(q)+q*(q/2+r*sp.diff(q,r))+w*B(q)
    S0=S(a,u);S1=A(b)+a*b+r*a*sp.diff(b,r)+u*B(b)+v*B(a)
    S2=b*b/2+r*b*sp.diff(b,r)+v*B(b)
    q=a+e*b;full=S(q,u+e*v)
    Np=-full/2+(2*sp.diff(q,r)+r*sp.diff(q,r,2))/e
    checks=dict(epsilon_polynomial_split=sp.expand(full-S0-e*S1-e**2*S2)==0)
    h,cd,ba,bb,ra,rb=sp.symbols("hb Cd ba bb ra rb")
    rB=4+h*cd;cw=(rB*rB-16)/(2*h)
    checks["P1_geometry_difference_factored_before_boxing"]=sp.expand(
        bb+rB*rb-ba-4*ra-((bb-ba)+4*(rb-ra)+h*cd*rb))==0
    checks["m15_geometry_difference_factored_before_boxing"]=sp.expand(
        rB*rB*rb-16*ra-(16*(rb-ra)+2*h*cw*rb))==0
    for k in (0,1):
        checks["Qaxis_radial_cancellation_Z"+str(k)]=sp.diff(a,r,z,k)==0
        pressure=e*Np-sp.diff(e*(b+r*sp.diff(b,r)),r)+e*full/2
        m15=-e**2*r*Np-sp.diff(-e**2*r*r*sp.diff(b,r),r)-e**2*r*full/2
        checks["P1_viscous_FTC_Z"+str(k)]=sp.cancel(sp.expand(sp.diff(pressure,z,k)))==0
        checks["m15_viscous_FTC_Z"+str(k)]=sp.cancel(sp.expand(sp.diff(m15,z,k)))==0
    require(all(checks.values()),"Physical pressure primitive algebra failed")
    return dict(passed=True,checks=checks)


def build(root,parent=None):
    root=Path(root);require(sha(root/PARENT)==PARENT_SHA256,"Accepted selected-word parent changed")
    saved=load(root/PARENT)
    if parent is None:parent=saved
    require(digest(parent)==digest(saved) and parent["checks"]["passed"] is True,"Same passed parent required")
    require(parent["actual_native_collar_transport"]["native_defining_leading_function_equals_canonical_collar"] is True and
        parent["actual_n1_cutoff_history_identity"]["N1p_units"]["passed"] is True and
        parent["actual_C_inner_m15_source_binding"]["actual_C_inner_m15_source_function_bound"] is True,
        "Actual native/canonical leading source and physical N1p moment units required")
    paths=(PARENT,MEAN,MEAN_CHECK,MEAN_SOURCE,DOMAIN,GRAPH,RAW_SOURCE,OP_SOURCE,CUT_SOURCE)
    pins={}
    for path in paths:
        expected=PARENT_SHA256 if path==PARENT else parent["input_hashes"][path[len(BASE):]]
        require(sha(root/path)==expected,"Changed pressure-primitive source: "+path)
        pins[path]=dict(path=path,expected_sha256=expected,passed=True)
    data,receipt,domain,graph=load(root/MEAN),load(root/MEAN_CHECK),load(root/DOMAIN),load(root/GRAPH)
    gate="current_original_leading_own_mean_and_radial_Q_epsilon_factors_installed"
    require(receipt["all_passed"] is True and receipt[gate] is True and
        data["source_family"]==receipt["source_family"]==parent["source_family"],
        "Accepted actual leading mean/Q source required")
    require(receipt["input_hashes"][Path(MEAN).name]==pins[MEAN]["expected_sha256"] and
        receipt["checks"]["uniform_original_tube_and_actual_inverse_width_orders_checked"] is True and
        receipt["checks"]["radial_leading_Uz_mean_Q_derivatives_retain_epsilon_pressure_sectors"] is True,
        "The accepted receipt must bind these actual endpoint rows and complex bounds")
    A,B=data["samples"]["actual_first_inlet"],data["samples"]["selected_actual_inner_exit"]
    require(receipt["checks"]["original_axis_and_identical_core_actual_first_inlet_factors"] is True and
        receipt["checks"]["genuine_three_n1_forcings_execute_from_exposed_leading_rows"] is True and
        A["Z"]==B["Z"]==
        parent["actual_selected_first_word_enclosure"]["selected_Z_interval"],
        "Same actual axial interval required")
    require(A["leading_Qaxis"]==B["leading_Qaxis"] and
        A["leading_radial_recovery_identity"]==B["leading_radial_recovery_identity"]==
        "Q=Qaxis+epsilon*Qpsi; Qpsi=[2Z*Psi-(1-delta)Z*meanPsi-d*meanPsi_Z]/L" and
        A["actual_own_mean_history_not_comparison_mean"] is True and
        B["actual_own_mean_history_not_comparison_mean"] is True,
        "Same own mean and exact radial-constant Qaxis source required")
    c=MPIntervalContext();c.prec=PRECISION_BITS
    z=interval(c,A["Z"])
    require(z.a>-1 and z.b<1,"Saved real axial interval inside (-1,1) required")
    s=interval(c,B["geometry"]["phase"])
    source_s=parent["actual_selected_first_word_enclosure"]["selected_coordinate_rows"][-1]["s_exact_selected_enclosure"]
    require(s==interval(c,source_s) and interval(c,A["geometry"]["phase"])==0,
        "Actual selected sc/4 endpoint required")
    nodes=graph["exact_graph_nodes"];params=graph["lower_collar_bindings"]["parameters"]
    require(params["delta"]==25220 and params["epsilon"]==35565 and params["Lambda"]==35563 and
        nodes[35563]==dict(operation="analytic_unary",name="exp",argument=35562) and
        nodes[35564]==dict(operation="negative",argument=35562) and
        nodes[35565]==dict(operation="analytic_unary",name="exp",argument=35564),
        "Actual epsilon*Lambda=1 source graph required")
    visited=set();delta=graph_scalar(c,nodes,25220,visited)
    require(delta>0 and delta<c.mpf("1e-200"),"Actual positive small delta required")
    geometry=parent["actual_native_collar_transport"]["domain"]
    require(geometry["passed"] and "0<hb<1/2" in geometry["rho_bound_proof"] and s>0 and s<c.mpf(1)/4096,
        "Actual positive width and microscopic phase required")
    guards=dict(
        factor=ast_guard(root,MEAN_SOURCE,"factored_leading_rows",(
            "qa=(u*(z*(1+physical.delta))-d*op.dz(u))/L",
            "qpsi=[(2*z*P-z*M*(1-physical.delta)-d*op.dz(M))/L",
            "q=[op.truncate(qa,n-1)+eps(op.truncate(Q,n-1))")),
        complex_bounds=ast_guard(root,MEAN_SOURCE,"complex_mean_bounds",(
            "uaxis=lift(u);qa=lift(((1+dt)*z*u+4*d)/L)",
            "q=[qa+eps(qpsi[0]),eps(qpsi[1]),eps(qpsi[2])]",
            "Qpsi_and_Q_bound_tube_radius=eta/4",
            "norm_bounds_not_selected_comparison_field=True,all_original_pressure_atoms_retained=True")),
        pressure=ast_guard(root,RAW_SOURCE,"raw_system",(
            "Np=a.add(a.scale('-1/2',S),a.scale(parameters['Lambda'],a.add(a.scale(2,c1(Q[1])),a.scale(r,c1(Q[2])))))",
            "g4=g.mul(c(2),eps,x,Np.value)")),
        m15=ast_guard(root,CUT_SOURCE,"signed_rho_integrands",(
            "m15=a.add(a.scale(a.mul(u,U),g.mul(g.constant(2),epsilon)),a.neg(a.scale(pressure,g.mul(eps2,r))))")),
    )
    axis=axis_jet(c,A,delta)
    norms=inertial_norms(c,B,domain);S1,S2=norm_polynomials(c,norms,delta)
    # rhoB-4=hb*Cd; integral_4^rhoB rho drho=hb*Cw.
    # exp(x)<=1/(1-x) gives outward exprel bounds for actual hb, without
    # selecting a numerical width or materializing exp(loghb).
    Cd=c.make_mpf(((4*s)._mpi_[0],(4*s/(1-s/2))._mpi_[1]))
    Cw=c.make_mpf(((16*s)._mpi_[0],(16*s/(1-s))._mpi_[1]))
    rows={};axis_phase_normalized={}
    for k in (0,1):
        ba,bb=poly(c,A["leading_Qpsi_rho_rows"][0],k),poly(c,B["leading_Qpsi_rho_rows"][0],k)
        ra,rb=poly(c,A["leading_Qpsi_rho_rows"][1],k),poly(c,B["leading_Qpsi_rho_rows"][1],k)
        db,dr=add(bb,scale(ba,-1)),add(rb,scale(ra,-1))
        viscP=shift(add(db,scale(dr,4)),(0,0,0,-1))
        viscP=add(viscP,shift(scale(rb,Cd),(1,0,0,-1)))
        viscM=shift(scale(dr,-16),(0,0,0,-2))
        viscM=add(viscM,shift(scale(rb,-2*Cw),(1,0,0,-2)))
        axisP={(1,0,0,-1):-Cd*axis[k]/2}
        axisM={(1,0,0,-2):Cw*axis[k]/2}
        axis_phase_normalized[str(k)]=dict(P1_W0=record(next(iter(axisP.values()))/s),
            m15_leading_N1p=record(next(iter(axisM.values()))/s))
        factor=1 if k==0 else 1/norms["gap"]
        pError=add(shift(scale(S1,Cd*factor/2),(1,0,0,-2)),
            shift(scale(S2,Cd*factor/2),(1,0,0,-3)))
        mError=add(shift(scale(S1,Cw*factor/2),(1,0,0,-3)),
            shift(scale(S2,Cw*factor/2),(1,0,0,-4)))
        rows[str(k)]=dict(
            P1_W0_local=dict(viscous_endpoint=export(viscP),axis_inertial=export(axisP),
                signed_approximation=export(add(viscP,axisP)),
                remaining_inertial_error_radius=export(pError,error=True)),
            m15_leading_N1p_local=dict(viscous_endpoint=export(viscM),axis_inertial=export(axisM),
                signed_approximation=export(add(viscM,axisM)),
                remaining_inertial_error_radius=export(mError,error=True)))
    algebra=exact_algebra()
    return dict(passed=True,source_family=parent["source_family"],parent_report=PARENT,
        parent_report_sha256=PARENT_SHA256,parent_mapping_sha256=digest(parent),
        source_pins=pins,source_AST=guards,exact_algebra=algebra,ordinary_Z_orders=[0,1],
        selected_Z_interval=A["Z"],selected_phase=record(s),ordinary_Z=rows,
        axis_inertial_phase_normalized_coefficients=axis_phase_normalized,
        exact_pressure_parameter_subgraph_nodes=sorted(visited),delta_enclosure=record(delta),
        radial_constant_Qaxis_canceled_as_same_function_before_boxing=True,
        own_mean_core_history_retained=True,
        geometry_difference_factored_before_endpoint_boxing=True,
        saved_endpoint_source_uncertainty_retained=True,
        source_log_basis_names=parent["actual_selected_first_word_enclosure"]["source_log_basis_names"],
        source_log_bases=parent["actual_selected_first_word_enclosure"]["source_log_bases"],
        radius_factors=dict(rhoB_minus4_over_hb=record(Cd),
            integral_rho_drho_over_hb=record(Cw),complex_S1_S2_derivative_gap=record(norms["gap"])),
        inertial_volume_norms=dict(S1=export(S1,error=True),S2=export(S2,error=True)),
        inertial_volume_bound_proof=dict(
            split="Q=a+epsilon*b, U=u+epsilon*v; S=Saxis+epsilon*S1+epsilon^2*S2",
            S1="A(b)+a*b+rho*a*b_rho+u*B(b)+v*B(a)",
            S2="b^2/2+rho*b*b_rho+v*B(b)",
            A="(q+Z*(1-delta)*q_Z/2+rho*q_rho)/L",
            B="(-2Z*q+(1-Z^2)*q_Z-2Z*rho*q_rho)/L",
            source_coefficient_tube="Qpsi eta/4, Psi eta/2; first q_Z loss eta/8",
            S1_S2_norm_tube="eta/8 on the real radial interval",
            ordinary_Z1_error="Cauchy on entire S1/S2 with gap eta/16",
            positive_polynomial_coefficients_are_error_bounds_not_signed_values=True),
        actual_signed_endpoint_and_axis_approximation_evaluated=True,
        actual_leading_N1p_integral_enclosed=True,
        full_W0_pressure_increment_enclosed=True,
        m15_component_is_only_leading_N1p=True,full_C_inner_m15_evaluated=False,
        residual_inertial_volume_signed_values_evaluated=False,
        endpoint_sign_or_useful_relative_accuracy_certified=False,
        all_Z_functional_coverage=False,main_source_exponential_factors_materialized=False,
        physical_N1p_no_Pstar_normalization_used=True,
        canonical_C_inner_source_handles=parent["actual_C_inner_m15_source_binding"]["canonical_C0_C3_handles"],
        source_component="-epsilon^2*rho*N1p in actual D15_rho; same component, not an added ledger term",
        P1_is_distinct_from_m15=True,
        old_error_caps_removed=False,six_term_ledger_modified=False,
        independent_terminal_matching_closed=False,global_five_moment_closure=False,
        accepted_ancestors_replayed=False,source_executed=False,source_only=False,
        precision_bits=PRECISION_BITS)
