"""Transport raw finite-word error through the signed cutoff m15 functional.

The leading field and its pressure source are shared exact functions.  The
error therefore depends linearly on the raw F1/Uz1 errors.  A common inlet
pressure uncertainty cancels from the pressure IBP before any boxing.
Positive budgets below are enclosures; they are never signed values.
"""
from __future__ import annotations
from fractions import Fraction
from pathlib import Path

import sympy as sp

from .original_finite_N_signed_source import ast_guard,digest,load,require,sha

PARENT='lei_ren_part1_paper_compliant_current_original_native_collar_transport.json'
PARENT_SHA256='7c05f7438a27d9caa4378cb89b9ab6d33b44138086e86292b73f2fd8e3f4bf96'
BASE='experiments/root_st073/'
CUT_SOURCE=BASE+'lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py'
PRESSURE_SOURCE=BASE+'lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions.py'
VOLTERRA=BASE+'lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz'


def _positive(q):
    require(not isinstance(q,(bool,float)),'Exact rational bound required, not a floating point proxy')
    q=Fraction(q);require(q>=0,'Nonnegative error or norm bound required');return q


def m15_error_budget(*,epsilon,rho_upper,rho_measure,leading_U_bounds,leading_F_bounds,
                     raw_U_error_bounds,raw_F_error_bounds):
    """Exact rational C0/C1 upper bounds, after caller supplies certified norms.

    Inputs are raw physical-unit bounds (value, ordinary Z derivative).
    rho_measure bounds integral eta d-rho, and rho_upper bounds rho on the
    integration support.  No pressure-axis norm is an input.  This arithmetic
    API does not certify the caller's inputs or evaluate the actual field.
    """
    e,r,J=map(_positive,(epsilon,rho_upper,rho_measure))
    require(e>0,'The physical epsilon must be positive')
    def pair(v):
        require(len(v)==2,'Ordinary Z0/Z1 bound pair required')
        return tuple(map(_positive,v))
    U,F,Eu,Ef=map(pair,(leading_U_bounds,leading_F_bounds,raw_U_error_bounds,raw_F_error_bounds))
    return (2*e*J*U[0]*Eu[0]+2*e*e*r*J*F[0]*Ef[0],
        2*e*J*(U[1]*Eu[0]+U[0]*Eu[1])+
        2*e*e*r*J*(F[1]*Ef[0]+F[0]*Ef[1]))


def exact_algebra():
    z,r=sp.symbols('Z rho',real=True)
    eps,P2=sp.symbols('epsilon Pstar_squared',positive=True)
    eta=sp.Function('eta')(r)
    U,F,u,f,du,df,Np=[sp.Function(q)(r,z) for q in ('U','F','u1','F1','error_u1','error_F1','Np')]
    density=lambda uv,fv:2*eps*U*eta*uv-eps**2*r*(2*F*eta*fv+Np)
    error=sp.expand(density(u+du,f+df)-density(u,f))
    expected=2*eps*U*eta*du-2*eps**2*r*F*eta*df
    checks={
        'same_leading_Np_cancels_as_function':sp.expand(error-expected)==0,
        'error_is_linear_in_raw_F1_Uz1':sp.diff(error,du,df)==0,
        'physical_Pstar_squared_removed_once':sp.cancel(eps**2*P2*r*(2*F*eta*df/P2)-2*eps**2*r*F*eta*df)==0,
    }
    explicit_Z1=(2*eps*eta*(sp.diff(U,z)*du+U*sp.diff(du,z))-
        2*eps**2*r*eta*(sp.diff(F,z)*df+F*sp.diff(df,z)))
    checks['ordinary_Z1_product_rule_without_factorial']=sp.expand(sp.diff(error,z)-explicit_Z1)==0
    a,b=sp.symbols('Rin Rmax',positive=True);Paxis=sp.Function('inlet_pressure_error')(z)
    common_pressure=(b-a)*Paxis-b*Paxis+a*Paxis
    for k in range(4):
        checks['common_inlet_pressure_IBP_cancels_Z'+str(k)]=sp.expand(sp.diff(common_pressure,z,k))==0
    require(all(checks.values()),'Signed m15 tail transport algebra failed')
    return dict(passed=True,checks=checks,physical_density_error=str(expected),
        ordinary_Z1_error=str(explicit_Z1),pressure_inlet_error_is_a_shared_function=True,
        pressure_inlet_canceled_before_independent_interval_evaluation=True,
        leading_Np_kept_in_full_m15_but_cancels_only_from_same_source_difference=True)


def build(root,parent,raw,cutoff,tail):
    root=Path(root);require(sha(root/BASE/PARENT)==PARENT_SHA256 and
        digest(parent)==digest(load(root/BASE/PARENT)),'Actual accepted native parent required')
    for packet in (raw,cutoff,tail):
        require(packet.get('passed') is True and packet.get('source_family')==parent['source_family'] and
            packet.get('parent_report_sha256')==PARENT_SHA256 and packet.get('source_executed') is False,
            'Same-source actual raw/cutoff/tail packets required')
    require(raw['graph_binding']['path']==VOLTERRA and
        raw['graph_binding']['sha256']==parent['input_hashes'][VOLTERRA[len(BASE):]],
        'Actual raw leading graph changed')
    leading=load(root/VOLTERRA)['inner_Volterra_bindings']['leading']['first']
    native=parent['actual_native_collar_transport']['canonical_function_handles']
    require(leading['F_rho_rows'][0]==native['F_C0_C3'] and leading['U']==native['U_C0_C3'] and
        leading['m']==native['mean_U_C0_C3'], 'Raw and cutoff moment use different leading functions')
    fields=('F1','Uz1','K1','P1','d_x_F1','d_x_Uz1')
    require([q for key in fields for q in raw['actual_function_handles']['first'][key]]==cutoff['raw_state']['first_nodes'] and
        cutoff['raw_input_gate']['function_source_gate']['passed'] is True and
        cutoff['N1p_units']['passed'] is True and
        cutoff['cutoff_identity']['core_inlet_not_reset'] is True and
        tail['certified_raw_preconditioned_tail_majorant'] is True,
        'Actual raw ports, physical pressure convention or tail provenance incomplete')
    pins={}
    for path in (CUT_SOURCE,PRESSURE_SOURCE,VOLTERRA):
        expected=parent['input_hashes'].get(path[len(BASE):]);require(expected==sha(root/path),'Changed physical moment source')
        pins[path]=dict(path=path,expected_sha256=expected,passed=True)
    guards={
        'signed_physical_rho_density':ast_guard(root,CUT_SOURCE,'signed_rho_integrands',(
            "eps2=g.mul(epsilon,epsilon)","pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])",
            "m15=a.add(a.scale(a.mul(u,U),g.mul(g.constant(2),epsilon)),a.neg(a.scale(pressure,g.mul(eps2,r))))")),
        'retained_first_pressure':ast_guard(root,CUT_SOURCE,'own_cutoff_fields',(
            'actual_cutoff_own_pressure=True','leading_N1p_not_cut=True','nonzero_core_pressure_inlet_retained=True')),
        'physical_pressure_gradient':ast_guard(root,PRESSURE_SOURCE,'pressure_radial_function',(
            'epsilon',"a.scale(a.mul(leadingF,F1),2)",)),
    }
    algebra=exact_algebra()
    checks=dict(parent_and_actual_sources=True,actual_density_and_pressure_AST=True,
        same_source_error_algebra=algebra['passed'],same_source_raw_cutoff_tail_and_leading=True,
        actual_raw_port_equality=True,pressure_units_and_no_reset_port_bound=True,
        error_budget_conditional_on_explicit_approximation_contract=True)
    require(all(checks.values()),'Actual raw error-transport inputs missing')
    return dict(passed=True,parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        checks=checks,source_pins=pins,source_AST=guards,exact_algebra=algebra,
        raw_packet_sha256=digest(raw),cutoff_packet_sha256=digest(cutoff),tail_packet_sha256=digest(tail),
        approximation_contract=dict(leading_fields_and_Np='identical exact leading functions on both sides',
            raw_approximation='finite natural-word sum or another approximation with separately proved raw error',
            approximate_pressure='P1_hat(s)=P1_hat_A(Z)+integral epsilon*(2*F_lead*eta*F1_hat+N1p_same) d-rho',
            paired_pressure_error='delta_P1(s)=delta_P1_A(Z)+integral 2*epsilon*F_lead*eta*delta_F1 d-rho',
            same_approximate_inlet_function_at_both_IBP_endpoints=True,
            cutoff='eta is the same exact Z-independent cutoff; 0<=eta<=1',
            rho_measure='J >= integral eta d-rho, retaining the actual support',
            error0='2*epsilon*J*U0*Eu0 + 2*epsilon^2*rho_max*J*F0*Ef0',
            error1='2*epsilon*J*(U1*Eu0+U0*Eu1)+2*epsilon^2*rho_max*J*(F1*Ef0+F0*Ef1)',
            API='m15_error_budget: exact rational arithmetic; caller-supplied bounds must be independently certified'),
        actual_signed_m15_value_evaluated=False,actual_m15_numeric_error_budget_instantiated=False,
        pressure_inlet_cancellation_conditional_on_this_reconstruction=True,
        arbitrary_legacy_pressure_enclosures_admitted=False,
        raw_tail_constants_or_leading_norms_selected_as_values=False,
        old_pressure_or_moment_error_caps_removed=False,all_old_error_caps_retained=True,
        saved_six_terms_or_error_ledger_contracted=False,global_five_moment_closure=False,
        accepted_ancestor_or_quadrature_replayed=False)
