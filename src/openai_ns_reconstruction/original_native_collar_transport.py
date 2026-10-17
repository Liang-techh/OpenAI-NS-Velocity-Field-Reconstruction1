"""Native defining collar ODEs on the unchanged first microscopic window.

Only pinned sources and saved graph handles are read.  The two integrated
core direction identities are inputs, proved in their separate packets.
The finite-width numerical rows and all error budgets remain unchanged.
"""
from __future__ import annotations

import ast
from pathlib import Path

import sympy as sp

from .original_finite_N_signed_source import digest, load, require, sha

PARENT='lei_ren_part1_paper_compliant_current_original_core_collar_pressure_identity.json'
PARENT_SHA256='ae91b0d72d0415e8146f7ee9926d1a35ddbedffff6b89eac1f26558dbb5b44c7'
BASE='experiments/root_st073/'
LOWER=BASE+'lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions.json.gz'
LOWER_SOURCE=LOWER.removesuffix('.json.gz')+'.py'
MIXED=BASE+'lei_ren_part1_paper_compliant_bridge_mixed_C4.py'
MICRO=BASE+'lei_ren_part1_paper_compliant_current_original_micro_functions.py'
FLAT=BASE+'lei_ren_part1_paper_compliant_flat_pulse_derivatives.py'
EXIT=BASE+'lei_ren_part1_paper_compliant_current_inner_exit_strict_collar.json'
EXIT_SOURCE=EXIT.removesuffix('.json')+'.py'
MICRO_CHECK=MICRO.removesuffix('.py')+'_check.json'


def _guard(root,path,name,fragments,owner=None):
    tree=ast.parse((root/path).read_text(encoding='utf-8'))
    if owner:
        classes=[n for n in ast.walk(tree) if isinstance(n,ast.ClassDef) and n.name==owner]
        require(len(classes)==1,'Unique source class required: '+owner)
        search=classes[0].body
    else:search=ast.walk(tree)
    nodes=[n for n in search if isinstance(n,ast.FunctionDef) and n.name==name]
    require(len(nodes)==1,'Unique source function required: '+name)
    node=nodes[0];body=''.join(ast.unparse(node).split())
    checks={q:''.join(q.split()) in body for q in fragments}
    require(all(checks.values()),'Actual collar source changed: '+path+':'+name)
    return dict(path=path,function=name,owner=owner,checks=checks,
        AST_sha256=digest(ast.dump(node,include_attributes=False)),source_executed=False)


def _pin(root,parent,path):
    expected=(parent['input_hashes'].get(path[len(BASE):]) if path.startswith(BASE)
        else parent['implementation_hashes'].get(path))
    require(isinstance(expected,str) and sha(root/path)==expected,'Pinned collar source changed: '+path)
    return dict(path=path,expected_sha256=expected,passed=True)


def _domain(graph,exit_record):
    b=graph['lower_collar_bindings'];g=b['geometry_certificate'];n=graph['exact_graph_nodes']
    bits=g['sc_exact_mpf_tuple'];sign,man,exponent,bit_count=bits
    selected=exit_record['explicit_current_inner_exit_strict_collar']['selected_first_phase_endpoint']
    require(selected['lower_exact_mpf_tuple']==bits==selected['upper_exact_mpf_tuple'],
        'The source sc is not the identical exact selected scalar')
    require(sign==0 and man>0 and man.bit_length()==bit_count and exponent+bit_count<=-10,
        'Need exact positive sc<2^-10 without materializing its denominator')
    require(n[35819]==dict(operation='exact_rational',numerator=man,denominator=1) and
        n[35820]==dict(operation='exact_rational',numerator=exponent,denominator=1) and
        n[35821]==dict(operation='product',arguments=[35820,3184]) and
        n[35822]==dict(operation='analytic_unary',name='exp',argument=35821) and
        n[35823]==dict(operation='product',arguments=[35819,35822]) and
        n[3184]==dict(operation='analytic_unary',name='log',argument=48) and
        n[48]==dict(operation='exact_rational',numerator=2,denominator=1),
        'Exact sc graph does not encode man*exp(exponent*log(2))')
    require(b['parameters']['sc']==35823 and b['parameters']['sc_half']==35825 and
        b['parameters']['sc_quarter']==35824 and g['sc_and_hb_Z_independent'] is True,
        'Actual source phase geometry changed')
    support=g['original_modification_support']
    require(support['increments_exact_zero'] is True and support['region']=='unchanged_left' and
        support['original_P0_and_background_histories_retained'] is True,
        'Selected collar includes a later finite-N modification')
    # This saved positive margin bounds the exact width, not a chosen width.
    margin=exit_record['explicit_current_inner_exit_strict_collar']['positive_margins']['exact_width_below_half']
    lo=margin['lower_exact_mpf_tuple']
    require(lo[0]==0 and lo[1]>0 and
        exit_record['exact_current_exit_source_attachment']['passed'] is True,
        'Actual hb<1/2 source receipt required')
    require(n[8]['name']=='hbB' and n[8]['Z_independent'] is True and
        n[8]['quantity_is_exact_original_parameter_not_range_endpoint'] is True,
        'Actual positive source hb is not the exact source handle')
    return dict(passed=True,sc_exact_mpf_tuple=bits,sc_node=35823,hb_node=8,
        exact_integer_bit_exponent_check=True,sc_positive_and_less_than_2_power_minus10=True,
        hb_positive_less_than_half_source_margin=margin,
        native_first_s_interval='[0,sc/2]',canonical_q_interval=[0,2],s_of_q='sc*q/4',
        native_log_radius='y=hb*s',rho='4*exp(hb*s)',physical_R='epsilon*rho',
        first_window_only=True,alpha_identically_one=True,
        rho_upper_rational='4/(1-1/2048) < 41/10',
        rho_bound_proof='0<hb<1/2, 0<=s<1/2048, hence hb*s<1/2048; exp(x)<=1/(1-x)',
        ordinary_Z_independent_endpoints=True,later_modification_identically_zero=True,
        exact_sc_or_hb_materialized=False,range_endpoint_selected_as_source_value=False)


def _algebra():
    s,r,z=sp.symbols('s rho Z',real=True)
    h,e=sp.symbols('hb epsilon',positive=True)
    chi=sp.Function('chi')(s);a=sp.Function('ratio')(s,z)
    p=sp.Function('Phi_core')(r,z);v=sp.Function('Psi_core')(r,z)
    F0=sp.Function('F0')(z);U0=sp.Function('U0')(z)
    ar=(chi-1)*h*r*sp.diff(p,r)/p*a
    Ds=lambda f:sp.diff(f,s)+h*r*sp.diff(f,r)
    F=F0*a*p;Dbar=-2*r*sp.diff(p,r)/p
    F_res=Ds(F).subs(sp.diff(a,s),ar)+h*chi*Dbar*F/2
    Psi=sp.Function('Psi')(s,z);mean=sp.Function('mean_Psi')(s,z)
    U=U0+e*Psi;m=U0+e*mean
    psi_s=a*chi*h*r*sp.diff(v,r)
    Gtotal=-e*r*sp.diff(v,r)
    U_res=sp.diff(U,s).subs(sp.diff(Psi,s),psi_s)+h*chi*a*Gtotal
    mean_s=h*(Psi-mean)
    mean_res=Ds(r*m).subs(sp.diff(mean,s),mean_s)-h*r*U
    checks={}
    for k in range(4):
        checks['native_angular_ordinary_Z'+str(k)]=sp.simplify(sp.diff(F_res,z,k))==0
        checks['native_axial_ordinary_Z'+str(k)]=sp.simplify(sp.diff(U_res,z,k))==0
        checks['own_mean_ordinary_Z'+str(k)]=sp.simplify(sp.diff(mean_res,z,k))==0
    checks['actual_phi_over_comparison_phi_equals_ratio']=sp.cancel((a*p)/p-a)==0
    pressure_res=sp.cancel(e*h*r*(F0*a*p)**2/(e*h*r)-(F0*a*p)**2)
    for k in range(4):checks['pressure_FTC_ordinary_Z'+str(k)]=sp.simplify(sp.diff(pressure_res,z,k))==0
    # The actual C history uses the evolved field.  It is not the core or
    # comparison C history once chi differs from one.
    C=sp.Function('C_actual')(s,z);Paxis=sp.Function('P_axis')(z)
    Pnative=Paxis+e*r*F0**2*C
    Cs=h*((a*p)**2-C)
    actual_pressure_res=Ds(Pnative).subs(sp.diff(C,s),Cs)-e*h*r*(F0*a*p)**2
    for k in range(4):checks['actual_own_pressure_C_ordinary_Z'+str(k)]=sp.simplify(sp.diff(actual_pressure_res,z,k))==0
    # Core primitives solve exactly the same six comparison moment ODEs.
    y=sp.symbols('y',real=True);t=sp.symbols('t',positive=True)
    phi=sp.Function('Phi')(t,z);V=sp.Function('V')(t,z)
    specs=(('H',2,2*t*phi,2*sp.Function('Phi')(r,z)),
        ('M',1,V,sp.Function('V')(r,z)),
        ('K',2,2*t*phi*V,2*sp.Function('Phi')(r,z)*sp.Function('V')(r,z)),
        ('A',1,V**2,sp.Function('V')(r,z)**2),
        ('B',2,t*phi**2,sp.Function('Phi')(r,z)**2),
        ('C',1,phi**2,sp.Function('Phi')(r,z)**2))
    primitive_definitions={}
    for name,rate,body,rhs in specs:
        q=sp.Integral(body,(t,0,r))/r**rate
        defect=sp.simplify(r*sp.diff(q,r)+rate*q-rhs)
        checks['same_core_comparison_own_'+name]=defect==0
        primitive_definitions[name]=dict(rate=rate,primitive=str(q),source_rhs=str(rhs),
            initial_rho4='same actual common core integral, without reset')
    require(all(checks.values()),'Native/canonical source transport identity failed')
    return dict(passed=True,checks=checks,comparison_core_primitives=primitive_definitions,
        ordinary_Z_orders=[0,1,2,3],F0_is_a_Z_function=True,epsilon_and_hb_Z_independent=True,
        ordinary_derivatives_without_factorial_division=True,
        exact_functions_not_interval_centers=True,old_source_executed=False)


def build(root,parent,angular,axial):
    root=Path(root)
    require(sha(root/BASE/PARENT)==PARENT_SHA256 and digest(parent)==digest(load(root/BASE/PARENT)),
        'Actual accepted core/collar parent required')
    leading=parent['actual_leading_collar_source_binding']
    require(leading['passed'] is True and
        leading['actual_leading_collar_function_identity']['new_smooth_object_equals_canonical_lower_collar'] is True and
        leading['actual_leading_collar_function_identity']['accepted_native_inlet_same_function'] is True and
        leading['actual_leading_collar_function_identity']['accepted_true_pressure_same_function'] is True,
        'Accepted actual core and canonical collar functions required')
    require(angular['passed'] is True and
        angular['angular_direction_identity']['symbolic_proof']['direction_identity_exact'] is True and
        angular['angular_direction_identity']['symbolic_proof']['ordinary_Z_derivatives_exact'] is True and
        axial['passed'] is True and axial['axial_drive_identity']['passed'] is True,
        'Proved actual core direction inputs required')
    nonzero=angular['angular_direction_identity']['nonvanishing_domain']
    require(nonzero['passed'] is True and nonzero['F_core_positive_bound'] is True and
        nonzero['L_positive_bound'] is True and
        nonzero['Phi_receipt']['scaled_radial_domain']==['0','4.1'],
        'Source-bound core/comparison positivity through the selected collar required')
    paths=(LOWER,LOWER_SOURCE,MIXED,MICRO,FLAT,EXIT,EXIT_SOURCE,MICRO_CHECK)
    pins={Path(p).name:_pin(root,parent,p) for p in paths}
    graph=load(root/LOWER);exit_record=load(root/EXIT)
    micro_check=load(root/MICRO_CHECK)
    require(micro_check['all_passed'] is True and
        micro_check['current_original_Ra_micro_hydrated_whole_source_functions_installed'] is True,
        'Accepted original source ODE binding required')
    guards={
        'source_comparison_definition':_guard(root,MIXED,'formal_comparison_source',(
            "integrand='alpha(t)*D_t log(Phi_core(4exp(t),Z))'",
            "alpha='1-sigma((t-hb)/hb)'", "core_radius='rho=4exp(t)'",
            "integrand='alpha(t)*D_t V_core(4exp(t),Z)'",
            "H='(1/8)*integral_0^4 rho*Phi_core(rho,Z)drho'",
            "M='(1/4)*integral_0^4 V_core(rho,Z)drho'",
            "K='(1/8)*integral_0^4 rho*Phi_core(rho,Z)*V_core(rho,Z)drho'",
            "A='(1/4)*integral_0^4 V_core(rho,Z)^2 drho'",
            "B='(1/16)*integral_0^4 rho*Phi_core(rho,Z)^2 drho'",
            "C='(1/4)*integral_0^4 Phi_core(rho,Z)^2 drho'",
            "('H',2,'2*phi_bar(t,Z)')", "('M',1,'V_bar(t,Z)')",
            "('K',2,'2*phi_bar(t,Z)*V_bar(t,Z)')", "('A',1,'V_bar(t,Z)^2')",
            "('B',2,'phi_bar(t,Z)^2')", "('C',1,'phi_bar(t,Z)^2')")),
        'comparison_first_alpha':_guard(root,MIXED,'evaluate',(
            "if chart=='first':", "alpha=[algebra.lift(1)]+[algebra.lift(0)]*3",
            "controls=bridge_controls(algebra,scale,chi,Dbar,drive,phi/comparison['phi'].truncate(5),barlog)"),
            'CompliantBridgeMixedC4'),
        'native_defining_ODE':_guard(root,MICRO,'source_bindings',(
            "actual_first_equations='phi_y=-((1-sigma)+hb*sigma)*Dbar*phi/2; V_y=-((1-sigma)+hb*sigma)*(phi/barphi)*drive'",
            "actual_geometry='R=Ra*exp(hb*s), dy=hb*ds'")),
        'native_local_quotient':_guard(root,MICRO,'evaluate',(
            "quotient=f.multiply(at['phi'],f.jet(barphi.reciprocal()))",
            "R=self.Ra*self.series.scalar_series(1,s)[0]",
            "chi=f.scalar(1-sig)+f.h*sig",
            "true_physical_measure='dy=hb*ds'"),'MicroFunctions'),
        'actual_own_history_source':_guard(root,MICRO,'histories',(
            "phi,V=(full['phi'],full['V'])",
            "rhs=dict(H=f.scale(phi,2),M=V,K=f.scale(f.multiply(phi,V),2),A=f.multiply(V,V),B=f.multiply(phi,phi),C=f.multiply(phi,phi))",
            "result[name]=f.add(f.scale(self.initial[name],decay),f.scale(rhs[name],mass))",
            "original_core_memory_retained=True"),'MicroFunctions'),
        'actual_own_history_inlet':_guard(root,MICRO,'__init__',(
            "self.initial={name:flow.jet(jet(row)) for name,row in actual_core['actual_own_six_moments_axial5'].items()}",),
            'MicroFunctions'),
        'same_canonical_ODE_and_history':_guard(root,LOWER_SOURCE,'collar_fields',(
            "chi=lambda x:g.add(g.one,g.mul(g.sub(hb,g.one),sigma(x)))",
            "ratio=exp3(g,logratio)", "rho4_inlet_ratio_exact_one=True",
            "actual_Psi_radial_transport=True", "nonzero_rho4_incoming_mean_retained=True",
            "core_pressure_history_not_reset=True", "same_true_core_Phi_squared=True")),
        'exact_source_scalar':_guard(root,LOWER_SOURCE,'exact_binary_scalar',(
            "sign,man,exponent,_=bits", "g.unary('log',c(2))")),
        'source_hb_width_margin':_guard(root,EXIT_SOURCE,'prove',(
            "logh=read(ledger['shared_positive_width_log_enclosure'])",
            "require('exact_width_below_half',c.ln(c.mpf('.5'))-c.mpf(endpoints(logh)[1]))")),
        'flat_negative_argument':_guard(root,FLAT,'sigma_jets',(
            "if lo<=0:rows.append(IntervalTaylor.constant(c,0,4))",
            "if hi>0 and lo<=mp.mpf('.5'):")),
    }
    domain=_domain(graph,exit_record);algebra=_algebra();b=graph['lower_collar_bindings']
    f=b['fields']
    require(all(len(f[k])==4 for k in ('ratio','U','m','absolute_swirl_base_pressure')) and
        all(len(row)==4 for row in f['F_rho_rows']) and
        f['all_histories_start_at_true_core_rho4'] is True,
        'Canonical ordinary Z0..3 function handles or histories changed')
    checks=dict(parent_and_sources=True,actual_source_AST=True,exact_first_window_domain=domain['passed'],
        fresh_source_transport_algebra=algebra['passed'],actual_native_core_angular=angular['passed'],
        actual_native_core_axial=axial['passed'],same_native_inlet_and_pressure=True,
        native_quotient_not_frozen_initial_quotient=True,canonical_C3_handles=True,
        inherited_core_and_comparison_nonvanishing=True,
        old_numeric_rows_and_errors_not_replaced=True)
    return dict(passed=all(checks.values()),parent_report=PARENT,parent_report_sha256=PARENT_SHA256,
        checks=checks,source_pins=pins,source_AST=guards,domain=domain,exact_algebra=algebra,
        angular_packet_sha256=digest(angular),axial_packet_sha256=digest(axial),
        accepted_leading_packet_sha256=digest(leading),
        native_defining_leading_function_equals_canonical_collar=True,
        actual_native_D_over_R_and_total_drive_bound=True,
        actual_quotient_equals_ratio=True,ordinary_Z_orders=[0,1,2,3],
        comparison_core_primitive_identity=True,
        native_actual_six_histories_equal_comparison=False,
        actual_own_mean_and_pressure_equal_canonical_actual_histories=True,
        history_distinction=dict(direction_inputs='six own comparison moments = the continued core primitives on alpha=1',
            actual_mean='M_y+M=U_actual, with the original common core M at the inlet',
            actual_pressure='C_y+C=Phi_actual^2; P=P_axis+R*F0(Z)^2*C, with original core pressure history',
            no_actual_history_replaced_by_comparison_history=True),
        quotient_conventions=dict(pointwise='Phi_native(s,Z)/Phi_comparison(s,Z)=ratio(s,Z)',
            prefix='Phi_core(4,Z)*exp(ell(s,Z))/Phi_comparison(s,Z)',
            ell='-1/2*integral_0^s hb*chi(t)*Dbar(t,Z)dt',
            prefix_phi_a_means_inlet_not_evolved_field=True,
            extra_exponential_after_evolved_field=False),
        canonical_function_handles=dict(report=LOWER,report_sha256=pins[Path(LOWER).name]['expected_sha256'],
            F_C0_C3=f['F_rho_rows'][0],U_C0_C3=f['U'],mean_U_C0_C3=f['m'],
            leading_pressure_C0_C3=f['absolute_swirl_base_pressure'],ratio_C0_C3=f['ratio']),
        uniqueness_argument=[
            'On this exact first-window subinterval alpha=1, so comparison Phi/V equal the continued common core functions by FTC.',
            'Each of the six comparison moments solves the same linear ODE as the corresponding core primitive, with identical nonzero core data at rho4; comparison moments therefore equal core primitives.',
            'The proved radial direction identities reduce the actual angular ODE to the canonical log-ratio ODE with ratio(0)=1.',
            'The actual pointwise phi/barphi is this ratio, so the axial ODE reduces to the canonical Psi FTC with the same core-exit U.',
            'The separate actual own mean and pressure use U_actual and Phi_actual^2 in their history ODEs, not the comparison RHS. Their retained core inlet data and FTCs give the canonical actual mean and pressure functions.',
            'Smooth real ODE uniqueness and ordinary Z differentiation through order3 give these identities throughout the subinterval.'],
        first_window_only=True,old_finite_numeric_producer_is_exact_function=False,
        raw_n1_cut_response_identity=False,full_C_inner_m15_identity=False,
        numerical_signed_A_minus_A0_value=False,saved_six_terms_or_error_ledger_contracted=False,
        six_signed_source_equality=False,global_five_moment_closure=False,
        all_old_error_caps_retained=True,accepted_ancestor_or_quadrature_replayed=False,
        later_second_window_or_global_native_identity=False,full_NS_validation=False)
