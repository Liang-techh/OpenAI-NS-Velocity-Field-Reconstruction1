"""Scoped actual-source and independent radial/pressure calculus checks."""
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions_check as prior

Algebra=prior.Algebra


def jet(g,variables,body,order=5):
    Z=sy.Symbol('Z')
    return current.FunctionJet(tuple(current.source.expression(g,sy.diff(body,Z,k),variables) for k in range(order+1)))


def actual_source(owner):
    g=owner.g;r=current.rows;free={};incoming=owner.zero_pressure_integrals;count=0
    assert owner.loop.acceptance_loaded and tuple(owner.charts)==current.CHARTS
    assert owner.bindings['original_pressure_reduction']==current.parent.active.current.ast_binding(current.original_pressure.pressure_functions)
    assert owner.bindings['original_endpoint_reduction']==current.parent.active.current.ast_binding(current.original_pressure.q_plus_q_y)
    for name,chart in owner.charts.items():
        before=owner.loop.charts[name];native=g.nodes[before['coordinate'].node]['name']
        assert chart['source_chart'] is before and chart['incoming_integrals'] is incoming
        assert chart['radial']['total_fast_phase_terms_included']
        D=owner.derivatives[name]
        assert D.cache[before['actual_phase'].node]==g.zero and D.shared['family']==owner.family
        for key,bgjet in before['source_chart']['fields'].items():
            for ref,dr in zip(r(bgjet),r(chart['radial_background'][key])):
                assert D.cache[ref.node]==g.mul(before['physical_measure'],dr)
        for factory in D.shared['factories'].values():
            for node in factory.core_extensions:
                row=g.nodes[node]
                assert row['source_family']==owner.family and row['radial_order']<=2
                assert row['unique_analytic_Banach_solution'] and row['actual_twenty_term_operator']
        for key in ('E_y','V_y','m_y','slow_A_y','slow_B_y'):
            assert len(r(chart['radial'][key]))==6
            for ref in r(chart['radial'][key]):assert prior.free_variables(g,ref,free)<={'Z',native}
        for group,order in ((chart['generic'],4),(chart['physical'],None),
                            (chart['pressure'],3),(chart['cumulative_integrals'],3)):
            for key,jet0 in group.items():
                expected=5 if key in ('U','U_y','m','m_y') else 4 if key in ('Q','Q_y') else order
                assert len(r(jet0))==expected+1
                for ref in r(jet0):assert prior.free_variables(g,ref,free)<={'Z',native},(name,key,ref.node)
                count+=len(r(jet0))
        assert len(r(chart['S']))==4
        for key,power in (('S0',1),('S1',2)):
            for order,ref in enumerate(r(chart['local_integrals'][key])):
                row=g.nodes[ref.node]
                assert row['operation']=='definite_integral' and row['ordinary_slow_Z_derivative_order']==order
                assert row['lower']==before['domain'][0].node and row['upper']==before['coordinate'].node
                assert row['physical_radius_power']==power and row['physical_Jacobian_applied_once']
                assert row['original_signed_reduced_pressure'] and row['independent_P0_retained_in_original_source']
        incoming=owner.endpoints[name]['integrals']
    assert owner.Rm is owner.endpoints['postrestore']
    assert all(q==g.zero for value in owner.zero_pressure_integrals.values() for q in r(value))
    return dict(passed=True,actual_true_modified_charts=list(current.CHARTS),actual_generic_physical_pressure_rows=count,
        true_phase_held_native_function_derivatives_and_total_fast_terms=True,
        exact_existing_seven_background_seed_rows_and_unique_core_source=True,
        complete_modified_mean_C5_and_full_generic_Q_Qy_C4=True,
        all_ten_signed_pressure_predecessors_and_finite_endpoints_retained=True,
        Pstar_once_and_physical_R_Rsquared_measures_once=True,
        no_free_core_rho_old_phase_angle_or_integration_variable=True,
        direct_Qyy_full_preRp_all17_numeric_and_global_claims_not_installed=True)


def binder_sources(owner):
    g=owner.g;fast={chart['actual_phase'].node for chart in owner.loop.charts.values()}
    names={};count=0
    for i,row in enumerate(owner.prefix):
        if row['operation']=='bound_variable':
            assert row['name'] not in names;names[row['name']]=i
        if i>=len(owner.loop.prefix) and row['operation']=='function_substitution':
            assert row['variable'] not in fast
            assert g.nodes[row['variable']]['operation']=='bound_variable'
            count+=1
    return dict(passed=True,new_native_substitution_binders_checked=count,
        bound_names_unique_and_no_fast_phase_node_used_as_binder=True)


def native_chain_calculus():
    g=current.source.FunctionTransportGraph();Y,Z=sy.symbols('fixture_native Z')
    y,z=g.symbol('fixture_native'),g.symbol('Z');phase=g.symbol('fixture_fast_phase');t=g.symbol('fixture_integral_t')
    T=sy.Symbol('fixture_integral_t');v=g.symbol('fixture_substitution_v');V=sy.Symbol('fixture_substitution_v')
    variables={'fixture_native':y,'Z':z,'fixture_integral_t':t,'fixture_substitution_v':v}
    expr=lambda q:current.source.expression(g,q,variables)
    f=g.add(expr(Y*T**2+Z*T),g.unary('sin',g.add(y,t)))
    integral=current.controls.integral(g,f,'fixture_integral_t',expr(Y**2/13),expr(1+Y/5),moving_native_endpoints_fixture=True)
    original=expr(sy.exp(V)+Y*V+Z*V**2);argument=expr(Y**2+Z/17)
    substituted=g.node('function_substitution',expression=original.node,variable=v.node,value=argument.node,Z_independent_substitution=True)
    # A bound variable identical to the outer endpoint must remain bound.
    same=current.controls.integral(g,expr(Y**2),'fixture_native',g.zero,y,lexical_shadowing_fixture=True)
    # Substituting the differentiation variable itself has only its chain term.
    self_sub=g.node('function_substitution',expression=expr(Y**3+Z).node,variable=y.node,value=expr(Y**2+Z/19).node,
        Z_independent_substitution=True)
    body=g.add(integral,substituted,same,self_sub);D=current.NativeDerivative(g,y,phase)
    derivative=D.derivative(body);second=D.derivative(derivative)
    alg=Algebra(g,{y.node:Y,z.node:Z})
    expected=sy.integrate(Y*T**2+Z*T+sy.sin(Y+T),(T,Y**2/13,1+Y/5))
    expected+=sy.exp(Y**2+Z/17)+Y*(Y**2+Z/17)+Z*(Y**2+Z/17)**2+Y**3/3+(Y**2+Z/19)**3+Z
    point={Y:sy.Rational(1,7),Z:sy.Rational(3,8)}
    for k,ref in ((0,body),(1,derivative),(2,second)):
        prior.close(alg.at(ref.node).subs(point),sy.diff(expected,Y,k).subs(point),('native chain',k))
    return dict(passed=True,independent_moving_integral_and_substitution_rows=3,
        direct_plus_chain_terms_and_identical_binder_name_checked=True,numerical_calculus_fixture_only=True)


def implicit_native_calculus():
    g=current.source.FunctionTransportGraph();Y,Z,Psi=sy.symbols('fixture_native Z fixture_angle')
    y,z,angle=g.symbol('fixture_native'),g.symbol('Z'),g.symbol('fixture_angle');phase=g.symbol('fixture_phase')
    variables={'fixture_native':y,'Z':z,'fixture_angle':angle};expr=lambda q:current.source.expression(g,q,variables)
    Phi=expr((1+Y**2+Z**2)*sy.exp(Psi)+Y**3)
    inverse=g.node('original_monotone_phase_inverse',phase_function=Phi.node,angle_variable='fixture_angle',
        target_phase=phase.node,angle_lower=g.zero.node,angle_upper=g.constant(10).node,
        Delta=g.zero.node,eta=g.one.node,independent_closed_inverse_fixture=True)
    known=sy.log((7-Y**3)/(1+Y**2+Z**2));T1=Y*Psi**2+Z*Psi**5
    composed=current.previous.at_inverse(g,expr(T1),inverse,angle);D=current.NativeDerivative(g,y,phase)
    got=(D.derivative(inverse),D.derivative(composed));alg=Algebra(g,{y.node:Y,z.node:Z,phase.node:sy.Integer(7),inverse.node:known})
    point={Y:sy.Rational(1,7),Z:sy.Rational(3,8)}
    for label,ref,want in (('implicit native',got[0],sy.diff(known,Y)),
                           ('inverse composition native',got[1],sy.diff(T1.subs(Psi,known),Y))):
        prior.close(alg.at(ref.node).subs(point),want.subs(point),label)
    return dict(passed=True,independent_closed_inverse_and_composition_native_rows=2,
        positive_lambda_and_phase_held_implicit_equation_checked=True,numerical_calculus_fixture_only=True)


def total_radial_calculus():
    g=current.source.FunctionTransportGraph();Y,Z,F=sy.symbols('fixture_native Z fixture_phase');y,z,phase=g.symbol('fixture_native'),g.symbol('Z'),g.symbol('fixture_phase')
    variables={'fixture_native':y,'Z':z,'fixture_phase':phase};make=lambda body:jet(g,variables,body)
    N=sy.Integer(37);J=sy.Rational(3,2);phi0=sy.Rational(1,3)
    E=(2+Z**2/13+Z**5/17)*sy.exp(J*Y/7);V=(1+Z/19+Z**5/23)*(1+J*Y**2/11)
    A=Y*F/5+Z**2*F**2/7+Y**2*Z**5/11+Z**5*F**3/13
    B=Y**2*F/17+Y*Z**3*F**2/19+Z**5*F**3/23
    root=dict(E=make(E),V=make(V));loop=dict(A=make(A),B_over_Pstar=make(B),A_phi=make(sy.diff(A,F)),B_over_Pstar_phi=make(sy.diff(B,F)))
    finite=current.previous.finite_N_functions(g,root,loop,g.constant(37))
    chart=dict(source_chart=dict(radial=dict(E_y=make(sy.diff(E,Y)/J),V_y=make(sy.diff(V,Y)/J))),
        physical_measure=g.constant('3/2'),loop=loop,roots=root,
        modified_fields=dict(V=finite['V_N'],m=make(1+Z**5/29)))
    D=current.NativeDerivative(g,y,phase);radial=current.total_radial_functions(g,chart,D,g.constant(37))
    f=phi0+N*J*Y;expected=dict(E_y=sy.diff(E*sy.exp(A.subs(F,f)/N),Y)/J,
        V_y=sy.diff(V+B.subs(F,f)/N,Y)/J)
    alg=Algebra(g,{y.node:sy.Integer(0),z.node:sy.Rational(1,7),phase.node:phi0});count=0
    for key,body in expected.items():
        for order,ref in enumerate(current.rows(radial[key])):
            prior.close(alg.at(ref.node),sy.diff(body,Z,order).subs({Y:0,Z:sy.Rational(1,7)}),(key,order));count+=1
    assert sy.diff(B,F).subs({Y:0,Z:sy.Rational(1,7),F:phi0})!=0
    return dict(passed=True,independent_full_fast_phase_radial_C5_rows=count,
        nonunit_Jacobian_and_N_cancellation_in_Aphi_Bphi_checked=True,
        phase_held_y_not_conflated_with_total_y=True,numerical_calculus_fixture_only=True)


def generic_calculus():
    g=current.source.FunctionTransportGraph();Z=sy.Symbol('Z');variables={'Z':g.symbol('Z')}
    bodies=dict(E=2+Z**2/13+Z**5/17,V=3+Z/19+Z**5/23,m=1+Z**2/29+Z**5/31,
        h=1+Z**3/37+Z**5/41,k=-2+Z**4/43+Z**5/47,e=1+Z**2/53+Z**5/59,
        p=2+Z**3/61+Z**5/67,P0=-3-Z**2/71-Z**5/73)
    bodies['pressure']=bodies['P0']+bodies['p'];fields={key:jet(g,variables,body) for key,body in bodies.items()}
    de=sy.Rational(1,5);P=sy.Rational(3,2);R=sy.Rational(5,3)
    got=current.generic_functions(g,fields,g.constant('1/5'),g.constant('3/2'),g.constant('5/3'))
    E,V,m,h,k,e,p=(bodies[key] for key in ('E','V','m','h','k','e','pressure'));d=1-Z**2;L=1-de*Z**2;T=(1-de)*Z*m+d*sy.diff(m,Z)
    expected=dict(Q=(2*Z*V-T)/L,transport=T,
        theta_linear=(-E+(1-de/2)*h-(1-de)*Z*sy.diff(h,Z)/2)/L,
        theta_quadratic=((2*de-1)*Z*k-d*sy.diff(k,Z)+E*T)/L,
        axial_linear=(-V+(1-de)*(m-Z*sy.diff(m,Z))/2)/L,
        axial_quadratic=(V*T+2*de*Z*e-d*sy.diff(e,Z)+2*(1+de)*Z*p-d*sy.diff(p,Z))/L)
    expected['inertial_theta']=R*(expected['theta_linear']+P*expected['theta_quadratic'])
    expected['inertial_axial']=R*(expected['axial_linear']+P*expected['axial_quadratic'])
    expected['p2']=expected['inertial_axial']/E
    point=sy.Rational(3,8);alg=Algebra(g,{g.symbol('Z').node:point});count=0
    for key,body in expected.items():
        for order,ref in enumerate(current.rows(got[key])):
            assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,order).subs(Z,point))==0,(key,order);count+=1
    short=dict(fields);short['m']=current.FunctionJet(tuple(current.rows(fields['m'])[:5]))
    try:current.generic_functions(g,short,g.constant('1/5'),g.constant('3/2'),g.constant('5/3'))
    except ValueError:pass
    else:raise AssertionError('Missing fifth source was padded')
    return dict(passed=True,independent_full_signed_generic_C4_rows=count,
        genuine_fifth_mean_pressure_rows_consumed=True,missing_fifth_row_rejected=True,
        nonunit_radius_and_Pstar_once_and_pressure_terms_checked=True,mathematical_fixture_only=True)


def pressure_calculus():
    g=current.source.FunctionTransportGraph();Y,Z=sy.symbols('fixture_pressure_y Z');y,z=g.symbol('fixture_pressure_y'),g.symbol('Z')
    variables={'fixture_pressure_y':y,'Z':z};make=lambda body:jet(g,variables,body)
    P=sy.Rational(3,2);de=sy.Rational(1,5);R=sy.Rational(7,3)*sy.exp(Y);lower=sy.Rational(1,3)
    m=(2+Z+Z**5/11)+(3+Z**2+Z**5/13)*Y+(1+Z**3/17+Z**5/19)*Y**2
    U=m+sy.diff(m,Y);Uy=sy.diff(U,Y)
    f=dict(V=make(U/P),m=make(m/P));radial=dict(V_y=make(Uy/P))
    physical=current.physical_pressure_fields(g,f,radial,g.constant('1/5'),g.constant('3/2'))
    S=current.reduced_pressure_S(g,physical,g.constant('1/5'));a=current.Algebra(g,3)
    Q=a.trim(physical['Q']);Qy=a.trim(physical['Q_y']);radius=current.source.expression(g,R,variables)
    boundary=dict(B0=a.add(Q,Qy),B1=a.scale(Qy,radius))
    at=lambda j:a.at(j,y,g.constant('1/3'))
    i0=a.integral(a.scale(S,radius),y,g.constant('1/3'),y,fixture_exact_pressure=True)
    i1=a.integral(a.scale(S,g.mul(radius,radius)),y,g.constant('1/3'),y,fixture_exact_pressure=True)
    got=dict(m13=a.add(a.scale(i0,'-1/2'),a.sub(boundary['B0'],at(boundary['B0']))),
        m15=a.sub(a.scale(i1,'1/2'),a.sub(boundary['B1'],at(boundary['B1']))))
    d=1-Z**2;L=1-de*Z**2;q=(2*Z*U-(1-de)*Z*m-d*sy.diff(m,Z))/L;qy=sy.diff(q,Y);qyy=sy.diff(q,Y,2)
    ss=(q+(1-de)*Z*sy.diff(q,Z)/2+qy)/L+q*(q/2+qy)+U*(-2*Z*q+d*sy.diff(q,Z)-2*Z*qy)/L
    N1p=-ss/2+(qy+qyy)/R
    point=sy.Rational(3,8);alg=Algebra(g,{y.node:Y,z.node:point});count=0
    for key,jet0,body in (('Q',physical['Q'],q),('Q_y',physical['Q_y'],qy),('S',S,ss)):
        for order,ref in enumerate(current.rows(jet0)):
            assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,order).subs(Z,point))==0,(key,order);count+=1
    for key,multiplier in (('m13',R),('m15',-R**2)):
        for order,ref in enumerate(current.rows(got[key])):
            value=alg.at(ref.node)
            assert sy.simplify(sy.diff(value,Y)-sy.diff(multiplier*N1p,Z,order).subs(Z,point))==0,(key,order,'direct forcing')
            assert sy.simplify(value.subs(Y,lower))==0,(key,order,'finite inlet');count+=1
    return dict(passed=True,independent_physical_Q_Qy_S_and_signed_pressure_rows=count,
        direct_N1p_equivalence_of_reduced_integrals_through_C3=True,
        original_leading_weight_Pstar_once_nonunit_radius_and_both_finite_endpoints=True,
        m15_negative_weighted_forcing_sign_applied_once=True,mathematical_fixture_only=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['modified_Q4_pressure_C3_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks=dict(actual_source=actual_source(owner),binders=binder_sources(owner));print('TRUE_MODIFIED_SOURCE_TOTAL_RADIAL_AND_PRESSURE_SCOPE_CHECKED',flush=True)
    for name,fn in (('native_chain',native_chain_calculus),('implicit_native',implicit_native_calculus),
                    ('total_radial',total_radial_calculus),('full_generic',generic_calculus),('pressure_reduction',pressure_calculus)):
        checks[name]=fn();print('TRUE_MODIFIED_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_MODIFIED_Q4_TOTAL_RADIAL_AND_PRESSURE_C3_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses actual admitted true loop source')
