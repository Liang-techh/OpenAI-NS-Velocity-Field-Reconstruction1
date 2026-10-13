"""Direct physical N1p integrals, true C1 rows and selected join identities."""
import json
from pathlib import Path
import time
import scipy.integrate
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions_check as prior_check


def fixture_symbolic(g,ref,cache):
    """Exact arithmetic on the admitted fixture domain L>=4/5.

    Rational cancellation before integration avoids artificial branches
    at Z=+/-sqrt(5), where the pressure formula is undefined and which lie
    outside the entire admitted axial interval [-1,1].
    """
    if ref.node in cache:return cache[ref.node]
    row=g.nodes[ref.node];op=row['operation']
    at=lambda i:fixture_symbolic(g,current.source.FunctionRef(g,i),cache)
    if op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
    elif op=='bound_variable':value=sy.Symbol(row['name'])
    elif op=='sum':value=sy.Add(*(at(i) for i in row['arguments']))
    elif op=='product':value=sy.Mul(*(at(i) for i in row['arguments']))
    elif op=='negative':value=-at(row['argument'])
    elif op=='positive_quotient':value=at(row['numerator'])/at(row['denominator'])
    elif op=='analytic_unary':value=getattr(sy,row['name'])(at(row['argument']))
    elif op=='function_substitution':value=at(row['expression']).subs(at(row['variable']),at(row['value']))
    elif op=='definite_integral':
        assert row.get('original_complete_pulse_chart') in ('main','gap','end')
        body=sy.cancel(at(row['integrand']))
        value=sy.integrate(body,(sy.Symbol(row['variable']),at(row['lower']),at(row['upper'])),conds='none')
        assert not value.has(sy.Piecewise)
    else:raise ValueError('Fixture cannot substitute a source/cap node: '+op)
    cache[ref.node]=value
    return value


def exact_selected_rows():
    A0,A1,D0,D1,I1,I2,P1,P2,mu,ap=sy.symbols('A0 A1 D0 D1 I1 I2 P1 P2 mu ap')
    det=A1*D0-A0*D1
    inverse=lambda r1,r2:sy.Matrix([(A1*r2-D1*r1)/det,(D0*r1-A0*r2)/det])
    u=inverse(-mu*I1,-(I2-I1));v=inverse(-mu*P1,-(P2-P1));C=u+ap*v
    assert sy.cancel((sy.Matrix([[A0,A1]])*C)[0]+mu*(I1+ap*P1))==0
    assert sy.cancel((sy.Matrix([[A0+mu*D0,A1+mu*D1]])*C)[0]+mu*(I2+ap*P2))==0
    # The explicit absolute determinant equals A1 D0-A0 D1 in the exact
    # original normalized beta measure; the finite second row is retained.
    J1,J2=sy.symbols('J1 J2');lam1=sy.Rational(1,2)-mu;lam2=sy.Rational(1,2)-2*mu
    first=[sy.exp(center*lam1)*J1 for center in (-3,-1)]
    second=[sy.exp(center*lam2)*J2 for center in (-3,-1)]
    divided=[(b-a)/mu for a,b in zip(first,second)]
    determinant=first[1]*divided[0]-first[0]*divided[1]
    explicit=sy.exp(-2+5*mu)*J1*J2*(sy.exp(2*mu)-1)/mu
    assert sy.simplify(determinant-explicit)==0
    return dict(passed=True,pointwise_selected_C_u_plus_ap_v_row_identities_exact=True,
        finite_mu_second_row_and_explicit_absolute_determinant_proved=True,
        original_two_moment_joins_are_function_identities=True)


def direct_physical_fixture():
    """Independent original forcing, without using reduced-S to define it."""
    g=current.source.FunctionTransportGraph();Z,R=sy.symbols('Z fixture_physical_R')
    z=g.symbol('Z');c=g.constant;delta=sy.Rational(1,5);mu=sy.Rational(3,20)
    unit=sy.Rational(3,2);lo,first,second,hi=map(sy.Rational,('37/10','43/10','61/10','74/10'))
    m=unit*(2+Z+Z**2+(1-2*Z+Z**3)*R);U=m+R*sy.diff(m,R)
    d,L=1-Z**2,1-delta*Z**2
    Q=sy.cancel((2*Z*U-(1-delta)*Z*m-d*sy.diff(m,Z))/L);Qy=R*sy.diff(Q,R)
    Uy=R*sy.diff(U,R)
    fields={};edges=dict(main=(lo,first),gap=(first,second),end=(second,hi))
    for chart,(left,right) in edges.items():
        x=g.symbol('pressure_fixture_'+chart)
        if chart=='end':
            radius=g.mul(c(str(hi)),g.unary('exp',x))
            endpoints=(g.unary('log',c(str(left/hi))),g.zero)
        else:
            radius=g.mul(c(str(lo)),g.unary('exp',g.quotient(x,c(str(mu)),'positive fixture mu')))
            endpoints=tuple(g.mul(c(str(mu)),g.unary('log',c(str(edge/lo)))) for edge in (left,right))
        variables={'Z':z,str(R):radius}
        bind=lambda expr:current.source.expression(g,sy.sympify(expr),variables)
        jet=lambda expr:current.own.c3.C3Function(*(bind(sy.diff(expr,Z,n)) for n in range(4)))
        c2=lambda expr:current.own.c2.C2Function(*(bind(sy.diff(expr,Z,n)) for n in range(3)))
        fields[chart]=dict(query=x,radius=radius,U=jet(U),U_y=jet(Uy),Q=c2(Q),Q_y=c2(Qy),
            domain=dict(lower=endpoints[0],upper=endpoints[1]))
    functions,pieces,bindings=current.pulse_pressure_functions(g,c(str(delta)),c(str(mu)),fields)
    # Build the original radial differential forcing independently.
    T=(Q+(1-delta)*Z*sy.diff(Q,Z)/2+R*sy.diff(Q,R))/L
    convective=Q*(Q/2+R*sy.diff(Q,R))
    axial=U*(-2*Z*Q+d*sy.diff(Q,Z)-2*Z*R*sy.diff(Q,R))/L
    N1p=-(T+convective+axial)/2+2*sy.diff(Q,R)+R*sy.diff(Q,R,2)
    direct=lambda left,right:dict(
        m13=sy.integrate(sy.Poly(sy.cancel(N1p),R).as_expr(),(R,left,right)),
        m15=-sy.integrate(sy.Poly(sy.cancel(R*N1p),R).as_expr(),(R,left,right)))
    cache={}
    symbolic=lambda ref:fixture_symbolic(g,ref,cache)
    full=direct(lo,hi)
    for chart,functions_at_chart in pieces.items():
        want=direct(*edges[chart])
        for key,pair in functions_at_chart.items():
            assert sy.cancel(symbolic(pair.value)-want[key])==0,(chart,key,'value')
            assert sy.cancel(symbolic(pair.Z)-sy.diff(want[key],Z))==0,(chart,key,'Z')
    for key,pair in functions.items():
        assert sy.cancel(symbolic(pair.value)-full[key])==0,(key,'complete value')
        assert sy.cancel(symbolic(pair.Z)-sy.diff(full[key],Z))==0,(key,'complete Z')
        assert sy.cancel(symbolic(pair.value)-sum(symbolic(pieces[chart][key].value) for chart in pieces))==0
        assert sy.cancel(symbolic(pair.Z)-sum(symbolic(pieces[chart][key].Z) for chart in pieces))==0
        assert sy.cancel(full[key])!=0
    expected_products=dict(I0_q_Z=sy.diff(Q,Z),I0_U_q_Z=U*sy.diff(Q,Z),
        I0_U_squared=U**2,I0_U_U_Z=U*sy.diff(U,Z),I1_q=Q,I1_q_Z=sy.diff(Q,Z),
        I1_q_squared=Q**2,I1_U_q_Z=U*sy.diff(Q,Z),I1_U_squared=U**2,I1_U_U_Z=U*sy.diff(U,Z))
    for chart,values in bindings['products'].items():
        x=fields[chart]['query'];radius=symbolic(fields[chart]['radius'])
        for key,pair in values.items():
            want=expected_products[key].subs(R,radius)
            assert sy.cancel(symbolic(pair.value)-want)==0
            assert sy.cancel(symbolic(pair.Z)-sy.diff(want,Z))==0
            integral=bindings['integrals'][chart][key]
            for n,ref in enumerate((integral.value,integral.Z)):
                node=g.nodes[ref.node]
                assert node['variable']==g.nodes[x.node]['name'] and node['ordinary_slow_Z_derivative_order']==n
                assert node['endpoint_independent_of_Z']
    # Neither finite viscous boundary vanishes in this nonunit fixture.
    for end in bindings['global_finite_endpoints']:
        assert sy.cancel(symbolic(end['q_plus_q_y'].value))!=0
        assert sy.cancel(symbolic(end['R_q_y'].value))!=0
    return dict(passed=True,independent_mathematical_fixture=True,original_source_evaluation=False,
        admitted_axial_domain='[-1,1]; 1-Z^2/5>=4/5; poles +/-sqrt(5) excluded',
        nonunit_Pstar=str(unit),finite_nonzero_mu=str(mu),delta=str(delta),physical_R_endpoints=[str(lo),str(hi)],
        all_three_chart_jacobians_checked_symbolically=True,
        all_ten_actual_integrands_and_true_first_Z_rows_checked=True,
        direct_original_N1p_integration_equals_complete_m13_m15_functions=True,
        direct_first_Z_derivatives_checked_as_symbolic_functions=True,
        nonzero_signed_finite_viscous_endpoints_retained=True,
        internal_endpoints_telescope_as_functions=True,m15_minus_sign_applied_once=True)


def actual_source_checks(owner):
    g=owner.g;zero=current.source.C1Function(g.zero,g.zero)
    complete=owner.own.jets.complete
    selected,_=complete.Gamma.pulse.original_coefficient_functions(complete.functions['weighted'])
    for key,pair in selected.items():
        assert pair==complete.pulse_functions[key],key
        assert pair.value==owner.own.jets.functions['pulse'][key].value
        assert pair.Z==owner.own.jets.functions['pulse'][key].Z
    low,high=owner.bindings['global_finite_endpoints']
    assert low['q_plus_q_y']==zero
    assert low['R_q'].value!=g.zero and low['R_q_y'].value!=g.zero
    assert all(pair==zero for pair in high.values())
    cert=owner.endpoint_certificate
    assert cert['exact_second_row']=='A+mu*D'
    assert cert['terminal_beta_centers']==[-3,-1] and cert['terminal_s']==0
    for chart,f in owner.own.fields.items():
        expected=current.physical.nonradial_integrands(g,f['U'],f['Q'])
        for key,pair in owner.bindings['products'][chart].items():
            assert pair.value==expected[key].value and pair.Z==expected[key].Z
            for n,ref in enumerate((owner.bindings['integrals'][chart][key].value,
                                    owner.bindings['integrals'][chart][key].Z)):
                node=g.nodes[ref.node]
                assert node['original_complete_pulse_chart']==chart
                assert node['ordinary_slow_Z_derivative_order']==n
                assert node['lower']==f['domain']['lower'].node and node['upper']==f['domain']['upper'].node
                assert node['no_finite_history_or_chart_cutoff'] and node['physical_Pstar_not_reapplied']
    gap=owner.bindings['products']['gap']
    assert all(gap[key].value==g.zero for key in ('I0_U_q_Z','I0_U_squared','I0_U_U_Z'))
    assert all(gap[key].value!=g.zero for key in ('I0_q_Z','I1_q','I1_q_squared'))
    return dict(passed=True,actual_selected_C1_handles_replayed_from_original_selector=True,
        actual_selected_C3_value_and_first_Z_handles_preserved=True,
        source_native_nonzero_Rp_Q_and_R_Qy_retained=True,
        source_flat_startup_Q_plus_Qy_zero_and_empty_end_zero_proved=True,
        ten_products_per_chart_rebuilt_from_actual_own_U_Q=True,
        both_whole_chart_integral_derivative_rows_and_exact_domains_preserved=True,
        gap_nonzero_pressure_history_terms_retained=True)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1PulsePressureFunctions(require_checked=False)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['exact_graph_nodes']==owner.g.nodes and raw['pulse_pressure_bindings']==owner.snapshot
    assert not any(raw[key] for key in current.OPEN)
    checks=dict(direct_physical=direct_physical_fixture(),selected_rows=exact_selected_rows(),
        actual_source=actual_source_checks(owner))
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        signed_function_graph_only=True,uniform_original_pressure_numerical_accuracy_claimed=False,
        full_five_defect_functions_and_actual_repair_not_installed=True,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name),
            Path(prior_check.__file__).name:current.sha(Path(prior_check.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_COMPLETE_PULSE_PRESSURE_C1_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()
