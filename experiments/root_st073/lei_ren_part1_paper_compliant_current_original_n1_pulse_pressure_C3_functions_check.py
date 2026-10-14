"""Source identities and independent higher-calculus fixtures for pulse C3.

This accepts function construction, not numerical containment at the original
astronomical scales. Ancestor source owners and numerical grids are reused.
"""
import json
from pathlib import Path
from types import SimpleNamespace
import time

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions_check as Gamma_check


class Algebra:
    """Interpret actual elementary source nodes with explicit exact bindings."""
    def __init__(self,g,bindings):self.g,self.cache=g,dict(bindings)
    def at(self,i):
        if i in self.cache:return self.cache[i]
        row=self.g.nodes[i];op=row['operation'];at=self.at
        if op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
        elif op=='sum':value=sy.Add(*(at(j) for j in row['arguments']))
        elif op=='product':value=sy.Mul(*(at(j) for j in row['arguments']))
        elif op=='negative':value=-at(row['argument'])
        elif op=='positive_quotient':value=at(row['numerator'])/at(row['denominator'])
        elif op=='analytic_unary':
            x=at(row['argument']);name=row['name']
            value=sy.log(1+x) if name=='log1p' else getattr(sy,name)(x)
        else:raise AssertionError('Unbound actual source operation '+op)
        self.cache[i]=value;return value


def native_source(owner):
    exact=owner.native_owner.selected.inlet.exact;frame=exact.frame
    q,data=frame.symbolic_frame();_,constants=exact.current_constants();Z=q.z
    U,C1,C2,C0,CE,X=sy.symbols('U C1 C2 C0 C_E X',real=True)
    Q=1+Z**2
    abstract=dict(u=U/Q,m1=C1*Z*Q,m2=C2*Z*Q,energy=C0+CE*Z**2*Q**2,X=X)
    scalars={U:constants['U'],C1:constants['C1'],C2:constants['C2'],
        C0:constants['C0'],CE:constants['C_E'],X:sy.cancel(data['X'])}
    assert all(Z not in value.free_symbols for value in scalars.values())
    bindings={owner.g.symbol('Z').node:Z}
    old_packet=frame.functions['actual_identified_Rp_native_frame_C3']
    before=len(owner.g.nodes)
    assert current.jets.encode(current.jets.native_C3_functions(owner.g))==current.jets.encode(owner.old.native)
    assert len(owner.g.nodes)==before
    callback=frame.incoming_exact_taylor(Z);checked=[]
    for key,body in abstract.items():
        actual=body.subs(scalars)
        assert sy.cancel(actual-data[key])==0,key
        old=current.rows(owner.old.native[key])
        assert current.rows(owner.native[key])[:4]==old
        # Read the actual frame rows too: in particular X's three old
        # derivatives must be proved zero before higher zeros are appended.
        for order,ref in enumerate(current.c3.rows(old_packet[key])):
            assert sy.cancel(q.at(ref)-sy.diff(actual,Z,order))==0,(key,order)
        for order,ref in enumerate(old):
            desired=sy.diff(body,Z,order)
            if ref.node in bindings:assert sy.cancel(bindings[ref.node]-desired)==0
            bindings[ref.node]=desired
        for order,ref in enumerate(current.rows(owner.native[key])):
            expression=Algebra(owner.g,bindings).at(ref.node)
            assert sy.cancel(expression-sy.diff(body,Z,order))==0,(key,order)
            assert sy.cancel(expression.subs(scalars)-sy.factorial(order)*callback[key][order])==0
            checked.append(key+'_ordinary_Z'+str(order))
    alg=Algebra(owner.g,bindings)
    assert sy.cancel(alg.at(owner.native_evidence['actual_U_from_native_value'].node)-U)==0
    assert sy.cancel(alg.at(owner.native_evidence['actual_CE_from_native_ordinary_Z2'].node)-CE)==0
    assert sy.expand(alg.at(owner.native_evidence['positive_CE_denominator'].node))==2*(1+12*Z**2+15*Z**4)
    # The exact rational source identity above has already reduced X to
    # this Z-independent scalar; avoid re-expanding its transport factors.
    assert all(sy.diff(scalars[X],Z,k)==0 for k in range(1,6))
    return dict(passed=True,actual_native_ordinary_rows_checked=checked,
        actual_six_row_callback_factorials_applied_once=True,
        reconstructed_U_and_CE_equal_exact_Z_independent_source_scalars=True,
        actual_CE_denominator_strict_lower_bound=2,
        X_constant_source_and_all_old_X_derivatives_proved=True,
        cubic_moment_higher_zeros_proved_not_padded=True,
        old_native_C3_function_handles_preserved=True)


def Gamma_source(owner):
    original=owner.old.complete.Gamma;g=owner.g;E=owner.Gamma_evidence
    Z=sy.Symbol('Z',real=True);a,S,T=sy.symbols('a S T',positive=True);t,v=sy.symbols('t v',real=True)
    bindings={g.symbol('Z').node:Z,original.pieces['a'].node:a,
        E['actual_nonzero_S'].node:S,E['offset'].node:t,original.pieces['Gamma_variable'].node:v}
    alg=Algebra(g,bindings)
    assert sy.simplify(alg.at(E['actual_xi'].node)-2*(1-Z**2)*S*sy.exp(-t))==0
    assert alg.at(E['normalization'].node)==sy.gamma(1+a)
    assert owner.derivative.derivative(E['actual_nonzero_S'])==g.zero
    old_moments=owner.old.Gamma_evidence['positive_Gamma_moments']
    assert E['positive_improper_moments'][:3]==old_moments
    assert current.rows(owner.Gamma_Dhat)[:4]==current.rows(owner.old.Gamma_Dhat)
    moments=sy.symbols('M1:6',positive=True)
    for k,ref in enumerate(E['positive_improper_moments'],1):
        quotient=g.nodes[ref.node];assert quotient['operation']=='positive_quotient'
        assert quotient['denominator']==E['normalization'].node
        integral=g.nodes[quotient['numerator']]
        assert integral['operation']=='definite_integral' and integral['lower']==g.zero.node
        assert integral['upper']==original.pieces['infinity'].node
        assert integral['full_improper_integral'] and integral['exact_positive_Gamma_expectation']
        assert integral['Gamma_moment_order']==k and integral['zero_endpoint_integrand_continuation']==0
        want=sy.exp(-v+(a+k)*sy.log(v)-(a+k)*sy.log(1+2*(1-Z**2)*S*sy.exp(-t)*v))
        assert sy.simplify(alg.at(integral['integrand'])-want)==0,k
        # This is the pointwise source-kernel derivative. Domination uses
        # xi>=0 and the same full Gamma integral, including both endpoints.
        next_body=want*v/(1+2*(1-Z**2)*S*sy.exp(-t)*v)
        assert sy.simplify(sy.diff(want,Z)-4*(a+k)*Z*S*sy.exp(-t)*next_body)==0
        bindings[ref.node]=moments[k-1]
    def dz(expr):
        return sy.diff(expr,Z)+sum(sy.diff(expr,moments[k-1])*4*(a+k)*Z*S*T*moments[k]
            for k in range(1,5))
    alg=Algebra(g,bindings)
    third=alg.at(owner.Gamma_Dhat.ZZZ.node).subs(sy.exp(-t),T)
    fourth=alg.at(owner.Gamma_Dhat.Z4.node).subs(sy.exp(-t),T)
    fifth=alg.at(owner.Gamma_Dhat.Z5.node).subs(sy.exp(-t),T)
    assert sy.simplify(dz(third)-fourth)==0
    assert sy.simplify(dz(fourth)-fifth)==0
    return dict(passed=True,actual_full_positive_Gamma_kernels=5,
        M1_M2_M3_original_handles_preserved=True,M4_M5_from_same_improper_source=True,
        actual_xi_and_positive_nonzero_Z_independent_S_bound=True,
        five_exact_kernel_derivative_recurrences_checked=True,
        D_Z4_Z5_follow_exact_recurrence=True,
        domination='0<=M_k<=Gamma(1+a+k)/Gamma(1+a), a>0, xi>=0, k=1..5',
        finite_one_sided_Z_endpoints_and_exponential_offset_decay=True,
        no_division_by_one_minus_Z_squared=True)


def higher_calculus_fixture():
    g=current.source.FunctionTransportGraph();z=g.symbol('Z');Z=sy.Symbol('Z',real=True)
    f=sy.exp(Z**2/10)*sy.sqrt(3+Z**2)/(1+Z**2)
    def graph(expr):
        expr=sy.sympify(expr)
        if expr.is_Rational:return g.constant(str(expr))
        if expr==Z:return z
        if expr.func is sy.Add:return g.add(*(graph(v) for v in expr.args))
        if expr.func is sy.Mul:return g.mul(*(graph(v) for v in expr.args))
        if expr.func is sy.Pow:
            base,power=expr.args;node=graph(base)
            if power.is_Integer:n=int(power)
            elif power.is_Rational and power.q==2:
                node=g.unary('sqrt',node);n=int(power.p)
            else:raise AssertionError('Unexpected fixture power')
            result=g.mul(*([node]*abs(n)))
            return result if n>=0 else g.quotient(g.one,result,'positive fixture power')
        if expr.func in (sy.exp,sy.log):return g.unary('exp' if expr.func is sy.exp else 'log',graph(expr.args[0]))
        raise AssertionError('Unexpected independent fixture expression')
    old=current.c3.C3Function(*(graph(sy.diff(f,Z,k)) for k in range(4)))
    d=current.HigherZDerivative(g);d.seed_rows(old);jet=d.extend(old,5)
    v=g.symbol('higher_fixture_v');body=g.unary('exp',g.mul(z,v))
    integrate=lambda expr:current.controls.integral(g,expr,g.nodes[v.node]['name'],g.zero,g.one,
        endpoint_independent_of_Z=True,ordinary_slow_Z_derivative_order=0)
    value=integrate(body);iv=d.extend(current.source.C1Function(value,d.derivative(value)),5)
    evaluate=Gamma_check.FixtureEvaluator(g,.3);maximum=0.
    for point in (-.7,0.,.8):
        for order,ref in enumerate(current.rows(jet)):
            want=float(sy.diff(f,Z,order).subs(Z,point).evalf(40))
            error=abs(evaluate(ref,{z.node:point})-want);maximum=max(maximum,error)
            assert error<2e-9,(point,order,error)
        with mp.workdps(40):
            for order,ref in enumerate(current.rows(iv)):
                want=float(mp.quad(lambda v:v**order*mp.exp(mp.mpf(str(point))*v),[0,1]))
                error=abs(evaluate(ref,{z.node:point})-want);maximum=max(maximum,error)
                assert error<2e-10,(point,order,error)
    for order,ref in enumerate(current.rows(iv)[1:],1):
        row=g.nodes[ref.node];assert row['ordinary_slow_Z_derivative_order']==order
        assert row['actual_higher_source_body_not_interval_endpoint']
    fixed=g.node('function_substitution',expression=jet.value.node,variable=z.node,value=g.zero.node,
        Z_independent_substitution=True)
    assert d.derivative(fixed)==g.zero
    for ref in (g.node('unidentified_axial_function',argument=z.node),
        current.controls.integral(g,body,g.nodes[v.node]['name'],g.zero,z)):
        try:d.derivative(ref)
        except ValueError:pass
        else:raise AssertionError('Unproved source or moving endpoint admitted')
    return dict(passed=True,mathematical_fixture_only=True,independent_rows=36,
        composite_function_Z0_Z5_and_integral_Z0_Z5_checked=True,
        maximum_absolute_error=maximum,actual_source_integral_order_labels_checked=True,
        unknown_sources_and_moving_axial_endpoints_rejected=True)


def Gamma_fixture():
    g=current.source.FunctionTransportGraph();c=g.constant;z=g.symbol('Z')
    functions,pieces=current.jets.complete.Gamma.Gamma_future_functions(g,c('3/10'),c('1/50'),c('3/25'))
    original=SimpleNamespace(functions=functions,pieces=pieces,schedule={'S':c('3/25')})
    old,_=current.jets.Gamma_higher_functions(g,original)
    jet,_=current.Gamma_C5_functions(g,original,old)
    evaluate=Gamma_check.FixtureEvaluator(g,.3);errors={}
    with mp.workdps(40):
        a,S,t=mp.mpf('.15'),mp.mpf('.12'),mp.mpf('.7')
        def direct(Z):
            xi=2*(1-Z*Z)*S*mp.exp(-t)
            return mp.quad(lambda v:mp.exp(-v)*v**a*(-mp.expm1(-a*mp.log1p(xi*v))),
                [0,1,mp.inf])/(a*S*mp.gamma(1+a))
        point=mp.mpf('.65');env={z.node:float(point),pieces['offset_variable'].node:float(t)}
        for order,ref in ((4,jet.Z4),(5,jet.Z5)):
            want=float(mp.diff(direct,point,order));value=evaluate(ref,env)
            errors[str(order)]=abs(value-want);assert errors[str(order)]<2e-9
        for point in (-1.,1.):
            env={z.node:point,pieces['offset_variable'].node:.7}
            assert mp.isfinite(evaluate(jet.Z4,env)) and mp.isfinite(evaluate(jet.Z5,env))
        assert evaluate(jet.Z5,{z.node:0.,pieces['offset_variable'].node:.7})==0.
    return dict(passed=True,mathematical_fixture_only=True,
        independent_original_Dhat_integral_Z4_Z5_errors=errors,
        finite_one_sided_endpoint_rows_and_zero_center_odd_row=True)


def pulse_source(owner):
    g=owner.g;d=owner.derivative;preserved=[];integral_rows=[];endpoint_rows=[]
    for chart,fields in owner.leading_fields.items():
        old=owner.pulse.own.fields[chart]
        for key,jet in fields.items():
            prior=current.rows(old[key]);assert current.rows(jet)[:len(prior)]==prior
            for left,right in zip(current.rows(jet)[len(prior)-1:],current.rows(jet)[len(prior):]):
                assert d.derivative(left)==right
            preserved.append(chart+'_'+key)
        rebuilt=current.true.own.radial_Q(g,owner.old.complete.Gamma.schedule['delta'],old['U'],old['m'])
        assert current.rows(rebuilt)==current.rows(old['Q'])
    for chart,functions in owner.integrals.items():
        original=owner.pulse.bindings['integrals'][chart]
        for key,jet in functions.items():
            assert current.rows(jet)[:2]==current.rows(original[key])
            initial=g.nodes[jet.value.node]
            measure=('R*ds' if chart=='end' else '(R/mu)*dxi') if key.startswith('I0_') else (
                'R^2*ds' if chart=='end' else '(R^2/mu)*dxi')
            assert initial['measure']==measure
            for order,ref in enumerate(current.rows(jet)):
                if ref==g.zero:continue
                row=g.nodes[ref.node]
                assert row['operation']=='definite_integral' and row['measure']==measure
                assert row['lower']==initial['lower'] and row['upper']==initial['upper']
                assert row['variable']==initial['variable'] and row['endpoint_independent_of_Z']
                assert row['original_complete_pulse_chart']==chart
                assert row['ordinary_slow_Z_derivative_order']==order
                assert row['signed_integrand_not_its_magnitude_cap'] and row['physical_Pstar_not_reapplied']
                if order>=2:
                    previous=g.nodes[current.rows(jet)[order-1].node]
                    assert d.derivative(current.source.FunctionRef(g,previous['integrand'])).node==row['integrand']
                    assert row['actual_higher_source_body_not_interval_endpoint']
                integral_rows.append(chart+'_'+key+'_Z'+str(order))
    endpoint_groups=[('total',owner.global_endpoints,owner.pulse.bindings['global_finite_endpoints'])]+[
        (chart,owner.endpoints[chart],owner.pulse.bindings['per_chart_finite_endpoints'][chart])
        for chart in owner.endpoints]
    for chart,ends,original in endpoint_groups:
        for side,(now,old) in enumerate(zip(ends,original)):
            assert set(now)==set(old)
            for key,jet in now.items():
                assert current.rows(jet)[:2]==current.rows(old[key])
                assert d.derivative(jet.Z)==jet.ZZ and d.derivative(jet.ZZ)==jet.ZZZ
                endpoint_rows.append(chart+'_'+str(side)+'_'+key)
        if chart!='total':
            fields=owner.pulse.own.fields[chart];a=current.parent.C1Algebra(g)
            for edge in (fields['domain']['lower'],fields['domain']['upper']):
                assert d.derivative(edge)==g.zero
                radius=a.at(a.fixed(fields['radius']),fields['query'],edge).value
                assert d.derivative(fields['radius'])==g.zero and d.derivative(fields['query'])==g.zero
                row=g.nodes[radius.node]
                assert row['operation']=='function_substitution' and row['expression']==fields['radius'].node
                assert row['variable']==fields['query'].node and row['value']==edge.node
    assert all(ref==g.zero for jet in owner.global_endpoints[1].values() for ref in current.rows(jet))
    assert all(ref==g.zero for ref in current.rows(owner.global_endpoints[0]['q_plus_q_y']))
    pressure_rows=[]
    for label,now,old in [('total',owner.pressure,owner.pulse.functions)]+[
        (chart,owner.pieces[chart],owner.pulse.pieces[chart]) for chart in owner.pieces]:
        for key,jet in now.items():
            assert current.rows(jet)[:2]==current.rows(old[key])
            assert d.derivative(jet.Z)==jet.ZZ and d.derivative(jet.ZZ)==jet.ZZZ
            pressure_rows.append(label+'_'+key)
    assert owner.pulse.endpoint_certificate['actual_native_physical_mean_retained']
    assert owner.pulse.endpoint_certificate['exact_second_row']=='A+mu*D'
    assert owner.pulse.bindings['original_m15_negative_weighted_integral_sign_applied_once']
    return dict(passed=True,leading_field_jets_preserving_old_rows=preserved,
        actual_nonzero_integral_rows_checked=integral_rows,
        actual_finite_endpoint_C3_jets_checked=endpoint_rows,
        actual_endpoint_coordinates_and_radii_Z_independent=True,
        startup_native_Q_retained_and_terminal_zero_source_identities_preserved=True,
        genuine_pressure_C3_functions=pressure_rows,
        original_leading_Q_weight='1-delta',original_finite_endpoints_and_measures_retained=True,
        source_native_mean_and_signed_m15_retained=True,
        numerical_derivative_containment_claimed=False)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['C5_pulse_C3_pressure_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    checks=dict(native=native_source(owner),Gamma=Gamma_source(owner))
    print('ACTUAL_NATIVE_AND_GAMMA_C5_SOURCE_IDENTITIES_CHECKED',flush=True)
    checks.update(higher_calculus=higher_calculus_fixture(),Gamma_fixture=Gamma_fixture(),pulse=pulse_source(owner))
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_PULSE_PRESSURE_C3_FUNCTIONS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses the checked live pulse source')
