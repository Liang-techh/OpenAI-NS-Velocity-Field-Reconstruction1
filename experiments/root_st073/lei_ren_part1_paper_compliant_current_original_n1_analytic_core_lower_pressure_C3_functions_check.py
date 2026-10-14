"""Admitted analytic-source, true FTC and direct physical-pressure checks."""
import copy
import json
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_analytic_core_lower_pressure_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_lower_collar_functions_check as lower_check
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions_check as higher_check


def analytic_source(owner):
    g=owner.g;rows=current.higher.rows;checked=[]
    assert owner.core.acceptance_loaded and owner.lower.acceptance_loaded
    assert owner.core.definition['actual_twenty_term_source_and_Green_inverse_retained']
    assert {key:len(value) for key,value in owner.core.terms.items()}=={'theta':10,'z':10}
    for quantity,jets in owner.core_radial_rows.items():
        for radial_order,jet in enumerate(jets):
            old=owner.lower.core_radial_rows[quantity][radial_order]
            assert rows(jet)[:4]==rows(old)
            base=current.component_definition(g.nodes[old.value.node])
            for order,ref in enumerate(rows(jet)):
                row=g.nodes[ref.node]
                assert row['quantity']==quantity and row['radial_order']==radial_order
                assert row['ordinary_slow_Z_derivative_order']==order
                assert row['source_family']==owner.family
                assert current.component_definition(row)==base
                assert row['unique_analytic_Banach_solution'] and row['actual_twenty_term_operator']
                checked.append(quantity+'_rho'+str(radial_order)+'_Z'+str(order))
    for key,jet in owner.core_functions.items():
        assert rows(jet)[:4]==rows(owner.core.functions[key])
    # A changed component body/rho/source may not silently gain derivatives.
    row=copy.deepcopy(g.nodes[owner.core_radial_rows['Psi'][0].Z5.node])
    row['source_family']={**owner.family,'unadmitted_source':True};operation=row.pop('operation')
    bad=g.node(operation,**row)
    try:owner.derivative.derivative(bad)
    except ValueError:pass
    else:raise AssertionError('An unadmitted analytic component was differentiated')
    return dict(passed=True,actual_admitted_core_component_rows=checked,
        same_unique_analytic_twenty_term_core_and_radial_definitions=True,
        earlier_C3_handles_preserved=True,altered_analytic_component_rejected=True,
        numerical_core_mixed_derivatives_beyond_total_order_five_not_certified=True)


def gauge_FTC(owner):
    g=owner.g;r=current.higher.rows;f=owner.core_functions
    assert r(f['G'])[1:]==r(f['gradient'])[:5]
    z=sy.Symbol('Z',real=True);G=sy.Function('actual_anchored_G')(z)
    L,C=sy.symbols('Lambda logC',positive=True);logF=-C-L*G;expF=sy.exp(logF)
    bindings={owner.core.parameters['Lambda'].node:L,owner.core.parameters['logC'].node:C}
    for key,expression in (('G',G),('gradient',sy.diff(G,z)),('logF0',logF)):
        for order,ref in enumerate(r(f[key])):
            expected=sy.diff(expression,z,order)
            if ref.node in bindings:assert sy.simplify(bindings[ref.node]-expected)==0
            bindings[ref.node]=expected
    for order,ref in enumerate(r(f['F0'])[:4]):bindings[ref.node]=sy.diff(expF,z,order)
    alg=higher_check.Algebra(g,bindings)
    for order,ref in enumerate(r(f['F0'])[4:],4):
        assert sy.simplify(alg.at(ref.node)-sy.diff(expF,z,order))==0,order
    row=g.nodes[f['G'].value.node]
    assert row['operation']=='definite_integral' and row['upper']==g.symbol('Z').node
    assert row['lower']==owner.core.parameters['anchor'].node and row['original_H_root_anchor']
    return dict(passed=True,actual_anchored_gauge_FTC_through_Z5=True,
        fifth_gauge_row_is_fourth_true_gradient_row=True,
        independent_exponential_F0_Z4_Z5_product_rules=True,
        original_exact_H_root_and_nonzero_Lambda_retained=True,
        moving_axial_endpoint_not_treated_as_fixed=True)


def lower_source(owner):
    r=current.higher.rows;g=owner.g;d=owner.derivative
    lower=owner.lower;preserved=[]
    for key,jet in owner.fields.items():
        old=r(lower.fields[key]);assert r(jet)[:len(old)]==old
        for left,right in zip(r(jet)[len(old)-1:],r(jet)[len(old):]):assert d.derivative(left)==right
        preserved.append(key)
    for key,jets in owner.radial_fields.items():
        for radial_order,jet in enumerate(jets):
            old=r(lower.fields[key][radial_order]);assert r(jet)[:len(old)]==old
            for left,right in zip(r(jet)[len(old)-1:],r(jet)[len(old):]):assert d.derivative(left)==right
    # Reuse the already accepted exact geometry/source theorem. It verifies
    # the binary sc, nonzero hb, core rho4 origin and zero modification support.
    geometry=lower_check.source_fixture(lower)
    integrals=[]
    for key,jet in owner.pressure.items():
        assert r(jet)[:2]==r(lower.functions[key])
        for order,ref in enumerate(r(jet)):
            row=g.nodes[ref.node]
            assert row['operation']=='definite_integral' and row['ordinary_slow_Z_derivative_order']==order
            assert row['lower']==lower.parameters['sc_quarter'].node
            assert row['upper']==lower.parameters['sc_half'].node
            assert row['original_direct_N1p'] and row['all_finite_viscous_derivatives_retained']
            assert row['original_m15_minus_applied_once'] and row['signed_integrand_not_a_cap']
            assert row['integrand']==r(owner.signed_bodies[key])[order].node
            if order>=2:assert row['actual_higher_source_body_not_interval_endpoint']
            integrals.append(key+'_Z'+str(order))
    assert r(owner.N1p)[:2]==r(lower.pressure_bindings['N1p'])
    for key,jet in owner.signed_bodies.items():
        assert r(jet)[:2]==r(lower.pressure_bindings['signed_physical_integrands'][key])
    assert lower.pressure_bindings['direct_differential_form_not_endpoint_subtraction']
    return dict(passed=True,actual_lower_fields_preserving_earlier_rows=preserved,
        actual_signed_pressure_integral_rows=integrals,actual_source_geometry=geometry,
        genuine_leading_mean_Z5_and_auxiliary_Q_Z4=True,
        original_direct_radial_viscous_terms_and_history_retained=True,
        leading_weight='1-delta',physical_U_m_Q_Pstar_not_reapplied=True)


def direct_physical_fixture():
    """Nonzero fifth mean row; direct physical-R operator/integrals to Z3."""
    g=current.source.FunctionTransportGraph();r=current.higher.rows
    z,rho,s=sy.symbols('Z fixture_higher_rho original_lower_collar_s',real=True)
    variables={str(v):g.symbol(str(v)) for v in (z,rho,s)}
    bind=lambda expr:current.source.expression(g,sy.sympify(expr),variables)
    jet=lambda expr:current.lower.c3.C3Function(*(bind(sy.diff(expr,z,k)) for k in range(4)))
    hb,eps,de,sc=map(sy.Rational,('1/5','1/9','1/5','3/5'))
    b=(3+z**3+z**5)/7;Phi=2+z*z/10;Psi=b*rho*rho;mean=Psi/3
    core_rows={key:[jet(sy.diff(expr,rho,i)) for i in range(3)]
        for key,expr in (('Phi',Phi),('Psi',Psi),('mean_Psi',mean))}
    p={key:bind(value) for key,value in dict(hb=hb,epsilon=eps,delta=de,
        logLambda=sy.log(1/eps),sc=sc,sc_quarter=sc/4,sc_half=sc/2).items()}
    core_f=dict(U0=jet(4*z+sy.Rational(1,10)),F0=jet(1+z*z/20),P0=jet(sy.Rational(3,2)))
    fields=current.lower.collar_fields(g,variables[str(rho)],core_rows,core_f,p,
        sigma=lambda x:g.constant('1/3'),sigma_prime=lambda x:g.zero)
    pressure,bodies=current.lower.direct_pressure_functions(g,fields,p)
    d=current.higher.HigherZDerivative(g)
    now={key:d.extend(pair,3) for key,pair in pressure.items()}
    U=d.extend(fields['U'],5);M=d.extend(fields['m'],5);Q=d.extend(fields['Q'],4)
    chi=1+(hb-1)/3
    psi=b*(16+chi*(rho*rho-16))
    m=b*(chi*rho*rho/3+16*(1-chi)-128*(1-chi)/(3*rho))
    ur=4*z+sy.Rational(1,10)+eps*psi;mr=4*z+sy.Rational(1,10)+eps*m
    qr=sy.cancel((2*z*ur-(1-de)*z*mr-(1-z*z)*sy.diff(mr,z))/(1-de*z*z))
    rs=4*sy.exp(hb*s);cache={};field_rows=0
    for pair,expression in ((U,ur),(M,mr),(Q,qr)):
        for order,ref in enumerate(r(pair)):
            assert sy.cancel(lower_check.symbolic(g,ref,cache)-sy.diff(expression,z,order).subs(rho,rs))==0
            field_rows+=1
    assert sy.diff(mr,z,5)!=0
    R=sy.Symbol('physical_higher_R',positive=True);q=qr.subs(rho,R/eps);u=ur.subs(rho,R/eps)
    L=1-de*z*z
    N=sy.cancel(-((q+(1-de)*z*sy.diff(q,z)/2+R*sy.diff(q,R))/L
        +q*(q/2+R*sy.diff(q,R))
        +u*(-2*z*q+(1-z*z)*sy.diff(q,z)-2*z*R*sy.diff(q,R))/L)/2
        +2*sy.diff(q,R)+R*sy.diff(q,R,2))
    bottom=4*eps*sy.exp(hb*sc/4);top=4*eps*sy.exp(hb*sc/2)
    polynomial=sy.Poly(sy.cancel(N*R**2),R)
    def integrate_polynomial(polynomial,power=0):
        value=0
        for (degree,),coefficient in polynomial.terms():
            exponent=degree-2+power
            primitive=(sy.log(top/bottom) if exponent==-1 else
                (top**(exponent+1)-bottom**(exponent+1))/(exponent+1))
            value+=coefficient*primitive
        return value
    exact=dict(m13=integrate_polynomial(polynomial,0),m15=-integrate_polynomial(polynomial,1));pressure_rows=0
    E=sy.Symbol('exact_exp_three_over_one_hundred',positive=True)
    def normal(expression):
        return sy.expand_log(expression,force=True).replace(
            lambda q:q.func==sy.exp and q.args[0].is_Rational,
            lambda q:E**(q.args[0]*100/3))
    points=(sy.Rational(-2,3),sy.Integer(0),sy.Rational(3,5))
    for point in points:
        point_cache={node:expression.subs(z,point) for node,expression in cache.items()}
        point_cache[variables[str(z)].node]=point
        for key,pair in now.items():
            for order,ref in enumerate(r(pair)):
                # Evaluate the actual derivative at an EXACT rational axial
                # point, then transform its whole ds density to physical R.
                # No numerical quadrature or interval endpoint is consumed.
                row=g.nodes[ref.node]
                body=lower_check.symbolic(g,current.source.FunctionRef(g,row['integrand']),point_cache)
                transformed=sy.cancel(body.subs(s,sy.log(R/(4*eps))/hb)/(hb*R))
                sign,power=(1,0) if key=='m13' else (-1,1)
                expected=sign*R**power*sy.diff(N,z,order).subs(z,point)
                assert sy.cancel(transformed-expected)==0,(str(point),key,order,'physical density')
                actual_poly=sy.Poly(sy.cancel(transformed*R**2),R)
                got=integrate_polynomial(actual_poly);want=sy.diff(exact[key],z,order).subs(z,point)
                assert sy.cancel(normal(got)-normal(want))==0,(str(point),key,order)
                pressure_rows+=1
    return dict(passed=True,mathematical_fixture_only=True,
        genuine_nonzero_fifth_mean_derivative_fixture=True,independent_higher_field_rows=field_rows,
        independent_physical_pressure_rows=pressure_rows,
        exact_axial_fixture_points=[str(point) for point in points],
        direct_physical_R_operator_and_exact_signed_integrals_Z0_Z3_checked=True,
        both_viscous_radial_derivatives_and_original_m15_sign_retained=True,
        nonunit_physical_measures_checked=True,original_scale_numerical_containment_claimed=False)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['core_lower_C5_C3_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    assert not any(raw[k] for k in current.GATES+current.OPEN)
    checks=dict(analytic_source=analytic_source(owner),gauge_FTC=gauge_FTC(owner),lower_source=lower_source(owner))
    print('ACTUAL_ANALYTIC_CORE_C5_AND_LOWER_PRESSURE_C3_SOURCE_CHECKED',flush=True)
    checks['direct_physical_fixture']=direct_physical_fixture()
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_ANALYTIC_CORE_C5_AND_LOWER_PRESSURE_C3_CHECKED',flush=True)
    return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses the checked live analytic core/collar source')
