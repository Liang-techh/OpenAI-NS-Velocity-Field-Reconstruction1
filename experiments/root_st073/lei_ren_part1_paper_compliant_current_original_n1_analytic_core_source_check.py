"""Independent coefficient/forcing identities and original core admission."""
from fractions import Fraction
import json
import math
from pathlib import Path
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_analytic_core_source as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def contains(got,want):return ep(got)[0]<=ep(want)[0] and ep(got)[1]>=ep(want)[1]


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported analytic core source admitted')


def exact_model_identities():
    coefficients=[Fraction((-1)**n,2**n*math.factorial(n)*math.factorial(n+1)) for n in range(50)]
    count=0
    for n in range(49):
        assert 4*(n+1)*(n+2)*coefficients[n+1]+2*coefficients[n]==0
        assert 4*(n+1)**2*(n+2)*coefficients[n+1]+2*n*coefficients[n]==-2*coefficients[n]
        count+=2
    assert coefficients[1]==Fraction(-1,4)
    return dict(independent_exact_Bessel_and_chi_derivative_ODE_rows=count,
                axis_value_and_x_derivative_zero_but_rho_derivative_minus_quarter=True)


def forcing_fixture(c):
    """Differentiate the physical F0*Phi first, independently of the adapter."""
    z,r=s.symbols('z rho');delta=s.Rational(1,11);lam=s.Integer(7)
    grad=s.Rational(1,3)+z/5+z*z/7
    amplitude=s.exp(-lam*(z/3+z*z/10+z**3/21))
    phi=1+r*z+r*r*z*z+z**3+r**3/17
    d=1-z*z;L=1-delta*z*z;a=-2-delta;b=-3
    def axial(value,weight):return (weight*z*value+d*s.diff(value,z)-2*r*z*s.diff(value,r))/L
    expected=s.cancel(s.expand(-2*axial(axial(amplitude*phi,a),b)/(amplitude*lam**2)))
    count=0
    precise=type(c)();precise.dps=2*c.dps
    for z0,r0 in ((s.Rational(0),s.Rational(0)),(s.Rational(1,5),s.Rational(3,2)),(s.Rational(-2,7),s.Rational(4))):
        rational=lambda ctx,v:ctx.mpf(int(v.p))/int(v.q)
        def jet(expr,order):
            return op.Jet(c,[rational(c,s.cancel(s.diff(expr,z,k).subs({z:z0,r:r0})/math.factorial(k))) for k in range(order+1)])
        algebra=op.N1LogAlgebra(c,(c.mpf(0),c.mpf(0),c.mpf(0),c.ln(7)),[])
        rows=[algebra.lift(jet(s.diff(phi,r,i),3-i)) for i in range(3)]
        sectors=current.normalized_theta_sectors(algebra,rows,algebra.lift(jet(z,3)),rational(c,r0),
                                                  rational(c,delta),algebra.lift(jet(grad,3)))
        actual=sum(sectors.values(),algebra.lift(0,1))
        assert set(actual.terms)<={(0,0,0,0),(0,0,0,-1),(0,0,0,-2)}
        for k in (0,1):
            want=rational(precise,s.cancel(s.diff(expected,z,k).subs({z:z0,r:r0})))
            got=sum((row[k]*c.mpf(7)**powers[3] for powers,row in actual.terms.items()),c.mpf(0))
            assert contains(got,want),(z0,r0,k);count+=1
    return dict(independent_physical_source_normalization_derivative_rows=count,
                full_m_and_L_derivatives_and_radial_rho_coordinate_retained=True)


def factorial_tail_fixture(c):
    count=0
    for ch in (0,c.mpf(1)/5,1):
        chi=op.Jet(c,[ch,c.mpf(1)/3])
        for radius in (0,2,4):
            low,tail=current.model_chi_derivatives(c,chi,c.mpf(radius),8)
            high,high_tail=current.model_chi_derivatives(c,chi,c.mpf(radius),96)
            for field in current.FIELDS:
                for k in (0,1):
                    got=low[field][k]+c.mpf([-1,1])*tail[field][k]
                    want=high[field][k]+c.mpf([-1,1])*high_tail[field][k]
                    assert contains(got,want),(ch,radius,field,k);count+=1
    # Independent exact coefficient ratios, including the extra n factors.
    for n in range(9,30):
        c0=Fraction(1,2**n*math.factorial(n-1)*math.factorial(n+1))
        c1=Fraction(1,2**(n+1)*math.factorial(n)*math.factorial(n+2))
        assert c1/c0==Fraction(1,2*n*(n+2))
        assert (n+1)*c1/(n*c0)==Fraction(n+1,2*n*n*(n+2))
    return dict(factorial_tail_containment_rows=count,independent_model_tail_ratio_rows=42)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**exact_model_identities(),**forcing_fixture(c),**factorial_tail_fixture(c))
    for name,packet in packets.items():
        view=owner.report(packet)
        assert all(not view[k] for k in current.OPEN)
        assert view['core_initial_source_integral_only']
        assert view['remaining_normalized_B0_B1_feedback_not_solved']
        assert view['actual_first_collar_integral_still_open']
        assert view['original_nonlinear_field_enclosed_not_replaced_by_model']
        values,errors=owner._packets[id(packet)][3:5]
        for value in values.values():assert set(value.terms)<={current.ZERO}
        for error in errors.values():
            assert set(error.terms)<={(0,0,0,0),(0,0,0,-1),(0,0,0,-2)}
            assert all(ep(v)[0]>=0 for row in error.terms.values() for v in row.coefficients)
        if name=='core_axis':
            assert all(not value.terms for value in (*values.values(),*errors.values()))
            assert all(row['initial_integral_exact_zero'] for f in view['normalized_initial_integral_budgets'].values() for row in f.values())
        else:
            assert view['passing_normalized_initial_width_rows']==4
            for label in ('normalized_initial_integral_budgets','physical_initial_integral_budgets'):
                assert all(row['relative_initial_integral_width_certified'] for f in view[label].values() for row in f.values())
    checks['actual_original_core_initial_integral_four_width_rows_pass_at_both_nonzero_queries']=True
    checks['actual_axis_exact_zero_and_original_error_sectors_retained']=True
    checks['physical_F0_Lambda_weight_and_ordinary_derivative_errors_transferred']=True
    checks['actual_collar_or_coupled_completion_not_promoted']=True
    checks['unsupported_endpoint_or_model_inputs_rejected']=all((
        rejected(lambda:owner.evaluate('.371','4.1')),
        rejected(lambda:owner.evaluate('.371',c.mpf([0,4]))),
        rejected(lambda:current.model_chi_derivatives(c,op.Jet(c,[1,0]),c.mpf(4),True)),
        rejected(lambda:owner.report(current.OriginalN1AnalyticCoreSourcePacket()))))
    entry=owner._packets[id(packets['core_exit'])];old=entry[1]['core_initial_source_integral_only']
    entry[1]['core_initial_source_integral_only']=False
    try:checks['false_scope_mutation_rejected']=rejected(lambda:owner.report(packets['core_exit']))
    finally:entry[1]['core_initial_source_integral_only']=old
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
                **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_N1_ANALYTIC_CORE_SOURCE',checks['independent_exact_Bessel_and_chi_derivative_ODE_rows'],
          checks['independent_physical_source_normalization_derivative_rows'],checks['factorial_tail_containment_rows'],flush=True)
    return result
