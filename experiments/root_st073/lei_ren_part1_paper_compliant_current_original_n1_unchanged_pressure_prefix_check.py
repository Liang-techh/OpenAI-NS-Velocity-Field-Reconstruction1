"""Independent source equations, physical forcing, and prefix integral checks."""
import json
from pathlib import Path
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_unchanged_pressure_prefix as current
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import FactoredAlgebra


def exact_collar_difference_fixture():
    r,z=s.symbols('r Z',real=True);A=1+z+z*z
    corephi=s.exp(A*r);corepsi=(1+2*z)*r+(2-z)*r*r
    chi=1+(r-4)/7;ratio=s.exp(A*(r-4)**2/14)
    p=s.Function('actual_Psi')(r);m=s.Function('actual_meanPsi')(r)
    pcmean=(1+2*z)*r/2+(2-z)*r*r/3
    dpsi0=p-corepsi;dmean0=m-pcmean
    got=current.collar_derivative_corrections(
        [corepsi,s.diff(corepsi,r),s.diff(corepsi,r,2)],
        [corephi,s.diff(corephi,r)],chi,s.diff(chi,r),ratio,dpsi0,dmean0,r)
    actual_p1=ratio*chi*s.diff(corepsi,r)
    actual_m1=(p-m)/r
    actual_m2=(actual_p1-2*actual_m1)/r
    expected=dict(dpsi1=actual_p1-s.diff(corepsi,r),
        dpsi2=s.diff(actual_p1,r)-s.diff(corepsi,r,2),
        dmean1=actual_m1-s.diff(pcmean,r),dmean2=actual_m2-s.diff(pcmean,r,2))
    for name,value in expected.items():
        assert s.simplify(got[name]-value)==0,name
        assert s.simplify(s.diff(got[name],z)-s.diff(value,z))==0,name+' Z1'
    assert s.simplify(s.diff(ratio,r)/ratio-(chi-1)*s.diff(corephi,r)/corephi)==0
    return dict(exact_difference_and_Z1_rows=8,nonconstant_core_phi=True,
                nonconstant_chi_and_original_ratio_equation=True,
                original_mean_FTC_derivative_equations_verified=True)


def flat_tail_theorem():
    x=s.Symbol('x',positive=True);D=1/x**2-1/(1-x)**2
    assert s.simplify(s.diff(D,x)+2/x**3+2/(1-x)**3)==0
    # On 0<x<=1/2, (1-x)^-2<=4 and (1-x)^-3<=x^-3.
    # sigma=e^-D/(1+e^-D); sigma'=sigma(1-sigma)*(-D').
    # The majorants exp(4-x^-2)*x^-3j increase on this left interval.
    for j in (0,1):
        derivative=s.diff(4-1/x**2-3*j*s.log(x),x)
        assert s.simplify(derivative-(2-3*j*x*x)/x**3)==0
        assert 2-s.Rational(3*j,4)>0
    with mp.workdps(80):
        for point in ('0.1','0.25','0.5'):
            q=mp.mpf(point);E=mp.exp(-1/q**2+1/(1-q)**2)
            sigma=E/(1+E);sigma1=sigma*(1-sigma)*(2/q**3+2/(1-q)**3)
            assert 0<sigma<=mp.exp(4-1/q**2)
            assert 0<sigma1<=4*mp.exp(4-1/q**2)/q**3
    return dict(original_sigma_and_first_derivative_majorants=True,
                monotone_left_tail_for_orders=[0,1],nonzero_sigma_not_replaced_by_zero=True,
                Gronwall_bound='|exp(w)-1|<=2|w| for |w|<=1/2')


def _contains(c,enclosure,point):
    lo,hi=current.ep(enclosure)
    return lo<=point<=hi


def independent_physical_forcing_fixture():
    """Resolve only small fixture scales, then compare the physical R operator."""
    c=MPIntervalContext();c.dps=100
    logs=(c.ln(c.mpf('0.125')),c.mpf(0),c.ln(3),c.ln(7))
    algebra=FactoredAlgebra(c,logs,[])
    z=current.op.Jet.variable(c,c.mpf('.371'),5);rho=c.mpf(['4','4.02'])
    a=algebra.lift(z+z*z)+algebra.shift(algebra.lift(1-z),(0,0,1,0))
    b=algebra.lift(z*z/7+c.mpf('0.25'))
    psi=[a*rho+b*rho*rho,a+2*b*rho,2*b]
    meanpsi=[a*rho/2+b*rho*rho/3,a/2+2*b*rho/3,2*b/3]
    axis=dict(z=z,u=-4*z+c.mpf('.25'),d=1-z*z,L=1-z*z*c.mpf('.1'))
    rows=dict(Psi=[current.op.truncate(v,n) for v,n in zip(psi,(5,4,3))],
        meanPsi=[current.op.truncate(v,n) for v,n in zip(meanpsi,(5,4,3))],axis=axis)
    fields=current.mean.factored_leading_rows(algebra,rows,SimpleNamespace(delta=c.mpf('.1')))
    f=[algebra.lift(1+z*z),algebra.lift(0,4),algebra.lift(0,3)]
    raw=current.op.known_forcing(algebra,dict(F=f,Uz=fields['Uz'],Q=fields['Q'],mean=fields['mean']),
        algebra.lift(z),rho,c.mpf('.1'))['N1p']
    resolved=algebra.resolve(raw)
    R,Z=s.symbols('R Z',positive=True);delta=s.Rational(1,10);Lambda=s.Integer(7)
    rr=Lambda*R;A=Z+Z*Z+3*(1-Z);B=Z*Z/7+s.Rational(1,4)
    U=-4*Z+s.Rational(1,4)+(A*rr+B*rr*rr)/Lambda
    M=-4*Z+s.Rational(1,4)+(A*rr/2+B*rr*rr/3)/Lambda
    L=1-delta*Z*Z;d=1-Z*Z
    Q=(2*Z*U-(1-delta)*Z*M-d*s.diff(M,Z))/L
    S=(Q+(1-delta)*Z*s.diff(Q,Z)/2+R*s.diff(Q,R))/L+Q*(Q/2+R*s.diff(Q,R))
    S+=U*(-2*Z*Q+d*s.diff(Q,Z)-2*Z*R*s.diff(Q,R))/L
    N=-S/2+2*s.diff(Q,R)+R*s.diff(Q,R,2)
    rows_checked=0
    for radial in ('4','4.011','4.02'):
        for k in (0,1):
            exact=s.diff(N,Z,k).subs({Z:s.Rational(371,1000),R:s.Rational(radial)/7})
            with mp.workdps(110):
                point=mp.mpf(str(s.N(exact,105)))
            assert _contains(c,resolved[k],point),(radial,k)
            rows_checked+=1
    return dict(independent_physical_R_pressure_and_Z1_rows=rows_checked,
                nonconstant_U_and_own_mean=True,Pi_and_Lambda_scales_nonunit=True,
                actual_source_scales_not_resolved_in_production=True)


def physical_measure_fixture():
    q,h,L,a,b=s.symbols('q hb Lambda a b',positive=True)
    r=4*s.exp(h*q);R=r/L
    assert s.simplify(s.diff(R,q)-h*r/L)==0
    assert s.simplify(R*s.diff(R,q)-h*r*r/L**2)==0
    class Algebra:
        def shift(self,value,powers):
            return value*h**powers[0]*L**powers[3]
        def width(self,value,power=1):
            return value*h**power
    N,rho=s.symbols('N rho',real=True)
    got=current.physical_pressure_range_integrals(Algebra(),N,rho,a,b)
    assert s.simplify(got['m13']-N*h*rho*(b-a)/L)==0
    assert s.simplify(got['m15']+N*h*rho*rho*(b-a)/L**2)==0
    return dict(original_hb_in_both_physical_measures_once=True,
                m15_minus_weighted_pressure_sign=True)


@current.source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph()
    fixture=dict(collar=exact_collar_difference_fixture(),sigma=flat_tail_theorem(),
                 forcing=independent_physical_forcing_fixture(),measure=physical_measure_fixture())
    assert len(packets)==3
    samples={}
    for name,packet in packets.items():
        view=owner.report(packet)
        assert view['source_family']==owner.family
        assert view['original_modification_support']['increments_exact_zero']
        assert view['query']['lower_is_actual_Rin'] and view['query']['upper_is_actual_rminus']
        assert view['direct_same_source_pressure_integral_not_endpoint_subtraction']
        assert view['original_hb_measure_retained_once']
        assert view['absolute_radius_hb_and_N_not_materialized']
        assert view['whole_Z_precision_not_certified']
        for proof in view['actual_source_reference_and_error_proof']['original_sigma_bounds']:
            assert current.ep(proof['log_sigma_over_hb_squared_upper'])[1]<current.ep(owner.c.ln(proof['nonzero_coefficient_upper']))[0]
            assert proof['source_hb_power']==2 and proof['source_value_not_materialized_or_replaced']
        for field,rows in view['signed_source_scaled_pressure_budgets'].items():
            for key,row in rows.items():
                assert row['relative_partial_moment_width_certified'],(name,field,key)
                assert row['sign'] in ('positive','negative')
                assert current.ep(row['total_width_over_abs_lower'])[1]<current.ep(owner.c.mpf('.001'))[0]
                assert row['actual_leading_source_and_collar_error_enclosure_used']
                scaled=view['signed_source_scaled_pressure_values'][field][key]
                assert current.ep(scaled['directed_reference_coefficient_norm'])[0]>0
                assert scaled['source_bound_normalization_not_selected_value']
        assert all(not view[key] for key in current.OPEN)
        samples[name]=dict(Z=str(current.ep(view['Z'])[0]),pressure_rows={field:{key:dict(sign=row['sign'],
            relative_width_upper=str(current.ep(row['total_width_over_abs_lower'])[1])) for key,row in rows.items()}
            for field,rows in view['signed_source_scaled_pressure_budgets'].items()})
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=fixture,
        actual_samples=samples,numeric_scope='three named original Z points, pressure values and Z1',
        original_comparison_errors_and_pressure_atoms_retained=True,
        whole_Z_relative_accuracy_certified=False,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_ACTUAL_UNCHANGED_PRESSURE_PREFIX',len(packets),flush=True)
    return result


if __name__=='__main__':
    raise SystemExit('Use run(owner, packets)')
