"""Physical Z2, source normalization and continuous Taylor-model checks."""
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_unchanged_axial_pressure_models as current
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import FactoredAlgebra


def contains(enclosure,value):
    lo,hi=current.ep(enclosure)
    if hasattr(value,'_mpi_'):a,b=current.ep(value)
    else:a=b=value
    return lo<=a and b<=hi


def independent_physical_Z2_fixture():
    c=MPIntervalContext();c.dps=110
    algebra=FactoredAlgebra(c,(c.ln(c.mpf('.125')),c.mpf(0),c.ln(3),c.ln(7)),[])
    z=current.op.Jet.variable(c,c.mpf(['-.4','.4']),5);rho=c.mpf(['4','4.02'])
    a=algebra.lift(z+z*z)+algebra.shift(algebra.lift(1-z),(0,0,1,0))
    b=algebra.lift(z*z/7+c.mpf('.25'))
    psi=[a*rho+b*rho*rho,a+2*b*rho,2*b]
    mean=[a*rho/2+b*rho*rho/3,a/2+2*b*rho/3,2*b/3]
    axis=dict(z=z,u=1+z+z*z,d=1-z*z,L=1-z*z*c.mpf('.1'))
    rows=dict(Psi=[current.op.truncate(v,n) for v,n in zip(psi,(5,4,3))],
        meanPsi=[current.op.truncate(v,n) for v,n in zip(mean,(5,4,3))],axis=axis)
    fields=current.prior.mean.factored_leading_rows(algebra,rows,SimpleNamespace(delta=c.mpf('.1')))
    F=[algebra.lift(1+z*z),algebra.lift(0,4),algebra.lift(0,3)]
    forcing=current.op.known_forcing(algebra,dict(F=F,Uz=fields['Uz'],Q=fields['Q'],mean=fields['mean']),
        algebra.lift(z),rho,c.mpf('.1'))['N1p']
    resolved=algebra.resolve(forcing)
    R,Z=s.symbols('R Z',real=True);rr=7*R;delta=s.Rational(1,10)
    A=Z+Z*Z+3*(1-Z);B=Z*Z/7+s.Rational(1,4)
    U=1+Z+Z*Z+(A*rr+B*rr*rr)/7
    M=1+Z+Z*Z+(A*rr/2+B*rr*rr/3)/7
    L=1-delta*Z*Z;d=1-Z*Z
    q=(2*Z*U-(1-delta)*Z*M-d*s.diff(M,Z))/L
    S=(q+(1-delta)*Z*s.diff(q,Z)/2+R*s.diff(q,R))/L+q*(q/2+R*s.diff(q,R))
    S+=U*(-2*Z*q+d*s.diff(q,Z)-2*Z*R*s.diff(q,R))/L
    N=-S/2+2*s.diff(q,R)+R*s.diff(q,R,2)
    count=0
    with mp.workdps(120):
        for zv in ('-.4','0','.4'):
            for rv in ('4','4.02'):
                for k in range(3):
                    point=mp.mpf(str(s.N(s.diff(N,Z,k).subs({Z:s.Rational(zv),R:s.Rational(rv)/7}),115)))
                    assert contains(resolved[k]*math.factorial(k),point),(zv,rv,k)
                    count+=1
    return dict(independent_physical_R_Z0_Z1_Z2_rows=count,nonunit_Pi_Lambda=True,
        nonlinear_U_and_own_mean=True,ordinary_Z2_factorial_verified=True)


def independent_Taylor_and_normalization_fixture():
    c=MPIntervalContext();c.dps=100;count=0
    # f=Z^3-Z has actual zero crossings. Verify value and derivative enclosures.
    for i in range(16):
        left=c.mpf(-1)+c.mpf(2)*i/16;right=c.mpf(-1)+c.mpf(2)*(i+1)/16
        center=(left+right)/2;curvature=3*c.mpf([current.ep(left)[0],current.ep(right)[1]])
        coefficients=[center**3-center,3*center**2-1]
        radius,_,_=current.affine_error(c,coefficients,(right-left)/2,curvature)
        for t in (0,c.mpf('.25'),c.mpf('.5'),c.mpf('.75'),1):
            z=left+(right-left)*t;offset=z-center
            assert contains(coefficients[0]+coefficients[1]*offset+c.mpf([-1,1])*radius,z**3-z)
            assert contains(coefficients[1]+2*curvature*offset,3*z*z-1)
            count+=2
    algebra=FactoredAlgebra(c,(c.ln(2),c.mpf(0),c.ln(3),c.ln(5)),[])
    value=algebra.lift(current.op.Jet(c,[c.mpf(['-2','1']),c.mpf(3),c.mpf(-4)]))
    value+=algebra.shift(algebra.lift(current.op.Jet(c,[c.mpf(1),c.mpf(-2),c.mpf(1)])),(1,0,0,-1))
    resolved=algebra.resolve(value)
    for k in range(3):
        got=current.normalized_row(value,k,(0,0,0,0),c.mpf(7))
        assert contains(got,resolved[k]/7)
    return dict(zero_crossing_cubic_value_and_derivative_checks=count,
        independent_nonunit_source_normalizations=3,no_relative_division_at_zeros=True)


@current.source_precision
def run(owner,data):
    began=time.monotonic();owner.assert_graph();c=owner.c
    assert data['source_family']==owner.family and data['input_hashes']==owner.hashes
    checks=dict(physical=independent_physical_Z2_fixture(),Taylor=independent_Taylor_and_normalization_fixture())
    assert data['exact_Z_domain']==['-1','1'] and data['continuous_domain_including_zeros']
    count=len(data['cells']);assert count>=4 and count&(count-1)==0
    assert not data[current.GATE] and all(not data[k] for k in current.OPEN)
    target=c.mpf(data['normalized_absolute_error_target']);maxerror=c.mpf(0)
    for index,cell in enumerate(data['cells']):
        left,right=[c.mpf(v) for v in cell['exact_Z_cell']];center=c.mpf(cell['exact_center'])
        assert left._mpi_==(c.mpf(-1)+2*c.mpf(index)/count)._mpi_
        assert right._mpi_==(c.mpf(-1)+2*c.mpf(index+1)/count)._mpi_
        assert center._mpi_==((left+right)/2)._mpi_
        for name,model in cell['models'].items():
            scale=data['normalizations'][name]
            assert scale['reference_source_exponents'][1]==0
            assert current.ep(current.read_interval(c,scale['directed_reference_coefficient_norm']))[0]>0
            coefficients=[current.read_interval(c,v) for v in model['affine_coefficients_Z0_Z1']]
            curvature=current.read_interval(c,scale['full_domain_Z2_Taylor_coefficient'])
            expected=current.affine_error(c,coefficients,(right-left)/2,curvature)
            for key,value in zip(('Taylor_remainder_abs_upper','coefficient_ambiguity_abs_upper','total_absolute_model_error_upper'),expected):
                assert current.read_interval(c,model[key])._mpi_==value._mpi_,(index,name,key)
            assert current.ep(expected[2])[1]<=current.ep(target)[0]
            maxerror=current.domain.upper(c,max(maxerror,expected[2]))
    full,_=current.original_values(owner,['-1','1'])
    for name,value in full.items():
        scale=data['normalizations'][name];ref=scale['reference_source_exponents'];M=current.read_interval(c,scale['directed_reference_coefficient_norm'])
        for k,key in enumerate(('full_domain_value_range','full_domain_first_derivative_range','full_domain_Z2_Taylor_coefficient')):
            assert current.normalized_row(value,k,ref,M)._mpi_==current.read_interval(c,scale[key])._mpi_
    points={}
    for query in ('-1','-.639','0','.211','.371','.639','1'):
        direct,_=current.original_values(owner,query);view=current.evaluate_models(c,data,query)
        for name,value in direct.items():
            scale=data['normalizations'][name];ref=scale['reference_source_exponents'];M=current.read_interval(c,scale['directed_reference_coefficient_norm'])
            for k,key in enumerate(('normalized_values','normalized_first_derivatives')):
                assert contains(view[key][name],current.normalized_row(value,k,ref,M)),(query,name,k)
            assert contains(current.read_interval(c,scale['full_domain_Z2_Taylor_coefficient']),
                current.normalized_row(value,2,ref,M)),(query,name,'Z2')
        if query=='.371':
            packet=owner.prefix.evaluate(query)
            original=owner.prefix._packets[id(packet)][2][3]
            assert all(current.op.snapshot(current.op.truncate(value,1))==current.op.snapshot(original[name])
                for name,value in direct.items())
        points[query]=True
    crossing=current.evaluate_models(c,data,['.211','.371'])
    assert all(current.ep(v)[0]<0<current.ep(v)[1] for v in crossing['normalized_values'].values())
    checks['coverage_and_error']=dict(exact_domain=['-1','1'],cells=count,pressure_value_models=2*count,
        maximum_uniform_absolute_model_error_upper=current.op.snapshot(maxerror),named_direct_queries=points,
        signed_zero_crossing_interval_enclosed=True,first_derivative_accuracy_not_claimed=True)
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
        scope='Rin to rminus only; both pressure values uniformly on real Z in [-1,1]',
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_UNCHANGED_FULL_Z_ABSOLUTE_PRESSURE_MODELS',count,flush=True)
    return result
