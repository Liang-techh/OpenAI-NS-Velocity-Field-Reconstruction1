"""Evaluate full original n1 axial Z2 with separate source/kernel/word errors."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from openai_ns_reconstruction.original_n1_core_axial2 import (
    OriginalN1CoreAxial2,forcing_difference_z2,kernel_tail_z2)
from openai_ns_reconstruction.original_n1_core_volterra_numeric import OriginalN1CoreVolterraNumeric,DifferenceJet
from openai_ns_reconstruction.original_core_numeric import _bounds

DIRECTORY=ROOT/'experiments/root_st073'
NAME='lei_ren_part1_paper_compliant_current_original_n1_core_axial2.json'


def raw_source_identities(solver):
    r,z=sy.symbols('rho Z');R=solver.R
    de,e,P2=sy.Rational(1,17),sy.Rational(1,3),sy.Integer(5)
    A=1+z*r/7+r*r/11
    U=2+z/3+r*z*z/5+r*r*z**3/9
    da=r*(1+z+z*z)/13+r**3*z/17
    du=r*r*(1+z+z**4)/19
    average=lambda v:sy.cancel(sy.integrate(v,r)/r)
    at={r:sy.Rational(1,2),z:sy.Rational(1,5)}
    def local(value):
        return R([[sy.diff(value,r,i,z,k).subs(at)/(math.factorial(i)*math.factorial(k))
                   for k in range(5)] for i in range(3)],2,4,sy.sympify)
    pair=lambda value,error=0:DifferenceJet(local(value),local(error))
    result=forcing_difference_z2(pair(A,da),pair(U,du),pair(average(U),average(du)),
        pair(z),pair(r),pair(de),pair(e),pair(P2),pair(e*(2+2*z)))
    d,L=1-z*z,1-de*z*z;F0=sy.exp(-z*z-2*z)
    T=lambda weight,value:(weight*z*value+d*sy.diff(value,z)-2*z*r*sy.diff(value,r))/L
    def raw(a,u):
        mean=average(u);q=(2*z*u-(1-de)*z*mean-d*sy.diff(mean,z))/L
        flow=(q+(1-de)*z*sy.diff(q,z)/2+r*sy.diff(q,r))/L
        flow+=q*(q/2+r*sy.diff(q,r))+u*T(-2,q)
        np=-flow/(2*P2)+(2*sy.diff(q,r)+r*sy.diff(q,r,2))/(e*P2)
        return dict(f=-2*e*e*T(-3,T(-2-de,F0*a))/F0,
                    u=-2*T(-2,T(-1-de,u))-4*z*r*e*P2*np/L,p=np)
    old,new=raw(A,U),raw(A+da,U+du);count=0
    for name,row in result.items():
        for k in range(3):
            value=lambda expr:sy.diff(expr,z,k).subs(at)/math.factorial(k)
            assert sy.simplify(row.base.rows[0][k]-value(old[name]))==0,(name,k,'base')
            assert sy.simplify(row.error.rows[0][k]-value(new[name]-old[name]))==0,(name,k,'change')
            count+=2
    return count


def positive_kernel_tail_checks(solver):
    c=mp.mp.clone();c.dps=170;iv,J=solver.iv,solver.leading.taylor
    chi=J(iv,[1,iv.mpf(1)/7,iv.mpf(1)/11])
    count=0
    for rho in (1,4):
        for K in (8,16,24):
            bound=kernel_tail_z2(iv,chi,str(rho),K)
            q,x=c.mpf(rho)/2,c.sqrt(rho)
            tail=[c.mpf(0)]*3;tailx=[c.mpf(0)]*3
            for n in range(K+1,K+101):
                t=q**n/(math.factorial(n)*math.factorial(n+1))
                coeff=[c.mpf(1),c.mpf(n)/7,c.mpf(n)/11+c.mpf(n*(n-1))/98]
                for k in range(3):tail[k]+=x*t*coeff[k]/2;tailx[k]+=n*t*coeff[k]
            for name,reference in (('f',tail),('fx',tailx)):
                for k in range(3):
                    assert reference[k]<_bounds(c,bound[name][k])[1],(rho,K,name,k)
                    count+=1
    return count


def compare(a,b):
    c=mp.mp.clone();c.dps=190
    delta=abs(c.mpf(a['approximation'])-c.mpf(b['approximation']))
    error=c.mpf(a['absolute_error_upper'])+c.mpf(b['absolute_error_upper'])
    assert delta<=error,(delta,error)


def solve(field,degree,kernel):
    result=[]
    for Z in ('-.3','0','.3'):
        for rho in ('0','1','2','4'):result.append(field.evaluate(Z,rho,degree=degree,kernel_degree=kernel))
        print('Original n1 core Z2: precision='+str(field.precision)+', Z='+Z,flush=True)
    return result


def run():
    start=time.perf_counter()
    coarse=OriginalN1CoreAxial2(precision=55,quadrature_nodes=96)
    exact=raw_source_identities(coarse);kernels=positive_kernel_tail_checks(coarse)
    print('Independent Z2 raw-source and positive kernel-tail checks pass',flush=True)
    first=solve(coarse,14,18)
    fine=OriginalN1CoreAxial2(precision=85,quadrature_nodes=192)
    second=solve(fine,20,28);comparisons=legacy_count=0
    relative={name:fine.ctx.mpf(0) for name in fine.fields}
    for old,new in zip(first,second):
        for name,rows in new['normalized_original_n1'].items():
            for k,(a,b) in enumerate(zip(rows,old['normalized_original_n1'][name])):
                compare(a,b);comparisons+=1
                if k==2 and name in relative and fine.ctx.mpf(a['approximation']):
                    error=fine.ctx.mpf(a['absolute_error_upper'])/abs(fine.ctx.mpf(a['approximation']))
                    assert error<fine.ctx.mpf('1e-12'),(new['Z'],new['rho'],name,error)
                    relative[name]=max(relative[name],error)
        if new['rho']=='0':
            assert all(row['original_source_enclosure_exact_mpf_tuples']==[[0,0,0,0],[0,0,0,0]]
                for rows in new['normalized_original_n1'].values() for row in rows)
        if new['rho']=='4':
            legacy=OriginalN1CoreVolterraNumeric.evaluate(fine,new['Z'],'4',degree=20,kernel_degree=28)
            for name in fine.fields:
                for a,b in zip(new['normalized_original_n1'][name][:2],legacy['normalized_original_n1'][name]):
                    compare(a,b);legacy_count+=1
            for name in ('mean_u','q'):
                compare(new['normalized_original_n1'][name][0],legacy['normalized_original_n1'][name]);legacy_count+=1
    rejected=0
    for Z,rho,kwargs in (('1','2',{}),('.3','4.1',{}),('.3','-1',{}),
                         ('.3','2',dict(degree=7)),('.3','2',dict(kernel_degree=7))):
        try:fine.evaluate(Z,rho,**kwargs)
        except ValueError:rejected+=1
        else:raise AssertionError(('Unsupported Z2 request admitted',Z,rho,kwargs))
    implementation=dict(fine.code_hashes)
    implementation[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=fine.family,checks=dict(passed=True,
        independent_raw_source_Z0_Z1_Z2_difference_identities=exact,
        independent_positive_kernel_Z2_tail_checks=kernels,
        original_core_packets_per_level=len(second),full_error_precision_degree_kernel_comparisons=comparisons,
        accepted_Z1_interface_comparisons=legacy_count,exact_axis_all_Z_coefficients_zero=True,
        maximum_sampled_Z2_Taylor_coefficient_relative_errors={name:fine.ctx.nstr(v,35) for name,v in relative.items()},
        rejected_unsupported_requests=rejected,feedback_Z2_ordinary_tail_divided_by_two=True,
        source_kernel_and_feedback_errors_separate=True,unknown_Z2_feedback_iterate_not_zero_padded=True),
        levels=[dict(precision=55,source_degree=14,kernel_degree=18,samples=first),
                dict(precision=85,source_degree=20,kernel_degree=28,samples=second)],
        accepted_all_word_tail=fine.tail,input_hashes=fine.hashes,implementation_hashes=implementation,
        limitations=dict(core_only=True,full_n1_state_Taylor_order=2,mean_u_Taylor_order=2,q_Taylor_order=1,
            nominal_feedback_iterate_computed=False,collar_Z2_transport=False,
            physical_mixed_derivatives=False,global_moment_repair=False,
            temporal_n_ge_2_recursion=False,full_corrected_NS=False),
        execution_seconds=round(time.perf_counter()-start,3))


def inspect(path):
    report=json.loads(path.read_text(encoding='utf-8'));assert report['checks']['passed']
    for name,digest in report['input_hashes'].items():
        assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in report['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true')
    parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    report=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=report['checks'],execution_seconds=report['execution_seconds']),indent=2),flush=True)
