"""Focused independent paper-equation and same-source input checks."""
import copy
import json
import math
from pathlib import Path
import time

import sympy as s
import lei_ren_part1_paper_compliant_current_original_n1_inner_operator as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def contains(got,want):
    a,b=current.ep(got);c,d=current.ep(want)
    return a<=c and b>=d


def overlaps(a,b):
    x,y=current.ep(a);u,v=current.ep(b)
    return max(x,u)<=min(y,v)


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported source/scope was admitted')


def exact_fixture(c):
    """Differentiate direct paper equations with a nonconstant F0 factor."""
    r,z=s.symbols('rho Z');delta=s.Rational(1,17);lam=s.Integer(3);eps=1/lam
    f=(1+z/5)*(1+r/7+r*z/11+r*r/13)
    w=2+z/3+r*z/5+r*r*(1+z)/19
    q=1+z*z/7+r*z/9+r*r/17
    L=1-delta*z*z;d=1-z*z
    za=lambda a,h:(a*z*h+d*s.diff(h,z)-2*z*r*s.diff(h,r))/L
    theta=-za(-3,za(-2-delta,f));axial=-za(-2,za(-1-delta,w))
    timepart=(q+(1-delta)*z*s.diff(q,z)/2+r*s.diff(q,r))/L
    pressure=-(timepart+q*(q/2+r*s.diff(q,r))+w*za(-2,q))/2+lam*(2*s.diff(q,r)+r*s.diff(q,r,2))
    X=2;zz=s.Rational(7,19);rr=s.Integer(4)
    rational=lambda v:c.mpf(int(v.p))/int(v.q)
    algebra=current.N1LogAlgebra(c,(c.mpf(0),)*4,[])
    def lift_expression(expr,order):
        return algebra.lift(current.Jet(c,[rational(s.diff(expr,z,k).subs({r:rr,z:zz}))/math.factorial(k) for k in range(order+1)]))
    rows={name:[lift_expression(s.diff(expr,r,i),3-i) for i in range(3)] for name,expr in (('F',f),('Uz',w),('Q',q))}
    zj=algebra.lift(current.Jet.variable(c,rational(zz),3));de=rational(delta)
    forcing=current.known_forcing(algebra,rows,zj,c.mpf(int(rr)),de)
    matrix=current.matrix(algebra,rows,forcing,zj,c.mpf(int(rr)),de)
    # These small rational bases belong ONLY to this exact synthetic fixture.
    # No production source logarithm is exponentiated.
    bases=(c.mpf(1),c.mpf(1),c.mpf(1),c.mpf(int(lam)))
    count=0
    def check(value,expression,order):
        nonlocal count
        for k in range(order+1):
            got=sum((jet[k]*math.factorial(k)*math.prod(bases[i]**p for i,p in enumerate(key))
                     for key,jet in value.terms.items()),c.mpf(0))
            want=rational(s.diff(expression,z,k).subs({r:rr,z:zz}))
            assert contains(got,want),(k,current.ep(got),current.ep(want))
            count+=1
    for name,expr in (('N1theta',theta),('N1z',axial),('N1p',pressure)):check(forcing[name],expr,1)
    F,U,K,P,YF,YU=s.symbols('F U K P YF YU');Fz,Uz,Kz,Pz=s.symbols('Fz Uz Kz Pz')
    H=(1-delta)*z/2+d*w;B=q+(1-2*z*w)/L;E=z*w-s.Rational(1,2)
    vn=((1-delta)*z*U-(1+delta)*z*K-d*(Uz+Kz))/L
    px=eps*(4*X*f*F+2*X*pressure)
    direct={1:YF,2:YU,3:-YU,4:px,
      5:eps*(X*B*YF+2*H*Fz/L+2*(q+(-2+delta)*E/L)*F+2*(f+r*s.diff(f,r))*vn+2*za(-2-delta,f)*U+2*theta),
      6:eps*(X*B*YU+2*H*Uz/L+2*(-1+delta)*E*U/L+2*r*s.diff(w,r)*vn+2*za(-1-delta,w)*U
             +2*(-2*z*P+d*Pz-z*X*px)/L+2*axial)}
    unknown=(F,U,K,P,YF,YU);zunknown=(Fz,Uz,Kz,Pz,s.Integer(0),s.Integer(0))
    for label,variables in (('B0',unknown),('B1',zunknown)):
        for row in range(1,7):
            for col,var in enumerate(variables,1):
                expr=s.Integer(0) if var==0 else s.diff(direct[row],var)
                check(matrix[label].get(f'{row},{col}',algebra.lift(0)),expr,2)
    for row in (4,5,6):check(matrix['g'][str(row)],direct[row].subs({var:0 for var in unknown+zunknown if var!=0}),1)
    return dict(independent_direct_paper_operator_and_forcing_rows=count,
        nonconstant_amplitude_ordinary_Z_derivatives_checked=True,pressure_gradient_feedback_checked=True)


def exact_bridge_identities():
    r=s.symbols('rho',positive=True);phi,w,m,pc,wc,chi=(s.Function(n)(r) for n in ('phi','w','m','pc','wc','chi'))
    A=s.diff(pc,r)/pc;ratio=phi/pc
    pr=chi*phi*A;wr=chi*ratio*s.diff(wc,r);mr=(w-m)/r
    replacements={s.diff(phi,r):pr,s.diff(w,r):wr,s.diff(m,r):mr}
    wants=(phi*(s.diff(chi,r)*A+chi*chi*A*A+chi*s.diff(A,r)),
           ratio*(s.diff(chi,r)*s.diff(wc,r)+chi*(chi-1)*A*s.diff(wc,r)+chi*s.diff(wc,r,2)),
           (wr-2*mr)/r)
    for expr,want in zip((pr,wr,mr),wants):assert s.simplify(s.diff(expr,r).subs(replacements)-want)==0
    return dict(exact_prescribed_bridge_and_own_mean_FTC_identities=3)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks=dict(**exact_fixture(c),**exact_bridge_identities());rows=0
    for name,packet in packets.items():
        view=owner.report(packet);entry=owner._packets[id(packet)]
        assert view['actual_own_mean_not_comparison_mean'] and not view['actual_point_moment_history_recovered']
        assert all(not view[key] for key in current.OPEN)
        if view['geometry']['chart']=='first':
            raw=owner.source.actual.packet(view['Z'],view['geometry']['phase'],'first')
            for k in range(3):
                assert overlaps(entry[2]['mean'][0].terms[(0,0,0,0)][k],raw['actual_own_six_moments_axial5']['M'][k])
                assert overlaps(entry[2]['Q'][0].terms[(0,0,0,0)][k],raw['actual_radial_Q_axial4'][k])
                rows+=2
            assert all(key[0]>=0 for label in ('B0','B1') for value in entry[4][label].values() for key in value.terms)
            assert all(key[0]>=-1 for value in entry[4]['g'].values() for key in value.terms)
        rejected(lambda:entry[5].resolve(entry[2]['F'][0]))
    checks['actual_own_mean_and_Q_source_consistency_rows']=rows
    checks['interval_consistency_is_not_exact_point_history_proof']=True
    checks['guards_rejected']={name:rejected(fn) for name,fn in (
        ('core_comparison_extension',lambda:owner.evaluate('.371','4.01','core')),
        ('outside_common_inner_collar',lambda:owner.evaluate('.371','.51','first')),
        ('foreign_chart',lambda:owner.evaluate('.371','0','macro')),
        ('foreign_packet',lambda:owner.report(current.OriginalN1InnerOperatorPacket())),
        ('missing_Z_source_derivative',lambda:current.truncate(owner._packets[id(packets['first_interior'])][3]['N1p'],2)))}
    frame=owner.frame('.371');oldlogs=frame[1].logs;frame[1].logs=(c.mpf(0),)*4
    try:checks['amplitude_frame_mutation_rejected']=rejected(lambda:owner.frame('.371'))
    finally:frame[1].logs=oldlogs
    packet=packets['first_interior'];entry=owner._packets[id(packet)];original=entry[3]['N1p']
    entry[3]['N1p']=original+1
    try:checks['live_math_packet_mutation_rejected']=rejected(lambda:owner.report(packet))
    finally:entry[3]['N1p']=original
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    result=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.snapshot(result),separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print('PASS_ACTUAL_N1_OPERATOR_INPUT_ENCLOSURES',checks['independent_direct_paper_operator_and_forcing_rows'],rows,flush=True)
    return result
