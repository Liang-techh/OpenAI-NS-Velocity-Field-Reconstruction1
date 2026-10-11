"""Independent symbolic source/matrix derivatives and regular axis limits."""
import json
import math
from pathlib import Path
import time

import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_second_radial_function as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=current.op;ep=current.ep


def rejected(fn):
    try:fn()
    except ValueError:return True
    raise AssertionError('Unsupported second radial packet admitted')


def overlaps(left,right,order=1):
    c=left.ctx
    for power in left.terms.keys()|right.terms.keys():
        for k in range(order+1):
            a=left.terms[power][k] if power in left.terms else c.mpf(0)
            b=right.terms[power][k] if power in right.terms else c.mpf(0)
            assert ep(a)[0]<=ep(b)[1] and ep(a)[1]>=ep(b)[0],(power,k)


def equal_coefficients(left,right,order=1):
    """Ignore source sectors whose retained coefficients are all zero."""
    c=left.ctx
    for power in left.terms.keys()|right.terms.keys():
        for k in range(order+1):
            a=left.terms[power][k] if power in left.terms else c.mpf(0)
            b=right.terms[power][k] if power in right.terms else c.mpf(0)
            assert a==b,(power,k)


def rational(c,value):
    value=s.cancel(value)
    return c.mpf(int(s.numer(value)))/int(s.denom(value))


def resolved(value,k):
    c=value.ctx
    return sum((row[k]*math.factorial(k)*c.exp(sum((value.algebra.logs[i]*p[i] for i in range(4)),c.mpf(0)))
                for p,row in value.terms.items()),c.mpf(0))


@source_precision
def independent_source_and_matrix_fixture(c):
    z,r,x=s.symbols('z rho x',positive=True);z0=s.Rational(1,4);sub={z:z0,r:s.Integer(1)}
    dt=s.Rational(1,5);lam=s.Integer(3);Pi=s.Integer(2);eps=1/lam
    d=1-z*z;L=1-dt*z*z;gradient=s.Rational(3,10)+z/5
    phi=1+r*(1+z/2)+r*r*(1+z*z)/7+r**3*(1-z)/11
    psi=r*(1+z)+r*r*z*z/5+r**3/9
    meanpsi=r*(1+z)/2+r*r*z*z/15+r**3/36
    u=4*z+s.Rational(3,10);w=u+eps*psi;wm=u+eps*meanpsi
    q=(2*z*w-(1-dt)*z*wm-d*s.diff(wm,z))/L
    F0=5*s.exp(-lam*(s.Rational(3,10)*(z-z0)+(z*z-z0*z0)/10))
    T=lambda a,f:(a*z*f+d*s.diff(f,z)-2*z*r*s.diff(f,r))/L
    Np=-((q+(1-dt)*z*s.diff(q,z)/2+r*s.diff(q,r))/L+
         q*(q/2+r*s.diff(q,r))+w*(-2*z*q+d*s.diff(q,z)-2*z*r*s.diff(q,r))/L)/2
    Np+=lam*(2*s.diff(q,r)+r*s.diff(q,r,2))
    expected_sources={'4_over_x':2*Np/Pi,
        '5':s.cancel(-2*T(-3,T(-2-dt,F0*phi))/(F0*lam**2)),
        '6':-2*T(-2,T(-1-dt,w))-eps*4*z*r*Np/L}
    algebra=op.N1LogAlgebra(c,(c.ln(c.mpf('.2')),c.ln(5),c.ln(2),c.ln(3)),[])
    jet=lambda f,n:op.Jet(c,[rational(c,s.diff(f,z,k).subs(sub))/math.factorial(k) for k in range(n+1)])
    grid=lambda f,grades:[algebra.lift(jet(s.diff(f,r,i),n)) for i,n in enumerate(grades)]
    rows={'Phi':grid(phi,(5,4,3,2)),'Psi':grid(psi,(5,4,3,2)),
          'meanPsi':grid(meanpsi,(5,4,3,2))}
    fields={'Uz':grid(w,(5,4,3,2)),'Q':grid(q,(4,3,2,1))}
    grad=jet(gradient,5);chi_expr=s.Rational(2,5)+z/10;chi=jet(chi_expr,2)
    sources=current.source_rho_pairs(algebra,rows,fields,c.mpf('.25'),c.mpf(1),c.mpf('.2'),grad)
    count=0
    for key,expr in expected_sources.items():
        for radial,v in ((0,sources[key].v),(1,sources[key].r)):
            for k in (0,1):
                got=resolved(v,k);want=rational(c,s.diff(expr,r,radial,z,k).subs(sub))
                assert ep(got)[0]<=ep(want)[1] and ep(got)[1]>=ep(want)[0],('source',key,radial,k)
                count+=1
    rr=s.sqrt(r);fr=s.diff(phi,r);wr=s.diff(w,r)
    H=z*(1-dt)/2+d*w;E=z*w-s.Rational(1,2);B=q+(1-2*z*w)/L;Cf=phi+r*fr
    Ta=T(-2-dt,phi);Zw=T(-1-dt,w);m=gradient*d/L;gauge=-2*eps*psi*d*gradient/L
    B0={'4,1':F0**2*4*rr*phi,
        '5,1':eps*2*(q+E*(-2+dt)/L)+gauge,
        '5,2':eps**2*(2*Ta-Cf*z*2*(-1+dt)/L)-eps*2*m*phi,
        '5,3':-eps**2*Cf*z*2*(1+dt)/L,'5,5':eps*B*rr,
        '6,1':-eps*F0**2*8*z*r*phi/L,
        '6,2':eps*2*(E*(-1+dt)/L+Zw-wr*z*(-1+dt)*r/L),
        '6,3':-eps*wr*z*2*(1+dt)*r/L,'6,4':-eps*4*z/L,'6,6':eps*B*rr}
    B1={'5,1':eps*2*H/L,'5,2':-eps**2*2*d*Cf/L,'5,3':-eps**2*2*d*Cf/L,
        '6,2':eps*2*(H-d*r*wr)/L,'6,3':-eps*2*d*r*wr/L,'6,4':eps*2*d/L}
    scales=(1,eps,eps,eps*Pi,1,eps)
    expected={}
    for label,entries in (('B0',B0),('B1',B1)):
        expected[label]={}
        for key,expr in entries.items():
            i,j=(int(v)-1 for v in key.split(','));expected[label][key]=expr*scales[j]/scales[i]
    for key,sgn in (('1,5',1),('2,6',1),('3,6',-1)):expected['B0'][key]=s.Integer(sgn)
    expected['B0']['5,1']-=2*chi_expr
    matrices=current.matrix_x_pairs(algebra,rows,fields,c.mpf('.25'),c.mpf(1),c.mpf('.2'),grad,chi)
    for label,entries in expected.items():
        assert set(entries)==set(matrices[label])
        for key,expr in entries.items():
            for radial,v in ((0,matrices[label][key].v),(1,matrices[label][key].r)):
                derivative=s.diff(expr.subs(r,x*x),x,radial).subs(x,1)
                for k in (0,1):
                    got=resolved(v,k);want=rational(c,s.diff(derivative,z,k).subs(z,z0))
                    assert ep(got)[0]<=ep(want)[1] and ep(got)[1]>=ep(want)[0],(label,key,radial,k)
                    count+=1
    return dict(independent_original_source_matrix_value_radial_and_Z1_fixture_rows=count,
                fixture_checks_physical_F0_Z_gauge_and_constant_natural_conjugation=True)


def collar_and_axis_identities():
    r=s.symbols('rho',positive=True);u=s.Function('u')(r);A=s.Function('A')(r)
    ratio=s.Function('ratio')(r);S=s.Function('S')(r)
    C=(u-1)*A;B=s.diff(u,r)+u*(u-1)*A
    psi1=ratio*u*s.diff(S,r)
    psi2=s.expand(s.diff(psi1,r).subs(s.diff(ratio,r),ratio*C))
    psi3=s.expand(s.diff(psi2,r).subs(s.diff(ratio,r),ratio*C))
    want=ratio*(u*s.diff(S,r,3)+2*B*s.diff(S,r,2)+(s.diff(B,r)+C*B)*s.diff(S,r))
    assert s.expand(psi3-want)==0
    P=s.Function('P')(r);phi=s.Function('phi')(r);Cp=u*A
    phi3=phi*(Cp**3+3*Cp*s.diff(Cp,r)+s.diff(Cp,r,2))
    expanded=phi*(u**3*A**3+3*u*s.diff(u,r)*A**2+3*u*u*A*s.diff(A,r)+
                   s.diff(u,r,2)*A+2*s.diff(u,r)*s.diff(A,r)+u*s.diff(A,r,2))
    assert s.expand(phi3-expanded)==0
    hb=s.symbols('hb',positive=True);t=s.symbols('t');sig=s.Function('sigma')
    phase=s.log(r/4)/hb;cut=1-(1-hb)*sig(phase)
    expected=(1-hb)*(hb*s.diff(sig(t),t).subs(t,phase)-s.diff(sig(t),t,2).subs(t,phase))/(hb**2*r*r)
    assert s.simplify(s.diff(cut,r,2)-expected)==0
    x=s.symbols('x',positive=True);a1,b1,a2,b2,a3,b3,a4,b4=s.symbols('a1 b1 a2 b2 a3 b3 a4 b4')
    F=a1*x*x+b1*x**4;U=a2*x*x+b2*x**4;K=a3*x*x+b3*x**4;P=a4*x*x+b4*x**4
    assert s.expand(s.diff(F,x,2)+3*s.diff(F,x)/x).coeff(x,2)==24*b1
    assert s.expand(s.diff(U,x,2)+s.diff(U,x)/x).coeff(x,2)==16*b2
    assert s.expand(s.diff(K,x)+2*K/x+s.diff(U,x)).coeff(x,3)==6*b3+4*b2
    assert s.expand(s.diff(P,x)).coeff(x,3)==4*b4
    assert [s.diff(v,x,2).subs(x,0) for v in (F,U,K,P)]==[2*a1,2*a2,2*a3,2*a4]
    return dict(actual_collar_third_derivative_and_hb_inverse_squared_chain_rule_verified=True,
                genuine_axis_even_rho2_factors_12_8_K_relation_and_pressure_verified=True,
                finite_axis_xx_limits_derived_from_regular_parity=True)


@source_precision
def run(owner,packets):
    began=time.monotonic();owner.assert_graph();c=owner.c
    checks={**independent_source_and_matrix_fixture(c),**collar_and_axis_identities()}
    for name,p in packets.items():
        v=owner.report(p);rows,fields,third,sources,matrices,srho,mx,gx,dxx,drr,recovery2=owner._packets[id(p)][3]
        pp=owner._packets[id(p)][2];W,dx,dr,_,_,system=owner.regular._packets[id(pp)][5]
        assert all(not v[k] for k in current.OPEN)
        assert all(third[k].order==2 for k in ('Phi','Psi','meanPsi')) and fields['Q'][3].order==1
        assert all(row.order==1 for row in dxx+drr)
        assert recovery2['mean1'].order==1 and recovery2['Q1'].order==recovery2['V1'].order==0
        for key in ('5','6'):overlaps(op.truncate(sources[key].v,1),op.truncate(system['g'][key],1))
        overlaps(op.truncate(sources['4_over_x'].v,1),op.truncate(system['g4_over_x'],1))
        for label in ('B0','B1'):
            assert set(matrices[label])==set(system[label])
            for key,pair in matrices[label].items():overlaps(op.truncate(pair.v,1),op.truncate(system[label][key],1))
        if name=='core_axis':
            assert len(drr)==4 and len(v['physical_second_rho_function_jets'])==4
            assert all(not d.terms for d in dxx[4:])
            expected=(system['g']['5']/4,system['g']['6']/2,-system['g']['6']/4,system['g4_over_x'])
            for got,want in zip(dxx,expected):assert op.export(got)==op.export(op.truncate(want,1))
        else:
            assert len(drr)==6
            equal_coefficients(dxx[0],op.truncate(dx[4],1))
            equal_coefficients(dxx[1],op.truncate(dx[5],1))
            for got,xx,vx in zip(drr,dxx,dx):
                assert op.export(got)==op.export(op.truncate((xx-op.truncate(vx,1)/system['x'])/(4*system['rho']),1))
        for field,xx,scale in zip(current.parent.FIELDS,dxx,current.natural.WEIGHTS):
            assert v['physical_second_x_function_jets'][field]==op.export(current.weighted.transform(xx,system['gradient'],scale,1))
    ce=owner.report(packets['core_exit']);fi=owner.report(packets['actual_first_inlet'])
    for key in ('leading_normalized_Phi_Psi_meanPsi_rho3','leading_physical_Q_rho3','normalized_source_rho1',
                'normalized_B0_B1_x1','normalized_source_x1','normalized_second_x_function_jets','normalized_second_rho_function_jets',
                'normalized_own_mean_Q_V_rho2','physical_own_mean_Q_V_rho2'):
        assert ce[key]==fi[key],key
    checks['new_leading_rho3_and_Q_rho3_grades_present_at_all_original_queries']=True
    checks['radial_derivative_source_values_agree_with_original_accepted_operator']=True
    checks['six_equation_second_x_and_rho_mixed_Z1_execute_with_full_function_errors']=True
    checks['core_actual_inlet_third_sources_and_second_coefficient_jets_equal']=True
    checks['physical_second_radial_Z1_adds_F0_gauge_once']=True
    checks['own_mean_rho2_Z1_and_Q_V_rho2_value_recovered_without_new_fit']=True
    entry=owner._packets[id(packets['actual_first_interior'])]
    old=entry[3][3]['5'].r;entry[3][3]['5'].r=old+1
    try:checks['hidden_radial_source_pair_mutation_rejected']=rejected(lambda:owner.report(packets['actual_first_interior']))
    finally:entry[3][3]['5'].r=old
    old=entry[1]['full_axis_Taylor_remainder_derivative_joins_own_repair_n2_still_open']
    entry[1]['full_axis_Taylor_remainder_derivative_joins_own_repair_n2_still_open']=False
    try:checks['false_full_radial_scope_mutation_rejected']=rejected(lambda:owner.report(packets['actual_first_interior']))
    finally:entry[1]['full_axis_Taylor_remainder_derivative_joins_own_repair_n2_still_open']=old
    checks['foreign_packet_and_unsupported_coordinate_rejected']=all((
        rejected(lambda:owner.report(current.OriginalN1SecondRadialPacket())),
        rejected(lambda:owner.evaluate('.371','4.1','core'))))
    hashes=dict(owner.hashes)
    for name in (current.NAME,Path(__file__).name):hashes[name]=current.sha(name)
    receipt=dict(all_passed=True,source_family=owner.family,input_hashes=hashes,checks=checks,
                 **{current.GATE:True},**dict.fromkeys(current.OPEN,False),execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(op.snapshot(receipt),separators=(',',':'))+'\n',encoding='utf8',newline='\n')
    print('PASS_N1_SECOND_RADIAL_FUNCTION',checks['independent_original_source_matrix_value_radial_and_Z1_fixture_rows'],flush=True)
    return receipt
