"""Focused weighted-density, source-inlet, local-history and m15 unit checks."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace,FunctionType

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_native_pulse_m15 import OriginalNativePulseM15Numeric
from openai_ns_reconstruction.original_native_pulse_pressure import OriginalNativePulsePressureNumeric
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from lei_ren_part1_paper_compliant_axial_pulse_field import gp,decay_integral,raw_beta
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_native_pulse_m15.json'


def independent_units():
    """Execute both actual unit recipes symbolically, with the original p=2."""
    from lei_ren_part1_paper_compliant_current_original_n1_five_defect_values import pressure_unit_log
    lp,lu,lc,xi=sy.symbols('logP logU0 logRc xi',real=True)
    lm=-4*lp-sy.log(1000);mu=sy.exp(lm);rp=lc-60*lm-2
    ctx=SimpleNamespace(mpf=lambda x:sy.Rational(x) if isinstance(x,(str,int)) else x,
                        log=sy.log,ln=sy.log,exp=sy.exp,sqrt=sy.sqrt,pi=sy.pi)
    n=SimpleNamespace(amplitude=SimpleNamespace(lp_iv=lp),
        quiet=SimpleNamespace(numeric=SimpleNamespace(outer=SimpleNamespace(logmui=lm))),
        constants=dict(U=sy.exp(lu)),_row=lambda x:x,_rows=lambda v:list(v.coefficients))
    f=object.__new__(OriginalNativePulseM15Numeric);f.iv,f.native,f.mu=ctx,n,mu
    f.J=lambda iv,rows:SimpleNamespace(coefficients=rows);f.restore=lambda x:x
    p=dict(coordinate=dict(Rp_global_y_log_terms=[rp-sy.log(110)]))
    old=FunctionType(pressure_unit_log.__code__,dict(pressure_unit_log.__globals__,interval=lambda ctx,x:x))
    original=old(ctx,dict(logP=lp,logUp=lu,logRc=lc),2)
    saddle=sy.exp(-sy.Rational(2,3)*(sy.log(4)+lm));finite=sy.Rational(2,3)*sy.log(4)-sy.log(6)/2
    logE=-1/mu-3*saddle-sy.Rational(4,3)*lm+finite;T=13/mu
    logRv=rp+T;logF=lp+logE-(sy.Rational(1,2)+mu)*T
    identities=[sum(f.units(p,'fixed_inlet'))-original,
        sum(f.units(p,'end_linear'))-(2*logRv+logF),
        sum(f.units(p,'end_quadratic'))-(2*logRv+2*logF),
        sum(f.units(p,'end_viscous'))-(logRv+logF),
        sum(f.units(p,'end_linear_gap_query',((xi+sy.Rational(11,2))/mu,)))-(2*logRv+logF+(xi-13)/mu)]
    for value in identities:assert sy.simplify(sy.expand_log(value,force=True))==0
    return len(identities)


def independent_weighted_density():
    """Directly integrate -R*P1_y on main, +R*P1_y on remaining end.

    The reference uses the original Q_y/Q_yy density, not S1 endpoint
    reduction or the new cell matrices. Finite mu=.02 exposes I1 rates.
    """
    c=mp.mp.clone();c.dps=48;iv=MPIntervalContext();iv.dps=85
    e=object.__new__(OriginalNativePulsePressureNumeric)
    e.ctx,e.iv,e.J,e.cells=c,iv,J,128
    e.mu,e.delta,e.normal,e.ell=iv.mpf('.02'),iv.mpf('.001'),iv.mpf('.7'),iv.mpf('.15')
    e.a,e.b=iv.mpf('.48'),iv.mpf('.52');e._kernels,e._bumps,e._quadratics={},{},{}
    e.gp,e.raw_beta=gp,raw_beta
    def data(Z):
        z=J.variable(iv,iv.mpf(Z),4)
        ap=1+z*iv.mpf('.2')+z*z*iv.mpf('.05')
        C=[z*iv.mpf('.1')-iv.mpf('.7'),z*iv.mpf('.2')+z*z*iv.mpf('.1')+iv.mpf('.5')]
        return dict(selected_ap_Taylor=ap,selected_scaled_end_coefficient_Taylor=C),None,(1+z*z).reciprocal()*iv.mpf('.9'),None,None
    n=SimpleNamespace(_data=data,decay_integral=decay_integral,constants=dict(U=iv.mpf('.9'),C1=iv.mpf('.03')))
    e.native=n
    m=object.__new__(OriginalNativePulseM15Numeric)
    m.pressure,m.native,m.ctx,m.iv,m.J=e,n,c,iv,J
    m.mu,m.delta,m.a0,m.a1,m.c,m.cells=e.mu,e.delta,e.a,iv.mpf('1.48'),iv.mpf('.96'),128
    m._k1,m._main_matrix,m._end_matrix={},{},{}
    m.inlet=lambda Z,packet:[]
    m.term=lambda packet,name,value,kind,extra=():dict(jet=value,kind=kind,extra=extra)
    rules={}
    def quad(f,l,r,order):
        if r<=l:return c.mpf(0)
        if order not in rules:rules[order]=c.gauss_quadrature(order,'legendre')
        xs,ws=rules[order];h=(r-l)/2;mid=(r+l)/2
        return h*sum((ws[k]*f(mid+h*xs[k]) for k in range(order)),c.mpf(0))
    def sigma(x):
        if x<=0:return c.mpf(0)
        if x>=1:return c.mpf(1)
        return 1/(1+c.exp(1/x**2-1/(1-x)**2))
    def sigma_prime(x):
        if x<=0 or x>=1:return c.mpf(0)
        s=sigma(x);return s*(1-s)*(2/x**3+2/(1-x)**3)
    def main_shape(x,order):
        primitive=x-c.mpf('.01') if x>=c.mpf('.02') else quad(lambda v:sigma(50*v),0,x,order)
        return primitive,c.mpf('.02')*sigma(50*x),c.mpf('.02')**2*50*sigma_prime(50*x)
    def beta(v,center):
        ell=c.mpf('.15');r=(v-center)/ell
        if abs(r)>=1:return (c.mpf(0),)*3
        f=c.exp(-1/(1-r*r))/(ell*c.mpf('.7'));d=-2*r/(ell*(1-r*r)**2)
        dd=-2/(ell**2*(1-r*r)**2)-8*r*r/(ell**2*(1-r*r)**3)
        return f,f*d,f*(d*d+dd)
    def density(z,v,rows,mean_normal,part=None):
        de=c.mpf('.001');theta=lambda w:c.mpf('.9')/(1+w*w)*c.exp(-c.mpf('.52')*v)
        U=lambda w:theta(w)*rows[0](w)
        b=c.mpf('.52')
        Uy=lambda w:theta(w)*(rows[1](w)-b*rows[0](w))
        Uyy=lambda w:theta(w)*(rows[2](w)-2*b*rows[1](w)+b*b*rows[0](w))
        mean=lambda w:theta(w)*mean_normal(w)
        op=lambda f,p,w:(2*w*f(w)-(1-de)*w*p(w)-(1-w*w)*c.diff(p,w))/(1-de*w*w)
        Q=lambda w:op(U,mean,w)
        Qy=lambda w:op(Uy,lambda q:U(q)-mean(q),w)
        Qyy=lambda w:op(Uyy,lambda q:Uy(q)-U(q)+mean(q),w)
        q,qy,qz=Q(z),Qy(z),c.diff(Q,z);L=1-de*z*z
        linear=(q+(1-de)*z*qz/2+qy)/L
        quadratic=q*(q/2+qy)+U(z)*(-2*z*q+(1-z*z)*qz-2*z*qy)/L
        pieces=[-c.exp(v)*linear/2,-c.exp(v)*quadratic/2,qy+Qyy(z)]
        return sum(pieces) if part is None else pieces[part]
    def reference(chart,x,Z,order,part=None):
        if chart=='end':
            lo,hi=x,c.mpf(0)
            edges=[c.mpf(center)+c.mpf('.15')*j/2 for center in (-3,-1) for j in range(-2,3)]
        else:
            lo,hi=c.mpf(0),x;edges=[c.mpf('.5'),c.mpf(1)]
        cuts=sorted(set([lo,hi]+[v for v in edges if lo<v<hi]));totals=[c.mpf(0),c.mpf(0)]
        if order not in rules:rules[order]=c.gauss_quadrature(order,'legendre')
        xs,ws=rules[order]
        for left,right in zip(cuts,cuts[1:]):
            for k in range(order):
                v=(left+right)/2+(right-left)*xs[k]/2;weight=(right-left)*ws[k]/2
                if chart=='end':
                    beta_rows=[beta(v,j) for j in (-3,-1)];H=[]
                    for center in (-3,-1):
                        a,b=max(v,c.mpf(center)-c.mpf('.15')),c.mpf(center)+c.mpf('.15')
                        pieces=sorted(set([a,b]+[c.mpf(center)+c.mpf('.15')*j/2 for j in (-1,0,1) if a<c.mpf(center)+c.mpf('.15')*j/2<b]))
                        H.append(-sum((quad(lambda w:c.exp(c.mpf('.48')*(w-v))*beta(w,center)[0],l,r,order) for l,r in zip(pieces,pieces[1:])),c.mpf(0)) if a<b else c.mpf(0))
                    C=lambda w:(c.mpf('-.7')+c.mpf('.1')*w,c.mpf('.5')+c.mpf('.2')*w+c.mpf('.1')*w*w)
                    rows=[lambda w,j=j:sum(C(w)[i]*beta_rows[i][j] for i in range(2)) for j in range(3)]
                    mean=lambda w:sum(C(w)[i]*H[i] for i in range(2))
                    sign=1
                else:
                    shapes=main_shape(c.mpf('.02')*v,order)
                    pieces=sorted(set([c.mpf(0),v]+[q for q in (c.mpf('.5'),c.mpf(1)) if 0<q<v]))
                    flow=sum((quad(lambda w:c.exp(-c.mpf('.48')*(v-w))*sigma(w),l,r,order) for l,r in zip(pieces,pieces[1:])),c.mpf(0))
                    k0=shapes[0]/c.mpf('.48')-c.mpf('.02')*flow/c.mpf('.48')
                    ap=lambda w:1+c.mpf('.2')*w+c.mpf('.05')*w*w
                    rows=[lambda w,j=j:ap(w)*shapes[j] for j in range(3)]
                    mean=lambda w:c.mpf('.03')*w*(1+w*w)*c.exp(-c.mpf('.48')*v)+ap(w)*k0
                    sign=-1
                value=lambda z:sign*c.exp(v)*density(z,v,rows,mean,part)
                totals[0]+=weight*value(Z);totals[1]+=weight*c.diff(value,Z)
        return totals
    checked=0;convergence=[]
    with mp.workdps(115):
        cases=[('end','-4'),('end','-3'),('end','-1'),('main','1'),('main','25')]
        for chart,text in cases:
            if chart=='end':terms=m.end('.2',text,{})[0]
            else:terms=m.main('.2','entrance',text,{})
            got=J.constant(iv,0,1)
            for term in terms:
                exponent=sum(term['extra'],iv.mpf(0));got+=term['jet']*iv.exp(exponent)
            coarse,fine=reference(chart,c.mpf(text),c.mpf('.2'),64),reference(chart,c.mpf(text),c.mpf('.2'),96)
            for row in (0,1):
                difference=abs(coarse[row]-fine[row])/(1+abs(fine[row]));convergence.append(c.nstr(difference,12))
                assert difference<c.mpf('2e-8'),(chart,text,row,'reference convergence',difference)
                lo,hi=_bounds(c,got[row])
                if not lo<=fine[row]<=hi:
                    detail=[(j,_bounds(c,terms[j]['jet'][row]),reference(chart,c.mpf(text),c.mpf('.2'),64,j)[row]) for j in range(3)] if chart=='end' else None
                    raise AssertionError((chart,text,row,lo,fine[row],hi,detail))
                checked+=1
        assert _bounds(c,m.end_matrix('-4')[2][3])[0]>0
        assert _bounds(c,m.end_matrix('-2')[3][3])[0]>0
    return dict(defining_weighted_N1p_C0_C1_comparisons=checked,reference_Gauss_orders=[64,96],
                reference_relative_convergence=convergence,nonzero_off_support_history_products=2)


def run():
    start=time.perf_counter();unit_identities=independent_units();independent=independent_weighted_density()
    print('Independent weighted N1p C0/C1, I1 rates and off-support histories pass',flush=True)
    f=OriginalNativePulseM15Numeric();samples=[];source=seams=history=zeros=units=rejected=0
    c=mp.mp.clone();c.dps=250
    def pair(row):return tuple(c.make_mpf(tuple(x)) for x in row['source_enclosure_exact_mpf_tuples'])
    points=[('entrance','0'),('entrance','1'),('main','1'),('main','10.5'),('main','11'),('gap','11'),('gap','12'),
            ('gap_end','-1000'),('gap_end','-4'),('end','-4'),('end','-3'),('end','-2'),('end','-1'),('end','-.5'),('end','0')]
    for Z in ('-.3','0','.3'):
        for chart,x in points:
            p=f.evaluate(Z,chart,x);samples.append(p)
            assert all(p['contract'].values());source+=len(p['contract'])
            assert p['selected_integer']==f.native.inlets[Z]['selected_integer'];source+=1
            assert not p['limits']['independent_preRp_numeric_seam'];source+=1
            for term in p['m15_C0_C1_terms']:
                assert len(term['coefficient_C0_C1'])==2;source+=1
            if chart=='entrance' and x=='0':
                rows=f.complete['whole_Z']['normalized_true_pulse_C1_enclosures']['m15']
                for row,original in zip(p['m15_C0_C1_terms'][0]['coefficient_C0_C1'],rows):
                    a,b=pair(row);ll,hh=_bounds(c,f.native._decode(original));assert max(a,-hh)<=min(b,-ll);seams+=1
            if chart in ('entrance','main') and x not in ('0','11'):
                t=next(t for t in p['m15_C0_C1_terms'] if t['name']=='main_weighted_original_mean_log')
                assert pair(t['coefficient_C0_C1'][1])==(0,0);history+=1
                assert pair(t['coefficient_C0_C1'][0])[1]<0;history+=1
            if chart=='end' and x in ('-.5','0'):
                for term in p['m15_C0_C1_terms']:
                    for row in term['coefficient_C0_C1']:assert pair(row)==(0,0);zeros+=1
            if chart=='gap_end' and x=='-1000':
                term=next(t for t in p['m15_C0_C1_terms'] if t['name']=='gap_weighted_linear_at_query')
                assert pair(term['log_scale_terms'][-1])==(-1000,-1000);units+=1
            if chart=='gap':
                ts=p['m15_C0_C1_terms'];k=len(p['coordinate']['Rp_global_y_log_terms'])+2
                a,b=pair(next(t for t in ts if t['name']=='gap_weighted_linear_at_query')['log_scale_terms'][k])
                ll,hh=pair(next(t for t in ts if t['name']=='gap_weighted_linear_at_minus4')['log_scale_terms'][k])
                num=33 if x=='11' else 35;assert max(37*a,num*ll)<=min(37*b,num*hh);units+=1
        get=lambda chart,x:next(p for p in samples if p['Z']==Z and p['native_pulse_chart']==chart and p['native_pulse_coordinate']==x)
        assert get('main','11')['m15_C0_C1_terms']==get('gap','11')['m15_C0_C1_terms'];seams+=1
        assert get('gap_end','-4')['m15_C0_C1_terms']==get('end','-4')['m15_C0_C1_terms'];seams+=1
        print('Actual repaired weighted m15 reaches Rv: Z='+Z,flush=True)
    for Z,chart,x in [('1','end','0'),('nan','main','1'),('.3','main','12'),('.3','end','1'),('.3','unknown','0')]:
        try:f.evaluate(Z,chart,x)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid m15 query admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent=independent,original_source_unit_identities=unit_identities,actual_weighted_m15_packets=len(samples),source_C1_and_scope_conditions=source,
            terminal_derived_inlet_and_actual_chart_seams=seams,original_nonzero_mean_log_conditions=history,
            terminal_and_empty_support_zeros=zeros,weighted_atomic_unit_conditions=units,rejected_requests=rejected),
        execution_seconds=round(time.perf_counter()-start,3),limits=dict(weighted_m15_C0_C1=True,
            actual_terminal_derived_Rp_inlet=True,independent_preRp_numeric_seam=False,resolved_signed_point_centers=False,
            global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,stress_cone=False,
            genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-report',action='store_true');p.add_argument('--inspect-report',action='store_true')
    args=p.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
