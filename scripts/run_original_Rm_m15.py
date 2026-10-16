"""Integrate original active-patch m15 to Rh, with full inherited pressure errors."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_Rm_m15 import OriginalRmM15
from openai_ns_reconstruction.original_Rm_moment_kernels import OriginalRmMomentKernels
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_Rm_m15_numeric.json'


def independent_physical_integral():
    """Finite polynomial-bump family; direct S1 integration uses independent operators."""
    c=mp.mp.clone();c.dps=95;iv=MPIntervalContext();iv.dps=120
    z,xx=sy.symbols('Z x');delta=sy.Rational(1,17);unit=sy.Rational(7,100);radius=sy.Integer(3)
    centers=(Fraction(5,4),Fraction(7,4));edges=[sy.Rational(1),sy.Rational(49,40),sy.Rational(51,40),sy.Rational(69,40),sy.Rational(71,40),sy.Integer(2)]
    H=[(1+z+z**3)/20,(2-z+z*z)/30]
    number=lambda value:c.mpf(str(value.p))/int(value.q)
    def sym_cdf(i,mid):
        a,b=sy.Rational(centers[i]-Fraction(1,40)),sy.Rational(centers[i]+Fraction(1,40))
        if mid<=a:return sy.Integer(0)
        if mid>=b:return sy.Integer(1)
        s=(xx-a)/(b-a);return 10*s**3-15*s**4+6*s**5
    with mp.workdps(95):nodes,weights=mp.gauss_quadrature(40,'legendre')
    def quad(fun,a,b):return sum((c.mpf(w)*fun((a+b)/2+(b-a)*c.mpf(node)/2)*(b-a)/2 for node,w in zip(nodes,weights)),c.mpf(0))
    def cdf(i,x):
        a=c.mpf(centers[i].numerator)/centers[i].denominator-c.mpf('.025');s=(x-a)/c.mpf('.05')
        if s<=0:return c.mpf(0)
        if s>=1:return c.mpf(1)
        return 10*s**3-15*s**4+6*s**5
    def beta(i,x):
        a=c.mpf(centers[i].numerator)/centers[i].denominator-c.mpf('.025');s=(x-a)/c.mpf('.05')
        if not 0<s<1:return c.mpf(0),c.mpf(0)
        return 600*s*s*(1-s)**2,24000*s*(1-s)*(1-2*s)
    f=object.__new__(OriginalRmM15);f.ctx,f.iv,f.J=c,iv,J;f._polynomials={}
    f.o=SimpleNamespace(core=SimpleNamespace(delta=c.mpf(1)/17,delta_iv=iv.mpf(1)/17))
    def geo(x):return dict(x_point=c.mpf(x),x=iv.mpf(x))
    f.rows=SimpleNamespace(geometry=geo)
    k=object.__new__(OriginalRmMomentKernels);k.ctx,k.iv,k.audit=c,iv,c;k._packets={};k.rows=f.rows;k.centers=centers
    k.cdf=lambda i,q:(cdf(i,c.mpf(q.numerator)/q.denominator),iv.mpf(str(cdf(i,c.mpf(q.numerator)/q.denominator))))
    k.tail=lambda i,q:(1-k.cdf(i,q)[0],1-k.cdf(i,q)[1])
    k.beta=lambda ctx,i,x,directed:tuple(ctx.mpf(str(v)) for v in beta(i,c.mpf(x)))
    def moments(i,q):
        end=c.mpf(q.numerator)/q.denominator;cuts=sorted({c.mpf(1),end}|{number(v) for v in edges[1:-1] if number(v)<end})
        funcs=dict(A=lambda x:x*beta(i,x)[0],B=lambda x:x*beta(i,x)[0]**2,
            J=lambda x:1-cdf(i,x),L=lambda x:(1-cdf(i,x))/x,S=lambda x:(1-cdf(i,x))**2/x)
        values={name:sum((quad(fun,a,b) for a,b in zip(cuts,cuts[1:])),c.mpf(0)) for name,fun in funcs.items()}
        return dict(point=values,truth={name:iv.mpf(str(v)) for name,v in values.items()},cells=0)
    k.moments=moments;f.kernels=k
    fixture=0
    for Z in ('-.3','.23'):
        zz=sy.Rational(Z)
        def jet(expr,ctx):return J(ctx,[ctx.mpf(str(sy.N(sy.diff(expr,z,n).subs(z,zz),95)))/math.factorial(n) for n in range(4)])
        hp,ht=[jet(h,c) for h in H],[jet(h,iv) for h in H]
        data=dict(normalized_point=hp,normalized_truth=ht,defects_truth=[-ht[0]-ht[1]])
        f.rows.controls=SimpleNamespace(coefficients=lambda Z:data)
        for end in ('1.26','1.6','2'):
            xp=c.mpf(end);kern=f.integral(Z,end)
            up,uf=f.endpoint_F(Z,end,False);low,lf=f.endpoint_F(Z,'1',False)
            direct=[c.mpf(0),c.mpf(0)]
            for left,right in zip(edges,edges[1:]):
                a=number(left);b=min(number(right),xp)
                if b<=a:continue
                I=[sym_cdf(i,(left+right)/2) for i in range(2)];bs=[sy.diff(v,xx) for v in I]
                U=4*z+unit*sum(H[i]*bs[i] for i in range(2))
                M=4*z-unit*sum(H[i]*(1-I[i]) for i in range(2))/xx
                A,D=1-delta*z*z,1-z*z
                Q=(2*z*U-(1-delta)*z*M-D*sy.diff(M,z))/A;Qy=xx*sy.diff(Q,xx)
                flow=((Q+(1-delta)*z*sy.diff(Q,z)/2+Qy)/A+Q*(Q/2+Qy)
                    +U*(-2*z*Q+D*sy.diff(Q,z)-2*z*Qy)/A)
                for n in (0,1):
                    fun=sy.lambdify(xx,sy.diff(radius**2*xx*flow/2,z,n).subs(z,zz),'mpmath',cse=True)
                    with mp.workdps(95):direct[n]+=quad(fun,a,b)
            for n in (0,1):
                e=c.mpf('0.07');want=direct[n]-3*e*(xp*uf['Q_y'][n]-lf['Q_y'][n])
                got=9*sum((e**d*(xp*xp*up[(d,0)][n]-low[(d,0)][n]+kern['point'][d][n]) for d in range(3)),c.mpf(0))
                got-=3*e*(xp*uf['Q_y'][n]-lf['Q_y'][n])
                assert abs(got-want)<c.mpf('1e-65')*max(1,abs(want)),(Z,end,n,str(got-want));fixture+=1
    return dict(passed=True,direct_independent_physical_C0_C1_integrals=fixture,
        finite_Rm=3,finite_t='.07',independent_polynomial_normalized_bumps=True,
        direct_S1_operators_include_Qy=True,original_half_S1_normalization_and_RQy_sign=True)


def run():
    start=time.perf_counter();fixture=independent_physical_integral()
    print('Independent finite physical patch integrals and native shape reductions pass',flush=True)
    f=OriginalRmM15(precision=85,cells=32,partial_cells=128);audit=mp.mp.clone();audit.dps=190;samples=[]
    counts=dict(exact156term_inlet_replays=0,serialized_increment_errors=0,empty_increment_rows=0,
        physical_pressure_source_seams=0,complete_prefix_error_conditions=0,lower_finite_N_endpoint_conditions=0,
        positive_scalar_kernel_conditions=0,exact_terminal_RQy_rows=0,source_mass_branch_diagnostics=0,
        refined_kernel_overlaps=0,strictly_refined_tail_kernel_widths=0)
    def enclosure(row):return tuple(audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples'])
    def loghigh(row):return audit.make_mpf(tuple(row['log_absolute_upper']['upper_exact_mpf_tuple']))
    for Z in ('-.3','0','.3'):
        for x in ('1','1.25','1.5','1.75','1.775','Rh'):
            out=f.evaluate(Z,x);samples.append(out);assert f.evaluate(Z,x) is out and all(out['contract'].values())
            assert len(out['m15_C0_C1_terms'])==168
            assert out['m15_C0_C1_terms'][:156]==out['incoming_m15_terms']==f.inlet(Z)['m15_C0_C1_terms']
            counts['exact156term_inlet_replays']+=1
            for term in out['m15_C0_C1_terms'][156:]:
                for value in term['coefficient_C0_C1']:
                    lo,hi=enclosure(value);p=audit.mpf(value['approximation'])
                    assert max(abs(lo-p),abs(hi-p))<=audit.mpf(value['absolute_error_upper']);counts['serialized_increment_errors']+=1
                    if out['exact_empty_increment']:
                        assert lo==hi==p==0;counts['empty_increment_rows']+=1
                assert term['actual_function_binding']['scale_degree_is_original_fixed_axial_t_not_epsilon']
            _,old,_=f.rows.background(Z,x);v=f.normalized_fields(Z,x,True)
            z=f.J.variable(f.iv,f.iv.mpf(Z),3);t=f.rows.patch.t
            for name,values in (('U',4*z+v['U']*t),('m',4*z+v['m']*t),('Q',v['Q0']+v['Q']*t),('Q_y',v['Q_y']*t)):
                for n in range(values.order+1):
                    a,b=_bounds(audit,values[n]);c,d=_bounds(audit,old[name][n]);assert max(a,c)<=min(b,d),(Z,x,name,n)
                    counts['physical_pressure_source_seams']+=1
            for value in out['exact_mass_row_branch_diagnostic']:
                lo,hi=enclosure(value);assert lo<=0<=hi;counts['source_mass_branch_diagnostics']+=1
            for code,value in out['scalar_kernel_receipt']['kernels'].items():
                lo,hi=enclosure(value);assert lo>=0 and hi>=lo;counts['positive_scalar_kernel_conditions']+=1
            if not out['exact_empty_increment']:
                endpoint=f.rows.evaluate(Z,x)['finite_N_history_error_logs'];prefix=out['actual_lower_finite_N_endpoint_budget']['prefix_error_packet']
                assert all(prefix['contract'].values())
                for name,rows in endpoint.items():
                    for n,row in enumerate(rows):
                        assert loghigh(prefix['finite_N_prefix_history_logs'][name][n])>=loghigh(row),(Z,x,name,n)
                        counts['complete_prefix_error_conditions']+=1
                lower=out['actual_lower_finite_N_endpoint_budget']
                for kind in ('endpoint','viscous'):
                    for row in lower[kind]:assert not row['exact_zero'];counts['lower_finite_N_endpoint_conditions']+=1
            if x in ('1.775','Rh'):
                term=next(t for t in out['m15_C0_C1_terms'] if t['name']=='Rm_patch_RQy_upper_t1')
                for row in term['coefficient_C0_C1']:
                    assert enclosure(row)==(0,0);counts['exact_terminal_RQy_rows']+=1
        print('Original Rm cumulative signed m15: Z='+Z+',156 ->168 terms',flush=True)
    refined=OriginalRmMomentKernels(f.rows,cells=64,partial_cells=256)
    for x in ('1.25','Rh'):
        low=f.kernels.evaluate(x);high=refined.evaluate(x)
        for code,box in high['truth'].items():
            a,b=_bounds(audit,low['truth'][code]);c,d=_bounds(audit,box);assert max(a,c)<=min(b,d),(x,code)
            counts['refined_kernel_overlaps']+=1
            if code in (18,54,36) and d-c<b-a:counts['strictly_refined_tail_kernel_widths']+=1
    assert counts['strictly_refined_tail_kernel_widths']>0
    rejected=0
    for Z,x in (('1','1'),('0','0'),('0','3')):
        try:f.evaluate(Z,x)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid patch moment query admitted')
    code=dict(f.code_hashes);name='scripts/run_original_Rm_m15.py';code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_physical_integral=fixture,actual_original_packets=len(samples),
            **counts,rejected_invalid_queries=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(Rm_Rh_cumulative_signed_m15_with_complete_absolute_errors=True,
            resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
            independent_preRp_signed_seam=False,global_five_moment_numeric_join=False,common_global_N=False,
            exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    data=json.loads(path.read_text(encoding='utf-8'));assert data['checks']['passed']
    for base,group in ((DIRECTORY,data['input_hashes']),(ROOT,data['implementation_hashes'])):
        for name,digest in group.items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;data=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
