"""Evaluate original changing O2 axial pressure and cumulative m15."""
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
from openai_ns_reconstruction.original_axial_m15 import OriginalAxialM15
from openai_ns_reconstruction.original_axial_pressure_rows import OriginalAxialPressureRows
from openai_ns_reconstruction.original_axial_rate_memory import OriginalAxialRateMemory
from openai_ns_reconstruction.original_Rm_finite_N_numeric import RATES
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_axial_m15_numeric.json'


def independent_fixtures():
    c=mp.mp.clone();c.dps=100;iv=MPIntervalContext();iv.dps=130
    z,de,b,m=sy.symbols('z de b m');A,D=1-de*z*z,1-z*z;U=4*z*b;Q=8*z*z*b/A-4*m
    F=Q/(2*A)+Q*Q/4-z*z*U*U/A**2
    G=-Q/(2*A)+z*(1-de)*sy.diff(Q,z)/(4*A)-Q*Q/4+D*U*sy.diff(Q,z)/(2*A)+z*z*(1-de)*U*U/A**2+z*D*U*sy.diff(U,z)/A**2
    f=object.__new__(OriginalAxialM15);f.ctx,f.iv,f.J=c,iv,J
    f.o=SimpleNamespace(core=SimpleNamespace(delta=c.mpf(1)/17,delta_iv=iv.mpf(1)/17))
    coefficients=0
    for Z in ('-.3','0','.23'):
        fp,gp,_=f.coefficients(Z,c)
        for poly,actual in ((F,fp),(G,gp)):
            for mode,powers in (('b',(1,0)),('bm',(1,1)),('m2',(0,2)),('m',(0,1)),('b2',(2,0))):
                target=sy.Poly(sy.expand(poly),b,m).coeff_monomial(b**powers[0]*m**powers[1])
                for n in (0,1):
                    want=c.mpf(str(sy.N(sy.diff(target,z,n).subs({z:sy.Rational(Z),de:sy.Rational(1,17)}),105)))
                    got=actual.get(mode,J.constant(c,0,1))[n]
                    assert abs(got-want)<c.mpf('1e-90');coefficients+=1
    # Independent direct physical quadrature with a nonconstant polynomial
    # cutoff having b(0)=1,b_y(0)=0 and its exact relaxation mean.
    y,k=sy.symbols('y k');bb=1-k*y*y;mm=1-k*(y*y-2*y+2-2*sy.exp(-y))
    FF,GG,QQ=(v.subs({b:bb,m:mm}) for v in (F,G,Q))
    density=sy.exp(2*y)*(2*FF+sy.diff(FF,y)+GG)-sy.exp(y)*(sy.diff(QQ,y)+sy.diff(QQ,y,2))
    physical=0
    for Z in ('-.3','0','.23'):
        fp,gp,vp=f.coefficients(Z,c)
        for endpoint in ('.07','.8'):
            L=c.mpf(endpoint);kk=c.mpf('.2');beta=lambda x:1-kk*x*x
            mean=lambda x:1-kk*(x*x-2*x+2-2*c.exp(-x))
            B,M=beta(L),mean(L);by=-2*kk*L
            ib=c.quad(lambda x:c.exp(2*x)*beta(x),[0,L])
            ib2=c.quad(lambda x:c.exp(2*x)*beta(x)**2,[0,L])
            im2=c.quad(lambda x:c.exp(2*x)*mean(x)**2,[0,L])
            ibm=(c.exp(2*L)*M*M-1)/2;im=c.exp(2*L)*M-1-ib
            assembled=(fp['b']*B+fp['bm']*B*M+fp['m2']*M*M+fp['m']*M)*c.exp(2*L)-sum(fp.values(),J.constant(c,0,1))
            assembled+=gp['b2']*ib2+gp['bm']*ibm+gp['m2']*im2+gp['b']*ib+gp['m']*im
            assembled+=(vp['by']*by+vp['b']*B+vp['m']*M)*c.exp(L)
            for n in (0,1):
                expr=sy.diff(density,z,n).subs({z:sy.Rational(Z),de:sy.Rational(1,17),k:sy.Rational(1,5)})
                direct=sy.lambdify(y,expr,modules=[{'exp':c.exp,'mpf':c.mpf},'mpmath'])
                want=c.quad(direct,[0,L]);assert abs(assembled[n]-want)<c.mpf('1e-85'),(Z,endpoint,n,c.nstr(assembled[n]-want,12));physical+=1
    # Exercise the production CDF square primitive with a positive cubic
    # source over three cells, against a separate nested direct quadrature.
    engine=OriginalAxialRateMemory(c,iv,J);L=c.mpf('.7');poly=lambda x:c.mpf('.8')+x*c.mpf('.03')+x*x*c.mpf('.01')+x**3*c.mpf('.002')
    cells=[]
    for a,bv in zip(map(c.mpf,('0','.13','.4')),map(c.mpf,('.13','.4','.7'))):
        mid=(a+bv)/2;co=[c.diff(poly,mid,n)/math.factorial(n) for n in range(4)]
        mass=c.quad(lambda x:c.exp(-x)*poly(x),[a,bv]);radius=c.mpf('1e-95')
        cells.append(dict(left=a,right=bv,coefficients=[iv.mpf(v) for v in co],point_coefficients=co,
            fourth_absolute_cap=iv.mpf(0),point=mass,truth=iv.mpf([mass-radius,mass+radius])))
    distribution=dict(ordinary=True,anchor=c.mpf(0),scale=c.mpf(1),q_origin=c.mpf(0),left=c.mpf(0),right=L,
        actual_distance=iv.mpf(L),point_distance=L,source_cells=cells,
        normalized_source_integral=sum((row['truth'] for row in cells),iv.mpf(0)),
        normalized_homogeneous=iv.exp(-iv.mpf(L)),normalized_point_homogeneous=c.exp(-L),
        normalized_left_tail=iv.mpf(0),normalized_right_tail=iv.mpf(0),normalized_endpoint_mass_error=iv.mpf(0),
        lower_endpoint_clipped=True,upper_endpoint_clipped=True)
    engine.evaluate=lambda *args:dict(distribution=distribution)
    square=engine.squared_mean_integral('.25')
    # Integrate the analytic source primitive independently with tanh-sinh.
    tail=lambda q:c.exp(-L)+c.quad(lambda x:c.exp(-x)*poly(x),[q,L])
    want=c.quad(lambda q:tail(q)**2,[0,L]);lo,hi=_bounds(c,square['truth'])
    assert lo<=want<=hi and abs(square['point']-want)<c.mpf('1e-85')
    # Native s cells are integrated in physical d, with all four ordinary rows.
    r=object.__new__(OriginalAxialPressureRows);r.ctx,r.iv,r.J=c,iv,J;r.chart='O2_axial';r._histories={};r.owner=SimpleNamespace(logN=iv.log(17))
    incoming={name:[iv.log(iv.mpf((i+1)*(n+1))/19) for n in range(4)] for i,name in enumerate(RATES)}
    r.incoming_history_logs=lambda Z:incoming
    intervals=[(Fraction(0),Fraction(1,100)),(Fraction(1,100),Fraction(1,50)),(Fraction(1,50),Fraction(1))]
    groups=[]
    for j,(a,bv) in enumerate(intervals):
        cell=dict(actual_density_ranges=dict(N_scaled={name:[iv.log(iv.mpf((j+1)*(i+1)*(n+1))/7) for n in range(4)] for i,name in enumerate(RATES)}))
        groups.append((a,bv,[cell,cell]))
    r.utility=SimpleNamespace(_cover=lambda Z,q:([['-1','0'],['0','1']],groups,[groups[1][2][0]]))
    r.cover_owner=SimpleNamespace(_log_cover=lambda values:iv.mpf(max(_bounds(c,v)[1] for v in values)))
    end=c.expm1(c.mpf('.6'));out=r.history_logs('0','.015');transport=0
    for i,(name,rate) in enumerate(RATES.items()):
        rr=c.mpf(rate)
        for n in range(4):
            expected=c.mpf((i+1)*(n+1))/19*c.exp(-rr*end)
            for j,(a,bv) in enumerate(intervals):
                left=c.expm1(40*c.mpf(a.numerator)/a.denominator)
                right=min(end,c.expm1(40*c.mpf(bv.numerator)/bv.denominator))
                if right<=left:continue
                driver=c.mpf((j+1)*(i+1)*(n+1))/(7*17)
                expected+=c.quad(lambda x:driver*c.exp(-rr*(end-x)),[left,right])
            bound=c.exp(_bounds(c,out['history'][name][n])[1]);sup=c.exp(_bounds(c,out['prefix_history'][name][n])[1])
            assert bound>=expected and sup>=bound;transport+=1
    return dict(passed=True,independent_F_G_ordinary_coefficients=coefficients,
        direct_changing_source_physical_C0_C1_integrals=physical,independent_cubic_future_CDF_square=1,
        physical_exponential_width_C3_Duhamel_rows=transport,finite_fixture_N=17)


def run():
    start=time.perf_counter();fixture=independent_fixtures()
    print('Independent axial pressure/CDF/physical-width fixtures pass',flush=True)
    f=OriginalAxialM15(precision=85);a=mp.mp.clone();a.dps=190;samples=[];pressure=[]
    counts=dict(exact_inlet_replays=0,serialized_increment_rows=0,exact_empty_rows=0,pressure_rows=0,
        C3_history_rows=0,prefix_error_conditions=0,lower_endpoint_conditions=0,native_width_conditions=0,
        original_archive_overlaps=0,terminal_nonzero_Qy_conditions=0,positive_kernel_conditions=0)
    def bounds(row):return tuple(a.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples'])
    def highlog(row):return a.make_mpf(tuple(row['log_absolute_upper']['upper_exact_mpf_tuple']))
    path=DIRECTORY/(PREFIX+'compliant_current_original_outer_axial_numeric.json');raw=path.read_bytes();old=json.loads(raw)
    assert old['checks']['passed'] and old['source_family']==f.family;f.hashes[path.name]=hashlib.sha256(raw).hexdigest()
    archive={(p['Z'],p['native_axial']):p for p in old['samples']}
    for Z in ('-.3','0','.3'):
        for value in ('0','.001','.5','.99','.99999','1'):
            out=f.evaluate(Z,value);p=f.rows.evaluate(Z,value);samples.append(out);pressure.append(p)
            assert len(out['m15_C0_C1_terms'])==195
            assert out['m15_C0_C1_terms'][:176]==out['incoming_m15_terms']==f.inlet(Z)['m15_C0_C1_terms'];counts['exact_inlet_replays']+=1
            assert all(out['contract'].values()) and all(p['contract'].values())
            for term in out['m15_C0_C1_terms'][176:]:
                for row in term['coefficient_C0_C1']:
                    lo,hi=bounds(row);point=a.mpf(row['approximation'])
                    assert max(abs(lo-point),abs(hi-point))<=a.mpf(row['absolute_error_upper']);counts['serialized_increment_rows']+=1
                    if value=='0':assert lo==hi==point==0;counts['exact_empty_rows']+=1
            pp,tt,geo=f.rows.background(Z,value)
            for name,field in p['physical_pressure_rows'].items():
                for n,row in enumerate(field['signed_background_ordinary_Z']):
                    assert row['source_enclosure_exact_mpf_tuples']==[list(v) for v in (tt[name][n]*math.factorial(n))._mpi_]
                    counts['pressure_rows']+=1
            assert bounds(p['physical_pressure_rows']['m']['signed_background_ordinary_Z'][3])==(0,0)
            for name,rows in p['finite_N_history_error_logs'].items():
                assert len(rows)==4 and all(not row['exact_zero'] for row in rows);counts['C3_history_rows']+=4
            prefix,_=f.rows.error_data(Z,value,prefix=True);end,_=f.rows.error_data(Z,value)
            for name in ('U','m','Q','U_y','Q_y'):
                for v,w in zip(f.rows.log_rows(prefix[name]),f.rows.log_rows(end[name])):
                    assert highlog(v)>=highlog(w);counts['prefix_error_conditions']+=1
            history=f.rows.history_logs(Z,value);width=sum((step['physical_log_width'] for step in history['steps']),f.iv.mpf(0))
            lo,hi=_bounds(a,width);want=a.expm1(40*a.mpf(value));assert lo<=want<=hi;counts['native_width_conditions']+=1
            if value=='0':assert p['finite_N_history_error_logs']==p['finite_N_incoming_history_error_logs']
            else:
                budget=out['actual_lower_finite_N_endpoint_budget']
                for kind in ('endpoint','viscous'):
                    assert all(not row['exact_zero'] for row in budget[kind]);counts['lower_endpoint_conditions']+=2
                kernels=f.kernels(value)
                for mode in ('beta','beta2'):
                    dist=kernels[mode]['distribution'];lo,hi=_bounds(a,dist['normalized_source_integral']);assert lo>0 and hi>=lo
                    counts['positive_kernel_conditions']+=1
                lo,hi=_bounds(a,kernels['mean2']['truth']);assert lo>0 and hi>=lo;counts['positive_kernel_conditions']+=1
            if value=='1':
                mean=f.rows.scalars(value)['truth']['mean'];lo,hi=_bounds(a,tt['Q_y'][0]-4*mean);assert lo<=0<=hi
                assert bounds(p['physical_pressure_rows']['Q_y']['signed_background_ordinary_Z'][0])[0]>0
                assert bounds(p['physical_pressure_rows']['U']['signed_background_ordinary_Z'][0])==(0,0)
                counts['terminal_nonzero_Qy_conditions']+=1
            oldp=archive.get((Z,value))
            if oldp:
                for name,oldname in (('U','Uz'),('m','Mz_over_R'),('Q','Ur_over_sqrt_R_over2')):
                    for n,row in enumerate(oldp['physical_similarity'][oldname]):
                        lo,hi=bounds(row);ll,hh=_bounds(a,tt[name][n]);ll*=math.factorial(n);hh*=math.factorial(n)
                        assert max(lo,ll)<=min(hi,hh),(Z,value,name,n);counts['original_archive_overlaps']+=1
            print('Original O2_axial pressure/m15: Z='+Z+', s='+value+', 176 ->195 terms',flush=True)
    rejected=0
    for Z,value in (('1','0'),('0','-.1'),('0','1.1')):
        try:f.evaluate(Z,value)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid axial pressure query admitted')
    name='scripts/run_original_axial_m15.py';f.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,pressure_row_samples=pressure,input_hashes=f.hashes,
        implementation_hashes=f.code_hashes,checks=dict(passed=True,independent_fixtures=fixture,
            actual_moment_packets=len(samples),actual_pressure_packets=len(pressure),**counts,rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start,3),limits=copy_limits(samples[-1]['limits']))


def copy_limits(values):return dict(values)


def inspect(path):
    data=json.loads(path.read_text(encoding='utf-8'));assert data['checks']['passed']
    for base,group in ((DIRECTORY,data['input_hashes']),(ROOT,data['implementation_hashes'])):
        for name,digest in group.items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;data=inspect(path) if args.inspect_report else run()
    if args.write_report:
        temp=path.with_suffix('.json.tmp');temp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');temp.replace(path)
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
