"""Directed original O2 slope profiles and actual selected-N physical values."""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import scipy
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_slope_scalar import OriginalSlopeScalar
from openai_ns_reconstruction.original_outer_slope_numeric import OriginalOuterSlopeVelocity
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from run_original_Rm_patch_velocity import inlet_physical_seam
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_outer_slope_numeric.json'


def exact_forcing_checks():
    c=mp.mp.clone();c.dps=65;iv=MPIntervalContext();iv.dps=85;count=0
    for sig in (0,1):
        s=OriginalSlopeScalar(c,iv,cells=64)
        s.sigma_jets=lambda ctx,x:IntervalTaylor.constant(ctx,sig,4)
        got=s.evaluate('.5');x=c.mpf('.5');a=c.mpf('.1')-c.mpf('.6')*sig
        F=c.exp(a*x);H=c.mpf(5)/8*c.exp(-c.mpf('1.5')*x)+(F-c.exp(-c.mpf('1.5')*x))/(c.mpf('1.5')+a)
        G=(c.mpf(5)/6+x)*c.exp(-x) if sig else c.mpf(5)/6*c.exp(-x)+(c.exp(2*a*x)-c.exp(-x))/(1+2*a)
        C=5+c.expm1(2*a*x)/(2*a)
        for name,exact in zip(s.names,(F,H,G,C)):
            lo,hi=_bounds(c,got['directed'][name]);assert lo<=exact<=hi,(sig,name);count+=1
    return count


def independent_scalar_solution():
    def sigma(x):
        if x<=0:return 0.0
        if x>=1:return 1.0
        odds=1/(1-x)**2-1/x**2
        return 1/(1+math.exp(-odds)) if odds>=0 else math.exp(odds)/(1+math.exp(odds))
    def rhs(x,v):
        F,H,G,C=v
        return ((.1-.6*sigma(x))*F,F-1.5*H,F*F-G,F*F)
    solution=solve_ivp(rhs,(0,1),(1,5/8,5/6,5),method='DOP853',rtol=3e-13,atol=3e-14,dense_output=True)
    assert solution.success
    return solution


def run():
    start=time.perf_counter();manufactured=exact_forcing_checks();oracle=independent_scalar_solution()
    print('Directed Taylor fifth-coefficient remainder contains exact constant-forcing/resonance solutions',flush=True)
    f=OriginalOuterSlopeVelocity();samples=[];numeric=physical=history=pressure=measure=0
    path=DIRECTORY/(PREFIX+'compliant_current_original_outer_reference_numeric.json')
    reference=json.loads(path.read_text(encoding='utf-8'))
    scalar_checks=0;max_width=f.ctx.mpf(0);max_point_difference=0.0
    for x in ('0','.25','.5','.75','1'):
        data=f.numeric.scalar.evaluate(x);expected=oracle.sol(float(x))
        for i,name in enumerate(f.numeric.scalar.names):
            lo,hi=_bounds(f.ctx,data['directed'][name]);max_width=max(max_width,hi-lo)
            tolerance=f.ctx.mpf('5e-12')
            value=f.ctx.mpf(float(expected[i]))
            assert lo-tolerance<=value<=hi+tolerance,(x,name,f.ctx.nstr(lo,20),f.ctx.nstr(hi,20),expected[i])
            max_point_difference=max(max_point_difference,abs(float(data['point'][name])-float(expected[i])));scalar_checks+=1
        print('Directed original slope scalar profiles: x='+x,flush=True)
    endpoint=f.numeric.scalar.evaluate('1');lo,hi=_bounds(f.ctx,endpoint['directed']['F'])
    exact=f.iv.exp(-f.iv.mpf(1)/5);ll,hh=_bounds(f.ctx,exact);assert max(lo,ll)<=min(hi,hh)
    for Z in ('-.3','0','.3','.17','.23'):
        for x in ('0','.25','.5','.75','1'):
            p=f.evaluate(Z,x,log_tau='-2',azimuth='.25');samples.append(p);inlet=f.numeric._inlets[Z]
            assert p['selected_integer']==inlet['selected_integer'];history+=1
            assert p['positive_order_n1_own_history']==inlet['positive_order_n1_own_history'];history+=1
            assert p['independent_P0_over_Pstar_squared']==inlet['independent_P0_over_Pstar_squared'];pressure+=1
            for k in range(3):
                P=f._restore(p['generic_leading_fields']['pressure'][k]);p0=f._restore(p['independent_P0_over_Pstar_squared'][k])
                own=f._restore(p['generic_leading_fields']['p'][k]);a,b=_bounds(f.ctx,P-p0-own);assert a<=0<=b;pressure+=1
            total=sum((Fraction(*s['exact_end'])-Fraction(*s['exact_left']) for s in p['physical_transport_steps']),Fraction(0))
            assert total==Fraction(x);measure+=1
            for s in p['physical_transport_steps']:
                decay=s['own_weights']['p']['decay'];assert decay['lower_exact_mpf_tuple']==decay['upper_exact_mpf_tuple']==[0,1,0,1];measure+=1
            if x=='0':
                for group in ('generic_modified_fields','physical_similarity','finite_N_history_errors'):
                    for name,rows in p[group].items():
                        for a,b in zip(rows,inlet[group][name]):overlap(a,b);numeric+=1
                old=next(a for a in reference['samples'] if a['Z']==Z and a['native_offset']=='0' and a['physical']['time']['kind']=='log_tau')
                physical+=inlet_physical_seam(f,old['physical'],p['physical'])
                assert p['finite_N_local_errors']==inlet['finite_N_local_errors'];history+=1
        print('Actual original O2 slope physical values: Z='+Z,flush=True)
    scaled=0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'1',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
        normal=next(a for a in samples if a['Z']==Z and a['native_slope']=='1')
        for name,terms in p['physical']['cylindrical'].items():
            for a,b in zip(terms,normal['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);scaled+=1
    boundary=0
    for x,expected in [('0.24'+'9'*138,1),('0.25'+'0'*137+'1',2)]:
        p=f.numeric.evaluate('.3',x);assert len(p['physical_transport_steps'])==expected;boundary+=1
    rejected=0
    for x,kwargs in (('-.1',{}),('1.1',{}),('nan',{}),('.5',dict(log_tau='1')),('.5',dict(azimuth='nan'))):
        try:f.evaluate('.3',x,**kwargs)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid O2 slope request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,exact_forcing_directed_Taylor_checks=manufactured,independent_DOP853_scalar_values=scalar_checks,
            independent_float_oracle_tolerance='5e-12; consistency only, not the source of directed bounds',scipy_version=scipy.__version__,
            max_scalar_directed_width=f.ctx.nstr(max_width,15),max_scalar_point_difference_from_independent_float=max_point_difference,
            endpoint_F_exact_symmetry_identity=True,actual_original_O2_slope_physical_packets=len(samples),Rref_numeric_seam_coefficients=numeric,
            Rref_physical_sum_seams=physical,N_and_nonzero_n1_history_checks=history,physical_measure_and_pressure_memory_checks=measure,
            P0_unit_and_once_only_checks=pressure,scaled_time_checks=scaled,exact_close_boundary_checks=boundary,rejected_requests=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(through_original_O2_slope=True,beyond_O2_slope=False,physical_Cartesian_values=True,native_coordinates_required=True,
            finite_N_values_are_point_plus_directed_error=True,resolved_signed_finite_N_centers=False,
            radial_fast_phase_derivatives=False,full_n1_pressure=False,reduced_m13_m15_integrals_continued=False,
            physical_gradients=False,direct_physical_inverse=False,global_n1_high_order_joins=False,global_N_admission=False,
            global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for n,d in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==d,n
    for n,d in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==d,n
    return r


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-report',action='store_true');p.add_argument('--inspect-report',action='store_true')
    args=p.parse_args();path=DIRECTORY/NAME;r=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=r['checks'],execution_seconds=r['execution_seconds']),indent=2),flush=True)
