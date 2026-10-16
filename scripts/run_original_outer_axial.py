"""Original O2 axial physical values with resolved positive endpoint memory."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_outer_axial_numeric import OriginalOuterAxialVelocity
from run_original_Rm_patch_velocity import inlet_physical_seam
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_outer_axial_numeric.json'


def independent_memory_integral(f,s,power,result):
    """Independent tanh-sinh quadrature of the original scalar integrand.

    Peak parameters choose coordinates only; this evaluates the defining
    sigmoid directly, rather than the interval Taylor recurrence.
    """
    c=mp.mp.clone();c.dps=65;proof=result['proof'];sp=c.mpf(s);D=c.expm1(40*sp);U=c.exp(40*sp)
    origin=c.mpf(proof['q_origin']);scale=c.mpf(proof['q_scale']);anchor=c.mpf(proof['normalized_log_anchor'])
    left,right=c.mpf(proof['integration_left']),c.mpf(proof['integration_right'])
    ordinary=proof['method'].startswith('exact exponential')
    def integrand(x):
        q=x if ordinary else origin+scale*x
        r=1-c.ln(U-q)/40
        if r>=1:logbeta=c.mpf(0)
        elif r<=0:return c.mpf(0)
        else:
            odds=1/(1-r)**2-1/r**2
            logbeta=-c.log1p(c.exp(-odds)) if odds>0 else odds-c.log1p(c.exp(odds))
        return c.exp(-q+power*logbeta-anchor)
    split=[left,right] if ordinary else [left,c.mpf(0),right]
    value=c.exp(-D)+c.exp(anchor)*scale*c.quad(integrand,split)
    lo,hi=_bounds(c,result['directed']);tol=abs(value)*c.mpf('1e-35')
    assert lo-tol<=value<=hi+tol,(s,power,c.nstr(c.ln(value),30))
    return c.nstr(abs(value-c.mpf(result['point']))/value,15)


def run():
    start=time.perf_counter();f=OriginalOuterAxialVelocity();samples=[]
    shape_checks=0;max_relative=f.ctx.mpf(0);endpoint_logs={};oracle=[]
    for s in ('0','.25','.5','.75','.9999993','.9999999','1'):
        results={k:f.numeric.memory.evaluate(s,k) for k in (1,2)}
        for k,r in results.items():
            lo,hi=_bounds(f.ctx,r['directed']);assert lo>0
            max_relative=max(max_relative,(hi-lo)/lo);shape_checks+=1
            if s=='1':endpoint_logs[str(k)]=f.ctx.nstr(f.ctx.ln(r['point']),30)
            if s in ('.25','.5','.9999993','1'):oracle.append(independent_memory_integral(f,s,k,r))
        a,b=_bounds(f.ctx,results[1]['directed']);aa,bb=_bounds(f.ctx,results[2]['directed'])
        assert aa<=b and a>=bb if s in ('1','.9999999') else aa<=b
        print('Original positive axial memories: s='+s,flush=True)
    reference=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_outer_slope_numeric.json')).read_text(encoding='utf-8'))
    numerical=physical=history=pressure=measure=endpoint=0
    def restore_step(row):return f.iv.mpf([mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))])
    for Z in ('-.3','0','.3','.17','.23'):
        for s in ('0','.25','.5','.75','1'):
            p=f.evaluate(Z,s,log_tau='-2',azimuth='.25');samples.append(p);inlet=f.numeric._inlets[Z]
            assert p['selected_integer']==inlet['selected_integer'];history+=1
            assert p['positive_order_n1_own_history']==inlet['positive_order_n1_own_history'];history+=1
            assert p['independent_P0_over_Pstar_squared']==inlet['independent_P0_over_Pstar_squared'];pressure+=1
            for k in range(3):
                value=f._restore(p['generic_leading_fields']['pressure'][k])-f._restore(p['generic_leading_fields']['p'][k])-f._restore(p['independent_P0_over_Pstar_squared'][k])
                a,b=_bounds(f.ctx,value);assert a<=0<=b;pressure+=1
            widths=sum((restore_step(v['physical_log_width']) for v in p['physical_transport_steps']),f.iv.mpf(0))
            a,b=_bounds(f.ctx,widths-f.iv.expm1(40*f.iv.mpf(s)));assert a<=0<=b;measure+=1
            for v in p['physical_transport_steps']:
                weights=v['own_weights']['p'];decay=weights['decay']
                assert decay['lower_exact_mpf_tuple']==decay['upper_exact_mpf_tuple']==[0,1,0,1];measure+=1
                a,b=_bounds(f.ctx,restore_step(weights['mass'])-restore_step(v['physical_log_width']));assert a<=0<=b;measure+=1
            if s=='0':
                for group in ('generic_modified_fields','physical_similarity','finite_N_history_errors'):
                    for name,rows in p[group].items():
                        for a,b in zip(rows,inlet[group][name]):overlap(a,b);numerical+=1
                old=next(a for a in reference['samples'] if a['Z']==Z and a['native_slope']=='1' and a['physical']['time']['kind']=='log_tau')
                physical+=inlet_physical_seam(f,old['physical'],p['physical'])
                assert p['finite_N_local_errors']==inlet['finite_N_local_errors'];history+=1
            if s=='1':
                assert p['axial_profiles']['beta']['approximation']=='0.0';endpoint+=1
                for n in ('M','W'):
                    assert _bounds(f.ctx,f._restore(p['axial_profiles'][n]))[0]>0;endpoint+=1
                if Z=='.3':
                    assert _bounds(f.ctx,f._restore(p['physical_similarity']['Ur_over_sqrt_R_over2'][0]))[1]<0;endpoint+=1
        print('Actual original O2 axial physical values: Z='+Z,flush=True)
    for s in ('.9999993','.9999999'):samples.append(f.evaluate('.3',s,log_tau='-2',azimuth='.25'))
    scaled=0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'1',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
        normal=next(a for a in samples if a['Z']==Z and a['native_axial']=='1')
        for name,terms in p['physical']['cylindrical'].items():
            for a,b in zip(terms,normal['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);scaled+=1
    rejected=0
    for s,kwargs in (('-.1',{}),('1.1',{}),('nan',{}),('.5',dict(log_tau='1')),('.5',dict(azimuth='nan'))):
        try:f.evaluate('.3',s,**kwargs)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original axial request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,positive_original_memory_checks=shape_checks,independent_tanh_sinh_memory_comparisons=len(oracle),
            independent_memory_relative_point_differences=oracle,max_memory_relative_interval_width=f.ctx.nstr(max_relative,15),
            endpoint_log_positive_memories=endpoint_logs,actual_original_O2_axial_physical_packets=len(samples),
            slope_numeric_seam_coefficients=numerical,slope_physical_sum_seams=physical,N_and_nonzero_n1_history_checks=history,
            physical_exponential_measure_and_pressure_memory_checks=measure,P0_unit_and_once_only_checks=pressure,
            endpoint_leading_zero_with_positive_history_checks=endpoint,scaled_time_checks=scaled,rejected_requests=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(through_original_O2_axial=True,beyond_O2_axial=False,physical_Cartesian_values=True,native_coordinates_required=True,
            finite_N_values_are_point_plus_directed_error=True,resolved_signed_finite_N_centers=False,radial_fast_phase_derivatives=False,
            full_n1_pressure=False,reduced_m13_m15_integrals_continued=False,physical_gradients=False,direct_physical_inverse=False,
            global_n1_high_order_joins=False,global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


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
