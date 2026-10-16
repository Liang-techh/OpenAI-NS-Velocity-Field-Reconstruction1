"""Original O2 buffer exact five-history and small-sector physical continuation."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_core_log_amplitude import OriginalCoreLogAmplitude
from openai_ns_reconstruction.original_outer_buffer_numeric import OriginalOuterBufferNumeric,OriginalOuterBufferVelocity
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from run_original_Rm_patch_velocity import inlet_physical_seam
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_outer_buffer_numeric.json'


def independent_five_ODEs():
    c=mp.mp.clone();c.dps=70;iv=MPIntervalContext();iv.dps=90
    a=object.__new__(OriginalCoreLogAmplitude);a.ctx=c;a.iv=iv;a.precision=65
    f=object.__new__(OriginalOuterBufferNumeric);f.ctx=c;f.iv=iv;f.J=J
    restore=lambda r:iv.mpf([mp.make_mpf(tuple(v)) for v in r['source_enclosure_exact_mpf_tuples']])
    f.patch=SimpleNamespace(leading=SimpleNamespace(defect=SimpleNamespace(logPp=c.ln(3),logP=iv.ln(3)),
        reference=SimpleNamespace(previous=SimpleNamespace(_restore=restore))))
    row=lambda v:a.summarize(c.mpf(v),iv.mpf(v))
    inlet=dict(original_slope_scalar_inlet={n:row(v) for n,v in zip(('F','H','G','C'),('.8','.63','.9','5.9'))},
        axial_profiles={n:row(v) for n,v in zip(('M','W','physical_log_distance'),('.3','.2','1.7'))})
    count=radial=0
    for Z in ('-.3','0','.3'):
        z=c.mpf(Z);A=1/(1+z*z);B=4*z/3;D=c.mpf('1.7');F1,H1,G1,C1=map(c.mpf,('.8','.63','.9','5.9'))
        E0=A*F1*c.exp(-D/2)
        initial=dict(m=B*c.mpf('.3'),h=A*(F1*c.exp(-D/2)+(H1-F1)*c.exp(-c.mpf('1.5')*D)),
            k=A*B*F1*c.exp(-D/2)*(c.mpf('.3')+(H1/F1-1)*c.exp(-D)),
            e=B*B*c.mpf('.2')-A*A*(G1+F1*F1*D)*c.exp(-D)/2,p=A*A*(C1+F1*F1*(1-c.exp(-D)))/2)
        for b in ('0','.5','11'):
            x=c.mpf(b);truth,sectors,Q=f._leading(Z,b,inlet,iv)
            def density(u):
                E=E0*c.exp(-u/2)
                return dict(m=c.mpf(0),h=E,k=c.mpf(0),e=-E*E/2,p=E*E/2)
            for name,rate in dict(m=1,h=c.mpf('1.5'),k=c.mpf('1.5'),e=1,p=0).items():
                value=initial[name]*c.exp(-rate*x)+c.quad(lambda u:c.exp(-rate*(x-u))*density(u)[name],[0,x])
                lo,hi=_bounds(c,truth[name][0]);tol=c.mpf('1e-60')*max(1,abs(value))
                assert lo-tol<=value<=hi+tol,(Z,b,name);count+=1
            # Independently recover Q from m and its actual Z derivative.
            delta=iv.mpf('.11');zz=J.variable(iv,iv.mpf(Z),2);mass=truth['m'];mz=J(iv,[mass[1],2*mass[2]])
            generic=(-zz*mass*(1-delta)-(1-zz*zz)*mz)/(1-zz*zz*delta)*3
            for k in range(2):
                lo,hi=_bounds(c,generic[k]-Q[k]);assert lo<=0<=hi;radial+=1
    return count,radial


def run():
    start=time.perf_counter();ode,radial=independent_five_ODEs()
    print('Independent raw five-density ODEs and shared-streamfunction Q recovery pass',flush=True)
    f=OriginalOuterBufferVelocity();samples=[];numerical=physical=history=pressure=measure=sectors=0
    reference=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_outer_axial_numeric.json')).read_text(encoding='utf-8'))
    for Z in ('-.3','0','.3','.17','.23'):
        for b in ('0','.5','1','5','8','11'):
            p=f.evaluate(Z,b,log_tau='-2',azimuth='.25');samples.append(p);inlet=f.numeric._inlets[Z]
            assert p['selected_integer']==inlet['selected_integer'];history+=1
            assert p['positive_order_n1_own_history']==inlet['positive_order_n1_own_history'];history+=1
            assert p['independent_P0_over_Pstar_squared']==inlet['independent_P0_over_Pstar_squared'];pressure+=1
            for k in range(3):
                value=f._restore(p['generic_leading_fields']['pressure'][k])-f._restore(p['generic_leading_fields']['p'][k])-f._restore(p['independent_P0_over_Pstar_squared'][k])
                a,bound=_bounds(f.ctx,value);assert a<=0<=bound;pressure+=1
            width=sum((Fraction(*s['exact_end'])-Fraction(*s['exact_left']) for s in p['physical_transport_steps']),Fraction(0))
            assert width==Fraction(b);measure+=1
            for s in p['physical_transport_steps']:
                decay=s['own_weights']['p']['decay'];assert decay['lower_exact_mpf_tuple']==decay['upper_exact_mpf_tuple']==[0,1,0,1];measure+=1
            assert _bounds(f.ctx,f._restore(p['leading_small_sectors']['h_minus_E'][0]))[1]<0;sectors+=1
            assert _bounds(f.ctx,f._restore(p['leading_small_sectors']['p_tail'][0]))[1]<0;sectors+=1
            if b!='0':
                assert _bounds(f.ctx,f._restore(p['leading_small_sectors']['p_buffer_increment'][0]))[0]>0;sectors+=1
            if b=='0':
                for group in ('generic_modified_fields','physical_similarity','finite_N_history_errors'):
                    for name,rows in p[group].items():
                        for a,bb in zip(rows,inlet[group][name]):overlap(a,bb);numerical+=1
                old=next(a for a in reference['samples'] if a['Z']==Z and a['native_axial']=='1' and a['physical']['time']['kind']=='log_tau')
                physical+=inlet_physical_seam(f,old['physical'],p['physical'])
                assert p['finite_N_local_errors']==inlet['finite_N_local_errors'];history+=1
            assert _bounds(f.ctx,f._restore(p['leading_Q_exact_Z_independent'][1]))==(f.ctx.mpf(0),f.ctx.mpf(0));radial+=1
        print('Actual original O2 buffer fields and retained sectors: Z='+Z,flush=True)
    tiny=f.numeric.evaluate('.3','1e-140')
    assert _bounds(f.ctx,f._restore(tiny['leading_small_sectors']['p_buffer_increment'][0]))[0]>0;sectors+=1
    boundary=0
    for b,count in [('0.'+'9'*140,1),('1.'+'0'*139+'1',2),('4.'+'9'*140,2),('5.'+'0'*139+'1',3)]:
        p=f.numeric.evaluate('.3',b);assert len(p['physical_transport_steps'])==count;boundary+=1
    scaled=0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'11',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
        normal=next(a for a in samples if a['Z']==Z and a['native_buffer']=='11')
        for name,terms in p['physical']['cylindrical'].items():
            for a,b in zip(terms,normal['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);scaled+=1
    rejected=0
    for b,kwargs in (('-.1',{}),('11.1',{}),('nan',{}),('1',dict(log_tau='1')),('1',dict(azimuth='nan'))):
        try:f.evaluate('.3',b,**kwargs)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original buffer request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_raw_five_ODE_values=ode,streamfunction_Q_identity_checks=radial,
            actual_original_O2_buffer_physical_packets=len(samples),axial_numeric_seam_coefficients=numerical,axial_physical_sum_seams=physical,
            N_and_nonzero_n1_history_checks=history,physical_measure_and_pressure_memory_checks=measure,
            P0_unit_and_once_only_checks=pressure,retained_nonzero_small_sector_checks=sectors,
            tiny_pressure_increment_at_b_1e_minus140=True,exact_close_boundary_checks=boundary,scaled_time_checks=scaled,rejected_requests=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(through_original_O2_buffer=True,beyond_O2_buffer=False,physical_Cartesian_values=True,native_coordinates_required=True,
            finite_N_values_are_point_plus_directed_error=True,resolved_signed_finite_N_centers=False,radial_fast_phase_derivatives=False,
            full_n1_pressure=False,reduced_m13_m15_integrals_continued=False,relative_pressure_cancellation_accuracy=False,
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
