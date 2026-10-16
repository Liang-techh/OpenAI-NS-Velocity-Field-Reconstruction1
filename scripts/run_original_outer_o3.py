"""Focused original O3 transition/power continuation, with nonzero mu sectors."""
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

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_o3_scalar import OriginalO3Scalar
from openai_ns_reconstruction.original_outer_o3_numeric import OriginalOuterO3Velocity
from run_original_Rm_patch_velocity import inlet_physical_seam
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_outer_o3_numeric.json'


def independent_scalar_checks():
    c=mp.mp.clone();c.dps=70;iv=MPIntervalContext();iv.dps=90;mu=c.mpf('.001')
    simple=OriginalO3Scalar(c,iv,mu,iv.mpf('.001'),cells=64)
    simple.sigma_jets=lambda ctx,x:[ctx.mpf('.5')]+[ctx.mpf(0)]*4
    constant=0
    for x in ('.125','.5','1'):
        t=c.mpf(x);got=simple.transition(x)
        expected=dict(J=t/2,U1=c.expm1(-mu*t/2)/mu,U2=c.expm1(-mu*t)/mu,
            C1=c.quad(lambda u:c.exp(u)*c.expm1(-mu*u/2)/mu,[0,t]),
            C0=c.quad(lambda u:c.expm1(-mu*u)/mu,[0,t]),
            Cm=c.quad(lambda u:c.exp(-u)*c.expm1(-mu*u)/mu,[0,t]))
        for n,v in expected.items():
            lo,hi=_bounds(c,got['directed'][n]);assert lo-c.mpf('1e-60')<=v<=hi+c.mpf('1e-60'),(x,n);constant+=1
    power=0
    for x in ('.125','1','2'):
        t=c.mpf(x);got=simple.power(x)
        for name,alpha,factor in (('C1',1,1),('C0',0,2),('Cm',-1,2)):
            v=c.quad(lambda u:c.exp(alpha*u)*c.expm1(-factor*mu*u/2)/mu,[0,1])
            v+=c.quad(lambda u:c.exp(alpha*(1+u))*c.expm1(-factor*mu*(c.mpf('.5')+u))/mu,[0,t])
            lo,hi=_bounds(c,got['directed'][name]);assert lo<=v<=hi,(x,name);power+=1
    # Independent DOP853 integration of the defining scalar densities for
    # the actual original sigmoid, with a representable nonzero test mu.
    import scipy
    from scipy.integrate import solve_ivp
    from scipy.special import expit
    def rhs(t,y):
        sigma=0. if t<=0 else 1. if t>=1 else expit(1/(1-t)**2-1/t**2)
        U1=math.expm1(-.001*y[0])/.001;U2=math.expm1(-.002*y[0])/.001
        return [sigma,math.exp(t)*U1,U2,math.exp(-t)*U2]
    targets=[.25,.5,.75,1.];result=solve_ivp(rhs,[0,1],[0.,0.,0.,0.],method='DOP853',t_eval=targets,rtol=3e-14,atol=3e-15)
    assert result.success;actual=OriginalO3Scalar(c,iv,mu,iv.mpf('.001'));count=0;difference=c.mpf(0)
    for j,x in enumerate(targets):
        got=actual.transition(str(x))
        for i,name in enumerate(('J','C1','C0','Cm')):
            v=c.mpf(float(result.y[i,j]));lo,hi=_bounds(c,got['directed'][name]);tol=c.mpf('5e-12')
            assert lo-tol<=v<=hi+tol,(x,name);difference=max(difference,abs(got['point'][name]-v));count+=1
    return dict(independent_constant_sigmoid_integrals=constant,independent_power_defining_integrals=power,
        independent_DOP853_original_sigmoid_values=count,independent_float_tolerance='5e-12; consistency only',
        max_scalar_point_difference=c.nstr(difference,16),scipy_version=scipy.__version__)


def run():
    start=time.perf_counter();independent=independent_scalar_checks()
    print('Independent normalized transition and original power-flow integrals pass',flush=True)
    f=OriginalOuterO3Velocity();samples=[];numerical=physical=history=pressure=measure=sectors=quiet=0
    accepted=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_outer_buffer_numeric.json')).read_text(encoding='utf-8'))
    widths=[]
    for Z in ('-.3','0','.3','.17','.23'):
        for chart,values in (('O3_transition',('0','.125','.25','.5','.75','1')),('O3_power',('0','.5','1','1.5','2'))):
            for x in values:
                p=f.evaluate(Z,chart,x,log_tau='-2',azimuth='.25');samples.append(p)
                inlet=f.numeric._inlets[Z] if chart=='O3_transition' else f.numeric.evaluate(Z,'O3_transition','1')
                assert p['selected_integer']==inlet['selected_integer'];history+=1
                assert p['positive_order_n1_own_history']==inlet['positive_order_n1_own_history'];history+=1
                assert p['independent_P0_over_Pstar_squared']==inlet['independent_P0_over_Pstar_squared'];pressure+=1
                for k in range(3):
                    v=f._restore(p['generic_leading_fields']['pressure'][k])-f._restore(p['generic_leading_fields']['p'][k])-f._restore(p['independent_P0_over_Pstar_squared'][k])
                    lo,hi=_bounds(f.ctx,v);assert lo<=0<=hi;pressure+=1
                width=sum((Fraction(*s['exact_end'])-Fraction(*s['exact_left']) for s in p['physical_transport_steps']),Fraction(0))
                assert width==Fraction(x);measure+=1
                for step in p['physical_transport_steps']:
                    decay=step['own_weights']['p']['decay'];assert decay['lower_exact_mpf_tuple']==decay['upper_exact_mpf_tuple']==[0,1,0,1];measure+=1
                if x=='0':
                    for group in ('generic_modified_fields','physical_similarity','finite_N_history_errors'):
                        for name,rows in p[group].items():
                            for a,b in zip(rows,inlet[group][name]):overlap(a,b);numerical+=1
                    old=next(a for a in accepted['samples'] if a['Z']==Z and a['native_buffer']=='11' and a['physical']['time']['kind']=='log_tau') if chart=='O3_transition' else next(a for a in samples if a['Z']==Z and a['native_O3_chart']=='O3_transition' and a['native_O3']=='1')
                    physical+=inlet_physical_seam(f,old['physical'],p['physical'])
                if chart=='O3_power' or x!='0':
                    for name,sign in (('E',-1),('h',-1),('e',1),('p',-1)):
                        lo,hi=_bounds(f.ctx,f._restore(p['leading_small_sectors']['mu_correction_'+name][0]))
                        assert hi<0 if sign<0 else lo>0,(Z,chart,x,name);sectors+=1
                assert _bounds(f.ctx,f._restore(p['original_mu']['mu']))[0]>0;sectors+=1
                for row in p['normalized_mu_profiles'].values():
                    lo,hi=_bounds(f.ctx,f._restore(row));widths.append(hi-lo)
                if chart=='O3_power':
                    for name,rows in p['finite_N_local_errors'].items():
                        for row in rows:assert _bounds(f.ctx,f._restore(row))==(0,0);quiet+=1
                    for name,rows in p['finite_N_history_errors'].items():
                        rate=dict(m=1,h='1.5',k='1.5',e=1,p=0)[name]
                        for a,b in zip(rows,inlet['finite_N_history_errors'][name]):
                            v=f._restore(a)*f.iv.exp(f.iv.mpf(rate)*f.iv.mpf(x))-f._restore(b)
                            lo,hi=_bounds(f.ctx,v);assert lo<=0<=hi;quiet+=1
                    # Incoming nonzero intervals survive even though local A/B vanish.
                    assert any(_bounds(f.ctx,f._restore(v))[0]!=_bounds(f.ctx,f._restore(v))[1] for v in p['finite_N_history_errors']['h']);quiet+=1
        print('Actual original O3 transition/power through Rc: Z='+Z,flush=True)
    scaled=0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'O3_power','2',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
        normal=next(a for a in samples if a['Z']==Z and a['native_O3_chart']=='O3_power' and a['native_O3']=='2')
        for name,terms in p['physical']['cylindrical'].items():
            for a,b in zip(terms,normal['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);scaled+=1
    boundary=0
    for chart,x,count in [('O3_transition','0.24'+'9'*138,1),('O3_transition','0.25'+'0'*137+'1',2),('O3_power','0.'+'9'*140,1),('O3_power','1.'+'0'*139+'1',2)]:
        p=f.numeric.evaluate('.3',chart,x);assert len(p['physical_transport_steps'])==count;boundary+=1
    rejected=0
    for chart,x,kwargs in [('O3_transition','-.1',{}),('O3_transition','1.1',{}),('O3_power','2.1',{}),('O3_power','nan',{}),('bad','0',{}),('O3_power','1',dict(log_tau='1')),('O3_power','1',dict(azimuth='nan'))]:
        try:f.evaluate('.3',chart,x,**kwargs)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid O3 request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,**independent,actual_original_O3_physical_packets=len(samples),buffer_and_transition_numeric_seam_coefficients=numerical,
            buffer_and_transition_physical_sum_seams=physical,N_and_nonzero_n1_history_checks=history,
            physical_measure_and_pressure_memory_checks=measure,P0_unit_and_once_only_checks=pressure,
            nonzero_mu_and_signed_sector_checks=sectors,quiet_power_zero_local_with_nonzero_history_checks=quiet,
            max_tested_normalized_profile_interval_width=f.ctx.nstr(max(widths),16),exact_close_boundary_checks=boundary,
            scaled_time_checks=scaled,rejected_requests=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(through_original_Rc=True,beyond_Rc=False,physical_Cartesian_values=True,native_coordinates_required=True,
            finite_N_values_are_point_plus_directed_error=True,resolved_signed_finite_N_centers=False,radial_fast_phase_derivatives=False,
            full_n1_pressure=False,reduced_m13_m15_integrals_continued=False,relative_pressure_cancellation_accuracy=False,
            admissible_stress_cone_margin=False,physical_gradients=False,direct_physical_inverse=False,global_n1_high_order_joins=False,
            global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


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
