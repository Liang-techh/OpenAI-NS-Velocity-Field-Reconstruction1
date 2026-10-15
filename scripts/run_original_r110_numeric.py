"""Evaluate original-N second/post-power values and six own histories to R110."""
import argparse
import hashlib
import importlib
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_r110_numeric import OriginalSecondSwitchNumeric,OriginalPostPowerNumeric,post_power_histories,RATES
from openai_ns_reconstruction.original_core_numeric import _bounds
from run_original_second_micro_numeric import overlap

NAME='lei_ren_part1_paper_compliant_current_original_r110_numeric.json'


def independent_checks():
    c=mp.mp.clone();c.dps=65;iv=MPIntervalContext();iv.dps=110
    module=importlib.import_module('lei_ren_part1_paper_compliant_actual_bridge_integrals');J=module.IntervalTaylor
    h=iv.mpf(1)/101;hp=c.mpf(1)/101;z=c.mpf('.3')
    b=object.__new__(module.CompliantActualBridgeIntegrals);b.ctx=iv;b.cap=h;b.weight=iv.mpf('.1');b.pulse_cache={}
    funcs=[lambda z:2+z,lambda z:-1+z*z,lambda z:1-z+z**3]
    def jet(fn,ctx=iv,order=3):
        return J(ctx,[ctx.mpf(c.nstr(c.diff(fn,z,k)/math.factorial(k),65)) for k in range(order+1)])
    modes=[jet(lambda z,fn=fn:fn(z)/10000) for fn in funcs]
    p=dict(R0=iv.mpf(2),modes=dict(D_over_R=modes))
    leading=SimpleNamespace(iv=iv,bridge=b,hcap=h,logh=iv.ln(h),_prepared={'.3':p},
        _transfer=lambda v:iv.mpf([mp.make_mpf(x) for x in v._mpi_]),_width_error=lambda v:h*abs(v))
    fake=object.__new__(OriginalSecondSwitchNumeric);fake.leading=leading
    fake.transport=SimpleNamespace(ctx=c,iv=iv,phase_iv=iv.mpf(1)/32);fake._second_switch_increments={}
    fake.Y=iv.ln(50);fake.Ypoint=c.ln(50);fake.post_log=iv.ln(iv.mpf(11)/10)
    fake.post_log_point=c.ln(c.mpf(11)/10);fake.post_length=fake.post_log-2*h;fake.logN=iv.ln(17)
    post=object.__new__(OriginalPostPowerNumeric);post.__dict__=dict(fake.__dict__)
    def sigma(t):
        if t<=0:return c.mpf(0)
        if t>=1:return c.mpf(1)
        odds=1/(1-t)**2-1/t**2
        if odds>1000:return c.mpf(1)
        if odds<-1000:return c.mpf(0)
        return 1/(1+c.exp(-odds))
    def direct_log(z,s):
        D=lambda t:sum((fn(z)/10000*c.mpf(2)**j*100**(1-j)*c.exp((1-j)*hp*(1+t)) for j,fn in enumerate(funcs)),c.mpf(0))
        return -hp*hp/2*c.quad(lambda t:(1-sigma(t))*D(t),[0,s/2,s])-c.mpf(2)/5*hp*c.quad(sigma,[0,s/2,s])
    checks=0
    for phase in ('.4','1'):
        s=c.mpf(phase);bound,evidence=fake._second_increment('.3',iv.mpf(phase))
        for k in range(4):
            value=abs(c.diff(lambda z:c.expm1(direct_log(z,s)),z,k)/math.factorial(k))
            assert value<=_bounds(c,bound*iv.exp(bound)/b.weight**k)[1];checks+=1
        assert evidence['native_shear_width_power']==2 and evidence['native_power_width_power']==1;checks+=1
        for view,actual in ((fake,4*c.exp(c.ln(50)+hp*(1+s))),
                            (post,4*c.exp(c.ln(50)+c.ln(c.mpf(11)/10)*s+2*hp*(1-s)))):
            rho,_,_=view._downstream_geometry(iv.mpf(phase));lo,hi=_bounds(c,rho)
            assert lo-c.mpf('1e-55')<=actual<=hi+c.mpf('1e-55');checks+=1
    # Independently integrate all six defining ODEs, with nonzero inlet and
    # genuine axial Taylor coefficients, rather than repeat closed weights.
    initials={name:(lambda z,n=n:1+n+z+z*z/(n+1)) for n,name in enumerate(RATES)}
    phi=lambda z:1+z*z;psi=lambda z:3-z
    for y in (c.mpf('.07'),c.mpf('.4')):
        result=post_power_histories(c,{name:jet(fn,c,2) for name,fn in initials.items()},jet(phi,c,2),jet(psi,c,2),y)
        for name,rate in RATES.items():
            def exact(z):
                def source(t):
                    p=phi(z)*c.exp(-c.mpf(2)*t/5);v=psi(z)
                    return dict(H=2*p,Mpsi=v,Kpsi=2*p*v,Apsi=v*v,B=p*p,C=p*p)[name]
                return initials[name](z)*c.exp(-rate*y)+c.quad(lambda t:c.exp(-rate*(y-t))*source(t),[0,y])
            for k in range(3):
                value=c.diff(exact,z,k)/math.factorial(k)
                assert abs(value-result[name][k])<c.mpf('1e-54');checks+=1
    rows=[dict(exact_zero=True)]*4
    def cell(a,b):return dict(exact_left=a,exact_right=b,actual_source_and_primitive_range=dict(A=rows,B=rows),actual_density_ranges=dict(N_scaled={name:rows for name in ('m','h','k','e','p')}))
    windows=[dict(chart=name,cells=[cell(a,b)]) for name,a,b in (
        ('first_micro',['sc_half'],[1,1]),('second_micro',[1,1],[2,1]),('frozen_macro',[0,1],[1,1]),
        ('first_switch',[0,1],[1,1]),('second_switch',[0,1],[1,1]),('post_power',[0,1],[1,1]))]
    fake.ranges=post.ranges=[dict(exact_Z_cell=['0','.5'],actual_17_chart_C3_transports=windows)]
    for view in (fake,post):
        local,prefix,_=view._covers('.3',iv.mpf('.4'));current=[p for p in prefix if p[3]==view.chart]
        assert len(current)==1 and _bounds(c,current[0][1])[1]>=c.mpf('.4');checks+=1
        measure=iv.ln(view.post_log) if view.chart=='post_power' else leading.logh
        assert current[0][2]._mpi_==measure._mpi_;checks+=1
    return checks


def run():
    began=time.perf_counter();independent=independent_checks()
    print('Independent second controls, six ODEs, geometry and measures pass',flush=True)
    f=OriginalSecondSwitchNumeric();post=f.post_power_field();samples=[];proofs=0
    for Z in ('-.3','0','.3'):
        for view in (f,post):
            for s in ('0','.5','1'):
                p=view.evaluate(Z,s,log_tau='-2',azimuth='.25');samples.append(p)
                assert p['source_cover']['local_source_chart']==view.chart
                assert [len(p['normalized_modified_leading'][n]) for n in ('Phi','Psi','Mpsi','Qpsi')]==[4,4,3,2]
                assert set(p['unmodified_own_six'])==set(RATES)
                assert all(len(v)==3 for v in p['unmodified_own_six'].values())
                assert p['unmodified_six_contract']['selected_N_five_histories_remain_separate']
                assert not p['limits']['global_field'] and not p['limits']['actual_mean_Z3_numeric']
                for row in p['small_error_proofs']:assert row['numerical_absolute_error_cap']=='1e-1000';proofs+=1
        print('Actual second/post-power and six histories: Z='+Z,flush=True)
    for view in (f,post):samples.append(view.evaluate('.17','.37',log_tau='-3',azimuth='1.1'))
    reference=json.loads((DIRECTORY/'lei_ren_part1_paper_compliant_current_original_first_switch_numeric.json').read_text(encoding='utf-8'))
    pick=lambda Z,chart,s:next(p for p in samples if p['Z']==Z and p['source_cover']['local_source_chart']==chart and p['coordinate']['value']==s)
    seams=history=six_seams=geometry=mean_checks=0
    for Z in ('-.3','0','.3'):
        old=next(p for p in reference['samples'] if p['Z']==Z and p['coordinate']['value']=='1')
        second0,second1=pick(Z,'second_switch','0'),pick(Z,'second_switch','1')
        post0,exit=pick(Z,'post_power','0'),pick(Z,'post_power','1')
        for left,right in ((old,second0),(second1,post0)):
            for name,rows in left['normalized_modified_leading'].items():
                for a,b in zip(rows,right['normalized_modified_leading'][name]):overlap(a,b);seams+=1
            for name,rows in left['physical']['cylindrical'].items():
                for a,b in zip(rows,right['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);seams+=1
            assert left['original_five_correction_history_bounds']==right['original_five_correction_history_bounds'];history+=20
        for name,rows in second1['unmodified_own_six'].items():
            for a,b in zip(rows,post0['unmodified_own_six'][name]):overlap(a,b);six_seams+=1
        for p in [p for p in samples if p['Z']==Z]:
            for a,b in zip(p['unmodified_own_six']['Mpsi'],p['normalized_modified_leading']['Mpsi']):overlap(a,b);mean_checks+=1
            own=p['physical']['n1_own_history'];assert own['actual_terminal_primitive']==old['physical']['n1_own_history']['actual_terminal_primitive'];history+=1
        rho=post._numeric(Z,'1',False)[3];lo,hi=_bounds(f.transport.ctx,rho*f.leading.core.epsilon_iv)
        assert lo<=110<=hi and exit['physical']['native']['R110_endpoint_exact'];geometry+=1
        assert exit['normalized_modified_leading']['Phi'][0]['approximation']!=second1['normalized_modified_leading']['Phi'][0]['approximation'];geometry+=1
    coarse=OriginalSecondSwitchNumeric(precision=60,amplitude_precision=90).post_power_field()
    old=coarse.evaluate('.3','1',log_tau='-2',azimuth='.25');new=pick('.3','post_power','1');comparisons=0
    for group in ('normalized_modified_leading','unmodified_own_six'):
        for name,rows in new[group].items():
            for a,b in zip(rows,old[group][name]):overlap(a,b);comparisons+=1
    for name,rows in new['physical']['cylindrical'].items():
        for a,b in zip(rows,old['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);comparisons+=1
    rejected=0
    for view in (f,post):
        for Z,s,kwargs in (('1','.5',{}),('.3','-.01',{}),('.3','1.01',{}),('.3','nan',{}),('.3','1',{'log_tau':'1'}),('.3','1',{'azimuth':'nan'})):
            try:view.evaluate(Z,s,**kwargs)
            except (ValueError,ArithmeticError):rejected+=1
            else:raise AssertionError('Unsupported downstream request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,checks=dict(passed=True,independent_second_controls_six_ODE_geometry_measure_checks=independent,
        actual_second_and_post_power_packets=len(samples),selected_N_small_error_proofs=proofs,
        first_second_and_post_native_physical_seams=seams,unmodified_six_history_join_coefficients=six_seams,
        correction_and_n1_terminal_memory_checks=history,unmodified_modified_own_mean_overlaps=mean_checks,
        actual_R110_geometry_and_nonconstant_swirl_checks=geometry,cross_precision_native_physical_and_six_history_overlaps=comparisons,
        rejected_unsupported_requests=rejected),samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        limitations=dict(Phi_Psi_Z_order=3,modified_own_mean_Z_order=2,Qpsi_Z_order=1,unmodified_six_Z_order=2,
            modified_all_six_merged=False,resolved_signed_correction_integrals=False,resolved_fast_phase=False,
            mean_Z3_and_QZZ_numeric=False,total_radial_derivatives=False,long_reshape=False,global_N_admission=False,
            global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False),execution_seconds=round(time.perf_counter()-began,3))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for name,digest in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    r=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=r['checks'],execution_seconds=r['execution_seconds']),indent=2),flush=True)
