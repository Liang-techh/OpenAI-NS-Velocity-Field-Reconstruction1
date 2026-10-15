"""Evaluate selected-N second-micro values, Z3 velocities and retained histories."""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
DIRECTORY=ROOT/'experiments/root_st073';sys.path.insert(0,str(DIRECTORY))
from openai_ns_reconstruction.original_second_micro_numeric import OriginalSecondMicroNumeric
from openai_ns_reconstruction.original_core_numeric import _bounds

NAME='lei_ren_part1_paper_compliant_current_original_second_micro_numeric.json'


def independent_control_and_memory_checks():
    c=mp.mp.clone();c.dps=90;iv=MPIntervalContext();iv.dps=130
    module=importlib.import_module('lei_ren_part1_paper_compliant_actual_bridge_integrals')
    J=module.IntervalTaylor
    bridge=object.__new__(module.CompliantActualBridgeIntegrals)
    bridge.ctx=iv;bridge.cap=iv.mpf(0);bridge.logh=-iv.ln(7);bridge.weight=iv.mpf('.1')
    bridge.bridge=SimpleNamespace(r=iv.mpf(4));bridge.pressure_log=iv.ln(11);bridge.swirl_log=iv.ln(13)
    bridge.pulse_masses=lambda s:(iv.mpf('.5'),iv.mpf('.5'))
    direction=lambda D,drives:dict(D_over_R=J.constant(iv,D,2),**{
        'drive_'+name:J.constant(iv,value,2) for name,value in zip(('hydro','pressure','swirl'),drives)})
    p=dict(data=dict(phi0=J.constant(iv,1,2)),angular_prefix_norm_per_hb=iv.mpf(1),
        micros=[dict(direction=direction(2,(5,11,17)),quotient=J.constant(iv,1,2)),
                dict(direction=direction(3,(7,13,19)),quotient=J.constant(iv,iv.mpf(6)/5,2))])
    checks=0;h=c.mpf(1)/7
    for text in ('1','1.25','1.6','2'):
        angular,axial=bridge.contributions(p,'second',text)
        nominal=lambda terms:sum((c.mpf(_bounds(c,t['normalized_axial_coefficients'][0])[1])*h**t['width_power'] for t in terms),c.mpf(0))
        sigma=lambda u:3*u*u-2*u*u*u
        mass=c.quad(lambda u:h*(1-sigma(u)+h*sigma(u)),[0,1])
        expected=-2*(2*mass+3*h*h*(c.mpf(text)-1))
        assert abs(nominal(angular)-expected)<c.mpf('1e-75');checks+=1
        for name,d1,d2,power in (('hydro',5,7,1),('pressure',11,13,1),('swirl',17,19,2)):
            expected=-c.mpf(4)**power*(d1*mass+c.mpf(6)/5*d2*h*h*(c.mpf(text)-1))
            assert abs(nominal(axial[name])-expected)<c.mpf('1e-75');checks+=1
    # Genuine third differentiation of a finite exponential, independent of
    # the production Taylor coefficient bound formula.
    for z in ('.2','.3','.4'):
        x=c.mpf(z);N=17;A=lambda u:1+u+u*u+u**3
        coeff=[iv.mpf(c.nstr(A(x),90))/N,iv.mpf(c.nstr(1+2*x+3*x*x,90))/N,
               iv.mpf(c.nstr(1+3*x,90))/N,iv.mpf(1)/N]
        bound=OriginalSecondMicroNumeric._third_exponential_increment_bound(coeff)
        value=abs(c.diff(lambda u:c.exp(A(u)/N)-1,x,3)/6)
        assert value<=_bounds(c,bound)[1];checks+=1
    fake=object.__new__(OriginalSecondMicroNumeric)
    fake.transport=SimpleNamespace(ctx=c,iv=iv,phase_iv=iv.mpf(1)/32)
    fake.leading=SimpleNamespace(core=SimpleNamespace(lp_iv=iv.mpf(0),ll_iv=iv.mpf(0)),logh=-iv.ln(7))
    fake.logN=iv.ln(17)
    def cap(v):
        if not v:return dict(exact_zero=True)
        q=iv.ln(v)
        return dict(exact_zero=False,log_absolute_upper=dict(lower_exact_mpf_tuple=list(q._mpi_[0]),upper_exact_mpf_tuple=list(q._mpi_[1])))
    rows=[cap(v) for v in (3,3,2,0)]
    def cell(left,right):return dict(exact_left=left,exact_right=right,actual_source_and_primitive_range=dict(A=rows,B=rows),
        actual_density_ranges=dict(N_scaled={name:rows for name in ('m','h','k','e','p')}))
    first=[cell(a,b) for a,b in ((['sc_half'],[1,4]),([1,4],[1,2]),([1,2],[3,4]),([3,4],[1,1]))]
    second=[cell(a,b) for a,b in (([1,1],[5,4]),([5,4],[3,2]),([3,2],[7,4]),([7,4],[2,1]))]
    fake.ranges=[dict(exact_Z_cell=['0','.5'],actual_17_chart_C3_transports=[dict(chart='first_micro',cells=first),dict(chart='second_micro',cells=second)])]
    for text in ('1','1.5','2'):
        local,prefix,_=fake._covers('.3',iv.mpf(text));_,_,H=fake._correction_logs(local,prefix,False)
        total=sum((p[1] for p in prefix),iv.mpf(0))
        assert _bounds(c,total-(iv.mpf(text)-iv.mpf(1)/16))[1]>=0;checks+=1
        assert len(prefix)>=4
        for name,rate in (('m',1),('h',c.mpf('1.5')),('k',c.mpf('1.5')),('e',1),('p',0)):
            expected=c.quad(lambda t:c.exp(-rate*(c.mpf(text)-t)/7)*c.mpf('1.39')*c.cos(17*(t-c.mpf(1)/16)/7)/119,
                            [c.mpf(1)/16,1,c.mpf(text)])
            assert abs(expected)<=c.exp(_bounds(c,H[name][0])[1]);checks+=1
            if text=='1':assert expected!=0
    return checks


def overlap(a,b):
    c=mp.mp.clone();c.dps=180
    assert abs(c.mpf(a['approximation'])-c.mpf(b['approximation']))<=c.mpf(a['absolute_error_upper'])+c.mpf(b['absolute_error_upper'])


def run():
    began=time.perf_counter();independent=independent_control_and_memory_checks()
    print('Independent coupled-control, third-exponential and incoming-memory checks pass',flush=True)
    field=OriginalSecondMicroNumeric();samples=[];proofs=history_checks=0
    for Z in ('-.3','0','.3'):
        for phase in ('1','1.5','2'):
            p=field.evaluate(Z,phase,log_tau='-2',azimuth='.25')
            v=p['normalized_modified_leading']
            assert [len(v[k]) for k in ('Phi','Psi','Mpsi','Qpsi')]==[4,4,3,2]
            assert p['source_cover']['cumulative_history_charts']==['first_micro','second_micro']
            assert p['source_cover']['prefix_source_cells']>=4
            assert p['limits']['actual_velocity_Z3_numeric'] and not p['limits']['actual_mean_Z3_numeric']
            assert not p['limits']['Qpsi_ZZ_numeric'] and not p['physical']['physical_gradient_or_vorticity_available']
            third=p['normalized_third_derivative_log_error_bounds']
            assert mp.make_mpf(tuple(third['delta_Psi']['log_absolute_upper']['upper_exact_mpf_tuple']))<0
            assert mp.make_mpf(tuple(third['delta_Mpsi']['log_absolute_upper']['lower_exact_mpf_tuple']))>0
            for row in p['small_error_proofs']:assert row['numerical_absolute_error_cap']=='1e-1000';proofs+=1
            samples.append(p)
        print('Actual selected-N second-micro values/Z3: Z='+Z,flush=True)
    samples.append(field.evaluate('.17','1.37',log_tau='-3',azimuth='1.1'))
    reference=json.loads((DIRECTORY/'lei_ren_part1_paper_compliant_current_original_first_micro_numeric.json').read_text(encoding='utf-8'))
    assert reference['checks']['passed'] and reference['source_family']==field.family
    seams=0
    for Z in ('-.3','0','.3'):
        first=next(p for p in reference['samples'] if p['Z']==Z and p['coordinate']==dict(kind='exact_decimal_phase',value='1'))
        second=next(p for p in samples if p['Z']==Z and p['coordinate']==dict(kind='exact_decimal_phase',value='1'))
        for name,rows in first['normalized_modified_leading'].items():
            for a,b in zip(rows,second['normalized_modified_leading'][name]):overlap(a,b);seams+=1
        for name,rows in first['physical']['cylindrical'].items():
            for a,b in zip(rows,second['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);seams+=1
        assert first['original_five_correction_history_bounds']==second['original_five_correction_history_bounds'];history_checks+=20
        exit=next(p for p in samples if p['Z']==Z and p['coordinate']==dict(kind='exact_decimal_phase',value='2'))
        assert exit['physical']['n1_own_history']['actual_terminal_primitive']==first['physical']['n1_own_history']['actual_terminal_primitive'];history_checks+=1
    coarse=OriginalSecondMicroNumeric(precision=60,amplitude_precision=90)
    old=coarse.evaluate('.3','2',log_tau='-2',azimuth='.25')
    new=next(p for p in samples if p['Z']=='.3' and p['coordinate']==dict(kind='exact_decimal_phase',value='2'))
    comparisons=0
    for name,rows in new['normalized_modified_leading'].items():
        for a,b in zip(rows,old['normalized_modified_leading'][name]):overlap(a,b);comparisons+=1
    for name,rows in new['physical']['cylindrical'].items():
        for a,b in zip(rows,old['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);comparisons+=1
    rejected=0
    for Z,s,kwargs in (('1','1.5',{}),('.3','.99',{}),('.3','2.01',{}),('.3','nan',{}),('.3','2',{'log_tau':'1'}),('.3','2',{'azimuth':'nan'})):
        try:field.evaluate(Z,s,**kwargs)
        except (ValueError,ArithmeticError):rejected+=1
        else:raise AssertionError('Unsupported second-micro request admitted')
    try:field.selected('.3','.5')
    except ValueError:rejected+=1
    else:raise AssertionError('First-chart selected coordinate admitted on second chart')
    code=dict(field.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=field.family,checks=dict(passed=True,
        independent_coupled_control_third_exponential_and_memory_comparisons=independent,
        actual_second_micro_packets=len(samples),selected_N_small_correction_error_proofs=proofs,
        first_second_native_and_physical_seam_overlaps=seams,all_five_incoming_and_n1_terminal_history_checks=history_checks,
        local_velocity_Z3_available_while_incoming_mean_Z3_unresolved=True,
        cross_precision_native_physical_coefficient_overlaps=comparisons,rejected_unsupported_requests=rejected),
        samples=samples,input_hashes=field.hashes,implementation_hashes=code,
        limitations=dict(Phi_Psi_Z_order=3,own_mean_Z_order=2,Qpsi_Z_order=1,
            total_physical_radial_derivatives=False,resolved_fast_phase_point=False,
            resolved_signed_correction_integrals=False,mean_Z3_and_QZZ_numeric=False,
            complete_Cartesian_C2_join=False,global_N_admission=False,global_field=False,
            genuine_temporal_n_ge_2=False,full_corrected_NS=False),execution_seconds=round(time.perf_counter()-began,3))


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
