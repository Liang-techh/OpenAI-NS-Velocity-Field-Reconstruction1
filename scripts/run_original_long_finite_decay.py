"""Check retained finite-long histories at formerly rejected ordinary distances."""
import argparse
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
from openai_ns_reconstruction.original_long_numeric import OriginalLongNumeric
from openai_ns_reconstruction.original_long_finite_decay import OriginalLongFiniteDecayNumeric,OriginalLongFiniteDecayVelocity
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_long_finite_decay.json'


def independent_scaled_history_checks():
    """Homogeneous ODE invariants after canceling enormous source factors."""
    from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
    import lei_ren_part1_paper_compliant_long_reshape_profiles as kernels
    c=mp.mp.clone();c.dps=65;iv=MPIntervalContext();iv.dps=115
    f=object.__new__(OriginalLongFiniteDecayNumeric);f.ctx=c;f.iv=iv;f.J=J;f.kernels=kernels;f.T=iv.mpf('1e100')
    f.leading=SimpleNamespace(_abs=lambda v:abs(v),_jet_summary=lambda p,t:[],
        _symmetric=lambda v:iv.mpf([-_bounds(c,v)[1],_bounds(c,v)[1]]))
    f.owner=SimpleNamespace(_cap=lambda *args:(_ for _ in ()).throw(AssertionError('No direct-branch absolute cap expected')))
    z=J.variable(c,c.mpf('.3'),2);zi=J.variable(iv,iv.mpf('.3'),2)
    scales=(c.exp(c.mpf('2e18')),iv.exp(iv.mpf('2e18')))
    def inputs(ctx,Z,scale):
        zero=J.constant(ctx,0,2);phi=1+Z*Z
        return dict(values=dict(Phi=phi,Psi=zero),B=zero,
            six=dict(H=phi,Mpsi=(1+Z*Z)*scale,Kpsi=(2-Z+Z*Z)*phi*scale,
                Apsi=(3+Z+Z*Z)*scale*scale,B=phi*phi/2,C=phi*phi))
    inp=[inputs(c,z,scales[0]),inputs(iv,zi,scales[1])];zero=(J.constant(c,0,2),J.constant(iv,0,2))
    expected=dict(Mpsi=1+zi*zi,Apsi=3+zi+zi*zi,Kpsi_over_Phi=2-zi+zi*zi)
    checks=0
    for raw in ('9999','10000','10001','1e18','1e30'):
        yp,y=c.mpf(raw),iv.mpf(raw)
        p,t,_,_=f._six(inp,zero,yp/c.mpf('1e100'),y/f.T,yp,y,[])
        for name,rate,scale_log in (('Mpsi',1,'2e18'),('Apsi',1,'4e18'),('Kpsi_over_Phi','1.6','2e18')):
            factor=iv.exp(iv.mpf(rate)*y-iv.mpf(scale_log))
            actual=t[name]*factor
            for k in range(3):
                lo,hi=_bounds(c,actual[k]);ll,hh=_bounds(c,expected[name][k])
                assert max(lo,ll)<=min(hi,hh),(raw,name,k);checks+=1
            assert p[name][0]!=0;checks+=1
    return checks


def run():
    started=time.perf_counter();independent=independent_scaled_history_checks()
    print('Independent scaled homogeneous-history invariants and nonzero centers pass',flush=True)
    f=OriginalLongFiniteDecayVelocity();samples=[];proofs=physical=0
    ys=('0','1','9999','10000','10001','1e17','1e18','1e30','1e31')
    for Z in ('-.3','0','.3'):
        for y in ys:
            p=f.at_log_radius(Z,y,log_tau='-2',azimuth='.25');samples.append(p)
            assert p['finite_decay_numeric_extension']['old_absolute_cap_not_weakened']
            assert p['finite_decay_numeric_extension']['direct_finite_exponential_used']==(f.ctx.mpf(y)<=f.ctx.mpf('1e30'))
            assert p['physical']['n1_own_history']['y_cancellation_before_enclosure'];physical+=1
            assert not p['limits']['global_field'] and not p['limits']['genuine_temporal_n_ge_2']
            proofs+=sum(isinstance(row,dict) and row.get('numerical_absolute_error_cap')=='1e-1000' for row in p['small_error_proofs'])
        print('Actual long ordinary/log-distance histories: Z='+Z,flush=True)
    samples.extend(f.evaluate('.3',q,log_tau='-2',azimuth='.25') for q in ('.5','1'))
    samples.append(f.at_log_radius('.3','1e18',log_tau=None,scaled_time='1e18',azimuth='1.7'))
    # Reproduce the old rejection on the identical already-loaded actual
    # source family. No new legacy constructor or historical suite is run.
    legacy=object.__new__(OriginalLongNumeric);legacy.__dict__=dict(f.numeric.__dict__)
    try:legacy.at_log_radius('.3','10000')
    except ArithmeticError as error:
        assert 'Mpsi_inlet_memory_Z0' in str(error);old_failure=str(error)
    else:raise AssertionError('Original finite-distance rejection no longer reproduced')
    old=json.loads((DIRECTORY/(PREFIX+'compliant_current_original_long_velocity.json')).read_text(encoding='utf-8'))
    seams=history=0
    for Z in ('-.3','0','.3'):
        before=next(p for p in old['samples'] if p['Z']==Z and p['coordinate']['value']=='0')
        after=next(p for p in samples if p['Z']==Z and p['coordinate']['value']=='0')
        for group in ('normalized_modified_leading','unmodified_six_current_amplitude_units'):
            for name,rows in before[group].items():
                for a,b in zip(rows,after[group][name]):overlap(a,b);seams+=1
        assert before['selected_integer']==after['selected_integer'];history+=1
        assert before['original_five_correction_history_bounds']==after['original_five_correction_history_bounds'];history+=1
        assert before['physical']['n1_own_history']['actual_terminal_primitive']==after['physical']['n1_own_history']['actual_terminal_primitive'];history+=1
    for p in samples:
        assert p['selected_integer']==samples[0]['selected_integer'];history+=1
        if p['coordinate']['kind']=='exact_decimal_q':
            before=next(q for q in old['samples'] if q['Z']=='.3' and q['coordinate']['kind']=='exact_decimal_q' and q['coordinate']['value']==p['coordinate']['value'])
            for n,rows in p['unmodified_six_current_amplitude_units'].items():
                for a,b in zip(rows,before['unmodified_six_current_amplitude_units'][n]):overlap(a,b);seams+=1
    coarse=OriginalLongFiniteDecayVelocity(precision=60,amplitude_precision=90)
    before=coarse.at_log_radius('.3','1e18',log_tau='-2',azimuth='.25')
    after=next(p for p in samples if p['Z']=='.3' and p['coordinate']['value']=='1e18' and p['physical']['time']['kind']=='log_tau')
    comparisons=0
    for group in ('normalized_modified_leading','unmodified_six_current_amplitude_units'):
        for n,rows in before[group].items():
            for a,b in zip(rows,after[group][n]):overlap(a,b);comparisons+=1
    for n,rows in before['physical']['cylindrical'].items():
        for a,b in zip(rows,after['physical']['cylindrical'][n]):
            overlap(a['coefficient'],b['coefficient']);comparisons+=1
            for aa,bb in zip(a['log_scale_terms'],b['log_scale_terms']):overlap(aa,bb);comparisons+=1
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,checks=dict(passed=True,independent_scaled_homogeneous_history_invariants=independent,
        actual_original_finite_decay_packets=len(samples),legacy_actual_source_failure_reproduced=True,
        direct_and_combined_log_history_packets=physical,selected_N_small_error_proofs=proofs,
        accepted_inlet_and_large_window_seam_coefficients=seams,original_integer_history_and_primitive_checks=history,
        cross_precision_coefficient_and_log_comparisons=comparisons),
        legacy_rejection=old_failure,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        execution_seconds=round(time.perf_counter()-started,3),
        limits=dict(native_coordinates_required=True,ordinary_long_decay_range_extended=True,
            global_N_admission=False,radial_or_fast_phase_derivatives=False,pressure_exported=False,
            actual_patch_connected=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    r=json.loads(path.read_text(encoding='utf-8'));assert r['checks']['passed']
    for n,d in r['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==d,n
    for n,d in r['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==d,n
    return r


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    result=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],execution_seconds=result['execution_seconds']),indent=2),flush=True)
