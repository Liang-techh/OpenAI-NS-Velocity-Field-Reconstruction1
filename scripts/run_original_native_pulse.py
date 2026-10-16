"""Actual Rp callback and original selected main/gap/end physical velocity."""
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
from openai_ns_reconstruction.original_native_pulse_numeric import OriginalNativePulseVelocity
from lei_ren_part1_paper_compliant_axial_pulse_field import CompliantAxialPulseField,gp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_native_pulse_numeric.json'


def independent_shapes_and_flows():
    c=mp.mp.clone();c.dps=85;iv=MPIntervalContext();iv.dps=85;checks=0
    def sigma(x):
        if x<=0:return c.mpf(0)
        if x>=1:return c.mpf(1)
        odds=-1/x**2+1/(1-x)**2
        return 1/(1+c.exp(-odds))
    def primitive(x):
        if x>=c.mpf('.02'):return x-c.mpf('.01')
        return x*c.quad(lambda v:sigma(50*x*v),[0,c.mpf('.5'),1])
    for text in ('0','.001','.01','.02','.5','10.5','11'):
        x=c.mpf(text);truth=primitive(x)*sigma(11-x)
        got=gp(iv,iv.mpf(text))['value'];lo,hi=_bounds(c,got)
        assert lo-c.mpf('1e-75')<=truth<=hi+c.mpf('1e-75'),text;checks+=1
    f=object.__new__(CompliantAxialPulseField);f.ctx=iv;f.mu=iv.mpf('.02');f.delta=iv.mpf('.001')
    f.rate=1-f.mu;f.prate=1+2*f.mu;f.Xp=iv.mpf('.9')
    selected=dict(selected_ap_Taylor=J.constant(iv,iv.mpf('1.01'),1))
    inlet=dict(Mp_over_Pstar_squared=[iv.mpf('.4'),iv.mpf(0)])
    f.data=lambda Z:(selected,inlet,J.constant(iv,1,1),J.constant(iv,iv.mpf('.1'),1),
                     [J.constant(iv,iv.mpf(v),1) for v in ('.2','-.3')])
    datum=SimpleNamespace(normalized_jets=lambda Z,n:dict(normalized_pressure_coefficients=[iv.mpf('.7'),iv.mpf(0)]))
    f.selection=SimpleNamespace(K=iv.mpf(1),future=SimpleNamespace(angular=SimpleNamespace(initial=SimpleNamespace(datum=datum))))
    for text in ('.5','1'):
        xi=c.mpf(text);t=xi/c.mpf('.02');packet=f.main('0',text,cells=64)
        for row,initial in ((1,c.mpf('.2')),(2,c.mpf('-.3'))):
            lam=c.mpf('.5')-row*c.mpf('.02')
            # Fubini reduces the independent convolution to one sigma
            # integral; it does not use the producer's rectangle kernels.
            start=c.quad(lambda r:c.exp(lam*r/c.mpf('.02'))*sigma(50*r),[0,c.mpf('.01'),c.mpf('.02')])
            bulk=c.mpf('.02')/lam*(c.exp(lam*t)-c.exp(lam))
            kernel=(primitive(xi)-c.exp(-lam*t)*(start+bulk))/lam
            truth=initial*c.exp(-lam*t)+c.mpf('1.01')*kernel
            name='Mz_over_R_Utheta' if row==1 else 'Mtheta_z_over_sqrt2_R_3half_Utheta_squared'
            lo,hi=_bounds(c,packet[name][0]);assert lo<=truth<=hi,(text,row);checks+=1
    return checks


def run():
    start=time.perf_counter();independent=independent_shapes_and_flows()
    print('Independent original gp and Fubini convolution checks pass',flush=True)
    f=OriginalNativePulseVelocity(cells=128);c,iv=f.ctx,f.iv;samples=[]
    callbacks=seams=nonzero=history=terminal=physical=C5=0
    points=[('entrance','0'),('entrance','1'),('main','1'),('main','10.5'),('main','11'),
            ('gap','12'),('gap_end','-1000'),('end','-4'),('end','-3'),('end','-2'),('end','-1'),('end','0')]
    for Z in ('-.3','0','.3'):
        for chart,x in points:
            p=f.evaluate(Z,chart,x,log_tau='-2',azimuth='.25');samples.append(p)
            assert p['native_selected_inlet_C5_source_callback'];callbacks+=1
            assert p['selected_integer']==f.numeric.inlets[Z]['selected_integer'];callbacks+=1
            assert p['physical']['pressure_orders_included']==[0];physical+=1
            assert not p['physical']['full_n1_pressure_included'];physical+=1
            assert p['physical']['interior_P1_not_substituted_by_zero'];physical+=1
            assert p['physical']['radial_half_power_cancelled_before_boxing'];physical+=1
            assert p['physical']['leading_pressure_includes_same_P0_once'];physical+=1
            if chart=='entrance' and x=='0':
                data=f.numeric._data(Z);incoming=f.numeric.inlets[Z]['native_pulse_input_frame']
                assert set(p['native_callback_frame_C5'])==set(incoming)
                for name,rows in p['native_callback_frame_C5'].items():
                    assert len(rows)==6;C5+=1
                    for a,b in zip(rows,incoming[name]):overlap(a,b);seams+=1
                for name,jet in (('u',data[2]),('energy',data[3]),('m1',data[4][0]),('m2',data[4][1])):
                    rows=f.numeric._rows(jet)
                    for a,b in zip(rows,incoming[name]):overlap(a,b);seams+=1
                for name in ('Uz_over_Utheta',):
                    assert _bounds(c,f._restore(p['normalized_fields'][name][0]))==(0,0);seams+=1
            if chart=='main' and x=='1':
                assert _bounds(c,f._restore(p['normalized_fields']['Uz_over_Utheta'][0]))[0]>0;nonzero+=1
            if chart=='end' and x in ('-3','-1'):
                lo,hi=_bounds(c,f._restore(p['normalized_fields']['Uz_over_Utheta'][0]));assert (hi<0 if x=='-3' else lo>0);nonzero+=1
                assert len(p['end_source_log_terms'])==4;physical+=1
            if chart in ('gap','gap_end'):
                assert _bounds(c,f._restore(p['normalized_fields']['Uz_over_Utheta'][0]))==(0,0);history+=1
                lo,hi=_bounds(c,f._restore(p['normalized_fields']['Mz_over_R_Utheta'][0]));assert lo<hi and (lo!=0 or hi!=0);history+=1
            if chart=='end' and x=='0':
                assert p['terminal_linear_zero_by_source_identity'];terminal+=1
                for name in ('Uz_over_Utheta','Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared','Q_over_Utheta'):
                    for row in p['normalized_fields'][name]:assert _bounds(c,f._restore(row))==(0,0);terminal+=1
            assert p['incoming_P1_at_Rp']==f.numeric.inlets[Z]['reserved_n1_pressure'];history+=1
        print('Actual selected pulse reaches Rv: Z='+Z,flush=True)
    scaled=0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'end','-3',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        normal=next(v for v in samples if v['Z']==Z and v['native_pulse_chart']=='end' and v['native_pulse_coordinate']=='-3')
        for group in ('ur','utheta','uz'):
            for a,b in zip(p['physical']['cylindrical'][group],normal['physical']['cylindrical'][group]):overlap(a['coefficient'],b['coefficient']);scaled+=1
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
    rejected=0
    for Z,chart,x in (('1','main','1'),('nan','main','1'),('.3','unknown','0'),('.3','main','12'),('.3','end','1')):
        try:f.numeric.evaluate(Z,chart,x)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original pulse query admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_original_gp_and_linear_convolutions=independent,actual_selected_pulse_physical_packets=len(samples),
            current_native_callback_and_same_N_conditions=callbacks,actual_Rp_inlet_coefficient_seams=seams,
            nonzero_original_main_and_end_velocity_conditions=nonzero,gap_history_and_retained_incoming_P1_conditions=history,
            actual_selected_terminal_linear_zero_conditions=terminal,physical_pressure_scope_and_atomic_end_units=physical,
            scaled_time_conditions=scaled,actual_eight_C5_callback_conditions=C5,rejected_requests=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(original_main_gap_end_velocity_values=True,leading_pressure=True,interior_P1=False,
            signed_coefficients_are_function_enclosures=True,terminal_zeros_from_source_identity_not_independent_signed_residual=True,
            physical_gradients=False,relative_pressure_cancellation_accuracy=False,common_global_N=False,
            exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
