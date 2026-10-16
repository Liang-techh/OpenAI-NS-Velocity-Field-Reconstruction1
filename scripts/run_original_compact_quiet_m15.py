"""Evaluate actual compact C3 pressure and exact quiet m15 continuation to Rp."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_compact_pressure_rows import OriginalCompactPressureRows
from openai_ns_reconstruction.original_compact_quiet_m15 import OriginalCompactQuietM15
from openai_ns_reconstruction.original_Rm_finite_N_numeric import RATES
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_compact_quiet_m15_numeric.json'


def independent_fixtures():
    c=mp.mp.clone();c.dps=100;iv=MPIntervalContext();iv.dps=130
    f=object.__new__(OriginalCompactPressureRows);f.ctx,f.iv,f.audit,f.J=c,iv,c,J
    f.chart='compact';f._histories={};f.L=iv.log(2)
    f.supports=[(f.L*(iv.mpf(v)-iv.mpf(1)/40),f.L*(iv.mpf(v)+iv.mpf(1)/40)) for v in ('.2','.5','.8')]
    f._matched=lambda Z:[dict(exact_Z_cell=['-1','1'])];f.geometry=lambda x:dict(local=iv.log(iv.mpf(x)))
    f.certificates={n:dict(synthetic_terminal_identity=True) for n in RATES}
    incoming={n:[iv.log(iv.mpf('1.25')*(i+1)*(k+1)/17) for k in range(4)] for i,n in enumerate(RATES)}
    f.incoming_history_logs=lambda Z:incoming
    width=c.log(2)/20;left=c.log(2)*(c.mpf('.2')-c.mpf(1)/40)
    index={n:i for i,n in enumerate(RATES)}
    def cap(frames,group,name,k):
        raw=iv.mpf((index[name]+1)*(k+1))/17
        return iv.log(raw*(iv.mpf('1.1') if group=='actual_incoming_correction_C3_bounds' else 2/(iv.log(2)/20)))
    f.cap=cap;transport=0;terminal=0
    for x in ('1','1.13','1.15','1.2','1.8','2'):
        out=f.history_logs('0',x);s=c.log(c.mpf(x));t=max(c.mpf(0),min(c.mpf(1),(s-left)/width))
        switch=6*t**5-15*t**4+10*t**3
        for i,(name,rate) in enumerate(RATES.items()):
            for k in range(4):
                actual=c.mpf((i+1)*(k+1))/17*c.exp(-c.mpf(rate)*s)*(1-switch)
                v=out['history'][name][k];bound=c.mpf(0) if v is None else c.exp(_bounds(c,v)[1])
                sup=c.exp(_bounds(c,out['prefix_history'][name][k])[1]);assert bound>=abs(actual) and sup>=bound
                transport+=1
                if x=='2':assert v is None and sup>0;terminal+=1
    # An actual mean perturbation and its U=M_y+M velocity generate Q/Qy.
    # Bound those ordinary derivatives via the production C3 recovery.
    f.o=SimpleNamespace(core=SimpleNamespace(delta_iv=iv.mpf(1)/17));pressure=0
    s=left+width*c.mpf('.35');t=(s-left)/width
    sw=6*t**5-15*t**4+10*t**3
    swy=(30*t**4-60*t**3+30*t**2)/width
    swyy=(120*t**3-180*t**2+60*t)/width**2
    for Z in ('-.3','0','.23'):
        z=J.variable(c,c.mpf(Z),3);de=c.mpf(1)/17;A=1-z*z*de;D=1-z*z
        amp=c.mpf('.13')+z*c.mpf('.07')+z*z*c.mpf('.02')+z**3*c.mpf('.01')
        M=amp*(c.exp(-s)*(1-sw));U=amp*(-c.exp(-s)*swy);Uy=amp*(c.exp(-s)*(swy-swyy))
        dz=lambda v:J(c,[(k+1)*v[k+1] for k in range(v.order)])
        Q=((2*z*U-z*(1-de)*M-D*dz(M))/A).truncate(2)
        Qy=((2*z*Uy-z*(1-de)*(U-M)-D*(dz(U)-dz(M)))/A).truncate(1)
        def logs(v,order):
            return [None if v[k]==0 else iv.log(iv.mpf(abs(v[k]*math.factorial(k)))*iv.mpf('1.001')) for k in range(order+1)]
        errors=f.finish_errors(Z,dict(U=logs(U,3),m=logs(M,3),U_y=logs(Uy,1)))
        assert errors['m'][3] is not None
        for name,field in (('Q',Q),('Q_y',Qy)):
            for k in range(field.order+1):
                bound=c.exp(_bounds(c,errors[name][k])[1]);assert bound>=abs(field[k]*math.factorial(k));pressure+=1
    return dict(passed=True,independent_supported_C3_Duhamel_rows=transport,
        exact_terminal_rows_with_positive_prefix=terminal,independent_actual_QZZ_QyZ_error_rows=pressure,
        finite_fixture_N=17,already_N_scaled_source_not_divided_again=True)


def run():
    start=time.perf_counter();fixture=independent_fixtures();print('Independent supported C3 transport and physical QZZ/QyZ fixtures pass',flush=True)
    compact=OriginalCompactQuietM15(precision=85);quiet=compact.next_field();engines=(compact,quiet)
    a=mp.mp.clone();a.dps=190;samples=[];pressure=[]
    bounds=lambda row:tuple(a.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples'])
    highlog=lambda row:-a.inf if row['exact_zero'] else a.make_mpf(tuple(row['log_absolute_upper']['upper_exact_mpf_tuple']))
    counts=dict(exact_inlet_replays=0,serialized_increment_rows=0,exact_empty_rows=0,pressure_rows=0,
        C3_history_rows=0,terminal_relative_zero_rows=0,quiet_exact_error_rows=0,positive_absolute_mean_conditions=0,
        prefix_error_conditions=0,nonzero_terminal_prefix_conditions=0,original_archive_overlaps=0,
        compact_profile_y_units_conditions=0,physical_width_conditions=0)
    archives={}
    for chart,suffix,key in (('compact','compact_numeric','native_compact_x'),('quiet','quiet_numeric','native_quiet_fraction')):
        path=DIRECTORY/(PREFIX+'compliant_current_original_'+suffix+'.json');raw=path.read_bytes();old=json.loads(raw)
        assert old['checks']['passed'] and old['source_family']==compact.family
        archives[chart]={(p['Z'],p[key]):p for p in old['samples']}
        for engine in engines:
            engine.hashes[path.name]=hashlib.sha256(raw).hexdigest()
            for name,digest in old['implementation_hashes'].items():
                assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
                assert engine.code_hashes.get(name,digest)==digest;engine.code_hashes[name]=digest
    for Z in ('-.3','0','.3'):
        values=('1','1.15','1.5','2') if Z=='0' else ('1','1.15','2')
        for engine,coordinates in ((compact,values),(quiet,('0','.5','1'))):
            for value in coordinates:
                out=engine.evaluate(Z,value);p=engine.rows.evaluate(Z,value);samples.append(out);pressure.append(p)
                count=engine.COUNTS[engine.chart];extra=4 if engine.chart=='compact' else 1
                assert len(out['m15_C0_C1_terms'])==count+extra
                assert out['m15_C0_C1_terms'][:count]==out['incoming_m15_terms']==engine.inlet(Z)['m15_C0_C1_terms'];counts['exact_inlet_replays']+=1
                assert all(out['contract'].values()) and all(p['contract'].values())
                for term in out['m15_C0_C1_terms'][count:]:
                    for row in term['coefficient_C0_C1']:
                        lo,hi=bounds(row);point=a.mpf(row['approximation']);assert max(abs(lo-point),abs(hi-point))<=a.mpf(row['absolute_error_upper']);counts['serialized_increment_rows']+=1
                        if out['exact_empty_increment']:assert lo==hi==point==0;counts['exact_empty_rows']+=1
                pp,tt,geo=engine.rows.background(Z,value)
                for name,field in p['physical_pressure_rows'].items():
                    for k,row in enumerate(field['signed_background_ordinary_Z']):
                        assert row['source_enclosure_exact_mpf_tuples']==[list(v) for v in (tt[name][k]*math.factorial(k))._mpi_];counts['pressure_rows']+=1
                assert _bounds(a,engine.rows.scalars(value)['truth']['mean'])[0]>0 and _bounds(a,tt['Q_y'][0])[0]>0
                counts['positive_absolute_mean_conditions']+=1
                terminal=engine.chart=='compact' and value=='2'
                for name,rows in p['finite_N_history_error_logs'].items():
                    assert len(rows)==4;counts['C3_history_rows']+=4
                    if terminal:assert all(row['exact_zero'] for row in rows);counts['terminal_relative_zero_rows']+=4
                    if engine.chart=='quiet':assert all(row['exact_zero'] for row in rows);counts['quiet_exact_error_rows']+=4
                end,_=engine.rows.error_data(Z,value);prefix,_=engine.rows.error_data(Z,value,prefix=True)
                for name in ('U','m','Q','U_y','Q_y'):
                    for v,w in zip(engine.rows.log_rows(prefix[name]),engine.rows.log_rows(end[name])):
                        assert highlog(v)>=highlog(w);counts['prefix_error_conditions']+=1
                        if engine.chart=='quiet':assert v['exact_zero'] and w['exact_zero']
                if terminal:
                    assert all(v is None for name in ('U','m','Q','U_y','Q_y') for v in end[name])
                    assert all(v is not None for v in prefix['m']);counts['nonzero_terminal_prefix_conditions']+=4
                if engine.chart=='compact':
                    assert _bounds(a,geo['local'])[0]<=a.log(a.mpf(value))<=_bounds(a,geo['local'])[1];counts['physical_width_conditions']+=1
                    h=engine.rows.history_logs(Z,value);active=bool(set(h['local_bumps'])&{0,2})
                    component=p['finite_N_axial_radial_error_components']['actual_physical_deltaV_y_C1']
                    for k,row in enumerate(component):
                        if active:
                            source=engine.rows.cover_owner._log_cover([fr['actual_selected_N_profile_y0_y1_y2_C3_bounds'][1]['delta_V'][k] for fr in h['source_frames']])
                            assert highlog(row)==_bounds(a,source+engine.rows.o.core.lp_iv)[1]
                        else:assert row['exact_zero']
                        counts['compact_profile_y_units_conditions']+=1
                else:
                    want=(-60*(a.log(a.mpf('.001'))-4*(a.exp(40)+11))-2-a.log(2))*a.mpf(value)
                    assert _bounds(a,geo['local'])[0]<=want<=_bounds(a,geo['local'])[1];counts['physical_width_conditions']+=1
                    assert out['incoming_m15_terms']==compact.evaluate(Z,'2')['m15_C0_C1_terms']
                old=archives[engine.chart].get((Z,value))
                if old:
                    for name,oldname in (('U','Uz'),('m','Mz_over_R'),('Q','Ur_over_sqrt_R_over2')):
                        for k,row in enumerate(old['physical_similarity'][oldname]):
                            lo,hi=bounds(row);ll,hh=_bounds(a,tt[name][k]);ll*=math.factorial(k);hh*=math.factorial(k)
                            assert max(lo,ll)<=min(hi,hh),(Z,engine.chart,value,name,k);counts['original_archive_overlaps']+=1
            print('Original '+engine.chart+' pressure/m15: Z='+Z+', '+str(count)+' ->'+str(count+extra)+' terms',flush=True)
    rejected=0
    for engine,Z,value in ((compact,'1','1'),(compact,'0','.9'),(compact,'0','2.1'),(quiet,'0','-.1'),(quiet,'0','1.1')):
        try:engine.evaluate(Z,value)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid compact/quiet query admitted')
    hashes={};code={}
    for engine in engines:hashes.update(engine.hashes);code.update(engine.code_hashes)
    name='scripts/run_original_compact_quiet_m15.py';code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=compact.family,samples=samples,pressure_row_samples=pressure,input_hashes=hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_fixtures=fixture,actual_moment_packets=len(samples),actual_pressure_packets=len(pressure),
            **counts,rejected_invalid_queries=rejected),execution_seconds=round(time.perf_counter()-start,3),limits=dict(samples[-1]['limits']))


def inspect(path):
    data=json.loads(path.read_text(encoding='utf-8'));assert data['checks']['passed']
    for base,group in ((DIRECTORY,data['input_hashes']),(ROOT,data['implementation_hashes'])):
        for name,digest in group.items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;data=inspect(path) if args.inspect_report else run()
    if args.write_report:
        temp=path.with_suffix('.json.tmp');temp.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');temp.replace(path)
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
