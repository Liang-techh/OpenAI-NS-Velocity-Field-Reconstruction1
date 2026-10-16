"""Evaluate original outer-reference/O2-slope pressure and cumulative m15."""
import argparse
import ast
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
from openai_ns_reconstruction.original_outer_pressure_rows import OriginalOuterPressureRows
from openai_ns_reconstruction.original_outer_m15 import OriginalOuterM15
from openai_ns_reconstruction.original_Rm_finite_N_numeric import RATES
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_outer_m15_numeric.json'


def independent_fixtures():
    c=mp.mp.clone();c.dps=100;iv=MPIntervalContext();iv.dps=130
    leading=0
    for filename in ('original_outer_reference_numeric.py','original_outer_slope_numeric.py'):
        path=ROOT/'src/openai_ns_reconstruction'/filename;tree=ast.parse(path.read_text(encoding='utf-8'))
        node=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='_leading')
        scope={};exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),scope)
        fake=SimpleNamespace(J=J,ctx=c,iv=iv,patch=SimpleNamespace(leading=SimpleNamespace(
            defect=SimpleNamespace(logPp=c.log(17),logP=iv.log(17)))))
        coordinate='-.7' if 'reference' in filename else {name:c.mpf(value) for name,value in zip(('F','H','G','C'),('.8','.7','1.1','3.2'))}
        for Z in ('-.3','0','.23'):
            out=scope['_leading'](fake,Z,coordinate,c)
            for name in ('V','m'):
                want=4*J.variable(c,c.mpf(Z),2)
                for n in range(3):assert abs(out[name][n]*17-want[n])<c.mpf('1e-95');leading+=1
    z,de=sy.symbols('Z delta');A,D=1-de*z*z,1-z*z;U=4*z;Q=(2*z*U-(1-de)*z*U-4*D)/A
    direct=((Q+(1-de)*z*sy.diff(Q,z)/2)/A+Q*Q/2+U*(-2*z*Q+D*sy.diff(Q,z))/A)/2
    f=object.__new__(OriginalOuterM15);f.ctx,f.iv,f.J=c,iv,J;f._densities={}
    rows=object.__new__(OriginalOuterPressureRows);rows.ctx,rows.iv,rows.J=c,iv,J;rows.chart='Rh_reference';rows._backgrounds={}
    rows.reference=SimpleNamespace(logCp=c.mpf(0),logC=iv.mpf(0),logPp=c.mpf(0),logP=iv.mpf(0))
    rows.o=SimpleNamespace(core=SimpleNamespace(delta=c.mpf(1)/17,delta_iv=iv.mpf(1)/17));f.rows=rows;f.o=rows.o;f.chart=rows.chart
    physical=0
    for Z in ('-.3','0','.23'):
        base=f.base_density(Z)['point']
        for length in ('.0001','.7','5'):
            L=c.mpf(length)
            for n in (0,1):
                want=c.mpf(str(sy.N(sy.diff(direct,z,n).subs({z:sy.Rational(Z),de:sy.Rational(1,17)}),100)))
                direct_integral=c.quad(lambda s:want*c.exp(2*s),[0,L])
                got=base[n]*c.expm1(2*L)/2
                assert abs(got-direct_integral)<c.mpf('1e-90')*max(1,abs(got));physical+=1
    # Exercise the production four-row log transport with finite, independently
    # integrated sources. This fixture never exponentiates original huge logs.
    r=object.__new__(OriginalOuterPressureRows);r.ctx,r.iv,r.J=c,iv,J;r.chart='Rh_reference';r._histories={}
    r.owner=SimpleNamespace(logN=iv.log(17))
    incoming={name:[iv.log(iv.mpf((i+1)*(n+1))/19) for n in range(4)] for i,name in enumerate(RATES)}
    r.incoming_history_logs=lambda Z:incoming
    intervals=[(Fraction(-5),Fraction(-9,2)),(Fraction(-9,2),Fraction(-7,2)),(Fraction(-7,2),Fraction(0))]
    groups=[]
    for j,(a,b) in enumerate(intervals):
        cell=dict(actual_density_ranges=dict(N_scaled={name:[iv.log(iv.mpf((j+1)*(i+1)*(n+1))/7) for n in range(4)] for i,name in enumerate(RATES)}))
        groups.append((a,b,[cell,cell]))
    r.utility=SimpleNamespace(_cover=lambda Z,q:([['-1','0'],['0','1']],groups,[groups[1][2][0]]))
    r.cover_owner=SimpleNamespace(_log_cover=lambda values:iv.mpf(max(_bounds(c,v)[1] for v in values)))
    query='-4.2';out=r.history_logs('0',query);transport=0
    for i,(name,rate) in enumerate(RATES.items()):
        rr=c.mpf(rate);end=c.mpf(query)
        for n in range(4):
            want=c.mpf((i+1)*(n+1))/19*c.exp(-rr*(end+5))
            for j,(a,b) in enumerate(intervals):
                left=c.mpf(a.numerator)/a.denominator;right=min(end,c.mpf(b.numerator)/b.denominator)
                if right<=left:continue
                driver=c.mpf((j+1)*(i+1)*(n+1))/(7*17)
                want+=c.quad(lambda s:driver*c.exp(-rr*(end-s)),[left,right])
            bound=c.exp(_bounds(c,out['history'][name][n])[1]);sup=c.exp(_bounds(c,out['prefix_history'][name][n])[1])
            assert bound>=want and sup>=bound,(name,n);transport+=1
    return dict(passed=True,original_leading_source_coefficients=leading,independent_physical_C0_C1_integrals=physical,
        independent_positive_C3_Duhamel_rows=transport,finite_fixture_N=17,
        overlapping_Z_covers_not_independent_radial_cells=True,ordinary_error_rows_not_Taylor_divided=True)


def run():
    start=time.perf_counter();fixture=independent_fixtures()
    print('Original leading-source and independent physical/Duhamel fixtures pass',flush=True)
    f=OriginalOuterM15(precision=85);slope=f.next_field();audit=mp.mp.clone();audit.dps=190;samples=[];pressure=[]
    counts=dict(exact_cumulative_inlet_replays=0,serialized_increment_errors=0,exact_empty_rows=0,
        serialized_pressure_rows=0,genuine_C3_history_rows=0,leading_zero_UyQy_conditions=0,
        native_width_conditions=0,positive_kernel_conditions=0,complete_prefix_error_conditions=0,
        lower_finite_N_endpoint_conditions=0,original_saved_pressure_overlaps=0,exact_history_seams=0,
        current_chart_slow_fast_unit_conditions=0)
    def bounds(row):return tuple(audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples'])
    def highlog(row):return audit.make_mpf(tuple(row['log_absolute_upper']['upper_exact_mpf_tuple']))
    archives={}
    for chart,suffix,key in (('Rh_reference','outer_reference_numeric','native_offset'),('O2_slope','outer_slope_numeric','native_slope')):
        path=DIRECTORY/(PREFIX+'compliant_current_original_'+suffix+'.json');raw=path.read_bytes();data=json.loads(raw)
        assert data['checks']['passed'] and data['source_family']==f.family
        archives[chart]={(p['Z'],p[key]):p for p in data['samples']};f.hashes[path.name]=hashlib.sha256(raw).hexdigest();slope.hashes[path.name]=hashlib.sha256(raw).hexdigest()
    for Z in ('-.3','0','.3'):
        for engine,coordinates in ((f,('-5','-4.5','-2.5','0')),(slope,('0','.25','.5','1'))):
            for value in coordinates:
                out=engine.evaluate(Z,value);samples.append(out);rows=engine.rows;p=rows.evaluate(Z,value);pressure.append(p)
                count=168 if engine.chart=='Rh_reference' else 172
                assert len(out['m15_C0_C1_terms'])==count+4
                assert out['m15_C0_C1_terms'][:count]==out['incoming_m15_terms']==engine.inlet(Z)['m15_C0_C1_terms']
                counts['exact_cumulative_inlet_replays']+=1
                assert engine.evaluate(Z,value) is out and all(out['contract'].values()) and all(p['contract'].values())
                for term in out['m15_C0_C1_terms'][count:]:
                    for row in term['coefficient_C0_C1']:
                        lo,hi=bounds(row);point=audit.mpf(row['approximation'])
                        assert max(abs(lo-point),abs(hi-point))<=audit.mpf(row['absolute_error_upper']);counts['serialized_increment_errors']+=1
                        if out['exact_empty_increment']:assert lo==hi==point==0;counts['exact_empty_rows']+=1
                point,truth,geo=rows.background(Z,value)
                for name,field in p['physical_pressure_rows'].items():
                    for n,row in enumerate(field['signed_background_ordinary_Z']):
                        assert row['source_enclosure_exact_mpf_tuples']==[list(v) for v in (truth[name][n]*math.factorial(n))._mpi_]
                        counts['serialized_pressure_rows']+=1
                    if name in ('U_y','Q_y'):
                        for row in field['signed_background_ordinary_Z']:assert bounds(row)==(0,0);counts['leading_zero_UyQy_conditions']+=1
                for name,values in p['finite_N_history_error_logs'].items():
                    assert len(values)==4 and not any(row['exact_zero'] for row in values);counts['genuine_C3_history_rows']+=4
                if out['exact_empty_increment']:
                    assert p['finite_N_history_error_logs']==p['finite_N_incoming_history_error_logs'];counts['exact_history_seams']+=1
                history=rows.history_logs(Z,value);width=sum((step['physical_log_width'] for step in history['steps']),rows.iv.mpf(0))
                a,b=_bounds(audit,width);want=audit.mpf(value)-audit.mpf(rows.DOMAINS[engine.chart][0].numerator)/rows.DOMAINS[engine.chart][0].denominator
                assert a<=want<=b;counts['native_width_conditions']+=1
                seen=set()
                for step in history['steps']:
                    key=str(step['exact_left']),str(step['exact_end']);assert key not in seen;seen.add(key)
                a,b=bounds(out['positive_native_kernel']);assert a>=0 and b>=a;counts['positive_kernel_conditions']+=1
                prefix,_=rows.error_data(Z,value,prefix=True);endpoint,_=rows.error_data(Z,value)
                for name in ('U','m','Q','U_y','Q_y'):
                    for a,b in zip(rows.log_rows(prefix[name]),rows.log_rows(endpoint[name])):
                        assert highlog(a)>=highlog(b),(Z,engine.chart,value,name);counts['complete_prefix_error_conditions']+=1
                if not out['exact_empty_increment']:
                    lower=out['actual_lower_finite_N_endpoint_budget']
                    for kind in ('endpoint','viscous'):
                        for row in lower[kind]:assert not row['exact_zero'];counts['lower_finite_N_endpoint_conditions']+=1
                old=archives[engine.chart].get((Z,value))
                if old is not None:
                    for name,oldname in (('U','Uz'),('m','Mz_over_R'),('Q','Ur_over_sqrt_R_over2')):
                        for n,row in enumerate(old['physical_similarity'][oldname]):
                            a,b=bounds(row);c,d=_bounds(audit,truth[name][n]);c*=math.factorial(n);d*=math.factorial(n)
                            assert max(a,c)<=min(b,d),(Z,engine.chart,value,name,n);counts['original_saved_pressure_overlaps']+=1
                for n,key in enumerate(('B_phi','B_Z_phi')):
                    source=rows.read_log(rows.loop['phase_and_first_slow_phase_log_bounds'][key])
                    value_row=p['finite_N_axial_radial_error_components']['fast_B_phi'][n]
                    if source is None:assert value_row['exact_zero']
                    else:assert highlog(value_row)==_bounds(audit,source+rows.o.core.lp_iv)[1]
                    counts['current_chart_slow_fast_unit_conditions']+=1
            print('Original '+engine.chart+' pressure/m15: Z='+Z+', '+str(count)+' ->'+str(count+4)+' terms',flush=True)
    rejected=0
    for engine,Z,value in ((f,'1','-5'),(f,'0','-6'),(f,'0','.1'),(slope,'0','-1'),(slope,'0','1.1')):
        try:engine.evaluate(Z,value)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid outer pressure query admitted')
    hashes={**f.hashes,**slope.hashes};code={**f.code_hashes,**slope.code_hashes}
    name='scripts/run_original_outer_m15.py';code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,pressure_row_samples=pressure,input_hashes=hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_source_and_transport=fixture,actual_moment_packets=len(samples),
            actual_pressure_packets=len(pressure),**counts,rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(original_cumulative_signed_m15_to_Rref_and_O2_slope_with_complete_absolute_errors=True,
            resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
            O2_axial_changing_source_pressure_integrated=False,independent_preRp_signed_seam=False,
            global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
            stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    data=json.loads(path.read_text(encoding='utf-8'));assert data['checks']['passed']
    for base,group in ((DIRECTORY,data['input_hashes']),(ROOT,data['implementation_hashes'])):
        for name,digest in group.items():assert hashlib.sha256((base/name).read_bytes()).hexdigest()==digest,name
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;data=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=data['checks'],execution_seconds=data['execution_seconds']),indent=2),flush=True)
