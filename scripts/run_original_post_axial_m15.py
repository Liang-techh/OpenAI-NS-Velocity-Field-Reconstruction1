"""Evaluate original pressure/m15 through O2 buffer and O3 to Rc."""
import argparse
import ast
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
from openai_ns_reconstruction.original_post_axial_m15 import OriginalPostAxialM15
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_post_axial_m15_numeric.json'


def independent_fixtures():
    c=mp.mp.clone();c.dps=100;iv=MPIntervalContext();iv.dps=130
    methods={}
    for filename in ('original_outer_buffer_numeric.py','original_outer_o3_numeric.py'):
        path=ROOT/'src/openai_ns_reconstruction'/filename;tree=ast.parse(path.read_text(encoding='utf-8'))
        node=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='_leading')
        scope={};exec(compile(ast.Module(body=[node],type_ignores=[]),str(path),'exec'),scope);methods[filename]=scope['_leading']
    buffer=SimpleNamespace(ctx=c,iv=iv,J=J,patch=SimpleNamespace(leading=SimpleNamespace(
        reference=SimpleNamespace(previous=SimpleNamespace(_restore=lambda row:iv.mpf(row['approximation']))),
        defect=SimpleNamespace(logPp=c.log(17),logP=iv.log(17)))))
    buffer._leading=lambda *args:methods['original_outer_buffer_numeric.py'](buffer,*args)
    inlet=dict(original_slope_scalar_inlet={n:dict(approximation=v) for n,v in zip(('F','H','G','C'),('.8','.7','1.1','3.2'))},
        axial_profiles={n:dict(approximation=v) for n,v in zip(('M','W','physical_log_distance'),('.125','.0625','.7'))})
    o3=SimpleNamespace(ctx=c,iv=iv,buffer=buffer,mu=c.mpf(1)/101,mui=iv.mpf(1)/101,
        _jet=lambda row,directed:J.constant(iv if directed else c,(iv if directed else c).mpf(row['approximation']),2))
    o3inlet=dict(original_slope_scalar_inlet=inlet['original_slope_scalar_inlet'],
        original_axial_profile_inlet=inlet['axial_profiles'],generic_leading_fields=dict(E=dict(approximation='.3')))
    scalar={label:{n:ctx.mpf(v) for n,v in zip(('U1','C1','C0','Cm'),('.3','.2','.1','.4'))}
        for label,ctx in (('point',c),('directed',iv))}
    source_rows=0;nonzero_mu=0
    for Z in ('-.3','0','.23'):
        for chart,values,origin in (('O2_buffer',('0','5.5','11'),0),('O3_transition',('0','.5','1'),11),('O3_power',('0','1','2'),12)):
            for x in values:
                if chart=='O2_buffer':fields,_,Q=buffer._leading(Z,x,inlet,c)
                else:
                    fields,sectors,Q=methods['original_outer_o3_numeric.py'](o3,Z,x,chart,o3inlet,scalar,c)
                    assert sectors['mu_correction_E'][0]!=0;nonzero_mu+=1
                tau=origin+c.mpf(x);z=J.variable(c,c.mpf(Z),2);mean=c.mpf('.125')*c.exp(-tau)
                for name,want in (('V',J.constant(c,0,2)),('m',4*z*mean)):
                    for k in range(3):assert abs(fields[name][k]*17-want[k])<c.mpf('1e-90');source_rows+=1
                for k in range(2):assert abs(Q[k]-(-4*mean if k==0 else 0))<c.mpf('1e-90');source_rows+=1
    # Derive the actual half-S1 density independently with finite delta and
    # nonzero mean. This catches an extra 1/2 or wrong O3 radius anchor.
    z,y,mstar,de=sy.symbols('z y mstar de');A=1-de*z*z;mean=mstar*sy.exp(-y);Q=-4*mean
    F=Q/(2*A)+Q*Q/4;G=-Q/(2*A)-Q*Q/4
    density=sy.exp(2*y)*(2*F+sy.diff(F,y)+G)-sy.exp(y)*(sy.diff(Q,y)+sy.diff(Q,y,2))
    assert sy.simplify(density+4*mstar*mstar)==0
    physical=0
    for Z in ('-.3','0','.23'):
        for origin,length in ((0,'5.5'),(11,'1'),(12,'2')):
            L=c.mpf(length)
            for k in (0,1):
                expr=sy.diff(density,z,k).subs({z:sy.Rational(Z),mstar:sy.Rational(1,8),de:sy.Rational(1,17)})
                fun=sy.lambdify(y,expr,modules=[{'exp':c.exp,'mpf':c.mpf},'mpmath'])
                want=c.quad(fun,[origin,origin+L]);got=-4*c.mpf('.125')**2*L if k==0 else c.mpf(0)
                assert abs(want-got)<c.mpf('1e-90');physical+=1
    return dict(passed=True,original_buffer_O3_leading_rows=source_rows,
        genuinely_nonzero_mu_source_cases=nonzero_mu,independent_physical_C0_C1_integrals=physical,
        exact_half_S1_density_identity=True,global_radius_mean_anchor_verified=True)


def run():
    start=time.perf_counter();fixture=independent_fixtures();print('Independent original buffer/O3 source and pressure fixtures pass',flush=True)
    buffer=OriginalPostAxialM15(precision=85);transition=buffer.next_field();power=transition.next_field()
    engines=(buffer,transition,power);a=mp.mp.clone();a.dps=190;samples=[];pressure=[]
    counts=dict(exact_inlet_replays=0,serialized_increment_rows=0,exact_empty_rows=0,pressure_rows=0,
        C3_history_rows=0,exact_history_seams=0,prefix_error_conditions=0,lower_endpoint_conditions=0,
        native_width_conditions=0,original_archive_overlaps=0,nonzero_RQy_identity_conditions=0,
        leading_Z1_zero_conditions=0,actual_fast_unit_conditions=0)
    bounds=lambda row:tuple(a.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples'])
    highlog=lambda row:-a.inf if row['exact_zero'] else a.make_mpf(tuple(row['log_absolute_upper']['upper_exact_mpf_tuple']))
    archives={}
    for suffix in ('outer_buffer','outer_o3'):
        path=DIRECTORY/(PREFIX+'compliant_current_original_'+suffix+'_numeric.json');raw=path.read_bytes();old=json.loads(raw)
        assert old['checks']['passed'] and old['source_family']==buffer.family
        for engine in engines:
            engine.hashes[path.name]=hashlib.sha256(raw).hexdigest()
            for name,digest in old['implementation_hashes'].items():
                assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
                assert engine.code_hashes.get(name,digest)==digest;engine.code_hashes[name]=digest
        for p in old['samples']:
            chart='O2_buffer' if suffix=='outer_buffer' else p['native_O3_chart']
            coord=p['native_buffer' if suffix=='outer_buffer' else 'native_O3'];archives[chart,p['Z'],coord]=p
    for Z in ('-.3','0','.3'):
        for engine,values in ((buffer,('0','5.5','11')),(transition,('0','1')),(power,('0','2'))):
            for value in values:
                out=engine.evaluate(Z,value);p=engine.rows.evaluate(Z,value);samples.append(out);pressure.append(p)
                count=engine.COUNTS[engine.chart];assert len(out['m15_C0_C1_terms'])==count+4
                assert out['m15_C0_C1_terms'][:count]==out['incoming_m15_terms']==engine.inlet(Z)['m15_C0_C1_terms'];counts['exact_inlet_replays']+=1
                assert all(out['contract'].values()) and all(p['contract'].values())
                assert engine.evaluate(Z,value) is out
                for term in out['m15_C0_C1_terms'][count:]:
                    for row in term['coefficient_C0_C1']:
                        lo,hi=bounds(row);point=a.mpf(row['approximation']);assert max(abs(lo-point),abs(hi-point))<=a.mpf(row['absolute_error_upper'])
                        counts['serialized_increment_rows']+=1
                        if value=='0':assert lo==hi==point==0;counts['exact_empty_rows']+=1
                leading=out['m15_C0_C1_terms'][count]
                assert bounds(leading['coefficient_C0_C1'][1])==(0,0);counts['leading_Z1_zero_conditions']+=1
                pp,tt,geo=engine.rows.background(Z,value)
                for name,field in p['physical_pressure_rows'].items():
                    for k,row in enumerate(field['signed_background_ordinary_Z']):
                        assert row['source_enclosure_exact_mpf_tuples']==[list(v) for v in (tt[name][k]*math.factorial(k))._mpi_];counts['pressure_rows']+=1
                mean=engine.rows.scalars(value)['truth']['mean'];assert _bounds(a,mean)[0]>0
                qy=tt['Q_y'][0];lo,hi=_bounds(a,qy+tt['Q'][0]);assert lo<=0<=hi and _bounds(a,qy)[0]>0
                relative=qy/(4*engine.rows.mean_truth)*engine.iv.exp(engine.rows.ORIGINS[engine.chart]+engine.iv.mpf(value))
                lo,hi=_bounds(a,relative);assert lo<=1<=hi;counts['nonzero_RQy_identity_conditions']+=1
                for name,values4 in p['finite_N_history_error_logs'].items():
                    assert len(values4)==4 and all(not row['exact_zero'] for row in values4);counts['C3_history_rows']+=4
                if value=='0':assert p['finite_N_history_error_logs']==p['finite_N_incoming_history_error_logs'];counts['exact_history_seams']+=1
                else:
                    for kind in ('endpoint','viscous'):
                        budget=out['actual_lower_finite_N_endpoint_budget'][kind];assert all(not row['exact_zero'] for row in budget);counts['lower_endpoint_conditions']+=2
                prefix,_=engine.rows.error_data(Z,value,prefix=True);endpoint,_=engine.rows.error_data(Z,value)
                for name in ('U','m','Q','U_y','Q_y'):
                    for v,w in zip(engine.rows.log_rows(prefix[name]),engine.rows.log_rows(endpoint[name])):
                        assert highlog(v)>=highlog(w);counts['prefix_error_conditions']+=1
                h=engine.rows.history_logs(Z,value);width=sum((step['physical_log_width'] for step in h['steps']),engine.iv.mpf(0))
                lo,hi=_bounds(a,width);assert lo<=a.mpf(value)<=hi;counts['native_width_conditions']+=1
                old=archives.get((engine.chart,Z,value))
                if old:
                    for name,oldname in (('U','Uz'),('m','Mz_over_R'),('Q','Ur_over_sqrt_R_over2')):
                        for k,row in enumerate(old['physical_similarity'][oldname]):
                            lo,hi=bounds(row);ll,hh=_bounds(a,tt[name][k]);ll*=math.factorial(k);hh*=math.factorial(k)
                            assert max(lo,ll)<=min(hi,hh),(Z,engine.chart,value,name,k);counts['original_archive_overlaps']+=1
                for k,key in enumerate(('B_phi','B_Z_phi')):
                    source=engine.rows.read_log(engine.rows.loop['phase_and_first_slow_phase_log_bounds'][key]);row=p['finite_N_axial_radial_error_components']['fast_B_phi'][k]
                    if source is None:assert row['exact_zero']
                    else:assert highlog(row)==_bounds(a,source+engine.rows.o.core.lp_iv)[1]
                    counts['actual_fast_unit_conditions']+=1
            print('Original '+engine.chart+' pressure/m15: Z='+Z+', '+str(count)+' ->'+str(count+4)+' terms',flush=True)
    rejected=0
    for engine,Z,value in ((buffer,'1','0'),(buffer,'0','11.1'),(transition,'0','1.1'),(power,'0','2.1'),(power,'0','-.1')):
        try:engine.evaluate(Z,value)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid post-axial query admitted')
    hashes={};code={}
    for engine in engines:hashes.update(engine.hashes);code.update(engine.code_hashes)
    name='scripts/run_original_post_axial_m15.py';code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=buffer.family,samples=samples,pressure_row_samples=pressure,input_hashes=hashes,implementation_hashes=code,
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
