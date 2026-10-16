"""Actual compact n1 primitives and physical orders0/1, pressure kept separate."""
import argparse
import ast
import hashlib
import inspect as inspectlib
import json
from pathlib import Path
import sys
import textwrap
import time
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_reserved_n1_numeric import OriginalReservedN1Numeric
from openai_ns_reconstruction.original_compact_n1_velocity import OriginalCompactN1Velocity
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_compact_n1_numeric.json'


def independent_primitives():
    c=mp.mp.clone();c.dps=70;iv=MPIntervalContext();iv.dps=90
    f=object.__new__(OriginalReservedN1Numeric);f.ctx,f.iv,f.J=c,iv,J;f.ell=iv.log(2)/400
    f.cells=256;f.beta_jets=beta_jets;f._integrals={};count=0;maximum=c.mpf(0)
    for sign in (-1,1):
        for r in ('-1','-.75','0','.75','1'):
            got=f._integral_to(c.mpf(r),sign)
            truth=c.quad(lambda t:c.exp(-1/(1-t*t)+sign*c.log(2)*t/400) if abs(t)<1 else c.mpf(0),[-1,c.mpf(r)])
            lo,hi=_bounds(c,got);assert lo-c.mpf('1e-60')<=truth<=hi+c.mpf('1e-60'),(r,sign);count+=1
            maximum=max(maximum,hi-lo)
    return count,c.nstr(maximum,16)


def independent_beta_unit():
    import sympy as sy
    import lei_ren_part1_paper_compliant_current_original_n1_five_defect_values as original
    Rc,lp,lu=sy.symbols('logRc logP logUp',real=True);lm=-4*lp-sy.log(1000)
    ctx=SimpleNamespace(mpf=lambda v:sy.Rational(str(v)),log=sy.log,ln=sy.log,exp=sy.exp,sqrt=sy.sqrt,pi=sy.pi)
    namespace=dict(original.pressure_unit_log.__globals__);namespace['interval']=lambda c,v:v
    tree=ast.parse(textwrap.dedent(inspectlib.getsource(original.pressure_unit_log)));exec(compile(tree,'<independent original pressure unit>','exec'),namespace)
    target=namespace['pressure_unit_log'](ctx,dict(logP=lp,logUp=lu,logRc=Rc),2)-lp-(-lp/2+sy.Rational(1,2))-Rc
    f=object.__new__(OriginalReservedN1Numeric);f.ctx=f.iv=ctx;f.logUp=lu;f._mid=lambda v:v;f._source_row=lambda v:v
    f.background=SimpleNamespace(patch=SimpleNamespace(leading=SimpleNamespace(reference=SimpleNamespace(_summary=lambda p,t:t),defect=SimpleNamespace(logP=lp))),outer=SimpleNamespace(logmui=lm))
    actual=sum(f.beta_unit_atoms(dict(coordinate=dict(global_y_log_terms=[Rc-sy.log(110),0,0,0,sy.Symbol('local_s')]))))
    assert sy.simplify(sy.expand(actual-target))==0
    return True


def run():
    start=time.perf_counter();independent,width=independent_primitives();unit=independent_beta_unit()
    print('Independent weighted-beta primitives and original physical beta-unit identity pass',flush=True)
    f=OriginalCompactN1Velocity();c,iv=f.ctx,f.iv;samples=[];source=mass=active=physical=seams=terminal=0
    # Exact decimal physical-radius ratios approximate source bump centers;
    # all computations use their directed log(x), not a rounded chosen phase.
    radial=[('inlet','1')]
    for j in (12,13,14,15,16):
        radial.append(('center'+str(j),c.nstr(c.exp(c.log(2)*j/40),95)))
    radial.extend([('partial12',c.nstr(c.exp(c.log(2)*(c.mpf(12)/40+c.mpf(1)/800)),95)),
        ('partial16',c.nstr(c.exp(c.log(2)*(c.mpf(16)/40-c.mpf(1)/800)),95)),('suffix','1.9'),('terminal','2')])
    for Z in ('-.3','0','.3'):
        for label,x in radial:
            p=f.evaluate(Z,x,log_tau='-2',azimuth='.25');p['sample_label']=label;samples.append(p);n=p['reserved_n1']
            assert p['physical']['coefficient_orders_included']==[0,1];physical+=1
            assert not p['physical']['n1_extra_Lambda_squared_factor_applied'];physical+=1
            assert p['physical']['n1_swirl_is_Utheta1_not_F1'];physical+=1
            assert not p['physical']['full_n1_pressure_included'];physical+=1
            assert n['contract']['same_actual_source_terminal_mass_identity'];source+=1
            assert not n['limits']['Q1_Z_derivative'] and not n['limits']['resolved_signed_coefficient_functions'];source+=1
            assert p['selected_integer']==f.numeric._inlets[Z]['selected_integer'];source+=1
            for k in (0,1):
                row=f._restore(n['M1_C0_C1'][k]);fw=f._restore(n['forward_M1_C0_C1'][k]);bw=f._restore(n['backward_terminal_M1_C0_C1'][k])
                lo,hi=_bounds(c,row);fl,fh=_bounds(c,fw);bl,bh=_bounds(c,bw)
                assert max(fl,bl)<=lo<=hi<=min(fh,bh);mass+=1
            z=iv.mpf(Z);delta=iv.mpf(f.amplitude.delta_iv)
            RUz=f._restore(n['R_times_Uz1_C0_C1'][0]);M=f._restore(n['M1_C0_C1'][0]);Mz=f._restore(n['M1_C0_C1'][1])
            v=f._restore(n['R_times_Q1'][0])-(2*z*RUz-(1+delta)*z*M-(1-z*z)*Mz)/(1-delta*z*z)
            lo,hi=_bounds(c,v);assert lo<=0<=hi;mass+=1
            if label=='inlet':
                incoming=p['incoming_n1_history_not_extended']['actual_terminal_primitive']['actual_primitive']
                for k in (0,1):
                    v=f._restore(n['M1_C0_C1'][k])-f._restore(incoming[k])*iv.exp(-2*f.amplitude.ll_iv)
                    lo,hi=_bounds(c,v);assert lo<=0<=hi;seams+=1
                assert n['active_bumps']==[];seams+=1
            elif label in ('suffix','terminal'):
                assert n['active_bumps']==[];terminal+=1
                for name in ('M1_C0_C1','R_times_Q1','Uz1_times_Rc_C0_C1','Utheta1_normalized_C0_C1'):
                    for row in n[name]:assert _bounds(c,f._restore(row))==(0,0);terminal+=1
                assert _bounds(c,f._restore(p['generic_leading_fields']['E'][0]))[0]>0;terminal+=1
            else:
                index=int(label.removeprefix('center').removeprefix('partial'))-12
                assert n['active_bumps']==[index];active+=1
                name='Uz1_times_Rc_C0_C1' if index in (0,4) else 'Utheta1_normalized_C0_C1'
                lo,hi=_bounds(c,f._restore(n[name][0]));assert lo<hi and (lo!=0 or hi!=0);active+=1
        print('Actual compact n1 velocity and own mass: Z='+Z,flush=True)
    scaled=0
    for Z in ('-.3','0','.3'):
        label,x=next(v for v in radial if v[0]=='center12')
        p=f.evaluate(Z,x,log_tau=None,scaled_time='1e18',azimuth='1.1');p['sample_label']='scaled_center12';samples.append(p)
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
        normal=next(a for a in samples if a['Z']==Z and a['sample_label']=='center12')
        for group in ('ur','utheta','uz'):
            for a,b in zip(p['physical']['cylindrical'][group],normal['physical']['cylindrical'][group]):overlap(a['coefficient'],b['coefficient']);scaled+=1
    rejected=0
    for Z,x in (('1','1.2'),('nan','1.2'),('.3','0.9'),('.3','2.1')):
        try:f.reserved.evaluate(Z,x)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid reserved n1 request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_weighted_beta_defining_integrals=independent,max_tested_beta_integral_width=width,
            independent_original_beta_physical_unit_identity=unit,actual_compact_leading_n1_physical_packets=len(samples),
            source_N_and_scope_conditions=source,mass_recovery_and_intersection_conditions=mass,actual_support_conditions=active,
            physical_units_orders_and_pressure_scope_conditions=physical,Rc_normalized_to_physical_mass_seams=seams,
            source_identity_zero_suffix_and_nonzero_background_conditions=terminal,scaled_time_conditions=scaled,rejected_requests=rejected),
        execution_seconds=round(time.perf_counter()-start,3),limits=dict(compact_physical_temporal_orders=[0,1],
            actual_reserved_velocity_values=True,mass_C0_C1=True,Q1_value_only=True,Q1_Z=False,n1_pressure=False,
            coefficient_enclosures_not_resolved_signed_functions=True,terminal_mass_zero_uses_checked_source_identity=True,
            independent_resolved_terminal_residual=False,global_n1_joins=False,common_global_N=False,physical_gradients=False,
            stress_cone=False,exact_heat_exterior=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for n,h in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/n).read_bytes()).hexdigest()==h,n
    for n,h in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
    return d


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write-report',action='store_true');p.add_argument('--inspect-report',action='store_true')
    args=p.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
