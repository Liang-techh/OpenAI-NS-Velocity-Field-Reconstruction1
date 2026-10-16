"""Build original second-micro pressure continuation with both endpoints."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from types import SimpleNamespace

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073'),str(ROOT/'scripts')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_second_micro_m15 import OriginalSecondMicroM15
from openai_ns_reconstruction.original_n1_inner_m15_numeric import leading_regular_frame,normalized_pressure_source
from openai_ns_reconstruction.original_first_collar_axial2 import OriginalFirstCollarAxial2
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J
from run_original_first_micro_m15 import independent_S1

NAME=PREFIX+'compliant_current_original_second_micro_m15_numeric.json'


def constant_control_source():
    ctx=SimpleNamespace(mpf=sy.sympify)
    jet=lambda a:J(ctx,[sy.Rational(str(v)) for v in a])
    const=lambda v:J.constant(ctx,sy.sympify(v),5)
    z=J.variable(ctx,sy.Rational(1,5),5);rho=const(4);x=const(2)
    de,e,P2=sy.Rational(1,17),sy.Rational(1,3),sy.Integer(5);h=sy.Symbol('hb',positive=True)
    comparison=dict(Phi=[jet([2,.1,.03,.01,.001,.0001]),jet([-.2,.04,.01,.001,.0001,0]),jet([.08,-.01,.001,0,0,0])],
        Psi=[jet([.3,.1,.01,.001,0,0]),jet([-.4,.02,.001,0,0,0]),jet([.07,.01,0,0,0,0])])
    values=dict(Phi=jet([2.1,.11,.03,.01,.001,.0001]),Psi=comparison['Psi'][0],Mpsi=jet([.2,.03,.002,0,0,0]))
    f=leading_regular_frame(values,comparison,z,rho,x,de,e,P2,const(1),h,sy.Rational(1,19))
    f['g']=jet([.2,.02,.001,0,0,0])
    owner=SimpleNamespace(leading=SimpleNamespace(_derivative=lambda q:J(q.ctx,[(k+1)*q[k+1] for k in range(q.order)])))
    actual=OriginalFirstCollarAxial2._forcing(owner,f)['p']/(2*x)
    normalized=normalized_pressure_source(f);assert set(actual.rows)=={0};assert set(normalized.rows)=={0}
    A=comparison['Phi'][1]/comparison['Phi'][0]
    Ar=comparison['Phi'][2]/comparison['Phi'][0]-A*A;ratio=values['Phi']/comparison['Phi'][0]
    tests=0
    for k in (0,1):
        assert sy.cancel(actual.rows[0][k]-normalized.rows[0][k])==0;tests+=1
        assert sy.cancel(f['Phi'][2].rows[0][k]-(values['Phi']*(h*h*A*A+h*Ar))[k])==0;tests+=1
        assert sy.cancel(f['U'][1].rows[0][k]/e-(ratio*comparison['Psi'][1]*h)[k])==0;tests+=1
        expected=ratio*(A*comparison['Psi'][1]*(h*(h-1))+comparison['Psi'][2]*h)
        assert sy.cancel(f['U'][2].rows[0][k]/e-expected[k])==0;tests+=1
    return dict(passed=True,actual_constant_hb_radial_and_pressure_C0_C1_identities=tests,
        no_sigma_prime_or_inverse_hb_Laurent_row=True,effective_sigma_one_means_chi_hb=True)


def run():
    start=time.perf_counter();source=constant_control_source();S1=independent_S1()
    print('Actual constant-hb second pressure source and reduced-S1 physical identities pass',flush=True)
    f=OriginalSecondMicroM15(precision=80);audit=mp.mp.clone();audit.dps=180
    samples=[];inlets=zeros=errors=frames=lower=controls=0
    def check(Z,coordinate):
        nonlocal inlets,zeros,errors,frames,lower,controls
        p=f.evaluate(Z,coordinate);samples.append(p)
        assert all(p['contract'].values()) and len(p['m15_C0_C1_terms'])==17
        assert p['m15_C0_C1_terms'][:13]==p['incoming_first_micro_m15_terms'];inlets+=1
        assert p['actual_selected_integer']==f.inlet(Z)['actual_selected_integer']
        assert not any('singular' in t['name'] for t in p['m15_C0_C1_terms'][13:]);controls+=1
        metadata=p['second_micro_full_prefix_cover']
        assert metadata['generic_derivative_chart']=='bridge_second'
        assert metadata['cumulative_history_charts']==['first_micro','second_micro'];frames+=1
        if coordinate=='2':assert metadata['uniform_prefix_source_cells']==8*len(metadata['exact_Z_cells']);frames+=1
        for term in p['m15_C0_C1_terms'][13:]:
            for row in term['coefficient_C0_C1']:
                ends=[audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']]
                assert max(abs(v-audit.mpf(row['approximation'])) for v in ends)<=audit.mpf(row['absolute_error_upper']);errors+=1
                if p['exact_empty_increment']:assert ends[0]==ends[1]==audit.mpf(row['approximation'])==0;zeros+=1
        if not p['exact_empty_increment']:
            for name in ('endpoint','viscous'):
                for k,row in enumerate(p['incoming_finite_N_endpoint_error_logs'][name]):
                    if not row['exact_zero']:
                        total=p['finite_N_error_logs_before_physical_measures'][name][k]
                        assert not total['exact_zero']
                        x=audit.make_mpf(tuple(row['log_absolute_upper']['upper_exact_mpf_tuple']))
                        y=audit.make_mpf(tuple(total['log_absolute_upper']['upper_exact_mpf_tuple']))
                        assert y>=x;lower+=1
        assert not p['limits']['resolved_finite_N_signed_centers']
        return p
    for Z in ('-.3','0','.3'):
        for coordinate in ('1','1.25','1.5','1.75','2'):check(Z,coordinate)
        print('Original cumulative m15 reaches second_micro endpoint s=2: Z='+Z,flush=True)
    check('.3','1.875')
    assert lower>0
    # Whole-prefix leading field errors explicitly include the completed
    # first window even at the empty second-moment increment.
    for Z in ('-.3','0','.3'):
        v=f.rows.inputs(Z,f.iv.mpf(1))
        assert v['prefix']=='full first control window plus second suffix s-1'
        assert _bounds(audit,v['delta']['Phi'])[1]>0;controls+=1
    rejected=0
    for args in [('1','2'),('.3','.9'),('.3','2.1'),('nan','2')]:
        try:f.evaluate(*args)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original second-micro pressure query admitted')
    code=dict(f.code_hashes)
    for name in ('scripts/run_original_first_micro_m15.py','scripts/run_original_second_micro_m15.py'):
        code[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,actual_constant_control_source=source,independent_actual_S1=S1,
            actual_original_second_micro_packets=len(samples),exact_complete_first_inlet_replays=inlets,
            serialized_increment_error_conditions=errors,exact_empty_second_increment_rows=zeros,
            cumulative_source_frame_conditions=frames,nonzero_lower_endpoint_budget_conditions=lower,
            constant_control_and_first_window_memory_conditions=controls,rejected_invalid_queries=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(second_micro_leading_signed_m15_with_complete_factored_finite_N_errors=True,
            resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
            remaining_active_compact_signed_pressure=False,independent_preRp_signed_seam=False,
            global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
            stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


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
