"""Evaluate actual first-micro signed m15 with complete factored errors."""
import argparse
import hashlib
import importlib
import itertools
import json
from pathlib import Path
import sys
import time
from types import FunctionType,SimpleNamespace

import mpmath as mp
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073'),str(ROOT/'scripts')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_first_micro_m15 import OriginalFirstMicroM15
from run_original_n1_inner_m15_numeric import independent_original_source

NAME=PREFIX+'compliant_current_original_first_micro_m15_numeric.json'


class SymbolicAlgebra:
    def __init__(self,g):pass
    def fixed(self,v):return sy.sympify(v)
    def add(self,*v):return sum(v)
    def mul(self,*v):return sy.prod(v)
    def neg(self,v):return -v
    def sub(self,a,b):return a-b
    def div(self,a,b):return a/b
    def scale(self,a,b):return sy.sympify(a)*b


def independent_S1():
    """Execute the accepted S1 body and compare against physical density."""
    module=importlib.import_module(PREFIX+'compliant_current_original_n1_preRp_moment_functions')
    scope=dict(module.__dict__,C1Algebra=SymbolicAlgebra,
               source=SimpleNamespace(C1Function=lambda value,Z:value))
    fn=FunctionType(module.pressure_functions.__code__,scope)
    r,z=sy.symbols('R Z',positive=True);de=sy.Rational(1,17)
    U=1+z*z+r*(sy.Rational(1,5)+z/10)+r*r*z/30
    m=1+z*z+r*(sy.Rational(1,5)+z/10)/2+r*r*z/90+(sy.Rational(1,10)+z*z/20)/r
    assert sy.cancel(r*sy.diff(m,r)-(U-m))==0
    d,L=1-z*z,1-de*z*z
    Q=(2*z*U-(1-de)*z*m-d*sy.diff(m,z))/L;Qy=r*sy.diff(Q,r)
    flow=(Q+(1-de)*z*sy.diff(Q,z)/2+Qy)/L+Q*(Q/2+Qy)
    flow+=U*(-2*z*Q+d*sy.diff(Q,z)-2*z*Qy)/L
    integrands=dict(I1_q=Q,I1_q_Z=sy.diff(Q,z),I1_q_squared=Q*Q,
        I1_U_q_Z=U*sy.diff(Q,z),I1_U_squared=U*U,I1_U_U_Z=U*sy.diff(U,z))
    for key in ('I0_q_Z','I0_U_q_Z','I0_U_squared','I0_U_U_Z'):integrands[key]=sy.Integer(0)
    endpoint=dict(R_q=r*Q,R_squared_q=r*r*Q,R_q_squared=r*Q*Q,
        R_squared_q_squared=r*r*Q*Q,R_U_squared=r*U*U,R_squared_U_squared=r*r*U*U,
        q_plus_q_y=Q+Qy,R_q_y=r*Qy)
    g=SimpleNamespace(symbol=lambda name:z,one=sy.Integer(1))
    # Differentiate each actual endpoint body; integrate r times each I1
    # density. This is the fundamental theorem applied to the actual source.
    zero={k:sy.Integer(0) for k in endpoint}
    source=fn(g,de,{k:r*v if k.startswith('I1') else v for k,v in integrands.items()},
              zero,{k:sy.diff(v,r) for k,v in endpoint.items()})['m15']
    density=r*flow/2-r*(2*sy.diff(Q,r)+r*sy.diff(Q,r,2))
    assert sy.cancel(source-density)==0
    assert sy.cancel(sy.diff(source,z)-sy.diff(density,z))==0
    count=0
    with mp.workdps(85):
        c=mp.mp;nodes,weights=c.gauss_quadrature(64,'legendre')
        ep=c.mpf('.3');h=c.mpf('.08');a=c.mpf('.02');b=c.mpf('.24')
        integrators={k:sy.lambdify((r,z),v,'mpmath') for k,v in integrands.items() if k.startswith('I1')}
        integratorsZ={k:sy.lambdify((r,z),sy.diff(v,z),'mpmath') for k,v in integrands.items() if k.startswith('I1')}
        endpoint_functions={k:sy.lambdify((r,z),v,'mpmath') for k,v in endpoint.items()}
        endpointZ={k:sy.lambdify((r,z),sy.diff(v,z),'mpmath') for k,v in endpoint.items()}
        density_fn=[sy.lambdify((r,z),v,'mpmath') for v in (density,sy.diff(density,z))]
        R=lambda s:ep*4*c.exp(h*s)
        for Z in ('-.3','0','.3'):
            zz=c.mpf(Z)
            for k,funcs in enumerate((integrators,integratorsZ)):
                values={key:c.mpf(0) for key in integrands}
                direct=c.mpf(0)
                for node,w in zip(nodes,weights):
                    s=a+(node+1)*(b-a)/2;rad=R(s);ds=w*(b-a)/2
                    for key,f in funcs.items():values[key]+=f(rad,zz)*rad*rad*h*ds
                    direct+=density_fn[k](rad,zz)*rad*h*ds
                ends=endpoint_functions if k==0 else endpointZ
                lo={key:f(R(a),zz) for key,f in ends.items()};hi={key:f(R(b),zz) for key,f in ends.items()}
                if k==0:
                    actual=fn(SimpleNamespace(symbol=lambda name:zz,one=c.mpf(1)),c.mpf(1)/17,values,lo,hi)['m15']
                else:
                    # The derivative also differentiates the Z coefficients.
                    # Build the entire source primitive at symbolic Z first.
                    I={key:sy.Symbol(key) for key in values};low={key:sy.Symbol('lo_'+key) for key in lo};high={key:sy.Symbol('hi_'+key) for key in hi}
                    expr=fn(g,de,I,low,high)['m15'];derivative=sy.diff(expr,z)
                    term=sy.Integer(0)
                    for key,v in I.items():term+=sy.diff(expr,v)*values[key]
                    for key,v in low.items():term+=sy.diff(expr,v)*lo[key]
                    for key,v in high.items():term+=sy.diff(expr,v)*hi[key]
                    # Coefficient derivatives use independent value integrals.
                    values0={key:c.mpf(0) for key in values}
                    for node,w in zip(nodes,weights):
                        s=a+(node+1)*(b-a)/2;rad=R(s)
                        for key,f in integrators.items():values0[key]+=f(rad,zz)*rad*rad*h*w*(b-a)/2
                    subs={z:zz,**{I[key]:v for key,v in values0.items()},
                          **{low[key]:f(R(a),zz) for key,f in endpoint_functions.items()},
                          **{high[key]:f(R(b),zz) for key,f in endpoint_functions.items()}}
                    actual=c.mpf(str(sy.N((term+derivative).subs(subs),80)))
                assert abs(c.mpf(actual)-direct)<c.mpf('1e-60');count+=1
    return dict(passed=True,actual_accepted_S1_source_density_identities=2,
                independent_physical_native_measure_value_Z1_integrals=count)


def product_errors(f):
    c,iv=f.ctx,f.iv
    a,b=[c.mpf('1.3'),c.mpf('-.2')],[c.mpf('.7'),c.mpf('.4')]
    ea,eb=[c.mpf('.03'),c.mpf('.04')],[c.mpf('.02'),c.mpf('.01')]
    logs=lambda values:[iv.log(abs(iv.mpf(v))) for v in values]
    bound=f.delta_product(logs(a),logs(b),logs(ea),logs(eb));n=0
    for signs in itertools.product((-1,1),repeat=4):
        x,y=[ea[k]*signs[k] for k in (0,1)],[eb[k]*signs[k+2] for k in (0,1)]
        actual=[(a[0]+x[0])*(b[0]+y[0])-a[0]*b[0],
            (a[1]+x[1])*(b[0]+y[0])+(a[0]+x[0])*(b[1]+y[1])-a[1]*b[0]-a[0]*b[1]]
        for k in (0,1):assert abs(actual[k])<=c.exp(_bounds(c,bound[k])[1]);n+=1
    return n


def run():
    start=time.perf_counter();original=independent_original_source();S1=independent_S1()
    print('Original Laurent pressure and actual reduced-S1 physical C0/C1 identities pass',flush=True)
    f=OriginalFirstMicroM15(cells=8,precision=80);products=product_errors(f)
    coarse=OriginalFirstMicroM15(cells=4,rows=f.rows);audit=mp.mp.clone();audit.dps=180
    samples=[];inlets=serialized=empty=atoms=finite=refinements=0
    for Z in ('-.3','0','.3'):
        for coordinate,selected in (('.5',True),('1',True),('.125',False),('.25',False)):
            p=f.evaluate(Z,coordinate,selected_sc_multiple=selected);samples.append(p)
            assert all(p['contract'].values());assert len(p['m15_C0_C1_terms'])==13
            assert p['m15_C0_C1_terms'][:6]==p['incoming_q2_m15_terms'];inlets+=1
            for term in p['m15_C0_C1_terms'][6:]:
                for row in term['coefficient_C0_C1']:
                    ends=[audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']]
                    a,e=audit.mpf(row['approximation']),audit.mpf(row['absolute_error_upper'])
                    assert max(abs(v-a) for v in ends)<=e;serialized+=1
                    if p['exact_empty_increment']:assert ends[0]==ends[1]==a==0;empty+=1
                if 'finite_N_' in term['name']:
                    assert term['actual_function_binding']['error_enclosure_only'];finite+=1
            if not p['exact_empty_increment']:
                lower=p['m15_C0_C1_terms'][8]
                assert lower['name']=='first_micro_singular_lower'
                assert lower['actual_function_binding']['endpoint_sign']==-1
                assert lower['log_scale_terms'][1]['original_nonzero_sigma_log_atom'];atoms+=1
            assert not p['limits']['resolved_finite_N_signed_centers']
        q=coarse.evaluate(Z,'.25');samples.append(q)
        fine=f.evaluate(Z,'.25')
        for k in (0,1):
            a=fine['m15_C0_C1_terms'][6]['coefficient_C0_C1'][k];b=q['m15_C0_C1_terms'][6]['coefficient_C0_C1'][k]
            assert abs(audit.mpf(a['approximation'])-audit.mpf(b['approximation']))<=audit.mpf(a['absolute_error_upper'])+audit.mpf(b['absolute_error_upper']);refinements+=1
        print('Original cumulative m15 integrated through first_micro 1/4: Z='+Z,flush=True)
    rejected=0
    for args in [('1','.25'),('.3','.3'),('.3','0'),('nan','.25')]:
        try:f.evaluate(*args)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original first-micro m15 query admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_actual_Laurent_pressure_source=original,independent_actual_S1=S1,
            actual_original_packets=len(samples),retained_exact_inlet_partitions=inlets,
            serialized_increment_errors=serialized,exact_empty_increment_rows=empty,
            retained_nonzero_lower_sigma_atoms=atoms,finite_N_error_only_sectors=finite,
            independent_quadratic_product_error_sign_cases=products,directed_cell_refinement_C0_C1=refinements,
            rejected_invalid_queries=rejected),execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(first_micro_m15_leading_signed_with_complete_factored_finite_N_errors=True,
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
