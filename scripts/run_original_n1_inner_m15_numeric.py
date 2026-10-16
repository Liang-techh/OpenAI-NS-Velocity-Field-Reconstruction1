"""Original signed m15 through the selected cutoff and lower collar."""
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

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_n1_inner_m15_numeric import (
    OriginalN1InnerM15Numeric,leading_regular_frame,normalized_pressure_source,Laurent)
from openai_ns_reconstruction.original_first_collar_axial2 import OriginalFirstCollarAxial2
from lei_ren_part1_paper_interval_taylor import IntervalTaylor as J

NAME=PREFIX+'compliant_current_original_n1_inner_m15_numeric.json'


def independent_original_source():
    import sympy as sy
    ctx=SimpleNamespace(mpf=sy.sympify)
    jet=lambda a:J(ctx,[sy.Rational(str(v)) for v in a])
    const=lambda a:J.constant(ctx,sy.Rational(str(a)),5)
    z=J.variable(ctx,sy.Rational(1,5),5);r=const(4);x=const(2)
    de,e,P2=map(sy.Rational,('1/17','1/3','5'))
    hb,sp=sy.symbols('original_h sigma_prime',real=True)
    cmp=dict(Phi=[jet([2,.1,.03,.01,.001,.0001]),jet([-.2,.04,.01,.001,.0001,0]),jet([.08,-.01,.001,0,0,0])],
             Psi=[jet([.3,.1,.01,.001,0,0]),jet([-.4,.02,.001,0,0,0]),jet([.07,.01,0,0,0,0])])
    val=dict(Phi=jet([2.1,.11,.03,.01,.001,.0001]),Psi=cmp['Psi'][0],Mpsi=jet([.2,.03,.002,0,0,0]))
    v=leading_regular_frame(val,cmp,z,r,x,de,e,P2,const('.07'),hb,sy.Rational(1,19))
    v.update(g=jet([.2,.02,.001,0,0,0]),S=jet([.03,.01,.001,0,0,0]))
    full=dict(v);full['U']=list(v['U']);full['Q']=list(v['Q'])
    ratio=val['Phi']/cmp['Phi'][0]
    inverse=-(ratio*(1-hb)*cmp['Psi'][1]/r)*sp*e
    full['U'][2]=Laurent({0:v['U'][2].rows[0],-1:inverse})
    d,L=1-z*z,1-z*z*de
    full['Q'][2]=(full['U'][2]*(2*z)-full['M'][2]*(z*(1-de))-full['M'][2].dz()*d)/L
    owner=SimpleNamespace(leading=SimpleNamespace(_derivative=lambda q:J(q.ctx,[(k+1)*q[k+1] for k in range(q.order)])))
    original=OriginalFirstCollarAxial2._forcing(owner,full)['p']/(2*x)
    reduced=normalized_pressure_source(v)
    A=ratio*cmp['Psi'][1]*z*r*r*(2*(1-hb))/L
    checks=0
    for k in (0,1):
        assert sy.cancel(original.rows[0][k]-reduced.rows[0][k])==0;checks+=1
        assert sy.cancel(-P2*r[k*0]**2*original.rows[-1][k]/sp-A[k])==0;checks+=1
    # Execute the actual cutoff coupling expression, retaining its units.
    import openai_ns_reconstruction.original_n1_cutoff_numeric as source
    tree=ast.parse(textwrap.dedent(inspectlib.getsource(source.OriginalN1CutoffNumeric._evaluate)))
    expr=next(n.value for n in ast.walk(tree) if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id=='coupling' for t in n.targets))
    ff=jet([.2,-.03,.001,0,0,0])
    actual=eval(compile(ast.Expression(expr),'<actual cut-F pressure source>','eval'),{},
        dict(raw=dict(S=v['S'],Phi=v['Phi'],e=e,P2=P2),rt=dict(f=ff)))
    expected=(v['S']*v['Phi'][0].rows[0]*ff*(2/(e**3*P2))).truncate(2)
    for k in (0,1):assert sy.cancel(actual[k]-expected[k])==0;checks+=1
    return dict(actual_regular_Np_and_singular_Laurent_C0_C1=4,actual_cutF_units_C0_C1=2,
        total=checks,entire_one_minus_h_in_singular_coefficient=True)


def independent_flat_and_stieltjes():
    c=mp.mp.clone();c.dps=75;iv=MPIntervalContext();iv.dps=105
    f=object.__new__(OriginalN1InnerM15Numeric);f.ctx=c;f.iv=iv
    phase=c.mpf(1)/64
    f.cut=SimpleNamespace(phase_iv=iv.mpf(1)/64,phase_tuple=list(phase._mpf_))
    f.o=SimpleNamespace(_summary=lambda p,t:dict(approximation=c.nstr(p,75),bits=[list(v) for v in t._mpi_]))
    flat=0
    for q in ('.5','1','2'):
        atom,k=f.sigma_atom(q);t=phase*c.mpf(q)
        odds=1/(1-t)**2-1/t**2;sig=c.exp(odds)/(1+c.exp(odds))
        lo,hi=_bounds(c,k);assert lo<=sig/c.exp(odds)<=hi;flat+=1
        a,b=[c.make_mpf(tuple(v)) for v in atom['bits']]
        assert a-c.mpf('1e-65')<=odds<=b+c.mpf('1e-65');flat+=1
    def sigma(t):
        return c.mpf(0) if t==0 else 1/(1+c.exp(1/t**2-1/(1-t)**2))
    t=c.mpf('.12');h=c.mpf('.02');sig=sigma(t)
    nodes,weights=c.gauss_quadrature(96,'legendre');tests=0
    for z in (c.mpf('-.3'),c.mpf(0),c.mpf('.3')):
        integral=[c.mpf(0)]*2;A0=1+z*z
        for node,w in zip(nodes,weights):
            s=(node+1)*t/2;ss=sigma(s)
            ds=ss*(1-ss)*(2/s**3+2/(1-s)**3)
            integral[0]+=(A0+h*(2+z+s))*ds*w*t/2
            integral[1]+=(2*z+h)*ds*w*t/2
        for k,(base,E) in enumerate(((A0,2+abs(z)+t),(2*z,c.mpf(1)))):
            assert abs(integral[k]-base*sig)<=h*E*sig+c.mpf('1e-60');tests+=1
    return dict(flat_log_and_kappa_checks=flat,direct_positive_dsigma_C0_C1_comparisons=tests,
        sigma_prime_integrated_before_hyper_exponential_evaluation=True)


def run():
    start=time.perf_counter();source=independent_original_source();flat=independent_flat_and_stieltjes()
    print('Actual regular/singular Laurent source, cut-F units and d-sigma identities pass',flush=True)
    f=OriginalN1InnerM15Numeric(cells=16,precision=80)
    coarse=OriginalN1InnerM15Numeric(cells=8,cutoff=f.cut,core=f.core)
    audit=mp.mp.clone();audit.dps=180
    def bounds(row):return tuple(audit.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples'])
    samples=[];conditions=zeros=seams=history=atoms=overlaps=0
    for Z in ('-.3','0','.3'):
        for q in ('0','.5','.625','.75','1','1.5','2'):
            p=f.evaluate(Z,q);samples.append(p)
            assert all(p['contract'].values());conditions+=len(p['contract'])
            assert p['m15_C0_C1_terms'][:2]==p['core_inlet']['m15_core_C0_C1_terms'];seams+=1
            for term in p['m15_C0_C1_terms'][2:]:
                for row in term['coefficient_C0_C1']:
                    lo,hi=bounds(row);a=audit.mpf(row['approximation']);e=audit.mpf(row['absolute_error_upper'])
                    assert max(abs(lo-a),abs(hi-a))<=e;conditions+=1
                    if q=='0':assert lo==hi==a==0;zeros+=1
            if q!='0':
                singular=p['m15_C0_C1_terms'][4]
                lo,hi=bounds(singular['log_scale_terms'][-1]);assert audit.isfinite(lo) and hi<-audit.mpf('1e100');atoms+=1
                assert singular['actual_function_binding']['actual_phase']=='s=q*(sc/4)';atoms+=1
            if q in ('1','1.5','2'):
                at=f.evaluate(Z,'.75')['m15_C0_C1_terms'][2]
                assert p['m15_C0_C1_terms'][2]==at;history+=1
        lower=f.evaluate(Z,'2');old=f.evaluate(Z,'1')
        a=lower['m15_C0_C1_terms'][3]['coefficient_C0_C1'][0]['approximation']
        b=old['m15_C0_C1_terms'][3]['coefficient_C0_C1'][0]['approximation']
        assert audit.mpf(a)!=audit.mpf(b);history+=1
        small=coarse.evaluate(Z,'2');samples.append(small)
        for a,b in zip(small['m15_C0_C1_terms'][2:],lower['m15_C0_C1_terms'][2:]):
            for x,y in zip(a['coefficient_C0_C1'],b['coefficient_C0_C1']):
                p0,p1=bounds(x),bounds(y);assert max(p0[0],p1[0])<=min(p0[1],p1[1]);overlaps+=1
        print('Signed core/first/lower m15 cumulative history reaches rminus: Z='+Z,flush=True)
    rejected=0
    for Z,q in [('1','1'),('nan','0'),('.3','-.1'),('.3','2.1')]:
        try:f.evaluate(Z,q)
        except ValueError:rejected+=1
        else:raise AssertionError('Outside selected inner/lower chart admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,independent_actual_source=source,independent_flat_and_stieltjes=flat,
            actual_signed_core_first_lower_packets=len(samples),source_and_serialized_error_conditions=conditions,
            unchanged_signed_core_inlet_conditions=seams,exact_empty_increment_zeros=zeros,
            preserved_flat_log_atom_conditions=atoms,velocity_stop_and_nonzero_pressure_history_conditions=history,
            directed_cell_refinement_C0_C1_comparisons=overlaps,rejected_requests=rejected),
        execution_seconds=round(time.perf_counter()-start,3),
        limits=dict(signed_core_first_lower_m15_C0_C1=True,
            signed_inner_lower_boundaries_numerically_available=True,active_compact_signed_pressure=False,
            independent_preRp_signed_numeric_seam=False,global_five_moment_numeric_join=False,
            common_global_N=False,exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,
            full_corrected_NS=False))


def inspect(path):
    d=json.loads(path.read_text(encoding='utf-8'));assert d['checks']['passed']
    for name,digest in d['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in d['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return d


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME;d=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(d,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=d['checks'],execution_seconds=d['execution_seconds']),indent=2),flush=True)
